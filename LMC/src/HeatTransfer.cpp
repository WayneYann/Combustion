//
// "Divu_Type" means S, where divergence U = S
// "Dsdt_Type" means pd S/pd t, where S is as above
// "Ydot_Type" means -omega_l/rho, i.e., the mass rate of decrease of species l due
//             to kinetics divided by rho
//
#include <winstd.H>
#include <unistd.h>

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstdio>
#include <cfloat>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <ctime>

using std::cout;
using std::endl;
using std::cerr;

#include <Geometry.H>
#include <Extrapolater.H>
#include <BoxDomain.H>
#include <ParmParse.H>
#include <ErrorList.H>
#include <HeatTransfer.H>
#include <HEATTRANSFER_F.H>
#include <DIFFUSION_F.H>
#include <MultiGrid.H>
#include <ArrayLim.H>
#include <SPACE.H>
#include <Interpolater.H>
#include <ccse-mpi.H>
#include <Utility.H>
#include <BLProfiler.H>

#if defined(BL_USE_NEWMECH) || defined(BL_USE_VELOCITY)
#include <DataServices.H>
#include <AmrData.H>
#endif

#include <PROB_F.H>
#include <NAVIERSTOKES_F.H>
#include <DERIVE_F.H>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <buildInfo.H>

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITSCOMP(fab,fabdat,fablo,fabhi,comp)  \
const int* fablo = (fab).loVect();                    \
const int* fabhi = (fab).hiVect();                    \
const Real* fabdat = (fab).dataPtr(comp);

#ifdef BL_USE_FLOAT
#  define Real_MIN FLT_MIN
#  define Real_MAX FLT_MAX
#else
#  define Real_MIN DBL_MIN
#  define Real_MAX DBL_MAX
#endif

#define GEOM_GROW   1
#define PRESS_GROW  1
#define DIVU_GROW   1
#define DSDT_GROW   1
#define bogus_value 1.e20
#define DQRAD_GROW  1
#define YDOT_GROW   1

const int LinOp_grow = 1;

static const std::string typical_values_filename("typical_values.fab");

static const Real typical_RhoH_value_default = -1.e10;

#define SHOWVALARR(val)                        \
{                                              \
    std::cout << #val << " = ";                \
    for (int i=0;i<val.size();++i)             \
    {                                          \
        std::cout << val[i] << " " ;           \
    }                                          \
    std::cout << std::endl;                    \
}                                             
#define SHOWVALARRA(val) { SHOWVALARR(val); BoxLib::Abort();}
#define SHOWVAL(val) { std::cout << #val << " = " << val << std::endl;}
#define SHOWVALA(val) { SHOWVAL(val); BoxLib::Abort();}

namespace
{
    bool initialized  = false;

}
//
// Set all default values in Initialize()!!!
//
namespace
{
    std::set<std::string> ShowMF_Sets;
    std::string           ShowMF_Dir;
    bool                  ShowMF_Verbose;
    bool                  ShowMF_Check_Nans;
    bool                  do_not_use_funccount;
    bool                  do_active_control;
    bool                  do_active_control_temp;
    Real                  temp_control;
    Real                  crse_dt;
    bool                  benchmarking;
}

int  HeatTransfer::num_divu_iters;
int  HeatTransfer::init_once_done;
int  HeatTransfer::do_OT_radiation;
int  HeatTransfer::do_heat_sink;
int  HeatTransfer::RhoH;
int  HeatTransfer::do_diffuse_sync;
int  HeatTransfer::do_reflux_visc;
int  HeatTransfer::dpdt_option;
int  HeatTransfer::Ydot_Type;
int  HeatTransfer::FuncCount_Type;
int  HeatTransfer::divu_ceiling;
Real HeatTransfer::divu_dt_factor;
Real HeatTransfer::min_rho_divu_ceiling;
int  HeatTransfer::have_trac;
int  HeatTransfer::have_rhort;
int  HeatTransfer::Trac;
int  HeatTransfer::RhoRT;
int  HeatTransfer::first_spec;
int  HeatTransfer::last_spec;
int  HeatTransfer::nspecies;
int  HeatTransfer::floor_species;
int  HeatTransfer::do_set_rho_to_species_sum;
Real HeatTransfer::rgas;
Real HeatTransfer::prandtl;
Real HeatTransfer::schmidt;
Real HeatTransfer::constant_mu_val;
Real HeatTransfer::constant_rhoD_val;
Real HeatTransfer::constant_lambda_val;
int  HeatTransfer::unity_Le;
Real HeatTransfer::htt_tempmin;
Real HeatTransfer::htt_tempmax;
Real HeatTransfer::htt_hmixTYP;
int  HeatTransfer::siegel_test;
int  HeatTransfer::zeroBndryVisc;
int  HeatTransfer::do_add_nonunityLe_corr_to_rhoh_adv_flux;
int  HeatTransfer::do_check_divudt;
int  HeatTransfer::hack_nochem;
int  HeatTransfer::hack_nospecdiff;
int  HeatTransfer::hack_noavgdivu;
int  HeatTransfer::use_tranlib;
Real HeatTransfer::trac_diff_coef;
Real HeatTransfer::P1atm_MKS;
bool HeatTransfer::plot_reactions;
bool HeatTransfer::plot_consumption;
bool HeatTransfer::plot_heat_release;
Real HeatTransfer::new_T_threshold;
int  HeatTransfer::reset_typical_vals_int=-1;
std::map<std::string,Real> HeatTransfer::typical_values_FileVals;

std::string                                HeatTransfer::turbFile;
ChemDriver*                                HeatTransfer::chemSolve;
std::map<std::string, Array<std::string> > HeatTransfer::auxDiag_names;

Array<Real> HeatTransfer::typical_values;

void
HeatTransfer::Initialize ()
{
    if (initialized) return;

    NavierStokesBase::Initialize();

    //
    // Set all default values here!!!
    //
    ShowMF_Verbose         = true;
    ShowMF_Check_Nans      = true;
    do_not_use_funccount   = false;
    do_active_control      = false;
    do_active_control_temp = false;
    temp_control           = -1;
    benchmarking           = false;
    crse_dt                = -1;
    
    HeatTransfer::num_divu_iters            = 1;
    HeatTransfer::init_once_done            = 0;
    HeatTransfer::do_OT_radiation           = 0;
    HeatTransfer::do_heat_sink              = 0;
    HeatTransfer::RhoH                      = -1;
    HeatTransfer::do_diffuse_sync           = 1;
    HeatTransfer::do_reflux_visc            = 1;
    HeatTransfer::dpdt_option               = 2;
    HeatTransfer::Ydot_Type                 = -1;
    HeatTransfer::FuncCount_Type            = -1;
    HeatTransfer::divu_ceiling              = 0;
    HeatTransfer::divu_dt_factor            = .5;
    HeatTransfer::min_rho_divu_ceiling      = -1.e20;
    HeatTransfer::have_trac                 = 0;
    HeatTransfer::have_rhort                = 0;
    HeatTransfer::Trac                      = -1;
    HeatTransfer::RhoRT                     = -1;
    HeatTransfer::first_spec                = -1;
    HeatTransfer::last_spec                 = -2;
    HeatTransfer::nspecies                  = 0;
    HeatTransfer::floor_species             = 0;
    HeatTransfer::do_set_rho_to_species_sum = 1;
    HeatTransfer::rgas                      = -1.0;
    HeatTransfer::prandtl                   = .7;
    HeatTransfer::schmidt                   = .7;
    HeatTransfer::constant_mu_val           = -1;
    HeatTransfer::constant_rhoD_val         = -1;
    HeatTransfer::constant_lambda_val       = -1;
    HeatTransfer::unity_Le                  = 1;
    HeatTransfer::htt_tempmin               = 298.0;
    HeatTransfer::htt_tempmax               = 40000.;
    HeatTransfer::htt_hmixTYP               = -1.;
    HeatTransfer::siegel_test               = 0;
    HeatTransfer::zeroBndryVisc             = 0;
    HeatTransfer::chemSolve                 = 0;
    HeatTransfer::do_check_divudt           = 1;
    HeatTransfer::hack_nochem               = 0;
    HeatTransfer::hack_nospecdiff           = 0;
    HeatTransfer::hack_noavgdivu            = 0;
    HeatTransfer::use_tranlib               = 0;
    HeatTransfer::trac_diff_coef            = 0.0;
    HeatTransfer::P1atm_MKS                 = -1.0;
    HeatTransfer::turbFile                  = "";
    HeatTransfer::plot_reactions            = false;
    HeatTransfer::plot_consumption          = true;
    HeatTransfer::plot_heat_release         = true;
    HeatTransfer::new_T_threshold           = -1;  // On new AMR level, max change in lower bound for T, not used if <=0

    HeatTransfer::do_add_nonunityLe_corr_to_rhoh_adv_flux = 1;

    HeatTransfer::reset_typical_vals_int    = -1;
    HeatTransfer::typical_values_FileVals.clear();

    ParmParse pp("ns");

    pp.query("benchmarking",benchmarking);

    pp.query("do_diffuse_sync",do_diffuse_sync);
    BL_ASSERT(do_diffuse_sync == 0 || do_diffuse_sync == 1);
    pp.query("do_reflux_visc",do_reflux_visc);
    BL_ASSERT(do_reflux_visc == 0 || do_reflux_visc == 1);
    pp.query("dpdt_option",dpdt_option);
    BL_ASSERT(dpdt_option >= 0 && dpdt_option <= 2);
    pp.query("do_active_control",do_active_control);
    pp.query("do_active_control_temp",do_active_control_temp);
    pp.query("temp_control",temp_control);

    if (do_active_control_temp && temp_control <= 0)
        BoxLib::Error("temp_control MUST be set with do_active_control_temp");

    verbose = pp.contains("v");

    pp.query("divu_ceiling",divu_ceiling);
    BL_ASSERT(divu_ceiling >= 0 && divu_ceiling <= 3);
    pp.query("divu_dt_factor",divu_dt_factor);
    BL_ASSERT(divu_dt_factor>0 && divu_dt_factor <= 1.0);
    pp.query("min_rho_divu_ceiling",min_rho_divu_ceiling);
    if (divu_ceiling) BL_ASSERT(min_rho_divu_ceiling >= 0.0);

    pp.query("htt_tempmin",htt_tempmin);
    pp.query("htt_tempmax",htt_tempmax);

    pp.query("floor_species",floor_species);
    BL_ASSERT(floor_species == 0 || floor_species == 1);

    pp.query("do_set_rho_to_species_sum",do_set_rho_to_species_sum);

    pp.query("num_divu_iters",num_divu_iters);

    pp.query("do_not_use_funccount",do_not_use_funccount);

    pp.query("schmidt",schmidt);
    pp.query("prandtl",prandtl);
    pp.query("unity_Le",unity_Le);
    unity_Le = unity_Le ? 1 : 0;
    if (unity_Le)
    {
	schmidt = prandtl;
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::read_params: Le=1, setting Sc = Pr" << '\n';
    }

    pp.query("constant_mu_val",constant_mu_val);
    pp.query("constant_rhoD_val",constant_rhoD_val);
    pp.query("constant_lambda_val",constant_lambda_val);
    if (constant_mu_val != -1)
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::read_params: using constant_mu_val = " 
                      << constant_mu_val << '\n';
    }
    if (constant_rhoD_val != -1)
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::read_params: using constant_rhoD_val = " 
                      << constant_rhoD_val << '\n';
    }
    if (constant_lambda_val != -1)
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::read_params: using constant_lambda_val = " 
                      << constant_lambda_val << '\n';
    }

    pp.query("do_add_nonunityLe_corr_to_rhoh_adv_flux", do_add_nonunityLe_corr_to_rhoh_adv_flux);
    pp.query("hack_nochem",hack_nochem);
    pp.query("hack_nospecdiff",hack_nospecdiff);
    pp.query("hack_noavgdivu",hack_noavgdivu);
    pp.query("do_check_divudt",do_check_divudt);
    pp.query("do_OT_radiation",do_OT_radiation);
    do_OT_radiation = (do_OT_radiation ? 1 : 0);
    pp.query("do_heat_sink",do_heat_sink);
    do_heat_sink = (do_heat_sink ? 1 : 0);

    pp.query("use_tranlib",use_tranlib);
    if (use_tranlib == 1) {
      chemSolve->SetTransport(ChemDriver::CD_TRANLIB);
      if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "HeatTransfer::read_params: Using Tranlib transport " << '\n';
    }
    else {
      chemSolve->SetTransport(ChemDriver::CD_EG);
      if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "HeatTransfer::read_params: Using EGLib transport " << '\n';
    }
    chemSolve = new ChemDriver();

    pp.query("turbFile",turbFile);

    pp.query("siegel_test",siegel_test);
    pp.query("zeroBndryVisc",zeroBndryVisc);
    //
    // Set variability/visc for velocities.
    //
    if (variable_vel_visc != 1)
      BoxLib::Error("HeatTransfer::read_params() -- must use variable viscosity");
    //
    // Set variability/visc for tracer
    //
    if (variable_scal_diff != 1)
      BoxLib::Error("HeatTransfer::read_params() -- must use variable scalar diffusivity");
    //
    // Read in scalar value and use it as tracer.
    //
    pp.query("scal_diff_coefs",trac_diff_coef);

    for (int i = 0; i < visc_coef.size(); i++)
        visc_coef[i] = bogus_value;

    // Useful for debugging
    if (int nsv=pp.countval("ShowMF_Sets"))
    {
        Array<std::string> ShowMF_set_names(nsv);
        pp.getarr("ShowMF_Sets",ShowMF_set_names);
        for (int i=0; i<nsv; ++i) {
            ShowMF_Sets.insert(ShowMF_set_names[i]);
        }
        ShowMF_Dir="."; pp.query("ShowMF_Dir",ShowMF_Dir);
        pp.query("ShowMF_Verbose",ShowMF_Verbose);
        pp.query("ShowMF_Check_Nans",ShowMF_Check_Nans);
    }

#ifdef PARTICLES
    read_particle_params ();
#endif
        
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "\nDumping ParmParse table:\n \n";
        ParmParse::dumpTable(std::cout);
        std::cout << "\n... done dumping ParmParse table.\n" << '\n';
    }

    BoxLib::ExecOnFinalize(HeatTransfer::Finalize);

    initialized = true;
}

void
HeatTransfer::Finalize ()
{
    initialized = false;
}

static
void
FabMinMax (FArrayBox& fab,
           const Box& box,
           Real       fmin,
           Real       fmax,
           int        sComp,
           int        nComp)
{
    BL_ASSERT(fab.box().contains(box));
    BL_ASSERT(sComp + nComp <= fab.nComp());

    const int* lo     = box.loVect();
    const int* hi     = box.hiVect();
    Real*      fabdat = fab.dataPtr(sComp);
    
    FORT_FABMINMAX(lo, hi,
                   fabdat, ARLIM(fab.loVect()), ARLIM(fab.hiVect()),
                   &fmin, &fmax, &nComp);
}

static
std::ostream&
levWrite(std::ostream& os, int level, std::string& message)
{
    if (ParallelDescriptor::IOProcessor())
    {
        for (int lev=0; lev<=level; ++lev) {
            //os << '\t';
            os << "  ";
        }
        os << message;
    }
    return os;
}

void
showMF(const std::string&   mySet,
       const MultiFab&      mf,
       const std::string&   name,
       int                  lev = -1,
       int                  iter = -1) // Default value = no append 2nd integer
{
    if (ShowMF_Sets.count(mySet)>0)
    {
        std::string DebugDir(ShowMF_Dir);
        if (ParallelDescriptor::IOProcessor())
            if (!BoxLib::UtilCreateDirectory(DebugDir, 0755))
                BoxLib::CreateDirectoryFailed(DebugDir);
        ParallelDescriptor::Barrier();

        std::string junkname = name;
        if (lev>=0) {
            junkname = BoxLib::Concatenate(junkname+"_",lev,1);
        }
        if (iter>=0) {
            junkname = BoxLib::Concatenate(junkname+"_",iter,1);
        }
        junkname = DebugDir + "/" + junkname;

        if (ShowMF_Verbose>0 && ParallelDescriptor::IOProcessor()) {
            cout << "   ******************************  Debug: writing " << junkname << endl;
        }

        if (ShowMF_Check_Nans)
        {
            for (MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                BL_ASSERT(!mf[mfi].contains_nan(mfi.validbox(),0,mf.nComp()));
            }
        }
        VisMF::Write(mf,junkname);
    }
}

HeatTransfer::FPLoc 
HeatTransfer::fpi_phys_loc (int p_bc)
{
    //
    // Location of data that FillPatchIterator returns at physical boundaries
    //
    if (p_bc == EXT_DIR || p_bc == HOEXTRAP || p_bc == FOEXTRAP)
    {
        return HT_Edge;
    }
    return HT_Center;
}
    
void
HeatTransfer::center_to_edge_fancy (const FArrayBox& cfab,
                                    FArrayBox&       efab,
                                    const Box&       ccBox,
                                    int              sComp,
                                    int              dComp,
                                    int              nComp,
                                    const Box&       domain,
                                    const FPLoc&     bc_lo,
                                    const FPLoc&     bc_hi)
{
    const Box&      ebox = efab.box();
    const IndexType ixt  = ebox.ixType();

    BL_ASSERT(!(ixt.cellCentered()) && !(ixt.nodeCentered()));

    int dir = -1;
    for (int d = 0; d < BL_SPACEDIM; d++)
        if (ixt.test(d))
            dir = d;

    BL_ASSERT(BoxLib::grow(ccBox,-BoxLib::BASISV(dir)).contains(BoxLib::enclosedCells(ebox)));
    BL_ASSERT(sComp+nComp <= cfab.nComp() && dComp+nComp <= efab.nComp());
    //
    // Exclude unnecessary cc->ec calcs
    //
    Box ccVBox = ccBox;
    if (bc_lo!=HT_Center)
        ccVBox.setSmall(dir,std::max(domain.smallEnd(dir),ccVBox.smallEnd(dir)));
    if (bc_hi!=HT_Center)
        ccVBox.setBig(dir,std::min(domain.bigEnd(dir),ccVBox.bigEnd(dir)));
    //
    // Shift cell-centered data to edges
    //
    const int isharm = def_harm_avg_cen2edge?1:0;
    FORT_CEN2EDG(ccVBox.loVect(),ccVBox.hiVect(),
                 ARLIM(cfab.loVect()),ARLIM(cfab.hiVect()),cfab.dataPtr(sComp),
                 ARLIM(efab.loVect()),ARLIM(efab.hiVect()),efab.dataPtr(dComp),
                 &nComp, &dir, &isharm);
    //
    // Fix boundary...i.e. fill-patched data in cfab REALLY lives on edges
    //
    if ( !(domain.contains(ccBox)) )
    {
        if (bc_lo==HT_Edge)
        {
            BoxList gCells = BoxLib::boxDiff(ccBox,domain);
            if (gCells.ok())
            {
                const int inc = +1;
                FArrayBox ovlpFab;
                for (BoxList::const_iterator bli = gCells.begin(), end = gCells.end();
                     bli != end;
                     ++bli)
                {
                    if (bc_lo == HT_Edge)
                    {
                        ovlpFab.resize(*bli,nComp);
                        ovlpFab.copy(cfab,sComp,0,nComp);
                        ovlpFab.shiftHalf(dir,inc);
                        efab.copy(ovlpFab,0,dComp,nComp);
                    }
                }
            }
        }
        if (bc_hi==HT_Edge)
        {
            BoxList gCells = BoxLib::boxDiff(ccBox,domain);
            if (gCells.ok())
            {
                const int inc = -1;
                FArrayBox ovlpFab;
                for (BoxList::const_iterator bli = gCells.begin(), end = gCells.end();
                     bli != end;
                     ++bli)
                {
                    if (bc_hi == HT_Edge)
                    {
                        ovlpFab.resize(*bli,nComp);
                        ovlpFab.copy(cfab,sComp,0,nComp);
                        ovlpFab.shiftHalf(dir,inc);
                        efab.copy(ovlpFab,0,dComp,nComp);
                    }
                }
            }
        }
    }
}    

void
HeatTransfer::variableCleanUp ()
{
    NavierStokesBase::variableCleanUp();

    delete chemSolve;
    chemSolve = 0;

    ShowMF_Sets.clear();
    auxDiag_names.clear();
    typical_values.clear();
}

HeatTransfer::HeatTransfer ()
{
    if (!init_once_done)
        init_once();

    if (!do_temp)
        BoxLib::Abort("do_temp MUST be true");

    if (!have_divu)
        BoxLib::Abort("have_divu MUST be true");

    if (!have_dsdt)
        BoxLib::Abort("have_dsdt MUST be true");

    EdgeState              = 0;
    SpecDiffusionFluxn     = 0;
    SpecDiffusionFluxnp1   = 0;
    FillPatchedOldState_ok = true;
    FillPatchedNewState_ok = true;
}

HeatTransfer::HeatTransfer (Amr&            papa,
                            int             lev,
                            const Geometry& level_geom,
                            const BoxArray& bl,
                            Real            time)
    :
    NavierStokesBase(papa,lev,level_geom,bl,time),
    //
    // Make room for all components except velocities in aux_boundary_data_old.
    //
    aux_boundary_data_old(bl,Godunov::hypgrow(),desc_lst[State_Type].nComp()-BL_SPACEDIM,level_geom),
    //
    // Only save RhoH in aux_boundary_data_new in component 0.
    //
    aux_boundary_data_new(bl,LinOp_grow,1,level_geom),
    FillPatchedOldState_ok(true),
    FillPatchedNewState_ok(true)
{
    if (!init_once_done)
        init_once();

    if (!do_temp)
        BoxLib::Abort("do_temp MUST be true");

    if (!have_divu)
        BoxLib::Abort("have_divu MUST be true");

    if (!have_dsdt)
        BoxLib::Abort("have_dsdt MUST be true");

    define_data();
}

HeatTransfer::~HeatTransfer ()
{
    ;
}

void
HeatTransfer::define_data ()
{
    const int nGrow       = 0;
    const int nEdgeStates = desc_lst[State_Type].nComp();

    EdgeState = (raii_fbs.push_back(new FluxBoxes(this, nEdgeStates, nGrow)))->get();
    
    if (nspecies>0 && !unity_Le)
    {
	SpecDiffusionFluxn   = (raii_fbs.push_back(
				    new FluxBoxes(this, nspecies, nGrow)))->get();
	SpecDiffusionFluxnp1 = (raii_fbs.push_back(
				    new FluxBoxes(this, nspecies, nGrow)))->get();
	spec_diffusion_flux_computed.resize(nspecies,HT_None);
    }

    for (std::map<std::string,Array<std::string> >::iterator it = auxDiag_names.begin(), end = auxDiag_names.end();
         it != end; ++it)
    {
        auxDiag[it->first] = raii_mfs.push_back(new MultiFab(grids,it->second.size(),0));
        auxDiag[it->first]->setVal(0);
    }

}

void
HeatTransfer::init_once ()
{
    //
    // Computes the static variables unique to HeatTransfer.
    // Check that (some) things are set up correctly.
    //
    int dummy_State_Type;

    int have_temp = isStateVariable("temp", dummy_State_Type, Temp);

    have_temp = have_temp && State_Type == dummy_State_Type;
    have_temp = have_temp && isStateVariable("rhoh", dummy_State_Type, RhoH);
    have_temp = have_temp && State_Type == dummy_State_Type;

    have_trac = isStateVariable("tracer", dummy_State_Type, Trac);
    have_trac = have_trac && State_Type == dummy_State_Type;

    have_rhort = isStateVariable("RhoRT", dummy_State_Type, RhoRT);
    have_rhort = have_rhort && State_Type == dummy_State_Type;

    if (!have_temp)
        BoxLib::Abort("HeatTransfer::init_once(): RhoH & Temp must both be the state");
    
    if (!have_rhort && verbose && ParallelDescriptor::IOProcessor())
        BoxLib::Warning("HeatTransfer::init_once(): RhoRT being stored in the Tracer slot");
    
    if (Temp < RhoH)
        BoxLib::Abort("HeatTransfer::init_once(): must have RhoH < Temp");
    //
    // Temperature must be non-conservative, rho*h must be conservative.
    //
    if (advectionType[Temp] == Conservative)
        BoxLib::Abort("HeatTransfer::init_once(): Temp must be non-conservative");

    if (advectionType[RhoH] != Conservative)
        BoxLib::Abort("HeatTransfer::init_once(): RhoH must be conservative");
    //
    // Species checks.
    //
    BL_ASSERT(Temp > RhoH && RhoH > Density);
    //
    // Here we want to count relative to Density instead of relative
    // to RhoH, so we can put RhoH after the species instead of before.  
    // This logic should work in both cases.
    //
    first_spec =  Density + 1;
    last_spec  = first_spec + getChemSolve().numSpecies() - 1;
    
    for (int i = first_spec; i <= last_spec; i++)
        if (advectionType[i] != Conservative)
            BoxLib::Error("HeatTransfer::init_once: species must be conservative");
    
    int diffuse_spec = is_diffusive[first_spec];
    for (int i = first_spec+1; i <= last_spec; i++)
        if (is_diffusive[i] != diffuse_spec)
            BoxLib::Error("HeatTransfer::init_once: Le != 1; diffuse");
    //
    // Load integer pointers into Fortran common, reqd for proper ICs.
    //
    const int density = (int)Density;

    FORT_SET_SCAL_NUMB(&density, &Temp, &Trac, &RhoH, &first_spec, &last_spec);
    //
    // Load constants into Fortran common to compute viscosities, etc.
    //
    const int var_visc = (constant_mu_val     == -1 ? 1 : 0);
    const int var_cond = (constant_lambda_val == -1 ? 1 : 0);
    const int var_diff = (constant_rhoD_val   == -1 ? 1 : 0);
    
    FORT_SET_HT_VISC_COMMON(&var_visc, &constant_mu_val,
                            &var_cond, &constant_lambda_val,
                            &var_diff, &constant_rhoD_val,
                            &prandtl,  &schmidt, &unity_Le);

    // initialize default typical values
    FORT_INIT_TYPVALS_COMMON();
    //
    // make space for typical values
    //
    typical_values.resize(NUM_STATE,-1); // -ve means don't use for anything
    typical_values[RhoH] = typical_RhoH_value_default;

    ParmParse pp("ht");

    const Array<std::string>& speciesNames = getChemSolve().speciesNames();

    for (int i=0; i<nspecies; ++i) {
      const std::string ppStr = std::string("typValY_") + speciesNames[i];
      if (pp.countval(ppStr.c_str())>0) {
        pp.get(ppStr.c_str(),typical_values_FileVals[speciesNames[i]]);
      }
    }
    std::string otherKeys[4] = {"Temp", "RhoH", "Vel", "Trac"};
    for (int i=0; i<4; ++i) {
      const std::string ppStr(std::string("typVal_")+otherKeys[i]);
      if (pp.countval(ppStr.c_str())>0) {
        pp.get(ppStr.c_str(),typical_values_FileVals[otherKeys[i]]);
      }
    }
    //
    // Get universal gas constant from Fortran.
    //
    rgas = getChemSolve().getRuniversal();
    P1atm_MKS = getChemSolve().getP1atm_MKS();

    if (rgas <= 0.0)
    {
        std::cerr << "HeatTransfer::init_once(): bad rgas: " << rgas << '\n';
        BoxLib::Abort();
    }
    if (P1atm_MKS <= 0.0)
    {
        std::cerr << "HeatTransfer::init_once(): bad P1atm_MKS: " << P1atm_MKS << '\n';
        BoxLib::Abort();
    }
    //
    // Chemistry.
    //
    int ydot_good = Ydot_Type >= 0 && Ydot_Type <desc_lst.size()
        && Ydot_Type != Divu_Type
        && Ydot_Type != Dsdt_Type
        && Ydot_Type != State_Type;
    
    if (!ydot_good)
        BoxLib::Error("HeatTransfer::init_once(): need Ydot_Type if do_chemistry");
    
    const StateDescriptor& ydot_cell = desc_lst[Ydot_Type];
    int nydot = ydot_cell.nComp();
    if (nydot < nspecies)
        BoxLib::Error("HeatTransfer::init_once(): Ydot_Type needs nspecies components");
    //
    // Enforce Le = 1, unless !unity_Le
    //
    if (unity_Le && (schmidt != prandtl) )
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "**************** WARNING ***************\n"
                      << "HeatTransfer::init_once() : currently must have"
                      << "equal Schmidt and Prandtl numbers unless !unity_Le.\n"
                      << "Setting Schmidt = Prandtl\n"
                      << "**************** WARNING ***************\n";
    
        schmidt = prandtl;
    }
    //
    // We are done.
    //
    num_state_type = desc_lst.size();

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "HeatTransfer::init_once(): num_state_type = " << num_state_type << '\n';

    pp.query("plot_reactions",plot_reactions);
    if (plot_reactions)
    {
        auxDiag_names["REACTIONS"].resize(getChemSolve().numReactions());
        for (int i = 0; i < auxDiag_names["REACTIONS"].size(); ++i)
            auxDiag_names["REACTIONS"][i] = BoxLib::Concatenate("R",i+1);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "***** Make sure to increase amr.regrid_int !!!!!" << '\n';
    }

    pp.query("plot_consumption",plot_consumption);
    pp.query("plot_auxDiags",plot_consumption); // This is for backward comptibility - FIXME
    if (plot_consumption)
    {
        auxDiag_names["CONSUMPTION"].resize(consumptionName.size());
        for (int j=0; j<consumptionName.size(); ++j)
        {
            auxDiag_names["CONSUMPTION"][j] = consumptionName[j] + "_ConsumptionRate";
        }
    }

    pp.query("plot_heat_release",plot_heat_release);
    if (plot_heat_release)
    {
        auxDiag_names["HEATRELEASE"].resize(1);
        auxDiag_names["HEATRELEASE"][0] = "HeatRelease";
    }
    
    pp.query("new_T_threshold",new_T_threshold);

#ifdef BL_COMM_PROFILING
    auxDiag_names["COMMPROF"].resize(3);
    auxDiag_names["COMMPROF"][0] = "mpiRank";
    auxDiag_names["COMMPROF"][1] = "proximityRank";
    auxDiag_names["COMMPROF"][2] = "proximityOrder";
#endif

    init_once_done = 1;
}

void
HeatTransfer::restart (Amr&          papa,
                       std::istream& is,
                       bool          bReadSpecial)
{

    NavierStokesBase::restart(papa,is,bReadSpecial);
    //
    // Make room for all components except velocities in aux_boundary_data_old.
    //
    aux_boundary_data_old.initialize(grids,Godunov::hypgrow(),desc_lst[State_Type].nComp()-BL_SPACEDIM,Geom());
    //
    // Only save RhoH in aux_boundary_data_new in component 0.
    //
    aux_boundary_data_new.initialize(grids,LinOp_grow,1,Geom());

    FillPatchedOldState_ok = true;
    FillPatchedNewState_ok = true;

    define_data();

    // Deal with typical values
    set_typical_values(true);
}

void
HeatTransfer::set_typical_values (bool restart)
{
    if (level==0) {

        int nComp = typical_values.size();

        if (restart) {

            BL_ASSERT(nComp==NUM_STATE);
            for (int i=0; i<nComp; ++i) {
                typical_values[i] = -1;
            }
            
            if (ParallelDescriptor::IOProcessor())
            {
                const std::string tvfile = parent->theRestartFile() + "/" + typical_values_filename;
                std::ifstream tvis;
                tvis.open(tvfile.c_str(),std::ios::in|std::ios::binary);
                
                if (tvis.good()) {
                    FArrayBox tvfab;
                    tvfab.readFrom(tvis);
                    if (tvfab.nComp() != typical_values.size()) {
                        BoxLib::Abort("Typical values file has wrong number of components");
                    }
                    for (int i=0; i<typical_values.size(); ++i) {
                        typical_values[i] = tvfab.dataPtr()[i];
                    }
                }
            }

	    ParallelDescriptor::ReduceRealMax(typical_values.dataPtr(),nComp); //FIXME: better way?
        }
	else { // not restart

	  // Check fortan common values, override values set above if fortran values > 0
	  Array<Real> tvTmp(nComp,0);
	  FORT_GETTYPICALVALS(tvTmp.dataPtr(), &nComp);
	  ParallelDescriptor::ReduceRealMax(tvTmp.dataPtr(),nComp);
        
	  for (int i=0; i<nComp; ++i) {
            if (tvTmp[i]>0) {
	      typical_values[i] = tvTmp[i];
            }
	  }
	}

        // If typVals specified in inputs, these take precedence componentwise
        for (std::map<std::string,Real>::const_iterator it=typical_values_FileVals.begin(), 
               End=typical_values_FileVals.end(); it!=End; ++it) {
          int idx = getChemSolve().index(it->first);
          if (idx>=0) {
            typical_values[first_spec+idx] = it->second;
          } else {
            if (it->first == "Temp") {
              typical_values[Temp] = it->second;
            }
            else if (it->first == "RhoH") {
              typical_values[RhoH] = it->second;
            }
            else if (it->first == "Trac") {
              typical_values[Trac] = it->second;
            }
            else if (it->first == "Vel") {
              for (int d=0; d<BL_SPACEDIM; ++d) {
                typical_values[d] = it->second;
              }
            }
          }
        }

        FORT_SETTYPICALVALS(typical_values.dataPtr(), &nComp);

        if (ParallelDescriptor::IOProcessor())
        {
            cout << "Typical vals: " << '\n';
            cout << "\tVelocity: ";
            for (int i=0; i<BL_SPACEDIM; ++i) {
                cout << typical_values[i] << " ";
            }
            cout << '\n';
            cout << "\tDensity: " << typical_values[Density] << '\n';
            cout << "\tTemp: "    << typical_values[Temp]    << '\n';
            cout << "\tRhoH: "    << typical_values[RhoH]    << '\n';
            const Array<std::string>& names = getChemSolve().speciesNames();
            for (int i=0; i<nspecies; ++i) {
                cout << "\tY_" << names[i] << ": " << typical_values[first_spec+i] << '\n';
            }
        }
    }
}

void
HeatTransfer::reset_typical_values(const MultiFab& S)
{
  // NOTE: Assumes that this level has valid data everywhere
  int nComp = typical_values.size();
  BL_ASSERT(nComp = S.nComp());
  for (int i=0; i<nComp; ++i) {
    Real thisMax = S.max(i);
    Real thisMin = S.min(i);
    Real newVal = std::abs(thisMax - thisMin);
    if (newVal > 0) {
      if ( (i>=first_spec && i<=last_spec) ) {
        typical_values[i] = newVal / typical_values[Density];
      }
      else {
        typical_values[i] = newVal;
      }
    }
  }

  // If typVals specified in inputs, these take precedence componentwise
  if (parent->levelSteps(0) == 0)
  {
    for (std::map<std::string,Real>::const_iterator it=typical_values_FileVals.begin(), 
	   End=typical_values_FileVals.end(); it!=End; ++it) {
      int idx = getChemSolve().index(it->first);
      if (idx>=0) {
	typical_values[first_spec+idx] = it->second;
      } else {
	if (it->first == "Temp") {
	  typical_values[Temp] = it->second;
	}
	else if (it->first == "RhoH") {
	  typical_values[RhoH] = it->second;
	}
	else if (it->first == "Trac") {
	  typical_values[Trac] = it->second;
	}
	else if (it->first == "Vel") {
	  for (int d=0; d<BL_SPACEDIM; ++d) {
	    typical_values[d] = it->second;
	  }
	}
      }
    }
  }

  FORT_SETTYPICALVALS(typical_values.dataPtr(), &nComp);

  if (ParallelDescriptor::IOProcessor()) {
    cout << "New typical vals: " << '\n';
    cout << "\tVelocity: ";
    for (int i=0; i<BL_SPACEDIM; ++i) {
      cout << typical_values[i] << " ";
    }
    cout << '\n';
    cout << "\tDensity: " << typical_values[Density] << '\n';
    cout << "\tTemp: "    << typical_values[Temp]    << '\n';
    cout << "\tRhoH: "    << typical_values[RhoH]    << '\n';
    const Array<std::string>& names = getChemSolve().speciesNames();
    for (int i=0; i<nspecies; ++i) {
      cout << "\tY_" << names[i] << ": " << typical_values[first_spec+i] << '\n';
    }
  }
}

Real
HeatTransfer::estTimeStep ()
{
    Real estdt = NavierStokesBase::estTimeStep();

    if (fixed_dt > 0.0 || !divu_ceiling)
        //
        // The n-s function did the right thing in this case.
        //
        return estdt;

    Real dt, ns_estdt = estdt, divu_dt = 1.0e20;

    const int   n_grow   = 1;
    const Real  cur_time = state[State_Type].curTime();
    const Real* dx       = geom.CellSize();
    MultiFab*   dsdt     = getDsdt(0,cur_time);
    MultiFab*   divu     = getDivCond(0,cur_time);

    for (FillPatchIterator U_fpi(*this,*divu,n_grow,cur_time,State_Type,Xvel,BL_SPACEDIM);
         U_fpi.isValid();
         ++U_fpi)
    {
        const int        i   = U_fpi.index();
        FArrayBox&       U   = U_fpi();
        const FArrayBox& Rho = rho_ctime[U_fpi];
        const int*       lo  = grids[i].loVect();
        const int*       hi  = grids[i].hiVect();

        DEF_CLIMITS((*divu)[U_fpi],sdat,slo,shi);
        DEF_CLIMITS(Rho,rhodat,rholo,rhohi);
        DEF_CLIMITS(U,vel,ulo,uhi);

        DEF_CLIMITS(volume[i],vol,v_lo,v_hi);

        DEF_CLIMITS(area[0][i],areax,ax_lo,ax_hi);
        DEF_CLIMITS(area[1][i],areay,ay_lo,ay_hi);
#if (BL_SPACEDIM==3)
        DEF_CLIMITS(area[2][i],areaz,az_lo,az_hi)
#endif
        FORT_EST_DIVU_DT(divu_ceiling,&divu_dt_factor,
                         dx,sdat,ARLIM(slo),ARLIM(shi),
                         (*dsdt)[U_fpi].dataPtr(),
                         rhodat,ARLIM(rholo),ARLIM(rhohi),
                         vel,ARLIM(ulo),ARLIM(uhi),
                         vol,ARLIM(v_lo),ARLIM(v_hi),
                         areax,ARLIM(ax_lo),ARLIM(ax_hi),
                         areay,ARLIM(ay_lo),ARLIM(ay_hi),
#if (BL_SPACEDIM==3) 
                         areaz,ARLIM(az_lo),ARLIM(az_hi),
#endif 
                         lo,hi,&dt,&min_rho_divu_ceiling);

        divu_dt = std::min(divu_dt,dt);
    }

    delete divu;
    delete dsdt;

    ParallelDescriptor::ReduceRealMin(divu_dt);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "HeatTransfer::estTimeStep(): estdt, divu_dt = " 
                  << estdt << ", " << divu_dt << '\n';
    }

    estdt = std::min(estdt, divu_dt);

    if (estdt < ns_estdt && verbose && ParallelDescriptor::IOProcessor())
        std::cout << "HeatTransfer::estTimeStep(): timestep reduced from " 
                  << ns_estdt << " to " << estdt << '\n';

    return estdt;
}

void
HeatTransfer::checkTimeStep (Real dt)
{
    if (fixed_dt > 0.0 || !divu_ceiling) 
        return;

    const int   n_grow    = 1;
    const Real  cur_time  = state[State_Type].curTime();
    const Real* dx        = geom.CellSize();
    MultiFab*   dsdt      = getDsdt(0,cur_time);
    MultiFab*   divu      = getDivCond(0,cur_time);

    for (FillPatchIterator U_fpi(*this,*divu,n_grow,cur_time,State_Type,Xvel,BL_SPACEDIM);
         U_fpi.isValid();
         ++U_fpi)
    {
        const int        i   = U_fpi.index();
        FArrayBox&       U   = U_fpi();
        const FArrayBox& Rho = rho_ctime[U_fpi];
        const int*       lo  = grids[i].loVect();
        const int*       hi  = grids[i].hiVect();

        DEF_LIMITS((*divu)[U_fpi],sdat,slo,shi);
        DEF_CLIMITS(Rho,rhodat,rholo,rhohi);
        DEF_CLIMITS(U,vel,ulo,uhi);

        DEF_CLIMITS(volume[i],vol,v_lo,v_hi);
        DEF_CLIMITS(area[0][i],areax,ax_lo,ax_hi);
        DEF_CLIMITS(area[1][i],areay,ay_lo,ay_hi);

#if (BL_SPACEDIM==3)
        DEF_CLIMITS(area[2][i],areaz,az_lo,az_hi);
#endif
        FORT_CHECK_DIVU_DT(divu_ceiling,&divu_dt_factor,
                           dx,sdat,ARLIM(slo),ARLIM(shi),
			   (*dsdt)[U_fpi].dataPtr(),
                           rhodat,ARLIM(rholo),ARLIM(rhohi),
                           vel,ARLIM(ulo),ARLIM(uhi),
                           vol,ARLIM(v_lo),ARLIM(v_hi),
                           areax,ARLIM(ax_lo),ARLIM(ax_hi),
                           areay,ARLIM(ay_lo),ARLIM(ay_hi),
#if (BL_SPACEDIM==3) 
                           areaz,ARLIM(az_lo),ARLIM(az_hi),
#endif 
                           lo,hi,&dt,&min_rho_divu_ceiling);
    }

    delete dsdt;
    delete divu;
}

void
HeatTransfer::setTimeLevel (Real time,
                            Real dt_old,
                            Real dt_new)
{
    NavierStokesBase::setTimeLevel(time, dt_old, dt_new);    

    state[Ydot_Type].setTimeLevel(time,dt_old,dt_new);

    state[FuncCount_Type].setTimeLevel(time,dt_old,dt_new);
}

//
// This (minus the NEWMECH stuff) is copied from NavierStokes.cpp
//

void
HeatTransfer::initData ()
{
    //
    // Initialize the state and the pressure.
    //
    int         ns       = NUM_STATE - BL_SPACEDIM;
    const Real* dx       = geom.CellSize();
    MultiFab&   S_new    = get_new_data(State_Type);
    MultiFab&   P_new    = get_new_data(Press_Type);
    const Real  cur_time = state[State_Type].curTime();

    ParmParse pp("ht");

#ifdef BL_USE_NEWMECH
    //
    // This code has a few drawbacks.  It assumes that the physical
    // domain size of the current problem is the same as that of the
    // one that generated the pltfile.  It also assumes that the pltfile
    // has at least as many levels as does the current problem.  If
    // either of these are false this code is likely to core dump.
    //

    std::string pltfile;
    pp.query("pltfile", pltfile);
    if (pltfile.empty())
        BoxLib::Abort("You must specify `pltfile'");
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "initData: reading data from: " << pltfile << '\n';

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(pltfile, fileType);

    if (!dataServices.AmrDataOk())
        //
        // This calls ParallelDescriptor::EndParallel() and exit()
        //
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    
    AmrData&                  amrData     = dataServices.AmrDataRef();
    const int                 nspecies    = getChemSolve().numSpecies();
    const Array<std::string>& names       = getChemSolve().speciesNames();   
    Array<std::string>        plotnames   = amrData.PlotVarNames();

    // Determine how species passed in
    bool has_mass, has_mole;
    pp.query("has_mass",has_mass);
    pp.query("has_mole",has_mole);
    std::string specStr;
    if (has_mole) {
        specStr = "X";
    }
    else if (has_mass) {
        specStr = "Y";
    } else {
        BoxLib::Abort("must declare whether restart plotfile has mass or mole fractions");
    }        

    int idT = -1, idX = -1, idSpec = -1;
    const std::string pltSpecName = specStr + "("+names[0]+")";
    for (int i = 0; i < plotnames.size(); ++i)
    {
        if (plotnames[i] == "temp")       idT = i;
        if (plotnames[i] == "x_velocity") idX = i;
        if (plotnames[i] == pltSpecName)  idSpec = i;
    }
    //
    // In the plotfile the mass fractions directly follow the velocities.
    //
    BL_ASSERT(idT>=0 && idX>=0 && idSpec>=0);

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        amrData.FillVar(S_new, level, plotnames[idX+i], Xvel+i);
        amrData.FlushGrids(idX+i);
    }
    amrData.FillVar(S_new, level, plotnames[idT], Temp);
    amrData.FlushGrids(idT);

    for (int i = 0; i < nspecies; i++)
    {
        amrData.FillVar(S_new, level, plotnames[idSpec+i], first_spec+i);
        amrData.FlushGrids(idSpec+i);
    }

    if (has_mole) {

        // assume ok to do in place
        ChemDriver& cd = getChemSolve();
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
            FArrayBox& fab = S_new[mfi];
            const Box& box = mfi.validbox();
            cd.moleFracToMassFrac(fab,fab,box,first_spec,first_spec);
        }
    }


    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "initData: finished init from pltfile" << '\n';
#endif

    for (MFIter snewmfi(S_new); snewmfi.isValid(); ++snewmfi)
    {
        BL_ASSERT(grids[snewmfi.index()] == snewmfi.validbox());

        P_new[snewmfi].setVal(0);

        const int  i       = snewmfi.index();
        RealBox    gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
        const int* lo      = snewmfi.validbox().loVect();
        const int* hi      = snewmfi.validbox().hiVect();
        const int* s_lo    = S_new[snewmfi].loVect();
        const int* s_hi    = S_new[snewmfi].hiVect();
        const int* p_lo    = P_new[snewmfi].loVect();
        const int* p_hi    = P_new[snewmfi].hiVect();

#ifdef BL_USE_NEWMECH
        FORT_INITDATANEWMECH (&level,&cur_time,lo,hi,&ns,
                              S_new[snewmfi].dataPtr(Xvel),
                              S_new[snewmfi].dataPtr(BL_SPACEDIM),
                              ARLIM(s_lo), ARLIM(s_hi),
                              P_new[snewmfi].dataPtr(),
                              ARLIM(p_lo), ARLIM(p_hi),
                              dx,gridloc.lo(),gridloc.hi() );
#else
        FORT_INITDATA (&level,&cur_time,lo,hi,&ns,
                       S_new[snewmfi].dataPtr(Xvel),
                       S_new[snewmfi].dataPtr(BL_SPACEDIM),
                       ARLIM(s_lo), ARLIM(s_hi),
                       P_new[snewmfi].dataPtr(),
                       ARLIM(p_lo), ARLIM(p_hi),
                       dx,gridloc.lo(),gridloc.hi() );
#endif
    }

#ifdef BL_USE_VELOCITY
    //
    // We want to add the velocity from the supplied plotfile
    // to what we already put into the velocity field via FORT_INITDATA.
    //
    // This code has a few drawbacks.  It assumes that the physical
    // domain size of the current problem is the same as that of the
    // one that generated the pltfile.  It also assumes that the pltfile
    // has at least as many levels (with the same refinement ratios) as does
    // the current problem.  If either of these are false this code is
    // likely to core dump.
    //

    std::string velocity_plotfile;
    pp.query("velocity_plotfile", velocity_plotfile);

    if (!velocity_plotfile.empty())
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "initData: reading data from: " << velocity_plotfile << '\n';

        DataServices::SetBatchMode();
        Amrvis::FileType fileType(Amrvis::NEWPLT);
        DataServices dataServices(velocity_plotfile, fileType);

        if (!dataServices.AmrDataOk())
            //
            // This calls ParallelDescriptor::EndParallel() and exit()
            //
            DataServices::Dispatch(DataServices::ExitRequest, NULL);

        AmrData&                  amrData   = dataServices.AmrDataRef();
        Array<std::string>        plotnames = amrData.PlotVarNames();

        if (amrData.FinestLevel() < level)
            BoxLib::Abort("initData: not enough levels in plotfile");

        if (amrData.ProbDomain()[level] != Domain())
            BoxLib::Abort("initData: problem domains do not match");
    
        int idX = -1;
        for (int i = 0; i < plotnames.size(); ++i)
            if (plotnames[i] == "x_velocity") idX = i;

        if (idX == -1)
            BoxLib::Abort("Could not find velocity fields in supplied velocity_plotfile");

        MultiFab tmp(S_new.boxArray(), 1, 0);
        for (int i = 0; i < BL_SPACEDIM; i++)
        {
            amrData.FillVar(tmp, level, plotnames[idX+i], 0);
            for (MFIter mfi(tmp); mfi.isValid(); ++mfi)
                S_new[mfi].plus(tmp[mfi], tmp[mfi].box(), 0, Xvel+i, 1);
            amrData.FlushGrids(idX+i);
        }

        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "initData: finished init from velocity_plotfile" << '\n';
    }
#endif /*BL_USE_VELOCITY*/

    make_rho_prev_time();
    make_rho_curr_time();
    //
    // Initialize other types.
    //
    initDataOtherTypes();
    //
    // Initialize divU and dSdt.
    //
    if (have_divu)
    {
        const Real dt       = 1.0;
        const Real dtin     = -1.0; // Dummy value denotes initialization.
        const Real cur_time = state[Divu_Type].curTime();
        MultiFab&  Divu_new = get_new_data(Divu_Type);

        state[State_Type].setTimeLevel(cur_time,dt,dt);

        calc_divu(cur_time,dtin,Divu_new);

        if (have_dsdt)
            get_new_data(Dsdt_Type).setVal(0);
    }

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point) 
    {
        get_new_data(Dpdt_Type).setVal(0);
    }

    is_first_step_after_regrid = false;
    old_intersect_new          = grids;

#ifdef PARTICLES
    NavierStokesBase::initParticleData();
#endif
}

void
HeatTransfer::initDataOtherTypes ()
{
    const Real cur_time  = state[State_Type].curTime();
    {
        MultiFab rhoh(grids,1,0);
        compute_rhohmix(cur_time,rhoh);
        get_new_data(State_Type).copy(rhoh,0,RhoH,1);
    }
    //
    // Set up diffusivities, viscosities (need for initial divu compute)
    //
    // Assume always variable diffusivity.
    // Assume always variable viscosity.
    //
    calcDiffusivity(cur_time,true);
    //
    // Assume that by now, S_new has "good" data
    //
    get_new_data(Ydot_Type).setVal(0);

    get_new_data(FuncCount_Type).setVal(1);

    setThermoPress(cur_time);
}

void
HeatTransfer::init (AmrLevel& old)
{
    NavierStokesBase::init(old);

    HeatTransfer* oldht    = (HeatTransfer*) &old;
    const Real    cur_time = oldht->state[State_Type].curTime();
    //
    // Get best ydot data.
    //
    MultiFab& Ydot = get_new_data(Ydot_Type);

    for (FillPatchIterator fpi(*oldht,Ydot,Ydot.nGrow(),cur_time,Ydot_Type,0,nspecies);
         fpi.isValid();
         ++fpi)
    {
        Ydot[fpi].copy(fpi());
    }

    RhoH_to_Temp(get_new_data(State_Type));

    MultiFab& FuncCount = get_new_data(FuncCount_Type);

    for (FillPatchIterator fpi(*oldht,FuncCount,FuncCount.nGrow(),cur_time,FuncCount_Type,0,1);
         fpi.isValid();
         ++fpi)
    {
        FuncCount[fpi].copy(fpi());
    }
}

//
// Inits the data on a new level that did not exist before regridding.
//
void
HeatTransfer::init ()
{
    NavierStokesBase::init();
 
    HeatTransfer& old      = getLevel(level-1);
    const Real    cur_time = old.state[State_Type].curTime();
    //
    // Get best ydot data.
    //
    FillCoarsePatch(get_new_data(Ydot_Type),0,cur_time,Ydot_Type,0,nspecies);

    RhoH_to_Temp(get_new_data(State_Type));

    if (new_T_threshold>0)
    {
        MultiFab& crse = old.get_new_data(State_Type);
        MultiFab& fine = get_new_data(State_Type);

        RhoH_to_Temp(crse,0); // Make sure T is current
        Real min_T_crse = crse.min(Temp);
        Real min_T_fine = min_T_crse * std::min(1.0, new_T_threshold);

        const int* ratio = crse_ratio.getVect();
        int Tcomp = (int)Temp;
        int Rcomp = (int)Density;
        int n_tmp = std::max( (int)Density, std::max( (int)Temp, std::max( last_spec, RhoH) ) );
        Array<Real> tmp(n_tmp);
        int num_cells_hacked = 0;
        for (MFIter mfi(fine); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab = fine[mfi];
            const Box& box = mfi.validbox();
            num_cells_hacked += 
                FORT_CONSERVATIVE_T_FLOOR(box.loVect(), box.hiVect(),
                                          fab.dataPtr(), ARLIM(fab.loVect()), ARLIM(fab.hiVect()),
                                          &min_T_fine, &Tcomp, &Rcomp, &first_spec, &last_spec, &RhoH,
                                          ratio, tmp.dataPtr(), &n_tmp);
        }

        ParallelDescriptor::ReduceIntSum(num_cells_hacked);

        if (num_cells_hacked > 0)
        {

            Real old_min = fine.min(Temp);
            RhoH_to_Temp(get_new_data(State_Type));
            Real new_min = fine.min(Temp);

            if (verbose && ParallelDescriptor::IOProcessor())
            {
                std::cout << "...level data adjusted to reduce new extrema ("
                          << num_cells_hacked
                          << " cells affected), new min = "
                          << new_min
                          << " (old min = "
                          << old_min
                          << ")\n";
            }
        }
    }

    FillCoarsePatch(get_new_data(FuncCount_Type),0,cur_time,FuncCount_Type,0,1);
}

void
HeatTransfer::post_timestep (int crse_iteration)
{
    NavierStokesBase::post_timestep(crse_iteration);
    
    if (plot_reactions && level == 0)
    {
        const int Ndiag = auxDiag["REACTIONS"]->nComp();
        //
        // Multiply by the inverse of the coarse timestep.
        //
        const Real factor = 1.0 / crse_dt;

        for (int i = parent->finestLevel(); i >= 0; --i)
            getLevel(i).auxDiag["REACTIONS"]->mult(factor);

        for (int i = parent->finestLevel(); i > 0; --i)
        {
            HeatTransfer& clev = getLevel(i-1);
            HeatTransfer& flev = getLevel(i);

            MultiFab& Ydot_crse = *(clev.auxDiag["REACTIONS"]);
            MultiFab& Ydot_fine = *(flev.auxDiag["REACTIONS"]);

            BoxLib::average_down(Ydot_fine,Ydot_crse,
				 flev.geom, clev.geom,
                                 0,Ndiag,parent->refRatio(i-1));
        }
    }
}
 
void
HeatTransfer::post_restart ()
{
    NavierStokesBase::post_restart();

    Real dummy  = 0;
    int MyProc  = ParallelDescriptor::MyProc();
    int step    = parent->levelSteps(0);
    int restart = 1;

    if (do_active_control)
    {
        int usetemp = 0;
        FORT_ACTIVECONTROL(&dummy,&dummy,&crse_dt,&MyProc,&step,&restart,&usetemp);
    }
    else if (do_active_control_temp)
    {
        int usetemp = 1;
        FORT_ACTIVECONTROL(&dummy,&dummy,&crse_dt,&MyProc,&step,&restart,&usetemp);
    }
}

void
HeatTransfer::postCoarseTimeStep (Real cumtime)
{
    //
    // postCoarseTimeStep() is only called by level 0.
    //
    BL_ASSERT(level == 0);
}

void
HeatTransfer::post_regrid (int lbase,
                           int new_finest)
{
    NavierStokesBase::post_regrid(lbase, new_finest);
    //
    // FIXME: This may be necessary regardless, unless the interpolation
    //        to fine from coarse data preserves rho=sum(rho.Y)
    //
    if (do_set_rho_to_species_sum)
    {
        const int nGrow = 0;
        if (parent->levelSteps(0)>0 && level>lbase)
            set_rho_to_species_sum(get_new_data(State_Type),0,nGrow,0);
    }
}

void
HeatTransfer::checkPoint (const std::string& dir,
                          std::ostream&      os,
                          VisMF::How         how,
                          bool               dump_old)
{
    BL_PROFILE("HeatTransfer::checkPoint()");

    NavierStokesBase::checkPoint(dir,os,how,dump_old);

    if (level == 0)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            const std::string tvfile = dir + "/" + typical_values_filename;
            std::ofstream tvos;
            tvos.open(tvfile.c_str(),std::ios::out|std::ios::trunc|std::ios::binary);
            if (!tvos.good())
                BoxLib::FileOpenFailed(tvfile);
            Box tvbox(IntVect(),(NUM_STATE-1)*BoxLib::BASISV(0));
            int nComp = typical_values.size();
            FArrayBox tvfab(tvbox,nComp);
            for (int i=0; i<nComp; ++i) {
                tvfab.dataPtr()[i] = typical_values[i];
            }
            tvfab.writeOn(tvos);
        }
    }
}

void
HeatTransfer::post_init (Real stop_time)
{
    if (level > 0)
        //
        // Nothing to sync up at level > 0.
        //
        return;

    BL_PROFILE("HeatTransfer::post_init()");

    const Real cur_time     = state[State_Type].curTime();
    const int  finest_level = parent->finestLevel();
    Real        dt_init     = 0.0;

    Array<Real> dt_save(finest_level+1);
    Array<int>  nc_save(finest_level+1);
    Real        dt_init2 = 0.0;
    Array<Real> dt_save2(finest_level+1);
    Array<int>  nc_save2(finest_level+1);
    //
    // Load typical values for each state component
    //
    set_typical_values(false);
    //
    // Ensure state is consistent, i.e. velocity field satisfies initial
    // estimate of constraint, coarse levels are fine level averages, pressure
    // is zero.
    //
    post_init_state();
    //
    // Estimate the initial timestepping.
    //
    post_init_estDT(dt_init,nc_save,dt_save,stop_time);
    //
    // Better estimate needs dt to estimate divu
    //
    const bool do_iter        = do_init_proj && projector;
    const int  init_divu_iter = do_iter ? num_divu_iters : 0;

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "doing num_divu_iters = " << num_divu_iters << '\n';

    for (int iter = 0; iter < init_divu_iter; ++iter)
    {
        //
        // Update species destruction rates in each level but not state.
        //
        if (nspecies > 0)
        {
            for (int k = 0; k <= finest_level; k++)
            {
                MultiFab& S_new = getLevel(k).get_new_data(State_Type);
                //
                // Don't update S_new in this strang_chem() call ...
                //
                MultiFab S_tmp(S_new.boxArray(),S_new.nComp(),0);

                S_tmp.copy(S_new);  // Parallel copy

                getLevel(k).strang_chem(S_tmp,dt_save[k],HT_EstimateYdotNew);
            }
        }
        //
        // Recompute the velocity to obey constraint with chemistry and
        // divqrad and then average that down.
        //
        if (nspecies > 0)
        {
            for (int k = 0; k <= finest_level; k++)
            {
                MultiFab&  Divu_new = getLevel(k).get_new_data(Divu_Type);
                getLevel(k).calc_divu(cur_time,dt_save[k],Divu_new);
            }
            if (!hack_noavgdivu)
            {
                for (int k = finest_level-1; k >= 0; k--)
                {
                    HeatTransfer&   fine_lev = getLevel(k+1);
                    const Geometry& fgeom    = fine_lev.geom;
                    
                    HeatTransfer&   crse_lev = getLevel(k);
                    const Geometry& cgeom    = crse_lev.geom;
                    const IntVect&  fratio   = crse_lev.fine_ratio;
                    
                    MultiFab& Divu_crse = crse_lev.get_new_data(Divu_Type);
                    MultiFab& Divu_fine = fine_lev.get_new_data(Divu_Type);

                    BoxLib::average_down(Divu_fine,Divu_crse,
					 fgeom, cgeom,
                                         0,1,fratio);
                }
            }
            //
            // Recompute the initial velocity field based on this new constraint
            //
            const Real divu_time = state[Divu_Type].curTime();

            int havedivu = 1;

            projector->initialVelocityProject(0,divu_time,havedivu);
            //
            // Average down the new velocity
            //
            for (int k = finest_level-1; k >= 0; k--)
            {
                const Geometry& fgeom   = getLevel(k+1).geom;
                const Geometry& cgeom   = getLevel(k  ).geom;
                MultiFab&       S_fine  = getLevel(k+1).get_new_data(State_Type);
                MultiFab&       S_crse  = getLevel(k  ).get_new_data(State_Type);
                IntVect&        fratio  = getLevel(k  ).fine_ratio;
                
                BoxLib::average_down(S_fine,S_crse,
				     fgeom, cgeom,
                                     Xvel,BL_SPACEDIM,fratio);
            }
        }
        //
        // Estimate the initial timestepping again, using new velocity
        // (necessary?) Need to pass space to save dt, nc, but these are
        // hacked, just pass something.
        //
	// reset Ncycle to nref...
        //
	parent->setNCycle(nc_save);
        post_init_estDT(dt_init2, nc_save2, dt_save2, stop_time);
	//
	// Compute dt_init,dt_save as the minimum of the values computed
	// in the calls to post_init_estDT
	// Then setTimeLevel and dt_level to these values.
	//
	dt_init = std::min(dt_init,dt_init2);
	Array<Real> dt_level(finest_level+1,dt_init);

	parent->setDtLevel(dt_level);
	for (int k = 0; k <= finest_level; k++)
        {
	    dt_save[k] = std::min(dt_save[k],dt_save2[k]);
	    getLevel(k).setTimeLevel(cur_time,dt_init,dt_init);
        }
    }
    //
    // Initialize the pressure by iterating the initial timestep.
    //
    post_init_press(dt_init, nc_save, dt_save);
    //
    // Compute the initial estimate of conservation.
    //
    if (sum_interval > 0)
        sum_integrated_quantities();

}

void
HeatTransfer::sum_integrated_quantities ()
{
    const int finest_level = parent->finestLevel();
    const Real time        = state[State_Type].curTime();

    Real mass = 0.0;
    for (int lev = 0; lev <= finest_level; lev++)
	mass += getLevel(lev).volWgtSum("density",time);

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "TIME= " << time << " MASS= " << mass;

    if (getChemSolve().index(fuelName) >= 0)
    {
        int MyProc  = ParallelDescriptor::MyProc();
        int step    = parent->levelSteps(0);
        int restart = 0;

        if (do_active_control)
        {
            Real fuelmass = 0.0;
            std::string fuel = "rho.Y(" + fuelName + ")";
            for (int lev = 0; lev <= finest_level; lev++)
                fuelmass += getLevel(lev).volWgtSum(fuel,time);

            if (verbose && ParallelDescriptor::IOProcessor())
                std::cout << " FUELMASS= " << fuelmass;

            int usetemp = 0;

            FORT_ACTIVECONTROL(&fuelmass,&time,&crse_dt,&MyProc,&step,&restart,&usetemp);
        }
        else if (do_active_control_temp)
        {
            const int   DM     = BL_SPACEDIM-1;
            const Real* dx     = geom.CellSize();
            const Real* problo = geom.ProbLo();
            Real        hival  = geom.ProbHi(DM);
            const Real  time   = state[State_Type].curTime();
            MultiFab&   mf     = get_new_data(State_Type);

            for (FillPatchIterator Tfpi(*this,mf,1,time,State_Type,Temp,1);
                 Tfpi.isValid();
                 ++Tfpi)
            {
                const FArrayBox& fab  = Tfpi();
                const Box&       vbox = Tfpi.validbox();

                for (IntVect iv = vbox.smallEnd(); iv <= vbox.bigEnd(); vbox.next(iv))
                {
                    const Real T_hi = fab(iv);

                    if (T_hi > temp_control)
                    {
                        const Real hi = problo[DM] + (iv[DM] + .5) * dx[DM];

                        if (hi < hival)
                        {
                            hival = hi;

                            IntVect lo_iv = iv;

                            lo_iv[DM] -= 1;

                            const Real T_lo = fab(lo_iv);

                            if (T_lo < temp_control)
                            {
                                const Real lo    = problo[DM] + (lo_iv[DM] + .5) * dx[DM];
                                const Real slope = (T_hi - T_lo) / (hi - lo);

                                hival = (temp_control - T_lo) / slope + lo;
                            }
                        }
                    }
                }
            }

            ParallelDescriptor::ReduceRealMin(hival);

            if (verbose && ParallelDescriptor::IOProcessor())
                std::cout << " HIVAL= " << hival;

            int usetemp = 1;

            FORT_ACTIVECONTROL(&hival,&time,&crse_dt,&MyProc,&step,&restart,&usetemp);
        }
    }

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << '\n';

    Real rho_h    = 0.0;
    Real rho_temp = 0.0;

    for (int lev = 0; lev <= finest_level; lev++)
    {
        rho_temp += getLevel(lev).volWgtSum("rho_temp",time);
        rho_h     = getLevel(lev).volWgtSum("rhoh",time);
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "TIME= " << time << " LEV= " << lev << " RHOH= " << rho_h << '\n';
    }
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "TIME= " << time << " RHO*TEMP= " << rho_temp << '\n';

    if (verbose) {

        int old_prec = std::cout.precision(15);
        int nc=BL_SPACEDIM+1;
        Array<Real> vmin(nc), vmax(nc);
        Array<int> comps(nc);
        for (int n=0; n<BL_SPACEDIM; ++n) {
            comps[n] = Xvel+n;
        }
        comps[nc-1] = Temp;
        
        for (int lev = 0; lev <= finest_level; lev++)
        {
            MultiFab& newstate = getLevel(lev).get_data(State_Type,time);
            for (int n=0; n<=BL_SPACEDIM; ++n) {
                Real thismin = newstate.min(comps[n],0);
                Real thismax = newstate.max(comps[n],0);

                if (lev==0) {
                    vmin[n] = thismin;
                    vmax[n] = thismax;
                } else {
                    vmin[n] = std::min(thismin,vmin[n]);
                    vmax[n] = std::max(thismax,vmax[n]);
                }
            }
        }

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "min,max temp = " << vmin[nc-1] << ", " << vmax[nc-1] << '\n';
            for (int n=0; n<BL_SPACEDIM; ++n) {
                std::string str = (n==0 ? "xvel" : (n==1 ? "yvel" : "zvel") );
                std::cout << "min,max "<< str << "  = " << vmin[Xvel+n] << ", " << vmax[Xvel+n] << '\n';
            }
        }

        if (nspecies > 0)
        {
            Real min_sum, max_sum;
            for (int lev = 0; lev <= finest_level; lev++) {
                MultiFab* mf = getLevel(lev).derive("rhominsumrhoY",time,0);
                Real this_min = mf->min(0);
                Real this_max = mf->max(0);
                if (lev==0) {
                    min_sum = this_min;
                    max_sum = this_max;
                } else {
                    min_sum = std::min(this_min,min_sum);
                    max_sum = std::max(this_max,max_sum);
                }
                delete mf;
            }
            
            if (ParallelDescriptor::IOProcessor()) {
                std::cout << "min,max rho-sum rho Y_l = "
                          << min_sum << ", " << max_sum << '\n';
            }
            
            for (int lev = 0; lev <= finest_level; lev++)
            {
                MultiFab* mf = getLevel(lev).derive("sumYdot",time,0);
                Real this_min = mf->min(0);
                Real this_max = mf->max(0);
                if (lev==0) {
                    min_sum = this_min;
                    max_sum = this_max;
                } else {
                    min_sum = std::min(this_min,min_sum);
                    max_sum = std::max(this_max,max_sum);
                }
                delete mf;
            }
            
            if (ParallelDescriptor::IOProcessor()) {
                std::cout << "min,max sum Ydot = "
                          << min_sum << ", " << max_sum << '\n';
            }
        }       

        std::cout << std::setprecision(old_prec);
    }
}
    
void
HeatTransfer::post_init_press (Real&        dt_init,
			       Array<int>&  nc_save,
			       Array<Real>& dt_save)
{
    BL_PROFILE("HeatTransfer::post_init_press()");

    const int  nState          = desc_lst[State_Type].nComp();
    const int  nGrow           = 0;
    const Real cur_time        = state[State_Type].curTime();
    const int  finest_level    = parent->finestLevel();
    NavierStokesBase::initial_iter = true;
    //
    // Make space to save a copy of the initial State_Type state data
    //
    PArray<MultiFab> saved_state(finest_level+1,PArrayManage);

    if (init_iter > 0)
        for (int k = 0; k <= finest_level; k++)
            saved_state.set(k,new MultiFab(getLevel(k).grids,nState,nGrow));
    //
    // Iterate over the advance function.
    //
    for (int iter = 0; iter < init_iter; iter++)
    {
	//
	// Squirrel away copy of pre-advance State_Type state
	//
        for (int k = 0; k <= finest_level; k++)
	    MultiFab::Copy(saved_state[k],
                           getLevel(k).get_new_data(State_Type),
			   0,
                           0,
                           nState,
                           nGrow);

        for (int k = 0; k <= finest_level; k++ )
        {
            getLevel(k).advance(cur_time,dt_init,1,1);
        }
        //
        // This constructs a guess at P, also sets p_old == p_new.
        //
        MultiFab** sig = new MultiFab*[finest_level+1];

        for (int k = 0; k <= finest_level; k++)
        {
            sig[k] = &(getLevel(k).get_rho_half_time());
        }

        if (projector)
        {
            int havedivu = 1;
            projector->initialSyncProject(0,sig,parent->dtLevel(0),cur_time,
                                          havedivu);
        }
        delete [] sig;

        for (int k = finest_level-1; k>= 0; k--)
        {
            getLevel(k).avgDown();
        }
        for (int k = 0; k <= finest_level; k++)
        {
            //
            // Reset state variables to initial time, but
            // do not pressure variable, only pressure time.
            //
            getLevel(k).resetState(cur_time, dt_init, dt_init);
        }
	//
	// For State_Type state, restore state we saved above
	//
        for (int k = 0; k <= finest_level; k++)
	    MultiFab::Copy(getLevel(k).get_new_data(State_Type),
                           saved_state[k],
			   0,
                           0,
                           nState,
                           nGrow);

        NavierStokesBase::initial_iter = false;
    }

    if (init_iter <= 0)
        NavierStokesBase::initial_iter = false; // Just being compulsive -- rbp.

    NavierStokesBase::initial_step = false;
    //
    // Re-instate timestep.
    //
    for (int k = 0; k <= finest_level; k++)
    {
        getLevel(k).setTimeLevel(cur_time,dt_save[k],dt_save[k]);
        if (getLevel(k).state[Press_Type].descriptor()->timeType() ==
            StateDescriptor::Point)
          getLevel(k).state[Press_Type].setNewTimeLevel(cur_time+.5*dt_init);
    }
    parent->setDtLevel(dt_save);
    parent->setNCycle(nc_save);
}

//
// Reset the time levels to time (time) and timestep dt.
// This is done at the start of the timestep in the pressure
// iteration section.
//

void
HeatTransfer::resetState (Real time,
                          Real dt_old,
                          Real dt_new)
{
    NavierStokesBase::resetState(time,dt_old,dt_new);

    state[Ydot_Type].reset();
    state[Ydot_Type].setTimeLevel(time,dt_old,dt_new);

    state[FuncCount_Type].reset();
    state[FuncCount_Type].setTimeLevel(time,dt_old,dt_new);
}

void
HeatTransfer::avgDown ()
{
    if (level == parent->finestLevel()) return;

    HeatTransfer&   fine_lev = getLevel(level+1);
    const Geometry& fgeom    = fine_lev.geom;
    MultiFab&       S_crse   = get_new_data(State_Type);
    MultiFab&       S_fine   = fine_lev.get_new_data(State_Type);

    BoxLib::average_down(S_fine,S_crse,
			 fgeom, geom,
                         0,S_crse.nComp(),fine_ratio);
    //
    // Fill rho_ctime at the current and finer levels with the correct data.
    //
    for (int lev = level; lev <= parent->finestLevel(); lev++)
    {
        getLevel(lev).make_rho_curr_time();
    }
    //
    // Reset the temperature
    //
    RhoH_to_Temp(S_crse);
    //
    // Now average down pressure over time n-(n+1) interval.
    //
    MultiFab&       P_crse      = get_new_data(Press_Type);
    MultiFab&       P_fine_init = fine_lev.get_new_data(Press_Type);
    MultiFab&       P_fine_avg  = fine_lev.p_avg;
    MultiFab&       P_fine      = initial_step ? P_fine_init : P_fine_avg;
    const BoxArray& P_fgrids    = fine_lev.state[Press_Type].boxArray();

    BoxArray crse_P_fine_BA = P_fgrids;

    crse_P_fine_BA.coarsen(fine_ratio);

    {
        MultiFab crse_P_fine(crse_P_fine_BA,1,0);

        for (MFIter mfi(P_fine); mfi.isValid(); ++mfi)
        {
            const int i = mfi.index();

            injectDown(crse_P_fine_BA[i],crse_P_fine[mfi],P_fine[mfi],fine_ratio);
        }

        P_crse.copy(crse_P_fine);  // Parallel copy
	const Geometry& cgeom = parent->Geom(level);
	cgeom.PeriodicCopy(P_crse, crse_P_fine);
    }
    //
    // Next average down divu and dSdT at new time.
    //
    if (hack_noavgdivu) 
    {
        //
        // Now that state averaged down, recompute divu (don't avgDown,
        // since that will give a very different value, and screw up mac)
        //
        StateData& divuSD = get_state_data(Divu_Type);// should be const...
        const Real time   = divuSD.curTime();
        const Real dt     = time - divuSD.prevTime();
        calc_divu(time,dt,divuSD.newData());
    }
    else
    {
        MultiFab& Divu_crse = get_new_data(Divu_Type);
        MultiFab& Divu_fine = fine_lev.get_new_data(Divu_Type);
            
        BoxLib::average_down(Divu_fine,Divu_crse,
                             parent->Geom(level+1),parent->Geom(level),
                             0,1,fine_ratio);
    }

    if (hack_noavgdivu)
    {
        //
        //  Now that have good divu, recompute dsdt (don't avgDown,
        //   since that will give a very different value, and screw up mac)
        //
        StateData& dsdtSD = get_state_data(Dsdt_Type);// should be const...
        MultiFab&  dsdt   = dsdtSD.newData();

        if (get_state_data(Divu_Type).hasOldData())
        {
            const Real time = dsdtSD.curTime();
            const Real dt   = time - dsdtSD.prevTime();
            calc_dsdt(time,dt,dsdt);
        }
        else
        {
            dsdt.setVal(0.0);
        }
    }
    else
    {
        MultiFab& Dsdt_crse = get_new_data(Dsdt_Type);
        MultiFab& Dsdt_fine = fine_lev.get_new_data(Dsdt_Type);
            
        BoxLib::average_down(Dsdt_fine,Dsdt_crse,
                             parent->Geom(level+1),parent->Geom(level),
                             0,1,fine_ratio);
    }
}

void
HeatTransfer::scalar_diffusion_update (Real dt,
                                       int  first_scalar, 
                                       int  last_scalar,
                                       int  corrector)
{
    BL_PROFILE("HeatTransfer::scalar_diffusion_update()");
    //
    // Build single component edge-centered array of MultiFabs for fluxes
    //
    FluxBoxes fb_fluxSCn  (this);
    FluxBoxes fb_fluxSCnp1(this);
    MultiFab** fluxSCn   =   fb_fluxSCn.get();
    MultiFab** fluxSCnp1 = fb_fluxSCnp1.get();
    //
    // Set diffusion solve mode.
    //
    Diffusion::SolveMode solve_mode = Diffusion::ONEPASS;
    //
    // Do implicit c-n solve for each scalar.
    //
    const MultiFab& RhoHalftime = get_rho_half_time();

    for (int sigma = first_scalar; sigma <= last_scalar; sigma++)
    {
        if (is_diffusive[sigma])
        {
            int rho_flag = 0;
	    
            MultiFab* delta_rhs = 0;
	    MultiFab* alpha = 0;
	    FluxBoxes fb_betan, fb_betanp1;

            diffuse_scalar_setup(dt,sigma,rho_flag,delta_rhs,alpha,fb_betan,fb_betanp1);
            //
            // Note: dt taken care of in diffuse_scalar.
            //
	    const int dataComp = 0; // Component of dR, alpha, betas to use.

            diffusion->diffuse_scalar(dt,sigma,be_cn_theta,RhoHalftime,rho_flag,fluxSCn,
                                      fluxSCnp1,0,delta_rhs,dataComp,alpha,dataComp,
                                      fb_betan.get(),fb_betanp1.get(),dataComp,solve_mode);
	    //
	    // If corrector, increment the viscous flux registers.
            // Assume corrector called only ONCE!.
	    //
	    if (do_reflux && corrector)
	    {
                for (int d = 0; d < BL_SPACEDIM; d++)
                {
                    for (MFIter fmfi(*fluxSCn[d]); fmfi.isValid(); ++fmfi)
                    {
                        (*fluxSCnp1[d])[fmfi].plus((*fluxSCn[d])[fmfi]);

                        if (level > 0)
                            getViscFluxReg().FineAdd((*fluxSCnp1[d])[fmfi],d,fmfi.index(),0,sigma,1,dt);
                    }
                    if (level < parent->finestLevel())
                    {
                        getLevel(level+1).getViscFluxReg().CrseInit(*fluxSCnp1[d],d,0,sigma,1,-dt);
                    }
                }
	    }

	    delete delta_rhs;
	    delete alpha;
        }  
    }    
}

void
HeatTransfer::differential_spec_diffusion_update (Real dt,
						  int  corrector)
{
    BL_PROFILE("HeatTransfer::differential_spec_diffusion_update()");

    if (verbose && benchmarking) ParallelDescriptor::Barrier();

    const Real strt_time = ParallelDescriptor::second();

    if (hack_nospecdiff)
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "... HACK!!! skipping spec diffusion " << '\n';

        if (!corrector)
            MultiFab::Copy(get_new_data(State_Type),get_old_data(State_Type),first_spec,first_spec,nspecies,0);

        for (int d = 0; d < BL_SPACEDIM; ++d)
        {
            SpecDiffusionFluxn[d]->setVal(0,0,nspecies);
            SpecDiffusionFluxnp1[d]->setVal(0,0,nspecies);
        }
        for (int comp = 0; comp < nspecies; ++comp)
            spec_diffusion_flux_computed[comp] = HT_Diffusion;

        return;
    }
    //
    // Build single component edge-centered array of MultiFabs for fluxes
    //
    MultiFab** fluxSCn;
    MultiFab** fluxSCnp1;
    const int nGrow   = 0;
    const int nCompSC = 1;
    const int sCompSC = 0;
    FluxBoxes fb_fluxSCn  (this, nCompSC, nGrow);
    FluxBoxes fb_fluxSCnp1(this, nCompSC, nGrow);
    fluxSCn   = fb_fluxSCn.get();
    fluxSCnp1 = fb_fluxSCnp1.get();
    //
    // Set diffusion solve mode
    //
    Diffusion::SolveMode solve_mode = Diffusion::PREDICTOR;
    //
    // Do implicit c-n solve for each scalar...but dont reflux.
    // Save the fluxes, coeffs and source term, we need 'em for later
    //
    const int sCompY = first_spec;
    const int nCompY = nspecies;
    MultiFab* alpha = 0; // Allocate lazily
    MultiFab delta_rhs(grids, nCompY, nGrow);
    FluxBoxes fb_betan  (this, nCompY, nGrow);
    FluxBoxes fb_betanp1(this, nCompY, nGrow);
    MultiFab **betan   = fb_betan.get();
    MultiFab **betanp1 = fb_betanp1.get();
    Array<int> rho_flag(nCompY,0);
    
    const MultiFab& RhoHalftime = get_rho_half_time();

    for (int sigma = 0; sigma < nCompY; ++sigma)
    {
	const int state_ind = sCompY + sigma;

	MultiFab* alphaSC;
	MultiFab* delta_rhsSC;
	FluxBoxes fb_betanSC, fb_betanp1SC;

	// delta_rhsSC is zero coming out of this
	diffuse_scalar_setup(dt, state_ind, rho_flag[sigma], delta_rhsSC,
			     alphaSC, fb_betanSC, fb_betanp1SC);
	
	MultiFab** betanSC = fb_betanSC.get();
	MultiFab** betanp1SC = fb_betanp1SC.get();

	if (state_ind == sCompY && alphaSC)
	{
	    alpha = new MultiFab(grids, nCompY, nGrow);
	}
	else
	{
	    if ((!alphaSC) ^ !alpha)
		BoxLib::Error("All diff-diffusion must be of same form");
	}   
	//    
	// Nab a copy of the coeffs, call diffuser
	//
	if (alphaSC)
	    MultiFab::Copy(*alpha,*alphaSC,sCompSC,sigma,nCompSC,nGrow);
	
	for (int d=0; d<BL_SPACEDIM; ++d)
	{
	    if (betanSC)
		MultiFab::Copy(  *betan[d],  *betanSC[d],sCompSC,sigma,nCompSC,nGrow);
	    if (betanp1SC)
		MultiFab::Copy(*betanp1[d],*betanp1SC[d],sCompSC,sigma,nCompSC,nGrow);
	}
	//
	// Make sure we've got a place for delta_rhs...predictor will dump any
	// explicit updates taken before this routine into delta_rhs for later.
	//
	if (delta_rhsSC)
	{
	    MultiFab::Copy(delta_rhs,*delta_rhsSC,sCompSC,sigma,nCompSC,nGrow);
	}
	else
	{
	    delta_rhs.setVal(0,sigma,1);
	}

	delete delta_rhsSC;
	delete alphaSC;

	diffusion->diffuse_scalar(dt,state_ind,be_cn_theta,RhoHalftime,rho_flag[sigma],
                                  fluxSCn,fluxSCnp1,0,&delta_rhs,sigma,alpha,sigma,
                                  betan,betanp1,sigma,solve_mode);
	//
	// Pull fluxes into flux array
	//
	for (int d = 0; d < BL_SPACEDIM; ++d)
        {
	    MultiFab::Copy(*SpecDiffusionFluxn[d],  *fluxSCn[d],  sCompSC,sigma,nCompSC,nGrow);
	    MultiFab::Copy(*SpecDiffusionFluxnp1[d],*fluxSCnp1[d],sCompSC,sigma,nCompSC,nGrow);
        }
	spec_diffusion_flux_computed[sigma] = HT_Diffusion;
    }
    fb_fluxSCn.clear();
    fb_fluxSCnp1.clear();
    //
    // Modify update/fluxes to preserve flux sum = 0, compute new update and
    // leave modified fluxes in level data.  Do this in two stages, first for
    // the explicit fluxes, then the implicit ones (so send in rhs=0 for the
    // second one...).
    //
    const int  dataComp  = 0;
    const Real prev_time = state[State_Type].prevTime();
    const Real cur_time  = state[State_Type].curTime();

    // conservatively correct fluxes at time n, then set new = old + (dt/2)*fluxes^old
    // the 1/2 was added in differential_spec_diffusion_update in compflux (b)
    adjust_spec_diffusion_update(get_new_data(State_Type),&get_old_data(State_Type),
				 sCompY,dt,prev_time,rho_flag,RhoHalftime,dataComp,
                                 &delta_rhs,alpha,betan);

    fb_betan.clear();
    delta_rhs.clear();

    // conservatively correct fluxes at time n+1, then set new = new + (dt/2)*fluxes^new
    // the 1/2 was added in differential_spec_diffusion_update in compflux (b)
    adjust_spec_diffusion_update(get_new_data(State_Type),&get_new_data(State_Type),
				 sCompY,dt,cur_time,rho_flag,RhoHalftime,dataComp,0,
                                 alpha,betanp1);

    fb_betanp1.clear();
    if (alpha)
	delete alpha;
    //
    // Now do reflux with new, improved fluxes
    //
    if (do_reflux && corrector)
    {
        FArrayBox fluxtot;

        for (int d = 0; d < BL_SPACEDIM; d++)
        {
            MultiFab fluxes;

            if (level < parent->finestLevel())
                fluxes.define(SpecDiffusionFluxn[d]->boxArray(), nCompY, 0, Fab_allocate);

            for (MFIter fmfi(*SpecDiffusionFluxn[d]); fmfi.isValid(); ++fmfi)
            {
                const Box& ebox = (*SpecDiffusionFluxn[d])[fmfi].box();

                fluxtot.resize(ebox,nCompY);
                fluxtot.copy((*SpecDiffusionFluxn[d])[fmfi], ebox,0,ebox,0,nCompY);
                fluxtot.plus((*SpecDiffusionFluxnp1[d])[fmfi],ebox,0,0,nCompY);

                if (level < parent->finestLevel())
                    fluxes[fmfi].copy(fluxtot);

                if (level > 0)
                    getViscFluxReg().FineAdd(fluxtot,d,fmfi.index(),0,sCompY,nCompY,dt);
            }

            if (level < parent->finestLevel())
            {
                getLevel(level+1).getViscFluxReg().CrseInit(fluxes,d,0,sCompY,nCompY,-dt);
            }
        }
    }

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::differential_spec_diffusion_update(): time: " << run_time << '\n';
    }
}

void
HeatTransfer::adjust_spec_diffusion_fluxes (Real time)
{
    //
    // In this function we explicitly adjust the species diffusion fluxes so that their sum
    // is zero.  The fluxes are class member data, either SpecDiffusionFluxn or 
    // SpecDiffusionFluxnp1, depending on time
    //
    const TimeLevel whichTime = which_time(State_Type,time);
    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);    
    MultiFab* const * flux = (whichTime == AmrOldTime) ? SpecDiffusionFluxn : SpecDiffusionFluxnp1;

    MultiFab& S = get_data(State_Type,time);
    //
    // Fill grow cells in the state for RhoY, Temp
    // For Dirichlet physical boundary grow cells, this data will live on the cell face, otherwise
    // it will live at cell centers.
    //
    const int nGrow = 1;
    BL_ASSERT(S.nGrow()>=nGrow);
    for (FillPatchIterator Tfpi(*this,S,nGrow,time,State_Type,Temp,1),
             Yfpi(*this,S,nGrow,time,State_Type,first_spec,nspecies);
         Yfpi.isValid() && Tfpi.isValid();
         ++Yfpi, ++Tfpi)
    {
        const Box& vbox = Yfpi.validbox();
        FArrayBox& fab = S[Yfpi];
        BoxList gcells = BoxLib::boxDiff(Box(vbox).grow(nGrow),vbox);
        for (BoxList::const_iterator it = gcells.begin(), end = gcells.end(); it != end; ++it)
        {
            const Box& gbox = *it;
            fab.copy(Tfpi(),gbox,0,gbox,Temp,1);
            fab.copy(Yfpi(),gbox,0,gbox,first_spec,nspecies);
        }
    }
    //
    // Get boundary info for Y (assume all Ys have the same boundary type.
    //
    const BCRec& Ybc = get_desc_lst()[State_Type].getBC(first_spec);
    // 
    // The following REPAIR_FLUX routine modifies the fluxes of all the species
    // to ensure that they sum to zero.  It requires the RhoY on valid + 1 grow (and, at least
    // as of this writing actually used the values of RhoY on edges that it gets by arithmetic
    // averaging.
    //
    const Box& domain = geom.Domain();
    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        const Box& box   = mfi.validbox();
        FArrayBox& Y = S[mfi];
        int sCompY=first_spec;

        for (int d =0; d < BL_SPACEDIM; ++d)
        {
            FArrayBox& f = (*flux[d])[mfi];
            FORT_REPAIR_FLUX(box.loVect(), box.hiVect(), domain.loVect(), domain.hiVect(),
                             f.dataPtr(),      ARLIM(f.loVect()),ARLIM(f.hiVect()),
                             Y.dataPtr(sCompY),ARLIM(Y.loVect()),ARLIM(Y.hiVect()),
                             &d, Ybc.vect());
        }
    }
}

void
HeatTransfer::adjust_spec_diffusion_update (MultiFab&              Phi_new,
					    const MultiFab*        Phi_old,
					    int                    sCompS,
					    Real                   dt,
					    Real                   time,
					    const Array<int>&      rho_flag,
					    const MultiFab&        rho_half,
					    int                    dataComp,
					    const MultiFab*        delta_rhs, 
					    const MultiFab*        alpha, 
					    const MultiFab* const* betanp1)
{
    BL_PROFILE("HeatTransfer::adjust_spec_diffusion_update()");
    //
    // Here, we're going to compute an update using fluxes computed from
    // the old state and a guess for the new one.  These fluxes are modified
    // arbitrarily however to preserve that the sum over all species of the
    // diffusive fluxes is zero.  We do this by marching though each edge 
    // and adjusting the flux of the downstream state that is the "dominant"
    // species.  Improved fluxes are left in SpecDiffusionFlux.
    //
    BL_ASSERT(!alpha ^ (alpha&&alpha->nComp()  == nspecies));
    BL_ASSERT(betanp1 && betanp1[0]->nComp()   == nspecies);
    BL_ASSERT(!delta_rhs || delta_rhs->nComp() == nspecies);

    const TimeLevel whichTime = which_time(State_Type,time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
    
    MultiFab* const * flux = (whichTime == AmrOldTime) ? SpecDiffusionFluxn : SpecDiffusionFluxnp1;

    adjust_spec_diffusion_fluxes(time);
    //
    // Reset Phi_new using "repaired" fluxes.  Here, we assume the following
    // arrangement of terms:
    //
    //     Phi_new = Phi_old + (dt/A)(RHS-Div(flux))
    //            A = rhonph.alpha for rho_flag == 1
    //            A = 1            for rho_flag == 2
    //     
    //
    FArrayBox update,  tmp;

    for (MFIter mfi(Phi_new); mfi.isValid(); ++mfi)
    {
	const Box& box   = mfi.validbox();
	//
	// Do same trick as in Diffusion class, so as to deal with r-z correctly
	// (i.e. multiply though by cell-centered r.dx.dz [integrating over vol].
	//  Since fluxes are already extensive, merely mult flux by dx, then scale
	//  flux divergence by vol...that should do it)
	//
	// Here then, the update becomes: (dt/A)[RHS - (1/V)[Div(flux.dx)]]
	//
	// In RECOMP_UPDATE, we compute (-1/V).Div(flux.dx)
        //
        // Also, by putting update into a separate fab, and copying into afterward
        // we allow that Phi_old = &Phi_new, i.e. that we are updating an existing
        // state
	//
        update.resize(box,nspecies);

	// update = -div(flux)/vol
	FORT_RECOMP_UPDATE(box.loVect(), box.hiVect(),
			   update.dataPtr(),
			   ARLIM(update.loVect()),       ARLIM(update.hiVect()),
			   (*flux[0])[mfi].dataPtr(),
                           ARLIM((*flux[0])[mfi].loVect()),
                           ARLIM((*flux[0])[mfi].hiVect()),
			   (*flux[1])[mfi].dataPtr(),
			   ARLIM((*flux[1])[mfi].loVect()),
                           ARLIM((*flux[1])[mfi].hiVect()),
#if BL_SPACEDIM == 3
			   (*flux[2])[mfi].dataPtr(),
			   ARLIM((*flux[2])[mfi].loVect()),
                           ARLIM((*flux[2])[mfi].hiVect()),
#endif
			   volume[mfi].dataPtr(),
			   ARLIM(volume[mfi].loVect()),ARLIM(volume[mfi].hiVect()),
			   &nspecies);
	
	if (delta_rhs)
	    update.plus((*delta_rhs)[mfi],box,dataComp,0,nspecies);

	update.mult(dt,box,0,nspecies);

	if (alpha)
	    update.divide((*alpha)[mfi],box,dataComp,0,nspecies);

        tmp.resize(box,1);
        tmp.copy(rho_half[mfi],0,0,1);
        tmp.invert(1);

	for (int i = 0; i < nspecies; ++i)
	    if (rho_flag[i] == 1)
		update.mult(tmp,0,i,1);

	if (Phi_old)
	    update.plus((*Phi_old)[mfi],box,sCompS,0,nspecies);

        Phi_new[mfi].copy(update,box,0,box,sCompS,nspecies);
    }
}

void
HeatTransfer::diffuse_scalar_setup (Real        dt,
                                    int         sigma,
                                    int&        rho_flag, 
                                    MultiFab*&  delta_rhs,
                                    MultiFab*&  alpha, 
                                    FluxBoxes&  fb_betan,
                                    FluxBoxes&  fb_betanp1)
{
    //
    // Do setup for implicit c-n solve for an arbitrary scalar.
    //
    // Note: should be ok for variable c_p.
    //
    const Real prev_time = state[State_Type].prevTime();

    NavierStokesBase::diffuse_scalar_setup(sigma, rho_flag);

    alpha     = 0;
    delta_rhs = 0;
   
    if (sigma == Temp)
    {
        rho_flag = 1;
    }
    else if (sigma == RhoH || (sigma >= first_spec && sigma <= last_spec))
    {
        rho_flag = 2;
    }

    if (sigma == RhoH)
    {
        diffuse_rhoh_setup(prev_time,dt,delta_rhs); 
    }
    else if (sigma == Temp)
    {
        diffuse_temp_setup(prev_time,dt,delta_rhs,alpha); 
    }
    else if (sigma >= first_spec && sigma <= last_spec)
    {
        diffuse_spec_setup(sigma,prev_time,dt,delta_rhs); 
    }

    MultiFab** betan   = fb_betan.define  (this);
    MultiFab** betanp1 = fb_betanp1.define(this);

    getDiffusivity(betan, prev_time, sigma, 0, 1);
    getDiffusivity(betanp1, prev_time+dt, sigma, 0, 1);
}

void
HeatTransfer::diffuse_spec_setup (int        istate,
                                  Real       time,
                                  Real       dt, 
                                  MultiFab*& delta_rhs)
{
    //
    // Chemistry split, no source terms
    //
    delta_rhs = new MultiFab(grids,1,0);
    delta_rhs->setVal(0);
}

void
HeatTransfer::compute_OT_radloss (Real      time,
                                  int       nGrow,
                                  MultiFab& dqrad)
{
    //
    // Get optically thin radiation losses (+ve when energy LOST).
    //
    BL_ASSERT(do_OT_radiation || do_heat_sink);
    BL_ASSERT(nGrow <= dqrad.nGrow());
    BL_ASSERT(dqrad.boxArray() == grids);

    const Real* dx = geom.CellSize();

    Real p_amb;
    FORT_GETPAMB(&p_amb);

    const Real Patm = p_amb / P1atm_MKS;
    const Real T_bg = 298.0;

    FArrayBox tmp;

    for (FillPatchIterator fpi(*this,dqrad,nGrow,time,State_Type,0,NUM_STATE);
         fpi.isValid();
         ++fpi)
    {
        FArrayBox& S   = fpi();
        const Box& box = S.box();

        tmp.resize(box,1);
        tmp.copy(S,Density,0,1);
        tmp.invert(1);
        for (int spec = first_spec; spec <= last_spec; spec++)
            S.mult(tmp,0,spec,1);

        FArrayBox& dqrad_fab = dqrad[fpi];

        tmp.resize(box,nspecies);

        getChemSolve().massFracToMoleFrac(tmp,S,box,first_spec,0);

        if (do_OT_radiation)
        {
            getChemSolve().getOTradLoss_TDF(dqrad_fab,S,tmp,Patm,T_bg,box,0,Temp,0);
        }
        else
        {
            dqrad_fab.setVal(0.0);
        }
        //
        // Add arbitrary heat sink.
        //
        if (do_heat_sink)
        {
            FArrayBox rad(box,1);
            FORT_RADLOSS(box.loVect(),box.hiVect(),
                         rad.dataPtr(),         ARLIM(rad.loVect()),ARLIM(rad.hiVect()),
                         S.dataPtr(Temp),       ARLIM(S.loVect()),  ARLIM(S.hiVect()),
                         S.dataPtr(first_spec), ARLIM(S.loVect()),  ARLIM(S.hiVect()),
                         dx, &Patm, &time);
            dqrad_fab.plus(rad,0,0,1);
        }
    }
}

void
HeatTransfer::diffuse_rhoh_setup (Real       time,
                                  Real       dt, 
                                  MultiFab*& delta_rhs)
{
    //
    // Do set-up for implicit c-n solve for rho*h using Le=1 equation.
    //
    BL_ASSERT(delta_rhs==0);
    delta_rhs = new MultiFab (grids,1,0);
    delta_rhs->setVal(0);
    const int nGrow = 0;

    if (do_OT_radiation || do_heat_sink)
    {
        MultiFab dqrad(grids,1,nGrow);

        compute_OT_radloss(time,nGrow,dqrad);
        for (MFIter mfi(dqrad); mfi.isValid(); ++mfi)
        {
            dqrad[mfi].mult(1.0-be_cn_theta);
            (*delta_rhs)[mfi].minus(dqrad[mfi]);
        }

        compute_OT_radloss(time+dt,nGrow,dqrad);
        for (MFIter mfi(dqrad); mfi.isValid(); ++mfi)
        {
            dqrad[mfi].mult(be_cn_theta);
            (*delta_rhs)[mfi].minus(dqrad[mfi]);
        }
    }
}

void
HeatTransfer::diffuse_temp_setup (Real       prev_time,
                                  Real       dt, 
                                  MultiFab*& delta_rhs,
                                  MultiFab*& alpha)
{
    //
    // Do set-up for implicit c-n solve for T.
    //
    BL_ASSERT(delta_rhs==0);
    delta_rhs = new MultiFab (grids,1,0);
    delta_rhs->setVal(0);
    const int nGrow = 0;

    if (do_OT_radiation || do_heat_sink)
    {
        MultiFab dqrad(grids,1,nGrow);

        compute_OT_radloss(prev_time,nGrow,dqrad);
        for (MFIter mfi(dqrad); mfi.isValid(); ++mfi)
        {
            dqrad[mfi].mult(1.0-be_cn_theta);
            (*delta_rhs)[mfi].minus(dqrad[mfi]);
        }

        compute_OT_radloss(prev_time+dt,nGrow,dqrad);
        for (MFIter mfi(dqrad); mfi.isValid(); ++mfi)
        {
            dqrad[mfi].mult(be_cn_theta);
            (*delta_rhs)[mfi].minus(dqrad[mfi]);
        }
    }
    //
    // Increment rhs by (+ sum_l rho D grad Y_l dot grad h_l)
    // Note: this way ensures isothermal preservation
    //
    MultiFab rdgydgh(grids,1,0);
    compute_rhoDgradYgradH(prev_time, rdgydgh);
    rdgydgh.mult(1.0-be_cn_theta,0);
    
    for (MFIter mfi(*delta_rhs); mfi.isValid(); ++mfi)
    {
        (*delta_rhs)[mfi].plus(rdgydgh[mfi],0,0,1);
    }
    compute_rhoDgradYgradH(prev_time+dt, rdgydgh);
    rdgydgh.mult(be_cn_theta,0);
    
    for (MFIter mfi(*delta_rhs); mfi.isValid(); ++mfi)
    {
        (*delta_rhs)[mfi].plus(rdgydgh[mfi],0,0,1);
    }
    rdgydgh.clear();
    //
    // alpha = c_p^(n+1/2)
    // Note: rho accounted for in diffusion box
    //
    BL_ASSERT(alpha==0);
    alpha = new MultiFab (grids,1,0);
    compute_cp(prev_time+dt/2,*alpha);
}

void
HeatTransfer::velocity_diffusion_update (Real dt)
{
    BL_PROFILE("HeatTransfer::velocity_diffusion_update()");

    const Real strt_time = ParallelDescriptor::second();
    //
    // Do implicit c-n solve for velocity
    // compute the viscous forcing
    // do following except at initial iteration--rbp, per jbb
    //
    if (is_diffusive[Xvel])
    {
        int rho_flag;
        if (do_mom_diff == 0)
        {
           rho_flag = 1;
        }
        else
        {
           rho_flag = 3;
        }

        MultiFab *delta_rhs = 0;
	FluxBoxes fb_betan, fb_betanp1;

        diffuse_velocity_setup(dt, delta_rhs, fb_betan, fb_betanp1);

        diffusion->diffuse_velocity(dt,be_cn_theta,get_rho_half_time(),rho_flag,
                                    delta_rhs,fb_betan.get(),fb_betanp1.get());

	delete delta_rhs;
    }

    if (verbose)
    {
        Real run_time    = ParallelDescriptor::second() - strt_time;
        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "HeatTransfer::velocity_diffusion_update(): lev: "
                      << level
                      << ", time: " << run_time << '\n';
        }
    }
}
    
void
HeatTransfer::diffuse_velocity_setup (Real        dt,
                                      MultiFab*&  delta_rhs,
                                      FluxBoxes&  fb_betan,
                                      FluxBoxes&  fb_betanp1)
{
    //
    // Do setup for implicit c-n solve for velocity.
    //
    BL_ASSERT(delta_rhs==0);
    const Real time = state[State_Type].prevTime();
    //
    // Assume always variable viscosity.
    //
    MultiFab** betan = fb_betan.define(this);
    MultiFab** betanp1 = fb_betanp1.define(this);

    getViscosity(betan, time);
    getViscosity(betanp1, time+dt);
    //
    // Do the following only if it is a non-trivial operation.
    //
    if (S_in_vel_diffusion)
    {
        //
        // Include div mu S*I terms in rhs
        //  (i.e. make nonzero delta_rhs to add into RHS):
        //
        // The scalar and tensor solvers incorporate the relevant pieces of
        //  of Div(tau), provided the flow is divergence-free.  Howeever, if
        //  Div(U) =/= 0, there is an additional piece not accounted for,
        //  which is of the form A.Div(U).  For constant viscosity, Div(tau)_i
        //  = Lapacian(U_i) + mu/3 d[Div(U)]/dx_i.  For mu not constant,
        //  Div(tau)_i = d[ mu(du_i/dx_j + du_j/dx_i) ]/dx_i - 2mu/3 d[Div(U)]/dx_i
        //
        // As a convenience, we treat this additional term as a "source" in
        // the diffusive solve, computing Div(U) in the "normal" way we
        // always do--via a call to calc_divu.  This routine computes delta_rhs
        // if necessary, and stores it as an auxilliary rhs to the viscous solves.
        // This is a little strange, but probably not bad.
        //
        delta_rhs = new MultiFab(grids,BL_SPACEDIM,0);
        delta_rhs->setVal(0);
      
        MultiFab divmusi(grids,BL_SPACEDIM,0);
        //
	// Assume always variable viscosity.
        //
	diffusion->compute_divmusi(time,betan,divmusi);
	divmusi.mult((-2./3.)*(1.0-be_cn_theta),0,BL_SPACEDIM,0);
        (*delta_rhs).plus(divmusi,0,BL_SPACEDIM,0);
        //
	// Assume always variable viscosity.
        //
	diffusion->compute_divmusi(time+dt,betanp1,divmusi);
	divmusi.mult((-2./3.)*be_cn_theta,0,BL_SPACEDIM,0);
        (*delta_rhs).plus(divmusi,0,BL_SPACEDIM,0);
    }
}

void
HeatTransfer::getViscTerms (MultiFab& visc_terms,
                            int       src_comp, 
                            int       num_comp,
                            Real      time)
{
    BL_PROFILE("HeatTransfer::getViscTerms()");
    //
    // Load "viscous" terms, starting from component = 0.
    //
    // JFG: for species, this procedure returns the *negative* of the divergence of 
    // of the diffusive fluxes.  specifically, in the mixture averaged case, the
    // diffusive flux vector for species k is
    //
    //       j_k = - rho D_k,mix grad Y_k
    //
    // so the divergence of the flux, div dot j_k, has a negative in it.  instead 
    // this procedure returns - div dot j_k to remove the negative.
    //
    // note the fluxes used in the code are extensive, that is, scaled by the areas
    // of the cell edges.  the calculation of the divergence is the sum of un-divided
    // differences of the extensive fluxes, all divided by volume.  so the effect is 
    // to give the true divided difference approximation to the divergence of the
    // intensive flux.
    //
    const int  last_comp = src_comp + num_comp - 1;
    int        icomp     = src_comp; // This is the current related state comp.
    int        load_comp = 0;        // Comp for result of current calculation.
    MultiFab** vel_visc  = 0;        // Potentially reused, raise scope
    FluxBoxes  fb_vel_visc;
    const int  nGrow     = visc_terms.nGrow();
    //
    // Get Div(tau) from the tensor operator, if velocity and have non-const viscosity
    //
    if (src_comp < BL_SPACEDIM)
    {
        if (src_comp != Xvel || num_comp < BL_SPACEDIM)
            BoxLib::Error("tensor v -> getViscTerms needs all v-components at once");

	vel_visc = fb_vel_visc.define(this);
        getViscosity(vel_visc, time);

        showMF("velVT",*viscn_cc,"velVT_viscn_cc",level);
        for (int dir=0; dir<BL_SPACEDIM; ++dir) {
            showMF("velVT",*(vel_visc[dir]),BoxLib::Concatenate("velVT_viscn_",dir,1),level);
        }

        diffusion->getTensorViscTerms(visc_terms,time,vel_visc,0);
        showMF("velVT",visc_terms,"velVT_visc_terms_1",level);

        icomp = load_comp = BL_SPACEDIM;
    }
    //
    // For Le != 1, need to get visc terms in a block
    //
    if (!unity_Le)
    {
	const int non_spec_comps = std::min(num_comp,
                                            std::max(0,first_spec-src_comp) + std::max(0,last_comp-last_spec));
	const bool has_spec = src_comp <= last_spec && last_comp >= first_spec;
	BL_ASSERT( !has_spec || (num_comp>=nspecies+non_spec_comps));
	if (has_spec)
	{
	    const int sCompY = first_spec - src_comp + load_comp;
	    compute_differential_diffusion_terms(visc_terms,sCompY,time);
        }
    }
    //
    // Now, do the rest.
    // Get Div(visc.Grad(state)) or Div(visc.Grad(state/rho)), as appropriate.
    //
    for ( ; icomp <= last_comp; load_comp++, icomp++)
    {
	const bool is_spec  = icomp >= first_spec && icomp <= last_spec;
	if ( !(is_spec && !unity_Le) )
	{
	    if (icomp == Temp)
	    {
		getTempViscTerms(visc_terms,Temp-load_comp,time);
	    }
	    else if (icomp == RhoH)
	    {
		;
	    }
	    else
	    {
		const int  rho_flag = Diffusion::set_rho_flag(diffusionType[icomp]);

		if (icomp == Density)
                {
                    visc_terms.setVal(0.0,load_comp,1);
                }
                else
                {
                    //
                    // Assume always variable viscosity / diffusivity.
                    //
		    FluxBoxes fb_beta(this);
		    MultiFab** beta = fb_beta.get();
                    getDiffusivity(beta, time, icomp, 0, 1);

                    diffusion->getViscTerms(visc_terms,icomp-load_comp,
                                            icomp,time,rho_flag,beta,0);
                }
	    }
	}
    }
    //
    // Add Div(u) term if desired, if this is velocity, and if Div(u) is nonzero
    // If const-visc, term is mu.Div(u)/3, else it's -Div(mu.Div(u).I)*2/3
    //
    if (src_comp < BL_SPACEDIM && S_in_vel_diffusion)
    {
        if (num_comp < BL_SPACEDIM)
            BoxLib::Error("getViscTerms() need all v-components at once");
    
        MultiFab divmusi(grids,BL_SPACEDIM,1);
#ifndef NDEBUG
        showMF("velVT",get_old_data(Divu_Type),"velVT_divu",level);
#endif
        //
	// Assume always using variable viscosity.
        //
	diffusion->compute_divmusi(time,vel_visc,divmusi); // pre-computed visc above
	divmusi.mult((-2./3.),0,BL_SPACEDIM,0);
        showMF("velVT",divmusi,"velVT_divmusi",level);
        visc_terms.plus(divmusi,Xvel,BL_SPACEDIM,0);
        showMF("velVT",visc_terms,"velVT_visc_terms_2",level);
    }
    //
    // Ensure consistent grow cells
    //    
    if (nGrow > 0)
    {
        visc_terms.FillBoundary(0,num_comp, geom.periodicity());
	Extrapolater::FirstOrderExtrap(visc_terms, geom, 0, num_comp);
    }
}


void dumpProfileFab(const FArrayBox& fab,
                    const std::string file)
{
    const Box& box = fab.box();
    int imid = (int)(0.5*(box.smallEnd()[0] + box.bigEnd()[0]));
    IntVect iv1 = box.smallEnd(); iv1[0] = imid;
    IntVect iv2 = box.bigEnd(); iv2[0] = imid;
    iv1[1] = 115; iv2[1]=125;
    Box pb(iv1,iv2);
    std::cout << "dumping over " << pb << std::endl;

    std::ofstream osf(file.c_str());
    for (IntVect iv = pb.smallEnd(); iv <= pb.bigEnd(); pb.next(iv))
    {
        osf << iv[1] << " " << (iv[1]+0.5)*3.5/256 << " ";
        for (int n=0; n<fab.nComp(); ++n) 
            osf << fab(iv,n) << " ";
        osf << std::endl;
    }
    BoxLib::Abort();
}

void dumpProfile(const MultiFab& mf,
                 const std::string file)
{
    const FArrayBox& fab = mf[0];
    dumpProfileFab(fab,file);
}


Real MFnorm(const MultiFab& mf,
            int             sComp,
            int             nComp)
{
    const int p = 0; // max norm
    Real normtot = 0.0;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        normtot = std::max(normtot,mf[mfi].norm(mfi.validbox(),p,sComp,nComp));
    ParallelDescriptor::ReduceRealMax(normtot);
    return normtot;
}

void
HeatTransfer::compute_differential_diffusion_terms (MultiFab& visc_terms,
						    int       sComp,
                                                    Real      time)
{
    BL_PROFILE("HeatTransfer::compute_differential_diffusion_terms()");
    if (hack_nospecdiff)
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "... HACK!!! zeroing spec diffusion terms " << '\n';
        visc_terms.setVal(0.0,sComp,nspecies);
        return;            
    }
    //
    // NOTE: This routine does not fill grow cells
    //
    BL_ASSERT(visc_terms.boxArray() == grids);
    BL_ASSERT(visc_terms.nComp() >= sComp+nspecies);

    const int        nGrowOp = 1; // Required by the operator to compute first-cut fluxes
    const Real*      dx      = geom.CellSize();

    const Array<int> rho_flag(nspecies,2); // Hardwired, for now.

    MultiFab s_tmp(grids,1,nGrowOp);

    const TimeLevel whichTime = which_time(State_Type,time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
    BL_ASSERT(first_spec == Density+1);
    //
    // Create and fill a full MultiFab of all species at this level and below.
    //
    FArrayBox tmp;

    MultiFab rho_and_species(grids,nspecies+1,nGrowOp);

    for (FillPatchIterator fpi(*this,rho_and_species,nGrowOp,time,State_Type,Density,nspecies+1);
         fpi.isValid();
         ++fpi)
    {
        FArrayBox& fab = rho_and_species[fpi];

        fab.copy(fpi(),0,0,nspecies+1);

        tmp.resize(fab.box(),1);
        tmp.copy(fab,0,0,1);
        tmp.invert(1);

        for (int comp = 0; comp < nspecies; ++comp) 
            if (rho_flag[comp] == 2)
                fab.mult(tmp,0,comp+1,1);
    }

    tmp.clear();

    MultiFab rho_and_species_crse;

    if (level > 0) 
    {
        int nGrow = 1;
        HeatTransfer& coarser  = *(HeatTransfer*) &(parent->getLevel(level-1));

        rho_and_species_crse.define(coarser.grids,nspecies+1,nGrowOp,Fab_allocate);

        for (FillPatchIterator fpi(coarser,rho_and_species_crse,nGrow,time,State_Type,Density,nspecies+1);
             fpi.isValid();
             ++fpi)
        {
            FArrayBox& fab = rho_and_species_crse[fpi];

            fab.copy(fpi(),0,0,nspecies+1);

            tmp.resize(fab.box(),1);
            tmp.copy(fab,0,0,1);
            tmp.invert(1);

            for (int comp = 0; comp < nspecies; ++comp) 
                if (rho_flag[comp] == 2)
                    fab.mult(tmp,0,comp+1,1);
        }
    }

    FluxBoxes fb_flux(this, 1, 0);
    FluxBoxes fb_beta(this, nspecies, 0);
    MultiFab **flux = fb_flux.get();
    MultiFab **beta = fb_beta.get();

    getDiffusivity(beta, time, first_spec, 0, nspecies);

    for (int comp = 0; comp < nspecies; ++comp)
    {
        const int state_ind = first_spec + comp;
	ViscBndry visc_bndry;
	diffusion->getBndryDataGivenS(visc_bndry,rho_and_species,rho_and_species_crse,
                                      state_ind,comp+1,1);

	ABecLaplacian visc_op(visc_bndry,dx);
	visc_op.setScalars(0.0,1.0);
	
	for (int d = 0; d < BL_SPACEDIM; d++)
	{
            MultiFab bcoeffs(area[d].boxArray(),1,0);
	    MultiFab::Copy(bcoeffs, area[d], 0, 0, 1, 0);
	    for (MFIter mfi(bcoeffs); mfi.isValid(); ++mfi)
		bcoeffs[mfi].mult((*beta[d])[mfi],comp,0,1);
	    visc_op.bCoefficients(bcoeffs,d);
	}

	MultiFab::Copy(s_tmp,rho_and_species,comp+1,0,1,0);
	
        visc_op.compFlux(D_DECL(*flux[0],*flux[1],*flux[2]),s_tmp);
        MultiFab* const* fluxKeep = (whichTime == AmrOldTime) ? SpecDiffusionFluxn : SpecDiffusionFluxnp1;
	for (int d = 0; d < BL_SPACEDIM; ++d)
	    MultiFab::Copy(*fluxKeep[d],*flux[d],0,comp,1,0);
	spec_diffusion_flux_computed[comp] = HT_ExplicitDiffusion;
    }

    s_tmp.clear();
    rho_and_species.clear();
    rho_and_species_crse.clear();
    //
    // Modify update/fluxes to preserve flux sum = 0, compute -Div(flux)
    // (use dt = 1.0,  since the routine actually updates does "state"-=Div(flux)*dt)
    //
    const Real      dt        = 1.0;
    const MultiFab* old_state = 0;
    const MultiFab* delta_rhs = 0;
    const MultiFab* alpha     = 0;
    const int       dataComp  = 0;
    adjust_spec_diffusion_update(visc_terms,old_state,sComp,dt,time,rho_flag,
				 get_rho_half_time(),dataComp,delta_rhs,alpha,beta);
}

void
HeatTransfer::getTempViscTerms (MultiFab& visc_terms,
                                int       src_comp, 
                                Real      time)
{
    //
    // If only one species,
    // this computes 1/c_p (div lambda grad T + div q_rad)
    // if more than one species,
    //               1/c_p (div lambda grad T + div q_rad + 
    //               sum_l rho D grad Y_l dot grad h_l)
    //
    // NOTE: This routine does not fill grow cells
    //
    BL_ASSERT(visc_terms.boxArray()==grids);

    const int nGrow     = 0;
    const int rho_flag  = 1;

    FluxBoxes fb_beta(this, 1, nGrow);
    MultiFab** beta = fb_beta.get();
    //
    // + div lambda grad T + Q
    //
    getDiffusivity(beta, time, Temp, 0, 1);
    diffusion->getViscTerms(visc_terms,src_comp,Temp,time,rho_flag,beta,0);
    fb_beta.clear();
    add_heat_sources(visc_terms,Temp-src_comp,time,nGrow,1.0);
    MultiFab delta_visc_terms(grids,1,nGrow);
    //
    // + sum_l rho D grad Y_l dot grad h_l
    //
    compute_rhoDgradYgradH(time,delta_visc_terms);
    //
    // Add to visc terms, then divide whole mess by c_p 
    //
    FArrayBox cp, spec, tmp;

    MultiFab& S = get_data(State_Type,time);

    for (MFIter dvt_mfi(delta_visc_terms); dvt_mfi.isValid(); ++dvt_mfi)
    {
        const Box& box = dvt_mfi.validbox();

        cp.resize(box,1);
        spec.resize(box,nspecies);

        spec.copy(S[dvt_mfi],box,first_spec,box,0,nspecies);

        tmp.resize(box,1);
        tmp.copy(S[dvt_mfi],Density,0,1);
        tmp.invert(1);

        for (int i = 0; i < nspecies; ++i)
            spec.mult(tmp,0,i,1);

        const int yComp   = 0;
        const int sCompCp = 0;
        getChemSolve().getCpmixGivenTY(cp,S[dvt_mfi],spec,box,Temp,yComp,sCompCp);
        
        visc_terms[dvt_mfi].plus(delta_visc_terms[dvt_mfi],box,0,Temp-src_comp,1);
        visc_terms[dvt_mfi].divide(cp,0,Temp-src_comp,1);
    }
}


void
HeatTransfer::add_heat_sources(MultiFab& sum,
                               int       dComp,
                               Real      time,
                               int       nGrow,
                               Real      scale)
{
    //
    // - div q rad
    //
    if (do_OT_radiation || do_heat_sink)
    {
        BL_ASSERT(sum.boxArray() == grids);
        MultiFab dqrad(grids,1,nGrow);
        compute_OT_radloss(time,nGrow,dqrad);
        for (MFIter mfi(sum); mfi.isValid(); ++mfi)
        {
            FArrayBox& sumFab = sum[mfi];
            FArrayBox& Qfab = dqrad[mfi];
            if (scale != 1)
                Qfab.mult(scale);
            sumFab.minus(Qfab,mfi.validbox(),0,dComp,1);
        }
    }    
}

void
HeatTransfer::set_rho_to_species_sum (MultiFab& S,
                                      int       strtcomp, 
                                      int       nghost_in,
                                      int       minzero)
{
    set_rho_to_species_sum(S, strtcomp, S, strtcomp, nghost_in, minzero);
}

//
// This function
//       sets the density component in S_out to the sum of 
//       the species components in S_in
// if minzero = 1, the species components are first "max-ed" w/ zero
//  thus changing S_in values    
//
// s_in_strt is the state component corresponding to the
// 0-th component of S_in. It is otherwise assumed that 
// that the components of S_in "align" with those in S_out.
//
void
HeatTransfer::set_rho_to_species_sum (MultiFab& S_in, 
                                      int       s_in_strt,
                                      MultiFab& S_out,
                                      int       s_out_strt,  
                                      int       nghost_in,
                                      int       minzero)

{
    const BoxArray& sgrids = S_in.boxArray();

    BL_ASSERT(sgrids == S_out.boxArray());

    const int s_first_spec = first_spec - s_in_strt;
    const int s_last_spec  = last_spec  - s_in_strt;
    const int s_num_spec   = last_spec - first_spec + 1;
    const int s_density    = Density-s_out_strt;
    const int nghost       = std::min(S_in.nGrow(), std::min(S_out.nGrow(), nghost_in));

    if (minzero)
    {
        for (MFIter mfi(S_in); mfi.isValid(); ++mfi)
        {
            const int i   = mfi.index();
            Box       box = BoxLib::grow(sgrids[i], nghost);
            FabMinMax(S_in[mfi], box, 0.0, Real_MAX, s_first_spec, s_num_spec);
        }
    }

    BL_ASSERT(s_density >= 0);

    S_out.setVal(0, s_density, 1, nghost);

    for (MFIter mfi(S_in); mfi.isValid(); ++mfi)
    {
        const int i   = mfi.index();
        Box       box = BoxLib::grow(sgrids[i],nghost);

        for (int spec = s_first_spec; spec <= s_last_spec; spec++)
        {
            S_out[mfi].plus(S_in[mfi],box,spec,s_density,1);
        }
    }

    make_rho_curr_time();
}

void
HeatTransfer::scale_species (MultiFab& S,
                             int       strtcomp, 
                             int       minzero)
    //
    // This function 
    //       scales the species components in S so
    //       that they add up to the density component.
    //
    // if minzero = 1, the species components are first "max-ed" w/ zero
    //
    // strtcomp is the state component corresponding to the
    // 0-th component of S. It is otherwise assumed that 
    // that the components of S "align" with those in the state.
    //
{
    const BoxArray& sgrids = S.boxArray();
    const int s_density    = Density-strtcomp;

    BL_ASSERT(s_density >= 0);

    const int s_first_spec = first_spec - strtcomp;
    const int s_last_spec  = last_spec  - strtcomp;
    const int nghost       = S.nGrow();

    if (minzero)
    {
        //
        // max (rho*Y)_l w/ zero
        //
        for (MFIter mfi(S); mfi.isValid(); ++mfi)
        {
            const int k   = mfi.index();
            Box       grd = BoxLib::grow(sgrids[k],nghost);
            const int vol = grd.volume();

            for (int spec = s_first_spec; spec <= s_last_spec; spec++)
            {
                if (S[mfi].min(spec) < 0.0)
                {
                    Real* sdat = S[mfi].dataPtr(spec);

                    for (int i = 0; i < vol; i++) 
                        sdat[i] = std::max(0.0,sdat[i]);
                }
            }
        }
    }

    MultiFab factor(sgrids,1,nghost);
    //
    // compute sum_l (rho*Y)_l
    //
    factor.setVal(0,nghost);

    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        for (int spec = s_first_spec; spec <= s_last_spec; spec++)
        {
            factor[mfi].plus(S[mfi],spec,0,1);
        }
    }
    //
    // Compute rho/sum_l (rho*Y)_l and mulitply (rho*Y)_l by this factor
    //
    // We go through the various min/max contortions because S
    // may represent a change to the state instead of the state itself
    //
    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        const int k      = mfi.index();
        Real  min_factor = factor[mfi].min();
        Real  max_factor = factor[mfi].max();
        Real  min_rho    = S[mfi].min(s_density);
        Real  max_rho    = S[mfi].max(s_density);

        if ((min_factor>0.0||max_factor<0.0) && (min_rho>0.0||max_rho<0.0))
        {
            //
            // Anything potentially worrisome is non-zero --> easy case.
            //
            factor[mfi].invert(1.0,0,1);
            factor[mfi].mult(S[mfi],s_density,0,1);
            for (int spec = s_first_spec; spec <= s_last_spec; spec++)
                S[mfi].mult(factor[mfi],0,spec,1);
        }
        else
        {
            //
            //  We are here only if S represents change to the
            //  state, not the state itself
            //
            //  rho has a zero value anr/or sum has a zero value
            //
            //  guiding principal: rho is correct and the species are wrong
            //
            Real* sumdat = factor[mfi].dataPtr();
            Real* rhodat = S[mfi].dataPtr(s_density);
            Real rhoi, sumi, factori;

            Box grd = BoxLib::grow(sgrids[k],nghost);
            int vol = grd.volume();

            for (int i = 0; i < vol; i++)
            { 
                sumi = sumdat[i];
                rhoi = rhodat[i];

                if (sumi != 0.0 && rhoi != 0.0)
                {
                    //
                    // Case 1: both non-zero --> do what is done above.
                    //
                    factori = rhoi/sumi;
                    for (int spec = s_first_spec; spec <= s_last_spec; spec++) 
                        S[mfi].dataPtr(spec)[i] *= factori;
                }
                else if (sumi == 0.0 && rhoi == 0.0)
                {
                    //
                    // Case 2: both zero --> do nothing
                    //
                }
                else
                {
                    //
                    // Cases 3 and 4: one is zero, but not both
                    // --> shift non-zero species so that sum is zero
                    //
                    int nonzero=0;
                    for (int spec = s_first_spec; spec <= s_last_spec; spec++) 
                        if (S[mfi].dataPtr(spec)[i] != 0.0)
                            nonzero++;
                    Real shift = (rhoi-sumi)/nonzero;
                    for (int spec = s_first_spec; spec <= s_last_spec; spec++) 
                        if (S[mfi].dataPtr(spec)[i] != 0.0) 
                            S[mfi].dataPtr(spec)[i] += shift; 
                }
            }
        }
    }
}

void
HeatTransfer::temperature_stats (MultiFab& S)
{
    if (verbose)
    {
        //
        // Calculate some minimums and maximums.
        //
        Real tdhmin[3] = { 1.0e30, 1.0e30, 1.0e30};
        Real tdhmax[3] = {-1.0e30,-1.0e30,-1.0e30};

        for (MFIter S_mfi(S); S_mfi.isValid(); ++S_mfi)
        {
            const Box& bx = S_mfi.validbox();

            tdhmin[0] = std::min(tdhmin[0],S[S_mfi].min(bx,Temp));
            tdhmin[1] = std::min(tdhmin[1],S[S_mfi].min(bx,Density));
            tdhmin[2] = std::min(tdhmin[2],S[S_mfi].min(bx,RhoH));

            tdhmax[0] = std::max(tdhmax[0],S[S_mfi].max(bx,Temp));
            tdhmax[1] = std::max(tdhmax[1],S[S_mfi].max(bx,Density));
            tdhmax[2] = std::max(tdhmax[2],S[S_mfi].max(bx,RhoH));
        }

        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceRealMin(tdhmin,3,IOProc);
        ParallelDescriptor::ReduceRealMax(tdhmax,3,IOProc);

        FArrayBox   Y, tmp;
        bool        aNegY = false;
        Array<Real> minY(nspecies,1.e30);

        for (MFIter S_mfi(S); S_mfi.isValid(); ++S_mfi)
        {
            Y.resize(S_mfi.validbox(),1);

            tmp.resize(S_mfi.validbox(),1);
            tmp.copy(S[S_mfi],Density,0,1);
            tmp.invert(1);

            for (int i = 0; i < nspecies; ++i)
            {
                Y.copy(S[S_mfi],first_spec+i,0,1);
                Y.mult(tmp,0,0,1);
                minY[i] = std::min(minY[i],Y.min(0));
            }
        }

        ParallelDescriptor::ReduceRealMin(minY.dataPtr(),nspecies,IOProc);

        for (int i = 0; i < nspecies; ++i)
        {
            if (minY[i] < 0) aNegY = true;
        }

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "  Min,max temp = " << tdhmin[0] << ", " << tdhmax[0] << '\n';
            std::cout << "  Min,max rho  = " << tdhmin[1] << ", " << tdhmax[1] << '\n';
            std::cout << "  Min,max rhoh = " << tdhmin[2] << ", " << tdhmax[2] << '\n';
            if (aNegY)
            {
                const Array<std::string>& names = HeatTransfer::getChemSolve().speciesNames();
                std::cout << "  Species w/min < 0: ";
                for (int i = 0; i < nspecies; ++i)
                    if (minY[i] < 0)
                        std::cout << "Y(" << names[i] << ") [" << minY[i] << "]  ";
                std::cout << '\n';
            }
        }
    }
}

void
HeatTransfer::compute_rhoRT (const MultiFab& S,
                             MultiFab&       p, 
                             int             pComp,
                             const MultiFab* temp)
{
    BL_ASSERT(pComp<p.nComp());

    const BoxArray& sgrids = S.boxArray();
    int nCompY = last_spec - first_spec + 1;

    FArrayBox state, tmp;

    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        const int  i      = mfi.index();
	const Box& box    = sgrids[i];
	const int  sCompR = 0;
	const int  sCompT = 1;
	const int  sCompY = 2;
	
        state.resize(box,nCompY+2);
	BL_ASSERT(S[mfi].box().contains(box));
	state.copy(S[mfi],box,Density,box,sCompR,1);

	if (temp)
	{
	    BL_ASSERT(temp->boxArray()[i].contains(box));
	    state.copy((*temp)[mfi],box,0,box,sCompT,1);
	}
        else
        {
	    state.copy(S[mfi],box,Temp,box,sCompT,1);
	}
	state.copy(S[mfi],box,first_spec,box,sCompY,nCompY);

        tmp.resize(box,1);
        tmp.copy(state,sCompR,0,1);
        tmp.invert(1);

        for (int k = 0; k < nCompY; k++)
            state.mult(tmp,0,sCompY+k,1);

	getChemSolve().getPGivenRTY(p[mfi],state,state,state,box,sCompR,sCompT,sCompY,pComp);
    }
}
			   
//
// Setup for the advance function.
//

#ifndef NDEBUG
#if defined(BL_OSF1)
#if defined(BL_USE_DOUBLE)
const Real BL_BOGUS      = DBL_QNAN;
#else
const Real BL_BOGUS      = FLT_QNAN;
#endif
#else
const Real BL_BOGUS      = 1.e30;
#endif
#endif

void
HeatTransfer::set_htt_hmixTYP ()
{
    const int finest_level = parent->finestLevel();

    // set typical value for hmix, needed for TfromHY solves if not provided explicitly
    if (typical_values[RhoH]==typical_RhoH_value_default)
    {
        htt_hmixTYP = 0;
        std::vector< std::pair<int,Box> > isects;
        for (int k = 0; k <= finest_level; k++)
        {
            AmrLevel& ht = getLevel(k);
            const MultiFab& S = ht.get_new_data(State_Type);
            const BoxArray& ba = ht.boxArray();
            MultiFab hmix(ba,1,0,Fab_allocate);
            MultiFab::Copy(hmix,S,RhoH,0,1,0);
            MultiFab::Divide(hmix,S,Density,0,1,0);
            if (k!=finest_level) 
            {
                AmrLevel& htf = getLevel(k+1);
                BoxArray baf = BoxArray(htf.boxArray()).coarsen(parent->refRatio(k));
                for (MFIter mfi(hmix); mfi.isValid(); ++mfi)
                {
                    baf.intersections(ba[mfi.index()],isects);
                    for (int i = 0; i < isects.size(); i++)
                    {
                        hmix[mfi].setVal(0,isects[i].second,0,1);
                    }
                }
            }
            htt_hmixTYP = std::max(htt_hmixTYP,hmix.norm0(0));
        }
        ParallelDescriptor::ReduceRealMax(htt_hmixTYP);
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "setting htt_hmixTYP(via domain scan) = " << htt_hmixTYP << '\n';
    }
    else
    {
        htt_hmixTYP = typical_values[RhoH];
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "setting htt_hmixTYP(from user input) = " << htt_hmixTYP << '\n';
    }
}

void
HeatTransfer::advance_setup (Real time,
                             Real dt,
                             int  iteration,
                             int  ncycle)
{
    NavierStokesBase::advance_setup(time, dt, iteration, ncycle);
    //
    // Make sure the new state has values so that c-n works
    // in the predictor (of the predictor)--rbp.
    //
    for (int k = 0; k < num_state_type; k++)
    {
        MultiFab& nstate = get_new_data(k);
        MultiFab& ostate = get_old_data(k);

        MultiFab::Copy(nstate,ostate,0,0,nstate.nComp(),nstate.nGrow());
    }
    if (level == 0)
        set_htt_hmixTYP();

    make_rho_curr_time();

#ifndef NDEBUG
    aux_boundary_data_old.setVal(BL_BOGUS);
    aux_boundary_data_new.setVal(BL_BOGUS);
#endif
    //
    // Set a dumbbell flag to help avoid stupid mistakes
    //
    for (int i = 0; i < spec_diffusion_flux_computed.size(); ++i)
	spec_diffusion_flux_computed[i] = HT_None;

    if (plot_reactions && level == 0)
    {
        for (int i = parent->finestLevel(); i >= 0; --i)
            getLevel(i).auxDiag["REACTIONS"]->setVal(0);
    }
}

void
HeatTransfer::setThermoPress(Real time)
{
    const TimeLevel whichTime = which_time(State_Type,time);
    
    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
    
    MultiFab& S = (whichTime == AmrOldTime) ? get_old_data(State_Type) : get_new_data(State_Type);
    
    const int pComp = (have_rhort ? RhoRT : Trac);

    compute_rhoRT (S,S,pComp);
}

#if 1
Real
HeatTransfer::predict_velocity (Real  dt,
                                Real& comp_cfl)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... predict edge velocities\n";
    //
    // Get simulation parameters.
    //
    const int   nComp          = BL_SPACEDIM;
    const Real* dx             = geom.CellSize();
    const Real  prev_time      = state[State_Type].prevTime();
    const Real  prev_pres_time = state[Press_Type].prevTime();
    //
    // Compute viscous terms at level n.
    // Ensure reasonable values in 1 grow cell.  Here, do extrap for
    // c-f/phys boundary, since we have no interpolator fn, also,
    // preserve extrap for corners at periodic/non-periodic intersections.
    //
    MultiFab visc_terms(grids,nComp,1);

    if (be_cn_theta != 1.0)
    {
	getViscTerms(visc_terms,Xvel,nComp,prev_time);
    }
    else
    {
	visc_terms.setVal(0);
    }
    //
    // Set up the timestep estimation.
    //
    Real cflgrid,u_max[3];
    Real cflmax = 1.0e-10;
    comp_cfl    = (level == 0) ? cflmax : comp_cfl;

    FArrayBox tforces;

    Array<int> bndry[BL_SPACEDIM];

    MultiFab Gp(grids,BL_SPACEDIM,1);

    getGradP(Gp, prev_pres_time);
    
    FArrayBox null_fab;

    showMF("mac",Gp,"pv_Gp",level);
    showMF("mac",get_old_data(State_Type),"pv_S_old",level);
    showMF("mac",rho_ptime,"pv_rho_old",level);
    showMF("mac",visc_terms,"pv_visc_terms",level);

#ifndef NDEBUG
    for (int dir=0; dir<BL_SPACEDIM; ++dir) {
        u_mac[dir].setVal(0.);
    }
    MultiFab Force(grids,BL_SPACEDIM,Godunov::hypgrow());
    Force.setVal(0);
#endif

    for (FillPatchIterator U_fpi(*this,visc_terms,Godunov::hypgrow(),prev_time,State_Type,Xvel,BL_SPACEDIM)
#ifdef MOREGENGETFORCE
	     , S_fpi(*this,visc_terms,1,prev_time,State_Type,Density,NUM_SCALARS);
	 S_fpi.isValid() && U_fpi.isValid();
	 ++S_fpi, ++U_fpi
#else
             ; U_fpi.isValid();
	 ++U_fpi
#endif
	)
    {
        const int i = U_fpi.index();

#ifdef GENGETFORCE
        getForce(tforces,i,1,Xvel,BL_SPACEDIM,prev_time,rho_ptime[U_fpi]);
#elif MOREGENGETFORCE
	if (ParallelDescriptor::IOProcessor() && getForceVerbose)
	    std::cout << "---" << std::endl << "A - Predict velocity:" << std::endl << " Calling getForce..." << std::endl;
        getForce(tforces,i,1,Xvel,BL_SPACEDIM,prev_time,U_fpi(),S_fpi(),0);
#else
	getForce(tforces,i,1,Xvel,BL_SPACEDIM,rho_ptime[U_fpi]);
#endif		 

        //
        // Test velocities, rho and cfl.
        //
        cflgrid  = godunov->test_u_rho(U_fpi(),rho_ptime[U_fpi],grids[i],dx,dt,u_max);
        cflmax   = std::max(cflgrid,cflmax);
        comp_cfl = std::max(cflgrid,comp_cfl);
        //
        // Compute the total forcing.
        //
        godunov->Sum_tf_gp_visc(tforces,visc_terms[U_fpi],Gp[U_fpi],rho_ptime[U_fpi]);

#ifndef NDEBUG
        Force[U_fpi].copy(tforces,0,0,BL_SPACEDIM);
#endif

        D_TERM(bndry[0] = getBCArray(State_Type,i,0,1);,
               bndry[1] = getBCArray(State_Type,i,1,1);,
               bndry[2] = getBCArray(State_Type,i,2,1);)

        godunov->Setup(grids[i], dx, dt, 1,
                       null_fab, bndry[0].dataPtr(),
                       null_fab, bndry[1].dataPtr(),
#if (BL_SPACEDIM == 3)                         
                       null_fab, bndry[2].dataPtr(),
#endif
                       U_fpi(), rho_ptime[U_fpi], tforces);

        godunov->ComputeUmac(grids[i], dx, dt, 
                             u_mac[0][U_fpi], bndry[0].dataPtr(),
                             u_mac[1][U_fpi], bndry[1].dataPtr(),
#if (BL_SPACEDIM == 3)
                             u_mac[2][U_fpi], bndry[2].dataPtr(),
#endif
                             U_fpi(), tforces);
    }

    showMF("mac",u_mac[0],"pv_umac0",level);
    showMF("mac",u_mac[1],"pv_umac1",level);
#if BL_SPACEDIM==3
    showMF("mac",u_mac[2],"pv_umac2",level);
#endif
#ifndef NDEBUG
    showMF("mac",Force,"pv_force",level);
#endif

    Real tempdt = std::min(change_max,cfl/cflmax);

    ParallelDescriptor::ReduceRealMin(tempdt);

    return dt*tempdt;
}
#endif

Real
HeatTransfer::advance (Real time,
                       Real dt,
                       int  iteration,
                       int  ncycle)
{
    BL_PROFILE("HeatTransfer::advance()");

    if (level == 0)
    {
        crse_dt = dt;
        int thisLevelStep = parent->levelSteps(0);
        FORT_SET_COMMON(&time,&thisLevelStep);
    }

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "Advancing level "    << level
                  << " : starting time = " << time
                  << " with dt = "         << dt << '\n';
    }

    advance_setup(time,dt,iteration,ncycle);

    if (level==0 && reset_typical_vals_int>0)
    {
      int L0_steps = parent->levelSteps(0);
      if (L0_steps>0 && L0_steps%reset_typical_vals_int==0)
      {
        reset_typical_values(get_old_data(State_Type));
      }
    }

    if (do_check_divudt)
        checkTimeStep(dt);
    
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& S_old = get_old_data(State_Type);

    //
    // Reset flags that fill patched state data is good
    //
    FillPatchedOldState_ok = true;
    FillPatchedNewState_ok = true;
    //
    // Compute traced states for normal comp of velocity at half time level.
    //
    Real dt_test = 0.0, dummy = 0.0;    
    dt_test = predict_velocity(dt,dummy);
    
    showMF("mac",u_mac[0],"adv_umac0",level);
    showMF("mac",u_mac[1],"adv_umac1",level);
#if BL_SPACEDIM==3
    showMF("mac",u_mac[2],"adv_umac2",level);
#endif
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "HeatTransfer::advance(): at start of time step\n";

    temperature_stats(S_old);

    const Real prev_time = state[State_Type].prevTime();
    const Real cur_time  = state[State_Type].curTime();
#if 0
    setThermoPress(prev_time);
#endif
    //
    // Do MAC projection and update edge velocities.
    //
    if (do_mac_proj) 
    {
        const int havedivu = 1;
        MultiFab mac_rhs(grids,1,0);
        create_mac_rhs(mac_rhs,0,time,dt);
        showMF("mac",mac_rhs,"mac_rhs",level);
        mac_project(time,dt,S_old,&mac_rhs,havedivu,umac_n_grow,true);
    }

    if (do_mom_diff == 0)
        velocity_advection(dt);
    //
    // Set old-time boundary data for RhoH
    // 
    set_overdetermined_boundary_cells(time);
    //
    // Advance the old state for a Strang-split dt/2.  Include grow cells in
    // advance, and squirrel these away for diffusion and Godunov guys to
    // access for overwriting non-advanced fill-patched grow cell data.
    //
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... advancing chem\n";

    BL_ASSERT(S_new.boxArray() == S_old.boxArray());

    const int nComp = NUM_STATE - BL_SPACEDIM;
    Array<int> consumptionComps; // Put in scope for work below
    if (plot_consumption)
    {
        //
        // Save off a copy of the pre-chem state
        //
        consumptionComps.resize(consumptionName.size());
        for (int j=0; j<consumptionComps.size(); ++j)
        {
            consumptionComps[j] = getChemSolve().index(consumptionName[j]) + first_spec;
            auxDiag["CONSUMPTION"]->copy(S_old,consumptionComps[j],j,1);
        }
    }

    MultiFab Qtmp; // Put in scope for work below
    if (plot_heat_release)
    {
        //
        // Save off a copy of the pre-chem state
        //
        Qtmp.define(grids,getChemSolve().numSpecies(),0,Fab_allocate);
        MultiFab::Copy(Qtmp,S_old,first_spec,0,Qtmp.nComp(),0);
    }
    //
    // Build a MultiFab parallel to "fabs".  Force it to have the
    // same distribution as aux_boundary_data_old.  This'll cut out a
    // parallel copy.  It doesn't happen by default (like you might think)
    // due to not being recached appropriately after regridding.
    //
    {
        MultiFab tmpFABs;

        tmpFABs.define(aux_boundary_data_old.equivBoxArray(),
                       NUM_STATE,
                       0,
                       aux_boundary_data_old.DistributionMap(),
                       Fab_allocate);

        const int ngrow = aux_boundary_data_old.nGrow();

        {
            MultiFab tmpS_old(S_old.boxArray(), NUM_STATE, ngrow);

	    FillPatch(*this,tmpS_old,ngrow,prev_time,State_Type,0,NUM_STATE,0);

            tmpFABs.copy(tmpS_old,0,0,NUM_STATE,ngrow,0);
        }

        strang_chem(S_old,  dt,HT_LeaveYdotAlone);
        strang_chem(tmpFABs,dt,HT_LeaveYdotAlone,ngrow);

        aux_boundary_data_old.copyFrom(tmpFABs,BL_SPACEDIM,0,nComp);
    }
    //
    // Find change due to first Strang step.
    //
    if (plot_consumption)
    {
        MultiFab tmp(auxDiag["CONSUMPTION"]->boxArray(),consumptionComps.size(),0);
        tmp.setVal(0);
        for (int j=0; j<consumptionComps.size(); ++j)
        {
            tmp.copy(S_old,consumptionComps[j],j,1);
        }
        for (MFIter mfi(*auxDiag["CONSUMPTION"]); mfi.isValid(); ++mfi)
        {
            (*auxDiag["CONSUMPTION"])[mfi].minus(tmp[mfi],0,0,consumptionComps.size());
        }
    }
    if (plot_heat_release)
    {
        for (MFIter mfi(Qtmp); mfi.isValid(); ++mfi)
        {
            Qtmp[mfi].minus(S_old[mfi],first_spec,0,Qtmp.nComp());
        }
    }

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "HeatTransfer::advance(): after first Strang step\n";

    temperature_stats(S_old);
    //
    // Activate hook in FillPatch hack to get "better" OLD data now.
    //
    FillPatchedOldState_ok = false;
    //
    // Compute tn coeffs based on chem-advance tn data
    //  (these are used in the Godunov extrapolation)
    //
    const int nScalDiffs = NUM_STATE-BL_SPACEDIM-1;
    calcDiffusivity(prev_time);
    //
    // Godunov-extrapolate states to cell edges
    // computes actual states, NOT fluxes, into EdgeStates
    //
    compute_edge_states(dt);
    //
    // Compute advective fluxes divergences, where possible
    // NOTE: If Le!=1 cannot do RhoH advection until after spec update fills
    //       spec diffusion fluxes.  Since it is never a bad idea to do this
    //       later (i.e. S_old does not change), we do it later always.
    //
    const int first_scalar = Density;
    const int last_scalar = first_scalar + NUM_SCALARS - 1;
    bool do_adv_reflux = do_reflux;
    //
    //  Load the advective flux registers into aofs.
    //
    if (RhoH > first_scalar)
	scalar_advection(dt,first_scalar,RhoH-1,do_adv_reflux);  // density and species
    if (RhoH < last_scalar)
        scalar_advection(dt,RhoH+1,last_scalar,do_adv_reflux);  // temperature and tracers
    //
    // Copy old-time boundary RhoH into estimate for new-time RhoH.
    //
    aux_boundary_data_new.copy(aux_boundary_data_old,RhoH-BL_SPACEDIM,0,1);
    //
    // Save rho used in rho-states, needed for replacing with new one
    //  NOTE: WE LOAD/USE GROW CELLS HERE SO WE CAN FIX BOUNDARY DATA AS WELL
    //
    MultiFab Rho_hold(grids,1,LinOp_grow);

    BL_ASSERT(LinOp_grow == 1);

    for (MFIter mfi(rho_ctime); mfi.isValid(); ++mfi)
    {
        const Box& box = BoxLib::grow(mfi.validbox(),LinOp_grow);
        Rho_hold[mfi].copy(rho_ctime[mfi],box,0,box,0,1);
    }
    //
    // Compute new and half-time densities.
    //
    const int rho_corr = 1;
    scalar_update(dt,Density,Density,rho_corr);
    //
    // Set saved rho at current time.
    //
    make_rho_curr_time();
    //
    // Reset rho-states to contain new rho.
    //
    reset_rho_in_rho_states(Rho_hold,cur_time,first_scalar+1,NUM_SCALARS-1);

    Rho_hold.clear();
    //
    // Compute the update to momentum
    //
    if (do_mom_diff == 1)
	momentum_advection(dt,do_reflux);

    //
    // Update energy and species.
    //
    {
        //
        // Predictor
        //
        int corrector = 0;
        //
        // Set tnp1 coeffs to tn values in first round of predictor
        //
        MultiFab::Copy(*diffnp1_cc,*diffn_cc,0,0,nScalDiffs,diffn_cc->nGrow());

	// Eq (13) in DayBell
	// Here, predict n+1 coeffs using n coeffs
	// results in Ttilde^np1*
        temp_update(dt,corrector);         
        temperature_stats(S_new);

	// compute coefficients before Eq (14a) in DayBell
	// results in D_m^{n+1,star}
        calcDiffusivity(cur_time);

	// Eq (14) and (14a) in DayBell, results in corrected Gamma^{n+1,star}
        spec_update(time,dt,corrector);
        
        set_overdetermined_boundary_cells(time + dt); // RhoH BC's to see new Y's at n+1
        //
        // Activate hook in FillPatch hack to get "better" NEW data now.
        //
        FillPatchedNewState_ok = false;
        
        do_adv_reflux = false;

	// we appear to be missing the call to calc_diffusivity that appears in
	// DayBell after equation (14a)

        scalar_advection(dt,RhoH,RhoH,do_adv_reflux); // Get aofs for RhoH now

        rhoh_update(time,dt,corrector);

        RhoH_to_Temp(S_new);

        temperature_stats(S_new);
        //
        // Corrector
        //
        corrector = 1;
        calcDiffusivity(cur_time);
        tracer_update(dt,corrector);
        spec_update(time,dt,corrector);
        
        set_overdetermined_boundary_cells(time + dt);// RhoH BC's to see new Y's at n+1
        
        do_adv_reflux = do_reflux;
        scalar_advection(dt,RhoH,RhoH,do_adv_reflux); // Get aofs for RhoH now
        
        rhoh_update(time,dt,corrector);
        RhoH_to_Temp(S_new); 
    }

    //
    // Second half of Strang-split chemistry (first half done in
    // compute_edge_states) This takes new-time data, and returns new-time
    // data, as well as providing a predicted Ydot for the velocity
    // constraint.  We write the result over the new state, but only care
    // about stuff in the valid region.
    //
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "HeatTransfer::advance(): after scalar_update\n";

    temperature_stats(S_new);

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... advancing chem\n";
    //
    // Adjust chemistry diagnostic before and after reactions.
    //
    if (plot_consumption)
    {
        for (MFIter mfi(*auxDiag["CONSUMPTION"]); mfi.isValid(); ++mfi)
        {
            for (int j=0; j<consumptionComps.size(); ++j)
            {
                (*auxDiag["CONSUMPTION"])[mfi].plus(S_new[mfi],consumptionComps[j],j,1);
            }
        }
    }
    if (plot_heat_release)
    {
        for (MFIter mfi(Qtmp); mfi.isValid(); ++mfi)
        {
            Qtmp[mfi].plus(S_new[mfi],first_spec,0,Qtmp.nComp());
        }
    }

    strang_chem(S_new,dt,HT_EstimateYdotNew);

    if (plot_consumption)
    {
        for (MFIter mfi(*auxDiag["CONSUMPTION"]); mfi.isValid(); ++mfi)
        {
            for (int j=0; j<consumptionComps.size(); ++j)
            {
                (*auxDiag["CONSUMPTION"])[mfi].minus(S_new[mfi],consumptionComps[j],j,1);
            }
            (*auxDiag["CONSUMPTION"])[mfi].mult(1.0/dt);
        }
    }
    if (plot_heat_release)
    {
        FArrayBox enthi, T;
        for (MFIter mfi(Qtmp); mfi.isValid(); ++mfi)
        {
            Qtmp[mfi].minus(S_new[mfi],first_spec,0,Qtmp.nComp());
            Qtmp[mfi].mult(1.0/dt);

            const Box& box = mfi.validbox();
            T.resize(mfi.validbox(),1);
            T.setVal(298.15);

            enthi.resize(mfi.validbox(),Qtmp.nComp());
            getChemSolve().getHGivenT(enthi,T,box,0,0);

            // Form heat release
            (*auxDiag["HEATRELEASE"])[mfi].setVal(0.);
            for (int j=0; j<Qtmp.nComp(); ++j)
            {
                Qtmp[mfi].mult(enthi,j,j,1);
                (*auxDiag["HEATRELEASE"])[mfi].plus(Qtmp[mfi],j,0,1);
            }
        }
    }

#ifdef BL_COMM_PROFILING
    for (MFIter mfi(*auxDiag["COMMPROF"]); mfi.isValid(); ++mfi)
    {
      int rank(ParallelDescriptor::MyProc());
      (*auxDiag["COMMPROF"])[mfi].setVal(rank, 0);
      (*auxDiag["COMMPROF"])[mfi].setVal(DistributionMapping::ProximityMap(rank),   1);
      (*auxDiag["COMMPROF"])[mfi].setVal(DistributionMapping::ProximityOrder(rank), 2);
    }
#endif

    //
    // Deactivate hook in FillPatch so that old data really is old data again.
    //
    FillPatchedOldState_ok = true;

    if (verbose && ParallelDescriptor::IOProcessor())
	std::cout << "HeatTransfer::advance(): after second Strang-split step\n";
    
    temperature_stats(S_new);
    //
    // S appears in rhs of the velocity update, so we better do it now.
    // (be sure to use most recent version of state to get
    // viscosity/diffusivity).
    //
    calcDiffusivity(cur_time,true);
    //
    // Set the dependent value of RhoRT to be the thermodynamic pressure.  By keeping this in
    // the state, we can use the average down stuff to be sure that RhoRT_avg is avg(RhoRT),
    // not ave(Rho)avg(R)avg(T), which seems to give the p-relax stuff in the mac Rhs troubles.
    //
    setThermoPress(cur_time);

    calc_divu(time+dt, dt, get_new_data(Divu_Type));

    if (!NavierStokesBase::initial_step && level != parent->finestLevel())
    {
        //
        // Set new divu to old div + dt*dsdt_old where covered by fine.
        //
        BoxArray crsndgrids = getLevel(level+1).grids;

        crsndgrids.coarsen(fine_ratio);
            
        MultiFab& divu_new = get_new_data(Divu_Type);
        MultiFab& divu_old = get_old_data(Divu_Type);
        MultiFab& dsdt_old = get_old_data(Dsdt_Type);

        std::vector< std::pair<int,Box> > isects;

        for (MFIter mfi(divu_new); mfi.isValid();++mfi)
        {
            crsndgrids.intersections(mfi.validbox(),isects);

            for (int i = 0, N = isects.size(); i < N; i++)
            {
                const Box& ovlp = isects[i].second;

                divu_new[mfi].copy(dsdt_old[mfi],ovlp,0,ovlp,0,1);
                divu_new[mfi].mult(dt,ovlp,0,1);
                divu_new[mfi].plus(divu_old[mfi],ovlp,0,0,1);
            }
        }
    }
        
    calc_dsdt(time, dt, get_new_data(Dsdt_Type));

    if (NavierStokesBase::initial_step)
        MultiFab::Copy(get_old_data(Dsdt_Type),get_new_data(Dsdt_Type),0,0,1,0);
    //
    // Add the advective and other terms to get velocity (or momentum) at t^{n+1}.
    //
    velocity_update(dt);
    //
    // Increment rho average.
    //
    if (!initial_step)
    {
        if (level > 0)
        {
            Real alpha = 1.0/Real(ncycle);
            if (iteration == ncycle)
                alpha = 0.5/Real(ncycle);
            incrRhoAvg(alpha);
        }
        //
        // Do a level project to update the pressure and velocity fields.
        //
        level_projector(dt,time,iteration);

        if (level > 0 && iteration == 1) p_avg.setVal(0);
    }

#ifdef PARTICLES
    if (theNSPC() != 0)
    {
        theNSPC()->AdvectWithUmac(u_mac, level, dt);
    }
#endif
    //
    // Deactivate hook in FillPatch so that new data really is new data again.
    //
    FillPatchedNewState_ok = true;

    advance_cleanup(iteration,ncycle);
    //
    // Update estimate for allowable time step.
    //
    dt_test = std::min(dt_test, estTimeStep());

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "HeatTransfer::advance(): at end of time step\n";

    temperature_stats(S_new);

    return dt_test;
}

void
HeatTransfer::create_mac_rhs (MultiFab& mac_rhs, int nGrow, Real time, Real dt)
{
    NavierStokesBase::create_mac_rhs(mac_rhs,nGrow,time,dt);

    showMF("mac",mac_rhs,"mac_rhs0_",level);

    if (dt > 0.0) 
    {
        MultiFab  dpdt(grids,1,0);
        calc_dpdt(time,dt,dpdt,u_mac);

        for (MFIter mfi(dpdt); mfi.isValid(); ++mfi)
            mac_rhs[mfi].plus(dpdt[mfi], grids[mfi.index()], 0,0,1);

        if (nGrow > 0)
            mac_rhs.FillBoundary();
    }

    showMF("mac",mac_rhs,"mac_rhs1_",level);
}

void
HeatTransfer::reset_rho_in_rho_states (const MultiFab& rho,
                                       Real            time,
                                       const int       sComp,
                                       const int       nComp)
{
    //
    // Divide the given rho from the states with diffusion terms of the
    // form Laplacian_SoverRho and multiply by the new Rho.
    //
    BL_ASSERT(rho.boxArray() == grids);
    //
    // Only do the valid region.
    //
    MultiFab& S = get_data(State_Type,time);

    BL_ASSERT(sComp + nComp <= S.nComp());

    FArrayBox tmp;

    for (MFIter Smfi(S); Smfi.isValid(); ++Smfi)
    {
        tmp.resize(Smfi.validbox(),1);
        tmp.copy(rho[Smfi],0,0,1);
        tmp.invert(1);

        for (int comp = sComp; comp < sComp + nComp; ++comp)
        {
            if (is_diffusive[comp] && diffusionType[comp]==Laplacian_SoverRho)
            {
                S[Smfi].mult(tmp,     Smfi.validbox(), 0,       comp, 1);
                S[Smfi].mult(S[Smfi], Smfi.validbox(), Density, comp, 1);
            }
        }
    }
}

void
HeatTransfer::set_overdetermined_boundary_cells (Real time)
{
    BL_ASSERT(first_spec == Density+1);

    const TimeLevel whichTime = which_time(State_Type,time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    AuxBoundaryData& rhoh_data = (whichTime == AmrOldTime) ? aux_boundary_data_old : aux_boundary_data_new;

    const int nGrow = (whichTime == AmrOldTime) ? Godunov::hypgrow() : LinOp_grow;
    //
    // Build a MultiFab parallel to State with appropriate # of ghost
    //
    MultiFab tmpS(grids,1,nGrow);

    const BoxArray& rhoh_BA = rhoh_data.equivBoxArray();

    const int sCompT = 0, sCompY = 1, sCompH = 0;
    //
    // A non-allocated MultiFab on grids for FPI below.
    //
    FArrayBox tmp, rhoh;

    MultiFab t(grids,1,0,Fab_noallocate);

    std::vector< std::pair<int,Box> > isects;

    for (FillPatchIterator T_fpi(*this,t,nGrow,time,State_Type,Temp,1),
             RhoY_fpi(*this,t,nGrow,time,State_Type,Density,nspecies+1);
         T_fpi.isValid() && RhoY_fpi.isValid();
         ++T_fpi, ++RhoY_fpi)
    {
        FArrayBox& RhoY = RhoY_fpi();

        BL_ASSERT(RhoY.box() == tmpS[RhoY_fpi].box());

        tmp.resize(RhoY.box(),1);
        tmp.copy(RhoY,0,0,1);
        tmp.invert(1);

        rhoh_BA.intersections(RhoY.box(),isects);

        for (int i = 0, N = isects.size(); i < N; i++)
        {
            const Box& isect = isects[i].second;

            for (int j = 1; j < nspecies+1; j++)
                RhoY.mult(tmp,isect,0,j,1);

            rhoh.resize(isect,1);
            getChemSolve().getHmixGivenTY(rhoh,T_fpi(),RhoY,isect,sCompT,sCompY,sCompH);
            rhoh.mult(RhoY,0,0,1);

            tmpS[T_fpi].copy(rhoh);
        }
    }

    const int RhoHcomp = (whichTime == AmrOldTime) ? RhoH-BL_SPACEDIM : 0;

    rhoh_data.copyFrom(tmpS,0,RhoHcomp,1,nGrow); // Parallel copy.
}

DistributionMapping
HeatTransfer::getFuncCountDM (const BoxArray& bxba, int ngrow)
{
    //
    // Sometimes "mf" is the valid region of the State.
    // Sometimes it's the region covered by AuxBoundaryData.
    // When ngrow>0 were doing AuxBoundaryData with nGrow()==ngrow.
    //
    DistributionMapping rr;
    rr.RoundRobinProcessorMap(bxba.size(),ParallelDescriptor::NProcs());

    MultiFab fctmpnew;
    fctmpnew.define(bxba, 1, 0, rr, Fab_allocate);
    fctmpnew.setVal(1);

    const MultiFab& FC = get_new_data(FuncCount_Type);
    fctmpnew.copy(FC,0,0,1,std::min(ngrow,FC.nGrow()),0);

    int count = 0;
    Array<long> vwrk(bxba.size());
    for (MFIter mfi(fctmpnew); mfi.isValid(); ++mfi)
        vwrk[count++] = static_cast<long>(fctmpnew[mfi].sum(0));

    fctmpnew.clear();

#if BL_USE_MPI
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    Array<int> nmtags(ParallelDescriptor::NProcs(),0);
    Array<int> offset(ParallelDescriptor::NProcs(),0);

    for (int i = 0; i < vwrk.size(); i++)
        nmtags[rr.ProcessorMap()[i]]++;

    BL_ASSERT(nmtags[ParallelDescriptor::MyProc()] == count);

    for (int i = 1; i < offset.size(); i++)
        offset[i] = offset[i-1] + nmtags[i-1];

    Array<long> vwrktmp = vwrk;

    BL_COMM_PROFILE(BLProfiler::Gatherv, vwrktmp.size() * sizeof(long),
                    ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());

    MPI_Gatherv(vwrk.dataPtr(),
                count,
                ParallelDescriptor::Mpi_typemap<long>::type(),
                vwrktmp.dataPtr(),
                nmtags.dataPtr(),
                offset.dataPtr(),
                ParallelDescriptor::Mpi_typemap<long>::type(),
                IOProc,
                ParallelDescriptor::Communicator());

    BL_COMM_PROFILE(BLProfiler::Gatherv, vwrktmp.size() * sizeof(long),
                    ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // We must now assemble vwrk in the proper order.
        //
        std::vector< std::vector<int> > table(ParallelDescriptor::NProcs());

        for (int i = 0; i < vwrk.size(); i++)
            table[rr.ProcessorMap()[i]].push_back(i);

        int idx = 0;
        for (int i = 0; i < table.size(); i++)
            for (int j = 0; j < table[i].size(); j++)
                vwrk[table[i][j]] = vwrktmp[idx++]; 
    }
    //
    // Send the properly-ordered vwrk to all processors.
    //
    ParallelDescriptor::Bcast(vwrk.dataPtr(), vwrk.size(), IOProc);
#endif

    DistributionMapping res;
    //
    // This call doesn't invoke the MinimizeCommCosts() stuff.
    //
    res.KnapSackProcessorMap(vwrk,ParallelDescriptor::NProcs());

    return res;
}

void
HeatTransfer::strang_chem (MultiFab&  mf,
                           Real       dt,
                           YdotAction Ydot_action,
                           int        ngrow)
{
    BL_PROFILE("HeatTransfer::strang_chem()");
    //
    // Sometimes "mf" is the valid region of the State.
    // Sometimes it's the region covered by AuxBoundaryData.
    // When ngrow>0 we're doing AuxBoundaryData with nGrow()==ngrow.
    //
    if (verbose && benchmarking) ParallelDescriptor::Barrier();

    const Real strt_time = ParallelDescriptor::second();
    //
    // I intend that this function be called just prior to the Godunov
    // extrapolation step for species and temperature (i.e. with FillPatched
    // data in the valid and grow cells), or just after other processes to
    // finish the chem evolve  Here, we:
    //
    //  (1) Carry out the Strang-split chemistry advance (for half time step).
    //
    //  (2) [potentially] Estimate, or improve the value of ydot.
    //       For improving ydot, it is assumed that ydot presently holds the
    //       effective changes in mass fraction, in terms of a rate computed
    //       over the second half of the timestep prior.  We can center that
    //       estimate by averaging with the effective rate over this first half
    // Note:
    //   The dt passed in is the full time step for this level ... and
    //   mf is in State_Type ordering, but starts at the scalars.
    //
    const int rho_comp  = Density; // mf and State_Type completely aligned here
    const int dCompYdot = 0;       // first component of Ydot corres. to first_spec
    const int ycomp     = first_spec - Density + rho_comp;
    const int Tcomp     = Temp - Density + rho_comp;

    MultiFab junk, *ydot_tmp = 0;

    if (Ydot_action == HT_ImproveYdotOld)
    {
        MultiFab& ydot_old = get_old_data(Ydot_Type);
  	junk.define(ydot_old.boxArray(),nspecies,0,Fab_allocate);
	ydot_tmp = &junk;
    }
    else if (Ydot_action == HT_EstimateYdotNew)
    {
	ydot_tmp = &get_new_data(Ydot_Type);
    }

    if (hack_nochem)
    {
        if (ydot_tmp)
        {
            ydot_tmp->setVal(0,dCompYdot,nspecies);
        }
    }
    else
    {
        Real p_amb;
        FORT_GETPAMB(&p_amb);
        const Real Patm = p_amb / P1atm_MKS;

        {
            FArrayBox tmp;
            for (MFIter Smfi(mf); Smfi.isValid(); ++Smfi)
            {
                tmp.resize(Smfi.validbox(),1);
                tmp.copy(mf[Smfi],rho_comp,0,1);
                tmp.invert(1);

                for (int comp = 0; comp < nspecies; ++comp)
                    mf[Smfi].mult(tmp,0,ycomp+comp,1);
            }
        }

        if (ydot_tmp) 
            ydot_tmp->copy(mf,ycomp,dCompYdot,nspecies);

        FArrayBox* chemDiag = 0;

        if (do_not_use_funccount)
        {
	    MultiFab tmp;

            tmp.define(mf.boxArray(), 1, 0, mf.DistributionMap(), Fab_allocate);

            for (MFIter Smfi(mf); Smfi.isValid(); ++Smfi)
            {
                FArrayBox& fb = mf[Smfi];
                const Box& bx = Smfi.validbox();
		FArrayBox& fc = tmp[Smfi];

                if (plot_reactions &&
                    BoxLib::intersect(mf.boxArray(),auxDiag["REACTIONS"]->boxArray()).size() != 0)
                {
                    chemDiag = &( (*auxDiag["REACTIONS"])[Smfi] );
                }

                bool ok = getChemSolve().solveTransient(fb,fb,fb,fb,fc,bx,ycomp,Tcomp,0.5*dt,Patm,chemDiag);

		if (!ok) {
		  BoxLib::Abort("ChemDriver::solveTransient failed");
		}

            }
            //
            // When ngrow>0 this does NOT properly update FuncCount_Type since parallel
            // copy()s do not touch ghost cells.  We'll ignore this since we're not using
            // the FuncCount_Type anyway.
            //
	    get_new_data(FuncCount_Type).copy(tmp);
        }
        else
        {
            //
            // Let's chop the grids up a bit.
            // We want to try and level out the chemistry work.
            //
            const int NProcs = ParallelDescriptor::NProcs();
            BoxArray  ba     = mf.boxArray();
            bool      done   = (ba.size() >= 3*NProcs);

            for (int cnt = 1; !done; cnt *= 2)
            {
                const int ChunkSize = parent->maxGridSize(level)/cnt;

                if (ChunkSize < 16)
                    //
                    // Don't let grids get too small. 
                    //
                    break;

                IntVect chunk(D_DECL(ChunkSize,ChunkSize,ChunkSize));

                for (int j = BL_SPACEDIM-1; j >=0  && ba.size() < 3*NProcs; j--)
                {
                    chunk[j] /= 2;
                    ba.maxSize(chunk);
                    if (ba.size() >= 3*NProcs) done = true;
                }
            }

            DistributionMapping dm = getFuncCountDM(ba,ngrow);

            MultiFab tmp, fcnCntTemp;

            tmp.define(ba, mf.nComp(), 0, dm, Fab_allocate);
            fcnCntTemp.define(ba, 1, 0, dm, Fab_allocate);

            MultiFab diagTemp;
            const bool do_diag = plot_reactions && BoxLib::intersect(ba,auxDiag["REACTIONS"]->boxArray()).size() != 0;
            if (do_diag)
            {
                diagTemp.define(ba, auxDiag["REACTIONS"]->nComp(), 0, dm, Fab_allocate);
                diagTemp.copy(*auxDiag["REACTIONS"]); // Parallel copy
            }

            if (verbose && ParallelDescriptor::IOProcessor())
                std::cout << "*** strang_chem: FABs in tmp MF: " << tmp.size() << '\n';

            tmp.copy(mf); // Parallel copy.

            for (MFIter Smfi(tmp); Smfi.isValid(); ++Smfi)
            {
                FArrayBox& fb = tmp[Smfi];
                const Box& bx = Smfi.validbox();
                FArrayBox& fc = fcnCntTemp[Smfi];
                chemDiag = (do_diag ? &(diagTemp[Smfi]) : 0);

                bool ok = getChemSolve().solveTransient(fb,fb,fb,fb,fc,bx,ycomp,Tcomp,0.5*dt,Patm,chemDiag);

		if (!ok) {
		  BoxLib::Abort("ChemDriver::solveTransient failed");
		}

            }

            mf.copy(tmp); // Parallel copy.

            if (do_diag)
            {
                auxDiag["REACTIONS"]->copy(diagTemp); // Parallel copy
            }

	    MultiFab& FC = get_new_data(FuncCount_Type);
	    FC.copy(fcnCntTemp,0,0,1,0,std::min(ngrow,FC.nGrow()));
        }

        if (ydot_tmp)
        {
            for (MFIter Smfi(mf); Smfi.isValid(); ++Smfi)
            {
                (*ydot_tmp)[Smfi].minus(mf[Smfi], Smfi.validbox(), ycomp, 0, nspecies);
                (*ydot_tmp)[Smfi].mult(1/(0.5*dt), Smfi.validbox(), 0, nspecies);
            }
        }

        for (MFIter Smfi(mf); Smfi.isValid(); ++Smfi) {
#ifdef DO_JBB_HACK_POST
            const Box& box = Smfi.validbox();
            getChemSolve().getHmixGivenTY(mf[Smfi],mf[Smfi],mf[Smfi],box,Temp,first_spec,RhoH);
            const Real Patm = p_amb / P1atm_MKS;
            getChemSolve().getRhoGivenPTY(mf[Smfi],Patm,mf[Smfi],mf[Smfi],box,Temp,first_spec,Density);
            mf[Smfi].mult(mf[Smfi],box,Density,RhoH,1);
#endif
            for (int comp = 0; comp < nspecies; ++comp) {
                mf[Smfi].mult(mf[Smfi],Smfi.validbox(),rho_comp,ycomp+comp,1);
            }
        }

        if (Ydot_action == HT_ImproveYdotOld)
        {
            BL_ASSERT(ydot_tmp != 0);

            MultiFab& ydot_old = get_old_data(Ydot_Type);

            for (MFIter Ymfi(*ydot_tmp); Ymfi.isValid(); ++Ymfi)
            {
                ydot_old[Ymfi].plus((*ydot_tmp)[Ymfi],Ymfi.validbox(),0,dCompYdot,nspecies);
                ydot_old[Ymfi].mult(0.5,Ymfi.validbox(),dCompYdot,nspecies);
            }
        }
    }

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::strang_chem time: " << run_time << '\n';
    }
}

void
HeatTransfer::compute_edge_states (Real              dt,
                                   std::vector<int>* state_comps_to_compute)
{
    BL_PROFILE("HeatTransfer::compute_edge_states()");
    //
    // Compute edge states, store internally.  Do this to
    // avoid recomputing these, and to allow inter-equation consistency.  Note,
    // storage order in EdgeState same as in State_Type.
    // NOTE: Ordering is important here, must do rho.Y and Temp BEFORE RhoH and
    //       Density, but then it doesn't matter.
    //
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... computing edge states\n";
    //
    // Get simulation parameters.
    //
    const Real* dx             = geom.CellSize();
    const Real  prev_time      = state[State_Type].prevTime();
    const Real  prev_pres_time = state[Press_Type].prevTime();
    //
    // NOTE: Massive memory bloat here...how many MultiFab can we waste??
    // Get viscous forcing on the ones we need, don't bother allocating others
    //
    const int nState = desc_lst[State_Type].nComp();
    PArray<MultiFab> visc_terms(nState,PArrayManage);

    const int use_forces_in_trans = godunov->useForcesInTrans();
    //
    // "do_predict" states are predicted normally, after special states
    //
    std::vector<int> do_predict(nState,true);

    if (do_mom_diff != 1) 
    {
        for (int d=0; d<BL_SPACEDIM; ++d)
            do_predict[Xvel+d] = false;
    } 
    for (int sigma=first_spec; sigma<=last_spec; ++sigma)
        do_predict[sigma] = false; // Will get these in a special way
    
    do_predict[RhoH] = false; // Current thought is that we always predict Temp, RhoY and eval Rho, RhoH

    if (do_set_rho_to_species_sum)
    {
        do_predict[Density] = false;
    }

    //
    // This logic and the associated array passed in allows the computation
    // of the edge states to be shut off for specific components.  This is
    // intended to allow special components, such as RhoK and RhoEps in
    // TurbHT be treated differently.  This logic tries to insure that
    // all components with interdependencies are turned on at the same time.
    //
    std::vector<int> compute_comp(nState, true);

    if (state_comps_to_compute != 0)
    {
        BL_ASSERT(state_comps_to_compute->size() == nState);

        for (int cmp = 0; cmp < nState; cmp++)
            compute_comp[cmp] = (*state_comps_to_compute)[cmp];

        if (compute_comp[Density] || compute_comp[Temp] ||
            compute_comp[RhoH]    || compute_comp[first_spec])
        {
            BL_ASSERT(compute_comp[Density]);
            BL_ASSERT(compute_comp[Temp]);
            BL_ASSERT(compute_comp[RhoH]);

            for (int sigma=first_spec; sigma<=last_spec; ++sigma)
                BL_ASSERT(compute_comp[sigma]);
        }
    }
    //
    // If !do_predict, but will need visc terms, get them explicity here
    //
    const int nGrowF = 1;

    MultiFab divu_fp(grids,1,nGrowF);

    create_mac_rhs(divu_fp,nGrowF,prev_time,dt);

    MultiFab Gp;

    if (use_forces_in_trans || (do_mom_diff == 1))
    {
        visc_terms.set(Xvel, new MultiFab(grids,BL_SPACEDIM,nGrowF));
        getViscTerms(visc_terms[Xvel],Xvel,BL_SPACEDIM,prev_time);

        Gp.define(grids,BL_SPACEDIM,1,Fab_allocate);
        getGradP(Gp, prev_pres_time);
    }

    if (compute_comp[first_spec])
    {
	visc_terms.set(first_spec, new MultiFab(grids,nspecies,nGrowF));
	getViscTerms(visc_terms[first_spec],first_spec,nspecies,prev_time);
    }
    //
    // Get all the normal visc terms for everything but velocity
    //
    for (int sigma=BL_SPACEDIM; sigma<nState; ++sigma)
    {
        if (do_predict[sigma] && compute_comp[sigma])
        {
            BL_ASSERT( sigma < first_spec || sigma > last_spec );
            visc_terms.set(sigma, new MultiFab(grids,1,nGrowF));
            if (be_cn_theta == 1.0)
            {
                visc_terms[sigma].setVal(0.0,0,1,nGrowF);
            }
            else
            {
                getViscTerms(visc_terms[sigma],sigma,1,prev_time);
            }
        }
    }
    //
    // Loop on grids, and compute edge fluxes
    //
    FArrayBox edge[BL_SPACEDIM],Rho,U,state,tforces,tvelforces,vel,h,spec;
    //
    // FillPatch'd state data.
    //
    for (FillPatchIterator S_fpi(*this,divu_fp,Godunov::hypgrow(),prev_time,State_Type,0,nState);
         S_fpi.isValid();
         ++S_fpi)
    {
        //
        // Gonna need this array on a per-grid basis
        //
        std::vector<int> this_edge_state_computed(nState,false);

        const int i = S_fpi.index();

        Rho.resize(S_fpi().box(),1);
        U.resize(S_fpi().box(),BL_SPACEDIM);

        Rho.copy(S_fpi(),Density,0,1);
        U.copy(S_fpi(),Xvel,0,BL_SPACEDIM);
        //
        // Get the spec forces based on CC data (forces on EC data in getViscTerms)
        //
        if (use_forces_in_trans || (do_mom_diff == 1))
        {
            NavierStokesBase::getForce(tvelforces,i,nGrowF,Xvel,BL_SPACEDIM,
#ifdef GENGETFORCE
				       prev_time,
#endif		 
				       Rho);
            godunov->Sum_tf_gp_visc(tvelforces,visc_terms[Xvel][S_fpi],Gp[S_fpi],Rho);
        }
        //
        // Set up the workspace for the godunov Box (also resize "edge" for later)
        //
        Array<int> u_bc[BL_SPACEDIM];
        D_TERM(u_bc[0] = getBCArray(State_Type,i,0,1);,
               u_bc[1] = getBCArray(State_Type,i,1,1);,
               u_bc[2] = getBCArray(State_Type,i,2,1);)

        godunov->Setup(grids[i], dx, dt, 0,
                       edge[0], u_bc[0].dataPtr(),
                       edge[1], u_bc[1].dataPtr(),
#if (BL_SPACEDIM == 3)
                       edge[2], u_bc[2].dataPtr(),
#endif
                       U, Rho, tvelforces);

        const int velpred = 0; // Already have edge velocities for transverse derivative

        if (do_mom_diff == 1) 
        {
            vel.resize(S_fpi().box(),BL_SPACEDIM);

            vel.copy(S_fpi(),0,0,BL_SPACEDIM);
            //
            // Loop over the velocity components.
            //
            for (int comp = 0 ; comp < BL_SPACEDIM ; comp++ )
            {
                if (predict_mom_together == 1) 
                {
                    vel.mult(Rho,S_fpi().box(),S_fpi().box(),0,comp,1);
                    tvelforces.mult(Rho,tvelforces.box(),tvelforces.box(),0,comp,1);
                }
                Array<int> bc = getBCArray(State_Type,i,comp,1);

                int iconserv_dummy = 0;
                FArrayBox divu_dummy;

                godunov->edge_states(grids[i], dx, dt, velpred,
                                     u_mac[0][S_fpi], edge[0],
                                     u_mac[1][S_fpi], edge[1],
#if (BL_SPACEDIM == 3)             
                                     u_mac[2][S_fpi], edge[2],
#endif
                                     U,vel,tvelforces,divu_dummy,
                                     comp,comp,bc.dataPtr(),
                                     iconserv_dummy,PRE_MAC);

                for (int d=0; d<BL_SPACEDIM; ++d)
                    (*EdgeState[d])[S_fpi].copy(edge[d],0,comp,1);

                this_edge_state_computed[comp] = true;
            }
        }
        //
        // Get spec edge states.
        // FIXME: Fab copy reqd, force sum below pulls state and forces from same comp
        //
        if (compute_comp[first_spec])
        {
            spec.resize(S_fpi().box(),nspecies);
            spec.copy(S_fpi(),first_spec,0,nspecies);

            NavierStokesBase::getForce(tforces,i,nGrowF,first_spec,nspecies,
#ifdef GENGETFORCE
				       prev_time,
#endif		 
				       Rho);

            for (int comp = 0 ; comp < nspecies ; comp++)
            {
                int state_ind = first_spec + comp;
                int use_conserv_diff = 
                      (advectionType[state_ind] == Conservative) ? true : false;
                Array<int> bc = getBCArray(State_Type,i,state_ind,1);

                AdvectionScheme adv_scheme = FPU;
                if (adv_scheme == PRE_MAC)
                {
                    godunov->Sum_tf_divu_visc(spec, tforces, comp, 1,
                                              visc_terms[first_spec][S_fpi], comp,
                                              divu_fp[S_fpi], Rho, use_conserv_diff);
                    
                    int iconserv_dummy = 0;
                    godunov->edge_states(grids[i], dx, dt, velpred,
                                         u_mac[0][S_fpi], edge[0],
                                         u_mac[1][S_fpi], edge[1],
#if (BL_SPACEDIM==3)
                                         u_mac[2][S_fpi], edge[2],
#endif
                                         U,spec,tforces,divu_fp[S_fpi],
                                         comp,state_ind,bc.dataPtr(),
                                         iconserv_dummy,PRE_MAC);
                }
                else
                {
                    FArrayBox junkDivu(tforces.box(),1);
                    junkDivu.setVal(0.);
                    godunov->Sum_tf_divu_visc(spec, tforces, comp, 1,
                                              visc_terms[first_spec][S_fpi], comp,
                                              junkDivu, Rho, use_conserv_diff);
                    
                    godunov->edge_states(grids[i], dx, dt, velpred,
                                         u_mac[0][S_fpi], edge[0],
                                         u_mac[1][S_fpi], edge[1],
#if (BL_SPACEDIM==3)
                                         u_mac[2][S_fpi], edge[2],
#endif
                                         U,spec,tforces,divu_fp[S_fpi],
                                         comp,state_ind,bc.dataPtr(), 
                                         use_conserv_diff,FPU);
                }

                for (int d=0; d<BL_SPACEDIM; ++d)
                    (*EdgeState[d])[S_fpi].copy(edge[d],0,state_ind,1);

                this_edge_state_computed[state_ind] = true;
            }
        }
        //
        // Get density edge states
        //
        if (compute_comp[Density])
        {
            if (do_set_rho_to_species_sum)
            {
                for (int d=0; d<BL_SPACEDIM; ++d)
                {
                    (*EdgeState[d])[S_fpi].setVal(0.0,edge[d].box(),Density,1);
                    for (int sigma=first_spec; sigma<=last_spec; ++sigma)
                        (*EdgeState[d])[S_fpi].plus((*EdgeState[d])[S_fpi],
                                                edge[d].box(), sigma,Density,1);
                }
                this_edge_state_computed[Density] = true;
            }
            else
            {
                BoxLib::Error("No code yet for rho != sum(rho.Y)");
            }

            if (do_mom_diff == 1 && predict_mom_together == 0)
               for (int icomp = 0; icomp < BL_SPACEDIM; icomp++)
                  for (int d=0; d<BL_SPACEDIM; ++d)
                    (*EdgeState[d])[S_fpi].mult((*EdgeState[d])[S_fpi],(*EdgeState[d])[S_fpi].box(),
                                            (*EdgeState[d])[S_fpi].box(),Density,icomp,1);
        }

        if (compute_comp[Temp])
        {
            //
            // Get Temp edge states via extrap.
            //
            const int comp = 0;
            const int state_ind = Temp;
            int use_conserv_diff = (advectionType[state_ind] == Conservative) ? true : false;
            state.resize(S_fpi().box(),1);
            state.copy(S_fpi(),state_ind,0,1);
            FArrayBox& vt = visc_terms[state_ind][S_fpi];
            int vtComp = 0;

            NavierStokesBase::getForce(tforces,i,nGrowF,state_ind,1,
#ifdef GENGETFORCE
				       prev_time,
#endif		 
				       Rho);
            Array<int> bc = getBCArray(State_Type,i,state_ind,1);

            AdvectionScheme adv_scheme = FPU;

            if (adv_scheme == PRE_MAC)
            {
                godunov->Sum_tf_divu_visc(state, tforces,  comp, 1,
                                          vt, vtComp,
                                          divu_fp[S_fpi], Rho, use_conserv_diff);

                int iconserv_dummy = 0;
                godunov->edge_states(grids[i], dx, dt, velpred,
                                     u_mac[0][S_fpi], edge[0],
                                     u_mac[1][S_fpi], edge[1],
#if (BL_SPACEDIM==3)
                                     u_mac[2][S_fpi], edge[2],
#endif
                                     U, state, tforces, divu_fp[S_fpi],
                                     comp, state_ind, bc.dataPtr(),
                                     iconserv_dummy, PRE_MAC);

            }
            else
            {
                FArrayBox junkDivu(tforces.box(),1);
                junkDivu.setVal(0);
                godunov->Sum_tf_divu_visc(state, tforces,  comp, 1,
                                          vt, vtComp,
                                          junkDivu, Rho, use_conserv_diff);

                godunov->edge_states(grids[i], dx, dt, velpred,
                                     u_mac[0][S_fpi], edge[0],
                                     u_mac[1][S_fpi], edge[1],
#if (BL_SPACEDIM==3)
                                     u_mac[2][S_fpi], edge[2],
#endif
                                     U, state, tforces, divu_fp[S_fpi],
                                     comp, state_ind, bc.dataPtr(), 
                                     use_conserv_diff, FPU);
            }

            for (int d=0; d<BL_SPACEDIM; ++d)
                (*EdgeState[d])[S_fpi].copy(edge[d],0,state_ind,1);

            this_edge_state_computed[state_ind] = true;
        }

        if (compute_comp[RhoH])
        {
            //
            // Set rhoh on edges = sum(rho.Y.H)
            //
            for (int d=0; d<BL_SPACEDIM; ++d)
            {
                (*EdgeState[d])[S_fpi].setVal(0.0,edge[d].box(),RhoH,1);
                h.resize(edge[d].box(),nspecies);
                getChemSolve().getHGivenT(h,(*EdgeState[d])[S_fpi],
                                          edge[d].box(),Temp,0);
                h.mult((*EdgeState[d])[S_fpi],edge[d].box(),first_spec,0,
                       nspecies);
                
                (*EdgeState[d])[S_fpi].setVal(0.0,edge[d].box(),RhoH,1);
                for (int comp=0; comp<nspecies; ++comp)
                    (*EdgeState[d])[S_fpi].plus(h,edge[d].box(),comp,RhoH,1);
            }
            this_edge_state_computed[RhoH] = true;
        }
        //
        // Now do the rest as normal
        //
        state.resize(S_fpi().box(),1);

        for (int state_ind=0; state_ind<nState; ++state_ind)
        {
            if (do_predict[state_ind]                &&
                !this_edge_state_computed[state_ind] &&
                compute_comp[state_ind])
            {
                int use_conserv_diff =
                    (advectionType[state_ind] == Conservative) ? true : false;
                //
                // Do it the old-fashioned way.
                //
                state.copy(S_fpi(),state_ind,0,1);
                const int comp = 0;
                NavierStokesBase::getForce(tforces,i,nGrowF,state_ind,1,
#ifdef GENGETFORCE
					   prev_time,
#endif		 
					   Rho);
                godunov->Sum_tf_divu_visc(state, tforces, comp, 1,
                                          visc_terms[state_ind][S_fpi], 0,
                                          divu_fp[S_fpi], Rho,
                                          use_conserv_diff);
                Array<int> bc = getBCArray(State_Type,i,state_ind,1);
                int iconserv_dummy = 0;
                godunov->edge_states(grids[i], dx, dt, velpred,
                                     u_mac[0][S_fpi], edge[0],
                                     u_mac[1][S_fpi], edge[1],
#if (BL_SPACEDIM==3)
                                     u_mac[2][S_fpi], edge[2],
#endif
                                     U,state,tforces,divu_fp[S_fpi],
                                     comp,state_ind,bc.dataPtr(),
                                     iconserv_dummy,PRE_MAC);

                for (int d=0; d<BL_SPACEDIM; ++d)
                    (*EdgeState[d])[S_fpi].copy(edge[d],0,state_ind,1);

                this_edge_state_computed[state_ind] = true;
            }
        }
    }
}

void
HeatTransfer::momentum_advection (Real dt, bool do_adv_reflux)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... advect momentum\n";

    const int finest_level = parent->finestLevel();

    FArrayBox edge[BL_SPACEDIM];
    //
    // Compute the advective forcing for momentum.
    //
    const int use_conserv_diff = true;

    MultiFab fluxes[BL_SPACEDIM];

    if (do_reflux && level < parent->finestLevel())
    {
        for (int d=0; d<BL_SPACEDIM; ++d)
            fluxes[d].define((*EdgeState[d]).boxArray(), BL_SPACEDIM, 0, Fab_allocate);
    }

    for (MFIter AofS_mfi(*aofs); AofS_mfi.isValid(); ++AofS_mfi)
    {
        const int i = AofS_mfi.index();

        for (int d=0; d<BL_SPACEDIM; ++d)
            edge[d].resize(BoxLib::surroundingNodes(grids[i],d),1);

        for (int comp = 0 ; comp < BL_SPACEDIM ; comp++ )
        {
            // 
            // If here, edge states at n+1/2 have already been computed, get a copy
            // 
            for (int d=0; d<BL_SPACEDIM; ++d)
                edge[d].copy((*EdgeState[d])[AofS_mfi],edge[d].box(),comp,edge[d].box(),0,1);
 
            godunov->ComputeAofs(grids[i],
                                 area[0][AofS_mfi],u_mac[0][AofS_mfi],edge[0],
                                 area[1][AofS_mfi],u_mac[1][AofS_mfi],edge[1],
#if BL_SPACEDIM==3
                                 area[2][AofS_mfi],u_mac[2][AofS_mfi],edge[2],
#endif
                                 volume[AofS_mfi],(*aofs)[AofS_mfi],comp,
                                 use_conserv_diff);
            if (do_adv_reflux)
            {
                if (level < parent->finestLevel())
                {
                    for (int d=0; d<BL_SPACEDIM; ++d)
                        fluxes[d][AofS_mfi].copy(edge[d],0,comp,1);
                }

                if (level > 0)
                {
                    for (int d=0; d<BL_SPACEDIM; ++d)
                        advflux_reg->FineAdd(edge[d],d,i,0,comp,1,dt);
                }
            }
        }
    }

    D_TERM(edge[0].clear();, edge[1].clear();, edge[2].clear(););

    if (do_adv_reflux && level < finest_level)
    {
        for (int d=0; d<BL_SPACEDIM; ++d)
        {
            getAdvFluxReg(level+1).CrseInit(fluxes[d],d,0,0,BL_SPACEDIM,-dt);
        }
    }
}

void
HeatTransfer::scalar_advection (Real dt,
                                int  fscalar,
                                int  lscalar,
                                bool do_adv_reflux)
{
    BL_PROFILE("HeatTransfer::scalar_advection()");

    if (verbose && benchmarking) ParallelDescriptor::Barrier();

    const Real strt_time = ParallelDescriptor::second();
    //
    // Compute the advection flux divergences
    //
    const bool do_special_rhoh = nspecies>0 
        && do_set_rho_to_species_sum 
        && fscalar>=RhoH && lscalar<=RhoH
        && !unity_Le 
        && do_add_nonunityLe_corr_to_rhoh_adv_flux;
    FluxBoxes fb_fluxNULN;
    MultiFab** fluxNULN;
    //
    // If RhoH included, compute non-unity Lewis number flux addition.stuff
    // (using LinOp stuff which works on MultiFabs, so need to do this prior
    // to the subsequent MFIter loop)
    //
    // Note that this requires reasonable species and Temperature values
    //  in the state at np1, and pulls transport coeffs from the level. 
    //  It should be that the coeffs in the level have been computed with the
    //  temp and species in the state (to preserve isothermal flow over time step)
    //
    if (do_special_rhoh)
    {
        MultiFab Soln(grids,1,1);
        const Real prev_time = state[State_Type].prevTime();
        const Real cur_time  = state[State_Type].curTime();
        MultiFab& S_new = get_new_data(State_Type);
        MultiFab& S_old = get_old_data(State_Type);
        //
        // For each species, make a LinOp that can compute -(lambda/cp).A.Grad(Y_l)
        // (the ViscBndry construction knows how to fill bc's correctly, so
        //    re-construct linop/bndrydata for each species)
        //
        const Real a = 1.0;     // Passed around, but not used
        Real rhsscale;          //  -ditto-

	// 1 = assumes d(c)/dt ~ grad c
	// 2 = assumes d(rho*c)/dt ~ grad c
	// 3 = assumes rho*d(c)/dt ~ grad c
        const int rho_flag = 2; 

        MultiFab *alpha=0;      //  -ditto-
	FluxBoxes fb_fluxSC(this, 1, 0); // tmp single-component lambda/cp grad h
        MultiFab   **fluxSC = fb_fluxSC.get();
	FluxBoxes fb_fluxi(this, nspecies, 0); // all the species fluxes, "D"
	MultiFab   **fluxi = fb_fluxi.get(); 
	FluxBoxes fb_rhoh_visc(this, 1, 0); // lambda/cp
	MultiFab   **rhoh_visc = fb_rhoh_visc.get();
	fluxNULN = fb_fluxNULN.define(this, 1, 0); // single-component "DD" or "NULN" fluxes

        const int nGrow    = 1; // Size to grow fil-patched fab for T below
        const int dataComp = 0; // coeffs loaded into 0-comp for all species
        //
        // Initialize fluxNULN (NULN = non-unity Lewis number)
        //
        for (int d = 0; d < BL_SPACEDIM; ++d)
            fluxNULN[d]->setVal(0.0);
        //
        // Get the NULN flux contrib from n data
        //

	// put lambda/cp into rhoh_visc
        getDiffusivity(rhoh_visc, prev_time, RhoH, 0, 1);

	// rho^n+1/2, don't think this is needed
        const MultiFab& RhoHalftime = get_rho_half_time();

        for (int comp = 0; comp < nspecies; ++comp)
        {
            const Real b     = 1.0 - be_cn_theta;
            const int  sigma = first_spec + comp;
            //
            // Start by getting lambda/cp.Grad(Y) (note: neg of usual diff flux)
            //
            ViscBndry      visc_bndry;
            ABecLaplacian* visc_op;

            visc_op = diffusion->getViscOp(sigma,a,b,prev_time,visc_bndry,
                                           RhoHalftime,rho_flag,&rhsscale,
                                           rhoh_visc,dataComp,alpha,0);

            visc_op->maxOrder(diffusion->maxOrder());
            MultiFab::Copy(Soln,S_old,sigma,0,1,0);

            for (MFIter Smfi(Soln); Smfi.isValid(); ++Smfi)
                Soln[Smfi].divide(S_old[Smfi],Smfi.validbox(),Density,0,1);

	    // compute the lamba/cp grad Y term
            visc_op->compFlux(D_DECL(*fluxSC[0],*fluxSC[1],*fluxSC[2]),Soln);
            for (int d=0; d < BL_SPACEDIM; ++d)
                fluxSC[d]->mult(-b/geom.CellSize()[d]);
            //
            // Here, get fluxi = (lambda/cp - rho.D)Grad(Y)
            //                 = lambda/cp.Grad(Y) + SpecDiffFlux
            //
            for (int d = 0; d < BL_SPACEDIM; ++d)
            {
                for (MFIter SDF_mfi(*SpecDiffusionFluxn[d]); SDF_mfi.isValid(); ++SDF_mfi)
                {
                    const Box& ebox    = SDF_mfi.validbox();
                    FArrayBox& SDF_fab = (*SpecDiffusionFluxn[d])[SDF_mfi];
		    // add the Gamma_m term
                    (*fluxi[d])[SDF_mfi].copy(SDF_fab,ebox,comp,ebox,comp,1);
		    // add the lambda/cp grad Y term
                    (*fluxi[d])[SDF_mfi].plus((*fluxSC[d])[SDF_mfi],ebox,0,comp,1);
                }
            }
            delete visc_op;
        }

        FArrayBox eTemp, h;
        //
        // Multiply fluxi by h_i, and add to running total.
        //
        for (FillPatchIterator Told_fpi(*this,S_old,nGrow,prev_time,State_Type,Temp,1);
             Told_fpi.isValid();
             ++Told_fpi)
        {
            const Box& box = Told_fpi.validbox();

            for (int d = 0; d < BL_SPACEDIM; ++d)
            {
                const Box& ebox = BoxLib::surroundingNodes(box,d);
                eTemp.resize(ebox,1);
                FPLoc bc_lo = fpi_phys_loc(get_desc_lst()[State_Type].getBC(Temp).lo(d));
                FPLoc bc_hi = fpi_phys_loc(get_desc_lst()[State_Type].getBC(Temp).hi(d));
                
                center_to_edge_fancy(Told_fpi(),eTemp,
                                     BoxLib::grow(box,BoxLib::BASISV(d)),0,0,1,
                                     geom.Domain(),bc_lo,bc_hi);
                
                h.resize(ebox,nspecies);
                getChemSolve().getHGivenT(h,eTemp,ebox,0,0);
                (*fluxi[d])[Told_fpi].mult(h,ebox,0,0,nspecies);
                
                for (int comp=0; comp<nspecies; ++comp)
                    (*fluxNULN[d])[Told_fpi].plus((*fluxi[d])[Told_fpi],ebox,comp,0,1);
            }
        }
        //
        // Get the Le!=1 flux contrib from n+1 data.
        //
        getDiffusivity(rhoh_visc, cur_time, RhoH, 0, 1);

        for (int comp = 0; comp < nspecies; ++comp)
        {
            const Real b     = be_cn_theta;
            const int  sigma = first_spec + comp;
            //
            //  start by getting lambda/cp.Grad(Y) (note: neg of usual diff flux)
            //
            ViscBndry      visc_bndry;
            ABecLaplacian* visc_op;

            visc_op = diffusion->getViscOp(sigma,a,b,cur_time,visc_bndry,
                                           RhoHalftime,rho_flag,&rhsscale,
                                           rhoh_visc,dataComp,alpha,0);

            visc_op->maxOrder(diffusion->maxOrder());
            MultiFab::Copy(Soln,S_new,sigma,0,1,0);

            for (MFIter Smfi(Soln); Smfi.isValid(); ++Smfi)
                Soln[Smfi].divide(S_new[Smfi],Smfi.validbox(),Density,0,1);

	    // compute the lamba/cp grad Y term
            visc_op->compFlux(D_DECL(*fluxSC[0],*fluxSC[1],*fluxSC[2]),Soln);
            for (int d=0; d < BL_SPACEDIM; ++d)
                fluxSC[d]->mult(-b/geom.CellSize()[d]);
            //
            // Here, get fluxi = (lambda/cp - rho.D)Grad(Y)
            //                 = lambda/cp.Grad(Y) + SpecDiffFlux
            //
            for (int d = 0; d < BL_SPACEDIM; ++d)
            {
                MFIter SDF_mfi(*SpecDiffusionFluxnp1[d]);
                for ( ; SDF_mfi.isValid(); ++SDF_mfi)
                {
                    FArrayBox& SDF_fab = (*SpecDiffusionFluxnp1[d])[SDF_mfi];
                    const Box& ebox    = SDF_mfi.validbox();
		    // add the Gamma_m term
                    (*fluxi[d])[SDF_mfi].copy(SDF_fab,ebox,comp,ebox,comp,1);
		    // add the lambda/cp grad Y term
                    (*fluxi[d])[SDF_mfi].plus((*fluxSC[d])[SDF_mfi],ebox,0,comp,1);
                }
            }
            delete visc_op;
        }

        Soln.clear();
        fb_fluxSC.clear();
        fb_rhoh_visc.clear();
        //
        // Multiply fluxi by h_i, and add to running total
        //
        for (FillPatchIterator Tnew_fpi(*this,S_new,nGrow,cur_time,State_Type,Temp,1);
             Tnew_fpi.isValid();
             ++Tnew_fpi)
        {
            const Box& box = Tnew_fpi.validbox();

            for (int d = 0; d < BL_SPACEDIM; ++d)
            {
                const Box& ebox = BoxLib::surroundingNodes(box,d);
                eTemp.resize(ebox,1);
                FPLoc bc_lo = fpi_phys_loc(get_desc_lst()[State_Type].getBC(Temp).lo(d));
                FPLoc bc_hi = fpi_phys_loc(get_desc_lst()[State_Type].getBC(Temp).hi(d));
                
                center_to_edge_fancy(Tnew_fpi(),eTemp,BoxLib::grow(box,BoxLib::BASISV(d)),
                                     0,0,1,geom.Domain(),bc_lo,bc_hi);
                
                h.resize(ebox,nspecies);
                getChemSolve().getHGivenT(h,eTemp,ebox,0,0);
                (*fluxi[d])[Tnew_fpi].mult(h,ebox,0,0,nspecies);
                
                for (int comp = 0; comp < nspecies; ++comp)
                    (*fluxNULN[d])[Tnew_fpi].plus((*fluxi[d])[Tnew_fpi],ebox,comp,0,1);
            }
        }
    }

    FArrayBox edge[BL_SPACEDIM];

    const int nscalar = lscalar - fscalar + 1;

    MultiFab fluxes[BL_SPACEDIM];

    if (do_adv_reflux && level < parent->finestLevel())
    {
        for (int d=0; d<BL_SPACEDIM; ++d)
            fluxes[d].define((*EdgeState[d]).boxArray(), nscalar, 0, Fab_allocate);
    }

    for (MFIter AofS_mfi(*aofs); AofS_mfi.isValid(); ++AofS_mfi)
    {
        const int i = AofS_mfi.index();

        for (int d=0; d<BL_SPACEDIM; ++d)
            edge[d].resize(BoxLib::surroundingNodes(grids[i],d),1);

        for (int sigma=fscalar; sigma<=lscalar; ++sigma)
        {
            // 
            // If here, edge states at n+1/2 have already been computed, get a copy
            // 
            for (int d=0; d<BL_SPACEDIM; ++d)
                edge[d].copy((*EdgeState[d])[AofS_mfi],edge[d].box(),sigma,edge[d].box(),0,1);
            
            int use_conserv_diff = (advectionType[sigma] == Conservative) ? true : false;

	    // takes edge states, multiplies them by umac and area
	    // returns aofs
	    // note that "edge" gets converted to a "umac*area*edge"
            godunov->ComputeAofs(grids[i],
                                 area[0][AofS_mfi],u_mac[0][AofS_mfi],edge[0],
                                 area[1][AofS_mfi],u_mac[1][AofS_mfi],edge[1],
#if BL_SPACEDIM==3
                                 area[2][AofS_mfi],u_mac[2][AofS_mfi],edge[2],
#endif
                                 volume[AofS_mfi],(*aofs)[AofS_mfi],sigma,
                                 use_conserv_diff);
            //
            // Add divergence of fluxNULN to aofs[RhoH], and increment advective
            // going into flux registers
            //
            if (sigma==RhoH && do_special_rhoh)
            {
                if (do_adv_reflux)
                    for (int d=0; d<BL_SPACEDIM; ++d)
                        edge[d].plus((*fluxNULN[d])[AofS_mfi],edge[d].box(),0,0,1);
                
                FArrayBox& staten = (*aofs)[AofS_mfi];
                const FArrayBox& stateo = staten;
                const Box& box = AofS_mfi.validbox();
                const FArrayBox& vol = volume[AofS_mfi];
                const Real mult = 1.0; // no dt scaling of aofs, done in scl_adv_upd
                const int nComp = 1;
                
                FORT_INCRWEXTFLXDIV(box.loVect(), box.hiVect(),
                                    (*fluxNULN[0])[AofS_mfi].dataPtr(),
                                    ARLIM((*fluxNULN[0])[AofS_mfi].loVect()),
                                    ARLIM((*fluxNULN[0])[AofS_mfi].hiVect()),
                                    (*fluxNULN[1])[AofS_mfi].dataPtr(),
                                    ARLIM((*fluxNULN[1])[AofS_mfi].loVect()),
                                    ARLIM((*fluxNULN[1])[AofS_mfi].hiVect()),
#if BL_SPACEDIM == 3
                                    (*fluxNULN[2])[AofS_mfi].dataPtr(),
                                    ARLIM((*fluxNULN[2])[AofS_mfi].loVect()),
                                    ARLIM((*fluxNULN[2])[AofS_mfi].hiVect()),
#endif
                                    stateo.dataPtr(RhoH),
                                    ARLIM(stateo.loVect()), ARLIM(stateo.hiVect()),
                                    staten.dataPtr(RhoH),
                                    ARLIM(staten.loVect()), ARLIM(staten.hiVect()),
                                    vol.dataPtr(),
                                    ARLIM(vol.loVect()), ARLIM(vol.hiVect()),
                                    &nComp, &mult);
            }

	    // at this point, rhoh component of edge contains area-weighted
	    // (U^ADV*rho*h)^{n+1/2} - (1/2)sum(hm(gamma_m lamba/cp grad y)^{n+1,*}
	    //                       - (1/2)sum(hm(gamma_m lamba/cp grad y)^n

	    // rhoY component of edge contains area-weighted
	    // (U^ADV*rho*Y_m)^{n+1/2}
            if (do_adv_reflux)
            {
                if (level < parent->finestLevel())
                {
                    const int fcomp = sigma - fscalar;

                    for (int d=0; d<BL_SPACEDIM; ++d)
                        fluxes[d][AofS_mfi].copy(edge[d],0,fcomp,1);
                }
                if (level > 0)
                {
                    for (int d=0; d<BL_SPACEDIM; ++d)
                        advflux_reg->FineAdd(edge[d],d,i,0,sigma,1,dt);
                }
            }
        }
    }

    D_TERM(edge[0].clear();, edge[1].clear();, edge[2].clear(););

    if (do_adv_reflux && level < parent->finestLevel())
    {
        for (int d=0; d<BL_SPACEDIM; ++d)
        {
            getAdvFluxReg(level+1).CrseInit(fluxes[d],d,0,fscalar,nscalar,-dt);
        }
    }

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::scalar_advection(): time: " << run_time << '\n';
    }
}

void
HeatTransfer::rhoh_update (Real time,
                           Real dt,
                           int  corrector) 
{
    //
    // Do implicit c-n solve for RhoH
    //
    scalar_update(dt,RhoH,RhoH,corrector);
}

void
HeatTransfer::temp_update (Real dt,
                           int  corrector) 
{
    //
    // Do implicit c-n solve for temperature.
    //
    scalar_update(dt,Temp,Temp,corrector);
}

void
HeatTransfer::spec_update (Real time,
                           Real dt,
                           int  corrector) 
{
    //
    // Do implicit c-n solve for rho*Y_l, l=0,nspecies-1.
    //
    scalar_advection_update(dt, first_spec, last_spec);
	    
    if (unity_Le)
    {
        scalar_diffusion_update(dt, first_spec, last_spec, corrector);
    }
    else
    {
        differential_spec_diffusion_update(dt, corrector);
    }
    //
    // Enforce sum_l rho U Y_l equals rho.
    //
    if (floor_species)
        scale_species(get_new_data(State_Type),0,1);
}

void
HeatTransfer::tracer_update (Real dt,
                             int  corrector) 
{
    //
    // Update tracer.
    //
    if (have_trac)
        scalar_update(dt,Trac,Trac,corrector);
    if (have_rhort)
        scalar_update(dt,RhoRT,RhoRT,corrector);
}

void
HeatTransfer::scalar_update (Real dt,
                             int  first_scalar, 
                             int  last_scalar,
                             int  corrector)
{
    //
    // Do implicit c-n solve for an arbitrary scalar (i.e., not velocity).
    //
    scalar_advection_update(dt, first_scalar, last_scalar);
    scalar_diffusion_update(dt, first_scalar, last_scalar, corrector);
}

//
// sum_l rho D grad Y_l dot grad h_l on the valid region of the multifab.
//
void
HeatTransfer::compute_rhoDgradYgradH (Real      time,
                                      MultiFab& rdgydgh)
{
    //
    // FIXME: Shouldn't this really be the species diffusion fluxes
    //        dotted with grad(h_i) to be consistent? I dunno....
    //
    // Get edge-centered rho.D
    // (copy spec visc from internal database directly (before inflow "zeroing"))
    //
    FluxBoxes fb_beta(this, nspecies, 0);
    MultiFab** beta = fb_beta.get();
    getDiffusivity(beta, time, first_spec, 0, nspecies);
    //
    // Get result, using cell-centered Y,h and edge-centered rhoD
    //
    FArrayBox rdgydgh_spec_i, tmp, h;

    const Real* dx = geom.CellSize();

    rdgydgh.setVal(0,0);
    //
    // nspecies = number of species, ncomp is one greater due to fillpatching
    // density and species together.
    //
    int nspecies = last_spec - first_spec + 1;

    for (FillPatchIterator rho_and_species_fpi(*this,rdgydgh,1,time,State_Type,Density,nspecies+1),
             Temp_fpi(*this,rdgydgh,1,time,State_Type,Temp,1);
         rho_and_species_fpi.isValid() && Temp_fpi.isValid();
         ++rho_and_species_fpi, ++Temp_fpi)
    {
        const int  i               = rho_and_species_fpi.index();
        const int* lo              = grids[i].loVect();
        const int* hi              = grids[i].hiVect();
        FArrayBox& rho_and_species = rho_and_species_fpi();
        const Box& bx              = rho_and_species.box();

        rdgydgh_spec_i.resize(grids[i],1);
        DEF_LIMITS(rdgydgh_spec_i,prod,prodlo,prodhi);

        const  int* speclo  = rho_and_species.loVect();
        const  int* spechi  = rho_and_species.hiVect();

        tmp.resize(rho_and_species.box(),1);
        tmp.copy(rho_and_species,0,0,1);
        tmp.invert(1);

        h.resize(bx,nspecies);

	getChemSolve().getHGivenT(h,Temp_fpi(),bx,0,0);

        for (int spec = first_spec; spec <= last_spec; spec++) 
        {
            const int comp = spec - first_spec;

            DEF_CLIMITSCOMP(h,hdat,hlo,hhi,comp);
            DEF_CLIMITSCOMP((*beta[0])[rho_and_species_fpi],betax,betaxlo,betaxhi,comp);
            DEF_CLIMITSCOMP((*beta[1])[rho_and_species_fpi],betay,betaylo,betayhi,comp);
#if (BL_SPACEDIM==3)
            DEF_CLIMITSCOMP((*beta[2])[rho_and_species_fpi],betaz,betazlo,betazhi,comp);
#endif

            rho_and_species.mult(tmp,0,comp+1,1);

            const Real* specdat = rho_and_species.dataPtr(comp+1);

            FORT_COMPUTE_RHODGRADHDOTGRADY(dx,lo,hi,
                                           ARLIM(speclo),ARLIM(spechi),specdat,
                                           ARLIM(hlo),ARLIM(hhi),hdat,
                                           ARLIM(betaxlo),ARLIM(betaxhi),betax,
                                           ARLIM(betaylo),ARLIM(betayhi),betay,
#if (BL_SPACEDIM==3) 
                                           ARLIM(betazlo),ARLIM(betazhi),betaz,
#endif            
                                           ARLIM(prodlo),ARLIM(prodhi),prod);

            rdgydgh[rho_and_species_fpi].plus(rdgydgh_spec_i);
        }
    }
}

//
// An enum to clean up mac_sync...questionable usefullness
//
enum SYNC_SCHEME {ReAdvect, UseEdgeState, Other};

void
HeatTransfer::mac_sync ()
{
    BL_PROFILE("HeatTransfer::mac_sync()");

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... mac_sync\n";

    if (verbose && benchmarking) ParallelDescriptor::Barrier();

    const Real strt_time = ParallelDescriptor::second();

    int        sigma;
    const int  finest_level   = parent->finestLevel();
    const int  ngrids         = grids.size();
    const Real prev_time      = state[State_Type].prevTime();
    const Real cur_time       = state[State_Type].curTime();
    const Real prev_pres_time = state[Press_Type].prevTime();
    const Real dt             = parent->dtLevel(level);
    MultiFab&  Rh             = get_rho_half_time();
    //
    // will hold q^{n+1,p} * (delta rho)^sync for conserved quantities
    // as defined before Eq (18) in DayBell:2000.  Note that in the paper, 
    // Eq (18) is missing Y_m^{n+1,p} * (delta rho)^sync in the RHS
    // and Eq (19) is missing the h^{n+1,p} * (delta rho)^sync in the RHS
    //
    MultiFab*  DeltaSsync     = 0;
    //
    // Compute the corrective pressure used to compute U^{ADV,corr} in mac_sync_compute
    //
    mac_projector->mac_sync_solve(level,dt,Rh,fine_ratio);

    if (do_reflux)
    {
        MultiFab& S_new = get_new_data(State_Type);

	Array<SYNC_SCHEME> sync_scheme(NUM_STATE,ReAdvect);

        if (do_mom_diff == 1)
          for (int i=0; i<BL_SPACEDIM; ++i)
            sync_scheme[i] = UseEdgeState;

        for (int i=BL_SPACEDIM; i<NUM_STATE; ++i)
            sync_scheme[i] = UseEdgeState;
        
        Array<int> incr_sync(NUM_STATE,0);
        for (int i=0; i<sync_scheme.size(); ++i)
            if (sync_scheme[i] == ReAdvect)
                incr_sync[i] = 1;

	//
	// After solving for mac_sync_phi in mac_sync_solve(), we
	// can now do the sync advect step in mac_sync_compute().
	// This consists of two steps
	//
	// 1. compute U^{ADV,corr} as the gradient of mac_sync_phi
	// 2. add -D^MAC ( U^{ADV,corr} * rho * q)^{n+1/2} ) to flux registers
	//

	// velocities
        if (do_mom_diff == 0) 
        {
            mac_projector->mac_sync_compute(level,u_mac,Vsync,Ssync,Rh,
                                            (level > 0) ? &getAdvFluxReg(level) : 0,
                                            advectionType,prev_time,
                                            prev_pres_time,dt,NUM_STATE,
                                            be_cn_theta,
                                            modify_reflux_normal_vel,
                                            do_mom_diff,
                                            incr_sync);
        }
        else
        {
            for (int comp=0; comp<BL_SPACEDIM; ++comp)
            {
                if (sync_scheme[comp]==UseEdgeState)
                {
                    mac_projector->mac_sync_compute(level,Vsync,comp,
                                                    comp,EdgeState, comp,Rh,
                                                    (level>0 ? &getAdvFluxReg(level):0),
                                                    advectionType,modify_reflux_normal_vel,dt);
                }
            }
        }

	// scalars
        for (int comp=BL_SPACEDIM; comp<NUM_STATE; ++comp)
        {
            if (sync_scheme[comp]==UseEdgeState)
            {
                int s_ind = comp - BL_SPACEDIM;
                //
                // This routine does a sync advect step for a single 
                // scalar component. The half-time edge states are passed in.
                // This routine is useful when the edge states are computed
                // in a physics-class-specific manner. (For example, as they are
                // in the calculation of div rho U h = div U sum_l (rho Y)_l h_l(T)).
                //
                mac_projector->mac_sync_compute(level,Ssync,comp,s_ind,
                                                EdgeState,comp,Rh,
                                                (level>0 ? &getAdvFluxReg(level):0),
                                                advectionType,modify_reflux_normal_vel,dt);
            }
        }
        
        Ssync.mult(dt,Ssync.nGrow());

        sync_setup(DeltaSsync);

        //
        // For all conservative variables Q (other than density)
	// set DeltaSsync = q^{n+1,p} * (delta rho)^sync
	// subtract q^{n+1,p} * (delta rho)^sync from Ssync
        //
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            const int  i   = mfi.index();
            const Box& grd = grids[i];

            int iconserved = -1;

            for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
            {
                if (istate != Density && advectionType[istate] == Conservative)
                {
                    iconserved++;
             
		    FArrayBox delta_ssync;
                    delta_ssync.resize(grd,1);
                    delta_ssync.copy(S_new[mfi],grd,istate,grd,0,1); // delta_ssync = (rho*q)^{n+1,p}
                    delta_ssync.divide(S_new[mfi],grd,Density,0,1); // delta_ssync = q^{n+1,p}
                    FArrayBox& s_sync = Ssync[mfi]; // Ssync = RHS of Eq (18), (19) without the q^{n+1,p} * (delta rho)^sync terms
                    delta_ssync.mult(s_sync,grd,Density-BL_SPACEDIM,0,1); // delta_ssync = q^{n+1,p} * (delta rho)^sync
                    (*DeltaSsync)[mfi].copy(delta_ssync,grd,0,grd,iconserved,1); // DeltaSsync = q^{n+1,p} * (delta rho)^sync
                    s_sync.minus(delta_ssync,grd,0,istate-BL_SPACEDIM,1); // Ssync = Ssync - q^{n+1,p} * (delta rho)^sync
                }
            }
        }
        //
        // Now, increment density.
        //
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            const int i = mfi.index();
            S_new[mfi].plus(Ssync[mfi],grids[i],Density-BL_SPACEDIM,Density,1);
        }
        make_rho_curr_time();

        const int numscal = NUM_STATE - BL_SPACEDIM;
        //
        // Set do_diffuse_sync to 0 for debugging reasons only.
        //
        if (do_mom_diff == 1)
        {
            for (MFIter Vsyncmfi(Vsync); Vsyncmfi.isValid(); ++Vsyncmfi)
            {
                const int  i    = Vsyncmfi.index();
                const Box& vbox = rho_ctime.box(i);

                D_TERM(Vsync[Vsyncmfi].divide(rho_ctime[Vsyncmfi],vbox,0,Xvel,1);,
                       Vsync[Vsyncmfi].divide(rho_ctime[Vsyncmfi],vbox,0,Yvel,1);,
                       Vsync[Vsyncmfi].divide(rho_ctime[Vsyncmfi],vbox,0,Zvel,1););
            }
        }

        if (do_diffuse_sync)
        {
	    FluxBoxes fb_beta(this);
	    MultiFab** beta = fb_beta.get();
            if (is_diffusive[Xvel])
            {
                int rho_flag = (do_mom_diff == 0) ? 1 : 3;
                getViscosity(beta, cur_time);
                diffusion->diffuse_Vsync(Vsync,dt,be_cn_theta,Rh,rho_flag,beta);
            }
	    
	    if (!unity_Le 
		&& nspecies>0 
		&& do_add_nonunityLe_corr_to_rhoh_adv_flux) 
	    {
		//
		// Diffuse the species syncs such that sum(SpecDiffSyncFluxes) = 0
	        // After exiting, SpecDiffusionFluxnp1 should contain rhoD grad (delta Y)^sync
		// Also, Ssync for species should contain rho^{n+1} * (delta Y)^sync
		differential_spec_diffuse_sync(dt);

                MultiFab Soln(grids,1,1);
                const Real cur_time  = state[State_Type].curTime();
                const Real a = 1.0;     // Passed around, but not used
                Real rhsscale;          //  -ditto-
                const int rho_flag = 2; // FIXME: Messy assumption
                MultiFab *alpha=0;      //  -ditto-
		FluxBoxes fb_fluxSC(this, 1, 0);
                MultiFab **fluxSC = fb_fluxSC.get();
		FluxBoxes fb_fluxNULN(this, nspecies, 0);
		MultiFab **fluxNULN = fb_fluxNULN.get();
		FluxBoxes fb_rhoh_visc(this, 1, 0);
		MultiFab **rhoh_visc = fb_rhoh_visc.get();

                const int nGrow    = 1; // Size to grow fil-patched fab for T below
                const int dataComp = 0; // coeffs loaded into 0-comp for all species
                  
		// get lambda/cp
                getDiffusivity(rhoh_visc, cur_time, RhoH, 0, 1);

		// compute lambda/cp grad (delta Y_m^sync)
                for (int comp = 0; comp < nspecies; ++comp)
                {
                    const Real b     = be_cn_theta;
                    const int  sigma = first_spec + comp;
                    //
                    //  start by getting lambda/cp.Grad(delta Y^sync)
                    //   (note: neg of usual diff flux)
                    //
                    ViscBndry      visc_bndry;
                    ABecLaplacian* visc_op;

                    visc_op = diffusion->getViscOp(sigma,a,b,cur_time,
                                                   visc_bndry,Rh,
                                                   rho_flag,&rhsscale,
                                                   rhoh_visc,dataComp,alpha,0);

                    visc_op->maxOrder(diffusion->maxOrder());

		    // copy rho^{n+1} * (delta Y)^sync into Soln
                    MultiFab::Copy(Soln,Ssync,sigma-BL_SPACEDIM,0,1,0);

		    // divide Soln by rho
                    for (MFIter Smfi(Soln); Smfi.isValid(); ++Smfi)
		    {
		      Soln[Smfi].divide(S_new[Smfi],Smfi.validbox(),Density,0,1);
		    }

		    // compute lambda/cp.Grad(delta Y^sync) and weight
		    visc_op->compFlux(D_DECL(*fluxSC[0],*fluxSC[1],*fluxSC[2]),Soln);
                    for (int d = 0; d < BL_SPACEDIM; ++d)
                        fluxSC[d]->mult(-b/geom.CellSize()[d]);
                    //
                    // Here, get fluxNULN = (lambda/cp - rho.D)Grad(delta Ysync)
                    //                    = lambda/cp.Grad(delta Ysync) + SpecSyncDiffFlux
                    //
                    BL_ASSERT(spec_diffusion_flux_computed[comp]==HT_SyncDiffusion);

                    for (int d = 0; d < BL_SPACEDIM; ++d)
                    {
                        MFIter SDF_mfi(*SpecDiffusionFluxnp1[d]);

                        for ( ; SDF_mfi.isValid(); ++SDF_mfi)
                        {
                            FArrayBox& fluxSC_fab   = (*fluxSC[d])[SDF_mfi];
                            FArrayBox& fluxNULN_fab = (*fluxNULN[d])[SDF_mfi];
                            FArrayBox& SDF_fab = (*SpecDiffusionFluxnp1[d])[SDF_mfi];
                            const Box& ebox    = SDF_mfi.validbox();
			    // copy in (delta Gamma)^sync
                            fluxNULN_fab.copy(SDF_fab,ebox,comp,ebox,comp,1);
			    // add in (lambda/cp) grad (delta Y^sync)
                            fluxNULN_fab.plus(fluxSC_fab,ebox,0,comp,1);
                        }
                    }
                    delete visc_op;
                }

                Soln.clear();
		fb_fluxSC.clear();
		fb_rhoh_visc.clear();
                //
                // Multiply fluxi by h_i (let FLXDIV routine below sum up the fluxes)
                //
                FArrayBox eTemp, h;

                for (FillPatchIterator Tnew_fpi(*this,S_new,nGrow,cur_time,State_Type,Temp,1);
                     Tnew_fpi.isValid();
                     ++Tnew_fpi)
                {
                    const Box& box = Tnew_fpi.validbox();

                    for (int d = 0; d < BL_SPACEDIM; ++d)
                    {
                        const Box& ebox = BoxLib::surroundingNodes(box,d);
                        eTemp.resize(ebox,1);
                        FPLoc bc_lo = fpi_phys_loc(get_desc_lst()[State_Type].getBC(Temp).lo(d));
                        FPLoc bc_hi = fpi_phys_loc(get_desc_lst()[State_Type].getBC(Temp).hi(d));
                        
                        center_to_edge_fancy(Tnew_fpi(),eTemp,BoxLib::grow(box,BoxLib::BASISV(d)),
                                             0,0,1,geom.Domain(),bc_lo,bc_hi);
                        
                        h.resize(ebox,nspecies);
                        getChemSolve().getHGivenT(h,eTemp,ebox,0,0);

			// multiply fluxNULN by h_m
                        (*fluxNULN[d])[Tnew_fpi].mult(h,ebox,0,0,nspecies);
                    }
                }

		// add the NULN fluxes to the RHS of the (delta h)^sync diffusion solve
		// afterwards, the entire RHS should be ready
                for (MFIter Ssync_mfi(Ssync); Ssync_mfi.isValid(); ++Ssync_mfi)
                {
                    FArrayBox& syncn = Ssync[Ssync_mfi];
                    const FArrayBox& synco = Ssync[Ssync_mfi];
                    const Box& box = Ssync_mfi.validbox();

                    //
                    // Multiply by dt*dt, one to make it extensive, and one because
                    // Ssync multiplied above by dt, need same units here.
                    //
                    const Real mult = dt*dt;
                    const int sigmaRhoH = RhoH - BL_SPACEDIM; // RhoH comp in Ssync
		    
                    FORT_INCRWEXTFLXDIV(box.loVect(), box.hiVect(),
                                        (*fluxNULN[0])[Ssync_mfi].dataPtr(),
                                        ARLIM((*fluxNULN[0])[Ssync_mfi].loVect()),
                                        ARLIM((*fluxNULN[0])[Ssync_mfi].hiVect()),
                                        (*fluxNULN[1])[Ssync_mfi].dataPtr(),
                                        ARLIM((*fluxNULN[1])[Ssync_mfi].loVect()),
                                        ARLIM((*fluxNULN[1])[Ssync_mfi].hiVect()),
#if BL_SPACEDIM == 3
                                        (*fluxNULN[2])[Ssync_mfi].dataPtr(),
                                        ARLIM((*fluxNULN[2])[Ssync_mfi].loVect()),
                                        ARLIM((*fluxNULN[2])[Ssync_mfi].hiVect()),
#endif
                                        synco.dataPtr(sigmaRhoH),
                                        ARLIM(synco.loVect()),
                                        ARLIM(synco.hiVect()),
                                        syncn.dataPtr(sigmaRhoH),
                                        ARLIM(syncn.loVect()),
                                        ARLIM(syncn.hiVect()),
                                        volume[Ssync_mfi].dataPtr(),
                                        ARLIM(volume[Ssync_mfi].loVect()),
                                        ARLIM(volume[Ssync_mfi].hiVect()),
                                        &nspecies, &mult);
                }
            }

	    FluxBoxes fb_flux(this);
	    MultiFab **flux = fb_flux.get();

            for (sigma = 0; sigma < numscal; sigma++)
            {
                int rho_flag = 0;
                int do_viscsyncflux = do_reflux;
		const int state_ind = BL_SPACEDIM + sigma;
		//
		// To diffuse, or not?
		// (1) Density, no
		// (2) RhoH...if diffusive
		// (3) Trac...if diffusive
		// (4) Spec:
		//    (a) if Le==1, and spec diffusive
		//    (b) if Le!=1, do differential diffusion instead, done above
		// (5) Temp, no (set instead by RhoH to Temp)
                //
		const bool is_spec = state_ind<=last_spec && state_ind>=first_spec;
                int do_it
		    =  state_ind!=Density 
		    && state_ind!=Temp
		    && is_diffusive[state_ind]
		    && !(is_spec && !unity_Le);
		
		if (do_it && (is_spec || state_ind==RhoH))
		    rho_flag = 2;

                if (do_it)
                {
                    MultiFab* alpha = 0;
                    getDiffusivity(beta, cur_time, state_ind, 0, 1);
                    
		    // on entry, Ssync = RHS for (delta h)^sync diffusive solve
		    // on exit, Ssync = rho^{n+1} * (delta h)^sync
		    // on exit, flux = coeff * grad phi
		    diffusion->diffuse_Ssync(Ssync,sigma,dt,be_cn_theta,Rh,
					     rho_flag,flux,0,beta,0,alpha,0);
		    if (do_viscsyncflux && level > 0)
		    {
			for (MFIter mfi(Ssync); mfi.isValid(); ++mfi)
			{
			    const int i=mfi.index();
			    for (int d=0; d<BL_SPACEDIM; ++d)
                                getViscFluxReg().FineAdd((*flux[d])[mfi],d,i,0,state_ind,1,dt);
			}
		    }
                }
            }
        }
        //
        // For all conservative variables Q (other than density)
        // increment sync by (sync_for_rho)*q_presync.
	// Before this loop, Ssync holds rho^{n+1} (delta phi)^sync
	// DeltaSsync holds (delta rho)^sync phi^p
        //
        for (MFIter mfi(Ssync); mfi.isValid(); ++mfi)
        {
            const int i = mfi.index();

            int iconserved = -1;

            for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
            {
                if (istate != Density && advectionType[istate] == Conservative)
                {
                    iconserved++;

                    Ssync[mfi].plus((*DeltaSsync)[mfi],grids[i],iconserved,istate-BL_SPACEDIM,1);
                }
            }
        }
        sync_cleanup(DeltaSsync);
        //
        // Increment the state (for all but rho, since that was done above)
        //
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            const int i = mfi.index();

            for (int sigma = 0; sigma < numscal; sigma++)
            {
                if (!(BL_SPACEDIM+sigma == Density))
                {
                    S_new[mfi].plus(Ssync[mfi],grids[i],sigma,BL_SPACEDIM+sigma,1);
                }
            }
        }
        //
        // Recompute temperature and rho R T after the mac_sync.
        //
        RhoH_to_Temp(S_new);
        setThermoPress(cur_time);
        //
        // Get boundary conditions.
        //
        Real mult = 1.0;
        Array<int*>         sync_bc(grids.size());
        Array< Array<int> > sync_bc_array(grids.size());
        for (int i = 0; i < ngrids; i++)
        {
            sync_bc_array[i] = getBCArray(State_Type,i,Density,numscal);
            sync_bc[i]       = sync_bc_array[i].dataPtr();
        }
        //
        // Interpolate the sync correction to the finer levels.
        //
        IntVect ratio = IntVect::TheUnitVector();
        for (int lev = level+1; lev <= finest_level; lev++)
        {
            ratio                   *= parent->refRatio(lev-1);
            HeatTransfer& fine_level = getLevel(lev);
            MultiFab& S_new_lev      = fine_level.get_new_data(State_Type);
            //
            // New way of interpolating syncs to make sure mass is conserved
            // and to ensure freestream preservation for species & temperature.
            //
            const BoxArray& fine_grids = S_new_lev.boxArray();
            const int nghost           = S_new_lev.nGrow();
            MultiFab increment(fine_grids, numscal, nghost);
            increment.setVal(0,nghost);
            //
            // Note: we use the lincc_interp (which_interp==3) for density,
            // rho*h and rho*Y, cell_cons_interp for everything else. Doing
            // so is needed for freestream preservation of Y and T.  The setting
            // which_interp=5 calls the lincc_interp, but then follows with a
            // "protection" step to be sure that all the components but the
            // first and the last have their sync adjusted to try to preserve
            // positivity after the sync is applied.  The density sync is then
            // adjusted to be the sum of the species syncs.
            //
            // HACK note: Presently, the species mass syncs are redistributed 
            //            without consequence to the enthalpy sync.  Clearly
            //            the species carry enthalphy, so the enthalph sync should
            //            be adjusted as well.  Note yet sure how to do this correctly.
            //            Punt for now...
            //
            const SyncInterpType which_interp = CellConsProt_T;

            const int nComp = 2+nspecies;

            SyncInterp(Ssync, level, increment, lev, ratio, 
                       Density-BL_SPACEDIM, Density-BL_SPACEDIM, nComp, 1, mult, 
                       sync_bc.dataPtr(), which_interp, Density);

            if (have_trac)
                SyncInterp(Ssync, level, increment, lev, ratio, 
                           Trac-BL_SPACEDIM, Trac-BL_SPACEDIM, 1, 1, mult, 
                           sync_bc.dataPtr());

            if (have_rhort)
                SyncInterp(Ssync, level, increment, lev, ratio, 
                           RhoRT-BL_SPACEDIM, RhoRT-BL_SPACEDIM, 1, 1, mult, 
                           sync_bc.dataPtr());

            SyncInterp(Ssync, level, increment, lev, ratio, 
                       Temp-BL_SPACEDIM, Temp-BL_SPACEDIM, 1, 1, mult, 
                       sync_bc.dataPtr());

            if (do_set_rho_to_species_sum)
            {
                increment.setVal(0,Density-BL_SPACEDIM,1,0);

                for (int istate = first_spec; istate <= last_spec; istate++)
                { 
                    for (MFIter mfi(increment); mfi.isValid(); ++mfi)
                    {
                        int i = mfi.index();
                        increment[mfi].plus(increment[mfi],fine_grids[i],
					  istate-BL_SPACEDIM,Density-BL_SPACEDIM,1);
                    }
                }
            }

            for (MFIter mfi(increment); mfi.isValid(); ++mfi)
            {
                int i = mfi.index();
                S_new_lev[mfi].plus(increment[mfi],fine_grids[i],0,Density,numscal);
            }
            fine_level.make_rho_curr_time();
            fine_level.incrRhoAvg(increment,Density-BL_SPACEDIM,1.0);
            //
            // Recompute temperature and rho R T after interpolation of the mac_sync correction
            //   of the individual quantities rho, Y, T.
            //
            RhoH_to_Temp(S_new_lev);
            fine_level.setThermoPress(cur_time);
        }
        //
        // Average down Trac = rho R T after interpolation of the mac_sync correction
        //   of the individual quantities rho, Y, T.
        //
        for (int lev = finest_level-1; lev >= level; lev--)
        {
            HeatTransfer&   fine_lev = getLevel(lev+1);
            const Geometry& fgeom    = fine_lev.geom;
            
            HeatTransfer&   crse_lev = getLevel(lev);
            const Geometry& cgeom    = crse_lev.geom;
            const IntVect&  fratio   = crse_lev.fine_ratio;
            
            MultiFab& S_crse = crse_lev.get_new_data(State_Type);
            MultiFab& S_fine = fine_lev.get_new_data(State_Type);

            const int pComp = (have_rhort ? RhoRT : Trac);
            BoxLib::average_down(S_fine,S_crse,
				 fgeom, cgeom,
                                 pComp,1,fratio);
        }
    }

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer:mac_sync(): lev: " << level << ", time: " << run_time << '\n';
    }
}

void
HeatTransfer::differential_spec_diffuse_sync (Real dt)
{
    BL_PROFILE("HeatTransfer::differential_spec_diffuse_sync()");
    //
    // Diffuse the species syncs such that sum(SpecDiffSyncFluxes) = 0
    // After exiting, SpecDiffusionFluxnp1 should contain rhoD grad (delta Y)^sync
    // Also, Ssync for species should contain rho^{n+1} * (delta Y)^sync
    //
    if (hack_nospecdiff)
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "... HACK!!! skipping spec sync diffusion " << '\n';

        for (int d=0; d<BL_SPACEDIM; ++d)
            SpecDiffusionFluxnp1[d]->setVal(0.0,0,nspecies);

        for (int comp=0; comp<nspecies; ++comp)
            spec_diffusion_flux_computed[comp] = HT_SyncDiffusion;

        return;            
    }

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Doing differential sync diffusion ..." << '\n';

    if (verbose && benchmarking) ParallelDescriptor::Barrier();

    const Real strt_time = ParallelDescriptor::second();
    //
    // Do implicit c-n solve for each scalar...but dont reflux.
    // Save the fluxes, coeffs and source term, we need 'em for later
    // Actually, since Ssync multiplied by dt in mac_sync before
    // a call to this routine (I hope...), Ssync is in units of s,
    // not the "usual" ds/dt...lets convert it (divide by dt) so
    // we can use a generic flux adjustment function
    //
    const Real cur_time = state[State_Type].curTime();
    FluxBoxes fb_betanp1(this, nspecies, 0);
    MultiFab **betanp1 = fb_betanp1.get();
    getDiffusivity(betanp1, cur_time, first_spec, 0, nspecies); // rhoD

    MultiFab Rhs(grids,nspecies,0);
    const int spec_Ssync_sComp = first_spec - BL_SPACEDIM;

    // Rhs and Ssync contain the RHS of DayBell:2000 Eq (18)
    // with the additional -Y_m^{n+1,p} * (delta rho)^sync term
    // Copy this into Rhs; we will need this later since we overwrite SSync in the solves
    MultiFab::Copy(Rhs,Ssync,spec_Ssync_sComp,0,nspecies,0);
    //
    // Some standard settings
    //
    const Array<int> rho_flag(nspecies,2);
    const MultiFab* alpha = 0;
    FluxBoxes fb_fluxSC(this);
    MultiFab** fluxSC = fb_fluxSC.get();

    const MultiFab& RhoHalftime = get_rho_half_time();

    for (int sigma = 0; sigma < nspecies; ++sigma)
    {
        //
        // Here, we use Ssync as a source in units of s, as expected by diffuse_Ssync
        // (i.e., ds/dt ~ d(Ssync)/dt, vs. ds/dt ~ Rhs in diffuse_scalar).  This was
        // apparently done to mimic diffuse_Vsync, which does the same, because the
        // diffused result is an acceleration, not a velocity, req'd by the projection.
        //
	const int ssync_ind = first_spec + sigma - Density;
                    
	// on entry, Ssync = RHS for (delta Ytilde)^sync diffusive solve
	// on exit, Ssync = rho^{n+1} * (delta Ytilde)^sync
	// on exit, fluxSC = rhoD grad (delta Ytilde)^sync
	diffusion->diffuse_Ssync(Ssync,ssync_ind,dt,be_cn_theta,
				 RhoHalftime,rho_flag[sigma],fluxSC,0,
                                 betanp1,sigma,alpha,0);
	//
	// Pull fluxes into flux array
	// this is the rhoD grad (delta Ytilde)^sync terms in DayBell:2000 Eq (18)
	//
	for (int d=0; d<BL_SPACEDIM; ++d)
	{
	  MultiFab::Copy(*SpecDiffusionFluxnp1[d],*fluxSC[d],0,sigma,1,0);
	}

	spec_diffusion_flux_computed[sigma] = HT_SyncDiffusion;
    }
    //
    // Modify update/fluxes to preserve flux sum = 0
    // (Be sure to pass the "normal" looking Rhs to this generic function)
    //
    const int sCompS = first_spec - BL_SPACEDIM;
    const MultiFab* old_sync = 0;
    const int dataComp = 0; 
    Rhs.mult(1.0/dt,0,nspecies,0); // adjust_spec_diffusion_update needs Rhs in units of dsdt
    //
    // on entry, Ssync = rho^{n+1} * (delta Ytilde)^sync
    // on exit,  Ssync = rho^{n+1} * (delta Y)^sync
    // on entry, SpecDiffusionFluxnp1 = rhoD grad (delta Ytilde)^sync
    // on exit,  SpecDiffusionFluxnp1 = rhoD grad (delta Y)^sync
    // on entry, Rhs = (1/dt) * "RHS of diffusion solve"
    //
    adjust_spec_diffusion_update(Ssync,old_sync,sCompS,dt,cur_time,rho_flag,
                                 RhoHalftime,dataComp,&Rhs,alpha,betanp1);

    //
    // Do refluxing AFTER flux adjustment
    //
    if (do_reflux && level > 0)
    {
      for (int d=0; d<BL_SPACEDIM; ++d)
      {
	for (MFIter fmfi(*SpecDiffusionFluxnp1[d]); fmfi.isValid(); ++fmfi)
	{
	  const int i=fmfi.index();
	  getViscFluxReg().FineAdd((*SpecDiffusionFluxnp1[d])[fmfi],d,i,0,first_spec,nspecies,dt);
	}
      }
    }

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::differential_spec_diffuse_sync(): time: " << run_time << '\n';
    }
}

void
HeatTransfer::reflux ()
{
    if (level == parent->finestLevel()) return;

    BL_ASSERT(do_reflux);
    //
    // First do refluxing step.
    //
    FluxRegister& fr_adv  = getAdvFluxReg(level+1);
    FluxRegister& fr_visc = getViscFluxReg(level+1);
    const Real    dt_crse = parent->dtLevel(level);
    const Real    scale   = 1.0/dt_crse;
    //
    // It is important, for do_mom_diff == 0, to do the viscous
    // refluxing first, since this will be divided by rho_half
    // before the advective refluxing is added.  In the case of
    // do_mom_diff == 1, both components of the refluxing will
    // be divided by rho^(n+1) in NavierStokes::level_sync.
    //
    // take divergence of diffusive flux registers into cell-centered RHS.
    //
    fr_visc.Reflux(Vsync,volume,scale,0,0,BL_SPACEDIM,geom);
    if (do_reflux_visc)
        fr_visc.Reflux(Ssync,volume,scale,BL_SPACEDIM,0,NUM_STATE-BL_SPACEDIM,geom);

    const MultiFab& RhoHalftime = get_rho_half_time();

    if (do_mom_diff == 0) 
    {
        for (MFIter mfi(Vsync); mfi.isValid(); ++mfi)
        {
            const int i = mfi.index();

            D_TERM(Vsync[mfi].divide(RhoHalftime[mfi],grids[i],0,Xvel,1);,
                   Vsync[mfi].divide(RhoHalftime[mfi],grids[i],0,Yvel,1);,
                   Vsync[mfi].divide(RhoHalftime[mfi],grids[i],0,Zvel,1););
        }
    }

    FArrayBox tmp;
    //
    // For any variables that used non-conservative advective differencing,
    // divide the sync by rhohalf.
    //
    for (MFIter mfi(Ssync); mfi.isValid(); ++mfi)
    {
        const int i = mfi.index();

        tmp.resize(grids[i],1);
        tmp.copy(RhoHalftime[mfi],0,0,1);
        tmp.invert(1);

        for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
        {
            if (advectionType[istate] == NonConservative)
            {
                const int sigma = istate -  BL_SPACEDIM;

                Ssync[mfi].mult(tmp,0,sigma,1);
            }
        }
    }
    //
    // Take divergence of advective flux registers into cell-centered RHS.
    //
    fr_adv.Reflux(Vsync,volume,scale,0,0,BL_SPACEDIM,geom);
    fr_adv.Reflux(Ssync,volume,scale,BL_SPACEDIM,0,NUM_STATE-BL_SPACEDIM,geom);

    BoxArray baf = getLevel(level+1).boxArray();

    baf.coarsen(fine_ratio);
    //
    // This is necessary in order to zero out the contribution to any
    // coarse grid cells which underlie fine grid cells.
    //
    std::vector< std::pair<int,Box> > isects;

    for (MFIter mfi(Vsync); mfi.isValid(); ++mfi)
    {
        baf.intersections(grids[mfi.index()],isects);

        for (int i = 0, N = isects.size(); i < N; i++)
        {
            Vsync[mfi].setVal(0,isects[i].second,0,BL_SPACEDIM);
            Ssync[mfi].setVal(0,isects[i].second,0,NUM_STATE-BL_SPACEDIM);
        }
    }
}

void
HeatTransfer::set_preferred_boundary_values (MultiFab& S,
                                             int       state_index,
                                             int       src_comp,
                                             int       dst_comp,
                                             int       num_comp,
                                             Real      time) const
{
    //
    // Only do copy if request contains affected states,
    // and fillpatched data known to be no good.
    //
    if (state_index == State_Type)
    {
        const TimeLevel whichTime = which_time(State_Type,time);

        if (!FillPatchedOldState_ok && whichTime == AmrOldTime)
        {
            //
            // To get chem-advanced data instead of FP'd data at old time.
            //
            if (src_comp > BL_SPACEDIM)
            {
                aux_boundary_data_old.copyTo(S,src_comp-BL_SPACEDIM,dst_comp,num_comp);
            }
        }

        if (!FillPatchedNewState_ok && whichTime == AmrNewTime)
        {
            //
            // To get RhoH computed with current T instead of FP'd RhoH.
            //
            if (src_comp <= RhoH && src_comp + num_comp > RhoH)
            {
                aux_boundary_data_new.copyTo(S,0,RhoH-src_comp+dst_comp,1);
            }
        }
    }
}

void
HeatTransfer::calcViscosity (const Real time,
			     const Real dt,
			     const int  iteration,
			     const int  ncycle)
{
    const TimeLevel whichTime = which_time(State_Type, time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    compute_vel_visc(time, whichTime == AmrOldTime ? viscn_cc : viscnp1_cc);
}

void
HeatTransfer::calcDiffusivity (const Real time)
{
    calcDiffusivity(time,false);
}

void
HeatTransfer::calcDiffusivity (const Real time,
                               bool       do_VelVisc)
{
    BL_PROFILE("HeatTransfer::calcDiffusivity()");

    const TimeLevel whichTime = which_time(State_Type, time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    const int  dotemp         = 1;
    const int  vflag          = do_VelVisc;
    const int  nc_bcen        = nspecies+2;
    const int  nGrow          = 1;
    const int  offset         = BL_SPACEDIM + 1; // No diffusion coeff for vels or rho
    const int  num_comp       = NUM_STATE-offset;
    const int  last_comp      = offset + num_comp - 1;
    const bool has_spec       = offset < last_spec && last_comp > first_spec;
    const int  non_spec_comps = std::max(0,first_spec-offset) + std::max(0,last_comp-last_spec);
    MultiFab&  visc           = (whichTime == AmrOldTime) ? (*diffn_cc) : (*diffnp1_cc);
    MultiFab&  beta           = (whichTime == AmrOldTime) ? (*viscn_cc) : (*viscnp1_cc);

    BL_ASSERT(first_spec == Density+1);
    BL_ASSERT(nspecies > 0 && has_spec);
    BL_ASSERT(num_comp >= nspecies+non_spec_comps);

    MultiFab& S_new = get_new_data(State_Type);

    FArrayBox bcen, temp, rhospec, cpmix;

    Real p_amb;
    FORT_GETPAMB(&p_amb);

    for (FillPatchIterator Rho_and_spec_fpi(*this,S_new,nGrow,time,State_Type,Density,nspecies+1),
             Temp_fpi(*this,S_new,nGrow,time,State_Type,Temp,1);
         Rho_and_spec_fpi.isValid() && Temp_fpi.isValid();
         ++Rho_and_spec_fpi, ++Temp_fpi)
    {
        const int  idx = Rho_and_spec_fpi.index();
        const Box& gbx = BoxLib::grow(grids[idx],nGrow);
        //
        // Convert from RhoY_l to Y_l
        //
        temp.resize(gbx,1);
        bcen.resize(gbx,nc_bcen);
        rhospec.resize(gbx,nspecies+1);

        temp.copy(Rho_and_spec_fpi(),0,0,1);
        temp.invert(1);

        for (int n = 1; n < nspecies+1; n++)
            Rho_and_spec_fpi().mult(temp,0,n,1);

        rhospec.copy(Rho_and_spec_fpi(),0,0,nspecies+1);

        temp.copy(Temp_fpi(),0,0,1);

        FORT_SPECTEMPVISC(gbx.loVect(),gbx.hiVect(),
                          ARLIM(temp.loVect()),ARLIM(temp.hiVect()),
                          temp.dataPtr(),
                          ARLIM(rhospec.loVect()),ARLIM(rhospec.hiVect()),
                          rhospec.dataPtr(1),
                          ARLIM(bcen.loVect()),ARLIM(bcen.hiVect()),bcen.dataPtr(),
                          &nc_bcen, &P1atm_MKS, &dotemp, &vflag, &p_amb);

        visc[Rho_and_spec_fpi].copy(bcen,0,first_spec-offset,nspecies);
        visc[Rho_and_spec_fpi].copy(bcen,nspecies,Temp-offset,1);

        if (do_VelVisc)
            beta[Rho_and_spec_fpi].copy(bcen,nspecies+1,0,1);
        //
        // Now get the rest.
        //
        for (int icomp = offset; icomp <= last_comp; icomp++)
        {
            const bool is_spec = icomp >= first_spec && icomp <= last_spec;

            if (!is_spec)
            {
                if (icomp == RhoH)
                {
                    visc[Rho_and_spec_fpi].copy(visc[Rho_and_spec_fpi],Temp-offset,RhoH-offset,1);
                    cpmix.resize(gbx,1);
                    const int sCompT = 0, sCompY = 1, sCompCp = 0;
                    getChemSolve().getCpmixGivenTY(cpmix,temp,rhospec,gbx,sCompT,sCompY,sCompCp);
                    visc[Rho_and_spec_fpi].divide(cpmix,0,RhoH-offset,1);
                }
                else if (icomp == Trac || icomp == RhoRT)
                {
                    visc.setVal(trac_diff_coef, icomp-offset, 1, nGrow);
                }
            }
        }
    }
}

void
HeatTransfer::getViscosity (MultiFab*  beta[BL_SPACEDIM],
                            const Real time)
{
    const TimeLevel whichTime = which_time(State_Type, time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    MultiFab* visc = (whichTime == AmrOldTime) ? viscn_cc : viscnp1_cc;

    for (MFIter viscMfi(*visc); viscMfi.isValid(); ++viscMfi)
    {
        const int i = viscMfi.index();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            FPLoc bc_lo = fpi_phys_loc(get_desc_lst()[State_Type].getBC(Density).lo(dir));
            FPLoc bc_hi = fpi_phys_loc(get_desc_lst()[State_Type].getBC(Density).hi(dir));

            center_to_edge_fancy((*visc)[viscMfi],(*beta[dir])[viscMfi],
                                 BoxLib::grow(grids[i],BoxLib::BASISV(dir)), 0, 0, 1,
                                 geom.Domain(), bc_lo, bc_hi);
        }
    }
}


void
HeatTransfer::getDiffusivity (MultiFab*  beta[BL_SPACEDIM],
                              const Real time,
                              const int  state_comp,
                              const int  dst_comp,
                              const int  ncomp)
{
    BL_ASSERT(state_comp > Density);

    const TimeLevel whichTime = which_time(State_Type, time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    MultiFab* diff      = (whichTime == AmrOldTime) ? diffn_cc : diffnp1_cc;
    const int offset    = BL_SPACEDIM + 1; // No diffusion coeff for vels or rho
    int       diff_comp = state_comp - offset;

    for (MFIter diffMfi(*diff); diffMfi.isValid(); ++diffMfi)
    {
        const int i = diffMfi.index();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            FPLoc bc_lo = fpi_phys_loc(get_desc_lst()[State_Type].getBC(Density).lo(dir));
            FPLoc bc_hi = fpi_phys_loc(get_desc_lst()[State_Type].getBC(Density).hi(dir));

            center_to_edge_fancy((*diff)[diffMfi],(*beta[dir])[diffMfi],
                                 BoxLib::grow(grids[i],BoxLib::BASISV(dir)), diff_comp, 
                                 dst_comp, ncomp, geom.Domain(), bc_lo, bc_hi);
        }
    }

    if (siegel_test && (state_comp == Temp || state_comp == RhoH))
        beta[1]->setVal(0);

    if (zeroBndryVisc > 0)
        zeroBoundaryVisc(beta,time,state_comp,dst_comp,ncomp);
}

void
HeatTransfer::zeroBoundaryVisc (MultiFab*  beta[BL_SPACEDIM],
                                const Real time,
                                const int  state_comp,
                                const int  dst_comp,
                                const int  ncomp) const
{
    BL_ASSERT(state_comp > Density);

    const int isrz = (int) geom.IsRZ();
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        Box edom = BoxLib::surroundingNodes(geom.Domain(),dir);
        
        for (MFIter mfi(*(beta[dir])); mfi.isValid(); ++mfi)
        {
            FArrayBox& beta_fab = (*(beta[dir]))[mfi];
            const Box& ebox     = BoxLib::surroundingNodes(mfi.validbox(),dir);
            FORT_ZEROVISC(beta_fab.dataPtr(dst_comp),
                          ARLIM(beta_fab.loVect()), ARLIM(beta_fab.hiVect()),
                          ebox.loVect(),  ebox.hiVect(),
                          edom.loVect(),  edom.hiVect(),
                          geom.CellSize(), geom.ProbLo(), phys_bc.vect(),
                          &dir, &isrz, &state_comp, &ncomp);
        }
    }
}

void
HeatTransfer::compute_vel_visc (Real      time,
                                MultiFab* beta)
{
    BL_PROFILE("HeatTransfer::compute_vel_visc()");

    const int nGrow = beta->nGrow();

    BL_ASSERT(nGrow == 1);
    BL_ASSERT(first_spec == Density+1);

    FArrayBox tmp;

    MultiFab dummy(grids,1,0,Fab_noallocate);

    for (FillPatchIterator Temp_fpi(*this,dummy,nGrow,time,State_Type,Temp,1),
             Rho_and_spec_fpi(*this,dummy,nGrow,time,State_Type,Density,nspecies+1);
         Rho_and_spec_fpi.isValid() && Temp_fpi.isValid();
         ++Rho_and_spec_fpi, ++Temp_fpi)
    {
        const Box& box          = Rho_and_spec_fpi().box();
        FArrayBox& temp         = Temp_fpi();
        FArrayBox& rho_and_spec = Rho_and_spec_fpi();
        //
        // Convert from RhoY_l to Y_l
        //
        tmp.resize(box,1);
        tmp.copy(rho_and_spec,0,0,1);
        tmp.invert(1);

        for (int n = 1; n < nspecies+1; ++n)
            rho_and_spec.mult(tmp,0,n,1);

        FORT_VELVISC(box.loVect(),box.hiVect(),
                     ARLIM(temp.loVect()),ARLIM(temp.hiVect()),temp.dataPtr(),
                     ARLIM(rho_and_spec.loVect()),ARLIM(rho_and_spec.hiVect()),
                     rho_and_spec.dataPtr(1),
                     ARLIM(tmp.loVect()),ARLIM(tmp.hiVect()),tmp.dataPtr());

        (*beta)[Rho_and_spec_fpi].copy(tmp,0,0,1);
    }
}

void
HeatTransfer::calc_divu (Real      time,
                         Real      dt,
                         MultiFab& divu)
{
    //
    // Get Mwmix, cpmix and pressure
    //
    const int nGrow = 0;

    int sCompR, sCompT, sCompY, sCompCp, sCompMw;
    sCompR=sCompT=sCompY=sCompCp=sCompMw=0;	

    MultiFab rho(grids,1,nGrow), temp(grids,1,nGrow);

    const MultiFab& Rho_time = get_rho(time);

    for (FillPatchIterator Temp_fpi(*this,temp,nGrow,time,State_Type,Temp,1);
         Temp_fpi.isValid();
         ++Temp_fpi)
    {
        rho[Temp_fpi].copy(Rho_time[Temp_fpi],0,sCompR,1);
        temp[Temp_fpi].copy(Temp_fpi(),0,sCompT,1);
    }
    //
    // Note that state contains rho*species, so divide species by rho.
    //
    showMF("divu",rho,"divu_rho",level);
    showMF("divu",temp,"divu_temp",level);

    MultiFab species(grids,nspecies,nGrow);

    FArrayBox tmp;

    for (FillPatchIterator Spec_fpi(*this,species,nGrow,time,State_Type,first_spec,nspecies);
         Spec_fpi.isValid();
         ++Spec_fpi)
    {
        const int i = Spec_fpi.index();

        species[Spec_fpi].copy(Spec_fpi(),0,sCompY,nspecies);

        tmp.resize(grids[i],1);
        tmp.copy(rho[Spec_fpi],sCompR,0,1);
        tmp.invert(1);

        for (int ispecies = 0; ispecies < nspecies; ispecies++)
            species[Spec_fpi].mult(tmp,0,ispecies,1);
    }

    showMF("divu",species,"divu_species",level);

    MultiFab mwmix(grids,1,nGrow), cp(grids,1,nGrow);

    for (MFIter Rho_mfi(rho); Rho_mfi.isValid(); ++Rho_mfi)
    {
        const int  iGrid = Rho_mfi.index();
        const Box& box   = Rho_mfi.validbox();

        BL_ASSERT(box == grids[iGrid]);

        getChemSolve().getMwmixGivenY(mwmix[Rho_mfi],species[Rho_mfi],
                                      box,sCompY,sCompMw);
        getChemSolve().getCpmixGivenTY(cp[Rho_mfi],temp[Rho_mfi],species[Rho_mfi],
                                       box,sCompT,sCompY,sCompCp);
    }

    showMF("divu",cp,"divu_cp",level);
    showMF("divu",mwmix,"divu_mwmix",level);

    species.clear();
    //
    // divu = 1/T DT/dt + W_mix * sum (1/W)DY/Dt
    //
    // Compute rho*DT/Dt as
    //
    //   1/c_p (div lambda grad T + sum_l rho D grad h_l dot grad Y_l)
    //
    {
        MultiFab visc_terms(grids,1,1);
        getViscTerms(visc_terms,Temp,1,time);
        MultiFab::Copy(divu,visc_terms,0,0,1,nGrow);
    }
    for (MFIter Divu_mfi(divu); Divu_mfi.isValid(); ++Divu_mfi)
    {
        const int iGrid = Divu_mfi.index();
        divu[Divu_mfi].divide(rho[Divu_mfi],grids[iGrid],sCompR,0,1);
        divu[Divu_mfi].divide(temp[Divu_mfi],grids[iGrid],sCompT,0,1);
    }
    showMF("divu",divu,"divu_VT_T_over_rhoT",level);

    MultiFab delta_divu(grids,1,nGrow), spec_visc_terms(grids,nspecies,0);

    delta_divu.setVal(0.0);

    const Array<Real>& mwt = getChemSolve().speciesMolecWt();

    getViscTerms(spec_visc_terms,first_spec,nspecies,time);

    for (MFIter mfi(spec_visc_terms); mfi.isValid(); ++mfi)
    {
        const int iGrid = mfi.index();

        for (int comp = 0; comp < nspecies; ++comp)
        {
            spec_visc_terms[mfi].mult(1.0/mwt[comp],comp,1);
            delta_divu[mfi].plus(spec_visc_terms[mfi],grids[iGrid],comp,0,1);
        }
    }
    showMF("divu",delta_divu,"divu_sum_VT_Y_over_Wi",level);

    spec_visc_terms.clear();

    for (MFIter Divu_mfi(divu); Divu_mfi.isValid(); ++Divu_mfi)
    {
        const Box& box   = Divu_mfi.validbox();
        delta_divu[Divu_mfi].divide(rho[Divu_mfi],box,sCompR,0,1);
        delta_divu[Divu_mfi].mult(mwmix[Divu_mfi],box,0,0,1);
        divu[Divu_mfi].plus(delta_divu[Divu_mfi],box,0,0,1);
    }
    showMF("divu",divu,"divu_1",level);

    rho.clear();

    if (dt > 0.0)
    {
        //
        // Increment divu by
        //    sum_l (h_l/(c_p*T) - mw_mix/mw_l)*delta Y_l/dt 
        // (i.e., Y_l"dot")
        //
        FArrayBox h;

        const int sCompH = 0;

        for (FillPatchIterator Ydot_fpi(*this,delta_divu,0,time,Ydot_Type,0,nspecies);
             Ydot_fpi.isValid();
             ++Ydot_fpi)
        {
            const int i = Ydot_fpi.index();

            h.resize(BoxLib::grow(grids[i],nGrow),nspecies);

            getChemSolve().getHGivenT(h,temp[Ydot_fpi],grids[i],sCompT,sCompH);

            for (int istate = first_spec; istate <= last_spec; istate++)
            {
                const int ispec = istate-first_spec;

                delta_divu[Ydot_fpi].copy(h,ispec,0,1);
                delta_divu[Ydot_fpi].divide(cp[Ydot_fpi]);
                delta_divu[Ydot_fpi].divide(temp[Ydot_fpi]);
                delta_divu[Ydot_fpi].mult(Ydot_fpi(),ispec,0,1);
                divu[Ydot_fpi].plus(delta_divu[Ydot_fpi]);

                delta_divu[Ydot_fpi].copy(mwmix[Ydot_fpi],0,0,1);
                delta_divu[Ydot_fpi].divide(mwt[ispec]);
                delta_divu[Ydot_fpi].mult(Ydot_fpi(),ispec,0,1);
                divu[Ydot_fpi].minus(delta_divu[Ydot_fpi]);
            }
        }
    }

    showMF("divu",divu,"divu_2",level);
}

//
// Compute the Eulerian Dp/Dt for use in pressure relaxation.
//
void
HeatTransfer::calc_dpdt (Real      time,
                         Real      dt_,
                         MultiFab& dpdt,
                         MultiFab* u_mac)
{
    Real dt = crse_dt, p_amb, dpdt_factor;
    
    FORT_GETPAMB(&p_amb);
    FORT_GETDPDT(&dpdt_factor);
    
    if (dt <= 0.0 || dpdt_factor <= 0)
    {
        dpdt.setVal(0);
        return;
    }

    MultiFab& S_old = get_old_data(State_Type);
    
    const int pComp = (have_rhort ? RhoRT : Trac);
    const int nGrow = 1;
    
    int sCompR, sCompT, sCompY, sCompCp, sCompMw;
    sCompR=sCompT=sCompY=sCompCp=sCompMw=0;
    
    MultiFab temp(grids,1,nGrow);
    MultiFab rhoRT(grids,1,nGrow);
    
    const MultiFab& rho = get_rho(time);

    BL_ASSERT(rho.boxArray()  == S_old.boxArray());
    BL_ASSERT(temp.boxArray() == S_old.boxArray());
    BL_ASSERT(dpdt.boxArray() == S_old.boxArray());
    
    for (FillPatchIterator Temp_fpi(*this,temp,nGrow,time,State_Type,Temp,1);
         Temp_fpi.isValid();
         ++Temp_fpi)
    {
        temp[Temp_fpi].copy(Temp_fpi(),0,sCompT,1);
    }

    for (FillPatchIterator rhoRT_fpi(*this,rhoRT,nGrow,time,State_Type,pComp,1);
         rhoRT_fpi.isValid();
         ++rhoRT_fpi)
    {
        rhoRT[rhoRT_fpi].copy(rhoRT_fpi(),0,0,1);
    }
    //
    // Note that state contains rho*species, so divide species by rho.
    //
    MultiFab species(grids,nspecies,nGrow);

    FArrayBox tmp;

    for (FillPatchIterator Spec_fpi(*this,species,nGrow,time,State_Type,first_spec,nspecies);
         Spec_fpi.isValid();
         ++Spec_fpi)
    {
        const int  i  = Spec_fpi.index();
        const Box& bx = BoxLib::grow(grids[i],nGrow);

        species[Spec_fpi].copy(Spec_fpi(),0,sCompY,nspecies);

        tmp.resize(bx,1);
        tmp.copy(rho[Spec_fpi],sCompR,0,1);
        tmp.invert(1);

        for (int ispecies = 0; ispecies < nspecies; ispecies++)
            species[Spec_fpi].mult(tmp,0,ispecies,1);
    }

    tmp.clear();
    
    MultiFab mwmix(grids,1,nGrow), cp(grids,1,nGrow);
    
    for (MFIter Rho_mfi(rho); Rho_mfi.isValid(); ++Rho_mfi)
    {
        const int  idx = Rho_mfi.index();
        const Box& box = BoxLib::grow(Rho_mfi.validbox(),nGrow);
        
        BL_ASSERT(Rho_mfi.validbox() == grids[idx]);
        
        getChemSolve().getMwmixGivenY(mwmix[Rho_mfi],species[Rho_mfi],box,sCompY,sCompMw);
        getChemSolve().getCpmixGivenTY(cp[Rho_mfi],temp[Rho_mfi],species[Rho_mfi],box,sCompT,sCompY,sCompCp);
    }

    temp.clear();
    species.clear();
    //
    // Now do pressure relaxation.
    //
    const Real* dx = geom.CellSize();
    //
    //  Get gamma = c_p/(c_p-R)
    //
    MultiFab gamma(grids,1,nGrow);

    FArrayBox rmix;

    for (MFIter mfi(gamma); mfi.isValid(); ++mfi)
    {
        const int idx = mfi.index();
        gamma[mfi].copy(cp[mfi]);
        rmix.resize(BoxLib::grow(grids[idx],nGrow),1);
        rmix.setVal(rgas);
        rmix.divide(mwmix[mfi]);
        gamma[mfi].minus(rmix);
        gamma[mfi].divide(cp[mfi]);
        gamma[mfi].invert(1.0);
    }

    rmix.clear();
    cp.clear();
    mwmix.clear();

    FArrayBox ugradp, p_denom;
    
    for (MFIter mfi(dpdt); mfi.isValid(); ++mfi)
    {
        const int i = mfi.index();
        
        dpdt[mfi].copy(S_old[mfi],grids[i],pComp,grids[i],0,1);
        dpdt[mfi].plus(-p_amb);
        dpdt[mfi].mult(1.0/dt);

        ugradp.resize(grids[i],1);
        
        const int* lo = grids[i].loVect();
        const int* hi = grids[i].hiVect();

        FORT_COMPUTE_UGRADP(rhoRT[mfi].dataPtr(0), 
                            ARLIM(rhoRT[mfi].box().loVect()), 
                            ARLIM(rhoRT[mfi].box().hiVect()),
                            ugradp.dataPtr(), 
                            ARLIM(ugradp.box().loVect()), 
                            ARLIM(ugradp.box().hiVect()),
                            u_mac[0][mfi].dataPtr(),
                            ARLIM(u_mac[0][mfi].box().loVect()),
                            ARLIM(u_mac[0][mfi].box().hiVect()),
                            u_mac[1][mfi].dataPtr(),
                            ARLIM(u_mac[1][mfi].box().loVect()),
                            ARLIM(u_mac[1][mfi].box().hiVect()),
#if (BL_SPACEDIM == 3)
                            u_mac[2][mfi].dataPtr(),
                            ARLIM(u_mac[2][mfi].box().loVect()),
                            ARLIM(u_mac[2][mfi].box().hiVect()),
#endif
                            lo,hi,dx);
        //
        // Note that this is minus because the term we want to add to S is
        //  - Dp/Dt.   We already have (p-p0)/dt = -dp/dt, so now
        // we subtract ugradp to form -Dp/dt.  
        // (Note the dp/dt term is negative because
        //  p is the current value, p0 is the value we're trying to get to,
        //  so dp/dt = (p0 - p)/dt.)
        //
        dpdt[mfi].minus(ugradp,0,0,1);
        //
        // Make sure to divide by gamma *after* subtracting ugradp.
        //
        dpdt[mfi].divide(gamma[mfi]);

        if (dpdt_option == 0)
        {
            dpdt[mfi].divide(S_old[mfi],pComp,0,1);
        }
        else if (dpdt_option == 1)
        {
            dpdt[mfi].divide(p_amb);
        }
        else
        {
            const int ncomp = 1;
            p_denom.resize(grids[i],ncomp);
            p_denom.copy(S_old[mfi],grids[i],pComp,grids[i],0,ncomp);
            Real num_norm = dpdt[mfi].norm(0,0,1);
            FabMinMax(p_denom,grids[i],2*num_norm*Real_MIN,p_amb,0,ncomp);
            dpdt[mfi].divide(p_denom);
        }
        
        dpdt[mfi].mult(dpdt_factor);
    }
}

//
// Function to use if Divu_Type and Dsdt_Type are in the state.
//

void
HeatTransfer::calc_dsdt (Real      time,
                         Real      dt,
                         MultiFab& dsdt)
{
    MultiFab& Divu_new = get_new_data(Divu_Type);
    MultiFab& Divu_old = get_old_data(Divu_Type);

    dsdt.setVal(0);

    for (MFIter mfi(dsdt); mfi.isValid(); ++mfi)
    {
        const int k = mfi.index();
        dsdt[mfi].copy(Divu_new[mfi],grids[k]);
        dsdt[mfi].minus(Divu_old[mfi],grids[k],0,0,1);
        dsdt[mfi].divide(dt,grids[k],0,1);
    }
}

void
HeatTransfer::RhoH_to_Temp (MultiFab& S,
                            int       nGrow)

{
    //
    // Solve hmix = sum(hl(temp).Yl) for the temp field in S.  S is state-like.
    //
    // Be careful -- this is called for levels other than the current level.
    //
    // nGrow is number of grow cells wanted, assumes the rest of the state
    // trusted out there...must make temp same size as S.
    //
    const BoxArray& sgrids = S.boxArray();

    MultiFab temp(sgrids,1,S.nGrow());

    const int do_minmax = 1;

    RhoH_to_Temp(S,temp,nGrow,do_minmax);

    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        const int k   = mfi.index();
        Box       box = BoxLib::grow(sgrids[k],nGrow);
        S[mfi].copy(temp[mfi],box,0,box,Temp,1);
    }
}

//
// Taking an initial guess from S, solve hmix = sum(hl(temp).Yl) for temp.
// hmix and Yl are taken from S (and are assumed to be multiplied by rho
// (i.e. S holds rho.hmix and rho.Yl).
//
// Be careful -- this is called for levels other than the current level
// Note: we only do this on the valid region. The ghost cells
//       are "operated" on, but only in copying S to temp
//       and applying min/max to temp; no actual rhoh-to-temp conversion
//       is applied to the ghost cells...superceeded, see below
// Note: assumes that Temp slot of S contains a reasonable initial guess
//       for temp
// Note: S is unchanged by this function, but we use it for temporary
//       space, in particular, to hold primitives h and Y_l
// Note: nGrow is number of grow cells wanted, assumes the rest of the state
//       trusted out there, nGrow must be same in temp and S
// Note: no good reason for above restriction on nGrow, removed
//
// JFG: RhoH_to_Temp is the top of a rather long sequence of wrappers:
//
//         RhoH_to_Temp    in   HeatTransfer.cpp  calls
//         getTGivenHY     in     ChemDriver.cpp  calls
//         FORT_TfromHY    in  ChemDriver_?D.F    calls
//         FORT_TfromHYpt  in   ChemDriver_F.F    does the work
//

void
HeatTransfer::RhoH_to_Temp (MultiFab& S,
                            MultiFab& temp,
                            int       nGrow,
                            int       dominmax)
{
    BL_ASSERT(S.nGrow() >= nGrow  &&  temp.nGrow() >= nGrow);
    //
    //  If this hasn't been set yet, we cannnot do it correct here (scan multilevel),
    //   just wing it.
    //
    const Real htt_hmixTYP_SAVE = htt_hmixTYP; 
    if (htt_hmixTYP <= 0)
    {
        if (typical_values[RhoH]==typical_RhoH_value_default)
	{
            htt_hmixTYP = S.norm0(RhoH);
	}
        else
	{
            htt_hmixTYP = typical_values[RhoH];
	}        
        if (ParallelDescriptor::IOProcessor())
            std::cout << "setting htt_hmixTYP = " << htt_hmixTYP << '\n';
    }

    const BoxArray& sgrids = S.boxArray();

    const int sCompT    = 0;
    int       max_iters = 0;

    MultiFab::Copy(temp,S,Temp,sCompT,1,nGrow);

    FArrayBox tmp;

    const Real eps = getChemSolve().getHtoTerrMAX();
    Real errMAX = eps*htt_hmixTYP;

    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        const int k   = mfi.index();
        Box       box = BoxLib::grow(sgrids[k],nGrow);
        //
        // Convert rho*h to h and rho*Y to Y for this operation.
        //
        tmp.resize(box,1);
        tmp.copy(S[mfi],Density,0,1);
        tmp.invert(1);

        S[mfi].mult(tmp,0,RhoH,1);

        for (int spec = first_spec; spec <= last_spec; spec++)
            S[mfi].mult(tmp,0,spec,1);

        int iters = getChemSolve().getTGivenHY(temp[mfi],S[mfi],S[mfi],
                                               box,RhoH,first_spec,sCompT,errMAX);
        if (iters < 0)
            BoxLib::Error("HeatTransfer::RhoH_to_Temp: error in H->T");

        max_iters = std::max(max_iters,iters);

        if (dominmax)
            FabMinMax(S[mfi], box, htt_tempmin, htt_tempmax, Temp, 1);
        //
        // Convert back to rho*h and rho*Y
        //
        S[mfi].mult(S[mfi],box,Density,RhoH,1);
        for (int spec = first_spec; spec <= last_spec; spec++)
            S[mfi].mult(S[mfi],box,Density,spec,1);
    }

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceIntMax(max_iters,IOProc);

        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::RhoH_to_Temp: max_iters = " << max_iters << '\n';
    }
    // Reset it back
    htt_hmixTYP = htt_hmixTYP_SAVE;
}

#if 0
void
HeatTransfer::compute_cp_and_hmix (const MultiFab& S,
                                   MultiFab&       cp, 
                                   MultiFab&       hmix,
                                   MultiFab*       temp, 
                                   int             nGrow,
                                   int             calchmix,
                                   int             floor_spec)
{
    //
    // Notes:
    //  1) S has the same number of components as the state.
    //     However, he RhoH-th and the first-spec_th through last_spec-th
    //     components are assumed to have been converted to primitive
    //     form, i.e., they hold h and Y, not rho*h and rho*Y
    //  2) We assume that S, cp, and hmix all have the same set of boxes
    //     and the same number of ghost cells
    //  3) We only compute on the valid region regardless of the number of 
    //     ghost cells.
    //
    BL_ASSERT(S.nGrow() >= nGrow && cp.nGrow() >= nGrow);
    BL_ASSERT(S.nGrow() == cp.nGrow());

    const BoxArray& sgrids   = S.boxArray();
    bool            tmp_temp = temp != 0;
    int             nCompY   = last_spec - first_spec + 1;
    bool            tmp_spec = floor_spec && nCompY>0;
    if (tmp_spec)
	BoxLib::Error("Spec flooring not currently implemented");
    
    bool need_tmp_data = tmp_temp || tmp_spec;

    FArrayBox tmp;

    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
	const int iGrid = mfi.index();
        Box       box   = BoxLib::grow(sgrids[iGrid],nGrow);

	int sCompY = -1, sCompT = -1;

	if (need_tmp_data)
	{
	    int nComp = 1 + nCompY;
	    tmp.resize(box,nComp);
	    sCompT = 0;
	    sCompY = 1;
	    if (tmp_temp)
	    {
		tmp.copy((*temp)[mfi],0,sCompT,1);
	    }
            else
            {
		tmp.copy(S[mfi],Temp,sCompT,1);
	    }
	    tmp.copy(S[mfi],first_spec,sCompY,nCompY);
	}
        else
        {
	    sCompT = Temp;
	    sCompY = first_spec;
	}
	
	const FArrayBox& state = need_tmp_data ? tmp : S[mfi];

	const int sCompCp = 0;
	getChemSolve().getCpmixGivenTY(cp[mfi],state,state,box,sCompT,sCompY,sCompCp);

	if (calchmix)
	{
	    const int sCompH = 0;
	    getChemSolve().getHmixGivenTY(hmix[mfi],state,state,box,sCompT,sCompY,sCompH);
	}
    }
}
#endif

void
HeatTransfer::compute_cp (Real      time,
                          MultiFab& cp)
{
    const int nGrow   = cp.nGrow();
    const int sComp   = std::min(std::min((int)Density,(int)Temp),first_spec);
    const int eComp   = std::max(std::max((int)Density,(int)Temp),first_spec+nspecies-1);
    const int nComp   = eComp - sComp + 1;
    const int sCompR  = Density - sComp;
    const int sCompT  = Temp - sComp;
    const int sCompY  = first_spec - sComp;
    const int sCompCp = 0;

    FArrayBox tmp;

    for (FillPatchIterator fpi(*this,cp,nGrow,time,State_Type,sComp,nComp);
         fpi.isValid();
         ++fpi)
    {
        //
        // Convert rho*Y to Y for this operation
        //
	const int  iGrid    = fpi.index();
	const Box& validBox = BoxLib::grow(grids[iGrid],nGrow);
	FArrayBox& state    = fpi();

        tmp.resize(validBox,1);
        tmp.copy(state,sCompR,0,1);
        tmp.invert(1);

        for (int k = 0; k < nspecies; k++)
            state.mult(tmp,0,sCompY+k,1);
	
	getChemSolve().getCpmixGivenTY(cp[fpi],state,state,validBox,sCompT,sCompY,sCompCp);
    }
}

void
HeatTransfer::compute_cp (const FArrayBox& temp, 
                          const FArrayBox& species,
                          FArrayBox&       cp)
{
    BL_PROFILE("HeatTransfer::compute_cp()");

    const Box& box    = temp.box();
    const int nSpec   = last_spec - first_spec + 1;
    const int nComp   = nSpec + 1;
    const int sCompT  = 0;
    const int sCompY  = 1;
    const int sCompCp = 0;

    cp.resize(box,1);

    FArrayBox state(box,nComp);
    state.copy(temp,0,sCompT,1);
    state.copy(species,0,sCompY,nSpec);

    getChemSolve().getCpmixGivenTY(cp,state,state,box,sCompT,sCompY,sCompCp);
}

void
HeatTransfer::compute_rhohmix (Real      time,
                               MultiFab& rhohmix)
{
    BL_PROFILE("HeatTransfer::compute_rhohmix()");

    const int ngrow  = 0; // We only do this on the valid region
    const int sComp  = std::min(std::min((int)Density,(int)Temp),first_spec);
    const int eComp  = std::max(std::max((int)Density,(int)Temp),first_spec+nspecies-1);
    const int nComp  = eComp - sComp + 1;
    const int sCompR = Density - sComp;
    const int sCompT = Temp - sComp;
    const int sCompY = first_spec - sComp;
    const int sCompH = 0;

    FArrayBox tmp;

    for (FillPatchIterator fpi(*this,rhohmix,ngrow,time,State_Type,sComp,nComp);
         fpi.isValid();
         ++fpi)
    {
        //
        // Convert rho*Y to Y for this operation
        //
	const int  iGrid    = fpi.index();
	const Box& validBox = grids[iGrid];
	FArrayBox& state    = fpi();

        tmp.resize(validBox,1);
        tmp.copy(state,sCompR,0,1);
        tmp.invert(1);

        for (int k = 0; k < nspecies; k++)
            state.mult(tmp,0,sCompY+k,1);
	
	getChemSolve().getHmixGivenTY(rhohmix[fpi],state,state,validBox,
				      sCompT,sCompY,sCompH);
        //
        // Convert hmix to rho*hmix
        //
        rhohmix[fpi].mult(state,validBox,sCompR,sCompH,1);
    }
}
 
void
HeatTransfer::compute_h (Real      time,
                         MultiFab& h)
{
    BL_PROFILE("HeatTransfer::compute_h()");

    BL_ASSERT(h.nComp() == nspecies);
    BL_ASSERT(h.boxArray() == grids);

    const int nGrow = h.nGrow();
    const int sComp = Temp;
    const int nComp = 1;
    const int sCompT = 0;
    const int sCompH = 0;

    for (FillPatchIterator fpi(*this,h,nGrow,time,State_Type,sComp,nComp);
         fpi.isValid();
         ++fpi)
    {
	const Box& box = BoxLib::grow(grids[fpi.index()],nGrow);

	BL_ASSERT(box == h[fpi].box());

	getChemSolve().getHGivenT(h[fpi],fpi(),box,sCompT,sCompH);
    }
}

void
HeatTransfer::setPlotVariables ()
{
    AmrLevel::setPlotVariables();

    const Array<std::string>& names = getChemSolve().speciesNames();

    ParmParse pp("ht");

    bool plot_ydot,plot_rhoY,plot_massFrac,plot_moleFrac,plot_conc;
    plot_ydot=plot_rhoY=plot_massFrac=plot_moleFrac=plot_conc = false;

    if (pp.query("plot_massfrac",plot_massFrac))
    {
        if (plot_massFrac)
        {
            for (int i = 0; i < names.size(); i++)
            {
                const std::string name = "Y("+names[i]+")";
                parent->addDerivePlotVar(name);
            }
        }
        else
        {
            for (int i = 0; i < names.size(); i++)
            {
                const std::string name = "Y("+names[i]+")";
                parent->deleteDerivePlotVar(name);
            }
        }
    }

    if (pp.query("plot_molefrac",plot_moleFrac))
    {
        if (plot_moleFrac)
            parent->addDerivePlotVar("molefrac");
        else
            parent->deleteDerivePlotVar("molefrac");
    }

    if (pp.query("plot_concentration",plot_conc))
    {
        if (plot_conc)
            parent->addDerivePlotVar("concentration");
        else
            parent->deleteDerivePlotVar("concentration");
    }

    if (pp.query("plot_ydot",plot_ydot))
    {
        if (plot_ydot)
        {
            for (int i = 0; i < names.size(); i++)
            {
                const std::string name = "d[Y("+names[i]+")]/dt";
                parent->addStatePlotVar(name);
            }
        }
        else
        {
            for (int i = 0; i < names.size(); i++)
            {
                const std::string name = "d[Y("+names[i]+")]/dt";
                parent->deleteStatePlotVar(name);
            }
        }
    }
  
    if (pp.query("plot_rhoY",plot_rhoY))
    {
        if (plot_rhoY)
        {
            for (int i = 0; i < names.size(); i++)
            {
                const std::string name = "rho.Y("+names[i]+")";
                parent->addStatePlotVar(name);
            }
        }
        else
        {
            for (int i = 0; i < names.size(); i++)
            {
                const std::string name = "rho.Y("+names[i]+")";
                parent->deleteStatePlotVar(name);
            }
        }
    }

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "\nState Plot Vars: ";

        std::list<std::string>::const_iterator li = parent->statePlotVars().begin(), end = parent->statePlotVars().end();

        for ( ; li != end; ++li)
            std::cout << *li << ' ';
        std::cout << '\n';

        std::cout << "\nDerive Plot Vars: ";

        li  = parent->derivePlotVars().begin();
        end = parent->derivePlotVars().end();

        for ( ; li != end; ++li)
            std::cout << *li << ' ';
        std::cout << '\n';
    }
}

void
HeatTransfer::writePlotFile (const std::string& dir,
                             std::ostream&  os,
                             VisMF::How     how)
{
    if ( ! Amr::Plot_Files_Output() ) return;
    //
    // Note that this is really the same as its NavierStokes counterpart,
    // but in order to add diagnostic MultiFabs into the plotfile, code had
    // to be interspersed within this function.
    //
    int i, n;
    //
    // The list of indices of State to write to plotfile.
    // first component of pair is state_type,
    // second component of pair is component # within the state_type
    //
    std::vector<std::pair<int,int> > plot_var_map;
    for (int typ = 0; typ < desc_lst.size(); typ++)
        for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
            if (parent->isStatePlotVar(desc_lst[typ].name(comp)) &&
                desc_lst[typ].getType() == IndexType::TheCellType())
                plot_var_map.push_back(std::pair<int,int>(typ,comp));

    int num_derive = 0;
    std::list<std::string> derive_names;
    const std::list<DeriveRec>& dlist = derive_lst.dlist();

    for (std::list<DeriveRec>::const_iterator it = dlist.begin(), end = dlist.end();
         it != end;
         ++it)
    {
        if (parent->isDerivePlotVar(it->name()))
	{
            derive_names.push_back(it->name());
            num_derive += it->numDerive();
	}
    }

    int num_auxDiag = 0;
    for (std::map<std::string,MultiFab*>::const_iterator it = auxDiag.begin(), end = auxDiag.end();
         it != end;
         ++it)
    {
        num_auxDiag += it->second->nComp();
    }

    int n_data_items = plot_var_map.size() + num_derive + num_auxDiag;
    Real cur_time = state[State_Type].curTime();

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_data_items == 0)
            BoxLib::Error("Must specify at least one valid data item to plot");

        os << n_data_items << '\n';

	//
	// Names of variables -- first state, then derived
	//
	for (i =0; i < plot_var_map.size(); i++)
        {
	    int typ  = plot_var_map[i].first;
	    int comp = plot_var_map[i].second;
	    os << desc_lst[typ].name(comp) << '\n';
        }

	for (std::list<std::string>::const_iterator it = derive_names.begin(), end = derive_names.end();
             it != end;
             ++it)
        {
	    const DeriveRec* rec = derive_lst.get(*it);
	    for (i = 0; i < rec->numDerive(); i++)
                os << rec->variableName(i) << '\n';
        }
        //
        // Hack in additional diagnostics.
        //
        for (std::map<std::string,Array<std::string> >::const_iterator it = auxDiag_names.begin(), end = auxDiag_names.end();
             it != end;
             ++it)
        {
            for (i=0; i<it->second.size(); ++i)
                os << it->second[i] << '\n';
        }

        os << BL_SPACEDIM << '\n';
        os << parent->cumTime() << '\n';
        int f_lev = parent->finestLevel();
        os << f_lev << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbLo(i) << ' ';
        os << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbHi(i) << ' ';
        os << '\n';
        for (i = 0; i < f_lev; i++)
            os << parent->refRatio(i)[0] << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->Geom(i).Domain() << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->levelSteps(i) << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
        {
            for (int k = 0; k < BL_SPACEDIM; k++)
                os << parent->Geom(i).CellSize()[k] << ' ';
            os << '\n';
        }
        os << (int) Geometry::Coord() << '\n';
        os << "0\n"; // Write bndry data.


        // job_info file with details about the run
	std::ofstream jobInfoFile;
	std::string FullPathJobInfoFile = dir;
	std::string PrettyLine = "===============================================================================\n";

	FullPathJobInfoFile += "/job_info";
	jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

	// job information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Job Information\n";
	jobInfoFile << PrettyLine;
	
	jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
	jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif
	jobInfoFile << "\n\n";

        // plotfile information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Plotfile Information\n";
	jobInfoFile << PrettyLine;

	time_t now = time(0);

	// Convert now to tm struct for local timezone
	tm* localtm = localtime(&now);
	jobInfoFile   << "output data / time: " << asctime(localtm);

	char currentDir[FILENAME_MAX];
	if (getcwd(currentDir, FILENAME_MAX)) {
	  jobInfoFile << "output dir:         " << currentDir << "\n";
	}

	jobInfoFile << "\n\n";


        // build information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Build Information\n";
	jobInfoFile << PrettyLine;

	jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
	jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
	jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
	jobInfoFile << "BoxLib dir:    " << buildInfoGetBoxlibDir() << "\n";

	jobInfoFile << "\n";
	
	jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
	jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";
	jobInfoFile << "FCOMP:         " << buildInfoGetFcomp() << "\n";
	jobInfoFile << "FCOMP version: " << buildInfoGetFcompVersion() << "\n";

	jobInfoFile << "\n";

	for (int n = 1; n <= buildInfoGetNumModules(); n++) {
	    jobInfoFile << buildInfoGetModuleName(n) << ": " << buildInfoGetModuleVal(n) << "\n";
	}

	jobInfoFile << "\n";

	const char* githash1 = buildInfoGetGitHash(1);
	const char* githash2 = buildInfoGetGitHash(2);
	const char* githash3 = buildInfoGetGitHash(3);
	if (strlen(githash1) > 0) {
	  jobInfoFile << "LMC    git hash: " << githash1 << "\n";
	}
	if (strlen(githash2) > 0) {
	  jobInfoFile << "BoxLib git hash: " << githash2 << "\n";
	}
	if (strlen(githash3) > 0) {
	  jobInfoFile << "IAMR   git hash: " << githash3 << "\n";
	}

	jobInfoFile << "\n\n";


	// runtime parameters
	jobInfoFile << PrettyLine;
	jobInfoFile << " Inputs File Parameters\n";
	jobInfoFile << PrettyLine;
	
	ParmParse::dumpTable(jobInfoFile, true);

	jobInfoFile.close();
	

    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";

    std::string Level = BoxLib::Concatenate("Level_", level, 1);
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
        FullPath += '/';
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
            BoxLib::CreateDirectoryFailed(FullPath);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (i = 0; i < grids.size(); ++i)
        {
            RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
            for (n = 0; n < BL_SPACEDIM; n++)
                os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        if (n_data_items > 0)
        {
            std::string PathNameInHeader = Level;
            PathNameInHeader += BaseName;
            os << PathNameInHeader << '\n';
        }
    }
    //
    // We combine all of the multifabs -- state, derived, etc -- into one
    // multifab -- plotMF.
    // NOTE: we are assuming that each state variable has one component,
    // but a derived variable is allowed to have multiple components.
    int       cnt   = 0;
    int       ncomp = 1;
    const int nGrow = 0;
    MultiFab  plotMF(grids,n_data_items,nGrow);
    MultiFab* this_dat = 0;
    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (i = 0; i < plot_var_map.size(); i++)
    {
	int typ  = plot_var_map[i].first;
	int comp = plot_var_map[i].second;
	this_dat = &state[typ].newData();
	MultiFab::Copy(plotMF,*this_dat,comp,cnt,ncomp,nGrow);
	cnt+= ncomp;
    }
    //
    // Cull data from derived variables.
    // 
    Real plot_time;

    if (derive_names.size() > 0)
    {
	for (std::list<std::string>::const_iterator it = derive_names.begin(), end = derive_names.end();
             it != end;
             ++it) 
	{
            if (*it == "avg_pressure" ||
                *it == "gradpx"       ||
                *it == "gradpy"       ||
                *it == "gradpz") 
            {
                if (state[Press_Type].descriptor()->timeType() == 
                    StateDescriptor::Interval) 
                {
                    plot_time = cur_time;
                }
                else
                {
                    int f_lev = parent->finestLevel();
                    plot_time = getLevel(f_lev).state[Press_Type].curTime();
                }
            }
            else
            {
                plot_time = cur_time;
            } 
	    const DeriveRec* rec = derive_lst.get(*it);
	    ncomp = rec->numDerive();
	    MultiFab* derive_dat = derive(*it,plot_time,nGrow);
	    MultiFab::Copy(plotMF,*derive_dat,0,cnt,ncomp,nGrow);
	    delete derive_dat;
	    cnt += ncomp;
	}
    }
    //
    // Cull data from diagnostic multifabs.
    //
    for (std::map<std::string,MultiFab*>::const_iterator it = auxDiag.begin(), end = auxDiag.end();
         it != end;
         ++it)
    {
        int nComp = it->second->nComp();
        MultiFab::Copy(plotMF,*it->second,0,cnt,nComp,nGrow);
        cnt += nComp;
    }
    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(plotMF,TheFullPath,how);
}

MultiFab*
HeatTransfer::derive (const std::string& name,
                      Real               time,
                      int                ngrow)
{        
  BL_ASSERT(ngrow >= 0);
  
  MultiFab* mf = 0;
  const DeriveRec* rec = derive_lst.get(name);
  if (rec)
  {
    BoxArray dstBA(grids);
    mf = new MultiFab(dstBA, rec->numDerive(), ngrow);
    int dcomp = 0;
    derive(name,time,*mf,dcomp);
  }
  else
  {
    mf = AmrLevel::derive(name,time,ngrow);
  }

  if (mf==0) {
    std::string msg("HeatTransfer::derive(): unknown variable: ");
    msg += name;
    BoxLib::Error(msg.c_str());
  }
  return mf;
}
 

void
HeatTransfer::derive (const std::string& name,
                      Real               time,
                      MultiFab&          mf,
                      int                dcomp)
{
    if (name == "mean_progress_curvature")
    {
        // 
        // Smooth the temperature, then use smoothed T to compute mean progress curvature
        //

        // Assert because we do not know how to un-convert the destination
        //   and also, implicitly assume the convert that was used to build mf BA is trivial
        BL_ASSERT(mf.boxArray()[0].ixType()==IndexType::TheCellType());

        const DeriveRec* rec = derive_lst.get(name);

        Box srcB = mf.boxArray()[0];
        Box dstB = rec->boxMap()(srcB);

        // Find nGrowSRC
        int nGrowSRC = 0;
        for ( ; !srcB.contains(dstB); ++nGrowSRC)
        {
            srcB.grow(1);
        }
        BL_ASSERT(nGrowSRC);  // Need grow cells for this to work!

        MultiFab tmf(mf.boxArray(),1,nGrowSRC);
        
        for (FillPatchIterator fpi(*this,tmf,nGrowSRC,time,State_Type,Temp,1);
             fpi.isValid();
             ++fpi)
        {
            tmf[fpi].copy(fpi());
        }

        int num_smooth_pre = 3;
        
        for (int i=0; i<num_smooth_pre; ++i)
        {
            // Fix up fine-fine and periodic
            tmf.FillBoundary(0,1,geom.periodicity());
                        
            for (MFIter mfi(tmf); mfi.isValid(); ++mfi)
            {
                // 
                // Use result MultiFab for temporary container to hold smooth T field
                // 
                const Box& box = mfi.validbox();
                FORT_SMOOTH(box.loVect(),box.hiVect(),
                            tmf[mfi].dataPtr(),
                            ARLIM(tmf[mfi].loVect()),ARLIM(tmf[mfi].hiVect()),
                            mf[mfi].dataPtr(dcomp),
                            ARLIM(mf[mfi].loVect()),ARLIM(mf[mfi].hiVect()));

                // Set result back into slot for smoothed T, leave grow cells
                tmf[mfi].copy(mf[mfi],box,dcomp,box,0,1);
            }
        }

        tmf.FillBoundary(0,1,geom.periodicity());

        const Real* dx = geom.CellSize();
        
        FArrayBox nWork;
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const FArrayBox& Tg = tmf[mfi];
            FArrayBox& MC = mf[mfi];
            const Box& box = mfi.validbox();
            const Box& nodebox = BoxLib::surroundingNodes(box);
            nWork.resize(nodebox,BL_SPACEDIM);
            
            FORT_MCURVE(box.loVect(),box.hiVect(),
                        Tg.dataPtr(),ARLIM(Tg.loVect()),ARLIM(Tg.hiVect()),
                        MC.dataPtr(dcomp),ARLIM(MC.loVect()),ARLIM(MC.hiVect()),
                        nWork.dataPtr(),ARLIM(nWork.loVect()),ARLIM(nWork.hiVect()),
                        dx);
        }
    } 
    else
    {
#ifdef PARTICLES
        ParticleDerive(name,time,mf,dcomp);
#else
        AmrLevel::derive(name,time,mf,dcomp);
#endif
    }
}
