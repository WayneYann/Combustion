//
// "Divu_Type" means S, where divergence U = S
// "Dsdt_Type" means pd S/pd t, where S is as above
// "Ydot_Type" means -omega_l/rho, i.e., the mass rate of decrease of species l due
//             to kinetics divided by rho
//
#include <winstd.H>

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstdio>
#include <cfloat>
#include <fstream>
#include <vector>

#include <CoordSys.H>
#include <Geometry.H>
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
#include <Profiler.H>
#include <ccse-mpi.H>

#ifndef NDEBUG
void inspectFAB (FArrayBox& unfab);
void inspectFAB (MultiFab& unfab);
void inspectFAB (MultiFab& unfab, int dummy, int comp);
void inspectFAB (MultiFab* unfab, int i);
#endif

#include <PROB_F.H>
#include <NAVIERSTOKES_F.H>
#include <VISCOPERATOR_F.H>
#include <DERIVE_F.H>

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

#define GEOM_GROW   1
#define HYP_GROW    3
#define PRESS_GROW  1
#define DIVU_GROW   1
#define DSDT_GROW   1
#define bogus_value 1.e20
#define DQRAD_GROW  1
#define YDOT_GROW   1
const int LinOp_grow = 1;

//
// Initialization of static members.
//
int       HeatTransfer::num_divu_iters        = 1;
int       HeatTransfer::init_once_done        = 0;
int       HeatTransfer::do_OT_radiation       = 0;
int       HeatTransfer::do_heat_sink          = 0;
int       HeatTransfer::RhoH                  = -1;
int       HeatTransfer::do_diffuse_sync       = 1;
int       HeatTransfer::do_reflux_visc        = 1;
int       HeatTransfer::dpdt_option           = 0;
int       HeatTransfer::Ydot_Type             = -1;
int       HeatTransfer::FuncCount_Type           = -1;
int       HeatTransfer::divu_ceiling          = 0;
Real      HeatTransfer::divu_dt_factor        = .5;
Real      HeatTransfer::min_rho_divu_ceiling  = -1.e20;
int       HeatTransfer::have_trac                 = 0;
int       HeatTransfer::have_rhort                = 0;
int       HeatTransfer::Trac                      = -1;
int       HeatTransfer::RhoRT                     = -1;
int       HeatTransfer::first_spec                = -1;
int       HeatTransfer::last_spec                 = -2;
int       HeatTransfer::nspecies                  = 0;
int       HeatTransfer::floor_species             = 0;
int       HeatTransfer::do_set_rho_to_species_sum = 1;
Real      HeatTransfer::rgas                      = -1.0;
Real      HeatTransfer::prandtl                   = .7;
Real      HeatTransfer::schmidt                   = .7;
Real      HeatTransfer::constant_mu_val           = -1;
Real      HeatTransfer::constant_rhoD_val         = -1;
Real      HeatTransfer::constant_lambda_val       = -1;
int       HeatTransfer::unity_Le                  = 1;
Real      HeatTransfer::htt_tempmin               = 298.0;
Real      HeatTransfer::htt_tempmax               = 40000.;
int       HeatTransfer::siegel_test               = 0;
int       HeatTransfer::zeroBndryVisc             = 0;
ChemDriver* HeatTransfer::chemSolve            = 0;
int       HeatTransfer::do_add_nonunityLe_corr_to_rhoh_adv_flux = 1;
int       HeatTransfer::do_check_divudt           = 1;
int       HeatTransfer::hack_nochem               = 0;
int       HeatTransfer::hack_nospecdiff           = 0;
int       HeatTransfer::hack_noavgdivu            = 0;
Real      HeatTransfer::trac_diff_coef            = 0.0;
Real      HeatTransfer::P1atm_MKS                 = -1.0;
std::string   HeatTransfer::turbFile                  ="";
ChemDriver::Chem_Evolve HeatTransfer::chem_integrator = ChemDriver::CKD_Vode;
Array<std::string> HeatTransfer::auxDiag_names(0);
bool      HeatTransfer::plot_auxDiags             = false;

static int  max_grid_size_chem   = 16;
static bool do_not_use_funccount = false;
static bool do_active_control    = false;
static Real crse_dt = -1;

#ifdef BL_USE_FLOAT
#  define Real_MIN FLT_MIN
#  define Real_MAX FLT_MAX
#else
#  define Real_MIN DBL_MIN
#  define Real_MAX DBL_MAX
#endif

#if defined(BL_PLOT_REACS) && defined(BL_PLOT_CONSUMPTION)
#error "Only one of BL_PLOT_REACS or BL_PLOT_CONSUMPTION may be defined"
#endif

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
    BL_PROFILE("HeatTransfer::center_to_edge_fancy()");

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
                for (BoxList::iterator bli = gCells.begin();
                     bli != gCells.end();
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
                for (BoxList::iterator bli = gCells.begin();
                     bli != gCells.end();
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
    NavierStokes::variableCleanUp();

    if (ParallelDescriptor::IOProcessor())
        std::cout << "Deleting chemSolver in variableCleanUp ... " << std::flush;

    delete chemSolve;

    if (ParallelDescriptor::IOProcessor())
        std::cout << "done\n";

    chemSolve = 0;
}

void
HeatTransfer::read_params ()
{
    //
    // Read parameters from input file and command line.
    //
    NavierStokes::read_params();

    ParmParse pp("ns");

    pp.query("do_diffuse_sync",do_diffuse_sync);
    BL_ASSERT(do_diffuse_sync == 0 || do_diffuse_sync == 1);
    pp.query("do_reflux_visc",do_reflux_visc);
    BL_ASSERT(do_reflux_visc == 0 || do_reflux_visc == 1);
    pp.query("dpdt_option",dpdt_option);
    BL_ASSERT(dpdt_option >= 0 && dpdt_option <= 2);
    pp.query("max_grid_size_chem",max_grid_size_chem);
    BL_ASSERT(max_grid_size_chem > 0);
    pp.query("do_active_control",do_active_control);

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
        if (ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::read_params: Le=1, setting Sc = Pr" << std::endl;
    }

    pp.query("constant_mu_val",constant_mu_val);
    pp.query("constant_rhoD_val",constant_rhoD_val);
    pp.query("constant_lambda_val",constant_lambda_val);
    if (constant_mu_val != -1)
    {
        if (ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::read_params: using constant_mu_val = " 
                      << constant_mu_val << std::endl;
    }
    if (constant_rhoD_val != -1)
    {
        if (ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::read_params: using constant_rhoD_val = " 
                      << constant_rhoD_val << std::endl;
    }
    if (constant_lambda_val != -1)
    {
        if (ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::read_params: using constant_lambda_val = " 
                      << constant_lambda_val << std::endl;
    }

    pp.query("do_add_nonunityLe_corr_to_rhoh_adv_flux",
             do_add_nonunityLe_corr_to_rhoh_adv_flux);
    pp.query("hack_nochem",hack_nochem);
    pp.query("hack_nospecdiff",hack_nospecdiff);
    pp.query("hack_noavgdivu",hack_noavgdivu);
    pp.query("do_check_divudt",do_check_divudt);

    pp.query("do_OT_radiation",do_OT_radiation);
    do_OT_radiation = (do_OT_radiation ? 1 : 0);
    pp.query("do_heat_sink",do_heat_sink);
    do_heat_sink = (do_heat_sink ? 1 : 0);

    int use_chemeq2 = 0; pp.query("use_chemeq2",use_chemeq2);
    chem_integrator = (use_chemeq2 ? ChemDriver::CKD_ChemEQ2 : ChemDriver::CKD_Vode );

    std::string tranfile=""; pp.query("tranfile",tranfile);
    chemSolve = new ChemDriver(tranfile);

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
    auxDiag                = 0;
}

HeatTransfer::HeatTransfer (Amr&            papa,
                            int             lev,
                            const Geometry& level_geom,
                            const BoxArray& bl,
                            Real            time)
    :
    NavierStokes(papa,lev,level_geom,bl,time),
    //
    // Make room for all components except velocities in aux_boundary_data_old.
    //
    aux_boundary_data_old(bl,HYP_GROW,desc_lst[State_Type].nComp()-BL_SPACEDIM,level_geom),
    //
    // Only save Density & RhoH in aux_boundary_data_new in components 0 & 1.
    //
    aux_boundary_data_new(bl,LinOp_grow,2,level_geom),
    FillPatchedOldState_ok(true)
{
    if (!init_once_done)
        init_once();

    if (!do_temp)
        BoxLib::Abort("do_temp MUST be true");

    if (!have_divu)
        BoxLib::Abort("have_divu MUST be true");

    if (!have_dsdt)
        BoxLib::Abort("have_dsdt MUST be true");

    const int nGrow       = 0;
    const int nEdgeStates = desc_lst[State_Type].nComp();
    diffusion->allocFluxBoxesLevel(EdgeState,nGrow,nEdgeStates);
    if (nspecies>0 && !unity_Le)
    {
	diffusion->allocFluxBoxesLevel(SpecDiffusionFluxn,  nGrow,nspecies);
	diffusion->allocFluxBoxesLevel(SpecDiffusionFluxnp1,nGrow,nspecies);
	spec_diffusion_flux_computed.resize(nspecies,HT_None);
    }

    auxDiag = 0;

    if (plot_auxDiags)
    {
        BL_ASSERT(auxDiag_names.size() > 0);
        auxDiag = new MultiFab(grids,auxDiag_names.size(),0);
        auxDiag->setVal(0);
    }
}

HeatTransfer::~HeatTransfer ()
{
    diffusion->removeFluxBoxesLevel(EdgeState);    
    if (nspecies>0 && !unity_Le)    
    {
	diffusion->removeFluxBoxesLevel(SpecDiffusionFluxn);    
	diffusion->removeFluxBoxesLevel(SpecDiffusionFluxnp1);    
    }
    delete auxDiag;
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
    
    if (!have_rhort && ParallelDescriptor::IOProcessor())
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
    //
    // Get universal gas constant from Fortran.
    //
    rgas = getChemSolve().getRuniversal();
    P1atm_MKS = getChemSolve().getP1atm_MKS();

    if (rgas <= 0.0)
    {
        std::cout << "HeatTransfer::init_once(): bad rgas: " << rgas << '\n';
        BoxLib::Abort();
    }
    if (P1atm_MKS <= 0.0)
    {
        std::cout << "HeatTransfer::init_once(): bad P1atm_MKS: " << P1atm_MKS << '\n';
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
        if (ParallelDescriptor::IOProcessor())
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

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "HeatTransfer::init_once(): num_state_type = "
                  << num_state_type << '\n';
    }
    ParmParse pp("ht");

    pp.query("plot_auxDiags",plot_auxDiags);

    if (plot_auxDiags)
    {
#ifdef BL_PLOT_REACS
        auxDiag_names.resize(getChemSolve().numReactions());
        for (int i = 0; i < auxDiag_names.size(); ++i)
            auxDiag_names[i] = BoxLib::Concatenate("R",i+1);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "***** Make sure to increase amr.regrid_int !!!!!" << std::endl;
#elif defined(BL_PLOT_CONSUMPTION)
        auxDiag_names.resize(1);
        auxDiag_names[0] = "FuelConsumptionRate";
#endif
    }

    init_once_done = 1;
}

void
HeatTransfer::restart (Amr&          papa,
                       std::istream& is,
                       bool          bReadSpecial)
{

    NavierStokes::restart(papa,is,bReadSpecial);
    //
    // Make room for all components except velocities in aux_boundary_data_old.
    //
    aux_boundary_data_old.initialize(grids,HYP_GROW,desc_lst[State_Type].nComp()-BL_SPACEDIM,Geom());
    //
    // Only save Density & RhoH in aux_boundary_data_new in components 0 & 1.
    //
    aux_boundary_data_new.initialize(grids,LinOp_grow,2,Geom());

    FillPatchedOldState_ok = true;

    set_overdetermined_boundary_cells(state[State_Type].curTime());

    BL_ASSERT(EdgeState == 0);
    const int nGrow       = 0;
    const int nEdgeStates = desc_lst[State_Type].nComp();
    diffusion->allocFluxBoxesLevel(EdgeState,nGrow,nEdgeStates);
    
    if (nspecies>0 && !unity_Le)
    {
	BL_ASSERT(SpecDiffusionFluxn == 0);
	BL_ASSERT(SpecDiffusionFluxnp1 == 0);
	diffusion->allocFluxBoxesLevel(SpecDiffusionFluxn,  nGrow,nspecies);
	diffusion->allocFluxBoxesLevel(SpecDiffusionFluxnp1,nGrow,nspecies);
	spec_diffusion_flux_computed.resize(nspecies,HT_None);
    }

    auxDiag = 0;

    if (plot_auxDiags)
    {
        BL_ASSERT(auxDiag_names.size() > 0);
        auxDiag = new MultiFab(grids,auxDiag_names.size(),0);
        auxDiag->setVal(0);
    }
}

Real
HeatTransfer::estTimeStep ()
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::estTimeStep()");

    Real estdt = NavierStokes::estTimeStep();

    if (fixed_dt > 0.0 || !divu_ceiling)
        //
        // The n-s function did the right thing in this case.
        //
        return estdt;

    Real dt, ns_estdt = estdt, divu_dt = 1.0e20;

    const int         n_grow   = 1;
    const Real        cur_time = state[State_Type].curTime();
    const Real*       dx       = geom.CellSize();
    MultiFab*         dsdt     = getDsdt(0,cur_time);
    MultiFab*         divu     = getDivCond(0,cur_time);

    MFIter  dsdtmfi(*dsdt);
    MFIter  divumfi(*divu);
    MFIter  Rho_mfi(*rho_ctime);

    for (FillPatchIterator U_fpi(*this,*divu,n_grow,cur_time,State_Type,Xvel,BL_SPACEDIM);
         dsdtmfi.isValid() && divumfi.isValid() && Rho_mfi.isValid() && U_fpi.isValid();
         ++dsdtmfi, ++divumfi, ++Rho_mfi, ++U_fpi)
    {
        BL_ASSERT(dsdtmfi.index() == divumfi.index());

        const int        i   = dsdtmfi.index();
        FArrayBox&       U   = U_fpi();
        const FArrayBox& Rho = (*rho_ctime)[Rho_mfi];
        const int*       lo  = grids[i].loVect();
        const int*       hi  = grids[i].hiVect();

        DEF_CLIMITS((*divu)[Rho_mfi],sdat,slo,shi);
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
                         (*dsdt)[Rho_mfi].dataPtr(),
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

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "HeatTransfer::estTimeStep(): estdt, divu_dt = " 
                  << estdt << ", " << divu_dt << '\n';
    }

    estdt = std::min(estdt, divu_dt);

    if (estdt < ns_estdt && ParallelDescriptor::IOProcessor())
    {
        std::cout << "HeatTransfer::estTimeStep() : timestep reduced from " 
                  << ns_estdt << " to " << estdt << '\n';
    }

    return estdt;
}

void
HeatTransfer::checkTimeStep (Real dt)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::checkTimeStep()");

    if (fixed_dt > 0.0 || !divu_ceiling) 
        return;

    const int   n_grow    = 1;
    const Real  cur_time  = state[State_Type].curTime();
    const Real* dx        = geom.CellSize();
    MultiFab*   dsdt      = getDsdt(0,cur_time);
    MultiFab*   divu      = getDivCond(0,cur_time);

    MFIter dsdtmfi(*dsdt);
    MFIter divumfi(*divu);

    for (FillPatchIterator U_fpi(*this,*divu,n_grow,cur_time,State_Type,Xvel,BL_SPACEDIM);
         dsdtmfi.isValid() && divumfi.isValid() && U_fpi.isValid();
         ++dsdtmfi, ++divumfi, ++U_fpi)
    {
        BL_ASSERT(dsdtmfi.index() == divumfi.index());

        const int        i   = divumfi.index();
        FArrayBox&       U   = U_fpi();
        const FArrayBox& Rho = (*rho_ctime)[i];
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
    NavierStokes::setTimeLevel(time, dt_old, dt_new);    

    state[Ydot_Type].setTimeLevel(time,dt_old,dt_new);

    state[FuncCount_Type].setTimeLevel(time,dt_old,dt_new);
}

void
HeatTransfer::initDataOtherTypes ()
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::initDataOtherTypes()");

    const Real* dx = geom.CellSize();

    int add_turb = 0, box_offset[BL_SPACEDIM];
    Real turb_size[BL_SPACEDIM], turb_scale;
    //
    // FORT_CHECK_TURB should be located in PROB_?D.F.
    // If it is not, merely comment out the call below and everything
    // should be ok.
    //
    FORT_CHECK_TURB(&add_turb,turb_size,box_offset,&turb_scale,dx);
  
    if (add_turb) 
    {
        MultiFab& S_new = get_new_data(State_Type);

        const int numVelComp = BL_SPACEDIM;
        const int numVelGrow = 0;
        //
        // Read in data about turbulent velocity field.
        //
        Real fileProbSize[BL_SPACEDIM];
        int vel_size[BL_SPACEDIM];
        {
            std::string dxFileName = turbFile + "/Size";
            std::ifstream is(dxFileName.c_str(),std::ios::in);
            if (is.fail())
            {
                std::string msg("HeatTransfer: cannot open velocity file: ");
                msg += turbFile;
                BoxLib::Error(msg.c_str());
            }
            for (int i = 0; i < BL_SPACEDIM; i++)
                is >> fileProbSize[i];
            for (int i = 0; i < BL_SPACEDIM; i++)
                is >> vel_size[i];
        }
        //
        // Make sure turbulent velocity field is correct size.
        //
        if (D_TERM(fileProbSize[0] != turb_size[0],
                   || fileProbSize[1] != turb_size[1],
                   || fileProbSize[2] != turb_size[2]))
        {
            if (ParallelDescriptor::IOProcessor())
                std::cout << "warning turb field is wrong size" << '\n';
        }
        //
        // Read in turbulent velocity field ... store in turbVelFile.
        //
        std::string MFPath = turbFile +"/MultiFab";
        MultiFab* turbVelFile = new MultiFab();
        VisMF::Read(*turbVelFile,MFPath);
        //
        // Compute refinement ratio between turbulence file data and this level.
        // 
        bool    do_avgdown;
        IntVect avgdown_ratio(IntVect::TheZeroVector());
        IntVect refine_ratio(IntVect::TheZeroVector());

        for (int i = 0;i < BL_SPACEDIM; i++)
        {
            Real dxFile = fileProbSize[i]/Real(vel_size[i]);
            if (dxFile > geom.CellSize(i))
            {
                refine_ratio[i] = int(dxFile/geom.CellSize(i));
                do_avgdown = false;
            }
            else
            {
                avgdown_ratio[i] = int(geom.CellSize(i)/dxFile);
                do_avgdown = true;
            }
        }
        //
        // Check for consistency in ratios.
        //
        if (do_avgdown) 
        {
            if (D_TERM(avgdown_ratio[0] < 1, 
                       || avgdown_ratio[1] < 1, 
                       || avgdown_ratio[2] < 1)
                || refine_ratio != IntVect::TheZeroVector())
                BoxLib::Error("bad sizing");
        }
        else
        {
            if (D_TERM(refine_ratio[0] < 1, 
                       || refine_ratio[1] < 1, 
                       || refine_ratio[2] < 1)
                || avgdown_ratio != IntVect::TheZeroVector())
                BoxLib::Error("bad sizing");
        }

        MultiFab* gridTurbVel;

        if (do_avgdown) 
        {
            //
            // Avg down data from fine level to this level.
            // The case where they are on the same level falls through here.
            //
            if (ParallelDescriptor::IOProcessor())
                std::cout << "averaging down turb data by factor of " 
                          << avgdown_ratio << std::endl;

            BoxArray fileTurbBA(turbVelFile->boxArray());
            BoxArray crseTurbBA(fileTurbBA);
            crseTurbBA.coarsen(avgdown_ratio);
            gridTurbVel = new MultiFab(crseTurbBA,numVelComp,numVelGrow);
	  
            Box fineDomain = BoxLib::refine(geom.Domain(),avgdown_ratio);
            Geometry fineGeom(fineDomain);
            MultiFab* crseVol = new MultiFab();
            MultiFab* fineVol = new MultiFab();
            geom.GetVolume(*crseVol,crseTurbBA,numVelGrow);
            fineGeom.GetVolume(*fineVol,fileTurbBA,numVelGrow);
	  
            for (MFIter finemfi(*turbVelFile); finemfi.isValid(); ++finemfi)
            {
                const int i       = finemfi.index();
                const int f_level = 0;
                const int c_level = 0;
                const Box ovlp    = (*gridTurbVel).box(i);
	      
                NavierStokes::avgDown_doit((*turbVelFile)[finemfi],
                                           (*gridTurbVel)[finemfi],
                                           (*fineVol)[finemfi],
                                           (*crseVol)[finemfi],
                                           f_level,c_level,
                                           ovlp,Xvel,numVelComp,avgdown_ratio);
            }

            delete crseVol;
            delete fineVol;
        }
        else
        {
            //
            // Interpolate data from this level to finer level.
            //
            if (ParallelDescriptor::IOProcessor())
                std::cout << "interpolating turb data by factor of " 
                          << refine_ratio << std::endl;

            if (turbVelFile->nGrow() == 0)
                BoxLib::Error("need ghost cells to interp -- use new turb fcn");
            BoxArray fileTurbBA(turbVelFile->boxArray());
            BoxArray gridTurbBA(fileTurbBA);
            gridTurbBA.refine(refine_ratio);
            gridTurbVel = new MultiFab(gridTurbBA,numVelComp,numVelGrow);

            Box fineDomain = BoxLib::refine(geom.Domain(),refine_ratio);
            Geometry fineGeom(fineDomain);
            const Geometry& crseGeom = geom;
            Array<BCRec> bcr(numVelComp);
            Interpolater* mapper = &pc_interp;
            const Box& pdomain = state[State_Type].getDomain();
            const StateDescriptor&  desc    = desc_lst[State_Type];

            for (MFIter crsemfi(*turbVelFile); crsemfi.isValid(); ++crsemfi)
            {
                const int  i    = crsemfi.index();
                const Box& dbox = gridTurbBA[i];
                BoxLib::setBC(dbox,pdomain,Xvel,Xvel,numVelComp,desc.getBCs(),bcr);
                mapper->interp((*turbVelFile)[crsemfi],Xvel,
                               (*gridTurbVel)[crsemfi],Xvel,
                               numVelComp, (*gridTurbVel).box(i),
                               refine_ratio,
                               crseGeom,fineGeom,
                               bcr);
            }
        }

        delete turbVelFile;
        //
        // Shift grid velocity by offset.
        //
        BoxArray gridTurbBA = gridTurbVel->boxArray();
        BoxArray shiftedGridTurbBA(gridTurbBA);
        for (int dir=0; dir < BL_SPACEDIM; dir++)
            shiftedGridTurbBA.shift(dir,box_offset[dir]);
        MultiFab* shiftedGridTurbVel
            = new MultiFab(shiftedGridTurbBA,numVelComp,numVelGrow);

        for (MFIter gridmfi(*gridTurbVel); gridmfi.isValid(); ++gridmfi)
        {

            (*shiftedGridTurbVel)[gridmfi].copy((*gridTurbVel)[gridmfi],
                                                gridmfi.validbox(),Xvel,
                                                (*shiftedGridTurbVel).box(gridmfi.index()),
                                                Xvel,numVelComp);
        }

        delete gridTurbVel;
        //
        // Copy vel into boxarray used by S_new,
        // multiply by turbulence scale and then add to S_new.
        //
        MultiFab* Vel = new MultiFab(S_new.boxArray(),numVelComp,numVelGrow);
        Vel->setVal(0.0);
        Vel->copy(*shiftedGridTurbVel,Xvel,Xvel,numVelComp);

        Vel->mult(turb_scale,Xvel,numVelComp,numVelGrow);

        S_new.plus(*Vel,Xvel,numVelComp,numVelGrow);

        delete shiftedGridTurbVel;
        delete Vel;
    }

    const Real cur_time  = state[State_Type].curTime();
    {
        MultiFab rhoh(grids,1,0);
        compute_rhohmix(cur_time,rhoh);
        get_new_data(State_Type).copy(rhoh,0,RhoH,1);
    }
    //
    // Set up diffusivities, viscosities (need for initial divu compute)
    //
    const Real dt        = -1.0;
    const int  iteration = 0;
    const int  ncycle    = 0;
    //
    // Assume always variable diffusivity.
    // Assume always variable viscosity.
    //
    const int offset   = BL_SPACEDIM + 1;
    const int num_diff = NUM_STATE-offset;
    calcDiffusivity(cur_time, dt, iteration, ncycle, offset, num_diff, true);
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
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::init(AmrLEvel&)");

    NavierStokes::init(old);

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
        Ydot[fpi.index()].copy(fpi(),0,0,nspecies);
    }

    RhoH_to_Temp(get_new_data(State_Type));

    MultiFab& FuncCount = get_new_data(FuncCount_Type);

    for (FillPatchIterator fpi(*oldht,FuncCount,FuncCount.nGrow(),cur_time,FuncCount_Type,0,1);
         fpi.isValid();
         ++fpi)
    {
        FuncCount[fpi.index()].copy(fpi());
    }
}

//
// Inits the data on a new level that did not exist before regridding.
//

void
HeatTransfer::init ()
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::init()");

    NavierStokes::init();
 
    HeatTransfer& old      = getLevel(level-1);
    const Real    cur_time = old.state[State_Type].curTime();
    //
    // Get best ydot data.
    //
    FillCoarsePatch(get_new_data(Ydot_Type),0,cur_time,Ydot_Type,0,nspecies);

    RhoH_to_Temp(get_new_data(State_Type));

    FillCoarsePatch(get_new_data(FuncCount_Type),0,cur_time,FuncCount_Type,0,1);
}

void
HeatTransfer::post_timestep (int crse_iteration)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::post_timestep()");

    NavierStokes::post_timestep(crse_iteration);

    if (plot_auxDiags && level == 0)
    {
        const int Ndiag = auxDiag->nComp();

#ifdef BL_PLOT_REACS
        //
        // Multiply by the inverse of the coarse timestep.
        //
        const Real factor = 1.0 / crse_dt;

        for (int i = parent->finestLevel(); i >= 0; --i)
            getLevel(i).auxDiag->mult(factor);
#endif

        for (int i = parent->finestLevel(); i > 0; --i)
        {
            HeatTransfer& clev = getLevel(i-1);
            HeatTransfer& flev = getLevel(i);

            MultiFab& Ydot_crse = *(clev.auxDiag);
            MultiFab& Ydot_fine = *(flev.auxDiag);
            
            NavierStokes::avgDown(clev.boxArray(),flev.boxArray(),
                                  Ydot_crse,Ydot_fine,
                                  clev.volume,flev.volume,
                                  i-1,i,0,Ndiag,parent->refRatio(i-1));
        }
    }
}

void
HeatTransfer::post_restart ()
{
    NavierStokes::post_restart();

    Real dummy  = 0;
    int MyProc  = ParallelDescriptor::MyProc();
    int step    = parent->levelSteps(0);
    int restart = 1;

    if (do_active_control)
        FORT_ACTIVECONTROL(&dummy,&dummy,&crse_dt,&MyProc,&step,&restart);
}

void
HeatTransfer::post_regrid (int lbase,
                           int new_finest)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::post_regrid()");

    NavierStokes::post_regrid(lbase, new_finest);
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
HeatTransfer::post_init (Real stop_time)
{
    if (level > 0)
        //
        // Nothing to sync up at level > 0.
        //
        return;

    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::post_init()");

    const Real cur_time     = state[State_Type].curTime();
    const int  finest_level = parent->finestLevel();
    Real        dt_init     = 0.0;

    Array<Real> dt_save(finest_level+1);
    Array<int>  nc_save(finest_level+1);
    Real        dt_init2 = 0.0;
    Array<Real> dt_save2(finest_level+1);
    Array<int>  nc_save2(finest_level+1);
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
      std::cout << "doing num_divu_iters = " << num_divu_iters << std::endl;

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
 
                S_tmp.copy(S_new);
 
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
                    const BoxArray& fgrids   = fine_lev.grids;
                    MultiFab&       fvolume  = fine_lev.volume;
                    
                    HeatTransfer&   crse_lev = getLevel(k);
                    const BoxArray& cgrids   = crse_lev.grids;
                    MultiFab&       cvolume  = crse_lev.volume;
                    const IntVect&  fratio   = crse_lev.fine_ratio;
                    
                    MultiFab& Divu_crse = crse_lev.get_new_data(Divu_Type);
                    MultiFab& Divu_fine = fine_lev.get_new_data(Divu_Type);
                    
                    crse_lev.NavierStokes::avgDown(cgrids,fgrids,
                                                   Divu_crse,Divu_fine,
                                                   cvolume,fvolume,
                                                   k,k+1,0,1,fratio);
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
                const BoxArray& fgrids  = getLevel(k+1).grids;
                const BoxArray& cgrids  = getLevel(k).grids;
                MultiFab&       fvolume = getLevel(k+1).volume;
                MultiFab&       cvolume = getLevel(k).volume;
                MultiFab&       S_fine  = getLevel(k+1).get_new_data(State_Type);
                MultiFab&       S_crse  = getLevel(k).get_new_data(State_Type);
                IntVect&        fratio  = getLevel(k).fine_ratio;
                
                getLevel(k).NavierStokes::avgDown(cgrids,fgrids,
                                                  S_crse,S_fine,
                                                  cvolume,fvolume,
                                                  k,k+1,Xvel,BL_SPACEDIM,
                                                  fratio);
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
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::sum_integrated_quantities()");

    const int finest_level = parent->finestLevel();
    const Real time        = state[State_Type].curTime();

    Real mass = 0.0;
    for (int lev = 0; lev <= finest_level; lev++)
	mass += getLevel(lev).volWgtSum("density",time);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "TIME= " << time << " MASS= " << mass;

    if (getChemSolve().index(fuelName) >= 0)
    {
        Real fuelmass = 0.0;
        std::string fuel = "rho.Y(" + fuelName + ")";
        for (int lev = 0; lev <= finest_level; lev++)
            fuelmass += getLevel(lev).volWgtSum(fuel,time);

        if (ParallelDescriptor::IOProcessor())
            std::cout << " FUELMASS= " << fuelmass;

        int MyProc  = ParallelDescriptor::MyProc();
        int step    = parent->levelSteps(0);
        int restart = 0;

        if (do_active_control)
            FORT_ACTIVECONTROL(&fuelmass,&time,&crse_dt,&MyProc,&step,&restart);
    }

    if (ParallelDescriptor::IOProcessor())
        std::cout << '\n';

    Real rho_h    = 0.0;
    Real rho_temp = 0.0;

    for (int lev = 0; lev <= finest_level; lev++)
    {
        rho_temp += getLevel(lev).volWgtSum("rho_temp",time);
        rho_h     = getLevel(lev).volWgtSum("rhoh",time);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "TIME= " << time << " LEV= " << lev << " RHOH= " << rho_h << '\n';
    }
    if (ParallelDescriptor::IOProcessor())
        std::cout << "TIME= " << time << " RHO*TEMP= " << rho_temp << '\n';

    int old_prec = std::cout.precision(25);

    Real min_temp = 1.0e20;
    Real max_temp = -1.0e20;
    Real min_xvel = 1.0e20;
    Real max_xvel = -1.0e20;
    Real min_yvel = 1.0e20;
    Real max_yvel = -1.0e20;
#if (BL_SPACEDIM == 3)
    Real min_zvel = 1.0e20;
    Real max_zvel = -1.0e20;
#endif

    for (int lev = 0; lev <= finest_level; lev++)
    {
	MultiFab& newstate = getLevel(lev).get_data(State_Type,time);

        min_temp = std::min(min_temp,newstate.min(Temp,0));
        max_temp = std::max(max_temp,newstate.max(Temp,0));
	min_xvel = std::min(min_xvel,newstate.min(Xvel,0));
	max_xvel = std::max(max_xvel,newstate.max(Xvel,0));
	min_yvel = std::min(min_yvel,newstate.min(Yvel,0));
	max_yvel = std::max(max_yvel,newstate.max(Yvel,0));
#if (BL_SPACEDIM == 3)
	min_zvel = std::min(min_zvel,newstate.min(Zvel,0));
	max_zvel = std::max(max_zvel,newstate.max(Zvel,0));
#endif
    }

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "min,max temp = " << min_temp << ", " << max_temp << '\n';

        std::cout << "min,max xvel = " << min_xvel << ", " << max_xvel << '\n';
        std::cout << "min,max yvel = " << min_yvel << ", " << max_yvel << '\n';
#if (BL_SPACEDIM == 3)
        std::cout << "min,max zvel = " << min_zvel << ", " << max_zvel << '\n';
#endif
    }

    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    if (nspecies > 0)
    {
	Real min_sum = 1.0e20;
	Real max_sum = -1.0e20;

	for (int lev = 0; lev <= finest_level; lev++)
	{
            MultiFab* mf = getLevel(lev).derive("rhominsumrhoY",time,0);

            for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
	    {
		min_sum = std::min(min_sum,(*mf)[mfi].min());
		max_sum = std::max(max_sum,(*mf)[mfi].max());
	    }

            delete mf;
	}

        ParallelDescriptor::ReduceRealMax(max_sum,IOProc);
        ParallelDescriptor::ReduceRealMin(min_sum,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "min,max rho-sum rho Y_l = " << min_sum << ", " << max_sum << '\n';
    }
    
    Real min_sum  = 1.0e20;
    Real max_sum  = -1.0e20;

    for (int lev = 0; lev <= finest_level; lev++)
    {
        MultiFab* mf = getLevel(lev).derive("sumYdot",time,0);

        for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
        {
            min_sum = std::min(min_sum,(*mf)[mfi].min());
            max_sum = std::max(max_sum,(*mf)[mfi].max());
        }

        delete mf;
    }

    ParallelDescriptor::ReduceRealMax(max_sum,IOProc);
    ParallelDescriptor::ReduceRealMin(min_sum,IOProc);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "min,max sum Ydot = " << min_sum << ", " << max_sum << '\n';

    std::cout << std::setprecision(old_prec);
}

void
HeatTransfer::post_init_press (Real&        dt_init,
			       Array<int>&  nc_save,
			       Array<Real>& dt_save)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::post_init_press()");

    const int  nState          = desc_lst[State_Type].nComp();
    const int  nGrow           = 0;
    const Real cur_time        = state[State_Type].curTime();
    const int  finest_level    = parent->finestLevel();
    MultiFab&  P_new           = get_new_data(Press_Type);
    MultiFab&  P_old           = get_old_data(Press_Type);
    NavierStokes::initial_iter = true;
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
            sig[k] = getLevel(k).get_rho_half_time();
        }

        if (projector)
        {
            int havedivu = 1;
            projector->initialSyncProject(0,sig,parent->dtLevel(0),cur_time,
                                          dt_init,havedivu);
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

        NavierStokes::initial_iter = false;
    }

    if (init_iter <= 0)
        NavierStokes::initial_iter = false; // Just being compulsive -- rbp.

    NavierStokes::initial_step = false;
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
    NavierStokes::resetState(time,dt_old,dt_new);

    state[Ydot_Type].reset();
    state[Ydot_Type].setTimeLevel(time,dt_old,dt_new);

    state[FuncCount_Type].reset();
    state[FuncCount_Type].setTimeLevel(time,dt_old,dt_new);
}

void
HeatTransfer::avgDown ()
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::avgDown()");

    if (level == parent->finestLevel()) return;

    HeatTransfer&   fine_lev = getLevel(level+1);
    const BoxArray& fgrids   = fine_lev.grids;
    MultiFab&       fvolume  = fine_lev.volume;
    MultiFab&       S_crse   = get_new_data(State_Type);
    MultiFab&       S_fine   = fine_lev.get_new_data(State_Type);

    NavierStokes::avgDown(grids,fgrids,S_crse,S_fine,volume,fvolume,
                          level,level+1,0,S_crse.nComp(),fine_ratio);
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
    MultiFab&       P_fine_avg  = *fine_lev.p_avg;
    MultiFab&       P_fine      = initial_step ? P_fine_init : P_fine_avg;
    const BoxArray& P_cgrids    = state[Press_Type].boxArray();
    const BoxArray& P_fgrids    = fine_lev.state[Press_Type].boxArray();

    BoxArray crse_P_fine_BA(P_fgrids.size());

    for (int i = 0; i < P_fgrids.size(); ++i)
    {
        crse_P_fine_BA.set(i,BoxLib::coarsen(P_fgrids[i],fine_ratio));
    }

    {
        MultiFab crse_P_fine(crse_P_fine_BA,1,0);

        for (MFIter mfi(P_fine); mfi.isValid(); ++mfi)
        {
            const int i = mfi.index();

            injectDown(crse_P_fine_BA[i],crse_P_fine[i],P_fine[i],fine_ratio);
        }

        P_crse.copy(crse_P_fine);
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
            
        NavierStokes::avgDown(grids,fgrids,
                              Divu_crse,Divu_fine,
                              volume,fvolume,
                              level,level+1,0,1,fine_ratio);
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
            
        NavierStokes::avgDown(grids,fgrids,
                              Dsdt_crse,Dsdt_fine,
                              volume,fvolume,
                              level,level+1,0,1,fine_ratio);
    }
}

void
HeatTransfer::scalar_diffusion_update (Real dt,
                                       int  first_scalar, 
                                       int  last_scalar,
                                       int  corrector)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::scalar_diffusion_update()");

    const Real strt_time = ParallelDescriptor::second();
    //
    // Build single component edge-centered array of MultiFabs for fluxes
    //
    MultiFab** fluxSCn;
    MultiFab** fluxSCnp1;
    const int nGrow = 0;
    const int nComp = 1;
    diffusion->allocFluxBoxesLevel(fluxSCn  ,nGrow,nComp);
    diffusion->allocFluxBoxesLevel(fluxSCnp1,nGrow,nComp);
    //
    // Set diffusion solve mode.
    //
    Diffusion::SolveMode solve_mode = Diffusion::ONEPASS;
    //
    // Do implicit c-n solve for each scalar.
    //
    const MultiFab* Rh = get_rho_half_time();

    for (int sigma = first_scalar; sigma <= last_scalar; sigma++)
    {
        if (is_diffusive[sigma])
        {
            int rho_flag = 0;
	    
            MultiFab *delta_rhs, *alpha, **betan, **betanp1;

            diffuse_scalar_setup(dt,sigma,&rho_flag,delta_rhs,alpha,betan,betanp1);
            //
            // Note: dt taken care of in diffuse_scalar.
            //
	    const int dataComp = 0; // Component of dR, alpha, betas to use.

            diffusion->diffuse_scalar(dt,sigma,be_cn_theta,Rh,rho_flag,fluxSCn,
                                      fluxSCnp1,dataComp,delta_rhs,alpha,betan,
                                      betanp1,solve_mode);
	    //
	    // If corrector, increment the viscous flux registers (assume
	    //  corrector called only ONCE!).
	    //
	    if (do_reflux && corrector)
	    {
                FArrayBox fluxtot;

		for (int d = 0; d < BL_SPACEDIM; d++)
		{
		    for (MFIter fmfi(*fluxSCn[d]); fmfi.isValid(); ++fmfi)
		    {
                        const Box& ebox = (*fluxSCn[d])[fmfi].box();
                        fluxtot.resize(ebox,nComp);
                        fluxtot.copy((*fluxSCn[d])[fmfi],ebox,0,ebox,0,nComp);
                        fluxtot.plus((*fluxSCnp1[d])[fmfi],ebox,0,0,nComp);
			if (level < parent->finestLevel())
			    getLevel(level+1).getViscFluxReg().CrseInit(fluxtot,ebox,
									d,0,sigma,
									nComp,-dt);
			
			if (level > 0)
			    getViscFluxReg().FineAdd(fluxtot,d,fmfi.index(),
						     0,sigma,nComp,dt);
		    }
		}
		if (level < parent->finestLevel())
		    getLevel(level+1).getViscFluxReg().CrseInitFinish();
	    }
	    //
	    // Clean up memory, etc
	    //
            diffuse_cleanup(delta_rhs, betan, betanp1, alpha);
        }  
    }    
    diffusion->removeFluxBoxesLevel(fluxSCn);
    diffusion->removeFluxBoxesLevel(fluxSCnp1);

    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_time = ParallelDescriptor::second() - strt_time;

    ParallelDescriptor::ReduceRealMax(run_time,IOProc);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "HeatTransfer::scalar_diffusion_update(): time: " << run_time << std::endl;
}

void
HeatTransfer::differential_spec_diffusion_update (Real dt,
						  int  corrector)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::differential_spec_diffusion_update()");

    const Real strt_time = ParallelDescriptor::second();

    if (hack_nospecdiff)
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "... HACK!!! skipping spec diffusion " << std::endl;

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
    diffusion->allocFluxBoxesLevel(fluxSCn,  nGrow,nCompSC);
    diffusion->allocFluxBoxesLevel(fluxSCnp1,nGrow,nCompSC);
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
    MultiFab **betan, **betanp1;
    diffusion->allocFluxBoxesLevel(betan  ,nGrow,nCompY);
    diffusion->allocFluxBoxesLevel(betanp1,nGrow,nCompY);
    Array<int> rho_flag(nCompY,0);
    
    MultiFab *alphaSC, *delta_rhsSC, **betanSC, **betanp1SC;

    const MultiFab* Rh = get_rho_half_time();

    for (int sigma = 0; sigma < nCompY; ++sigma)
    {
	const int state_ind = sCompY + sigma;

	diffuse_scalar_setup(dt, state_ind, &rho_flag[sigma], delta_rhsSC,
			     alphaSC, betanSC, betanp1SC);
	
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
	//
	// Clean up single-component stuff, then diffuse the scalar
	//
	diffuse_cleanup(delta_rhsSC, betanSC, betanp1SC, alphaSC);

	diffusion->diffuse_scalar(dt,state_ind,be_cn_theta,Rh,rho_flag[sigma],
                                  fluxSCn,fluxSCnp1,sigma,&delta_rhs,alpha,
                                  betan,betanp1,solve_mode);
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
    //
    // Modify update/fluxes to preserve flux sum = 0, compute new update and
    // leave modified fluxes in level data.  Do this in two stages, first for
    // the explicit fluxes, then the implicit ones (so send in rhs=0 for the
    // second one...).
    //
    const int  dataComp  = 0;
    const Real prev_time = state[State_Type].prevTime();
    const Real cur_time  = state[State_Type].curTime();

    adjust_spec_diffusion_update(get_new_data(State_Type),&get_old_data(State_Type),
				 sCompY,dt,prev_time,rho_flag,Rh,dataComp,
                                 &delta_rhs,alpha,betan);
    adjust_spec_diffusion_update(get_new_data(State_Type),&get_new_data(State_Type),
				 sCompY,dt,cur_time,rho_flag,Rh,dataComp,0,
                                 alpha,betanp1);
    //
    // Now do reflux with new, improved fluxes
    //
    if (do_reflux && corrector)
    {
        FArrayBox fluxtot;

	for (int d = 0; d < BL_SPACEDIM; d++)
	{
	    for (MFIter fmfi(*SpecDiffusionFluxn[d]); fmfi.isValid(); ++fmfi)
	    {
                const Box& ebox = (*SpecDiffusionFluxn[d])[fmfi].box();

                fluxtot.resize(ebox,nCompY);
                fluxtot.copy((*SpecDiffusionFluxn[d])[fmfi], ebox,0,ebox,0,nCompY);
                fluxtot.plus((*SpecDiffusionFluxnp1[d])[fmfi],ebox,0,0,nCompY);

		if (level < parent->finestLevel())
		    getLevel(level+1).getViscFluxReg().CrseInit(fluxtot,ebox,d,0,sCompY,nCompY,-dt);

		if (level > 0)
		    getViscFluxReg().FineAdd(fluxtot,d,fmfi.index(),0,sCompY,nCompY,dt);
	    }
	}

	if (level < parent->finestLevel())
	    getLevel(level+1).getViscFluxReg().CrseInitFinish();
    }
    //
    // Clean up memory
    //
    if (alpha)
	delete alpha;
    diffusion->removeFluxBoxesLevel(fluxSCn);
    diffusion->removeFluxBoxesLevel(fluxSCnp1);
    diffusion->removeFluxBoxesLevel(betan);
    diffusion->removeFluxBoxesLevel(betanp1);

    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_time = ParallelDescriptor::second() - strt_time;

    ParallelDescriptor::ReduceRealMax(run_time,IOProc);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "HeatTransfer::differential_spec_diffusion_update(): time: " << run_time << std::endl;
}

void
HeatTransfer::make_rho_prev_time ()
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::make_rho_prev_time()");

    const Real prev_time = state[State_Type].prevTime();
    //
    // We are assuming that the state is loaded with Density then species.
    //
    int ncomp = nspecies + 1;

    for (FillPatchIterator fpi(*this,*rho_ptime,1,prev_time,State_Type,Density,ncomp);
         fpi.isValid();
         ++fpi)
    {
        (*rho_ptime)[fpi.index()].copy(fpi());
    }
}

void
HeatTransfer::make_rho_curr_time ()
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::make_rho_curr_time()");

    const Real curr_time = state[State_Type].curTime();
    //
    // We are assuming that the state is loaded with Density then species.
    //
    int ncomp = nspecies + 1;

    for (FillPatchIterator fpi(*this,*rho_ctime,1,curr_time,State_Type,Density,ncomp);
         fpi.isValid();
         ++fpi)
    {
        (*rho_ctime)[fpi.index()].copy(fpi());
    }
}

void
HeatTransfer::adjust_spec_diffusion_update (MultiFab&              Phi_new,
					    const MultiFab*        Phi_old,
					    int                    sCompS,
					    Real                   dt,
					    Real                   time,
					    const Array<int>&      rho_flag,
					    const MultiFab*        rho_half,
					    int                    dataComp,
					    const MultiFab*        delta_rhs, 
					    const MultiFab*        alpha, 
					    const MultiFab* const* betanp1)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::adjust_spec_diffusion_update()");
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

    const int nGrowOp = 1;
    MultiFab rho_and_species(grids,nspecies+1,nGrowOp);
    //
    // Create and fill a full MultiFab of all species at this level and below.
    //
    FArrayBox tmp;

    for (FillPatchIterator fpi(*this,rho_and_species,nGrowOp,time,State_Type,Density,nspecies+1);
         fpi.isValid();
         ++fpi)
    {
        FArrayBox& rho_and_spec = rho_and_species[fpi];

        rho_and_spec.copy(fpi(),0,0,nspecies+1);

        tmp.resize(rho_and_spec.box(),1);
        tmp.copy(rho_and_spec,0,0,1);
        tmp.invert(1);

        for (int comp = 0; comp < nspecies; ++comp) 
            if (rho_flag[comp] == 2)
                rho_and_spec.mult(tmp,0,comp+1,1);
    }

    MultiFab rho_and_species_crse;

    if (level > 0) 
    {
        const int     nGrow   = 1;
        HeatTransfer& coarser = *(HeatTransfer*) &(parent->getLevel(level-1));

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

    const Real a = 1.0;
    const Real b = -dt;

    for (int comp = 0; comp < nspecies; ++comp)
    {
	const int state_ind = first_spec + comp;

        ViscBndry  visc_bndry;
        diffusion->getBndryDataGivenS(visc_bndry,rho_and_species,rho_and_species_crse,
                                      state_ind,comp+1,1,time,rho_flag[comp]);
        bool bndry_already_filled = true;

	Real           rhsscale;
        ABecLaplacian* visc_op;
	visc_op = diffusion->getViscOp(state_ind,a,b,time,visc_bndry,
                                       rho_half,rho_flag[comp],&rhsscale,
                                       dataComp+comp,betanp1,alpha,bndry_already_filled);
	visc_op->maxOrder(diffusion->maxOrder());

	rho_and_species.setBndry(bogus_value,comp+1,1); // Ensure computable corners

	visc_op->applyBC(rho_and_species,comp+1,1);

	delete visc_op;

	if (rho_flag[comp] == 2)
        {
	    for (MFIter Smfi(rho_and_species); Smfi.isValid(); ++Smfi)
            {
                FArrayBox& fab = rho_and_species[Smfi];
     		fab.mult(fab,fab.box(),0,comp+1,1);
            }
        }
    }
    //
    // "Repair" fluxes to ensure they sum to zero, use rho_and_species to pick dominant comp
    //  Note that rho_and_species is now rho and rho.Y, not rho and Y.
    //
    for (int d =0; d < BL_SPACEDIM; ++d)
    {
        for (MFIter mfi(rho_and_species); mfi.isValid(); ++mfi)
        {
            FArrayBox& state = rho_and_species[mfi];
            const Box& box   = mfi.validbox();

            FArrayBox& fab = (*flux[d])[mfi.index()];
            FORT_REPAIR_FLUX(box.loVect(), box.hiVect(),
                             fab.dataPtr(),  ARLIM(fab.loVect()),  ARLIM(fab.hiVect()),
                             state.dataPtr(1),ARLIM(state.loVect()),ARLIM(state.hiVect()),
                             &nspecies, &d);
        }
    }
    //
    // Reset Phi_new using "repaired" fluxes.  Here, we assume the following
    // arrangement of terms:
    //
    //     Phi_new = Phi_old + (dt/A)(RHS-Div(flux))
    //            A = rhonph.alpha for rho_flag == 1
    //            A = 1            for rho_flag == 2
    //     
    //
    FArrayBox update;

    for (MFIter mfi(Phi_new); mfi.isValid(); ++mfi)
    {
	int        iGrid = mfi.index();
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
	FORT_RECOMP_UPDATE(box.loVect(), box.hiVect(),
			   update.dataPtr(),
			   ARLIM(update.loVect()),       ARLIM(update.hiVect()),
			   (*flux[0])[iGrid].dataPtr(),
                           ARLIM((*flux[0])[iGrid].loVect()),
                           ARLIM((*flux[0])[iGrid].hiVect()),
			   (*flux[1])[iGrid].dataPtr(),
			   ARLIM((*flux[1])[iGrid].loVect()),
                           ARLIM((*flux[1])[iGrid].hiVect()),
#if BL_SPACEDIM == 3
			   (*flux[2])[iGrid].dataPtr(),
			   ARLIM((*flux[2])[iGrid].loVect()),
                           ARLIM((*flux[2])[iGrid].hiVect()),
#endif
			   volume[iGrid].dataPtr(),
			   ARLIM(volume[iGrid].loVect()),ARLIM(volume[iGrid].hiVect()),
			   &nspecies);
	
	if (delta_rhs)
	    update.plus((*delta_rhs)[iGrid],box,dataComp,0,nspecies);

	update.mult(dt,box,0,nspecies);

	if (alpha)
	    update.divide((*alpha)[iGrid],box,dataComp,0,nspecies);

        tmp.resize(box,1);
        tmp.copy((*rho_half)[iGrid],0,0,1);
        tmp.invert(1);

	for (int i = 0; i < nspecies; ++i)
	    if (rho_flag[i] == 1)
		update.mult(tmp,0,i,1);

	if (Phi_old)
	    update.plus((*Phi_old)[iGrid],box,sCompS,0,nspecies);

        Phi_new[mfi].copy(update,box,0,box,sCompS,nspecies);
    }
}

void
HeatTransfer::diffuse_scalar_setup (Real        dt,
                                    int         sigma,
                                    int*        rho_flag, 
                                    MultiFab*&  delta_rhs,
                                    MultiFab*&  alpha, 
                                    MultiFab**& betan,
                                    MultiFab**& betanp1)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::diffuse_scalar_setup()");
    //
    // Do setup for implicit c-n solve for an arbitrary scalar.
    //
    // Note: should be ok for variable c_p.
    //
    const Real prev_time = state[State_Type].prevTime();

    NavierStokes::diffuse_scalar_setup(dt, sigma, rho_flag, 
                                       delta_rhs, alpha, betan, betanp1);
    alpha     = 0;
    delta_rhs = 0;
    betan     = 0;
    betanp1   = 0;
   
    if (sigma == Temp)
    {
        (*rho_flag) = 1;
    }
    else if (sigma == RhoH || sigma >= first_spec && sigma <= last_spec)
    {
        (*rho_flag) = 2;
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

    diffusion->allocFluxBoxesLevel(betan);
    diffusion->allocFluxBoxesLevel(betanp1);
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
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::compute_OT_radloss()");
    //
    // Get optically thin radiation losses (+ve when energy LOST).
    //
    BL_ASSERT(do_OT_radiation || do_heat_sink);
    BL_ASSERT(nGrow <= dqrad.nGrow());
    BL_ASSERT(dqrad.boxArray() == grids);

    const Real* dx = geom.CellSize();

    Real p_amb, dpdt_factor;
    FORT_GETPAMB(&p_amb, &dpdt_factor);

    const Real Patm = p_amb / P1atm_MKS;
    const Real T_bg = 298.0;

    FArrayBox X,tmp,rad;

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

        X.resize(box,nspecies);
        getChemSolve().massFracToMoleFrac(X,S,box,first_spec,0);

        FArrayBox& dqrad_fab = dqrad[fpi];

        if (do_OT_radiation)
        {
            getChemSolve().getOTradLoss_TDF(dqrad_fab,S,X,Patm,T_bg,box,0,Temp,0);
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
            rad.resize(box,1);
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
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::diffuse_temp_setup()");
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
    //
    // alpha = c_p^(n+1/2)
    // Note: rho accounted for in diffusion box
    //
    BL_ASSERT(alpha==0);
    alpha = new MultiFab (grids,1,0);
    compute_cp(prev_time+dt/2,*alpha);
}

void
HeatTransfer::diffuse_cleanup (MultiFab*&  delta_rhs,
                               MultiFab**& betan,
                               MultiFab**& betanp1,
                               MultiFab*&  alpha)
{
    delete delta_rhs;
    delete alpha;
    alpha = delta_rhs = 0;

    diffusion->removeFluxBoxesLevel(betan);
    diffusion->removeFluxBoxesLevel(betanp1);
}

void
HeatTransfer::diffuse_cleanup (MultiFab*&  delta_rhs,
                               MultiFab**& betan,
                               MultiFab**& betanp1)
{
    MultiFab* alpha = 0;
    diffuse_cleanup(delta_rhs, betan, betanp1, alpha);
}

void
HeatTransfer::velocity_diffusion_update (Real dt)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::velocity_diffusion_update()");
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

        MultiFab *delta_rhs = 0, **betan = 0, **betanp1 = 0;

        diffuse_velocity_setup(dt, delta_rhs, betan, betanp1);

        diffusion->diffuse_velocity(dt,be_cn_theta,get_rho_half_time(),rho_flag,
                                    delta_rhs,betan,betanp1);

        diffuse_cleanup(delta_rhs, betan, betanp1);
    }
}
    
void
HeatTransfer::diffuse_velocity_setup (Real        dt,
                                      MultiFab*&  delta_rhs,
                                      MultiFab**& betan,
                                      MultiFab**& betanp1)
{
    //
    // Do setup for implicit c-n solve for velocity.
    //
    BL_ASSERT(delta_rhs==0);
    BL_ASSERT(betan==0);
    BL_ASSERT(betanp1==0);
    const Real time = state[State_Type].prevTime();
    //
    // Assume always variable viscosity.
    //
    diffusion->allocFluxBoxesLevel(betan);
    diffusion->allocFluxBoxesLevel(betanp1);
    
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
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::getViscTerms()");
    //
    // Load "viscous" terms, starting from component = 0.
    //
    const int  last_comp = src_comp + num_comp - 1;
    int        icomp     = src_comp; // This is the current related state comp.
    int        load_comp = 0;        // Comp for result of current calculation.
    MultiFab** vel_visc  = 0;        // Potentially reused, raise scope
    const int  nGrow     = visc_terms.nGrow();
    //
    // Get Div(tau) from the tensor operator, if velocity and have non-const viscosity
    //
    if (src_comp < BL_SPACEDIM)
    {
        if (src_comp != Xvel || num_comp < BL_SPACEDIM)
            BoxLib::Error("tensor v -> getViscTerms needs all v-components at once");

        diffusion->allocFluxBoxesLevel(vel_visc);
        getViscosity(vel_visc, time);
        diffusion->getTensorViscTerms(visc_terms,time,0,vel_visc);
        icomp = load_comp = BL_SPACEDIM;
    }
    //
    // For Le != 1, need to get visc terms in a block
    //
    if (!unity_Le)
    {
	const int non_spec_comps = std::max(0,first_spec-src_comp) + std::max(0,last_comp-last_spec);
	const bool has_spec = src_comp <= last_spec && last_comp >= first_spec;
	BL_ASSERT( !has_spec || (num_comp>=nspecies+non_spec_comps));
	if (has_spec)
	{
	    const int sComp = first_spec - src_comp;
	    compute_differential_diffusion_terms(visc_terms,sComp,time);
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
		getRhoHViscTerms(visc_terms,RhoH-load_comp,time);
	    }
	    else
	    {
		const int  rho_flag = Diffusion::set_rho_flag(diffusionType[icomp]);
		MultiFab** beta     = 0;

		if (icomp == Density)
                {
                    visc_terms.setVal(0.0,load_comp,1);
                }
                else
                {
                    //
                    // Assume always variable viscosity / diffusivity.
                    //
                    diffusion->allocFluxBoxesLevel(beta);
                    getDiffusivity(beta, time, icomp, 0, 1);

                    diffusion->getViscTerms(visc_terms,icomp-load_comp,
                                            icomp,time,rho_flag,0,beta);
                    
                    diffusion->removeFluxBoxesLevel(beta);
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
        //
	// Assume always using variable viscosity.
        //
	diffusion->compute_divmusi(time,vel_visc,divmusi); // pre-computed visc above
	divmusi.mult((-2./3.),0,BL_SPACEDIM,0);
        visc_terms.plus(divmusi,Xvel,BL_SPACEDIM,0);
    }
    //
    // Clean up your mess ...
    //
    if (vel_visc)
        diffusion->removeFluxBoxesLevel(vel_visc);
    //
    // Ensure consistent grow cells
    //    
    if (nGrow > 0)
    {
        for (MFIter mfi(visc_terms); mfi.isValid(); ++mfi)
        {
            FArrayBox& vt  = visc_terms[mfi];
            const Box& box = mfi.validbox();
            FORT_VISCEXTRAP(vt.dataPtr(),ARLIM(vt.loVect()),ARLIM(vt.hiVect()),
                            box.loVect(),box.hiVect(),&num_comp);
        }
        visc_terms.FillBoundary(0,num_comp);
        //
        // Note: this is a special periodic fill in that we want to
        // preserve the extrapolated grow values when periodic --
        // usually we preserve only valid data.  The scheme relies on
        // the fact that there is good data in the "non-periodic" grow cells.
        // ("good" data produced via VISCEXTRAP above)
        //
        geom.FillPeriodicBoundary(visc_terms,0,num_comp,true);
    }
}

void
HeatTransfer::compute_differential_diffusion_terms (MultiFab& visc_terms,
						    int       sComp,
                                                    Real      time)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::compute_differential_diffusion_terms()");

    if (hack_nospecdiff)
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "... HACK!!! zeroing spec diffusion terms " << std::endl;
        visc_terms.setVal(0.0,sComp,nspecies);
        return;            
    }
    //
    // NOTE: This routine does not fill grow cells
    //
    BL_ASSERT(visc_terms.boxArray() == grids);
    BL_ASSERT(visc_terms.nComp() >= sComp+nspecies);

    MultiFab& S = get_data(State_Type,time);
    const Real* dx = geom.CellSize();

    const Array<int> rho_flag(nspecies,2); // Hardwired, for now.
    const int nGrowOp = 1; // Required by the operator to compute first-cut fluxes
    MultiFab s_tmp(grids,1,nGrowOp);
    //
    // FIXME: Should have this in level data
    //
    MultiFab **beta;
    diffusion->allocFluxBoxesLevel(beta,0,nspecies);
    getDiffusivity(beta, time, first_spec, 0, nspecies);

    MultiFab **flux;
    diffusion->allocFluxBoxesLevel(flux,0,1);

    const TimeLevel whichTime = which_time(State_Type,time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
    BL_ASSERT(first_spec == Density+1);
    //
    // Create and fill a full MultiFab of all species at this level and below.
    //
    MultiFab rho_and_species(grids,nspecies+1,nGrowOp);

    FArrayBox tmp;

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

    for (int comp = 0; comp < nspecies; ++comp)
    {
        const int state_ind = first_spec + comp;
	ViscBndry visc_bndry;
	diffusion->getBndryDataGivenS(visc_bndry,rho_and_species,rho_and_species_crse,
                                      state_ind,comp+1,1,time,rho_flag[comp]);

	ABecLaplacian visc_op(visc_bndry,dx);
	visc_op.setScalars(0.0,1.0);
	
	for (int d = 0; d < BL_SPACEDIM; d++)
	{
	    MultiFab bcoeffs(area[d].boxArray(),1,0);
	    MultiFab::Copy(bcoeffs,area[d],0,0,1,0);
	    for (MFIter mfi(bcoeffs); mfi.isValid(); ++mfi)
		bcoeffs[mfi].mult((*beta[d])[mfi.index()],comp,0,1);
	    visc_op.bCoefficients(bcoeffs,d);
	}

	MultiFab::Copy(s_tmp,rho_and_species,comp+1,0,1,0);
	
        visc_op.compFlux(D_DECL(*flux[0],*flux[1],*flux[2]),s_tmp);
        MultiFab* const* fluxKeep = (whichTime == AmrOldTime) ? SpecDiffusionFluxn : SpecDiffusionFluxnp1;
	for (int d = 0; d < BL_SPACEDIM; ++d)
	    MultiFab::Copy(*fluxKeep[d],*flux[d],0,comp,1,0);
	spec_diffusion_flux_computed[comp] = HT_ExplicitDiffusion;
    }
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
    diffusion->removeFluxBoxesLevel(beta);
    diffusion->removeFluxBoxesLevel(flux);
}

void
HeatTransfer::getTempViscTerms (MultiFab& visc_terms,
                                int       src_comp, 
                                Real      time)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::getTempViscTerms()");
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

    MultiFab** beta;
    diffusion->allocFluxBoxesLevel(beta);
    //
    // + div lambda grad T
    //
    getDiffusivity(beta, time, Temp, 0, 1);

    diffusion->getViscTerms(visc_terms,src_comp,Temp,time,rho_flag,0,beta);
    //
    // (see note at top of file on sign convention at top of file)
    // - div q rad
    //
    if (do_OT_radiation || do_heat_sink)
    {
        MultiFab dqrad(grids,1,nGrow);
        compute_OT_radloss(time,nGrow,dqrad);
        for (MFIter mfi(visc_terms); mfi.isValid(); ++mfi)
        {
            const int i = mfi.index();
            visc_terms[i].minus(dqrad[i],grids[i],0,Temp-src_comp,1);
        }
    }

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
        getChemSolve().getCpmixGivenTY(cp,S[dvt_mfi],spec,box,
                                       Temp,yComp,sCompCp);
        
        visc_terms[dvt_mfi].plus(delta_visc_terms[dvt_mfi],box,0,Temp-src_comp,1);
        visc_terms[dvt_mfi].divide(cp,0,Temp-src_comp,1);
    }
    //
    // Clean up.
    //
    diffusion->removeFluxBoxesLevel(beta);
}

void
HeatTransfer::getRhoHViscTerms (MultiFab& visc_terms,
                                int       src_comp, 
                                Real      time)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::getRhoHViscTerms()");
    //
    // NOTE: This routine does not fill grow cells
    //
    BL_ASSERT(visc_terms.boxArray()==grids);

    MultiFab** beta;
    const int nGrow = 0;
    diffusion->allocFluxBoxesLevel(beta);
    
    const int rhoh_rho_flag = 2;
    getDiffusivity(beta, time, RhoH, 0, 1);

    diffusion->getViscTerms(visc_terms,src_comp,RhoH,time,rhoh_rho_flag,0,beta);

    if (do_OT_radiation || do_heat_sink)
    {
        MultiFab dqrad(grids,1,nGrow);
        compute_OT_radloss(time,nGrow,dqrad);
        for (MFIter mfi(visc_terms); mfi.isValid(); ++mfi)
        {
            const int i = mfi.index();
            visc_terms[i].minus(dqrad[i],grids[i],0,RhoH-src_comp,1);
        }
    }
    //
    // Clean up.
    //
    diffusion->removeFluxBoxesLevel(beta);
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
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::set_rho_to_species_sum()");

    const BoxArray& sgrids = S_in.boxArray();

    BL_ASSERT(sgrids == S_out.boxArray());

    const int ngrids       = sgrids.size();
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
            FabMinMax(S_in[i], box, 0.0, Real_MAX, s_first_spec, s_num_spec);
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
            S_out[i].plus(S_in[i],box,spec,s_density,1);
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
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::scale_species()");

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
                if (S[k].min(spec) < 0.0)
                {
                    Real* sdat = S[k].dataPtr(spec);

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
        const int k = mfi.index();

        for (int spec = s_first_spec; spec <= s_last_spec; spec++)
        {
            factor[k].plus(S[k],spec,0,1);
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
        Real  min_factor = factor[k].min();
        Real  max_factor = factor[k].max();
        Real  min_rho    = S[k].min(s_density);
        Real  max_rho    = S[k].max(s_density);

        if ((min_factor>0.0||max_factor<0.0) && (min_rho>0.0||max_rho<0.0))
        {
            //
            // Anything potentially worrisome is non-zero --> easy case.
            //
            factor[k].invert(1.0,0,1);
            factor[k].mult(S[k],s_density,0,1);
            for (int spec = s_first_spec; spec <= s_last_spec; spec++)
                S[k].mult(factor[k],0,spec,1);
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
            Real* sumdat = factor[k].dataPtr();
            Real* rhodat = S[k].dataPtr(s_density);
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
                        S[k].dataPtr(spec)[i] *= factori;
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
                        if (S[k].dataPtr(spec)[i] != 0.0)
                            nonzero++;
                    Real shift = (rhoi-sumi)/nonzero;
                    for (int spec = s_first_spec; spec <= s_last_spec; spec++) 
                        if (S[k].dataPtr(spec)[i] != 0.0) 
                            S[k].dataPtr(spec)[i] += shift; 
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
        BL_PROFILE(BL_PROFILE_THIS_NAME() + "::temperature_stats()");
        //
        // Calculate some minimums and maximums.
        //
        Real tmin, dmin, hmin, tmax, dmax, hmax;

        tmin = dmin = hmin =  1.0e30;
        tmax = dmax = hmax = -1.0e30;

        for (MFIter S_mfi(S); S_mfi.isValid(); ++S_mfi)
        {
            const Box& bx = S_mfi.validbox();

            tmin = std::min(tmin,S[S_mfi].min(bx,Temp));
            tmax = std::max(tmax,S[S_mfi].max(bx,Temp));
            dmin = std::min(dmin,S[S_mfi].min(bx,Density));
            dmax = std::max(dmax,S[S_mfi].max(bx,Density));
            hmin = std::min(hmin,S[S_mfi].min(bx,RhoH));
            hmax = std::max(hmax,S[S_mfi].max(bx,RhoH));
        }

        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceRealMin(tmin,IOProc);
        ParallelDescriptor::ReduceRealMin(dmin,IOProc);
        ParallelDescriptor::ReduceRealMin(hmin,IOProc);

        ParallelDescriptor::ReduceRealMax(tmax,IOProc);
        ParallelDescriptor::ReduceRealMax(dmax,IOProc);
        ParallelDescriptor::ReduceRealMax(hmax,IOProc);

        FArrayBox   Y;
        bool        aNegY = false;
        Array<Real> minY(nspecies,1.e30);
        FArrayBox   tmp;

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

        for (int i = 0; i < nspecies; ++i)
        {
            ParallelDescriptor::ReduceRealMin(minY[i],IOProc);

            if (minY[i] < 0) aNegY = true;
        }

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "  Min,max temp = " << tmin << ", " << tmax << '\n';
            std::cout << "  Min,max rho  = " << dmin << ", " << dmax << '\n';
            std::cout << "  Min,max rhoh = " << hmin << ", " << hmax << '\n';
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
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::compute_rhoRT()");

    BL_ASSERT(pComp<p.nComp());

    const BoxArray& sgrids = S.boxArray();
    int nCompY = last_spec - first_spec + 1;

    FArrayBox state,tmp;

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
	    state.copy((*temp)[i],box,0,box,sCompT,1);
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

	getChemSolve().getPGivenRTY(p[i],state,state,state,
				    box,sCompR,sCompT,sCompY,pComp);
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
HeatTransfer::advance_setup (Real time,
                             Real dt,
                             int  iteration,
                             int  ncycle)
{
    NavierStokes::advance_setup(time, dt, iteration, ncycle);
    //
    // Make sure the new state has values so that c-n works
    // in the predictor (of the predictor)--rbp.
    //
    for (int k = 0; k < num_state_type; k++)
    {
        MultiFab& nstate = get_new_data(k);
        MultiFab& ostate = get_old_data(k);

        MultiFab::Copy(nstate,ostate,0,0,nstate.nComp(),0);
    }

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

Real
HeatTransfer::advance (Real time,
                       Real dt,
                       int  iteration,
                       int  ncycle)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::advance(Real,...");

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
                  << " with dt = "         << dt << std::endl;
    }

    advance_setup(time,dt,iteration,ncycle);

    if (do_check_divudt)
        checkTimeStep(dt);
    
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& S_old = get_old_data(State_Type);
    //
    // Reset flag that fill patched state data is good
    //
    FillPatchedOldState_ok = true;
    //
    // Compute traced states for normal comp of velocity at half time level.
    //
    Real dt_test = 0.0, dummy = 0.0;    
    dt_test = predict_velocity(dt,dummy);
    
    if (ParallelDescriptor::IOProcessor())
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
        int havedivu = 1;
        mac_rhs = create_mac_rhs(time,dt);
        mac_project(time,dt,S_old,mac_rhs,havedivu);
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
    //
    // Save off a copy of the pre-chem state
    //
    const int fuelComp = getChemSolve().index(fuelName) + first_spec;
#ifdef BL_PLOT_CONSUMPTION
    if (plot_auxDiags)
        auxDiag->copy(S_old,fuelComp,0,1);
#endif
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
        {
            BoxList bl;

            for (int i = 0; i < S_old.size(); i++)
                bl.push_back(BoxLib::grow(S_old.boxArray()[i],aux_boundary_data_old.nGrow()));

            BoxArray ba(bl);
            //
            // This MF is guaranteed to cover tmpFABs & valid region of S_old.
            //
            MultiFab tmpS_old(ba,NUM_STATE,0);
            //
            // Note that S_old & tmpS_old have the same distribution.
            //
            const int ngrow = aux_boundary_data_old.nGrow();

            for (FillPatchIterator fpi(*this,S_old,ngrow,prev_time,State_Type,0,NUM_STATE);
                 fpi.isValid();
                 ++fpi)
            {
                tmpS_old[fpi.index()].copy(fpi());
            }

            tmpFABs.copy(tmpS_old);
        }

        strang_chem(S_old,  dt,HT_LeaveYdotAlone);
        strang_chem(tmpFABs,dt,HT_LeaveYdotAlone);

        aux_boundary_data_old.copyFrom(tmpFABs,BL_SPACEDIM,0,nComp);
    }
    //
    // Find change due to first Strang step.
    //
#ifdef BL_PLOT_CONSUMPTION
    if (plot_auxDiags)
    {
        MultiFab tmp(auxDiag->boxArray(),1,0);
        tmp.setVal(0);
        tmp.copy(S_old,fuelComp,0,1);
        for (MFIter mfi(*auxDiag); mfi.isValid(); ++mfi)
            (*auxDiag)[mfi].minus(tmp[mfi],0,0,1);
    }
#endif

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "HeatTransfer::advance(): after first Strang step\n";

    temperature_stats(S_old);
    //
    // Activate hook in FillPatch hack to get better data now.
    //
    FillPatchedOldState_ok = false;
    //
    // Compute tn coeffs based on chem-advance tn data
    //  (these are used in the Godunov extrapolation)
    //
    const int nScalDiffs = NUM_STATE-BL_SPACEDIM-1;
    calcDiffusivity(prev_time,dt,iteration,ncycle,Density+1,nScalDiffs);
    //
    // Godunov-extrapolate states to cell edges
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
    bool do_adv_reflux = true;
    if (RhoH > first_scalar)
	scalar_advection(dt,first_scalar,RhoH-1,do_adv_reflux);
    if (RhoH < last_scalar)
	scalar_advection(dt,RhoH+1,last_scalar,do_adv_reflux);
    //
    // Copy old-time boundary Density & RhoH into estimate for new-time RhoH.
    //
    aux_boundary_data_new.copy(aux_boundary_data_old,Density-BL_SPACEDIM,0,1);
    aux_boundary_data_new.copy(aux_boundary_data_old,RhoH-BL_SPACEDIM,   1,1);
    //
    // Save rho used in rho-states, needed for replacing with new one
    //  NOTE: WE LOAD/USE GROW CELLS HERE SO WE CAN FIX BOUNDARY DATA AS WELL
    //
    MultiFab Rho_hold(grids,1,LinOp_grow);

    BL_ASSERT(LinOp_grow == 1);

    for (MFIter mfi(*rho_ctime); mfi.isValid(); ++mfi)
    {
        const Box box = BoxLib::grow(mfi.validbox(),LinOp_grow);
        Rho_hold[mfi.index()].copy((*rho_ctime)[mfi],box,0,box,0,1);
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
    //
    // Compute the update to momentum
    //
    if (do_mom_diff == 1)
      momentum_advection(dt,do_reflux);
    //
    // Predictor
    //
    int corrector = 0;
    //
    // Set tnp1 coeffs to tn values in first round of predictor
    //
    BL_PROFILE_TIMER(ptimer, BL_PROFILE_THIS_NAME() + "::advance()" + "-predictor");
    BL_PROFILE_START(ptimer);

    MultiFab::Copy(*diffnp1_cc,*diffn_cc,0,0,nScalDiffs,diffn_cc->nGrow());
    temp_update(dt,corrector);         // Here, predict n+1 coeffs using n coeffs
    temperature_stats(S_new);

    calcDiffusivity(cur_time,dt,iteration,ncycle,Density+1,nScalDiffs);
    spec_update(time,dt,corrector);
    
    set_overdetermined_boundary_cells(time + dt); // RhoH BC's to see new Y's at n+1

    do_adv_reflux = false;
    scalar_advection(dt,RhoH,RhoH,do_adv_reflux); // Get aofs for RhoH now

    rhoh_update(time,dt,corrector);
    RhoH_to_Temp(S_new);
    temperature_stats(S_new);

    BL_PROFILE_STOP(ptimer);
    //
    // Corrector
    //
    BL_PROFILE_TIMER(ctimer, BL_PROFILE_THIS_NAME() + "::advance()" + "-corrector");
    BL_PROFILE_START(ctimer);

    corrector = 1;
    calcDiffusivity(cur_time,dt,iteration,ncycle,Density+1,nScalDiffs);
    tracer_update(dt,corrector);
    spec_update(time,dt,corrector);

    set_overdetermined_boundary_cells(time + dt);// RhoH BC's to see new Y's at n+1

    do_adv_reflux = true;
    scalar_advection(dt,RhoH,RhoH,do_adv_reflux); // Get aofs for RhoH now
    
    rhoh_update(time,dt,corrector);
    RhoH_to_Temp(S_new); 

    BL_PROFILE_STOP(ctimer);
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
#ifdef BL_PLOT_CONSUMPTION
    if (plot_auxDiags)
    {
        MultiFab tmp(auxDiag->boxArray(),1,0);
        tmp.setVal(0);
        tmp.copy(S_new,fuelComp,0,1);
        for (MFIter mfi(*auxDiag); mfi.isValid(); ++mfi)
            (*auxDiag)[mfi].plus(tmp[mfi],0,0,1);
    }
#endif

    strang_chem(S_new,dt,HT_EstimateYdotNew);

#ifdef BL_PLOT_CONSUMPTION
    if (plot_auxDiags)
    {
        MultiFab tmp(auxDiag->boxArray(),1,0);
        tmp.setVal(0);
        tmp.copy(S_new,fuelComp,0,1);
        for (MFIter mfi(*auxDiag); mfi.isValid(); ++mfi)
        {
            (*auxDiag)[mfi].minus(tmp[mfi],0,0,1);
            (*auxDiag)[mfi].mult(1.0/dt);
        }
    }
#endif
    //
    //  HACK!!  What are we really supposed to do here?
    //  Deactivate hook in FillPatch hack so that old data really is old data again
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
    calcDiffusivity(cur_time,dt,iteration,ncycle,Density+1,nScalDiffs,true);
    //
    // Set the dependent value of RhoRT to be the thermodynamic pressure.  By keeping this in
    //   the state, we can use the average down stuff to be sure that RhoRT_avg is avg(RhoRT),
    //   not ave(Rho)avg(R)avg(T), which seems to give the p-relax stuff in the mac Rhs troubles.
    //
    setThermoPress(cur_time);

    calc_divu(time+dt, dt, get_new_data(Divu_Type));

    if (!NavierStokes::initial_step && level != parent->finestLevel())
    {
        //
        // Set new divu to old div + dt*dsdt_old where covered by fine.
        //
        BoxArray crsndgrids = getLevel(level+1).grids;
        crsndgrids.coarsen(fine_ratio);
            
        MultiFab& divu_new = get_new_data(Divu_Type);
        MultiFab& divu_old = get_old_data(Divu_Type);
        MultiFab& dsdt_old = get_old_data(Dsdt_Type);
            
        for (MFIter mfi(divu_new); mfi.isValid();++mfi)
        {
            for (int i=0; i<crsndgrids.size(); ++i)
            {
                const Box ovlp = crsndgrids[i] & mfi.validbox();
                if (ovlp.ok())
                {
                    divu_new[mfi].copy(dsdt_old[mfi],ovlp,0,ovlp,0,1);
                    divu_new[mfi].mult(dt,ovlp,0,1);
                    divu_new[mfi].plus(divu_old[mfi],ovlp,0,0,1);
                }
            }
        }
    }
        
    calc_dsdt(time, dt, get_new_data(Dsdt_Type));
    if (NavierStokes::initial_step)
        MultiFab::Copy(get_old_data(Dsdt_Type),get_new_data(Dsdt_Type),0,0,1,0);
    //
    // Add the advective and other terms to get velocity (or momentum) at t^{n+1}.
    //
    velocity_update(dt);
    
    advance_cleanup(dt,iteration,ncycle);
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

        if (level > 0 && iteration == 1) p_avg->setVal(0);
    }
    //
    // Update estimate for allowable time step.
    //
    dt_test = std::min(dt_test, estTimeStep());
    if (ParallelDescriptor::IOProcessor())
        std::cout << "HeatTransfer::advance(): at end of time step\n";

    temperature_stats(S_new);
   
    return dt_test;
}

MultiFab*
HeatTransfer::create_mac_rhs (Real time, Real dt)
{
   MultiFab*    dsdt = getDsdt(0,time);
   MultiFab* mac_rhs = getDivCond(0,time);
   for (MFIter dsdtmfi(*dsdt); dsdtmfi.isValid(); ++dsdtmfi)
   {
       (*dsdt)[dsdtmfi].mult(.5*dt);
       (*mac_rhs)[dsdtmfi].plus((*dsdt)[dsdtmfi]);
   }

   delete dsdt;

   if (dt > 0.0) 
   {
       MultiFab  dpdt(grids,1,0);
       calc_dpdt(time,dt,dpdt,u_mac);

       for (MFIter mfi(dpdt); mfi.isValid(); ++mfi)
           (*mac_rhs)[mfi].plus(dpdt[mfi], grids[mfi.index()], 0,0,1);
   }

   return mac_rhs;
}

void
HeatTransfer::reset_rho_in_rho_states (const MultiFab& rho,
                                       Real            time,
                                       const int       sComp,
                                       const int       nComp)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::reset_rho_in_rho_states()");
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
    //
    // Now do the same for AuxBoundaryData.
    // This routine should only be called at new time.
    //
    const TimeLevel whichTime = which_time(State_Type,time);

    BL_ASSERT(whichTime == AmrNewTime);

    const int nGrow = LinOp_grow;

    BL_ASSERT(rho.nGrow() >= nGrow);
    //
    // Make a multifab of "rho" on the boxes of the boundary data.
    // Force them to have the same DistributonMapping().
    //
    MultiFab tmpRho;

    tmpRho.define(aux_boundary_data_new.equivBoxArray(),
                  1,
                  0,
                  aux_boundary_data_new.DistributionMap(),
                  Fab_allocate);
    BoxList bl;
    for (int i = 0; i < rho.size(); i++)
        bl.push_back(BoxLib::grow(rho.boxArray()[i],nGrow));
    //
    // This MF is guaranteed to cover tmpRho.
    //
    BoxArray bat(bl);
    MultiFab rhoGrow(bat,1,0);

    for (MFIter rmfi(rhoGrow); rmfi.isValid(); ++rmfi)
        rhoGrow[rmfi].copy(rho[rmfi],rho[rmfi].box());

    tmpRho.copy(rhoGrow);  // Parallel copy.

    for (MFIter rmfi(tmpRho); rmfi.isValid(); ++rmfi)
    {
        BL_ASSERT(rmfi.validbox() == aux_boundary_data_new[rmfi].box());

        tmp.resize(rmfi.validbox(),1);
        tmp.copy(tmpRho[rmfi],0,0,1);
        tmp.invert(1);

        BL_ASSERT(is_diffusive[RhoH]);
        BL_ASSERT(diffusionType[RhoH] == Laplacian_SoverRho);

        aux_boundary_data_new[rmfi].mult(tmp,0,1,1);
        aux_boundary_data_new[rmfi].mult(aux_boundary_data_new[rmfi],0,1,1);
    }
}

void
HeatTransfer::set_overdetermined_boundary_cells (Real time)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::set_overdetermined_boundary_cells()");

    const TimeLevel whichTime = which_time(State_Type,time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    AuxBoundaryData& rhoh_data = (whichTime == AmrOldTime) ? aux_boundary_data_old : aux_boundary_data_new;

    const int nGrow = (whichTime == AmrOldTime) ? HYP_GROW : LinOp_grow;
    //
    // Build a MultiFab parallel to State with appropriate # of ghost
    // cells built into the FABs themselves to cover rhoh_data.
    //
    BoxList bl;
    for (int i = 0; i < grids.size(); i++)
        bl.push_back(BoxLib::grow(grids[i],nGrow));

    BoxArray _bl(bl);
    MultiFab tmpS(_bl,1,0);

    BL_ASSERT(first_spec == Density+1);

    FArrayBox rhoh, tmp;

    const int sCompT = 0;
    const int sCompY = 1;
    const int sCompH = 0;
    //
    // A non-allocated MultiFab on grids for FPI below.
    //
    MultiFab t(grids,1,0,Fab_noallocate);

    for (FillPatchIterator T_fpi(*this,t,nGrow,time,State_Type,Temp,1),
             RhoY_fpi(*this,t,nGrow,time,State_Type,Density,nspecies+1);
         T_fpi.isValid() && RhoY_fpi.isValid();
         ++T_fpi, ++RhoY_fpi)
    {
        FArrayBox& RhoY = RhoY_fpi();

        BL_ASSERT(RhoY.box() == tmpS[RhoY_fpi.index()].box());

        tmp.resize(RhoY.box(),1);
        tmp.copy(RhoY,0,0,1);
        tmp.invert(1);

        const BoxArray& rhoh_BA = rhoh_data.equivBoxArray();

        for (int i = 0; i < rhoh_BA.size(); i++)
        {
            const Box isect = RhoY.box() & rhoh_BA[i];

            if (isect.ok())
            {
                for (int i = 1; i < nspecies+1; ++i)
                    RhoY.mult(tmp,isect,0,i,1);

                rhoh.resize(isect,1);
                getChemSolve().getHmixGivenTY(rhoh,T_fpi(),RhoY,isect,sCompT,sCompY,sCompH);
                rhoh.mult(RhoY,0,0,1);

                tmpS[T_fpi.index()].copy(rhoh);
            }
        }
    }

    const int RhoHcomp = (whichTime == AmrOldTime) ? RhoH-BL_SPACEDIM : 1;

    rhoh_data.copyFrom(tmpS,0,RhoHcomp,1); // Parallel copy.
}

DistributionMapping
HeatTransfer::getFuncCountDM (const BoxArray& bxba)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::getFuncCountDM()");

    DistributionMapping rr;
    rr.RoundRobinProcessorMap(bxba.size(),ParallelDescriptor::NProcs());

    MultiFab fctmpnew;
    fctmpnew.define(bxba, 1, 0, rr, Fab_allocate);
    fctmpnew.setVal(1);
    fctmpnew.copy(get_new_data(FuncCount_Type));  // Parallel copy.

    int count = 0;
    Array<long> vwrk(bxba.size());
    for (MFIter mfi(fctmpnew); mfi.isValid(); ++mfi)
        vwrk[count++] = static_cast<long>(fctmpnew[mfi].sum(0));

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

    MPI_Gatherv(vwrk.dataPtr(),
                count,
                ParallelDescriptor::Mpi_typemap<long>::type(),
                vwrktmp.dataPtr(),
                nmtags.dataPtr(),
                offset.dataPtr(),
                ParallelDescriptor::Mpi_typemap<long>::type(),
                IOProc,
                ParallelDescriptor::Communicator());

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
HeatTransfer::strang_chem (MultiFab&  state,
                           Real       dt,
                           YdotAction Ydot_action,
                           int        ngrow)
{
    //
    // Note: we only update the valid region of the "state" MultiFab.
    //
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::strang_chem(MultiFab&,...");

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
    //   state is in State_Type ordering, but starts at the scalars.
    //
    const int rho_comp  = Density; // state and State_Type completely aligned here
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
        Real p_amb, dpdt_factor;
        FORT_GETPAMB(&p_amb, &dpdt_factor);
        const Real Patm = p_amb / P1atm_MKS;

        //if (chem_integrator == ChemDriver::CKD_Vode)
        //{
            FArrayBox tmp;
            for (MFIter Smfi(state); Smfi.isValid(); ++Smfi)
            {

                tmp.resize(Smfi.validbox(),1);
                tmp.copy(state[Smfi],rho_comp,0,1);
                tmp.invert(1);

                for (int comp = 0; comp < nspecies; ++comp)
                    state[Smfi].mult(tmp,0,ycomp+comp,1);
            }
        //}

        if (ydot_tmp) 
            ydot_tmp->copy(state,ycomp,dCompYdot,nspecies);

        FArrayBox* chemDiag = 0;

        if (do_not_use_funccount)
        {
	    MultiFab tmp(state.boxArray(), 1, 0);

            for (MFIter Smfi(state); Smfi.isValid(); ++Smfi)
            {
                FArrayBox& fb = state[Smfi];
                const Box& bx = Smfi.validbox();
		FArrayBox& fc = tmp[Smfi];
#ifdef BL_PLOT_REACS
                if (plot_auxDiags && BoxLib::intersect(state.boxArray(),auxDiag->boxArray()).size() != 0)
                    chemDiag = &((*auxDiag)[Smfi]);
#endif
                getChemSolve().solveTransient(fb,fb,fb,fb,fc,bx,ycomp,Tcomp,0.5*dt,Patm,chem_integrator,chemDiag);
            }

	    get_new_data(FuncCount_Type).copy(tmp);
        }
        else
        {
            BoxList bl(state.boxArray());
            bl.minimize();
            bl.maxSize(max_grid_size_chem);

            if (bl.size() < 2*ParallelDescriptor::NProcs() && max_grid_size_chem >= 16)
                //
                // Let's chop the grids up a bit more.
                // We want to try and level out the chemistry work.
                //
                bl.maxSize(max_grid_size_chem/2);

            const BoxArray ba(bl);

            DistributionMapping dm = getFuncCountDM(ba);

            MultiFab tmp, fcnCntTemp;

            tmp.define(ba, state.nComp(), 0, dm, Fab_allocate);

            fcnCntTemp.define(ba, 1, 0, dm, Fab_allocate);

#ifdef BL_PLOT_REACS
            MultiFab diagTemp;
            const bool do_diag = plot_auxDiags && BoxLib::intersect(ba,auxDiag->boxArray()).size() != 0;
            if (do_diag)
            {
                diagTemp.define(ba, auxDiag->nComp(), 0, dm, Fab_allocate);
                diagTemp.copy(*auxDiag); // Parallel copy
            }
#endif
            if (verbose && ParallelDescriptor::IOProcessor())
                std::cout << "*** strang_chem: FABs in tmp MF: " << tmp.size() << std::endl;

            tmp.copy(state); // Parallel copy.

            for (MFIter Smfi(tmp); Smfi.isValid(); ++Smfi)
            {
                FArrayBox& fb = tmp[Smfi];
                const Box& bx = Smfi.validbox();
                FArrayBox& fc = fcnCntTemp[Smfi];
#ifdef BL_PLOT_REACS
                chemDiag = (do_diag ? &(diagTemp[Smfi]) : 0);
#endif
                getChemSolve().solveTransient(fb,fb,fb,fb,fc,bx,ycomp,Tcomp,0.5*dt,Patm,chem_integrator,chemDiag);
            }

            state.copy(tmp); // Parallel copy.

#ifdef BL_PLOT_REACS
            if (do_diag) auxDiag->copy(diagTemp); // Parallel copy
#endif
            get_new_data(FuncCount_Type).copy(fcnCntTemp); // Parallel copy.
        }

        if (ydot_tmp)
        {
            MultiFab state_tmp(ydot_tmp->boxArray(),nspecies,0);

            state_tmp.copy(state,ycomp,0,nspecies);
            
            for (MFIter Smfi(state_tmp); Smfi.isValid(); ++Smfi)
            {
                (*ydot_tmp)[Smfi].minus(state_tmp[Smfi]);
                (*ydot_tmp)[Smfi].mult(1/(0.5*dt),Smfi.validbox(),0,nspecies);
            }
        }

        //if (chem_integrator == ChemDriver::CKD_Vode)
            for (MFIter Smfi(state); Smfi.isValid(); ++Smfi)
                for (int comp = 0; comp < nspecies; ++comp)
                    state[Smfi].mult(state[Smfi],Smfi.validbox(),rho_comp,ycomp+comp,1);

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

    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_time = ParallelDescriptor::second() - strt_time;

    ParallelDescriptor::ReduceRealMax(run_time,IOProc);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "HeatTransfer::strang_chem time: " << run_time << std::endl;
}

static
BoxArray
GetBndryCells(const BoxArray& ba,
              int             n_grow,
              const Geometry& geom)
{
    const BoxList blgrids = BoxList(ba);

    BoxDomain bd;

    for (int i = 0; i < ba.size(); ++i)
    {
	BoxList gCells = BoxLib::boxDiff(BoxLib::grow(ba[i],n_grow),ba[i]);

	for (BoxList::iterator bli = gCells.begin(); bli != gCells.end(); ++bli)
	    bd.add(BoxLib::complementIn(*bli,blgrids));
    }

    BoxList bl;

    Array<IntVect> pshifts(27);

    for (BoxDomain::const_iterator bdi = bd.begin(); bdi != bd.end(); ++bdi)
    {
        bl.push_back(*bdi);
        //
        // Add in periodic ghost cells shifted to valid region.
        //
        geom.periodicShift(geom.Domain(), *bdi, pshifts);

        for (int i = 0; i < pshifts.size(); i++)
        {
            Box shftbox = *bdi + pshifts[i];

            bl.push_back(geom.Domain() & shftbox);
        }
    }

    bl.simplify();

    return BoxArray(bl);
}

void
HeatTransfer::compute_edge_states (Real               dt,
                                   std::vector<bool>* state_comps_to_compute)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::compute_edge_states()");
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
    std::vector<bool> do_predict(nState,true);
    if (do_mom_diff == 1) 
    {
      for (int d=0; d<BL_SPACEDIM; ++d)
        do_predict[Xvel+d] = true;
    }
    else
    {
      for (int d=0; d<BL_SPACEDIM; ++d)
        do_predict[Xvel+d] = false;
    } 
    for (int sigma=first_spec; sigma<=last_spec; ++sigma)
        do_predict[sigma] = false;
    if (do_set_rho_to_species_sum)
    {
        do_predict[Density] = false;
        do_predict[RhoH]    = false;
    }
    //
    // This logic and the associated array passed in allows the computation
    // of the edge states to be shut off for specific components.  This is
    // intended to allow special components, such as RhoK and RhoEps in
    // TurbHT be treated differently.  This logic tries to insure that
    // all components with interdependencies are turned on at the same time.
    //
    std::vector<bool> compute_comp(nState, true);
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
    if (use_forces_in_trans || (do_mom_diff == 1))
    {
        visc_terms.set(Xvel, new MultiFab(grids,BL_SPACEDIM,nGrowF));
        getViscTerms(visc_terms[Xvel],Xvel,BL_SPACEDIM,prev_time);
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
    FArrayBox edge[BL_SPACEDIM],tforces,tvelforces,Rho,U,spec,h,state,edgeVals;
    //
    // FillPatch'd state data.
    //
    MultiFab* divu_fp = NavierStokes::create_mac_rhs_grown(nGrowF,prev_time,dt);

    MultiFab Gp(grids,BL_SPACEDIM,1);

    getGradP(Gp, prev_pres_time);

    for (FillPatchIterator S_fpi(*this,*divu_fp,HYP_GROW,prev_time,State_Type,0,nState);
         S_fpi.isValid();
         ++S_fpi)
    {
        //
        // Gonna need this array on a per-grid basis
        //
        std::vector<bool> this_edge_state_computed(nState,false);

        const int i = S_fpi.index();

        Rho.resize(S_fpi().box(),1);
        Rho.copy(S_fpi(),Density,0,1);

        U.resize(S_fpi().box(),BL_SPACEDIM);
        U.copy(S_fpi(),Xvel,0,BL_SPACEDIM);
        //
        // Get the spec forces based on CC data (forces on EC data in getViscTerms)
        //
        if (use_forces_in_trans || (do_mom_diff == 1))
        {
            NavierStokes::getForce(tvelforces,i,nGrowF,Xvel,BL_SPACEDIM,Rho);
            godunov->Sum_tf_gp_visc(tvelforces,visc_terms[Xvel][i],Gp[i],Rho);
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

        FArrayBox vel;

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
                                     u_mac[0][i], edge[0],
                                     u_mac[1][i], edge[1],
#if (BL_SPACEDIM == 3)             
                                     u_mac[2][i], edge[2],
#endif
                                     U,vel,tvelforces,divu_dummy,
                                     comp,comp,bc.dataPtr(),
                                     iconserv_dummy,PRE_MAC);

                for (int d=0; d<BL_SPACEDIM; ++d)
                    (*EdgeState[d])[i].copy(edge[d],0,comp,1);

                this_edge_state_computed[comp] = true;
            }
        }
        //
        // Get spec edge states
        // FIXME: Fab copy reqd, force sum below pulls state and forces from same comp
        //
        if (compute_comp[first_spec])
        {
            spec.resize(S_fpi().box(),nspecies);
            spec.copy(S_fpi(),first_spec,0,nspecies);

            NavierStokes::getForce(tforces,i,nGrowF,first_spec,nspecies,Rho);

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
                                              visc_terms[first_spec][i], comp,
                                              (*divu_fp)[i], Rho, use_conserv_diff);
                    
                    int iconserv_dummy = 0;
                    godunov->edge_states(grids[i], dx, dt, velpred,
                                         u_mac[0][i], edge[0],
                                         u_mac[1][i], edge[1],
#if (BL_SPACEDIM==3)
                                         u_mac[2][i], edge[2],
#endif
                                         U,spec,tforces,(*divu_fp)[i],
                                         comp,state_ind,bc.dataPtr(),
                                         iconserv_dummy,PRE_MAC);
                }
                else
                {
                    FArrayBox junkDivu(tforces.box(),1);
                    junkDivu.setVal(0.);
                    godunov->Sum_tf_divu_visc(spec, tforces, comp, 1,
                                              visc_terms[first_spec][i], comp,
                                              junkDivu, Rho, use_conserv_diff);
                    
                    godunov->edge_states(grids[i], dx, dt, velpred,
                                         u_macG[0][i], edge[0],
                                         u_macG[1][i], edge[1],
#if (BL_SPACEDIM==3)
                                         u_macG[2][i], edge[2],
#endif
                                         U,spec,tforces,(*divu_fp)[i],
                                         comp,state_ind,bc.dataPtr(), 
                                         use_conserv_diff,FPU);
                }

                for (int d=0; d<BL_SPACEDIM; ++d)
                    (*EdgeState[d])[i].copy(edge[d],0,state_ind,1);

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
                    (*EdgeState[d])[i].setVal(0.0,edge[d].box(),Density,1);
                    for (int sigma=first_spec; sigma<=last_spec; ++sigma)
                        (*EdgeState[d])[i].plus((*EdgeState[d])[i],
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
                    (*EdgeState[d])[i].mult((*EdgeState[d])[i],(*EdgeState[d])[i].box(),
                                            (*EdgeState[d])[i].box(),Density,icomp,1);
        }

        if (compute_comp[Temp])
        {
            //
            // Get Temp edge states via extrap.
            // NOTE: This will provide a guess for H->T solve.
            //
            const int comp = 0;
            const int state_ind = Temp;
            int use_conserv_diff = 
                      (advectionType[state_ind] == Conservative) ? true : false;
            state.resize(S_fpi().box(),1);
            state.copy(S_fpi(),state_ind,0,1);

            NavierStokes::getForce(tforces,i,nGrowF,state_ind,1,Rho);
            Array<int> bc = getBCArray(State_Type,i,state_ind,1);

            AdvectionScheme adv_scheme = FPU;

            if (adv_scheme == PRE_MAC)
            {
                godunov->Sum_tf_divu_visc(state, tforces,  comp, 1,
                                          visc_terms[state_ind][i], 0,
                                          (*divu_fp)[i], Rho, use_conserv_diff);

                int iconserv_dummy = 0;
                godunov->edge_states(grids[i], dx, dt, velpred,
                                     u_mac[0][i], edge[0],
                                     u_mac[1][i], edge[1],
#if (BL_SPACEDIM==3)
                                     u_mac[2][i], edge[2],
#endif
                                     U, state, tforces, (*divu_fp)[i],
                                     comp, state_ind, bc.dataPtr(),
                                     iconserv_dummy, PRE_MAC);

            } else {

                FArrayBox junkDivu(tforces.box(),1);
                junkDivu.setVal(0.);
                godunov->Sum_tf_divu_visc(state, tforces,  comp, 1,
                                          visc_terms[state_ind][i], 0,
                                          junkDivu, Rho, use_conserv_diff);

                godunov->edge_states(grids[i], dx, dt, velpred,
                                     u_macG[0][i], edge[0],
                                     u_macG[1][i], edge[1],
#if (BL_SPACEDIM==3)
                                     u_macG[2][i], edge[2],
#endif
                                     U, state, tforces, (*divu_fp)[i],
                                     comp, state_ind, bc.dataPtr(), 
                                     use_conserv_diff, FPU);
            }

            for (int d=0; d<BL_SPACEDIM; ++d)
                (*EdgeState[d])[i].copy(edge[d],0,state_ind,1);

            this_edge_state_computed[state_ind] = true;
        }

        if (compute_comp[RhoH])
        {
            //
            // Set rhoh on edges = sum(rho.Y.H)
            //
            for (int d=0; d<BL_SPACEDIM; ++d)
            {
                (*EdgeState[d])[i].setVal(0.0,edge[d].box(),RhoH,1);
                h.resize(edge[d].box(),nspecies);
                getChemSolve().getHGivenT(h,(*EdgeState[d])[i],
                                          edge[d].box(),Temp,0);
                h.mult((*EdgeState[d])[i],edge[d].box(),first_spec,0,
                       nspecies);
                
                (*EdgeState[d])[i].setVal(0.0,edge[d].box(),RhoH,1);
                for (int comp=0; comp<nspecies; ++comp)
                    (*EdgeState[d])[i].plus(h,edge[d].box(),comp,RhoH,1);
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
                NavierStokes::getForce(tforces,i,nGrowF,state_ind,1,Rho);
                godunov->Sum_tf_divu_visc(state, tforces, comp, 1,
                                          visc_terms[state_ind][i], 0,
                                          (*divu_fp)[i], Rho,
                                          use_conserv_diff);
                Array<int> bc = getBCArray(State_Type,i,state_ind,1);
                int iconserv_dummy = 0;
                godunov->edge_states(grids[i], dx, dt, velpred,
                                     u_mac[0][i], edge[0],
                                     u_mac[1][i], edge[1],
#if (BL_SPACEDIM==3)
                                     u_mac[2][i], edge[2],
#endif
                                     U,state,tforces,(*divu_fp)[i],
                                     comp,state_ind,bc.dataPtr(),
                                     iconserv_dummy,PRE_MAC);

                for (int d=0; d<BL_SPACEDIM; ++d)
                    (*EdgeState[d])[i].copy(edge[d],0,state_ind,1);

                this_edge_state_computed[state_ind] = true;
            }
        }
    }

    delete divu_fp;
}

void
HeatTransfer::momentum_advection (Real dt, bool do_adv_reflux)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... advect momentum\n";

    const int   finest_level   = parent->finestLevel();
    Box edge_bx[BL_SPACEDIM];
    FArrayBox edge[BL_SPACEDIM], h[BL_SPACEDIM], flux;
    //
    // Compute the advective forcing for momentum.
    //
    int use_conserv_diff = true;
    for (MFIter AofS_mfi(*aofs); AofS_mfi.isValid(); ++AofS_mfi)
    {
        const int i = AofS_mfi.index();

        for (int d=0; d<BL_SPACEDIM; ++d)
        {
            edge_bx[d] = BoxLib::surroundingNodes(grids[i],d);
            edge[d].resize(edge_bx[d],1);
        }

        for (int comp = 0 ; comp < BL_SPACEDIM ; comp++ )
        {
            // 
            // If here, edge states at n+1/2 have already been computed, get a copy
            // 
            for (int d=0; d<BL_SPACEDIM; ++d)
                edge[d].copy((*EdgeState[d])[i],edge_bx[d],comp,edge_bx[d],0,1);
 
            godunov->ComputeAofs(grids[i],
                                 area[0][i],u_mac[0][i],edge[0],
                                 area[1][i],u_mac[1][i],edge[1],
#if BL_SPACEDIM==3
                                 area[2][i],u_mac[2][i],edge[2],
#endif
                                 volume[i],(*aofs)[i],comp,
                                 use_conserv_diff);
            //
            // Get fluxes for diagnostics and refluxing.
            //
            pullFluxes(i, comp, 1, edge[0], edge[1], edge[2], dt);
        }
    }
    //
    // pullFluxes() contains CrseInit() calls -- complete the process.
    //
    if (do_adv_reflux && level < finest_level)
        getAdvFluxReg(level+1).CrseInitFinish();
}

void
HeatTransfer::scalar_advection (Real dt,
                                int  fscalar,
                                int  lscalar,
                                bool do_adv_reflux)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::scalar_advection()");
    //
    // Compute the advection flux divergences
    //
    Box edge_bx[BL_SPACEDIM];
    FArrayBox edge[BL_SPACEDIM], h[BL_SPACEDIM], flux;
    const bool do_special_rhoh = nspecies>0 
        && do_set_rho_to_species_sum 
        && fscalar>=RhoH && lscalar<=RhoH
        && !unity_Le 
        && do_add_nonunityLe_corr_to_rhoh_adv_flux;
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
        const int rho_flag = 2; // FIXME: Messy assumption
        MultiFab *alpha=0;      //  -ditto-
        MultiFab **fluxSC, **fluxi, **rhoh_visc;
        diffusion->allocFluxBoxesLevel(fluxSC,0,1);
        diffusion->allocFluxBoxesLevel(fluxi,0,nspecies);
        diffusion->allocFluxBoxesLevel(fluxNULN,0,1);
        diffusion->allocFluxBoxesLevel(rhoh_visc,0,1);
        FArrayBox eTemp, h;
        const int dataComp = 0; // coeffs loaded into 0-comp for all species
        const int nGrow = 1; // Size to grow fil-patched fab for T below
        //
        // Initialize fluxNULN (NULN = non-unity Lewis number)
        //
        for (int d = 0; d < BL_SPACEDIM; ++d)
            fluxNULN[d]->setVal(0.0);
        //
        // Get the NULN flux contrib from n data
        //
        getDiffusivity(rhoh_visc, prev_time, RhoH, 0, 1);

        const MultiFab* Rh = get_rho_half_time();
        //
        // TODO -- try to invert this loop?
        //
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
                                           Rh,rho_flag,&rhsscale,dataComp,
                                           rhoh_visc,alpha);

            visc_op->maxOrder(diffusion->maxOrder());
            MultiFab::Copy(Soln,S_old,sigma,0,1,0);

            for (MFIter Smfi(Soln); Smfi.isValid(); ++Smfi)
                Soln[Smfi].divide(S_old[Smfi],Smfi.validbox(),Density,0,1);

            visc_op->compFlux(D_DECL(*fluxSC[0],*fluxSC[1],*fluxSC[2]),Soln);
            for (int d=0; d < BL_SPACEDIM; ++d)
                fluxSC[d]->mult(-b/geom.CellSize()[d]);
            //
            // Here, get fluxi = (lambda/cp - rho.D)Grad(Y)
            //                 = lambda/cp.Grad(Y) + SpecDiffFlux
            //
            for (int d = 0; d < BL_SPACEDIM; ++d)
            {
                for (MFIter SDF_mfi(*SpecDiffusionFluxn[d]);
                     SDF_mfi.isValid();
                     ++SDF_mfi)
                {
                    const Box& ebox    = SDF_mfi.validbox();
                    FArrayBox& SDF_fab = (*SpecDiffusionFluxn[d])[SDF_mfi];
                    (*fluxi[d])[SDF_mfi].copy(SDF_fab,ebox,comp,ebox,comp,1);
                    (*fluxi[d])[SDF_mfi].plus((*fluxSC[d])[SDF_mfi],ebox,0,comp,1);
                }
            }
            delete visc_op;
        }
        //
        // Multiply fluxi by h_i, and add to running total.
        //
        for (FillPatchIterator Told_fpi(*this,S_old,nGrow,prev_time,State_Type,Temp,1);
             Told_fpi.isValid();
             ++Told_fpi)
        {
            const int i    = Told_fpi.index();
            const Box& box = Told_fpi.validbox();

            for (int d = 0; d < BL_SPACEDIM; ++d)
            {
                const Box ebox = BoxLib::surroundingNodes(box,d);
                eTemp.resize(ebox,1);
                FPLoc bc_lo =
                    fpi_phys_loc(get_desc_lst()[State_Type].getBC(Temp).lo(d));
                FPLoc bc_hi = 
                    fpi_phys_loc(get_desc_lst()[State_Type].getBC(Temp).hi(d));
                
                center_to_edge_fancy(Told_fpi(),eTemp,
                                     BoxLib::grow(box,BoxLib::BASISV(d)),0,0,1,
                                     geom.Domain(),bc_lo,bc_hi);
                
                h.resize(ebox,nspecies);
                getChemSolve().getHGivenT(h,eTemp,ebox,0,0);
                (*fluxi[d])[i].mult(h,ebox,0,0,nspecies);
                
                for (int comp=0; comp<nspecies; ++comp)
                    (*fluxNULN[d])[i].plus((*fluxi[d])[i],ebox,comp,0,1);
            }
        }
        //
        // Get the Le!=1 flux contrib from n+1 data.
        //
        getDiffusivity(rhoh_visc, cur_time, RhoH, 0, 1);
        //
        // TODO -- try to invert this loop?
        //
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
                                           Rh,rho_flag,&rhsscale,dataComp,
                                           rhoh_visc,alpha);

            visc_op->maxOrder(diffusion->maxOrder());
            MultiFab::Copy(Soln,S_new,sigma,0,1,0);

            for (MFIter Smfi(Soln); Smfi.isValid(); ++Smfi)
                Soln[Smfi].divide(S_new[Smfi],Smfi.validbox(),Density,0,1);

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
                    (*fluxi[d])[SDF_mfi].copy(SDF_fab,ebox,comp,ebox,comp,1);
                    (*fluxi[d])[SDF_mfi].plus((*fluxSC[d])[SDF_mfi],ebox,0,comp,1);
                }
            }
            delete visc_op;
        }
        //
        // Multiply fluxi by h_i, and add to running total
        //
        for (FillPatchIterator Tnew_fpi(*this,S_new,nGrow,cur_time,State_Type,Temp,1);
             Tnew_fpi.isValid();
             ++Tnew_fpi)
        {
            const int i    = Tnew_fpi.index();
            const Box& box = Tnew_fpi.validbox();

            for (int d = 0; d < BL_SPACEDIM; ++d)
            {
                const Box ebox = BoxLib::surroundingNodes(box,d);
                eTemp.resize(ebox,1);
                FPLoc bc_lo =
                    fpi_phys_loc(get_desc_lst()[State_Type].getBC(Temp).lo(d));
                FPLoc bc_hi = 
                    fpi_phys_loc(get_desc_lst()[State_Type].getBC(Temp).hi(d));
                
                center_to_edge_fancy(Tnew_fpi(),eTemp,BoxLib::grow(box,BoxLib::BASISV(d)),
                                     0,0,1,geom.Domain(),bc_lo,bc_hi);
                
                h.resize(ebox,nspecies);
                getChemSolve().getHGivenT(h,eTemp,ebox,0,0);
                (*fluxi[d])[i].mult(h,ebox,0,0,nspecies);
                
                for (int comp = 0; comp < nspecies; ++comp)
                    (*fluxNULN[d])[i].plus((*fluxi[d])[i],ebox,comp,0,1);
            }
        }
        //
        // Clean up
        //
        diffusion->removeFluxBoxesLevel(fluxSC);
        diffusion->removeFluxBoxesLevel(fluxi);
        diffusion->removeFluxBoxesLevel(rhoh_visc);
    }
    
    for (MFIter AofS_mfi(*aofs); AofS_mfi.isValid(); ++AofS_mfi)
    {
        const int i = AofS_mfi.index();
        for (int d=0; d<BL_SPACEDIM; ++d)
        {
            edge_bx[d] = BoxLib::surroundingNodes(grids[i],d);
            edge[d].resize(edge_bx[d],1);
        }
        for (int sigma=fscalar; sigma<=lscalar; ++sigma)
        {
            // 
            // If here, edge states at n+1/2 have already been computed, get a copy
            // 
            for (int d=0; d<BL_SPACEDIM; ++d)
                edge[d].copy((*EdgeState[d])[i],edge_bx[d],sigma,edge_bx[d],0,1);
            
            int use_conserv_diff = 
                      (advectionType[sigma] == Conservative) ? true : false;

            godunov->ComputeAofs(grids[i],
                                 area[0][i],u_mac[0][i],edge[0],
                                 area[1][i],u_mac[1][i],edge[1],
#if BL_SPACEDIM==3
                                 area[2][i],u_mac[2][i],edge[2],
#endif
                                 volume[i],(*aofs)[i],sigma,
                                 use_conserv_diff);
            //
            // Add divergence of fluxNULN to aofs[RhoH], and increment advective
            //  going into flux registers
            //
            if (sigma==RhoH && do_special_rhoh)
            {
                if (do_adv_reflux)
                    for (int d=0; d<BL_SPACEDIM; ++d)
                        edge[d].plus((*fluxNULN[d])[i],edge_bx[d],0,0,1);
                
                FArrayBox& staten = (*aofs)[i];
                const FArrayBox& stateo = staten;
                const Box& box = AofS_mfi.validbox();
                const FArrayBox& vol = volume[i];
                const Real mult = 1.0; // no dt scaling of aofs, done in scl_adv_upd
                const int nComp = 1;
                
                FORT_INCRWEXTFLXDIV(box.loVect(), box.hiVect(),
                                    (*fluxNULN[0])[i].dataPtr(),
                                    ARLIM((*fluxNULN[0])[i].loVect()),
                                    ARLIM((*fluxNULN[0])[i].hiVect()),
                                    (*fluxNULN[1])[i].dataPtr(),
                                    ARLIM((*fluxNULN[1])[i].loVect()),
                                    ARLIM((*fluxNULN[1])[i].hiVect()),
#if BL_SPACEDIM == 3
                                    (*fluxNULN[2])[i].dataPtr(),
                                    ARLIM((*fluxNULN[2])[i].loVect()),
                                    ARLIM((*fluxNULN[2])[i].hiVect()),
#endif
                                    stateo.dataPtr(RhoH),
                                    ARLIM(stateo.loVect()), ARLIM(stateo.hiVect()),
                                    staten.dataPtr(RhoH),
                                    ARLIM(staten.loVect()), ARLIM(staten.hiVect()),
                                    vol.dataPtr(),
                                    ARLIM(vol.loVect()), ARLIM(vol.hiVect()),
                                    &nComp, &mult);
            }
            //
            // Do refluxing
            //
            if (do_adv_reflux)
                pullFluxes(i,sigma,1,edge[0],edge[1],edge[2],dt);
        }
    }
    //
    // pullFluxes() contains CrseInit() calls. Got to complete the process.
    //
    if (do_reflux && level < parent->finestLevel())
        getAdvFluxReg(level+1).CrseInitFinish();
    //
    // Cleanup memoroy , if necessary
    //
    if (do_special_rhoh)
        diffusion->removeFluxBoxesLevel(fluxNULN);
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
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::spec_update()");
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
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::compute_rhoDgradYgradH()");
    //
    // Get cell-centered h_i
    //
    MultiFab h(grids,nspecies,1);
    compute_h(time,h);
    //
    // Get edge-centered rho.D
    // (copy spec visc from internal database directly (before inflow "zeroing"))
    //
    MultiFab** beta;
    diffusion->allocFluxBoxesLevel(beta,0,nspecies);
    getDiffusivity(beta, time, first_spec, 0, nspecies);
    //
    // Get result, using cell-centered Y,h and edge-centered rhoD
    //
    FArrayBox rdgydgh_spec_i,tmp;

    const Real* dx = geom.CellSize();

    rdgydgh.setVal(0,0);
    //
    // nspecies = number of species, ncomp is one greater due to fillpatching
    // density and species together.
    //
    int nspecies = last_spec - first_spec + 1;
    int ncomp    = nspecies + 1;

    for (FillPatchIterator fpi(*this,h,1,time,State_Type,Density,ncomp);
          fpi.isValid();
          ++fpi)
    {
            const int  i               = fpi.index();
            const int* lo              = grids[i].loVect();
            const int* hi              = grids[i].hiVect();
            FArrayBox& rho_and_species = fpi();

            rdgydgh_spec_i.resize(grids[i],1);
            DEF_LIMITS(rdgydgh_spec_i,prod,prodlo,prodhi);

            const  int* speclo  = rho_and_species.loVect();
            const  int* spechi  = rho_and_species.hiVect();

            tmp.resize(rho_and_species.box(),1);
            tmp.copy(rho_and_species,0,0,1);
            tmp.invert(1);

            for (int spec = first_spec; spec <= last_spec; spec++) 
            {
                const int comp = spec - first_spec;

                DEF_CLIMITSCOMP(h[i],hdat,hlo,hhi,comp);
                DEF_CLIMITSCOMP((*beta[0])[i],betax,betaxlo,betaxhi,comp);
                DEF_CLIMITSCOMP((*beta[1])[i],betay,betaylo,betayhi,comp);
#if (BL_SPACEDIM==3)
                DEF_CLIMITSCOMP((*beta[2])[i],betaz,betazlo,betazhi,comp);
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

                rdgydgh[i].plus(rdgydgh_spec_i);
            }
    }
  
    diffusion->removeFluxBoxesLevel(beta);
}

//
// An enum to clean up mac_sync...questionable usefullness
//
enum SYNC_SCHEME {ReAdvect, UseEdgeState, Other};

void
HeatTransfer::mac_sync ()
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::mac_sync()");

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... mac_sync\n";

    const Real strt_time = ParallelDescriptor::second();

    int        sigma;
    const int  finest_level   = parent->finestLevel();
    const int  ngrids         = grids.size();
    const Real prev_time      = state[State_Type].prevTime();
    const Real cur_time       = state[State_Type].curTime();
    const Real prev_pres_time = state[Press_Type].prevTime();
    const Real dt             = parent->dtLevel(level);
    MultiFab*  DeltaSsync     = 0; // hold (Delta rho)*q for conserved quantities
    MultiFab*  Rh             = get_rho_half_time();

    sync_setup(DeltaSsync);
    //
    // Compute the correction velocity.
    //
    mac_projector->mac_sync_solve(level,dt,Rh,fine_ratio);
    //
    // Update coarse grid state by adding correction from mac_sync solve.
    //
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

        for (int comp=BL_SPACEDIM; comp<NUM_STATE; ++comp)
        {
            if (sync_scheme[comp]==UseEdgeState)
            {
                int s_ind = comp - BL_SPACEDIM;
                mac_projector->mac_sync_compute(level,Ssync,comp,s_ind,
                                                EdgeState,comp,Rh,
                                                (level>0 ? &getAdvFluxReg(level):0),
                                                advectionType,modify_reflux_normal_vel,dt);
            }
        }
        
        Ssync->mult(dt,Ssync->nGrow());
        //
        // For all conservative variables Q (other than density)
        // express Q as rho*q and increment sync by -(sync_for_rho)*q 
        //
        FArrayBox delta_ssync;
            
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
                        
                    delta_ssync.resize(grd,1);
                    delta_ssync.copy(S_new[i],grd,istate,grd,0,1);
                    delta_ssync.divide(S_new[i],grd,Density,0,1);
                    FArrayBox& s_sync = (*Ssync)[i];
                    delta_ssync.mult(s_sync,grd,Density-BL_SPACEDIM,0,1);
                    (*DeltaSsync)[i].copy(delta_ssync,grd,0,grd,iconserved,1);
                    s_sync.minus(delta_ssync,grd,0,istate-BL_SPACEDIM,1);
                }
            }
        }
        //
        // Now, increment density.
        //
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            const int i = mfi.index();
            S_new[i].plus((*Ssync)[i],grids[i],Density-BL_SPACEDIM,Density,1);
        }
        make_rho_curr_time();

        const int numscal = NUM_STATE - BL_SPACEDIM;
        //
        // Set do_diffuse_sync to 0 for debugging reasons only.
        //
        if (do_mom_diff == 1)
        {
            for (MFIter Vsyncmfi(*Vsync); Vsyncmfi.isValid(); ++Vsyncmfi)
            {
                const int i    = Vsyncmfi.index();
                const Box vbox = (*rho_ctime).box(i);

                D_TERM((*Vsync)[i].divide((*rho_ctime)[i],vbox,0,Xvel,1);,
                       (*Vsync)[i].divide((*rho_ctime)[i],vbox,0,Yvel,1);,
                       (*Vsync)[i].divide((*rho_ctime)[i],vbox,0,Zvel,1););
            }
        }

        if (do_diffuse_sync)
        {
	    MultiFab** beta;
	    diffusion->allocFluxBoxesLevel(beta);
            if (is_diffusive[Xvel])
            {
                int rho_flag = (do_mom_diff == 0) ? 1 : 3;

                getViscosity(beta, cur_time);

                diffusion->diffuse_Vsync(Vsync,dt,be_cn_theta,Rh,rho_flag,beta);
            }
	    
	    if (!unity_Le && nspecies>0 && do_add_nonunityLe_corr_to_rhoh_adv_flux)
	    {
		//
		// Diffuse the species syncs such that sum(SpecDiffSyncFluxes) = 0
		//
		differential_spec_diffuse_sync(dt);

                MultiFab Soln(grids,1,1);
                const Real cur_time  = state[State_Type].curTime();
                const Real a = 1.0;     // Passed around, but not used
                Real rhsscale;          //  -ditto-
                const int rho_flag = 2; // FIXME: Messy assumption
                MultiFab *alpha=0;      //  -ditto-
                MultiFab **fluxSC, **fluxNULN, **rhoh_visc;
                diffusion->allocFluxBoxesLevel(fluxSC,0,1);
                diffusion->allocFluxBoxesLevel(fluxNULN,0,nspecies);
                diffusion->allocFluxBoxesLevel(rhoh_visc,0,1);
                FArrayBox eTemp, h;
                const int dataComp = 0; // coeffs loaded into 0-comp for all species
                const int nGrow = 1; // Size to grow fil-patched fab for T below
                  
                getDiffusivity(rhoh_visc, cur_time, RhoH, 0, 1);
                //
                // TODO -- try to invert this loop?
                //
                for (int comp = 0; comp < nspecies; ++comp)
                {
                    const Real b     = be_cn_theta;
                    const int  sigma = first_spec + comp;
                    //
                    //  start by getting lambda/cp.Grad(Ysync)
                    //   (note: neg of usual diff flux)
                    //
                    ViscBndry      visc_bndry;
                    ABecLaplacian* visc_op;

                    visc_op = diffusion->getViscOp(sigma,a,b,cur_time,
                                                   visc_bndry,Rh,
                                                   rho_flag,&rhsscale,dataComp,
                                                   rhoh_visc,alpha);

                    visc_op->maxOrder(diffusion->maxOrder());

                    MultiFab::Copy(Soln,*Ssync,sigma-BL_SPACEDIM,0,1,0);

                    for (MFIter Smfi(Soln); Smfi.isValid(); ++Smfi)
                        Soln[Smfi].divide(S_new[Smfi],Smfi.validbox(),Density,0,1);

		    visc_op->compFlux(D_DECL(*fluxSC[0],*fluxSC[1],*fluxSC[2]),Soln);
                    for (int d = 0; d < BL_SPACEDIM; ++d)
                        fluxSC[d]->mult(-b/geom.CellSize()[d]);
                    //
                    // Here, get fluxNULN = (lambda/cp - rho.D)Grad(Ysync)
                    //                    = lambda/cp.Grad(Ysync) + SpecSyncDiffFlux
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
                            fluxNULN_fab.copy(SDF_fab,ebox,comp,ebox,comp,1);
                            fluxNULN_fab.plus(fluxSC_fab,ebox,0,comp,1);
                        }
                    }
                    delete visc_op;
                }
                //
                // Multiply fluxi by h_i (let FLXDIV routine below sum up the fluxes)
                //
                for (FillPatchIterator Tnew_fpi(*this,S_new,nGrow,cur_time,State_Type,Temp,1);
                     Tnew_fpi.isValid();
                     ++Tnew_fpi)
                {
                    const int i    = Tnew_fpi.index();
                    const Box& box = Tnew_fpi.validbox();

                    for (int d = 0; d < BL_SPACEDIM; ++d)
                    {
                        const Box ebox = BoxLib::surroundingNodes(box,d);
                        eTemp.resize(ebox,1);
                        FPLoc bc_lo =
                            fpi_phys_loc(get_desc_lst()[State_Type].getBC(Temp).lo(d));
                        FPLoc bc_hi = 
                            fpi_phys_loc(get_desc_lst()[State_Type].getBC(Temp).hi(d));
                        
                        center_to_edge_fancy(Tnew_fpi(),eTemp,BoxLib::grow(box,BoxLib::BASISV(d)),
                                             0,0,1,geom.Domain(),bc_lo,bc_hi);
                        
                        h.resize(ebox,nspecies);
                        getChemSolve().getHGivenT(h,eTemp,ebox,0,0);
                        (*fluxNULN[d])[i].mult(h,ebox,0,0,nspecies);
                    }
                }

                for (MFIter Ssync_mfi(*Ssync); Ssync_mfi.isValid(); ++Ssync_mfi)
                {
                    const int i      = Ssync_mfi.index();
                    FArrayBox& syncn = (*Ssync)[i];
                    const FArrayBox& synco = (*Ssync)[i];
                    const Box& box = Ssync_mfi.validbox();
                    const FArrayBox& vol = volume[i];
                    //
                    // Multiply by dt*dt, one to make it extensive, and one because
                    // Ssync multiplied above by dt, need same units here.
                    //
                    const Real mult = dt*dt;
                    const int sigmaRhoH = RhoH - BL_SPACEDIM; // RhoH comp in Ssync
		    
                    FORT_INCRWEXTFLXDIV(box.loVect(), box.hiVect(),
                                        (*fluxNULN[0])[i].dataPtr(),
                                        ARLIM((*fluxNULN[0])[i].loVect()),
                                        ARLIM((*fluxNULN[0])[i].hiVect()),
                                        (*fluxNULN[1])[i].dataPtr(),
                                        ARLIM((*fluxNULN[1])[i].loVect()),
                                        ARLIM((*fluxNULN[1])[i].hiVect()),
#if BL_SPACEDIM == 3
                                        (*fluxNULN[2])[i].dataPtr(),
                                        ARLIM((*fluxNULN[2])[i].loVect()),
                                        ARLIM((*fluxNULN[2])[i].hiVect()),
#endif
                                        synco.dataPtr(sigmaRhoH),
                                        ARLIM(synco.loVect()),
                                        ARLIM(synco.hiVect()),
                                        syncn.dataPtr(sigmaRhoH),
                                        ARLIM(syncn.loVect()),
                                        ARLIM(syncn.hiVect()),
                                        vol.dataPtr(),
                                        ARLIM(vol.loVect()),
                                        ARLIM(vol.hiVect()),
                                        &nspecies, &mult);
                }
                diffusion->removeFluxBoxesLevel(fluxSC);
                diffusion->removeFluxBoxesLevel(fluxNULN);
                diffusion->removeFluxBoxesLevel(rhoh_visc);
            }

	    MultiFab **flux;
            diffusion->allocFluxBoxesLevel(flux);

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
                    
		    diffusion->diffuse_Ssync(Ssync,sigma,dt,be_cn_theta,Rh,
                                             rho_flag,flux,0,beta,alpha);
		    if (do_viscsyncflux && level > 0)
		    {
			for (MFIter mfi(*Ssync); mfi.isValid(); ++mfi)
			{
			    const int i=mfi.index();
			    for (int d=0; d<BL_SPACEDIM; ++d)
                                getViscFluxReg().FineAdd((*flux[d])[i],d,i,
                                                         0,state_ind,1,dt);
			}
		    }
                }
            }
            diffusion->removeFluxBoxesLevel(flux);
	    diffusion->removeFluxBoxesLevel(beta);
        }
        //
        // For all conservative variables Q (other than density)
        // increment sync by (sync_for_rho)*q_presync.
        //
        for (MFIter mfi(*Ssync); mfi.isValid(); ++mfi)
        {
            const int i = mfi.index();

            int iconserved = -1;

            for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
            {
                if (istate != Density && advectionType[istate] == Conservative)
                {
                    iconserved++;

                    (*Ssync)[i].plus((*DeltaSsync)[i],grids[i],
                                     iconserved,istate-BL_SPACEDIM,1);
                }
            }
        }
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
                    S_new[i].plus((*Ssync)[i],grids[i],sigma,BL_SPACEDIM+sigma,1);
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
            MultiFab increment(fine_grids, numscal, S_new_lev.nGrow());
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

            SyncInterp(*Ssync, level, increment, lev, ratio, 
                       Density-BL_SPACEDIM, Density-BL_SPACEDIM, nComp, 1, mult, 
                       sync_bc.dataPtr(), which_interp, Density);

            if (have_trac)
                SyncInterp(*Ssync, level, increment, lev, ratio, 
                           Trac-BL_SPACEDIM, Trac-BL_SPACEDIM, 1, 1, mult, 
                           sync_bc.dataPtr());

            if (have_rhort)
                SyncInterp(*Ssync, level, increment, lev, ratio, 
                           RhoRT-BL_SPACEDIM, RhoRT-BL_SPACEDIM, 1, 1, mult, 
                           sync_bc.dataPtr());

            SyncInterp(*Ssync, level, increment, lev, ratio, 
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
                        increment[i].plus(increment[i],fine_grids[i],
					  istate-BL_SPACEDIM,Density-BL_SPACEDIM,1);
                    }
                }
            }

            for (MFIter mfi(increment); mfi.isValid(); ++mfi)
            {
                int i = mfi.index();
                S_new_lev[i].plus(increment[i],fine_grids[i],0,Density,numscal);
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
            const BoxArray& fgrids   = fine_lev.grids;
            MultiFab&       fvolume  = fine_lev.volume;
            
            HeatTransfer&   crse_lev = getLevel(lev);
            const BoxArray& cgrids   = crse_lev.grids;
            MultiFab&       cvolume  = crse_lev.volume;
            const IntVect&  fratio   = crse_lev.fine_ratio;
            
            MultiFab& S_crse = crse_lev.get_new_data(State_Type);
            MultiFab& S_fine = fine_lev.get_new_data(State_Type);

            const int pComp = (have_rhort ? RhoRT : Trac);
            crse_lev.NavierStokes::avgDown(cgrids,fgrids,
                                           S_crse,S_fine,
                                           cvolume,fvolume,
                                           lev,lev+1,pComp,1,fratio);
        }
    }

    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_time = ParallelDescriptor::second() - strt_time;

    ParallelDescriptor::ReduceRealMax(run_time,IOProc);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "HeatTransfer:mac_sync(): lev: "
                  << level
                  << ", time: " << run_time << std::endl;
    }

    sync_cleanup(DeltaSsync);
}

void
HeatTransfer::differential_spec_diffuse_sync(Real dt)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::differential_spec_diffuse_sync()");

    if (hack_nospecdiff)
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "... HACK!!! skipping spec sync diffusion " << std::endl;

        for (int d=0; d<BL_SPACEDIM; ++d)
            SpecDiffusionFluxnp1[d]->setVal(0.0,0,nspecies);

        for (int comp=0; comp<nspecies; ++comp)
            spec_diffusion_flux_computed[comp] = HT_SyncDiffusion;

        return;            
    }

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "...doing differential sync diffusion..." << std::endl;
    //
    // Do implicit c-n solve for each scalar...but dont reflux.
    // Save the fluxes, coeffs and source term, we need 'em for later
    // Actually, since Ssync multiplied by dt in mac_sync before
    // a call to this routine (I hope...), Ssync is in units of s,
    // not the "usual" ds/dt...lets convert it (divide by dt) so
    // we can use a generic flux adjustment function
    //
    const Real cur_time = state[State_Type].curTime();
    MultiFab **betanp1;
    diffusion->allocFluxBoxesLevel(betanp1,0,nspecies);
    getDiffusivity(betanp1, cur_time, first_spec, 0, nspecies);

    MultiFab Rhs(grids,nspecies,0);
    const int spec_Ssync_sComp = first_spec - BL_SPACEDIM;
    MultiFab::Copy(Rhs,*Ssync,spec_Ssync_sComp,0,nspecies,0);
    Rhs.mult(1.0/dt,0,nspecies,0); // Make Rhs in units of ds/dt again...
    //
    // Some standard settings
    //
    const Array<int> rho_flag(nspecies,2);
    const MultiFab* alpha = 0;
    MultiFab** fluxSC;
    diffusion->allocFluxBoxesLevel(fluxSC,0,1);

    const MultiFab* Rh = get_rho_half_time();

    for (int sigma = 0; sigma < nspecies; ++sigma)
    {
        //
        // Here, we use Ssync as a source in units of s, as expected by diffuse_Ssync
        // (i.e., ds/dt ~ d(Ssync)/dt, vs. ds/dt ~ Rhs in diffuse_scalar).  This was
        // apparently done to mimic diffuse_Vsync, which does the same, because the
        // diffused result is an acceleration, not a velocity, req'd by the projection.
        //
	const int ssync_ind = first_spec + sigma - Density;
	diffusion->diffuse_Ssync(Ssync,ssync_ind,dt,be_cn_theta,
				 Rh,rho_flag[sigma],fluxSC,
                                 sigma,betanp1,alpha);
	//
	// Pull fluxes into flux array
	//
	for (int d=0; d<BL_SPACEDIM; ++d)
	    MultiFab::Copy(*SpecDiffusionFluxnp1[d],*fluxSC[d],0,sigma,1,0);
	spec_diffusion_flux_computed[sigma] = HT_SyncDiffusion;
    }
    diffusion->removeFluxBoxesLevel(fluxSC);
    //
    // Modify update/fluxes to preserve flux sum = 0
    // (Be sure to pass the "normal" looking Rhs to this generic function)
    //
    const int sCompS = first_spec - BL_SPACEDIM;
    const MultiFab* old_sync = 0;
    const int dataComp = 0; 
    adjust_spec_diffusion_update(*Ssync,old_sync,sCompS,dt,cur_time,rho_flag,
                                 Rh,dataComp,&Rhs,alpha,betanp1);
    //
    // Clean up memory
    //
    diffusion->removeFluxBoxesLevel(betanp1);
    //
    // Do refluxing AFTER flux adjustment
    //
    if (do_reflux)
    {
	for (int d=0; d<BL_SPACEDIM; ++d)
	{
	    for (MFIter fmfi(*SpecDiffusionFluxnp1[d]); fmfi.isValid(); ++fmfi)
	    {
                FArrayBox& fmfi_fab = (*SpecDiffusionFluxnp1[d])[fmfi];
		if (level < parent->finestLevel())
		    getLevel(level+1).getViscFluxReg().CrseInit(fmfi_fab,
								fmfi_fab.box(),
								d,0,first_spec,
								nspecies,-dt);
		if (level > 0)
		    getViscFluxReg().FineAdd(fmfi_fab,d,fmfi.index(),
					     0,first_spec,nspecies,dt);
	    }
	}
	if (level < parent->finestLevel())
	    getLevel(level+1).getViscFluxReg().CrseInitFinish();
    }
    
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "...finished differential sync diffusion..." << std::endl;    
}

void
HeatTransfer::reflux ()
{
    if (level == parent->finestLevel())
        return;

    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::reflux()");

    BL_ASSERT(do_reflux);

    MultiFab& S_crse = get_new_data(State_Type);
    //
    // First do refluxing step.
    //
    FluxRegister& fr_adv  = getAdvFluxReg(level+1);
    FluxRegister& fr_visc = getViscFluxReg(level+1);
    const Real    dt_crse = parent->dtLevel(level);
    const Real    scale   = 1.0/dt_crse;
    //
    // It is important, for do_mom_diff == 0, to do the viscous
    //   refluxing first, since this will be divided by rho_half
    //   before the advective refluxing is added.  In the case of
    //   do_mom_diff == 1, both components of the refluxing will
    //   be divided by rho^(n+1) in NavierStokes::level_sync.
    //
    fr_visc.Reflux(*Vsync,volume,scale,0,0,BL_SPACEDIM,geom);
    //
    // Set do_reflux_visc to 0 for debugging reasons only.
    //
    if (do_reflux_visc)
        fr_visc.Reflux(*Ssync,volume,scale,BL_SPACEDIM,0,NUM_STATE-BL_SPACEDIM,geom);

    const MultiFab* Rh = get_rho_half_time();

    if (do_mom_diff == 0) 
    {
        for (MFIter mfi(*Vsync); mfi.isValid(); ++mfi)
        {
            const int i = mfi.index();

            D_TERM((*Vsync)[i].divide((*Rh)[i],grids[i],0,Xvel,1);,
                   (*Vsync)[i].divide((*Rh)[i],grids[i],0,Yvel,1);,
                   (*Vsync)[i].divide((*Rh)[i],grids[i],0,Zvel,1););
        }
    }

    FArrayBox tmp;

    for (MFIter mfi(*Ssync); mfi.isValid(); ++mfi)
    {
        const int i = mfi.index();

        tmp.resize(grids[i],1);
        tmp.copy((*Rh)[i],0,0,1);
        tmp.invert(1);

        for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
        {
            if (advectionType[istate] == NonConservative)
            {
                const int sigma = istate -  BL_SPACEDIM;

                (*Ssync)[i].mult(tmp,0,sigma,1);
            }
        }
    }

    fr_adv.Reflux(*Vsync,volume,scale,0,0,BL_SPACEDIM,geom);
    fr_adv.Reflux(*Ssync,volume,scale,BL_SPACEDIM,0,NUM_STATE-BL_SPACEDIM,geom);

    const BoxArray& fine_boxes = getLevel(level+1).boxArray();
    const int nfine            = fine_boxes.size();
    //
    // This is necessary in order to zero out the contribution to any
    // coarse grid cells which underlie fine grid cells.
    //
    for (int kf = 0; kf < nfine; kf++)
    {
        Box bf = BoxLib::coarsen(fine_boxes[kf],fine_ratio);

        for (MFIter mfi(*Vsync); mfi.isValid(); ++mfi)
        {
            const int k  = mfi.index();
            const Box bx = bf & grids[k];

            if (bx.ok())
            {
                (*Vsync)[k].setVal(0,bx,0,BL_SPACEDIM);
                (*Ssync)[k].setVal(0,bx,0,NUM_STATE-BL_SPACEDIM);
            }
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
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::set_preferred_boundary_values()");
    //
    // Only do copy if request contains affected states,
    // and fillpatched data known to be no good.
    //
    if (state_index == State_Type)
    {
        const TimeLevel whichTime = which_time(State_Type,time);
        //
        // To get chem-advanced data instead of FP'd data at old time
        //
        if (!FillPatchedOldState_ok && whichTime == AmrOldTime && src_comp > BL_SPACEDIM)
        {
            aux_boundary_data_old.copyTo(S,src_comp-BL_SPACEDIM,dst_comp,num_comp);
        }
        //
        // To get RhoH computed with current T, Y data instead of FP'd RhoH
        // Note: "advance" is to make sure that RhoH is correct as needed.
        //
        if (src_comp <= RhoH && src_comp + num_comp > RhoH)
        {
            const AuxBoundaryData& data = (whichTime == AmrOldTime) ? aux_boundary_data_old : aux_boundary_data_new;

            const int RhoHcomp = (whichTime == AmrOldTime) ? RhoH-BL_SPACEDIM : 1;

            data.copyTo(S,RhoHcomp,RhoH-src_comp+dst_comp,1);
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
HeatTransfer::calcDiffusivity (const Real time,
			       const Real dt,
			       const int  iteration,
			       const int  ncycle,
			       const int  src_comp,
			       const int  num_comp)
{
    calcDiffusivity(time,dt,iteration,ncycle,src_comp,num_comp,false);
}

void
HeatTransfer::calcDiffusivity (const Real time,
			       const Real dt,
			       const int  iteration,
			       const int  ncycle,
			       const int  src_comp,
			       const int  num_comp,
                               bool       do_VelVisc)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::calcDiffusivity()");

    const TimeLevel whichTime = which_time(State_Type, time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    const int  nGrow           = 1;
    const int  offset          = BL_SPACEDIM + 1; // No diffusion coeff for vels or rho
    const int  last_comp       = src_comp + num_comp - 1;
    const bool has_spec        = src_comp < last_spec && last_comp > first_spec;
    const int  non_spec_comps  = std::max(0,first_spec-src_comp) + std::max(0,last_comp-last_spec);
    MultiFab*  temp            = new MultiFab(grids,1,nGrow);
    MultiFab*  rhospec         = new MultiFab(grids,nspecies+1,nGrow);
    MultiFab*  visc            = (whichTime == AmrOldTime) ? diffn_cc : diffnp1_cc;

    BL_ASSERT( !has_spec || (num_comp >= nspecies+non_spec_comps) );

    FArrayBox tmp, bcen;

    for (FillPatchIterator Rho_and_spec_fpi(*this,*rhospec,nGrow,time,State_Type,Density,nspecies+1),
             Temp_fpi(*this,*rhospec,nGrow,time,State_Type,Temp,nGrow);
         Rho_and_spec_fpi.isValid() && Temp_fpi.isValid();
         ++Rho_and_spec_fpi, ++Temp_fpi)
    {
        const int idx = Rho_and_spec_fpi.index();
        const Box box = BoxLib::grow(grids[idx],nGrow);
        //
        // Convert from RhoY_l to Y_l
        //
        tmp.resize(box,1);
        tmp.copy(Rho_and_spec_fpi(),0,0,1);
        tmp.invert(1);

        for (int n = 1; n < nspecies+1; n++)
            Rho_and_spec_fpi().mult(tmp,0,n,1);

        (*rhospec)[idx].copy(Rho_and_spec_fpi(),0,0,nspecies+1);

        (*temp)[idx].copy(Temp_fpi(),0,0,1);
    }

    BL_ASSERT(nspecies > 0 && has_spec);

    if (nspecies > 0 && has_spec)
    {
        BL_ASSERT(first_spec == Density+1);

        for (MFIter rho_and_spec(*rhospec); rho_and_spec.isValid(); ++rho_and_spec)
        {
            const int  idx     = rho_and_spec.index();
            const Box  box     = BoxLib::grow(grids[idx],nGrow);
            const int  vflag   = do_VelVisc;
            FArrayBox& tempFab = (*temp)[idx];

            BL_ASSERT(box == (*rhospec)[idx].box());

            const int nc_bcen = nspecies+2;
            int       dotemp  = 1;
            bcen.resize(box,nc_bcen);

            FORT_SPECTEMPVISC(box.loVect(),box.hiVect(),
                              ARLIM(tempFab.loVect()),ARLIM(tempFab.hiVect()),
                              tempFab.dataPtr(),
                              ARLIM((*rhospec)[idx].loVect()),ARLIM((*rhospec)[idx].hiVect()),
                              (*rhospec)[idx].dataPtr(1),
                              ARLIM(bcen.loVect()),ARLIM(bcen.hiVect()),bcen.dataPtr(),
                              &nc_bcen, &P1atm_MKS, &dotemp, &vflag);

            (*visc)[idx].copy(bcen,0,first_spec-offset,nspecies);
            (*visc)[idx].copy(bcen,nspecies,Temp-offset,1);

            if (do_VelVisc)
            {
                MultiFab* beta = (whichTime==AmrOldTime) ? viscn_cc : viscnp1_cc;

                (*beta)[idx].copy(bcen,nspecies+1,0,1);
            }
        }
    }
    //
    // Now get the rest.
    //
    for (int icomp = src_comp; icomp <= last_comp; icomp++)
    {
        const bool is_spec = icomp >= first_spec && icomp <= last_spec;

        if (!is_spec)
        {
            if (icomp == RhoH)
            {
                //
                // We assume that Temp has already been done.
                //
                BL_ASSERT(first_spec == Density+1);

                FArrayBox cpmix;

                for (MFIter rho_and_species(*rhospec); rho_and_species.isValid(); ++rho_and_species)
                {
                    const int idx = rho_and_species.index();
                    const Box box = BoxLib::grow(grids[idx],nGrow);
                    (*visc)[idx].copy((*visc)[idx],Temp-offset,RhoH-offset,1);
                    cpmix.resize(box,1);
                    const int sCompT = 0, sCompY = 1, sCompCp = 0;
                    getChemSolve().getCpmixGivenTY(cpmix,(*temp)[idx],(*rhospec)[idx],box,sCompT,sCompY,sCompCp);
                    (*visc)[idx].divide(cpmix,0,RhoH-offset,1);
                }
            }
            else if (icomp == Trac || icomp == RhoRT)
            {
                visc->setVal(trac_diff_coef, icomp-offset, 1, nGrow);
            }
        }
    }

    delete rhospec;
    delete temp;
}

void
HeatTransfer::getViscosity (MultiFab*  beta[BL_SPACEDIM],
                            const Real time)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::getViscosity()");

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

            center_to_edge_fancy((*visc)[viscMfi],(*beta[dir])[i],
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
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::getDiffusivity()");

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

            center_to_edge_fancy((*diff)[diffMfi],(*beta[dir])[i],
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
            const Box ebox      = BoxLib::surroundingNodes(mfi.validbox(),dir);
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
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::compute_vel_visc()");

    const int nGrow = beta->nGrow();

    BL_ASSERT(nGrow == 1);
    BL_ASSERT(first_spec == Density+1);

    MultiFab dummy(grids,1,0,Fab_noallocate);
    FArrayBox bcen,tmp;

    for (FillPatchIterator Temp_fpi(*this,dummy,nGrow,time,State_Type,Temp,1),
             Rho_and_spec_fpi(*this,dummy,nGrow,time,State_Type,Density,nspecies+1);
         Rho_and_spec_fpi.isValid() && Temp_fpi.isValid();
         ++Rho_and_spec_fpi, ++Temp_fpi)
    {
        const int  i            = Rho_and_spec_fpi.index();
        Box        box          = BoxLib::grow(grids[i],nGrow);
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

        bcen.resize(box,1);

        FORT_VELVISC(box.loVect(),box.hiVect(),
                     ARLIM(temp.loVect()),ARLIM(temp.hiVect()),temp.dataPtr(),
                     ARLIM(rho_and_spec.loVect()),ARLIM(rho_and_spec.hiVect()),
                     rho_and_spec.dataPtr(1),
                     ARLIM(bcen.loVect()),ARLIM(bcen.hiVect()),bcen.dataPtr());

        (*beta)[i].copy(bcen,0,0,1);
    }
}

void
HeatTransfer::calc_divu (Real      time,
                         Real      dt,
                         MultiFab& divu)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::calc_divu()");
    //
    // Get Mwmix, cpmix and pressure
    //
    const int nGrow = 0;

    int sCompR, sCompT, sCompY, sCompP, sCompCp, sCompMw;
    sCompR=sCompT=sCompY=sCompP=sCompCp=sCompMw=0;	

    MultiFab rho(grids,1,nGrow);
    MultiFab temp(grids,1,nGrow);

    const MultiFab& Rho_time = get_rho(time);
    MFIter          Rho_mfi(Rho_time);

    for (FillPatchIterator Temp_fpi(*this,temp,nGrow,time,State_Type,Temp,1);
         Rho_mfi.isValid() && Temp_fpi.isValid();
         ++Rho_mfi, ++Temp_fpi)
    {
        const int i = Rho_mfi.index();

        rho[i].copy(Rho_time[i],0,sCompR,1);
        temp[i].copy(Temp_fpi(),0,sCompT,1);
    }
    //
    // Note that state contains rho*species, so divide species by rho.
    //
    MultiFab species;

    species.define(grids,nspecies,nGrow,Fab_allocate);

    FArrayBox tmp;

    for (FillPatchIterator Spec_fpi(*this,species,nGrow,time,State_Type,first_spec,nspecies);
         Spec_fpi.isValid();
         ++Spec_fpi)
    {
        const int i = Spec_fpi.index();

        species[i].copy(Spec_fpi(),0,sCompY,nspecies);

        tmp.resize(grids[i],1);
        tmp.copy(rho[i],sCompR,0,1);
        tmp.invert(1);

        for (int ispecies = 0; ispecies < nspecies; ispecies++)
            species[i].mult(tmp,0,ispecies,1);
    }

    MultiFab mwmix(grids,1,nGrow);
    MultiFab cp(grids,1,nGrow);
    MultiFab p(grids,1,nGrow);
	    
    for (MFIter Rho_mfi(rho); Rho_mfi.isValid(); ++Rho_mfi)
    {
        const int  iGrid = Rho_mfi.index();
        const Box& box   = Rho_mfi.validbox();

        BL_ASSERT(box == grids[iGrid]);

        getChemSolve().getMwmixGivenY(mwmix[iGrid],species[iGrid],
                                      box,sCompY,sCompMw);
        getChemSolve().getCpmixGivenTY(cp[iGrid],temp[iGrid],species[iGrid],
                                       box,sCompT,sCompY,sCompCp);
        getChemSolve().getPGivenRTY(p[iGrid],
                                    rho[iGrid],temp[iGrid],species[iGrid],
                                    box,sCompR,sCompT,sCompY,sCompP);
    }
    //
    // divu = 1/T DT/dt + W_mix * sum (1/W)DY/Dt
    //
    MultiFab visc_terms(grids,1,1);
    MultiFab delta_divu(grids,1,nGrow);
    //
    // Compute rho*DT/Dt as
    //
    //   1/c_p (div lambda grad T + sum_l rho D grad h_l dot grad Y_l)
    //
    getViscTerms(visc_terms,Temp,1,time);
    MultiFab::Copy(divu,visc_terms,0,0,1,nGrow);

    for (MFIter Divu_mfi(divu); Divu_mfi.isValid(); ++Divu_mfi)
    {
        const int iGrid = Divu_mfi.index();
        divu[iGrid].divide(rho[iGrid],grids[iGrid],sCompR,0,1);
        divu[iGrid].divide(temp[iGrid],grids[iGrid],sCompT,0,1);
    }

    delta_divu.setVal(0.0);

    const Array<Real> mwt = getChemSolve().speciesMolecWt();
    MultiFab spec_visc_terms(grids,nspecies,0);
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
	    
    for (MFIter Divu_mfi(divu); Divu_mfi.isValid(); ++Divu_mfi)
    {
        const int  iGrid = Divu_mfi.index();
        const Box& box   = Divu_mfi.validbox();
        delta_divu[iGrid].divide(rho[iGrid],box,sCompR,0,1);
        delta_divu[iGrid].mult(mwmix[iGrid],box,0,0,1);
        divu[iGrid].plus(delta_divu[iGrid],box,0,0,1);
    }

    if (dt > 0.0)
    {
        //
        // Increment divu by
        //    sum_l (h_l/(c_p*T) - mw_mix/mw_l)*delta Y_l/dt 
        // (i.e., Y_l"dot")
        //
        const Array<Real> mwt = getChemSolve().speciesMolecWt();
        MultiFab h(grids,nspecies,nGrow);
        const int sCompH = 0;
        for (MFIter mfi(h); mfi.isValid(); ++mfi)
        {
            getChemSolve().getHGivenT(h[mfi],temp[mfi],mfi.validbox(),sCompT,sCompH);
        }

        for (FillPatchIterator Ydot_fpi(*this,delta_divu,0,time,Ydot_Type,0,nspecies);
             Ydot_fpi.isValid();
             ++Ydot_fpi)
        {
            for (int istate = first_spec; istate <= last_spec; istate++)
            {
                const int ispec = istate-first_spec;
                const int i     = Ydot_fpi.index();

                delta_divu[i].copy(h[i],ispec,0,1);
                delta_divu[i].divide(cp[i]);
                delta_divu[i].divide(temp[i]);
                delta_divu[i].mult(Ydot_fpi(),ispec,0,1);
                divu[i].plus(delta_divu[i]);

                delta_divu[i].copy(mwmix[i],0,0,1);
                delta_divu[i].divide(mwt[ispec]);
                delta_divu[i].mult(Ydot_fpi(),ispec,0,1);
                divu[i].minus(delta_divu[i]);
            }
        }
    }
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
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::calc_dpdt()");

    Real dt = crse_dt;

    Real p_amb, dpdt_factor;
    
    FORT_GETPAMB(&p_amb, &dpdt_factor);
    
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
        temp[Temp_fpi.index()].copy(Temp_fpi(),0,sCompT,1);
    }

    for (FillPatchIterator rhoRT_fpi(*this,rhoRT,nGrow,time,State_Type,pComp,1);
         rhoRT_fpi.isValid();
         ++rhoRT_fpi)
    {
        rhoRT[rhoRT_fpi.index()].copy(rhoRT_fpi(),0,0,1);
    }
    //
    // Note that state contains rho*species, so divide species by rho.
    //
    MultiFab species;

    species.define(grids,nspecies,nGrow,Fab_allocate);

    FArrayBox tmp;

    for (FillPatchIterator Spec_fpi(*this,species,nGrow,time,State_Type,first_spec,nspecies);
         Spec_fpi.isValid();
         ++Spec_fpi)
    {
        const int i  = Spec_fpi.index();
        const Box bx = BoxLib::grow(grids[i],nGrow);

        species[i].copy(Spec_fpi(),0,sCompY,nspecies);

        tmp.resize(bx,1);
        tmp.copy(rho[i],sCompR,0,1);
        tmp.invert(1);

        for (int ispecies = 0; ispecies < nspecies; ispecies++)
            species[i].mult(tmp,0,ispecies,1);
    }
    
    MultiFab mwmix(grids,1,nGrow);
    MultiFab cp(grids,1,nGrow);
    
    for (MFIter Rho_mfi(rho); Rho_mfi.isValid(); ++Rho_mfi)
    {
        const int idx = Rho_mfi.index();
        const Box box = BoxLib::grow(Rho_mfi.validbox(),nGrow);
        
        BL_ASSERT(Rho_mfi.validbox() == grids[idx]);
        
        getChemSolve().getMwmixGivenY(mwmix[idx],species[idx],box,sCompY,sCompMw);
        getChemSolve().getCpmixGivenTY(cp[idx],temp[idx],species[idx],box,sCompT,sCompY,sCompCp);
    }
    //
    // Now do pressure relaxation.
    //
    const Real* dx = geom.CellSize();
    FArrayBox p_denom, ugradp, rmix;
    //
    //  Get gamma = c_p/(c_p-R)
    //
    MultiFab gamma(grids,1,nGrow);
    
    for (MFIter mfi(gamma); mfi.isValid(); ++mfi)
    {
        const int idx = mfi.index();
        gamma[idx].copy(cp[idx]);
        rmix.resize(BoxLib::grow(grids[idx],nGrow),1);
        rmix.setVal(rgas);
        rmix.divide(mwmix[idx]);
        gamma[idx].minus(rmix);
        gamma[idx].divide(cp[idx]);
        gamma[idx].invert(1.0);
    }
    
    for (MFIter mfi(dpdt); mfi.isValid(); ++mfi)
    {
        const int i = mfi.index();
        
        dpdt[i].copy(S_old[i],grids[i],pComp,grids[i],0,1);
        dpdt[i].plus(-p_amb);
        dpdt[i].mult(1.0/dt);

        ugradp.resize(grids[i],1);
        
        const int* lo = grids[i].loVect();
        const int* hi = grids[i].hiVect();

        FORT_COMPUTE_UGRADP(rhoRT[i].dataPtr(0), 
                            ARLIM(rhoRT[i].box().loVect()), 
                            ARLIM(rhoRT[i].box().hiVect()),
                            ugradp.dataPtr(), 
                            ARLIM(ugradp.box().loVect()), 
                            ARLIM(ugradp.box().hiVect()),
                            u_mac[0][i].dataPtr(),
                            ARLIM(u_mac[0][i].box().loVect()),
                            ARLIM(u_mac[0][i].box().hiVect()),
                            u_mac[1][i].dataPtr(),
                            ARLIM(u_mac[1][i].box().loVect()),
                            ARLIM(u_mac[1][i].box().hiVect()),
#if (BL_SPACEDIM == 3)
                            u_mac[2][i].dataPtr(),
                            ARLIM(u_mac[2][i].box().loVect()),
                            ARLIM(u_mac[2][i].box().hiVect()),
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
        dpdt[i].minus(ugradp,0,0,1);
        //
        // Make sure to divide by gamma *after* subtracting ugradp.
        //
        dpdt[i].divide(gamma[i]);

        if (dpdt_option == 0)
        {
            dpdt[i].divide(S_old[i],pComp,0,1);
        }
        else if (dpdt_option == 1)
        {
            dpdt[i].divide(p_amb);
        }
        else
        {
            const int ncomp = 1;
            p_denom.resize(grids[i],ncomp);
            p_denom.copy(S_old[i],grids[i],pComp,grids[i],0,ncomp);
            Real num_norm = dpdt[i].norm(0,0,1);
            FabMinMax(p_denom,grids[i],2*num_norm*Real_MIN,p_amb,0,ncomp);
            dpdt[i].divide(p_denom);
        }
        
        dpdt[i].mult(dpdt_factor);
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
        dsdt[k].copy(Divu_new[k],grids[k]);
        dsdt[k].minus(Divu_old[k],grids[k],0,0,1);
        dsdt[k].divide(dt,grids[k],0,1);
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
        S[k].copy(temp[k],box,0,box,Temp,1);
    }
}

//
// Taking an initial guess from S, solve hmix = sum(hl(temp).Yl) for temp.
// hmix and Yl are taken from S (and are assumed to be multiplied by rho
// (i.e. S hold rho.hmix and rho.Yl).
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

void
HeatTransfer::RhoH_to_Temp (MultiFab& S,
                            MultiFab& temp,
                            int       nGrow,
                            int       dominmax)
{
    BL_PROFILE("HeatTransfer::RhoH_to_Temp()");

    BL_ASSERT(S.nGrow() >= nGrow  &&  temp.nGrow() >= nGrow);

    const BoxArray& sgrids = S.boxArray();

    const int sCompT    = 0;
    int       max_iters = 0;

    MultiFab::Copy(temp,S,Temp,sCompT,1,nGrow);

    FArrayBox tmp;

    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        const int k   = mfi.index();
        Box       box = BoxLib::grow(sgrids[k],nGrow);
        //
        // Convert rho*h to h and rho*Y to Y for this operation.
        //
        tmp.resize(box,1);
        tmp.copy(S[k],Density,0,1);
        tmp.invert(1);

        S[k].mult(tmp,0,RhoH,1);

        for (int spec = first_spec; spec <= last_spec; spec++)
            S[k].mult(tmp,0,spec,1);

        int iters = getChemSolve().getTGivenHY(temp[k],S[k],S[k],
                                               box,RhoH,first_spec,sCompT);
        if (iters < 0)
            BoxLib::Error("HeatTransfer::RhoH_to_Temp: error in H->T");

        max_iters = std::max(max_iters,iters);

        if (dominmax)
            FabMinMax(S[k], box, htt_tempmin, htt_tempmax, Temp, 1);
        //
        // Convert back to rho*h and rho*Y
        //
        S[k].mult(S[k],box,Density,RhoH,1);
        for (int spec = first_spec; spec <= last_spec; spec++)
            S[k].mult(S[k],box,Density,spec,1);
    }

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceIntMax(max_iters,IOProc);

        if (verbose && ParallelDescriptor::IOProcessor())
        {
            std::cout << "HeatTransfer::RhoH_to_Temp: max_iters = " << max_iters << '\n';
        }
    }
}

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
		tmp.copy((*temp)[iGrid],0,sCompT,1);
	    }
            else
            {
		tmp.copy(S[iGrid],Temp,sCompT,1);
	    }
	    tmp.copy(S[iGrid],first_spec,sCompY,nCompY);
	}
        else
        {
	    sCompT = Temp;
	    sCompY = first_spec;
	}
	
	const FArrayBox& state = need_tmp_data ? tmp : S[iGrid];

	const int sCompCp = 0;
	getChemSolve().getCpmixGivenTY(cp[iGrid],state,state,box,sCompT,sCompY,sCompCp);

	if (calchmix)
	{
	    const int sCompH = 0;
	    getChemSolve().getHmixGivenTY(hmix[iGrid],state,state,
					  box,sCompT,sCompY,sCompH);
	}
    }
}

void
HeatTransfer::compute_cp (Real      time,
                          MultiFab& cp)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::compute_cp()");

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
	const Box& fabBox   = fpi().box();
	FArrayBox& state    = fpi();

        tmp.resize(validBox,1);
        tmp.copy(state,sCompR,0,1);
        tmp.invert(1);

        for (int k = 0; k < nspecies; k++)
            state.mult(tmp,0,sCompY+k,1);
	
	getChemSolve().getCpmixGivenTY(cp[iGrid],state,state,validBox,
				       sCompT,sCompY,sCompCp);
    }
}

void
HeatTransfer::compute_cp (const FArrayBox& temp, 
                          const FArrayBox& species,
                          FArrayBox&       cp)
{
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
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::compute_rhohmix()");

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
	const Box& fabBox   = fpi().box();
	FArrayBox& state    = fpi();

        tmp.resize(validBox,1);
        tmp.copy(state,sCompR,0,1);
        tmp.invert(1);

        for (int k = 0; k < nspecies; k++)
            state.mult(tmp,0,sCompY+k,1);
	
	getChemSolve().getHmixGivenTY(rhohmix[iGrid],state,state,validBox,
				      sCompT,sCompY,sCompH);
        //
        // Convert hmix to rho*hmix
        //
        rhohmix[iGrid].mult(state,validBox,sCompR,sCompH,1);
    }
}
 
void
HeatTransfer::compute_h (Real      time,
                         MultiFab& h)
{
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
	const Box box = BoxLib::grow(grids[fpi.index()],nGrow);

	BL_ASSERT(box == h[fpi.index()].box());

	getChemSolve().getHGivenT(h[fpi.index()],fpi(),box,sCompT,sCompH);
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

        std::list<std::string>::const_iterator li = parent->statePlotVars().begin();

        for ( ; li != parent->statePlotVars().end(); ++li)
            std::cout << *li << ' ';
        std::cout << '\n';

        std::cout << "\nDerive Plot Vars: ";

        li = parent->derivePlotVars().begin();

        for ( ; li != parent->derivePlotVars().end(); ++li)
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

    for (std::list<DeriveRec>::const_iterator it = dlist.begin();
         it != dlist.end();
         ++it)
    {
        if (parent->isDerivePlotVar(it->name()))
	{
            derive_names.push_back(it->name());
            num_derive += it->numDerive();
	}
    }

    int num_auxDiag = 0;
    if (plot_auxDiags)
        num_auxDiag = auxDiag_names.size();
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

	for (std::list<std::string>::const_iterator it = derive_names.begin();
             it != derive_names.end();
             ++it)
        {
	    const DeriveRec* rec = derive_lst.get(*it);
	    for (i = 0; i < rec->numDerive(); i++)
                os << rec->variableName(i) << '\n';
        }
        //
        // Hack in additional diagnostics.
        //
        if (plot_auxDiags)
            for (i=0; i<auxDiag_names.size(); ++i)
                os << auxDiag_names[i] << '\n';

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
        os << (int) CoordSys::Coord() << '\n';
        os << "0\n"; // Write bndry data.
    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";
    char buf[64];
    sprintf(buf, "Level_%d", level);
    std::string Level = buf;
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
            for (n = 0; n < BL_SPACEDIM; n++)
                os << grid_loc[i].lo(n) << ' ' << grid_loc[i].hi(n) << '\n';
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
	for (std::list<std::string>::const_iterator it = derive_names.begin();
             it != derive_names.end();
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
    if (plot_auxDiags)
        MultiFab::Copy(plotMF,*auxDiag,0,cnt,num_auxDiag,nGrow);
    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(plotMF,TheFullPath,how,true);
}

AuxBoundaryData::AuxBoundaryData ()
    :
    m_ngrow(0),
    m_initialized(false)
{}

AuxBoundaryData::AuxBoundaryData (const BoxArray& ba,
                                  int             n_grow,
                                  int             n_comp,
                                  const Geometry& geom)
    :
    m_ngrow(n_grow),
    m_initialized(false)
{
    initialize(ba,n_grow,n_comp,geom);
}

void
AuxBoundaryData::copy (const AuxBoundaryData& src,
                       int                    src_comp,
                       int                    dst_comp,
                       int                    num_comp)
{
    BL_ASSERT(m_initialized);
    BL_ASSERT(src_comp + num_comp <= src.m_fabs.nComp());
    BL_ASSERT(dst_comp + num_comp <= m_fabs.nComp());

    m_fabs.copy(src.m_fabs,src_comp,dst_comp,num_comp);
}

AuxBoundaryData::AuxBoundaryData (const AuxBoundaryData& rhs)
    :
    m_fabs(rhs.m_fabs.boxArray(),rhs.m_fabs.nComp(),0),
    m_ngrow(rhs.m_ngrow)
{
    m_fabs.copy(rhs.m_fabs,0,0,rhs.m_fabs.nComp());

    m_initialized = true;
}

void
AuxBoundaryData::initialize (const BoxArray& ba,
			     int             n_grow,
			     int             n_comp,
                             const Geometry& geom)
{
    BL_ASSERT(!m_initialized);

    m_ngrow = n_grow;

    const BoxList blgrids = BoxList(ba);

    BoxDomain bd;

    for (int i = 0; i < ba.size(); ++i)
    {
	BoxList gCells = BoxLib::boxDiff(BoxLib::grow(ba[i],n_grow),ba[i]);

	for (BoxList::iterator bli = gCells.begin();
             bli != gCells.end();
             ++bli)
        {
	    bd.add(BoxLib::complementIn(*bli,blgrids));
        }
    }

    if (geom.isAnyPeriodic())
    {
        //
        // Remove any intersections with periodically shifted valid region.
        //
        Box dmn = geom.Domain();

        for (int d = 0; d < BL_SPACEDIM; d++)
            if (!geom.isPeriodic(d)) 
                dmn.grow(d,n_grow);

        BoxList bl;

        for (BoxDomain::const_iterator bdi = bd.begin(); bdi != bd.end(); ++bdi)
        {
            Box isect = *bdi & dmn;

            if (isect.ok())
                bl.push_back(isect);
        }

	bd.clear();

	bd.add(bl);
    }

    bd.simplify();

    BL_ASSERT(bd.ok());

    BoxArray nba(bd.size());

    int cnt = 0;
    for (BoxDomain::const_iterator bdi = bd.begin(); bdi != bd.end(); ++bdi, ++cnt)
    {
        nba.set(cnt,*bdi);
    }

    m_fabs.define(nba,n_comp,0,Fab_allocate);

    m_initialized = true;
}

void
AuxBoundaryData::copyTo (MultiFab& mf,
                         int       src_comp,
                         int       dst_comp,
                         int       num_comp) const
{
    BL_ASSERT(m_initialized);

    mf.copy(m_fabs,src_comp,dst_comp,num_comp);
}

void
AuxBoundaryData::copyFrom (const MultiFab& mf,
                           int       src_comp,
                           int       dst_comp,
                           int       num_comp)
{
    BL_ASSERT(m_initialized);

    m_fabs.copy(mf,src_comp,dst_comp,num_comp);
}

#ifndef NDEBUG
//
// Debugging routines.
//
static int lastOne[2] = {0, 0};
static long Nresets = 0;
int whichBox(const BoxArray& ba, const IntVect& idx, int toggle)
{
    // check first the soln for last time, then do full search
    if (ba[lastOne[toggle]].contains(idx)) return lastOne[toggle];
    for (int i=0; i<ba.size(); ++i)
    {
	if (ba[i].contains(idx))
	{
	    lastOne[toggle] = i;
	    Nresets++;
	    return lastOne[toggle];
	}
    }
    return -1;
}

IntVect myReflect(const IntVect& in, int loc, int dir)
{
  IntVect tmpInt(in);
  tmpInt.setVal(dir,2*loc-in[dir]-1);
  return tmpInt;
}

void
flipFab(FArrayBox* out,
	const FArrayBox* in,
	int dir1,
	int dir2)
{
   const Box& box = in->box();
   IntVect lo(box.smallEnd());
   IntVect hi(box.bigEnd());
   std::swap(lo[dir1],lo[dir2]);
   std::swap(hi[dir1],hi[dir2]);
   Box flipped_box(lo,hi);
   int nComp = in->nComp();
   out->resize(flipped_box,nComp);
   
   for (IntVect fidx = flipped_box.smallEnd(); fidx <= flipped_box.bigEnd();
        box.next(fidx))
   {
     IntVect idx(fidx);
     std::swap(idx[dir1],idx[dir2]);
     for (int k = 0; k < nComp; ++k) 
       (*out)(fidx,k) = (*in)(idx,k); 
   } 
}

void 
flipMF(MultiFab*       out,
       const MultiFab* in,
       int dir1,
       int dir2)
{
  const BoxArray& grids = in->boxArray();
  int len = grids.size();
  BoxArray ba(len);
  for (int i=0; i<len; ++i)
    {
        const Box& gbox = grids[i];
	IntVect lo(gbox.smallEnd());
	IntVect hi(gbox.bigEnd());
	std::swap(lo[dir1],lo[dir2]);
	std::swap(hi[dir1],hi[dir2]);
        Box flipped_box(lo,hi);
        ba.set(i,flipped_box);
    }
    int nGrow = in->nGrow();
    int nComp = in->nComp();
    out->define(ba,nComp,nGrow,Fab_allocate);
    BL_ASSERT(ParallelDescriptor::NProcs()==1);
    for (int i=0; i<out->size(); ++i)
    {
        FArrayBox& fab = (*out)[i];
        const Box& fb = fab.box();
        for (IntVect fidx=fb.smallEnd(); fidx <= fb.bigEnd(); fb.next(fidx))
        {
	    IntVect idx(fidx);
	    std::swap(idx[dir1],idx[dir2]);
            for (int k=0; k<nComp; ++k)
                (*out)[i](fidx,k) = (*in)[i](idx,k);
        }
    }
}

template <class Op>
void
reflectFab(FArrayBox*       out,
	   const FArrayBox* in,
	   int              loc,
	   Op               oper,
	   int              dir)
{
    //
    // Build reflected boxes.
    //
    const Box& box = in->box();
    IntVect lo(box.smallEnd());
    IntVect hi(box.bigEnd());
    lo.setVal(dir,loc-(  box.bigEnd(dir)-loc+1));
    hi.setVal(dir,loc-(box.smallEnd(dir)-loc+1));
    Box reflected_box(lo,hi);
    reflected_box &= box;
    int nComp = in->nComp();
    out->resize(reflected_box,nComp);

    for (IntVect idx=reflected_box.smallEnd(); idx <= reflected_box.bigEnd();
	 reflected_box.next(idx))
    {
	IntVect ridx = myReflect(idx,loc,dir);
	for (int k=0; k<nComp; ++k)
	    (*out)(idx,k) = oper((*in)(ridx,k),(*in)(idx,k));
    }
}

template <class Op>
void
reflectMF (MultiFab*       out,
	   const MultiFab* in,
	   int             loc,
	   Op              oper,
	   int             dir)
{
    // Build reflected boxes
    const BoxArray& grids = in->boxArray();
    int len = grids.size();
    BoxList reflected_bl;
    for (int i=0; i<len; ++i)
    {
	const Box& gbox = grids[i];
	IntVect lo(gbox.smallEnd());
	IntVect hi(gbox.bigEnd());
	lo.setVal(dir,loc-(  gbox.bigEnd(dir)-loc+1));
	hi.setVal(dir,loc-(gbox.smallEnd(dir)-loc+1));
	BoxList ovlps = BoxLib::intersect(BoxList(grids),Box(lo,hi));
	reflected_bl.join(ovlps);
    }
    int nComp = in->nComp();
    out->define(BoxArray(reflected_bl),nComp,0,Fab_allocate);
    BL_ASSERT(ParallelDescriptor::NProcs()==1);
    for (int i=0; i<out->size(); ++i)
    {
	FArrayBox& fab = (*out)[i];
	const Box& fb = fab.box();
	for (IntVect idx=fb.smallEnd(); idx <= fb.bigEnd(); fb.next(idx))
	{
	    IntVect ridx = myReflect(idx,loc,dir);
	    int j = whichBox(grids,idx,0);
	    int rj = whichBox(grids,ridx,1);
	    BL_ASSERT(j>=0 && rj>=0);
	    for (int k=0; k<nComp; ++k)
		fab(idx,k) = oper((*in)[rj](ridx,k),(*in)[j](idx,k));
	}
    }
}

class MFinMemory
{
public:
    MFinMemory() {}
    ~MFinMemory() {}

    bool store(const MultiFab* mf,
	       int             loc)
    {
	if (mf == 0)
	    return false;
	if (mfp.size() <= loc)
	    mfp.resize(loc+1);
	if (mfp.defined(loc))
	    mfp.clear(loc);
	mfp.set(loc, new MultiFab(mf->boxArray(),mf->nComp(),mf->nGrow()));
	MultiFab::Copy(mfp[loc],*mf,0,0,mf->nComp(),mf->nGrow());
	return true;
    }

    const MultiFab* recall(int loc)
    {
	return (mfp.size() <= loc  ?  0  : &mfp[loc]);
    }

private:
    PArray<MultiFab> mfp;
};

static MFinMemory mfInMemory;

#if BL_SPACEDIM==2
static std::string AmrVis = "amrvis2d";
#else
static std::string AmrVis = "amrvis3d";
#endif

extern "C"
{
    bool dumpFab(const FArrayBox* fab,
		 const char* file)
      {
	std::ofstream os(file);
	os.precision(20);
	os << *fab << std::endl;
	return true;
      }

    bool dumpMF(const MultiFab* mf,
		 const char* file)
      {
	std::ofstream os(file);
	os.precision(20);
	for (int i=0; i < mf->size() ; i++)
	  os << (*mf)[i] << std::endl;
	return true;
      }

    bool dumpFabcomp(const FArrayBox* fab,
		     const char* file,
		     int sComp,
		     int nComp)
      {
	BL_ASSERT (nComp>0);
	BL_ASSERT(fab->nComp()>=sComp+nComp);
	FArrayBox junk(fab->box(),nComp);
	junk.copy(*fab,sComp,0,nComp);
	return dumpFab(&junk,file);
      }

    bool dumpMFcomp(const MultiFab* mf,
		    const char* file,
		    int sComp,
		    int nComp)
      {
	BL_ASSERT (nComp>0);
	BL_ASSERT(mf->nComp()>=sComp+nComp);
	MultiFab junk(mf->boxArray(),nComp,mf->nGrow());
	MultiFab::Copy(junk,*mf,sComp,0,nComp,mf->nGrow());
	return dumpMF(&junk,file);
      }

    bool writeFab (const FArrayBox* fab,
		   const char*      file)
	{
	    std::ofstream os(file);
	    fab->writeOn(os);
	    return true;
	}
    
    bool writeMF (const MultiFab* mf,
		  const char*     file)
	{
	    std::string FullPath = file;
	    if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
		FullPath += '/';
	    
	    if (ParallelDescriptor::IOProcessor())
		if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
		    BoxLib::CreateDirectoryFailed(FullPath);

            ParallelDescriptor::Barrier();
	    
	    static const std::string MultiFabBaseName("MultiFab");
	    FullPath += MultiFabBaseName;
	    
	    VisMF::Write(*mf,FullPath,VisMF::OneFilePerCPU);
	    
	    return true;
	}

    bool writeMFcomp(const MultiFab* mf,
		     const char*     file,
		     int             sComp,
		     int             nComp)
	{
	    BL_ASSERT(nComp>0);
	    BL_ASSERT(mf->nComp()>=sComp+nComp);
	    const int nGrow = 0;
	    MultiFab junk(mf->boxArray(),nComp,nGrow);
	    MultiFab::Copy(junk,*mf,sComp,0,nComp,nGrow);
	    return writeMF(&junk,file);
	}

    bool amrvisFab(const FArrayBox* fab)
	{
	    std::string buf("amrvis.fab");
	    writeFab(fab,buf.c_str());
	    buf = AmrVis + std::string(" -fab ") + buf + std::string(" &");
	    system(buf.c_str());
	    return true;
	}

    bool amrvisFabBox(const FArrayBox* fab,
                      const Box&       box)
	{
            FArrayBox tmp(box,fab->nComp());
            tmp.copy(*fab);
	    std::string buf("amrvis.fab");
	    writeFab(&tmp,buf.c_str());
	    buf = AmrVis + std::string(" -fab ") + buf + std::string(" &");
	    system(buf.c_str());
	    return true;
	}
    
    bool amrvisFabDiff(const FArrayBox* fab1, const FArrayBox* fab2)
	{
	    const Box box = fab1->box() & fab2->box();
	    const int nComp = std::min(fab1->nComp(),fab2->nComp());
	    if (box.ok() && nComp > 0)
	    {
		FArrayBox fab(box,nComp);
		fab.copy(*fab1,box,0,box,0,nComp);
		fab.minus(*fab2,box,0,0,nComp);
		return amrvisFab(&fab);
	    }
	    else
	    {
		BoxLib::Warning("..........fabdiff: Bad intersection");
	    }
	    return true;
	}
    
    bool amrvisMF(const MultiFab* mf)
	{
	    std::string buf("amrvis.mfab");
	    writeMF(mf,buf.c_str());
	    buf = AmrVis + std::string(" -multifab ") + buf + std::string("/MultiFab &");
	    system(buf.c_str());
	    return true;
	}

    bool amrvisMFcomp(const MultiFab* mf,
		      int             sComp,
		      int             nComp = 1)
	{
	    const int nGrow = 0;
	    MultiFab junk(mf->boxArray(),nComp,nGrow);
	    for (MFIter mfi(junk); mfi.isValid(); ++mfi)
	    {
		const Box box = BoxLib::grow(mfi.validbox(),nGrow);
		junk[mfi].copy((*mf)[mfi],box,sComp,box,0,nComp);
	    }
	    return amrvisMF(&junk);
	}

    bool amrvisFabComp(const FArrayBox* fab,
		       int              sComp,
		       int              nComp = 1,
		       int              nGrow = 0)
	{
	    const Box box = BoxLib::grow(fab->box(),-nGrow);
	    FArrayBox junk(box,nComp);
	    junk.copy(*fab,box,sComp,box,0,nComp);
	    return amrvisFab(&junk);
	}

    bool amrvisMFdiff(const MultiFab* mf1, const MultiFab* mf2, int sComp, int nComp)
	{
	    const int nGrow = 0;
	    BL_ASSERT(nComp>0);
	    BL_ASSERT(mf1->boxArray() == mf2->boxArray());
	    BL_ASSERT (mf1->nComp() >= sComp+nComp);
	    BL_ASSERT (mf2->nComp() >= sComp+nComp);
	    MultiFab junk(mf1->boxArray(),nComp,nGrow);
	    for (MFIter mfi(junk); mfi.isValid(); ++mfi)
	    {
		const Box box = BoxLib::grow(mf1->boxArray()[mfi.index()],nGrow);
		junk[mfi].copy((*mf1)[mfi],box,0,box,0,nComp);
		junk[mfi].minus((*mf2)[mfi],box,0,0,nComp);
	    }
	    return amrvisMF(&junk);
	}
    
    bool amrvisXsymFab(const FArrayBox* fab,
                       int              loc)
	{
	    int dir = 0;
	    FArrayBox junk;
	    reflectFab(&junk,fab,loc,std::minus<Real>(),dir);
	    return amrvisFab(&junk);
	}
    
    bool amrvisXantiSymFab(const FArrayBox* fab,
                           int              loc)
	{
	    int dir = 0;
	    FArrayBox junk;
	    reflectFab(&junk,fab,loc,std::plus<Real>(),dir);
	    return amrvisFab(&junk);
	}
    
    bool amrvisXsymMF(const MultiFab* mf,
		      int             loc)
	{
	    int dir = 0;
	    MultiFab junk;
	    reflectMF(&junk,mf,loc,std::minus<Real>(),dir);
	    return amrvisMF(&junk);
	}

    bool amrvisXantiSymMF(const MultiFab* mf,
                          int             loc)
	{
	    int dir = 0;
	    MultiFab junk;
	    reflectMF(&junk,mf,loc,std::plus<Real>(),dir);
	    return amrvisMF(&junk);
	}


    bool amrvisYsymFab(const FArrayBox* fab,
                       int              loc)
	{
	    int dir = 1;
	    FArrayBox junk;
	    reflectFab(&junk,fab,loc,std::minus<Real>(),dir);
	    return amrvisFab(&junk);
	}
    
    bool amrvisYantiSymFab(const FArrayBox* fab,
                           int              loc)
	{
	    int dir = 1;
	    FArrayBox junk;
	    reflectFab(&junk,fab,loc,std::plus<Real>(),dir);
	    return amrvisFab(&junk);
	}
    
    bool amrvisYsymMF(const MultiFab* mf,
		      int             loc)
	{
	    int dir = 1;
	    MultiFab junk;
	    reflectMF(&junk,mf,loc,std::minus<Real>(),dir);
	    return amrvisMF(&junk);
	}

    bool amrvisYantiSymMF(const MultiFab* mf,
                          int             loc)
	{
	    int dir = 1;
	    MultiFab junk;
	    reflectMF(&junk,mf,loc,std::plus<Real>(),dir);
	    return amrvisMF(&junk);
	}


    bool amrvisZsymFab(const FArrayBox* fab,
                       int              loc)
	{
	    int dir = 2;
	    FArrayBox junk;
	    reflectFab(&junk,fab,loc,std::minus<Real>(),dir);
	    return amrvisFab(&junk);
	}
    
    bool amrvisZantiSymFab(const FArrayBox* fab,
                           int              loc)
	{
	    int dir = 2;
	    FArrayBox junk;
	    reflectFab(&junk,fab,loc,std::plus<Real>(),dir);
	    return amrvisFab(&junk);
	}
    
    bool amrvisZsymMF(const MultiFab* mf,
		      int             loc)
	{
	    int dir = 2;
	    MultiFab junk;
	    reflectMF(&junk,mf,loc,std::minus<Real>(),dir);
	    return amrvisMF(&junk);
	}

    bool amrvisZantiSymMF(const MultiFab* mf,
                          int             loc)
	{
	    int dir = 2;
	    MultiFab junk;
	    reflectMF(&junk,mf,loc,std::plus<Real>(),dir);
	    return amrvisMF(&junk);
	}

    bool amrvisSumMF(const MultiFab* mf,
		     int             sComp,
		     int             nComp)
	{
	    const int nGrow = 0;
	    MultiFab junk(mf->boxArray(),1,nGrow);
	    junk.setVal(0.0);
	    for (MFIter mfi(junk); mfi.isValid(); ++mfi)
	    {
		const Box box = BoxLib::grow(mfi.validbox(),nGrow);
		for (int n=0; n<nComp; n++)
		    junk[mfi].plus((*mf)[mfi],box,sComp+n,0,1);
	    }
	    return amrvisMF(&junk);
	}

    bool amrvisMFsumY(const MultiFab* mf)
	{
	    // Specific to my hack of HT...
	    const int nComp = HeatTransfer::getChemSolve().numSpecies();
	    const int nGrow = 0;
	    const int sComp = BL_SPACEDIM+2;
	    const int rhoComp = BL_SPACEDIM;
	    MultiFab junk(mf->boxArray(),1,nGrow);
	    junk.setVal(0.0);
	    for (MFIter mfi(junk); mfi.isValid(); ++mfi)
	    {
		const Box box = BoxLib::grow(mfi.validbox(),nGrow);
		for (int n=0; n<nComp; n++)
		    junk[mfi].plus((*mf)[mfi],box,sComp+n,0,1);
		junk[mfi].divide((*mf)[mfi],box,rhoComp,0,1);
	    }
	    return amrvisMF(&junk);
	}

    bool amrvisMFY(const MultiFab* mf)
	{
	    // Specific to my hack of HT...
	    const int nComp = HeatTransfer::getChemSolve().numSpecies();
	    const int nGrow = 0;
	    const int sComp = BL_SPACEDIM+2;
	    const int rhoComp = BL_SPACEDIM;
	    MultiFab junk(mf->boxArray(),nComp,nGrow);
	    MultiFab::Copy(junk,*mf,sComp,0,nComp,nGrow);
	    for (MFIter mfi(junk); mfi.isValid(); ++mfi)
	    {
		const Box box = BoxLib::grow(mfi.validbox(),nGrow);
		for (int n=0; n<nComp; n++)
		    junk[mfi].divide((*mf)[mfi],box,rhoComp,n,1);
	    }
	    return amrvisMF(&junk);
	}

    bool amrvisFabSumY(const FArrayBox* fab,
		       const int        nGrow)
	{
	    // Specific to my hack of HT...
	    const int nComp = HeatTransfer::getChemSolve().numSpecies();
	    const int sComp = BL_SPACEDIM+2;
	    const int rhoComp = BL_SPACEDIM;
	    const Box box = BoxLib::grow(fab->box(),-nGrow);
	    FArrayBox junk(box,1);
	    junk.setVal(0.0);
	    for (int n=0; n<nComp; ++n)
		junk.plus((*fab),box,sComp+n,0,1);
	    junk.divide((*fab),box,rhoComp,0,1);
	    return amrvisFab(&junk);
	}
	    

    bool amrvisSumFab(const FArrayBox* fab,
		      int              sComp,
		      int              nComp)
	{
	    FArrayBox junk(fab->box(),1);
	    junk.setVal(0.0);
	    for (int n=0; n<nComp; n++)
		junk.plus(*fab,sComp+n,0,1);
	    return amrvisFab(&junk);
	}

    bool amrvisFlipFabinXY(const FArrayBox* fab)
      {
	FArrayBox junk;
	int dir1 = 0;
	int dir2 = 1;
	flipFab(&junk,fab,dir1,dir2);
	return amrvisFab(&junk);
      }

    bool amrvisFlipMFinXY(const MultiFab* mf)
      {
	MultiFab junk;
	int dir1 = 0;
	int dir2 = 1;
	flipMF(&junk, mf,dir1,dir2);
	return amrvisMF(&junk);
      }

    bool writeFlipFabinXY(const FArrayBox* fab,
			  const char* file)
      {
	FArrayBox junk;
	int dir1 = 0;
	int dir2 = 1;
	flipFab(&junk,fab,dir1,dir2);
	return writeFab(&junk,file);
      }

    bool writeFlipMFinXY(const MultiFab* mf,
			 const char* file)
      {
	MultiFab junk;
	int dir1 = 0;
	int dir2 = 1;
	flipMF(&junk, mf,dir1,dir2);
	return writeMF(&junk,file);
      }

    bool dumpFlipFabinXY(const FArrayBox* fab,
			const char* file)
      {
	FArrayBox junk;
	int dir1 = 0;
	int dir2 = 1;
	flipFab(&junk, fab,dir1,dir2);
	return dumpFab(&junk,file);
      }

    bool dumpFlipMFinXY(const MultiFab* mf,
			const char* file)
      {
	MultiFab junk;
	int dir1 = 0;
	int dir2 = 1;
	flipMF(&junk, mf,dir1,dir2);
	return dumpMF(&junk,file);
      }


    bool amrvisFlipFabinXZ(const FArrayBox* fab)
      {
	FArrayBox junk;
	int dir1 = 0;
	int dir2 = 2;
	flipFab(&junk,fab,dir1,dir2);
	return amrvisFab(&junk);
      }

    bool amrvisFlipMFinXZ(const MultiFab* mf)
      {
	MultiFab junk;
	int dir1 = 0;
	int dir2 = 2;
	flipMF(&junk, mf,dir1,dir2);
	return amrvisMF(&junk);
      }

    bool writeFlipFabinXZ(const FArrayBox* fab,
			  const char* file)
      {
	FArrayBox junk;
	int dir1 = 0;
	int dir2 = 2;
	flipFab(&junk,fab,dir1,dir2);
	return writeFab(&junk,file);
      }

    bool writeFlipMFinXZ(const MultiFab* mf,
			 const char* file)
      {
	MultiFab junk;
	int dir1 = 0;
	int dir2 = 2;
	flipMF(&junk, mf,dir1,dir2);
	return writeMF(&junk,file);
      }

    bool dumpFlipFabinXZ(const FArrayBox* fab,
			const char* file)
      {
	FArrayBox junk;
	int dir1 = 0;
	int dir2 = 2;
	flipFab(&junk, fab,dir1,dir2);
	return dumpFab(&junk,file);
      }

    bool dumpFlipMFinXZ(const MultiFab* mf,
			const char* file)
      {
	MultiFab junk;
	int dir1 = 0;
	int dir2 = 2;
	flipMF(&junk, mf,dir1,dir2);
	return dumpMF(&junk,file);
      }


    bool amrvisFlipFabinYZ(const FArrayBox* fab)
      {
	FArrayBox junk;
	int dir1 = 1;
	int dir2 = 2;
	flipFab(&junk,fab,dir1,dir2);
	return amrvisFab(&junk);
      }

    bool amrvisFlipMFinYZ(const MultiFab* mf)
      {
	MultiFab junk;
	int dir1 = 1;
	int dir2 = 2;
	flipMF(&junk, mf,dir1,dir2);
	return amrvisMF(&junk);
      }

    bool writeFlipFabinYZ(const FArrayBox* fab,
			  const char* file)
      {
	FArrayBox junk;
	int dir1 = 1;
	int dir2 = 2;
	flipFab(&junk,fab,dir1,dir2);
	return writeFab(&junk,file);
      }

    bool writeFlipMFinYZ(const MultiFab* mf,
			 const char* file)
      {
	MultiFab junk;
	int dir1 = 1;
	int dir2 = 1;
	flipMF(&junk, mf,dir1,dir2);
	return writeMF(&junk,file);
      }

    bool dumpFlipFabinYZ(const FArrayBox* fab,
			const char* file)
      {
	FArrayBox junk;
	int dir1 = 1;
	int dir2 = 2;
	flipFab(&junk, fab,dir1,dir2);
	return dumpFab(&junk,file);
      }

    bool dumpFlipMFinYZ(const MultiFab* mf,
			const char* file)
      {
	MultiFab junk;
	int dir1 = 1;
	int dir2 = 2;
	flipMF(&junk, mf,dir1,dir2);
	return dumpMF(&junk,file);
      }


    bool storeMF(const MultiFab* mf,
		 int             loc)
	{
	    mfInMemory.store(mf,loc);
	    return true;
	}

    const MultiFab* recallMF(int loc)
	{
	    return mfInMemory.recall(loc);
	}
}

//
// Need this so linker wont toss out the uncalled function.
//
void NeverCalledHT()
{
    Box garbage_box;
    dumpFab(NULL,"DORK");
    dumpMF(NULL,"DORK");
    dumpFabcomp(NULL,"DORK",0,0);
    dumpMFcomp(NULL,"DORK",0,0);
    writeMFcomp(NULL,"DORK",0,0);
    writeFab( NULL, "DORK" );
    writeMF( NULL, "DORK" );
    amrvisFab( NULL );
    amrvisFabBox( NULL, garbage_box);
    amrvisFabDiff( NULL, NULL );
    amrvisMF( NULL );
    amrvisMFdiff(NULL, NULL, 0, 0 );
    amrvisMFcomp( NULL, 0, 0 );
    amrvisXsymFab( NULL, 0 );
    amrvisXsymMF( NULL, 0 );
    amrvisYsymFab( NULL, 0 );
    amrvisYsymMF( NULL, 0 );
    amrvisZsymFab( NULL, 0 );
    amrvisZsymMF( NULL, 0 );
    amrvisSumMF( NULL, 0, 0);
    amrvisMFsumY( NULL );
    amrvisMFY( NULL );
    amrvisFabSumY( NULL, 0 );
    amrvisSumFab( NULL, 0, 0);
    amrvisFlipFabinXY(NULL);
    amrvisFlipMFinXY(NULL);
    writeFlipFabinXY(NULL,"DORK");
    writeFlipMFinXY(NULL,"DORK");
    dumpFlipFabinXY(NULL,"DORK");
    dumpFlipMFinXY(NULL,"DORK");
    amrvisFlipFabinXZ(NULL);
    amrvisFlipMFinXZ(NULL);
    writeFlipFabinXZ(NULL,"DORK");
    writeFlipMFinXZ(NULL,"DORK");
    dumpFlipFabinXZ(NULL,"DORK");
    dumpFlipMFinXZ(NULL,"DORK");
    amrvisFlipFabinYZ(NULL);
    amrvisFlipMFinYZ(NULL);
    writeFlipFabinYZ(NULL,"DORK");
    writeFlipMFinYZ(NULL,"DORK");
    dumpFlipFabinYZ(NULL,"DORK");
    dumpFlipMFinYZ(NULL,"DORK");
    storeMF( NULL, 0 );
    recallMF(0);
}

#endif /*NDEBUG*/
