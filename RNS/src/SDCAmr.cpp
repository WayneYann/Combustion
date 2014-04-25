/*
 * Multilevel SDC + AMR controller.
 *
 * When RNS is compiled with USE_SDCLIB=TRUE, the time-stepping is
 * done using the multi-level SDC (MLSDC) algorithm, with IMEX
 * sweepers on each level (hydrodynamics are explicit, chemistry is
 * implicit).  The MLSDC algorithm is implemented in C in SDCLib (in
 * SDCLib multi-level SDC is called multi-grid SDC).
 *
 * The interface between SDCLib and RNS is (mostly) contained in
 * SDCAmr (derived from Amr) and SDCAmrEncap.
 *
 * Note that in MLSDC, there is no concept of "sub-cycling" fine
 * levels.  As such, the "timeStep" method defined in SDCAmr is only
 * called on the coarsest level, and it essentially passes control to
 * SDCLib to advance the solution.
 *
 * SDCLib handles interpolation and restriction in time.  SDCLib also
 * takes care of allocating multifabs at each SDC node.  As such the
 * traditional Amr concepts of "old data", "new data", and defining an
 * "advance" routine no longer apply.  In order to reuse as much of
 * the existing Amr code as possible, SDCAmr puts the solution
 * obtained from SDCLib into the "new data" state (this means all the
 * logic to write plotfiles etc still works), but never defines "old
 * data" or any "mid data".
 *
 * Some notes:
 *
 * 1. Since we never use 'old data', any calls to fill patch iterator
 *    will always use 'new data', regardless of any 'time' information
 *    set on the state data.
 *
 * Known issues:
 *
 * 1. The SDC hierarchy is currently not depth limited.
 *
 * 2. We're using Gauss-Lobatto nodes, but Gauss-Radau would probably
 *    be better for chemistry.
 *
 * 3. mlsdc_amr_interpolate won't work for correcton at wall boundary.
 */

#include <SDCAmr.H>
#include <MultiFab.H>
#include <ParmParse.H>
#include <StateDescriptor.H>
#include <AmrLevel.H>
#include <Interpolater.H>
#include <FabArray.H>
#include <iomanip>

#include "RNS.H"
#include "RNS_F.H"

#ifdef BL_USE_ARRAYVIEW
#include <ArrayView.H>
#endif

#ifdef USE_COLOROUTPUT
// only works on some systems
#define RESETCOLOR       "\033[0m"
#define BOLDFONT         "\033[1m"
#define REDCOLOR         "\033[31m"      /* Red */
#define GREENCOLOR       "\033[32m"      /* Green */
#define YELLOWCOLOR      "\033[33m"      /* Yellow */
#define BLUECOLOR        "\033[34m"      /* Blue */
#define MAGENTACOLOR     "\033[35m"      /* Magenta */
#define CYANCOLOR        "\033[36m"      /* Cyan */
#else
#define RESETCOLOR       ""
#define BOLDFONT         ""
#define REDCOLOR         ""
#define GREENCOLOR       ""
#define YELLOWCOLOR      ""
#define BLUECOLOR        ""
#define MAGENTACOLOR     ""
#define CYANCOLOR        ""
#endif

using namespace std;

BEGIN_EXTERN_C

/*
 * Spatial interpolation between MultiFabs.  Called by SDCLib.
 */
void mlsdc_amr_interpolate(void *Fp, void *Gp, sdc_state *state, void *ctxF, void *ctxG)
{
  BL_PROFILE("MLSDC_AMR_INTERPOLATE()");

  bool isCorrection = state->kind == SDC_CORRECTION;

  RNSEncap& F      = *((RNSEncap*) Fp);
  RNSEncap& G      = *((RNSEncap*) Gp);
  MultiFab& UF     = *F.U;
  MultiFab& UG     = *G.U;
  RNS&      levelF = *((RNS*) ctxF);
  RNS&      levelG = *((RNS*) ctxG);

  const IntVect         ratio = levelG.fineRatio();
  const DescriptorList& dl    = levelF.get_desc_lst();
  const Array<BCRec>&   bcs   = dl[0].getBCs();
  const int             ncomp = dl[0].nComp();
  Interpolater&         map   = *dl[0].interp();
  Array<BCRec>          bcr(ncomp);

  const Geometry& geomG = levelG.Geom();

  RNS_SETNAN(UF);

  // make a coarse version (UC) of the fine multifab (UF)
  BoxArray ba_C(UF.size());
  for (int i=0; i<ba_C.size(); i++)
    ba_C.set(i, map.CoarseBox(UF.fabbox(i), ratio));

  MultiFab UC(ba_C, ncomp, 0);
  RNS_SETNAN(UC);

  bool touch = false;
  bool touch_periodic = false;
  const Box& crse_domain_box = levelG.Domain();
  if (geomG.isAnyPeriodic()) {
      for (int i=0; i < ba_C.size(); i++) {
	  if (! crse_domain_box.contains(ba_C[i])) {
	      touch = true;
	      for (int idim=0; idim<BL_SPACEDIM; idim++) {
		  if (geomG.isPeriodic(i)   ||
		      ba_C[i].bigEnd(idim) > crse_domain_box.bigEnd(idim) ||
		      ba_C[i].smallEnd(idim) < crse_domain_box.smallEnd(idim) )
		  {
		      touch_periodic = true;
		      break;
		  }
	      }
	  }
	  if (touch_periodic) break;
      }
  }
  else {
      for (int i=0; i < ba_C.size(); i++) {
	  if (! crse_domain_box.contains(ba_C[i])) {
	      touch = true;
	      break;
	  }
      }
  }

  if (!touch) {
    // If level F does not touch physical boundaries, then AMR levels are
    // properly nested so that the valid rgions of UC are contained inside
    // the valid regions of UG.  So FabArray::copy is all we need.
    UC.copy(UG);
  }
  else if (touch_periodic) {

      // This case is more complicated because level F might touch only one
      // of the periodic boundaries.
      int ng_C = 0;
      {
	  Box box_C = ba_C.minimalBox();
	  for (int idim=0; idim < BL_SPACEDIM; idim++) {
	      int gap_hi = box_C.bigEnd(idim) - crse_domain_box.bigEnd(idim);
	      int gap_lo = crse_domain_box.smallEnd(idim) - box_C.smallEnd(idim);
	      ng_C = std::max(ng_C, gap_hi);
	      ng_C = std::max(ng_C, gap_lo);
	  }
      }
      int ng_G = UG.nGrow();
      const BoxArray& ba_G = UG.boxArray();

      MultiFab* UG_safe;
      MultiFab UGG;

      if (ng_C > ng_G) {
	  UGG.define(ba_G, ncomp, ng_C, Fab_allocate);
	  RNS_SETNAN(UGG);
	  MultiFab::Copy(UGG, UG, 0, 0, ncomp, 0);
	  UG_safe = &UGG;
      }
      else {
	  UG_safe = &UG;
      }

      if (isCorrection) UG_safe->setBndry(0.0);

      levelG.fill_boundary(*UG_safe, state->t, RNS::use_FillBoundary, isCorrection);

      // We cannot do FabArray::copy() directly on UG because it copies only form
      // valid regions.  So we need to make the ghost cells of UG valid.
      BoxArray ba_G2(UG.size());
      for (int i=0; i<ba_G2.size(); i++) {
	  ba_G2.set(i, BoxLib::grow(ba_G[i],ng_C));
      }

      MultiFab UG2(ba_G2, ncomp, 0);
      for (MFIter mfi(UG2); mfi.isValid(); ++mfi)
      {
	  int i = mfi.index();
	  UG2[i].copy((*UG_safe)[i]);  // Fab to Fab copy
      }

      UC.copy(UG2);

  }
  else {
      if (isCorrection) {
#ifdef NDEBUG
	  UC.setVal(0.0);
	  UC.copy(UG);
#else
	  UC.copy(UG);
	  const Box& crse_domain_box = levelG.Domain();
	  for (MFIter mfi(UC); mfi.isValid(); ++mfi) {
	      UC[mfi].setComplement(0.0, crse_domain_box, 0, ncomp);
	  }
#endif
      }
      else {
	  UC.copy(UG);
	  levelG.fill_boundary(UC, state->t, RNS::set_PhysBoundary, isCorrection);
      }
  }

  RNS_ASSERTNONAN(UC);

  // now that UF is completely contained within UC, cycle through each
  // FAB in UF and interpolate from the corresponding FAB in UC
// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
  for (MFIter mfi(UF); mfi.isValid(); ++mfi) {
    BoxLib::setBC(UF[mfi].box(), levelF.Domain(), 0, 0, ncomp, bcs, bcr);
    Geometry fine_geom(UF[mfi].box());
    Geometry crse_geom(UC[mfi].box());

    map.interp(UC[mfi], 0, UF[mfi], 0, ncomp, UF[mfi].box(), ratio,
               crse_geom, fine_geom, bcr, 0, 0);
  }

  if (isCorrection) UF.setBndry(0.0);
  levelF.fill_boundary(UF, state->t, RNS::set_PhysBoundary, isCorrection);
  RNS_ASSERTNONAN(UF);
}

/*
 * Spatial restriction between solutions and integrals of function
 * evals.  Called by SDCLib.
 */
void mlsdc_amr_restrict(void *Fp, void *Gp, sdc_state *state, void *ctxF, void *ctxG)
{
  BL_PROFILE("MLSDC_AMR_RESTRICT()");

  RNSEncap& F      = *((RNSEncap*) Fp);
  RNSEncap& G      = *((RNSEncap*) Gp);
  MultiFab& UF     = *F.U;
  MultiFab& UG     = *G.U;
  RNS&      levelF = *((RNS*) ctxF);
  RNS&      levelG = *((RNS*) ctxG);

  BL_ASSERT(G.type==SDC_SOLUTION || G.type==SDC_TAU);

  if (G.type == SDC_TAU) {
    SDCAmr&       amr        = *levelF.getSDCAmr();
    IntVect       crse_ratio = amr.refRatio(levelG.Level());
    FluxRegister& flxF       = *F.fine_flux;
    FluxRegister& flxG       = *G.crse_flux;

    FluxRegister flx(UF.boxArray(), crse_ratio, levelF.Level(), UF.nComp());
    for (OrientationIter face; face; ++face)
      for (FabSetIter bfsi(flxF[face()]); bfsi.isValid(); ++bfsi) {
        flx[face()][bfsi].copy(flxF[face()][bfsi]);
        flx[face()][bfsi].saxpy(1.0, flxG[face()][bfsi]);
      }

    flx.Reflux(UG, levelG.Volume(), 1.0, 0, 0, UG.nComp(), levelG.Geom());
  }

  levelG.avgDown(UG, UF);
  if (G.type == SDC_SOLUTION)
    levelG.fill_boundary(UG, state->t, RNS::use_FillBoundary);

#ifndef NDEBUG
  BL_ASSERT(UG.contains_nan() == false);
#endif
}

END_EXTERN_C


void SDCAmr::final_integrate(double t, double dt, int niter)
{
  int nlevs = mg.nlevels;
  for (int l=0; l<nlevs; l++) {
    sdc_sweeper* swp    = mg.sweepers[l];
    swp->sweep(swp, t, dt, niter, SDC_SWEEP_NOEVAL);
  }
}

/*
 * Take one SDC+AMR time step.
 *
 * This is only called on the coarsest SDC+AMR level.
 */
void SDCAmr::timeStep(int level, Real time,
		      int /* iteration */, int /* niter */,
		      Real stop_time)
{
  BL_PROFILE("SDCAmr::timeStep()");

  BL_ASSERT(level == 0);

  // regrid
  int lev_top = min(finest_level, max_level-1);
  for (int i=level; i<=lev_top; i++) {
    const int post_regrid_flag = 1;
    const int old_finest = finest_level;

    // grids_changed = grids_changed || regrid(i,time);
    regrid(i,time);
    amr_level[0].computeNewDt(finest_level, sub_cycle, n_cycle, ref_ratio,
  			      dt_min, dt_level, stop_time, post_regrid_flag);

    for (int k=i; k<=finest_level; k++)
      level_count[k] = 0;

    if (old_finest > finest_level)
      lev_top = min(finest_level, max_level-1);
  }

  // if (mg.sweepers[0] == NULL) grids_changed = true;
  // if (grids_changed)
  rebuild_mlsdc();

  double dt = dt_level[0];
  for (int lev=1; lev<first_refinement_level; lev++)
    dt /= trat;

  if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
    cout << "MLSDC advancing with dt: " << dt << " (" << dt_level[0] << ")" << endl;
  }

  // set intial conditions
  for (int lev=0; lev<=finest_level; lev++) {
    MultiFab& Unew = getLevel(lev).get_new_data(0);
    RNSEncap& Q0   = *((RNSEncap*) mg.sweepers[lev]->nset->Q[0]);
    MultiFab& U0   = *Q0.U;
    MultiFab::Copy(U0, Unew, 0, 0, U0.nComp(), U0.nGrow());
    getLevel(lev).get_state_data(0).setTimeLevel(time+dt, dt, dt);
  }

  // fill fine boundaries using coarse data
  for (int lev=0; lev<=finest_level; lev++) {
    RNS&      rns  = *dynamic_cast<RNS*>(&getLevel(lev));
    RNSEncap& Q0   = *((RNSEncap*) mg.sweepers[lev]->nset->Q[0]);
    MultiFab& U0   = *Q0.U;
    rns.fill_boundary(U0, time, RNS::use_FillCoarsePatch);
  }

  BL_PROFILE_VAR("SDCAmr::timeStep-iters", sdc_iters);

  Array< Array<Real> > r0p(finest_level+1, Array<Real>(3));
  Array< Array<Real> > r2p(finest_level+1, Array<Real>(3));

  sdc_mg_spread_per_level(&mg, time, dt);

  for (int k=0; k<max_iters; k++) {
    int flags = SDC_MG_MIXEDINTERP | SDC_SWEEP_MONITOR;
    if (k==max_iters-1) flags |= SDC_MG_HALFSWEEP;
    if (k==0)           flags |= SDC_SWEEP_FIRST;

    sdc_mg_sweep(&mg, time, dt, k, flags);

    if (verbose > 0) {

      int len=string(" iter: 0, level: 0,   rho").size();

      int ncomps = 1;
      if (RNS::fuelID >= 0) ncomps++;
      if (RNS::flameTracID >= 0) ncomps++;
      Array<int> comps(ncomps);
      Array<string> names(ncomps);
      Array<string> colors(ncomps);

      int ic=0;
      comps[ic] = 0;
      names[ic] = "rho";
      colors[ic] = REDCOLOR;
      if (RNS::fuelID >= 0) {
	ic++;
        comps[ic] = RNS::FirstSpec + RNS::fuelID;
	string vname = "rho*Y(" + RNS::fuelName + ")";
	string space(len-vname.size(), ' ');
	names[ic] = space + vname;
	colors[ic] = GREENCOLOR;
      }
      if (RNS::flameTracID >= 0) {
	ic++;
        comps[ic] = RNS::FirstSpec + RNS::flameTracID;
	string vname = "rho*Y(" + RNS::flameTracName + ")";
	string space(len-vname.size(), ' ');
	names[ic] = space + vname;
	colors[ic] = BLUECOLOR;
      }

      for (int lev=0; lev<=finest_level; lev++) {
	RNSEncap* Rencap = (RNSEncap*) encaps[lev]->create(SDC_INTEGRAL, encaps[lev]->ctx);
        MultiFab& R      = *Rencap->U;

	sdc_sweeper_residual(mg.sweepers[lev], dt, Rencap);

	Array<Real> r0(R.norm0(comps));
	Array<Real> r2(R.norm2(comps));

	encaps[lev]->destroy(Rencap);

	if (ParallelDescriptor::IOProcessor()) {
	  std::ios_base::fmtflags ff = cout.flags();
	  int oldprec = cout.precision(2);
	  cout << scientific;
	  if (level%2 == 1) cout << BOLDFONT;
	  if (k == 0) {
	      cout << " iter: " << k << ", level: " << lev << ",   ";
	      for (int ic=0; ic<ncomps; ic++) {
		  cout << colors[ic] << names[ic]
		       << " res norm0: " << r0[ic] << "         "
		       <<    "  norm2: " << r2[ic] << endl;
	      }
	  }
	  else {
	      cout << " iter: " << k << ", level: " << lev << ",   ";
	      for (int ic=0; ic<ncomps; ic++) {
		  cout << colors[ic] << names[ic]
		       << " res norm0: " << r0[ic] << " " << r0p[lev][ic]/(r0[ic]+1.e-80)
		       <<    ", norm2: " << r2[ic] << " " << r2p[lev][ic]/(r2[ic]+1.e-80) << endl;
	      }
	  }
	  cout << RESETCOLOR;
	  cout.precision(oldprec);
	  cout.flags(ff);
	}

	for (int ic=0; ic<ncomps; ic++) {
	    r0p[lev][ic] = r0[ic];
	    r2p[lev][ic] = r2[ic];
	}
      }
    }
  }

  final_integrate(time, dt, max_iters);

  BL_PROFILE_VAR_STOP(sdc_iters);

  // copy final solution from sdclib to new_data
  for (int lev=0; lev<=finest_level; lev++) {
    int       nnodes = mg.sweepers[lev]->nset->nnodes;
    MultiFab& Unew   = getLevel(lev).get_new_data(0);
    RNSEncap& Qend   = *((RNSEncap*) mg.sweepers[lev]->nset->Q[nnodes-1]);
    MultiFab& Uend   = *Qend.U;

    MultiFab::Copy(Unew, Uend, 0, 0, Uend.nComp(), 0);

    RNS& rns = *dynamic_cast<RNS*>(&getLevel(lev));
    rns.post_update(Unew);
  }

  for (int lev = finest_level-1; lev>= 0; lev--) {
    RNS& rns = *dynamic_cast<RNS*>(&getLevel(lev));
    rns.avgDown();
  }

  level_steps[level]++;
  level_count[level]++;
}


/*
 * Build single SDC level.
 */
sdc_sweeper* SDCAmr::build_level(int lev)
{
  BL_PROFILE("SDCAmr::build_level()");

  if (finest_level - max_trefs > 0)
    first_refinement_level = finest_level - max_trefs + 1;
  else
    first_refinement_level = 1;

  int nnodes = nnodes0;
  for (int l=first_refinement_level; l<=lev; l++)
      nnodes = (nnodes-1)*trat + 1;

  double nodes[3] = { 0.0, 0.5, 1.0 };
  int nrepeat     = (nnodes-1)/2;
  sdc_imex* imex = sdc_imex_create(nodes, 3, nrepeat, 0,
				   sdc_f1eval, sdc_f2eval, sdc_f2comp);

  sdc_imex_setup(imex, NULL, NULL);
  sdc_hooks_add(imex->hooks, SDC_HOOK_POST_STEP, sdc_poststep_hook);

  if (lev < first_refinement_level)
    mg.nsweeps[lev] = nsweeps[0];
  else
    mg.nsweeps[lev] = nsweeps[lev-first_refinement_level+1];

  return (sdc_sweeper*) imex;
}

/*
 * Rebuild MLSDC hierarchy.
 */
void SDCAmr::rebuild_mlsdc()
{
  BL_PROFILE("SDCAmr::rebuild_mlsdc()");

  // reset previous and clear sweepers etc
  sdc_mg_reset(&mg);
  destroy_mlsdc();

  // rebuild
  for (int lev=0; lev<=finest_level; lev++) {
    encaps[lev] = build_encap(lev);
    sweepers[lev] = build_level(lev);
    sweepers[lev]->nset->ctx = &getLevel(lev);
    sweepers[lev]->nset->encap = encaps[lev];
    sdc_mg_add_level(&mg, sweepers[lev], mlsdc_amr_interpolate, mlsdc_amr_restrict);
  }

  sdc_mg_setup(&mg, 0);
  sdc_mg_allocate(&mg);

  if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
    cout << "Rebuilt MLSDC: " << mg.nlevels << ", nnodes: ";
    for (int lev=0; lev<=finest_level; lev++)
      cout << sweepers[lev]->nset->nnodes << " ";
    cout << endl;
    if (verbose > 2)
      sdc_mg_print(&mg, 2);
  }
}

/*
 * Initialize SDC multigrid sweeper, set parameters.
 */
SDCAmr::SDCAmr ()
{
  ParmParse ppsdc("mlsdc");
  if (!ppsdc.query("max_iters", max_iters)) max_iters = 4;
  if (!ppsdc.query("max_trefs", max_trefs)) max_trefs = 2;
  if (!ppsdc.query("nnodes0",   nnodes0))   nnodes0 = 3;
  if (!ppsdc.query("trat",      trat))      trat = 2;

  if (verbose > 2)
    sdc_log_set_stdout(SDC_LOG_DEBUG);
  else if (verbose > 1)
    sdc_log_set_stdout(SDC_LOG_INFO);

  sdc_mg_build(&mg, max_level+1);
  sdc_hooks_add(mg.hooks, SDC_HOOK_POST_TRANS, sdc_poststep_hook);

  sweepers.resize(max_level+1);
  encaps.resize(max_level+1);
  nsweeps.resize(max_trefs+1);

  if (!ppsdc.queryarr("nsweeps", nsweeps)) {
    nsweeps[0] = 2;
    for (int l=1; l<max_trefs+1; l++)
      nsweeps[l] = 1;
  }


  for (int i=0; i<=max_level; i++)
    sweepers[i] = NULL;

  if (max_level > 0)
    for (int i=0; i<=max_level; i++)
      if (blockingFactor(i) < 4)
	BoxLib::Abort("For AMR runs, set blocking_factor to at least 4.");
}

SDCAmr::~SDCAmr()
{
  destroy_mlsdc();
  sdc_mg_destroy(&mg);
}
