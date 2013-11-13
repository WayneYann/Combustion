/*
 * Multilevel SDC + AMR controller.
 */

#include <SDCAmr.H>
#include <MultiFab.H>
#include <ParmParse.H>
#include <StateDescriptor.H>
#include <AmrLevel.H>
#include <Interpolater.H>
#include <FabArray.H>

#include "RNS.H"
#include "RNS_F.H"

#ifdef BL_USE_ARRAYVIEW
#include <ArrayView.H>
#endif

using namespace std;

BEGIN_EXTERN_C

/*
 * Spatial interpolation between MultiFabs.  Called by SDCLib.
 */
void mlsdc_amr_interpolate(void *F, void *G, sdc_state *state, void *ctxF, void *ctxG)
{
  MultiFab& UF     = *((MultiFab*) F);
  MultiFab& UG     = *((MultiFab*) G);
  RNS&      levelF = *((RNS*) ctxF);
  RNS&      levelG = *((RNS*) ctxG);

  const IntVect         ratio = levelG.fineRatio();
  const DescriptorList& dl    = levelF.get_desc_lst();
  const Array<BCRec>&   bcs   = dl[0].getBCs();
  const int             ncomp = dl[0].nComp();
  Interpolater&         map   = *dl[0].interp();
  Array<BCRec>          bcr(ncomp);

  RNS_SETNAN(UF);

  // make a coarse version (UC) of the fine multifab (UF)
  BoxArray crseba(UF.size());
  for (int i=0; i<crseba.size(); i++)
    crseba.set(i, map.CoarseBox(UF.fabbox(i), ratio));
  MultiFab UC(crseba, ncomp, 0);

  RNS_SETNAN(UC);
  UC.copy(UG);
  levelG.fill_boundary(UC, state->t, RNS::use_FillBoundary);
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


  levelF.fill_boundary(UF, state->t, RNS::set_PhysBoundary);
  RNS_ASSERTNONAN(UF);
}


/*
 * Spatial restriction between MultiFabs.  Called by SDCLib.
 */
void mlsdc_amr_restrict(void *F, void *G, sdc_state *state, void *ctxF, void *ctxG)
{
  MultiFab& UF      = *((MultiFab*) F);
  MultiFab& UG      = *((MultiFab*) G);
  RNS&      levelF  = *((RNS*) ctxF);
  RNS&      levelG  = *((RNS*) ctxG);

  if (state->kind == SDC_SOLUTION)
    levelF.fill_boundary(UF, state->t, RNS::use_FillBoundary);

  levelG.avgDown(UG, UF);

  if (state->kind == SDC_SOLUTION)
    levelG.fill_boundary(UG, state->t, RNS::use_FillBoundary);
}


END_EXTERN_C


void SDCAmr::timeStep(int level,
		      Real time,
		      int /* iteration */,
		      int /* niter */,
		      Real stop_time )
{
  BL_ASSERT(level == 0);
  double dt = dt_level[0];

  // build sdc hierarchy
  if (sweepers[0] == NULL) rebuild_mlsdc();

  // vector<MultiFab*> Unew, U0, Uend;
  // vector<RNS*>      rns;

  // for (int lev=0; lev<=finest_level; lev++) {
  //   Unew[lev] = &getLevel(lev).get_new_data(0);
  //   // U0.push_back(*((MultiFab*) mg.sweepers[lev]->nset->Q[0]));
  // }


  // regrid
  int lev_top = min(finest_level, max_level-1);
  for (int i=level; i<=lev_top; i++) {
    const int post_regrid_flag = 1;
    const int old_finest = finest_level;

    regrid(i,time);
    amr_level[0].computeNewDt(finest_level, sub_cycle, n_cycle, ref_ratio,
  			      dt_min, dt_level, stop_time, post_regrid_flag);

    for (int k=i; k<=finest_level; k++) level_count[k] = 0;
    if (old_finest > finest_level) lev_top = min(finest_level, max_level-1);
  }

  _time_step ++;

  // note: since we never use 'old data', any calls to fill patch
  // iterator will always use 'new data'.

  // set intial conditions and times
  for (int lev=0; lev<=finest_level; lev++) {
    MultiFab& Unew = getLevel(lev).get_new_data(0);
    MultiFab& U0   = *((MultiFab*) mg.sweepers[lev]->nset->Q[0]);
    MultiFab::Copy(U0, Unew, 0, 0, U0.nComp(), U0.nGrow());
    getLevel(lev).get_state_data(0).setTimeLevel(time+dt, dt, dt);
  }

  for (int lev=0; lev<=finest_level; lev++) {
    RNS&      rns  = *dynamic_cast<RNS*>(&getLevel(lev));
    MultiFab& U0   = *((MultiFab*) mg.sweepers[lev]->nset->Q[0]);
    rns.fill_boundary(U0, time, RNS::use_FillCoarsePatch);
  }

  // spread and iterate (XXX: spread from qend if step>0)
  if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
    cout << "MLSDC advancing with dt: " << dt << endl;
  }

  sdc_mg_spread(&mg, time, dt, 0);
  for (int k=0; k<max_iters; k++) {
    sdc_mg_sweep(&mg, time, dt, (k==max_iters-1) ? SDC_MG_LAST_SWEEP : 0);

    // echo residuals
    if (verbose > 0) {
      for (int lev=0; lev<=finest_level; lev++) {
        int       nnodes = mg.sweepers[lev]->nset->nnodes;
        MultiFab& R      = *((MultiFab*) mg.sweepers[lev]->nset->R[nnodes-2]);
	double    r0     = R.norm0();
	double    r2     = R.norm2();

	if (ParallelDescriptor::IOProcessor()) {
	  cout << "MLSDC iter: " << k << ", level: " << lev
	       << ", res norm0: " << scientific << r0 << ", res norm2: " << r2 << endl;
	}
      }
    }
  }

  // copy final solution from SDCLib to 'new data'
  for (int lev=0; lev<=finest_level; lev++) {
    int       nnodes = mg.sweepers[lev]->nset->nnodes;
    MultiFab& Unew   = getLevel(lev).get_new_data(0);
    MultiFab& Uend   = *((MultiFab*) mg.sweepers[lev]->nset->Q[nnodes-1]);
    MultiFab::Copy(Unew, Uend, 0, 0, Uend.nComp(), 0);
  }

  level_steps[level]++;
  level_count[level]++;
}

sdc_sweeper* SDCAmr::build_mlsdc_level(int lev)
{
  int first_refinement_level, nnodes;

  if (finest_level - max_trefs > 1)
    first_refinement_level = finest_level - max_trefs;
  else
    first_refinement_level = 1;

  if (lev < first_refinement_level)
    nnodes = nnodes0;
  else
    nnodes = 1 + (nnodes0 - 1) * ((int) pow(trat, lev-first_refinement_level+1));

  sdc_nodes* nodes = sdc_nodes_create(nnodes, SDC_UNIFORM);
  sdc_imex*  imex  = sdc_imex_create(nodes, sdc_f1eval, sdc_f2eval, sdc_f2comp);

  // XXX: for fine levels, need to make the integration tables etc local only

  sdc_imex_setup(imex, NULL, NULL);
  sdc_hooks_add(imex->hooks, SDC_HOOK_POST_STEP, sdc_poststep_hook);
  sdc_nodes_destroy(nodes);
  return (sdc_sweeper*) imex;
}

void SDCAmr::rebuild_mlsdc()
{
  // reset previous and clear sweepers etc
  sdc_mg_reset(&mg);
  for (unsigned int lev=0; lev<=max_level; lev++) {
    if (sweepers[lev] != NULL) {
      sweepers[lev]->destroy(sweepers[lev]); sweepers[lev] = NULL;
      destroy_encap(lev);
    }
  }

  // rebuild
  for (int lev=0; lev<=finest_level; lev++) {
    encaps[lev] = build_encap(lev);
    sweepers[lev] = build_mlsdc_level(lev);
    sweepers[lev]->nset->ctx = &getLevel(lev);
    sweepers[lev]->nset->encap = encaps[lev];
    sdc_mg_add_level(&mg, sweepers[lev], mlsdc_amr_interpolate, mlsdc_amr_restrict);
  }
  sdc_mg_setup(&mg);
  sdc_mg_allocate(&mg);

  // XXX: for fine levels, need to make the interpolation matrices local only

  if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
    cout << "Rebuilt MLSDC: " << mg.nlevels << ", nnodes: ";
    for (int lev=0; lev<=finest_level; lev++)
      cout << sweepers[lev]->nset->nnodes << " ";
    cout << endl;
  }
}

void SDCAmr::regrid(int  lbase, Real time, bool initial)
{
  Amr::regrid(lbase, time, initial);
  rebuild_mlsdc();
}

SDCAmr::SDCAmr ()
{
  ParmParse ppsdc("mlsdc");
  if (!ppsdc.query("max_iters", max_iters)) max_iters = 8;
  if (!ppsdc.query("max_trefs", max_trefs)) max_trefs = 3;
  if (!ppsdc.query("nnodes0",   nnodes0))   nnodes0 = 3;
  if (!ppsdc.query("trat",      trat))      trat = 2;

  // sdc_log_set_stdout(SDC_LOG_DEBUG);
  sdc_mg_build(&mg, max_level+1);
  sdc_hooks_add(mg.hooks, SDC_HOOK_POST_TRANS, sdc_poststep_hook);

  sweepers.resize(max_level+1);
  encaps.resize(max_level+1);

  for (int i=0; i<=max_level; i++) sweepers[i] = NULL;

  if (max_level > 0) {
      for (int i=0; i<=max_level; i++) {
	  if (blockingFactor(i) < 4) {
	      BoxLib::Abort("For AMR runs, set blocking_factor to at least 4.");
	  }
      }
  }
}

SDCAmr::~SDCAmr()
{
  sdc_mg_destroy(&mg);
}
