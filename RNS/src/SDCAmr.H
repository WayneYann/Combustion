
#ifdef USE_SDCLIB
#ifndef _SDCAMR_H_
#define _SDCAMR_H_

#include <MultiFab.H>
#include <FluxRegister.H>
#include <Amr.H>
#include <RNS.H>
#include <sdc.h>

#if defined(NDEBUG)
#define RNS_SETNAN(MF)
#define RNS_ASSERTNONAN(MF)
#else
#define RNS_SETNAN(MF) MF.setVal(NAN, MF.nGrow());
#define RNS_ASSERTNONAN(MF) BL_ASSERT(MF.contains_nan() == false);
#endif

BEGIN_EXTERN_C
void sdc_f1eval(void*, void*, double, sdc_state*, void*);
void sdc_f2eval(void*, void*, double, sdc_state*, void*);
void sdc_f2comp(void*, void*, double, double, void*, sdc_state*, void*);
void sdc_poststep_hook(void*, sdc_state*, void*);
END_EXTERN_C

using namespace std;

/*
 * MultiFab encapsulation context for SDCLib.
 */

typedef struct {
  const BoxArray *ba;
  IntVect crse_ratio;
  int ncomp, ngrow, level, finest;
  RNS *rns;
} RNSEncapCtx;

typedef struct {
  MultiFab* U;			// solution or function eval
  FluxRegister* fine_flux;	// current levels flux register
  FluxRegister* crse_flux;	// coarse levels flux register
  RNS* rns;
  int type;
} RNSEncap;

/*
 * Multigrid SDC+AMR class.
 */

class SDCAmr : public Amr
{
public:

  int max_iters;                // maximum number of sdc iterations
  int max_trefs;                // maximum number of time refinements
  int nnodes0;			// number of nodes on the coarsest level
  int trat;			// time refinement ratio
  int ho_imex;
  int first_refinement_level;   // first amr level that has time refinement

  sdc_mg mg;
  vector<sdc_sweeper*> sweepers;
  vector<sdc_encap*>   encaps;
  vector<int>          nsweeps;

  SDCAmr();
  ~SDCAmr();

  sdc_encap*   build_encap(int lev);
  sdc_sweeper* build_level(int lev);
  void         rebuild_mlsdc();
  void         destroy_mlsdc();
  void         final_integrate(double t, double dt, int niter);

  virtual void init (Real strt_time, Real stop_time);

  virtual void timeStep(int level,
			Real time,
			int iteration,
			int niter ,
			Real stop_time);

};

void dgp_send_mf(MultiFab& U, int level, int comp, int wait);

#endif
#endif
