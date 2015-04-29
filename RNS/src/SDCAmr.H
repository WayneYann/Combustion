
#ifdef USE_SDCLIB
#ifndef _SDCAMR_H_
#define _SDCAMR_H_

#include <MLSDCAmr.H>
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
 * RNS MLSDC+AMR class.
 */

class SDCAmr : public MLSDCAmr
{
  int ho_imex;
public:
  SDCAmr();
  virtual sdc_sweeper* BuildLevel(int lev);
  virtual void coarseTimeStep (Real stop_time);
};

#endif
#endif
