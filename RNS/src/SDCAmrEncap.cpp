
#include <MultiFab.H>
#include <SDCAmr.H>
#include <StateDescriptor.H>
#include <AmrLevel.H>

#include <cassert>

BEGIN_EXTERN_C

void *mf_encap_create(int type, void *encap_ctx)
{
  mf_encap* ctx = (mf_encap*) encap_ctx;
  if (type == SDC_SOLUTION || type == SDC_WORK) {
    return new MultiFab(*ctx->ba, ctx->ncomp, ctx->ngrow);
  } else {
    return new MultiFab(*ctx->ba, ctx->ncomp, 0);
  }
}

void mf_encap_destroy(void *sol)
{
  MultiFab* mf = (MultiFab*) sol;
  delete mf;
}

void mf_encap_setval(void *sol, sdc_dtype val)
{
  MultiFab& mf = *((MultiFab*) sol);
  mf.setVal(val, mf.nGrow());
}

void mf_encap_copy(void *dstp, const void *srcp)
{
  MultiFab& dst = *((MultiFab*) dstp);
  MultiFab& src = *((MultiFab*) srcp);
  for (MFIter mfi(dst); mfi.isValid(); ++mfi)
    dst[mfi].copy(src[mfi]);
}

void mf_encap_saxpy(void *yp, sdc_dtype a, void *xp)
{
  MultiFab& y = *((MultiFab*) yp);
  MultiFab& x = *((MultiFab*) xp);
  for (MFIter mfi(y); mfi.isValid(); ++mfi)
    y[mfi].saxpy(a, x[mfi]);
}

END_EXTERN_C


sdc_encap* SDCAmr::build_encap(int lev)
{
  const DescriptorList& dl = getLevel(lev).get_desc_lst();

  assert(dl.size() == 1);	// valid for RNS

  mf_encap* ctx = new mf_encap;
  ctx->ba    = &boxArray(lev);
  ctx->ncomp = dl[0].nComp();
  ctx->ngrow = dl[0].nExtra();

  sdc_encap* encap = new sdc_encap;
  encap->create  = mf_encap_create;
  encap->destroy = mf_encap_destroy;
  encap->setval  = mf_encap_setval;
  encap->copy    = mf_encap_copy;
  encap->saxpy   = mf_encap_saxpy;
  encap->ctx     = ctx;

  return encap;
}
