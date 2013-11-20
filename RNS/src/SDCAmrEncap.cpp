
#include <MultiFab.H>
#include <SDCAmr.H>
#include <StateDescriptor.H>
#include <AmrLevel.H>

#include <cassert>

BEGIN_EXTERN_C

void *mf_encap_create(int type, void *encap_ctx)
{
  RNSEncapCtx* ctx   = (RNSEncapCtx*) encap_ctx;
  RNSEncap*    encap = new RNSEncap;
  encap->flux = 0;
  if (type == SDC_SOLUTION || type == SDC_WORK) {
    encap->U    = new MultiFab(*ctx->ba, ctx->ncomp, ctx->ngrow);
  } else {
    encap->U    = new MultiFab(*ctx->ba, ctx->ncomp, 0);
    if (ctx->level > 0)
      encap->flux = new FluxRegister(*ctx->ba, ctx->rr, ctx->level, ctx->ncomp);
  }
  return encap;
}

void mf_encap_destroy(void *Qptr)
{
  RNSEncap* Q = (RNSEncap*) Qptr;
  delete Q->U;
  if (Q->flux != NULL)
    delete Q->flux;
  delete Q;
}

void mf_encap_setval(void *Qptr, sdc_dtype val)
{
  RNSEncap& Q = *((RNSEncap*) Qptr);
  MultiFab& U = *Q.U;
  U.setVal(val, U.nGrow());
}

void mf_encap_copy(void *dstp, const void *srcp)
{
  RNSEncap& Qdst = *((RNSEncap*) dstp);
  RNSEncap& Qsrc = *((RNSEncap*) srcp);
  MultiFab& Udst = *Qdst.U;
  MultiFab& Usrc = *Qsrc.U;
  MultiFab::Copy(Udst, Usrc, 0, 0, Udst.nComp(), Udst.nGrow());
#ifndef NDEBUG
  BL_ASSERT(Usrc.contains_nan() == false);
  BL_ASSERT(Udst.contains_nan() == false);
#endif
}

void mf_encap_saxpy(void *yp, sdc_dtype a, void *xp)
{
  RNSEncap& Qy = *((RNSEncap*) yp);
  RNSEncap& Qx = *((RNSEncap*) xp);
  MultiFab& Uy = *Qy.U;
  MultiFab& Ux = *Qx.U;
// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
  BL_ASSERT(Uy.boxArray() == Ux.boxArray());

  for (MFIter mfi(Uy); mfi.isValid(); ++mfi)
    Uy[mfi].saxpy(a, Ux[mfi]);
}

END_EXTERN_C


sdc_encap* SDCAmr::build_encap(int lev)
{
  const DescriptorList& dl = getLevel(lev).get_desc_lst();
  assert(dl.size() == 1);	// valid for RNS

  RNSEncapCtx* ctx = new RNSEncapCtx;
  ctx->ba    = &boxArray(lev);
  ctx->rr    = refRatio(lev);
  ctx->ncomp = dl[0].nComp();
  ctx->ngrow = dl[0].nExtra();
  ctx->level = lev;

  sdc_encap* encap = new sdc_encap;
  encap->create  = mf_encap_create;
  encap->destroy = mf_encap_destroy;
  encap->setval  = mf_encap_setval;
  encap->copy    = mf_encap_copy;
  encap->saxpy   = mf_encap_saxpy;
  encap->ctx     = ctx;

  return encap;
}

void SDCAmr::destroy_encap(int lev)
{
  // XXX: memory leaks here?
  delete (RNSEncapCtx*) encaps[lev]->ctx;
  delete encaps[lev];
}

