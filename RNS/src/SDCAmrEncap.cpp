/*
 * RNS encapsulation for SDCLib.
 *
 * Notes:
 *   - State/solution encaps are created with grow/ghost cells.
 *   - Function evaluation encaps are created without grow/ghost cells.
 *   - Integral encaps are created without grow/ghost cells.
 */

#include <MultiFab.H>
#include <SDCAmr.H>
#include <StateDescriptor.H>
#include <AmrLevel.H>

#include <cassert>

BEGIN_EXTERN_C

void mf_encap_setval(void *Qptr, sdc_dtype val);


void *mf_encap_create(int type, void *encap_ctx)
{
  RNSEncapCtx* ctx   = (RNSEncapCtx*) encap_ctx;
  RNSEncap*    encap = new RNSEncap;

  BoxArray ba(*ctx->ba);

  encap->flux = 0;
  encap->type = type;
  encap->rns  = ctx->rns;

  switch (type) {
  case SDC_SOLUTION:
  case SDC_WORK:
    encap->U = new MultiFab(ba, ctx->ncomp, ctx->ngrow);
    mf_encap_setval(encap, 0.0);
    break;
  case SDC_FEVAL:
    encap->U = new MultiFab(ba, ctx->ncomp, 0);
    if (ctx->level > 0)
      encap->flux = new FluxRegister(ba, ctx->crse_ratio, ctx->level, ctx->ncomp);
    mf_encap_setval(encap, 0.0);
    break;
  case SDC_INTEGRAL:
  case SDC_TAU:
    encap->U    = new MultiFab(ba, ctx->ncomp, 0);
    if (ctx->level > 0)
      encap->flux = new FluxRegister(ba, ctx->crse_ratio, ctx->level, ctx->ncomp);
    mf_encap_setval(encap, 0.0);
    break;
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

  if (Q.flux) {
    FluxRegister& F = *Q.flux;
    for (OrientationIter face; face; ++face)
      for (FabSetIter bfsi(F[face()]); bfsi.isValid(); ++bfsi)
        F[face()][bfsi].setVal(val);
  }
}

void mf_encap_copy(void *dstp, const void *srcp)
{
  RNSEncap& Qdst = *((RNSEncap*) dstp);
  RNSEncap& Qsrc = *((RNSEncap*) srcp);
  MultiFab& Udst = *Qdst.U;
  MultiFab& Usrc = *Qsrc.U;

  for (MFIter mfi(Udst); mfi.isValid(); ++mfi)
    Udst[mfi].copy(Usrc[mfi]);

  // MultiFab::Copy(Udst, Usrc, 0, 0, Udst.nComp(), Udst.nGrow());

  if (Qdst.flux && Qsrc.flux) {
    // XXX: should test this
    FluxRegister& Fdst = *Qdst.flux;
    FluxRegister& Fsrc = *Qsrc.flux;
    for (OrientationIter face; face; ++face)
      for (FabSetIter bfsi(Fdst[face()]); bfsi.isValid(); ++bfsi)
        Fdst[face()][bfsi].copy(Fsrc[face()][bfsi]);
  }

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

  if ((Qy.type==SDC_TAU) && (Qx.flux!=NULL)) {
    FluxRegister& Fx  = *Qx.flux;
    FluxRegister& Fy  = *Qy.flux;
    for (OrientationIter face; face; ++face)
      for (FabSetIter bfsi(Fy[face()]); bfsi.isValid(); ++bfsi) {
        Fy[face()][bfsi].saxpy(a, Fx[face()][bfsi]);
	// cout << "SAXPY ";
        // cout << Fy[face()][bfsi].norm(0, 0, Uy.nComp()) << " ";
        // cout << Fx[face()][bfsi].norm(0, 0, Ux.nComp()) << endl;
      }
  }
}

END_EXTERN_C


sdc_encap* SDCAmr::build_encap(int lev)
{
  const DescriptorList& dl = getLevel(lev).get_desc_lst();
  assert(dl.size() == 1);       // valid for RNS

  RNSEncapCtx* ctx = new RNSEncapCtx;
  ctx->level  = lev;
  ctx->ba     = &boxArray(lev);
  ctx->rns    = dynamic_cast<RNS*>(&getLevel(lev));
  ctx->finest = lev == finest_level;
  ctx->ncomp  = dl[0].nComp();
  ctx->ngrow  = dl[0].nExtra();
  if (lev > 0)
    ctx->crse_ratio = refRatio(lev-1);
  if (lev < finest_level)
    ctx->fine_ba = &boxArray(lev+1);

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
