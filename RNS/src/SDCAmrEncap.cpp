/*
 * RNS encapsulation for SDCLib.
 *
 * Notes:
 *   - State/solution encaps are created with grow/ghost cells.
 *   - Function evaluation encaps are created without grow/ghost cells.
 *   - Integral encaps are created without grow/ghost cells.
 *
 * XXX: Note that the FEVAL encapsulations have flux registers, and
 * since we're using the IMEX sweeper, both the "explicit" feval and
 * "implicit" feval will have flux registers, but this isn't
 * necessary.  Matt should clean this up sometime.
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

  encap->fine_flux = 0;
  encap->crse_flux = 0;
  encap->type = type;
  encap->rns  = ctx->rns;

  switch (type) {
  case SDC_SOLUTION:
  case SDC_WORK:
    encap->U = new MultiFab(ba, ctx->ncomp, ctx->ngrow);
    mf_encap_setval(encap, 0.0);
    break;
  case SDC_FEVAL:
  case SDC_INTEGRAL:
  case SDC_TAU:
    encap->U = new MultiFab(ba, ctx->ncomp, 0);
    if (ctx->level > 0)
      encap->fine_flux = new FluxRegister(ba, ctx->crse_ratio, ctx->level, ctx->ncomp);
    if (! ctx->finest) {
      SDCAmr&   amr  = *encap->rns->getSDCAmr();
      AmrLevel& rnsF = amr.getLevel(ctx->level+1);
      encap->crse_flux = new FluxRegister(rnsF.boxArray(), amr.refRatio(ctx->level), rnsF.Level(), ctx->ncomp);
    }
    mf_encap_setval(encap, 0.0);
    break;
  }

  return encap;
}

void mf_encap_destroy(void *Qptr)
{
  RNSEncap* Q = (RNSEncap*) Qptr;
  delete Q->U;
  if (Q->fine_flux != NULL) delete Q->fine_flux;
  if (Q->crse_flux != NULL) delete Q->crse_flux;
  delete Q;
}

void mf_encap_setval_flux(FluxRegister& dst, sdc_dtype val)
{
  for (OrientationIter face; face; ++face)
    for (FabSetIter bfsi(dst[face()]); bfsi.isValid(); ++bfsi)
      dst[face()][bfsi].setVal(val);
}

void mf_encap_setval(void *Qptr, sdc_dtype val)
{
  RNSEncap& Q = *((RNSEncap*) Qptr);
  MultiFab& U = *Q.U;
  U.setVal(val, U.nGrow());

  if (Q.fine_flux) mf_encap_setval_flux(*Q.fine_flux, val);
  if (Q.crse_flux) mf_encap_setval_flux(*Q.crse_flux, val);
}

void mf_encap_copy_flux(FluxRegister& dst, FluxRegister& src)
{
  for (OrientationIter face; face; ++face)
    for (FabSetIter bfsi(dst[face()]); bfsi.isValid(); ++bfsi)
      dst[face()][bfsi].copy(src[face()][bfsi]);
}

void mf_encap_copy(void *dstp, const void *srcp, int flags)
{
  RNSEncap& Qdst = *((RNSEncap*) dstp);
  RNSEncap& Qsrc = *((RNSEncap*) srcp);
  MultiFab& Udst = *Qdst.U;
  MultiFab& Usrc = *Qsrc.U;

  int nghost = 0;
  if (flags & SDC_COPY_GHOST) {
      int ngsrc = Usrc.nGrow();
      int ngdst = Udst.nGrow();
      nghost = (ngdst < ngsrc) ? ngdst : ngsrc;
  }
  MultiFab::Copy(Udst, Usrc, 0, 0, Usrc.nComp(), nghost);

  if (Qdst.fine_flux && Qsrc.fine_flux) mf_encap_copy_flux(*Qdst.fine_flux, *Qsrc.fine_flux);
  if (Qdst.crse_flux && Qsrc.crse_flux) mf_encap_copy_flux(*Qdst.crse_flux, *Qsrc.crse_flux);

#ifndef NDEBUG
  BL_ASSERT(Usrc.contains_nan() == false);
  BL_ASSERT(Udst.contains_nan() == false);
#endif
}

void mf_encap_saxpy_flux(FluxRegister& y, sdc_dtype a, FluxRegister& x)
{
  for (OrientationIter face; face; ++face)
    for (FabSetIter bfsi(y[face()]); bfsi.isValid(); ++bfsi)
      y[face()][bfsi].saxpy(a, x[face()][bfsi]);
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

  if ((Qy.type==SDC_TAU) && (Qx.fine_flux!=NULL)) mf_encap_saxpy_flux(*Qy.fine_flux, a, *Qx.fine_flux);
  if ((Qy.type==SDC_TAU) && (Qx.crse_flux!=NULL)) mf_encap_saxpy_flux(*Qy.crse_flux, a, *Qx.crse_flux);
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

  sdc_encap* encap = new sdc_encap;
  encap->create  = mf_encap_create;
  encap->destroy = mf_encap_destroy;
  encap->setval  = mf_encap_setval;
  encap->copy    = mf_encap_copy;
  encap->saxpy   = mf_encap_saxpy;
  encap->ctx     = ctx;

  return encap;
}

void SDCAmr::destroy_mlsdc()
{
  for (unsigned int lev=0; lev<=max_level; lev++) {
    if (sweepers[lev] != NULL) {
      sweepers[lev]->destroy(sweepers[lev]);
      sweepers[lev] = NULL;
      delete (RNSEncapCtx*) encaps[lev]->ctx;
      delete encaps[lev];
    }
  }
}
