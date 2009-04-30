#include <winstd.H>

#include "DDOp.H"
#include "DDOp_F.H"
#include "LO_F.H"

// Note: The ratio between the original grids and the coarse-fine boundary data
//       registers can be anything, but the ratio between adjacent DD operators
//       generated for multigrid will always be two.
const IntVect MGIV = IntVect(D_DECL(2,2,2));

DDOp::DDOp (const ChemDriver& ckd)
    : ckdriver(ckd), coarser(0)
{}

DDOp::DDOp (const BoxArray&   grids,
            const Box&        box,
            const ChemDriver& ckd,
            const IntVect&    ratio,
            int               mgLevel)
    : ckdriver(ckd), coarser(0)
{
    define(grids,box,ratio,mgLevel);
}

DDOp::~DDOp ()
{
    if (coarser)
        delete coarser;
}

bool can_coarsen(const BoxArray& ba)
{
    // Ratio between levels here will always be MGIV
    for (int i = 0; i < ba.size(); ++i)
    {
        Box tmp = ba[i];
        Box ctmp  = BoxLib::coarsen(ba[i],MGIV);
        Box rctmp = BoxLib::refine(ctmp,MGIV);
        if (tmp != rctmp || ctmp.numPts() == 1)
            return false;
    }
    return true;
}

bool
DDOp::coarser_exists(int level) const
{
    BL_ASSERT(level >= 0);
    if (level == 0)
    {
        if (coarser == 0)
            return false;
        return true;
    }
    return coarser->coarser_exists(--level);
}

void
DDOp::define (const BoxArray& _grids,
              const Box&      box,
              const IntVect&  ratio,
              int             mgLevel)
{
    grids = _grids;
    Geometry geom(box);
    Tbd.define(grids,1,geom,mgLevel);
    Ybd.define(grids,ckdriver.numSpecies(),geom,mgLevel);
    cfRatio = ratio;
    const int gGrow = 0;
    geom.GetVolume(volume,grids,gGrow);
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
        geom.GetFaceArea(area[dir],grids,dir,gGrow);

    // Generate coarser one (ratio = MGIV), if possible
    if (can_coarsen(grids))
    {
        const BoxArray cGrids = BoxArray(grids).coarsen(MGIV);
        const Box cBox = Box(box).coarsen(MGIV);
        BL_ASSERT(Box(cBox).refine(MGIV) == box);
        coarser = new DDOp(cGrids,cBox,ckdriver,cfRatio,mgLevel+1);
    }
}

void
DDOp::center_to_edge (const FArrayBox& cfab,
                      FArrayBox&       efab,
                      const Box&       ccBox,
                      int              sComp,
                      int              dComp,
                      int              nComp) const
{
    // Compute data on edges between each pair of cc values in the dir direction
    const Box&      ebox = efab.box();
    const IndexType ixt  = ebox.ixType();

    BL_ASSERT(!(ixt.cellCentered()) && !(ixt.nodeCentered()));

    int dir = -1;
    for (int d = 0; d < BL_SPACEDIM; d++)
        if (ixt.test(d))
            dir = d;

    BL_ASSERT(BoxLib::grow(ccBox,-BoxLib::BASISV(dir)).contains(BoxLib::enclosedCells(ebox)));
    BL_ASSERT(sComp+nComp <= cfab.nComp() && dComp+nComp <= efab.nComp());

    FORT_DDC2E(ccBox.loVect(),ccBox.hiVect(),
               ARLIM(cfab.loVect()),ARLIM(cfab.hiVect()),cfab.dataPtr(sComp),
               ARLIM(efab.loVect()),ARLIM(efab.hiVect()),efab.dataPtr(dComp),
               &nComp, &dir);
}    

void
coarsenBndryData(const MultiFab&      fmf,
                 int                  fmf_sc,
                 const BndryRegister* fbr,
                 int                  fbr_sc,
                 MultiFab&            cmf,
                 int                  cmf_sc,
                 BndryRegister*       cbr,
                 int                  cbr_sc,
                 int                  nc)
{
    // coarsen boundary data, always by MGIV
    BL_ASSERT(fmf.boxArray() == BoxArray(cmf.boxArray()).refine(MGIV));
    for (MFIter mfi(fmf); mfi.isValid(); ++mfi)
    {
        const FArrayBox& ffab = fmf[mfi];
        FArrayBox& cfab = cmf[mfi];
        const Box& gbox = mfi.validbox();
        for (OrientationIter oitr; oitr; ++oitr)
        {
            const Box fbox = BoxLib::adjCell(gbox,oitr(),1);
            const int face = (int)oitr();
            FORT_CRSNCCBND(fbox.loVect(), fbox.hiVect(),
                           ffab.dataPtr(fmf_sc), ARLIM(ffab.loVect()), ARLIM(ffab.hiVect()),
                           cfab.dataPtr(cmf_sc), ARLIM(cfab.loVect()), ARLIM(cfab.hiVect()),
                           &nc, &face, MGIV.getVect());
        }
    }

    if (fbr && cbr)
    {
        const BoxArray& fgrids = fbr->boxes();
        const BoxArray& cgrids = cbr->boxes();
        BL_ASSERT(fgrids == BoxArray(cgrids).refine(MGIV));
        for (OrientationIter oitr; oitr; ++oitr)
        {
            const FabSet& ffs = (*fbr)[oitr()];
            for (FabSetIter ffsi(ffs); ffsi.isValid(); ++ffsi)
            {
                const Box fbox = BoxLib::adjCell(fgrids[ffsi.index()],oitr(),1);
                const int face = (int)oitr();
                const FArrayBox& ffab = ffs[ffsi];
                FArrayBox& cfab = (*cbr)[oitr()][ffsi];
                FORT_CRSNCCBND(fbox.loVect(), fbox.hiVect(),
                               ffab.dataPtr(fbr_sc), ARLIM(ffab.loVect()), ARLIM(ffab.hiVect()),
                               cfab.dataPtr(cbr_sc), ARLIM(cfab.loVect()), ARLIM(cfab.hiVect()),
                               &nc, &face, MGIV.getVect());
            }
        }
    }
}                 

void
DDOp::setBoundaryData(const MultiFab&      fineT,
                      int                  fStartT,
                      const MultiFab&      fineY,
                      int                  fStartY,
                      const BndryRegister* cbrT,
                      int                  cStartT,
                      const BndryRegister* cbrY,
                      int                  cStartY,
                      const BCRec&         bcT,
                      const BCRec&         bcY)
{
    BL_ASSERT(fineT.boxArray() == grids);
    BL_ASSERT(fineY.boxArray() == grids);
    const int Nspec = ckdriver.numSpecies();

    if (coarser)  // if so, then it will be MGIV coarser
    {
        const int newTmfComp = Nspec;
        const int newTbrComp = Nspec;
        const int newYmfComp = 0;
        const int newYbrComp = 0;
        const BoxArray cBA = BoxArray(grids).coarsen(MGIV);
        MultiFab newMF(cBA,Nspec+1,1);
#ifndef NDEBUG
        newMF.setVal(-1.0);
#endif
        BoxArray ccBA = BoxArray(cBA).coarsen(MGIV);
        BndryRegister* newBR = 0;
        if (cbrT && cbrY && coarser->coarser)
        {
            newBR = new BndryRegister();
            newBR->setBoxes(ccBA);
            for (OrientationIter fi; fi; ++fi)
                newBR->define(fi(),IndexType::TheCellType(),0,1,1,Nspec+1);
#ifndef NDEBUG
            newBR->setVal(-1.0);
#endif
        }
        coarsenBndryData(fineT,fStartT,cbrT,cStartT,
                         newMF,newTmfComp,newBR,newTbrComp,1);
        coarsenBndryData(fineY,fStartY,cbrY,cStartY,
                         newMF,newYmfComp,newBR,newYbrComp,Nspec);
        coarser->setBoundaryData(newMF,newTmfComp,newMF,newYmfComp,
                                 newBR,newTbrComp,newBR,newYbrComp,bcT,bcY);
    }

    IntVect ratio(cfRatio); // To avoid const problems
    if (cbrT == 0)
    {
        Tbd.setBndryValues(fineT,fStartT,0,1,bcT);
    }
    else
    {
        Tbd.setBndryValues(const_cast<BndryRegister&>(*cbrT),
                           cStartT,fineT,fStartT,0,1,ratio,bcT);
    }

    if (cbrY == 0)
    {
        Ybd.setBndryValues(fineY,fStartY,0,Nspec,bcY);
    }
    else
    {
        Ybd.setBndryValues(const_cast<BndryRegister&>(*cbrY),
                           cStartY,fineY,fStartY,0,Nspec,ratio,bcY);
    }
}

void
DDOp::setGrowCells(MultiFab& T,
                   int       compT,
                   MultiFab& Y,
                   int       compY) const
{
    BL_ASSERT(T.nGrow() >= 1);
    BL_ASSERT(Y.nGrow() >= 1);
    const int Nspec = ckdriver.numSpecies();

    BL_ASSERT(T.nComp() > compT);
    BL_ASSERT(Y.nComp() >= compY + Nspec);
    BL_ASSERT(&T != &Y || (compT<compY || compT >= compY+Nspec));

    BL_ASSERT(T.boxArray() == grids);
    BL_ASSERT(Y.boxArray() == grids);

    T.FillBoundary(compT,1);
    Tbd.getGeom().FillPeriodicBoundary(T,compT,1);
    Y.FillBoundary(compY,Nspec);
    Ybd.getGeom().FillPeriodicBoundary(Y,compY,Nspec);

    const int flagbc  = 1;
    const int flagden = 0; // Use LinOp's bc interpolator, but don't save the coeff
    const int maxorder = 3;
    Real* dummy = 0;
    Box dumbox(IntVect(D_DECL(0,0,0)),IntVect(D_DECL(0,0,0)));
    const Real* dx = Tbd.getGeom().CellSize();
    for (OrientationIter oitr; oitr; ++oitr)
    {
        const Orientation&      face = oitr();
        const int              iFace = (int)face;

        const Array<Array<BoundCond> >& Tbc = Tbd.bndryConds(face);
        const Array<Real>&      Tloc = Tbd.bndryLocs(face);
        const FabSet&            Tfs = Tbd.bndryValues(face);
        const int                Tnc = 1;
        
        const Array<Array<BoundCond> >&  Ybc = Ybd.bndryConds(face);
        const Array<Real>&      Yloc = Ybd.bndryLocs(face);
        const FabSet&            Yfs = Ybd.bndryValues(face);
        const int                Ync = Nspec;

        const int comp = 0;
        for (MFIter mfi(T); mfi.isValid(); ++mfi)
        {
            const int   idx = mfi.index();
            const Box& vbox = mfi.validbox();
            BL_ASSERT(grids[idx] == vbox);

            FArrayBox& Tfab = T[mfi];
            const FArrayBox& Tb = Tfs[mfi];

            const Mask& Tm  = Tbd.bndryMasks(face)[idx];
            const Real Tbcl = Tloc[idx];
            const int Tbct  = Tbc[idx][comp];

            FORT_APPLYBC(&flagden, &flagbc, &maxorder,
                         Tfab.dataPtr(compT), ARLIM(Tfab.loVect()), ARLIM(Tfab.hiVect()),
                         &iFace, &Tbct, &Tbcl,
                         Tb.dataPtr(), ARLIM(Tb.loVect()), ARLIM(Tb.hiVect()),
                         Tm.dataPtr(), ARLIM(Tm.loVect()), ARLIM(Tm.hiVect()),
                         dummy, ARLIM(dumbox.loVect()), ARLIM(dumbox.hiVect()),
                         vbox.loVect(),vbox.hiVect(), &Tnc, dx);

            FArrayBox& Yfab = Y[mfi];
            const FArrayBox& Yb = Yfs[mfi];

            const Mask& Ym  = Ybd.bndryMasks(face)[idx];
            const Real Ybcl = Yloc[idx];
            const int Ybct  = Ybc[idx][comp];

            FORT_APPLYBC(&flagden, &flagbc, &maxorder,
                         Yfab.dataPtr(compY), ARLIM(Yfab.loVect()), ARLIM(Yfab.hiVect()),
                         &iFace, &Ybct, &Ybcl,
                         Yb.dataPtr(), ARLIM(Yb.loVect()), ARLIM(Yb.hiVect()),
                         Ym.dataPtr(), ARLIM(Ym.loVect()), ARLIM(Ym.hiVect()),
                         dummy, ARLIM(dumbox.loVect()), ARLIM(dumbox.hiVect()),
                         vbox.loVect(),vbox.hiVect(), &Ync, dx);
        }
    }
}

void
DDOp::applyOp(MultiFab&         outH,
              int               dCompH,
              MultiFab&         outY,
              int               dCompY,
              const MultiFab&   inT,
              int               sCompT,
              const MultiFab&   inY,
              int               sCompY,
              PArray<MultiFab>& fluxH,
              int               dCompFH,
              PArray<MultiFab>& fluxY,
              int               dCompFY,
              DD_ApForTorRH     whichApp,
              int               level) const
{
    BL_ASSERT(level >= 0);
    if (level > 0)
    {
        coarser->applyOp(outH,dCompH,outY,dCompY,inT,sCompT,inY,sCompY,
                         fluxH,dCompFH,fluxY,dCompFY,whichApp,level-1);
        return;
    }

    const int Nspec = ckdriver.numSpecies();
    const int nGrow = 1;

    BL_ASSERT(outH.nComp() > dCompH);
    BL_ASSERT(outY.nComp() >= dCompY + Nspec);
    BL_ASSERT(&inT != &inY || (sCompT<sCompY || sCompT >= sCompY+Nspec));
    BL_ASSERT(&outH != &outY || (dCompH<dCompY || dCompH >= dCompY+Nspec));
    BL_ASSERT(outH.boxArray() == grids);
    BL_ASSERT(outY.boxArray() == grids);
    BL_ASSERT(inT.boxArray() == grids);
    BL_ASSERT(inY.boxArray() == grids);
    BL_ASSERT(inT.nGrow() >= nGrow);
    BL_ASSERT(inY.nGrow() >= nGrow);

    // Need grow cells in X,T to compute forcing, and to get Ye for evaluating coeffs
    // Promise to change only the grow cells in T,Y
    setGrowCells(const_cast<MultiFab&>(inT),sCompT,const_cast<MultiFab&>(inY),sCompY);

    // Get diffusion flux divergence
    const Real* dx = Tbd.getGeom().CellSize();
    FArrayBox de, et, Xc, TYe, Hie;
    const int sCompde = 0;
    const int sCompXc = 0;
    const int sCompTe = Nspec;
    const int sCompYe = 0;
    const int do_add_enth_flux = (whichApp==DD_RhoH ? 1 : 0);

    // Allocate/initialize output multifabs
    outH.setVal(0.0,dCompH,1);
    outY.setVal(0.0,dCompY,Nspec);

    for (MFIter mfi(inY); mfi.isValid(); ++mfi)
    {
        FArrayBox& outHc = outH[mfi];
        FArrayBox& outYc = outY[mfi];
        const FArrayBox& Yc = inY[mfi];
        const FArrayBox& Tc = inT[mfi];
        const Box& box = mfi.validbox();
        const Box gbox = BoxLib::grow(box,nGrow);

        // Get cc mole frac in valid+grow so that we can compute grad(X)
        // NOTE: No worries for corner cells, since C2E only fills surroundingNodes(box,dir)
        Xc.resize(gbox,Nspec);
        ckdriver.massFracToMoleFrac(Xc,Yc,gbox,sCompY,sCompXc);

        for (int i=0; i<BL_SPACEDIM; ++i)
        {
            const Box ebox = BoxLib::surroundingNodes(box,i);
            de.resize(ebox,Nspec+1); // force on edge
            TYe.resize(ebox,Nspec+1); // state on edge

            // Get ec values from cc values
            center_to_edge(Tc,TYe,gbox,sCompT,sCompTe,1);
            center_to_edge(Yc,TYe,gbox,sCompY,sCompYe,Nspec);

            // cc X,T in and ec de out
            FORT_DIFFFORCE(box.loVect(),box.hiVect(),
                           de.dataPtr(sCompde),ARLIM(de.loVect()),ARLIM(de.hiVect()),
                           Xc.dataPtr(sCompXc),ARLIM(Xc.loVect()),ARLIM(Xc.hiVect()),
                           Tc.dataPtr(sCompT), ARLIM(Tc.loVect()),ARLIM(Tc.hiVect()),
                           dx, &i);

            // Get h_i(Te) on edge temperature
            Hie.resize(ebox,Nspec);
            ckdriver.getHGivenT(Hie,TYe,ebox,sCompTe,0);

            // Get diffusion fluxes on edges based on edge state/forces
            FArrayBox& Hfl = fluxH[i][mfi];
            FArrayBox& Yfl = fluxY[i][mfi];
            FORT_FLUX(box.loVect(),box.hiVect(),
                      Hfl.dataPtr(dCompFH),ARLIM(Hfl.loVect()),ARLIM(Hfl.hiVect()),
                      Yfl.dataPtr(dCompFY),ARLIM(Yfl.loVect()),ARLIM(Yfl.hiVect()),
                      de.dataPtr(),        ARLIM(de.loVect()), ARLIM(de.hiVect()),
                      TYe.dataPtr(sCompYe),ARLIM(TYe.loVect()),ARLIM(TYe.hiVect()),
                      TYe.dataPtr(sCompTe),ARLIM(TYe.loVect()),ARLIM(TYe.hiVect()),
                      Hie.dataPtr(),       ARLIM(Hie.loVect()),ARLIM(Hie.hiVect()),
                      &i, &do_add_enth_flux);

            // If whichApp==DD_Temp, this is being used for the Temperature equation.
            //   We do not add hi.fluxi to the energy flux, but we do add
            //     -(grad(T).sum(cpi.fluxi))/cpmix to the total result.
            // HACK: need to implement

            // Multiply fluxes times edge areas
            Hfl.mult(area[i][mfi],0,dCompFH,1);
            for (int n=0; n<Nspec; ++n)
                Yfl.mult(area[i][mfi],0,dCompFY+n,1);

            // Now get flux divergences ( really, Div(flux.Area.dx)/Vol -> dx cancels, not needed )
            const int ncH = 1;
            FORT_INCRDIV(box.loVect(),box.hiVect(),
                         outHc.dataPtr(dCompH),ARLIM(outHc.loVect()),ARLIM(outHc.hiVect()),
                         Hfl.dataPtr(dCompFH), ARLIM(Hfl.loVect()),  ARLIM(Hfl.hiVect()),
                         &i, &ncH);

            FORT_INCRDIV(box.loVect(),box.hiVect(),
                         outYc.dataPtr(dCompY),ARLIM(outYc.loVect()),ARLIM(outYc.hiVect()),
                         Yfl.dataPtr(dCompFY), ARLIM(Yfl.loVect()),  ARLIM(Yfl.hiVect()),
                         &i, &Nspec);
        }

        const FArrayBox& v = volume[mfi];
        outHc.divide(v,0,dCompH,1);
        for (int n=0; n<Nspec; ++n)
            outYc.divide(v,0,dCompY+n,1);
    }
    // Flux returned is extensive (i.e. flux.Area)
}

void
DDOp::setRelax(MultiFab&         lambda,
               int               dCompH,
               const MultiFab&   inT,
               int               sCompT,
               const MultiFab&   inY,
               int               sCompY,
               int               level) const
{
    BL_ASSERT(level >= 0);
    if (level > 0)
    {
        coarser->setRelax(lambda,dCompH,inT,sCompT,inY,sCompY,level-1);
        return;
    }

    const int Nspec = ckdriver.numSpecies();
    const int nGrow = 1;

    BL_ASSERT(lambda.nComp() > dCompH);
    BL_ASSERT(&inT != &inY || (sCompT<sCompY || sCompT >= sCompY+Nspec));
    BL_ASSERT(lambda.boxArray() == grids);
    BL_ASSERT(inT.boxArray() == grids);
    BL_ASSERT(inY.boxArray() == grids);
    BL_ASSERT(inT.nGrow() >= nGrow);
    BL_ASSERT(inY.nGrow() >= nGrow);

    // Need grow cells in X,T to get edge values of Y,T for evaluating relax factor
    // Promise to change only the grow cells in T,Y
    setGrowCells(const_cast<MultiFab&>(inT),sCompT,const_cast<MultiFab&>(inY),sCompY);

    // Get relax factor
    FArrayBox TYe;
    const int sCompTe = Nspec;
    const int sCompYe = 0;
    const Real* dx = Tbd.getGeom().CellSize();
    for (MFIter mfi(lambda); mfi.isValid(); ++mfi)
    {
        FArrayBox& lam = lambda[mfi];
        const FArrayBox& Yc = inY[mfi];
        const FArrayBox& Tc = inT[mfi];

        const Box& box = mfi.validbox();
        const Box gbox = BoxLib::grow(box,nGrow);

        for (int i=0; i<BL_SPACEDIM; ++i)
        {
            Box ebox = BoxLib::surroundingNodes(box,i);
            TYe.resize(ebox,Nspec+1);

            // Get edge values from center values of input data
            center_to_edge(Tc,TYe,gbox,0,sCompTe,1);
            center_to_edge(Yc,TYe,gbox,0,sCompYe,Nspec);

            // Compute thermal conductivity using edge-based data
            FORT_THERM(box.loVect(),box.hiVect(),
                       lam.dataPtr(dCompH), ARLIM(lam.loVect()),ARLIM(lam.hiVect()),
                       TYe.dataPtr(sCompYe),ARLIM(TYe.loVect()),ARLIM(TYe.hiVect()),
                       TYe.dataPtr(sCompTe),ARLIM(TYe.loVect()),ARLIM(TYe.hiVect()),
                       &i);
        }
        
        // Compute lambda / cpmix using cell-centered input data
        FORT_CPSCALE(box.loVect(),box.hiVect(),
                     lam.dataPtr(dCompH),ARLIM(lam.loVect()),ARLIM(lam.hiVect()),
                     Yc.dataPtr(sCompY), ARLIM(Yc.loVect()), ARLIM(Yc.hiVect()),
                     Tc.dataPtr(sCompT), ARLIM(Tc.loVect()), ARLIM(Tc.hiVect()));

    }
}

void
DDOp::average (MultiFab&       mfC,
               int             dCompC,
               const MultiFab& mfF,
               int             sCompF,
               int             nComp)
{
    BL_ASSERT(mfC.nComp()>=dCompC+nComp);
    BL_ASSERT(mfF.nComp()>=sCompF+nComp);
    for (MFIter mfi(mfC); mfi.isValid(); ++mfi)
    {
        FArrayBox& C = mfC[mfi];
        const FArrayBox& F = mfF[mfi];

        const Box& cbox = mfi.validbox();
        FORT_DDCCAVG(C.dataPtr(dCompC),ARLIM(C.loVect()), ARLIM(C.hiVect()),
                     F.dataPtr(sCompF),ARLIM(F.loVect()), ARLIM(F.hiVect()),
                     cbox.loVect(), cbox.hiVect(), &nComp, MGIV.getVect());
    }
}

void
DDOp::interpolate (MultiFab&       mfF,
                   int             dCompF,
                   const MultiFab& mfC,
                   int             sCompC,
                   int             nComp)
{
    BL_ASSERT(mfF.nComp()>=dCompF+nComp);
    BL_ASSERT(mfC.nComp()>=sCompC+nComp);
    for (MFIter mfi(mfF); mfi.isValid(); ++mfi)
    {
        FArrayBox& F = mfF[mfi];
        const FArrayBox& C = mfC[mfi];
        const Box cbox = BoxLib::refine(mfi.validbox(),MGIV);
        FORT_DDCCINT(F.dataPtr(dCompF),ARLIM(F.loVect()), ARLIM(F.hiVect()),
                     C.dataPtr(sCompC),ARLIM(C.loVect()), ARLIM(C.hiVect()),
                     cbox.loVect(), cbox.hiVect(), &nComp, MGIV.getVect());
    }
}
