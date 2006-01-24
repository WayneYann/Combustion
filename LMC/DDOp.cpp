#include "DDOp.H"
#include "DDOp_F.H"
#include "LO_F.H"

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
DDOp::define (const BoxArray&      _grids,
              const Box&           box,
              const IntVect&       ratio,
              int                  mgLevel)
{
    grids = _grids;
    Geometry geom(box);
    Tbd.define(grids,1,geom,mgLevel);
    Ybd.define(grids,ckdriver.numSpecies(),geom,mgLevel);
    Xbd.define(grids,ckdriver.numSpecies(),geom,mgLevel);
    cfRatio = ratio;
    const int gGrow = 0;
    geom.GetVolume(volume,grids,gGrow);
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
        geom.GetFaceArea(area[dir],grids,dir,gGrow);

    const int nGrowOp = 1;
    Soln.define(grids,ckdriver.numSpecies()+1,nGrowOp,Fab_allocate);

    // Generate coarser one, if possible
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
                           &nc, &face);
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
                               &nc, &face);
            }
        }
    }
}                 

void
DDOp::setBoundaryData_Mass(const MultiFab&      fineT,
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
    if (coarser)
    {
        const int nspecies = ckdriver.numSpecies();
        const int newTmfComp = nspecies;
        const int newTbrComp = nspecies;
        const int newYmfComp = 0;
        const int newYbrComp = 0;
        const BoxArray cBA = BoxArray(grids).coarsen(MGIV);
        MultiFab newMF(cBA,nspecies+1,1);
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
                newBR->define(fi(),IndexType::TheCellType(),0,1,1,nspecies+1);
#ifndef NDEBUG
            newBR->setVal(-1.0);
#endif
        }
        coarsenBndryData(fineT,fStartT,cbrT,cStartT,
                         newMF,newTmfComp,newBR,newTbrComp,1);
        coarsenBndryData(fineY,fStartY,cbrY,cStartY,
                         newMF,newYmfComp,newBR,newYbrComp,nspecies);
        coarser->setBoundaryData_Mass(newMF,newTmfComp,newMF,newYmfComp,
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

    const int Nspec = ckdriver.numSpecies();
    const int nGrow = 1;
    MultiFab fineX(grids,Nspec,nGrow);
    const int fStartX = 0;
    for (MFIter mfi(fineX); mfi.isValid(); ++mfi)
    {
        const Box gbox = BoxLib::grow(mfi.validbox(),nGrow);
        ckdriver.massFracToMoleFrac(fineX[mfi],fineY[mfi],gbox,fStartY,fStartX);
    }

    if (cbrY == 0)
    {
        Ybd.setBndryValues(fineY,fStartY,0,Nspec,bcY);
        Xbd.setBndryValues(fineX,fStartX,0,Nspec,bcY);
    }
    else
    {
        const BoxArray cbrGrids = BoxArray(grids).coarsen(cfRatio);
        BndryRegister cbrX(cbrGrids,0,1,1,Nspec);
        const int cStartX = 0;
        for (OrientationIter oitr; oitr; ++oitr)
        {
            FabSet& fsX = cbrX[oitr()];
            const FabSet& fsY = (*cbrY)[oitr()];
            for (FabSetIter fsi(fsX); fsi.isValid(); ++fsi)
            {
                ckdriver.massFracToMoleFrac(fsX[fsi],fsY[fsi],fsX[fsi].box(),cStartY,cStartX);
            }
        }
        Ybd.setBndryValues(const_cast<BndryRegister&>(*cbrY),
                           cStartY,fineY,fStartY,0,Nspec,ratio,bcY);
        Xbd.setBndryValues(cbrX,cStartX,fineX,fStartX,0,Nspec,ratio,bcY);
    }
}

void
DDOp::setBoundaryData_Mole(const MultiFab&      fineT,
                           int                  fStartT,
                           const MultiFab&      fineX,
                           int                  fStartX,
                           const BndryRegister* cbrT,
                           int                  cStartT,
                           const BndryRegister* cbrX,
                           int                  cStartX,
                           const BCRec&         bcT,
                           const BCRec&         bcX)
{
    BL_ASSERT(fineT.boxArray() == grids);
    BL_ASSERT(fineX.boxArray() == grids);
    if (coarser)
    {
        const int nspecies = ckdriver.numSpecies();
        const int newTmfComp = nspecies;
        const int newTbrComp = nspecies;
        const int newXmfComp = 0;
        const int newXbrComp = 0;
        const BoxArray cBA = BoxArray(grids).coarsen(MGIV);
        MultiFab newMF(cBA,nspecies+1,1);
        BoxArray ccBA = BoxArray(cBA).coarsen(MGIV);
        BndryRegister* newBR = 0;
        if (cbrT && cbrX && coarser->coarser)
        {
            newBR = new BndryRegister();
            newBR->setBoxes(ccBA);
            for (OrientationIter fi; fi; ++fi)
                newBR->define(fi(),IndexType::TheCellType(),0,1,1,nspecies+1);
        }
        coarsenBndryData(fineT,fStartT,cbrT,cStartT,
                         newMF,newTmfComp,newBR,newTbrComp,1);
        coarsenBndryData(fineX,fStartX,cbrX,cStartX,
                         newMF,newXmfComp,newBR,newXbrComp,nspecies);
        coarser->setBoundaryData_Mole(newMF,newTmfComp,newMF,newXmfComp,
                                      newBR,newTbrComp,newBR,newXbrComp,bcT,bcX);
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

    const int Nspec = ckdriver.numSpecies();
    const int nGrow = 1;
    MultiFab fineY(grids,Nspec,nGrow);
    const int fStartY = 0;
    for (MFIter mfi(fineY); mfi.isValid(); ++mfi)
    {
        const Box gbox = BoxLib::grow(mfi.validbox(),nGrow);
        FArrayBox& Y = fineY[mfi];
        const FArrayBox& X = fineX[mfi];
        ckdriver.moleFracToMassFrac(Y,X,gbox,fStartX,fStartY);
    }

    if (cbrX == 0)
    {
        Xbd.setBndryValues(fineX,fStartX,0,Nspec,bcX);
        Ybd.setBndryValues(fineY,fStartY,0,Nspec,bcX);
    }
    else
    {
        const BoxArray cbrGrids = BoxArray(grids).coarsen(cfRatio);
        BndryRegister cbrY(cbrGrids,0,1,1,Nspec);
        const int cStartY = 0;
        for (OrientationIter oitr; oitr; ++oitr)
        {
            FabSet& fsY = cbrY[oitr()];
            const FabSet& fsX = (*cbrX)[oitr()];
            for (FabSetIter fsi(fsY); fsi.isValid(); ++fsi)
            {
                FArrayBox& Y = fsY[fsi];
                const FArrayBox& X = fsY[fsi];
                ckdriver.moleFracToMassFrac(Y,X,X.box(),cStartX,cStartY);
            }
        }
        Xbd.setBndryValues(const_cast<BndryRegister&>(*cbrX),
                           cStartX,fineX,fStartX,0,Nspec,ratio,bcX);
        Ybd.setBndryValues(cbrY,cStartY,fineY,fStartY,0,Nspec,ratio,bcX);
    }
}

void
DDOp::setGrowCells(MultiFab& T,
                   int       compT,
                   MultiFab& X,
                   int       compX) const
{
    BL_ASSERT(T.nGrow() >= 1);
    BL_ASSERT(X.nGrow() >= 1);
    const int Nspec = ckdriver.numSpecies();

    BL_ASSERT(T.nComp() > compT);
    BL_ASSERT(X.nComp() >= compX + Nspec);
    BL_ASSERT(&T != &X || (compT<compX || compT >= compX+Nspec));

    BL_ASSERT(T.boxArray() == grids);
    BL_ASSERT(X.boxArray() == grids);

    T.FillBoundary(compT,1);
    X.FillBoundary(compX,Nspec);

    const int flagbc  = 1;
    const int flagden = 0; // Use LinOp's bc interpolator, but don't save the coeff
    const int maxorder = 3;
    Real* dummy;
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
        
        const Array<Array<BoundCond> >&  Xbc = Xbd.bndryConds(face);
        const Array<Real>&      Xloc = Xbd.bndryLocs(face);
        const FabSet&            Xfs = Xbd.bndryValues(face);
        const int                Xnc = Nspec;

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

            FArrayBox& Xfab = X[mfi];
            const FArrayBox& Xb = Xfs[mfi];

            const Mask& Xm  = Xbd.bndryMasks(face)[idx];
            const Real Xbcl = Xloc[idx];
            const int Xbct  = Xbc[idx][comp];

            FORT_APPLYBC(&flagden, &flagbc, &maxorder,
                         Xfab.dataPtr(compX), ARLIM(Xfab.loVect()), ARLIM(Xfab.hiVect()),
                         &iFace, &Xbct, &Xbcl,
                         Xb.dataPtr(), ARLIM(Xb.loVect()), ARLIM(Xb.hiVect()),
                         Xm.dataPtr(), ARLIM(Xm.loVect()), ARLIM(Xm.hiVect()),
                         dummy, ARLIM(dumbox.loVect()), ARLIM(dumbox.hiVect()),
                         vbox.loVect(),vbox.hiVect(), &Xnc, dx);
        }
    }
}

void
DDOp::cellToEdge(MultiFab&       Te,
                 int             dCompTe,
                 MultiFab&       Ye,
                 int             dCompYe,
                 const MultiFab& Tc,
                 int             sCompTc,
                 const MultiFab& Yc,
                 int             sCompYc,
                 int             nGrow,
                 int             dir) const
{
    // Assumes nGrow grow cells properly filled for Tc and Yc
    BL_ASSERT(Tc.nGrow() >= 1);
    BL_ASSERT(Yc.nGrow() >= 1);
    const int Nspec = ckdriver.numSpecies();
    BL_ASSERT(Tc.boxArray() == grids);
    BL_ASSERT(Yc.boxArray() == grids);

    BL_ASSERT(Tc.nComp() > sCompTc);
    BL_ASSERT(Yc.nComp() >= sCompYc + Nspec);
    BL_ASSERT(&Tc != &Yc || (sCompTc<sCompYc || sCompTc >= sCompYc+Nspec));
    BL_ASSERT(Te.nComp() > dCompTe);
    BL_ASSERT(Ye.nComp() >= dCompYe + Nspec);
    BL_ASSERT(&Te != &Ye || (dCompTe<dCompYe || dCompTe >= dCompYe+Nspec));

    for (MFIter mfi(Tc); mfi.isValid(); ++mfi)
    {
        const Box gbox = BoxLib::grow(mfi.validbox(),nGrow);
        center_to_edge(T[mfi],Te[mfi],gbox,0,dCompTe,1);
        center_to_edge(Y[mfi],Ye[mfi],gbox,0,dCompYe,Nspec);
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
              int               level) const
{
    BL_ASSERT(level >= 0);
    if (level > 0)
    {
        coarser->applyOp(outH,dCompH,outY,dCompY,inT,sCompT,inY,sCompY,
                         fluxH,dCompFH,fluxY,dCompFY,level-1);
        return;
    }

    const int Nspec = ckdriver.numSpecies();

    BL_ASSERT(outH.nComp() > dCompH);
    BL_ASSERT(outY.nComp() >= dCompY + Nspec);
    BL_ASSERT(&inT != &inY || (sCompT<sCompY || sCompT >= sCompY+Nspec));
    BL_ASSERT(&outH != &outY || (dCompH<dCompY || dCompH >= dCompY+Nspec));
    BL_ASSERT(outH.boxArray() == grids);
    BL_ASSERT(outY.boxArray() == grids);
    BL_ASSERT(inT.boxArray() == grids);
    BL_ASSERT(inY.boxArray() == grids);

    // Need grow cells in X,T to compute forcing, and edge values of Y,T for evaluating coeffs
    // Promise to change only the grow cells in T,Y
    setGrowCells(const_cast<MultiFab&>(inT),sCompT,const_cast<MultiFab&>(inY),sCompY);

    // Get CC mole fractions over valid+nGrow
    const int nGrow = 1;
    MultiFab inX(grids,Nspec,nGrow);
    const int sCompX = 0;
    for (MFIter mfi(inX); mfi.isValid(); ++mfi)
    {
        const Box gbox = BoxLib::grow(mfi.validbox(),nGrow);
        ckdriver.massFracToMoleFrac(inX[mfi],inY[mfi],gbox,sCompY,sCompX);
    }

    // Compute edge-based T and Y for evaluating diffusion fluxes
    const int sCompTe = Nspec;
    const int sCompYe = 0;
    PArray<MultiFab> TYe(BL_SPACEDIM,PArrayManage);
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        TYe.set(dir,new MultiFab(BoxArray(grids).surroundingNodes(dir),
                                 Nspec+1,0,Fab_allocate));
        cellToEdge(TYe[dir],sCompTe,TYe[dir],sCompYe,inT,sCompT,inY,sCompY,dir);
    }

    // Get fluxes
    const Real* dx = Tbd.getGeom().CellSize();
    FArrayBox de, He, X, Y;
    const int sCompd = 0;
    for (MFIter mfi(inY); mfi.isValid(); ++mfi)
    {
        const FArrayBox& Y = inY[mfi];
        const FArrayBox& X = inX[mfi];
        const FArrayBox& T = inT[mfi];

        for (int i=0; i<BL_SPACEDIM; ++i)
        {
            const int idx = mfi.index();
            const Box& box = grids[idx];
            const Box gbox = BoxLib::grow(box,nGrow);
            const Box ebox = BoxLib::surroundingNodes(box,i);
            de.resize(ebox,Nspec+1);

            // CC X,T in and EC de out
            FORT_DIFFFORCE(box.loVect(),box.hiVect(),
                           de.dataPtr(sCompd),ARLIM(de.loVect()),ARLIM(de.hiVect()),
                           X.dataPtr(sCompX), ARLIM(X.loVect()), ARLIM(X.hiVect()),
                           T.dataPtr(sCompT), ARLIM(T.loVect()), ARLIM(T.hiVect()),
                           dx, &i);

            // Get h_i(Te) on edge temperature
            He.resize(ebox,Nspec);
            ckdriver.getHGivenT(He,TYe[i][mfi],ebox,sCompTe,0);

            // Get diffusion fluxes on edges based on edge state/forces
            FORT_FLUX(box.loVect(),box.hiVect(),
                      fluxH[i][mfi].dataPtr(dCompFH),ARLIM(fluxH[i][mfi].loVect()),ARLIM(fluxH[i][mfi].hiVect()),
                      fluxY[i][mfi].dataPtr(dCompFY),ARLIM(fluxY[i][mfi].loVect()),ARLIM(fluxY[i][mfi].hiVect()),
                      de.dataPtr(),ARLIM(de.loVect()),ARLIM(de.hiVect()),
                      TYe[i][mfi].dataPtr(sCompYe), ARLIM(TYe[i][mfi].loVect()),ARLIM(TYe[i][mfi].hiVect()),
                      TYe[i][mfi].dataPtr(sCompTe), ARLIM(TYe[i][mfi].loVect()),ARLIM(TYe[i][mfi].hiVect()),
                      He.dataPtr(),ARLIM(He.loVect()),ARLIM(He.hiVect()),
                      &i);

            // Multiply fluxes times edge areas
            fluxH[i][mfi].mult(area[i][mfi],0,dCompFH,1);
            for (int n=0; n<Nspec; ++n)
                fluxY[i][mfi].mult(area[i][mfi],0,dCompFY+n,1);
        }
    }

    // Compute flux divergence ( really, Div(flux.Area.dx)/Vol )
    outH.setVal(0.0,dCompH,1);
    outY.setVal(0.0,dCompY,Nspec);
    Array<Real> h1(BL_SPACEDIM,1.0);
    for (MFIter mfi(outH); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        FArrayBox Hfab = outH[mfi];
        FArrayBox Yfab = outY[mfi];

        for (int i=0; i<BL_SPACEDIM; ++i)
        {
            const FArrayBox& Hfl = fluxH[i][mfi];
            const FArrayBox& Yfl = fluxY[i][mfi];

            const int ncH = 1;
            FORT_INCRDIV(box.loVect(),box.hiVect(),
                         Hfab.dataPtr(dCompH), ARLIM(Hfab.loVect()), ARLIM(Hfab.hiVect()),
                         Hfl.dataPtr(dCompFH), ARLIM(Hfl.loVect()),  ARLIM(Hfl.hiVect()),
                         h1.dataPtr(), &i, &ncH);

            FORT_INCRDIV(box.loVect(),box.hiVect(),
                         Yfab.dataPtr(dCompY), ARLIM(Yfab.loVect()), ARLIM(Yfab.hiVect()),
                         Yfl.dataPtr(dCompFY), ARLIM(Yfl.loVect()),  ARLIM(Yfl.hiVect()),
                         h1.dataPtr(), &i, &Nspec);
        }

        const FArrayBox& v = volume[mfi];
        Hfab.divide(v,0,dCompH,1);
        for (int n=0; n<Nspec; ++n)
            Yfab.divide(v,0,dCompY+n,1);
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

    BL_ASSERT(lambda.nComp() > dCompH);
    BL_ASSERT(&inT != &inY || (sCompT<sCompY || sCompT >= sCompY+Nspec));
    BL_ASSERT(lambda.boxArray() == grids);
    BL_ASSERT(inT.boxArray() == grids);
    BL_ASSERT(inY.boxArray() == grids);

    // Need grow cells in X,T to get edge values of Y,T for evaluating relax factor
    // Promise to change only the grow cells in T,Y
    setGrowCells(const_cast<MultiFab&>(inT),sCompT,const_cast<MultiFab&>(inY),sCompY);

    // Get relax factor
    FARrayBox TYe;
    const int sCompTe = Nspec;
    const int sCompYe = 0;
    const Real* dx = Tbd.getGeom().CellSize();
    for (MFIter mfi(lambda); mfi.isValid(); ++mfi)
    {
        FArrayBox& lam = lambda[mfi];
        const FArrayBox Yc = inY[mfi];
        const FArrayBox Tc = inT[mfi];

        const Box& box = mfi.validbox();
        const Box gbox = BoxLib::grow(box,nGrow);

        for (int i=0; i<BL_SPACEDIM; ++i)
        {
            Box ebox = BoxLib::surroundingNodes(box,i);
            TYe.resize(ebox,Nspec+1);

            // Get edge values from center values of input data
            center_to_edge(Tc,TYe,gbox,0,dCompTe,1);
            center_to_edge(Yc,TYe,gbox,0,dCompYe,Nspec);

            // Compute thermal conductivity using edge-based data
            FORT_THERM(box.loVect(),box.hiVect(),
                       lam.dataPtr(dCompH),ARLIM(lam.loVect()),ARLIM(lam.hiVect()),
                       TYe.dataPtr(sCompYe),ARLIM(Ye.loVect()), ARLIM(Ye.hiVect()),
                       TYe.dataPtr(sCompTe),ARLIM(Te.loVect()), ARLIM(Te.hiVect()),
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
                     cbox.loVect(), cbox.hiVect(), &nComp);
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
                     cbox.loVect(), cbox.hiVect(), &nComp);
    }
}
