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
    for (MultiFabIterator Ymfi(fineY); Ymfi.isValid(); ++Ymfi)
    {
        const Box gbox = BoxLib::grow(Ymfi.validbox(),nGrow);
        DependentMultiFabIterator Xmfi(Ymfi, fineX);
        ckdriver.moleFracToMassFrac(Ymfi(),Xmfi(),gbox,fStartX,fStartY);
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
            for (FabSetIterator Yfsi(cbrY[oitr()]); Yfsi.isValid(); ++Yfsi)
            {
                DependentFabSetIterator Xfsi(Yfsi, (*cbrX)[oitr()]);
                ckdriver.moleFracToMassFrac(Yfsi(),Xfsi(),Xfsi().box(),cStartX,cStartY);
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
    const Real xInt = -0.5; // Grow cell location wrt face, in dx units, +ve inward
    Real* dummy;
    Box dumbox(IntVect(D_DECL(0,0,0)),IntVect(D_DECL(0,0,0)));
    const Real* dx = Tbd.getGeom().CellSize();
    for (OrientationIter oitr; oitr; ++oitr)
    {
        const Orientation&      face = oitr();
        const int              iFace = (int)face;

        const Array<BoundCond>&  Tbc = Tbd.bndryConds(face);
        const Array<Real>&      Tloc = Tbd.bndryLocs(face);
        const FabSet&            Tfs = Tbd.bndryValues(face);
        const int                Tnc = 1;
        
        const Array<BoundCond>&  Xbc = Xbd.bndryConds(face);
        const Array<Real>&      Xloc = Xbd.bndryLocs(face);
        const FabSet&            Xfs = Xbd.bndryValues(face);
        const int                Xnc = Nspec;

        for (MultiFabIterator Tmfi(T); Tmfi.isValid(); ++Tmfi)
        {
            const int   idx = Tmfi.index();
            const Box& vbox = Tmfi.validbox();
            BL_ASSERT(grids[idx] == vbox);

            DependentFabSetIterator Tfsfsi(Tmfi, Tfs);
            const Mask& Tm  = Tbd.bndryMasks(face)[idx];
            const Real Tbcl = Tloc[idx];
            const int Tbct  = Tbc[idx];

            FORT_APPLYBC(&flagden, &flagbc, &maxorder,
                         Tmfi().dataPtr(compT),
                         ARLIM(Tmfi().loVect()), ARLIM(Tmfi().hiVect()),
                         Tmfi().dataPtr(compT),
                         ARLIM(Tmfi().loVect()), ARLIM(Tmfi().hiVect()),
                         &iFace, &Tbct, &Tbcl,
                         Tfsfsi().dataPtr(),
                         ARLIM(Tfsfsi().loVect()), ARLIM(Tfsfsi().hiVect()),
                         Tm.dataPtr(), ARLIM(Tm.loVect()), ARLIM(Tm.hiVect()),
                         dummy, ARLIM(dumbox.loVect()), ARLIM(dumbox.hiVect()),
                         vbox.loVect(),vbox.hiVect(), &Tnc, dx, &xInt);

            DependentMultiFabIterator Xmfi(Tmfi,X);
            DependentFabSetIterator Xfsfsi(Tmfi, Xfs);
            const Mask& Xm  = Xbd.bndryMasks(face)[idx];
            const Real Xbcl = Xloc[idx];
            const int Xbct  = Xbc[idx];

            FORT_APPLYBC(&flagden, &flagbc, &maxorder,
                         Xmfi().dataPtr(compX),
                         ARLIM(Xmfi().loVect()), ARLIM(Xmfi().hiVect()),
                         Xmfi().dataPtr(compX),
                         ARLIM(Xmfi().loVect()), ARLIM(Xmfi().hiVect()),
                         &iFace, &Xbct, &Xbcl,
                         Xfsfsi().dataPtr(),
                         ARLIM(Xfsfsi().loVect()), ARLIM(Xfsfsi().hiVect()),
                         Xm.dataPtr(), ARLIM(Xm.loVect()), ARLIM(Xm.hiVect()),
                         dummy, ARLIM(dumbox.loVect()), ARLIM(dumbox.hiVect()),
                         vbox.loVect(),vbox.hiVect(), &Xnc, dx, &xInt);
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
                 int             dir) const
{
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

    const int nGrow = 1;
    MultiFab T(grids,1,nGrow,Fab_allocate);
    MultiFab Y(grids,Nspec,nGrow,Fab_allocate);
    for (MultiFabIterator Tmfi(T); Tmfi.isValid(); ++Tmfi)
    {
        DependentMultiFabIterator Ymfi(Tmfi, Y);
        DependentMultiFabIterator Tcmfi(Tmfi, Tc);
        DependentMultiFabIterator Ycmfi(Tmfi, Yc);
        const Box& box = Tmfi.validbox();
        Tmfi().copy(Tcmfi(),box,sCompTc,box,0,1);
        Ymfi().copy(Ycmfi(),box,sCompYc,box,0,Nspec);
    }
    
    const Real bogusVal = -1.e20;
    T.setBndry(bogusVal,0,1);
    Y.setBndry(bogusVal,0,Nspec);
    T.FillBoundary(0,1);
    Y.FillBoundary(0,Nspec);

    // Simple averaging will get interior edges, and fine-fine edges correct
    // The stuff below will fix up c-f and phys bc's
    for (MultiFabIterator Tmfi(T); Tmfi.isValid(); ++Tmfi)
    {
        DependentMultiFabIterator Ymfi(Tmfi, Y);
        DependentMultiFabIterator Temfi(Tmfi, Te);
        DependentMultiFabIterator Yemfi(Tmfi, Ye);
        const Box& gbox = BoxLib::grow(Tmfi.validbox(),nGrow);
        
        center_to_edge(Tmfi(),Temfi(),gbox,0,dCompTe,1);
        center_to_edge(Ymfi(),Yemfi(),gbox,0,dCompYe,Nspec);
    }

    // Now, use LinOp's interp stuff to get edge values on c-f and phys bc
    const int flagbc  = 1;
    const int flagden = 0; // Use LinOp's bc interpolator, but don't save the coeff
    const int maxorder = 3;
    const Real xInt = 0.0; // Grow cell location wrt face, in dx units, +ve inward
    Real* dummy;
    Box dumbox(IntVect(D_DECL(0,0,0)),IntVect(D_DECL(0,0,0)));
    const Real* dx = Tbd.getGeom().CellSize();
    Orientation face(dir,Orientation::low);
    for (int side=0; side<2; ++side, face = face.flip())
    {
        const int              iFace = (int)face;

        const Array<BoundCond>&  Tbc = Tbd.bndryConds(face);
        const Array<Real>&      Tloc = Tbd.bndryLocs(face);
        const FabSet&            Tfs = Tbd.bndryValues(face);
        const int                Tnc = 1;
        
        const Array<BoundCond>&  Ybc = Ybd.bndryConds(face);
        const Array<Real>&      Yloc = Ybd.bndryLocs(face);
        const FabSet&            Yfs = Ybd.bndryValues(face);
        const int                Ync = Nspec;

        for (MultiFabIterator Tmfi(T); Tmfi.isValid(); ++Tmfi)
        {
            const int idx = Tmfi.index();
            const Box& vbox = Tmfi.validbox();
            BL_ASSERT(grids[idx] == vbox);
            const Box gbox = BoxLib::grow(vbox,nGrow);

            DependentFabSetIterator Tfsfsi(Tmfi, Tfs);
            DependentMultiFabIterator Temfi(Tmfi, Te);
            const Mask& Tm  = Tbd.bndryMasks(face)[idx];
            const Real Tbcl = Tloc[idx];
            const int Tbct  = Tbc[idx];

            // Shift edge data to look like cell data for interpolator
            IntVect iv(D_DECL(0,0,0));
            if (face.isLow())
                iv = -BoxLib::BASISV(dir);

            Temfi().shift(iv);
            FORT_APPLYBC(&flagden, &flagbc, &maxorder,
                         Tmfi().dataPtr(),ARLIM(Tmfi().loVect()), ARLIM(Tmfi().hiVect()),
                         Temfi().dataPtr(dCompTe),
                         ARLIM(Temfi().loVect()), ARLIM(Temfi().hiVect()),
                         &iFace, &Tbct, &Tbcl,
                         Tfsfsi().dataPtr(),
                         ARLIM(Tfsfsi().loVect()), ARLIM(Tfsfsi().hiVect()),
                         Tm.dataPtr(), ARLIM(Tm.loVect()), ARLIM(Tm.hiVect()),
                         dummy, ARLIM(dumbox.loVect()), ARLIM(dumbox.hiVect()),
                         vbox.loVect(),vbox.hiVect(), &Tnc, dx, &xInt);
            Temfi().shift(-iv);

            DependentMultiFabIterator Ymfi(Tmfi,Y);
            DependentMultiFabIterator Yemfi(Tmfi,Ye);
            DependentFabSetIterator Yfsfsi(Tmfi, Yfs);
            const Mask& Ym  = Ybd.bndryMasks(face)[idx];
            const Real Ybcl = Yloc[idx];
            const int Ybct  = Ybc[idx];
            
            Yemfi().shift(iv);
            FORT_APPLYBC(&flagden, &flagbc, &maxorder,
                         Ymfi().dataPtr(),
                         ARLIM(Ymfi().loVect()), ARLIM(Ymfi().hiVect()),
                         Yemfi().dataPtr(dCompYe),
                         ARLIM(Yemfi().loVect()), ARLIM(Yemfi().hiVect()),
                         &iFace, &Ybct, &Ybcl,
                         Yfsfsi().dataPtr(),
                         ARLIM(Yfsfsi().loVect()), ARLIM(Yfsfsi().hiVect()),
                         Ym.dataPtr(), ARLIM(Ym.loVect()), ARLIM(Ym.hiVect()),
                         dummy, ARLIM(dumbox.loVect()), ARLIM(dumbox.hiVect()),
                         vbox.loVect(),vbox.hiVect(), &Ync, dx, &xInt);
            Yemfi().shift(-iv);
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

    // Get CC mole fractions over valid (for computing forces)
    const int nGrow = 1;
    MultiFab inX(grids,Nspec,nGrow);
    const int sCompX = 0;
    for (MultiFabIterator Xmfi(inX); Xmfi.isValid(); ++Xmfi)
    {
        DependentMultiFabIterator Ymfi(Xmfi, inY);
        ckdriver.massFracToMoleFrac(Xmfi(),Ymfi(),Xmfi.validbox(),sCompY,sCompX);
    }

    // Need edge-based T and Y for evaluating diffusion fluxes below
    const int sCompTe = Nspec;
    const int sCompYe = 0;
    PArray<MultiFab> TYe(BL_SPACEDIM,PArrayManage);
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        TYe.set(dir,new MultiFab(BoxArray(grids).surroundingNodes(dir),
                                 Nspec+1,0,Fab_allocate));
        cellToEdge(TYe[dir],sCompTe,TYe[dir],sCompYe,inT,sCompT,inY,sCompY,dir);
    }

    // Need grow cells in X,T to compute diffusion force below
    // Promise to change only the grow cells in T
    setGrowCells(const_cast<MultiFab&>(inT),sCompT,inX,sCompX);
    
    // Get fluxes
    const Real* dx = Tbd.getGeom().CellSize();
    FArrayBox de, He, X, Y;
    const int sCompd = 0;
    for (ConstMultiFabIterator Ymfi(inY); Ymfi.isValid(); ++Ymfi)
    {
        for (int i=0; i<BL_SPACEDIM; ++i)
        {
            const int idx = Ymfi.index();
            const Box& box = grids[idx];
            const Box gbox = BoxLib::grow(box,nGrow);
            const Box ebox = surroundingNodes(box,i);

            const FArrayBox& Y = Ymfi();
            const FArrayBox& X = inX[idx];
            const FArrayBox& T = inT[idx];
            de.resize(ebox,Nspec+1);

            // CC X,T in and EC de out
            FORT_DIFFFORCE(box.loVect(),box.hiVect(),
                           de.dataPtr(sCompd),ARLIM(de.loVect()),ARLIM(de.hiVect()),
                           X.dataPtr(sCompX), ARLIM(X.loVect()), ARLIM(X.hiVect()),
                           T.dataPtr(sCompT), ARLIM(T.loVect()), ARLIM(T.hiVect()),
                           dx, &i);

            // Get h_i(Te) on edge temperature
            He.resize(ebox,Nspec);
            ckdriver.getHGivenT(He,TYe[i][idx],ebox,sCompTe,0);

            // Get diffusion fluxes on edges based on edge state/forces
            FORT_FLUX(box.loVect(),box.hiVect(),
                      fluxH[i][idx].dataPtr(dCompFH),
                      ARLIM(fluxH[i][idx].loVect()),ARLIM(fluxH[i][idx].hiVect()),
                      fluxY[i][idx].dataPtr(dCompFY),
                      ARLIM(fluxY[i][idx].loVect()),ARLIM(fluxY[i][idx].hiVect()),
                      de.dataPtr(),ARLIM(de.loVect()),ARLIM(de.hiVect()),
                      TYe[i][idx].dataPtr(sCompYe),
                      ARLIM(TYe[i][idx].loVect()),ARLIM(TYe[i][idx].hiVect()),
                      TYe[i][idx].dataPtr(sCompTe),
                      ARLIM(TYe[i][idx].loVect()),ARLIM(TYe[i][idx].hiVect()),
                      He.dataPtr(),ARLIM(He.loVect()),ARLIM(He.hiVect()),
                      &i);

            // Multiply fluxes times edge areas
            ConstDependentMultiFabIterator ai(Ymfi, area[i]);
            fluxH[i][idx].mult(ai(),0,dCompFH,1);
            for (int n=0; n<Nspec; ++n)
                fluxY[i][idx].mult(ai(),0,dCompFY+n,1);
        }
    }

    // Compute flux divergence ( really, Div(flux.Area.dx)/Vol )
    outH.setVal(0.0,dCompH,1);
    outY.setVal(0.0,dCompY,Nspec);
    Array<Real> h1(BL_SPACEDIM,1.0);
    for (MultiFabIterator outHmfi(outH); outHmfi.isValid(); ++outHmfi)
    {
        DependentMultiFabIterator outYmfi(outHmfi, outY);
        const int idx = outHmfi.index();
        const Box& box = outHmfi.validbox();

        for (int i=0; i<BL_SPACEDIM; ++i)
        {
            const int ncH = 1;
            FORT_INCRDIV(box.loVect(),box.hiVect(),
                         outHmfi().dataPtr(dCompH),
                         ARLIM(outHmfi().loVect()),ARLIM(outHmfi().hiVect()),
                         fluxH[i][idx].dataPtr(dCompFH),
                         ARLIM(fluxH[i][idx].loVect()),ARLIM(fluxH[i][idx].hiVect()),
                         h1.dataPtr(), &i, &ncH);

            FORT_INCRDIV(box.loVect(),box.hiVect(),
                         outYmfi().dataPtr(dCompY),
                         ARLIM(outYmfi().loVect()),ARLIM(outYmfi().hiVect()),
                         fluxY[i][idx].dataPtr(dCompFY),
                         ARLIM(fluxY[i][idx].loVect()),ARLIM(fluxY[i][idx].hiVect()),
                         h1.dataPtr(), &i, &Nspec);
        }

        DependentMultiFabIterator vmfi(outHmfi, volume);
        outHmfi().divide(vmfi(),0,dCompH,1);
        for (int n=0; n<Nspec; ++n)
            outYmfi().divide(vmfi(),0,dCompY+n,1);
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

    // Get mole fractions over valid
    const int nGrow = 1;
    MultiFab inX(grids,Nspec,nGrow);
    const int sCompX = 0;
    for (MultiFabIterator Xmfi(inX); Xmfi.isValid(); ++Xmfi)
    {
        DependentMultiFabIterator Ymfi(Xmfi, inY);
        ckdriver.massFracToMoleFrac(Xmfi(),Ymfi(),Xmfi.validbox(),sCompY,sCompX);
    }

    // Need edge-based T and Y for evaluating diffusion fluxes below
    const int sCompTe = Nspec;
    const int sCompYe = 0;
    PArray<MultiFab> TYe(BL_SPACEDIM,PArrayManage);
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        TYe.set(dir,new MultiFab(BoxArray(grids).surroundingNodes(dir),
                                 Nspec+1,0,Fab_allocate));
        cellToEdge(TYe[dir],sCompTe,TYe[dir],sCompYe,inT,sCompT,inY,sCompY,dir);
    }

    // Need grow cells in X,T to compute diffusion force below
    // Promise to change only the grow cells in T
    setGrowCells(const_cast<MultiFab&>(inT),sCompT,inX,sCompX);
    
    // Get fluxes
    const Real* dx = Tbd.getGeom().CellSize();
    for (MultiFabIterator Lmfi(lambda); Lmfi.isValid(); ++Lmfi)
    {
        const int idx = Lmfi.index();
        const Box& box = grids[idx];
        for (int i=0; i<BL_SPACEDIM; ++i)
        {
            // Compute thermal conductivity
            FORT_THERM(box.loVect(),box.hiVect(),
                       Lmfi().dataPtr(dCompH),      ARLIM(Lmfi().loVect()),     ARLIM(Lmfi().hiVect()),
                       TYe[i][idx].dataPtr(sCompYe),ARLIM(TYe[i][idx].loVect()),ARLIM(TYe[i][idx].hiVect()),
                       TYe[i][idx].dataPtr(sCompTe),ARLIM(TYe[i][idx].loVect()),ARLIM(TYe[i][idx].hiVect()),
                       &i);

        }

        DependentMultiFabIterator Ymfi(Lmfi, inY);
        const FArrayBox& Y = Ymfi();
        const FArrayBox& T = inT[idx];
        
        // Compute lambda / cpmix
        FORT_CPSCALE(box.loVect(),box.hiVect(),
                     Lmfi().dataPtr(dCompH),ARLIM(Lmfi().loVect()),ARLIM(Lmfi().hiVect()),
                     Y.dataPtr(sCompY),     ARLIM(Y.loVect()), ARLIM(Y.hiVect()),
                     T.dataPtr(sCompT),     ARLIM(T.loVect()), ARLIM(T.hiVect()));

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
    for (MultiFabIterator C_mfi(mfC); C_mfi.isValid(); ++C_mfi)
    {
        DependentMultiFabIterator F_mfi(C_mfi, mfF);
        const Box& box = C_mfi.validbox();
        BL_ASSERT(BoxLib::refine(box,MGIV) == F_mfi.validbox());
        FORT_DDCCAVG(C_mfi().dataPtr(dCompC),
                     ARLIM(C_mfi().loVect()), ARLIM(C_mfi().hiVect()),
                     F_mfi().dataPtr(sCompF),
                     ARLIM(F_mfi().loVect()), ARLIM(F_mfi().hiVect()),
                     box.loVect(), box.hiVect(), &nComp);
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
    for (MultiFabIterator F_mfi(mfF); F_mfi.isValid(); ++F_mfi)
    {
        DependentMultiFabIterator C_mfi(F_mfi, mfC);
        const Box& box = C_mfi.validbox();
        BL_ASSERT(BoxLib::refine(box,MGIV) == F_mfi.validbox());
        FORT_DDCCINT(F_mfi().dataPtr(dCompF),
                     ARLIM(F_mfi().loVect()), ARLIM(F_mfi().hiVect()),
                     C_mfi().dataPtr(sCompC),
                     ARLIM(C_mfi().loVect()), ARLIM(C_mfi().hiVect()),
                     box.loVect(), box.hiVect(), &nComp);
    }
}
