#include <winstd.H>

#include "DDOp.H"
#include "DDOp_F.H"
#include "LO_F.H"

#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
#include "BoxLib_Data_Dump.H"


// Note: The ratio between the original grids and the coarse-fine boundary data
//       registers can be anything, but the ratio between adjacent DD operators
//       generated for multigrid will always be two.
const IntVect MGIV = IntVect(D_DECL(2,2,2));
DDOp::DD_Model DDOp::transport_model = DDOp::DD_Model_NumModels; // ...must set prior to first ctr call
static std::string DD_FULL_NAME("DD_Model_Full"); // ParmParsed name for full
static std::string DD_MIX_NAME("DD_Model_MixAvg"); // ParmParsed name for mixture-averaged

ChemDriver* DDOp::chem = 0; // ...must set prior to first ctr call
int DDOp::maxorder = 3;
int DDOp::mgLevelsMAX = -1;
int DDOp::mgLevels = 0;
bool DDOp::first_define = true;

DDOp::DDOp ()
{
}

DDOp::DDOp (const BoxArray&   grids,
            const Box&        box,
            const IntVect&    amrRatio)
    : Tbd(0), Ybd(0)
{
    BL_ASSERT(ddOps.size()==0);
    ddOps.resize(1);
    ddOps.set(0,this);
    define(grids,box,amrRatio,0,this);
}

DDOp::DDOp (const BoxArray&   grids,
            const Box&        box,
            const IntVect&    amrRatio,
            int               mgLevel,
            DDOp*             parent)
    : Tbd(0), Ybd(0)
{
    define(grids,box,amrRatio,mgLevel,parent);
}

DDOp::~DDOp ()
{
    if (mg_parent==this) {
        BL_ASSERT(mg_level==0);
        delete Tbd;
        delete Ybd;
        for (int i=1; i<mgLevels; ++i)
            delete ddOps.remove(i);
    }
}

const BoxArray&
DDOp::boxArray(int level) const
{
    BL_ASSERT( (level>=0)  &&  (level<mgLevels) );
    BL_ASSERT(mg_level==0);
    return ddOps[level].grids;
}

void
DDOp::ensure_valid_transport_is_set() const
{
    std::string id;
    if (transport_model==DD_Model_Full)
    {
        id = DD_FULL_NAME;
    }
    else if (transport_model==DD_Model_MixAvg) 
    {
        id = DD_MIX_NAME;
    }
    else
    {
        BoxLib::Abort("Must set the static DDOp::transport_model");
    }
    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "DDOp transport model: " << id << std::endl;
    }
}

bool can_coarsen(const BoxArray& ba)
{
    // Ratio between levels here will always be MGIV
    for (int i = 0; i < ba.size(); ++i)
    {
        Box tmp = ba[i];
        Box ctmp  = BoxLib::coarsen(ba[i],MGIV);
        Box rctmp = BoxLib::refine(ctmp,MGIV);
        if (tmp != rctmp || ctmp.numPts() <= 1)
            return false;
    }
    return true;
}

bool can_coarsen(const Box& box)
{
    // Ratio between levels here will always be MGIV
    Box ctmp  = BoxLib::coarsen(box,MGIV);
    Box rctmp = BoxLib::refine(ctmp,MGIV);
    if (box != rctmp || ctmp.numPts() <= 1)
        return false;
    return true;
}

void
DDOp::define (const BoxArray& _grids,
              const Box&      box,
              const IntVect&  amrRatio)
{
    define(_grids,box,amrRatio,0,this);
}

void
DDOp::define (const BoxArray& _grids,
              const Box&      box,
              const IntVect&  amrRatio,
              int             mgLevel,
              DDOp*           parent)
{
    mg_parent = parent;
    mg_level = mgLevel;

    if (!chem) {
        BoxLib::Abort("ChemDriver must be set prior to first DDOp ctr call");
    }
    if (first_define) {
        ensure_valid_transport_is_set();
        first_define = false;
    }
    amr_ratio = amrRatio;
    const IntVect mg_ratio = MGIV;
    grids = _grids;
    geom.define(box);

    int Nspec = chem->numSpecies();

    if (mgLevel==0) 
    {
        mgLevels = 1;
        BL_ASSERT(ddOps.size()==0);
        ddOps.resize(1,PArrayNoManage);
        ddOps.set(mgLevels-1,this);

        Box cbox(box);
        BoxArray fgrids, cgrids;
        bool not_done = can_coarsen(grids)  &&  ( (mgLevelsMAX<0) || (mgLevels < mgLevelsMAX) );
        while (not_done)
        {
            fgrids = ddOps[mgLevels-1].grids;
            cgrids = BoxArray(fgrids).coarsen(mg_ratio);
            cbox.coarsen(mg_ratio);

            mgLevels++;
            ddOps.resize(mgLevels);
            ddOps.set(mgLevels-1,new DDOp(cgrids,cbox,amr_ratio,mgLevels-1,parent));

            not_done = can_coarsen(cgrids)  &&  can_coarsen(cbox) 
                &&  ( (mgLevelsMAX<0) || (mgLevels < mgLevelsMAX) );
        }

        Tbd = new DDBndry(grids,1,geom,mgLevels);
        Ybd = new DDBndry(grids,Nspec,geom,mgLevels);
    }

    dx.resize(BL_SPACEDIM);
    for (int i=0; i<BL_SPACEDIM; ++i)
        dx[i] = geom.CellSize()[i];
    
    const int gGrow = 0;
    geom.GetVolume(volInv,grids,gGrow);
    for (MFIter mfi(volInv); mfi.isValid(); ++mfi) {
        volInv[mfi].invert(1);
    }
    area.resize(BL_SPACEDIM,PArrayManage);
    for (int dir = 0; dir < BL_SPACEDIM; dir++) {
        area.set(dir,new MultiFab());
        geom.GetFaceArea(area[dir],grids,dir,gGrow);
    }
    /*  Not yet...
    flux.resize(BL_SPACEDIM,PArrayManage);
    for (int dir = 0; dir < BL_SPACEDIM; dir++) {
        BoxArray eba=BoxArray(grids).surroundingNodes(dir);
        flux.set(dir,new MultiFab(eba,Nspec+1,0));
    }
    */
    int model_DD0_MA1 = (transport_model==DD_Model_Full ? 0 : 1);
    int nComp = FORT_DDNCOEFS(model_DD0_MA1);
    int nGrow = 1;
    coefs.define(grids,nComp,nGrow,Fab_allocate);

    // Weights used to generate grow cell data
    stencilWeight.setBoxes(grids);
    for (OrientationIter face; face; ++face)
    {
        int in_rad_st = 1;
        int out_rad_st = 0;
        int extent_rad_st = 0;
        int ncomp_st = 2;
        stencilWeight.define(face(),IndexType::TheCellType(),in_rad_st,out_rad_st,extent_rad_st,ncomp_st);
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
               cfab.dataPtr(sComp),ARLIM(cfab.loVect()),ARLIM(cfab.hiVect()),
               efab.dataPtr(dComp),ARLIM(efab.loVect()),ARLIM(efab.hiVect()),
               &nComp, &dir);
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
                      const BCRec&         bcY,
                      int                  max_order)
{
    BL_ASSERT(mg_level==0);
    BL_ASSERT(fineY.boxArray() == grids);
    const int Nspec = chem->numSpecies();
    BL_ASSERT( (cbrT==0 ? cbrY==0 : cbrY!=0 ) );

    if (cbrT==0) {

        Tbd->setBndryValues(fineT,fStartT,0,1,bcT);
        Ybd->setBndryValues(fineY,fStartY,0,Nspec,bcY);

    } else {

        Tbd->setBndryValues(const_cast<BndryRegister&>(*cbrT),
                           cStartT,fineT,fStartT,0,1,amr_ratio,bcT,max_order);
        Ybd->setBndryValues(const_cast<BndryRegister&>(*cbrY),
                           cStartY,fineY,fStartY,0,Nspec,amr_ratio,bcY,max_order);
    }
}

void
DDOp::setGrowCells(MultiFab& YT,
                   int       compT,
                   int       compY)
{
    BL_ASSERT(YT.nGrow() >= 1);
    const int Nspec = chem->numSpecies();

    BL_ASSERT(YT.nComp() > compT);
    BL_ASSERT(YT.nComp() >= compY + Nspec);

    BL_ASSERT(YT.boxArray() == grids);

    const int flagbc  = 1;
    const int flagden = 1;

    DDBndry& Tbndry = (*mg_parent->Tbd);
    DDBndry& Ybndry = (*mg_parent->Ybd);

    for (OrientationIter oitr; oitr; ++oitr)
    {
        const Orientation&      face = oitr();
        const int              iFace = (int)face;

        const Array<Array<BoundCond> >& Tbc = Tbndry.bndryConds(face);
        const Array<Real>&      Tloc = Tbndry.bndryLocs(face);
        const FabSet&            Tfs = Tbndry.bndryValues(face,mg_level);
        FabSet&                  Tst = stencilWeight[face];
        const int                Tnc = 1;
        
        const Array<Array<BoundCond> >&  Ybc = Ybndry.bndryConds(face);
        const Array<Real>&      Yloc = Ybndry.bndryLocs(face);
        const FabSet&            Yfs = Ybndry.bndryValues(face,mg_level);
        FabSet&                  Yst = stencilWeight[face];
        const int                Ync = Nspec;

        const int comp = 0;
        for (MFIter mfi(YT); mfi.isValid(); ++mfi)
        {
            const int   idx = mfi.index();
            const Box& vbox = mfi.validbox();
            BL_ASSERT(grids[idx] == vbox);

            FArrayBox& Tfab = YT[mfi];
            FArrayBox& Tstfab = Tst[mfi];
            const FArrayBox& Tb = Tfs[mfi];

            const Mask& Tm  = Tbndry.bndryMasks(face,mg_level)[idx];
            const Real Tbcl = Tloc[idx];
            const int Tbct  = Tbc[idx][comp];

            FORT_APPLYBC(&flagden, &flagbc, &maxorder,
                         Tfab.dataPtr(compT), ARLIM(Tfab.loVect()), ARLIM(Tfab.hiVect()),
                         &iFace, &Tbct, &Tbcl,
                         Tb.dataPtr(), ARLIM(Tb.loVect()), ARLIM(Tb.hiVect()),
                         Tm.dataPtr(), ARLIM(Tm.loVect()), ARLIM(Tm.hiVect()),
                         Tstfab.dataPtr(1), ARLIM(Tstfab.loVect()), ARLIM(Tstfab.hiVect()),
                         vbox.loVect(),vbox.hiVect(), &Tnc, dx.dataPtr());

            FArrayBox& Yfab = YT[mfi];
            FArrayBox& Ystfab = Yst[mfi];
            const FArrayBox& Yb = Yfs[mfi];

            const Mask& Ym  = Ybndry.bndryMasks(face,mg_level)[idx];
            const Real Ybcl = Yloc[idx];
            const int Ybct  = Ybc[idx][comp];

            FORT_APPLYBC(&flagden, &flagbc, &maxorder,
                         Yfab.dataPtr(compY), ARLIM(Yfab.loVect()), ARLIM(Yfab.hiVect()),
                         &iFace, &Ybct, &Ybcl,
                         Yb.dataPtr(), ARLIM(Yb.loVect()), ARLIM(Yb.hiVect()),
                         Ym.dataPtr(), ARLIM(Ym.loVect()), ARLIM(Ym.hiVect()),
                         Ystfab.dataPtr(0), ARLIM(Ystfab.loVect()), ARLIM(Ystfab.hiVect()),
                         vbox.loVect(),vbox.hiVect(), &Ync, dx.dataPtr());
        }
    }


    YT.FillBoundary(compT,1);
    geom.FillPeriodicBoundary(YT,compT,1,true);
    YT.FillBoundary(compY,Nspec);
    geom.FillPeriodicBoundary(YT,compY,Nspec,true);
}

void
DDOp::setCoefficients(const MFIter&    mfi,
                      const FArrayBox& inYT,
                      const FArrayBox& cpi)
{
    int nGrow = 1;
    int Nspec = chem->numSpecies();
    BL_ASSERT(inYT.nComp()>=Nspec+1);

    const Box gbox = Box(grids[mfi.index()]).grow(nGrow);
    BL_ASSERT(inYT.box().contains(gbox));
    BL_ASSERT(cpi.box().contains(gbox));

    int Full0_Mix1 = (transport_model == DD_Model_Full ? 0 : 1);
    FORT_DDCOEFS(gbox.loVect(), gbox.hiVect(),
                 coefs[mfi].dataPtr(), ARLIM(coefs[mfi].loVect()), ARLIM(coefs[mfi].hiVect()),
                 inYT.dataPtr(), ARLIM(inYT.loVect()),  ARLIM(inYT.hiVect()),
                 cpi.dataPtr(), ARLIM(cpi.loVect()),  ARLIM(cpi.hiVect()),
                 Full0_Mix1);
}

bool
DDOp::thisIsParent() const
{
    return ( (mg_level==0)  &&  (this == &(ddOps[0]))  );
}

void
DDOp::applyOp(MultiFab&         outYH,
              const MultiFab&   inYT,
              PArray<MultiFab>& fluxYH,
              DD_ApForTorRH     whichApp,
              bool              updateCoefs,
              int               level,
              bool              getAlpha,
              MultiFab*         alpha)
{
    BL_ASSERT(level >= 0  &&  ((mgLevelsMAX<0) || (level<mgLevelsMAX)) );
    BL_ASSERT(thisIsParent());
    ddOps[level].applyOp_DoIt(outYH,inYT,fluxYH,whichApp,updateCoefs,getAlpha,alpha);
}

void
DDOp::applyOp_DoIt(MultiFab&         outYH,
                   const MultiFab&   inYT,
                   PArray<MultiFab>& fluxYH,
                   DD_ApForTorRH     whichApp,
                   bool              updateCoefs,
                   bool              getAlpha,
                   MultiFab*         alpha)
{
    const int Nspec = chem->numSpecies();
    const int nGrow = 1;

    int nc = Nspec+1;
    BL_ASSERT(outYH.nComp() >= nc);
    BL_ASSERT(inYT.nComp() >= nc);
    BL_ASSERT(outYH.boxArray() == grids);
    BL_ASSERT(inYT.boxArray() == grids);
    BL_ASSERT(inYT.nGrow() >= nGrow);

    // Need grow cells in X,T to compute forcing
    // Promise to change only the grow cells in T,Y
    int sCompY = 0;
    int sCompT = Nspec;

    setGrowCells(const_cast<MultiFab&>(inYT),sCompT,sCompY);

    if (true)
    {
        for (MFIter mfi(inYT); mfi.isValid(); ++mfi)
        {
            const FArrayBox& fab = inYT[mfi];
            Box              bx  = BoxLib::grow(mfi.validbox(),1);

            for (IntVect p = bx.smallEnd(); p <= bx.bigEnd(); bx.next(p))
            {
                if (isnan(fab(p)))
                {
                    std::cout << "Got a nan at " << p << '\n';
                }
            }
        }
    }

    // Initialize output
    outYH.setVal(0,0,nc);
    for (int d=0; d<BL_SPACEDIM; ++d) {
        fluxYH[d].setVal(0,0,nc);
    }

    const IntVect iv=IntVect::TheZeroVector();
    FArrayBox dum(Box(iv,iv),1);
    if (getAlpha) {
        BL_ASSERT(alpha->nGrow()>=1);
        alpha->setVal(0);
    }

    int for_T0_H1 = (whichApp==DD_Temp ? 0 : 1);
    int Full0_Mix1 = (transport_model == DD_Model_Full ? 0 : 1);

    FArrayBox FcpDTc, CPic, Xc;
    FArrayBox Hic(Box(iv,iv),1);
    FArrayBox FcpDTe(Box(iv,iv),1);

    if (0 && updateCoefs && ParallelDescriptor::IOProcessor()) {
        std::cout << "DDOp::apply: Setting coefficients " << std::endl;
    }

    // HACK
    bool do_abort = false;
    if (0 && mg_level==3)
    {
        std::string str = "DDOpMARC";
        cout << "writing DDOp special" << endl;
        (*mg_parent).Write(str);
    }
    for (MFIter mfi(inYT); mfi.isValid(); ++mfi)
    {
        FArrayBox& outYHc = outYH[mfi];
        const FArrayBox& YTc = inYT[mfi];
        const FArrayBox& c = coefs[mfi];
        const Box& box = mfi.validbox();

        // Actually only need this if for_T0_H1 == 1
        FcpDTc.resize(box,1);
        FcpDTc.setVal(0);

        Box gbox = Box(box).grow(nGrow);
        CPic.resize(gbox,Nspec);
        chem->getCpGivenT(CPic,YTc,gbox,sCompT,0);

        // Actually only need this if for_T0_H1 == 1
        Hic.resize(gbox,Nspec);
        chem->getHGivenT(Hic,YTc,gbox,Nspec,0);

        Xc.resize(gbox,Nspec);
        chem->massFracToMoleFrac(Xc,YTc,gbox,sCompY,0);

        if (updateCoefs) {
            setCoefficients(mfi,YTc,CPic);
        }

        int fillAlpha = getAlpha;
        FArrayBox& alfc = (fillAlpha ? (*alpha)[mfi] : dum);

        IntVect pp(D_DECL(8,0,8));
        if (0 && box.contains(pp))
        {
            IntVect ivs(D_DECL(6,0,8));
            IntVect ive(D_DECL(9,1,9));

            if (box==Box(ivs,ive))
            {
                std::ofstream osf;
                osf.open("junk");
                YTc.writeOn(osf);
                osf.close();

                cout << "************************************************found it" << endl;
                cout << "************************************************level: " << mg_level << endl;

                int model_DD0_MA1 = (transport_model==DD_Model_Full ? 0 : 1);
                int nCoef = FORT_DDNCOEFS(model_DD0_MA1);
                int nState = Nspec+1;
                IntVect iv(D_DECL(8,0,8));
                Box ibox(IntVect(D_DECL(-1,-1,-1)),IntVect(D_DECL(+1,+1,+1)));
                for (IntVect idx=ibox.smallEnd(); idx<=ibox.bigEnd(); ibox.next(idx))
                {
                    IntVect ivp=iv+idx;
                    cout << "pt: " << ivp << endl;
                    cout << "YT: ";
                    for (int n=0; n<nState; ++n) cout << YTc(ivp,n) << " (" << n << ") ";
                    cout << endl;
                    cout << "coef: ";
                    for (int n=0; n<nCoef; ++n) cout << c(ivp,n) << " (" << n << ") ";
                    cout << endl;
                }
                do_abort = true;
            }
        }
        
        for (int dir=0; dir<BL_SPACEDIM; ++dir) {

            const Box ebox = BoxLib::surroundingNodes(box,dir);

            // Actually only need this if for_T0_H1 == 1
            FcpDTe.resize(ebox,Nspec);
            
            // Returns fluxes (and Fn.cpn.gradT if for T)
            FArrayBox& fe = fluxYH[dir][mfi];
            const FArrayBox& ae = area[dir][mfi];

            FORT_DDFLUX(box.loVect(), box.hiVect(), dx.dataPtr(), &dir,
                        fe.dataPtr(), ARLIM(fe.loVect()), ARLIM(fe.hiVect()),
                        FcpDTe.dataPtr(), ARLIM(FcpDTe.loVect()),  ARLIM(FcpDTe.hiVect()),
                        YTc.dataPtr(), ARLIM(YTc.loVect()),  ARLIM(YTc.hiVect()),
                        Xc.dataPtr(), ARLIM(Xc.loVect()),  ARLIM(Xc.hiVect()),
                        c.dataPtr(), ARLIM(c.loVect()),  ARLIM(c.hiVect()),
                        CPic.dataPtr(), ARLIM(CPic.loVect()),  ARLIM(CPic.hiVect()),
                        ae.dataPtr(), ARLIM(ae.loVect()),  ARLIM(ae.hiVect()),
                        &for_T0_H1, Hic.dataPtr(), ARLIM(Hic.loVect()), ARLIM(Hic.hiVect()),
                        &fillAlpha, alfc.dataPtr(), ARLIM(alfc.loVect()), ARLIM(alfc.hiVect()),
                        Full0_Mix1);

            // If for T, increment running sum on cell centers with -avg(F.cp.gT) across faces (in FcpDTe).
            // If for H, q was incremented with +He.Fe inside DDFLUX, nothing more to do
            if (for_T0_H1 == 0) {                
                const int diff0_avg1 = 1;
                const Real a = -1;
                const int oc = 1;
                FORT_DDETC(box.loVect(),box.hiVect(),
                           FcpDTc.dataPtr(),ARLIM(FcpDTc.loVect()),ARLIM(FcpDTc.hiVect()),
                           FcpDTe.dataPtr(),ARLIM(FcpDTe.loVect()),ARLIM(FcpDTe.hiVect()),
                           &a, &dir, &oc, &diff0_avg1);
            }
            
            // Now form -Div(F.Area) add to running total
            const int diff0_avg1 = 0;
            const Real a = -1.;
            FORT_DDETC(box.loVect(),box.hiVect(),
                       outYHc.dataPtr(),ARLIM(outYHc.loVect()),ARLIM(outYHc.hiVect()),
                       fe.dataPtr(),ARLIM(fe.loVect()),ARLIM(fe.hiVect()), &a, &dir, &nc, &diff0_avg1);            
        }
        // Form -(1/Vol) Div(F.Area)
        for (int n=0; n<nc; ++n) {
            outYHc.mult(volInv[mfi],0,n,1);
        }

        // For rho.DT/Dt, form -(1/Vol) Div(F.Area) - avg(F.cp.DT)
        if (for_T0_H1 == 0) {
            outYHc.plus(FcpDTc,0,Nspec,1);
        }

        // Build 1/cpb (CPic <- CPic.Y, CPic_0=sum(CPic_n)
        CPic.mult(YTc,0,0,Nspec);
        for (int n=1; n<Nspec; ++n) {
            CPic.plus(CPic,n,0,1);
        }
        CPic.invert(1,0,1);
        
        // For rho.DT/Dt, form (1/Cp) [ -(1/Vol) Div(F.Area) - avg(F.cp.DT) ]
        if (for_T0_H1 == 0) {
            outYHc.mult(CPic,0,Nspec,1);
        }

        if (getAlpha) {
            for (int n=0; n<nc; ++n) {
                alfc.mult(volInv[mfi],0,n,1);
            }
            alfc.mult(CPic,0,Nspec,1);

            // Modify coefficient based on what was needed to fill adjacent grow cells
            for (OrientationIter oitr; oitr; ++oitr) {
                const Orientation face = oitr();
                const FabSet& stfs = stencilWeight[face];
                int dir = face.coordDir();
                int shiftCells = ( face.isLow() ? -1 : +1 ); 

                const FArrayBox& src = stfs[mfi];
                const Box& srcBox = src.box();
                const Box dstBox = Box(srcBox).shift(dir,shiftCells); 
                    
                // Do T component
                {
                    int sComp = 1; // Stencil coef for T in slot 1
                    int dComp = Nspec;
                    alfc.mult(src,srcBox,dstBox,sComp,dComp,1);
                    alfc.plus(alfc,dstBox,srcBox,dComp,dComp,1);
                }
                // Do Y component
                {
                    int sComp = 0; // stencil coef for all Y in slot 0
                    for (int i=0; i<Nspec; ++i) {
                        int dComp = i;
                        alfc.mult(src,srcBox,dstBox,sComp,dComp,1);
                        alfc.plus(alfc,dstBox,srcBox,dComp,dComp,1);
                    }
                }
            }
        }
    }
    ParallelDescriptor::Barrier();
    if (do_abort)
        BoxLib::Abort();
    ParallelDescriptor::Barrier();
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
    for (MFIter mfi(mfC); mfi.isValid(); ++mfi)
    {
        FArrayBox& F = mfF[mfi];
        const FArrayBox& C = mfC[mfi];
        const Box cbox = mfi.validbox();
        FORT_DDCCINT(F.dataPtr(dCompF),ARLIM(F.loVect()), ARLIM(F.hiVect()),
                     C.dataPtr(sCompC),ARLIM(C.loVect()), ARLIM(C.hiVect()),
                     cbox.loVect(), cbox.hiVect(), &nComp, MGIV.getVect());
    }
}
#include "Utility.H"
#ifdef WIN32
static std::string sep = "\\";
#else
static std::string sep = "/";
#endif

void
DDOp::Write(std::string& outfile) const 
{
    if (mg_parent!=this)
    {
        BoxLib::Abort("DDOp::Write: only the parent can write the DDOp structure");
    }

    // Set filenames
    std::string DDOpNAME = outfile + sep + "DDOp" + sep + "MG_Level_";
    std::string DDBndryTNAME = outfile + sep + "DDBndryT" + sep + "MG_Level_";
    std::string DDBndryYNAME = outfile + sep + "DDBndryY" + sep + "MG_Level_";
    Array<std::string> DDOpNAMElev(ddOps.size());
    Array<std::string> DDBndryTNAMElev(ddOps.size());
    Array<std::string> DDBndryYNAMElev(ddOps.size());

    for (int i=0;i<ddOps.size(); ++i)
    {
        DDOpNAMElev[i] = BoxLib::Concatenate(DDOpNAME,i,2);
        DDBndryTNAMElev[i] = BoxLib::Concatenate(DDBndryTNAME,i,2);
        DDBndryYNAMElev[i] = BoxLib::Concatenate(DDBndryYNAME,i,2);
    }


    if (ParallelDescriptor::IOProcessor()) {
        if (!BoxLib::UtilCreateDirectory(outfile, 0755))
            BoxLib::CreateDirectoryFailed(outfile);

        std::string hdrName("Header");
        hdrName = outfile + sep + hdrName;
        std::ofstream ofs;
        ofs.open(hdrName.c_str());
        ofs << BL_SPACEDIM << '\n';
        ofs << geom.Domain() << '\n';
        ofs << grids << '\n';
        ofs << amr_ratio << '\n';
        std::string modelName = (transport_model==DD_Model_Full ? DD_FULL_NAME : DD_MIX_NAME);
        ofs << modelName << '\n';
        ofs << maxorder << '\n';
        ofs << mgLevels << '\n';
        ofs << mgLevelsMAX << '\n';    
        ofs.close();
    }
    ParallelDescriptor::Barrier();

    for (int i=0;i<ddOps.size(); ++i)
    {
        ddOps[i].WriteSub(DDOpNAMElev[i]);
        Tbd->Write(DDBndryTNAMElev[i]);
        Ybd->Write(DDBndryYNAMElev[i]);
    }
}

void
DDOp::WriteSub(std::string& outfile, int level) const 
{
    if (mg_parent!=this)
    {
        BoxLib::Abort("DDOp::Write: only the parent can write");
    }
    if (level>=0  && level<mgLevels)
    {
        BL_ASSERT(ddOps.defined(level));
        ddOps[level].WriteSub(outfile);
    }
}

void
DDOp::WriteSub(std::string& outfile) const 
{
    if (ParallelDescriptor::IOProcessor()) {
        if (!BoxLib::UtilCreateDirectory(outfile, 0755))
            BoxLib::CreateDirectoryFailed(outfile);
        
        std::ofstream osf;
        std::string hdrName = "Header";
        hdrName = outfile + sep + hdrName;
        std::ofstream ofs;
        ofs.open(hdrName.c_str());
        ofs << grids << '\n';
        ofs.close();
    }
    ParallelDescriptor::Barrier();

    VisMF::Write(volInv,outfile+sep+"volInv");
    VisMF::Write(coefs,outfile+sep+"coefs");
    int mindigits = 1;
    char buf[32];
    for (int i=0; i<BL_SPACEDIM; ++i) {
        sprintf(buf, "%0*d",  mindigits, i);
        VisMF::Write(area[i],outfile+sep+"area_"+std::string(buf));
    }
    for (OrientationIter oitr; oitr; ++oitr) {
        sprintf(buf, "%0*d",  mindigits, (int)oitr());
        stencilWeight[oitr()].write(outfile+sep+"stencilWeight"+std::string(buf));
    }
}
