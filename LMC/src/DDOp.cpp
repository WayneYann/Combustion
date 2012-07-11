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
DDOp::DD_Average DDOp::average_normal = DD_Arithmetic; // Used to move transport coeffs evaluated at cc to edges
DDOp::DD_Average DDOp::average_tangential = DD_Arithmetic; // Used to coarsen edge-based transport coeffients
int DDOp::transport_coefs_nComp = -1; // Invalid until valid transport model identified

DDOp::DDOp ()
    :
    Tbd(0),
    Ybd(0),
    mg_parent(0)
{
}

DDOp::DDOp (const BoxArray&   grids,
            const Box&        box,
            const IntVect&    amrRatio)
    : Tbd(0), Ybd(0)
{
    BL_ASSERT(ddOps.size()==0);
    ddOps.resize(1,PArrayNoManage); // Kinda weird, since ddOps[0] is me.  Will explicitly manage in dtr
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
        for (int i=1; i<ddOps.size(); ++i) {
            delete ddOps.remove(i);
        }
        ddOps.remove(0); // returns "this", don't delete
        ddOps.clear();
    }
    volInv.clear();
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
            ddOps.resize(mgLevels,PArrayNoManage);
            ddOps.set(mgLevels-1,new DDOp(cgrids,cbox,amr_ratio,mgLevels-1,parent));

            not_done = can_coarsen(cgrids)  &&  can_coarsen(cbox) 
                &&  ( (mgLevelsMAX<0) || (mgLevels < mgLevelsMAX) );
        }

        Tbd = new DDBndry(grids,1,geom,mgLevels);
        Ybd = new DDBndry(grids,Nspec,geom,mgLevels);
    }

    dx.resize(BL_SPACEDIM,PArrayManage);
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
    transport_coefs_nComp = FORT_DDNCOEFS(&model_DD0_MA1);
    transport_coefs.resize(BL_SPACEDIM,PArrayManage);
    for (int dir = 0; dir < BL_SPACEDIM; dir++) {
        BoxArray egrids = BoxArray(grids).surroundingNodes(dir);
        transport_coefs.set(dir,new MultiFab(egrids,transport_coefs_nComp,0,Fab_allocate));
    }
    cpi.define(grids,Nspec,0,Fab_allocate);

    // Weights used to generate grow cell data, only two components here because all Ys have to have the same bc
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
DDOp::setGrowCells(MultiFab& T,
                   int       sCompT,
                   MultiFab& Y,
                   int       sCompY)
{
    BL_ASSERT(T.nGrow() >= 1);
    BL_ASSERT(Y.nGrow() >= 1);
    const int Nspec = chem->numSpecies();

    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.nComp() >= sCompY + Nspec);

    BL_ASSERT(T.boxArray() == grids);
    BL_ASSERT(Y.boxArray() == grids);

    const int flagbc  = 1;
    const int flagden = 1;

    DDBndry& Tbndry = (*mg_parent->Tbd);
    DDBndry& Ybndry = (*mg_parent->Ybd);

    T.FillBoundary(sCompT,1);
    geom.FillPeriodicBoundary(T,sCompT,1);
    Y.FillBoundary(sCompY,Nspec);
    geom.FillPeriodicBoundary(Y,sCompY,Nspec);

    for (MFIter mfi(T); mfi.isValid(); ++mfi)
    {
        BL_ASSERT(!T[mfi].contains_nan(mfi.validbox(),sCompT,1));
        BL_ASSERT(!Y[mfi].contains_nan(mfi.validbox(),sCompY,Nspec));
    }

    for (OrientationIter oitr; oitr; ++oitr)
    {
        const Orientation&      face = oitr();
        const int              iFace = (int)face;

        const FabSet&            Tfs = Tbndry.bndryValues(face,mg_level);
        FabSet&                  Tst = stencilWeight[face];
        const int                Tnc = 1;

        const FabSet&            Yfs = Ybndry.bndryValues(face,mg_level);
        FabSet&                  Yst = stencilWeight[face];
        const int                Ync = Nspec;

        const int comp = 0;

        for (MFIter mfi(T); mfi.isValid(); ++mfi)
        {
            const int   idx = mfi.index();
            const Box& vbox = mfi.validbox();
            BL_ASSERT(grids[idx] == vbox);

            FArrayBox& Tfab = T[mfi];
            FArrayBox& Tstfab = Tst[mfi];
            const FArrayBox& Tb = Tfs[mfi];

            BL_ASSERT(!Tb.contains_nan());

            const Mask& Tm  = Tbndry.bndryMasks(face,mg_level)[idx];
            const Real Tbcl = Tbndry.bndryLocs(face,idx);
            const int  Tbct = Tbndry.bndryConds(face,idx)[comp];

            FORT_APPLYBC(&flagden, &flagbc, &maxorder,
                         Tfab.dataPtr(sCompT), ARLIM(Tfab.loVect()), ARLIM(Tfab.hiVect()),
                         &iFace, &Tbct, &Tbcl,
                         Tb.dataPtr(), ARLIM(Tb.loVect()), ARLIM(Tb.hiVect()),
                         Tm.dataPtr(), ARLIM(Tm.loVect()), ARLIM(Tm.hiVect()),
                         Tstfab.dataPtr(1), ARLIM(Tstfab.loVect()), ARLIM(Tstfab.hiVect()),
                         vbox.loVect(),vbox.hiVect(), &Tnc, dx.dataPtr());

            if (Tfab.contains_nan(Tb.box(),sCompT,1))
            {
                std::cout << "Tfab contains NaNs: " << Tfab.box() << ' ' << Tb.box() << '\n';

                std::cout << "Tbcl: " << Tbcl << ", Tbct: " << Tbct << '\n';

                for (IntVect p = Tb.box().smallEnd(); p <= Tb.box().bigEnd(); Tb.box().next(p))
                {
                    //if (std::isnan(Tfab(p,sCompT)))
                        //std::cout << "T isnan @ p = " << p << '\n';
                }
            }

            BL_ASSERT(!Tfab.contains_nan(Tb.box(),sCompT,1));

            FArrayBox& Yfab = Y[mfi];
            FArrayBox& Ystfab = Yst[mfi];
            const FArrayBox& Yb = Yfs[mfi];

            BL_ASSERT(!Yb.contains_nan());

            const Mask& Ym  = Ybndry.bndryMasks(face,mg_level)[idx];
            const Real Ybcl = Ybndry.bndryLocs(face,idx);
            const int Ybct  = Ybndry.bndryConds(face,idx)[comp];

            FORT_APPLYBC(&flagden, &flagbc, &maxorder,
                         Yfab.dataPtr(sCompY), ARLIM(Yfab.loVect()), ARLIM(Yfab.hiVect()),
                         &iFace, &Ybct, &Ybcl,
                         Yb.dataPtr(), ARLIM(Yb.loVect()), ARLIM(Yb.hiVect()),
                         Ym.dataPtr(), ARLIM(Ym.loVect()), ARLIM(Ym.hiVect()),
                         Ystfab.dataPtr(0), ARLIM(Ystfab.loVect()), ARLIM(Ystfab.hiVect()),
                         vbox.loVect(),vbox.hiVect(), &Ync, dx.dataPtr());

            BL_ASSERT(!Yfab.contains_nan(Yb.box(),sCompY,Ync));
        }
    } // orientation

    T.FillBoundary(sCompT,1);
    geom.FillPeriodicBoundary(T,sCompT,1);
    Y.FillBoundary(sCompY,Nspec);
    geom.FillPeriodicBoundary(Y,sCompY,Nspec);
}

void
DDOp::setCoefficients(const MultiFab& T,
                      int             sCompT,
                      const MultiFab& Y,
                      int             sCompY)
{
    // Only to be called at mg_level=0 currently, crse data generated automatically
    //
    // This computes transport coefs (and cp) on cell centers, including grow cells
    // and then generates edge-based coef data for mg_level=0.  Then coarser data is 
    // generated by averaging edge data.
    BL_ASSERT(thisIsParent());

    int nGrow = 1;
    int Nspec = chem->numSpecies();
    BL_ASSERT(sCompT<T.nComp());
    BL_ASSERT(sCompY+Nspec<Y.nComp());
    BL_ASSERT(T.boxArray()==grids);
    BL_ASSERT(Y.boxArray()==grids);
    BL_ASSERT(T.nGrow()>=nGrow);
    BL_ASSERT(Y.nGrow()>=nGrow);
    BL_ASSERT(transport_coefs_nComp>0);

    // Do const_cast...promise to change only grow cells
#ifndef NDEBUG
    const_cast<MultiFab&>(T).setBndry(0.,sCompT,1);
    const_cast<MultiFab&>(Y).setBndry(0.,sCompY,Nspec);
#endif
    setGrowCells(const_cast<MultiFab&>(T),sCompT,const_cast<MultiFab&>(Y),sCompY);

    FArrayBox tcpfab, tcoefs;
    for (MFIter mfi(T); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const Box gbox = Box(box).grow(nGrow);
        int Full0_Mix1 = (transport_model == DD_Model_Full ? 0 : 1);
        const FArrayBox& Tfab = T[mfi];
        const FArrayBox& Yfab = Y[mfi];
        
        tcpfab.resize(gbox,Nspec);
        BL_ASSERT(transport_coefs_nComp>0);
        tcoefs.resize(gbox,transport_coefs_nComp);

        BL_ASSERT(!Tfab.contains_nan(gbox,sCompT,1));
        chem->getCpGivenT(tcpfab,Tfab,gbox,sCompT,0);
        BL_ASSERT(!tcpfab.contains_nan(gbox,0,Nspec));
        cpi[mfi].copy(tcpfab); // will only need valid-region cp afterward

        BL_ASSERT(!Tfab.contains_nan(gbox,sCompT,1));
        BL_ASSERT(!Yfab.contains_nan(gbox,sCompY,Nspec));
        FORT_DDCOEFS(gbox.loVect(), gbox.hiVect(),
                     tcoefs.dataPtr(), ARLIM(tcoefs.loVect()), ARLIM(tcoefs.hiVect()),
                     Tfab.dataPtr(sCompT), ARLIM(Tfab.loVect()),   ARLIM(Tfab.hiVect()),
                     Yfab.dataPtr(sCompY), ARLIM(Yfab.loVect()),   ARLIM(Yfab.hiVect()),
                     tcpfab.dataPtr(), ARLIM(tcpfab.loVect()), ARLIM(tcpfab.hiVect()),
                     &Full0_Mix1);
        BL_ASSERT(!tcoefs.contains_nan(gbox,0,transport_coefs_nComp));

        for (int dir=0; dir<BL_SPACEDIM; ++dir)
        {
            int avgNormal_H0_A1 = (average_normal==DD_Harmonic ? 0 : 1);
            center_to_edge(tcoefs,transport_coefs[dir][mfi],box,0,0,transport_coefs_nComp,avgNormal_H0_A1);
        }
    }

    // Now, coarsen data
    for (int lev=1; lev<mgLevels; ++lev)
    {
        BL_ASSERT(ddOps.size() > lev);
        const DDOp& fine = ddOps[lev-1];
        DDOp&       crse = ddOps[lev];
        average(crse.cpi,0,fine.cpi,0,fine.cpi.nComp());
        for (int dir=0; dir<BL_SPACEDIM; ++dir)
        {
            int avgTang_H0_A1 = (average_tangential==DD_Harmonic ? 0 : 1);
            average_ec(crse.transport_coefs[dir],0,fine.transport_coefs[dir],0,transport_coefs_nComp,dir,avgTang_H0_A1);
        }
    }
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
              int               level,
              bool              getAlpha,
              MultiFab*         alpha)
{
    BL_ASSERT(level >= 0  &&  ((mgLevelsMAX<0) || (level<mgLevelsMAX)) );
    BL_ASSERT(thisIsParent());
    ddOps[level].applyOp_DoIt(outYH,inYT,fluxYH,whichApp,getAlpha,alpha);
}

void
DDOp::applyOp_DoIt(MultiFab&         outYH,
                   const MultiFab&   inYT,
                   PArray<MultiFab>& fluxYH,
                   DD_ApForTorRH     whichApp,
                   bool              getAlpha,
                   MultiFab*         alpha)
{
    //
    //  Compute - (1/V)Div(A.Fi) for Y components and -(1/V)Div(A.Q) for H
    //   (if whichApp==Temp, compute -(1/V)Div(Q-Fi.hi) - Fi.cpi.Grad(T) )
    //  And return diagnonal coefficient of the discrete operator.
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
    // FIXME: Should probably pass these in, as they assume the layout of inYT
    //   Also then enforce the layout of outYH, fluxYH and alpha
    int sCompY = 0;
    int sCompT = Nspec;

    setGrowCells(const_cast<MultiFab&>(inYT),sCompT,const_cast<MultiFab&>(inYT),sCompY);

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

    FArrayBox sumFcpDTc, Xc, F_dTe, F_dTc;
    FArrayBox Hic(Box(iv,iv),1);

    bool bad_skip=false;
    for (MFIter mfi(inYT); mfi.isValid(); ++mfi)
    {
        FArrayBox& outYHc = outYH[mfi];
        const FArrayBox& YTc = inYT[mfi];
        const Box& box = mfi.validbox();
        Box gbox = Box(box).grow(nGrow);

#ifndef NDEBUG
        BoxArray myRegion(2*BL_SPACEDIM+1);
        myRegion.set(0,box);
        int cnt=0;
        for (OrientationIter oitr; oitr; ++oitr)
        {
            myRegion.set(++cnt,BoxLib::adjCell(box,oitr(),inYT.nGrow()));
        }
        BoxArray corners = BoxLib::complementIn(gbox,myRegion);
        //BoxArray corners = BoxLib::complementIn(gbox,grids);
#endif

        Hic.resize(gbox,Nspec);
#ifndef NDEBUG
        for (int i=0; i<corners.size(); ++i) {
            const_cast<FArrayBox&>(YTc).setVal(0,corners[i],sCompT,1);
            const_cast<FArrayBox&>(YTc).setVal(0,corners[i],sCompY,Nspec);
        }
#endif
        BL_ASSERT(!YTc.contains_nan(gbox,sCompT,1));
        BL_ASSERT(!YTc.contains_nan(gbox,sCompY,Nspec));
        chem->getHGivenT(Hic,YTc,gbox,sCompT,0);
#ifndef NDEBUG
        for (int i=0; i<corners.size(); ++i)
            Hic.setVal(0,corners[i],0,Nspec);
#endif
        // FIXME
#ifndef NDEBUG
        if (Hic.contains_nan(gbox,0,Nspec)) {
            cout << "Hic" << endl;
            cout << "box,id: " << box << " " << mfi.index() << endl;
            cout << "gbox: " << gbox << endl;
            BoxArray bat = BoxLib::complementIn(gbox,outYH.boxArray());
            cout << "uncovered: " << bat;
            for (int k=0; k<bat.size(); ++k) {
                FArrayBox tf(bat[k],Nspec);
                tf.copy(Hic);
                if (tf.contains_nan()) {
                    FArrayBox tt(tf.box(),1);
                    tt.copy(YTc,sCompT,0,1);
                    cout << tt << endl;
                    cout << tf << endl;
                }
                else
                {
                    cout << bat[k] << " is ok" << endl;
                }
            }
            bad_skip = true;
            continue;
        }
#endif
        BL_ASSERT(!Hic.contains_nan(gbox,0,Nspec));

        Xc.resize(gbox,Nspec);
        chem->massFracToMoleFrac(Xc,YTc,gbox,sCompY,0);
#ifndef NDEBUG
        for (int i=0; i<corners.size(); ++i)
            Xc.setVal(0,corners[i],0,Nspec);
#endif
        BL_ASSERT(!Xc.contains_nan(gbox,0,Xc.nComp()));

        int fillAlpha = getAlpha;
        FArrayBox& alfc = (fillAlpha ? (*alpha)[mfi] : dum);

        F_dTc.resize(box,Nspec);
        F_dTc.setVal(0);

        for (int dir=0; dir<BL_SPACEDIM; ++dir) {

            const Box ebox = BoxLib::surroundingNodes(box,dir);

            // Actually only need this if for_T0_H1 == 1, but resize here so dataPtr call below wont fail...
            F_dTe.resize(ebox,Nspec);
            
            FArrayBox& fe = fluxYH[dir][mfi];
            const FArrayBox& ae = area[dir][mfi];

            const FArrayBox& c = transport_coefs[dir][mfi];
            BL_ASSERT(!c.contains_nan(box,0,c.nComp()));

            // Compute Fi.Area, and Q.Area for RhoH,  or (Q-Fi.hi).Area and Fi.cp.Grad(T) for Temp
            FORT_DDFLUX(box.loVect(), box.hiVect(), dx.dataPtr(), &dir,
                        fe.dataPtr(), ARLIM(fe.loVect()), ARLIM(fe.hiVect()),
                        F_dTe.dataPtr(), ARLIM(F_dTe.loVect()),  ARLIM(F_dTe.hiVect()),
                        YTc.dataPtr(), ARLIM(YTc.loVect()),  ARLIM(YTc.hiVect()),
                        Xc.dataPtr(), ARLIM(Xc.loVect()),  ARLIM(Xc.hiVect()),
                        c.dataPtr(), ARLIM(c.loVect()),  ARLIM(c.hiVect()),
                        ae.dataPtr(), ARLIM(ae.loVect()),  ARLIM(ae.hiVect()),
                        &for_T0_H1, Hic.dataPtr(), ARLIM(Hic.loVect()), ARLIM(Hic.hiVect()),
                        &fillAlpha, alfc.dataPtr(), ARLIM(alfc.loVect()), ARLIM(alfc.hiVect()),
                        &Full0_Mix1);

            BL_ASSERT(!fe.contains_nan(ebox,0,fe.nComp()));
            
            BL_ASSERT(!F_dTe.contains_nan(ebox,0,for_T0_H1 == 0 ? F_dTe.nComp() : 1));

            // If for T, increment running sum on cell centers with -avg(F.grad(T)) across faces (in F_dTe).
            // If for H, Q.Area was incremented with sum(Fi.hi).Area inside DDFLUX, nothing more to do
            if (for_T0_H1 == 0) {
                const int diff0_avg1 = 1;
                const Real a = -1;
                FORT_DDETC(box.loVect(),box.hiVect(),
                           F_dTc.dataPtr(),ARLIM(F_dTc.loVect()),ARLIM(F_dTc.hiVect()),
                           F_dTe.dataPtr(),ARLIM(F_dTe.loVect()),ARLIM(F_dTe.hiVect()),
                           &a, &dir, &Nspec, &diff0_avg1);
            }
            
            // Now form -Div(F.Area), add to running total
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
        if (for_T0_H1 == 0)
        {
            F_dTc.mult(cpi[mfi],0,0,Nspec);

            // Compute Sum(F.cp.dT) 
            for (int i=1; i<Nspec; ++i) {
                F_dTc.plus(F_dTc,i,0,1);
            }

            // Increment result
            outYHc.plus(F_dTc,0,Nspec,1);
        }

        if (getAlpha) {

            // Scale alpha by -(1/Vol) to be consistent with construction above 
            alfc.mult(-1,0,nc);
            for (int n=0; n<nc; ++n) {
                alfc.mult(volInv[mfi],0,n,1);
            }
            // For RhoH equation, scale lambda by 1/cpmix_cc (use Xc container since its already about the correct size)
            if (for_T0_H1 == 1) {
                FArrayBox& cpmix = Xc;
                cpmix.copy(cpi[mfi],box,0,box,0,Nspec);
                cpmix.mult(YTc,box,box,sCompY,0,Nspec);
                for (int n=1; n<Nspec; ++n)
                    cpmix.plus(cpmix,box,box,n,0,1);
                cpmix.invert(1,box,0,1);
                alfc.mult(cpmix,box,box,0,sCompT,1);
            }

            // Modify coefficient based on what was needed to fill adjacent grow cells
            for (OrientationIter oitr; oitr; ++oitr) {
                const Orientation face = oitr();
                if (geom.isPeriodic(oitr().coordDir())) continue;

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
            BL_ASSERT(!alfc.contains_nan(box,0,Nspec));
            BL_ASSERT(!alfc.contains_nan(box,Nspec,1));
        }
#if 0
        if (outYHc.contains_nan(box,0,Nspec))
        {
            for (int i=0; i<Nspec; ++i) {
                if (outYHc.contains_nan(box,i,1))
                    cout << "component " << i << " contains a nan" << endl;
            }
        }
#endif
        BL_ASSERT(!outYHc.contains_nan(box,0,Nspec));
        BL_ASSERT(!outYHc.contains_nan(box,Nspec,1));
    }
    if (bad_skip) BoxLib::Abort();
}

void
DDOp::average_ec(MultiFab&       mfC,
                 int             dCompC,
                 const MultiFab& mfF,
                 int             sCompF,
                 int             nComp,
                 int             dir,
                 int             avgTang_H0_A1)
{
    BL_ASSERT(mfC.nComp()>=dCompC+nComp);
    BL_ASSERT(mfF.nComp()>=sCompF+nComp);
    const int* ratio = MGIV.getVect();
    for (MFIter mfi(mfC); mfi.isValid(); ++mfi)
    {
        FArrayBox& C = mfC[mfi];
        const FArrayBox& F = mfF[mfi];
        const Box& ebox = mfi.validbox();

        FORT_DDECAVG(ebox.loVect(),ebox.hiVect(),
                     C.dataPtr(), ARLIM(C.loVect()), ARLIM(C.hiVect()),
                     F.dataPtr(), ARLIM(F.loVect()), ARLIM(F.hiVect()),
                     &nComp, ratio, &dir, &avgTang_H0_A1);
    }
}

void
DDOp::center_to_edge (const FArrayBox& cfab,
                      FArrayBox&       efab,
                      const Box&       ccBox,
                      int              sComp,
                      int              dComp,
                      int              nComp,
                      int              avgNormal_H0_A1)
{
    // Compute data on edges between each pair of cc values in the dir direction
    const Box&      ebox = efab.box();
    const IndexType ixt  = ebox.ixType();

    BL_ASSERT(!(ixt.cellCentered()) && !(ixt.nodeCentered()));

    int dir = -1;
    for (int d = 0; d < BL_SPACEDIM; d++)
        if (ixt.test(d))
            dir = d;

    BL_ASSERT(ccBox.contains(BoxLib::enclosedCells(ebox)));
    BL_ASSERT(cfab.box().contains(Box(ccBox).grow(dir,1)));
    BL_ASSERT(sComp+nComp <= cfab.nComp() && dComp+nComp <= efab.nComp());

    FORT_DDC2E(ccBox.loVect(),ccBox.hiVect(),
               cfab.dataPtr(sComp),ARLIM(cfab.loVect()),ARLIM(cfab.hiVect()),
               efab.dataPtr(dComp),ARLIM(efab.loVect()),ARLIM(efab.hiVect()),
               &nComp, &dir, &avgNormal_H0_A1);
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
    const int* ratio = MGIV.getVect();
    for (MFIter mfi(mfC); mfi.isValid(); ++mfi)
    {
        FArrayBox& C = mfC[mfi];
        const FArrayBox& F = mfF[mfi];
        const Box& cbox = mfi.validbox();
        BL_ASSERT(cbox == BoxLib::coarsen(mfF.boxArray()[mfi.index()],MGIV));
        FORT_DDCCAVG(cbox.loVect(), cbox.hiVect(), 
                     C.dataPtr(dCompC),ARLIM(C.loVect()), ARLIM(C.hiVect()),
                     F.dataPtr(sCompF),ARLIM(F.loVect()), ARLIM(F.hiVect()),
                     &nComp, ratio);
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
        const Box& cbox = mfi.validbox();
        BL_ASSERT(F.box().contains(BoxLib::refine(cbox,MGIV)));
        FORT_DDCCINT(cbox.loVect(), cbox.hiVect(),
                     F.dataPtr(dCompF),ARLIM(F.loVect()), ARLIM(F.hiVect()),
                     C.dataPtr(sCompC),ARLIM(C.loVect()), ARLIM(C.hiVect()),
                     &nComp, MGIV.getVect());
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
    Array<std::string> DDOpNAMElev(ddOps.size());

    for (int i=0;i<ddOps.size(); ++i)
    {
        DDOpNAMElev[i] = BoxLib::Concatenate(DDOpNAME,i,2);
    }


    bool verbose = true;
    if (ParallelDescriptor::IOProcessor()) {
        if (!BoxLib::UtilCreateDirectory(outfile, 0755))
            BoxLib::CreateDirectoryFailed(outfile);

        std::string hdrName("Header");
        hdrName = outfile + sep + hdrName;
        std::ofstream ofs;
        if (verbose)
            cout << "DDOp::Write " << ": writing Header" << endl;
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
        if (verbose && ParallelDescriptor::IOProcessor())
            cout << "DDOp::Write: writing mg_level " << i << " to " << DDOpNAMElev[i] << endl;
        ddOps[i].WriteSub(DDOpNAMElev[i]);
    }

    std::string DDBndryTNAME = outfile + sep + "DDBndryT";
    if (verbose && ParallelDescriptor::IOProcessor())
        cout << "DDOp::Write: writing DDBndryT to " << DDBndryTNAME << endl;
    Tbd->Write(DDBndryTNAME);

    std::string DDBndryYNAME = outfile + sep + "DDBndryY";
    if (verbose && ParallelDescriptor::IOProcessor())
        cout << "DDOp::Write: writing DDBndryY to " << DDBndryYNAME << endl;
    Ybd->Write(DDBndryYNAME);
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
    bool verbose = true;
    if (verbose && ParallelDescriptor::IOProcessor()) {
        if (!BoxLib::UtilCreateDirectory(outfile, 0755))
            BoxLib::CreateDirectoryFailed(outfile);
        
        std::ofstream osf;
        std::string hdrName = "Header";
        hdrName = outfile + sep + hdrName;
        std::ofstream ofs;
        if (verbose)
            cout << "DDOp::WriteSub mg_level " << mg_level << ": writing Header" << endl;
        ofs.open(hdrName.c_str());
        ofs << grids << '\n';
        ofs.close();
    }
    ParallelDescriptor::Barrier();

    if (verbose && ParallelDescriptor::IOProcessor())
        cout << "DDOp::WriteSub mg_level " << mg_level << ": writing volInv" << endl;
    VisMF::Write(volInv,outfile+sep+"volInv");
    if (verbose && ParallelDescriptor::IOProcessor())
        cout << "DDOp::WriteSub mg_level " << mg_level << ": writing cpi" << endl;
    VisMF::Write(cpi,outfile+sep+"cpi");

    for (int i=0; i<BL_SPACEDIM; ++i) {
        std::string dim = BoxLib::Concatenate("", i, 1);
        if (verbose && ParallelDescriptor::IOProcessor())
            cout << "DDOp::WriteSub mg_level " << mg_level << ": writing area " << i << endl;
        VisMF::Write(area[i],outfile+sep+"area_"+dim);
        if (verbose && ParallelDescriptor::IOProcessor())
            cout << "DDOp::WriteSub mg_level " << mg_level << ": writing transport_coefs " << i << endl;
        VisMF::Write(transport_coefs[i],outfile+sep+"transport_coefs_"+dim);
    }
    for (OrientationIter oitr; oitr; ++oitr) {
        std::string face = BoxLib::Concatenate("", oitr(), 1);
        if (verbose && ParallelDescriptor::IOProcessor())
            cout << "DDOp::WriteSub mg_level " << mg_level << ": writing stencil weights " << oitr() << endl;
        stencilWeight[oitr()].write(outfile+sep+"stencilWeight"+face);
    }
}
