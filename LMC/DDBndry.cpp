//BL_COPYRIGHT_NOTICE

//
// $Id: DDBndry.cpp,v 1.6 2011-02-02 01:36:44 marc Exp $
//
#include <winstd.H>
#include <cmath>

#include <LO_BCTYPES.H>
#include <ArrayLim.H>
#include <DDBndry.H>
#include <DDBndry_F.H>

using std::cout;
using std::endl;
#include <fstream>

const IntVect MGIV = IntVect(D_DECL(2,2,2));
DDBndry::DDBndry (const BoxArray& grids,
                  int             ncomp,
                  const Geometry& geom,
                  int             numMGlevels)
{
    define(grids,ncomp,geom,numMGlevels);
}

void
DDBndry::define (const BoxArray& grids,
                 int             ncomp,
                 const Geometry& geom,
                 int             numMGlevels)
{
    BL_ASSERT(numMGlevels > 0);
    ViscBndry::define(grids,ncomp,geom);
    num_mg_levels = numMGlevels;
    
    bnd_vals.resize(num_mg_levels);
    geoms.resize(num_mg_levels);

    geoms.set(0,new Geometry(geom));
    BoxArray cgrids = grids;
    BoxArray fgrids;

    for (int clev=1; clev<num_mg_levels; ++clev) 
    {
        fgrids = cgrids;
        cgrids.coarsen(MGIV);

        Box cdomain = Box(geoms[clev-1].Domain()).coarsen(MGIV);
        geoms.set(clev,new Geometry(cdomain));
        BL_ASSERT(BoxArray(cgrids).refine(MGIV)==fgrids);
        BL_ASSERT(Box(cdomain).refine(MGIV)==geoms[clev-1].Domain());
        bnd_vals.set(clev,new BndryData(cgrids,ncomp,geoms[clev]));
    }
}

void
DDBndry::setBndryValues (::BndryRegister&  crse,
                         int             c_start,
                         const MultiFab& fine,
                         int             f_start,
                         int             bnd_start,
                         int             num_comp,
                         int             ratio, 
                         const BCRec&    phys_bc,
                         int             max_order)
{
    IntVect iv = ratio*IntVect::TheUnitVector();
    setBndryValues(crse,c_start,fine,f_start,bnd_start,num_comp,iv,phys_bc,max_order);
}

void
coarsen_data_in_fabsets(FabSet&            fine,
                        int                sComp,
                        FabSet&            crse,
                        int                dComp,
                        int                nComp,
                        const Orientation& face)
{
    //
    // Too messy to work out indices here.  Try something else (note mult below,
    //  pretty clever, I think...)  Should probably find a way to call DDOp::coarsen
    //  here rather than having two fortran routines that are essentially identical.
    //  
    FArrayBox f_fab, c_fab;
    const int dir = face.coordDir();
    for (FabSetIter mfi(crse); mfi.isValid(); ++mfi)
    {
        const FArrayBox& ffab = fine[mfi];
        FArrayBox& cfab = crse[mfi];
        Box c_bnd_bx = cfab.box();
        Box f_bnd_bx = BoxLib::refine(c_bnd_bx,MGIV);
        
        f_fab.resize(f_bnd_bx,nComp);
        f_fab.setVal(0.);
        Box ovlp = f_bnd_bx & ffab.box();
        f_fab.copy(ffab,ovlp,sComp,ovlp,0,nComp);
        f_fab.mult(MGIV[dir]);
        
        c_fab.resize(c_bnd_bx,nComp);
        FORT_DDBNDCCAVG(c_fab.dataPtr(), ARLIM(c_bnd_bx.loVect()), ARLIM(c_bnd_bx.hiVect()),
                        f_fab.dataPtr(), ARLIM(f_bnd_bx.loVect()), ARLIM(f_bnd_bx.hiVect()),
                        c_bnd_bx.loVect(), c_bnd_bx.hiVect(), &nComp, MGIV.getVect());
        cfab.copy(c_fab,0,dComp,nComp);
    }
}                 

void
DDBndry::setBndryValues (::BndryRegister&  crse,
                         int             c_start,
                         const MultiFab& fine,
                         int             f_start,
                         int             bnd_start,
                         int             num_comp,
                         IntVect&        ratio, 
                         const BCRec&    phys_bc,
                         int             max_order)
{
    //
    // Check that boxarrays are identical.
    //
    BL_ASSERT(grids.size());
    BL_ASSERT(grids == fine.boxArray());
    //
    // Set amr structure, which will be the same for all mg levels
    //
    BoxArray amr_crse_grids = BoxArray(grids).coarsen(ratio);
    if (BoxArray(amr_crse_grids).refine(ratio) != grids)
    {
        BoxLib::Abort("Unaceptable ratio in DDBndry::setBndryValues");
    }

    // Set BR data at level 0, as before
    InterpBndryData::setBndryValues(crse,c_start,fine,f_start,bnd_start,
                                    num_comp,ratio,phys_bc,max_order);

    // Now simply average data in BRs
    for (int mg_lev=1; mg_lev<num_mg_levels; ++mg_lev)
    {
        for (OrientationIter oitr; oitr; ++oitr)
        {
            Orientation face = oitr();
            FabSet& ffs = (mg_lev==1 ? bndry[face] : bnd_vals[mg_lev-1][face]);
            FabSet& cfs = bnd_vals[mg_lev][face];
            coarsen_data_in_fabsets(ffs,bnd_start,cfs,bnd_start,num_comp,face);
        }
    }
}

void
DDBndry::setBndryValues (const MultiFab& mf,
                         int             mf_start,
                         int             bnd_start,
                         int             num_comp,
                         const BCRec&    bc)
{
    // Set BR data at level 0, as before
    InterpBndryData::setBndryValues(mf,mf_start,bnd_start,num_comp,bc);
    
    // Now simply average data in BRs
    for (int mg_lev=1; mg_lev<num_mg_levels; ++mg_lev)
    {
        for (OrientationIter oitr; oitr; ++oitr)
        {
            Orientation face = oitr();
            FabSet& ffs = (mg_lev==1 ? bndry[face] : bnd_vals[mg_lev-1][face]);
            FabSet& cfs = bnd_vals[mg_lev][face];
            coarsen_data_in_fabsets(ffs,bnd_start,cfs,bnd_start,num_comp,face);
        }
    }
}

const FabSet&
DDBndry::bndryValues (const Orientation& face,int level) const
{
    BL_ASSERT(level>=0);
    if (level==0) 
    {
        return ViscBndry::bndryValues(face);
    }
    BL_ASSERT(level<bnd_vals.size());
    return bnd_vals[level][face];
}

const PArray<Mask>&
DDBndry::bndryMasks (const Orientation& face, int level) const
{
    BL_ASSERT(level>=0);
    if (level==0) 
    {
        return ViscBndry::bndryMasks(face);
    }
    BL_ASSERT(level<bnd_vals.size());
    return bnd_vals[level].bndryMasks(face);
}

#include "Utility.H"
#ifdef WIN32
static std::string sep = "\\";
#else
static std::string sep = "/";
#endif
void
DDBndry::Write(std::string& outfile) const
{
    if (ParallelDescriptor::IOProcessor()) {
        if (!BoxLib::UtilCreateDirectory(outfile, 0755))
            BoxLib::CreateDirectoryFailed(outfile);
    }
    ParallelDescriptor::Barrier();

    for (int lev=0; lev<bnd_vals.size(); ++lev)
    {
        int mindigits = 2;
        for (OrientationIter oitr; oitr; ++oitr) {
            char buf[32];
            sprintf(buf, "%0*d",  mindigits, (int)oitr());
            const FabSet& vals = bndryValues(oitr(),lev);
            vals.write(outfile+sep+"bndVals_"+std::string(buf));
        }
    }
}
