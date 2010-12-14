#ifndef _DDOP_H_
#define _DDOP_H_

#include "DDBndry.H"
#include "ChemDriver.H"

class DDOp
{
public:

    DDOp (const ChemDriver& ckd);

    DDOp (const BoxArray&   grids,
          const Box&        box,
          const ChemDriver& ckd,
          const IntVect&    ratio,
          int               mgLevel = 0);

    ~DDOp ();

    void define (const BoxArray& grids,
                 const Box&      box,
                 const IntVect&  ratio,
                 int             mgLevel = 0);
        
    void setBoundaryData(const MultiFab&      fineT,
                         int                  fStartT,
                         const MultiFab&      fineX,
                         int                  fStartX,
                         const BndryRegister* cbrT,
                         int                  cStartT,
                         const BndryRegister* cbrX,
                         int                  cStartX,
                         const BCRec&         bcT,
                         const BCRec&         bcX);

    void setGrowCells(MultiFab& YT,
                      int       compT,
                      int       compY) const;

    enum DD_ApForTorRH {DD_Temp=0, DD_RhoH, DD_NumAppTypes};
    enum DD_Model {DD_Model_Full=0, DD_Model_MixAvg, DD_Model_NumModels};

    void applyOp(MultiFab&         outYH,
                 const MultiFab&   inYT,
                 PArray<MultiFab>& fluxYH,
                 DD_ApForTorRH     whichApp,
                 int               level = 0,
                 bool              getAlpha = false,
                 MultiFab*         alpha = 0,
                 bool              updateCoefs = true);

    void setCoefficients(const MFIter&    mfi,
                         const FArrayBox& inYT,
                         const FArrayBox& cpi);
    
    static void average(MultiFab&       mfC,
                        int             dCompC,
                        const MultiFab& mfF,
                        int             sCompF,
                        int             nComp);

    static void interpolate(MultiFab&       mfF,
                            int             dCompF,
                            const MultiFab& mfC,
                            int             sCompC,
                            int             nComp);

    bool coarser_exists(int level) const;
    const Box& domain() const;
    const BoxArray& boxArray() const;
    const DDBndry& TBndry() const;
    const DDBndry& YBndry() const;

private:

    void center_to_edge (const FArrayBox& cfab,
                         FArrayBox&       efab,
                         const Box&       ccBox,
                         int              sComp,
                         int              dComp,
                         int              nComp) const;

    const ChemDriver& ckdriver;
    BoxArray grids;
    DDBndry Tbd, Ybd;
    IntVect cfRatio;
    MultiFab volume;
    MultiFab area[BL_SPACEDIM];
    MultiFab coefs;
    DDOp* coarser;
    static DD_Model transport_model;
};

inline
const DDBndry&
DDOp::TBndry() const
{
    return Tbd;
}

inline
const DDBndry&
DDOp::YBndry() const
{
    return Ybd;
}

inline
const Box&
DDOp::domain() const
{
    BL_ASSERT(Tbd.getGeom().Domain()==Ybd.getGeom().Domain());
    return Tbd.getGeom().Domain();
}

inline
const BoxArray&
DDOp::boxArray() const
{
    return grids;
}

#endif
