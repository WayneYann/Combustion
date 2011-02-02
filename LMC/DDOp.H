#ifndef _DDOP_H_
#define _DDOP_H_

#include "DDBndry.H"
#include "ChemDriver.H"

class DDOp
{
public:

    DDOp ();

    DDOp (const BoxArray&   grids,
          const Box&        box,
          const IntVect&    amrRatio);

    ~DDOp ();

    void define (const BoxArray& grids,
                 const Box&      box,
                 const IntVect&  amrRatio);
        
    void setBoundaryData(const MultiFab&      fineT,
                         int                  fStartT,
                         const MultiFab&      fineX,
                         int                  fStartX,
                         const BndryRegister* cbrT,
                         int                  cStartT,
                         const BndryRegister* cbrX,
                         int                  cStartX,
                         const BCRec&         bcT,
                         const BCRec&         bcX,
                         int                  max_order = InterpBndryData::IBD_max_order_DEF);
    
    enum DD_ApForTorRH {DD_Temp=0, DD_RhoH, DD_NumAppTypes};
    enum DD_Model {DD_Model_Full=0, DD_Model_MixAvg, DD_Model_NumModels};

    static void set_transport_model(const DD_Model& model);
    static void set_chem_driver(ChemDriver& chem);
    static void set_mgLevelsMAX(int mgLevelsMAX);

    void applyOp(MultiFab&         outYH,
                 const MultiFab&   inYT,
                 PArray<MultiFab>& fluxYH,
                 DD_ApForTorRH     whichApp,
                 bool              updateCoefs = true,
                 int               level = 0,
                 bool              getAlpha = false,
                 MultiFab*         alpha = 0);

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

    static int numLevels();
    const Box& domain() const {return geom.Domain();}
    const BoxArray& boxArray(int level) const;
    void Write(std::string& outfile) const;
    void WriteSub(std::string& outfile, int level) const;

protected:

    void WriteSub(std::string& outfile) const;

    DDOp (const BoxArray& grids,
          const Box&      box,
          const IntVect&  amrRatio,
          int             mgLevel,
          DDOp*           parent);

    void define (const BoxArray& grids,
                 const Box&      box,
                 const IntVect&  amrRatio,
                 int             mgLevel,
                 DDOp*           parent);
        
    bool thisIsParent() const;

    void applyOp_DoIt(MultiFab&         outYH,
                      const MultiFab&   inYT,
                      PArray<MultiFab>& fluxYH,
                      DD_ApForTorRH     whichApp,
                      bool              updateCoefs,
                      bool              getAlpha,
                      MultiFab*         alpha);

    void setGrowCells(MultiFab& YT,
                      int       compT,
                      int       compY);

    void ensure_valid_transport_is_set() const;

    void center_to_edge (const FArrayBox& cfab,
                         FArrayBox&       efab,
                         const Box&       ccBox,
                         int              sComp,
                         int              dComp,
                         int              nComp) const;

    static ChemDriver* chem;
    BoxArray grids;
    DDBndry *Tbd, *Ybd; // only reeally defined at mg_level=0
    BndryRegister stencilWeight;
    IntVect amr_ratio;
    MultiFab volInv;
    PArray<MultiFab> area;
    PArray<MultiFab> flux;
    MultiFab coefs;
    DDOp* mg_parent;
    PArray<DDOp> ddOps;
    Array<Real> dx;
    Geometry geom;
    int mg_level;
    static DD_Model transport_model;
    static int maxorder;
    static int mgLevelsMAX;
    static int mgLevels;
    static bool first_define;
};

inline
void
DDOp::set_transport_model(const DD_Model& model)
{
    transport_model = model;
}

inline
void
DDOp::set_chem_driver(ChemDriver& chemdriver)
{
    chem = &chemdriver;
}

inline
void
DDOp::set_mgLevelsMAX(int mgLevelsMAX_)
{
    if (mgLevelsMAX > mgLevels) {
        BoxLib::Abort("DDOp::set_mgLevelsMAX: Cannot set mgLevelsMAX less than number presently used");
    }
    mgLevelsMAX = mgLevelsMAX_;
}

inline
int
DDOp::numLevels()
{
    return mgLevels;
}

#endif
