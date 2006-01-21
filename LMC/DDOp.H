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
        
    void setBoundaryData_Mass(const MultiFab&      fineT,
                              int                  fStartT,
                              const MultiFab&      fineY,
                              int                  fStartY,
                              const BndryRegister* cbrT,
                              int                  cStartT,
                              const BndryRegister* cbrY,
                              int                  cStartY,
                              const BCRec&         bcT,
                              const BCRec&         bcY);

    void setBoundaryData_Mole(const MultiFab&      fineT,
                              int                  fStartT,
                              const MultiFab&      fineX,
                              int                  fStartX,
                              const BndryRegister* cbrT,
                              int                  cStartT,
                              const BndryRegister* cbrX,
                              int                  cStartX,
                              const BCRec&         bcT,
                              const BCRec&         bcX);

    void setGrowCells(MultiFab& T,
                      int       compT,
                      MultiFab& X,
                      int       compX) const;

    void cellToEdge(MultiFab&       Te,
                    int             dCompTe,
                    MultiFab&       Ye,
                    int             dCompYe,
                    const MultiFab& Tc,
                    int             sCompTc,
                    const MultiFab& Yc,
                    int             sCompYc,
                    int             dir) const;

    void applyOp(MultiFab&         outH,
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
                 int               level = 0) const;

     void setRelax(MultiFab&         lambda,
                   int               dCompH,
                   const MultiFab&   inT,
                   int               sCompT,
                   const MultiFab&   inY,
                   int               sCompY,
                   int               level = 0) const;

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
    const DDBndry& XBndry() const;

private:

    void center_to_edge (const FArrayBox& cfab,
                         FArrayBox&       efab,
                         const Box&       ccBox,
                         int              sComp,
                         int              dComp,
                         int              nComp) const;

    const ChemDriver& ckdriver;
    BoxArray grids;
    DDBndry Tbd, Ybd, Xbd;
    IntVect cfRatio;
    MultiFab volume;
    MultiFab area[BL_SPACEDIM];
    DDOp* coarser;
    MultiFab Soln;
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
const DDBndry&
DDOp::XBndry() const
{
    return Xbd;
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
