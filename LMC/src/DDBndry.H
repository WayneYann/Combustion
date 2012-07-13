
#ifndef _DDBndry_H_
#define _DDBndry_H_

#include <ViscBndry.H>
#include <BC_TYPES.H>

class DDBndry
    :
    public ViscBndry
{
public:
    DDBndry()
        : ViscBndry()  {}

    DDBndry (const BoxArray& grids,
             int             ncomp,
             const Geometry& geom,
             int             numMGlevels);

    virtual void define  (const BoxArray& grids,
                          int             ncomp,
                          const Geometry& geom,
                          int             numMGlevels);

    virtual void setBndryValues (const MultiFab& mf,
                                 int             mf_start,
                                 int             bnd_start,
                                 int             num_comp,
                                 const BCRec&    bc);

    virtual void setBndryValues (::BndryRegister&  crse,
                                 int             c_start,
                                 const MultiFab& fine,
                                 int             f_start,
                                 int             bnd_start,
                                 int             num_comp,
                                 int             ratio, 
                                 const BCRec&    phys_bc,
                                 int             max_order = InterpBndryData::IBD_max_order_DEF);    

    virtual void setBndryValues (::BndryRegister&  crse,
                                 int             c_start,
                                 const MultiFab& fine,
                                 int             f_start,
                                 int             bnd_start,
                                 int             num_comp,
                                 IntVect&        ratio, 
                                 const BCRec&    phys_bc,
                                 int             max_order = InterpBndryData::IBD_max_order_DEF);

    const Geometry& getGeom(int level) const {return geoms[level];}
    const FabSet& bndryValues (const Orientation& face,int level) const;
    //
    // A useful typedef.
    //
    typedef BndryData::MaskTuple MaskTuple;

    const MaskTuple& bndryMasks (int igrid, int level) const;

    void Write(std::string& outfile) const;

protected:

    // Hide these, since there is now a level associated with each
    const Geometry& getGeom() const;
    const FabSet& bndryValues (const Orientation& face) const;

    int num_mg_levels;
    PArray<BndryData> bnd_vals;
    PArray<Geometry> geoms;
};

#endif /*_DDBndry_H_*/
