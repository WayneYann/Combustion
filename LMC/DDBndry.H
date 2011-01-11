//BL_COPYRIGHT_NOTICE

#ifndef _DDBndry_H_
#define _DDBndry_H_

//
// $Id: DDBndry.H,v 1.3 2011-01-11 03:18:46 marc Exp $
//

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
             int             mgLevel = 0);

    void define (const BoxArray& grids,
                 int             ncomp,
                 const Geometry& geom,
                 int             mgLevel);

    virtual void setBndryConds (const BCRec& phys_bc,
                                int          ratio);

    virtual void setBndryConds (const BCRec&       phys_bc,
                                /*const*/ IntVect& ratio,
                                int                comp = 0);

    int mgLevel() const {return mg_level;}

protected:
    int mg_level;
};

#endif /*_DDBndry_H_*/
