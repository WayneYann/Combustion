//BL_COPYRIGHT_NOTICE

//
// $Id: DDBndry.cpp,v 1.4 2007-03-27 17:28:05 lijewski Exp $
//
#include <winstd.H>
#include <cmath>

#include <LO_BCTYPES.H>
#include <DDBndry.H>

DDBndry::DDBndry (const BoxArray& grids,
                  int             ncomp,
                  const Geometry& geom,
                  int             mgLevel)
{
    define(grids,ncomp,geom,mgLevel);
}

void
DDBndry::define (const BoxArray& grids,
                 int             ncomp,
                 const Geometry& geom,
                 int             mgLevel)
{
    mg_level = mgLevel;
    ViscBndry::define(grids,ncomp,geom);
}

void
DDBndry::setBndryConds (const BCRec& phys_bc,
                        int          ratio)
{
    IntVect ratio_vect = ratio * IntVect::TheUnitVector();
    setBndryConds(phys_bc, ratio_vect);
}

void
DDBndry::setBndryConds (const BCRec&       bc,
                        /*const*/ IntVect& ratio,
                        int                comp)
{
    //
    //  NOTE: ALL BCLOC VALUES ARE DEFINED AS A LENGTH IN PHYSICAL
    //        DIMENSIONS *RELATIVE* TO THE FACE, NOT IN ABSOLUTE PHYSICAL SPACE
    //
    const Real* dx    = getGeom().CellSize();
    const Box& domain = getGeom().Domain();

    ViscBndry::setBndryConds(bc,ratio,comp);

    // Here, just overwrite the c-f boundary instructions
    for (OrientationIter fi; fi; ++fi)
    {
        const Orientation& face = fi();
        const int dir           = face.coordDir();
        const Real delta        = dx[dir]*ratio[dir]/pow(float(2),float(mg_level));
        const int p_bc          = (face.isLow() ? bc.lo(dir) : bc.hi(dir));

        for (int i = 0; i < boxes().size(); i++)
        {
            if (!(domain[fi()] == boxes()[i][fi()] && !getGeom().isPeriodic(dir)))
            {
                bcond[face][i][0] = LO_DIRICHLET;
                bcloc[face][i] = 0.5*delta;
            }
        }
    }
}

