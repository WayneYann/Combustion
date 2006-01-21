//BL_COPYRIGHT_NOTICE

//
// $Id: DDBndry.cpp,v 1.1 2006-01-21 02:12:38 marc Exp $
//

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
DDBndry::setBndryConds (const BCRec&       bc,
                        /*const*/ IntVect& ratio)
{
    //
    //  NOTE: ALL BCLOC VALUES ARE DEFINED AS A LENGTH IN PHYSICAL
    //        DIMENSIONS *RELATIVE* TO THE FACE, NOT IN ABSOLUTE PHYSICAL SPACE
    //
    const Real* dx    = getGeom().CellSize();
    const Box& domain = getGeom().Domain();

    ViscBndry::setBndryConds(bc,ratio);

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
                //setBoundCond(face,i,LO_DIRICHLET);
                //setBoundLoc(face,i,0.5*delta);
            }
        }
    }
}

