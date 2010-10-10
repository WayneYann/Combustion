
#include "LevelBld.H"
#include "CNSReact.H"

class CNSReactBld
    :
    public LevelBld
{
    virtual void variableSetUp ();
    virtual void variableCleanUp ();
    virtual AmrLevel *operator() ();
    virtual AmrLevel *operator() (Amr&            papa,
                                  int             lev,
                                  const Geometry& level_geom,
                                  const BoxArray& ba,
                                  Real            time);
};

CNSReactBld CNSReact_bld;

LevelBld*
getLevelBld ()
{
    return &CNSReact_bld;
}

void
CNSReactBld::variableSetUp ()
{
    CNSReact::variableSetUp();
}

void
CNSReactBld::variableCleanUp ()
{
    CNSReact::variableCleanUp();
}

AmrLevel*
CNSReactBld::operator() ()
{
    return new CNSReact;
}

AmrLevel*
CNSReactBld::operator() (Amr&            papa,
                       int             lev,
                       const Geometry& level_geom,
                       const BoxArray& ba,
                       Real            time)
{
    return new CNSReact(papa, lev, level_geom, ba, time);
}
