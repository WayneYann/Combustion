
#include "LevelBld.H"
#include "RNS.H"
#include <cmath>

class RNSBld
    :
    public LevelBld
{
    virtual void variableSetUp () override;
    virtual void variableCleanUp () override;
    virtual AmrLevel *operator() () override;
    virtual AmrLevel *operator() (Amr&            papa,
                                  int             lev,
                                  const Geometry& level_geom,
                                  const BoxArray& ba,
                                  Real            time) override;
};

RNSBld RNS_bld;

LevelBld*
getLevelBld ()
{
    return &RNS_bld;
}

void
RNSBld::variableSetUp ()
{
    RNS::variableSetUp();
}

void
RNSBld::variableCleanUp ()
{
    RNS::variableCleanUp();
}

AmrLevel*
RNSBld::operator() ()
{
    return new RNS;
}

AmrLevel*
RNSBld::operator() (Amr&            papa,
		    int             lev,
		    const Geometry& level_geom,
		    const BoxArray& ba,
		    Real            time)
{
    return new RNS(papa, lev, level_geom, ba, time);
}
