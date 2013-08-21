
#include "LevelBld.H"
#include "RNS.H"
#include <cmath>

#ifdef USE_SDCLIB
#include "SDCLevelBld.H"
#endif

class RNSBld
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



#ifdef USE_SDCLIB
class RNSSDCBld
    :
    public SDCLevelBld
{
  sdc_sweeper_t* build (SDCAmr& papa,
                        int     lev,
                        Real    time)
  {
    sdc_nodes_t nodes;
    int nnodes0 = 3;
    int trat    = 2;
    int nnodes  = 1 + (nnodes0 - 1) * ((int) pow(trat, lev));
    sdc_nodes_build(&nodes, nnodes, SDC_GAUSS_LOBATTO);
    sdc_imex_t* imex = sdc_imex_create(&nodes, NULL, NULL, NULL);
    return (sdc_sweeper_t*) imex;
  }
};
#endif

