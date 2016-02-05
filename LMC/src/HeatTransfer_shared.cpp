
#include <HeatTransfer.H>

#ifdef PARTICLES
namespace
{
    bool do_curvature_sample = false;
}

void
HeatTransfer::read_particle_params ()
{
    ParmParse ppht("ht");
    ppht.query("do_curvature_sample", do_curvature_sample);
}

int
HeatTransfer::timestamp_num_extras ()
{
    return do_curvature_sample ? 1 : 0;
}

void
HeatTransfer::timestamp_add_extras (int lev,
				    Real time,
				    MultiFab& mf)
{
    if (do_curvature_sample)
    {
	AmrLevel& amr_level = parent->getLevel(lev);
	int cComp = mf.nComp()-1;

	amr_level.derive("mean_progress_curvature", time, mf, cComp);

	mf.setBndry(0,cComp,1);
    }
}

#endif /*PARTICLES*/

