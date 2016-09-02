
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

void
HeatTransfer::compute_rhohmix (Real      time,
                               MultiFab& rhohmix)
{
    const Real strt_time = ParallelDescriptor::second();

    const int ngrow  = 0; // We only do this on the valid region
    const int sComp  = std::min(std::min((int)Density,(int)Temp),first_spec);
    const int eComp  = std::max(std::max((int)Density,(int)Temp),first_spec+nspecies-1);
    const int nComp  = eComp - sComp + 1;
    const int sCompR = Density - sComp;
    const int sCompT = Temp - sComp;
    const int sCompY = first_spec - sComp;
    const int sCompH = 0;

    FArrayBox tmp;

    for (FillPatchIterator fpi(*this,rhohmix,ngrow,time,State_Type,sComp,nComp);
         fpi.isValid();
         ++fpi)
    {
        //
        // Convert rho*Y to Y for this operation
        //
	const int  iGrid    = fpi.index();
	const Box& validBox = grids[iGrid];
	FArrayBox& state    = fpi();

        tmp.resize(validBox,1);
        tmp.copy(state,sCompR,0,1);
        tmp.invert(1);

        for (int k = 0; k < nspecies; k++)
            state.mult(tmp,0,sCompY+k,1);
	
	getChemSolve().getHmixGivenTY(rhohmix[fpi],state,state,validBox,
				      sCompT,sCompY,sCompH);
        //
        // Convert hmix to rho*hmix
        //
        rhohmix[fpi].mult(state,validBox,sCompR,sCompH,1);
    }

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::compute_rhohmix(): lev: " << level << ", time: " << run_time << '\n';
    }
}
