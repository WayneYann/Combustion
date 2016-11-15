
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

    FillPatchIterator fpi(*this,rhohmix,ngrow,time,State_Type,sComp,nComp);
    MultiFab& statemf = fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	FArrayBox tmp;

	for (MFIter mfi(rhohmix, true); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.tilebox();
	    FArrayBox& state  = statemf[mfi];
	    FArrayBox& rhmfab = rhohmix[mfi];

	    //
	    // Convert rho*Y to Y for this operation
	    //
	    tmp.resize(bx,1);
	    tmp.copy(state,sCompR,0,1);
	    tmp.invert(1.0);
	    
	    for (int k = 0; k < nspecies; k++) {
		state.mult(tmp,0,sCompY+k,1);
	    }
	    
	    getChemSolve().getHmixGivenTY(rhohmix[mfi],state,state,bx,
					  sCompT,sCompY,sCompH);
	    //
	    // Convert hmix to rho*hmix
	    //
	    rhohmix[mfi].mult(state,bx,sCompR,sCompH,1);
	}
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
