
#include <winstd.H>

#include <new>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cmath>

#ifndef WIN32
#include <unistd.h>
#endif

#include <CArena.H>
#include <REAL.H>
#include <Utility.H>
#include <IntVect.H>
#include <Box.H>
#include <Amr.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <AmrLevel.H>

#ifdef USE_SDCLIB
#include <SDCAmr.H>
#endif

int main (int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    BL_PROFILE_VAR("main()", pmain);

    BL_PROFILE_VAR("main-init", pmain_init);

    Real dRunTime1 = ParallelDescriptor::second();

    std::cout << std::setprecision(10);

    int  max_step;
    Real strt_time;
    Real stop_time;
    Real walltime_limit;
    int walltime_int;
    ParmParse pp; 

    max_step  = -1;
    strt_time =  0.0;
    stop_time = -1.0;
    walltime_limit = -1.0;
    walltime_int = 10;
    
    pp.query("max_step",max_step);
    pp.query("strt_time",strt_time);
    pp.query("stop_time",stop_time);
    pp.query("walltime_limit",walltime_limit);
    pp.query("walltime_int",walltime_int);
    
    if (strt_time < 0.0)
	BoxLib::Abort("MUST SPECIFY a non-negative strt_time"); 

    if (max_step < 0 && stop_time < 0.0)
	BoxLib::Abort("Exiting because neither max_step nor stop_time is non-negative.");

#ifndef USE_SDCLIB
    Amr* amrptr = new Amr;
#else
    Amr* amrptr = new SDCAmr();
#endif

    amrptr->init(strt_time,stop_time);

    // If we set the regrid_on_restart flag and if we are *not* going to take
    //    a time step then we want to go ahead and regrid here.
    if ( amrptr->RegridOnRestart() && 
	 ( (amrptr->levelSteps(0) >= max_step) ||
	   (amrptr->cumTime() >= stop_time) ) )
    {
	//
	// Regrid only!
	//
	amrptr->RegridOnly(amrptr->cumTime());
    }

    BL_PROFILE_VAR_STOP(pmain_init);

    int istep = 0;
    int walltime_limit_reached = 0;
    while ( amrptr->okToContinue()                            &&
	    (amrptr->levelSteps(0) < max_step || max_step < 0) &&
	    (amrptr->cumTime() < stop_time || stop_time < 0.0) &&
	    !walltime_limit_reached )
    {
	//
	// Do a timestep.
	//
	amrptr->coarseTimeStep(stop_time);

	istep++;
	if (walltime_limit > 0 && istep % walltime_int == 0)
	{
	    Real walltime_used = (ParallelDescriptor::second() - dRunTime1) / 3600.0;
	    walltime_limit_reached = walltime_used > walltime_limit;
	    ParallelDescriptor::Bcast(&walltime_limit_reached, 1);
	}
    }
  
    // Write final checkpoint and plotfile
    if (amrptr->stepOfLastCheckPoint() < amrptr->levelSteps(0)) 
    {
	amrptr->checkPoint();
    }
  
    if (amrptr->stepOfLastPlotFile() < amrptr->levelSteps(0)) 
    {
	amrptr->writePlotFile();
    }
    
    delete amrptr;
    
    //
    // This MUST follow the above delete as ~Amr() may dump files to disk.
    //
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    
    Real dRunTime2 = ParallelDescriptor::second() - dRunTime1;
    
    ParallelDescriptor::ReduceRealMax(dRunTime2,IOProc);
    
    if (ParallelDescriptor::IOProcessor()) 
	std::cout << "Run time = " << dRunTime2 << std::endl;
  
    if (CArena* arena = dynamic_cast<CArena*>(BoxLib::The_Arena())) 
    {
	//
	// A barrier to make sure our output follows that of RunStats.
	//
	ParallelDescriptor::Barrier();
	//
	// We're using a CArena -- output some FAB memory stats.
	//
	// This'll output total # of bytes of heap space in the Arena.
	//
	// It's actually the high water mark of heap space required by FABs.
	//
	char buf[256];
	
	sprintf(buf,
		"CPU(%d): Heap Space (bytes) used by Coalescing FAB Arena: %ld",
		ParallelDescriptor::MyProc(),
		arena->heap_space_used());
	
	std::cout << buf << std::endl;
    }

    BL_PROFILE_VAR_STOP(pmain); 
    
    BoxLib::Finalize();
    return 0;
}
