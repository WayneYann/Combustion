#include <iostream>
#include <iomanip>
#include <cmath>

#include <RNS.H>
#include <RNS_F.H>

void 
RNS::sum_conserved_variables()
{
    int finest_level = parent->finestLevel();
    Real time        = state[State_Type].curTime();

    Array<Real> s(Eden+1, 0.0);

    for (int lev = 0; lev <= finest_level; lev++)
    {
	RNS& rns = getLevel(lev);

	rns.volWgtSumCons(s);
    }
    
    ParallelDescriptor::ReduceRealSum(&s[0], s.size(), ParallelDescriptor::IOProcessorNumber());    

    if (ParallelDescriptor::IOProcessor())
    {
	int oldprec = std::cout.precision(16);
	std::cout << "At t = " << time << std::endl;
	std::cout << "   Total Mass       is " << s[Density] << std::endl;
	std::cout << "   Total x-Momentum is " << s[Xmom] << std::endl;
#if (BL_SPACEDIM >= 2)
	std::cout << "   Total y-Momentum is " << s[Ymom] << std::endl;
#endif
#if (BL_SPACEDIM >= 3)
	std::cout << "   Total z-Momentum is " << s[Zmom] << std::endl;
#endif
	std::cout << "   Total Energy     is " << s[Eden] << "\n";
	std::cout.precision(oldprec);
    }
}


void
RNS::volWgtSumCons(Array<Real>& s)
{
    const MultiFab& Unew = get_new_data(State_Type);

    BoxArray baf;
    if (level < parent->finestLevel()) {
	baf = parent->boxArray(level+1);
	baf.coarsen(parent->refRatio(level));
    }

    for (MFIter mfi(Unew); mfi.isValid(); ++mfi)
    {
	int i = mfi.index();
        const Box& bx = mfi.validbox();

	FArrayBox msk(bx,1);
	msk.setVal(1.0); 
	// Now mask off finer level
	if (level < parent->finestLevel()) {
	    std::vector< std::pair<int,Box> > isects = baf.intersections(bx);
	    for (int ii = 0; ii < isects.size(); ii++) {
		msk.setVal(0.0, isects[ii].second, 0);
	    }
	}

	BL_FORT_PROC_CALL(RNS_SUM_CONS, rns_sum_cons)
	    (BL_TO_FORTRAN(Unew[i]),
	     BL_TO_FORTRAN(msk),
	     BL_TO_FORTRAN(volume[i]),
	     &s[0]);
    }
}
