#include <winstd.H>

#include "RNS.H"
#include "RNS_F.H"

using std::string;

// xxxxx need to add flux register for AMR
Real
RNS::advance (Real time,
	      Real dt,
	      int  iteration,
	      int  ncycle)
{
    for (int k = 0; k < NUM_STATE_TYPE; k++) 
    {
        state[k].allocOldData();
        state[k].allocMidData();
        state[k].swapTimeLevels(dt);
    }

    BL_ASSERT( state[State_Type].sizeMidData() == 1 );

    MultiFab Uprime(grids,NUM_STATE,0);

    MultiFab& Unew = get_new_data(State_Type);
    MultiFab& Uold = get_old_data(State_Type);
    MultiFab& Umid = get_mid_data(State_Type, 0);

    dUdt(Uold, Uprime, time);
    update_rk(Umid, Uold, 0.5*dt, Uprime); // Umid = Uold + 0.5*dt*Uprime
    post_update(Umid);

    bool doFillpatch = (level == 0) ? false : true; // no valid mid data on level 0
    dUdt(Umid, Uprime, time+0.5*dt, doFillpatch);
    update_rk(Unew, Uold, dt, Uprime); // Unew = Uold + dt*Uprime
    post_update(Unew);

    if (Unew.contains_nan(0,NUM_STATE,0))
    {
	for (int i=0; i<NUM_STATE; i++)
	{
	    if (Unew.contains_nan(i, 1, 0))
	    {
		std::cout << "RNS::advance: Testing component i for NaNs: " << i << std::endl;
                BoxLib::Abort("RNS::advance: Has NaNs in this component.");
	    }
	}
    }
    
    return dt;
}

// xxxxx need to add flux register for AMR
void
RNS::dUdt(MultiFab& U, MultiFab& Uprime, Real time, bool do_fillpatch)
{

    if (do_fillpatch)
    {
	for (FillPatchIterator fpi(*this, U, NUM_GROW, time, State_Type, 0, NUM_STATE); 
	     fpi.isValid(); ++fpi) 
	{
	    U[fpi].copy(fpi());
	}
    }
    else
    {
	U.FillBoundary();
	geom.FillPeriodicBoundary(U, true);
	for (MFIter mfi(U); mfi.isValid(); ++mfi)
	{
	    setPhysBoundaryValues(U[mfi], 
				  State_Type,
				  time,
				  0,
				  0,
				  NUM_STATE);
	}
    }

    
    const Real *dx = geom.CellSize();

    for (MFIter mfi(Uprime); mfi.isValid(); ++mfi)
    {
	int i = mfi.index();
	const Box& bx = mfi.validbox();

	BL_FORT_PROC_CALL(RNS_DUDT,rns_dudt)
	    (bx.loVect(), bx.hiVect(),
	     BL_TO_FORTRAN(U[i]),
	     BL_TO_FORTRAN(Uprime[i]),
	     dx);
    }
}


// Compute U1 = U2 + c Uprime.
void 
RNS::update_rk(MultiFab& U1, const MultiFab& U2, Real c, const MultiFab& Uprime)
{
    // not meant for performance
    MultiFab::Copy(U1, Uprime, 0, 0, NUM_STATE, 0);
    U1.mult(c);
    U1.plus(U2, 0, NUM_STATE, 0);
}


void
RNS::post_update(MultiFab& U)
{
    for (MFIter mfi(U); mfi.isValid(); ++mfi)
    {
	const int   i = mfi.index();
	const Box& bx = mfi.validbox();
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();

//	    BL_FORT_PROC_CALL(RNS_ENFORCE_CONSISTENT_Y, rns_enforce_consistent_y)
//		(lo, hi, BL_TO_FORTRAN(U[i]));

	BL_FORT_PROC_CALL(RNS_COMPUTE_TEMP, rns_compute_temp)
	    (lo, hi, BL_TO_FORTRAN(U[i]));	
    }
}
