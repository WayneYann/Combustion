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
        state[k].swapTimeLevels(dt);
    }

    MultiFab Uprime(grids,NUM_STATE,0);

    MultiFab& Unew = get_new_data(State_Type);
    MultiFab& Uold = get_old_data(State_Type);

    bool doFillpatch;

    if (level == 0)
    {
	doFillpatch = false;
	MultiFab::Copy(Unew, Uold, 0, 0, NUM_STATE, 0);
    }
    else
    {
	doFillpatch = true;
    }
    fill_boundary(Unew, time, doFillpatch);
    
    // do half-dt chemistry
    advance_chemistry(Unew, 0.5*dt);
    post_update(Unew);

    // Step 1 of RK2
    doFillpatch = false;
    dUdt(Unew, Uprime, time, doFillpatch);
    update_rk(Unew, Uold, 0.5*dt, Uprime); // Unew = Uold + 0.5*dt*Uprime
    post_update(Unew);

    // Step 2 of RK2
    doFillpatch = (level == 0) ? false : true;
    dUdt(Unew, Uprime, time+0.5*dt, doFillpatch);
    update_rk(Unew, Uold, dt, Uprime); // Unew = Uold + dt*Uprime
    post_update(Unew);

    // do another half-dt chemistry
    advance_chemistry(Unew, 0.5*dt);
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


void
RNS::fill_boundary(MultiFab& U, Real time, bool do_fillpatch)
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
}


// xxxxx need to add flux register for AMR
void
RNS::dUdt(MultiFab& U, MultiFab& Uprime, Real time, bool do_fillpatch)
{
    fill_boundary(U, time, do_fillpatch);
    
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
    for (MFIter mfi(U1); mfi.isValid(); ++mfi)
    {
	const int   i = mfi.index();
	const Box& bx = mfi.validbox();
	
	U1[i].copy(U2[i], bx);
	U1[i].saxpy(c, Uprime[i]);
    }
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

	BL_FORT_PROC_CALL(RNS_ENFORCE_CONSISTENT_Y, rns_enforce_consistent_y)
	    (lo, hi, BL_TO_FORTRAN(U[i]));

	BL_FORT_PROC_CALL(RNS_COMPUTE_TEMP, rns_compute_temp)
	    (lo, hi, BL_TO_FORTRAN(U[i]));	
    }
}


void
RNS::advance_chemistry(MultiFab& U, Real dt)
{
    if (chemSolve->isNull) return;

    for (MFIter mfi(U); mfi.isValid(); ++mfi)
    {
	const int   i = mfi.index();
	const Box& bx = mfi.validbox();
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();

	BL_FORT_PROC_CALL(RNS_ADVCHEM, rns_advchem)
	    (lo, hi, BL_TO_FORTRAN(U[i]), dt);
    }
}
