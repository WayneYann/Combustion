#include <winstd.H>

#include "RNS.H"
#include "RNS_F.H"

using std::string;

Real
RNS::advance (Real time,
	      Real dt,
	      int  iteration,
	      int  ncycle)
{
    if (use_sdc)
    {
	advance_sdc(time,dt,iteration,ncycle);
    }
    else
    {
	advance_rk(time,dt,iteration,ncycle);
    }
    return dt;
}

void
RNS::advance_sdc (Real time,
		  Real dt,
		  int  iteration,
		  int  ncycle)
{
    
}

void
RNS::advance_rk (Real time,
		 Real dt,
		 int  iteration,
		 int  ncycle)
{
    BL_ASSERT( state[State_Type].sizeMidData() == 1 );

    MultiFab Uprime(grids,NUM_STATE,0);

    MultiFab& Unew = get_new_data(State_Type);
    MultiFab& Uold = get_old_data(State_Type);
    MultiFab& Umid = get_mid_data(State_Type, 0);

    dUdt(Uold, Uprime, time);
    update_rk(Umid, Uold, 0.5*dt, Uprime); // Umid = Uold + 0.5*dt*Uprime
    // reset species mass fraction here

    dUdt(Umid, Uprime, time+0.5*dt);
    update_rk(Unew, Uold, dt, Uprime); // Unew = Uold + dt*Uprime
    // reset species mass fraction here
}

// xxxxx need to add flux register for AMR
void
RNS::dUdt(MultiFab& U, MultiFab& Uprime, Real time)
{
    for (FillPatchIterator fpi(*this, U, NUM_GROW, time, State_Type, 0, NUM_STATE); 
	 fpi.isValid(); ++fpi) 
    {
	U[fpi].copy(fpi());
    }

    // for testing
    Uprime.setVal(0.0);
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
