Real
HeatTransfer::rk_step_selection (const Real time)
{
/*
  estimate the time step for the Runge-Kutta implementation of the diffusion
  update.  this estimation is based on the diffusion coefficients of the 
  species into the mixture, independent of whether or not the "operator" is 
  based on the mixture averaged or the multicomponent formulation with Soret 
  and Dufour effects.

                                space_step ** 2
                  time_step  =  ---------------
                                D_mix * 2 * DIM

  the minimum of this formula is taken over all cells and species

  Real time     ! INPUT time, either prev_time or cur_time
  returned Real ! OUTPUT the time step
*/

    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::rk_step_selection()");

    // check that time is either prev_time or cur_time
    const TimeLevel whichTime = which_time(State_Type, time);
    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    // get some constants used as dimensions
    int ncomps = NUM_STATE;
    // const int nspecies = getChemSolve().numSpecies(); // this is globally defined

    // get some constants used in the calculation of temperature
    Real maximum_error = getChemSolve().getHtoTerrMAX ();
    int maximum_iterations = getChemSolve().getHtoTiterMAX ();

    // get locations of things in the state and convert to fortran indexing
    int index_of_firstY = Density + 1;
    int index_of_lastY = index_of_firstY + nspecies - 1;
    int index_of_rho = Density;
    int index_of_rhoH = RhoH;
    int index_of_T = Temp;
    int fort_index_of_firstY = 1 + index_of_firstY;
    int fort_index_of_lastY  = 1 + index_of_lastY;
    int fort_index_of_rho    = 1 + index_of_rho;
    int fort_index_of_rhoH   = 1 + index_of_rhoH;
    int fort_index_of_T      = 1 + index_of_T;

    // get the state at the desired time
    MultiFab& state 
	= (whichTime == AmrOldTime) 
	? get_old_data(State_Type) 
	: get_new_data(State_Type);

    // since it is possible that a processor has no boxes on this
    // level, the variable that holds the smallest time step found
    // here must be initialized with a large value so that the 
    // ParallelDescriptor::ReduceRealMin will function propoerly.
    Real smallest_time_step = 1.0e20;

    // loop over fabs in the state
    for (MFIter state_mfi(state); state_mfi.isValid(); ++state_mfi)
    {
	// get index of the present box
        const int idx = state_mfi.index();

/*
c     arguments are alphabetical, mostly:
c
c     lo, hi             ! INPUT limits of valid region of the box
c     dx                 ! INPUT physical dimensions of grid cells
c     index_of_firstY    ! INPUT index of rho Y for the first species in the state
c     index_of_lastY     ! INPUT index of rho Y for the last species in the state
c     index_of_rho       ! INPUT index of rho in the state
c     index_of_rhoH      ! INPUT index of rho H in the state
c     index_of_T         ! INPUT index of T in the state
c     maximum_error      ! INPUT maximum error in calculation of T
c     maximum_iterations ! INPUT maximum iterations in calculation of T
c     ncomps             ! INPUT total number of components in the state
c     nspecies           ! INPUT total number of species in the state
c     smallest_time_step ! OUTPUT smallest time step over all cells and species
c     state, DIMS(state) ! INPUT all variables in the state
*/

#define DATA_AND_LIMITS(foo) foo.dataPtr(),foo.loVect()[0],foo.loVect()[1],foo.hiVect()[0],foo.hiVect()[1]

	FORT_RK_STEP_SELECTION
	    (grids[idx].loVect(), grids[idx].hiVect(),
	     geom.CellSize(),
	     &fort_index_of_firstY,
	     &fort_index_of_lastY,
	     &fort_index_of_rho,
	     &fort_index_of_rhoH,
	     &fort_index_of_T,
	     &maximum_error,
	     &maximum_iterations,
	     &ncomps,
	     &nspecies,
	     &smallest_time_step,
	     DATA_AND_LIMITS(state[state_mfi])
	     );
    }
    // reconcile the time step over all the processors
    ParallelDescriptor::ReduceRealMin (smallest_time_step);

    // apply the safety factor
    Real time_step = rk_time_step_multiplier * smallest_time_step;

    // return the scaled time step
    return time_step;
}
