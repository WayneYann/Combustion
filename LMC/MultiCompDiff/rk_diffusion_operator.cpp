void
HeatTransfer::rk_diffusion_operator (const Real time,
				     const Real dt,
				     MultiFab **& flux_for_H,
				     MultiFab **& flux_for_Y,
				     MultiFab *& update_for_H,
				     MultiFab *& update_for_Y)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::rk_diffusion_operator()");

    // check that time is either prev_time or cur_time
    const TimeLevel whichTime = which_time(State_Type, time);
    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    if (ParallelDescriptor::IOProcessor())
    {
	if (whichTime == AmrOldTime)
	    std::cout << "JFG: time is AmrOldTime\n" << std::flush;
	else
	    std::cout << "JFG: time is AmrNewTime\n" << std::flush;
    }

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

    // allocate OUTPUT multifabs for fluxes and updates
    int ngrow = 0;
    diffusion->allocFluxBoxesLevel(flux_for_H,ngrow,1);
    diffusion->allocFluxBoxesLevel(flux_for_Y,ngrow,nspecies);
    update_for_H = new MultiFab(grids,1,ngrow);
    update_for_Y = new MultiFab(grids,nspecies,ngrow);

    // loop over fabs in the state at the specified time
    MultiFab dummy (grids,1,0,Fab_noallocate);
    MFIter update_for_H_mfi(*update_for_H);
    MFIter update_for_Y_mfi(*update_for_Y);
    MFIter xflux_for_H_mfi(*flux_for_H[0]);
    MFIter yflux_for_H_mfi(*flux_for_H[1]);
    MFIter xflux_for_Y_mfi(*flux_for_Y[0]);
    MFIter yflux_for_Y_mfi(*flux_for_Y[1]);
    ngrow = 1;
    for (FillPatchIterator state_fpi (*this, dummy, ngrow, time, State_Type, 0, ncomps);
         state_fpi.isValid();
         ++state_fpi, 
	     ++update_for_H_mfi, 
	     ++update_for_Y_mfi,
	     ++xflux_for_H_mfi,
	     ++yflux_for_H_mfi,
	     ++xflux_for_Y_mfi,
	     ++yflux_for_Y_mfi)
    {
        BL_ASSERT (
	    update_for_H_mfi.isValid() &&
	    update_for_Y_mfi.isValid() &&
	    xflux_for_H_mfi.isValid() &&
	    yflux_for_H_mfi.isValid() &&
	    xflux_for_Y_mfi.isValid() &&
	    yflux_for_Y_mfi.isValid()
	    );

	// get index of the present box
        const int idx = state_fpi.index();

	// get boundary condition array for all components
	Array<int> bc = getBCArray (State_Type, idx, 0, ncomps);

        // print some stuff
	const int* lo_vect = grids[idx].loVect();
	const int* hi_vect = grids[idx].hiVect();
	if (ParallelDescriptor::IOProcessor())
	    std::cout << "JFG: in rk_diffusion\n" 
		      << "state_fpi.index() = " << state_fpi.index() << "\n"
		      << "lo_vect[0] = " << lo_vect[0] << "\n"
		      << std::flush;

/*
c     arguments are alphabetical, mostly:
c
c     domain_lo, domain_hi,             ! INPUT limits of valid region of the domain
c     lo, hi,                           ! INPUT limits of valid region of the box
c     areax, DIMS(areax),               ! INPUT areas of the faces perpendicular to x axis
c     areay, DIMS(areay),               ! INPUT areas of the faces perpendicular to y axis
c     bc,                               ! INPUT boundary condition array for all comps
c     dt,                               ! INPUT timestep
c     dx,                               ! INPUT physical dimensions of grid cells
c     index_of_firstY,                  ! INPUT index of rho Y for the first species in the state
c     index_of_lastY,                   ! INPUT index of rho Y for the last species in the state
c     index_of_rho,                     ! INPUT index of rho in the state
c     index_of_rhoH,                    ! INPUT index of rho H in the state
c     index_of_T,                       ! INPUT index of T in the state
c     maximum_error,                    ! INPUT maximum error in calculation of T
c     maximum_iterations,               ! INPUT maximum iterations in calculation of T
c     ncomps,                           ! INPUT total number of components in the state
c     nspecies,                         ! INPUT total number of species in the state
c     state, DIMS(state),               ! INPUT all variables in the state
c     update_for_H, DIMS(update_for_H), ! OUTPUT divergences of the fluxes
c     update_for_Y, DIMS(update_for_Y), ! OUTPUT divergences of the fluxes
c     volume, DIMS(volume),             ! INPUT volumes of the cells
c     xflux_for_H, DIMS(xflux_for_H),   ! OUTPUT x fluxes for enthalpy
c     xflux_for_Y, DIMS(xflux_for_Y),   ! OUTPUT x fluxes for species
c     yflux_for_H, DIMS(yflux_for_H),   ! OUTPUT y fluxes for enthalpy
c     yflux_for_Y, DIMS(yflux_for_Y),   ! OUTPUT y fluxes for species
*/

#define DATA_AND_LIMITS(foo) foo.dataPtr(),foo.loVect()[0],foo.loVect()[1],foo.hiVect()[0],foo.hiVect()[1]

	FORT_RK_DIFFUSION (geom.Domain().loVect(), geom.Domain().hiVect(), 
			   grids[idx].loVect(), grids[idx].hiVect(),
			   DATA_AND_LIMITS(area[0][idx]),
			   DATA_AND_LIMITS(area[1][idx]),
			   bc.dataPtr(),
			   &dt,
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
			   DATA_AND_LIMITS(state_fpi()),
			   DATA_AND_LIMITS((*update_for_H)[update_for_H_mfi]),
			   DATA_AND_LIMITS((*update_for_Y)[update_for_Y_mfi]),
			   DATA_AND_LIMITS(volume[idx]),
			   DATA_AND_LIMITS((*flux_for_H[0])[xflux_for_H_mfi]),
			   DATA_AND_LIMITS((*flux_for_Y[0])[xflux_for_Y_mfi]),
			   DATA_AND_LIMITS((*flux_for_H[1])[yflux_for_H_mfi]),
			   DATA_AND_LIMITS((*flux_for_Y[1])[yflux_for_Y_mfi]));
    }
}
