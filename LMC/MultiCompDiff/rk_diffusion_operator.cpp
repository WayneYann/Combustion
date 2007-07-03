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

    // get some constants used as dimensions
    int ncomps = NUM_STATE;
    // const int nspecies = getChemSolve().numSpecies(); // this is globally defined

    // get locations of things in the state
    int index_of_firstY = Density + 1;
    int index_of_lastY = index_of_firstY + nspecies - 1;
    int index_of_rho = Density;
    int index_of_rhoH = RhoH;
    int index_of_T = Temp;

    // allocate INPUT multifab for state
    int ngrow = 1;
    MultiFab state (grids, ncomps, ngrow);

    // allocate OUTPUT multifabs for fluxes and updates
    ngrow = 0;
    diffusion->allocFluxBoxesLevel(flux_for_H,ngrow,1);
    diffusion->allocFluxBoxesLevel(flux_for_Y,ngrow,nspecies);
    update_for_H = new MultiFab(grids,1,ngrow);
    update_for_Y = new MultiFab(grids,nspecies,ngrow);

    // loop over boxes
    MFIter update_for_H_mfi(*update_for_H);
    MFIter update_for_Y_mfi(*update_for_Y);
    MFIter xflux_for_H_mfi(*flux_for_H[0]);
    MFIter yflux_for_H_mfi(*flux_for_H[1]);
    MFIter xflux_for_Y_mfi(*flux_for_Y[0]);
    MFIter yflux_for_Y_mfi(*flux_for_Y[1]);
    for (FillPatchIterator state_fpi (*this, state, ngrow, time, State_Type, 0, ncomps);
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
c     index_of_T                        ! INPUT index of T in the state
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
			   &index_of_firstY,
			   &index_of_lastY,
			   &index_of_rho,
			   &index_of_rhoH,
			   &index_of_T,
			   &ncomps,
			   &nspecies,
			   DATA_AND_LIMITS(state[idx]),
			   DATA_AND_LIMITS((*update_for_H)[update_for_H_mfi]),
			   DATA_AND_LIMITS((*update_for_Y)[update_for_Y_mfi]),
			   DATA_AND_LIMITS(volume[idx]),
			   DATA_AND_LIMITS((*flux_for_H[0])[xflux_for_H_mfi]),
			   DATA_AND_LIMITS((*flux_for_Y[0])[xflux_for_Y_mfi]),
			   DATA_AND_LIMITS((*flux_for_H[1])[yflux_for_H_mfi]),
			   DATA_AND_LIMITS((*flux_for_Y[1])[yflux_for_Y_mfi]));
    }
}
