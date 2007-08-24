void
HeatTransfer::rk_diffusion_operator (const Real time,
				     const Real scale_factor,
				     MultiFab *& div_of_flux_for_H,
				     MultiFab *& div_of_flux_for_Y,
				     MultiFab **& flux_for_H,
				     MultiFab **& flux_for_Y)
{
/*
  evaluate the "operator" (extensive fluxes and their divergences) for the 
  Runge-Kutta implementation of the diffusion update using the mixture
  averaged or multicomponent formulation.
  
  the divergences are scaled by 1/vol so they do approximate the analytic
  divergence.  when used in the Runge-Kutta formula their additional scaling,
  given by scale_factor, should be set to - dt.  the minus occurs because
  the standard ODE is written y_prime = f(y), that is, the divergences are
  moved to the opposite side of the equation from the time derivative.  when
  used in getViscTerms the scale_factor should be set to -1.

  Real time                       ! INPUT time, either prev_time or cur_time
  Real scale_factor               ! INPUT scale factor
  MultiFab *& div_of_flux_for_H   ! OUTPUT divergence of the flux for rho H weighted by scale_factor / vol
  MultiFab *& div_of_flux_for_Y   ! OUTPUT divergence of the flux for rho Y weighted by scale_factor / vol
  MultiFab **& flux_for_H         ! OUTPUT extensive rho H flux
  MultiFab **& flux_for_Y         ! OUTPUT extensive rho Y flux
*/

    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::rk_diffusion_operator()");

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

    // allocate OUTPUT multifabs for fluxes and updates
    int ngrow = 0;
    diffusion->allocFluxBoxesLevel(flux_for_H,ngrow,1);
    diffusion->allocFluxBoxesLevel(flux_for_Y,ngrow,nspecies);
    div_of_flux_for_H = new MultiFab(grids,1,ngrow);
    div_of_flux_for_Y = new MultiFab(grids,nspecies,ngrow);

    // loop over fabs in the state at the specified time
    MultiFab dummy (grids,1,0,Fab_noallocate);
    MFIter div_of_flux_for_H_mfi(*div_of_flux_for_H);
    MFIter div_of_flux_for_Y_mfi(*div_of_flux_for_Y);
    MFIter xflux_for_H_mfi(*flux_for_H[0]);
    MFIter yflux_for_H_mfi(*flux_for_H[1]);
    MFIter xflux_for_Y_mfi(*flux_for_Y[0]);
    MFIter yflux_for_Y_mfi(*flux_for_Y[1]);
    ngrow = 1;
    for (FillPatchIterator state_fpi (*this, dummy, ngrow, time, State_Type, 0, ncomps);
         state_fpi.isValid();
         ++state_fpi, 
	     ++div_of_flux_for_H_mfi, 
	     ++div_of_flux_for_Y_mfi,
	     ++xflux_for_H_mfi,
	     ++yflux_for_H_mfi,
	     ++xflux_for_Y_mfi,
	     ++yflux_for_Y_mfi)
    {
        BL_ASSERT (
	    div_of_flux_for_H_mfi.isValid() &&
	    div_of_flux_for_Y_mfi.isValid() &&
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
c   * areax, DIMS(areax),               ! INPUT areas of the faces perendicular to x axis
c   * areay, DIMS(areay),               ! INPUT areas of the faces perpendicular to y axis
c     bc,                               ! INPUT boundary condition array for all comps
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
c     scale_factor,                     ! INPUT scale_factor
c     state, DIMS(state),               ! INPUT all variables in the state
c   * volume, DIMS(volume),             ! INPUT volumes of the cells
c     div_of_flux_for_H, DIMS(div_of_flux_for_H), ! OUTPUT divergence of the flux for rho H
c     div_of_flux_for_Y, DIMS(div_of_flux_for_Y), ! OUTPUT divergences of the fluxes for rho Y
c     xflux_for_H, DIMS(xflux_for_H),   ! OUTPUT extensive x fluxes for rho H
c     xflux_for_Y, DIMS(xflux_for_Y),   ! OUTPUT extensive x fluxes for rho Y
c     yflux_for_H, DIMS(yflux_for_H),   ! OUTPUT extensive y fluxes for rho H
c     yflux_for_Y, DIMS(yflux_for_Y),   ! OUTPUT extensive y fluxes for rho Y
c
c     * these arguments are not used
*/

#define DATA_AND_LIMITS(foo) foo.dataPtr(),foo.loVect()[0],foo.loVect()[1],foo.hiVect()[0],foo.hiVect()[1]

	if (!rk_mixture_averaged)
	{
	    // multicomponent is the default
	    FORT_RK_MULTICOMPONENT
		(geom.Domain().loVect(), geom.Domain().hiVect(), 
		 grids[idx].loVect(), grids[idx].hiVect(),
		 DATA_AND_LIMITS(area[0][idx]),
		 DATA_AND_LIMITS(area[1][idx]),
		 bc.dataPtr(),
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
		 &scale_factor,
		 DATA_AND_LIMITS(state_fpi()),
		 DATA_AND_LIMITS(volume[idx]),
		 DATA_AND_LIMITS((*div_of_flux_for_H)[div_of_flux_for_H_mfi]),
		 DATA_AND_LIMITS((*div_of_flux_for_Y)[div_of_flux_for_Y_mfi]),
		 DATA_AND_LIMITS((*flux_for_H[0])[xflux_for_H_mfi]),
		 DATA_AND_LIMITS((*flux_for_Y[0])[xflux_for_Y_mfi]),
		 DATA_AND_LIMITS((*flux_for_H[1])[yflux_for_H_mfi]),
		 DATA_AND_LIMITS((*flux_for_Y[1])[yflux_for_Y_mfi]));
	}
	else
	{
	    // do mixture averaged
	    FORT_RK_MIXTURE_AVERAGED
		(geom.Domain().loVect(), geom.Domain().hiVect(), 
		 grids[idx].loVect(), grids[idx].hiVect(),
		 DATA_AND_LIMITS(area[0][idx]),
		 DATA_AND_LIMITS(area[1][idx]),
		 bc.dataPtr(),
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
		 &scale_factor,
		 DATA_AND_LIMITS(state_fpi()),
		 DATA_AND_LIMITS(volume[idx]),
		 DATA_AND_LIMITS((*div_of_flux_for_H)[div_of_flux_for_H_mfi]),
		 DATA_AND_LIMITS((*div_of_flux_for_Y)[div_of_flux_for_Y_mfi]),
		 DATA_AND_LIMITS((*flux_for_H[0])[xflux_for_H_mfi]),
		 DATA_AND_LIMITS((*flux_for_Y[0])[xflux_for_Y_mfi]),
		 DATA_AND_LIMITS((*flux_for_H[1])[yflux_for_H_mfi]),
		 DATA_AND_LIMITS((*flux_for_Y[1])[yflux_for_Y_mfi]));
	    std::cout << "                            x sum = " << ((*flux_for_H[0])[xflux_for_H_mfi]).norm (1, 0, 1) << std::endl;
	    std::cout << "                            y sum = " << ((*flux_for_H[1])[yflux_for_H_mfi]).norm (1, 0, 1) << std::endl;
	    
	}

        // BoxLib::Abort("JFG: stopping here for now");
    }
}
