// this is some of Joe's debug stuff
// define this to get the orginal visc terms
// comment this out to get the new operator
// #define ORIGINAL_GETVISCTERMS

void
HeatTransfer::getViscTerms (MultiFab& visc_terms,
                            int       src_comp, 
                            int       num_comp,
                            Real      time)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::getViscTerms()");
    //
    // Load "viscous" terms, starting from component = 0.
    //
    // JFG: for species, this procedure returns the *negative* of the divergence of 
    // of the diffusive fluxes.  specifically, in the mixture averaged case, the
    // diffusive flux vector for species k is
    //
    //       j_k = - rho D_k,mix grad Y_k
    //
    // so the divergence of the flux is div dot j_k.  instead this procedure returns 
    // - div dot j_k.
    //
    // note the fluxes used in the code are extensive, that is, scaled by the areas
    // of the cell edges.  the calculation of the divergence is the sum of un-divided
    // differences of the extensive fluxes, all divided by volume.  so the effect is 
    // to give the true divided difference approximation to the divergence of the
    // intensive flux.
    //
    const int  last_comp = src_comp + num_comp - 1;
    int        icomp     = src_comp; // This is the current related state comp.
    int        load_comp = 0;        // Comp for result of current calculation.
    MultiFab** vel_visc  = 0;        // Potentially reused, raise scope
    const int  nGrow     = visc_terms.nGrow();
    //
    // Get Div(tau) from the tensor operator, if velocity and have non-const viscosity
    //
    if (src_comp < BL_SPACEDIM)
    {
        if (src_comp != Xvel || num_comp < BL_SPACEDIM)
            BoxLib::Error("tensor v -> getViscTerms needs all v-components at once");

        diffusion->allocFluxBoxesLevel(vel_visc);
        getViscosity(vel_visc, time);
        diffusion->getTensorViscTerms(visc_terms,time,0,vel_visc);
        icomp = load_comp = BL_SPACEDIM;
    }
    //
    // For Le != 1, need to get visc terms in a block
    //
    if (!unity_Le)
    {
	const int non_spec_comps = std::min(num_comp,
                                            std::max(0,first_spec-src_comp) + std::max(0,last_comp-last_spec));
	const bool has_spec = src_comp <= last_spec && last_comp >= first_spec;
	BL_ASSERT( !has_spec || (num_comp>=nspecies+non_spec_comps));
	if (has_spec)
	{
	    const int sCompY = first_spec - src_comp + load_comp;
            const int sCompT = Temp - src_comp + load_comp;
            if (do_mcdd)
            {
		// code for Marc's multicomponent diffusion

                BL_ASSERT(sCompT > 0 && sCompT < visc_terms.nComp());
                compute_mcdd_visc_terms(visc_terms,sCompY,visc_terms,sCompT,time,nGrow,DDOp::DD_Temp);
            }
#ifndef ORIGINAL_GETVISCTERMS
            else if (do_rk_diffusion)
	    {
		// code for Joe's Runge-Kutta diffusion

		if (ParallelDescriptor::IOProcessor())
		    std::cout << "JFG: at top of getViscTerms do_rk_diffusion block\n" << std::flush;

		// apply the diffusion operator to the old or new state to get updates
		// associated with that state.  only the species updates are used here, 
		// but space is still needed for the others.
		MultiFab* div_of_flux_for_H;
		MultiFab* div_of_flux_for_Y;
		MultiFab** flux_for_H;
		MultiFab** flux_for_Y;

		// per the comments above, pass - 1.0 as the scaling argument.
		rk_diffusion_operator (time,
				       - 1.0,
				       div_of_flux_for_H,
				       div_of_flux_for_Y,
				       flux_for_H,
				       flux_for_Y);

		// examination of visc_terms.nGrow() reveals getViscTerms is called twice,
		// once with 1 growth cell, and once with 0.  the original implementation,
		// in the following else block, compute_differential_diffusion_terms 
                // contains a remark that it does not fill growth cells, so we also do
		// not fill them here.  note the arguments for MultiFab:: procedures are 
		// (dst, src, srccomp, dstcomp, ncomp, nghost);
		MultiFab::Copy (visc_terms, *div_of_flux_for_Y, 0, sCompY, nspecies, 0);

		// delete the space for fluxes and updates
		delete div_of_flux_for_H;
		delete div_of_flux_for_Y;
		diffusion->removeFluxBoxesLevel (flux_for_H);
		diffusion->removeFluxBoxesLevel (flux_for_Y);

		if (ParallelDescriptor::IOProcessor())
		    std::cout << "JFG: at bottom of getViscTerms do_rk_diffusion block\n" << std::flush;
	    }
#endif
	    else {
		// code for the original implementation

		compute_differential_diffusion_terms(visc_terms,sCompY,time);
            }
        }
    }
    //
    // Now, do the rest.
    // Get Div(visc.Grad(state)) or Div(visc.Grad(state/rho)), as appropriate.
    //
    for ( ; icomp <= last_comp; load_comp++, icomp++)
    {
	const bool is_spec  = icomp >= first_spec && icomp <= last_spec;
	if ( !(is_spec && !unity_Le) )
	{
	    if (icomp == Temp)
	    {
                if (!do_mcdd) // Because in this case, was done above
                    getTempViscTerms(visc_terms,Temp-load_comp,time);
	    }
	    else if (icomp == RhoH)
	    {
                if (do_mcdd)
                    // What to do here?  
                    BoxLib::Abort("do we really want to get RhoH VT when do_mcdd?");
	    }
	    else
	    {
		const int  rho_flag = Diffusion::set_rho_flag(diffusionType[icomp]);
		MultiFab** beta     = 0;

		if (icomp == Density)
                {
                    visc_terms.setVal(0.0,load_comp,1);
                }
                else
                {
                    //
                    // Assume always variable viscosity / diffusivity.
                    //
                    diffusion->allocFluxBoxesLevel(beta);
                    getDiffusivity(beta, time, icomp, 0, 1);

                    diffusion->getViscTerms(visc_terms,icomp-load_comp,
                                            icomp,time,rho_flag,0,beta);
                    
                    diffusion->removeFluxBoxesLevel(beta);
                }
	    }
	}
    }
    //
    // Add Div(u) term if desired, if this is velocity, and if Div(u) is nonzero
    // If const-visc, term is mu.Div(u)/3, else it's -Div(mu.Div(u).I)*2/3
    //
    if (src_comp < BL_SPACEDIM && S_in_vel_diffusion)
    {
        if (num_comp < BL_SPACEDIM)
            BoxLib::Error("getViscTerms() need all v-components at once");
    
        MultiFab divmusi(grids,BL_SPACEDIM,1);
        //
	// Assume always using variable viscosity.
        //
	diffusion->compute_divmusi(time,vel_visc,divmusi); // pre-computed visc above
	divmusi.mult((-2./3.),0,BL_SPACEDIM,0);
        visc_terms.plus(divmusi,Xvel,BL_SPACEDIM,0);
    }
    //
    // Clean up your mess ...
    //
    if (vel_visc)
        diffusion->removeFluxBoxesLevel(vel_visc);
    //
    // Ensure consistent grow cells
    //    
    if (nGrow > 0)
    {
        for (MFIter mfi(visc_terms); mfi.isValid(); ++mfi)
        {
            FArrayBox& vt  = visc_terms[mfi];
            const Box& box = mfi.validbox();
            FORT_VISCEXTRAP(vt.dataPtr(),ARLIM(vt.loVect()),ARLIM(vt.hiVect()),
                            box.loVect(),box.hiVect(),&num_comp);
        }
        visc_terms.FillBoundary(0,num_comp);
        //
        // Note: this is a special periodic fill in that we want to
        // preserve the extrapolated grow values when periodic --
        // usually we preserve only valid data.  The scheme relies on
        // the fact that there is good data in the "non-periodic" grow cells.
        // ("good" data produced via VISCEXTRAP above)
        //
        geom.FillPeriodicBoundary(visc_terms,0,num_comp,true);
    }
}
