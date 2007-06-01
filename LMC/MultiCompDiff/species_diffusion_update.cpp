// formerly differential_spec_diffusion_update

void
HeatTransfer::species_diffusion_update (Real dt,
					int  corrector)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::species_diffusion_update()");

    const Real strt_time = ParallelDescriptor::second();

    if (hack_nospecdiff)
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "... HACK!!! skipping spec diffusion " << std::endl;

        if (!corrector)
            MultiFab::Copy(get_new_data(State_Type),get_old_data(State_Type),first_spec,first_spec,nspecies,0);

        for (int d = 0; d < BL_SPACEDIM; ++d)
        {
            SpecDiffusionFluxn[d]->setVal(0,0,nspecies);
            SpecDiffusionFluxnp1[d]->setVal(0,0,nspecies);
        }
        for (int comp = 0; comp < nspecies; ++comp)
            spec_diffusion_flux_computed[comp] = HT_Diffusion;

        return;
    }
    //
    // Build single component edge-centered array of MultiFabs for fluxes
    //
    MultiFab** fluxSCn;
    MultiFab** fluxSCnp1;
    const int nGrow   = 0;
    const int nCompSC = 1;
    const int sCompSC = 0;
    diffusion->allocFluxBoxesLevel(fluxSCn,  nGrow,nCompSC);
    diffusion->allocFluxBoxesLevel(fluxSCnp1,nGrow,nCompSC);
    //
    // Set diffusion solve mode
    //
    Diffusion::SolveMode solve_mode = Diffusion::PREDICTOR;
    //
    // Do implicit c-n solve for each scalar...but dont reflux.
    // Save the fluxes, coeffs and source term, we need 'em for later
    //
    const int sCompY = first_spec;
    const int nCompY = nspecies;
    MultiFab* alpha = 0; // Allocate lazily
    MultiFab delta_rhs(grids, nCompY, nGrow);
    MultiFab **betan, **betanp1;
    diffusion->allocFluxBoxesLevel(betan  ,nGrow,nCompY);
    diffusion->allocFluxBoxesLevel(betanp1,nGrow,nCompY);
    Array<int> rho_flag(nCompY,0);
    
    MultiFab *alphaSC, *delta_rhsSC, **betanSC, **betanp1SC;

    const MultiFab* Rh = get_rho_half_time();

    for (int sigma = 0; sigma < nCompY; ++sigma)
    {
	const int state_ind = sCompY + sigma;

	diffuse_scalar_setup(dt, state_ind, &rho_flag[sigma], delta_rhsSC,
			     alphaSC, betanSC, betanp1SC);
	
	if (state_ind == sCompY && alphaSC)
	{
	    alpha = new MultiFab(grids, nCompY, nGrow);
	}
	else
	{
	    if ((!alphaSC) ^ !alpha)
		BoxLib::Error("All diff-diffusion must be of same form");
	}   
	//    
	// Nab a copy of the coeffs, call diffuser
	//
	if (alphaSC)
	    MultiFab::Copy(*alpha,*alphaSC,sCompSC,sigma,nCompSC,nGrow);
	
	for (int d=0; d<BL_SPACEDIM; ++d)
	{
	    if (betanSC)
		MultiFab::Copy(  *betan[d],  *betanSC[d],sCompSC,sigma,nCompSC,nGrow);
	    if (betanp1SC)
		MultiFab::Copy(*betanp1[d],*betanp1SC[d],sCompSC,sigma,nCompSC,nGrow);
	}
	//
	// Make sure we've got a place for delta_rhs...predictor will dump any
	// explicit updates taken before this routine into delta_rhs for later.
	//
	if (delta_rhsSC)
	{
	    MultiFab::Copy(delta_rhs,*delta_rhsSC,sCompSC,sigma,nCompSC,nGrow);
	}
	else
	{
	    delta_rhs.setVal(0,sigma,1);
	}
	//
	// Clean up single-component stuff, then diffuse the scalar
	//
	diffuse_cleanup(delta_rhsSC, betanSC, betanp1SC, alphaSC);

	diffusion->diffuse_scalar(dt,state_ind,be_cn_theta,Rh,rho_flag[sigma],
                                  fluxSCn,fluxSCnp1,sigma,&delta_rhs,alpha,
                                  betan,betanp1,solve_mode);
	//
	// Pull fluxes into flux array
	//
	for (int d = 0; d < BL_SPACEDIM; ++d)
        {
	    MultiFab::Copy(*SpecDiffusionFluxn[d],  *fluxSCn[d],  sCompSC,sigma,nCompSC,nGrow);
	    MultiFab::Copy(*SpecDiffusionFluxnp1[d],*fluxSCnp1[d],sCompSC,sigma,nCompSC,nGrow);
        }
	spec_diffusion_flux_computed[sigma] = HT_Diffusion;
    }
    //
    // Modify update/fluxes to preserve flux sum = 0, compute new update and
    // leave modified fluxes in level data.  Do this in two stages, first for
    // the explicit fluxes, then the implicit ones (so send in rhs=0 for the
    // second one...).
    //
    const int  dataComp  = 0;
    const Real prev_time = state[State_Type].prevTime();
    const Real cur_time  = state[State_Type].curTime();

    adjust_spec_diffusion_update(get_new_data(State_Type),&get_old_data(State_Type),
				 sCompY,dt,prev_time,rho_flag,Rh,dataComp,
                                 &delta_rhs,alpha,betan);
    adjust_spec_diffusion_update(get_new_data(State_Type),&get_new_data(State_Type),
				 sCompY,dt,cur_time,rho_flag,Rh,dataComp,0,
                                 alpha,betanp1);
    //
    // Now do reflux with new, improved fluxes
    //
    if (do_reflux && corrector)
    {
        FArrayBox fluxtot;

	for (int d = 0; d < BL_SPACEDIM; d++)
	{
	    for (MFIter fmfi(*SpecDiffusionFluxn[d]); fmfi.isValid(); ++fmfi)
	    {
                const Box& ebox = (*SpecDiffusionFluxn[d])[fmfi].box();

                fluxtot.resize(ebox,nCompY);
                fluxtot.copy((*SpecDiffusionFluxn[d])[fmfi], ebox,0,ebox,0,nCompY);
                fluxtot.plus((*SpecDiffusionFluxnp1[d])[fmfi],ebox,0,0,nCompY);

		if (level < parent->finestLevel())
		    getLevel(level+1).getViscFluxReg().CrseInit(fluxtot,ebox,d,0,sCompY,nCompY,-dt);

		if (level > 0)
		    getViscFluxReg().FineAdd(fluxtot,d,fmfi.index(),0,sCompY,nCompY,dt);
	    }
	}

	if (level < parent->finestLevel())
	    getLevel(level+1).getViscFluxReg().CrseInitFinish();
    }
    //
    // Clean up memory
    //
    if (alpha)
	delete alpha;
    diffusion->removeFluxBoxesLevel(fluxSCn);
    diffusion->removeFluxBoxesLevel(fluxSCnp1);
    diffusion->removeFluxBoxesLevel(betan);
    diffusion->removeFluxBoxesLevel(betanp1);

    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_time = ParallelDescriptor::second() - strt_time;

    ParallelDescriptor::ReduceRealMax(run_time,IOProc);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "HeatTransfer::species_diffusion_update(): time: " << run_time << std::endl;
}
