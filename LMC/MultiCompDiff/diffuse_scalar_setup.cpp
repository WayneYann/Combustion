void
HeatTransfer::diffuse_scalar_setup (Real        dt,
                                    int         sigma,
                                    int*        rho_flag, 
                                    MultiFab*&  delta_rhs,
                                    MultiFab*&  alpha, 
                                    MultiFab**& betan,
                                    MultiFab**& betanp1)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::diffuse_scalar_setup()");
    //
    // Do setup for implicit c-n solve for an arbitrary scalar.
    //
    // Note: should be ok for variable c_p.
    //
    const Real prev_time = state[State_Type].prevTime();

    NavierStokes::diffuse_scalar_setup(dt, sigma, rho_flag, 
                                       delta_rhs, alpha, betan, betanp1);
    alpha     = 0;
    delta_rhs = 0;
    betan     = 0;
    betanp1   = 0;
   
    if (sigma == Temp)
    {
        (*rho_flag) = 1;
        diffuse_temp_setup(prev_time,dt,delta_rhs,alpha); 
    }
    else if (sigma == RhoH)
    {
        (*rho_flag) = 2;
        diffuse_rhoh_setup(prev_time,dt,delta_rhs); 
    }
    else if (sigma >= first_spec && sigma <= last_spec)
    {
        (*rho_flag) = 2;
        // formerly: diffuse_spec_setup(sigma,prev_time,dt,delta_rhs); 
	// Chemistry split, no source terms
	delta_rhs = new MultiFab(grids,1,0);
	delta_rhs->setVal(0);
    }

    diffusion->allocFluxBoxesLevel(betan);
    diffusion->allocFluxBoxesLevel(betanp1);
    getDiffusivity(betan, prev_time, sigma, 0, 1);
    getDiffusivity(betanp1, prev_time+dt, sigma, 0, 1);
}
