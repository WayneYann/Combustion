void
HeatTransfer::calc_divu (Real      time,
                         Real      dt,
                         MultiFab& divu)
{
    if (ParallelDescriptor::IOProcessor())
	std::cout << "JFG: at top of calc_divu\n" << std::flush;

    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::calc_divu()");
    //
    // Get Mwmix, cpmix and pressure
    //
    const int nGrow = 0;

    int sCompR, sCompT, sCompY, sCompP, sCompCp, sCompMw;
    sCompR=sCompT=sCompY=sCompP=sCompCp=sCompMw=0;	

    //
    // mcdd: get Y,T visc terms together, use in place of individual calls below
    //
    MultiFab mcViscTerms;
    const int vtCompT = nspecies; // T terms stacked on Y terms
    const int vtCompY = 0;
    if (do_mcdd)
    {
        mcViscTerms.define(grids,nspecies+1,nGrow,Fab_allocate);
        compute_mcdd_visc_terms(mcViscTerms,vtCompY,mcViscTerms,vtCompT,time,nGrow,DDOp::DD_Temp);
    }

    MultiFab rho(grids,1,nGrow);
    MultiFab temp(grids,1,nGrow);
    
    const MultiFab& Rho_time = get_rho(time);
    MFIter          Rho_mfi(Rho_time);
    
    for (FillPatchIterator Temp_fpi(*this,temp,nGrow,time,State_Type,Temp,1);
         Rho_mfi.isValid() && Temp_fpi.isValid();
         ++Rho_mfi, ++Temp_fpi)
    {
        const int i = Rho_mfi.index();
        
        rho[i].copy(Rho_time[i],0,sCompR,1);
        temp[i].copy(Temp_fpi(),0,sCompT,1);
    }
    //
    // Note that state contains rho*species, so divide species by rho.
    //
    MultiFab species;
    
    species.define(grids,nspecies,nGrow,Fab_allocate);
    
    FArrayBox tmp;
    
    for (FillPatchIterator Spec_fpi(*this,species,nGrow,time,State_Type,first_spec,nspecies);
         Spec_fpi.isValid();
         ++Spec_fpi)
    {
        const int i = Spec_fpi.index();
        
        species[i].copy(Spec_fpi(),0,sCompY,nspecies);
        
        tmp.resize(grids[i],1);
        tmp.copy(rho[i],sCompR,0,1);
        tmp.invert(1);
        
        for (int ispecies = 0; ispecies < nspecies; ispecies++)
            species[i].mult(tmp,0,ispecies,1);
    }
    
    MultiFab mwmix(grids,1,nGrow);
    MultiFab cp(grids,1,nGrow);
    MultiFab p(grids,1,nGrow);
    
    for (MFIter Rho_mfi(rho); Rho_mfi.isValid(); ++Rho_mfi)
    {
        const int  iGrid = Rho_mfi.index();
        const Box& box   = Rho_mfi.validbox();
        
        BL_ASSERT(box == grids[iGrid]);
        
        getChemSolve().getMwmixGivenY(mwmix[iGrid],species[iGrid],
                                      box,sCompY,sCompMw);
        getChemSolve().getCpmixGivenTY(cp[iGrid],temp[iGrid],species[iGrid],
                                       box,sCompT,sCompY,sCompCp);
        getChemSolve().getPGivenRTY(p[iGrid],
                                    rho[iGrid],temp[iGrid],species[iGrid],
                                    box,sCompR,sCompT,sCompY,sCompP);
    }
    //
    // divu = 1/T DT/dt + W_mix * sum (1/W)DY/Dt
    //
    MultiFab visc_terms(grids,1,1);
    MultiFab delta_divu(grids,1,nGrow);
    //
    // Compute rho*DT/Dt as
    //
    //   1/c_p (div lambda grad T + sum_l rho D grad h_l dot grad Y_l)
    //
    
    if (do_mcdd)
    {
        MultiFab::Copy(divu,mcViscTerms,vtCompT,0,1,nGrow);
    }
    else
    {
        getViscTerms(visc_terms,Temp,1,time);
        MultiFab::Copy(divu,visc_terms,0,0,1,nGrow);
    }

    for (MFIter Divu_mfi(divu); Divu_mfi.isValid(); ++Divu_mfi)
    {
        const int iGrid = Divu_mfi.index();
        divu[iGrid].divide(rho[iGrid],grids[iGrid],sCompR,0,1);
        divu[iGrid].divide(temp[iGrid],grids[iGrid],sCompT,0,1);
    }

    delta_divu.setVal(0.0);

    const Array<Real> mwt = getChemSolve().speciesMolecWt();
    MultiFab spec_visc_terms(grids,nspecies,0);
    if (do_mcdd)
    {
        MultiFab::Copy(spec_visc_terms,mcViscTerms,vtCompY,0,nspecies,0);
    }
    else
    {
        getViscTerms(spec_visc_terms,first_spec,nspecies,time);
    }

    for (MFIter mfi(spec_visc_terms); mfi.isValid(); ++mfi)
    {
        const int iGrid = mfi.index();
        for (int comp = 0; comp < nspecies; ++comp)
        {
            spec_visc_terms[mfi].mult(1.0/mwt[comp],comp,1);
            delta_divu[mfi].plus(spec_visc_terms[mfi],grids[iGrid],comp,0,1);
        }
    }
	    
    for (MFIter Divu_mfi(divu); Divu_mfi.isValid(); ++Divu_mfi)
    {
        const int  iGrid = Divu_mfi.index();
        const Box& box   = Divu_mfi.validbox();
        delta_divu[iGrid].divide(rho[iGrid],box,sCompR,0,1);
        delta_divu[iGrid].mult(mwmix[iGrid],box,0,0,1);
        divu[iGrid].plus(delta_divu[iGrid],box,0,0,1);
    }

    if (dt > 0.0)
    {
        //
        // Increment divu by
        //    sum_l (h_l/(c_p*T) - mw_mix/mw_l)*delta Y_l/dt 
        // (i.e., Y_l"dot")
        //
        const Array<Real> mwt = getChemSolve().speciesMolecWt();
        MultiFab h(grids,nspecies,nGrow);
        const int sCompH = 0;
        for (MFIter mfi(h); mfi.isValid(); ++mfi)
        {
            getChemSolve().getHGivenT(h[mfi],temp[mfi],mfi.validbox(),sCompT,sCompH);
        }

        for (FillPatchIterator Ydot_fpi(*this,delta_divu,0,time,Ydot_Type,0,nspecies);
             Ydot_fpi.isValid();
             ++Ydot_fpi)
        {
            for (int istate = first_spec; istate <= last_spec; istate++)
            {
                const int ispec = istate-first_spec;
                const int i     = Ydot_fpi.index();

                delta_divu[i].copy(h[i],ispec,0,1);
                delta_divu[i].divide(cp[i]);
                delta_divu[i].divide(temp[i]);
                delta_divu[i].mult(Ydot_fpi(),ispec,0,1);
                divu[i].plus(delta_divu[i]);

                delta_divu[i].copy(mwmix[i],0,0,1);
                delta_divu[i].divide(mwt[ispec]);
                delta_divu[i].mult(Ydot_fpi(),ispec,0,1);
                divu[i].minus(delta_divu[i]);
            }
        }
    }
    if (ParallelDescriptor::IOProcessor())
	std::cout << "JFG: at bottom of calc_divu\n" << std::flush;
}
