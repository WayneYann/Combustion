#include <winstd.H>

#include "CNSReact.H"
#include "CNSReact_F.H"

#ifdef DIFFUSION
#include "Diffusion.H"
#endif

using std::string;

Real
CNSReact::advance (Real time,
                 Real dt,
                 int  iteration,
                 int  ncycle)
{
    Real dt_new = dt;

    dt_new = advance_hydro(time,dt,iteration,ncycle);
 
    Real cur_time = state[State_Type].curTime();
    set_special_tagging_flag(cur_time);

    return dt_new;
}

Real
CNSReact::advance_hydro (Real time,
                       Real dt,
                       int  iteration,
                       int  ncycle)
{
    u_gdnv = new MultiFab[BL_SPACEDIM];
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
	BoxArray edge_grids(grids);
	edge_grids.surroundingNodes(dir);
	u_gdnv[dir].define(edge_grids,1,1,Fab_allocate);
	u_gdnv[dir].setVal(1.e40);
    }
    
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level) {
        //
        // Set reflux registers to zero.
        //
        getFluxReg(level+1).setVal(0.0);
    }

    for (int k = 0; k < NUM_STATE_TYPE; k++) {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    MultiFab grav_vec_old;

    // It's possible for interpolation to create very small negative values for
    //   species so we make sure here that all species are non-negative after this point
    enforce_nonnegative_species(S_old);

    const Real prev_time = state[State_Type].prevTime();

#ifdef REACTIONS
    react_first_half_dt(S_old,time,dt);
#endif

    Real dt_new = dt;

    Real cur_time = state[State_Type].curTime();
        
        //
        // Get pointers to Flux registers, or set pointer to zero if not there.
        //
        FluxRegister *fine    = 0;
        FluxRegister *current = 0;
        
        if (do_reflux && level < finest_level)
            fine = &getFluxReg(level+1);
        if (do_reflux && level > 0)
            current = &getFluxReg(level);

        AmrLevel &levelData = *this;
        Geometry g = levelData.Geom();
        MultiFab levelVolume;
        g.GetVolume(levelVolume,grids,NUM_GROW);
        
        MultiFab levelArea[BL_SPACEDIM];
        for (int i = 0; i < BL_SPACEDIM ; i++)
            g.GetFaceArea(levelArea[i],grids,i,NUM_GROW);

        // Integrate, looping through the grids.
        FArrayBox divu, grid_volume;
        FArrayBox dloga, area[BL_SPACEDIM];
        FArrayBox flux[BL_SPACEDIM];
        
        const Real *dx = geom.CellSize();
        Real courno    = -1.0e+20;

        MultiFab fluxes[BL_SPACEDIM];

        if (do_reflux && fine)
        {
            for (int j = 0; j < BL_SPACEDIM; j++)
            {
                BoxArray ba = S_new.boxArray();
                ba.surroundingNodes(j);
                fluxes[j].define(ba, NUM_STATE, 0, Fab_allocate);
            }
        }

       // Define the gravity vector so we can pass this to ca_umdrv.
       MultiFab grav_vector(grids,BL_SPACEDIM,3);
       grav_vector.setVal(0.);

       MultiFab ext_src_old(grids,NUM_STATE,1,Fab_allocate);
       ext_src_old.setVal(0.0);

       if (add_ext_src)
          getOldSource(prev_time,dt,ext_src_old);

#ifdef DIFFUSION
       MultiFab OldTempDiffTerm(grids,1,1);
       add_diffusion_to_old_source(ext_src_old,OldTempDiffTerm,prev_time);
#endif
       ext_src_old.FillBoundary();
        
       for (FillPatchIterator fpi(*this, S_new, NUM_GROW,
                                   time, State_Type, 0, NUM_STATE);
             fpi.isValid();
             ++fpi)
        {
            int mfiindex = fpi.index();

            Box bx(fpi.UngrownBox());

            Box bx_g4(BoxLib::grow(bx,4));
            // Create FAB for extended grid values (including boundaries) and fill.
            FArrayBox &state = fpi();
            FArrayBox &stateout = S_new[fpi];
            
            grid_volume.resize(bx_g4,1);
            grid_volume.copy(levelVolume[fpi]);

            for (int i = 0; i < BL_SPACEDIM ; i++) {
                area[i].resize(BoxLib::surroundingNodes(bx_g4,i));
                area[i].copy(levelArea[i][mfiindex]);
            }
            
#if (BL_SPACEDIM <=2)
            dloga.resize(bx_g4);
            dloga.copy(dLogArea[0][mfiindex]);
#endif

            // Allocate fabs for fluxes.
            for (int i = 0; i < BL_SPACEDIM ; i++)  
                flux[i].resize(BoxLib::surroundingNodes(bx,i),NUM_STATE);

            const int*  domain_lo = geom.Domain().loVect();
            const int*  domain_hi = geom.Domain().hiVect();

            Real cflLoc = -1.0e+20;
            int is_finest_level = 0;
            if (level == finest_level) is_finest_level = 1;
            BL_FORT_PROC_CALL(CA_UMDRV,ca_umdrv)
                (&is_finest_level,&time,
                 bx.loVect(), bx.hiVect(),
                 domain_lo, domain_hi,
                 BL_TO_FORTRAN(state), BL_TO_FORTRAN(stateout),
		 BL_TO_FORTRAN(u_gdnv[0][fpi]),
#if (BL_SPACEDIM >= 2)
		 BL_TO_FORTRAN(u_gdnv[1][fpi]),
#endif
#if (BL_SPACEDIM == 3)
		 BL_TO_FORTRAN(u_gdnv[2][fpi]),
#endif
                 BL_TO_FORTRAN(ext_src_old[fpi]),
                 BL_TO_FORTRAN(grav_vector[fpi]), 
                 dx, &dt,
                 D_DECL(BL_TO_FORTRAN(flux[0]), 
                        BL_TO_FORTRAN(flux[1]), 
                        BL_TO_FORTRAN(flux[2])), 
                 D_DECL(BL_TO_FORTRAN(area[0]), 
                        BL_TO_FORTRAN(area[1]), 
                        BL_TO_FORTRAN(area[2])), 
#if (BL_SPACEDIM < 3) 
                 BL_TO_FORTRAN(dloga), 
#endif
                 BL_TO_FORTRAN(grid_volume), 
                 &cflLoc,verbose);

            if (do_reflux)
            {
                if (fine)
                {
                    for (int i = 0; i < BL_SPACEDIM ; i++)
                        fluxes[i][mfiindex].copy(flux[i]);
                }
                if (current)
                {
                    for (int i = 0; i < BL_SPACEDIM ; i++)
                        current->FineAdd(flux[i],i,mfiindex,0,0,NUM_STATE,1);
                }
            }
            courno = std::max(courno,cflLoc);
        }

        if (do_reflux && fine)
        {
            for (int i = 0; i < BL_SPACEDIM ; i++)
                fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1);
        }

        ParallelDescriptor::ReduceRealMax(courno);

        if (courno > 1.0) {
           if (ParallelDescriptor::IOProcessor()) 
              std::cout << "OOPS -- EFFECTIVE CFL AT THIS LEVEL " << level << " IS " << courno << std::endl;
           BoxLib::Abort("CFL is too high at this level -- go back to a checkpoint and restart with lower cfl number");
        }
        
        dt_new = dt/courno;

        if (add_ext_src)  
        {
           getOldSource(prev_time,dt,ext_src_old);
           ext_src_old.mult(-0.5*dt);

           // Must compute new temperature in case it is needed in the source term evaluation
           computeTemp(S_new);

           // Compute source at new time (no ghost cells needed)
           MultiFab ext_src_new(grids,NUM_STATE,0,Fab_allocate);
           ext_src_new.setVal(0.0);

           getNewSource(prev_time,cur_time,dt,ext_src_new);
           ext_src_new.mult(0.5*dt);

           // Subtract off half of the old source term, and add half of the new.
           MultiFab::Add(S_new,ext_src_old,0,0,S_new.nComp(),0);
           MultiFab::Add(S_new,ext_src_new,0,0,S_new.nComp(),0);
           computeTemp(S_new);
        }

#ifdef DIFFUSION
     time_center_diffusion(S_new, OldTempDiffTerm, cur_time, dt);
#endif

     reset_internal_energy(S_new);

#ifdef REACTIONS
    react_second_half_dt(S_new,cur_time,dt);
#endif

    delete [] u_gdnv;

    return dt_new;
}
