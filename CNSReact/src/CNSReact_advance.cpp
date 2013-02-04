#include <winstd.H>

#include "CNSReact.H"
#include "CNSReact_F.H"

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

    const Real prev_time = state[State_Type].prevTime();
    const Real  cur_time = state[State_Type].curTime();

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    if (S_old.contains_nan(Density,S_old.nComp(),0))
    {
        for (int i = 0; i < S_old.nComp(); i++)
        {
            if (S_old.contains_nan(Density+i,1,0))
            {
                std::cout << "Testing component i for NaNs: " << i << std::endl;
                BoxLib::Abort("S_old has NaNs in this component::advance_hydro()");
            }
        }
    }

    // It's possible for interpolation to create very small negative values for
    //   species so we make sure here that all species are non-negative after this point
    enforce_nonnegative_species(S_old);

    Real dt_new = dt;

        
    //
    // Get pointers to Flux registers, or set pointer to zero if not there.
    //
    FluxRegister *fine    = 0;
    FluxRegister *current = 0;
        
    if (do_reflux && level < finest_level) {
      fine = &getFluxReg(level+1);
    }
    if (do_reflux && level > 0) {
      current = &getFluxReg(level);
    }

    AmrLevel &levelData = *this;
    Geometry g = levelData.Geom();
    MultiFab levelVolume;
    g.GetVolume(levelVolume,grids,NUM_GROW);
        
    MultiFab levelArea[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM ; i++) {
      g.GetFaceArea(levelArea[i],grids,i,NUM_GROW);
    }

    // Integrate, looping through the grids.
    FArrayBox divu, grid_volume;
    FArrayBox dloga, area[BL_SPACEDIM];
    FArrayBox flux[BL_SPACEDIM];
        
    const Real *dx = geom.CellSize();
    Real courno    = -1.0e+20;

    MultiFab fluxes[BL_SPACEDIM];

    if (do_reflux && fine) {
      for (int j = 0; j < BL_SPACEDIM; j++) {
	BoxArray ba = S_new.boxArray();
	ba.surroundingNodes(j);
	fluxes[j].define(ba, NUM_STATE, 0, Fab_allocate);
      }
    }

    for (FillPatchIterator fpi(*this, S_new, NUM_GROW,
			       time, State_Type, 0, NUM_STATE);
	 fpi.isValid(); ++fpi) {

      int mfiindex = fpi.index();

      Box bx(fpi.UngrownBox());

      Box bx_g(BoxLib::grow(bx,NUM_GROW));
      // Create FAB for extended grid values (including boundaries) and fill.
      FArrayBox &state = fpi();
      FArrayBox &stateout = S_new[fpi];
            
      grid_volume.resize(bx_g,1);
      grid_volume.copy(levelVolume[fpi]);

      for (int i = 0; i < BL_SPACEDIM ; i++) {
	area[i].resize(BoxLib::surroundingNodes(bx_g,i));
	area[i].copy(levelArea[i][mfiindex]);
      }
            
#if (BL_SPACEDIM <=2)
      dloga.resize(bx_g);
      dloga.copy(dLogArea[0][mfiindex]);
#endif

      // Allocate fabs for fluxes.
      Box bx_g1(BoxLib::grow(bx,1));
      for (int i = 0; i < BL_SPACEDIM ; i++) {
	flux[i].resize(BoxLib::surroundingNodes(bx_g1,i),NUM_STATE);
      }

      const int*  domain_lo = geom.Domain().loVect();
      const int*  domain_hi = geom.Domain().hiVect();

      Real cflLoc = -1.0e+20;
      int is_finest_level = 0;
      if (level == finest_level) is_finest_level = 1;

      BL_FORT_PROC_CALL(CNS_UMDRV,cns_umdrv)
	(bx.loVect(), bx.hiVect(),
	 BL_TO_FORTRAN(state), BL_TO_FORTRAN(stateout),
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

      if (do_reflux) {
	if (fine) {
	  for (int i = 0; i < BL_SPACEDIM ; i++) {
                        fluxes[i][mfiindex].copy(flux[i]);
	  }
	}
	if (current) {
	  for (int i = 0; i < BL_SPACEDIM ; i++) {
	    current->FineAdd(flux[i],i,mfiindex,0,0,NUM_STATE,1);
	  }
	}
      }

      courno = std::max(courno,cflLoc);

    } // end for(fpi...)

    if (do_reflux && fine) {
      for (int i = 0; i < BL_SPACEDIM ; i++) {
	fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1);
      }
    }

    ParallelDescriptor::ReduceRealMax(courno);

    if (courno > 1.0) {
      if (ParallelDescriptor::IOProcessor()) 
	std::cout << "OOPS -- EFFECTIVE CFL AT THIS LEVEL " << level << " IS " << courno << std::endl;
      BoxLib::Abort("CFL is too high at this level -- go back to a checkpoint and restart with lower cfl number");
    }
        
    dt_new = dt/courno;

    if (S_new.contains_nan(Density,S_new.nComp(),0))
    {
        for (int i = 0; i < S_new.nComp(); i++)
        {
            if (S_new.contains_nan(Density + i, 1, 0))
            {
                std::cout << "Testing component i for NaNs: " << i << std::endl;
                BoxLib::Abort("S_new has NaNs in this component::advance_hydro()");
            }
        }
    }

    reset_internal_energy(S_new);

    return dt_new;
}
