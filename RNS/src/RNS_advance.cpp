#include <winstd.H>

#include "RNS.H"
#include "RNS_F.H"

// #include <ArrayView.H>

using std::string;

Real
RNS::advance (Real time,
	      Real dt,
	      int  iteration,
	      int  ncycle)
{
    for (int k = 0; k < NUM_STATE_TYPE; k++)
    {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

    MultiFab& Unew = get_new_data(State_Type);
    MultiFab& Uold = get_old_data(State_Type);

    int fill_boundary_type;

    if (level == 0)
    {
	MultiFab::Copy(Unew, Uold, 0, 0, NUM_STATE, 0);
	fill_boundary_type = use_FillBoundary;
    }
    else
    {
	fill_boundary_type = use_FillPatchIterator;
    }
    fill_boundary(Unew, time, fill_boundary_type);

    if (! chemSolve->isNull)
    {
	// do half-dt chemistry
	advance_chemistry(Unew, 0.5*dt);

	// After chemistry, fill ghosts cells with the post-chemistry state of current level
	fill_boundary(Unew, time, use_FillBoundary);
    }

    // Advance Advection & Diffusion
    advance_AD(Unew, time, dt);

    if (! chemSolve->isNull)
    {
	// fill boundary for chemistry
	fill_boundary_type = use_FillCoarsePatch;
	fill_boundary(Unew, time+dt, fill_boundary_type);

	// do another half-dt chemistry
	advance_chemistry(Unew, 0.5*dt);
    }

    if (Unew.contains_nan(0,NUM_STATE,0))
    {
	for (int i=0; i<NUM_STATE; i++)
	{
	    if (Unew.contains_nan(i, 1, 0))
	    {
		std::cout << "RNS::advance: Testing component i for NaNs: " << i << std::endl;
                BoxLib::Abort("RNS::advance: Has NaNs in this component.");
	    }
	}
    }

    return dt;
}


void
RNS::fill_boundary(MultiFab& U, Real time, int type_in)
{
    if (type_in == no_fill) return;

    int type = type_in;
    if (level == 0 && type_in == use_FillCoarsePatch) type = use_FillBoundary;

    switch (type)
    {
    case use_FillPatchIterator:

	for (FillPatchIterator fpi(*this, U, NUM_GROW, time, State_Type, 0, NUM_STATE);
	     fpi.isValid(); ++fpi)
	{
	    U[fpi].copy(fpi());
	}

	break;

    case use_FillCoarsePatch:  // so that valid region of U will not be touched
    {
	const Box& domain_box = geom.Domain();
	BoxArray grids_g(grids);
	for (int ibox=0; ibox<grids_g.size(); ibox++)
	{
	    const Box b = BoxLib::grow(grids_g[ibox], NUM_GROW) & domain_box;
	    grids_g.set(ibox, b);
	}

	MultiFab Utmp(grids_g, NUM_STATE, 0);
	FillCoarsePatch(Utmp, 0, time, State_Type, 0, NUM_STATE);

	for (MFIter mfi(U); mfi.isValid(); ++mfi)
	{
	    int i = mfi.index();

	    const Box& vbox = grids[i];
	    const Box& gbox = Utmp[i].box();
	    const BoxArray& ba = BoxLib::boxComplement(gbox, vbox);

	    for (int ibox=0; ibox<ba.size(); ibox++)
	    {
		U[i].copy(Utmp[i], ba[ibox]);
	    }
	}
    }

    // no break; so it will go to next case and call FillBoundary

    case use_FillBoundary:

	U.FillBoundary();
	geom.FillPeriodicBoundary(U, true);

	// no break; so it will go to next case and set physical boundaries

    case set_PhysBoundary:

	for (MFIter mfi(U); mfi.isValid(); ++mfi)
	{
	    setPhysBoundaryValues(U[mfi],
				  State_Type,
				  time,
				  0,
				  0,
				  NUM_STATE);
	}

	break;
    }
}


void
RNS::dUdt_AD(MultiFab& U, MultiFab& Uprime, Real time, int fill_boundary_type,
	     FluxRegister* fine, FluxRegister* current, Real dt)
{
    FArrayBox  flux[BL_SPACEDIM];
    MultiFab fluxes[BL_SPACEDIM];
    if (do_reflux && fine)
    {
	for (int idim = 0; idim < BL_SPACEDIM; idim++)
	{
	    BoxArray ba = U.boxArray();
	    ba.surroundingNodes(idim);
	    fluxes[idim].define(ba, NUM_STATE, 0, Fab_allocate);
	}
    }

    fill_boundary(U, time, fill_boundary_type);

    const Real *dx = geom.CellSize();

    for (MFIter mfi(Uprime); mfi.isValid(); ++mfi)
    {
	int i = mfi.index();
	const Box& bx = mfi.validbox();

	for (int idim = 0; idim < BL_SPACEDIM ; idim++)
	{
	    flux[idim].resize(BoxLib::surroundingNodes(bx,idim),NUM_STATE);
	}

	BL_FORT_PROC_CALL(RNS_DUDT_AD,rns_dudt_ad)
	    (bx.loVect(), bx.hiVect(),
	     BL_TO_FORTRAN(U[i]),
	     BL_TO_FORTRAN(Uprime[i]),
	     D_DECL(BL_TO_FORTRAN(flux[0]),
		    BL_TO_FORTRAN(flux[1]),
		    BL_TO_FORTRAN(flux[2])),
	     dx);

	if (do_reflux)
	{
	    if (fine)
	    {
		for (int idim = 0; idim < BL_SPACEDIM ; idim++)
		{
		    fluxes[idim][i].copy(flux[idim]);
		}
	    }

	    if (current)
	    {
		for (int idim = 0; idim < BL_SPACEDIM ; idim++)
		{
		    current->FineAdd(flux[idim],area[idim][i],idim,i,0,0,NUM_STATE,dt);
		}
	    }
	}
    }

    if (do_reflux && fine)
    {
	for (int idim = 0; idim < BL_SPACEDIM ; idim++)
	{
	    fine->CrseInit(fluxes[idim],area[idim],idim,0,0,NUM_STATE,-dt);
	}
    }
}


void
RNS::dUdt_chemistry(const MultiFab& U, MultiFab& Uprime)
{
    BL_ASSERT( ! chemSolve->isNull );

    for (MFIter mfi(U); mfi.isValid(); ++mfi)
    {
	const int   i = mfi.index();
	const Box& bx = mfi.validbox();
	const int* lo = bx.loVect();
	const int* hi = bx.hiVect();

	BL_FORT_PROC_CALL(RNS_DUDT_CHEM, rns_dudt_chem)
	    (lo, hi, BL_TO_FORTRAN(U[i]), BL_TO_FORTRAN(Uprime[i]));
    }
}


// Compute U1 = U2 + c Uprime.
void
RNS::update_rk(MultiFab& U1, const MultiFab& U2, Real c, const MultiFab& Uprime)
{
    for (MFIter mfi(U1); mfi.isValid(); ++mfi)
    {
	const int   i = mfi.index();
	const Box& bx = mfi.validbox();

	U1[i].copy(U2[i], bx);
	U1[i].saxpy(c, Uprime[i]);
    }
}

// Compute U = a Ua + b Ub + c Uprime.
void
RNS::update_rk(MultiFab& U, Real a, const MultiFab& Ua, Real b, const MultiFab& Ub,
	       Real c, const MultiFab& Uprime)
{
    for (MFIter mfi(U); mfi.isValid(); ++mfi)
    {
	const int   i = mfi.index();

	U[i].setVal(0.0);
	U[i].saxpy(a, Ua[i]);
	U[i].saxpy(b, Ub[i]);
	U[i].saxpy(c, Uprime[i]);
    }
}


void
RNS::post_update(MultiFab& U)
{
    for (MFIter mfi(U); mfi.isValid(); ++mfi)
    {
	const int   i = mfi.index();
	const Box& bx = mfi.validbox();
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();

	BL_FORT_PROC_CALL(RNS_ENFORCE_CONSISTENT_Y, rns_enforce_consistent_y)
	    (lo, hi, BL_TO_FORTRAN(U[i]));

	BL_FORT_PROC_CALL(RNS_COMPUTE_TEMP, rns_compute_temp)
	    (lo, hi, BL_TO_FORTRAN(U[i]));
    }
}


void
RNS::advance_chemistry(MultiFab& U, Real dt)
{
    BL_ASSERT( ! chemSolve->isNull );

    for (MFIter mfi(U); mfi.isValid(); ++mfi)
    {
	const int   i = mfi.index();
	const Box& bx = mfi.validbox();
	const int* lo = bx.loVect();
	const int* hi = bx.hiVect();

	BL_FORT_PROC_CALL(RNS_ADVCHEM, rns_advchem)
	    (lo, hi, BL_TO_FORTRAN(U[i]), dt);
    }

    post_update(U);
}


void
RNS::advance_AD(MultiFab& Unew, Real time, Real dt)
{
    MultiFab Uprime(grids,NUM_STATE,0);
    MultiFab U0(grids,NUM_STATE,0);
    MultiFab::Copy(U0, Unew, 0, 0, NUM_STATE, 0);

    if (RK_order == 2)
    {
	int finest_level = parent->finestLevel();

	//
	// Get pointers to Flux registers, or set pointer to zero if not there.
	//
	FluxRegister *fine    = 0;
	FluxRegister *current = 0;

	if (do_reflux && level < finest_level) {
	    fine = &getFluxReg(level+1);
	    fine->setVal(0.0);
	}
	if (do_reflux && level > 0) {
	    current = &getFluxReg(level);
	}

	// Step 1 of RK2
	dUdt_AD(Unew, Uprime, time, no_fill);
	update_rk(Unew, U0, 0.5*dt, Uprime); // Unew = U0 + 0.5*dt*Uprime
	post_update(Unew);

	// Step 2 of RK2
	int fill_boundary_type = use_FillCoarsePatch;
	dUdt_AD(Unew, Uprime, time+0.5*dt, fill_boundary_type, fine, current, dt);
	update_rk(Unew, U0, dt, Uprime); // Unew = U0 + dt*Uprime
	post_update(Unew);
    }
    else if (RK_order == 3)
    {
	BL_ASSERT(level==0);

	int fill_boundary_type = use_FillBoundary;
	MultiFab Utmp(grids,NUM_STATE,NUM_GROW);

	// Step 1 of RK3
	dUdt_AD(Unew, Uprime, time, no_fill);
	update_rk(Unew, U0, dt, Uprime);
	post_update(Unew);

	// Step 2 of RK3
	dUdt_AD(Unew, Uprime, time+(1./3.)*dt, fill_boundary_type);
	update_rk(Utmp, 0.75, U0, 0.25, Unew, 0.25*dt, Uprime);
	post_update(Utmp);

	// Step 3 of RK3
	dUdt_AD(Utmp, Uprime, time+(2./3.)*dt, fill_boundary_type);
	update_rk(Unew, 1./3., U0, 2./3., Utmp, (2./3.)*dt, Uprime);
	post_update(Unew);
    }
}

#ifdef USE_SDCLIB
BEGIN_EXTERN_C

//
// Compute dU_{AD}/dt.
//
// XXX: it might be interesting to track the magnitude of reflux
// registers.
//
void sdc_f1eval(void *F, void *Q, double t, sdc_state *state, void *ctx)
{
  RNS&      rns    = *((RNS*) ctx);
  MultiFab& U      = *((MultiFab*) Q);
  MultiFab& Uprime = *((MultiFab*) F);

  if (rns.verbose > 1 && ParallelDescriptor::IOProcessor()) {
    cout << "MLSDC evaluating adv/diff:  level: " << rns.Level() << ", node: " << state->node << endl;
  }

  rns.dUdt_AD(U, Uprime, t, RNS::use_FillBoundary, 0, 0, 0.0);
}

//
// Compute dU_R/dt.
//
// Note:
//
//   * Uprime doesn't have ghost cells
//
//   * Calling advance_chemistry with dt=0.0 puts dU_R/dt into tmp.
//
void sdc_f2eval(void *F, void *Q, double t, sdc_state *state, void *ctx)
{
  RNS&      rns    = *((RNS*) ctx);
  MultiFab& U      = *((MultiFab*) Q);
  MultiFab& Uprime = *((MultiFab*) F);

  Uprime.setVal(0.0);

  if (rns.chemSolve->isNull) return;

  if (rns.verbose > 1 && ParallelDescriptor::IOProcessor()) {
    cout << "MLSDC evaluating chemistry: level: " << rns.Level() << ", node: " << state->node << endl;
  }

  rns.fill_boundary(U, state->t, RNS::use_FillBoundary);
  BL_ASSERT(U.contains_nan() == false);
  rns.dUdt_chemistry(U, Uprime);
  BL_ASSERT(Uprime.contains_nan() == false);
}

//
// Solve U - dt dU_R/dt (U) = RHS for U.
//
// This advances chemistry from 0 to dt using RHS as the initial
// condition to obtain U.  Then, dU_R/dt is set to (U - RHS) / dt.
//
// XXX: it might be interesting to track the difference between Uprime
// as calculated above and calling f2eval...
//
void sdc_f2comp(void *F, void *Q, double t, double dt, void *RHS, sdc_state *state, void *ctx)
{
  RNS&      rns    = *((RNS*) ctx);
  MultiFab& U      = *((MultiFab*) Q);
  MultiFab& Uprime = *((MultiFab*) F);
  MultiFab& Urhs   = *((MultiFab*) RHS);

  BL_ASSERT(Urhs.contains_nan() == false);
  MultiFab::Copy(U, Urhs, 0, 0, U.nComp(), U.nGrow());

  if (rns.chemSolve->isNull) {
   Uprime.setVal(0.0);
    return;
  }

  if (rns.verbose > 1 && ParallelDescriptor::IOProcessor()) {
    cout << "MLSDC advancing  chemistry: level: " << rns.Level() << ", node: " << state->node << endl;
  }

  rns.fill_boundary(U, state->t, RNS::use_FillBoundary);

  BL_ASSERT(U.contains_nan() == false);
  rns.advance_chemistry(U, dt);

  Uprime.copy(U);
  Uprime.minus(Urhs, 0, Uprime.nComp(), 0);
  Uprime.mult(1./dt);
}

void sdc_poststep_hook(void *Q, sdc_state *state, void *ctx)
{
  RNS&      rns    = *((RNS*) ctx);
  MultiFab& U      = *((MultiFab*) Q);

  rns.post_update(U);
}

END_EXTERN_C
#endif
