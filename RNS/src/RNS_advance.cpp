#include <winstd.H>

#include "RNS.H"
#include "RNS_F.H"

#ifdef USE_SDCLIB
#include "SDCAmr.H"
#endif

// #include <ArrayView.H>

using std::string;

Real
RNS::advance (Real time,
	      Real dt,
	      int  iteration,
	      int  ncycle)
{
    if (level == 0) {
	for (int lev=0; lev<parent->finestLevel(); lev++) {
	    RNS& rns = *dynamic_cast<RNS*>(&getLevel(lev));
	    rns.zeroChemStatus();
	}
    }

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

    if (! ChemDriver::isNull())
    {
	// do half-dt chemistry
	advance_chemistry(Unew, 0.5*dt);

	// After chemistry, fill ghosts cells with the post-chemistry state of current level
	fill_boundary(Unew, time, use_FillBoundary);
    }

    // Advance Advection & Diffusion
    advance_AD(Unew, time, dt, iteration, ncycle);

    if (! ChemDriver::isNull())
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
RNS::fill_boundary(MultiFab& U, Real time, int type_in, bool isCorrection, bool isFEval)
{
    BL_PROFILE("RNS::fill_boundary()");

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
	BoxArray grids_g(grids);
	for (int ibox=0; ibox<grids_g.size(); ibox++)
	{
	    const Box b = BoxLib::grow(grids_g[ibox], NUM_GROW);
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

	// no break; so it will go to next case and set physical boundaries

    case set_PhysBoundary:

	geom.FillPeriodicBoundary(U, true);

	if (isCorrection) break;

	if (isFEval) BL_FORT_PROC_CALL(SET_SDC_BOUNDARY_FLAG, set_sdc_boundary_flag)();

	for (MFIter mfi(U); mfi.isValid(); ++mfi)
	{
	    setPhysBoundaryValues(U[mfi],
				  State_Type,
				  time,
				  0,
				  0,
				  NUM_STATE);
	}

	if (isFEval) BL_FORT_PROC_CALL(UNSET_SDC_BOUNDARY_FLAG, unset_sdc_boundary_flag)();

	break;
    }
}


#ifndef USE_SDCLIB
void
RNS::fill_rk_boundary(MultiFab& U, Real time, Real dt, int stage, int iteration, int ncycle)
{
    // all boundaries must be periodic!

    BL_ASSERT(level > 0);
    BL_ASSERT(RK_order > 2);

    static MultiFab* U0;
    static MultiFab* k1;
    static MultiFab* k2;
    static MultiFab* k3;
    static MultiFab* k4;

    const int ncomp = U.nComp();
    const int ngrow = U.nGrow();

    if (iteration == 1 && stage == 0) {

	U0 = new MultiFab(grids, ncomp, ngrow);
	k1 = new MultiFab(grids, ncomp, ngrow);
	k2 = new MultiFab(grids, ncomp, ngrow);
	k3 = new MultiFab(grids, ncomp, ngrow);
	k4 = (RK_order == 4) ? new MultiFab(grids, ncomp, ngrow) : 0;

	U0->setVal(0.0);
	k1->setVal(0.0);
	k2->setVal(0.0);
	k3->setVal(0.0);
	if (k4) k4->setVal(0.0);

	RNS& levelG = *dynamic_cast<RNS*>(&getLevel(level-1));

	const IntVect         ratio = levelG.fineRatio();
	const DescriptorList& dl    = get_desc_lst();
	const Array<BCRec>&   bcs   = dl[0].getBCs();
	Interpolater&         map   = *dl[0].interp();
	Array<BCRec>          bcr(ncomp);
	
	const Geometry& geomG = levelG.Geom();
	
	// make a coarse version of the fine multifab
	BoxArray ba_C(U.size());
	for (int i=0; i<ba_C.size(); i++) {
	    ba_C.set(i, map.CoarseBox(U.fabbox(i), ratio));
	}
	
	MultiFab UC(ba_C, ncomp, 0);
	
	bool touch = false;
	bool touch_periodic = false;
	const Box& crse_domain_box = levelG.Domain();
	if (geomG.isAnyPeriodic()) {
	    for (int i=0; i < ba_C.size(); i++) {
		if (! crse_domain_box.contains(ba_C[i])) {
		    touch = true;
		    for (int idim=0; idim<BL_SPACEDIM; idim++) {
			if (geomG.isPeriodic(i)   ||
			    ba_C[i].bigEnd(idim) > crse_domain_box.bigEnd(idim) ||
			    ba_C[i].smallEnd(idim) < crse_domain_box.smallEnd(idim) )
			{
			    touch_periodic = true;
			    break;
			}
		    }
		}
		if (touch_periodic) break;
	    }
	}
	else {
	    for (int i=0; i < ba_C.size(); i++) {
		if (! crse_domain_box.contains(ba_C[i])) {
	      touch = true;
	      break;
		}
	    }
	}

	MultiFab* UG;
	MultiFab* UF;
	
	for (int ii=0; ii<RK_order+1; ii++) {
	    switch (ii) {
	    case 0:
		UG = &(levelG.get_old_data(0));
		UF = U0;
		break;
	    case 1:
		UG = &(levelG.getRKk(0));
		UF = k1;
		break;
	    case 2:
		UG = &(levelG.getRKk(1));
		UF = k2;
		break;
	    case 3:
		UG = &(levelG.getRKk(2));
		UF = k3;
		break;
	    case 4:
		UG = &(levelG.getRKk(3));
		UF = k4;
		break;
	    }

	    if (!touch) {
		// If level F does not touch physical boundaries, then AMR levels are
		// properly nested so that the valid rgions of UC are contained inside
		// the valid regions of UG.  So FabArray::copy is all we need.
		UC.copy(*UG);
	    }
	    else if (touch_periodic) {
		// This case is more complicated because level F might touch only one
		// of the periodic boundaries.
		int ng_C = 0;
		{
		    Box box_C = ba_C.minimalBox();
		    for (int idim=0; idim < BL_SPACEDIM; idim++) {
			int gap_hi = box_C.bigEnd(idim) - crse_domain_box.bigEnd(idim);
			int gap_lo = crse_domain_box.smallEnd(idim) - box_C.smallEnd(idim);
			ng_C = std::max(ng_C, gap_hi);
			ng_C = std::max(ng_C, gap_lo);
		    }
		}
		int ng_G = UG->nGrow();
		const BoxArray& ba_G = UG->boxArray();
		
		MultiFab* UG_safe;
		MultiFab UGG;
		
		if (ng_C > ng_G) {
		    UGG.define(ba_G, ncomp, ng_C, Fab_allocate);
		    MultiFab::Copy(UGG, *UG, 0, 0, ncomp, 0);
		    UG_safe = &UGG;
		}
		else {
		    UG_safe = UG;
		}
	    
		levelG.fill_boundary(*UG_safe, time, RNS::use_FillBoundary);

		// We cannot do FabArray::copy() directly on UG because it copies only form
		// valid regions.  So we need to make the ghost cells of UG valid.
		BoxArray ba_G2(UG->size());
		for (int i=0; i<ba_G2.size(); i++) {
		    ba_G2.set(i, BoxLib::grow(ba_G[i],ng_C));
		}
		
		MultiFab UG2(ba_G2, ncomp, 0);
		for (MFIter mfi(UG2); mfi.isValid(); ++mfi)
		{
		    int i = mfi.index();
		    UG2[i].copy((*UG_safe)[i]);  // Fab to Fab copy
		}
		
		UC.copy(UG2);
	    }
	    else {
		UC.copy(*UG);
		levelG.fill_boundary(UC, time, RNS::set_PhysBoundary);
	    }

	    // now that UF is completely contained within UC, cycle through each
	    // FAB in UF and interpolate from the corresponding FAB in UC
	    for (MFIter mfi(*UF); mfi.isValid(); ++mfi) {
		BoxLib::setBC((*UF)[mfi].box(), Domain(), 0, 0, ncomp, bcs, bcr);
		Geometry fine_geom((*UF)[mfi].box());
		Geometry crse_geom(UC[mfi].box());
	    
		map.interp(UC[mfi], 0, (*UF)[mfi], 0, ncomp, (*UF)[mfi].box(), ratio,
			   crse_geom, fine_geom, bcr, 0, 0);
	    }
	}
    }

    Real dtdt = 1.0/ncycle;
    Real xsi0 = (iteration-1.0)*dtdt;

    for (MFIter mfi(U); mfi.isValid(); ++mfi) 
    {
	int i = mfi.index();

	const Box& vbox = grids[i];
	const Box& gbox = U[i].box();
	const BoxArray& ba = BoxLib::boxComplement(gbox, vbox);

	for (int ibox=0; ibox<ba.size(); ibox++)
	{
	    if (RK_order == 3) {
		BL_FORT_PROC_CALL(RNS_FILL_RK3_BNDRY, rns_fill_rk3_bndry)
		    (ba[ibox].loVect(), ba[ibox].hiVect(),
		     BL_TO_FORTRAN(U[i]),
		     BL_TO_FORTRAN((*U0)[i]),
		     BL_TO_FORTRAN((*k1)[i]),
		     BL_TO_FORTRAN((*k2)[i]),
		     BL_TO_FORTRAN((*k3)[i]),
		     dtdt, xsi0, stage);
	    }
	    else {
		BL_FORT_PROC_CALL(RNS_FILL_RK4_BNDRY, rns_fill_rk4_bndry)
		    (ba[ibox].loVect(), ba[ibox].hiVect(),
		     BL_TO_FORTRAN(U[i]),
		     BL_TO_FORTRAN((*U0)[i]),
		     BL_TO_FORTRAN((*k1)[i]),
		     BL_TO_FORTRAN((*k2)[i]),
		     BL_TO_FORTRAN((*k3)[i]),
		     BL_TO_FORTRAN((*k4)[i]),
		     dtdt, xsi0, stage);
	    }
	}
    }

    fill_boundary(U, time, RNS::use_FillBoundary);

    if (iteration == ncycle && stage+1 == RK_order) {
	delete U0;
	delete k1;
	delete k2;
	delete k3;
	if (k4) delete k4;
    }
}
#endif


void
RNS::dUdt_AD(MultiFab& U, MultiFab& Uprime, Real time, int fill_boundary_type,
	     FluxRegister* fine, FluxRegister* current, Real dt, bool partialUpdate)
{
    BL_PROFILE("RNS::dUdt_AD()");

    if (partialUpdate && touchFine.empty()) buildTouchFine();

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

#ifndef NDEBUG
    BL_ASSERT(U.contains_nan() == false);
#endif

    const Real *dx = geom.CellSize();

    int iteration=-1;
    BL_FORT_PROC_CALL(RNS_PASSINFO,rns_passinfo)(level,iteration,time);

    for (MFIter mfi(Uprime); mfi.isValid(); ++mfi)
    {
	int i = mfi.index();

	if (partialUpdate && !touchFine[i]) {
	    if (do_reflux && fine) {
		for (int idim = 0; idim < BL_SPACEDIM ; idim++) {
		    fluxes[idim][i].setVal(0.0);
		}
	    }
	    continue;
	}

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
RNS::dUdt_chemistry(const MultiFab& U, MultiFab& Uprime, bool partialUpdate)
{
    BL_PROFILE("RNS::dUdt_chemistry()");

    BL_ASSERT( ! ChemDriver::isNull() );

    if (partialUpdate && touchFine.empty()) buildTouchFine();

    int iteration=-1;
    Real time=-1.0;
    BL_FORT_PROC_CALL(RNS_PASSINFO,rns_passinfo)(level,iteration,time);

    for (MFIter mfi(U); mfi.isValid(); ++mfi)
    {
	const int i = mfi.index();

	if (partialUpdate && !touchFine[i]) continue;

	const Box& bx = mfi.validbox();
	const int* lo = bx.loVect();
	const int* hi = bx.hiVect();

	BL_FORT_PROC_CALL(RNS_DUDT_CHEM, rns_dudt_chem)
	    (lo, hi, BL_TO_FORTRAN(U[i]), BL_TO_FORTRAN(Uprime[i]),
	     BL_TO_FORTRAN((*chemstatus)[i]));
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
    BL_PROFILE("RNS::advance_chemistry()");

    BL_ASSERT( ! ChemDriver::isNull() );

    int iteration=-1;
    Real time=-1.;
    BL_FORT_PROC_CALL(RNS_PASSINFO,rns_passinfo)(level,iteration,time);

    for (MFIter mfi(U); mfi.isValid(); ++mfi)
    {
	const int   i = mfi.index();
	const Box& bx = mfi.validbox();
	const int* lo = bx.loVect();
	const int* hi = bx.hiVect();

	BL_FORT_PROC_CALL(RNS_ADVCHEM, rns_advchem)
	    (lo, hi, BL_TO_FORTRAN(U[i]), BL_TO_FORTRAN((*chemstatus)[i]), dt);
    }

    post_update(U);
}


void
RNS::advance_chemistry(MultiFab& U, const MultiFab& Uguess, Real dt)
{
    BL_PROFILE("RNS::advance_chemistry() w/ guess");

    BL_ASSERT( ! ChemDriver::isNull() );
    BL_ASSERT( Uguess.nGrow() == 0 );

    int iteration=-1;
    Real time=-1.0;
    BL_FORT_PROC_CALL(RNS_PASSINFO,rns_passinfo)(level,iteration,time);

    for (MFIter mfi(U); mfi.isValid(); ++mfi)
    {
	const int   i = mfi.index();
	const Box& bx = mfi.validbox();
	const int* lo = bx.loVect();
	const int* hi = bx.hiVect();

	BL_FORT_PROC_CALL(RNS_ADVCHEM2, rns_advchem2)
	    (lo, hi, BL_TO_FORTRAN(U[i]), BL_TO_FORTRAN((*chemstatus)[i]),
	     BL_TO_FORTRAN(Uguess[i]), dt);
    }

    post_update(U);
}


void
RNS::advance_AD(MultiFab& Unew, Real time, Real dt, int iteration, int ncycle)
{
    MultiFab U0(grids,NUM_STATE,0);
    MultiFab::Copy(U0, Unew, 0, 0, NUM_STATE, 0);

    int finest_level = parent->finestLevel();

    FluxRegister *fine    = 0;
    FluxRegister *current = 0;
    FluxRegister *fine_RK = 0;

    if (do_reflux && level < finest_level) {
	fine = &getFluxReg(level+1);
	fine->setVal(0.0);
	if (RK_order > 2) {
	    fine_RK = &getFluxRegRK(level+1);
	}
    }
    if (do_reflux && level > 0) {
	current = &getFluxReg(level);
    }

    if (RK_order == 2)
    {
	MultiFab Uprime(grids,NUM_STATE,0);

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
#ifndef USE_SDCLIB
    else if (RK_order == 3)
    {
	int fill_boundary_type, stage;

	// Step 1 of RK3
	stage = 0;
	if (level > 0) {
	    fill_rk_boundary(Unew, time, dt, stage, iteration, ncycle);
	}
	if (fine_RK) {
	    fine_RK->setVal(0.0);
	}
	dUdt_AD(Unew, RK_k[stage], time, no_fill, fine_RK, current, dt/6.);
	if (fine_RK) {	
	    (*fine) += (*fine_RK);
	}
	RK_k[stage].mult(dt);
	update_rk(Unew, U0, 1.0, RK_k[stage]);
	post_update(Unew);

	// Step 2 of RK3
	stage = 1;
	if (level > 0) {
	    fill_rk_boundary(Unew, time+dt, dt, stage, iteration, ncycle);
	    fill_boundary_type = no_fill;
	}
	else {
	    fill_boundary_type = use_FillBoundary;
	}
	if (fine_RK) {
	    fine_RK->setVal(0.0);
	}
	dUdt_AD(Unew, RK_k[stage], time+dt, fill_boundary_type, fine_RK, current, dt/6.);
	if (fine_RK) {
	    (*fine) += (*fine_RK);
	}
	RK_k[stage].mult(dt);
	for (MFIter mfi(Unew); mfi.isValid(); ++mfi)
	{
	    const int i = mfi.index();
	    const Box& bx = mfi.validbox();
	    
	    Unew[i].copy(U0[i], bx);
	    Unew[i].saxpy(0.25, RK_k[0][i]);
	    Unew[i].saxpy(0.25, RK_k[1][i]);
	}
	post_update(Unew);

	// Step 3 of RK3
	stage = 2;
	if (level > 0) {
	    fill_rk_boundary(Unew, time+0.5*dt, dt, stage, iteration, ncycle);
	    fill_boundary_type = no_fill;
	}
	else {
	    fill_boundary_type = use_FillBoundary;
	}	
	if (fine_RK) {
	    fine_RK->setVal(0.0);
	}
	dUdt_AD(Unew, RK_k[stage], time+0.5*dt, fill_boundary_type, fine_RK, current, dt*(2./3.));
	if (fine_RK) {
	    (*fine) += (*fine_RK);
	}
	RK_k[stage].mult(dt);
	for (MFIter mfi(Unew); mfi.isValid(); ++mfi)
	{
	    const int i = mfi.index();
	    const Box& bx = mfi.validbox();
	    
	    Unew[i].copy(U0[i], bx);
	    Unew[i].saxpy(1./6., RK_k[0][i]);
	    Unew[i].saxpy(1./6., RK_k[1][i]);
	    Unew[i].saxpy(2./3., RK_k[2][i]);
	}
	post_update(Unew);
    }
    else if (RK_order == 4)
    {
	int fill_boundary_type, stage;

	// Step 1 of RK4
	stage = 0;
	if (level > 0) {
	    fill_rk_boundary(Unew, time, dt, stage, iteration, ncycle);
	}
	if (fine_RK) {
	    fine_RK->setVal(0.0);
	}
	dUdt_AD(Unew, RK_k[stage], time, no_fill, fine_RK, current, dt/6.);
	if (fine_RK) {	
	    (*fine) += (*fine_RK);
	}
	RK_k[stage].mult(dt);
	update_rk(Unew, U0, 0.5, RK_k[stage]);
	post_update(Unew);

	// Step 2 of RK4
	stage = 1;
	if (level > 0) {
	    fill_rk_boundary(Unew, time+0.5*dt, dt, stage, iteration, ncycle);
	    fill_boundary_type = no_fill;
	}
	else {
	    fill_boundary_type = use_FillBoundary;
	}
	if (fine_RK) {
	    fine_RK->setVal(0.0);
	}
	dUdt_AD(Unew, RK_k[stage], time+0.5*dt, fill_boundary_type, fine_RK, current, dt/3.);
	if (fine_RK) {
	    (*fine) += (*fine_RK);
	}
	RK_k[stage].mult(dt);
	update_rk(Unew, U0, 0.5, RK_k[stage]);
	post_update(Unew);

	// Step 3 of RK4
	stage = 2;
	if (level > 0) {
	    fill_rk_boundary(Unew, time+0.5*dt, dt, stage, iteration, ncycle);
	    fill_boundary_type = no_fill;
	}
	else {
	    fill_boundary_type = use_FillBoundary;
	}	
	if (fine_RK) {
	    fine_RK->setVal(0.0);
	}
	dUdt_AD(Unew, RK_k[stage], time+0.5*dt, fill_boundary_type, fine_RK, current, dt/3.);
	if (fine_RK) {
	    (*fine) += (*fine_RK);
	}
	RK_k[stage].mult(dt);
	update_rk(Unew, U0, 1.0, RK_k[stage]);
	post_update(Unew);

	// Step 4 of RK4
	stage = 3;
	if (level > 0) {
	    fill_rk_boundary(Unew, time+dt, dt, stage, iteration, ncycle);
	    fill_boundary_type = no_fill;
	}
	else {
	    fill_boundary_type = use_FillBoundary;
	}	
	if (fine_RK) {
	    fine_RK->setVal(0.0);
	}
	dUdt_AD(Unew, RK_k[stage], time+dt, fill_boundary_type, fine_RK, current, dt/6.);
	if (fine_RK) {
	    (*fine) += (*fine_RK);
	}
	RK_k[stage].mult(dt);
	for (MFIter mfi(Unew); mfi.isValid(); ++mfi)
	{
	    const int i = mfi.index();
	    const Box& bx = mfi.validbox();

	    Unew[i].copy(U0[i], bx);
	    Unew[i].saxpy(1./6., RK_k[0][i]);
	    Unew[i].saxpy(1./3., RK_k[1][i]);
	    Unew[i].saxpy(1./3., RK_k[2][i]);
	    Unew[i].saxpy(1./6., RK_k[3][i]);
	}
	post_update(Unew);	
    }
#endif
}

#ifdef USE_SDCLIB
BEGIN_EXTERN_C

//
// Compute dU_{AD}/dt.
//
void sdc_f1eval(void *Fp, void *Qp, double t, sdc_state *state, void *ctx)
{
  RNS&      rns    = *((RNS*) ctx);
  RNSEncap& Q      = *((RNSEncap*) Qp);
  RNSEncap& F      = *((RNSEncap*) Fp);
  MultiFab& U      = *Q.U;
  MultiFab& Uprime = *F.U;

  if (rns.verbose > 1 && ParallelDescriptor::IOProcessor()) {
    cout << "MLSDC evaluating adv/diff:"
	 << "  level: " << rns.Level() << ", node: " << state->node << endl;
  }

  bool partialUpdate;
  if (state->flags & SDC_POSTRESTRICT) {
      partialUpdate = true;
  }
  else {
      partialUpdate = false;
  }

  if (F.fine_flux != NULL) {
    F.fine_flux->setVal(0.0);
    partialUpdate = false;  // partial update doesn't work when there is fine flux register
                            // (i.e., this is a fine level)
  }

  rns.dUdt_AD(U, Uprime, t, RNS::use_FillBoundary, F.crse_flux, F.fine_flux, 1.0,
	      partialUpdate);

#ifndef NDEBUG
  BL_ASSERT(U.contains_nan() == false);
  BL_ASSERT(Uprime.contains_nan() == false);
#endif
}

//
// Compute dU_R/dt.
//
// Note:
//
//   * Uprime doesn't have ghost cells
//
void sdc_f2eval(void *Fp, void *Qp, double t, sdc_state *state, void *ctx)
{
  RNS&      rns    = *((RNS*) ctx);
  RNSEncap& Q      = *((RNSEncap*) Qp);
  RNSEncap& F      = *((RNSEncap*) Fp);
  MultiFab& U      = *Q.U;
  MultiFab& Uprime = *F.U;
  Real dt = state->dt;

  if (ChemDriver::isNull() || !RNS::do_chemistry) {
      Uprime.setVal(0.0);
      return;
  }

  if (rns.verbose > 1 && ParallelDescriptor::IOProcessor()) {
    cout << "MLSDC evaluating chemistry:"
	 << " level: " << rns.Level() << ", node: " << state->node
	 << ", dt = " << dt << endl;
  }

  BL_FORT_PROC_CALL(RNS_PASSINFO, rns_passinfo)(rns.Level(),state->iter,t);

  rns.fill_boundary(U, state->t, RNS::use_FillBoundary);
  BL_ASSERT(U.contains_nan() == false);

  bool partialUpdate;
  if (state->flags & SDC_POSTRESTRICT) {
      partialUpdate = true;
  }
  else {
      partialUpdate = false;
  }

  rns.dUdt_chemistry(U, Uprime, partialUpdate);
  BL_ASSERT(Uprime.contains_nan() == false);
}

//
// Solve U - dt dU_R/dt (U) = RHS for U.
//
// This advances chemistry from 0 to dt using RHS as the initial
// condition to obtain U.  Then, dU_R/dt is either set to (U - RHS) / dt
// or computed by calling dUdt_chmistry.
//
void sdc_f2comp(void *Fp, void *Qp, double t, double dt, void *RHSp, sdc_state *state, void *ctx)
{
  static int first = 1;
  RNS&      rns    = *((RNS*) ctx);
  RNSEncap& Q      = *((RNSEncap*) Qp);
  RNSEncap& F      = *((RNSEncap*) Fp);
  RNSEncap& RHS    = *((RNSEncap*) RHSp);
  MultiFab& U      = *Q.U;
  MultiFab& Uprime = *F.U;
  MultiFab& Urhs   = *RHS.U;

  BL_ASSERT(Urhs.contains_nan() == false);

  if (ChemDriver::isNull() || !RNS::do_chemistry) {
      MultiFab::Copy(U, Urhs, 0, 0, U.nComp(), 0);
      Uprime.setVal(0.0);
      return;
  }

  BL_FORT_PROC_CALL(RNS_PASSINFO, rns_passinfo)(rns.Level(),state->iter,t);

  if (first) {
      const SDCAmr* sdcamr = rns.getSDCAmr();
      if (rns.check_imex_order(sdcamr->ho_imex)) {
	  BoxLib::Warning("\n*** ho_imex is incompatible with chem_solver");
      }
      first = 0;
  }

  if (rns.verbose > 1 && ParallelDescriptor::IOProcessor()) {
      cout << "MLSDC advancing  chemistry:"
	   << " level: " << rns.Level() << ", node: " << state->node << endl;
  }

  MultiFab Uguess;
  bool Uguess_defined=false;
  if ( rns.f2comp_timer[state->node] >= rns.f2comp_nbdf ) {
      Uguess.define(U.boxArray(), U.nComp(), 0, Fab_allocate);
      MultiFab::Copy(Uguess, U, 0, 0, U.nComp(), 0);
      Uguess_defined = true;
  }

  MultiFab::Copy(U, Urhs, 0, 0, U.nComp(), 0);

  rns.fill_boundary(U, state->t, RNS::use_FillBoundary);

  BL_ASSERT(U.contains_nan() == false);

  if (Uguess_defined) {
      rns.advance_chemistry(U, Uguess, dt);
  }
  else {
      rns.advance_chemistry(U, dt);
  }

  if (rns.f2comp_simple_dUdt) {
      Uprime.copy(U);
      Uprime.minus(Urhs, 0, Uprime.nComp(), 0);
      Uprime.mult(1./dt);
  }
  else {
      rns.fill_boundary(U, state->t+dt, RNS::use_FillBoundary);
      rns.dUdt_chemistry(U, Uprime);
  }

  rns.f2comp_timer[state->node]++;
}

void sdc_poststep_hook(void *Qp, sdc_state *state, void *ctx)
{
  RNS&      rns = *((RNS*) ctx);
  RNSEncap& Q   = *((RNSEncap*) Qp);
  MultiFab& U   = *Q.U;

  rns.post_update(U);
}

END_EXTERN_C
#endif
