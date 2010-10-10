#include <winstd.H>

#include "CNSReact.H"
#include "CNSReact_F.H"

using std::string;

#ifdef REACTIONS
void
CNSReact::react_first_half_dt(MultiFab& S_old, Real time, Real dt) 
{
    // Make sure to zero these even if do_react == 0.
    MultiFab& ReactMF_old = get_old_data(Reactions_Type);
    ReactMF_old.setVal(0.);
    MultiFab& ReactMF = get_new_data(Reactions_Type);
    ReactMF.setVal(0.);
    if (do_react == 1)
    {
        strang_chem(S_old,ReactMF,time,dt);
        reset_internal_energy(S_old);
    }
}

void
CNSReact::react_second_half_dt(MultiFab& S_new, Real cur_time, Real dt) 
{
    if (do_react == 1) 
    {
        MultiFab& ReactMF = get_new_data(Reactions_Type);

        strang_chem(S_new,ReactMF,cur_time,dt);
        ReactMF.mult(1.0/dt);
        reset_internal_energy(S_new);
    }
}

void
CNSReact::strang_chem (MultiFab&  state,
                     MultiFab&  React_mf,
                     Real       time,
                     Real       dt)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::strang_chem(MultiFab&,...");
    const Real strt_time = ParallelDescriptor::second();

    for (MFIter Smfi(state); Smfi.isValid(); ++Smfi)
    {
        FArrayBox& fb   = state[Smfi];
        const Box& bx   = Smfi.validbox();
        reactState(fb, fb, React_mf[Smfi], bx, time, 0.5*dt);
    }

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

//      if (ParallelDescriptor::IOProcessor()) 
//          std::cout << "strang_chem time = " << run_time << '\n';
    }
}

void
CNSReact::reactState(FArrayBox&        Snew,
                   FArrayBox&        Sold,
                   FArrayBox&        ReactionTerms,
                   const Box&        box,
                   Real              time,
                   Real              dt_react)
{
    BL_FORT_PROC_CALL(CA_REACT_STATE,ca_react_state)
                     (box.loVect(), box.hiVect(), 
                     BL_TO_FORTRAN(Sold),
                     BL_TO_FORTRAN(Snew),
                     BL_TO_FORTRAN(ReactionTerms),
                     time,dt_react);
}
#endif
