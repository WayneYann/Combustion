// #include <winstd.H>

// #include "CNSReact.H"
// #include "CNSReact_F.H"
// #include "Diffusion.H"

// using std::string;

// void
// CNSReact::add_diffusion_to_old_source (MultiFab& ext_src_old, MultiFab& OldTempDiffTerm, Real prev_time)
// {
//     // Define an explicit temperature update.
//     OldTempDiffTerm.setVal(0.);
//     if (diffuse_temp == 1) {
//       getTempDiffusionTerm(prev_time,OldTempDiffTerm);
//       MultiFab::Add(ext_src_old,OldTempDiffTerm,0,Eden,1,0);
//       MultiFab::Add(ext_src_old,OldTempDiffTerm,0,Eint,1,0);
//     }
//     geom.FillPeriodicBoundary(ext_src_old,0,NUM_STATE);
// }

// void
// CNSReact::time_center_diffusion(MultiFab& S_new, MultiFab& OldTempDiffTerm, Real cur_time, Real dt)
// {
//         // Correct the temperature update so that it will be time-centered.
//         MultiFab NewTempDiffTerm(grids,1,1);
//         NewTempDiffTerm.setVal(0.);
//         if (diffuse_temp == 1) {
//            getTempDiffusionTerm(cur_time,NewTempDiffTerm);
//            NewTempDiffTerm.mult( 0.5*dt);
//            OldTempDiffTerm.mult(-0.5*dt);
//            // Subtract off half of the old source term, and add half of the new.
//            MultiFab::Add(S_new,OldTempDiffTerm,0,Eden,1,0);
//            MultiFab::Add(S_new,OldTempDiffTerm,0,Eint,1,0);
//            MultiFab::Add(S_new,NewTempDiffTerm,0,Eden,1,0);
//            MultiFab::Add(S_new,NewTempDiffTerm,0,Eint,1,0);
//            computeTemp(S_new);
//         }
// }

// void
// CNSReact::full_diffusion_update (MultiFab& S_new, Real prev_time, Real cur_time, Real dt)
// {
//         if (diffuse_temp == 1) {
//            // Define an explicit temperature update.
//            MultiFab OldTempDiffTerm(grids,1,1);
//            OldTempDiffTerm.setVal(0.);
//            getTempDiffusionTerm(prev_time,OldTempDiffTerm);
//            OldTempDiffTerm.mult(dt);
//            MultiFab::Add(S_new,OldTempDiffTerm,0,Eden,1,0);
//            MultiFab::Add(S_new,OldTempDiffTerm,0,Eint,1,0);
//            computeTemp(S_new);

//            // Correct the temperature update so that it will be time-centered.
//            MultiFab NewTempDiffTerm(grids,1,1);
//            NewTempDiffTerm.setVal(0.);
//            getTempDiffusionTerm(cur_time,NewTempDiffTerm);
//            NewTempDiffTerm.mult( 0.5*dt);
//            OldTempDiffTerm.mult(-0.5*dt);
//            // Subtract off half of the old source term, and add half of the new.
//            MultiFab::Add(S_new,OldTempDiffTerm,0,Eden,1,0);
//            MultiFab::Add(S_new,OldTempDiffTerm,0,Eint,1,0);
//            MultiFab::Add(S_new,NewTempDiffTerm,0,Eden,1,0);
//            MultiFab::Add(S_new,NewTempDiffTerm,0,Eint,1,0);
//            computeTemp(S_new);
//         }
// }
