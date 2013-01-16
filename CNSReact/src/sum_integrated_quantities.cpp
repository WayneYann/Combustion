#include <iomanip>

#include <CNSReact.H>
#include <CNSReact_F.H>

void
CNSReact::sum_integrated_quantities ()
{
    if (verbose <= 0) return;

    int finest_level = parent->finestLevel();
    Real dt_crse     = parent->dtLevel(0);
    Real time        = state[State_Type].curTime();
    Real mass        = 0.0;
    Real xmom        = 0.0;
    Real rho_E       = 0.0;
    Real com_xloc    = 0.0;
    Real com_xvel    = 0.0;
#if (BL_SPACEDIM>=2)
    Real ymom        = 0.0;
    Real com_yloc    = 0.0;
    Real com_yvel    = 0.0;
#endif
#if (BL_SPACEDIM==3)
    Real zmom        = 0.0;
    Real com_zloc    = 0.0;
    Real com_zvel    = 0.0;
#endif

    for (int lev = 0; lev <= finest_level; lev++)
    {
        CNSReact& cns_lev = getLevel(lev);

        mass     += cns_lev.volWgtSum("density", time);
        xmom     += cns_lev.volWgtSum("xmom", time);
#if (BL_SPACEDIM == 2)
       if (Geometry::IsRZ()) 
          xmom = 0.;
#endif

#if (BL_SPACEDIM>=2)
       ymom     += cns_lev.volWgtSum("ymom", time);
#endif
#if (BL_SPACEDIM==3)
       zmom     += cns_lev.volWgtSum("zmom", time);
#endif
        rho_E    += cns_lev.volWgtSum("rho_E", time);
    }

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        std::cout << '\n';
        std::cout << "TIME= " << time << " MASS        = "   << mass  << '\n';
        std::cout << "TIME= " << time << " XMOM        = "   << xmom     << '\n';
#if (BL_SPACEDIM>=2)
        std::cout << "TIME= " << time << " YMOM        = "   << ymom     << '\n';
#endif
#if (BL_SPACEDIM==3)
        std::cout << "TIME= " << time << " ZMOM        = "   << zmom     << '\n';
#endif
        std::cout << "TIME= " << time << " RHO*E       = "   << rho_E     << '\n';
	std::cout<<'\n';
    }
}

Real
CNSReact::volWgtSum (const std::string& name,
                   Real               time)
{
    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf      = derive(name,time,0);

    BL_ASSERT(mf != 0);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }
        Real s;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();
#if(BL_SPACEDIM < 3) 
        const Real* rad = radius[mfi.index()].dataPtr();
        int irlo        = lo[0]-radius_grow;
        int irhi        = hi[0]+radius_grow;
#endif

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //
#if(BL_SPACEDIM == 1) 
	BL_FORT_PROC_CALL(CNS_SUMMASS,cns_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s,rad,irlo,irhi);
#elif(BL_SPACEDIM == 2)
	BL_FORT_PROC_CALL(CNS_SUMMASS,cns_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s,rad,irlo,irhi);
#elif(BL_SPACEDIM == 3)
	BL_FORT_PROC_CALL(CNS_SUMMASS,cns_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s);
#endif
        sum += s;
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}
