#include <winstd.H>

#ifndef WIN32
#include <unistd.h>
#endif

#include <iomanip>
#include <iostream>
#include <string>
#include <ctime>

#include <Utility.H>
#include <ParmParse.H>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "RNS.H"
#include "RNS_F.H"

#include "buildInfo.H"

using std::string;

// I/O routines for RNS

void
RNS::restart (Amr&     papa,
	      istream& is,
	      bool     bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

    buildMetrics();
    
    BL_ASSERT(flux_reg == 0);
    if (level > 0 && do_reflux) 
    {
	flux_reg = new FluxRegister(grids,crse_ratio,level,NUM_STATE);
    }

    BL_ASSERT(chemstatus == 0);
    if (! ChemDriver::isNull()) {
      chemstatus = new MultiFab(grids,1,1);
      chemstatus->setVal(0.0,1);
    }

    BL_ASSERT(RK_k == 0);
    BL_ASSERT(flux_reg_RK == 0);
#ifndef USE_SDCLIB
    if (RK_order > 2) {
	RK_k = new MultiFab[RK_order];
	for (int i=0; i<RK_order; i++) {
	    RK_k[i].define(grids,NUM_STATE,0,Fab_allocate);
	}
	if (flux_reg) 
	{
	    flux_reg_RK = new FluxRegister(grids,crse_ratio,level,NUM_STATE);	    
	}
    }
#endif

    // get the elapsed CPU time to now;
    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
	// get ellapsed CPU time
	std::ifstream CPUFile;
	std::string FullPathCPUFile = parent->theRestartFile();
	FullPathCPUFile += "/CPUtime";
	CPUFile.open(FullPathCPUFile.c_str(), std::ios::in);	
  
	CPUFile >> previousCPUTimeUsed;
	CPUFile.close();

	std::cout << "read CPU time: " << previousCPUTimeUsed << "\n";
    }

    if (level == 0)
    {
	// get problem-specific stuff -- note all processors do this,
	// eliminating the need for a broadcast
	std::string dir = parent->theRestartFile();
	
	char * dir_for_pass = new char[dir.size() + 1];
	std::copy(dir.begin(), dir.end(), dir_for_pass);
	dir_for_pass[dir.size()] = '\0';
	
	int len = dir.size();
	
	Array<int> int_dir_name(len);
	for (int j = 0; j < len; j++)
	    int_dir_name[j] = (int) dir_for_pass[j];
	
	BL_FORT_PROC_CALL(PROBLEM_RESTART, problem_restart)(int_dir_name.dataPtr(), &len);
	
	delete [] dir_for_pass;
      }
}

void
RNS::checkPoint(const std::string& dir,
		std::ostream&  os,
		VisMF::How     how,
		bool dump_old_default)
{
    AmrLevel::checkPoint(dir, os, how, dump_old);

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
	// store ellapsed CPU time
	std::ofstream CPUFile;
	std::string FullPathCPUFile = dir;
	FullPathCPUFile += "/CPUtime";
	CPUFile.open(FullPathCPUFile.c_str(), std::ios::out);	
  
	CPUFile << std::setprecision(15) << getCPUTime();
	CPUFile.close();

	// store any problem-specific stuff
	char * dir_for_pass = new char[dir.size() + 1];
	std::copy(dir.begin(), dir.end(), dir_for_pass);
	dir_for_pass[dir.size()] = '\0';
	
	int len = dir.size();
	
	Array<int> int_dir_name(len);
	for (int j = 0; j < len; j++)
	    int_dir_name[j] = (int) dir_for_pass[j];

	BL_FORT_PROC_CALL(PROBLEM_CHECKPOINT, problem_checkpoint)(int_dir_name.dataPtr(), &len); 

	delete [] dir_for_pass;
    }
}

std::string
RNS::thePlotFileType () const
{
    //
    // Increment this whenever the writePlotFile() format changes.
    //
    static const std::string the_plot_file_type("HyperCLaw-V1.1");
    
    return the_plot_file_type;
}

void
RNS::setPlotVariables ()
{
    const Array<std::string>& spec_names = chemSolve->speciesNames();
    int inext = 0;

    if (plot_cons)
    {
	icomp_cons = 0;
	inext += NUM_STATE;
	plot_names.push_back("<rho>");
	plot_names.push_back("<xmom>");
#if (BL_SPACEDIM >= 2)
	plot_names.push_back("<ymom>");
#endif
#if (BL_SPACEDIM == 3)
	plot_names.push_back("<zmom>");
#endif
	plot_names.push_back("<rhoE>");
	plot_names.push_back("<T>");
	for (int i=0; i<NumSpec; i++)
	{
	    plot_names.push_back("<rho.Y(" + spec_names[i] + ")>");
	}
    }

    if (plot_prim)
    {
	icomp_prim = inext;
	inext += NUM_STATE;
	plot_names.push_back("density");
	plot_names.push_back("x_vel");
#if (BL_SPACEDIM >= 2)
	plot_names.push_back("y_vel");
#endif
#if (BL_SPACEDIM == 3)
	plot_names.push_back("z_vel");
#endif
	plot_names.push_back("pressure");
	plot_names.push_back("temperature");
	for (int i=0; i<NumSpec; i++)
	{
	    plot_names.push_back("Y(" + spec_names[i] + ")");
	}
    }

    if (plot_primplus)
    {
	plot_primplus = 0;

	if (plot_magvel)
	{
	    plot_primplus = 1;
	    icomp_magvel = inext;
	    inext++;
	    plot_names.push_back("magvel");
	}
	
	if (plot_Mach)
	{
	    plot_primplus = 1;
	    icomp_Mach = inext;
	    inext++;
	    plot_names.push_back("Mach");
	}
	
	if (plot_divu)
	{
	    plot_primplus = 1;
	    icomp_divu = inext;
	    inext++;
	    plot_names.push_back("divu");
	}
	
#if (BL_SPACEDIM >= 2)
	if (plot_magvort)
	{
	    plot_primplus = 1;
	    icomp_magvort = inext;
	    inext++;
	    plot_names.push_back("magvort");
	}
#endif
	
	if (plot_X)
	{
	    plot_primplus = 1;
	    icomp_X = inext;
	    inext += NumSpec;
	    for (int i=0; i<NumSpec; i++)
	    {
		plot_names.push_back("X(" + spec_names[i] + ")");
	    }	
	}
	
	if (plot_omegadot)
	{
	    plot_primplus = 1;
	    icomp_omegadot = inext;
	    inext += NumSpec;
	    for (int i=0; i<NumSpec; i++)
	    {
		plot_names.push_back("rho*omgdot(" + spec_names[i] + ")");
	    }	
	}
	
	if (plot_dYdt)
	{
	    plot_primplus = 1;
	    icomp_dYdt = inext;
	    inext += NumSpec;
	    for (int i=0; i<NumSpec; i++)
	    {
		plot_names.push_back("Ydot(" + spec_names[i] + ")");
	    }		
	}
	
	if (plot_heatRelease)
	{
	    plot_primplus = 1;
	    icomp_heatRelease = inext;
	    inext++;
	    plot_names.push_back("HeatRelease");	
	}
	
	if (plot_fuelConsumption)
	{
	    if (fuelID >= 0)
	    {
		plot_primplus = 1;
		icomp_fuelConsumption = inext;
		inext++;
		plot_names.push_back("FuelConsumptionRate");
	    }
	    else
	    {
//	    BoxLib::Warning("plot_fuelConsumption is true, but fuelName is not set correctly.");
		plot_fuelConsumption = 0;
	    }
	}
    }
}

void
RNS::writePlotFile (const std::string& dir,
                       ostream&       os,
                       VisMF::How     how)
{
    int i, n;
    Real cur_time = state[State_Type].curTime();
    int n_plot_vars = plot_names.size();

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_plot_vars == 0)
            BoxLib::Error("Must specify at least one valid data item to plot");

        os << n_plot_vars << '\n';

	//
	// Names of variables 
	//
	for (std::vector<string>::iterator it = plot_names.begin();
	     it != plot_names.end(); ++it)
	{
	    os << *it << '\n';
        }

        os << BL_SPACEDIM << '\n';
        os << parent->cumTime() << '\n';
        int f_lev = parent->finestLevel();
        os << f_lev << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbLo(i) << ' ';
        os << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbHi(i) << ' ';
        os << '\n';
        for (i = 0; i < f_lev; i++)
            os << parent->refRatio(i)[0] << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->Geom(i).Domain() << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->levelSteps(i) << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
        {
            for (int k = 0; k < BL_SPACEDIM; k++)
                os << parent->Geom(i).CellSize()[k] << ' ';
            os << '\n';
        }
        os << (int) Geometry::Coord() << '\n';
        os << "0\n"; // Write bndry data.


        // job_info file with details about the run
	std::ofstream jobInfoFile;
	std::string FullPathJobInfoFile = dir;
	FullPathJobInfoFile += "/job_info";
	jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);	

	std::string PrettyLine = "===============================================================================\n";
	std::string OtherLine = "--------------------------------------------------------------------------------\n";
	std::string SkipSpace = "        ";

	// job information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Job Information\n";
	jobInfoFile << PrettyLine;
	
	jobInfoFile << "job name: " << job_name << "\n\n";

	jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
	jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif

	jobInfoFile << "\n";
	jobInfoFile << "CPU time used since start of simulation (CPU-hours): " <<
	  getCPUTime()/3600.0;

	jobInfoFile << "\n\n";

        // plotfile information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Plotfile Information\n";
	jobInfoFile << PrettyLine;

	time_t now = time(0);

	// Convert now to tm struct for local timezone
	tm* localtm = localtime(&now);
	jobInfoFile   << "output data / time: " << asctime(localtm);

	char currentDir[FILENAME_MAX];
	if (getcwd(currentDir, FILENAME_MAX)) {
	  jobInfoFile << "output dir:         " << currentDir << "\n";
	}

	jobInfoFile << "\n\n";

        // build information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Build Information\n";
	jobInfoFile << PrettyLine;

	jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
	jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
	jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
	jobInfoFile << "BoxLib dir:    " << buildInfoGetBoxlibDir() << "\n";

	jobInfoFile << "\n";
	
	jobInfoFile << "COMP:  " << buildInfoGetComp() << "\n";
	jobInfoFile << "FCOMP: " << buildInfoGetFcomp() << "\n";

	jobInfoFile << "\n";

	jobInfoFile << "Chemistry: " << buildInfoGetAux(1) << "\n";

	jobInfoFile << "\n";

	const char* githash1 = buildInfoGetGitHash(1);
	const char* githash2 = buildInfoGetGitHash(2);
	const char* githash3 = buildInfoGetGitHash(3);
	if (strlen(githash1) > 0) {
	  jobInfoFile << "Combustion git hash: " << githash1 << "\n";
	}
	if (strlen(githash2) > 0) {
	  jobInfoFile << "BoxLib     git hash: " << githash2 << "\n";
	}
	if (strlen(githash3) > 0) {	
	  jobInfoFile << "SDCLib     git hash: " << githash3 << "\n";
	}

	jobInfoFile << "\n\n";

	// grid information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Grid Information\n";
	jobInfoFile << PrettyLine;

	for (i = 0; i <= f_lev; i++)
	  {
	    jobInfoFile << "level: " << i << "\n";
	    jobInfoFile << "  number of boxes = " << parent->numGrids(i) << "\n";
	    jobInfoFile << "  maximum zones   = ";
	    for (n = 0; n < BL_SPACEDIM; n++)
	      {
		jobInfoFile << parent->Geom(i).Domain().length(n) << " ";
		//jobInfoFile << parent->Geom(i).ProbHi(n) << " ";
	      }
	    jobInfoFile << "\n\n";
	  }

	jobInfoFile << "\n";

	// runtime parameters
	jobInfoFile << PrettyLine;
	jobInfoFile << " Inputs File Parameters\n";
	jobInfoFile << PrettyLine;
	
	ParmParse::dumpTable(jobInfoFile, true);

	jobInfoFile.close();
    }

    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";
    char buf[64];
    sprintf(buf, "Level_%d", level);
    std::string Level = buf;
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if (!FullPath.empty() && FullPath[FullPath.size()-1] != '/')
        FullPath += '/';
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
            BoxLib::CreateDirectoryFailed(FullPath);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (i = 0; i < grids.size(); ++i)
        {
            RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
            for (n = 0; n < BL_SPACEDIM; n++)
                os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
	std::string PathNameInHeader = Level;
	PathNameInHeader += BaseName;
	os << PathNameInHeader << '\n';
    }

    int ngrow = 0;
    MultiFab plotMF(grids,n_plot_vars,ngrow);

    ngrow = (plot_divu || plot_magvort) ? 2: 1;
    for (FillPatchIterator fpi(*this, plotMF, ngrow, cur_time, State_Type, 0, NUM_STATE); 
	 fpi.isValid(); ++fpi) 
    {
	int i = fpi.index();

	if (plot_cons)
	{
	    plotMF[fpi].copy(fpi(), 0, icomp_cons, NUM_STATE);
	}

	if (plot_prim || plot_primplus)
	{
	    const Box& bx  = grids[i];
	    const Real* dx = geom.CellSize();
	    const Box& bxp = BoxLib::grow(bx,ngrow-1);

	    FArrayBox prim(bxp, NUM_STATE);
	
	    BL_FORT_PROC_CALL(RNS_CTOPRIM,rns_ctoprim)
		(bxp.loVect(), bxp.hiVect(),
		 BL_TO_FORTRAN(fpi()),
		 BL_TO_FORTRAN(prim));

	    if (plot_prim)
	    {
		plotMF[fpi].copy(prim, 0, icomp_prim, NUM_STATE);
	    }

	    if (plot_primplus)
	    {
		BL_FORT_PROC_CALL(RNS_MAKEPLOTVAR,rns_makeplotvar)
		    (bx.loVect(), bx.hiVect(), dx,
		     BL_TO_FORTRAN(prim),
		     BL_TO_FORTRAN(plotMF[fpi]),
		     n_plot_vars, icomp_magvel, icomp_Mach, icomp_divu, icomp_magvort, 
		     icomp_X, icomp_omegadot, icomp_dYdt, icomp_heatRelease, 
		     icomp_fuelConsumption, fuelID);
	    }
	}
    }

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(plotMF,TheFullPath,how,true);
}

