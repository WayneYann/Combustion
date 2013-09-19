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

    if (plot_cons)
    {
	plot_names.push_back("density");
	plot_names.push_back("xmom");
#if (BL_SPACEDIM >= 2)
	plot_names.push_back("ymom");
#endif
#if (BL_SPACEDIM == 3)
	plot_names.push_back("zmom");
#endif
	plot_names.push_back("rho_E");
	plot_names.push_back("Temp");
	for (int i=0; i<NumSpec; i++)
	{
	    plot_names.push_back("rho.Y(" + spec_names[i] + ")");
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

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (plot_names.size() == 0)
            BoxLib::Error("Must specify at least one valid data item to plot");

        os << plot_names.size() << '\n';

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
    MultiFab plotMF(grids,plot_names.size(),ngrow);

    ngrow = 1;
    for (FillPatchIterator fpi(*this, plotMF, ngrow, cur_time, State_Type, 0, NUM_STATE); 
	 fpi.isValid(); ++fpi) 
    {
	if (plot_cons)
	{
	    plotMF[fpi].copy(fpi());
	}

	
    }

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(plotMF,TheFullPath,how,true);
}

