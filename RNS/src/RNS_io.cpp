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
    AmrLevel::setPlotVariables();
    
    const Array<std::string>& names = chemSolve->speciesNames();
    
    ParmParse pp("rns");
    
    bool plot_rhoY,plot_massFrac,plot_moleFrac,plot_conc;
    plot_rhoY=plot_massFrac=plot_moleFrac=plot_conc = false;
    
    if (pp.query("plot_massfrac",plot_massFrac))
    {
        if (plot_massFrac)
        {
            for (int i = 0; i < names.size(); i++)
            {
                const std::string name = "Y("+names[i]+")";
                parent->addDerivePlotVar(name);
            }
        }
        else
        {
            for (int i = 0; i < names.size(); i++)
            {
                const std::string name = "Y("+names[i]+")";
                parent->deleteDerivePlotVar(name);
            }
        }
    }

    if (pp.query("plot_molefrac",plot_moleFrac))
    {
        if (plot_moleFrac)
	{
            parent->addDerivePlotVar("molefrac");
	}
        else
	{
            parent->deleteDerivePlotVar("molefrac");
	}
    }
    
    if (pp.query("plot_concentration",plot_conc))
    {
        if (plot_conc)
	{
            parent->addDerivePlotVar("concentration");
	}
        else
	{
            parent->deleteDerivePlotVar("concentration");
	}
    }
    
    if (pp.query("plot_rhoY",plot_rhoY))
    {
        if (plot_rhoY)
        {
            for (int i = 0; i < names.size(); i++)
            {
                const std::string name = "rho.Y("+names[i]+")";
                parent->addStatePlotVar(name);
            }
        }
        else
        {
            for (int i = 0; i < names.size(); i++)
            {
                const std::string name = "rho.Y("+names[i]+")";
                parent->deleteStatePlotVar(name);
            }
        }
    }

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "\nState Plot Vars: ";

        std::list<std::string>::const_iterator li = parent->statePlotVars().begin(), end = parent->statePlotVars().end();

        for ( ; li != end; ++li)
            std::cout << *li << ' ';
        std::cout << '\n';

        std::cout << "\nDerive Plot Vars: ";

        li  = parent->derivePlotVars().begin();
        end = parent->derivePlotVars().end();

        for ( ; li != end; ++li)
            std::cout << *li << ' ';
        std::cout << '\n';
    }

}

void
RNS::writePlotFile (const std::string& dir,
                       ostream&       os,
                       VisMF::How     how)
{
    int i, n, ncomp;
    //
    // The list of indices of State to write to plotfile.
    // first component of pair is state_type,
    // second component of pair is component # within the state_type
    //
    std::vector<std::pair<int,int> > plot_var_map;
    for (int typ = 0; typ < desc_lst.size(); typ++)
    {
        for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
	{
            if (parent->isStatePlotVar(desc_lst[typ].name(comp)) &&
                desc_lst[typ].getType() == IndexType::TheCellType())
	    {
                plot_var_map.push_back(std::pair<int,int>(typ,comp));
	    }
	}
    }

    int num_derive = 0;
    std::list<std::string> derive_names;
    const std::list<DeriveRec>& dlist = derive_lst.dlist();

    for (std::list<DeriveRec>::const_iterator it = dlist.begin();
	 it != dlist.end();
	 ++it)
    {
        if (parent->isDerivePlotVar(it->name()))
	{
            derive_names.push_back(it->name());
            num_derive += it->numDerive();
	}
    }

    int n_data_items = plot_var_map.size() + num_derive;

    Real cur_time = state[State_Type].curTime();

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_data_items == 0)
            BoxLib::Error("Must specify at least one valid data item to plot");

        os << n_data_items << '\n';

	//
	// Names of variables -- first state, then derived
	//
	for (i =0; i < plot_var_map.size(); i++)
        {
	    int typ = plot_var_map[i].first;
	    int comp = plot_var_map[i].second;
	    os << desc_lst[typ].name(comp) << '\n';
        }

	for ( std::list<std::string>::iterator it = derive_names.begin();
	      it != derive_names.end(); ++it)
        {
	    const DeriveRec* rec = derive_lst.get(*it);
	    for (i = 0; i < rec->numDerive(); i++)
                os << rec->variableName(i) << '\n';
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
        if (n_data_items > 0)
        {
            std::string PathNameInHeader = Level;
            PathNameInHeader += BaseName;
            os << PathNameInHeader << '\n';
        }
    }
    //
    // We combine all of the multifabs -- state, derived, etc -- into one
    // multifab -- plotMF.
    // NOTE: we are assuming that each state variable has one component,
    // but a derived variable is allowed to have multiple components.
    int       cnt   = 0;
    const int nGrow = 0;
    MultiFab  plotMF(grids,n_data_items,nGrow);
    MultiFab* this_dat = 0;
    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (i = 0; i < plot_var_map.size(); i++)
    {
	int typ  = plot_var_map[i].first;
	int comp = plot_var_map[i].second;
	this_dat = &state[typ].newData();
	MultiFab::Copy(plotMF,*this_dat,comp,cnt,1,nGrow);
	cnt++;
    }
    //
    // Cull data from derived variables.
    // 
    if (derive_names.size() > 0)
    {
	for (std::list<std::string>::iterator it = derive_names.begin();
	     it != derive_names.end(); ++it) 
	{
	    const DeriveRec* rec = derive_lst.get(*it);
	    ncomp = rec->numDerive();
	    MultiFab* derive_dat = derive(*it,cur_time,nGrow);
	    MultiFab::Copy(plotMF,*derive_dat,0,cnt,ncomp,nGrow);
	    delete derive_dat;
	    cnt += ncomp;
	}
    }

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(plotMF,TheFullPath,how,true);
}

