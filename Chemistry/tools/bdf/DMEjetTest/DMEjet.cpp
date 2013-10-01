
#include <winstd.H>

#include <new>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#ifndef WIN32
#include <unistd.h>
#endif

#include <Utility.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <FArrayBox.H>

#include <ChemDriver.H>

#include "DMEjet_F.H"

static int         verbose = 0;
static std::string fabin("DMEjet.fab");
static std::string fabout_vode("vode.fab");
static std::string fabout_bdf("bdf.fab");
static Real stop_time = 4.e-8;
static Real dt        = 4.e-10;

int main (int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    if (argc==1)
    {
	std::cout << "Since no inputs file is provided, I am going to use the default parameters specified in DMEjet.cpp\n" << std::endl;
    }

    {
	ParmParse pp;
	pp.query("verbose", verbose);
	pp.query("fabin", fabin);
	pp.query("fabout_vode", fabout_vode);
	pp.query("fabout_bdf", fabout_bdf);
	pp.query("stop_time", stop_time);
	pp.query("dt", dt);
    }

    FArrayBox inData;
    {
	std::ifstream is;
	is.open(fabin.c_str(), std::ios::in);
	if (is)
	{
	    inData.readFrom(is);
	    is.close();
	}
	else
	{
	    std::cerr << "Error: cannot open file " << fabin <<std::endl;
	    BoxLib::Abort();
	}
    }
    
    const Box& box  = inData.box();
    const int ncomp = inData.nComp();
    const int* lo = box.loVect();
    const int* hi = box.hiVect();

    for (int i=0; i<2; i++)
    {
	int use_vode = (i==0) ? 1 : 0;

	FArrayBox outData(box, ncomp);
	ChemDriver chemdriver(use_vode);
	if (chemdriver.numSpecies()+2 != ncomp)
	{
	    std::cerr << "Error: # of species do not match " << chemdriver.numSpecies()
		      << " " << ncomp-2 <<std::endl;
	    BoxLib::Abort();	    
	}

	if (use_vode)
	{
	    std::cout << "VODE: ..." << std::endl;
	}
	else
	{
	    std::cout << "BDF: ..." << std::endl;
	}

	double strt_time = ParallelDescriptor::second();

	BL_FORT_PROC_CALL(CHEMSOLVE, chemsolve)
	    (lo, hi, BL_TO_FORTRAN(inData), BL_TO_FORTRAN(outData), dt, verbose, use_vode);

	double stop_time = ParallelDescriptor::second();

	std::cout << " ... solve time: " << stop_time-strt_time << "\n" << std::endl;

	std::string fabout = (use_vode) ? fabout_vode : fabout_bdf;
	std::ofstream os;
	os.open(fabout.c_str(), std::ios::trunc);
	if (os)
	{
	    outData.writeOn(os);
	    os.close();
	}
	else
	{
	    std::cerr << "Error: cannot open file " << fabout <<std::endl;
	    BoxLib::Abort();
	}

	chemdriver.reset();
    }

    BoxLib::Finalize();
    return 0;
}
