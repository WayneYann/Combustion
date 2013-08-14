#include <iostream>
#include <fstream>

#include "Utility.H"
#include "ParallelDescriptor.H"
#include "ChemDriver.H"
#include "ParmParse.H"

static Real Patm_DEF    = 1;
static Real dt_DEF      = 1.e-5;
static bool verbose_DEF = false;

static
void 
print_usage (int,
             char* argv[])
{
   std::cerr << "usage:\n";
   std::cerr << argv[0] << " pmf_file=<input fab file name> [options] \n";
   exit(1);
}

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);

    if (argc<2) print_usage(argc,argv);
    
    // Parse command line
    ParmParse pp;
    bool verbose=verbose_DEF; pp.query("verbose",verbose);
    
    ChemDriver cd;
    if (verbose)
        cd.set_verbose_vode();

    // Read fab containing pmf solution
    std::string pmf_file=""; pp.get("pmf_file",pmf_file);
    std::ifstream is;
    is.open(pmf_file.c_str());
    FArrayBox ostate;
    ostate.readFrom(is);
    is.close();

    // Simple check to see if number of species is same between compiled mech and fab file
    const Box& box = ostate.box();
    const int nSpec = cd.numSpecies();
    const int nComp = nSpec + 4;
    if (nComp != ostate.nComp()) {
      std::cout << "pmf file is not compatible with the mechanism compiled into this code" << '\n';
      BoxLib::Abort();
    }
    FArrayBox nstate(box,nComp);

    const int sCompY = 4; // An assumption...
    const int sCompT = 1; // Another assumption...
    FArrayBox funcCnt(box,1);

    funcCnt.setVal(0);
    nstate.copy(ostate,box,sCompT,box,sCompT,1);

    //for (int i=0; i<cd.numReactions(); ++i) {
    for (int i=0; i<1; ++i) {
      ChemDriver::Parameter param(i,ChemDriver::FWD_A);
      std::cout << param << std::endl;
      param = 4;
      std::cout << param << std::endl;
      std::cout << ChemDriver::Parameter(i,ChemDriver::FWD_A) << std::endl;
      param.ResetToDefault();
      std::cout << ChemDriver::Parameter(i,ChemDriver::FWD_A) << std::endl;
      param = 6;
      std::cout << param << std::endl;
      cd.ResetAllParamsToDefault();
      std::cout << param << std::endl;

      std::cout << ChemDriver::Parameter(i,ChemDriver::FWD_BETA) << std::endl;
      std::cout << ChemDriver::Parameter(i,ChemDriver::FWD_EA) << std::endl;
      std::cout << ChemDriver::Parameter(i,ChemDriver::LOW_A) << std::endl;
      std::cout << ChemDriver::Parameter(i,ChemDriver::LOW_BETA) << std::endl;
      std::cout << ChemDriver::Parameter(i,ChemDriver::LOW_EA) << std::endl;
      std::cout << ChemDriver::Parameter(i,ChemDriver::REV_A) << std::endl;
      std::cout << ChemDriver::Parameter(i,ChemDriver::REV_BETA) << std::endl;
      std::cout << ChemDriver::Parameter(i,ChemDriver::REV_EA) << std::endl;
      std::cout << ChemDriver::Parameter(i,ChemDriver::TROE_A) << std::endl;
      std::cout << ChemDriver::Parameter(i,ChemDriver::TROE_TS) << std::endl;
      std::cout << ChemDriver::Parameter(i,ChemDriver::TROE_TSS) << std::endl;
      std::cout << ChemDriver::Parameter(i,ChemDriver::TROE_TSSS) << std::endl;
      std::cout << ChemDriver::Parameter(i,ChemDriver::SRI_A) << std::endl;
      std::cout << ChemDriver::Parameter(i,ChemDriver::SRI_B) << std::endl;
      std::cout << ChemDriver::Parameter(i,ChemDriver::SRI_C) << std::endl;
      std::cout << ChemDriver::Parameter(i,ChemDriver::SRI_D) << std::endl;
      std::cout << ChemDriver::Parameter(i,ChemDriver::SRI_E) << std::endl;
      std::cout << std::endl;
    }
    BoxLib::Finalize();
}

