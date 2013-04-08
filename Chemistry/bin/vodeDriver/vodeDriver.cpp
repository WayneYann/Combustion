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
   std::cerr << "\tOptions:       Patm = <pressure, in atmospheres>              [DEFAULT = " << Patm_DEF << "]\n";
   std::cerr << "\t                 dt = <time interval, in seconds>             [DEFAULT = " << dt_DEF << "]\n";
   std::cerr << "\t        fabfile_out = <output fab file name, null->no output> [DEFAULT = \"\"""]\n";
   std::cerr << "\t            verbose = <0,1>                                   [DEFAULT = " << verbose_DEF << "]\n";
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

    Real Patm = Patm_DEF; pp.query("Patm",Patm);
    Real dt = dt_DEF; pp.query("dt",dt);

    funcCnt.setVal(0);
    nstate.copy(ostate,box,sCompT,box,sCompT,1);

#ifdef LMC_SDC
    const int sCompR = 3;
    ostate.mult(1.e3,sCompR,1);
    const int sCompRH = 2; // Replace cgs velocity with MKS rho.Hmix
    cd.getHmixGivenTY(ostate,ostate,ostate,box,sCompT,sCompY,sCompRH);
    ostate.mult(ostate,sCompR,sCompRH,1);
    for (int i=0; i<nSpec; ++i) {
      ostate.mult(ostate,sCompR,sCompY+i,1);
    }
    FArrayBox c_0(box,nSpec+1); c_0.setVal(0);
    FArrayBox I_R(box,nSpec+1);
    cd.solveTransient_sdc(nstate,nstate,nstate,ostate,ostate,ostate,
                          c_0,I_R,funcCnt,box,sCompY,sCompRH,sCompT,dt,Patm,0,true);
#else
    cd.solveTransient(nstate,nstate,ostate,ostate,funcCnt,
                      box,sCompY,sCompT,dt,Patm);
#endif

    std::cout << " ... total function evals: " << funcCnt.norm(1) << '\n';
    std::cout << " ... max evals at a point: " << funcCnt.norm(0) << '\n';

    std::string fabfile_out=""; pp.query("fabfile_out",fabfile_out);
    if (fabfile_out != "") {
      std::ofstream os;
      os.open(fabfile_out.c_str());
      nstate.writeOn(os);
      os.close();
    }
    
    BoxLib::Finalize();
}

