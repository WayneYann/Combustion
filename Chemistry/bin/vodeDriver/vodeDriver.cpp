#include <iostream>
#include <fstream>
using std::cout;
using std::endl;
using std::string;

#include "Utility.H"
#include "ParallelDescriptor.H"
#include "ChemDriver.H"
#include "ParmParse.H"

bool verbose_DEF=true;

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);
    
    // Parse command line
    ParmParse pp;
    bool verbose=verbose_DEF; pp.query("verbose",verbose);
    
    ChemDriver cd;
    if (verbose)
        cd.set_verbose_vode();

    const Array<std::string>& names = cd.speciesNames();
    Array<Real> mwt = cd.speciesMolecWt();
    const int nSpec = cd.numSpecies();

    std::string pmf_file=""; pp.get("pmf_file",pmf_file);
    std::ifstream is;
    is.open(pmf_file.c_str());
    FArrayBox ostate;
    ostate.readFrom(is);
    is.close();

    const Box& box = ostate.box();
    const int nComp = nSpec + 4;
    if (nComp != ostate.nComp()) {
      cout << "pmf file is not compatible with the mechanism compiled into this code" << endl;
      BoxLib::Abort();
    }
    FArrayBox nstate(box,nComp);

    const int sCompY = 4;
    const int sCompT = 1;
    FArrayBox funcCnt(box,1);

    Real Patm = 1.0; pp.query("Patm",Patm);
    Real dt = 1.0e-5; pp.query("dt",dt);

    funcCnt.setVal(0);
    cd.solveTransient(nstate,nstate,ostate,ostate,funcCnt,
                      box,sCompY,sCompT,dt,Patm);

    cout << " ... total function evals: " << funcCnt.norm(1) << endl;
    cout << " ... max evals at a point: " << funcCnt.norm(0) << endl;

    std::string write=""; pp.query("write",write);
    if (write != "") {
      std::ofstream os;
      os.open(write.c_str());
      nstate.writeOn(os);
      os.close();
    }
    
    BoxLib::Finalize();
}

    
