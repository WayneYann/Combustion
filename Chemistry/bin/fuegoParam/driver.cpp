#include <iostream>

#include <Observation.H>

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " pmf_file=<input fab file name> [options] \n";
  std::cerr << "\tOptions:       Patm = <pressure, in atmospheres> \n";
  std::cerr << "\t                 dt = <time interval, in seconds> \n";
  std::cerr << "\t              Tfile = <T search value, in K> \n";
  exit(1);
}

// This function is in the f90 file
extern "C" {
  void test_observation(void(*)(int, const Real*,Real*),int, const Real*, Real*);
}

int
main (int   argc,
      char* argv[])
{
  BoxLib::Initialize(argc,argv);

  if (argc<2) print_usage(argc,argv);
    
  Observation ctx;

  // Populate UserContext with parameters, initialize
  // array of Reals with default values (returned by AddParameter)
  Array<Real> pdata;
  pdata.push_back(ctx.AddParameter(8,ChemDriver::FWD_BETA));
  pdata.push_back(ctx.AddParameter(8,ChemDriver::FWD_EA));
  pdata.push_back(ctx.AddParameter(8,ChemDriver::LOW_A));
  pdata.push_back(ctx.AddParameter(8,ChemDriver::LOW_BETA));
  pdata.push_back(ctx.AddParameter(8,ChemDriver::LOW_EA));
  pdata.push_back(ctx.AddParameter(8,ChemDriver::TROE_A));
  pdata.push_back(ctx.AddParameter(8,ChemDriver::TROE_TS));
  pdata.push_back(ctx.AddParameter(8,ChemDriver::TROE_TSSS));
  int cnt = ctx.NumParameters();

  // Call observation function from C++
  Real ret;
  observation_function(cnt,pdata.dataPtr(),&ret);
  std::cout << "Observation (C++): " << ret << std::endl;

  // Call observation function via function ptr passed to f90 routine
  test_observation(&observation_function,cnt,pdata.dataPtr(),&ret);
  std::cout << "Observation (F90): " << ret << std::endl;

  BoxLib::Finalize();
}

