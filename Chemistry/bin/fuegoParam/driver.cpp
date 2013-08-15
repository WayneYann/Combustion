#include <iostream>
#include <fstream>

#include "Utility.H"
#include "ParallelDescriptor.H"
#include "ChemDriver.H"
#include "ParmParse.H"

static Real Patm_DEF    = 1;
static Real dt_DEF      = 1.e-8;

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " pmf_file=<input fab file name> [options] \n";
  std::cerr << "\tOptions:       Patm = <pressure, in atmospheres>              [DEFAULT = " << Patm_DEF << "]\n";
  std::cerr << "\t                 dt = <time interval, in seconds>             [DEFAULT = " << dt_DEF << "]\n";
  exit(1);
}

typedef ChemDriver::Parameter Param;
typedef PArray<Param> PPArray;


struct UserContext
{
  PPArray active_params; // The set of active parameters
  FArrayBox s_init, s_final, C_0, I_R, funcCnt; 
  int sCompY, sCompT, sCompR, sCompRH;
  Real Patm, dt;
  ChemDriver* cd;
};

// The observation function
Real ObservationFunction(int         n,
                         const Real* pvals,
                         void*       ctx);

// Helper Functions
void Evolve(void* ptr);
void InitializeAuxilaryData(ParmParse& pp,
                            void*      ptr);

int
main (int   argc,
      char* argv[])
{
  BoxLib::Initialize(argc,argv);

  if (argc<2) print_usage(argc,argv);
    
  ChemDriver cd;

  // Build user context with active parameters and auxiliary data
  UserContext ctx;
  ctx.cd = &cd;
  ctx.sCompT  = 1;
  ctx.sCompRH = 2;
  ctx.sCompR  = 3;
  ctx.sCompY  = 4;

  ctx.active_params.resize(9,PArrayManage);
  int cnt = 0;
  ctx.active_params.set(cnt++,new Param(8,ChemDriver::FWD_A));
  ctx.active_params.set(cnt++,new Param(8,ChemDriver::FWD_BETA));
  ctx.active_params.set(cnt++,new Param(8,ChemDriver::FWD_EA));
  ctx.active_params.set(cnt++,new Param(8,ChemDriver::LOW_A));
  ctx.active_params.set(cnt++,new Param(8,ChemDriver::LOW_BETA));
  ctx.active_params.set(cnt++,new Param(8,ChemDriver::LOW_EA));
  ctx.active_params.set(cnt++,new Param(8,ChemDriver::TROE_A));
  ctx.active_params.set(cnt++,new Param(8,ChemDriver::TROE_TS));
  ctx.active_params.set(cnt++,new Param(8,ChemDriver::TROE_TSSS));

  // Initialize active parameters to default values
  Array<Real> pdata(cnt);
  for (int i=0; i<cnt; ++i) {
    pdata[i] = ctx.active_params[i].DefaultValue();
  }

  ParmParse pp;
  void* user = (void*)(&ctx);
  InitializeAuxilaryData(pp,user);

  // Call observation function
  const Real* pvals = pdata.dataPtr();
  Real ret = ObservationFunction(cnt,pvals,user);
  std::cout << "Total time for observation: " << ret << std::endl;

  BoxLib::Finalize();
}

Real ObservationFunction(int         n,
                         const Real* pvals,
                         void*       ptr)
{
  // Returns CPU time to evolve...useless, but a placeholder
  UserContext* ctx = (UserContext*)(ptr);

  // Load active parameters
  for (int i=0; i<n; ++i) {
    ctx->active_params[i] = pvals[i];
  }

  double strt_time = ParallelDescriptor::second();
  ctx->funcCnt.setVal(0);
  Evolve(ptr);
  double evolve_time = ParallelDescriptor::second() - strt_time;
  
  std::cout << " ... total function evals: " << ctx->funcCnt.norm(1) << '\n';
  std::cout << " ... max evals at a point: " << ctx->funcCnt.norm(0) << '\n';

  return evolve_time;
}

void Evolve(void* ptr)
{
  UserContext* ctx = (UserContext*)(ptr);
  FArrayBox& Fcnt = ctx->funcCnt;
  const Box& box = Fcnt.box();

#ifdef LMC_SDC
  FArrayBox& rYold = ctx->s_init;
  FArrayBox& rYnew = ctx->s_final;
  FArrayBox& rHold = ctx->s_init;
  FArrayBox& rHnew = ctx->s_final;
  FArrayBox& Told  = ctx->s_init;
  FArrayBox& Tnew  = ctx->s_final;
  FArrayBox& C_0   = ctx->C_0;
  FArrayBox& I_R   = ctx->I_R;
  FArrayBox* diag = 0;
  ctx->cd->solveTransient_sdc(rYnew,rHnew,Tnew,rYold,rHold,Told,C_0,I_R,
                              Fcnt,box,ctx->sCompY,ctx->sCompRH,ctx->sCompT,
                              ctx->dt,ctx->Patm,0,true);
#else
  FArrayBox& Yold = ctx->s_init;
  FArrayBox& Ynew = ctx->s_final;
  FArrayBox& Told = ctx->s_init;
  FArrayBox& Tnew = ctx->s_final;
  cd->solveTransient(Ynew,Tnew,Yold,Told,Fcnt,box,
                     ctx->sCompY,ctx->sCompT,ctx->dt,ctx->Patm);
#endif
}

void InitializeAuxilaryData(ParmParse& pp,
                            void*      ptr)
{
  UserContext* ctx = (UserContext*)(ptr);

  std::string pmf_file=""; pp.get("pmf_file",pmf_file);
  std::ifstream is;
  is.open(pmf_file.c_str());
  ctx->s_init.readFrom(is);
  is.close();

  // Simple check to see if number of species is same between compiled mech and fab file
  const Box& box = ctx->s_init.box();
  const int nSpec = ctx->cd->numSpecies();
  const int nComp = nSpec + 4;
  if (nComp != ctx->s_init.nComp()) {
    std::cout << "pmf file is not compatible with the mechanism compiled into this code" << '\n';
    BoxLib::Abort();
  }

  ctx->Patm = Patm_DEF; pp.query("Patm",ctx->Patm);
  ctx->dt = dt_DEF; pp.query("dt",ctx->dt);
  ctx->funcCnt.resize(box,1);
  ctx->s_final.resize(box,ctx->s_init.nComp());

#ifdef LMC_SDC
  ctx->s_init.mult(1.e3,ctx->sCompR,1);
  ctx->cd->getHmixGivenTY(ctx->s_init,ctx->s_init,ctx->s_init,
                         box,ctx->sCompT,ctx->sCompY,ctx->sCompRH);
  ctx->s_init.mult(ctx->s_init,ctx->sCompR,ctx->sCompRH,1);
  for (int i=0; i<nSpec; ++i) {
    ctx->s_init.mult(ctx->s_init,ctx->sCompR,ctx->sCompY+i,1);
  }
  ctx->C_0.resize(box,nSpec+1); ctx->C_0.setVal(0);
  ctx->I_R.resize(box,nSpec+1);
#endif
}
