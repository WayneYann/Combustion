#include <fstream>

#include <Utility.H>
#include <ParallelDescriptor.H>
#include <ParmParse.H>

#include <Observation.H>

static Real Patm_DEF = 1;
static Real dt_DEF   = 1.e-6;
static Real Tfile_DEF = 1000;

static Observation* the_obs_ptr = 0;

extern "C" {
  void observation_function(int num_vals, const Real* pvals, Real* y)
  {
    BL_ASSERT(the_obs_ptr);
    Observation& ctx = *the_obs_ptr;
    
    // Load active parameters
    for (int i=0; i<num_vals; ++i) {
      ctx.Parameter(i) = pvals[i];
    }

    *y = ctx.Evolve();
  }
}

Real
Observation::Evolve()
{
  // Evolve the initial state for time interval = dt(s), at pressure Patm (atm)
  // Return final value of evolved state
  Reset();
  const Box& box = funcCnt.box();

#ifdef LMC_SDC
  FArrayBox& rYold = s_init;
  FArrayBox& rYnew = s_final;
  FArrayBox& rHold = s_init;
  FArrayBox& rHnew = s_final;
  FArrayBox& Told  = s_init;
  FArrayBox& Tnew  = s_final;
  FArrayBox* diag = 0;
  cd->solveTransient_sdc(rYnew,rHnew,Tnew,rYold,rHold,Told,C_0,I_R,
                         funcCnt,box,sCompY,sCompRH,sCompT,
                         dt,Patm,diag,true);
#else
  FArrayBox& Yold = s_init;
  FArrayBox& Ynew = s_final;
  FArrayBox& Told = s_init;
  FArrayBox& Tnew = s_final;
  cd->solveTransient(Ynew,Tnew,Yold,Told,funcCnt,box,
                          sCompY,sCompT,dt,Patm);
#endif

  return FinalValue();
}

Real
Observation::FinalValue() const
{
  // Return the final temperature of the cell that was evolved
  return s_final(s_final.box().smallEnd(),sCompT);
}

Observation::Observation()
  : cd(0)
{
  if (the_obs_ptr != 0) {
    BoxLib::Abort("Currently, may only have one Observation object constructed at a time");
  }
  the_obs_ptr = this;
  cd = new ChemDriver();

  ParmParse pp;
  std::string pmf_file=""; pp.get("pmf_file",pmf_file);
  std::ifstream is;
  is.open(pmf_file.c_str());
  FArrayBox fileFAB;
  fileFAB.readFrom(is);
  is.close();

  // Simple check to see if number of species is same between compiled mech and fab file
  const Box& box = fileFAB.box();
  const int nSpec = cd->numSpecies();
  const int nComp = nSpec + 4;
  if (nComp != fileFAB.nComp()) {
    std::cout << "pmf file is not compatible with the mechanism compiled into this code" << '\n';
    BoxLib::Abort();
  }
  sCompT  = 1;
  sCompRH = 2;
  sCompR  = 3;
  sCompY  = 4;

  // Find location
  bool found = false;
  IntVect iv=box.smallEnd();
  Real Tfile = Tfile_DEF; pp.query("Tfile",Tfile);
  for (IntVect End=box.bigEnd(); iv<=End && !found; box.next(iv)) {
    if (fileFAB(iv,sCompT)>=Tfile) found = true;
  }

  Box bx(iv,iv);
  Patm = Patm_DEF; pp.query("Patm",Patm);
  dt = dt_DEF; pp.query("dt",dt);
  s_init.resize(bx,fileFAB.nComp()); s_init.copy(fileFAB);
  funcCnt.resize(bx,1);
  
#ifdef LMC_SDC
  s_init.mult(1.e3,sCompR,1);
  cd->getHmixGivenTY(s_init,s_init,s_init,bx,sCompT,sCompY,sCompRH);
  s_init.mult(s_init,sCompR,sCompRH,1);
  for (int i=0; i<nSpec; ++i) {
    s_init.mult(s_init,sCompR,sCompY+i,1);
  }
  C_0.resize(bx,nSpec+1); C_0.setVal(0);
  I_R.resize(bx,nSpec+1); I_R.setVal(0);
#endif

  s_final.resize(bx,s_init.nComp());
  s_final.copy(s_init);
}

Observation::~Observation()
{
  delete cd;
  active_params.clear();
  the_obs_ptr = 0;
}

Real
Observation::AddParameter(int reaction, 
                          const ChemDriver::REACTION_PARAMETER& rp)
{
  int len = active_params.size();
  active_params.resize(len+1,PArrayManage);
  active_params.set(len, new ChemDriver::Parameter(reaction,rp));
  return active_params[len].DefaultValue();
}

ChemDriver::Parameter&
Observation::Parameter(int i)
{
  return active_params[i];
}

const ChemDriver::Parameter&
Observation::Parameter(int i) const
{
  return active_params[i];
}

int
Observation::NumParameters() const
{
  return active_params.size();
}

void
Observation::Reset()
{
  funcCnt.setVal(0);
}
