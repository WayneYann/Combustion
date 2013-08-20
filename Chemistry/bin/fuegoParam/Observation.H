#ifndef Observation_H_
#define Observation_H_

#include <ChemDriver.H>

struct Observation
{
public:
  Observation(); // Initialize simulated observation
  ~Observation();

  // Add parameter to active set, return default value
  Real AddParameter(int reaction, const ChemDriver::REACTION_PARAMETER& rp);

  // Retrieve Parameter datastructure
  ChemDriver::Parameter& Parameter(int i);
  const ChemDriver::Parameter& Parameter(int i) const;

  // Return length of current active set
  int NumParameters() const;

  // Perform simulation to create "observation" and return value,
  // called by observation_function after pushing new parameter values
  Real Evolve();

private:
  // Reset internal data back to state of initialization, in preparation for
  // next call to observation_function
  void Reset();

  // Compute observation from final data after evolve
  Real FinalValue() const;

  PArray<ChemDriver::Parameter> active_params; // The set of active parameters
  FArrayBox s_init, s_final, C_0, I_R, funcCnt; 
  int sCompY, sCompT, sCompR, sCompRH;
  Real Patm, dt;
  ChemDriver* cd;
};


// The observation function
extern "C" {
  void observation_function(int num_vals, const Real* pvals, Real* y);
}
#endif
