#include <cstdio>

#include "min_func.hh"

int main() {
  const double dt = 0.005; // Time step
  const double T  = 4; //Simulation time
  //const int    gap = 100; // # of steps between obs
  const int    gap = 10; // # of steps between obs
  const double xo = 3.6314, yo = 6.6136,zo = 10.6044; // initial conditions 
  const double sD = 0.1; // standard deviation of observation noise
  const double so = 9.5, ro = 27, bo = 2.5; // means of prior
  const double ss = 0.4*so, sr = 0.4*ro, sb = 0.4*bo; // standard deviations  
  const double st = 10, rt = 28, bt = 8./3.; // true parameter values

  min_func mf(dt,T,gap,xo,yo,zo,sD,so,ro,bo,ss,sr,sb,st,rt,bt);
  mf.info();

  // Generate true state trajectory
  double *x = new double [3*mf.Steps];   
  mf.LorenzTrajectory(x,st,bt,rt);
  mf.Save_all("True.txt",x,mf.Steps);
  
  // Generate the data
  mf.genData(st,bt,rt);  
  
  // Mean of prior
  mf.LorenzTrajectory(x,so,bo,ro);
  mf.Save_all("Prior.txt",x,mf.Steps);

  // Generate a cloud of trajectories from prior
  mf.TrajectoryCloud(100);
 
  // Optimization 
  double sMAP=so,rMAP=ro,bMAP=bo;
  mf.minimize(sMAP,bMAP,rMAP);
  mf.LorenzTrajectory(x,sMAP,bMAP,rMAP);    
  mf.Save_all("MAP.txt",x,mf.Steps); 
  // Output result of minimization
  printf("Maximum likelihood parameters \n");
  printf("%g \n",sMAP);
  printf("%g \n",rMAP);
  printf("%g \n",bMAP);

  // Sampling
  mf.sampling(1000,sMAP,bMAP,rMAP);
}

