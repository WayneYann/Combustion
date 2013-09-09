#ifndef MIN_FUNC_HH
#define MIN_FUNC_HH

#include <cstdio>
#include <cstdlib>
#include <math.h>

#include "minpack.h"

class min_func {

public:

  const double dt,T;         
  const int gap;
  const double xo,yo,zo;
  const double sD,so,ro,bo;       
  const double ss,sr,sb;
  const double st,rt,bt;
  const int Steps;		
  const int nObs;
  double *xData,*yData,*zData;
		
  
  min_func(double dt_,double T_,int gap_,double xo_,double yo_,double zo_,double sD_,double so_,double ro_,double bo_,double ss_,double sr_,double sb_,double st_, double rt_, double bt_); 
 ~min_func();
  double drand();              
  double randn();              
  void info();
  void ff(double *dx,double *x, double s,double b,double r);
  void LorenzTrajectory(double *x,double s,double b,double r);
  void Lorenz(double *x,double s,double b,double r);
  void TrajectoryCloud(int NOS);
  void genData(double s,double b,double r);
  void Save(const char *Filename,double *x,int length);
  void Save_all(const char *fn,double *x,int length);
  double funcF(double *X);
  void FCN(const int *NP,const double *X,double *FVEC,int *IFLAGP);
  void minimize(double &s, double &b, double &r);
  double minArray(double *a,int length);
  double maxArray(double *a,int length);
  void hess(double *H,double s,double b, double r);
  void grad(const double *X,double *gradF);
  void sqrtH(double *L, double *H);
  void Av(double *X,double *A, double *v);
  double abs(double x);
  double funcFo(double *pSample,double phi,double *Hess);
  double vAv(double *X,double *A, double *v);
  void sampling(int NOS,double sMAP,double bMAP,double rMAP);
  double der_cfd(double *X,int N, int i);
  double der_ffd(const double *X,int N, int K);
private:

};

void global_fcn(const int *NP,const double *X,double *FVEC,int *IFLAGP);

#endif
