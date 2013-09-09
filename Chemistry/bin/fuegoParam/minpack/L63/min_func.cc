#include "min_func.hh"

// Tell the compiler about the existence of the required LAPACK functions
extern "C" {
  //Chris's template
  //void dgesv_(int *n,int *nrhs,double *a,int *lda,int *ipivot,double *b,int *ldb,int *info);
 
  // eigenvalues
  void dsyev_(char *JOBZ,char *UPLO,int *N, double *A,int *LDA, double *W, double *WORK, int *LWORK, int *INFO); 
}

int cdsyev(char jobz, char uplo, int n, double *a, int lda, double *w)
{
  int info, lwork;
  double *work;
  double work_query;

  // Ask how much space we need.
  lwork = -1; // query
  dsyev_(&jobz, &uplo, &n, a, &lda, w, &work_query, &lwork, &info);

  lwork = work_query;
  work = new double[lwork];
  dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
  delete work;

  return info;
}



min_func::min_func(double dt_,double T_,int gap_,double xo_,double yo_,double zo_,double sD_,double so_,double ro_,double bo_,double ss_,double sr_,double sb_,double st_,double rt_, double bt_)
  : dt(dt_),T(T_),gap(gap_),xo(xo_),yo(yo_),zo(zo_),sD(sD_), so(so_), ro(ro_),
	bo(bo_), ss(ss_), sr(sr_), sb(sb_), st(st_), rt(rt_), bt(bt_),
	Steps(T/dt), nObs(Steps/gap+1), xData(new double[3*nObs]), yData(xData+nObs), zData(yData+nObs) {
}

min_func::~min_func() {
	delete [] xData;
}

double min_func::drand(){
  return (rand()+1.0)/(RAND_MAX+1.0);
}

double min_func::randn(){
  double pi;
  pi =  3.14159265358979323846;
  return sqrt(-2*log(drand())) * cos(2*pi*drand());
}

void min_func::info() {
  puts("Configuration");
  printf("Time step %g\n",dt);
  printf("Integration time %g\n",T);	
  printf("Gap %i\n",gap);
  printf("Initial conditions\n %g\n %g\n %g\n",xo,yo,zo);
  printf("Standard deviation of obs %g\n",sD);
  printf("Prior means\n %g\n %g\n %g\n",so,ro,bo);
  printf("Prior standard deviation\n %g\n %g\n %g\n",ss,sr,sb);
  printf("True parameter values\n %g\n %g\n %g\n",st,rt,bt);
}

void min_func::LorenzTrajectory(double *x,double s,double b,double r){
  int J,K;
  double *k1 = new double [3];
  double *xk1 = new double [3];
  double *k2 = new double [3];
  double *xk2 = new double [3];
  double *k3 = new double [3];
  double *xk3 = new double [3];
  double *k4 = new double [3];
  double *xn = new double [3];  

  double *y=x+Steps, *z=y+Steps;
  xn[0]=xo;xn[1]=yo;xn[2]=zo;
  for(J=0;J<Steps;J++){
    // move                                                                                                                                                           
    ff(k1,xn,s,b,r);
    for(K=0;K<3;K++) xk1[K]=xn[K]+dt*0.5*k1[K];
    ff(k2,xk1,s,b,r);
    for(K=0;K<3;K++) xk2[K]=xn[K]+dt*0.5*k2[K];
    ff(k3,xk2,s,b,r);
    for(K=0;K<3;K++) xk2[K]=xn[K]+dt*k3[K];
    ff(k4,xk3,s,b,r);
    
    // save                                                                                       
    x[J] = xn[0];
    y[J] = xn[1];
    z[J] = xn[2];
    
    // update
    for(K=0;K<3;K++) xn[K]=xn[K]+(k1[K]+2*k2[K]+2*k3[K]+k4[K])*dt/6;
  }
}

void min_func::TrajectoryCloud(int NOS){
  int ii = 0;
  double *X = new double [3*Steps];
  char buf[256];
  for(ii=0;ii<NOS;ii++){
    LorenzTrajectory(X,so+ss*randn(),bo+sb*randn(),ro+sr*randn());
    sprintf(buf,"Sample%d.txt",ii);
    Save_all(buf,X,Steps);
  };

}

void min_func::ff(double *dx,double *x, double s,double b,double r){
  dx[0]=s*(x[1]-x[0]);
  dx[1]=x[0]*(r-x[2])-x[1];
  dx[2]=x[0]*x[1]-b*x[2];
}

void min_func::Lorenz(double *x,double s,double b,double r){
  int J,K;
  double *k1 = new double [3];
  double *xk1 = new double [3];
  double *k2 = new double [3];
  double *xk2 = new double [3];
  double *k3 = new double [3];
  double *xk3 = new double [3];
  double *k4 = new double [3];
  double *xn = new double [3];
  xn[0]=xo;xn[1]=yo;xn[2]=zo;
  double *y=x+nObs, *z=y+nObs;

  for(J=0;J<Steps;J++){
    // move                  
    ff(k1,xn,s,b,r);
    for(K=0;K<3;K++) xk1[K]=xn[K]+dt*0.5*k1[K];
    ff(k2,xk1,s,b,r);
    for(K=0;K<3;K++) xk2[K]=xn[K]+dt*0.5*k2[K];
    ff(k3,xk2,s,b,r);
    for(K=0;K<3;K++) xk2[K]=xn[K]+dt*k3[K];
    ff(k4,xk3,s,b,r);

    // save          
    if(J%gap==0){
      x[J/gap] = xn[0];
      y[J/gap] = xn[1];
      z[J/gap] = xn[2];
    }
    // update                                                                                                                       
    for(K=0;K<3;K++) xn[K]=xn[K]+(k1[K]+2*k2[K]+2*k3[K]+k4[K])*dt/6;
  }
  x[nObs-1] = xn[0];
  y[nObs-1] = xn[1];
  z[nObs-1] = xn[2];
}

void min_func::genData(double s,double b, double r){
  double *x = new double [3*nObs];
  Lorenz(x,s,b,r);
  double *y=x+nObs, *z=y+nObs;
  int ii = 0;
  for(ii=0;ii<nObs;ii++){
    x[ii]+=sD*randn();
    xData[ii]=x[ii];
    y[ii]+=sD*randn();
    yData[ii]=y[ii];
    z[ii]+=sD*randn();
    zData[ii]=z[ii];
  }
  Save_all("Data.txt",x,nObs);
}

void min_func::Save_all(const char *fn,double *x,int length) {
	char buf[256];
	sprintf(buf,"x%s",fn);
	Save(buf,x,length);
	sprintf(buf,"y%s",fn);
	Save(buf,x+length,length);
	sprintf(buf,"z%s",fn);
	Save(buf,x+2*length,length);
}

void min_func::Save(const char *Filename, double *x,int length){
 FILE *fp;
 int J;
 fp = fopen(Filename,"w");
 for(J=0;J<length;J++) fprintf(fp,"%g\n",x[J]);
 fclose(fp);
}


double min_func::funcF(double *X){
  double s = X[0];
  double b = X[1];
  double r = X[2];
  double F = 0;
  // Contribution from prior
  F += 0.5*( ((s-so)/ss)*((s-so)/ss) + ((b-bo)/sb)*((b-bo)/sb) + ((r-ro)/sr)*((r-ro)/sr) );
  // Contribution from data
  double *x = new double [3*nObs];
  Lorenz(x,s,b,r);
  double *y=x+nObs, *z=y+nObs;
  int ii=0;
  for(ii=0;ii<nObs;ii++){
    F +=0.5*((xData[ii]-x[ii])/sD)*((xData[ii]-x[ii])/sD)+0.5*((yData[ii]-y[ii])/sD)*((yData[ii]-y[ii])/sD)+0.5*((zData[ii]-z[ii])/sD)*((zData[ii]-z[ii])/sD); 
  } 
  return F;
}

double min_func::funcFo(double *pSample,double phi,double *Hess){
  double Fo = phi+vAv(pSample,Hess,pSample);
  return Fo;
}

void min_func::FCN(const int *NP,const double *X,double *FVEC,int *IFLAGP){
 int ii;	
 double *gradF = new double[3];	
 grad(X,gradF);	
 for(ii=0;ii<3;ii++) FVEC[ii]=gradF[ii];
}

void min_func::grad(const double *X,double *gradF){
// Finite difference	
int ii,N=3;;
  for(ii=0;ii<N;ii++){  
   gradF[ii] = der_ffd(X,N,ii); 
  } 

// Centered finite differences	
/*int ii,N=3;;
  for(ii=0;ii<N;ii++){  
   gradF[ii] = der_cfd(X,N,ii); 
  }
*/  
}

// Compute the derivative of the function funcF with respect to 
// the Kth variable (centered finite differences)
double min_func::der_cfd(double *X,int N, int K){
  int ii;
  double x1,x2,dx,h=1.5e-8;
  double fp;
  double xdX[3],xdX1[3],xdX2[3];
  for(ii=0;ii<N;ii++){
     xdX[ii]  = X[ii];
     xdX1[ii] = X[ii];
     xdX2[ii] = X[ii];
  }
    xdX1[K] = xdX[K]+h*xdX[K];
    x1 = funcF(xdX1);
    xdX2[K] = xdX[K]-h*xdX[K];
    x2 = funcF(xdX2);
    dx = xdX1[K]-xdX2[K]; 
    fp = (x1-x2)/dx;
    return fp;
}	

// Compute the derivative of the function funcF with respect to 
// the Kth variable (forward finite differences)
double min_func::der_ffd(const double *X,int N, int K){
  int ii;
  double x1,x2,dx,h=1.5e-8;
  double xdX[3],xdX1[3];
  double fp;
  for(ii=0;ii<N;ii++){
     xdX[ii] = X[ii];
     xdX1[ii]= X[ii];
  }
  xdX1[K] = xdX[K]+h*xdX[K];
  x1 = funcF(xdX1);
  x2 = funcF(xdX);
  dx = xdX1[K]-xdX[K]; 
  fp = (x1-x2)/dx;
  return fp;
}	


void min_func::hess(double *H,double s,double b, double r){
  int ii;
  double X[3]={s,b,r};
  double f1,f2,dx,h=1.5e-8;
  double X2[3];
  
  for(ii=0;ii<3;ii++) X2[ii]=X[ii];

  f2    = der_cfd(X,3,0);  
  X2[0] = X[0]+h*X[0];
  f1    = der_cfd(X2,3,0);
  dx    = X2[0]-X[0];
  X2[0] = X[0]; 
  H[0]  = (f1-f2)/dx;

  f2    = der_cfd(X,3,0);  
  X2[1] = X[1]+h*X[1];
  f1    = der_cfd(X2,3,0);
  dx    = X2[1]-X[1];
  X2[1] = X[1]; 
  H[1]  = (f1-f2)/dx;

  f2    = der_cfd(X,3,0);  
  X2[2] = X[2]+h*X[2];
  f1    = der_cfd(X2,3,0);
  dx    = X2[2]-X[2];
  X2[2] = X[2]; 
  H[2]  = (f1-f2)/dx;





  f2    = der_cfd(X,3,1);  
  X2[0] = X[0]+h*X[0];
  f1    = der_cfd(X2,3,1);
  dx    = X2[0]-X[0];
  X2[0] = X[1]; 
  H[3]  = (f1-f2)/dx;

  f2    = der_cfd(X,3,1);  
  X2[1] = X[1]+h*X[1];
  f1    = der_cfd(X2,3,1);
  dx    = X2[1]-X[1];
  X2[1] = X[1]; 
  H[4]  = (f1-f2)/dx;

  f2    = der_cfd(X,3,1);  
  X2[2] = X[2]+h*X[2];
  f1    = der_cfd(X2,3,1);
  dx    = X2[2]-X[2];
  X2[2] = X[2]; 
  H[5]  = (f1-f2)/dx;




  f2    = der_cfd(X,3,2);  
  X2[0] = X[0]+h*X[0];
  f1    = der_cfd(X2,3,2);
  dx    = X2[0]-X[0];
  X2[0] = X[0];
  H[6]  = (f1-f2)/dx;

  f2    = der_cfd(X,3,2);  
  X2[1] = X[1]+h*X[1];
  f1    = der_cfd(X2,3,2);
  dx    = X2[1]-X[1];
  X2[1] = X[1];
  H[7]  = (f1-f2)/dx;

  f2    = der_cfd(X,3,2);  
  X2[2] = X[2]+h*X[2];
  f1    = der_cfd(X2,3,2);
  dx    = X2[2]-X[2];
  X2[2] = X[2];
  H[8]  = (f1-f2)/dx;

  printf("H()= %g %g %g\n",H[0],H[1],H[2]);
  printf("     %g %g %g\n",H[3],H[4],H[5]);
  printf("     %g %g %g\n",H[6],H[7],H[8]);
}

void min_func::sqrtH(double *L,double *H){
  int ii;
  double W[3];

  int INFO = cdsyev('V', 'U', 3, H, 3, W);
  if(INFO==0){
    printf("Eigenvalue calculation ok\n");
  }
  else{
    printf("Eigenvalue calculation not ok\n");     
  }

 
  printf("%i\n",INFO);
  puts("\n Eigenvectors\n");
  printf("     %g %g %g\n",H[0],H[1],H[2]);
  printf("     %g %g %g\n",H[3],H[4],H[5]);
  printf("     %g %g %g\n",H[6],H[7],H[8]);
  puts("\n Eigenvalues\n");
  printf(" %g %g %g\n",W[0],W[1],W[2]);
 

  // (pseudo) inverting
  for(ii=0;ii<3;ii++){
    if(abs(W[ii])>1e5){
	W[ii]=0;
      }
    else{
      W[ii]=1/sqrt(W[ii]); 
      }
  }

  // matrix square root of inverse Hessian
  L[0]=H[0]*W[0];
  L[1]=H[3]*W[1];
  L[2]=H[6]*W[2];
  L[3]=H[1]*W[0];
  L[4]=H[4]*W[1];
  L[5]=H[7]*W[2];
  L[6]=H[2]*W[0];
  L[7]=H[5]*W[1];
  L[8]=H[8]*W[2];
}

double min_func::abs(double x){
  if(x<0) x=-x;
  return x;
}

void min_func::Av(double *X,double *A, double *v){
  X[1]  = A[0]*v[0] + A[1]*v[1] + A[2]*v[2];
  X[2]  = A[3]*v[0] + A[4]*v[1] + A[5]*v[2];
  X[3]  = A[6]*v[0] + A[7]*v[1] + A[8]*v[2];
}

double min_func::vAv(double *X,double *A, double *v){
  int out;
  Av(X,A,v);
  out = X[0]*v[0]+X[1]*v[1]+X[2]*v[2];
  return out;
}

min_func *global_mfp;

void min_func::minimize(double &s, double &b, double &r){
 
     int N=3,INFO,LWA=180;                                                                                  
     double TOL=1e-7,FNORM;                                                                                  
     double X[3],FVEC[3],WA[180];                                                                       
     X[0] = s;
     X[1] = b;
     X[2] = r;                                                                         

     global_mfp=this;
     hybrd1_(global_fcn,&N,X,FVEC,&TOL,&INFO,WA,&LWA);                                                         
     FNORM=enorm_(&N,FVEC);                                                                             

     if(INFO==1){
       printf("Optimization successful, res = %g\n",FNORM); 
     }
     else{
       printf("Problems with optimization. INFO = %i\n",INFO);
     }
     s = X[0];
     b = X[1];
     r = X[2];
}

double min_func::minArray(double *a,int length){
  int ii;
  double min=a[0];
  for(ii=0;ii<length;ii++){
    if(a[ii]<min) min=a[ii];
  }
  return min;
}

double min_func::maxArray(double *a,int length){
  int ii;
  double max=a[0];
  for(ii=0;ii<length;ii++){
    if(a[ii]>max) max=a[ii];
  }
  return max;
}

void min_func::sampling(int NOS,double sMAP,double bMAP,double rMAP){
  int     ii,jj;
  double *H = new double[9];
  double *Hess = new double[9];
  double *L = new double[9];
  double *xi = new double[NOS];
  double *pp = new double[NOS];
  double *ps = new double[NOS];
  double *pb = new double[NOS];
  double *pr = new double[NOS];
  double *x  = new double[Steps]; 
  char   buf[256];
  double *w = new double[NOS];
  double *pSample = new double[3];
  double phi;
  
  // compute Hessian
  hess(H,sMAP,bMAP,rMAP);
  // save Hessian for weights                                                                                                               
  for(ii=0;ii<9;ii++) Hess[ii]=H[ii];
  // compute square root (using evals)
  sqrtH(L,H);// H are now evecs of Hessian                                                                                              
  // generate samples and save them + their trajectories
  for(ii=0;ii<NOS;ii++){
    for(jj=1;jj<NOS;jj++) xi[jj]=randn();
    Av(pp,L,xi);
    ps[ii] = sMAP+pp[0];
    pb[ii] = bMAP+pp[1];
    pr[ii] = rMAP+pp[2];
    LorenzTrajectory(x,ps[ii],pb[ii],pr[ii]);
    sprintf(buf,"pSample%d.txt",ii);
    Save_all(buf,x,Steps);
  }

  // Weighting                                                                                                                            
  pSample[0]=sMAP;
  pSample[1]=bMAP;
  pSample[2]=rMAP;
  phi = funcF(pSample); //Minimum                                                                                                  
  // compute weights
  for(ii=0;ii<NOS;ii++){
    pSample[0]=ps[ii];
    pSample[1]=pb[ii];
    pSample[2]=pr[ii];
    w[ii] = -funcF(pSample);
    w[ii]+=  funcFo(pSample,phi,Hess);
  }
  // subtract max weight                                                                                                                    
  double wMax = maxArray(w,NOS);
  for(ii=0;ii<NOS;ii++){
    w[ii]-= wMax;
    w[ii] = exp(w[ii]);
  }
  // Normalize weights                                                                                                                      
  double wSum = 0;
  for(ii=0;ii<NOS;ii++) wSum+=w[ii];
  for(ii=0;ii<NOS;ii++) w[ii]=w[ii]/wSum;
  // save weights
  Save("weights.txt",w,NOS);
}

//void global_fcn(int *NP,double *X,double *FVEC,int *IFLAGP) {
void global_fcn(const int *NP,const double *X,double *FVEC,int *IFLAGP) {
	global_mfp->FCN(NP,X,FVEC,IFLAGP);
}
