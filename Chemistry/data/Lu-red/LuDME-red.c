
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#if defined(BL_FORT_USE_UPPERCASE)
#define FORT_CKWYP DSCKWYP
#define CKINDX CKINDX
#define CKINIT CKINIT
#define CKXNUM CKXNUM
#define CKSYME CKSYME
#define CKSYMS CKSYMS
#define CKRP CKRP
#define CKPX CKPX
#define CKPY CKPY
#define CKPC CKPC
#define CKRHOX CKRHOX
#define CKRHOY CKRHOY
#define CKRHOC CKRHOC
#define CKWT CKWT
#define CKMMWY CKMMWY
#define CKMMWX CKMMWX
#define CKMMWC CKMMWC
#define CKYTX CKYTX
#define CKYTCP CKYTCP
#define CKYTCR CKYTCR
#define CKXTY CKXTY
#define CKXTCP CKXTCP
#define CKXTCR CKXTCR
#define CKCTX CKCTX
#define CKCTY CKCTY
#define CKCPOR CKCPOR
#define CKHORT CKHORT
#define CKSOR CKSOR
#define CKCVML CKCVML
#define CKCPML CKCPML
#define CKUML CKUML
#define CKHML CKHML
#define CKGML CKGML
#define CKAML CKAML
#define CKSML CKSML
#define CKCVMS CKCVMS
#define CKCPMS CKCPMS
#define CKUMS CKUMS
#define CKHMS CKHMS
#define CKGMS CKGMS
#define CKAMS CKAMS
#define CKSMS CKSMS
#define CKCPBL CKCPBL
#define CKCPBS CKCPBS
#define CKCVBL CKCVBL
#define CKCVBS CKCVBS
#define CKHBML CKHBML
#define CKHBMS CKHBMS
#define CKUBML CKUBML
#define CKUBMS CKUBMS
#define CKSBML CKSBML
#define CKSBMS CKSBMS
#define CKGBML CKGBML
#define CKGBMS CKGBMS
#define CKABML CKABML
#define CKABMS CKABMS
#define CKWC CKWC
#define CKWYP CKWYP
#define CKWXP CKWXP
#define CKWYR CKWYR
#define CKWXR CKWXR
#define CKQC CKQC
#define CKKFKR CKKFKR
#define CKQYP CKQYP
#define CKQXP CKQXP
#define CKQYR CKQYR
#define CKQXR CKQXR
#define CKNU CKNU
#define CKNCF CKNCF
#define CKABE CKABE
#define CKEQC CKEQC
#define CKEQYP CKEQYP
#define CKEQXP CKEQXP
#define CKEQYR CKEQYR
#define CKEQXR CKEQXR
#define DWDOT DWDOT
#define VCKHMS VCKHMS
#define VCKPY VCKPY
#define VCKWYR VCKWYR
#define VCKYTX VCKYTX
#define GET_T_GIVEN_EY GET_T_GIVEN_EY
#elif defined(BL_FORT_USE_LOWERCASE)
#define FORT_CKWYP dsckwyp
#define CKINDX ckindx
#define CKINIT ckinit
#define CKXNUM ckxnum
#define CKSYME cksyme
#define CKSYMS cksyms
#define CKRP ckrp
#define CKPX ckpx
#define CKPY ckpy
#define CKPC ckpc
#define CKRHOX ckrhox
#define CKRHOY ckrhoy
#define CKRHOC ckrhoc
#define CKWT ckwt
#define CKMMWY ckmmwy
#define CKMMWX ckmmwx
#define CKMMWC ckmmwc
#define CKYTX ckytx
#define CKYTCP ckytcp
#define CKYTCR ckytcr
#define CKXTY ckxty
#define CKXTCP ckxtcp
#define CKXTCR ckxtcr
#define CKCTX ckctx
#define CKCTY ckcty
#define CKCPOR ckcpor
#define CKHORT ckhort
#define CKSOR cksor
#define CKCVML ckcvml
#define CKCPML ckcpml
#define CKUML ckuml
#define CKHML ckhml
#define CKGML ckgml
#define CKAML ckaml
#define CKSML cksml
#define CKCVMS ckcvms
#define CKCPMS ckcpms
#define CKUMS ckums
#define CKHMS ckhms
#define CKGMS ckgms
#define CKAMS ckams
#define CKSMS cksms
#define CKCPBL ckcpbl
#define CKCPBS ckcpbs
#define CKCVBL ckcvbl
#define CKCVBS ckcvbs
#define CKHBML ckhbml
#define CKHBMS ckhbms
#define CKUBML ckubml
#define CKUBMS ckubms
#define CKSBML cksbml
#define CKSBMS cksbms
#define CKGBML ckgbml
#define CKGBMS ckgbms
#define CKABML ckabml
#define CKABMS ckabms
#define CKWC ckwc
#define CKWYP ckwyp
#define CKWXP ckwxp
#define CKWYR ckwyr
#define CKWXR ckwxr
#define CKQC ckqc
#define CKKFKR ckkfkr
#define CKQYP ckqyp
#define CKQXP ckqxp
#define CKQYR ckqyr
#define CKQXR ckqxr
#define CKNU cknu
#define CKNCF ckncf
#define CKABE ckabe
#define CKEQC ckeqc
#define CKEQYP ckeqyp
#define CKEQXP ckeqxp
#define CKEQYR ckeqyr
#define CKEQXR ckeqxr
#define DWDOT dwdot
#define VCKHMS vckhms
#define VCKPY vckpy
#define VCKWYR vckwyr
#define VCKYTX vckytx
#define GET_T_GIVEN_EY get_t_given_ey
#elif defined(BL_FORT_USE_UNDERSCORE)
#define FORT_CKWYP dsckwyp_
#define CKINDX ckindx_
#define CKINIT ckinit_
#define CKXNUM ckxnum_
#define CKSYME cksyme_
#define CKSYMS cksyms_
#define CKRP ckrp_
#define CKPX ckpx_
#define CKPY ckpy_
#define CKPC ckpc_
#define CKRHOX ckrhox_
#define CKRHOY ckrhoy_
#define CKRHOC ckrhoc_
#define CKWT ckwt_
#define CKMMWY ckmmwy_
#define CKMMWX ckmmwx_
#define CKMMWC ckmmwc_
#define CKYTX ckytx_
#define CKYTCP ckytcp_
#define CKYTCR ckytcr_
#define CKXTY ckxty_
#define CKXTCP ckxtcp_
#define CKXTCR ckxtcr_
#define CKCTX ckctx_
#define CKCTY ckcty_
#define CKCPOR ckcpor_
#define CKHORT ckhort_
#define CKSOR cksor_
#define CKCVML ckcvml_
#define CKCPML ckcpml_
#define CKUML ckuml_
#define CKHML ckhml_
#define CKGML ckgml_
#define CKAML ckaml_
#define CKSML cksml_
#define CKCVMS ckcvms_
#define CKCPMS ckcpms_
#define CKUMS ckums_
#define CKHMS ckhms_
#define CKGMS ckgms_
#define CKAMS ckams_
#define CKSMS cksms_
#define CKCPBL ckcpbl_
#define CKCPBS ckcpbs_
#define CKCVBL ckcvbl_
#define CKCVBS ckcvbs_
#define CKHBML ckhbml_
#define CKHBMS ckhbms_
#define CKUBML ckubml_
#define CKUBMS ckubms_
#define CKSBML cksbml_
#define CKSBMS cksbms_
#define CKGBML ckgbml_
#define CKGBMS ckgbms_
#define CKABML ckabml_
#define CKABMS ckabms_
#define CKWC ckwc_
#define CKWYP ckwyp_
#define CKWXP ckwxp_
#define CKWYR ckwyr_
#define CKWXR ckwxr_
#define CKQC ckqc_
#define CKKFKR ckkfkr_
#define CKQYP ckqyp_
#define CKQXP ckqxp_
#define CKQYR ckqyr_
#define CKQXR ckqxr_
#define CKNU cknu_
#define CKNCF ckncf_
#define CKABE ckabe_
#define CKEQC ckeqc_
#define CKEQYP ckeqyp_
#define CKEQXP ckeqxp_
#define CKEQYR ckeqyr_
#define CKEQXR ckeqxr_
#define DWDOT dwdot_
#define VCKHMS vckhms_
#define VCKPY vckpy_
#define VCKWYR vckwyr_
#define VCKYTX vckytx_
#define GET_T_GIVEN_EY get_t_given_ey_
#endif

/*function declarations */
void FORT_CKWYP(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  wdot);
void molecularWeight(double * restrict  wt);
void gibbs(double * restrict  species, double * restrict  tc);
void helmholtz(double * restrict  species, double * restrict  tc);
void speciesInternalEnergy(double * restrict  species, double * restrict  tc);
void speciesEnthalpy(double * restrict  species, double * restrict  tc);
void speciesEntropy(double * restrict  species, double * restrict  tc);
void cp_R(double * restrict  species, double * restrict  tc);
void cv_R(double * restrict  species, double * restrict  tc);
void equilibriumConstants(double * restrict  kc, double * restrict  g_RT, double T);
void productionRate(double * restrict  wdot, double * restrict  sc, double T);
void progressRate(double * restrict  qdot, double * restrict  speciesConc, double T);
void progressRateFR(double * restrict  q_f, double * restrict  q_r, double * restrict  speciesConc, double T);
void CKINDX(int * iwrk, double * restrict rwrk, int * mm, int * kk, int * ii, int * nfit );
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double * restrict  rval, int * kerr, int lenline);
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * restrict  rval, int * kerr, int lenline, int lenkray);
void CKSYME(int * kname, int * lenkname);
void CKSYMS(int * kname, int * lenkname);
void CKRP(int * ickwrk, double * restrict  rckwrk, double * restrict  ru, double * restrict  ruc, double * restrict  pa);
void CKPX(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict rwrk, double * restrict  P);
void CKPY(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  P);
void CKPC(double * restrict  rho, double * restrict  T, double * restrict  c, int * iwrk, double * restrict rwrk, double * restrict  P);
void CKRHOX(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict rwrk, double * restrict  rho);
void CKRHOY(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  rho);
void CKRHOC(double * restrict  P, double * restrict  T, double * restrict  c, int * iwrk, double * restrict rwrk, double * restrict  rho);
void CKWT(int * iwrk, double * restrict rwrk, double * restrict  wt);
void CKMMWY(double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  wtm);
void CKMMWX(double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  wtm);
void CKMMWC(double * restrict  c, int * iwrk, double * restrict  rwrk, double * restrict  wtm);
void CKYTX(double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  x);
void CKYTCP(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  c);
void CKYTCR(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  c);
void CKXTY(double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  y);
void CKXTCP(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  c);
void CKXTCR(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  c);
void CKCTX(double * restrict  c, int * iwrk, double * restrict  rwrk, double * restrict  x);
void CKCTY(double * restrict  c, int * iwrk, double * restrict  rwrk, double * restrict  y);
void CKCPOR(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  cpor);
void CKHORT(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  hort);
void CKSOR(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  sor);
void CKCVML(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  cvml);
void CKCPML(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  cvml);
void CKUML(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  uml);
void CKHML(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  uml);
void CKGML(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  gml);
void CKAML(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  aml);
void CKSML(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  sml);
void CKCVMS(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  cvms);
void CKCPMS(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  cvms);
void CKUMS(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  ums);
void CKHMS(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  ums);
void CKGMS(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  gms);
void CKAMS(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  ams);
void CKSMS(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  sms);
void CKCPBL(double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  cpbl);
void CKCPBS(double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  cpbs);
void CKCVBL(double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  cpbl);
void CKCVBS(double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  cpbs);
void CKHBML(double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  hbml);
void CKHBMS(double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  hbms);
void CKUBML(double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  ubml);
void CKUBMS(double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  ubms);
void CKSBML(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  sbml);
void CKSBMS(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  sbms);
void CKGBML(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  gbml);
void CKGBMS(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  gbms);
void CKABML(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  abml);
void CKABMS(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  abms);
void CKWC(double * restrict  T, double * restrict  C, int * iwrk, double * restrict rwrk, double * restrict  wdot);
void CKWYP(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  wdot);
void CKWXP(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict rwrk, double * restrict  wdot);
void CKWYR(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  wdot);
void CKWXR(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict rwrk, double * restrict  wdot);
void CKQC(double * restrict  T, double * restrict  C, int * iwrk, double * restrict rwrk, double * restrict  qdot);
void CKKFKR(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict rwrk, double * restrict  q_f, double * restrict  q_r);
void CKQYP(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  qdot);
void CKQXP(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict rwrk, double * restrict  qdot);
void CKQYR(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  qdot);
void CKQXR(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict rwrk, double * restrict  qdot);
void CKNU(int * kdim, int * iwrk, double * restrict rwrk, int * nuki);
void CKNCF(int * mdim, int * iwrk, double * restrict rwrk, int * ncf);
void CKABE(int * iwrk, double * restrict rwrk, double * restrict  a, double * restrict  b, double * restrict  e );
void CKEQC(double * restrict  T, double * restrict  C , int * iwrk, double * restrict rwrk, double * restrict  eqcon );
void CKEQYP(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  eqcon);
void CKEQXP(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict rwrk, double * restrict  eqcon);
void CKEQYR(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  eqcon);
void CKEQXR(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict rwrk, double * restrict  eqcon);
void DWDOT(double * restrict  J, double * restrict  sc, double * T, int * consP);
void aJacobian(double * restrict J, double * restrict sc, double T, int consP);
void dcvpRdT(double * restrict  species, double * restrict  tc);
void GET_T_GIVEN_EY(double * restrict  e, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  t, int *ierr);
/*vector version */
void vproductionRate(int npt, double * restrict wdot, double * restrict c, double * restrict T);
void VCKHMS(int * restrict np, double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  ums);
void VCKPY(int * restrict np, double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  P);
void VCKWYR(int * restrict np, double * restrict rho, double * restrict T,
            double * restrict y, int * restrict iwrk, double * restrict rwrk,
            double * restrict wdot);
void VCKYTX(int * restrict np, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  x);

/* Inverse molecular weights */
static const double imw[30] = {
    1.0 / 1.007970,  /*H */
    1.0 / 2.015940,  /*H2 */
    1.0 / 15.035060,  /*CH3 */
    1.0 / 15.999400,  /*O */
    1.0 / 16.043030,  /*CH4 */
    1.0 / 17.007370,  /*OH */
    1.0 / 18.015340,  /*H2O */
    1.0 / 26.038240,  /*C2H2 */
    1.0 / 28.010550,  /*CO */
    1.0 / 28.054180,  /*C2H4 */
    1.0 / 29.062150,  /*C2H5 */
    1.0 / 30.026490,  /*CH2O */
    1.0 / 30.070120,  /*C2H6 */
    1.0 / 31.034460,  /*CH3O */
    1.0 / 31.998800,  /*O2 */
    1.0 / 33.006770,  /*HO2 */
    1.0 / 34.014740,  /*H2O2 */
    1.0 / 44.009950,  /*CO2 */
    1.0 / 44.053580,  /*CH3HCO */
    1.0 / 46.025890,  /*HCOOH */
    1.0 / 46.069520,  /*CH3OCH3 */
    1.0 / 59.045010,  /*CH3OCO */
    1.0 / 60.052980,  /*CH3OCHO */
    1.0 / 62.068920,  /*CH3OCH2OH */
    1.0 / 75.044410,  /*OCH2OCHO */
    1.0 / 75.044410,  /*HOCH2OCO */
    1.0 / 77.060350,  /*CH3OCH2O2 */
    1.0 / 92.051780,  /*HO2CH2OCHO */
    1.0 / 109.059150,  /*O2CH2OCH2O2H */
    1.0 / 28.013400};  /*N2 */


/*A few mechanism parameters */
void CKINDX(int * iwrk, double * restrict  rwrk, int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 4;
    *kk = 30;
    *ii = 3;
    *nfit = -1; /*Why do you need this anyway ?  */
}


/* ckxnum... for parsing strings  */
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double * restrict  rval, int * kerr, int lenline )
{
    int n,i; /*Loop Counters */
    char *p; /*String Tokens */
    char cstr[1000];
    /* Strip Comments  */
    for (i=0; i<lenline; ++i) {
        if (line[i]=='!') {
            cstr[i] = '\0';
            break;
        }
        cstr[i] = line[i];
    }

    p = strtok(cstr," ");
    if (!p) {
        *nval = 0;
        *kerr = 1;
        return;
    }
    for (n=0; n<*nexp; ++n) {
        rval[n] = atof(p);
        p = strtok(NULL, " ");
        if (!p) break;
    }
    *nval = n+1;
    if (*nval < *nexp) *kerr = 1;
    return;
}


/* cksnum... for parsing strings  */
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * restrict  rval, int * kerr, int lenline, int lenkray)
{
    /*Not done yet ... */
}


/* Returns the char strings of element names */
void CKSYME(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*4; i++) {
        kname[i] = ' ';
    }

    /* C  */
    kname[ 0*lenkname + 0 ] = 'C';
    kname[ 0*lenkname + 1 ] = ' ';

    /* H  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = ' ';

    /* O  */
    kname[ 2*lenkname + 0 ] = 'O';
    kname[ 2*lenkname + 1 ] = ' ';

    /* N  */
    kname[ 3*lenkname + 0 ] = 'N';
    kname[ 3*lenkname + 1 ] = ' ';

}


/* Returns the char strings of species names */
void CKSYMS(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*30; i++) {
        kname[i] = ' ';
    }

    /* H  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = ' ';

    /* H2  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = '2';
    kname[ 1*lenkname + 2 ] = ' ';

    /* CH3  */
    kname[ 2*lenkname + 0 ] = 'C';
    kname[ 2*lenkname + 1 ] = 'H';
    kname[ 2*lenkname + 2 ] = '3';
    kname[ 2*lenkname + 3 ] = ' ';

    /* O  */
    kname[ 3*lenkname + 0 ] = 'O';
    kname[ 3*lenkname + 1 ] = ' ';

    /* CH4  */
    kname[ 4*lenkname + 0 ] = 'C';
    kname[ 4*lenkname + 1 ] = 'H';
    kname[ 4*lenkname + 2 ] = '4';
    kname[ 4*lenkname + 3 ] = ' ';

    /* OH  */
    kname[ 5*lenkname + 0 ] = 'O';
    kname[ 5*lenkname + 1 ] = 'H';
    kname[ 5*lenkname + 2 ] = ' ';

    /* H2O  */
    kname[ 6*lenkname + 0 ] = 'H';
    kname[ 6*lenkname + 1 ] = '2';
    kname[ 6*lenkname + 2 ] = 'O';
    kname[ 6*lenkname + 3 ] = ' ';

    /* C2H2  */
    kname[ 7*lenkname + 0 ] = 'C';
    kname[ 7*lenkname + 1 ] = '2';
    kname[ 7*lenkname + 2 ] = 'H';
    kname[ 7*lenkname + 3 ] = '2';
    kname[ 7*lenkname + 4 ] = ' ';

    /* CO  */
    kname[ 8*lenkname + 0 ] = 'C';
    kname[ 8*lenkname + 1 ] = 'O';
    kname[ 8*lenkname + 2 ] = ' ';

    /* C2H4  */
    kname[ 9*lenkname + 0 ] = 'C';
    kname[ 9*lenkname + 1 ] = '2';
    kname[ 9*lenkname + 2 ] = 'H';
    kname[ 9*lenkname + 3 ] = '4';
    kname[ 9*lenkname + 4 ] = ' ';

    /* C2H5  */
    kname[ 10*lenkname + 0 ] = 'C';
    kname[ 10*lenkname + 1 ] = '2';
    kname[ 10*lenkname + 2 ] = 'H';
    kname[ 10*lenkname + 3 ] = '5';
    kname[ 10*lenkname + 4 ] = ' ';

    /* CH2O  */
    kname[ 11*lenkname + 0 ] = 'C';
    kname[ 11*lenkname + 1 ] = 'H';
    kname[ 11*lenkname + 2 ] = '2';
    kname[ 11*lenkname + 3 ] = 'O';
    kname[ 11*lenkname + 4 ] = ' ';

    /* C2H6  */
    kname[ 12*lenkname + 0 ] = 'C';
    kname[ 12*lenkname + 1 ] = '2';
    kname[ 12*lenkname + 2 ] = 'H';
    kname[ 12*lenkname + 3 ] = '6';
    kname[ 12*lenkname + 4 ] = ' ';

    /* CH3O  */
    kname[ 13*lenkname + 0 ] = 'C';
    kname[ 13*lenkname + 1 ] = 'H';
    kname[ 13*lenkname + 2 ] = '3';
    kname[ 13*lenkname + 3 ] = 'O';
    kname[ 13*lenkname + 4 ] = ' ';

    /* O2  */
    kname[ 14*lenkname + 0 ] = 'O';
    kname[ 14*lenkname + 1 ] = '2';
    kname[ 14*lenkname + 2 ] = ' ';

    /* HO2  */
    kname[ 15*lenkname + 0 ] = 'H';
    kname[ 15*lenkname + 1 ] = 'O';
    kname[ 15*lenkname + 2 ] = '2';
    kname[ 15*lenkname + 3 ] = ' ';

    /* H2O2  */
    kname[ 16*lenkname + 0 ] = 'H';
    kname[ 16*lenkname + 1 ] = '2';
    kname[ 16*lenkname + 2 ] = 'O';
    kname[ 16*lenkname + 3 ] = '2';
    kname[ 16*lenkname + 4 ] = ' ';

    /* CO2  */
    kname[ 17*lenkname + 0 ] = 'C';
    kname[ 17*lenkname + 1 ] = 'O';
    kname[ 17*lenkname + 2 ] = '2';
    kname[ 17*lenkname + 3 ] = ' ';

    /* CH3HCO  */
    kname[ 18*lenkname + 0 ] = 'C';
    kname[ 18*lenkname + 1 ] = 'H';
    kname[ 18*lenkname + 2 ] = '3';
    kname[ 18*lenkname + 3 ] = 'H';
    kname[ 18*lenkname + 4 ] = 'C';
    kname[ 18*lenkname + 5 ] = 'O';
    kname[ 18*lenkname + 6 ] = ' ';

    /* HCOOH  */
    kname[ 19*lenkname + 0 ] = 'H';
    kname[ 19*lenkname + 1 ] = 'C';
    kname[ 19*lenkname + 2 ] = 'O';
    kname[ 19*lenkname + 3 ] = 'O';
    kname[ 19*lenkname + 4 ] = 'H';
    kname[ 19*lenkname + 5 ] = ' ';

    /* CH3OCH3  */
    kname[ 20*lenkname + 0 ] = 'C';
    kname[ 20*lenkname + 1 ] = 'H';
    kname[ 20*lenkname + 2 ] = '3';
    kname[ 20*lenkname + 3 ] = 'O';
    kname[ 20*lenkname + 4 ] = 'C';
    kname[ 20*lenkname + 5 ] = 'H';
    kname[ 20*lenkname + 6 ] = '3';
    kname[ 20*lenkname + 7 ] = ' ';

    /* CH3OCO  */
    kname[ 21*lenkname + 0 ] = 'C';
    kname[ 21*lenkname + 1 ] = 'H';
    kname[ 21*lenkname + 2 ] = '3';
    kname[ 21*lenkname + 3 ] = 'O';
    kname[ 21*lenkname + 4 ] = 'C';
    kname[ 21*lenkname + 5 ] = 'O';
    kname[ 21*lenkname + 6 ] = ' ';

    /* CH3OCHO  */
    kname[ 22*lenkname + 0 ] = 'C';
    kname[ 22*lenkname + 1 ] = 'H';
    kname[ 22*lenkname + 2 ] = '3';
    kname[ 22*lenkname + 3 ] = 'O';
    kname[ 22*lenkname + 4 ] = 'C';
    kname[ 22*lenkname + 5 ] = 'H';
    kname[ 22*lenkname + 6 ] = 'O';
    kname[ 22*lenkname + 7 ] = ' ';

    /* CH3OCH2OH  */
    kname[ 23*lenkname + 0 ] = 'C';
    kname[ 23*lenkname + 1 ] = 'H';
    kname[ 23*lenkname + 2 ] = '3';
    kname[ 23*lenkname + 3 ] = 'O';
    kname[ 23*lenkname + 4 ] = 'C';
    kname[ 23*lenkname + 5 ] = 'H';
    kname[ 23*lenkname + 6 ] = '2';
    kname[ 23*lenkname + 7 ] = 'O';
    kname[ 23*lenkname + 8 ] = 'H';
    kname[ 23*lenkname + 9 ] = ' ';

    /* OCH2OCHO  */
    kname[ 24*lenkname + 0 ] = 'O';
    kname[ 24*lenkname + 1 ] = 'C';
    kname[ 24*lenkname + 2 ] = 'H';
    kname[ 24*lenkname + 3 ] = '2';
    kname[ 24*lenkname + 4 ] = 'O';
    kname[ 24*lenkname + 5 ] = 'C';
    kname[ 24*lenkname + 6 ] = 'H';
    kname[ 24*lenkname + 7 ] = 'O';
    kname[ 24*lenkname + 8 ] = ' ';

    /* HOCH2OCO  */
    kname[ 25*lenkname + 0 ] = 'H';
    kname[ 25*lenkname + 1 ] = 'O';
    kname[ 25*lenkname + 2 ] = 'C';
    kname[ 25*lenkname + 3 ] = 'H';
    kname[ 25*lenkname + 4 ] = '2';
    kname[ 25*lenkname + 5 ] = 'O';
    kname[ 25*lenkname + 6 ] = 'C';
    kname[ 25*lenkname + 7 ] = 'O';
    kname[ 25*lenkname + 8 ] = ' ';

    /* CH3OCH2O2  */
    kname[ 26*lenkname + 0 ] = 'C';
    kname[ 26*lenkname + 1 ] = 'H';
    kname[ 26*lenkname + 2 ] = '3';
    kname[ 26*lenkname + 3 ] = 'O';
    kname[ 26*lenkname + 4 ] = 'C';
    kname[ 26*lenkname + 5 ] = 'H';
    kname[ 26*lenkname + 6 ] = '2';
    kname[ 26*lenkname + 7 ] = 'O';
    kname[ 26*lenkname + 8 ] = '2';
    kname[ 26*lenkname + 9 ] = ' ';

    /* HO2CH2OCHO  */
    kname[ 27*lenkname + 0 ] = 'H';
    kname[ 27*lenkname + 1 ] = 'O';
    kname[ 27*lenkname + 2 ] = '2';
    kname[ 27*lenkname + 3 ] = 'C';
    kname[ 27*lenkname + 4 ] = 'H';
    kname[ 27*lenkname + 5 ] = '2';
    kname[ 27*lenkname + 6 ] = 'O';
    kname[ 27*lenkname + 7 ] = 'C';
    kname[ 27*lenkname + 8 ] = 'H';
    kname[ 27*lenkname + 9 ] = 'O';
    kname[ 27*lenkname + 10 ] = ' ';

    /* O2CH2OCH2O2H  */
    kname[ 28*lenkname + 0 ] = 'O';
    kname[ 28*lenkname + 1 ] = '2';
    kname[ 28*lenkname + 2 ] = 'C';
    kname[ 28*lenkname + 3 ] = 'H';
    kname[ 28*lenkname + 4 ] = '2';
    kname[ 28*lenkname + 5 ] = 'O';
    kname[ 28*lenkname + 6 ] = 'C';
    kname[ 28*lenkname + 7 ] = 'H';
    kname[ 28*lenkname + 8 ] = '2';
    kname[ 28*lenkname + 9 ] = 'O';
    kname[ 28*lenkname + 10 ] = '2';
    kname[ 28*lenkname + 11 ] = 'H';
    kname[ 28*lenkname + 12 ] = ' ';

    /* N2  */
    kname[ 29*lenkname + 0 ] = 'N';
    kname[ 29*lenkname + 1 ] = '2';
    kname[ 29*lenkname + 2 ] = ' ';

}


/* Returns R, Rc, Patm */
void CKRP(int * ickwrk, double * restrict  rckwrk, double * restrict  ru, double * restrict  ruc, double * restrict  pa)
{
     *ru  = 8.31451e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}


/*Compute P = rhoRT/W(x) */
void CKPX(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  P)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.035060; /*CH3 */
    XW += x[3]*15.999400; /*O */
    XW += x[4]*16.043030; /*CH4 */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*26.038240; /*C2H2 */
    XW += x[8]*28.010550; /*CO */
    XW += x[9]*28.054180; /*C2H4 */
    XW += x[10]*29.062150; /*C2H5 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*30.070120; /*C2H6 */
    XW += x[13]*31.034460; /*CH3O */
    XW += x[14]*31.998800; /*O2 */
    XW += x[15]*33.006770; /*HO2 */
    XW += x[16]*34.014740; /*H2O2 */
    XW += x[17]*44.009950; /*CO2 */
    XW += x[18]*44.053580; /*CH3HCO */
    XW += x[19]*46.025890; /*HCOOH */
    XW += x[20]*46.069520; /*CH3OCH3 */
    XW += x[21]*59.045010; /*CH3OCO */
    XW += x[22]*60.052980; /*CH3OCHO */
    XW += x[23]*62.068920; /*CH3OCH2OH */
    XW += x[24]*75.044410; /*OCH2OCHO */
    XW += x[25]*75.044410; /*HOCH2OCO */
    XW += x[26]*77.060350; /*CH3OCH2O2 */
    XW += x[27]*92.051780; /*HO2CH2OCHO */
    XW += x[28]*109.059150; /*O2CH2OCH2O2H */
    XW += x[29]*28.013400; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void CKPY(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*CH3 */
    YOW += y[3]*imw[3]; /*O */
    YOW += y[4]*imw[4]; /*CH4 */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*H2O */
    YOW += y[7]*imw[7]; /*C2H2 */
    YOW += y[8]*imw[8]; /*CO */
    YOW += y[9]*imw[9]; /*C2H4 */
    YOW += y[10]*imw[10]; /*C2H5 */
    YOW += y[11]*imw[11]; /*CH2O */
    YOW += y[12]*imw[12]; /*C2H6 */
    YOW += y[13]*imw[13]; /*CH3O */
    YOW += y[14]*imw[14]; /*O2 */
    YOW += y[15]*imw[15]; /*HO2 */
    YOW += y[16]*imw[16]; /*H2O2 */
    YOW += y[17]*imw[17]; /*CO2 */
    YOW += y[18]*imw[18]; /*CH3HCO */
    YOW += y[19]*imw[19]; /*HCOOH */
    YOW += y[20]*imw[20]; /*CH3OCH3 */
    YOW += y[21]*imw[21]; /*CH3OCO */
    YOW += y[22]*imw[22]; /*CH3OCHO */
    YOW += y[23]*imw[23]; /*CH3OCH2OH */
    YOW += y[24]*imw[24]; /*OCH2OCHO */
    YOW += y[25]*imw[25]; /*HOCH2OCO */
    YOW += y[26]*imw[26]; /*CH3OCH2O2 */
    YOW += y[27]*imw[27]; /*HO2CH2OCHO */
    YOW += y[28]*imw[28]; /*O2CH2OCH2O2H */
    YOW += y[29]*imw[29]; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void VCKPY(int * restrict np, double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  P)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<30; n++) {
        for (int i=0; i<(*np); i++) {
            YOW[i] += y[n*(*np)+i] * imw[n];
        }
    }

    for (int i=0; i<(*np); i++) {
        P[i] = rho[i] * 8.31451e+07 * T[i] * YOW[i]; /*P = rho*R*T/W */
    }

    return;
}


/*Compute P = rhoRT/W(c) */
void CKPC(double * restrict  rho, double * restrict  T, double * restrict  c, int * iwrk, double * restrict  rwrk, double * restrict  P)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*1.007970; /*H */
    W += c[1]*2.015940; /*H2 */
    W += c[2]*15.035060; /*CH3 */
    W += c[3]*15.999400; /*O */
    W += c[4]*16.043030; /*CH4 */
    W += c[5]*17.007370; /*OH */
    W += c[6]*18.015340; /*H2O */
    W += c[7]*26.038240; /*C2H2 */
    W += c[8]*28.010550; /*CO */
    W += c[9]*28.054180; /*C2H4 */
    W += c[10]*29.062150; /*C2H5 */
    W += c[11]*30.026490; /*CH2O */
    W += c[12]*30.070120; /*C2H6 */
    W += c[13]*31.034460; /*CH3O */
    W += c[14]*31.998800; /*O2 */
    W += c[15]*33.006770; /*HO2 */
    W += c[16]*34.014740; /*H2O2 */
    W += c[17]*44.009950; /*CO2 */
    W += c[18]*44.053580; /*CH3HCO */
    W += c[19]*46.025890; /*HCOOH */
    W += c[20]*46.069520; /*CH3OCH3 */
    W += c[21]*59.045010; /*CH3OCO */
    W += c[22]*60.052980; /*CH3OCHO */
    W += c[23]*62.068920; /*CH3OCH2OH */
    W += c[24]*75.044410; /*OCH2OCHO */
    W += c[25]*75.044410; /*HOCH2OCO */
    W += c[26]*77.060350; /*CH3OCH2O2 */
    W += c[27]*92.051780; /*HO2CH2OCHO */
    W += c[28]*109.059150; /*O2CH2OCH2O2H */
    W += c[29]*28.013400; /*N2 */

    for (id = 0; id < 30; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31451e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.035060; /*CH3 */
    XW += x[3]*15.999400; /*O */
    XW += x[4]*16.043030; /*CH4 */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*26.038240; /*C2H2 */
    XW += x[8]*28.010550; /*CO */
    XW += x[9]*28.054180; /*C2H4 */
    XW += x[10]*29.062150; /*C2H5 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*30.070120; /*C2H6 */
    XW += x[13]*31.034460; /*CH3O */
    XW += x[14]*31.998800; /*O2 */
    XW += x[15]*33.006770; /*HO2 */
    XW += x[16]*34.014740; /*H2O2 */
    XW += x[17]*44.009950; /*CO2 */
    XW += x[18]*44.053580; /*CH3HCO */
    XW += x[19]*46.025890; /*HCOOH */
    XW += x[20]*46.069520; /*CH3OCH3 */
    XW += x[21]*59.045010; /*CH3OCO */
    XW += x[22]*60.052980; /*CH3OCHO */
    XW += x[23]*62.068920; /*CH3OCH2OH */
    XW += x[24]*75.044410; /*OCH2OCHO */
    XW += x[25]*75.044410; /*HOCH2OCO */
    XW += x[26]*77.060350; /*CH3OCH2O2 */
    XW += x[27]*92.051780; /*HO2CH2OCHO */
    XW += x[28]*109.059150; /*O2CH2OCH2O2H */
    XW += x[29]*28.013400; /*N2 */
    *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  rho)
{
    double YOW = 0;
    double tmp[30];

    for (int i = 0; i < 30; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 30; i++)
    {
        YOW += tmp[i];
    }

    *rho = *P / (83145100 * (*T) * YOW);/*rho = P*W/(R*T) */
    return;
}


/*Compute rho = P*W(c)/(R*T) */
void CKRHOC(double * restrict  P, double * restrict  T, double * restrict  c, int * iwrk, double * restrict  rwrk, double * restrict  rho)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*1.007970; /*H */
    W += c[1]*2.015940; /*H2 */
    W += c[2]*15.035060; /*CH3 */
    W += c[3]*15.999400; /*O */
    W += c[4]*16.043030; /*CH4 */
    W += c[5]*17.007370; /*OH */
    W += c[6]*18.015340; /*H2O */
    W += c[7]*26.038240; /*C2H2 */
    W += c[8]*28.010550; /*CO */
    W += c[9]*28.054180; /*C2H4 */
    W += c[10]*29.062150; /*C2H5 */
    W += c[11]*30.026490; /*CH2O */
    W += c[12]*30.070120; /*C2H6 */
    W += c[13]*31.034460; /*CH3O */
    W += c[14]*31.998800; /*O2 */
    W += c[15]*33.006770; /*HO2 */
    W += c[16]*34.014740; /*H2O2 */
    W += c[17]*44.009950; /*CO2 */
    W += c[18]*44.053580; /*CH3HCO */
    W += c[19]*46.025890; /*HCOOH */
    W += c[20]*46.069520; /*CH3OCH3 */
    W += c[21]*59.045010; /*CH3OCO */
    W += c[22]*60.052980; /*CH3OCHO */
    W += c[23]*62.068920; /*CH3OCH2OH */
    W += c[24]*75.044410; /*OCH2OCHO */
    W += c[25]*75.044410; /*HOCH2OCO */
    W += c[26]*77.060350; /*CH3OCH2O2 */
    W += c[27]*92.051780; /*HO2CH2OCHO */
    W += c[28]*109.059150; /*O2CH2OCH2O2H */
    W += c[29]*28.013400; /*N2 */

    for (id = 0; id < 30; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.31451e+07); /*rho = PW/(R*T) */

    return;
}


/*get molecular weight for all species */
void CKWT(int * iwrk, double * restrict  rwrk, double * restrict  wt)
{
    molecularWeight(wt);
}


/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWY(double * restrict y, int * iwrk, double * restrict  rwrk, double * restrict  wtm)
{
    double YOW = 0;
    double tmp[30];

    for (int i = 0; i < 30; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 30; i++)
    {
        YOW += tmp[i];
    }

    *wtm = 1.0 / YOW;
    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWX(double * restrict x, int * iwrk, double * restrict  rwrk, double * restrict  wtm)
{
    double XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.035060; /*CH3 */
    XW += x[3]*15.999400; /*O */
    XW += x[4]*16.043030; /*CH4 */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*26.038240; /*C2H2 */
    XW += x[8]*28.010550; /*CO */
    XW += x[9]*28.054180; /*C2H4 */
    XW += x[10]*29.062150; /*C2H5 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*30.070120; /*C2H6 */
    XW += x[13]*31.034460; /*CH3O */
    XW += x[14]*31.998800; /*O2 */
    XW += x[15]*33.006770; /*HO2 */
    XW += x[16]*34.014740; /*H2O2 */
    XW += x[17]*44.009950; /*CO2 */
    XW += x[18]*44.053580; /*CH3HCO */
    XW += x[19]*46.025890; /*HCOOH */
    XW += x[20]*46.069520; /*CH3OCH3 */
    XW += x[21]*59.045010; /*CH3OCO */
    XW += x[22]*60.052980; /*CH3OCHO */
    XW += x[23]*62.068920; /*CH3OCH2OH */
    XW += x[24]*75.044410; /*OCH2OCHO */
    XW += x[25]*75.044410; /*HOCH2OCO */
    XW += x[26]*77.060350; /*CH3OCH2O2 */
    XW += x[27]*92.051780; /*HO2CH2OCHO */
    XW += x[28]*109.059150; /*O2CH2OCH2O2H */
    XW += x[29]*28.013400; /*N2 */
    *wtm = XW;

    return;
}


/*given c[species]: molar concentration */
/*returns mean molecular weight (gm/mole) */
void CKMMWC(double * restrict c, int * iwrk, double * restrict  rwrk, double * restrict  wtm)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*1.007970; /*H */
    W += c[1]*2.015940; /*H2 */
    W += c[2]*15.035060; /*CH3 */
    W += c[3]*15.999400; /*O */
    W += c[4]*16.043030; /*CH4 */
    W += c[5]*17.007370; /*OH */
    W += c[6]*18.015340; /*H2O */
    W += c[7]*26.038240; /*C2H2 */
    W += c[8]*28.010550; /*CO */
    W += c[9]*28.054180; /*C2H4 */
    W += c[10]*29.062150; /*C2H5 */
    W += c[11]*30.026490; /*CH2O */
    W += c[12]*30.070120; /*C2H6 */
    W += c[13]*31.034460; /*CH3O */
    W += c[14]*31.998800; /*O2 */
    W += c[15]*33.006770; /*HO2 */
    W += c[16]*34.014740; /*H2O2 */
    W += c[17]*44.009950; /*CO2 */
    W += c[18]*44.053580; /*CH3HCO */
    W += c[19]*46.025890; /*HCOOH */
    W += c[20]*46.069520; /*CH3OCH3 */
    W += c[21]*59.045010; /*CH3OCO */
    W += c[22]*60.052980; /*CH3OCHO */
    W += c[23]*62.068920; /*CH3OCH2OH */
    W += c[24]*75.044410; /*OCH2OCHO */
    W += c[25]*75.044410; /*HOCH2OCO */
    W += c[26]*77.060350; /*CH3OCH2O2 */
    W += c[27]*92.051780; /*HO2CH2OCHO */
    W += c[28]*109.059150; /*O2CH2OCH2O2H */
    W += c[29]*28.013400; /*N2 */

    for (id = 0; id < 30; ++id) {
        sumC += c[id];
    }
    /* CK provides no guard against divison by zero */
    *wtm = W/sumC;

    return;
}


/*convert y[species] (mass fracs) to x[species] (mole fracs) */
void CKYTX(double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  x)
{
    double YOW = 0;
    double tmp[30];

    for (int i = 0; i < 30; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 30; i++)
    {
        YOW += tmp[i];
    }

    double YOWINV = 1.0/YOW;

    for (int i = 0; i < 30; i++)
    {
        x[i] = y[i]*imw[i]*YOWINV;
    }
    return;
}


/*convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs) */
void VCKYTX(int * restrict np, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  x)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<30; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<30; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCP(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  c)
{
    double YOW = 0;
    double PWORT;

    /*Compute inverse of mean molecular wt first */
    for (int i = 0; i < 30; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 30; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 30; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCR(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  c)
{
    for (int i = 0; i < 30; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
void CKXTY(double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.035060; /*CH3 */
    XW += x[3]*15.999400; /*O */
    XW += x[4]*16.043030; /*CH4 */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*26.038240; /*C2H2 */
    XW += x[8]*28.010550; /*CO */
    XW += x[9]*28.054180; /*C2H4 */
    XW += x[10]*29.062150; /*C2H5 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*30.070120; /*C2H6 */
    XW += x[13]*31.034460; /*CH3O */
    XW += x[14]*31.998800; /*O2 */
    XW += x[15]*33.006770; /*HO2 */
    XW += x[16]*34.014740; /*H2O2 */
    XW += x[17]*44.009950; /*CO2 */
    XW += x[18]*44.053580; /*CH3HCO */
    XW += x[19]*46.025890; /*HCOOH */
    XW += x[20]*46.069520; /*CH3OCH3 */
    XW += x[21]*59.045010; /*CH3OCO */
    XW += x[22]*60.052980; /*CH3OCHO */
    XW += x[23]*62.068920; /*CH3OCH2OH */
    XW += x[24]*75.044410; /*OCH2OCHO */
    XW += x[25]*75.044410; /*HOCH2OCO */
    XW += x[26]*77.060350; /*CH3OCH2O2 */
    XW += x[27]*92.051780; /*HO2CH2OCHO */
    XW += x[28]*109.059150; /*O2CH2OCH2O2H */
    XW += x[29]*28.013400; /*N2 */
    /*Now compute conversion */
    y[0] = x[0]*1.007970/XW; 
    y[1] = x[1]*2.015940/XW; 
    y[2] = x[2]*15.035060/XW; 
    y[3] = x[3]*15.999400/XW; 
    y[4] = x[4]*16.043030/XW; 
    y[5] = x[5]*17.007370/XW; 
    y[6] = x[6]*18.015340/XW; 
    y[7] = x[7]*26.038240/XW; 
    y[8] = x[8]*28.010550/XW; 
    y[9] = x[9]*28.054180/XW; 
    y[10] = x[10]*29.062150/XW; 
    y[11] = x[11]*30.026490/XW; 
    y[12] = x[12]*30.070120/XW; 
    y[13] = x[13]*31.034460/XW; 
    y[14] = x[14]*31.998800/XW; 
    y[15] = x[15]*33.006770/XW; 
    y[16] = x[16]*34.014740/XW; 
    y[17] = x[17]*44.009950/XW; 
    y[18] = x[18]*44.053580/XW; 
    y[19] = x[19]*46.025890/XW; 
    y[20] = x[20]*46.069520/XW; 
    y[21] = x[21]*59.045010/XW; 
    y[22] = x[22]*60.052980/XW; 
    y[23] = x[23]*62.068920/XW; 
    y[24] = x[24]*75.044410/XW; 
    y[25] = x[25]*75.044410/XW; 
    y[26] = x[26]*77.060350/XW; 
    y[27] = x[27]*92.051780/XW; 
    y[28] = x[28]*109.059150/XW; 
    y[29] = x[29]*28.013400/XW; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 30; ++id) {
        c[id] = x[id]*PORT;
    }

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCR(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  c)
{
    int id; /*loop counter */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.035060; /*CH3 */
    XW += x[3]*15.999400; /*O */
    XW += x[4]*16.043030; /*CH4 */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*26.038240; /*C2H2 */
    XW += x[8]*28.010550; /*CO */
    XW += x[9]*28.054180; /*C2H4 */
    XW += x[10]*29.062150; /*C2H5 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*30.070120; /*C2H6 */
    XW += x[13]*31.034460; /*CH3O */
    XW += x[14]*31.998800; /*O2 */
    XW += x[15]*33.006770; /*HO2 */
    XW += x[16]*34.014740; /*H2O2 */
    XW += x[17]*44.009950; /*CO2 */
    XW += x[18]*44.053580; /*CH3HCO */
    XW += x[19]*46.025890; /*HCOOH */
    XW += x[20]*46.069520; /*CH3OCH3 */
    XW += x[21]*59.045010; /*CH3OCO */
    XW += x[22]*60.052980; /*CH3OCHO */
    XW += x[23]*62.068920; /*CH3OCH2OH */
    XW += x[24]*75.044410; /*OCH2OCHO */
    XW += x[25]*75.044410; /*HOCH2OCO */
    XW += x[26]*77.060350; /*CH3OCH2O2 */
    XW += x[27]*92.051780; /*HO2CH2OCHO */
    XW += x[28]*109.059150; /*O2CH2OCH2O2H */
    XW += x[29]*28.013400; /*N2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 30; ++id) {
        c[id] = x[id]*ROW;
    }

    return;
}


/*convert c[species] (molar conc) to x[species] (mole fracs) */
void CKCTX(double * restrict  c, int * iwrk, double * restrict  rwrk, double * restrict  x)
{
    int id; /*loop counter */
    double sumC = 0; 

    /*compute sum of c  */
    for (id = 0; id < 30; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    for (id = 0; id < 30; ++id) {
        x[id] = c[id]/sumC;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
void CKCTY(double * restrict  c, int * iwrk, double * restrict  rwrk, double * restrict  y)
{
    double CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*1.007970; /*H */
    CW += c[1]*2.015940; /*H2 */
    CW += c[2]*15.035060; /*CH3 */
    CW += c[3]*15.999400; /*O */
    CW += c[4]*16.043030; /*CH4 */
    CW += c[5]*17.007370; /*OH */
    CW += c[6]*18.015340; /*H2O */
    CW += c[7]*26.038240; /*C2H2 */
    CW += c[8]*28.010550; /*CO */
    CW += c[9]*28.054180; /*C2H4 */
    CW += c[10]*29.062150; /*C2H5 */
    CW += c[11]*30.026490; /*CH2O */
    CW += c[12]*30.070120; /*C2H6 */
    CW += c[13]*31.034460; /*CH3O */
    CW += c[14]*31.998800; /*O2 */
    CW += c[15]*33.006770; /*HO2 */
    CW += c[16]*34.014740; /*H2O2 */
    CW += c[17]*44.009950; /*CO2 */
    CW += c[18]*44.053580; /*CH3HCO */
    CW += c[19]*46.025890; /*HCOOH */
    CW += c[20]*46.069520; /*CH3OCH3 */
    CW += c[21]*59.045010; /*CH3OCO */
    CW += c[22]*60.052980; /*CH3OCHO */
    CW += c[23]*62.068920; /*CH3OCH2OH */
    CW += c[24]*75.044410; /*OCH2OCHO */
    CW += c[25]*75.044410; /*HOCH2OCO */
    CW += c[26]*77.060350; /*CH3OCH2O2 */
    CW += c[27]*92.051780; /*HO2CH2OCHO */
    CW += c[28]*109.059150; /*O2CH2OCH2O2H */
    CW += c[29]*28.013400; /*N2 */
    /*Now compute conversion */
    y[0] = c[0]*1.007970/CW; 
    y[1] = c[1]*2.015940/CW; 
    y[2] = c[2]*15.035060/CW; 
    y[3] = c[3]*15.999400/CW; 
    y[4] = c[4]*16.043030/CW; 
    y[5] = c[5]*17.007370/CW; 
    y[6] = c[6]*18.015340/CW; 
    y[7] = c[7]*26.038240/CW; 
    y[8] = c[8]*28.010550/CW; 
    y[9] = c[9]*28.054180/CW; 
    y[10] = c[10]*29.062150/CW; 
    y[11] = c[11]*30.026490/CW; 
    y[12] = c[12]*30.070120/CW; 
    y[13] = c[13]*31.034460/CW; 
    y[14] = c[14]*31.998800/CW; 
    y[15] = c[15]*33.006770/CW; 
    y[16] = c[16]*34.014740/CW; 
    y[17] = c[17]*44.009950/CW; 
    y[18] = c[18]*44.053580/CW; 
    y[19] = c[19]*46.025890/CW; 
    y[20] = c[20]*46.069520/CW; 
    y[21] = c[21]*59.045010/CW; 
    y[22] = c[22]*60.052980/CW; 
    y[23] = c[23]*62.068920/CW; 
    y[24] = c[24]*75.044410/CW; 
    y[25] = c[25]*75.044410/CW; 
    y[26] = c[26]*77.060350/CW; 
    y[27] = c[27]*92.051780/CW; 
    y[28] = c[28]*109.059150/CW; 
    y[29] = c[29]*28.013400/CW; 

    return;
}


/*get Cp/R as a function of T  */
/*for all species (Eq 19) */
void CKCPOR(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  cpor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpor, tc);
}


/*get H/RT as a function of T  */
/*for all species (Eq 20) */
void CKHORT(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  hort)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEnthalpy(hort, tc);
}


/*get S/R as a function of T  */
/*for all species (Eq 21) */
void CKSOR(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  sor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sor, tc);
}


/*get specific heat at constant volume as a function  */
/*of T for all species (molar units) */
void CKCVML(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  cvml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 30; ++id) {
        cvml[id] *= 8.31451e+07;
    }
}


/*get specific heat at constant pressure as a  */
/*function of T for all species (molar units) */
void CKCPML(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  cpml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 30; ++id) {
        cpml[id] *= 8.31451e+07;
    }
}


/*get internal energy as a function  */
/*of T for all species (molar units) */
void CKUML(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  uml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 30; ++id) {
        uml[id] *= RT;
    }
}


/*get enthalpy as a function  */
/*of T for all species (molar units) */
void CKHML(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  hml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 30; ++id) {
        hml[id] *= RT;
    }
}


/*get standard-state Gibbs energy as a function  */
/*of T for all species (molar units) */
void CKGML(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  gml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 30; ++id) {
        gml[id] *= RT;
    }
}


/*get standard-state Helmholtz free energy as a  */
/*function of T for all species (molar units) */
void CKAML(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  aml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(aml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 30; ++id) {
        aml[id] *= RT;
    }
}


/*Returns the standard-state entropies in molar units */
void CKSML(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  sml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 30; ++id) {
        sml[id] *= 8.31451e+07;
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
void CKCVMS(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  cvms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 8.248767324424338e+07; /*H */
    cvms[1] *= 4.124383662212169e+07; /*H2 */
    cvms[2] *= 5.530081023953346e+06; /*CH3 */
    cvms[3] *= 5.196763628636074e+06; /*O */
    cvms[4] *= 5.182630712527496e+06; /*CH4 */
    cvms[5] *= 4.888768810227566e+06; /*OH */
    cvms[6] *= 4.615239012974499e+06; /*H2O */
    cvms[7] *= 3.193192012977835e+06; /*C2H2 */
    cvms[8] *= 2.968349425484326e+06; /*CO */
    cvms[9] *= 2.963733033722604e+06; /*C2H4 */
    cvms[10] *= 2.860941121011349e+06; /*C2H5 */
    cvms[11] *= 2.769058254894261e+06; /*CH2O */
    cvms[12] *= 2.765040511976673e+06; /*C2H6 */
    cvms[13] *= 2.679121853578248e+06; /*CH3O */
    cvms[14] *= 2.598381814318037e+06; /*O2 */
    cvms[15] *= 2.519031701678171e+06; /*HO2 */
    cvms[16] *= 2.444384405113783e+06; /*H2O2 */
    cvms[17] *= 1.889234139098090e+06; /*CO2 */
    cvms[18] *= 1.887363070152301e+06; /*CH3HCO */
    cvms[19] *= 1.806485436783514e+06; /*HCOOH */
    cvms[20] *= 1.804774610197805e+06; /*CH3OCH3 */
    cvms[21] *= 1.408164720439543e+06; /*CH3OCO */
    cvms[22] *= 1.384529127447131e+06; /*CH3OCHO */
    cvms[23] *= 1.339560926789124e+06; /*CH3OCH2OH */
    cvms[24] *= 1.107945282000351e+06; /*OCH2OCHO */
    cvms[25] *= 1.107945282000351e+06; /*HOCH2OCO */
    cvms[26] *= 1.078960840432207e+06; /*CH3OCH2O2 */
    cvms[27] *= 9.032427183917572e+05; /*HO2CH2OCHO */
    cvms[28] *= 7.623853661063744e+05; /*O2CH2OCH2O2H */
    cvms[29] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
void CKCPMS(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 8.248767324424338e+07; /*H */
    cpms[1] *= 4.124383662212169e+07; /*H2 */
    cpms[2] *= 5.530081023953346e+06; /*CH3 */
    cpms[3] *= 5.196763628636074e+06; /*O */
    cpms[4] *= 5.182630712527496e+06; /*CH4 */
    cpms[5] *= 4.888768810227566e+06; /*OH */
    cpms[6] *= 4.615239012974499e+06; /*H2O */
    cpms[7] *= 3.193192012977835e+06; /*C2H2 */
    cpms[8] *= 2.968349425484326e+06; /*CO */
    cpms[9] *= 2.963733033722604e+06; /*C2H4 */
    cpms[10] *= 2.860941121011349e+06; /*C2H5 */
    cpms[11] *= 2.769058254894261e+06; /*CH2O */
    cpms[12] *= 2.765040511976673e+06; /*C2H6 */
    cpms[13] *= 2.679121853578248e+06; /*CH3O */
    cpms[14] *= 2.598381814318037e+06; /*O2 */
    cpms[15] *= 2.519031701678171e+06; /*HO2 */
    cpms[16] *= 2.444384405113783e+06; /*H2O2 */
    cpms[17] *= 1.889234139098090e+06; /*CO2 */
    cpms[18] *= 1.887363070152301e+06; /*CH3HCO */
    cpms[19] *= 1.806485436783514e+06; /*HCOOH */
    cpms[20] *= 1.804774610197805e+06; /*CH3OCH3 */
    cpms[21] *= 1.408164720439543e+06; /*CH3OCO */
    cpms[22] *= 1.384529127447131e+06; /*CH3OCHO */
    cpms[23] *= 1.339560926789124e+06; /*CH3OCH2OH */
    cpms[24] *= 1.107945282000351e+06; /*OCH2OCHO */
    cpms[25] *= 1.107945282000351e+06; /*HOCH2OCO */
    cpms[26] *= 1.078960840432207e+06; /*CH3OCH2O2 */
    cpms[27] *= 9.032427183917572e+05; /*HO2CH2OCHO */
    cpms[28] *= 7.623853661063744e+05; /*O2CH2OCH2O2H */
    cpms[29] *= 2.968047434442088e+06; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    for (int i = 0; i < 30; i++)
    {
        ums[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void CKHMS(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hms, tc);
    for (int i = 0; i < 30; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int * restrict np, double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  hms)
{
    double tc[5], h[30];

    for (int i=0; i<(*np); i++) {
        tc[0] = 0.0;
        tc[1] = T[i];
        tc[2] = T[i]*T[i];
        tc[3] = T[i]*T[i]*T[i];
        tc[4] = T[i]*T[i]*T[i]*T[i];

        speciesEnthalpy(h, tc);

        hms[0*(*np)+i] = h[0];
        hms[1*(*np)+i] = h[1];
        hms[2*(*np)+i] = h[2];
        hms[3*(*np)+i] = h[3];
        hms[4*(*np)+i] = h[4];
        hms[5*(*np)+i] = h[5];
        hms[6*(*np)+i] = h[6];
        hms[7*(*np)+i] = h[7];
        hms[8*(*np)+i] = h[8];
        hms[9*(*np)+i] = h[9];
        hms[10*(*np)+i] = h[10];
        hms[11*(*np)+i] = h[11];
        hms[12*(*np)+i] = h[12];
        hms[13*(*np)+i] = h[13];
        hms[14*(*np)+i] = h[14];
        hms[15*(*np)+i] = h[15];
        hms[16*(*np)+i] = h[16];
        hms[17*(*np)+i] = h[17];
        hms[18*(*np)+i] = h[18];
        hms[19*(*np)+i] = h[19];
        hms[20*(*np)+i] = h[20];
        hms[21*(*np)+i] = h[21];
        hms[22*(*np)+i] = h[22];
        hms[23*(*np)+i] = h[23];
        hms[24*(*np)+i] = h[24];
        hms[25*(*np)+i] = h[25];
        hms[26*(*np)+i] = h[26];
        hms[27*(*np)+i] = h[27];
        hms[28*(*np)+i] = h[28];
        hms[29*(*np)+i] = h[29];
    }

    for (int n=0; n<30; n++) {
        for (int i=0; i<(*np); i++) {
            hms[n*(*np)+i] *= 8.31451e+07 * T[i] * imw[n];
        }
    }
}


/*Returns gibbs in mass units (Eq 31.) */
void CKGMS(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  gms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gms, tc);
    for (int i = 0; i < 30; i++)
    {
        gms[i] *= RT*imw[i];
    }
}


/*Returns helmholtz in mass units (Eq 32.) */
void CKAMS(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  ams)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(ams, tc);
    for (int i = 0; i < 30; i++)
    {
        ams[i] *= RT*imw[i];
    }
}


/*Returns the entropies in mass units (Eq 28.) */
void CKSMS(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  sms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 8.248767324424338e+07; /*H */
    sms[1] *= 4.124383662212169e+07; /*H2 */
    sms[2] *= 5.530081023953346e+06; /*CH3 */
    sms[3] *= 5.196763628636074e+06; /*O */
    sms[4] *= 5.182630712527496e+06; /*CH4 */
    sms[5] *= 4.888768810227566e+06; /*OH */
    sms[6] *= 4.615239012974499e+06; /*H2O */
    sms[7] *= 3.193192012977835e+06; /*C2H2 */
    sms[8] *= 2.968349425484326e+06; /*CO */
    sms[9] *= 2.963733033722604e+06; /*C2H4 */
    sms[10] *= 2.860941121011349e+06; /*C2H5 */
    sms[11] *= 2.769058254894261e+06; /*CH2O */
    sms[12] *= 2.765040511976673e+06; /*C2H6 */
    sms[13] *= 2.679121853578248e+06; /*CH3O */
    sms[14] *= 2.598381814318037e+06; /*O2 */
    sms[15] *= 2.519031701678171e+06; /*HO2 */
    sms[16] *= 2.444384405113783e+06; /*H2O2 */
    sms[17] *= 1.889234139098090e+06; /*CO2 */
    sms[18] *= 1.887363070152301e+06; /*CH3HCO */
    sms[19] *= 1.806485436783514e+06; /*HCOOH */
    sms[20] *= 1.804774610197805e+06; /*CH3OCH3 */
    sms[21] *= 1.408164720439543e+06; /*CH3OCO */
    sms[22] *= 1.384529127447131e+06; /*CH3OCHO */
    sms[23] *= 1.339560926789124e+06; /*CH3OCH2OH */
    sms[24] *= 1.107945282000351e+06; /*OCH2OCHO */
    sms[25] *= 1.107945282000351e+06; /*HOCH2OCO */
    sms[26] *= 1.078960840432207e+06; /*CH3OCH2O2 */
    sms[27] *= 9.032427183917572e+05; /*HO2CH2OCHO */
    sms[28] *= 7.623853661063744e+05; /*O2CH2OCH2O2H */
    sms[29] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double * restrict T, double * restrict x, int * iwrk, double * restrict  rwrk, double * restrict  cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[30]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 30; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
void CKCPBS(double * restrict T, double * restrict y, int * iwrk, double * restrict  rwrk, double * restrict  cpbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[30], tresult[30]; /* temporary storage */
    cp_R(cpor, tc);
    for (int i = 0; i < 30; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 30; i++)
    {
        result += tresult[i];
    }

    *cpbs = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double * restrict T, double * restrict x, int * iwrk, double * restrict  rwrk, double * restrict  cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[30]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 30; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
void CKCVBS(double * restrict T, double * restrict y, int * iwrk, double * restrict  rwrk, double * restrict  cvbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[30]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H */
    result += cvor[1]*y[1]*imw[1]; /*H2 */
    result += cvor[2]*y[2]*imw[2]; /*CH3 */
    result += cvor[3]*y[3]*imw[3]; /*O */
    result += cvor[4]*y[4]*imw[4]; /*CH4 */
    result += cvor[5]*y[5]*imw[5]; /*OH */
    result += cvor[6]*y[6]*imw[6]; /*H2O */
    result += cvor[7]*y[7]*imw[7]; /*C2H2 */
    result += cvor[8]*y[8]*imw[8]; /*CO */
    result += cvor[9]*y[9]*imw[9]; /*C2H4 */
    result += cvor[10]*y[10]*imw[10]; /*C2H5 */
    result += cvor[11]*y[11]*imw[11]; /*CH2O */
    result += cvor[12]*y[12]*imw[12]; /*C2H6 */
    result += cvor[13]*y[13]*imw[13]; /*CH3O */
    result += cvor[14]*y[14]*imw[14]; /*O2 */
    result += cvor[15]*y[15]*imw[15]; /*HO2 */
    result += cvor[16]*y[16]*imw[16]; /*H2O2 */
    result += cvor[17]*y[17]*imw[17]; /*CO2 */
    result += cvor[18]*y[18]*imw[18]; /*CH3HCO */
    result += cvor[19]*y[19]*imw[19]; /*HCOOH */
    result += cvor[20]*y[20]*imw[20]; /*CH3OCH3 */
    result += cvor[21]*y[21]*imw[21]; /*CH3OCO */
    result += cvor[22]*y[22]*imw[22]; /*CH3OCHO */
    result += cvor[23]*y[23]*imw[23]; /*CH3OCH2OH */
    result += cvor[24]*y[24]*imw[24]; /*OCH2OCHO */
    result += cvor[25]*y[25]*imw[25]; /*HOCH2OCO */
    result += cvor[26]*y[26]*imw[26]; /*CH3OCH2O2 */
    result += cvor[27]*y[27]*imw[27]; /*HO2CH2OCHO */
    result += cvor[28]*y[28]*imw[28]; /*O2CH2OCH2O2H */
    result += cvor[29]*y[29]*imw[29]; /*N2 */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double * restrict T, double * restrict x, int * iwrk, double * restrict  rwrk, double * restrict  hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[30]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 30; ++id) {
        result += x[id]*hml[id];
    }

    *hbml = result * RT;
}


/*Returns mean enthalpy of mixture in mass units */
void CKHBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict  rwrk, double * restrict  hbms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[30], tmp[30]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    int id;
    for (id = 0; id < 30; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 30; ++id) {
        result += tmp[id];
    }

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double * restrict T, double * restrict x, int * iwrk, double * restrict  rwrk, double * restrict  ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[30]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 30; ++id) {
        result += x[id]*uml[id];
    }

    *ubml = result * RT;
}


/*get mean internal energy in mass units */
void CKUBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict  rwrk, double * restrict  ubms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double ums[30]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H */
    result += y[1]*ums[1]*imw[1]; /*H2 */
    result += y[2]*ums[2]*imw[2]; /*CH3 */
    result += y[3]*ums[3]*imw[3]; /*O */
    result += y[4]*ums[4]*imw[4]; /*CH4 */
    result += y[5]*ums[5]*imw[5]; /*OH */
    result += y[6]*ums[6]*imw[6]; /*H2O */
    result += y[7]*ums[7]*imw[7]; /*C2H2 */
    result += y[8]*ums[8]*imw[8]; /*CO */
    result += y[9]*ums[9]*imw[9]; /*C2H4 */
    result += y[10]*ums[10]*imw[10]; /*C2H5 */
    result += y[11]*ums[11]*imw[11]; /*CH2O */
    result += y[12]*ums[12]*imw[12]; /*C2H6 */
    result += y[13]*ums[13]*imw[13]; /*CH3O */
    result += y[14]*ums[14]*imw[14]; /*O2 */
    result += y[15]*ums[15]*imw[15]; /*HO2 */
    result += y[16]*ums[16]*imw[16]; /*H2O2 */
    result += y[17]*ums[17]*imw[17]; /*CO2 */
    result += y[18]*ums[18]*imw[18]; /*CH3HCO */
    result += y[19]*ums[19]*imw[19]; /*HCOOH */
    result += y[20]*ums[20]*imw[20]; /*CH3OCH3 */
    result += y[21]*ums[21]*imw[21]; /*CH3OCO */
    result += y[22]*ums[22]*imw[22]; /*CH3OCHO */
    result += y[23]*ums[23]*imw[23]; /*CH3OCH2OH */
    result += y[24]*ums[24]*imw[24]; /*OCH2OCHO */
    result += y[25]*ums[25]*imw[25]; /*HOCH2OCO */
    result += y[26]*ums[26]*imw[26]; /*CH3OCH2O2 */
    result += y[27]*ums[27]*imw[27]; /*HO2CH2OCHO */
    result += y[28]*ums[28]*imw[28]; /*O2CH2OCH2O2H */
    result += y[29]*ums[29]*imw[29]; /*N2 */

    *ubms = result * RT;
}


/*get mixture entropy in molar units */
void CKSBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict  rwrk, double * restrict  sbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[30]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 30; ++id) {
        result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
    }

    *sbml = result * 8.31451e+07;
}


/*get mixture entropy in mass units */
void CKSBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict  rwrk, double * restrict  sbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[30]; /* temporary storage */
    double x[30]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*CH3 */
    YOW += y[3]*imw[3]; /*O */
    YOW += y[4]*imw[4]; /*CH4 */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*H2O */
    YOW += y[7]*imw[7]; /*C2H2 */
    YOW += y[8]*imw[8]; /*CO */
    YOW += y[9]*imw[9]; /*C2H4 */
    YOW += y[10]*imw[10]; /*C2H5 */
    YOW += y[11]*imw[11]; /*CH2O */
    YOW += y[12]*imw[12]; /*C2H6 */
    YOW += y[13]*imw[13]; /*CH3O */
    YOW += y[14]*imw[14]; /*O2 */
    YOW += y[15]*imw[15]; /*HO2 */
    YOW += y[16]*imw[16]; /*H2O2 */
    YOW += y[17]*imw[17]; /*CO2 */
    YOW += y[18]*imw[18]; /*CH3HCO */
    YOW += y[19]*imw[19]; /*HCOOH */
    YOW += y[20]*imw[20]; /*CH3OCH3 */
    YOW += y[21]*imw[21]; /*CH3OCO */
    YOW += y[22]*imw[22]; /*CH3OCHO */
    YOW += y[23]*imw[23]; /*CH3OCH2OH */
    YOW += y[24]*imw[24]; /*OCH2OCHO */
    YOW += y[25]*imw[25]; /*HOCH2OCO */
    YOW += y[26]*imw[26]; /*CH3OCH2O2 */
    YOW += y[27]*imw[27]; /*HO2CH2OCHO */
    YOW += y[28]*imw[28]; /*O2CH2OCH2O2H */
    YOW += y[29]*imw[29]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(2.015940*YOW); 
    x[2] = y[2]/(15.035060*YOW); 
    x[3] = y[3]/(15.999400*YOW); 
    x[4] = y[4]/(16.043030*YOW); 
    x[5] = y[5]/(17.007370*YOW); 
    x[6] = y[6]/(18.015340*YOW); 
    x[7] = y[7]/(26.038240*YOW); 
    x[8] = y[8]/(28.010550*YOW); 
    x[9] = y[9]/(28.054180*YOW); 
    x[10] = y[10]/(29.062150*YOW); 
    x[11] = y[11]/(30.026490*YOW); 
    x[12] = y[12]/(30.070120*YOW); 
    x[13] = y[13]/(31.034460*YOW); 
    x[14] = y[14]/(31.998800*YOW); 
    x[15] = y[15]/(33.006770*YOW); 
    x[16] = y[16]/(34.014740*YOW); 
    x[17] = y[17]/(44.009950*YOW); 
    x[18] = y[18]/(44.053580*YOW); 
    x[19] = y[19]/(46.025890*YOW); 
    x[20] = y[20]/(46.069520*YOW); 
    x[21] = y[21]/(59.045010*YOW); 
    x[22] = y[22]/(60.052980*YOW); 
    x[23] = y[23]/(62.068920*YOW); 
    x[24] = y[24]/(75.044410*YOW); 
    x[25] = y[25]/(75.044410*YOW); 
    x[26] = y[26]/(77.060350*YOW); 
    x[27] = y[27]/(92.051780*YOW); 
    x[28] = y[28]/(109.059150*YOW); 
    x[29] = y[29]/(28.013400*YOW); 
    speciesEntropy(sor, tc);
    /*Perform computation in Eq 42 and 43 */
    result += x[0]*(sor[0]-log((x[0]+1e-100))-logPratio);
    result += x[1]*(sor[1]-log((x[1]+1e-100))-logPratio);
    result += x[2]*(sor[2]-log((x[2]+1e-100))-logPratio);
    result += x[3]*(sor[3]-log((x[3]+1e-100))-logPratio);
    result += x[4]*(sor[4]-log((x[4]+1e-100))-logPratio);
    result += x[5]*(sor[5]-log((x[5]+1e-100))-logPratio);
    result += x[6]*(sor[6]-log((x[6]+1e-100))-logPratio);
    result += x[7]*(sor[7]-log((x[7]+1e-100))-logPratio);
    result += x[8]*(sor[8]-log((x[8]+1e-100))-logPratio);
    result += x[9]*(sor[9]-log((x[9]+1e-100))-logPratio);
    result += x[10]*(sor[10]-log((x[10]+1e-100))-logPratio);
    result += x[11]*(sor[11]-log((x[11]+1e-100))-logPratio);
    result += x[12]*(sor[12]-log((x[12]+1e-100))-logPratio);
    result += x[13]*(sor[13]-log((x[13]+1e-100))-logPratio);
    result += x[14]*(sor[14]-log((x[14]+1e-100))-logPratio);
    result += x[15]*(sor[15]-log((x[15]+1e-100))-logPratio);
    result += x[16]*(sor[16]-log((x[16]+1e-100))-logPratio);
    result += x[17]*(sor[17]-log((x[17]+1e-100))-logPratio);
    result += x[18]*(sor[18]-log((x[18]+1e-100))-logPratio);
    result += x[19]*(sor[19]-log((x[19]+1e-100))-logPratio);
    result += x[20]*(sor[20]-log((x[20]+1e-100))-logPratio);
    result += x[21]*(sor[21]-log((x[21]+1e-100))-logPratio);
    result += x[22]*(sor[22]-log((x[22]+1e-100))-logPratio);
    result += x[23]*(sor[23]-log((x[23]+1e-100))-logPratio);
    result += x[24]*(sor[24]-log((x[24]+1e-100))-logPratio);
    result += x[25]*(sor[25]-log((x[25]+1e-100))-logPratio);
    result += x[26]*(sor[26]-log((x[26]+1e-100))-logPratio);
    result += x[27]*(sor[27]-log((x[27]+1e-100))-logPratio);
    result += x[28]*(sor[28]-log((x[28]+1e-100))-logPratio);
    result += x[29]*(sor[29]-log((x[29]+1e-100))-logPratio);
    /*Scale by R/W */
    *sbms = result * 8.31451e+07 * YOW;
}


/*Returns mean gibbs free energy in molar units */
void CKGBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict  rwrk, double * restrict  gbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double gort[30]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 30; ++id) {
        result += x[id]*(gort[id]+log((x[id]+1e-100))+logPratio);
    }

    *gbml = result * RT;
}


/*Returns mixture gibbs free energy in mass units */
void CKGBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict  rwrk, double * restrict  gbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double gort[30]; /* temporary storage */
    double x[30]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*CH3 */
    YOW += y[3]*imw[3]; /*O */
    YOW += y[4]*imw[4]; /*CH4 */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*H2O */
    YOW += y[7]*imw[7]; /*C2H2 */
    YOW += y[8]*imw[8]; /*CO */
    YOW += y[9]*imw[9]; /*C2H4 */
    YOW += y[10]*imw[10]; /*C2H5 */
    YOW += y[11]*imw[11]; /*CH2O */
    YOW += y[12]*imw[12]; /*C2H6 */
    YOW += y[13]*imw[13]; /*CH3O */
    YOW += y[14]*imw[14]; /*O2 */
    YOW += y[15]*imw[15]; /*HO2 */
    YOW += y[16]*imw[16]; /*H2O2 */
    YOW += y[17]*imw[17]; /*CO2 */
    YOW += y[18]*imw[18]; /*CH3HCO */
    YOW += y[19]*imw[19]; /*HCOOH */
    YOW += y[20]*imw[20]; /*CH3OCH3 */
    YOW += y[21]*imw[21]; /*CH3OCO */
    YOW += y[22]*imw[22]; /*CH3OCHO */
    YOW += y[23]*imw[23]; /*CH3OCH2OH */
    YOW += y[24]*imw[24]; /*OCH2OCHO */
    YOW += y[25]*imw[25]; /*HOCH2OCO */
    YOW += y[26]*imw[26]; /*CH3OCH2O2 */
    YOW += y[27]*imw[27]; /*HO2CH2OCHO */
    YOW += y[28]*imw[28]; /*O2CH2OCH2O2H */
    YOW += y[29]*imw[29]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(2.015940*YOW); 
    x[2] = y[2]/(15.035060*YOW); 
    x[3] = y[3]/(15.999400*YOW); 
    x[4] = y[4]/(16.043030*YOW); 
    x[5] = y[5]/(17.007370*YOW); 
    x[6] = y[6]/(18.015340*YOW); 
    x[7] = y[7]/(26.038240*YOW); 
    x[8] = y[8]/(28.010550*YOW); 
    x[9] = y[9]/(28.054180*YOW); 
    x[10] = y[10]/(29.062150*YOW); 
    x[11] = y[11]/(30.026490*YOW); 
    x[12] = y[12]/(30.070120*YOW); 
    x[13] = y[13]/(31.034460*YOW); 
    x[14] = y[14]/(31.998800*YOW); 
    x[15] = y[15]/(33.006770*YOW); 
    x[16] = y[16]/(34.014740*YOW); 
    x[17] = y[17]/(44.009950*YOW); 
    x[18] = y[18]/(44.053580*YOW); 
    x[19] = y[19]/(46.025890*YOW); 
    x[20] = y[20]/(46.069520*YOW); 
    x[21] = y[21]/(59.045010*YOW); 
    x[22] = y[22]/(60.052980*YOW); 
    x[23] = y[23]/(62.068920*YOW); 
    x[24] = y[24]/(75.044410*YOW); 
    x[25] = y[25]/(75.044410*YOW); 
    x[26] = y[26]/(77.060350*YOW); 
    x[27] = y[27]/(92.051780*YOW); 
    x[28] = y[28]/(109.059150*YOW); 
    x[29] = y[29]/(28.013400*YOW); 
    gibbs(gort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(gort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(gort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(gort[2]+log((x[2]+1e-100))+logPratio);
    result += x[3]*(gort[3]+log((x[3]+1e-100))+logPratio);
    result += x[4]*(gort[4]+log((x[4]+1e-100))+logPratio);
    result += x[5]*(gort[5]+log((x[5]+1e-100))+logPratio);
    result += x[6]*(gort[6]+log((x[6]+1e-100))+logPratio);
    result += x[7]*(gort[7]+log((x[7]+1e-100))+logPratio);
    result += x[8]*(gort[8]+log((x[8]+1e-100))+logPratio);
    result += x[9]*(gort[9]+log((x[9]+1e-100))+logPratio);
    result += x[10]*(gort[10]+log((x[10]+1e-100))+logPratio);
    result += x[11]*(gort[11]+log((x[11]+1e-100))+logPratio);
    result += x[12]*(gort[12]+log((x[12]+1e-100))+logPratio);
    result += x[13]*(gort[13]+log((x[13]+1e-100))+logPratio);
    result += x[14]*(gort[14]+log((x[14]+1e-100))+logPratio);
    result += x[15]*(gort[15]+log((x[15]+1e-100))+logPratio);
    result += x[16]*(gort[16]+log((x[16]+1e-100))+logPratio);
    result += x[17]*(gort[17]+log((x[17]+1e-100))+logPratio);
    result += x[18]*(gort[18]+log((x[18]+1e-100))+logPratio);
    result += x[19]*(gort[19]+log((x[19]+1e-100))+logPratio);
    result += x[20]*(gort[20]+log((x[20]+1e-100))+logPratio);
    result += x[21]*(gort[21]+log((x[21]+1e-100))+logPratio);
    result += x[22]*(gort[22]+log((x[22]+1e-100))+logPratio);
    result += x[23]*(gort[23]+log((x[23]+1e-100))+logPratio);
    result += x[24]*(gort[24]+log((x[24]+1e-100))+logPratio);
    result += x[25]*(gort[25]+log((x[25]+1e-100))+logPratio);
    result += x[26]*(gort[26]+log((x[26]+1e-100))+logPratio);
    result += x[27]*(gort[27]+log((x[27]+1e-100))+logPratio);
    result += x[28]*(gort[28]+log((x[28]+1e-100))+logPratio);
    result += x[29]*(gort[29]+log((x[29]+1e-100))+logPratio);
    /*Scale by RT/W */
    *gbms = result * RT * YOW;
}


/*Returns mean helmholtz free energy in molar units */
void CKABML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict  rwrk, double * restrict  abml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double aort[30]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 30; ++id) {
        result += x[id]*(aort[id]+log((x[id]+1e-100))+logPratio);
    }

    *abml = result * RT;
}


/*Returns mixture helmholtz free energy in mass units */
void CKABMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict  rwrk, double * restrict  abms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double aort[30]; /* temporary storage */
    double x[30]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*CH3 */
    YOW += y[3]*imw[3]; /*O */
    YOW += y[4]*imw[4]; /*CH4 */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*H2O */
    YOW += y[7]*imw[7]; /*C2H2 */
    YOW += y[8]*imw[8]; /*CO */
    YOW += y[9]*imw[9]; /*C2H4 */
    YOW += y[10]*imw[10]; /*C2H5 */
    YOW += y[11]*imw[11]; /*CH2O */
    YOW += y[12]*imw[12]; /*C2H6 */
    YOW += y[13]*imw[13]; /*CH3O */
    YOW += y[14]*imw[14]; /*O2 */
    YOW += y[15]*imw[15]; /*HO2 */
    YOW += y[16]*imw[16]; /*H2O2 */
    YOW += y[17]*imw[17]; /*CO2 */
    YOW += y[18]*imw[18]; /*CH3HCO */
    YOW += y[19]*imw[19]; /*HCOOH */
    YOW += y[20]*imw[20]; /*CH3OCH3 */
    YOW += y[21]*imw[21]; /*CH3OCO */
    YOW += y[22]*imw[22]; /*CH3OCHO */
    YOW += y[23]*imw[23]; /*CH3OCH2OH */
    YOW += y[24]*imw[24]; /*OCH2OCHO */
    YOW += y[25]*imw[25]; /*HOCH2OCO */
    YOW += y[26]*imw[26]; /*CH3OCH2O2 */
    YOW += y[27]*imw[27]; /*HO2CH2OCHO */
    YOW += y[28]*imw[28]; /*O2CH2OCH2O2H */
    YOW += y[29]*imw[29]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(2.015940*YOW); 
    x[2] = y[2]/(15.035060*YOW); 
    x[3] = y[3]/(15.999400*YOW); 
    x[4] = y[4]/(16.043030*YOW); 
    x[5] = y[5]/(17.007370*YOW); 
    x[6] = y[6]/(18.015340*YOW); 
    x[7] = y[7]/(26.038240*YOW); 
    x[8] = y[8]/(28.010550*YOW); 
    x[9] = y[9]/(28.054180*YOW); 
    x[10] = y[10]/(29.062150*YOW); 
    x[11] = y[11]/(30.026490*YOW); 
    x[12] = y[12]/(30.070120*YOW); 
    x[13] = y[13]/(31.034460*YOW); 
    x[14] = y[14]/(31.998800*YOW); 
    x[15] = y[15]/(33.006770*YOW); 
    x[16] = y[16]/(34.014740*YOW); 
    x[17] = y[17]/(44.009950*YOW); 
    x[18] = y[18]/(44.053580*YOW); 
    x[19] = y[19]/(46.025890*YOW); 
    x[20] = y[20]/(46.069520*YOW); 
    x[21] = y[21]/(59.045010*YOW); 
    x[22] = y[22]/(60.052980*YOW); 
    x[23] = y[23]/(62.068920*YOW); 
    x[24] = y[24]/(75.044410*YOW); 
    x[25] = y[25]/(75.044410*YOW); 
    x[26] = y[26]/(77.060350*YOW); 
    x[27] = y[27]/(92.051780*YOW); 
    x[28] = y[28]/(109.059150*YOW); 
    x[29] = y[29]/(28.013400*YOW); 
    helmholtz(aort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(aort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(aort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(aort[2]+log((x[2]+1e-100))+logPratio);
    result += x[3]*(aort[3]+log((x[3]+1e-100))+logPratio);
    result += x[4]*(aort[4]+log((x[4]+1e-100))+logPratio);
    result += x[5]*(aort[5]+log((x[5]+1e-100))+logPratio);
    result += x[6]*(aort[6]+log((x[6]+1e-100))+logPratio);
    result += x[7]*(aort[7]+log((x[7]+1e-100))+logPratio);
    result += x[8]*(aort[8]+log((x[8]+1e-100))+logPratio);
    result += x[9]*(aort[9]+log((x[9]+1e-100))+logPratio);
    result += x[10]*(aort[10]+log((x[10]+1e-100))+logPratio);
    result += x[11]*(aort[11]+log((x[11]+1e-100))+logPratio);
    result += x[12]*(aort[12]+log((x[12]+1e-100))+logPratio);
    result += x[13]*(aort[13]+log((x[13]+1e-100))+logPratio);
    result += x[14]*(aort[14]+log((x[14]+1e-100))+logPratio);
    result += x[15]*(aort[15]+log((x[15]+1e-100))+logPratio);
    result += x[16]*(aort[16]+log((x[16]+1e-100))+logPratio);
    result += x[17]*(aort[17]+log((x[17]+1e-100))+logPratio);
    result += x[18]*(aort[18]+log((x[18]+1e-100))+logPratio);
    result += x[19]*(aort[19]+log((x[19]+1e-100))+logPratio);
    result += x[20]*(aort[20]+log((x[20]+1e-100))+logPratio);
    result += x[21]*(aort[21]+log((x[21]+1e-100))+logPratio);
    result += x[22]*(aort[22]+log((x[22]+1e-100))+logPratio);
    result += x[23]*(aort[23]+log((x[23]+1e-100))+logPratio);
    result += x[24]*(aort[24]+log((x[24]+1e-100))+logPratio);
    result += x[25]*(aort[25]+log((x[25]+1e-100))+logPratio);
    result += x[26]*(aort[26]+log((x[26]+1e-100))+logPratio);
    result += x[27]*(aort[27]+log((x[27]+1e-100))+logPratio);
    result += x[28]*(aort[28]+log((x[28]+1e-100))+logPratio);
    result += x[29]*(aort[29]+log((x[29]+1e-100))+logPratio);
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double * restrict  T, double * restrict  C, int * iwrk, double * restrict  rwrk, double * restrict  wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 30; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 30; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  wdot)
{
    FORT_CKWYP(P, T, y, iwrk, rwrk, wdot);
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  wdot)
{
    double y[30];
    CKXTY(x, iwrk, rwrk, y);
    FORT_CKWYP(P, T, y, iwrk, rwrk, wdot);
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  wdot)
{
    double pres;
    CKPY(rho, T, y, iwrk, rwrk, &pres);
    FORT_CKWYP(&pres, T, y, iwrk, rwrk, wdot);
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int * restrict np, double * restrict rho, double * restrict T,
	    double * restrict y, int * restrict iwrk, double * restrict rwrk,
	    double * restrict wdot)
{
    printf("VCKWYR not supported!\n");
    exit(1);
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  wdot)
{
    double y[30], pres;
    CKXTY(x, iwrk, rwrk, y);
    CKPY(rho, T, y, iwrk, rwrk, &pres);
    FORT_CKWYP(&pres, T, y, iwrk, rwrk, wdot);
}


/*Returns the rate of progress for each reaction */
void CKQC(double * restrict  T, double * restrict  C, int * iwrk, double * restrict  rwrk, double * restrict  qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 30; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 30; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 3; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  q_f, double * restrict  q_r)
{
    int id; /*loop counter */
    double c[30]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 30; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  qdot)
{
    int id; /*loop counter */
    double c[30]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*CH3 */
    YOW += y[3]*imw[3]; /*O */
    YOW += y[4]*imw[4]; /*CH4 */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*H2O */
    YOW += y[7]*imw[7]; /*C2H2 */
    YOW += y[8]*imw[8]; /*CO */
    YOW += y[9]*imw[9]; /*C2H4 */
    YOW += y[10]*imw[10]; /*C2H5 */
    YOW += y[11]*imw[11]; /*CH2O */
    YOW += y[12]*imw[12]; /*C2H6 */
    YOW += y[13]*imw[13]; /*CH3O */
    YOW += y[14]*imw[14]; /*O2 */
    YOW += y[15]*imw[15]; /*HO2 */
    YOW += y[16]*imw[16]; /*H2O2 */
    YOW += y[17]*imw[17]; /*CO2 */
    YOW += y[18]*imw[18]; /*CH3HCO */
    YOW += y[19]*imw[19]; /*HCOOH */
    YOW += y[20]*imw[20]; /*CH3OCH3 */
    YOW += y[21]*imw[21]; /*CH3OCO */
    YOW += y[22]*imw[22]; /*CH3OCHO */
    YOW += y[23]*imw[23]; /*CH3OCH2OH */
    YOW += y[24]*imw[24]; /*OCH2OCHO */
    YOW += y[25]*imw[25]; /*HOCH2OCO */
    YOW += y[26]*imw[26]; /*CH3OCH2O2 */
    YOW += y[27]*imw[27]; /*HO2CH2OCHO */
    YOW += y[28]*imw[28]; /*O2CH2OCH2O2H */
    YOW += y[29]*imw[29]; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]*imw[0]; 
    c[1] = PWORT * y[1]*imw[1]; 
    c[2] = PWORT * y[2]*imw[2]; 
    c[3] = PWORT * y[3]*imw[3]; 
    c[4] = PWORT * y[4]*imw[4]; 
    c[5] = PWORT * y[5]*imw[5]; 
    c[6] = PWORT * y[6]*imw[6]; 
    c[7] = PWORT * y[7]*imw[7]; 
    c[8] = PWORT * y[8]*imw[8]; 
    c[9] = PWORT * y[9]*imw[9]; 
    c[10] = PWORT * y[10]*imw[10]; 
    c[11] = PWORT * y[11]*imw[11]; 
    c[12] = PWORT * y[12]*imw[12]; 
    c[13] = PWORT * y[13]*imw[13]; 
    c[14] = PWORT * y[14]*imw[14]; 
    c[15] = PWORT * y[15]*imw[15]; 
    c[16] = PWORT * y[16]*imw[16]; 
    c[17] = PWORT * y[17]*imw[17]; 
    c[18] = PWORT * y[18]*imw[18]; 
    c[19] = PWORT * y[19]*imw[19]; 
    c[20] = PWORT * y[20]*imw[20]; 
    c[21] = PWORT * y[21]*imw[21]; 
    c[22] = PWORT * y[22]*imw[22]; 
    c[23] = PWORT * y[23]*imw[23]; 
    c[24] = PWORT * y[24]*imw[24]; 
    c[25] = PWORT * y[25]*imw[25]; 
    c[26] = PWORT * y[26]*imw[26]; 
    c[27] = PWORT * y[27]*imw[27]; 
    c[28] = PWORT * y[28]*imw[28]; 
    c[29] = PWORT * y[29]*imw[29]; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  qdot)
{
    int id; /*loop counter */
    double c[30]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 30; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  qdot)
{
    int id; /*loop counter */
    double c[30]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 
    c[3] = 1e6 * (*rho) * y[3]*imw[3]; 
    c[4] = 1e6 * (*rho) * y[4]*imw[4]; 
    c[5] = 1e6 * (*rho) * y[5]*imw[5]; 
    c[6] = 1e6 * (*rho) * y[6]*imw[6]; 
    c[7] = 1e6 * (*rho) * y[7]*imw[7]; 
    c[8] = 1e6 * (*rho) * y[8]*imw[8]; 
    c[9] = 1e6 * (*rho) * y[9]*imw[9]; 
    c[10] = 1e6 * (*rho) * y[10]*imw[10]; 
    c[11] = 1e6 * (*rho) * y[11]*imw[11]; 
    c[12] = 1e6 * (*rho) * y[12]*imw[12]; 
    c[13] = 1e6 * (*rho) * y[13]*imw[13]; 
    c[14] = 1e6 * (*rho) * y[14]*imw[14]; 
    c[15] = 1e6 * (*rho) * y[15]*imw[15]; 
    c[16] = 1e6 * (*rho) * y[16]*imw[16]; 
    c[17] = 1e6 * (*rho) * y[17]*imw[17]; 
    c[18] = 1e6 * (*rho) * y[18]*imw[18]; 
    c[19] = 1e6 * (*rho) * y[19]*imw[19]; 
    c[20] = 1e6 * (*rho) * y[20]*imw[20]; 
    c[21] = 1e6 * (*rho) * y[21]*imw[21]; 
    c[22] = 1e6 * (*rho) * y[22]*imw[22]; 
    c[23] = 1e6 * (*rho) * y[23]*imw[23]; 
    c[24] = 1e6 * (*rho) * y[24]*imw[24]; 
    c[25] = 1e6 * (*rho) * y[25]*imw[25]; 
    c[26] = 1e6 * (*rho) * y[26]*imw[26]; 
    c[27] = 1e6 * (*rho) * y[27]*imw[27]; 
    c[28] = 1e6 * (*rho) * y[28]*imw[28]; 
    c[29] = 1e6 * (*rho) * y[29]*imw[29]; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  qdot)
{
    int id; /*loop counter */
    double c[30]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.035060; /*CH3 */
    XW += x[3]*15.999400; /*O */
    XW += x[4]*16.043030; /*CH4 */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*26.038240; /*C2H2 */
    XW += x[8]*28.010550; /*CO */
    XW += x[9]*28.054180; /*C2H4 */
    XW += x[10]*29.062150; /*C2H5 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*30.070120; /*C2H6 */
    XW += x[13]*31.034460; /*CH3O */
    XW += x[14]*31.998800; /*O2 */
    XW += x[15]*33.006770; /*HO2 */
    XW += x[16]*34.014740; /*H2O2 */
    XW += x[17]*44.009950; /*CO2 */
    XW += x[18]*44.053580; /*CH3HCO */
    XW += x[19]*46.025890; /*HCOOH */
    XW += x[20]*46.069520; /*CH3OCH3 */
    XW += x[21]*59.045010; /*CH3OCO */
    XW += x[22]*60.052980; /*CH3OCHO */
    XW += x[23]*62.068920; /*CH3OCH2OH */
    XW += x[24]*75.044410; /*OCH2OCHO */
    XW += x[25]*75.044410; /*HOCH2OCO */
    XW += x[26]*77.060350; /*CH3OCH2O2 */
    XW += x[27]*92.051780; /*HO2CH2OCHO */
    XW += x[28]*109.059150; /*O2CH2OCH2O2H */
    XW += x[29]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 30; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the stoichiometric coefficients */
/*of the reaction mechanism. (Eq 50) */
void CKNU(int * kdim, int * iwrk, double * restrict  rwrk, int * nuki)
{
    int id; /*loop counter */
    int kd = (*kdim); 
    /*Zero nuki */
    for (id = 0; id < 30 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: H + O2 <=> O + OH */
    nuki[ 0 * kd + 0 ] += -1 ;
    nuki[ 14 * kd + 0 ] += -1 ;
    nuki[ 3 * kd + 0 ] += +1 ;
    nuki[ 5 * kd + 0 ] += +1 ;

    /*reaction 2: O + H2 <=> H + OH */
    nuki[ 3 * kd + 1 ] += -1 ;
    nuki[ 1 * kd + 1 ] += -1 ;
    nuki[ 0 * kd + 1 ] += +1 ;
    nuki[ 5 * kd + 1 ] += +1 ;

    /*reaction 3: H2 + OH <=> H2O + H */
    nuki[ 1 * kd + 2 ] += -1 ;
    nuki[ 5 * kd + 2 ] += -1 ;
    nuki[ 6 * kd + 2 ] += +1 ;
    nuki[ 0 * kd + 2 ] += +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * restrict  rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < 4 * 30; ++ id) {
         ncf[id] = 0; 
    }

    /*H */
    ncf[ 0 * kd + 1 ] = 1; /*H */

    /*H2 */
    ncf[ 1 * kd + 1 ] = 2; /*H */

    /*CH3 */
    ncf[ 2 * kd + 0 ] = 1; /*C */
    ncf[ 2 * kd + 1 ] = 3; /*H */

    /*O */
    ncf[ 3 * kd + 2 ] = 1; /*O */

    /*CH4 */
    ncf[ 4 * kd + 0 ] = 1; /*C */
    ncf[ 4 * kd + 1 ] = 4; /*H */

    /*OH */
    ncf[ 5 * kd + 2 ] = 1; /*O */
    ncf[ 5 * kd + 1 ] = 1; /*H */

    /*H2O */
    ncf[ 6 * kd + 1 ] = 2; /*H */
    ncf[ 6 * kd + 2 ] = 1; /*O */

    /*C2H2 */
    ncf[ 7 * kd + 0 ] = 2; /*C */
    ncf[ 7 * kd + 1 ] = 2; /*H */

    /*CO */
    ncf[ 8 * kd + 0 ] = 1; /*C */
    ncf[ 8 * kd + 2 ] = 1; /*O */

    /*C2H4 */
    ncf[ 9 * kd + 0 ] = 2; /*C */
    ncf[ 9 * kd + 1 ] = 4; /*H */

    /*C2H5 */
    ncf[ 10 * kd + 0 ] = 2; /*C */
    ncf[ 10 * kd + 1 ] = 5; /*H */

    /*CH2O */
    ncf[ 11 * kd + 0 ] = 1; /*C */
    ncf[ 11 * kd + 1 ] = 2; /*H */
    ncf[ 11 * kd + 2 ] = 1; /*O */

    /*C2H6 */
    ncf[ 12 * kd + 0 ] = 2; /*C */
    ncf[ 12 * kd + 1 ] = 6; /*H */

    /*CH3O */
    ncf[ 13 * kd + 0 ] = 1; /*C */
    ncf[ 13 * kd + 1 ] = 3; /*H */
    ncf[ 13 * kd + 2 ] = 1; /*O */

    /*O2 */
    ncf[ 14 * kd + 2 ] = 2; /*O */

    /*HO2 */
    ncf[ 15 * kd + 1 ] = 1; /*H */
    ncf[ 15 * kd + 2 ] = 2; /*O */

    /*H2O2 */
    ncf[ 16 * kd + 1 ] = 2; /*H */
    ncf[ 16 * kd + 2 ] = 2; /*O */

    /*CO2 */
    ncf[ 17 * kd + 0 ] = 1; /*C */
    ncf[ 17 * kd + 2 ] = 2; /*O */

    /*CH3HCO */
    ncf[ 18 * kd + 0 ] = 2; /*C */
    ncf[ 18 * kd + 1 ] = 4; /*H */
    ncf[ 18 * kd + 2 ] = 1; /*O */

    /*HCOOH */
    ncf[ 19 * kd + 0 ] = 1; /*C */
    ncf[ 19 * kd + 1 ] = 2; /*H */
    ncf[ 19 * kd + 2 ] = 2; /*O */

    /*CH3OCH3 */
    ncf[ 20 * kd + 0 ] = 2; /*C */
    ncf[ 20 * kd + 1 ] = 6; /*H */
    ncf[ 20 * kd + 2 ] = 1; /*O */

    /*CH3OCO */
    ncf[ 21 * kd + 0 ] = 2; /*C */
    ncf[ 21 * kd + 1 ] = 3; /*H */
    ncf[ 21 * kd + 2 ] = 2; /*O */

    /*CH3OCHO */
    ncf[ 22 * kd + 0 ] = 2; /*C */
    ncf[ 22 * kd + 1 ] = 4; /*H */
    ncf[ 22 * kd + 2 ] = 2; /*O */

    /*CH3OCH2OH */
    ncf[ 23 * kd + 0 ] = 2; /*C */
    ncf[ 23 * kd + 1 ] = 6; /*H */
    ncf[ 23 * kd + 2 ] = 2; /*O */

    /*OCH2OCHO */
    ncf[ 24 * kd + 0 ] = 2; /*C */
    ncf[ 24 * kd + 1 ] = 3; /*H */
    ncf[ 24 * kd + 2 ] = 3; /*O */

    /*HOCH2OCO */
    ncf[ 25 * kd + 0 ] = 2; /*C */
    ncf[ 25 * kd + 1 ] = 3; /*H */
    ncf[ 25 * kd + 2 ] = 3; /*O */

    /*CH3OCH2O2 */
    ncf[ 26 * kd + 0 ] = 2; /*C */
    ncf[ 26 * kd + 1 ] = 5; /*H */
    ncf[ 26 * kd + 2 ] = 3; /*O */

    /*HO2CH2OCHO */
    ncf[ 27 * kd + 0 ] = 2; /*C */
    ncf[ 27 * kd + 1 ] = 4; /*H */
    ncf[ 27 * kd + 2 ] = 4; /*O */

    /*O2CH2OCH2O2H */
    ncf[ 28 * kd + 0 ] = 2; /*C */
    ncf[ 28 * kd + 1 ] = 5; /*H */
    ncf[ 28 * kd + 2 ] = 5; /*O */

    /*N2 */
    ncf[ 29 * kd + 3 ] = 2; /*N */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE(int * iwrk, double * restrict  rwrk, double * restrict  a, double * restrict  b, double * restrict  e)
{

    /*reaction 1: H + O2 <=> O + OH */
    a[0] = 0;
    b[0] = 0;
    e[0] = 0;

    /*reaction 2: O + H2 <=> H + OH */
    a[1] = 0;
    b[1] = 0;
    e[1] = 0;

    /*reaction 3: H2 + OH <=> H2O + H */
    a[2] = 0;
    b[2] = 0;
    e[2] = 0;

    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double * restrict  T, double * restrict  C, int * iwrk, double * restrict  rwrk, double * restrict  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[30]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> O + OH */
    /*eqcon[0] *= 1;  */

    /*reaction 2: O + H2 <=> H + OH */
    /*eqcon[1] *= 1;  */

    /*reaction 3: H2 + OH <=> H2O + H */
    /*eqcon[2] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[30]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> O + OH */
    /*eqcon[0] *= 1;  */

    /*reaction 2: O + H2 <=> H + OH */
    /*eqcon[1] *= 1;  */

    /*reaction 3: H2 + OH <=> H2O + H */
    /*eqcon[2] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[30]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> O + OH */
    /*eqcon[0] *= 1;  */

    /*reaction 2: O + H2 <=> H + OH */
    /*eqcon[1] *= 1;  */

    /*reaction 3: H2 + OH <=> H2O + H */
    /*eqcon[2] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[30]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> O + OH */
    /*eqcon[0] *= 1;  */

    /*reaction 2: O + H2 <=> H + OH */
    /*eqcon[1] *= 1;  */

    /*reaction 3: H2 + OH <=> H2O + H */
    /*eqcon[2] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[30]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> O + OH */
    /*eqcon[0] *= 1;  */

    /*reaction 2: O + H2 <=> H + OH */
    /*eqcon[1] *= 1;  */

    /*reaction 3: H2 + OH <=> H2O + H */
    /*eqcon[2] *= 1;  */
}

static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[3];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[3];
#ifdef _OPENMP
#pragma omp threadprivate(Kc_save)
#endif

/*compute the production rate for each species */
void productionRate(double * restrict  wdot, double * restrict  sc, double T)
{
    double qdot;

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[30];                /*Gibbs free energy */
    double Kc;                      /*equilibrium constant */
    double k_f;                     /*forward reaction rate */
    double k_r;                     /*reverse reaction rate */
    double q_f;                     /*forward progress rate */
    double q_r;                     /*reverse progress rate */
    double phi_f;                   /*forward phase space factor */
    double phi_r;                   /*reverse phase space factor */
    double alpha;                   /*enhancement */
    double redP;                    /*reduced pressure */
    double logPred;                 /*log of above */
    double F;                       /*fallof rate enhancement */

    double F_troe;                  /*TROE intermediate */
    double logFcent;                /*TROE intermediate */
    double troe;                    /*TROE intermediate */
    double troe_c;                  /*TROE intermediate */
    double troe_n;                  /*TROE intermediate */

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */

    double invT = 1.0 / tc[1];

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;
    double refCinv = 1 / refC;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 30; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*zero out wdot */
    for (id = 0; id < 30; ++id) {
        wdot[id] = 0.0;
    }

    if (T != T_save)
    {
        T_save = T;

        k_f_save[0] = 1e-06 * 0;
        k_f_save[1] = 1e-06 * 0;
        k_f_save[2] = 1e-06 * 0;

        Kc_save[0] = exp((g_RT[0] + g_RT[14]) - (g_RT[3] + g_RT[5]));
        Kc_save[1] = exp((g_RT[3] + g_RT[1]) - (g_RT[0] + g_RT[5]));
        Kc_save[2] = exp((g_RT[1] + g_RT[5]) - (g_RT[6] + g_RT[0]));
    }

    /*reaction 1: H + O2 <=> O + OH */
    phi_f = sc[0]*sc[14];
    k_f = k_f_save[0];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = Kc_save[0];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 2: O + H2 <=> H + OH */
    phi_f = sc[3]*sc[1];
    k_f = k_f_save[1];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[5];
    Kc = Kc_save[1];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 3: H2 + OH <=> H2O + H */
    phi_f = sc[1]*sc[5];
    k_f = k_f_save[2];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = Kc_save[2];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[0] += 1 * qdot;

    return;
}


/*compute the production rate for each species */
void vproductionRate(int npt, double * restrict wdot, double * restrict sc, double * restrict T)
{
    double k_f_s[3][npt], Kc_s[3][npt], mixture[npt], g_RT[30*npt];
    double tc[5*npt], invT[npt];

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        tc[0*npt+i] = log(T[i]);
        tc[1*npt+i] = T[i];
        tc[2*npt+i] = T[i]*T[i];
        tc[3*npt+i] = T[i]*T[i]*T[i];
        tc[4*npt+i] = T[i]*T[i]*T[i]*T[i];
        invT[i] = 1.0 / T[i];
    }

#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        k_f_s[0][i] = 1e-06 * 0.0;
        k_f_s[1][i] = 1e-06 * 0.0;
        k_f_s[2][i] = 1e-06 * 0.0;
    }

    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        double tg[5], g[30];
        tg[0] = tc[0*npt+i];
        tg[1] = tc[1*npt+i];
        tg[2] = tc[2*npt+i];
        tg[3] = tc[3*npt+i];
        tg[4] = tc[4*npt+i];

        gibbs(g, tg);

        g_RT[0*npt+i] = g[0];
        g_RT[1*npt+i] = g[1];
        g_RT[2*npt+i] = g[2];
        g_RT[3*npt+i] = g[3];
        g_RT[4*npt+i] = g[4];
        g_RT[5*npt+i] = g[5];
        g_RT[6*npt+i] = g[6];
        g_RT[7*npt+i] = g[7];
        g_RT[8*npt+i] = g[8];
        g_RT[9*npt+i] = g[9];
        g_RT[10*npt+i] = g[10];
        g_RT[11*npt+i] = g[11];
        g_RT[12*npt+i] = g[12];
        g_RT[13*npt+i] = g[13];
        g_RT[14*npt+i] = g[14];
        g_RT[15*npt+i] = g[15];
        g_RT[16*npt+i] = g[16];
        g_RT[17*npt+i] = g[17];
        g_RT[18*npt+i] = g[18];
        g_RT[19*npt+i] = g[19];
        g_RT[20*npt+i] = g[20];
        g_RT[21*npt+i] = g[21];
        g_RT[22*npt+i] = g[22];
        g_RT[23*npt+i] = g[23];
        g_RT[24*npt+i] = g[24];
        g_RT[25*npt+i] = g[25];
        g_RT[26*npt+i] = g[26];
        g_RT[27*npt+i] = g[27];
        g_RT[28*npt+i] = g[28];
        g_RT[29*npt+i] = g[29];
    }

#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        double refC = 101325. / 8.31451 / T[i];

        Kc_s[0][i] = exp((g_RT[0*npt+i] + g_RT[14*npt+i]) - (g_RT[3*npt+i] + g_RT[5*npt+i]));
        Kc_s[1][i] = exp((g_RT[3*npt+i] + g_RT[1*npt+i]) - (g_RT[0*npt+i] + g_RT[5*npt+i]));
        Kc_s[2][i] = exp((g_RT[1*npt+i] + g_RT[5*npt+i]) - (g_RT[6*npt+i] + g_RT[0*npt+i]));
    }

    for (int i=0; i<npt; i++) {
        mixture[i] = 0.0;
    }

    for (int n=0; n<30; n++) {
        for (int i=0; i<npt; i++) {
            mixture[i] += sc[n*npt+i];
            wdot[n*npt+i] = 0.0;
        }
    }

#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        double alpha, redP, F, logPred, logFcent;
        double troe_c, troe_n, troe, F_troe;
        double X, F_src;

        /*reaction 1: H + O2 <=> O + OH */
        phi_f = sc[0*npt+i]*sc[14*npt+i];
        k_f = k_f_s[0][i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[5*npt+i];
        Kc = Kc_s[0][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= 1 * qdot;
        wdot[14*npt+i] -= 1 * qdot;
        wdot[3*npt+i] += 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;

        /*reaction 2: O + H2 <=> H + OH */
        phi_f = sc[3*npt+i]*sc[1*npt+i];
        k_f = k_f_s[1][i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[5*npt+i];
        Kc = Kc_s[1][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= 1 * qdot;
        wdot[1*npt+i] -= 1 * qdot;
        wdot[0*npt+i] += 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;

        /*reaction 3: H2 + OH <=> H2O + H */
        phi_f = sc[1*npt+i]*sc[5*npt+i];
        k_f = k_f_s[2][i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[0*npt+i];
        Kc = Kc_s[2][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= 1 * qdot;
        wdot[5*npt+i] -= 1 * qdot;
        wdot[6*npt+i] += 1 * qdot;
        wdot[0*npt+i] += 1 * qdot;
    }
}

/*compute the reaction Jacobian */
void DWDOT(double * restrict  J, double * restrict  sc, double * Tp, int * consP)
{
    printf("DWDOT not supported!\n");
    exit(1);
}

/*compute the reaction Jacobian */
void aJacobian(double * restrict J, double * restrict sc, double T, int consP)
{
    for (int i=0; i<961; i++) {
        J[i] = 0.0;
    }

    double wdot[30];
    for (int k=0; k<30; k++) {
        wdot[k] = 0.0;
    }

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];
    double invT2 = invT * invT;

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;
    double refCinv = 1 / refC;

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int k = 0; k < 30; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[30];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[30];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[30];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    /*reaction 1: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[14];
    k_f = 1e-06 * 0;
    dlnkfdT = 0.0;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp((g_RT[0] + g_RT[14]) - (g_RT[3] + g_RT[5]));
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[14]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[3] += q; /* O */
    wdot[5] += q; /* OH */
    wdot[14] -= q; /* O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[14];
    J[0] += -1 * dqdci;           /* dwdot[H]/d[H] */
    J[3] += 1 * dqdci;            /* dwdot[O]/d[H] */
    J[5] += 1 * dqdci;            /* dwdot[OH]/d[H] */
    J[14] += -1 * dqdci;          /* dwdot[O2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[5];
    J[93] += -1 * dqdci;          /* dwdot[H]/d[O] */
    J[96] += 1 * dqdci;           /* dwdot[O]/d[O] */
    J[98] += 1 * dqdci;           /* dwdot[OH]/d[O] */
    J[107] += -1 * dqdci;         /* dwdot[O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[3];
    J[155] += -1 * dqdci;         /* dwdot[H]/d[OH] */
    J[158] += 1 * dqdci;          /* dwdot[O]/d[OH] */
    J[160] += 1 * dqdci;          /* dwdot[OH]/d[OH] */
    J[169] += -1 * dqdci;         /* dwdot[O2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[0];
    J[434] += -1 * dqdci;         /* dwdot[H]/d[O2] */
    J[437] += 1 * dqdci;          /* dwdot[O]/d[O2] */
    J[439] += 1 * dqdci;          /* dwdot[OH]/d[O2] */
    J[448] += -1 * dqdci;         /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[930] += -1 * dqdT;          /* dwdot[H]/dT */
    J[933] += 1 * dqdT;           /* dwdot[O]/dT */
    J[935] += 1 * dqdT;           /* dwdot[OH]/dT */
    J[944] += -1 * dqdT;          /* dwdot[O2]/dT */

    /*reaction 2: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1e-06 * 0;
    dlnkfdT = 0.0;
    /* reverse */
    phi_r = sc[0]*sc[5];
    Kc = exp((g_RT[1] + g_RT[3]) - (g_RT[0] + g_RT[5]));
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[0] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[3] -= q; /* O */
    wdot[5] += q; /* OH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5];
    J[0] += 1 * dqdci;            /* dwdot[H]/d[H] */
    J[1] += -1 * dqdci;           /* dwdot[H2]/d[H] */
    J[3] += -1 * dqdci;           /* dwdot[O]/d[H] */
    J[5] += 1 * dqdci;            /* dwdot[OH]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[3];
    J[31] += 1 * dqdci;           /* dwdot[H]/d[H2] */
    J[32] += -1 * dqdci;          /* dwdot[H2]/d[H2] */
    J[34] += -1 * dqdci;          /* dwdot[O]/d[H2] */
    J[36] += 1 * dqdci;           /* dwdot[OH]/d[H2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[1];
    J[93] += 1 * dqdci;           /* dwdot[H]/d[O] */
    J[94] += -1 * dqdci;          /* dwdot[H2]/d[O] */
    J[96] += -1 * dqdci;          /* dwdot[O]/d[O] */
    J[98] += 1 * dqdci;           /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0];
    J[155] += 1 * dqdci;          /* dwdot[H]/d[OH] */
    J[156] += -1 * dqdci;         /* dwdot[H2]/d[OH] */
    J[158] += -1 * dqdci;         /* dwdot[O]/d[OH] */
    J[160] += 1 * dqdci;          /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[930] += 1 * dqdT;           /* dwdot[H]/dT */
    J[931] += -1 * dqdT;          /* dwdot[H2]/dT */
    J[933] += -1 * dqdT;          /* dwdot[O]/dT */
    J[935] += 1 * dqdT;           /* dwdot[OH]/dT */

    /*reaction 3: H2 + OH <=> H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[5];
    k_f = 1e-06 * 0;
    dlnkfdT = 0.0;
    /* reverse */
    phi_r = sc[0]*sc[6];
    Kc = exp((g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[6]));
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[5]) + (h_RT[0] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* H2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[6];
    J[0] += 1 * dqdci;            /* dwdot[H]/d[H] */
    J[1] += -1 * dqdci;           /* dwdot[H2]/d[H] */
    J[5] += -1 * dqdci;           /* dwdot[OH]/d[H] */
    J[6] += 1 * dqdci;            /* dwdot[H2O]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[5];
    J[31] += 1 * dqdci;           /* dwdot[H]/d[H2] */
    J[32] += -1 * dqdci;          /* dwdot[H2]/d[H2] */
    J[36] += -1 * dqdci;          /* dwdot[OH]/d[H2] */
    J[37] += 1 * dqdci;           /* dwdot[H2O]/d[H2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[1];
    J[155] += 1 * dqdci;          /* dwdot[H]/d[OH] */
    J[156] += -1 * dqdci;         /* dwdot[H2]/d[OH] */
    J[160] += -1 * dqdci;         /* dwdot[OH]/d[OH] */
    J[161] += 1 * dqdci;          /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[0];
    J[186] += 1 * dqdci;          /* dwdot[H]/d[H2O] */
    J[187] += -1 * dqdci;         /* dwdot[H2]/d[H2O] */
    J[191] += -1 * dqdci;         /* dwdot[OH]/d[H2O] */
    J[192] += 1 * dqdci;          /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[930] += 1 * dqdT;           /* dwdot[H]/dT */
    J[931] += -1 * dqdT;          /* dwdot[H2]/dT */
    J[935] += -1 * dqdT;          /* dwdot[OH]/dT */
    J[936] += 1 * dqdT;           /* dwdot[H2O]/dT */

    double c_R[30], dcRdT[30], e_RT[30];
    double * eh_RT;
    if (consP) {
        cp_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        eh_RT = &h_RT[0];
    }
    else {
        cv_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        speciesInternalEnergy(e_RT, tc);
        eh_RT = &e_RT[0];
    }

    double cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 30; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[930+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 30; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 30; ++m) {
            dehmixdc += eh_RT[m]*J[k*31+m];
        }
        J[k*31+30] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[960] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}


/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void dcvpRdT(double * restrict  species, double * restrict  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1376 kelvin */
    if (T < 1376) {
        /*species 19: HCOOH */
        species[19] =
            +1.63363016e-02
            -2.12514842e-05 * tc[1]
            +9.96398931e-09 * tc[2]
            -1.60870441e-12 * tc[3];
    } else {
        /*species 19: HCOOH */
        species[19] =
            +5.14289368e-03
            -3.64477026e-06 * tc[1]
            +8.69157489e-10 * tc[2]
            -6.83568796e-14 * tc[3];
    }

    /*species with midpoint at T=1475 kelvin */
    if (T < 1475) {
        /*species 24: OCH2OCHO */
        species[24] =
            +1.58839723e-02
            +7.07081094e-07 * tc[1]
            -1.83137077e-08 * tc[2]
            +7.78647204e-12 * tc[3];
    } else {
        /*species 24: OCH2OCHO */
        species[24] =
            +8.11262659e-03
            -5.82712924e-06 * tc[1]
            +1.40202115e-09 * tc[2]
            -1.10950210e-13 * tc[3];
    }

    /*species with midpoint at T=710 kelvin */
    if (T < 710) {
        /*species 20: CH3OCH3 */
        species[20] =
            -5.39434751e-03
            +1.29894550e-04 * tc[1]
            -2.41519595e-07 * tc[2]
            +1.30989607e-10 * tc[3];
    } else {
        /*species 20: CH3OCH3 */
        species[20] =
            +2.69173263e-02
            -2.77749554e-05 * tc[1]
            +1.04254524e-08 * tc[2]
            -1.36682714e-12 * tc[3];
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 1: H2 */
        species[1] =
            +8.24944200e-04
            -1.62860300e-06 * tc[1]
            -2.84263020e-10 * tc[2]
            +1.65394880e-12 * tc[3];
        /*species 2: CH3 */
        species[2] =
            +2.12659790e-03
            +1.09167766e-05 * tc[1]
            -1.98543009e-08 * tc[2]
            +9.86282960e-12 * tc[3];
        /*species 3: O */
        species[3] =
            -1.63816600e-03
            +4.84206400e-06 * tc[1]
            -4.80852900e-09 * tc[2]
            +1.55627840e-12 * tc[3];
        /*species 4: CH4 */
        species[4] =
            +1.74766800e-02
            -5.56681800e-05 * tc[1]
            +9.14912400e-08 * tc[2]
            -4.89572400e-11 * tc[3];
        /*species 5: OH */
        species[5] =
            -3.22544939e-03
            +1.30552938e-05 * tc[1]
            -1.73956093e-08 * tc[2]
            +8.24949516e-12 * tc[3];
        /*species 6: H2O */
        species[6] =
            +3.47498200e-03
            -1.27093920e-05 * tc[1]
            +2.09057430e-08 * tc[2]
            -1.00263520e-11 * tc[3];
        /*species 7: C2H2 */
        species[7] =
            +2.33615629e-02
            -7.10343630e-05 * tc[1]
            +8.40457311e-08 * tc[2]
            -3.40029190e-11 * tc[3];
        /*species 8: CO */
        species[8] =
            +1.51194100e-03
            -7.76351000e-06 * tc[1]
            +1.67458320e-08 * tc[2]
            -9.89980400e-12 * tc[3];
        /*species 9: C2H4 */
        species[9] =
            -7.57052247e-03
            +1.14198058e-04 * tc[1]
            -2.07476626e-07 * tc[2]
            +1.07953749e-10 * tc[3];
        /*species 10: C2H5 */
        species[10] =
            -4.18336380e-03
            +9.94145400e-05 * tc[1]
            -1.79717622e-07 * tc[2]
            +9.21939120e-11 * tc[3];
        /*species 12: C2H6 */
        species[12] =
            +1.54946700e-02
            +1.15610140e-05 * tc[1]
            -3.77349600e-08 * tc[2]
            +1.83450680e-11 * tc[3];
        /*species 13: CH3O */
        species[13] =
            +7.21659500e-03
            +1.06769440e-05 * tc[1]
            -2.21329080e-08 * tc[2]
            +8.30244400e-12 * tc[3];
        /*species 14: O2 */
        species[14] =
            +1.12748600e-03
            -1.15123000e-06 * tc[1]
            +3.94163100e-09 * tc[2]
            -3.50742160e-12 * tc[3];
        /*species 15: HO2 */
        species[15] =
            -4.74912051e-03
            +4.23165782e-05 * tc[1]
            -7.28291682e-08 * tc[2]
            +3.71690050e-11 * tc[3];
        /*species 16: H2O2 */
        species[16] =
            +6.56922600e-03
            -2.97002600e-07 * tc[1]
            -1.38774180e-08 * tc[2]
            +9.88606000e-12 * tc[3];
        /*species 17: CO2 */
        species[17] =
            +9.92207200e-03
            -2.08182200e-05 * tc[1]
            +2.06000610e-08 * tc[2]
            -8.46912000e-12 * tc[3];
        /*species 18: CH3HCO */
        species[18] =
            -3.19328580e-03
            +9.50698420e-05 * tc[1]
            -1.72375833e-07 * tc[2]
            +8.77244480e-11 * tc[3];
        /*species 29: N2 */
        species[29] =
            +1.40824000e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77942000e-12 * tc[3];
    } else {
        /*species 0: H */
        species[0] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 1: H2 */
        species[1] =
            +7.00064400e-04
            -1.12676580e-07 * tc[1]
            -2.76947340e-11 * tc[2]
            +6.33100800e-15 * tc[3];
        /*species 2: CH3 */
        species[2] =
            +5.79785200e-03
            -3.95116000e-06 * tc[1]
            +9.21893700e-10 * tc[2]
            -7.16696640e-14 * tc[3];
        /*species 3: O */
        species[3] =
            -2.75506200e-05
            -6.20560600e-09 * tc[1]
            +1.36532010e-11 * tc[2]
            -1.74722080e-15 * tc[3];
        /*species 4: CH4 */
        species[4] =
            +1.02372400e-02
            -7.75025800e-06 * tc[1]
            +2.03567550e-09 * tc[2]
            -1.80136920e-13 * tc[3];
        /*species 5: OH */
        species[5] =
            +1.05650448e-03
            -5.18165516e-07 * tc[1]
            +9.15656022e-11 * tc[2]
            -5.32783504e-15 * tc[3];
        /*species 6: H2O */
        species[6] =
            +3.05629300e-03
            -1.74605200e-06 * tc[1]
            +3.60298800e-10 * tc[2]
            -2.55664720e-14 * tc[3];
        /*species 7: C2H2 */
        species[7] =
            +5.96166664e-03
            -4.74589704e-06 * tc[1]
            +1.40223651e-09 * tc[2]
            -1.44494085e-13 * tc[3];
        /*species 8: CO */
        species[8] =
            +1.44268900e-03
            -1.12616560e-06 * tc[1]
            +3.05574300e-10 * tc[2]
            -2.76438080e-14 * tc[3];
        /*species 9: C2H4 */
        species[9] =
            +1.46454151e-02
            -1.34215583e-05 * tc[1]
            +4.41668769e-09 * tc[2]
            -5.02824244e-13 * tc[3];
        /*species 10: C2H5 */
        species[10] =
            +1.24338930e-02
            -8.82782380e-06 * tc[1]
            +2.11962306e-09 * tc[2]
            -1.68140544e-13 * tc[3];
        /*species 12: C2H6 */
        species[12] =
            +1.38404300e-02
            -9.11451800e-06 * tc[1]
            +2.01749010e-09 * tc[2]
            -1.43926440e-13 * tc[3];
        /*species 13: CH3O */
        species[13] =
            +7.87149700e-03
            -5.31276800e-06 * tc[1]
            +1.18332930e-09 * tc[2]
            -8.45046400e-14 * tc[3];
        /*species 14: O2 */
        species[14] =
            +6.13519700e-04
            -2.51768400e-07 * tc[1]
            +5.32584300e-11 * tc[2]
            -4.54574000e-15 * tc[3];
        /*species 15: HO2 */
        species[15] =
            +2.23982013e-03
            -1.26731630e-06 * tc[1]
            +3.42739110e-10 * tc[2]
            -4.31634140e-14 * tc[3];
        /*species 16: H2O2 */
        species[16] =
            +4.33613600e-03
            -2.94937800e-06 * tc[1]
            +7.04671200e-10 * tc[2]
            -5.72661600e-14 * tc[3];
        /*species 17: CO2 */
        species[17] =
            +3.14016900e-03
            -2.55682200e-06 * tc[1]
            +7.18199100e-10 * tc[2]
            -6.67613200e-14 * tc[3];
        /*species 18: CH3HCO */
        species[18] =
            +1.17230590e-02
            -8.45262740e-06 * tc[1]
            +2.05117353e-09 * tc[2]
            -1.63939452e-13 * tc[3];
        /*species 29: N2 */
        species[29] =
            +1.48797700e-03
            -1.13695220e-06 * tc[1]
            +3.02911200e-10 * tc[2]
            -2.70134040e-14 * tc[3];
    }

    /*species with midpoint at T=1387 kelvin */
    if (T < 1387) {
        /*species 27: HO2CH2OCHO */
        species[27] =
            +4.02952392e-02
            -6.60218592e-05 * tc[1]
            +4.03080351e-08 * tc[2]
            -8.74406320e-12 * tc[3];
    } else {
        /*species 27: HO2CH2OCHO */
        species[27] =
            +8.52683511e-03
            -6.08227000e-06 * tc[1]
            +1.45679072e-09 * tc[2]
            -1.14926534e-13 * tc[3];
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 26: CH3OCH2O2 */
        species[26] =
            +3.68877454e-02
            -5.65123110e-05 * tc[1]
            +3.47191599e-08 * tc[2]
            -7.88521880e-12 * tc[3];
    } else {
        /*species 26: CH3OCH2O2 */
        species[26] =
            +1.18705986e-02
            -8.15813064e-06 * tc[1]
            +1.90593243e-09 * tc[2]
            -1.47771147e-13 * tc[3];
    }

    /*species with midpoint at T=1200 kelvin */
    if (T < 1200) {
        /*species 11: CH2O */
        species[11] =
            +4.92614230e-03
            +1.65652988e-06 * tc[1]
            -1.65114588e-09 * tc[2]
            -1.58441304e-12 * tc[3];
    } else {
        /*species 11: CH2O */
        species[11] =
            +2.86780160e-03
            -4.75652660e-07 * tc[1]
            -4.83339090e-10 * tc[2]
            +1.14266940e-13 * tc[3];
    }

    /*species with midpoint at T=1362 kelvin */
    if (T < 1362) {
        /*species 21: CH3OCO */
        species[21] =
            +2.43434884e-02
            -3.31191120e-05 * tc[1]
            +1.37561223e-08 * tc[2]
            -1.32718283e-12 * tc[3];
    } else {
        /*species 21: CH3OCO */
        species[21] =
            +4.53544950e-03
            -3.30192728e-06 * tc[1]
            +8.01591831e-10 * tc[2]
            -6.38307452e-14 * tc[3];
    }

    /*species with midpoint at T=1603 kelvin */
    if (T < 1603) {
        /*species 25: HOCH2OCO */
        species[25] =
            +1.28768359e-02
            +4.08838836e-06 * tc[1]
            -1.83046476e-08 * tc[2]
            +7.19282236e-12 * tc[3];
    } else {
        /*species 25: HOCH2OCO */
        species[25] =
            +8.17663898e-03
            -5.84068042e-06 * tc[1]
            +1.40008685e-09 * tc[2]
            -1.10510729e-13 * tc[3];
    }

    /*species with midpoint at T=1686 kelvin */
    if (T < 1686) {
        /*species 22: CH3OCHO */
        species[22] =
            +2.03760048e-02
            -1.36955408e-05 * tc[1]
            -2.18455861e-09 * tc[2]
            +2.24852086e-12 * tc[3];
    } else {
        /*species 22: CH3OCHO */
        species[22] =
            +1.15503122e-02
            -8.55564972e-06 * tc[1]
            +2.10759918e-09 * tc[2]
            -1.69733421e-13 * tc[3];
    }

    /*species with midpoint at T=1402 kelvin */
    if (T < 1402) {
        /*species 28: O2CH2OCH2O2H */
        species[28] =
            +5.83226232e-02
            -1.10651956e-04 * tc[1]
            +7.79431620e-08 * tc[2]
            -1.90856402e-11 * tc[3];
    } else {
        /*species 28: O2CH2OCH2O2H */
        species[28] =
            +1.04394841e-02
            -7.21165878e-06 * tc[1]
            +1.69137853e-09 * tc[2]
            -1.31522886e-13 * tc[3];
    }

    /*species with midpoint at T=2014 kelvin */
    if (T < 2014) {
        /*species 23: CH3OCH2OH */
        species[23] =
            +2.44325751e-02
            -1.73396957e-05 * tc[1]
            -1.77995798e-10 * tc[2]
            +1.74560001e-12 * tc[3];
    } else {
        /*species 23: CH3OCH2OH */
        species[23] =
            +1.53602372e-02
            -1.08200758e-05 * tc[1]
            +2.58172034e-09 * tc[2]
            -2.03527901e-13 * tc[3];
    }
    return;
}


/*compute the progress rate for each reaction */
void progressRate(double * restrict  qdot, double * restrict  sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[30];                /*Gibbs free energy */
    double Kc;                      /*equilibrium constant */
    double k_f;                     /*forward reaction rate */
    double k_r;                     /*reverse reaction rate */
    double q_f;                     /*forward progress rate */
    double q_r;                     /*reverse progress rate */
    double phi_f;                   /*forward phase space factor */
    double phi_r;                   /*reverse phase space factor */
    double alpha;                   /*enhancement */
    double redP;                    /*reduced pressure */
    double logPred;                 /*log of above */
    double F;                       /*fallof rate enhancement */

    double F_troe;                  /*TROE intermediate */
    double logFcent;                /*TROE intermediate */
    double troe;                    /*TROE intermediate */
    double troe_c;                  /*TROE intermediate */
    double troe_n;                  /*TROE intermediate */

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */

    double invT = 1.0 / tc[1];

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;
    double refCinv = 1 / refC;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 30; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    if (T != T_save)
    {
        T_save = T;

        k_f_save[0] = 1e-06 * 0;
        k_f_save[1] = 1e-06 * 0;
        k_f_save[2] = 1e-06 * 0;

        Kc_save[0] = exp((g_RT[0] + g_RT[14]) - (g_RT[3] + g_RT[5]));
        Kc_save[1] = exp((g_RT[3] + g_RT[1]) - (g_RT[0] + g_RT[5]));
        Kc_save[2] = exp((g_RT[1] + g_RT[5]) - (g_RT[6] + g_RT[0]));
    }

    /*reaction 1: H + O2 <=> O + OH */
    phi_f = sc[0]*sc[14];
    k_f = k_f_save[0];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = Kc_save[0];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[0] = q_f - q_r;

    /*reaction 2: O + H2 <=> H + OH */
    phi_f = sc[3]*sc[1];
    k_f = k_f_save[1];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[5];
    Kc = Kc_save[1];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[1] = q_f - q_r;

    /*reaction 3: H2 + OH <=> H2O + H */
    phi_f = sc[1]*sc[5];
    k_f = k_f_save[2];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = Kc_save[2];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[2] = q_f - q_r;

    return;
}


/*compute the progress rate for each reaction */
void progressRateFR(double * restrict  q_f, double * restrict  q_r, double * restrict  sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[30];                /*Gibbs free energy */
    double Kc;                      /*equilibrium constant */
    double k_f;                     /*forward reaction rate */
    double k_r;                     /*reverse reaction rate */
    double phi_f;                   /*forward phase space factor */
    double phi_r;                   /*reverse phase space factor */
    double alpha;                   /*enhancement */
    double redP;                    /*reduced pressure */
    double logPred;                 /*log of above */
    double F;                       /*fallof rate enhancement */

    double F_troe;                  /*TROE intermediate */
    double logFcent;                /*TROE intermediate */
    double troe;                    /*TROE intermediate */
    double troe_c;                  /*TROE intermediate */
    double troe_n;                  /*TROE intermediate */

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */

    double invT = 1.0 / tc[1];

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 30; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    if (T != T_save)
    {
        T_save = T;

        k_f_save[0] = 1e-06 * 0;
        k_f_save[1] = 1e-06 * 0;
        k_f_save[2] = 1e-06 * 0;

        Kc_save[0] = exp((g_RT[0] + g_RT[14]) - (g_RT[3] + g_RT[5]));
        Kc_save[1] = exp((g_RT[3] + g_RT[1]) - (g_RT[0] + g_RT[5]));
        Kc_save[2] = exp((g_RT[1] + g_RT[5]) - (g_RT[6] + g_RT[0]));
    }

    /*reaction 1: H + O2 <=> O + OH */
    phi_f = sc[0]*sc[14];
    k_f = k_f_save[0];
    q_f[0] = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = Kc_save[0];
    k_r = k_f / Kc;
    q_r[0] = phi_r * k_r;

    /*reaction 2: O + H2 <=> H + OH */
    phi_f = sc[3]*sc[1];
    k_f = k_f_save[1];
    q_f[1] = phi_f * k_f;
    phi_r = sc[0]*sc[5];
    Kc = Kc_save[1];
    k_r = k_f / Kc;
    q_r[1] = phi_r * k_r;

    /*reaction 3: H2 + OH <=> H2O + H */
    phi_f = sc[1]*sc[5];
    k_f = k_f_save[2];
    q_f[2] = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = Kc_save[2];
    k_r = k_f / Kc;
    q_r[2] = phi_r * k_r;

    return;
}


/*compute the equilibrium constants for each reaction */
void equilibriumConstants(double * restrict kc, double * restrict  g_RT, double T)
{
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    /*reaction 1: H + O2 <=> O + OH */
    kc[0] = exp((g_RT[0] + g_RT[14]) - (g_RT[3] + g_RT[5]));

    /*reaction 2: O + H2 <=> H + OH */
    kc[1] = exp((g_RT[3] + g_RT[1]) - (g_RT[0] + g_RT[5]));

    /*reaction 3: H2 + OH <=> H2O + H */
    kc[2] = exp((g_RT[1] + g_RT[5]) - (g_RT[6] + g_RT[0]));

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void gibbs(double * restrict  species, double * restrict  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1376 kelvin */
    if (T < 1376) {
        /*species 19: HCOOH */
        species[19] =
            -4.646165040000000e+04 * invT
            -1.585309795000000e+01
            -1.435481850000000e+00 * tc[0]
            -8.168150799999999e-03 * tc[1]
            +1.770957016666667e-06 * tc[2]
            -2.767774808333333e-10 * tc[3]
            +2.010880515000000e-14 * tc[4];
    } else {
        /*species 19: HCOOH */
        species[19] =
            -4.839954000000000e+04 * invT
            +1.799780993000000e+01
            -6.687330130000000e+00 * tc[0]
            -2.571446840000000e-03 * tc[1]
            +3.037308550000000e-07 * tc[2]
            -2.414326358333334e-11 * tc[3]
            +8.544609950000001e-16 * tc[4];
    }

    /*species with midpoint at T=1475 kelvin */
    if (T < 1475) {
        /*species 24: OCH2OCHO */
        species[24] =
            -4.022427920000000e+04 * invT
            -9.195499099999997e-01
            -5.196908370000000e+00 * tc[0]
            -7.941986150000000e-03 * tc[1]
            -5.892342450000000e-08 * tc[2]
            +5.087141025000000e-10 * tc[3]
            -9.733090050000001e-14 * tc[4];
    } else {
        /*species 24: OCH2OCHO */
        species[24] =
            -4.336472310000000e+04 * invT
            +4.539257250000000e+01
            -1.202339160000000e+01 * tc[0]
            -4.056313295000000e-03 * tc[1]
            +4.855941033333334e-07 * tc[2]
            -3.894503200000000e-11 * tc[3]
            +1.386877625000000e-15 * tc[4];
    }

    /*species with midpoint at T=710 kelvin */
    if (T < 710) {
        /*species 20: CH3OCH3 */
        species[20] =
            -2.397554550000000e+04 * invT
            +6.317929965999999e+00
            -5.680974470000000e+00 * tc[0]
            +2.697173755000000e-03 * tc[1]
            -1.082454583333333e-05 * tc[2]
            +6.708877650000000e-09 * tc[3]
            -1.637370090000000e-12 * tc[4];
    } else {
        /*species 20: CH3OCH3 */
        species[20] =
            -2.341209750000000e+04 * invT
            -1.938662045400000e+01
            -8.308155460000000e-01 * tc[0]
            -1.345866315000000e-02 * tc[1]
            +2.314579616666667e-06 * tc[2]
            -2.895958991666667e-10 * tc[3]
            +1.708533920000000e-14 * tc[4];
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.547163000000000e+04 * invT
            +2.960117600000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -1.012521000000000e+03 * invT
            +6.592218000000000e+00
            -3.298124000000000e+00 * tc[0]
            -4.124721000000000e-04 * tc[1]
            +1.357169166666667e-07 * tc[2]
            +7.896194999999999e-12 * tc[3]
            -2.067436000000000e-14 * tc[4];
        /*species 2: CH3 */
        species[2] =
            +1.642271600000000e+04 * invT
            +1.983644300000000e+00
            -3.657179700000000e+00 * tc[0]
            -1.063298950000000e-03 * tc[1]
            -9.097313833333333e-07 * tc[2]
            +5.515083583333334e-10 * tc[3]
            -1.232853700000000e-13 * tc[4];
        /*species 3: O */
        species[3] =
            +2.914764000000000e+04 * invT
            -1.756599999999997e-02
            -2.946429000000000e+00 * tc[0]
            +8.190830000000000e-04 * tc[1]
            -4.035053333333333e-07 * tc[2]
            +1.335702500000000e-10 * tc[3]
            -1.945348000000000e-14 * tc[4];
        /*species 4: CH4 */
        species[4] =
            -9.825228999999999e+03 * invT
            -1.294344850000000e+01
            -7.787415000000000e-01 * tc[0]
            -8.738340000000001e-03 * tc[1]
            +4.639015000000000e-06 * tc[2]
            -2.541423333333333e-09 * tc[3]
            +6.119655000000000e-13 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.346309130000000e+03 * invT
            +4.815738570000000e+00
            -4.125305610000000e+00 * tc[0]
            +1.612724695000000e-03 * tc[1]
            -1.087941151666667e-06 * tc[2]
            +4.832113691666666e-10 * tc[3]
            -1.031186895000000e-13 * tc[4];
        /*species 6: H2O */
        species[6] =
            -3.020811000000000e+04 * invT
            +7.966090000000001e-01
            -3.386842000000000e+00 * tc[0]
            -1.737491000000000e-03 * tc[1]
            +1.059116000000000e-06 * tc[2]
            -5.807150833333333e-10 * tc[3]
            +1.253294000000000e-13 * tc[4];
        /*species 7: C2H2 */
        species[7] =
            +2.642898070000000e+04 * invT
            -1.313102400600000e+01
            -8.086810940000000e-01 * tc[0]
            -1.168078145000000e-02 * tc[1]
            +5.919530250000000e-06 * tc[2]
            -2.334603641666667e-09 * tc[3]
            +4.250364870000000e-13 * tc[4];
        /*species 8: CO */
        species[8] =
            -1.431054000000000e+04 * invT
            -1.586445000000000e+00
            -3.262452000000000e+00 * tc[0]
            -7.559705000000000e-04 * tc[1]
            +6.469591666666667e-07 * tc[2]
            -4.651620000000000e-10 * tc[3]
            +1.237475500000000e-13 * tc[4];
        /*species 9: C2H4 */
        species[9] =
            +5.089775930000000e+03 * invT
            -1.381294799999999e-01
            -3.959201480000000e+00 * tc[0]
            +3.785261235000000e-03 * tc[1]
            -9.516504866666667e-06 * tc[2]
            +5.763239608333333e-09 * tc[3]
            -1.349421865000000e-12 * tc[4];
        /*species 10: C2H5 */
        species[10] =
            +1.284171400000000e+04 * invT
            -4.041656000000007e-01
            -4.305858000000000e+00 * tc[0]
            +2.091681900000000e-03 * tc[1]
            -8.284545000000000e-06 * tc[2]
            +4.992156166666666e-09 * tc[3]
            -1.152423900000000e-12 * tc[4];
        /*species 12: C2H6 */
        species[12] =
            -1.123918000000000e+04 * invT
            -1.296975100000000e+01
            -1.462539000000000e+00 * tc[0]
            -7.747335000000000e-03 * tc[1]
            -9.634178333333333e-07 * tc[2]
            +1.048193333333333e-09 * tc[3]
            -2.293133500000000e-13 * tc[4];
        /*species 13: CH3O */
        species[13] =
            +9.786011000000000e+02 * invT
            -1.104597600000000e+01
            -2.106204000000000e+00 * tc[0]
            -3.608297500000000e-03 * tc[1]
            -8.897453333333333e-07 * tc[2]
            +6.148030000000000e-10 * tc[3]
            -1.037805500000000e-13 * tc[4];
        /*species 14: O2 */
        species[14] =
            -1.005249000000000e+03 * invT
            -2.821802000000000e+00
            -3.212936000000000e+00 * tc[0]
            -5.637430000000000e-04 * tc[1]
            +9.593583333333333e-08 * tc[2]
            -1.094897500000000e-10 * tc[3]
            +4.384277000000000e-14 * tc[4];
        /*species 15: HO2 */
        species[15] =
            +2.948080400000000e+02 * invT
            +5.851355599999999e-01
            -4.301798010000000e+00 * tc[0]
            +2.374560255000000e-03 * tc[1]
            -3.526381516666666e-06 * tc[2]
            +2.023032450000000e-09 * tc[3]
            -4.646125620000001e-13 * tc[4];
        /*species 16: H2O2 */
        species[16] =
            -1.766315000000000e+04 * invT
            -3.396609000000000e+00
            -3.388754000000000e+00 * tc[0]
            -3.284613000000000e-03 * tc[1]
            +2.475021666666666e-08 * tc[2]
            +3.854838333333333e-10 * tc[3]
            -1.235757500000000e-13 * tc[4];
        /*species 17: CO2 */
        species[17] =
            -4.837314000000000e+04 * invT
            -7.912765000000000e+00
            -2.275725000000000e+00 * tc[0]
            -4.961036000000000e-03 * tc[1]
            +1.734851666666667e-06 * tc[2]
            -5.722239166666667e-10 * tc[3]
            +1.058640000000000e-13 * tc[4];
        /*species 18: CH3HCO */
        species[18] =
            -2.157287800000000e+04 * invT
            +6.264436000000000e-01
            -4.729459500000000e+00 * tc[0]
            +1.596642900000000e-03 * tc[1]
            -7.922486833333334e-06 * tc[2]
            +4.788217583333333e-09 * tc[3]
            -1.096555600000000e-12 * tc[4];
        /*species 29: N2 */
        species[29] =
            -1.020900000000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041200000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427500000000e-13 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.547163000000000e+04 * invT
            +2.960117600000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -8.350340000000000e+02 * invT
            +4.346533000000000e+00
            -2.991423000000000e+00 * tc[0]
            -3.500322000000000e-04 * tc[1]
            +9.389715000000000e-09 * tc[2]
            +7.692981666666667e-13 * tc[3]
            -7.913760000000000e-17 * tc[4];
        /*species 2: CH3 */
        species[2] =
            +1.650951300000000e+04 * invT
            -1.744359300000000e+00
            -2.978120600000000e+00 * tc[0]
            -2.898926000000000e-03 * tc[1]
            +3.292633333333333e-07 * tc[2]
            -2.560815833333334e-11 * tc[3]
            +8.958708000000000e-16 * tc[4];
        /*species 3: O */
        species[3] =
            +2.923080000000000e+04 * invT
            -2.378248000000000e+00
            -2.542060000000000e+00 * tc[0]
            +1.377531000000000e-05 * tc[1]
            +5.171338333333333e-10 * tc[2]
            -3.792555833333334e-13 * tc[3]
            +2.184026000000000e-17 * tc[4];
        /*species 4: CH4 */
        species[4] =
            -1.008079000000000e+04 * invT
            -7.939916000000000e+00
            -1.683479000000000e+00 * tc[0]
            -5.118620000000000e-03 * tc[1]
            +6.458548333333333e-07 * tc[2]
            -5.654654166666667e-11 * tc[3]
            +2.251711500000000e-15 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.683628750000000e+03 * invT
            -2.836911870000000e+00
            -2.864728860000000e+00 * tc[0]
            -5.282522400000000e-04 * tc[1]
            +4.318045966666667e-08 * tc[2]
            -2.543488950000000e-12 * tc[3]
            +6.659793800000000e-17 * tc[4];
        /*species 6: H2O */
        species[6] =
            -2.989921000000000e+04 * invT
            -4.190671000000000e+00
            -2.672146000000000e+00 * tc[0]
            -1.528146500000000e-03 * tc[1]
            +1.455043333333333e-07 * tc[2]
            -1.000830000000000e-11 * tc[3]
            +3.195809000000000e-16 * tc[4];
        /*species 7: C2H2 */
        species[7] =
            +2.593599920000000e+04 * invT
            +5.377850850000001e+00
            -4.147569640000000e+00 * tc[0]
            -2.980833320000000e-03 * tc[1]
            +3.954914200000000e-07 * tc[2]
            -3.895101425000000e-11 * tc[3]
            +1.806176065000000e-15 * tc[4];
        /*species 8: CO */
        species[8] =
            -1.426835000000000e+04 * invT
            -3.083140000000000e+00
            -3.025078000000000e+00 * tc[0]
            -7.213445000000000e-04 * tc[1]
            +9.384713333333334e-08 * tc[2]
            -8.488174999999999e-12 * tc[3]
            +3.455476000000000e-16 * tc[4];
        /*species 9: C2H4 */
        species[9] =
            +4.939886140000000e+03 * invT
            -8.269258140000002e+00
            -2.036111160000000e+00 * tc[0]
            -7.322707550000000e-03 * tc[1]
            +1.118463191666667e-06 * tc[2]
            -1.226857691666667e-10 * tc[3]
            +6.285303050000000e-15 * tc[4];
        /*species 10: C2H5 */
        species[10] =
            +1.205645500000000e+04 * invT
            +3.441855570000000e+00
            -4.287881400000000e+00 * tc[0]
            -6.216946500000000e-03 * tc[1]
            +7.356519833333333e-07 * tc[2]
            -5.887841833333334e-11 * tc[3]
            +2.101756800000000e-15 * tc[4];
        /*species 12: C2H6 */
        species[12] =
            -1.271779000000000e+04 * invT
            +1.006544500000000e+01
            -4.825938000000000e+00 * tc[0]
            -6.920215000000000e-03 * tc[1]
            +7.595431666666667e-07 * tc[2]
            -5.604139166666666e-11 * tc[3]
            +1.799080500000000e-15 * tc[4];
        /*species 13: CH3O */
        species[13] =
            +1.278325000000000e+02 * invT
            +8.412250000000001e-01
            -3.770800000000000e+00 * tc[0]
            -3.935748500000000e-03 * tc[1]
            +4.427306666666667e-07 * tc[2]
            -3.287025833333333e-11 * tc[3]
            +1.056308000000000e-15 * tc[4];
        /*species 14: O2 */
        species[14] =
            -1.233930000000000e+03 * invT
            +5.084119999999999e-01
            -3.697578000000000e+00 * tc[0]
            -3.067598500000000e-04 * tc[1]
            +2.098070000000000e-08 * tc[2]
            -1.479400833333333e-12 * tc[3]
            +5.682175000000001e-17 * tc[4];
        /*species 15: HO2 */
        species[15] =
            +1.118567130000000e+02 * invT
            +2.321087500000001e-01
            -4.017210900000000e+00 * tc[0]
            -1.119910065000000e-03 * tc[1]
            +1.056096916666667e-07 * tc[2]
            -9.520530833333334e-12 * tc[3]
            +5.395426750000000e-16 * tc[4];
        /*species 16: H2O2 */
        species[16] =
            -1.800696000000000e+04 * invT
            +4.072030000000000e+00
            -4.573167000000000e+00 * tc[0]
            -2.168068000000000e-03 * tc[1]
            +2.457815000000000e-07 * tc[2]
            -1.957420000000000e-11 * tc[3]
            +7.158270000000000e-16 * tc[4];
        /*species 17: CO2 */
        species[17] =
            -4.896696000000000e+04 * invT
            +5.409018900000000e+00
            -4.453623000000000e+00 * tc[0]
            -1.570084500000000e-03 * tc[1]
            +2.130685000000000e-07 * tc[2]
            -1.994997500000000e-11 * tc[3]
            +8.345165000000000e-16 * tc[4];
        /*species 18: CH3HCO */
        species[18] =
            -2.259312200000000e+04 * invT
            +8.884902499999999e+00
            -5.404110800000000e+00 * tc[0]
            -5.861529500000000e-03 * tc[1]
            +7.043856166666666e-07 * tc[2]
            -5.697704250000000e-11 * tc[3]
            +2.049243150000000e-15 * tc[4];
        /*species 29: N2 */
        species[29] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439885000000000e-04 * tc[1]
            +9.474601666666666e-08 * tc[2]
            -8.414199999999999e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
    }

    /*species with midpoint at T=1387 kelvin */
    if (T < 1387) {
        /*species 27: HO2CH2OCHO */
        species[27] =
            -5.806299340000000e+04 * invT
            -1.177278217000000e+01
            -3.479357030000000e+00 * tc[0]
            -2.014761960000000e-02 * tc[1]
            +5.501821600000000e-06 * tc[2]
            -1.119667641666667e-09 * tc[3]
            +1.093007900000000e-13 * tc[4];
    } else {
        /*species 27: HO2CH2OCHO */
        species[27] =
            -6.239596080000000e+04 * invT
            +7.035084370000000e+01
            -1.645842980000000e+01 * tc[0]
            -4.263417555000000e-03 * tc[1]
            +5.068558333333334e-07 * tc[2]
            -4.046640900000000e-11 * tc[3]
            +1.436581670000000e-15 * tc[4];
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 26: CH3OCH2O2 */
        species[26] =
            -1.949409400000000e+04 * invT
            -1.693606398000000e+01
            -2.210296120000000e+00 * tc[0]
            -1.844387270000000e-02 * tc[1]
            +4.709359250000000e-06 * tc[2]
            -9.644211083333333e-10 * tc[3]
            +9.856523500000001e-14 * tc[4];
    } else {
        /*species 26: CH3OCH2O2 */
        species[26] =
            -2.296792380000000e+04 * invT
            +4.779898740000000e+01
            -1.242497290000000e+01 * tc[0]
            -5.935299300000000e-03 * tc[1]
            +6.798442200000000e-07 * tc[2]
            -5.294256741666666e-11 * tc[3]
            +1.847139335000000e-15 * tc[4];
    }

    /*species with midpoint at T=1200 kelvin */
    if (T < 1200) {
        /*species 11: CH2O */
        species[11] =
            -1.497079300000000e+04 * invT
            -6.773498699999999e+00
            -2.696261200000000e+00 * tc[0]
            -2.463071150000000e-03 * tc[1]
            -1.380441566666667e-07 * tc[2]
            +4.586516333333333e-11 * tc[3]
            +1.980516300000000e-14 * tc[4];
    } else {
        /*species 11: CH2O */
        species[11] =
            -1.623017300000000e+04 * invT
            +1.026957180000000e+01
            -5.148190500000000e+00 * tc[0]
            -1.433900800000000e-03 * tc[1]
            +3.963772166666667e-08 * tc[2]
            +1.342608583333333e-11 * tc[3]
            -1.428336750000000e-15 * tc[4];
    }

    /*species with midpoint at T=1362 kelvin */
    if (T < 1362) {
        /*species 21: CH3OCO */
        species[21] =
            -2.144048290000000e+04 * invT
            -1.275344461000000e+01
            -3.941991590000000e+00 * tc[0]
            -1.217174420000000e-02 * tc[1]
            +2.759926000000000e-06 * tc[2]
            -3.821145091666667e-10 * tc[3]
            +1.658978540000000e-14 * tc[4];
    } else {
        /*species 21: CH3OCO */
        species[21] =
            -2.466164000000000e+04 * invT
            +4.587916509999999e+01
            -1.308776000000000e+01 * tc[0]
            -2.267724750000000e-03 * tc[1]
            +2.751606066666667e-07 * tc[2]
            -2.226643975000000e-11 * tc[3]
            +7.978843150000000e-16 * tc[4];
    }

    /*species with midpoint at T=1603 kelvin */
    if (T < 1603) {
        /*species 25: HOCH2OCO */
        species[25] =
            -4.395261830000000e+04 * invT
            +3.541263520000000e+00
            -6.081808010000000e+00 * tc[0]
            -6.438417950000000e-03 * tc[1]
            -3.406990300000000e-07 * tc[2]
            +5.084624341666666e-10 * tc[3]
            -8.991027950000001e-14 * tc[4];
    } else {
        /*species 25: HOCH2OCO */
        species[25] =
            -4.655757430000000e+04 * invT
            +3.997726560000000e+01
            -1.137373910000000e+01 * tc[0]
            -4.088319490000000e-03 * tc[1]
            +4.867233683333334e-07 * tc[2]
            -3.889130133333334e-11 * tc[3]
            +1.381384115000000e-15 * tc[4];
    }

    /*species with midpoint at T=1686 kelvin */
    if (T < 1686) {
        /*species 22: CH3OCHO */
        species[22] =
            -4.418551670000000e+04 * invT
            -9.448074070000001e+00
            -3.088397830000000e+00 * tc[0]
            -1.018800240000000e-02 * tc[1]
            +1.141295066666667e-06 * tc[2]
            +6.068218358333333e-11 * tc[3]
            -2.810651080000000e-14 * tc[4];
    } else {
        /*species 22: CH3OCHO */
        species[22] =
            -4.643647690000000e+04 * invT
            +2.762138298000000e+01
            -8.691235180000000e+00 * tc[0]
            -5.775156100000000e-03 * tc[1]
            +7.129708100000000e-07 * tc[2]
            -5.854442158333334e-11 * tc[3]
            +2.121667760000000e-15 * tc[4];
    }

    /*species with midpoint at T=1402 kelvin */
    if (T < 1402) {
        /*species 28: O2CH2OCH2O2H */
        species[28] =
            -3.276287420000000e+04 * invT
            -2.242509499000000e+01
            -1.996405510000000e+00 * tc[0]
            -2.916131160000000e-02 * tc[1]
            +9.220996300000001e-06 * tc[2]
            -2.165087833333333e-09 * tc[3]
            +2.385705025000000e-13 * tc[4];
    } else {
        /*species 28: O2CH2OCH2O2H */
        species[28] =
            -3.792070550000000e+04 * invT
            +8.438853190000000e+01
            -1.920380460000000e+01 * tc[0]
            -5.219742050000000e-03 * tc[1]
            +6.009715650000000e-07 * tc[2]
            -4.698273691666667e-11 * tc[3]
            +1.644036070000000e-15 * tc[4];
    }

    /*species with midpoint at T=2014 kelvin */
    if (T < 2014) {
        /*species 23: CH3OCH2OH */
        species[23] =
            -4.544888990000000e+04 * invT
            -9.892604739999999e+00
            -3.158518760000000e+00 * tc[0]
            -1.221628755000000e-02 * tc[1]
            +1.444974640000000e-06 * tc[2]
            +4.944327733333333e-12 * tc[3]
            -2.182000015000000e-14 * tc[4];
    } else {
        /*species 23: CH3OCH2OH */
        species[23] =
            -4.766071150000000e+04 * invT
            +2.673248590000000e+01
            -8.709815700000000e+00 * tc[0]
            -7.680118600000000e-03 * tc[1]
            +9.016729800000000e-07 * tc[2]
            -7.171445383333333e-11 * tc[3]
            +2.544098760000000e-15 * tc[4];
    }
    return;
}


/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void helmholtz(double * restrict  species, double * restrict  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1376 kelvin */
    if (T < 1376) {
        /*species 19: HCOOH */
        species[19] =
            -4.64616504e+04 * invT
            -1.68530979e+01
            -1.43548185e+00 * tc[0]
            -8.16815080e-03 * tc[1]
            +1.77095702e-06 * tc[2]
            -2.76777481e-10 * tc[3]
            +2.01088051e-14 * tc[4];
    } else {
        /*species 19: HCOOH */
        species[19] =
            -4.83995400e+04 * invT
            +1.69978099e+01
            -6.68733013e+00 * tc[0]
            -2.57144684e-03 * tc[1]
            +3.03730855e-07 * tc[2]
            -2.41432636e-11 * tc[3]
            +8.54460995e-16 * tc[4];
    }

    /*species with midpoint at T=1475 kelvin */
    if (T < 1475) {
        /*species 24: OCH2OCHO */
        species[24] =
            -4.02242792e+04 * invT
            -1.91954991e+00
            -5.19690837e+00 * tc[0]
            -7.94198615e-03 * tc[1]
            -5.89234245e-08 * tc[2]
            +5.08714103e-10 * tc[3]
            -9.73309005e-14 * tc[4];
    } else {
        /*species 24: OCH2OCHO */
        species[24] =
            -4.33647231e+04 * invT
            +4.43925725e+01
            -1.20233916e+01 * tc[0]
            -4.05631329e-03 * tc[1]
            +4.85594103e-07 * tc[2]
            -3.89450320e-11 * tc[3]
            +1.38687762e-15 * tc[4];
    }

    /*species with midpoint at T=710 kelvin */
    if (T < 710) {
        /*species 20: CH3OCH3 */
        species[20] =
            -2.39755455e+04 * invT
            +5.31792997e+00
            -5.68097447e+00 * tc[0]
            +2.69717376e-03 * tc[1]
            -1.08245458e-05 * tc[2]
            +6.70887765e-09 * tc[3]
            -1.63737009e-12 * tc[4];
    } else {
        /*species 20: CH3OCH3 */
        species[20] =
            -2.34120975e+04 * invT
            -2.03866205e+01
            -8.30815546e-01 * tc[0]
            -1.34586631e-02 * tc[1]
            +2.31457962e-06 * tc[2]
            -2.89595899e-10 * tc[3]
            +1.70853392e-14 * tc[4];
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.54716300e+04 * invT
            +1.96011760e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -1.01252100e+03 * invT
            +5.59221800e+00
            -3.29812400e+00 * tc[0]
            -4.12472100e-04 * tc[1]
            +1.35716917e-07 * tc[2]
            +7.89619500e-12 * tc[3]
            -2.06743600e-14 * tc[4];
        /*species 2: CH3 */
        species[2] =
            +1.64227160e+04 * invT
            +9.83644300e-01
            -3.65717970e+00 * tc[0]
            -1.06329895e-03 * tc[1]
            -9.09731383e-07 * tc[2]
            +5.51508358e-10 * tc[3]
            -1.23285370e-13 * tc[4];
        /*species 3: O */
        species[3] =
            +2.91476400e+04 * invT
            -1.01756600e+00
            -2.94642900e+00 * tc[0]
            +8.19083000e-04 * tc[1]
            -4.03505333e-07 * tc[2]
            +1.33570250e-10 * tc[3]
            -1.94534800e-14 * tc[4];
        /*species 4: CH4 */
        species[4] =
            -9.82522900e+03 * invT
            -1.39434485e+01
            -7.78741500e-01 * tc[0]
            -8.73834000e-03 * tc[1]
            +4.63901500e-06 * tc[2]
            -2.54142333e-09 * tc[3]
            +6.11965500e-13 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.34630913e+03 * invT
            +3.81573857e+00
            -4.12530561e+00 * tc[0]
            +1.61272470e-03 * tc[1]
            -1.08794115e-06 * tc[2]
            +4.83211369e-10 * tc[3]
            -1.03118689e-13 * tc[4];
        /*species 6: H2O */
        species[6] =
            -3.02081100e+04 * invT
            -2.03391000e-01
            -3.38684200e+00 * tc[0]
            -1.73749100e-03 * tc[1]
            +1.05911600e-06 * tc[2]
            -5.80715083e-10 * tc[3]
            +1.25329400e-13 * tc[4];
        /*species 7: C2H2 */
        species[7] =
            +2.64289807e+04 * invT
            -1.41310240e+01
            -8.08681094e-01 * tc[0]
            -1.16807815e-02 * tc[1]
            +5.91953025e-06 * tc[2]
            -2.33460364e-09 * tc[3]
            +4.25036487e-13 * tc[4];
        /*species 8: CO */
        species[8] =
            -1.43105400e+04 * invT
            -2.58644500e+00
            -3.26245200e+00 * tc[0]
            -7.55970500e-04 * tc[1]
            +6.46959167e-07 * tc[2]
            -4.65162000e-10 * tc[3]
            +1.23747550e-13 * tc[4];
        /*species 9: C2H4 */
        species[9] =
            +5.08977593e+03 * invT
            -1.13812948e+00
            -3.95920148e+00 * tc[0]
            +3.78526124e-03 * tc[1]
            -9.51650487e-06 * tc[2]
            +5.76323961e-09 * tc[3]
            -1.34942187e-12 * tc[4];
        /*species 10: C2H5 */
        species[10] =
            +1.28417140e+04 * invT
            -1.40416560e+00
            -4.30585800e+00 * tc[0]
            +2.09168190e-03 * tc[1]
            -8.28454500e-06 * tc[2]
            +4.99215617e-09 * tc[3]
            -1.15242390e-12 * tc[4];
        /*species 12: C2H6 */
        species[12] =
            -1.12391800e+04 * invT
            -1.39697510e+01
            -1.46253900e+00 * tc[0]
            -7.74733500e-03 * tc[1]
            -9.63417833e-07 * tc[2]
            +1.04819333e-09 * tc[3]
            -2.29313350e-13 * tc[4];
        /*species 13: CH3O */
        species[13] =
            +9.78601100e+02 * invT
            -1.20459760e+01
            -2.10620400e+00 * tc[0]
            -3.60829750e-03 * tc[1]
            -8.89745333e-07 * tc[2]
            +6.14803000e-10 * tc[3]
            -1.03780550e-13 * tc[4];
        /*species 14: O2 */
        species[14] =
            -1.00524900e+03 * invT
            -3.82180200e+00
            -3.21293600e+00 * tc[0]
            -5.63743000e-04 * tc[1]
            +9.59358333e-08 * tc[2]
            -1.09489750e-10 * tc[3]
            +4.38427700e-14 * tc[4];
        /*species 15: HO2 */
        species[15] =
            +2.94808040e+02 * invT
            -4.14864440e-01
            -4.30179801e+00 * tc[0]
            +2.37456025e-03 * tc[1]
            -3.52638152e-06 * tc[2]
            +2.02303245e-09 * tc[3]
            -4.64612562e-13 * tc[4];
        /*species 16: H2O2 */
        species[16] =
            -1.76631500e+04 * invT
            -4.39660900e+00
            -3.38875400e+00 * tc[0]
            -3.28461300e-03 * tc[1]
            +2.47502167e-08 * tc[2]
            +3.85483833e-10 * tc[3]
            -1.23575750e-13 * tc[4];
        /*species 17: CO2 */
        species[17] =
            -4.83731400e+04 * invT
            -8.91276500e+00
            -2.27572500e+00 * tc[0]
            -4.96103600e-03 * tc[1]
            +1.73485167e-06 * tc[2]
            -5.72223917e-10 * tc[3]
            +1.05864000e-13 * tc[4];
        /*species 18: CH3HCO */
        species[18] =
            -2.15728780e+04 * invT
            -3.73556400e-01
            -4.72945950e+00 * tc[0]
            +1.59664290e-03 * tc[1]
            -7.92248683e-06 * tc[2]
            +4.78821758e-09 * tc[3]
            -1.09655560e-12 * tc[4];
        /*species 29: N2 */
        species[29] =
            -1.02090000e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120000e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242750e-13 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.54716300e+04 * invT
            +1.96011760e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -8.35034000e+02 * invT
            +3.34653300e+00
            -2.99142300e+00 * tc[0]
            -3.50032200e-04 * tc[1]
            +9.38971500e-09 * tc[2]
            +7.69298167e-13 * tc[3]
            -7.91376000e-17 * tc[4];
        /*species 2: CH3 */
        species[2] =
            +1.65095130e+04 * invT
            -2.74435930e+00
            -2.97812060e+00 * tc[0]
            -2.89892600e-03 * tc[1]
            +3.29263333e-07 * tc[2]
            -2.56081583e-11 * tc[3]
            +8.95870800e-16 * tc[4];
        /*species 3: O */
        species[3] =
            +2.92308000e+04 * invT
            -3.37824800e+00
            -2.54206000e+00 * tc[0]
            +1.37753100e-05 * tc[1]
            +5.17133833e-10 * tc[2]
            -3.79255583e-13 * tc[3]
            +2.18402600e-17 * tc[4];
        /*species 4: CH4 */
        species[4] =
            -1.00807900e+04 * invT
            -8.93991600e+00
            -1.68347900e+00 * tc[0]
            -5.11862000e-03 * tc[1]
            +6.45854833e-07 * tc[2]
            -5.65465417e-11 * tc[3]
            +2.25171150e-15 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.68362875e+03 * invT
            -3.83691187e+00
            -2.86472886e+00 * tc[0]
            -5.28252240e-04 * tc[1]
            +4.31804597e-08 * tc[2]
            -2.54348895e-12 * tc[3]
            +6.65979380e-17 * tc[4];
        /*species 6: H2O */
        species[6] =
            -2.98992100e+04 * invT
            -5.19067100e+00
            -2.67214600e+00 * tc[0]
            -1.52814650e-03 * tc[1]
            +1.45504333e-07 * tc[2]
            -1.00083000e-11 * tc[3]
            +3.19580900e-16 * tc[4];
        /*species 7: C2H2 */
        species[7] =
            +2.59359992e+04 * invT
            +4.37785085e+00
            -4.14756964e+00 * tc[0]
            -2.98083332e-03 * tc[1]
            +3.95491420e-07 * tc[2]
            -3.89510143e-11 * tc[3]
            +1.80617607e-15 * tc[4];
        /*species 8: CO */
        species[8] =
            -1.42683500e+04 * invT
            -4.08314000e+00
            -3.02507800e+00 * tc[0]
            -7.21344500e-04 * tc[1]
            +9.38471333e-08 * tc[2]
            -8.48817500e-12 * tc[3]
            +3.45547600e-16 * tc[4];
        /*species 9: C2H4 */
        species[9] =
            +4.93988614e+03 * invT
            -9.26925814e+00
            -2.03611116e+00 * tc[0]
            -7.32270755e-03 * tc[1]
            +1.11846319e-06 * tc[2]
            -1.22685769e-10 * tc[3]
            +6.28530305e-15 * tc[4];
        /*species 10: C2H5 */
        species[10] =
            +1.20564550e+04 * invT
            +2.44185557e+00
            -4.28788140e+00 * tc[0]
            -6.21694650e-03 * tc[1]
            +7.35651983e-07 * tc[2]
            -5.88784183e-11 * tc[3]
            +2.10175680e-15 * tc[4];
        /*species 12: C2H6 */
        species[12] =
            -1.27177900e+04 * invT
            +9.06544500e+00
            -4.82593800e+00 * tc[0]
            -6.92021500e-03 * tc[1]
            +7.59543167e-07 * tc[2]
            -5.60413917e-11 * tc[3]
            +1.79908050e-15 * tc[4];
        /*species 13: CH3O */
        species[13] =
            +1.27832500e+02 * invT
            -1.58775000e-01
            -3.77080000e+00 * tc[0]
            -3.93574850e-03 * tc[1]
            +4.42730667e-07 * tc[2]
            -3.28702583e-11 * tc[3]
            +1.05630800e-15 * tc[4];
        /*species 14: O2 */
        species[14] =
            -1.23393000e+03 * invT
            -4.91588000e-01
            -3.69757800e+00 * tc[0]
            -3.06759850e-04 * tc[1]
            +2.09807000e-08 * tc[2]
            -1.47940083e-12 * tc[3]
            +5.68217500e-17 * tc[4];
        /*species 15: HO2 */
        species[15] =
            +1.11856713e+02 * invT
            -7.67891250e-01
            -4.01721090e+00 * tc[0]
            -1.11991006e-03 * tc[1]
            +1.05609692e-07 * tc[2]
            -9.52053083e-12 * tc[3]
            +5.39542675e-16 * tc[4];
        /*species 16: H2O2 */
        species[16] =
            -1.80069600e+04 * invT
            +3.07203000e+00
            -4.57316700e+00 * tc[0]
            -2.16806800e-03 * tc[1]
            +2.45781500e-07 * tc[2]
            -1.95742000e-11 * tc[3]
            +7.15827000e-16 * tc[4];
        /*species 17: CO2 */
        species[17] =
            -4.89669600e+04 * invT
            +4.40901890e+00
            -4.45362300e+00 * tc[0]
            -1.57008450e-03 * tc[1]
            +2.13068500e-07 * tc[2]
            -1.99499750e-11 * tc[3]
            +8.34516500e-16 * tc[4];
        /*species 18: CH3HCO */
        species[18] =
            -2.25931220e+04 * invT
            +7.88490250e+00
            -5.40411080e+00 * tc[0]
            -5.86152950e-03 * tc[1]
            +7.04385617e-07 * tc[2]
            -5.69770425e-11 * tc[3]
            +2.04924315e-15 * tc[4];
        /*species 29: N2 */
        species[29] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988500e-04 * tc[1]
            +9.47460167e-08 * tc[2]
            -8.41420000e-12 * tc[3]
            +3.37667550e-16 * tc[4];
    }

    /*species with midpoint at T=1387 kelvin */
    if (T < 1387) {
        /*species 27: HO2CH2OCHO */
        species[27] =
            -5.80629934e+04 * invT
            -1.27727822e+01
            -3.47935703e+00 * tc[0]
            -2.01476196e-02 * tc[1]
            +5.50182160e-06 * tc[2]
            -1.11966764e-09 * tc[3]
            +1.09300790e-13 * tc[4];
    } else {
        /*species 27: HO2CH2OCHO */
        species[27] =
            -6.23959608e+04 * invT
            +6.93508437e+01
            -1.64584298e+01 * tc[0]
            -4.26341756e-03 * tc[1]
            +5.06855833e-07 * tc[2]
            -4.04664090e-11 * tc[3]
            +1.43658167e-15 * tc[4];
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 26: CH3OCH2O2 */
        species[26] =
            -1.94940940e+04 * invT
            -1.79360640e+01
            -2.21029612e+00 * tc[0]
            -1.84438727e-02 * tc[1]
            +4.70935925e-06 * tc[2]
            -9.64421108e-10 * tc[3]
            +9.85652350e-14 * tc[4];
    } else {
        /*species 26: CH3OCH2O2 */
        species[26] =
            -2.29679238e+04 * invT
            +4.67989874e+01
            -1.24249729e+01 * tc[0]
            -5.93529930e-03 * tc[1]
            +6.79844220e-07 * tc[2]
            -5.29425674e-11 * tc[3]
            +1.84713933e-15 * tc[4];
    }

    /*species with midpoint at T=1200 kelvin */
    if (T < 1200) {
        /*species 11: CH2O */
        species[11] =
            -1.49707930e+04 * invT
            -7.77349870e+00
            -2.69626120e+00 * tc[0]
            -2.46307115e-03 * tc[1]
            -1.38044157e-07 * tc[2]
            +4.58651633e-11 * tc[3]
            +1.98051630e-14 * tc[4];
    } else {
        /*species 11: CH2O */
        species[11] =
            -1.62301730e+04 * invT
            +9.26957180e+00
            -5.14819050e+00 * tc[0]
            -1.43390080e-03 * tc[1]
            +3.96377217e-08 * tc[2]
            +1.34260858e-11 * tc[3]
            -1.42833675e-15 * tc[4];
    }

    /*species with midpoint at T=1362 kelvin */
    if (T < 1362) {
        /*species 21: CH3OCO */
        species[21] =
            -2.14404829e+04 * invT
            -1.37534446e+01
            -3.94199159e+00 * tc[0]
            -1.21717442e-02 * tc[1]
            +2.75992600e-06 * tc[2]
            -3.82114509e-10 * tc[3]
            +1.65897854e-14 * tc[4];
    } else {
        /*species 21: CH3OCO */
        species[21] =
            -2.46616400e+04 * invT
            +4.48791651e+01
            -1.30877600e+01 * tc[0]
            -2.26772475e-03 * tc[1]
            +2.75160607e-07 * tc[2]
            -2.22664398e-11 * tc[3]
            +7.97884315e-16 * tc[4];
    }

    /*species with midpoint at T=1603 kelvin */
    if (T < 1603) {
        /*species 25: HOCH2OCO */
        species[25] =
            -4.39526183e+04 * invT
            +2.54126352e+00
            -6.08180801e+00 * tc[0]
            -6.43841795e-03 * tc[1]
            -3.40699030e-07 * tc[2]
            +5.08462434e-10 * tc[3]
            -8.99102795e-14 * tc[4];
    } else {
        /*species 25: HOCH2OCO */
        species[25] =
            -4.65575743e+04 * invT
            +3.89772656e+01
            -1.13737391e+01 * tc[0]
            -4.08831949e-03 * tc[1]
            +4.86723368e-07 * tc[2]
            -3.88913013e-11 * tc[3]
            +1.38138412e-15 * tc[4];
    }

    /*species with midpoint at T=1686 kelvin */
    if (T < 1686) {
        /*species 22: CH3OCHO */
        species[22] =
            -4.41855167e+04 * invT
            -1.04480741e+01
            -3.08839783e+00 * tc[0]
            -1.01880024e-02 * tc[1]
            +1.14129507e-06 * tc[2]
            +6.06821836e-11 * tc[3]
            -2.81065108e-14 * tc[4];
    } else {
        /*species 22: CH3OCHO */
        species[22] =
            -4.64364769e+04 * invT
            +2.66213830e+01
            -8.69123518e+00 * tc[0]
            -5.77515610e-03 * tc[1]
            +7.12970810e-07 * tc[2]
            -5.85444216e-11 * tc[3]
            +2.12166776e-15 * tc[4];
    }

    /*species with midpoint at T=1402 kelvin */
    if (T < 1402) {
        /*species 28: O2CH2OCH2O2H */
        species[28] =
            -3.27628742e+04 * invT
            -2.34250950e+01
            -1.99640551e+00 * tc[0]
            -2.91613116e-02 * tc[1]
            +9.22099630e-06 * tc[2]
            -2.16508783e-09 * tc[3]
            +2.38570502e-13 * tc[4];
    } else {
        /*species 28: O2CH2OCH2O2H */
        species[28] =
            -3.79207055e+04 * invT
            +8.33885319e+01
            -1.92038046e+01 * tc[0]
            -5.21974205e-03 * tc[1]
            +6.00971565e-07 * tc[2]
            -4.69827369e-11 * tc[3]
            +1.64403607e-15 * tc[4];
    }

    /*species with midpoint at T=2014 kelvin */
    if (T < 2014) {
        /*species 23: CH3OCH2OH */
        species[23] =
            -4.54488899e+04 * invT
            -1.08926047e+01
            -3.15851876e+00 * tc[0]
            -1.22162875e-02 * tc[1]
            +1.44497464e-06 * tc[2]
            +4.94432773e-12 * tc[3]
            -2.18200002e-14 * tc[4];
    } else {
        /*species 23: CH3OCH2OH */
        species[23] =
            -4.76607115e+04 * invT
            +2.57324859e+01
            -8.70981570e+00 * tc[0]
            -7.68011860e-03 * tc[1]
            +9.01672980e-07 * tc[2]
            -7.17144538e-11 * tc[3]
            +2.54409876e-15 * tc[4];
    }
    return;
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cv_R(double * restrict  species, double * restrict  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1376 kelvin */
    if (T < 1376) {
        /*species 19: HCOOH */
        species[19] =
            +4.35481850e-01
            +1.63363016e-02 * tc[1]
            -1.06257421e-05 * tc[2]
            +3.32132977e-09 * tc[3]
            -4.02176103e-13 * tc[4];
    } else {
        /*species 19: HCOOH */
        species[19] =
            +5.68733013e+00
            +5.14289368e-03 * tc[1]
            -1.82238513e-06 * tc[2]
            +2.89719163e-10 * tc[3]
            -1.70892199e-14 * tc[4];
    }

    /*species with midpoint at T=1475 kelvin */
    if (T < 1475) {
        /*species 24: OCH2OCHO */
        species[24] =
            +4.19690837e+00
            +1.58839723e-02 * tc[1]
            +3.53540547e-07 * tc[2]
            -6.10456923e-09 * tc[3]
            +1.94661801e-12 * tc[4];
    } else {
        /*species 24: OCH2OCHO */
        species[24] =
            +1.10233916e+01
            +8.11262659e-03 * tc[1]
            -2.91356462e-06 * tc[2]
            +4.67340384e-10 * tc[3]
            -2.77375525e-14 * tc[4];
    }

    /*species with midpoint at T=710 kelvin */
    if (T < 710) {
        /*species 20: CH3OCH3 */
        species[20] =
            +4.68097447e+00
            -5.39434751e-03 * tc[1]
            +6.49472750e-05 * tc[2]
            -8.05065318e-08 * tc[3]
            +3.27474018e-11 * tc[4];
    } else {
        /*species 20: CH3OCH3 */
        species[20] =
            -1.69184454e-01
            +2.69173263e-02 * tc[1]
            -1.38874777e-05 * tc[2]
            +3.47515079e-09 * tc[3]
            -3.41706784e-13 * tc[4];
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +2.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 2: CH3 */
        species[2] =
            +2.65717970e+00
            +2.12659790e-03 * tc[1]
            +5.45838830e-06 * tc[2]
            -6.61810030e-09 * tc[3]
            +2.46570740e-12 * tc[4];
        /*species 3: O */
        species[3] =
            +1.94642900e+00
            -1.63816600e-03 * tc[1]
            +2.42103200e-06 * tc[2]
            -1.60284300e-09 * tc[3]
            +3.89069600e-13 * tc[4];
        /*species 4: CH4 */
        species[4] =
            -2.21258500e-01
            +1.74766800e-02 * tc[1]
            -2.78340900e-05 * tc[2]
            +3.04970800e-08 * tc[3]
            -1.22393100e-11 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 6: H2O */
        species[6] =
            +2.38684200e+00
            +3.47498200e-03 * tc[1]
            -6.35469600e-06 * tc[2]
            +6.96858100e-09 * tc[3]
            -2.50658800e-12 * tc[4];
        /*species 7: C2H2 */
        species[7] =
            -1.91318906e-01
            +2.33615629e-02 * tc[1]
            -3.55171815e-05 * tc[2]
            +2.80152437e-08 * tc[3]
            -8.50072974e-12 * tc[4];
        /*species 8: CO */
        species[8] =
            +2.26245200e+00
            +1.51194100e-03 * tc[1]
            -3.88175500e-06 * tc[2]
            +5.58194400e-09 * tc[3]
            -2.47495100e-12 * tc[4];
        /*species 9: C2H4 */
        species[9] =
            +2.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 10: C2H5 */
        species[10] =
            +3.30585800e+00
            -4.18336380e-03 * tc[1]
            +4.97072700e-05 * tc[2]
            -5.99058740e-08 * tc[3]
            +2.30484780e-11 * tc[4];
        /*species 12: C2H6 */
        species[12] =
            +4.62539000e-01
            +1.54946700e-02 * tc[1]
            +5.78050700e-06 * tc[2]
            -1.25783200e-08 * tc[3]
            +4.58626700e-12 * tc[4];
        /*species 13: CH3O */
        species[13] =
            +1.10620400e+00
            +7.21659500e-03 * tc[1]
            +5.33847200e-06 * tc[2]
            -7.37763600e-09 * tc[3]
            +2.07561100e-12 * tc[4];
        /*species 14: O2 */
        species[14] =
            +2.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
        /*species 15: HO2 */
        species[15] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 16: H2O2 */
        species[16] =
            +2.38875400e+00
            +6.56922600e-03 * tc[1]
            -1.48501300e-07 * tc[2]
            -4.62580600e-09 * tc[3]
            +2.47151500e-12 * tc[4];
        /*species 17: CO2 */
        species[17] =
            +1.27572500e+00
            +9.92207200e-03 * tc[1]
            -1.04091100e-05 * tc[2]
            +6.86668700e-09 * tc[3]
            -2.11728000e-12 * tc[4];
        /*species 18: CH3HCO */
        species[18] =
            +3.72945950e+00
            -3.19328580e-03 * tc[1]
            +4.75349210e-05 * tc[2]
            -5.74586110e-08 * tc[3]
            +2.19311120e-11 * tc[4];
        /*species 29: N2 */
        species[29] =
            +2.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +1.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 2: CH3 */
        species[2] =
            +1.97812060e+00
            +5.79785200e-03 * tc[1]
            -1.97558000e-06 * tc[2]
            +3.07297900e-10 * tc[3]
            -1.79174160e-14 * tc[4];
        /*species 3: O */
        species[3] =
            +1.54206000e+00
            -2.75506200e-05 * tc[1]
            -3.10280300e-09 * tc[2]
            +4.55106700e-12 * tc[3]
            -4.36805200e-16 * tc[4];
        /*species 4: CH4 */
        species[4] =
            +6.83479000e-01
            +1.02372400e-02 * tc[1]
            -3.87512900e-06 * tc[2]
            +6.78558500e-10 * tc[3]
            -4.50342300e-14 * tc[4];
        /*species 5: OH */
        species[5] =
            +1.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 6: H2O */
        species[6] =
            +1.67214600e+00
            +3.05629300e-03 * tc[1]
            -8.73026000e-07 * tc[2]
            +1.20099600e-10 * tc[3]
            -6.39161800e-15 * tc[4];
        /*species 7: C2H2 */
        species[7] =
            +3.14756964e+00
            +5.96166664e-03 * tc[1]
            -2.37294852e-06 * tc[2]
            +4.67412171e-10 * tc[3]
            -3.61235213e-14 * tc[4];
        /*species 8: CO */
        species[8] =
            +2.02507800e+00
            +1.44268900e-03 * tc[1]
            -5.63082800e-07 * tc[2]
            +1.01858100e-10 * tc[3]
            -6.91095200e-15 * tc[4];
        /*species 9: C2H4 */
        species[9] =
            +1.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 10: C2H5 */
        species[10] =
            +3.28788140e+00
            +1.24338930e-02 * tc[1]
            -4.41391190e-06 * tc[2]
            +7.06541020e-10 * tc[3]
            -4.20351360e-14 * tc[4];
        /*species 12: C2H6 */
        species[12] =
            +3.82593800e+00
            +1.38404300e-02 * tc[1]
            -4.55725900e-06 * tc[2]
            +6.72496700e-10 * tc[3]
            -3.59816100e-14 * tc[4];
        /*species 13: CH3O */
        species[13] =
            +2.77080000e+00
            +7.87149700e-03 * tc[1]
            -2.65638400e-06 * tc[2]
            +3.94443100e-10 * tc[3]
            -2.11261600e-14 * tc[4];
        /*species 14: O2 */
        species[14] =
            +2.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
        /*species 15: HO2 */
        species[15] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 16: H2O2 */
        species[16] =
            +3.57316700e+00
            +4.33613600e-03 * tc[1]
            -1.47468900e-06 * tc[2]
            +2.34890400e-10 * tc[3]
            -1.43165400e-14 * tc[4];
        /*species 17: CO2 */
        species[17] =
            +3.45362300e+00
            +3.14016900e-03 * tc[1]
            -1.27841100e-06 * tc[2]
            +2.39399700e-10 * tc[3]
            -1.66903300e-14 * tc[4];
        /*species 18: CH3HCO */
        species[18] =
            +4.40411080e+00
            +1.17230590e-02 * tc[1]
            -4.22631370e-06 * tc[2]
            +6.83724510e-10 * tc[3]
            -4.09848630e-14 * tc[4];
        /*species 29: N2 */
        species[29] =
            +1.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }

    /*species with midpoint at T=1387 kelvin */
    if (T < 1387) {
        /*species 27: HO2CH2OCHO */
        species[27] =
            +2.47935703e+00
            +4.02952392e-02 * tc[1]
            -3.30109296e-05 * tc[2]
            +1.34360117e-08 * tc[3]
            -2.18601580e-12 * tc[4];
    } else {
        /*species 27: HO2CH2OCHO */
        species[27] =
            +1.54584298e+01
            +8.52683511e-03 * tc[1]
            -3.04113500e-06 * tc[2]
            +4.85596908e-10 * tc[3]
            -2.87316334e-14 * tc[4];
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 26: CH3OCH2O2 */
        species[26] =
            +1.21029612e+00
            +3.68877454e-02 * tc[1]
            -2.82561555e-05 * tc[2]
            +1.15730533e-08 * tc[3]
            -1.97130470e-12 * tc[4];
    } else {
        /*species 26: CH3OCH2O2 */
        species[26] =
            +1.14249729e+01
            +1.18705986e-02 * tc[1]
            -4.07906532e-06 * tc[2]
            +6.35310809e-10 * tc[3]
            -3.69427867e-14 * tc[4];
    }

    /*species with midpoint at T=1200 kelvin */
    if (T < 1200) {
        /*species 11: CH2O */
        species[11] =
            +1.69626120e+00
            +4.92614230e-03 * tc[1]
            +8.28264940e-07 * tc[2]
            -5.50381960e-10 * tc[3]
            -3.96103260e-13 * tc[4];
    } else {
        /*species 11: CH2O */
        species[11] =
            +4.14819050e+00
            +2.86780160e-03 * tc[1]
            -2.37826330e-07 * tc[2]
            -1.61113030e-10 * tc[3]
            +2.85667350e-14 * tc[4];
    }

    /*species with midpoint at T=1362 kelvin */
    if (T < 1362) {
        /*species 21: CH3OCO */
        species[21] =
            +2.94199159e+00
            +2.43434884e-02 * tc[1]
            -1.65595560e-05 * tc[2]
            +4.58537411e-09 * tc[3]
            -3.31795708e-13 * tc[4];
    } else {
        /*species 21: CH3OCO */
        species[21] =
            +1.20877600e+01
            +4.53544950e-03 * tc[1]
            -1.65096364e-06 * tc[2]
            +2.67197277e-10 * tc[3]
            -1.59576863e-14 * tc[4];
    }

    /*species with midpoint at T=1603 kelvin */
    if (T < 1603) {
        /*species 25: HOCH2OCO */
        species[25] =
            +5.08180801e+00
            +1.28768359e-02 * tc[1]
            +2.04419418e-06 * tc[2]
            -6.10154921e-09 * tc[3]
            +1.79820559e-12 * tc[4];
    } else {
        /*species 25: HOCH2OCO */
        species[25] =
            +1.03737391e+01
            +8.17663898e-03 * tc[1]
            -2.92034021e-06 * tc[2]
            +4.66695616e-10 * tc[3]
            -2.76276823e-14 * tc[4];
    }

    /*species with midpoint at T=1686 kelvin */
    if (T < 1686) {
        /*species 22: CH3OCHO */
        species[22] =
            +2.08839783e+00
            +2.03760048e-02 * tc[1]
            -6.84777040e-06 * tc[2]
            -7.28186203e-10 * tc[3]
            +5.62130216e-13 * tc[4];
    } else {
        /*species 22: CH3OCHO */
        species[22] =
            +7.69123518e+00
            +1.15503122e-02 * tc[1]
            -4.27782486e-06 * tc[2]
            +7.02533059e-10 * tc[3]
            -4.24333552e-14 * tc[4];
    }

    /*species with midpoint at T=1402 kelvin */
    if (T < 1402) {
        /*species 28: O2CH2OCH2O2H */
        species[28] =
            +9.96405510e-01
            +5.83226232e-02 * tc[1]
            -5.53259778e-05 * tc[2]
            +2.59810540e-08 * tc[3]
            -4.77141005e-12 * tc[4];
    } else {
        /*species 28: O2CH2OCH2O2H */
        species[28] =
            +1.82038046e+01
            +1.04394841e-02 * tc[1]
            -3.60582939e-06 * tc[2]
            +5.63792843e-10 * tc[3]
            -3.28807214e-14 * tc[4];
    }

    /*species with midpoint at T=2014 kelvin */
    if (T < 2014) {
        /*species 23: CH3OCH2OH */
        species[23] =
            +2.15851876e+00
            +2.44325751e-02 * tc[1]
            -8.66984784e-06 * tc[2]
            -5.93319328e-11 * tc[3]
            +4.36400003e-13 * tc[4];
    } else {
        /*species 23: CH3OCH2OH */
        species[23] =
            +7.70981570e+00
            +1.53602372e-02 * tc[1]
            -5.41003788e-06 * tc[2]
            +8.60573446e-10 * tc[3]
            -5.08819752e-14 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cp_R(double * restrict  species, double * restrict  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1376 kelvin */
    if (T < 1376) {
        /*species 19: HCOOH */
        species[19] =
            +1.43548185e+00
            +1.63363016e-02 * tc[1]
            -1.06257421e-05 * tc[2]
            +3.32132977e-09 * tc[3]
            -4.02176103e-13 * tc[4];
    } else {
        /*species 19: HCOOH */
        species[19] =
            +6.68733013e+00
            +5.14289368e-03 * tc[1]
            -1.82238513e-06 * tc[2]
            +2.89719163e-10 * tc[3]
            -1.70892199e-14 * tc[4];
    }

    /*species with midpoint at T=1475 kelvin */
    if (T < 1475) {
        /*species 24: OCH2OCHO */
        species[24] =
            +5.19690837e+00
            +1.58839723e-02 * tc[1]
            +3.53540547e-07 * tc[2]
            -6.10456923e-09 * tc[3]
            +1.94661801e-12 * tc[4];
    } else {
        /*species 24: OCH2OCHO */
        species[24] =
            +1.20233916e+01
            +8.11262659e-03 * tc[1]
            -2.91356462e-06 * tc[2]
            +4.67340384e-10 * tc[3]
            -2.77375525e-14 * tc[4];
    }

    /*species with midpoint at T=710 kelvin */
    if (T < 710) {
        /*species 20: CH3OCH3 */
        species[20] =
            +5.68097447e+00
            -5.39434751e-03 * tc[1]
            +6.49472750e-05 * tc[2]
            -8.05065318e-08 * tc[3]
            +3.27474018e-11 * tc[4];
    } else {
        /*species 20: CH3OCH3 */
        species[20] =
            +8.30815546e-01
            +2.69173263e-02 * tc[1]
            -1.38874777e-05 * tc[2]
            +3.47515079e-09 * tc[3]
            -3.41706784e-13 * tc[4];
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +3.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 2: CH3 */
        species[2] =
            +3.65717970e+00
            +2.12659790e-03 * tc[1]
            +5.45838830e-06 * tc[2]
            -6.61810030e-09 * tc[3]
            +2.46570740e-12 * tc[4];
        /*species 3: O */
        species[3] =
            +2.94642900e+00
            -1.63816600e-03 * tc[1]
            +2.42103200e-06 * tc[2]
            -1.60284300e-09 * tc[3]
            +3.89069600e-13 * tc[4];
        /*species 4: CH4 */
        species[4] =
            +7.78741500e-01
            +1.74766800e-02 * tc[1]
            -2.78340900e-05 * tc[2]
            +3.04970800e-08 * tc[3]
            -1.22393100e-11 * tc[4];
        /*species 5: OH */
        species[5] =
            +4.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 6: H2O */
        species[6] =
            +3.38684200e+00
            +3.47498200e-03 * tc[1]
            -6.35469600e-06 * tc[2]
            +6.96858100e-09 * tc[3]
            -2.50658800e-12 * tc[4];
        /*species 7: C2H2 */
        species[7] =
            +8.08681094e-01
            +2.33615629e-02 * tc[1]
            -3.55171815e-05 * tc[2]
            +2.80152437e-08 * tc[3]
            -8.50072974e-12 * tc[4];
        /*species 8: CO */
        species[8] =
            +3.26245200e+00
            +1.51194100e-03 * tc[1]
            -3.88175500e-06 * tc[2]
            +5.58194400e-09 * tc[3]
            -2.47495100e-12 * tc[4];
        /*species 9: C2H4 */
        species[9] =
            +3.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 10: C2H5 */
        species[10] =
            +4.30585800e+00
            -4.18336380e-03 * tc[1]
            +4.97072700e-05 * tc[2]
            -5.99058740e-08 * tc[3]
            +2.30484780e-11 * tc[4];
        /*species 12: C2H6 */
        species[12] =
            +1.46253900e+00
            +1.54946700e-02 * tc[1]
            +5.78050700e-06 * tc[2]
            -1.25783200e-08 * tc[3]
            +4.58626700e-12 * tc[4];
        /*species 13: CH3O */
        species[13] =
            +2.10620400e+00
            +7.21659500e-03 * tc[1]
            +5.33847200e-06 * tc[2]
            -7.37763600e-09 * tc[3]
            +2.07561100e-12 * tc[4];
        /*species 14: O2 */
        species[14] =
            +3.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
        /*species 15: HO2 */
        species[15] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 16: H2O2 */
        species[16] =
            +3.38875400e+00
            +6.56922600e-03 * tc[1]
            -1.48501300e-07 * tc[2]
            -4.62580600e-09 * tc[3]
            +2.47151500e-12 * tc[4];
        /*species 17: CO2 */
        species[17] =
            +2.27572500e+00
            +9.92207200e-03 * tc[1]
            -1.04091100e-05 * tc[2]
            +6.86668700e-09 * tc[3]
            -2.11728000e-12 * tc[4];
        /*species 18: CH3HCO */
        species[18] =
            +4.72945950e+00
            -3.19328580e-03 * tc[1]
            +4.75349210e-05 * tc[2]
            -5.74586110e-08 * tc[3]
            +2.19311120e-11 * tc[4];
        /*species 29: N2 */
        species[29] =
            +3.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +2.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 2: CH3 */
        species[2] =
            +2.97812060e+00
            +5.79785200e-03 * tc[1]
            -1.97558000e-06 * tc[2]
            +3.07297900e-10 * tc[3]
            -1.79174160e-14 * tc[4];
        /*species 3: O */
        species[3] =
            +2.54206000e+00
            -2.75506200e-05 * tc[1]
            -3.10280300e-09 * tc[2]
            +4.55106700e-12 * tc[3]
            -4.36805200e-16 * tc[4];
        /*species 4: CH4 */
        species[4] =
            +1.68347900e+00
            +1.02372400e-02 * tc[1]
            -3.87512900e-06 * tc[2]
            +6.78558500e-10 * tc[3]
            -4.50342300e-14 * tc[4];
        /*species 5: OH */
        species[5] =
            +2.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 6: H2O */
        species[6] =
            +2.67214600e+00
            +3.05629300e-03 * tc[1]
            -8.73026000e-07 * tc[2]
            +1.20099600e-10 * tc[3]
            -6.39161800e-15 * tc[4];
        /*species 7: C2H2 */
        species[7] =
            +4.14756964e+00
            +5.96166664e-03 * tc[1]
            -2.37294852e-06 * tc[2]
            +4.67412171e-10 * tc[3]
            -3.61235213e-14 * tc[4];
        /*species 8: CO */
        species[8] =
            +3.02507800e+00
            +1.44268900e-03 * tc[1]
            -5.63082800e-07 * tc[2]
            +1.01858100e-10 * tc[3]
            -6.91095200e-15 * tc[4];
        /*species 9: C2H4 */
        species[9] =
            +2.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 10: C2H5 */
        species[10] =
            +4.28788140e+00
            +1.24338930e-02 * tc[1]
            -4.41391190e-06 * tc[2]
            +7.06541020e-10 * tc[3]
            -4.20351360e-14 * tc[4];
        /*species 12: C2H6 */
        species[12] =
            +4.82593800e+00
            +1.38404300e-02 * tc[1]
            -4.55725900e-06 * tc[2]
            +6.72496700e-10 * tc[3]
            -3.59816100e-14 * tc[4];
        /*species 13: CH3O */
        species[13] =
            +3.77080000e+00
            +7.87149700e-03 * tc[1]
            -2.65638400e-06 * tc[2]
            +3.94443100e-10 * tc[3]
            -2.11261600e-14 * tc[4];
        /*species 14: O2 */
        species[14] =
            +3.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
        /*species 15: HO2 */
        species[15] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 16: H2O2 */
        species[16] =
            +4.57316700e+00
            +4.33613600e-03 * tc[1]
            -1.47468900e-06 * tc[2]
            +2.34890400e-10 * tc[3]
            -1.43165400e-14 * tc[4];
        /*species 17: CO2 */
        species[17] =
            +4.45362300e+00
            +3.14016900e-03 * tc[1]
            -1.27841100e-06 * tc[2]
            +2.39399700e-10 * tc[3]
            -1.66903300e-14 * tc[4];
        /*species 18: CH3HCO */
        species[18] =
            +5.40411080e+00
            +1.17230590e-02 * tc[1]
            -4.22631370e-06 * tc[2]
            +6.83724510e-10 * tc[3]
            -4.09848630e-14 * tc[4];
        /*species 29: N2 */
        species[29] =
            +2.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }

    /*species with midpoint at T=1387 kelvin */
    if (T < 1387) {
        /*species 27: HO2CH2OCHO */
        species[27] =
            +3.47935703e+00
            +4.02952392e-02 * tc[1]
            -3.30109296e-05 * tc[2]
            +1.34360117e-08 * tc[3]
            -2.18601580e-12 * tc[4];
    } else {
        /*species 27: HO2CH2OCHO */
        species[27] =
            +1.64584298e+01
            +8.52683511e-03 * tc[1]
            -3.04113500e-06 * tc[2]
            +4.85596908e-10 * tc[3]
            -2.87316334e-14 * tc[4];
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 26: CH3OCH2O2 */
        species[26] =
            +2.21029612e+00
            +3.68877454e-02 * tc[1]
            -2.82561555e-05 * tc[2]
            +1.15730533e-08 * tc[3]
            -1.97130470e-12 * tc[4];
    } else {
        /*species 26: CH3OCH2O2 */
        species[26] =
            +1.24249729e+01
            +1.18705986e-02 * tc[1]
            -4.07906532e-06 * tc[2]
            +6.35310809e-10 * tc[3]
            -3.69427867e-14 * tc[4];
    }

    /*species with midpoint at T=1200 kelvin */
    if (T < 1200) {
        /*species 11: CH2O */
        species[11] =
            +2.69626120e+00
            +4.92614230e-03 * tc[1]
            +8.28264940e-07 * tc[2]
            -5.50381960e-10 * tc[3]
            -3.96103260e-13 * tc[4];
    } else {
        /*species 11: CH2O */
        species[11] =
            +5.14819050e+00
            +2.86780160e-03 * tc[1]
            -2.37826330e-07 * tc[2]
            -1.61113030e-10 * tc[3]
            +2.85667350e-14 * tc[4];
    }

    /*species with midpoint at T=1362 kelvin */
    if (T < 1362) {
        /*species 21: CH3OCO */
        species[21] =
            +3.94199159e+00
            +2.43434884e-02 * tc[1]
            -1.65595560e-05 * tc[2]
            +4.58537411e-09 * tc[3]
            -3.31795708e-13 * tc[4];
    } else {
        /*species 21: CH3OCO */
        species[21] =
            +1.30877600e+01
            +4.53544950e-03 * tc[1]
            -1.65096364e-06 * tc[2]
            +2.67197277e-10 * tc[3]
            -1.59576863e-14 * tc[4];
    }

    /*species with midpoint at T=1603 kelvin */
    if (T < 1603) {
        /*species 25: HOCH2OCO */
        species[25] =
            +6.08180801e+00
            +1.28768359e-02 * tc[1]
            +2.04419418e-06 * tc[2]
            -6.10154921e-09 * tc[3]
            +1.79820559e-12 * tc[4];
    } else {
        /*species 25: HOCH2OCO */
        species[25] =
            +1.13737391e+01
            +8.17663898e-03 * tc[1]
            -2.92034021e-06 * tc[2]
            +4.66695616e-10 * tc[3]
            -2.76276823e-14 * tc[4];
    }

    /*species with midpoint at T=1686 kelvin */
    if (T < 1686) {
        /*species 22: CH3OCHO */
        species[22] =
            +3.08839783e+00
            +2.03760048e-02 * tc[1]
            -6.84777040e-06 * tc[2]
            -7.28186203e-10 * tc[3]
            +5.62130216e-13 * tc[4];
    } else {
        /*species 22: CH3OCHO */
        species[22] =
            +8.69123518e+00
            +1.15503122e-02 * tc[1]
            -4.27782486e-06 * tc[2]
            +7.02533059e-10 * tc[3]
            -4.24333552e-14 * tc[4];
    }

    /*species with midpoint at T=1402 kelvin */
    if (T < 1402) {
        /*species 28: O2CH2OCH2O2H */
        species[28] =
            +1.99640551e+00
            +5.83226232e-02 * tc[1]
            -5.53259778e-05 * tc[2]
            +2.59810540e-08 * tc[3]
            -4.77141005e-12 * tc[4];
    } else {
        /*species 28: O2CH2OCH2O2H */
        species[28] =
            +1.92038046e+01
            +1.04394841e-02 * tc[1]
            -3.60582939e-06 * tc[2]
            +5.63792843e-10 * tc[3]
            -3.28807214e-14 * tc[4];
    }

    /*species with midpoint at T=2014 kelvin */
    if (T < 2014) {
        /*species 23: CH3OCH2OH */
        species[23] =
            +3.15851876e+00
            +2.44325751e-02 * tc[1]
            -8.66984784e-06 * tc[2]
            -5.93319328e-11 * tc[3]
            +4.36400003e-13 * tc[4];
    } else {
        /*species 23: CH3OCH2OH */
        species[23] =
            +8.70981570e+00
            +1.53602372e-02 * tc[1]
            -5.41003788e-06 * tc[2]
            +8.60573446e-10 * tc[3]
            -5.08819752e-14 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesInternalEnergy(double * restrict  species, double * restrict  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1376 kelvin */
    if (T < 1376) {
        /*species 19: HCOOH */
        species[19] =
            +4.35481850e-01
            +8.16815080e-03 * tc[1]
            -3.54191403e-06 * tc[2]
            +8.30332443e-10 * tc[3]
            -8.04352206e-14 * tc[4]
            -4.64616504e+04 * invT;
    } else {
        /*species 19: HCOOH */
        species[19] =
            +5.68733013e+00
            +2.57144684e-03 * tc[1]
            -6.07461710e-07 * tc[2]
            +7.24297908e-11 * tc[3]
            -3.41784398e-15 * tc[4]
            -4.83995400e+04 * invT;
    }

    /*species with midpoint at T=1475 kelvin */
    if (T < 1475) {
        /*species 24: OCH2OCHO */
        species[24] =
            +4.19690837e+00
            +7.94198615e-03 * tc[1]
            +1.17846849e-07 * tc[2]
            -1.52614231e-09 * tc[3]
            +3.89323602e-13 * tc[4]
            -4.02242792e+04 * invT;
    } else {
        /*species 24: OCH2OCHO */
        species[24] =
            +1.10233916e+01
            +4.05631329e-03 * tc[1]
            -9.71188207e-07 * tc[2]
            +1.16835096e-10 * tc[3]
            -5.54751050e-15 * tc[4]
            -4.33647231e+04 * invT;
    }

    /*species with midpoint at T=710 kelvin */
    if (T < 710) {
        /*species 20: CH3OCH3 */
        species[20] =
            +4.68097447e+00
            -2.69717376e-03 * tc[1]
            +2.16490917e-05 * tc[2]
            -2.01266330e-08 * tc[3]
            +6.54948036e-12 * tc[4]
            -2.39755455e+04 * invT;
    } else {
        /*species 20: CH3OCH3 */
        species[20] =
            -1.69184454e-01
            +1.34586631e-02 * tc[1]
            -4.62915923e-06 * tc[2]
            +8.68787697e-10 * tc[3]
            -6.83413568e-14 * tc[4]
            -2.34120975e+04 * invT;
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +2.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 * invT;
        /*species 2: CH3 */
        species[2] =
            +2.65717970e+00
            +1.06329895e-03 * tc[1]
            +1.81946277e-06 * tc[2]
            -1.65452507e-09 * tc[3]
            +4.93141480e-13 * tc[4]
            +1.64227160e+04 * invT;
        /*species 3: O */
        species[3] =
            +1.94642900e+00
            -8.19083000e-04 * tc[1]
            +8.07010667e-07 * tc[2]
            -4.00710750e-10 * tc[3]
            +7.78139200e-14 * tc[4]
            +2.91476400e+04 * invT;
        /*species 4: CH4 */
        species[4] =
            -2.21258500e-01
            +8.73834000e-03 * tc[1]
            -9.27803000e-06 * tc[2]
            +7.62427000e-09 * tc[3]
            -2.44786200e-12 * tc[4]
            -9.82522900e+03 * invT;
        /*species 5: OH */
        species[5] =
            +3.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.34630913e+03 * invT;
        /*species 6: H2O */
        species[6] =
            +2.38684200e+00
            +1.73749100e-03 * tc[1]
            -2.11823200e-06 * tc[2]
            +1.74214525e-09 * tc[3]
            -5.01317600e-13 * tc[4]
            -3.02081100e+04 * invT;
        /*species 7: C2H2 */
        species[7] =
            -1.91318906e-01
            +1.16807815e-02 * tc[1]
            -1.18390605e-05 * tc[2]
            +7.00381092e-09 * tc[3]
            -1.70014595e-12 * tc[4]
            +2.64289807e+04 * invT;
        /*species 8: CO */
        species[8] =
            +2.26245200e+00
            +7.55970500e-04 * tc[1]
            -1.29391833e-06 * tc[2]
            +1.39548600e-09 * tc[3]
            -4.94990200e-13 * tc[4]
            -1.43105400e+04 * invT;
        /*species 9: C2H4 */
        species[9] =
            +2.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 * invT;
        /*species 10: C2H5 */
        species[10] =
            +3.30585800e+00
            -2.09168190e-03 * tc[1]
            +1.65690900e-05 * tc[2]
            -1.49764685e-08 * tc[3]
            +4.60969560e-12 * tc[4]
            +1.28417140e+04 * invT;
        /*species 12: C2H6 */
        species[12] =
            +4.62539000e-01
            +7.74733500e-03 * tc[1]
            +1.92683567e-06 * tc[2]
            -3.14458000e-09 * tc[3]
            +9.17253400e-13 * tc[4]
            -1.12391800e+04 * invT;
        /*species 13: CH3O */
        species[13] =
            +1.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122200e-13 * tc[4]
            +9.78601100e+02 * invT;
        /*species 14: O2 */
        species[14] =
            +2.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 * invT;
        /*species 15: HO2 */
        species[15] =
            +3.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 16: H2O2 */
        species[16] =
            +2.38875400e+00
            +3.28461300e-03 * tc[1]
            -4.95004333e-08 * tc[2]
            -1.15645150e-09 * tc[3]
            +4.94303000e-13 * tc[4]
            -1.76631500e+04 * invT;
        /*species 17: CO2 */
        species[17] =
            +1.27572500e+00
            +4.96103600e-03 * tc[1]
            -3.46970333e-06 * tc[2]
            +1.71667175e-09 * tc[3]
            -4.23456000e-13 * tc[4]
            -4.83731400e+04 * invT;
        /*species 18: CH3HCO */
        species[18] =
            +3.72945950e+00
            -1.59664290e-03 * tc[1]
            +1.58449737e-05 * tc[2]
            -1.43646527e-08 * tc[3]
            +4.38622240e-12 * tc[4]
            -2.15728780e+04 * invT;
        /*species 29: N2 */
        species[29] =
            +2.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
    } else {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +1.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 * invT;
        /*species 2: CH3 */
        species[2] =
            +1.97812060e+00
            +2.89892600e-03 * tc[1]
            -6.58526667e-07 * tc[2]
            +7.68244750e-11 * tc[3]
            -3.58348320e-15 * tc[4]
            +1.65095130e+04 * invT;
        /*species 3: O */
        species[3] =
            +1.54206000e+00
            -1.37753100e-05 * tc[1]
            -1.03426767e-09 * tc[2]
            +1.13776675e-12 * tc[3]
            -8.73610400e-17 * tc[4]
            +2.92308000e+04 * invT;
        /*species 4: CH4 */
        species[4] =
            +6.83479000e-01
            +5.11862000e-03 * tc[1]
            -1.29170967e-06 * tc[2]
            +1.69639625e-10 * tc[3]
            -9.00684600e-15 * tc[4]
            -1.00807900e+04 * invT;
        /*species 5: OH */
        species[5] =
            +1.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.68362875e+03 * invT;
        /*species 6: H2O */
        species[6] =
            +1.67214600e+00
            +1.52814650e-03 * tc[1]
            -2.91008667e-07 * tc[2]
            +3.00249000e-11 * tc[3]
            -1.27832360e-15 * tc[4]
            -2.98992100e+04 * invT;
        /*species 7: C2H2 */
        species[7] =
            +3.14756964e+00
            +2.98083332e-03 * tc[1]
            -7.90982840e-07 * tc[2]
            +1.16853043e-10 * tc[3]
            -7.22470426e-15 * tc[4]
            +2.59359992e+04 * invT;
        /*species 8: CO */
        species[8] =
            +2.02507800e+00
            +7.21344500e-04 * tc[1]
            -1.87694267e-07 * tc[2]
            +2.54645250e-11 * tc[3]
            -1.38219040e-15 * tc[4]
            -1.42683500e+04 * invT;
        /*species 9: C2H4 */
        species[9] =
            +1.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 * invT;
        /*species 10: C2H5 */
        species[10] =
            +3.28788140e+00
            +6.21694650e-03 * tc[1]
            -1.47130397e-06 * tc[2]
            +1.76635255e-10 * tc[3]
            -8.40702720e-15 * tc[4]
            +1.20564550e+04 * invT;
        /*species 12: C2H6 */
        species[12] =
            +3.82593800e+00
            +6.92021500e-03 * tc[1]
            -1.51908633e-06 * tc[2]
            +1.68124175e-10 * tc[3]
            -7.19632200e-15 * tc[4]
            -1.27177900e+04 * invT;
        /*species 13: CH3O */
        species[13] =
            +2.77080000e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832500e+02 * invT;
        /*species 14: O2 */
        species[14] =
            +2.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 * invT;
        /*species 15: HO2 */
        species[15] =
            +3.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 16: H2O2 */
        species[16] =
            +3.57316700e+00
            +2.16806800e-03 * tc[1]
            -4.91563000e-07 * tc[2]
            +5.87226000e-11 * tc[3]
            -2.86330800e-15 * tc[4]
            -1.80069600e+04 * invT;
        /*species 17: CO2 */
        species[17] =
            +3.45362300e+00
            +1.57008450e-03 * tc[1]
            -4.26137000e-07 * tc[2]
            +5.98499250e-11 * tc[3]
            -3.33806600e-15 * tc[4]
            -4.89669600e+04 * invT;
        /*species 18: CH3HCO */
        species[18] =
            +4.40411080e+00
            +5.86152950e-03 * tc[1]
            -1.40877123e-06 * tc[2]
            +1.70931128e-10 * tc[3]
            -8.19697260e-15 * tc[4]
            -2.25931220e+04 * invT;
        /*species 29: N2 */
        species[29] =
            +1.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }

    /*species with midpoint at T=1387 kelvin */
    if (T < 1387) {
        /*species 27: HO2CH2OCHO */
        species[27] =
            +2.47935703e+00
            +2.01476196e-02 * tc[1]
            -1.10036432e-05 * tc[2]
            +3.35900293e-09 * tc[3]
            -4.37203160e-13 * tc[4]
            -5.80629934e+04 * invT;
    } else {
        /*species 27: HO2CH2OCHO */
        species[27] =
            +1.54584298e+01
            +4.26341756e-03 * tc[1]
            -1.01371167e-06 * tc[2]
            +1.21399227e-10 * tc[3]
            -5.74632668e-15 * tc[4]
            -6.23959608e+04 * invT;
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 26: CH3OCH2O2 */
        species[26] =
            +1.21029612e+00
            +1.84438727e-02 * tc[1]
            -9.41871850e-06 * tc[2]
            +2.89326332e-09 * tc[3]
            -3.94260940e-13 * tc[4]
            -1.94940940e+04 * invT;
    } else {
        /*species 26: CH3OCH2O2 */
        species[26] =
            +1.14249729e+01
            +5.93529930e-03 * tc[1]
            -1.35968844e-06 * tc[2]
            +1.58827702e-10 * tc[3]
            -7.38855734e-15 * tc[4]
            -2.29679238e+04 * invT;
    }

    /*species with midpoint at T=1200 kelvin */
    if (T < 1200) {
        /*species 11: CH2O */
        species[11] =
            +1.69626120e+00
            +2.46307115e-03 * tc[1]
            +2.76088313e-07 * tc[2]
            -1.37595490e-10 * tc[3]
            -7.92206520e-14 * tc[4]
            -1.49707930e+04 * invT;
    } else {
        /*species 11: CH2O */
        species[11] =
            +4.14819050e+00
            +1.43390080e-03 * tc[1]
            -7.92754433e-08 * tc[2]
            -4.02782575e-11 * tc[3]
            +5.71334700e-15 * tc[4]
            -1.62301730e+04 * invT;
    }

    /*species with midpoint at T=1362 kelvin */
    if (T < 1362) {
        /*species 21: CH3OCO */
        species[21] =
            +2.94199159e+00
            +1.21717442e-02 * tc[1]
            -5.51985200e-06 * tc[2]
            +1.14634353e-09 * tc[3]
            -6.63591416e-14 * tc[4]
            -2.14404829e+04 * invT;
    } else {
        /*species 21: CH3OCO */
        species[21] =
            +1.20877600e+01
            +2.26772475e-03 * tc[1]
            -5.50321213e-07 * tc[2]
            +6.67993193e-11 * tc[3]
            -3.19153726e-15 * tc[4]
            -2.46616400e+04 * invT;
    }

    /*species with midpoint at T=1603 kelvin */
    if (T < 1603) {
        /*species 25: HOCH2OCO */
        species[25] =
            +5.08180801e+00
            +6.43841795e-03 * tc[1]
            +6.81398060e-07 * tc[2]
            -1.52538730e-09 * tc[3]
            +3.59641118e-13 * tc[4]
            -4.39526183e+04 * invT;
    } else {
        /*species 25: HOCH2OCO */
        species[25] =
            +1.03737391e+01
            +4.08831949e-03 * tc[1]
            -9.73446737e-07 * tc[2]
            +1.16673904e-10 * tc[3]
            -5.52553646e-15 * tc[4]
            -4.65575743e+04 * invT;
    }

    /*species with midpoint at T=1686 kelvin */
    if (T < 1686) {
        /*species 22: CH3OCHO */
        species[22] =
            +2.08839783e+00
            +1.01880024e-02 * tc[1]
            -2.28259013e-06 * tc[2]
            -1.82046551e-10 * tc[3]
            +1.12426043e-13 * tc[4]
            -4.41855167e+04 * invT;
    } else {
        /*species 22: CH3OCHO */
        species[22] =
            +7.69123518e+00
            +5.77515610e-03 * tc[1]
            -1.42594162e-06 * tc[2]
            +1.75633265e-10 * tc[3]
            -8.48667104e-15 * tc[4]
            -4.64364769e+04 * invT;
    }

    /*species with midpoint at T=1402 kelvin */
    if (T < 1402) {
        /*species 28: O2CH2OCH2O2H */
        species[28] =
            +9.96405510e-01
            +2.91613116e-02 * tc[1]
            -1.84419926e-05 * tc[2]
            +6.49526350e-09 * tc[3]
            -9.54282010e-13 * tc[4]
            -3.27628742e+04 * invT;
    } else {
        /*species 28: O2CH2OCH2O2H */
        species[28] =
            +1.82038046e+01
            +5.21974205e-03 * tc[1]
            -1.20194313e-06 * tc[2]
            +1.40948211e-10 * tc[3]
            -6.57614428e-15 * tc[4]
            -3.79207055e+04 * invT;
    }

    /*species with midpoint at T=2014 kelvin */
    if (T < 2014) {
        /*species 23: CH3OCH2OH */
        species[23] =
            +2.15851876e+00
            +1.22162875e-02 * tc[1]
            -2.88994928e-06 * tc[2]
            -1.48329832e-11 * tc[3]
            +8.72800006e-14 * tc[4]
            -4.54488899e+04 * invT;
    } else {
        /*species 23: CH3OCH2OH */
        species[23] =
            +7.70981570e+00
            +7.68011860e-03 * tc[1]
            -1.80334596e-06 * tc[2]
            +2.15143362e-10 * tc[3]
            -1.01763950e-14 * tc[4]
            -4.76607115e+04 * invT;
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEnthalpy(double * restrict  species, double * restrict  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1376 kelvin */
    if (T < 1376) {
        /*species 19: HCOOH */
        species[19] =
            +1.43548185e+00
            +8.16815080e-03 * tc[1]
            -3.54191403e-06 * tc[2]
            +8.30332443e-10 * tc[3]
            -8.04352206e-14 * tc[4]
            -4.64616504e+04 * invT;
    } else {
        /*species 19: HCOOH */
        species[19] =
            +6.68733013e+00
            +2.57144684e-03 * tc[1]
            -6.07461710e-07 * tc[2]
            +7.24297908e-11 * tc[3]
            -3.41784398e-15 * tc[4]
            -4.83995400e+04 * invT;
    }

    /*species with midpoint at T=1475 kelvin */
    if (T < 1475) {
        /*species 24: OCH2OCHO */
        species[24] =
            +5.19690837e+00
            +7.94198615e-03 * tc[1]
            +1.17846849e-07 * tc[2]
            -1.52614231e-09 * tc[3]
            +3.89323602e-13 * tc[4]
            -4.02242792e+04 * invT;
    } else {
        /*species 24: OCH2OCHO */
        species[24] =
            +1.20233916e+01
            +4.05631329e-03 * tc[1]
            -9.71188207e-07 * tc[2]
            +1.16835096e-10 * tc[3]
            -5.54751050e-15 * tc[4]
            -4.33647231e+04 * invT;
    }

    /*species with midpoint at T=710 kelvin */
    if (T < 710) {
        /*species 20: CH3OCH3 */
        species[20] =
            +5.68097447e+00
            -2.69717376e-03 * tc[1]
            +2.16490917e-05 * tc[2]
            -2.01266330e-08 * tc[3]
            +6.54948036e-12 * tc[4]
            -2.39755455e+04 * invT;
    } else {
        /*species 20: CH3OCH3 */
        species[20] =
            +8.30815546e-01
            +1.34586631e-02 * tc[1]
            -4.62915923e-06 * tc[2]
            +8.68787697e-10 * tc[3]
            -6.83413568e-14 * tc[4]
            -2.34120975e+04 * invT;
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +3.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 * invT;
        /*species 2: CH3 */
        species[2] =
            +3.65717970e+00
            +1.06329895e-03 * tc[1]
            +1.81946277e-06 * tc[2]
            -1.65452507e-09 * tc[3]
            +4.93141480e-13 * tc[4]
            +1.64227160e+04 * invT;
        /*species 3: O */
        species[3] =
            +2.94642900e+00
            -8.19083000e-04 * tc[1]
            +8.07010667e-07 * tc[2]
            -4.00710750e-10 * tc[3]
            +7.78139200e-14 * tc[4]
            +2.91476400e+04 * invT;
        /*species 4: CH4 */
        species[4] =
            +7.78741500e-01
            +8.73834000e-03 * tc[1]
            -9.27803000e-06 * tc[2]
            +7.62427000e-09 * tc[3]
            -2.44786200e-12 * tc[4]
            -9.82522900e+03 * invT;
        /*species 5: OH */
        species[5] =
            +4.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.34630913e+03 * invT;
        /*species 6: H2O */
        species[6] =
            +3.38684200e+00
            +1.73749100e-03 * tc[1]
            -2.11823200e-06 * tc[2]
            +1.74214525e-09 * tc[3]
            -5.01317600e-13 * tc[4]
            -3.02081100e+04 * invT;
        /*species 7: C2H2 */
        species[7] =
            +8.08681094e-01
            +1.16807815e-02 * tc[1]
            -1.18390605e-05 * tc[2]
            +7.00381092e-09 * tc[3]
            -1.70014595e-12 * tc[4]
            +2.64289807e+04 * invT;
        /*species 8: CO */
        species[8] =
            +3.26245200e+00
            +7.55970500e-04 * tc[1]
            -1.29391833e-06 * tc[2]
            +1.39548600e-09 * tc[3]
            -4.94990200e-13 * tc[4]
            -1.43105400e+04 * invT;
        /*species 9: C2H4 */
        species[9] =
            +3.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 * invT;
        /*species 10: C2H5 */
        species[10] =
            +4.30585800e+00
            -2.09168190e-03 * tc[1]
            +1.65690900e-05 * tc[2]
            -1.49764685e-08 * tc[3]
            +4.60969560e-12 * tc[4]
            +1.28417140e+04 * invT;
        /*species 12: C2H6 */
        species[12] =
            +1.46253900e+00
            +7.74733500e-03 * tc[1]
            +1.92683567e-06 * tc[2]
            -3.14458000e-09 * tc[3]
            +9.17253400e-13 * tc[4]
            -1.12391800e+04 * invT;
        /*species 13: CH3O */
        species[13] =
            +2.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122200e-13 * tc[4]
            +9.78601100e+02 * invT;
        /*species 14: O2 */
        species[14] =
            +3.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 * invT;
        /*species 15: HO2 */
        species[15] =
            +4.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 16: H2O2 */
        species[16] =
            +3.38875400e+00
            +3.28461300e-03 * tc[1]
            -4.95004333e-08 * tc[2]
            -1.15645150e-09 * tc[3]
            +4.94303000e-13 * tc[4]
            -1.76631500e+04 * invT;
        /*species 17: CO2 */
        species[17] =
            +2.27572500e+00
            +4.96103600e-03 * tc[1]
            -3.46970333e-06 * tc[2]
            +1.71667175e-09 * tc[3]
            -4.23456000e-13 * tc[4]
            -4.83731400e+04 * invT;
        /*species 18: CH3HCO */
        species[18] =
            +4.72945950e+00
            -1.59664290e-03 * tc[1]
            +1.58449737e-05 * tc[2]
            -1.43646527e-08 * tc[3]
            +4.38622240e-12 * tc[4]
            -2.15728780e+04 * invT;
        /*species 29: N2 */
        species[29] =
            +3.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
    } else {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +2.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 * invT;
        /*species 2: CH3 */
        species[2] =
            +2.97812060e+00
            +2.89892600e-03 * tc[1]
            -6.58526667e-07 * tc[2]
            +7.68244750e-11 * tc[3]
            -3.58348320e-15 * tc[4]
            +1.65095130e+04 * invT;
        /*species 3: O */
        species[3] =
            +2.54206000e+00
            -1.37753100e-05 * tc[1]
            -1.03426767e-09 * tc[2]
            +1.13776675e-12 * tc[3]
            -8.73610400e-17 * tc[4]
            +2.92308000e+04 * invT;
        /*species 4: CH4 */
        species[4] =
            +1.68347900e+00
            +5.11862000e-03 * tc[1]
            -1.29170967e-06 * tc[2]
            +1.69639625e-10 * tc[3]
            -9.00684600e-15 * tc[4]
            -1.00807900e+04 * invT;
        /*species 5: OH */
        species[5] =
            +2.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.68362875e+03 * invT;
        /*species 6: H2O */
        species[6] =
            +2.67214600e+00
            +1.52814650e-03 * tc[1]
            -2.91008667e-07 * tc[2]
            +3.00249000e-11 * tc[3]
            -1.27832360e-15 * tc[4]
            -2.98992100e+04 * invT;
        /*species 7: C2H2 */
        species[7] =
            +4.14756964e+00
            +2.98083332e-03 * tc[1]
            -7.90982840e-07 * tc[2]
            +1.16853043e-10 * tc[3]
            -7.22470426e-15 * tc[4]
            +2.59359992e+04 * invT;
        /*species 8: CO */
        species[8] =
            +3.02507800e+00
            +7.21344500e-04 * tc[1]
            -1.87694267e-07 * tc[2]
            +2.54645250e-11 * tc[3]
            -1.38219040e-15 * tc[4]
            -1.42683500e+04 * invT;
        /*species 9: C2H4 */
        species[9] =
            +2.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 * invT;
        /*species 10: C2H5 */
        species[10] =
            +4.28788140e+00
            +6.21694650e-03 * tc[1]
            -1.47130397e-06 * tc[2]
            +1.76635255e-10 * tc[3]
            -8.40702720e-15 * tc[4]
            +1.20564550e+04 * invT;
        /*species 12: C2H6 */
        species[12] =
            +4.82593800e+00
            +6.92021500e-03 * tc[1]
            -1.51908633e-06 * tc[2]
            +1.68124175e-10 * tc[3]
            -7.19632200e-15 * tc[4]
            -1.27177900e+04 * invT;
        /*species 13: CH3O */
        species[13] =
            +3.77080000e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832500e+02 * invT;
        /*species 14: O2 */
        species[14] =
            +3.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 * invT;
        /*species 15: HO2 */
        species[15] =
            +4.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 16: H2O2 */
        species[16] =
            +4.57316700e+00
            +2.16806800e-03 * tc[1]
            -4.91563000e-07 * tc[2]
            +5.87226000e-11 * tc[3]
            -2.86330800e-15 * tc[4]
            -1.80069600e+04 * invT;
        /*species 17: CO2 */
        species[17] =
            +4.45362300e+00
            +1.57008450e-03 * tc[1]
            -4.26137000e-07 * tc[2]
            +5.98499250e-11 * tc[3]
            -3.33806600e-15 * tc[4]
            -4.89669600e+04 * invT;
        /*species 18: CH3HCO */
        species[18] =
            +5.40411080e+00
            +5.86152950e-03 * tc[1]
            -1.40877123e-06 * tc[2]
            +1.70931128e-10 * tc[3]
            -8.19697260e-15 * tc[4]
            -2.25931220e+04 * invT;
        /*species 29: N2 */
        species[29] =
            +2.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }

    /*species with midpoint at T=1387 kelvin */
    if (T < 1387) {
        /*species 27: HO2CH2OCHO */
        species[27] =
            +3.47935703e+00
            +2.01476196e-02 * tc[1]
            -1.10036432e-05 * tc[2]
            +3.35900293e-09 * tc[3]
            -4.37203160e-13 * tc[4]
            -5.80629934e+04 * invT;
    } else {
        /*species 27: HO2CH2OCHO */
        species[27] =
            +1.64584298e+01
            +4.26341756e-03 * tc[1]
            -1.01371167e-06 * tc[2]
            +1.21399227e-10 * tc[3]
            -5.74632668e-15 * tc[4]
            -6.23959608e+04 * invT;
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 26: CH3OCH2O2 */
        species[26] =
            +2.21029612e+00
            +1.84438727e-02 * tc[1]
            -9.41871850e-06 * tc[2]
            +2.89326332e-09 * tc[3]
            -3.94260940e-13 * tc[4]
            -1.94940940e+04 * invT;
    } else {
        /*species 26: CH3OCH2O2 */
        species[26] =
            +1.24249729e+01
            +5.93529930e-03 * tc[1]
            -1.35968844e-06 * tc[2]
            +1.58827702e-10 * tc[3]
            -7.38855734e-15 * tc[4]
            -2.29679238e+04 * invT;
    }

    /*species with midpoint at T=1200 kelvin */
    if (T < 1200) {
        /*species 11: CH2O */
        species[11] =
            +2.69626120e+00
            +2.46307115e-03 * tc[1]
            +2.76088313e-07 * tc[2]
            -1.37595490e-10 * tc[3]
            -7.92206520e-14 * tc[4]
            -1.49707930e+04 * invT;
    } else {
        /*species 11: CH2O */
        species[11] =
            +5.14819050e+00
            +1.43390080e-03 * tc[1]
            -7.92754433e-08 * tc[2]
            -4.02782575e-11 * tc[3]
            +5.71334700e-15 * tc[4]
            -1.62301730e+04 * invT;
    }

    /*species with midpoint at T=1362 kelvin */
    if (T < 1362) {
        /*species 21: CH3OCO */
        species[21] =
            +3.94199159e+00
            +1.21717442e-02 * tc[1]
            -5.51985200e-06 * tc[2]
            +1.14634353e-09 * tc[3]
            -6.63591416e-14 * tc[4]
            -2.14404829e+04 * invT;
    } else {
        /*species 21: CH3OCO */
        species[21] =
            +1.30877600e+01
            +2.26772475e-03 * tc[1]
            -5.50321213e-07 * tc[2]
            +6.67993193e-11 * tc[3]
            -3.19153726e-15 * tc[4]
            -2.46616400e+04 * invT;
    }

    /*species with midpoint at T=1603 kelvin */
    if (T < 1603) {
        /*species 25: HOCH2OCO */
        species[25] =
            +6.08180801e+00
            +6.43841795e-03 * tc[1]
            +6.81398060e-07 * tc[2]
            -1.52538730e-09 * tc[3]
            +3.59641118e-13 * tc[4]
            -4.39526183e+04 * invT;
    } else {
        /*species 25: HOCH2OCO */
        species[25] =
            +1.13737391e+01
            +4.08831949e-03 * tc[1]
            -9.73446737e-07 * tc[2]
            +1.16673904e-10 * tc[3]
            -5.52553646e-15 * tc[4]
            -4.65575743e+04 * invT;
    }

    /*species with midpoint at T=1686 kelvin */
    if (T < 1686) {
        /*species 22: CH3OCHO */
        species[22] =
            +3.08839783e+00
            +1.01880024e-02 * tc[1]
            -2.28259013e-06 * tc[2]
            -1.82046551e-10 * tc[3]
            +1.12426043e-13 * tc[4]
            -4.41855167e+04 * invT;
    } else {
        /*species 22: CH3OCHO */
        species[22] =
            +8.69123518e+00
            +5.77515610e-03 * tc[1]
            -1.42594162e-06 * tc[2]
            +1.75633265e-10 * tc[3]
            -8.48667104e-15 * tc[4]
            -4.64364769e+04 * invT;
    }

    /*species with midpoint at T=1402 kelvin */
    if (T < 1402) {
        /*species 28: O2CH2OCH2O2H */
        species[28] =
            +1.99640551e+00
            +2.91613116e-02 * tc[1]
            -1.84419926e-05 * tc[2]
            +6.49526350e-09 * tc[3]
            -9.54282010e-13 * tc[4]
            -3.27628742e+04 * invT;
    } else {
        /*species 28: O2CH2OCH2O2H */
        species[28] =
            +1.92038046e+01
            +5.21974205e-03 * tc[1]
            -1.20194313e-06 * tc[2]
            +1.40948211e-10 * tc[3]
            -6.57614428e-15 * tc[4]
            -3.79207055e+04 * invT;
    }

    /*species with midpoint at T=2014 kelvin */
    if (T < 2014) {
        /*species 23: CH3OCH2OH */
        species[23] =
            +3.15851876e+00
            +1.22162875e-02 * tc[1]
            -2.88994928e-06 * tc[2]
            -1.48329832e-11 * tc[3]
            +8.72800006e-14 * tc[4]
            -4.54488899e+04 * invT;
    } else {
        /*species 23: CH3OCH2OH */
        species[23] =
            +8.70981570e+00
            +7.68011860e-03 * tc[1]
            -1.80334596e-06 * tc[2]
            +2.15143362e-10 * tc[3]
            -1.01763950e-14 * tc[4]
            -4.76607115e+04 * invT;
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEntropy(double * restrict  species, double * restrict  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1376 kelvin */
    if (T < 1376) {
        /*species 19: HCOOH */
        species[19] =
            +1.43548185e+00 * tc[0]
            +1.63363016e-02 * tc[1]
            -5.31287105e-06 * tc[2]
            +1.10710992e-09 * tc[3]
            -1.00544026e-13 * tc[4]
            +1.72885798e+01 ;
    } else {
        /*species 19: HCOOH */
        species[19] =
            +6.68733013e+00 * tc[0]
            +5.14289368e-03 * tc[1]
            -9.11192565e-07 * tc[2]
            +9.65730543e-11 * tc[3]
            -4.27230498e-15 * tc[4]
            -1.13104798e+01 ;
    }

    /*species with midpoint at T=1475 kelvin */
    if (T < 1475) {
        /*species 24: OCH2OCHO */
        species[24] =
            +5.19690837e+00 * tc[0]
            +1.58839723e-02 * tc[1]
            +1.76770274e-07 * tc[2]
            -2.03485641e-09 * tc[3]
            +4.86654503e-13 * tc[4]
            +6.11645828e+00 ;
    } else {
        /*species 24: OCH2OCHO */
        species[24] =
            +1.20233916e+01 * tc[0]
            +8.11262659e-03 * tc[1]
            -1.45678231e-06 * tc[2]
            +1.55780128e-10 * tc[3]
            -6.93438812e-15 * tc[4]
            -3.33691809e+01 ;
    }

    /*species with midpoint at T=710 kelvin */
    if (T < 710) {
        /*species 20: CH3OCH3 */
        species[20] =
            +5.68097447e+00 * tc[0]
            -5.39434751e-03 * tc[1]
            +3.24736375e-05 * tc[2]
            -2.68355106e-08 * tc[3]
            +8.18685045e-12 * tc[4]
            -6.36955496e-01 ;
    } else {
        /*species 20: CH3OCH3 */
        species[20] =
            +8.30815546e-01 * tc[0]
            +2.69173263e-02 * tc[1]
            -6.94373885e-06 * tc[2]
            +1.15838360e-09 * tc[3]
            -8.54266960e-14 * tc[4]
            +2.02174360e+01 ;
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117600e-01 ;
        /*species 1: H2 */
        species[1] =
            +3.29812400e+00 * tc[0]
            +8.24944200e-04 * tc[1]
            -4.07150750e-07 * tc[2]
            -3.15847800e-11 * tc[3]
            +1.03371800e-13 * tc[4]
            -3.29409400e+00 ;
        /*species 2: CH3 */
        species[2] =
            +3.65717970e+00 * tc[0]
            +2.12659790e-03 * tc[1]
            +2.72919415e-06 * tc[2]
            -2.20603343e-09 * tc[3]
            +6.16426850e-13 * tc[4]
            +1.67353540e+00 ;
        /*species 3: O */
        species[3] =
            +2.94642900e+00 * tc[0]
            -1.63816600e-03 * tc[1]
            +1.21051600e-06 * tc[2]
            -5.34281000e-10 * tc[3]
            +9.72674000e-14 * tc[4]
            +2.96399500e+00 ;
        /*species 4: CH4 */
        species[4] =
            +7.78741500e-01 * tc[0]
            +1.74766800e-02 * tc[1]
            -1.39170450e-05 * tc[2]
            +1.01656933e-08 * tc[3]
            -3.05982750e-12 * tc[4]
            +1.37221900e+01 ;
        /*species 5: OH */
        species[5] =
            +4.12530561e+00 * tc[0]
            -3.22544939e-03 * tc[1]
            +3.26382346e-06 * tc[2]
            -1.93284548e-09 * tc[3]
            +5.15593447e-13 * tc[4]
            -6.90432960e-01 ;
        /*species 6: H2O */
        species[6] =
            +3.38684200e+00 * tc[0]
            +3.47498200e-03 * tc[1]
            -3.17734800e-06 * tc[2]
            +2.32286033e-09 * tc[3]
            -6.26647000e-13 * tc[4]
            +2.59023300e+00 ;
        /*species 7: C2H2 */
        species[7] =
            +8.08681094e-01 * tc[0]
            +2.33615629e-02 * tc[1]
            -1.77585907e-05 * tc[2]
            +9.33841457e-09 * tc[3]
            -2.12518243e-12 * tc[4]
            +1.39397051e+01 ;
        /*species 8: CO */
        species[8] =
            +3.26245200e+00 * tc[0]
            +1.51194100e-03 * tc[1]
            -1.94087750e-06 * tc[2]
            +1.86064800e-09 * tc[3]
            -6.18737750e-13 * tc[4]
            +4.84889700e+00 ;
        /*species 9: C2H4 */
        species[9] =
            +3.95920148e+00 * tc[0]
            -7.57052247e-03 * tc[1]
            +2.85495146e-05 * tc[2]
            -2.30529584e-08 * tc[3]
            +6.74710933e-12 * tc[4]
            +4.09733096e+00 ;
        /*species 10: C2H5 */
        species[10] =
            +4.30585800e+00 * tc[0]
            -4.18336380e-03 * tc[1]
            +2.48536350e-05 * tc[2]
            -1.99686247e-08 * tc[3]
            +5.76211950e-12 * tc[4]
            +4.71002360e+00 ;
        /*species 12: C2H6 */
        species[12] =
            +1.46253900e+00 * tc[0]
            +1.54946700e-02 * tc[1]
            +2.89025350e-06 * tc[2]
            -4.19277333e-09 * tc[3]
            +1.14656675e-12 * tc[4]
            +1.44322900e+01 ;
        /*species 13: CH3O */
        species[13] =
            +2.10620400e+00 * tc[0]
            +7.21659500e-03 * tc[1]
            +2.66923600e-06 * tc[2]
            -2.45921200e-09 * tc[3]
            +5.18902750e-13 * tc[4]
            +1.31521800e+01 ;
        /*species 14: O2 */
        species[14] =
            +3.21293600e+00 * tc[0]
            +1.12748600e-03 * tc[1]
            -2.87807500e-07 * tc[2]
            +4.37959000e-10 * tc[3]
            -2.19213850e-13 * tc[4]
            +6.03473800e+00 ;
        /*species 15: HO2 */
        species[15] =
            +4.30179801e+00 * tc[0]
            -4.74912051e-03 * tc[1]
            +1.05791445e-05 * tc[2]
            -8.09212980e-09 * tc[3]
            +2.32306281e-12 * tc[4]
            +3.71666245e+00 ;
        /*species 16: H2O2 */
        species[16] =
            +3.38875400e+00 * tc[0]
            +6.56922600e-03 * tc[1]
            -7.42506500e-08 * tc[2]
            -1.54193533e-09 * tc[3]
            +6.17878750e-13 * tc[4]
            +6.78536300e+00 ;
        /*species 17: CO2 */
        species[17] =
            +2.27572500e+00 * tc[0]
            +9.92207200e-03 * tc[1]
            -5.20455500e-06 * tc[2]
            +2.28889567e-09 * tc[3]
            -5.29320000e-13 * tc[4]
            +1.01884900e+01 ;
        /*species 18: CH3HCO */
        species[18] =
            +4.72945950e+00 * tc[0]
            -3.19328580e-03 * tc[1]
            +2.37674605e-05 * tc[2]
            -1.91528703e-08 * tc[3]
            +5.48277800e-12 * tc[4]
            +4.10301590e+00 ;
        /*species 29: N2 */
        species[29] =
            +3.29867700e+00 * tc[0]
            +1.40824000e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213750e-13 * tc[4]
            +3.95037200e+00 ;
    } else {
        /*species 0: H */
        species[0] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117600e-01 ;
        /*species 1: H2 */
        species[1] =
            +2.99142300e+00 * tc[0]
            +7.00064400e-04 * tc[1]
            -2.81691450e-08 * tc[2]
            -3.07719267e-12 * tc[3]
            +3.95688000e-16 * tc[4]
            -1.35511000e+00 ;
        /*species 2: CH3 */
        species[2] =
            +2.97812060e+00 * tc[0]
            +5.79785200e-03 * tc[1]
            -9.87790000e-07 * tc[2]
            +1.02432633e-10 * tc[3]
            -4.47935400e-15 * tc[4]
            +4.72247990e+00 ;
        /*species 3: O */
        species[3] =
            +2.54206000e+00 * tc[0]
            -2.75506200e-05 * tc[1]
            -1.55140150e-09 * tc[2]
            +1.51702233e-12 * tc[3]
            -1.09201300e-16 * tc[4]
            +4.92030800e+00 ;
        /*species 4: CH4 */
        species[4] =
            +1.68347900e+00 * tc[0]
            +1.02372400e-02 * tc[1]
            -1.93756450e-06 * tc[2]
            +2.26186167e-10 * tc[3]
            -1.12585575e-14 * tc[4]
            +9.62339500e+00 ;
        /*species 5: OH */
        species[5] =
            +2.86472886e+00 * tc[0]
            +1.05650448e-03 * tc[1]
            -1.29541379e-07 * tc[2]
            +1.01739558e-11 * tc[3]
            -3.32989690e-16 * tc[4]
            +5.70164073e+00 ;
        /*species 6: H2O */
        species[6] =
            +2.67214600e+00 * tc[0]
            +3.05629300e-03 * tc[1]
            -4.36513000e-07 * tc[2]
            +4.00332000e-11 * tc[3]
            -1.59790450e-15 * tc[4]
            +6.86281700e+00 ;
        /*species 7: C2H2 */
        species[7] =
            +4.14756964e+00 * tc[0]
            +5.96166664e-03 * tc[1]
            -1.18647426e-06 * tc[2]
            +1.55804057e-10 * tc[3]
            -9.03088033e-15 * tc[4]
            -1.23028121e+00 ;
        /*species 8: CO */
        species[8] =
            +3.02507800e+00 * tc[0]
            +1.44268900e-03 * tc[1]
            -2.81541400e-07 * tc[2]
            +3.39527000e-11 * tc[3]
            -1.72773800e-15 * tc[4]
            +6.10821800e+00 ;
        /*species 9: C2H4 */
        species[9] =
            +2.03611116e+00 * tc[0]
            +1.46454151e-02 * tc[1]
            -3.35538958e-06 * tc[2]
            +4.90743077e-10 * tc[3]
            -3.14265152e-14 * tc[4]
            +1.03053693e+01 ;
        /*species 10: C2H5 */
        species[10] =
            +4.28788140e+00 * tc[0]
            +1.24338930e-02 * tc[1]
            -2.20695595e-06 * tc[2]
            +2.35513673e-10 * tc[3]
            -1.05087840e-14 * tc[4]
            +8.46025830e-01 ;
        /*species 12: C2H6 */
        species[12] =
            +4.82593800e+00 * tc[0]
            +1.38404300e-02 * tc[1]
            -2.27862950e-06 * tc[2]
            +2.24165567e-10 * tc[3]
            -8.99540250e-15 * tc[4]
            -5.23950700e+00 ;
        /*species 13: CH3O */
        species[13] =
            +3.77080000e+00 * tc[0]
            +7.87149700e-03 * tc[1]
            -1.32819200e-06 * tc[2]
            +1.31481033e-10 * tc[3]
            -5.28154000e-15 * tc[4]
            +2.92957500e+00 ;
        /*species 14: O2 */
        species[14] =
            +3.69757800e+00 * tc[0]
            +6.13519700e-04 * tc[1]
            -6.29421000e-08 * tc[2]
            +5.91760333e-12 * tc[3]
            -2.84108750e-16 * tc[4]
            +3.18916600e+00 ;
        /*species 15: HO2 */
        species[15] =
            +4.01721090e+00 * tc[0]
            +2.23982013e-03 * tc[1]
            -3.16829075e-07 * tc[2]
            +3.80821233e-11 * tc[3]
            -2.69771337e-15 * tc[4]
            +3.78510215e+00 ;
        /*species 16: H2O2 */
        species[16] =
            +4.57316700e+00 * tc[0]
            +4.33613600e-03 * tc[1]
            -7.37344500e-07 * tc[2]
            +7.82968000e-11 * tc[3]
            -3.57913500e-15 * tc[4]
            +5.01137000e-01 ;
        /*species 17: CO2 */
        species[17] =
            +4.45362300e+00 * tc[0]
            +3.14016900e-03 * tc[1]
            -6.39205500e-07 * tc[2]
            +7.97999000e-11 * tc[3]
            -4.17258250e-15 * tc[4]
            -9.55395900e-01 ;
        /*species 18: CH3HCO */
        species[18] =
            +5.40411080e+00 * tc[0]
            +1.17230590e-02 * tc[1]
            -2.11315685e-06 * tc[2]
            +2.27908170e-10 * tc[3]
            -1.02462158e-14 * tc[4]
            -3.48079170e+00 ;
        /*species 29: N2 */
        species[29] =
            +2.92664000e+00 * tc[0]
            +1.48797700e-03 * tc[1]
            -2.84238050e-07 * tc[2]
            +3.36568000e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
    }

    /*species with midpoint at T=1387 kelvin */
    if (T < 1387) {
        /*species 27: HO2CH2OCHO */
        species[27] =
            +3.47935703e+00 * tc[0]
            +4.02952392e-02 * tc[1]
            -1.65054648e-05 * tc[2]
            +4.47867057e-09 * tc[3]
            -5.46503950e-13 * tc[4]
            +1.52521392e+01 ;
    } else {
        /*species 27: HO2CH2OCHO */
        species[27] =
            +1.64584298e+01 * tc[0]
            +8.52683511e-03 * tc[1]
            -1.52056750e-06 * tc[2]
            +1.61865636e-10 * tc[3]
            -7.18290835e-15 * tc[4]
            -5.38924139e+01 ;
    }

    /*species with midpoint at T=1389 kelvin */
    if (T < 1389) {
        /*species 26: CH3OCH2O2 */
        species[26] =
            +2.21029612e+00 * tc[0]
            +3.68877454e-02 * tc[1]
            -1.41280778e-05 * tc[2]
            +3.85768443e-09 * tc[3]
            -4.92826175e-13 * tc[4]
            +1.91463601e+01 ;
    } else {
        /*species 26: CH3OCH2O2 */
        species[26] =
            +1.24249729e+01 * tc[0]
            +1.18705986e-02 * tc[1]
            -2.03953266e-06 * tc[2]
            +2.11770270e-10 * tc[3]
            -9.23569667e-15 * tc[4]
            -3.53740145e+01 ;
    }

    /*species with midpoint at T=1200 kelvin */
    if (T < 1200) {
        /*species 11: CH2O */
        species[11] =
            +2.69626120e+00 * tc[0]
            +4.92614230e-03 * tc[1]
            +4.14132470e-07 * tc[2]
            -1.83460653e-10 * tc[3]
            -9.90258150e-14 * tc[4]
            +9.46975990e+00 ;
    } else {
        /*species 11: CH2O */
        species[11] =
            +5.14819050e+00 * tc[0]
            +2.86780160e-03 * tc[1]
            -1.18913165e-07 * tc[2]
            -5.37043433e-11 * tc[3]
            +7.14168375e-15 * tc[4]
            -5.12138130e+00 ;
    }

    /*species with midpoint at T=1362 kelvin */
    if (T < 1362) {
        /*species 21: CH3OCO */
        species[21] =
            +3.94199159e+00 * tc[0]
            +2.43434884e-02 * tc[1]
            -8.27977800e-06 * tc[2]
            +1.52845804e-09 * tc[3]
            -8.29489270e-14 * tc[4]
            +1.66954362e+01 ;
    } else {
        /*species 21: CH3OCO */
        species[21] =
            +1.30877600e+01 * tc[0]
            +4.53544950e-03 * tc[1]
            -8.25481820e-07 * tc[2]
            +8.90657590e-11 * tc[3]
            -3.98942157e-15 * tc[4]
            -3.27914051e+01 ;
    }

    /*species with midpoint at T=1603 kelvin */
    if (T < 1603) {
        /*species 25: HOCH2OCO */
        species[25] =
            +6.08180801e+00 * tc[0]
            +1.28768359e-02 * tc[1]
            +1.02209709e-06 * tc[2]
            -2.03384974e-09 * tc[3]
            +4.49551398e-13 * tc[4]
            +2.54054449e+00 ;
    } else {
        /*species 25: HOCH2OCO */
        species[25] =
            +1.13737391e+01 * tc[0]
            +8.17663898e-03 * tc[1]
            -1.46017011e-06 * tc[2]
            +1.55565205e-10 * tc[3]
            -6.90692057e-15 * tc[4]
            -2.86035265e+01 ;
    }

    /*species with midpoint at T=1686 kelvin */
    if (T < 1686) {
        /*species 22: CH3OCHO */
        species[22] =
            +3.08839783e+00 * tc[0]
            +2.03760048e-02 * tc[1]
            -3.42388520e-06 * tc[2]
            -2.42728734e-10 * tc[3]
            +1.40532554e-13 * tc[4]
            +1.25364719e+01 ;
    } else {
        /*species 22: CH3OCHO */
        species[22] =
            +8.69123518e+00 * tc[0]
            +1.15503122e-02 * tc[1]
            -2.13891243e-06 * tc[2]
            +2.34177686e-10 * tc[3]
            -1.06083388e-14 * tc[4]
            -1.89301478e+01 ;
    }

    /*species with midpoint at T=1402 kelvin */
    if (T < 1402) {
        /*species 28: O2CH2OCH2O2H */
        species[28] =
            +1.99640551e+00 * tc[0]
            +5.83226232e-02 * tc[1]
            -2.76629889e-05 * tc[2]
            +8.66035133e-09 * tc[3]
            -1.19285251e-12 * tc[4]
            +2.44215005e+01 ;
    } else {
        /*species 28: O2CH2OCH2O2H */
        species[28] =
            +1.92038046e+01 * tc[0]
            +1.04394841e-02 * tc[1]
            -1.80291469e-06 * tc[2]
            +1.87930948e-10 * tc[3]
            -8.22018035e-15 * tc[4]
            -6.51847273e+01 ;
    }

    /*species with midpoint at T=2014 kelvin */
    if (T < 2014) {
        /*species 23: CH3OCH2OH */
        species[23] =
            +3.15851876e+00 * tc[0]
            +2.44325751e-02 * tc[1]
            -4.33492392e-06 * tc[2]
            -1.97773109e-11 * tc[3]
            +1.09100001e-13 * tc[4]
            +1.30511235e+01 ;
    } else {
        /*species 23: CH3OCH2OH */
        species[23] =
            +8.70981570e+00 * tc[0]
            +1.53602372e-02 * tc[1]
            -2.70501894e-06 * tc[2]
            +2.86857815e-10 * tc[3]
            -1.27204938e-14 * tc[4]
            -1.80226702e+01 ;
    }
    return;
}


/*save molecular weights into array */
void molecularWeight(double * restrict  wt)
{
    wt[0] = 1.007970; /*H */
    wt[1] = 2.015940; /*H2 */
    wt[2] = 15.035060; /*CH3 */
    wt[3] = 15.999400; /*O */
    wt[4] = 16.043030; /*CH4 */
    wt[5] = 17.007370; /*OH */
    wt[6] = 18.015340; /*H2O */
    wt[7] = 26.038240; /*C2H2 */
    wt[8] = 28.010550; /*CO */
    wt[9] = 28.054180; /*C2H4 */
    wt[10] = 29.062150; /*C2H5 */
    wt[11] = 30.026490; /*CH2O */
    wt[12] = 30.070120; /*C2H6 */
    wt[13] = 31.034460; /*CH3O */
    wt[14] = 31.998800; /*O2 */
    wt[15] = 33.006770; /*HO2 */
    wt[16] = 34.014740; /*H2O2 */
    wt[17] = 44.009950; /*CO2 */
    wt[18] = 44.053580; /*CH3HCO */
    wt[19] = 46.025890; /*HCOOH */
    wt[20] = 46.069520; /*CH3OCH3 */
    wt[21] = 59.045010; /*CH3OCO */
    wt[22] = 60.052980; /*CH3OCHO */
    wt[23] = 62.068920; /*CH3OCH2OH */
    wt[24] = 75.044410; /*OCH2OCHO */
    wt[25] = 75.044410; /*HOCH2OCO */
    wt[26] = 77.060350; /*CH3OCH2O2 */
    wt[27] = 92.051780; /*HO2CH2OCHO */
    wt[28] = 109.059150; /*O2CH2OCH2O2H */
    wt[29] = 28.013400; /*N2 */

    return;
}
/* get temperature given internal energy in mass units and mass fracs */
void GET_T_GIVEN_EY(double * restrict  e, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const double tol  = 1.e-12;
#else
    const int maxiter = 200;
    const double tol  = 1.e-6;
#endif
    double ein  = *e;
    double tmin = 250;/*max lower bound for thermo def */
    double tmax = 3500;/*min upper bound for thermo def */
    double e1,emin,emax,cv,t1,dt;
    int i;/* loop counter */
    CKUBMS(&tmin, y, iwrk, rwrk, &emin);
    CKUBMS(&tmax, y, iwrk, rwrk, &emax);
    if (ein < emin) {
        /*Linear Extrapolation below tmin */
        CKCVBS(&tmin, y, iwrk, rwrk, &cv);
        *t = tmin - (emin-ein)/cv;
        *ierr = 1;
        return;
    }
    if (ein > emax) {
        /*Linear Extrapolation above tmax */
        CKCVBS(&tmax, y, iwrk, rwrk, &cv);
        *t = tmax - (emax-ein)/cv;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKUBMS(&t1,y,iwrk,rwrk,&e1);
        CKCVBS(&t1,y,iwrk,rwrk,&cv);
        dt = (ein - e1) / cv;
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
    }
    *t = t1;
    *ierr = 0;
    return;
}

/* End of file  */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENIMC EGTRANSETLENIMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENIMC egtransetlenimc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENIMC egtransetlenimc_
#endif
void egtransetLENIMC(int* LENIMC) {
  *LENIMC =          122;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENRMC EGTRANSETLENRMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENRMC egtransetlenrmc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENRMC egtransetlenrmc_
#endif
void egtransetLENRMC(int* LENRMC) {
  *LENRMC =        18150;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNO EGTRANSETNO
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNO egtransetno
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNO egtransetno_
#endif
void egtransetNO(int* NO) {
  *NO =            4;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetKK EGTRANSETKK
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetKK egtransetkk
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetKK egtransetkk_
#endif
void egtransetKK(int* KK) {
  *KK =           30;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLITE EGTRANSETNLITE
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLITE egtransetnlite
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLITE egtransetnlite_
#endif
void egtransetNLITE(int* NLITE) {
  *NLITE =            2;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetPATM EGTRANSETPATM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetPATM egtransetpatm
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetPATM egtransetpatm_
#endif
void egtransetPATM(double* PATM) {
  *PATM =   0.1013250000000000E+07;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetWT EGTRANSETWT
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetWT egtransetwt
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetWT egtransetwt_
#endif
void egtransetWT(double* WT) {
  WT[           0] =   0.1007969975471497E+01;
  WT[           1] =   0.2015939950942993E+01;
  WT[           2] =   0.1503506028652191E+02;
  WT[           3] =   0.1599940013885498E+02;
  WT[           4] =   0.1604303026199341E+02;
  WT[           5] =   0.1700737011432648E+02;
  WT[           6] =   0.1801534008979797E+02;
  WT[           7] =   0.2603824067115784E+02;
  WT[           8] =   0.2801055049896240E+02;
  WT[           9] =   0.2805418062210083E+02;
  WT[          10] =   0.2906215059757233E+02;
  WT[          11] =   0.3002649044990540E+02;
  WT[          12] =   0.3007012057304382E+02;
  WT[          13] =   0.3103446042537689E+02;
  WT[          14] =   0.3199880027770996E+02;
  WT[          15] =   0.3300677025318146E+02;
  WT[          16] =   0.3401474022865295E+02;
  WT[          17] =   0.4400995063781738E+02;
  WT[          18] =   0.4405358076095581E+02;
  WT[          19] =   0.4602589058876038E+02;
  WT[          20] =   0.4606952071189880E+02;
  WT[          21] =   0.5904501092433929E+02;
  WT[          22] =   0.6005298089981079E+02;
  WT[          23] =   0.6206892085075378E+02;
  WT[          24] =   0.7504441106319427E+02;
  WT[          25] =   0.7504441106319427E+02;
  WT[          26] =   0.7706035101413727E+02;
  WT[          27] =   0.9205178117752075E+02;
  WT[          28] =   0.1090591512918472E+03;
  WT[          29] =   0.2801339912414551E+02;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetEPS EGTRANSETEPS
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetEPS egtranseteps
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetEPS egtranseteps_
#endif
void egtransetEPS(double* EPS) {
  EPS[           0] =   0.1450000000000000E+03;
  EPS[           1] =   0.3800000000000000E+02;
  EPS[           2] =   0.1440000000000000E+03;
  EPS[           3] =   0.8000000000000000E+02;
  EPS[           4] =   0.1414000000000000E+03;
  EPS[           5] =   0.8000000000000000E+02;
  EPS[           6] =   0.5724000000000000E+03;
  EPS[           7] =   0.2653000000000000E+03;
  EPS[           8] =   0.9809999999999999E+02;
  EPS[           9] =   0.2384000000000000E+03;
  EPS[          10] =   0.2475000000000000E+03;
  EPS[          11] =   0.4980000000000000E+03;
  EPS[          12] =   0.2475000000000000E+03;
  EPS[          13] =   0.4170000000000000E+03;
  EPS[          14] =   0.1074000000000000E+03;
  EPS[          15] =   0.1074000000000000E+03;
  EPS[          16] =   0.1074000000000000E+03;
  EPS[          17] =   0.2440000000000000E+03;
  EPS[          18] =   0.4360000000000000E+03;
  EPS[          19] =   0.4706000000000000E+03;
  EPS[          20] =   0.3294000000000000E+03;
  EPS[          21] =   0.4065000000000000E+03;
  EPS[          22] =   0.4065000000000000E+03;
  EPS[          23] =   0.3294000000000000E+03;
  EPS[          24] =   0.3294000000000000E+03;
  EPS[          25] =   0.3294000000000000E+03;
  EPS[          26] =   0.3294000000000000E+03;
  EPS[          27] =   0.3294000000000000E+03;
  EPS[          28] =   0.3294000000000000E+03;
  EPS[          29] =   0.9753000000000000E+02;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetSIG EGTRANSETSIG
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetSIG egtransetsig
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetSIG egtransetsig_
#endif
void egtransetSIG(double* SIG) {
  SIG[           0] =   0.2050000000000000E+01;
  SIG[           1] =   0.2920000000000000E+01;
  SIG[           2] =   0.3800000000000000E+01;
  SIG[           3] =   0.2750000000000000E+01;
  SIG[           4] =   0.3746000000000000E+01;
  SIG[           5] =   0.2750000000000000E+01;
  SIG[           6] =   0.2605000000000000E+01;
  SIG[           7] =   0.3721000000000000E+01;
  SIG[           8] =   0.3650000000000000E+01;
  SIG[           9] =   0.3496000000000000E+01;
  SIG[          10] =   0.4350000000000000E+01;
  SIG[          11] =   0.3590000000000000E+01;
  SIG[          12] =   0.4350000000000000E+01;
  SIG[          13] =   0.3690000000000000E+01;
  SIG[          14] =   0.3458000000000000E+01;
  SIG[          15] =   0.3458000000000000E+01;
  SIG[          16] =   0.3458000000000000E+01;
  SIG[          17] =   0.3763000000000000E+01;
  SIG[          18] =   0.3970000000000000E+01;
  SIG[          19] =   0.4410000000000000E+01;
  SIG[          20] =   0.4624000000000000E+01;
  SIG[          21] =   0.4709000000000000E+01;
  SIG[          22] =   0.4709000000000000E+01;
  SIG[          23] =   0.4624000000000000E+01;
  SIG[          24] =   0.4624000000000000E+01;
  SIG[          25] =   0.4624000000000000E+01;
  SIG[          26] =   0.4624000000000000E+01;
  SIG[          27] =   0.4624000000000000E+01;
  SIG[          28] =   0.4624000000000000E+01;
  SIG[          29] =   0.3621000000000000E+01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetDIP EGTRANSETDIP
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetDIP egtransetdip
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetDIP egtransetdip_
#endif
void egtransetDIP(double* DIP) {
  DIP[           0] =   0.0000000000000000E+00;
  DIP[           1] =   0.0000000000000000E+00;
  DIP[           2] =   0.0000000000000000E+00;
  DIP[           3] =   0.0000000000000000E+00;
  DIP[           4] =   0.0000000000000000E+00;
  DIP[           5] =   0.0000000000000000E+00;
  DIP[           6] =   0.1844000000000000E+01;
  DIP[           7] =   0.0000000000000000E+00;
  DIP[           8] =   0.0000000000000000E+00;
  DIP[           9] =   0.0000000000000000E+00;
  DIP[          10] =   0.0000000000000000E+00;
  DIP[          11] =   0.0000000000000000E+00;
  DIP[          12] =   0.0000000000000000E+00;
  DIP[          13] =   0.1700000000000000E+01;
  DIP[          14] =   0.0000000000000000E+00;
  DIP[          15] =   0.0000000000000000E+00;
  DIP[          16] =   0.0000000000000000E+00;
  DIP[          17] =   0.0000000000000000E+00;
  DIP[          18] =   0.0000000000000000E+00;
  DIP[          19] =   0.0000000000000000E+00;
  DIP[          20] =   0.0000000000000000E+00;
  DIP[          21] =   0.0000000000000000E+00;
  DIP[          22] =   0.0000000000000000E+00;
  DIP[          23] =   0.0000000000000000E+00;
  DIP[          24] =   0.0000000000000000E+00;
  DIP[          25] =   0.0000000000000000E+00;
  DIP[          26] =   0.0000000000000000E+00;
  DIP[          27] =   0.0000000000000000E+00;
  DIP[          28] =   0.0000000000000000E+00;
  DIP[          29] =   0.0000000000000000E+00;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetPOL EGTRANSETPOL
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetPOL egtransetpol
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetPOL egtransetpol_
#endif
void egtransetPOL(double* POL) {
  POL[           0] =   0.0000000000000000E+00;
  POL[           1] =   0.7900000000000000E+00;
  POL[           2] =   0.0000000000000000E+00;
  POL[           3] =   0.0000000000000000E+00;
  POL[           4] =   0.2600000000000000E+01;
  POL[           5] =   0.0000000000000000E+00;
  POL[           6] =   0.0000000000000000E+00;
  POL[           7] =   0.0000000000000000E+00;
  POL[           8] =   0.1950000000000000E+01;
  POL[           9] =   0.0000000000000000E+00;
  POL[          10] =   0.0000000000000000E+00;
  POL[          11] =   0.0000000000000000E+00;
  POL[          12] =   0.0000000000000000E+00;
  POL[          13] =   0.0000000000000000E+00;
  POL[          14] =   0.1600000000000000E+01;
  POL[          15] =   0.0000000000000000E+00;
  POL[          16] =   0.0000000000000000E+00;
  POL[          17] =   0.2650000000000000E+01;
  POL[          18] =   0.0000000000000000E+00;
  POL[          19] =   0.0000000000000000E+00;
  POL[          20] =   0.0000000000000000E+00;
  POL[          21] =   0.0000000000000000E+00;
  POL[          22] =   0.0000000000000000E+00;
  POL[          23] =   0.0000000000000000E+00;
  POL[          24] =   0.0000000000000000E+00;
  POL[          25] =   0.0000000000000000E+00;
  POL[          26] =   0.0000000000000000E+00;
  POL[          27] =   0.0000000000000000E+00;
  POL[          28] =   0.0000000000000000E+00;
  POL[          29] =   0.1760000000000000E+01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetZROT EGTRANSETZROT
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetZROT egtransetzrot
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetZROT egtransetzrot_
#endif
void egtransetZROT(double* ZROT) {
  ZROT[           0] =   0.0000000000000000E+00;
  ZROT[           1] =   0.2800000000000000E+03;
  ZROT[           2] =   0.0000000000000000E+00;
  ZROT[           3] =   0.0000000000000000E+00;
  ZROT[           4] =   0.1300000000000000E+02;
  ZROT[           5] =   0.0000000000000000E+00;
  ZROT[           6] =   0.4000000000000000E+01;
  ZROT[           7] =   0.2500000000000000E+01;
  ZROT[           8] =   0.1800000000000000E+01;
  ZROT[           9] =   0.1500000000000000E+01;
  ZROT[          10] =   0.1500000000000000E+01;
  ZROT[          11] =   0.2000000000000000E+01;
  ZROT[          12] =   0.1500000000000000E+01;
  ZROT[          13] =   0.2000000000000000E+01;
  ZROT[          14] =   0.3800000000000000E+01;
  ZROT[          15] =   0.1000000000000000E+01;
  ZROT[          16] =   0.3800000000000000E+01;
  ZROT[          17] =   0.2100000000000000E+01;
  ZROT[          18] =   0.2000000000000000E+01;
  ZROT[          19] =   0.1500000000000000E+01;
  ZROT[          20] =   0.1000000000000000E+01;
  ZROT[          21] =   0.1000000000000000E+01;
  ZROT[          22] =   0.1000000000000000E+01;
  ZROT[          23] =   0.1000000000000000E+01;
  ZROT[          24] =   0.1000000000000000E+01;
  ZROT[          25] =   0.1000000000000000E+01;
  ZROT[          26] =   0.1000000000000000E+01;
  ZROT[          27] =   0.1000000000000000E+01;
  ZROT[          28] =   0.1000000000000000E+01;
  ZROT[          29] =   0.4000000000000000E+01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLIN EGTRANSETNLIN
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLIN egtransetnlin
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLIN egtransetnlin_
#endif
void egtransetNLIN(int* NLIN) {
  NLIN[           0] =            0;
  NLIN[           1] =            1;
  NLIN[           2] =            1;
  NLIN[           3] =            0;
  NLIN[           4] =            2;
  NLIN[           5] =            1;
  NLIN[           6] =            2;
  NLIN[           7] =            1;
  NLIN[           8] =            1;
  NLIN[           9] =            2;
  NLIN[          10] =            2;
  NLIN[          11] =            2;
  NLIN[          12] =            2;
  NLIN[          13] =            2;
  NLIN[          14] =            1;
  NLIN[          15] =            2;
  NLIN[          16] =            2;
  NLIN[          17] =            1;
  NLIN[          18] =            2;
  NLIN[          19] =            2;
  NLIN[          20] =            2;
  NLIN[          21] =            2;
  NLIN[          22] =            2;
  NLIN[          23] =            2;
  NLIN[          24] =            2;
  NLIN[          25] =            2;
  NLIN[          26] =            2;
  NLIN[          27] =            2;
  NLIN[          28] =            2;
  NLIN[          29] =            1;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFLAM EGTRANSETCOFLAM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFLAM egtransetcoflam
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFLAM egtransetcoflam_
#endif
void egtransetCOFLAM(double* COFLAM) {
  COFLAM[           0] =  -0.8554281253040258E+00;
  COFLAM[           1] =   0.3652691486294877E+01;
  COFLAM[           2] =  -0.3980303020607808E+00;
  COFLAM[           3] =   0.1757072886401263E-01;
  COFLAM[           4] =   0.9132344040594148E+01;
  COFLAM[           5] =  -0.4346961503128693E+00;
  COFLAM[           6] =   0.1126754364016490E+00;
  COFLAM[           7] =  -0.2531803165785582E-02;
  COFLAM[           8] =   0.1329540069466844E+02;
  COFLAM[           9] =  -0.4317035712196896E+01;
  COFLAM[          10] =   0.8574796350806355E+00;
  COFLAM[          11] =  -0.4516104542171843E-01;
  COFLAM[          12] =   0.1685025383414630E+01;
  COFLAM[          13] =   0.1929024677614472E+01;
  COFLAM[          14] =  -0.1738657445088021E+00;
  COFLAM[          15] =   0.7841476915125366E-02;
  COFLAM[          16] =   0.1295934751845324E+02;
  COFLAM[          17] =  -0.4856398414574292E+01;
  COFLAM[          18] =   0.1029040711986845E+01;
  COFLAM[          19] =  -0.5698710041428736E-01;
  COFLAM[          20] =   0.1500752966747290E+02;
  COFLAM[          21] =  -0.3630974904682330E+01;
  COFLAM[          22] =   0.5926158891396912E+00;
  COFLAM[          23] =  -0.2628250121345284E-01;
  COFLAM[          24] =   0.2350117885411275E+02;
  COFLAM[          25] =  -0.9060443719032351E+01;
  COFLAM[          26] =   0.1548701551684282E+01;
  COFLAM[          27] =  -0.7721043889928444E-01;
  COFLAM[          28] =  -0.9069304423379325E+01;
  COFLAM[          29] =   0.5065458511298031E+01;
  COFLAM[          30] =  -0.4576172762886371E+00;
  COFLAM[          31] =   0.1590924413166833E-01;
  COFLAM[          32] =   0.1157879307056280E+02;
  COFLAM[          33] =  -0.3021050899869008E+01;
  COFLAM[          34] =   0.5821675560201971E+00;
  COFLAM[          35] =  -0.2934671843536578E-01;
  COFLAM[          36] =  -0.1342414135164640E+02;
  COFLAM[          37] =   0.6083698271393525E+01;
  COFLAM[          38] =  -0.4766239458521935E+00;
  COFLAM[          39] =   0.1178953116268849E-01;
  COFLAM[          40] =  -0.9751824445950614E+01;
  COFLAM[          41] =   0.4413378685593287E+01;
  COFLAM[          42] =  -0.2468062806272802E+00;
  COFLAM[          43] =   0.1377122333476894E-02;
  COFLAM[          44] =   0.3562399283829174E+01;
  COFLAM[          45] =  -0.1581301349807987E+01;
  COFLAM[          46] =   0.6203935337751392E+00;
  COFLAM[          47] =  -0.4004666880047965E-01;
  COFLAM[          48] =  -0.1392568895481882E+02;
  COFLAM[          49] =   0.6085567690300339E+01;
  COFLAM[          50] =  -0.4671887217581250E+00;
  COFLAM[          51] =   0.1122436679253518E-01;
  COFLAM[          52] =  -0.6138049360593907E+01;
  COFLAM[          53] =   0.2471259577153958E+01;
  COFLAM[          54] =   0.6476679484655880E-01;
  COFLAM[          55] =  -0.1455106906367778E-01;
  COFLAM[          56] =  -0.1938382757224525E+01;
  COFLAM[          57] =   0.2891812367108134E+01;
  COFLAM[          58] =  -0.2711200250304408E+00;
  COFLAM[          59] =   0.1152522945287249E-01;
  COFLAM[          60] =  -0.1130096296226325E+01;
  COFLAM[          61] =   0.2340066562944334E+01;
  COFLAM[          62] =  -0.1632055932993225E+00;
  COFLAM[          63] =   0.5799980518407136E-02;
  COFLAM[          64] =   0.8981260119898151E+00;
  COFLAM[          65] =   0.1325300734490483E+01;
  COFLAM[          66] =   0.1818510064860679E-01;
  COFLAM[          67] =  -0.4467526168120878E-02;
  COFLAM[          68] =  -0.1154097849909565E+02;
  COFLAM[          69] =   0.5968337094767509E+01;
  COFLAM[          70] =  -0.5826401624646567E+00;
  COFLAM[          71] =   0.2109963828235269E-01;
  COFLAM[          72] =  -0.9884154801418376E+01;
  COFLAM[          73] =   0.4139058806234985E+01;
  COFLAM[          74] =  -0.1749630791503862E+00;
  COFLAM[          75] =  -0.3280828993011904E-02;
  COFLAM[          76] =  -0.1364058777845794E+02;
  COFLAM[          77] =   0.5702536702255589E+01;
  COFLAM[          78] =  -0.4144369477580017E+00;
  COFLAM[          79] =   0.8355378454356676E-02;
  COFLAM[          80] =  -0.8307846828394501E+01;
  COFLAM[          81] =   0.3455145940791863E+01;
  COFLAM[          82] =  -0.7676228923521952E-01;
  COFLAM[          83] =  -0.8016848317299876E-02;
  COFLAM[          84] =  -0.1792571228314050E+02;
  COFLAM[          85] =   0.8090046753015994E+01;
  COFLAM[          86] =  -0.8013688267705129E+00;
  COFLAM[          87] =   0.2806493657228643E-01;
  COFLAM[          88] =  -0.8639135073712978E+01;
  COFLAM[          89] =   0.3577784819715171E+01;
  COFLAM[          90] =  -0.1011234435899089E+00;
  COFLAM[          91] =  -0.6856167407634219E-02;
  COFLAM[          92] =  -0.8201973292350056E+01;
  COFLAM[          93] =   0.3614059165453200E+01;
  COFLAM[          94] =  -0.1246288473368713E+00;
  COFLAM[          95] =  -0.4811623494744587E-02;
  COFLAM[          96] =  -0.9028894095695616E+01;
  COFLAM[          97] =   0.4130186445296287E+01;
  COFLAM[          98] =  -0.2208907562042679E+00;
  COFLAM[          99] =   0.1809715093822095E-03;
  COFLAM[         100] =  -0.4747765397850146E+01;
  COFLAM[         101] =   0.2392750481960765E+01;
  COFLAM[         102] =   0.1131073948410975E-01;
  COFLAM[         103] =  -0.1010837059251814E-01;
  COFLAM[         104] =  -0.1494705931999075E+02;
  COFLAM[         105] =   0.6751990761547317E+01;
  COFLAM[         106] =  -0.6002258006930086E+00;
  COFLAM[         107] =   0.1859844416221185E-01;
  COFLAM[         108] =  -0.1839997799853145E+02;
  COFLAM[         109] =   0.8432845949842360E+01;
  COFLAM[         110] =  -0.8585690252027328E+00;
  COFLAM[         111] =   0.3119497600749326E-01;
  COFLAM[         112] =  -0.2249644940665438E+02;
  COFLAM[         113] =   0.1027211233031001E+02;
  COFLAM[         114] =  -0.1126822567074276E+01;
  COFLAM[         115] =   0.4418820857328377E-01;
  COFLAM[         116] =   0.1293008947455327E+02;
  COFLAM[         117] =  -0.3528394988071363E+01;
  COFLAM[         118] =   0.6455858113801279E+00;
  COFLAM[         119] =  -0.3194427362954873E-01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFETA EGTRANSETCOFETA
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFETA egtransetcofeta
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFETA egtransetcofeta_
#endif
void egtransetCOFETA(double* COFETA) {
  COFETA[           0] =  -0.2040534341454124E+02;
  COFETA[           1] =   0.3652691486295009E+01;
  COFETA[           2] =  -0.3980303020608006E+00;
  COFETA[           3] =   0.1757072886401361E-01;
  COFETA[           4] =  -0.1384035451131115E+02;
  COFETA[           5] =   0.1003491325708292E+01;
  COFETA[           6] =  -0.5016044554561987E-01;
  COFETA[           7] =   0.2330995224065489E-02;
  COFETA[           8] =  -0.2022833518474935E+02;
  COFETA[           9] =   0.3630408718312077E+01;
  COFETA[          10] =  -0.3952256976799788E+00;
  COFETA[          11] =   0.1745288966941275E-01;
  COFETA[          12] =  -0.1510027705857393E+02;
  COFETA[          13] =   0.1929024677614498E+01;
  COFETA[          14] =  -0.1738657445088048E+00;
  COFETA[          15] =   0.7841476915125451E-02;
  COFETA[          16] =  -0.2000457400541702E+02;
  COFETA[          17] =   0.3569542093003574E+01;
  COFETA[          18] =  -0.3874920392785122E+00;
  COFETA[          19] =   0.1712461411247914E-01;
  COFETA[          20] =  -0.1506972928055982E+02;
  COFETA[          21] =   0.1929024677614455E+01;
  COFETA[          22] =  -0.1738657445087985E+00;
  COFETA[          23] =   0.7841476915125141E-02;
  COFETA[          24] =  -0.1055754979620539E+02;
  COFETA[          25] =  -0.1377850378695579E+01;
  COFETA[          26] =   0.4213981638352425E+00;
  COFETA[          27] =  -0.2414423055966423E-01;
  COFETA[          28] =  -0.2474505588723548E+02;
  COFETA[          29] =   0.5290034654913134E+01;
  COFETA[          30] =  -0.5878529487588758E+00;
  COFETA[          31] =   0.2486206312507097E-01;
  COFETA[          32] =  -0.1661561262594180E+02;
  COFETA[          33] =   0.2400975158113569E+01;
  COFETA[          34] =  -0.2357717790312281E+00;
  COFETA[          35] =   0.1054820948438182E-01;
  COFETA[          36] =  -0.2393711445663787E+02;
  COFETA[          37] =   0.5100701997082546E+01;
  COFETA[          38] =  -0.5700811865820409E+00;
  COFETA[          39] =   0.2436923618342013E-01;
  COFETA[          40] =  -0.2452950079636804E+02;
  COFETA[          41] =   0.5145720188115567E+01;
  COFETA[          42] =  -0.5734515028626993E+00;
  COFETA[          43] =   0.2441369795441310E-01;
  COFETA[          44] =  -0.1985295977766010E+02;
  COFETA[          45] =   0.2703811616927011E+01;
  COFETA[          46] =  -0.1672355189641908E+00;
  COFETA[          47] =   0.3212257118449540E-02;
  COFETA[          48] =  -0.2451245312341208E+02;
  COFETA[          49] =   0.5145720188115497E+01;
  COFETA[          50] =  -0.5734515028626889E+00;
  COFETA[          51] =   0.2441369795441260E-01;
  COFETA[          52] =  -0.1997384799268116E+02;
  COFETA[          53] =   0.2861750522819142E+01;
  COFETA[          54] =  -0.2024899143538781E+00;
  COFETA[          55] =   0.5362464436156140E-02;
  COFETA[          56] =  -0.1715809053469514E+02;
  COFETA[          57] =   0.2678088349030608E+01;
  COFETA[          58] =  -0.2721592407921913E+00;
  COFETA[          59] =   0.1214173232611473E-01;
  COFETA[          60] =  -0.1714258339027710E+02;
  COFETA[          61] =   0.2678088349030630E+01;
  COFETA[          62] =  -0.2721592407921950E+00;
  COFETA[          63] =   0.1214173232611492E-01;
  COFETA[          64] =  -0.1712754275668106E+02;
  COFETA[          65] =   0.2678088349030581E+01;
  COFETA[          66] =  -0.2721592407921878E+00;
  COFETA[          67] =   0.1214173232611457E-01;
  COFETA[          68] =  -0.2397057295682966E+02;
  COFETA[          69] =   0.5130426196036950E+01;
  COFETA[          70] =  -0.5724284704186094E+00;
  COFETA[          71] =   0.2440888721969576E-01;
  COFETA[          72] =  -0.2231997819119366E+02;
  COFETA[          73] =   0.3866553929075860E+01;
  COFETA[          74] =  -0.3419581218331932E+00;
  COFETA[          75] =   0.1173156132419994E-01;
  COFETA[          76] =  -0.2125098599057449E+02;
  COFETA[          77] =   0.3262696346070605E+01;
  COFETA[          78] =  -0.2503683629086272E+00;
  COFETA[          79] =   0.7235695952542563E-02;
  COFETA[          80] =  -0.2502945492144446E+02;
  COFETA[          81] =   0.5159446447609293E+01;
  COFETA[          82] =  -0.5500619928927599E+00;
  COFETA[          83] =   0.2237585130462668E-01;
  COFETA[          84] =  -0.2348004890799276E+02;
  COFETA[          85] =   0.4341458455321099E+01;
  COFETA[          86] =  -0.4151006205826759E+00;
  COFETA[          87] =   0.1536292773840759E-01;
  COFETA[          88] =  -0.2347158533923254E+02;
  COFETA[          89] =   0.4341458455321031E+01;
  COFETA[          90] =  -0.4151006205826658E+00;
  COFETA[          91] =   0.1536292773840708E-01;
  COFETA[          92] =  -0.2488040801201845E+02;
  COFETA[          93] =   0.5159446447609339E+01;
  COFETA[          94] =  -0.5500619928927664E+00;
  COFETA[          95] =   0.2237585130462699E-01;
  COFETA[          96] =  -0.2478549066625603E+02;
  COFETA[          97] =   0.5159446447609294E+01;
  COFETA[          98] =  -0.5500619928927601E+00;
  COFETA[          99] =   0.2237585130462669E-01;
  COFETA[         100] =  -0.2478549066625603E+02;
  COFETA[         101] =   0.5159446447609294E+01;
  COFETA[         102] =  -0.5500619928927601E+00;
  COFETA[         103] =   0.2237585130462669E-01;
  COFETA[         104] =  -0.2477223626202213E+02;
  COFETA[         105] =   0.5159446447609334E+01;
  COFETA[         106] =  -0.5500619928927658E+00;
  COFETA[         107] =   0.2237585130462696E-01;
  COFETA[         108] =  -0.2468335508035342E+02;
  COFETA[         109] =   0.5159446447609318E+01;
  COFETA[         110] =  -0.5500619928927633E+00;
  COFETA[         111] =   0.2237585130462684E-01;
  COFETA[         112] =  -0.2459858550543459E+02;
  COFETA[         113] =   0.5159446447609233E+01;
  COFETA[         114] =  -0.5500619928927513E+00;
  COFETA[         115] =   0.2237585130462627E-01;
  COFETA[         116] =  -0.1656563666406150E+02;
  COFETA[         117] =   0.2388167035581858E+01;
  COFETA[         118] =  -0.2341208182867099E+00;
  COFETA[         119] =   0.1047727172770704E-01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFD EGTRANSETCOFD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFD egtransetcofd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFD egtransetcofd_
#endif
void egtransetCOFD(double* COFD) {
  COFD[           0] =  -0.1519084593013382E+02;
  COFD[           1] =   0.4384306309333000E+01;
  COFD[           2] =  -0.3557550068895214E+00;
  COFD[           3] =   0.1548014912218549E-01;
  COFD[           4] =  -0.1201910564637032E+02;
  COFD[           5] =   0.3033139269467014E+01;
  COFD[           6] =  -0.1860150369684279E+00;
  COFD[           7] =   0.8361628520974778E-02;
  COFD[           8] =  -0.1801636950310948E+02;
  COFD[           9] =   0.5090088857460223E+01;
  COFD[          10] =  -0.4457552498826564E+00;
  COFD[          11] =   0.1929141580461785E-01;
  COFD[          12] =  -0.1531287529990816E+02;
  COFD[          13] =   0.4270888061293578E+01;
  COFD[          14] =  -0.3482328696611331E+00;
  COFD[          15] =   0.1545097630860032E-01;
  COFD[          16] =  -0.1790236652645725E+02;
  COFD[          17] =   0.5054241512059871E+01;
  COFD[          18] =  -0.4413163862361035E+00;
  COFD[          19] =   0.1910841720545356E-01;
  COFD[          20] =  -0.1533072430565730E+02;
  COFD[          21] =   0.4277567720643019E+01;
  COFD[          22] =  -0.3491407301790579E+00;
  COFD[          23] =   0.1549202933974097E-01;
  COFD[          24] =  -0.1651315089501800E+02;
  COFD[          25] =   0.4220701654162998E+01;
  COFD[          26] =  -0.2825947916575933E+00;
  COFD[          27] =   0.1018065117940972E-01;
  COFD[          28] =  -0.1699415309027476E+02;
  COFD[          29] =   0.4447105495943485E+01;
  COFD[          30] =  -0.3323553153514281E+00;
  COFD[          31] =   0.1310173916235249E-01;
  COFD[          32] =  -0.1719132602269712E+02;
  COFD[          33] =   0.4859367175853364E+01;
  COFD[          34] =  -0.4240580935572132E+00;
  COFD[          35] =   0.1870696403708276E-01;
  COFD[          36] =  -0.1718718424688124E+02;
  COFD[          37] =   0.4609957970156615E+01;
  COFD[          38] =  -0.3605879229749922E+00;
  COFD[          39] =   0.1462630202481900E-01;
  COFD[          40] =  -0.1739069348322545E+02;
  COFD[          41] =   0.4517486135284712E+01;
  COFD[          42] =  -0.3398734156166323E+00;
  COFD[          43] =   0.1334270297473222E-01;
  COFD[          44] =  -0.1539084008249366E+02;
  COFD[          45] =   0.3565006575777335E+01;
  COFD[          46] =  -0.1867580729432450E+00;
  COFD[          47] =   0.5536884066371989E-02;
  COFD[          48] =  -0.1738939651119507E+02;
  COFD[          49] =   0.4516024567959285E+01;
  COFD[          50] =  -0.3395710901746170E+00;
  COFD[          51] =   0.1332458515958147E-01;
  COFD[          52] =  -0.1561125740910447E+02;
  COFD[          53] =   0.3698144089263720E+01;
  COFD[          54] =  -0.2103080207621357E+00;
  COFD[          55] =   0.6792204226221277E-02;
  COFD[          56] =  -0.1761107854705398E+02;
  COFD[          57] =   0.5048254268555682E+01;
  COFD[          58] =  -0.4490719655111177E+00;
  COFD[          59] =   0.1981332179135272E-01;
  COFD[          60] =  -0.1762111333528786E+02;
  COFD[          61] =   0.5052180357519856E+01;
  COFD[          62] =  -0.4496001366665920E+00;
  COFD[          63] =   0.1983696581079775E-01;
  COFD[          64] =  -0.1763070827156596E+02;
  COFD[          65] =   0.5055937247475734E+01;
  COFD[          66] =  -0.4501055595390027E+00;
  COFD[          67] =   0.1985959199754691E-01;
  COFD[          68] =  -0.1725034548774479E+02;
  COFD[          69] =   0.4566921339972978E+01;
  COFD[          70] =  -0.3509192271095338E+00;
  COFD[          71] =   0.1402685086906688E-01;
  COFD[          72] =  -0.1508763665602342E+02;
  COFD[          73] =   0.3400141263568119E+01;
  COFD[          74] =  -0.1645508303328382E+00;
  COFD[          75] =   0.4517926211356998E-02;
  COFD[          76] =  -0.1434856733752978E+02;
  COFD[          77] =   0.2972259604459018E+01;
  COFD[          78] =  -0.9821869542826608E-01;
  COFD[          79] =   0.1198325317137719E-02;
  COFD[          80] =  -0.1569651926124499E+02;
  COFD[          81] =   0.3625575306654659E+01;
  COFD[          82] =  -0.1996382427440662E+00;
  COFD[          83] =   0.6236602581187985E-02;
  COFD[          84] =  -0.1433143401213129E+02;
  COFD[          85] =   0.2944305087428476E+01;
  COFD[          86] =  -0.9451478525573859E-01;
  COFD[          87] =   0.9946176585968385E-03;
  COFD[          88] =  -0.1432407507363763E+02;
  COFD[          89] =   0.2940787058428872E+01;
  COFD[          90] =  -0.9397416078466805E-01;
  COFD[          91] =   0.9675574576547308E-03;
  COFD[          92] =  -0.1561208690963572E+02;
  COFD[          93] =   0.3582789091852060E+01;
  COFD[          94] =  -0.1928143728969187E+00;
  COFD[          95] =   0.5885365651518750E-02;
  COFD[          96] =  -0.1556317061681839E+02;
  COFD[          97] =   0.3558114950801765E+01;
  COFD[          98] =  -0.1888799785819696E+00;
  COFD[          99] =   0.5682884908013199E-02;
  COFD[         100] =  -0.1556317061681839E+02;
  COFD[         101] =   0.3558114950801765E+01;
  COFD[         102] =  -0.1888799785819696E+00;
  COFD[         103] =   0.5682884908013199E-02;
  COFD[         104] =  -0.1555662334580070E+02;
  COFD[         105] =   0.3554818568932333E+01;
  COFD[         106] =  -0.1883544064628437E+00;
  COFD[         107] =   0.5655838508786585E-02;
  COFD[         108] =  -0.1551443723602216E+02;
  COFD[         109] =   0.3533613016857426E+01;
  COFD[         110] =  -0.1849737082823091E+00;
  COFD[         111] =   0.5481875490852130E-02;
  COFD[         112] =  -0.1547687771509115E+02;
  COFD[         113] =   0.3514781617610363E+01;
  COFD[         114] =  -0.1819719591701729E+00;
  COFD[         115] =   0.5327428541928384E-02;
  COFD[         116] =  -0.1711567017943538E+02;
  COFD[         117] =   0.4833264586465570E+01;
  COFD[         118] =  -0.4205908440542898E+00;
  COFD[         119] =   0.1855345683322859E-01;
  COFD[         120] =  -0.1201910564637032E+02;
  COFD[         121] =   0.3033139269467014E+01;
  COFD[         122] =  -0.1860150369684279E+00;
  COFD[         123] =   0.8361628520974778E-02;
  COFD[         124] =  -0.1031825736753964E+02;
  COFD[         125] =   0.2192412640005562E+01;
  COFD[         126] =  -0.7539036147594932E-01;
  COFD[         127] =   0.3511225190832653E-02;
  COFD[         128] =  -0.1358415095584097E+02;
  COFD[         129] =   0.3239845970896248E+01;
  COFD[         130] =  -0.2150500406377150E+00;
  COFD[         131] =   0.9715225271588263E-02;
  COFD[         132] =  -0.1093372460559999E+02;
  COFD[         133] =   0.2306133612968357E+01;
  COFD[         134] =  -0.8736427899565385E-01;
  COFD[         135] =   0.3897526119633928E-02;
  COFD[         136] =  -0.1353824314289343E+02;
  COFD[         137] =   0.3227352939332374E+01;
  COFD[         138] =  -0.2134245159235849E+00;
  COFD[         139] =   0.9644648995332943E-02;
  COFD[         140] =  -0.1093673058791694E+02;
  COFD[         141] =   0.2306118133758978E+01;
  COFD[         142] =  -0.8736238957367691E-01;
  COFD[         143] =   0.3897448280154413E-02;
  COFD[         144] =  -0.1818497883283899E+02;
  COFD[         145] =   0.5043624597310085E+01;
  COFD[         146] =  -0.4377066809210735E+00;
  COFD[         147] =   0.1887643866477002E-01;
  COFD[         148] =  -0.1586051006549287E+02;
  COFD[         149] =   0.4116349686729230E+01;
  COFD[         150] =  -0.3299412099375534E+00;
  COFD[         151] =   0.1473720935464579E-01;
  COFD[         152] =  -0.1268296956267843E+02;
  COFD[         153] =   0.2903428401599904E+01;
  COFD[         154] =  -0.1707667302851327E+00;
  COFD[         155] =   0.7771111802295727E-02;
  COFD[         156] =  -0.1521342233375189E+02;
  COFD[         157] =   0.3885324867997948E+01;
  COFD[         158] =  -0.2989244982570487E+00;
  COFD[         159] =   0.1334849796046345E-01;
  COFD[         160] =  -0.1609618892523237E+02;
  COFD[         161] =   0.4147159660004064E+01;
  COFD[         162] =  -0.3347661201338916E+00;
  COFD[         163] =   0.1498106836636175E-01;
  COFD[         164] =  -0.1747771970561251E+02;
  COFD[         165] =   0.4673159679696371E+01;
  COFD[         166] =  -0.3935109715079382E+00;
  COFD[         167] =   0.1711386921544412E-01;
  COFD[         168] =  -0.1610856624128618E+02;
  COFD[         169] =   0.4151906359557744E+01;
  COFD[         170] =  -0.3354185801849658E+00;
  COFD[         171] =   0.1501089101157305E-01;
  COFD[         172] =  -0.1797854335218850E+02;
  COFD[         173] =   0.4905020138101169E+01;
  COFD[         174] =  -0.4293084201751438E+00;
  COFD[         175] =   0.1891302999394448E-01;
  COFD[         176] =  -0.1270707915711888E+02;
  COFD[         177] =   0.2929619774221139E+01;
  COFD[         178] =  -0.1739670601708100E+00;
  COFD[         179] =   0.7901198812704345E-02;
  COFD[         180] =  -0.1271124360389947E+02;
  COFD[         181] =   0.2931047949505553E+01;
  COFD[         182] =  -0.1741679825526774E+00;
  COFD[         183] =   0.7910585398301936E-02;
  COFD[         184] =  -0.1271518733681371E+02;
  COFD[         185] =   0.2932402241140174E+01;
  COFD[         186] =  -0.1743585082612658E+00;
  COFD[         187] =   0.7919486191765343E-02;
  COFD[         188] =  -0.1564739122681956E+02;
  COFD[         189] =   0.4025825276486324E+01;
  COFD[         190] =  -0.3181195147607377E+00;
  COFD[         191] =   0.1422136468533644E-01;
  COFD[         192] =  -0.1770170994272074E+02;
  COFD[         193] =   0.4743686490643380E+01;
  COFD[         194] =  -0.4058992902065018E+00;
  COFD[         195] =   0.1778646492198370E-01;
  COFD[         196] =  -0.1856440584784442E+02;
  COFD[         197] =   0.5034796976697628E+01;
  COFD[         198] =  -0.4442640448549661E+00;
  COFD[         199] =   0.1947279036625017E-01;
  COFD[         200] =  -0.1779799510057660E+02;
  COFD[         201] =   0.4774820362539309E+01;
  COFD[         202] =  -0.4169337844744687E+00;
  COFD[         203] =   0.1856826481796923E-01;
  COFD[         204] =  -0.1866728292989077E+02;
  COFD[         205] =   0.5078383565507850E+01;
  COFD[         206] =  -0.4538464035977413E+00;
  COFD[         207] =   0.2005643240392637E-01;
  COFD[         208] =  -0.1867240731108900E+02;
  COFD[         209] =   0.5080367429995977E+01;
  COFD[         210] =  -0.4541125321914433E+00;
  COFD[         211] =   0.2006831547910542E-01;
  COFD[         212] =  -0.1789359789035752E+02;
  COFD[         213] =   0.4812277632492618E+01;
  COFD[         214] =  -0.4220338205719207E+00;
  COFD[         215] =   0.1879929285184377E-01;
  COFD[         216] =  -0.1794643219743921E+02;
  COFD[         217] =   0.4833075314023510E+01;
  COFD[         218] =  -0.4248631135846847E+00;
  COFD[         219] =   0.1892735769307228E-01;
  COFD[         220] =  -0.1794643219743921E+02;
  COFD[         221] =   0.4833075314023510E+01;
  COFD[         222] =  -0.4248631135846847E+00;
  COFD[         223] =   0.1892735769307228E-01;
  COFD[         224] =  -0.1795336719760172E+02;
  COFD[         225] =   0.4835810632619655E+01;
  COFD[         226] =  -0.4252350727772528E+00;
  COFD[         227] =   0.1894418782785674E-01;
  COFD[         228] =  -0.1799730159578534E+02;
  COFD[         229] =   0.4853169597207414E+01;
  COFD[         230] =  -0.4275947086398001E+00;
  COFD[         231] =   0.1905091769934767E-01;
  COFD[         232] =  -0.1803535719272616E+02;
  COFD[         233] =   0.4868249132378355E+01;
  COFD[         234] =  -0.4296431131249400E+00;
  COFD[         235] =   0.1914351307721480E-01;
  COFD[         236] =  -0.1266066516595179E+02;
  COFD[         237] =   0.2898076057146818E+01;
  COFD[         238] =  -0.1700453759318706E+00;
  COFD[         239] =   0.7738690296212197E-02;
  COFD[         240] =  -0.1801636950310948E+02;
  COFD[         241] =   0.5090088857460223E+01;
  COFD[         242] =  -0.4457552498826564E+00;
  COFD[         243] =   0.1929141580461785E-01;
  COFD[         244] =  -0.1358415095584097E+02;
  COFD[         245] =   0.3239845970896248E+01;
  COFD[         246] =  -0.2150500406377150E+00;
  COFD[         247] =   0.9715225271588263E-02;
  COFD[         248] =  -0.1772990792056323E+02;
  COFD[         249] =   0.4367699379901269E+01;
  COFD[         250] =  -0.3537277912983985E+00;
  COFD[         251] =   0.1539755051214736E-01;
  COFD[         252] =  -0.1561687122378518E+02;
  COFD[         253] =   0.3703886367486989E+01;
  COFD[         254] =  -0.2712410029515650E+00;
  COFD[         255] =   0.1197280722366742E-01;
  COFD[         256] =  -0.1766459414605863E+02;
  COFD[         257] =   0.4343242779501765E+01;
  COFD[         258] =  -0.3506769601296933E+00;
  COFD[         259] =   0.1527024825917116E-01;
  COFD[         260] =  -0.1562210625534589E+02;
  COFD[         261] =   0.3699988818188416E+01;
  COFD[         262] =  -0.2707081263160739E+00;
  COFD[         263] =   0.1194858184298051E-01;
  COFD[         264] =  -0.2042453573019913E+02;
  COFD[         265] =   0.5222435968901098E+01;
  COFD[         266] =  -0.4324987976330752E+00;
  COFD[         267] =   0.1751958312757915E-01;
  COFD[         268] =  -0.1938215107014084E+02;
  COFD[         269] =   0.4860644015548545E+01;
  COFD[         270] =  -0.4060689337372705E+00;
  COFD[         271] =   0.1718274802992044E-01;
  COFD[         272] =  -0.1661214673917502E+02;
  COFD[         273] =   0.3930381940303202E+01;
  COFD[         274] =  -0.2994310262293675E+00;
  COFD[         275] =   0.1314286546846328E-01;
  COFD[         276] =  -0.1906497740611485E+02;
  COFD[         277] =   0.4777590849182513E+01;
  COFD[         278] =  -0.3975472982609866E+00;
  COFD[         279] =   0.1690531426791657E-01;
  COFD[         280] =  -0.1938165462851885E+02;
  COFD[         281] =   0.4799455794296509E+01;
  COFD[         282] =  -0.3987831540937672E+00;
  COFD[         283] =   0.1688981189405884E-01;
  COFD[         284] =  -0.2045023591963550E+02;
  COFD[         285] =   0.5108641161759671E+01;
  COFD[         286] =  -0.4199387136777197E+00;
  COFD[         287] =   0.1704933199723114E-01;
  COFD[         288] =  -0.1938533704705245E+02;
  COFD[         289] =   0.4797779503867271E+01;
  COFD[         290] =  -0.3984320673303173E+00;
  COFD[         291] =   0.1686862693067165E-01;
  COFD[         292] =  -0.2012688925964889E+02;
  COFD[         293] =   0.5018247199767933E+01;
  COFD[         294] =  -0.4131414445091446E+00;
  COFD[         295] =   0.1694646037906594E-01;
  COFD[         296] =  -0.1707094794267795E+02;
  COFD[         297] =   0.4119810624341889E+01;
  COFD[         298] =  -0.3242844616145015E+00;
  COFD[         299] =   0.1423184749658152E-01;
  COFD[         300] =  -0.1709003104221693E+02;
  COFD[         301] =   0.4125627998285617E+01;
  COFD[         302] =  -0.3250670331906144E+00;
  COFD[         303] =   0.1426686739847163E-01;
  COFD[         304] =  -0.1710868501481421E+02;
  COFD[         305] =   0.4131370391002061E+01;
  COFD[         306] =  -0.3258395192522007E+00;
  COFD[         307] =   0.1430143624455670E-01;
  COFD[         308] =  -0.1923136095279614E+02;
  COFD[         309] =   0.4771652981887103E+01;
  COFD[         310] =  -0.3944968044729790E+00;
  COFD[         311] =   0.1667283377652949E-01;
  COFD[         312] =  -0.2002012999655281E+02;
  COFD[         313] =   0.4895179416250934E+01;
  COFD[         314] =  -0.3929898448063785E+00;
  COFD[         315] =   0.1589904718416698E-01;
  COFD[         316] =  -0.2017297081385243E+02;
  COFD[         317] =   0.4879127539439356E+01;
  COFD[         318] =  -0.3874272290357592E+00;
  COFD[         319] =   0.1551838396554711E-01;
  COFD[         320] =  -0.1979190800756452E+02;
  COFD[         321] =   0.4812635185197916E+01;
  COFD[         322] =  -0.3895191141694820E+00;
  COFD[         323] =   0.1603137873834851E-01;
  COFD[         324] =  -0.1984047911980011E+02;
  COFD[         325] =   0.4733570408392962E+01;
  COFD[         326] =  -0.3701539096262219E+00;
  COFD[         327] =   0.1481934816743204E-01;
  COFD[         328] =  -0.1983092143620168E+02;
  COFD[         329] =   0.4728281143697850E+01;
  COFD[         330] =  -0.3693402228450711E+00;
  COFD[         331] =   0.1477857564328594E-01;
  COFD[         332] =  -0.1969361932498264E+02;
  COFD[         333] =   0.4748312088352801E+01;
  COFD[         334] =  -0.3792560078874427E+00;
  COFD[         335] =   0.1550280308751674E-01;
  COFD[         336] =  -0.1963494149489594E+02;
  COFD[         337] =   0.4711113691258022E+01;
  COFD[         338] =  -0.3733233453642334E+00;
  COFD[         339] =   0.1519736886874179E-01;
  COFD[         340] =  -0.1963494149489594E+02;
  COFD[         341] =   0.4711113691258022E+01;
  COFD[         342] =  -0.3733233453642334E+00;
  COFD[         343] =   0.1519736886874179E-01;
  COFD[         344] =  -0.1962710518342632E+02;
  COFD[         345] =   0.4706196970979815E+01;
  COFD[         346] =  -0.3725393175856263E+00;
  COFD[         347] =   0.1515701000479143E-01;
  COFD[         348] =  -0.1957709787054268E+02;
  COFD[         349] =   0.4675065521077943E+01;
  COFD[         350] =  -0.3675756867404701E+00;
  COFD[         351] =   0.1490152940712903E-01;
  COFD[         352] =  -0.1953370018155996E+02;
  COFD[         353] =   0.4648360822131973E+01;
  COFD[         354] =  -0.3633186948567744E+00;
  COFD[         355] =   0.1468245773755453E-01;
  COFD[         356] =  -0.1658229952875101E+02;
  COFD[         357] =   0.3922184822536529E+01;
  COFD[         358] =  -0.2983959485115160E+00;
  COFD[         359] =   0.1309927190981370E-01;
  COFD[         360] =  -0.1531287529990816E+02;
  COFD[         361] =   0.4270888061293578E+01;
  COFD[         362] =  -0.3482328696611331E+00;
  COFD[         363] =   0.1545097630860032E-01;
  COFD[         364] =  -0.1093372460559999E+02;
  COFD[         365] =   0.2306133612968357E+01;
  COFD[         366] =  -0.8736427899565385E-01;
  COFD[         367] =   0.3897526119633928E-02;
  COFD[         368] =  -0.1561687122378518E+02;
  COFD[         369] =   0.3703886367486989E+01;
  COFD[         370] =  -0.2712410029515650E+00;
  COFD[         371] =   0.1197280722366742E-01;
  COFD[         372] =  -0.1356847565002614E+02;
  COFD[         373] =   0.3062175626208075E+01;
  COFD[         374] =  -0.1889691489381823E+00;
  COFD[         375] =   0.8454139853458284E-02;
  COFD[         376] =  -0.1556217965348612E+02;
  COFD[         377] =   0.3683628994574988E+01;
  COFD[         378] =  -0.2686556693855442E+00;
  COFD[         379] =   0.1186253607235680E-01;
  COFD[         380] =  -0.1358466752355666E+02;
  COFD[         381] =   0.3062671492876047E+01;
  COFD[         382] =  -0.1890385291581307E+00;
  COFD[         383] =   0.8457362732499073E-02;
  COFD[         384] =  -0.1906508600184398E+02;
  COFD[         385] =   0.4989055420187643E+01;
  COFD[         386] =  -0.4192832417141421E+00;
  COFD[         387] =   0.1761870656062123E-01;
  COFD[         388] =  -0.1784670232120704E+02;
  COFD[         389] =   0.4483333996808572E+01;
  COFD[         390] =  -0.3684301675626975E+00;
  COFD[         391] =   0.1601948073360996E-01;
  COFD[         392] =  -0.1478198899238034E+02;
  COFD[         393] =   0.3371092581023257E+01;
  COFD[         394] =  -0.2298567856863946E+00;
  COFD[         395] =   0.1025846107153025E-01;
  COFD[         396] =  -0.1741286728352899E+02;
  COFD[         397] =   0.4346824031835780E+01;
  COFD[         398] =  -0.3515350630281370E+00;
  COFD[         399] =   0.1532053929906590E-01;
  COFD[         400] =  -0.1797996041020607E+02;
  COFD[         401] =   0.4462878081531843E+01;
  COFD[         402] =  -0.3662487521791470E+00;
  COFD[         403] =   0.1594206335130083E-01;
  COFD[         404] =  -0.1914986770399175E+02;
  COFD[         405] =   0.4864899621234062E+01;
  COFD[         406] =  -0.4039220297382060E+00;
  COFD[         407] =   0.1697258735720344E-01;
  COFD[         408] =  -0.1800786212441529E+02;
  COFD[         409] =   0.4471520139789249E+01;
  COFD[         410] =  -0.3673588813728833E+00;
  COFD[         411] =   0.1598943652325460E-01;
  COFD[         412] =  -0.1878355293596396E+02;
  COFD[         413] =   0.4745970603113269E+01;
  COFD[         414] =  -0.3924372163577848E+00;
  COFD[         415] =   0.1663659729501466E-01;
  COFD[         416] =  -0.1511165027002312E+02;
  COFD[         417] =   0.3515441252015493E+01;
  COFD[         418] =  -0.2491198761221894E+00;
  COFD[         419] =   0.1111625360314578E-01;
  COFD[         420] =  -0.1513524421958843E+02;
  COFD[         421] =   0.3523297948981071E+01;
  COFD[         422] =  -0.2502118585975575E+00;
  COFD[         423] =   0.1116667641067510E-01;
  COFD[         424] =  -0.1515824683455671E+02;
  COFD[         425] =   0.3531011422568503E+01;
  COFD[         426] =  -0.2512839088019835E+00;
  COFD[         427] =   0.1121617798528330E-01;
  COFD[         428] =  -0.1793673861488134E+02;
  COFD[         429] =   0.4489875914439151E+01;
  COFD[         430] =  -0.3697556337020395E+00;
  COFD[         431] =   0.1609253879718407E-01;
  COFD[         432] =  -0.1900626883865293E+02;
  COFD[         433] =   0.4754407177756363E+01;
  COFD[         434] =  -0.3901998067876299E+00;
  COFD[         435] =   0.1639341623231018E-01;
  COFD[         436] =  -0.1917251019075018E+02;
  COFD[         437] =   0.4736217002347487E+01;
  COFD[         438] =  -0.3839903059289852E+00;
  COFD[         439] =   0.1596104139772183E-01;
  COFD[         440] =  -0.1891691641938918E+02;
  COFD[         441] =   0.4708474963146238E+01;
  COFD[         442] =  -0.3908198429896564E+00;
  COFD[         443] =   0.1669173693263477E-01;
  COFD[         444] =  -0.1935099040011764E+02;
  COFD[         445] =   0.4796859780440744E+01;
  COFD[         446] =  -0.3951546730488872E+00;
  COFD[         447] =   0.1657370960389658E-01;
  COFD[         448] =  -0.1935470245750297E+02;
  COFD[         449] =   0.4797074097380636E+01;
  COFD[         450] =  -0.3950970944572844E+00;
  COFD[         451] =   0.1656735836140260E-01;
  COFD[         452] =  -0.1905867186145966E+02;
  COFD[         453] =   0.4745150284504969E+01;
  COFD[         454] =  -0.3946441379637312E+00;
  COFD[         455] =   0.1681441351814178E-01;
  COFD[         456] =  -0.1913891562399630E+02;
  COFD[         457] =   0.4766338982641832E+01;
  COFD[         458] =  -0.3968432892885747E+00;
  COFD[         459] =   0.1688440639323740E-01;
  COFD[         460] =  -0.1913891562399630E+02;
  COFD[         461] =   0.4766338982641832E+01;
  COFD[         462] =  -0.3968432892885747E+00;
  COFD[         463] =   0.1688440639323740E-01;
  COFD[         464] =  -0.1914947266287874E+02;
  COFD[         465] =   0.4769147094975692E+01;
  COFD[         466] =  -0.3971341787719369E+00;
  COFD[         467] =   0.1689363434236904E-01;
  COFD[         468] =  -0.1921620409589607E+02;
  COFD[         469] =   0.4787003103790272E+01;
  COFD[         470] =  -0.3989806507176665E+00;
  COFD[         471] =   0.1695203644218339E-01;
  COFD[         472] =  -0.1927347440067483E+02;
  COFD[         473] =   0.4802472345441742E+01;
  COFD[         474] =  -0.4005754779946568E+00;
  COFD[         475] =   0.1700221644405138E-01;
  COFD[         476] =  -0.1474819919893963E+02;
  COFD[         477] =   0.3361311502126538E+01;
  COFD[         478] =  -0.2285465363406641E+00;
  COFD[         479] =   0.1019990809984005E-01;
  COFD[         480] =  -0.1790236652645725E+02;
  COFD[         481] =   0.5054241512059871E+01;
  COFD[         482] =  -0.4413163862361035E+00;
  COFD[         483] =   0.1910841720545356E-01;
  COFD[         484] =  -0.1353824314289343E+02;
  COFD[         485] =   0.3227352939332374E+01;
  COFD[         486] =  -0.2134245159235849E+00;
  COFD[         487] =   0.9644648995332943E-02;
  COFD[         488] =  -0.1766459414605863E+02;
  COFD[         489] =   0.4343242779501765E+01;
  COFD[         490] =  -0.3506769601296933E+00;
  COFD[         491] =   0.1527024825917116E-01;
  COFD[         492] =  -0.1556217965348612E+02;
  COFD[         493] =   0.3683628994574988E+01;
  COFD[         494] =  -0.2686556693855442E+00;
  COFD[         495] =   0.1186253607235680E-01;
  COFD[         496] =  -0.1760280207090635E+02;
  COFD[         497] =   0.4320037087778842E+01;
  COFD[         498] =  -0.3477962090917468E+00;
  COFD[         499] =   0.1515062790890059E-01;
  COFD[         500] =  -0.1556532980903788E+02;
  COFD[         501] =   0.3678661452138608E+01;
  COFD[         502] =  -0.2679777825820168E+00;
  COFD[         503] =   0.1183177506717497E-01;
  COFD[         504] =  -0.2052802769939090E+02;
  COFD[         505] =   0.5190809253464677E+01;
  COFD[         506] =  -0.4197156907449909E+00;
  COFD[         507] =   0.1661902490158712E-01;
  COFD[         508] =  -0.1933923251029205E+02;
  COFD[         509] =   0.4846414430750782E+01;
  COFD[         510] =  -0.4047681110056587E+00;
  COFD[         511] =   0.1714937694388643E-01;
  COFD[         512] =  -0.1654504240089281E+02;
  COFD[         513] =   0.3902990496709335E+01;
  COFD[         514] =  -0.2959726841289695E+00;
  COFD[         515] =   0.1299731318228341E-01;
  COFD[         516] =  -0.1902332581359150E+02;
  COFD[         517] =   0.4763347082494895E+01;
  COFD[         518] =  -0.3961877942189693E+00;
  COFD[         519] =   0.1686706463721130E-01;
  COFD[         520] =  -0.1933548765861729E+02;
  COFD[         521] =   0.4783226123251535E+01;
  COFD[         522] =  -0.3972356460765860E+00;
  COFD[         523] =   0.1684656833081372E-01;
  COFD[         524] =  -0.2047278277859408E+02;
  COFD[         525] =   0.5123414930187060E+01;
  COFD[         526] =  -0.4228706903648477E+00;
  COFD[         527] =   0.1721819830097161E-01;
  COFD[         528] =  -0.1933986861083130E+02;
  COFD[         529] =   0.4781769764498772E+01;
  COFD[         530] =  -0.3969195763875204E+00;
  COFD[         531] =   0.1682718826455844E-01;
  COFD[         532] =  -0.2034321636694622E+02;
  COFD[         533] =   0.5089562300279502E+01;
  COFD[         534] =  -0.4211891334428361E+00;
  COFD[         535] =   0.1724909778947470E-01;
  COFD[         536] =  -0.1696131964138645E+02;
  COFD[         537] =   0.4074167708212651E+01;
  COFD[         538] =  -0.3182576642303228E+00;
  COFD[         539] =   0.1396612551510442E-01;
  COFD[         540] =  -0.1698006191382479E+02;
  COFD[         541] =   0.4079755891690069E+01;
  COFD[         542] =  -0.3190093940831369E+00;
  COFD[         543] =   0.1399976445173436E-01;
  COFD[         544] =  -0.1699844377766270E+02;
  COFD[         545] =   0.4085300971938842E+01;
  COFD[         546] =  -0.3197553288341113E+00;
  COFD[         547] =   0.1403314440725374E-01;
  COFD[         548] =  -0.1920918070133062E+02;
  COFD[         549] =   0.4764713252836613E+01;
  COFD[         550] =  -0.3942405534249802E+00;
  COFD[         551] =   0.1668923743050566E-01;
  COFD[         552] =  -0.2005721601235365E+02;
  COFD[         553] =   0.4914479851418061E+01;
  COFD[         554] =  -0.3965430355140246E+00;
  COFD[         555] =   0.1609667930913030E-01;
  COFD[         556] =  -0.2013382472396569E+02;
  COFD[         557] =   0.4863179495148787E+01;
  COFD[         558] =  -0.3856326932017908E+00;
  COFD[         559] =   0.1544705776346803E-01;
  COFD[         560] =  -0.1978759220318865E+02;
  COFD[         561] =   0.4813522654682605E+01;
  COFD[         562] =  -0.3904695959382222E+00;
  COFD[         563] =   0.1610776277485814E-01;
  COFD[         564] =  -0.1988931959990085E+02;
  COFD[         565] =   0.4756594033154599E+01;
  COFD[         566] =  -0.3742321497803292E+00;
  COFD[         567] =   0.1504181447689004E-01;
  COFD[         568] =  -0.1987978855175104E+02;
  COFD[         569] =   0.4751262797957020E+01;
  COFD[         570] =  -0.3734104385798870E+00;
  COFD[         571] =   0.1500058297946865E-01;
  COFD[         572] =  -0.1969237276879183E+02;
  COFD[         573] =   0.4749599566916338E+01;
  COFD[         574] =  -0.3802354707820949E+00;
  COFD[         575] =   0.1557944508694219E-01;
  COFD[         576] =  -0.1963449815268955E+02;
  COFD[         577] =   0.4712178662790742E+01;
  COFD[         578] =  -0.3742472943897919E+00;
  COFD[         579] =   0.1527044453467252E-01;
  COFD[         580] =  -0.1963449815268955E+02;
  COFD[         581] =   0.4712178662790742E+01;
  COFD[         582] =  -0.3742472943897919E+00;
  COFD[         583] =   0.1527044453467252E-01;
  COFD[         584] =  -0.1962672344426228E+02;
  COFD[         585] =   0.4707211623266654E+01;
  COFD[         586] =  -0.3734526020103864E+00;
  COFD[         587] =   0.1522944338733201E-01;
  COFD[         588] =  -0.1957689044936644E+02;
  COFD[         589] =   0.4675660284789568E+01;
  COFD[         590] =  -0.3684053408172555E+00;
  COFD[         591] =   0.1496906925902828E-01;
  COFD[         592] =  -0.1953337113047812E+02;
  COFD[         593] =   0.4648468363013314E+01;
  COFD[         594] =  -0.3640564212695835E+00;
  COFD[         595] =   0.1474476403759242E-01;
  COFD[         596] =  -0.1651845590523054E+02;
  COFD[         597] =   0.3896132878391479E+01;
  COFD[         598] =  -0.2951170694006329E+00;
  COFD[         599] =   0.1296170338310411E-01;
  COFD[         600] =  -0.1533072430565730E+02;
  COFD[         601] =   0.4277567720643019E+01;
  COFD[         602] =  -0.3491407301790579E+00;
  COFD[         603] =   0.1549202933974097E-01;
  COFD[         604] =  -0.1093673058791694E+02;
  COFD[         605] =   0.2306118133758978E+01;
  COFD[         606] =  -0.8736238957367691E-01;
  COFD[         607] =   0.3897448280154413E-02;
  COFD[         608] =  -0.1562210625534589E+02;
  COFD[         609] =   0.3699988818188416E+01;
  COFD[         610] =  -0.2707081263160739E+00;
  COFD[         611] =   0.1194858184298051E-01;
  COFD[         612] =  -0.1358466752355666E+02;
  COFD[         613] =   0.3062671492876047E+01;
  COFD[         614] =  -0.1890385291581307E+00;
  COFD[         615] =   0.8457362732499073E-02;
  COFD[         616] =  -0.1556532980903788E+02;
  COFD[         617] =   0.3678661452138608E+01;
  COFD[         618] =  -0.2679777825820168E+00;
  COFD[         619] =   0.1183177506717497E-01;
  COFD[         620] =  -0.1359902342804023E+02;
  COFD[         621] =   0.3062175626208111E+01;
  COFD[         622] =  -0.1889691489381874E+00;
  COFD[         623] =   0.8454139853458534E-02;
  COFD[         624] =  -0.1908171031786331E+02;
  COFD[         625] =   0.4989633286602976E+01;
  COFD[         626] =  -0.4194071192530971E+00;
  COFD[         627] =   0.1762632910887997E-01;
  COFD[         628] =  -0.1783539626494398E+02;
  COFD[         629] =   0.4471504634175316E+01;
  COFD[         630] =  -0.3669099771101654E+00;
  COFD[         631] =   0.1595459643293670E-01;
  COFD[         632] =  -0.1476985219999908E+02;
  COFD[         633] =   0.3357751196108267E+01;
  COFD[         634] =  -0.2280046831004087E+00;
  COFD[         635] =   0.1017303667696902E-01;
  COFD[         636] =  -0.1740489593482223E+02;
  COFD[         637] =   0.4336092039978253E+01;
  COFD[         638] =  -0.3501568532068898E+00;
  COFD[         639] =   0.1526174747325281E-01;
  COFD[         640] =  -0.1796101043339682E+02;
  COFD[         641] =   0.4447717004791725E+01;
  COFD[         642] =  -0.3643006119047739E+00;
  COFD[         643] =   0.1585890337276808E-01;
  COFD[         644] =  -0.1918177070837826E+02;
  COFD[         645] =   0.4871873199573439E+01;
  COFD[         646] =  -0.4051602536432213E+00;
  COFD[         647] =   0.1704106539271711E-01;
  COFD[         648] =  -0.1798854490941215E+02;
  COFD[         649] =   0.4456122648699593E+01;
  COFD[         650] =  -0.3653808041282116E+00;
  COFD[         651] =   0.1590501766674265E-01;
  COFD[         652] =  -0.1880409131690931E+02;
  COFD[         653] =   0.4747896766189557E+01;
  COFD[         654] =  -0.3929513558466455E+00;
  COFD[         655] =   0.1667071927589917E-01;
  COFD[         656] =  -0.1509620843583920E+02;
  COFD[         657] =   0.3500304707531062E+01;
  COFD[         658] =  -0.2470160084537919E+00;
  COFD[         659] =   0.1101910392929726E-01;
  COFD[         660] =  -0.1511946285859442E+02;
  COFD[         661] =   0.3507927470352835E+01;
  COFD[         662] =  -0.2480755290441053E+00;
  COFD[         663] =   0.1106802953504252E-01;
  COFD[         664] =  -0.1514219804803715E+02;
  COFD[         665] =   0.3515441252015498E+01;
  COFD[         666] =  -0.2491198761221899E+00;
  COFD[         667] =   0.1111625360314580E-01;
  COFD[         668] =  -0.1792454027766008E+02;
  COFD[         669] =   0.4476326984714260E+01;
  COFD[         670] =  -0.3680191278135618E+00;
  COFD[         671] =   0.1601860672433683E-01;
  COFD[         672] =  -0.1903172737849495E+02;
  COFD[         673] =   0.4757798023561699E+01;
  COFD[         674] =  -0.3909686059341044E+00;
  COFD[         675] =   0.1644144176618614E-01;
  COFD[         676] =  -0.1920854681460308E+02;
  COFD[         677] =   0.4744336922244572E+01;
  COFD[         678] =  -0.3854716437444286E+00;
  COFD[         679] =   0.1604421161000185E-01;
  COFD[         680] =  -0.1891630498772998E+02;
  COFD[         681] =   0.4700565744006908E+01;
  COFD[         682] =  -0.3899919372649346E+00;
  COFD[         683] =   0.1666500811895143E-01;
  COFD[         684] =  -0.1936781090512805E+02;
  COFD[         685] =   0.4796060596924874E+01;
  COFD[         686] =  -0.3953630801104229E+00;
  COFD[         687] =   0.1659683769348035E-01;
  COFD[         688] =  -0.1937166591561948E+02;
  COFD[         689] =   0.4796286056188764E+01;
  COFD[         690] =  -0.3953052681051102E+00;
  COFD[         691] =   0.1659039955935203E-01;
  COFD[         692] =  -0.1906167945299586E+02;
  COFD[         693] =   0.4737944806110550E+01;
  COFD[         694] =  -0.3938946014462452E+00;
  COFD[         695] =   0.1679046692522489E-01;
  COFD[         696] =  -0.1914452665813124E+02;
  COFD[         697] =   0.4759724497948794E+01;
  COFD[         698] =  -0.3961575775878822E+00;
  COFD[         699] =   0.1686262510597603E-01;
  COFD[         700] =  -0.1914452665813124E+02;
  COFD[         701] =   0.4759724497948794E+01;
  COFD[         702] =  -0.3961575775878822E+00;
  COFD[         703] =   0.1686262510597603E-01;
  COFD[         704] =  -0.1915545032909366E+02;
  COFD[         705] =   0.4762618776775587E+01;
  COFD[         706] =  -0.3964577122504647E+00;
  COFD[         707] =   0.1687216358936520E-01;
  COFD[         708] =  -0.1922460653113258E+02;
  COFD[         709] =   0.4781057331709897E+01;
  COFD[         710] =  -0.3983664361805447E+00;
  COFD[         711] =   0.1693264372808513E-01;
  COFD[         712] =  -0.1928407792404587E+02;
  COFD[         713] =   0.4797069911008538E+01;
  COFD[         714] =  -0.4000190517915601E+00;
  COFD[         715] =   0.1698473884996324E-01;
  COFD[         716] =  -0.1473661449825940E+02;
  COFD[         717] =   0.3348204558833826E+01;
  COFD[         718] =  -0.2267271723233180E+00;
  COFD[         719] =   0.1011600240359858E-01;
  COFD[         720] =  -0.1651315089501800E+02;
  COFD[         721] =   0.4220701654162998E+01;
  COFD[         722] =  -0.2825947916575933E+00;
  COFD[         723] =   0.1018065117940972E-01;
  COFD[         724] =  -0.1818497883283899E+02;
  COFD[         725] =   0.5043624597310085E+01;
  COFD[         726] =  -0.4377066809210735E+00;
  COFD[         727] =   0.1887643866477002E-01;
  COFD[         728] =  -0.2042453573019913E+02;
  COFD[         729] =   0.5222435968901098E+01;
  COFD[         730] =  -0.4324987976330752E+00;
  COFD[         731] =   0.1751958312757915E-01;
  COFD[         732] =  -0.1906508600184398E+02;
  COFD[         733] =   0.4989055420187643E+01;
  COFD[         734] =  -0.4192832417141421E+00;
  COFD[         735] =   0.1761870656062123E-01;
  COFD[         736] =  -0.2052802769939090E+02;
  COFD[         737] =   0.5190809253464677E+01;
  COFD[         738] =  -0.4197156907449909E+00;
  COFD[         739] =   0.1661902490158712E-01;
  COFD[         740] =  -0.1908171031786331E+02;
  COFD[         741] =   0.4989633286602976E+01;
  COFD[         742] =  -0.4194071192530971E+00;
  COFD[         743] =   0.1762632910887997E-01;
  COFD[         744] =  -0.1180242969654581E+02;
  COFD[         745] =   0.8887910786873973E+00;
  COFD[         746] =   0.2461702200576275E+00;
  COFD[         747] =  -0.1607058638443385E-01;
  COFD[         748] =  -0.2000931867510009E+02;
  COFD[         749] =   0.4769188096212780E+01;
  COFD[         750] =  -0.3467567852966127E+00;
  COFD[         751] =   0.1276719834602697E-01;
  COFD[         752] =  -0.2031234004461199E+02;
  COFD[         753] =   0.5186312454417227E+01;
  COFD[         754] =  -0.4313742876779729E+00;
  COFD[         755] =   0.1759974844922337E-01;
  COFD[         756] =  -0.2016248856160786E+02;
  COFD[         757] =   0.4899874571765012E+01;
  COFD[         758] =  -0.3689624059605574E+00;
  COFD[         759] =   0.1393109509915711E-01;
  COFD[         760] =  -0.2002167887304272E+02;
  COFD[         761] =   0.4706367098965421E+01;
  COFD[         762] =  -0.3398295196679849E+00;
  COFD[         763] =   0.1251278469477959E-01;
  COFD[         764] =  -0.1757470287826831E+02;
  COFD[         765] =   0.3509436887842453E+01;
  COFD[         766] =  -0.1496102927025162E+00;
  COFD[         767] =   0.2961974572337444E-02;
  COFD[         768] =  -0.1999142375296304E+02;
  COFD[         769] =   0.4690231520616589E+01;
  COFD[         770] =  -0.3375026912192592E+00;
  COFD[         771] =   0.1240317371013226E-01;
  COFD[         772] =  -0.1619587514529141E+02;
  COFD[         773] =   0.2843235032975332E+01;
  COFD[         774] =  -0.5104314541884596E-01;
  COFD[         775] =  -0.1707937587486378E-02;
  COFD[         776] =  -0.2029800570223752E+02;
  COFD[         777] =   0.5167991782317953E+01;
  COFD[         778] =  -0.4260628585538284E+00;
  COFD[         779] =   0.1725358239585169E-01;
  COFD[         780] =  -0.1969226895876598E+02;
  COFD[         781] =   0.4975023535625081E+01;
  COFD[         782] =  -0.4062730945223436E+00;
  COFD[         783] =   0.1660121367605663E-01;
  COFD[         784] =  -0.1967561834414118E+02;
  COFD[         785] =   0.4964896032502449E+01;
  COFD[         786] =  -0.4047371090295282E+00;
  COFD[         787] =   0.1652517237194395E-01;
  COFD[         788] =  -0.1924815690645235E+02;
  COFD[         789] =   0.4356731855020855E+01;
  COFD[         790] =  -0.2835834611821879E+00;
  COFD[         791] =   0.9647123135613308E-02;
  COFD[         792] =  -0.1835811548485455E+02;
  COFD[         793] =   0.3820061187148583E+01;
  COFD[         794] =  -0.1987078411533762E+00;
  COFD[         795] =   0.5458528950815168E-02;
  COFD[         796] =  -0.1814291631525713E+02;
  COFD[         797] =   0.3646821154166707E+01;
  COFD[         798] =  -0.1727110481572554E+00;
  COFD[         799] =   0.4220560947533801E-02;
  COFD[         800] =  -0.1920855889011546E+02;
  COFD[         801] =   0.4192213485908265E+01;
  COFD[         802] =  -0.2588520836689803E+00;
  COFD[         803] =   0.8501867421401927E-02;
  COFD[         804] =  -0.1874280769185066E+02;
  COFD[         805] =   0.3905545522430255E+01;
  COFD[         806] =  -0.2138630507050396E+00;
  COFD[         807] =   0.6308839293777503E-02;
  COFD[         808] =  -0.1873753222748401E+02;
  COFD[         809] =   0.3902621516519401E+01;
  COFD[         810] =  -0.2134883091855596E+00;
  COFD[         811] =   0.6294056333331442E-02;
  COFD[         812] =  -0.1898993862716132E+02;
  COFD[         813] =   0.4082522815084645E+01;
  COFD[         814] =  -0.2435514237071018E+00;
  COFD[         815] =   0.7811876908468279E-02;
  COFD[         816] =  -0.1885646323568239E+02;
  COFD[         817] =   0.4016884863488960E+01;
  COFD[         818] =  -0.2343948477217928E+00;
  COFD[         819] =   0.7398933738679844E-02;
  COFD[         820] =  -0.1885646323568239E+02;
  COFD[         821] =   0.4016884863488960E+01;
  COFD[         822] =  -0.2343948477217928E+00;
  COFD[         823] =   0.7398933738679844E-02;
  COFD[         824] =  -0.1883841467755942E+02;
  COFD[         825] =   0.4008070394200236E+01;
  COFD[         826] =  -0.2331651697525292E+00;
  COFD[         827] =   0.7343476142094791E-02;
  COFD[         828] =  -0.1872172887795907E+02;
  COFD[         829] =   0.3951398222423000E+01;
  COFD[         830] =  -0.2252586833510302E+00;
  COFD[         831] =   0.6986886080269861E-02;
  COFD[         832] =  -0.1861792788931437E+02;
  COFD[         833] =   0.3901413915830852E+01;
  COFD[         834] =  -0.2182846281641017E+00;
  COFD[         835] =   0.6672322090977868E-02;
  COFD[         836] =  -0.2024331300258610E+02;
  COFD[         837] =   0.5167058816761965E+01;
  COFD[         838] =  -0.4292618020154789E+00;
  COFD[         839] =   0.1752039422653532E-01;
  COFD[         840] =  -0.1699415309027476E+02;
  COFD[         841] =   0.4447105495943485E+01;
  COFD[         842] =  -0.3323553153514281E+00;
  COFD[         843] =   0.1310173916235249E-01;
  COFD[         844] =  -0.1586051006549287E+02;
  COFD[         845] =   0.4116349686729230E+01;
  COFD[         846] =  -0.3299412099375534E+00;
  COFD[         847] =   0.1473720935464579E-01;
  COFD[         848] =  -0.1938215107014084E+02;
  COFD[         849] =   0.4860644015548545E+01;
  COFD[         850] =  -0.4060689337372705E+00;
  COFD[         851] =   0.1718274802992044E-01;
  COFD[         852] =  -0.1784670232120704E+02;
  COFD[         853] =   0.4483333996808572E+01;
  COFD[         854] =  -0.3684301675626975E+00;
  COFD[         855] =   0.1601948073360996E-01;
  COFD[         856] =  -0.1933923251029205E+02;
  COFD[         857] =   0.4846414430750782E+01;
  COFD[         858] =  -0.4047681110056587E+00;
  COFD[         859] =   0.1714937694388643E-01;
  COFD[         860] =  -0.1783539626494398E+02;
  COFD[         861] =   0.4471504634175316E+01;
  COFD[         862] =  -0.3669099771101654E+00;
  COFD[         863] =   0.1595459643293670E-01;
  COFD[         864] =  -0.2000931867510009E+02;
  COFD[         865] =   0.4769188096212780E+01;
  COFD[         866] =  -0.3467567852966127E+00;
  COFD[         867] =   0.1276719834602697E-01;
  COFD[         868] =  -0.2082280203393144E+02;
  COFD[         869] =   0.5221216866940531E+01;
  COFD[         870] =  -0.4376943707328301E+00;
  COFD[         871] =   0.1795212687201138E-01;
  COFD[         872] =  -0.1850939531817895E+02;
  COFD[         873] =   0.4544732744781582E+01;
  COFD[         874] =  -0.3730497114680474E+00;
  COFD[         875] =   0.1608218736056522E-01;
  COFD[         876] =  -0.2050228473364046E+02;
  COFD[         877] =   0.5141702811583297E+01;
  COFD[         878] =  -0.4303114088412218E+00;
  COFD[         879] =   0.1774700289277886E-01;
  COFD[         880] =  -0.2079667220948008E+02;
  COFD[         881] =   0.5153157167644873E+01;
  COFD[         882] =  -0.4303605765762956E+00;
  COFD[         883] =   0.1769056906056647E-01;
  COFD[         884] =  -0.2091089956165765E+02;
  COFD[         885] =   0.5028244114071345E+01;
  COFD[         886] =  -0.3884228452359734E+00;
  COFD[         887] =   0.1487815240691972E-01;
  COFD[         888] =  -0.2079503005458607E+02;
  COFD[         889] =   0.5148706542754562E+01;
  COFD[         890] =  -0.4296850962181173E+00;
  COFD[         891] =   0.1765708156263098E-01;
  COFD[         892] =  -0.2101915140926890E+02;
  COFD[         893] =   0.5126566140193276E+01;
  COFD[         894] =  -0.4085741082041421E+00;
  COFD[         895] =   0.1603488286694611E-01;
  COFD[         896] =  -0.1881310947360568E+02;
  COFD[         897] =   0.4658928999542997E+01;
  COFD[         898] =  -0.3866746793655438E+00;
  COFD[         899] =   0.1662299447760606E-01;
  COFD[         900] =  -0.1882047781345688E+02;
  COFD[         901] =   0.4658986213777973E+01;
  COFD[         902] =  -0.3866620595057007E+00;
  COFD[         903] =   0.1662154903565528E-01;
  COFD[         904] =  -0.1882775755549848E+02;
  COFD[         905] =   0.4659119842910808E+01;
  COFD[         906] =  -0.3866568941685569E+00;
  COFD[         907] =   0.1662031309488973E-01;
  COFD[         908] =  -0.2058394974760212E+02;
  COFD[         909] =   0.5090753092957280E+01;
  COFD[         910] =  -0.4215256831575097E+00;
  COFD[         911] =   0.1727466009556999E-01;
  COFD[         912] =  -0.2100066296086839E+02;
  COFD[         913] =   0.5035858183028025E+01;
  COFD[         914] =  -0.3937618393670701E+00;
  COFD[         915] =   0.1527466676375227E-01;
  COFD[         916] =  -0.2092466715580316E+02;
  COFD[         917] =   0.4921376438675481E+01;
  COFD[         918] =  -0.3746202469475856E+00;
  COFD[         919] =   0.1427858739868015E-01;
  COFD[         920] =  -0.2097636580694235E+02;
  COFD[         921] =   0.5051780969690284E+01;
  COFD[         922] =  -0.4052813410201178E+00;
  COFD[         923] =   0.1613327373268888E-01;
  COFD[         924] =  -0.2080529541825178E+02;
  COFD[         925] =   0.4870282456696042E+01;
  COFD[         926] =  -0.3717177531779383E+00;
  COFD[         927] =   0.1428959829981116E-01;
  COFD[         928] =  -0.2079012230483711E+02;
  COFD[         929] =   0.4862345324915401E+01;
  COFD[         930] =  -0.3705511888834099E+00;
  COFD[         931] =   0.1423346818325797E-01;
  COFD[         932] =  -0.2073725116067806E+02;
  COFD[         933] =   0.4921491851668374E+01;
  COFD[         934] =  -0.3859047434126549E+00;
  COFD[         935] =   0.1518991301474521E-01;
  COFD[         936] =  -0.2057011439802693E+02;
  COFD[         937] =   0.4834011551777849E+01;
  COFD[         938] =  -0.3728956950726914E+00;
  COFD[         939] =   0.1455662572760375E-01;
  COFD[         940] =  -0.2057011439802693E+02;
  COFD[         941] =   0.4834011551777849E+01;
  COFD[         942] =  -0.3728956950726914E+00;
  COFD[         943] =   0.1455662572760375E-01;
  COFD[         944] =  -0.2054662946443170E+02;
  COFD[         945] =   0.4821857326037926E+01;
  COFD[         946] =  -0.3710883189732943E+00;
  COFD[         947] =   0.1446864565064257E-01;
  COFD[         948] =  -0.2039108134211626E+02;
  COFD[         949] =   0.4741987400654799E+01;
  COFD[         950] =  -0.3592117190236912E+00;
  COFD[         951] =   0.1389053577261777E-01;
  COFD[         952] =  -0.2024909562164469E+02;
  COFD[         953] =   0.4669841294589520E+01;
  COFD[         954] =  -0.3484841342800297E+00;
  COFD[         955] =   0.1336838998616518E-01;
  COFD[         956] =  -0.1848414582884554E+02;
  COFD[         957] =   0.4538825155691820E+01;
  COFD[         958] =  -0.3723671333664587E+00;
  COFD[         959] =   0.1605610876285614E-01;
  COFD[         960] =  -0.1719132602269712E+02;
  COFD[         961] =   0.4859367175853364E+01;
  COFD[         962] =  -0.4240580935572132E+00;
  COFD[         963] =   0.1870696403708276E-01;
  COFD[         964] =  -0.1268296956267843E+02;
  COFD[         965] =   0.2903428401599904E+01;
  COFD[         966] =  -0.1707667302851327E+00;
  COFD[         967] =   0.7771111802295727E-02;
  COFD[         968] =  -0.1661214673917502E+02;
  COFD[         969] =   0.3930381940303202E+01;
  COFD[         970] =  -0.2994310262293675E+00;
  COFD[         971] =   0.1314286546846328E-01;
  COFD[         972] =  -0.1478198899238034E+02;
  COFD[         973] =   0.3371092581023257E+01;
  COFD[         974] =  -0.2298567856863946E+00;
  COFD[         975] =   0.1025846107153025E-01;
  COFD[         976] =  -0.1654504240089281E+02;
  COFD[         977] =   0.3902990496709335E+01;
  COFD[         978] =  -0.2959726841289695E+00;
  COFD[         979] =   0.1299731318228341E-01;
  COFD[         980] =  -0.1476985219999908E+02;
  COFD[         981] =   0.3357751196108267E+01;
  COFD[         982] =  -0.2280046831004087E+00;
  COFD[         983] =   0.1017303667696902E-01;
  COFD[         984] =  -0.2031234004461199E+02;
  COFD[         985] =   0.5186312454417227E+01;
  COFD[         986] =  -0.4313742876779729E+00;
  COFD[         987] =   0.1759974844922337E-01;
  COFD[         988] =  -0.1850939531817895E+02;
  COFD[         989] =   0.4544732744781582E+01;
  COFD[         990] =  -0.3730497114680474E+00;
  COFD[         991] =   0.1608218736056522E-01;
  COFD[         992] =  -0.1551400389810888E+02;
  COFD[         993] =   0.3472777644794381E+01;
  COFD[         994] =  -0.2416432013010369E+00;
  COFD[         995] =   0.1070799930256162E-01;
  COFD[         996] =  -0.1819896532390013E+02;
  COFD[         997] =   0.4458091795954766E+01;
  COFD[         998] =  -0.3635634287967068E+00;
  COFD[         999] =   0.1574446016490575E-01;
  COFD[        1000] =  -0.1859114968467777E+02;
  COFD[        1001] =   0.4515459784417187E+01;
  COFD[        1002] =  -0.3706448354090633E+00;
  COFD[        1003] =   0.1603614788210589E-01;
  COFD[        1004] =  -0.2009778046496907E+02;
  COFD[        1005] =   0.5037181430781024E+01;
  COFD[        1006] =  -0.4240362270814515E+00;
  COFD[        1007] =   0.1776546650877937E-01;
  COFD[        1008] =  -0.1860391914467572E+02;
  COFD[        1009] =   0.4517177293840912E+01;
  COFD[        1010] =  -0.3708581736449252E+00;
  COFD[        1011] =   0.1604490864407949E-01;
  COFD[        1012] =  -0.2003301781547066E+02;
  COFD[        1013] =   0.5028625403684593E+01;
  COFD[        1014] =  -0.4259384280972174E+00;
  COFD[        1015] =   0.1797007932483027E-01;
  COFD[        1016] =  -0.1582244348422915E+02;
  COFD[        1017] =   0.3599211636009877E+01;
  COFD[        1018] =  -0.2581898369286856E+00;
  COFD[        1019] =   0.1143096575439008E-01;
  COFD[        1020] =  -0.1583284545979817E+02;
  COFD[        1021] =   0.3600569993944359E+01;
  COFD[        1022] =  -0.2583773478248265E+00;
  COFD[        1023] =   0.1143956672273590E-01;
  COFD[        1024] =  -0.1584365035133258E+02;
  COFD[        1025] =   0.3602236757402370E+01;
  COFD[        1026] =  -0.2586073726747737E+00;
  COFD[        1027] =   0.1145011549760728E-01;
  COFD[        1028] =  -0.1853079633836986E+02;
  COFD[        1029] =   0.4516962151189352E+01;
  COFD[        1030] =  -0.3708427592293742E+00;
  COFD[        1031] =   0.1604368344993687E-01;
  COFD[        1032] =  -0.1996442254756272E+02;
  COFD[        1033] =   0.4938191170899541E+01;
  COFD[        1034] =  -0.4135433397485079E+00;
  COFD[        1035] =   0.1740206709551280E-01;
  COFD[        1036] =  -0.2012556972106850E+02;
  COFD[        1037] =   0.4927762527501130E+01;
  COFD[        1038] =  -0.4091613023934828E+00;
  COFD[        1039] =   0.1708636787303381E-01;
  COFD[        1040] =  -0.1954575230865094E+02;
  COFD[        1041] =   0.4765735338615951E+01;
  COFD[        1042] =  -0.3970172796487949E+00;
  COFD[        1043] =   0.1692142125782123E-01;
  COFD[        1044] =  -0.1997997031335933E+02;
  COFD[        1045] =   0.4852495827753419E+01;
  COFD[        1046] =  -0.4017712513240834E+00;
  COFD[        1047] =   0.1685558673482876E-01;
  COFD[        1048] =  -0.1997892794374273E+02;
  COFD[        1049] =   0.4850385203230214E+01;
  COFD[        1050] =  -0.4014003482747929E+00;
  COFD[        1051] =   0.1683520757557910E-01;
  COFD[        1052] =  -0.1962170244526679E+02;
  COFD[        1053] =   0.4768861169263409E+01;
  COFD[        1054] =  -0.3963823359686274E+00;
  COFD[        1055] =   0.1684681197314066E-01;
  COFD[        1056] =  -0.1966500520249660E+02;
  COFD[        1057] =   0.4770450920673444E+01;
  COFD[        1058] =  -0.3958664514070914E+00;
  COFD[        1059] =   0.1679182031481241E-01;
  COFD[        1060] =  -0.1966500520249660E+02;
  COFD[        1061] =   0.4770450920673444E+01;
  COFD[        1062] =  -0.3958664514070914E+00;
  COFD[        1063] =   0.1679182031481241E-01;
  COFD[        1064] =  -0.1967071253473714E+02;
  COFD[        1065] =   0.4770649303948582E+01;
  COFD[        1066] =  -0.3957907059871365E+00;
  COFD[        1067] =   0.1678395401774585E-01;
  COFD[        1068] =  -0.1970677705449042E+02;
  COFD[        1069] =   0.4771839592787130E+01;
  COFD[        1070] =  -0.3952722585184899E+00;
  COFD[        1071] =   0.1673110821116964E-01;
  COFD[        1072] =  -0.1973759080924860E+02;
  COFD[        1073] =   0.4772765644230508E+01;
  COFD[        1074] =  -0.3947767037628411E+00;
  COFD[        1075] =   0.1668185546520601E-01;
  COFD[        1076] =  -0.1549017806948574E+02;
  COFD[        1077] =   0.3466859355125401E+01;
  COFD[        1078] =  -0.2408856054712291E+00;
  COFD[        1079] =   0.1067561900190753E-01;
  COFD[        1080] =  -0.1718718424688124E+02;
  COFD[        1081] =   0.4609957970156615E+01;
  COFD[        1082] =  -0.3605879229749922E+00;
  COFD[        1083] =   0.1462630202481900E-01;
  COFD[        1084] =  -0.1521342233375189E+02;
  COFD[        1085] =   0.3885324867997948E+01;
  COFD[        1086] =  -0.2989244982570487E+00;
  COFD[        1087] =   0.1334849796046345E-01;
  COFD[        1088] =  -0.1906497740611485E+02;
  COFD[        1089] =   0.4777590849182513E+01;
  COFD[        1090] =  -0.3975472982609866E+00;
  COFD[        1091] =   0.1690531426791657E-01;
  COFD[        1092] =  -0.1741286728352899E+02;
  COFD[        1093] =   0.4346824031835780E+01;
  COFD[        1094] =  -0.3515350630281370E+00;
  COFD[        1095] =   0.1532053929906590E-01;
  COFD[        1096] =  -0.1902332581359150E+02;
  COFD[        1097] =   0.4763347082494895E+01;
  COFD[        1098] =  -0.3961877942189693E+00;
  COFD[        1099] =   0.1686706463721130E-01;
  COFD[        1100] =  -0.1740489593482223E+02;
  COFD[        1101] =   0.4336092039978253E+01;
  COFD[        1102] =  -0.3501568532068898E+00;
  COFD[        1103] =   0.1526174747325281E-01;
  COFD[        1104] =  -0.2016248856160786E+02;
  COFD[        1105] =   0.4899874571765012E+01;
  COFD[        1106] =  -0.3689624059605574E+00;
  COFD[        1107] =   0.1393109509915711E-01;
  COFD[        1108] =  -0.2050228473364046E+02;
  COFD[        1109] =   0.5141702811583297E+01;
  COFD[        1110] =  -0.4303114088412218E+00;
  COFD[        1111] =   0.1774700289277886E-01;
  COFD[        1112] =  -0.1819896532390013E+02;
  COFD[        1113] =   0.4458091795954766E+01;
  COFD[        1114] =  -0.3635634287967068E+00;
  COFD[        1115] =   0.1574446016490575E-01;
  COFD[        1116] =  -0.2033351260815833E+02;
  COFD[        1117] =   0.5124525436391270E+01;
  COFD[        1118] =  -0.4314296217915469E+00;
  COFD[        1119] =   0.1792928453096072E-01;
  COFD[        1120] =  -0.2057313309604061E+02;
  COFD[        1121] =   0.5110289634910290E+01;
  COFD[        1122] =  -0.4279140715878681E+00;
  COFD[        1123] =   0.1770783981675040E-01;
  COFD[        1124] =  -0.2094985654740211E+02;
  COFD[        1125] =   0.5103559016459362E+01;
  COFD[        1126] =  -0.4028197832671332E+00;
  COFD[        1127] =   0.1567912907759180E-01;
  COFD[        1128] =  -0.2057389670179838E+02;
  COFD[        1129] =   0.5106758677684447E+01;
  COFD[        1130] =  -0.4273719773584432E+00;
  COFD[        1131] =   0.1768071705566059E-01;
  COFD[        1132] =  -0.2108393478877553E+02;
  COFD[        1133] =   0.5212547607346515E+01;
  COFD[        1134] =  -0.4246122376091541E+00;
  COFD[        1135] =   0.1692225360937287E-01;
  COFD[        1136] =  -0.1839743677478215E+02;
  COFD[        1137] =   0.4528923888489936E+01;
  COFD[        1138] =  -0.3712220054808055E+00;
  COFD[        1139] =   0.1601226224608377E-01;
  COFD[        1140] =  -0.1840518682221029E+02;
  COFD[        1141] =   0.4529051601311094E+01;
  COFD[        1142] =  -0.3712235453019769E+00;
  COFD[        1143] =   0.1601165440687174E-01;
  COFD[        1144] =  -0.1841287661408400E+02;
  COFD[        1145] =   0.4529271857206017E+01;
  COFD[        1146] =  -0.3712347871837591E+00;
  COFD[        1147] =   0.1601136037617702E-01;
  COFD[        1148] =  -0.2041077552923072E+02;
  COFD[        1149] =   0.5069647730132703E+01;
  COFD[        1150] =  -0.4221575179904544E+00;
  COFD[        1151] =   0.1743626998603495E-01;
  COFD[        1152] =  -0.2096286670043905E+02;
  COFD[        1153] =   0.5073824718604745E+01;
  COFD[        1154] =  -0.4026978389039572E+00;
  COFD[        1155] =   0.1581329096487812E-01;
  COFD[        1156] =  -0.2089621715369740E+02;
  COFD[        1157] =   0.4962017913186494E+01;
  COFD[        1158] =  -0.3839107017164409E+00;
  COFD[        1159] =   0.1483109113820474E-01;
  COFD[        1160] =  -0.2090423070644263E+02;
  COFD[        1161] =   0.5069554235969105E+01;
  COFD[        1162] =  -0.4110387472810522E+00;
  COFD[        1163] =   0.1651243898718931E-01;
  COFD[        1164] =  -0.2072976841817859E+02;
  COFD[        1165] =   0.4886126782915337E+01;
  COFD[        1166] =  -0.3771688345094954E+00;
  COFD[        1167] =   0.1464949867829970E-01;
  COFD[        1168] =  -0.2071293986412085E+02;
  COFD[        1169] =   0.4877335837263783E+01;
  COFD[        1170] =  -0.3758684263056035E+00;
  COFD[        1171] =   0.1458649350932613E-01;
  COFD[        1172] =  -0.2065920762233094E+02;
  COFD[        1173] =   0.4934466495930646E+01;
  COFD[        1174] =  -0.3907995086461343E+00;
  COFD[        1175] =   0.1551999187366136E-01;
  COFD[        1176] =  -0.2048311817894017E+02;
  COFD[        1177] =   0.4841596720970757E+01;
  COFD[        1178] =  -0.3768868905348668E+00;
  COFD[        1179] =   0.1483785654563742E-01;
  COFD[        1180] =  -0.2048311817894017E+02;
  COFD[        1181] =   0.4841596720970757E+01;
  COFD[        1182] =  -0.3768868905348668E+00;
  COFD[        1183] =   0.1483785654563742E-01;
  COFD[        1184] =  -0.2045813151569588E+02;
  COFD[        1185] =   0.4828580179065280E+01;
  COFD[        1186] =  -0.3749369856707143E+00;
  COFD[        1187] =   0.1474225738632903E-01;
  COFD[        1188] =  -0.2029132185536470E+02;
  COFD[        1189] =   0.4742426481711147E+01;
  COFD[        1190] =  -0.3620314088800772E+00;
  COFD[        1191] =   0.1410955509763312E-01;
  COFD[        1192] =  -0.2013725011091987E+02;
  COFD[        1193] =   0.4663752728745078E+01;
  COFD[        1194] =  -0.3502469131795646E+00;
  COFD[        1195] =   0.1353185410858909E-01;
  COFD[        1196] =  -0.1816875575125756E+02;
  COFD[        1197] =   0.4450099332111185E+01;
  COFD[        1198] =  -0.3625779431233763E+00;
  COFD[        1199] =   0.1570387818274669E-01;
  COFD[        1200] =  -0.1739069348322545E+02;
  COFD[        1201] =   0.4517486135284712E+01;
  COFD[        1202] =  -0.3398734156166323E+00;
  COFD[        1203] =   0.1334270297473222E-01;
  COFD[        1204] =  -0.1609618892523237E+02;
  COFD[        1205] =   0.4147159660004064E+01;
  COFD[        1206] =  -0.3347661201338916E+00;
  COFD[        1207] =   0.1498106836636175E-01;
  COFD[        1208] =  -0.1938165462851885E+02;
  COFD[        1209] =   0.4799455794296509E+01;
  COFD[        1210] =  -0.3987831540937672E+00;
  COFD[        1211] =   0.1688981189405884E-01;
  COFD[        1212] =  -0.1797996041020607E+02;
  COFD[        1213] =   0.4462878081531843E+01;
  COFD[        1214] =  -0.3662487521791470E+00;
  COFD[        1215] =   0.1594206335130083E-01;
  COFD[        1216] =  -0.1933548765861729E+02;
  COFD[        1217] =   0.4783226123251535E+01;
  COFD[        1218] =  -0.3972356460765860E+00;
  COFD[        1219] =   0.1684656833081372E-01;
  COFD[        1220] =  -0.1796101043339682E+02;
  COFD[        1221] =   0.4447717004791725E+01;
  COFD[        1222] =  -0.3643006119047739E+00;
  COFD[        1223] =   0.1585890337276808E-01;
  COFD[        1224] =  -0.2002167887304272E+02;
  COFD[        1225] =   0.4706367098965421E+01;
  COFD[        1226] =  -0.3398295196679849E+00;
  COFD[        1227] =   0.1251278469477959E-01;
  COFD[        1228] =  -0.2079667220948008E+02;
  COFD[        1229] =   0.5153157167644873E+01;
  COFD[        1230] =  -0.4303605765762956E+00;
  COFD[        1231] =   0.1769056906056647E-01;
  COFD[        1232] =  -0.1859114968467777E+02;
  COFD[        1233] =   0.4515459784417187E+01;
  COFD[        1234] =  -0.3706448354090633E+00;
  COFD[        1235] =   0.1603614788210589E-01;
  COFD[        1236] =  -0.2057313309604061E+02;
  COFD[        1237] =   0.5110289634910290E+01;
  COFD[        1238] =  -0.4279140715878681E+00;
  COFD[        1239] =   0.1770783981675040E-01;
  COFD[        1240] =  -0.2086882696434525E+02;
  COFD[        1241] =   0.5134268558863930E+01;
  COFD[        1242] =  -0.4303935480868216E+00;
  COFD[        1243] =   0.1779289806072862E-01;
  COFD[        1244] =  -0.2116982137731672E+02;
  COFD[        1245] =   0.5082529061705947E+01;
  COFD[        1246] =  -0.3985695163760528E+00;
  COFD[        1247] =   0.1543772091972775E-01;
  COFD[        1248] =  -0.2087667755513102E+02;
  COFD[        1249] =   0.5133988787671894E+01;
  COFD[        1250] =  -0.4303508137436654E+00;
  COFD[        1251] =   0.1779076924798608E-01;
  COFD[        1252] =  -0.2127808242822088E+02;
  COFD[        1253] =   0.5181120604040352E+01;
  COFD[        1254] =  -0.4187531251865307E+00;
  COFD[        1255] =   0.1659700305087440E-01;
  COFD[        1256] =  -0.1876367440516571E+02;
  COFD[        1257] =   0.4572319637339744E+01;
  COFD[        1258] =  -0.3763258628262886E+00;
  COFD[        1259] =   0.1621138071251859E-01;
  COFD[        1260] =  -0.1876945901632086E+02;
  COFD[        1261] =   0.4571736685656203E+01;
  COFD[        1262] =  -0.3762545867091497E+00;
  COFD[        1263] =   0.1620851017345105E-01;
  COFD[        1264] =  -0.1877525427053113E+02;
  COFD[        1265] =   0.4571270559057186E+01;
  COFD[        1266] =  -0.3761952939001092E+00;
  COFD[        1267] =   0.1620601231486099E-01;
  COFD[        1268] =  -0.2074265405272435E+02;
  COFD[        1269] =   0.5101856595916880E+01;
  COFD[        1270] =  -0.4258231875144630E+00;
  COFD[        1271] =   0.1757612519891745E-01;
  COFD[        1272] =  -0.2133138415589590E+02;
  COFD[        1273] =   0.5124135795112915E+01;
  COFD[        1274] =  -0.4090120016728032E+00;
  COFD[        1275] =   0.1608355890403254E-01;
  COFD[        1276] =  -0.2136316488507793E+02;
  COFD[        1277] =   0.5059984882792701E+01;
  COFD[        1278] =  -0.3970696812921098E+00;
  COFD[        1279] =   0.1542543509204860E-01;
  COFD[        1280] =  -0.2137421544016685E+02;
  COFD[        1281] =   0.5174419154841422E+01;
  COFD[        1282] =  -0.4256066415182231E+00;
  COFD[        1283] =   0.1719123893975486E-01;
  COFD[        1284] =  -0.2130885185428343E+02;
  COFD[        1285] =   0.5038362567088167E+01;
  COFD[        1286] =  -0.3984885360207525E+00;
  COFD[        1287] =   0.1564530982215940E-01;
  COFD[        1288] =  -0.2129718866568306E+02;
  COFD[        1289] =   0.5031879745702688E+01;
  COFD[        1290] =  -0.3975318057376973E+00;
  COFD[        1291] =   0.1559908037190491E-01;
  COFD[        1292] =  -0.2120949634659705E+02;
  COFD[        1293] =   0.5075665220886366E+01;
  COFD[        1294] =  -0.4108520560652987E+00;
  COFD[        1295] =   0.1646976074078464E-01;
  COFD[        1296] =  -0.2107933630762049E+02;
  COFD[        1297] =   0.5003615183707528E+01;
  COFD[        1298] =  -0.4000889881542739E+00;
  COFD[        1299] =   0.1594355520326905E-01;
  COFD[        1300] =  -0.2107933630762049E+02;
  COFD[        1301] =   0.5003615183707528E+01;
  COFD[        1302] =  -0.4000889881542739E+00;
  COFD[        1303] =   0.1594355520326905E-01;
  COFD[        1304] =  -0.2106042187580892E+02;
  COFD[        1305] =   0.4993353000413181E+01;
  COFD[        1306] =  -0.3985560772762201E+00;
  COFD[        1307] =   0.1586861625332952E-01;
  COFD[        1308] =  -0.2093223365539976E+02;
  COFD[        1309] =   0.4924726826746843E+01;
  COFD[        1310] =  -0.3883055489377778E+00;
  COFD[        1311] =   0.1536753032747617E-01;
  COFD[        1312] =  -0.2081181618605448E+02;
  COFD[        1313] =   0.4861337517486841E+01;
  COFD[        1314] =  -0.3788379332137031E+00;
  COFD[        1315] =   0.1490475778206881E-01;
  COFD[        1316] =  -0.1856460854973782E+02;
  COFD[        1317] =   0.4508650490412063E+01;
  COFD[        1318] =  -0.3698234379228469E+00;
  COFD[        1319] =   0.1600311768013340E-01;
  COFD[        1320] =  -0.1539084008249366E+02;
  COFD[        1321] =   0.3565006575777335E+01;
  COFD[        1322] =  -0.1867580729432450E+00;
  COFD[        1323] =   0.5536884066371989E-02;
  COFD[        1324] =  -0.1747771970561251E+02;
  COFD[        1325] =   0.4673159679696371E+01;
  COFD[        1326] =  -0.3935109715079382E+00;
  COFD[        1327] =   0.1711386921544412E-01;
  COFD[        1328] =  -0.2045023591963550E+02;
  COFD[        1329] =   0.5108641161759671E+01;
  COFD[        1330] =  -0.4199387136777197E+00;
  COFD[        1331] =   0.1704933199723114E-01;
  COFD[        1332] =  -0.1914986770399175E+02;
  COFD[        1333] =   0.4864899621234062E+01;
  COFD[        1334] =  -0.4039220297382060E+00;
  COFD[        1335] =   0.1697258735720344E-01;
  COFD[        1336] =  -0.2047278277859408E+02;
  COFD[        1337] =   0.5123414930187060E+01;
  COFD[        1338] =  -0.4228706903648477E+00;
  COFD[        1339] =   0.1721819830097161E-01;
  COFD[        1340] =  -0.1918177070837826E+02;
  COFD[        1341] =   0.4871873199573439E+01;
  COFD[        1342] =  -0.4051602536432213E+00;
  COFD[        1343] =   0.1704106539271711E-01;
  COFD[        1344] =  -0.1757470287826831E+02;
  COFD[        1345] =   0.3509436887842453E+01;
  COFD[        1346] =  -0.1496102927025162E+00;
  COFD[        1347] =   0.2961974572337444E-02;
  COFD[        1348] =  -0.2091089956165765E+02;
  COFD[        1349] =   0.5028244114071345E+01;
  COFD[        1350] =  -0.3884228452359734E+00;
  COFD[        1351] =   0.1487815240691972E-01;
  COFD[        1352] =  -0.2009778046496907E+02;
  COFD[        1353] =   0.5037181430781024E+01;
  COFD[        1354] =  -0.4240362270814515E+00;
  COFD[        1355] =   0.1776546650877937E-01;
  COFD[        1356] =  -0.2094985654740211E+02;
  COFD[        1357] =   0.5103559016459362E+01;
  COFD[        1358] =  -0.4028197832671332E+00;
  COFD[        1359] =   0.1567912907759180E-01;
  COFD[        1360] =  -0.2116982137731672E+02;
  COFD[        1361] =   0.5082529061705947E+01;
  COFD[        1362] =  -0.3985695163760528E+00;
  COFD[        1363] =   0.1543772091972775E-01;
  COFD[        1364] =  -0.1887260515065304E+02;
  COFD[        1365] =   0.3923352162852601E+01;
  COFD[        1366] =  -0.2115865709699200E+00;
  COFD[        1367] =   0.5939285843445600E-02;
  COFD[        1368] =  -0.2117357921779836E+02;
  COFD[        1369] =   0.5080400399832968E+01;
  COFD[        1370] =  -0.3982616990582741E+00;
  COFD[        1371] =   0.1542316631896301E-01;
  COFD[        1372] =  -0.1982650591715289E+02;
  COFD[        1373] =   0.4379424471801094E+01;
  COFD[        1374] =  -0.2814553408580349E+00;
  COFD[        1375] =   0.9373962429419786E-02;
  COFD[        1376] =  -0.2032689236579443E+02;
  COFD[        1377] =   0.5114221260960875E+01;
  COFD[        1378] =  -0.4320136109530544E+00;
  COFD[        1379] =   0.1803190494311217E-01;
  COFD[        1380] =  -0.2033308000404864E+02;
  COFD[        1381] =   0.5113689932626547E+01;
  COFD[        1382] =  -0.4319400850127593E+00;
  COFD[        1383] =   0.1802856188319298E-01;
  COFD[        1384] =  -0.2033818516187319E+02;
  COFD[        1385] =   0.5112803625509553E+01;
  COFD[        1386] =  -0.4318111908846191E+00;
  COFD[        1387] =   0.1802241389102129E-01;
  COFD[        1388] =  -0.2106203504941583E+02;
  COFD[        1389] =   0.5066201505062354E+01;
  COFD[        1390] =  -0.3966848046994042E+00;
  COFD[        1391] =   0.1536587788350550E-01;
  COFD[        1392] =  -0.1974614299047013E+02;
  COFD[        1393] =   0.4263552495214677E+01;
  COFD[        1394] =  -0.2638401037897016E+00;
  COFD[        1395] =   0.8519318791865089E-02;
  COFD[        1396] =  -0.1947441791465671E+02;
  COFD[        1397] =   0.4068997090767398E+01;
  COFD[        1398] =  -0.2342259510016974E+00;
  COFD[        1399] =   0.7076017538245195E-02;
  COFD[        1400] =  -0.2074011798729127E+02;
  COFD[        1401] =   0.4716308251797323E+01;
  COFD[        1402] =  -0.3368455089670934E+00;
  COFD[        1403] =   0.1222481812086106E-01;
  COFD[        1404] =  -0.2009226460847481E+02;
  COFD[        1405] =   0.4335080602770606E+01;
  COFD[        1406] =  -0.2764812336210347E+00;
  COFD[        1407] =   0.9208446946513201E-02;
  COFD[        1408] =  -0.2008729443847752E+02;
  COFD[        1409] =   0.4331784646722137E+01;
  COFD[        1410] =  -0.2760303040339949E+00;
  COFD[        1411] =   0.9188681380243065E-02;
  COFD[        1412] =  -0.2060180540567790E+02;
  COFD[        1413] =   0.4632891986456493E+01;
  COFD[        1414] =  -0.3249863297154995E+00;
  COFD[        1415] =   0.1167639716644466E-01;
  COFD[        1416] =  -0.2049312404884751E+02;
  COFD[        1417] =   0.4573043093616539E+01;
  COFD[        1418] =  -0.3164763709052368E+00;
  COFD[        1419] =   0.1128280291029378E-01;
  COFD[        1420] =  -0.2049312404884751E+02;
  COFD[        1421] =   0.4573043093616539E+01;
  COFD[        1422] =  -0.3164763709052368E+00;
  COFD[        1423] =   0.1128280291029378E-01;
  COFD[        1424] =  -0.2047732010042510E+02;
  COFD[        1425] =   0.4564542740763746E+01;
  COFD[        1426] =  -0.3152676686066944E+00;
  COFD[        1427] =   0.1122689874475350E-01;
  COFD[        1428] =  -0.2036999878611920E+02;
  COFD[        1429] =   0.4507734427051779E+01;
  COFD[        1430] =  -0.3071897892401887E+00;
  COFD[        1431] =   0.1085329061385159E-01;
  COFD[        1432] =  -0.2026866351460656E+02;
  COFD[        1433] =   0.4455185153044745E+01;
  COFD[        1434] =  -0.2997175567738866E+00;
  COFD[        1435] =   0.1050770872018541E-01;
  COFD[        1436] =  -0.2007323696432831E+02;
  COFD[        1437] =   0.5032033200707714E+01;
  COFD[        1438] =  -0.4235009227577490E+00;
  COFD[        1439] =   0.1774762291235686E-01;
  COFD[        1440] =  -0.1738939651119507E+02;
  COFD[        1441] =   0.4516024567959285E+01;
  COFD[        1442] =  -0.3395710901746170E+00;
  COFD[        1443] =   0.1332458515958147E-01;
  COFD[        1444] =  -0.1610856624128618E+02;
  COFD[        1445] =   0.4151906359557744E+01;
  COFD[        1446] =  -0.3354185801849658E+00;
  COFD[        1447] =   0.1501089101157305E-01;
  COFD[        1448] =  -0.1938533704705245E+02;
  COFD[        1449] =   0.4797779503867271E+01;
  COFD[        1450] =  -0.3984320673303173E+00;
  COFD[        1451] =   0.1686862693067165E-01;
  COFD[        1452] =  -0.1800786212441529E+02;
  COFD[        1453] =   0.4471520139789249E+01;
  COFD[        1454] =  -0.3673588813728833E+00;
  COFD[        1455] =   0.1598943652325460E-01;
  COFD[        1456] =  -0.1933986861083130E+02;
  COFD[        1457] =   0.4781769764498772E+01;
  COFD[        1458] =  -0.3969195763875204E+00;
  COFD[        1459] =   0.1682718826455844E-01;
  COFD[        1460] =  -0.1798854490941215E+02;
  COFD[        1461] =   0.4456122648699593E+01;
  COFD[        1462] =  -0.3653808041282116E+00;
  COFD[        1463] =   0.1590501766674265E-01;
  COFD[        1464] =  -0.1999142375296304E+02;
  COFD[        1465] =   0.4690231520616589E+01;
  COFD[        1466] =  -0.3375026912192592E+00;
  COFD[        1467] =   0.1240317371013226E-01;
  COFD[        1468] =  -0.2079503005458607E+02;
  COFD[        1469] =   0.5148706542754562E+01;
  COFD[        1470] =  -0.4296850962181173E+00;
  COFD[        1471] =   0.1765708156263098E-01;
  COFD[        1472] =  -0.1860391914467572E+02;
  COFD[        1473] =   0.4517177293840912E+01;
  COFD[        1474] =  -0.3708581736449252E+00;
  COFD[        1475] =   0.1604490864407949E-01;
  COFD[        1476] =  -0.2057389670179838E+02;
  COFD[        1477] =   0.5106758677684447E+01;
  COFD[        1478] =  -0.4273719773584432E+00;
  COFD[        1479] =   0.1768071705566059E-01;
  COFD[        1480] =  -0.2087667755513102E+02;
  COFD[        1481] =   0.5133988787671894E+01;
  COFD[        1482] =  -0.4303508137436654E+00;
  COFD[        1483] =   0.1779076924798608E-01;
  COFD[        1484] =  -0.2117357921779836E+02;
  COFD[        1485] =   0.5080400399832968E+01;
  COFD[        1486] =  -0.3982616990582741E+00;
  COFD[        1487] =   0.1542316631896301E-01;
  COFD[        1488] =  -0.2088587463730100E+02;
  COFD[        1489] =   0.5134268558863905E+01;
  COFD[        1490] =  -0.4303935480868180E+00;
  COFD[        1491] =   0.1779289806072845E-01;
  COFD[        1492] =  -0.2128386918932766E+02;
  COFD[        1493] =   0.5179813730552705E+01;
  COFD[        1494] =  -0.4185624846862125E+00;
  COFD[        1495] =   0.1658789623841498E-01;
  COFD[        1496] =  -0.1877458843620260E+02;
  COFD[        1497] =   0.4573081324845698E+01;
  COFD[        1498] =  -0.3764161851618111E+00;
  COFD[        1499] =   0.1621488416322674E-01;
  COFD[        1500] =  -0.1878015806354195E+02;
  COFD[        1501] =   0.4572383256671336E+01;
  COFD[        1502] =  -0.3763335089270711E+00;
  COFD[        1503] =   0.1621168231729516E-01;
  COFD[        1504] =  -0.1878574706059010E+02;
  COFD[        1505] =   0.4571806260778968E+01;
  COFD[        1506] =  -0.3762632297433025E+00;
  COFD[        1507] =   0.1620886476912668E-01;
  COFD[        1508] =  -0.2076217207534043E+02;
  COFD[        1509] =   0.5106213575785350E+01;
  COFD[        1510] =  -0.4264909663519035E+00;
  COFD[        1511] =   0.1760950164404084E-01;
  COFD[        1512] =  -0.2135679402100113E+02;
  COFD[        1513] =   0.5130945196700601E+01;
  COFD[        1514] =  -0.4100141138210215E+00;
  COFD[        1515] =   0.1613184415858697E-01;
  COFD[        1516] =  -0.2139377655814223E+02;
  COFD[        1517] =   0.5069003650557725E+01;
  COFD[        1518] =  -0.3983908652693064E+00;
  COFD[        1519] =   0.1548879283162612E-01;
  COFD[        1520] =  -0.2140524845935924E+02;
  COFD[        1521] =   0.5183807431960805E+01;
  COFD[        1522] =  -0.4270095585158289E+00;
  COFD[        1523] =   0.1725985135601222E-01;
  COFD[        1524] =  -0.2134888342401024E+02;
  COFD[        1525] =   0.5051235704256541E+01;
  COFD[        1526] =  -0.4003883872653834E+00;
  COFD[        1527] =   0.1573711341614223E-01;
  COFD[        1528] =  -0.2133756080522132E+02;
  COFD[        1529] =   0.5044876892625164E+01;
  COFD[        1530] =  -0.3994499299174794E+00;
  COFD[        1531] =   0.1569176539047349E-01;
  COFD[        1532] =  -0.2124842946919059E+02;
  COFD[        1533] =   0.5088113178440255E+01;
  COFD[        1534] =  -0.4127116950602789E+00;
  COFD[        1535] =   0.1656068511393674E-01;
  COFD[        1536] =  -0.2112052652732195E+02;
  COFD[        1537] =   0.5016777693101416E+01;
  COFD[        1538] =  -0.4020551643132008E+00;
  COFD[        1539] =   0.1603967677731761E-01;
  COFD[        1540] =  -0.2112052652732195E+02;
  COFD[        1541] =   0.5016777693101416E+01;
  COFD[        1542] =  -0.4020551643132008E+00;
  COFD[        1543] =   0.1603967677731761E-01;
  COFD[        1544] =  -0.2110177556093154E+02;
  COFD[        1545] =   0.5006548246650174E+01;
  COFD[        1546] =  -0.4005271173561529E+00;
  COFD[        1547] =   0.1596497410669515E-01;
  COFD[        1548] =  -0.2097388276834499E+02;
  COFD[        1549] =   0.4937790584419703E+01;
  COFD[        1550] =  -0.3902567887483665E+00;
  COFD[        1551] =   0.1546291085732362E-01;
  COFD[        1552] =  -0.2085271434976424E+02;
  COFD[        1553] =   0.4873828197990187E+01;
  COFD[        1554] =  -0.3807034496066365E+00;
  COFD[        1555] =   0.1499594029854510E-01;
  COFD[        1556] =  -0.1857749863267182E+02;
  COFD[        1557] =   0.4510414228568457E+01;
  COFD[        1558] =  -0.3700425246017726E+00;
  COFD[        1559] =   0.1601211519431907E-01;
  COFD[        1560] =  -0.1561125740910447E+02;
  COFD[        1561] =   0.3698144089263720E+01;
  COFD[        1562] =  -0.2103080207621357E+00;
  COFD[        1563] =   0.6792204226221277E-02;
  COFD[        1564] =  -0.1797854335218850E+02;
  COFD[        1565] =   0.4905020138101169E+01;
  COFD[        1566] =  -0.4293084201751438E+00;
  COFD[        1567] =   0.1891302999394448E-01;
  COFD[        1568] =  -0.2012688925964889E+02;
  COFD[        1569] =   0.5018247199767933E+01;
  COFD[        1570] =  -0.4131414445091446E+00;
  COFD[        1571] =   0.1694646037906594E-01;
  COFD[        1572] =  -0.1878355293596396E+02;
  COFD[        1573] =   0.4745970603113269E+01;
  COFD[        1574] =  -0.3924372163577848E+00;
  COFD[        1575] =   0.1663659729501466E-01;
  COFD[        1576] =  -0.2034321636694622E+02;
  COFD[        1577] =   0.5089562300279502E+01;
  COFD[        1578] =  -0.4211891334428361E+00;
  COFD[        1579] =   0.1724909778947470E-01;
  COFD[        1580] =  -0.1880409131690931E+02;
  COFD[        1581] =   0.4747896766189557E+01;
  COFD[        1582] =  -0.3929513558466455E+00;
  COFD[        1583] =   0.1667071927589917E-01;
  COFD[        1584] =  -0.1619587514529141E+02;
  COFD[        1585] =   0.2843235032975332E+01;
  COFD[        1586] =  -0.5104314541884596E-01;
  COFD[        1587] =  -0.1707937587486378E-02;
  COFD[        1588] =  -0.2101915140926890E+02;
  COFD[        1589] =   0.5126566140193276E+01;
  COFD[        1590] =  -0.4085741082041421E+00;
  COFD[        1591] =   0.1603488286694611E-01;
  COFD[        1592] =  -0.2003301781547066E+02;
  COFD[        1593] =   0.5028625403684593E+01;
  COFD[        1594] =  -0.4259384280972174E+00;
  COFD[        1595] =   0.1797007932483027E-01;
  COFD[        1596] =  -0.2108393478877553E+02;
  COFD[        1597] =   0.5212547607346515E+01;
  COFD[        1598] =  -0.4246122376091541E+00;
  COFD[        1599] =   0.1692225360937287E-01;
  COFD[        1600] =  -0.2127808242822088E+02;
  COFD[        1601] =   0.5181120604040352E+01;
  COFD[        1602] =  -0.4187531251865307E+00;
  COFD[        1603] =   0.1659700305087440E-01;
  COFD[        1604] =  -0.1982650591715289E+02;
  COFD[        1605] =   0.4379424471801094E+01;
  COFD[        1606] =  -0.2814553408580349E+00;
  COFD[        1607] =   0.9373962429419786E-02;
  COFD[        1608] =  -0.2128386918932766E+02;
  COFD[        1609] =   0.5179813730552705E+01;
  COFD[        1610] =  -0.4185624846862125E+00;
  COFD[        1611] =   0.1658789623841498E-01;
  COFD[        1612] =  -0.1906244143057827E+02;
  COFD[        1613] =   0.4049325847301449E+01;
  COFD[        1614] =  -0.2359825286728792E+00;
  COFD[        1615] =   0.7307807686876559E-02;
  COFD[        1616] =  -0.2016013281728147E+02;
  COFD[        1617] =   0.5063566265537697E+01;
  COFD[        1618] =  -0.4281515792120217E+00;
  COFD[        1619] =   0.1797180686456706E-01;
  COFD[        1620] =  -0.1994316533928780E+02;
  COFD[        1621] =   0.4989293390342089E+01;
  COFD[        1622] =  -0.4200959859596722E+00;
  COFD[        1623] =   0.1768643162186107E-01;
  COFD[        1624] =  -0.1994945707378163E+02;
  COFD[        1625] =   0.4988909811915903E+01;
  COFD[        1626] =  -0.4200448560138315E+00;
  COFD[        1627] =   0.1768418875431581E-01;
  COFD[        1628] =  -0.2111193727907456E+02;
  COFD[        1629] =   0.5115004803331022E+01;
  COFD[        1630] =  -0.4073081193984433E+00;
  COFD[        1631] =   0.1598596435486438E-01;
  COFD[        1632] =  -0.2039176723952609E+02;
  COFD[        1633] =   0.4589563967495875E+01;
  COFD[        1634] =  -0.3155845121772854E+00;
  COFD[        1635] =   0.1111636872703185E-01;
  COFD[        1636] =  -0.2021372325579300E+02;
  COFD[        1637] =   0.4433392392216290E+01;
  COFD[        1638] =  -0.2910494454945646E+00;
  COFD[        1639] =   0.9894154374316580E-02;
  COFD[        1640] =  -0.2108918057162479E+02;
  COFD[        1641] =   0.4917141448841289E+01;
  COFD[        1642] =  -0.3711916302984747E+00;
  COFD[        1643] =   0.1402561629509121E-01;
  COFD[        1644] =  -0.2056103945202781E+02;
  COFD[        1645] =   0.4582953188069841E+01;
  COFD[        1646] =  -0.3167741504985752E+00;
  COFD[        1647] =   0.1125103765063593E-01;
  COFD[        1648] =  -0.2055209768497027E+02;
  COFD[        1649] =   0.4577812802798165E+01;
  COFD[        1650] =  -0.3160440880823138E+00;
  COFD[        1651] =   0.1121730347594009E-01;
  COFD[        1652] =  -0.2092002080559684E+02;
  COFD[        1653] =   0.4818457997428789E+01;
  COFD[        1654] =  -0.3569233511421860E+00;
  COFD[        1655] =   0.1335173475837706E-01;
  COFD[        1656] =  -0.2078560813380206E+02;
  COFD[        1657] =   0.4746137595541170E+01;
  COFD[        1658] =  -0.3464661176921034E+00;
  COFD[        1659] =   0.1285782723518177E-01;
  COFD[        1660] =  -0.2078560813380206E+02;
  COFD[        1661] =   0.4746137595541170E+01;
  COFD[        1662] =  -0.3464661176921034E+00;
  COFD[        1663] =   0.1285782723518177E-01;
  COFD[        1664] =  -0.2076598860041296E+02;
  COFD[        1665] =   0.4735797568586775E+01;
  COFD[        1666] =  -0.3449709850285935E+00;
  COFD[        1667] =   0.1278721123756670E-01;
  COFD[        1668] =  -0.2063237685170529E+02;
  COFD[        1669] =   0.4666357783222749E+01;
  COFD[        1670] =  -0.3349302789942741E+00;
  COFD[        1671] =   0.1231299384240708E-01;
  COFD[        1672] =  -0.2050575086267360E+02;
  COFD[        1673] =   0.4601709555855408E+01;
  COFD[        1674] =  -0.3255825756010016E+00;
  COFD[        1675] =   0.1187152855312509E-01;
  COFD[        1676] =  -0.2000088490858467E+02;
  COFD[        1677] =   0.5021789747731799E+01;
  COFD[        1678] =  -0.4253480595637520E+00;
  COFD[        1679] =   0.1795667045275812E-01;
  COFD[        1680] =  -0.1761107854705398E+02;
  COFD[        1681] =   0.5048254268555682E+01;
  COFD[        1682] =  -0.4490719655111177E+00;
  COFD[        1683] =   0.1981332179135272E-01;
  COFD[        1684] =  -0.1270707915711888E+02;
  COFD[        1685] =   0.2929619774221139E+01;
  COFD[        1686] =  -0.1739670601708100E+00;
  COFD[        1687] =   0.7901198812704345E-02;
  COFD[        1688] =  -0.1707094794267795E+02;
  COFD[        1689] =   0.4119810624341889E+01;
  COFD[        1690] =  -0.3242844616145015E+00;
  COFD[        1691] =   0.1423184749658152E-01;
  COFD[        1692] =  -0.1511165027002312E+02;
  COFD[        1693] =   0.3515441252015493E+01;
  COFD[        1694] =  -0.2491198761221894E+00;
  COFD[        1695] =   0.1111625360314578E-01;
  COFD[        1696] =  -0.1696131964138645E+02;
  COFD[        1697] =   0.4074167708212651E+01;
  COFD[        1698] =  -0.3182576642303228E+00;
  COFD[        1699] =   0.1396612551510442E-01;
  COFD[        1700] =  -0.1509620843583920E+02;
  COFD[        1701] =   0.3500304707531062E+01;
  COFD[        1702] =  -0.2470160084537919E+00;
  COFD[        1703] =   0.1101910392929726E-01;
  COFD[        1704] =  -0.2029800570223752E+02;
  COFD[        1705] =   0.5167991782317953E+01;
  COFD[        1706] =  -0.4260628585538284E+00;
  COFD[        1707] =   0.1725358239585169E-01;
  COFD[        1708] =  -0.1881310947360568E+02;
  COFD[        1709] =   0.4658928999542997E+01;
  COFD[        1710] =  -0.3866746793655438E+00;
  COFD[        1711] =   0.1662299447760606E-01;
  COFD[        1712] =  -0.1582244348422915E+02;
  COFD[        1713] =   0.3599211636009877E+01;
  COFD[        1714] =  -0.2581898369286856E+00;
  COFD[        1715] =   0.1143096575439008E-01;
  COFD[        1716] =  -0.1839743677478215E+02;
  COFD[        1717] =   0.4528923888489936E+01;
  COFD[        1718] =  -0.3712220054808055E+00;
  COFD[        1719] =   0.1601226224608377E-01;
  COFD[        1720] =  -0.1876367440516571E+02;
  COFD[        1721] =   0.4572319637339744E+01;
  COFD[        1722] =  -0.3763258628262886E+00;
  COFD[        1723] =   0.1621138071251859E-01;
  COFD[        1724] =  -0.2032689236579443E+02;
  COFD[        1725] =   0.5114221260960875E+01;
  COFD[        1726] =  -0.4320136109530544E+00;
  COFD[        1727] =   0.1803190494311217E-01;
  COFD[        1728] =  -0.1877458843620260E+02;
  COFD[        1729] =   0.4573081324845698E+01;
  COFD[        1730] =  -0.3764161851618111E+00;
  COFD[        1731] =   0.1621488416322674E-01;
  COFD[        1732] =  -0.2016013281728147E+02;
  COFD[        1733] =   0.5063566265537697E+01;
  COFD[        1734] =  -0.4281515792120217E+00;
  COFD[        1735] =   0.1797180686456706E-01;
  COFD[        1736] =  -0.1605076758771229E+02;
  COFD[        1737] =   0.3688226344868140E+01;
  COFD[        1738] =  -0.2690947493495641E+00;
  COFD[        1739] =   0.1187504521154671E-01;
  COFD[        1740] =  -0.1605887758547131E+02;
  COFD[        1741] =   0.3688399632843664E+01;
  COFD[        1742] =  -0.2691183502617316E+00;
  COFD[        1743] =   0.1187611401222301E-01;
  COFD[        1744] =  -0.1606742328871227E+02;
  COFD[        1745] =   0.3688898239072069E+01;
  COFD[        1746] =  -0.2691862577212903E+00;
  COFD[        1747] =   0.1187918929470605E-01;
  COFD[        1748] =  -0.1868514242925209E+02;
  COFD[        1749] =   0.4565512377700616E+01;
  COFD[        1750] =  -0.3754538805580755E+00;
  COFD[        1751] =   0.1617305773699574E-01;
  COFD[        1752] =  -0.2009700773592217E+02;
  COFD[        1753] =   0.4972105108130217E+01;
  COFD[        1754] =  -0.4156366087533929E+00;
  COFD[        1755] =   0.1739973132656829E-01;
  COFD[        1756] =  -0.2038667713911414E+02;
  COFD[        1757] =   0.5013528836767304E+01;
  COFD[        1758] =  -0.4183960149034390E+00;
  COFD[        1759] =   0.1741045522425236E-01;
  COFD[        1760] =  -0.1967283644204365E+02;
  COFD[        1761] =   0.4798018532274038E+01;
  COFD[        1762] =  -0.3991636476917164E+00;
  COFD[        1763] =   0.1693109854880472E-01;
  COFD[        1764] =  -0.2012864081733589E+02;
  COFD[        1765] =   0.4891098053091741E+01;
  COFD[        1766] =  -0.4047491719429138E+00;
  COFD[        1767] =   0.1690406746445060E-01;
  COFD[        1768] =  -0.2012555233929497E+02;
  COFD[        1769] =   0.4887988169322795E+01;
  COFD[        1770] =  -0.4042349632777354E+00;
  COFD[        1771] =   0.1687690382084792E-01;
  COFD[        1772] =  -0.1971435410948289E+02;
  COFD[        1773] =   0.4784277383958439E+01;
  COFD[        1774] =  -0.3961230823099747E+00;
  COFD[        1775] =   0.1674302776354957E-01;
  COFD[        1776] =  -0.1973273884329878E+02;
  COFD[        1777] =   0.4773632942733207E+01;
  COFD[        1778] =  -0.3938198339615520E+00;
  COFD[        1779] =   0.1660198213420950E-01;
  COFD[        1780] =  -0.1973273884329878E+02;
  COFD[        1781] =   0.4773632942733207E+01;
  COFD[        1782] =  -0.3938198339615520E+00;
  COFD[        1783] =   0.1660198213420950E-01;
  COFD[        1784] =  -0.1973489604720685E+02;
  COFD[        1785] =   0.4772084197088009E+01;
  COFD[        1786] =  -0.3934868007358818E+00;
  COFD[        1787] =   0.1658164649850388E-01;
  COFD[        1788] =  -0.1974717304559367E+02;
  COFD[        1789] =   0.4761539000992064E+01;
  COFD[        1790] =  -0.3912295069504163E+00;
  COFD[        1791] =   0.1644410435383097E-01;
  COFD[        1792] =  -0.1975589425049316E+02;
  COFD[        1793] =   0.4751520609448664E+01;
  COFD[        1794] =  -0.3890982795753707E+00;
  COFD[        1795] =   0.1631462582834976E-01;
  COFD[        1796] =  -0.1578873902175925E+02;
  COFD[        1797] =   0.3589258101555060E+01;
  COFD[        1798] =  -0.2568718883472193E+00;
  COFD[        1799] =   0.1137271971440306E-01;
  COFD[        1800] =  -0.1762111333528786E+02;
  COFD[        1801] =   0.5052180357519856E+01;
  COFD[        1802] =  -0.4496001366665920E+00;
  COFD[        1803] =   0.1983696581079775E-01;
  COFD[        1804] =  -0.1271124360389947E+02;
  COFD[        1805] =   0.2931047949505553E+01;
  COFD[        1806] =  -0.1741679825526774E+00;
  COFD[        1807] =   0.7910585398301936E-02;
  COFD[        1808] =  -0.1709003104221693E+02;
  COFD[        1809] =   0.4125627998285617E+01;
  COFD[        1810] =  -0.3250670331906144E+00;
  COFD[        1811] =   0.1426686739847163E-01;
  COFD[        1812] =  -0.1513524421958843E+02;
  COFD[        1813] =   0.3523297948981071E+01;
  COFD[        1814] =  -0.2502118585975575E+00;
  COFD[        1815] =   0.1116667641067510E-01;
  COFD[        1816] =  -0.1698006191382479E+02;
  COFD[        1817] =   0.4079755891690069E+01;
  COFD[        1818] =  -0.3190093940831369E+00;
  COFD[        1819] =   0.1399976445173436E-01;
  COFD[        1820] =  -0.1511946285859442E+02;
  COFD[        1821] =   0.3507927470352835E+01;
  COFD[        1822] =  -0.2480755290441053E+00;
  COFD[        1823] =   0.1106802953504252E-01;
  COFD[        1824] =  -0.1969226895876598E+02;
  COFD[        1825] =   0.4975023535625081E+01;
  COFD[        1826] =  -0.4062730945223436E+00;
  COFD[        1827] =   0.1660121367605663E-01;
  COFD[        1828] =  -0.1882047781345688E+02;
  COFD[        1829] =   0.4658986213777973E+01;
  COFD[        1830] =  -0.3866620595057007E+00;
  COFD[        1831] =   0.1662154903565528E-01;
  COFD[        1832] =  -0.1583284545979817E+02;
  COFD[        1833] =   0.3600569993944359E+01;
  COFD[        1834] =  -0.2583773478248265E+00;
  COFD[        1835] =   0.1143956672273590E-01;
  COFD[        1836] =  -0.1840518682221029E+02;
  COFD[        1837] =   0.4529051601311094E+01;
  COFD[        1838] =  -0.3712235453019769E+00;
  COFD[        1839] =   0.1601165440687174E-01;
  COFD[        1840] =  -0.1876945901632086E+02;
  COFD[        1841] =   0.4571736685656203E+01;
  COFD[        1842] =  -0.3762545867091497E+00;
  COFD[        1843] =   0.1620851017345105E-01;
  COFD[        1844] =  -0.2033308000404864E+02;
  COFD[        1845] =   0.5113689932626547E+01;
  COFD[        1846] =  -0.4319400850127593E+00;
  COFD[        1847] =   0.1802856188319298E-01;
  COFD[        1848] =  -0.1878015806354195E+02;
  COFD[        1849] =   0.4572383256671336E+01;
  COFD[        1850] =  -0.3763335089270711E+00;
  COFD[        1851] =   0.1621168231729516E-01;
  COFD[        1852] =  -0.1994316533928780E+02;
  COFD[        1853] =   0.4989293390342089E+01;
  COFD[        1854] =  -0.4200959859596722E+00;
  COFD[        1855] =   0.1768643162186107E-01;
  COFD[        1856] =  -0.1605887758547131E+02;
  COFD[        1857] =   0.3688399632843664E+01;
  COFD[        1858] =  -0.2691183502617316E+00;
  COFD[        1859] =   0.1187611401222301E-01;
  COFD[        1860] =  -0.1606627473213039E+02;
  COFD[        1861] =   0.3688226344868143E+01;
  COFD[        1862] =  -0.2690947493495643E+00;
  COFD[        1863] =   0.1187504521154673E-01;
  COFD[        1864] =  -0.1607413035852147E+02;
  COFD[        1865] =   0.3688389366224249E+01;
  COFD[        1866] =  -0.2691169520025368E+00;
  COFD[        1867] =   0.1187605069009256E-01;
  COFD[        1868] =  -0.1868972720602865E+02;
  COFD[        1869] =   0.4564003101584411E+01;
  COFD[        1870] =  -0.3752914154261894E+00;
  COFD[        1871] =   0.1616756912994680E-01;
  COFD[        1872] =  -0.2011340049596478E+02;
  COFD[        1873] =   0.4975835496600746E+01;
  COFD[        1874] =  -0.4162401315456316E+00;
  COFD[        1875] =   0.1743113787867068E-01;
  COFD[        1876] =  -0.2040768203667939E+02;
  COFD[        1877] =   0.5019349654721370E+01;
  COFD[        1878] =  -0.4193216095365674E+00;
  COFD[        1879] =   0.1745804721026501E-01;
  COFD[        1880] =  -0.1968298220975192E+02;
  COFD[        1881] =   0.4799139981296912E+01;
  COFD[        1882] =  -0.3994189723094190E+00;
  COFD[        1883] =   0.1694708582407344E-01;
  COFD[        1884] =  -0.2014955127994030E+02;
  COFD[        1885] =   0.4896710570985273E+01;
  COFD[        1886] =  -0.4056774618307388E+00;
  COFD[        1887] =   0.1695311601042129E-01;
  COFD[        1888] =  -0.2014664412250720E+02;
  COFD[        1889] =   0.4893661499235972E+01;
  COFD[        1890] =  -0.4051731103333032E+00;
  COFD[        1891] =   0.1692646556995869E-01;
  COFD[        1892] =  -0.1972631705838112E+02;
  COFD[        1893] =   0.4785911211303439E+01;
  COFD[        1894] =  -0.3964794306558241E+00;
  COFD[        1895] =   0.1676492861574991E-01;
  COFD[        1896] =  -0.1974560293529215E+02;
  COFD[        1897] =   0.4775430480336665E+01;
  COFD[        1898] =  -0.3942069403173756E+00;
  COFD[        1899] =   0.1662563585220353E-01;
  COFD[        1900] =  -0.1974560293529215E+02;
  COFD[        1901] =   0.4775430480336665E+01;
  COFD[        1902] =  -0.3942069403173756E+00;
  COFD[        1903] =   0.1662563585220353E-01;
  COFD[        1904] =  -0.1974787011577460E+02;
  COFD[        1905] =   0.4773895101413700E+01;
  COFD[        1906] =  -0.3938762513976486E+00;
  COFD[        1907] =   0.1660542834670177E-01;
  COFD[        1908] =  -0.1976079148787423E+02;
  COFD[        1909] =   0.4763387289184921E+01;
  COFD[        1910] =  -0.3916239954434353E+00;
  COFD[        1911] =   0.1646810838530097E-01;
  COFD[        1912] =  -0.1976999225259944E+02;
  COFD[        1913] =   0.4753334674389479E+01;
  COFD[        1914] =  -0.3894833752885402E+00;
  COFD[        1915] =   0.1633799793889408E-01;
  COFD[        1916] =  -0.1579929195048582E+02;
  COFD[        1917] =   0.3590679879226839E+01;
  COFD[        1918] =  -0.2570681340587627E+00;
  COFD[        1919] =   0.1138172060899465E-01;
  COFD[        1920] =  -0.1763070827156596E+02;
  COFD[        1921] =   0.5055937247475734E+01;
  COFD[        1922] =  -0.4501055595390027E+00;
  COFD[        1923] =   0.1985959199754691E-01;
  COFD[        1924] =  -0.1271518733681371E+02;
  COFD[        1925] =   0.2932402241140174E+01;
  COFD[        1926] =  -0.1743585082612658E+00;
  COFD[        1927] =   0.7919486191765343E-02;
  COFD[        1928] =  -0.1710868501481421E+02;
  COFD[        1929] =   0.4131370391002061E+01;
  COFD[        1930] =  -0.3258395192522007E+00;
  COFD[        1931] =   0.1430143624455670E-01;
  COFD[        1932] =  -0.1515824683455671E+02;
  COFD[        1933] =   0.3531011422568503E+01;
  COFD[        1934] =  -0.2512839088019835E+00;
  COFD[        1935] =   0.1121617798528330E-01;
  COFD[        1936] =  -0.1699844377766270E+02;
  COFD[        1937] =   0.4085300971938842E+01;
  COFD[        1938] =  -0.3197553288341113E+00;
  COFD[        1939] =   0.1403314440725374E-01;
  COFD[        1940] =  -0.1514219804803715E+02;
  COFD[        1941] =   0.3515441252015498E+01;
  COFD[        1942] =  -0.2491198761221899E+00;
  COFD[        1943] =   0.1111625360314580E-01;
  COFD[        1944] =  -0.1967561834414118E+02;
  COFD[        1945] =   0.4964896032502449E+01;
  COFD[        1946] =  -0.4047371090295282E+00;
  COFD[        1947] =   0.1652517237194395E-01;
  COFD[        1948] =  -0.1882775755549848E+02;
  COFD[        1949] =   0.4659119842910808E+01;
  COFD[        1950] =  -0.3866568941685569E+00;
  COFD[        1951] =   0.1662031309488973E-01;
  COFD[        1952] =  -0.1584365035133258E+02;
  COFD[        1953] =   0.3602236757402370E+01;
  COFD[        1954] =  -0.2586073726747737E+00;
  COFD[        1955] =   0.1145011549760728E-01;
  COFD[        1956] =  -0.1841287661408400E+02;
  COFD[        1957] =   0.4529271857206017E+01;
  COFD[        1958] =  -0.3712347871837591E+00;
  COFD[        1959] =   0.1601136037617702E-01;
  COFD[        1960] =  -0.1877525427053113E+02;
  COFD[        1961] =   0.4571270559057186E+01;
  COFD[        1962] =  -0.3761952939001092E+00;
  COFD[        1963] =   0.1620601231486099E-01;
  COFD[        1964] =  -0.2033818516187319E+02;
  COFD[        1965] =   0.5112803625509553E+01;
  COFD[        1966] =  -0.4318111908846191E+00;
  COFD[        1967] =   0.1802241389102129E-01;
  COFD[        1968] =  -0.1878574706059010E+02;
  COFD[        1969] =   0.4571806260778968E+01;
  COFD[        1970] =  -0.3762632297433025E+00;
  COFD[        1971] =   0.1620886476912668E-01;
  COFD[        1972] =  -0.1994945707378163E+02;
  COFD[        1973] =   0.4988909811915903E+01;
  COFD[        1974] =  -0.4200448560138315E+00;
  COFD[        1975] =   0.1768418875431581E-01;
  COFD[        1976] =  -0.1606742328871227E+02;
  COFD[        1977] =   0.3688898239072069E+01;
  COFD[        1978] =  -0.2691862577212903E+00;
  COFD[        1979] =   0.1187918929470605E-01;
  COFD[        1980] =  -0.1607413035852147E+02;
  COFD[        1981] =   0.3688389366224249E+01;
  COFD[        1982] =  -0.2691169520025368E+00;
  COFD[        1983] =   0.1187605069009256E-01;
  COFD[        1984] =  -0.1608131536572632E+02;
  COFD[        1985] =   0.3688226344868147E+01;
  COFD[        1986] =  -0.2690947493495650E+00;
  COFD[        1987] =   0.1187504521154676E-01;
  COFD[        1988] =  -0.1869430855359004E+02;
  COFD[        1989] =   0.4562621031385030E+01;
  COFD[        1990] =  -0.3751419042160075E+00;
  COFD[        1991] =   0.1616247820310002E-01;
  COFD[        1992] =  -0.2012868226172504E+02;
  COFD[        1993] =   0.4979206995573055E+01;
  COFD[        1994] =  -0.4167863191869534E+00;
  COFD[        1995] =   0.1745958878239685E-01;
  COFD[        1996] =  -0.2042738004900579E+02;
  COFD[        1997] =   0.5024724404363650E+01;
  COFD[        1998] =  -0.4201768197596846E+00;
  COFD[        1999] =   0.1750204200208273E-01;
  COFD[        2000] =  -0.1969261804237637E+02;
  COFD[        2001] =   0.4800160734138452E+01;
  COFD[        2002] =  -0.3996532469352685E+00;
  COFD[        2003] =   0.1696180432302862E-01;
  COFD[        2004] =  -0.2016949900503337E+02;
  COFD[        2005] =   0.4902039892776266E+01;
  COFD[        2006] =  -0.4065592769263701E+00;
  COFD[        2007] =   0.1699972276094747E-01;
  COFD[        2008] =  -0.2016677937256390E+02;
  COFD[        2009] =   0.4899055461722977E+01;
  COFD[        2010] =  -0.4060654125033917E+00;
  COFD[        2011] =   0.1697361869912191E-01;
  COFD[        2012] =  -0.1973776103143215E+02;
  COFD[        2013] =   0.4787458937602214E+01;
  COFD[        2014] =  -0.3968179100773139E+00;
  COFD[        2015] =   0.1678575643964539E-01;
  COFD[        2016] =  -0.1975795254936833E+02;
  COFD[        2017] =   0.4777156336087093E+01;
  COFD[        2018] =  -0.3945792342207480E+00;
  COFD[        2019] =   0.1664840211599068E-01;
  COFD[        2020] =  -0.1975795254936833E+02;
  COFD[        2021] =   0.4777156336087093E+01;
  COFD[        2022] =  -0.3945792342207480E+00;
  COFD[        2023] =   0.1664840211599068E-01;
  COFD[        2024] =  -0.1976033064762442E+02;
  COFD[        2025] =   0.4775636407093550E+01;
  COFD[        2026] =  -0.3942513288166545E+00;
  COFD[        2027] =   0.1662834933434690E-01;
  COFD[        2028] =  -0.1977390351786128E+02;
  COFD[        2029] =   0.4765179530799094E+01;
  COFD[        2030] =  -0.3920069566617906E+00;
  COFD[        2031] =   0.1649142352328586E-01;
  COFD[        2032] =  -0.1978359074953881E+02;
  COFD[        2033] =   0.4755104128079081E+01;
  COFD[        2034] =  -0.3898593304886394E+00;
  COFD[        2035] =   0.1636082491100651E-01;
  COFD[        2036] =  -0.1581023736073806E+02;
  COFD[        2037] =   0.3592405682636056E+01;
  COFD[        2038] =  -0.2573062871053550E+00;
  COFD[        2039] =   0.1139264138769127E-01;
  COFD[        2040] =  -0.1725034548774479E+02;
  COFD[        2041] =   0.4566921339972978E+01;
  COFD[        2042] =  -0.3509192271095338E+00;
  COFD[        2043] =   0.1402685086906688E-01;
  COFD[        2044] =  -0.1564739122681956E+02;
  COFD[        2045] =   0.4025825276486324E+01;
  COFD[        2046] =  -0.3181195147607377E+00;
  COFD[        2047] =   0.1422136468533644E-01;
  COFD[        2048] =  -0.1923136095279614E+02;
  COFD[        2049] =   0.4771652981887103E+01;
  COFD[        2050] =  -0.3944968044729790E+00;
  COFD[        2051] =   0.1667283377652949E-01;
  COFD[        2052] =  -0.1793673861488134E+02;
  COFD[        2053] =   0.4489875914439151E+01;
  COFD[        2054] =  -0.3697556337020395E+00;
  COFD[        2055] =   0.1609253879718407E-01;
  COFD[        2056] =  -0.1920918070133062E+02;
  COFD[        2057] =   0.4764713252836613E+01;
  COFD[        2058] =  -0.3942405534249802E+00;
  COFD[        2059] =   0.1668923743050566E-01;
  COFD[        2060] =  -0.1792454027766008E+02;
  COFD[        2061] =   0.4476326984714260E+01;
  COFD[        2062] =  -0.3680191278135618E+00;
  COFD[        2063] =   0.1601860672433683E-01;
  COFD[        2064] =  -0.1924815690645235E+02;
  COFD[        2065] =   0.4356731855020855E+01;
  COFD[        2066] =  -0.2835834611821879E+00;
  COFD[        2067] =   0.9647123135613308E-02;
  COFD[        2068] =  -0.2058394974760212E+02;
  COFD[        2069] =   0.5090753092957280E+01;
  COFD[        2070] =  -0.4215256831575097E+00;
  COFD[        2071] =   0.1727466009556999E-01;
  COFD[        2072] =  -0.1853079633836986E+02;
  COFD[        2073] =   0.4516962151189352E+01;
  COFD[        2074] =  -0.3708427592293742E+00;
  COFD[        2075] =   0.1604368344993687E-01;
  COFD[        2076] =  -0.2041077552923072E+02;
  COFD[        2077] =   0.5069647730132703E+01;
  COFD[        2078] =  -0.4221575179904544E+00;
  COFD[        2079] =   0.1743626998603495E-01;
  COFD[        2080] =  -0.2074265405272435E+02;
  COFD[        2081] =   0.5101856595916880E+01;
  COFD[        2082] =  -0.4258231875144630E+00;
  COFD[        2083] =   0.1757612519891745E-01;
  COFD[        2084] =  -0.2106203504941583E+02;
  COFD[        2085] =   0.5066201505062354E+01;
  COFD[        2086] =  -0.3966848046994042E+00;
  COFD[        2087] =   0.1536587788350550E-01;
  COFD[        2088] =  -0.2076217207534043E+02;
  COFD[        2089] =   0.5106213575785350E+01;
  COFD[        2090] =  -0.4264909663519035E+00;
  COFD[        2091] =   0.1760950164404084E-01;
  COFD[        2092] =  -0.2111193727907456E+02;
  COFD[        2093] =   0.5115004803331022E+01;
  COFD[        2094] =  -0.4073081193984433E+00;
  COFD[        2095] =   0.1598596435486438E-01;
  COFD[        2096] =  -0.1868514242925209E+02;
  COFD[        2097] =   0.4565512377700616E+01;
  COFD[        2098] =  -0.3754538805580755E+00;
  COFD[        2099] =   0.1617305773699574E-01;
  COFD[        2100] =  -0.1868972720602865E+02;
  COFD[        2101] =   0.4564003101584411E+01;
  COFD[        2102] =  -0.3752914154261894E+00;
  COFD[        2103] =   0.1616756912994680E-01;
  COFD[        2104] =  -0.1869430855359004E+02;
  COFD[        2105] =   0.4562621031385030E+01;
  COFD[        2106] =  -0.3751419042160075E+00;
  COFD[        2107] =   0.1616247820310002E-01;
  COFD[        2108] =  -0.2074682089064003E+02;
  COFD[        2109] =   0.5125913401864712E+01;
  COFD[        2110] =  -0.4300036496008145E+00;
  COFD[        2111] =   0.1780187794631543E-01;
  COFD[        2112] =  -0.2141485527132566E+02;
  COFD[        2113] =   0.5180861875662476E+01;
  COFD[        2114] =  -0.4178144385214836E+00;
  COFD[        2115] =   0.1652354263065661E-01;
  COFD[        2116] =  -0.2148234754894915E+02;
  COFD[        2117] =   0.5128240108050982E+01;
  COFD[        2118] =  -0.4075464684376063E+00;
  COFD[        2119] =   0.1594434352315902E-01;
  COFD[        2120] =  -0.2148607955277492E+02;
  COFD[        2121] =   0.5238740013363649E+01;
  COFD[        2122] =  -0.4356664393291084E+00;
  COFD[        2123] =   0.1769834275551969E-01;
  COFD[        2124] =  -0.2152239257799227E+02;
  COFD[        2125] =   0.5141731472937461E+01;
  COFD[        2126] =  -0.4142263111871738E+00;
  COFD[        2127] =   0.1642254345298395E-01;
  COFD[        2128] =  -0.2151480952421122E+02;
  COFD[        2129] =   0.5136706180179550E+01;
  COFD[        2130] =  -0.4134839036187021E+00;
  COFD[        2131] =   0.1638663893119980E-01;
  COFD[        2132] =  -0.2141313429386006E+02;
  COFD[        2133] =   0.5174902392506038E+01;
  COFD[        2134] =  -0.4261176538552951E+00;
  COFD[        2135] =   0.1723100350853119E-01;
  COFD[        2136] =  -0.2131043155757551E+02;
  COFD[        2137] =   0.5111234850420700E+01;
  COFD[        2138] =  -0.4165974605140776E+00;
  COFD[        2139] =   0.1676520530944528E-01;
  COFD[        2140] =  -0.2131043155757551E+02;
  COFD[        2141] =   0.5111234850420700E+01;
  COFD[        2142] =  -0.4165974605140776E+00;
  COFD[        2143] =   0.1676520530944528E-01;
  COFD[        2144] =  -0.2129337978312540E+02;
  COFD[        2145] =   0.5101265807720457E+01;
  COFD[        2146] =  -0.4151069111142752E+00;
  COFD[        2147] =   0.1669228251038974E-01;
  COFD[        2148] =  -0.2116627201427422E+02;
  COFD[        2149] =   0.5029554692101149E+01;
  COFD[        2150] =  -0.4043854700274245E+00;
  COFD[        2151] =   0.1616778720354987E-01;
  COFD[        2152] =  -0.2103101981945703E+02;
  COFD[        2153] =   0.4956251882737971E+01;
  COFD[        2154] =  -0.3934269877726836E+00;
  COFD[        2155] =   0.1563174832696796E-01;
  COFD[        2156] =  -0.1850244542172855E+02;
  COFD[        2157] =   0.4509546506902997E+01;
  COFD[        2158] =  -0.3699250247632380E+00;
  COFD[        2159] =   0.1600569104209367E-01;
  COFD[        2160] =  -0.1508763665602342E+02;
  COFD[        2161] =   0.3400141263568119E+01;
  COFD[        2162] =  -0.1645508303328382E+00;
  COFD[        2163] =   0.4517926211356998E-02;
  COFD[        2164] =  -0.1770170994272074E+02;
  COFD[        2165] =   0.4743686490643380E+01;
  COFD[        2166] =  -0.4058992902065018E+00;
  COFD[        2167] =   0.1778646492198370E-01;
  COFD[        2168] =  -0.2002012999655281E+02;
  COFD[        2169] =   0.4895179416250934E+01;
  COFD[        2170] =  -0.3929898448063785E+00;
  COFD[        2171] =   0.1589904718416698E-01;
  COFD[        2172] =  -0.1900626883865293E+02;
  COFD[        2173] =   0.4754407177756363E+01;
  COFD[        2174] =  -0.3901998067876299E+00;
  COFD[        2175] =   0.1639341623231018E-01;
  COFD[        2176] =  -0.2005721601235365E+02;
  COFD[        2177] =   0.4914479851418061E+01;
  COFD[        2178] =  -0.3965430355140246E+00;
  COFD[        2179] =   0.1609667930913030E-01;
  COFD[        2180] =  -0.1903172737849495E+02;
  COFD[        2181] =   0.4757798023561699E+01;
  COFD[        2182] =  -0.3909686059341044E+00;
  COFD[        2183] =   0.1644144176618614E-01;
  COFD[        2184] =  -0.1835811548485455E+02;
  COFD[        2185] =   0.3820061187148583E+01;
  COFD[        2186] =  -0.1987078411533762E+00;
  COFD[        2187] =   0.5458528950815168E-02;
  COFD[        2188] =  -0.2100066296086839E+02;
  COFD[        2189] =   0.5035858183028025E+01;
  COFD[        2190] =  -0.3937618393670701E+00;
  COFD[        2191] =   0.1527466676375227E-01;
  COFD[        2192] =  -0.1996442254756272E+02;
  COFD[        2193] =   0.4938191170899541E+01;
  COFD[        2194] =  -0.4135433397485079E+00;
  COFD[        2195] =   0.1740206709551280E-01;
  COFD[        2196] =  -0.2096286670043905E+02;
  COFD[        2197] =   0.5073824718604745E+01;
  COFD[        2198] =  -0.4026978389039572E+00;
  COFD[        2199] =   0.1581329096487812E-01;
  COFD[        2200] =  -0.2133138415589590E+02;
  COFD[        2201] =   0.5124135795112915E+01;
  COFD[        2202] =  -0.4090120016728032E+00;
  COFD[        2203] =   0.1608355890403254E-01;
  COFD[        2204] =  -0.1974614299047013E+02;
  COFD[        2205] =   0.4263552495214677E+01;
  COFD[        2206] =  -0.2638401037897016E+00;
  COFD[        2207] =   0.8519318791865089E-02;
  COFD[        2208] =  -0.2135679402100113E+02;
  COFD[        2209] =   0.5130945196700601E+01;
  COFD[        2210] =  -0.4100141138210215E+00;
  COFD[        2211] =   0.1613184415858697E-01;
  COFD[        2212] =  -0.2039176723952609E+02;
  COFD[        2213] =   0.4589563967495875E+01;
  COFD[        2214] =  -0.3155845121772854E+00;
  COFD[        2215] =   0.1111636872703185E-01;
  COFD[        2216] =  -0.2009700773592217E+02;
  COFD[        2217] =   0.4972105108130217E+01;
  COFD[        2218] =  -0.4156366087533929E+00;
  COFD[        2219] =   0.1739973132656829E-01;
  COFD[        2220] =  -0.2011340049596478E+02;
  COFD[        2221] =   0.4975835496600746E+01;
  COFD[        2222] =  -0.4162401315456316E+00;
  COFD[        2223] =   0.1743113787867068E-01;
  COFD[        2224] =  -0.2012868226172504E+02;
  COFD[        2225] =   0.4979206995573055E+01;
  COFD[        2226] =  -0.4167863191869534E+00;
  COFD[        2227] =   0.1745958878239685E-01;
  COFD[        2228] =  -0.2141485527132566E+02;
  COFD[        2229] =   0.5180861875662476E+01;
  COFD[        2230] =  -0.4178144385214836E+00;
  COFD[        2231] =   0.1652354263065661E-01;
  COFD[        2232] =  -0.2050229661591078E+02;
  COFD[        2233] =   0.4549054655103792E+01;
  COFD[        2234] =  -0.3083705960322782E+00;
  COFD[        2235] =   0.1072673844162690E-01;
  COFD[        2236] =  -0.2035193547368082E+02;
  COFD[        2237] =   0.4405149516051482E+01;
  COFD[        2238] =  -0.2854877489789628E+00;
  COFD[        2239] =   0.9575581860730991E-02;
  COFD[        2240] =  -0.2130548041809759E+02;
  COFD[        2241] =   0.4924201371605050E+01;
  COFD[        2242] =  -0.3707322745788748E+00;
  COFD[        2243] =   0.1394920260822761E-01;
  COFD[        2244] =  -0.2084204817296711E+02;
  COFD[        2245] =   0.4614418443612200E+01;
  COFD[        2246] =  -0.3198657824249118E+00;
  COFD[        2247] =   0.1134133188499373E-01;
  COFD[        2248] =  -0.2083926067165808E+02;
  COFD[        2249] =   0.4611660894654445E+01;
  COFD[        2250] =  -0.3194768328381069E+00;
  COFD[        2251] =   0.1132352504214021E-01;
  COFD[        2252] =  -0.2126292802833702E+02;
  COFD[        2253] =   0.4875913544147601E+01;
  COFD[        2254] =  -0.3637700066706523E+00;
  COFD[        2255] =   0.1362121204586473E-01;
  COFD[        2256] =  -0.2118199216189331E+02;
  COFD[        2257] =   0.4823957828991331E+01;
  COFD[        2258] =  -0.3562746590271360E+00;
  COFD[        2259] =   0.1326788764265950E-01;
  COFD[        2260] =  -0.2118199216189331E+02;
  COFD[        2261] =   0.4823957828991331E+01;
  COFD[        2262] =  -0.3562746590271360E+00;
  COFD[        2263] =   0.1326788764265950E-01;
  COFD[        2264] =  -0.2116807005556548E+02;
  COFD[        2265] =   0.4815698824842710E+01;
  COFD[        2266] =  -0.3550830672893983E+00;
  COFD[        2267] =   0.1321171158250588E-01;
  COFD[        2268] =  -0.2106221702527866E+02;
  COFD[        2269] =   0.4755729339634630E+01;
  COFD[        2270] =  -0.3464303724989257E+00;
  COFD[        2271] =   0.1280377425307382E-01;
  COFD[        2272] =  -0.2094723649561429E+02;
  COFD[        2273] =   0.4693781529535505E+01;
  COFD[        2274] =  -0.3374918691050086E+00;
  COFD[        2275] =   0.1238235308953681E-01;
  COFD[        2276] =  -0.1994083366515423E+02;
  COFD[        2277] =   0.4933053831289947E+01;
  COFD[        2278] =  -0.4129818889961173E+00;
  COFD[        2279] =   0.1738193532007086E-01;
  COFD[        2280] =  -0.1434856733752978E+02;
  COFD[        2281] =   0.2972259604459018E+01;
  COFD[        2282] =  -0.9821869542826608E-01;
  COFD[        2283] =   0.1198325317137719E-02;
  COFD[        2284] =  -0.1856440584784442E+02;
  COFD[        2285] =   0.5034796976697628E+01;
  COFD[        2286] =  -0.4442640448549661E+00;
  COFD[        2287] =   0.1947279036625017E-01;
  COFD[        2288] =  -0.2017297081385243E+02;
  COFD[        2289] =   0.4879127539439356E+01;
  COFD[        2290] =  -0.3874272290357592E+00;
  COFD[        2291] =   0.1551838396554711E-01;
  COFD[        2292] =  -0.1917251019075018E+02;
  COFD[        2293] =   0.4736217002347487E+01;
  COFD[        2294] =  -0.3839903059289852E+00;
  COFD[        2295] =   0.1596104139772183E-01;
  COFD[        2296] =  -0.2013382472396569E+02;
  COFD[        2297] =   0.4863179495148787E+01;
  COFD[        2298] =  -0.3856326932017908E+00;
  COFD[        2299] =   0.1544705776346803E-01;
  COFD[        2300] =  -0.1920854681460308E+02;
  COFD[        2301] =   0.4744336922244572E+01;
  COFD[        2302] =  -0.3854716437444286E+00;
  COFD[        2303] =   0.1604421161000185E-01;
  COFD[        2304] =  -0.1814291631525713E+02;
  COFD[        2305] =   0.3646821154166707E+01;
  COFD[        2306] =  -0.1727110481572554E+00;
  COFD[        2307] =   0.4220560947533801E-02;
  COFD[        2308] =  -0.2092466715580316E+02;
  COFD[        2309] =   0.4921376438675481E+01;
  COFD[        2310] =  -0.3746202469475856E+00;
  COFD[        2311] =   0.1427858739868015E-01;
  COFD[        2312] =  -0.2012556972106850E+02;
  COFD[        2313] =   0.4927762527501130E+01;
  COFD[        2314] =  -0.4091613023934828E+00;
  COFD[        2315] =   0.1708636787303381E-01;
  COFD[        2316] =  -0.2089621715369740E+02;
  COFD[        2317] =   0.4962017913186494E+01;
  COFD[        2318] =  -0.3839107017164409E+00;
  COFD[        2319] =   0.1483109113820474E-01;
  COFD[        2320] =  -0.2136316488507793E+02;
  COFD[        2321] =   0.5059984882792701E+01;
  COFD[        2322] =  -0.3970696812921098E+00;
  COFD[        2323] =   0.1542543509204860E-01;
  COFD[        2324] =  -0.1947441791465671E+02;
  COFD[        2325] =   0.4068997090767398E+01;
  COFD[        2326] =  -0.2342259510016974E+00;
  COFD[        2327] =   0.7076017538245195E-02;
  COFD[        2328] =  -0.2139377655814223E+02;
  COFD[        2329] =   0.5069003650557725E+01;
  COFD[        2330] =  -0.3983908652693064E+00;
  COFD[        2331] =   0.1548879283162612E-01;
  COFD[        2332] =  -0.2021372325579300E+02;
  COFD[        2333] =   0.4433392392216290E+01;
  COFD[        2334] =  -0.2910494454945646E+00;
  COFD[        2335] =   0.9894154374316580E-02;
  COFD[        2336] =  -0.2038667713911414E+02;
  COFD[        2337] =   0.5013528836767304E+01;
  COFD[        2338] =  -0.4183960149034390E+00;
  COFD[        2339] =   0.1741045522425236E-01;
  COFD[        2340] =  -0.2040768203667939E+02;
  COFD[        2341] =   0.5019349654721370E+01;
  COFD[        2342] =  -0.4193216095365674E+00;
  COFD[        2343] =   0.1745804721026501E-01;
  COFD[        2344] =  -0.2042738004900579E+02;
  COFD[        2345] =   0.5024724404363650E+01;
  COFD[        2346] =  -0.4201768197596846E+00;
  COFD[        2347] =   0.1750204200208273E-01;
  COFD[        2348] =  -0.2148234754894915E+02;
  COFD[        2349] =   0.5128240108050982E+01;
  COFD[        2350] =  -0.4075464684376063E+00;
  COFD[        2351] =   0.1594434352315902E-01;
  COFD[        2352] =  -0.2035193547368082E+02;
  COFD[        2353] =   0.4405149516051482E+01;
  COFD[        2354] =  -0.2854877489789628E+00;
  COFD[        2355] =   0.9575581860730991E-02;
  COFD[        2356] =  -0.2015827532495969E+02;
  COFD[        2357] =   0.4245811407071466E+01;
  COFD[        2358] =  -0.2605338071793272E+00;
  COFD[        2359] =   0.8331438357037944E-02;
  COFD[        2360] =  -0.2137644847507669E+02;
  COFD[        2361] =   0.4879361403534794E+01;
  COFD[        2362] =  -0.3619130719444060E+00;
  COFD[        2363] =   0.1345545271628123E-01;
  COFD[        2364] =  -0.2080211055309488E+02;
  COFD[        2365] =   0.4522132770796719E+01;
  COFD[        2366] =  -0.3043845860308388E+00;
  COFD[        2367] =   0.1053714557960114E-01;
  COFD[        2368] =  -0.2080221542804032E+02;
  COFD[        2369] =   0.4520604274531649E+01;
  COFD[        2370] =  -0.3041726314312171E+00;
  COFD[        2371] =   0.1052767250427024E-01;
  COFD[        2372] =  -0.2139170339962327E+02;
  COFD[        2373] =   0.4855952270590682E+01;
  COFD[        2374] =  -0.3585715335315449E+00;
  COFD[        2375] =   0.1330006698014156E-01;
  COFD[        2376] =  -0.2135235844472485E+02;
  COFD[        2377] =   0.4822077076791630E+01;
  COFD[        2378] =  -0.3537269179793018E+00;
  COFD[        2379] =   0.1307429400877126E-01;
  COFD[        2380] =  -0.2135235844472485E+02;
  COFD[        2381] =   0.4822077076791630E+01;
  COFD[        2382] =  -0.3537269179793018E+00;
  COFD[        2383] =   0.1307429400877126E-01;
  COFD[        2384] =  -0.2134436611864803E+02;
  COFD[        2385] =   0.4816403941921898E+01;
  COFD[        2386] =  -0.3529153872118850E+00;
  COFD[        2387] =   0.1303646422215331E-01;
  COFD[        2388] =  -0.2127828879489648E+02;
  COFD[        2389] =   0.4773820066846246E+01;
  COFD[        2390] =  -0.3468230916877841E+00;
  COFD[        2391] =   0.1275243316292433E-01;
  COFD[        2392] =  -0.2120031254217025E+02;
  COFD[        2393] =   0.4728097335132442E+01;
  COFD[        2394] =  -0.3402810017723648E+00;
  COFD[        2395] =   0.1244740133248284E-01;
  COFD[        2396] =  -0.2013110520151531E+02;
  COFD[        2397] =   0.4933748223113413E+01;
  COFD[        2398] =  -0.4099898956624159E+00;
  COFD[        2399] =   0.1712216655141095E-01;
  COFD[        2400] =  -0.1569651926124499E+02;
  COFD[        2401] =   0.3625575306654659E+01;
  COFD[        2402] =  -0.1996382427440662E+00;
  COFD[        2403] =   0.6236602581187985E-02;
  COFD[        2404] =  -0.1779799510057660E+02;
  COFD[        2405] =   0.4774820362539309E+01;
  COFD[        2406] =  -0.4169337844744687E+00;
  COFD[        2407] =   0.1856826481796923E-01;
  COFD[        2408] =  -0.1979190800756452E+02;
  COFD[        2409] =   0.4812635185197916E+01;
  COFD[        2410] =  -0.3895191141694820E+00;
  COFD[        2411] =   0.1603137873834851E-01;
  COFD[        2412] =  -0.1891691641938918E+02;
  COFD[        2413] =   0.4708474963146238E+01;
  COFD[        2414] =  -0.3908198429896564E+00;
  COFD[        2415] =   0.1669173693263477E-01;
  COFD[        2416] =  -0.1978759220318865E+02;
  COFD[        2417] =   0.4813522654682605E+01;
  COFD[        2418] =  -0.3904695959382222E+00;
  COFD[        2419] =   0.1610776277485814E-01;
  COFD[        2420] =  -0.1891630498772998E+02;
  COFD[        2421] =   0.4700565744006908E+01;
  COFD[        2422] =  -0.3899919372649346E+00;
  COFD[        2423] =   0.1666500811895143E-01;
  COFD[        2424] =  -0.1920855889011546E+02;
  COFD[        2425] =   0.4192213485908265E+01;
  COFD[        2426] =  -0.2588520836689803E+00;
  COFD[        2427] =   0.8501867421401927E-02;
  COFD[        2428] =  -0.2097636580694235E+02;
  COFD[        2429] =   0.5051780969690284E+01;
  COFD[        2430] =  -0.4052813410201178E+00;
  COFD[        2431] =   0.1613327373268888E-01;
  COFD[        2432] =  -0.1954575230865094E+02;
  COFD[        2433] =   0.4765735338615951E+01;
  COFD[        2434] =  -0.3970172796487949E+00;
  COFD[        2435] =   0.1692142125782123E-01;
  COFD[        2436] =  -0.2090423070644263E+02;
  COFD[        2437] =   0.5069554235969105E+01;
  COFD[        2438] =  -0.4110387472810522E+00;
  COFD[        2439] =   0.1651243898718931E-01;
  COFD[        2440] =  -0.2137421544016685E+02;
  COFD[        2441] =   0.5174419154841422E+01;
  COFD[        2442] =  -0.4256066415182231E+00;
  COFD[        2443] =   0.1719123893975486E-01;
  COFD[        2444] =  -0.2074011798729127E+02;
  COFD[        2445] =   0.4716308251797323E+01;
  COFD[        2446] =  -0.3368455089670934E+00;
  COFD[        2447] =   0.1222481812086106E-01;
  COFD[        2448] =  -0.2140524845935924E+02;
  COFD[        2449] =   0.5183807431960805E+01;
  COFD[        2450] =  -0.4270095585158289E+00;
  COFD[        2451] =   0.1725985135601222E-01;
  COFD[        2452] =  -0.2108918057162479E+02;
  COFD[        2453] =   0.4917141448841289E+01;
  COFD[        2454] =  -0.3711916302984747E+00;
  COFD[        2455] =   0.1402561629509121E-01;
  COFD[        2456] =  -0.1967283644204365E+02;
  COFD[        2457] =   0.4798018532274038E+01;
  COFD[        2458] =  -0.3991636476917164E+00;
  COFD[        2459] =   0.1693109854880472E-01;
  COFD[        2460] =  -0.1968298220975192E+02;
  COFD[        2461] =   0.4799139981296912E+01;
  COFD[        2462] =  -0.3994189723094190E+00;
  COFD[        2463] =   0.1694708582407344E-01;
  COFD[        2464] =  -0.1969261804237637E+02;
  COFD[        2465] =   0.4800160734138452E+01;
  COFD[        2466] =  -0.3996532469352685E+00;
  COFD[        2467] =   0.1696180432302862E-01;
  COFD[        2468] =  -0.2148607955277492E+02;
  COFD[        2469] =   0.5238740013363649E+01;
  COFD[        2470] =  -0.4356664393291084E+00;
  COFD[        2471] =   0.1769834275551969E-01;
  COFD[        2472] =  -0.2130548041809759E+02;
  COFD[        2473] =   0.4924201371605050E+01;
  COFD[        2474] =  -0.3707322745788748E+00;
  COFD[        2475] =   0.1394920260822761E-01;
  COFD[        2476] =  -0.2137644847507669E+02;
  COFD[        2477] =   0.4879361403534794E+01;
  COFD[        2478] =  -0.3619130719444060E+00;
  COFD[        2479] =   0.1345545271628123E-01;
  COFD[        2480] =  -0.2176129979912939E+02;
  COFD[        2481] =   0.5158607600017159E+01;
  COFD[        2482] =  -0.4139137284684214E+00;
  COFD[        2483] =   0.1631382826546863E-01;
  COFD[        2484] =  -0.2163570087754292E+02;
  COFD[        2485] =   0.4991338208642202E+01;
  COFD[        2486] =  -0.3826379744946652E+00;
  COFD[        2487] =   0.1458812354473534E-01;
  COFD[        2488] =  -0.2163339664854972E+02;
  COFD[        2489] =   0.4988686988670290E+01;
  COFD[        2490] =  -0.3822530252797285E+00;
  COFD[        2491] =   0.1456984260990830E-01;
  COFD[        2492] =  -0.2176867529242300E+02;
  COFD[        2493] =   0.5131070812581271E+01;
  COFD[        2494] =  -0.4098687890615916E+00;
  COFD[        2495] =   0.1611928094862400E-01;
  COFD[        2496] =  -0.2171041726819223E+02;
  COFD[        2497] =   0.5087936731469250E+01;
  COFD[        2498] =  -0.4035328710130172E+00;
  COFD[        2499] =   0.1581455637665996E-01;
  COFD[        2500] =  -0.2171041726819223E+02;
  COFD[        2501] =   0.5087936731469250E+01;
  COFD[        2502] =  -0.4035328710130172E+00;
  COFD[        2503] =   0.1581455637665996E-01;
  COFD[        2504] =  -0.2169907962284337E+02;
  COFD[        2505] =   0.5080635121967584E+01;
  COFD[        2506] =  -0.4024603628256680E+00;
  COFD[        2507] =   0.1576297563539379E-01;
  COFD[        2508] =  -0.2160710054610370E+02;
  COFD[        2509] =   0.5025481182776788E+01;
  COFD[        2510] =  -0.3943591536682922E+00;
  COFD[        2511] =   0.1537337185442625E-01;
  COFD[        2512] =  -0.2150038725326603E+02;
  COFD[        2513] =   0.4965856294754335E+01;
  COFD[        2514] =  -0.3856015603443829E+00;
  COFD[        2515] =   0.1495222567994429E-01;
  COFD[        2516] =  -0.1952277059880253E+02;
  COFD[        2517] =   0.4760388219539090E+01;
  COFD[        2518] =  -0.3964022009449412E+00;
  COFD[        2519] =   0.1689791166807166E-01;
  COFD[        2520] =  -0.1433143401213129E+02;
  COFD[        2521] =   0.2944305087428476E+01;
  COFD[        2522] =  -0.9451478525573859E-01;
  COFD[        2523] =   0.9946176585968385E-03;
  COFD[        2524] =  -0.1866728292989077E+02;
  COFD[        2525] =   0.5078383565507850E+01;
  COFD[        2526] =  -0.4538464035977413E+00;
  COFD[        2527] =   0.2005643240392637E-01;
  COFD[        2528] =  -0.1984047911980011E+02;
  COFD[        2529] =   0.4733570408392962E+01;
  COFD[        2530] =  -0.3701539096262219E+00;
  COFD[        2531] =   0.1481934816743204E-01;
  COFD[        2532] =  -0.1935099040011764E+02;
  COFD[        2533] =   0.4796859780440744E+01;
  COFD[        2534] =  -0.3951546730488872E+00;
  COFD[        2535] =   0.1657370960389658E-01;
  COFD[        2536] =  -0.1988931959990085E+02;
  COFD[        2537] =   0.4756594033154599E+01;
  COFD[        2538] =  -0.3742321497803292E+00;
  COFD[        2539] =   0.1504181447689004E-01;
  COFD[        2540] =  -0.1936781090512805E+02;
  COFD[        2541] =   0.4796060596924874E+01;
  COFD[        2542] =  -0.3953630801104229E+00;
  COFD[        2543] =   0.1659683769348035E-01;
  COFD[        2544] =  -0.1874280769185066E+02;
  COFD[        2545] =   0.3905545522430255E+01;
  COFD[        2546] =  -0.2138630507050396E+00;
  COFD[        2547] =   0.6308839293777503E-02;
  COFD[        2548] =  -0.2080529541825178E+02;
  COFD[        2549] =   0.4870282456696042E+01;
  COFD[        2550] =  -0.3717177531779383E+00;
  COFD[        2551] =   0.1428959829981116E-01;
  COFD[        2552] =  -0.1997997031335933E+02;
  COFD[        2553] =   0.4852495827753419E+01;
  COFD[        2554] =  -0.4017712513240834E+00;
  COFD[        2555] =   0.1685558673482876E-01;
  COFD[        2556] =  -0.2072976841817859E+02;
  COFD[        2557] =   0.4886126782915337E+01;
  COFD[        2558] =  -0.3771688345094954E+00;
  COFD[        2559] =   0.1464949867829970E-01;
  COFD[        2560] =  -0.2130885185428343E+02;
  COFD[        2561] =   0.5038362567088167E+01;
  COFD[        2562] =  -0.3984885360207525E+00;
  COFD[        2563] =   0.1564530982215940E-01;
  COFD[        2564] =  -0.2009226460847481E+02;
  COFD[        2565] =   0.4335080602770606E+01;
  COFD[        2566] =  -0.2764812336210347E+00;
  COFD[        2567] =   0.9208446946513201E-02;
  COFD[        2568] =  -0.2134888342401024E+02;
  COFD[        2569] =   0.5051235704256541E+01;
  COFD[        2570] =  -0.4003883872653834E+00;
  COFD[        2571] =   0.1573711341614223E-01;
  COFD[        2572] =  -0.2056103945202781E+02;
  COFD[        2573] =   0.4582953188069841E+01;
  COFD[        2574] =  -0.3167741504985752E+00;
  COFD[        2575] =   0.1125103765063593E-01;
  COFD[        2576] =  -0.2012864081733589E+02;
  COFD[        2577] =   0.4891098053091741E+01;
  COFD[        2578] =  -0.4047491719429138E+00;
  COFD[        2579] =   0.1690406746445060E-01;
  COFD[        2580] =  -0.2014955127994030E+02;
  COFD[        2581] =   0.4896710570985273E+01;
  COFD[        2582] =  -0.4056774618307388E+00;
  COFD[        2583] =   0.1695311601042129E-01;
  COFD[        2584] =  -0.2016949900503337E+02;
  COFD[        2585] =   0.4902039892776266E+01;
  COFD[        2586] =  -0.4065592769263701E+00;
  COFD[        2587] =   0.1699972276094747E-01;
  COFD[        2588] =  -0.2152239257799227E+02;
  COFD[        2589] =   0.5141731472937461E+01;
  COFD[        2590] =  -0.4142263111871738E+00;
  COFD[        2591] =   0.1642254345298395E-01;
  COFD[        2592] =  -0.2084204817296711E+02;
  COFD[        2593] =   0.4614418443612200E+01;
  COFD[        2594] =  -0.3198657824249118E+00;
  COFD[        2595] =   0.1134133188499373E-01;
  COFD[        2596] =  -0.2080211055309488E+02;
  COFD[        2597] =   0.4522132770796719E+01;
  COFD[        2598] =  -0.3043845860308388E+00;
  COFD[        2599] =   0.1053714557960114E-01;
  COFD[        2600] =  -0.2163570087754292E+02;
  COFD[        2601] =   0.4991338208642202E+01;
  COFD[        2602] =  -0.3826379744946652E+00;
  COFD[        2603] =   0.1458812354473534E-01;
  COFD[        2604] =  -0.2142272784698728E+02;
  COFD[        2605] =   0.4785909395121438E+01;
  COFD[        2606] =  -0.3464660102393260E+00;
  COFD[        2607] =   0.1265829609478411E-01;
  COFD[        2608] =  -0.2142679763402700E+02;
  COFD[        2609] =   0.4785846726597538E+01;
  COFD[        2610] =  -0.3464570913102326E+00;
  COFD[        2611] =   0.1265788303675542E-01;
  COFD[        2612] =  -0.2176051031769040E+02;
  COFD[        2613] =   0.5011965066548862E+01;
  COFD[        2614] =  -0.3856361645401574E+00;
  COFD[        2615] =   0.1473067061515274E-01;
  COFD[        2616] =  -0.2177389401397997E+02;
  COFD[        2617] =   0.4998448046198806E+01;
  COFD[        2618] =  -0.3836787642667724E+00;
  COFD[        2619] =   0.1463798128403765E-01;
  COFD[        2620] =  -0.2177389401397997E+02;
  COFD[        2621] =   0.4998448046198806E+01;
  COFD[        2622] =  -0.3836787642667724E+00;
  COFD[        2623] =   0.1463798128403765E-01;
  COFD[        2624] =  -0.2177183976325961E+02;
  COFD[        2625] =   0.4994984064373654E+01;
  COFD[        2626] =  -0.3831767216659396E+00;
  COFD[        2627] =   0.1461418638286886E-01;
  COFD[        2628] =  -0.2173580178381468E+02;
  COFD[        2629] =   0.4962896685693639E+01;
  COFD[        2630] =  -0.3785249617010963E+00;
  COFD[        2631] =   0.1439364791769246E-01;
  COFD[        2632] =  -0.2167081740518373E+02;
  COFD[        2633] =   0.4920315656610881E+01;
  COFD[        2634] =  -0.3723507532504461E+00;
  COFD[        2635] =   0.1410087595979374E-01;
  COFD[        2636] =  -0.1995821009052950E+02;
  COFD[        2637] =   0.4847729338512834E+01;
  COFD[        2638] =  -0.4012454442759492E+00;
  COFD[        2639] =   0.1683646904857131E-01;
  COFD[        2640] =  -0.1432407507363763E+02;
  COFD[        2641] =   0.2940787058428872E+01;
  COFD[        2642] =  -0.9397416078466805E-01;
  COFD[        2643] =   0.9675574576547308E-03;
  COFD[        2644] =  -0.1867240731108900E+02;
  COFD[        2645] =   0.5080367429995977E+01;
  COFD[        2646] =  -0.4541125321914433E+00;
  COFD[        2647] =   0.2006831547910542E-01;
  COFD[        2648] =  -0.1983092143620168E+02;
  COFD[        2649] =   0.4728281143697850E+01;
  COFD[        2650] =  -0.3693402228450711E+00;
  COFD[        2651] =   0.1477857564328594E-01;
  COFD[        2652] =  -0.1935470245750297E+02;
  COFD[        2653] =   0.4797074097380636E+01;
  COFD[        2654] =  -0.3950970944572844E+00;
  COFD[        2655] =   0.1656735836140260E-01;
  COFD[        2656] =  -0.1987978855175104E+02;
  COFD[        2657] =   0.4751262797957020E+01;
  COFD[        2658] =  -0.3734104385798870E+00;
  COFD[        2659] =   0.1500058297946865E-01;
  COFD[        2660] =  -0.1937166591561948E+02;
  COFD[        2661] =   0.4796286056188764E+01;
  COFD[        2662] =  -0.3953052681051102E+00;
  COFD[        2663] =   0.1659039955935203E-01;
  COFD[        2664] =  -0.1873753222748401E+02;
  COFD[        2665] =   0.3902621516519401E+01;
  COFD[        2666] =  -0.2134883091855596E+00;
  COFD[        2667] =   0.6294056333331442E-02;
  COFD[        2668] =  -0.2079012230483711E+02;
  COFD[        2669] =   0.4862345324915401E+01;
  COFD[        2670] =  -0.3705511888834099E+00;
  COFD[        2671] =   0.1423346818325797E-01;
  COFD[        2672] =  -0.1997892794374273E+02;
  COFD[        2673] =   0.4850385203230214E+01;
  COFD[        2674] =  -0.4014003482747929E+00;
  COFD[        2675] =   0.1683520757557910E-01;
  COFD[        2676] =  -0.2071293986412085E+02;
  COFD[        2677] =   0.4877335837263783E+01;
  COFD[        2678] =  -0.3758684263056035E+00;
  COFD[        2679] =   0.1458649350932613E-01;
  COFD[        2680] =  -0.2129718866568306E+02;
  COFD[        2681] =   0.5031879745702688E+01;
  COFD[        2682] =  -0.3975318057376973E+00;
  COFD[        2683] =   0.1559908037190491E-01;
  COFD[        2684] =  -0.2008729443847752E+02;
  COFD[        2685] =   0.4331784646722137E+01;
  COFD[        2686] =  -0.2760303040339949E+00;
  COFD[        2687] =   0.9188681380243065E-02;
  COFD[        2688] =  -0.2133756080522132E+02;
  COFD[        2689] =   0.5044876892625164E+01;
  COFD[        2690] =  -0.3994499299174794E+00;
  COFD[        2691] =   0.1569176539047349E-01;
  COFD[        2692] =  -0.2055209768497027E+02;
  COFD[        2693] =   0.4577812802798165E+01;
  COFD[        2694] =  -0.3160440880823138E+00;
  COFD[        2695] =   0.1121730347594009E-01;
  COFD[        2696] =  -0.2012555233929497E+02;
  COFD[        2697] =   0.4887988169322795E+01;
  COFD[        2698] =  -0.4042349632777354E+00;
  COFD[        2699] =   0.1687690382084792E-01;
  COFD[        2700] =  -0.2014664412250720E+02;
  COFD[        2701] =   0.4893661499235972E+01;
  COFD[        2702] =  -0.4051731103333032E+00;
  COFD[        2703] =   0.1692646556995869E-01;
  COFD[        2704] =  -0.2016677937256390E+02;
  COFD[        2705] =   0.4899055461722977E+01;
  COFD[        2706] =  -0.4060654125033917E+00;
  COFD[        2707] =   0.1697361869912191E-01;
  COFD[        2708] =  -0.2151480952421122E+02;
  COFD[        2709] =   0.5136706180179550E+01;
  COFD[        2710] =  -0.4134839036187021E+00;
  COFD[        2711] =   0.1638663893119980E-01;
  COFD[        2712] =  -0.2083926067165808E+02;
  COFD[        2713] =   0.4611660894654445E+01;
  COFD[        2714] =  -0.3194768328381069E+00;
  COFD[        2715] =   0.1132352504214021E-01;
  COFD[        2716] =  -0.2080221542804032E+02;
  COFD[        2717] =   0.4520604274531649E+01;
  COFD[        2718] =  -0.3041726314312171E+00;
  COFD[        2719] =   0.1052767250427024E-01;
  COFD[        2720] =  -0.2163339664854972E+02;
  COFD[        2721] =   0.4988686988670290E+01;
  COFD[        2722] =  -0.3822530252797285E+00;
  COFD[        2723] =   0.1456984260990830E-01;
  COFD[        2724] =  -0.2142679763402700E+02;
  COFD[        2725] =   0.4785846726597538E+01;
  COFD[        2726] =  -0.3464570913102326E+00;
  COFD[        2727] =   0.1265788303675542E-01;
  COFD[        2728] =  -0.2143119141574728E+02;
  COFD[        2729] =   0.4785909395121408E+01;
  COFD[        2730] =  -0.3464660102393218E+00;
  COFD[        2731] =   0.1265829609478391E-01;
  COFD[        2732] =  -0.2176529489902795E+02;
  COFD[        2733] =   0.5012168456615743E+01;
  COFD[        2734] =  -0.3856653477006033E+00;
  COFD[        2735] =   0.1473203875352443E-01;
  COFD[        2736] =  -0.2178318492664675E+02;
  COFD[        2737] =   0.5000465951311985E+01;
  COFD[        2738] =  -0.3839711955741957E+00;
  COFD[        2739] =   0.1465183998943084E-01;
  COFD[        2740] =  -0.2178318492664675E+02;
  COFD[        2741] =   0.5000465951311985E+01;
  COFD[        2742] =  -0.3839711955741957E+00;
  COFD[        2743] =   0.1465183998943084E-01;
  COFD[        2744] =  -0.2178172115101771E+02;
  COFD[        2745] =   0.4997238443372161E+01;
  COFD[        2746] =  -0.3835034599913829E+00;
  COFD[        2747] =   0.1462967284123923E-01;
  COFD[        2748] =  -0.2174924769698764E+02;
  COFD[        2749] =   0.4966564768133120E+01;
  COFD[        2750] =  -0.3790567913063227E+00;
  COFD[        2751] =   0.1441886477883868E-01;
  COFD[        2752] =  -0.2168687631153109E+02;
  COFD[        2753] =   0.4924990034391955E+01;
  COFD[        2754] =  -0.3730285599341068E+00;
  COFD[        2755] =   0.1413301763976827E-01;
  COFD[        2756] =  -0.1995722688303096E+02;
  COFD[        2757] =   0.4845641046786140E+01;
  COFD[        2758] =  -0.4008772701599357E+00;
  COFD[        2759] =   0.1681620052332664E-01;
  COFD[        2760] =  -0.1561208690963572E+02;
  COFD[        2761] =   0.3582789091852060E+01;
  COFD[        2762] =  -0.1928143728969187E+00;
  COFD[        2763] =   0.5885365651518750E-02;
  COFD[        2764] =  -0.1789359789035752E+02;
  COFD[        2765] =   0.4812277632492618E+01;
  COFD[        2766] =  -0.4220338205719207E+00;
  COFD[        2767] =   0.1879929285184377E-01;
  COFD[        2768] =  -0.1969361932498264E+02;
  COFD[        2769] =   0.4748312088352801E+01;
  COFD[        2770] =  -0.3792560078874427E+00;
  COFD[        2771] =   0.1550280308751674E-01;
  COFD[        2772] =  -0.1905867186145966E+02;
  COFD[        2773] =   0.4745150284504969E+01;
  COFD[        2774] =  -0.3946441379637312E+00;
  COFD[        2775] =   0.1681441351814178E-01;
  COFD[        2776] =  -0.1969237276879183E+02;
  COFD[        2777] =   0.4749599566916338E+01;
  COFD[        2778] =  -0.3802354707820949E+00;
  COFD[        2779] =   0.1557944508694219E-01;
  COFD[        2780] =  -0.1906167945299586E+02;
  COFD[        2781] =   0.4737944806110550E+01;
  COFD[        2782] =  -0.3938946014462452E+00;
  COFD[        2783] =   0.1679046692522489E-01;
  COFD[        2784] =  -0.1898993862716132E+02;
  COFD[        2785] =   0.4082522815084645E+01;
  COFD[        2786] =  -0.2435514237071018E+00;
  COFD[        2787] =   0.7811876908468279E-02;
  COFD[        2788] =  -0.2073725116067806E+02;
  COFD[        2789] =   0.4921491851668374E+01;
  COFD[        2790] =  -0.3859047434126549E+00;
  COFD[        2791] =   0.1518991301474521E-01;
  COFD[        2792] =  -0.1962170244526679E+02;
  COFD[        2793] =   0.4768861169263409E+01;
  COFD[        2794] =  -0.3963823359686274E+00;
  COFD[        2795] =   0.1684681197314066E-01;
  COFD[        2796] =  -0.2065920762233094E+02;
  COFD[        2797] =   0.4934466495930646E+01;
  COFD[        2798] =  -0.3907995086461343E+00;
  COFD[        2799] =   0.1551999187366136E-01;
  COFD[        2800] =  -0.2120949634659705E+02;
  COFD[        2801] =   0.5075665220886366E+01;
  COFD[        2802] =  -0.4108520560652987E+00;
  COFD[        2803] =   0.1646976074078464E-01;
  COFD[        2804] =  -0.2060180540567790E+02;
  COFD[        2805] =   0.4632891986456493E+01;
  COFD[        2806] =  -0.3249863297154995E+00;
  COFD[        2807] =   0.1167639716644466E-01;
  COFD[        2808] =  -0.2124842946919059E+02;
  COFD[        2809] =   0.5088113178440255E+01;
  COFD[        2810] =  -0.4127116950602789E+00;
  COFD[        2811] =   0.1656068511393674E-01;
  COFD[        2812] =  -0.2092002080559684E+02;
  COFD[        2813] =   0.4818457997428789E+01;
  COFD[        2814] =  -0.3569233511421860E+00;
  COFD[        2815] =   0.1335173475837706E-01;
  COFD[        2816] =  -0.1971435410948289E+02;
  COFD[        2817] =   0.4784277383958439E+01;
  COFD[        2818] =  -0.3961230823099747E+00;
  COFD[        2819] =   0.1674302776354957E-01;
  COFD[        2820] =  -0.1972631705838112E+02;
  COFD[        2821] =   0.4785911211303439E+01;
  COFD[        2822] =  -0.3964794306558241E+00;
  COFD[        2823] =   0.1676492861574991E-01;
  COFD[        2824] =  -0.1973776103143215E+02;
  COFD[        2825] =   0.4787458937602214E+01;
  COFD[        2826] =  -0.3968179100773139E+00;
  COFD[        2827] =   0.1678575643964539E-01;
  COFD[        2828] =  -0.2141313429386006E+02;
  COFD[        2829] =   0.5174902392506038E+01;
  COFD[        2830] =  -0.4261176538552951E+00;
  COFD[        2831] =   0.1723100350853119E-01;
  COFD[        2832] =  -0.2126292802833702E+02;
  COFD[        2833] =   0.4875913544147601E+01;
  COFD[        2834] =  -0.3637700066706523E+00;
  COFD[        2835] =   0.1362121204586473E-01;
  COFD[        2836] =  -0.2139170339962327E+02;
  COFD[        2837] =   0.4855952270590682E+01;
  COFD[        2838] =  -0.3585715335315449E+00;
  COFD[        2839] =   0.1330006698014156E-01;
  COFD[        2840] =  -0.2176867529242300E+02;
  COFD[        2841] =   0.5131070812581271E+01;
  COFD[        2842] =  -0.4098687890615916E+00;
  COFD[        2843] =   0.1611928094862400E-01;
  COFD[        2844] =  -0.2176051031769040E+02;
  COFD[        2845] =   0.5011965066548862E+01;
  COFD[        2846] =  -0.3856361645401574E+00;
  COFD[        2847] =   0.1473067061515274E-01;
  COFD[        2848] =  -0.2176529489902795E+02;
  COFD[        2849] =   0.5012168456615743E+01;
  COFD[        2850] =  -0.3856653477006033E+00;
  COFD[        2851] =   0.1473203875352443E-01;
  COFD[        2852] =  -0.2191034670855553E+02;
  COFD[        2853] =   0.5158607600017168E+01;
  COFD[        2854] =  -0.4139137284684229E+00;
  COFD[        2855] =   0.1631382826546871E-01;
  COFD[        2856] =  -0.2193017488088049E+02;
  COFD[        2857] =   0.5147264054224069E+01;
  COFD[        2858] =  -0.4122474414807874E+00;
  COFD[        2859] =   0.1623368509412096E-01;
  COFD[        2860] =  -0.2193017488088049E+02;
  COFD[        2861] =   0.5147264054224069E+01;
  COFD[        2862] =  -0.4122474414807874E+00;
  COFD[        2863] =   0.1623368509412096E-01;
  COFD[        2864] =  -0.2192865333155479E+02;
  COFD[        2865] =   0.5143921941237063E+01;
  COFD[        2866] =  -0.4117565109010827E+00;
  COFD[        2867] =   0.1621007304720648E-01;
  COFD[        2868] =  -0.2189359568671224E+02;
  COFD[        2869] =   0.5111398712052350E+01;
  COFD[        2870] =  -0.4069791541191179E+00;
  COFD[        2871] =   0.1598030296830180E-01;
  COFD[        2872] =  -0.2182529370017760E+02;
  COFD[        2873] =   0.5066380027973548E+01;
  COFD[        2874] =  -0.4003664957816425E+00;
  COFD[        2875] =   0.1566227519627694E-01;
  COFD[        2876] =  -0.1959997385926940E+02;
  COFD[        2877] =   0.4764004414585760E+01;
  COFD[        2878] =  -0.3958298367826015E+00;
  COFD[        2879] =   0.1682596038127772E-01;
  COFD[        2880] =  -0.1556317061681839E+02;
  COFD[        2881] =   0.3558114950801765E+01;
  COFD[        2882] =  -0.1888799785819696E+00;
  COFD[        2883] =   0.5682884908013199E-02;
  COFD[        2884] =  -0.1794643219743921E+02;
  COFD[        2885] =   0.4833075314023510E+01;
  COFD[        2886] =  -0.4248631135846847E+00;
  COFD[        2887] =   0.1892735769307228E-01;
  COFD[        2888] =  -0.1963494149489594E+02;
  COFD[        2889] =   0.4711113691258022E+01;
  COFD[        2890] =  -0.3733233453642334E+00;
  COFD[        2891] =   0.1519736886874179E-01;
  COFD[        2892] =  -0.1913891562399630E+02;
  COFD[        2893] =   0.4766338982641832E+01;
  COFD[        2894] =  -0.3968432892885747E+00;
  COFD[        2895] =   0.1688440639323740E-01;
  COFD[        2896] =  -0.1963449815268955E+02;
  COFD[        2897] =   0.4712178662790742E+01;
  COFD[        2898] =  -0.3742472943897919E+00;
  COFD[        2899] =   0.1527044453467252E-01;
  COFD[        2900] =  -0.1914452665813124E+02;
  COFD[        2901] =   0.4759724497948794E+01;
  COFD[        2902] =  -0.3961575775878822E+00;
  COFD[        2903] =   0.1686262510597603E-01;
  COFD[        2904] =  -0.1885646323568239E+02;
  COFD[        2905] =   0.4016884863488960E+01;
  COFD[        2906] =  -0.2343948477217928E+00;
  COFD[        2907] =   0.7398933738679844E-02;
  COFD[        2908] =  -0.2057011439802693E+02;
  COFD[        2909] =   0.4834011551777849E+01;
  COFD[        2910] =  -0.3728956950726914E+00;
  COFD[        2911] =   0.1455662572760375E-01;
  COFD[        2912] =  -0.1966500520249660E+02;
  COFD[        2913] =   0.4770450920673444E+01;
  COFD[        2914] =  -0.3958664514070914E+00;
  COFD[        2915] =   0.1679182031481241E-01;
  COFD[        2916] =  -0.2048311817894017E+02;
  COFD[        2917] =   0.4841596720970757E+01;
  COFD[        2918] =  -0.3768868905348668E+00;
  COFD[        2919] =   0.1483785654563742E-01;
  COFD[        2920] =  -0.2107933630762049E+02;
  COFD[        2921] =   0.5003615183707528E+01;
  COFD[        2922] =  -0.4000889881542739E+00;
  COFD[        2923] =   0.1594355520326905E-01;
  COFD[        2924] =  -0.2049312404884751E+02;
  COFD[        2925] =   0.4573043093616539E+01;
  COFD[        2926] =  -0.3164763709052368E+00;
  COFD[        2927] =   0.1128280291029378E-01;
  COFD[        2928] =  -0.2112052652732195E+02;
  COFD[        2929] =   0.5016777693101416E+01;
  COFD[        2930] =  -0.4020551643132008E+00;
  COFD[        2931] =   0.1603967677731761E-01;
  COFD[        2932] =  -0.2078560813380206E+02;
  COFD[        2933] =   0.4746137595541170E+01;
  COFD[        2934] =  -0.3464661176921034E+00;
  COFD[        2935] =   0.1285782723518177E-01;
  COFD[        2936] =  -0.1973273884329878E+02;
  COFD[        2937] =   0.4773632942733207E+01;
  COFD[        2938] =  -0.3938198339615520E+00;
  COFD[        2939] =   0.1660198213420950E-01;
  COFD[        2940] =  -0.1974560293529215E+02;
  COFD[        2941] =   0.4775430480336665E+01;
  COFD[        2942] =  -0.3942069403173756E+00;
  COFD[        2943] =   0.1662563585220353E-01;
  COFD[        2944] =  -0.1975795254936833E+02;
  COFD[        2945] =   0.4777156336087093E+01;
  COFD[        2946] =  -0.3945792342207480E+00;
  COFD[        2947] =   0.1664840211599068E-01;
  COFD[        2948] =  -0.2131043155757551E+02;
  COFD[        2949] =   0.5111234850420700E+01;
  COFD[        2950] =  -0.4165974605140776E+00;
  COFD[        2951] =   0.1676520530944528E-01;
  COFD[        2952] =  -0.2118199216189331E+02;
  COFD[        2953] =   0.4823957828991331E+01;
  COFD[        2954] =  -0.3562746590271360E+00;
  COFD[        2955] =   0.1326788764265950E-01;
  COFD[        2956] =  -0.2135235844472485E+02;
  COFD[        2957] =   0.4822077076791630E+01;
  COFD[        2958] =  -0.3537269179793018E+00;
  COFD[        2959] =   0.1307429400877126E-01;
  COFD[        2960] =  -0.2171041726819223E+02;
  COFD[        2961] =   0.5087936731469250E+01;
  COFD[        2962] =  -0.4035328710130172E+00;
  COFD[        2963] =   0.1581455637665996E-01;
  COFD[        2964] =  -0.2177389401397997E+02;
  COFD[        2965] =   0.4998448046198806E+01;
  COFD[        2966] =  -0.3836787642667724E+00;
  COFD[        2967] =   0.1463798128403765E-01;
  COFD[        2968] =  -0.2178318492664675E+02;
  COFD[        2969] =   0.5000465951311985E+01;
  COFD[        2970] =  -0.3839711955741957E+00;
  COFD[        2971] =   0.1465183998943084E-01;
  COFD[        2972] =  -0.2193017488088049E+02;
  COFD[        2973] =   0.5147264054224069E+01;
  COFD[        2974] =  -0.4122474414807874E+00;
  COFD[        2975] =   0.1623368509412096E-01;
  COFD[        2976] =  -0.2200526405431778E+02;
  COFD[        2977] =   0.5158607600017140E+01;
  COFD[        2978] =  -0.4139137284684187E+00;
  COFD[        2979] =   0.1631382826546851E-01;
  COFD[        2980] =  -0.2200526405431778E+02;
  COFD[        2981] =   0.5158607600017140E+01;
  COFD[        2982] =  -0.4139137284684187E+00;
  COFD[        2983] =   0.1631382826546851E-01;
  COFD[        2984] =  -0.2201134713900483E+02;
  COFD[        2985] =   0.5158384050126705E+01;
  COFD[        2986] =  -0.4138808904484991E+00;
  COFD[        2987] =   0.1631224885111492E-01;
  COFD[        2988] =  -0.2202438945187234E+02;
  COFD[        2989] =   0.5145495392374555E+01;
  COFD[        2990] =  -0.4119876385887463E+00;
  COFD[        2991] =   0.1622118947228123E-01;
  COFD[        2992] =  -0.2199457139050309E+02;
  COFD[        2993] =   0.5115943807127820E+01;
  COFD[        2994] =  -0.4076467799460342E+00;
  COFD[        2995] =   0.1601241241828947E-01;
  COFD[        2996] =  -0.1964411600052203E+02;
  COFD[        2997] =   0.4765927467119864E+01;
  COFD[        2998] =  -0.3953572323512016E+00;
  COFD[        2999] =   0.1677284140586059E-01;
  COFD[        3000] =  -0.1556317061681839E+02;
  COFD[        3001] =   0.3558114950801765E+01;
  COFD[        3002] =  -0.1888799785819696E+00;
  COFD[        3003] =   0.5682884908013199E-02;
  COFD[        3004] =  -0.1794643219743921E+02;
  COFD[        3005] =   0.4833075314023510E+01;
  COFD[        3006] =  -0.4248631135846847E+00;
  COFD[        3007] =   0.1892735769307228E-01;
  COFD[        3008] =  -0.1963494149489594E+02;
  COFD[        3009] =   0.4711113691258022E+01;
  COFD[        3010] =  -0.3733233453642334E+00;
  COFD[        3011] =   0.1519736886874179E-01;
  COFD[        3012] =  -0.1913891562399630E+02;
  COFD[        3013] =   0.4766338982641832E+01;
  COFD[        3014] =  -0.3968432892885747E+00;
  COFD[        3015] =   0.1688440639323740E-01;
  COFD[        3016] =  -0.1963449815268955E+02;
  COFD[        3017] =   0.4712178662790742E+01;
  COFD[        3018] =  -0.3742472943897919E+00;
  COFD[        3019] =   0.1527044453467252E-01;
  COFD[        3020] =  -0.1914452665813124E+02;
  COFD[        3021] =   0.4759724497948794E+01;
  COFD[        3022] =  -0.3961575775878822E+00;
  COFD[        3023] =   0.1686262510597603E-01;
  COFD[        3024] =  -0.1885646323568239E+02;
  COFD[        3025] =   0.4016884863488960E+01;
  COFD[        3026] =  -0.2343948477217928E+00;
  COFD[        3027] =   0.7398933738679844E-02;
  COFD[        3028] =  -0.2057011439802693E+02;
  COFD[        3029] =   0.4834011551777849E+01;
  COFD[        3030] =  -0.3728956950726914E+00;
  COFD[        3031] =   0.1455662572760375E-01;
  COFD[        3032] =  -0.1966500520249660E+02;
  COFD[        3033] =   0.4770450920673444E+01;
  COFD[        3034] =  -0.3958664514070914E+00;
  COFD[        3035] =   0.1679182031481241E-01;
  COFD[        3036] =  -0.2048311817894017E+02;
  COFD[        3037] =   0.4841596720970757E+01;
  COFD[        3038] =  -0.3768868905348668E+00;
  COFD[        3039] =   0.1483785654563742E-01;
  COFD[        3040] =  -0.2107933630762049E+02;
  COFD[        3041] =   0.5003615183707528E+01;
  COFD[        3042] =  -0.4000889881542739E+00;
  COFD[        3043] =   0.1594355520326905E-01;
  COFD[        3044] =  -0.2049312404884751E+02;
  COFD[        3045] =   0.4573043093616539E+01;
  COFD[        3046] =  -0.3164763709052368E+00;
  COFD[        3047] =   0.1128280291029378E-01;
  COFD[        3048] =  -0.2112052652732195E+02;
  COFD[        3049] =   0.5016777693101416E+01;
  COFD[        3050] =  -0.4020551643132008E+00;
  COFD[        3051] =   0.1603967677731761E-01;
  COFD[        3052] =  -0.2078560813380206E+02;
  COFD[        3053] =   0.4746137595541170E+01;
  COFD[        3054] =  -0.3464661176921034E+00;
  COFD[        3055] =   0.1285782723518177E-01;
  COFD[        3056] =  -0.1973273884329878E+02;
  COFD[        3057] =   0.4773632942733207E+01;
  COFD[        3058] =  -0.3938198339615520E+00;
  COFD[        3059] =   0.1660198213420950E-01;
  COFD[        3060] =  -0.1974560293529215E+02;
  COFD[        3061] =   0.4775430480336665E+01;
  COFD[        3062] =  -0.3942069403173756E+00;
  COFD[        3063] =   0.1662563585220353E-01;
  COFD[        3064] =  -0.1975795254936833E+02;
  COFD[        3065] =   0.4777156336087093E+01;
  COFD[        3066] =  -0.3945792342207480E+00;
  COFD[        3067] =   0.1664840211599068E-01;
  COFD[        3068] =  -0.2131043155757551E+02;
  COFD[        3069] =   0.5111234850420700E+01;
  COFD[        3070] =  -0.4165974605140776E+00;
  COFD[        3071] =   0.1676520530944528E-01;
  COFD[        3072] =  -0.2118199216189331E+02;
  COFD[        3073] =   0.4823957828991331E+01;
  COFD[        3074] =  -0.3562746590271360E+00;
  COFD[        3075] =   0.1326788764265950E-01;
  COFD[        3076] =  -0.2135235844472485E+02;
  COFD[        3077] =   0.4822077076791630E+01;
  COFD[        3078] =  -0.3537269179793018E+00;
  COFD[        3079] =   0.1307429400877126E-01;
  COFD[        3080] =  -0.2171041726819223E+02;
  COFD[        3081] =   0.5087936731469250E+01;
  COFD[        3082] =  -0.4035328710130172E+00;
  COFD[        3083] =   0.1581455637665996E-01;
  COFD[        3084] =  -0.2177389401397997E+02;
  COFD[        3085] =   0.4998448046198806E+01;
  COFD[        3086] =  -0.3836787642667724E+00;
  COFD[        3087] =   0.1463798128403765E-01;
  COFD[        3088] =  -0.2178318492664675E+02;
  COFD[        3089] =   0.5000465951311985E+01;
  COFD[        3090] =  -0.3839711955741957E+00;
  COFD[        3091] =   0.1465183998943084E-01;
  COFD[        3092] =  -0.2193017488088049E+02;
  COFD[        3093] =   0.5147264054224069E+01;
  COFD[        3094] =  -0.4122474414807874E+00;
  COFD[        3095] =   0.1623368509412096E-01;
  COFD[        3096] =  -0.2200526405431778E+02;
  COFD[        3097] =   0.5158607600017140E+01;
  COFD[        3098] =  -0.4139137284684187E+00;
  COFD[        3099] =   0.1631382826546851E-01;
  COFD[        3100] =  -0.2200526405431778E+02;
  COFD[        3101] =   0.5158607600017140E+01;
  COFD[        3102] =  -0.4139137284684187E+00;
  COFD[        3103] =   0.1631382826546851E-01;
  COFD[        3104] =  -0.2201134713900483E+02;
  COFD[        3105] =   0.5158384050126705E+01;
  COFD[        3106] =  -0.4138808904484991E+00;
  COFD[        3107] =   0.1631224885111492E-01;
  COFD[        3108] =  -0.2202438945187234E+02;
  COFD[        3109] =   0.5145495392374555E+01;
  COFD[        3110] =  -0.4119876385887463E+00;
  COFD[        3111] =   0.1622118947228123E-01;
  COFD[        3112] =  -0.2199457139050309E+02;
  COFD[        3113] =   0.5115943807127820E+01;
  COFD[        3114] =  -0.4076467799460342E+00;
  COFD[        3115] =   0.1601241241828947E-01;
  COFD[        3116] =  -0.1964411600052203E+02;
  COFD[        3117] =   0.4765927467119864E+01;
  COFD[        3118] =  -0.3953572323512016E+00;
  COFD[        3119] =   0.1677284140586059E-01;
  COFD[        3120] =  -0.1555662334580070E+02;
  COFD[        3121] =   0.3554818568932333E+01;
  COFD[        3122] =  -0.1883544064628437E+00;
  COFD[        3123] =   0.5655838508786585E-02;
  COFD[        3124] =  -0.1795336719760172E+02;
  COFD[        3125] =   0.4835810632619655E+01;
  COFD[        3126] =  -0.4252350727772528E+00;
  COFD[        3127] =   0.1894418782785674E-01;
  COFD[        3128] =  -0.1962710518342632E+02;
  COFD[        3129] =   0.4706196970979815E+01;
  COFD[        3130] =  -0.3725393175856263E+00;
  COFD[        3131] =   0.1515701000479143E-01;
  COFD[        3132] =  -0.1914947266287874E+02;
  COFD[        3133] =   0.4769147094975692E+01;
  COFD[        3134] =  -0.3971341787719369E+00;
  COFD[        3135] =   0.1689363434236904E-01;
  COFD[        3136] =  -0.1962672344426228E+02;
  COFD[        3137] =   0.4707211623266654E+01;
  COFD[        3138] =  -0.3734526020103864E+00;
  COFD[        3139] =   0.1522944338733201E-01;
  COFD[        3140] =  -0.1915545032909366E+02;
  COFD[        3141] =   0.4762618776775587E+01;
  COFD[        3142] =  -0.3964577122504647E+00;
  COFD[        3143] =   0.1687216358936520E-01;
  COFD[        3144] =  -0.1883841467755942E+02;
  COFD[        3145] =   0.4008070394200236E+01;
  COFD[        3146] =  -0.2331651697525292E+00;
  COFD[        3147] =   0.7343476142094791E-02;
  COFD[        3148] =  -0.2054662946443170E+02;
  COFD[        3149] =   0.4821857326037926E+01;
  COFD[        3150] =  -0.3710883189732943E+00;
  COFD[        3151] =   0.1446864565064257E-01;
  COFD[        3152] =  -0.1967071253473714E+02;
  COFD[        3153] =   0.4770649303948582E+01;
  COFD[        3154] =  -0.3957907059871365E+00;
  COFD[        3155] =   0.1678395401774585E-01;
  COFD[        3156] =  -0.2045813151569588E+02;
  COFD[        3157] =   0.4828580179065280E+01;
  COFD[        3158] =  -0.3749369856707143E+00;
  COFD[        3159] =   0.1474225738632903E-01;
  COFD[        3160] =  -0.2106042187580892E+02;
  COFD[        3161] =   0.4993353000413181E+01;
  COFD[        3162] =  -0.3985560772762201E+00;
  COFD[        3163] =   0.1586861625332952E-01;
  COFD[        3164] =  -0.2047732010042510E+02;
  COFD[        3165] =   0.4564542740763746E+01;
  COFD[        3166] =  -0.3152676686066944E+00;
  COFD[        3167] =   0.1122689874475350E-01;
  COFD[        3168] =  -0.2110177556093154E+02;
  COFD[        3169] =   0.5006548246650174E+01;
  COFD[        3170] =  -0.4005271173561529E+00;
  COFD[        3171] =   0.1596497410669515E-01;
  COFD[        3172] =  -0.2076598860041296E+02;
  COFD[        3173] =   0.4735797568586775E+01;
  COFD[        3174] =  -0.3449709850285935E+00;
  COFD[        3175] =   0.1278721123756670E-01;
  COFD[        3176] =  -0.1973489604720685E+02;
  COFD[        3177] =   0.4772084197088009E+01;
  COFD[        3178] =  -0.3934868007358818E+00;
  COFD[        3179] =   0.1658164649850388E-01;
  COFD[        3180] =  -0.1974787011577460E+02;
  COFD[        3181] =   0.4773895101413700E+01;
  COFD[        3182] =  -0.3938762513976486E+00;
  COFD[        3183] =   0.1660542834670177E-01;
  COFD[        3184] =  -0.1976033064762442E+02;
  COFD[        3185] =   0.4775636407093550E+01;
  COFD[        3186] =  -0.3942513288166545E+00;
  COFD[        3187] =   0.1662834933434690E-01;
  COFD[        3188] =  -0.2129337978312540E+02;
  COFD[        3189] =   0.5101265807720457E+01;
  COFD[        3190] =  -0.4151069111142752E+00;
  COFD[        3191] =   0.1669228251038974E-01;
  COFD[        3192] =  -0.2116807005556548E+02;
  COFD[        3193] =   0.4815698824842710E+01;
  COFD[        3194] =  -0.3550830672893983E+00;
  COFD[        3195] =   0.1321171158250588E-01;
  COFD[        3196] =  -0.2134436611864803E+02;
  COFD[        3197] =   0.4816403941921898E+01;
  COFD[        3198] =  -0.3529153872118850E+00;
  COFD[        3199] =   0.1303646422215331E-01;
  COFD[        3200] =  -0.2169907962284337E+02;
  COFD[        3201] =   0.5080635121967584E+01;
  COFD[        3202] =  -0.4024603628256680E+00;
  COFD[        3203] =   0.1576297563539379E-01;
  COFD[        3204] =  -0.2177183976325961E+02;
  COFD[        3205] =   0.4994984064373654E+01;
  COFD[        3206] =  -0.3831767216659396E+00;
  COFD[        3207] =   0.1461418638286886E-01;
  COFD[        3208] =  -0.2178172115101771E+02;
  COFD[        3209] =   0.4997238443372161E+01;
  COFD[        3210] =  -0.3835034599913829E+00;
  COFD[        3211] =   0.1462967284123923E-01;
  COFD[        3212] =  -0.2192865333155479E+02;
  COFD[        3213] =   0.5143921941237063E+01;
  COFD[        3214] =  -0.4117565109010827E+00;
  COFD[        3215] =   0.1621007304720648E-01;
  COFD[        3216] =  -0.2201134713900483E+02;
  COFD[        3217] =   0.5158384050126705E+01;
  COFD[        3218] =  -0.4138808904484991E+00;
  COFD[        3219] =   0.1631224885111492E-01;
  COFD[        3220] =  -0.2201134713900483E+02;
  COFD[        3221] =   0.5158384050126705E+01;
  COFD[        3222] =  -0.4138808904484991E+00;
  COFD[        3223] =   0.1631224885111492E-01;
  COFD[        3224] =  -0.2201851845855177E+02;
  COFD[        3225] =   0.5158607600017143E+01;
  COFD[        3226] =  -0.4139137284684192E+00;
  COFD[        3227] =   0.1631382826546854E-01;
  COFD[        3228] =  -0.2203870123235601E+02;
  COFD[        3229] =   0.5148647764622056E+01;
  COFD[        3230] =  -0.4124506981529867E+00;
  COFD[        3231] =   0.1624346105339284E-01;
  COFD[        3232] =  -0.2201498076448029E+02;
  COFD[        3233] =   0.5121575031494952E+01;
  COFD[        3234] =  -0.4084739492202148E+00;
  COFD[        3235] =   0.1605219530909675E-01;
  COFD[        3236] =  -0.1964994017551725E+02;
  COFD[        3237] =   0.4766172408467972E+01;
  COFD[        3238] =  -0.3952875571309086E+00;
  COFD[        3239] =   0.1676523887567621E-01;
  COFD[        3240] =  -0.1551443723602216E+02;
  COFD[        3241] =   0.3533613016857426E+01;
  COFD[        3242] =  -0.1849737082823091E+00;
  COFD[        3243] =   0.5481875490852130E-02;
  COFD[        3244] =  -0.1799730159578534E+02;
  COFD[        3245] =   0.4853169597207414E+01;
  COFD[        3246] =  -0.4275947086398001E+00;
  COFD[        3247] =   0.1905091769934767E-01;
  COFD[        3248] =  -0.1957709787054268E+02;
  COFD[        3249] =   0.4675065521077943E+01;
  COFD[        3250] =  -0.3675756867404701E+00;
  COFD[        3251] =   0.1490152940712903E-01;
  COFD[        3252] =  -0.1921620409589607E+02;
  COFD[        3253] =   0.4787003103790272E+01;
  COFD[        3254] =  -0.3989806507176665E+00;
  COFD[        3255] =   0.1695203644218339E-01;
  COFD[        3256] =  -0.1957689044936644E+02;
  COFD[        3257] =   0.4675660284789568E+01;
  COFD[        3258] =  -0.3684053408172555E+00;
  COFD[        3259] =   0.1496906925902828E-01;
  COFD[        3260] =  -0.1922460653113258E+02;
  COFD[        3261] =   0.4781057331709897E+01;
  COFD[        3262] =  -0.3983664361805447E+00;
  COFD[        3263] =   0.1693264372808513E-01;
  COFD[        3264] =  -0.1872172887795907E+02;
  COFD[        3265] =   0.3951398222423000E+01;
  COFD[        3266] =  -0.2252586833510302E+00;
  COFD[        3267] =   0.6986886080269861E-02;
  COFD[        3268] =  -0.2039108134211626E+02;
  COFD[        3269] =   0.4741987400654799E+01;
  COFD[        3270] =  -0.3592117190236912E+00;
  COFD[        3271] =   0.1389053577261777E-01;
  COFD[        3272] =  -0.1970677705449042E+02;
  COFD[        3273] =   0.4771839592787130E+01;
  COFD[        3274] =  -0.3952722585184899E+00;
  COFD[        3275] =   0.1673110821116964E-01;
  COFD[        3276] =  -0.2029132185536470E+02;
  COFD[        3277] =   0.4742426481711147E+01;
  COFD[        3278] =  -0.3620314088800772E+00;
  COFD[        3279] =   0.1410955509763312E-01;
  COFD[        3280] =  -0.2093223365539976E+02;
  COFD[        3281] =   0.4924726826746843E+01;
  COFD[        3282] =  -0.3883055489377778E+00;
  COFD[        3283] =   0.1536753032747617E-01;
  COFD[        3284] =  -0.2036999878611920E+02;
  COFD[        3285] =   0.4507734427051779E+01;
  COFD[        3286] =  -0.3071897892401887E+00;
  COFD[        3287] =   0.1085329061385159E-01;
  COFD[        3288] =  -0.2097388276834499E+02;
  COFD[        3289] =   0.4937790584419703E+01;
  COFD[        3290] =  -0.3902567887483665E+00;
  COFD[        3291] =   0.1546291085732362E-01;
  COFD[        3292] =  -0.2063237685170529E+02;
  COFD[        3293] =   0.4666357783222749E+01;
  COFD[        3294] =  -0.3349302789942741E+00;
  COFD[        3295] =   0.1231299384240708E-01;
  COFD[        3296] =  -0.1974717304559367E+02;
  COFD[        3297] =   0.4761539000992064E+01;
  COFD[        3298] =  -0.3912295069504163E+00;
  COFD[        3299] =   0.1644410435383097E-01;
  COFD[        3300] =  -0.1976079148787423E+02;
  COFD[        3301] =   0.4763387289184921E+01;
  COFD[        3302] =  -0.3916239954434353E+00;
  COFD[        3303] =   0.1646810838530097E-01;
  COFD[        3304] =  -0.1977390351786128E+02;
  COFD[        3305] =   0.4765179530799094E+01;
  COFD[        3306] =  -0.3920069566617906E+00;
  COFD[        3307] =   0.1649142352328586E-01;
  COFD[        3308] =  -0.2116627201427422E+02;
  COFD[        3309] =   0.5029554692101149E+01;
  COFD[        3310] =  -0.4043854700274245E+00;
  COFD[        3311] =   0.1616778720354987E-01;
  COFD[        3312] =  -0.2106221702527866E+02;
  COFD[        3313] =   0.4755729339634630E+01;
  COFD[        3314] =  -0.3464303724989257E+00;
  COFD[        3315] =   0.1280377425307382E-01;
  COFD[        3316] =  -0.2127828879489648E+02;
  COFD[        3317] =   0.4773820066846246E+01;
  COFD[        3318] =  -0.3468230916877841E+00;
  COFD[        3319] =   0.1275243316292433E-01;
  COFD[        3320] =  -0.2160710054610370E+02;
  COFD[        3321] =   0.5025481182776788E+01;
  COFD[        3322] =  -0.3943591536682922E+00;
  COFD[        3323] =   0.1537337185442625E-01;
  COFD[        3324] =  -0.2173580178381468E+02;
  COFD[        3325] =   0.4962896685693639E+01;
  COFD[        3326] =  -0.3785249617010963E+00;
  COFD[        3327] =   0.1439364791769246E-01;
  COFD[        3328] =  -0.2174924769698764E+02;
  COFD[        3329] =   0.4966564768133120E+01;
  COFD[        3330] =  -0.3790567913063227E+00;
  COFD[        3331] =   0.1441886477883868E-01;
  COFD[        3332] =  -0.2189359568671224E+02;
  COFD[        3333] =   0.5111398712052350E+01;
  COFD[        3334] =  -0.4069791541191179E+00;
  COFD[        3335] =   0.1598030296830180E-01;
  COFD[        3336] =  -0.2202438945187234E+02;
  COFD[        3337] =   0.5145495392374555E+01;
  COFD[        3338] =  -0.4119876385887463E+00;
  COFD[        3339] =   0.1622118947228123E-01;
  COFD[        3340] =  -0.2202438945187234E+02;
  COFD[        3341] =   0.5145495392374555E+01;
  COFD[        3342] =  -0.4119876385887463E+00;
  COFD[        3343] =   0.1622118947228123E-01;
  COFD[        3344] =  -0.2203870123235601E+02;
  COFD[        3345] =   0.5148647764622056E+01;
  COFD[        3346] =  -0.4124506981529867E+00;
  COFD[        3347] =   0.1624346105339284E-01;
  COFD[        3348] =  -0.2210739964022042E+02;
  COFD[        3349] =   0.5158607600017130E+01;
  COFD[        3350] =  -0.4139137284684174E+00;
  COFD[        3351] =   0.1631382826546844E-01;
  COFD[        3352] =  -0.2212770151166406E+02;
  COFD[        3353] =   0.5149540198412157E+01;
  COFD[        3354] =  -0.4125817900731084E+00;
  COFD[        3355] =   0.1624976613901457E-01;
  COFD[        3356] =  -0.1968677474693097E+02;
  COFD[        3357] =   0.4767670071867062E+01;
  COFD[        3358] =  -0.3948092641573855E+00;
  COFD[        3359] =   0.1671414151551012E-01;
  COFD[        3360] =  -0.1547687771509115E+02;
  COFD[        3361] =   0.3514781617610363E+01;
  COFD[        3362] =  -0.1819719591701729E+00;
  COFD[        3363] =   0.5327428541928384E-02;
  COFD[        3364] =  -0.1803535719272616E+02;
  COFD[        3365] =   0.4868249132378355E+01;
  COFD[        3366] =  -0.4296431131249400E+00;
  COFD[        3367] =   0.1914351307721480E-01;
  COFD[        3368] =  -0.1953370018155996E+02;
  COFD[        3369] =   0.4648360822131973E+01;
  COFD[        3370] =  -0.3633186948567744E+00;
  COFD[        3371] =   0.1468245773755453E-01;
  COFD[        3372] =  -0.1927347440067483E+02;
  COFD[        3373] =   0.4802472345441742E+01;
  COFD[        3374] =  -0.4005754779946568E+00;
  COFD[        3375] =   0.1700221644405138E-01;
  COFD[        3376] =  -0.1953337113047812E+02;
  COFD[        3377] =   0.4648468363013314E+01;
  COFD[        3378] =  -0.3640564212695835E+00;
  COFD[        3379] =   0.1474476403759242E-01;
  COFD[        3380] =  -0.1928407792404587E+02;
  COFD[        3381] =   0.4797069911008538E+01;
  COFD[        3382] =  -0.4000190517915601E+00;
  COFD[        3383] =   0.1698473884996324E-01;
  COFD[        3384] =  -0.1861792788931437E+02;
  COFD[        3385] =   0.3901413915830852E+01;
  COFD[        3386] =  -0.2182846281641017E+00;
  COFD[        3387] =   0.6672322090977868E-02;
  COFD[        3388] =  -0.2024909562164469E+02;
  COFD[        3389] =   0.4669841294589520E+01;
  COFD[        3390] =  -0.3484841342800297E+00;
  COFD[        3391] =   0.1336838998616518E-01;
  COFD[        3392] =  -0.1973759080924860E+02;
  COFD[        3393] =   0.4772765644230508E+01;
  COFD[        3394] =  -0.3947767037628411E+00;
  COFD[        3395] =   0.1668185546520601E-01;
  COFD[        3396] =  -0.2013725011091987E+02;
  COFD[        3397] =   0.4663752728745078E+01;
  COFD[        3398] =  -0.3502469131795646E+00;
  COFD[        3399] =   0.1353185410858909E-01;
  COFD[        3400] =  -0.2081181618605448E+02;
  COFD[        3401] =   0.4861337517486841E+01;
  COFD[        3402] =  -0.3788379332137031E+00;
  COFD[        3403] =   0.1490475778206881E-01;
  COFD[        3404] =  -0.2026866351460656E+02;
  COFD[        3405] =   0.4455185153044745E+01;
  COFD[        3406] =  -0.2997175567738866E+00;
  COFD[        3407] =   0.1050770872018541E-01;
  COFD[        3408] =  -0.2085271434976424E+02;
  COFD[        3409] =   0.4873828197990187E+01;
  COFD[        3410] =  -0.3807034496066365E+00;
  COFD[        3411] =   0.1499594029854510E-01;
  COFD[        3412] =  -0.2050575086267360E+02;
  COFD[        3413] =   0.4601709555855408E+01;
  COFD[        3414] =  -0.3255825756010016E+00;
  COFD[        3415] =   0.1187152855312509E-01;
  COFD[        3416] =  -0.1975589425049316E+02;
  COFD[        3417] =   0.4751520609448664E+01;
  COFD[        3418] =  -0.3890982795753707E+00;
  COFD[        3419] =   0.1631462582834976E-01;
  COFD[        3420] =  -0.1976999225259944E+02;
  COFD[        3421] =   0.4753334674389479E+01;
  COFD[        3422] =  -0.3894833752885402E+00;
  COFD[        3423] =   0.1633799793889408E-01;
  COFD[        3424] =  -0.1978359074953881E+02;
  COFD[        3425] =   0.4755104128079081E+01;
  COFD[        3426] =  -0.3898593304886394E+00;
  COFD[        3427] =   0.1636082491100651E-01;
  COFD[        3428] =  -0.2103101981945703E+02;
  COFD[        3429] =   0.4956251882737971E+01;
  COFD[        3430] =  -0.3934269877726836E+00;
  COFD[        3431] =   0.1563174832696796E-01;
  COFD[        3432] =  -0.2094723649561429E+02;
  COFD[        3433] =   0.4693781529535505E+01;
  COFD[        3434] =  -0.3374918691050086E+00;
  COFD[        3435] =   0.1238235308953681E-01;
  COFD[        3436] =  -0.2120031254217025E+02;
  COFD[        3437] =   0.4728097335132442E+01;
  COFD[        3438] =  -0.3402810017723648E+00;
  COFD[        3439] =   0.1244740133248284E-01;
  COFD[        3440] =  -0.2150038725326603E+02;
  COFD[        3441] =   0.4965856294754335E+01;
  COFD[        3442] =  -0.3856015603443829E+00;
  COFD[        3443] =   0.1495222567994429E-01;
  COFD[        3444] =  -0.2167081740518373E+02;
  COFD[        3445] =   0.4920315656610881E+01;
  COFD[        3446] =  -0.3723507532504461E+00;
  COFD[        3447] =   0.1410087595979374E-01;
  COFD[        3448] =  -0.2168687631153109E+02;
  COFD[        3449] =   0.4924990034391955E+01;
  COFD[        3450] =  -0.3730285599341068E+00;
  COFD[        3451] =   0.1413301763976827E-01;
  COFD[        3452] =  -0.2182529370017760E+02;
  COFD[        3453] =   0.5066380027973548E+01;
  COFD[        3454] =  -0.4003664957816425E+00;
  COFD[        3455] =   0.1566227519627694E-01;
  COFD[        3456] =  -0.2199457139050309E+02;
  COFD[        3457] =   0.5115943807127820E+01;
  COFD[        3458] =  -0.4076467799460342E+00;
  COFD[        3459] =   0.1601241241828947E-01;
  COFD[        3460] =  -0.2199457139050309E+02;
  COFD[        3461] =   0.5115943807127820E+01;
  COFD[        3462] =  -0.4076467799460342E+00;
  COFD[        3463] =   0.1601241241828947E-01;
  COFD[        3464] =  -0.2201498076448029E+02;
  COFD[        3465] =   0.5121575031494952E+01;
  COFD[        3466] =  -0.4084739492202148E+00;
  COFD[        3467] =   0.1605219530909675E-01;
  COFD[        3468] =  -0.2212770151166406E+02;
  COFD[        3469] =   0.5149540198412157E+01;
  COFD[        3470] =  -0.4125817900731084E+00;
  COFD[        3471] =   0.1624976613901457E-01;
  COFD[        3472] =  -0.2219216921513908E+02;
  COFD[        3473] =   0.5158607600017143E+01;
  COFD[        3474] =  -0.4139137284684192E+00;
  COFD[        3475] =   0.1631382826546853E-01;
  COFD[        3476] =  -0.1971828836771711E+02;
  COFD[        3477] =   0.4768875742807525E+01;
  COFD[        3478] =  -0.3943502742168009E+00;
  COFD[        3479] =   0.1666648262004003E-01;
  COFD[        3480] =  -0.1711567017943538E+02;
  COFD[        3481] =   0.4833264586465570E+01;
  COFD[        3482] =  -0.4205908440542898E+00;
  COFD[        3483] =   0.1855345683322859E-01;
  COFD[        3484] =  -0.1266066516595179E+02;
  COFD[        3485] =   0.2898076057146818E+01;
  COFD[        3486] =  -0.1700453759318706E+00;
  COFD[        3487] =   0.7738690296212197E-02;
  COFD[        3488] =  -0.1658229952875101E+02;
  COFD[        3489] =   0.3922184822536529E+01;
  COFD[        3490] =  -0.2983959485115160E+00;
  COFD[        3491] =   0.1309927190981370E-01;
  COFD[        3492] =  -0.1474819919893963E+02;
  COFD[        3493] =   0.3361311502126538E+01;
  COFD[        3494] =  -0.2285465363406641E+00;
  COFD[        3495] =   0.1019990809984005E-01;
  COFD[        3496] =  -0.1651845590523054E+02;
  COFD[        3497] =   0.3896132878391479E+01;
  COFD[        3498] =  -0.2951170694006329E+00;
  COFD[        3499] =   0.1296170338310411E-01;
  COFD[        3500] =  -0.1473661449825940E+02;
  COFD[        3501] =   0.3348204558833826E+01;
  COFD[        3502] =  -0.2267271723233180E+00;
  COFD[        3503] =   0.1011600240359858E-01;
  COFD[        3504] =  -0.2024331300258610E+02;
  COFD[        3505] =   0.5167058816761965E+01;
  COFD[        3506] =  -0.4292618020154789E+00;
  COFD[        3507] =   0.1752039422653532E-01;
  COFD[        3508] =  -0.1848414582884554E+02;
  COFD[        3509] =   0.4538825155691820E+01;
  COFD[        3510] =  -0.3723671333664587E+00;
  COFD[        3511] =   0.1605610876285614E-01;
  COFD[        3512] =  -0.1549017806948574E+02;
  COFD[        3513] =   0.3466859355125401E+01;
  COFD[        3514] =  -0.2408856054712291E+00;
  COFD[        3515] =   0.1067561900190753E-01;
  COFD[        3516] =  -0.1816875575125756E+02;
  COFD[        3517] =   0.4450099332111185E+01;
  COFD[        3518] =  -0.3625779431233763E+00;
  COFD[        3519] =   0.1570387818274669E-01;
  COFD[        3520] =  -0.1856460854973782E+02;
  COFD[        3521] =   0.4508650490412063E+01;
  COFD[        3522] =  -0.3698234379228469E+00;
  COFD[        3523] =   0.1600311768013340E-01;
  COFD[        3524] =  -0.2007323696432831E+02;
  COFD[        3525] =   0.5032033200707714E+01;
  COFD[        3526] =  -0.4235009227577490E+00;
  COFD[        3527] =   0.1774762291235686E-01;
  COFD[        3528] =  -0.1857749863267182E+02;
  COFD[        3529] =   0.4510414228568457E+01;
  COFD[        3530] =  -0.3700425246017726E+00;
  COFD[        3531] =   0.1601211519431907E-01;
  COFD[        3532] =  -0.2000088490858467E+02;
  COFD[        3533] =   0.5021789747731799E+01;
  COFD[        3534] =  -0.4253480595637520E+00;
  COFD[        3535] =   0.1795667045275812E-01;
  COFD[        3536] =  -0.1578873902175925E+02;
  COFD[        3537] =   0.3589258101555060E+01;
  COFD[        3538] =  -0.2568718883472193E+00;
  COFD[        3539] =   0.1137271971440306E-01;
  COFD[        3540] =  -0.1579929195048582E+02;
  COFD[        3541] =   0.3590679879226839E+01;
  COFD[        3542] =  -0.2570681340587627E+00;
  COFD[        3543] =   0.1138172060899465E-01;
  COFD[        3544] =  -0.1581023736073806E+02;
  COFD[        3545] =   0.3592405682636056E+01;
  COFD[        3546] =  -0.2573062871053550E+00;
  COFD[        3547] =   0.1139264138769127E-01;
  COFD[        3548] =  -0.1850244542172855E+02;
  COFD[        3549] =   0.4509546506902997E+01;
  COFD[        3550] =  -0.3699250247632380E+00;
  COFD[        3551] =   0.1600569104209367E-01;
  COFD[        3552] =  -0.1994083366515423E+02;
  COFD[        3553] =   0.4933053831289947E+01;
  COFD[        3554] =  -0.4129818889961173E+00;
  COFD[        3555] =   0.1738193532007086E-01;
  COFD[        3556] =  -0.2013110520151531E+02;
  COFD[        3557] =   0.4933748223113413E+01;
  COFD[        3558] =  -0.4099898956624159E+00;
  COFD[        3559] =   0.1712216655141095E-01;
  COFD[        3560] =  -0.1952277059880253E+02;
  COFD[        3561] =   0.4760388219539090E+01;
  COFD[        3562] =  -0.3964022009449412E+00;
  COFD[        3563] =   0.1689791166807166E-01;
  COFD[        3564] =  -0.1995821009052950E+02;
  COFD[        3565] =   0.4847729338512834E+01;
  COFD[        3566] =  -0.4012454442759492E+00;
  COFD[        3567] =   0.1683646904857131E-01;
  COFD[        3568] =  -0.1995722688303096E+02;
  COFD[        3569] =   0.4845641046786140E+01;
  COFD[        3570] =  -0.4008772701599357E+00;
  COFD[        3571] =   0.1681620052332664E-01;
  COFD[        3572] =  -0.1959997385926940E+02;
  COFD[        3573] =   0.4764004414585760E+01;
  COFD[        3574] =  -0.3958298367826015E+00;
  COFD[        3575] =   0.1682596038127772E-01;
  COFD[        3576] =  -0.1964411600052203E+02;
  COFD[        3577] =   0.4765927467119864E+01;
  COFD[        3578] =  -0.3953572323512016E+00;
  COFD[        3579] =   0.1677284140586059E-01;
  COFD[        3580] =  -0.1964411600052203E+02;
  COFD[        3581] =   0.4765927467119864E+01;
  COFD[        3582] =  -0.3953572323512016E+00;
  COFD[        3583] =   0.1677284140586059E-01;
  COFD[        3584] =  -0.1964994017551725E+02;
  COFD[        3585] =   0.4766172408467972E+01;
  COFD[        3586] =  -0.3952875571309086E+00;
  COFD[        3587] =   0.1676523887567621E-01;
  COFD[        3588] =  -0.1968677474693097E+02;
  COFD[        3589] =   0.4767670071867062E+01;
  COFD[        3590] =  -0.3948092641573855E+00;
  COFD[        3591] =   0.1671414151551012E-01;
  COFD[        3592] =  -0.1971828836771711E+02;
  COFD[        3593] =   0.4768875742807525E+01;
  COFD[        3594] =  -0.3943502742168009E+00;
  COFD[        3595] =   0.1666648262004003E-01;
  COFD[        3596] =  -0.1546612761154138E+02;
  COFD[        3597] =   0.3460858520880728E+01;
  COFD[        3598] =  -0.2401164641793465E+00;
  COFD[        3599] =   0.1064270570979806E-01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetKTDIF EGTRANSETKTDIF
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetKTDIF egtransetktdif
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetKTDIF egtransetktdif_
#endif
void egtransetKTDIF(int* KTDIF) {
  KTDIF[           0] =            1;
  KTDIF[           1] =            2;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFTD EGTRANSETCOFTD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFTD egtransetcoftd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFTD egtransetcoftd_
#endif
void egtransetCOFTD(double* COFTD) {
  COFTD[           0] =   0.0000000000000000E+00;
  COFTD[           1] =   0.0000000000000000E+00;
  COFTD[           2] =   0.0000000000000000E+00;
  COFTD[           3] =   0.0000000000000000E+00;
  COFTD[           4] =   0.1441522353365304E+00;
  COFTD[           5] =   0.7999936505724083E-04;
  COFTD[           6] =  -0.4897074975189896E-07;
  COFTD[           7] =   0.9142773649236618E-11;
  COFTD[           8] =   0.1000392102974500E+00;
  COFTD[           9] =   0.6504687084211580E-03;
  COFTD[          10] =  -0.3417790254177462E-06;
  COFTD[          11] =   0.5627791733432544E-10;
  COFTD[          12] =   0.2352832283423578E+00;
  COFTD[          13] =   0.4656706350120388E-03;
  COFTD[          14] =  -0.2609398458577575E-06;
  COFTD[          15] =   0.4492718580048779E-10;
  COFTD[          16] =   0.1051242233178346E+00;
  COFTD[          17] =   0.6506660047178345E-03;
  COFTD[          18] =  -0.3425645641009918E-06;
  COFTD[          19] =   0.5648041610405423E-10;
  COFTD[          20] =   0.2370534621354604E+00;
  COFTD[          21] =   0.4691742672103970E-03;
  COFTD[          22] =  -0.2629031159826994E-06;
  COFTD[          23] =   0.4526521083989819E-10;
  COFTD[          24] =  -0.1743526120296758E+00;
  COFTD[          25] =   0.8622469282596517E-03;
  COFTD[          26] =  -0.3795455129739271E-06;
  COFTD[          27] =   0.5602621271895701E-10;
  COFTD[          28] =  -0.3814697966940322E-01;
  COFTD[          29] =   0.8398335469283959E-03;
  COFTD[          30] =  -0.4116889425506645E-06;
  COFTD[          31] =   0.6491249948882154E-10;
  COFTD[          32] =   0.2001200095273127E+00;
  COFTD[          33] =   0.5647937464189248E-03;
  COFTD[          34] =  -0.3094455087464756E-06;
  COFTD[          35] =   0.5241393748626840E-10;
  COFTD[          36] =  -0.1420994083273805E-01;
  COFTD[          37] =   0.8238121586424461E-03;
  COFTD[          38] =  -0.4089955430834602E-06;
  COFTD[          39] =   0.6498993528275147E-10;
  COFTD[          40] =  -0.2281049585009957E-01;
  COFTD[          41] =   0.8334704591750301E-03;
  COFTD[          42] =  -0.4119691401467454E-06;
  COFTD[          43] =   0.6528594139428977E-10;
  COFTD[          44] =  -0.1613574734420234E+00;
  COFTD[          45] =   0.9059203176635061E-03;
  COFTD[          46] =  -0.4078791780205072E-06;
  COFTD[          47] =   0.6106263270843215E-10;
  COFTD[          48] =  -0.2286365866852431E-01;
  COFTD[          49] =   0.8354129701565902E-03;
  COFTD[          50] =  -0.4129292876480630E-06;
  COFTD[          51] =   0.6543809874636368E-10;
  COFTD[          52] =  -0.1312444259932660E+00;
  COFTD[          53] =   0.9039014427379127E-03;
  COFTD[          54] =  -0.4178315332949202E-06;
  COFTD[          55] =   0.6357257058426966E-10;
  COFTD[          56] =   0.1798404103253150E+00;
  COFTD[          57] =   0.6017229463172420E-03;
  COFTD[          58] =  -0.3264339193657546E-06;
  COFTD[          59] =   0.5491123423829806E-10;
  COFTD[          60] =   0.1801870770056730E+00;
  COFTD[          61] =   0.6028828485656729E-03;
  COFTD[          62] =  -0.3270631648338786E-06;
  COFTD[          63] =   0.5501708305866783E-10;
  COFTD[          64] =   0.1805137892277234E+00;
  COFTD[          65] =   0.6039759857559981E-03;
  COFTD[          66] =  -0.3276561903444169E-06;
  COFTD[          67] =   0.5511683912195231E-10;
  COFTD[          68] =  -0.2003084371216687E-01;
  COFTD[          69] =   0.8504401712698103E-03;
  COFTD[          70] =  -0.4210644955951721E-06;
  COFTD[          71] =   0.6679597524516173E-10;
  COFTD[          72] =  -0.1419517543438646E+00;
  COFTD[          73] =   0.9234297382553384E-03;
  COFTD[          74] =  -0.4241404020095680E-06;
  COFTD[          75] =   0.6428102343789282E-10;
  COFTD[          76] =  -0.1556824608103752E+00;
  COFTD[          77] =   0.9271592371665584E-03;
  COFTD[          78] =  -0.4210744576330030E-06;
  COFTD[          79] =   0.6337596666091549E-10;
  COFTD[          80] =  -0.8712265411232842E-01;
  COFTD[          81] =   0.9023015642605635E-03;
  COFTD[          82] =  -0.4305715034288798E-06;
  COFTD[          83] =   0.6678968571533750E-10;
  COFTD[          84] =  -0.1305225264594721E+00;
  COFTD[          85] =   0.9307619077421273E-03;
  COFTD[          86] =  -0.4317978939742757E-06;
  COFTD[          87] =   0.6584179186431713E-10;
  COFTD[          88] =  -0.1305973678476730E+00;
  COFTD[          89] =   0.9312956049907893E-03;
  COFTD[          90] =  -0.4320454861308468E-06;
  COFTD[          91] =   0.6587954543252040E-10;
  COFTD[          92] =  -0.8811128910364355E-01;
  COFTD[          93] =   0.9125405418059032E-03;
  COFTD[          94] =  -0.4354574663152285E-06;
  COFTD[          95] =   0.6754758985668746E-10;
  COFTD[          96] =  -0.8860760312277161E-01;
  COFTD[          97] =   0.9176807079359012E-03;
  COFTD[          98] =  -0.4379103148374198E-06;
  COFTD[          99] =   0.6792807249569135E-10;
  COFTD[         100] =  -0.8860760312277161E-01;
  COFTD[         101] =   0.9176807079359012E-03;
  COFTD[         102] =  -0.4379103148374198E-06;
  COFTD[         103] =   0.6792807249569135E-10;
  COFTD[         104] =  -0.8866990566937646E-01;
  COFTD[         105] =   0.9183259555564174E-03;
  COFTD[         106] =  -0.4382182221369902E-06;
  COFTD[         107] =   0.6797583467132073E-10;
  COFTD[         108] =  -0.8904854192934855E-01;
  COFTD[         109] =   0.9222473706365612E-03;
  COFTD[         110] =  -0.4400894918471452E-06;
  COFTD[         111] =   0.6826610356936570E-10;
  COFTD[         112] =  -0.8935321462130538E-01;
  COFTD[         113] =   0.9254027686136021E-03;
  COFTD[         114] =  -0.4415952239711987E-06;
  COFTD[         115] =   0.6849967075747767E-10;
  COFTD[         116] =   0.2015217542366260E+00;
  COFTD[         117] =   0.5627441294507066E-03;
  COFTD[         118] =  -0.3085192624064295E-06;
  COFTD[         119] =   0.5228060239663243E-10;
  COFTD[         120] =  -0.1441522353365304E+00;
  COFTD[         121] =  -0.7999936505724083E-04;
  COFTD[         122] =   0.4897074975189896E-07;
  COFTD[         123] =  -0.9142773649236618E-11;
  COFTD[         124] =   0.0000000000000000E+00;
  COFTD[         125] =   0.0000000000000000E+00;
  COFTD[         126] =   0.0000000000000000E+00;
  COFTD[         127] =   0.0000000000000000E+00;
  COFTD[         128] =   0.3311912921591189E+00;
  COFTD[         129] =   0.1813267307504234E-03;
  COFTD[         130] =  -0.1110964046786237E-06;
  COFTD[         131] =   0.2076359836954698E-10;
  COFTD[         132] =   0.4066826061623353E+00;
  COFTD[         133] =   0.3847052789609525E-04;
  COFTD[         134] =  -0.2548469497314269E-07;
  COFTD[         135] =   0.5863025071898190E-11;
  COFTD[         136] =   0.3395573510152542E+00;
  COFTD[         137] =   0.1793350524644130E-03;
  COFTD[         138] =  -0.1101357188562818E-06;
  COFTD[         139] =   0.2064272631295983E-10;
  COFTD[         140] =   0.4128957304916273E+00;
  COFTD[         141] =   0.3905826430087406E-04;
  COFTD[         142] =  -0.2587403933152670E-07;
  COFTD[         143] =   0.5952597881665513E-11;
  COFTD[         144] =   0.2274730213366172E-01;
  COFTD[         145] =   0.6730786443559050E-03;
  COFTD[         146] =  -0.3409357669584640E-06;
  COFTD[         147] =   0.5484991439904711E-10;
  COFTD[         148] =   0.2580669416192538E+00;
  COFTD[         149] =   0.4050726264183527E-03;
  COFTD[         150] =  -0.2305874642262042E-06;
  COFTD[         151] =   0.4018638759886466E-10;
  COFTD[         152] =   0.4306056706094065E+00;
  COFTD[         153] =   0.9359619868800346E-04;
  COFTD[         154] =  -0.6039837349070145E-07;
  COFTD[         155] =   0.1231151900773763E-10;
  COFTD[         156] =   0.2829745044907168E+00;
  COFTD[         157] =   0.3730329798247335E-03;
  COFTD[         158] =  -0.2149591816261514E-06;
  COFTD[         159] =   0.3783551881490937E-10;
  COFTD[         160] =   0.2767260753066081E+00;
  COFTD[         161] =   0.3877928501945458E-03;
  COFTD[         162] =  -0.2225046014869723E-06;
  COFTD[         163] =   0.3902511769192158E-10;
  COFTD[         164] =   0.1226934830045553E+00;
  COFTD[         165] =   0.6212781891511368E-03;
  COFTD[         166] =  -0.3299652330704521E-06;
  COFTD[         167] =   0.5471615878539925E-10;
  COFTD[         168] =   0.2780220079882950E+00;
  COFTD[         169] =   0.3896089184047235E-03;
  COFTD[         170] =  -0.2235466102119295E-06;
  COFTD[         171] =   0.3920787576908362E-10;
  COFTD[         172] =   0.1403144190321918E+00;
  COFTD[         173] =   0.6012660048574510E-03;
  COFTD[         174] =  -0.3219150906946780E-06;
  COFTD[         175] =   0.5366790096442973E-10;
  COFTD[         176] =   0.4265800666136179E+00;
  COFTD[         177] =   0.1204072844552778E-03;
  COFTD[         178] =  -0.7672988784371693E-07;
  COFTD[         179] =   0.1520903572120744E-10;
  COFTD[         180] =   0.4282310128428898E+00;
  COFTD[         181] =   0.1208732836141846E-03;
  COFTD[         182] =  -0.7702684714613741E-07;
  COFTD[         183] =   0.1526789759061949E-10;
  COFTD[         184] =   0.4297895875586662E+00;
  COFTD[         185] =   0.1213132098175750E-03;
  COFTD[         186] =  -0.7730719138277033E-07;
  COFTD[         187] =   0.1532346610021863E-10;
  COFTD[         188] =   0.2931916398206432E+00;
  COFTD[         189] =   0.4014300382755085E-03;
  COFTD[         190] =  -0.2307057834480938E-06;
  COFTD[         191] =   0.4051766204108762E-10;
  COFTD[         192] =   0.1599912935342585E+00;
  COFTD[         193] =   0.6054913477163982E-03;
  COFTD[         194] =  -0.3262695984490298E-06;
  COFTD[         195] =   0.5463067913108842E-10;
  COFTD[         196] =   0.1422359682361410E+00;
  COFTD[         197] =   0.6329414071426726E-03;
  COFTD[         198] =  -0.3382378176090615E-06;
  COFTD[         199] =   0.5631821549456155E-10;
  COFTD[         200] =   0.2270093817830355E+00;
  COFTD[         201] =   0.5111594095278553E-03;
  COFTD[         202] =  -0.2840231885195104E-06;
  COFTD[         203] =   0.4859641263479976E-10;
  COFTD[         204] =   0.1808705686939265E+00;
  COFTD[         205] =   0.5957384609631979E-03;
  COFTD[         206] =  -0.3234737254417002E-06;
  COFTD[         207] =   0.5444656545597373E-10;
  COFTD[         208] =   0.1810782285725225E+00;
  COFTD[         209] =   0.5964224361249499E-03;
  COFTD[         210] =  -0.3238451098799718E-06;
  COFTD[         211] =   0.5450907627381549E-10;
  COFTD[         212] =   0.2321964078855084E+00;
  COFTD[         213] =   0.5228390906886971E-03;
  COFTD[         214] =  -0.2905129453788386E-06;
  COFTD[         215] =   0.4970681106346105E-10;
  COFTD[         216] =   0.2348213595204900E+00;
  COFTD[         217] =   0.5287497218583759E-03;
  COFTD[         218] =  -0.2937971582480208E-06;
  COFTD[         219] =   0.5026874040663625E-10;
  COFTD[         220] =   0.2348213595204900E+00;
  COFTD[         221] =   0.5287497218583759E-03;
  COFTD[         222] =  -0.2937971582480208E-06;
  COFTD[         223] =   0.5026874040663625E-10;
  COFTD[         224] =   0.2351518692719966E+00;
  COFTD[         225] =   0.5294939341376058E-03;
  COFTD[         226] =  -0.2942106761067208E-06;
  COFTD[         227] =   0.5033949337789185E-10;
  COFTD[         228] =   0.2371653158537335E+00;
  COFTD[         229] =   0.5340276329554817E-03;
  COFTD[         230] =  -0.2967298033496871E-06;
  COFTD[         231] =   0.5077051644898994E-10;
  COFTD[         232] =   0.2387914779012326E+00;
  COFTD[         233] =   0.5376892791194756E-03;
  COFTD[         234] =  -0.2987643788643753E-06;
  COFTD[         235] =   0.5111863264247342E-10;
  COFTD[         236] =   0.4313313898218810E+00;
  COFTD[         237] =   0.9205368774727828E-04;
  COFTD[         238] =  -0.5945097223307069E-07;
  COFTD[         239] =   0.1214380118036443E-10;
};




#if 0




\\
\\
\\  This is the mechanism file
\\
\\

ELEMENTS
C H O N 
END

SPECIES
H	H2	CH3	O	CH4
OH	H2O	C2H2	CO	C2H4
C2H5	CH2O	C2H6	CH3O	O2
HO2	H2O2	CO2	CH3HCO	HCOOH
CH3OCH3	CH3OCO	CH3OCHO	CH3OCH2OH
OCH2OCHO	HOCH2OCO	CH3OCH2O2
HO2CH2OCHO	O2CH2OCH2O2H	N2	
END

REACTIONS
!1
H+O2 = O+OH                                 0.0	       0.0   0.0
!2
O+H2 = H+OH                                 0.0	       0.0   0.0 
!3
H2+OH = H2O+H                               0.0	       0.0   0.0
END

\\
\\
\\  This is the therm file
\\
\\
THERMO ALL
0300.00  1000.00  5000.00
! From GRI mechanism
HO2               L 5/89H   1O   2   00   00G   200.000  3500.000  1000.000    1
 4.01721090E+00 2.23982013E-03-6.33658150E-07 1.14246370E-10-1.07908535E-14    2
 1.11856713E+02 3.78510215E+00 4.30179801E+00-4.74912051E-03 2.11582891E-05    3
-2.42763894E-08 9.29225124E-12 2.94808040E+02 3.71666245E+00 1.00021620E+04    4
CH2O            NIST/98 C   1H   2O   1    0G   300.000  4000.000 1200.00      1
 0.51481905E+01 0.28678016E-02-0.23782633E-06-0.16111303E-09 0.28566735E-13    2
-0.16230173E+05-0.51213813E+01 0.26962612E+01 0.49261423E-02 0.82826494E-06    3
-0.55038196E-09-0.39610326E-12-0.14970793E+05 0.94697599E+01                   4
AR                120186AR  1               G  0300.00   5000.00  1000.00      1
 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-0.07453750E+04 0.04366001E+02 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-0.07453750E+04 0.04366001E+02                   4
C2H6              121686C   2H   6          G  0300.00   4000.00  1000.00      1
 0.04825938E+02 0.01384043E+00-0.04557259E-04 0.06724967E-08-0.03598161E-12    2
-0.01271779E+06-0.05239507E+02 0.01462539E+02 0.01549467E+00 0.05780507E-04    3
-0.01257832E-06 0.04586267E-10-0.01123918E+06 0.01443229E+03                   4
CH3O              121686C   1H   3O   1     G  0300.00   3000.00  1000.00      1
 0.03770800E+02 0.07871497E-01-0.02656384E-04 0.03944431E-08-0.02112616E-12    2
 0.01278325E+04 0.02929575E+02 0.02106204E+02 0.07216595E-01 0.05338472E-04    3
-0.07377636E-07 0.02075611E-10 0.09786011E+04 0.01315218E+03                   4
CH3OH             121686C   1H   4O   1     G  0300.00   5000.00  1000.00      1
 0.04029061E+02 0.09376593E-01-0.03050254E-04 0.04358793E-08-0.02224723E-12    2
-0.02615791E+06 0.02378196E+02 0.02660115E+02 0.07341508E-01 0.07170051E-04    3
-0.08793194E-07 0.02390570E-10-0.02535348E+06 0.01123263E+03                   4
CH4               121286C   1H   4          G  0300.00   5000.00  1000.00      1
 0.01683479E+02 0.01023724E+00-0.03875129E-04 0.06785585E-08-0.04503423E-12    2
-0.01008079E+06 0.09623395E+02 0.07787415E+01 0.01747668E+00-0.02783409E-03    3
 0.03049708E-06-0.01223931E-09-0.09825229E+05 0.01372219E+03                   4
CO                121286C   1O   1          G  0300.00   5000.00  1000.00      1
 0.03025078E+02 0.01442689E-01-0.05630828E-05 0.01018581E-08-0.06910952E-13    2
-0.01426835E+06 0.06108218E+02 0.03262452E+02 0.01511941E-01-0.03881755E-04    3
 0.05581944E-07-0.02474951E-10-0.01431054E+06 0.04848897E+02                   4
CO2               121286C   1O   2          G  0300.00   5000.00  1000.00      1
 0.04453623E+02 0.03140169E-01-0.01278411E-04 0.02393997E-08-0.01669033E-12    2
-0.04896696E+06-0.09553959E+01 0.02275725E+02 0.09922072E-01-0.01040911E-03    3
 0.06866687E-07-0.02117280E-10-0.04837314E+06 0.01018849E+03                   4
H                 120186H   1               G  0300.00   5000.00  1000.00      1
 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
 0.02547163E+06-0.04601176E+01 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 0.02547163E+06-0.04601176E+01                   4
H2                121286H   2               G  0300.00   5000.00  1000.00      1
 0.02991423E+02 0.07000644E-02-0.05633829E-06-0.09231578E-10 0.01582752E-13    2
-0.08350340E+04-0.01355110E+02 0.03298124E+02 0.08249442E-02-0.08143015E-05    3
-0.09475434E-09 0.04134872E-11-0.01012521E+05-0.03294094E+02                   4
H2O                20387H   2O   1          G  0300.00   5000.00  1000.00      1
 0.02672146E+02 0.03056293E-01-0.08730260E-05 0.01200996E-08-0.06391618E-13    2
-0.02989921E+06 0.06862817E+02 0.03386842E+02 0.03474982E-01-0.06354696E-04    3
 0.06968581E-07-0.02506588E-10-0.03020811E+06 0.02590233E+02                   4
H2O2              120186H   2O   2          G  0300.00   5000.00  1000.00      1
 0.04573167E+02 0.04336136E-01-0.01474689E-04 0.02348904E-08-0.01431654E-12    2
-0.01800696E+06 0.05011370E+01 0.03388754E+02 0.06569226E-01-0.01485013E-05    3
-0.04625806E-07 0.02471515E-10-0.01766315E+06 0.06785363E+02                   4
HCO               121286H   1C   1O   1     G  0300.00   5000.00  1000.00      1
 0.03557271E+02 0.03345573E-01-0.01335006E-04 0.02470573E-08-0.01713851E-12    2
 0.03916324E+05 0.05552299E+02 0.02898330E+02 0.06199147E-01-0.09623084E-04    3
 0.01089825E-06-0.04574885E-10 0.04159922E+05 0.08983614E+02                   4
O                 120186O   1               G  0300.00   5000.00  1000.00      1
 0.02542060E+02-0.02755062E-03-0.03102803E-07 0.04551067E-10-0.04368052E-14    2
 0.02923080E+06 0.04920308E+02 0.02946429E+02-0.01638166E-01 0.02421032E-04    3
-0.01602843E-07 0.03890696E-11 0.02914764E+06 0.02963995E+02                   4
O2                121386O   2               G  0300.00   5000.00  1000.00      1
 0.03697578E+02 0.06135197E-02-0.01258842E-05 0.01775281E-09-0.01136435E-13    2
-0.01233930E+05 0.03189166E+02 0.03212936E+02 0.01127486E-01-0.05756150E-05    3
 0.01313877E-07-0.08768554E-11-0.01005249E+05 0.06034738E+02                   4
N2                121286N   2               G  0300.00   5000.00  1000.00      1
 0.02926640E+02 0.01487977E-01-0.05684761E-05 0.01009704E-08-0.06753351E-13    2
-0.09227977E+04 0.05980528E+02 0.03298677E+02 0.01408240E-01-0.03963222E-04    3
 0.05641515E-07-0.02444855E-10-0.01020900E+05 0.03950372E+02                   4
!
! Henry Curran
HOCH2O     8/ 4/99 THERMC   1H   3O   2    0G   300.000  5000.000 1452.000     1
 6.39521515E+00 7.43673043E-03-2.50422354E-06 3.84879712E-10-2.21778689E-14    2
-2.47500385E+04-7.29290847E+00 4.11183145E+00 7.53850697E-03 3.77337370E-06    3
-5.38746005E-09 1.45615887E-12-2.34414546E+04 6.81381989E+00                   4
!
C2H3              L 2/92C   2H   3   00   00G   200.000  3500.000  1000.000    1
 3.01672400E+00 1.03302292E-02-4.68082349E-06 1.01763288E-09-8.62607041E-14    2
 3.46128739E+04 7.78732378E+00 3.21246645E+00 1.51479162E-03 2.59209412E-05    3
-3.57657847E-08 1.47150873E-11 3.48598468E+04 8.51054025E+00 1.05750490E+04    4
C2H4              L 1/91C   2H   4   00   00G   200.000  3500.000  1000.000    1
 2.03611116E+00 1.46454151E-02-6.71077915E-06 1.47222923E-09-1.25706061E-13    2
 4.93988614E+03 1.03053693E+01 3.95920148E+00-7.57052247E-03 5.70990292E-05    3
-6.91588753E-08 2.69884373E-11 5.08977593E+03 4.09733096E+00 1.05186890E+04    4
!
! A. Burcat
CH2HCO            T04/83O   1H   3C   2    0G   300.     5000.                 1
 0.59756699E+01 0.81305914E-02-0.27436245E-05 0.40703041E-09-0.21760171E-13    2
 0.49032178E+03-0.50320879E+01 0.34090624E+01 0.10738574E-01 0.18914925E-05    3
-0.71585831E-08 0.28673851E-11 0.15214766E+04 0.95714535E+01                   4
CH2CO             L 5/90C   2H   2O   1   00G   200.000  3500.000  1000.000    1
 4.51129732E+00 9.00359745E-03-4.16939635E-06 9.23345882E-10-7.94838201E-14    2
-7.55105311E+03 6.32247205E-01 2.13583630E+00 1.81188721E-02-1.73947474E-05    3
 9.34397568E-09-2.01457615E-12-7.04291804E+03 1.22156480E+01 1.17977430E+04    4
CH3HCO            L 8/88C   2H   4O   1    0G   200.000  6000.0    1000.0      1
 0.54041108E+01 0.11723059E-01-0.42263137E-05 0.68372451E-09-0.40984863E-13    2
-0.22593122E+05-0.34807917E+01 0.47294595E+01-0.31932858E-02 0.47534921E-04    3
-0.57458611E-07 0.21931112E-10-0.21572878E+05 0.41030159E+01-0.19987949E+05    4
C2H2              L 1/91C   2H   2   00   00G   200.000  3500.000  1000.000    1
 4.14756964E+00 5.96166664E-03-2.37294852E-06 4.67412171E-10-3.61235213E-14    2
 2.59359992E+04-1.23028121E+00 8.08681094E-01 2.33615629E-02-3.55171815E-05    3
 2.80152437E-08-8.50072974E-12 2.64289807E+04 1.39397051E+01 1.00058390E+04    4
HCCO              SRIC91H   1C   2O   1     G  0300.00   4000.00  1000.00      1
 0.56282058E+01 0.40853401E-02-0.15934547E-05 0.28626052E-09-0.19407832E-13    2
 0.19327215E+05-0.39302595E+01 0.22517214E+01 0.17655021E-01-0.23729101E-04    3
 0.17275759E-07-0.50664811E-11 0.20059449E+05 0.12490417E+02                   4
C2H               L 1/91C   2H   1   00   00G   200.000  3500.000  1000.000    1
 3.16780652E+00 4.75221902E-03-1.83787077E-06 3.04190252E-10-1.77232770E-14    2
 6.71210650E+04 6.63589475E+00 2.88965733E+00 1.34099611E-02-2.84769501E-05    3
 2.94791045E-08-1.09331511E-11 6.68393932E+04 6.22296438E+00 1.04544720E+04    4
CH2               L S/93C   1H   2   00   00G   200.000  3500.000  1000.000    1
 2.87410113E+00 3.65639292E-03-1.40894597E-06 2.60179549E-10-1.87727567E-14    2
 4.62636040E+04 6.17119324E+00 3.76267867E+00 9.68872143E-04 2.79489841E-06    3
-3.85091153E-09 1.68741719E-12 4.60040401E+04 1.56253185E+00 1.00274170E+04    4
CH2(S)            L S/93C   1H   2   00   00G   200.000  3500.000  1000.000    1
 2.29203842E+00 4.65588637E-03-2.01191947E-06 4.17906000E-10-3.39716365E-14    2
 5.09259997E+04 8.62650169E+00 4.19860411E+00-2.36661419E-03 8.23296220E-06    3
-6.68815981E-09 1.94314737E-12 5.04968163E+04-7.69118967E-01 9.93967200E+03    4
HCCOH              SRI91C   2O   1H   20   0G   300.000  5000.000   1000.0     1
 0.59238291E+01 0.67923600E-02-0.25658564E-05 0.44987841E-09-0.29940101E-13    2
 0.72646260E+04-0.76017742E+01 0.12423733E+01 0.31072201E-01-0.50866864E-04    3
 0.43137131E-07-0.14014594E-10 0.80316143E+04 0.13874319E+02                   4
CH3CO             T 9/92C   2H   3O   1    0G   200.000  6000.0    1000.0      1
 0.59447731E+01 0.78667205E-02-0.28865882E-05 0.47270875E-09-0.28599861E-13    2
-0.37873075E+04-0.50136751E+01 0.41634257E+01-0.23261610E-03 0.34267820E-04    3
-0.44105227E-07 0.17275612E-10-0.26574529E+04 0.73468280E+01-0.12027167E+04    4
!
! A. Burcat
CH3CH2O           T11/82O   1C   2H   5    0G   300.000  5000.0   1000.        1  
 0.60114346E+01 0.12165219E-01-0.40449604E-05 0.59076588E-09-0.30969595E-13    2
-0.49366992E+04-0.67901798E+01 0.17302504E+01 0.16908489E-01 0.39996221E-05    3
-0.13711180E-07 0.57643603E-11-0.32922483E+04 0.17336115E+02-0.20138288E+04    4
C2H4OH            T 4/83H   5C   2O   1    0G   300.000  5000.    1000.        1  
 0.75944014E+01 0.93229339E-02-0.30303854E-05 0.43216319E-09-0.21970039E-13    2
-0.57727852E+04-0.13955572E+02 0.14019508E+01 0.21543175E-01-0.22326512E-05    3
-0.14464092E-07 0.80488420E-11-0.38464519E+04 0.19148981E+02-0.25154820E+04    4
CH3CHOH           T 4/83C   2O   1H   5    0G   300.000  5000.    1000.        1  
 0.67665424E+01 0.11634436E-01-0.37790651E-05 0.53828875E-09-0.27315345E-13    2
-0.56092969E+04-0.93980442E+01 0.24813328E+01 0.16790036E-01 0.37755499E-05    3
-0.13923497E-07 0.60095193E-11-0.40120054E+04 0.14581622E+02-0.25172860E+04    4
CH                121286C   1H   1    0    0G   300.000  5000.000 1365.000    01
 2.21342975E+00 2.15716841E-03-5.73074267E-07 6.56809882E-11-2.67729778E-15    2
 7.09278928E+04 9.20307557E+00 3.50778851E+00-4.60811343E-04 1.32973337E-06    3
-5.15427373E-10 5.83223892E-14 7.04334301E+04 2.09978249E+00                   4
HOC2H4O2                C   2H   5O   3    0G   300.000  5000.000 1391.000    01
 1.00941573E+01 1.23879015E-02-3.73811683E-06 5.46874551E-10-3.09943951E-14    2
-2.37710522E+04-2.00956526E+01 4.44209543E+00 2.52880383E-02-1.51605275E-05    3
 5.24921198E-09-7.91470852E-13-2.17507126E+04 1.04122371E+01                   4
!~~~~~~~~~~~~~~~~~ Notes on THERMO data for species below ~~~~~~~~~~~~~~~~~~~
! 
! OH dHf adjusted to 8.91 kcal/mol (Ruscic et al., 2002)
! CH3 dHf appears to be consistant with Ruscic et al., 1999 already
!
! CH2OH thermo was fit directly to Johnson&Hudgens (1996) table,
! temperatures above 2000 K were extrapolated towards Cp_inf= 12.5 R:
! Cp/Cp_inf = a + (1-a)*exp(-(b/T)^n), where
! a = 0.4238, b = 696.36 K, n = 1.3807
!
! Equilibrium constants:
! C2H5OH = CH3 + CH2OH
! A            n            E (kcal/mol)
! 7.82635E+09 -1.93207E+00  8.92141E+04
!
! C2H5OH = C2H5 + OH
! A            n            E (kcal/mol)
! 1.21865E+08 -1.55343E+00  9.47874E+04
!
CH3               IU0702C  1.H  3.   0.   0.G   200.000  6000.000 1000.        1
 0.29781206E+01 0.57978520E-02-0.19755800E-05 0.30729790E-09-0.17917416E-13    2
 0.16509513E+05 0.47224799E+01 0.36571797E+01 0.21265979E-02 0.54583883E-05    3
-0.66181003E-08 0.24657074E-11 0.16422716E+05 0.16735354E+01 0.17643935E+05    4
CH2OH              JH/96C   1H   3O   1N   0G   250.000  3000.000  750.00      1
 0.37469103E+01 0.88646121E-02-0.42580722E-05 0.10088040E-08-0.94501561E-13    2
-0.36664824E+04 0.54281095E+01 0.46119792E+01-0.31203760E-02 0.35531680E-04    3
-0.49379398E-07 0.22027247E-10-0.36040734E+04 0.28351399E+01                   4
C2H5              T12/91C   2H   5    0    0G   200.000  6000.000 1000.        1
 0.42878814E+01 0.12433893E-01-0.44139119E-05 0.70654102E-09-0.42035136E-13    2
 0.12056455E+05 0.84602583E+00 0.43058580E+01-0.41833638E-02 0.49707270E-04    3
-0.59905874E-07 0.23048478E-10 0.12841714E+05 0.47100236E+01 0.14271225E+05    4
! Previous value prior to Ruscic update
!OH         8/12/99 THERMH   1O   1    0    0G   300.000  5000.000 1357.000    01
! 2.62599754E+00 1.31992406E-03-3.59724670E-07 4.25630800E-11-1.82048016E-15    2
! 4.12085374E+03 7.10667307E+00 3.43586219E+00 2.02235804E-04-1.13546412E-07    3
! 2.42445149E-10-7.43651031E-14 3.74321252E+03 2.45014127E+00                   4
OH                S 9/01O   1H   1    0    0G   200.000  6000.000 1000.        1
 2.86472886E+00 1.05650448E-03-2.59082758E-07 3.05218674E-11-1.33195876E-15    2
 3.68362875E+03 5.70164073E+00 4.12530561E+00-3.22544939E-03 6.52764691E-06    3
-5.79853643E-09 2.06237379E-12 3.34630913E+03-6.90432960E-01 4.51532273E+03    4
C2H5OH            L 8/88C   2H   6O   1    0G   200.000  6000.000 1000.        1
 0.65624365E+01 0.15204222E-01-0.53896795E-05 0.86225011E-09-0.51289787E-13    2
-0.31525621E+05-0.94730202E+01 0.48586957E+01-0.37401726E-02 0.69555378E-04    3
-0.88654796E-07 0.35168835E-10-0.29996132E+05 0.48018545E+01-0.28257829E+05    4
!
CH3OCH3           AK0904C   2H   6O   1    0G   270.000  3000.000  710.00      1
 8.30815546E-01 2.69173263E-02-1.38874777E-05 3.47515079E-09-3.41706784E-13    2
-2.34120975E+04 2.02174360E+01 5.68097447E+00-5.39434751E-03 6.49472750E-05    3
-8.05065318E-08 3.27474018E-11-2.39755455E+04-6.36955496E-01                   4
!
CH3OCH2    7/20/98 THERMC   2H   5O   1    0G   300.000  5000.000 1376.000    21
 8.17137842E+00 1.10086181E-02-3.82352277E-06 5.99637202E-10-3.50317513E-14    2
-3.41941605E+03-1.78650856E+01 2.91327415E+00 2.03364659E-02-9.59712342E-06    3
 2.07478525E-09-1.71343362E-13-1.18844240E+03 1.16066817E+01                   4
CH3OCO     4/20/99 THERMC   2H   3O   2    0G   300.000  5000.000 1362.000    21
 1.30877600E+01 4.53544950E-03-1.65096364E-06 2.67197277E-10-1.59576863E-14    2
-2.46616400E+04-3.27914051E+01 3.94199159E+00 2.43434884E-02-1.65595560E-05    3
 4.58537411E-09-3.31795708E-13-2.14404829E+04 1.66954362E+01                   4
CH3OCHO    4/20/99 THERMC   2H   4O   2    0G   300.000  5000.000 1686.000    21
 8.69123518E+00 1.15503122E-02-4.27782486E-06 7.02533059E-10-4.24333552E-14    2
-4.64364769E+04-1.89301478E+01 3.08839783E+00 2.03760048E-02-6.84777040E-06    3
-7.28186203E-10 5.62130216E-13-4.41855167E+04 1.25364719E+01                   4
OCHO       2/14/95 THERMC   1H   1O   2    0G   300.000  5000.000 1690.000    01
 6.12628782E+00 3.75602932E-03-1.42010352E-06 2.36429200E-10-1.44167651E-14    2
-2.17698466E+04-8.01574694E+00 1.35213452E+00 1.50082004E-02-1.09896141E-05    3
 3.73679840E-09-4.81014498E-13-2.02253647E+04 1.74373147E+01                   4
!
HE                120186HE  1               G  0300.00   5000.00  1000.00      1
 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-0.07453750E+04 0.09153489E+01 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-0.07453750E+04 0.09153488E+01                   4
!
! Species involved in low temperature DME oxidation (from LLNL database)
CH3OCH2O2  7/20/98 THERMC   2H   5O   3    0G   300.000  5000.000 1389.000    31
 1.24249729E+01 1.18705986E-02-4.07906532E-06 6.35310809E-10-3.69427867E-14    2
-2.29679238E+04-3.53740145E+01 2.21029612E+00 3.68877454E-02-2.82561555E-05    3
 1.15730533E-08-1.97130470E-12-1.94940940E+04 1.91463601E+01                   4
CH2OCH2O2H 7/20/98 THERMC   2H   5O   3    0G   300.000  5000.000 1393.000    41
 1.51191783E+01 9.23718883E-03-3.19127505E-06 4.99114678E-10-2.91162488E-14    2
-1.84114867E+04-4.85706618E+01 2.52895507E+00 4.24128290E-02-3.73406386E-05    3
 1.66639333E-08-2.96443312E-12-1.44293306E+04 1.76899251E+01                   4
CH3OCH2O2H 7/20/98 THERMC   2H   6O   3    0G   300.000  5000.000 1392.000    41
 1.49370964E+01 1.19465829E-02-4.12746359E-06 6.45422590E-10-3.76427939E-14    2
-4.11001068E+04-4.99552737E+01 1.19855761E+00 4.59060764E-02-3.66252420E-05    3
 1.49318970E-08-2.46057445E-12-3.65363161E+04 2.31339904E+01                   4
CH3OCH2O   2/ 9/96 THERMC   2H   5O   2    0G   300.000  5000.000 2012.000    21
 8.60261845E+00 1.35772195E-02-4.84661602E-06 7.77766193E-10-4.62633624E-14    2
-2.13762444E+04-1.75775023E+01 3.25889339E+00 2.22146359E-02-7.78556340E-06    3
-2.41484158E-10 4.51914496E-13-1.92377212E+04 1.23680069E+01                   4
CH3OCH2OH  2/ 9/96 THERMC   2H   6O   2    0G   300.000  5000.000 2014.000    31
 8.70981570e+00 1.53602372e-02-5.41003788e-06 8.60573446e-10-5.08819752e-14    2
-4.76607115e+04-1.80226702e+01 3.15851876e+00 2.44325751e-02-8.66984784e-06    3
-5.93319328e-11 4.36400003e-13-4.54488899e+04 1.30511235e+01                   4
O2CH2OCH2O2H 7/20/98 TRMC   2H   5O   5    0G   300.000  5000.000 1402.000    51
 1.92038046E+01 1.04394841E-02-3.60582939E-06 5.63792843E-10-3.28807214E-14    2
-3.79207055E+04-6.51847273E+01 1.99640551E+00 5.83226232E-02-5.53259778E-05    3
 2.59810540E-08-4.77141005E-12-3.27628742E+04 2.44215005E+01                   4
HO2CH2OCHO 8/31/99 THERMC   2H   4O   4    0G   300.000  5000.000 1387.000    41
 1.64584298E+01 8.52683511E-03-3.04113500E-06 4.85596908E-10-2.87316334E-14    2
-6.23959608E+04-5.38924139E+01 3.47935703E+00 4.02952392E-02-3.30109296E-05    3
 1.34360117E-08-2.18601580E-12-5.80629934E+04 1.52521392E+01                   4
OCH2OCHO   4/ 9/98 THERMC   2H   3O   3    0G   300.000  5000.000 1475.000    21
 1.20233916E+01 8.11262659E-03-2.91356462E-06 4.67340384E-10-2.77375525E-14    2
-4.33647231E+04-3.33691809E+01 5.19690837E+00 1.58839723E-02 3.53540547E-07    3
-6.10456923E-09 1.94661801E-12-4.02242792E+04 6.11645828E+00                   4
HOCH2OCO   8/31/99 THERMC   2H   3O   3    0G   300.000  5000.000 1603.000    31
 1.13737391E+01 8.17663898E-03-2.92034021E-06 4.66695616E-10-2.76276823E-14    2
-4.65575743E+04-2.86035265E+01 6.08180801E+00 1.28768359E-02 2.04419418E-06    3
-6.10154921E-09 1.79820559E-12-4.39526183E+04 2.54054449E+00                   4
HCOOH      7/ 1/99 THERMC   1H   2O   2    0G   300.000  5000.000 1376.000    11
 6.68733013e+00 5.14289368e-03-1.82238513e-06 2.89719163e-10-1.70892199e-14    2
-4.83995400e+04-1.13104798e+01 1.43548185e+00 1.63363016e-02-1.06257421e-05    3
 3.32132977e-09-4.02176103e-13-4.64616504e+04 1.72885798e+01                   4
END

\\
\\
\\  This is the tran file
\\
\\
AR                 0   136.500     3.330     0.000     0.000     0.000
AS                 0  1045.500     4.580     0.000     0.000     0.000 ! MEC
ASH                1   199.300     4.215     0.000     0.000     1.000 ! MEC
ASH2               2   229.600     4.180     0.000     0.000     1.000 ! MEC
C                  0    71.400     3.298     0.000     0.000     0.000 ! *
C2                 1    97.530     3.621     0.000     1.760     4.000
C2O                1   232.400     3.828     0.000     0.000     1.000 ! *
CN2                1   232.400     3.828     0.000     0.000     1.000 ! OIS
C2H                1   265.300     3.721     0.000     0.000     2.500 ! NMM
C2H2               1   265.300     3.721     0.000     0.000     2.500 ! NMM
C2H2OH             2   224.700     4.162     0.000     0.000     1.000 ! *
C2H3               2   265.300     3.721     0.000     0.000     1.000 ! NMM
C2H4               2   238.400     3.496     0.000     0.000     1.500 ! NMM
C2H5               2   247.500     4.350     0.000     0.000     1.500 ! NMM
C2H5OH             2   470.600     4.410     0.000     0.000     1.500 ! NMM
CH3CHOH            2   470.600     4.410     0.000     0.000     1.500
CH3CH2O            2   470.600     4.410     0.000     0.000     1.500 ! NMM
C2H4OH             2   470.600     4.410     0.000     0.000     1.500
PC2H4OH            2   470.600     4.410     0.000     0.000     1.500 ! NMM
SC2H4OH            2   470.600     4.410     0.000     0.000     1.500 ! NMM
HCOOH              2   470.600     4.410     0.000     0.000     1.500
C2H6               2   247.500     4.350     0.000     0.000     1.500 ! NMM
C2N                1   232.400     3.828     0.000     0.000     1.000 ! OIS
C2N2               1   349.000     4.361     0.000     0.000     1.000 ! OIS
C3H2               2   209.000     4.100     0.000     0.000     1.000 ! *
C3H2(S)            2   209.000     4.100     0.000     0.000     1.000 ! *
C3H3               1   324.800     4.290     0.000     0.000     1.000 ! NMM
C4H3               1   357.000     4.720     0.000     0.000     1.000 ! NMM
C3H4O              2   443.200     4.120     0.000     0.000     1.000 ! NMM
ACETONE            2   443.200     4.120     0.000     0.000     1.000
CH2CHCH2O          2   443.200     4.120     0.000     0.000     1.000
HOC2H4O2           2   443.200     4.120     0.000     0.000     1.000
CH3COCH2           2   443.200     4.120     0.000     0.000     1.000
CHCHCHO            2   443.200     4.120     0.000     0.000     1.000 ! NMM
HCCCHO             2   443.200     4.120     0.000     0.000     1.000 ! NMM
HCCCO              2   443.200     4.120     0.000     0.000     1.000 ! NMM
H2CCHCO            2   443.200     4.120     0.000     0.000     1.000 ! NMM
CH3CCO             2   443.200     4.120     0.000     0.000     1.000 ! NMM
CH3CHCO            2   443.200     4.120     0.000     0.000     1.000 ! NMM
CH2CHCO            2   443.200     4.120     0.000     0.000     1.000 ! NMM
C2H3CO             2   443.200     4.120     0.000     0.000     1.000 ! NMM
C2H5CHO            2   424.600     4.820     0.000     0.000     1.000 ! NMM
CH2CH2CHO          2   424.600     4.820     0.000     0.000     1.000 ! NMM
CH2CHCHO           2   424.600     4.820     0.000     0.000     1.000
C2H5CO             2   424.600     4.820     0.000     0.000     1.000 ! NMM
CH3COCH3           2   435.500     4.860     0.000     0.000     1.000 ! NMM
CH3COCH2           2   435.500     4.860     0.000     0.000     1.000 ! NMM
AC3H4              1   324.800     4.290     0.000     0.000     1.000 ! NMM
PC3H4              1   324.800     4.290     0.000     0.000     1.000 ! NMM
C3H4C              2   324.800     4.290     0.000     0.000     1.000 ! NMM
C3H6               2   307.800     4.140     0.000     0.000     1.000 ! NMM
C3H6OH             2   487.900     4.820     0.000     0.000     1.000 ! NMM
C3H6O              2   411.000     4.820     0.000     0.000     1.000 ! NMM
C3H5O              2   411.000     4.820     0.000     0.000     1.000 ! NMM
C3H7               2   303.400     4.810     0.000     0.000     1.000 ! NMM
C4H6               2   357.000     4.720     0.000     0.000     1.000 ! NMM
IC3H7              2   303.400     4.810     0.000     0.000     1.000 ! NMM
NC3H7              2   303.400     4.810     0.000     0.000     1.000 ! NMM
C3H8               2   303.400     4.810     0.000     0.000     1.000 ! NMM
C4H                1   357.000     4.720     0.000     0.000     1.000 ! NMM
C4H2               1   357.000     4.720     0.000     0.000     1.000 ! NMM
C4H2OH             2   224.700     4.162     0.000     0.000     1.000 ! *
CH3CHCCH           2   355.000     4.650     0.000     0.000     1.000 ! NMM
IC4H7              2   355.000     4.650     0.000     0.000     1.000 ! NMM
C4H7               2   355.000     4.650     0.000     0.000     1.000 ! NMM
C4H8               2   355.000     4.650     0.000     0.000     1.000 ! NMM
C4H8-1             2   355.000     4.650     0.000     0.000     1.000 ! NMM
C4H8-2             2   355.000     4.650     0.000     0.000     1.000 ! NMM
IC4H8              2   355.000     4.650     0.000     0.000     1.000 ! NMM
PC4H9              2   352.000     5.240     0.000     0.000     1.000 ! NMM
C4H9               2   352.000     5.240     0.000     0.000     1.000 ! NMM
SC4H9              2   352.000     5.240     0.000     0.000     1.000 ! NMM
TC4H9              2   352.000     5.240     0.000     0.000     1.000 ! NMM
IC4H9              2   352.000     5.240     0.000     0.000     1.000 ! NMM
C4H10              2   352.000     5.240     0.000     0.000     1.000 ! NMM
IC4H10             2   352.000     5.240     0.000     0.000     1.000 ! NMM
C5H2               1   408.000     5.200     0.000     0.000     1.000 ! NMM
C5H3               1   408.000     5.200     0.000     0.000     1.000 ! NMM
C5H5               1   408.000     5.200     0.000     0.000     1.000 ! NMM
C5H6               1   408.000     5.200     0.000     0.000     1.000 ! NMM
C5H7               2   408.000     5.200     0.000     0.000     1.000 ! NMM
CYC5H7             2   408.000     5.200     0.000     0.000     1.000 !
C5H8               2   408.000     5.200     0.000     0.000     1.000 ! NMM
C6H2               1   408.000     5.200     0.000     0.000     1.000 ! NMM
C6H4               2   412.300     5.349     0.000     0.000     1.000 ! JAM
C6H5               2   412.300     5.349     0.000     0.000     1.000 ! JAM
C6H5(L)            2   412.300     5.349     0.000     0.000     1.000 ! JAM
C6H5OH             2   450.000     5.500     0.000     0.000     1.000 ! NMM
C6H5O              2   450.000     5.500     0.000     0.000     1.000 ! JAM
C6H4O2             2   450.000     5.500     0.000     0.000     1.000
C5H4O              2   450.000     5.500     0.000     0.000     1.000 ! NMM
C5H4OH             2   450.000     5.500     0.000     0.000     1.000 ! NMM
C5H5O              2   450.000     5.500     0.000     0.000     1.000 ! NMM
C5H5OH             2   450.000     5.500     0.000     0.000     1.000 ! NMM
C6H5C2H            2   468.500     5.230     0.000     0.000     1.000 ! NMM
LC6H5              2   426.300     5.510     0.000     0.000     1.000
C6H6               2   468.500     5.230     0.000     10.30     1.000 ! NMM
C6H7               2   468.500     5.230     0.000     0.000     1.000 ! NMM
CYC6H7             2   468.500     5.230     0.000     0.000     1.000 ! NMM
CYC6H8             2   468.500     5.230     0.000     0.000     1.000 ! NMM
C6H5CH2            2   495.300     5.680     0.000     0.000     1.000 ! NMM
C6H5CH3            2   495.300     5.680     0.430     12.30     1.000 ! NMM
C6H5CO             2   622.400     5.530     0.000     0.000     1.000 ! NMM
C6H5CHO            2   622.400     5.530     0.000     0.000     1.000 ! NMM
C6H5CH2OH          2   622.400     5.530     0.000     0.000     1.000 ! NMM
OC6H4CH3           2   621.100     5.640     0.000     0.000     1.000 ! NMM
HOC6H4CH3          2   621.100     5.640     0.000     0.000     1.000 ! NMM
XYLYLENE           2   523.600     6.182     0.000     0.000     1.000
XYLYLRAD           2   523.600     6.182     0.000     0.000     1.000
C6H5C2H5           2   523.600     5.960     0.000     0.000     1.000 ! NMM
C6H9               2   426.300     5.510     0.000     0.000     1.000 ! NMM
C6H10              2   426.300     5.510     0.000     0.000     1.000 ! NMM
C8H14              2   494.000     6.170     0.000     0.000     1.000 ! NMM
IC8H14             2   494.000     6.170     0.000     0.000     1.000 ! NMM
C6H5C2H3           2   546.200     6.000     0.130     15.00     1.000 ! NMM
C6H5CHCH           2   546.200     6.000     0.000     0.000     1.000 ! NMM
C6H5CCH2           2   546.200     6.000     0.000     0.000     1.000
C6H5C2H            2   534.300     5.710     0.770     0.000     1.000 ! NMM
C6H4C2H3           2   546.200     6.000     0.000     0.000     1.000 ! NMM
C6H4C2H            2   534.300     5.710     0.000     0.000     1.000 ! NMM
C6H5CCO            2   588.200     5.940     0.000     0.000     1.0001
C10H7              2   630.400     6.180     0.000     0.000     1.000 ! NMM
C10H7O             2   630.400     6.180     0.000     0.000     1.000 ! NMM
C10H8              2   630.400     6.180     0.000     16.50     1.000 ! NMM
C10H9              2   630.400     6.180     0.000     0.000     1.000 ! NMM
C10H10             2   630.400     6.180     0.000     0.000     1.000 ! NMM
C10H7CH2           2   660.00      6.350     0.000     0.000     1.000 ! NMM
C10H7OH            2   663.45      6.362     0.000     0.000     1.000
C10H7CH3           2   660.0       6.350     0.000     0.000     1.000 ! NMM
FLRNTHN            2   812.3       7.170     0.000     0.000     1.000 ! NMM
ACEPHEN            2   812.3       7.170     0.000     0.000     1.000
ANTHRACN           2   772.0       6.960     0.000     25.40     1.000
CH3INDENE          2   625.0       6.150     0.000     0.000     1.000
CH3INDENYL         2   625.0       6.150     0.000     0.000     1.000
PHNTHRN            2   772.0       6.960     0.000     38.80     1.000
PYRENE             2   834.9       7.240     0.000     0.000     1.000 ! NMM
PYRENYL            2   834.9       7.240     0.000     0.000     1.000 ! NMM
DHPYRENE           2   834.9       7.240     0.000     0.000     1.000
BENZOAP            2   832.5       7.550     1.400     0.000     1.000 ! NMM
BENZOGHI           2   832.5       7.550     0.000     0.000     1.000
CPENTACD           2   832.5       7.550     0.000     0.000     1.000
INDENE             2   588.6       5.960     0.650     0.000     1.000 ! NMM
INDENYL            2   588.6       5.960     0.000     0.000     1.000 ! NMM
CH3FLRNE           2   712.6       6.890     0.000     0.000     1.000
CH3FLRNL           2   712.6       6.890     0.000     0.000     1.000
FLRENE             2   712.6       6.890     0.000     0.000     1.000
BIBENZYL           2   783.800     6.640     0.000     0.000     1.000 ! NMM
STILBENE           2   772.0       6.960     0.000     0.000     1.000 ! NMM
STILBNRD           2   772.0       6.960     0.000     0.000     1.000 ! NMM
ANTHRACN           2   772.0       6.960     0.000     0.000     1.000 ! NMM
DHANTHRN           2   772.0       6.960     0.000     0.000     1.000
CH                 1    80.000     2.750     0.000     0.000     0.000
CH2                1   144.000     3.800     0.000     0.000     0.000
CH2(S)             1   144.000     3.800     0.000     0.000     0.000
CH2CHCCH           2   373.700     4.790     0.000     0.000     1.000 ! NMM
CH2CHCCH2          2   373.700     4.790     0.000     0.000     1.000 ! NMM
CH3CCCH2           2   357.100     4.720     0.000     0.000     1.000
CH3CHCCH2          2   357.100     4.720     0.000     0.000     1.000
CH2CH2CCH          2   373.700     4.790     0.000     0.000     1.000 ! NMM
CH2CHCH2           2   316.000     4.220     0.000     0.000     1.000 ! NMM
AC3H5              2   316.000     4.220     0.000     0.000     1.000
PC3H5              2   316.000     4.220     0.000     0.000     1.000
CH2CHCHCH          2   357.100     4.720     0.000     0.000     1.000 ! NMM
CH2CHCHCH2         2   357.100     4.720     0.000     0.000     1.000 ! NMM
CH2CO              2   436.000     3.970     0.000     0.000     2.000
CH2O               2   498.000     3.590     0.000     0.000     2.000
HCOH               2   498.000     3.590     0.000     0.000     1.000
H2CO               2   498.000     3.590     0.000     0.000     2.000
CH2OH              2   417.000     3.690     1.700     0.000     2.000
CH2HCO             2   436.000     3.970     0.000     0.000     2.000
CHOCHO             1   440.200     4.010     0.000     0.000     2.000 ! NMM
CHOCO              1   440.200     4.010     0.000     0.000     2.000 ! NMM
CH3                1   144.000     3.800     0.000     0.000     0.000
CH3CC              2   252.000     4.760     0.000     0.000     1.000 ! JAM
CH3CHCCH           2   373.700     4.790     0.000     0.000     1.000 ! NMM
CH3CCCH2           2   357.100     4.720     0.000     0.000     1.000 ! NMM
CH3CCCH3           2   357.100     4.720     0.000     0.000     1.000 ! NMM
CH3CCH2            2   316.000     4.220     0.000     0.000     1.000 ! NMM
TC3H5              2   316.000     4.220     0.000     0.000     1.000
CH3CHCH            2   316.000     4.220     0.000     0.000     1.000 ! NMM
SC3H5              2   316.000     4.220     0.000     0.000     1.000
CH3CH2CCH          2   357.100     4.720     0.000     0.000     1.000 ! NMM
CH3HCO             2   436.000     3.970     0.000     0.000     2.000
CH3CO              2   436.000     3.970     0.000     0.000     2.000
CH3O               2   417.000     3.690     1.700     0.000     2.000
CH3OH              2   481.800     3.626     0.000     0.000     1.000 ! SVE
CH4                2   141.400     3.746     0.000     2.600    13.000
CH4O               2   417.000     3.690     1.700     0.000     2.000
CN                 1    75.000     3.856     0.000     0.000     1.000 ! OIS
CNC                1   232.400     3.828     0.000     0.000     1.000 ! OIS
CNN                1   232.400     3.828     0.000     0.000     1.000 ! OIS
CO                 1    98.100     3.650     0.000     1.950     1.800
CO2                1   244.000     3.763     0.000     2.650     2.100
F                  0    80.000     2.750     0.000     0.000     0.000
F2                 1   125.700     3.301     0.000     1.600     3.800
H                  0   145.000     2.050     0.000     0.000     0.000
GAH                1   335.500     4.240     0.000     0.000     1.000 ! MEC
H2C4O              2   357.000     5.180     0.000     0.000     1.000 ! JAM
H2                 1    38.000     2.920     0.000     0.790   280.000
H2CCC              2   265.300     3.721     0.000     0.000     1.000 ! *
H2CCC(S)           2   265.300     3.721     0.000     0.000     1.000 ! *
H2CCCH             2   252.000     4.760     0.000     0.000     1.000 ! JAM
H2CCCCH            2   357.100     4.720     0.000     0.000     1.000 ! NMM
H2CCCCH2           2   357.100     4.720     0.000     0.000     1.000 ! NMM
H2CCCCCH           1   408.000     5.200     0.000     0.000     1.000 ! NMM
H2CN               1   569.000     3.630     0.000     0.000     1.000 ! OS/JM
H2NO               2   116.700     3.492     0.000     0.000     1.000 ! JAM
H2O                2   572.400     2.605     1.844     0.000     4.000
H2O2               2   107.400     3.458     0.000     0.000     3.800
H2S                2   301.000     3.600     0.000     0.000     1.000 ! OIS
HC2N2              1   349.000     4.361     0.000     0.000     1.000 ! OIS
HCCHCCH            2   357.100     4.720     0.000     0.000     1.000 ! NMM
HCCO               2   150.000     2.500     0.000     0.000     1.000 ! *
HCCOH              2   436.000     3.970     0.000     0.000     2.000
HCCCHCCH           1   408.000     5.200     0.000     0.000     1.000 ! NMM
HCN                1   569.000     3.630     0.000     0.000     1.000 ! OIS
HCO                2   498.000     3.590     0.000     0.000     0.000
HCO+               1   498.000     3.590     0.000     0.000     0.000
HE                 0    10.200     2.576     0.000     0.000     0.000 ! *
HF                 1   330.000     3.148     1.920     2.460     1.000 ! SV/MEC
HF0                1   352.000     2.490     1.730     0.000     5.000
HF1                1   352.000     2.490     1.730     0.000     5.000
HF2                1   352.000     2.490     1.730     0.000     5.000
HF3                1   352.000     2.490     1.730     0.000     5.000
HF4                1   352.000     2.490     1.730     0.000     5.000
HF5                1   352.000     2.490     1.730     0.000     5.000
HF6                1   352.000     2.490     1.730     0.000     5.000
HF7                1   352.000     2.490     1.730     0.000     5.000
HF8                1   352.000     2.490     1.730     0.000     5.000
HCNO               2   232.400     3.828     0.000     0.000     1.000 ! JAM
HOCN               2   232.400     3.828     0.000     0.000     1.000 ! JAM
HNCO               2   232.400     3.828     0.000     0.000     1.000 ! OIS
HNNO               2   232.400     3.828     0.000     0.000     1.000 ! *
HNO                2   116.700     3.492     0.000     0.000     1.000 ! *
HNOH               2   116.700     3.492     0.000     0.000     1.000 ! JAM
HO2                2   107.400     3.458     0.000     0.000     1.000 ! *
HSO2               2   252.000     4.290     0.000     0.000     1.000 ! OIS
N                  0    71.400     3.298     0.000     0.000     0.000 ! *
N2                 1    97.530     3.621     0.000     1.760     4.000
N2H2               2    71.400     3.798     0.000     0.000     1.000 ! *
N2H3               2   200.000     3.900     0.000     0.000     1.000 ! *
N2H4               2   205.000     4.230     0.000     4.260     1.500
N2O                1   232.400     3.828     0.000     0.000     1.000 ! *
NCN                1   232.400     3.828     0.000     0.000     1.000 ! OIS
NCO                1   232.400     3.828     0.000     0.000     1.000 ! OIS
NH                 1    80.000     2.650     0.000     0.000     4.000
NH2                2    80.000     2.650     0.000     2.260     4.000
NH3                2   481.000     2.920     1.470     0.000    10.000
NNH                2    71.400     3.798     0.000     0.000     1.000 ! *
NO                 1    97.530     3.621     0.000     1.760     4.000
NCNO               2   232.400     3.828     0.000     0.000     1.000 ! OIS
NO2                2   200.000     3.500     0.000     0.000     1.000 ! *
O                  0    80.000     2.750     0.000     0.000     0.000
O2                 1   107.400     3.458     0.000     1.600     3.800
O3                 2   180.000     4.100     0.000     0.000     2.000
OH                 1    80.000     2.750     0.000     0.000     0.000
S                  0   847.000     3.839     0.000     0.000     0.000 ! OIS
S2                 1   847.000     3.900     0.000     0.000     1.000 ! OIS
SH                 1   847.000     3.900     0.000     0.000     1.000 ! OIS
SO                 1   301.000     3.993     0.000     0.000     1.000 ! OIS
SO2                2   252.000     4.290     0.000     0.000     1.000 ! OIS
SO3                2   378.400     4.175     0.000     0.000     1.000 ! OIS
SIH4               2   207.6       4.084     0.000     0.000     1.000 ! MEC
SIH3               2   170.3       3.943     0.000     0.000     1.000 ! MEC
SIH2               2   133.1       3.803     0.000     0.000     1.000 ! MEC
SIH                1    95.8       3.662     0.000     0.000     1.000 ! MEC
SI                 0  3036.        2.910     0.000     0.000     0.000 ! MEC
SI2H6              2   301.3       4.828     0.000     0.000     1.000 ! MEC
SI2H5              2   306.9       4.717     0.000     0.000     1.000 ! MEC
SI2H4              2   312.6       4.601     0.000     0.000     1.000 ! MEC
SI2H3              2   318.2       4.494     0.000     0.000     1.000 ! MEC
SI2H2              2   323.8       4.383     0.000     0.000     1.000 ! MEC
SI2                1  3036.        3.280     0.000     0.000     1.000 ! MEC
SI3                2  3036.        3.550     0.000     0.000     1.000 ! MEC
SIF3               2   309.6       4.359     0.000     0.000     1.000 ! MEC
SIF3NH2            2   231.0       4.975     0.000     0.000     1.000 ! MEC
SIF4               2   171.9       4.880     0.000     0.000     1.000 ! SVE
SIHF3              2   180.8       4.681     0.000     0.000     1.000 ! MEC
H2SISIH2           2   312.6       4.601     0.000     0.000     1.000 ! MEC
H3SISIH            2   312.6       4.601     0.000     0.000     1.000 ! MEC
SI3H8              2   331.2       5.562     0.000     0.000     1.000 ! MEC
ASH3               2   259.8       4.145     0.000     0.000     1.000 ! MEC
AS2                1   1045.5      5.510     0.000     0.000     1.000 ! MEC
GAME3              2   378.2       5.52      0.000     0.000     1.000 ! MEC
GAME2              2   675.8       5.22      0.000     0.000     1.000 ! MEC
GAME               2   972.7       4.92      0.000     0.000     1.000 ! MEC
GA                 0  2961.8       4.62      0.000     0.000     0.000 ! MEC
K                  0   850.        4.25      0.000     0.000     1.000 ! SINGH
KOH                2  1213.        4.52      0.000     0.000     1.000 ! SINGH
KO2                2  1213.        4.69      0.000     0.000     1.000 ! SINGH
KH                 1    93.3       3.542     0.000     0.000     1.000 ! SINGH
K+                 0   850.        4.25      0.000     0.000     1.000 ! SINGH
E                  0   850.        425.      0.000     0.000     1.000 ! SINGH
KCL                1  1989.        4.186     0.000     0.000     1.000 ! SINGH
CL                 0   130.8       3.613     0.000     0.000     1.000 ! SINGH
CL-                0   130.8       3.613     0.000     0.000     1.000 ! SINGH
HCL                1   344.7       3.339     1.084     0.000     1.000 ! SINGH
KO                 1   383.0       3.812     0.000     0.000     1.000 ! SINGH
CH3OCH3            2   329.400     4.624     0.000     0.000     1.000 !loc_est
CH3OCH2            2   329.400     4.624     0.000     0.000     1.000 !=CH3OCH3
CH3OCH2O           2   470.900     4.862     0.000     0.000     1.000 !loc_est
CH3OCHO            2   406.500     4.709     0.000     0.000     1.000 !loc_est
CH3OCO             2   406.500     4.709     0.000     0.000     1.000 !=CH3OCHO
OCHO               2   485.400     4.410     0.000     0.000     1.000 !loc_est
HOCH2O             2   481.800     3.626     0.000     0.000     1.000 !=CH3OH
CH3OCH2OH          2   329.400     4.624     0.000     0.000     1.000 !=CH3OCH3
OCH2OCHO           2   329.400     4.624     0.000     0.000     1.000 !=CH3OCH3
HOCH2OCO           2   329.400     4.624     0.000     0.000     1.000 !=CH3OCH3
CH3OCH2O2          2   329.400     4.624     0.000     0.000     1.000 !=CH3OCH3
CH2OCH2O2H         2   329.400     4.624     0.000     0.000     1.000 !=CH3OCH3
CH3OCH2O2H         2   329.400     4.624     0.000     0.000     1.000 !=CH3OCH3
HO2CH2OCHO         2   329.400     4.624     0.000     0.000     1.000 !=CH3OCH3
O2CH2OCH2O2H       2   329.400     4.624     0.000     0.000     1.000 !=CH3OCH3

#endif
