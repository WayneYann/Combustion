
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#if defined(BL_FORT_USE_UPPERCASE)
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
#elif defined(BL_FORT_USE_LOWERCASE)
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
#elif defined(BL_FORT_USE_UNDERSCORE)
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
#endif

/*function declarations */
extern "C" {
void molecularWeight(double * wt);
void gibbs(double * species, double * tc);
void helmholtz(double * species, double * tc);
void speciesInternalEnergy(double * species, double * tc);
void speciesEnthalpy(double * species, double * tc);
void speciesEntropy(double * species, double * tc);
void cp_R(double * species, double * tc);
void cv_R(double * species, double * tc);
void equilibriumConstants(double * kc, double * g_RT, double T);
void productionRate(double * wdot, double * sc, double T);
void progressRate(double * qdot, double * speciesConc, double T);
void progressRateFR(double * q_f, double * q_r, double * speciesConc, double T);
void CKINDX(int * iwrk, double *rwrk, int * mm, int * kk, int * ii, int * nfit );
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double * rval, int * kerr, int lenline);
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * rval, int * kerr, int lenline, int lenkray);
void CKSYME(int * kname, int * lenkname);
void CKSYMS(int * kname, int * lenkname);
void CKRP(int * ickwrk, double * rckwrk, double * ru, double * ruc, double * pa);
void CKPX(double * rho, double * T, double * x, int * iwrk, double *rwrk, double * P);
void CKPY(double * rho, double * T, double * y, int * iwrk, double *rwrk, double * P);
void CKPC(double * rho, double * T, double * c, int * iwrk, double *rwrk, double * P);
void CKRHOX(double * P, double * T, double * x, int * iwrk, double *rwrk, double * rho);
void CKRHOY(double * P, double * T, double * y, int * iwrk, double *rwrk, double * rho);
void CKRHOC(double * P, double * T, double * c, int * iwrk, double *rwrk, double * rho);
void CKWT(int * iwrk, double *rwrk, double * wt);
void CKMMWY(double * y, int * iwrk, double * rwrk, double * wtm);
void CKMMWX(double * x, int * iwrk, double * rwrk, double * wtm);
void CKMMWC(double * c, int * iwrk, double * rwrk, double * wtm);
void CKYTX(double * y, int * iwrk, double * rwrk, double * x);
void CKYTCP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * c);
void CKYTCR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * c);
void CKXTY(double * x, int * iwrk, double * rwrk, double * y);
void CKXTCP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * c);
void CKXTCR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * c);
void CKCTX(double * c, int * iwrk, double * rwrk, double * x);
void CKCTY(double * c, int * iwrk, double * rwrk, double * y);
void CKCPOR(double * T, int * iwrk, double * rwrk, double * cpor);
void CKHORT(double * T, int * iwrk, double * rwrk, double * hort);
void CKSOR(double * T, int * iwrk, double * rwrk, double * sor);
void CKCVML(double * T, int * iwrk, double * rwrk, double * cvml);
void CKCPML(double * T, int * iwrk, double * rwrk, double * cvml);
void CKUML(double * T, int * iwrk, double * rwrk, double * uml);
void CKHML(double * T, int * iwrk, double * rwrk, double * uml);
void CKGML(double * T, int * iwrk, double * rwrk, double * gml);
void CKAML(double * T, int * iwrk, double * rwrk, double * aml);
void CKSML(double * T, int * iwrk, double * rwrk, double * sml);
void CKCVMS(double * T, int * iwrk, double * rwrk, double * cvms);
void CKCPMS(double * T, int * iwrk, double * rwrk, double * cvms);
void CKUMS(double * T, int * iwrk, double * rwrk, double * ums);
void CKHMS(double * T, int * iwrk, double * rwrk, double * ums);
void CKGMS(double * T, int * iwrk, double * rwrk, double * gms);
void CKAMS(double * T, int * iwrk, double * rwrk, double * ams);
void CKSMS(double * T, int * iwrk, double * rwrk, double * sms);
void CKCPBL(double * T, double * x, int * iwrk, double * rwrk, double * cpbl);
void CKCPBS(double * T, double * y, int * iwrk, double * rwrk, double * cpbs);
void CKCVBL(double * T, double * x, int * iwrk, double * rwrk, double * cpbl);
void CKCVBS(double * T, double * y, int * iwrk, double * rwrk, double * cpbs);
void CKHBML(double * T, double * x, int * iwrk, double * rwrk, double * hbml);
void CKHBMS(double * T, double * y, int * iwrk, double * rwrk, double * hbms);
void CKUBML(double * T, double * x, int * iwrk, double * rwrk, double * ubml);
void CKUBMS(double * T, double * y, int * iwrk, double * rwrk, double * ubms);
void CKSBML(double * P, double * T, double * x, int * iwrk, double * rwrk, double * sbml);
void CKSBMS(double * P, double * T, double * y, int * iwrk, double * rwrk, double * sbms);
void CKGBML(double * P, double * T, double * x, int * iwrk, double * rwrk, double * gbml);
void CKGBMS(double * P, double * T, double * y, int * iwrk, double * rwrk, double * gbms);
void CKABML(double * P, double * T, double * x, int * iwrk, double * rwrk, double * abml);
void CKABMS(double * P, double * T, double * y, int * iwrk, double * rwrk, double * abms);
void CKWC(double * T, double * C, int * iwrk, double *rwrk, double * wdot);
void CKWYP(double * P, double * T, double * y, int * iwrk, double *rwrk, double * wdot);
void CKWXP(double * P, double * T, double * x, int * iwrk, double *rwrk, double * wdot);
void CKWYR(double * rho, double * T, double * y, int * iwrk, double *rwrk, double * wdot);
void CKWXR(double * rho, double * T, double * x, int * iwrk, double *rwrk, double * wdot);
void CKQC(double * T, double * C, int * iwrk, double *rwrk, double * qdot);
void CKKFKR(double * P, double * T, double * x, int * iwrk, double *rwrk, double * q_f, double * q_r);
void CKQYP(double * P, double * T, double * y, int * iwrk, double *rwrk, double * qdot);
void CKQXP(double * P, double * T, double * x, int * iwrk, double *rwrk, double * qdot);
void CKQYR(double * rho, double * T, double * y, int * iwrk, double *rwrk, double * qdot);
void CKQXR(double * rho, double * T, double * x, int * iwrk, double *rwrk, double * qdot);
void CKNU(int * kdim, int * iwrk, double *rwrk, int * nuki);
void CKNCF(int * mdim, int * iwrk, double *rwrk, int * ncf);
void CKABE(int * iwrk, double *rwrk, double * a, double * b, double * e );
void CKEQC(double * T, double * C , int * iwrk, double *rwrk, double * eqcon );
void CKEQYP(double * P, double * T, double * y, int * iwrk, double *rwrk, double * eqcon);
void CKEQXP(double * P, double * T, double * x, int * iwrk, double *rwrk, double * eqcon);
void CKEQYR(double * rho, double * T, double * y, int * iwrk, double *rwrk, double * eqcon);
void CKEQXR(double * rho, double * T, double * x, int * iwrk, double *rwrk, double * eqcon);
int  feeytt_(double * e, double * y, int * iwrk, double *rwrk, double * t);
void fephity_(double * phi, int * iwrk, double *rwrk, double * y);
void feytphi_(double * y, int * iwrk, double *rwrk, double * phi);
void fectyr_(double * c, double * rho, int * iwrk, double *rwrk, double * y);
void fecvrhs_(double * time, double * phi, double * phidot, double * rckwrk, int * ickwrk);
int fecvdim_();
void fezndrhs_(double * time, double * z, double * zdot, double * rckwrk, int * ickwrk);
int feznddim_();
char* femechfile_();
char* fesymname_(int sn);
int fesymnum_(const char* s1);
}


/*A few mechanism parameters */
void CKINDX(int * iwrk, double * rwrk, int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 3;
    *kk = 15;
    *ii = 58;
    *nfit = -1; /*Why do you need this anyway ?  */
}


/*Dummy ckinit */
void fginit_(int * leniwk, int * lenrwk, int * lencwk, int * linc, int * lout, int * ickwrk, double * rckwrk, char * cckwrk )
{
    if ((*lout) != 0) {
        /* printf(" ***       Congratulations       *** \n"); */
        /* printf(" * You are using the Fuego Library * \n"); */
        /* printf(" *****    Say NO to cklib.f    ***** \n"); */
    }
}


/* ckxnum... for parsing strings  */
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double * rval, int * kerr, int lenline )
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
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * rval, int * kerr, int lenline, int lenkray)
{
    /*Not done yet ... */
}


/* Returns the char strings of element names */
void CKSYME(int * kname, int * plenkname)
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*3; i++) {
        kname[i] = ' ';
    }

    /* O  */
    kname[ 0*lenkname + 0 ] = 'O';
    kname[ 0*lenkname + 1 ] = ' ';

    /* H  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = ' ';

    /* N  */
    kname[ 2*lenkname + 0 ] = 'N';
    kname[ 2*lenkname + 1 ] = ' ';
}


/* Returns the char strings of species names */
void CKSYMS(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*15; i++) {
        kname[i] = ' ';
    }

    /* H  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = ' ';

    /* O  */
    kname[ 1*lenkname + 0 ] = 'O';
    kname[ 1*lenkname + 1 ] = ' ';

    /* OH  */
    kname[ 2*lenkname + 0 ] = 'O';
    kname[ 2*lenkname + 1 ] = 'H';
    kname[ 2*lenkname + 2 ] = ' ';

    /* H2  */
    kname[ 3*lenkname + 0 ] = 'H';
    kname[ 3*lenkname + 1 ] = '2';
    kname[ 3*lenkname + 2 ] = ' ';

    /* O2  */
    kname[ 4*lenkname + 0 ] = 'O';
    kname[ 4*lenkname + 1 ] = '2';
    kname[ 4*lenkname + 2 ] = ' ';

    /* HO2  */
    kname[ 5*lenkname + 0 ] = 'H';
    kname[ 5*lenkname + 1 ] = 'O';
    kname[ 5*lenkname + 2 ] = '2';
    kname[ 5*lenkname + 3 ] = ' ';

    /* H2O  */
    kname[ 6*lenkname + 0 ] = 'H';
    kname[ 6*lenkname + 1 ] = '2';
    kname[ 6*lenkname + 2 ] = 'O';
    kname[ 6*lenkname + 3 ] = ' ';

    /* H2O2  */
    kname[ 7*lenkname + 0 ] = 'H';
    kname[ 7*lenkname + 1 ] = '2';
    kname[ 7*lenkname + 2 ] = 'O';
    kname[ 7*lenkname + 3 ] = '2';
    kname[ 7*lenkname + 4 ] = ' ';

    /* NO  */
    kname[ 8*lenkname + 0 ] = 'N';
    kname[ 8*lenkname + 1 ] = 'O';
    kname[ 8*lenkname + 2 ] = ' ';

    /* NO2  */
    kname[ 9*lenkname + 0 ] = 'N';
    kname[ 9*lenkname + 1 ] = 'O';
    kname[ 9*lenkname + 2 ] = '2';
    kname[ 9*lenkname + 3 ] = ' ';

    /* N2O  */
    kname[ 10*lenkname + 0 ] = 'N';
    kname[ 10*lenkname + 1 ] = '2';
    kname[ 10*lenkname + 2 ] = 'O';
    kname[ 10*lenkname + 3 ] = ' ';

    /* NH  */
    kname[ 11*lenkname + 0 ] = 'N';
    kname[ 11*lenkname + 1 ] = 'H';
    kname[ 11*lenkname + 2 ] = ' ';

    /* N  */
    kname[ 12*lenkname + 0 ] = 'N';
    kname[ 12*lenkname + 1 ] = ' ';

    /* NNH  */
    kname[ 13*lenkname + 0 ] = 'N';
    kname[ 13*lenkname + 1 ] = 'N';
    kname[ 13*lenkname + 2 ] = 'H';
    kname[ 13*lenkname + 3 ] = ' ';

    /* N2  */
    kname[ 14*lenkname + 0 ] = 'N';
    kname[ 14*lenkname + 1 ] = '2';
    kname[ 14*lenkname + 2 ] = ' ';

}


/* Returns R, Rc, Patm */
void CKRP(int * ickwrk, double * rckwrk, double * ru, double * ruc, double * pa)
{
     *ru  = 8.314e+07; 
     *ruc = 1.987; 
     *pa  = 1.01325e+06; 
}


/*Compute P = rhoRT/W(x) */
void CKPX(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * P)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*17.007370; /*OH */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*31.998800; /*O2 */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*30.006100; /*NO */
    XW += x[9]*46.005500; /*NO2 */
    XW += x[10]*44.012800; /*N2O */
    XW += x[11]*15.014670; /*NH */
    XW += x[12]*14.006700; /*N */
    XW += x[13]*29.021370; /*NNH */
    XW += x[14]*28.013400; /*N2 */
    *P = *rho * 8.314e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void CKPY(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]/1.007970; /*H */
    YOW += y[1]/15.999400; /*O */
    YOW += y[2]/17.007370; /*OH */
    YOW += y[3]/2.015940; /*H2 */
    YOW += y[4]/31.998800; /*O2 */
    YOW += y[5]/33.006770; /*HO2 */
    YOW += y[6]/18.015340; /*H2O */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/30.006100; /*NO */
    YOW += y[9]/46.005500; /*NO2 */
    YOW += y[10]/44.012800; /*N2O */
    YOW += y[11]/15.014670; /*NH */
    YOW += y[12]/14.006700; /*N */
    YOW += y[13]/29.021370; /*NNH */
    YOW += y[14]/28.013400; /*N2 */
    *P = *rho * 8.314e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(c) */
void CKPC(double * rho, double * T, double * c, int * iwrk, double * rwrk, double * P)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*1.007970; /*H */
    W += c[1]*15.999400; /*O */
    W += c[2]*17.007370; /*OH */
    W += c[3]*2.015940; /*H2 */
    W += c[4]*31.998800; /*O2 */
    W += c[5]*33.006770; /*HO2 */
    W += c[6]*18.015340; /*H2O */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*30.006100; /*NO */
    W += c[9]*46.005500; /*NO2 */
    W += c[10]*44.012800; /*N2O */
    W += c[11]*15.014670; /*NH */
    W += c[12]*14.006700; /*N */
    W += c[13]*29.021370; /*NNH */
    W += c[14]*28.013400; /*N2 */

    for (id = 0; id < 15; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.314e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double * P, double * T, double * x, int * iwrk, double * rwrk, double * rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*17.007370; /*OH */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*31.998800; /*O2 */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*30.006100; /*NO */
    XW += x[9]*46.005500; /*NO2 */
    XW += x[10]*44.012800; /*N2O */
    XW += x[11]*15.014670; /*NH */
    XW += x[12]*14.006700; /*N */
    XW += x[13]*29.021370; /*NNH */
    XW += x[14]*28.013400; /*N2 */
    *rho = *P * XW / (8.314e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double * P, double * T, double * y, int * iwrk, double * rwrk, double * rho)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]/1.007970; /*H */
    YOW += y[1]/15.999400; /*O */
    YOW += y[2]/17.007370; /*OH */
    YOW += y[3]/2.015940; /*H2 */
    YOW += y[4]/31.998800; /*O2 */
    YOW += y[5]/33.006770; /*HO2 */
    YOW += y[6]/18.015340; /*H2O */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/30.006100; /*NO */
    YOW += y[9]/46.005500; /*NO2 */
    YOW += y[10]/44.012800; /*N2O */
    YOW += y[11]/15.014670; /*NH */
    YOW += y[12]/14.006700; /*N */
    YOW += y[13]/29.021370; /*NNH */
    YOW += y[14]/28.013400; /*N2 */
    *rho = *P / (8.314e+07 * (*T) * YOW); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(c)/(R*T) */
void CKRHOC(double * P, double * T, double * c, int * iwrk, double * rwrk, double * rho)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*1.007970; /*H */
    W += c[1]*15.999400; /*O */
    W += c[2]*17.007370; /*OH */
    W += c[3]*2.015940; /*H2 */
    W += c[4]*31.998800; /*O2 */
    W += c[5]*33.006770; /*HO2 */
    W += c[6]*18.015340; /*H2O */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*30.006100; /*NO */
    W += c[9]*46.005500; /*NO2 */
    W += c[10]*44.012800; /*N2O */
    W += c[11]*15.014670; /*NH */
    W += c[12]*14.006700; /*N */
    W += c[13]*29.021370; /*NNH */
    W += c[14]*28.013400; /*N2 */

    for (id = 0; id < 15; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.314e+07); /*rho = PW/(R*T) */

    return;
}


/*get molecular weight for all species */
void CKWT(int * iwrk, double * rwrk, double * wt)
{
    molecularWeight(wt);
}


/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWY(double *y, int * iwrk, double * rwrk, double * wtm)
{
    double YOW = 0;/* see Eq 3 in CK Manual */
    YOW += y[0]/1.007970; /*H */
    YOW += y[1]/15.999400; /*O */
    YOW += y[2]/17.007370; /*OH */
    YOW += y[3]/2.015940; /*H2 */
    YOW += y[4]/31.998800; /*O2 */
    YOW += y[5]/33.006770; /*HO2 */
    YOW += y[6]/18.015340; /*H2O */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/30.006100; /*NO */
    YOW += y[9]/46.005500; /*NO2 */
    YOW += y[10]/44.012800; /*N2O */
    YOW += y[11]/15.014670; /*NH */
    YOW += y[12]/14.006700; /*N */
    YOW += y[13]/29.021370; /*NNH */
    YOW += y[14]/28.013400; /*N2 */
    *wtm = 1.0 / YOW;

    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWX(double *x, int * iwrk, double * rwrk, double * wtm)
{
    double XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*17.007370; /*OH */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*31.998800; /*O2 */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*30.006100; /*NO */
    XW += x[9]*46.005500; /*NO2 */
    XW += x[10]*44.012800; /*N2O */
    XW += x[11]*15.014670; /*NH */
    XW += x[12]*14.006700; /*N */
    XW += x[13]*29.021370; /*NNH */
    XW += x[14]*28.013400; /*N2 */
    *wtm = XW;

    return;
}


/*given c[species]: molar concentration */
/*returns mean molecular weight (gm/mole) */
void CKMMWC(double *c, int * iwrk, double * rwrk, double * wtm)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*1.007970; /*H */
    W += c[1]*15.999400; /*O */
    W += c[2]*17.007370; /*OH */
    W += c[3]*2.015940; /*H2 */
    W += c[4]*31.998800; /*O2 */
    W += c[5]*33.006770; /*HO2 */
    W += c[6]*18.015340; /*H2O */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*30.006100; /*NO */
    W += c[9]*46.005500; /*NO2 */
    W += c[10]*44.012800; /*N2O */
    W += c[11]*15.014670; /*NH */
    W += c[12]*14.006700; /*N */
    W += c[13]*29.021370; /*NNH */
    W += c[14]*28.013400; /*N2 */

    for (id = 0; id < 15; ++id) {
        sumC += c[id];
    }
    /* CK provides no guard against divison by zero */
    *wtm = W/sumC;

    return;
}


/*convert y[species] (mass fracs) to x[species] (mole fracs) */
void CKYTX(double * y, int * iwrk, double * rwrk, double * x)
{
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/1.007970; /*H */
    YOW += y[1]/15.999400; /*O */
    YOW += y[2]/17.007370; /*OH */
    YOW += y[3]/2.015940; /*H2 */
    YOW += y[4]/31.998800; /*O2 */
    YOW += y[5]/33.006770; /*HO2 */
    YOW += y[6]/18.015340; /*H2O */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/30.006100; /*NO */
    YOW += y[9]/46.005500; /*NO2 */
    YOW += y[10]/44.012800; /*N2O */
    YOW += y[11]/15.014670; /*NH */
    YOW += y[12]/14.006700; /*N */
    YOW += y[13]/29.021370; /*NNH */
    YOW += y[14]/28.013400; /*N2 */
    /*Now compute conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(15.999400*YOW); 
    x[2] = y[2]/(17.007370*YOW); 
    x[3] = y[3]/(2.015940*YOW); 
    x[4] = y[4]/(31.998800*YOW); 
    x[5] = y[5]/(33.006770*YOW); 
    x[6] = y[6]/(18.015340*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(30.006100*YOW); 
    x[9] = y[9]/(46.005500*YOW); 
    x[10] = y[10]/(44.012800*YOW); 
    x[11] = y[11]/(15.014670*YOW); 
    x[12] = y[12]/(14.006700*YOW); 
    x[13] = y[13]/(29.021370*YOW); 
    x[14] = y[14]/(28.013400*YOW); 

    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * c)
{
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/1.007970; /*H */
    YOW += y[1]/15.999400; /*O */
    YOW += y[2]/17.007370; /*OH */
    YOW += y[3]/2.015940; /*H2 */
    YOW += y[4]/31.998800; /*O2 */
    YOW += y[5]/33.006770; /*HO2 */
    YOW += y[6]/18.015340; /*H2O */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/30.006100; /*NO */
    YOW += y[9]/46.005500; /*NO2 */
    YOW += y[10]/44.012800; /*N2O */
    YOW += y[11]/15.014670; /*NH */
    YOW += y[12]/14.006700; /*N */
    YOW += y[13]/29.021370; /*NNH */
    YOW += y[14]/28.013400; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.314e+07 * (*T)); 
    /*Now compute conversion */
    c[0] = PWORT * y[0]/1.007970; 
    c[1] = PWORT * y[1]/15.999400; 
    c[2] = PWORT * y[2]/17.007370; 
    c[3] = PWORT * y[3]/2.015940; 
    c[4] = PWORT * y[4]/31.998800; 
    c[5] = PWORT * y[5]/33.006770; 
    c[6] = PWORT * y[6]/18.015340; 
    c[7] = PWORT * y[7]/34.014740; 
    c[8] = PWORT * y[8]/30.006100; 
    c[9] = PWORT * y[9]/46.005500; 
    c[10] = PWORT * y[10]/44.012800; 
    c[11] = PWORT * y[11]/15.014670; 
    c[12] = PWORT * y[12]/14.006700; 
    c[13] = PWORT * y[13]/29.021370; 
    c[14] = PWORT * y[14]/28.013400; 

    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * c)
{
    /*See Eq 8 (Temperature not used) */
    c[0] = (*rho) * y[0]/1.007970; 
    c[1] = (*rho) * y[1]/15.999400; 
    c[2] = (*rho) * y[2]/17.007370; 
    c[3] = (*rho) * y[3]/2.015940; 
    c[4] = (*rho) * y[4]/31.998800; 
    c[5] = (*rho) * y[5]/33.006770; 
    c[6] = (*rho) * y[6]/18.015340; 
    c[7] = (*rho) * y[7]/34.014740; 
    c[8] = (*rho) * y[8]/30.006100; 
    c[9] = (*rho) * y[9]/46.005500; 
    c[10] = (*rho) * y[10]/44.012800; 
    c[11] = (*rho) * y[11]/15.014670; 
    c[12] = (*rho) * y[12]/14.006700; 
    c[13] = (*rho) * y[13]/29.021370; 
    c[14] = (*rho) * y[14]/28.013400; 

    return;
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
void CKXTY(double * x, int * iwrk, double * rwrk, double * y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*17.007370; /*OH */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*31.998800; /*O2 */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*30.006100; /*NO */
    XW += x[9]*46.005500; /*NO2 */
    XW += x[10]*44.012800; /*N2O */
    XW += x[11]*15.014670; /*NH */
    XW += x[12]*14.006700; /*N */
    XW += x[13]*29.021370; /*NNH */
    XW += x[14]*28.013400; /*N2 */
    /*Now compute conversion */
    y[0] = x[0]*1.007970/XW; 
    y[1] = x[1]*15.999400/XW; 
    y[2] = x[2]*17.007370/XW; 
    y[3] = x[3]*2.015940/XW; 
    y[4] = x[4]*31.998800/XW; 
    y[5] = x[5]*33.006770/XW; 
    y[6] = x[6]*18.015340/XW; 
    y[7] = x[7]*34.014740/XW; 
    y[8] = x[8]*30.006100/XW; 
    y[9] = x[9]*46.005500/XW; 
    y[10] = x[10]*44.012800/XW; 
    y[11] = x[11]*15.014670/XW; 
    y[12] = x[12]*14.006700/XW; 
    y[13] = x[13]*29.021370/XW; 
    y[14] = x[14]*28.013400/XW; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.314e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 15; ++id) {
        c[id] = x[id]*PORT;
    }

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * c)
{
    int id; /*loop counter */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*17.007370; /*OH */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*31.998800; /*O2 */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*30.006100; /*NO */
    XW += x[9]*46.005500; /*NO2 */
    XW += x[10]*44.012800; /*N2O */
    XW += x[11]*15.014670; /*NH */
    XW += x[12]*14.006700; /*N */
    XW += x[13]*29.021370; /*NNH */
    XW += x[14]*28.013400; /*N2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 15; ++id) {
        c[id] = x[id]*ROW;
    }

    return;
}


/*convert c[species] (molar conc) to x[species] (mole fracs) */
void CKCTX(double * c, int * iwrk, double * rwrk, double * x)
{
    int id; /*loop counter */
    double sumC = 0; 

    /*compute sum of c  */
    for (id = 0; id < 15; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    for (id = 0; id < 15; ++id) {
        x[id] = c[id]/sumC;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
void CKCTY(double * c, int * iwrk, double * rwrk, double * y)
{
    double CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*1.007970; /*H */
    CW += c[1]*15.999400; /*O */
    CW += c[2]*17.007370; /*OH */
    CW += c[3]*2.015940; /*H2 */
    CW += c[4]*31.998800; /*O2 */
    CW += c[5]*33.006770; /*HO2 */
    CW += c[6]*18.015340; /*H2O */
    CW += c[7]*34.014740; /*H2O2 */
    CW += c[8]*30.006100; /*NO */
    CW += c[9]*46.005500; /*NO2 */
    CW += c[10]*44.012800; /*N2O */
    CW += c[11]*15.014670; /*NH */
    CW += c[12]*14.006700; /*N */
    CW += c[13]*29.021370; /*NNH */
    CW += c[14]*28.013400; /*N2 */
    /*Now compute conversion */
    y[0] = c[0]*1.007970/CW; 
    y[1] = c[1]*15.999400/CW; 
    y[2] = c[2]*17.007370/CW; 
    y[3] = c[3]*2.015940/CW; 
    y[4] = c[4]*31.998800/CW; 
    y[5] = c[5]*33.006770/CW; 
    y[6] = c[6]*18.015340/CW; 
    y[7] = c[7]*34.014740/CW; 
    y[8] = c[8]*30.006100/CW; 
    y[9] = c[9]*46.005500/CW; 
    y[10] = c[10]*44.012800/CW; 
    y[11] = c[11]*15.014670/CW; 
    y[12] = c[12]*14.006700/CW; 
    y[13] = c[13]*29.021370/CW; 
    y[14] = c[14]*28.013400/CW; 

    return;
}


/*get Cp/R as a function of T  */
/*for all species (Eq 19) */
void CKCPOR(double *T, int * iwrk, double * rwrk, double * cpor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpor, tc);
}


/*get H/RT as a function of T  */
/*for all species (Eq 20) */
void CKHORT(double *T, int * iwrk, double * rwrk, double * hort)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEnthalpy(hort, tc);
}


/*get S/R as a function of T  */
/*for all species (Eq 21) */
void CKSOR(double *T, int * iwrk, double * rwrk, double * sor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sor, tc);
}


/*get specific heat at constant volume as a function  */
/*of T for all species (molar units) */
void CKCVML(double *T, int * iwrk, double * rwrk, double * cvml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 15; ++id) {
        cvml[id] *= 8.314e+07;
    }
}


/*get specific heat at constant pressure as a  */
/*function of T for all species (molar units) */
void CKCPML(double *T, int * iwrk, double * rwrk, double * cpml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 15; ++id) {
        cpml[id] *= 8.314e+07;
    }
}


/*get internal energy as a function  */
/*of T for all species (molar units) */
void CKUML(double *T, int * iwrk, double * rwrk, double * uml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 15; ++id) {
        uml[id] *= RT;
    }
}


/*get enthalpy as a function  */
/*of T for all species (molar units) */
void CKHML(double *T, int * iwrk, double * rwrk, double * hml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 15; ++id) {
        hml[id] *= RT;
    }
}


/*get standard-state Gibbs energy as a function  */
/*of T for all species (molar units) */
void CKGML(double *T, int * iwrk, double * rwrk, double * gml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    gibbs(gml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 15; ++id) {
        gml[id] *= RT;
    }
}


/*get standard-state Helmholtz free energy as a  */
/*function of T for all species (molar units) */
void CKAML(double *T, int * iwrk, double * rwrk, double * aml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    helmholtz(aml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 15; ++id) {
        aml[id] *= RT;
    }
}


/*Returns the standard-state entropies in molar units */
void CKSML(double *T, int * iwrk, double * rwrk, double * sml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 15; ++id) {
        sml[id] *= 8.314e+07;
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
void CKCVMS(double *T, int * iwrk, double * rwrk, double * cvms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 82482613.569848; /*H */
    cvms[1] *= 5196444.866683; /*O */
    cvms[2] *= 4888468.940230; /*OH */
    cvms[3] *= 41241306.784924; /*H2 */
    cvms[4] *= 2598222.433341; /*O2 */
    cvms[5] *= 2518877.187922; /*HO2 */
    cvms[6] *= 4614955.920899; /*H2O */
    cvms[7] *= 2444234.470115; /*H2O2 */
    cvms[8] *= 2770769.943445; /*NO */
    cvms[9] *= 1807175.229049; /*NO2 */
    cvms[10] *= 1888995.928457; /*N2O */
    cvms[11] *= 5537251.234959; /*NH */
    cvms[12] *= 5935730.757423; /*N */
    cvms[13] *= 2864785.501167; /*NNH */
    cvms[14] *= 2967865.378712; /*N2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
void CKCPMS(double *T, int * iwrk, double * rwrk, double * cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 82482613.569848; /*H */
    cpms[1] *= 5196444.866683; /*O */
    cpms[2] *= 4888468.940230; /*OH */
    cpms[3] *= 41241306.784924; /*H2 */
    cpms[4] *= 2598222.433341; /*O2 */
    cpms[5] *= 2518877.187922; /*HO2 */
    cpms[6] *= 4614955.920899; /*H2O */
    cpms[7] *= 2444234.470115; /*H2O2 */
    cpms[8] *= 2770769.943445; /*NO */
    cpms[9] *= 1807175.229049; /*NO2 */
    cpms[10] *= 1888995.928457; /*N2O */
    cpms[11] *= 5537251.234959; /*NH */
    cpms[12] *= 5935730.757423; /*N */
    cpms[13] *= 2864785.501167; /*NNH */
    cpms[14] *= 2967865.378712; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double *T, int * iwrk, double * rwrk, double * ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    ums[0] *= RT/1.007970; /*H */
    ums[1] *= RT/15.999400; /*O */
    ums[2] *= RT/17.007370; /*OH */
    ums[3] *= RT/2.015940; /*H2 */
    ums[4] *= RT/31.998800; /*O2 */
    ums[5] *= RT/33.006770; /*HO2 */
    ums[6] *= RT/18.015340; /*H2O */
    ums[7] *= RT/34.014740; /*H2O2 */
    ums[8] *= RT/30.006100; /*NO */
    ums[9] *= RT/46.005500; /*NO2 */
    ums[10] *= RT/44.012800; /*N2O */
    ums[11] *= RT/15.014670; /*NH */
    ums[12] *= RT/14.006700; /*N */
    ums[13] *= RT/29.021370; /*NNH */
    ums[14] *= RT/28.013400; /*N2 */
}


/*Returns enthalpy in mass units (Eq 27.) */
void CKHMS(double *T, int * iwrk, double * rwrk, double * hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    speciesEnthalpy(hms, tc);
    hms[0] *= RT/1.007970; /*H */
    hms[1] *= RT/15.999400; /*O */
    hms[2] *= RT/17.007370; /*OH */
    hms[3] *= RT/2.015940; /*H2 */
    hms[4] *= RT/31.998800; /*O2 */
    hms[5] *= RT/33.006770; /*HO2 */
    hms[6] *= RT/18.015340; /*H2O */
    hms[7] *= RT/34.014740; /*H2O2 */
    hms[8] *= RT/30.006100; /*NO */
    hms[9] *= RT/46.005500; /*NO2 */
    hms[10] *= RT/44.012800; /*N2O */
    hms[11] *= RT/15.014670; /*NH */
    hms[12] *= RT/14.006700; /*N */
    hms[13] *= RT/29.021370; /*NNH */
    hms[14] *= RT/28.013400; /*N2 */
}


/*Returns gibbs in mass units (Eq 31.) */
void CKGMS(double *T, int * iwrk, double * rwrk, double * gms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    gibbs(gms, tc);
    gms[0] *= RT/1.007970; /*H */
    gms[1] *= RT/15.999400; /*O */
    gms[2] *= RT/17.007370; /*OH */
    gms[3] *= RT/2.015940; /*H2 */
    gms[4] *= RT/31.998800; /*O2 */
    gms[5] *= RT/33.006770; /*HO2 */
    gms[6] *= RT/18.015340; /*H2O */
    gms[7] *= RT/34.014740; /*H2O2 */
    gms[8] *= RT/30.006100; /*NO */
    gms[9] *= RT/46.005500; /*NO2 */
    gms[10] *= RT/44.012800; /*N2O */
    gms[11] *= RT/15.014670; /*NH */
    gms[12] *= RT/14.006700; /*N */
    gms[13] *= RT/29.021370; /*NNH */
    gms[14] *= RT/28.013400; /*N2 */
}


/*Returns helmholtz in mass units (Eq 32.) */
void CKAMS(double *T, int * iwrk, double * rwrk, double * ams)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    helmholtz(ams, tc);
    ams[0] *= RT/1.007970; /*H */
    ams[1] *= RT/15.999400; /*O */
    ams[2] *= RT/17.007370; /*OH */
    ams[3] *= RT/2.015940; /*H2 */
    ams[4] *= RT/31.998800; /*O2 */
    ams[5] *= RT/33.006770; /*HO2 */
    ams[6] *= RT/18.015340; /*H2O */
    ams[7] *= RT/34.014740; /*H2O2 */
    ams[8] *= RT/30.006100; /*NO */
    ams[9] *= RT/46.005500; /*NO2 */
    ams[10] *= RT/44.012800; /*N2O */
    ams[11] *= RT/15.014670; /*NH */
    ams[12] *= RT/14.006700; /*N */
    ams[13] *= RT/29.021370; /*NNH */
    ams[14] *= RT/28.013400; /*N2 */
}


/*Returns the entropies in mass units (Eq 28.) */
void CKSMS(double *T, int * iwrk, double * rwrk, double * sms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 82482613.569848; /*H */
    sms[1] *= 5196444.866683; /*O */
    sms[2] *= 4888468.940230; /*OH */
    sms[3] *= 41241306.784924; /*H2 */
    sms[4] *= 2598222.433341; /*O2 */
    sms[5] *= 2518877.187922; /*HO2 */
    sms[6] *= 4614955.920899; /*H2O */
    sms[7] *= 2444234.470115; /*H2O2 */
    sms[8] *= 2770769.943445; /*NO */
    sms[9] *= 1807175.229049; /*NO2 */
    sms[10] *= 1888995.928457; /*N2O */
    sms[11] *= 5537251.234959; /*NH */
    sms[12] *= 5935730.757423; /*N */
    sms[13] *= 2864785.501167; /*NNH */
    sms[14] *= 2967865.378712; /*N2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double *T, double *x, int * iwrk, double * rwrk, double * cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[15]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 15; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.314e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
void CKCPBS(double *T, double *y, int * iwrk, double * rwrk, double * cpbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[15]; /* temporary storage */
    cp_R(cpor, tc);
    /*multiply by y/molecularweight */
    result += cpor[0]*y[0]/1.00797; /*H */
    result += cpor[1]*y[1]/15.9994; /*O */
    result += cpor[2]*y[2]/17.0074; /*OH */
    result += cpor[3]*y[3]/2.01594; /*H2 */
    result += cpor[4]*y[4]/31.9988; /*O2 */
    result += cpor[5]*y[5]/33.0068; /*HO2 */
    result += cpor[6]*y[6]/18.0153; /*H2O */
    result += cpor[7]*y[7]/34.0147; /*H2O2 */
    result += cpor[8]*y[8]/30.0061; /*NO */
    result += cpor[9]*y[9]/46.0055; /*NO2 */
    result += cpor[10]*y[10]/44.0128; /*N2O */
    result += cpor[11]*y[11]/15.0147; /*NH */
    result += cpor[12]*y[12]/14.0067; /*N */
    result += cpor[13]*y[13]/29.0214; /*NNH */
    result += cpor[14]*y[14]/28.0134; /*N2 */

    *cpbs = result * 8.314e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double *T, double *x, int * iwrk, double * rwrk, double * cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[15]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 15; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.314e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
void CKCVBS(double *T, double *y, int * iwrk, double * rwrk, double * cvbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[15]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]/1.00797; /*H */
    result += cvor[1]*y[1]/15.9994; /*O */
    result += cvor[2]*y[2]/17.0074; /*OH */
    result += cvor[3]*y[3]/2.01594; /*H2 */
    result += cvor[4]*y[4]/31.9988; /*O2 */
    result += cvor[5]*y[5]/33.0068; /*HO2 */
    result += cvor[6]*y[6]/18.0153; /*H2O */
    result += cvor[7]*y[7]/34.0147; /*H2O2 */
    result += cvor[8]*y[8]/30.0061; /*NO */
    result += cvor[9]*y[9]/46.0055; /*NO2 */
    result += cvor[10]*y[10]/44.0128; /*N2O */
    result += cvor[11]*y[11]/15.0147; /*NH */
    result += cvor[12]*y[12]/14.0067; /*N */
    result += cvor[13]*y[13]/29.0214; /*NNH */
    result += cvor[14]*y[14]/28.0134; /*N2 */

    *cvbs = result * 8.314e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double *T, double *x, int * iwrk, double * rwrk, double * hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[15]; /* temporary storage */
    double RT = 8.314e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 15; ++id) {
        result += x[id]*hml[id];
    }

    *hbml = result * RT;
}


/*Returns mean enthalpy of mixture in mass units */
void CKHBMS(double *T, double *y, int * iwrk, double * rwrk, double * hbms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[15]; /* temporary storage */
    double RT = 8.314e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*hml[0]/1.007970; /*H */
    result += y[1]*hml[1]/15.999400; /*O */
    result += y[2]*hml[2]/17.007370; /*OH */
    result += y[3]*hml[3]/2.015940; /*H2 */
    result += y[4]*hml[4]/31.998800; /*O2 */
    result += y[5]*hml[5]/33.006770; /*HO2 */
    result += y[6]*hml[6]/18.015340; /*H2O */
    result += y[7]*hml[7]/34.014740; /*H2O2 */
    result += y[8]*hml[8]/30.006100; /*NO */
    result += y[9]*hml[9]/46.005500; /*NO2 */
    result += y[10]*hml[10]/44.012800; /*N2O */
    result += y[11]*hml[11]/15.014670; /*NH */
    result += y[12]*hml[12]/14.006700; /*N */
    result += y[13]*hml[13]/29.021370; /*NNH */
    result += y[14]*hml[14]/28.013400; /*N2 */

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double *T, double *x, int * iwrk, double * rwrk, double * ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[15]; /* temporary energy array */
    double RT = 8.314e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 15; ++id) {
        result += x[id]*uml[id];
    }

    *ubml = result * RT;
}


/*get mean internal energy in mass units */
void CKUBMS(double *T, double *y, int * iwrk, double * rwrk, double * ubms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double ums[15]; /* temporary energy array */
    double RT = 8.314e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]/1.007970; /*H */
    result += y[1]*ums[1]/15.999400; /*O */
    result += y[2]*ums[2]/17.007370; /*OH */
    result += y[3]*ums[3]/2.015940; /*H2 */
    result += y[4]*ums[4]/31.998800; /*O2 */
    result += y[5]*ums[5]/33.006770; /*HO2 */
    result += y[6]*ums[6]/18.015340; /*H2O */
    result += y[7]*ums[7]/34.014740; /*H2O2 */
    result += y[8]*ums[8]/30.006100; /*NO */
    result += y[9]*ums[9]/46.005500; /*NO2 */
    result += y[10]*ums[10]/44.012800; /*N2O */
    result += y[11]*ums[11]/15.014670; /*NH */
    result += y[12]*ums[12]/14.006700; /*N */
    result += y[13]*ums[13]/29.021370; /*NNH */
    result += y[14]*ums[14]/28.013400; /*N2 */

    *ubms = result * RT;
}


/*get mixture entropy in molar units */
void CKSBML(double *P, double *T, double *x, int * iwrk, double * rwrk, double * sbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[15]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 15; ++id) {
        result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
    }

    *sbml = result * 8.314e+07;
}


/*get mixture entropy in mass units */
void CKSBMS(double *P, double *T, double *y, int * iwrk, double * rwrk, double * sbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[15]; /* temporary storage */
    double x[15]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/1.007970; /*H */
    YOW += y[1]/15.999400; /*O */
    YOW += y[2]/17.007370; /*OH */
    YOW += y[3]/2.015940; /*H2 */
    YOW += y[4]/31.998800; /*O2 */
    YOW += y[5]/33.006770; /*HO2 */
    YOW += y[6]/18.015340; /*H2O */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/30.006100; /*NO */
    YOW += y[9]/46.005500; /*NO2 */
    YOW += y[10]/44.012800; /*N2O */
    YOW += y[11]/15.014670; /*NH */
    YOW += y[12]/14.006700; /*N */
    YOW += y[13]/29.021370; /*NNH */
    YOW += y[14]/28.013400; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(15.999400*YOW); 
    x[2] = y[2]/(17.007370*YOW); 
    x[3] = y[3]/(2.015940*YOW); 
    x[4] = y[4]/(31.998800*YOW); 
    x[5] = y[5]/(33.006770*YOW); 
    x[6] = y[6]/(18.015340*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(30.006100*YOW); 
    x[9] = y[9]/(46.005500*YOW); 
    x[10] = y[10]/(44.012800*YOW); 
    x[11] = y[11]/(15.014670*YOW); 
    x[12] = y[12]/(14.006700*YOW); 
    x[13] = y[13]/(29.021370*YOW); 
    x[14] = y[14]/(28.013400*YOW); 
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
    /*Scale by R/W */
    *sbms = result * 8.314e+07 * YOW;
}


/*Returns mean gibbs free energy in molar units */
void CKGBML(double *P, double *T, double *x, int * iwrk, double * rwrk, double * gbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    double gort[15]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 15; ++id) {
        result += x[id]*(gort[id]+log((x[id]+1e-100))+logPratio);
    }

    *gbml = result * RT;
}


/*Returns mixture gibbs free energy in mass units */
void CKGBMS(double *P, double *T, double *y, int * iwrk, double * rwrk, double * gbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    double gort[15]; /* temporary storage */
    double x[15]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/1.007970; /*H */
    YOW += y[1]/15.999400; /*O */
    YOW += y[2]/17.007370; /*OH */
    YOW += y[3]/2.015940; /*H2 */
    YOW += y[4]/31.998800; /*O2 */
    YOW += y[5]/33.006770; /*HO2 */
    YOW += y[6]/18.015340; /*H2O */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/30.006100; /*NO */
    YOW += y[9]/46.005500; /*NO2 */
    YOW += y[10]/44.012800; /*N2O */
    YOW += y[11]/15.014670; /*NH */
    YOW += y[12]/14.006700; /*N */
    YOW += y[13]/29.021370; /*NNH */
    YOW += y[14]/28.013400; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(15.999400*YOW); 
    x[2] = y[2]/(17.007370*YOW); 
    x[3] = y[3]/(2.015940*YOW); 
    x[4] = y[4]/(31.998800*YOW); 
    x[5] = y[5]/(33.006770*YOW); 
    x[6] = y[6]/(18.015340*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(30.006100*YOW); 
    x[9] = y[9]/(46.005500*YOW); 
    x[10] = y[10]/(44.012800*YOW); 
    x[11] = y[11]/(15.014670*YOW); 
    x[12] = y[12]/(14.006700*YOW); 
    x[13] = y[13]/(29.021370*YOW); 
    x[14] = y[14]/(28.013400*YOW); 
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
    /*Scale by RT/W */
    *gbms = result * RT * YOW;
}


/*Returns mean helmholtz free energy in molar units */
void CKABML(double *P, double *T, double *x, int * iwrk, double * rwrk, double * abml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    double aort[15]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 15; ++id) {
        result += x[id]*(aort[id]+log((x[id]+1e-100))+logPratio);
    }

    *abml = result * RT;
}


/*Returns mixture helmholtz free energy in mass units */
void CKABMS(double *P, double *T, double *y, int * iwrk, double * rwrk, double * abms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    double aort[15]; /* temporary storage */
    double x[15]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/1.007970; /*H */
    YOW += y[1]/15.999400; /*O */
    YOW += y[2]/17.007370; /*OH */
    YOW += y[3]/2.015940; /*H2 */
    YOW += y[4]/31.998800; /*O2 */
    YOW += y[5]/33.006770; /*HO2 */
    YOW += y[6]/18.015340; /*H2O */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/30.006100; /*NO */
    YOW += y[9]/46.005500; /*NO2 */
    YOW += y[10]/44.012800; /*N2O */
    YOW += y[11]/15.014670; /*NH */
    YOW += y[12]/14.006700; /*N */
    YOW += y[13]/29.021370; /*NNH */
    YOW += y[14]/28.013400; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(15.999400*YOW); 
    x[2] = y[2]/(17.007370*YOW); 
    x[3] = y[3]/(2.015940*YOW); 
    x[4] = y[4]/(31.998800*YOW); 
    x[5] = y[5]/(33.006770*YOW); 
    x[6] = y[6]/(18.015340*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(30.006100*YOW); 
    x[9] = y[9]/(46.005500*YOW); 
    x[10] = y[10]/(44.012800*YOW); 
    x[11] = y[11]/(15.014670*YOW); 
    x[12] = y[12]/(14.006700*YOW); 
    x[13] = y[13]/(29.021370*YOW); 
    x[14] = y[14]/(28.013400*YOW); 
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
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double * T, double * C, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 15; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 15; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[15]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/1.007970; /*H */
    YOW += y[1]/15.999400; /*O */
    YOW += y[2]/17.007370; /*OH */
    YOW += y[3]/2.015940; /*H2 */
    YOW += y[4]/31.998800; /*O2 */
    YOW += y[5]/33.006770; /*HO2 */
    YOW += y[6]/18.015340; /*H2O */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/30.006100; /*NO */
    YOW += y[9]/46.005500; /*NO2 */
    YOW += y[10]/44.012800; /*N2O */
    YOW += y[11]/15.014670; /*NH */
    YOW += y[12]/14.006700; /*N */
    YOW += y[13]/29.021370; /*NNH */
    YOW += y[14]/28.013400; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.314e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]/1.007970; 
    c[1] = PWORT * y[1]/15.999400; 
    c[2] = PWORT * y[2]/17.007370; 
    c[3] = PWORT * y[3]/2.015940; 
    c[4] = PWORT * y[4]/31.998800; 
    c[5] = PWORT * y[5]/33.006770; 
    c[6] = PWORT * y[6]/18.015340; 
    c[7] = PWORT * y[7]/34.014740; 
    c[8] = PWORT * y[8]/30.006100; 
    c[9] = PWORT * y[9]/46.005500; 
    c[10] = PWORT * y[10]/44.012800; 
    c[11] = PWORT * y[11]/15.014670; 
    c[12] = PWORT * y[12]/14.006700; 
    c[13] = PWORT * y[13]/29.021370; 
    c[14] = PWORT * y[14]/28.013400; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 15; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[15]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.314e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 15; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 15; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[15]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]/1.007970; 
    c[1] = 1e6 * (*rho) * y[1]/15.999400; 
    c[2] = 1e6 * (*rho) * y[2]/17.007370; 
    c[3] = 1e6 * (*rho) * y[3]/2.015940; 
    c[4] = 1e6 * (*rho) * y[4]/31.998800; 
    c[5] = 1e6 * (*rho) * y[5]/33.006770; 
    c[6] = 1e6 * (*rho) * y[6]/18.015340; 
    c[7] = 1e6 * (*rho) * y[7]/34.014740; 
    c[8] = 1e6 * (*rho) * y[8]/30.006100; 
    c[9] = 1e6 * (*rho) * y[9]/46.005500; 
    c[10] = 1e6 * (*rho) * y[10]/44.012800; 
    c[11] = 1e6 * (*rho) * y[11]/15.014670; 
    c[12] = 1e6 * (*rho) * y[12]/14.006700; 
    c[13] = 1e6 * (*rho) * y[13]/29.021370; 
    c[14] = 1e6 * (*rho) * y[14]/28.013400; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 15; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[15]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*17.007370; /*OH */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*31.998800; /*O2 */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*30.006100; /*NO */
    XW += x[9]*46.005500; /*NO2 */
    XW += x[10]*44.012800; /*N2O */
    XW += x[11]*15.014670; /*NH */
    XW += x[12]*14.006700; /*N */
    XW += x[13]*29.021370; /*NNH */
    XW += x[14]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 15; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 15; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double * T, double * C, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 15; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 15; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 58; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKKFKR(double * P, double * T, double * x, int * iwrk, double * rwrk, double * q_f, double * q_r)
{
    int id; /*loop counter */
    double c[15]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.314e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 15; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    for (id = 0; id < 58; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[15]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/1.007970; /*H */
    YOW += y[1]/15.999400; /*O */
    YOW += y[2]/17.007370; /*OH */
    YOW += y[3]/2.015940; /*H2 */
    YOW += y[4]/31.998800; /*O2 */
    YOW += y[5]/33.006770; /*HO2 */
    YOW += y[6]/18.015340; /*H2O */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/30.006100; /*NO */
    YOW += y[9]/46.005500; /*NO2 */
    YOW += y[10]/44.012800; /*N2O */
    YOW += y[11]/15.014670; /*NH */
    YOW += y[12]/14.006700; /*N */
    YOW += y[13]/29.021370; /*NNH */
    YOW += y[14]/28.013400; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.314e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]/1.007970; 
    c[1] = PWORT * y[1]/15.999400; 
    c[2] = PWORT * y[2]/17.007370; 
    c[3] = PWORT * y[3]/2.015940; 
    c[4] = PWORT * y[4]/31.998800; 
    c[5] = PWORT * y[5]/33.006770; 
    c[6] = PWORT * y[6]/18.015340; 
    c[7] = PWORT * y[7]/34.014740; 
    c[8] = PWORT * y[8]/30.006100; 
    c[9] = PWORT * y[9]/46.005500; 
    c[10] = PWORT * y[10]/44.012800; 
    c[11] = PWORT * y[11]/15.014670; 
    c[12] = PWORT * y[12]/14.006700; 
    c[13] = PWORT * y[13]/29.021370; 
    c[14] = PWORT * y[14]/28.013400; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 58; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[15]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.314e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 15; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 58; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[15]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]/1.007970; 
    c[1] = 1e6 * (*rho) * y[1]/15.999400; 
    c[2] = 1e6 * (*rho) * y[2]/17.007370; 
    c[3] = 1e6 * (*rho) * y[3]/2.015940; 
    c[4] = 1e6 * (*rho) * y[4]/31.998800; 
    c[5] = 1e6 * (*rho) * y[5]/33.006770; 
    c[6] = 1e6 * (*rho) * y[6]/18.015340; 
    c[7] = 1e6 * (*rho) * y[7]/34.014740; 
    c[8] = 1e6 * (*rho) * y[8]/30.006100; 
    c[9] = 1e6 * (*rho) * y[9]/46.005500; 
    c[10] = 1e6 * (*rho) * y[10]/44.012800; 
    c[11] = 1e6 * (*rho) * y[11]/15.014670; 
    c[12] = 1e6 * (*rho) * y[12]/14.006700; 
    c[13] = 1e6 * (*rho) * y[13]/29.021370; 
    c[14] = 1e6 * (*rho) * y[14]/28.013400; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 58; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[15]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*17.007370; /*OH */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*31.998800; /*O2 */
    XW += x[5]*33.006770; /*HO2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*30.006100; /*NO */
    XW += x[9]*46.005500; /*NO2 */
    XW += x[10]*44.012800; /*N2O */
    XW += x[11]*15.014670; /*NH */
    XW += x[12]*14.006700; /*N */
    XW += x[13]*29.021370; /*NNH */
    XW += x[14]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 15; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 58; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the stoichiometric coefficients */
/*of the reaction mechanism. (Eq 50) */
void CKNU(int * kdim, int * iwrk, double * rwrk, int * nuki)
{
    int id; /*loop counter */
    int kd = (*kdim); 
    /*Zero nuki */
    for (id = 0; id < 15 * 58; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: H + O2 <=> O + OH */
    nuki[ 0 * kd + 0 ] = -1 ;
    nuki[ 4 * kd + 0 ] = -1 ;
    nuki[ 1 * kd + 0 ] = +1 ;
    nuki[ 2 * kd + 0 ] = +1 ;

    /*reaction 2: H + H + M <=> H2 + M */
    nuki[ 0 * kd + 1 ] = -1 ;
    nuki[ 0 * kd + 1 ] = -1 ;
    nuki[ 3 * kd + 1 ] = +1 ;

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    nuki[ 0 * kd + 2 ] = -1 ;
    nuki[ 0 * kd + 2 ] = -1 ;
    nuki[ 14 * kd + 2 ] = -1 ;
    nuki[ 3 * kd + 2 ] = +1 ;
    nuki[ 14 * kd + 2 ] = +1 ;

    /*reaction 4: H + H + H2 <=> H2 + H2 */
    nuki[ 0 * kd + 3 ] = -1 ;
    nuki[ 0 * kd + 3 ] = -1 ;
    nuki[ 3 * kd + 3 ] = -1 ;
    nuki[ 3 * kd + 3 ] = +1 ;
    nuki[ 3 * kd + 3 ] = +1 ;

    /*reaction 5: H + H + H2O <=> H2 + H2O */
    nuki[ 0 * kd + 4 ] = -1 ;
    nuki[ 0 * kd + 4 ] = -1 ;
    nuki[ 6 * kd + 4 ] = -1 ;
    nuki[ 3 * kd + 4 ] = +1 ;
    nuki[ 6 * kd + 4 ] = +1 ;

    /*reaction 6: H + O + M <=> OH + M */
    nuki[ 0 * kd + 5 ] = -1 ;
    nuki[ 1 * kd + 5 ] = -1 ;
    nuki[ 2 * kd + 5 ] = +1 ;

    /*reaction 7: H + O2 (+M) <=> HO2 (+M) */
    nuki[ 0 * kd + 6 ] = -1 ;
    nuki[ 4 * kd + 6 ] = -1 ;
    nuki[ 5 * kd + 6 ] = +1 ;

    /*reaction 8: H + O2 (+N2) <=> HO2 (+N2) */
    nuki[ 0 * kd + 7 ] = -1 ;
    nuki[ 4 * kd + 7 ] = -1 ;
    nuki[ 5 * kd + 7 ] = +1 ;

    /*reaction 9: O + O + M <=> O2 + M */
    nuki[ 1 * kd + 8 ] = -1 ;
    nuki[ 1 * kd + 8 ] = -1 ;
    nuki[ 4 * kd + 8 ] = +1 ;

    /*reaction 10: O + H2 <=> OH + H */
    nuki[ 1 * kd + 9 ] = -1 ;
    nuki[ 3 * kd + 9 ] = -1 ;
    nuki[ 2 * kd + 9 ] = +1 ;
    nuki[ 0 * kd + 9 ] = +1 ;

    /*reaction 11: O + H2 <=> OH + H */
    nuki[ 1 * kd + 10 ] = -1 ;
    nuki[ 3 * kd + 10 ] = -1 ;
    nuki[ 2 * kd + 10 ] = +1 ;
    nuki[ 0 * kd + 10 ] = +1 ;

    /*reaction 12: OH + OH <=> O + H2O */
    nuki[ 2 * kd + 11 ] = -1 ;
    nuki[ 2 * kd + 11 ] = -1 ;
    nuki[ 1 * kd + 11 ] = +1 ;
    nuki[ 6 * kd + 11 ] = +1 ;

    /*reaction 13: OH + H + M <=> H2O + M */
    nuki[ 2 * kd + 12 ] = -1 ;
    nuki[ 0 * kd + 12 ] = -1 ;
    nuki[ 6 * kd + 12 ] = +1 ;

    /*reaction 14: OH + H2 <=> H + H2O */
    nuki[ 2 * kd + 13 ] = -1 ;
    nuki[ 3 * kd + 13 ] = -1 ;
    nuki[ 0 * kd + 13 ] = +1 ;
    nuki[ 6 * kd + 13 ] = +1 ;

    /*reaction 15: H2 + O2 <=> HO2 + H */
    nuki[ 3 * kd + 14 ] = -1 ;
    nuki[ 4 * kd + 14 ] = -1 ;
    nuki[ 5 * kd + 14 ] = +1 ;
    nuki[ 0 * kd + 14 ] = +1 ;

    /*reaction 16: HO2 + H <=> OH + OH */
    nuki[ 5 * kd + 15 ] = -1 ;
    nuki[ 0 * kd + 15 ] = -1 ;
    nuki[ 2 * kd + 15 ] = +1 ;
    nuki[ 2 * kd + 15 ] = +1 ;

    /*reaction 17: HO2 + H <=> H2O + O */
    nuki[ 5 * kd + 16 ] = -1 ;
    nuki[ 0 * kd + 16 ] = -1 ;
    nuki[ 6 * kd + 16 ] = +1 ;
    nuki[ 1 * kd + 16 ] = +1 ;

    /*reaction 18: HO2 + O <=> OH + O2 */
    nuki[ 5 * kd + 17 ] = -1 ;
    nuki[ 1 * kd + 17 ] = -1 ;
    nuki[ 2 * kd + 17 ] = +1 ;
    nuki[ 4 * kd + 17 ] = +1 ;

    /*reaction 19: HO2 + OH <=> H2O + O2 */
    nuki[ 5 * kd + 18 ] = -1 ;
    nuki[ 2 * kd + 18 ] = -1 ;
    nuki[ 6 * kd + 18 ] = +1 ;
    nuki[ 4 * kd + 18 ] = +1 ;

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    nuki[ 5 * kd + 19 ] = -1 ;
    nuki[ 2 * kd + 19 ] = -1 ;
    nuki[ 6 * kd + 19 ] = +1 ;
    nuki[ 4 * kd + 19 ] = +1 ;

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    nuki[ 5 * kd + 20 ] = -1 ;
    nuki[ 2 * kd + 20 ] = -1 ;
    nuki[ 6 * kd + 20 ] = +1 ;
    nuki[ 4 * kd + 20 ] = +1 ;

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    nuki[ 5 * kd + 21 ] = -1 ;
    nuki[ 5 * kd + 21 ] = -1 ;
    nuki[ 7 * kd + 21 ] = +1 ;
    nuki[ 4 * kd + 21 ] = +1 ;

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    nuki[ 5 * kd + 22 ] = -1 ;
    nuki[ 5 * kd + 22 ] = -1 ;
    nuki[ 7 * kd + 22 ] = +1 ;
    nuki[ 4 * kd + 22 ] = +1 ;

    /*reaction 24: H2O2 (+M) <=> OH + OH (+M) */
    nuki[ 7 * kd + 23 ] = -1 ;
    nuki[ 2 * kd + 23 ] = +1 ;
    nuki[ 2 * kd + 23 ] = +1 ;

    /*reaction 25: H2O2 + H <=> H2O + OH */
    nuki[ 7 * kd + 24 ] = -1 ;
    nuki[ 0 * kd + 24 ] = -1 ;
    nuki[ 6 * kd + 24 ] = +1 ;
    nuki[ 2 * kd + 24 ] = +1 ;

    /*reaction 26: H2O2 + H <=> HO2 + H2 */
    nuki[ 7 * kd + 25 ] = -1 ;
    nuki[ 0 * kd + 25 ] = -1 ;
    nuki[ 5 * kd + 25 ] = +1 ;
    nuki[ 3 * kd + 25 ] = +1 ;

    /*reaction 27: H2O2 + O <=> HO2 + OH */
    nuki[ 7 * kd + 26 ] = -1 ;
    nuki[ 1 * kd + 26 ] = -1 ;
    nuki[ 5 * kd + 26 ] = +1 ;
    nuki[ 2 * kd + 26 ] = +1 ;

    /*reaction 28: H2O2 + OH <=> H2O + HO2 */
    nuki[ 7 * kd + 27 ] = -1 ;
    nuki[ 2 * kd + 27 ] = -1 ;
    nuki[ 6 * kd + 27 ] = +1 ;
    nuki[ 5 * kd + 27 ] = +1 ;

    /*reaction 29: H2O2 + OH <=> H2O + HO2 */
    nuki[ 7 * kd + 28 ] = -1 ;
    nuki[ 2 * kd + 28 ] = -1 ;
    nuki[ 6 * kd + 28 ] = +1 ;
    nuki[ 5 * kd + 28 ] = +1 ;

    /*reaction 30: NO + O (+M) <=> NO2 (+M) */
    nuki[ 8 * kd + 29 ] = -1 ;
    nuki[ 1 * kd + 29 ] = -1 ;
    nuki[ 9 * kd + 29 ] = +1 ;

    /*reaction 31: NO + HO2 <=> NO2 + OH */
    nuki[ 8 * kd + 30 ] = -1 ;
    nuki[ 5 * kd + 30 ] = -1 ;
    nuki[ 9 * kd + 30 ] = +1 ;
    nuki[ 2 * kd + 30 ] = +1 ;

    /*reaction 32: NO2 + H <=> NO + OH */
    nuki[ 9 * kd + 31 ] = -1 ;
    nuki[ 0 * kd + 31 ] = -1 ;
    nuki[ 8 * kd + 31 ] = +1 ;
    nuki[ 2 * kd + 31 ] = +1 ;

    /*reaction 33: NO2 + O <=> NO + O2 */
    nuki[ 9 * kd + 32 ] = -1 ;
    nuki[ 1 * kd + 32 ] = -1 ;
    nuki[ 8 * kd + 32 ] = +1 ;
    nuki[ 4 * kd + 32 ] = +1 ;

    /*reaction 34: NO2 + NO2 <=> NO + NO + O2 */
    nuki[ 9 * kd + 33 ] = -1 ;
    nuki[ 9 * kd + 33 ] = -1 ;
    nuki[ 8 * kd + 33 ] = +1 ;
    nuki[ 8 * kd + 33 ] = +1 ;
    nuki[ 4 * kd + 33 ] = +1 ;

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    nuki[ 10 * kd + 34 ] = -1 ;
    nuki[ 14 * kd + 34 ] = +1 ;
    nuki[ 1 * kd + 34 ] = +1 ;

    /*reaction 36: N2O + H <=> N2 + OH */
    nuki[ 10 * kd + 35 ] = -1 ;
    nuki[ 0 * kd + 35 ] = -1 ;
    nuki[ 14 * kd + 35 ] = +1 ;
    nuki[ 2 * kd + 35 ] = +1 ;

    /*reaction 37: N2O + H <=> N2 + OH */
    nuki[ 10 * kd + 36 ] = -1 ;
    nuki[ 0 * kd + 36 ] = -1 ;
    nuki[ 14 * kd + 36 ] = +1 ;
    nuki[ 2 * kd + 36 ] = +1 ;

    /*reaction 38: N2O + O <=> NO + NO */
    nuki[ 10 * kd + 37 ] = -1 ;
    nuki[ 1 * kd + 37 ] = -1 ;
    nuki[ 8 * kd + 37 ] = +1 ;
    nuki[ 8 * kd + 37 ] = +1 ;

    /*reaction 39: N2O + O <=> N2 + O2 */
    nuki[ 10 * kd + 38 ] = -1 ;
    nuki[ 1 * kd + 38 ] = -1 ;
    nuki[ 14 * kd + 38 ] = +1 ;
    nuki[ 4 * kd + 38 ] = +1 ;

    /*reaction 40: NH + H <=> N + H2 */
    nuki[ 11 * kd + 39 ] = -1 ;
    nuki[ 0 * kd + 39 ] = -1 ;
    nuki[ 12 * kd + 39 ] = +1 ;
    nuki[ 3 * kd + 39 ] = +1 ;

    /*reaction 41: NH + O <=> NO + H */
    nuki[ 11 * kd + 40 ] = -1 ;
    nuki[ 1 * kd + 40 ] = -1 ;
    nuki[ 8 * kd + 40 ] = +1 ;
    nuki[ 0 * kd + 40 ] = +1 ;

    /*reaction 42: NH + OH <=> N + H2O */
    nuki[ 11 * kd + 41 ] = -1 ;
    nuki[ 2 * kd + 41 ] = -1 ;
    nuki[ 12 * kd + 41 ] = +1 ;
    nuki[ 6 * kd + 41 ] = +1 ;

    /*reaction 43: NH + O2 <=> NO + OH */
    nuki[ 11 * kd + 42 ] = -1 ;
    nuki[ 4 * kd + 42 ] = -1 ;
    nuki[ 8 * kd + 42 ] = +1 ;
    nuki[ 2 * kd + 42 ] = +1 ;

    /*reaction 44: NH + NO <=> N2O + H */
    nuki[ 11 * kd + 43 ] = -1 ;
    nuki[ 8 * kd + 43 ] = -1 ;
    nuki[ 10 * kd + 43 ] = +1 ;
    nuki[ 0 * kd + 43 ] = +1 ;

    /*reaction 45: NH + NO <=> N2O + H */
    nuki[ 11 * kd + 44 ] = -1 ;
    nuki[ 8 * kd + 44 ] = -1 ;
    nuki[ 10 * kd + 44 ] = +1 ;
    nuki[ 0 * kd + 44 ] = +1 ;

    /*reaction 46: NH + NO <=> N2 + OH */
    nuki[ 11 * kd + 45 ] = -1 ;
    nuki[ 8 * kd + 45 ] = -1 ;
    nuki[ 14 * kd + 45 ] = +1 ;
    nuki[ 2 * kd + 45 ] = +1 ;

    /*reaction 47: NH + NO2 <=> N2O + OH */
    nuki[ 11 * kd + 46 ] = -1 ;
    nuki[ 9 * kd + 46 ] = -1 ;
    nuki[ 10 * kd + 46 ] = +1 ;
    nuki[ 2 * kd + 46 ] = +1 ;

    /*reaction 48: N + OH <=> NO + H */
    nuki[ 12 * kd + 47 ] = -1 ;
    nuki[ 2 * kd + 47 ] = -1 ;
    nuki[ 8 * kd + 47 ] = +1 ;
    nuki[ 0 * kd + 47 ] = +1 ;

    /*reaction 49: N + O2 <=> NO + O */
    nuki[ 12 * kd + 48 ] = -1 ;
    nuki[ 4 * kd + 48 ] = -1 ;
    nuki[ 8 * kd + 48 ] = +1 ;
    nuki[ 1 * kd + 48 ] = +1 ;

    /*reaction 50: N + NO <=> N2 + O */
    nuki[ 12 * kd + 49 ] = -1 ;
    nuki[ 8 * kd + 49 ] = -1 ;
    nuki[ 14 * kd + 49 ] = +1 ;
    nuki[ 1 * kd + 49 ] = +1 ;

    /*reaction 51: NNH <=> N2 + H */
    nuki[ 13 * kd + 50 ] = -1 ;
    nuki[ 14 * kd + 50 ] = +1 ;
    nuki[ 0 * kd + 50 ] = +1 ;

    /*reaction 52: NNH + H <=> N2 + H2 */
    nuki[ 13 * kd + 51 ] = -1 ;
    nuki[ 0 * kd + 51 ] = -1 ;
    nuki[ 14 * kd + 51 ] = +1 ;
    nuki[ 3 * kd + 51 ] = +1 ;

    /*reaction 53: NNH + O <=> N2O + H */
    nuki[ 13 * kd + 52 ] = -1 ;
    nuki[ 1 * kd + 52 ] = -1 ;
    nuki[ 10 * kd + 52 ] = +1 ;
    nuki[ 0 * kd + 52 ] = +1 ;

    /*reaction 54: NNH + O <=> N2 + OH */
    nuki[ 13 * kd + 53 ] = -1 ;
    nuki[ 1 * kd + 53 ] = -1 ;
    nuki[ 14 * kd + 53 ] = +1 ;
    nuki[ 2 * kd + 53 ] = +1 ;

    /*reaction 55: NNH + O <=> NH + NO */
    nuki[ 13 * kd + 54 ] = -1 ;
    nuki[ 1 * kd + 54 ] = -1 ;
    nuki[ 11 * kd + 54 ] = +1 ;
    nuki[ 8 * kd + 54 ] = +1 ;

    /*reaction 56: NNH + OH <=> N2 + H2O */
    nuki[ 13 * kd + 55 ] = -1 ;
    nuki[ 2 * kd + 55 ] = -1 ;
    nuki[ 14 * kd + 55 ] = +1 ;
    nuki[ 6 * kd + 55 ] = +1 ;

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    nuki[ 13 * kd + 56 ] = -1 ;
    nuki[ 4 * kd + 56 ] = -1 ;
    nuki[ 14 * kd + 56 ] = +1 ;
    nuki[ 5 * kd + 56 ] = +1 ;

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    nuki[ 13 * kd + 57 ] = -1 ;
    nuki[ 4 * kd + 57 ] = -1 ;
    nuki[ 14 * kd + 57 ] = +1 ;
    nuki[ 0 * kd + 57 ] = +1 ;
    nuki[ 4 * kd + 57 ] = +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < kd * 15; ++ id) {
         ncf[id] = 0; 
    }

    /*H */
    ncf[ 0 * kd + 1 ] = 1; /*H */

    /*O */
    ncf[ 1 * kd + 0 ] = 1; /*O */

    /*OH */
    ncf[ 2 * kd + 0 ] = 1; /*O */
    ncf[ 2 * kd + 1 ] = 1; /*H */

    /*H2 */
    ncf[ 3 * kd + 1 ] = 2; /*H */

    /*O2 */
    ncf[ 4 * kd + 0 ] = 2; /*O */

    /*HO2 */
    ncf[ 5 * kd + 1 ] = 1; /*H */
    ncf[ 5 * kd + 0 ] = 2; /*O */

    /*H2O */
    ncf[ 6 * kd + 1 ] = 2; /*H */
    ncf[ 6 * kd + 0 ] = 1; /*O */

    /*H2O2 */
    ncf[ 7 * kd + 1 ] = 2; /*H */
    ncf[ 7 * kd + 0 ] = 2; /*O */

    /*NO */
    ncf[ 8 * kd + 2 ] = 1; /*N */
    ncf[ 8 * kd + 0 ] = 1; /*O */

    /*NO2 */
    ncf[ 9 * kd + 2 ] = 1; /*N */
    ncf[ 9 * kd + 0 ] = 2; /*O */

    /*N2O */
    ncf[ 10 * kd + 2 ] = 2; /*N */
    ncf[ 10 * kd + 0 ] = 1; /*O */

    /*NH */
    ncf[ 11 * kd + 2 ] = 1; /*N */
    ncf[ 11 * kd + 1 ] = 1; /*H */

    /*N */
    ncf[ 12 * kd + 2 ] = 1; /*N */

    /*NNH */
    ncf[ 13 * kd + 2 ] = 2; /*N */
    ncf[ 13 * kd + 1 ] = 1; /*H */

    /*N2 */
    ncf[ 14 * kd + 2 ] = 2; /*N */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE(int * iwrk, double * rwrk, double * a, double * b, double * e)
{

    /*reaction 1: H + O2 <=> O + OH */
    a[0] = 3.6e+15;
    b[0] = -0.41;
    e[0] = 16600;

    /*reaction 2: H + H + M <=> H2 + M */
    a[1] = 7e+17;
    b[1] = -1;
    e[1] = 0;

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    a[2] = 5.4e+18;
    b[2] = -1.3;
    e[2] = 0;

    /*reaction 4: H + H + H2 <=> H2 + H2 */
    a[3] = 1e+17;
    b[3] = -0.6;
    e[3] = 0;

    /*reaction 5: H + H + H2O <=> H2 + H2O */
    a[4] = 1e+19;
    b[4] = -1;
    e[4] = 0;

    /*reaction 6: H + O + M <=> OH + M */
    a[5] = 6.2e+16;
    b[5] = -0.6;
    e[5] = 0;

    /*reaction 7: H + O2 (+M) <=> HO2 (+M) */
    a[6] = 1.5e+12;
    b[6] = 0.6;
    e[6] = 0;

    /*reaction 8: H + O2 (+N2) <=> HO2 (+N2) */
    a[7] = 1.5e+12;
    b[7] = 0.6;
    e[7] = 0;

    /*reaction 9: O + O + M <=> O2 + M */
    a[8] = 1.9e+13;
    b[8] = 0;
    e[8] = -1788;

    /*reaction 10: O + H2 <=> OH + H */
    a[9] = 3.8e+12;
    b[9] = 0;
    e[9] = 7948;

    /*reaction 11: O + H2 <=> OH + H */
    a[10] = 8.8e+14;
    b[10] = 0;
    e[10] = 19175;

    /*reaction 12: OH + OH <=> O + H2O */
    a[11] = 4300;
    b[11] = 2.7;
    e[11] = -1822;

    /*reaction 13: OH + H + M <=> H2O + M */
    a[12] = 4.5e+22;
    b[12] = -2;
    e[12] = 0;

    /*reaction 14: OH + H2 <=> H + H2O */
    a[13] = 2.1e+08;
    b[13] = 1.52;
    e[13] = 3449;

    /*reaction 15: H2 + O2 <=> HO2 + H */
    a[14] = 740000;
    b[14] = 2.433;
    e[14] = 53502;

    /*reaction 16: HO2 + H <=> OH + OH */
    a[15] = 8.4e+13;
    b[15] = 0;
    e[15] = 400;

    /*reaction 17: HO2 + H <=> H2O + O */
    a[16] = 1.4e+12;
    b[16] = 0;
    e[16] = 0;

    /*reaction 18: HO2 + O <=> OH + O2 */
    a[17] = 1.6e+13;
    b[17] = 0;
    e[17] = -445;

    /*reaction 19: HO2 + OH <=> H2O + O2 */
    a[18] = 3.6e+21;
    b[18] = -2.1;
    e[18] = 9000;

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    a[19] = 2e+15;
    b[19] = -0.6;
    e[19] = 0;

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    a[20] = -2.2e+96;
    b[20] = -24;
    e[20] = 49000;

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    a[21] = 1.9e+11;
    b[21] = 0;
    e[21] = -1408;

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    a[22] = 1e+14;
    b[22] = 0;
    e[22] = 11034;

    /*reaction 24: H2O2 (+M) <=> OH + OH (+M) */
    a[23] = 4e+11;
    b[23] = 0;
    e[23] = 37137;

    /*reaction 25: H2O2 + H <=> H2O + OH */
    a[24] = 1e+13;
    b[24] = 0;
    e[24] = 3580;

    /*reaction 26: H2O2 + H <=> HO2 + H2 */
    a[25] = 1.7e+12;
    b[25] = 0;
    e[25] = 3760;

    /*reaction 27: H2O2 + O <=> HO2 + OH */
    a[26] = 9.6e+06;
    b[26] = 2;
    e[26] = 3970;

    /*reaction 28: H2O2 + OH <=> H2O + HO2 */
    a[27] = 1.9e+12;
    b[27] = 0;
    e[27] = 427;

    /*reaction 29: H2O2 + OH <=> H2O + HO2 */
    a[28] = 1.6e+18;
    b[28] = 0;
    e[28] = 29410;

    /*reaction 30: NO + O (+M) <=> NO2 (+M) */
    a[29] = 1.3e+15;
    b[29] = -0.75;
    e[29] = 0;

    /*reaction 31: NO + HO2 <=> NO2 + OH */
    a[30] = 2.1e+12;
    b[30] = 0;
    e[30] = -497;

    /*reaction 32: NO2 + H <=> NO + OH */
    a[31] = 1.3e+14;
    b[31] = 0;
    e[31] = 362;

    /*reaction 33: NO2 + O <=> NO + O2 */
    a[32] = 1.1e+14;
    b[32] = -0.52;
    e[32] = 0;

    /*reaction 34: NO2 + NO2 <=> NO + NO + O2 */
    a[33] = 4.5e+12;
    b[33] = 0;
    e[33] = 27599;

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    a[34] = 1.3e+12;
    b[34] = 0;
    e[34] = 62570;

    /*reaction 36: N2O + H <=> N2 + OH */
    a[35] = 3.3e+10;
    b[35] = 0;
    e[35] = 4729;

    /*reaction 37: N2O + H <=> N2 + OH */
    a[36] = 4.4e+14;
    b[36] = 0;
    e[36] = 19254;

    /*reaction 38: N2O + O <=> NO + NO */
    a[37] = 9.2e+13;
    b[37] = 0;
    e[37] = 27679;

    /*reaction 39: N2O + O <=> N2 + O2 */
    a[38] = 3.7e+12;
    b[38] = 0;
    e[38] = 15936;

    /*reaction 40: NH + H <=> N + H2 */
    a[39] = 3e+13;
    b[39] = 0;
    e[39] = 0;

    /*reaction 41: NH + O <=> NO + H */
    a[40] = 9.2e+13;
    b[40] = 0;
    e[40] = 0;

    /*reaction 42: NH + OH <=> N + H2O */
    a[41] = 5e+11;
    b[41] = 0.5;
    e[41] = 2000;

    /*reaction 43: NH + O2 <=> NO + OH */
    a[42] = 1.3e+06;
    b[42] = 1.5;
    e[42] = 100;

    /*reaction 44: NH + NO <=> N2O + H */
    a[43] = 2.9e+14;
    b[43] = -0.4;
    e[43] = 0;

    /*reaction 45: NH + NO <=> N2O + H */
    a[44] = -2.2e+13;
    b[44] = -0.23;
    e[44] = 0;

    /*reaction 46: NH + NO <=> N2 + OH */
    a[45] = 2.2e+13;
    b[45] = -0.23;
    e[45] = 0;

    /*reaction 47: NH + NO2 <=> N2O + OH */
    a[46] = 1e+13;
    b[46] = 0;
    e[46] = 0;

    /*reaction 48: N + OH <=> NO + H */
    a[47] = 3.8e+13;
    b[47] = 0;
    e[47] = 0;

    /*reaction 49: N + O2 <=> NO + O */
    a[48] = 6.4e+09;
    b[48] = 1;
    e[48] = 6280;

    /*reaction 50: N + NO <=> N2 + O */
    a[49] = 2.1e+13;
    b[49] = 0;
    e[49] = 0;

    /*reaction 51: NNH <=> N2 + H */
    a[50] = 6.5e+07;
    b[50] = 0;
    e[50] = 0;

    /*reaction 52: NNH + H <=> N2 + H2 */
    a[51] = 1e+14;
    b[51] = 0;
    e[51] = 0;

    /*reaction 53: NNH + O <=> N2O + H */
    a[52] = 1e+14;
    b[52] = 0;
    e[52] = 0;

    /*reaction 54: NNH + O <=> N2 + OH */
    a[53] = 8e+13;
    b[53] = 0;
    e[53] = 0;

    /*reaction 55: NNH + O <=> NH + NO */
    a[54] = 5e+13;
    b[54] = 0;
    e[54] = 0;

    /*reaction 56: NNH + OH <=> N2 + H2O */
    a[55] = 5e+13;
    b[55] = 0;
    e[55] = 0;

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    a[56] = 2e+14;
    b[56] = 0;
    e[56] = 0;

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    a[57] = 5e+13;
    b[57] = 0;
    e[57] = 0;

    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double * T, double * C, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[15]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> O + OH */
    /*eqcon[0] *= 1;  */

    /*reaction 2: H + H + M <=> H2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + H2 <=> H2 + H2 */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + H + H2O <=> H2 + H2O */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + O + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + O2 (+M) <=> HO2 (+M) */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 (+N2) <=> HO2 (+N2) */
    eqcon[7] *= 1e+06; 

    /*reaction 9: O + O + M <=> O2 + M */
    eqcon[8] *= 1e+06; 

    /*reaction 10: O + H2 <=> OH + H */
    /*eqcon[9] *= 1;  */

    /*reaction 11: O + H2 <=> OH + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: OH + H + M <=> H2O + M */
    eqcon[12] *= 1e+06; 

    /*reaction 14: OH + H2 <=> H + H2O */
    /*eqcon[13] *= 1;  */

    /*reaction 15: H2 + O2 <=> HO2 + H */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + H <=> OH + OH */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + H <=> H2O + O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + O <=> OH + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + OH <=> H2O + O2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[23] *= 1e-06; 

    /*reaction 25: H2O2 + H <=> H2O + OH */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + H <=> HO2 + H2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + O <=> HO2 + OH */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[28] *= 1;  */

    /*reaction 30: NO + O (+M) <=> NO2 (+M) */
    eqcon[29] *= 1e+06; 

    /*reaction 31: NO + HO2 <=> NO2 + OH */
    /*eqcon[30] *= 1;  */

    /*reaction 32: NO2 + H <=> NO + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: NO2 + O <=> NO + O2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: NO2 + NO2 <=> NO + NO + O2 */
    eqcon[33] *= 1e-06; 

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    eqcon[34] *= 1e-06; 

    /*reaction 36: N2O + H <=> N2 + OH */
    /*eqcon[35] *= 1;  */

    /*reaction 37: N2O + H <=> N2 + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: N2O + O <=> NO + NO */
    /*eqcon[37] *= 1;  */

    /*reaction 39: N2O + O <=> N2 + O2 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: NH + H <=> N + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: NH + O <=> NO + H */
    /*eqcon[40] *= 1;  */

    /*reaction 42: NH + OH <=> N + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: NH + O2 <=> NO + OH */
    /*eqcon[42] *= 1;  */

    /*reaction 44: NH + NO <=> N2O + H */
    /*eqcon[43] *= 1;  */

    /*reaction 45: NH + NO <=> N2O + H */
    /*eqcon[44] *= 1;  */

    /*reaction 46: NH + NO <=> N2 + OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: NH + NO2 <=> N2O + OH */
    /*eqcon[46] *= 1;  */

    /*reaction 48: N + OH <=> NO + H */
    /*eqcon[47] *= 1;  */

    /*reaction 49: N + O2 <=> NO + O */
    /*eqcon[48] *= 1;  */

    /*reaction 50: N + NO <=> N2 + O */
    /*eqcon[49] *= 1;  */

    /*reaction 51: NNH <=> N2 + H */
    eqcon[50] *= 1e-06; 

    /*reaction 52: NNH + H <=> N2 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: NNH + O <=> N2O + H */
    /*eqcon[52] *= 1;  */

    /*reaction 54: NNH + O <=> N2 + OH */
    /*eqcon[53] *= 1;  */

    /*reaction 55: NNH + O <=> NH + NO */
    /*eqcon[54] *= 1;  */

    /*reaction 56: NNH + OH <=> N2 + H2O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    eqcon[57] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[15]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> O + OH */
    /*eqcon[0] *= 1;  */

    /*reaction 2: H + H + M <=> H2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + H2 <=> H2 + H2 */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + H + H2O <=> H2 + H2O */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + O + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + O2 (+M) <=> HO2 (+M) */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 (+N2) <=> HO2 (+N2) */
    eqcon[7] *= 1e+06; 

    /*reaction 9: O + O + M <=> O2 + M */
    eqcon[8] *= 1e+06; 

    /*reaction 10: O + H2 <=> OH + H */
    /*eqcon[9] *= 1;  */

    /*reaction 11: O + H2 <=> OH + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: OH + H + M <=> H2O + M */
    eqcon[12] *= 1e+06; 

    /*reaction 14: OH + H2 <=> H + H2O */
    /*eqcon[13] *= 1;  */

    /*reaction 15: H2 + O2 <=> HO2 + H */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + H <=> OH + OH */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + H <=> H2O + O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + O <=> OH + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + OH <=> H2O + O2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[23] *= 1e-06; 

    /*reaction 25: H2O2 + H <=> H2O + OH */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + H <=> HO2 + H2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + O <=> HO2 + OH */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[28] *= 1;  */

    /*reaction 30: NO + O (+M) <=> NO2 (+M) */
    eqcon[29] *= 1e+06; 

    /*reaction 31: NO + HO2 <=> NO2 + OH */
    /*eqcon[30] *= 1;  */

    /*reaction 32: NO2 + H <=> NO + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: NO2 + O <=> NO + O2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: NO2 + NO2 <=> NO + NO + O2 */
    eqcon[33] *= 1e-06; 

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    eqcon[34] *= 1e-06; 

    /*reaction 36: N2O + H <=> N2 + OH */
    /*eqcon[35] *= 1;  */

    /*reaction 37: N2O + H <=> N2 + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: N2O + O <=> NO + NO */
    /*eqcon[37] *= 1;  */

    /*reaction 39: N2O + O <=> N2 + O2 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: NH + H <=> N + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: NH + O <=> NO + H */
    /*eqcon[40] *= 1;  */

    /*reaction 42: NH + OH <=> N + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: NH + O2 <=> NO + OH */
    /*eqcon[42] *= 1;  */

    /*reaction 44: NH + NO <=> N2O + H */
    /*eqcon[43] *= 1;  */

    /*reaction 45: NH + NO <=> N2O + H */
    /*eqcon[44] *= 1;  */

    /*reaction 46: NH + NO <=> N2 + OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: NH + NO2 <=> N2O + OH */
    /*eqcon[46] *= 1;  */

    /*reaction 48: N + OH <=> NO + H */
    /*eqcon[47] *= 1;  */

    /*reaction 49: N + O2 <=> NO + O */
    /*eqcon[48] *= 1;  */

    /*reaction 50: N + NO <=> N2 + O */
    /*eqcon[49] *= 1;  */

    /*reaction 51: NNH <=> N2 + H */
    eqcon[50] *= 1e-06; 

    /*reaction 52: NNH + H <=> N2 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: NNH + O <=> N2O + H */
    /*eqcon[52] *= 1;  */

    /*reaction 54: NNH + O <=> N2 + OH */
    /*eqcon[53] *= 1;  */

    /*reaction 55: NNH + O <=> NH + NO */
    /*eqcon[54] *= 1;  */

    /*reaction 56: NNH + OH <=> N2 + H2O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    eqcon[57] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[15]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> O + OH */
    /*eqcon[0] *= 1;  */

    /*reaction 2: H + H + M <=> H2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + H2 <=> H2 + H2 */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + H + H2O <=> H2 + H2O */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + O + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + O2 (+M) <=> HO2 (+M) */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 (+N2) <=> HO2 (+N2) */
    eqcon[7] *= 1e+06; 

    /*reaction 9: O + O + M <=> O2 + M */
    eqcon[8] *= 1e+06; 

    /*reaction 10: O + H2 <=> OH + H */
    /*eqcon[9] *= 1;  */

    /*reaction 11: O + H2 <=> OH + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: OH + H + M <=> H2O + M */
    eqcon[12] *= 1e+06; 

    /*reaction 14: OH + H2 <=> H + H2O */
    /*eqcon[13] *= 1;  */

    /*reaction 15: H2 + O2 <=> HO2 + H */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + H <=> OH + OH */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + H <=> H2O + O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + O <=> OH + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + OH <=> H2O + O2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[23] *= 1e-06; 

    /*reaction 25: H2O2 + H <=> H2O + OH */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + H <=> HO2 + H2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + O <=> HO2 + OH */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[28] *= 1;  */

    /*reaction 30: NO + O (+M) <=> NO2 (+M) */
    eqcon[29] *= 1e+06; 

    /*reaction 31: NO + HO2 <=> NO2 + OH */
    /*eqcon[30] *= 1;  */

    /*reaction 32: NO2 + H <=> NO + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: NO2 + O <=> NO + O2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: NO2 + NO2 <=> NO + NO + O2 */
    eqcon[33] *= 1e-06; 

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    eqcon[34] *= 1e-06; 

    /*reaction 36: N2O + H <=> N2 + OH */
    /*eqcon[35] *= 1;  */

    /*reaction 37: N2O + H <=> N2 + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: N2O + O <=> NO + NO */
    /*eqcon[37] *= 1;  */

    /*reaction 39: N2O + O <=> N2 + O2 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: NH + H <=> N + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: NH + O <=> NO + H */
    /*eqcon[40] *= 1;  */

    /*reaction 42: NH + OH <=> N + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: NH + O2 <=> NO + OH */
    /*eqcon[42] *= 1;  */

    /*reaction 44: NH + NO <=> N2O + H */
    /*eqcon[43] *= 1;  */

    /*reaction 45: NH + NO <=> N2O + H */
    /*eqcon[44] *= 1;  */

    /*reaction 46: NH + NO <=> N2 + OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: NH + NO2 <=> N2O + OH */
    /*eqcon[46] *= 1;  */

    /*reaction 48: N + OH <=> NO + H */
    /*eqcon[47] *= 1;  */

    /*reaction 49: N + O2 <=> NO + O */
    /*eqcon[48] *= 1;  */

    /*reaction 50: N + NO <=> N2 + O */
    /*eqcon[49] *= 1;  */

    /*reaction 51: NNH <=> N2 + H */
    eqcon[50] *= 1e-06; 

    /*reaction 52: NNH + H <=> N2 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: NNH + O <=> N2O + H */
    /*eqcon[52] *= 1;  */

    /*reaction 54: NNH + O <=> N2 + OH */
    /*eqcon[53] *= 1;  */

    /*reaction 55: NNH + O <=> NH + NO */
    /*eqcon[54] *= 1;  */

    /*reaction 56: NNH + OH <=> N2 + H2O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    eqcon[57] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[15]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> O + OH */
    /*eqcon[0] *= 1;  */

    /*reaction 2: H + H + M <=> H2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + H2 <=> H2 + H2 */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + H + H2O <=> H2 + H2O */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + O + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + O2 (+M) <=> HO2 (+M) */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 (+N2) <=> HO2 (+N2) */
    eqcon[7] *= 1e+06; 

    /*reaction 9: O + O + M <=> O2 + M */
    eqcon[8] *= 1e+06; 

    /*reaction 10: O + H2 <=> OH + H */
    /*eqcon[9] *= 1;  */

    /*reaction 11: O + H2 <=> OH + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: OH + H + M <=> H2O + M */
    eqcon[12] *= 1e+06; 

    /*reaction 14: OH + H2 <=> H + H2O */
    /*eqcon[13] *= 1;  */

    /*reaction 15: H2 + O2 <=> HO2 + H */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + H <=> OH + OH */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + H <=> H2O + O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + O <=> OH + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + OH <=> H2O + O2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[23] *= 1e-06; 

    /*reaction 25: H2O2 + H <=> H2O + OH */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + H <=> HO2 + H2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + O <=> HO2 + OH */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[28] *= 1;  */

    /*reaction 30: NO + O (+M) <=> NO2 (+M) */
    eqcon[29] *= 1e+06; 

    /*reaction 31: NO + HO2 <=> NO2 + OH */
    /*eqcon[30] *= 1;  */

    /*reaction 32: NO2 + H <=> NO + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: NO2 + O <=> NO + O2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: NO2 + NO2 <=> NO + NO + O2 */
    eqcon[33] *= 1e-06; 

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    eqcon[34] *= 1e-06; 

    /*reaction 36: N2O + H <=> N2 + OH */
    /*eqcon[35] *= 1;  */

    /*reaction 37: N2O + H <=> N2 + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: N2O + O <=> NO + NO */
    /*eqcon[37] *= 1;  */

    /*reaction 39: N2O + O <=> N2 + O2 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: NH + H <=> N + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: NH + O <=> NO + H */
    /*eqcon[40] *= 1;  */

    /*reaction 42: NH + OH <=> N + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: NH + O2 <=> NO + OH */
    /*eqcon[42] *= 1;  */

    /*reaction 44: NH + NO <=> N2O + H */
    /*eqcon[43] *= 1;  */

    /*reaction 45: NH + NO <=> N2O + H */
    /*eqcon[44] *= 1;  */

    /*reaction 46: NH + NO <=> N2 + OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: NH + NO2 <=> N2O + OH */
    /*eqcon[46] *= 1;  */

    /*reaction 48: N + OH <=> NO + H */
    /*eqcon[47] *= 1;  */

    /*reaction 49: N + O2 <=> NO + O */
    /*eqcon[48] *= 1;  */

    /*reaction 50: N + NO <=> N2 + O */
    /*eqcon[49] *= 1;  */

    /*reaction 51: NNH <=> N2 + H */
    eqcon[50] *= 1e-06; 

    /*reaction 52: NNH + H <=> N2 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: NNH + O <=> N2O + H */
    /*eqcon[52] *= 1;  */

    /*reaction 54: NNH + O <=> N2 + OH */
    /*eqcon[53] *= 1;  */

    /*reaction 55: NNH + O <=> NH + NO */
    /*eqcon[54] *= 1;  */

    /*reaction 56: NNH + OH <=> N2 + H2O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    eqcon[57] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[15]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> O + OH */
    /*eqcon[0] *= 1;  */

    /*reaction 2: H + H + M <=> H2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + H2 <=> H2 + H2 */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + H + H2O <=> H2 + H2O */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + O + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + O2 (+M) <=> HO2 (+M) */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 (+N2) <=> HO2 (+N2) */
    eqcon[7] *= 1e+06; 

    /*reaction 9: O + O + M <=> O2 + M */
    eqcon[8] *= 1e+06; 

    /*reaction 10: O + H2 <=> OH + H */
    /*eqcon[9] *= 1;  */

    /*reaction 11: O + H2 <=> OH + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: OH + H + M <=> H2O + M */
    eqcon[12] *= 1e+06; 

    /*reaction 14: OH + H2 <=> H + H2O */
    /*eqcon[13] *= 1;  */

    /*reaction 15: H2 + O2 <=> HO2 + H */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + H <=> OH + OH */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + H <=> H2O + O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + O <=> OH + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + OH <=> H2O + O2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[23] *= 1e-06; 

    /*reaction 25: H2O2 + H <=> H2O + OH */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + H <=> HO2 + H2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + O <=> HO2 + OH */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H2O2 + OH <=> H2O + HO2 */
    /*eqcon[28] *= 1;  */

    /*reaction 30: NO + O (+M) <=> NO2 (+M) */
    eqcon[29] *= 1e+06; 

    /*reaction 31: NO + HO2 <=> NO2 + OH */
    /*eqcon[30] *= 1;  */

    /*reaction 32: NO2 + H <=> NO + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: NO2 + O <=> NO + O2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: NO2 + NO2 <=> NO + NO + O2 */
    eqcon[33] *= 1e-06; 

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    eqcon[34] *= 1e-06; 

    /*reaction 36: N2O + H <=> N2 + OH */
    /*eqcon[35] *= 1;  */

    /*reaction 37: N2O + H <=> N2 + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: N2O + O <=> NO + NO */
    /*eqcon[37] *= 1;  */

    /*reaction 39: N2O + O <=> N2 + O2 */
    /*eqcon[38] *= 1;  */

    /*reaction 40: NH + H <=> N + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: NH + O <=> NO + H */
    /*eqcon[40] *= 1;  */

    /*reaction 42: NH + OH <=> N + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: NH + O2 <=> NO + OH */
    /*eqcon[42] *= 1;  */

    /*reaction 44: NH + NO <=> N2O + H */
    /*eqcon[43] *= 1;  */

    /*reaction 45: NH + NO <=> N2O + H */
    /*eqcon[44] *= 1;  */

    /*reaction 46: NH + NO <=> N2 + OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: NH + NO2 <=> N2O + OH */
    /*eqcon[46] *= 1;  */

    /*reaction 48: N + OH <=> NO + H */
    /*eqcon[47] *= 1;  */

    /*reaction 49: N + O2 <=> NO + O */
    /*eqcon[48] *= 1;  */

    /*reaction 50: N + NO <=> N2 + O */
    /*eqcon[49] *= 1;  */

    /*reaction 51: NNH <=> N2 + H */
    eqcon[50] *= 1e-06; 

    /*reaction 52: NNH + H <=> N2 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: NNH + O <=> N2O + H */
    /*eqcon[52] *= 1;  */

    /*reaction 54: NNH + O <=> N2 + OH */
    /*eqcon[53] *= 1;  */

    /*reaction 55: NNH + O <=> NH + NO */
    /*eqcon[54] *= 1;  */

    /*reaction 56: NNH + OH <=> N2 + H2O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    eqcon[57] *= 1e-06; 
}

/*compute the production rate for each species */
void productionRate(double * wdot, double * sc, double T)
{
    double qdot;

    static double T_old = -1, k_f_old[58], Kc_old[58];

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[15];                /*Gibbs free energy */
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

    double X;                       /*SRI intermediate */
    double F_sri;                   /*SRI intermediate */
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.314 / T;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 15; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*zero out wdot */
    for (id = 0; id < 15; ++id) {
        wdot[id] = 0.0;
    }

    if (T != T_old)
    {
        T_old = T;

        k_f_old[0] = 1e-06 * 3.6e+15*exp(-0.41*tc[0]-8354.3/tc[1]);
        k_f_old[1] = 1e-12 * 7e+17*exp(-1*tc[0]);
        k_f_old[2] = 1e-12 * 5.4e+18*exp(-1.3*tc[0]);
        k_f_old[3] = 1e-12 * 1e+17*exp(-0.6*tc[0]);
        k_f_old[4] = 1e-12 * 1e+19*exp(-1*tc[0]);
        k_f_old[5] = 1e-12 * 6.2e+16*exp(-0.6*tc[0]);
        k_f_old[6] = 1e-06 * 1.5e+12*exp(0.6*tc[0]);
        k_f_old[7] = 1e-06 * 1.5e+12*exp(0.6*tc[0]);
        k_f_old[8] = 1e-12 * 1.9e+13*exp(+899.849/tc[1]);
        k_f_old[9] = 1e-06 * 3.8e+12*exp(-4000/tc[1]);
        k_f_old[10] = 1e-06 * 8.8e+14*exp(-9650.23/tc[1]);
        k_f_old[11] = 1e-06 * 4300*exp(2.7*tc[0]+916.96/tc[1]);
        k_f_old[12] = 1e-12 * 4.5e+22*exp(-2*tc[0]);
        k_f_old[13] = 1e-06 * 2.1e+08*exp(1.52*tc[0]-1735.78/tc[1]);
        k_f_old[14] = 1e-06 * 740000*exp(2.433*tc[0]-26926/tc[1]);
        k_f_old[15] = 1e-06 * 8.4e+13*exp(-201.309/tc[1]);
        k_f_old[16] = 1e-06 * 1.4e+12;
        k_f_old[17] = 1e-06 * 1.6e+13*exp(+223.956/tc[1]);
        k_f_old[18] = 1e-06 * 3.6e+21*exp(-2.1*tc[0]-4529.44/tc[1]);
        k_f_old[19] = 1e-06 * 2e+15*exp(-0.6*tc[0]);
        k_f_old[20] = 1e-06 * -2.2e+96*exp(-24*tc[0]-24660.3/tc[1]);
        k_f_old[21] = 1e-06 * 1.9e+11*exp(+708.606/tc[1]);
        k_f_old[22] = 1e-06 * 1e+14*exp(-5553.1/tc[1]);
        k_f_old[23] = 1 * 4e+11*exp(-18690/tc[1]);
        k_f_old[24] = 1e-06 * 1e+13*exp(-1801.71/tc[1]);
        k_f_old[25] = 1e-06 * 1.7e+12*exp(-1892.3/tc[1]);
        k_f_old[26] = 1e-06 * 9.6e+06*exp(2*tc[0]-1997.99/tc[1]);
        k_f_old[27] = 1e-06 * 1.9e+12*exp(-214.897/tc[1]);
        k_f_old[28] = 1e-06 * 1.6e+18*exp(-14801.2/tc[1]);
        k_f_old[29] = 1e-06 * 1.3e+15*exp(-0.75*tc[0]);
        k_f_old[30] = 1e-06 * 2.1e+12*exp(+250.126/tc[1]);
        k_f_old[31] = 1e-06 * 1.3e+14*exp(-182.184/tc[1]);
        k_f_old[32] = 1e-06 * 1.1e+14*exp(-0.52*tc[0]);
        k_f_old[33] = 1e-06 * 4.5e+12*exp(-13889.8/tc[1]);
        k_f_old[34] = 1 * 1.3e+12*exp(-31489.7/tc[1]);
        k_f_old[35] = 1e-06 * 3.3e+10*exp(-2379.97/tc[1]);
        k_f_old[36] = 1e-06 * 4.4e+14*exp(-9689.98/tc[1]);
        k_f_old[37] = 1e-06 * 9.2e+13*exp(-13930/tc[1]);
        k_f_old[38] = 1e-06 * 3.7e+12*exp(-8020.13/tc[1]);
        k_f_old[39] = 1e-06 * 3e+13;
        k_f_old[40] = 1e-06 * 9.2e+13;
        k_f_old[41] = 1e-06 * 5e+11*exp(0.5*tc[0]-1006.54/tc[1]);
        k_f_old[42] = 1e-06 * 1.3e+06*exp(1.5*tc[0]-50.3271/tc[1]);
        k_f_old[43] = 1e-06 * 2.9e+14*exp(-0.4*tc[0]);
        k_f_old[44] = 1e-06 * -2.2e+13*exp(-0.23*tc[0]);
        k_f_old[45] = 1e-06 * 2.2e+13*exp(-0.23*tc[0]);
        k_f_old[46] = 1e-06 * 1e+13;
        k_f_old[47] = 1e-06 * 3.8e+13;
        k_f_old[48] = 1e-06 * 6.4e+09*exp(1*tc[0]-3160.54/tc[1]);
        k_f_old[49] = 1e-06 * 2.1e+13;
        k_f_old[50] = 1 * 6.5e+07;
        k_f_old[51] = 1e-06 * 1e+14;
        k_f_old[52] = 1e-06 * 1e+14;
        k_f_old[53] = 1e-06 * 8e+13;
        k_f_old[54] = 1e-06 * 5e+13;
        k_f_old[55] = 1e-06 * 5e+13;
        k_f_old[56] = 1e-06 * 2e+14;
        k_f_old[57] = 1e-06 * 5e+13;

        Kc_old[0] = exp((g_RT[0] + g_RT[4]) - (g_RT[1] + g_RT[2]));
        Kc_old[1] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0]) - (g_RT[3]));
        Kc_old[2] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[14]) - (g_RT[3] + g_RT[14]));
        Kc_old[3] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[3]) - (g_RT[3] + g_RT[3]));
        Kc_old[4] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[6]) - (g_RT[3] + g_RT[6]));
        Kc_old[5] = 1.0 / (refC) * exp((g_RT[0] + g_RT[1]) - (g_RT[2]));
        Kc_old[6] = 1.0 / (refC) * exp((g_RT[0] + g_RT[4]) - (g_RT[5]));
        Kc_old[7] = 1.0 / (refC) * exp((g_RT[0] + g_RT[4]) - (g_RT[5]));
        Kc_old[8] = 1.0 / (refC) * exp((g_RT[1] + g_RT[1]) - (g_RT[4]));
        Kc_old[9] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[0]));
        Kc_old[10] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[0]));
        Kc_old[11] = exp((g_RT[2] + g_RT[2]) - (g_RT[1] + g_RT[6]));
        Kc_old[12] = 1.0 / (refC) * exp((g_RT[2] + g_RT[0]) - (g_RT[6]));
        Kc_old[13] = exp((g_RT[2] + g_RT[3]) - (g_RT[0] + g_RT[6]));
        Kc_old[14] = exp((g_RT[3] + g_RT[4]) - (g_RT[5] + g_RT[0]));
        Kc_old[15] = exp((g_RT[5] + g_RT[0]) - (g_RT[2] + g_RT[2]));
        Kc_old[16] = exp((g_RT[5] + g_RT[0]) - (g_RT[6] + g_RT[1]));
        Kc_old[17] = exp((g_RT[5] + g_RT[1]) - (g_RT[2] + g_RT[4]));
        Kc_old[18] = exp((g_RT[5] + g_RT[2]) - (g_RT[6] + g_RT[4]));
        Kc_old[19] = exp((g_RT[5] + g_RT[2]) - (g_RT[6] + g_RT[4]));
        Kc_old[20] = exp((g_RT[5] + g_RT[2]) - (g_RT[6] + g_RT[4]));
        Kc_old[21] = exp((g_RT[5] + g_RT[5]) - (g_RT[7] + g_RT[4]));
        Kc_old[22] = exp((g_RT[5] + g_RT[5]) - (g_RT[7] + g_RT[4]));
        Kc_old[23] = refC * exp((g_RT[7]) - (g_RT[2] + g_RT[2]));
        Kc_old[24] = exp((g_RT[7] + g_RT[0]) - (g_RT[6] + g_RT[2]));
        Kc_old[25] = exp((g_RT[7] + g_RT[0]) - (g_RT[5] + g_RT[3]));
        Kc_old[26] = exp((g_RT[7] + g_RT[1]) - (g_RT[5] + g_RT[2]));
        Kc_old[27] = exp((g_RT[7] + g_RT[2]) - (g_RT[6] + g_RT[5]));
        Kc_old[28] = exp((g_RT[7] + g_RT[2]) - (g_RT[6] + g_RT[5]));
        Kc_old[29] = 1.0 / (refC) * exp((g_RT[8] + g_RT[1]) - (g_RT[9]));
        Kc_old[30] = exp((g_RT[8] + g_RT[5]) - (g_RT[9] + g_RT[2]));
        Kc_old[31] = exp((g_RT[9] + g_RT[0]) - (g_RT[8] + g_RT[2]));
        Kc_old[32] = exp((g_RT[9] + g_RT[1]) - (g_RT[8] + g_RT[4]));
        Kc_old[33] = refC * exp((g_RT[9] + g_RT[9]) - (g_RT[8] + g_RT[8] + g_RT[4]));
        Kc_old[34] = refC * exp((g_RT[10]) - (g_RT[14] + g_RT[1]));
        Kc_old[35] = exp((g_RT[10] + g_RT[0]) - (g_RT[14] + g_RT[2]));
        Kc_old[36] = exp((g_RT[10] + g_RT[0]) - (g_RT[14] + g_RT[2]));
        Kc_old[37] = exp((g_RT[10] + g_RT[1]) - (g_RT[8] + g_RT[8]));
        Kc_old[38] = exp((g_RT[10] + g_RT[1]) - (g_RT[14] + g_RT[4]));
        Kc_old[39] = exp((g_RT[11] + g_RT[0]) - (g_RT[12] + g_RT[3]));
        Kc_old[40] = exp((g_RT[11] + g_RT[1]) - (g_RT[8] + g_RT[0]));
        Kc_old[41] = exp((g_RT[11] + g_RT[2]) - (g_RT[12] + g_RT[6]));
        Kc_old[42] = exp((g_RT[11] + g_RT[4]) - (g_RT[8] + g_RT[2]));
        Kc_old[43] = exp((g_RT[11] + g_RT[8]) - (g_RT[10] + g_RT[0]));
        Kc_old[44] = exp((g_RT[11] + g_RT[8]) - (g_RT[10] + g_RT[0]));
        Kc_old[45] = exp((g_RT[11] + g_RT[8]) - (g_RT[14] + g_RT[2]));
        Kc_old[46] = exp((g_RT[11] + g_RT[9]) - (g_RT[10] + g_RT[2]));
        Kc_old[47] = exp((g_RT[12] + g_RT[2]) - (g_RT[8] + g_RT[0]));
        Kc_old[48] = exp((g_RT[12] + g_RT[4]) - (g_RT[8] + g_RT[1]));
        Kc_old[49] = exp((g_RT[12] + g_RT[8]) - (g_RT[14] + g_RT[1]));
        Kc_old[50] = refC * exp((g_RT[13]) - (g_RT[14] + g_RT[0]));
        Kc_old[51] = exp((g_RT[13] + g_RT[0]) - (g_RT[14] + g_RT[3]));
        Kc_old[52] = exp((g_RT[13] + g_RT[1]) - (g_RT[10] + g_RT[0]));
        Kc_old[53] = exp((g_RT[13] + g_RT[1]) - (g_RT[14] + g_RT[2]));
        Kc_old[54] = exp((g_RT[13] + g_RT[1]) - (g_RT[11] + g_RT[8]));
        Kc_old[55] = exp((g_RT[13] + g_RT[2]) - (g_RT[14] + g_RT[6]));
        Kc_old[56] = exp((g_RT[13] + g_RT[4]) - (g_RT[14] + g_RT[5]));
        Kc_old[57] = refC * exp((g_RT[13] + g_RT[4]) - (g_RT[14] + g_RT[0] + g_RT[4]));
    }

    /*reaction 1: H + O2 <=> O + OH */
    phi_f = sc[0]*sc[4];
    k_f = k_f_old[0];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[2];
    Kc = Kc_old[0];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 2: H + H + M <=> H2 + M */
    phi_f = sc[0]*sc[0];
    alpha = mixture + -1*sc[14] + -1*sc[6] + -1*sc[3];
    k_f = alpha * k_f_old[1];
    q_f = phi_f * k_f;
    phi_r = sc[3];
    Kc = Kc_old[1];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    phi_f = sc[0]*sc[0]*sc[14];
    k_f = k_f_old[2];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[14];
    Kc = Kc_old[2];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 4: H + H + H2 <=> H2 + H2 */
    phi_f = sc[0]*sc[0]*sc[3];
    k_f = k_f_old[3];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[3];
    Kc = Kc_old[3];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 5: H + H + H2O <=> H2 + H2O */
    phi_f = sc[0]*sc[0]*sc[6];
    k_f = k_f_old[4];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[6];
    Kc = Kc_old[4];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 6: H + O + M <=> OH + M */
    phi_f = sc[0]*sc[1];
    alpha = mixture + 4*sc[6];
    k_f = alpha * k_f_old[5];
    q_f = phi_f * k_f;
    phi_r = sc[2];
    Kc = Kc_old[5];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 7: H + O2 (+M) <=> HO2 (+M) */
    phi_f = sc[0]*sc[4];
    alpha = mixture + -1*sc[14] + 10*sc[6] + sc[3] + -0.22*sc[4];
    k_f = k_f_old[6];
    redP = 1e-12 * alpha / k_f * 3.5e+16*exp(-0.41*tc[0]+561.651/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-1e-30))+ (0.5*exp(T/-1e+30)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = Kc_old[6];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 8: H + O2 (+N2) <=> HO2 (+N2) */
    phi_f = sc[0]*sc[4];
    alpha = sc[14];
    k_f = k_f_old[7];
    redP = 1e-12 * alpha / k_f * 6.37e+20*exp(-1.72*tc[0]-261.701/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2*exp(T/-1e-30))+ (0.8*exp(T/-1e+30)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = Kc_old[7];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 9: O + O + M <=> O2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + 0.5*sc[14] + 0.5*sc[4] + 9*sc[6];
    k_f = alpha * k_f_old[8];
    q_f = phi_f * k_f;
    phi_r = sc[4];
    Kc = Kc_old[8];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 10: O + H2 <=> OH + H */
    phi_f = sc[1]*sc[3];
    k_f = k_f_old[9];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[0];
    Kc = Kc_old[9];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 11: O + H2 <=> OH + H */
    phi_f = sc[1]*sc[3];
    k_f = k_f_old[10];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[0];
    Kc = Kc_old[10];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 12: OH + OH <=> O + H2O */
    phi_f = sc[2]*sc[2];
    k_f = k_f_old[11];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[6];
    Kc = Kc_old[11];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 13: OH + H + M <=> H2O + M */
    phi_f = sc[2]*sc[0];
    alpha = mixture + -0.27*sc[3] + 11*sc[6];
    k_f = alpha * k_f_old[12];
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = Kc_old[12];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 14: OH + H2 <=> H + H2O */
    phi_f = sc[2]*sc[3];
    k_f = k_f_old[13];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[6];
    Kc = Kc_old[13];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 15: H2 + O2 <=> HO2 + H */
    phi_f = sc[3]*sc[4];
    k_f = k_f_old[14];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[0];
    Kc = Kc_old[14];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 16: HO2 + H <=> OH + OH */
    phi_f = sc[5]*sc[0];
    k_f = k_f_old[15];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[2];
    Kc = Kc_old[15];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 17: HO2 + H <=> H2O + O */
    phi_f = sc[5]*sc[0];
    k_f = k_f_old[16];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[1];
    Kc = Kc_old[16];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 18: HO2 + O <=> OH + O2 */
    phi_f = sc[5]*sc[1];
    k_f = k_f_old[17];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = Kc_old[17];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 19: HO2 + OH <=> H2O + O2 */
    phi_f = sc[5]*sc[2];
    k_f = k_f_old[18];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[4];
    Kc = Kc_old[18];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    phi_f = sc[5]*sc[2];
    k_f = k_f_old[19];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[4];
    Kc = Kc_old[19];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    phi_f = sc[5]*sc[2];
    k_f = k_f_old[20];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[4];
    Kc = Kc_old[20];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[5]*sc[5];
    k_f = k_f_old[21];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[4];
    Kc = Kc_old[21];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[5]*sc[5];
    k_f = k_f_old[22];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[4];
    Kc = Kc_old[22];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 24: H2O2 (+M) <=> OH + OH (+M) */
    phi_f = sc[7];
    alpha = mixture + 11*sc[6] + 1.5*sc[3];
    k_f = k_f_old[23];
    redP = 1e-6 * alpha / k_f * 2.291e+16*exp(-21959.3689/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-1e-30))+ (0.5*exp(T/-1e+30))+ (exp(-1e+30/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[2];
    Kc = Kc_old[23];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 25: H2O2 + H <=> H2O + OH */
    phi_f = sc[7]*sc[0];
    k_f = k_f_old[24];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[2];
    Kc = Kc_old[24];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 26: H2O2 + H <=> HO2 + H2 */
    phi_f = sc[7]*sc[0];
    k_f = k_f_old[25];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[3];
    Kc = Kc_old[25];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 27: H2O2 + O <=> HO2 + OH */
    phi_f = sc[7]*sc[1];
    k_f = k_f_old[26];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[2];
    Kc = Kc_old[26];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 28: H2O2 + OH <=> H2O + HO2 */
    phi_f = sc[7]*sc[2];
    k_f = k_f_old[27];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = Kc_old[27];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 29: H2O2 + OH <=> H2O + HO2 */
    phi_f = sc[7]*sc[2];
    k_f = k_f_old[28];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = Kc_old[28];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 30: NO + O (+M) <=> NO2 (+M) */
    phi_f = sc[8]*sc[1];
    alpha = mixture;
    k_f = k_f_old[29];
    redP = 1e-12 * alpha / k_f * 4.72e+24*exp(-2.87*tc[0]-780.07/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.12*exp(T/-1000))+ (0.88*exp(T/-10000))+ (exp(-1e+30/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[9];
    Kc = Kc_old[29];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 31: NO + HO2 <=> NO2 + OH */
    phi_f = sc[8]*sc[5];
    k_f = k_f_old[30];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[2];
    Kc = Kc_old[30];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 32: NO2 + H <=> NO + OH */
    phi_f = sc[9]*sc[0];
    k_f = k_f_old[31];
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[2];
    Kc = Kc_old[31];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 33: NO2 + O <=> NO + O2 */
    phi_f = sc[9]*sc[1];
    k_f = k_f_old[32];
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[4];
    Kc = Kc_old[32];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 34: NO2 + NO2 <=> NO + NO + O2 */
    phi_f = sc[9]*sc[9];
    k_f = k_f_old[33];
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[8]*sc[4];
    Kc = Kc_old[33];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[9] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    phi_f = sc[10];
    alpha = mixture + 0.7*sc[14] + 0.4*sc[4] + 11*sc[6];
    k_f = k_f_old[34];
    redP = 1e-6 * alpha / k_f * 4e+14*exp(-28485.2/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[1];
    Kc = Kc_old[34];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 36: N2O + H <=> N2 + OH */
    phi_f = sc[10]*sc[0];
    k_f = k_f_old[35];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[2];
    Kc = Kc_old[35];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 37: N2O + H <=> N2 + OH */
    phi_f = sc[10]*sc[0];
    k_f = k_f_old[36];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[2];
    Kc = Kc_old[36];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 38: N2O + O <=> NO + NO */
    phi_f = sc[10]*sc[1];
    k_f = k_f_old[37];
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[8];
    Kc = Kc_old[37];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 39: N2O + O <=> N2 + O2 */
    phi_f = sc[10]*sc[1];
    k_f = k_f_old[38];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[4];
    Kc = Kc_old[38];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 40: NH + H <=> N + H2 */
    phi_f = sc[11]*sc[0];
    k_f = k_f_old[39];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[3];
    Kc = Kc_old[39];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 41: NH + O <=> NO + H */
    phi_f = sc[11]*sc[1];
    k_f = k_f_old[40];
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = Kc_old[40];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 42: NH + OH <=> N + H2O */
    phi_f = sc[11]*sc[2];
    k_f = k_f_old[41];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[6];
    Kc = Kc_old[41];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 43: NH + O2 <=> NO + OH */
    phi_f = sc[11]*sc[4];
    k_f = k_f_old[42];
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[2];
    Kc = Kc_old[42];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 44: NH + NO <=> N2O + H */
    phi_f = sc[11]*sc[8];
    k_f = k_f_old[43];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[0];
    Kc = Kc_old[43];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 45: NH + NO <=> N2O + H */
    phi_f = sc[11]*sc[8];
    k_f = k_f_old[44];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[0];
    Kc = Kc_old[44];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 46: NH + NO <=> N2 + OH */
    phi_f = sc[11]*sc[8];
    k_f = k_f_old[45];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[2];
    Kc = Kc_old[45];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 47: NH + NO2 <=> N2O + OH */
    phi_f = sc[11]*sc[9];
    k_f = k_f_old[46];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[2];
    Kc = Kc_old[46];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[9] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 48: N + OH <=> NO + H */
    phi_f = sc[12]*sc[2];
    k_f = k_f_old[47];
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = Kc_old[47];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 49: N + O2 <=> NO + O */
    phi_f = sc[12]*sc[4];
    k_f = k_f_old[48];
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[1];
    Kc = Kc_old[48];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 50: N + NO <=> N2 + O */
    phi_f = sc[12]*sc[8];
    k_f = k_f_old[49];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[1];
    Kc = Kc_old[49];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 51: NNH <=> N2 + H */
    phi_f = sc[13];
    k_f = k_f_old[50];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[0];
    Kc = Kc_old[50];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 52: NNH + H <=> N2 + H2 */
    phi_f = sc[13]*sc[0];
    k_f = k_f_old[51];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[3];
    Kc = Kc_old[51];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 53: NNH + O <=> N2O + H */
    phi_f = sc[13]*sc[1];
    k_f = k_f_old[52];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[0];
    Kc = Kc_old[52];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 54: NNH + O <=> N2 + OH */
    phi_f = sc[13]*sc[1];
    k_f = k_f_old[53];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[2];
    Kc = Kc_old[53];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 55: NNH + O <=> NH + NO */
    phi_f = sc[13]*sc[1];
    k_f = k_f_old[54];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[8];
    Kc = Kc_old[54];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 56: NNH + OH <=> N2 + H2O */
    phi_f = sc[13]*sc[2];
    k_f = k_f_old[55];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[6];
    Kc = Kc_old[55];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    phi_f = sc[13]*sc[4];
    k_f = k_f_old[56];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[5];
    Kc = Kc_old[56];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    phi_f = sc[13]*sc[4];
    k_f = k_f_old[57];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[0]*sc[4];
    Kc = Kc_old[57];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[4] += 1 * qdot;

    return;
}

/*compute the progress rate for each reaction */
void progressRate(double * qdot, double * sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[15];                /*Gibbs free energy */
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

    double X;                       /*SRI intermediate */
    double F_sri;                   /*SRI intermediate */
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.314 / T;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 15; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*reaction 1: H + O2 <=> O + OH */
    phi_f = sc[0]*sc[4];
    k_f = 1e-06 * 3.6e+15*exp(-0.41*tc[0]-8354.3/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[2];
    Kc = exp((g_RT[0] + g_RT[4]) - (g_RT[1] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[0] = q_f - q_r;

    /*reaction 2: H + H + M <=> H2 + M */
    phi_f = sc[0]*sc[0];
    alpha = mixture + -1*sc[14] + -1*sc[6] + -1*sc[3];
    k_f = 1e-12 * alpha * 7e+17*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[3];
    Kc = 1.0 / (refC) * exp((g_RT[0] + g_RT[0]) - (g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[1] = q_f - q_r;

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    phi_f = sc[0]*sc[0]*sc[14];
    k_f = 1e-12 * 5.4e+18*exp(-1.3*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[14];
    Kc = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[14]) - (g_RT[3] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[2] = q_f - q_r;

    /*reaction 4: H + H + H2 <=> H2 + H2 */
    phi_f = sc[0]*sc[0]*sc[3];
    k_f = 1e-12 * 1e+17*exp(-0.6*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[3];
    Kc = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[3]) - (g_RT[3] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[3] = q_f - q_r;

    /*reaction 5: H + H + H2O <=> H2 + H2O */
    phi_f = sc[0]*sc[0]*sc[6];
    k_f = 1e-12 * 1e+19*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[6];
    Kc = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[6]) - (g_RT[3] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[4] = q_f - q_r;

    /*reaction 6: H + O + M <=> OH + M */
    phi_f = sc[0]*sc[1];
    alpha = mixture + 4*sc[6];
    k_f = 1e-12 * alpha * 6.2e+16*exp(-0.6*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[2];
    Kc = 1.0 / (refC) * exp((g_RT[0] + g_RT[1]) - (g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[5] = q_f - q_r;

    /*reaction 7: H + O2 (+M) <=> HO2 (+M) */
    phi_f = sc[0]*sc[4];
    alpha = mixture + -1*sc[14] + 10*sc[6] + sc[3] + -0.22*sc[4];
    k_f = 1e-06 * 1.5e+12*exp(0.6*tc[0]);
    redP = 1e-12 * alpha / k_f * 3.5e+16*exp(-0.41*tc[0]+561.651/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-1e-30))+ (0.5*exp(T/-1e+30)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = 1.0 / (refC) * exp((g_RT[0] + g_RT[4]) - (g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[6] = q_f - q_r;

    /*reaction 8: H + O2 (+N2) <=> HO2 (+N2) */
    phi_f = sc[0]*sc[4];
    alpha = sc[14];
    k_f = 1e-06 * 1.5e+12*exp(0.6*tc[0]);
    redP = 1e-12 * alpha / k_f * 6.37e+20*exp(-1.72*tc[0]-261.701/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2*exp(T/-1e-30))+ (0.8*exp(T/-1e+30)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = 1.0 / (refC) * exp((g_RT[0] + g_RT[4]) - (g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[7] = q_f - q_r;

    /*reaction 9: O + O + M <=> O2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + 0.5*sc[14] + 0.5*sc[4] + 9*sc[6];
    k_f = 1e-12 * alpha * 1.9e+13*exp(+899.849/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[1]) - (g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[8] = q_f - q_r;

    /*reaction 10: O + H2 <=> OH + H */
    phi_f = sc[1]*sc[3];
    k_f = 1e-06 * 3.8e+12*exp(-4000/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[0];
    Kc = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[9] = q_f - q_r;

    /*reaction 11: O + H2 <=> OH + H */
    phi_f = sc[1]*sc[3];
    k_f = 1e-06 * 8.8e+14*exp(-9650.23/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[0];
    Kc = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[10] = q_f - q_r;

    /*reaction 12: OH + OH <=> O + H2O */
    phi_f = sc[2]*sc[2];
    k_f = 1e-06 * 4300*exp(2.7*tc[0]+916.96/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[6];
    Kc = exp((g_RT[2] + g_RT[2]) - (g_RT[1] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[11] = q_f - q_r;

    /*reaction 13: OH + H + M <=> H2O + M */
    phi_f = sc[2]*sc[0];
    alpha = mixture + -0.27*sc[3] + 11*sc[6];
    k_f = 1e-12 * alpha * 4.5e+22*exp(-2*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = 1.0 / (refC) * exp((g_RT[2] + g_RT[0]) - (g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[12] = q_f - q_r;

    /*reaction 14: OH + H2 <=> H + H2O */
    phi_f = sc[2]*sc[3];
    k_f = 1e-06 * 2.1e+08*exp(1.52*tc[0]-1735.78/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[6];
    Kc = exp((g_RT[2] + g_RT[3]) - (g_RT[0] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[13] = q_f - q_r;

    /*reaction 15: H2 + O2 <=> HO2 + H */
    phi_f = sc[3]*sc[4];
    k_f = 1e-06 * 740000*exp(2.433*tc[0]-26926/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[0];
    Kc = exp((g_RT[3] + g_RT[4]) - (g_RT[5] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[14] = q_f - q_r;

    /*reaction 16: HO2 + H <=> OH + OH */
    phi_f = sc[5]*sc[0];
    k_f = 1e-06 * 8.4e+13*exp(-201.309/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[2];
    Kc = exp((g_RT[5] + g_RT[0]) - (g_RT[2] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[15] = q_f - q_r;

    /*reaction 17: HO2 + H <=> H2O + O */
    phi_f = sc[5]*sc[0];
    k_f = 1e-06 * 1.4e+12;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[1];
    Kc = exp((g_RT[5] + g_RT[0]) - (g_RT[6] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[16] = q_f - q_r;

    /*reaction 18: HO2 + O <=> OH + O2 */
    phi_f = sc[5]*sc[1];
    k_f = 1e-06 * 1.6e+13*exp(+223.956/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = exp((g_RT[5] + g_RT[1]) - (g_RT[2] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[17] = q_f - q_r;

    /*reaction 19: HO2 + OH <=> H2O + O2 */
    phi_f = sc[5]*sc[2];
    k_f = 1e-06 * 3.6e+21*exp(-2.1*tc[0]-4529.44/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[4];
    Kc = exp((g_RT[5] + g_RT[2]) - (g_RT[6] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[18] = q_f - q_r;

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    phi_f = sc[5]*sc[2];
    k_f = 1e-06 * 2e+15*exp(-0.6*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[4];
    Kc = exp((g_RT[5] + g_RT[2]) - (g_RT[6] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[19] = q_f - q_r;

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    phi_f = sc[5]*sc[2];
    k_f = 1e-06 * -2.2e+96*exp(-24*tc[0]-24660.3/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[4];
    Kc = exp((g_RT[5] + g_RT[2]) - (g_RT[6] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[20] = q_f - q_r;

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[5]*sc[5];
    k_f = 1e-06 * 1.9e+11*exp(+708.606/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[4];
    Kc = exp((g_RT[5] + g_RT[5]) - (g_RT[7] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[21] = q_f - q_r;

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[5]*sc[5];
    k_f = 1e-06 * 1e+14*exp(-5553.1/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[4];
    Kc = exp((g_RT[5] + g_RT[5]) - (g_RT[7] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[22] = q_f - q_r;

    /*reaction 24: H2O2 (+M) <=> OH + OH (+M) */
    phi_f = sc[7];
    alpha = mixture + 11*sc[6] + 1.5*sc[3];
    k_f = 1 * 4e+11*exp(-18690/tc[1]);
    redP = 1e-6 * alpha / k_f * 2.291e+16*exp(-21959.3689/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-1e-30))+ (0.5*exp(T/-1e+30))+ (exp(-1e+30/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[2];
    Kc = refC * exp((g_RT[7]) - (g_RT[2] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[23] = q_f - q_r;

    /*reaction 25: H2O2 + H <=> H2O + OH */
    phi_f = sc[7]*sc[0];
    k_f = 1e-06 * 1e+13*exp(-1801.71/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[2];
    Kc = exp((g_RT[7] + g_RT[0]) - (g_RT[6] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[24] = q_f - q_r;

    /*reaction 26: H2O2 + H <=> HO2 + H2 */
    phi_f = sc[7]*sc[0];
    k_f = 1e-06 * 1.7e+12*exp(-1892.3/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[3];
    Kc = exp((g_RT[7] + g_RT[0]) - (g_RT[5] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[25] = q_f - q_r;

    /*reaction 27: H2O2 + O <=> HO2 + OH */
    phi_f = sc[7]*sc[1];
    k_f = 1e-06 * 9.6e+06*exp(2*tc[0]-1997.99/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[2];
    Kc = exp((g_RT[7] + g_RT[1]) - (g_RT[5] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[26] = q_f - q_r;

    /*reaction 28: H2O2 + OH <=> H2O + HO2 */
    phi_f = sc[7]*sc[2];
    k_f = 1e-06 * 1.9e+12*exp(-214.897/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = exp((g_RT[7] + g_RT[2]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[27] = q_f - q_r;

    /*reaction 29: H2O2 + OH <=> H2O + HO2 */
    phi_f = sc[7]*sc[2];
    k_f = 1e-06 * 1.6e+18*exp(-14801.2/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = exp((g_RT[7] + g_RT[2]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[28] = q_f - q_r;

    /*reaction 30: NO + O (+M) <=> NO2 (+M) */
    phi_f = sc[8]*sc[1];
    alpha = mixture;
    k_f = 1e-06 * 1.3e+15*exp(-0.75*tc[0]);
    redP = 1e-12 * alpha / k_f * 4.72e+24*exp(-2.87*tc[0]-780.07/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.12*exp(T/-1000))+ (0.88*exp(T/-10000))+ (exp(-1e+30/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[9];
    Kc = 1.0 / (refC) * exp((g_RT[8] + g_RT[1]) - (g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[29] = q_f - q_r;

    /*reaction 31: NO + HO2 <=> NO2 + OH */
    phi_f = sc[8]*sc[5];
    k_f = 1e-06 * 2.1e+12*exp(+250.126/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[2];
    Kc = exp((g_RT[8] + g_RT[5]) - (g_RT[9] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[30] = q_f - q_r;

    /*reaction 32: NO2 + H <=> NO + OH */
    phi_f = sc[9]*sc[0];
    k_f = 1e-06 * 1.3e+14*exp(-182.184/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[2];
    Kc = exp((g_RT[9] + g_RT[0]) - (g_RT[8] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[31] = q_f - q_r;

    /*reaction 33: NO2 + O <=> NO + O2 */
    phi_f = sc[9]*sc[1];
    k_f = 1e-06 * 1.1e+14*exp(-0.52*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[4];
    Kc = exp((g_RT[9] + g_RT[1]) - (g_RT[8] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[32] = q_f - q_r;

    /*reaction 34: NO2 + NO2 <=> NO + NO + O2 */
    phi_f = sc[9]*sc[9];
    k_f = 1e-06 * 4.5e+12*exp(-13889.8/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[8]*sc[4];
    Kc = refC * exp((g_RT[9] + g_RT[9]) - (g_RT[8] + g_RT[8] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[33] = q_f - q_r;

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    phi_f = sc[10];
    alpha = mixture + 0.7*sc[14] + 0.4*sc[4] + 11*sc[6];
    k_f = 1 * 1.3e+12*exp(-31489.7/tc[1]);
    redP = 1e-6 * alpha / k_f * 4e+14*exp(-28485.2/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[1];
    Kc = refC * exp((g_RT[10]) - (g_RT[14] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[34] = q_f - q_r;

    /*reaction 36: N2O + H <=> N2 + OH */
    phi_f = sc[10]*sc[0];
    k_f = 1e-06 * 3.3e+10*exp(-2379.97/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[2];
    Kc = exp((g_RT[10] + g_RT[0]) - (g_RT[14] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[35] = q_f - q_r;

    /*reaction 37: N2O + H <=> N2 + OH */
    phi_f = sc[10]*sc[0];
    k_f = 1e-06 * 4.4e+14*exp(-9689.98/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[2];
    Kc = exp((g_RT[10] + g_RT[0]) - (g_RT[14] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[36] = q_f - q_r;

    /*reaction 38: N2O + O <=> NO + NO */
    phi_f = sc[10]*sc[1];
    k_f = 1e-06 * 9.2e+13*exp(-13930/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[8];
    Kc = exp((g_RT[10] + g_RT[1]) - (g_RT[8] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[37] = q_f - q_r;

    /*reaction 39: N2O + O <=> N2 + O2 */
    phi_f = sc[10]*sc[1];
    k_f = 1e-06 * 3.7e+12*exp(-8020.13/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[4];
    Kc = exp((g_RT[10] + g_RT[1]) - (g_RT[14] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[38] = q_f - q_r;

    /*reaction 40: NH + H <=> N + H2 */
    phi_f = sc[11]*sc[0];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[3];
    Kc = exp((g_RT[11] + g_RT[0]) - (g_RT[12] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[39] = q_f - q_r;

    /*reaction 41: NH + O <=> NO + H */
    phi_f = sc[11]*sc[1];
    k_f = 1e-06 * 9.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = exp((g_RT[11] + g_RT[1]) - (g_RT[8] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[40] = q_f - q_r;

    /*reaction 42: NH + OH <=> N + H2O */
    phi_f = sc[11]*sc[2];
    k_f = 1e-06 * 5e+11*exp(0.5*tc[0]-1006.54/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[6];
    Kc = exp((g_RT[11] + g_RT[2]) - (g_RT[12] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[41] = q_f - q_r;

    /*reaction 43: NH + O2 <=> NO + OH */
    phi_f = sc[11]*sc[4];
    k_f = 1e-06 * 1.3e+06*exp(1.5*tc[0]-50.3271/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[2];
    Kc = exp((g_RT[11] + g_RT[4]) - (g_RT[8] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[42] = q_f - q_r;

    /*reaction 44: NH + NO <=> N2O + H */
    phi_f = sc[11]*sc[8];
    k_f = 1e-06 * 2.9e+14*exp(-0.4*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[0];
    Kc = exp((g_RT[11] + g_RT[8]) - (g_RT[10] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[43] = q_f - q_r;

    /*reaction 45: NH + NO <=> N2O + H */
    phi_f = sc[11]*sc[8];
    k_f = 1e-06 * -2.2e+13*exp(-0.23*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[0];
    Kc = exp((g_RT[11] + g_RT[8]) - (g_RT[10] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[44] = q_f - q_r;

    /*reaction 46: NH + NO <=> N2 + OH */
    phi_f = sc[11]*sc[8];
    k_f = 1e-06 * 2.2e+13*exp(-0.23*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[2];
    Kc = exp((g_RT[11] + g_RT[8]) - (g_RT[14] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[45] = q_f - q_r;

    /*reaction 47: NH + NO2 <=> N2O + OH */
    phi_f = sc[11]*sc[9];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[2];
    Kc = exp((g_RT[11] + g_RT[9]) - (g_RT[10] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[46] = q_f - q_r;

    /*reaction 48: N + OH <=> NO + H */
    phi_f = sc[12]*sc[2];
    k_f = 1e-06 * 3.8e+13;
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = exp((g_RT[12] + g_RT[2]) - (g_RT[8] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[47] = q_f - q_r;

    /*reaction 49: N + O2 <=> NO + O */
    phi_f = sc[12]*sc[4];
    k_f = 1e-06 * 6.4e+09*exp(1*tc[0]-3160.54/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[1];
    Kc = exp((g_RT[12] + g_RT[4]) - (g_RT[8] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[48] = q_f - q_r;

    /*reaction 50: N + NO <=> N2 + O */
    phi_f = sc[12]*sc[8];
    k_f = 1e-06 * 2.1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[1];
    Kc = exp((g_RT[12] + g_RT[8]) - (g_RT[14] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[49] = q_f - q_r;

    /*reaction 51: NNH <=> N2 + H */
    phi_f = sc[13];
    k_f = 1 * 6.5e+07;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[0];
    Kc = refC * exp((g_RT[13]) - (g_RT[14] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[50] = q_f - q_r;

    /*reaction 52: NNH + H <=> N2 + H2 */
    phi_f = sc[13]*sc[0];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[3];
    Kc = exp((g_RT[13] + g_RT[0]) - (g_RT[14] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[51] = q_f - q_r;

    /*reaction 53: NNH + O <=> N2O + H */
    phi_f = sc[13]*sc[1];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[0];
    Kc = exp((g_RT[13] + g_RT[1]) - (g_RT[10] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[52] = q_f - q_r;

    /*reaction 54: NNH + O <=> N2 + OH */
    phi_f = sc[13]*sc[1];
    k_f = 1e-06 * 8e+13;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[2];
    Kc = exp((g_RT[13] + g_RT[1]) - (g_RT[14] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[53] = q_f - q_r;

    /*reaction 55: NNH + O <=> NH + NO */
    phi_f = sc[13]*sc[1];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[8];
    Kc = exp((g_RT[13] + g_RT[1]) - (g_RT[11] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[54] = q_f - q_r;

    /*reaction 56: NNH + OH <=> N2 + H2O */
    phi_f = sc[13]*sc[2];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[6];
    Kc = exp((g_RT[13] + g_RT[2]) - (g_RT[14] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[55] = q_f - q_r;

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    phi_f = sc[13]*sc[4];
    k_f = 1e-06 * 2e+14;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[5];
    Kc = exp((g_RT[13] + g_RT[4]) - (g_RT[14] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[56] = q_f - q_r;

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    phi_f = sc[13]*sc[4];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[0]*sc[4];
    Kc = refC * exp((g_RT[13] + g_RT[4]) - (g_RT[14] + g_RT[0] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[57] = q_f - q_r;

    return;
}


/*compute the progress rate for each reaction */
void progressRateFR(double * q_f, double * q_r, double * sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[15];                /*Gibbs free energy */
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

    double X;                       /*SRI intermediate */
    double F_sri;                   /*SRI intermediate */
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.314 / T;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 15; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*reaction 1: H + O2 <=> O + OH */
    phi_f = sc[0]*sc[4];
    k_f = 1e-06 * 3.6e+15*exp(-0.41*tc[0]-8354.3/tc[1]);
    q_f[0] = phi_f * k_f;
    phi_r = sc[1]*sc[2];
    Kc = exp((g_RT[0] + g_RT[4]) - (g_RT[1] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[0] = phi_r * k_r;

    /*reaction 2: H + H + M <=> H2 + M */
    phi_f = sc[0]*sc[0];
    alpha = mixture + -1*sc[14] + -1*sc[6] + -1*sc[3];
    k_f = 1e-12 * alpha * 7e+17*exp(-1*tc[0]);
    q_f[1] = phi_f * k_f;
    phi_r = sc[3];
    Kc = 1.0 / (refC) * exp((g_RT[0] + g_RT[0]) - (g_RT[3]));
    k_r = k_f / Kc;
    q_r[1] = phi_r * k_r;

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    phi_f = sc[0]*sc[0]*sc[14];
    k_f = 1e-12 * 5.4e+18*exp(-1.3*tc[0]);
    q_f[2] = phi_f * k_f;
    phi_r = sc[3]*sc[14];
    Kc = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[14]) - (g_RT[3] + g_RT[14]));
    k_r = k_f / Kc;
    q_r[2] = phi_r * k_r;

    /*reaction 4: H + H + H2 <=> H2 + H2 */
    phi_f = sc[0]*sc[0]*sc[3];
    k_f = 1e-12 * 1e+17*exp(-0.6*tc[0]);
    q_f[3] = phi_f * k_f;
    phi_r = sc[3]*sc[3];
    Kc = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[3]) - (g_RT[3] + g_RT[3]));
    k_r = k_f / Kc;
    q_r[3] = phi_r * k_r;

    /*reaction 5: H + H + H2O <=> H2 + H2O */
    phi_f = sc[0]*sc[0]*sc[6];
    k_f = 1e-12 * 1e+19*exp(-1*tc[0]);
    q_f[4] = phi_f * k_f;
    phi_r = sc[3]*sc[6];
    Kc = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[6]) - (g_RT[3] + g_RT[6]));
    k_r = k_f / Kc;
    q_r[4] = phi_r * k_r;

    /*reaction 6: H + O + M <=> OH + M */
    phi_f = sc[0]*sc[1];
    alpha = mixture + 4*sc[6];
    k_f = 1e-12 * alpha * 6.2e+16*exp(-0.6*tc[0]);
    q_f[5] = phi_f * k_f;
    phi_r = sc[2];
    Kc = 1.0 / (refC) * exp((g_RT[0] + g_RT[1]) - (g_RT[2]));
    k_r = k_f / Kc;
    q_r[5] = phi_r * k_r;

    /*reaction 7: H + O2 (+M) <=> HO2 (+M) */
    phi_f = sc[0]*sc[4];
    alpha = mixture + -1*sc[14] + 10*sc[6] + sc[3] + -0.22*sc[4];
    k_f = 1e-06 * 1.5e+12*exp(0.6*tc[0]);
    redP = 1e-12 * alpha / k_f * 3.5e+16*exp(-0.41*tc[0]+561.651/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-1e-30))+ (0.5*exp(T/-1e+30)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[6] = phi_f * k_f;
    phi_r = sc[5];
    Kc = 1.0 / (refC) * exp((g_RT[0] + g_RT[4]) - (g_RT[5]));
    k_r = k_f / Kc;
    q_r[6] = phi_r * k_r;

    /*reaction 8: H + O2 (+N2) <=> HO2 (+N2) */
    phi_f = sc[0]*sc[4];
    alpha = sc[14];
    k_f = 1e-06 * 1.5e+12*exp(0.6*tc[0]);
    redP = 1e-12 * alpha / k_f * 6.37e+20*exp(-1.72*tc[0]-261.701/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2*exp(T/-1e-30))+ (0.8*exp(T/-1e+30)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[7] = phi_f * k_f;
    phi_r = sc[5];
    Kc = 1.0 / (refC) * exp((g_RT[0] + g_RT[4]) - (g_RT[5]));
    k_r = k_f / Kc;
    q_r[7] = phi_r * k_r;

    /*reaction 9: O + O + M <=> O2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + 0.5*sc[14] + 0.5*sc[4] + 9*sc[6];
    k_f = 1e-12 * alpha * 1.9e+13*exp(+899.849/tc[1]);
    q_f[8] = phi_f * k_f;
    phi_r = sc[4];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[1]) - (g_RT[4]));
    k_r = k_f / Kc;
    q_r[8] = phi_r * k_r;

    /*reaction 10: O + H2 <=> OH + H */
    phi_f = sc[1]*sc[3];
    k_f = 1e-06 * 3.8e+12*exp(-4000/tc[1]);
    q_f[9] = phi_f * k_f;
    phi_r = sc[2]*sc[0];
    Kc = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[9] = phi_r * k_r;

    /*reaction 11: O + H2 <=> OH + H */
    phi_f = sc[1]*sc[3];
    k_f = 1e-06 * 8.8e+14*exp(-9650.23/tc[1]);
    q_f[10] = phi_f * k_f;
    phi_r = sc[2]*sc[0];
    Kc = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[10] = phi_r * k_r;

    /*reaction 12: OH + OH <=> O + H2O */
    phi_f = sc[2]*sc[2];
    k_f = 1e-06 * 4300*exp(2.7*tc[0]+916.96/tc[1]);
    q_f[11] = phi_f * k_f;
    phi_r = sc[1]*sc[6];
    Kc = exp((g_RT[2] + g_RT[2]) - (g_RT[1] + g_RT[6]));
    k_r = k_f / Kc;
    q_r[11] = phi_r * k_r;

    /*reaction 13: OH + H + M <=> H2O + M */
    phi_f = sc[2]*sc[0];
    alpha = mixture + -0.27*sc[3] + 11*sc[6];
    k_f = 1e-12 * alpha * 4.5e+22*exp(-2*tc[0]);
    q_f[12] = phi_f * k_f;
    phi_r = sc[6];
    Kc = 1.0 / (refC) * exp((g_RT[2] + g_RT[0]) - (g_RT[6]));
    k_r = k_f / Kc;
    q_r[12] = phi_r * k_r;

    /*reaction 14: OH + H2 <=> H + H2O */
    phi_f = sc[2]*sc[3];
    k_f = 1e-06 * 2.1e+08*exp(1.52*tc[0]-1735.78/tc[1]);
    q_f[13] = phi_f * k_f;
    phi_r = sc[0]*sc[6];
    Kc = exp((g_RT[2] + g_RT[3]) - (g_RT[0] + g_RT[6]));
    k_r = k_f / Kc;
    q_r[13] = phi_r * k_r;

    /*reaction 15: H2 + O2 <=> HO2 + H */
    phi_f = sc[3]*sc[4];
    k_f = 1e-06 * 740000*exp(2.433*tc[0]-26926/tc[1]);
    q_f[14] = phi_f * k_f;
    phi_r = sc[5]*sc[0];
    Kc = exp((g_RT[3] + g_RT[4]) - (g_RT[5] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[14] = phi_r * k_r;

    /*reaction 16: HO2 + H <=> OH + OH */
    phi_f = sc[5]*sc[0];
    k_f = 1e-06 * 8.4e+13*exp(-201.309/tc[1]);
    q_f[15] = phi_f * k_f;
    phi_r = sc[2]*sc[2];
    Kc = exp((g_RT[5] + g_RT[0]) - (g_RT[2] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[15] = phi_r * k_r;

    /*reaction 17: HO2 + H <=> H2O + O */
    phi_f = sc[5]*sc[0];
    k_f = 1e-06 * 1.4e+12;
    q_f[16] = phi_f * k_f;
    phi_r = sc[6]*sc[1];
    Kc = exp((g_RT[5] + g_RT[0]) - (g_RT[6] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[16] = phi_r * k_r;

    /*reaction 18: HO2 + O <=> OH + O2 */
    phi_f = sc[5]*sc[1];
    k_f = 1e-06 * 1.6e+13*exp(+223.956/tc[1]);
    q_f[17] = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = exp((g_RT[5] + g_RT[1]) - (g_RT[2] + g_RT[4]));
    k_r = k_f / Kc;
    q_r[17] = phi_r * k_r;

    /*reaction 19: HO2 + OH <=> H2O + O2 */
    phi_f = sc[5]*sc[2];
    k_f = 1e-06 * 3.6e+21*exp(-2.1*tc[0]-4529.44/tc[1]);
    q_f[18] = phi_f * k_f;
    phi_r = sc[6]*sc[4];
    Kc = exp((g_RT[5] + g_RT[2]) - (g_RT[6] + g_RT[4]));
    k_r = k_f / Kc;
    q_r[18] = phi_r * k_r;

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    phi_f = sc[5]*sc[2];
    k_f = 1e-06 * 2e+15*exp(-0.6*tc[0]);
    q_f[19] = phi_f * k_f;
    phi_r = sc[6]*sc[4];
    Kc = exp((g_RT[5] + g_RT[2]) - (g_RT[6] + g_RT[4]));
    k_r = k_f / Kc;
    q_r[19] = phi_r * k_r;

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    phi_f = sc[5]*sc[2];
    k_f = 1e-06 * -2.2e+96*exp(-24*tc[0]-24660.3/tc[1]);
    q_f[20] = phi_f * k_f;
    phi_r = sc[6]*sc[4];
    Kc = exp((g_RT[5] + g_RT[2]) - (g_RT[6] + g_RT[4]));
    k_r = k_f / Kc;
    q_r[20] = phi_r * k_r;

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[5]*sc[5];
    k_f = 1e-06 * 1.9e+11*exp(+708.606/tc[1]);
    q_f[21] = phi_f * k_f;
    phi_r = sc[7]*sc[4];
    Kc = exp((g_RT[5] + g_RT[5]) - (g_RT[7] + g_RT[4]));
    k_r = k_f / Kc;
    q_r[21] = phi_r * k_r;

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[5]*sc[5];
    k_f = 1e-06 * 1e+14*exp(-5553.1/tc[1]);
    q_f[22] = phi_f * k_f;
    phi_r = sc[7]*sc[4];
    Kc = exp((g_RT[5] + g_RT[5]) - (g_RT[7] + g_RT[4]));
    k_r = k_f / Kc;
    q_r[22] = phi_r * k_r;

    /*reaction 24: H2O2 (+M) <=> OH + OH (+M) */
    phi_f = sc[7];
    alpha = mixture + 11*sc[6] + 1.5*sc[3];
    k_f = 1 * 4e+11*exp(-18690/tc[1]);
    redP = 1e-6 * alpha / k_f * 2.291e+16*exp(-21959.3689/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-1e-30))+ (0.5*exp(T/-1e+30))+ (exp(-1e+30/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[23] = phi_f * k_f;
    phi_r = sc[2]*sc[2];
    Kc = refC * exp((g_RT[7]) - (g_RT[2] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[23] = phi_r * k_r;

    /*reaction 25: H2O2 + H <=> H2O + OH */
    phi_f = sc[7]*sc[0];
    k_f = 1e-06 * 1e+13*exp(-1801.71/tc[1]);
    q_f[24] = phi_f * k_f;
    phi_r = sc[6]*sc[2];
    Kc = exp((g_RT[7] + g_RT[0]) - (g_RT[6] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[24] = phi_r * k_r;

    /*reaction 26: H2O2 + H <=> HO2 + H2 */
    phi_f = sc[7]*sc[0];
    k_f = 1e-06 * 1.7e+12*exp(-1892.3/tc[1]);
    q_f[25] = phi_f * k_f;
    phi_r = sc[5]*sc[3];
    Kc = exp((g_RT[7] + g_RT[0]) - (g_RT[5] + g_RT[3]));
    k_r = k_f / Kc;
    q_r[25] = phi_r * k_r;

    /*reaction 27: H2O2 + O <=> HO2 + OH */
    phi_f = sc[7]*sc[1];
    k_f = 1e-06 * 9.6e+06*exp(2*tc[0]-1997.99/tc[1]);
    q_f[26] = phi_f * k_f;
    phi_r = sc[5]*sc[2];
    Kc = exp((g_RT[7] + g_RT[1]) - (g_RT[5] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[26] = phi_r * k_r;

    /*reaction 28: H2O2 + OH <=> H2O + HO2 */
    phi_f = sc[7]*sc[2];
    k_f = 1e-06 * 1.9e+12*exp(-214.897/tc[1]);
    q_f[27] = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = exp((g_RT[7] + g_RT[2]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r[27] = phi_r * k_r;

    /*reaction 29: H2O2 + OH <=> H2O + HO2 */
    phi_f = sc[7]*sc[2];
    k_f = 1e-06 * 1.6e+18*exp(-14801.2/tc[1]);
    q_f[28] = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = exp((g_RT[7] + g_RT[2]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r[28] = phi_r * k_r;

    /*reaction 30: NO + O (+M) <=> NO2 (+M) */
    phi_f = sc[8]*sc[1];
    alpha = mixture;
    k_f = 1e-06 * 1.3e+15*exp(-0.75*tc[0]);
    redP = 1e-12 * alpha / k_f * 4.72e+24*exp(-2.87*tc[0]-780.07/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.12*exp(T/-1000))+ (0.88*exp(T/-10000))+ (exp(-1e+30/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[29] = phi_f * k_f;
    phi_r = sc[9];
    Kc = 1.0 / (refC) * exp((g_RT[8] + g_RT[1]) - (g_RT[9]));
    k_r = k_f / Kc;
    q_r[29] = phi_r * k_r;

    /*reaction 31: NO + HO2 <=> NO2 + OH */
    phi_f = sc[8]*sc[5];
    k_f = 1e-06 * 2.1e+12*exp(+250.126/tc[1]);
    q_f[0] = phi_f * k_f;
    phi_r = sc[9]*sc[2];
    Kc = exp((g_RT[8] + g_RT[5]) - (g_RT[9] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[30] = phi_r * k_r;

    /*reaction 32: NO2 + H <=> NO + OH */
    phi_f = sc[9]*sc[0];
    k_f = 1e-06 * 1.3e+14*exp(-182.184/tc[1]);
    q_f[31] = phi_f * k_f;
    phi_r = sc[8]*sc[2];
    Kc = exp((g_RT[9] + g_RT[0]) - (g_RT[8] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[31] = phi_r * k_r;

    /*reaction 33: NO2 + O <=> NO + O2 */
    phi_f = sc[9]*sc[1];
    k_f = 1e-06 * 1.1e+14*exp(-0.52*tc[0]);
    q_f[32] = phi_f * k_f;
    phi_r = sc[8]*sc[4];
    Kc = exp((g_RT[9] + g_RT[1]) - (g_RT[8] + g_RT[4]));
    k_r = k_f / Kc;
    q_r[32] = phi_r * k_r;

    /*reaction 34: NO2 + NO2 <=> NO + NO + O2 */
    phi_f = sc[9]*sc[9];
    k_f = 1e-06 * 4.5e+12*exp(-13889.8/tc[1]);
    q_f[33] = phi_f * k_f;
    phi_r = sc[8]*sc[8]*sc[4];
    Kc = refC * exp((g_RT[9] + g_RT[9]) - (g_RT[8] + g_RT[8] + g_RT[4]));
    k_r = k_f / Kc;
    q_r[33] = phi_r * k_r;

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    phi_f = sc[10];
    alpha = mixture + 0.7*sc[14] + 0.4*sc[4] + 11*sc[6];
    k_f = 1 * 1.3e+12*exp(-31489.7/tc[1]);
    redP = 1e-6 * alpha / k_f * 4e+14*exp(-28485.2/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f[34] = phi_f * k_f;
    phi_r = sc[14]*sc[1];
    Kc = refC * exp((g_RT[10]) - (g_RT[14] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[34] = phi_r * k_r;

    /*reaction 36: N2O + H <=> N2 + OH */
    phi_f = sc[10]*sc[0];
    k_f = 1e-06 * 3.3e+10*exp(-2379.97/tc[1]);
    q_f[35] = phi_f * k_f;
    phi_r = sc[14]*sc[2];
    Kc = exp((g_RT[10] + g_RT[0]) - (g_RT[14] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[35] = phi_r * k_r;

    /*reaction 37: N2O + H <=> N2 + OH */
    phi_f = sc[10]*sc[0];
    k_f = 1e-06 * 4.4e+14*exp(-9689.98/tc[1]);
    q_f[36] = phi_f * k_f;
    phi_r = sc[14]*sc[2];
    Kc = exp((g_RT[10] + g_RT[0]) - (g_RT[14] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[36] = phi_r * k_r;

    /*reaction 38: N2O + O <=> NO + NO */
    phi_f = sc[10]*sc[1];
    k_f = 1e-06 * 9.2e+13*exp(-13930/tc[1]);
    q_f[37] = phi_f * k_f;
    phi_r = sc[8]*sc[8];
    Kc = exp((g_RT[10] + g_RT[1]) - (g_RT[8] + g_RT[8]));
    k_r = k_f / Kc;
    q_r[37] = phi_r * k_r;

    /*reaction 39: N2O + O <=> N2 + O2 */
    phi_f = sc[10]*sc[1];
    k_f = 1e-06 * 3.7e+12*exp(-8020.13/tc[1]);
    q_f[38] = phi_f * k_f;
    phi_r = sc[14]*sc[4];
    Kc = exp((g_RT[10] + g_RT[1]) - (g_RT[14] + g_RT[4]));
    k_r = k_f / Kc;
    q_r[38] = phi_r * k_r;

    /*reaction 40: NH + H <=> N + H2 */
    phi_f = sc[11]*sc[0];
    k_f = 1e-06 * 3e+13;
    q_f[39] = phi_f * k_f;
    phi_r = sc[12]*sc[3];
    Kc = exp((g_RT[11] + g_RT[0]) - (g_RT[12] + g_RT[3]));
    k_r = k_f / Kc;
    q_r[39] = phi_r * k_r;

    /*reaction 41: NH + O <=> NO + H */
    phi_f = sc[11]*sc[1];
    k_f = 1e-06 * 9.2e+13;
    q_f[40] = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = exp((g_RT[11] + g_RT[1]) - (g_RT[8] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[40] = phi_r * k_r;

    /*reaction 42: NH + OH <=> N + H2O */
    phi_f = sc[11]*sc[2];
    k_f = 1e-06 * 5e+11*exp(0.5*tc[0]-1006.54/tc[1]);
    q_f[41] = phi_f * k_f;
    phi_r = sc[12]*sc[6];
    Kc = exp((g_RT[11] + g_RT[2]) - (g_RT[12] + g_RT[6]));
    k_r = k_f / Kc;
    q_r[41] = phi_r * k_r;

    /*reaction 43: NH + O2 <=> NO + OH */
    phi_f = sc[11]*sc[4];
    k_f = 1e-06 * 1.3e+06*exp(1.5*tc[0]-50.3271/tc[1]);
    q_f[42] = phi_f * k_f;
    phi_r = sc[8]*sc[2];
    Kc = exp((g_RT[11] + g_RT[4]) - (g_RT[8] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[42] = phi_r * k_r;

    /*reaction 44: NH + NO <=> N2O + H */
    phi_f = sc[11]*sc[8];
    k_f = 1e-06 * 2.9e+14*exp(-0.4*tc[0]);
    q_f[43] = phi_f * k_f;
    phi_r = sc[10]*sc[0];
    Kc = exp((g_RT[11] + g_RT[8]) - (g_RT[10] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[43] = phi_r * k_r;

    /*reaction 45: NH + NO <=> N2O + H */
    phi_f = sc[11]*sc[8];
    k_f = 1e-06 * -2.2e+13*exp(-0.23*tc[0]);
    q_f[44] = phi_f * k_f;
    phi_r = sc[10]*sc[0];
    Kc = exp((g_RT[11] + g_RT[8]) - (g_RT[10] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[44] = phi_r * k_r;

    /*reaction 46: NH + NO <=> N2 + OH */
    phi_f = sc[11]*sc[8];
    k_f = 1e-06 * 2.2e+13*exp(-0.23*tc[0]);
    q_f[45] = phi_f * k_f;
    phi_r = sc[14]*sc[2];
    Kc = exp((g_RT[11] + g_RT[8]) - (g_RT[14] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[45] = phi_r * k_r;

    /*reaction 47: NH + NO2 <=> N2O + OH */
    phi_f = sc[11]*sc[9];
    k_f = 1e-06 * 1e+13;
    q_f[46] = phi_f * k_f;
    phi_r = sc[10]*sc[2];
    Kc = exp((g_RT[11] + g_RT[9]) - (g_RT[10] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[46] = phi_r * k_r;

    /*reaction 48: N + OH <=> NO + H */
    phi_f = sc[12]*sc[2];
    k_f = 1e-06 * 3.8e+13;
    q_f[47] = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = exp((g_RT[12] + g_RT[2]) - (g_RT[8] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[47] = phi_r * k_r;

    /*reaction 49: N + O2 <=> NO + O */
    phi_f = sc[12]*sc[4];
    k_f = 1e-06 * 6.4e+09*exp(1*tc[0]-3160.54/tc[1]);
    q_f[48] = phi_f * k_f;
    phi_r = sc[8]*sc[1];
    Kc = exp((g_RT[12] + g_RT[4]) - (g_RT[8] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[48] = phi_r * k_r;

    /*reaction 50: N + NO <=> N2 + O */
    phi_f = sc[12]*sc[8];
    k_f = 1e-06 * 2.1e+13;
    q_f[49] = phi_f * k_f;
    phi_r = sc[14]*sc[1];
    Kc = exp((g_RT[12] + g_RT[8]) - (g_RT[14] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[49] = phi_r * k_r;

    /*reaction 51: NNH <=> N2 + H */
    phi_f = sc[13];
    k_f = 1 * 6.5e+07;
    q_f[50] = phi_f * k_f;
    phi_r = sc[14]*sc[0];
    Kc = refC * exp((g_RT[13]) - (g_RT[14] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[50] = phi_r * k_r;

    /*reaction 52: NNH + H <=> N2 + H2 */
    phi_f = sc[13]*sc[0];
    k_f = 1e-06 * 1e+14;
    q_f[51] = phi_f * k_f;
    phi_r = sc[14]*sc[3];
    Kc = exp((g_RT[13] + g_RT[0]) - (g_RT[14] + g_RT[3]));
    k_r = k_f / Kc;
    q_r[51] = phi_r * k_r;

    /*reaction 53: NNH + O <=> N2O + H */
    phi_f = sc[13]*sc[1];
    k_f = 1e-06 * 1e+14;
    q_f[52] = phi_f * k_f;
    phi_r = sc[10]*sc[0];
    Kc = exp((g_RT[13] + g_RT[1]) - (g_RT[10] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[52] = phi_r * k_r;

    /*reaction 54: NNH + O <=> N2 + OH */
    phi_f = sc[13]*sc[1];
    k_f = 1e-06 * 8e+13;
    q_f[53] = phi_f * k_f;
    phi_r = sc[14]*sc[2];
    Kc = exp((g_RT[13] + g_RT[1]) - (g_RT[14] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[53] = phi_r * k_r;

    /*reaction 55: NNH + O <=> NH + NO */
    phi_f = sc[13]*sc[1];
    k_f = 1e-06 * 5e+13;
    q_f[54] = phi_f * k_f;
    phi_r = sc[11]*sc[8];
    Kc = exp((g_RT[13] + g_RT[1]) - (g_RT[11] + g_RT[8]));
    k_r = k_f / Kc;
    q_r[54] = phi_r * k_r;

    /*reaction 56: NNH + OH <=> N2 + H2O */
    phi_f = sc[13]*sc[2];
    k_f = 1e-06 * 5e+13;
    q_f[55] = phi_f * k_f;
    phi_r = sc[14]*sc[6];
    Kc = exp((g_RT[13] + g_RT[2]) - (g_RT[14] + g_RT[6]));
    k_r = k_f / Kc;
    q_r[55] = phi_r * k_r;

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    phi_f = sc[13]*sc[4];
    k_f = 1e-06 * 2e+14;
    q_f[56] = phi_f * k_f;
    phi_r = sc[14]*sc[5];
    Kc = exp((g_RT[13] + g_RT[4]) - (g_RT[14] + g_RT[5]));
    k_r = k_f / Kc;
    q_r[56] = phi_r * k_r;

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    phi_f = sc[13]*sc[4];
    k_f = 1e-06 * 5e+13;
    q_f[57] = phi_f * k_f;
    phi_r = sc[14]*sc[0]*sc[4];
    Kc = refC * exp((g_RT[13] + g_RT[4]) - (g_RT[14] + g_RT[0] + g_RT[4]));
    k_r = k_f / Kc;
    q_r[57] = phi_r * k_r;

    return;
}


/*compute the equilibrium constants for each reaction */
void equilibriumConstants(double *kc, double * g_RT, double T)
{
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.314 / T;

    /*reaction 1: H + O2 <=> O + OH */
    kc[0] = exp((g_RT[0] + g_RT[4]) - (g_RT[1] + g_RT[2]));

    /*reaction 2: H + H + M <=> H2 + M */
    kc[1] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0]) - (g_RT[3]));

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    kc[2] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[14]) - (g_RT[3] + g_RT[14]));

    /*reaction 4: H + H + H2 <=> H2 + H2 */
    kc[3] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[3]) - (g_RT[3] + g_RT[3]));

    /*reaction 5: H + H + H2O <=> H2 + H2O */
    kc[4] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[6]) - (g_RT[3] + g_RT[6]));

    /*reaction 6: H + O + M <=> OH + M */
    kc[5] = 1.0 / (refC) * exp((g_RT[0] + g_RT[1]) - (g_RT[2]));

    /*reaction 7: H + O2 (+M) <=> HO2 (+M) */
    kc[6] = 1.0 / (refC) * exp((g_RT[0] + g_RT[4]) - (g_RT[5]));

    /*reaction 8: H + O2 (+N2) <=> HO2 (+N2) */
    kc[7] = 1.0 / (refC) * exp((g_RT[0] + g_RT[4]) - (g_RT[5]));

    /*reaction 9: O + O + M <=> O2 + M */
    kc[8] = 1.0 / (refC) * exp((g_RT[1] + g_RT[1]) - (g_RT[4]));

    /*reaction 10: O + H2 <=> OH + H */
    kc[9] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[0]));

    /*reaction 11: O + H2 <=> OH + H */
    kc[10] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[0]));

    /*reaction 12: OH + OH <=> O + H2O */
    kc[11] = exp((g_RT[2] + g_RT[2]) - (g_RT[1] + g_RT[6]));

    /*reaction 13: OH + H + M <=> H2O + M */
    kc[12] = 1.0 / (refC) * exp((g_RT[2] + g_RT[0]) - (g_RT[6]));

    /*reaction 14: OH + H2 <=> H + H2O */
    kc[13] = exp((g_RT[2] + g_RT[3]) - (g_RT[0] + g_RT[6]));

    /*reaction 15: H2 + O2 <=> HO2 + H */
    kc[14] = exp((g_RT[3] + g_RT[4]) - (g_RT[5] + g_RT[0]));

    /*reaction 16: HO2 + H <=> OH + OH */
    kc[15] = exp((g_RT[5] + g_RT[0]) - (g_RT[2] + g_RT[2]));

    /*reaction 17: HO2 + H <=> H2O + O */
    kc[16] = exp((g_RT[5] + g_RT[0]) - (g_RT[6] + g_RT[1]));

    /*reaction 18: HO2 + O <=> OH + O2 */
    kc[17] = exp((g_RT[5] + g_RT[1]) - (g_RT[2] + g_RT[4]));

    /*reaction 19: HO2 + OH <=> H2O + O2 */
    kc[18] = exp((g_RT[5] + g_RT[2]) - (g_RT[6] + g_RT[4]));

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    kc[19] = exp((g_RT[5] + g_RT[2]) - (g_RT[6] + g_RT[4]));

    /*reaction 21: HO2 + OH <=> H2O + O2 */
    kc[20] = exp((g_RT[5] + g_RT[2]) - (g_RT[6] + g_RT[4]));

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    kc[21] = exp((g_RT[5] + g_RT[5]) - (g_RT[7] + g_RT[4]));

    /*reaction 23: HO2 + HO2 <=> H2O2 + O2 */
    kc[22] = exp((g_RT[5] + g_RT[5]) - (g_RT[7] + g_RT[4]));

    /*reaction 24: H2O2 (+M) <=> OH + OH (+M) */
    kc[23] = refC * exp((g_RT[7]) - (g_RT[2] + g_RT[2]));

    /*reaction 25: H2O2 + H <=> H2O + OH */
    kc[24] = exp((g_RT[7] + g_RT[0]) - (g_RT[6] + g_RT[2]));

    /*reaction 26: H2O2 + H <=> HO2 + H2 */
    kc[25] = exp((g_RT[7] + g_RT[0]) - (g_RT[5] + g_RT[3]));

    /*reaction 27: H2O2 + O <=> HO2 + OH */
    kc[26] = exp((g_RT[7] + g_RT[1]) - (g_RT[5] + g_RT[2]));

    /*reaction 28: H2O2 + OH <=> H2O + HO2 */
    kc[27] = exp((g_RT[7] + g_RT[2]) - (g_RT[6] + g_RT[5]));

    /*reaction 29: H2O2 + OH <=> H2O + HO2 */
    kc[28] = exp((g_RT[7] + g_RT[2]) - (g_RT[6] + g_RT[5]));

    /*reaction 30: NO + O (+M) <=> NO2 (+M) */
    kc[29] = 1.0 / (refC) * exp((g_RT[8] + g_RT[1]) - (g_RT[9]));

    /*reaction 31: NO + HO2 <=> NO2 + OH */
    kc[30] = exp((g_RT[8] + g_RT[5]) - (g_RT[9] + g_RT[2]));

    /*reaction 32: NO2 + H <=> NO + OH */
    kc[31] = exp((g_RT[9] + g_RT[0]) - (g_RT[8] + g_RT[2]));

    /*reaction 33: NO2 + O <=> NO + O2 */
    kc[32] = exp((g_RT[9] + g_RT[1]) - (g_RT[8] + g_RT[4]));

    /*reaction 34: NO2 + NO2 <=> NO + NO + O2 */
    kc[33] = refC * exp((g_RT[9] + g_RT[9]) - (g_RT[8] + g_RT[8] + g_RT[4]));

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    kc[34] = refC * exp((g_RT[10]) - (g_RT[14] + g_RT[1]));

    /*reaction 36: N2O + H <=> N2 + OH */
    kc[35] = exp((g_RT[10] + g_RT[0]) - (g_RT[14] + g_RT[2]));

    /*reaction 37: N2O + H <=> N2 + OH */
    kc[36] = exp((g_RT[10] + g_RT[0]) - (g_RT[14] + g_RT[2]));

    /*reaction 38: N2O + O <=> NO + NO */
    kc[37] = exp((g_RT[10] + g_RT[1]) - (g_RT[8] + g_RT[8]));

    /*reaction 39: N2O + O <=> N2 + O2 */
    kc[38] = exp((g_RT[10] + g_RT[1]) - (g_RT[14] + g_RT[4]));

    /*reaction 40: NH + H <=> N + H2 */
    kc[39] = exp((g_RT[11] + g_RT[0]) - (g_RT[12] + g_RT[3]));

    /*reaction 41: NH + O <=> NO + H */
    kc[40] = exp((g_RT[11] + g_RT[1]) - (g_RT[8] + g_RT[0]));

    /*reaction 42: NH + OH <=> N + H2O */
    kc[41] = exp((g_RT[11] + g_RT[2]) - (g_RT[12] + g_RT[6]));

    /*reaction 43: NH + O2 <=> NO + OH */
    kc[42] = exp((g_RT[11] + g_RT[4]) - (g_RT[8] + g_RT[2]));

    /*reaction 44: NH + NO <=> N2O + H */
    kc[43] = exp((g_RT[11] + g_RT[8]) - (g_RT[10] + g_RT[0]));

    /*reaction 45: NH + NO <=> N2O + H */
    kc[44] = exp((g_RT[11] + g_RT[8]) - (g_RT[10] + g_RT[0]));

    /*reaction 46: NH + NO <=> N2 + OH */
    kc[45] = exp((g_RT[11] + g_RT[8]) - (g_RT[14] + g_RT[2]));

    /*reaction 47: NH + NO2 <=> N2O + OH */
    kc[46] = exp((g_RT[11] + g_RT[9]) - (g_RT[10] + g_RT[2]));

    /*reaction 48: N + OH <=> NO + H */
    kc[47] = exp((g_RT[12] + g_RT[2]) - (g_RT[8] + g_RT[0]));

    /*reaction 49: N + O2 <=> NO + O */
    kc[48] = exp((g_RT[12] + g_RT[4]) - (g_RT[8] + g_RT[1]));

    /*reaction 50: N + NO <=> N2 + O */
    kc[49] = exp((g_RT[12] + g_RT[8]) - (g_RT[14] + g_RT[1]));

    /*reaction 51: NNH <=> N2 + H */
    kc[50] = refC * exp((g_RT[13]) - (g_RT[14] + g_RT[0]));

    /*reaction 52: NNH + H <=> N2 + H2 */
    kc[51] = exp((g_RT[13] + g_RT[0]) - (g_RT[14] + g_RT[3]));

    /*reaction 53: NNH + O <=> N2O + H */
    kc[52] = exp((g_RT[13] + g_RT[1]) - (g_RT[10] + g_RT[0]));

    /*reaction 54: NNH + O <=> N2 + OH */
    kc[53] = exp((g_RT[13] + g_RT[1]) - (g_RT[14] + g_RT[2]));

    /*reaction 55: NNH + O <=> NH + NO */
    kc[54] = exp((g_RT[13] + g_RT[1]) - (g_RT[11] + g_RT[8]));

    /*reaction 56: NNH + OH <=> N2 + H2O */
    kc[55] = exp((g_RT[13] + g_RT[2]) - (g_RT[14] + g_RT[6]));

    /*reaction 57: NNH + O2 <=> N2 + HO2 */
    kc[56] = exp((g_RT[13] + g_RT[4]) - (g_RT[14] + g_RT[5]));

    /*reaction 58: NNH + O2 <=> N2 + H + O2 */
    kc[57] = refC * exp((g_RT[13] + g_RT[4]) - (g_RT[14] + g_RT[0] + g_RT[4]));

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void gibbs(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    static double T_old = -1, species_old[15];

    if (T == T_old)
    {
        for (int i = 0; i < 15; i++)
            species[i] = species_old[i];
        return;
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.54736600e+04 / tc[1]
            +2.94668285e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 1: O */
        species[1] =
            +2.91222592e+04 / tc[1]
            +1.11633364e+00
            -3.16826710e+00 * tc[0]
            +1.63965942e-03 * tc[1]
            -1.10717733e-06 * tc[2]
            +5.10672187e-10 * tc[3]
            -1.05632985e-13 * tc[4];
        /*species 2: OH */
        species[2] =
            +3.37165248e+03 / tc[1]
            +4.09579830e+00
            -3.99198424e+00 * tc[0]
            +1.20053327e-03 * tc[1]
            -7.69440055e-07 * tc[2]
            +3.23263588e-10 * tc[3]
            -6.81597510e-14 * tc[4];
        /*species 3: H2 */
        species[3] =
            -9.17924130e+02 / tc[1]
            +1.66130072e+00
            -2.34430290e+00 * tc[0]
            -3.99021240e-03 * tc[1]
            +3.24631950e-06 * tc[2]
            -1.67974725e-09 * tc[3]
            +3.68801445e-13 * tc[4];
        /*species 4: O2 */
        species[4] =
            -1.06394356e+03 / tc[1]
            +1.24780630e-01
            -3.78245636e+00 * tc[0]
            +1.49836707e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774590e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +2.63190983e+02 / tc[1]
            +5.85910600e-01
            -4.30178800e+00 * tc[0]
            +2.37451005e-03 * tc[1]
            -3.52632550e-06 * tc[2]
            +2.02299675e-09 * tc[3]
            -4.64603350e-13 * tc[4];
        /*species 6: H2O */
        species[6] =
            -3.02937260e+04 / tc[1]
            +5.04764421e+00
            -4.19863520e+00 * tc[0]
            +1.01820085e-03 * tc[1]
            -1.08672360e-06 * tc[2]
            +4.57327242e-10 * tc[3]
            -8.85984000e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.76843601e+04 / tc[1]
            +1.04141933e+00
            -4.31515149e+00 * tc[0]
            +4.23695311e-04 * tc[1]
            -2.94007205e-06 * tc[2]
            +1.88969120e-09 * tc[3]
            -4.54475079e-13 * tc[4];
        /*species 8: NO */
        species[8] =
            +9.81823786e+03 / tc[1]
            +1.93798944e+00
            -4.21859896e+00 * tc[0]
            +2.31994062e-03 * tc[1]
            -1.84071748e-06 * tc[2]
            +7.78379589e-10 * tc[3]
            -1.40277437e-13 * tc[4];
        /*species 9: NO2 */
        species[9] =
            +2.87409757e+03 / tc[1]
            -2.36796070e+00
            -3.94403120e+00 * tc[0]
            +7.92714500e-04 * tc[1]
            -2.77630200e-06 * tc[2]
            +1.70628550e-09 * tc[3]
            -3.91752820e-13 * tc[4];
        /*species 10: N2O */
        species[10] =
            +8.76510000e+03 / tc[1]
            -6.96816400e+00
            -2.54305800e+00 * tc[0]
            -4.74609650e-03 * tc[1]
            +1.63212917e-06 * tc[2]
            -5.21987083e-10 * tc[3]
            +9.50913000e-14 * tc[4];
        /*species 11: NH */
        species[11] =
            +4.18942940e+04 / tc[1]
            +1.64458070e+00
            -3.49290840e+00 * tc[0]
            -1.55895985e-04 * tc[1]
            +2.48174733e-07 * tc[2]
            -2.06803683e-10 * tc[3]
            +5.17848350e-14 * tc[4];
        /*species 12: N */
        species[12] =
            +5.60989000e+04 / tc[1]
            -1.66449500e+00
            -2.50307100e+00 * tc[0]
            +1.09000900e-05 * tc[1]
            -9.03421500e-09 * tc[2]
            +4.70630000e-12 * tc[3]
            -1.04995200e-15 * tc[4];
        /*species 13: NNH */
        species[13] =
            +2.83334700e+04 / tc[1]
            -2.89049300e+00
            -3.50134400e+00 * tc[0]
            -1.02679350e-03 * tc[1]
            -1.19506833e-07 * tc[2]
            -4.10112333e-11 * tc[3]
            +4.83558500e-14 * tc[4];
        /*species 14: N2 */
        species[14] =
            -1.04697628e+03 / tc[1]
            +5.63534900e-01
            -3.53100528e+00 * tc[0]
            +6.18304940e-05 * tc[1]
            +8.38332388e-08 * tc[2]
            -2.02942177e-10 * tc[3]
            +7.04406175e-14 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.54736600e+04 / tc[1]
            +2.94668285e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 1: O */
        species[1] =
            +2.92260120e+04 / tc[1]
            -2.37865760e+00
            -2.54363697e+00 * tc[0]
            +1.36581243e-05 * tc[1]
            +6.98382533e-10 * tc[2]
            -4.12901538e-13 * tc[3]
            +2.39776847e-17 * tc[4];
        /*species 2: OH */
        species[2] =
            +3.70056220e+03 / tc[1]
            -3.00660061e+00
            -2.83853033e+00 * tc[0]
            -5.53706445e-04 * tc[1]
            +4.90000348e-08 * tc[2]
            -3.50582274e-12 * tc[3]
            +1.21144945e-16 * tc[4];
        /*species 3: H2 */
        species[3] =
            -8.13055820e+02 / tc[1]
            +3.95714690e+00
            -2.93283050e+00 * tc[0]
            -4.13299010e-04 * tc[1]
            +2.44000950e-08 * tc[2]
            -1.28415425e-12 * tc[3]
            +3.44398075e-17 * tc[4];
        /*species 4: O2 */
        species[4] =
            -1.21597725e+03 / tc[1]
            +2.45598990e-01
            -3.66096083e+00 * tc[0]
            -3.28182761e-04 * tc[1]
            +2.35249142e-08 * tc[2]
            -1.71498048e-12 * tc[3]
            +6.49566240e-17 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +3.02010736e+01 / tc[1]
            +1.21529210e+00
            -4.17226590e+00 * tc[0]
            -9.40604900e-04 * tc[1]
            +5.77154950e-08 * tc[2]
            -1.62237633e-12 * tc[3]
            -8.80457650e-18 * tc[4];
        /*species 6: H2O */
        species[6] =
            -2.98858940e+04 / tc[1]
            -4.20551110e+00
            -2.67703890e+00 * tc[0]
            -1.48659080e-03 * tc[1]
            +1.28961482e-07 * tc[2]
            -7.86945950e-12 * tc[3]
            +2.13449955e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.79847939e+04 / tc[1]
            +3.91480339e+00
            -4.57977305e+00 * tc[0]
            -2.02663002e-03 * tc[1]
            +2.16407883e-07 * tc[2]
            -1.65176167e-11 * tc[3]
            +5.69843960e-16 * tc[4];
        /*species 8: NO */
        species[8] =
            +9.89456954e+03 / tc[1]
            -3.10829235e+00
            -3.26071234e+00 * tc[0]
            -5.95505675e-04 * tc[1]
            +7.15204410e-08 * tc[2]
            -5.78734552e-12 * tc[3]
            +2.01647840e-16 * tc[4];
        /*species 9: NO2 */
        species[9] =
            +2.29397777e+03 / tc[1]
            +5.00217095e+00
            -4.88475400e+00 * tc[0]
            -1.08619775e-03 * tc[1]
            +1.38011515e-07 * tc[2]
            -1.31229250e-11 * tc[3]
            +5.25544750e-16 * tc[4];
        /*species 10: N2O */
        species[10] =
            +8.16581100e+03 / tc[1]
            +6.37622700e+00
            -4.71897700e+00 * tc[0]
            -1.43685700e-03 * tc[1]
            +1.99582667e-07 * tc[2]
            -1.87546000e-11 * tc[3]
            +7.87668500e-16 * tc[4];
        /*species 11: NH */
        species[11] =
            +4.21345140e+04 / tc[1]
            -2.95708690e+00
            -2.78369290e+00 * tc[0]
            -6.64921450e-04 * tc[1]
            +7.07967450e-08 * tc[2]
            -6.52904200e-12 * tc[3]
            +2.75222350e-16 * tc[4];
        /*species 12: N */
        species[12] =
            +5.61160400e+04 / tc[1]
            -1.99849000e+00
            -2.45026800e+00 * tc[0]
            -5.33073000e-05 * tc[1]
            +1.24422283e-08 * tc[2]
            -1.56637667e-12 * tc[3]
            +5.12992000e-17 * tc[4];
        /*species 13: NNH */
        species[13] =
            +2.78802900e+04 / tc[1]
            +3.51105320e+00
            -4.41534200e+00 * tc[0]
            -8.07194000e-04 * tc[1]
            +2.72149000e-08 * tc[2]
            +7.13320500e-12 * tc[3]
            -8.07395500e-16 * tc[4];
        /*species 14: N2 */
        species[14] =
            -9.23948688e+02 / tc[1]
            -2.91931125e+00
            -2.95257637e+00 * tc[0]
            -6.98450200e-04 * tc[1]
            +8.21052672e-08 * tc[2]
            -6.55008496e-12 * tc[3]
            +2.30377602e-16 * tc[4];
    }
    T_old = T;
    for (int i = 0; i < 15; i++)
        species_old[i] = species[i];
    return;
}


/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void helmholtz(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.54736600e+04 / tc[1]
            +1.94668285e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 1: O */
        species[1] =
            +2.91222592e+04 / tc[1]
            +1.16333640e-01
            -3.16826710e+00 * tc[0]
            +1.63965942e-03 * tc[1]
            -1.10717733e-06 * tc[2]
            +5.10672187e-10 * tc[3]
            -1.05632985e-13 * tc[4];
        /*species 2: OH */
        species[2] =
            +3.37165248e+03 / tc[1]
            +3.09579830e+00
            -3.99198424e+00 * tc[0]
            +1.20053327e-03 * tc[1]
            -7.69440055e-07 * tc[2]
            +3.23263588e-10 * tc[3]
            -6.81597510e-14 * tc[4];
        /*species 3: H2 */
        species[3] =
            -9.17924130e+02 / tc[1]
            +6.61300720e-01
            -2.34430290e+00 * tc[0]
            -3.99021240e-03 * tc[1]
            +3.24631950e-06 * tc[2]
            -1.67974725e-09 * tc[3]
            +3.68801445e-13 * tc[4];
        /*species 4: O2 */
        species[4] =
            -1.06394356e+03 / tc[1]
            -8.75219370e-01
            -3.78245636e+00 * tc[0]
            +1.49836707e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774590e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +2.63190983e+02 / tc[1]
            -4.14089400e-01
            -4.30178800e+00 * tc[0]
            +2.37451005e-03 * tc[1]
            -3.52632550e-06 * tc[2]
            +2.02299675e-09 * tc[3]
            -4.64603350e-13 * tc[4];
        /*species 6: H2O */
        species[6] =
            -3.02937260e+04 / tc[1]
            +4.04764421e+00
            -4.19863520e+00 * tc[0]
            +1.01820085e-03 * tc[1]
            -1.08672360e-06 * tc[2]
            +4.57327242e-10 * tc[3]
            -8.85984000e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.76843601e+04 / tc[1]
            +4.14193300e-02
            -4.31515149e+00 * tc[0]
            +4.23695311e-04 * tc[1]
            -2.94007205e-06 * tc[2]
            +1.88969120e-09 * tc[3]
            -4.54475079e-13 * tc[4];
        /*species 8: NO */
        species[8] =
            +9.81823786e+03 / tc[1]
            +9.37989440e-01
            -4.21859896e+00 * tc[0]
            +2.31994062e-03 * tc[1]
            -1.84071748e-06 * tc[2]
            +7.78379589e-10 * tc[3]
            -1.40277437e-13 * tc[4];
        /*species 9: NO2 */
        species[9] =
            +2.87409757e+03 / tc[1]
            -3.36796070e+00
            -3.94403120e+00 * tc[0]
            +7.92714500e-04 * tc[1]
            -2.77630200e-06 * tc[2]
            +1.70628550e-09 * tc[3]
            -3.91752820e-13 * tc[4];
        /*species 10: N2O */
        species[10] =
            +8.76510000e+03 / tc[1]
            -7.96816400e+00
            -2.54305800e+00 * tc[0]
            -4.74609650e-03 * tc[1]
            +1.63212917e-06 * tc[2]
            -5.21987083e-10 * tc[3]
            +9.50913000e-14 * tc[4];
        /*species 11: NH */
        species[11] =
            +4.18942940e+04 / tc[1]
            +6.44580700e-01
            -3.49290840e+00 * tc[0]
            -1.55895985e-04 * tc[1]
            +2.48174733e-07 * tc[2]
            -2.06803683e-10 * tc[3]
            +5.17848350e-14 * tc[4];
        /*species 12: N */
        species[12] =
            +5.60989000e+04 / tc[1]
            -2.66449500e+00
            -2.50307100e+00 * tc[0]
            +1.09000900e-05 * tc[1]
            -9.03421500e-09 * tc[2]
            +4.70630000e-12 * tc[3]
            -1.04995200e-15 * tc[4];
        /*species 13: NNH */
        species[13] =
            +2.83334700e+04 / tc[1]
            -3.89049300e+00
            -3.50134400e+00 * tc[0]
            -1.02679350e-03 * tc[1]
            -1.19506833e-07 * tc[2]
            -4.10112333e-11 * tc[3]
            +4.83558500e-14 * tc[4];
        /*species 14: N2 */
        species[14] =
            -1.04697628e+03 / tc[1]
            -4.36465100e-01
            -3.53100528e+00 * tc[0]
            +6.18304940e-05 * tc[1]
            +8.38332388e-08 * tc[2]
            -2.02942177e-10 * tc[3]
            +7.04406175e-14 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.54736600e+04 / tc[1]
            +1.94668285e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 1: O */
        species[1] =
            +2.92260120e+04 / tc[1]
            -3.37865760e+00
            -2.54363697e+00 * tc[0]
            +1.36581243e-05 * tc[1]
            +6.98382533e-10 * tc[2]
            -4.12901538e-13 * tc[3]
            +2.39776847e-17 * tc[4];
        /*species 2: OH */
        species[2] =
            +3.70056220e+03 / tc[1]
            -4.00660061e+00
            -2.83853033e+00 * tc[0]
            -5.53706445e-04 * tc[1]
            +4.90000348e-08 * tc[2]
            -3.50582274e-12 * tc[3]
            +1.21144945e-16 * tc[4];
        /*species 3: H2 */
        species[3] =
            -8.13055820e+02 / tc[1]
            +2.95714690e+00
            -2.93283050e+00 * tc[0]
            -4.13299010e-04 * tc[1]
            +2.44000950e-08 * tc[2]
            -1.28415425e-12 * tc[3]
            +3.44398075e-17 * tc[4];
        /*species 4: O2 */
        species[4] =
            -1.21597725e+03 / tc[1]
            -7.54401010e-01
            -3.66096083e+00 * tc[0]
            -3.28182761e-04 * tc[1]
            +2.35249142e-08 * tc[2]
            -1.71498048e-12 * tc[3]
            +6.49566240e-17 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +3.02010736e+01 / tc[1]
            +2.15292100e-01
            -4.17226590e+00 * tc[0]
            -9.40604900e-04 * tc[1]
            +5.77154950e-08 * tc[2]
            -1.62237633e-12 * tc[3]
            -8.80457650e-18 * tc[4];
        /*species 6: H2O */
        species[6] =
            -2.98858940e+04 / tc[1]
            -5.20551110e+00
            -2.67703890e+00 * tc[0]
            -1.48659080e-03 * tc[1]
            +1.28961482e-07 * tc[2]
            -7.86945950e-12 * tc[3]
            +2.13449955e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.79847939e+04 / tc[1]
            +2.91480339e+00
            -4.57977305e+00 * tc[0]
            -2.02663002e-03 * tc[1]
            +2.16407883e-07 * tc[2]
            -1.65176167e-11 * tc[3]
            +5.69843960e-16 * tc[4];
        /*species 8: NO */
        species[8] =
            +9.89456954e+03 / tc[1]
            -4.10829235e+00
            -3.26071234e+00 * tc[0]
            -5.95505675e-04 * tc[1]
            +7.15204410e-08 * tc[2]
            -5.78734552e-12 * tc[3]
            +2.01647840e-16 * tc[4];
        /*species 9: NO2 */
        species[9] =
            +2.29397777e+03 / tc[1]
            +4.00217095e+00
            -4.88475400e+00 * tc[0]
            -1.08619775e-03 * tc[1]
            +1.38011515e-07 * tc[2]
            -1.31229250e-11 * tc[3]
            +5.25544750e-16 * tc[4];
        /*species 10: N2O */
        species[10] =
            +8.16581100e+03 / tc[1]
            +5.37622700e+00
            -4.71897700e+00 * tc[0]
            -1.43685700e-03 * tc[1]
            +1.99582667e-07 * tc[2]
            -1.87546000e-11 * tc[3]
            +7.87668500e-16 * tc[4];
        /*species 11: NH */
        species[11] =
            +4.21345140e+04 / tc[1]
            -3.95708690e+00
            -2.78369290e+00 * tc[0]
            -6.64921450e-04 * tc[1]
            +7.07967450e-08 * tc[2]
            -6.52904200e-12 * tc[3]
            +2.75222350e-16 * tc[4];
        /*species 12: N */
        species[12] =
            +5.61160400e+04 / tc[1]
            -2.99849000e+00
            -2.45026800e+00 * tc[0]
            -5.33073000e-05 * tc[1]
            +1.24422283e-08 * tc[2]
            -1.56637667e-12 * tc[3]
            +5.12992000e-17 * tc[4];
        /*species 13: NNH */
        species[13] =
            +2.78802900e+04 / tc[1]
            +2.51105320e+00
            -4.41534200e+00 * tc[0]
            -8.07194000e-04 * tc[1]
            +2.72149000e-08 * tc[2]
            +7.13320500e-12 * tc[3]
            -8.07395500e-16 * tc[4];
        /*species 14: N2 */
        species[14] =
            -9.23948688e+02 / tc[1]
            -3.91931125e+00
            -2.95257637e+00 * tc[0]
            -6.98450200e-04 * tc[1]
            +8.21052672e-08 * tc[2]
            -6.55008496e-12 * tc[3]
            +2.30377602e-16 * tc[4];
    }
    return;
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cv_R(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: O */
        species[1] =
            +2.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 2: OH */
        species[2] =
            +2.99198424e+00
            -2.40106655e-03 * tc[1]
            +4.61664033e-06 * tc[2]
            -3.87916306e-09 * tc[3]
            +1.36319502e-12 * tc[4];
        /*species 3: H2 */
        species[3] =
            +1.34430290e+00
            +7.98042480e-03 * tc[1]
            -1.94779170e-05 * tc[2]
            +2.01569670e-08 * tc[3]
            -7.37602890e-12 * tc[4];
        /*species 4: O2 */
        species[4] =
            +2.78245636e+00
            -2.99673415e-03 * tc[1]
            +9.84730200e-06 * tc[2]
            -9.68129508e-09 * tc[3]
            +3.24372836e-12 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +3.30178800e+00
            -4.74902010e-03 * tc[1]
            +2.11579530e-05 * tc[2]
            -2.42759610e-08 * tc[3]
            +9.29206700e-12 * tc[4];
        /*species 6: H2O */
        species[6] =
            +3.19863520e+00
            -2.03640170e-03 * tc[1]
            +6.52034160e-06 * tc[2]
            -5.48792690e-09 * tc[3]
            +1.77196800e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.31515149e+00
            -8.47390622e-04 * tc[1]
            +1.76404323e-05 * tc[2]
            -2.26762944e-08 * tc[3]
            +9.08950158e-12 * tc[4];
        /*species 8: NO */
        species[8] =
            +3.21859896e+00
            -4.63988124e-03 * tc[1]
            +1.10443049e-05 * tc[2]
            -9.34055507e-09 * tc[3]
            +2.80554874e-12 * tc[4];
        /*species 9: NO2 */
        species[9] =
            +2.94403120e+00
            -1.58542900e-03 * tc[1]
            +1.66578120e-05 * tc[2]
            -2.04754260e-08 * tc[3]
            +7.83505640e-12 * tc[4];
        /*species 10: N2O */
        species[10] =
            +1.54305800e+00
            +9.49219300e-03 * tc[1]
            -9.79277500e-06 * tc[2]
            +6.26384500e-09 * tc[3]
            -1.90182600e-12 * tc[4];
        /*species 11: NH */
        species[11] =
            +2.49290840e+00
            +3.11791970e-04 * tc[1]
            -1.48904840e-06 * tc[2]
            +2.48164420e-09 * tc[3]
            -1.03569670e-12 * tc[4];
        /*species 12: N */
        species[12] =
            +1.50307100e+00
            -2.18001800e-05 * tc[1]
            +5.42052900e-08 * tc[2]
            -5.64756000e-11 * tc[3]
            +2.09990400e-14 * tc[4];
        /*species 13: NNH */
        species[13] =
            +2.50134400e+00
            +2.05358700e-03 * tc[1]
            +7.17041000e-07 * tc[2]
            +4.92134800e-10 * tc[3]
            -9.67117000e-13 * tc[4];
        /*species 14: N2 */
        species[14] =
            +2.53100528e+00
            -1.23660988e-04 * tc[1]
            -5.02999433e-07 * tc[2]
            +2.43530612e-09 * tc[3]
            -1.40881235e-12 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: O */
        species[1] =
            +1.54363697e+00
            -2.73162486e-05 * tc[1]
            -4.19029520e-09 * tc[2]
            +4.95481845e-12 * tc[3]
            -4.79553694e-16 * tc[4];
        /*species 2: OH */
        species[2] =
            +1.83853033e+00
            +1.10741289e-03 * tc[1]
            -2.94000209e-07 * tc[2]
            +4.20698729e-11 * tc[3]
            -2.42289890e-15 * tc[4];
        /*species 3: H2 */
        species[3] =
            +1.93283050e+00
            +8.26598020e-04 * tc[1]
            -1.46400570e-07 * tc[2]
            +1.54098510e-11 * tc[3]
            -6.88796150e-16 * tc[4];
        /*species 4: O2 */
        species[4] =
            +2.66096083e+00
            +6.56365523e-04 * tc[1]
            -1.41149485e-07 * tc[2]
            +2.05797658e-11 * tc[3]
            -1.29913248e-15 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +3.17226590e+00
            +1.88120980e-03 * tc[1]
            -3.46292970e-07 * tc[2]
            +1.94685160e-11 * tc[3]
            +1.76091530e-16 * tc[4];
        /*species 6: H2O */
        species[6] =
            +1.67703890e+00
            +2.97318160e-03 * tc[1]
            -7.73768890e-07 * tc[2]
            +9.44335140e-11 * tc[3]
            -4.26899910e-15 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.57977305e+00
            +4.05326003e-03 * tc[1]
            -1.29844730e-06 * tc[2]
            +1.98211400e-10 * tc[3]
            -1.13968792e-14 * tc[4];
        /*species 8: NO */
        species[8] =
            +2.26071234e+00
            +1.19101135e-03 * tc[1]
            -4.29122646e-07 * tc[2]
            +6.94481463e-11 * tc[3]
            -4.03295681e-15 * tc[4];
        /*species 9: NO2 */
        species[9] =
            +3.88475400e+00
            +2.17239550e-03 * tc[1]
            -8.28069090e-07 * tc[2]
            +1.57475100e-10 * tc[3]
            -1.05108950e-14 * tc[4];
        /*species 10: N2O */
        species[10] =
            +3.71897700e+00
            +2.87371400e-03 * tc[1]
            -1.19749600e-06 * tc[2]
            +2.25055200e-10 * tc[3]
            -1.57533700e-14 * tc[4];
        /*species 11: NH */
        species[11] =
            +1.78369290e+00
            +1.32984290e-03 * tc[1]
            -4.24780470e-07 * tc[2]
            +7.83485040e-11 * tc[3]
            -5.50444700e-15 * tc[4];
        /*species 12: N */
        species[12] =
            +1.45026800e+00
            +1.06614600e-04 * tc[1]
            -7.46533700e-08 * tc[2]
            +1.87965200e-11 * tc[3]
            -1.02598400e-15 * tc[4];
        /*species 13: NNH */
        species[13] =
            +3.41534200e+00
            +1.61438800e-03 * tc[1]
            -1.63289400e-07 * tc[2]
            -8.55984600e-11 * tc[3]
            +1.61479100e-14 * tc[4];
        /*species 14: N2 */
        species[14] =
            +1.95257637e+00
            +1.39690040e-03 * tc[1]
            -4.92631603e-07 * tc[2]
            +7.86010195e-11 * tc[3]
            -4.60755204e-15 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cp_R(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: O */
        species[1] =
            +3.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 2: OH */
        species[2] =
            +3.99198424e+00
            -2.40106655e-03 * tc[1]
            +4.61664033e-06 * tc[2]
            -3.87916306e-09 * tc[3]
            +1.36319502e-12 * tc[4];
        /*species 3: H2 */
        species[3] =
            +2.34430290e+00
            +7.98042480e-03 * tc[1]
            -1.94779170e-05 * tc[2]
            +2.01569670e-08 * tc[3]
            -7.37602890e-12 * tc[4];
        /*species 4: O2 */
        species[4] =
            +3.78245636e+00
            -2.99673415e-03 * tc[1]
            +9.84730200e-06 * tc[2]
            -9.68129508e-09 * tc[3]
            +3.24372836e-12 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +4.30178800e+00
            -4.74902010e-03 * tc[1]
            +2.11579530e-05 * tc[2]
            -2.42759610e-08 * tc[3]
            +9.29206700e-12 * tc[4];
        /*species 6: H2O */
        species[6] =
            +4.19863520e+00
            -2.03640170e-03 * tc[1]
            +6.52034160e-06 * tc[2]
            -5.48792690e-09 * tc[3]
            +1.77196800e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.31515149e+00
            -8.47390622e-04 * tc[1]
            +1.76404323e-05 * tc[2]
            -2.26762944e-08 * tc[3]
            +9.08950158e-12 * tc[4];
        /*species 8: NO */
        species[8] =
            +4.21859896e+00
            -4.63988124e-03 * tc[1]
            +1.10443049e-05 * tc[2]
            -9.34055507e-09 * tc[3]
            +2.80554874e-12 * tc[4];
        /*species 9: NO2 */
        species[9] =
            +3.94403120e+00
            -1.58542900e-03 * tc[1]
            +1.66578120e-05 * tc[2]
            -2.04754260e-08 * tc[3]
            +7.83505640e-12 * tc[4];
        /*species 10: N2O */
        species[10] =
            +2.54305800e+00
            +9.49219300e-03 * tc[1]
            -9.79277500e-06 * tc[2]
            +6.26384500e-09 * tc[3]
            -1.90182600e-12 * tc[4];
        /*species 11: NH */
        species[11] =
            +3.49290840e+00
            +3.11791970e-04 * tc[1]
            -1.48904840e-06 * tc[2]
            +2.48164420e-09 * tc[3]
            -1.03569670e-12 * tc[4];
        /*species 12: N */
        species[12] =
            +2.50307100e+00
            -2.18001800e-05 * tc[1]
            +5.42052900e-08 * tc[2]
            -5.64756000e-11 * tc[3]
            +2.09990400e-14 * tc[4];
        /*species 13: NNH */
        species[13] =
            +3.50134400e+00
            +2.05358700e-03 * tc[1]
            +7.17041000e-07 * tc[2]
            +4.92134800e-10 * tc[3]
            -9.67117000e-13 * tc[4];
        /*species 14: N2 */
        species[14] =
            +3.53100528e+00
            -1.23660988e-04 * tc[1]
            -5.02999433e-07 * tc[2]
            +2.43530612e-09 * tc[3]
            -1.40881235e-12 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: O */
        species[1] =
            +2.54363697e+00
            -2.73162486e-05 * tc[1]
            -4.19029520e-09 * tc[2]
            +4.95481845e-12 * tc[3]
            -4.79553694e-16 * tc[4];
        /*species 2: OH */
        species[2] =
            +2.83853033e+00
            +1.10741289e-03 * tc[1]
            -2.94000209e-07 * tc[2]
            +4.20698729e-11 * tc[3]
            -2.42289890e-15 * tc[4];
        /*species 3: H2 */
        species[3] =
            +2.93283050e+00
            +8.26598020e-04 * tc[1]
            -1.46400570e-07 * tc[2]
            +1.54098510e-11 * tc[3]
            -6.88796150e-16 * tc[4];
        /*species 4: O2 */
        species[4] =
            +3.66096083e+00
            +6.56365523e-04 * tc[1]
            -1.41149485e-07 * tc[2]
            +2.05797658e-11 * tc[3]
            -1.29913248e-15 * tc[4];
        /*species 5: HO2 */
        species[5] =
            +4.17226590e+00
            +1.88120980e-03 * tc[1]
            -3.46292970e-07 * tc[2]
            +1.94685160e-11 * tc[3]
            +1.76091530e-16 * tc[4];
        /*species 6: H2O */
        species[6] =
            +2.67703890e+00
            +2.97318160e-03 * tc[1]
            -7.73768890e-07 * tc[2]
            +9.44335140e-11 * tc[3]
            -4.26899910e-15 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.57977305e+00
            +4.05326003e-03 * tc[1]
            -1.29844730e-06 * tc[2]
            +1.98211400e-10 * tc[3]
            -1.13968792e-14 * tc[4];
        /*species 8: NO */
        species[8] =
            +3.26071234e+00
            +1.19101135e-03 * tc[1]
            -4.29122646e-07 * tc[2]
            +6.94481463e-11 * tc[3]
            -4.03295681e-15 * tc[4];
        /*species 9: NO2 */
        species[9] =
            +4.88475400e+00
            +2.17239550e-03 * tc[1]
            -8.28069090e-07 * tc[2]
            +1.57475100e-10 * tc[3]
            -1.05108950e-14 * tc[4];
        /*species 10: N2O */
        species[10] =
            +4.71897700e+00
            +2.87371400e-03 * tc[1]
            -1.19749600e-06 * tc[2]
            +2.25055200e-10 * tc[3]
            -1.57533700e-14 * tc[4];
        /*species 11: NH */
        species[11] =
            +2.78369290e+00
            +1.32984290e-03 * tc[1]
            -4.24780470e-07 * tc[2]
            +7.83485040e-11 * tc[3]
            -5.50444700e-15 * tc[4];
        /*species 12: N */
        species[12] =
            +2.45026800e+00
            +1.06614600e-04 * tc[1]
            -7.46533700e-08 * tc[2]
            +1.87965200e-11 * tc[3]
            -1.02598400e-15 * tc[4];
        /*species 13: NNH */
        species[13] =
            +4.41534200e+00
            +1.61438800e-03 * tc[1]
            -1.63289400e-07 * tc[2]
            -8.55984600e-11 * tc[3]
            +1.61479100e-14 * tc[4];
        /*species 14: N2 */
        species[14] =
            +2.95257637e+00
            +1.39690040e-03 * tc[1]
            -4.92631603e-07 * tc[2]
            +7.86010195e-11 * tc[3]
            -4.60755204e-15 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesInternalEnergy(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54736600e+04 / tc[1];
        /*species 1: O */
        species[1] =
            +2.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 / tc[1];
        /*species 2: OH */
        species[2] =
            +2.99198424e+00
            -1.20053327e-03 * tc[1]
            +1.53888011e-06 * tc[2]
            -9.69790765e-10 * tc[3]
            +2.72639004e-13 * tc[4]
            +3.37165248e+03 / tc[1];
        /*species 3: H2 */
        species[3] =
            +1.34430290e+00
            +3.99021240e-03 * tc[1]
            -6.49263900e-06 * tc[2]
            +5.03924175e-09 * tc[3]
            -1.47520578e-12 * tc[4]
            -9.17924130e+02 / tc[1];
        /*species 4: O2 */
        species[4] =
            +2.78245636e+00
            -1.49836707e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745672e-13 * tc[4]
            -1.06394356e+03 / tc[1];
        /*species 5: HO2 */
        species[5] =
            +3.30178800e+00
            -2.37451005e-03 * tc[1]
            +7.05265100e-06 * tc[2]
            -6.06899025e-09 * tc[3]
            +1.85841340e-12 * tc[4]
            +2.63190983e+02 / tc[1];
        /*species 6: H2O */
        species[6] =
            +3.19863520e+00
            -1.01820085e-03 * tc[1]
            +2.17344720e-06 * tc[2]
            -1.37198172e-09 * tc[3]
            +3.54393600e-13 * tc[4]
            -3.02937260e+04 / tc[1];
        /*species 7: H2O2 */
        species[7] =
            +3.31515149e+00
            -4.23695311e-04 * tc[1]
            +5.88014410e-06 * tc[2]
            -5.66907360e-09 * tc[3]
            +1.81790032e-12 * tc[4]
            -1.76843601e+04 / tc[1];
        /*species 8: NO */
        species[8] =
            +3.21859896e+00
            -2.31994062e-03 * tc[1]
            +3.68143497e-06 * tc[2]
            -2.33513877e-09 * tc[3]
            +5.61109748e-13 * tc[4]
            +9.81823786e+03 / tc[1];
        /*species 9: NO2 */
        species[9] =
            +2.94403120e+00
            -7.92714500e-04 * tc[1]
            +5.55260400e-06 * tc[2]
            -5.11885650e-09 * tc[3]
            +1.56701128e-12 * tc[4]
            +2.87409757e+03 / tc[1];
        /*species 10: N2O */
        species[10] =
            +1.54305800e+00
            +4.74609650e-03 * tc[1]
            -3.26425833e-06 * tc[2]
            +1.56596125e-09 * tc[3]
            -3.80365200e-13 * tc[4]
            +8.76510000e+03 / tc[1];
        /*species 11: NH */
        species[11] =
            +2.49290840e+00
            +1.55895985e-04 * tc[1]
            -4.96349467e-07 * tc[2]
            +6.20411050e-10 * tc[3]
            -2.07139340e-13 * tc[4]
            +4.18942940e+04 / tc[1];
        /*species 12: N */
        species[12] =
            +1.50307100e+00
            -1.09000900e-05 * tc[1]
            +1.80684300e-08 * tc[2]
            -1.41189000e-11 * tc[3]
            +4.19980800e-15 * tc[4]
            +5.60989000e+04 / tc[1];
        /*species 13: NNH */
        species[13] =
            +2.50134400e+00
            +1.02679350e-03 * tc[1]
            +2.39013667e-07 * tc[2]
            +1.23033700e-10 * tc[3]
            -1.93423400e-13 * tc[4]
            +2.83334700e+04 / tc[1];
        /*species 14: N2 */
        species[14] =
            +2.53100528e+00
            -6.18304940e-05 * tc[1]
            -1.67666478e-07 * tc[2]
            +6.08826530e-10 * tc[3]
            -2.81762470e-13 * tc[4]
            -1.04697628e+03 / tc[1];
    } else {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54736600e+04 / tc[1];
        /*species 1: O */
        species[1] =
            +1.54363697e+00
            -1.36581243e-05 * tc[1]
            -1.39676507e-09 * tc[2]
            +1.23870461e-12 * tc[3]
            -9.59107388e-17 * tc[4]
            +2.92260120e+04 / tc[1];
        /*species 2: OH */
        species[2] =
            +1.83853033e+00
            +5.53706445e-04 * tc[1]
            -9.80000697e-08 * tc[2]
            +1.05174682e-11 * tc[3]
            -4.84579780e-16 * tc[4]
            +3.70056220e+03 / tc[1];
        /*species 3: H2 */
        species[3] =
            +1.93283050e+00
            +4.13299010e-04 * tc[1]
            -4.88001900e-08 * tc[2]
            +3.85246275e-12 * tc[3]
            -1.37759230e-16 * tc[4]
            -8.13055820e+02 / tc[1];
        /*species 4: O2 */
        species[4] =
            +2.66096083e+00
            +3.28182761e-04 * tc[1]
            -4.70498283e-08 * tc[2]
            +5.14494145e-12 * tc[3]
            -2.59826496e-16 * tc[4]
            -1.21597725e+03 / tc[1];
        /*species 5: HO2 */
        species[5] =
            +3.17226590e+00
            +9.40604900e-04 * tc[1]
            -1.15430990e-07 * tc[2]
            +4.86712900e-12 * tc[3]
            +3.52183060e-17 * tc[4]
            +3.02010736e+01 / tc[1];
        /*species 6: H2O */
        species[6] =
            +1.67703890e+00
            +1.48659080e-03 * tc[1]
            -2.57922963e-07 * tc[2]
            +2.36083785e-11 * tc[3]
            -8.53799820e-16 * tc[4]
            -2.98858940e+04 / tc[1];
        /*species 7: H2O2 */
        species[7] =
            +3.57977305e+00
            +2.02663002e-03 * tc[1]
            -4.32815767e-07 * tc[2]
            +4.95528500e-11 * tc[3]
            -2.27937584e-15 * tc[4]
            -1.79847939e+04 / tc[1];
        /*species 8: NO */
        species[8] =
            +2.26071234e+00
            +5.95505675e-04 * tc[1]
            -1.43040882e-07 * tc[2]
            +1.73620366e-11 * tc[3]
            -8.06591362e-16 * tc[4]
            +9.89456954e+03 / tc[1];
        /*species 9: NO2 */
        species[9] =
            +3.88475400e+00
            +1.08619775e-03 * tc[1]
            -2.76023030e-07 * tc[2]
            +3.93687750e-11 * tc[3]
            -2.10217900e-15 * tc[4]
            +2.29397777e+03 / tc[1];
        /*species 10: N2O */
        species[10] =
            +3.71897700e+00
            +1.43685700e-03 * tc[1]
            -3.99165333e-07 * tc[2]
            +5.62638000e-11 * tc[3]
            -3.15067400e-15 * tc[4]
            +8.16581100e+03 / tc[1];
        /*species 11: NH */
        species[11] =
            +1.78369290e+00
            +6.64921450e-04 * tc[1]
            -1.41593490e-07 * tc[2]
            +1.95871260e-11 * tc[3]
            -1.10088940e-15 * tc[4]
            +4.21345140e+04 / tc[1];
        /*species 12: N */
        species[12] =
            +1.45026800e+00
            +5.33073000e-05 * tc[1]
            -2.48844567e-08 * tc[2]
            +4.69913000e-12 * tc[3]
            -2.05196800e-16 * tc[4]
            +5.61160400e+04 / tc[1];
        /*species 13: NNH */
        species[13] =
            +3.41534200e+00
            +8.07194000e-04 * tc[1]
            -5.44298000e-08 * tc[2]
            -2.13996150e-11 * tc[3]
            +3.22958200e-15 * tc[4]
            +2.78802900e+04 / tc[1];
        /*species 14: N2 */
        species[14] =
            +1.95257637e+00
            +6.98450200e-04 * tc[1]
            -1.64210534e-07 * tc[2]
            +1.96502549e-11 * tc[3]
            -9.21510408e-16 * tc[4]
            -9.23948688e+02 / tc[1];
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEnthalpy(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54736600e+04 / tc[1];
        /*species 1: O */
        species[1] =
            +3.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 / tc[1];
        /*species 2: OH */
        species[2] =
            +3.99198424e+00
            -1.20053327e-03 * tc[1]
            +1.53888011e-06 * tc[2]
            -9.69790765e-10 * tc[3]
            +2.72639004e-13 * tc[4]
            +3.37165248e+03 / tc[1];
        /*species 3: H2 */
        species[3] =
            +2.34430290e+00
            +3.99021240e-03 * tc[1]
            -6.49263900e-06 * tc[2]
            +5.03924175e-09 * tc[3]
            -1.47520578e-12 * tc[4]
            -9.17924130e+02 / tc[1];
        /*species 4: O2 */
        species[4] =
            +3.78245636e+00
            -1.49836707e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745672e-13 * tc[4]
            -1.06394356e+03 / tc[1];
        /*species 5: HO2 */
        species[5] =
            +4.30178800e+00
            -2.37451005e-03 * tc[1]
            +7.05265100e-06 * tc[2]
            -6.06899025e-09 * tc[3]
            +1.85841340e-12 * tc[4]
            +2.63190983e+02 / tc[1];
        /*species 6: H2O */
        species[6] =
            +4.19863520e+00
            -1.01820085e-03 * tc[1]
            +2.17344720e-06 * tc[2]
            -1.37198172e-09 * tc[3]
            +3.54393600e-13 * tc[4]
            -3.02937260e+04 / tc[1];
        /*species 7: H2O2 */
        species[7] =
            +4.31515149e+00
            -4.23695311e-04 * tc[1]
            +5.88014410e-06 * tc[2]
            -5.66907360e-09 * tc[3]
            +1.81790032e-12 * tc[4]
            -1.76843601e+04 / tc[1];
        /*species 8: NO */
        species[8] =
            +4.21859896e+00
            -2.31994062e-03 * tc[1]
            +3.68143497e-06 * tc[2]
            -2.33513877e-09 * tc[3]
            +5.61109748e-13 * tc[4]
            +9.81823786e+03 / tc[1];
        /*species 9: NO2 */
        species[9] =
            +3.94403120e+00
            -7.92714500e-04 * tc[1]
            +5.55260400e-06 * tc[2]
            -5.11885650e-09 * tc[3]
            +1.56701128e-12 * tc[4]
            +2.87409757e+03 / tc[1];
        /*species 10: N2O */
        species[10] =
            +2.54305800e+00
            +4.74609650e-03 * tc[1]
            -3.26425833e-06 * tc[2]
            +1.56596125e-09 * tc[3]
            -3.80365200e-13 * tc[4]
            +8.76510000e+03 / tc[1];
        /*species 11: NH */
        species[11] =
            +3.49290840e+00
            +1.55895985e-04 * tc[1]
            -4.96349467e-07 * tc[2]
            +6.20411050e-10 * tc[3]
            -2.07139340e-13 * tc[4]
            +4.18942940e+04 / tc[1];
        /*species 12: N */
        species[12] =
            +2.50307100e+00
            -1.09000900e-05 * tc[1]
            +1.80684300e-08 * tc[2]
            -1.41189000e-11 * tc[3]
            +4.19980800e-15 * tc[4]
            +5.60989000e+04 / tc[1];
        /*species 13: NNH */
        species[13] =
            +3.50134400e+00
            +1.02679350e-03 * tc[1]
            +2.39013667e-07 * tc[2]
            +1.23033700e-10 * tc[3]
            -1.93423400e-13 * tc[4]
            +2.83334700e+04 / tc[1];
        /*species 14: N2 */
        species[14] =
            +3.53100528e+00
            -6.18304940e-05 * tc[1]
            -1.67666478e-07 * tc[2]
            +6.08826530e-10 * tc[3]
            -2.81762470e-13 * tc[4]
            -1.04697628e+03 / tc[1];
    } else {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54736600e+04 / tc[1];
        /*species 1: O */
        species[1] =
            +2.54363697e+00
            -1.36581243e-05 * tc[1]
            -1.39676507e-09 * tc[2]
            +1.23870461e-12 * tc[3]
            -9.59107388e-17 * tc[4]
            +2.92260120e+04 / tc[1];
        /*species 2: OH */
        species[2] =
            +2.83853033e+00
            +5.53706445e-04 * tc[1]
            -9.80000697e-08 * tc[2]
            +1.05174682e-11 * tc[3]
            -4.84579780e-16 * tc[4]
            +3.70056220e+03 / tc[1];
        /*species 3: H2 */
        species[3] =
            +2.93283050e+00
            +4.13299010e-04 * tc[1]
            -4.88001900e-08 * tc[2]
            +3.85246275e-12 * tc[3]
            -1.37759230e-16 * tc[4]
            -8.13055820e+02 / tc[1];
        /*species 4: O2 */
        species[4] =
            +3.66096083e+00
            +3.28182761e-04 * tc[1]
            -4.70498283e-08 * tc[2]
            +5.14494145e-12 * tc[3]
            -2.59826496e-16 * tc[4]
            -1.21597725e+03 / tc[1];
        /*species 5: HO2 */
        species[5] =
            +4.17226590e+00
            +9.40604900e-04 * tc[1]
            -1.15430990e-07 * tc[2]
            +4.86712900e-12 * tc[3]
            +3.52183060e-17 * tc[4]
            +3.02010736e+01 / tc[1];
        /*species 6: H2O */
        species[6] =
            +2.67703890e+00
            +1.48659080e-03 * tc[1]
            -2.57922963e-07 * tc[2]
            +2.36083785e-11 * tc[3]
            -8.53799820e-16 * tc[4]
            -2.98858940e+04 / tc[1];
        /*species 7: H2O2 */
        species[7] =
            +4.57977305e+00
            +2.02663002e-03 * tc[1]
            -4.32815767e-07 * tc[2]
            +4.95528500e-11 * tc[3]
            -2.27937584e-15 * tc[4]
            -1.79847939e+04 / tc[1];
        /*species 8: NO */
        species[8] =
            +3.26071234e+00
            +5.95505675e-04 * tc[1]
            -1.43040882e-07 * tc[2]
            +1.73620366e-11 * tc[3]
            -8.06591362e-16 * tc[4]
            +9.89456954e+03 / tc[1];
        /*species 9: NO2 */
        species[9] =
            +4.88475400e+00
            +1.08619775e-03 * tc[1]
            -2.76023030e-07 * tc[2]
            +3.93687750e-11 * tc[3]
            -2.10217900e-15 * tc[4]
            +2.29397777e+03 / tc[1];
        /*species 10: N2O */
        species[10] =
            +4.71897700e+00
            +1.43685700e-03 * tc[1]
            -3.99165333e-07 * tc[2]
            +5.62638000e-11 * tc[3]
            -3.15067400e-15 * tc[4]
            +8.16581100e+03 / tc[1];
        /*species 11: NH */
        species[11] =
            +2.78369290e+00
            +6.64921450e-04 * tc[1]
            -1.41593490e-07 * tc[2]
            +1.95871260e-11 * tc[3]
            -1.10088940e-15 * tc[4]
            +4.21345140e+04 / tc[1];
        /*species 12: N */
        species[12] =
            +2.45026800e+00
            +5.33073000e-05 * tc[1]
            -2.48844567e-08 * tc[2]
            +4.69913000e-12 * tc[3]
            -2.05196800e-16 * tc[4]
            +5.61160400e+04 / tc[1];
        /*species 13: NNH */
        species[13] =
            +4.41534200e+00
            +8.07194000e-04 * tc[1]
            -5.44298000e-08 * tc[2]
            -2.13996150e-11 * tc[3]
            +3.22958200e-15 * tc[4]
            +2.78802900e+04 / tc[1];
        /*species 14: N2 */
        species[14] =
            +2.95257637e+00
            +6.98450200e-04 * tc[1]
            -1.64210534e-07 * tc[2]
            +1.96502549e-11 * tc[3]
            -9.21510408e-16 * tc[4]
            -9.23948688e+02 / tc[1];
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEntropy(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.46682850e-01 ;
        /*species 1: O */
        species[1] =
            +3.16826710e+00 * tc[0]
            -3.27931884e-03 * tc[1]
            +3.32153198e-06 * tc[2]
            -2.04268875e-09 * tc[3]
            +5.28164927e-13 * tc[4]
            +2.05193346e+00 ;
        /*species 2: OH */
        species[2] =
            +3.99198424e+00 * tc[0]
            -2.40106655e-03 * tc[1]
            +2.30832017e-06 * tc[2]
            -1.29305435e-09 * tc[3]
            +3.40798755e-13 * tc[4]
            -1.03814059e-01 ;
        /*species 3: H2 */
        species[3] =
            +2.34430290e+00 * tc[0]
            +7.98042480e-03 * tc[1]
            -9.73895850e-06 * tc[2]
            +6.71898900e-09 * tc[3]
            -1.84400722e-12 * tc[4]
            +6.83002180e-01 ;
        /*species 4: O2 */
        species[4] =
            +3.78245636e+00 * tc[0]
            -2.99673415e-03 * tc[1]
            +4.92365100e-06 * tc[2]
            -3.22709836e-09 * tc[3]
            +8.10932090e-13 * tc[4]
            +3.65767573e+00 ;
        /*species 5: HO2 */
        species[5] =
            +4.30178800e+00 * tc[0]
            -4.74902010e-03 * tc[1]
            +1.05789765e-05 * tc[2]
            -8.09198700e-09 * tc[3]
            +2.32301675e-12 * tc[4]
            +3.71587740e+00 ;
        /*species 6: H2O */
        species[6] =
            +4.19863520e+00 * tc[0]
            -2.03640170e-03 * tc[1]
            +3.26017080e-06 * tc[2]
            -1.82930897e-09 * tc[3]
            +4.42992000e-13 * tc[4]
            -8.49009010e-01 ;
        /*species 7: H2O2 */
        species[7] =
            +4.31515149e+00 * tc[0]
            -8.47390622e-04 * tc[1]
            +8.82021615e-06 * tc[2]
            -7.55876480e-09 * tc[3]
            +2.27237539e-12 * tc[4]
            +3.27373216e+00 ;
        /*species 8: NO */
        species[8] =
            +4.21859896e+00 * tc[0]
            -4.63988124e-03 * tc[1]
            +5.52215245e-06 * tc[2]
            -3.11351836e-09 * tc[3]
            +7.01387185e-13 * tc[4]
            +2.28060952e+00 ;
        /*species 9: NO2 */
        species[9] =
            +3.94403120e+00 * tc[0]
            -1.58542900e-03 * tc[1]
            +8.32890600e-06 * tc[2]
            -6.82514200e-09 * tc[3]
            +1.95876410e-12 * tc[4]
            +6.31199190e+00 ;
        /*species 10: N2O */
        species[10] =
            +2.54305800e+00 * tc[0]
            +9.49219300e-03 * tc[1]
            -4.89638750e-06 * tc[2]
            +2.08794833e-09 * tc[3]
            -4.75456500e-13 * tc[4]
            +9.51122200e+00 ;
        /*species 11: NH */
        species[11] =
            +3.49290840e+00 * tc[0]
            +3.11791970e-04 * tc[1]
            -7.44524200e-07 * tc[2]
            +8.27214733e-10 * tc[3]
            -2.58924175e-13 * tc[4]
            +1.84832770e+00 ;
        /*species 12: N */
        species[12] =
            +2.50307100e+00 * tc[0]
            -2.18001800e-05 * tc[1]
            +2.71026450e-08 * tc[2]
            -1.88252000e-11 * tc[3]
            +5.24976000e-15 * tc[4]
            +4.16756600e+00 ;
        /*species 13: NNH */
        species[13] =
            +3.50134400e+00 * tc[0]
            +2.05358700e-03 * tc[1]
            +3.58520500e-07 * tc[2]
            +1.64044933e-10 * tc[3]
            -2.41779250e-13 * tc[4]
            +6.39183700e+00 ;
        /*species 14: N2 */
        species[14] =
            +3.53100528e+00 * tc[0]
            -1.23660988e-04 * tc[1]
            -2.51499717e-07 * tc[2]
            +8.11768707e-10 * tc[3]
            -3.52203088e-13 * tc[4]
            +2.96747038e+00 ;
    } else {
        /*species 0: H */
        species[0] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.46682850e-01 ;
        /*species 1: O */
        species[1] =
            +2.54363697e+00 * tc[0]
            -2.73162486e-05 * tc[1]
            -2.09514760e-09 * tc[2]
            +1.65160615e-12 * tc[3]
            -1.19888423e-16 * tc[4]
            +4.92229457e+00 ;
        /*species 2: OH */
        species[2] =
            +2.83853033e+00 * tc[0]
            +1.10741289e-03 * tc[1]
            -1.47000104e-07 * tc[2]
            +1.40232910e-11 * tc[3]
            -6.05724725e-16 * tc[4]
            +5.84513094e+00 ;
        /*species 3: H2 */
        species[3] =
            +2.93283050e+00 * tc[0]
            +8.26598020e-04 * tc[1]
            -7.32002850e-08 * tc[2]
            +5.13661700e-12 * tc[3]
            -1.72199037e-16 * tc[4]
            -1.02431640e+00 ;
        /*species 4: O2 */
        species[4] =
            +3.66096083e+00 * tc[0]
            +6.56365523e-04 * tc[1]
            -7.05747425e-08 * tc[2]
            +6.85992193e-12 * tc[3]
            -3.24783120e-16 * tc[4]
            +3.41536184e+00 ;
        /*species 5: HO2 */
        species[5] =
            +4.17226590e+00 * tc[0]
            +1.88120980e-03 * tc[1]
            -1.73146485e-07 * tc[2]
            +6.48950533e-12 * tc[3]
            +4.40228825e-17 * tc[4]
            +2.95697380e+00 ;
        /*species 6: H2O */
        species[6] =
            +2.67703890e+00 * tc[0]
            +2.97318160e-03 * tc[1]
            -3.86884445e-07 * tc[2]
            +3.14778380e-11 * tc[3]
            -1.06724977e-15 * tc[4]
            +6.88255000e+00 ;
        /*species 7: H2O2 */
        species[7] =
            +4.57977305e+00 * tc[0]
            +4.05326003e-03 * tc[1]
            -6.49223650e-07 * tc[2]
            +6.60704667e-11 * tc[3]
            -2.84921980e-15 * tc[4]
            +6.64969660e-01 ;
        /*species 8: NO */
        species[8] =
            +3.26071234e+00 * tc[0]
            +1.19101135e-03 * tc[1]
            -2.14561323e-07 * tc[2]
            +2.31493821e-11 * tc[3]
            -1.00823920e-15 * tc[4]
            +6.36900469e+00 ;
        /*species 9: NO2 */
        species[9] =
            +4.88475400e+00 * tc[0]
            +2.17239550e-03 * tc[1]
            -4.14034545e-07 * tc[2]
            +5.24917000e-11 * tc[3]
            -2.62772375e-15 * tc[4]
            -1.17416951e-01 ;
        /*species 10: N2O */
        species[10] =
            +4.71897700e+00 * tc[0]
            +2.87371400e-03 * tc[1]
            -5.98748000e-07 * tc[2]
            +7.50184000e-11 * tc[3]
            -3.93834250e-15 * tc[4]
            -1.65725000e+00 ;
        /*species 11: NH */
        species[11] =
            +2.78369290e+00 * tc[0]
            +1.32984290e-03 * tc[1]
            -2.12390235e-07 * tc[2]
            +2.61161680e-11 * tc[3]
            -1.37611175e-15 * tc[4]
            +5.74077980e+00 ;
        /*species 12: N */
        species[12] =
            +2.45026800e+00 * tc[0]
            +1.06614600e-04 * tc[1]
            -3.73266850e-08 * tc[2]
            +6.26550667e-12 * tc[3]
            -2.56496000e-16 * tc[4]
            +4.44875800e+00 ;
        /*species 13: NNH */
        species[13] =
            +4.41534200e+00 * tc[0]
            +1.61438800e-03 * tc[1]
            -8.16447000e-08 * tc[2]
            -2.85328200e-11 * tc[3]
            +4.03697750e-15 * tc[4]
            +9.04288800e-01 ;
        /*species 14: N2 */
        species[14] =
            +2.95257637e+00 * tc[0]
            +1.39690040e-03 * tc[1]
            -2.46315802e-07 * tc[2]
            +2.62003398e-11 * tc[3]
            -1.15188801e-15 * tc[4]
            +5.87188762e+00 ;
    }
    return;
}


/*save molecular weights into array */
void molecularWeight(double * wt)
{
    wt[0] = 1.007970; /*H */
    wt[1] = 15.999400; /*O */
    wt[2] = 17.007370; /*OH */
    wt[3] = 2.015940; /*H2 */
    wt[4] = 31.998800; /*O2 */
    wt[5] = 33.006770; /*HO2 */
    wt[6] = 18.015340; /*H2O */
    wt[7] = 34.014740; /*H2O2 */
    wt[8] = 30.006100; /*NO */
    wt[9] = 46.005500; /*NO2 */
    wt[10] = 44.012800; /*N2O */
    wt[11] = 15.014670; /*NH */
    wt[12] = 14.006700; /*N */
    wt[13] = 29.021370; /*NNH */
    wt[14] = 28.013400; /*N2 */

    return;
}


/*get temperature given internal energy in mass units and mass fracs */
int feeytt_(double * e, double * y, int * iwrk, double * rwrk, double * t)
{
    const int maxiter = 50;
    const double tol  = 0.001;
    double ein  = *e;
    double tmin = 300; // max lower bound for thermo def
    double tmax = 4000; // min upper bound for thermo def
    double e1,emin,emax,cv,t1,dt;
    int i; // loop counter
    CKUBMS(&tmin, y, iwrk, rwrk, &emin);
    CKUBMS(&tmax, y, iwrk, rwrk, &emax);
    if (ein < emin) {
        /*Linear Extrapolation below tmin */
        CKCVBS(&tmin, y, iwrk, rwrk, &cv);
        *t = tmin - (emin-ein)/cv;
        return 1;
    }
    if (ein > emax) {
        /*Linear Extrapolation above tmax */
        CKCVBS(&tmax, y, iwrk, rwrk, &cv);
        *t = tmax - (emax-ein)/cv;
        return 1;
    }
    t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);
    for (i = 0; i < maxiter; ++i) {
        CKUBMS(&t1,y,iwrk,rwrk,&e1);
        CKCVBS(&t1,y,iwrk,rwrk,&cv);
        dt = (ein - e1) / cv;
        if (dt > 100) { dt = 100; }
        else if (dt < -100) { dt = -100; }
        else if (fabs(dt) < tol) break;
        t1 += dt;
    }
    *t = t1;
    return 0;
}


/*convert phi[species] (specific mole nums) to y[species] (mass fracs) */
void fephity_(double * phi, int * iwrk, double * rwrk, double * y)
{
    double XW  = 0; 
    int id; /*loop counter */
    /*Compute mean molecular wt first */
    y[0] = phi[0]*1.007970;   XW += y[0]; /*H */
    y[1] = phi[1]*15.999400;   XW += y[1]; /*O */
    y[2] = phi[2]*17.007370;   XW += y[2]; /*OH */
    y[3] = phi[3]*2.015940;   XW += y[3]; /*H2 */
    y[4] = phi[4]*31.998800;   XW += y[4]; /*O2 */
    y[5] = phi[5]*33.006770;   XW += y[5]; /*HO2 */
    y[6] = phi[6]*18.015340;   XW += y[6]; /*H2O */
    y[7] = phi[7]*34.014740;   XW += y[7]; /*H2O2 */
    y[8] = phi[8]*30.006100;   XW += y[8]; /*NO */
    y[9] = phi[9]*46.005500;   XW += y[9]; /*NO2 */
    y[10] = phi[10]*44.012800;   XW += y[10]; /*N2O */
    y[11] = phi[11]*15.014670;   XW += y[11]; /*NH */
    y[12] = phi[12]*14.006700;   XW += y[12]; /*N */
    y[13] = phi[13]*29.021370;   XW += y[13]; /*NNH */
    y[14] = phi[14]*28.013400;   XW += y[14]; /*N2 */
    for (id = 0; id < 15; ++id) {
        y[id] = y[id]/XW;
    }

    return;
}


/*convert y[species] (mass fracs) to phi[species] (specific mole num) */
void feytphi_(double * y, int * iwrk, double * rwrk, double * phi)
{
    phi[0] = y[0]/ 1.00797000e-03; /*H (wt in kg) */
    phi[1] = y[1]/ 1.59994000e-02; /*O (wt in kg) */
    phi[2] = y[2]/ 1.70073700e-02; /*OH (wt in kg) */
    phi[3] = y[3]/ 2.01594000e-03; /*H2 (wt in kg) */
    phi[4] = y[4]/ 3.19988000e-02; /*O2 (wt in kg) */
    phi[5] = y[5]/ 3.30067700e-02; /*HO2 (wt in kg) */
    phi[6] = y[6]/ 1.80153400e-02; /*H2O (wt in kg) */
    phi[7] = y[7]/ 3.40147400e-02; /*H2O2 (wt in kg) */
    phi[8] = y[8]/ 3.00061000e-02; /*NO (wt in kg) */
    phi[9] = y[9]/ 4.60055000e-02; /*NO2 (wt in kg) */
    phi[10] = y[10]/ 4.40128000e-02; /*N2O (wt in kg) */
    phi[11] = y[11]/ 1.50146700e-02; /*NH (wt in kg) */
    phi[12] = y[12]/ 1.40067000e-02; /*N (wt in kg) */
    phi[13] = y[13]/ 2.90213700e-02; /*NNH (wt in kg) */
    phi[14] = y[14]/ 2.80134000e-02; /*N2 (wt in kg) */

    return;
}


/*reverse of ytcr, useful for rate computations */
void fectyr_(double * c, double * rho, int * iwrk, double * rwrk, double * y)
{
    y[0] = c[0] * 1.007970 / (*rho); 
    y[1] = c[1] * 15.999400 / (*rho); 
    y[2] = c[2] * 17.007370 / (*rho); 
    y[3] = c[3] * 2.015940 / (*rho); 
    y[4] = c[4] * 31.998800 / (*rho); 
    y[5] = c[5] * 33.006770 / (*rho); 
    y[6] = c[6] * 18.015340 / (*rho); 
    y[7] = c[7] * 34.014740 / (*rho); 
    y[8] = c[8] * 30.006100 / (*rho); 
    y[9] = c[9] * 46.005500 / (*rho); 
    y[10] = c[10] * 44.012800 / (*rho); 
    y[11] = c[11] * 15.014670 / (*rho); 
    y[12] = c[12] * 14.006700 / (*rho); 
    y[13] = c[13] * 29.021370 / (*rho); 
    y[14] = c[14] * 28.013400 / (*rho); 

    return;
}


/*ddebdf compatible right hand side of CV burner */
/*rwrk[0] and rwrk[1] should contain rho and ene respectively */
/*working variable phi contains specific mole numbers */
void fecvrhs_(double * time, double * phi, double * phidot, double * rwrk, int * iwrk)
{
    double rho,ene; /*CV Parameters */
    double y[15], wdot[15]; /*temporary storage */
    int i; /*Loop counter */
    double temperature,pressure; /*temporary var */
    rho = rwrk[0];
    ene = rwrk[1];
    fephity_(phi, iwrk, rwrk, y);
    feeytt_(&ene, y, iwrk, rwrk, &temperature);
    CKPY(&rho, &temperature,  y, iwrk, rwrk, &pressure);
    CKWYP(&pressure, &temperature,  y, iwrk, rwrk, wdot);
    for (i=0; i<15; ++i) phidot[i] = wdot[i] / (rho/1000.0); 

    return;
}


/*returns the dimensionality of the cv burner (number of species) */
int fecvdim_()
{
    return 15;
}


/*ddebdf compatible right hand side of ZND solver */
/*rwrk[0] : scaling factor for pressure */
/*rwrk[1] : preshock density (g/cc)  */
/*rwrk[2] : detonation velocity (cm/s)  */
/*solution vector: [P; rho; y0 ... ylast]  */
void fezndrhs_(double * time, double * z, double * zdot, double * rwrk, int * iwrk)
{
    double psc,rho1,udet; /*ZND Parameters */
    double wt[15], hms[15], wdot[15]; /*temporary storage */
    int i; /*Loop counter */
    /*temporary variables */
    double ru, T, uvel, wtm, p, rho, gam, son, xm, sum, drdy, eta, cp, cv ;
    double *y; /*mass frac pointer */

    ru = 8.314e+07;

    psc = rwrk[0];
    rho1 = rwrk[1];
    udet = rwrk[2];

    p = z[0] * psc;
    rho = z[1];

    y = &z[3];

    CKMMWY(y, 0, 0, &wtm);

    T = p * wtm / rho / ru;

    uvel = (rho1 * udet)/ rho;

    CKCPBS(&T, y, 0, 0, &cp);
    CKCVBS(&T, y, 0, 0, &cv);
    gam = cp/cv;

    son = sqrt(fabs(gam*ru*T/wtm));
    xm = uvel/son;

    CKHMS(&T, 0, 0, hms);
    CKWT(0, 0, wt);
    CKWYP(&p, &T, y, 0, 0, wdot);

    sum = 0.0;
    for (i=0; i<15; ++i) {
        zdot[i+3] = wdot[i] * wt[i] / rho;
        drdy = -rho * wtm / wt[i];
        sum += -( drdy + rho * hms[i]/ (cp*T) ) * zdot[i+3];
    }

    eta = 1.0 - xm*xm;
    zdot[0] = -(uvel*uvel/eta/psc)*sum;
    zdot[1] = -sum/eta;
    zdot[2] = uvel;

    return;
}


/*returns the dimensionality of the ZND solver (3+number of species) */
int feznddim_()
{
    return 18;
}


/*returns the name of the source mechanism file  */
char* femechfile_()
{
    return "";
}


/*returns the species number */
int fesymnum_(const char* s1)
{
    if (strcmp(s1, "H")==0) return 0; 
    if (strcmp(s1, "O")==0) return 1; 
    if (strcmp(s1, "OH")==0) return 2; 
    if (strcmp(s1, "H2")==0) return 3; 
    if (strcmp(s1, "O2")==0) return 4; 
    if (strcmp(s1, "HO2")==0) return 5; 
    if (strcmp(s1, "H2O")==0) return 6; 
    if (strcmp(s1, "H2O2")==0) return 7; 
    if (strcmp(s1, "NO")==0) return 8; 
    if (strcmp(s1, "NO2")==0) return 9; 
    if (strcmp(s1, "N2O")==0) return 10; 
    if (strcmp(s1, "NH")==0) return 11; 
    if (strcmp(s1, "N")==0) return 12; 
    if (strcmp(s1, "NNH")==0) return 13; 
    if (strcmp(s1, "N2")==0) return 14; 
    /*species name not found */
    return -1;
}


/*returns the species name */
char* fesymname_(int sn)
{
    if (sn==0) return "H"; 
    if (sn==1) return "O"; 
    if (sn==2) return "OH"; 
    if (sn==3) return "H2"; 
    if (sn==4) return "O2"; 
    if (sn==5) return "HO2"; 
    if (sn==6) return "H2O"; 
    if (sn==7) return "H2O2"; 
    if (sn==8) return "NO"; 
    if (sn==9) return "NO2"; 
    if (sn==10) return "N2O"; 
    if (sn==11) return "NH"; 
    if (sn==12) return "N"; 
    if (sn==13) return "NNH"; 
    if (sn==14) return "N2"; 
    /*species name not found */
    return "NOTFOUND";
}

/* End of file  */
