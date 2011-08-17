/*  -*- C -*-  */
/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *                                -*- author -*-
 *                             -*- organization -*-
 *                    (C) -*- years -*-  All Rights Reserved
 *
 * <LicenseText>
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

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
    *mm = 4;
    *kk = 10;
    *ii = 31;
    *nfit = -1; /*Why do you need this anyway ?  */
}


/*Dummy ckinit */
void fginit_(int * leniwk, int * lenrwk, int * lencwk, int * linc, int * lout, int * ickwrk, double * rckwrk, char * cckwrk )
{
    if ((*lout) != 0) {
        printf(" ***       Congratulations       *** \n");
        printf(" * You are using the Fuego Library * \n");
        printf(" *****    Say NO to cklib.f    ***** \n");
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
void CKSYME(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*4; i++) {
        kname[i] = ' ';
    }

    /* O  */
    kname[ 0*lenkname + 0 ] = 'O';
    kname[ 0*lenkname + 1 ] = ' ';

    /* N  */
    kname[ 1*lenkname + 0 ] = 'N';
    kname[ 1*lenkname + 1 ] = ' ';

    /* AR  */
    kname[ 2*lenkname + 0 ] = 'A';
    kname[ 2*lenkname + 1 ] = 'R';
    kname[ 2*lenkname + 2 ] = ' ';

    /* H  */
    kname[ 3*lenkname + 0 ] = 'H';
    kname[ 3*lenkname + 1 ] = ' ';

}


/* Returns the char strings of species names */
void CKSYMS(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*10; i++) {
        kname[i] = ' ';
    }

    /* H  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = ' ';

    /* H2  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = '2';
    kname[ 1*lenkname + 2 ] = ' ';

    /* O  */
    kname[ 2*lenkname + 0 ] = 'O';
    kname[ 2*lenkname + 1 ] = ' ';

    /* O2  */
    kname[ 3*lenkname + 0 ] = 'O';
    kname[ 3*lenkname + 1 ] = '2';
    kname[ 3*lenkname + 2 ] = ' ';

    /* H2O  */
    kname[ 4*lenkname + 0 ] = 'H';
    kname[ 4*lenkname + 1 ] = '2';
    kname[ 4*lenkname + 2 ] = 'O';
    kname[ 4*lenkname + 3 ] = ' ';

    /* OH  */
    kname[ 5*lenkname + 0 ] = 'O';
    kname[ 5*lenkname + 1 ] = 'H';
    kname[ 5*lenkname + 2 ] = ' ';

    /* H2O2  */
    kname[ 6*lenkname + 0 ] = 'H';
    kname[ 6*lenkname + 1 ] = '2';
    kname[ 6*lenkname + 2 ] = 'O';
    kname[ 6*lenkname + 3 ] = '2';
    kname[ 6*lenkname + 4 ] = ' ';

    /* HO2  */
    kname[ 7*lenkname + 0 ] = 'H';
    kname[ 7*lenkname + 1 ] = 'O';
    kname[ 7*lenkname + 2 ] = '2';
    kname[ 7*lenkname + 3 ] = ' ';

    /* AR  */
    kname[ 8*lenkname + 0 ] = 'A';
    kname[ 8*lenkname + 1 ] = 'R';
    kname[ 8*lenkname + 2 ] = ' ';

    /* N2  */
    kname[ 9*lenkname + 0 ] = 'N';
    kname[ 9*lenkname + 1 ] = '2';
    kname[ 9*lenkname + 2 ] = ' ';

}


/* Returns R, Rc, Patm */
void CKRP(int * ickwrk, double * rckwrk, double * ru, double * ruc, double * pa)
{
     *ru  = 8.31451e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}


/*Compute P = rhoRT/W(x) */
void CKPX(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * P)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*34.014740; /*H2O2 */
    XW += x[7]*33.006770; /*HO2 */
    XW += x[8]*39.948000; /*AR */
    XW += x[9]*28.013400; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void CKPY(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]/1.007970; /*H */
    YOW += y[1]/2.015940; /*H2 */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/18.015340; /*H2O */
    YOW += y[5]/17.007370; /*OH */
    YOW += y[6]/34.014740; /*H2O2 */
    YOW += y[7]/33.006770; /*HO2 */
    YOW += y[8]/39.948000; /*AR */
    YOW += y[9]/28.013400; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) * YOW; /*P = rho*R*T/W */

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
    W += c[1]*2.015940; /*H2 */
    W += c[2]*15.999400; /*O */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*17.007370; /*OH */
    W += c[6]*34.014740; /*H2O2 */
    W += c[7]*33.006770; /*HO2 */
    W += c[8]*39.948000; /*AR */
    W += c[9]*28.013400; /*N2 */

    for (id = 0; id < 10; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31451e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double * P, double * T, double * x, int * iwrk, double * rwrk, double * rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*34.014740; /*H2O2 */
    XW += x[7]*33.006770; /*HO2 */
    XW += x[8]*39.948000; /*AR */
    XW += x[9]*28.013400; /*N2 */
    *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double * P, double * T, double * y, int * iwrk, double * rwrk, double * rho)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]/1.007970; /*H */
    YOW += y[1]/2.015940; /*H2 */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/18.015340; /*H2O */
    YOW += y[5]/17.007370; /*OH */
    YOW += y[6]/34.014740; /*H2O2 */
    YOW += y[7]/33.006770; /*HO2 */
    YOW += y[8]/39.948000; /*AR */
    YOW += y[9]/28.013400; /*N2 */
    *rho = *P / (8.31451e+07 * (*T) * YOW); /*rho = P*W/(R*T) */

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
    W += c[1]*2.015940; /*H2 */
    W += c[2]*15.999400; /*O */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*17.007370; /*OH */
    W += c[6]*34.014740; /*H2O2 */
    W += c[7]*33.006770; /*HO2 */
    W += c[8]*39.948000; /*AR */
    W += c[9]*28.013400; /*N2 */

    for (id = 0; id < 10; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.31451e+07); /*rho = PW/(R*T) */

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
    YOW += y[1]/2.015940; /*H2 */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/18.015340; /*H2O */
    YOW += y[5]/17.007370; /*OH */
    YOW += y[6]/34.014740; /*H2O2 */
    YOW += y[7]/33.006770; /*HO2 */
    YOW += y[8]/39.948000; /*AR */
    YOW += y[9]/28.013400; /*N2 */
    *wtm = 1.0 / YOW;

    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWX(double *x, int * iwrk, double * rwrk, double * wtm)
{
    double XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*34.014740; /*H2O2 */
    XW += x[7]*33.006770; /*HO2 */
    XW += x[8]*39.948000; /*AR */
    XW += x[9]*28.013400; /*N2 */
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
    W += c[1]*2.015940; /*H2 */
    W += c[2]*15.999400; /*O */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*17.007370; /*OH */
    W += c[6]*34.014740; /*H2O2 */
    W += c[7]*33.006770; /*HO2 */
    W += c[8]*39.948000; /*AR */
    W += c[9]*28.013400; /*N2 */

    for (id = 0; id < 10; ++id) {
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
    YOW += y[1]/2.015940; /*H2 */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/18.015340; /*H2O */
    YOW += y[5]/17.007370; /*OH */
    YOW += y[6]/34.014740; /*H2O2 */
    YOW += y[7]/33.006770; /*HO2 */
    YOW += y[8]/39.948000; /*AR */
    YOW += y[9]/28.013400; /*N2 */
    /*Now compute conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(2.015940*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(17.007370*YOW); 
    x[6] = y[6]/(34.014740*YOW); 
    x[7] = y[7]/(33.006770*YOW); 
    x[8] = y[8]/(39.948000*YOW); 
    x[9] = y[9]/(28.013400*YOW); 

    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * c)
{
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/1.007970; /*H */
    YOW += y[1]/2.015940; /*H2 */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/18.015340; /*H2O */
    YOW += y[5]/17.007370; /*OH */
    YOW += y[6]/34.014740; /*H2O2 */
    YOW += y[7]/33.006770; /*HO2 */
    YOW += y[8]/39.948000; /*AR */
    YOW += y[9]/28.013400; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*Now compute conversion */
    c[0] = PWORT * y[0]/1.007970; 
    c[1] = PWORT * y[1]/2.015940; 
    c[2] = PWORT * y[2]/15.999400; 
    c[3] = PWORT * y[3]/31.998800; 
    c[4] = PWORT * y[4]/18.015340; 
    c[5] = PWORT * y[5]/17.007370; 
    c[6] = PWORT * y[6]/34.014740; 
    c[7] = PWORT * y[7]/33.006770; 
    c[8] = PWORT * y[8]/39.948000; 
    c[9] = PWORT * y[9]/28.013400; 

    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * c)
{
    /*See Eq 8 (Temperature not used) */
    c[0] = (*rho) * y[0]/1.007970; 
    c[1] = (*rho) * y[1]/2.015940; 
    c[2] = (*rho) * y[2]/15.999400; 
    c[3] = (*rho) * y[3]/31.998800; 
    c[4] = (*rho) * y[4]/18.015340; 
    c[5] = (*rho) * y[5]/17.007370; 
    c[6] = (*rho) * y[6]/34.014740; 
    c[7] = (*rho) * y[7]/33.006770; 
    c[8] = (*rho) * y[8]/39.948000; 
    c[9] = (*rho) * y[9]/28.013400; 

    return;
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
void CKXTY(double * x, int * iwrk, double * rwrk, double * y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*34.014740; /*H2O2 */
    XW += x[7]*33.006770; /*HO2 */
    XW += x[8]*39.948000; /*AR */
    XW += x[9]*28.013400; /*N2 */
    /*Now compute conversion */
    y[0] = x[0]*1.007970/XW; 
    y[1] = x[1]*2.015940/XW; 
    y[2] = x[2]*15.999400/XW; 
    y[3] = x[3]*31.998800/XW; 
    y[4] = x[4]*18.015340/XW; 
    y[5] = x[5]*17.007370/XW; 
    y[6] = x[6]*34.014740/XW; 
    y[7] = x[7]*33.006770/XW; 
    y[8] = x[8]*39.948000/XW; 
    y[9] = x[9]*28.013400/XW; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 10; ++id) {
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
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*34.014740; /*H2O2 */
    XW += x[7]*33.006770; /*HO2 */
    XW += x[8]*39.948000; /*AR */
    XW += x[9]*28.013400; /*N2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 10; ++id) {
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
    for (id = 0; id < 10; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    for (id = 0; id < 10; ++id) {
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
    CW += c[1]*2.015940; /*H2 */
    CW += c[2]*15.999400; /*O */
    CW += c[3]*31.998800; /*O2 */
    CW += c[4]*18.015340; /*H2O */
    CW += c[5]*17.007370; /*OH */
    CW += c[6]*34.014740; /*H2O2 */
    CW += c[7]*33.006770; /*HO2 */
    CW += c[8]*39.948000; /*AR */
    CW += c[9]*28.013400; /*N2 */
    /*Now compute conversion */
    y[0] = c[0]*1.007970/CW; 
    y[1] = c[1]*2.015940/CW; 
    y[2] = c[2]*15.999400/CW; 
    y[3] = c[3]*31.998800/CW; 
    y[4] = c[4]*18.015340/CW; 
    y[5] = c[5]*17.007370/CW; 
    y[6] = c[6]*34.014740/CW; 
    y[7] = c[7]*33.006770/CW; 
    y[8] = c[8]*39.948000/CW; 
    y[9] = c[9]*28.013400/CW; 

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
    for (id = 0; id < 10; ++id) {
        cvml[id] *= 8.31451e+07;
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
    for (id = 0; id < 10; ++id) {
        cpml[id] *= 8.31451e+07;
    }
}


/*get internal energy as a function  */
/*of T for all species (molar units) */
void CKUML(double *T, int * iwrk, double * rwrk, double * uml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 10; ++id) {
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
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 10; ++id) {
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
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 10; ++id) {
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
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(aml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 10; ++id) {
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
    for (id = 0; id < 10; ++id) {
        sml[id] *= 8.31451e+07;
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
    cvms[0] *= 8.248767324424338e+07; /*H */
    cvms[1] *= 4.124383662212169e+07; /*H2 */
    cvms[2] *= 5.196763628636074e+06; /*O */
    cvms[3] *= 2.598381814318037e+06; /*O2 */
    cvms[4] *= 4.615239012974499e+06; /*H2O */
    cvms[5] *= 4.888768810227566e+06; /*OH */
    cvms[6] *= 2.444384405113783e+06; /*H2O2 */
    cvms[7] *= 2.519031701678171e+06; /*HO2 */
    cvms[8] *= 2.081333233203164e+06; /*AR */
    cvms[9] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
void CKCPMS(double *T, int * iwrk, double * rwrk, double * cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 8.248767324424338e+07; /*H */
    cpms[1] *= 4.124383662212169e+07; /*H2 */
    cpms[2] *= 5.196763628636074e+06; /*O */
    cpms[3] *= 2.598381814318037e+06; /*O2 */
    cpms[4] *= 4.615239012974499e+06; /*H2O */
    cpms[5] *= 4.888768810227566e+06; /*OH */
    cpms[6] *= 2.444384405113783e+06; /*H2O2 */
    cpms[7] *= 2.519031701678171e+06; /*HO2 */
    cpms[8] *= 2.081333233203164e+06; /*AR */
    cpms[9] *= 2.968047434442088e+06; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double *T, int * iwrk, double * rwrk, double * ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    ums[0] *= RT/1.007970; /*H */
    ums[1] *= RT/2.015940; /*H2 */
    ums[2] *= RT/15.999400; /*O */
    ums[3] *= RT/31.998800; /*O2 */
    ums[4] *= RT/18.015340; /*H2O */
    ums[5] *= RT/17.007370; /*OH */
    ums[6] *= RT/34.014740; /*H2O2 */
    ums[7] *= RT/33.006770; /*HO2 */
    ums[8] *= RT/39.948000; /*AR */
    ums[9] *= RT/28.013400; /*N2 */
}


/*Returns enthalpy in mass units (Eq 27.) */
void CKHMS(double *T, int * iwrk, double * rwrk, double * hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hms, tc);
    hms[0] *= RT/1.007970; /*H */
    hms[1] *= RT/2.015940; /*H2 */
    hms[2] *= RT/15.999400; /*O */
    hms[3] *= RT/31.998800; /*O2 */
    hms[4] *= RT/18.015340; /*H2O */
    hms[5] *= RT/17.007370; /*OH */
    hms[6] *= RT/34.014740; /*H2O2 */
    hms[7] *= RT/33.006770; /*HO2 */
    hms[8] *= RT/39.948000; /*AR */
    hms[9] *= RT/28.013400; /*N2 */
}


/*Returns gibbs in mass units (Eq 31.) */
void CKGMS(double *T, int * iwrk, double * rwrk, double * gms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gms, tc);
    gms[0] *= RT/1.007970; /*H */
    gms[1] *= RT/2.015940; /*H2 */
    gms[2] *= RT/15.999400; /*O */
    gms[3] *= RT/31.998800; /*O2 */
    gms[4] *= RT/18.015340; /*H2O */
    gms[5] *= RT/17.007370; /*OH */
    gms[6] *= RT/34.014740; /*H2O2 */
    gms[7] *= RT/33.006770; /*HO2 */
    gms[8] *= RT/39.948000; /*AR */
    gms[9] *= RT/28.013400; /*N2 */
}


/*Returns helmholtz in mass units (Eq 32.) */
void CKAMS(double *T, int * iwrk, double * rwrk, double * ams)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(ams, tc);
    ams[0] *= RT/1.007970; /*H */
    ams[1] *= RT/2.015940; /*H2 */
    ams[2] *= RT/15.999400; /*O */
    ams[3] *= RT/31.998800; /*O2 */
    ams[4] *= RT/18.015340; /*H2O */
    ams[5] *= RT/17.007370; /*OH */
    ams[6] *= RT/34.014740; /*H2O2 */
    ams[7] *= RT/33.006770; /*HO2 */
    ams[8] *= RT/39.948000; /*AR */
    ams[9] *= RT/28.013400; /*N2 */
}


/*Returns the entropies in mass units (Eq 28.) */
void CKSMS(double *T, int * iwrk, double * rwrk, double * sms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 8.248767324424338e+07; /*H */
    sms[1] *= 4.124383662212169e+07; /*H2 */
    sms[2] *= 5.196763628636074e+06; /*O */
    sms[3] *= 2.598381814318037e+06; /*O2 */
    sms[4] *= 4.615239012974499e+06; /*H2O */
    sms[5] *= 4.888768810227566e+06; /*OH */
    sms[6] *= 2.444384405113783e+06; /*H2O2 */
    sms[7] *= 2.519031701678171e+06; /*HO2 */
    sms[8] *= 2.081333233203164e+06; /*AR */
    sms[9] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double *T, double *x, int * iwrk, double * rwrk, double * cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[10]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 10; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
void CKCPBS(double *T, double *y, int * iwrk, double * rwrk, double * cpbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[10]; /* temporary storage */
    cp_R(cpor, tc);
    /*multiply by y/molecularweight */
    result += cpor[0]*y[0]/1.007970; /*H */
    result += cpor[1]*y[1]/2.015940; /*H2 */
    result += cpor[2]*y[2]/15.999400; /*O */
    result += cpor[3]*y[3]/31.998800; /*O2 */
    result += cpor[4]*y[4]/18.015340; /*H2O */
    result += cpor[5]*y[5]/17.007370; /*OH */
    result += cpor[6]*y[6]/34.014740; /*H2O2 */
    result += cpor[7]*y[7]/33.006770; /*HO2 */
    result += cpor[8]*y[8]/39.948000; /*AR */
    result += cpor[9]*y[9]/28.013400; /*N2 */

    *cpbs = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double *T, double *x, int * iwrk, double * rwrk, double * cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[10]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 10; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
void CKCVBS(double *T, double *y, int * iwrk, double * rwrk, double * cvbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[10]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]/1.007970; /*H */
    result += cvor[1]*y[1]/2.015940; /*H2 */
    result += cvor[2]*y[2]/15.999400; /*O */
    result += cvor[3]*y[3]/31.998800; /*O2 */
    result += cvor[4]*y[4]/18.015340; /*H2O */
    result += cvor[5]*y[5]/17.007370; /*OH */
    result += cvor[6]*y[6]/34.014740; /*H2O2 */
    result += cvor[7]*y[7]/33.006770; /*HO2 */
    result += cvor[8]*y[8]/39.948000; /*AR */
    result += cvor[9]*y[9]/28.013400; /*N2 */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double *T, double *x, int * iwrk, double * rwrk, double * hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[10]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 10; ++id) {
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
    double hml[10]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*hml[0]/1.007970; /*H */
    result += y[1]*hml[1]/2.015940; /*H2 */
    result += y[2]*hml[2]/15.999400; /*O */
    result += y[3]*hml[3]/31.998800; /*O2 */
    result += y[4]*hml[4]/18.015340; /*H2O */
    result += y[5]*hml[5]/17.007370; /*OH */
    result += y[6]*hml[6]/34.014740; /*H2O2 */
    result += y[7]*hml[7]/33.006770; /*HO2 */
    result += y[8]*hml[8]/39.948000; /*AR */
    result += y[9]*hml[9]/28.013400; /*N2 */

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double *T, double *x, int * iwrk, double * rwrk, double * ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[10]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 10; ++id) {
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
    double ums[10]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]/1.007970; /*H */
    result += y[1]*ums[1]/2.015940; /*H2 */
    result += y[2]*ums[2]/15.999400; /*O */
    result += y[3]*ums[3]/31.998800; /*O2 */
    result += y[4]*ums[4]/18.015340; /*H2O */
    result += y[5]*ums[5]/17.007370; /*OH */
    result += y[6]*ums[6]/34.014740; /*H2O2 */
    result += y[7]*ums[7]/33.006770; /*HO2 */
    result += y[8]*ums[8]/39.948000; /*AR */
    result += y[9]*ums[9]/28.013400; /*N2 */

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
    double sor[10]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 10; ++id) {
        result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
    }

    *sbml = result * 8.31451e+07;
}


/*get mixture entropy in mass units */
void CKSBMS(double *P, double *T, double *y, int * iwrk, double * rwrk, double * sbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[10]; /* temporary storage */
    double x[10]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/1.007970; /*H */
    YOW += y[1]/2.015940; /*H2 */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/18.015340; /*H2O */
    YOW += y[5]/17.007370; /*OH */
    YOW += y[6]/34.014740; /*H2O2 */
    YOW += y[7]/33.006770; /*HO2 */
    YOW += y[8]/39.948000; /*AR */
    YOW += y[9]/28.013400; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(2.015940*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(17.007370*YOW); 
    x[6] = y[6]/(34.014740*YOW); 
    x[7] = y[7]/(33.006770*YOW); 
    x[8] = y[8]/(39.948000*YOW); 
    x[9] = y[9]/(28.013400*YOW); 
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
    /*Scale by R/W */
    *sbms = result * 8.31451e+07 * YOW;
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
    double RT = 8.31451e+07*tT; /*R*T */
    double gort[10]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 10; ++id) {
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
    double RT = 8.31451e+07*tT; /*R*T */
    double gort[10]; /* temporary storage */
    double x[10]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/1.007970; /*H */
    YOW += y[1]/2.015940; /*H2 */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/18.015340; /*H2O */
    YOW += y[5]/17.007370; /*OH */
    YOW += y[6]/34.014740; /*H2O2 */
    YOW += y[7]/33.006770; /*HO2 */
    YOW += y[8]/39.948000; /*AR */
    YOW += y[9]/28.013400; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(2.015940*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(17.007370*YOW); 
    x[6] = y[6]/(34.014740*YOW); 
    x[7] = y[7]/(33.006770*YOW); 
    x[8] = y[8]/(39.948000*YOW); 
    x[9] = y[9]/(28.013400*YOW); 
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
    double RT = 8.31451e+07*tT; /*R*T */
    double aort[10]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 10; ++id) {
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
    double RT = 8.31451e+07*tT; /*R*T */
    double aort[10]; /* temporary storage */
    double x[10]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/1.007970; /*H */
    YOW += y[1]/2.015940; /*H2 */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/18.015340; /*H2O */
    YOW += y[5]/17.007370; /*OH */
    YOW += y[6]/34.014740; /*H2O2 */
    YOW += y[7]/33.006770; /*HO2 */
    YOW += y[8]/39.948000; /*AR */
    YOW += y[9]/28.013400; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(2.015940*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(17.007370*YOW); 
    x[6] = y[6]/(34.014740*YOW); 
    x[7] = y[7]/(33.006770*YOW); 
    x[8] = y[8]/(39.948000*YOW); 
    x[9] = y[9]/(28.013400*YOW); 
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
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double * T, double * C, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 10; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 10; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[10]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/1.007970; /*H */
    YOW += y[1]/2.015940; /*H2 */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/18.015340; /*H2O */
    YOW += y[5]/17.007370; /*OH */
    YOW += y[6]/34.014740; /*H2O2 */
    YOW += y[7]/33.006770; /*HO2 */
    YOW += y[8]/39.948000; /*AR */
    YOW += y[9]/28.013400; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]/1.007970; 
    c[1] = PWORT * y[1]/2.015940; 
    c[2] = PWORT * y[2]/15.999400; 
    c[3] = PWORT * y[3]/31.998800; 
    c[4] = PWORT * y[4]/18.015340; 
    c[5] = PWORT * y[5]/17.007370; 
    c[6] = PWORT * y[6]/34.014740; 
    c[7] = PWORT * y[7]/33.006770; 
    c[8] = PWORT * y[8]/39.948000; 
    c[9] = PWORT * y[9]/28.013400; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 10; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[10]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 10; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 10; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[10]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]/1.007970; 
    c[1] = 1e6 * (*rho) * y[1]/2.015940; 
    c[2] = 1e6 * (*rho) * y[2]/15.999400; 
    c[3] = 1e6 * (*rho) * y[3]/31.998800; 
    c[4] = 1e6 * (*rho) * y[4]/18.015340; 
    c[5] = 1e6 * (*rho) * y[5]/17.007370; 
    c[6] = 1e6 * (*rho) * y[6]/34.014740; 
    c[7] = 1e6 * (*rho) * y[7]/33.006770; 
    c[8] = 1e6 * (*rho) * y[8]/39.948000; 
    c[9] = 1e6 * (*rho) * y[9]/28.013400; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 10; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[10]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*34.014740; /*H2O2 */
    XW += x[7]*33.006770; /*HO2 */
    XW += x[8]*39.948000; /*AR */
    XW += x[9]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 10; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 10; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double * T, double * C, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 10; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 10; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 31; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double * P, double * T, double * x, int * iwrk, double * rwrk, double * q_f, double * q_r)
{
    int id; /*loop counter */
    double c[10]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 10; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 31; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[10]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/1.007970; /*H */
    YOW += y[1]/2.015940; /*H2 */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/18.015340; /*H2O */
    YOW += y[5]/17.007370; /*OH */
    YOW += y[6]/34.014740; /*H2O2 */
    YOW += y[7]/33.006770; /*HO2 */
    YOW += y[8]/39.948000; /*AR */
    YOW += y[9]/28.013400; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]/1.007970; 
    c[1] = PWORT * y[1]/2.015940; 
    c[2] = PWORT * y[2]/15.999400; 
    c[3] = PWORT * y[3]/31.998800; 
    c[4] = PWORT * y[4]/18.015340; 
    c[5] = PWORT * y[5]/17.007370; 
    c[6] = PWORT * y[6]/34.014740; 
    c[7] = PWORT * y[7]/33.006770; 
    c[8] = PWORT * y[8]/39.948000; 
    c[9] = PWORT * y[9]/28.013400; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 31; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[10]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 10; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 31; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[10]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]/1.007970; 
    c[1] = 1e6 * (*rho) * y[1]/2.015940; 
    c[2] = 1e6 * (*rho) * y[2]/15.999400; 
    c[3] = 1e6 * (*rho) * y[3]/31.998800; 
    c[4] = 1e6 * (*rho) * y[4]/18.015340; 
    c[5] = 1e6 * (*rho) * y[5]/17.007370; 
    c[6] = 1e6 * (*rho) * y[6]/34.014740; 
    c[7] = 1e6 * (*rho) * y[7]/33.006770; 
    c[8] = 1e6 * (*rho) * y[8]/39.948000; 
    c[9] = 1e6 * (*rho) * y[9]/28.013400; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 31; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[10]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*34.014740; /*H2O2 */
    XW += x[7]*33.006770; /*HO2 */
    XW += x[8]*39.948000; /*AR */
    XW += x[9]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 10; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 31; ++id) {
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
    for (id = 0; id < 10 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: H + H + M <=> H2 + M */
    nuki[ 0 * kd + 0 ] += -1 ;
    nuki[ 0 * kd + 0 ] += -1 ;
    nuki[ 1 * kd + 0 ] += +1 ;

    /*reaction 2: H + H + H2 <=> H2 + H2 */
    nuki[ 0 * kd + 1 ] += -1 ;
    nuki[ 0 * kd + 1 ] += -1 ;
    nuki[ 1 * kd + 1 ] += -1 ;
    nuki[ 1 * kd + 1 ] += +1 ;
    nuki[ 1 * kd + 1 ] += +1 ;

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    nuki[ 0 * kd + 2 ] += -1 ;
    nuki[ 0 * kd + 2 ] += -1 ;
    nuki[ 9 * kd + 2 ] += -1 ;
    nuki[ 1 * kd + 2 ] += +1 ;
    nuki[ 9 * kd + 2 ] += +1 ;

    /*reaction 4: H + H + H <=> H2 + H */
    nuki[ 0 * kd + 3 ] += -1 ;
    nuki[ 0 * kd + 3 ] += -1 ;
    nuki[ 0 * kd + 3 ] += -1 ;
    nuki[ 1 * kd + 3 ] += +1 ;
    nuki[ 0 * kd + 3 ] += +1 ;

    /*reaction 5: O + O + M <=> O2 + M */
    nuki[ 2 * kd + 4 ] += -1 ;
    nuki[ 2 * kd + 4 ] += -1 ;
    nuki[ 3 * kd + 4 ] += +1 ;

    /*reaction 6: O + H + M <=> OH + M */
    nuki[ 2 * kd + 5 ] += -1 ;
    nuki[ 0 * kd + 5 ] += -1 ;
    nuki[ 5 * kd + 5 ] += +1 ;

    /*reaction 7: H + OH + M <=> H2O + M */
    nuki[ 0 * kd + 6 ] += -1 ;
    nuki[ 5 * kd + 6 ] += -1 ;
    nuki[ 4 * kd + 6 ] += +1 ;

    /*reaction 8: H + O2 (+M) <=> HO2 (+M) */
    nuki[ 0 * kd + 7 ] += -1 ;
    nuki[ 3 * kd + 7 ] += -1 ;
    nuki[ 7 * kd + 7 ] += +1 ;

    /*reaction 9: H + O2 (+AR) <=> HO2 (+AR) */
    nuki[ 0 * kd + 8 ] += -1 ;
    nuki[ 3 * kd + 8 ] += -1 ;
    nuki[ 7 * kd + 8 ] += +1 ;

    /*reaction 10: H + O2 (+O2) <=> HO2 (+O2) */
    nuki[ 0 * kd + 9 ] += -1 ;
    nuki[ 3 * kd + 9 ] += -1 ;
    nuki[ 7 * kd + 9 ] += +1 ;

    /*reaction 11: H + O2 (+H2O) <=> HO2 (+H2O) */
    nuki[ 0 * kd + 10 ] += -1 ;
    nuki[ 3 * kd + 10 ] += -1 ;
    nuki[ 7 * kd + 10 ] += +1 ;

    /*reaction 12: OH + OH (+M) <=> H2O2 (+M) */
    nuki[ 5 * kd + 11 ] += -1 ;
    nuki[ 5 * kd + 11 ] += -1 ;
    nuki[ 6 * kd + 11 ] += +1 ;

    /*reaction 13: OH + OH (+H2O) <=> H2O2 (+H2O) */
    nuki[ 5 * kd + 12 ] += -1 ;
    nuki[ 5 * kd + 12 ] += -1 ;
    nuki[ 6 * kd + 12 ] += +1 ;

    /*reaction 14: O + H2 <=> OH + H */
    nuki[ 2 * kd + 13 ] += -1 ;
    nuki[ 1 * kd + 13 ] += -1 ;
    nuki[ 5 * kd + 13 ] += +1 ;
    nuki[ 0 * kd + 13 ] += +1 ;

    /*reaction 15: H + O2 <=> OH + O */
    nuki[ 0 * kd + 14 ] += -1 ;
    nuki[ 3 * kd + 14 ] += -1 ;
    nuki[ 5 * kd + 14 ] += +1 ;
    nuki[ 2 * kd + 14 ] += +1 ;

    /*reaction 16: H2 + OH <=> H2O + H */
    nuki[ 1 * kd + 15 ] += -1 ;
    nuki[ 5 * kd + 15 ] += -1 ;
    nuki[ 4 * kd + 15 ] += +1 ;
    nuki[ 0 * kd + 15 ] += +1 ;

    /*reaction 17: OH + OH <=> H2O + O */
    nuki[ 5 * kd + 16 ] += -1 ;
    nuki[ 5 * kd + 16 ] += -1 ;
    nuki[ 4 * kd + 16 ] += +1 ;
    nuki[ 2 * kd + 16 ] += +1 ;

    /*reaction 18: HO2 + O <=> OH + O2 */
    nuki[ 7 * kd + 17 ] += -1 ;
    nuki[ 2 * kd + 17 ] += -1 ;
    nuki[ 5 * kd + 17 ] += +1 ;
    nuki[ 3 * kd + 17 ] += +1 ;

    /*reaction 19: H + HO2 <=> OH + OH */
    nuki[ 0 * kd + 18 ] += -1 ;
    nuki[ 7 * kd + 18 ] += -1 ;
    nuki[ 5 * kd + 18 ] += +1 ;
    nuki[ 5 * kd + 18 ] += +1 ;

    /*reaction 20: H + HO2 <=> H2O + O */
    nuki[ 0 * kd + 19 ] += -1 ;
    nuki[ 7 * kd + 19 ] += -1 ;
    nuki[ 4 * kd + 19 ] += +1 ;
    nuki[ 2 * kd + 19 ] += +1 ;

    /*reaction 21: H2 + O2 <=> H + HO2 */
    nuki[ 1 * kd + 20 ] += -1 ;
    nuki[ 3 * kd + 20 ] += -1 ;
    nuki[ 0 * kd + 20 ] += +1 ;
    nuki[ 7 * kd + 20 ] += +1 ;

    /*reaction 22: H2 + O2 <=> OH + OH */
    nuki[ 1 * kd + 21 ] += -1 ;
    nuki[ 3 * kd + 21 ] += -1 ;
    nuki[ 5 * kd + 21 ] += +1 ;
    nuki[ 5 * kd + 21 ] += +1 ;

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    nuki[ 7 * kd + 22 ] += -1 ;
    nuki[ 5 * kd + 22 ] += -1 ;
    nuki[ 4 * kd + 22 ] += +1 ;
    nuki[ 3 * kd + 22 ] += +1 ;

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    nuki[ 7 * kd + 23 ] += -1 ;
    nuki[ 5 * kd + 23 ] += -1 ;
    nuki[ 4 * kd + 23 ] += +1 ;
    nuki[ 3 * kd + 23 ] += +1 ;

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    nuki[ 7 * kd + 24 ] += -1 ;
    nuki[ 7 * kd + 24 ] += -1 ;
    nuki[ 6 * kd + 24 ] += +1 ;
    nuki[ 3 * kd + 24 ] += +1 ;

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    nuki[ 7 * kd + 25 ] += -1 ;
    nuki[ 7 * kd + 25 ] += -1 ;
    nuki[ 6 * kd + 25 ] += +1 ;
    nuki[ 3 * kd + 25 ] += +1 ;

    /*reaction 27: H2O2 + H <=> HO2 + H2 */
    nuki[ 6 * kd + 26 ] += -1 ;
    nuki[ 0 * kd + 26 ] += -1 ;
    nuki[ 7 * kd + 26 ] += +1 ;
    nuki[ 1 * kd + 26 ] += +1 ;

    /*reaction 28: H2O2 + H <=> H2O + OH */
    nuki[ 6 * kd + 27 ] += -1 ;
    nuki[ 0 * kd + 27 ] += -1 ;
    nuki[ 4 * kd + 27 ] += +1 ;
    nuki[ 5 * kd + 27 ] += +1 ;

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    nuki[ 6 * kd + 28 ] += -1 ;
    nuki[ 2 * kd + 28 ] += -1 ;
    nuki[ 7 * kd + 28 ] += +1 ;
    nuki[ 5 * kd + 28 ] += +1 ;

    /*reaction 30: H2O2 + OH <=> HO2 + H2O */
    nuki[ 6 * kd + 29 ] += -1 ;
    nuki[ 5 * kd + 29 ] += -1 ;
    nuki[ 7 * kd + 29 ] += +1 ;
    nuki[ 4 * kd + 29 ] += +1 ;

    /*reaction 31: H2O2 + OH <=> HO2 + H2O */
    nuki[ 6 * kd + 30 ] += -1 ;
    nuki[ 5 * kd + 30 ] += -1 ;
    nuki[ 7 * kd + 30 ] += +1 ;
    nuki[ 4 * kd + 30 ] += +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < 4 * 10; ++ id) {
         ncf[id] = 0; 
    }

    /*H */
    ncf[ 0 * kd + 3 ] = 1; /*H */

    /*H2 */
    ncf[ 1 * kd + 3 ] = 2; /*H */

    /*O */
    ncf[ 2 * kd + 0 ] = 1; /*O */

    /*O2 */
    ncf[ 3 * kd + 0 ] = 2; /*O */

    /*H2O */
    ncf[ 4 * kd + 3 ] = 2; /*H */
    ncf[ 4 * kd + 0 ] = 1; /*O */

    /*OH */
    ncf[ 5 * kd + 0 ] = 1; /*O */
    ncf[ 5 * kd + 3 ] = 1; /*H */

    /*H2O2 */
    ncf[ 6 * kd + 3 ] = 2; /*H */
    ncf[ 6 * kd + 0 ] = 2; /*O */

    /*HO2 */
    ncf[ 7 * kd + 3 ] = 1; /*H */
    ncf[ 7 * kd + 0 ] = 2; /*O */

    /*AR */
    ncf[ 8 * kd + 2 ] = 1; /*AR */

    /*N2 */
    ncf[ 9 * kd + 1 ] = 2; /*N */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE(int * iwrk, double * rwrk, double * a, double * b, double * e)
{

    /*reaction 1: H + H + M <=> H2 + M */
    a[0] = 7e+17;
    b[0] = -1;
    e[0] = 0;

    /*reaction 2: H + H + H2 <=> H2 + H2 */
    a[1] = 1e+17;
    b[1] = -0.6;
    e[1] = 0;

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    a[2] = 5.4e+18;
    b[2] = -1.3;
    e[2] = 0;

    /*reaction 4: H + H + H <=> H2 + H */
    a[3] = 3.2e+15;
    b[3] = 0;
    e[3] = 0;

    /*reaction 5: O + O + M <=> O2 + M */
    a[4] = 1e+17;
    b[4] = -1;
    e[4] = 0;

    /*reaction 6: O + H + M <=> OH + M */
    a[5] = 6.75e+18;
    b[5] = -1;
    e[5] = 0;

    /*reaction 7: H + OH + M <=> H2O + M */
    a[6] = 2.2e+22;
    b[6] = -2;
    e[6] = 0;

    /*reaction 8: H + O2 (+M) <=> HO2 (+M) */
    a[7] = 4.66e+12;
    b[7] = 0.44;
    e[7] = 0;

    /*reaction 9: H + O2 (+AR) <=> HO2 (+AR) */
    a[8] = 4.66e+12;
    b[8] = 0.44;
    e[8] = 0;

    /*reaction 10: H + O2 (+O2) <=> HO2 (+O2) */
    a[9] = 4.66e+12;
    b[9] = 0.44;
    e[9] = 0;

    /*reaction 11: H + O2 (+H2O) <=> HO2 (+H2O) */
    a[10] = 9.06e+12;
    b[10] = 0.2;
    e[10] = 0;

    /*reaction 12: OH + OH (+M) <=> H2O2 (+M) */
    a[11] = 7.2e+13;
    b[11] = -0.37;
    e[11] = 0;

    /*reaction 13: OH + OH (+H2O) <=> H2O2 (+H2O) */
    a[12] = 7.2e+13;
    b[12] = -0.37;
    e[12] = 0;

    /*reaction 14: O + H2 <=> OH + H */
    a[13] = 50600;
    b[13] = 2.67;
    e[13] = 6290;

    /*reaction 15: H + O2 <=> OH + O */
    a[14] = 9.75e+13;
    b[14] = 0;
    e[14] = 14850;

    /*reaction 16: H2 + OH <=> H2O + H */
    a[15] = 1e+08;
    b[15] = 1.6;
    e[15] = 3300;

    /*reaction 17: OH + OH <=> H2O + O */
    a[16] = 35700;
    b[16] = 2.4;
    e[16] = -2110;

    /*reaction 18: HO2 + O <=> OH + O2 */
    a[17] = 1.63e+13;
    b[17] = 0;
    e[17] = -445;

    /*reaction 19: H + HO2 <=> OH + OH */
    a[18] = 1.7e+14;
    b[18] = 0;
    e[18] = 875;

    /*reaction 20: H + HO2 <=> H2O + O */
    a[19] = 3e+13;
    b[19] = 0;
    e[19] = 1720;

    /*reaction 21: H2 + O2 <=> H + HO2 */
    a[20] = 740000;
    b[20] = 2.4328;
    e[20] = 53500;

    /*reaction 22: H2 + O2 <=> OH + OH */
    a[21] = 2.04e+12;
    b[21] = 0.44;
    e[21] = 69155;

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    a[22] = 1e+13;
    b[22] = 0;
    e[22] = 0;

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    a[23] = 5.8e+13;
    b[23] = 0;
    e[23] = 3974;

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    a[24] = 1.03e+14;
    b[24] = 0;
    e[24] = 11040;

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    a[25] = 1.94e+11;
    b[25] = 0;
    e[25] = -1409;

    /*reaction 27: H2O2 + H <=> HO2 + H2 */
    a[26] = 1.7e+12;
    b[26] = 0;
    e[26] = 3755;

    /*reaction 28: H2O2 + H <=> H2O + OH */
    a[27] = 1e+13;
    b[27] = 0;
    e[27] = 3575;

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    a[28] = 2.8e+13;
    b[28] = 0;
    e[28] = 6400;

    /*reaction 30: H2O2 + OH <=> HO2 + H2O */
    a[29] = 2e+12;
    b[29] = 0;
    e[29] = 427;

    /*reaction 31: H2O2 + OH <=> HO2 + H2O */
    a[30] = 1.7e+18;
    b[30] = 0;
    e[30] = 29400;

    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double * T, double * C, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[10]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + H + M <=> H2 + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H + H + H2 <=> H2 + H2 */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + H <=> H2 + H */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + O + M <=> O2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + OH + M <=> H2O + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 (+M) <=> HO2 (+M) */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H + O2 (+AR) <=> HO2 (+AR) */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 (+O2) <=> HO2 (+O2) */
    eqcon[9] *= 1e+06; 

    /*reaction 11: H + O2 (+H2O) <=> HO2 (+H2O) */
    eqcon[10] *= 1e+06; 

    /*reaction 12: OH + OH (+M) <=> H2O2 (+M) */
    eqcon[11] *= 1e+06; 

    /*reaction 13: OH + OH (+H2O) <=> H2O2 (+H2O) */
    eqcon[12] *= 1e+06; 

    /*reaction 14: O + H2 <=> OH + H */
    /*eqcon[13] *= 1;  */

    /*reaction 15: H + O2 <=> OH + O */
    /*eqcon[14] *= 1;  */

    /*reaction 16: H2 + OH <=> H2O + H */
    /*eqcon[15] *= 1;  */

    /*reaction 17: OH + OH <=> H2O + O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + O <=> OH + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + HO2 <=> OH + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H + HO2 <=> H2O + O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2 + O2 <=> H + HO2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: H2 + O2 <=> OH + OH */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + H <=> HO2 + H2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + H <=> H2O + OH */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    /*eqcon[28] *= 1;  */

    /*reaction 30: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[29] *= 1;  */

    /*reaction 31: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[30] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[10]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + H + M <=> H2 + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H + H + H2 <=> H2 + H2 */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + H <=> H2 + H */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + O + M <=> O2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + OH + M <=> H2O + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 (+M) <=> HO2 (+M) */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H + O2 (+AR) <=> HO2 (+AR) */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 (+O2) <=> HO2 (+O2) */
    eqcon[9] *= 1e+06; 

    /*reaction 11: H + O2 (+H2O) <=> HO2 (+H2O) */
    eqcon[10] *= 1e+06; 

    /*reaction 12: OH + OH (+M) <=> H2O2 (+M) */
    eqcon[11] *= 1e+06; 

    /*reaction 13: OH + OH (+H2O) <=> H2O2 (+H2O) */
    eqcon[12] *= 1e+06; 

    /*reaction 14: O + H2 <=> OH + H */
    /*eqcon[13] *= 1;  */

    /*reaction 15: H + O2 <=> OH + O */
    /*eqcon[14] *= 1;  */

    /*reaction 16: H2 + OH <=> H2O + H */
    /*eqcon[15] *= 1;  */

    /*reaction 17: OH + OH <=> H2O + O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + O <=> OH + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + HO2 <=> OH + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H + HO2 <=> H2O + O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2 + O2 <=> H + HO2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: H2 + O2 <=> OH + OH */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + H <=> HO2 + H2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + H <=> H2O + OH */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    /*eqcon[28] *= 1;  */

    /*reaction 30: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[29] *= 1;  */

    /*reaction 31: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[30] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[10]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + H + M <=> H2 + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H + H + H2 <=> H2 + H2 */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + H <=> H2 + H */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + O + M <=> O2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + OH + M <=> H2O + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 (+M) <=> HO2 (+M) */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H + O2 (+AR) <=> HO2 (+AR) */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 (+O2) <=> HO2 (+O2) */
    eqcon[9] *= 1e+06; 

    /*reaction 11: H + O2 (+H2O) <=> HO2 (+H2O) */
    eqcon[10] *= 1e+06; 

    /*reaction 12: OH + OH (+M) <=> H2O2 (+M) */
    eqcon[11] *= 1e+06; 

    /*reaction 13: OH + OH (+H2O) <=> H2O2 (+H2O) */
    eqcon[12] *= 1e+06; 

    /*reaction 14: O + H2 <=> OH + H */
    /*eqcon[13] *= 1;  */

    /*reaction 15: H + O2 <=> OH + O */
    /*eqcon[14] *= 1;  */

    /*reaction 16: H2 + OH <=> H2O + H */
    /*eqcon[15] *= 1;  */

    /*reaction 17: OH + OH <=> H2O + O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + O <=> OH + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + HO2 <=> OH + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H + HO2 <=> H2O + O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2 + O2 <=> H + HO2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: H2 + O2 <=> OH + OH */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + H <=> HO2 + H2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + H <=> H2O + OH */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    /*eqcon[28] *= 1;  */

    /*reaction 30: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[29] *= 1;  */

    /*reaction 31: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[30] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[10]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + H + M <=> H2 + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H + H + H2 <=> H2 + H2 */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + H <=> H2 + H */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + O + M <=> O2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + OH + M <=> H2O + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 (+M) <=> HO2 (+M) */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H + O2 (+AR) <=> HO2 (+AR) */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 (+O2) <=> HO2 (+O2) */
    eqcon[9] *= 1e+06; 

    /*reaction 11: H + O2 (+H2O) <=> HO2 (+H2O) */
    eqcon[10] *= 1e+06; 

    /*reaction 12: OH + OH (+M) <=> H2O2 (+M) */
    eqcon[11] *= 1e+06; 

    /*reaction 13: OH + OH (+H2O) <=> H2O2 (+H2O) */
    eqcon[12] *= 1e+06; 

    /*reaction 14: O + H2 <=> OH + H */
    /*eqcon[13] *= 1;  */

    /*reaction 15: H + O2 <=> OH + O */
    /*eqcon[14] *= 1;  */

    /*reaction 16: H2 + OH <=> H2O + H */
    /*eqcon[15] *= 1;  */

    /*reaction 17: OH + OH <=> H2O + O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + O <=> OH + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + HO2 <=> OH + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H + HO2 <=> H2O + O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2 + O2 <=> H + HO2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: H2 + O2 <=> OH + OH */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + H <=> HO2 + H2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + H <=> H2O + OH */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    /*eqcon[28] *= 1;  */

    /*reaction 30: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[29] *= 1;  */

    /*reaction 31: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[30] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[10]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + H + M <=> H2 + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H + H + H2 <=> H2 + H2 */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + H <=> H2 + H */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + O + M <=> O2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + OH + M <=> H2O + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 (+M) <=> HO2 (+M) */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H + O2 (+AR) <=> HO2 (+AR) */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 (+O2) <=> HO2 (+O2) */
    eqcon[9] *= 1e+06; 

    /*reaction 11: H + O2 (+H2O) <=> HO2 (+H2O) */
    eqcon[10] *= 1e+06; 

    /*reaction 12: OH + OH (+M) <=> H2O2 (+M) */
    eqcon[11] *= 1e+06; 

    /*reaction 13: OH + OH (+H2O) <=> H2O2 (+H2O) */
    eqcon[12] *= 1e+06; 

    /*reaction 14: O + H2 <=> OH + H */
    /*eqcon[13] *= 1;  */

    /*reaction 15: H + O2 <=> OH + O */
    /*eqcon[14] *= 1;  */

    /*reaction 16: H2 + OH <=> H2O + H */
    /*eqcon[15] *= 1;  */

    /*reaction 17: OH + OH <=> H2O + O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + O <=> OH + O2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + HO2 <=> OH + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H + HO2 <=> H2O + O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2 + O2 <=> H + HO2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: H2 + O2 <=> OH + OH */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + H <=> HO2 + H2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + H <=> H2O + OH */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    /*eqcon[28] *= 1;  */

    /*reaction 30: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[29] *= 1;  */

    /*reaction 31: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[30] *= 1;  */
}

static double T_save = -1;
#ifdef BL_USE_OMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[31];
#ifdef BL_USE_OMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[31];
#ifdef BL_USE_OMP
#pragma omp threadprivate(Kc_save)
#endif

/*compute the production rate for each species */
void productionRate(double * wdot, double * sc, double T)
{
    double qdot;

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[10];                /*Gibbs free energy */
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

    double invT = 1.0 / tc[1];

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 10; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*zero out wdot */
    for (id = 0; id < 10; ++id) {
        wdot[id] = 0.0;
    }

    if (T != T_save)
    {
        T_save = T;

        k_f_save[0] = 1e-12 * 7e+17*exp(-1*tc[0]);
        k_f_save[1] = 1e-12 * 1e+17*exp(-0.6*tc[0]);
        k_f_save[2] = 1e-12 * 5.4e+18*exp(-1.3*tc[0]);
        k_f_save[3] = 1e-12 * 3.2e+15;
        k_f_save[4] = 1e-12 * 1e+17*exp(-1*tc[0]);
        k_f_save[5] = 1e-12 * 6.75e+18*exp(-1*tc[0]);
        k_f_save[6] = 1e-12 * 2.2e+22*exp(-2*tc[0]);
        k_f_save[7] = 1e-06 * 4.66e+12*exp(0.44*tc[0]);
        k_f_save[8] = 1e-06 * 4.66e+12*exp(0.44*tc[0]);
        k_f_save[9] = 1e-06 * 4.66e+12*exp(0.44*tc[0]);
        k_f_save[10] = 1e-06 * 9.06e+12*exp(0.2*tc[0]);
        k_f_save[11] = 1e-06 * 7.2e+13*exp(-0.37*tc[0]);
        k_f_save[12] = 1e-06 * 7.2e+13*exp(-0.37*tc[0]);
        k_f_save[13] = 1e-06 * 50600*exp(2.67*tc[0]-3165.2328279116868543*invT);
        k_f_save[14] = 1e-06 * 9.75e+13*exp(-7472.7674872000870891*invT);
        k_f_save[15] = 1e-06 * 1e+08*exp(1.6*tc[0]-1660.6149971555748834*invT);
        k_f_save[16] = 1e-06 * 35700*exp(2.4*tc[0]+1061.7871648479585929*invT);
        k_f_save[17] = 1e-06 * 1.63e+13*exp(+223.93141628310027613*invT);
        k_f_save[18] = 1e-06 * 1.7e+14*exp(-440.31458257912970566*invT);
        k_f_save[19] = 1e-06 * 3e+13*exp(-865.5326651841179455*invT);
        k_f_save[20] = 1e-06 * 740000*exp(2.4328*tc[0]-26922.091620552502718*invT);
        k_f_save[21] = 1e-06 * 2.04e+12*exp(0.44*tc[0]-34799.94852372539026*invT);
        k_f_save[22] = 1e-06 * 1e+13;
        k_f_save[23] = 1e-06 * 5.8e+13*exp(-1999.7830299079560064*invT);
        k_f_save[24] = 1e-06 * 1.03e+14*exp(-5555.5119904841048992*invT);
        k_f_save[25] = 1e-06 * 1.94e+11*exp(+709.03228211885004839*invT);
        k_f_save[26] = 1e-06 * 1.7e+12*exp(-1889.5785800967225896*invT);
        k_f_save[27] = 1e-06 * 1e+13*exp(-1798.9995802518728851*invT);
        k_f_save[28] = 1e-06 * 2.8e+13*exp(-3220.586661150206055*invT);
        k_f_save[29] = 1e-06 * 2e+12*exp(-214.8735162986153*invT);
        k_f_save[30] = 1e-06 * 1.7e+18*exp(-14794.569974658759747*invT);

        Kc_save[0] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0]) - (g_RT[1]));
        Kc_save[1] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[1]) - (g_RT[1] + g_RT[1]));
        Kc_save[2] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[9]) - (g_RT[1] + g_RT[9]));
        Kc_save[3] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[0]) - (g_RT[1] + g_RT[0]));
        Kc_save[4] = 1.0 / (refC) * exp((g_RT[2] + g_RT[2]) - (g_RT[3]));
        Kc_save[5] = 1.0 / (refC) * exp((g_RT[2] + g_RT[0]) - (g_RT[5]));
        Kc_save[6] = 1.0 / (refC) * exp((g_RT[0] + g_RT[5]) - (g_RT[4]));
        Kc_save[7] = 1.0 / (refC) * exp((g_RT[0] + g_RT[3]) - (g_RT[7]));
        Kc_save[8] = 1.0 / (refC) * exp((g_RT[0] + g_RT[3]) - (g_RT[7]));
        Kc_save[9] = 1.0 / (refC) * exp((g_RT[0] + g_RT[3]) - (g_RT[7]));
        Kc_save[10] = 1.0 / (refC) * exp((g_RT[0] + g_RT[3]) - (g_RT[7]));
        Kc_save[11] = 1.0 / (refC) * exp((g_RT[5] + g_RT[5]) - (g_RT[6]));
        Kc_save[12] = 1.0 / (refC) * exp((g_RT[5] + g_RT[5]) - (g_RT[6]));
        Kc_save[13] = exp((g_RT[2] + g_RT[1]) - (g_RT[5] + g_RT[0]));
        Kc_save[14] = exp((g_RT[0] + g_RT[3]) - (g_RT[5] + g_RT[2]));
        Kc_save[15] = exp((g_RT[1] + g_RT[5]) - (g_RT[4] + g_RT[0]));
        Kc_save[16] = exp((g_RT[5] + g_RT[5]) - (g_RT[4] + g_RT[2]));
        Kc_save[17] = exp((g_RT[7] + g_RT[2]) - (g_RT[5] + g_RT[3]));
        Kc_save[18] = exp((g_RT[0] + g_RT[7]) - (g_RT[5] + g_RT[5]));
        Kc_save[19] = exp((g_RT[0] + g_RT[7]) - (g_RT[4] + g_RT[2]));
        Kc_save[20] = exp((g_RT[1] + g_RT[3]) - (g_RT[0] + g_RT[7]));
        Kc_save[21] = exp((g_RT[1] + g_RT[3]) - (g_RT[5] + g_RT[5]));
        Kc_save[22] = exp((g_RT[7] + g_RT[5]) - (g_RT[4] + g_RT[3]));
        Kc_save[23] = exp((g_RT[7] + g_RT[5]) - (g_RT[4] + g_RT[3]));
        Kc_save[24] = exp((g_RT[7] + g_RT[7]) - (g_RT[6] + g_RT[3]));
        Kc_save[25] = exp((g_RT[7] + g_RT[7]) - (g_RT[6] + g_RT[3]));
        Kc_save[26] = exp((g_RT[6] + g_RT[0]) - (g_RT[7] + g_RT[1]));
        Kc_save[27] = exp((g_RT[6] + g_RT[0]) - (g_RT[4] + g_RT[5]));
        Kc_save[28] = exp((g_RT[6] + g_RT[2]) - (g_RT[7] + g_RT[5]));
        Kc_save[29] = exp((g_RT[6] + g_RT[5]) - (g_RT[7] + g_RT[4]));
        Kc_save[30] = exp((g_RT[6] + g_RT[5]) - (g_RT[7] + g_RT[4]));
    }

    /*reaction 1: H + H + M <=> H2 + M */
    phi_f = sc[0]*sc[0];
    alpha = mixture + -1*sc[1] + -1*sc[9] + -1*sc[0] + 13.3*sc[4];
    k_f = alpha * k_f_save[0];
    q_f = phi_f * k_f;
    phi_r = sc[1];
    Kc = Kc_save[0];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 2: H + H + H2 <=> H2 + H2 */
    phi_f = sc[0]*sc[0]*sc[1];
    k_f = k_f_save[1];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[1];
    Kc = Kc_save[1];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    phi_f = sc[0]*sc[0]*sc[9];
    k_f = k_f_save[2];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[9];
    Kc = Kc_save[2];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[9] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 4: H + H + H <=> H2 + H */
    phi_f = sc[0]*sc[0]*sc[0];
    k_f = k_f_save[3];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[0];
    Kc = Kc_save[3];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 5: O + O + M <=> O2 + M */
    phi_f = sc[2]*sc[2];
    alpha = mixture + 27.8*sc[2] + 7*sc[3] + sc[9] + 4*sc[4];
    k_f = alpha * k_f_save[4];
    q_f = phi_f * k_f;
    phi_r = sc[3];
    Kc = Kc_save[4];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 6: O + H + M <=> OH + M */
    phi_f = sc[2]*sc[0];
    alpha = mixture + 4*sc[4];
    k_f = alpha * k_f_save[5];
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = Kc_save[5];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 7: H + OH + M <=> H2O + M */
    phi_f = sc[0]*sc[5];
    alpha = mixture + 5.4*sc[4] + -0.62*sc[8];
    k_f = alpha * k_f_save[6];
    q_f = phi_f * k_f;
    phi_r = sc[4];
    Kc = Kc_save[6];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 8: H + O2 (+M) <=> HO2 (+M) */
    phi_f = sc[0]*sc[3];
    alpha = mixture + -1*sc[8] + -1*sc[4] + -1*sc[3] + 0.5*sc[1];
    k_f = k_f_save[7];
    redP = 1e-12 * alpha / k_f * 5.7e+19*exp(-1.4*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-100000))+ (0.5*exp(T/-10)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[7];
    Kc = Kc_save[7];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 9: H + O2 (+AR) <=> HO2 (+AR) */
    phi_f = sc[0]*sc[3];
    alpha = sc[8];
    k_f = k_f_save[8];
    redP = 1e-12 * alpha / k_f * 7.43e+18*exp(-1.2*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-10))+ (0.5*exp(T/-100000)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[7];
    Kc = Kc_save[8];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 10: H + O2 (+O2) <=> HO2 (+O2) */
    phi_f = sc[0]*sc[3];
    alpha = sc[3];
    k_f = k_f_save[9];
    redP = 1e-12 * alpha / k_f * 5.69e+18*exp(-1.094*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-100000))+ (0.5*exp(T/-10)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[7];
    Kc = Kc_save[9];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 11: H + O2 (+H2O) <=> HO2 (+H2O) */
    phi_f = sc[0]*sc[3];
    alpha = sc[4];
    k_f = k_f_save[10];
    redP = 1e-12 * alpha / k_f * 3.67e+19*exp(-1*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2*exp(T/-10))+ (0.8*exp(T/-100000)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[7];
    Kc = Kc_save[10];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 12: OH + OH (+M) <=> H2O2 (+M) */
    phi_f = sc[5]*sc[5];
    alpha = mixture + -1*sc[4];
    k_f = k_f_save[11];
    redP = 1e-12 * alpha / k_f * 2.2e+19*exp(-0.76*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-100000))+ (0.5*exp(T/-10)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = Kc_save[11];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 13: OH + OH (+H2O) <=> H2O2 (+H2O) */
    phi_f = sc[5]*sc[5];
    alpha = sc[4];
    k_f = k_f_save[12];
    redP = 1e-12 * alpha / k_f * 1.45e+18;
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = Kc_save[12];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 14: O + H2 <=> OH + H */
    phi_f = sc[2]*sc[1];
    k_f = k_f_save[13];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[0];
    Kc = Kc_save[13];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 15: H + O2 <=> OH + O */
    phi_f = sc[0]*sc[3];
    k_f = k_f_save[14];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[2];
    Kc = Kc_save[14];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 16: H2 + OH <=> H2O + H */
    phi_f = sc[1]*sc[5];
    k_f = k_f_save[15];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[0];
    Kc = Kc_save[15];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 17: OH + OH <=> H2O + O */
    phi_f = sc[5]*sc[5];
    k_f = k_f_save[16];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[2];
    Kc = Kc_save[16];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 18: HO2 + O <=> OH + O2 */
    phi_f = sc[7]*sc[2];
    k_f = k_f_save[17];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[3];
    Kc = Kc_save[17];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 19: H + HO2 <=> OH + OH */
    phi_f = sc[0]*sc[7];
    k_f = k_f_save[18];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[5];
    Kc = Kc_save[18];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 20: H + HO2 <=> H2O + O */
    phi_f = sc[0]*sc[7];
    k_f = k_f_save[19];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[2];
    Kc = Kc_save[19];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 21: H2 + O2 <=> H + HO2 */
    phi_f = sc[1]*sc[3];
    k_f = k_f_save[20];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[7];
    Kc = Kc_save[20];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 22: H2 + O2 <=> OH + OH */
    phi_f = sc[1]*sc[3];
    k_f = k_f_save[21];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[5];
    Kc = Kc_save[21];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    phi_f = sc[7]*sc[5];
    k_f = k_f_save[22];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = Kc_save[22];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    phi_f = sc[7]*sc[5];
    k_f = k_f_save[23];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = Kc_save[23];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[7]*sc[7];
    k_f = k_f_save[24];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = Kc_save[24];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[7]*sc[7];
    k_f = k_f_save[25];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = Kc_save[25];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 27: H2O2 + H <=> HO2 + H2 */
    phi_f = sc[6]*sc[0];
    k_f = k_f_save[26];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[1];
    Kc = Kc_save[26];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 28: H2O2 + H <=> H2O + OH */
    phi_f = sc[6]*sc[0];
    k_f = k_f_save[27];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[5];
    Kc = Kc_save[27];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    phi_f = sc[6]*sc[2];
    k_f = k_f_save[28];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[5];
    Kc = Kc_save[28];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 30: H2O2 + OH <=> HO2 + H2O */
    phi_f = sc[6]*sc[5];
    k_f = k_f_save[29];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[4];
    Kc = Kc_save[29];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 31: H2O2 + OH <=> HO2 + H2O */
    phi_f = sc[6]*sc[5];
    k_f = k_f_save[30];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[4];
    Kc = Kc_save[30];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[4] += 1 * qdot;

    return;
}


/*compute the progress rate for each reaction */
void progressRate(double * qdot, double * sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[10];                /*Gibbs free energy */
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

    double invT = 1.0 / tc[1];

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 10; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    if (T != T_save)
    {
        T_save = T;

        k_f_save[0] = 1e-12 * 7e+17*exp(-1*tc[0]);
        k_f_save[1] = 1e-12 * 1e+17*exp(-0.6*tc[0]);
        k_f_save[2] = 1e-12 * 5.4e+18*exp(-1.3*tc[0]);
        k_f_save[3] = 1e-12 * 3.2e+15;
        k_f_save[4] = 1e-12 * 1e+17*exp(-1*tc[0]);
        k_f_save[5] = 1e-12 * 6.75e+18*exp(-1*tc[0]);
        k_f_save[6] = 1e-12 * 2.2e+22*exp(-2*tc[0]);
        k_f_save[7] = 1e-06 * 4.66e+12*exp(0.44*tc[0]);
        k_f_save[8] = 1e-06 * 4.66e+12*exp(0.44*tc[0]);
        k_f_save[9] = 1e-06 * 4.66e+12*exp(0.44*tc[0]);
        k_f_save[10] = 1e-06 * 9.06e+12*exp(0.2*tc[0]);
        k_f_save[11] = 1e-06 * 7.2e+13*exp(-0.37*tc[0]);
        k_f_save[12] = 1e-06 * 7.2e+13*exp(-0.37*tc[0]);
        k_f_save[13] = 1e-06 * 50600*exp(2.67*tc[0]-3165.2328279116868543*invT);
        k_f_save[14] = 1e-06 * 9.75e+13*exp(-7472.7674872000870891*invT);
        k_f_save[15] = 1e-06 * 1e+08*exp(1.6*tc[0]-1660.6149971555748834*invT);
        k_f_save[16] = 1e-06 * 35700*exp(2.4*tc[0]+1061.7871648479585929*invT);
        k_f_save[17] = 1e-06 * 1.63e+13*exp(+223.93141628310027613*invT);
        k_f_save[18] = 1e-06 * 1.7e+14*exp(-440.31458257912970566*invT);
        k_f_save[19] = 1e-06 * 3e+13*exp(-865.5326651841179455*invT);
        k_f_save[20] = 1e-06 * 740000*exp(2.4328*tc[0]-26922.091620552502718*invT);
        k_f_save[21] = 1e-06 * 2.04e+12*exp(0.44*tc[0]-34799.94852372539026*invT);
        k_f_save[22] = 1e-06 * 1e+13;
        k_f_save[23] = 1e-06 * 5.8e+13*exp(-1999.7830299079560064*invT);
        k_f_save[24] = 1e-06 * 1.03e+14*exp(-5555.5119904841048992*invT);
        k_f_save[25] = 1e-06 * 1.94e+11*exp(+709.03228211885004839*invT);
        k_f_save[26] = 1e-06 * 1.7e+12*exp(-1889.5785800967225896*invT);
        k_f_save[27] = 1e-06 * 1e+13*exp(-1798.9995802518728851*invT);
        k_f_save[28] = 1e-06 * 2.8e+13*exp(-3220.586661150206055*invT);
        k_f_save[29] = 1e-06 * 2e+12*exp(-214.8735162986153*invT);
        k_f_save[30] = 1e-06 * 1.7e+18*exp(-14794.569974658759747*invT);

        Kc_save[0] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0]) - (g_RT[1]));
        Kc_save[1] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[1]) - (g_RT[1] + g_RT[1]));
        Kc_save[2] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[9]) - (g_RT[1] + g_RT[9]));
        Kc_save[3] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[0]) - (g_RT[1] + g_RT[0]));
        Kc_save[4] = 1.0 / (refC) * exp((g_RT[2] + g_RT[2]) - (g_RT[3]));
        Kc_save[5] = 1.0 / (refC) * exp((g_RT[2] + g_RT[0]) - (g_RT[5]));
        Kc_save[6] = 1.0 / (refC) * exp((g_RT[0] + g_RT[5]) - (g_RT[4]));
        Kc_save[7] = 1.0 / (refC) * exp((g_RT[0] + g_RT[3]) - (g_RT[7]));
        Kc_save[8] = 1.0 / (refC) * exp((g_RT[0] + g_RT[3]) - (g_RT[7]));
        Kc_save[9] = 1.0 / (refC) * exp((g_RT[0] + g_RT[3]) - (g_RT[7]));
        Kc_save[10] = 1.0 / (refC) * exp((g_RT[0] + g_RT[3]) - (g_RT[7]));
        Kc_save[11] = 1.0 / (refC) * exp((g_RT[5] + g_RT[5]) - (g_RT[6]));
        Kc_save[12] = 1.0 / (refC) * exp((g_RT[5] + g_RT[5]) - (g_RT[6]));
        Kc_save[13] = exp((g_RT[2] + g_RT[1]) - (g_RT[5] + g_RT[0]));
        Kc_save[14] = exp((g_RT[0] + g_RT[3]) - (g_RT[5] + g_RT[2]));
        Kc_save[15] = exp((g_RT[1] + g_RT[5]) - (g_RT[4] + g_RT[0]));
        Kc_save[16] = exp((g_RT[5] + g_RT[5]) - (g_RT[4] + g_RT[2]));
        Kc_save[17] = exp((g_RT[7] + g_RT[2]) - (g_RT[5] + g_RT[3]));
        Kc_save[18] = exp((g_RT[0] + g_RT[7]) - (g_RT[5] + g_RT[5]));
        Kc_save[19] = exp((g_RT[0] + g_RT[7]) - (g_RT[4] + g_RT[2]));
        Kc_save[20] = exp((g_RT[1] + g_RT[3]) - (g_RT[0] + g_RT[7]));
        Kc_save[21] = exp((g_RT[1] + g_RT[3]) - (g_RT[5] + g_RT[5]));
        Kc_save[22] = exp((g_RT[7] + g_RT[5]) - (g_RT[4] + g_RT[3]));
        Kc_save[23] = exp((g_RT[7] + g_RT[5]) - (g_RT[4] + g_RT[3]));
        Kc_save[24] = exp((g_RT[7] + g_RT[7]) - (g_RT[6] + g_RT[3]));
        Kc_save[25] = exp((g_RT[7] + g_RT[7]) - (g_RT[6] + g_RT[3]));
        Kc_save[26] = exp((g_RT[6] + g_RT[0]) - (g_RT[7] + g_RT[1]));
        Kc_save[27] = exp((g_RT[6] + g_RT[0]) - (g_RT[4] + g_RT[5]));
        Kc_save[28] = exp((g_RT[6] + g_RT[2]) - (g_RT[7] + g_RT[5]));
        Kc_save[29] = exp((g_RT[6] + g_RT[5]) - (g_RT[7] + g_RT[4]));
        Kc_save[30] = exp((g_RT[6] + g_RT[5]) - (g_RT[7] + g_RT[4]));
    }

    /*reaction 1: H + H + M <=> H2 + M */
    phi_f = sc[0]*sc[0];
    alpha = mixture + -1*sc[1] + -1*sc[9] + -1*sc[0] + 13.3*sc[4];
    k_f = alpha * k_f_save[0];
    q_f = phi_f * k_f;
    phi_r = sc[1];
    Kc = Kc_save[0];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[0] = q_f - q_r;

    /*reaction 2: H + H + H2 <=> H2 + H2 */
    phi_f = sc[0]*sc[0]*sc[1];
    k_f = k_f_save[1];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[1];
    Kc = Kc_save[1];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[1] = q_f - q_r;

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    phi_f = sc[0]*sc[0]*sc[9];
    k_f = k_f_save[2];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[9];
    Kc = Kc_save[2];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[2] = q_f - q_r;

    /*reaction 4: H + H + H <=> H2 + H */
    phi_f = sc[0]*sc[0]*sc[0];
    k_f = k_f_save[3];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[0];
    Kc = Kc_save[3];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[3] = q_f - q_r;

    /*reaction 5: O + O + M <=> O2 + M */
    phi_f = sc[2]*sc[2];
    alpha = mixture + 27.8*sc[2] + 7*sc[3] + sc[9] + 4*sc[4];
    k_f = alpha * k_f_save[4];
    q_f = phi_f * k_f;
    phi_r = sc[3];
    Kc = Kc_save[4];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[4] = q_f - q_r;

    /*reaction 6: O + H + M <=> OH + M */
    phi_f = sc[2]*sc[0];
    alpha = mixture + 4*sc[4];
    k_f = alpha * k_f_save[5];
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = Kc_save[5];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[5] = q_f - q_r;

    /*reaction 7: H + OH + M <=> H2O + M */
    phi_f = sc[0]*sc[5];
    alpha = mixture + 5.4*sc[4] + -0.62*sc[8];
    k_f = alpha * k_f_save[6];
    q_f = phi_f * k_f;
    phi_r = sc[4];
    Kc = Kc_save[6];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[6] = q_f - q_r;

    /*reaction 8: H + O2 (+M) <=> HO2 (+M) */
    phi_f = sc[0]*sc[3];
    alpha = mixture + -1*sc[8] + -1*sc[4] + -1*sc[3] + 0.5*sc[1];
    k_f = k_f_save[7];
    redP = 1e-12 * alpha / k_f * 5.7e+19*exp(-1.4*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-100000))+ (0.5*exp(T/-10)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[7];
    Kc = Kc_save[7];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[7] = q_f - q_r;

    /*reaction 9: H + O2 (+AR) <=> HO2 (+AR) */
    phi_f = sc[0]*sc[3];
    alpha = sc[8];
    k_f = k_f_save[8];
    redP = 1e-12 * alpha / k_f * 7.43e+18*exp(-1.2*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-10))+ (0.5*exp(T/-100000)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[7];
    Kc = Kc_save[8];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[8] = q_f - q_r;

    /*reaction 10: H + O2 (+O2) <=> HO2 (+O2) */
    phi_f = sc[0]*sc[3];
    alpha = sc[3];
    k_f = k_f_save[9];
    redP = 1e-12 * alpha / k_f * 5.69e+18*exp(-1.094*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-100000))+ (0.5*exp(T/-10)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[7];
    Kc = Kc_save[9];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[9] = q_f - q_r;

    /*reaction 11: H + O2 (+H2O) <=> HO2 (+H2O) */
    phi_f = sc[0]*sc[3];
    alpha = sc[4];
    k_f = k_f_save[10];
    redP = 1e-12 * alpha / k_f * 3.67e+19*exp(-1*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2*exp(T/-10))+ (0.8*exp(T/-100000)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[7];
    Kc = Kc_save[10];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[10] = q_f - q_r;

    /*reaction 12: OH + OH (+M) <=> H2O2 (+M) */
    phi_f = sc[5]*sc[5];
    alpha = mixture + -1*sc[4];
    k_f = k_f_save[11];
    redP = 1e-12 * alpha / k_f * 2.2e+19*exp(-0.76*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-100000))+ (0.5*exp(T/-10)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = Kc_save[11];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[11] = q_f - q_r;

    /*reaction 13: OH + OH (+H2O) <=> H2O2 (+H2O) */
    phi_f = sc[5]*sc[5];
    alpha = sc[4];
    k_f = k_f_save[12];
    redP = 1e-12 * alpha / k_f * 1.45e+18;
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = Kc_save[12];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[12] = q_f - q_r;

    /*reaction 14: O + H2 <=> OH + H */
    phi_f = sc[2]*sc[1];
    k_f = k_f_save[13];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[0];
    Kc = Kc_save[13];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[13] = q_f - q_r;

    /*reaction 15: H + O2 <=> OH + O */
    phi_f = sc[0]*sc[3];
    k_f = k_f_save[14];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[2];
    Kc = Kc_save[14];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[14] = q_f - q_r;

    /*reaction 16: H2 + OH <=> H2O + H */
    phi_f = sc[1]*sc[5];
    k_f = k_f_save[15];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[0];
    Kc = Kc_save[15];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[15] = q_f - q_r;

    /*reaction 17: OH + OH <=> H2O + O */
    phi_f = sc[5]*sc[5];
    k_f = k_f_save[16];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[2];
    Kc = Kc_save[16];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[16] = q_f - q_r;

    /*reaction 18: HO2 + O <=> OH + O2 */
    phi_f = sc[7]*sc[2];
    k_f = k_f_save[17];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[3];
    Kc = Kc_save[17];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[17] = q_f - q_r;

    /*reaction 19: H + HO2 <=> OH + OH */
    phi_f = sc[0]*sc[7];
    k_f = k_f_save[18];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[5];
    Kc = Kc_save[18];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[18] = q_f - q_r;

    /*reaction 20: H + HO2 <=> H2O + O */
    phi_f = sc[0]*sc[7];
    k_f = k_f_save[19];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[2];
    Kc = Kc_save[19];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[19] = q_f - q_r;

    /*reaction 21: H2 + O2 <=> H + HO2 */
    phi_f = sc[1]*sc[3];
    k_f = k_f_save[20];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[7];
    Kc = Kc_save[20];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[20] = q_f - q_r;

    /*reaction 22: H2 + O2 <=> OH + OH */
    phi_f = sc[1]*sc[3];
    k_f = k_f_save[21];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[5];
    Kc = Kc_save[21];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[21] = q_f - q_r;

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    phi_f = sc[7]*sc[5];
    k_f = k_f_save[22];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = Kc_save[22];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[22] = q_f - q_r;

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    phi_f = sc[7]*sc[5];
    k_f = k_f_save[23];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = Kc_save[23];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[23] = q_f - q_r;

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[7]*sc[7];
    k_f = k_f_save[24];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = Kc_save[24];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[24] = q_f - q_r;

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[7]*sc[7];
    k_f = k_f_save[25];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = Kc_save[25];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[25] = q_f - q_r;

    /*reaction 27: H2O2 + H <=> HO2 + H2 */
    phi_f = sc[6]*sc[0];
    k_f = k_f_save[26];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[1];
    Kc = Kc_save[26];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[26] = q_f - q_r;

    /*reaction 28: H2O2 + H <=> H2O + OH */
    phi_f = sc[6]*sc[0];
    k_f = k_f_save[27];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[5];
    Kc = Kc_save[27];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[27] = q_f - q_r;

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    phi_f = sc[6]*sc[2];
    k_f = k_f_save[28];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[5];
    Kc = Kc_save[28];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[28] = q_f - q_r;

    /*reaction 30: H2O2 + OH <=> HO2 + H2O */
    phi_f = sc[6]*sc[5];
    k_f = k_f_save[29];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[4];
    Kc = Kc_save[29];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[29] = q_f - q_r;

    /*reaction 31: H2O2 + OH <=> HO2 + H2O */
    phi_f = sc[6]*sc[5];
    k_f = k_f_save[30];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[4];
    Kc = Kc_save[30];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[30] = q_f - q_r;

    return;
}


/*compute the progress rate for each reaction */
void progressRateFR(double * q_f, double * q_r, double * sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[10];                /*Gibbs free energy */
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

    double invT = 1.0 / tc[1];

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 10; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    if (T != T_save)
    {
        T_save = T;

        k_f_save[0] = 1e-12 * 7e+17*exp(-1*tc[0]);
        k_f_save[1] = 1e-12 * 1e+17*exp(-0.6*tc[0]);
        k_f_save[2] = 1e-12 * 5.4e+18*exp(-1.3*tc[0]);
        k_f_save[3] = 1e-12 * 3.2e+15;
        k_f_save[4] = 1e-12 * 1e+17*exp(-1*tc[0]);
        k_f_save[5] = 1e-12 * 6.75e+18*exp(-1*tc[0]);
        k_f_save[6] = 1e-12 * 2.2e+22*exp(-2*tc[0]);
        k_f_save[7] = 1e-06 * 4.66e+12*exp(0.44*tc[0]);
        k_f_save[8] = 1e-06 * 4.66e+12*exp(0.44*tc[0]);
        k_f_save[9] = 1e-06 * 4.66e+12*exp(0.44*tc[0]);
        k_f_save[10] = 1e-06 * 9.06e+12*exp(0.2*tc[0]);
        k_f_save[11] = 1e-06 * 7.2e+13*exp(-0.37*tc[0]);
        k_f_save[12] = 1e-06 * 7.2e+13*exp(-0.37*tc[0]);
        k_f_save[13] = 1e-06 * 50600*exp(2.67*tc[0]-3165.2328279116868543*invT);
        k_f_save[14] = 1e-06 * 9.75e+13*exp(-7472.7674872000870891*invT);
        k_f_save[15] = 1e-06 * 1e+08*exp(1.6*tc[0]-1660.6149971555748834*invT);
        k_f_save[16] = 1e-06 * 35700*exp(2.4*tc[0]+1061.7871648479585929*invT);
        k_f_save[17] = 1e-06 * 1.63e+13*exp(+223.93141628310027613*invT);
        k_f_save[18] = 1e-06 * 1.7e+14*exp(-440.31458257912970566*invT);
        k_f_save[19] = 1e-06 * 3e+13*exp(-865.5326651841179455*invT);
        k_f_save[20] = 1e-06 * 740000*exp(2.4328*tc[0]-26922.091620552502718*invT);
        k_f_save[21] = 1e-06 * 2.04e+12*exp(0.44*tc[0]-34799.94852372539026*invT);
        k_f_save[22] = 1e-06 * 1e+13;
        k_f_save[23] = 1e-06 * 5.8e+13*exp(-1999.7830299079560064*invT);
        k_f_save[24] = 1e-06 * 1.03e+14*exp(-5555.5119904841048992*invT);
        k_f_save[25] = 1e-06 * 1.94e+11*exp(+709.03228211885004839*invT);
        k_f_save[26] = 1e-06 * 1.7e+12*exp(-1889.5785800967225896*invT);
        k_f_save[27] = 1e-06 * 1e+13*exp(-1798.9995802518728851*invT);
        k_f_save[28] = 1e-06 * 2.8e+13*exp(-3220.586661150206055*invT);
        k_f_save[29] = 1e-06 * 2e+12*exp(-214.8735162986153*invT);
        k_f_save[30] = 1e-06 * 1.7e+18*exp(-14794.569974658759747*invT);

        Kc_save[0] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0]) - (g_RT[1]));
        Kc_save[1] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[1]) - (g_RT[1] + g_RT[1]));
        Kc_save[2] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[9]) - (g_RT[1] + g_RT[9]));
        Kc_save[3] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[0]) - (g_RT[1] + g_RT[0]));
        Kc_save[4] = 1.0 / (refC) * exp((g_RT[2] + g_RT[2]) - (g_RT[3]));
        Kc_save[5] = 1.0 / (refC) * exp((g_RT[2] + g_RT[0]) - (g_RT[5]));
        Kc_save[6] = 1.0 / (refC) * exp((g_RT[0] + g_RT[5]) - (g_RT[4]));
        Kc_save[7] = 1.0 / (refC) * exp((g_RT[0] + g_RT[3]) - (g_RT[7]));
        Kc_save[8] = 1.0 / (refC) * exp((g_RT[0] + g_RT[3]) - (g_RT[7]));
        Kc_save[9] = 1.0 / (refC) * exp((g_RT[0] + g_RT[3]) - (g_RT[7]));
        Kc_save[10] = 1.0 / (refC) * exp((g_RT[0] + g_RT[3]) - (g_RT[7]));
        Kc_save[11] = 1.0 / (refC) * exp((g_RT[5] + g_RT[5]) - (g_RT[6]));
        Kc_save[12] = 1.0 / (refC) * exp((g_RT[5] + g_RT[5]) - (g_RT[6]));
        Kc_save[13] = exp((g_RT[2] + g_RT[1]) - (g_RT[5] + g_RT[0]));
        Kc_save[14] = exp((g_RT[0] + g_RT[3]) - (g_RT[5] + g_RT[2]));
        Kc_save[15] = exp((g_RT[1] + g_RT[5]) - (g_RT[4] + g_RT[0]));
        Kc_save[16] = exp((g_RT[5] + g_RT[5]) - (g_RT[4] + g_RT[2]));
        Kc_save[17] = exp((g_RT[7] + g_RT[2]) - (g_RT[5] + g_RT[3]));
        Kc_save[18] = exp((g_RT[0] + g_RT[7]) - (g_RT[5] + g_RT[5]));
        Kc_save[19] = exp((g_RT[0] + g_RT[7]) - (g_RT[4] + g_RT[2]));
        Kc_save[20] = exp((g_RT[1] + g_RT[3]) - (g_RT[0] + g_RT[7]));
        Kc_save[21] = exp((g_RT[1] + g_RT[3]) - (g_RT[5] + g_RT[5]));
        Kc_save[22] = exp((g_RT[7] + g_RT[5]) - (g_RT[4] + g_RT[3]));
        Kc_save[23] = exp((g_RT[7] + g_RT[5]) - (g_RT[4] + g_RT[3]));
        Kc_save[24] = exp((g_RT[7] + g_RT[7]) - (g_RT[6] + g_RT[3]));
        Kc_save[25] = exp((g_RT[7] + g_RT[7]) - (g_RT[6] + g_RT[3]));
        Kc_save[26] = exp((g_RT[6] + g_RT[0]) - (g_RT[7] + g_RT[1]));
        Kc_save[27] = exp((g_RT[6] + g_RT[0]) - (g_RT[4] + g_RT[5]));
        Kc_save[28] = exp((g_RT[6] + g_RT[2]) - (g_RT[7] + g_RT[5]));
        Kc_save[29] = exp((g_RT[6] + g_RT[5]) - (g_RT[7] + g_RT[4]));
        Kc_save[30] = exp((g_RT[6] + g_RT[5]) - (g_RT[7] + g_RT[4]));
    }

    /*reaction 1: H + H + M <=> H2 + M */
    phi_f = sc[0]*sc[0];
    alpha = mixture + -1*sc[1] + -1*sc[9] + -1*sc[0] + 13.3*sc[4];
    k_f = alpha * k_f_save[0];
    q_f[0] = phi_f * k_f;
    phi_r = sc[1];
    Kc = Kc_save[0];
    k_r = k_f / Kc;
    q_r[0] = phi_r * k_r;

    /*reaction 2: H + H + H2 <=> H2 + H2 */
    phi_f = sc[0]*sc[0]*sc[1];
    k_f = k_f_save[1];
    q_f[1] = phi_f * k_f;
    phi_r = sc[1]*sc[1];
    Kc = Kc_save[1];
    k_r = k_f / Kc;
    q_r[1] = phi_r * k_r;

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    phi_f = sc[0]*sc[0]*sc[9];
    k_f = k_f_save[2];
    q_f[2] = phi_f * k_f;
    phi_r = sc[1]*sc[9];
    Kc = Kc_save[2];
    k_r = k_f / Kc;
    q_r[2] = phi_r * k_r;

    /*reaction 4: H + H + H <=> H2 + H */
    phi_f = sc[0]*sc[0]*sc[0];
    k_f = k_f_save[3];
    q_f[3] = phi_f * k_f;
    phi_r = sc[1]*sc[0];
    Kc = Kc_save[3];
    k_r = k_f / Kc;
    q_r[3] = phi_r * k_r;

    /*reaction 5: O + O + M <=> O2 + M */
    phi_f = sc[2]*sc[2];
    alpha = mixture + 27.8*sc[2] + 7*sc[3] + sc[9] + 4*sc[4];
    k_f = alpha * k_f_save[4];
    q_f[4] = phi_f * k_f;
    phi_r = sc[3];
    Kc = Kc_save[4];
    k_r = k_f / Kc;
    q_r[4] = phi_r * k_r;

    /*reaction 6: O + H + M <=> OH + M */
    phi_f = sc[2]*sc[0];
    alpha = mixture + 4*sc[4];
    k_f = alpha * k_f_save[5];
    q_f[5] = phi_f * k_f;
    phi_r = sc[5];
    Kc = Kc_save[5];
    k_r = k_f / Kc;
    q_r[5] = phi_r * k_r;

    /*reaction 7: H + OH + M <=> H2O + M */
    phi_f = sc[0]*sc[5];
    alpha = mixture + 5.4*sc[4] + -0.62*sc[8];
    k_f = alpha * k_f_save[6];
    q_f[6] = phi_f * k_f;
    phi_r = sc[4];
    Kc = Kc_save[6];
    k_r = k_f / Kc;
    q_r[6] = phi_r * k_r;

    /*reaction 8: H + O2 (+M) <=> HO2 (+M) */
    phi_f = sc[0]*sc[3];
    alpha = mixture + -1*sc[8] + -1*sc[4] + -1*sc[3] + 0.5*sc[1];
    k_f = k_f_save[7];
    redP = 1e-12 * alpha / k_f * 5.7e+19*exp(-1.4*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-100000))+ (0.5*exp(T/-10)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[7] = phi_f * k_f;
    phi_r = sc[7];
    Kc = Kc_save[7];
    k_r = k_f / Kc;
    q_r[7] = phi_r * k_r;

    /*reaction 9: H + O2 (+AR) <=> HO2 (+AR) */
    phi_f = sc[0]*sc[3];
    alpha = sc[8];
    k_f = k_f_save[8];
    redP = 1e-12 * alpha / k_f * 7.43e+18*exp(-1.2*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-10))+ (0.5*exp(T/-100000)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[8] = phi_f * k_f;
    phi_r = sc[7];
    Kc = Kc_save[8];
    k_r = k_f / Kc;
    q_r[8] = phi_r * k_r;

    /*reaction 10: H + O2 (+O2) <=> HO2 (+O2) */
    phi_f = sc[0]*sc[3];
    alpha = sc[3];
    k_f = k_f_save[9];
    redP = 1e-12 * alpha / k_f * 5.69e+18*exp(-1.094*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-100000))+ (0.5*exp(T/-10)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[9] = phi_f * k_f;
    phi_r = sc[7];
    Kc = Kc_save[9];
    k_r = k_f / Kc;
    q_r[9] = phi_r * k_r;

    /*reaction 11: H + O2 (+H2O) <=> HO2 (+H2O) */
    phi_f = sc[0]*sc[3];
    alpha = sc[4];
    k_f = k_f_save[10];
    redP = 1e-12 * alpha / k_f * 3.67e+19*exp(-1*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2*exp(T/-10))+ (0.8*exp(T/-100000)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[10] = phi_f * k_f;
    phi_r = sc[7];
    Kc = Kc_save[10];
    k_r = k_f / Kc;
    q_r[10] = phi_r * k_r;

    /*reaction 12: OH + OH (+M) <=> H2O2 (+M) */
    phi_f = sc[5]*sc[5];
    alpha = mixture + -1*sc[4];
    k_f = k_f_save[11];
    redP = 1e-12 * alpha / k_f * 2.2e+19*exp(-0.76*tc[0]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-100000))+ (0.5*exp(T/-10)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[11] = phi_f * k_f;
    phi_r = sc[6];
    Kc = Kc_save[11];
    k_r = k_f / Kc;
    q_r[11] = phi_r * k_r;

    /*reaction 13: OH + OH (+H2O) <=> H2O2 (+H2O) */
    phi_f = sc[5]*sc[5];
    alpha = sc[4];
    k_f = k_f_save[12];
    redP = 1e-12 * alpha / k_f * 1.45e+18;
    F = redP / (1 + redP);
    k_f *= F;
    q_f[12] = phi_f * k_f;
    phi_r = sc[6];
    Kc = Kc_save[12];
    k_r = k_f / Kc;
    q_r[12] = phi_r * k_r;

    /*reaction 14: O + H2 <=> OH + H */
    phi_f = sc[2]*sc[1];
    k_f = k_f_save[13];
    q_f[13] = phi_f * k_f;
    phi_r = sc[5]*sc[0];
    Kc = Kc_save[13];
    k_r = k_f / Kc;
    q_r[13] = phi_r * k_r;

    /*reaction 15: H + O2 <=> OH + O */
    phi_f = sc[0]*sc[3];
    k_f = k_f_save[14];
    q_f[14] = phi_f * k_f;
    phi_r = sc[5]*sc[2];
    Kc = Kc_save[14];
    k_r = k_f / Kc;
    q_r[14] = phi_r * k_r;

    /*reaction 16: H2 + OH <=> H2O + H */
    phi_f = sc[1]*sc[5];
    k_f = k_f_save[15];
    q_f[15] = phi_f * k_f;
    phi_r = sc[4]*sc[0];
    Kc = Kc_save[15];
    k_r = k_f / Kc;
    q_r[15] = phi_r * k_r;

    /*reaction 17: OH + OH <=> H2O + O */
    phi_f = sc[5]*sc[5];
    k_f = k_f_save[16];
    q_f[16] = phi_f * k_f;
    phi_r = sc[4]*sc[2];
    Kc = Kc_save[16];
    k_r = k_f / Kc;
    q_r[16] = phi_r * k_r;

    /*reaction 18: HO2 + O <=> OH + O2 */
    phi_f = sc[7]*sc[2];
    k_f = k_f_save[17];
    q_f[17] = phi_f * k_f;
    phi_r = sc[5]*sc[3];
    Kc = Kc_save[17];
    k_r = k_f / Kc;
    q_r[17] = phi_r * k_r;

    /*reaction 19: H + HO2 <=> OH + OH */
    phi_f = sc[0]*sc[7];
    k_f = k_f_save[18];
    q_f[18] = phi_f * k_f;
    phi_r = sc[5]*sc[5];
    Kc = Kc_save[18];
    k_r = k_f / Kc;
    q_r[18] = phi_r * k_r;

    /*reaction 20: H + HO2 <=> H2O + O */
    phi_f = sc[0]*sc[7];
    k_f = k_f_save[19];
    q_f[19] = phi_f * k_f;
    phi_r = sc[4]*sc[2];
    Kc = Kc_save[19];
    k_r = k_f / Kc;
    q_r[19] = phi_r * k_r;

    /*reaction 21: H2 + O2 <=> H + HO2 */
    phi_f = sc[1]*sc[3];
    k_f = k_f_save[20];
    q_f[20] = phi_f * k_f;
    phi_r = sc[0]*sc[7];
    Kc = Kc_save[20];
    k_r = k_f / Kc;
    q_r[20] = phi_r * k_r;

    /*reaction 22: H2 + O2 <=> OH + OH */
    phi_f = sc[1]*sc[3];
    k_f = k_f_save[21];
    q_f[21] = phi_f * k_f;
    phi_r = sc[5]*sc[5];
    Kc = Kc_save[21];
    k_r = k_f / Kc;
    q_r[21] = phi_r * k_r;

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    phi_f = sc[7]*sc[5];
    k_f = k_f_save[22];
    q_f[22] = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = Kc_save[22];
    k_r = k_f / Kc;
    q_r[22] = phi_r * k_r;

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    phi_f = sc[7]*sc[5];
    k_f = k_f_save[23];
    q_f[23] = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = Kc_save[23];
    k_r = k_f / Kc;
    q_r[23] = phi_r * k_r;

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[7]*sc[7];
    k_f = k_f_save[24];
    q_f[24] = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = Kc_save[24];
    k_r = k_f / Kc;
    q_r[24] = phi_r * k_r;

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[7]*sc[7];
    k_f = k_f_save[25];
    q_f[25] = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = Kc_save[25];
    k_r = k_f / Kc;
    q_r[25] = phi_r * k_r;

    /*reaction 27: H2O2 + H <=> HO2 + H2 */
    phi_f = sc[6]*sc[0];
    k_f = k_f_save[26];
    q_f[26] = phi_f * k_f;
    phi_r = sc[7]*sc[1];
    Kc = Kc_save[26];
    k_r = k_f / Kc;
    q_r[26] = phi_r * k_r;

    /*reaction 28: H2O2 + H <=> H2O + OH */
    phi_f = sc[6]*sc[0];
    k_f = k_f_save[27];
    q_f[27] = phi_f * k_f;
    phi_r = sc[4]*sc[5];
    Kc = Kc_save[27];
    k_r = k_f / Kc;
    q_r[27] = phi_r * k_r;

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    phi_f = sc[6]*sc[2];
    k_f = k_f_save[28];
    q_f[28] = phi_f * k_f;
    phi_r = sc[7]*sc[5];
    Kc = Kc_save[28];
    k_r = k_f / Kc;
    q_r[28] = phi_r * k_r;

    /*reaction 30: H2O2 + OH <=> HO2 + H2O */
    phi_f = sc[6]*sc[5];
    k_f = k_f_save[29];
    q_f[29] = phi_f * k_f;
    phi_r = sc[7]*sc[4];
    Kc = Kc_save[29];
    k_r = k_f / Kc;
    q_r[29] = phi_r * k_r;

    /*reaction 31: H2O2 + OH <=> HO2 + H2O */
    phi_f = sc[6]*sc[5];
    k_f = k_f_save[30];
    q_f[30] = phi_f * k_f;
    phi_r = sc[7]*sc[4];
    Kc = Kc_save[30];
    k_r = k_f / Kc;
    q_r[30] = phi_r * k_r;

    return;
}


/*compute the equilibrium constants for each reaction */
void equilibriumConstants(double *kc, double * g_RT, double T)
{
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    /*reaction 1: H + H + M <=> H2 + M */
    kc[0] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0]) - (g_RT[1]));

    /*reaction 2: H + H + H2 <=> H2 + H2 */
    kc[1] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[1]) - (g_RT[1] + g_RT[1]));

    /*reaction 3: H + H + N2 <=> H2 + N2 */
    kc[2] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[9]) - (g_RT[1] + g_RT[9]));

    /*reaction 4: H + H + H <=> H2 + H */
    kc[3] = 1.0 / (refC) * exp((g_RT[0] + g_RT[0] + g_RT[0]) - (g_RT[1] + g_RT[0]));

    /*reaction 5: O + O + M <=> O2 + M */
    kc[4] = 1.0 / (refC) * exp((g_RT[2] + g_RT[2]) - (g_RT[3]));

    /*reaction 6: O + H + M <=> OH + M */
    kc[5] = 1.0 / (refC) * exp((g_RT[2] + g_RT[0]) - (g_RT[5]));

    /*reaction 7: H + OH + M <=> H2O + M */
    kc[6] = 1.0 / (refC) * exp((g_RT[0] + g_RT[5]) - (g_RT[4]));

    /*reaction 8: H + O2 (+M) <=> HO2 (+M) */
    kc[7] = 1.0 / (refC) * exp((g_RT[0] + g_RT[3]) - (g_RT[7]));

    /*reaction 9: H + O2 (+AR) <=> HO2 (+AR) */
    kc[8] = 1.0 / (refC) * exp((g_RT[0] + g_RT[3]) - (g_RT[7]));

    /*reaction 10: H + O2 (+O2) <=> HO2 (+O2) */
    kc[9] = 1.0 / (refC) * exp((g_RT[0] + g_RT[3]) - (g_RT[7]));

    /*reaction 11: H + O2 (+H2O) <=> HO2 (+H2O) */
    kc[10] = 1.0 / (refC) * exp((g_RT[0] + g_RT[3]) - (g_RT[7]));

    /*reaction 12: OH + OH (+M) <=> H2O2 (+M) */
    kc[11] = 1.0 / (refC) * exp((g_RT[5] + g_RT[5]) - (g_RT[6]));

    /*reaction 13: OH + OH (+H2O) <=> H2O2 (+H2O) */
    kc[12] = 1.0 / (refC) * exp((g_RT[5] + g_RT[5]) - (g_RT[6]));

    /*reaction 14: O + H2 <=> OH + H */
    kc[13] = exp((g_RT[2] + g_RT[1]) - (g_RT[5] + g_RT[0]));

    /*reaction 15: H + O2 <=> OH + O */
    kc[14] = exp((g_RT[0] + g_RT[3]) - (g_RT[5] + g_RT[2]));

    /*reaction 16: H2 + OH <=> H2O + H */
    kc[15] = exp((g_RT[1] + g_RT[5]) - (g_RT[4] + g_RT[0]));

    /*reaction 17: OH + OH <=> H2O + O */
    kc[16] = exp((g_RT[5] + g_RT[5]) - (g_RT[4] + g_RT[2]));

    /*reaction 18: HO2 + O <=> OH + O2 */
    kc[17] = exp((g_RT[7] + g_RT[2]) - (g_RT[5] + g_RT[3]));

    /*reaction 19: H + HO2 <=> OH + OH */
    kc[18] = exp((g_RT[0] + g_RT[7]) - (g_RT[5] + g_RT[5]));

    /*reaction 20: H + HO2 <=> H2O + O */
    kc[19] = exp((g_RT[0] + g_RT[7]) - (g_RT[4] + g_RT[2]));

    /*reaction 21: H2 + O2 <=> H + HO2 */
    kc[20] = exp((g_RT[1] + g_RT[3]) - (g_RT[0] + g_RT[7]));

    /*reaction 22: H2 + O2 <=> OH + OH */
    kc[21] = exp((g_RT[1] + g_RT[3]) - (g_RT[5] + g_RT[5]));

    /*reaction 23: HO2 + OH <=> H2O + O2 */
    kc[22] = exp((g_RT[7] + g_RT[5]) - (g_RT[4] + g_RT[3]));

    /*reaction 24: HO2 + OH <=> H2O + O2 */
    kc[23] = exp((g_RT[7] + g_RT[5]) - (g_RT[4] + g_RT[3]));

    /*reaction 25: HO2 + HO2 <=> H2O2 + O2 */
    kc[24] = exp((g_RT[7] + g_RT[7]) - (g_RT[6] + g_RT[3]));

    /*reaction 26: HO2 + HO2 <=> H2O2 + O2 */
    kc[25] = exp((g_RT[7] + g_RT[7]) - (g_RT[6] + g_RT[3]));

    /*reaction 27: H2O2 + H <=> HO2 + H2 */
    kc[26] = exp((g_RT[6] + g_RT[0]) - (g_RT[7] + g_RT[1]));

    /*reaction 28: H2O2 + H <=> H2O + OH */
    kc[27] = exp((g_RT[6] + g_RT[0]) - (g_RT[4] + g_RT[5]));

    /*reaction 29: H2O2 + O <=> HO2 + OH */
    kc[28] = exp((g_RT[6] + g_RT[2]) - (g_RT[7] + g_RT[5]));

    /*reaction 30: H2O2 + OH <=> HO2 + H2O */
    kc[29] = exp((g_RT[6] + g_RT[5]) - (g_RT[7] + g_RT[4]));

    /*reaction 31: H2O2 + OH <=> HO2 + H2O */
    kc[30] = exp((g_RT[6] + g_RT[5]) - (g_RT[7] + g_RT[4]));

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void gibbs(double * species, double * tc)
{

    /*temperature */
    double T = tc[1], invT = 1.0 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.547366000000000e+04 * invT
            +2.946682850000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -9.179241300000000e+02 * invT
            +1.661300720000000e+00
            -2.344302900000000e+00 * tc[0]
            -3.990212400000000e-03 * tc[1]
            +3.246319500000000e-06 * tc[2]
            -1.679747250000000e-09 * tc[3]
            +3.688014450000000e-13 * tc[4];
        /*species 2: O */
        species[2] =
            +2.912225920000000e+04 * invT
            +1.116333640000000e+00
            -3.168267100000000e+00 * tc[0]
            +1.639659420000000e-03 * tc[1]
            -1.107177326666667e-06 * tc[2]
            +5.106721866666666e-10 * tc[3]
            -1.056329855000000e-13 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.063943560000000e+03 * invT
            +1.247806300000001e-01
            -3.782456360000000e+00 * tc[0]
            +1.498367075000000e-03 * tc[1]
            -1.641217000000000e-06 * tc[2]
            +8.067745900000000e-10 * tc[3]
            -1.621864180000000e-13 * tc[4];
        /*species 4: H2O */
        species[4] =
            -3.029372600000000e+04 * invT
            +5.047644210000000e+00
            -4.198635200000000e+00 * tc[0]
            +1.018200850000000e-03 * tc[1]
            -1.086723600000000e-06 * tc[2]
            +4.573272416666666e-10 * tc[3]
            -8.859840000000001e-14 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.368898360000000e+03 * invT
            +4.095982717000000e+00
            -3.991984240000000e+00 * tc[0]
            +1.200533275000000e-03 * tc[1]
            -7.694400550000000e-07 * tc[2]
            +3.232635883333333e-10 * tc[3]
            -6.815975100000000e-14 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            -1.770674370000000e+04 * invT
            +1.041418300000000e+00
            -4.315151490000000e+00 * tc[0]
            +4.236953110000000e-04 * tc[1]
            -2.940072050000000e-06 * tc[2]
            +1.889691200000000e-09 * tc[3]
            -4.544750790000000e-13 * tc[4];
        /*species 7: HO2 */
        species[7] =
            +2.948087600000000e+02 * invT
            +5.850870000000001e-01
            -4.301788000000000e+00 * tc[0]
            +2.374510050000000e-03 * tc[1]
            -3.526325500000000e-06 * tc[2]
            +2.022996750000000e-09 * tc[3]
            -4.646033500000000e-13 * tc[4];
        /*species 8: AR */
        species[8] =
            -7.453750000000000e+02 * invT
            -1.879674900000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 9: N2 */
        species[9] =
            -1.046976280000000e+03 * invT
            +5.635349000000001e-01
            -3.531005280000000e+00 * tc[0]
            +6.183049400000000e-05 * tc[1]
            +8.383323883333334e-08 * tc[2]
            -2.029421766666667e-10 * tc[3]
            +7.044061750000001e-14 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.547366000000000e+04 * invT
            +2.946682850000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -8.130558200000000e+02 * invT
            +3.957146900000000e+00
            -2.932830500000000e+00 * tc[0]
            -4.132990100000000e-04 * tc[1]
            +2.440009500000000e-08 * tc[2]
            -1.284154250000000e-12 * tc[3]
            +3.443980750000000e-17 * tc[4];
        /*species 2: O */
        species[2] =
            +2.922601200000000e+04 * invT
            -2.378657600000000e+00
            -2.543636970000000e+00 * tc[0]
            +1.365812430000000e-05 * tc[1]
            +6.983825333333333e-10 * tc[2]
            -4.129015375000000e-13 * tc[3]
            +2.397768470000000e-17 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.215977250000000e+03 * invT
            +2.455989900000000e-01
            -3.660960830000000e+00 * tc[0]
            -3.281827615000000e-04 * tc[1]
            +2.352491416666667e-08 * tc[2]
            -1.714980483333333e-12 * tc[3]
            +6.495662400000000e-17 * tc[4];
        /*species 4: H2O */
        species[4] =
            -2.988589400000000e+04 * invT
            -4.205511100000001e+00
            -2.677038900000000e+00 * tc[0]
            -1.486590800000000e-03 * tc[1]
            +1.289614816666667e-07 * tc[2]
            -7.869459500000001e-12 * tc[3]
            +2.134499550000000e-16 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.697808080000000e+03 * invT
            -3.006416189999999e+00
            -2.838530330000000e+00 * tc[0]
            -5.537064450000000e-04 * tc[1]
            +4.900003483333333e-08 * tc[2]
            -3.505822741666666e-12 * tc[3]
            +1.211449450000000e-16 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            -1.800717750000000e+04 * invT
            +3.914802356000000e+00
            -4.579773050000000e+00 * tc[0]
            -2.026630015000000e-03 * tc[1]
            +2.164078833333333e-07 * tc[2]
            -1.651761666666667e-11 * tc[3]
            +5.698439600000000e-16 * tc[4];
        /*species 7: HO2 */
        species[7] =
            +6.181885100000000e+01 * invT
            +1.214468500000000e+00
            -4.172265900000000e+00 * tc[0]
            -9.406049000000000e-04 * tc[1]
            +5.771549500000000e-08 * tc[2]
            -1.622376333333333e-12 * tc[3]
            -8.804576500000000e-18 * tc[4];
        /*species 8: AR */
        species[8] =
            -7.453750000000000e+02 * invT
            -1.879674900000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 9: N2 */
        species[9] =
            -9.239486879999999e+02 * invT
            -2.919311250000000e+00
            -2.952576370000000e+00 * tc[0]
            -6.984502000000001e-04 * tc[1]
            +8.210526716666667e-08 * tc[2]
            -6.550084958333333e-12 * tc[3]
            +2.303776020000000e-16 * tc[4];
    }
    return;
}


/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void helmholtz(double * species, double * tc)
{

    /*temperature */
    double T = tc[1], invT = 1.0 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.54736600e+04 * invT
            +1.94668285e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -9.17924130e+02 * invT
            +6.61300720e-01
            -2.34430290e+00 * tc[0]
            -3.99021240e-03 * tc[1]
            +3.24631950e-06 * tc[2]
            -1.67974725e-09 * tc[3]
            +3.68801445e-13 * tc[4];
        /*species 2: O */
        species[2] =
            +2.91222592e+04 * invT
            +1.16333640e-01
            -3.16826710e+00 * tc[0]
            +1.63965942e-03 * tc[1]
            -1.10717733e-06 * tc[2]
            +5.10672187e-10 * tc[3]
            -1.05632985e-13 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.06394356e+03 * invT
            -8.75219370e-01
            -3.78245636e+00 * tc[0]
            +1.49836707e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774590e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 4: H2O */
        species[4] =
            -3.02937260e+04 * invT
            +4.04764421e+00
            -4.19863520e+00 * tc[0]
            +1.01820085e-03 * tc[1]
            -1.08672360e-06 * tc[2]
            +4.57327242e-10 * tc[3]
            -8.85984000e-14 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.36889836e+03 * invT
            +3.09598272e+00
            -3.99198424e+00 * tc[0]
            +1.20053327e-03 * tc[1]
            -7.69440055e-07 * tc[2]
            +3.23263588e-10 * tc[3]
            -6.81597510e-14 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            -1.77067437e+04 * invT
            +4.14183000e-02
            -4.31515149e+00 * tc[0]
            +4.23695311e-04 * tc[1]
            -2.94007205e-06 * tc[2]
            +1.88969120e-09 * tc[3]
            -4.54475079e-13 * tc[4];
        /*species 7: HO2 */
        species[7] =
            +2.94808760e+02 * invT
            -4.14913000e-01
            -4.30178800e+00 * tc[0]
            +2.37451005e-03 * tc[1]
            -3.52632550e-06 * tc[2]
            +2.02299675e-09 * tc[3]
            -4.64603350e-13 * tc[4];
        /*species 8: AR */
        species[8] =
            -7.45375000e+02 * invT
            -2.87967490e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 9: N2 */
        species[9] =
            -1.04697628e+03 * invT
            -4.36465100e-01
            -3.53100528e+00 * tc[0]
            +6.18304940e-05 * tc[1]
            +8.38332388e-08 * tc[2]
            -2.02942177e-10 * tc[3]
            +7.04406175e-14 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.54736600e+04 * invT
            +1.94668285e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -8.13055820e+02 * invT
            +2.95714690e+00
            -2.93283050e+00 * tc[0]
            -4.13299010e-04 * tc[1]
            +2.44000950e-08 * tc[2]
            -1.28415425e-12 * tc[3]
            +3.44398075e-17 * tc[4];
        /*species 2: O */
        species[2] =
            +2.92260120e+04 * invT
            -3.37865760e+00
            -2.54363697e+00 * tc[0]
            +1.36581243e-05 * tc[1]
            +6.98382533e-10 * tc[2]
            -4.12901538e-13 * tc[3]
            +2.39776847e-17 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.21597725e+03 * invT
            -7.54401010e-01
            -3.66096083e+00 * tc[0]
            -3.28182761e-04 * tc[1]
            +2.35249142e-08 * tc[2]
            -1.71498048e-12 * tc[3]
            +6.49566240e-17 * tc[4];
        /*species 4: H2O */
        species[4] =
            -2.98858940e+04 * invT
            -5.20551110e+00
            -2.67703890e+00 * tc[0]
            -1.48659080e-03 * tc[1]
            +1.28961482e-07 * tc[2]
            -7.86945950e-12 * tc[3]
            +2.13449955e-16 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.69780808e+03 * invT
            -4.00641619e+00
            -2.83853033e+00 * tc[0]
            -5.53706445e-04 * tc[1]
            +4.90000348e-08 * tc[2]
            -3.50582274e-12 * tc[3]
            +1.21144945e-16 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            -1.80071775e+04 * invT
            +2.91480236e+00
            -4.57977305e+00 * tc[0]
            -2.02663002e-03 * tc[1]
            +2.16407883e-07 * tc[2]
            -1.65176167e-11 * tc[3]
            +5.69843960e-16 * tc[4];
        /*species 7: HO2 */
        species[7] =
            +6.18188510e+01 * invT
            +2.14468500e-01
            -4.17226590e+00 * tc[0]
            -9.40604900e-04 * tc[1]
            +5.77154950e-08 * tc[2]
            -1.62237633e-12 * tc[3]
            -8.80457650e-18 * tc[4];
        /*species 8: AR */
        species[8] =
            -7.45375000e+02 * invT
            -2.87967490e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 9: N2 */
        species[9] =
            -9.23948688e+02 * invT
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
    double T = tc[1], invT = 1.0 / T;

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
            +1.34430290e+00
            +7.98042480e-03 * tc[1]
            -1.94779170e-05 * tc[2]
            +2.01569670e-08 * tc[3]
            -7.37602890e-12 * tc[4];
        /*species 2: O */
        species[2] =
            +2.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 3: O2 */
        species[3] =
            +2.78245636e+00
            -2.99673415e-03 * tc[1]
            +9.84730200e-06 * tc[2]
            -9.68129508e-09 * tc[3]
            +3.24372836e-12 * tc[4];
        /*species 4: H2O */
        species[4] =
            +3.19863520e+00
            -2.03640170e-03 * tc[1]
            +6.52034160e-06 * tc[2]
            -5.48792690e-09 * tc[3]
            +1.77196800e-12 * tc[4];
        /*species 5: OH */
        species[5] =
            +2.99198424e+00
            -2.40106655e-03 * tc[1]
            +4.61664033e-06 * tc[2]
            -3.87916306e-09 * tc[3]
            +1.36319502e-12 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            +3.31515149e+00
            -8.47390622e-04 * tc[1]
            +1.76404323e-05 * tc[2]
            -2.26762944e-08 * tc[3]
            +9.08950158e-12 * tc[4];
        /*species 7: HO2 */
        species[7] =
            +3.30178800e+00
            -4.74902010e-03 * tc[1]
            +2.11579530e-05 * tc[2]
            -2.42759610e-08 * tc[3]
            +9.29206700e-12 * tc[4];
        /*species 8: AR */
        species[8] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 9: N2 */
        species[9] =
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
        /*species 1: H2 */
        species[1] =
            +1.93283050e+00
            +8.26598020e-04 * tc[1]
            -1.46400570e-07 * tc[2]
            +1.54098510e-11 * tc[3]
            -6.88796150e-16 * tc[4];
        /*species 2: O */
        species[2] =
            +1.54363697e+00
            -2.73162486e-05 * tc[1]
            -4.19029520e-09 * tc[2]
            +4.95481845e-12 * tc[3]
            -4.79553694e-16 * tc[4];
        /*species 3: O2 */
        species[3] =
            +2.66096083e+00
            +6.56365523e-04 * tc[1]
            -1.41149485e-07 * tc[2]
            +2.05797658e-11 * tc[3]
            -1.29913248e-15 * tc[4];
        /*species 4: H2O */
        species[4] =
            +1.67703890e+00
            +2.97318160e-03 * tc[1]
            -7.73768890e-07 * tc[2]
            +9.44335140e-11 * tc[3]
            -4.26899910e-15 * tc[4];
        /*species 5: OH */
        species[5] =
            +1.83853033e+00
            +1.10741289e-03 * tc[1]
            -2.94000209e-07 * tc[2]
            +4.20698729e-11 * tc[3]
            -2.42289890e-15 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            +3.57977305e+00
            +4.05326003e-03 * tc[1]
            -1.29844730e-06 * tc[2]
            +1.98211400e-10 * tc[3]
            -1.13968792e-14 * tc[4];
        /*species 7: HO2 */
        species[7] =
            +3.17226590e+00
            +1.88120980e-03 * tc[1]
            -3.46292970e-07 * tc[2]
            +1.94685160e-11 * tc[3]
            +1.76091530e-16 * tc[4];
        /*species 8: AR */
        species[8] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 9: N2 */
        species[9] =
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
    double T = tc[1], invT = 1.0 / T;

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
            +2.34430290e+00
            +7.98042480e-03 * tc[1]
            -1.94779170e-05 * tc[2]
            +2.01569670e-08 * tc[3]
            -7.37602890e-12 * tc[4];
        /*species 2: O */
        species[2] =
            +3.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00
            -2.99673415e-03 * tc[1]
            +9.84730200e-06 * tc[2]
            -9.68129508e-09 * tc[3]
            +3.24372836e-12 * tc[4];
        /*species 4: H2O */
        species[4] =
            +4.19863520e+00
            -2.03640170e-03 * tc[1]
            +6.52034160e-06 * tc[2]
            -5.48792690e-09 * tc[3]
            +1.77196800e-12 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.99198424e+00
            -2.40106655e-03 * tc[1]
            +4.61664033e-06 * tc[2]
            -3.87916306e-09 * tc[3]
            +1.36319502e-12 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            +4.31515149e+00
            -8.47390622e-04 * tc[1]
            +1.76404323e-05 * tc[2]
            -2.26762944e-08 * tc[3]
            +9.08950158e-12 * tc[4];
        /*species 7: HO2 */
        species[7] =
            +4.30178800e+00
            -4.74902010e-03 * tc[1]
            +2.11579530e-05 * tc[2]
            -2.42759610e-08 * tc[3]
            +9.29206700e-12 * tc[4];
        /*species 8: AR */
        species[8] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 9: N2 */
        species[9] =
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
        /*species 1: H2 */
        species[1] =
            +2.93283050e+00
            +8.26598020e-04 * tc[1]
            -1.46400570e-07 * tc[2]
            +1.54098510e-11 * tc[3]
            -6.88796150e-16 * tc[4];
        /*species 2: O */
        species[2] =
            +2.54363697e+00
            -2.73162486e-05 * tc[1]
            -4.19029520e-09 * tc[2]
            +4.95481845e-12 * tc[3]
            -4.79553694e-16 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.66096083e+00
            +6.56365523e-04 * tc[1]
            -1.41149485e-07 * tc[2]
            +2.05797658e-11 * tc[3]
            -1.29913248e-15 * tc[4];
        /*species 4: H2O */
        species[4] =
            +2.67703890e+00
            +2.97318160e-03 * tc[1]
            -7.73768890e-07 * tc[2]
            +9.44335140e-11 * tc[3]
            -4.26899910e-15 * tc[4];
        /*species 5: OH */
        species[5] =
            +2.83853033e+00
            +1.10741289e-03 * tc[1]
            -2.94000209e-07 * tc[2]
            +4.20698729e-11 * tc[3]
            -2.42289890e-15 * tc[4];
        /*species 6: H2O2 */
        species[6] =
            +4.57977305e+00
            +4.05326003e-03 * tc[1]
            -1.29844730e-06 * tc[2]
            +1.98211400e-10 * tc[3]
            -1.13968792e-14 * tc[4];
        /*species 7: HO2 */
        species[7] =
            +4.17226590e+00
            +1.88120980e-03 * tc[1]
            -3.46292970e-07 * tc[2]
            +1.94685160e-11 * tc[3]
            +1.76091530e-16 * tc[4];
        /*species 8: AR */
        species[8] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 9: N2 */
        species[9] =
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
    double T = tc[1], invT = 1.0 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54736600e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +1.34430290e+00
            +3.99021240e-03 * tc[1]
            -6.49263900e-06 * tc[2]
            +5.03924175e-09 * tc[3]
            -1.47520578e-12 * tc[4]
            -9.17924130e+02 * invT;
        /*species 2: O */
        species[2] =
            +2.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +2.78245636e+00
            -1.49836707e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745672e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +3.19863520e+00
            -1.01820085e-03 * tc[1]
            +2.17344720e-06 * tc[2]
            -1.37198172e-09 * tc[3]
            +3.54393600e-13 * tc[4]
            -3.02937260e+04 * invT;
        /*species 5: OH */
        species[5] =
            +2.99198424e+00
            -1.20053327e-03 * tc[1]
            +1.53888011e-06 * tc[2]
            -9.69790765e-10 * tc[3]
            +2.72639004e-13 * tc[4]
            +3.36889836e+03 * invT;
        /*species 6: H2O2 */
        species[6] =
            +3.31515149e+00
            -4.23695311e-04 * tc[1]
            +5.88014410e-06 * tc[2]
            -5.66907360e-09 * tc[3]
            +1.81790032e-12 * tc[4]
            -1.77067437e+04 * invT;
        /*species 7: HO2 */
        species[7] =
            +3.30178800e+00
            -2.37451005e-03 * tc[1]
            +7.05265100e-06 * tc[2]
            -6.06899025e-09 * tc[3]
            +1.85841340e-12 * tc[4]
            +2.94808760e+02 * invT;
        /*species 8: AR */
        species[8] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 9: N2 */
        species[9] =
            +2.53100528e+00
            -6.18304940e-05 * tc[1]
            -1.67666478e-07 * tc[2]
            +6.08826530e-10 * tc[3]
            -2.81762470e-13 * tc[4]
            -1.04697628e+03 * invT;
    } else {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54736600e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +1.93283050e+00
            +4.13299010e-04 * tc[1]
            -4.88001900e-08 * tc[2]
            +3.85246275e-12 * tc[3]
            -1.37759230e-16 * tc[4]
            -8.13055820e+02 * invT;
        /*species 2: O */
        species[2] =
            +1.54363697e+00
            -1.36581243e-05 * tc[1]
            -1.39676507e-09 * tc[2]
            +1.23870461e-12 * tc[3]
            -9.59107388e-17 * tc[4]
            +2.92260120e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +2.66096083e+00
            +3.28182761e-04 * tc[1]
            -4.70498283e-08 * tc[2]
            +5.14494145e-12 * tc[3]
            -2.59826496e-16 * tc[4]
            -1.21597725e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +1.67703890e+00
            +1.48659080e-03 * tc[1]
            -2.57922963e-07 * tc[2]
            +2.36083785e-11 * tc[3]
            -8.53799820e-16 * tc[4]
            -2.98858940e+04 * invT;
        /*species 5: OH */
        species[5] =
            +1.83853033e+00
            +5.53706445e-04 * tc[1]
            -9.80000697e-08 * tc[2]
            +1.05174682e-11 * tc[3]
            -4.84579780e-16 * tc[4]
            +3.69780808e+03 * invT;
        /*species 6: H2O2 */
        species[6] =
            +3.57977305e+00
            +2.02663002e-03 * tc[1]
            -4.32815767e-07 * tc[2]
            +4.95528500e-11 * tc[3]
            -2.27937584e-15 * tc[4]
            -1.80071775e+04 * invT;
        /*species 7: HO2 */
        species[7] =
            +3.17226590e+00
            +9.40604900e-04 * tc[1]
            -1.15430990e-07 * tc[2]
            +4.86712900e-12 * tc[3]
            +3.52183060e-17 * tc[4]
            +6.18188510e+01 * invT;
        /*species 8: AR */
        species[8] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 9: N2 */
        species[9] =
            +1.95257637e+00
            +6.98450200e-04 * tc[1]
            -1.64210534e-07 * tc[2]
            +1.96502549e-11 * tc[3]
            -9.21510408e-16 * tc[4]
            -9.23948688e+02 * invT;
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEnthalpy(double * species, double * tc)
{

    /*temperature */
    double T = tc[1], invT = 1.0 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54736600e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +2.34430290e+00
            +3.99021240e-03 * tc[1]
            -6.49263900e-06 * tc[2]
            +5.03924175e-09 * tc[3]
            -1.47520578e-12 * tc[4]
            -9.17924130e+02 * invT;
        /*species 2: O */
        species[2] =
            +3.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00
            -1.49836707e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745672e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +4.19863520e+00
            -1.01820085e-03 * tc[1]
            +2.17344720e-06 * tc[2]
            -1.37198172e-09 * tc[3]
            +3.54393600e-13 * tc[4]
            -3.02937260e+04 * invT;
        /*species 5: OH */
        species[5] =
            +3.99198424e+00
            -1.20053327e-03 * tc[1]
            +1.53888011e-06 * tc[2]
            -9.69790765e-10 * tc[3]
            +2.72639004e-13 * tc[4]
            +3.36889836e+03 * invT;
        /*species 6: H2O2 */
        species[6] =
            +4.31515149e+00
            -4.23695311e-04 * tc[1]
            +5.88014410e-06 * tc[2]
            -5.66907360e-09 * tc[3]
            +1.81790032e-12 * tc[4]
            -1.77067437e+04 * invT;
        /*species 7: HO2 */
        species[7] =
            +4.30178800e+00
            -2.37451005e-03 * tc[1]
            +7.05265100e-06 * tc[2]
            -6.06899025e-09 * tc[3]
            +1.85841340e-12 * tc[4]
            +2.94808760e+02 * invT;
        /*species 8: AR */
        species[8] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 9: N2 */
        species[9] =
            +3.53100528e+00
            -6.18304940e-05 * tc[1]
            -1.67666478e-07 * tc[2]
            +6.08826530e-10 * tc[3]
            -2.81762470e-13 * tc[4]
            -1.04697628e+03 * invT;
    } else {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54736600e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +2.93283050e+00
            +4.13299010e-04 * tc[1]
            -4.88001900e-08 * tc[2]
            +3.85246275e-12 * tc[3]
            -1.37759230e-16 * tc[4]
            -8.13055820e+02 * invT;
        /*species 2: O */
        species[2] =
            +2.54363697e+00
            -1.36581243e-05 * tc[1]
            -1.39676507e-09 * tc[2]
            +1.23870461e-12 * tc[3]
            -9.59107388e-17 * tc[4]
            +2.92260120e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +3.66096083e+00
            +3.28182761e-04 * tc[1]
            -4.70498283e-08 * tc[2]
            +5.14494145e-12 * tc[3]
            -2.59826496e-16 * tc[4]
            -1.21597725e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +2.67703890e+00
            +1.48659080e-03 * tc[1]
            -2.57922963e-07 * tc[2]
            +2.36083785e-11 * tc[3]
            -8.53799820e-16 * tc[4]
            -2.98858940e+04 * invT;
        /*species 5: OH */
        species[5] =
            +2.83853033e+00
            +5.53706445e-04 * tc[1]
            -9.80000697e-08 * tc[2]
            +1.05174682e-11 * tc[3]
            -4.84579780e-16 * tc[4]
            +3.69780808e+03 * invT;
        /*species 6: H2O2 */
        species[6] =
            +4.57977305e+00
            +2.02663002e-03 * tc[1]
            -4.32815767e-07 * tc[2]
            +4.95528500e-11 * tc[3]
            -2.27937584e-15 * tc[4]
            -1.80071775e+04 * invT;
        /*species 7: HO2 */
        species[7] =
            +4.17226590e+00
            +9.40604900e-04 * tc[1]
            -1.15430990e-07 * tc[2]
            +4.86712900e-12 * tc[3]
            +3.52183060e-17 * tc[4]
            +6.18188510e+01 * invT;
        /*species 8: AR */
        species[8] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 9: N2 */
        species[9] =
            +2.95257637e+00
            +6.98450200e-04 * tc[1]
            -1.64210534e-07 * tc[2]
            +1.96502549e-11 * tc[3]
            -9.21510408e-16 * tc[4]
            -9.23948688e+02 * invT;
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEntropy(double * species, double * tc)
{

    /*temperature */
    double T = tc[1], invT = 1.0 / T;

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
        /*species 1: H2 */
        species[1] =
            +2.34430290e+00 * tc[0]
            +7.98042480e-03 * tc[1]
            -9.73895850e-06 * tc[2]
            +6.71898900e-09 * tc[3]
            -1.84400722e-12 * tc[4]
            +6.83002180e-01 ;
        /*species 2: O */
        species[2] =
            +3.16826710e+00 * tc[0]
            -3.27931884e-03 * tc[1]
            +3.32153198e-06 * tc[2]
            -2.04268875e-09 * tc[3]
            +5.28164927e-13 * tc[4]
            +2.05193346e+00 ;
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00 * tc[0]
            -2.99673415e-03 * tc[1]
            +4.92365100e-06 * tc[2]
            -3.22709836e-09 * tc[3]
            +8.10932090e-13 * tc[4]
            +3.65767573e+00 ;
        /*species 4: H2O */
        species[4] =
            +4.19863520e+00 * tc[0]
            -2.03640170e-03 * tc[1]
            +3.26017080e-06 * tc[2]
            -1.82930897e-09 * tc[3]
            +4.42992000e-13 * tc[4]
            -8.49009010e-01 ;
        /*species 5: OH */
        species[5] =
            +3.99198424e+00 * tc[0]
            -2.40106655e-03 * tc[1]
            +2.30832017e-06 * tc[2]
            -1.29305435e-09 * tc[3]
            +3.40798755e-13 * tc[4]
            -1.03998477e-01 ;
        /*species 6: H2O2 */
        species[6] =
            +4.31515149e+00 * tc[0]
            -8.47390622e-04 * tc[1]
            +8.82021615e-06 * tc[2]
            -7.55876480e-09 * tc[3]
            +2.27237539e-12 * tc[4]
            +3.27373319e+00 ;
        /*species 7: HO2 */
        species[7] =
            +4.30178800e+00 * tc[0]
            -4.74902010e-03 * tc[1]
            +1.05789765e-05 * tc[2]
            -8.09198700e-09 * tc[3]
            +2.32301675e-12 * tc[4]
            +3.71670100e+00 ;
        /*species 8: AR */
        species[8] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.37967490e+00 ;
        /*species 9: N2 */
        species[9] =
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
        /*species 1: H2 */
        species[1] =
            +2.93283050e+00 * tc[0]
            +8.26598020e-04 * tc[1]
            -7.32002850e-08 * tc[2]
            +5.13661700e-12 * tc[3]
            -1.72199037e-16 * tc[4]
            -1.02431640e+00 ;
        /*species 2: O */
        species[2] =
            +2.54363697e+00 * tc[0]
            -2.73162486e-05 * tc[1]
            -2.09514760e-09 * tc[2]
            +1.65160615e-12 * tc[3]
            -1.19888423e-16 * tc[4]
            +4.92229457e+00 ;
        /*species 3: O2 */
        species[3] =
            +3.66096083e+00 * tc[0]
            +6.56365523e-04 * tc[1]
            -7.05747425e-08 * tc[2]
            +6.85992193e-12 * tc[3]
            -3.24783120e-16 * tc[4]
            +3.41536184e+00 ;
        /*species 4: H2O */
        species[4] =
            +2.67703890e+00 * tc[0]
            +2.97318160e-03 * tc[1]
            -3.86884445e-07 * tc[2]
            +3.14778380e-11 * tc[3]
            -1.06724977e-15 * tc[4]
            +6.88255000e+00 ;
        /*species 5: OH */
        species[5] =
            +2.83853033e+00 * tc[0]
            +1.10741289e-03 * tc[1]
            -1.47000104e-07 * tc[2]
            +1.40232910e-11 * tc[3]
            -6.05724725e-16 * tc[4]
            +5.84494652e+00 ;
        /*species 6: H2O2 */
        species[6] =
            +4.57977305e+00 * tc[0]
            +4.05326003e-03 * tc[1]
            -6.49223650e-07 * tc[2]
            +6.60704667e-11 * tc[3]
            -2.84921980e-15 * tc[4]
            +6.64970694e-01 ;
        /*species 7: HO2 */
        species[7] =
            +4.17226590e+00 * tc[0]
            +1.88120980e-03 * tc[1]
            -1.73146485e-07 * tc[2]
            +6.48950533e-12 * tc[3]
            +4.40228825e-17 * tc[4]
            +2.95779740e+00 ;
        /*species 8: AR */
        species[8] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.37967490e+00 ;
        /*species 9: N2 */
        species[9] =
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
    wt[1] = 2.015940; /*H2 */
    wt[2] = 15.999400; /*O */
    wt[3] = 31.998800; /*O2 */
    wt[4] = 18.015340; /*H2O */
    wt[5] = 17.007370; /*OH */
    wt[6] = 34.014740; /*H2O2 */
    wt[7] = 33.006770; /*HO2 */
    wt[8] = 39.948000; /*AR */
    wt[9] = 28.013400; /*N2 */

    return;
}


/*get temperature given internal energy in mass units and mass fracs */
int feeytt_(double * e, double * y, int * iwrk, double * rwrk, double * t)
{
    const int maxiter = 50;
    const double tol  = 0.001;
    double ein  = *e;
    double tmin = 200; // max lower bound for thermo def
    double tmax = 6000; // min upper bound for thermo def
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
    y[1] = phi[1]*2.015940;   XW += y[1]; /*H2 */
    y[2] = phi[2]*15.999400;   XW += y[2]; /*O */
    y[3] = phi[3]*31.998800;   XW += y[3]; /*O2 */
    y[4] = phi[4]*18.015340;   XW += y[4]; /*H2O */
    y[5] = phi[5]*17.007370;   XW += y[5]; /*OH */
    y[6] = phi[6]*34.014740;   XW += y[6]; /*H2O2 */
    y[7] = phi[7]*33.006770;   XW += y[7]; /*HO2 */
    y[8] = phi[8]*39.948000;   XW += y[8]; /*AR */
    y[9] = phi[9]*28.013400;   XW += y[9]; /*N2 */
    for (id = 0; id < 10; ++id) {
        y[id] = y[id]/XW;
    }

    return;
}


/*convert y[species] (mass fracs) to phi[species] (specific mole num) */
void feytphi_(double * y, int * iwrk, double * rwrk, double * phi)
{
    phi[0] = y[0]/ 1.00797000e-03; /*H (wt in kg) */
    phi[1] = y[1]/ 2.01594000e-03; /*H2 (wt in kg) */
    phi[2] = y[2]/ 1.59994000e-02; /*O (wt in kg) */
    phi[3] = y[3]/ 3.19988000e-02; /*O2 (wt in kg) */
    phi[4] = y[4]/ 1.80153400e-02; /*H2O (wt in kg) */
    phi[5] = y[5]/ 1.70073700e-02; /*OH (wt in kg) */
    phi[6] = y[6]/ 3.40147400e-02; /*H2O2 (wt in kg) */
    phi[7] = y[7]/ 3.30067700e-02; /*HO2 (wt in kg) */
    phi[8] = y[8]/ 3.99480000e-02; /*AR (wt in kg) */
    phi[9] = y[9]/ 2.80134000e-02; /*N2 (wt in kg) */

    return;
}


/*reverse of ytcr, useful for rate computations */
void fectyr_(double * c, double * rho, int * iwrk, double * rwrk, double * y)
{
    y[0] = c[0] * 1.007970 / (*rho); 
    y[1] = c[1] * 2.015940 / (*rho); 
    y[2] = c[2] * 15.999400 / (*rho); 
    y[3] = c[3] * 31.998800 / (*rho); 
    y[4] = c[4] * 18.015340 / (*rho); 
    y[5] = c[5] * 17.007370 / (*rho); 
    y[6] = c[6] * 34.014740 / (*rho); 
    y[7] = c[7] * 33.006770 / (*rho); 
    y[8] = c[8] * 39.948000 / (*rho); 
    y[9] = c[9] * 28.013400 / (*rho); 

    return;
}


/*ddebdf compatible right hand side of CV burner */
/*rwrk[0] and rwrk[1] should contain rho and ene respectively */
/*working variable phi contains specific mole numbers */
void fecvrhs_(double * time, double * phi, double * phidot, double * rwrk, int * iwrk)
{
    double rho,ene; /*CV Parameters */
    double y[10], wdot[10]; /*temporary storage */
    int i; /*Loop counter */
    double temperature,pressure; /*temporary var */
    rho = rwrk[0];
    ene = rwrk[1];
    fephity_(phi, iwrk, rwrk, y);
    feeytt_(&ene, y, iwrk, rwrk, &temperature);
    CKPY(&rho, &temperature,  y, iwrk, rwrk, &pressure);
    CKWYP(&pressure, &temperature,  y, iwrk, rwrk, wdot);
    for (i=0; i<10; ++i) phidot[i] = wdot[i] / (rho/1000.0); 

    return;
}


/*returns the dimensionality of the cv burner (number of species) */
int fecvdim_()
{
    return 10;
}


/*ddebdf compatible right hand side of ZND solver */
/*rwrk[0] : scaling factor for pressure */
/*rwrk[1] : preshock density (g/cc)  */
/*rwrk[2] : detonation velocity (cm/s)  */
/*solution vector: [P; rho; y0 ... ylast]  */
void fezndrhs_(double * time, double * z, double * zdot, double * rwrk, int * iwrk)
{
    double psc,rho1,udet; /*ZND Parameters */
    double wt[10], hms[10], wdot[10]; /*temporary storage */
    int i; /*Loop counter */
    /*temporary variables */
    double ru, T, uvel, wtm, p, rho, gam, son, xm, sum, drdy, eta, cp, cv ;
    double *y; /*mass frac pointer */

    ru = 8.31451e+07;

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
    for (i=0; i<10; ++i) {
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
    return 13;
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
    if (strcmp(s1, "H2")==0) return 1; 
    if (strcmp(s1, "O")==0) return 2; 
    if (strcmp(s1, "O2")==0) return 3; 
    if (strcmp(s1, "H2O")==0) return 4; 
    if (strcmp(s1, "OH")==0) return 5; 
    if (strcmp(s1, "H2O2")==0) return 6; 
    if (strcmp(s1, "HO2")==0) return 7; 
    if (strcmp(s1, "AR")==0) return 8; 
    if (strcmp(s1, "N2")==0) return 9; 
    /*species name not found */
    return -1;
}


/*returns the species name */
char* fesymname_(int sn)
{
    if (sn==0) return "H"; 
    if (sn==1) return "H2"; 
    if (sn==2) return "O"; 
    if (sn==3) return "O2"; 
    if (sn==4) return "H2O"; 
    if (sn==5) return "OH"; 
    if (sn==6) return "H2O2"; 
    if (sn==7) return "HO2"; 
    if (sn==8) return "AR"; 
    if (sn==9) return "N2"; 
    /*species name not found */
    return "NOTFOUND";
}

/* End of file  */



#if 0



\\
\\
\\  This is the mechanism file
\\
\\
ELEM  O N AR H
END
SPECIES ! structure, source of thermo-data, source of transport 
H  ! burcat, chemkin
H2 ! burcat, chemkin
O   O2    H2O    OH    H2O2    HO2
AR ! burcat, chemkin 
N2 ! burcat, chemkin 

END
REACTIONS

!*********************************************************************
!
! A.KONNOV's detailed reaction mechanism   h/o 
!
!*********************************************************************
H+H+M=H2+M                     7.000E+17     -1.0         0.0
 H2/0.0/ N2/0.0/ H/0.0/ H2O/14.3/ ! CO/3.0/ CO2/3.0/ 
H+H+H2=H2+H2                   1.000E+17     -0.6         0.0 
H+H+N2=H2+N2                   5.400E+18     -1.3         0.0 
H+H+H=H2+H                     3.200E+15      0.0         0.0 
O+O+M=O2+M                     1.000E+17     -1.0         0.0 
 O/28.8/ O2/8.0/ N2/2.0/ H2O/5.0/ ! NO/2.0/ N/2.0/ 
O+H+M=OH+M                     6.750E+18     -1.0         0.0 
 H2O/5.0/ 
H+OH+M=H2O+M                   2.200E+22     -2.0         0.0 
    H2O/6.4/ AR/0.38/ ! CO2/1.9/ 
H+O2(+M)=HO2(+M)               4.660E+12      0.44        0.0 
    LOW /5.70E+19 -1.4 0.0/ 
    TROE /0.5 100000 10/
 AR/0.0/ H2O/0.0/ O2/0.0/ H2/1.5/ ! CO2/2.4/ CH4/3.5/ !HE/0.57/ 
H+O2(+AR)=HO2(+AR)             4.660E+12      0.44        0.0
    LOW /7.430E+18 -1.2  0.0/
    TROE /0.5 10 100000/ 
H+O2(+O2)=HO2(+O2)             4.660E+12      0.44        0.0
    LOW /5.690E+18 -1.094 0.0/ 
    TROE /0.5 100000 10/
H+O2(+H2O)=HO2(+H2O)           9.060E+12      0.2         0.0
    LOW /3.670E+19 -1.0 0.0/ 
    TROE /0.8 10 100000/ 
OH+OH(+M)=H2O2(+M)             7.200E+13     -0.37        0.0 
    LOW /2.2E+19 -0.76 0.0/ 
    TROE /0.5 100000 10/
    H2O/0.0/
OH+OH(+H2O)=H2O2(+H2O)         7.200E+13     -0.37        0.0 
    LOW /1.45E+18 0.0 0.0/ 
O+H2=OH+H                      5.060E+04      2.67     6290.0 
H+O2=OH+O                      9.750E+13      0.0     14850.0 
H2+OH=H2O+H                    1.000E+08      1.6      3300.0 
OH+OH=H2O+O                    3.570E+04      2.4     -2110.0 
HO2+O=OH+O2                    1.630E+13      0.0      -445.0 
H+HO2=OH+OH                    1.700E+14      0.0       875.0 
H+HO2=H2O+O                    3.000E+13      0.0      1720.0 
H2+O2=H+HO2                    7.400E+05     2.4328   53500.0 
H2+O2=OH+OH                    2.040E+12      0.44    69155.0 
HO2+OH=H2O+O2                  1.000E+13      0.0         0.0 
    DUPLICATE
HO2+OH=H2O+O2                  5.800E+13      0.0      3974.0 
    DUPLICATE
HO2+HO2=H2O2+O2                1.030E+14      0.0     11040.0 
    DUPLICATE
HO2+HO2=H2O2+O2                1.940E+11      0.0     -1409.0 
    DUPLICATE
H2O2+H=HO2+H2                  1.700E+12      0.0      3755.0 
H2O2+H=H2O+OH                  1.000E+13      0.0      3575.0 
H2O2+O=HO2+OH                  2.800E+13      0.0      6400.0 
H2O2+OH=HO2+H2O                2.000E+12      0.0       427.0 
    DUPLICATE
H2O2+OH=HO2+H2O                1.700E+18      0.0     29400.0 
    DUPLICATE 
END

\\
\\
\\  This is the therm file
\\
\\
THERMO ALL
 300.000  1000.000  5000.000
H                 L 6/94H   10   00   00   0G   200.000  6000.00  1000.0       1
 0.25000000E+01 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
 0.25473660E+05-0.44668285E+00 0.25000000E+01 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 0.25473660E+05-0.44668285E+00 0.26219035E+05    4
H2  REF ELEMENT   RUS 78H   20   00   00   0G   200.000  6000.00  1000.0       1
 0.29328305E+01 0.82659802E-03-0.14640057E-06 0.15409851E-10-0.68879615E-15    2
-0.81305582E+03-0.10243164E+01 0.23443029E+01 0.79804248E-02-0.19477917E-04    3
 0.20156967E-07-0.73760289E-11-0.91792413E+03 0.68300218E+00 0.00000000E+00    4
O                 L 1/90O   10   00   00   0G   200.000  6000.00  1000.0       1
 2.54363697E+00-2.73162486E-05-4.19029520E-09 4.95481845E-12-4.79553694E-16    2
 2.92260120E+04 4.92229457E+00 3.16826710E+00-3.27931884E-03 6.64306396E-06    3
-6.12806624E-09 2.11265971E-12 2.91222592E+04 2.05193346E+00 2.99687009E+04    4
O2 REF ELEMENT    RUS 89O   20   00   00   0G   200.000  6000.00  1000.0       1
 3.66096083E+00 6.56365523E-04-1.41149485E-07 2.05797658E-11-1.29913248E-15    2
-1.21597725E+03 3.41536184E+00 3.78245636E+00-2.99673415E-03 9.84730200E-06    3
-9.68129508E-09 3.24372836E-12-1.06394356E+03 3.65767573E+00 0.00000000E+00    4
OH                      O   1H   10   00   0G   200.000  6000.00  1000.0       1
 2.83853033E+00 1.10741289E-03-2.94000209E-07 4.20698729E-11-2.42289890E-15    2
 3.69780808E+03 5.84494652E+00 3.99198424E+00-2.40106655E-03 4.61664033E-06    3
-3.87916306E-09 1.36319502E-12 3.36889836E+03-1.03998477E-01 4.48615380E+03    4
HO2               L 5/89H   1O   20   00   0G   200.000  6000.00  1000.0       1
 0.41722659E+01 0.18812098E-02-0.34629297E-06 0.19468516E-10 0.17609153E-15    2
 0.61818851E+02 0.29577974E+01 0.43017880E+01-0.47490201E-02 0.21157953E-04    3
-0.24275961E-07 0.92920670E-11 0.29480876E+03 0.37167010E+01 0.15096500E+04    4
H2O               L 5/89H   2O   10   00   0G   200.000  6000.00  1000.0       1
 0.26770389E+01 0.29731816E-02-0.77376889E-06 0.94433514E-10-0.42689991E-14    2
-0.29885894E+05 0.68825500E+01 0.41986352E+01-0.20364017E-02 0.65203416E-05    3
-0.54879269E-08 0.17719680E-11-0.30293726E+05-0.84900901E+00-0.29084817E+05    4
H2O2              L 2/93H   2O   20   00   0G   200.000  6000.00  1000.0       1
 4.57977305E+00 4.05326003E-03-1.29844730E-06 1.98211400E-10-1.13968792E-14    2
-1.80071775E+04 6.64970694E-01 4.31515149E+00-8.47390622E-04 1.76404323E-05    3
-2.26762944E-08 9.08950158E-12-1.77067437E+04 3.27373319E+00-1.63425145E+04    4
AR REF ELEMENT    L 6/88AR  1    0    0    0G   200.000  6000.00  1000.0       1
 0.25000000E+01 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-0.74537500E+03 0.43796749E+01 0.25000000E+01 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-0.74537500E+03 0.43796749E+01 0.00000000E+00    4
N2  REF ELEMENT   8/02  N   2    0    0    0G   200.000  6000.00  1000.0       1
 2.95257637E+00 1.39690040E-03-4.92631603E-07 7.86010195E-11-4.60755204E-15    2
-9.23948688E+02 5.87188762E+00 3.53100528E+00-1.23660988E-04-5.02999433E-07    3
 2.43530612E-09-1.40881235E-12-1.04697628E+03 2.96747038E+00 0.00000000E+00    4
\\
\\
\\  This is the tran file
\\
\\
H                  0   145.000     2.050     0.000     0.000     0.000
H2                 1    38.000     2.920     0.000     0.790   280.000
O                  0    80.000     2.750     0.000     0.000     0.000
OH                 1    80.000     2.750     0.000     0.000     0.000
H2O                2   572.400     2.605     1.844     0.000     4.000
O2                 1   107.400     3.458     0.000     1.600     3.800
HO2                2   107.400     3.458     0.000     0.000     1.000  !(*)
H2O2               2   107.400     3.458     0.000     0.000     3.800
AR                 0   136.500     3.330     0.000     0.000     0.000
N2                 1    97.530     3.621     0.000     1.760     4.000

#endif
