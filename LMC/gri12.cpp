/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *                                   Marc Day
 *                    Lawrence Berkeley National Laboratory
 *                      (C) 1998-2003  All Rights Reserved
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
    *mm = 5;
    *kk = 32;
    *ii = 177;
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

    for (i=0; i<lenkname*5; i++) {
        kname[i] = ' ';
    }

    /* O  */
    kname[ 0*lenkname + 0 ] = 'O';
    kname[ 0*lenkname + 1 ] = ' ';

    /* H  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = ' ';

    /* C  */
    kname[ 2*lenkname + 0 ] = 'C';
    kname[ 2*lenkname + 1 ] = ' ';

    /* N  */
    kname[ 3*lenkname + 0 ] = 'N';
    kname[ 3*lenkname + 1 ] = ' ';

    /* AR  */
    kname[ 4*lenkname + 0 ] = 'A';
    kname[ 4*lenkname + 1 ] = 'R';
}


/* Returns the char strings of species names */
void CKSYMS(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*32; i++) {
        kname[i] = ' ';
    }

    /* H2  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = '2';
    kname[ 0*lenkname + 2 ] = ' ';

    /* H  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = ' ';

    /* O  */
    kname[ 2*lenkname + 0 ] = 'O';
    kname[ 2*lenkname + 1 ] = ' ';

    /* O2  */
    kname[ 3*lenkname + 0 ] = 'O';
    kname[ 3*lenkname + 1 ] = '2';
    kname[ 3*lenkname + 2 ] = ' ';

    /* OH  */
    kname[ 4*lenkname + 0 ] = 'O';
    kname[ 4*lenkname + 1 ] = 'H';
    kname[ 4*lenkname + 2 ] = ' ';

    /* H2O  */
    kname[ 5*lenkname + 0 ] = 'H';
    kname[ 5*lenkname + 1 ] = '2';
    kname[ 5*lenkname + 2 ] = 'O';
    kname[ 5*lenkname + 3 ] = ' ';

    /* HO2  */
    kname[ 6*lenkname + 0 ] = 'H';
    kname[ 6*lenkname + 1 ] = 'O';
    kname[ 6*lenkname + 2 ] = '2';
    kname[ 6*lenkname + 3 ] = ' ';

    /* H2O2  */
    kname[ 7*lenkname + 0 ] = 'H';
    kname[ 7*lenkname + 1 ] = '2';
    kname[ 7*lenkname + 2 ] = 'O';
    kname[ 7*lenkname + 3 ] = '2';
    kname[ 7*lenkname + 4 ] = ' ';

    /* C  */
    kname[ 8*lenkname + 0 ] = 'C';
    kname[ 8*lenkname + 1 ] = ' ';

    /* CH  */
    kname[ 9*lenkname + 0 ] = 'C';
    kname[ 9*lenkname + 1 ] = 'H';
    kname[ 9*lenkname + 2 ] = ' ';

    /* CH2  */
    kname[ 10*lenkname + 0 ] = 'C';
    kname[ 10*lenkname + 1 ] = 'H';
    kname[ 10*lenkname + 2 ] = '2';
    kname[ 10*lenkname + 3 ] = ' ';

    /* CH2(S)  */
    kname[ 11*lenkname + 0 ] = 'C';
    kname[ 11*lenkname + 1 ] = 'H';
    kname[ 11*lenkname + 2 ] = '2';
    kname[ 11*lenkname + 3 ] = '(';
    kname[ 11*lenkname + 4 ] = 'S';
    kname[ 11*lenkname + 5 ] = ')';
    kname[ 11*lenkname + 6 ] = ' ';

    /* CH3  */
    kname[ 12*lenkname + 0 ] = 'C';
    kname[ 12*lenkname + 1 ] = 'H';
    kname[ 12*lenkname + 2 ] = '3';
    kname[ 12*lenkname + 3 ] = ' ';

    /* CH4  */
    kname[ 13*lenkname + 0 ] = 'C';
    kname[ 13*lenkname + 1 ] = 'H';
    kname[ 13*lenkname + 2 ] = '4';
    kname[ 13*lenkname + 3 ] = ' ';

    /* CO  */
    kname[ 14*lenkname + 0 ] = 'C';
    kname[ 14*lenkname + 1 ] = 'O';
    kname[ 14*lenkname + 2 ] = ' ';

    /* CO2  */
    kname[ 15*lenkname + 0 ] = 'C';
    kname[ 15*lenkname + 1 ] = 'O';
    kname[ 15*lenkname + 2 ] = '2';
    kname[ 15*lenkname + 3 ] = ' ';

    /* HCO  */
    kname[ 16*lenkname + 0 ] = 'H';
    kname[ 16*lenkname + 1 ] = 'C';
    kname[ 16*lenkname + 2 ] = 'O';
    kname[ 16*lenkname + 3 ] = ' ';

    /* CH2O  */
    kname[ 17*lenkname + 0 ] = 'C';
    kname[ 17*lenkname + 1 ] = 'H';
    kname[ 17*lenkname + 2 ] = '2';
    kname[ 17*lenkname + 3 ] = 'O';
    kname[ 17*lenkname + 4 ] = ' ';

    /* CH2OH  */
    kname[ 18*lenkname + 0 ] = 'C';
    kname[ 18*lenkname + 1 ] = 'H';
    kname[ 18*lenkname + 2 ] = '2';
    kname[ 18*lenkname + 3 ] = 'O';
    kname[ 18*lenkname + 4 ] = 'H';
    kname[ 18*lenkname + 5 ] = ' ';

    /* CH3O  */
    kname[ 19*lenkname + 0 ] = 'C';
    kname[ 19*lenkname + 1 ] = 'H';
    kname[ 19*lenkname + 2 ] = '3';
    kname[ 19*lenkname + 3 ] = 'O';
    kname[ 19*lenkname + 4 ] = ' ';

    /* CH3OH  */
    kname[ 20*lenkname + 0 ] = 'C';
    kname[ 20*lenkname + 1 ] = 'H';
    kname[ 20*lenkname + 2 ] = '3';
    kname[ 20*lenkname + 3 ] = 'O';
    kname[ 20*lenkname + 4 ] = 'H';
    kname[ 20*lenkname + 5 ] = ' ';

    /* C2H  */
    kname[ 21*lenkname + 0 ] = 'C';
    kname[ 21*lenkname + 1 ] = '2';
    kname[ 21*lenkname + 2 ] = 'H';
    kname[ 21*lenkname + 3 ] = ' ';

    /* C2H2  */
    kname[ 22*lenkname + 0 ] = 'C';
    kname[ 22*lenkname + 1 ] = '2';
    kname[ 22*lenkname + 2 ] = 'H';
    kname[ 22*lenkname + 3 ] = '2';
    kname[ 22*lenkname + 4 ] = ' ';

    /* C2H3  */
    kname[ 23*lenkname + 0 ] = 'C';
    kname[ 23*lenkname + 1 ] = '2';
    kname[ 23*lenkname + 2 ] = 'H';
    kname[ 23*lenkname + 3 ] = '3';
    kname[ 23*lenkname + 4 ] = ' ';

    /* C2H4  */
    kname[ 24*lenkname + 0 ] = 'C';
    kname[ 24*lenkname + 1 ] = '2';
    kname[ 24*lenkname + 2 ] = 'H';
    kname[ 24*lenkname + 3 ] = '4';
    kname[ 24*lenkname + 4 ] = ' ';

    /* C2H5  */
    kname[ 25*lenkname + 0 ] = 'C';
    kname[ 25*lenkname + 1 ] = '2';
    kname[ 25*lenkname + 2 ] = 'H';
    kname[ 25*lenkname + 3 ] = '5';
    kname[ 25*lenkname + 4 ] = ' ';

    /* C2H6  */
    kname[ 26*lenkname + 0 ] = 'C';
    kname[ 26*lenkname + 1 ] = '2';
    kname[ 26*lenkname + 2 ] = 'H';
    kname[ 26*lenkname + 3 ] = '6';
    kname[ 26*lenkname + 4 ] = ' ';

    /* HCCO  */
    kname[ 27*lenkname + 0 ] = 'H';
    kname[ 27*lenkname + 1 ] = 'C';
    kname[ 27*lenkname + 2 ] = 'C';
    kname[ 27*lenkname + 3 ] = 'O';
    kname[ 27*lenkname + 4 ] = ' ';

    /* CH2CO  */
    kname[ 28*lenkname + 0 ] = 'C';
    kname[ 28*lenkname + 1 ] = 'H';
    kname[ 28*lenkname + 2 ] = '2';
    kname[ 28*lenkname + 3 ] = 'C';
    kname[ 28*lenkname + 4 ] = 'O';
    kname[ 28*lenkname + 5 ] = ' ';

    /* HCCOH  */
    kname[ 29*lenkname + 0 ] = 'H';
    kname[ 29*lenkname + 1 ] = 'C';
    kname[ 29*lenkname + 2 ] = 'C';
    kname[ 29*lenkname + 3 ] = 'O';
    kname[ 29*lenkname + 4 ] = 'H';
    kname[ 29*lenkname + 5 ] = ' ';

    /* N2  */
    kname[ 30*lenkname + 0 ] = 'N';
    kname[ 30*lenkname + 1 ] = '2';
    kname[ 30*lenkname + 2 ] = ' ';

    /* AR  */
    kname[ 31*lenkname + 0 ] = 'A';
    kname[ 31*lenkname + 1 ] = 'R';
    kname[ 31*lenkname + 2 ] = ' ';

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
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*12.011000; /*C */
    XW += x[9]*13.018970; /*CH */
    XW += x[10]*14.026940; /*CH2 */
    XW += x[11]*14.026940; /*CH2(S) */
    XW += x[12]*15.034910; /*CH3 */
    XW += x[13]*16.042880; /*CH4 */
    XW += x[14]*28.010400; /*CO */
    XW += x[15]*44.009800; /*CO2 */
    XW += x[16]*29.018370; /*HCO */
    XW += x[17]*30.026340; /*CH2O */
    XW += x[18]*31.034310; /*CH2OH */
    XW += x[19]*31.034310; /*CH3O */
    XW += x[20]*32.042280; /*CH3OH */
    XW += x[21]*25.029970; /*C2H */
    XW += x[22]*26.037940; /*C2H2 */
    XW += x[23]*27.045910; /*C2H3 */
    XW += x[24]*28.053880; /*C2H4 */
    XW += x[25]*29.061850; /*C2H5 */
    XW += x[26]*30.069820; /*C2H6 */
    XW += x[27]*41.029370; /*HCCO */
    XW += x[28]*42.037340; /*CH2CO */
    XW += x[29]*42.037340; /*HCCOH */
    XW += x[30]*28.013400; /*N2 */
    XW += x[31]*39.948000; /*AR */
    *P = *rho * 8.314e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void CKPY(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]/2.015940; /*H2 */
    YOW += y[1]/1.007970; /*H */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/18.015340; /*H2O */
    YOW += y[6]/33.006770; /*HO2 */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/12.011000; /*C */
    YOW += y[9]/13.018970; /*CH */
    YOW += y[10]/14.026940; /*CH2 */
    YOW += y[11]/14.026940; /*CH2(S) */
    YOW += y[12]/15.034910; /*CH3 */
    YOW += y[13]/16.042880; /*CH4 */
    YOW += y[14]/28.010400; /*CO */
    YOW += y[15]/44.009800; /*CO2 */
    YOW += y[16]/29.018370; /*HCO */
    YOW += y[17]/30.026340; /*CH2O */
    YOW += y[18]/31.034310; /*CH2OH */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/32.042280; /*CH3OH */
    YOW += y[21]/25.029970; /*C2H */
    YOW += y[22]/26.037940; /*C2H2 */
    YOW += y[23]/27.045910; /*C2H3 */
    YOW += y[24]/28.053880; /*C2H4 */
    YOW += y[25]/29.061850; /*C2H5 */
    YOW += y[26]/30.069820; /*C2H6 */
    YOW += y[27]/41.029370; /*HCCO */
    YOW += y[28]/42.037340; /*CH2CO */
    YOW += y[29]/42.037340; /*HCCOH */
    YOW += y[30]/28.013400; /*N2 */
    YOW += y[31]/39.948000; /*AR */
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
    W += c[0]*2.015940; /*H2 */
    W += c[1]*1.007970; /*H */
    W += c[2]*15.999400; /*O */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*17.007370; /*OH */
    W += c[5]*18.015340; /*H2O */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*12.011000; /*C */
    W += c[9]*13.018970; /*CH */
    W += c[10]*14.026940; /*CH2 */
    W += c[11]*14.026940; /*CH2(S) */
    W += c[12]*15.034910; /*CH3 */
    W += c[13]*16.042880; /*CH4 */
    W += c[14]*28.010400; /*CO */
    W += c[15]*44.009800; /*CO2 */
    W += c[16]*29.018370; /*HCO */
    W += c[17]*30.026340; /*CH2O */
    W += c[18]*31.034310; /*CH2OH */
    W += c[19]*31.034310; /*CH3O */
    W += c[20]*32.042280; /*CH3OH */
    W += c[21]*25.029970; /*C2H */
    W += c[22]*26.037940; /*C2H2 */
    W += c[23]*27.045910; /*C2H3 */
    W += c[24]*28.053880; /*C2H4 */
    W += c[25]*29.061850; /*C2H5 */
    W += c[26]*30.069820; /*C2H6 */
    W += c[27]*41.029370; /*HCCO */
    W += c[28]*42.037340; /*CH2CO */
    W += c[29]*42.037340; /*HCCOH */
    W += c[30]*28.013400; /*N2 */
    W += c[31]*39.948000; /*AR */

    for (id = 0; id < 32; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.314e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double * P, double * T, double * x, int * iwrk, double * rwrk, double * rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*12.011000; /*C */
    XW += x[9]*13.018970; /*CH */
    XW += x[10]*14.026940; /*CH2 */
    XW += x[11]*14.026940; /*CH2(S) */
    XW += x[12]*15.034910; /*CH3 */
    XW += x[13]*16.042880; /*CH4 */
    XW += x[14]*28.010400; /*CO */
    XW += x[15]*44.009800; /*CO2 */
    XW += x[16]*29.018370; /*HCO */
    XW += x[17]*30.026340; /*CH2O */
    XW += x[18]*31.034310; /*CH2OH */
    XW += x[19]*31.034310; /*CH3O */
    XW += x[20]*32.042280; /*CH3OH */
    XW += x[21]*25.029970; /*C2H */
    XW += x[22]*26.037940; /*C2H2 */
    XW += x[23]*27.045910; /*C2H3 */
    XW += x[24]*28.053880; /*C2H4 */
    XW += x[25]*29.061850; /*C2H5 */
    XW += x[26]*30.069820; /*C2H6 */
    XW += x[27]*41.029370; /*HCCO */
    XW += x[28]*42.037340; /*CH2CO */
    XW += x[29]*42.037340; /*HCCOH */
    XW += x[30]*28.013400; /*N2 */
    XW += x[31]*39.948000; /*AR */
    *rho = *P * XW / (8.314e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double * P, double * T, double * y, int * iwrk, double * rwrk, double * rho)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]/2.015940; /*H2 */
    YOW += y[1]/1.007970; /*H */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/18.015340; /*H2O */
    YOW += y[6]/33.006770; /*HO2 */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/12.011000; /*C */
    YOW += y[9]/13.018970; /*CH */
    YOW += y[10]/14.026940; /*CH2 */
    YOW += y[11]/14.026940; /*CH2(S) */
    YOW += y[12]/15.034910; /*CH3 */
    YOW += y[13]/16.042880; /*CH4 */
    YOW += y[14]/28.010400; /*CO */
    YOW += y[15]/44.009800; /*CO2 */
    YOW += y[16]/29.018370; /*HCO */
    YOW += y[17]/30.026340; /*CH2O */
    YOW += y[18]/31.034310; /*CH2OH */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/32.042280; /*CH3OH */
    YOW += y[21]/25.029970; /*C2H */
    YOW += y[22]/26.037940; /*C2H2 */
    YOW += y[23]/27.045910; /*C2H3 */
    YOW += y[24]/28.053880; /*C2H4 */
    YOW += y[25]/29.061850; /*C2H5 */
    YOW += y[26]/30.069820; /*C2H6 */
    YOW += y[27]/41.029370; /*HCCO */
    YOW += y[28]/42.037340; /*CH2CO */
    YOW += y[29]/42.037340; /*HCCOH */
    YOW += y[30]/28.013400; /*N2 */
    YOW += y[31]/39.948000; /*AR */
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
    W += c[0]*2.015940; /*H2 */
    W += c[1]*1.007970; /*H */
    W += c[2]*15.999400; /*O */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*17.007370; /*OH */
    W += c[5]*18.015340; /*H2O */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*12.011000; /*C */
    W += c[9]*13.018970; /*CH */
    W += c[10]*14.026940; /*CH2 */
    W += c[11]*14.026940; /*CH2(S) */
    W += c[12]*15.034910; /*CH3 */
    W += c[13]*16.042880; /*CH4 */
    W += c[14]*28.010400; /*CO */
    W += c[15]*44.009800; /*CO2 */
    W += c[16]*29.018370; /*HCO */
    W += c[17]*30.026340; /*CH2O */
    W += c[18]*31.034310; /*CH2OH */
    W += c[19]*31.034310; /*CH3O */
    W += c[20]*32.042280; /*CH3OH */
    W += c[21]*25.029970; /*C2H */
    W += c[22]*26.037940; /*C2H2 */
    W += c[23]*27.045910; /*C2H3 */
    W += c[24]*28.053880; /*C2H4 */
    W += c[25]*29.061850; /*C2H5 */
    W += c[26]*30.069820; /*C2H6 */
    W += c[27]*41.029370; /*HCCO */
    W += c[28]*42.037340; /*CH2CO */
    W += c[29]*42.037340; /*HCCOH */
    W += c[30]*28.013400; /*N2 */
    W += c[31]*39.948000; /*AR */

    for (id = 0; id < 32; ++id) {
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
    YOW += y[0]/2.015940; /*H2 */
    YOW += y[1]/1.007970; /*H */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/18.015340; /*H2O */
    YOW += y[6]/33.006770; /*HO2 */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/12.011000; /*C */
    YOW += y[9]/13.018970; /*CH */
    YOW += y[10]/14.026940; /*CH2 */
    YOW += y[11]/14.026940; /*CH2(S) */
    YOW += y[12]/15.034910; /*CH3 */
    YOW += y[13]/16.042880; /*CH4 */
    YOW += y[14]/28.010400; /*CO */
    YOW += y[15]/44.009800; /*CO2 */
    YOW += y[16]/29.018370; /*HCO */
    YOW += y[17]/30.026340; /*CH2O */
    YOW += y[18]/31.034310; /*CH2OH */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/32.042280; /*CH3OH */
    YOW += y[21]/25.029970; /*C2H */
    YOW += y[22]/26.037940; /*C2H2 */
    YOW += y[23]/27.045910; /*C2H3 */
    YOW += y[24]/28.053880; /*C2H4 */
    YOW += y[25]/29.061850; /*C2H5 */
    YOW += y[26]/30.069820; /*C2H6 */
    YOW += y[27]/41.029370; /*HCCO */
    YOW += y[28]/42.037340; /*CH2CO */
    YOW += y[29]/42.037340; /*HCCOH */
    YOW += y[30]/28.013400; /*N2 */
    YOW += y[31]/39.948000; /*AR */
    *wtm = 1.0 / YOW;

    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWX(double *x, int * iwrk, double * rwrk, double * wtm)
{
    double XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*12.011000; /*C */
    XW += x[9]*13.018970; /*CH */
    XW += x[10]*14.026940; /*CH2 */
    XW += x[11]*14.026940; /*CH2(S) */
    XW += x[12]*15.034910; /*CH3 */
    XW += x[13]*16.042880; /*CH4 */
    XW += x[14]*28.010400; /*CO */
    XW += x[15]*44.009800; /*CO2 */
    XW += x[16]*29.018370; /*HCO */
    XW += x[17]*30.026340; /*CH2O */
    XW += x[18]*31.034310; /*CH2OH */
    XW += x[19]*31.034310; /*CH3O */
    XW += x[20]*32.042280; /*CH3OH */
    XW += x[21]*25.029970; /*C2H */
    XW += x[22]*26.037940; /*C2H2 */
    XW += x[23]*27.045910; /*C2H3 */
    XW += x[24]*28.053880; /*C2H4 */
    XW += x[25]*29.061850; /*C2H5 */
    XW += x[26]*30.069820; /*C2H6 */
    XW += x[27]*41.029370; /*HCCO */
    XW += x[28]*42.037340; /*CH2CO */
    XW += x[29]*42.037340; /*HCCOH */
    XW += x[30]*28.013400; /*N2 */
    XW += x[31]*39.948000; /*AR */
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
    W += c[0]*2.015940; /*H2 */
    W += c[1]*1.007970; /*H */
    W += c[2]*15.999400; /*O */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*17.007370; /*OH */
    W += c[5]*18.015340; /*H2O */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*12.011000; /*C */
    W += c[9]*13.018970; /*CH */
    W += c[10]*14.026940; /*CH2 */
    W += c[11]*14.026940; /*CH2(S) */
    W += c[12]*15.034910; /*CH3 */
    W += c[13]*16.042880; /*CH4 */
    W += c[14]*28.010400; /*CO */
    W += c[15]*44.009800; /*CO2 */
    W += c[16]*29.018370; /*HCO */
    W += c[17]*30.026340; /*CH2O */
    W += c[18]*31.034310; /*CH2OH */
    W += c[19]*31.034310; /*CH3O */
    W += c[20]*32.042280; /*CH3OH */
    W += c[21]*25.029970; /*C2H */
    W += c[22]*26.037940; /*C2H2 */
    W += c[23]*27.045910; /*C2H3 */
    W += c[24]*28.053880; /*C2H4 */
    W += c[25]*29.061850; /*C2H5 */
    W += c[26]*30.069820; /*C2H6 */
    W += c[27]*41.029370; /*HCCO */
    W += c[28]*42.037340; /*CH2CO */
    W += c[29]*42.037340; /*HCCOH */
    W += c[30]*28.013400; /*N2 */
    W += c[31]*39.948000; /*AR */

    for (id = 0; id < 32; ++id) {
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
    YOW += y[0]/2.015940; /*H2 */
    YOW += y[1]/1.007970; /*H */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/18.015340; /*H2O */
    YOW += y[6]/33.006770; /*HO2 */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/12.011000; /*C */
    YOW += y[9]/13.018970; /*CH */
    YOW += y[10]/14.026940; /*CH2 */
    YOW += y[11]/14.026940; /*CH2(S) */
    YOW += y[12]/15.034910; /*CH3 */
    YOW += y[13]/16.042880; /*CH4 */
    YOW += y[14]/28.010400; /*CO */
    YOW += y[15]/44.009800; /*CO2 */
    YOW += y[16]/29.018370; /*HCO */
    YOW += y[17]/30.026340; /*CH2O */
    YOW += y[18]/31.034310; /*CH2OH */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/32.042280; /*CH3OH */
    YOW += y[21]/25.029970; /*C2H */
    YOW += y[22]/26.037940; /*C2H2 */
    YOW += y[23]/27.045910; /*C2H3 */
    YOW += y[24]/28.053880; /*C2H4 */
    YOW += y[25]/29.061850; /*C2H5 */
    YOW += y[26]/30.069820; /*C2H6 */
    YOW += y[27]/41.029370; /*HCCO */
    YOW += y[28]/42.037340; /*CH2CO */
    YOW += y[29]/42.037340; /*HCCOH */
    YOW += y[30]/28.013400; /*N2 */
    YOW += y[31]/39.948000; /*AR */
    /*Now compute conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(12.011000*YOW); 
    x[9] = y[9]/(13.018970*YOW); 
    x[10] = y[10]/(14.026940*YOW); 
    x[11] = y[11]/(14.026940*YOW); 
    x[12] = y[12]/(15.034910*YOW); 
    x[13] = y[13]/(16.042880*YOW); 
    x[14] = y[14]/(28.010400*YOW); 
    x[15] = y[15]/(44.009800*YOW); 
    x[16] = y[16]/(29.018370*YOW); 
    x[17] = y[17]/(30.026340*YOW); 
    x[18] = y[18]/(31.034310*YOW); 
    x[19] = y[19]/(31.034310*YOW); 
    x[20] = y[20]/(32.042280*YOW); 
    x[21] = y[21]/(25.029970*YOW); 
    x[22] = y[22]/(26.037940*YOW); 
    x[23] = y[23]/(27.045910*YOW); 
    x[24] = y[24]/(28.053880*YOW); 
    x[25] = y[25]/(29.061850*YOW); 
    x[26] = y[26]/(30.069820*YOW); 
    x[27] = y[27]/(41.029370*YOW); 
    x[28] = y[28]/(42.037340*YOW); 
    x[29] = y[29]/(42.037340*YOW); 
    x[30] = y[30]/(28.013400*YOW); 
    x[31] = y[31]/(39.948000*YOW); 

    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * c)
{
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/2.015940; /*H2 */
    YOW += y[1]/1.007970; /*H */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/18.015340; /*H2O */
    YOW += y[6]/33.006770; /*HO2 */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/12.011000; /*C */
    YOW += y[9]/13.018970; /*CH */
    YOW += y[10]/14.026940; /*CH2 */
    YOW += y[11]/14.026940; /*CH2(S) */
    YOW += y[12]/15.034910; /*CH3 */
    YOW += y[13]/16.042880; /*CH4 */
    YOW += y[14]/28.010400; /*CO */
    YOW += y[15]/44.009800; /*CO2 */
    YOW += y[16]/29.018370; /*HCO */
    YOW += y[17]/30.026340; /*CH2O */
    YOW += y[18]/31.034310; /*CH2OH */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/32.042280; /*CH3OH */
    YOW += y[21]/25.029970; /*C2H */
    YOW += y[22]/26.037940; /*C2H2 */
    YOW += y[23]/27.045910; /*C2H3 */
    YOW += y[24]/28.053880; /*C2H4 */
    YOW += y[25]/29.061850; /*C2H5 */
    YOW += y[26]/30.069820; /*C2H6 */
    YOW += y[27]/41.029370; /*HCCO */
    YOW += y[28]/42.037340; /*CH2CO */
    YOW += y[29]/42.037340; /*HCCOH */
    YOW += y[30]/28.013400; /*N2 */
    YOW += y[31]/39.948000; /*AR */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.314e+07 * (*T)); 
    /*Now compute conversion */
    c[0] = PWORT * y[0]/2.015940; 
    c[1] = PWORT * y[1]/1.007970; 
    c[2] = PWORT * y[2]/15.999400; 
    c[3] = PWORT * y[3]/31.998800; 
    c[4] = PWORT * y[4]/17.007370; 
    c[5] = PWORT * y[5]/18.015340; 
    c[6] = PWORT * y[6]/33.006770; 
    c[7] = PWORT * y[7]/34.014740; 
    c[8] = PWORT * y[8]/12.011000; 
    c[9] = PWORT * y[9]/13.018970; 
    c[10] = PWORT * y[10]/14.026940; 
    c[11] = PWORT * y[11]/14.026940; 
    c[12] = PWORT * y[12]/15.034910; 
    c[13] = PWORT * y[13]/16.042880; 
    c[14] = PWORT * y[14]/28.010400; 
    c[15] = PWORT * y[15]/44.009800; 
    c[16] = PWORT * y[16]/29.018370; 
    c[17] = PWORT * y[17]/30.026340; 
    c[18] = PWORT * y[18]/31.034310; 
    c[19] = PWORT * y[19]/31.034310; 
    c[20] = PWORT * y[20]/32.042280; 
    c[21] = PWORT * y[21]/25.029970; 
    c[22] = PWORT * y[22]/26.037940; 
    c[23] = PWORT * y[23]/27.045910; 
    c[24] = PWORT * y[24]/28.053880; 
    c[25] = PWORT * y[25]/29.061850; 
    c[26] = PWORT * y[26]/30.069820; 
    c[27] = PWORT * y[27]/41.029370; 
    c[28] = PWORT * y[28]/42.037340; 
    c[29] = PWORT * y[29]/42.037340; 
    c[30] = PWORT * y[30]/28.013400; 
    c[31] = PWORT * y[31]/39.948000; 

    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * c)
{
    /*See Eq 8 (Temperature not used) */
    c[0] = (*rho) * y[0]/2.015940; 
    c[1] = (*rho) * y[1]/1.007970; 
    c[2] = (*rho) * y[2]/15.999400; 
    c[3] = (*rho) * y[3]/31.998800; 
    c[4] = (*rho) * y[4]/17.007370; 
    c[5] = (*rho) * y[5]/18.015340; 
    c[6] = (*rho) * y[6]/33.006770; 
    c[7] = (*rho) * y[7]/34.014740; 
    c[8] = (*rho) * y[8]/12.011000; 
    c[9] = (*rho) * y[9]/13.018970; 
    c[10] = (*rho) * y[10]/14.026940; 
    c[11] = (*rho) * y[11]/14.026940; 
    c[12] = (*rho) * y[12]/15.034910; 
    c[13] = (*rho) * y[13]/16.042880; 
    c[14] = (*rho) * y[14]/28.010400; 
    c[15] = (*rho) * y[15]/44.009800; 
    c[16] = (*rho) * y[16]/29.018370; 
    c[17] = (*rho) * y[17]/30.026340; 
    c[18] = (*rho) * y[18]/31.034310; 
    c[19] = (*rho) * y[19]/31.034310; 
    c[20] = (*rho) * y[20]/32.042280; 
    c[21] = (*rho) * y[21]/25.029970; 
    c[22] = (*rho) * y[22]/26.037940; 
    c[23] = (*rho) * y[23]/27.045910; 
    c[24] = (*rho) * y[24]/28.053880; 
    c[25] = (*rho) * y[25]/29.061850; 
    c[26] = (*rho) * y[26]/30.069820; 
    c[27] = (*rho) * y[27]/41.029370; 
    c[28] = (*rho) * y[28]/42.037340; 
    c[29] = (*rho) * y[29]/42.037340; 
    c[30] = (*rho) * y[30]/28.013400; 
    c[31] = (*rho) * y[31]/39.948000; 

    return;
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
void CKXTY(double * x, int * iwrk, double * rwrk, double * y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*12.011000; /*C */
    XW += x[9]*13.018970; /*CH */
    XW += x[10]*14.026940; /*CH2 */
    XW += x[11]*14.026940; /*CH2(S) */
    XW += x[12]*15.034910; /*CH3 */
    XW += x[13]*16.042880; /*CH4 */
    XW += x[14]*28.010400; /*CO */
    XW += x[15]*44.009800; /*CO2 */
    XW += x[16]*29.018370; /*HCO */
    XW += x[17]*30.026340; /*CH2O */
    XW += x[18]*31.034310; /*CH2OH */
    XW += x[19]*31.034310; /*CH3O */
    XW += x[20]*32.042280; /*CH3OH */
    XW += x[21]*25.029970; /*C2H */
    XW += x[22]*26.037940; /*C2H2 */
    XW += x[23]*27.045910; /*C2H3 */
    XW += x[24]*28.053880; /*C2H4 */
    XW += x[25]*29.061850; /*C2H5 */
    XW += x[26]*30.069820; /*C2H6 */
    XW += x[27]*41.029370; /*HCCO */
    XW += x[28]*42.037340; /*CH2CO */
    XW += x[29]*42.037340; /*HCCOH */
    XW += x[30]*28.013400; /*N2 */
    XW += x[31]*39.948000; /*AR */
    /*Now compute conversion */
    y[0] = x[0]*2.015940/XW; 
    y[1] = x[1]*1.007970/XW; 
    y[2] = x[2]*15.999400/XW; 
    y[3] = x[3]*31.998800/XW; 
    y[4] = x[4]*17.007370/XW; 
    y[5] = x[5]*18.015340/XW; 
    y[6] = x[6]*33.006770/XW; 
    y[7] = x[7]*34.014740/XW; 
    y[8] = x[8]*12.011000/XW; 
    y[9] = x[9]*13.018970/XW; 
    y[10] = x[10]*14.026940/XW; 
    y[11] = x[11]*14.026940/XW; 
    y[12] = x[12]*15.034910/XW; 
    y[13] = x[13]*16.042880/XW; 
    y[14] = x[14]*28.010400/XW; 
    y[15] = x[15]*44.009800/XW; 
    y[16] = x[16]*29.018370/XW; 
    y[17] = x[17]*30.026340/XW; 
    y[18] = x[18]*31.034310/XW; 
    y[19] = x[19]*31.034310/XW; 
    y[20] = x[20]*32.042280/XW; 
    y[21] = x[21]*25.029970/XW; 
    y[22] = x[22]*26.037940/XW; 
    y[23] = x[23]*27.045910/XW; 
    y[24] = x[24]*28.053880/XW; 
    y[25] = x[25]*29.061850/XW; 
    y[26] = x[26]*30.069820/XW; 
    y[27] = x[27]*41.029370/XW; 
    y[28] = x[28]*42.037340/XW; 
    y[29] = x[29]*42.037340/XW; 
    y[30] = x[30]*28.013400/XW; 
    y[31] = x[31]*39.948000/XW; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.314e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 32; ++id) {
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
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*12.011000; /*C */
    XW += x[9]*13.018970; /*CH */
    XW += x[10]*14.026940; /*CH2 */
    XW += x[11]*14.026940; /*CH2(S) */
    XW += x[12]*15.034910; /*CH3 */
    XW += x[13]*16.042880; /*CH4 */
    XW += x[14]*28.010400; /*CO */
    XW += x[15]*44.009800; /*CO2 */
    XW += x[16]*29.018370; /*HCO */
    XW += x[17]*30.026340; /*CH2O */
    XW += x[18]*31.034310; /*CH2OH */
    XW += x[19]*31.034310; /*CH3O */
    XW += x[20]*32.042280; /*CH3OH */
    XW += x[21]*25.029970; /*C2H */
    XW += x[22]*26.037940; /*C2H2 */
    XW += x[23]*27.045910; /*C2H3 */
    XW += x[24]*28.053880; /*C2H4 */
    XW += x[25]*29.061850; /*C2H5 */
    XW += x[26]*30.069820; /*C2H6 */
    XW += x[27]*41.029370; /*HCCO */
    XW += x[28]*42.037340; /*CH2CO */
    XW += x[29]*42.037340; /*HCCOH */
    XW += x[30]*28.013400; /*N2 */
    XW += x[31]*39.948000; /*AR */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 32; ++id) {
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
    for (id = 0; id < 32; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    for (id = 0; id < 32; ++id) {
        x[id] = c[id]/sumC;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
void CKCTY(double * c, int * iwrk, double * rwrk, double * y)
{
    double CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*2.015940; /*H2 */
    CW += c[1]*1.007970; /*H */
    CW += c[2]*15.999400; /*O */
    CW += c[3]*31.998800; /*O2 */
    CW += c[4]*17.007370; /*OH */
    CW += c[5]*18.015340; /*H2O */
    CW += c[6]*33.006770; /*HO2 */
    CW += c[7]*34.014740; /*H2O2 */
    CW += c[8]*12.011000; /*C */
    CW += c[9]*13.018970; /*CH */
    CW += c[10]*14.026940; /*CH2 */
    CW += c[11]*14.026940; /*CH2(S) */
    CW += c[12]*15.034910; /*CH3 */
    CW += c[13]*16.042880; /*CH4 */
    CW += c[14]*28.010400; /*CO */
    CW += c[15]*44.009800; /*CO2 */
    CW += c[16]*29.018370; /*HCO */
    CW += c[17]*30.026340; /*CH2O */
    CW += c[18]*31.034310; /*CH2OH */
    CW += c[19]*31.034310; /*CH3O */
    CW += c[20]*32.042280; /*CH3OH */
    CW += c[21]*25.029970; /*C2H */
    CW += c[22]*26.037940; /*C2H2 */
    CW += c[23]*27.045910; /*C2H3 */
    CW += c[24]*28.053880; /*C2H4 */
    CW += c[25]*29.061850; /*C2H5 */
    CW += c[26]*30.069820; /*C2H6 */
    CW += c[27]*41.029370; /*HCCO */
    CW += c[28]*42.037340; /*CH2CO */
    CW += c[29]*42.037340; /*HCCOH */
    CW += c[30]*28.013400; /*N2 */
    CW += c[31]*39.948000; /*AR */
    /*Now compute conversion */
    y[0] = c[0]*2.015940/CW; 
    y[1] = c[1]*1.007970/CW; 
    y[2] = c[2]*15.999400/CW; 
    y[3] = c[3]*31.998800/CW; 
    y[4] = c[4]*17.007370/CW; 
    y[5] = c[5]*18.015340/CW; 
    y[6] = c[6]*33.006770/CW; 
    y[7] = c[7]*34.014740/CW; 
    y[8] = c[8]*12.011000/CW; 
    y[9] = c[9]*13.018970/CW; 
    y[10] = c[10]*14.026940/CW; 
    y[11] = c[11]*14.026940/CW; 
    y[12] = c[12]*15.034910/CW; 
    y[13] = c[13]*16.042880/CW; 
    y[14] = c[14]*28.010400/CW; 
    y[15] = c[15]*44.009800/CW; 
    y[16] = c[16]*29.018370/CW; 
    y[17] = c[17]*30.026340/CW; 
    y[18] = c[18]*31.034310/CW; 
    y[19] = c[19]*31.034310/CW; 
    y[20] = c[20]*32.042280/CW; 
    y[21] = c[21]*25.029970/CW; 
    y[22] = c[22]*26.037940/CW; 
    y[23] = c[23]*27.045910/CW; 
    y[24] = c[24]*28.053880/CW; 
    y[25] = c[25]*29.061850/CW; 
    y[26] = c[26]*30.069820/CW; 
    y[27] = c[27]*41.029370/CW; 
    y[28] = c[28]*42.037340/CW; 
    y[29] = c[29]*42.037340/CW; 
    y[30] = c[30]*28.013400/CW; 
    y[31] = c[31]*39.948000/CW; 

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
    for (id = 0; id < 32; ++id) {
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
    for (id = 0; id < 32; ++id) {
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
    for (id = 0; id < 32; ++id) {
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
    for (id = 0; id < 32; ++id) {
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
    for (id = 0; id < 32; ++id) {
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
    for (id = 0; id < 32; ++id) {
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
    for (id = 0; id < 32; ++id) {
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
    cvms[0] *= 41241306.784924; /*H2 */
    cvms[1] *= 82482613.569848; /*H */
    cvms[2] *= 5196444.866683; /*O */
    cvms[3] *= 2598222.433341; /*O2 */
    cvms[4] *= 4888468.940230; /*OH */
    cvms[5] *= 4614955.920899; /*H2O */
    cvms[6] *= 2518877.187922; /*HO2 */
    cvms[7] *= 2444234.470115; /*H2O2 */
    cvms[8] *= 6921988.177504; /*C */
    cvms[9] *= 6386065.871570; /*CH */
    cvms[10] *= 5927165.867966; /*CH2 */
    cvms[11] *= 5927165.867966; /*CH2(S) */
    cvms[12] *= 5529796.985815; /*CH3 */
    cvms[13] *= 5182361.271792; /*CH4 */
    cvms[14] *= 2968183.246223; /*CO */
    cvms[15] *= 1889124.694954; /*CO2 */
    cvms[16] *= 2865081.670680; /*HCO */
    cvms[17] *= 2768902.237169; /*CH2O */
    cvms[18] *= 2678970.468491; /*CH2OH */
    cvms[19] *= 2678970.468491; /*CH3O */
    cvms[20] *= 2594696.756910; /*CH3OH */
    cvms[21] *= 3321618.044289; /*C2H */
    cvms[22] *= 3193032.935785; /*C2H2 */
    cvms[23] *= 3074032.265877; /*C2H3 */
    cvms[24] *= 2963582.933983; /*C2H4 */
    cvms[25] *= 2860795.166171; /*C2H5 */
    cvms[26] *= 2764898.492908; /*C2H6 */
    cvms[27] *= 2026353.317148; /*HCCO */
    cvms[28] *= 1977765.481831; /*CH2CO */
    cvms[29] *= 1977765.481831; /*HCCOH */
    cvms[30] *= 2967865.378712; /*N2 */
    cvms[31] *= 2081205.567237; /*AR */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
void CKCPMS(double *T, int * iwrk, double * rwrk, double * cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 41241306.784924; /*H2 */
    cpms[1] *= 82482613.569848; /*H */
    cpms[2] *= 5196444.866683; /*O */
    cpms[3] *= 2598222.433341; /*O2 */
    cpms[4] *= 4888468.940230; /*OH */
    cpms[5] *= 4614955.920899; /*H2O */
    cpms[6] *= 2518877.187922; /*HO2 */
    cpms[7] *= 2444234.470115; /*H2O2 */
    cpms[8] *= 6921988.177504; /*C */
    cpms[9] *= 6386065.871570; /*CH */
    cpms[10] *= 5927165.867966; /*CH2 */
    cpms[11] *= 5927165.867966; /*CH2(S) */
    cpms[12] *= 5529796.985815; /*CH3 */
    cpms[13] *= 5182361.271792; /*CH4 */
    cpms[14] *= 2968183.246223; /*CO */
    cpms[15] *= 1889124.694954; /*CO2 */
    cpms[16] *= 2865081.670680; /*HCO */
    cpms[17] *= 2768902.237169; /*CH2O */
    cpms[18] *= 2678970.468491; /*CH2OH */
    cpms[19] *= 2678970.468491; /*CH3O */
    cpms[20] *= 2594696.756910; /*CH3OH */
    cpms[21] *= 3321618.044289; /*C2H */
    cpms[22] *= 3193032.935785; /*C2H2 */
    cpms[23] *= 3074032.265877; /*C2H3 */
    cpms[24] *= 2963582.933983; /*C2H4 */
    cpms[25] *= 2860795.166171; /*C2H5 */
    cpms[26] *= 2764898.492908; /*C2H6 */
    cpms[27] *= 2026353.317148; /*HCCO */
    cpms[28] *= 1977765.481831; /*CH2CO */
    cpms[29] *= 1977765.481831; /*HCCOH */
    cpms[30] *= 2967865.378712; /*N2 */
    cpms[31] *= 2081205.567237; /*AR */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double *T, int * iwrk, double * rwrk, double * ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    ums[0] *= RT/2.015940; /*H2 */
    ums[1] *= RT/1.007970; /*H */
    ums[2] *= RT/15.999400; /*O */
    ums[3] *= RT/31.998800; /*O2 */
    ums[4] *= RT/17.007370; /*OH */
    ums[5] *= RT/18.015340; /*H2O */
    ums[6] *= RT/33.006770; /*HO2 */
    ums[7] *= RT/34.014740; /*H2O2 */
    ums[8] *= RT/12.011000; /*C */
    ums[9] *= RT/13.018970; /*CH */
    ums[10] *= RT/14.026940; /*CH2 */
    ums[11] *= RT/14.026940; /*CH2(S) */
    ums[12] *= RT/15.034910; /*CH3 */
    ums[13] *= RT/16.042880; /*CH4 */
    ums[14] *= RT/28.010400; /*CO */
    ums[15] *= RT/44.009800; /*CO2 */
    ums[16] *= RT/29.018370; /*HCO */
    ums[17] *= RT/30.026340; /*CH2O */
    ums[18] *= RT/31.034310; /*CH2OH */
    ums[19] *= RT/31.034310; /*CH3O */
    ums[20] *= RT/32.042280; /*CH3OH */
    ums[21] *= RT/25.029970; /*C2H */
    ums[22] *= RT/26.037940; /*C2H2 */
    ums[23] *= RT/27.045910; /*C2H3 */
    ums[24] *= RT/28.053880; /*C2H4 */
    ums[25] *= RT/29.061850; /*C2H5 */
    ums[26] *= RT/30.069820; /*C2H6 */
    ums[27] *= RT/41.029370; /*HCCO */
    ums[28] *= RT/42.037340; /*CH2CO */
    ums[29] *= RT/42.037340; /*HCCOH */
    ums[30] *= RT/28.013400; /*N2 */
    ums[31] *= RT/39.948000; /*AR */
}


/*Returns enthalpy in mass units (Eq 27.) */
void CKHMS(double *T, int * iwrk, double * rwrk, double * hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    speciesEnthalpy(hms, tc);
    hms[0] *= RT/2.015940; /*H2 */
    hms[1] *= RT/1.007970; /*H */
    hms[2] *= RT/15.999400; /*O */
    hms[3] *= RT/31.998800; /*O2 */
    hms[4] *= RT/17.007370; /*OH */
    hms[5] *= RT/18.015340; /*H2O */
    hms[6] *= RT/33.006770; /*HO2 */
    hms[7] *= RT/34.014740; /*H2O2 */
    hms[8] *= RT/12.011000; /*C */
    hms[9] *= RT/13.018970; /*CH */
    hms[10] *= RT/14.026940; /*CH2 */
    hms[11] *= RT/14.026940; /*CH2(S) */
    hms[12] *= RT/15.034910; /*CH3 */
    hms[13] *= RT/16.042880; /*CH4 */
    hms[14] *= RT/28.010400; /*CO */
    hms[15] *= RT/44.009800; /*CO2 */
    hms[16] *= RT/29.018370; /*HCO */
    hms[17] *= RT/30.026340; /*CH2O */
    hms[18] *= RT/31.034310; /*CH2OH */
    hms[19] *= RT/31.034310; /*CH3O */
    hms[20] *= RT/32.042280; /*CH3OH */
    hms[21] *= RT/25.029970; /*C2H */
    hms[22] *= RT/26.037940; /*C2H2 */
    hms[23] *= RT/27.045910; /*C2H3 */
    hms[24] *= RT/28.053880; /*C2H4 */
    hms[25] *= RT/29.061850; /*C2H5 */
    hms[26] *= RT/30.069820; /*C2H6 */
    hms[27] *= RT/41.029370; /*HCCO */
    hms[28] *= RT/42.037340; /*CH2CO */
    hms[29] *= RT/42.037340; /*HCCOH */
    hms[30] *= RT/28.013400; /*N2 */
    hms[31] *= RT/39.948000; /*AR */
}


/*Returns gibbs in mass units (Eq 31.) */
void CKGMS(double *T, int * iwrk, double * rwrk, double * gms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    gibbs(gms, tc);
    gms[0] *= RT/2.015940; /*H2 */
    gms[1] *= RT/1.007970; /*H */
    gms[2] *= RT/15.999400; /*O */
    gms[3] *= RT/31.998800; /*O2 */
    gms[4] *= RT/17.007370; /*OH */
    gms[5] *= RT/18.015340; /*H2O */
    gms[6] *= RT/33.006770; /*HO2 */
    gms[7] *= RT/34.014740; /*H2O2 */
    gms[8] *= RT/12.011000; /*C */
    gms[9] *= RT/13.018970; /*CH */
    gms[10] *= RT/14.026940; /*CH2 */
    gms[11] *= RT/14.026940; /*CH2(S) */
    gms[12] *= RT/15.034910; /*CH3 */
    gms[13] *= RT/16.042880; /*CH4 */
    gms[14] *= RT/28.010400; /*CO */
    gms[15] *= RT/44.009800; /*CO2 */
    gms[16] *= RT/29.018370; /*HCO */
    gms[17] *= RT/30.026340; /*CH2O */
    gms[18] *= RT/31.034310; /*CH2OH */
    gms[19] *= RT/31.034310; /*CH3O */
    gms[20] *= RT/32.042280; /*CH3OH */
    gms[21] *= RT/25.029970; /*C2H */
    gms[22] *= RT/26.037940; /*C2H2 */
    gms[23] *= RT/27.045910; /*C2H3 */
    gms[24] *= RT/28.053880; /*C2H4 */
    gms[25] *= RT/29.061850; /*C2H5 */
    gms[26] *= RT/30.069820; /*C2H6 */
    gms[27] *= RT/41.029370; /*HCCO */
    gms[28] *= RT/42.037340; /*CH2CO */
    gms[29] *= RT/42.037340; /*HCCOH */
    gms[30] *= RT/28.013400; /*N2 */
    gms[31] *= RT/39.948000; /*AR */
}


/*Returns helmholtz in mass units (Eq 32.) */
void CKAMS(double *T, int * iwrk, double * rwrk, double * ams)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    helmholtz(ams, tc);
    ams[0] *= RT/2.015940; /*H2 */
    ams[1] *= RT/1.007970; /*H */
    ams[2] *= RT/15.999400; /*O */
    ams[3] *= RT/31.998800; /*O2 */
    ams[4] *= RT/17.007370; /*OH */
    ams[5] *= RT/18.015340; /*H2O */
    ams[6] *= RT/33.006770; /*HO2 */
    ams[7] *= RT/34.014740; /*H2O2 */
    ams[8] *= RT/12.011000; /*C */
    ams[9] *= RT/13.018970; /*CH */
    ams[10] *= RT/14.026940; /*CH2 */
    ams[11] *= RT/14.026940; /*CH2(S) */
    ams[12] *= RT/15.034910; /*CH3 */
    ams[13] *= RT/16.042880; /*CH4 */
    ams[14] *= RT/28.010400; /*CO */
    ams[15] *= RT/44.009800; /*CO2 */
    ams[16] *= RT/29.018370; /*HCO */
    ams[17] *= RT/30.026340; /*CH2O */
    ams[18] *= RT/31.034310; /*CH2OH */
    ams[19] *= RT/31.034310; /*CH3O */
    ams[20] *= RT/32.042280; /*CH3OH */
    ams[21] *= RT/25.029970; /*C2H */
    ams[22] *= RT/26.037940; /*C2H2 */
    ams[23] *= RT/27.045910; /*C2H3 */
    ams[24] *= RT/28.053880; /*C2H4 */
    ams[25] *= RT/29.061850; /*C2H5 */
    ams[26] *= RT/30.069820; /*C2H6 */
    ams[27] *= RT/41.029370; /*HCCO */
    ams[28] *= RT/42.037340; /*CH2CO */
    ams[29] *= RT/42.037340; /*HCCOH */
    ams[30] *= RT/28.013400; /*N2 */
    ams[31] *= RT/39.948000; /*AR */
}


/*Returns the entropies in mass units (Eq 28.) */
void CKSMS(double *T, int * iwrk, double * rwrk, double * sms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 41241306.784924; /*H2 */
    sms[1] *= 82482613.569848; /*H */
    sms[2] *= 5196444.866683; /*O */
    sms[3] *= 2598222.433341; /*O2 */
    sms[4] *= 4888468.940230; /*OH */
    sms[5] *= 4614955.920899; /*H2O */
    sms[6] *= 2518877.187922; /*HO2 */
    sms[7] *= 2444234.470115; /*H2O2 */
    sms[8] *= 6921988.177504; /*C */
    sms[9] *= 6386065.871570; /*CH */
    sms[10] *= 5927165.867966; /*CH2 */
    sms[11] *= 5927165.867966; /*CH2(S) */
    sms[12] *= 5529796.985815; /*CH3 */
    sms[13] *= 5182361.271792; /*CH4 */
    sms[14] *= 2968183.246223; /*CO */
    sms[15] *= 1889124.694954; /*CO2 */
    sms[16] *= 2865081.670680; /*HCO */
    sms[17] *= 2768902.237169; /*CH2O */
    sms[18] *= 2678970.468491; /*CH2OH */
    sms[19] *= 2678970.468491; /*CH3O */
    sms[20] *= 2594696.756910; /*CH3OH */
    sms[21] *= 3321618.044289; /*C2H */
    sms[22] *= 3193032.935785; /*C2H2 */
    sms[23] *= 3074032.265877; /*C2H3 */
    sms[24] *= 2963582.933983; /*C2H4 */
    sms[25] *= 2860795.166171; /*C2H5 */
    sms[26] *= 2764898.492908; /*C2H6 */
    sms[27] *= 2026353.317148; /*HCCO */
    sms[28] *= 1977765.481831; /*CH2CO */
    sms[29] *= 1977765.481831; /*HCCOH */
    sms[30] *= 2967865.378712; /*N2 */
    sms[31] *= 2081205.567237; /*AR */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double *T, double *x, int * iwrk, double * rwrk, double * cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[32]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 32; ++id) {
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
    double cpor[32]; /* temporary storage */
    cp_R(cpor, tc);
    /*multiply by y/molecularweight */
    result += cpor[0]*y[0]/2.01594; /*H2 */
    result += cpor[1]*y[1]/1.00797; /*H */
    result += cpor[2]*y[2]/15.9994; /*O */
    result += cpor[3]*y[3]/31.9988; /*O2 */
    result += cpor[4]*y[4]/17.0074; /*OH */
    result += cpor[5]*y[5]/18.0153; /*H2O */
    result += cpor[6]*y[6]/33.0068; /*HO2 */
    result += cpor[7]*y[7]/34.0147; /*H2O2 */
    result += cpor[8]*y[8]/12.011; /*C */
    result += cpor[9]*y[9]/13.019; /*CH */
    result += cpor[10]*y[10]/14.0269; /*CH2 */
    result += cpor[11]*y[11]/14.0269; /*CH2(S) */
    result += cpor[12]*y[12]/15.0349; /*CH3 */
    result += cpor[13]*y[13]/16.0429; /*CH4 */
    result += cpor[14]*y[14]/28.0104; /*CO */
    result += cpor[15]*y[15]/44.0098; /*CO2 */
    result += cpor[16]*y[16]/29.0184; /*HCO */
    result += cpor[17]*y[17]/30.0263; /*CH2O */
    result += cpor[18]*y[18]/31.0343; /*CH2OH */
    result += cpor[19]*y[19]/31.0343; /*CH3O */
    result += cpor[20]*y[20]/32.0423; /*CH3OH */
    result += cpor[21]*y[21]/25.03; /*C2H */
    result += cpor[22]*y[22]/26.0379; /*C2H2 */
    result += cpor[23]*y[23]/27.0459; /*C2H3 */
    result += cpor[24]*y[24]/28.0539; /*C2H4 */
    result += cpor[25]*y[25]/29.0618; /*C2H5 */
    result += cpor[26]*y[26]/30.0698; /*C2H6 */
    result += cpor[27]*y[27]/41.0294; /*HCCO */
    result += cpor[28]*y[28]/42.0373; /*CH2CO */
    result += cpor[29]*y[29]/42.0373; /*HCCOH */
    result += cpor[30]*y[30]/28.0134; /*N2 */
    result += cpor[31]*y[31]/39.948; /*AR */

    *cpbs = result * 8.314e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double *T, double *x, int * iwrk, double * rwrk, double * cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[32]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 32; ++id) {
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
    double cvor[32]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]/2.01594; /*H2 */
    result += cvor[1]*y[1]/1.00797; /*H */
    result += cvor[2]*y[2]/15.9994; /*O */
    result += cvor[3]*y[3]/31.9988; /*O2 */
    result += cvor[4]*y[4]/17.0074; /*OH */
    result += cvor[5]*y[5]/18.0153; /*H2O */
    result += cvor[6]*y[6]/33.0068; /*HO2 */
    result += cvor[7]*y[7]/34.0147; /*H2O2 */
    result += cvor[8]*y[8]/12.011; /*C */
    result += cvor[9]*y[9]/13.019; /*CH */
    result += cvor[10]*y[10]/14.0269; /*CH2 */
    result += cvor[11]*y[11]/14.0269; /*CH2(S) */
    result += cvor[12]*y[12]/15.0349; /*CH3 */
    result += cvor[13]*y[13]/16.0429; /*CH4 */
    result += cvor[14]*y[14]/28.0104; /*CO */
    result += cvor[15]*y[15]/44.0098; /*CO2 */
    result += cvor[16]*y[16]/29.0184; /*HCO */
    result += cvor[17]*y[17]/30.0263; /*CH2O */
    result += cvor[18]*y[18]/31.0343; /*CH2OH */
    result += cvor[19]*y[19]/31.0343; /*CH3O */
    result += cvor[20]*y[20]/32.0423; /*CH3OH */
    result += cvor[21]*y[21]/25.03; /*C2H */
    result += cvor[22]*y[22]/26.0379; /*C2H2 */
    result += cvor[23]*y[23]/27.0459; /*C2H3 */
    result += cvor[24]*y[24]/28.0539; /*C2H4 */
    result += cvor[25]*y[25]/29.0618; /*C2H5 */
    result += cvor[26]*y[26]/30.0698; /*C2H6 */
    result += cvor[27]*y[27]/41.0294; /*HCCO */
    result += cvor[28]*y[28]/42.0373; /*CH2CO */
    result += cvor[29]*y[29]/42.0373; /*HCCOH */
    result += cvor[30]*y[30]/28.0134; /*N2 */
    result += cvor[31]*y[31]/39.948; /*AR */

    *cvbs = result * 8.314e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double *T, double *x, int * iwrk, double * rwrk, double * hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[32]; /* temporary storage */
    double RT = 8.314e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 32; ++id) {
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
    double hml[32]; /* temporary storage */
    double RT = 8.314e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*hml[0]/2.015940; /*H2 */
    result += y[1]*hml[1]/1.007970; /*H */
    result += y[2]*hml[2]/15.999400; /*O */
    result += y[3]*hml[3]/31.998800; /*O2 */
    result += y[4]*hml[4]/17.007370; /*OH */
    result += y[5]*hml[5]/18.015340; /*H2O */
    result += y[6]*hml[6]/33.006770; /*HO2 */
    result += y[7]*hml[7]/34.014740; /*H2O2 */
    result += y[8]*hml[8]/12.011000; /*C */
    result += y[9]*hml[9]/13.018970; /*CH */
    result += y[10]*hml[10]/14.026940; /*CH2 */
    result += y[11]*hml[11]/14.026940; /*CH2(S) */
    result += y[12]*hml[12]/15.034910; /*CH3 */
    result += y[13]*hml[13]/16.042880; /*CH4 */
    result += y[14]*hml[14]/28.010400; /*CO */
    result += y[15]*hml[15]/44.009800; /*CO2 */
    result += y[16]*hml[16]/29.018370; /*HCO */
    result += y[17]*hml[17]/30.026340; /*CH2O */
    result += y[18]*hml[18]/31.034310; /*CH2OH */
    result += y[19]*hml[19]/31.034310; /*CH3O */
    result += y[20]*hml[20]/32.042280; /*CH3OH */
    result += y[21]*hml[21]/25.029970; /*C2H */
    result += y[22]*hml[22]/26.037940; /*C2H2 */
    result += y[23]*hml[23]/27.045910; /*C2H3 */
    result += y[24]*hml[24]/28.053880; /*C2H4 */
    result += y[25]*hml[25]/29.061850; /*C2H5 */
    result += y[26]*hml[26]/30.069820; /*C2H6 */
    result += y[27]*hml[27]/41.029370; /*HCCO */
    result += y[28]*hml[28]/42.037340; /*CH2CO */
    result += y[29]*hml[29]/42.037340; /*HCCOH */
    result += y[30]*hml[30]/28.013400; /*N2 */
    result += y[31]*hml[31]/39.948000; /*AR */

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double *T, double *x, int * iwrk, double * rwrk, double * ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[32]; /* temporary energy array */
    double RT = 8.314e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 32; ++id) {
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
    double ums[32]; /* temporary energy array */
    double RT = 8.314e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]/2.015940; /*H2 */
    result += y[1]*ums[1]/1.007970; /*H */
    result += y[2]*ums[2]/15.999400; /*O */
    result += y[3]*ums[3]/31.998800; /*O2 */
    result += y[4]*ums[4]/17.007370; /*OH */
    result += y[5]*ums[5]/18.015340; /*H2O */
    result += y[6]*ums[6]/33.006770; /*HO2 */
    result += y[7]*ums[7]/34.014740; /*H2O2 */
    result += y[8]*ums[8]/12.011000; /*C */
    result += y[9]*ums[9]/13.018970; /*CH */
    result += y[10]*ums[10]/14.026940; /*CH2 */
    result += y[11]*ums[11]/14.026940; /*CH2(S) */
    result += y[12]*ums[12]/15.034910; /*CH3 */
    result += y[13]*ums[13]/16.042880; /*CH4 */
    result += y[14]*ums[14]/28.010400; /*CO */
    result += y[15]*ums[15]/44.009800; /*CO2 */
    result += y[16]*ums[16]/29.018370; /*HCO */
    result += y[17]*ums[17]/30.026340; /*CH2O */
    result += y[18]*ums[18]/31.034310; /*CH2OH */
    result += y[19]*ums[19]/31.034310; /*CH3O */
    result += y[20]*ums[20]/32.042280; /*CH3OH */
    result += y[21]*ums[21]/25.029970; /*C2H */
    result += y[22]*ums[22]/26.037940; /*C2H2 */
    result += y[23]*ums[23]/27.045910; /*C2H3 */
    result += y[24]*ums[24]/28.053880; /*C2H4 */
    result += y[25]*ums[25]/29.061850; /*C2H5 */
    result += y[26]*ums[26]/30.069820; /*C2H6 */
    result += y[27]*ums[27]/41.029370; /*HCCO */
    result += y[28]*ums[28]/42.037340; /*CH2CO */
    result += y[29]*ums[29]/42.037340; /*HCCOH */
    result += y[30]*ums[30]/28.013400; /*N2 */
    result += y[31]*ums[31]/39.948000; /*AR */

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
    double sor[32]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 32; ++id) {
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
    double sor[32]; /* temporary storage */
    double x[32]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/2.015940; /*H2 */
    YOW += y[1]/1.007970; /*H */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/18.015340; /*H2O */
    YOW += y[6]/33.006770; /*HO2 */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/12.011000; /*C */
    YOW += y[9]/13.018970; /*CH */
    YOW += y[10]/14.026940; /*CH2 */
    YOW += y[11]/14.026940; /*CH2(S) */
    YOW += y[12]/15.034910; /*CH3 */
    YOW += y[13]/16.042880; /*CH4 */
    YOW += y[14]/28.010400; /*CO */
    YOW += y[15]/44.009800; /*CO2 */
    YOW += y[16]/29.018370; /*HCO */
    YOW += y[17]/30.026340; /*CH2O */
    YOW += y[18]/31.034310; /*CH2OH */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/32.042280; /*CH3OH */
    YOW += y[21]/25.029970; /*C2H */
    YOW += y[22]/26.037940; /*C2H2 */
    YOW += y[23]/27.045910; /*C2H3 */
    YOW += y[24]/28.053880; /*C2H4 */
    YOW += y[25]/29.061850; /*C2H5 */
    YOW += y[26]/30.069820; /*C2H6 */
    YOW += y[27]/41.029370; /*HCCO */
    YOW += y[28]/42.037340; /*CH2CO */
    YOW += y[29]/42.037340; /*HCCOH */
    YOW += y[30]/28.013400; /*N2 */
    YOW += y[31]/39.948000; /*AR */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(12.011000*YOW); 
    x[9] = y[9]/(13.018970*YOW); 
    x[10] = y[10]/(14.026940*YOW); 
    x[11] = y[11]/(14.026940*YOW); 
    x[12] = y[12]/(15.034910*YOW); 
    x[13] = y[13]/(16.042880*YOW); 
    x[14] = y[14]/(28.010400*YOW); 
    x[15] = y[15]/(44.009800*YOW); 
    x[16] = y[16]/(29.018370*YOW); 
    x[17] = y[17]/(30.026340*YOW); 
    x[18] = y[18]/(31.034310*YOW); 
    x[19] = y[19]/(31.034310*YOW); 
    x[20] = y[20]/(32.042280*YOW); 
    x[21] = y[21]/(25.029970*YOW); 
    x[22] = y[22]/(26.037940*YOW); 
    x[23] = y[23]/(27.045910*YOW); 
    x[24] = y[24]/(28.053880*YOW); 
    x[25] = y[25]/(29.061850*YOW); 
    x[26] = y[26]/(30.069820*YOW); 
    x[27] = y[27]/(41.029370*YOW); 
    x[28] = y[28]/(42.037340*YOW); 
    x[29] = y[29]/(42.037340*YOW); 
    x[30] = y[30]/(28.013400*YOW); 
    x[31] = y[31]/(39.948000*YOW); 
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
    result += x[30]*(sor[30]-log((x[30]+1e-100))-logPratio);
    result += x[31]*(sor[31]-log((x[31]+1e-100))-logPratio);
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
    double gort[32]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 32; ++id) {
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
    double gort[32]; /* temporary storage */
    double x[32]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/2.015940; /*H2 */
    YOW += y[1]/1.007970; /*H */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/18.015340; /*H2O */
    YOW += y[6]/33.006770; /*HO2 */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/12.011000; /*C */
    YOW += y[9]/13.018970; /*CH */
    YOW += y[10]/14.026940; /*CH2 */
    YOW += y[11]/14.026940; /*CH2(S) */
    YOW += y[12]/15.034910; /*CH3 */
    YOW += y[13]/16.042880; /*CH4 */
    YOW += y[14]/28.010400; /*CO */
    YOW += y[15]/44.009800; /*CO2 */
    YOW += y[16]/29.018370; /*HCO */
    YOW += y[17]/30.026340; /*CH2O */
    YOW += y[18]/31.034310; /*CH2OH */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/32.042280; /*CH3OH */
    YOW += y[21]/25.029970; /*C2H */
    YOW += y[22]/26.037940; /*C2H2 */
    YOW += y[23]/27.045910; /*C2H3 */
    YOW += y[24]/28.053880; /*C2H4 */
    YOW += y[25]/29.061850; /*C2H5 */
    YOW += y[26]/30.069820; /*C2H6 */
    YOW += y[27]/41.029370; /*HCCO */
    YOW += y[28]/42.037340; /*CH2CO */
    YOW += y[29]/42.037340; /*HCCOH */
    YOW += y[30]/28.013400; /*N2 */
    YOW += y[31]/39.948000; /*AR */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(12.011000*YOW); 
    x[9] = y[9]/(13.018970*YOW); 
    x[10] = y[10]/(14.026940*YOW); 
    x[11] = y[11]/(14.026940*YOW); 
    x[12] = y[12]/(15.034910*YOW); 
    x[13] = y[13]/(16.042880*YOW); 
    x[14] = y[14]/(28.010400*YOW); 
    x[15] = y[15]/(44.009800*YOW); 
    x[16] = y[16]/(29.018370*YOW); 
    x[17] = y[17]/(30.026340*YOW); 
    x[18] = y[18]/(31.034310*YOW); 
    x[19] = y[19]/(31.034310*YOW); 
    x[20] = y[20]/(32.042280*YOW); 
    x[21] = y[21]/(25.029970*YOW); 
    x[22] = y[22]/(26.037940*YOW); 
    x[23] = y[23]/(27.045910*YOW); 
    x[24] = y[24]/(28.053880*YOW); 
    x[25] = y[25]/(29.061850*YOW); 
    x[26] = y[26]/(30.069820*YOW); 
    x[27] = y[27]/(41.029370*YOW); 
    x[28] = y[28]/(42.037340*YOW); 
    x[29] = y[29]/(42.037340*YOW); 
    x[30] = y[30]/(28.013400*YOW); 
    x[31] = y[31]/(39.948000*YOW); 
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
    result += x[30]*(gort[30]+log((x[30]+1e-100))+logPratio);
    result += x[31]*(gort[31]+log((x[31]+1e-100))+logPratio);
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
    double aort[32]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 32; ++id) {
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
    double aort[32]; /* temporary storage */
    double x[32]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/2.015940; /*H2 */
    YOW += y[1]/1.007970; /*H */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/18.015340; /*H2O */
    YOW += y[6]/33.006770; /*HO2 */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/12.011000; /*C */
    YOW += y[9]/13.018970; /*CH */
    YOW += y[10]/14.026940; /*CH2 */
    YOW += y[11]/14.026940; /*CH2(S) */
    YOW += y[12]/15.034910; /*CH3 */
    YOW += y[13]/16.042880; /*CH4 */
    YOW += y[14]/28.010400; /*CO */
    YOW += y[15]/44.009800; /*CO2 */
    YOW += y[16]/29.018370; /*HCO */
    YOW += y[17]/30.026340; /*CH2O */
    YOW += y[18]/31.034310; /*CH2OH */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/32.042280; /*CH3OH */
    YOW += y[21]/25.029970; /*C2H */
    YOW += y[22]/26.037940; /*C2H2 */
    YOW += y[23]/27.045910; /*C2H3 */
    YOW += y[24]/28.053880; /*C2H4 */
    YOW += y[25]/29.061850; /*C2H5 */
    YOW += y[26]/30.069820; /*C2H6 */
    YOW += y[27]/41.029370; /*HCCO */
    YOW += y[28]/42.037340; /*CH2CO */
    YOW += y[29]/42.037340; /*HCCOH */
    YOW += y[30]/28.013400; /*N2 */
    YOW += y[31]/39.948000; /*AR */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(12.011000*YOW); 
    x[9] = y[9]/(13.018970*YOW); 
    x[10] = y[10]/(14.026940*YOW); 
    x[11] = y[11]/(14.026940*YOW); 
    x[12] = y[12]/(15.034910*YOW); 
    x[13] = y[13]/(16.042880*YOW); 
    x[14] = y[14]/(28.010400*YOW); 
    x[15] = y[15]/(44.009800*YOW); 
    x[16] = y[16]/(29.018370*YOW); 
    x[17] = y[17]/(30.026340*YOW); 
    x[18] = y[18]/(31.034310*YOW); 
    x[19] = y[19]/(31.034310*YOW); 
    x[20] = y[20]/(32.042280*YOW); 
    x[21] = y[21]/(25.029970*YOW); 
    x[22] = y[22]/(26.037940*YOW); 
    x[23] = y[23]/(27.045910*YOW); 
    x[24] = y[24]/(28.053880*YOW); 
    x[25] = y[25]/(29.061850*YOW); 
    x[26] = y[26]/(30.069820*YOW); 
    x[27] = y[27]/(41.029370*YOW); 
    x[28] = y[28]/(42.037340*YOW); 
    x[29] = y[29]/(42.037340*YOW); 
    x[30] = y[30]/(28.013400*YOW); 
    x[31] = y[31]/(39.948000*YOW); 
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
    result += x[30]*(aort[30]+log((x[30]+1e-100))+logPratio);
    result += x[31]*(aort[31]+log((x[31]+1e-100))+logPratio);
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double * T, double * C, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 32; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 32; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[32]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/2.015940; /*H2 */
    YOW += y[1]/1.007970; /*H */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/18.015340; /*H2O */
    YOW += y[6]/33.006770; /*HO2 */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/12.011000; /*C */
    YOW += y[9]/13.018970; /*CH */
    YOW += y[10]/14.026940; /*CH2 */
    YOW += y[11]/14.026940; /*CH2(S) */
    YOW += y[12]/15.034910; /*CH3 */
    YOW += y[13]/16.042880; /*CH4 */
    YOW += y[14]/28.010400; /*CO */
    YOW += y[15]/44.009800; /*CO2 */
    YOW += y[16]/29.018370; /*HCO */
    YOW += y[17]/30.026340; /*CH2O */
    YOW += y[18]/31.034310; /*CH2OH */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/32.042280; /*CH3OH */
    YOW += y[21]/25.029970; /*C2H */
    YOW += y[22]/26.037940; /*C2H2 */
    YOW += y[23]/27.045910; /*C2H3 */
    YOW += y[24]/28.053880; /*C2H4 */
    YOW += y[25]/29.061850; /*C2H5 */
    YOW += y[26]/30.069820; /*C2H6 */
    YOW += y[27]/41.029370; /*HCCO */
    YOW += y[28]/42.037340; /*CH2CO */
    YOW += y[29]/42.037340; /*HCCOH */
    YOW += y[30]/28.013400; /*N2 */
    YOW += y[31]/39.948000; /*AR */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.314e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]/2.015940; 
    c[1] = PWORT * y[1]/1.007970; 
    c[2] = PWORT * y[2]/15.999400; 
    c[3] = PWORT * y[3]/31.998800; 
    c[4] = PWORT * y[4]/17.007370; 
    c[5] = PWORT * y[5]/18.015340; 
    c[6] = PWORT * y[6]/33.006770; 
    c[7] = PWORT * y[7]/34.014740; 
    c[8] = PWORT * y[8]/12.011000; 
    c[9] = PWORT * y[9]/13.018970; 
    c[10] = PWORT * y[10]/14.026940; 
    c[11] = PWORT * y[11]/14.026940; 
    c[12] = PWORT * y[12]/15.034910; 
    c[13] = PWORT * y[13]/16.042880; 
    c[14] = PWORT * y[14]/28.010400; 
    c[15] = PWORT * y[15]/44.009800; 
    c[16] = PWORT * y[16]/29.018370; 
    c[17] = PWORT * y[17]/30.026340; 
    c[18] = PWORT * y[18]/31.034310; 
    c[19] = PWORT * y[19]/31.034310; 
    c[20] = PWORT * y[20]/32.042280; 
    c[21] = PWORT * y[21]/25.029970; 
    c[22] = PWORT * y[22]/26.037940; 
    c[23] = PWORT * y[23]/27.045910; 
    c[24] = PWORT * y[24]/28.053880; 
    c[25] = PWORT * y[25]/29.061850; 
    c[26] = PWORT * y[26]/30.069820; 
    c[27] = PWORT * y[27]/41.029370; 
    c[28] = PWORT * y[28]/42.037340; 
    c[29] = PWORT * y[29]/42.037340; 
    c[30] = PWORT * y[30]/28.013400; 
    c[31] = PWORT * y[31]/39.948000; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 32; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[32]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.314e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 32; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 32; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[32]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]/2.015940; 
    c[1] = 1e6 * (*rho) * y[1]/1.007970; 
    c[2] = 1e6 * (*rho) * y[2]/15.999400; 
    c[3] = 1e6 * (*rho) * y[3]/31.998800; 
    c[4] = 1e6 * (*rho) * y[4]/17.007370; 
    c[5] = 1e6 * (*rho) * y[5]/18.015340; 
    c[6] = 1e6 * (*rho) * y[6]/33.006770; 
    c[7] = 1e6 * (*rho) * y[7]/34.014740; 
    c[8] = 1e6 * (*rho) * y[8]/12.011000; 
    c[9] = 1e6 * (*rho) * y[9]/13.018970; 
    c[10] = 1e6 * (*rho) * y[10]/14.026940; 
    c[11] = 1e6 * (*rho) * y[11]/14.026940; 
    c[12] = 1e6 * (*rho) * y[12]/15.034910; 
    c[13] = 1e6 * (*rho) * y[13]/16.042880; 
    c[14] = 1e6 * (*rho) * y[14]/28.010400; 
    c[15] = 1e6 * (*rho) * y[15]/44.009800; 
    c[16] = 1e6 * (*rho) * y[16]/29.018370; 
    c[17] = 1e6 * (*rho) * y[17]/30.026340; 
    c[18] = 1e6 * (*rho) * y[18]/31.034310; 
    c[19] = 1e6 * (*rho) * y[19]/31.034310; 
    c[20] = 1e6 * (*rho) * y[20]/32.042280; 
    c[21] = 1e6 * (*rho) * y[21]/25.029970; 
    c[22] = 1e6 * (*rho) * y[22]/26.037940; 
    c[23] = 1e6 * (*rho) * y[23]/27.045910; 
    c[24] = 1e6 * (*rho) * y[24]/28.053880; 
    c[25] = 1e6 * (*rho) * y[25]/29.061850; 
    c[26] = 1e6 * (*rho) * y[26]/30.069820; 
    c[27] = 1e6 * (*rho) * y[27]/41.029370; 
    c[28] = 1e6 * (*rho) * y[28]/42.037340; 
    c[29] = 1e6 * (*rho) * y[29]/42.037340; 
    c[30] = 1e6 * (*rho) * y[30]/28.013400; 
    c[31] = 1e6 * (*rho) * y[31]/39.948000; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 32; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[32]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*12.011000; /*C */
    XW += x[9]*13.018970; /*CH */
    XW += x[10]*14.026940; /*CH2 */
    XW += x[11]*14.026940; /*CH2(S) */
    XW += x[12]*15.034910; /*CH3 */
    XW += x[13]*16.042880; /*CH4 */
    XW += x[14]*28.010400; /*CO */
    XW += x[15]*44.009800; /*CO2 */
    XW += x[16]*29.018370; /*HCO */
    XW += x[17]*30.026340; /*CH2O */
    XW += x[18]*31.034310; /*CH2OH */
    XW += x[19]*31.034310; /*CH3O */
    XW += x[20]*32.042280; /*CH3OH */
    XW += x[21]*25.029970; /*C2H */
    XW += x[22]*26.037940; /*C2H2 */
    XW += x[23]*27.045910; /*C2H3 */
    XW += x[24]*28.053880; /*C2H4 */
    XW += x[25]*29.061850; /*C2H5 */
    XW += x[26]*30.069820; /*C2H6 */
    XW += x[27]*41.029370; /*HCCO */
    XW += x[28]*42.037340; /*CH2CO */
    XW += x[29]*42.037340; /*HCCOH */
    XW += x[30]*28.013400; /*N2 */
    XW += x[31]*39.948000; /*AR */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 32; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 32; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double * T, double * C, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 32; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 32; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 177; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[32]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/2.015940; /*H2 */
    YOW += y[1]/1.007970; /*H */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/18.015340; /*H2O */
    YOW += y[6]/33.006770; /*HO2 */
    YOW += y[7]/34.014740; /*H2O2 */
    YOW += y[8]/12.011000; /*C */
    YOW += y[9]/13.018970; /*CH */
    YOW += y[10]/14.026940; /*CH2 */
    YOW += y[11]/14.026940; /*CH2(S) */
    YOW += y[12]/15.034910; /*CH3 */
    YOW += y[13]/16.042880; /*CH4 */
    YOW += y[14]/28.010400; /*CO */
    YOW += y[15]/44.009800; /*CO2 */
    YOW += y[16]/29.018370; /*HCO */
    YOW += y[17]/30.026340; /*CH2O */
    YOW += y[18]/31.034310; /*CH2OH */
    YOW += y[19]/31.034310; /*CH3O */
    YOW += y[20]/32.042280; /*CH3OH */
    YOW += y[21]/25.029970; /*C2H */
    YOW += y[22]/26.037940; /*C2H2 */
    YOW += y[23]/27.045910; /*C2H3 */
    YOW += y[24]/28.053880; /*C2H4 */
    YOW += y[25]/29.061850; /*C2H5 */
    YOW += y[26]/30.069820; /*C2H6 */
    YOW += y[27]/41.029370; /*HCCO */
    YOW += y[28]/42.037340; /*CH2CO */
    YOW += y[29]/42.037340; /*HCCOH */
    YOW += y[30]/28.013400; /*N2 */
    YOW += y[31]/39.948000; /*AR */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.314e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]/2.015940; 
    c[1] = PWORT * y[1]/1.007970; 
    c[2] = PWORT * y[2]/15.999400; 
    c[3] = PWORT * y[3]/31.998800; 
    c[4] = PWORT * y[4]/17.007370; 
    c[5] = PWORT * y[5]/18.015340; 
    c[6] = PWORT * y[6]/33.006770; 
    c[7] = PWORT * y[7]/34.014740; 
    c[8] = PWORT * y[8]/12.011000; 
    c[9] = PWORT * y[9]/13.018970; 
    c[10] = PWORT * y[10]/14.026940; 
    c[11] = PWORT * y[11]/14.026940; 
    c[12] = PWORT * y[12]/15.034910; 
    c[13] = PWORT * y[13]/16.042880; 
    c[14] = PWORT * y[14]/28.010400; 
    c[15] = PWORT * y[15]/44.009800; 
    c[16] = PWORT * y[16]/29.018370; 
    c[17] = PWORT * y[17]/30.026340; 
    c[18] = PWORT * y[18]/31.034310; 
    c[19] = PWORT * y[19]/31.034310; 
    c[20] = PWORT * y[20]/32.042280; 
    c[21] = PWORT * y[21]/25.029970; 
    c[22] = PWORT * y[22]/26.037940; 
    c[23] = PWORT * y[23]/27.045910; 
    c[24] = PWORT * y[24]/28.053880; 
    c[25] = PWORT * y[25]/29.061850; 
    c[26] = PWORT * y[26]/30.069820; 
    c[27] = PWORT * y[27]/41.029370; 
    c[28] = PWORT * y[28]/42.037340; 
    c[29] = PWORT * y[29]/42.037340; 
    c[30] = PWORT * y[30]/28.013400; 
    c[31] = PWORT * y[31]/39.948000; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 177; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[32]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.314e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 32; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 177; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[32]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]/2.015940; 
    c[1] = 1e6 * (*rho) * y[1]/1.007970; 
    c[2] = 1e6 * (*rho) * y[2]/15.999400; 
    c[3] = 1e6 * (*rho) * y[3]/31.998800; 
    c[4] = 1e6 * (*rho) * y[4]/17.007370; 
    c[5] = 1e6 * (*rho) * y[5]/18.015340; 
    c[6] = 1e6 * (*rho) * y[6]/33.006770; 
    c[7] = 1e6 * (*rho) * y[7]/34.014740; 
    c[8] = 1e6 * (*rho) * y[8]/12.011000; 
    c[9] = 1e6 * (*rho) * y[9]/13.018970; 
    c[10] = 1e6 * (*rho) * y[10]/14.026940; 
    c[11] = 1e6 * (*rho) * y[11]/14.026940; 
    c[12] = 1e6 * (*rho) * y[12]/15.034910; 
    c[13] = 1e6 * (*rho) * y[13]/16.042880; 
    c[14] = 1e6 * (*rho) * y[14]/28.010400; 
    c[15] = 1e6 * (*rho) * y[15]/44.009800; 
    c[16] = 1e6 * (*rho) * y[16]/29.018370; 
    c[17] = 1e6 * (*rho) * y[17]/30.026340; 
    c[18] = 1e6 * (*rho) * y[18]/31.034310; 
    c[19] = 1e6 * (*rho) * y[19]/31.034310; 
    c[20] = 1e6 * (*rho) * y[20]/32.042280; 
    c[21] = 1e6 * (*rho) * y[21]/25.029970; 
    c[22] = 1e6 * (*rho) * y[22]/26.037940; 
    c[23] = 1e6 * (*rho) * y[23]/27.045910; 
    c[24] = 1e6 * (*rho) * y[24]/28.053880; 
    c[25] = 1e6 * (*rho) * y[25]/29.061850; 
    c[26] = 1e6 * (*rho) * y[26]/30.069820; 
    c[27] = 1e6 * (*rho) * y[27]/41.029370; 
    c[28] = 1e6 * (*rho) * y[28]/42.037340; 
    c[29] = 1e6 * (*rho) * y[29]/42.037340; 
    c[30] = 1e6 * (*rho) * y[30]/28.013400; 
    c[31] = 1e6 * (*rho) * y[31]/39.948000; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 177; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[32]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*12.011000; /*C */
    XW += x[9]*13.018970; /*CH */
    XW += x[10]*14.026940; /*CH2 */
    XW += x[11]*14.026940; /*CH2(S) */
    XW += x[12]*15.034910; /*CH3 */
    XW += x[13]*16.042880; /*CH4 */
    XW += x[14]*28.010400; /*CO */
    XW += x[15]*44.009800; /*CO2 */
    XW += x[16]*29.018370; /*HCO */
    XW += x[17]*30.026340; /*CH2O */
    XW += x[18]*31.034310; /*CH2OH */
    XW += x[19]*31.034310; /*CH3O */
    XW += x[20]*32.042280; /*CH3OH */
    XW += x[21]*25.029970; /*C2H */
    XW += x[22]*26.037940; /*C2H2 */
    XW += x[23]*27.045910; /*C2H3 */
    XW += x[24]*28.053880; /*C2H4 */
    XW += x[25]*29.061850; /*C2H5 */
    XW += x[26]*30.069820; /*C2H6 */
    XW += x[27]*41.029370; /*HCCO */
    XW += x[28]*42.037340; /*CH2CO */
    XW += x[29]*42.037340; /*HCCOH */
    XW += x[30]*28.013400; /*N2 */
    XW += x[31]*39.948000; /*AR */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 32; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 177; ++id) {
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
    for (id = 0; id < 32 * 177; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: 2 O + M <=> O2 + M */
    nuki[ 2 * kd + 0 ] += -2 ;
    nuki[ 3 * kd + 0 ] += +1 ;

    /*reaction 2: O + H + M <=> OH + M */
    nuki[ 2 * kd + 1 ] += -1 ;
    nuki[ 1 * kd + 1 ] += -1 ;
    nuki[ 4 * kd + 1 ] += +1 ;

    /*reaction 3: O + H2 <=> H + OH */
    nuki[ 2 * kd + 2 ] += -1 ;
    nuki[ 0 * kd + 2 ] += -1 ;
    nuki[ 1 * kd + 2 ] += +1 ;
    nuki[ 4 * kd + 2 ] += +1 ;

    /*reaction 4: O + HO2 <=> OH + O2 */
    nuki[ 2 * kd + 3 ] += -1 ;
    nuki[ 6 * kd + 3 ] += -1 ;
    nuki[ 4 * kd + 3 ] += +1 ;
    nuki[ 3 * kd + 3 ] += +1 ;

    /*reaction 5: O + H2O2 <=> OH + HO2 */
    nuki[ 2 * kd + 4 ] += -1 ;
    nuki[ 7 * kd + 4 ] += -1 ;
    nuki[ 4 * kd + 4 ] += +1 ;
    nuki[ 6 * kd + 4 ] += +1 ;

    /*reaction 6: O + CH <=> H + CO */
    nuki[ 2 * kd + 5 ] += -1 ;
    nuki[ 9 * kd + 5 ] += -1 ;
    nuki[ 1 * kd + 5 ] += +1 ;
    nuki[ 14 * kd + 5 ] += +1 ;

    /*reaction 7: O + CH2 <=> H + HCO */
    nuki[ 2 * kd + 6 ] += -1 ;
    nuki[ 10 * kd + 6 ] += -1 ;
    nuki[ 1 * kd + 6 ] += +1 ;
    nuki[ 16 * kd + 6 ] += +1 ;

    /*reaction 8: O + CH2(S) <=> H2 + CO */
    nuki[ 2 * kd + 7 ] += -1 ;
    nuki[ 11 * kd + 7 ] += -1 ;
    nuki[ 0 * kd + 7 ] += +1 ;
    nuki[ 14 * kd + 7 ] += +1 ;

    /*reaction 9: O + CH2(S) <=> H + HCO */
    nuki[ 2 * kd + 8 ] += -1 ;
    nuki[ 11 * kd + 8 ] += -1 ;
    nuki[ 1 * kd + 8 ] += +1 ;
    nuki[ 16 * kd + 8 ] += +1 ;

    /*reaction 10: O + CH3 <=> H + CH2O */
    nuki[ 2 * kd + 9 ] += -1 ;
    nuki[ 12 * kd + 9 ] += -1 ;
    nuki[ 1 * kd + 9 ] += +1 ;
    nuki[ 17 * kd + 9 ] += +1 ;

    /*reaction 11: O + CH4 <=> OH + CH3 */
    nuki[ 2 * kd + 10 ] += -1 ;
    nuki[ 13 * kd + 10 ] += -1 ;
    nuki[ 4 * kd + 10 ] += +1 ;
    nuki[ 12 * kd + 10 ] += +1 ;

    /*reaction 12: O + CO + M <=> CO2 + M */
    nuki[ 2 * kd + 11 ] += -1 ;
    nuki[ 14 * kd + 11 ] += -1 ;
    nuki[ 15 * kd + 11 ] += +1 ;

    /*reaction 13: O + HCO <=> OH + CO */
    nuki[ 2 * kd + 12 ] += -1 ;
    nuki[ 16 * kd + 12 ] += -1 ;
    nuki[ 4 * kd + 12 ] += +1 ;
    nuki[ 14 * kd + 12 ] += +1 ;

    /*reaction 14: O + HCO <=> H + CO2 */
    nuki[ 2 * kd + 13 ] += -1 ;
    nuki[ 16 * kd + 13 ] += -1 ;
    nuki[ 1 * kd + 13 ] += +1 ;
    nuki[ 15 * kd + 13 ] += +1 ;

    /*reaction 15: O + CH2O <=> OH + HCO */
    nuki[ 2 * kd + 14 ] += -1 ;
    nuki[ 17 * kd + 14 ] += -1 ;
    nuki[ 4 * kd + 14 ] += +1 ;
    nuki[ 16 * kd + 14 ] += +1 ;

    /*reaction 16: O + CH2OH <=> OH + CH2O */
    nuki[ 2 * kd + 15 ] += -1 ;
    nuki[ 18 * kd + 15 ] += -1 ;
    nuki[ 4 * kd + 15 ] += +1 ;
    nuki[ 17 * kd + 15 ] += +1 ;

    /*reaction 17: O + CH3O <=> OH + CH2O */
    nuki[ 2 * kd + 16 ] += -1 ;
    nuki[ 19 * kd + 16 ] += -1 ;
    nuki[ 4 * kd + 16 ] += +1 ;
    nuki[ 17 * kd + 16 ] += +1 ;

    /*reaction 18: O + CH3OH <=> OH + CH2OH */
    nuki[ 2 * kd + 17 ] += -1 ;
    nuki[ 20 * kd + 17 ] += -1 ;
    nuki[ 4 * kd + 17 ] += +1 ;
    nuki[ 18 * kd + 17 ] += +1 ;

    /*reaction 19: O + CH3OH <=> OH + CH3O */
    nuki[ 2 * kd + 18 ] += -1 ;
    nuki[ 20 * kd + 18 ] += -1 ;
    nuki[ 4 * kd + 18 ] += +1 ;
    nuki[ 19 * kd + 18 ] += +1 ;

    /*reaction 20: O + C2H <=> CH + CO */
    nuki[ 2 * kd + 19 ] += -1 ;
    nuki[ 21 * kd + 19 ] += -1 ;
    nuki[ 9 * kd + 19 ] += +1 ;
    nuki[ 14 * kd + 19 ] += +1 ;

    /*reaction 21: O + C2H2 <=> H + HCCO */
    nuki[ 2 * kd + 20 ] += -1 ;
    nuki[ 22 * kd + 20 ] += -1 ;
    nuki[ 1 * kd + 20 ] += +1 ;
    nuki[ 27 * kd + 20 ] += +1 ;

    /*reaction 22: O + C2H2 <=> OH + C2H */
    nuki[ 2 * kd + 21 ] += -1 ;
    nuki[ 22 * kd + 21 ] += -1 ;
    nuki[ 4 * kd + 21 ] += +1 ;
    nuki[ 21 * kd + 21 ] += +1 ;

    /*reaction 23: O + C2H2 <=> CO + CH2 */
    nuki[ 2 * kd + 22 ] += -1 ;
    nuki[ 22 * kd + 22 ] += -1 ;
    nuki[ 14 * kd + 22 ] += +1 ;
    nuki[ 10 * kd + 22 ] += +1 ;

    /*reaction 24: O + C2H3 <=> H + CH2CO */
    nuki[ 2 * kd + 23 ] += -1 ;
    nuki[ 23 * kd + 23 ] += -1 ;
    nuki[ 1 * kd + 23 ] += +1 ;
    nuki[ 28 * kd + 23 ] += +1 ;

    /*reaction 25: O + C2H4 <=> CH3 + HCO */
    nuki[ 2 * kd + 24 ] += -1 ;
    nuki[ 24 * kd + 24 ] += -1 ;
    nuki[ 12 * kd + 24 ] += +1 ;
    nuki[ 16 * kd + 24 ] += +1 ;

    /*reaction 26: O + C2H5 <=> CH3 + CH2O */
    nuki[ 2 * kd + 25 ] += -1 ;
    nuki[ 25 * kd + 25 ] += -1 ;
    nuki[ 12 * kd + 25 ] += +1 ;
    nuki[ 17 * kd + 25 ] += +1 ;

    /*reaction 27: O + C2H6 <=> OH + C2H5 */
    nuki[ 2 * kd + 26 ] += -1 ;
    nuki[ 26 * kd + 26 ] += -1 ;
    nuki[ 4 * kd + 26 ] += +1 ;
    nuki[ 25 * kd + 26 ] += +1 ;

    /*reaction 28: O + HCCO <=> H + 2 CO */
    nuki[ 2 * kd + 27 ] += -1 ;
    nuki[ 27 * kd + 27 ] += -1 ;
    nuki[ 1 * kd + 27 ] += +1 ;
    nuki[ 14 * kd + 27 ] += +2 ;

    /*reaction 29: O + CH2CO <=> OH + HCCO */
    nuki[ 2 * kd + 28 ] += -1 ;
    nuki[ 28 * kd + 28 ] += -1 ;
    nuki[ 4 * kd + 28 ] += +1 ;
    nuki[ 27 * kd + 28 ] += +1 ;

    /*reaction 30: O + CH2CO <=> CH2 + CO2 */
    nuki[ 2 * kd + 29 ] += -1 ;
    nuki[ 28 * kd + 29 ] += -1 ;
    nuki[ 10 * kd + 29 ] += +1 ;
    nuki[ 15 * kd + 29 ] += +1 ;

    /*reaction 31: O2 + CO <=> O + CO2 */
    nuki[ 3 * kd + 30 ] += -1 ;
    nuki[ 14 * kd + 30 ] += -1 ;
    nuki[ 2 * kd + 30 ] += +1 ;
    nuki[ 15 * kd + 30 ] += +1 ;

    /*reaction 32: O2 + CH2O <=> HO2 + HCO */
    nuki[ 3 * kd + 31 ] += -1 ;
    nuki[ 17 * kd + 31 ] += -1 ;
    nuki[ 6 * kd + 31 ] += +1 ;
    nuki[ 16 * kd + 31 ] += +1 ;

    /*reaction 33: H + O2 + M <=> HO2 + M */
    nuki[ 1 * kd + 32 ] += -1 ;
    nuki[ 3 * kd + 32 ] += -1 ;
    nuki[ 6 * kd + 32 ] += +1 ;

    /*reaction 34: H + 2 O2 <=> HO2 + O2 */
    nuki[ 1 * kd + 33 ] += -1 ;
    nuki[ 3 * kd + 33 ] += -2 ;
    nuki[ 6 * kd + 33 ] += +1 ;
    nuki[ 3 * kd + 33 ] += +1 ;

    /*reaction 35: H + O2 + H2O <=> HO2 + H2O */
    nuki[ 1 * kd + 34 ] += -1 ;
    nuki[ 3 * kd + 34 ] += -1 ;
    nuki[ 5 * kd + 34 ] += -1 ;
    nuki[ 6 * kd + 34 ] += +1 ;
    nuki[ 5 * kd + 34 ] += +1 ;

    /*reaction 36: H + O2 + N2 <=> HO2 + N2 */
    nuki[ 1 * kd + 35 ] += -1 ;
    nuki[ 3 * kd + 35 ] += -1 ;
    nuki[ 30 * kd + 35 ] += -1 ;
    nuki[ 6 * kd + 35 ] += +1 ;
    nuki[ 30 * kd + 35 ] += +1 ;

    /*reaction 37: H + O2 + AR <=> HO2 + AR */
    nuki[ 1 * kd + 36 ] += -1 ;
    nuki[ 3 * kd + 36 ] += -1 ;
    nuki[ 31 * kd + 36 ] += -1 ;
    nuki[ 6 * kd + 36 ] += +1 ;
    nuki[ 31 * kd + 36 ] += +1 ;

    /*reaction 38: H + O2 <=> O + OH */
    nuki[ 1 * kd + 37 ] += -1 ;
    nuki[ 3 * kd + 37 ] += -1 ;
    nuki[ 2 * kd + 37 ] += +1 ;
    nuki[ 4 * kd + 37 ] += +1 ;

    /*reaction 39: 2 H + M <=> H2 + M */
    nuki[ 1 * kd + 38 ] += -2 ;
    nuki[ 0 * kd + 38 ] += +1 ;

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    nuki[ 1 * kd + 39 ] += -2 ;
    nuki[ 0 * kd + 39 ] += -1 ;
    nuki[ 0 * kd + 39 ] += +2 ;

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    nuki[ 1 * kd + 40 ] += -2 ;
    nuki[ 5 * kd + 40 ] += -1 ;
    nuki[ 0 * kd + 40 ] += +1 ;
    nuki[ 5 * kd + 40 ] += +1 ;

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    nuki[ 1 * kd + 41 ] += -2 ;
    nuki[ 15 * kd + 41 ] += -1 ;
    nuki[ 0 * kd + 41 ] += +1 ;
    nuki[ 15 * kd + 41 ] += +1 ;

    /*reaction 43: H + OH + M <=> H2O + M */
    nuki[ 1 * kd + 42 ] += -1 ;
    nuki[ 4 * kd + 42 ] += -1 ;
    nuki[ 5 * kd + 42 ] += +1 ;

    /*reaction 44: H + HO2 <=> O + H2O */
    nuki[ 1 * kd + 43 ] += -1 ;
    nuki[ 6 * kd + 43 ] += -1 ;
    nuki[ 2 * kd + 43 ] += +1 ;
    nuki[ 5 * kd + 43 ] += +1 ;

    /*reaction 45: H + HO2 <=> O2 + H2 */
    nuki[ 1 * kd + 44 ] += -1 ;
    nuki[ 6 * kd + 44 ] += -1 ;
    nuki[ 3 * kd + 44 ] += +1 ;
    nuki[ 0 * kd + 44 ] += +1 ;

    /*reaction 46: H + HO2 <=> 2 OH */
    nuki[ 1 * kd + 45 ] += -1 ;
    nuki[ 6 * kd + 45 ] += -1 ;
    nuki[ 4 * kd + 45 ] += +2 ;

    /*reaction 47: H + H2O2 <=> HO2 + H2 */
    nuki[ 1 * kd + 46 ] += -1 ;
    nuki[ 7 * kd + 46 ] += -1 ;
    nuki[ 6 * kd + 46 ] += +1 ;
    nuki[ 0 * kd + 46 ] += +1 ;

    /*reaction 48: H + H2O2 <=> OH + H2O */
    nuki[ 1 * kd + 47 ] += -1 ;
    nuki[ 7 * kd + 47 ] += -1 ;
    nuki[ 4 * kd + 47 ] += +1 ;
    nuki[ 5 * kd + 47 ] += +1 ;

    /*reaction 49: H + CH <=> C + H2 */
    nuki[ 1 * kd + 48 ] += -1 ;
    nuki[ 9 * kd + 48 ] += -1 ;
    nuki[ 8 * kd + 48 ] += +1 ;
    nuki[ 0 * kd + 48 ] += +1 ;

    /*reaction 50: H + CH2 (+M) <=> CH3 (+M) */
    nuki[ 1 * kd + 49 ] += -1 ;
    nuki[ 10 * kd + 49 ] += -1 ;
    nuki[ 12 * kd + 49 ] += +1 ;

    /*reaction 51: H + CH2(S) <=> CH + H2 */
    nuki[ 1 * kd + 50 ] += -1 ;
    nuki[ 11 * kd + 50 ] += -1 ;
    nuki[ 9 * kd + 50 ] += +1 ;
    nuki[ 0 * kd + 50 ] += +1 ;

    /*reaction 52: H + CH3 (+M) <=> CH4 (+M) */
    nuki[ 1 * kd + 51 ] += -1 ;
    nuki[ 12 * kd + 51 ] += -1 ;
    nuki[ 13 * kd + 51 ] += +1 ;

    /*reaction 53: H + CH4 <=> CH3 + H2 */
    nuki[ 1 * kd + 52 ] += -1 ;
    nuki[ 13 * kd + 52 ] += -1 ;
    nuki[ 12 * kd + 52 ] += +1 ;
    nuki[ 0 * kd + 52 ] += +1 ;

    /*reaction 54: H + HCO (+M) <=> CH2O (+M) */
    nuki[ 1 * kd + 53 ] += -1 ;
    nuki[ 16 * kd + 53 ] += -1 ;
    nuki[ 17 * kd + 53 ] += +1 ;

    /*reaction 55: H + HCO <=> H2 + CO */
    nuki[ 1 * kd + 54 ] += -1 ;
    nuki[ 16 * kd + 54 ] += -1 ;
    nuki[ 0 * kd + 54 ] += +1 ;
    nuki[ 14 * kd + 54 ] += +1 ;

    /*reaction 56: H + CH2O (+M) <=> CH2OH (+M) */
    nuki[ 1 * kd + 55 ] += -1 ;
    nuki[ 17 * kd + 55 ] += -1 ;
    nuki[ 18 * kd + 55 ] += +1 ;

    /*reaction 57: H + CH2O (+M) <=> CH3O (+M) */
    nuki[ 1 * kd + 56 ] += -1 ;
    nuki[ 17 * kd + 56 ] += -1 ;
    nuki[ 19 * kd + 56 ] += +1 ;

    /*reaction 58: H + CH2O <=> HCO + H2 */
    nuki[ 1 * kd + 57 ] += -1 ;
    nuki[ 17 * kd + 57 ] += -1 ;
    nuki[ 16 * kd + 57 ] += +1 ;
    nuki[ 0 * kd + 57 ] += +1 ;

    /*reaction 59: H + CH2OH (+M) <=> CH3OH (+M) */
    nuki[ 1 * kd + 58 ] += -1 ;
    nuki[ 18 * kd + 58 ] += -1 ;
    nuki[ 20 * kd + 58 ] += +1 ;

    /*reaction 60: H + CH2OH <=> H2 + CH2O */
    nuki[ 1 * kd + 59 ] += -1 ;
    nuki[ 18 * kd + 59 ] += -1 ;
    nuki[ 0 * kd + 59 ] += +1 ;
    nuki[ 17 * kd + 59 ] += +1 ;

    /*reaction 61: H + CH2OH <=> OH + CH3 */
    nuki[ 1 * kd + 60 ] += -1 ;
    nuki[ 18 * kd + 60 ] += -1 ;
    nuki[ 4 * kd + 60 ] += +1 ;
    nuki[ 12 * kd + 60 ] += +1 ;

    /*reaction 62: H + CH2OH <=> CH2(S) + H2O */
    nuki[ 1 * kd + 61 ] += -1 ;
    nuki[ 18 * kd + 61 ] += -1 ;
    nuki[ 11 * kd + 61 ] += +1 ;
    nuki[ 5 * kd + 61 ] += +1 ;

    /*reaction 63: H + CH3O (+M) <=> CH3OH (+M) */
    nuki[ 1 * kd + 62 ] += -1 ;
    nuki[ 19 * kd + 62 ] += -1 ;
    nuki[ 20 * kd + 62 ] += +1 ;

    /*reaction 64: H + CH3O <=> H + CH2OH */
    nuki[ 1 * kd + 63 ] += -1 ;
    nuki[ 19 * kd + 63 ] += -1 ;
    nuki[ 1 * kd + 63 ] += +1 ;
    nuki[ 18 * kd + 63 ] += +1 ;

    /*reaction 65: H + CH3O <=> H2 + CH2O */
    nuki[ 1 * kd + 64 ] += -1 ;
    nuki[ 19 * kd + 64 ] += -1 ;
    nuki[ 0 * kd + 64 ] += +1 ;
    nuki[ 17 * kd + 64 ] += +1 ;

    /*reaction 66: H + CH3O <=> OH + CH3 */
    nuki[ 1 * kd + 65 ] += -1 ;
    nuki[ 19 * kd + 65 ] += -1 ;
    nuki[ 4 * kd + 65 ] += +1 ;
    nuki[ 12 * kd + 65 ] += +1 ;

    /*reaction 67: H + CH3O <=> CH2(S) + H2O */
    nuki[ 1 * kd + 66 ] += -1 ;
    nuki[ 19 * kd + 66 ] += -1 ;
    nuki[ 11 * kd + 66 ] += +1 ;
    nuki[ 5 * kd + 66 ] += +1 ;

    /*reaction 68: H + CH3OH <=> CH2OH + H2 */
    nuki[ 1 * kd + 67 ] += -1 ;
    nuki[ 20 * kd + 67 ] += -1 ;
    nuki[ 18 * kd + 67 ] += +1 ;
    nuki[ 0 * kd + 67 ] += +1 ;

    /*reaction 69: H + CH3OH <=> CH3O + H2 */
    nuki[ 1 * kd + 68 ] += -1 ;
    nuki[ 20 * kd + 68 ] += -1 ;
    nuki[ 19 * kd + 68 ] += +1 ;
    nuki[ 0 * kd + 68 ] += +1 ;

    /*reaction 70: H + C2H (+M) <=> C2H2 (+M) */
    nuki[ 1 * kd + 69 ] += -1 ;
    nuki[ 21 * kd + 69 ] += -1 ;
    nuki[ 22 * kd + 69 ] += +1 ;

    /*reaction 71: H + C2H2 (+M) <=> C2H3 (+M) */
    nuki[ 1 * kd + 70 ] += -1 ;
    nuki[ 22 * kd + 70 ] += -1 ;
    nuki[ 23 * kd + 70 ] += +1 ;

    /*reaction 72: H + C2H3 (+M) <=> C2H4 (+M) */
    nuki[ 1 * kd + 71 ] += -1 ;
    nuki[ 23 * kd + 71 ] += -1 ;
    nuki[ 24 * kd + 71 ] += +1 ;

    /*reaction 73: H + C2H3 <=> H2 + C2H2 */
    nuki[ 1 * kd + 72 ] += -1 ;
    nuki[ 23 * kd + 72 ] += -1 ;
    nuki[ 0 * kd + 72 ] += +1 ;
    nuki[ 22 * kd + 72 ] += +1 ;

    /*reaction 74: H + C2H4 (+M) <=> C2H5 (+M) */
    nuki[ 1 * kd + 73 ] += -1 ;
    nuki[ 24 * kd + 73 ] += -1 ;
    nuki[ 25 * kd + 73 ] += +1 ;

    /*reaction 75: H + C2H4 <=> C2H3 + H2 */
    nuki[ 1 * kd + 74 ] += -1 ;
    nuki[ 24 * kd + 74 ] += -1 ;
    nuki[ 23 * kd + 74 ] += +1 ;
    nuki[ 0 * kd + 74 ] += +1 ;

    /*reaction 76: H + C2H5 (+M) <=> C2H6 (+M) */
    nuki[ 1 * kd + 75 ] += -1 ;
    nuki[ 25 * kd + 75 ] += -1 ;
    nuki[ 26 * kd + 75 ] += +1 ;

    /*reaction 77: H + C2H5 <=> H2 + C2H4 */
    nuki[ 1 * kd + 76 ] += -1 ;
    nuki[ 25 * kd + 76 ] += -1 ;
    nuki[ 0 * kd + 76 ] += +1 ;
    nuki[ 24 * kd + 76 ] += +1 ;

    /*reaction 78: H + C2H6 <=> C2H5 + H2 */
    nuki[ 1 * kd + 77 ] += -1 ;
    nuki[ 26 * kd + 77 ] += -1 ;
    nuki[ 25 * kd + 77 ] += +1 ;
    nuki[ 0 * kd + 77 ] += +1 ;

    /*reaction 79: H + HCCO <=> CH2(S) + CO */
    nuki[ 1 * kd + 78 ] += -1 ;
    nuki[ 27 * kd + 78 ] += -1 ;
    nuki[ 11 * kd + 78 ] += +1 ;
    nuki[ 14 * kd + 78 ] += +1 ;

    /*reaction 80: H + CH2CO <=> HCCO + H2 */
    nuki[ 1 * kd + 79 ] += -1 ;
    nuki[ 28 * kd + 79 ] += -1 ;
    nuki[ 27 * kd + 79 ] += +1 ;
    nuki[ 0 * kd + 79 ] += +1 ;

    /*reaction 81: H + CH2CO <=> CH3 + CO */
    nuki[ 1 * kd + 80 ] += -1 ;
    nuki[ 28 * kd + 80 ] += -1 ;
    nuki[ 12 * kd + 80 ] += +1 ;
    nuki[ 14 * kd + 80 ] += +1 ;

    /*reaction 82: H + HCCOH <=> H + CH2CO */
    nuki[ 1 * kd + 81 ] += -1 ;
    nuki[ 29 * kd + 81 ] += -1 ;
    nuki[ 1 * kd + 81 ] += +1 ;
    nuki[ 28 * kd + 81 ] += +1 ;

    /*reaction 83: H2 + CO (+M) <=> CH2O (+M) */
    nuki[ 0 * kd + 82 ] += -1 ;
    nuki[ 14 * kd + 82 ] += -1 ;
    nuki[ 17 * kd + 82 ] += +1 ;

    /*reaction 84: OH + H2 <=> H + H2O */
    nuki[ 4 * kd + 83 ] += -1 ;
    nuki[ 0 * kd + 83 ] += -1 ;
    nuki[ 1 * kd + 83 ] += +1 ;
    nuki[ 5 * kd + 83 ] += +1 ;

    /*reaction 85: 2 OH (+M) <=> H2O2 (+M) */
    nuki[ 4 * kd + 84 ] += -2 ;
    nuki[ 7 * kd + 84 ] += +1 ;

    /*reaction 86: 2 OH <=> O + H2O */
    nuki[ 4 * kd + 85 ] += -2 ;
    nuki[ 2 * kd + 85 ] += +1 ;
    nuki[ 5 * kd + 85 ] += +1 ;

    /*reaction 87: OH + HO2 <=> O2 + H2O */
    nuki[ 4 * kd + 86 ] += -1 ;
    nuki[ 6 * kd + 86 ] += -1 ;
    nuki[ 3 * kd + 86 ] += +1 ;
    nuki[ 5 * kd + 86 ] += +1 ;

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    nuki[ 4 * kd + 87 ] += -1 ;
    nuki[ 7 * kd + 87 ] += -1 ;
    nuki[ 6 * kd + 87 ] += +1 ;
    nuki[ 5 * kd + 87 ] += +1 ;

    /*reaction 89: OH + H2O2 <=> HO2 + H2O */
    nuki[ 4 * kd + 88 ] += -1 ;
    nuki[ 7 * kd + 88 ] += -1 ;
    nuki[ 6 * kd + 88 ] += +1 ;
    nuki[ 5 * kd + 88 ] += +1 ;

    /*reaction 90: OH + C <=> H + CO */
    nuki[ 4 * kd + 89 ] += -1 ;
    nuki[ 8 * kd + 89 ] += -1 ;
    nuki[ 1 * kd + 89 ] += +1 ;
    nuki[ 14 * kd + 89 ] += +1 ;

    /*reaction 91: OH + CH <=> H + HCO */
    nuki[ 4 * kd + 90 ] += -1 ;
    nuki[ 9 * kd + 90 ] += -1 ;
    nuki[ 1 * kd + 90 ] += +1 ;
    nuki[ 16 * kd + 90 ] += +1 ;

    /*reaction 92: OH + CH2 <=> H + CH2O */
    nuki[ 4 * kd + 91 ] += -1 ;
    nuki[ 10 * kd + 91 ] += -1 ;
    nuki[ 1 * kd + 91 ] += +1 ;
    nuki[ 17 * kd + 91 ] += +1 ;

    /*reaction 93: OH + CH2 <=> CH + H2O */
    nuki[ 4 * kd + 92 ] += -1 ;
    nuki[ 10 * kd + 92 ] += -1 ;
    nuki[ 9 * kd + 92 ] += +1 ;
    nuki[ 5 * kd + 92 ] += +1 ;

    /*reaction 94: OH + CH2(S) <=> H + CH2O */
    nuki[ 4 * kd + 93 ] += -1 ;
    nuki[ 11 * kd + 93 ] += -1 ;
    nuki[ 1 * kd + 93 ] += +1 ;
    nuki[ 17 * kd + 93 ] += +1 ;

    /*reaction 95: OH + CH3 (+M) <=> CH3OH (+M) */
    nuki[ 4 * kd + 94 ] += -1 ;
    nuki[ 12 * kd + 94 ] += -1 ;
    nuki[ 20 * kd + 94 ] += +1 ;

    /*reaction 96: OH + CH3 <=> CH2 + H2O */
    nuki[ 4 * kd + 95 ] += -1 ;
    nuki[ 12 * kd + 95 ] += -1 ;
    nuki[ 10 * kd + 95 ] += +1 ;
    nuki[ 5 * kd + 95 ] += +1 ;

    /*reaction 97: OH + CH3 <=> CH2(S) + H2O */
    nuki[ 4 * kd + 96 ] += -1 ;
    nuki[ 12 * kd + 96 ] += -1 ;
    nuki[ 11 * kd + 96 ] += +1 ;
    nuki[ 5 * kd + 96 ] += +1 ;

    /*reaction 98: OH + CH4 <=> CH3 + H2O */
    nuki[ 4 * kd + 97 ] += -1 ;
    nuki[ 13 * kd + 97 ] += -1 ;
    nuki[ 12 * kd + 97 ] += +1 ;
    nuki[ 5 * kd + 97 ] += +1 ;

    /*reaction 99: OH + CO <=> H + CO2 */
    nuki[ 4 * kd + 98 ] += -1 ;
    nuki[ 14 * kd + 98 ] += -1 ;
    nuki[ 1 * kd + 98 ] += +1 ;
    nuki[ 15 * kd + 98 ] += +1 ;

    /*reaction 100: OH + HCO <=> H2O + CO */
    nuki[ 4 * kd + 99 ] += -1 ;
    nuki[ 16 * kd + 99 ] += -1 ;
    nuki[ 5 * kd + 99 ] += +1 ;
    nuki[ 14 * kd + 99 ] += +1 ;

    /*reaction 101: OH + CH2O <=> HCO + H2O */
    nuki[ 4 * kd + 100 ] += -1 ;
    nuki[ 17 * kd + 100 ] += -1 ;
    nuki[ 16 * kd + 100 ] += +1 ;
    nuki[ 5 * kd + 100 ] += +1 ;

    /*reaction 102: OH + CH2OH <=> H2O + CH2O */
    nuki[ 4 * kd + 101 ] += -1 ;
    nuki[ 18 * kd + 101 ] += -1 ;
    nuki[ 5 * kd + 101 ] += +1 ;
    nuki[ 17 * kd + 101 ] += +1 ;

    /*reaction 103: OH + CH3O <=> H2O + CH2O */
    nuki[ 4 * kd + 102 ] += -1 ;
    nuki[ 19 * kd + 102 ] += -1 ;
    nuki[ 5 * kd + 102 ] += +1 ;
    nuki[ 17 * kd + 102 ] += +1 ;

    /*reaction 104: OH + CH3OH <=> CH2OH + H2O */
    nuki[ 4 * kd + 103 ] += -1 ;
    nuki[ 20 * kd + 103 ] += -1 ;
    nuki[ 18 * kd + 103 ] += +1 ;
    nuki[ 5 * kd + 103 ] += +1 ;

    /*reaction 105: OH + CH3OH <=> CH3O + H2O */
    nuki[ 4 * kd + 104 ] += -1 ;
    nuki[ 20 * kd + 104 ] += -1 ;
    nuki[ 19 * kd + 104 ] += +1 ;
    nuki[ 5 * kd + 104 ] += +1 ;

    /*reaction 106: OH + C2H <=> H + HCCO */
    nuki[ 4 * kd + 105 ] += -1 ;
    nuki[ 21 * kd + 105 ] += -1 ;
    nuki[ 1 * kd + 105 ] += +1 ;
    nuki[ 27 * kd + 105 ] += +1 ;

    /*reaction 107: OH + C2H2 <=> H + CH2CO */
    nuki[ 4 * kd + 106 ] += -1 ;
    nuki[ 22 * kd + 106 ] += -1 ;
    nuki[ 1 * kd + 106 ] += +1 ;
    nuki[ 28 * kd + 106 ] += +1 ;

    /*reaction 108: OH + C2H2 <=> H + HCCOH */
    nuki[ 4 * kd + 107 ] += -1 ;
    nuki[ 22 * kd + 107 ] += -1 ;
    nuki[ 1 * kd + 107 ] += +1 ;
    nuki[ 29 * kd + 107 ] += +1 ;

    /*reaction 109: OH + C2H2 <=> C2H + H2O */
    nuki[ 4 * kd + 108 ] += -1 ;
    nuki[ 22 * kd + 108 ] += -1 ;
    nuki[ 21 * kd + 108 ] += +1 ;
    nuki[ 5 * kd + 108 ] += +1 ;

    /*reaction 110: OH + C2H2 <=> CH3 + CO */
    nuki[ 4 * kd + 109 ] += -1 ;
    nuki[ 22 * kd + 109 ] += -1 ;
    nuki[ 12 * kd + 109 ] += +1 ;
    nuki[ 14 * kd + 109 ] += +1 ;

    /*reaction 111: OH + C2H3 <=> H2O + C2H2 */
    nuki[ 4 * kd + 110 ] += -1 ;
    nuki[ 23 * kd + 110 ] += -1 ;
    nuki[ 5 * kd + 110 ] += +1 ;
    nuki[ 22 * kd + 110 ] += +1 ;

    /*reaction 112: OH + C2H4 <=> C2H3 + H2O */
    nuki[ 4 * kd + 111 ] += -1 ;
    nuki[ 24 * kd + 111 ] += -1 ;
    nuki[ 23 * kd + 111 ] += +1 ;
    nuki[ 5 * kd + 111 ] += +1 ;

    /*reaction 113: OH + C2H6 <=> C2H5 + H2O */
    nuki[ 4 * kd + 112 ] += -1 ;
    nuki[ 26 * kd + 112 ] += -1 ;
    nuki[ 25 * kd + 112 ] += +1 ;
    nuki[ 5 * kd + 112 ] += +1 ;

    /*reaction 114: OH + CH2CO <=> HCCO + H2O */
    nuki[ 4 * kd + 113 ] += -1 ;
    nuki[ 28 * kd + 113 ] += -1 ;
    nuki[ 27 * kd + 113 ] += +1 ;
    nuki[ 5 * kd + 113 ] += +1 ;

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    nuki[ 6 * kd + 114 ] += -2 ;
    nuki[ 3 * kd + 114 ] += +1 ;
    nuki[ 7 * kd + 114 ] += +1 ;

    /*reaction 116: 2 HO2 <=> O2 + H2O2 */
    nuki[ 6 * kd + 115 ] += -2 ;
    nuki[ 3 * kd + 115 ] += +1 ;
    nuki[ 7 * kd + 115 ] += +1 ;

    /*reaction 117: HO2 + CH2 <=> OH + CH2O */
    nuki[ 6 * kd + 116 ] += -1 ;
    nuki[ 10 * kd + 116 ] += -1 ;
    nuki[ 4 * kd + 116 ] += +1 ;
    nuki[ 17 * kd + 116 ] += +1 ;

    /*reaction 118: HO2 + CH3 <=> O2 + CH4 */
    nuki[ 6 * kd + 117 ] += -1 ;
    nuki[ 12 * kd + 117 ] += -1 ;
    nuki[ 3 * kd + 117 ] += +1 ;
    nuki[ 13 * kd + 117 ] += +1 ;

    /*reaction 119: HO2 + CH3 <=> OH + CH3O */
    nuki[ 6 * kd + 118 ] += -1 ;
    nuki[ 12 * kd + 118 ] += -1 ;
    nuki[ 4 * kd + 118 ] += +1 ;
    nuki[ 19 * kd + 118 ] += +1 ;

    /*reaction 120: HO2 + CO <=> OH + CO2 */
    nuki[ 6 * kd + 119 ] += -1 ;
    nuki[ 14 * kd + 119 ] += -1 ;
    nuki[ 4 * kd + 119 ] += +1 ;
    nuki[ 15 * kd + 119 ] += +1 ;

    /*reaction 121: HO2 + CH2O <=> HCO + H2O2 */
    nuki[ 6 * kd + 120 ] += -1 ;
    nuki[ 17 * kd + 120 ] += -1 ;
    nuki[ 16 * kd + 120 ] += +1 ;
    nuki[ 7 * kd + 120 ] += +1 ;

    /*reaction 122: C + O2 <=> O + CO */
    nuki[ 8 * kd + 121 ] += -1 ;
    nuki[ 3 * kd + 121 ] += -1 ;
    nuki[ 2 * kd + 121 ] += +1 ;
    nuki[ 14 * kd + 121 ] += +1 ;

    /*reaction 123: C + CH2 <=> H + C2H */
    nuki[ 8 * kd + 122 ] += -1 ;
    nuki[ 10 * kd + 122 ] += -1 ;
    nuki[ 1 * kd + 122 ] += +1 ;
    nuki[ 21 * kd + 122 ] += +1 ;

    /*reaction 124: C + CH3 <=> H + C2H2 */
    nuki[ 8 * kd + 123 ] += -1 ;
    nuki[ 12 * kd + 123 ] += -1 ;
    nuki[ 1 * kd + 123 ] += +1 ;
    nuki[ 22 * kd + 123 ] += +1 ;

    /*reaction 125: CH + O2 <=> O + HCO */
    nuki[ 9 * kd + 124 ] += -1 ;
    nuki[ 3 * kd + 124 ] += -1 ;
    nuki[ 2 * kd + 124 ] += +1 ;
    nuki[ 16 * kd + 124 ] += +1 ;

    /*reaction 126: CH + H2 <=> H + CH2 */
    nuki[ 9 * kd + 125 ] += -1 ;
    nuki[ 0 * kd + 125 ] += -1 ;
    nuki[ 1 * kd + 125 ] += +1 ;
    nuki[ 10 * kd + 125 ] += +1 ;

    /*reaction 127: CH + H2O <=> H + CH2O */
    nuki[ 9 * kd + 126 ] += -1 ;
    nuki[ 5 * kd + 126 ] += -1 ;
    nuki[ 1 * kd + 126 ] += +1 ;
    nuki[ 17 * kd + 126 ] += +1 ;

    /*reaction 128: CH + CH2 <=> H + C2H2 */
    nuki[ 9 * kd + 127 ] += -1 ;
    nuki[ 10 * kd + 127 ] += -1 ;
    nuki[ 1 * kd + 127 ] += +1 ;
    nuki[ 22 * kd + 127 ] += +1 ;

    /*reaction 129: CH + CH3 <=> H + C2H3 */
    nuki[ 9 * kd + 128 ] += -1 ;
    nuki[ 12 * kd + 128 ] += -1 ;
    nuki[ 1 * kd + 128 ] += +1 ;
    nuki[ 23 * kd + 128 ] += +1 ;

    /*reaction 130: CH + CH4 <=> H + C2H4 */
    nuki[ 9 * kd + 129 ] += -1 ;
    nuki[ 13 * kd + 129 ] += -1 ;
    nuki[ 1 * kd + 129 ] += +1 ;
    nuki[ 24 * kd + 129 ] += +1 ;

    /*reaction 131: CH + CO (+M) <=> HCCO (+M) */
    nuki[ 9 * kd + 130 ] += -1 ;
    nuki[ 14 * kd + 130 ] += -1 ;
    nuki[ 27 * kd + 130 ] += +1 ;

    /*reaction 132: CH + CO2 <=> HCO + CO */
    nuki[ 9 * kd + 131 ] += -1 ;
    nuki[ 15 * kd + 131 ] += -1 ;
    nuki[ 16 * kd + 131 ] += +1 ;
    nuki[ 14 * kd + 131 ] += +1 ;

    /*reaction 133: CH + CH2O <=> H + CH2CO */
    nuki[ 9 * kd + 132 ] += -1 ;
    nuki[ 17 * kd + 132 ] += -1 ;
    nuki[ 1 * kd + 132 ] += +1 ;
    nuki[ 28 * kd + 132 ] += +1 ;

    /*reaction 134: CH + HCCO <=> CO + C2H2 */
    nuki[ 9 * kd + 133 ] += -1 ;
    nuki[ 27 * kd + 133 ] += -1 ;
    nuki[ 14 * kd + 133 ] += +1 ;
    nuki[ 22 * kd + 133 ] += +1 ;

    /*reaction 135: CH2 + O2 <=> OH + HCO */
    nuki[ 10 * kd + 134 ] += -1 ;
    nuki[ 3 * kd + 134 ] += -1 ;
    nuki[ 4 * kd + 134 ] += +1 ;
    nuki[ 16 * kd + 134 ] += +1 ;

    /*reaction 136: CH2 + H2 <=> H + CH3 */
    nuki[ 10 * kd + 135 ] += -1 ;
    nuki[ 0 * kd + 135 ] += -1 ;
    nuki[ 1 * kd + 135 ] += +1 ;
    nuki[ 12 * kd + 135 ] += +1 ;

    /*reaction 137: 2 CH2 <=> H2 + C2H2 */
    nuki[ 10 * kd + 136 ] += -2 ;
    nuki[ 0 * kd + 136 ] += +1 ;
    nuki[ 22 * kd + 136 ] += +1 ;

    /*reaction 138: CH2 + CH3 <=> H + C2H4 */
    nuki[ 10 * kd + 137 ] += -1 ;
    nuki[ 12 * kd + 137 ] += -1 ;
    nuki[ 1 * kd + 137 ] += +1 ;
    nuki[ 24 * kd + 137 ] += +1 ;

    /*reaction 139: CH2 + CH4 <=> 2 CH3 */
    nuki[ 10 * kd + 138 ] += -1 ;
    nuki[ 13 * kd + 138 ] += -1 ;
    nuki[ 12 * kd + 138 ] += +2 ;

    /*reaction 140: CH2 + CO (+M) <=> CH2CO (+M) */
    nuki[ 10 * kd + 139 ] += -1 ;
    nuki[ 14 * kd + 139 ] += -1 ;
    nuki[ 28 * kd + 139 ] += +1 ;

    /*reaction 141: CH2 + HCCO <=> C2H3 + CO */
    nuki[ 10 * kd + 140 ] += -1 ;
    nuki[ 27 * kd + 140 ] += -1 ;
    nuki[ 23 * kd + 140 ] += +1 ;
    nuki[ 14 * kd + 140 ] += +1 ;

    /*reaction 142: CH2(S) + N2 <=> CH2 + N2 */
    nuki[ 11 * kd + 141 ] += -1 ;
    nuki[ 30 * kd + 141 ] += -1 ;
    nuki[ 10 * kd + 141 ] += +1 ;
    nuki[ 30 * kd + 141 ] += +1 ;

    /*reaction 143: CH2(S) + AR <=> CH2 + AR */
    nuki[ 11 * kd + 142 ] += -1 ;
    nuki[ 31 * kd + 142 ] += -1 ;
    nuki[ 10 * kd + 142 ] += +1 ;
    nuki[ 31 * kd + 142 ] += +1 ;

    /*reaction 144: CH2(S) + O2 <=> H + OH + CO */
    nuki[ 11 * kd + 143 ] += -1 ;
    nuki[ 3 * kd + 143 ] += -1 ;
    nuki[ 1 * kd + 143 ] += +1 ;
    nuki[ 4 * kd + 143 ] += +1 ;
    nuki[ 14 * kd + 143 ] += +1 ;

    /*reaction 145: CH2(S) + O2 <=> CO + H2O */
    nuki[ 11 * kd + 144 ] += -1 ;
    nuki[ 3 * kd + 144 ] += -1 ;
    nuki[ 14 * kd + 144 ] += +1 ;
    nuki[ 5 * kd + 144 ] += +1 ;

    /*reaction 146: CH2(S) + H2 <=> CH3 + H */
    nuki[ 11 * kd + 145 ] += -1 ;
    nuki[ 0 * kd + 145 ] += -1 ;
    nuki[ 12 * kd + 145 ] += +1 ;
    nuki[ 1 * kd + 145 ] += +1 ;

    /*reaction 147: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    nuki[ 11 * kd + 146 ] += -1 ;
    nuki[ 5 * kd + 146 ] += -1 ;
    nuki[ 20 * kd + 146 ] += +1 ;

    /*reaction 148: CH2(S) + H2O <=> CH2 + H2O */
    nuki[ 11 * kd + 147 ] += -1 ;
    nuki[ 5 * kd + 147 ] += -1 ;
    nuki[ 10 * kd + 147 ] += +1 ;
    nuki[ 5 * kd + 147 ] += +1 ;

    /*reaction 149: CH2(S) + CH3 <=> H + C2H4 */
    nuki[ 11 * kd + 148 ] += -1 ;
    nuki[ 12 * kd + 148 ] += -1 ;
    nuki[ 1 * kd + 148 ] += +1 ;
    nuki[ 24 * kd + 148 ] += +1 ;

    /*reaction 150: CH2(S) + CH4 <=> 2 CH3 */
    nuki[ 11 * kd + 149 ] += -1 ;
    nuki[ 13 * kd + 149 ] += -1 ;
    nuki[ 12 * kd + 149 ] += +2 ;

    /*reaction 151: CH2(S) + CO <=> CH2 + CO */
    nuki[ 11 * kd + 150 ] += -1 ;
    nuki[ 14 * kd + 150 ] += -1 ;
    nuki[ 10 * kd + 150 ] += +1 ;
    nuki[ 14 * kd + 150 ] += +1 ;

    /*reaction 152: CH2(S) + CO2 <=> CH2 + CO2 */
    nuki[ 11 * kd + 151 ] += -1 ;
    nuki[ 15 * kd + 151 ] += -1 ;
    nuki[ 10 * kd + 151 ] += +1 ;
    nuki[ 15 * kd + 151 ] += +1 ;

    /*reaction 153: CH2(S) + CO2 <=> CO + CH2O */
    nuki[ 11 * kd + 152 ] += -1 ;
    nuki[ 15 * kd + 152 ] += -1 ;
    nuki[ 14 * kd + 152 ] += +1 ;
    nuki[ 17 * kd + 152 ] += +1 ;

    /*reaction 154: CH2(S) + C2H6 <=> CH3 + C2H5 */
    nuki[ 11 * kd + 153 ] += -1 ;
    nuki[ 26 * kd + 153 ] += -1 ;
    nuki[ 12 * kd + 153 ] += +1 ;
    nuki[ 25 * kd + 153 ] += +1 ;

    /*reaction 155: CH3 + O2 <=> O + CH3O */
    nuki[ 12 * kd + 154 ] += -1 ;
    nuki[ 3 * kd + 154 ] += -1 ;
    nuki[ 2 * kd + 154 ] += +1 ;
    nuki[ 19 * kd + 154 ] += +1 ;

    /*reaction 156: CH3 + O2 <=> OH + CH2O */
    nuki[ 12 * kd + 155 ] += -1 ;
    nuki[ 3 * kd + 155 ] += -1 ;
    nuki[ 4 * kd + 155 ] += +1 ;
    nuki[ 17 * kd + 155 ] += +1 ;

    /*reaction 157: CH3 + H2O2 <=> HO2 + CH4 */
    nuki[ 12 * kd + 156 ] += -1 ;
    nuki[ 7 * kd + 156 ] += -1 ;
    nuki[ 6 * kd + 156 ] += +1 ;
    nuki[ 13 * kd + 156 ] += +1 ;

    /*reaction 158: 2 CH3 (+M) <=> C2H6 (+M) */
    nuki[ 12 * kd + 157 ] += -2 ;
    nuki[ 26 * kd + 157 ] += +1 ;

    /*reaction 159: 2 CH3 <=> H + C2H5 */
    nuki[ 12 * kd + 158 ] += -2 ;
    nuki[ 1 * kd + 158 ] += +1 ;
    nuki[ 25 * kd + 158 ] += +1 ;

    /*reaction 160: CH3 + HCO <=> CH4 + CO */
    nuki[ 12 * kd + 159 ] += -1 ;
    nuki[ 16 * kd + 159 ] += -1 ;
    nuki[ 13 * kd + 159 ] += +1 ;
    nuki[ 14 * kd + 159 ] += +1 ;

    /*reaction 161: CH3 + CH2O <=> HCO + CH4 */
    nuki[ 12 * kd + 160 ] += -1 ;
    nuki[ 17 * kd + 160 ] += -1 ;
    nuki[ 16 * kd + 160 ] += +1 ;
    nuki[ 13 * kd + 160 ] += +1 ;

    /*reaction 162: CH3 + CH3OH <=> CH2OH + CH4 */
    nuki[ 12 * kd + 161 ] += -1 ;
    nuki[ 20 * kd + 161 ] += -1 ;
    nuki[ 18 * kd + 161 ] += +1 ;
    nuki[ 13 * kd + 161 ] += +1 ;

    /*reaction 163: CH3 + CH3OH <=> CH3O + CH4 */
    nuki[ 12 * kd + 162 ] += -1 ;
    nuki[ 20 * kd + 162 ] += -1 ;
    nuki[ 19 * kd + 162 ] += +1 ;
    nuki[ 13 * kd + 162 ] += +1 ;

    /*reaction 164: CH3 + C2H4 <=> C2H3 + CH4 */
    nuki[ 12 * kd + 163 ] += -1 ;
    nuki[ 24 * kd + 163 ] += -1 ;
    nuki[ 23 * kd + 163 ] += +1 ;
    nuki[ 13 * kd + 163 ] += +1 ;

    /*reaction 165: CH3 + C2H6 <=> C2H5 + CH4 */
    nuki[ 12 * kd + 164 ] += -1 ;
    nuki[ 26 * kd + 164 ] += -1 ;
    nuki[ 25 * kd + 164 ] += +1 ;
    nuki[ 13 * kd + 164 ] += +1 ;

    /*reaction 166: HCO + H2O <=> H + CO + H2O */
    nuki[ 16 * kd + 165 ] += -1 ;
    nuki[ 5 * kd + 165 ] += -1 ;
    nuki[ 1 * kd + 165 ] += +1 ;
    nuki[ 14 * kd + 165 ] += +1 ;
    nuki[ 5 * kd + 165 ] += +1 ;

    /*reaction 167: HCO + M <=> H + CO + M */
    nuki[ 16 * kd + 166 ] += -1 ;
    nuki[ 1 * kd + 166 ] += +1 ;
    nuki[ 14 * kd + 166 ] += +1 ;

    /*reaction 168: HCO + O2 <=> HO2 + CO */
    nuki[ 16 * kd + 167 ] += -1 ;
    nuki[ 3 * kd + 167 ] += -1 ;
    nuki[ 6 * kd + 167 ] += +1 ;
    nuki[ 14 * kd + 167 ] += +1 ;

    /*reaction 169: CH2OH + O2 <=> HO2 + CH2O */
    nuki[ 18 * kd + 168 ] += -1 ;
    nuki[ 3 * kd + 168 ] += -1 ;
    nuki[ 6 * kd + 168 ] += +1 ;
    nuki[ 17 * kd + 168 ] += +1 ;

    /*reaction 170: CH3O + O2 <=> HO2 + CH2O */
    nuki[ 19 * kd + 169 ] += -1 ;
    nuki[ 3 * kd + 169 ] += -1 ;
    nuki[ 6 * kd + 169 ] += +1 ;
    nuki[ 17 * kd + 169 ] += +1 ;

    /*reaction 171: C2H + O2 <=> HCO + CO */
    nuki[ 21 * kd + 170 ] += -1 ;
    nuki[ 3 * kd + 170 ] += -1 ;
    nuki[ 16 * kd + 170 ] += +1 ;
    nuki[ 14 * kd + 170 ] += +1 ;

    /*reaction 172: C2H + H2 <=> H + C2H2 */
    nuki[ 21 * kd + 171 ] += -1 ;
    nuki[ 0 * kd + 171 ] += -1 ;
    nuki[ 1 * kd + 171 ] += +1 ;
    nuki[ 22 * kd + 171 ] += +1 ;

    /*reaction 173: C2H3 + O2 <=> HCO + CH2O */
    nuki[ 23 * kd + 172 ] += -1 ;
    nuki[ 3 * kd + 172 ] += -1 ;
    nuki[ 16 * kd + 172 ] += +1 ;
    nuki[ 17 * kd + 172 ] += +1 ;

    /*reaction 174: C2H4 (+M) <=> H2 + C2H2 (+M) */
    nuki[ 24 * kd + 173 ] += -1 ;
    nuki[ 0 * kd + 173 ] += +1 ;
    nuki[ 22 * kd + 173 ] += +1 ;

    /*reaction 175: C2H5 + O2 <=> HO2 + C2H4 */
    nuki[ 25 * kd + 174 ] += -1 ;
    nuki[ 3 * kd + 174 ] += -1 ;
    nuki[ 6 * kd + 174 ] += +1 ;
    nuki[ 24 * kd + 174 ] += +1 ;

    /*reaction 176: HCCO + O2 <=> OH + 2 CO */
    nuki[ 27 * kd + 175 ] += -1 ;
    nuki[ 3 * kd + 175 ] += -1 ;
    nuki[ 4 * kd + 175 ] += +1 ;
    nuki[ 14 * kd + 175 ] += +2 ;

    /*reaction 177: 2 HCCO <=> 2 CO + C2H2 */
    nuki[ 27 * kd + 176 ] += -2 ;
    nuki[ 14 * kd + 176 ] += +2 ;
    nuki[ 22 * kd + 176 ] += +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < 5 * 32; ++ id) {
         ncf[id] = 0; 
    }

    /*H2 */
    ncf[ 0 * kd + 1 ] = 2; /*H */

    /*H */
    ncf[ 1 * kd + 1 ] = 1; /*H */

    /*O */
    ncf[ 2 * kd + 0 ] = 1; /*O */

    /*O2 */
    ncf[ 3 * kd + 0 ] = 2; /*O */

    /*OH */
    ncf[ 4 * kd + 0 ] = 1; /*O */
    ncf[ 4 * kd + 1 ] = 1; /*H */

    /*H2O */
    ncf[ 5 * kd + 1 ] = 2; /*H */
    ncf[ 5 * kd + 0 ] = 1; /*O */

    /*HO2 */
    ncf[ 6 * kd + 1 ] = 1; /*H */
    ncf[ 6 * kd + 0 ] = 2; /*O */

    /*H2O2 */
    ncf[ 7 * kd + 1 ] = 2; /*H */
    ncf[ 7 * kd + 0 ] = 2; /*O */

    /*C */
    ncf[ 8 * kd + 2 ] = 1; /*C */

    /*CH */
    ncf[ 9 * kd + 2 ] = 1; /*C */
    ncf[ 9 * kd + 1 ] = 1; /*H */

    /*CH2 */
    ncf[ 10 * kd + 2 ] = 1; /*C */
    ncf[ 10 * kd + 1 ] = 2; /*H */

    /*CH2(S) */
    ncf[ 11 * kd + 2 ] = 1; /*C */
    ncf[ 11 * kd + 1 ] = 2; /*H */

    /*CH3 */
    ncf[ 12 * kd + 2 ] = 1; /*C */
    ncf[ 12 * kd + 1 ] = 3; /*H */

    /*CH4 */
    ncf[ 13 * kd + 2 ] = 1; /*C */
    ncf[ 13 * kd + 1 ] = 4; /*H */

    /*CO */
    ncf[ 14 * kd + 2 ] = 1; /*C */
    ncf[ 14 * kd + 0 ] = 1; /*O */

    /*CO2 */
    ncf[ 15 * kd + 2 ] = 1; /*C */
    ncf[ 15 * kd + 0 ] = 2; /*O */

    /*HCO */
    ncf[ 16 * kd + 1 ] = 1; /*H */
    ncf[ 16 * kd + 2 ] = 1; /*C */
    ncf[ 16 * kd + 0 ] = 1; /*O */

    /*CH2O */
    ncf[ 17 * kd + 1 ] = 2; /*H */
    ncf[ 17 * kd + 2 ] = 1; /*C */
    ncf[ 17 * kd + 0 ] = 1; /*O */

    /*CH2OH */
    ncf[ 18 * kd + 2 ] = 1; /*C */
    ncf[ 18 * kd + 1 ] = 3; /*H */
    ncf[ 18 * kd + 0 ] = 1; /*O */

    /*CH3O */
    ncf[ 19 * kd + 2 ] = 1; /*C */
    ncf[ 19 * kd + 1 ] = 3; /*H */
    ncf[ 19 * kd + 0 ] = 1; /*O */

    /*CH3OH */
    ncf[ 20 * kd + 2 ] = 1; /*C */
    ncf[ 20 * kd + 1 ] = 4; /*H */
    ncf[ 20 * kd + 0 ] = 1; /*O */

    /*C2H */
    ncf[ 21 * kd + 2 ] = 2; /*C */
    ncf[ 21 * kd + 1 ] = 1; /*H */

    /*C2H2 */
    ncf[ 22 * kd + 2 ] = 2; /*C */
    ncf[ 22 * kd + 1 ] = 2; /*H */

    /*C2H3 */
    ncf[ 23 * kd + 2 ] = 2; /*C */
    ncf[ 23 * kd + 1 ] = 3; /*H */

    /*C2H4 */
    ncf[ 24 * kd + 2 ] = 2; /*C */
    ncf[ 24 * kd + 1 ] = 4; /*H */

    /*C2H5 */
    ncf[ 25 * kd + 2 ] = 2; /*C */
    ncf[ 25 * kd + 1 ] = 5; /*H */

    /*C2H6 */
    ncf[ 26 * kd + 2 ] = 2; /*C */
    ncf[ 26 * kd + 1 ] = 6; /*H */

    /*HCCO */
    ncf[ 27 * kd + 1 ] = 1; /*H */
    ncf[ 27 * kd + 2 ] = 2; /*C */
    ncf[ 27 * kd + 0 ] = 1; /*O */

    /*CH2CO */
    ncf[ 28 * kd + 2 ] = 2; /*C */
    ncf[ 28 * kd + 1 ] = 2; /*H */
    ncf[ 28 * kd + 0 ] = 1; /*O */

    /*HCCOH */
    ncf[ 29 * kd + 2 ] = 2; /*C */
    ncf[ 29 * kd + 0 ] = 1; /*O */
    ncf[ 29 * kd + 1 ] = 2; /*H */

    /*N2 */
    ncf[ 30 * kd + 3 ] = 2; /*N */

    /*AR */
    ncf[ 31 * kd + 4 ] = 1; /*AR */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE(int * iwrk, double * rwrk, double * a, double * b, double * e)
{

    /*reaction 1: 2 O + M <=> O2 + M */
    a[0] = 1.2e+17;
    b[0] = -1;
    e[0] = 0;

    /*reaction 2: O + H + M <=> OH + M */
    a[1] = 5e+17;
    b[1] = -1;
    e[1] = 0;

    /*reaction 3: O + H2 <=> H + OH */
    a[2] = 50000;
    b[2] = 2.67;
    e[2] = 6290;

    /*reaction 4: O + HO2 <=> OH + O2 */
    a[3] = 2e+13;
    b[3] = 0;
    e[3] = 0;

    /*reaction 5: O + H2O2 <=> OH + HO2 */
    a[4] = 9.63e+06;
    b[4] = 2;
    e[4] = 4000;

    /*reaction 6: O + CH <=> H + CO */
    a[5] = 5.7e+13;
    b[5] = 0;
    e[5] = 0;

    /*reaction 7: O + CH2 <=> H + HCO */
    a[6] = 8e+13;
    b[6] = 0;
    e[6] = 0;

    /*reaction 8: O + CH2(S) <=> H2 + CO */
    a[7] = 1.5e+13;
    b[7] = 0;
    e[7] = 0;

    /*reaction 9: O + CH2(S) <=> H + HCO */
    a[8] = 1.5e+13;
    b[8] = 0;
    e[8] = 0;

    /*reaction 10: O + CH3 <=> H + CH2O */
    a[9] = 8.43e+13;
    b[9] = 0;
    e[9] = 0;

    /*reaction 11: O + CH4 <=> OH + CH3 */
    a[10] = 1.02e+09;
    b[10] = 1.5;
    e[10] = 8600;

    /*reaction 12: O + CO + M <=> CO2 + M */
    a[11] = 6.02e+14;
    b[11] = 0;
    e[11] = 3000;

    /*reaction 13: O + HCO <=> OH + CO */
    a[12] = 3e+13;
    b[12] = 0;
    e[12] = 0;

    /*reaction 14: O + HCO <=> H + CO2 */
    a[13] = 3e+13;
    b[13] = 0;
    e[13] = 0;

    /*reaction 15: O + CH2O <=> OH + HCO */
    a[14] = 3.9e+13;
    b[14] = 0;
    e[14] = 3540;

    /*reaction 16: O + CH2OH <=> OH + CH2O */
    a[15] = 1e+13;
    b[15] = 0;
    e[15] = 0;

    /*reaction 17: O + CH3O <=> OH + CH2O */
    a[16] = 1e+13;
    b[16] = 0;
    e[16] = 0;

    /*reaction 18: O + CH3OH <=> OH + CH2OH */
    a[17] = 388000;
    b[17] = 2.5;
    e[17] = 3100;

    /*reaction 19: O + CH3OH <=> OH + CH3O */
    a[18] = 130000;
    b[18] = 2.5;
    e[18] = 5000;

    /*reaction 20: O + C2H <=> CH + CO */
    a[19] = 5e+13;
    b[19] = 0;
    e[19] = 0;

    /*reaction 21: O + C2H2 <=> H + HCCO */
    a[20] = 1.02e+07;
    b[20] = 2;
    e[20] = 1900;

    /*reaction 22: O + C2H2 <=> OH + C2H */
    a[21] = 4.6e+19;
    b[21] = -1.41;
    e[21] = 28950;

    /*reaction 23: O + C2H2 <=> CO + CH2 */
    a[22] = 1.02e+07;
    b[22] = 2;
    e[22] = 1900;

    /*reaction 24: O + C2H3 <=> H + CH2CO */
    a[23] = 3e+13;
    b[23] = 0;
    e[23] = 0;

    /*reaction 25: O + C2H4 <=> CH3 + HCO */
    a[24] = 1.92e+07;
    b[24] = 1.83;
    e[24] = 220;

    /*reaction 26: O + C2H5 <=> CH3 + CH2O */
    a[25] = 1.32e+14;
    b[25] = 0;
    e[25] = 0;

    /*reaction 27: O + C2H6 <=> OH + C2H5 */
    a[26] = 8.98e+07;
    b[26] = 1.92;
    e[26] = 5690;

    /*reaction 28: O + HCCO <=> H + 2 CO */
    a[27] = 1e+14;
    b[27] = 0;
    e[27] = 0;

    /*reaction 29: O + CH2CO <=> OH + HCCO */
    a[28] = 1e+13;
    b[28] = 0;
    e[28] = 8000;

    /*reaction 30: O + CH2CO <=> CH2 + CO2 */
    a[29] = 1.75e+12;
    b[29] = 0;
    e[29] = 1350;

    /*reaction 31: O2 + CO <=> O + CO2 */
    a[30] = 2.5e+12;
    b[30] = 0;
    e[30] = 47800;

    /*reaction 32: O2 + CH2O <=> HO2 + HCO */
    a[31] = 1e+14;
    b[31] = 0;
    e[31] = 40000;

    /*reaction 33: H + O2 + M <=> HO2 + M */
    a[32] = 2.8e+18;
    b[32] = -0.86;
    e[32] = 0;

    /*reaction 34: H + 2 O2 <=> HO2 + O2 */
    a[33] = 3e+20;
    b[33] = -1.72;
    e[33] = 0;

    /*reaction 35: H + O2 + H2O <=> HO2 + H2O */
    a[34] = 9.38e+18;
    b[34] = -0.76;
    e[34] = 0;

    /*reaction 36: H + O2 + N2 <=> HO2 + N2 */
    a[35] = 3.75e+20;
    b[35] = -1.72;
    e[35] = 0;

    /*reaction 37: H + O2 + AR <=> HO2 + AR */
    a[36] = 7e+17;
    b[36] = -0.8;
    e[36] = 0;

    /*reaction 38: H + O2 <=> O + OH */
    a[37] = 8.3e+13;
    b[37] = 0;
    e[37] = 14413;

    /*reaction 39: 2 H + M <=> H2 + M */
    a[38] = 1e+18;
    b[38] = -1;
    e[38] = 0;

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    a[39] = 9e+16;
    b[39] = -0.6;
    e[39] = 0;

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    a[40] = 6e+19;
    b[40] = -1.25;
    e[40] = 0;

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    a[41] = 5.5e+20;
    b[41] = -2;
    e[41] = 0;

    /*reaction 43: H + OH + M <=> H2O + M */
    a[42] = 2.2e+22;
    b[42] = -2;
    e[42] = 0;

    /*reaction 44: H + HO2 <=> O + H2O */
    a[43] = 3.97e+12;
    b[43] = 0;
    e[43] = 671;

    /*reaction 45: H + HO2 <=> O2 + H2 */
    a[44] = 2.8e+13;
    b[44] = 0;
    e[44] = 1068;

    /*reaction 46: H + HO2 <=> 2 OH */
    a[45] = 1.34e+14;
    b[45] = 0;
    e[45] = 635;

    /*reaction 47: H + H2O2 <=> HO2 + H2 */
    a[46] = 1.21e+07;
    b[46] = 2;
    e[46] = 5200;

    /*reaction 48: H + H2O2 <=> OH + H2O */
    a[47] = 1e+13;
    b[47] = 0;
    e[47] = 3600;

    /*reaction 49: H + CH <=> C + H2 */
    a[48] = 1.1e+14;
    b[48] = 0;
    e[48] = 0;

    /*reaction 50: H + CH2 (+M) <=> CH3 (+M) */
    a[49] = 2.5e+16;
    b[49] = -0.8;
    e[49] = 0;

    /*reaction 51: H + CH2(S) <=> CH + H2 */
    a[50] = 3e+13;
    b[50] = 0;
    e[50] = 0;

    /*reaction 52: H + CH3 (+M) <=> CH4 (+M) */
    a[51] = 1.27e+16;
    b[51] = -0.63;
    e[51] = 383;

    /*reaction 53: H + CH4 <=> CH3 + H2 */
    a[52] = 6.6e+08;
    b[52] = 1.62;
    e[52] = 10840;

    /*reaction 54: H + HCO (+M) <=> CH2O (+M) */
    a[53] = 1.09e+12;
    b[53] = 0.48;
    e[53] = -260;

    /*reaction 55: H + HCO <=> H2 + CO */
    a[54] = 7.34e+13;
    b[54] = 0;
    e[54] = 0;

    /*reaction 56: H + CH2O (+M) <=> CH2OH (+M) */
    a[55] = 5.4e+11;
    b[55] = 0.454;
    e[55] = 3600;

    /*reaction 57: H + CH2O (+M) <=> CH3O (+M) */
    a[56] = 5.4e+11;
    b[56] = 0.454;
    e[56] = 2600;

    /*reaction 58: H + CH2O <=> HCO + H2 */
    a[57] = 2.3e+10;
    b[57] = 1.05;
    e[57] = 3275;

    /*reaction 59: H + CH2OH (+M) <=> CH3OH (+M) */
    a[58] = 1.8e+13;
    b[58] = 0;
    e[58] = 0;

    /*reaction 60: H + CH2OH <=> H2 + CH2O */
    a[59] = 2e+13;
    b[59] = 0;
    e[59] = 0;

    /*reaction 61: H + CH2OH <=> OH + CH3 */
    a[60] = 1.2e+13;
    b[60] = 0;
    e[60] = 0;

    /*reaction 62: H + CH2OH <=> CH2(S) + H2O */
    a[61] = 6e+12;
    b[61] = 0;
    e[61] = 0;

    /*reaction 63: H + CH3O (+M) <=> CH3OH (+M) */
    a[62] = 5e+13;
    b[62] = 0;
    e[62] = 0;

    /*reaction 64: H + CH3O <=> H + CH2OH */
    a[63] = 3.4e+06;
    b[63] = 1.6;
    e[63] = 0;

    /*reaction 65: H + CH3O <=> H2 + CH2O */
    a[64] = 2e+13;
    b[64] = 0;
    e[64] = 0;

    /*reaction 66: H + CH3O <=> OH + CH3 */
    a[65] = 3.2e+13;
    b[65] = 0;
    e[65] = 0;

    /*reaction 67: H + CH3O <=> CH2(S) + H2O */
    a[66] = 1.6e+13;
    b[66] = 0;
    e[66] = 0;

    /*reaction 68: H + CH3OH <=> CH2OH + H2 */
    a[67] = 1.7e+07;
    b[67] = 2.1;
    e[67] = 4870;

    /*reaction 69: H + CH3OH <=> CH3O + H2 */
    a[68] = 4.2e+06;
    b[68] = 2.1;
    e[68] = 4870;

    /*reaction 70: H + C2H (+M) <=> C2H2 (+M) */
    a[69] = 1e+17;
    b[69] = -1;
    e[69] = 0;

    /*reaction 71: H + C2H2 (+M) <=> C2H3 (+M) */
    a[70] = 5.6e+12;
    b[70] = 0;
    e[70] = 2400;

    /*reaction 72: H + C2H3 (+M) <=> C2H4 (+M) */
    a[71] = 6.08e+12;
    b[71] = 0.27;
    e[71] = 280;

    /*reaction 73: H + C2H3 <=> H2 + C2H2 */
    a[72] = 3e+13;
    b[72] = 0;
    e[72] = 0;

    /*reaction 74: H + C2H4 (+M) <=> C2H5 (+M) */
    a[73] = 1.08e+12;
    b[73] = 0.454;
    e[73] = 1820;

    /*reaction 75: H + C2H4 <=> C2H3 + H2 */
    a[74] = 1.325e+06;
    b[74] = 2.53;
    e[74] = 12240;

    /*reaction 76: H + C2H5 (+M) <=> C2H6 (+M) */
    a[75] = 5.21e+17;
    b[75] = -0.99;
    e[75] = 1580;

    /*reaction 77: H + C2H5 <=> H2 + C2H4 */
    a[76] = 2e+12;
    b[76] = 0;
    e[76] = 0;

    /*reaction 78: H + C2H6 <=> C2H5 + H2 */
    a[77] = 1.15e+08;
    b[77] = 1.9;
    e[77] = 7530;

    /*reaction 79: H + HCCO <=> CH2(S) + CO */
    a[78] = 1e+14;
    b[78] = 0;
    e[78] = 0;

    /*reaction 80: H + CH2CO <=> HCCO + H2 */
    a[79] = 5e+13;
    b[79] = 0;
    e[79] = 8000;

    /*reaction 81: H + CH2CO <=> CH3 + CO */
    a[80] = 1.13e+13;
    b[80] = 0;
    e[80] = 3428;

    /*reaction 82: H + HCCOH <=> H + CH2CO */
    a[81] = 1e+13;
    b[81] = 0;
    e[81] = 0;

    /*reaction 83: H2 + CO (+M) <=> CH2O (+M) */
    a[82] = 4.3e+07;
    b[82] = 1.5;
    e[82] = 79600;

    /*reaction 84: OH + H2 <=> H + H2O */
    a[83] = 2.16e+08;
    b[83] = 1.51;
    e[83] = 3430;

    /*reaction 85: 2 OH (+M) <=> H2O2 (+M) */
    a[84] = 7.4e+13;
    b[84] = -0.37;
    e[84] = 0;

    /*reaction 86: 2 OH <=> O + H2O */
    a[85] = 35700;
    b[85] = 2.4;
    e[85] = -2110;

    /*reaction 87: OH + HO2 <=> O2 + H2O */
    a[86] = 2.9e+13;
    b[86] = 0;
    e[86] = -500;

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    a[87] = 1.75e+12;
    b[87] = 0;
    e[87] = 320;

    /*reaction 89: OH + H2O2 <=> HO2 + H2O */
    a[88] = 5.8e+14;
    b[88] = 0;
    e[88] = 9560;

    /*reaction 90: OH + C <=> H + CO */
    a[89] = 5e+13;
    b[89] = 0;
    e[89] = 0;

    /*reaction 91: OH + CH <=> H + HCO */
    a[90] = 3e+13;
    b[90] = 0;
    e[90] = 0;

    /*reaction 92: OH + CH2 <=> H + CH2O */
    a[91] = 2e+13;
    b[91] = 0;
    e[91] = 0;

    /*reaction 93: OH + CH2 <=> CH + H2O */
    a[92] = 1.13e+07;
    b[92] = 2;
    e[92] = 3000;

    /*reaction 94: OH + CH2(S) <=> H + CH2O */
    a[93] = 3e+13;
    b[93] = 0;
    e[93] = 0;

    /*reaction 95: OH + CH3 (+M) <=> CH3OH (+M) */
    a[94] = 6.3e+13;
    b[94] = 0;
    e[94] = 0;

    /*reaction 96: OH + CH3 <=> CH2 + H2O */
    a[95] = 5.6e+07;
    b[95] = 1.6;
    e[95] = 5420;

    /*reaction 97: OH + CH3 <=> CH2(S) + H2O */
    a[96] = 2.501e+13;
    b[96] = 0;
    e[96] = 0;

    /*reaction 98: OH + CH4 <=> CH3 + H2O */
    a[97] = 1e+08;
    b[97] = 1.6;
    e[97] = 3120;

    /*reaction 99: OH + CO <=> H + CO2 */
    a[98] = 4.76e+07;
    b[98] = 1.228;
    e[98] = 70;

    /*reaction 100: OH + HCO <=> H2O + CO */
    a[99] = 5e+13;
    b[99] = 0;
    e[99] = 0;

    /*reaction 101: OH + CH2O <=> HCO + H2O */
    a[100] = 3.43e+09;
    b[100] = 1.18;
    e[100] = -447;

    /*reaction 102: OH + CH2OH <=> H2O + CH2O */
    a[101] = 5e+12;
    b[101] = 0;
    e[101] = 0;

    /*reaction 103: OH + CH3O <=> H2O + CH2O */
    a[102] = 5e+12;
    b[102] = 0;
    e[102] = 0;

    /*reaction 104: OH + CH3OH <=> CH2OH + H2O */
    a[103] = 1.44e+06;
    b[103] = 2;
    e[103] = -840;

    /*reaction 105: OH + CH3OH <=> CH3O + H2O */
    a[104] = 6.3e+06;
    b[104] = 2;
    e[104] = 1500;

    /*reaction 106: OH + C2H <=> H + HCCO */
    a[105] = 2e+13;
    b[105] = 0;
    e[105] = 0;

    /*reaction 107: OH + C2H2 <=> H + CH2CO */
    a[106] = 0.000218;
    b[106] = 4.5;
    e[106] = -1000;

    /*reaction 108: OH + C2H2 <=> H + HCCOH */
    a[107] = 504000;
    b[107] = 2.3;
    e[107] = 13500;

    /*reaction 109: OH + C2H2 <=> C2H + H2O */
    a[108] = 3.37e+07;
    b[108] = 2;
    e[108] = 14000;

    /*reaction 110: OH + C2H2 <=> CH3 + CO */
    a[109] = 0.000483;
    b[109] = 4;
    e[109] = -2000;

    /*reaction 111: OH + C2H3 <=> H2O + C2H2 */
    a[110] = 5e+12;
    b[110] = 0;
    e[110] = 0;

    /*reaction 112: OH + C2H4 <=> C2H3 + H2O */
    a[111] = 3.6e+06;
    b[111] = 2;
    e[111] = 2500;

    /*reaction 113: OH + C2H6 <=> C2H5 + H2O */
    a[112] = 3.54e+06;
    b[112] = 2.12;
    e[112] = 870;

    /*reaction 114: OH + CH2CO <=> HCCO + H2O */
    a[113] = 7.5e+12;
    b[113] = 0;
    e[113] = 2000;

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    a[114] = 1.3e+11;
    b[114] = 0;
    e[114] = -1630;

    /*reaction 116: 2 HO2 <=> O2 + H2O2 */
    a[115] = 4.2e+14;
    b[115] = 0;
    e[115] = 12000;

    /*reaction 117: HO2 + CH2 <=> OH + CH2O */
    a[116] = 2e+13;
    b[116] = 0;
    e[116] = 0;

    /*reaction 118: HO2 + CH3 <=> O2 + CH4 */
    a[117] = 1e+12;
    b[117] = 0;
    e[117] = 0;

    /*reaction 119: HO2 + CH3 <=> OH + CH3O */
    a[118] = 2e+13;
    b[118] = 0;
    e[118] = 0;

    /*reaction 120: HO2 + CO <=> OH + CO2 */
    a[119] = 1.5e+14;
    b[119] = 0;
    e[119] = 23600;

    /*reaction 121: HO2 + CH2O <=> HCO + H2O2 */
    a[120] = 1e+12;
    b[120] = 0;
    e[120] = 8000;

    /*reaction 122: C + O2 <=> O + CO */
    a[121] = 5.8e+13;
    b[121] = 0;
    e[121] = 576;

    /*reaction 123: C + CH2 <=> H + C2H */
    a[122] = 5e+13;
    b[122] = 0;
    e[122] = 0;

    /*reaction 124: C + CH3 <=> H + C2H2 */
    a[123] = 5e+13;
    b[123] = 0;
    e[123] = 0;

    /*reaction 125: CH + O2 <=> O + HCO */
    a[124] = 3.3e+13;
    b[124] = 0;
    e[124] = 0;

    /*reaction 126: CH + H2 <=> H + CH2 */
    a[125] = 1.107e+08;
    b[125] = 1.79;
    e[125] = 1670;

    /*reaction 127: CH + H2O <=> H + CH2O */
    a[126] = 5.71e+12;
    b[126] = 0;
    e[126] = -755;

    /*reaction 128: CH + CH2 <=> H + C2H2 */
    a[127] = 4e+13;
    b[127] = 0;
    e[127] = 0;

    /*reaction 129: CH + CH3 <=> H + C2H3 */
    a[128] = 3e+13;
    b[128] = 0;
    e[128] = 0;

    /*reaction 130: CH + CH4 <=> H + C2H4 */
    a[129] = 6e+13;
    b[129] = 0;
    e[129] = 0;

    /*reaction 131: CH + CO (+M) <=> HCCO (+M) */
    a[130] = 5e+13;
    b[130] = 0;
    e[130] = 0;

    /*reaction 132: CH + CO2 <=> HCO + CO */
    a[131] = 3.4e+12;
    b[131] = 0;
    e[131] = 690;

    /*reaction 133: CH + CH2O <=> H + CH2CO */
    a[132] = 9.46e+13;
    b[132] = 0;
    e[132] = -515;

    /*reaction 134: CH + HCCO <=> CO + C2H2 */
    a[133] = 5e+13;
    b[133] = 0;
    e[133] = 0;

    /*reaction 135: CH2 + O2 <=> OH + HCO */
    a[134] = 1.32e+13;
    b[134] = 0;
    e[134] = 1500;

    /*reaction 136: CH2 + H2 <=> H + CH3 */
    a[135] = 500000;
    b[135] = 2;
    e[135] = 7230;

    /*reaction 137: 2 CH2 <=> H2 + C2H2 */
    a[136] = 3.2e+13;
    b[136] = 0;
    e[136] = 0;

    /*reaction 138: CH2 + CH3 <=> H + C2H4 */
    a[137] = 4e+13;
    b[137] = 0;
    e[137] = 0;

    /*reaction 139: CH2 + CH4 <=> 2 CH3 */
    a[138] = 2.46e+06;
    b[138] = 2;
    e[138] = 8270;

    /*reaction 140: CH2 + CO (+M) <=> CH2CO (+M) */
    a[139] = 8.1e+11;
    b[139] = 0.5;
    e[139] = 4510;

    /*reaction 141: CH2 + HCCO <=> C2H3 + CO */
    a[140] = 3e+13;
    b[140] = 0;
    e[140] = 0;

    /*reaction 142: CH2(S) + N2 <=> CH2 + N2 */
    a[141] = 1.5e+13;
    b[141] = 0;
    e[141] = 600;

    /*reaction 143: CH2(S) + AR <=> CH2 + AR */
    a[142] = 9e+12;
    b[142] = 0;
    e[142] = 600;

    /*reaction 144: CH2(S) + O2 <=> H + OH + CO */
    a[143] = 2.8e+13;
    b[143] = 0;
    e[143] = 0;

    /*reaction 145: CH2(S) + O2 <=> CO + H2O */
    a[144] = 1.2e+13;
    b[144] = 0;
    e[144] = 0;

    /*reaction 146: CH2(S) + H2 <=> CH3 + H */
    a[145] = 7e+13;
    b[145] = 0;
    e[145] = 0;

    /*reaction 147: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    a[146] = 2e+13;
    b[146] = 0;
    e[146] = 0;

    /*reaction 148: CH2(S) + H2O <=> CH2 + H2O */
    a[147] = 3e+13;
    b[147] = 0;
    e[147] = 0;

    /*reaction 149: CH2(S) + CH3 <=> H + C2H4 */
    a[148] = 1.2e+13;
    b[148] = 0;
    e[148] = -570;

    /*reaction 150: CH2(S) + CH4 <=> 2 CH3 */
    a[149] = 1.6e+13;
    b[149] = 0;
    e[149] = -570;

    /*reaction 151: CH2(S) + CO <=> CH2 + CO */
    a[150] = 9e+12;
    b[150] = 0;
    e[150] = 0;

    /*reaction 152: CH2(S) + CO2 <=> CH2 + CO2 */
    a[151] = 7e+12;
    b[151] = 0;
    e[151] = 0;

    /*reaction 153: CH2(S) + CO2 <=> CO + CH2O */
    a[152] = 1.4e+13;
    b[152] = 0;
    e[152] = 0;

    /*reaction 154: CH2(S) + C2H6 <=> CH3 + C2H5 */
    a[153] = 4e+13;
    b[153] = 0;
    e[153] = -550;

    /*reaction 155: CH3 + O2 <=> O + CH3O */
    a[154] = 2.675e+13;
    b[154] = 0;
    e[154] = 28800;

    /*reaction 156: CH3 + O2 <=> OH + CH2O */
    a[155] = 3.6e+10;
    b[155] = 0;
    e[155] = 8940;

    /*reaction 157: CH3 + H2O2 <=> HO2 + CH4 */
    a[156] = 24500;
    b[156] = 2.47;
    e[156] = 5180;

    /*reaction 158: 2 CH3 (+M) <=> C2H6 (+M) */
    a[157] = 2.12e+16;
    b[157] = -0.97;
    e[157] = 620;

    /*reaction 159: 2 CH3 <=> H + C2H5 */
    a[158] = 4.99e+12;
    b[158] = 0.1;
    e[158] = 10600;

    /*reaction 160: CH3 + HCO <=> CH4 + CO */
    a[159] = 2.648e+13;
    b[159] = 0;
    e[159] = 0;

    /*reaction 161: CH3 + CH2O <=> HCO + CH4 */
    a[160] = 3320;
    b[160] = 2.81;
    e[160] = 5860;

    /*reaction 162: CH3 + CH3OH <=> CH2OH + CH4 */
    a[161] = 3e+07;
    b[161] = 1.5;
    e[161] = 9940;

    /*reaction 163: CH3 + CH3OH <=> CH3O + CH4 */
    a[162] = 1e+07;
    b[162] = 1.5;
    e[162] = 9940;

    /*reaction 164: CH3 + C2H4 <=> C2H3 + CH4 */
    a[163] = 227000;
    b[163] = 2;
    e[163] = 9200;

    /*reaction 165: CH3 + C2H6 <=> C2H5 + CH4 */
    a[164] = 6.14e+06;
    b[164] = 1.74;
    e[164] = 10450;

    /*reaction 166: HCO + H2O <=> H + CO + H2O */
    a[165] = 2.244e+18;
    b[165] = -1;
    e[165] = 17000;

    /*reaction 167: HCO + M <=> H + CO + M */
    a[166] = 1.87e+17;
    b[166] = -1;
    e[166] = 17000;

    /*reaction 168: HCO + O2 <=> HO2 + CO */
    a[167] = 7.6e+12;
    b[167] = 0;
    e[167] = 400;

    /*reaction 169: CH2OH + O2 <=> HO2 + CH2O */
    a[168] = 1.8e+13;
    b[168] = 0;
    e[168] = 900;

    /*reaction 170: CH3O + O2 <=> HO2 + CH2O */
    a[169] = 4.28e-13;
    b[169] = 7.6;
    e[169] = -3530;

    /*reaction 171: C2H + O2 <=> HCO + CO */
    a[170] = 5e+13;
    b[170] = 0;
    e[170] = 1500;

    /*reaction 172: C2H + H2 <=> H + C2H2 */
    a[171] = 407000;
    b[171] = 2.4;
    e[171] = 200;

    /*reaction 173: C2H3 + O2 <=> HCO + CH2O */
    a[172] = 3.98e+12;
    b[172] = 0;
    e[172] = -240;

    /*reaction 174: C2H4 (+M) <=> H2 + C2H2 (+M) */
    a[173] = 8e+12;
    b[173] = 0.44;
    e[173] = 88770;

    /*reaction 175: C2H5 + O2 <=> HO2 + C2H4 */
    a[174] = 8.4e+11;
    b[174] = 0;
    e[174] = 3875;

    /*reaction 176: HCCO + O2 <=> OH + 2 CO */
    a[175] = 1.6e+12;
    b[175] = 0;
    e[175] = 854;

    /*reaction 177: 2 HCCO <=> 2 CO + C2H2 */
    a[176] = 1e+13;
    b[176] = 0;
    e[176] = 0;

    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double * T, double * C, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[32]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 O + M <=> O2 + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: O + H + M <=> OH + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: O + H2 <=> H + OH */
    /*eqcon[2] *= 1;  */

    /*reaction 4: O + HO2 <=> OH + O2 */
    /*eqcon[3] *= 1;  */

    /*reaction 5: O + H2O2 <=> OH + HO2 */
    /*eqcon[4] *= 1;  */

    /*reaction 6: O + CH <=> H + CO */
    /*eqcon[5] *= 1;  */

    /*reaction 7: O + CH2 <=> H + HCO */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + CH2(S) <=> H2 + CO */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + CH2(S) <=> H + HCO */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + CH3 <=> H + CH2O */
    /*eqcon[9] *= 1;  */

    /*reaction 11: O + CH4 <=> OH + CH3 */
    /*eqcon[10] *= 1;  */

    /*reaction 12: O + CO + M <=> CO2 + M */
    eqcon[11] *= 1e+06; 

    /*reaction 13: O + HCO <=> OH + CO */
    /*eqcon[12] *= 1;  */

    /*reaction 14: O + HCO <=> H + CO2 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O + CH2O <=> OH + HCO */
    /*eqcon[14] *= 1;  */

    /*reaction 16: O + CH2OH <=> OH + CH2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: O + CH3O <=> OH + CH2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: O + CH3OH <=> OH + CH2OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: O + CH3OH <=> OH + CH3O */
    /*eqcon[18] *= 1;  */

    /*reaction 20: O + C2H <=> CH + CO */
    /*eqcon[19] *= 1;  */

    /*reaction 21: O + C2H2 <=> H + HCCO */
    /*eqcon[20] *= 1;  */

    /*reaction 22: O + C2H2 <=> OH + C2H */
    /*eqcon[21] *= 1;  */

    /*reaction 23: O + C2H2 <=> CO + CH2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: O + C2H3 <=> H + CH2CO */
    /*eqcon[23] *= 1;  */

    /*reaction 25: O + C2H4 <=> CH3 + HCO */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O + C2H5 <=> CH3 + CH2O */
    /*eqcon[25] *= 1;  */

    /*reaction 27: O + C2H6 <=> OH + C2H5 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: O + HCCO <=> H + 2 CO */
    eqcon[27] *= 1e-06; 

    /*reaction 29: O + CH2CO <=> OH + HCCO */
    /*eqcon[28] *= 1;  */

    /*reaction 30: O + CH2CO <=> CH2 + CO2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: O2 + CO <=> O + CO2 */
    /*eqcon[30] *= 1;  */

    /*reaction 32: O2 + CH2O <=> HO2 + HCO */
    /*eqcon[31] *= 1;  */

    /*reaction 33: H + O2 + M <=> HO2 + M */
    eqcon[32] *= 1e+06; 

    /*reaction 34: H + 2 O2 <=> HO2 + O2 */
    eqcon[33] *= 1e+06; 

    /*reaction 35: H + O2 + H2O <=> HO2 + H2O */
    eqcon[34] *= 1e+06; 

    /*reaction 36: H + O2 + N2 <=> HO2 + N2 */
    eqcon[35] *= 1e+06; 

    /*reaction 37: H + O2 + AR <=> HO2 + AR */
    eqcon[36] *= 1e+06; 

    /*reaction 38: H + O2 <=> O + OH */
    /*eqcon[37] *= 1;  */

    /*reaction 39: 2 H + M <=> H2 + M */
    eqcon[38] *= 1e+06; 

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    eqcon[39] *= 1e+06; 

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    eqcon[40] *= 1e+06; 

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    eqcon[41] *= 1e+06; 

    /*reaction 43: H + OH + M <=> H2O + M */
    eqcon[42] *= 1e+06; 

    /*reaction 44: H + HO2 <=> O + H2O */
    /*eqcon[43] *= 1;  */

    /*reaction 45: H + HO2 <=> O2 + H2 */
    /*eqcon[44] *= 1;  */

    /*reaction 46: H + HO2 <=> 2 OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + H2O2 <=> HO2 + H2 */
    /*eqcon[46] *= 1;  */

    /*reaction 48: H + H2O2 <=> OH + H2O */
    /*eqcon[47] *= 1;  */

    /*reaction 49: H + CH <=> C + H2 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[49] *= 1e+06; 

    /*reaction 51: H + CH2(S) <=> CH + H2 */
    /*eqcon[50] *= 1;  */

    /*reaction 52: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[51] *= 1e+06; 

    /*reaction 53: H + CH4 <=> CH3 + H2 */
    /*eqcon[52] *= 1;  */

    /*reaction 54: H + HCO (+M) <=> CH2O (+M) */
    eqcon[53] *= 1e+06; 

    /*reaction 55: H + HCO <=> H2 + CO */
    /*eqcon[54] *= 1;  */

    /*reaction 56: H + CH2O (+M) <=> CH2OH (+M) */
    eqcon[55] *= 1e+06; 

    /*reaction 57: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[56] *= 1e+06; 

    /*reaction 58: H + CH2O <=> HCO + H2 */
    /*eqcon[57] *= 1;  */

    /*reaction 59: H + CH2OH (+M) <=> CH3OH (+M) */
    eqcon[58] *= 1e+06; 

    /*reaction 60: H + CH2OH <=> H2 + CH2O */
    /*eqcon[59] *= 1;  */

    /*reaction 61: H + CH2OH <=> OH + CH3 */
    /*eqcon[60] *= 1;  */

    /*reaction 62: H + CH2OH <=> CH2(S) + H2O */
    /*eqcon[61] *= 1;  */

    /*reaction 63: H + CH3O (+M) <=> CH3OH (+M) */
    eqcon[62] *= 1e+06; 

    /*reaction 64: H + CH3O <=> H + CH2OH */
    /*eqcon[63] *= 1;  */

    /*reaction 65: H + CH3O <=> H2 + CH2O */
    /*eqcon[64] *= 1;  */

    /*reaction 66: H + CH3O <=> OH + CH3 */
    /*eqcon[65] *= 1;  */

    /*reaction 67: H + CH3O <=> CH2(S) + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: H + CH3OH <=> CH2OH + H2 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: H + CH3OH <=> CH3O + H2 */
    /*eqcon[68] *= 1;  */

    /*reaction 70: H + C2H (+M) <=> C2H2 (+M) */
    eqcon[69] *= 1e+06; 

    /*reaction 71: H + C2H2 (+M) <=> C2H3 (+M) */
    eqcon[70] *= 1e+06; 

    /*reaction 72: H + C2H3 (+M) <=> C2H4 (+M) */
    eqcon[71] *= 1e+06; 

    /*reaction 73: H + C2H3 <=> H2 + C2H2 */
    /*eqcon[72] *= 1;  */

    /*reaction 74: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[73] *= 1e+06; 

    /*reaction 75: H + C2H4 <=> C2H3 + H2 */
    /*eqcon[74] *= 1;  */

    /*reaction 76: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[75] *= 1e+06; 

    /*reaction 77: H + C2H5 <=> H2 + C2H4 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[77] *= 1;  */

    /*reaction 79: H + HCCO <=> CH2(S) + CO */
    /*eqcon[78] *= 1;  */

    /*reaction 80: H + CH2CO <=> HCCO + H2 */
    /*eqcon[79] *= 1;  */

    /*reaction 81: H + CH2CO <=> CH3 + CO */
    /*eqcon[80] *= 1;  */

    /*reaction 82: H + HCCOH <=> H + CH2CO */
    /*eqcon[81] *= 1;  */

    /*reaction 83: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[82] *= 1e+06; 

    /*reaction 84: OH + H2 <=> H + H2O */
    /*eqcon[83] *= 1;  */

    /*reaction 85: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[84] *= 1e+06; 

    /*reaction 86: 2 OH <=> O + H2O */
    /*eqcon[85] *= 1;  */

    /*reaction 87: OH + HO2 <=> O2 + H2O */
    /*eqcon[86] *= 1;  */

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[87] *= 1;  */

    /*reaction 89: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[88] *= 1;  */

    /*reaction 90: OH + C <=> H + CO */
    /*eqcon[89] *= 1;  */

    /*reaction 91: OH + CH <=> H + HCO */
    /*eqcon[90] *= 1;  */

    /*reaction 92: OH + CH2 <=> H + CH2O */
    /*eqcon[91] *= 1;  */

    /*reaction 93: OH + CH2 <=> CH + H2O */
    /*eqcon[92] *= 1;  */

    /*reaction 94: OH + CH2(S) <=> H + CH2O */
    /*eqcon[93] *= 1;  */

    /*reaction 95: OH + CH3 (+M) <=> CH3OH (+M) */
    eqcon[94] *= 1e+06; 

    /*reaction 96: OH + CH3 <=> CH2 + H2O */
    /*eqcon[95] *= 1;  */

    /*reaction 97: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[96] *= 1;  */

    /*reaction 98: OH + CH4 <=> CH3 + H2O */
    /*eqcon[97] *= 1;  */

    /*reaction 99: OH + CO <=> H + CO2 */
    /*eqcon[98] *= 1;  */

    /*reaction 100: OH + HCO <=> H2O + CO */
    /*eqcon[99] *= 1;  */

    /*reaction 101: OH + CH2O <=> HCO + H2O */
    /*eqcon[100] *= 1;  */

    /*reaction 102: OH + CH2OH <=> H2O + CH2O */
    /*eqcon[101] *= 1;  */

    /*reaction 103: OH + CH3O <=> H2O + CH2O */
    /*eqcon[102] *= 1;  */

    /*reaction 104: OH + CH3OH <=> CH2OH + H2O */
    /*eqcon[103] *= 1;  */

    /*reaction 105: OH + CH3OH <=> CH3O + H2O */
    /*eqcon[104] *= 1;  */

    /*reaction 106: OH + C2H <=> H + HCCO */
    /*eqcon[105] *= 1;  */

    /*reaction 107: OH + C2H2 <=> H + CH2CO */
    /*eqcon[106] *= 1;  */

    /*reaction 108: OH + C2H2 <=> H + HCCOH */
    /*eqcon[107] *= 1;  */

    /*reaction 109: OH + C2H2 <=> C2H + H2O */
    /*eqcon[108] *= 1;  */

    /*reaction 110: OH + C2H2 <=> CH3 + CO */
    /*eqcon[109] *= 1;  */

    /*reaction 111: OH + C2H3 <=> H2O + C2H2 */
    /*eqcon[110] *= 1;  */

    /*reaction 112: OH + C2H4 <=> C2H3 + H2O */
    /*eqcon[111] *= 1;  */

    /*reaction 113: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[112] *= 1;  */

    /*reaction 114: OH + CH2CO <=> HCCO + H2O */
    /*eqcon[113] *= 1;  */

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[114] *= 1;  */

    /*reaction 116: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[115] *= 1;  */

    /*reaction 117: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[116] *= 1;  */

    /*reaction 118: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[117] *= 1;  */

    /*reaction 119: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[118] *= 1;  */

    /*reaction 120: HO2 + CO <=> OH + CO2 */
    /*eqcon[119] *= 1;  */

    /*reaction 121: HO2 + CH2O <=> HCO + H2O2 */
    /*eqcon[120] *= 1;  */

    /*reaction 122: C + O2 <=> O + CO */
    /*eqcon[121] *= 1;  */

    /*reaction 123: C + CH2 <=> H + C2H */
    /*eqcon[122] *= 1;  */

    /*reaction 124: C + CH3 <=> H + C2H2 */
    /*eqcon[123] *= 1;  */

    /*reaction 125: CH + O2 <=> O + HCO */
    /*eqcon[124] *= 1;  */

    /*reaction 126: CH + H2 <=> H + CH2 */
    /*eqcon[125] *= 1;  */

    /*reaction 127: CH + H2O <=> H + CH2O */
    /*eqcon[126] *= 1;  */

    /*reaction 128: CH + CH2 <=> H + C2H2 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: CH + CH3 <=> H + C2H3 */
    /*eqcon[128] *= 1;  */

    /*reaction 130: CH + CH4 <=> H + C2H4 */
    /*eqcon[129] *= 1;  */

    /*reaction 131: CH + CO (+M) <=> HCCO (+M) */
    eqcon[130] *= 1e+06; 

    /*reaction 132: CH + CO2 <=> HCO + CO */
    /*eqcon[131] *= 1;  */

    /*reaction 133: CH + CH2O <=> H + CH2CO */
    /*eqcon[132] *= 1;  */

    /*reaction 134: CH + HCCO <=> CO + C2H2 */
    /*eqcon[133] *= 1;  */

    /*reaction 135: CH2 + O2 <=> OH + HCO */
    /*eqcon[134] *= 1;  */

    /*reaction 136: CH2 + H2 <=> H + CH3 */
    /*eqcon[135] *= 1;  */

    /*reaction 137: 2 CH2 <=> H2 + C2H2 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[137] *= 1;  */

    /*reaction 139: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[138] *= 1;  */

    /*reaction 140: CH2 + CO (+M) <=> CH2CO (+M) */
    eqcon[139] *= 1e+06; 

    /*reaction 141: CH2 + HCCO <=> C2H3 + CO */
    /*eqcon[140] *= 1;  */

    /*reaction 142: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[141] *= 1;  */

    /*reaction 143: CH2(S) + AR <=> CH2 + AR */
    /*eqcon[142] *= 1;  */

    /*reaction 144: CH2(S) + O2 <=> H + OH + CO */
    eqcon[143] *= 1e-06; 

    /*reaction 145: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[144] *= 1;  */

    /*reaction 146: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[145] *= 1;  */

    /*reaction 147: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    eqcon[146] *= 1e+06; 

    /*reaction 148: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[147] *= 1;  */

    /*reaction 149: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[148] *= 1;  */

    /*reaction 150: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[149] *= 1;  */

    /*reaction 151: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[150] *= 1;  */

    /*reaction 152: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[152] *= 1;  */

    /*reaction 154: CH2(S) + C2H6 <=> CH3 + C2H5 */
    /*eqcon[153] *= 1;  */

    /*reaction 155: CH3 + O2 <=> O + CH3O */
    /*eqcon[154] *= 1;  */

    /*reaction 156: CH3 + O2 <=> OH + CH2O */
    /*eqcon[155] *= 1;  */

    /*reaction 157: CH3 + H2O2 <=> HO2 + CH4 */
    /*eqcon[156] *= 1;  */

    /*reaction 158: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[157] *= 1e+06; 

    /*reaction 159: 2 CH3 <=> H + C2H5 */
    /*eqcon[158] *= 1;  */

    /*reaction 160: CH3 + HCO <=> CH4 + CO */
    /*eqcon[159] *= 1;  */

    /*reaction 161: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[160] *= 1;  */

    /*reaction 162: CH3 + CH3OH <=> CH2OH + CH4 */
    /*eqcon[161] *= 1;  */

    /*reaction 163: CH3 + CH3OH <=> CH3O + CH4 */
    /*eqcon[162] *= 1;  */

    /*reaction 164: CH3 + C2H4 <=> C2H3 + CH4 */
    /*eqcon[163] *= 1;  */

    /*reaction 165: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[164] *= 1;  */

    /*reaction 166: HCO + H2O <=> H + CO + H2O */
    eqcon[165] *= 1e-06; 

    /*reaction 167: HCO + M <=> H + CO + M */
    eqcon[166] *= 1e-06; 

    /*reaction 168: HCO + O2 <=> HO2 + CO */
    /*eqcon[167] *= 1;  */

    /*reaction 169: CH2OH + O2 <=> HO2 + CH2O */
    /*eqcon[168] *= 1;  */

    /*reaction 170: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[169] *= 1;  */

    /*reaction 171: C2H + O2 <=> HCO + CO */
    /*eqcon[170] *= 1;  */

    /*reaction 172: C2H + H2 <=> H + C2H2 */
    /*eqcon[171] *= 1;  */

    /*reaction 173: C2H3 + O2 <=> HCO + CH2O */
    /*eqcon[172] *= 1;  */

    /*reaction 174: C2H4 (+M) <=> H2 + C2H2 (+M) */
    eqcon[173] *= 1e-06; 

    /*reaction 175: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[174] *= 1;  */

    /*reaction 176: HCCO + O2 <=> OH + 2 CO */
    eqcon[175] *= 1e-06; 

    /*reaction 177: 2 HCCO <=> 2 CO + C2H2 */
    eqcon[176] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[32]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 O + M <=> O2 + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: O + H + M <=> OH + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: O + H2 <=> H + OH */
    /*eqcon[2] *= 1;  */

    /*reaction 4: O + HO2 <=> OH + O2 */
    /*eqcon[3] *= 1;  */

    /*reaction 5: O + H2O2 <=> OH + HO2 */
    /*eqcon[4] *= 1;  */

    /*reaction 6: O + CH <=> H + CO */
    /*eqcon[5] *= 1;  */

    /*reaction 7: O + CH2 <=> H + HCO */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + CH2(S) <=> H2 + CO */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + CH2(S) <=> H + HCO */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + CH3 <=> H + CH2O */
    /*eqcon[9] *= 1;  */

    /*reaction 11: O + CH4 <=> OH + CH3 */
    /*eqcon[10] *= 1;  */

    /*reaction 12: O + CO + M <=> CO2 + M */
    eqcon[11] *= 1e+06; 

    /*reaction 13: O + HCO <=> OH + CO */
    /*eqcon[12] *= 1;  */

    /*reaction 14: O + HCO <=> H + CO2 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O + CH2O <=> OH + HCO */
    /*eqcon[14] *= 1;  */

    /*reaction 16: O + CH2OH <=> OH + CH2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: O + CH3O <=> OH + CH2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: O + CH3OH <=> OH + CH2OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: O + CH3OH <=> OH + CH3O */
    /*eqcon[18] *= 1;  */

    /*reaction 20: O + C2H <=> CH + CO */
    /*eqcon[19] *= 1;  */

    /*reaction 21: O + C2H2 <=> H + HCCO */
    /*eqcon[20] *= 1;  */

    /*reaction 22: O + C2H2 <=> OH + C2H */
    /*eqcon[21] *= 1;  */

    /*reaction 23: O + C2H2 <=> CO + CH2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: O + C2H3 <=> H + CH2CO */
    /*eqcon[23] *= 1;  */

    /*reaction 25: O + C2H4 <=> CH3 + HCO */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O + C2H5 <=> CH3 + CH2O */
    /*eqcon[25] *= 1;  */

    /*reaction 27: O + C2H6 <=> OH + C2H5 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: O + HCCO <=> H + 2 CO */
    eqcon[27] *= 1e-06; 

    /*reaction 29: O + CH2CO <=> OH + HCCO */
    /*eqcon[28] *= 1;  */

    /*reaction 30: O + CH2CO <=> CH2 + CO2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: O2 + CO <=> O + CO2 */
    /*eqcon[30] *= 1;  */

    /*reaction 32: O2 + CH2O <=> HO2 + HCO */
    /*eqcon[31] *= 1;  */

    /*reaction 33: H + O2 + M <=> HO2 + M */
    eqcon[32] *= 1e+06; 

    /*reaction 34: H + 2 O2 <=> HO2 + O2 */
    eqcon[33] *= 1e+06; 

    /*reaction 35: H + O2 + H2O <=> HO2 + H2O */
    eqcon[34] *= 1e+06; 

    /*reaction 36: H + O2 + N2 <=> HO2 + N2 */
    eqcon[35] *= 1e+06; 

    /*reaction 37: H + O2 + AR <=> HO2 + AR */
    eqcon[36] *= 1e+06; 

    /*reaction 38: H + O2 <=> O + OH */
    /*eqcon[37] *= 1;  */

    /*reaction 39: 2 H + M <=> H2 + M */
    eqcon[38] *= 1e+06; 

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    eqcon[39] *= 1e+06; 

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    eqcon[40] *= 1e+06; 

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    eqcon[41] *= 1e+06; 

    /*reaction 43: H + OH + M <=> H2O + M */
    eqcon[42] *= 1e+06; 

    /*reaction 44: H + HO2 <=> O + H2O */
    /*eqcon[43] *= 1;  */

    /*reaction 45: H + HO2 <=> O2 + H2 */
    /*eqcon[44] *= 1;  */

    /*reaction 46: H + HO2 <=> 2 OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + H2O2 <=> HO2 + H2 */
    /*eqcon[46] *= 1;  */

    /*reaction 48: H + H2O2 <=> OH + H2O */
    /*eqcon[47] *= 1;  */

    /*reaction 49: H + CH <=> C + H2 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[49] *= 1e+06; 

    /*reaction 51: H + CH2(S) <=> CH + H2 */
    /*eqcon[50] *= 1;  */

    /*reaction 52: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[51] *= 1e+06; 

    /*reaction 53: H + CH4 <=> CH3 + H2 */
    /*eqcon[52] *= 1;  */

    /*reaction 54: H + HCO (+M) <=> CH2O (+M) */
    eqcon[53] *= 1e+06; 

    /*reaction 55: H + HCO <=> H2 + CO */
    /*eqcon[54] *= 1;  */

    /*reaction 56: H + CH2O (+M) <=> CH2OH (+M) */
    eqcon[55] *= 1e+06; 

    /*reaction 57: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[56] *= 1e+06; 

    /*reaction 58: H + CH2O <=> HCO + H2 */
    /*eqcon[57] *= 1;  */

    /*reaction 59: H + CH2OH (+M) <=> CH3OH (+M) */
    eqcon[58] *= 1e+06; 

    /*reaction 60: H + CH2OH <=> H2 + CH2O */
    /*eqcon[59] *= 1;  */

    /*reaction 61: H + CH2OH <=> OH + CH3 */
    /*eqcon[60] *= 1;  */

    /*reaction 62: H + CH2OH <=> CH2(S) + H2O */
    /*eqcon[61] *= 1;  */

    /*reaction 63: H + CH3O (+M) <=> CH3OH (+M) */
    eqcon[62] *= 1e+06; 

    /*reaction 64: H + CH3O <=> H + CH2OH */
    /*eqcon[63] *= 1;  */

    /*reaction 65: H + CH3O <=> H2 + CH2O */
    /*eqcon[64] *= 1;  */

    /*reaction 66: H + CH3O <=> OH + CH3 */
    /*eqcon[65] *= 1;  */

    /*reaction 67: H + CH3O <=> CH2(S) + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: H + CH3OH <=> CH2OH + H2 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: H + CH3OH <=> CH3O + H2 */
    /*eqcon[68] *= 1;  */

    /*reaction 70: H + C2H (+M) <=> C2H2 (+M) */
    eqcon[69] *= 1e+06; 

    /*reaction 71: H + C2H2 (+M) <=> C2H3 (+M) */
    eqcon[70] *= 1e+06; 

    /*reaction 72: H + C2H3 (+M) <=> C2H4 (+M) */
    eqcon[71] *= 1e+06; 

    /*reaction 73: H + C2H3 <=> H2 + C2H2 */
    /*eqcon[72] *= 1;  */

    /*reaction 74: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[73] *= 1e+06; 

    /*reaction 75: H + C2H4 <=> C2H3 + H2 */
    /*eqcon[74] *= 1;  */

    /*reaction 76: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[75] *= 1e+06; 

    /*reaction 77: H + C2H5 <=> H2 + C2H4 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[77] *= 1;  */

    /*reaction 79: H + HCCO <=> CH2(S) + CO */
    /*eqcon[78] *= 1;  */

    /*reaction 80: H + CH2CO <=> HCCO + H2 */
    /*eqcon[79] *= 1;  */

    /*reaction 81: H + CH2CO <=> CH3 + CO */
    /*eqcon[80] *= 1;  */

    /*reaction 82: H + HCCOH <=> H + CH2CO */
    /*eqcon[81] *= 1;  */

    /*reaction 83: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[82] *= 1e+06; 

    /*reaction 84: OH + H2 <=> H + H2O */
    /*eqcon[83] *= 1;  */

    /*reaction 85: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[84] *= 1e+06; 

    /*reaction 86: 2 OH <=> O + H2O */
    /*eqcon[85] *= 1;  */

    /*reaction 87: OH + HO2 <=> O2 + H2O */
    /*eqcon[86] *= 1;  */

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[87] *= 1;  */

    /*reaction 89: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[88] *= 1;  */

    /*reaction 90: OH + C <=> H + CO */
    /*eqcon[89] *= 1;  */

    /*reaction 91: OH + CH <=> H + HCO */
    /*eqcon[90] *= 1;  */

    /*reaction 92: OH + CH2 <=> H + CH2O */
    /*eqcon[91] *= 1;  */

    /*reaction 93: OH + CH2 <=> CH + H2O */
    /*eqcon[92] *= 1;  */

    /*reaction 94: OH + CH2(S) <=> H + CH2O */
    /*eqcon[93] *= 1;  */

    /*reaction 95: OH + CH3 (+M) <=> CH3OH (+M) */
    eqcon[94] *= 1e+06; 

    /*reaction 96: OH + CH3 <=> CH2 + H2O */
    /*eqcon[95] *= 1;  */

    /*reaction 97: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[96] *= 1;  */

    /*reaction 98: OH + CH4 <=> CH3 + H2O */
    /*eqcon[97] *= 1;  */

    /*reaction 99: OH + CO <=> H + CO2 */
    /*eqcon[98] *= 1;  */

    /*reaction 100: OH + HCO <=> H2O + CO */
    /*eqcon[99] *= 1;  */

    /*reaction 101: OH + CH2O <=> HCO + H2O */
    /*eqcon[100] *= 1;  */

    /*reaction 102: OH + CH2OH <=> H2O + CH2O */
    /*eqcon[101] *= 1;  */

    /*reaction 103: OH + CH3O <=> H2O + CH2O */
    /*eqcon[102] *= 1;  */

    /*reaction 104: OH + CH3OH <=> CH2OH + H2O */
    /*eqcon[103] *= 1;  */

    /*reaction 105: OH + CH3OH <=> CH3O + H2O */
    /*eqcon[104] *= 1;  */

    /*reaction 106: OH + C2H <=> H + HCCO */
    /*eqcon[105] *= 1;  */

    /*reaction 107: OH + C2H2 <=> H + CH2CO */
    /*eqcon[106] *= 1;  */

    /*reaction 108: OH + C2H2 <=> H + HCCOH */
    /*eqcon[107] *= 1;  */

    /*reaction 109: OH + C2H2 <=> C2H + H2O */
    /*eqcon[108] *= 1;  */

    /*reaction 110: OH + C2H2 <=> CH3 + CO */
    /*eqcon[109] *= 1;  */

    /*reaction 111: OH + C2H3 <=> H2O + C2H2 */
    /*eqcon[110] *= 1;  */

    /*reaction 112: OH + C2H4 <=> C2H3 + H2O */
    /*eqcon[111] *= 1;  */

    /*reaction 113: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[112] *= 1;  */

    /*reaction 114: OH + CH2CO <=> HCCO + H2O */
    /*eqcon[113] *= 1;  */

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[114] *= 1;  */

    /*reaction 116: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[115] *= 1;  */

    /*reaction 117: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[116] *= 1;  */

    /*reaction 118: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[117] *= 1;  */

    /*reaction 119: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[118] *= 1;  */

    /*reaction 120: HO2 + CO <=> OH + CO2 */
    /*eqcon[119] *= 1;  */

    /*reaction 121: HO2 + CH2O <=> HCO + H2O2 */
    /*eqcon[120] *= 1;  */

    /*reaction 122: C + O2 <=> O + CO */
    /*eqcon[121] *= 1;  */

    /*reaction 123: C + CH2 <=> H + C2H */
    /*eqcon[122] *= 1;  */

    /*reaction 124: C + CH3 <=> H + C2H2 */
    /*eqcon[123] *= 1;  */

    /*reaction 125: CH + O2 <=> O + HCO */
    /*eqcon[124] *= 1;  */

    /*reaction 126: CH + H2 <=> H + CH2 */
    /*eqcon[125] *= 1;  */

    /*reaction 127: CH + H2O <=> H + CH2O */
    /*eqcon[126] *= 1;  */

    /*reaction 128: CH + CH2 <=> H + C2H2 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: CH + CH3 <=> H + C2H3 */
    /*eqcon[128] *= 1;  */

    /*reaction 130: CH + CH4 <=> H + C2H4 */
    /*eqcon[129] *= 1;  */

    /*reaction 131: CH + CO (+M) <=> HCCO (+M) */
    eqcon[130] *= 1e+06; 

    /*reaction 132: CH + CO2 <=> HCO + CO */
    /*eqcon[131] *= 1;  */

    /*reaction 133: CH + CH2O <=> H + CH2CO */
    /*eqcon[132] *= 1;  */

    /*reaction 134: CH + HCCO <=> CO + C2H2 */
    /*eqcon[133] *= 1;  */

    /*reaction 135: CH2 + O2 <=> OH + HCO */
    /*eqcon[134] *= 1;  */

    /*reaction 136: CH2 + H2 <=> H + CH3 */
    /*eqcon[135] *= 1;  */

    /*reaction 137: 2 CH2 <=> H2 + C2H2 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[137] *= 1;  */

    /*reaction 139: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[138] *= 1;  */

    /*reaction 140: CH2 + CO (+M) <=> CH2CO (+M) */
    eqcon[139] *= 1e+06; 

    /*reaction 141: CH2 + HCCO <=> C2H3 + CO */
    /*eqcon[140] *= 1;  */

    /*reaction 142: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[141] *= 1;  */

    /*reaction 143: CH2(S) + AR <=> CH2 + AR */
    /*eqcon[142] *= 1;  */

    /*reaction 144: CH2(S) + O2 <=> H + OH + CO */
    eqcon[143] *= 1e-06; 

    /*reaction 145: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[144] *= 1;  */

    /*reaction 146: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[145] *= 1;  */

    /*reaction 147: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    eqcon[146] *= 1e+06; 

    /*reaction 148: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[147] *= 1;  */

    /*reaction 149: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[148] *= 1;  */

    /*reaction 150: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[149] *= 1;  */

    /*reaction 151: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[150] *= 1;  */

    /*reaction 152: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[152] *= 1;  */

    /*reaction 154: CH2(S) + C2H6 <=> CH3 + C2H5 */
    /*eqcon[153] *= 1;  */

    /*reaction 155: CH3 + O2 <=> O + CH3O */
    /*eqcon[154] *= 1;  */

    /*reaction 156: CH3 + O2 <=> OH + CH2O */
    /*eqcon[155] *= 1;  */

    /*reaction 157: CH3 + H2O2 <=> HO2 + CH4 */
    /*eqcon[156] *= 1;  */

    /*reaction 158: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[157] *= 1e+06; 

    /*reaction 159: 2 CH3 <=> H + C2H5 */
    /*eqcon[158] *= 1;  */

    /*reaction 160: CH3 + HCO <=> CH4 + CO */
    /*eqcon[159] *= 1;  */

    /*reaction 161: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[160] *= 1;  */

    /*reaction 162: CH3 + CH3OH <=> CH2OH + CH4 */
    /*eqcon[161] *= 1;  */

    /*reaction 163: CH3 + CH3OH <=> CH3O + CH4 */
    /*eqcon[162] *= 1;  */

    /*reaction 164: CH3 + C2H4 <=> C2H3 + CH4 */
    /*eqcon[163] *= 1;  */

    /*reaction 165: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[164] *= 1;  */

    /*reaction 166: HCO + H2O <=> H + CO + H2O */
    eqcon[165] *= 1e-06; 

    /*reaction 167: HCO + M <=> H + CO + M */
    eqcon[166] *= 1e-06; 

    /*reaction 168: HCO + O2 <=> HO2 + CO */
    /*eqcon[167] *= 1;  */

    /*reaction 169: CH2OH + O2 <=> HO2 + CH2O */
    /*eqcon[168] *= 1;  */

    /*reaction 170: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[169] *= 1;  */

    /*reaction 171: C2H + O2 <=> HCO + CO */
    /*eqcon[170] *= 1;  */

    /*reaction 172: C2H + H2 <=> H + C2H2 */
    /*eqcon[171] *= 1;  */

    /*reaction 173: C2H3 + O2 <=> HCO + CH2O */
    /*eqcon[172] *= 1;  */

    /*reaction 174: C2H4 (+M) <=> H2 + C2H2 (+M) */
    eqcon[173] *= 1e-06; 

    /*reaction 175: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[174] *= 1;  */

    /*reaction 176: HCCO + O2 <=> OH + 2 CO */
    eqcon[175] *= 1e-06; 

    /*reaction 177: 2 HCCO <=> 2 CO + C2H2 */
    eqcon[176] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[32]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 O + M <=> O2 + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: O + H + M <=> OH + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: O + H2 <=> H + OH */
    /*eqcon[2] *= 1;  */

    /*reaction 4: O + HO2 <=> OH + O2 */
    /*eqcon[3] *= 1;  */

    /*reaction 5: O + H2O2 <=> OH + HO2 */
    /*eqcon[4] *= 1;  */

    /*reaction 6: O + CH <=> H + CO */
    /*eqcon[5] *= 1;  */

    /*reaction 7: O + CH2 <=> H + HCO */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + CH2(S) <=> H2 + CO */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + CH2(S) <=> H + HCO */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + CH3 <=> H + CH2O */
    /*eqcon[9] *= 1;  */

    /*reaction 11: O + CH4 <=> OH + CH3 */
    /*eqcon[10] *= 1;  */

    /*reaction 12: O + CO + M <=> CO2 + M */
    eqcon[11] *= 1e+06; 

    /*reaction 13: O + HCO <=> OH + CO */
    /*eqcon[12] *= 1;  */

    /*reaction 14: O + HCO <=> H + CO2 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O + CH2O <=> OH + HCO */
    /*eqcon[14] *= 1;  */

    /*reaction 16: O + CH2OH <=> OH + CH2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: O + CH3O <=> OH + CH2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: O + CH3OH <=> OH + CH2OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: O + CH3OH <=> OH + CH3O */
    /*eqcon[18] *= 1;  */

    /*reaction 20: O + C2H <=> CH + CO */
    /*eqcon[19] *= 1;  */

    /*reaction 21: O + C2H2 <=> H + HCCO */
    /*eqcon[20] *= 1;  */

    /*reaction 22: O + C2H2 <=> OH + C2H */
    /*eqcon[21] *= 1;  */

    /*reaction 23: O + C2H2 <=> CO + CH2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: O + C2H3 <=> H + CH2CO */
    /*eqcon[23] *= 1;  */

    /*reaction 25: O + C2H4 <=> CH3 + HCO */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O + C2H5 <=> CH3 + CH2O */
    /*eqcon[25] *= 1;  */

    /*reaction 27: O + C2H6 <=> OH + C2H5 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: O + HCCO <=> H + 2 CO */
    eqcon[27] *= 1e-06; 

    /*reaction 29: O + CH2CO <=> OH + HCCO */
    /*eqcon[28] *= 1;  */

    /*reaction 30: O + CH2CO <=> CH2 + CO2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: O2 + CO <=> O + CO2 */
    /*eqcon[30] *= 1;  */

    /*reaction 32: O2 + CH2O <=> HO2 + HCO */
    /*eqcon[31] *= 1;  */

    /*reaction 33: H + O2 + M <=> HO2 + M */
    eqcon[32] *= 1e+06; 

    /*reaction 34: H + 2 O2 <=> HO2 + O2 */
    eqcon[33] *= 1e+06; 

    /*reaction 35: H + O2 + H2O <=> HO2 + H2O */
    eqcon[34] *= 1e+06; 

    /*reaction 36: H + O2 + N2 <=> HO2 + N2 */
    eqcon[35] *= 1e+06; 

    /*reaction 37: H + O2 + AR <=> HO2 + AR */
    eqcon[36] *= 1e+06; 

    /*reaction 38: H + O2 <=> O + OH */
    /*eqcon[37] *= 1;  */

    /*reaction 39: 2 H + M <=> H2 + M */
    eqcon[38] *= 1e+06; 

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    eqcon[39] *= 1e+06; 

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    eqcon[40] *= 1e+06; 

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    eqcon[41] *= 1e+06; 

    /*reaction 43: H + OH + M <=> H2O + M */
    eqcon[42] *= 1e+06; 

    /*reaction 44: H + HO2 <=> O + H2O */
    /*eqcon[43] *= 1;  */

    /*reaction 45: H + HO2 <=> O2 + H2 */
    /*eqcon[44] *= 1;  */

    /*reaction 46: H + HO2 <=> 2 OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + H2O2 <=> HO2 + H2 */
    /*eqcon[46] *= 1;  */

    /*reaction 48: H + H2O2 <=> OH + H2O */
    /*eqcon[47] *= 1;  */

    /*reaction 49: H + CH <=> C + H2 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[49] *= 1e+06; 

    /*reaction 51: H + CH2(S) <=> CH + H2 */
    /*eqcon[50] *= 1;  */

    /*reaction 52: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[51] *= 1e+06; 

    /*reaction 53: H + CH4 <=> CH3 + H2 */
    /*eqcon[52] *= 1;  */

    /*reaction 54: H + HCO (+M) <=> CH2O (+M) */
    eqcon[53] *= 1e+06; 

    /*reaction 55: H + HCO <=> H2 + CO */
    /*eqcon[54] *= 1;  */

    /*reaction 56: H + CH2O (+M) <=> CH2OH (+M) */
    eqcon[55] *= 1e+06; 

    /*reaction 57: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[56] *= 1e+06; 

    /*reaction 58: H + CH2O <=> HCO + H2 */
    /*eqcon[57] *= 1;  */

    /*reaction 59: H + CH2OH (+M) <=> CH3OH (+M) */
    eqcon[58] *= 1e+06; 

    /*reaction 60: H + CH2OH <=> H2 + CH2O */
    /*eqcon[59] *= 1;  */

    /*reaction 61: H + CH2OH <=> OH + CH3 */
    /*eqcon[60] *= 1;  */

    /*reaction 62: H + CH2OH <=> CH2(S) + H2O */
    /*eqcon[61] *= 1;  */

    /*reaction 63: H + CH3O (+M) <=> CH3OH (+M) */
    eqcon[62] *= 1e+06; 

    /*reaction 64: H + CH3O <=> H + CH2OH */
    /*eqcon[63] *= 1;  */

    /*reaction 65: H + CH3O <=> H2 + CH2O */
    /*eqcon[64] *= 1;  */

    /*reaction 66: H + CH3O <=> OH + CH3 */
    /*eqcon[65] *= 1;  */

    /*reaction 67: H + CH3O <=> CH2(S) + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: H + CH3OH <=> CH2OH + H2 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: H + CH3OH <=> CH3O + H2 */
    /*eqcon[68] *= 1;  */

    /*reaction 70: H + C2H (+M) <=> C2H2 (+M) */
    eqcon[69] *= 1e+06; 

    /*reaction 71: H + C2H2 (+M) <=> C2H3 (+M) */
    eqcon[70] *= 1e+06; 

    /*reaction 72: H + C2H3 (+M) <=> C2H4 (+M) */
    eqcon[71] *= 1e+06; 

    /*reaction 73: H + C2H3 <=> H2 + C2H2 */
    /*eqcon[72] *= 1;  */

    /*reaction 74: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[73] *= 1e+06; 

    /*reaction 75: H + C2H4 <=> C2H3 + H2 */
    /*eqcon[74] *= 1;  */

    /*reaction 76: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[75] *= 1e+06; 

    /*reaction 77: H + C2H5 <=> H2 + C2H4 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[77] *= 1;  */

    /*reaction 79: H + HCCO <=> CH2(S) + CO */
    /*eqcon[78] *= 1;  */

    /*reaction 80: H + CH2CO <=> HCCO + H2 */
    /*eqcon[79] *= 1;  */

    /*reaction 81: H + CH2CO <=> CH3 + CO */
    /*eqcon[80] *= 1;  */

    /*reaction 82: H + HCCOH <=> H + CH2CO */
    /*eqcon[81] *= 1;  */

    /*reaction 83: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[82] *= 1e+06; 

    /*reaction 84: OH + H2 <=> H + H2O */
    /*eqcon[83] *= 1;  */

    /*reaction 85: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[84] *= 1e+06; 

    /*reaction 86: 2 OH <=> O + H2O */
    /*eqcon[85] *= 1;  */

    /*reaction 87: OH + HO2 <=> O2 + H2O */
    /*eqcon[86] *= 1;  */

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[87] *= 1;  */

    /*reaction 89: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[88] *= 1;  */

    /*reaction 90: OH + C <=> H + CO */
    /*eqcon[89] *= 1;  */

    /*reaction 91: OH + CH <=> H + HCO */
    /*eqcon[90] *= 1;  */

    /*reaction 92: OH + CH2 <=> H + CH2O */
    /*eqcon[91] *= 1;  */

    /*reaction 93: OH + CH2 <=> CH + H2O */
    /*eqcon[92] *= 1;  */

    /*reaction 94: OH + CH2(S) <=> H + CH2O */
    /*eqcon[93] *= 1;  */

    /*reaction 95: OH + CH3 (+M) <=> CH3OH (+M) */
    eqcon[94] *= 1e+06; 

    /*reaction 96: OH + CH3 <=> CH2 + H2O */
    /*eqcon[95] *= 1;  */

    /*reaction 97: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[96] *= 1;  */

    /*reaction 98: OH + CH4 <=> CH3 + H2O */
    /*eqcon[97] *= 1;  */

    /*reaction 99: OH + CO <=> H + CO2 */
    /*eqcon[98] *= 1;  */

    /*reaction 100: OH + HCO <=> H2O + CO */
    /*eqcon[99] *= 1;  */

    /*reaction 101: OH + CH2O <=> HCO + H2O */
    /*eqcon[100] *= 1;  */

    /*reaction 102: OH + CH2OH <=> H2O + CH2O */
    /*eqcon[101] *= 1;  */

    /*reaction 103: OH + CH3O <=> H2O + CH2O */
    /*eqcon[102] *= 1;  */

    /*reaction 104: OH + CH3OH <=> CH2OH + H2O */
    /*eqcon[103] *= 1;  */

    /*reaction 105: OH + CH3OH <=> CH3O + H2O */
    /*eqcon[104] *= 1;  */

    /*reaction 106: OH + C2H <=> H + HCCO */
    /*eqcon[105] *= 1;  */

    /*reaction 107: OH + C2H2 <=> H + CH2CO */
    /*eqcon[106] *= 1;  */

    /*reaction 108: OH + C2H2 <=> H + HCCOH */
    /*eqcon[107] *= 1;  */

    /*reaction 109: OH + C2H2 <=> C2H + H2O */
    /*eqcon[108] *= 1;  */

    /*reaction 110: OH + C2H2 <=> CH3 + CO */
    /*eqcon[109] *= 1;  */

    /*reaction 111: OH + C2H3 <=> H2O + C2H2 */
    /*eqcon[110] *= 1;  */

    /*reaction 112: OH + C2H4 <=> C2H3 + H2O */
    /*eqcon[111] *= 1;  */

    /*reaction 113: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[112] *= 1;  */

    /*reaction 114: OH + CH2CO <=> HCCO + H2O */
    /*eqcon[113] *= 1;  */

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[114] *= 1;  */

    /*reaction 116: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[115] *= 1;  */

    /*reaction 117: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[116] *= 1;  */

    /*reaction 118: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[117] *= 1;  */

    /*reaction 119: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[118] *= 1;  */

    /*reaction 120: HO2 + CO <=> OH + CO2 */
    /*eqcon[119] *= 1;  */

    /*reaction 121: HO2 + CH2O <=> HCO + H2O2 */
    /*eqcon[120] *= 1;  */

    /*reaction 122: C + O2 <=> O + CO */
    /*eqcon[121] *= 1;  */

    /*reaction 123: C + CH2 <=> H + C2H */
    /*eqcon[122] *= 1;  */

    /*reaction 124: C + CH3 <=> H + C2H2 */
    /*eqcon[123] *= 1;  */

    /*reaction 125: CH + O2 <=> O + HCO */
    /*eqcon[124] *= 1;  */

    /*reaction 126: CH + H2 <=> H + CH2 */
    /*eqcon[125] *= 1;  */

    /*reaction 127: CH + H2O <=> H + CH2O */
    /*eqcon[126] *= 1;  */

    /*reaction 128: CH + CH2 <=> H + C2H2 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: CH + CH3 <=> H + C2H3 */
    /*eqcon[128] *= 1;  */

    /*reaction 130: CH + CH4 <=> H + C2H4 */
    /*eqcon[129] *= 1;  */

    /*reaction 131: CH + CO (+M) <=> HCCO (+M) */
    eqcon[130] *= 1e+06; 

    /*reaction 132: CH + CO2 <=> HCO + CO */
    /*eqcon[131] *= 1;  */

    /*reaction 133: CH + CH2O <=> H + CH2CO */
    /*eqcon[132] *= 1;  */

    /*reaction 134: CH + HCCO <=> CO + C2H2 */
    /*eqcon[133] *= 1;  */

    /*reaction 135: CH2 + O2 <=> OH + HCO */
    /*eqcon[134] *= 1;  */

    /*reaction 136: CH2 + H2 <=> H + CH3 */
    /*eqcon[135] *= 1;  */

    /*reaction 137: 2 CH2 <=> H2 + C2H2 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[137] *= 1;  */

    /*reaction 139: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[138] *= 1;  */

    /*reaction 140: CH2 + CO (+M) <=> CH2CO (+M) */
    eqcon[139] *= 1e+06; 

    /*reaction 141: CH2 + HCCO <=> C2H3 + CO */
    /*eqcon[140] *= 1;  */

    /*reaction 142: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[141] *= 1;  */

    /*reaction 143: CH2(S) + AR <=> CH2 + AR */
    /*eqcon[142] *= 1;  */

    /*reaction 144: CH2(S) + O2 <=> H + OH + CO */
    eqcon[143] *= 1e-06; 

    /*reaction 145: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[144] *= 1;  */

    /*reaction 146: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[145] *= 1;  */

    /*reaction 147: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    eqcon[146] *= 1e+06; 

    /*reaction 148: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[147] *= 1;  */

    /*reaction 149: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[148] *= 1;  */

    /*reaction 150: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[149] *= 1;  */

    /*reaction 151: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[150] *= 1;  */

    /*reaction 152: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[152] *= 1;  */

    /*reaction 154: CH2(S) + C2H6 <=> CH3 + C2H5 */
    /*eqcon[153] *= 1;  */

    /*reaction 155: CH3 + O2 <=> O + CH3O */
    /*eqcon[154] *= 1;  */

    /*reaction 156: CH3 + O2 <=> OH + CH2O */
    /*eqcon[155] *= 1;  */

    /*reaction 157: CH3 + H2O2 <=> HO2 + CH4 */
    /*eqcon[156] *= 1;  */

    /*reaction 158: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[157] *= 1e+06; 

    /*reaction 159: 2 CH3 <=> H + C2H5 */
    /*eqcon[158] *= 1;  */

    /*reaction 160: CH3 + HCO <=> CH4 + CO */
    /*eqcon[159] *= 1;  */

    /*reaction 161: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[160] *= 1;  */

    /*reaction 162: CH3 + CH3OH <=> CH2OH + CH4 */
    /*eqcon[161] *= 1;  */

    /*reaction 163: CH3 + CH3OH <=> CH3O + CH4 */
    /*eqcon[162] *= 1;  */

    /*reaction 164: CH3 + C2H4 <=> C2H3 + CH4 */
    /*eqcon[163] *= 1;  */

    /*reaction 165: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[164] *= 1;  */

    /*reaction 166: HCO + H2O <=> H + CO + H2O */
    eqcon[165] *= 1e-06; 

    /*reaction 167: HCO + M <=> H + CO + M */
    eqcon[166] *= 1e-06; 

    /*reaction 168: HCO + O2 <=> HO2 + CO */
    /*eqcon[167] *= 1;  */

    /*reaction 169: CH2OH + O2 <=> HO2 + CH2O */
    /*eqcon[168] *= 1;  */

    /*reaction 170: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[169] *= 1;  */

    /*reaction 171: C2H + O2 <=> HCO + CO */
    /*eqcon[170] *= 1;  */

    /*reaction 172: C2H + H2 <=> H + C2H2 */
    /*eqcon[171] *= 1;  */

    /*reaction 173: C2H3 + O2 <=> HCO + CH2O */
    /*eqcon[172] *= 1;  */

    /*reaction 174: C2H4 (+M) <=> H2 + C2H2 (+M) */
    eqcon[173] *= 1e-06; 

    /*reaction 175: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[174] *= 1;  */

    /*reaction 176: HCCO + O2 <=> OH + 2 CO */
    eqcon[175] *= 1e-06; 

    /*reaction 177: 2 HCCO <=> 2 CO + C2H2 */
    eqcon[176] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[32]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 O + M <=> O2 + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: O + H + M <=> OH + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: O + H2 <=> H + OH */
    /*eqcon[2] *= 1;  */

    /*reaction 4: O + HO2 <=> OH + O2 */
    /*eqcon[3] *= 1;  */

    /*reaction 5: O + H2O2 <=> OH + HO2 */
    /*eqcon[4] *= 1;  */

    /*reaction 6: O + CH <=> H + CO */
    /*eqcon[5] *= 1;  */

    /*reaction 7: O + CH2 <=> H + HCO */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + CH2(S) <=> H2 + CO */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + CH2(S) <=> H + HCO */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + CH3 <=> H + CH2O */
    /*eqcon[9] *= 1;  */

    /*reaction 11: O + CH4 <=> OH + CH3 */
    /*eqcon[10] *= 1;  */

    /*reaction 12: O + CO + M <=> CO2 + M */
    eqcon[11] *= 1e+06; 

    /*reaction 13: O + HCO <=> OH + CO */
    /*eqcon[12] *= 1;  */

    /*reaction 14: O + HCO <=> H + CO2 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O + CH2O <=> OH + HCO */
    /*eqcon[14] *= 1;  */

    /*reaction 16: O + CH2OH <=> OH + CH2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: O + CH3O <=> OH + CH2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: O + CH3OH <=> OH + CH2OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: O + CH3OH <=> OH + CH3O */
    /*eqcon[18] *= 1;  */

    /*reaction 20: O + C2H <=> CH + CO */
    /*eqcon[19] *= 1;  */

    /*reaction 21: O + C2H2 <=> H + HCCO */
    /*eqcon[20] *= 1;  */

    /*reaction 22: O + C2H2 <=> OH + C2H */
    /*eqcon[21] *= 1;  */

    /*reaction 23: O + C2H2 <=> CO + CH2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: O + C2H3 <=> H + CH2CO */
    /*eqcon[23] *= 1;  */

    /*reaction 25: O + C2H4 <=> CH3 + HCO */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O + C2H5 <=> CH3 + CH2O */
    /*eqcon[25] *= 1;  */

    /*reaction 27: O + C2H6 <=> OH + C2H5 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: O + HCCO <=> H + 2 CO */
    eqcon[27] *= 1e-06; 

    /*reaction 29: O + CH2CO <=> OH + HCCO */
    /*eqcon[28] *= 1;  */

    /*reaction 30: O + CH2CO <=> CH2 + CO2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: O2 + CO <=> O + CO2 */
    /*eqcon[30] *= 1;  */

    /*reaction 32: O2 + CH2O <=> HO2 + HCO */
    /*eqcon[31] *= 1;  */

    /*reaction 33: H + O2 + M <=> HO2 + M */
    eqcon[32] *= 1e+06; 

    /*reaction 34: H + 2 O2 <=> HO2 + O2 */
    eqcon[33] *= 1e+06; 

    /*reaction 35: H + O2 + H2O <=> HO2 + H2O */
    eqcon[34] *= 1e+06; 

    /*reaction 36: H + O2 + N2 <=> HO2 + N2 */
    eqcon[35] *= 1e+06; 

    /*reaction 37: H + O2 + AR <=> HO2 + AR */
    eqcon[36] *= 1e+06; 

    /*reaction 38: H + O2 <=> O + OH */
    /*eqcon[37] *= 1;  */

    /*reaction 39: 2 H + M <=> H2 + M */
    eqcon[38] *= 1e+06; 

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    eqcon[39] *= 1e+06; 

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    eqcon[40] *= 1e+06; 

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    eqcon[41] *= 1e+06; 

    /*reaction 43: H + OH + M <=> H2O + M */
    eqcon[42] *= 1e+06; 

    /*reaction 44: H + HO2 <=> O + H2O */
    /*eqcon[43] *= 1;  */

    /*reaction 45: H + HO2 <=> O2 + H2 */
    /*eqcon[44] *= 1;  */

    /*reaction 46: H + HO2 <=> 2 OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + H2O2 <=> HO2 + H2 */
    /*eqcon[46] *= 1;  */

    /*reaction 48: H + H2O2 <=> OH + H2O */
    /*eqcon[47] *= 1;  */

    /*reaction 49: H + CH <=> C + H2 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[49] *= 1e+06; 

    /*reaction 51: H + CH2(S) <=> CH + H2 */
    /*eqcon[50] *= 1;  */

    /*reaction 52: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[51] *= 1e+06; 

    /*reaction 53: H + CH4 <=> CH3 + H2 */
    /*eqcon[52] *= 1;  */

    /*reaction 54: H + HCO (+M) <=> CH2O (+M) */
    eqcon[53] *= 1e+06; 

    /*reaction 55: H + HCO <=> H2 + CO */
    /*eqcon[54] *= 1;  */

    /*reaction 56: H + CH2O (+M) <=> CH2OH (+M) */
    eqcon[55] *= 1e+06; 

    /*reaction 57: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[56] *= 1e+06; 

    /*reaction 58: H + CH2O <=> HCO + H2 */
    /*eqcon[57] *= 1;  */

    /*reaction 59: H + CH2OH (+M) <=> CH3OH (+M) */
    eqcon[58] *= 1e+06; 

    /*reaction 60: H + CH2OH <=> H2 + CH2O */
    /*eqcon[59] *= 1;  */

    /*reaction 61: H + CH2OH <=> OH + CH3 */
    /*eqcon[60] *= 1;  */

    /*reaction 62: H + CH2OH <=> CH2(S) + H2O */
    /*eqcon[61] *= 1;  */

    /*reaction 63: H + CH3O (+M) <=> CH3OH (+M) */
    eqcon[62] *= 1e+06; 

    /*reaction 64: H + CH3O <=> H + CH2OH */
    /*eqcon[63] *= 1;  */

    /*reaction 65: H + CH3O <=> H2 + CH2O */
    /*eqcon[64] *= 1;  */

    /*reaction 66: H + CH3O <=> OH + CH3 */
    /*eqcon[65] *= 1;  */

    /*reaction 67: H + CH3O <=> CH2(S) + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: H + CH3OH <=> CH2OH + H2 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: H + CH3OH <=> CH3O + H2 */
    /*eqcon[68] *= 1;  */

    /*reaction 70: H + C2H (+M) <=> C2H2 (+M) */
    eqcon[69] *= 1e+06; 

    /*reaction 71: H + C2H2 (+M) <=> C2H3 (+M) */
    eqcon[70] *= 1e+06; 

    /*reaction 72: H + C2H3 (+M) <=> C2H4 (+M) */
    eqcon[71] *= 1e+06; 

    /*reaction 73: H + C2H3 <=> H2 + C2H2 */
    /*eqcon[72] *= 1;  */

    /*reaction 74: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[73] *= 1e+06; 

    /*reaction 75: H + C2H4 <=> C2H3 + H2 */
    /*eqcon[74] *= 1;  */

    /*reaction 76: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[75] *= 1e+06; 

    /*reaction 77: H + C2H5 <=> H2 + C2H4 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[77] *= 1;  */

    /*reaction 79: H + HCCO <=> CH2(S) + CO */
    /*eqcon[78] *= 1;  */

    /*reaction 80: H + CH2CO <=> HCCO + H2 */
    /*eqcon[79] *= 1;  */

    /*reaction 81: H + CH2CO <=> CH3 + CO */
    /*eqcon[80] *= 1;  */

    /*reaction 82: H + HCCOH <=> H + CH2CO */
    /*eqcon[81] *= 1;  */

    /*reaction 83: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[82] *= 1e+06; 

    /*reaction 84: OH + H2 <=> H + H2O */
    /*eqcon[83] *= 1;  */

    /*reaction 85: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[84] *= 1e+06; 

    /*reaction 86: 2 OH <=> O + H2O */
    /*eqcon[85] *= 1;  */

    /*reaction 87: OH + HO2 <=> O2 + H2O */
    /*eqcon[86] *= 1;  */

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[87] *= 1;  */

    /*reaction 89: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[88] *= 1;  */

    /*reaction 90: OH + C <=> H + CO */
    /*eqcon[89] *= 1;  */

    /*reaction 91: OH + CH <=> H + HCO */
    /*eqcon[90] *= 1;  */

    /*reaction 92: OH + CH2 <=> H + CH2O */
    /*eqcon[91] *= 1;  */

    /*reaction 93: OH + CH2 <=> CH + H2O */
    /*eqcon[92] *= 1;  */

    /*reaction 94: OH + CH2(S) <=> H + CH2O */
    /*eqcon[93] *= 1;  */

    /*reaction 95: OH + CH3 (+M) <=> CH3OH (+M) */
    eqcon[94] *= 1e+06; 

    /*reaction 96: OH + CH3 <=> CH2 + H2O */
    /*eqcon[95] *= 1;  */

    /*reaction 97: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[96] *= 1;  */

    /*reaction 98: OH + CH4 <=> CH3 + H2O */
    /*eqcon[97] *= 1;  */

    /*reaction 99: OH + CO <=> H + CO2 */
    /*eqcon[98] *= 1;  */

    /*reaction 100: OH + HCO <=> H2O + CO */
    /*eqcon[99] *= 1;  */

    /*reaction 101: OH + CH2O <=> HCO + H2O */
    /*eqcon[100] *= 1;  */

    /*reaction 102: OH + CH2OH <=> H2O + CH2O */
    /*eqcon[101] *= 1;  */

    /*reaction 103: OH + CH3O <=> H2O + CH2O */
    /*eqcon[102] *= 1;  */

    /*reaction 104: OH + CH3OH <=> CH2OH + H2O */
    /*eqcon[103] *= 1;  */

    /*reaction 105: OH + CH3OH <=> CH3O + H2O */
    /*eqcon[104] *= 1;  */

    /*reaction 106: OH + C2H <=> H + HCCO */
    /*eqcon[105] *= 1;  */

    /*reaction 107: OH + C2H2 <=> H + CH2CO */
    /*eqcon[106] *= 1;  */

    /*reaction 108: OH + C2H2 <=> H + HCCOH */
    /*eqcon[107] *= 1;  */

    /*reaction 109: OH + C2H2 <=> C2H + H2O */
    /*eqcon[108] *= 1;  */

    /*reaction 110: OH + C2H2 <=> CH3 + CO */
    /*eqcon[109] *= 1;  */

    /*reaction 111: OH + C2H3 <=> H2O + C2H2 */
    /*eqcon[110] *= 1;  */

    /*reaction 112: OH + C2H4 <=> C2H3 + H2O */
    /*eqcon[111] *= 1;  */

    /*reaction 113: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[112] *= 1;  */

    /*reaction 114: OH + CH2CO <=> HCCO + H2O */
    /*eqcon[113] *= 1;  */

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[114] *= 1;  */

    /*reaction 116: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[115] *= 1;  */

    /*reaction 117: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[116] *= 1;  */

    /*reaction 118: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[117] *= 1;  */

    /*reaction 119: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[118] *= 1;  */

    /*reaction 120: HO2 + CO <=> OH + CO2 */
    /*eqcon[119] *= 1;  */

    /*reaction 121: HO2 + CH2O <=> HCO + H2O2 */
    /*eqcon[120] *= 1;  */

    /*reaction 122: C + O2 <=> O + CO */
    /*eqcon[121] *= 1;  */

    /*reaction 123: C + CH2 <=> H + C2H */
    /*eqcon[122] *= 1;  */

    /*reaction 124: C + CH3 <=> H + C2H2 */
    /*eqcon[123] *= 1;  */

    /*reaction 125: CH + O2 <=> O + HCO */
    /*eqcon[124] *= 1;  */

    /*reaction 126: CH + H2 <=> H + CH2 */
    /*eqcon[125] *= 1;  */

    /*reaction 127: CH + H2O <=> H + CH2O */
    /*eqcon[126] *= 1;  */

    /*reaction 128: CH + CH2 <=> H + C2H2 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: CH + CH3 <=> H + C2H3 */
    /*eqcon[128] *= 1;  */

    /*reaction 130: CH + CH4 <=> H + C2H4 */
    /*eqcon[129] *= 1;  */

    /*reaction 131: CH + CO (+M) <=> HCCO (+M) */
    eqcon[130] *= 1e+06; 

    /*reaction 132: CH + CO2 <=> HCO + CO */
    /*eqcon[131] *= 1;  */

    /*reaction 133: CH + CH2O <=> H + CH2CO */
    /*eqcon[132] *= 1;  */

    /*reaction 134: CH + HCCO <=> CO + C2H2 */
    /*eqcon[133] *= 1;  */

    /*reaction 135: CH2 + O2 <=> OH + HCO */
    /*eqcon[134] *= 1;  */

    /*reaction 136: CH2 + H2 <=> H + CH3 */
    /*eqcon[135] *= 1;  */

    /*reaction 137: 2 CH2 <=> H2 + C2H2 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[137] *= 1;  */

    /*reaction 139: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[138] *= 1;  */

    /*reaction 140: CH2 + CO (+M) <=> CH2CO (+M) */
    eqcon[139] *= 1e+06; 

    /*reaction 141: CH2 + HCCO <=> C2H3 + CO */
    /*eqcon[140] *= 1;  */

    /*reaction 142: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[141] *= 1;  */

    /*reaction 143: CH2(S) + AR <=> CH2 + AR */
    /*eqcon[142] *= 1;  */

    /*reaction 144: CH2(S) + O2 <=> H + OH + CO */
    eqcon[143] *= 1e-06; 

    /*reaction 145: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[144] *= 1;  */

    /*reaction 146: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[145] *= 1;  */

    /*reaction 147: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    eqcon[146] *= 1e+06; 

    /*reaction 148: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[147] *= 1;  */

    /*reaction 149: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[148] *= 1;  */

    /*reaction 150: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[149] *= 1;  */

    /*reaction 151: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[150] *= 1;  */

    /*reaction 152: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[152] *= 1;  */

    /*reaction 154: CH2(S) + C2H6 <=> CH3 + C2H5 */
    /*eqcon[153] *= 1;  */

    /*reaction 155: CH3 + O2 <=> O + CH3O */
    /*eqcon[154] *= 1;  */

    /*reaction 156: CH3 + O2 <=> OH + CH2O */
    /*eqcon[155] *= 1;  */

    /*reaction 157: CH3 + H2O2 <=> HO2 + CH4 */
    /*eqcon[156] *= 1;  */

    /*reaction 158: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[157] *= 1e+06; 

    /*reaction 159: 2 CH3 <=> H + C2H5 */
    /*eqcon[158] *= 1;  */

    /*reaction 160: CH3 + HCO <=> CH4 + CO */
    /*eqcon[159] *= 1;  */

    /*reaction 161: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[160] *= 1;  */

    /*reaction 162: CH3 + CH3OH <=> CH2OH + CH4 */
    /*eqcon[161] *= 1;  */

    /*reaction 163: CH3 + CH3OH <=> CH3O + CH4 */
    /*eqcon[162] *= 1;  */

    /*reaction 164: CH3 + C2H4 <=> C2H3 + CH4 */
    /*eqcon[163] *= 1;  */

    /*reaction 165: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[164] *= 1;  */

    /*reaction 166: HCO + H2O <=> H + CO + H2O */
    eqcon[165] *= 1e-06; 

    /*reaction 167: HCO + M <=> H + CO + M */
    eqcon[166] *= 1e-06; 

    /*reaction 168: HCO + O2 <=> HO2 + CO */
    /*eqcon[167] *= 1;  */

    /*reaction 169: CH2OH + O2 <=> HO2 + CH2O */
    /*eqcon[168] *= 1;  */

    /*reaction 170: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[169] *= 1;  */

    /*reaction 171: C2H + O2 <=> HCO + CO */
    /*eqcon[170] *= 1;  */

    /*reaction 172: C2H + H2 <=> H + C2H2 */
    /*eqcon[171] *= 1;  */

    /*reaction 173: C2H3 + O2 <=> HCO + CH2O */
    /*eqcon[172] *= 1;  */

    /*reaction 174: C2H4 (+M) <=> H2 + C2H2 (+M) */
    eqcon[173] *= 1e-06; 

    /*reaction 175: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[174] *= 1;  */

    /*reaction 176: HCCO + O2 <=> OH + 2 CO */
    eqcon[175] *= 1e-06; 

    /*reaction 177: 2 HCCO <=> 2 CO + C2H2 */
    eqcon[176] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[32]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 O + M <=> O2 + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: O + H + M <=> OH + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: O + H2 <=> H + OH */
    /*eqcon[2] *= 1;  */

    /*reaction 4: O + HO2 <=> OH + O2 */
    /*eqcon[3] *= 1;  */

    /*reaction 5: O + H2O2 <=> OH + HO2 */
    /*eqcon[4] *= 1;  */

    /*reaction 6: O + CH <=> H + CO */
    /*eqcon[5] *= 1;  */

    /*reaction 7: O + CH2 <=> H + HCO */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + CH2(S) <=> H2 + CO */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + CH2(S) <=> H + HCO */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + CH3 <=> H + CH2O */
    /*eqcon[9] *= 1;  */

    /*reaction 11: O + CH4 <=> OH + CH3 */
    /*eqcon[10] *= 1;  */

    /*reaction 12: O + CO + M <=> CO2 + M */
    eqcon[11] *= 1e+06; 

    /*reaction 13: O + HCO <=> OH + CO */
    /*eqcon[12] *= 1;  */

    /*reaction 14: O + HCO <=> H + CO2 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O + CH2O <=> OH + HCO */
    /*eqcon[14] *= 1;  */

    /*reaction 16: O + CH2OH <=> OH + CH2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: O + CH3O <=> OH + CH2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: O + CH3OH <=> OH + CH2OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: O + CH3OH <=> OH + CH3O */
    /*eqcon[18] *= 1;  */

    /*reaction 20: O + C2H <=> CH + CO */
    /*eqcon[19] *= 1;  */

    /*reaction 21: O + C2H2 <=> H + HCCO */
    /*eqcon[20] *= 1;  */

    /*reaction 22: O + C2H2 <=> OH + C2H */
    /*eqcon[21] *= 1;  */

    /*reaction 23: O + C2H2 <=> CO + CH2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: O + C2H3 <=> H + CH2CO */
    /*eqcon[23] *= 1;  */

    /*reaction 25: O + C2H4 <=> CH3 + HCO */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O + C2H5 <=> CH3 + CH2O */
    /*eqcon[25] *= 1;  */

    /*reaction 27: O + C2H6 <=> OH + C2H5 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: O + HCCO <=> H + 2 CO */
    eqcon[27] *= 1e-06; 

    /*reaction 29: O + CH2CO <=> OH + HCCO */
    /*eqcon[28] *= 1;  */

    /*reaction 30: O + CH2CO <=> CH2 + CO2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: O2 + CO <=> O + CO2 */
    /*eqcon[30] *= 1;  */

    /*reaction 32: O2 + CH2O <=> HO2 + HCO */
    /*eqcon[31] *= 1;  */

    /*reaction 33: H + O2 + M <=> HO2 + M */
    eqcon[32] *= 1e+06; 

    /*reaction 34: H + 2 O2 <=> HO2 + O2 */
    eqcon[33] *= 1e+06; 

    /*reaction 35: H + O2 + H2O <=> HO2 + H2O */
    eqcon[34] *= 1e+06; 

    /*reaction 36: H + O2 + N2 <=> HO2 + N2 */
    eqcon[35] *= 1e+06; 

    /*reaction 37: H + O2 + AR <=> HO2 + AR */
    eqcon[36] *= 1e+06; 

    /*reaction 38: H + O2 <=> O + OH */
    /*eqcon[37] *= 1;  */

    /*reaction 39: 2 H + M <=> H2 + M */
    eqcon[38] *= 1e+06; 

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    eqcon[39] *= 1e+06; 

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    eqcon[40] *= 1e+06; 

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    eqcon[41] *= 1e+06; 

    /*reaction 43: H + OH + M <=> H2O + M */
    eqcon[42] *= 1e+06; 

    /*reaction 44: H + HO2 <=> O + H2O */
    /*eqcon[43] *= 1;  */

    /*reaction 45: H + HO2 <=> O2 + H2 */
    /*eqcon[44] *= 1;  */

    /*reaction 46: H + HO2 <=> 2 OH */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + H2O2 <=> HO2 + H2 */
    /*eqcon[46] *= 1;  */

    /*reaction 48: H + H2O2 <=> OH + H2O */
    /*eqcon[47] *= 1;  */

    /*reaction 49: H + CH <=> C + H2 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[49] *= 1e+06; 

    /*reaction 51: H + CH2(S) <=> CH + H2 */
    /*eqcon[50] *= 1;  */

    /*reaction 52: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[51] *= 1e+06; 

    /*reaction 53: H + CH4 <=> CH3 + H2 */
    /*eqcon[52] *= 1;  */

    /*reaction 54: H + HCO (+M) <=> CH2O (+M) */
    eqcon[53] *= 1e+06; 

    /*reaction 55: H + HCO <=> H2 + CO */
    /*eqcon[54] *= 1;  */

    /*reaction 56: H + CH2O (+M) <=> CH2OH (+M) */
    eqcon[55] *= 1e+06; 

    /*reaction 57: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[56] *= 1e+06; 

    /*reaction 58: H + CH2O <=> HCO + H2 */
    /*eqcon[57] *= 1;  */

    /*reaction 59: H + CH2OH (+M) <=> CH3OH (+M) */
    eqcon[58] *= 1e+06; 

    /*reaction 60: H + CH2OH <=> H2 + CH2O */
    /*eqcon[59] *= 1;  */

    /*reaction 61: H + CH2OH <=> OH + CH3 */
    /*eqcon[60] *= 1;  */

    /*reaction 62: H + CH2OH <=> CH2(S) + H2O */
    /*eqcon[61] *= 1;  */

    /*reaction 63: H + CH3O (+M) <=> CH3OH (+M) */
    eqcon[62] *= 1e+06; 

    /*reaction 64: H + CH3O <=> H + CH2OH */
    /*eqcon[63] *= 1;  */

    /*reaction 65: H + CH3O <=> H2 + CH2O */
    /*eqcon[64] *= 1;  */

    /*reaction 66: H + CH3O <=> OH + CH3 */
    /*eqcon[65] *= 1;  */

    /*reaction 67: H + CH3O <=> CH2(S) + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: H + CH3OH <=> CH2OH + H2 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: H + CH3OH <=> CH3O + H2 */
    /*eqcon[68] *= 1;  */

    /*reaction 70: H + C2H (+M) <=> C2H2 (+M) */
    eqcon[69] *= 1e+06; 

    /*reaction 71: H + C2H2 (+M) <=> C2H3 (+M) */
    eqcon[70] *= 1e+06; 

    /*reaction 72: H + C2H3 (+M) <=> C2H4 (+M) */
    eqcon[71] *= 1e+06; 

    /*reaction 73: H + C2H3 <=> H2 + C2H2 */
    /*eqcon[72] *= 1;  */

    /*reaction 74: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[73] *= 1e+06; 

    /*reaction 75: H + C2H4 <=> C2H3 + H2 */
    /*eqcon[74] *= 1;  */

    /*reaction 76: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[75] *= 1e+06; 

    /*reaction 77: H + C2H5 <=> H2 + C2H4 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[77] *= 1;  */

    /*reaction 79: H + HCCO <=> CH2(S) + CO */
    /*eqcon[78] *= 1;  */

    /*reaction 80: H + CH2CO <=> HCCO + H2 */
    /*eqcon[79] *= 1;  */

    /*reaction 81: H + CH2CO <=> CH3 + CO */
    /*eqcon[80] *= 1;  */

    /*reaction 82: H + HCCOH <=> H + CH2CO */
    /*eqcon[81] *= 1;  */

    /*reaction 83: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[82] *= 1e+06; 

    /*reaction 84: OH + H2 <=> H + H2O */
    /*eqcon[83] *= 1;  */

    /*reaction 85: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[84] *= 1e+06; 

    /*reaction 86: 2 OH <=> O + H2O */
    /*eqcon[85] *= 1;  */

    /*reaction 87: OH + HO2 <=> O2 + H2O */
    /*eqcon[86] *= 1;  */

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[87] *= 1;  */

    /*reaction 89: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[88] *= 1;  */

    /*reaction 90: OH + C <=> H + CO */
    /*eqcon[89] *= 1;  */

    /*reaction 91: OH + CH <=> H + HCO */
    /*eqcon[90] *= 1;  */

    /*reaction 92: OH + CH2 <=> H + CH2O */
    /*eqcon[91] *= 1;  */

    /*reaction 93: OH + CH2 <=> CH + H2O */
    /*eqcon[92] *= 1;  */

    /*reaction 94: OH + CH2(S) <=> H + CH2O */
    /*eqcon[93] *= 1;  */

    /*reaction 95: OH + CH3 (+M) <=> CH3OH (+M) */
    eqcon[94] *= 1e+06; 

    /*reaction 96: OH + CH3 <=> CH2 + H2O */
    /*eqcon[95] *= 1;  */

    /*reaction 97: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[96] *= 1;  */

    /*reaction 98: OH + CH4 <=> CH3 + H2O */
    /*eqcon[97] *= 1;  */

    /*reaction 99: OH + CO <=> H + CO2 */
    /*eqcon[98] *= 1;  */

    /*reaction 100: OH + HCO <=> H2O + CO */
    /*eqcon[99] *= 1;  */

    /*reaction 101: OH + CH2O <=> HCO + H2O */
    /*eqcon[100] *= 1;  */

    /*reaction 102: OH + CH2OH <=> H2O + CH2O */
    /*eqcon[101] *= 1;  */

    /*reaction 103: OH + CH3O <=> H2O + CH2O */
    /*eqcon[102] *= 1;  */

    /*reaction 104: OH + CH3OH <=> CH2OH + H2O */
    /*eqcon[103] *= 1;  */

    /*reaction 105: OH + CH3OH <=> CH3O + H2O */
    /*eqcon[104] *= 1;  */

    /*reaction 106: OH + C2H <=> H + HCCO */
    /*eqcon[105] *= 1;  */

    /*reaction 107: OH + C2H2 <=> H + CH2CO */
    /*eqcon[106] *= 1;  */

    /*reaction 108: OH + C2H2 <=> H + HCCOH */
    /*eqcon[107] *= 1;  */

    /*reaction 109: OH + C2H2 <=> C2H + H2O */
    /*eqcon[108] *= 1;  */

    /*reaction 110: OH + C2H2 <=> CH3 + CO */
    /*eqcon[109] *= 1;  */

    /*reaction 111: OH + C2H3 <=> H2O + C2H2 */
    /*eqcon[110] *= 1;  */

    /*reaction 112: OH + C2H4 <=> C2H3 + H2O */
    /*eqcon[111] *= 1;  */

    /*reaction 113: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[112] *= 1;  */

    /*reaction 114: OH + CH2CO <=> HCCO + H2O */
    /*eqcon[113] *= 1;  */

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[114] *= 1;  */

    /*reaction 116: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[115] *= 1;  */

    /*reaction 117: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[116] *= 1;  */

    /*reaction 118: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[117] *= 1;  */

    /*reaction 119: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[118] *= 1;  */

    /*reaction 120: HO2 + CO <=> OH + CO2 */
    /*eqcon[119] *= 1;  */

    /*reaction 121: HO2 + CH2O <=> HCO + H2O2 */
    /*eqcon[120] *= 1;  */

    /*reaction 122: C + O2 <=> O + CO */
    /*eqcon[121] *= 1;  */

    /*reaction 123: C + CH2 <=> H + C2H */
    /*eqcon[122] *= 1;  */

    /*reaction 124: C + CH3 <=> H + C2H2 */
    /*eqcon[123] *= 1;  */

    /*reaction 125: CH + O2 <=> O + HCO */
    /*eqcon[124] *= 1;  */

    /*reaction 126: CH + H2 <=> H + CH2 */
    /*eqcon[125] *= 1;  */

    /*reaction 127: CH + H2O <=> H + CH2O */
    /*eqcon[126] *= 1;  */

    /*reaction 128: CH + CH2 <=> H + C2H2 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: CH + CH3 <=> H + C2H3 */
    /*eqcon[128] *= 1;  */

    /*reaction 130: CH + CH4 <=> H + C2H4 */
    /*eqcon[129] *= 1;  */

    /*reaction 131: CH + CO (+M) <=> HCCO (+M) */
    eqcon[130] *= 1e+06; 

    /*reaction 132: CH + CO2 <=> HCO + CO */
    /*eqcon[131] *= 1;  */

    /*reaction 133: CH + CH2O <=> H + CH2CO */
    /*eqcon[132] *= 1;  */

    /*reaction 134: CH + HCCO <=> CO + C2H2 */
    /*eqcon[133] *= 1;  */

    /*reaction 135: CH2 + O2 <=> OH + HCO */
    /*eqcon[134] *= 1;  */

    /*reaction 136: CH2 + H2 <=> H + CH3 */
    /*eqcon[135] *= 1;  */

    /*reaction 137: 2 CH2 <=> H2 + C2H2 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[137] *= 1;  */

    /*reaction 139: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[138] *= 1;  */

    /*reaction 140: CH2 + CO (+M) <=> CH2CO (+M) */
    eqcon[139] *= 1e+06; 

    /*reaction 141: CH2 + HCCO <=> C2H3 + CO */
    /*eqcon[140] *= 1;  */

    /*reaction 142: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[141] *= 1;  */

    /*reaction 143: CH2(S) + AR <=> CH2 + AR */
    /*eqcon[142] *= 1;  */

    /*reaction 144: CH2(S) + O2 <=> H + OH + CO */
    eqcon[143] *= 1e-06; 

    /*reaction 145: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[144] *= 1;  */

    /*reaction 146: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[145] *= 1;  */

    /*reaction 147: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    eqcon[146] *= 1e+06; 

    /*reaction 148: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[147] *= 1;  */

    /*reaction 149: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[148] *= 1;  */

    /*reaction 150: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[149] *= 1;  */

    /*reaction 151: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[150] *= 1;  */

    /*reaction 152: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[152] *= 1;  */

    /*reaction 154: CH2(S) + C2H6 <=> CH3 + C2H5 */
    /*eqcon[153] *= 1;  */

    /*reaction 155: CH3 + O2 <=> O + CH3O */
    /*eqcon[154] *= 1;  */

    /*reaction 156: CH3 + O2 <=> OH + CH2O */
    /*eqcon[155] *= 1;  */

    /*reaction 157: CH3 + H2O2 <=> HO2 + CH4 */
    /*eqcon[156] *= 1;  */

    /*reaction 158: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[157] *= 1e+06; 

    /*reaction 159: 2 CH3 <=> H + C2H5 */
    /*eqcon[158] *= 1;  */

    /*reaction 160: CH3 + HCO <=> CH4 + CO */
    /*eqcon[159] *= 1;  */

    /*reaction 161: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[160] *= 1;  */

    /*reaction 162: CH3 + CH3OH <=> CH2OH + CH4 */
    /*eqcon[161] *= 1;  */

    /*reaction 163: CH3 + CH3OH <=> CH3O + CH4 */
    /*eqcon[162] *= 1;  */

    /*reaction 164: CH3 + C2H4 <=> C2H3 + CH4 */
    /*eqcon[163] *= 1;  */

    /*reaction 165: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[164] *= 1;  */

    /*reaction 166: HCO + H2O <=> H + CO + H2O */
    eqcon[165] *= 1e-06; 

    /*reaction 167: HCO + M <=> H + CO + M */
    eqcon[166] *= 1e-06; 

    /*reaction 168: HCO + O2 <=> HO2 + CO */
    /*eqcon[167] *= 1;  */

    /*reaction 169: CH2OH + O2 <=> HO2 + CH2O */
    /*eqcon[168] *= 1;  */

    /*reaction 170: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[169] *= 1;  */

    /*reaction 171: C2H + O2 <=> HCO + CO */
    /*eqcon[170] *= 1;  */

    /*reaction 172: C2H + H2 <=> H + C2H2 */
    /*eqcon[171] *= 1;  */

    /*reaction 173: C2H3 + O2 <=> HCO + CH2O */
    /*eqcon[172] *= 1;  */

    /*reaction 174: C2H4 (+M) <=> H2 + C2H2 (+M) */
    eqcon[173] *= 1e-06; 

    /*reaction 175: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[174] *= 1;  */

    /*reaction 176: HCCO + O2 <=> OH + 2 CO */
    eqcon[175] *= 1e-06; 

    /*reaction 177: 2 HCCO <=> 2 CO + C2H2 */
    eqcon[176] *= 1e-06; 
}


/*compute the production rate for each species */
void productionRate(double * wdot, double * sc, double T)
{
    double qdot;

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[32];                /*Gibbs free energy */
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
    for (id = 0; id < 32; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*zero out wdot */
    for (id = 0; id < 32; ++id) {
        wdot[id] = 0.0;
    }

    /*reaction 1: 2 O + M <=> O2 + M */
    phi_f = sc[2]*sc[2];
    alpha = mixture + 1.4*sc[0] + 14.4*sc[5] + sc[13] + 0.75*sc[14] + 2.6*sc[15] + 2*sc[26] + -0.17*sc[31];
    k_f = 1e-12 * alpha * 1.2e+17*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[3];
    Kc = 1.0 / (refC) * exp((2 * g_RT[2]) - (g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 2 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 2: O + H + M <=> OH + M */
    phi_f = sc[2]*sc[1];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-12 * alpha * 5e+17*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[4];
    Kc = 1.0 / (refC) * exp((g_RT[2] + g_RT[1]) - (g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 3: O + H2 <=> H + OH */
    phi_f = sc[2]*sc[0];
    k_f = 1e-06 * 50000*exp(2.67*tc[0]-3165.58/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[4];
    Kc = exp((g_RT[2] + g_RT[0]) - (g_RT[1] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 4: O + HO2 <=> OH + O2 */
    phi_f = sc[2]*sc[6];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = exp((g_RT[2] + g_RT[6]) - (g_RT[4] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 5: O + H2O2 <=> OH + HO2 */
    phi_f = sc[2]*sc[7];
    k_f = 1e-06 * 9.63e+06*exp(2*tc[0]-2013.09/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[6];
    Kc = exp((g_RT[2] + g_RT[7]) - (g_RT[4] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 6: O + CH <=> H + CO */
    phi_f = sc[2]*sc[9];
    k_f = 1e-06 * 5.7e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = exp((g_RT[2] + g_RT[9]) - (g_RT[1] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[9] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 7: O + CH2 <=> H + HCO */
    phi_f = sc[2]*sc[10];
    k_f = 1e-06 * 8e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = exp((g_RT[2] + g_RT[10]) - (g_RT[1] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 8: O + CH2(S) <=> H2 + CO */
    phi_f = sc[2]*sc[11];
    k_f = 1e-06 * 1.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[14];
    Kc = exp((g_RT[2] + g_RT[11]) - (g_RT[0] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 9: O + CH2(S) <=> H + HCO */
    phi_f = sc[2]*sc[11];
    k_f = 1e-06 * 1.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = exp((g_RT[2] + g_RT[11]) - (g_RT[1] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 10: O + CH3 <=> H + CH2O */
    phi_f = sc[2]*sc[12];
    k_f = 1e-06 * 8.43e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = exp((g_RT[2] + g_RT[12]) - (g_RT[1] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 11: O + CH4 <=> OH + CH3 */
    phi_f = sc[2]*sc[13];
    k_f = 1e-06 * 1.02e+09*exp(1.5*tc[0]-4328.13/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = exp((g_RT[2] + g_RT[13]) - (g_RT[4] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 12: O + CO + M <=> CO2 + M */
    phi_f = sc[2]*sc[14];
    alpha = mixture + sc[0] + 5*sc[3] + 5*sc[5] + sc[13] + 0.5*sc[14] + 2.5*sc[15] + 2*sc[26] + -0.5*sc[31];
    k_f = 1e-12 * alpha * 6.02e+14*exp(-1509.81/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15];
    Kc = 1.0 / (refC) * exp((g_RT[2] + g_RT[14]) - (g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 13: O + HCO <=> OH + CO */
    phi_f = sc[2]*sc[16];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[14];
    Kc = exp((g_RT[2] + g_RT[16]) - (g_RT[4] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[16] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 14: O + HCO <=> H + CO2 */
    phi_f = sc[2]*sc[16];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[15];
    Kc = exp((g_RT[2] + g_RT[16]) - (g_RT[1] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[16] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 15: O + CH2O <=> OH + HCO */
    phi_f = sc[2]*sc[17];
    k_f = 1e-06 * 3.9e+13*exp(-1781.58/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[16];
    Kc = exp((g_RT[2] + g_RT[17]) - (g_RT[4] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 16: O + CH2OH <=> OH + CH2O */
    phi_f = sc[2]*sc[18];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = exp((g_RT[2] + g_RT[18]) - (g_RT[4] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 17: O + CH3O <=> OH + CH2O */
    phi_f = sc[2]*sc[19];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = exp((g_RT[2] + g_RT[19]) - (g_RT[4] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 18: O + CH3OH <=> OH + CH2OH */
    phi_f = sc[2]*sc[20];
    k_f = 1e-06 * 388000*exp(2.5*tc[0]-1560.14/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[18];
    Kc = exp((g_RT[2] + g_RT[20]) - (g_RT[4] + g_RT[18]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[18] += 1 * qdot;

    /*reaction 19: O + CH3OH <=> OH + CH3O */
    phi_f = sc[2]*sc[20];
    k_f = 1e-06 * 130000*exp(2.5*tc[0]-2516.36/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[19];
    Kc = exp((g_RT[2] + g_RT[20]) - (g_RT[4] + g_RT[19]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[19] += 1 * qdot;

    /*reaction 20: O + C2H <=> CH + CO */
    phi_f = sc[2]*sc[21];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[14];
    Kc = exp((g_RT[2] + g_RT[21]) - (g_RT[9] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[21] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 21: O + C2H2 <=> H + HCCO */
    phi_f = sc[2]*sc[22];
    k_f = 1e-06 * 1.02e+07*exp(2*tc[0]-956.215/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[27];
    Kc = exp((g_RT[2] + g_RT[22]) - (g_RT[1] + g_RT[27]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[27] += 1 * qdot;

    /*reaction 22: O + C2H2 <=> OH + C2H */
    phi_f = sc[2]*sc[22];
    k_f = 1e-06 * 4.6e+19*exp(-1.41*tc[0]-14569.7/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[21];
    Kc = exp((g_RT[2] + g_RT[22]) - (g_RT[4] + g_RT[21]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[21] += 1 * qdot;

    /*reaction 23: O + C2H2 <=> CO + CH2 */
    phi_f = sc[2]*sc[22];
    k_f = 1e-06 * 1.02e+07*exp(2*tc[0]-956.215/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[10];
    Kc = exp((g_RT[2] + g_RT[22]) - (g_RT[14] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 24: O + C2H3 <=> H + CH2CO */
    phi_f = sc[2]*sc[23];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = exp((g_RT[2] + g_RT[23]) - (g_RT[1] + g_RT[28]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[23] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[28] += 1 * qdot;

    /*reaction 25: O + C2H4 <=> CH3 + HCO */
    phi_f = sc[2]*sc[24];
    k_f = 1e-06 * 1.92e+07*exp(1.83*tc[0]-110.72/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[16];
    Kc = exp((g_RT[2] + g_RT[24]) - (g_RT[12] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[24] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 26: O + C2H5 <=> CH3 + CH2O */
    phi_f = sc[2]*sc[25];
    k_f = 1e-06 * 1.32e+14;
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[17];
    Kc = exp((g_RT[2] + g_RT[25]) - (g_RT[12] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[25] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 27: O + C2H6 <=> OH + C2H5 */
    phi_f = sc[2]*sc[26];
    k_f = 1e-06 * 8.98e+07*exp(1.92*tc[0]-2863.61/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[25];
    Kc = exp((g_RT[2] + g_RT[26]) - (g_RT[4] + g_RT[25]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[26] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[25] += 1 * qdot;

    /*reaction 28: O + HCCO <=> H + 2 CO */
    phi_f = sc[2]*sc[27];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14]*sc[14];
    Kc = refC * exp((g_RT[2] + g_RT[27]) - (g_RT[1] + 2 * g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[27] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[14] += 2 * qdot;

    /*reaction 29: O + CH2CO <=> OH + HCCO */
    phi_f = sc[2]*sc[28];
    k_f = 1e-06 * 1e+13*exp(-4026.17/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[27];
    Kc = exp((g_RT[2] + g_RT[28]) - (g_RT[4] + g_RT[27]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[28] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[27] += 1 * qdot;

    /*reaction 30: O + CH2CO <=> CH2 + CO2 */
    phi_f = sc[2]*sc[28];
    k_f = 1e-06 * 1.75e+12*exp(-679.416/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[15];
    Kc = exp((g_RT[2] + g_RT[28]) - (g_RT[10] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[28] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 31: O2 + CO <=> O + CO2 */
    phi_f = sc[3]*sc[14];
    k_f = 1e-06 * 2.5e+12*exp(-24056.4/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[15];
    Kc = exp((g_RT[3] + g_RT[14]) - (g_RT[2] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 32: O2 + CH2O <=> HO2 + HCO */
    phi_f = sc[3]*sc[17];
    k_f = 1e-06 * 1e+14*exp(-20130.9/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[16];
    Kc = exp((g_RT[3] + g_RT[17]) - (g_RT[6] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 33: H + O2 + M <=> HO2 + M */
    phi_f = sc[1]*sc[3];
    alpha = mixture + -1*sc[3] + -1*sc[5] + -0.25*sc[14] + 0.5*sc[15] + 0.5*sc[26] + -1*sc[30] + -1*sc[31];
    k_f = 1e-12 * alpha * 2.8e+18*exp(-0.86*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 34: H + 2 O2 <=> HO2 + O2 */
    phi_f = sc[1]*sc[3]*sc[3];
    k_f = 1e-12 * 3e+20*exp(-1.72*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = 1.0 / (refC) * exp((g_RT[1] + 2 * g_RT[3]) - (g_RT[6] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 2 * qdot;
    wdot[6] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 35: H + O2 + H2O <=> HO2 + H2O */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = 1e-12 * 9.38e+18*exp(-0.76*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 36: H + O2 + N2 <=> HO2 + N2 */
    phi_f = sc[1]*sc[3]*sc[30];
    k_f = 1e-12 * 3.75e+20*exp(-1.72*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[30];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[30]) - (g_RT[6] + g_RT[30]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[30] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[30] += 1 * qdot;

    /*reaction 37: H + O2 + AR <=> HO2 + AR */
    phi_f = sc[1]*sc[3]*sc[31];
    k_f = 1e-12 * 7e+17*exp(-0.8*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[31];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[31]) - (g_RT[6] + g_RT[31]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[31] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[31] += 1 * qdot;

    /*reaction 38: H + O2 <=> O + OH */
    phi_f = sc[1]*sc[3];
    k_f = 1e-06 * 8.3e+13*exp(-7253.65/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 39: 2 H + M <=> H2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[0] + -1*sc[5] + sc[13] + -1*sc[15] + 2*sc[26] + -0.37*sc[31];
    k_f = 1e-12 * alpha * 1e+18*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    phi_f = sc[1]*sc[1]*sc[0];
    k_f = 1e-12 * 9e+16*exp(-0.6*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[0];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[0] += 2 * qdot;

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    phi_f = sc[1]*sc[1]*sc[5];
    k_f = 1e-12 * 6e+19*exp(-1.25*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[5];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    phi_f = sc[1]*sc[1]*sc[15];
    k_f = 1e-12 * 5.5e+20*exp(-2*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[15];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[15]) - (g_RT[0] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 43: H + OH + M <=> H2O + M */
    phi_f = sc[1]*sc[4];
    alpha = mixture + -0.27*sc[0] + 2.65*sc[5] + sc[13] + 2*sc[26] + -0.62*sc[31];
    k_f = 1e-12 * alpha * 2.2e+22*exp(-2*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 44: H + HO2 <=> O + H2O */
    phi_f = sc[1]*sc[6];
    k_f = 1e-06 * 3.97e+12*exp(-337.695/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = exp((g_RT[1] + g_RT[6]) - (g_RT[2] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 45: H + HO2 <=> O2 + H2 */
    phi_f = sc[1]*sc[6];
    k_f = 1e-06 * 2.8e+13*exp(-537.494/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 46: H + HO2 <=> 2 OH */
    phi_f = sc[1]*sc[6];
    k_f = 1e-06 * 1.34e+14*exp(-319.577/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[4] += 2 * qdot;

    /*reaction 47: H + H2O2 <=> HO2 + H2 */
    phi_f = sc[1]*sc[7];
    k_f = 1e-06 * 1.21e+07*exp(2*tc[0]-2617.01/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = exp((g_RT[1] + g_RT[7]) - (g_RT[6] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 48: H + H2O2 <=> OH + H2O */
    phi_f = sc[1]*sc[7];
    k_f = 1e-06 * 1e+13*exp(-1811.78/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[5];
    Kc = exp((g_RT[1] + g_RT[7]) - (g_RT[4] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 49: H + CH <=> C + H2 */
    phi_f = sc[1]*sc[9];
    k_f = 1e-06 * 1.1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = exp((g_RT[1] + g_RT[9]) - (g_RT[8] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[9] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 50: H + CH2 (+M) <=> CH3 (+M) */
    phi_f = sc[1]*sc[10];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 2.5e+16*exp(-0.8*tc[0]);
    redP = 1e-12 * alpha / k_f * 3.2e+27*exp(-3.14*tc[0]-619.024/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.32*exp(T/-78))+ (0.68*exp(T/-1995))+ (exp(-5590/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[12];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 51: H + CH2(S) <=> CH + H2 */
    phi_f = sc[1]*sc[11];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = exp((g_RT[1] + g_RT[11]) - (g_RT[9] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 52: H + CH3 (+M) <=> CH4 (+M) */
    phi_f = sc[1]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 1.27e+16*exp(-0.63*tc[0]-192.753/tc[1]);
    redP = 1e-12 * alpha / k_f * 2.477e+33*exp(-4.76*tc[0]-1227.98/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.217*exp(T/-74))+ (0.783*exp(T/-2941))+ (exp(-6964/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[13];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[12]) - (g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 53: H + CH4 <=> CH3 + H2 */
    phi_f = sc[1]*sc[13];
    k_f = 1e-06 * 6.6e+08*exp(1.62*tc[0]-5455.46/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[0];
    Kc = exp((g_RT[1] + g_RT[13]) - (g_RT[12] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 54: H + HCO (+M) <=> CH2O (+M) */
    phi_f = sc[1]*sc[16];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 1.09e+12*exp(0.48*tc[0]+130.851/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.35e+24*exp(-2.57*tc[0]-717.162/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2176*exp(T/-271))+ (0.7824*exp(T/-2755))+ (exp(-6570/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[17];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[16]) - (g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[16] -= 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 55: H + HCO <=> H2 + CO */
    phi_f = sc[1]*sc[16];
    k_f = 1e-06 * 7.34e+13;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[14];
    Kc = exp((g_RT[1] + g_RT[16]) - (g_RT[0] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[16] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 56: H + CH2O (+M) <=> CH2OH (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1811.78/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.27e+32*exp(-4.82*tc[0]-3286.36/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2813*exp(T/-103))+ (0.7187*exp(T/-1291))+ (exp(-4160/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[18];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[18]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[18] += 1 * qdot;

    /*reaction 57: H + CH2O (+M) <=> CH3O (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1308.51/tc[1]);
    redP = 1e-12 * alpha / k_f * 2.2e+30*exp(-4.8*tc[0]-2798.19/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.242*exp(T/-94))+ (0.758*exp(T/-1555))+ (exp(-4200/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[19];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[19]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[19] += 1 * qdot;

    /*reaction 58: H + CH2O <=> HCO + H2 */
    phi_f = sc[1]*sc[17];
    k_f = 1e-06 * 2.3e+10*exp(1.05*tc[0]-1648.21/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[0];
    Kc = exp((g_RT[1] + g_RT[17]) - (g_RT[16] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 59: H + CH2OH (+M) <=> CH3OH (+M) */
    phi_f = sc[1]*sc[18];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = 1e-06 * 1.8e+13;
    redP = 1e-12 * alpha / k_f * 3e+31*exp(-4.8*tc[0]-1660.8/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2321*exp(T/-338))+ (0.7679*exp(T/-1812))+ (exp(-5081/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[20];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[18]) - (g_RT[20]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[20] += 1 * qdot;

    /*reaction 60: H + CH2OH <=> H2 + CH2O */
    phi_f = sc[1]*sc[18];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[17];
    Kc = exp((g_RT[1] + g_RT[18]) - (g_RT[0] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 61: H + CH2OH <=> OH + CH3 */
    phi_f = sc[1]*sc[18];
    k_f = 1e-06 * 1.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = exp((g_RT[1] + g_RT[18]) - (g_RT[4] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 62: H + CH2OH <=> CH2(S) + H2O */
    phi_f = sc[1]*sc[18];
    k_f = 1e-06 * 6e+12;
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = exp((g_RT[1] + g_RT[18]) - (g_RT[11] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 63: H + CH3O (+M) <=> CH3OH (+M) */
    phi_f = sc[1]*sc[19];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = 1e-06 * 5e+13;
    redP = 1e-12 * alpha / k_f * 8.6e+28*exp(-4*tc[0]-1522.4/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.1098*exp(T/-144))+ (0.8902*exp(T/-2838))+ (exp(-45569/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[20];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[19]) - (g_RT[20]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[20] += 1 * qdot;

    /*reaction 64: H + CH3O <=> H + CH2OH */
    phi_f = sc[1]*sc[19];
    k_f = 1e-06 * 3.4e+06*exp(1.6*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[18];
    Kc = exp((g_RT[1] + g_RT[19]) - (g_RT[1] + g_RT[18]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[18] += 1 * qdot;

    /*reaction 65: H + CH3O <=> H2 + CH2O */
    phi_f = sc[1]*sc[19];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[17];
    Kc = exp((g_RT[1] + g_RT[19]) - (g_RT[0] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 66: H + CH3O <=> OH + CH3 */
    phi_f = sc[1]*sc[19];
    k_f = 1e-06 * 3.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = exp((g_RT[1] + g_RT[19]) - (g_RT[4] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 67: H + CH3O <=> CH2(S) + H2O */
    phi_f = sc[1]*sc[19];
    k_f = 1e-06 * 1.6e+13;
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = exp((g_RT[1] + g_RT[19]) - (g_RT[11] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 68: H + CH3OH <=> CH2OH + H2 */
    phi_f = sc[1]*sc[20];
    k_f = 1e-06 * 1.7e+07*exp(2.1*tc[0]-2450.93/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[0];
    Kc = exp((g_RT[1] + g_RT[20]) - (g_RT[18] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 69: H + CH3OH <=> CH3O + H2 */
    phi_f = sc[1]*sc[20];
    k_f = 1e-06 * 4.2e+06*exp(2.1*tc[0]-2450.93/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[0];
    Kc = exp((g_RT[1] + g_RT[20]) - (g_RT[19] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[19] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 70: H + C2H (+M) <=> C2H2 (+M) */
    phi_f = sc[1]*sc[21];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 1e+17*exp(-1*tc[0]);
    redP = 1e-12 * alpha / k_f * 3.75e+33*exp(-4.8*tc[0]-956.215/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.3536*exp(T/-132))+ (0.6464*exp(T/-1315))+ (exp(-5566/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[22];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[21]) - (g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[21] -= 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 71: H + C2H2 (+M) <=> C2H3 (+M) */
    phi_f = sc[1]*sc[22];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 5.6e+12*exp(-1207.85/tc[1]);
    redP = 1e-12 * alpha / k_f * 3.8e+40*exp(-7.27*tc[0]-3633.62/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2493*exp(T/-98.5))+ (0.7507*exp(T/-1302))+ (exp(-4167/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[23];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[22]) - (g_RT[23]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[23] += 1 * qdot;

    /*reaction 72: H + C2H3 (+M) <=> C2H4 (+M) */
    phi_f = sc[1]*sc[23];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 6.08e+12*exp(0.27*tc[0]-140.916/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.4e+30*exp(-3.86*tc[0]-1670.86/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.218*exp(T/-207.5))+ (0.782*exp(T/-2663))+ (exp(-6095/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[24];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[23]) - (g_RT[24]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[23] -= 1 * qdot;
    wdot[24] += 1 * qdot;

    /*reaction 73: H + C2H3 <=> H2 + C2H2 */
    phi_f = sc[1]*sc[23];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = exp((g_RT[1] + g_RT[23]) - (g_RT[0] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[23] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 74: H + C2H4 (+M) <=> C2H5 (+M) */
    phi_f = sc[1]*sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 1.08e+12*exp(0.454*tc[0]-915.954/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.2e+42*exp(-7.62*tc[0]-3507.8/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.0247*exp(T/-210))+ (0.9753*exp(T/-984))+ (exp(-4374/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[25];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[24]) - (g_RT[25]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[24] -= 1 * qdot;
    wdot[25] += 1 * qdot;

    /*reaction 75: H + C2H4 <=> C2H3 + H2 */
    phi_f = sc[1]*sc[24];
    k_f = 1e-06 * 1.325e+06*exp(2.53*tc[0]-6160.04/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[0];
    Kc = exp((g_RT[1] + g_RT[24]) - (g_RT[23] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[24] -= 1 * qdot;
    wdot[23] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 76: H + C2H5 (+M) <=> C2H6 (+M) */
    phi_f = sc[1]*sc[25];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 5.21e+17*exp(-0.99*tc[0]-795.169/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.99e+41*exp(-7.08*tc[0]-3364.37/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.1578*exp(T/-125))+ (0.8422*exp(T/-2219))+ (exp(-6882/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[26];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[25]) - (g_RT[26]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[25] -= 1 * qdot;
    wdot[26] += 1 * qdot;

    /*reaction 77: H + C2H5 <=> H2 + C2H4 */
    phi_f = sc[1]*sc[25];
    k_f = 1e-06 * 2e+12;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[24];
    Kc = exp((g_RT[1] + g_RT[25]) - (g_RT[0] + g_RT[24]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[25] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[24] += 1 * qdot;

    /*reaction 78: H + C2H6 <=> C2H5 + H2 */
    phi_f = sc[1]*sc[26];
    k_f = 1e-06 * 1.15e+08*exp(1.9*tc[0]-3789.63/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[0];
    Kc = exp((g_RT[1] + g_RT[26]) - (g_RT[25] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[26] -= 1 * qdot;
    wdot[25] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 79: H + HCCO <=> CH2(S) + CO */
    phi_f = sc[1]*sc[27];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[14];
    Kc = exp((g_RT[1] + g_RT[27]) - (g_RT[11] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[27] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 80: H + CH2CO <=> HCCO + H2 */
    phi_f = sc[1]*sc[28];
    k_f = 1e-06 * 5e+13*exp(-4026.17/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[0];
    Kc = exp((g_RT[1] + g_RT[28]) - (g_RT[27] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[28] -= 1 * qdot;
    wdot[27] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 81: H + CH2CO <=> CH3 + CO */
    phi_f = sc[1]*sc[28];
    k_f = 1e-06 * 1.13e+13*exp(-1725.21/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[14];
    Kc = exp((g_RT[1] + g_RT[28]) - (g_RT[12] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[28] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 82: H + HCCOH <=> H + CH2CO */
    phi_f = sc[1]*sc[29];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = exp((g_RT[1] + g_RT[29]) - (g_RT[1] + g_RT[28]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[29] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[28] += 1 * qdot;

    /*reaction 83: H2 + CO (+M) <=> CH2O (+M) */
    phi_f = sc[0]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 4.3e+07*exp(1.5*tc[0]-40060.4/tc[1]);
    redP = 1e-12 * alpha / k_f * 5.07e+27*exp(-3.42*tc[0]-42450.9/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.068*exp(T/-197))+ (0.932*exp(T/-1540))+ (exp(-10300/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[17];
    Kc = 1.0 / (refC) * exp((g_RT[0] + g_RT[14]) - (g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 84: OH + H2 <=> H + H2O */
    phi_f = sc[4]*sc[0];
    k_f = 1e-06 * 2.16e+08*exp(1.51*tc[0]-1726.22/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    Kc = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 85: 2 OH (+M) <=> H2O2 (+M) */
    phi_f = sc[4]*sc[4];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 7.4e+13*exp(-0.37*tc[0]);
    redP = 1e-12 * alpha / k_f * 2.3e+18*exp(-0.9*tc[0]+855.561/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2654*exp(T/-94))+ (0.7346*exp(T/-1756))+ (exp(-5182/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[7];
    Kc = 1.0 / (refC) * exp((2 * g_RT[4]) - (g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 2 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 86: 2 OH <=> O + H2O */
    phi_f = sc[4]*sc[4];
    k_f = 1e-06 * 35700*exp(2.4*tc[0]+1061.9/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 2 * qdot;
    wdot[2] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 87: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = 1e-06 * 2.9e+13*exp(+251.636/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    phi_f = sc[4]*sc[7];
    k_f = 1e-06 * 1.75e+12*exp(-161.047/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 89: OH + H2O2 <=> HO2 + H2O */
    phi_f = sc[4]*sc[7];
    k_f = 1e-06 * 5.8e+14*exp(-4811.27/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 90: OH + C <=> H + CO */
    phi_f = sc[4]*sc[8];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = exp((g_RT[4] + g_RT[8]) - (g_RT[1] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 91: OH + CH <=> H + HCO */
    phi_f = sc[4]*sc[9];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = exp((g_RT[4] + g_RT[9]) - (g_RT[1] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[9] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 92: OH + CH2 <=> H + CH2O */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = exp((g_RT[4] + g_RT[10]) - (g_RT[1] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 93: OH + CH2 <=> CH + H2O */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 1.13e+07*exp(2*tc[0]-1509.81/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[5];
    Kc = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 94: OH + CH2(S) <=> H + CH2O */
    phi_f = sc[4]*sc[11];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = exp((g_RT[4] + g_RT[11]) - (g_RT[1] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 95: OH + CH3 (+M) <=> CH3OH (+M) */
    phi_f = sc[4]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = 1e-06 * 6.3e+13;
    redP = 1e-12 * alpha / k_f * 2.7e+38*exp(-6.3*tc[0]-1560.14/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.7895*exp(T/-83.5))+ (0.2105*exp(T/-5398))+ (exp(-8370/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[20];
    Kc = 1.0 / (refC) * exp((g_RT[4] + g_RT[12]) - (g_RT[20]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[20] += 1 * qdot;

    /*reaction 96: OH + CH3 <=> CH2 + H2O */
    phi_f = sc[4]*sc[12];
    k_f = 1e-06 * 5.6e+07*exp(1.6*tc[0]-2727.73/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[5];
    Kc = exp((g_RT[4] + g_RT[12]) - (g_RT[10] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 97: OH + CH3 <=> CH2(S) + H2O */
    phi_f = sc[4]*sc[12];
    k_f = 1e-06 * 2.501e+13;
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = exp((g_RT[4] + g_RT[12]) - (g_RT[11] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 98: OH + CH4 <=> CH3 + H2O */
    phi_f = sc[4]*sc[13];
    k_f = 1e-06 * 1e+08*exp(1.6*tc[0]-1570.21/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[5];
    Kc = exp((g_RT[4] + g_RT[13]) - (g_RT[12] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 99: OH + CO <=> H + CO2 */
    phi_f = sc[4]*sc[14];
    k_f = 1e-06 * 4.76e+07*exp(1.228*tc[0]-35.229/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[15];
    Kc = exp((g_RT[4] + g_RT[14]) - (g_RT[1] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 100: OH + HCO <=> H2O + CO */
    phi_f = sc[4]*sc[16];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[14];
    Kc = exp((g_RT[4] + g_RT[16]) - (g_RT[5] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[16] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 101: OH + CH2O <=> HCO + H2O */
    phi_f = sc[4]*sc[17];
    k_f = 1e-06 * 3.43e+09*exp(1.18*tc[0]+224.962/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[5];
    Kc = exp((g_RT[4] + g_RT[17]) - (g_RT[16] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 102: OH + CH2OH <=> H2O + CH2O */
    phi_f = sc[4]*sc[18];
    k_f = 1e-06 * 5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[17];
    Kc = exp((g_RT[4] + g_RT[18]) - (g_RT[5] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 103: OH + CH3O <=> H2O + CH2O */
    phi_f = sc[4]*sc[19];
    k_f = 1e-06 * 5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[17];
    Kc = exp((g_RT[4] + g_RT[19]) - (g_RT[5] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 104: OH + CH3OH <=> CH2OH + H2O */
    phi_f = sc[4]*sc[20];
    k_f = 1e-06 * 1.44e+06*exp(2*tc[0]+422.748/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[5];
    Kc = exp((g_RT[4] + g_RT[20]) - (g_RT[18] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 105: OH + CH3OH <=> CH3O + H2O */
    phi_f = sc[4]*sc[20];
    k_f = 1e-06 * 6.3e+06*exp(2*tc[0]-754.907/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[5];
    Kc = exp((g_RT[4] + g_RT[20]) - (g_RT[19] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[19] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 106: OH + C2H <=> H + HCCO */
    phi_f = sc[4]*sc[21];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[27];
    Kc = exp((g_RT[4] + g_RT[21]) - (g_RT[1] + g_RT[27]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[21] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[27] += 1 * qdot;

    /*reaction 107: OH + C2H2 <=> H + CH2CO */
    phi_f = sc[4]*sc[22];
    k_f = 1e-06 * 0.000218*exp(4.5*tc[0]+503.271/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[28]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[28] += 1 * qdot;

    /*reaction 108: OH + C2H2 <=> H + HCCOH */
    phi_f = sc[4]*sc[22];
    k_f = 1e-06 * 504000*exp(2.3*tc[0]-6794.16/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[29];
    Kc = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[29]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[29] += 1 * qdot;

    /*reaction 109: OH + C2H2 <=> C2H + H2O */
    phi_f = sc[4]*sc[22];
    k_f = 1e-06 * 3.37e+07*exp(2*tc[0]-7045.8/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[5];
    Kc = exp((g_RT[4] + g_RT[22]) - (g_RT[21] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[21] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 110: OH + C2H2 <=> CH3 + CO */
    phi_f = sc[4]*sc[22];
    k_f = 1e-06 * 0.000483*exp(4*tc[0]+1006.54/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[14];
    Kc = exp((g_RT[4] + g_RT[22]) - (g_RT[12] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 111: OH + C2H3 <=> H2O + C2H2 */
    phi_f = sc[4]*sc[23];
    k_f = 1e-06 * 5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[22];
    Kc = exp((g_RT[4] + g_RT[23]) - (g_RT[5] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[23] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 112: OH + C2H4 <=> C2H3 + H2O */
    phi_f = sc[4]*sc[24];
    k_f = 1e-06 * 3.6e+06*exp(2*tc[0]-1258.18/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[5];
    Kc = exp((g_RT[4] + g_RT[24]) - (g_RT[23] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[24] -= 1 * qdot;
    wdot[23] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 113: OH + C2H6 <=> C2H5 + H2O */
    phi_f = sc[4]*sc[26];
    k_f = 1e-06 * 3.54e+06*exp(2.12*tc[0]-437.846/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[5];
    Kc = exp((g_RT[4] + g_RT[26]) - (g_RT[25] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[26] -= 1 * qdot;
    wdot[25] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 114: OH + CH2CO <=> HCCO + H2O */
    phi_f = sc[4]*sc[28];
    k_f = 1e-06 * 7.5e+12*exp(-1006.54/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[5];
    Kc = exp((g_RT[4] + g_RT[28]) - (g_RT[27] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[28] -= 1 * qdot;
    wdot[27] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    phi_f = sc[6]*sc[6];
    k_f = 1e-06 * 1.3e+11*exp(+820.332/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 2 * qdot;
    wdot[3] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 116: 2 HO2 <=> O2 + H2O2 */
    phi_f = sc[6]*sc[6];
    k_f = 1e-06 * 4.2e+14*exp(-6039.26/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 2 * qdot;
    wdot[3] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 117: HO2 + CH2 <=> OH + CH2O */
    phi_f = sc[6]*sc[10];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = exp((g_RT[6] + g_RT[10]) - (g_RT[4] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 118: HO2 + CH3 <=> O2 + CH4 */
    phi_f = sc[6]*sc[12];
    k_f = 1e-06 * 1e+12;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[13];
    Kc = exp((g_RT[6] + g_RT[12]) - (g_RT[3] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 119: HO2 + CH3 <=> OH + CH3O */
    phi_f = sc[6]*sc[12];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[19];
    Kc = exp((g_RT[6] + g_RT[12]) - (g_RT[4] + g_RT[19]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[19] += 1 * qdot;

    /*reaction 120: HO2 + CO <=> OH + CO2 */
    phi_f = sc[6]*sc[14];
    k_f = 1e-06 * 1.5e+14*exp(-11877.2/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[15];
    Kc = exp((g_RT[6] + g_RT[14]) - (g_RT[4] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 121: HO2 + CH2O <=> HCO + H2O2 */
    phi_f = sc[6]*sc[17];
    k_f = 1e-06 * 1e+12*exp(-4026.17/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[7];
    Kc = exp((g_RT[6] + g_RT[17]) - (g_RT[16] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 122: C + O2 <=> O + CO */
    phi_f = sc[8]*sc[3];
    k_f = 1e-06 * 5.8e+13*exp(-289.884/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[14];
    Kc = exp((g_RT[8] + g_RT[3]) - (g_RT[2] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 123: C + CH2 <=> H + C2H */
    phi_f = sc[8]*sc[10];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[21];
    Kc = exp((g_RT[8] + g_RT[10]) - (g_RT[1] + g_RT[21]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[21] += 1 * qdot;

    /*reaction 124: C + CH3 <=> H + C2H2 */
    phi_f = sc[8]*sc[12];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = exp((g_RT[8] + g_RT[12]) - (g_RT[1] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 125: CH + O2 <=> O + HCO */
    phi_f = sc[9]*sc[3];
    k_f = 1e-06 * 3.3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[16];
    Kc = exp((g_RT[9] + g_RT[3]) - (g_RT[2] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 126: CH + H2 <=> H + CH2 */
    phi_f = sc[9]*sc[0];
    k_f = 1e-06 * 1.107e+08*exp(1.79*tc[0]-840.463/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[10];
    Kc = exp((g_RT[9] + g_RT[0]) - (g_RT[1] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 127: CH + H2O <=> H + CH2O */
    phi_f = sc[9]*sc[5];
    k_f = 1e-06 * 5.71e+12*exp(+379.97/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = exp((g_RT[9] + g_RT[5]) - (g_RT[1] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 128: CH + CH2 <=> H + C2H2 */
    phi_f = sc[9]*sc[10];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = exp((g_RT[9] + g_RT[10]) - (g_RT[1] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 129: CH + CH3 <=> H + C2H3 */
    phi_f = sc[9]*sc[12];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[23];
    Kc = exp((g_RT[9] + g_RT[12]) - (g_RT[1] + g_RT[23]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[23] += 1 * qdot;

    /*reaction 130: CH + CH4 <=> H + C2H4 */
    phi_f = sc[9]*sc[13];
    k_f = 1e-06 * 6e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = exp((g_RT[9] + g_RT[13]) - (g_RT[1] + g_RT[24]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[24] += 1 * qdot;

    /*reaction 131: CH + CO (+M) <=> HCCO (+M) */
    phi_f = sc[9]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 5e+13;
    redP = 1e-12 * alpha / k_f * 2.69e+28*exp(-3.74*tc[0]-974.333/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.4243*exp(T/-237))+ (0.5757*exp(T/-1652))+ (exp(-5069/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[27];
    Kc = 1.0 / (refC) * exp((g_RT[9] + g_RT[14]) - (g_RT[27]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[27] += 1 * qdot;

    /*reaction 132: CH + CO2 <=> HCO + CO */
    phi_f = sc[9]*sc[15];
    k_f = 1e-06 * 3.4e+12*exp(-347.257/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[14];
    Kc = exp((g_RT[9] + g_RT[15]) - (g_RT[16] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 133: CH + CH2O <=> H + CH2CO */
    phi_f = sc[9]*sc[17];
    k_f = 1e-06 * 9.46e+13*exp(+259.185/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = exp((g_RT[9] + g_RT[17]) - (g_RT[1] + g_RT[28]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[28] += 1 * qdot;

    /*reaction 134: CH + HCCO <=> CO + C2H2 */
    phi_f = sc[9]*sc[27];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[22];
    Kc = exp((g_RT[9] + g_RT[27]) - (g_RT[14] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[27] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 135: CH2 + O2 <=> OH + HCO */
    phi_f = sc[10]*sc[3];
    k_f = 1e-06 * 1.32e+13*exp(-754.907/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[16];
    Kc = exp((g_RT[10] + g_RT[3]) - (g_RT[4] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 136: CH2 + H2 <=> H + CH3 */
    phi_f = sc[10]*sc[0];
    k_f = 1e-06 * 500000*exp(2*tc[0]-3638.65/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = exp((g_RT[10] + g_RT[0]) - (g_RT[1] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 137: 2 CH2 <=> H2 + C2H2 */
    phi_f = sc[10]*sc[10];
    k_f = 1e-06 * 3.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = exp((2 * g_RT[10]) - (g_RT[0] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 2 * qdot;
    wdot[0] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 138: CH2 + CH3 <=> H + C2H4 */
    phi_f = sc[10]*sc[12];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = exp((g_RT[10] + g_RT[12]) - (g_RT[1] + g_RT[24]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[24] += 1 * qdot;

    /*reaction 139: CH2 + CH4 <=> 2 CH3 */
    phi_f = sc[10]*sc[13];
    k_f = 1e-06 * 2.46e+06*exp(2*tc[0]-4162.05/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[12];
    Kc = exp((g_RT[10] + g_RT[13]) - (2 * g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[12] += 2 * qdot;

    /*reaction 140: CH2 + CO (+M) <=> CH2CO (+M) */
    phi_f = sc[10]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 8.1e+11*exp(0.5*tc[0]-2269.75/tc[1]);
    redP = 1e-12 * alpha / k_f * 2.69e+33*exp(-5.11*tc[0]-3570.71/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.4093*exp(T/-275))+ (0.5907*exp(T/-1226))+ (exp(-5185/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[28];
    Kc = 1.0 / (refC) * exp((g_RT[10] + g_RT[14]) - (g_RT[28]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[28] += 1 * qdot;

    /*reaction 141: CH2 + HCCO <=> C2H3 + CO */
    phi_f = sc[10]*sc[27];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[14];
    Kc = exp((g_RT[10] + g_RT[27]) - (g_RT[23] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[27] -= 1 * qdot;
    wdot[23] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 142: CH2(S) + N2 <=> CH2 + N2 */
    phi_f = sc[11]*sc[30];
    k_f = 1e-06 * 1.5e+13*exp(-301.963/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[30];
    Kc = exp((g_RT[11] + g_RT[30]) - (g_RT[10] + g_RT[30]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[30] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[30] += 1 * qdot;

    /*reaction 143: CH2(S) + AR <=> CH2 + AR */
    phi_f = sc[11]*sc[31];
    k_f = 1e-06 * 9e+12*exp(-301.963/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[31];
    Kc = exp((g_RT[11] + g_RT[31]) - (g_RT[10] + g_RT[31]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[31] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[31] += 1 * qdot;

    /*reaction 144: CH2(S) + O2 <=> H + OH + CO */
    phi_f = sc[11]*sc[3];
    k_f = 1e-06 * 2.8e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[4]*sc[14];
    Kc = refC * exp((g_RT[11] + g_RT[3]) - (g_RT[1] + g_RT[4] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 145: CH2(S) + O2 <=> CO + H2O */
    phi_f = sc[11]*sc[3];
    k_f = 1e-06 * 1.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[5];
    Kc = exp((g_RT[11] + g_RT[3]) - (g_RT[14] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 146: CH2(S) + H2 <=> CH3 + H */
    phi_f = sc[11]*sc[0];
    k_f = 1e-06 * 7e+13;
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[1];
    Kc = exp((g_RT[11] + g_RT[0]) - (g_RT[12] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 147: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    phi_f = sc[11]*sc[5];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = 1e-06 * 2e+13;
    redP = 1e-12 * alpha / k_f * 2.7e+38*exp(-6.3*tc[0]-1560.14/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.8493*exp(T/-134))+ (0.1507*exp(T/-2383))+ (exp(-7265/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[20];
    Kc = 1.0 / (refC) * exp((g_RT[11] + g_RT[5]) - (g_RT[20]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[20] += 1 * qdot;

    /*reaction 148: CH2(S) + H2O <=> CH2 + H2O */
    phi_f = sc[11]*sc[5];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[5];
    Kc = exp((g_RT[11] + g_RT[5]) - (g_RT[10] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 149: CH2(S) + CH3 <=> H + C2H4 */
    phi_f = sc[11]*sc[12];
    k_f = 1e-06 * 1.2e+13*exp(+286.865/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = exp((g_RT[11] + g_RT[12]) - (g_RT[1] + g_RT[24]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[24] += 1 * qdot;

    /*reaction 150: CH2(S) + CH4 <=> 2 CH3 */
    phi_f = sc[11]*sc[13];
    k_f = 1e-06 * 1.6e+13*exp(+286.865/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[12];
    Kc = exp((g_RT[11] + g_RT[13]) - (2 * g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[12] += 2 * qdot;

    /*reaction 151: CH2(S) + CO <=> CH2 + CO */
    phi_f = sc[11]*sc[14];
    k_f = 1e-06 * 9e+12;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[14];
    Kc = exp((g_RT[11] + g_RT[14]) - (g_RT[10] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 152: CH2(S) + CO2 <=> CH2 + CO2 */
    phi_f = sc[11]*sc[15];
    k_f = 1e-06 * 7e+12;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[15];
    Kc = exp((g_RT[11] + g_RT[15]) - (g_RT[10] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 153: CH2(S) + CO2 <=> CO + CH2O */
    phi_f = sc[11]*sc[15];
    k_f = 1e-06 * 1.4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[17];
    Kc = exp((g_RT[11] + g_RT[15]) - (g_RT[14] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 154: CH2(S) + C2H6 <=> CH3 + C2H5 */
    phi_f = sc[11]*sc[26];
    k_f = 1e-06 * 4e+13*exp(+276.799/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[25];
    Kc = exp((g_RT[11] + g_RT[26]) - (g_RT[12] + g_RT[25]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[26] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[25] += 1 * qdot;

    /*reaction 155: CH3 + O2 <=> O + CH3O */
    phi_f = sc[12]*sc[3];
    k_f = 1e-06 * 2.675e+13*exp(-14494.2/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[19];
    Kc = exp((g_RT[12] + g_RT[3]) - (g_RT[2] + g_RT[19]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[19] += 1 * qdot;

    /*reaction 156: CH3 + O2 <=> OH + CH2O */
    phi_f = sc[12]*sc[3];
    k_f = 1e-06 * 3.6e+10*exp(-4499.25/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = exp((g_RT[12] + g_RT[3]) - (g_RT[4] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 157: CH3 + H2O2 <=> HO2 + CH4 */
    phi_f = sc[12]*sc[7];
    k_f = 1e-06 * 24500*exp(2.47*tc[0]-2606.95/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = exp((g_RT[12] + g_RT[7]) - (g_RT[6] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 158: 2 CH3 (+M) <=> C2H6 (+M) */
    phi_f = sc[12]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 2.12e+16*exp(-0.97*tc[0]-312.028/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.77e+50*exp(-9.67*tc[0]-3130.35/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.4675*exp(T/-151))+ (0.5325*exp(T/-1038))+ (exp(-4970/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[26];
    Kc = 1.0 / (refC) * exp((2 * g_RT[12]) - (g_RT[26]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 2 * qdot;
    wdot[26] += 1 * qdot;

    /*reaction 159: 2 CH3 <=> H + C2H5 */
    phi_f = sc[12]*sc[12];
    k_f = 1e-06 * 4.99e+12*exp(0.1*tc[0]-5334.68/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[25];
    Kc = exp((2 * g_RT[12]) - (g_RT[1] + g_RT[25]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 2 * qdot;
    wdot[1] += 1 * qdot;
    wdot[25] += 1 * qdot;

    /*reaction 160: CH3 + HCO <=> CH4 + CO */
    phi_f = sc[12]*sc[16];
    k_f = 1e-06 * 2.648e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[14];
    Kc = exp((g_RT[12] + g_RT[16]) - (g_RT[13] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[16] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 161: CH3 + CH2O <=> HCO + CH4 */
    phi_f = sc[12]*sc[17];
    k_f = 1e-06 * 3320*exp(2.81*tc[0]-2949.17/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[13];
    Kc = exp((g_RT[12] + g_RT[17]) - (g_RT[16] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 162: CH3 + CH3OH <=> CH2OH + CH4 */
    phi_f = sc[12]*sc[20];
    k_f = 1e-06 * 3e+07*exp(1.5*tc[0]-5002.52/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[13];
    Kc = exp((g_RT[12] + g_RT[20]) - (g_RT[18] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 163: CH3 + CH3OH <=> CH3O + CH4 */
    phi_f = sc[12]*sc[20];
    k_f = 1e-06 * 1e+07*exp(1.5*tc[0]-5002.52/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[13];
    Kc = exp((g_RT[12] + g_RT[20]) - (g_RT[19] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[19] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 164: CH3 + C2H4 <=> C2H3 + CH4 */
    phi_f = sc[12]*sc[24];
    k_f = 1e-06 * 227000*exp(2*tc[0]-4630.1/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[13];
    Kc = exp((g_RT[12] + g_RT[24]) - (g_RT[23] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[24] -= 1 * qdot;
    wdot[23] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 165: CH3 + C2H6 <=> C2H5 + CH4 */
    phi_f = sc[12]*sc[26];
    k_f = 1e-06 * 6.14e+06*exp(1.74*tc[0]-5259.18/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[13];
    Kc = exp((g_RT[12] + g_RT[26]) - (g_RT[25] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[26] -= 1 * qdot;
    wdot[25] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 166: HCO + H2O <=> H + CO + H2O */
    phi_f = sc[16]*sc[5];
    k_f = 1e-06 * 2.244e+18*exp(-1*tc[0]-8555.61/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14]*sc[5];
    Kc = refC * exp((g_RT[16] + g_RT[5]) - (g_RT[1] + g_RT[14] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 167: HCO + M <=> H + CO + M */
    phi_f = sc[16];
    alpha = mixture + sc[0] + -1*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = 1e-06 * alpha * 1.87e+17*exp(-1*tc[0]-8555.61/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = refC * exp((g_RT[16]) - (g_RT[1] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 168: HCO + O2 <=> HO2 + CO */
    phi_f = sc[16]*sc[3];
    k_f = 1e-06 * 7.6e+12*exp(-201.309/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[14];
    Kc = exp((g_RT[16] + g_RT[3]) - (g_RT[6] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 169: CH2OH + O2 <=> HO2 + CH2O */
    phi_f = sc[18]*sc[3];
    k_f = 1e-06 * 1.8e+13*exp(-452.944/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = exp((g_RT[18] + g_RT[3]) - (g_RT[6] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 170: CH3O + O2 <=> HO2 + CH2O */
    phi_f = sc[19]*sc[3];
    k_f = 1e-06 * 4.28e-13*exp(7.6*tc[0]+1776.55/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = exp((g_RT[19] + g_RT[3]) - (g_RT[6] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 171: C2H + O2 <=> HCO + CO */
    phi_f = sc[21]*sc[3];
    k_f = 1e-06 * 5e+13*exp(-754.907/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[14];
    Kc = exp((g_RT[21] + g_RT[3]) - (g_RT[16] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 172: C2H + H2 <=> H + C2H2 */
    phi_f = sc[21]*sc[0];
    k_f = 1e-06 * 407000*exp(2.4*tc[0]-100.654/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = exp((g_RT[21] + g_RT[0]) - (g_RT[1] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 173: C2H3 + O2 <=> HCO + CH2O */
    phi_f = sc[23]*sc[3];
    k_f = 1e-06 * 3.98e+12*exp(+120.785/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[17];
    Kc = exp((g_RT[23] + g_RT[3]) - (g_RT[16] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[23] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 174: C2H4 (+M) <=> H2 + C2H2 (+M) */
    phi_f = sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1 * 8e+12*exp(0.44*tc[0]-44675.4/tc[1]);
    redP = 1e-12 * alpha / k_f * 7e+50*exp(-9.31*tc[0]-50256.7/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2655*exp(T/-180))+ (0.7345*exp(T/-1035))+ (exp(-5417/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = refC * exp((g_RT[24]) - (g_RT[0] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[24] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 175: C2H5 + O2 <=> HO2 + C2H4 */
    phi_f = sc[25]*sc[3];
    k_f = 1e-06 * 8.4e+11*exp(-1950.18/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[24];
    Kc = exp((g_RT[25] + g_RT[3]) - (g_RT[6] + g_RT[24]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[25] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[24] += 1 * qdot;

    /*reaction 176: HCCO + O2 <=> OH + 2 CO */
    phi_f = sc[27]*sc[3];
    k_f = 1e-06 * 1.6e+12*exp(-429.794/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[14]*sc[14];
    Kc = refC * exp((g_RT[27] + g_RT[3]) - (g_RT[4] + 2 * g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[27] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[14] += 2 * qdot;

    /*reaction 177: 2 HCCO <=> 2 CO + C2H2 */
    phi_f = sc[27]*sc[27];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[14]*sc[22];
    Kc = refC * exp((2 * g_RT[27]) - (2 * g_RT[14] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[27] -= 2 * qdot;
    wdot[14] += 2 * qdot;
    wdot[22] += 1 * qdot;

    return;
}


/*compute the progress rate for each reaction */
void progressRate(double * qdot, double * sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[32];                /*Gibbs free energy */
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
    for (id = 0; id < 32; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*reaction 1: 2 O + M <=> O2 + M */
    phi_f = sc[2]*sc[2];
    alpha = mixture + 1.4*sc[0] + 14.4*sc[5] + sc[13] + 0.75*sc[14] + 2.6*sc[15] + 2*sc[26] + -0.17*sc[31];
    k_f = 1e-12 * alpha * 1.2e+17*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[3];
    Kc = 1.0 / (refC) * exp((2 * g_RT[2]) - (g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[0] = q_f - q_r;

    /*reaction 2: O + H + M <=> OH + M */
    phi_f = sc[2]*sc[1];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-12 * alpha * 5e+17*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[4];
    Kc = 1.0 / (refC) * exp((g_RT[2] + g_RT[1]) - (g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[1] = q_f - q_r;

    /*reaction 3: O + H2 <=> H + OH */
    phi_f = sc[2]*sc[0];
    k_f = 1e-06 * 50000*exp(2.67*tc[0]-3165.58/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[4];
    Kc = exp((g_RT[2] + g_RT[0]) - (g_RT[1] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[2] = q_f - q_r;

    /*reaction 4: O + HO2 <=> OH + O2 */
    phi_f = sc[2]*sc[6];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = exp((g_RT[2] + g_RT[6]) - (g_RT[4] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[3] = q_f - q_r;

    /*reaction 5: O + H2O2 <=> OH + HO2 */
    phi_f = sc[2]*sc[7];
    k_f = 1e-06 * 9.63e+06*exp(2*tc[0]-2013.09/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[6];
    Kc = exp((g_RT[2] + g_RT[7]) - (g_RT[4] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[4] = q_f - q_r;

    /*reaction 6: O + CH <=> H + CO */
    phi_f = sc[2]*sc[9];
    k_f = 1e-06 * 5.7e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = exp((g_RT[2] + g_RT[9]) - (g_RT[1] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[5] = q_f - q_r;

    /*reaction 7: O + CH2 <=> H + HCO */
    phi_f = sc[2]*sc[10];
    k_f = 1e-06 * 8e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = exp((g_RT[2] + g_RT[10]) - (g_RT[1] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[6] = q_f - q_r;

    /*reaction 8: O + CH2(S) <=> H2 + CO */
    phi_f = sc[2]*sc[11];
    k_f = 1e-06 * 1.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[14];
    Kc = exp((g_RT[2] + g_RT[11]) - (g_RT[0] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[7] = q_f - q_r;

    /*reaction 9: O + CH2(S) <=> H + HCO */
    phi_f = sc[2]*sc[11];
    k_f = 1e-06 * 1.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = exp((g_RT[2] + g_RT[11]) - (g_RT[1] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[8] = q_f - q_r;

    /*reaction 10: O + CH3 <=> H + CH2O */
    phi_f = sc[2]*sc[12];
    k_f = 1e-06 * 8.43e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = exp((g_RT[2] + g_RT[12]) - (g_RT[1] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[9] = q_f - q_r;

    /*reaction 11: O + CH4 <=> OH + CH3 */
    phi_f = sc[2]*sc[13];
    k_f = 1e-06 * 1.02e+09*exp(1.5*tc[0]-4328.13/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = exp((g_RT[2] + g_RT[13]) - (g_RT[4] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[10] = q_f - q_r;

    /*reaction 12: O + CO + M <=> CO2 + M */
    phi_f = sc[2]*sc[14];
    alpha = mixture + sc[0] + 5*sc[3] + 5*sc[5] + sc[13] + 0.5*sc[14] + 2.5*sc[15] + 2*sc[26] + -0.5*sc[31];
    k_f = 1e-12 * alpha * 6.02e+14*exp(-1509.81/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15];
    Kc = 1.0 / (refC) * exp((g_RT[2] + g_RT[14]) - (g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[11] = q_f - q_r;

    /*reaction 13: O + HCO <=> OH + CO */
    phi_f = sc[2]*sc[16];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[14];
    Kc = exp((g_RT[2] + g_RT[16]) - (g_RT[4] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[12] = q_f - q_r;

    /*reaction 14: O + HCO <=> H + CO2 */
    phi_f = sc[2]*sc[16];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[15];
    Kc = exp((g_RT[2] + g_RT[16]) - (g_RT[1] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[13] = q_f - q_r;

    /*reaction 15: O + CH2O <=> OH + HCO */
    phi_f = sc[2]*sc[17];
    k_f = 1e-06 * 3.9e+13*exp(-1781.58/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[16];
    Kc = exp((g_RT[2] + g_RT[17]) - (g_RT[4] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[14] = q_f - q_r;

    /*reaction 16: O + CH2OH <=> OH + CH2O */
    phi_f = sc[2]*sc[18];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = exp((g_RT[2] + g_RT[18]) - (g_RT[4] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[15] = q_f - q_r;

    /*reaction 17: O + CH3O <=> OH + CH2O */
    phi_f = sc[2]*sc[19];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = exp((g_RT[2] + g_RT[19]) - (g_RT[4] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[16] = q_f - q_r;

    /*reaction 18: O + CH3OH <=> OH + CH2OH */
    phi_f = sc[2]*sc[20];
    k_f = 1e-06 * 388000*exp(2.5*tc[0]-1560.14/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[18];
    Kc = exp((g_RT[2] + g_RT[20]) - (g_RT[4] + g_RT[18]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[17] = q_f - q_r;

    /*reaction 19: O + CH3OH <=> OH + CH3O */
    phi_f = sc[2]*sc[20];
    k_f = 1e-06 * 130000*exp(2.5*tc[0]-2516.36/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[19];
    Kc = exp((g_RT[2] + g_RT[20]) - (g_RT[4] + g_RT[19]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[18] = q_f - q_r;

    /*reaction 20: O + C2H <=> CH + CO */
    phi_f = sc[2]*sc[21];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[14];
    Kc = exp((g_RT[2] + g_RT[21]) - (g_RT[9] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[19] = q_f - q_r;

    /*reaction 21: O + C2H2 <=> H + HCCO */
    phi_f = sc[2]*sc[22];
    k_f = 1e-06 * 1.02e+07*exp(2*tc[0]-956.215/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[27];
    Kc = exp((g_RT[2] + g_RT[22]) - (g_RT[1] + g_RT[27]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[20] = q_f - q_r;

    /*reaction 22: O + C2H2 <=> OH + C2H */
    phi_f = sc[2]*sc[22];
    k_f = 1e-06 * 4.6e+19*exp(-1.41*tc[0]-14569.7/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[21];
    Kc = exp((g_RT[2] + g_RT[22]) - (g_RT[4] + g_RT[21]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[21] = q_f - q_r;

    /*reaction 23: O + C2H2 <=> CO + CH2 */
    phi_f = sc[2]*sc[22];
    k_f = 1e-06 * 1.02e+07*exp(2*tc[0]-956.215/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[10];
    Kc = exp((g_RT[2] + g_RT[22]) - (g_RT[14] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[22] = q_f - q_r;

    /*reaction 24: O + C2H3 <=> H + CH2CO */
    phi_f = sc[2]*sc[23];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = exp((g_RT[2] + g_RT[23]) - (g_RT[1] + g_RT[28]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[23] = q_f - q_r;

    /*reaction 25: O + C2H4 <=> CH3 + HCO */
    phi_f = sc[2]*sc[24];
    k_f = 1e-06 * 1.92e+07*exp(1.83*tc[0]-110.72/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[16];
    Kc = exp((g_RT[2] + g_RT[24]) - (g_RT[12] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[24] = q_f - q_r;

    /*reaction 26: O + C2H5 <=> CH3 + CH2O */
    phi_f = sc[2]*sc[25];
    k_f = 1e-06 * 1.32e+14;
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[17];
    Kc = exp((g_RT[2] + g_RT[25]) - (g_RT[12] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[25] = q_f - q_r;

    /*reaction 27: O + C2H6 <=> OH + C2H5 */
    phi_f = sc[2]*sc[26];
    k_f = 1e-06 * 8.98e+07*exp(1.92*tc[0]-2863.61/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[25];
    Kc = exp((g_RT[2] + g_RT[26]) - (g_RT[4] + g_RT[25]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[26] = q_f - q_r;

    /*reaction 28: O + HCCO <=> H + 2 CO */
    phi_f = sc[2]*sc[27];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14]*sc[14];
    Kc = refC * exp((g_RT[2] + g_RT[27]) - (g_RT[1] + 2 * g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[27] = q_f - q_r;

    /*reaction 29: O + CH2CO <=> OH + HCCO */
    phi_f = sc[2]*sc[28];
    k_f = 1e-06 * 1e+13*exp(-4026.17/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[27];
    Kc = exp((g_RT[2] + g_RT[28]) - (g_RT[4] + g_RT[27]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[28] = q_f - q_r;

    /*reaction 30: O + CH2CO <=> CH2 + CO2 */
    phi_f = sc[2]*sc[28];
    k_f = 1e-06 * 1.75e+12*exp(-679.416/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[15];
    Kc = exp((g_RT[2] + g_RT[28]) - (g_RT[10] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[29] = q_f - q_r;

    /*reaction 31: O2 + CO <=> O + CO2 */
    phi_f = sc[3]*sc[14];
    k_f = 1e-06 * 2.5e+12*exp(-24056.4/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[15];
    Kc = exp((g_RT[3] + g_RT[14]) - (g_RT[2] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[30] = q_f - q_r;

    /*reaction 32: O2 + CH2O <=> HO2 + HCO */
    phi_f = sc[3]*sc[17];
    k_f = 1e-06 * 1e+14*exp(-20130.9/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[16];
    Kc = exp((g_RT[3] + g_RT[17]) - (g_RT[6] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[31] = q_f - q_r;

    /*reaction 33: H + O2 + M <=> HO2 + M */
    phi_f = sc[1]*sc[3];
    alpha = mixture + -1*sc[3] + -1*sc[5] + -0.25*sc[14] + 0.5*sc[15] + 0.5*sc[26] + -1*sc[30] + -1*sc[31];
    k_f = 1e-12 * alpha * 2.8e+18*exp(-0.86*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[32] = q_f - q_r;

    /*reaction 34: H + 2 O2 <=> HO2 + O2 */
    phi_f = sc[1]*sc[3]*sc[3];
    k_f = 1e-12 * 3e+20*exp(-1.72*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = 1.0 / (refC) * exp((g_RT[1] + 2 * g_RT[3]) - (g_RT[6] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[33] = q_f - q_r;

    /*reaction 35: H + O2 + H2O <=> HO2 + H2O */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = 1e-12 * 9.38e+18*exp(-0.76*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[34] = q_f - q_r;

    /*reaction 36: H + O2 + N2 <=> HO2 + N2 */
    phi_f = sc[1]*sc[3]*sc[30];
    k_f = 1e-12 * 3.75e+20*exp(-1.72*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[30];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[30]) - (g_RT[6] + g_RT[30]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[35] = q_f - q_r;

    /*reaction 37: H + O2 + AR <=> HO2 + AR */
    phi_f = sc[1]*sc[3]*sc[31];
    k_f = 1e-12 * 7e+17*exp(-0.8*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[31];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[31]) - (g_RT[6] + g_RT[31]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[36] = q_f - q_r;

    /*reaction 38: H + O2 <=> O + OH */
    phi_f = sc[1]*sc[3];
    k_f = 1e-06 * 8.3e+13*exp(-7253.65/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[37] = q_f - q_r;

    /*reaction 39: 2 H + M <=> H2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[0] + -1*sc[5] + sc[13] + -1*sc[15] + 2*sc[26] + -0.37*sc[31];
    k_f = 1e-12 * alpha * 1e+18*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[38] = q_f - q_r;

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    phi_f = sc[1]*sc[1]*sc[0];
    k_f = 1e-12 * 9e+16*exp(-0.6*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[0];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[39] = q_f - q_r;

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    phi_f = sc[1]*sc[1]*sc[5];
    k_f = 1e-12 * 6e+19*exp(-1.25*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[5];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[40] = q_f - q_r;

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    phi_f = sc[1]*sc[1]*sc[15];
    k_f = 1e-12 * 5.5e+20*exp(-2*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[15];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[15]) - (g_RT[0] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[41] = q_f - q_r;

    /*reaction 43: H + OH + M <=> H2O + M */
    phi_f = sc[1]*sc[4];
    alpha = mixture + -0.27*sc[0] + 2.65*sc[5] + sc[13] + 2*sc[26] + -0.62*sc[31];
    k_f = 1e-12 * alpha * 2.2e+22*exp(-2*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[42] = q_f - q_r;

    /*reaction 44: H + HO2 <=> O + H2O */
    phi_f = sc[1]*sc[6];
    k_f = 1e-06 * 3.97e+12*exp(-337.695/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = exp((g_RT[1] + g_RT[6]) - (g_RT[2] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[43] = q_f - q_r;

    /*reaction 45: H + HO2 <=> O2 + H2 */
    phi_f = sc[1]*sc[6];
    k_f = 1e-06 * 2.8e+13*exp(-537.494/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[44] = q_f - q_r;

    /*reaction 46: H + HO2 <=> 2 OH */
    phi_f = sc[1]*sc[6];
    k_f = 1e-06 * 1.34e+14*exp(-319.577/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[45] = q_f - q_r;

    /*reaction 47: H + H2O2 <=> HO2 + H2 */
    phi_f = sc[1]*sc[7];
    k_f = 1e-06 * 1.21e+07*exp(2*tc[0]-2617.01/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = exp((g_RT[1] + g_RT[7]) - (g_RT[6] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[46] = q_f - q_r;

    /*reaction 48: H + H2O2 <=> OH + H2O */
    phi_f = sc[1]*sc[7];
    k_f = 1e-06 * 1e+13*exp(-1811.78/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[5];
    Kc = exp((g_RT[1] + g_RT[7]) - (g_RT[4] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[47] = q_f - q_r;

    /*reaction 49: H + CH <=> C + H2 */
    phi_f = sc[1]*sc[9];
    k_f = 1e-06 * 1.1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = exp((g_RT[1] + g_RT[9]) - (g_RT[8] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[48] = q_f - q_r;

    /*reaction 50: H + CH2 (+M) <=> CH3 (+M) */
    phi_f = sc[1]*sc[10];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 2.5e+16*exp(-0.8*tc[0]);
    redP = 1e-12 * alpha / k_f * 3.2e+27*exp(-3.14*tc[0]-619.024/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.32*exp(T/-78))+ (0.68*exp(T/-1995))+ (exp(-5590/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[12];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[49] = q_f - q_r;

    /*reaction 51: H + CH2(S) <=> CH + H2 */
    phi_f = sc[1]*sc[11];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = exp((g_RT[1] + g_RT[11]) - (g_RT[9] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[50] = q_f - q_r;

    /*reaction 52: H + CH3 (+M) <=> CH4 (+M) */
    phi_f = sc[1]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 1.27e+16*exp(-0.63*tc[0]-192.753/tc[1]);
    redP = 1e-12 * alpha / k_f * 2.477e+33*exp(-4.76*tc[0]-1227.98/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.217*exp(T/-74))+ (0.783*exp(T/-2941))+ (exp(-6964/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[13];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[12]) - (g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[51] = q_f - q_r;

    /*reaction 53: H + CH4 <=> CH3 + H2 */
    phi_f = sc[1]*sc[13];
    k_f = 1e-06 * 6.6e+08*exp(1.62*tc[0]-5455.46/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[0];
    Kc = exp((g_RT[1] + g_RT[13]) - (g_RT[12] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[52] = q_f - q_r;

    /*reaction 54: H + HCO (+M) <=> CH2O (+M) */
    phi_f = sc[1]*sc[16];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 1.09e+12*exp(0.48*tc[0]+130.851/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.35e+24*exp(-2.57*tc[0]-717.162/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2176*exp(T/-271))+ (0.7824*exp(T/-2755))+ (exp(-6570/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[17];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[16]) - (g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[53] = q_f - q_r;

    /*reaction 55: H + HCO <=> H2 + CO */
    phi_f = sc[1]*sc[16];
    k_f = 1e-06 * 7.34e+13;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[14];
    Kc = exp((g_RT[1] + g_RT[16]) - (g_RT[0] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[54] = q_f - q_r;

    /*reaction 56: H + CH2O (+M) <=> CH2OH (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1811.78/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.27e+32*exp(-4.82*tc[0]-3286.36/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2813*exp(T/-103))+ (0.7187*exp(T/-1291))+ (exp(-4160/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[18];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[18]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[55] = q_f - q_r;

    /*reaction 57: H + CH2O (+M) <=> CH3O (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1308.51/tc[1]);
    redP = 1e-12 * alpha / k_f * 2.2e+30*exp(-4.8*tc[0]-2798.19/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.242*exp(T/-94))+ (0.758*exp(T/-1555))+ (exp(-4200/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[19];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[19]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[56] = q_f - q_r;

    /*reaction 58: H + CH2O <=> HCO + H2 */
    phi_f = sc[1]*sc[17];
    k_f = 1e-06 * 2.3e+10*exp(1.05*tc[0]-1648.21/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[0];
    Kc = exp((g_RT[1] + g_RT[17]) - (g_RT[16] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[57] = q_f - q_r;

    /*reaction 59: H + CH2OH (+M) <=> CH3OH (+M) */
    phi_f = sc[1]*sc[18];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = 1e-06 * 1.8e+13;
    redP = 1e-12 * alpha / k_f * 3e+31*exp(-4.8*tc[0]-1660.8/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2321*exp(T/-338))+ (0.7679*exp(T/-1812))+ (exp(-5081/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[20];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[18]) - (g_RT[20]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[58] = q_f - q_r;

    /*reaction 60: H + CH2OH <=> H2 + CH2O */
    phi_f = sc[1]*sc[18];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[17];
    Kc = exp((g_RT[1] + g_RT[18]) - (g_RT[0] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[59] = q_f - q_r;

    /*reaction 61: H + CH2OH <=> OH + CH3 */
    phi_f = sc[1]*sc[18];
    k_f = 1e-06 * 1.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = exp((g_RT[1] + g_RT[18]) - (g_RT[4] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[60] = q_f - q_r;

    /*reaction 62: H + CH2OH <=> CH2(S) + H2O */
    phi_f = sc[1]*sc[18];
    k_f = 1e-06 * 6e+12;
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = exp((g_RT[1] + g_RT[18]) - (g_RT[11] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[61] = q_f - q_r;

    /*reaction 63: H + CH3O (+M) <=> CH3OH (+M) */
    phi_f = sc[1]*sc[19];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = 1e-06 * 5e+13;
    redP = 1e-12 * alpha / k_f * 8.6e+28*exp(-4*tc[0]-1522.4/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.1098*exp(T/-144))+ (0.8902*exp(T/-2838))+ (exp(-45569/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[20];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[19]) - (g_RT[20]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[62] = q_f - q_r;

    /*reaction 64: H + CH3O <=> H + CH2OH */
    phi_f = sc[1]*sc[19];
    k_f = 1e-06 * 3.4e+06*exp(1.6*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[18];
    Kc = exp((g_RT[1] + g_RT[19]) - (g_RT[1] + g_RT[18]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[63] = q_f - q_r;

    /*reaction 65: H + CH3O <=> H2 + CH2O */
    phi_f = sc[1]*sc[19];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[17];
    Kc = exp((g_RT[1] + g_RT[19]) - (g_RT[0] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[64] = q_f - q_r;

    /*reaction 66: H + CH3O <=> OH + CH3 */
    phi_f = sc[1]*sc[19];
    k_f = 1e-06 * 3.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = exp((g_RT[1] + g_RT[19]) - (g_RT[4] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[65] = q_f - q_r;

    /*reaction 67: H + CH3O <=> CH2(S) + H2O */
    phi_f = sc[1]*sc[19];
    k_f = 1e-06 * 1.6e+13;
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = exp((g_RT[1] + g_RT[19]) - (g_RT[11] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[66] = q_f - q_r;

    /*reaction 68: H + CH3OH <=> CH2OH + H2 */
    phi_f = sc[1]*sc[20];
    k_f = 1e-06 * 1.7e+07*exp(2.1*tc[0]-2450.93/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[0];
    Kc = exp((g_RT[1] + g_RT[20]) - (g_RT[18] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[67] = q_f - q_r;

    /*reaction 69: H + CH3OH <=> CH3O + H2 */
    phi_f = sc[1]*sc[20];
    k_f = 1e-06 * 4.2e+06*exp(2.1*tc[0]-2450.93/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[0];
    Kc = exp((g_RT[1] + g_RT[20]) - (g_RT[19] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[68] = q_f - q_r;

    /*reaction 70: H + C2H (+M) <=> C2H2 (+M) */
    phi_f = sc[1]*sc[21];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 1e+17*exp(-1*tc[0]);
    redP = 1e-12 * alpha / k_f * 3.75e+33*exp(-4.8*tc[0]-956.215/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.3536*exp(T/-132))+ (0.6464*exp(T/-1315))+ (exp(-5566/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[22];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[21]) - (g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[69] = q_f - q_r;

    /*reaction 71: H + C2H2 (+M) <=> C2H3 (+M) */
    phi_f = sc[1]*sc[22];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 5.6e+12*exp(-1207.85/tc[1]);
    redP = 1e-12 * alpha / k_f * 3.8e+40*exp(-7.27*tc[0]-3633.62/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2493*exp(T/-98.5))+ (0.7507*exp(T/-1302))+ (exp(-4167/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[23];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[22]) - (g_RT[23]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[70] = q_f - q_r;

    /*reaction 72: H + C2H3 (+M) <=> C2H4 (+M) */
    phi_f = sc[1]*sc[23];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 6.08e+12*exp(0.27*tc[0]-140.916/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.4e+30*exp(-3.86*tc[0]-1670.86/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.218*exp(T/-207.5))+ (0.782*exp(T/-2663))+ (exp(-6095/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[24];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[23]) - (g_RT[24]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[71] = q_f - q_r;

    /*reaction 73: H + C2H3 <=> H2 + C2H2 */
    phi_f = sc[1]*sc[23];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = exp((g_RT[1] + g_RT[23]) - (g_RT[0] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[72] = q_f - q_r;

    /*reaction 74: H + C2H4 (+M) <=> C2H5 (+M) */
    phi_f = sc[1]*sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 1.08e+12*exp(0.454*tc[0]-915.954/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.2e+42*exp(-7.62*tc[0]-3507.8/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.0247*exp(T/-210))+ (0.9753*exp(T/-984))+ (exp(-4374/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[25];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[24]) - (g_RT[25]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[73] = q_f - q_r;

    /*reaction 75: H + C2H4 <=> C2H3 + H2 */
    phi_f = sc[1]*sc[24];
    k_f = 1e-06 * 1.325e+06*exp(2.53*tc[0]-6160.04/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[0];
    Kc = exp((g_RT[1] + g_RT[24]) - (g_RT[23] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[74] = q_f - q_r;

    /*reaction 76: H + C2H5 (+M) <=> C2H6 (+M) */
    phi_f = sc[1]*sc[25];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 5.21e+17*exp(-0.99*tc[0]-795.169/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.99e+41*exp(-7.08*tc[0]-3364.37/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.1578*exp(T/-125))+ (0.8422*exp(T/-2219))+ (exp(-6882/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[26];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[25]) - (g_RT[26]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[75] = q_f - q_r;

    /*reaction 77: H + C2H5 <=> H2 + C2H4 */
    phi_f = sc[1]*sc[25];
    k_f = 1e-06 * 2e+12;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[24];
    Kc = exp((g_RT[1] + g_RT[25]) - (g_RT[0] + g_RT[24]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[76] = q_f - q_r;

    /*reaction 78: H + C2H6 <=> C2H5 + H2 */
    phi_f = sc[1]*sc[26];
    k_f = 1e-06 * 1.15e+08*exp(1.9*tc[0]-3789.63/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[0];
    Kc = exp((g_RT[1] + g_RT[26]) - (g_RT[25] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[77] = q_f - q_r;

    /*reaction 79: H + HCCO <=> CH2(S) + CO */
    phi_f = sc[1]*sc[27];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[14];
    Kc = exp((g_RT[1] + g_RT[27]) - (g_RT[11] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[78] = q_f - q_r;

    /*reaction 80: H + CH2CO <=> HCCO + H2 */
    phi_f = sc[1]*sc[28];
    k_f = 1e-06 * 5e+13*exp(-4026.17/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[0];
    Kc = exp((g_RT[1] + g_RT[28]) - (g_RT[27] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[79] = q_f - q_r;

    /*reaction 81: H + CH2CO <=> CH3 + CO */
    phi_f = sc[1]*sc[28];
    k_f = 1e-06 * 1.13e+13*exp(-1725.21/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[14];
    Kc = exp((g_RT[1] + g_RT[28]) - (g_RT[12] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[80] = q_f - q_r;

    /*reaction 82: H + HCCOH <=> H + CH2CO */
    phi_f = sc[1]*sc[29];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = exp((g_RT[1] + g_RT[29]) - (g_RT[1] + g_RT[28]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[81] = q_f - q_r;

    /*reaction 83: H2 + CO (+M) <=> CH2O (+M) */
    phi_f = sc[0]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 4.3e+07*exp(1.5*tc[0]-40060.4/tc[1]);
    redP = 1e-12 * alpha / k_f * 5.07e+27*exp(-3.42*tc[0]-42450.9/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.068*exp(T/-197))+ (0.932*exp(T/-1540))+ (exp(-10300/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[17];
    Kc = 1.0 / (refC) * exp((g_RT[0] + g_RT[14]) - (g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[82] = q_f - q_r;

    /*reaction 84: OH + H2 <=> H + H2O */
    phi_f = sc[4]*sc[0];
    k_f = 1e-06 * 2.16e+08*exp(1.51*tc[0]-1726.22/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    Kc = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[83] = q_f - q_r;

    /*reaction 85: 2 OH (+M) <=> H2O2 (+M) */
    phi_f = sc[4]*sc[4];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 7.4e+13*exp(-0.37*tc[0]);
    redP = 1e-12 * alpha / k_f * 2.3e+18*exp(-0.9*tc[0]+855.561/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2654*exp(T/-94))+ (0.7346*exp(T/-1756))+ (exp(-5182/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[7];
    Kc = 1.0 / (refC) * exp((2 * g_RT[4]) - (g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[84] = q_f - q_r;

    /*reaction 86: 2 OH <=> O + H2O */
    phi_f = sc[4]*sc[4];
    k_f = 1e-06 * 35700*exp(2.4*tc[0]+1061.9/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[85] = q_f - q_r;

    /*reaction 87: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = 1e-06 * 2.9e+13*exp(+251.636/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[86] = q_f - q_r;

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    phi_f = sc[4]*sc[7];
    k_f = 1e-06 * 1.75e+12*exp(-161.047/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[87] = q_f - q_r;

    /*reaction 89: OH + H2O2 <=> HO2 + H2O */
    phi_f = sc[4]*sc[7];
    k_f = 1e-06 * 5.8e+14*exp(-4811.27/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[88] = q_f - q_r;

    /*reaction 90: OH + C <=> H + CO */
    phi_f = sc[4]*sc[8];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = exp((g_RT[4] + g_RT[8]) - (g_RT[1] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[89] = q_f - q_r;

    /*reaction 91: OH + CH <=> H + HCO */
    phi_f = sc[4]*sc[9];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = exp((g_RT[4] + g_RT[9]) - (g_RT[1] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[90] = q_f - q_r;

    /*reaction 92: OH + CH2 <=> H + CH2O */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = exp((g_RT[4] + g_RT[10]) - (g_RT[1] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[91] = q_f - q_r;

    /*reaction 93: OH + CH2 <=> CH + H2O */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 1.13e+07*exp(2*tc[0]-1509.81/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[5];
    Kc = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[92] = q_f - q_r;

    /*reaction 94: OH + CH2(S) <=> H + CH2O */
    phi_f = sc[4]*sc[11];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = exp((g_RT[4] + g_RT[11]) - (g_RT[1] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[93] = q_f - q_r;

    /*reaction 95: OH + CH3 (+M) <=> CH3OH (+M) */
    phi_f = sc[4]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = 1e-06 * 6.3e+13;
    redP = 1e-12 * alpha / k_f * 2.7e+38*exp(-6.3*tc[0]-1560.14/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.7895*exp(T/-83.5))+ (0.2105*exp(T/-5398))+ (exp(-8370/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[20];
    Kc = 1.0 / (refC) * exp((g_RT[4] + g_RT[12]) - (g_RT[20]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[94] = q_f - q_r;

    /*reaction 96: OH + CH3 <=> CH2 + H2O */
    phi_f = sc[4]*sc[12];
    k_f = 1e-06 * 5.6e+07*exp(1.6*tc[0]-2727.73/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[5];
    Kc = exp((g_RT[4] + g_RT[12]) - (g_RT[10] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[95] = q_f - q_r;

    /*reaction 97: OH + CH3 <=> CH2(S) + H2O */
    phi_f = sc[4]*sc[12];
    k_f = 1e-06 * 2.501e+13;
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = exp((g_RT[4] + g_RT[12]) - (g_RT[11] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[96] = q_f - q_r;

    /*reaction 98: OH + CH4 <=> CH3 + H2O */
    phi_f = sc[4]*sc[13];
    k_f = 1e-06 * 1e+08*exp(1.6*tc[0]-1570.21/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[5];
    Kc = exp((g_RT[4] + g_RT[13]) - (g_RT[12] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[97] = q_f - q_r;

    /*reaction 99: OH + CO <=> H + CO2 */
    phi_f = sc[4]*sc[14];
    k_f = 1e-06 * 4.76e+07*exp(1.228*tc[0]-35.229/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[15];
    Kc = exp((g_RT[4] + g_RT[14]) - (g_RT[1] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[98] = q_f - q_r;

    /*reaction 100: OH + HCO <=> H2O + CO */
    phi_f = sc[4]*sc[16];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[14];
    Kc = exp((g_RT[4] + g_RT[16]) - (g_RT[5] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[99] = q_f - q_r;

    /*reaction 101: OH + CH2O <=> HCO + H2O */
    phi_f = sc[4]*sc[17];
    k_f = 1e-06 * 3.43e+09*exp(1.18*tc[0]+224.962/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[5];
    Kc = exp((g_RT[4] + g_RT[17]) - (g_RT[16] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[100] = q_f - q_r;

    /*reaction 102: OH + CH2OH <=> H2O + CH2O */
    phi_f = sc[4]*sc[18];
    k_f = 1e-06 * 5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[17];
    Kc = exp((g_RT[4] + g_RT[18]) - (g_RT[5] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[101] = q_f - q_r;

    /*reaction 103: OH + CH3O <=> H2O + CH2O */
    phi_f = sc[4]*sc[19];
    k_f = 1e-06 * 5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[17];
    Kc = exp((g_RT[4] + g_RT[19]) - (g_RT[5] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[102] = q_f - q_r;

    /*reaction 104: OH + CH3OH <=> CH2OH + H2O */
    phi_f = sc[4]*sc[20];
    k_f = 1e-06 * 1.44e+06*exp(2*tc[0]+422.748/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[5];
    Kc = exp((g_RT[4] + g_RT[20]) - (g_RT[18] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[103] = q_f - q_r;

    /*reaction 105: OH + CH3OH <=> CH3O + H2O */
    phi_f = sc[4]*sc[20];
    k_f = 1e-06 * 6.3e+06*exp(2*tc[0]-754.907/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[5];
    Kc = exp((g_RT[4] + g_RT[20]) - (g_RT[19] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[104] = q_f - q_r;

    /*reaction 106: OH + C2H <=> H + HCCO */
    phi_f = sc[4]*sc[21];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[27];
    Kc = exp((g_RT[4] + g_RT[21]) - (g_RT[1] + g_RT[27]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[105] = q_f - q_r;

    /*reaction 107: OH + C2H2 <=> H + CH2CO */
    phi_f = sc[4]*sc[22];
    k_f = 1e-06 * 0.000218*exp(4.5*tc[0]+503.271/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[28]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[106] = q_f - q_r;

    /*reaction 108: OH + C2H2 <=> H + HCCOH */
    phi_f = sc[4]*sc[22];
    k_f = 1e-06 * 504000*exp(2.3*tc[0]-6794.16/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[29];
    Kc = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[29]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[107] = q_f - q_r;

    /*reaction 109: OH + C2H2 <=> C2H + H2O */
    phi_f = sc[4]*sc[22];
    k_f = 1e-06 * 3.37e+07*exp(2*tc[0]-7045.8/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[5];
    Kc = exp((g_RT[4] + g_RT[22]) - (g_RT[21] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[108] = q_f - q_r;

    /*reaction 110: OH + C2H2 <=> CH3 + CO */
    phi_f = sc[4]*sc[22];
    k_f = 1e-06 * 0.000483*exp(4*tc[0]+1006.54/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[14];
    Kc = exp((g_RT[4] + g_RT[22]) - (g_RT[12] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[109] = q_f - q_r;

    /*reaction 111: OH + C2H3 <=> H2O + C2H2 */
    phi_f = sc[4]*sc[23];
    k_f = 1e-06 * 5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[22];
    Kc = exp((g_RT[4] + g_RT[23]) - (g_RT[5] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[110] = q_f - q_r;

    /*reaction 112: OH + C2H4 <=> C2H3 + H2O */
    phi_f = sc[4]*sc[24];
    k_f = 1e-06 * 3.6e+06*exp(2*tc[0]-1258.18/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[5];
    Kc = exp((g_RT[4] + g_RT[24]) - (g_RT[23] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[111] = q_f - q_r;

    /*reaction 113: OH + C2H6 <=> C2H5 + H2O */
    phi_f = sc[4]*sc[26];
    k_f = 1e-06 * 3.54e+06*exp(2.12*tc[0]-437.846/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[5];
    Kc = exp((g_RT[4] + g_RT[26]) - (g_RT[25] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[112] = q_f - q_r;

    /*reaction 114: OH + CH2CO <=> HCCO + H2O */
    phi_f = sc[4]*sc[28];
    k_f = 1e-06 * 7.5e+12*exp(-1006.54/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[5];
    Kc = exp((g_RT[4] + g_RT[28]) - (g_RT[27] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[113] = q_f - q_r;

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    phi_f = sc[6]*sc[6];
    k_f = 1e-06 * 1.3e+11*exp(+820.332/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[114] = q_f - q_r;

    /*reaction 116: 2 HO2 <=> O2 + H2O2 */
    phi_f = sc[6]*sc[6];
    k_f = 1e-06 * 4.2e+14*exp(-6039.26/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[115] = q_f - q_r;

    /*reaction 117: HO2 + CH2 <=> OH + CH2O */
    phi_f = sc[6]*sc[10];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = exp((g_RT[6] + g_RT[10]) - (g_RT[4] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[116] = q_f - q_r;

    /*reaction 118: HO2 + CH3 <=> O2 + CH4 */
    phi_f = sc[6]*sc[12];
    k_f = 1e-06 * 1e+12;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[13];
    Kc = exp((g_RT[6] + g_RT[12]) - (g_RT[3] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[117] = q_f - q_r;

    /*reaction 119: HO2 + CH3 <=> OH + CH3O */
    phi_f = sc[6]*sc[12];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[19];
    Kc = exp((g_RT[6] + g_RT[12]) - (g_RT[4] + g_RT[19]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[118] = q_f - q_r;

    /*reaction 120: HO2 + CO <=> OH + CO2 */
    phi_f = sc[6]*sc[14];
    k_f = 1e-06 * 1.5e+14*exp(-11877.2/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[15];
    Kc = exp((g_RT[6] + g_RT[14]) - (g_RT[4] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[119] = q_f - q_r;

    /*reaction 121: HO2 + CH2O <=> HCO + H2O2 */
    phi_f = sc[6]*sc[17];
    k_f = 1e-06 * 1e+12*exp(-4026.17/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[7];
    Kc = exp((g_RT[6] + g_RT[17]) - (g_RT[16] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[120] = q_f - q_r;

    /*reaction 122: C + O2 <=> O + CO */
    phi_f = sc[8]*sc[3];
    k_f = 1e-06 * 5.8e+13*exp(-289.884/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[14];
    Kc = exp((g_RT[8] + g_RT[3]) - (g_RT[2] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[121] = q_f - q_r;

    /*reaction 123: C + CH2 <=> H + C2H */
    phi_f = sc[8]*sc[10];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[21];
    Kc = exp((g_RT[8] + g_RT[10]) - (g_RT[1] + g_RT[21]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[122] = q_f - q_r;

    /*reaction 124: C + CH3 <=> H + C2H2 */
    phi_f = sc[8]*sc[12];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = exp((g_RT[8] + g_RT[12]) - (g_RT[1] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[123] = q_f - q_r;

    /*reaction 125: CH + O2 <=> O + HCO */
    phi_f = sc[9]*sc[3];
    k_f = 1e-06 * 3.3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[16];
    Kc = exp((g_RT[9] + g_RT[3]) - (g_RT[2] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[124] = q_f - q_r;

    /*reaction 126: CH + H2 <=> H + CH2 */
    phi_f = sc[9]*sc[0];
    k_f = 1e-06 * 1.107e+08*exp(1.79*tc[0]-840.463/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[10];
    Kc = exp((g_RT[9] + g_RT[0]) - (g_RT[1] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[125] = q_f - q_r;

    /*reaction 127: CH + H2O <=> H + CH2O */
    phi_f = sc[9]*sc[5];
    k_f = 1e-06 * 5.71e+12*exp(+379.97/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = exp((g_RT[9] + g_RT[5]) - (g_RT[1] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[126] = q_f - q_r;

    /*reaction 128: CH + CH2 <=> H + C2H2 */
    phi_f = sc[9]*sc[10];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = exp((g_RT[9] + g_RT[10]) - (g_RT[1] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[127] = q_f - q_r;

    /*reaction 129: CH + CH3 <=> H + C2H3 */
    phi_f = sc[9]*sc[12];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[23];
    Kc = exp((g_RT[9] + g_RT[12]) - (g_RT[1] + g_RT[23]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[128] = q_f - q_r;

    /*reaction 130: CH + CH4 <=> H + C2H4 */
    phi_f = sc[9]*sc[13];
    k_f = 1e-06 * 6e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = exp((g_RT[9] + g_RT[13]) - (g_RT[1] + g_RT[24]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[129] = q_f - q_r;

    /*reaction 131: CH + CO (+M) <=> HCCO (+M) */
    phi_f = sc[9]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 5e+13;
    redP = 1e-12 * alpha / k_f * 2.69e+28*exp(-3.74*tc[0]-974.333/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.4243*exp(T/-237))+ (0.5757*exp(T/-1652))+ (exp(-5069/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[27];
    Kc = 1.0 / (refC) * exp((g_RT[9] + g_RT[14]) - (g_RT[27]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[130] = q_f - q_r;

    /*reaction 132: CH + CO2 <=> HCO + CO */
    phi_f = sc[9]*sc[15];
    k_f = 1e-06 * 3.4e+12*exp(-347.257/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[14];
    Kc = exp((g_RT[9] + g_RT[15]) - (g_RT[16] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[131] = q_f - q_r;

    /*reaction 133: CH + CH2O <=> H + CH2CO */
    phi_f = sc[9]*sc[17];
    k_f = 1e-06 * 9.46e+13*exp(+259.185/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = exp((g_RT[9] + g_RT[17]) - (g_RT[1] + g_RT[28]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[132] = q_f - q_r;

    /*reaction 134: CH + HCCO <=> CO + C2H2 */
    phi_f = sc[9]*sc[27];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[22];
    Kc = exp((g_RT[9] + g_RT[27]) - (g_RT[14] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[133] = q_f - q_r;

    /*reaction 135: CH2 + O2 <=> OH + HCO */
    phi_f = sc[10]*sc[3];
    k_f = 1e-06 * 1.32e+13*exp(-754.907/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[16];
    Kc = exp((g_RT[10] + g_RT[3]) - (g_RT[4] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[134] = q_f - q_r;

    /*reaction 136: CH2 + H2 <=> H + CH3 */
    phi_f = sc[10]*sc[0];
    k_f = 1e-06 * 500000*exp(2*tc[0]-3638.65/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = exp((g_RT[10] + g_RT[0]) - (g_RT[1] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[135] = q_f - q_r;

    /*reaction 137: 2 CH2 <=> H2 + C2H2 */
    phi_f = sc[10]*sc[10];
    k_f = 1e-06 * 3.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = exp((2 * g_RT[10]) - (g_RT[0] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[136] = q_f - q_r;

    /*reaction 138: CH2 + CH3 <=> H + C2H4 */
    phi_f = sc[10]*sc[12];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = exp((g_RT[10] + g_RT[12]) - (g_RT[1] + g_RT[24]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[137] = q_f - q_r;

    /*reaction 139: CH2 + CH4 <=> 2 CH3 */
    phi_f = sc[10]*sc[13];
    k_f = 1e-06 * 2.46e+06*exp(2*tc[0]-4162.05/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[12];
    Kc = exp((g_RT[10] + g_RT[13]) - (2 * g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[138] = q_f - q_r;

    /*reaction 140: CH2 + CO (+M) <=> CH2CO (+M) */
    phi_f = sc[10]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 8.1e+11*exp(0.5*tc[0]-2269.75/tc[1]);
    redP = 1e-12 * alpha / k_f * 2.69e+33*exp(-5.11*tc[0]-3570.71/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.4093*exp(T/-275))+ (0.5907*exp(T/-1226))+ (exp(-5185/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[28];
    Kc = 1.0 / (refC) * exp((g_RT[10] + g_RT[14]) - (g_RT[28]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[139] = q_f - q_r;

    /*reaction 141: CH2 + HCCO <=> C2H3 + CO */
    phi_f = sc[10]*sc[27];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[14];
    Kc = exp((g_RT[10] + g_RT[27]) - (g_RT[23] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[140] = q_f - q_r;

    /*reaction 142: CH2(S) + N2 <=> CH2 + N2 */
    phi_f = sc[11]*sc[30];
    k_f = 1e-06 * 1.5e+13*exp(-301.963/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[30];
    Kc = exp((g_RT[11] + g_RT[30]) - (g_RT[10] + g_RT[30]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[141] = q_f - q_r;

    /*reaction 143: CH2(S) + AR <=> CH2 + AR */
    phi_f = sc[11]*sc[31];
    k_f = 1e-06 * 9e+12*exp(-301.963/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[31];
    Kc = exp((g_RT[11] + g_RT[31]) - (g_RT[10] + g_RT[31]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[142] = q_f - q_r;

    /*reaction 144: CH2(S) + O2 <=> H + OH + CO */
    phi_f = sc[11]*sc[3];
    k_f = 1e-06 * 2.8e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[4]*sc[14];
    Kc = refC * exp((g_RT[11] + g_RT[3]) - (g_RT[1] + g_RT[4] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[143] = q_f - q_r;

    /*reaction 145: CH2(S) + O2 <=> CO + H2O */
    phi_f = sc[11]*sc[3];
    k_f = 1e-06 * 1.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[5];
    Kc = exp((g_RT[11] + g_RT[3]) - (g_RT[14] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[144] = q_f - q_r;

    /*reaction 146: CH2(S) + H2 <=> CH3 + H */
    phi_f = sc[11]*sc[0];
    k_f = 1e-06 * 7e+13;
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[1];
    Kc = exp((g_RT[11] + g_RT[0]) - (g_RT[12] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[145] = q_f - q_r;

    /*reaction 147: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    phi_f = sc[11]*sc[5];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = 1e-06 * 2e+13;
    redP = 1e-12 * alpha / k_f * 2.7e+38*exp(-6.3*tc[0]-1560.14/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.8493*exp(T/-134))+ (0.1507*exp(T/-2383))+ (exp(-7265/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[20];
    Kc = 1.0 / (refC) * exp((g_RT[11] + g_RT[5]) - (g_RT[20]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[146] = q_f - q_r;

    /*reaction 148: CH2(S) + H2O <=> CH2 + H2O */
    phi_f = sc[11]*sc[5];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[5];
    Kc = exp((g_RT[11] + g_RT[5]) - (g_RT[10] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[147] = q_f - q_r;

    /*reaction 149: CH2(S) + CH3 <=> H + C2H4 */
    phi_f = sc[11]*sc[12];
    k_f = 1e-06 * 1.2e+13*exp(+286.865/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = exp((g_RT[11] + g_RT[12]) - (g_RT[1] + g_RT[24]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[148] = q_f - q_r;

    /*reaction 150: CH2(S) + CH4 <=> 2 CH3 */
    phi_f = sc[11]*sc[13];
    k_f = 1e-06 * 1.6e+13*exp(+286.865/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[12];
    Kc = exp((g_RT[11] + g_RT[13]) - (2 * g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[149] = q_f - q_r;

    /*reaction 151: CH2(S) + CO <=> CH2 + CO */
    phi_f = sc[11]*sc[14];
    k_f = 1e-06 * 9e+12;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[14];
    Kc = exp((g_RT[11] + g_RT[14]) - (g_RT[10] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[150] = q_f - q_r;

    /*reaction 152: CH2(S) + CO2 <=> CH2 + CO2 */
    phi_f = sc[11]*sc[15];
    k_f = 1e-06 * 7e+12;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[15];
    Kc = exp((g_RT[11] + g_RT[15]) - (g_RT[10] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[151] = q_f - q_r;

    /*reaction 153: CH2(S) + CO2 <=> CO + CH2O */
    phi_f = sc[11]*sc[15];
    k_f = 1e-06 * 1.4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[17];
    Kc = exp((g_RT[11] + g_RT[15]) - (g_RT[14] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[152] = q_f - q_r;

    /*reaction 154: CH2(S) + C2H6 <=> CH3 + C2H5 */
    phi_f = sc[11]*sc[26];
    k_f = 1e-06 * 4e+13*exp(+276.799/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[25];
    Kc = exp((g_RT[11] + g_RT[26]) - (g_RT[12] + g_RT[25]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[153] = q_f - q_r;

    /*reaction 155: CH3 + O2 <=> O + CH3O */
    phi_f = sc[12]*sc[3];
    k_f = 1e-06 * 2.675e+13*exp(-14494.2/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[19];
    Kc = exp((g_RT[12] + g_RT[3]) - (g_RT[2] + g_RT[19]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[154] = q_f - q_r;

    /*reaction 156: CH3 + O2 <=> OH + CH2O */
    phi_f = sc[12]*sc[3];
    k_f = 1e-06 * 3.6e+10*exp(-4499.25/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = exp((g_RT[12] + g_RT[3]) - (g_RT[4] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[155] = q_f - q_r;

    /*reaction 157: CH3 + H2O2 <=> HO2 + CH4 */
    phi_f = sc[12]*sc[7];
    k_f = 1e-06 * 24500*exp(2.47*tc[0]-2606.95/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = exp((g_RT[12] + g_RT[7]) - (g_RT[6] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[156] = q_f - q_r;

    /*reaction 158: 2 CH3 (+M) <=> C2H6 (+M) */
    phi_f = sc[12]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1e-06 * 2.12e+16*exp(-0.97*tc[0]-312.028/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.77e+50*exp(-9.67*tc[0]-3130.35/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.4675*exp(T/-151))+ (0.5325*exp(T/-1038))+ (exp(-4970/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[26];
    Kc = 1.0 / (refC) * exp((2 * g_RT[12]) - (g_RT[26]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[157] = q_f - q_r;

    /*reaction 159: 2 CH3 <=> H + C2H5 */
    phi_f = sc[12]*sc[12];
    k_f = 1e-06 * 4.99e+12*exp(0.1*tc[0]-5334.68/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[25];
    Kc = exp((2 * g_RT[12]) - (g_RT[1] + g_RT[25]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[158] = q_f - q_r;

    /*reaction 160: CH3 + HCO <=> CH4 + CO */
    phi_f = sc[12]*sc[16];
    k_f = 1e-06 * 2.648e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[14];
    Kc = exp((g_RT[12] + g_RT[16]) - (g_RT[13] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[159] = q_f - q_r;

    /*reaction 161: CH3 + CH2O <=> HCO + CH4 */
    phi_f = sc[12]*sc[17];
    k_f = 1e-06 * 3320*exp(2.81*tc[0]-2949.17/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[13];
    Kc = exp((g_RT[12] + g_RT[17]) - (g_RT[16] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[160] = q_f - q_r;

    /*reaction 162: CH3 + CH3OH <=> CH2OH + CH4 */
    phi_f = sc[12]*sc[20];
    k_f = 1e-06 * 3e+07*exp(1.5*tc[0]-5002.52/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[13];
    Kc = exp((g_RT[12] + g_RT[20]) - (g_RT[18] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[161] = q_f - q_r;

    /*reaction 163: CH3 + CH3OH <=> CH3O + CH4 */
    phi_f = sc[12]*sc[20];
    k_f = 1e-06 * 1e+07*exp(1.5*tc[0]-5002.52/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[13];
    Kc = exp((g_RT[12] + g_RT[20]) - (g_RT[19] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[162] = q_f - q_r;

    /*reaction 164: CH3 + C2H4 <=> C2H3 + CH4 */
    phi_f = sc[12]*sc[24];
    k_f = 1e-06 * 227000*exp(2*tc[0]-4630.1/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[13];
    Kc = exp((g_RT[12] + g_RT[24]) - (g_RT[23] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[163] = q_f - q_r;

    /*reaction 165: CH3 + C2H6 <=> C2H5 + CH4 */
    phi_f = sc[12]*sc[26];
    k_f = 1e-06 * 6.14e+06*exp(1.74*tc[0]-5259.18/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[13];
    Kc = exp((g_RT[12] + g_RT[26]) - (g_RT[25] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[164] = q_f - q_r;

    /*reaction 166: HCO + H2O <=> H + CO + H2O */
    phi_f = sc[16]*sc[5];
    k_f = 1e-06 * 2.244e+18*exp(-1*tc[0]-8555.61/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14]*sc[5];
    Kc = refC * exp((g_RT[16] + g_RT[5]) - (g_RT[1] + g_RT[14] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[165] = q_f - q_r;

    /*reaction 167: HCO + M <=> H + CO + M */
    phi_f = sc[16];
    alpha = mixture + sc[0] + -1*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = 1e-06 * alpha * 1.87e+17*exp(-1*tc[0]-8555.61/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = refC * exp((g_RT[16]) - (g_RT[1] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[166] = q_f - q_r;

    /*reaction 168: HCO + O2 <=> HO2 + CO */
    phi_f = sc[16]*sc[3];
    k_f = 1e-06 * 7.6e+12*exp(-201.309/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[14];
    Kc = exp((g_RT[16] + g_RT[3]) - (g_RT[6] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[167] = q_f - q_r;

    /*reaction 169: CH2OH + O2 <=> HO2 + CH2O */
    phi_f = sc[18]*sc[3];
    k_f = 1e-06 * 1.8e+13*exp(-452.944/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = exp((g_RT[18] + g_RT[3]) - (g_RT[6] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[168] = q_f - q_r;

    /*reaction 170: CH3O + O2 <=> HO2 + CH2O */
    phi_f = sc[19]*sc[3];
    k_f = 1e-06 * 4.28e-13*exp(7.6*tc[0]+1776.55/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = exp((g_RT[19] + g_RT[3]) - (g_RT[6] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[169] = q_f - q_r;

    /*reaction 171: C2H + O2 <=> HCO + CO */
    phi_f = sc[21]*sc[3];
    k_f = 1e-06 * 5e+13*exp(-754.907/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[14];
    Kc = exp((g_RT[21] + g_RT[3]) - (g_RT[16] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[170] = q_f - q_r;

    /*reaction 172: C2H + H2 <=> H + C2H2 */
    phi_f = sc[21]*sc[0];
    k_f = 1e-06 * 407000*exp(2.4*tc[0]-100.654/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = exp((g_RT[21] + g_RT[0]) - (g_RT[1] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[171] = q_f - q_r;

    /*reaction 173: C2H3 + O2 <=> HCO + CH2O */
    phi_f = sc[23]*sc[3];
    k_f = 1e-06 * 3.98e+12*exp(+120.785/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[17];
    Kc = exp((g_RT[23] + g_RT[3]) - (g_RT[16] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[172] = q_f - q_r;

    /*reaction 174: C2H4 (+M) <=> H2 + C2H2 (+M) */
    phi_f = sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = 1 * 8e+12*exp(0.44*tc[0]-44675.4/tc[1]);
    redP = 1e-12 * alpha / k_f * 7e+50*exp(-9.31*tc[0]-50256.7/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2655*exp(T/-180))+ (0.7345*exp(T/-1035))+ (exp(-5417/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = refC * exp((g_RT[24]) - (g_RT[0] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[173] = q_f - q_r;

    /*reaction 175: C2H5 + O2 <=> HO2 + C2H4 */
    phi_f = sc[25]*sc[3];
    k_f = 1e-06 * 8.4e+11*exp(-1950.18/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[24];
    Kc = exp((g_RT[25] + g_RT[3]) - (g_RT[6] + g_RT[24]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[174] = q_f - q_r;

    /*reaction 176: HCCO + O2 <=> OH + 2 CO */
    phi_f = sc[27]*sc[3];
    k_f = 1e-06 * 1.6e+12*exp(-429.794/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[14]*sc[14];
    Kc = refC * exp((g_RT[27] + g_RT[3]) - (g_RT[4] + 2 * g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[175] = q_f - q_r;

    /*reaction 177: 2 HCCO <=> 2 CO + C2H2 */
    phi_f = sc[27]*sc[27];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[14]*sc[22];
    Kc = refC * exp((2 * g_RT[27]) - (2 * g_RT[14] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[176] = q_f - q_r;

    return;
}


/*compute the equilibrium constants for each reaction */
void equilibriumConstants(double *kc, double * g_RT, double T)
{
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.314 / T;

    /*reaction 1: 2 O + M <=> O2 + M */
    kc[0] = 1.0 / (refC) * exp((2 * g_RT[2]) - (g_RT[3]));

    /*reaction 2: O + H + M <=> OH + M */
    kc[1] = 1.0 / (refC) * exp((g_RT[2] + g_RT[1]) - (g_RT[4]));

    /*reaction 3: O + H2 <=> H + OH */
    kc[2] = exp((g_RT[2] + g_RT[0]) - (g_RT[1] + g_RT[4]));

    /*reaction 4: O + HO2 <=> OH + O2 */
    kc[3] = exp((g_RT[2] + g_RT[6]) - (g_RT[4] + g_RT[3]));

    /*reaction 5: O + H2O2 <=> OH + HO2 */
    kc[4] = exp((g_RT[2] + g_RT[7]) - (g_RT[4] + g_RT[6]));

    /*reaction 6: O + CH <=> H + CO */
    kc[5] = exp((g_RT[2] + g_RT[9]) - (g_RT[1] + g_RT[14]));

    /*reaction 7: O + CH2 <=> H + HCO */
    kc[6] = exp((g_RT[2] + g_RT[10]) - (g_RT[1] + g_RT[16]));

    /*reaction 8: O + CH2(S) <=> H2 + CO */
    kc[7] = exp((g_RT[2] + g_RT[11]) - (g_RT[0] + g_RT[14]));

    /*reaction 9: O + CH2(S) <=> H + HCO */
    kc[8] = exp((g_RT[2] + g_RT[11]) - (g_RT[1] + g_RT[16]));

    /*reaction 10: O + CH3 <=> H + CH2O */
    kc[9] = exp((g_RT[2] + g_RT[12]) - (g_RT[1] + g_RT[17]));

    /*reaction 11: O + CH4 <=> OH + CH3 */
    kc[10] = exp((g_RT[2] + g_RT[13]) - (g_RT[4] + g_RT[12]));

    /*reaction 12: O + CO + M <=> CO2 + M */
    kc[11] = 1.0 / (refC) * exp((g_RT[2] + g_RT[14]) - (g_RT[15]));

    /*reaction 13: O + HCO <=> OH + CO */
    kc[12] = exp((g_RT[2] + g_RT[16]) - (g_RT[4] + g_RT[14]));

    /*reaction 14: O + HCO <=> H + CO2 */
    kc[13] = exp((g_RT[2] + g_RT[16]) - (g_RT[1] + g_RT[15]));

    /*reaction 15: O + CH2O <=> OH + HCO */
    kc[14] = exp((g_RT[2] + g_RT[17]) - (g_RT[4] + g_RT[16]));

    /*reaction 16: O + CH2OH <=> OH + CH2O */
    kc[15] = exp((g_RT[2] + g_RT[18]) - (g_RT[4] + g_RT[17]));

    /*reaction 17: O + CH3O <=> OH + CH2O */
    kc[16] = exp((g_RT[2] + g_RT[19]) - (g_RT[4] + g_RT[17]));

    /*reaction 18: O + CH3OH <=> OH + CH2OH */
    kc[17] = exp((g_RT[2] + g_RT[20]) - (g_RT[4] + g_RT[18]));

    /*reaction 19: O + CH3OH <=> OH + CH3O */
    kc[18] = exp((g_RT[2] + g_RT[20]) - (g_RT[4] + g_RT[19]));

    /*reaction 20: O + C2H <=> CH + CO */
    kc[19] = exp((g_RT[2] + g_RT[21]) - (g_RT[9] + g_RT[14]));

    /*reaction 21: O + C2H2 <=> H + HCCO */
    kc[20] = exp((g_RT[2] + g_RT[22]) - (g_RT[1] + g_RT[27]));

    /*reaction 22: O + C2H2 <=> OH + C2H */
    kc[21] = exp((g_RT[2] + g_RT[22]) - (g_RT[4] + g_RT[21]));

    /*reaction 23: O + C2H2 <=> CO + CH2 */
    kc[22] = exp((g_RT[2] + g_RT[22]) - (g_RT[14] + g_RT[10]));

    /*reaction 24: O + C2H3 <=> H + CH2CO */
    kc[23] = exp((g_RT[2] + g_RT[23]) - (g_RT[1] + g_RT[28]));

    /*reaction 25: O + C2H4 <=> CH3 + HCO */
    kc[24] = exp((g_RT[2] + g_RT[24]) - (g_RT[12] + g_RT[16]));

    /*reaction 26: O + C2H5 <=> CH3 + CH2O */
    kc[25] = exp((g_RT[2] + g_RT[25]) - (g_RT[12] + g_RT[17]));

    /*reaction 27: O + C2H6 <=> OH + C2H5 */
    kc[26] = exp((g_RT[2] + g_RT[26]) - (g_RT[4] + g_RT[25]));

    /*reaction 28: O + HCCO <=> H + 2 CO */
    kc[27] = refC * exp((g_RT[2] + g_RT[27]) - (g_RT[1] + 2 * g_RT[14]));

    /*reaction 29: O + CH2CO <=> OH + HCCO */
    kc[28] = exp((g_RT[2] + g_RT[28]) - (g_RT[4] + g_RT[27]));

    /*reaction 30: O + CH2CO <=> CH2 + CO2 */
    kc[29] = exp((g_RT[2] + g_RT[28]) - (g_RT[10] + g_RT[15]));

    /*reaction 31: O2 + CO <=> O + CO2 */
    kc[30] = exp((g_RT[3] + g_RT[14]) - (g_RT[2] + g_RT[15]));

    /*reaction 32: O2 + CH2O <=> HO2 + HCO */
    kc[31] = exp((g_RT[3] + g_RT[17]) - (g_RT[6] + g_RT[16]));

    /*reaction 33: H + O2 + M <=> HO2 + M */
    kc[32] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));

    /*reaction 34: H + 2 O2 <=> HO2 + O2 */
    kc[33] = 1.0 / (refC) * exp((g_RT[1] + 2 * g_RT[3]) - (g_RT[6] + g_RT[3]));

    /*reaction 35: H + O2 + H2O <=> HO2 + H2O */
    kc[34] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));

    /*reaction 36: H + O2 + N2 <=> HO2 + N2 */
    kc[35] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[30]) - (g_RT[6] + g_RT[30]));

    /*reaction 37: H + O2 + AR <=> HO2 + AR */
    kc[36] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[31]) - (g_RT[6] + g_RT[31]));

    /*reaction 38: H + O2 <=> O + OH */
    kc[37] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));

    /*reaction 39: 2 H + M <=> H2 + M */
    kc[38] = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    kc[39] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    kc[40] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    kc[41] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[15]) - (g_RT[0] + g_RT[15]));

    /*reaction 43: H + OH + M <=> H2O + M */
    kc[42] = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));

    /*reaction 44: H + HO2 <=> O + H2O */
    kc[43] = exp((g_RT[1] + g_RT[6]) - (g_RT[2] + g_RT[5]));

    /*reaction 45: H + HO2 <=> O2 + H2 */
    kc[44] = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));

    /*reaction 46: H + HO2 <=> 2 OH */
    kc[45] = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));

    /*reaction 47: H + H2O2 <=> HO2 + H2 */
    kc[46] = exp((g_RT[1] + g_RT[7]) - (g_RT[6] + g_RT[0]));

    /*reaction 48: H + H2O2 <=> OH + H2O */
    kc[47] = exp((g_RT[1] + g_RT[7]) - (g_RT[4] + g_RT[5]));

    /*reaction 49: H + CH <=> C + H2 */
    kc[48] = exp((g_RT[1] + g_RT[9]) - (g_RT[8] + g_RT[0]));

    /*reaction 50: H + CH2 (+M) <=> CH3 (+M) */
    kc[49] = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[12]));

    /*reaction 51: H + CH2(S) <=> CH + H2 */
    kc[50] = exp((g_RT[1] + g_RT[11]) - (g_RT[9] + g_RT[0]));

    /*reaction 52: H + CH3 (+M) <=> CH4 (+M) */
    kc[51] = 1.0 / (refC) * exp((g_RT[1] + g_RT[12]) - (g_RT[13]));

    /*reaction 53: H + CH4 <=> CH3 + H2 */
    kc[52] = exp((g_RT[1] + g_RT[13]) - (g_RT[12] + g_RT[0]));

    /*reaction 54: H + HCO (+M) <=> CH2O (+M) */
    kc[53] = 1.0 / (refC) * exp((g_RT[1] + g_RT[16]) - (g_RT[17]));

    /*reaction 55: H + HCO <=> H2 + CO */
    kc[54] = exp((g_RT[1] + g_RT[16]) - (g_RT[0] + g_RT[14]));

    /*reaction 56: H + CH2O (+M) <=> CH2OH (+M) */
    kc[55] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[18]));

    /*reaction 57: H + CH2O (+M) <=> CH3O (+M) */
    kc[56] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[19]));

    /*reaction 58: H + CH2O <=> HCO + H2 */
    kc[57] = exp((g_RT[1] + g_RT[17]) - (g_RT[16] + g_RT[0]));

    /*reaction 59: H + CH2OH (+M) <=> CH3OH (+M) */
    kc[58] = 1.0 / (refC) * exp((g_RT[1] + g_RT[18]) - (g_RT[20]));

    /*reaction 60: H + CH2OH <=> H2 + CH2O */
    kc[59] = exp((g_RT[1] + g_RT[18]) - (g_RT[0] + g_RT[17]));

    /*reaction 61: H + CH2OH <=> OH + CH3 */
    kc[60] = exp((g_RT[1] + g_RT[18]) - (g_RT[4] + g_RT[12]));

    /*reaction 62: H + CH2OH <=> CH2(S) + H2O */
    kc[61] = exp((g_RT[1] + g_RT[18]) - (g_RT[11] + g_RT[5]));

    /*reaction 63: H + CH3O (+M) <=> CH3OH (+M) */
    kc[62] = 1.0 / (refC) * exp((g_RT[1] + g_RT[19]) - (g_RT[20]));

    /*reaction 64: H + CH3O <=> H + CH2OH */
    kc[63] = exp((g_RT[1] + g_RT[19]) - (g_RT[1] + g_RT[18]));

    /*reaction 65: H + CH3O <=> H2 + CH2O */
    kc[64] = exp((g_RT[1] + g_RT[19]) - (g_RT[0] + g_RT[17]));

    /*reaction 66: H + CH3O <=> OH + CH3 */
    kc[65] = exp((g_RT[1] + g_RT[19]) - (g_RT[4] + g_RT[12]));

    /*reaction 67: H + CH3O <=> CH2(S) + H2O */
    kc[66] = exp((g_RT[1] + g_RT[19]) - (g_RT[11] + g_RT[5]));

    /*reaction 68: H + CH3OH <=> CH2OH + H2 */
    kc[67] = exp((g_RT[1] + g_RT[20]) - (g_RT[18] + g_RT[0]));

    /*reaction 69: H + CH3OH <=> CH3O + H2 */
    kc[68] = exp((g_RT[1] + g_RT[20]) - (g_RT[19] + g_RT[0]));

    /*reaction 70: H + C2H (+M) <=> C2H2 (+M) */
    kc[69] = 1.0 / (refC) * exp((g_RT[1] + g_RT[21]) - (g_RT[22]));

    /*reaction 71: H + C2H2 (+M) <=> C2H3 (+M) */
    kc[70] = 1.0 / (refC) * exp((g_RT[1] + g_RT[22]) - (g_RT[23]));

    /*reaction 72: H + C2H3 (+M) <=> C2H4 (+M) */
    kc[71] = 1.0 / (refC) * exp((g_RT[1] + g_RT[23]) - (g_RT[24]));

    /*reaction 73: H + C2H3 <=> H2 + C2H2 */
    kc[72] = exp((g_RT[1] + g_RT[23]) - (g_RT[0] + g_RT[22]));

    /*reaction 74: H + C2H4 (+M) <=> C2H5 (+M) */
    kc[73] = 1.0 / (refC) * exp((g_RT[1] + g_RT[24]) - (g_RT[25]));

    /*reaction 75: H + C2H4 <=> C2H3 + H2 */
    kc[74] = exp((g_RT[1] + g_RT[24]) - (g_RT[23] + g_RT[0]));

    /*reaction 76: H + C2H5 (+M) <=> C2H6 (+M) */
    kc[75] = 1.0 / (refC) * exp((g_RT[1] + g_RT[25]) - (g_RT[26]));

    /*reaction 77: H + C2H5 <=> H2 + C2H4 */
    kc[76] = exp((g_RT[1] + g_RT[25]) - (g_RT[0] + g_RT[24]));

    /*reaction 78: H + C2H6 <=> C2H5 + H2 */
    kc[77] = exp((g_RT[1] + g_RT[26]) - (g_RT[25] + g_RT[0]));

    /*reaction 79: H + HCCO <=> CH2(S) + CO */
    kc[78] = exp((g_RT[1] + g_RT[27]) - (g_RT[11] + g_RT[14]));

    /*reaction 80: H + CH2CO <=> HCCO + H2 */
    kc[79] = exp((g_RT[1] + g_RT[28]) - (g_RT[27] + g_RT[0]));

    /*reaction 81: H + CH2CO <=> CH3 + CO */
    kc[80] = exp((g_RT[1] + g_RT[28]) - (g_RT[12] + g_RT[14]));

    /*reaction 82: H + HCCOH <=> H + CH2CO */
    kc[81] = exp((g_RT[1] + g_RT[29]) - (g_RT[1] + g_RT[28]));

    /*reaction 83: H2 + CO (+M) <=> CH2O (+M) */
    kc[82] = 1.0 / (refC) * exp((g_RT[0] + g_RT[14]) - (g_RT[17]));

    /*reaction 84: OH + H2 <=> H + H2O */
    kc[83] = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));

    /*reaction 85: 2 OH (+M) <=> H2O2 (+M) */
    kc[84] = 1.0 / (refC) * exp((2 * g_RT[4]) - (g_RT[7]));

    /*reaction 86: 2 OH <=> O + H2O */
    kc[85] = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));

    /*reaction 87: OH + HO2 <=> O2 + H2O */
    kc[86] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    kc[87] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));

    /*reaction 89: OH + H2O2 <=> HO2 + H2O */
    kc[88] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));

    /*reaction 90: OH + C <=> H + CO */
    kc[89] = exp((g_RT[4] + g_RT[8]) - (g_RT[1] + g_RT[14]));

    /*reaction 91: OH + CH <=> H + HCO */
    kc[90] = exp((g_RT[4] + g_RT[9]) - (g_RT[1] + g_RT[16]));

    /*reaction 92: OH + CH2 <=> H + CH2O */
    kc[91] = exp((g_RT[4] + g_RT[10]) - (g_RT[1] + g_RT[17]));

    /*reaction 93: OH + CH2 <=> CH + H2O */
    kc[92] = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[5]));

    /*reaction 94: OH + CH2(S) <=> H + CH2O */
    kc[93] = exp((g_RT[4] + g_RT[11]) - (g_RT[1] + g_RT[17]));

    /*reaction 95: OH + CH3 (+M) <=> CH3OH (+M) */
    kc[94] = 1.0 / (refC) * exp((g_RT[4] + g_RT[12]) - (g_RT[20]));

    /*reaction 96: OH + CH3 <=> CH2 + H2O */
    kc[95] = exp((g_RT[4] + g_RT[12]) - (g_RT[10] + g_RT[5]));

    /*reaction 97: OH + CH3 <=> CH2(S) + H2O */
    kc[96] = exp((g_RT[4] + g_RT[12]) - (g_RT[11] + g_RT[5]));

    /*reaction 98: OH + CH4 <=> CH3 + H2O */
    kc[97] = exp((g_RT[4] + g_RT[13]) - (g_RT[12] + g_RT[5]));

    /*reaction 99: OH + CO <=> H + CO2 */
    kc[98] = exp((g_RT[4] + g_RT[14]) - (g_RT[1] + g_RT[15]));

    /*reaction 100: OH + HCO <=> H2O + CO */
    kc[99] = exp((g_RT[4] + g_RT[16]) - (g_RT[5] + g_RT[14]));

    /*reaction 101: OH + CH2O <=> HCO + H2O */
    kc[100] = exp((g_RT[4] + g_RT[17]) - (g_RT[16] + g_RT[5]));

    /*reaction 102: OH + CH2OH <=> H2O + CH2O */
    kc[101] = exp((g_RT[4] + g_RT[18]) - (g_RT[5] + g_RT[17]));

    /*reaction 103: OH + CH3O <=> H2O + CH2O */
    kc[102] = exp((g_RT[4] + g_RT[19]) - (g_RT[5] + g_RT[17]));

    /*reaction 104: OH + CH3OH <=> CH2OH + H2O */
    kc[103] = exp((g_RT[4] + g_RT[20]) - (g_RT[18] + g_RT[5]));

    /*reaction 105: OH + CH3OH <=> CH3O + H2O */
    kc[104] = exp((g_RT[4] + g_RT[20]) - (g_RT[19] + g_RT[5]));

    /*reaction 106: OH + C2H <=> H + HCCO */
    kc[105] = exp((g_RT[4] + g_RT[21]) - (g_RT[1] + g_RT[27]));

    /*reaction 107: OH + C2H2 <=> H + CH2CO */
    kc[106] = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[28]));

    /*reaction 108: OH + C2H2 <=> H + HCCOH */
    kc[107] = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[29]));

    /*reaction 109: OH + C2H2 <=> C2H + H2O */
    kc[108] = exp((g_RT[4] + g_RT[22]) - (g_RT[21] + g_RT[5]));

    /*reaction 110: OH + C2H2 <=> CH3 + CO */
    kc[109] = exp((g_RT[4] + g_RT[22]) - (g_RT[12] + g_RT[14]));

    /*reaction 111: OH + C2H3 <=> H2O + C2H2 */
    kc[110] = exp((g_RT[4] + g_RT[23]) - (g_RT[5] + g_RT[22]));

    /*reaction 112: OH + C2H4 <=> C2H3 + H2O */
    kc[111] = exp((g_RT[4] + g_RT[24]) - (g_RT[23] + g_RT[5]));

    /*reaction 113: OH + C2H6 <=> C2H5 + H2O */
    kc[112] = exp((g_RT[4] + g_RT[26]) - (g_RT[25] + g_RT[5]));

    /*reaction 114: OH + CH2CO <=> HCCO + H2O */
    kc[113] = exp((g_RT[4] + g_RT[28]) - (g_RT[27] + g_RT[5]));

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    kc[114] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));

    /*reaction 116: 2 HO2 <=> O2 + H2O2 */
    kc[115] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));

    /*reaction 117: HO2 + CH2 <=> OH + CH2O */
    kc[116] = exp((g_RT[6] + g_RT[10]) - (g_RT[4] + g_RT[17]));

    /*reaction 118: HO2 + CH3 <=> O2 + CH4 */
    kc[117] = exp((g_RT[6] + g_RT[12]) - (g_RT[3] + g_RT[13]));

    /*reaction 119: HO2 + CH3 <=> OH + CH3O */
    kc[118] = exp((g_RT[6] + g_RT[12]) - (g_RT[4] + g_RT[19]));

    /*reaction 120: HO2 + CO <=> OH + CO2 */
    kc[119] = exp((g_RT[6] + g_RT[14]) - (g_RT[4] + g_RT[15]));

    /*reaction 121: HO2 + CH2O <=> HCO + H2O2 */
    kc[120] = exp((g_RT[6] + g_RT[17]) - (g_RT[16] + g_RT[7]));

    /*reaction 122: C + O2 <=> O + CO */
    kc[121] = exp((g_RT[8] + g_RT[3]) - (g_RT[2] + g_RT[14]));

    /*reaction 123: C + CH2 <=> H + C2H */
    kc[122] = exp((g_RT[8] + g_RT[10]) - (g_RT[1] + g_RT[21]));

    /*reaction 124: C + CH3 <=> H + C2H2 */
    kc[123] = exp((g_RT[8] + g_RT[12]) - (g_RT[1] + g_RT[22]));

    /*reaction 125: CH + O2 <=> O + HCO */
    kc[124] = exp((g_RT[9] + g_RT[3]) - (g_RT[2] + g_RT[16]));

    /*reaction 126: CH + H2 <=> H + CH2 */
    kc[125] = exp((g_RT[9] + g_RT[0]) - (g_RT[1] + g_RT[10]));

    /*reaction 127: CH + H2O <=> H + CH2O */
    kc[126] = exp((g_RT[9] + g_RT[5]) - (g_RT[1] + g_RT[17]));

    /*reaction 128: CH + CH2 <=> H + C2H2 */
    kc[127] = exp((g_RT[9] + g_RT[10]) - (g_RT[1] + g_RT[22]));

    /*reaction 129: CH + CH3 <=> H + C2H3 */
    kc[128] = exp((g_RT[9] + g_RT[12]) - (g_RT[1] + g_RT[23]));

    /*reaction 130: CH + CH4 <=> H + C2H4 */
    kc[129] = exp((g_RT[9] + g_RT[13]) - (g_RT[1] + g_RT[24]));

    /*reaction 131: CH + CO (+M) <=> HCCO (+M) */
    kc[130] = 1.0 / (refC) * exp((g_RT[9] + g_RT[14]) - (g_RT[27]));

    /*reaction 132: CH + CO2 <=> HCO + CO */
    kc[131] = exp((g_RT[9] + g_RT[15]) - (g_RT[16] + g_RT[14]));

    /*reaction 133: CH + CH2O <=> H + CH2CO */
    kc[132] = exp((g_RT[9] + g_RT[17]) - (g_RT[1] + g_RT[28]));

    /*reaction 134: CH + HCCO <=> CO + C2H2 */
    kc[133] = exp((g_RT[9] + g_RT[27]) - (g_RT[14] + g_RT[22]));

    /*reaction 135: CH2 + O2 <=> OH + HCO */
    kc[134] = exp((g_RT[10] + g_RT[3]) - (g_RT[4] + g_RT[16]));

    /*reaction 136: CH2 + H2 <=> H + CH3 */
    kc[135] = exp((g_RT[10] + g_RT[0]) - (g_RT[1] + g_RT[12]));

    /*reaction 137: 2 CH2 <=> H2 + C2H2 */
    kc[136] = exp((2 * g_RT[10]) - (g_RT[0] + g_RT[22]));

    /*reaction 138: CH2 + CH3 <=> H + C2H4 */
    kc[137] = exp((g_RT[10] + g_RT[12]) - (g_RT[1] + g_RT[24]));

    /*reaction 139: CH2 + CH4 <=> 2 CH3 */
    kc[138] = exp((g_RT[10] + g_RT[13]) - (2 * g_RT[12]));

    /*reaction 140: CH2 + CO (+M) <=> CH2CO (+M) */
    kc[139] = 1.0 / (refC) * exp((g_RT[10] + g_RT[14]) - (g_RT[28]));

    /*reaction 141: CH2 + HCCO <=> C2H3 + CO */
    kc[140] = exp((g_RT[10] + g_RT[27]) - (g_RT[23] + g_RT[14]));

    /*reaction 142: CH2(S) + N2 <=> CH2 + N2 */
    kc[141] = exp((g_RT[11] + g_RT[30]) - (g_RT[10] + g_RT[30]));

    /*reaction 143: CH2(S) + AR <=> CH2 + AR */
    kc[142] = exp((g_RT[11] + g_RT[31]) - (g_RT[10] + g_RT[31]));

    /*reaction 144: CH2(S) + O2 <=> H + OH + CO */
    kc[143] = refC * exp((g_RT[11] + g_RT[3]) - (g_RT[1] + g_RT[4] + g_RT[14]));

    /*reaction 145: CH2(S) + O2 <=> CO + H2O */
    kc[144] = exp((g_RT[11] + g_RT[3]) - (g_RT[14] + g_RT[5]));

    /*reaction 146: CH2(S) + H2 <=> CH3 + H */
    kc[145] = exp((g_RT[11] + g_RT[0]) - (g_RT[12] + g_RT[1]));

    /*reaction 147: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    kc[146] = 1.0 / (refC) * exp((g_RT[11] + g_RT[5]) - (g_RT[20]));

    /*reaction 148: CH2(S) + H2O <=> CH2 + H2O */
    kc[147] = exp((g_RT[11] + g_RT[5]) - (g_RT[10] + g_RT[5]));

    /*reaction 149: CH2(S) + CH3 <=> H + C2H4 */
    kc[148] = exp((g_RT[11] + g_RT[12]) - (g_RT[1] + g_RT[24]));

    /*reaction 150: CH2(S) + CH4 <=> 2 CH3 */
    kc[149] = exp((g_RT[11] + g_RT[13]) - (2 * g_RT[12]));

    /*reaction 151: CH2(S) + CO <=> CH2 + CO */
    kc[150] = exp((g_RT[11] + g_RT[14]) - (g_RT[10] + g_RT[14]));

    /*reaction 152: CH2(S) + CO2 <=> CH2 + CO2 */
    kc[151] = exp((g_RT[11] + g_RT[15]) - (g_RT[10] + g_RT[15]));

    /*reaction 153: CH2(S) + CO2 <=> CO + CH2O */
    kc[152] = exp((g_RT[11] + g_RT[15]) - (g_RT[14] + g_RT[17]));

    /*reaction 154: CH2(S) + C2H6 <=> CH3 + C2H5 */
    kc[153] = exp((g_RT[11] + g_RT[26]) - (g_RT[12] + g_RT[25]));

    /*reaction 155: CH3 + O2 <=> O + CH3O */
    kc[154] = exp((g_RT[12] + g_RT[3]) - (g_RT[2] + g_RT[19]));

    /*reaction 156: CH3 + O2 <=> OH + CH2O */
    kc[155] = exp((g_RT[12] + g_RT[3]) - (g_RT[4] + g_RT[17]));

    /*reaction 157: CH3 + H2O2 <=> HO2 + CH4 */
    kc[156] = exp((g_RT[12] + g_RT[7]) - (g_RT[6] + g_RT[13]));

    /*reaction 158: 2 CH3 (+M) <=> C2H6 (+M) */
    kc[157] = 1.0 / (refC) * exp((2 * g_RT[12]) - (g_RT[26]));

    /*reaction 159: 2 CH3 <=> H + C2H5 */
    kc[158] = exp((2 * g_RT[12]) - (g_RT[1] + g_RT[25]));

    /*reaction 160: CH3 + HCO <=> CH4 + CO */
    kc[159] = exp((g_RT[12] + g_RT[16]) - (g_RT[13] + g_RT[14]));

    /*reaction 161: CH3 + CH2O <=> HCO + CH4 */
    kc[160] = exp((g_RT[12] + g_RT[17]) - (g_RT[16] + g_RT[13]));

    /*reaction 162: CH3 + CH3OH <=> CH2OH + CH4 */
    kc[161] = exp((g_RT[12] + g_RT[20]) - (g_RT[18] + g_RT[13]));

    /*reaction 163: CH3 + CH3OH <=> CH3O + CH4 */
    kc[162] = exp((g_RT[12] + g_RT[20]) - (g_RT[19] + g_RT[13]));

    /*reaction 164: CH3 + C2H4 <=> C2H3 + CH4 */
    kc[163] = exp((g_RT[12] + g_RT[24]) - (g_RT[23] + g_RT[13]));

    /*reaction 165: CH3 + C2H6 <=> C2H5 + CH4 */
    kc[164] = exp((g_RT[12] + g_RT[26]) - (g_RT[25] + g_RT[13]));

    /*reaction 166: HCO + H2O <=> H + CO + H2O */
    kc[165] = refC * exp((g_RT[16] + g_RT[5]) - (g_RT[1] + g_RT[14] + g_RT[5]));

    /*reaction 167: HCO + M <=> H + CO + M */
    kc[166] = refC * exp((g_RT[16]) - (g_RT[1] + g_RT[14]));

    /*reaction 168: HCO + O2 <=> HO2 + CO */
    kc[167] = exp((g_RT[16] + g_RT[3]) - (g_RT[6] + g_RT[14]));

    /*reaction 169: CH2OH + O2 <=> HO2 + CH2O */
    kc[168] = exp((g_RT[18] + g_RT[3]) - (g_RT[6] + g_RT[17]));

    /*reaction 170: CH3O + O2 <=> HO2 + CH2O */
    kc[169] = exp((g_RT[19] + g_RT[3]) - (g_RT[6] + g_RT[17]));

    /*reaction 171: C2H + O2 <=> HCO + CO */
    kc[170] = exp((g_RT[21] + g_RT[3]) - (g_RT[16] + g_RT[14]));

    /*reaction 172: C2H + H2 <=> H + C2H2 */
    kc[171] = exp((g_RT[21] + g_RT[0]) - (g_RT[1] + g_RT[22]));

    /*reaction 173: C2H3 + O2 <=> HCO + CH2O */
    kc[172] = exp((g_RT[23] + g_RT[3]) - (g_RT[16] + g_RT[17]));

    /*reaction 174: C2H4 (+M) <=> H2 + C2H2 (+M) */
    kc[173] = refC * exp((g_RT[24]) - (g_RT[0] + g_RT[22]));

    /*reaction 175: C2H5 + O2 <=> HO2 + C2H4 */
    kc[174] = exp((g_RT[25] + g_RT[3]) - (g_RT[6] + g_RT[24]));

    /*reaction 176: HCCO + O2 <=> OH + 2 CO */
    kc[175] = refC * exp((g_RT[27] + g_RT[3]) - (g_RT[4] + 2 * g_RT[14]));

    /*reaction 177: 2 HCCO <=> 2 CO + C2H2 */
    kc[176] = refC * exp((2 * g_RT[27]) - (2 * g_RT[14] + g_RT[22]));

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void gibbs(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            -9.17935173e+02 / tc[1]
            +1.66132088e+00
            -2.34433112e+00 * tc[0]
            -3.99026037e-03 * tc[1]
            +3.24635850e-06 * tc[2]
            -1.67976745e-09 * tc[3]
            +3.68805881e-13 * tc[4];
        /*species 1: H */
        species[1] =
            +2.54736599e+04 / tc[1]
            +2.94668285e+00
            -2.50000000e+00 * tc[0]
            -3.52666409e-13 * tc[1]
            +3.32653273e-16 * tc[2]
            -1.91734693e-19 * tc[3]
            +4.63866166e-23 * tc[4];
        /*species 2: O */
        species[2] =
            +2.91222592e+04 / tc[1]
            +1.11633364e+00
            -3.16826710e+00 * tc[0]
            +1.63965942e-03 * tc[1]
            -1.10717733e-06 * tc[2]
            +5.10672187e-10 * tc[3]
            -1.05632985e-13 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.06394356e+03 / tc[1]
            +1.24780630e-01
            -3.78245636e+00 * tc[0]
            +1.49836708e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774591e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.61508056e+03 / tc[1]
            +4.09594089e+00
            -3.99201543e+00 * tc[0]
            +1.20065876e-03 * tc[1]
            -7.69656402e-07 * tc[2]
            +3.23427778e-10 * tc[3]
            -6.82057350e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.02937267e+04 / tc[1]
            +5.04767277e+00
            -4.19864056e+00 * tc[0]
            +1.01821705e-03 * tc[1]
            -1.08673369e-06 * tc[2]
            +4.57330885e-10 * tc[3]
            -8.85989085e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +2.94808040e+02 / tc[1]
            +5.85135560e-01
            -4.30179801e+00 * tc[0]
            +2.37456025e-03 * tc[1]
            -3.52638152e-06 * tc[2]
            +2.02303245e-09 * tc[3]
            -4.64612562e-13 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.77025821e+04 / tc[1]
            +8.41061950e-01
            -4.27611269e+00 * tc[0]
            +2.71411208e-04 * tc[1]
            -2.78892835e-06 * tc[2]
            +1.79809011e-09 * tc[3]
            -4.31227182e-13 * tc[4];
        /*species 8: C */
        species[8] =
            +8.54438832e+04 / tc[1]
            -1.97706893e+00
            -2.55423955e+00 * tc[0]
            +1.60768862e-04 * tc[1]
            -1.22298707e-07 * tc[2]
            +6.10195741e-11 * tc[3]
            -1.33260723e-14 * tc[4];
        /*species 9: CH */
        species[9] =
            +7.07972934e+04 / tc[1]
            +1.40580557e+00
            -3.48981665e+00 * tc[0]
            -1.61917771e-04 * tc[1]
            +2.81498442e-07 * tc[2]
            -2.63514439e-10 * tc[3]
            +7.03045335e-14 * tc[4];
        /*species 10: CH2 */
        species[10] =
            +4.60040401e+04 / tc[1]
            +2.20014682e+00
            -3.76267867e+00 * tc[0]
            -4.84436072e-04 * tc[1]
            -4.65816402e-07 * tc[2]
            +3.20909294e-10 * tc[3]
            -8.43708595e-14 * tc[4];
        /*species 11: CH2(S) */
        species[11] =
            +5.04968163e+04 / tc[1]
            +4.96772308e+00
            -4.19860411e+00 * tc[0]
            +1.18330710e-03 * tc[1]
            -1.37216037e-06 * tc[2]
            +5.57346651e-10 * tc[3]
            -9.71573685e-14 * tc[4];
        /*species 12: CH3 */
        species[12] =
            +1.64449988e+04 / tc[1]
            +2.06902607e+00
            -3.67359040e+00 * tc[0]
            -1.00547588e-03 * tc[1]
            -9.55036427e-07 * tc[2]
            +5.72597854e-10 * tc[3]
            -1.27192867e-13 * tc[4];
        /*species 13: CH4 */
        species[13] =
            -1.02466476e+04 / tc[1]
            +9.79117989e+00
            -5.14987613e+00 * tc[0]
            +6.83548940e-03 * tc[1]
            -8.19667665e-06 * tc[2]
            +4.03952522e-09 * tc[3]
            -8.33469780e-13 * tc[4];
        /*species 14: CO */
        species[14] =
            -1.43440860e+04 / tc[1]
            +7.11241900e-02
            -3.57953347e+00 * tc[0]
            +3.05176840e-04 * tc[1]
            -1.69469055e-07 * tc[2]
            -7.55838237e-11 * tc[3]
            +4.52212249e-14 * tc[4];
        /*species 15: CO2 */
        species[15] =
            -4.83719697e+04 / tc[1]
            -7.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
        /*species 16: HCO */
        species[16] =
            +3.83956496e+03 / tc[1]
            +8.26813410e-01
            -4.22118584e+00 * tc[0]
            +1.62196266e-03 * tc[1]
            -2.29665743e-06 * tc[2]
            +1.10953411e-09 * tc[3]
            -2.16884432e-13 * tc[4];
        /*species 17: CH2O */
        species[17] =
            -1.43089567e+04 / tc[1]
            +4.19091025e+00
            -4.79372315e+00 * tc[0]
            +4.95416684e-03 * tc[1]
            -6.22033347e-06 * tc[2]
            +3.16071051e-09 * tc[3]
            -6.58863260e-13 * tc[4];
        /*species 18: CH2OH */
        species[18] =
            -3.19391367e+03 / tc[1]
            -1.60913325e+00
            -3.86388918e+00 * tc[0]
            -2.79836152e-03 * tc[1]
            -9.88786318e-07 * tc[2]
            +8.71100100e-10 * tc[3]
            -2.18483639e-13 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +9.78601100e+02 / tc[1]
            -1.10459730e+01
            -2.10620400e+00 * tc[0]
            -3.60829750e-03 * tc[1]
            -8.89745333e-07 * tc[2]
            +6.14803000e-10 * tc[3]
            -1.03780500e-13 * tc[4];
        /*species 20: CH3OH */
        species[20] =
            -2.56427656e+04 / tc[1]
            +7.21949405e+00
            -5.71539582e+00 * tc[0]
            +7.61545645e-03 * tc[1]
            -1.08740193e-05 * tc[2]
            +5.92339074e-09 * tc[3]
            -1.30676349e-12 * tc[4];
        /*species 21: C2H */
        species[21] =
            +6.68393932e+04 / tc[1]
            -3.33330705e+00
            -2.88965733e+00 * tc[0]
            -6.70498055e-03 * tc[1]
            +4.74615835e-06 * tc[2]
            -2.45659204e-09 * tc[3]
            +5.46657555e-13 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            +2.64289807e+04 / tc[1]
            -1.31310240e+01
            -8.08681094e-01 * tc[0]
            -1.16807815e-02 * tc[1]
            +5.91953025e-06 * tc[2]
            -2.33460364e-09 * tc[3]
            +4.25036487e-13 * tc[4];
        /*species 23: C2H3 */
        species[23] =
            +3.48598468e+04 / tc[1]
            -5.29807380e+00
            -3.21246645e+00 * tc[0]
            -7.57395810e-04 * tc[1]
            -4.32015687e-06 * tc[2]
            +2.98048206e-09 * tc[3]
            -7.35754365e-13 * tc[4];
        /*species 24: C2H4 */
        species[24] =
            +5.08977593e+03 / tc[1]
            -1.38129480e-01
            -3.95920148e+00 * tc[0]
            +3.78526124e-03 * tc[1]
            -9.51650487e-06 * tc[2]
            +5.76323961e-09 * tc[3]
            -1.34942187e-12 * tc[4];
        /*species 25: C2H5 */
        species[25] =
            +1.28416265e+04 / tc[1]
            -4.00743560e-01
            -4.30646568e+00 * tc[0]
            +2.09329446e-03 * tc[1]
            -8.28571345e-06 * tc[2]
            +4.99272172e-09 * tc[3]
            -1.15254502e-12 * tc[4];
        /*species 26: C2H6 */
        species[26] =
            -1.15222055e+04 / tc[1]
            +1.62460176e+00
            -4.29142492e+00 * tc[0]
            +2.75077135e-03 * tc[1]
            -9.99063813e-06 * tc[2]
            +5.90388571e-09 * tc[3]
            -1.34342886e-12 * tc[4];
        /*species 27: HCCO */
        species[27] =
            +2.00594490e+04 / tc[1]
            -1.02386956e+01
            -2.25172140e+00 * tc[0]
            -8.82751050e-03 * tc[1]
            +3.95485017e-06 * tc[2]
            -1.43964658e-09 * tc[3]
            +2.53324055e-13 * tc[4];
        /*species 28: CH2CO */
        species[28] =
            -7.04291804e+03 / tc[1]
            -1.00798117e+01
            -2.13583630e+00 * tc[0]
            -9.05943605e-03 * tc[1]
            +2.89912457e-06 * tc[2]
            -7.78664640e-10 * tc[3]
            +1.00728807e-13 * tc[4];
        /*species 29: HCCOH */
        species[29] =
            +8.03161430e+03 / tc[1]
            -1.26319457e+01
            -1.24237330e+00 * tc[0]
            -1.55361005e-02 * tc[1]
            +8.47781067e-06 * tc[2]
            -3.59476092e-09 * tc[3]
            +7.00729700e-13 * tc[4];
        /*species 30: N2 */
        species[30] =
            -1.02089990e+03 / tc[1]
            -6.51695000e-01
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
        /*species 31: AR */
        species[31] =
            -7.45375000e+02 / tc[1]
            -1.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -9.50158922e+02 / tc[1]
            +6.54230251e+00
            -3.33727920e+00 * tc[0]
            +2.47012365e-05 * tc[1]
            -8.32427963e-08 * tc[2]
            +1.49638662e-11 * tc[3]
            -1.00127688e-15 * tc[4];
        /*species 1: H */
        species[1] =
            +2.54736599e+04 / tc[1]
            +2.94668292e+00
            -2.50000001e+00 * tc[0]
            +1.15421486e-11 * tc[1]
            -2.69269913e-15 * tc[2]
            +3.94596029e-19 * tc[3]
            -2.49098679e-23 * tc[4];
        /*species 2: O */
        species[2] =
            +2.92175791e+04 / tc[1]
            -2.21491786e+00
            -2.56942078e+00 * tc[0]
            +4.29870569e-05 * tc[1]
            -6.99140982e-09 * tc[2]
            +8.34814992e-13 * tc[3]
            -6.14168455e-17 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.08845772e+03 / tc[1]
            -2.17069345e+00
            -3.28253784e+00 * tc[0]
            -7.41543770e-04 * tc[1]
            +1.26327778e-07 * tc[2]
            -1.74558796e-11 * tc[3]
            +1.08358897e-15 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.85865700e+03 / tc[1]
            -1.38380843e+00
            -3.09288767e+00 * tc[0]
            -2.74214858e-04 * tc[1]
            -2.10842047e-08 * tc[2]
            +7.32884630e-12 * tc[3]
            -5.87061880e-16 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.00042971e+04 / tc[1]
            -1.93277761e+00
            -3.03399249e+00 * tc[0]
            -1.08845902e-03 * tc[1]
            +2.73454197e-08 * tc[2]
            +8.08683225e-12 * tc[3]
            -8.41004960e-16 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.11856713e+02 / tc[1]
            +2.32108750e-01
            -4.01721090e+00 * tc[0]
            -1.11991006e-03 * tc[1]
            +1.05609692e-07 * tc[2]
            -9.52053083e-12 * tc[3]
            +5.39542675e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.78617877e+04 / tc[1]
            +1.24884623e+00
            -4.16500285e+00 * tc[0]
            -2.45415847e-03 * tc[1]
            +3.16898708e-07 * tc[2]
            -3.09321655e-11 * tc[3]
            +1.43954153e-15 * tc[4];
        /*species 8: C */
        species[8] =
            +8.54512953e+04 / tc[1]
            -2.30883485e+00
            -2.49266888e+00 * tc[0]
            -2.39944642e-05 * tc[1]
            +1.20722503e-08 * tc[2]
            -3.11909191e-12 * tc[3]
            +2.43638946e-16 * tc[4];
        /*species 9: CH */
        species[9] =
            +7.10124364e+04 / tc[1]
            -2.60651526e+00
            -2.87846473e+00 * tc[0]
            -4.85456840e-04 * tc[1]
            -2.40742758e-08 * tc[2]
            +1.08906541e-11 * tc[3]
            -8.80396915e-16 * tc[4];
        /*species 10: CH2 */
        species[10] =
            +4.62636040e+04 / tc[1]
            -3.29709211e+00
            -2.87410113e+00 * tc[0]
            -1.82819646e-03 * tc[1]
            +2.34824328e-07 * tc[2]
            -2.16816291e-11 * tc[3]
            +9.38637835e-16 * tc[4];
        /*species 11: CH2(S) */
        species[11] =
            +5.09259997e+04 / tc[1]
            -6.33446327e+00
            -2.29203842e+00 * tc[0]
            -2.32794318e-03 * tc[1]
            +3.35319912e-07 * tc[2]
            -3.48255000e-11 * tc[3]
            +1.69858182e-15 * tc[4];
        /*species 12: CH3 */
        species[12] =
            +1.67755843e+04 / tc[1]
            -6.19435407e+00
            -2.28571772e+00 * tc[0]
            -3.61995018e-03 * tc[1]
            +4.97857247e-07 * tc[2]
            -4.96403870e-11 * tc[3]
            +2.33577197e-15 * tc[4];
        /*species 13: CH4 */
        species[13] =
            -9.46834459e+03 / tc[1]
            -1.83624665e+01
            -7.48514950e-02 * tc[0]
            -6.69547335e-03 * tc[1]
            +9.55476348e-07 * tc[2]
            -1.01910446e-10 * tc[3]
            +5.09076150e-15 * tc[4];
        /*species 14: CO */
        species[14] =
            -1.41518724e+04 / tc[1]
            -5.10350211e+00
            -2.71518561e+00 * tc[0]
            -1.03126372e-03 * tc[1]
            +1.66470962e-07 * tc[2]
            -1.91710840e-11 * tc[3]
            +1.01823858e-15 * tc[4];
        /*species 15: CO2 */
        species[15] =
            -4.87591660e+04 / tc[1]
            +1.58582223e+00
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
        /*species 16: HCO */
        species[16] =
            +4.01191815e+03 / tc[1]
            -7.02617054e+00
            -2.77217438e+00 * tc[0]
            -2.47847763e-03 * tc[1]
            +4.14076022e-07 * tc[2]
            -4.90968148e-11 * tc[3]
            +2.66754356e-15 * tc[4];
        /*species 17: CH2O */
        species[17] =
            -1.39958323e+04 / tc[1]
            -1.18956329e+01
            -1.76069008e+00 * tc[0]
            -4.60000041e-03 * tc[1]
            +7.37098022e-07 * tc[2]
            -8.38676767e-11 * tc[3]
            +4.41927820e-15 * tc[4];
        /*species 18: CH2OH */
        species[18] =
            -3.24250627e+03 / tc[1]
            -2.11776646e+00
            -3.69266569e+00 * tc[0]
            -4.32288399e-03 * tc[1]
            +6.25168533e-07 * tc[2]
            -6.56028863e-11 * tc[3]
            +3.24277101e-15 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +1.27832520e+02 / tc[1]
            +8.41224000e-01
            -3.77079900e+00 * tc[0]
            -3.93574850e-03 * tc[1]
            +4.42730667e-07 * tc[2]
            -3.28702583e-11 * tc[3]
            +1.05630800e-15 * tc[4];
        /*species 20: CH3OH */
        species[20] =
            -2.53748747e+04 / tc[1]
            -1.27126544e+01
            -1.78970791e+00 * tc[0]
            -7.04691460e-03 * tc[1]
            +1.06083472e-06 * tc[2]
            -1.15142571e-10 * tc[3]
            +5.85301100e-15 * tc[4];
        /*species 21: C2H */
        species[21] =
            +6.71210650e+04 / tc[1]
            -3.46808823e+00
            -3.16780652e+00 * tc[0]
            -2.37610951e-03 * tc[1]
            +3.06311795e-07 * tc[2]
            -2.53491877e-11 * tc[3]
            +8.86163850e-16 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            +2.59359992e+04 / tc[1]
            +5.37785085e+00
            -4.14756964e+00 * tc[0]
            -2.98083332e-03 * tc[1]
            +3.95491420e-07 * tc[2]
            -3.89510143e-11 * tc[3]
            +1.80617607e-15 * tc[4];
        /*species 23: C2H3 */
        species[23] =
            +3.46128739e+04 / tc[1]
            -4.77059978e+00
            -3.01672400e+00 * tc[0]
            -5.16511460e-03 * tc[1]
            +7.80137248e-07 * tc[2]
            -8.48027400e-11 * tc[3]
            +4.31303520e-15 * tc[4];
        /*species 24: C2H4 */
        species[24] =
            +4.93988614e+03 / tc[1]
            -8.26925814e+00
            -2.03611116e+00 * tc[0]
            -7.32270755e-03 * tc[1]
            +1.11846319e-06 * tc[2]
            -1.22685769e-10 * tc[3]
            +6.28530305e-15 * tc[4];
        /*species 25: C2H5 */
        species[25] =
            +1.28575200e+04 / tc[1]
            -1.15077779e+01
            -1.95465642e+00 * tc[0]
            -8.69863610e-03 * tc[1]
            +1.33034445e-06 * tc[2]
            -1.46014741e-10 * tc[3]
            +7.48207880e-15 * tc[4];
        /*species 26: C2H6 */
        species[26] =
            -1.14263932e+04 / tc[1]
            -1.40437292e+01
            -1.07188150e+00 * tc[0]
            -1.08426339e-02 * tc[1]
            +1.67093445e-06 * tc[2]
            -1.84510001e-10 * tc[3]
            +9.50014450e-15 * tc[4];
        /*species 27: HCCO */
        species[27] =
            +1.93272150e+04 / tc[1]
            +9.55846530e+00
            -5.62820580e+00 * tc[0]
            -2.04267005e-03 * tc[1]
            +2.65575783e-07 * tc[2]
            -2.38550433e-11 * tc[3]
            +9.70391600e-16 * tc[4];
        /*species 28: CH2CO */
        species[28] =
            -7.55105311e+03 / tc[1]
            +3.87905011e+00
            -4.51129732e+00 * tc[0]
            -4.50179872e-03 * tc[1]
            +6.94899392e-07 * tc[2]
            -7.69454902e-11 * tc[3]
            +3.97419100e-15 * tc[4];
        /*species 29: HCCOH */
        species[29] =
            +7.26462600e+03 / tc[1]
            +1.35256033e+01
            -5.92382910e+00 * tc[0]
            -3.39618000e-03 * tc[1]
            +4.27642733e-07 * tc[2]
            -3.74898675e-11 * tc[3]
            +1.49700505e-15 * tc[4];
        /*species 30: N2 */
        species[30] =
            -9.22797700e+02 / tc[1]
            -3.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 31: AR */
        species[31] =
            -7.45375000e+02 / tc[1]
            -1.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
    }
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
        /*species 0: H2 */
        species[0] =
            -9.17935173e+02 / tc[1]
            +6.61320882e-01
            -2.34433112e+00 * tc[0]
            -3.99026037e-03 * tc[1]
            +3.24635850e-06 * tc[2]
            -1.67976745e-09 * tc[3]
            +3.68805881e-13 * tc[4];
        /*species 1: H */
        species[1] =
            +2.54736599e+04 / tc[1]
            +1.94668285e+00
            -2.50000000e+00 * tc[0]
            -3.52666409e-13 * tc[1]
            +3.32653273e-16 * tc[2]
            -1.91734693e-19 * tc[3]
            +4.63866166e-23 * tc[4];
        /*species 2: O */
        species[2] =
            +2.91222592e+04 / tc[1]
            +1.16333640e-01
            -3.16826710e+00 * tc[0]
            +1.63965942e-03 * tc[1]
            -1.10717733e-06 * tc[2]
            +5.10672187e-10 * tc[3]
            -1.05632985e-13 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.06394356e+03 / tc[1]
            -8.75219370e-01
            -3.78245636e+00 * tc[0]
            +1.49836708e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774591e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.61508056e+03 / tc[1]
            +3.09594089e+00
            -3.99201543e+00 * tc[0]
            +1.20065876e-03 * tc[1]
            -7.69656402e-07 * tc[2]
            +3.23427778e-10 * tc[3]
            -6.82057350e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.02937267e+04 / tc[1]
            +4.04767277e+00
            -4.19864056e+00 * tc[0]
            +1.01821705e-03 * tc[1]
            -1.08673369e-06 * tc[2]
            +4.57330885e-10 * tc[3]
            -8.85989085e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +2.94808040e+02 / tc[1]
            -4.14864440e-01
            -4.30179801e+00 * tc[0]
            +2.37456025e-03 * tc[1]
            -3.52638152e-06 * tc[2]
            +2.02303245e-09 * tc[3]
            -4.64612562e-13 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.77025821e+04 / tc[1]
            -1.58938050e-01
            -4.27611269e+00 * tc[0]
            +2.71411208e-04 * tc[1]
            -2.78892835e-06 * tc[2]
            +1.79809011e-09 * tc[3]
            -4.31227182e-13 * tc[4];
        /*species 8: C */
        species[8] =
            +8.54438832e+04 / tc[1]
            -2.97706893e+00
            -2.55423955e+00 * tc[0]
            +1.60768862e-04 * tc[1]
            -1.22298707e-07 * tc[2]
            +6.10195741e-11 * tc[3]
            -1.33260723e-14 * tc[4];
        /*species 9: CH */
        species[9] =
            +7.07972934e+04 / tc[1]
            +4.05805570e-01
            -3.48981665e+00 * tc[0]
            -1.61917771e-04 * tc[1]
            +2.81498442e-07 * tc[2]
            -2.63514439e-10 * tc[3]
            +7.03045335e-14 * tc[4];
        /*species 10: CH2 */
        species[10] =
            +4.60040401e+04 / tc[1]
            +1.20014682e+00
            -3.76267867e+00 * tc[0]
            -4.84436072e-04 * tc[1]
            -4.65816402e-07 * tc[2]
            +3.20909294e-10 * tc[3]
            -8.43708595e-14 * tc[4];
        /*species 11: CH2(S) */
        species[11] =
            +5.04968163e+04 / tc[1]
            +3.96772308e+00
            -4.19860411e+00 * tc[0]
            +1.18330710e-03 * tc[1]
            -1.37216037e-06 * tc[2]
            +5.57346651e-10 * tc[3]
            -9.71573685e-14 * tc[4];
        /*species 12: CH3 */
        species[12] =
            +1.64449988e+04 / tc[1]
            +1.06902607e+00
            -3.67359040e+00 * tc[0]
            -1.00547588e-03 * tc[1]
            -9.55036427e-07 * tc[2]
            +5.72597854e-10 * tc[3]
            -1.27192867e-13 * tc[4];
        /*species 13: CH4 */
        species[13] =
            -1.02466476e+04 / tc[1]
            +8.79117989e+00
            -5.14987613e+00 * tc[0]
            +6.83548940e-03 * tc[1]
            -8.19667665e-06 * tc[2]
            +4.03952522e-09 * tc[3]
            -8.33469780e-13 * tc[4];
        /*species 14: CO */
        species[14] =
            -1.43440860e+04 / tc[1]
            -9.28875810e-01
            -3.57953347e+00 * tc[0]
            +3.05176840e-04 * tc[1]
            -1.69469055e-07 * tc[2]
            -7.55838237e-11 * tc[3]
            +4.52212249e-14 * tc[4];
        /*species 15: CO2 */
        species[15] =
            -4.83719697e+04 / tc[1]
            -8.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
        /*species 16: HCO */
        species[16] =
            +3.83956496e+03 / tc[1]
            -1.73186590e-01
            -4.22118584e+00 * tc[0]
            +1.62196266e-03 * tc[1]
            -2.29665743e-06 * tc[2]
            +1.10953411e-09 * tc[3]
            -2.16884432e-13 * tc[4];
        /*species 17: CH2O */
        species[17] =
            -1.43089567e+04 / tc[1]
            +3.19091025e+00
            -4.79372315e+00 * tc[0]
            +4.95416684e-03 * tc[1]
            -6.22033347e-06 * tc[2]
            +3.16071051e-09 * tc[3]
            -6.58863260e-13 * tc[4];
        /*species 18: CH2OH */
        species[18] =
            -3.19391367e+03 / tc[1]
            -2.60913325e+00
            -3.86388918e+00 * tc[0]
            -2.79836152e-03 * tc[1]
            -9.88786318e-07 * tc[2]
            +8.71100100e-10 * tc[3]
            -2.18483639e-13 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +9.78601100e+02 / tc[1]
            -1.20459730e+01
            -2.10620400e+00 * tc[0]
            -3.60829750e-03 * tc[1]
            -8.89745333e-07 * tc[2]
            +6.14803000e-10 * tc[3]
            -1.03780500e-13 * tc[4];
        /*species 20: CH3OH */
        species[20] =
            -2.56427656e+04 / tc[1]
            +6.21949405e+00
            -5.71539582e+00 * tc[0]
            +7.61545645e-03 * tc[1]
            -1.08740193e-05 * tc[2]
            +5.92339074e-09 * tc[3]
            -1.30676349e-12 * tc[4];
        /*species 21: C2H */
        species[21] =
            +6.68393932e+04 / tc[1]
            -4.33330705e+00
            -2.88965733e+00 * tc[0]
            -6.70498055e-03 * tc[1]
            +4.74615835e-06 * tc[2]
            -2.45659204e-09 * tc[3]
            +5.46657555e-13 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            +2.64289807e+04 / tc[1]
            -1.41310240e+01
            -8.08681094e-01 * tc[0]
            -1.16807815e-02 * tc[1]
            +5.91953025e-06 * tc[2]
            -2.33460364e-09 * tc[3]
            +4.25036487e-13 * tc[4];
        /*species 23: C2H3 */
        species[23] =
            +3.48598468e+04 / tc[1]
            -6.29807380e+00
            -3.21246645e+00 * tc[0]
            -7.57395810e-04 * tc[1]
            -4.32015687e-06 * tc[2]
            +2.98048206e-09 * tc[3]
            -7.35754365e-13 * tc[4];
        /*species 24: C2H4 */
        species[24] =
            +5.08977593e+03 / tc[1]
            -1.13812948e+00
            -3.95920148e+00 * tc[0]
            +3.78526124e-03 * tc[1]
            -9.51650487e-06 * tc[2]
            +5.76323961e-09 * tc[3]
            -1.34942187e-12 * tc[4];
        /*species 25: C2H5 */
        species[25] =
            +1.28416265e+04 / tc[1]
            -1.40074356e+00
            -4.30646568e+00 * tc[0]
            +2.09329446e-03 * tc[1]
            -8.28571345e-06 * tc[2]
            +4.99272172e-09 * tc[3]
            -1.15254502e-12 * tc[4];
        /*species 26: C2H6 */
        species[26] =
            -1.15222055e+04 / tc[1]
            +6.24601760e-01
            -4.29142492e+00 * tc[0]
            +2.75077135e-03 * tc[1]
            -9.99063813e-06 * tc[2]
            +5.90388571e-09 * tc[3]
            -1.34342886e-12 * tc[4];
        /*species 27: HCCO */
        species[27] =
            +2.00594490e+04 / tc[1]
            -1.12386956e+01
            -2.25172140e+00 * tc[0]
            -8.82751050e-03 * tc[1]
            +3.95485017e-06 * tc[2]
            -1.43964658e-09 * tc[3]
            +2.53324055e-13 * tc[4];
        /*species 28: CH2CO */
        species[28] =
            -7.04291804e+03 / tc[1]
            -1.10798117e+01
            -2.13583630e+00 * tc[0]
            -9.05943605e-03 * tc[1]
            +2.89912457e-06 * tc[2]
            -7.78664640e-10 * tc[3]
            +1.00728807e-13 * tc[4];
        /*species 29: HCCOH */
        species[29] =
            +8.03161430e+03 / tc[1]
            -1.36319457e+01
            -1.24237330e+00 * tc[0]
            -1.55361005e-02 * tc[1]
            +8.47781067e-06 * tc[2]
            -3.59476092e-09 * tc[3]
            +7.00729700e-13 * tc[4];
        /*species 30: N2 */
        species[30] =
            -1.02089990e+03 / tc[1]
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
        /*species 31: AR */
        species[31] =
            -7.45375000e+02 / tc[1]
            -2.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -9.50158922e+02 / tc[1]
            +5.54230251e+00
            -3.33727920e+00 * tc[0]
            +2.47012365e-05 * tc[1]
            -8.32427963e-08 * tc[2]
            +1.49638662e-11 * tc[3]
            -1.00127688e-15 * tc[4];
        /*species 1: H */
        species[1] =
            +2.54736599e+04 / tc[1]
            +1.94668292e+00
            -2.50000001e+00 * tc[0]
            +1.15421486e-11 * tc[1]
            -2.69269913e-15 * tc[2]
            +3.94596029e-19 * tc[3]
            -2.49098679e-23 * tc[4];
        /*species 2: O */
        species[2] =
            +2.92175791e+04 / tc[1]
            -3.21491786e+00
            -2.56942078e+00 * tc[0]
            +4.29870569e-05 * tc[1]
            -6.99140982e-09 * tc[2]
            +8.34814992e-13 * tc[3]
            -6.14168455e-17 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.08845772e+03 / tc[1]
            -3.17069345e+00
            -3.28253784e+00 * tc[0]
            -7.41543770e-04 * tc[1]
            +1.26327778e-07 * tc[2]
            -1.74558796e-11 * tc[3]
            +1.08358897e-15 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.85865700e+03 / tc[1]
            -2.38380843e+00
            -3.09288767e+00 * tc[0]
            -2.74214858e-04 * tc[1]
            -2.10842047e-08 * tc[2]
            +7.32884630e-12 * tc[3]
            -5.87061880e-16 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.00042971e+04 / tc[1]
            -2.93277761e+00
            -3.03399249e+00 * tc[0]
            -1.08845902e-03 * tc[1]
            +2.73454197e-08 * tc[2]
            +8.08683225e-12 * tc[3]
            -8.41004960e-16 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.11856713e+02 / tc[1]
            -7.67891250e-01
            -4.01721090e+00 * tc[0]
            -1.11991006e-03 * tc[1]
            +1.05609692e-07 * tc[2]
            -9.52053083e-12 * tc[3]
            +5.39542675e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.78617877e+04 / tc[1]
            +2.48846230e-01
            -4.16500285e+00 * tc[0]
            -2.45415847e-03 * tc[1]
            +3.16898708e-07 * tc[2]
            -3.09321655e-11 * tc[3]
            +1.43954153e-15 * tc[4];
        /*species 8: C */
        species[8] =
            +8.54512953e+04 / tc[1]
            -3.30883485e+00
            -2.49266888e+00 * tc[0]
            -2.39944642e-05 * tc[1]
            +1.20722503e-08 * tc[2]
            -3.11909191e-12 * tc[3]
            +2.43638946e-16 * tc[4];
        /*species 9: CH */
        species[9] =
            +7.10124364e+04 / tc[1]
            -3.60651526e+00
            -2.87846473e+00 * tc[0]
            -4.85456840e-04 * tc[1]
            -2.40742758e-08 * tc[2]
            +1.08906541e-11 * tc[3]
            -8.80396915e-16 * tc[4];
        /*species 10: CH2 */
        species[10] =
            +4.62636040e+04 / tc[1]
            -4.29709211e+00
            -2.87410113e+00 * tc[0]
            -1.82819646e-03 * tc[1]
            +2.34824328e-07 * tc[2]
            -2.16816291e-11 * tc[3]
            +9.38637835e-16 * tc[4];
        /*species 11: CH2(S) */
        species[11] =
            +5.09259997e+04 / tc[1]
            -7.33446327e+00
            -2.29203842e+00 * tc[0]
            -2.32794318e-03 * tc[1]
            +3.35319912e-07 * tc[2]
            -3.48255000e-11 * tc[3]
            +1.69858182e-15 * tc[4];
        /*species 12: CH3 */
        species[12] =
            +1.67755843e+04 / tc[1]
            -7.19435407e+00
            -2.28571772e+00 * tc[0]
            -3.61995018e-03 * tc[1]
            +4.97857247e-07 * tc[2]
            -4.96403870e-11 * tc[3]
            +2.33577197e-15 * tc[4];
        /*species 13: CH4 */
        species[13] =
            -9.46834459e+03 / tc[1]
            -1.93624665e+01
            -7.48514950e-02 * tc[0]
            -6.69547335e-03 * tc[1]
            +9.55476348e-07 * tc[2]
            -1.01910446e-10 * tc[3]
            +5.09076150e-15 * tc[4];
        /*species 14: CO */
        species[14] =
            -1.41518724e+04 / tc[1]
            -6.10350211e+00
            -2.71518561e+00 * tc[0]
            -1.03126372e-03 * tc[1]
            +1.66470962e-07 * tc[2]
            -1.91710840e-11 * tc[3]
            +1.01823858e-15 * tc[4];
        /*species 15: CO2 */
        species[15] =
            -4.87591660e+04 / tc[1]
            +5.85822230e-01
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
        /*species 16: HCO */
        species[16] =
            +4.01191815e+03 / tc[1]
            -8.02617054e+00
            -2.77217438e+00 * tc[0]
            -2.47847763e-03 * tc[1]
            +4.14076022e-07 * tc[2]
            -4.90968148e-11 * tc[3]
            +2.66754356e-15 * tc[4];
        /*species 17: CH2O */
        species[17] =
            -1.39958323e+04 / tc[1]
            -1.28956329e+01
            -1.76069008e+00 * tc[0]
            -4.60000041e-03 * tc[1]
            +7.37098022e-07 * tc[2]
            -8.38676767e-11 * tc[3]
            +4.41927820e-15 * tc[4];
        /*species 18: CH2OH */
        species[18] =
            -3.24250627e+03 / tc[1]
            -3.11776646e+00
            -3.69266569e+00 * tc[0]
            -4.32288399e-03 * tc[1]
            +6.25168533e-07 * tc[2]
            -6.56028863e-11 * tc[3]
            +3.24277101e-15 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +1.27832520e+02 / tc[1]
            -1.58776000e-01
            -3.77079900e+00 * tc[0]
            -3.93574850e-03 * tc[1]
            +4.42730667e-07 * tc[2]
            -3.28702583e-11 * tc[3]
            +1.05630800e-15 * tc[4];
        /*species 20: CH3OH */
        species[20] =
            -2.53748747e+04 / tc[1]
            -1.37126544e+01
            -1.78970791e+00 * tc[0]
            -7.04691460e-03 * tc[1]
            +1.06083472e-06 * tc[2]
            -1.15142571e-10 * tc[3]
            +5.85301100e-15 * tc[4];
        /*species 21: C2H */
        species[21] =
            +6.71210650e+04 / tc[1]
            -4.46808823e+00
            -3.16780652e+00 * tc[0]
            -2.37610951e-03 * tc[1]
            +3.06311795e-07 * tc[2]
            -2.53491877e-11 * tc[3]
            +8.86163850e-16 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            +2.59359992e+04 / tc[1]
            +4.37785085e+00
            -4.14756964e+00 * tc[0]
            -2.98083332e-03 * tc[1]
            +3.95491420e-07 * tc[2]
            -3.89510143e-11 * tc[3]
            +1.80617607e-15 * tc[4];
        /*species 23: C2H3 */
        species[23] =
            +3.46128739e+04 / tc[1]
            -5.77059978e+00
            -3.01672400e+00 * tc[0]
            -5.16511460e-03 * tc[1]
            +7.80137248e-07 * tc[2]
            -8.48027400e-11 * tc[3]
            +4.31303520e-15 * tc[4];
        /*species 24: C2H4 */
        species[24] =
            +4.93988614e+03 / tc[1]
            -9.26925814e+00
            -2.03611116e+00 * tc[0]
            -7.32270755e-03 * tc[1]
            +1.11846319e-06 * tc[2]
            -1.22685769e-10 * tc[3]
            +6.28530305e-15 * tc[4];
        /*species 25: C2H5 */
        species[25] =
            +1.28575200e+04 / tc[1]
            -1.25077779e+01
            -1.95465642e+00 * tc[0]
            -8.69863610e-03 * tc[1]
            +1.33034445e-06 * tc[2]
            -1.46014741e-10 * tc[3]
            +7.48207880e-15 * tc[4];
        /*species 26: C2H6 */
        species[26] =
            -1.14263932e+04 / tc[1]
            -1.50437292e+01
            -1.07188150e+00 * tc[0]
            -1.08426339e-02 * tc[1]
            +1.67093445e-06 * tc[2]
            -1.84510001e-10 * tc[3]
            +9.50014450e-15 * tc[4];
        /*species 27: HCCO */
        species[27] =
            +1.93272150e+04 / tc[1]
            +8.55846530e+00
            -5.62820580e+00 * tc[0]
            -2.04267005e-03 * tc[1]
            +2.65575783e-07 * tc[2]
            -2.38550433e-11 * tc[3]
            +9.70391600e-16 * tc[4];
        /*species 28: CH2CO */
        species[28] =
            -7.55105311e+03 / tc[1]
            +2.87905011e+00
            -4.51129732e+00 * tc[0]
            -4.50179872e-03 * tc[1]
            +6.94899392e-07 * tc[2]
            -7.69454902e-11 * tc[3]
            +3.97419100e-15 * tc[4];
        /*species 29: HCCOH */
        species[29] =
            +7.26462600e+03 / tc[1]
            +1.25256033e+01
            -5.92382910e+00 * tc[0]
            -3.39618000e-03 * tc[1]
            +4.27642733e-07 * tc[2]
            -3.74898675e-11 * tc[3]
            +1.49700505e-15 * tc[4];
        /*species 30: N2 */
        species[30] =
            -9.22797700e+02 / tc[1]
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 31: AR */
        species[31] =
            -7.45375000e+02 / tc[1]
            -2.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
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
        /*species 0: H2 */
        species[0] =
            +1.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 1: H */
        species[1] =
            +1.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
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
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 4: OH */
        species[4] =
            +2.99201543e+00
            -2.40131752e-03 * tc[1]
            +4.61793841e-06 * tc[2]
            -3.88113333e-09 * tc[3]
            +1.36411470e-12 * tc[4];
        /*species 5: H2O */
        species[5] =
            +3.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 8: C */
        species[8] =
            +1.55423955e+00
            -3.21537724e-04 * tc[1]
            +7.33792245e-07 * tc[2]
            -7.32234889e-10 * tc[3]
            +2.66521446e-13 * tc[4];
        /*species 9: CH */
        species[9] =
            +2.48981665e+00
            +3.23835541e-04 * tc[1]
            -1.68899065e-06 * tc[2]
            +3.16217327e-09 * tc[3]
            -1.40609067e-12 * tc[4];
        /*species 10: CH2 */
        species[10] =
            +2.76267867e+00
            +9.68872143e-04 * tc[1]
            +2.79489841e-06 * tc[2]
            -3.85091153e-09 * tc[3]
            +1.68741719e-12 * tc[4];
        /*species 11: CH2(S) */
        species[11] =
            +3.19860411e+00
            -2.36661419e-03 * tc[1]
            +8.23296220e-06 * tc[2]
            -6.68815981e-09 * tc[3]
            +1.94314737e-12 * tc[4];
        /*species 12: CH3 */
        species[12] =
            +2.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 13: CH4 */
        species[13] =
            +4.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 14: CO */
        species[14] =
            +2.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 15: CO2 */
        species[15] =
            +1.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 16: HCO */
        species[16] =
            +3.22118584e+00
            -3.24392532e-03 * tc[1]
            +1.37799446e-05 * tc[2]
            -1.33144093e-08 * tc[3]
            +4.33768865e-12 * tc[4];
        /*species 17: CH2O */
        species[17] =
            +3.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 18: CH2OH */
        species[18] =
            +2.86388918e+00
            +5.59672304e-03 * tc[1]
            +5.93271791e-06 * tc[2]
            -1.04532012e-08 * tc[3]
            +4.36967278e-12 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +1.10620400e+00
            +7.21659500e-03 * tc[1]
            +5.33847200e-06 * tc[2]
            -7.37763600e-09 * tc[3]
            +2.07561000e-12 * tc[4];
        /*species 20: CH3OH */
        species[20] =
            +4.71539582e+00
            -1.52309129e-02 * tc[1]
            +6.52441155e-05 * tc[2]
            -7.10806889e-08 * tc[3]
            +2.61352698e-11 * tc[4];
        /*species 21: C2H */
        species[21] =
            +1.88965733e+00
            +1.34099611e-02 * tc[1]
            -2.84769501e-05 * tc[2]
            +2.94791045e-08 * tc[3]
            -1.09331511e-11 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            -1.91318906e-01
            +2.33615629e-02 * tc[1]
            -3.55171815e-05 * tc[2]
            +2.80152437e-08 * tc[3]
            -8.50072974e-12 * tc[4];
        /*species 23: C2H3 */
        species[23] =
            +2.21246645e+00
            +1.51479162e-03 * tc[1]
            +2.59209412e-05 * tc[2]
            -3.57657847e-08 * tc[3]
            +1.47150873e-11 * tc[4];
        /*species 24: C2H4 */
        species[24] =
            +2.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 25: C2H5 */
        species[25] =
            +3.30646568e+00
            -4.18658892e-03 * tc[1]
            +4.97142807e-05 * tc[2]
            -5.99126606e-08 * tc[3]
            +2.30509004e-11 * tc[4];
        /*species 26: C2H6 */
        species[26] =
            +3.29142492e+00
            -5.50154270e-03 * tc[1]
            +5.99438288e-05 * tc[2]
            -7.08466285e-08 * tc[3]
            +2.68685771e-11 * tc[4];
        /*species 27: HCCO */
        species[27] =
            +1.25172140e+00
            +1.76550210e-02 * tc[1]
            -2.37291010e-05 * tc[2]
            +1.72757590e-08 * tc[3]
            -5.06648110e-12 * tc[4];
        /*species 28: CH2CO */
        species[28] =
            +1.13583630e+00
            +1.81188721e-02 * tc[1]
            -1.73947474e-05 * tc[2]
            +9.34397568e-09 * tc[3]
            -2.01457615e-12 * tc[4];
        /*species 29: HCCOH */
        species[29] =
            +2.42373300e-01
            +3.10722010e-02 * tc[1]
            -5.08668640e-05 * tc[2]
            +4.31371310e-08 * tc[3]
            -1.40145940e-11 * tc[4];
        /*species 30: N2 */
        species[30] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 31: AR */
        species[31] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +2.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 1: H */
        species[1] =
            +1.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 2: O */
        species[2] =
            +1.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 3: O2 */
        species[3] =
            +2.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 4: OH */
        species[4] =
            +2.09288767e+00
            +5.48429716e-04 * tc[1]
            +1.26505228e-07 * tc[2]
            -8.79461556e-11 * tc[3]
            +1.17412376e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            +2.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 8: C */
        species[8] =
            +1.49266888e+00
            +4.79889284e-05 * tc[1]
            -7.24335020e-08 * tc[2]
            +3.74291029e-11 * tc[3]
            -4.87277893e-15 * tc[4];
        /*species 9: CH */
        species[9] =
            +1.87846473e+00
            +9.70913681e-04 * tc[1]
            +1.44445655e-07 * tc[2]
            -1.30687849e-10 * tc[3]
            +1.76079383e-14 * tc[4];
        /*species 10: CH2 */
        species[10] =
            +1.87410113e+00
            +3.65639292e-03 * tc[1]
            -1.40894597e-06 * tc[2]
            +2.60179549e-10 * tc[3]
            -1.87727567e-14 * tc[4];
        /*species 11: CH2(S) */
        species[11] =
            +1.29203842e+00
            +4.65588637e-03 * tc[1]
            -2.01191947e-06 * tc[2]
            +4.17906000e-10 * tc[3]
            -3.39716365e-14 * tc[4];
        /*species 12: CH3 */
        species[12] =
            +1.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 13: CH4 */
        species[13] =
            -9.25148505e-01
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 14: CO */
        species[14] =
            +1.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 15: CO2 */
        species[15] =
            +2.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 16: HCO */
        species[16] =
            +1.77217438e+00
            +4.95695526e-03 * tc[1]
            -2.48445613e-06 * tc[2]
            +5.89161778e-10 * tc[3]
            -5.33508711e-14 * tc[4];
        /*species 17: CH2O */
        species[17] =
            +7.60690080e-01
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 18: CH2OH */
        species[18] =
            +2.69266569e+00
            +8.64576797e-03 * tc[1]
            -3.75101120e-06 * tc[2]
            +7.87234636e-10 * tc[3]
            -6.48554201e-14 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +2.77079900e+00
            +7.87149700e-03 * tc[1]
            -2.65638400e-06 * tc[2]
            +3.94443100e-10 * tc[3]
            -2.11261600e-14 * tc[4];
        /*species 20: CH3OH */
        species[20] =
            +7.89707910e-01
            +1.40938292e-02 * tc[1]
            -6.36500835e-06 * tc[2]
            +1.38171085e-09 * tc[3]
            -1.17060220e-13 * tc[4];
        /*species 21: C2H */
        species[21] =
            +2.16780652e+00
            +4.75221902e-03 * tc[1]
            -1.83787077e-06 * tc[2]
            +3.04190252e-10 * tc[3]
            -1.77232770e-14 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            +3.14756964e+00
            +5.96166664e-03 * tc[1]
            -2.37294852e-06 * tc[2]
            +4.67412171e-10 * tc[3]
            -3.61235213e-14 * tc[4];
        /*species 23: C2H3 */
        species[23] =
            +2.01672400e+00
            +1.03302292e-02 * tc[1]
            -4.68082349e-06 * tc[2]
            +1.01763288e-09 * tc[3]
            -8.62607041e-14 * tc[4];
        /*species 24: C2H4 */
        species[24] =
            +1.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 25: C2H5 */
        species[25] =
            +9.54656420e-01
            +1.73972722e-02 * tc[1]
            -7.98206668e-06 * tc[2]
            +1.75217689e-09 * tc[3]
            -1.49641576e-13 * tc[4];
        /*species 26: C2H6 */
        species[26] =
            +7.18815000e-02
            +2.16852677e-02 * tc[1]
            -1.00256067e-05 * tc[2]
            +2.21412001e-09 * tc[3]
            -1.90002890e-13 * tc[4];
        /*species 27: HCCO */
        species[27] =
            +4.62820580e+00
            +4.08534010e-03 * tc[1]
            -1.59345470e-06 * tc[2]
            +2.86260520e-10 * tc[3]
            -1.94078320e-14 * tc[4];
        /*species 28: CH2CO */
        species[28] =
            +3.51129732e+00
            +9.00359745e-03 * tc[1]
            -4.16939635e-06 * tc[2]
            +9.23345882e-10 * tc[3]
            -7.94838201e-14 * tc[4];
        /*species 29: HCCOH */
        species[29] =
            +4.92382910e+00
            +6.79236000e-03 * tc[1]
            -2.56585640e-06 * tc[2]
            +4.49878410e-10 * tc[3]
            -2.99401010e-14 * tc[4];
        /*species 30: N2 */
        species[30] =
            +1.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 31: AR */
        species[31] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
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
        /*species 0: H2 */
        species[0] =
            +2.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 1: H */
        species[1] =
            +2.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
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
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.99201543e+00
            -2.40131752e-03 * tc[1]
            +4.61793841e-06 * tc[2]
            -3.88113333e-09 * tc[3]
            +1.36411470e-12 * tc[4];
        /*species 5: H2O */
        species[5] =
            +4.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 8: C */
        species[8] =
            +2.55423955e+00
            -3.21537724e-04 * tc[1]
            +7.33792245e-07 * tc[2]
            -7.32234889e-10 * tc[3]
            +2.66521446e-13 * tc[4];
        /*species 9: CH */
        species[9] =
            +3.48981665e+00
            +3.23835541e-04 * tc[1]
            -1.68899065e-06 * tc[2]
            +3.16217327e-09 * tc[3]
            -1.40609067e-12 * tc[4];
        /*species 10: CH2 */
        species[10] =
            +3.76267867e+00
            +9.68872143e-04 * tc[1]
            +2.79489841e-06 * tc[2]
            -3.85091153e-09 * tc[3]
            +1.68741719e-12 * tc[4];
        /*species 11: CH2(S) */
        species[11] =
            +4.19860411e+00
            -2.36661419e-03 * tc[1]
            +8.23296220e-06 * tc[2]
            -6.68815981e-09 * tc[3]
            +1.94314737e-12 * tc[4];
        /*species 12: CH3 */
        species[12] =
            +3.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 13: CH4 */
        species[13] =
            +5.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 14: CO */
        species[14] =
            +3.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 15: CO2 */
        species[15] =
            +2.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 16: HCO */
        species[16] =
            +4.22118584e+00
            -3.24392532e-03 * tc[1]
            +1.37799446e-05 * tc[2]
            -1.33144093e-08 * tc[3]
            +4.33768865e-12 * tc[4];
        /*species 17: CH2O */
        species[17] =
            +4.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 18: CH2OH */
        species[18] =
            +3.86388918e+00
            +5.59672304e-03 * tc[1]
            +5.93271791e-06 * tc[2]
            -1.04532012e-08 * tc[3]
            +4.36967278e-12 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +2.10620400e+00
            +7.21659500e-03 * tc[1]
            +5.33847200e-06 * tc[2]
            -7.37763600e-09 * tc[3]
            +2.07561000e-12 * tc[4];
        /*species 20: CH3OH */
        species[20] =
            +5.71539582e+00
            -1.52309129e-02 * tc[1]
            +6.52441155e-05 * tc[2]
            -7.10806889e-08 * tc[3]
            +2.61352698e-11 * tc[4];
        /*species 21: C2H */
        species[21] =
            +2.88965733e+00
            +1.34099611e-02 * tc[1]
            -2.84769501e-05 * tc[2]
            +2.94791045e-08 * tc[3]
            -1.09331511e-11 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            +8.08681094e-01
            +2.33615629e-02 * tc[1]
            -3.55171815e-05 * tc[2]
            +2.80152437e-08 * tc[3]
            -8.50072974e-12 * tc[4];
        /*species 23: C2H3 */
        species[23] =
            +3.21246645e+00
            +1.51479162e-03 * tc[1]
            +2.59209412e-05 * tc[2]
            -3.57657847e-08 * tc[3]
            +1.47150873e-11 * tc[4];
        /*species 24: C2H4 */
        species[24] =
            +3.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 25: C2H5 */
        species[25] =
            +4.30646568e+00
            -4.18658892e-03 * tc[1]
            +4.97142807e-05 * tc[2]
            -5.99126606e-08 * tc[3]
            +2.30509004e-11 * tc[4];
        /*species 26: C2H6 */
        species[26] =
            +4.29142492e+00
            -5.50154270e-03 * tc[1]
            +5.99438288e-05 * tc[2]
            -7.08466285e-08 * tc[3]
            +2.68685771e-11 * tc[4];
        /*species 27: HCCO */
        species[27] =
            +2.25172140e+00
            +1.76550210e-02 * tc[1]
            -2.37291010e-05 * tc[2]
            +1.72757590e-08 * tc[3]
            -5.06648110e-12 * tc[4];
        /*species 28: CH2CO */
        species[28] =
            +2.13583630e+00
            +1.81188721e-02 * tc[1]
            -1.73947474e-05 * tc[2]
            +9.34397568e-09 * tc[3]
            -2.01457615e-12 * tc[4];
        /*species 29: HCCOH */
        species[29] =
            +1.24237330e+00
            +3.10722010e-02 * tc[1]
            -5.08668640e-05 * tc[2]
            +4.31371310e-08 * tc[3]
            -1.40145940e-11 * tc[4];
        /*species 30: N2 */
        species[30] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 31: AR */
        species[31] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +3.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 1: H */
        species[1] =
            +2.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 2: O */
        species[2] =
            +2.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.09288767e+00
            +5.48429716e-04 * tc[1]
            +1.26505228e-07 * tc[2]
            -8.79461556e-11 * tc[3]
            +1.17412376e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            +3.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 8: C */
        species[8] =
            +2.49266888e+00
            +4.79889284e-05 * tc[1]
            -7.24335020e-08 * tc[2]
            +3.74291029e-11 * tc[3]
            -4.87277893e-15 * tc[4];
        /*species 9: CH */
        species[9] =
            +2.87846473e+00
            +9.70913681e-04 * tc[1]
            +1.44445655e-07 * tc[2]
            -1.30687849e-10 * tc[3]
            +1.76079383e-14 * tc[4];
        /*species 10: CH2 */
        species[10] =
            +2.87410113e+00
            +3.65639292e-03 * tc[1]
            -1.40894597e-06 * tc[2]
            +2.60179549e-10 * tc[3]
            -1.87727567e-14 * tc[4];
        /*species 11: CH2(S) */
        species[11] =
            +2.29203842e+00
            +4.65588637e-03 * tc[1]
            -2.01191947e-06 * tc[2]
            +4.17906000e-10 * tc[3]
            -3.39716365e-14 * tc[4];
        /*species 12: CH3 */
        species[12] =
            +2.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 13: CH4 */
        species[13] =
            +7.48514950e-02
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 14: CO */
        species[14] =
            +2.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 15: CO2 */
        species[15] =
            +3.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 16: HCO */
        species[16] =
            +2.77217438e+00
            +4.95695526e-03 * tc[1]
            -2.48445613e-06 * tc[2]
            +5.89161778e-10 * tc[3]
            -5.33508711e-14 * tc[4];
        /*species 17: CH2O */
        species[17] =
            +1.76069008e+00
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 18: CH2OH */
        species[18] =
            +3.69266569e+00
            +8.64576797e-03 * tc[1]
            -3.75101120e-06 * tc[2]
            +7.87234636e-10 * tc[3]
            -6.48554201e-14 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +3.77079900e+00
            +7.87149700e-03 * tc[1]
            -2.65638400e-06 * tc[2]
            +3.94443100e-10 * tc[3]
            -2.11261600e-14 * tc[4];
        /*species 20: CH3OH */
        species[20] =
            +1.78970791e+00
            +1.40938292e-02 * tc[1]
            -6.36500835e-06 * tc[2]
            +1.38171085e-09 * tc[3]
            -1.17060220e-13 * tc[4];
        /*species 21: C2H */
        species[21] =
            +3.16780652e+00
            +4.75221902e-03 * tc[1]
            -1.83787077e-06 * tc[2]
            +3.04190252e-10 * tc[3]
            -1.77232770e-14 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            +4.14756964e+00
            +5.96166664e-03 * tc[1]
            -2.37294852e-06 * tc[2]
            +4.67412171e-10 * tc[3]
            -3.61235213e-14 * tc[4];
        /*species 23: C2H3 */
        species[23] =
            +3.01672400e+00
            +1.03302292e-02 * tc[1]
            -4.68082349e-06 * tc[2]
            +1.01763288e-09 * tc[3]
            -8.62607041e-14 * tc[4];
        /*species 24: C2H4 */
        species[24] =
            +2.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 25: C2H5 */
        species[25] =
            +1.95465642e+00
            +1.73972722e-02 * tc[1]
            -7.98206668e-06 * tc[2]
            +1.75217689e-09 * tc[3]
            -1.49641576e-13 * tc[4];
        /*species 26: C2H6 */
        species[26] =
            +1.07188150e+00
            +2.16852677e-02 * tc[1]
            -1.00256067e-05 * tc[2]
            +2.21412001e-09 * tc[3]
            -1.90002890e-13 * tc[4];
        /*species 27: HCCO */
        species[27] =
            +5.62820580e+00
            +4.08534010e-03 * tc[1]
            -1.59345470e-06 * tc[2]
            +2.86260520e-10 * tc[3]
            -1.94078320e-14 * tc[4];
        /*species 28: CH2CO */
        species[28] =
            +4.51129732e+00
            +9.00359745e-03 * tc[1]
            -4.16939635e-06 * tc[2]
            +9.23345882e-10 * tc[3]
            -7.94838201e-14 * tc[4];
        /*species 29: HCCOH */
        species[29] =
            +5.92382910e+00
            +6.79236000e-03 * tc[1]
            -2.56585640e-06 * tc[2]
            +4.49878410e-10 * tc[3]
            -2.99401010e-14 * tc[4];
        /*species 30: N2 */
        species[30] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 31: AR */
        species[31] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
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
        /*species 0: H2 */
        species[0] =
            +1.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 / tc[1];
        /*species 1: H */
        species[1] =
            +1.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 / tc[1];
        /*species 2: O */
        species[2] =
            +2.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 / tc[1];
        /*species 3: O2 */
        species[3] =
            +2.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 / tc[1];
        /*species 4: OH */
        species[4] =
            +2.99201543e+00
            -1.20065876e-03 * tc[1]
            +1.53931280e-06 * tc[2]
            -9.70283332e-10 * tc[3]
            +2.72822940e-13 * tc[4]
            +3.61508056e+03 / tc[1];
        /*species 5: H2O */
        species[5] =
            +3.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 / tc[1];
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 / tc[1];
        /*species 7: H2O2 */
        species[7] =
            +3.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 / tc[1];
        /*species 8: C */
        species[8] =
            +1.55423955e+00
            -1.60768862e-04 * tc[1]
            +2.44597415e-07 * tc[2]
            -1.83058722e-10 * tc[3]
            +5.33042892e-14 * tc[4]
            +8.54438832e+04 / tc[1];
        /*species 9: CH */
        species[9] =
            +2.48981665e+00
            +1.61917771e-04 * tc[1]
            -5.62996883e-07 * tc[2]
            +7.90543317e-10 * tc[3]
            -2.81218134e-13 * tc[4]
            +7.07972934e+04 / tc[1];
        /*species 10: CH2 */
        species[10] =
            +2.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 / tc[1];
        /*species 11: CH2(S) */
        species[11] =
            +3.19860411e+00
            -1.18330710e-03 * tc[1]
            +2.74432073e-06 * tc[2]
            -1.67203995e-09 * tc[3]
            +3.88629474e-13 * tc[4]
            +5.04968163e+04 / tc[1];
        /*species 12: CH3 */
        species[12] =
            +2.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 / tc[1];
        /*species 13: CH4 */
        species[13] =
            +4.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 / tc[1];
        /*species 14: CO */
        species[14] =
            +2.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 / tc[1];
        /*species 15: CO2 */
        species[15] =
            +1.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 / tc[1];
        /*species 16: HCO */
        species[16] =
            +3.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 / tc[1];
        /*species 17: CH2O */
        species[17] =
            +3.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 / tc[1];
        /*species 18: CH2OH */
        species[18] =
            +2.86388918e+00
            +2.79836152e-03 * tc[1]
            +1.97757264e-06 * tc[2]
            -2.61330030e-09 * tc[3]
            +8.73934556e-13 * tc[4]
            -3.19391367e+03 / tc[1];
        /*species 19: CH3O */
        species[19] =
            +1.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122000e-13 * tc[4]
            +9.78601100e+02 / tc[1];
        /*species 20: CH3OH */
        species[20] =
            +4.71539582e+00
            -7.61545645e-03 * tc[1]
            +2.17480385e-05 * tc[2]
            -1.77701722e-08 * tc[3]
            +5.22705396e-12 * tc[4]
            -2.56427656e+04 / tc[1];
        /*species 21: C2H */
        species[21] =
            +1.88965733e+00
            +6.70498055e-03 * tc[1]
            -9.49231670e-06 * tc[2]
            +7.36977613e-09 * tc[3]
            -2.18663022e-12 * tc[4]
            +6.68393932e+04 / tc[1];
        /*species 22: C2H2 */
        species[22] =
            -1.91318906e-01
            +1.16807815e-02 * tc[1]
            -1.18390605e-05 * tc[2]
            +7.00381092e-09 * tc[3]
            -1.70014595e-12 * tc[4]
            +2.64289807e+04 / tc[1];
        /*species 23: C2H3 */
        species[23] =
            +2.21246645e+00
            +7.57395810e-04 * tc[1]
            +8.64031373e-06 * tc[2]
            -8.94144617e-09 * tc[3]
            +2.94301746e-12 * tc[4]
            +3.48598468e+04 / tc[1];
        /*species 24: C2H4 */
        species[24] =
            +2.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 / tc[1];
        /*species 25: C2H5 */
        species[25] =
            +3.30646568e+00
            -2.09329446e-03 * tc[1]
            +1.65714269e-05 * tc[2]
            -1.49781651e-08 * tc[3]
            +4.61018008e-12 * tc[4]
            +1.28416265e+04 / tc[1];
        /*species 26: C2H6 */
        species[26] =
            +3.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 / tc[1];
        /*species 27: HCCO */
        species[27] =
            +1.25172140e+00
            +8.82751050e-03 * tc[1]
            -7.90970033e-06 * tc[2]
            +4.31893975e-09 * tc[3]
            -1.01329622e-12 * tc[4]
            +2.00594490e+04 / tc[1];
        /*species 28: CH2CO */
        species[28] =
            +1.13583630e+00
            +9.05943605e-03 * tc[1]
            -5.79824913e-06 * tc[2]
            +2.33599392e-09 * tc[3]
            -4.02915230e-13 * tc[4]
            -7.04291804e+03 / tc[1];
        /*species 29: HCCOH */
        species[29] =
            +2.42373300e-01
            +1.55361005e-02 * tc[1]
            -1.69556213e-05 * tc[2]
            +1.07842828e-08 * tc[3]
            -2.80291880e-12 * tc[4]
            +8.03161430e+03 / tc[1];
        /*species 30: N2 */
        species[30] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 / tc[1];
        /*species 31: AR */
        species[31] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 / tc[1];
    } else {
        /*species 0: H2 */
        species[0] =
            +2.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 / tc[1];
        /*species 1: H */
        species[1] =
            +1.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 / tc[1];
        /*species 2: O */
        species[2] =
            +1.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 / tc[1];
        /*species 3: O2 */
        species[3] =
            +2.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 / tc[1];
        /*species 4: OH */
        species[4] =
            +2.09288767e+00
            +2.74214858e-04 * tc[1]
            +4.21684093e-08 * tc[2]
            -2.19865389e-11 * tc[3]
            +2.34824752e-15 * tc[4]
            +3.85865700e+03 / tc[1];
        /*species 5: H2O */
        species[5] =
            +2.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 / tc[1];
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 / tc[1];
        /*species 7: H2O2 */
        species[7] =
            +3.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 / tc[1];
        /*species 8: C */
        species[8] =
            +1.49266888e+00
            +2.39944642e-05 * tc[1]
            -2.41445007e-08 * tc[2]
            +9.35727573e-12 * tc[3]
            -9.74555786e-16 * tc[4]
            +8.54512953e+04 / tc[1];
        /*species 9: CH */
        species[9] =
            +1.87846473e+00
            +4.85456840e-04 * tc[1]
            +4.81485517e-08 * tc[2]
            -3.26719623e-11 * tc[3]
            +3.52158766e-15 * tc[4]
            +7.10124364e+04 / tc[1];
        /*species 10: CH2 */
        species[10] =
            +1.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 / tc[1];
        /*species 11: CH2(S) */
        species[11] =
            +1.29203842e+00
            +2.32794318e-03 * tc[1]
            -6.70639823e-07 * tc[2]
            +1.04476500e-10 * tc[3]
            -6.79432730e-15 * tc[4]
            +5.09259997e+04 / tc[1];
        /*species 12: CH3 */
        species[12] =
            +1.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 / tc[1];
        /*species 13: CH4 */
        species[13] =
            -9.25148505e-01
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 / tc[1];
        /*species 14: CO */
        species[14] =
            +1.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 / tc[1];
        /*species 15: CO2 */
        species[15] =
            +2.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 / tc[1];
        /*species 16: HCO */
        species[16] =
            +1.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 / tc[1];
        /*species 17: CH2O */
        species[17] =
            +7.60690080e-01
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 / tc[1];
        /*species 18: CH2OH */
        species[18] =
            +2.69266569e+00
            +4.32288399e-03 * tc[1]
            -1.25033707e-06 * tc[2]
            +1.96808659e-10 * tc[3]
            -1.29710840e-14 * tc[4]
            -3.24250627e+03 / tc[1];
        /*species 19: CH3O */
        species[19] =
            +2.77079900e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832520e+02 / tc[1];
        /*species 20: CH3OH */
        species[20] =
            +7.89707910e-01
            +7.04691460e-03 * tc[1]
            -2.12166945e-06 * tc[2]
            +3.45427713e-10 * tc[3]
            -2.34120440e-14 * tc[4]
            -2.53748747e+04 / tc[1];
        /*species 21: C2H */
        species[21] =
            +2.16780652e+00
            +2.37610951e-03 * tc[1]
            -6.12623590e-07 * tc[2]
            +7.60475630e-11 * tc[3]
            -3.54465540e-15 * tc[4]
            +6.71210650e+04 / tc[1];
        /*species 22: C2H2 */
        species[22] =
            +3.14756964e+00
            +2.98083332e-03 * tc[1]
            -7.90982840e-07 * tc[2]
            +1.16853043e-10 * tc[3]
            -7.22470426e-15 * tc[4]
            +2.59359992e+04 / tc[1];
        /*species 23: C2H3 */
        species[23] =
            +2.01672400e+00
            +5.16511460e-03 * tc[1]
            -1.56027450e-06 * tc[2]
            +2.54408220e-10 * tc[3]
            -1.72521408e-14 * tc[4]
            +3.46128739e+04 / tc[1];
        /*species 24: C2H4 */
        species[24] =
            +1.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 / tc[1];
        /*species 25: C2H5 */
        species[25] =
            +9.54656420e-01
            +8.69863610e-03 * tc[1]
            -2.66068889e-06 * tc[2]
            +4.38044223e-10 * tc[3]
            -2.99283152e-14 * tc[4]
            +1.28575200e+04 / tc[1];
        /*species 26: C2H6 */
        species[26] =
            +7.18815000e-02
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 / tc[1];
        /*species 27: HCCO */
        species[27] =
            +4.62820580e+00
            +2.04267005e-03 * tc[1]
            -5.31151567e-07 * tc[2]
            +7.15651300e-11 * tc[3]
            -3.88156640e-15 * tc[4]
            +1.93272150e+04 / tc[1];
        /*species 28: CH2CO */
        species[28] =
            +3.51129732e+00
            +4.50179872e-03 * tc[1]
            -1.38979878e-06 * tc[2]
            +2.30836470e-10 * tc[3]
            -1.58967640e-14 * tc[4]
            -7.55105311e+03 / tc[1];
        /*species 29: HCCOH */
        species[29] =
            +4.92382910e+00
            +3.39618000e-03 * tc[1]
            -8.55285467e-07 * tc[2]
            +1.12469603e-10 * tc[3]
            -5.98802020e-15 * tc[4]
            +7.26462600e+03 / tc[1];
        /*species 30: N2 */
        species[30] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 / tc[1];
        /*species 31: AR */
        species[31] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 / tc[1];
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
        /*species 0: H2 */
        species[0] =
            +2.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 / tc[1];
        /*species 1: H */
        species[1] =
            +2.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 / tc[1];
        /*species 2: O */
        species[2] =
            +3.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 / tc[1];
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 / tc[1];
        /*species 4: OH */
        species[4] =
            +3.99201543e+00
            -1.20065876e-03 * tc[1]
            +1.53931280e-06 * tc[2]
            -9.70283332e-10 * tc[3]
            +2.72822940e-13 * tc[4]
            +3.61508056e+03 / tc[1];
        /*species 5: H2O */
        species[5] =
            +4.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 / tc[1];
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 / tc[1];
        /*species 7: H2O2 */
        species[7] =
            +4.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 / tc[1];
        /*species 8: C */
        species[8] =
            +2.55423955e+00
            -1.60768862e-04 * tc[1]
            +2.44597415e-07 * tc[2]
            -1.83058722e-10 * tc[3]
            +5.33042892e-14 * tc[4]
            +8.54438832e+04 / tc[1];
        /*species 9: CH */
        species[9] =
            +3.48981665e+00
            +1.61917771e-04 * tc[1]
            -5.62996883e-07 * tc[2]
            +7.90543317e-10 * tc[3]
            -2.81218134e-13 * tc[4]
            +7.07972934e+04 / tc[1];
        /*species 10: CH2 */
        species[10] =
            +3.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 / tc[1];
        /*species 11: CH2(S) */
        species[11] =
            +4.19860411e+00
            -1.18330710e-03 * tc[1]
            +2.74432073e-06 * tc[2]
            -1.67203995e-09 * tc[3]
            +3.88629474e-13 * tc[4]
            +5.04968163e+04 / tc[1];
        /*species 12: CH3 */
        species[12] =
            +3.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 / tc[1];
        /*species 13: CH4 */
        species[13] =
            +5.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 / tc[1];
        /*species 14: CO */
        species[14] =
            +3.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 / tc[1];
        /*species 15: CO2 */
        species[15] =
            +2.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 / tc[1];
        /*species 16: HCO */
        species[16] =
            +4.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 / tc[1];
        /*species 17: CH2O */
        species[17] =
            +4.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 / tc[1];
        /*species 18: CH2OH */
        species[18] =
            +3.86388918e+00
            +2.79836152e-03 * tc[1]
            +1.97757264e-06 * tc[2]
            -2.61330030e-09 * tc[3]
            +8.73934556e-13 * tc[4]
            -3.19391367e+03 / tc[1];
        /*species 19: CH3O */
        species[19] =
            +2.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122000e-13 * tc[4]
            +9.78601100e+02 / tc[1];
        /*species 20: CH3OH */
        species[20] =
            +5.71539582e+00
            -7.61545645e-03 * tc[1]
            +2.17480385e-05 * tc[2]
            -1.77701722e-08 * tc[3]
            +5.22705396e-12 * tc[4]
            -2.56427656e+04 / tc[1];
        /*species 21: C2H */
        species[21] =
            +2.88965733e+00
            +6.70498055e-03 * tc[1]
            -9.49231670e-06 * tc[2]
            +7.36977613e-09 * tc[3]
            -2.18663022e-12 * tc[4]
            +6.68393932e+04 / tc[1];
        /*species 22: C2H2 */
        species[22] =
            +8.08681094e-01
            +1.16807815e-02 * tc[1]
            -1.18390605e-05 * tc[2]
            +7.00381092e-09 * tc[3]
            -1.70014595e-12 * tc[4]
            +2.64289807e+04 / tc[1];
        /*species 23: C2H3 */
        species[23] =
            +3.21246645e+00
            +7.57395810e-04 * tc[1]
            +8.64031373e-06 * tc[2]
            -8.94144617e-09 * tc[3]
            +2.94301746e-12 * tc[4]
            +3.48598468e+04 / tc[1];
        /*species 24: C2H4 */
        species[24] =
            +3.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 / tc[1];
        /*species 25: C2H5 */
        species[25] =
            +4.30646568e+00
            -2.09329446e-03 * tc[1]
            +1.65714269e-05 * tc[2]
            -1.49781651e-08 * tc[3]
            +4.61018008e-12 * tc[4]
            +1.28416265e+04 / tc[1];
        /*species 26: C2H6 */
        species[26] =
            +4.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 / tc[1];
        /*species 27: HCCO */
        species[27] =
            +2.25172140e+00
            +8.82751050e-03 * tc[1]
            -7.90970033e-06 * tc[2]
            +4.31893975e-09 * tc[3]
            -1.01329622e-12 * tc[4]
            +2.00594490e+04 / tc[1];
        /*species 28: CH2CO */
        species[28] =
            +2.13583630e+00
            +9.05943605e-03 * tc[1]
            -5.79824913e-06 * tc[2]
            +2.33599392e-09 * tc[3]
            -4.02915230e-13 * tc[4]
            -7.04291804e+03 / tc[1];
        /*species 29: HCCOH */
        species[29] =
            +1.24237330e+00
            +1.55361005e-02 * tc[1]
            -1.69556213e-05 * tc[2]
            +1.07842828e-08 * tc[3]
            -2.80291880e-12 * tc[4]
            +8.03161430e+03 / tc[1];
        /*species 30: N2 */
        species[30] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 / tc[1];
        /*species 31: AR */
        species[31] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 / tc[1];
    } else {
        /*species 0: H2 */
        species[0] =
            +3.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 / tc[1];
        /*species 1: H */
        species[1] =
            +2.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 / tc[1];
        /*species 2: O */
        species[2] =
            +2.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 / tc[1];
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 / tc[1];
        /*species 4: OH */
        species[4] =
            +3.09288767e+00
            +2.74214858e-04 * tc[1]
            +4.21684093e-08 * tc[2]
            -2.19865389e-11 * tc[3]
            +2.34824752e-15 * tc[4]
            +3.85865700e+03 / tc[1];
        /*species 5: H2O */
        species[5] =
            +3.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 / tc[1];
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 / tc[1];
        /*species 7: H2O2 */
        species[7] =
            +4.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 / tc[1];
        /*species 8: C */
        species[8] =
            +2.49266888e+00
            +2.39944642e-05 * tc[1]
            -2.41445007e-08 * tc[2]
            +9.35727573e-12 * tc[3]
            -9.74555786e-16 * tc[4]
            +8.54512953e+04 / tc[1];
        /*species 9: CH */
        species[9] =
            +2.87846473e+00
            +4.85456840e-04 * tc[1]
            +4.81485517e-08 * tc[2]
            -3.26719623e-11 * tc[3]
            +3.52158766e-15 * tc[4]
            +7.10124364e+04 / tc[1];
        /*species 10: CH2 */
        species[10] =
            +2.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 / tc[1];
        /*species 11: CH2(S) */
        species[11] =
            +2.29203842e+00
            +2.32794318e-03 * tc[1]
            -6.70639823e-07 * tc[2]
            +1.04476500e-10 * tc[3]
            -6.79432730e-15 * tc[4]
            +5.09259997e+04 / tc[1];
        /*species 12: CH3 */
        species[12] =
            +2.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 / tc[1];
        /*species 13: CH4 */
        species[13] =
            +7.48514950e-02
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 / tc[1];
        /*species 14: CO */
        species[14] =
            +2.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 / tc[1];
        /*species 15: CO2 */
        species[15] =
            +3.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 / tc[1];
        /*species 16: HCO */
        species[16] =
            +2.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 / tc[1];
        /*species 17: CH2O */
        species[17] =
            +1.76069008e+00
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 / tc[1];
        /*species 18: CH2OH */
        species[18] =
            +3.69266569e+00
            +4.32288399e-03 * tc[1]
            -1.25033707e-06 * tc[2]
            +1.96808659e-10 * tc[3]
            -1.29710840e-14 * tc[4]
            -3.24250627e+03 / tc[1];
        /*species 19: CH3O */
        species[19] =
            +3.77079900e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832520e+02 / tc[1];
        /*species 20: CH3OH */
        species[20] =
            +1.78970791e+00
            +7.04691460e-03 * tc[1]
            -2.12166945e-06 * tc[2]
            +3.45427713e-10 * tc[3]
            -2.34120440e-14 * tc[4]
            -2.53748747e+04 / tc[1];
        /*species 21: C2H */
        species[21] =
            +3.16780652e+00
            +2.37610951e-03 * tc[1]
            -6.12623590e-07 * tc[2]
            +7.60475630e-11 * tc[3]
            -3.54465540e-15 * tc[4]
            +6.71210650e+04 / tc[1];
        /*species 22: C2H2 */
        species[22] =
            +4.14756964e+00
            +2.98083332e-03 * tc[1]
            -7.90982840e-07 * tc[2]
            +1.16853043e-10 * tc[3]
            -7.22470426e-15 * tc[4]
            +2.59359992e+04 / tc[1];
        /*species 23: C2H3 */
        species[23] =
            +3.01672400e+00
            +5.16511460e-03 * tc[1]
            -1.56027450e-06 * tc[2]
            +2.54408220e-10 * tc[3]
            -1.72521408e-14 * tc[4]
            +3.46128739e+04 / tc[1];
        /*species 24: C2H4 */
        species[24] =
            +2.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 / tc[1];
        /*species 25: C2H5 */
        species[25] =
            +1.95465642e+00
            +8.69863610e-03 * tc[1]
            -2.66068889e-06 * tc[2]
            +4.38044223e-10 * tc[3]
            -2.99283152e-14 * tc[4]
            +1.28575200e+04 / tc[1];
        /*species 26: C2H6 */
        species[26] =
            +1.07188150e+00
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 / tc[1];
        /*species 27: HCCO */
        species[27] =
            +5.62820580e+00
            +2.04267005e-03 * tc[1]
            -5.31151567e-07 * tc[2]
            +7.15651300e-11 * tc[3]
            -3.88156640e-15 * tc[4]
            +1.93272150e+04 / tc[1];
        /*species 28: CH2CO */
        species[28] =
            +4.51129732e+00
            +4.50179872e-03 * tc[1]
            -1.38979878e-06 * tc[2]
            +2.30836470e-10 * tc[3]
            -1.58967640e-14 * tc[4]
            -7.55105311e+03 / tc[1];
        /*species 29: HCCOH */
        species[29] =
            +5.92382910e+00
            +3.39618000e-03 * tc[1]
            -8.55285467e-07 * tc[2]
            +1.12469603e-10 * tc[3]
            -5.98802020e-15 * tc[4]
            +7.26462600e+03 / tc[1];
        /*species 30: N2 */
        species[30] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 / tc[1];
        /*species 31: AR */
        species[31] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 / tc[1];
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
        /*species 0: H2 */
        species[0] =
            +2.34433112e+00 * tc[0]
            +7.98052075e-03 * tc[1]
            -9.73907550e-06 * tc[2]
            +6.71906980e-09 * tc[3]
            -1.84402940e-12 * tc[4]
            +6.83010238e-01 ;
        /*species 1: H */
        species[1] =
            +2.50000000e+00 * tc[0]
            +7.05332819e-13 * tc[1]
            -9.97959820e-16 * tc[2]
            +7.66938773e-19 * tc[3]
            -2.31933083e-22 * tc[4]
            -4.46682853e-01 ;
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
            -2.99673416e-03 * tc[1]
            +4.92365101e-06 * tc[2]
            -3.22709836e-09 * tc[3]
            +8.10932092e-13 * tc[4]
            +3.65767573e+00 ;
        /*species 4: OH */
        species[4] =
            +3.99201543e+00 * tc[0]
            -2.40131752e-03 * tc[1]
            +2.30896920e-06 * tc[2]
            -1.29371111e-09 * tc[3]
            +3.41028675e-13 * tc[4]
            -1.03925458e-01 ;
        /*species 5: H2O */
        species[5] =
            +4.19864056e+00 * tc[0]
            -2.03643410e-03 * tc[1]
            +3.26020105e-06 * tc[2]
            -1.82932354e-09 * tc[3]
            +4.42994543e-13 * tc[4]
            -8.49032208e-01 ;
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00 * tc[0]
            -4.74912051e-03 * tc[1]
            +1.05791445e-05 * tc[2]
            -8.09212980e-09 * tc[3]
            +2.32306281e-12 * tc[4]
            +3.71666245e+00 ;
        /*species 7: H2O2 */
        species[7] =
            +4.27611269e+00 * tc[0]
            -5.42822417e-04 * tc[1]
            +8.36678505e-06 * tc[2]
            -7.19236043e-09 * tc[3]
            +2.15613591e-12 * tc[4]
            +3.43505074e+00 ;
        /*species 8: C */
        species[8] =
            +2.55423955e+00 * tc[0]
            -3.21537724e-04 * tc[1]
            +3.66896122e-07 * tc[2]
            -2.44078296e-10 * tc[3]
            +6.66303615e-14 * tc[4]
            +4.53130848e+00 ;
        /*species 9: CH */
        species[9] =
            +3.48981665e+00 * tc[0]
            +3.23835541e-04 * tc[1]
            -8.44495325e-07 * tc[2]
            +1.05405776e-09 * tc[3]
            -3.51522668e-13 * tc[4]
            +2.08401108e+00 ;
        /*species 10: CH2 */
        species[10] =
            +3.76267867e+00 * tc[0]
            +9.68872143e-04 * tc[1]
            +1.39744921e-06 * tc[2]
            -1.28363718e-09 * tc[3]
            +4.21854298e-13 * tc[4]
            +1.56253185e+00 ;
        /*species 11: CH2(S) */
        species[11] =
            +4.19860411e+00 * tc[0]
            -2.36661419e-03 * tc[1]
            +4.11648110e-06 * tc[2]
            -2.22938660e-09 * tc[3]
            +4.85786843e-13 * tc[4]
            -7.69118967e-01 ;
        /*species 12: CH3 */
        species[12] =
            +3.67359040e+00 * tc[0]
            +2.01095175e-03 * tc[1]
            +2.86510928e-06 * tc[2]
            -2.29039142e-09 * tc[3]
            +6.35964335e-13 * tc[4]
            +1.60456433e+00 ;
        /*species 13: CH4 */
        species[13] =
            +5.14987613e+00 * tc[0]
            -1.36709788e-02 * tc[1]
            +2.45900299e-05 * tc[2]
            -1.61581009e-08 * tc[3]
            +4.16734890e-12 * tc[4]
            -4.64130376e+00 ;
        /*species 14: CO */
        species[14] =
            +3.57953347e+00 * tc[0]
            -6.10353680e-04 * tc[1]
            +5.08407165e-07 * tc[2]
            +3.02335295e-10 * tc[3]
            -2.26106125e-13 * tc[4]
            +3.50840928e+00 ;
        /*species 15: CO2 */
        species[15] =
            +2.35677352e+00 * tc[0]
            +8.98459677e-03 * tc[1]
            -3.56178134e-06 * tc[2]
            +8.19730073e-10 * tc[3]
            -3.59248870e-14 * tc[4]
            +9.90105222e+00 ;
        /*species 16: HCO */
        species[16] =
            +4.22118584e+00 * tc[0]
            -3.24392532e-03 * tc[1]
            +6.88997230e-06 * tc[2]
            -4.43813643e-09 * tc[3]
            +1.08442216e-12 * tc[4]
            +3.39437243e+00 ;
        /*species 17: CH2O */
        species[17] =
            +4.79372315e+00 * tc[0]
            -9.90833369e-03 * tc[1]
            +1.86610004e-05 * tc[2]
            -1.26428420e-08 * tc[3]
            +3.29431630e-12 * tc[4]
            +6.02812900e-01 ;
        /*species 18: CH2OH */
        species[18] =
            +3.86388918e+00 * tc[0]
            +5.59672304e-03 * tc[1]
            +2.96635895e-06 * tc[2]
            -3.48440040e-09 * tc[3]
            +1.09241820e-12 * tc[4]
            +5.47302243e+00 ;
        /*species 19: CH3O */
        species[19] =
            +2.10620400e+00 * tc[0]
            +7.21659500e-03 * tc[1]
            +2.66923600e-06 * tc[2]
            -2.45921200e-09 * tc[3]
            +5.18902500e-13 * tc[4]
            +1.31521770e+01 ;
        /*species 20: CH3OH */
        species[20] =
            +5.71539582e+00 * tc[0]
            -1.52309129e-02 * tc[1]
            +3.26220578e-05 * tc[2]
            -2.36935630e-08 * tc[3]
            +6.53381745e-12 * tc[4]
            -1.50409823e+00 ;
        /*species 21: C2H */
        species[21] =
            +2.88965733e+00 * tc[0]
            +1.34099611e-02 * tc[1]
            -1.42384751e-05 * tc[2]
            +9.82636817e-09 * tc[3]
            -2.73328777e-12 * tc[4]
            +6.22296438e+00 ;
        /*species 22: C2H2 */
        species[22] =
            +8.08681094e-01 * tc[0]
            +2.33615629e-02 * tc[1]
            -1.77585907e-05 * tc[2]
            +9.33841457e-09 * tc[3]
            -2.12518243e-12 * tc[4]
            +1.39397051e+01 ;
        /*species 23: C2H3 */
        species[23] =
            +3.21246645e+00 * tc[0]
            +1.51479162e-03 * tc[1]
            +1.29604706e-05 * tc[2]
            -1.19219282e-08 * tc[3]
            +3.67877182e-12 * tc[4]
            +8.51054025e+00 ;
        /*species 24: C2H4 */
        species[24] =
            +3.95920148e+00 * tc[0]
            -7.57052247e-03 * tc[1]
            +2.85495146e-05 * tc[2]
            -2.30529584e-08 * tc[3]
            +6.74710933e-12 * tc[4]
            +4.09733096e+00 ;
        /*species 25: C2H5 */
        species[25] =
            +4.30646568e+00 * tc[0]
            -4.18658892e-03 * tc[1]
            +2.48571403e-05 * tc[2]
            -1.99708869e-08 * tc[3]
            +5.76272510e-12 * tc[4]
            +4.70720924e+00 ;
        /*species 26: C2H6 */
        species[26] =
            +4.29142492e+00 * tc[0]
            -5.50154270e-03 * tc[1]
            +2.99719144e-05 * tc[2]
            -2.36155428e-08 * tc[3]
            +6.71714427e-12 * tc[4]
            +2.66682316e+00 ;
        /*species 27: HCCO */
        species[27] =
            +2.25172140e+00 * tc[0]
            +1.76550210e-02 * tc[1]
            -1.18645505e-05 * tc[2]
            +5.75858633e-09 * tc[3]
            -1.26662028e-12 * tc[4]
            +1.24904170e+01 ;
        /*species 28: CH2CO */
        species[28] =
            +2.13583630e+00 * tc[0]
            +1.81188721e-02 * tc[1]
            -8.69737370e-06 * tc[2]
            +3.11465856e-09 * tc[3]
            -5.03644037e-13 * tc[4]
            +1.22156480e+01 ;
        /*species 29: HCCOH */
        species[29] =
            +1.24237330e+00 * tc[0]
            +3.10722010e-02 * tc[1]
            -2.54334320e-05 * tc[2]
            +1.43790437e-08 * tc[3]
            -3.50364850e-12 * tc[4]
            +1.38743190e+01 ;
        /*species 30: N2 */
        species[30] =
            +3.29867700e+00 * tc[0]
            +1.40824040e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213500e-13 * tc[4]
            +3.95037200e+00 ;
        /*species 31: AR */
        species[31] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.36600000e+00 ;
    } else {
        /*species 0: H2 */
        species[0] =
            +3.33727920e+00 * tc[0]
            -4.94024731e-05 * tc[1]
            +2.49728389e-07 * tc[2]
            -5.98554647e-11 * tc[3]
            +5.00638440e-15 * tc[4]
            -3.20502331e+00 ;
        /*species 1: H */
        species[1] =
            +2.50000001e+00 * tc[0]
            -2.30842973e-11 * tc[1]
            +8.07809740e-15 * tc[2]
            -1.57838412e-18 * tc[3]
            +1.24549339e-22 * tc[4]
            -4.46682914e-01 ;
        /*species 2: O */
        species[2] =
            +2.56942078e+00 * tc[0]
            -8.59741137e-05 * tc[1]
            +2.09742295e-08 * tc[2]
            -3.33925997e-12 * tc[3]
            +3.07084227e-16 * tc[4]
            +4.78433864e+00 ;
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00 * tc[0]
            +1.48308754e-03 * tc[1]
            -3.78983334e-07 * tc[2]
            +6.98235183e-11 * tc[3]
            -5.41794485e-15 * tc[4]
            +5.45323129e+00 ;
        /*species 4: OH */
        species[4] =
            +3.09288767e+00 * tc[0]
            +5.48429716e-04 * tc[1]
            +6.32526140e-08 * tc[2]
            -2.93153852e-11 * tc[3]
            +2.93530940e-15 * tc[4]
            +4.47669610e+00 ;
        /*species 5: H2O */
        species[5] =
            +3.03399249e+00 * tc[0]
            +2.17691804e-03 * tc[1]
            -8.20362590e-08 * tc[2]
            -3.23473290e-11 * tc[3]
            +4.20502480e-15 * tc[4]
            +4.96677010e+00 ;
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00 * tc[0]
            +2.23982013e-03 * tc[1]
            -3.16829075e-07 * tc[2]
            +3.80821233e-11 * tc[3]
            -2.69771337e-15 * tc[4]
            +3.78510215e+00 ;
        /*species 7: H2O2 */
        species[7] =
            +4.16500285e+00 * tc[0]
            +4.90831694e-03 * tc[1]
            -9.50696125e-07 * tc[2]
            +1.23728662e-10 * tc[3]
            -7.19770763e-15 * tc[4]
            +2.91615662e+00 ;
        /*species 8: C */
        species[8] =
            +2.49266888e+00 * tc[0]
            +4.79889284e-05 * tc[1]
            -3.62167510e-08 * tc[2]
            +1.24763676e-11 * tc[3]
            -1.21819473e-15 * tc[4]
            +4.80150373e+00 ;
        /*species 9: CH */
        species[9] =
            +2.87846473e+00 * tc[0]
            +9.70913681e-04 * tc[1]
            +7.22228275e-08 * tc[2]
            -4.35626163e-11 * tc[3]
            +4.40198457e-15 * tc[4]
            +5.48497999e+00 ;
        /*species 10: CH2 */
        species[10] =
            +2.87410113e+00 * tc[0]
            +3.65639292e-03 * tc[1]
            -7.04472985e-07 * tc[2]
            +8.67265163e-11 * tc[3]
            -4.69318918e-15 * tc[4]
            +6.17119324e+00 ;
        /*species 11: CH2(S) */
        species[11] =
            +2.29203842e+00 * tc[0]
            +4.65588637e-03 * tc[1]
            -1.00595973e-06 * tc[2]
            +1.39302000e-10 * tc[3]
            -8.49290912e-15 * tc[4]
            +8.62650169e+00 ;
        /*species 12: CH3 */
        species[12] =
            +2.28571772e+00 * tc[0]
            +7.23990037e-03 * tc[1]
            -1.49357174e-06 * tc[2]
            +1.98561548e-10 * tc[3]
            -1.16788599e-14 * tc[4]
            +8.48007179e+00 ;
        /*species 13: CH4 */
        species[13] =
            +7.48514950e-02 * tc[0]
            +1.33909467e-02 * tc[1]
            -2.86642905e-06 * tc[2]
            +4.07641783e-10 * tc[3]
            -2.54538075e-14 * tc[4]
            +1.84373180e+01 ;
        /*species 14: CO */
        species[14] =
            +2.71518561e+00 * tc[0]
            +2.06252743e-03 * tc[1]
            -4.99412886e-07 * tc[2]
            +7.66843360e-11 * tc[3]
            -5.09119290e-15 * tc[4]
            +7.81868772e+00 ;
        /*species 15: CO2 */
        species[15] =
            +3.85746029e+00 * tc[0]
            +4.41437026e-03 * tc[1]
            -1.10740702e-06 * tc[2]
            +1.74496729e-10 * tc[3]
            -1.18021041e-14 * tc[4]
            +2.27163806e+00 ;
        /*species 16: HCO */
        species[16] =
            +2.77217438e+00 * tc[0]
            +4.95695526e-03 * tc[1]
            -1.24222806e-06 * tc[2]
            +1.96387259e-10 * tc[3]
            -1.33377178e-14 * tc[4]
            +9.79834492e+00 ;
        /*species 17: CH2O */
        species[17] =
            +1.76069008e+00 * tc[0]
            +9.20000082e-03 * tc[1]
            -2.21129406e-06 * tc[2]
            +3.35470707e-10 * tc[3]
            -2.20963910e-14 * tc[4]
            +1.36563230e+01 ;
        /*species 18: CH2OH */
        species[18] =
            +3.69266569e+00 * tc[0]
            +8.64576797e-03 * tc[1]
            -1.87550560e-06 * tc[2]
            +2.62411545e-10 * tc[3]
            -1.62138550e-14 * tc[4]
            +5.81043215e+00 ;
        /*species 19: CH3O */
        species[19] =
            +3.77079900e+00 * tc[0]
            +7.87149700e-03 * tc[1]
            -1.32819200e-06 * tc[2]
            +1.31481033e-10 * tc[3]
            -5.28154000e-15 * tc[4]
            +2.92957500e+00 ;
        /*species 20: CH3OH */
        species[20] =
            +1.78970791e+00 * tc[0]
            +1.40938292e-02 * tc[1]
            -3.18250418e-06 * tc[2]
            +4.60570283e-10 * tc[3]
            -2.92650550e-14 * tc[4]
            +1.45023623e+01 ;
        /*species 21: C2H */
        species[21] =
            +3.16780652e+00 * tc[0]
            +4.75221902e-03 * tc[1]
            -9.18935385e-07 * tc[2]
            +1.01396751e-10 * tc[3]
            -4.43081925e-15 * tc[4]
            +6.63589475e+00 ;
        /*species 22: C2H2 */
        species[22] =
            +4.14756964e+00 * tc[0]
            +5.96166664e-03 * tc[1]
            -1.18647426e-06 * tc[2]
            +1.55804057e-10 * tc[3]
            -9.03088033e-15 * tc[4]
            -1.23028121e+00 ;
        /*species 23: C2H3 */
        species[23] =
            +3.01672400e+00 * tc[0]
            +1.03302292e-02 * tc[1]
            -2.34041174e-06 * tc[2]
            +3.39210960e-10 * tc[3]
            -2.15651760e-14 * tc[4]
            +7.78732378e+00 ;
        /*species 24: C2H4 */
        species[24] =
            +2.03611116e+00 * tc[0]
            +1.46454151e-02 * tc[1]
            -3.35538958e-06 * tc[2]
            +4.90743077e-10 * tc[3]
            -3.14265152e-14 * tc[4]
            +1.03053693e+01 ;
        /*species 25: C2H5 */
        species[25] =
            +1.95465642e+00 * tc[0]
            +1.73972722e-02 * tc[1]
            -3.99103334e-06 * tc[2]
            +5.84058963e-10 * tc[3]
            -3.74103940e-14 * tc[4]
            +1.34624343e+01 ;
        /*species 26: C2H6 */
        species[26] =
            +1.07188150e+00 * tc[0]
            +2.16852677e-02 * tc[1]
            -5.01280335e-06 * tc[2]
            +7.38040003e-10 * tc[3]
            -4.75007225e-14 * tc[4]
            +1.51156107e+01 ;
        /*species 27: HCCO */
        species[27] =
            +5.62820580e+00 * tc[0]
            +4.08534010e-03 * tc[1]
            -7.96727350e-07 * tc[2]
            +9.54201733e-11 * tc[3]
            -4.85195800e-15 * tc[4]
            -3.93025950e+00 ;
        /*species 28: CH2CO */
        species[28] =
            +4.51129732e+00 * tc[0]
            +9.00359745e-03 * tc[1]
            -2.08469817e-06 * tc[2]
            +3.07781961e-10 * tc[3]
            -1.98709550e-14 * tc[4]
            +6.32247205e-01 ;
        /*species 29: HCCOH */
        species[29] =
            +5.92382910e+00 * tc[0]
            +6.79236000e-03 * tc[1]
            -1.28292820e-06 * tc[2]
            +1.49959470e-10 * tc[3]
            -7.48502525e-15 * tc[4]
            -7.60177420e+00 ;
        /*species 30: N2 */
        species[30] =
            +2.92664000e+00 * tc[0]
            +1.48797680e-03 * tc[1]
            -2.84238000e-07 * tc[2]
            +3.36567933e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
        /*species 31: AR */
        species[31] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.36600000e+00 ;
    }
    return;
}


/*save molecular weights into array */
void molecularWeight(double * wt)
{
    wt[0] = 2.015940; /*H2 */
    wt[1] = 1.007970; /*H */
    wt[2] = 15.999400; /*O */
    wt[3] = 31.998800; /*O2 */
    wt[4] = 17.007370; /*OH */
    wt[5] = 18.015340; /*H2O */
    wt[6] = 33.006770; /*HO2 */
    wt[7] = 34.014740; /*H2O2 */
    wt[8] = 12.011000; /*C */
    wt[9] = 13.018970; /*CH */
    wt[10] = 14.026940; /*CH2 */
    wt[11] = 14.026940; /*CH2(S) */
    wt[12] = 15.034910; /*CH3 */
    wt[13] = 16.042880; /*CH4 */
    wt[14] = 28.010400; /*CO */
    wt[15] = 44.009800; /*CO2 */
    wt[16] = 29.018370; /*HCO */
    wt[17] = 30.026340; /*CH2O */
    wt[18] = 31.034310; /*CH2OH */
    wt[19] = 31.034310; /*CH3O */
    wt[20] = 32.042280; /*CH3OH */
    wt[21] = 25.029970; /*C2H */
    wt[22] = 26.037940; /*C2H2 */
    wt[23] = 27.045910; /*C2H3 */
    wt[24] = 28.053880; /*C2H4 */
    wt[25] = 29.061850; /*C2H5 */
    wt[26] = 30.069820; /*C2H6 */
    wt[27] = 41.029370; /*HCCO */
    wt[28] = 42.037340; /*CH2CO */
    wt[29] = 42.037340; /*HCCOH */
    wt[30] = 28.013400; /*N2 */
    wt[31] = 39.948000; /*AR */

    return;
}


/*get temperature given internal energy in mass units and mass fracs */
int feeytt_(double * e, double * y, int * iwrk, double * rwrk, double * t)
{
    const int maxiter = 50;
    const double tol  = 0.001;
    double ein  = *e;
    double tmin = 300; // max lower bound for thermo def
    double tmax = 3000; // min upper bound for thermo def
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
    y[0] = phi[0]*2.015940;   XW += y[0]; /*H2 */
    y[1] = phi[1]*1.007970;   XW += y[1]; /*H */
    y[2] = phi[2]*15.999400;   XW += y[2]; /*O */
    y[3] = phi[3]*31.998800;   XW += y[3]; /*O2 */
    y[4] = phi[4]*17.007370;   XW += y[4]; /*OH */
    y[5] = phi[5]*18.015340;   XW += y[5]; /*H2O */
    y[6] = phi[6]*33.006770;   XW += y[6]; /*HO2 */
    y[7] = phi[7]*34.014740;   XW += y[7]; /*H2O2 */
    y[8] = phi[8]*12.011000;   XW += y[8]; /*C */
    y[9] = phi[9]*13.018970;   XW += y[9]; /*CH */
    y[10] = phi[10]*14.026940;   XW += y[10]; /*CH2 */
    y[11] = phi[11]*14.026940;   XW += y[11]; /*CH2(S) */
    y[12] = phi[12]*15.034910;   XW += y[12]; /*CH3 */
    y[13] = phi[13]*16.042880;   XW += y[13]; /*CH4 */
    y[14] = phi[14]*28.010400;   XW += y[14]; /*CO */
    y[15] = phi[15]*44.009800;   XW += y[15]; /*CO2 */
    y[16] = phi[16]*29.018370;   XW += y[16]; /*HCO */
    y[17] = phi[17]*30.026340;   XW += y[17]; /*CH2O */
    y[18] = phi[18]*31.034310;   XW += y[18]; /*CH2OH */
    y[19] = phi[19]*31.034310;   XW += y[19]; /*CH3O */
    y[20] = phi[20]*32.042280;   XW += y[20]; /*CH3OH */
    y[21] = phi[21]*25.029970;   XW += y[21]; /*C2H */
    y[22] = phi[22]*26.037940;   XW += y[22]; /*C2H2 */
    y[23] = phi[23]*27.045910;   XW += y[23]; /*C2H3 */
    y[24] = phi[24]*28.053880;   XW += y[24]; /*C2H4 */
    y[25] = phi[25]*29.061850;   XW += y[25]; /*C2H5 */
    y[26] = phi[26]*30.069820;   XW += y[26]; /*C2H6 */
    y[27] = phi[27]*41.029370;   XW += y[27]; /*HCCO */
    y[28] = phi[28]*42.037340;   XW += y[28]; /*CH2CO */
    y[29] = phi[29]*42.037340;   XW += y[29]; /*HCCOH */
    y[30] = phi[30]*28.013400;   XW += y[30]; /*N2 */
    y[31] = phi[31]*39.948000;   XW += y[31]; /*AR */
    for (id = 0; id < 32; ++id) {
        y[id] = y[id]/XW;
    }

    return;
}


/*convert y[species] (mass fracs) to phi[species] (specific mole num) */
void feytphi_(double * y, int * iwrk, double * rwrk, double * phi)
{
    phi[0] = y[0]/ 2.01594000e-03; /*H2 (wt in kg) */
    phi[1] = y[1]/ 1.00797000e-03; /*H (wt in kg) */
    phi[2] = y[2]/ 1.59994000e-02; /*O (wt in kg) */
    phi[3] = y[3]/ 3.19988000e-02; /*O2 (wt in kg) */
    phi[4] = y[4]/ 1.70073700e-02; /*OH (wt in kg) */
    phi[5] = y[5]/ 1.80153400e-02; /*H2O (wt in kg) */
    phi[6] = y[6]/ 3.30067700e-02; /*HO2 (wt in kg) */
    phi[7] = y[7]/ 3.40147400e-02; /*H2O2 (wt in kg) */
    phi[8] = y[8]/ 1.20110000e-02; /*C (wt in kg) */
    phi[9] = y[9]/ 1.30189700e-02; /*CH (wt in kg) */
    phi[10] = y[10]/ 1.40269400e-02; /*CH2 (wt in kg) */
    phi[11] = y[11]/ 1.40269400e-02; /*CH2(S) (wt in kg) */
    phi[12] = y[12]/ 1.50349100e-02; /*CH3 (wt in kg) */
    phi[13] = y[13]/ 1.60428800e-02; /*CH4 (wt in kg) */
    phi[14] = y[14]/ 2.80104000e-02; /*CO (wt in kg) */
    phi[15] = y[15]/ 4.40098000e-02; /*CO2 (wt in kg) */
    phi[16] = y[16]/ 2.90183700e-02; /*HCO (wt in kg) */
    phi[17] = y[17]/ 3.00263400e-02; /*CH2O (wt in kg) */
    phi[18] = y[18]/ 3.10343100e-02; /*CH2OH (wt in kg) */
    phi[19] = y[19]/ 3.10343100e-02; /*CH3O (wt in kg) */
    phi[20] = y[20]/ 3.20422800e-02; /*CH3OH (wt in kg) */
    phi[21] = y[21]/ 2.50299700e-02; /*C2H (wt in kg) */
    phi[22] = y[22]/ 2.60379400e-02; /*C2H2 (wt in kg) */
    phi[23] = y[23]/ 2.70459100e-02; /*C2H3 (wt in kg) */
    phi[24] = y[24]/ 2.80538800e-02; /*C2H4 (wt in kg) */
    phi[25] = y[25]/ 2.90618500e-02; /*C2H5 (wt in kg) */
    phi[26] = y[26]/ 3.00698200e-02; /*C2H6 (wt in kg) */
    phi[27] = y[27]/ 4.10293700e-02; /*HCCO (wt in kg) */
    phi[28] = y[28]/ 4.20373400e-02; /*CH2CO (wt in kg) */
    phi[29] = y[29]/ 4.20373400e-02; /*HCCOH (wt in kg) */
    phi[30] = y[30]/ 2.80134000e-02; /*N2 (wt in kg) */
    phi[31] = y[31]/ 3.99480000e-02; /*AR (wt in kg) */

    return;
}


/*reverse of ytcr, useful for rate computations */
void fectyr_(double * c, double * rho, int * iwrk, double * rwrk, double * y)
{
    y[0] = c[0] * 2.015940 / (*rho); 
    y[1] = c[1] * 1.007970 / (*rho); 
    y[2] = c[2] * 15.999400 / (*rho); 
    y[3] = c[3] * 31.998800 / (*rho); 
    y[4] = c[4] * 17.007370 / (*rho); 
    y[5] = c[5] * 18.015340 / (*rho); 
    y[6] = c[6] * 33.006770 / (*rho); 
    y[7] = c[7] * 34.014740 / (*rho); 
    y[8] = c[8] * 12.011000 / (*rho); 
    y[9] = c[9] * 13.018970 / (*rho); 
    y[10] = c[10] * 14.026940 / (*rho); 
    y[11] = c[11] * 14.026940 / (*rho); 
    y[12] = c[12] * 15.034910 / (*rho); 
    y[13] = c[13] * 16.042880 / (*rho); 
    y[14] = c[14] * 28.010400 / (*rho); 
    y[15] = c[15] * 44.009800 / (*rho); 
    y[16] = c[16] * 29.018370 / (*rho); 
    y[17] = c[17] * 30.026340 / (*rho); 
    y[18] = c[18] * 31.034310 / (*rho); 
    y[19] = c[19] * 31.034310 / (*rho); 
    y[20] = c[20] * 32.042280 / (*rho); 
    y[21] = c[21] * 25.029970 / (*rho); 
    y[22] = c[22] * 26.037940 / (*rho); 
    y[23] = c[23] * 27.045910 / (*rho); 
    y[24] = c[24] * 28.053880 / (*rho); 
    y[25] = c[25] * 29.061850 / (*rho); 
    y[26] = c[26] * 30.069820 / (*rho); 
    y[27] = c[27] * 41.029370 / (*rho); 
    y[28] = c[28] * 42.037340 / (*rho); 
    y[29] = c[29] * 42.037340 / (*rho); 
    y[30] = c[30] * 28.013400 / (*rho); 
    y[31] = c[31] * 39.948000 / (*rho); 

    return;
}


/*ddebdf compatible right hand side of CV burner */
/*rwrk[0] and rwrk[1] should contain rho and ene respectively */
/*working variable phi contains specific mole numbers */
void fecvrhs_(double * time, double * phi, double * phidot, double * rwrk, int * iwrk)
{
    double rho,ene; /*CV Parameters */
    double y[32], wdot[32]; /*temporary storage */
    int i; /*Loop counter */
    double temperature,pressure; /*temporary var */
    rho = rwrk[0];
    ene = rwrk[1];
    fephity_(phi, iwrk, rwrk, y);
    feeytt_(&ene, y, iwrk, rwrk, &temperature);
    CKPY(&rho, &temperature,  y, iwrk, rwrk, &pressure);
    CKWYP(&pressure, &temperature,  y, iwrk, rwrk, wdot);
    for (i=0; i<32; ++i) phidot[i] = wdot[i] / (rho/1000.0); 

    return;
}


/*returns the dimensionality of the cv burner (number of species) */
int fecvdim_()
{
    return 32;
}


/*ddebdf compatible right hand side of ZND solver */
/*rwrk[0] : scaling factor for pressure */
/*rwrk[1] : preshock density (g/cc)  */
/*rwrk[2] : detonation velocity (cm/s)  */
/*solution vector: [P; rho; y0 ... ylast]  */
void fezndrhs_(double * time, double * z, double * zdot, double * rwrk, int * iwrk)
{
    double psc,rho1,udet; /*ZND Parameters */
    double wt[32], hms[32], wdot[32]; /*temporary storage */
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
    for (i=0; i<32; ++i) {
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
    return 35;
}


/*returns the name of the source mechanism file  */
char* femechfile_()
{
    return "";
}


/*returns the species number */
int fesymnum_(const char* s1)
{
    if (strcmp(s1, "H2")==0) return 0; 
    if (strcmp(s1, "H")==0) return 1; 
    if (strcmp(s1, "O")==0) return 2; 
    if (strcmp(s1, "O2")==0) return 3; 
    if (strcmp(s1, "OH")==0) return 4; 
    if (strcmp(s1, "H2O")==0) return 5; 
    if (strcmp(s1, "HO2")==0) return 6; 
    if (strcmp(s1, "H2O2")==0) return 7; 
    if (strcmp(s1, "C")==0) return 8; 
    if (strcmp(s1, "CH")==0) return 9; 
    if (strcmp(s1, "CH2")==0) return 10; 
    if (strcmp(s1, "CH2(S)")==0) return 11; 
    if (strcmp(s1, "CH3")==0) return 12; 
    if (strcmp(s1, "CH4")==0) return 13; 
    if (strcmp(s1, "CO")==0) return 14; 
    if (strcmp(s1, "CO2")==0) return 15; 
    if (strcmp(s1, "HCO")==0) return 16; 
    if (strcmp(s1, "CH2O")==0) return 17; 
    if (strcmp(s1, "CH2OH")==0) return 18; 
    if (strcmp(s1, "CH3O")==0) return 19; 
    if (strcmp(s1, "CH3OH")==0) return 20; 
    if (strcmp(s1, "C2H")==0) return 21; 
    if (strcmp(s1, "C2H2")==0) return 22; 
    if (strcmp(s1, "C2H3")==0) return 23; 
    if (strcmp(s1, "C2H4")==0) return 24; 
    if (strcmp(s1, "C2H5")==0) return 25; 
    if (strcmp(s1, "C2H6")==0) return 26; 
    if (strcmp(s1, "HCCO")==0) return 27; 
    if (strcmp(s1, "CH2CO")==0) return 28; 
    if (strcmp(s1, "HCCOH")==0) return 29; 
    if (strcmp(s1, "N2")==0) return 30; 
    if (strcmp(s1, "AR")==0) return 31; 
    /*species name not found */
    return -1;
}


/*returns the species name */
char* fesymname_(int sn)
{
    if (sn==0) return "H2"; 
    if (sn==1) return "H"; 
    if (sn==2) return "O"; 
    if (sn==3) return "O2"; 
    if (sn==4) return "OH"; 
    if (sn==5) return "H2O"; 
    if (sn==6) return "HO2"; 
    if (sn==7) return "H2O2"; 
    if (sn==8) return "C"; 
    if (sn==9) return "CH"; 
    if (sn==10) return "CH2"; 
    if (sn==11) return "CH2(S)"; 
    if (sn==12) return "CH3"; 
    if (sn==13) return "CH4"; 
    if (sn==14) return "CO"; 
    if (sn==15) return "CO2"; 
    if (sn==16) return "HCO"; 
    if (sn==17) return "CH2O"; 
    if (sn==18) return "CH2OH"; 
    if (sn==19) return "CH3O"; 
    if (sn==20) return "CH3OH"; 
    if (sn==21) return "C2H"; 
    if (sn==22) return "C2H2"; 
    if (sn==23) return "C2H3"; 
    if (sn==24) return "C2H4"; 
    if (sn==25) return "C2H5"; 
    if (sn==26) return "C2H6"; 
    if (sn==27) return "HCCO"; 
    if (sn==28) return "CH2CO"; 
    if (sn==29) return "HCCOH"; 
    if (sn==30) return "N2"; 
    if (sn==31) return "AR"; 
    /*species name not found */
    return "NOTFOUND";
}

/* End of file  */
