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
    *kk = 53;
    *ii = 325;
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
void CKSYME(int * kname, int * plenkname)
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
    for (i=0; i<lenkname*53; i++) {
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

    /* N  */
    kname[ 30*lenkname + 0 ] = 'N';
    kname[ 30*lenkname + 1 ] = ' ';

    /* NH  */
    kname[ 31*lenkname + 0 ] = 'N';
    kname[ 31*lenkname + 1 ] = 'H';
    kname[ 31*lenkname + 2 ] = ' ';

    /* NH2  */
    kname[ 32*lenkname + 0 ] = 'N';
    kname[ 32*lenkname + 1 ] = 'H';
    kname[ 32*lenkname + 2 ] = '2';
    kname[ 32*lenkname + 3 ] = ' ';

    /* NH3  */
    kname[ 33*lenkname + 0 ] = 'N';
    kname[ 33*lenkname + 1 ] = 'H';
    kname[ 33*lenkname + 2 ] = '3';
    kname[ 33*lenkname + 3 ] = ' ';

    /* NNH  */
    kname[ 34*lenkname + 0 ] = 'N';
    kname[ 34*lenkname + 1 ] = 'N';
    kname[ 34*lenkname + 2 ] = 'H';
    kname[ 34*lenkname + 3 ] = ' ';

    /* NO  */
    kname[ 35*lenkname + 0 ] = 'N';
    kname[ 35*lenkname + 1 ] = 'O';
    kname[ 35*lenkname + 2 ] = ' ';

    /* NO2  */
    kname[ 36*lenkname + 0 ] = 'N';
    kname[ 36*lenkname + 1 ] = 'O';
    kname[ 36*lenkname + 2 ] = '2';
    kname[ 36*lenkname + 3 ] = ' ';

    /* N2O  */
    kname[ 37*lenkname + 0 ] = 'N';
    kname[ 37*lenkname + 1 ] = '2';
    kname[ 37*lenkname + 2 ] = 'O';
    kname[ 37*lenkname + 3 ] = ' ';

    /* HNO  */
    kname[ 38*lenkname + 0 ] = 'H';
    kname[ 38*lenkname + 1 ] = 'N';
    kname[ 38*lenkname + 2 ] = 'O';
    kname[ 38*lenkname + 3 ] = ' ';

    /* CN  */
    kname[ 39*lenkname + 0 ] = 'C';
    kname[ 39*lenkname + 1 ] = 'N';
    kname[ 39*lenkname + 2 ] = ' ';

    /* HCN  */
    kname[ 40*lenkname + 0 ] = 'H';
    kname[ 40*lenkname + 1 ] = 'C';
    kname[ 40*lenkname + 2 ] = 'N';
    kname[ 40*lenkname + 3 ] = ' ';

    /* H2CN  */
    kname[ 41*lenkname + 0 ] = 'H';
    kname[ 41*lenkname + 1 ] = '2';
    kname[ 41*lenkname + 2 ] = 'C';
    kname[ 41*lenkname + 3 ] = 'N';
    kname[ 41*lenkname + 4 ] = ' ';

    /* HCNN  */
    kname[ 42*lenkname + 0 ] = 'H';
    kname[ 42*lenkname + 1 ] = 'C';
    kname[ 42*lenkname + 2 ] = 'N';
    kname[ 42*lenkname + 3 ] = 'N';
    kname[ 42*lenkname + 4 ] = ' ';

    /* HCNO  */
    kname[ 43*lenkname + 0 ] = 'H';
    kname[ 43*lenkname + 1 ] = 'C';
    kname[ 43*lenkname + 2 ] = 'N';
    kname[ 43*lenkname + 3 ] = 'O';
    kname[ 43*lenkname + 4 ] = ' ';

    /* HOCN  */
    kname[ 44*lenkname + 0 ] = 'H';
    kname[ 44*lenkname + 1 ] = 'O';
    kname[ 44*lenkname + 2 ] = 'C';
    kname[ 44*lenkname + 3 ] = 'N';
    kname[ 44*lenkname + 4 ] = ' ';

    /* HNCO  */
    kname[ 45*lenkname + 0 ] = 'H';
    kname[ 45*lenkname + 1 ] = 'N';
    kname[ 45*lenkname + 2 ] = 'C';
    kname[ 45*lenkname + 3 ] = 'O';
    kname[ 45*lenkname + 4 ] = ' ';

    /* NCO  */
    kname[ 46*lenkname + 0 ] = 'N';
    kname[ 46*lenkname + 1 ] = 'C';
    kname[ 46*lenkname + 2 ] = 'O';
    kname[ 46*lenkname + 3 ] = ' ';

    /* N2  */
    kname[ 47*lenkname + 0 ] = 'N';
    kname[ 47*lenkname + 1 ] = '2';
    kname[ 47*lenkname + 2 ] = ' ';

    /* AR  */
    kname[ 48*lenkname + 0 ] = 'A';
    kname[ 48*lenkname + 1 ] = 'R';
    kname[ 48*lenkname + 2 ] = ' ';

    /* C3H7  */
    kname[ 49*lenkname + 0 ] = 'C';
    kname[ 49*lenkname + 1 ] = '3';
    kname[ 49*lenkname + 2 ] = 'H';
    kname[ 49*lenkname + 3 ] = '7';
    kname[ 49*lenkname + 4 ] = ' ';

    /* C3H8  */
    kname[ 50*lenkname + 0 ] = 'C';
    kname[ 50*lenkname + 1 ] = '3';
    kname[ 50*lenkname + 2 ] = 'H';
    kname[ 50*lenkname + 3 ] = '8';
    kname[ 50*lenkname + 4 ] = ' ';

    /* CH2CHO  */
    kname[ 51*lenkname + 0 ] = 'C';
    kname[ 51*lenkname + 1 ] = 'H';
    kname[ 51*lenkname + 2 ] = '2';
    kname[ 51*lenkname + 3 ] = 'C';
    kname[ 51*lenkname + 4 ] = 'H';
    kname[ 51*lenkname + 5 ] = 'O';
    kname[ 51*lenkname + 6 ] = ' ';

    /* CH3CHO  */
    kname[ 52*lenkname + 0 ] = 'C';
    kname[ 52*lenkname + 1 ] = 'H';
    kname[ 52*lenkname + 2 ] = '3';
    kname[ 52*lenkname + 3 ] = 'C';
    kname[ 52*lenkname + 4 ] = 'H';
    kname[ 52*lenkname + 5 ] = 'O';
    kname[ 52*lenkname + 6 ] = ' ';

}


/* Returns R, Rc, Patm */
void CKRP(int * ickwrk, double * rckwrk, double * ru, double * ruc, double * pa)
{
     *ru  = 8.31451e+07; 
     //*ruc = 1.987; 
     *ruc = (*ru)/4.184e7; 
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
    XW += x[30]*14.006700; /*N */
    XW += x[31]*15.014670; /*NH */
    XW += x[32]*16.022640; /*NH2 */
    XW += x[33]*17.030610; /*NH3 */
    XW += x[34]*29.021370; /*NNH */
    XW += x[35]*30.006100; /*NO */
    XW += x[36]*46.005500; /*NO2 */
    XW += x[37]*44.012800; /*N2O */
    XW += x[38]*31.014070; /*HNO */
    XW += x[39]*26.017700; /*CN */
    XW += x[40]*27.025670; /*HCN */
    XW += x[41]*28.033640; /*H2CN */
    XW += x[42]*41.032370; /*HCNN */
    XW += x[43]*43.025070; /*HCNO */
    XW += x[44]*43.025070; /*HOCN */
    XW += x[45]*43.025070; /*HNCO */
    XW += x[46]*42.017100; /*NCO */
    XW += x[47]*28.013400; /*N2 */
    XW += x[48]*39.948000; /*AR */
    XW += x[49]*43.088790; /*C3H7 */
    XW += x[50]*44.096760; /*C3H8 */
    XW += x[51]*43.045310; /*CH2CHO */
    XW += x[52]*44.053280; /*CH3CHO */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

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
    YOW += y[30]/14.006700; /*N */
    YOW += y[31]/15.014670; /*NH */
    YOW += y[32]/16.022640; /*NH2 */
    YOW += y[33]/17.030610; /*NH3 */
    YOW += y[34]/29.021370; /*NNH */
    YOW += y[35]/30.006100; /*NO */
    YOW += y[36]/46.005500; /*NO2 */
    YOW += y[37]/44.012800; /*N2O */
    YOW += y[38]/31.014070; /*HNO */
    YOW += y[39]/26.017700; /*CN */
    YOW += y[40]/27.025670; /*HCN */
    YOW += y[41]/28.033640; /*H2CN */
    YOW += y[42]/41.032370; /*HCNN */
    YOW += y[43]/43.025070; /*HCNO */
    YOW += y[44]/43.025070; /*HOCN */
    YOW += y[45]/43.025070; /*HNCO */
    YOW += y[46]/42.017100; /*NCO */
    YOW += y[47]/28.013400; /*N2 */
    YOW += y[48]/39.948000; /*AR */
    YOW += y[49]/43.088790; /*C3H7 */
    YOW += y[50]/44.096760; /*C3H8 */
    YOW += y[51]/43.045310; /*CH2CHO */
    YOW += y[52]/44.053280; /*CH3CHO */
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
    W += c[30]*14.006700; /*N */
    W += c[31]*15.014670; /*NH */
    W += c[32]*16.022640; /*NH2 */
    W += c[33]*17.030610; /*NH3 */
    W += c[34]*29.021370; /*NNH */
    W += c[35]*30.006100; /*NO */
    W += c[36]*46.005500; /*NO2 */
    W += c[37]*44.012800; /*N2O */
    W += c[38]*31.014070; /*HNO */
    W += c[39]*26.017700; /*CN */
    W += c[40]*27.025670; /*HCN */
    W += c[41]*28.033640; /*H2CN */
    W += c[42]*41.032370; /*HCNN */
    W += c[43]*43.025070; /*HCNO */
    W += c[44]*43.025070; /*HOCN */
    W += c[45]*43.025070; /*HNCO */
    W += c[46]*42.017100; /*NCO */
    W += c[47]*28.013400; /*N2 */
    W += c[48]*39.948000; /*AR */
    W += c[49]*43.088790; /*C3H7 */
    W += c[50]*44.096760; /*C3H8 */
    W += c[51]*43.045310; /*CH2CHO */
    W += c[52]*44.053280; /*CH3CHO */

    for (id = 0; id < 53; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31451e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

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
    XW += x[30]*14.006700; /*N */
    XW += x[31]*15.014670; /*NH */
    XW += x[32]*16.022640; /*NH2 */
    XW += x[33]*17.030610; /*NH3 */
    XW += x[34]*29.021370; /*NNH */
    XW += x[35]*30.006100; /*NO */
    XW += x[36]*46.005500; /*NO2 */
    XW += x[37]*44.012800; /*N2O */
    XW += x[38]*31.014070; /*HNO */
    XW += x[39]*26.017700; /*CN */
    XW += x[40]*27.025670; /*HCN */
    XW += x[41]*28.033640; /*H2CN */
    XW += x[42]*41.032370; /*HCNN */
    XW += x[43]*43.025070; /*HCNO */
    XW += x[44]*43.025070; /*HOCN */
    XW += x[45]*43.025070; /*HNCO */
    XW += x[46]*42.017100; /*NCO */
    XW += x[47]*28.013400; /*N2 */
    XW += x[48]*39.948000; /*AR */
    XW += x[49]*43.088790; /*C3H7 */
    XW += x[50]*44.096760; /*C3H8 */
    XW += x[51]*43.045310; /*CH2CHO */
    XW += x[52]*44.053280; /*CH3CHO */
    *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

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
    YOW += y[30]/14.006700; /*N */
    YOW += y[31]/15.014670; /*NH */
    YOW += y[32]/16.022640; /*NH2 */
    YOW += y[33]/17.030610; /*NH3 */
    YOW += y[34]/29.021370; /*NNH */
    YOW += y[35]/30.006100; /*NO */
    YOW += y[36]/46.005500; /*NO2 */
    YOW += y[37]/44.012800; /*N2O */
    YOW += y[38]/31.014070; /*HNO */
    YOW += y[39]/26.017700; /*CN */
    YOW += y[40]/27.025670; /*HCN */
    YOW += y[41]/28.033640; /*H2CN */
    YOW += y[42]/41.032370; /*HCNN */
    YOW += y[43]/43.025070; /*HCNO */
    YOW += y[44]/43.025070; /*HOCN */
    YOW += y[45]/43.025070; /*HNCO */
    YOW += y[46]/42.017100; /*NCO */
    YOW += y[47]/28.013400; /*N2 */
    YOW += y[48]/39.948000; /*AR */
    YOW += y[49]/43.088790; /*C3H7 */
    YOW += y[50]/44.096760; /*C3H8 */
    YOW += y[51]/43.045310; /*CH2CHO */
    YOW += y[52]/44.053280; /*CH3CHO */
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
    W += c[30]*14.006700; /*N */
    W += c[31]*15.014670; /*NH */
    W += c[32]*16.022640; /*NH2 */
    W += c[33]*17.030610; /*NH3 */
    W += c[34]*29.021370; /*NNH */
    W += c[35]*30.006100; /*NO */
    W += c[36]*46.005500; /*NO2 */
    W += c[37]*44.012800; /*N2O */
    W += c[38]*31.014070; /*HNO */
    W += c[39]*26.017700; /*CN */
    W += c[40]*27.025670; /*HCN */
    W += c[41]*28.033640; /*H2CN */
    W += c[42]*41.032370; /*HCNN */
    W += c[43]*43.025070; /*HCNO */
    W += c[44]*43.025070; /*HOCN */
    W += c[45]*43.025070; /*HNCO */
    W += c[46]*42.017100; /*NCO */
    W += c[47]*28.013400; /*N2 */
    W += c[48]*39.948000; /*AR */
    W += c[49]*43.088790; /*C3H7 */
    W += c[50]*44.096760; /*C3H8 */
    W += c[51]*43.045310; /*CH2CHO */
    W += c[52]*44.053280; /*CH3CHO */

    for (id = 0; id < 53; ++id) {
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
    YOW += y[30]/14.006700; /*N */
    YOW += y[31]/15.014670; /*NH */
    YOW += y[32]/16.022640; /*NH2 */
    YOW += y[33]/17.030610; /*NH3 */
    YOW += y[34]/29.021370; /*NNH */
    YOW += y[35]/30.006100; /*NO */
    YOW += y[36]/46.005500; /*NO2 */
    YOW += y[37]/44.012800; /*N2O */
    YOW += y[38]/31.014070; /*HNO */
    YOW += y[39]/26.017700; /*CN */
    YOW += y[40]/27.025670; /*HCN */
    YOW += y[41]/28.033640; /*H2CN */
    YOW += y[42]/41.032370; /*HCNN */
    YOW += y[43]/43.025070; /*HCNO */
    YOW += y[44]/43.025070; /*HOCN */
    YOW += y[45]/43.025070; /*HNCO */
    YOW += y[46]/42.017100; /*NCO */
    YOW += y[47]/28.013400; /*N2 */
    YOW += y[48]/39.948000; /*AR */
    YOW += y[49]/43.088790; /*C3H7 */
    YOW += y[50]/44.096760; /*C3H8 */
    YOW += y[51]/43.045310; /*CH2CHO */
    YOW += y[52]/44.053280; /*CH3CHO */
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
    XW += x[30]*14.006700; /*N */
    XW += x[31]*15.014670; /*NH */
    XW += x[32]*16.022640; /*NH2 */
    XW += x[33]*17.030610; /*NH3 */
    XW += x[34]*29.021370; /*NNH */
    XW += x[35]*30.006100; /*NO */
    XW += x[36]*46.005500; /*NO2 */
    XW += x[37]*44.012800; /*N2O */
    XW += x[38]*31.014070; /*HNO */
    XW += x[39]*26.017700; /*CN */
    XW += x[40]*27.025670; /*HCN */
    XW += x[41]*28.033640; /*H2CN */
    XW += x[42]*41.032370; /*HCNN */
    XW += x[43]*43.025070; /*HCNO */
    XW += x[44]*43.025070; /*HOCN */
    XW += x[45]*43.025070; /*HNCO */
    XW += x[46]*42.017100; /*NCO */
    XW += x[47]*28.013400; /*N2 */
    XW += x[48]*39.948000; /*AR */
    XW += x[49]*43.088790; /*C3H7 */
    XW += x[50]*44.096760; /*C3H8 */
    XW += x[51]*43.045310; /*CH2CHO */
    XW += x[52]*44.053280; /*CH3CHO */
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
    W += c[30]*14.006700; /*N */
    W += c[31]*15.014670; /*NH */
    W += c[32]*16.022640; /*NH2 */
    W += c[33]*17.030610; /*NH3 */
    W += c[34]*29.021370; /*NNH */
    W += c[35]*30.006100; /*NO */
    W += c[36]*46.005500; /*NO2 */
    W += c[37]*44.012800; /*N2O */
    W += c[38]*31.014070; /*HNO */
    W += c[39]*26.017700; /*CN */
    W += c[40]*27.025670; /*HCN */
    W += c[41]*28.033640; /*H2CN */
    W += c[42]*41.032370; /*HCNN */
    W += c[43]*43.025070; /*HCNO */
    W += c[44]*43.025070; /*HOCN */
    W += c[45]*43.025070; /*HNCO */
    W += c[46]*42.017100; /*NCO */
    W += c[47]*28.013400; /*N2 */
    W += c[48]*39.948000; /*AR */
    W += c[49]*43.088790; /*C3H7 */
    W += c[50]*44.096760; /*C3H8 */
    W += c[51]*43.045310; /*CH2CHO */
    W += c[52]*44.053280; /*CH3CHO */

    for (id = 0; id < 53; ++id) {
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
    YOW += y[30]/14.006700; /*N */
    YOW += y[31]/15.014670; /*NH */
    YOW += y[32]/16.022640; /*NH2 */
    YOW += y[33]/17.030610; /*NH3 */
    YOW += y[34]/29.021370; /*NNH */
    YOW += y[35]/30.006100; /*NO */
    YOW += y[36]/46.005500; /*NO2 */
    YOW += y[37]/44.012800; /*N2O */
    YOW += y[38]/31.014070; /*HNO */
    YOW += y[39]/26.017700; /*CN */
    YOW += y[40]/27.025670; /*HCN */
    YOW += y[41]/28.033640; /*H2CN */
    YOW += y[42]/41.032370; /*HCNN */
    YOW += y[43]/43.025070; /*HCNO */
    YOW += y[44]/43.025070; /*HOCN */
    YOW += y[45]/43.025070; /*HNCO */
    YOW += y[46]/42.017100; /*NCO */
    YOW += y[47]/28.013400; /*N2 */
    YOW += y[48]/39.948000; /*AR */
    YOW += y[49]/43.088790; /*C3H7 */
    YOW += y[50]/44.096760; /*C3H8 */
    YOW += y[51]/43.045310; /*CH2CHO */
    YOW += y[52]/44.053280; /*CH3CHO */
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
    x[30] = y[30]/(14.006700*YOW); 
    x[31] = y[31]/(15.014670*YOW); 
    x[32] = y[32]/(16.022640*YOW); 
    x[33] = y[33]/(17.030610*YOW); 
    x[34] = y[34]/(29.021370*YOW); 
    x[35] = y[35]/(30.006100*YOW); 
    x[36] = y[36]/(46.005500*YOW); 
    x[37] = y[37]/(44.012800*YOW); 
    x[38] = y[38]/(31.014070*YOW); 
    x[39] = y[39]/(26.017700*YOW); 
    x[40] = y[40]/(27.025670*YOW); 
    x[41] = y[41]/(28.033640*YOW); 
    x[42] = y[42]/(41.032370*YOW); 
    x[43] = y[43]/(43.025070*YOW); 
    x[44] = y[44]/(43.025070*YOW); 
    x[45] = y[45]/(43.025070*YOW); 
    x[46] = y[46]/(42.017100*YOW); 
    x[47] = y[47]/(28.013400*YOW); 
    x[48] = y[48]/(39.948000*YOW); 
    x[49] = y[49]/(43.088790*YOW); 
    x[50] = y[50]/(44.096760*YOW); 
    x[51] = y[51]/(43.045310*YOW); 
    x[52] = y[52]/(44.053280*YOW); 

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
    YOW += y[30]/14.006700; /*N */
    YOW += y[31]/15.014670; /*NH */
    YOW += y[32]/16.022640; /*NH2 */
    YOW += y[33]/17.030610; /*NH3 */
    YOW += y[34]/29.021370; /*NNH */
    YOW += y[35]/30.006100; /*NO */
    YOW += y[36]/46.005500; /*NO2 */
    YOW += y[37]/44.012800; /*N2O */
    YOW += y[38]/31.014070; /*HNO */
    YOW += y[39]/26.017700; /*CN */
    YOW += y[40]/27.025670; /*HCN */
    YOW += y[41]/28.033640; /*H2CN */
    YOW += y[42]/41.032370; /*HCNN */
    YOW += y[43]/43.025070; /*HCNO */
    YOW += y[44]/43.025070; /*HOCN */
    YOW += y[45]/43.025070; /*HNCO */
    YOW += y[46]/42.017100; /*NCO */
    YOW += y[47]/28.013400; /*N2 */
    YOW += y[48]/39.948000; /*AR */
    YOW += y[49]/43.088790; /*C3H7 */
    YOW += y[50]/44.096760; /*C3H8 */
    YOW += y[51]/43.045310; /*CH2CHO */
    YOW += y[52]/44.053280; /*CH3CHO */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
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
    c[30] = PWORT * y[30]/14.006700; 
    c[31] = PWORT * y[31]/15.014670; 
    c[32] = PWORT * y[32]/16.022640; 
    c[33] = PWORT * y[33]/17.030610; 
    c[34] = PWORT * y[34]/29.021370; 
    c[35] = PWORT * y[35]/30.006100; 
    c[36] = PWORT * y[36]/46.005500; 
    c[37] = PWORT * y[37]/44.012800; 
    c[38] = PWORT * y[38]/31.014070; 
    c[39] = PWORT * y[39]/26.017700; 
    c[40] = PWORT * y[40]/27.025670; 
    c[41] = PWORT * y[41]/28.033640; 
    c[42] = PWORT * y[42]/41.032370; 
    c[43] = PWORT * y[43]/43.025070; 
    c[44] = PWORT * y[44]/43.025070; 
    c[45] = PWORT * y[45]/43.025070; 
    c[46] = PWORT * y[46]/42.017100; 
    c[47] = PWORT * y[47]/28.013400; 
    c[48] = PWORT * y[48]/39.948000; 
    c[49] = PWORT * y[49]/43.088790; 
    c[50] = PWORT * y[50]/44.096760; 
    c[51] = PWORT * y[51]/43.045310; 
    c[52] = PWORT * y[52]/44.053280; 

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
    c[30] = (*rho) * y[30]/14.006700; 
    c[31] = (*rho) * y[31]/15.014670; 
    c[32] = (*rho) * y[32]/16.022640; 
    c[33] = (*rho) * y[33]/17.030610; 
    c[34] = (*rho) * y[34]/29.021370; 
    c[35] = (*rho) * y[35]/30.006100; 
    c[36] = (*rho) * y[36]/46.005500; 
    c[37] = (*rho) * y[37]/44.012800; 
    c[38] = (*rho) * y[38]/31.014070; 
    c[39] = (*rho) * y[39]/26.017700; 
    c[40] = (*rho) * y[40]/27.025670; 
    c[41] = (*rho) * y[41]/28.033640; 
    c[42] = (*rho) * y[42]/41.032370; 
    c[43] = (*rho) * y[43]/43.025070; 
    c[44] = (*rho) * y[44]/43.025070; 
    c[45] = (*rho) * y[45]/43.025070; 
    c[46] = (*rho) * y[46]/42.017100; 
    c[47] = (*rho) * y[47]/28.013400; 
    c[48] = (*rho) * y[48]/39.948000; 
    c[49] = (*rho) * y[49]/43.088790; 
    c[50] = (*rho) * y[50]/44.096760; 
    c[51] = (*rho) * y[51]/43.045310; 
    c[52] = (*rho) * y[52]/44.053280; 

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
    XW += x[30]*14.006700; /*N */
    XW += x[31]*15.014670; /*NH */
    XW += x[32]*16.022640; /*NH2 */
    XW += x[33]*17.030610; /*NH3 */
    XW += x[34]*29.021370; /*NNH */
    XW += x[35]*30.006100; /*NO */
    XW += x[36]*46.005500; /*NO2 */
    XW += x[37]*44.012800; /*N2O */
    XW += x[38]*31.014070; /*HNO */
    XW += x[39]*26.017700; /*CN */
    XW += x[40]*27.025670; /*HCN */
    XW += x[41]*28.033640; /*H2CN */
    XW += x[42]*41.032370; /*HCNN */
    XW += x[43]*43.025070; /*HCNO */
    XW += x[44]*43.025070; /*HOCN */
    XW += x[45]*43.025070; /*HNCO */
    XW += x[46]*42.017100; /*NCO */
    XW += x[47]*28.013400; /*N2 */
    XW += x[48]*39.948000; /*AR */
    XW += x[49]*43.088790; /*C3H7 */
    XW += x[50]*44.096760; /*C3H8 */
    XW += x[51]*43.045310; /*CH2CHO */
    XW += x[52]*44.053280; /*CH3CHO */
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
    y[30] = x[30]*14.006700/XW; 
    y[31] = x[31]*15.014670/XW; 
    y[32] = x[32]*16.022640/XW; 
    y[33] = x[33]*17.030610/XW; 
    y[34] = x[34]*29.021370/XW; 
    y[35] = x[35]*30.006100/XW; 
    y[36] = x[36]*46.005500/XW; 
    y[37] = x[37]*44.012800/XW; 
    y[38] = x[38]*31.014070/XW; 
    y[39] = x[39]*26.017700/XW; 
    y[40] = x[40]*27.025670/XW; 
    y[41] = x[41]*28.033640/XW; 
    y[42] = x[42]*41.032370/XW; 
    y[43] = x[43]*43.025070/XW; 
    y[44] = x[44]*43.025070/XW; 
    y[45] = x[45]*43.025070/XW; 
    y[46] = x[46]*42.017100/XW; 
    y[47] = x[47]*28.013400/XW; 
    y[48] = x[48]*39.948000/XW; 
    y[49] = x[49]*43.088790/XW; 
    y[50] = x[50]*44.096760/XW; 
    y[51] = x[51]*43.045310/XW; 
    y[52] = x[52]*44.053280/XW; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 53; ++id) {
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
    XW += x[30]*14.006700; /*N */
    XW += x[31]*15.014670; /*NH */
    XW += x[32]*16.022640; /*NH2 */
    XW += x[33]*17.030610; /*NH3 */
    XW += x[34]*29.021370; /*NNH */
    XW += x[35]*30.006100; /*NO */
    XW += x[36]*46.005500; /*NO2 */
    XW += x[37]*44.012800; /*N2O */
    XW += x[38]*31.014070; /*HNO */
    XW += x[39]*26.017700; /*CN */
    XW += x[40]*27.025670; /*HCN */
    XW += x[41]*28.033640; /*H2CN */
    XW += x[42]*41.032370; /*HCNN */
    XW += x[43]*43.025070; /*HCNO */
    XW += x[44]*43.025070; /*HOCN */
    XW += x[45]*43.025070; /*HNCO */
    XW += x[46]*42.017100; /*NCO */
    XW += x[47]*28.013400; /*N2 */
    XW += x[48]*39.948000; /*AR */
    XW += x[49]*43.088790; /*C3H7 */
    XW += x[50]*44.096760; /*C3H8 */
    XW += x[51]*43.045310; /*CH2CHO */
    XW += x[52]*44.053280; /*CH3CHO */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 53; ++id) {
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
    for (id = 0; id < 53; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    for (id = 0; id < 53; ++id) {
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
    CW += c[30]*14.006700; /*N */
    CW += c[31]*15.014670; /*NH */
    CW += c[32]*16.022640; /*NH2 */
    CW += c[33]*17.030610; /*NH3 */
    CW += c[34]*29.021370; /*NNH */
    CW += c[35]*30.006100; /*NO */
    CW += c[36]*46.005500; /*NO2 */
    CW += c[37]*44.012800; /*N2O */
    CW += c[38]*31.014070; /*HNO */
    CW += c[39]*26.017700; /*CN */
    CW += c[40]*27.025670; /*HCN */
    CW += c[41]*28.033640; /*H2CN */
    CW += c[42]*41.032370; /*HCNN */
    CW += c[43]*43.025070; /*HCNO */
    CW += c[44]*43.025070; /*HOCN */
    CW += c[45]*43.025070; /*HNCO */
    CW += c[46]*42.017100; /*NCO */
    CW += c[47]*28.013400; /*N2 */
    CW += c[48]*39.948000; /*AR */
    CW += c[49]*43.088790; /*C3H7 */
    CW += c[50]*44.096760; /*C3H8 */
    CW += c[51]*43.045310; /*CH2CHO */
    CW += c[52]*44.053280; /*CH3CHO */
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
    y[30] = c[30]*14.006700/CW; 
    y[31] = c[31]*15.014670/CW; 
    y[32] = c[32]*16.022640/CW; 
    y[33] = c[33]*17.030610/CW; 
    y[34] = c[34]*29.021370/CW; 
    y[35] = c[35]*30.006100/CW; 
    y[36] = c[36]*46.005500/CW; 
    y[37] = c[37]*44.012800/CW; 
    y[38] = c[38]*31.014070/CW; 
    y[39] = c[39]*26.017700/CW; 
    y[40] = c[40]*27.025670/CW; 
    y[41] = c[41]*28.033640/CW; 
    y[42] = c[42]*41.032370/CW; 
    y[43] = c[43]*43.025070/CW; 
    y[44] = c[44]*43.025070/CW; 
    y[45] = c[45]*43.025070/CW; 
    y[46] = c[46]*42.017100/CW; 
    y[47] = c[47]*28.013400/CW; 
    y[48] = c[48]*39.948000/CW; 
    y[49] = c[49]*43.088790/CW; 
    y[50] = c[50]*44.096760/CW; 
    y[51] = c[51]*43.045310/CW; 
    y[52] = c[52]*44.053280/CW; 

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
    for (id = 0; id < 53; ++id) {
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
    for (id = 0; id < 53; ++id) {
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
    for (id = 0; id < 53; ++id) {
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
    for (id = 0; id < 53; ++id) {
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
    for (id = 0; id < 53; ++id) {
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
    for (id = 0; id < 53; ++id) {
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
    for (id = 0; id < 53; ++id) {
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
    cvms[30] *= 5935730.757423; /*N */
    cvms[31] *= 5537251.234959; /*NH */
    cvms[32] *= 5188907.695611; /*NH2 */
    cvms[33] *= 4881798.127020; /*NH3 */
    cvms[34] *= 2864785.501167; /*NNH */
    cvms[35] *= 2770769.943445; /*NO */
    cvms[36] *= 1807175.229049; /*NO2 */
    cvms[37] *= 1888995.928457; /*N2O */
    cvms[38] *= 2680718.783442; /*HNO */
    cvms[39] *= 3195516.898112; /*CN */
    cvms[40] *= 3076334.462753; /*HCN */
    cvms[41] *= 2965722.610407; /*H2CN */
    cvms[42] *= 2026205.164362; /*HCNN */
    cvms[43] *= 1932361.760248; /*HCNO */
    cvms[44] *= 1932361.760248; /*HOCN */
    cvms[45] *= 1932361.760248; /*HNCO */
    cvms[46] *= 1978718.188547; /*NCO */
    cvms[47] *= 2967865.378712; /*N2 */
    cvms[48] *= 2081205.567237; /*AR */
    cvms[49] *= 1929504.170342; /*C3H7 */
    cvms[50] *= 1885399.290107; /*C3H8 */
    cvms[51] *= 1931453.159473; /*CH2CHO */
    cvms[52] *= 1887260.154068; /*CH3CHO */
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
    cpms[30] *= 5935730.757423; /*N */
    cpms[31] *= 5537251.234959; /*NH */
    cpms[32] *= 5188907.695611; /*NH2 */
    cpms[33] *= 4881798.127020; /*NH3 */
    cpms[34] *= 2864785.501167; /*NNH */
    cpms[35] *= 2770769.943445; /*NO */
    cpms[36] *= 1807175.229049; /*NO2 */
    cpms[37] *= 1888995.928457; /*N2O */
    cpms[38] *= 2680718.783442; /*HNO */
    cpms[39] *= 3195516.898112; /*CN */
    cpms[40] *= 3076334.462753; /*HCN */
    cpms[41] *= 2965722.610407; /*H2CN */
    cpms[42] *= 2026205.164362; /*HCNN */
    cpms[43] *= 1932361.760248; /*HCNO */
    cpms[44] *= 1932361.760248; /*HOCN */
    cpms[45] *= 1932361.760248; /*HNCO */
    cpms[46] *= 1978718.188547; /*NCO */
    cpms[47] *= 2967865.378712; /*N2 */
    cpms[48] *= 2081205.567237; /*AR */
    cpms[49] *= 1929504.170342; /*C3H7 */
    cpms[50] *= 1885399.290107; /*C3H8 */
    cpms[51] *= 1931453.159473; /*CH2CHO */
    cpms[52] *= 1887260.154068; /*CH3CHO */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double *T, int * iwrk, double * rwrk, double * ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
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
    ums[30] *= RT/14.006700; /*N */
    ums[31] *= RT/15.014670; /*NH */
    ums[32] *= RT/16.022640; /*NH2 */
    ums[33] *= RT/17.030610; /*NH3 */
    ums[34] *= RT/29.021370; /*NNH */
    ums[35] *= RT/30.006100; /*NO */
    ums[36] *= RT/46.005500; /*NO2 */
    ums[37] *= RT/44.012800; /*N2O */
    ums[38] *= RT/31.014070; /*HNO */
    ums[39] *= RT/26.017700; /*CN */
    ums[40] *= RT/27.025670; /*HCN */
    ums[41] *= RT/28.033640; /*H2CN */
    ums[42] *= RT/41.032370; /*HCNN */
    ums[43] *= RT/43.025070; /*HCNO */
    ums[44] *= RT/43.025070; /*HOCN */
    ums[45] *= RT/43.025070; /*HNCO */
    ums[46] *= RT/42.017100; /*NCO */
    ums[47] *= RT/28.013400; /*N2 */
    ums[48] *= RT/39.948000; /*AR */
    ums[49] *= RT/43.088790; /*C3H7 */
    ums[50] *= RT/44.096760; /*C3H8 */
    ums[51] *= RT/43.045310; /*CH2CHO */
    ums[52] *= RT/44.053280; /*CH3CHO */
}


/*Returns enthalpy in mass units (Eq 27.) */
void CKHMS(double *T, int * iwrk, double * rwrk, double * hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
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
    hms[30] *= RT/14.006700; /*N */
    hms[31] *= RT/15.014670; /*NH */
    hms[32] *= RT/16.022640; /*NH2 */
    hms[33] *= RT/17.030610; /*NH3 */
    hms[34] *= RT/29.021370; /*NNH */
    hms[35] *= RT/30.006100; /*NO */
    hms[36] *= RT/46.005500; /*NO2 */
    hms[37] *= RT/44.012800; /*N2O */
    hms[38] *= RT/31.014070; /*HNO */
    hms[39] *= RT/26.017700; /*CN */
    hms[40] *= RT/27.025670; /*HCN */
    hms[41] *= RT/28.033640; /*H2CN */
    hms[42] *= RT/41.032370; /*HCNN */
    hms[43] *= RT/43.025070; /*HCNO */
    hms[44] *= RT/43.025070; /*HOCN */
    hms[45] *= RT/43.025070; /*HNCO */
    hms[46] *= RT/42.017100; /*NCO */
    hms[47] *= RT/28.013400; /*N2 */
    hms[48] *= RT/39.948000; /*AR */
    hms[49] *= RT/43.088790; /*C3H7 */
    hms[50] *= RT/44.096760; /*C3H8 */
    hms[51] *= RT/43.045310; /*CH2CHO */
    hms[52] *= RT/44.053280; /*CH3CHO */
}


/*Returns gibbs in mass units (Eq 31.) */
void CKGMS(double *T, int * iwrk, double * rwrk, double * gms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
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
    gms[30] *= RT/14.006700; /*N */
    gms[31] *= RT/15.014670; /*NH */
    gms[32] *= RT/16.022640; /*NH2 */
    gms[33] *= RT/17.030610; /*NH3 */
    gms[34] *= RT/29.021370; /*NNH */
    gms[35] *= RT/30.006100; /*NO */
    gms[36] *= RT/46.005500; /*NO2 */
    gms[37] *= RT/44.012800; /*N2O */
    gms[38] *= RT/31.014070; /*HNO */
    gms[39] *= RT/26.017700; /*CN */
    gms[40] *= RT/27.025670; /*HCN */
    gms[41] *= RT/28.033640; /*H2CN */
    gms[42] *= RT/41.032370; /*HCNN */
    gms[43] *= RT/43.025070; /*HCNO */
    gms[44] *= RT/43.025070; /*HOCN */
    gms[45] *= RT/43.025070; /*HNCO */
    gms[46] *= RT/42.017100; /*NCO */
    gms[47] *= RT/28.013400; /*N2 */
    gms[48] *= RT/39.948000; /*AR */
    gms[49] *= RT/43.088790; /*C3H7 */
    gms[50] *= RT/44.096760; /*C3H8 */
    gms[51] *= RT/43.045310; /*CH2CHO */
    gms[52] *= RT/44.053280; /*CH3CHO */
}


/*Returns helmholtz in mass units (Eq 32.) */
void CKAMS(double *T, int * iwrk, double * rwrk, double * ams)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
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
    ams[30] *= RT/14.006700; /*N */
    ams[31] *= RT/15.014670; /*NH */
    ams[32] *= RT/16.022640; /*NH2 */
    ams[33] *= RT/17.030610; /*NH3 */
    ams[34] *= RT/29.021370; /*NNH */
    ams[35] *= RT/30.006100; /*NO */
    ams[36] *= RT/46.005500; /*NO2 */
    ams[37] *= RT/44.012800; /*N2O */
    ams[38] *= RT/31.014070; /*HNO */
    ams[39] *= RT/26.017700; /*CN */
    ams[40] *= RT/27.025670; /*HCN */
    ams[41] *= RT/28.033640; /*H2CN */
    ams[42] *= RT/41.032370; /*HCNN */
    ams[43] *= RT/43.025070; /*HCNO */
    ams[44] *= RT/43.025070; /*HOCN */
    ams[45] *= RT/43.025070; /*HNCO */
    ams[46] *= RT/42.017100; /*NCO */
    ams[47] *= RT/28.013400; /*N2 */
    ams[48] *= RT/39.948000; /*AR */
    ams[49] *= RT/43.088790; /*C3H7 */
    ams[50] *= RT/44.096760; /*C3H8 */
    ams[51] *= RT/43.045310; /*CH2CHO */
    ams[52] *= RT/44.053280; /*CH3CHO */
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
    sms[30] *= 5935730.757423; /*N */
    sms[31] *= 5537251.234959; /*NH */
    sms[32] *= 5188907.695611; /*NH2 */
    sms[33] *= 4881798.127020; /*NH3 */
    sms[34] *= 2864785.501167; /*NNH */
    sms[35] *= 2770769.943445; /*NO */
    sms[36] *= 1807175.229049; /*NO2 */
    sms[37] *= 1888995.928457; /*N2O */
    sms[38] *= 2680718.783442; /*HNO */
    sms[39] *= 3195516.898112; /*CN */
    sms[40] *= 3076334.462753; /*HCN */
    sms[41] *= 2965722.610407; /*H2CN */
    sms[42] *= 2026205.164362; /*HCNN */
    sms[43] *= 1932361.760248; /*HCNO */
    sms[44] *= 1932361.760248; /*HOCN */
    sms[45] *= 1932361.760248; /*HNCO */
    sms[46] *= 1978718.188547; /*NCO */
    sms[47] *= 2967865.378712; /*N2 */
    sms[48] *= 2081205.567237; /*AR */
    sms[49] *= 1929504.170342; /*C3H7 */
    sms[50] *= 1885399.290107; /*C3H8 */
    sms[51] *= 1931453.159473; /*CH2CHO */
    sms[52] *= 1887260.154068; /*CH3CHO */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double *T, double *x, int * iwrk, double * rwrk, double * cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[53]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 53; ++id) {
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
    double cpor[53]; /* temporary storage */
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
    result += cpor[30]*y[30]/14.0067; /*N */
    result += cpor[31]*y[31]/15.0147; /*NH */
    result += cpor[32]*y[32]/16.0226; /*NH2 */
    result += cpor[33]*y[33]/17.0306; /*NH3 */
    result += cpor[34]*y[34]/29.0214; /*NNH */
    result += cpor[35]*y[35]/30.0061; /*NO */
    result += cpor[36]*y[36]/46.0055; /*NO2 */
    result += cpor[37]*y[37]/44.0128; /*N2O */
    result += cpor[38]*y[38]/31.0141; /*HNO */
    result += cpor[39]*y[39]/26.0177; /*CN */
    result += cpor[40]*y[40]/27.0257; /*HCN */
    result += cpor[41]*y[41]/28.0336; /*H2CN */
    result += cpor[42]*y[42]/41.0324; /*HCNN */
    result += cpor[43]*y[43]/43.0251; /*HCNO */
    result += cpor[44]*y[44]/43.0251; /*HOCN */
    result += cpor[45]*y[45]/43.0251; /*HNCO */
    result += cpor[46]*y[46]/42.0171; /*NCO */
    result += cpor[47]*y[47]/28.0134; /*N2 */
    result += cpor[48]*y[48]/39.948; /*AR */
    result += cpor[49]*y[49]/43.0888; /*C3H7 */
    result += cpor[50]*y[50]/44.0968; /*C3H8 */
    result += cpor[51]*y[51]/43.0453; /*CH2CHO */
    result += cpor[52]*y[52]/44.0533; /*CH3CHO */

    *cpbs = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double *T, double *x, int * iwrk, double * rwrk, double * cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[53]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 53; ++id) {
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
    double cvor[53]; /* temporary storage */
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
    result += cvor[30]*y[30]/14.0067; /*N */
    result += cvor[31]*y[31]/15.0147; /*NH */
    result += cvor[32]*y[32]/16.0226; /*NH2 */
    result += cvor[33]*y[33]/17.0306; /*NH3 */
    result += cvor[34]*y[34]/29.0214; /*NNH */
    result += cvor[35]*y[35]/30.0061; /*NO */
    result += cvor[36]*y[36]/46.0055; /*NO2 */
    result += cvor[37]*y[37]/44.0128; /*N2O */
    result += cvor[38]*y[38]/31.0141; /*HNO */
    result += cvor[39]*y[39]/26.0177; /*CN */
    result += cvor[40]*y[40]/27.0257; /*HCN */
    result += cvor[41]*y[41]/28.0336; /*H2CN */
    result += cvor[42]*y[42]/41.0324; /*HCNN */
    result += cvor[43]*y[43]/43.0251; /*HCNO */
    result += cvor[44]*y[44]/43.0251; /*HOCN */
    result += cvor[45]*y[45]/43.0251; /*HNCO */
    result += cvor[46]*y[46]/42.0171; /*NCO */
    result += cvor[47]*y[47]/28.0134; /*N2 */
    result += cvor[48]*y[48]/39.948; /*AR */
    result += cvor[49]*y[49]/43.0888; /*C3H7 */
    result += cvor[50]*y[50]/44.0968; /*C3H8 */
    result += cvor[51]*y[51]/43.0453; /*CH2CHO */
    result += cvor[52]*y[52]/44.0533; /*CH3CHO */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double *T, double *x, int * iwrk, double * rwrk, double * hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[53]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 53; ++id) {
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
    double hml[53]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
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
    result += y[30]*hml[30]/14.006700; /*N */
    result += y[31]*hml[31]/15.014670; /*NH */
    result += y[32]*hml[32]/16.022640; /*NH2 */
    result += y[33]*hml[33]/17.030610; /*NH3 */
    result += y[34]*hml[34]/29.021370; /*NNH */
    result += y[35]*hml[35]/30.006100; /*NO */
    result += y[36]*hml[36]/46.005500; /*NO2 */
    result += y[37]*hml[37]/44.012800; /*N2O */
    result += y[38]*hml[38]/31.014070; /*HNO */
    result += y[39]*hml[39]/26.017700; /*CN */
    result += y[40]*hml[40]/27.025670; /*HCN */
    result += y[41]*hml[41]/28.033640; /*H2CN */
    result += y[42]*hml[42]/41.032370; /*HCNN */
    result += y[43]*hml[43]/43.025070; /*HCNO */
    result += y[44]*hml[44]/43.025070; /*HOCN */
    result += y[45]*hml[45]/43.025070; /*HNCO */
    result += y[46]*hml[46]/42.017100; /*NCO */
    result += y[47]*hml[47]/28.013400; /*N2 */
    result += y[48]*hml[48]/39.948000; /*AR */
    result += y[49]*hml[49]/43.088790; /*C3H7 */
    result += y[50]*hml[50]/44.096760; /*C3H8 */
    result += y[51]*hml[51]/43.045310; /*CH2CHO */
    result += y[52]*hml[52]/44.053280; /*CH3CHO */

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double *T, double *x, int * iwrk, double * rwrk, double * ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[53]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 53; ++id) {
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
    double ums[53]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
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
    result += y[30]*ums[30]/14.006700; /*N */
    result += y[31]*ums[31]/15.014670; /*NH */
    result += y[32]*ums[32]/16.022640; /*NH2 */
    result += y[33]*ums[33]/17.030610; /*NH3 */
    result += y[34]*ums[34]/29.021370; /*NNH */
    result += y[35]*ums[35]/30.006100; /*NO */
    result += y[36]*ums[36]/46.005500; /*NO2 */
    result += y[37]*ums[37]/44.012800; /*N2O */
    result += y[38]*ums[38]/31.014070; /*HNO */
    result += y[39]*ums[39]/26.017700; /*CN */
    result += y[40]*ums[40]/27.025670; /*HCN */
    result += y[41]*ums[41]/28.033640; /*H2CN */
    result += y[42]*ums[42]/41.032370; /*HCNN */
    result += y[43]*ums[43]/43.025070; /*HCNO */
    result += y[44]*ums[44]/43.025070; /*HOCN */
    result += y[45]*ums[45]/43.025070; /*HNCO */
    result += y[46]*ums[46]/42.017100; /*NCO */
    result += y[47]*ums[47]/28.013400; /*N2 */
    result += y[48]*ums[48]/39.948000; /*AR */
    result += y[49]*ums[49]/43.088790; /*C3H7 */
    result += y[50]*ums[50]/44.096760; /*C3H8 */
    result += y[51]*ums[51]/43.045310; /*CH2CHO */
    result += y[52]*ums[52]/44.053280; /*CH3CHO */

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
    double sor[53]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 53; ++id) {
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
    double sor[53]; /* temporary storage */
    double x[53]; /* need a ytx conversion */
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
    YOW += y[30]/14.006700; /*N */
    YOW += y[31]/15.014670; /*NH */
    YOW += y[32]/16.022640; /*NH2 */
    YOW += y[33]/17.030610; /*NH3 */
    YOW += y[34]/29.021370; /*NNH */
    YOW += y[35]/30.006100; /*NO */
    YOW += y[36]/46.005500; /*NO2 */
    YOW += y[37]/44.012800; /*N2O */
    YOW += y[38]/31.014070; /*HNO */
    YOW += y[39]/26.017700; /*CN */
    YOW += y[40]/27.025670; /*HCN */
    YOW += y[41]/28.033640; /*H2CN */
    YOW += y[42]/41.032370; /*HCNN */
    YOW += y[43]/43.025070; /*HCNO */
    YOW += y[44]/43.025070; /*HOCN */
    YOW += y[45]/43.025070; /*HNCO */
    YOW += y[46]/42.017100; /*NCO */
    YOW += y[47]/28.013400; /*N2 */
    YOW += y[48]/39.948000; /*AR */
    YOW += y[49]/43.088790; /*C3H7 */
    YOW += y[50]/44.096760; /*C3H8 */
    YOW += y[51]/43.045310; /*CH2CHO */
    YOW += y[52]/44.053280; /*CH3CHO */
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
    x[30] = y[30]/(14.006700*YOW); 
    x[31] = y[31]/(15.014670*YOW); 
    x[32] = y[32]/(16.022640*YOW); 
    x[33] = y[33]/(17.030610*YOW); 
    x[34] = y[34]/(29.021370*YOW); 
    x[35] = y[35]/(30.006100*YOW); 
    x[36] = y[36]/(46.005500*YOW); 
    x[37] = y[37]/(44.012800*YOW); 
    x[38] = y[38]/(31.014070*YOW); 
    x[39] = y[39]/(26.017700*YOW); 
    x[40] = y[40]/(27.025670*YOW); 
    x[41] = y[41]/(28.033640*YOW); 
    x[42] = y[42]/(41.032370*YOW); 
    x[43] = y[43]/(43.025070*YOW); 
    x[44] = y[44]/(43.025070*YOW); 
    x[45] = y[45]/(43.025070*YOW); 
    x[46] = y[46]/(42.017100*YOW); 
    x[47] = y[47]/(28.013400*YOW); 
    x[48] = y[48]/(39.948000*YOW); 
    x[49] = y[49]/(43.088790*YOW); 
    x[50] = y[50]/(44.096760*YOW); 
    x[51] = y[51]/(43.045310*YOW); 
    x[52] = y[52]/(44.053280*YOW); 
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
    result += x[32]*(sor[32]-log((x[32]+1e-100))-logPratio);
    result += x[33]*(sor[33]-log((x[33]+1e-100))-logPratio);
    result += x[34]*(sor[34]-log((x[34]+1e-100))-logPratio);
    result += x[35]*(sor[35]-log((x[35]+1e-100))-logPratio);
    result += x[36]*(sor[36]-log((x[36]+1e-100))-logPratio);
    result += x[37]*(sor[37]-log((x[37]+1e-100))-logPratio);
    result += x[38]*(sor[38]-log((x[38]+1e-100))-logPratio);
    result += x[39]*(sor[39]-log((x[39]+1e-100))-logPratio);
    result += x[40]*(sor[40]-log((x[40]+1e-100))-logPratio);
    result += x[41]*(sor[41]-log((x[41]+1e-100))-logPratio);
    result += x[42]*(sor[42]-log((x[42]+1e-100))-logPratio);
    result += x[43]*(sor[43]-log((x[43]+1e-100))-logPratio);
    result += x[44]*(sor[44]-log((x[44]+1e-100))-logPratio);
    result += x[45]*(sor[45]-log((x[45]+1e-100))-logPratio);
    result += x[46]*(sor[46]-log((x[46]+1e-100))-logPratio);
    result += x[47]*(sor[47]-log((x[47]+1e-100))-logPratio);
    result += x[48]*(sor[48]-log((x[48]+1e-100))-logPratio);
    result += x[49]*(sor[49]-log((x[49]+1e-100))-logPratio);
    result += x[50]*(sor[50]-log((x[50]+1e-100))-logPratio);
    result += x[51]*(sor[51]-log((x[51]+1e-100))-logPratio);
    result += x[52]*(sor[52]-log((x[52]+1e-100))-logPratio);
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
    double gort[53]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 53; ++id) {
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
    double gort[53]; /* temporary storage */
    double x[53]; /* need a ytx conversion */
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
    YOW += y[30]/14.006700; /*N */
    YOW += y[31]/15.014670; /*NH */
    YOW += y[32]/16.022640; /*NH2 */
    YOW += y[33]/17.030610; /*NH3 */
    YOW += y[34]/29.021370; /*NNH */
    YOW += y[35]/30.006100; /*NO */
    YOW += y[36]/46.005500; /*NO2 */
    YOW += y[37]/44.012800; /*N2O */
    YOW += y[38]/31.014070; /*HNO */
    YOW += y[39]/26.017700; /*CN */
    YOW += y[40]/27.025670; /*HCN */
    YOW += y[41]/28.033640; /*H2CN */
    YOW += y[42]/41.032370; /*HCNN */
    YOW += y[43]/43.025070; /*HCNO */
    YOW += y[44]/43.025070; /*HOCN */
    YOW += y[45]/43.025070; /*HNCO */
    YOW += y[46]/42.017100; /*NCO */
    YOW += y[47]/28.013400; /*N2 */
    YOW += y[48]/39.948000; /*AR */
    YOW += y[49]/43.088790; /*C3H7 */
    YOW += y[50]/44.096760; /*C3H8 */
    YOW += y[51]/43.045310; /*CH2CHO */
    YOW += y[52]/44.053280; /*CH3CHO */
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
    x[30] = y[30]/(14.006700*YOW); 
    x[31] = y[31]/(15.014670*YOW); 
    x[32] = y[32]/(16.022640*YOW); 
    x[33] = y[33]/(17.030610*YOW); 
    x[34] = y[34]/(29.021370*YOW); 
    x[35] = y[35]/(30.006100*YOW); 
    x[36] = y[36]/(46.005500*YOW); 
    x[37] = y[37]/(44.012800*YOW); 
    x[38] = y[38]/(31.014070*YOW); 
    x[39] = y[39]/(26.017700*YOW); 
    x[40] = y[40]/(27.025670*YOW); 
    x[41] = y[41]/(28.033640*YOW); 
    x[42] = y[42]/(41.032370*YOW); 
    x[43] = y[43]/(43.025070*YOW); 
    x[44] = y[44]/(43.025070*YOW); 
    x[45] = y[45]/(43.025070*YOW); 
    x[46] = y[46]/(42.017100*YOW); 
    x[47] = y[47]/(28.013400*YOW); 
    x[48] = y[48]/(39.948000*YOW); 
    x[49] = y[49]/(43.088790*YOW); 
    x[50] = y[50]/(44.096760*YOW); 
    x[51] = y[51]/(43.045310*YOW); 
    x[52] = y[52]/(44.053280*YOW); 
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
    result += x[32]*(gort[32]+log((x[32]+1e-100))+logPratio);
    result += x[33]*(gort[33]+log((x[33]+1e-100))+logPratio);
    result += x[34]*(gort[34]+log((x[34]+1e-100))+logPratio);
    result += x[35]*(gort[35]+log((x[35]+1e-100))+logPratio);
    result += x[36]*(gort[36]+log((x[36]+1e-100))+logPratio);
    result += x[37]*(gort[37]+log((x[37]+1e-100))+logPratio);
    result += x[38]*(gort[38]+log((x[38]+1e-100))+logPratio);
    result += x[39]*(gort[39]+log((x[39]+1e-100))+logPratio);
    result += x[40]*(gort[40]+log((x[40]+1e-100))+logPratio);
    result += x[41]*(gort[41]+log((x[41]+1e-100))+logPratio);
    result += x[42]*(gort[42]+log((x[42]+1e-100))+logPratio);
    result += x[43]*(gort[43]+log((x[43]+1e-100))+logPratio);
    result += x[44]*(gort[44]+log((x[44]+1e-100))+logPratio);
    result += x[45]*(gort[45]+log((x[45]+1e-100))+logPratio);
    result += x[46]*(gort[46]+log((x[46]+1e-100))+logPratio);
    result += x[47]*(gort[47]+log((x[47]+1e-100))+logPratio);
    result += x[48]*(gort[48]+log((x[48]+1e-100))+logPratio);
    result += x[49]*(gort[49]+log((x[49]+1e-100))+logPratio);
    result += x[50]*(gort[50]+log((x[50]+1e-100))+logPratio);
    result += x[51]*(gort[51]+log((x[51]+1e-100))+logPratio);
    result += x[52]*(gort[52]+log((x[52]+1e-100))+logPratio);
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
    double aort[53]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 53; ++id) {
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
    double aort[53]; /* temporary storage */
    double x[53]; /* need a ytx conversion */
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
    YOW += y[30]/14.006700; /*N */
    YOW += y[31]/15.014670; /*NH */
    YOW += y[32]/16.022640; /*NH2 */
    YOW += y[33]/17.030610; /*NH3 */
    YOW += y[34]/29.021370; /*NNH */
    YOW += y[35]/30.006100; /*NO */
    YOW += y[36]/46.005500; /*NO2 */
    YOW += y[37]/44.012800; /*N2O */
    YOW += y[38]/31.014070; /*HNO */
    YOW += y[39]/26.017700; /*CN */
    YOW += y[40]/27.025670; /*HCN */
    YOW += y[41]/28.033640; /*H2CN */
    YOW += y[42]/41.032370; /*HCNN */
    YOW += y[43]/43.025070; /*HCNO */
    YOW += y[44]/43.025070; /*HOCN */
    YOW += y[45]/43.025070; /*HNCO */
    YOW += y[46]/42.017100; /*NCO */
    YOW += y[47]/28.013400; /*N2 */
    YOW += y[48]/39.948000; /*AR */
    YOW += y[49]/43.088790; /*C3H7 */
    YOW += y[50]/44.096760; /*C3H8 */
    YOW += y[51]/43.045310; /*CH2CHO */
    YOW += y[52]/44.053280; /*CH3CHO */
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
    x[30] = y[30]/(14.006700*YOW); 
    x[31] = y[31]/(15.014670*YOW); 
    x[32] = y[32]/(16.022640*YOW); 
    x[33] = y[33]/(17.030610*YOW); 
    x[34] = y[34]/(29.021370*YOW); 
    x[35] = y[35]/(30.006100*YOW); 
    x[36] = y[36]/(46.005500*YOW); 
    x[37] = y[37]/(44.012800*YOW); 
    x[38] = y[38]/(31.014070*YOW); 
    x[39] = y[39]/(26.017700*YOW); 
    x[40] = y[40]/(27.025670*YOW); 
    x[41] = y[41]/(28.033640*YOW); 
    x[42] = y[42]/(41.032370*YOW); 
    x[43] = y[43]/(43.025070*YOW); 
    x[44] = y[44]/(43.025070*YOW); 
    x[45] = y[45]/(43.025070*YOW); 
    x[46] = y[46]/(42.017100*YOW); 
    x[47] = y[47]/(28.013400*YOW); 
    x[48] = y[48]/(39.948000*YOW); 
    x[49] = y[49]/(43.088790*YOW); 
    x[50] = y[50]/(44.096760*YOW); 
    x[51] = y[51]/(43.045310*YOW); 
    x[52] = y[52]/(44.053280*YOW); 
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
    result += x[32]*(aort[32]+log((x[32]+1e-100))+logPratio);
    result += x[33]*(aort[33]+log((x[33]+1e-100))+logPratio);
    result += x[34]*(aort[34]+log((x[34]+1e-100))+logPratio);
    result += x[35]*(aort[35]+log((x[35]+1e-100))+logPratio);
    result += x[36]*(aort[36]+log((x[36]+1e-100))+logPratio);
    result += x[37]*(aort[37]+log((x[37]+1e-100))+logPratio);
    result += x[38]*(aort[38]+log((x[38]+1e-100))+logPratio);
    result += x[39]*(aort[39]+log((x[39]+1e-100))+logPratio);
    result += x[40]*(aort[40]+log((x[40]+1e-100))+logPratio);
    result += x[41]*(aort[41]+log((x[41]+1e-100))+logPratio);
    result += x[42]*(aort[42]+log((x[42]+1e-100))+logPratio);
    result += x[43]*(aort[43]+log((x[43]+1e-100))+logPratio);
    result += x[44]*(aort[44]+log((x[44]+1e-100))+logPratio);
    result += x[45]*(aort[45]+log((x[45]+1e-100))+logPratio);
    result += x[46]*(aort[46]+log((x[46]+1e-100))+logPratio);
    result += x[47]*(aort[47]+log((x[47]+1e-100))+logPratio);
    result += x[48]*(aort[48]+log((x[48]+1e-100))+logPratio);
    result += x[49]*(aort[49]+log((x[49]+1e-100))+logPratio);
    result += x[50]*(aort[50]+log((x[50]+1e-100))+logPratio);
    result += x[51]*(aort[51]+log((x[51]+1e-100))+logPratio);
    result += x[52]*(aort[52]+log((x[52]+1e-100))+logPratio);
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double * T, double * C, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 53; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 53; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[53]; /*temporary storage */
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
    YOW += y[30]/14.006700; /*N */
    YOW += y[31]/15.014670; /*NH */
    YOW += y[32]/16.022640; /*NH2 */
    YOW += y[33]/17.030610; /*NH3 */
    YOW += y[34]/29.021370; /*NNH */
    YOW += y[35]/30.006100; /*NO */
    YOW += y[36]/46.005500; /*NO2 */
    YOW += y[37]/44.012800; /*N2O */
    YOW += y[38]/31.014070; /*HNO */
    YOW += y[39]/26.017700; /*CN */
    YOW += y[40]/27.025670; /*HCN */
    YOW += y[41]/28.033640; /*H2CN */
    YOW += y[42]/41.032370; /*HCNN */
    YOW += y[43]/43.025070; /*HCNO */
    YOW += y[44]/43.025070; /*HOCN */
    YOW += y[45]/43.025070; /*HNCO */
    YOW += y[46]/42.017100; /*NCO */
    YOW += y[47]/28.013400; /*N2 */
    YOW += y[48]/39.948000; /*AR */
    YOW += y[49]/43.088790; /*C3H7 */
    YOW += y[50]/44.096760; /*C3H8 */
    YOW += y[51]/43.045310; /*CH2CHO */
    YOW += y[52]/44.053280; /*CH3CHO */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
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
    c[30] = PWORT * y[30]/14.006700; 
    c[31] = PWORT * y[31]/15.014670; 
    c[32] = PWORT * y[32]/16.022640; 
    c[33] = PWORT * y[33]/17.030610; 
    c[34] = PWORT * y[34]/29.021370; 
    c[35] = PWORT * y[35]/30.006100; 
    c[36] = PWORT * y[36]/46.005500; 
    c[37] = PWORT * y[37]/44.012800; 
    c[38] = PWORT * y[38]/31.014070; 
    c[39] = PWORT * y[39]/26.017700; 
    c[40] = PWORT * y[40]/27.025670; 
    c[41] = PWORT * y[41]/28.033640; 
    c[42] = PWORT * y[42]/41.032370; 
    c[43] = PWORT * y[43]/43.025070; 
    c[44] = PWORT * y[44]/43.025070; 
    c[45] = PWORT * y[45]/43.025070; 
    c[46] = PWORT * y[46]/42.017100; 
    c[47] = PWORT * y[47]/28.013400; 
    c[48] = PWORT * y[48]/39.948000; 
    c[49] = PWORT * y[49]/43.088790; 
    c[50] = PWORT * y[50]/44.096760; 
    c[51] = PWORT * y[51]/43.045310; 
    c[52] = PWORT * y[52]/44.053280; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 53; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[53]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 53; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 53; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[53]; /*temporary storage */
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
    c[30] = 1e6 * (*rho) * y[30]/14.006700; 
    c[31] = 1e6 * (*rho) * y[31]/15.014670; 
    c[32] = 1e6 * (*rho) * y[32]/16.022640; 
    c[33] = 1e6 * (*rho) * y[33]/17.030610; 
    c[34] = 1e6 * (*rho) * y[34]/29.021370; 
    c[35] = 1e6 * (*rho) * y[35]/30.006100; 
    c[36] = 1e6 * (*rho) * y[36]/46.005500; 
    c[37] = 1e6 * (*rho) * y[37]/44.012800; 
    c[38] = 1e6 * (*rho) * y[38]/31.014070; 
    c[39] = 1e6 * (*rho) * y[39]/26.017700; 
    c[40] = 1e6 * (*rho) * y[40]/27.025670; 
    c[41] = 1e6 * (*rho) * y[41]/28.033640; 
    c[42] = 1e6 * (*rho) * y[42]/41.032370; 
    c[43] = 1e6 * (*rho) * y[43]/43.025070; 
    c[44] = 1e6 * (*rho) * y[44]/43.025070; 
    c[45] = 1e6 * (*rho) * y[45]/43.025070; 
    c[46] = 1e6 * (*rho) * y[46]/42.017100; 
    c[47] = 1e6 * (*rho) * y[47]/28.013400; 
    c[48] = 1e6 * (*rho) * y[48]/39.948000; 
    c[49] = 1e6 * (*rho) * y[49]/43.088790; 
    c[50] = 1e6 * (*rho) * y[50]/44.096760; 
    c[51] = 1e6 * (*rho) * y[51]/43.045310; 
    c[52] = 1e6 * (*rho) * y[52]/44.053280; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 53; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[53]; /*temporary storage */
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
    XW += x[30]*14.006700; /*N */
    XW += x[31]*15.014670; /*NH */
    XW += x[32]*16.022640; /*NH2 */
    XW += x[33]*17.030610; /*NH3 */
    XW += x[34]*29.021370; /*NNH */
    XW += x[35]*30.006100; /*NO */
    XW += x[36]*46.005500; /*NO2 */
    XW += x[37]*44.012800; /*N2O */
    XW += x[38]*31.014070; /*HNO */
    XW += x[39]*26.017700; /*CN */
    XW += x[40]*27.025670; /*HCN */
    XW += x[41]*28.033640; /*H2CN */
    XW += x[42]*41.032370; /*HCNN */
    XW += x[43]*43.025070; /*HCNO */
    XW += x[44]*43.025070; /*HOCN */
    XW += x[45]*43.025070; /*HNCO */
    XW += x[46]*42.017100; /*NCO */
    XW += x[47]*28.013400; /*N2 */
    XW += x[48]*39.948000; /*AR */
    XW += x[49]*43.088790; /*C3H7 */
    XW += x[50]*44.096760; /*C3H8 */
    XW += x[51]*43.045310; /*CH2CHO */
    XW += x[52]*44.053280; /*CH3CHO */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 53; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 53; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double * T, double * C, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 53; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 53; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 325; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[53]; /*temporary storage */
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
    YOW += y[30]/14.006700; /*N */
    YOW += y[31]/15.014670; /*NH */
    YOW += y[32]/16.022640; /*NH2 */
    YOW += y[33]/17.030610; /*NH3 */
    YOW += y[34]/29.021370; /*NNH */
    YOW += y[35]/30.006100; /*NO */
    YOW += y[36]/46.005500; /*NO2 */
    YOW += y[37]/44.012800; /*N2O */
    YOW += y[38]/31.014070; /*HNO */
    YOW += y[39]/26.017700; /*CN */
    YOW += y[40]/27.025670; /*HCN */
    YOW += y[41]/28.033640; /*H2CN */
    YOW += y[42]/41.032370; /*HCNN */
    YOW += y[43]/43.025070; /*HCNO */
    YOW += y[44]/43.025070; /*HOCN */
    YOW += y[45]/43.025070; /*HNCO */
    YOW += y[46]/42.017100; /*NCO */
    YOW += y[47]/28.013400; /*N2 */
    YOW += y[48]/39.948000; /*AR */
    YOW += y[49]/43.088790; /*C3H7 */
    YOW += y[50]/44.096760; /*C3H8 */
    YOW += y[51]/43.045310; /*CH2CHO */
    YOW += y[52]/44.053280; /*CH3CHO */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
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
    c[30] = PWORT * y[30]/14.006700; 
    c[31] = PWORT * y[31]/15.014670; 
    c[32] = PWORT * y[32]/16.022640; 
    c[33] = PWORT * y[33]/17.030610; 
    c[34] = PWORT * y[34]/29.021370; 
    c[35] = PWORT * y[35]/30.006100; 
    c[36] = PWORT * y[36]/46.005500; 
    c[37] = PWORT * y[37]/44.012800; 
    c[38] = PWORT * y[38]/31.014070; 
    c[39] = PWORT * y[39]/26.017700; 
    c[40] = PWORT * y[40]/27.025670; 
    c[41] = PWORT * y[41]/28.033640; 
    c[42] = PWORT * y[42]/41.032370; 
    c[43] = PWORT * y[43]/43.025070; 
    c[44] = PWORT * y[44]/43.025070; 
    c[45] = PWORT * y[45]/43.025070; 
    c[46] = PWORT * y[46]/42.017100; 
    c[47] = PWORT * y[47]/28.013400; 
    c[48] = PWORT * y[48]/39.948000; 
    c[49] = PWORT * y[49]/43.088790; 
    c[50] = PWORT * y[50]/44.096760; 
    c[51] = PWORT * y[51]/43.045310; 
    c[52] = PWORT * y[52]/44.053280; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 325; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[53]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 53; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 325; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[53]; /*temporary storage */
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
    c[30] = 1e6 * (*rho) * y[30]/14.006700; 
    c[31] = 1e6 * (*rho) * y[31]/15.014670; 
    c[32] = 1e6 * (*rho) * y[32]/16.022640; 
    c[33] = 1e6 * (*rho) * y[33]/17.030610; 
    c[34] = 1e6 * (*rho) * y[34]/29.021370; 
    c[35] = 1e6 * (*rho) * y[35]/30.006100; 
    c[36] = 1e6 * (*rho) * y[36]/46.005500; 
    c[37] = 1e6 * (*rho) * y[37]/44.012800; 
    c[38] = 1e6 * (*rho) * y[38]/31.014070; 
    c[39] = 1e6 * (*rho) * y[39]/26.017700; 
    c[40] = 1e6 * (*rho) * y[40]/27.025670; 
    c[41] = 1e6 * (*rho) * y[41]/28.033640; 
    c[42] = 1e6 * (*rho) * y[42]/41.032370; 
    c[43] = 1e6 * (*rho) * y[43]/43.025070; 
    c[44] = 1e6 * (*rho) * y[44]/43.025070; 
    c[45] = 1e6 * (*rho) * y[45]/43.025070; 
    c[46] = 1e6 * (*rho) * y[46]/42.017100; 
    c[47] = 1e6 * (*rho) * y[47]/28.013400; 
    c[48] = 1e6 * (*rho) * y[48]/39.948000; 
    c[49] = 1e6 * (*rho) * y[49]/43.088790; 
    c[50] = 1e6 * (*rho) * y[50]/44.096760; 
    c[51] = 1e6 * (*rho) * y[51]/43.045310; 
    c[52] = 1e6 * (*rho) * y[52]/44.053280; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 325; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[53]; /*temporary storage */
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
    XW += x[30]*14.006700; /*N */
    XW += x[31]*15.014670; /*NH */
    XW += x[32]*16.022640; /*NH2 */
    XW += x[33]*17.030610; /*NH3 */
    XW += x[34]*29.021370; /*NNH */
    XW += x[35]*30.006100; /*NO */
    XW += x[36]*46.005500; /*NO2 */
    XW += x[37]*44.012800; /*N2O */
    XW += x[38]*31.014070; /*HNO */
    XW += x[39]*26.017700; /*CN */
    XW += x[40]*27.025670; /*HCN */
    XW += x[41]*28.033640; /*H2CN */
    XW += x[42]*41.032370; /*HCNN */
    XW += x[43]*43.025070; /*HCNO */
    XW += x[44]*43.025070; /*HOCN */
    XW += x[45]*43.025070; /*HNCO */
    XW += x[46]*42.017100; /*NCO */
    XW += x[47]*28.013400; /*N2 */
    XW += x[48]*39.948000; /*AR */
    XW += x[49]*43.088790; /*C3H7 */
    XW += x[50]*44.096760; /*C3H8 */
    XW += x[51]*43.045310; /*CH2CHO */
    XW += x[52]*44.053280; /*CH3CHO */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 53; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 325; ++id) {
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
    for (id = 0; id < 53 * 325; ++ id) {
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

    /*reaction 12: O + CO (+M) <=> CO2 (+M) */
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
    nuki[ 47 * kd + 35 ] += -1 ;
    nuki[ 6 * kd + 35 ] += +1 ;
    nuki[ 47 * kd + 35 ] += +1 ;

    /*reaction 37: H + O2 + AR <=> HO2 + AR */
    nuki[ 1 * kd + 36 ] += -1 ;
    nuki[ 3 * kd + 36 ] += -1 ;
    nuki[ 48 * kd + 36 ] += -1 ;
    nuki[ 6 * kd + 36 ] += +1 ;
    nuki[ 48 * kd + 36 ] += +1 ;

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

    /*reaction 135: CH2 + O2 => OH + H + CO */
    nuki[ 10 * kd + 134 ] += -1 ;
    nuki[ 3 * kd + 134 ] += -1 ;
    nuki[ 4 * kd + 134 ] += +1 ;
    nuki[ 1 * kd + 134 ] += +1 ;
    nuki[ 14 * kd + 134 ] += +1 ;

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
    nuki[ 47 * kd + 141 ] += -1 ;
    nuki[ 10 * kd + 141 ] += +1 ;
    nuki[ 47 * kd + 141 ] += +1 ;

    /*reaction 143: CH2(S) + AR <=> CH2 + AR */
    nuki[ 11 * kd + 142 ] += -1 ;
    nuki[ 48 * kd + 142 ] += -1 ;
    nuki[ 10 * kd + 142 ] += +1 ;
    nuki[ 48 * kd + 142 ] += +1 ;

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

    /*reaction 178: N + NO <=> N2 + O */
    nuki[ 30 * kd + 177 ] += -1 ;
    nuki[ 35 * kd + 177 ] += -1 ;
    nuki[ 47 * kd + 177 ] += +1 ;
    nuki[ 2 * kd + 177 ] += +1 ;

    /*reaction 179: N + O2 <=> NO + O */
    nuki[ 30 * kd + 178 ] += -1 ;
    nuki[ 3 * kd + 178 ] += -1 ;
    nuki[ 35 * kd + 178 ] += +1 ;
    nuki[ 2 * kd + 178 ] += +1 ;

    /*reaction 180: N + OH <=> NO + H */
    nuki[ 30 * kd + 179 ] += -1 ;
    nuki[ 4 * kd + 179 ] += -1 ;
    nuki[ 35 * kd + 179 ] += +1 ;
    nuki[ 1 * kd + 179 ] += +1 ;

    /*reaction 181: N2O + O <=> N2 + O2 */
    nuki[ 37 * kd + 180 ] += -1 ;
    nuki[ 2 * kd + 180 ] += -1 ;
    nuki[ 47 * kd + 180 ] += +1 ;
    nuki[ 3 * kd + 180 ] += +1 ;

    /*reaction 182: N2O + O <=> 2 NO */
    nuki[ 37 * kd + 181 ] += -1 ;
    nuki[ 2 * kd + 181 ] += -1 ;
    nuki[ 35 * kd + 181 ] += +2 ;

    /*reaction 183: N2O + H <=> N2 + OH */
    nuki[ 37 * kd + 182 ] += -1 ;
    nuki[ 1 * kd + 182 ] += -1 ;
    nuki[ 47 * kd + 182 ] += +1 ;
    nuki[ 4 * kd + 182 ] += +1 ;

    /*reaction 184: N2O + OH <=> N2 + HO2 */
    nuki[ 37 * kd + 183 ] += -1 ;
    nuki[ 4 * kd + 183 ] += -1 ;
    nuki[ 47 * kd + 183 ] += +1 ;
    nuki[ 6 * kd + 183 ] += +1 ;

    /*reaction 185: N2O (+M) <=> N2 + O (+M) */
    nuki[ 37 * kd + 184 ] += -1 ;
    nuki[ 47 * kd + 184 ] += +1 ;
    nuki[ 2 * kd + 184 ] += +1 ;

    /*reaction 186: HO2 + NO <=> NO2 + OH */
    nuki[ 6 * kd + 185 ] += -1 ;
    nuki[ 35 * kd + 185 ] += -1 ;
    nuki[ 36 * kd + 185 ] += +1 ;
    nuki[ 4 * kd + 185 ] += +1 ;

    /*reaction 187: NO + O + M <=> NO2 + M */
    nuki[ 35 * kd + 186 ] += -1 ;
    nuki[ 2 * kd + 186 ] += -1 ;
    nuki[ 36 * kd + 186 ] += +1 ;

    /*reaction 188: NO2 + O <=> NO + O2 */
    nuki[ 36 * kd + 187 ] += -1 ;
    nuki[ 2 * kd + 187 ] += -1 ;
    nuki[ 35 * kd + 187 ] += +1 ;
    nuki[ 3 * kd + 187 ] += +1 ;

    /*reaction 189: NO2 + H <=> NO + OH */
    nuki[ 36 * kd + 188 ] += -1 ;
    nuki[ 1 * kd + 188 ] += -1 ;
    nuki[ 35 * kd + 188 ] += +1 ;
    nuki[ 4 * kd + 188 ] += +1 ;

    /*reaction 190: NH + O <=> NO + H */
    nuki[ 31 * kd + 189 ] += -1 ;
    nuki[ 2 * kd + 189 ] += -1 ;
    nuki[ 35 * kd + 189 ] += +1 ;
    nuki[ 1 * kd + 189 ] += +1 ;

    /*reaction 191: NH + H <=> N + H2 */
    nuki[ 31 * kd + 190 ] += -1 ;
    nuki[ 1 * kd + 190 ] += -1 ;
    nuki[ 30 * kd + 190 ] += +1 ;
    nuki[ 0 * kd + 190 ] += +1 ;

    /*reaction 192: NH + OH <=> HNO + H */
    nuki[ 31 * kd + 191 ] += -1 ;
    nuki[ 4 * kd + 191 ] += -1 ;
    nuki[ 38 * kd + 191 ] += +1 ;
    nuki[ 1 * kd + 191 ] += +1 ;

    /*reaction 193: NH + OH <=> N + H2O */
    nuki[ 31 * kd + 192 ] += -1 ;
    nuki[ 4 * kd + 192 ] += -1 ;
    nuki[ 30 * kd + 192 ] += +1 ;
    nuki[ 5 * kd + 192 ] += +1 ;

    /*reaction 194: NH + O2 <=> HNO + O */
    nuki[ 31 * kd + 193 ] += -1 ;
    nuki[ 3 * kd + 193 ] += -1 ;
    nuki[ 38 * kd + 193 ] += +1 ;
    nuki[ 2 * kd + 193 ] += +1 ;

    /*reaction 195: NH + O2 <=> NO + OH */
    nuki[ 31 * kd + 194 ] += -1 ;
    nuki[ 3 * kd + 194 ] += -1 ;
    nuki[ 35 * kd + 194 ] += +1 ;
    nuki[ 4 * kd + 194 ] += +1 ;

    /*reaction 196: NH + N <=> N2 + H */
    nuki[ 31 * kd + 195 ] += -1 ;
    nuki[ 30 * kd + 195 ] += -1 ;
    nuki[ 47 * kd + 195 ] += +1 ;
    nuki[ 1 * kd + 195 ] += +1 ;

    /*reaction 197: NH + H2O <=> HNO + H2 */
    nuki[ 31 * kd + 196 ] += -1 ;
    nuki[ 5 * kd + 196 ] += -1 ;
    nuki[ 38 * kd + 196 ] += +1 ;
    nuki[ 0 * kd + 196 ] += +1 ;

    /*reaction 198: NH + NO <=> N2 + OH */
    nuki[ 31 * kd + 197 ] += -1 ;
    nuki[ 35 * kd + 197 ] += -1 ;
    nuki[ 47 * kd + 197 ] += +1 ;
    nuki[ 4 * kd + 197 ] += +1 ;

    /*reaction 199: NH + NO <=> N2O + H */
    nuki[ 31 * kd + 198 ] += -1 ;
    nuki[ 35 * kd + 198 ] += -1 ;
    nuki[ 37 * kd + 198 ] += +1 ;
    nuki[ 1 * kd + 198 ] += +1 ;

    /*reaction 200: NH2 + O <=> OH + NH */
    nuki[ 32 * kd + 199 ] += -1 ;
    nuki[ 2 * kd + 199 ] += -1 ;
    nuki[ 4 * kd + 199 ] += +1 ;
    nuki[ 31 * kd + 199 ] += +1 ;

    /*reaction 201: NH2 + O <=> H + HNO */
    nuki[ 32 * kd + 200 ] += -1 ;
    nuki[ 2 * kd + 200 ] += -1 ;
    nuki[ 1 * kd + 200 ] += +1 ;
    nuki[ 38 * kd + 200 ] += +1 ;

    /*reaction 202: NH2 + H <=> NH + H2 */
    nuki[ 32 * kd + 201 ] += -1 ;
    nuki[ 1 * kd + 201 ] += -1 ;
    nuki[ 31 * kd + 201 ] += +1 ;
    nuki[ 0 * kd + 201 ] += +1 ;

    /*reaction 203: NH2 + OH <=> NH + H2O */
    nuki[ 32 * kd + 202 ] += -1 ;
    nuki[ 4 * kd + 202 ] += -1 ;
    nuki[ 31 * kd + 202 ] += +1 ;
    nuki[ 5 * kd + 202 ] += +1 ;

    /*reaction 204: NNH <=> N2 + H */
    nuki[ 34 * kd + 203 ] += -1 ;
    nuki[ 47 * kd + 203 ] += +1 ;
    nuki[ 1 * kd + 203 ] += +1 ;

    /*reaction 205: NNH + M <=> N2 + H + M */
    nuki[ 34 * kd + 204 ] += -1 ;
    nuki[ 47 * kd + 204 ] += +1 ;
    nuki[ 1 * kd + 204 ] += +1 ;

    /*reaction 206: NNH + O2 <=> HO2 + N2 */
    nuki[ 34 * kd + 205 ] += -1 ;
    nuki[ 3 * kd + 205 ] += -1 ;
    nuki[ 6 * kd + 205 ] += +1 ;
    nuki[ 47 * kd + 205 ] += +1 ;

    /*reaction 207: NNH + O <=> OH + N2 */
    nuki[ 34 * kd + 206 ] += -1 ;
    nuki[ 2 * kd + 206 ] += -1 ;
    nuki[ 4 * kd + 206 ] += +1 ;
    nuki[ 47 * kd + 206 ] += +1 ;

    /*reaction 208: NNH + O <=> NH + NO */
    nuki[ 34 * kd + 207 ] += -1 ;
    nuki[ 2 * kd + 207 ] += -1 ;
    nuki[ 31 * kd + 207 ] += +1 ;
    nuki[ 35 * kd + 207 ] += +1 ;

    /*reaction 209: NNH + H <=> H2 + N2 */
    nuki[ 34 * kd + 208 ] += -1 ;
    nuki[ 1 * kd + 208 ] += -1 ;
    nuki[ 0 * kd + 208 ] += +1 ;
    nuki[ 47 * kd + 208 ] += +1 ;

    /*reaction 210: NNH + OH <=> H2O + N2 */
    nuki[ 34 * kd + 209 ] += -1 ;
    nuki[ 4 * kd + 209 ] += -1 ;
    nuki[ 5 * kd + 209 ] += +1 ;
    nuki[ 47 * kd + 209 ] += +1 ;

    /*reaction 211: NNH + CH3 <=> CH4 + N2 */
    nuki[ 34 * kd + 210 ] += -1 ;
    nuki[ 12 * kd + 210 ] += -1 ;
    nuki[ 13 * kd + 210 ] += +1 ;
    nuki[ 47 * kd + 210 ] += +1 ;

    /*reaction 212: H + NO + M <=> HNO + M */
    nuki[ 1 * kd + 211 ] += -1 ;
    nuki[ 35 * kd + 211 ] += -1 ;
    nuki[ 38 * kd + 211 ] += +1 ;

    /*reaction 213: HNO + O <=> NO + OH */
    nuki[ 38 * kd + 212 ] += -1 ;
    nuki[ 2 * kd + 212 ] += -1 ;
    nuki[ 35 * kd + 212 ] += +1 ;
    nuki[ 4 * kd + 212 ] += +1 ;

    /*reaction 214: HNO + H <=> H2 + NO */
    nuki[ 38 * kd + 213 ] += -1 ;
    nuki[ 1 * kd + 213 ] += -1 ;
    nuki[ 0 * kd + 213 ] += +1 ;
    nuki[ 35 * kd + 213 ] += +1 ;

    /*reaction 215: HNO + OH <=> NO + H2O */
    nuki[ 38 * kd + 214 ] += -1 ;
    nuki[ 4 * kd + 214 ] += -1 ;
    nuki[ 35 * kd + 214 ] += +1 ;
    nuki[ 5 * kd + 214 ] += +1 ;

    /*reaction 216: HNO + O2 <=> HO2 + NO */
    nuki[ 38 * kd + 215 ] += -1 ;
    nuki[ 3 * kd + 215 ] += -1 ;
    nuki[ 6 * kd + 215 ] += +1 ;
    nuki[ 35 * kd + 215 ] += +1 ;

    /*reaction 217: CN + O <=> CO + N */
    nuki[ 39 * kd + 216 ] += -1 ;
    nuki[ 2 * kd + 216 ] += -1 ;
    nuki[ 14 * kd + 216 ] += +1 ;
    nuki[ 30 * kd + 216 ] += +1 ;

    /*reaction 218: CN + OH <=> NCO + H */
    nuki[ 39 * kd + 217 ] += -1 ;
    nuki[ 4 * kd + 217 ] += -1 ;
    nuki[ 46 * kd + 217 ] += +1 ;
    nuki[ 1 * kd + 217 ] += +1 ;

    /*reaction 219: CN + H2O <=> HCN + OH */
    nuki[ 39 * kd + 218 ] += -1 ;
    nuki[ 5 * kd + 218 ] += -1 ;
    nuki[ 40 * kd + 218 ] += +1 ;
    nuki[ 4 * kd + 218 ] += +1 ;

    /*reaction 220: CN + O2 <=> NCO + O */
    nuki[ 39 * kd + 219 ] += -1 ;
    nuki[ 3 * kd + 219 ] += -1 ;
    nuki[ 46 * kd + 219 ] += +1 ;
    nuki[ 2 * kd + 219 ] += +1 ;

    /*reaction 221: CN + H2 <=> HCN + H */
    nuki[ 39 * kd + 220 ] += -1 ;
    nuki[ 0 * kd + 220 ] += -1 ;
    nuki[ 40 * kd + 220 ] += +1 ;
    nuki[ 1 * kd + 220 ] += +1 ;

    /*reaction 222: NCO + O <=> NO + CO */
    nuki[ 46 * kd + 221 ] += -1 ;
    nuki[ 2 * kd + 221 ] += -1 ;
    nuki[ 35 * kd + 221 ] += +1 ;
    nuki[ 14 * kd + 221 ] += +1 ;

    /*reaction 223: NCO + H <=> NH + CO */
    nuki[ 46 * kd + 222 ] += -1 ;
    nuki[ 1 * kd + 222 ] += -1 ;
    nuki[ 31 * kd + 222 ] += +1 ;
    nuki[ 14 * kd + 222 ] += +1 ;

    /*reaction 224: NCO + OH <=> NO + H + CO */
    nuki[ 46 * kd + 223 ] += -1 ;
    nuki[ 4 * kd + 223 ] += -1 ;
    nuki[ 35 * kd + 223 ] += +1 ;
    nuki[ 1 * kd + 223 ] += +1 ;
    nuki[ 14 * kd + 223 ] += +1 ;

    /*reaction 225: NCO + N <=> N2 + CO */
    nuki[ 46 * kd + 224 ] += -1 ;
    nuki[ 30 * kd + 224 ] += -1 ;
    nuki[ 47 * kd + 224 ] += +1 ;
    nuki[ 14 * kd + 224 ] += +1 ;

    /*reaction 226: NCO + O2 <=> NO + CO2 */
    nuki[ 46 * kd + 225 ] += -1 ;
    nuki[ 3 * kd + 225 ] += -1 ;
    nuki[ 35 * kd + 225 ] += +1 ;
    nuki[ 15 * kd + 225 ] += +1 ;

    /*reaction 227: NCO + M <=> N + CO + M */
    nuki[ 46 * kd + 226 ] += -1 ;
    nuki[ 30 * kd + 226 ] += +1 ;
    nuki[ 14 * kd + 226 ] += +1 ;

    /*reaction 228: NCO + NO <=> N2O + CO */
    nuki[ 46 * kd + 227 ] += -1 ;
    nuki[ 35 * kd + 227 ] += -1 ;
    nuki[ 37 * kd + 227 ] += +1 ;
    nuki[ 14 * kd + 227 ] += +1 ;

    /*reaction 229: NCO + NO <=> N2 + CO2 */
    nuki[ 46 * kd + 228 ] += -1 ;
    nuki[ 35 * kd + 228 ] += -1 ;
    nuki[ 47 * kd + 228 ] += +1 ;
    nuki[ 15 * kd + 228 ] += +1 ;

    /*reaction 230: HCN + M <=> H + CN + M */
    nuki[ 40 * kd + 229 ] += -1 ;
    nuki[ 1 * kd + 229 ] += +1 ;
    nuki[ 39 * kd + 229 ] += +1 ;

    /*reaction 231: HCN + O <=> NCO + H */
    nuki[ 40 * kd + 230 ] += -1 ;
    nuki[ 2 * kd + 230 ] += -1 ;
    nuki[ 46 * kd + 230 ] += +1 ;
    nuki[ 1 * kd + 230 ] += +1 ;

    /*reaction 232: HCN + O <=> NH + CO */
    nuki[ 40 * kd + 231 ] += -1 ;
    nuki[ 2 * kd + 231 ] += -1 ;
    nuki[ 31 * kd + 231 ] += +1 ;
    nuki[ 14 * kd + 231 ] += +1 ;

    /*reaction 233: HCN + O <=> CN + OH */
    nuki[ 40 * kd + 232 ] += -1 ;
    nuki[ 2 * kd + 232 ] += -1 ;
    nuki[ 39 * kd + 232 ] += +1 ;
    nuki[ 4 * kd + 232 ] += +1 ;

    /*reaction 234: HCN + OH <=> HOCN + H */
    nuki[ 40 * kd + 233 ] += -1 ;
    nuki[ 4 * kd + 233 ] += -1 ;
    nuki[ 44 * kd + 233 ] += +1 ;
    nuki[ 1 * kd + 233 ] += +1 ;

    /*reaction 235: HCN + OH <=> HNCO + H */
    nuki[ 40 * kd + 234 ] += -1 ;
    nuki[ 4 * kd + 234 ] += -1 ;
    nuki[ 45 * kd + 234 ] += +1 ;
    nuki[ 1 * kd + 234 ] += +1 ;

    /*reaction 236: HCN + OH <=> NH2 + CO */
    nuki[ 40 * kd + 235 ] += -1 ;
    nuki[ 4 * kd + 235 ] += -1 ;
    nuki[ 32 * kd + 235 ] += +1 ;
    nuki[ 14 * kd + 235 ] += +1 ;

    /*reaction 237: H + HCN (+M) <=> H2CN (+M) */
    nuki[ 1 * kd + 236 ] += -1 ;
    nuki[ 40 * kd + 236 ] += -1 ;
    nuki[ 41 * kd + 236 ] += +1 ;

    /*reaction 238: H2CN + N <=> N2 + CH2 */
    nuki[ 41 * kd + 237 ] += -1 ;
    nuki[ 30 * kd + 237 ] += -1 ;
    nuki[ 47 * kd + 237 ] += +1 ;
    nuki[ 10 * kd + 237 ] += +1 ;

    /*reaction 239: C + N2 <=> CN + N */
    nuki[ 8 * kd + 238 ] += -1 ;
    nuki[ 47 * kd + 238 ] += -1 ;
    nuki[ 39 * kd + 238 ] += +1 ;
    nuki[ 30 * kd + 238 ] += +1 ;

    /*reaction 240: CH + N2 <=> HCN + N */
    nuki[ 9 * kd + 239 ] += -1 ;
    nuki[ 47 * kd + 239 ] += -1 ;
    nuki[ 40 * kd + 239 ] += +1 ;
    nuki[ 30 * kd + 239 ] += +1 ;

    /*reaction 241: CH + N2 (+M) <=> HCNN (+M) */
    nuki[ 9 * kd + 240 ] += -1 ;
    nuki[ 47 * kd + 240 ] += -1 ;
    nuki[ 42 * kd + 240 ] += +1 ;

    /*reaction 242: CH2 + N2 <=> HCN + NH */
    nuki[ 10 * kd + 241 ] += -1 ;
    nuki[ 47 * kd + 241 ] += -1 ;
    nuki[ 40 * kd + 241 ] += +1 ;
    nuki[ 31 * kd + 241 ] += +1 ;

    /*reaction 243: CH2(S) + N2 <=> NH + HCN */
    nuki[ 11 * kd + 242 ] += -1 ;
    nuki[ 47 * kd + 242 ] += -1 ;
    nuki[ 31 * kd + 242 ] += +1 ;
    nuki[ 40 * kd + 242 ] += +1 ;

    /*reaction 244: C + NO <=> CN + O */
    nuki[ 8 * kd + 243 ] += -1 ;
    nuki[ 35 * kd + 243 ] += -1 ;
    nuki[ 39 * kd + 243 ] += +1 ;
    nuki[ 2 * kd + 243 ] += +1 ;

    /*reaction 245: C + NO <=> CO + N */
    nuki[ 8 * kd + 244 ] += -1 ;
    nuki[ 35 * kd + 244 ] += -1 ;
    nuki[ 14 * kd + 244 ] += +1 ;
    nuki[ 30 * kd + 244 ] += +1 ;

    /*reaction 246: CH + NO <=> HCN + O */
    nuki[ 9 * kd + 245 ] += -1 ;
    nuki[ 35 * kd + 245 ] += -1 ;
    nuki[ 40 * kd + 245 ] += +1 ;
    nuki[ 2 * kd + 245 ] += +1 ;

    /*reaction 247: CH + NO <=> H + NCO */
    nuki[ 9 * kd + 246 ] += -1 ;
    nuki[ 35 * kd + 246 ] += -1 ;
    nuki[ 1 * kd + 246 ] += +1 ;
    nuki[ 46 * kd + 246 ] += +1 ;

    /*reaction 248: CH + NO <=> N + HCO */
    nuki[ 9 * kd + 247 ] += -1 ;
    nuki[ 35 * kd + 247 ] += -1 ;
    nuki[ 30 * kd + 247 ] += +1 ;
    nuki[ 16 * kd + 247 ] += +1 ;

    /*reaction 249: CH2 + NO <=> H + HNCO */
    nuki[ 10 * kd + 248 ] += -1 ;
    nuki[ 35 * kd + 248 ] += -1 ;
    nuki[ 1 * kd + 248 ] += +1 ;
    nuki[ 45 * kd + 248 ] += +1 ;

    /*reaction 250: CH2 + NO <=> OH + HCN */
    nuki[ 10 * kd + 249 ] += -1 ;
    nuki[ 35 * kd + 249 ] += -1 ;
    nuki[ 4 * kd + 249 ] += +1 ;
    nuki[ 40 * kd + 249 ] += +1 ;

    /*reaction 251: CH2 + NO <=> H + HCNO */
    nuki[ 10 * kd + 250 ] += -1 ;
    nuki[ 35 * kd + 250 ] += -1 ;
    nuki[ 1 * kd + 250 ] += +1 ;
    nuki[ 43 * kd + 250 ] += +1 ;

    /*reaction 252: CH2(S) + NO <=> H + HNCO */
    nuki[ 11 * kd + 251 ] += -1 ;
    nuki[ 35 * kd + 251 ] += -1 ;
    nuki[ 1 * kd + 251 ] += +1 ;
    nuki[ 45 * kd + 251 ] += +1 ;

    /*reaction 253: CH2(S) + NO <=> OH + HCN */
    nuki[ 11 * kd + 252 ] += -1 ;
    nuki[ 35 * kd + 252 ] += -1 ;
    nuki[ 4 * kd + 252 ] += +1 ;
    nuki[ 40 * kd + 252 ] += +1 ;

    /*reaction 254: CH2(S) + NO <=> H + HCNO */
    nuki[ 11 * kd + 253 ] += -1 ;
    nuki[ 35 * kd + 253 ] += -1 ;
    nuki[ 1 * kd + 253 ] += +1 ;
    nuki[ 43 * kd + 253 ] += +1 ;

    /*reaction 255: CH3 + NO <=> HCN + H2O */
    nuki[ 12 * kd + 254 ] += -1 ;
    nuki[ 35 * kd + 254 ] += -1 ;
    nuki[ 40 * kd + 254 ] += +1 ;
    nuki[ 5 * kd + 254 ] += +1 ;

    /*reaction 256: CH3 + NO <=> H2CN + OH */
    nuki[ 12 * kd + 255 ] += -1 ;
    nuki[ 35 * kd + 255 ] += -1 ;
    nuki[ 41 * kd + 255 ] += +1 ;
    nuki[ 4 * kd + 255 ] += +1 ;

    /*reaction 257: HCNN + O <=> CO + H + N2 */
    nuki[ 42 * kd + 256 ] += -1 ;
    nuki[ 2 * kd + 256 ] += -1 ;
    nuki[ 14 * kd + 256 ] += +1 ;
    nuki[ 1 * kd + 256 ] += +1 ;
    nuki[ 47 * kd + 256 ] += +1 ;

    /*reaction 258: HCNN + O <=> HCN + NO */
    nuki[ 42 * kd + 257 ] += -1 ;
    nuki[ 2 * kd + 257 ] += -1 ;
    nuki[ 40 * kd + 257 ] += +1 ;
    nuki[ 35 * kd + 257 ] += +1 ;

    /*reaction 259: HCNN + O2 <=> O + HCO + N2 */
    nuki[ 42 * kd + 258 ] += -1 ;
    nuki[ 3 * kd + 258 ] += -1 ;
    nuki[ 2 * kd + 258 ] += +1 ;
    nuki[ 16 * kd + 258 ] += +1 ;
    nuki[ 47 * kd + 258 ] += +1 ;

    /*reaction 260: HCNN + OH <=> H + HCO + N2 */
    nuki[ 42 * kd + 259 ] += -1 ;
    nuki[ 4 * kd + 259 ] += -1 ;
    nuki[ 1 * kd + 259 ] += +1 ;
    nuki[ 16 * kd + 259 ] += +1 ;
    nuki[ 47 * kd + 259 ] += +1 ;

    /*reaction 261: HCNN + H <=> CH2 + N2 */
    nuki[ 42 * kd + 260 ] += -1 ;
    nuki[ 1 * kd + 260 ] += -1 ;
    nuki[ 10 * kd + 260 ] += +1 ;
    nuki[ 47 * kd + 260 ] += +1 ;

    /*reaction 262: HNCO + O <=> NH + CO2 */
    nuki[ 45 * kd + 261 ] += -1 ;
    nuki[ 2 * kd + 261 ] += -1 ;
    nuki[ 31 * kd + 261 ] += +1 ;
    nuki[ 15 * kd + 261 ] += +1 ;

    /*reaction 263: HNCO + O <=> HNO + CO */
    nuki[ 45 * kd + 262 ] += -1 ;
    nuki[ 2 * kd + 262 ] += -1 ;
    nuki[ 38 * kd + 262 ] += +1 ;
    nuki[ 14 * kd + 262 ] += +1 ;

    /*reaction 264: HNCO + O <=> NCO + OH */
    nuki[ 45 * kd + 263 ] += -1 ;
    nuki[ 2 * kd + 263 ] += -1 ;
    nuki[ 46 * kd + 263 ] += +1 ;
    nuki[ 4 * kd + 263 ] += +1 ;

    /*reaction 265: HNCO + H <=> NH2 + CO */
    nuki[ 45 * kd + 264 ] += -1 ;
    nuki[ 1 * kd + 264 ] += -1 ;
    nuki[ 32 * kd + 264 ] += +1 ;
    nuki[ 14 * kd + 264 ] += +1 ;

    /*reaction 266: HNCO + H <=> H2 + NCO */
    nuki[ 45 * kd + 265 ] += -1 ;
    nuki[ 1 * kd + 265 ] += -1 ;
    nuki[ 0 * kd + 265 ] += +1 ;
    nuki[ 46 * kd + 265 ] += +1 ;

    /*reaction 267: HNCO + OH <=> NCO + H2O */
    nuki[ 45 * kd + 266 ] += -1 ;
    nuki[ 4 * kd + 266 ] += -1 ;
    nuki[ 46 * kd + 266 ] += +1 ;
    nuki[ 5 * kd + 266 ] += +1 ;

    /*reaction 268: HNCO + OH <=> NH2 + CO2 */
    nuki[ 45 * kd + 267 ] += -1 ;
    nuki[ 4 * kd + 267 ] += -1 ;
    nuki[ 32 * kd + 267 ] += +1 ;
    nuki[ 15 * kd + 267 ] += +1 ;

    /*reaction 269: HNCO + M <=> NH + CO + M */
    nuki[ 45 * kd + 268 ] += -1 ;
    nuki[ 31 * kd + 268 ] += +1 ;
    nuki[ 14 * kd + 268 ] += +1 ;

    /*reaction 270: HCNO + H <=> H + HNCO */
    nuki[ 43 * kd + 269 ] += -1 ;
    nuki[ 1 * kd + 269 ] += -1 ;
    nuki[ 1 * kd + 269 ] += +1 ;
    nuki[ 45 * kd + 269 ] += +1 ;

    /*reaction 271: HCNO + H <=> OH + HCN */
    nuki[ 43 * kd + 270 ] += -1 ;
    nuki[ 1 * kd + 270 ] += -1 ;
    nuki[ 4 * kd + 270 ] += +1 ;
    nuki[ 40 * kd + 270 ] += +1 ;

    /*reaction 272: HCNO + H <=> NH2 + CO */
    nuki[ 43 * kd + 271 ] += -1 ;
    nuki[ 1 * kd + 271 ] += -1 ;
    nuki[ 32 * kd + 271 ] += +1 ;
    nuki[ 14 * kd + 271 ] += +1 ;

    /*reaction 273: HOCN + H <=> H + HNCO */
    nuki[ 44 * kd + 272 ] += -1 ;
    nuki[ 1 * kd + 272 ] += -1 ;
    nuki[ 1 * kd + 272 ] += +1 ;
    nuki[ 45 * kd + 272 ] += +1 ;

    /*reaction 274: HCCO + NO <=> HCNO + CO */
    nuki[ 27 * kd + 273 ] += -1 ;
    nuki[ 35 * kd + 273 ] += -1 ;
    nuki[ 43 * kd + 273 ] += +1 ;
    nuki[ 14 * kd + 273 ] += +1 ;

    /*reaction 275: CH3 + N <=> H2CN + H */
    nuki[ 12 * kd + 274 ] += -1 ;
    nuki[ 30 * kd + 274 ] += -1 ;
    nuki[ 41 * kd + 274 ] += +1 ;
    nuki[ 1 * kd + 274 ] += +1 ;

    /*reaction 276: CH3 + N <=> HCN + H2 */
    nuki[ 12 * kd + 275 ] += -1 ;
    nuki[ 30 * kd + 275 ] += -1 ;
    nuki[ 40 * kd + 275 ] += +1 ;
    nuki[ 0 * kd + 275 ] += +1 ;

    /*reaction 277: NH3 + H <=> NH2 + H2 */
    nuki[ 33 * kd + 276 ] += -1 ;
    nuki[ 1 * kd + 276 ] += -1 ;
    nuki[ 32 * kd + 276 ] += +1 ;
    nuki[ 0 * kd + 276 ] += +1 ;

    /*reaction 278: NH3 + OH <=> NH2 + H2O */
    nuki[ 33 * kd + 277 ] += -1 ;
    nuki[ 4 * kd + 277 ] += -1 ;
    nuki[ 32 * kd + 277 ] += +1 ;
    nuki[ 5 * kd + 277 ] += +1 ;

    /*reaction 279: NH3 + O <=> NH2 + OH */
    nuki[ 33 * kd + 278 ] += -1 ;
    nuki[ 2 * kd + 278 ] += -1 ;
    nuki[ 32 * kd + 278 ] += +1 ;
    nuki[ 4 * kd + 278 ] += +1 ;

    /*reaction 280: NH + CO2 <=> HNO + CO */
    nuki[ 31 * kd + 279 ] += -1 ;
    nuki[ 15 * kd + 279 ] += -1 ;
    nuki[ 38 * kd + 279 ] += +1 ;
    nuki[ 14 * kd + 279 ] += +1 ;

    /*reaction 281: CN + NO2 <=> NCO + NO */
    nuki[ 39 * kd + 280 ] += -1 ;
    nuki[ 36 * kd + 280 ] += -1 ;
    nuki[ 46 * kd + 280 ] += +1 ;
    nuki[ 35 * kd + 280 ] += +1 ;

    /*reaction 282: NCO + NO2 <=> N2O + CO2 */
    nuki[ 46 * kd + 281 ] += -1 ;
    nuki[ 36 * kd + 281 ] += -1 ;
    nuki[ 37 * kd + 281 ] += +1 ;
    nuki[ 15 * kd + 281 ] += +1 ;

    /*reaction 283: N + CO2 <=> NO + CO */
    nuki[ 30 * kd + 282 ] += -1 ;
    nuki[ 15 * kd + 282 ] += -1 ;
    nuki[ 35 * kd + 282 ] += +1 ;
    nuki[ 14 * kd + 282 ] += +1 ;

    /*reaction 284: O + CH3 => H + H2 + CO */
    nuki[ 2 * kd + 283 ] += -1 ;
    nuki[ 12 * kd + 283 ] += -1 ;
    nuki[ 1 * kd + 283 ] += +1 ;
    nuki[ 0 * kd + 283 ] += +1 ;
    nuki[ 14 * kd + 283 ] += +1 ;

    /*reaction 285: O + C2H4 <=> H + CH2CHO */
    nuki[ 2 * kd + 284 ] += -1 ;
    nuki[ 24 * kd + 284 ] += -1 ;
    nuki[ 1 * kd + 284 ] += +1 ;
    nuki[ 51 * kd + 284 ] += +1 ;

    /*reaction 286: O + C2H5 <=> H + CH3CHO */
    nuki[ 2 * kd + 285 ] += -1 ;
    nuki[ 25 * kd + 285 ] += -1 ;
    nuki[ 1 * kd + 285 ] += +1 ;
    nuki[ 52 * kd + 285 ] += +1 ;

    /*reaction 287: OH + HO2 <=> O2 + H2O */
    nuki[ 4 * kd + 286 ] += -1 ;
    nuki[ 6 * kd + 286 ] += -1 ;
    nuki[ 3 * kd + 286 ] += +1 ;
    nuki[ 5 * kd + 286 ] += +1 ;

    /*reaction 288: OH + CH3 => H2 + CH2O */
    nuki[ 4 * kd + 287 ] += -1 ;
    nuki[ 12 * kd + 287 ] += -1 ;
    nuki[ 0 * kd + 287 ] += +1 ;
    nuki[ 17 * kd + 287 ] += +1 ;

    /*reaction 289: CH + H2 (+M) <=> CH3 (+M) */
    nuki[ 9 * kd + 288 ] += -1 ;
    nuki[ 0 * kd + 288 ] += -1 ;
    nuki[ 12 * kd + 288 ] += +1 ;

    /*reaction 290: CH2 + O2 => 2 H + CO2 */
    nuki[ 10 * kd + 289 ] += -1 ;
    nuki[ 3 * kd + 289 ] += -1 ;
    nuki[ 1 * kd + 289 ] += +2 ;
    nuki[ 15 * kd + 289 ] += +1 ;

    /*reaction 291: CH2 + O2 <=> O + CH2O */
    nuki[ 10 * kd + 290 ] += -1 ;
    nuki[ 3 * kd + 290 ] += -1 ;
    nuki[ 2 * kd + 290 ] += +1 ;
    nuki[ 17 * kd + 290 ] += +1 ;

    /*reaction 292: CH2 + CH2 => 2 H + C2H2 */
    nuki[ 10 * kd + 291 ] += -1 ;
    nuki[ 10 * kd + 291 ] += -1 ;
    nuki[ 1 * kd + 291 ] += +2 ;
    nuki[ 22 * kd + 291 ] += +1 ;

    /*reaction 293: CH2(S) + H2O => H2 + CH2O */
    nuki[ 11 * kd + 292 ] += -1 ;
    nuki[ 5 * kd + 292 ] += -1 ;
    nuki[ 0 * kd + 292 ] += +1 ;
    nuki[ 17 * kd + 292 ] += +1 ;

    /*reaction 294: C2H3 + O2 <=> O + CH2CHO */
    nuki[ 23 * kd + 293 ] += -1 ;
    nuki[ 3 * kd + 293 ] += -1 ;
    nuki[ 2 * kd + 293 ] += +1 ;
    nuki[ 51 * kd + 293 ] += +1 ;

    /*reaction 295: C2H3 + O2 <=> HO2 + C2H2 */
    nuki[ 23 * kd + 294 ] += -1 ;
    nuki[ 3 * kd + 294 ] += -1 ;
    nuki[ 6 * kd + 294 ] += +1 ;
    nuki[ 22 * kd + 294 ] += +1 ;

    /*reaction 296: O + CH3CHO <=> OH + CH2CHO */
    nuki[ 2 * kd + 295 ] += -1 ;
    nuki[ 52 * kd + 295 ] += -1 ;
    nuki[ 4 * kd + 295 ] += +1 ;
    nuki[ 51 * kd + 295 ] += +1 ;

    /*reaction 297: O + CH3CHO => OH + CH3 + CO */
    nuki[ 2 * kd + 296 ] += -1 ;
    nuki[ 52 * kd + 296 ] += -1 ;
    nuki[ 4 * kd + 296 ] += +1 ;
    nuki[ 12 * kd + 296 ] += +1 ;
    nuki[ 14 * kd + 296 ] += +1 ;

    /*reaction 298: O2 + CH3CHO => HO2 + CH3 + CO */
    nuki[ 3 * kd + 297 ] += -1 ;
    nuki[ 52 * kd + 297 ] += -1 ;
    nuki[ 6 * kd + 297 ] += +1 ;
    nuki[ 12 * kd + 297 ] += +1 ;
    nuki[ 14 * kd + 297 ] += +1 ;

    /*reaction 299: H + CH3CHO <=> CH2CHO + H2 */
    nuki[ 1 * kd + 298 ] += -1 ;
    nuki[ 52 * kd + 298 ] += -1 ;
    nuki[ 51 * kd + 298 ] += +1 ;
    nuki[ 0 * kd + 298 ] += +1 ;

    /*reaction 300: H + CH3CHO => CH3 + H2 + CO */
    nuki[ 1 * kd + 299 ] += -1 ;
    nuki[ 52 * kd + 299 ] += -1 ;
    nuki[ 12 * kd + 299 ] += +1 ;
    nuki[ 0 * kd + 299 ] += +1 ;
    nuki[ 14 * kd + 299 ] += +1 ;

    /*reaction 301: OH + CH3CHO => CH3 + H2O + CO */
    nuki[ 4 * kd + 300 ] += -1 ;
    nuki[ 52 * kd + 300 ] += -1 ;
    nuki[ 12 * kd + 300 ] += +1 ;
    nuki[ 5 * kd + 300 ] += +1 ;
    nuki[ 14 * kd + 300 ] += +1 ;

    /*reaction 302: HO2 + CH3CHO => CH3 + H2O2 + CO */
    nuki[ 6 * kd + 301 ] += -1 ;
    nuki[ 52 * kd + 301 ] += -1 ;
    nuki[ 12 * kd + 301 ] += +1 ;
    nuki[ 7 * kd + 301 ] += +1 ;
    nuki[ 14 * kd + 301 ] += +1 ;

    /*reaction 303: CH3 + CH3CHO => CH3 + CH4 + CO */
    nuki[ 12 * kd + 302 ] += -1 ;
    nuki[ 52 * kd + 302 ] += -1 ;
    nuki[ 12 * kd + 302 ] += +1 ;
    nuki[ 13 * kd + 302 ] += +1 ;
    nuki[ 14 * kd + 302 ] += +1 ;

    /*reaction 304: H + CH2CO (+M) <=> CH2CHO (+M) */
    nuki[ 1 * kd + 303 ] += -1 ;
    nuki[ 28 * kd + 303 ] += -1 ;
    nuki[ 51 * kd + 303 ] += +1 ;

    /*reaction 305: O + CH2CHO => H + CH2 + CO2 */
    nuki[ 2 * kd + 304 ] += -1 ;
    nuki[ 51 * kd + 304 ] += -1 ;
    nuki[ 1 * kd + 304 ] += +1 ;
    nuki[ 10 * kd + 304 ] += +1 ;
    nuki[ 15 * kd + 304 ] += +1 ;

    /*reaction 306: O2 + CH2CHO => OH + CO + CH2O */
    nuki[ 3 * kd + 305 ] += -1 ;
    nuki[ 51 * kd + 305 ] += -1 ;
    nuki[ 4 * kd + 305 ] += +1 ;
    nuki[ 14 * kd + 305 ] += +1 ;
    nuki[ 17 * kd + 305 ] += +1 ;

    /*reaction 307: O2 + CH2CHO => OH + 2 HCO */
    nuki[ 3 * kd + 306 ] += -1 ;
    nuki[ 51 * kd + 306 ] += -1 ;
    nuki[ 4 * kd + 306 ] += +1 ;
    nuki[ 16 * kd + 306 ] += +2 ;

    /*reaction 308: H + CH2CHO <=> CH3 + HCO */
    nuki[ 1 * kd + 307 ] += -1 ;
    nuki[ 51 * kd + 307 ] += -1 ;
    nuki[ 12 * kd + 307 ] += +1 ;
    nuki[ 16 * kd + 307 ] += +1 ;

    /*reaction 309: H + CH2CHO <=> CH2CO + H2 */
    nuki[ 1 * kd + 308 ] += -1 ;
    nuki[ 51 * kd + 308 ] += -1 ;
    nuki[ 28 * kd + 308 ] += +1 ;
    nuki[ 0 * kd + 308 ] += +1 ;

    /*reaction 310: OH + CH2CHO <=> H2O + CH2CO */
    nuki[ 4 * kd + 309 ] += -1 ;
    nuki[ 51 * kd + 309 ] += -1 ;
    nuki[ 5 * kd + 309 ] += +1 ;
    nuki[ 28 * kd + 309 ] += +1 ;

    /*reaction 311: OH + CH2CHO <=> HCO + CH2OH */
    nuki[ 4 * kd + 310 ] += -1 ;
    nuki[ 51 * kd + 310 ] += -1 ;
    nuki[ 16 * kd + 310 ] += +1 ;
    nuki[ 18 * kd + 310 ] += +1 ;

    /*reaction 312: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    nuki[ 12 * kd + 311 ] += -1 ;
    nuki[ 25 * kd + 311 ] += -1 ;
    nuki[ 50 * kd + 311 ] += +1 ;

    /*reaction 313: O + C3H8 <=> OH + C3H7 */
    nuki[ 2 * kd + 312 ] += -1 ;
    nuki[ 50 * kd + 312 ] += -1 ;
    nuki[ 4 * kd + 312 ] += +1 ;
    nuki[ 49 * kd + 312 ] += +1 ;

    /*reaction 314: H + C3H8 <=> C3H7 + H2 */
    nuki[ 1 * kd + 313 ] += -1 ;
    nuki[ 50 * kd + 313 ] += -1 ;
    nuki[ 49 * kd + 313 ] += +1 ;
    nuki[ 0 * kd + 313 ] += +1 ;

    /*reaction 315: OH + C3H8 <=> C3H7 + H2O */
    nuki[ 4 * kd + 314 ] += -1 ;
    nuki[ 50 * kd + 314 ] += -1 ;
    nuki[ 49 * kd + 314 ] += +1 ;
    nuki[ 5 * kd + 314 ] += +1 ;

    /*reaction 316: C3H7 + H2O2 <=> HO2 + C3H8 */
    nuki[ 49 * kd + 315 ] += -1 ;
    nuki[ 7 * kd + 315 ] += -1 ;
    nuki[ 6 * kd + 315 ] += +1 ;
    nuki[ 50 * kd + 315 ] += +1 ;

    /*reaction 317: CH3 + C3H8 <=> C3H7 + CH4 */
    nuki[ 12 * kd + 316 ] += -1 ;
    nuki[ 50 * kd + 316 ] += -1 ;
    nuki[ 49 * kd + 316 ] += +1 ;
    nuki[ 13 * kd + 316 ] += +1 ;

    /*reaction 318: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    nuki[ 12 * kd + 317 ] += -1 ;
    nuki[ 24 * kd + 317 ] += -1 ;
    nuki[ 49 * kd + 317 ] += +1 ;

    /*reaction 319: O + C3H7 <=> C2H5 + CH2O */
    nuki[ 2 * kd + 318 ] += -1 ;
    nuki[ 49 * kd + 318 ] += -1 ;
    nuki[ 25 * kd + 318 ] += +1 ;
    nuki[ 17 * kd + 318 ] += +1 ;

    /*reaction 320: H + C3H7 (+M) <=> C3H8 (+M) */
    nuki[ 1 * kd + 319 ] += -1 ;
    nuki[ 49 * kd + 319 ] += -1 ;
    nuki[ 50 * kd + 319 ] += +1 ;

    /*reaction 321: H + C3H7 <=> CH3 + C2H5 */
    nuki[ 1 * kd + 320 ] += -1 ;
    nuki[ 49 * kd + 320 ] += -1 ;
    nuki[ 12 * kd + 320 ] += +1 ;
    nuki[ 25 * kd + 320 ] += +1 ;

    /*reaction 322: OH + C3H7 <=> C2H5 + CH2OH */
    nuki[ 4 * kd + 321 ] += -1 ;
    nuki[ 49 * kd + 321 ] += -1 ;
    nuki[ 25 * kd + 321 ] += +1 ;
    nuki[ 18 * kd + 321 ] += +1 ;

    /*reaction 323: HO2 + C3H7 <=> O2 + C3H8 */
    nuki[ 6 * kd + 322 ] += -1 ;
    nuki[ 49 * kd + 322 ] += -1 ;
    nuki[ 3 * kd + 322 ] += +1 ;
    nuki[ 50 * kd + 322 ] += +1 ;

    /*reaction 324: HO2 + C3H7 => OH + C2H5 + CH2O */
    nuki[ 6 * kd + 323 ] += -1 ;
    nuki[ 49 * kd + 323 ] += -1 ;
    nuki[ 4 * kd + 323 ] += +1 ;
    nuki[ 25 * kd + 323 ] += +1 ;
    nuki[ 17 * kd + 323 ] += +1 ;

    /*reaction 325: CH3 + C3H7 <=> 2 C2H5 */
    nuki[ 12 * kd + 324 ] += -1 ;
    nuki[ 49 * kd + 324 ] += -1 ;
    nuki[ 25 * kd + 324 ] += +2 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < 5 * 53; ++ id) {
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

    /*N */
    ncf[ 30 * kd + 3 ] = 1; /*N */

    /*NH */
    ncf[ 31 * kd + 3 ] = 1; /*N */
    ncf[ 31 * kd + 1 ] = 1; /*H */

    /*NH2 */
    ncf[ 32 * kd + 3 ] = 1; /*N */
    ncf[ 32 * kd + 1 ] = 2; /*H */

    /*NH3 */
    ncf[ 33 * kd + 3 ] = 1; /*N */
    ncf[ 33 * kd + 1 ] = 3; /*H */

    /*NNH */
    ncf[ 34 * kd + 3 ] = 2; /*N */
    ncf[ 34 * kd + 1 ] = 1; /*H */

    /*NO */
    ncf[ 35 * kd + 3 ] = 1; /*N */
    ncf[ 35 * kd + 0 ] = 1; /*O */

    /*NO2 */
    ncf[ 36 * kd + 3 ] = 1; /*N */
    ncf[ 36 * kd + 0 ] = 2; /*O */

    /*N2O */
    ncf[ 37 * kd + 3 ] = 2; /*N */
    ncf[ 37 * kd + 0 ] = 1; /*O */

    /*HNO */
    ncf[ 38 * kd + 1 ] = 1; /*H */
    ncf[ 38 * kd + 3 ] = 1; /*N */
    ncf[ 38 * kd + 0 ] = 1; /*O */

    /*CN */
    ncf[ 39 * kd + 2 ] = 1; /*C */
    ncf[ 39 * kd + 3 ] = 1; /*N */

    /*HCN */
    ncf[ 40 * kd + 1 ] = 1; /*H */
    ncf[ 40 * kd + 2 ] = 1; /*C */
    ncf[ 40 * kd + 3 ] = 1; /*N */

    /*H2CN */
    ncf[ 41 * kd + 1 ] = 2; /*H */
    ncf[ 41 * kd + 2 ] = 1; /*C */
    ncf[ 41 * kd + 3 ] = 1; /*N */

    /*HCNN */
    ncf[ 42 * kd + 2 ] = 1; /*C */
    ncf[ 42 * kd + 3 ] = 2; /*N */
    ncf[ 42 * kd + 1 ] = 1; /*H */

    /*HCNO */
    ncf[ 43 * kd + 1 ] = 1; /*H */
    ncf[ 43 * kd + 3 ] = 1; /*N */
    ncf[ 43 * kd + 2 ] = 1; /*C */
    ncf[ 43 * kd + 0 ] = 1; /*O */

    /*HOCN */
    ncf[ 44 * kd + 1 ] = 1; /*H */
    ncf[ 44 * kd + 3 ] = 1; /*N */
    ncf[ 44 * kd + 2 ] = 1; /*C */
    ncf[ 44 * kd + 0 ] = 1; /*O */

    /*HNCO */
    ncf[ 45 * kd + 1 ] = 1; /*H */
    ncf[ 45 * kd + 3 ] = 1; /*N */
    ncf[ 45 * kd + 2 ] = 1; /*C */
    ncf[ 45 * kd + 0 ] = 1; /*O */

    /*NCO */
    ncf[ 46 * kd + 3 ] = 1; /*N */
    ncf[ 46 * kd + 2 ] = 1; /*C */
    ncf[ 46 * kd + 0 ] = 1; /*O */

    /*N2 */
    ncf[ 47 * kd + 3 ] = 2; /*N */

    /*AR */
    ncf[ 48 * kd + 4 ] = 1; /*AR */

    /*C3H7 */
    ncf[ 49 * kd + 2 ] = 3; /*C */
    ncf[ 49 * kd + 1 ] = 7; /*H */

    /*C3H8 */
    ncf[ 50 * kd + 2 ] = 3; /*C */
    ncf[ 50 * kd + 1 ] = 8; /*H */

    /*CH2CHO */
    ncf[ 51 * kd + 0 ] = 1; /*O */
    ncf[ 51 * kd + 1 ] = 3; /*H */
    ncf[ 51 * kd + 2 ] = 2; /*C */

    /*CH3CHO */
    ncf[ 52 * kd + 2 ] = 2; /*C */
    ncf[ 52 * kd + 1 ] = 4; /*H */
    ncf[ 52 * kd + 0 ] = 1; /*O */

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
    a[2] = 38700;
    b[2] = 2.7;
    e[2] = 6260;

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
    a[9] = 5.06e+13;
    b[9] = 0;
    e[9] = 0;

    /*reaction 11: O + CH4 <=> OH + CH3 */
    a[10] = 1.02e+09;
    b[10] = 1.5;
    e[10] = 8600;

    /*reaction 12: O + CO (+M) <=> CO2 (+M) */
    a[11] = 1.8e+10;
    b[11] = 0;
    e[11] = 2385;

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
    a[20] = 1.35e+07;
    b[20] = 2;
    e[20] = 1900;

    /*reaction 22: O + C2H2 <=> OH + C2H */
    a[21] = 4.6e+19;
    b[21] = -1.41;
    e[21] = 28950;

    /*reaction 23: O + C2H2 <=> CO + CH2 */
    a[22] = 6.94e+06;
    b[22] = 2;
    e[22] = 1900;

    /*reaction 24: O + C2H3 <=> H + CH2CO */
    a[23] = 3e+13;
    b[23] = 0;
    e[23] = 0;

    /*reaction 25: O + C2H4 <=> CH3 + HCO */
    a[24] = 1.25e+07;
    b[24] = 1.83;
    e[24] = 220;

    /*reaction 26: O + C2H5 <=> CH3 + CH2O */
    a[25] = 2.24e+13;
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
    a[33] = 2.08e+19;
    b[33] = -1.24;
    e[33] = 0;

    /*reaction 35: H + O2 + H2O <=> HO2 + H2O */
    a[34] = 1.126e+19;
    b[34] = -0.76;
    e[34] = 0;

    /*reaction 36: H + O2 + N2 <=> HO2 + N2 */
    a[35] = 2.6e+19;
    b[35] = -1.24;
    e[35] = 0;

    /*reaction 37: H + O2 + AR <=> HO2 + AR */
    a[36] = 7e+17;
    b[36] = -0.8;
    e[36] = 0;

    /*reaction 38: H + O2 <=> O + OH */
    a[37] = 2.65e+16;
    b[37] = -0.6707;
    e[37] = 17041;

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
    a[44] = 4.48e+13;
    b[44] = 0;
    e[44] = 1068;

    /*reaction 46: H + HO2 <=> 2 OH */
    a[45] = 8.4e+13;
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
    a[48] = 1.65e+14;
    b[48] = 0;
    e[48] = 0;

    /*reaction 50: H + CH2 (+M) <=> CH3 (+M) */
    a[49] = 6e+14;
    b[49] = 0;
    e[49] = 0;

    /*reaction 51: H + CH2(S) <=> CH + H2 */
    a[50] = 3e+13;
    b[50] = 0;
    e[50] = 0;

    /*reaction 52: H + CH3 (+M) <=> CH4 (+M) */
    a[51] = 1.39e+16;
    b[51] = -0.534;
    e[51] = 536;

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
    a[57] = 5.74e+07;
    b[57] = 1.9;
    e[57] = 2742;

    /*reaction 59: H + CH2OH (+M) <=> CH3OH (+M) */
    a[58] = 1.055e+12;
    b[58] = 0.5;
    e[58] = 86;

    /*reaction 60: H + CH2OH <=> H2 + CH2O */
    a[59] = 2e+13;
    b[59] = 0;
    e[59] = 0;

    /*reaction 61: H + CH2OH <=> OH + CH3 */
    a[60] = 1.65e+11;
    b[60] = 0.65;
    e[60] = -284;

    /*reaction 62: H + CH2OH <=> CH2(S) + H2O */
    a[61] = 3.28e+13;
    b[61] = -0.09;
    e[61] = 610;

    /*reaction 63: H + CH3O (+M) <=> CH3OH (+M) */
    a[62] = 2.43e+12;
    b[62] = 0.515;
    e[62] = 50;

    /*reaction 64: H + CH3O <=> H + CH2OH */
    a[63] = 4.15e+07;
    b[63] = 1.63;
    e[63] = 1924;

    /*reaction 65: H + CH3O <=> H2 + CH2O */
    a[64] = 2e+13;
    b[64] = 0;
    e[64] = 0;

    /*reaction 66: H + CH3O <=> OH + CH3 */
    a[65] = 1.5e+12;
    b[65] = 0.5;
    e[65] = -110;

    /*reaction 67: H + CH3O <=> CH2(S) + H2O */
    a[66] = 2.62e+14;
    b[66] = -0.23;
    e[66] = 1070;

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
    a[73] = 5.4e+11;
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
    a[86] = 1.45e+13;
    b[86] = 0;
    e[86] = -500;

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    a[87] = 2e+12;
    b[87] = 0;
    e[87] = 427;

    /*reaction 89: OH + H2O2 <=> HO2 + H2O */
    a[88] = 1.7e+18;
    b[88] = 0;
    e[88] = 29410;

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
    a[94] = 2.79e+18;
    b[94] = -1.43;
    e[94] = 1330;

    /*reaction 96: OH + CH3 <=> CH2 + H2O */
    a[95] = 5.6e+07;
    b[95] = 1.6;
    e[95] = 5420;

    /*reaction 97: OH + CH3 <=> CH2(S) + H2O */
    a[96] = 6.44e+17;
    b[96] = -1.34;
    e[96] = 1417;

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
    a[118] = 3.78e+13;
    b[118] = 0;
    e[118] = 0;

    /*reaction 120: HO2 + CO <=> OH + CO2 */
    a[119] = 1.5e+14;
    b[119] = 0;
    e[119] = 23600;

    /*reaction 121: HO2 + CH2O <=> HCO + H2O2 */
    a[120] = 5.6e+06;
    b[120] = 2;
    e[120] = 12000;

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
    a[124] = 6.71e+13;
    b[124] = 0;
    e[124] = 0;

    /*reaction 126: CH + H2 <=> H + CH2 */
    a[125] = 1.08e+14;
    b[125] = 0;
    e[125] = 3110;

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
    a[131] = 1.9e+14;
    b[131] = 0;
    e[131] = 15792;

    /*reaction 133: CH + CH2O <=> H + CH2CO */
    a[132] = 9.46e+13;
    b[132] = 0;
    e[132] = -515;

    /*reaction 134: CH + HCCO <=> CO + C2H2 */
    a[133] = 5e+13;
    b[133] = 0;
    e[133] = 0;

    /*reaction 135: CH2 + O2 => OH + H + CO */
    a[134] = 5e+12;
    b[134] = 0;
    e[134] = 1500;

    /*reaction 136: CH2 + H2 <=> H + CH3 */
    a[135] = 500000;
    b[135] = 2;
    e[135] = 7230;

    /*reaction 137: 2 CH2 <=> H2 + C2H2 */
    a[136] = 1.6e+15;
    b[136] = 0;
    e[136] = 11944;

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
    a[146] = 4.82e+17;
    b[146] = -1.16;
    e[146] = 1145;

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
    a[154] = 3.56e+13;
    b[154] = 0;
    e[154] = 30480;

    /*reaction 156: CH3 + O2 <=> OH + CH2O */
    a[155] = 2.31e+12;
    b[155] = 0;
    e[155] = 20315;

    /*reaction 157: CH3 + H2O2 <=> HO2 + CH4 */
    a[156] = 24500;
    b[156] = 2.47;
    e[156] = 5180;

    /*reaction 158: 2 CH3 (+M) <=> C2H6 (+M) */
    a[157] = 6.77e+16;
    b[157] = -1.18;
    e[157] = 654;

    /*reaction 159: 2 CH3 <=> H + C2H5 */
    a[158] = 6.84e+12;
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
    a[165] = 1.5e+18;
    b[165] = -1;
    e[165] = 17000;

    /*reaction 167: HCO + M <=> H + CO + M */
    a[166] = 1.87e+17;
    b[166] = -1;
    e[166] = 17000;

    /*reaction 168: HCO + O2 <=> HO2 + CO */
    a[167] = 1.345e+13;
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
    a[170] = 1e+13;
    b[170] = 0;
    e[170] = -755;

    /*reaction 172: C2H + H2 <=> H + C2H2 */
    a[171] = 5.68e+10;
    b[171] = 0.9;
    e[171] = 1993;

    /*reaction 173: C2H3 + O2 <=> HCO + CH2O */
    a[172] = 4.58e+16;
    b[172] = -1.39;
    e[172] = 1015;

    /*reaction 174: C2H4 (+M) <=> H2 + C2H2 (+M) */
    a[173] = 8e+12;
    b[173] = 0.44;
    e[173] = 86770;

    /*reaction 175: C2H5 + O2 <=> HO2 + C2H4 */
    a[174] = 8.4e+11;
    b[174] = 0;
    e[174] = 3875;

    /*reaction 176: HCCO + O2 <=> OH + 2 CO */
    a[175] = 3.2e+12;
    b[175] = 0;
    e[175] = 854;

    /*reaction 177: 2 HCCO <=> 2 CO + C2H2 */
    a[176] = 1e+13;
    b[176] = 0;
    e[176] = 0;

    /*reaction 178: N + NO <=> N2 + O */
    a[177] = 2.7e+13;
    b[177] = 0;
    e[177] = 355;

    /*reaction 179: N + O2 <=> NO + O */
    a[178] = 9e+09;
    b[178] = 1;
    e[178] = 6500;

    /*reaction 180: N + OH <=> NO + H */
    a[179] = 3.36e+13;
    b[179] = 0;
    e[179] = 385;

    /*reaction 181: N2O + O <=> N2 + O2 */
    a[180] = 1.4e+12;
    b[180] = 0;
    e[180] = 10810;

    /*reaction 182: N2O + O <=> 2 NO */
    a[181] = 2.9e+13;
    b[181] = 0;
    e[181] = 23150;

    /*reaction 183: N2O + H <=> N2 + OH */
    a[182] = 3.87e+14;
    b[182] = 0;
    e[182] = 18880;

    /*reaction 184: N2O + OH <=> N2 + HO2 */
    a[183] = 2e+12;
    b[183] = 0;
    e[183] = 21060;

    /*reaction 185: N2O (+M) <=> N2 + O (+M) */
    a[184] = 7.91e+10;
    b[184] = 0;
    e[184] = 56020;

    /*reaction 186: HO2 + NO <=> NO2 + OH */
    a[185] = 2.11e+12;
    b[185] = 0;
    e[185] = -480;

    /*reaction 187: NO + O + M <=> NO2 + M */
    a[186] = 1.06e+20;
    b[186] = -1.41;
    e[186] = 0;

    /*reaction 188: NO2 + O <=> NO + O2 */
    a[187] = 3.9e+12;
    b[187] = 0;
    e[187] = -240;

    /*reaction 189: NO2 + H <=> NO + OH */
    a[188] = 1.32e+14;
    b[188] = 0;
    e[188] = 360;

    /*reaction 190: NH + O <=> NO + H */
    a[189] = 4e+13;
    b[189] = 0;
    e[189] = 0;

    /*reaction 191: NH + H <=> N + H2 */
    a[190] = 3.2e+13;
    b[190] = 0;
    e[190] = 330;

    /*reaction 192: NH + OH <=> HNO + H */
    a[191] = 2e+13;
    b[191] = 0;
    e[191] = 0;

    /*reaction 193: NH + OH <=> N + H2O */
    a[192] = 2e+09;
    b[192] = 1.2;
    e[192] = 0;

    /*reaction 194: NH + O2 <=> HNO + O */
    a[193] = 461000;
    b[193] = 2;
    e[193] = 6500;

    /*reaction 195: NH + O2 <=> NO + OH */
    a[194] = 1.28e+06;
    b[194] = 1.5;
    e[194] = 100;

    /*reaction 196: NH + N <=> N2 + H */
    a[195] = 1.5e+13;
    b[195] = 0;
    e[195] = 0;

    /*reaction 197: NH + H2O <=> HNO + H2 */
    a[196] = 2e+13;
    b[196] = 0;
    e[196] = 13850;

    /*reaction 198: NH + NO <=> N2 + OH */
    a[197] = 2.16e+13;
    b[197] = -0.23;
    e[197] = 0;

    /*reaction 199: NH + NO <=> N2O + H */
    a[198] = 3.65e+14;
    b[198] = -0.45;
    e[198] = 0;

    /*reaction 200: NH2 + O <=> OH + NH */
    a[199] = 3e+12;
    b[199] = 0;
    e[199] = 0;

    /*reaction 201: NH2 + O <=> H + HNO */
    a[200] = 3.9e+13;
    b[200] = 0;
    e[200] = 0;

    /*reaction 202: NH2 + H <=> NH + H2 */
    a[201] = 4e+13;
    b[201] = 0;
    e[201] = 3650;

    /*reaction 203: NH2 + OH <=> NH + H2O */
    a[202] = 9e+07;
    b[202] = 1.5;
    e[202] = -460;

    /*reaction 204: NNH <=> N2 + H */
    a[203] = 3.3e+08;
    b[203] = 0;
    e[203] = 0;

    /*reaction 205: NNH + M <=> N2 + H + M */
    a[204] = 1.3e+14;
    b[204] = -0.11;
    e[204] = 4980;

    /*reaction 206: NNH + O2 <=> HO2 + N2 */
    a[205] = 5e+12;
    b[205] = 0;
    e[205] = 0;

    /*reaction 207: NNH + O <=> OH + N2 */
    a[206] = 2.5e+13;
    b[206] = 0;
    e[206] = 0;

    /*reaction 208: NNH + O <=> NH + NO */
    a[207] = 7e+13;
    b[207] = 0;
    e[207] = 0;

    /*reaction 209: NNH + H <=> H2 + N2 */
    a[208] = 5e+13;
    b[208] = 0;
    e[208] = 0;

    /*reaction 210: NNH + OH <=> H2O + N2 */
    a[209] = 2e+13;
    b[209] = 0;
    e[209] = 0;

    /*reaction 211: NNH + CH3 <=> CH4 + N2 */
    a[210] = 2.5e+13;
    b[210] = 0;
    e[210] = 0;

    /*reaction 212: H + NO + M <=> HNO + M */
    a[211] = 4.48e+19;
    b[211] = -1.32;
    e[211] = 740;

    /*reaction 213: HNO + O <=> NO + OH */
    a[212] = 2.5e+13;
    b[212] = 0;
    e[212] = 0;

    /*reaction 214: HNO + H <=> H2 + NO */
    a[213] = 9e+11;
    b[213] = 0.72;
    e[213] = 660;

    /*reaction 215: HNO + OH <=> NO + H2O */
    a[214] = 1.3e+07;
    b[214] = 1.9;
    e[214] = -950;

    /*reaction 216: HNO + O2 <=> HO2 + NO */
    a[215] = 1e+13;
    b[215] = 0;
    e[215] = 13000;

    /*reaction 217: CN + O <=> CO + N */
    a[216] = 7.7e+13;
    b[216] = 0;
    e[216] = 0;

    /*reaction 218: CN + OH <=> NCO + H */
    a[217] = 4e+13;
    b[217] = 0;
    e[217] = 0;

    /*reaction 219: CN + H2O <=> HCN + OH */
    a[218] = 8e+12;
    b[218] = 0;
    e[218] = 7460;

    /*reaction 220: CN + O2 <=> NCO + O */
    a[219] = 6.14e+12;
    b[219] = 0;
    e[219] = -440;

    /*reaction 221: CN + H2 <=> HCN + H */
    a[220] = 295000;
    b[220] = 2.45;
    e[220] = 2240;

    /*reaction 222: NCO + O <=> NO + CO */
    a[221] = 2.35e+13;
    b[221] = 0;
    e[221] = 0;

    /*reaction 223: NCO + H <=> NH + CO */
    a[222] = 5.4e+13;
    b[222] = 0;
    e[222] = 0;

    /*reaction 224: NCO + OH <=> NO + H + CO */
    a[223] = 2.5e+12;
    b[223] = 0;
    e[223] = 0;

    /*reaction 225: NCO + N <=> N2 + CO */
    a[224] = 2e+13;
    b[224] = 0;
    e[224] = 0;

    /*reaction 226: NCO + O2 <=> NO + CO2 */
    a[225] = 2e+12;
    b[225] = 0;
    e[225] = 20000;

    /*reaction 227: NCO + M <=> N + CO + M */
    a[226] = 3.1e+14;
    b[226] = 0;
    e[226] = 54050;

    /*reaction 228: NCO + NO <=> N2O + CO */
    a[227] = 1.9e+17;
    b[227] = -1.52;
    e[227] = 740;

    /*reaction 229: NCO + NO <=> N2 + CO2 */
    a[228] = 3.8e+18;
    b[228] = -2;
    e[228] = 800;

    /*reaction 230: HCN + M <=> H + CN + M */
    a[229] = 1.04e+29;
    b[229] = -3.3;
    e[229] = 126600;

    /*reaction 231: HCN + O <=> NCO + H */
    a[230] = 20300;
    b[230] = 2.64;
    e[230] = 4980;

    /*reaction 232: HCN + O <=> NH + CO */
    a[231] = 5070;
    b[231] = 2.64;
    e[231] = 4980;

    /*reaction 233: HCN + O <=> CN + OH */
    a[232] = 3.91e+09;
    b[232] = 1.58;
    e[232] = 26600;

    /*reaction 234: HCN + OH <=> HOCN + H */
    a[233] = 1.1e+06;
    b[233] = 2.03;
    e[233] = 13370;

    /*reaction 235: HCN + OH <=> HNCO + H */
    a[234] = 4400;
    b[234] = 2.26;
    e[234] = 6400;

    /*reaction 236: HCN + OH <=> NH2 + CO */
    a[235] = 160;
    b[235] = 2.56;
    e[235] = 9000;

    /*reaction 237: H + HCN (+M) <=> H2CN (+M) */
    a[236] = 3.3e+13;
    b[236] = 0;
    e[236] = 0;

    /*reaction 238: H2CN + N <=> N2 + CH2 */
    a[237] = 6e+13;
    b[237] = 0;
    e[237] = 400;

    /*reaction 239: C + N2 <=> CN + N */
    a[238] = 6.3e+13;
    b[238] = 0;
    e[238] = 46020;

    /*reaction 240: CH + N2 <=> HCN + N */
    a[239] = 3.12e+09;
    b[239] = 0.88;
    e[239] = 20130;

    /*reaction 241: CH + N2 (+M) <=> HCNN (+M) */
    a[240] = 3.1e+12;
    b[240] = 0.15;
    e[240] = 0;

    /*reaction 242: CH2 + N2 <=> HCN + NH */
    a[241] = 1e+13;
    b[241] = 0;
    e[241] = 74000;

    /*reaction 243: CH2(S) + N2 <=> NH + HCN */
    a[242] = 1e+11;
    b[242] = 0;
    e[242] = 65000;

    /*reaction 244: C + NO <=> CN + O */
    a[243] = 1.9e+13;
    b[243] = 0;
    e[243] = 0;

    /*reaction 245: C + NO <=> CO + N */
    a[244] = 2.9e+13;
    b[244] = 0;
    e[244] = 0;

    /*reaction 246: CH + NO <=> HCN + O */
    a[245] = 4.1e+13;
    b[245] = 0;
    e[245] = 0;

    /*reaction 247: CH + NO <=> H + NCO */
    a[246] = 1.62e+13;
    b[246] = 0;
    e[246] = 0;

    /*reaction 248: CH + NO <=> N + HCO */
    a[247] = 2.46e+13;
    b[247] = 0;
    e[247] = 0;

    /*reaction 249: CH2 + NO <=> H + HNCO */
    a[248] = 3.1e+17;
    b[248] = -1.38;
    e[248] = 1270;

    /*reaction 250: CH2 + NO <=> OH + HCN */
    a[249] = 2.9e+14;
    b[249] = -0.69;
    e[249] = 760;

    /*reaction 251: CH2 + NO <=> H + HCNO */
    a[250] = 3.8e+13;
    b[250] = -0.36;
    e[250] = 580;

    /*reaction 252: CH2(S) + NO <=> H + HNCO */
    a[251] = 3.1e+17;
    b[251] = -1.38;
    e[251] = 1270;

    /*reaction 253: CH2(S) + NO <=> OH + HCN */
    a[252] = 2.9e+14;
    b[252] = -0.69;
    e[252] = 760;

    /*reaction 254: CH2(S) + NO <=> H + HCNO */
    a[253] = 3.8e+13;
    b[253] = -0.36;
    e[253] = 580;

    /*reaction 255: CH3 + NO <=> HCN + H2O */
    a[254] = 9.6e+13;
    b[254] = 0;
    e[254] = 28800;

    /*reaction 256: CH3 + NO <=> H2CN + OH */
    a[255] = 1e+12;
    b[255] = 0;
    e[255] = 21750;

    /*reaction 257: HCNN + O <=> CO + H + N2 */
    a[256] = 2.2e+13;
    b[256] = 0;
    e[256] = 0;

    /*reaction 258: HCNN + O <=> HCN + NO */
    a[257] = 2e+12;
    b[257] = 0;
    e[257] = 0;

    /*reaction 259: HCNN + O2 <=> O + HCO + N2 */
    a[258] = 1.2e+13;
    b[258] = 0;
    e[258] = 0;

    /*reaction 260: HCNN + OH <=> H + HCO + N2 */
    a[259] = 1.2e+13;
    b[259] = 0;
    e[259] = 0;

    /*reaction 261: HCNN + H <=> CH2 + N2 */
    a[260] = 1e+14;
    b[260] = 0;
    e[260] = 0;

    /*reaction 262: HNCO + O <=> NH + CO2 */
    a[261] = 9.8e+07;
    b[261] = 1.41;
    e[261] = 8500;

    /*reaction 263: HNCO + O <=> HNO + CO */
    a[262] = 1.5e+08;
    b[262] = 1.57;
    e[262] = 44000;

    /*reaction 264: HNCO + O <=> NCO + OH */
    a[263] = 2.2e+06;
    b[263] = 2.11;
    e[263] = 11400;

    /*reaction 265: HNCO + H <=> NH2 + CO */
    a[264] = 2.25e+07;
    b[264] = 1.7;
    e[264] = 3800;

    /*reaction 266: HNCO + H <=> H2 + NCO */
    a[265] = 105000;
    b[265] = 2.5;
    e[265] = 13300;

    /*reaction 267: HNCO + OH <=> NCO + H2O */
    a[266] = 3.3e+07;
    b[266] = 1.5;
    e[266] = 3600;

    /*reaction 268: HNCO + OH <=> NH2 + CO2 */
    a[267] = 3.3e+06;
    b[267] = 1.5;
    e[267] = 3600;

    /*reaction 269: HNCO + M <=> NH + CO + M */
    a[268] = 1.18e+16;
    b[268] = 0;
    e[268] = 84720;

    /*reaction 270: HCNO + H <=> H + HNCO */
    a[269] = 2.1e+15;
    b[269] = -0.69;
    e[269] = 2850;

    /*reaction 271: HCNO + H <=> OH + HCN */
    a[270] = 2.7e+11;
    b[270] = 0.18;
    e[270] = 2120;

    /*reaction 272: HCNO + H <=> NH2 + CO */
    a[271] = 1.7e+14;
    b[271] = -0.75;
    e[271] = 2890;

    /*reaction 273: HOCN + H <=> H + HNCO */
    a[272] = 2e+07;
    b[272] = 2;
    e[272] = 2000;

    /*reaction 274: HCCO + NO <=> HCNO + CO */
    a[273] = 9e+12;
    b[273] = 0;
    e[273] = 0;

    /*reaction 275: CH3 + N <=> H2CN + H */
    a[274] = 6.1e+14;
    b[274] = -0.31;
    e[274] = 290;

    /*reaction 276: CH3 + N <=> HCN + H2 */
    a[275] = 3.7e+12;
    b[275] = 0.15;
    e[275] = -90;

    /*reaction 277: NH3 + H <=> NH2 + H2 */
    a[276] = 540000;
    b[276] = 2.4;
    e[276] = 9915;

    /*reaction 278: NH3 + OH <=> NH2 + H2O */
    a[277] = 5e+07;
    b[277] = 1.6;
    e[277] = 955;

    /*reaction 279: NH3 + O <=> NH2 + OH */
    a[278] = 9.4e+06;
    b[278] = 1.94;
    e[278] = 6460;

    /*reaction 280: NH + CO2 <=> HNO + CO */
    a[279] = 1e+13;
    b[279] = 0;
    e[279] = 14350;

    /*reaction 281: CN + NO2 <=> NCO + NO */
    a[280] = 6.16e+15;
    b[280] = -0.752;
    e[280] = 345;

    /*reaction 282: NCO + NO2 <=> N2O + CO2 */
    a[281] = 3.25e+12;
    b[281] = 0;
    e[281] = -705;

    /*reaction 283: N + CO2 <=> NO + CO */
    a[282] = 3e+12;
    b[282] = 0;
    e[282] = 11300;

    /*reaction 284: O + CH3 => H + H2 + CO */
    a[283] = 3.37e+13;
    b[283] = 0;
    e[283] = 0;

    /*reaction 285: O + C2H4 <=> H + CH2CHO */
    a[284] = 6.7e+06;
    b[284] = 1.83;
    e[284] = 220;

    /*reaction 286: O + C2H5 <=> H + CH3CHO */
    a[285] = 1.096e+14;
    b[285] = 0;
    e[285] = 0;

    /*reaction 287: OH + HO2 <=> O2 + H2O */
    a[286] = 5e+15;
    b[286] = 0;
    e[286] = 17330;

    /*reaction 288: OH + CH3 => H2 + CH2O */
    a[287] = 8e+09;
    b[287] = 0.5;
    e[287] = -1755;

    /*reaction 289: CH + H2 (+M) <=> CH3 (+M) */
    a[288] = 1.97e+12;
    b[288] = 0.43;
    e[288] = -370;

    /*reaction 290: CH2 + O2 => 2 H + CO2 */
    a[289] = 5.8e+12;
    b[289] = 0;
    e[289] = 1500;

    /*reaction 291: CH2 + O2 <=> O + CH2O */
    a[290] = 2.4e+12;
    b[290] = 0;
    e[290] = 1500;

    /*reaction 292: CH2 + CH2 => 2 H + C2H2 */
    a[291] = 2e+14;
    b[291] = 0;
    e[291] = 10989;

    /*reaction 293: CH2(S) + H2O => H2 + CH2O */
    a[292] = 6.82e+10;
    b[292] = 0.25;
    e[292] = -935;

    /*reaction 294: C2H3 + O2 <=> O + CH2CHO */
    a[293] = 3.03e+11;
    b[293] = 0.29;
    e[293] = 11;

    /*reaction 295: C2H3 + O2 <=> HO2 + C2H2 */
    a[294] = 1.337e+06;
    b[294] = 1.61;
    e[294] = -384;

    /*reaction 296: O + CH3CHO <=> OH + CH2CHO */
    a[295] = 5.84e+12;
    b[295] = 0;
    e[295] = 1808;

    /*reaction 297: O + CH3CHO => OH + CH3 + CO */
    a[296] = 5.84e+12;
    b[296] = 0;
    e[296] = 1808;

    /*reaction 298: O2 + CH3CHO => HO2 + CH3 + CO */
    a[297] = 3.01e+13;
    b[297] = 0;
    e[297] = 39150;

    /*reaction 299: H + CH3CHO <=> CH2CHO + H2 */
    a[298] = 2.05e+09;
    b[298] = 1.16;
    e[298] = 2405;

    /*reaction 300: H + CH3CHO => CH3 + H2 + CO */
    a[299] = 2.05e+09;
    b[299] = 1.16;
    e[299] = 2405;

    /*reaction 301: OH + CH3CHO => CH3 + H2O + CO */
    a[300] = 2.343e+10;
    b[300] = 0.73;
    e[300] = -1113;

    /*reaction 302: HO2 + CH3CHO => CH3 + H2O2 + CO */
    a[301] = 3.01e+12;
    b[301] = 0;
    e[301] = 11923;

    /*reaction 303: CH3 + CH3CHO => CH3 + CH4 + CO */
    a[302] = 2.72e+06;
    b[302] = 1.77;
    e[302] = 5920;

    /*reaction 304: H + CH2CO (+M) <=> CH2CHO (+M) */
    a[303] = 4.865e+11;
    b[303] = 0.422;
    e[303] = -1755;

    /*reaction 305: O + CH2CHO => H + CH2 + CO2 */
    a[304] = 1.5e+14;
    b[304] = 0;
    e[304] = 0;

    /*reaction 306: O2 + CH2CHO => OH + CO + CH2O */
    a[305] = 1.81e+10;
    b[305] = 0;
    e[305] = 0;

    /*reaction 307: O2 + CH2CHO => OH + 2 HCO */
    a[306] = 2.35e+10;
    b[306] = 0;
    e[306] = 0;

    /*reaction 308: H + CH2CHO <=> CH3 + HCO */
    a[307] = 2.2e+13;
    b[307] = 0;
    e[307] = 0;

    /*reaction 309: H + CH2CHO <=> CH2CO + H2 */
    a[308] = 1.1e+13;
    b[308] = 0;
    e[308] = 0;

    /*reaction 310: OH + CH2CHO <=> H2O + CH2CO */
    a[309] = 1.2e+13;
    b[309] = 0;
    e[309] = 0;

    /*reaction 311: OH + CH2CHO <=> HCO + CH2OH */
    a[310] = 3.01e+13;
    b[310] = 0;
    e[310] = 0;

    /*reaction 312: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    a[311] = 9.43e+12;
    b[311] = 0;
    e[311] = 0;

    /*reaction 313: O + C3H8 <=> OH + C3H7 */
    a[312] = 193000;
    b[312] = 2.68;
    e[312] = 3716;

    /*reaction 314: H + C3H8 <=> C3H7 + H2 */
    a[313] = 1.32e+06;
    b[313] = 2.54;
    e[313] = 6756;

    /*reaction 315: OH + C3H8 <=> C3H7 + H2O */
    a[314] = 3.16e+07;
    b[314] = 1.8;
    e[314] = 934;

    /*reaction 316: C3H7 + H2O2 <=> HO2 + C3H8 */
    a[315] = 378;
    b[315] = 2.72;
    e[315] = 1500;

    /*reaction 317: CH3 + C3H8 <=> C3H7 + CH4 */
    a[316] = 0.903;
    b[316] = 3.65;
    e[316] = 7154;

    /*reaction 318: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    a[317] = 2.55e+06;
    b[317] = 1.6;
    e[317] = 5700;

    /*reaction 319: O + C3H7 <=> C2H5 + CH2O */
    a[318] = 9.64e+13;
    b[318] = 0;
    e[318] = 0;

    /*reaction 320: H + C3H7 (+M) <=> C3H8 (+M) */
    a[319] = 3.613e+13;
    b[319] = 0;
    e[319] = 0;

    /*reaction 321: H + C3H7 <=> CH3 + C2H5 */
    a[320] = 4.06e+06;
    b[320] = 2.19;
    e[320] = 890;

    /*reaction 322: OH + C3H7 <=> C2H5 + CH2OH */
    a[321] = 2.41e+13;
    b[321] = 0;
    e[321] = 0;

    /*reaction 323: HO2 + C3H7 <=> O2 + C3H8 */
    a[322] = 2.55e+10;
    b[322] = 0.255;
    e[322] = -943;

    /*reaction 324: HO2 + C3H7 => OH + C2H5 + CH2O */
    a[323] = 2.41e+13;
    b[323] = 0;
    e[323] = 0;

    /*reaction 325: CH3 + C3H7 <=> 2 C2H5 */
    a[324] = 1.927e+13;
    b[324] = -0.32;
    e[324] = 0;

    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double * T, double * C, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[53]; /* temporary storage */

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

    /*reaction 12: O + CO (+M) <=> CO2 (+M) */
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

    /*reaction 135: CH2 + O2 => OH + H + CO */
    eqcon[134] *= 1e-06; 

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

    /*reaction 178: N + NO <=> N2 + O */
    /*eqcon[177] *= 1;  */

    /*reaction 179: N + O2 <=> NO + O */
    /*eqcon[178] *= 1;  */

    /*reaction 180: N + OH <=> NO + H */
    /*eqcon[179] *= 1;  */

    /*reaction 181: N2O + O <=> N2 + O2 */
    /*eqcon[180] *= 1;  */

    /*reaction 182: N2O + O <=> 2 NO */
    /*eqcon[181] *= 1;  */

    /*reaction 183: N2O + H <=> N2 + OH */
    /*eqcon[182] *= 1;  */

    /*reaction 184: N2O + OH <=> N2 + HO2 */
    /*eqcon[183] *= 1;  */

    /*reaction 185: N2O (+M) <=> N2 + O (+M) */
    eqcon[184] *= 1e-06; 

    /*reaction 186: HO2 + NO <=> NO2 + OH */
    /*eqcon[185] *= 1;  */

    /*reaction 187: NO + O + M <=> NO2 + M */
    eqcon[186] *= 1e+06; 

    /*reaction 188: NO2 + O <=> NO + O2 */
    /*eqcon[187] *= 1;  */

    /*reaction 189: NO2 + H <=> NO + OH */
    /*eqcon[188] *= 1;  */

    /*reaction 190: NH + O <=> NO + H */
    /*eqcon[189] *= 1;  */

    /*reaction 191: NH + H <=> N + H2 */
    /*eqcon[190] *= 1;  */

    /*reaction 192: NH + OH <=> HNO + H */
    /*eqcon[191] *= 1;  */

    /*reaction 193: NH + OH <=> N + H2O */
    /*eqcon[192] *= 1;  */

    /*reaction 194: NH + O2 <=> HNO + O */
    /*eqcon[193] *= 1;  */

    /*reaction 195: NH + O2 <=> NO + OH */
    /*eqcon[194] *= 1;  */

    /*reaction 196: NH + N <=> N2 + H */
    /*eqcon[195] *= 1;  */

    /*reaction 197: NH + H2O <=> HNO + H2 */
    /*eqcon[196] *= 1;  */

    /*reaction 198: NH + NO <=> N2 + OH */
    /*eqcon[197] *= 1;  */

    /*reaction 199: NH + NO <=> N2O + H */
    /*eqcon[198] *= 1;  */

    /*reaction 200: NH2 + O <=> OH + NH */
    /*eqcon[199] *= 1;  */

    /*reaction 201: NH2 + O <=> H + HNO */
    /*eqcon[200] *= 1;  */

    /*reaction 202: NH2 + H <=> NH + H2 */
    /*eqcon[201] *= 1;  */

    /*reaction 203: NH2 + OH <=> NH + H2O */
    /*eqcon[202] *= 1;  */

    /*reaction 204: NNH <=> N2 + H */
    eqcon[203] *= 1e-06; 

    /*reaction 205: NNH + M <=> N2 + H + M */
    eqcon[204] *= 1e-06; 

    /*reaction 206: NNH + O2 <=> HO2 + N2 */
    /*eqcon[205] *= 1;  */

    /*reaction 207: NNH + O <=> OH + N2 */
    /*eqcon[206] *= 1;  */

    /*reaction 208: NNH + O <=> NH + NO */
    /*eqcon[207] *= 1;  */

    /*reaction 209: NNH + H <=> H2 + N2 */
    /*eqcon[208] *= 1;  */

    /*reaction 210: NNH + OH <=> H2O + N2 */
    /*eqcon[209] *= 1;  */

    /*reaction 211: NNH + CH3 <=> CH4 + N2 */
    /*eqcon[210] *= 1;  */

    /*reaction 212: H + NO + M <=> HNO + M */
    eqcon[211] *= 1e+06; 

    /*reaction 213: HNO + O <=> NO + OH */
    /*eqcon[212] *= 1;  */

    /*reaction 214: HNO + H <=> H2 + NO */
    /*eqcon[213] *= 1;  */

    /*reaction 215: HNO + OH <=> NO + H2O */
    /*eqcon[214] *= 1;  */

    /*reaction 216: HNO + O2 <=> HO2 + NO */
    /*eqcon[215] *= 1;  */

    /*reaction 217: CN + O <=> CO + N */
    /*eqcon[216] *= 1;  */

    /*reaction 218: CN + OH <=> NCO + H */
    /*eqcon[217] *= 1;  */

    /*reaction 219: CN + H2O <=> HCN + OH */
    /*eqcon[218] *= 1;  */

    /*reaction 220: CN + O2 <=> NCO + O */
    /*eqcon[219] *= 1;  */

    /*reaction 221: CN + H2 <=> HCN + H */
    /*eqcon[220] *= 1;  */

    /*reaction 222: NCO + O <=> NO + CO */
    /*eqcon[221] *= 1;  */

    /*reaction 223: NCO + H <=> NH + CO */
    /*eqcon[222] *= 1;  */

    /*reaction 224: NCO + OH <=> NO + H + CO */
    eqcon[223] *= 1e-06; 

    /*reaction 225: NCO + N <=> N2 + CO */
    /*eqcon[224] *= 1;  */

    /*reaction 226: NCO + O2 <=> NO + CO2 */
    /*eqcon[225] *= 1;  */

    /*reaction 227: NCO + M <=> N + CO + M */
    eqcon[226] *= 1e-06; 

    /*reaction 228: NCO + NO <=> N2O + CO */
    /*eqcon[227] *= 1;  */

    /*reaction 229: NCO + NO <=> N2 + CO2 */
    /*eqcon[228] *= 1;  */

    /*reaction 230: HCN + M <=> H + CN + M */
    eqcon[229] *= 1e-06; 

    /*reaction 231: HCN + O <=> NCO + H */
    /*eqcon[230] *= 1;  */

    /*reaction 232: HCN + O <=> NH + CO */
    /*eqcon[231] *= 1;  */

    /*reaction 233: HCN + O <=> CN + OH */
    /*eqcon[232] *= 1;  */

    /*reaction 234: HCN + OH <=> HOCN + H */
    /*eqcon[233] *= 1;  */

    /*reaction 235: HCN + OH <=> HNCO + H */
    /*eqcon[234] *= 1;  */

    /*reaction 236: HCN + OH <=> NH2 + CO */
    /*eqcon[235] *= 1;  */

    /*reaction 237: H + HCN (+M) <=> H2CN (+M) */
    eqcon[236] *= 1e+06; 

    /*reaction 238: H2CN + N <=> N2 + CH2 */
    /*eqcon[237] *= 1;  */

    /*reaction 239: C + N2 <=> CN + N */
    /*eqcon[238] *= 1;  */

    /*reaction 240: CH + N2 <=> HCN + N */
    /*eqcon[239] *= 1;  */

    /*reaction 241: CH + N2 (+M) <=> HCNN (+M) */
    eqcon[240] *= 1e+06; 

    /*reaction 242: CH2 + N2 <=> HCN + NH */
    /*eqcon[241] *= 1;  */

    /*reaction 243: CH2(S) + N2 <=> NH + HCN */
    /*eqcon[242] *= 1;  */

    /*reaction 244: C + NO <=> CN + O */
    /*eqcon[243] *= 1;  */

    /*reaction 245: C + NO <=> CO + N */
    /*eqcon[244] *= 1;  */

    /*reaction 246: CH + NO <=> HCN + O */
    /*eqcon[245] *= 1;  */

    /*reaction 247: CH + NO <=> H + NCO */
    /*eqcon[246] *= 1;  */

    /*reaction 248: CH + NO <=> N + HCO */
    /*eqcon[247] *= 1;  */

    /*reaction 249: CH2 + NO <=> H + HNCO */
    /*eqcon[248] *= 1;  */

    /*reaction 250: CH2 + NO <=> OH + HCN */
    /*eqcon[249] *= 1;  */

    /*reaction 251: CH2 + NO <=> H + HCNO */
    /*eqcon[250] *= 1;  */

    /*reaction 252: CH2(S) + NO <=> H + HNCO */
    /*eqcon[251] *= 1;  */

    /*reaction 253: CH2(S) + NO <=> OH + HCN */
    /*eqcon[252] *= 1;  */

    /*reaction 254: CH2(S) + NO <=> H + HCNO */
    /*eqcon[253] *= 1;  */

    /*reaction 255: CH3 + NO <=> HCN + H2O */
    /*eqcon[254] *= 1;  */

    /*reaction 256: CH3 + NO <=> H2CN + OH */
    /*eqcon[255] *= 1;  */

    /*reaction 257: HCNN + O <=> CO + H + N2 */
    eqcon[256] *= 1e-06; 

    /*reaction 258: HCNN + O <=> HCN + NO */
    /*eqcon[257] *= 1;  */

    /*reaction 259: HCNN + O2 <=> O + HCO + N2 */
    eqcon[258] *= 1e-06; 

    /*reaction 260: HCNN + OH <=> H + HCO + N2 */
    eqcon[259] *= 1e-06; 

    /*reaction 261: HCNN + H <=> CH2 + N2 */
    /*eqcon[260] *= 1;  */

    /*reaction 262: HNCO + O <=> NH + CO2 */
    /*eqcon[261] *= 1;  */

    /*reaction 263: HNCO + O <=> HNO + CO */
    /*eqcon[262] *= 1;  */

    /*reaction 264: HNCO + O <=> NCO + OH */
    /*eqcon[263] *= 1;  */

    /*reaction 265: HNCO + H <=> NH2 + CO */
    /*eqcon[264] *= 1;  */

    /*reaction 266: HNCO + H <=> H2 + NCO */
    /*eqcon[265] *= 1;  */

    /*reaction 267: HNCO + OH <=> NCO + H2O */
    /*eqcon[266] *= 1;  */

    /*reaction 268: HNCO + OH <=> NH2 + CO2 */
    /*eqcon[267] *= 1;  */

    /*reaction 269: HNCO + M <=> NH + CO + M */
    eqcon[268] *= 1e-06; 

    /*reaction 270: HCNO + H <=> H + HNCO */
    /*eqcon[269] *= 1;  */

    /*reaction 271: HCNO + H <=> OH + HCN */
    /*eqcon[270] *= 1;  */

    /*reaction 272: HCNO + H <=> NH2 + CO */
    /*eqcon[271] *= 1;  */

    /*reaction 273: HOCN + H <=> H + HNCO */
    /*eqcon[272] *= 1;  */

    /*reaction 274: HCCO + NO <=> HCNO + CO */
    /*eqcon[273] *= 1;  */

    /*reaction 275: CH3 + N <=> H2CN + H */
    /*eqcon[274] *= 1;  */

    /*reaction 276: CH3 + N <=> HCN + H2 */
    /*eqcon[275] *= 1;  */

    /*reaction 277: NH3 + H <=> NH2 + H2 */
    /*eqcon[276] *= 1;  */

    /*reaction 278: NH3 + OH <=> NH2 + H2O */
    /*eqcon[277] *= 1;  */

    /*reaction 279: NH3 + O <=> NH2 + OH */
    /*eqcon[278] *= 1;  */

    /*reaction 280: NH + CO2 <=> HNO + CO */
    /*eqcon[279] *= 1;  */

    /*reaction 281: CN + NO2 <=> NCO + NO */
    /*eqcon[280] *= 1;  */

    /*reaction 282: NCO + NO2 <=> N2O + CO2 */
    /*eqcon[281] *= 1;  */

    /*reaction 283: N + CO2 <=> NO + CO */
    /*eqcon[282] *= 1;  */

    /*reaction 284: O + CH3 => H + H2 + CO */
    eqcon[283] *= 1e-06; 

    /*reaction 285: O + C2H4 <=> H + CH2CHO */
    /*eqcon[284] *= 1;  */

    /*reaction 286: O + C2H5 <=> H + CH3CHO */
    /*eqcon[285] *= 1;  */

    /*reaction 287: OH + HO2 <=> O2 + H2O */
    /*eqcon[286] *= 1;  */

    /*reaction 288: OH + CH3 => H2 + CH2O */
    /*eqcon[287] *= 1;  */

    /*reaction 289: CH + H2 (+M) <=> CH3 (+M) */
    eqcon[288] *= 1e+06; 

    /*reaction 290: CH2 + O2 => 2 H + CO2 */
    eqcon[289] *= 1e-06; 

    /*reaction 291: CH2 + O2 <=> O + CH2O */
    /*eqcon[290] *= 1;  */

    /*reaction 292: CH2 + CH2 => 2 H + C2H2 */
    eqcon[291] *= 1e-06; 

    /*reaction 293: CH2(S) + H2O => H2 + CH2O */
    /*eqcon[292] *= 1;  */

    /*reaction 294: C2H3 + O2 <=> O + CH2CHO */
    /*eqcon[293] *= 1;  */

    /*reaction 295: C2H3 + O2 <=> HO2 + C2H2 */
    /*eqcon[294] *= 1;  */

    /*reaction 296: O + CH3CHO <=> OH + CH2CHO */
    /*eqcon[295] *= 1;  */

    /*reaction 297: O + CH3CHO => OH + CH3 + CO */
    eqcon[296] *= 1e-06; 

    /*reaction 298: O2 + CH3CHO => HO2 + CH3 + CO */
    eqcon[297] *= 1e-06; 

    /*reaction 299: H + CH3CHO <=> CH2CHO + H2 */
    /*eqcon[298] *= 1;  */

    /*reaction 300: H + CH3CHO => CH3 + H2 + CO */
    eqcon[299] *= 1e-06; 

    /*reaction 301: OH + CH3CHO => CH3 + H2O + CO */
    eqcon[300] *= 1e-06; 

    /*reaction 302: HO2 + CH3CHO => CH3 + H2O2 + CO */
    eqcon[301] *= 1e-06; 

    /*reaction 303: CH3 + CH3CHO => CH3 + CH4 + CO */
    eqcon[302] *= 1e-06; 

    /*reaction 304: H + CH2CO (+M) <=> CH2CHO (+M) */
    eqcon[303] *= 1e+06; 

    /*reaction 305: O + CH2CHO => H + CH2 + CO2 */
    eqcon[304] *= 1e-06; 

    /*reaction 306: O2 + CH2CHO => OH + CO + CH2O */
    eqcon[305] *= 1e-06; 

    /*reaction 307: O2 + CH2CHO => OH + 2 HCO */
    eqcon[306] *= 1e-06; 

    /*reaction 308: H + CH2CHO <=> CH3 + HCO */
    /*eqcon[307] *= 1;  */

    /*reaction 309: H + CH2CHO <=> CH2CO + H2 */
    /*eqcon[308] *= 1;  */

    /*reaction 310: OH + CH2CHO <=> H2O + CH2CO */
    /*eqcon[309] *= 1;  */

    /*reaction 311: OH + CH2CHO <=> HCO + CH2OH */
    /*eqcon[310] *= 1;  */

    /*reaction 312: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    eqcon[311] *= 1e+06; 

    /*reaction 313: O + C3H8 <=> OH + C3H7 */
    /*eqcon[312] *= 1;  */

    /*reaction 314: H + C3H8 <=> C3H7 + H2 */
    /*eqcon[313] *= 1;  */

    /*reaction 315: OH + C3H8 <=> C3H7 + H2O */
    /*eqcon[314] *= 1;  */

    /*reaction 316: C3H7 + H2O2 <=> HO2 + C3H8 */
    /*eqcon[315] *= 1;  */

    /*reaction 317: CH3 + C3H8 <=> C3H7 + CH4 */
    /*eqcon[316] *= 1;  */

    /*reaction 318: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    eqcon[317] *= 1e+06; 

    /*reaction 319: O + C3H7 <=> C2H5 + CH2O */
    /*eqcon[318] *= 1;  */

    /*reaction 320: H + C3H7 (+M) <=> C3H8 (+M) */
    eqcon[319] *= 1e+06; 

    /*reaction 321: H + C3H7 <=> CH3 + C2H5 */
    /*eqcon[320] *= 1;  */

    /*reaction 322: OH + C3H7 <=> C2H5 + CH2OH */
    /*eqcon[321] *= 1;  */

    /*reaction 323: HO2 + C3H7 <=> O2 + C3H8 */
    /*eqcon[322] *= 1;  */

    /*reaction 324: HO2 + C3H7 => OH + C2H5 + CH2O */
    eqcon[323] *= 1e-06; 

    /*reaction 325: CH3 + C3H7 <=> 2 C2H5 */
    /*eqcon[324] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[53]; /* temporary storage */

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

    /*reaction 12: O + CO (+M) <=> CO2 (+M) */
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

    /*reaction 135: CH2 + O2 => OH + H + CO */
    eqcon[134] *= 1e-06; 

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

    /*reaction 178: N + NO <=> N2 + O */
    /*eqcon[177] *= 1;  */

    /*reaction 179: N + O2 <=> NO + O */
    /*eqcon[178] *= 1;  */

    /*reaction 180: N + OH <=> NO + H */
    /*eqcon[179] *= 1;  */

    /*reaction 181: N2O + O <=> N2 + O2 */
    /*eqcon[180] *= 1;  */

    /*reaction 182: N2O + O <=> 2 NO */
    /*eqcon[181] *= 1;  */

    /*reaction 183: N2O + H <=> N2 + OH */
    /*eqcon[182] *= 1;  */

    /*reaction 184: N2O + OH <=> N2 + HO2 */
    /*eqcon[183] *= 1;  */

    /*reaction 185: N2O (+M) <=> N2 + O (+M) */
    eqcon[184] *= 1e-06; 

    /*reaction 186: HO2 + NO <=> NO2 + OH */
    /*eqcon[185] *= 1;  */

    /*reaction 187: NO + O + M <=> NO2 + M */
    eqcon[186] *= 1e+06; 

    /*reaction 188: NO2 + O <=> NO + O2 */
    /*eqcon[187] *= 1;  */

    /*reaction 189: NO2 + H <=> NO + OH */
    /*eqcon[188] *= 1;  */

    /*reaction 190: NH + O <=> NO + H */
    /*eqcon[189] *= 1;  */

    /*reaction 191: NH + H <=> N + H2 */
    /*eqcon[190] *= 1;  */

    /*reaction 192: NH + OH <=> HNO + H */
    /*eqcon[191] *= 1;  */

    /*reaction 193: NH + OH <=> N + H2O */
    /*eqcon[192] *= 1;  */

    /*reaction 194: NH + O2 <=> HNO + O */
    /*eqcon[193] *= 1;  */

    /*reaction 195: NH + O2 <=> NO + OH */
    /*eqcon[194] *= 1;  */

    /*reaction 196: NH + N <=> N2 + H */
    /*eqcon[195] *= 1;  */

    /*reaction 197: NH + H2O <=> HNO + H2 */
    /*eqcon[196] *= 1;  */

    /*reaction 198: NH + NO <=> N2 + OH */
    /*eqcon[197] *= 1;  */

    /*reaction 199: NH + NO <=> N2O + H */
    /*eqcon[198] *= 1;  */

    /*reaction 200: NH2 + O <=> OH + NH */
    /*eqcon[199] *= 1;  */

    /*reaction 201: NH2 + O <=> H + HNO */
    /*eqcon[200] *= 1;  */

    /*reaction 202: NH2 + H <=> NH + H2 */
    /*eqcon[201] *= 1;  */

    /*reaction 203: NH2 + OH <=> NH + H2O */
    /*eqcon[202] *= 1;  */

    /*reaction 204: NNH <=> N2 + H */
    eqcon[203] *= 1e-06; 

    /*reaction 205: NNH + M <=> N2 + H + M */
    eqcon[204] *= 1e-06; 

    /*reaction 206: NNH + O2 <=> HO2 + N2 */
    /*eqcon[205] *= 1;  */

    /*reaction 207: NNH + O <=> OH + N2 */
    /*eqcon[206] *= 1;  */

    /*reaction 208: NNH + O <=> NH + NO */
    /*eqcon[207] *= 1;  */

    /*reaction 209: NNH + H <=> H2 + N2 */
    /*eqcon[208] *= 1;  */

    /*reaction 210: NNH + OH <=> H2O + N2 */
    /*eqcon[209] *= 1;  */

    /*reaction 211: NNH + CH3 <=> CH4 + N2 */
    /*eqcon[210] *= 1;  */

    /*reaction 212: H + NO + M <=> HNO + M */
    eqcon[211] *= 1e+06; 

    /*reaction 213: HNO + O <=> NO + OH */
    /*eqcon[212] *= 1;  */

    /*reaction 214: HNO + H <=> H2 + NO */
    /*eqcon[213] *= 1;  */

    /*reaction 215: HNO + OH <=> NO + H2O */
    /*eqcon[214] *= 1;  */

    /*reaction 216: HNO + O2 <=> HO2 + NO */
    /*eqcon[215] *= 1;  */

    /*reaction 217: CN + O <=> CO + N */
    /*eqcon[216] *= 1;  */

    /*reaction 218: CN + OH <=> NCO + H */
    /*eqcon[217] *= 1;  */

    /*reaction 219: CN + H2O <=> HCN + OH */
    /*eqcon[218] *= 1;  */

    /*reaction 220: CN + O2 <=> NCO + O */
    /*eqcon[219] *= 1;  */

    /*reaction 221: CN + H2 <=> HCN + H */
    /*eqcon[220] *= 1;  */

    /*reaction 222: NCO + O <=> NO + CO */
    /*eqcon[221] *= 1;  */

    /*reaction 223: NCO + H <=> NH + CO */
    /*eqcon[222] *= 1;  */

    /*reaction 224: NCO + OH <=> NO + H + CO */
    eqcon[223] *= 1e-06; 

    /*reaction 225: NCO + N <=> N2 + CO */
    /*eqcon[224] *= 1;  */

    /*reaction 226: NCO + O2 <=> NO + CO2 */
    /*eqcon[225] *= 1;  */

    /*reaction 227: NCO + M <=> N + CO + M */
    eqcon[226] *= 1e-06; 

    /*reaction 228: NCO + NO <=> N2O + CO */
    /*eqcon[227] *= 1;  */

    /*reaction 229: NCO + NO <=> N2 + CO2 */
    /*eqcon[228] *= 1;  */

    /*reaction 230: HCN + M <=> H + CN + M */
    eqcon[229] *= 1e-06; 

    /*reaction 231: HCN + O <=> NCO + H */
    /*eqcon[230] *= 1;  */

    /*reaction 232: HCN + O <=> NH + CO */
    /*eqcon[231] *= 1;  */

    /*reaction 233: HCN + O <=> CN + OH */
    /*eqcon[232] *= 1;  */

    /*reaction 234: HCN + OH <=> HOCN + H */
    /*eqcon[233] *= 1;  */

    /*reaction 235: HCN + OH <=> HNCO + H */
    /*eqcon[234] *= 1;  */

    /*reaction 236: HCN + OH <=> NH2 + CO */
    /*eqcon[235] *= 1;  */

    /*reaction 237: H + HCN (+M) <=> H2CN (+M) */
    eqcon[236] *= 1e+06; 

    /*reaction 238: H2CN + N <=> N2 + CH2 */
    /*eqcon[237] *= 1;  */

    /*reaction 239: C + N2 <=> CN + N */
    /*eqcon[238] *= 1;  */

    /*reaction 240: CH + N2 <=> HCN + N */
    /*eqcon[239] *= 1;  */

    /*reaction 241: CH + N2 (+M) <=> HCNN (+M) */
    eqcon[240] *= 1e+06; 

    /*reaction 242: CH2 + N2 <=> HCN + NH */
    /*eqcon[241] *= 1;  */

    /*reaction 243: CH2(S) + N2 <=> NH + HCN */
    /*eqcon[242] *= 1;  */

    /*reaction 244: C + NO <=> CN + O */
    /*eqcon[243] *= 1;  */

    /*reaction 245: C + NO <=> CO + N */
    /*eqcon[244] *= 1;  */

    /*reaction 246: CH + NO <=> HCN + O */
    /*eqcon[245] *= 1;  */

    /*reaction 247: CH + NO <=> H + NCO */
    /*eqcon[246] *= 1;  */

    /*reaction 248: CH + NO <=> N + HCO */
    /*eqcon[247] *= 1;  */

    /*reaction 249: CH2 + NO <=> H + HNCO */
    /*eqcon[248] *= 1;  */

    /*reaction 250: CH2 + NO <=> OH + HCN */
    /*eqcon[249] *= 1;  */

    /*reaction 251: CH2 + NO <=> H + HCNO */
    /*eqcon[250] *= 1;  */

    /*reaction 252: CH2(S) + NO <=> H + HNCO */
    /*eqcon[251] *= 1;  */

    /*reaction 253: CH2(S) + NO <=> OH + HCN */
    /*eqcon[252] *= 1;  */

    /*reaction 254: CH2(S) + NO <=> H + HCNO */
    /*eqcon[253] *= 1;  */

    /*reaction 255: CH3 + NO <=> HCN + H2O */
    /*eqcon[254] *= 1;  */

    /*reaction 256: CH3 + NO <=> H2CN + OH */
    /*eqcon[255] *= 1;  */

    /*reaction 257: HCNN + O <=> CO + H + N2 */
    eqcon[256] *= 1e-06; 

    /*reaction 258: HCNN + O <=> HCN + NO */
    /*eqcon[257] *= 1;  */

    /*reaction 259: HCNN + O2 <=> O + HCO + N2 */
    eqcon[258] *= 1e-06; 

    /*reaction 260: HCNN + OH <=> H + HCO + N2 */
    eqcon[259] *= 1e-06; 

    /*reaction 261: HCNN + H <=> CH2 + N2 */
    /*eqcon[260] *= 1;  */

    /*reaction 262: HNCO + O <=> NH + CO2 */
    /*eqcon[261] *= 1;  */

    /*reaction 263: HNCO + O <=> HNO + CO */
    /*eqcon[262] *= 1;  */

    /*reaction 264: HNCO + O <=> NCO + OH */
    /*eqcon[263] *= 1;  */

    /*reaction 265: HNCO + H <=> NH2 + CO */
    /*eqcon[264] *= 1;  */

    /*reaction 266: HNCO + H <=> H2 + NCO */
    /*eqcon[265] *= 1;  */

    /*reaction 267: HNCO + OH <=> NCO + H2O */
    /*eqcon[266] *= 1;  */

    /*reaction 268: HNCO + OH <=> NH2 + CO2 */
    /*eqcon[267] *= 1;  */

    /*reaction 269: HNCO + M <=> NH + CO + M */
    eqcon[268] *= 1e-06; 

    /*reaction 270: HCNO + H <=> H + HNCO */
    /*eqcon[269] *= 1;  */

    /*reaction 271: HCNO + H <=> OH + HCN */
    /*eqcon[270] *= 1;  */

    /*reaction 272: HCNO + H <=> NH2 + CO */
    /*eqcon[271] *= 1;  */

    /*reaction 273: HOCN + H <=> H + HNCO */
    /*eqcon[272] *= 1;  */

    /*reaction 274: HCCO + NO <=> HCNO + CO */
    /*eqcon[273] *= 1;  */

    /*reaction 275: CH3 + N <=> H2CN + H */
    /*eqcon[274] *= 1;  */

    /*reaction 276: CH3 + N <=> HCN + H2 */
    /*eqcon[275] *= 1;  */

    /*reaction 277: NH3 + H <=> NH2 + H2 */
    /*eqcon[276] *= 1;  */

    /*reaction 278: NH3 + OH <=> NH2 + H2O */
    /*eqcon[277] *= 1;  */

    /*reaction 279: NH3 + O <=> NH2 + OH */
    /*eqcon[278] *= 1;  */

    /*reaction 280: NH + CO2 <=> HNO + CO */
    /*eqcon[279] *= 1;  */

    /*reaction 281: CN + NO2 <=> NCO + NO */
    /*eqcon[280] *= 1;  */

    /*reaction 282: NCO + NO2 <=> N2O + CO2 */
    /*eqcon[281] *= 1;  */

    /*reaction 283: N + CO2 <=> NO + CO */
    /*eqcon[282] *= 1;  */

    /*reaction 284: O + CH3 => H + H2 + CO */
    eqcon[283] *= 1e-06; 

    /*reaction 285: O + C2H4 <=> H + CH2CHO */
    /*eqcon[284] *= 1;  */

    /*reaction 286: O + C2H5 <=> H + CH3CHO */
    /*eqcon[285] *= 1;  */

    /*reaction 287: OH + HO2 <=> O2 + H2O */
    /*eqcon[286] *= 1;  */

    /*reaction 288: OH + CH3 => H2 + CH2O */
    /*eqcon[287] *= 1;  */

    /*reaction 289: CH + H2 (+M) <=> CH3 (+M) */
    eqcon[288] *= 1e+06; 

    /*reaction 290: CH2 + O2 => 2 H + CO2 */
    eqcon[289] *= 1e-06; 

    /*reaction 291: CH2 + O2 <=> O + CH2O */
    /*eqcon[290] *= 1;  */

    /*reaction 292: CH2 + CH2 => 2 H + C2H2 */
    eqcon[291] *= 1e-06; 

    /*reaction 293: CH2(S) + H2O => H2 + CH2O */
    /*eqcon[292] *= 1;  */

    /*reaction 294: C2H3 + O2 <=> O + CH2CHO */
    /*eqcon[293] *= 1;  */

    /*reaction 295: C2H3 + O2 <=> HO2 + C2H2 */
    /*eqcon[294] *= 1;  */

    /*reaction 296: O + CH3CHO <=> OH + CH2CHO */
    /*eqcon[295] *= 1;  */

    /*reaction 297: O + CH3CHO => OH + CH3 + CO */
    eqcon[296] *= 1e-06; 

    /*reaction 298: O2 + CH3CHO => HO2 + CH3 + CO */
    eqcon[297] *= 1e-06; 

    /*reaction 299: H + CH3CHO <=> CH2CHO + H2 */
    /*eqcon[298] *= 1;  */

    /*reaction 300: H + CH3CHO => CH3 + H2 + CO */
    eqcon[299] *= 1e-06; 

    /*reaction 301: OH + CH3CHO => CH3 + H2O + CO */
    eqcon[300] *= 1e-06; 

    /*reaction 302: HO2 + CH3CHO => CH3 + H2O2 + CO */
    eqcon[301] *= 1e-06; 

    /*reaction 303: CH3 + CH3CHO => CH3 + CH4 + CO */
    eqcon[302] *= 1e-06; 

    /*reaction 304: H + CH2CO (+M) <=> CH2CHO (+M) */
    eqcon[303] *= 1e+06; 

    /*reaction 305: O + CH2CHO => H + CH2 + CO2 */
    eqcon[304] *= 1e-06; 

    /*reaction 306: O2 + CH2CHO => OH + CO + CH2O */
    eqcon[305] *= 1e-06; 

    /*reaction 307: O2 + CH2CHO => OH + 2 HCO */
    eqcon[306] *= 1e-06; 

    /*reaction 308: H + CH2CHO <=> CH3 + HCO */
    /*eqcon[307] *= 1;  */

    /*reaction 309: H + CH2CHO <=> CH2CO + H2 */
    /*eqcon[308] *= 1;  */

    /*reaction 310: OH + CH2CHO <=> H2O + CH2CO */
    /*eqcon[309] *= 1;  */

    /*reaction 311: OH + CH2CHO <=> HCO + CH2OH */
    /*eqcon[310] *= 1;  */

    /*reaction 312: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    eqcon[311] *= 1e+06; 

    /*reaction 313: O + C3H8 <=> OH + C3H7 */
    /*eqcon[312] *= 1;  */

    /*reaction 314: H + C3H8 <=> C3H7 + H2 */
    /*eqcon[313] *= 1;  */

    /*reaction 315: OH + C3H8 <=> C3H7 + H2O */
    /*eqcon[314] *= 1;  */

    /*reaction 316: C3H7 + H2O2 <=> HO2 + C3H8 */
    /*eqcon[315] *= 1;  */

    /*reaction 317: CH3 + C3H8 <=> C3H7 + CH4 */
    /*eqcon[316] *= 1;  */

    /*reaction 318: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    eqcon[317] *= 1e+06; 

    /*reaction 319: O + C3H7 <=> C2H5 + CH2O */
    /*eqcon[318] *= 1;  */

    /*reaction 320: H + C3H7 (+M) <=> C3H8 (+M) */
    eqcon[319] *= 1e+06; 

    /*reaction 321: H + C3H7 <=> CH3 + C2H5 */
    /*eqcon[320] *= 1;  */

    /*reaction 322: OH + C3H7 <=> C2H5 + CH2OH */
    /*eqcon[321] *= 1;  */

    /*reaction 323: HO2 + C3H7 <=> O2 + C3H8 */
    /*eqcon[322] *= 1;  */

    /*reaction 324: HO2 + C3H7 => OH + C2H5 + CH2O */
    eqcon[323] *= 1e-06; 

    /*reaction 325: CH3 + C3H7 <=> 2 C2H5 */
    /*eqcon[324] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[53]; /* temporary storage */

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

    /*reaction 12: O + CO (+M) <=> CO2 (+M) */
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

    /*reaction 135: CH2 + O2 => OH + H + CO */
    eqcon[134] *= 1e-06; 

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

    /*reaction 178: N + NO <=> N2 + O */
    /*eqcon[177] *= 1;  */

    /*reaction 179: N + O2 <=> NO + O */
    /*eqcon[178] *= 1;  */

    /*reaction 180: N + OH <=> NO + H */
    /*eqcon[179] *= 1;  */

    /*reaction 181: N2O + O <=> N2 + O2 */
    /*eqcon[180] *= 1;  */

    /*reaction 182: N2O + O <=> 2 NO */
    /*eqcon[181] *= 1;  */

    /*reaction 183: N2O + H <=> N2 + OH */
    /*eqcon[182] *= 1;  */

    /*reaction 184: N2O + OH <=> N2 + HO2 */
    /*eqcon[183] *= 1;  */

    /*reaction 185: N2O (+M) <=> N2 + O (+M) */
    eqcon[184] *= 1e-06; 

    /*reaction 186: HO2 + NO <=> NO2 + OH */
    /*eqcon[185] *= 1;  */

    /*reaction 187: NO + O + M <=> NO2 + M */
    eqcon[186] *= 1e+06; 

    /*reaction 188: NO2 + O <=> NO + O2 */
    /*eqcon[187] *= 1;  */

    /*reaction 189: NO2 + H <=> NO + OH */
    /*eqcon[188] *= 1;  */

    /*reaction 190: NH + O <=> NO + H */
    /*eqcon[189] *= 1;  */

    /*reaction 191: NH + H <=> N + H2 */
    /*eqcon[190] *= 1;  */

    /*reaction 192: NH + OH <=> HNO + H */
    /*eqcon[191] *= 1;  */

    /*reaction 193: NH + OH <=> N + H2O */
    /*eqcon[192] *= 1;  */

    /*reaction 194: NH + O2 <=> HNO + O */
    /*eqcon[193] *= 1;  */

    /*reaction 195: NH + O2 <=> NO + OH */
    /*eqcon[194] *= 1;  */

    /*reaction 196: NH + N <=> N2 + H */
    /*eqcon[195] *= 1;  */

    /*reaction 197: NH + H2O <=> HNO + H2 */
    /*eqcon[196] *= 1;  */

    /*reaction 198: NH + NO <=> N2 + OH */
    /*eqcon[197] *= 1;  */

    /*reaction 199: NH + NO <=> N2O + H */
    /*eqcon[198] *= 1;  */

    /*reaction 200: NH2 + O <=> OH + NH */
    /*eqcon[199] *= 1;  */

    /*reaction 201: NH2 + O <=> H + HNO */
    /*eqcon[200] *= 1;  */

    /*reaction 202: NH2 + H <=> NH + H2 */
    /*eqcon[201] *= 1;  */

    /*reaction 203: NH2 + OH <=> NH + H2O */
    /*eqcon[202] *= 1;  */

    /*reaction 204: NNH <=> N2 + H */
    eqcon[203] *= 1e-06; 

    /*reaction 205: NNH + M <=> N2 + H + M */
    eqcon[204] *= 1e-06; 

    /*reaction 206: NNH + O2 <=> HO2 + N2 */
    /*eqcon[205] *= 1;  */

    /*reaction 207: NNH + O <=> OH + N2 */
    /*eqcon[206] *= 1;  */

    /*reaction 208: NNH + O <=> NH + NO */
    /*eqcon[207] *= 1;  */

    /*reaction 209: NNH + H <=> H2 + N2 */
    /*eqcon[208] *= 1;  */

    /*reaction 210: NNH + OH <=> H2O + N2 */
    /*eqcon[209] *= 1;  */

    /*reaction 211: NNH + CH3 <=> CH4 + N2 */
    /*eqcon[210] *= 1;  */

    /*reaction 212: H + NO + M <=> HNO + M */
    eqcon[211] *= 1e+06; 

    /*reaction 213: HNO + O <=> NO + OH */
    /*eqcon[212] *= 1;  */

    /*reaction 214: HNO + H <=> H2 + NO */
    /*eqcon[213] *= 1;  */

    /*reaction 215: HNO + OH <=> NO + H2O */
    /*eqcon[214] *= 1;  */

    /*reaction 216: HNO + O2 <=> HO2 + NO */
    /*eqcon[215] *= 1;  */

    /*reaction 217: CN + O <=> CO + N */
    /*eqcon[216] *= 1;  */

    /*reaction 218: CN + OH <=> NCO + H */
    /*eqcon[217] *= 1;  */

    /*reaction 219: CN + H2O <=> HCN + OH */
    /*eqcon[218] *= 1;  */

    /*reaction 220: CN + O2 <=> NCO + O */
    /*eqcon[219] *= 1;  */

    /*reaction 221: CN + H2 <=> HCN + H */
    /*eqcon[220] *= 1;  */

    /*reaction 222: NCO + O <=> NO + CO */
    /*eqcon[221] *= 1;  */

    /*reaction 223: NCO + H <=> NH + CO */
    /*eqcon[222] *= 1;  */

    /*reaction 224: NCO + OH <=> NO + H + CO */
    eqcon[223] *= 1e-06; 

    /*reaction 225: NCO + N <=> N2 + CO */
    /*eqcon[224] *= 1;  */

    /*reaction 226: NCO + O2 <=> NO + CO2 */
    /*eqcon[225] *= 1;  */

    /*reaction 227: NCO + M <=> N + CO + M */
    eqcon[226] *= 1e-06; 

    /*reaction 228: NCO + NO <=> N2O + CO */
    /*eqcon[227] *= 1;  */

    /*reaction 229: NCO + NO <=> N2 + CO2 */
    /*eqcon[228] *= 1;  */

    /*reaction 230: HCN + M <=> H + CN + M */
    eqcon[229] *= 1e-06; 

    /*reaction 231: HCN + O <=> NCO + H */
    /*eqcon[230] *= 1;  */

    /*reaction 232: HCN + O <=> NH + CO */
    /*eqcon[231] *= 1;  */

    /*reaction 233: HCN + O <=> CN + OH */
    /*eqcon[232] *= 1;  */

    /*reaction 234: HCN + OH <=> HOCN + H */
    /*eqcon[233] *= 1;  */

    /*reaction 235: HCN + OH <=> HNCO + H */
    /*eqcon[234] *= 1;  */

    /*reaction 236: HCN + OH <=> NH2 + CO */
    /*eqcon[235] *= 1;  */

    /*reaction 237: H + HCN (+M) <=> H2CN (+M) */
    eqcon[236] *= 1e+06; 

    /*reaction 238: H2CN + N <=> N2 + CH2 */
    /*eqcon[237] *= 1;  */

    /*reaction 239: C + N2 <=> CN + N */
    /*eqcon[238] *= 1;  */

    /*reaction 240: CH + N2 <=> HCN + N */
    /*eqcon[239] *= 1;  */

    /*reaction 241: CH + N2 (+M) <=> HCNN (+M) */
    eqcon[240] *= 1e+06; 

    /*reaction 242: CH2 + N2 <=> HCN + NH */
    /*eqcon[241] *= 1;  */

    /*reaction 243: CH2(S) + N2 <=> NH + HCN */
    /*eqcon[242] *= 1;  */

    /*reaction 244: C + NO <=> CN + O */
    /*eqcon[243] *= 1;  */

    /*reaction 245: C + NO <=> CO + N */
    /*eqcon[244] *= 1;  */

    /*reaction 246: CH + NO <=> HCN + O */
    /*eqcon[245] *= 1;  */

    /*reaction 247: CH + NO <=> H + NCO */
    /*eqcon[246] *= 1;  */

    /*reaction 248: CH + NO <=> N + HCO */
    /*eqcon[247] *= 1;  */

    /*reaction 249: CH2 + NO <=> H + HNCO */
    /*eqcon[248] *= 1;  */

    /*reaction 250: CH2 + NO <=> OH + HCN */
    /*eqcon[249] *= 1;  */

    /*reaction 251: CH2 + NO <=> H + HCNO */
    /*eqcon[250] *= 1;  */

    /*reaction 252: CH2(S) + NO <=> H + HNCO */
    /*eqcon[251] *= 1;  */

    /*reaction 253: CH2(S) + NO <=> OH + HCN */
    /*eqcon[252] *= 1;  */

    /*reaction 254: CH2(S) + NO <=> H + HCNO */
    /*eqcon[253] *= 1;  */

    /*reaction 255: CH3 + NO <=> HCN + H2O */
    /*eqcon[254] *= 1;  */

    /*reaction 256: CH3 + NO <=> H2CN + OH */
    /*eqcon[255] *= 1;  */

    /*reaction 257: HCNN + O <=> CO + H + N2 */
    eqcon[256] *= 1e-06; 

    /*reaction 258: HCNN + O <=> HCN + NO */
    /*eqcon[257] *= 1;  */

    /*reaction 259: HCNN + O2 <=> O + HCO + N2 */
    eqcon[258] *= 1e-06; 

    /*reaction 260: HCNN + OH <=> H + HCO + N2 */
    eqcon[259] *= 1e-06; 

    /*reaction 261: HCNN + H <=> CH2 + N2 */
    /*eqcon[260] *= 1;  */

    /*reaction 262: HNCO + O <=> NH + CO2 */
    /*eqcon[261] *= 1;  */

    /*reaction 263: HNCO + O <=> HNO + CO */
    /*eqcon[262] *= 1;  */

    /*reaction 264: HNCO + O <=> NCO + OH */
    /*eqcon[263] *= 1;  */

    /*reaction 265: HNCO + H <=> NH2 + CO */
    /*eqcon[264] *= 1;  */

    /*reaction 266: HNCO + H <=> H2 + NCO */
    /*eqcon[265] *= 1;  */

    /*reaction 267: HNCO + OH <=> NCO + H2O */
    /*eqcon[266] *= 1;  */

    /*reaction 268: HNCO + OH <=> NH2 + CO2 */
    /*eqcon[267] *= 1;  */

    /*reaction 269: HNCO + M <=> NH + CO + M */
    eqcon[268] *= 1e-06; 

    /*reaction 270: HCNO + H <=> H + HNCO */
    /*eqcon[269] *= 1;  */

    /*reaction 271: HCNO + H <=> OH + HCN */
    /*eqcon[270] *= 1;  */

    /*reaction 272: HCNO + H <=> NH2 + CO */
    /*eqcon[271] *= 1;  */

    /*reaction 273: HOCN + H <=> H + HNCO */
    /*eqcon[272] *= 1;  */

    /*reaction 274: HCCO + NO <=> HCNO + CO */
    /*eqcon[273] *= 1;  */

    /*reaction 275: CH3 + N <=> H2CN + H */
    /*eqcon[274] *= 1;  */

    /*reaction 276: CH3 + N <=> HCN + H2 */
    /*eqcon[275] *= 1;  */

    /*reaction 277: NH3 + H <=> NH2 + H2 */
    /*eqcon[276] *= 1;  */

    /*reaction 278: NH3 + OH <=> NH2 + H2O */
    /*eqcon[277] *= 1;  */

    /*reaction 279: NH3 + O <=> NH2 + OH */
    /*eqcon[278] *= 1;  */

    /*reaction 280: NH + CO2 <=> HNO + CO */
    /*eqcon[279] *= 1;  */

    /*reaction 281: CN + NO2 <=> NCO + NO */
    /*eqcon[280] *= 1;  */

    /*reaction 282: NCO + NO2 <=> N2O + CO2 */
    /*eqcon[281] *= 1;  */

    /*reaction 283: N + CO2 <=> NO + CO */
    /*eqcon[282] *= 1;  */

    /*reaction 284: O + CH3 => H + H2 + CO */
    eqcon[283] *= 1e-06; 

    /*reaction 285: O + C2H4 <=> H + CH2CHO */
    /*eqcon[284] *= 1;  */

    /*reaction 286: O + C2H5 <=> H + CH3CHO */
    /*eqcon[285] *= 1;  */

    /*reaction 287: OH + HO2 <=> O2 + H2O */
    /*eqcon[286] *= 1;  */

    /*reaction 288: OH + CH3 => H2 + CH2O */
    /*eqcon[287] *= 1;  */

    /*reaction 289: CH + H2 (+M) <=> CH3 (+M) */
    eqcon[288] *= 1e+06; 

    /*reaction 290: CH2 + O2 => 2 H + CO2 */
    eqcon[289] *= 1e-06; 

    /*reaction 291: CH2 + O2 <=> O + CH2O */
    /*eqcon[290] *= 1;  */

    /*reaction 292: CH2 + CH2 => 2 H + C2H2 */
    eqcon[291] *= 1e-06; 

    /*reaction 293: CH2(S) + H2O => H2 + CH2O */
    /*eqcon[292] *= 1;  */

    /*reaction 294: C2H3 + O2 <=> O + CH2CHO */
    /*eqcon[293] *= 1;  */

    /*reaction 295: C2H3 + O2 <=> HO2 + C2H2 */
    /*eqcon[294] *= 1;  */

    /*reaction 296: O + CH3CHO <=> OH + CH2CHO */
    /*eqcon[295] *= 1;  */

    /*reaction 297: O + CH3CHO => OH + CH3 + CO */
    eqcon[296] *= 1e-06; 

    /*reaction 298: O2 + CH3CHO => HO2 + CH3 + CO */
    eqcon[297] *= 1e-06; 

    /*reaction 299: H + CH3CHO <=> CH2CHO + H2 */
    /*eqcon[298] *= 1;  */

    /*reaction 300: H + CH3CHO => CH3 + H2 + CO */
    eqcon[299] *= 1e-06; 

    /*reaction 301: OH + CH3CHO => CH3 + H2O + CO */
    eqcon[300] *= 1e-06; 

    /*reaction 302: HO2 + CH3CHO => CH3 + H2O2 + CO */
    eqcon[301] *= 1e-06; 

    /*reaction 303: CH3 + CH3CHO => CH3 + CH4 + CO */
    eqcon[302] *= 1e-06; 

    /*reaction 304: H + CH2CO (+M) <=> CH2CHO (+M) */
    eqcon[303] *= 1e+06; 

    /*reaction 305: O + CH2CHO => H + CH2 + CO2 */
    eqcon[304] *= 1e-06; 

    /*reaction 306: O2 + CH2CHO => OH + CO + CH2O */
    eqcon[305] *= 1e-06; 

    /*reaction 307: O2 + CH2CHO => OH + 2 HCO */
    eqcon[306] *= 1e-06; 

    /*reaction 308: H + CH2CHO <=> CH3 + HCO */
    /*eqcon[307] *= 1;  */

    /*reaction 309: H + CH2CHO <=> CH2CO + H2 */
    /*eqcon[308] *= 1;  */

    /*reaction 310: OH + CH2CHO <=> H2O + CH2CO */
    /*eqcon[309] *= 1;  */

    /*reaction 311: OH + CH2CHO <=> HCO + CH2OH */
    /*eqcon[310] *= 1;  */

    /*reaction 312: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    eqcon[311] *= 1e+06; 

    /*reaction 313: O + C3H8 <=> OH + C3H7 */
    /*eqcon[312] *= 1;  */

    /*reaction 314: H + C3H8 <=> C3H7 + H2 */
    /*eqcon[313] *= 1;  */

    /*reaction 315: OH + C3H8 <=> C3H7 + H2O */
    /*eqcon[314] *= 1;  */

    /*reaction 316: C3H7 + H2O2 <=> HO2 + C3H8 */
    /*eqcon[315] *= 1;  */

    /*reaction 317: CH3 + C3H8 <=> C3H7 + CH4 */
    /*eqcon[316] *= 1;  */

    /*reaction 318: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    eqcon[317] *= 1e+06; 

    /*reaction 319: O + C3H7 <=> C2H5 + CH2O */
    /*eqcon[318] *= 1;  */

    /*reaction 320: H + C3H7 (+M) <=> C3H8 (+M) */
    eqcon[319] *= 1e+06; 

    /*reaction 321: H + C3H7 <=> CH3 + C2H5 */
    /*eqcon[320] *= 1;  */

    /*reaction 322: OH + C3H7 <=> C2H5 + CH2OH */
    /*eqcon[321] *= 1;  */

    /*reaction 323: HO2 + C3H7 <=> O2 + C3H8 */
    /*eqcon[322] *= 1;  */

    /*reaction 324: HO2 + C3H7 => OH + C2H5 + CH2O */
    eqcon[323] *= 1e-06; 

    /*reaction 325: CH3 + C3H7 <=> 2 C2H5 */
    /*eqcon[324] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[53]; /* temporary storage */

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

    /*reaction 12: O + CO (+M) <=> CO2 (+M) */
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

    /*reaction 135: CH2 + O2 => OH + H + CO */
    eqcon[134] *= 1e-06; 

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

    /*reaction 178: N + NO <=> N2 + O */
    /*eqcon[177] *= 1;  */

    /*reaction 179: N + O2 <=> NO + O */
    /*eqcon[178] *= 1;  */

    /*reaction 180: N + OH <=> NO + H */
    /*eqcon[179] *= 1;  */

    /*reaction 181: N2O + O <=> N2 + O2 */
    /*eqcon[180] *= 1;  */

    /*reaction 182: N2O + O <=> 2 NO */
    /*eqcon[181] *= 1;  */

    /*reaction 183: N2O + H <=> N2 + OH */
    /*eqcon[182] *= 1;  */

    /*reaction 184: N2O + OH <=> N2 + HO2 */
    /*eqcon[183] *= 1;  */

    /*reaction 185: N2O (+M) <=> N2 + O (+M) */
    eqcon[184] *= 1e-06; 

    /*reaction 186: HO2 + NO <=> NO2 + OH */
    /*eqcon[185] *= 1;  */

    /*reaction 187: NO + O + M <=> NO2 + M */
    eqcon[186] *= 1e+06; 

    /*reaction 188: NO2 + O <=> NO + O2 */
    /*eqcon[187] *= 1;  */

    /*reaction 189: NO2 + H <=> NO + OH */
    /*eqcon[188] *= 1;  */

    /*reaction 190: NH + O <=> NO + H */
    /*eqcon[189] *= 1;  */

    /*reaction 191: NH + H <=> N + H2 */
    /*eqcon[190] *= 1;  */

    /*reaction 192: NH + OH <=> HNO + H */
    /*eqcon[191] *= 1;  */

    /*reaction 193: NH + OH <=> N + H2O */
    /*eqcon[192] *= 1;  */

    /*reaction 194: NH + O2 <=> HNO + O */
    /*eqcon[193] *= 1;  */

    /*reaction 195: NH + O2 <=> NO + OH */
    /*eqcon[194] *= 1;  */

    /*reaction 196: NH + N <=> N2 + H */
    /*eqcon[195] *= 1;  */

    /*reaction 197: NH + H2O <=> HNO + H2 */
    /*eqcon[196] *= 1;  */

    /*reaction 198: NH + NO <=> N2 + OH */
    /*eqcon[197] *= 1;  */

    /*reaction 199: NH + NO <=> N2O + H */
    /*eqcon[198] *= 1;  */

    /*reaction 200: NH2 + O <=> OH + NH */
    /*eqcon[199] *= 1;  */

    /*reaction 201: NH2 + O <=> H + HNO */
    /*eqcon[200] *= 1;  */

    /*reaction 202: NH2 + H <=> NH + H2 */
    /*eqcon[201] *= 1;  */

    /*reaction 203: NH2 + OH <=> NH + H2O */
    /*eqcon[202] *= 1;  */

    /*reaction 204: NNH <=> N2 + H */
    eqcon[203] *= 1e-06; 

    /*reaction 205: NNH + M <=> N2 + H + M */
    eqcon[204] *= 1e-06; 

    /*reaction 206: NNH + O2 <=> HO2 + N2 */
    /*eqcon[205] *= 1;  */

    /*reaction 207: NNH + O <=> OH + N2 */
    /*eqcon[206] *= 1;  */

    /*reaction 208: NNH + O <=> NH + NO */
    /*eqcon[207] *= 1;  */

    /*reaction 209: NNH + H <=> H2 + N2 */
    /*eqcon[208] *= 1;  */

    /*reaction 210: NNH + OH <=> H2O + N2 */
    /*eqcon[209] *= 1;  */

    /*reaction 211: NNH + CH3 <=> CH4 + N2 */
    /*eqcon[210] *= 1;  */

    /*reaction 212: H + NO + M <=> HNO + M */
    eqcon[211] *= 1e+06; 

    /*reaction 213: HNO + O <=> NO + OH */
    /*eqcon[212] *= 1;  */

    /*reaction 214: HNO + H <=> H2 + NO */
    /*eqcon[213] *= 1;  */

    /*reaction 215: HNO + OH <=> NO + H2O */
    /*eqcon[214] *= 1;  */

    /*reaction 216: HNO + O2 <=> HO2 + NO */
    /*eqcon[215] *= 1;  */

    /*reaction 217: CN + O <=> CO + N */
    /*eqcon[216] *= 1;  */

    /*reaction 218: CN + OH <=> NCO + H */
    /*eqcon[217] *= 1;  */

    /*reaction 219: CN + H2O <=> HCN + OH */
    /*eqcon[218] *= 1;  */

    /*reaction 220: CN + O2 <=> NCO + O */
    /*eqcon[219] *= 1;  */

    /*reaction 221: CN + H2 <=> HCN + H */
    /*eqcon[220] *= 1;  */

    /*reaction 222: NCO + O <=> NO + CO */
    /*eqcon[221] *= 1;  */

    /*reaction 223: NCO + H <=> NH + CO */
    /*eqcon[222] *= 1;  */

    /*reaction 224: NCO + OH <=> NO + H + CO */
    eqcon[223] *= 1e-06; 

    /*reaction 225: NCO + N <=> N2 + CO */
    /*eqcon[224] *= 1;  */

    /*reaction 226: NCO + O2 <=> NO + CO2 */
    /*eqcon[225] *= 1;  */

    /*reaction 227: NCO + M <=> N + CO + M */
    eqcon[226] *= 1e-06; 

    /*reaction 228: NCO + NO <=> N2O + CO */
    /*eqcon[227] *= 1;  */

    /*reaction 229: NCO + NO <=> N2 + CO2 */
    /*eqcon[228] *= 1;  */

    /*reaction 230: HCN + M <=> H + CN + M */
    eqcon[229] *= 1e-06; 

    /*reaction 231: HCN + O <=> NCO + H */
    /*eqcon[230] *= 1;  */

    /*reaction 232: HCN + O <=> NH + CO */
    /*eqcon[231] *= 1;  */

    /*reaction 233: HCN + O <=> CN + OH */
    /*eqcon[232] *= 1;  */

    /*reaction 234: HCN + OH <=> HOCN + H */
    /*eqcon[233] *= 1;  */

    /*reaction 235: HCN + OH <=> HNCO + H */
    /*eqcon[234] *= 1;  */

    /*reaction 236: HCN + OH <=> NH2 + CO */
    /*eqcon[235] *= 1;  */

    /*reaction 237: H + HCN (+M) <=> H2CN (+M) */
    eqcon[236] *= 1e+06; 

    /*reaction 238: H2CN + N <=> N2 + CH2 */
    /*eqcon[237] *= 1;  */

    /*reaction 239: C + N2 <=> CN + N */
    /*eqcon[238] *= 1;  */

    /*reaction 240: CH + N2 <=> HCN + N */
    /*eqcon[239] *= 1;  */

    /*reaction 241: CH + N2 (+M) <=> HCNN (+M) */
    eqcon[240] *= 1e+06; 

    /*reaction 242: CH2 + N2 <=> HCN + NH */
    /*eqcon[241] *= 1;  */

    /*reaction 243: CH2(S) + N2 <=> NH + HCN */
    /*eqcon[242] *= 1;  */

    /*reaction 244: C + NO <=> CN + O */
    /*eqcon[243] *= 1;  */

    /*reaction 245: C + NO <=> CO + N */
    /*eqcon[244] *= 1;  */

    /*reaction 246: CH + NO <=> HCN + O */
    /*eqcon[245] *= 1;  */

    /*reaction 247: CH + NO <=> H + NCO */
    /*eqcon[246] *= 1;  */

    /*reaction 248: CH + NO <=> N + HCO */
    /*eqcon[247] *= 1;  */

    /*reaction 249: CH2 + NO <=> H + HNCO */
    /*eqcon[248] *= 1;  */

    /*reaction 250: CH2 + NO <=> OH + HCN */
    /*eqcon[249] *= 1;  */

    /*reaction 251: CH2 + NO <=> H + HCNO */
    /*eqcon[250] *= 1;  */

    /*reaction 252: CH2(S) + NO <=> H + HNCO */
    /*eqcon[251] *= 1;  */

    /*reaction 253: CH2(S) + NO <=> OH + HCN */
    /*eqcon[252] *= 1;  */

    /*reaction 254: CH2(S) + NO <=> H + HCNO */
    /*eqcon[253] *= 1;  */

    /*reaction 255: CH3 + NO <=> HCN + H2O */
    /*eqcon[254] *= 1;  */

    /*reaction 256: CH3 + NO <=> H2CN + OH */
    /*eqcon[255] *= 1;  */

    /*reaction 257: HCNN + O <=> CO + H + N2 */
    eqcon[256] *= 1e-06; 

    /*reaction 258: HCNN + O <=> HCN + NO */
    /*eqcon[257] *= 1;  */

    /*reaction 259: HCNN + O2 <=> O + HCO + N2 */
    eqcon[258] *= 1e-06; 

    /*reaction 260: HCNN + OH <=> H + HCO + N2 */
    eqcon[259] *= 1e-06; 

    /*reaction 261: HCNN + H <=> CH2 + N2 */
    /*eqcon[260] *= 1;  */

    /*reaction 262: HNCO + O <=> NH + CO2 */
    /*eqcon[261] *= 1;  */

    /*reaction 263: HNCO + O <=> HNO + CO */
    /*eqcon[262] *= 1;  */

    /*reaction 264: HNCO + O <=> NCO + OH */
    /*eqcon[263] *= 1;  */

    /*reaction 265: HNCO + H <=> NH2 + CO */
    /*eqcon[264] *= 1;  */

    /*reaction 266: HNCO + H <=> H2 + NCO */
    /*eqcon[265] *= 1;  */

    /*reaction 267: HNCO + OH <=> NCO + H2O */
    /*eqcon[266] *= 1;  */

    /*reaction 268: HNCO + OH <=> NH2 + CO2 */
    /*eqcon[267] *= 1;  */

    /*reaction 269: HNCO + M <=> NH + CO + M */
    eqcon[268] *= 1e-06; 

    /*reaction 270: HCNO + H <=> H + HNCO */
    /*eqcon[269] *= 1;  */

    /*reaction 271: HCNO + H <=> OH + HCN */
    /*eqcon[270] *= 1;  */

    /*reaction 272: HCNO + H <=> NH2 + CO */
    /*eqcon[271] *= 1;  */

    /*reaction 273: HOCN + H <=> H + HNCO */
    /*eqcon[272] *= 1;  */

    /*reaction 274: HCCO + NO <=> HCNO + CO */
    /*eqcon[273] *= 1;  */

    /*reaction 275: CH3 + N <=> H2CN + H */
    /*eqcon[274] *= 1;  */

    /*reaction 276: CH3 + N <=> HCN + H2 */
    /*eqcon[275] *= 1;  */

    /*reaction 277: NH3 + H <=> NH2 + H2 */
    /*eqcon[276] *= 1;  */

    /*reaction 278: NH3 + OH <=> NH2 + H2O */
    /*eqcon[277] *= 1;  */

    /*reaction 279: NH3 + O <=> NH2 + OH */
    /*eqcon[278] *= 1;  */

    /*reaction 280: NH + CO2 <=> HNO + CO */
    /*eqcon[279] *= 1;  */

    /*reaction 281: CN + NO2 <=> NCO + NO */
    /*eqcon[280] *= 1;  */

    /*reaction 282: NCO + NO2 <=> N2O + CO2 */
    /*eqcon[281] *= 1;  */

    /*reaction 283: N + CO2 <=> NO + CO */
    /*eqcon[282] *= 1;  */

    /*reaction 284: O + CH3 => H + H2 + CO */
    eqcon[283] *= 1e-06; 

    /*reaction 285: O + C2H4 <=> H + CH2CHO */
    /*eqcon[284] *= 1;  */

    /*reaction 286: O + C2H5 <=> H + CH3CHO */
    /*eqcon[285] *= 1;  */

    /*reaction 287: OH + HO2 <=> O2 + H2O */
    /*eqcon[286] *= 1;  */

    /*reaction 288: OH + CH3 => H2 + CH2O */
    /*eqcon[287] *= 1;  */

    /*reaction 289: CH + H2 (+M) <=> CH3 (+M) */
    eqcon[288] *= 1e+06; 

    /*reaction 290: CH2 + O2 => 2 H + CO2 */
    eqcon[289] *= 1e-06; 

    /*reaction 291: CH2 + O2 <=> O + CH2O */
    /*eqcon[290] *= 1;  */

    /*reaction 292: CH2 + CH2 => 2 H + C2H2 */
    eqcon[291] *= 1e-06; 

    /*reaction 293: CH2(S) + H2O => H2 + CH2O */
    /*eqcon[292] *= 1;  */

    /*reaction 294: C2H3 + O2 <=> O + CH2CHO */
    /*eqcon[293] *= 1;  */

    /*reaction 295: C2H3 + O2 <=> HO2 + C2H2 */
    /*eqcon[294] *= 1;  */

    /*reaction 296: O + CH3CHO <=> OH + CH2CHO */
    /*eqcon[295] *= 1;  */

    /*reaction 297: O + CH3CHO => OH + CH3 + CO */
    eqcon[296] *= 1e-06; 

    /*reaction 298: O2 + CH3CHO => HO2 + CH3 + CO */
    eqcon[297] *= 1e-06; 

    /*reaction 299: H + CH3CHO <=> CH2CHO + H2 */
    /*eqcon[298] *= 1;  */

    /*reaction 300: H + CH3CHO => CH3 + H2 + CO */
    eqcon[299] *= 1e-06; 

    /*reaction 301: OH + CH3CHO => CH3 + H2O + CO */
    eqcon[300] *= 1e-06; 

    /*reaction 302: HO2 + CH3CHO => CH3 + H2O2 + CO */
    eqcon[301] *= 1e-06; 

    /*reaction 303: CH3 + CH3CHO => CH3 + CH4 + CO */
    eqcon[302] *= 1e-06; 

    /*reaction 304: H + CH2CO (+M) <=> CH2CHO (+M) */
    eqcon[303] *= 1e+06; 

    /*reaction 305: O + CH2CHO => H + CH2 + CO2 */
    eqcon[304] *= 1e-06; 

    /*reaction 306: O2 + CH2CHO => OH + CO + CH2O */
    eqcon[305] *= 1e-06; 

    /*reaction 307: O2 + CH2CHO => OH + 2 HCO */
    eqcon[306] *= 1e-06; 

    /*reaction 308: H + CH2CHO <=> CH3 + HCO */
    /*eqcon[307] *= 1;  */

    /*reaction 309: H + CH2CHO <=> CH2CO + H2 */
    /*eqcon[308] *= 1;  */

    /*reaction 310: OH + CH2CHO <=> H2O + CH2CO */
    /*eqcon[309] *= 1;  */

    /*reaction 311: OH + CH2CHO <=> HCO + CH2OH */
    /*eqcon[310] *= 1;  */

    /*reaction 312: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    eqcon[311] *= 1e+06; 

    /*reaction 313: O + C3H8 <=> OH + C3H7 */
    /*eqcon[312] *= 1;  */

    /*reaction 314: H + C3H8 <=> C3H7 + H2 */
    /*eqcon[313] *= 1;  */

    /*reaction 315: OH + C3H8 <=> C3H7 + H2O */
    /*eqcon[314] *= 1;  */

    /*reaction 316: C3H7 + H2O2 <=> HO2 + C3H8 */
    /*eqcon[315] *= 1;  */

    /*reaction 317: CH3 + C3H8 <=> C3H7 + CH4 */
    /*eqcon[316] *= 1;  */

    /*reaction 318: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    eqcon[317] *= 1e+06; 

    /*reaction 319: O + C3H7 <=> C2H5 + CH2O */
    /*eqcon[318] *= 1;  */

    /*reaction 320: H + C3H7 (+M) <=> C3H8 (+M) */
    eqcon[319] *= 1e+06; 

    /*reaction 321: H + C3H7 <=> CH3 + C2H5 */
    /*eqcon[320] *= 1;  */

    /*reaction 322: OH + C3H7 <=> C2H5 + CH2OH */
    /*eqcon[321] *= 1;  */

    /*reaction 323: HO2 + C3H7 <=> O2 + C3H8 */
    /*eqcon[322] *= 1;  */

    /*reaction 324: HO2 + C3H7 => OH + C2H5 + CH2O */
    eqcon[323] *= 1e-06; 

    /*reaction 325: CH3 + C3H7 <=> 2 C2H5 */
    /*eqcon[324] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[53]; /* temporary storage */

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

    /*reaction 12: O + CO (+M) <=> CO2 (+M) */
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

    /*reaction 135: CH2 + O2 => OH + H + CO */
    eqcon[134] *= 1e-06; 

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

    /*reaction 178: N + NO <=> N2 + O */
    /*eqcon[177] *= 1;  */

    /*reaction 179: N + O2 <=> NO + O */
    /*eqcon[178] *= 1;  */

    /*reaction 180: N + OH <=> NO + H */
    /*eqcon[179] *= 1;  */

    /*reaction 181: N2O + O <=> N2 + O2 */
    /*eqcon[180] *= 1;  */

    /*reaction 182: N2O + O <=> 2 NO */
    /*eqcon[181] *= 1;  */

    /*reaction 183: N2O + H <=> N2 + OH */
    /*eqcon[182] *= 1;  */

    /*reaction 184: N2O + OH <=> N2 + HO2 */
    /*eqcon[183] *= 1;  */

    /*reaction 185: N2O (+M) <=> N2 + O (+M) */
    eqcon[184] *= 1e-06; 

    /*reaction 186: HO2 + NO <=> NO2 + OH */
    /*eqcon[185] *= 1;  */

    /*reaction 187: NO + O + M <=> NO2 + M */
    eqcon[186] *= 1e+06; 

    /*reaction 188: NO2 + O <=> NO + O2 */
    /*eqcon[187] *= 1;  */

    /*reaction 189: NO2 + H <=> NO + OH */
    /*eqcon[188] *= 1;  */

    /*reaction 190: NH + O <=> NO + H */
    /*eqcon[189] *= 1;  */

    /*reaction 191: NH + H <=> N + H2 */
    /*eqcon[190] *= 1;  */

    /*reaction 192: NH + OH <=> HNO + H */
    /*eqcon[191] *= 1;  */

    /*reaction 193: NH + OH <=> N + H2O */
    /*eqcon[192] *= 1;  */

    /*reaction 194: NH + O2 <=> HNO + O */
    /*eqcon[193] *= 1;  */

    /*reaction 195: NH + O2 <=> NO + OH */
    /*eqcon[194] *= 1;  */

    /*reaction 196: NH + N <=> N2 + H */
    /*eqcon[195] *= 1;  */

    /*reaction 197: NH + H2O <=> HNO + H2 */
    /*eqcon[196] *= 1;  */

    /*reaction 198: NH + NO <=> N2 + OH */
    /*eqcon[197] *= 1;  */

    /*reaction 199: NH + NO <=> N2O + H */
    /*eqcon[198] *= 1;  */

    /*reaction 200: NH2 + O <=> OH + NH */
    /*eqcon[199] *= 1;  */

    /*reaction 201: NH2 + O <=> H + HNO */
    /*eqcon[200] *= 1;  */

    /*reaction 202: NH2 + H <=> NH + H2 */
    /*eqcon[201] *= 1;  */

    /*reaction 203: NH2 + OH <=> NH + H2O */
    /*eqcon[202] *= 1;  */

    /*reaction 204: NNH <=> N2 + H */
    eqcon[203] *= 1e-06; 

    /*reaction 205: NNH + M <=> N2 + H + M */
    eqcon[204] *= 1e-06; 

    /*reaction 206: NNH + O2 <=> HO2 + N2 */
    /*eqcon[205] *= 1;  */

    /*reaction 207: NNH + O <=> OH + N2 */
    /*eqcon[206] *= 1;  */

    /*reaction 208: NNH + O <=> NH + NO */
    /*eqcon[207] *= 1;  */

    /*reaction 209: NNH + H <=> H2 + N2 */
    /*eqcon[208] *= 1;  */

    /*reaction 210: NNH + OH <=> H2O + N2 */
    /*eqcon[209] *= 1;  */

    /*reaction 211: NNH + CH3 <=> CH4 + N2 */
    /*eqcon[210] *= 1;  */

    /*reaction 212: H + NO + M <=> HNO + M */
    eqcon[211] *= 1e+06; 

    /*reaction 213: HNO + O <=> NO + OH */
    /*eqcon[212] *= 1;  */

    /*reaction 214: HNO + H <=> H2 + NO */
    /*eqcon[213] *= 1;  */

    /*reaction 215: HNO + OH <=> NO + H2O */
    /*eqcon[214] *= 1;  */

    /*reaction 216: HNO + O2 <=> HO2 + NO */
    /*eqcon[215] *= 1;  */

    /*reaction 217: CN + O <=> CO + N */
    /*eqcon[216] *= 1;  */

    /*reaction 218: CN + OH <=> NCO + H */
    /*eqcon[217] *= 1;  */

    /*reaction 219: CN + H2O <=> HCN + OH */
    /*eqcon[218] *= 1;  */

    /*reaction 220: CN + O2 <=> NCO + O */
    /*eqcon[219] *= 1;  */

    /*reaction 221: CN + H2 <=> HCN + H */
    /*eqcon[220] *= 1;  */

    /*reaction 222: NCO + O <=> NO + CO */
    /*eqcon[221] *= 1;  */

    /*reaction 223: NCO + H <=> NH + CO */
    /*eqcon[222] *= 1;  */

    /*reaction 224: NCO + OH <=> NO + H + CO */
    eqcon[223] *= 1e-06; 

    /*reaction 225: NCO + N <=> N2 + CO */
    /*eqcon[224] *= 1;  */

    /*reaction 226: NCO + O2 <=> NO + CO2 */
    /*eqcon[225] *= 1;  */

    /*reaction 227: NCO + M <=> N + CO + M */
    eqcon[226] *= 1e-06; 

    /*reaction 228: NCO + NO <=> N2O + CO */
    /*eqcon[227] *= 1;  */

    /*reaction 229: NCO + NO <=> N2 + CO2 */
    /*eqcon[228] *= 1;  */

    /*reaction 230: HCN + M <=> H + CN + M */
    eqcon[229] *= 1e-06; 

    /*reaction 231: HCN + O <=> NCO + H */
    /*eqcon[230] *= 1;  */

    /*reaction 232: HCN + O <=> NH + CO */
    /*eqcon[231] *= 1;  */

    /*reaction 233: HCN + O <=> CN + OH */
    /*eqcon[232] *= 1;  */

    /*reaction 234: HCN + OH <=> HOCN + H */
    /*eqcon[233] *= 1;  */

    /*reaction 235: HCN + OH <=> HNCO + H */
    /*eqcon[234] *= 1;  */

    /*reaction 236: HCN + OH <=> NH2 + CO */
    /*eqcon[235] *= 1;  */

    /*reaction 237: H + HCN (+M) <=> H2CN (+M) */
    eqcon[236] *= 1e+06; 

    /*reaction 238: H2CN + N <=> N2 + CH2 */
    /*eqcon[237] *= 1;  */

    /*reaction 239: C + N2 <=> CN + N */
    /*eqcon[238] *= 1;  */

    /*reaction 240: CH + N2 <=> HCN + N */
    /*eqcon[239] *= 1;  */

    /*reaction 241: CH + N2 (+M) <=> HCNN (+M) */
    eqcon[240] *= 1e+06; 

    /*reaction 242: CH2 + N2 <=> HCN + NH */
    /*eqcon[241] *= 1;  */

    /*reaction 243: CH2(S) + N2 <=> NH + HCN */
    /*eqcon[242] *= 1;  */

    /*reaction 244: C + NO <=> CN + O */
    /*eqcon[243] *= 1;  */

    /*reaction 245: C + NO <=> CO + N */
    /*eqcon[244] *= 1;  */

    /*reaction 246: CH + NO <=> HCN + O */
    /*eqcon[245] *= 1;  */

    /*reaction 247: CH + NO <=> H + NCO */
    /*eqcon[246] *= 1;  */

    /*reaction 248: CH + NO <=> N + HCO */
    /*eqcon[247] *= 1;  */

    /*reaction 249: CH2 + NO <=> H + HNCO */
    /*eqcon[248] *= 1;  */

    /*reaction 250: CH2 + NO <=> OH + HCN */
    /*eqcon[249] *= 1;  */

    /*reaction 251: CH2 + NO <=> H + HCNO */
    /*eqcon[250] *= 1;  */

    /*reaction 252: CH2(S) + NO <=> H + HNCO */
    /*eqcon[251] *= 1;  */

    /*reaction 253: CH2(S) + NO <=> OH + HCN */
    /*eqcon[252] *= 1;  */

    /*reaction 254: CH2(S) + NO <=> H + HCNO */
    /*eqcon[253] *= 1;  */

    /*reaction 255: CH3 + NO <=> HCN + H2O */
    /*eqcon[254] *= 1;  */

    /*reaction 256: CH3 + NO <=> H2CN + OH */
    /*eqcon[255] *= 1;  */

    /*reaction 257: HCNN + O <=> CO + H + N2 */
    eqcon[256] *= 1e-06; 

    /*reaction 258: HCNN + O <=> HCN + NO */
    /*eqcon[257] *= 1;  */

    /*reaction 259: HCNN + O2 <=> O + HCO + N2 */
    eqcon[258] *= 1e-06; 

    /*reaction 260: HCNN + OH <=> H + HCO + N2 */
    eqcon[259] *= 1e-06; 

    /*reaction 261: HCNN + H <=> CH2 + N2 */
    /*eqcon[260] *= 1;  */

    /*reaction 262: HNCO + O <=> NH + CO2 */
    /*eqcon[261] *= 1;  */

    /*reaction 263: HNCO + O <=> HNO + CO */
    /*eqcon[262] *= 1;  */

    /*reaction 264: HNCO + O <=> NCO + OH */
    /*eqcon[263] *= 1;  */

    /*reaction 265: HNCO + H <=> NH2 + CO */
    /*eqcon[264] *= 1;  */

    /*reaction 266: HNCO + H <=> H2 + NCO */
    /*eqcon[265] *= 1;  */

    /*reaction 267: HNCO + OH <=> NCO + H2O */
    /*eqcon[266] *= 1;  */

    /*reaction 268: HNCO + OH <=> NH2 + CO2 */
    /*eqcon[267] *= 1;  */

    /*reaction 269: HNCO + M <=> NH + CO + M */
    eqcon[268] *= 1e-06; 

    /*reaction 270: HCNO + H <=> H + HNCO */
    /*eqcon[269] *= 1;  */

    /*reaction 271: HCNO + H <=> OH + HCN */
    /*eqcon[270] *= 1;  */

    /*reaction 272: HCNO + H <=> NH2 + CO */
    /*eqcon[271] *= 1;  */

    /*reaction 273: HOCN + H <=> H + HNCO */
    /*eqcon[272] *= 1;  */

    /*reaction 274: HCCO + NO <=> HCNO + CO */
    /*eqcon[273] *= 1;  */

    /*reaction 275: CH3 + N <=> H2CN + H */
    /*eqcon[274] *= 1;  */

    /*reaction 276: CH3 + N <=> HCN + H2 */
    /*eqcon[275] *= 1;  */

    /*reaction 277: NH3 + H <=> NH2 + H2 */
    /*eqcon[276] *= 1;  */

    /*reaction 278: NH3 + OH <=> NH2 + H2O */
    /*eqcon[277] *= 1;  */

    /*reaction 279: NH3 + O <=> NH2 + OH */
    /*eqcon[278] *= 1;  */

    /*reaction 280: NH + CO2 <=> HNO + CO */
    /*eqcon[279] *= 1;  */

    /*reaction 281: CN + NO2 <=> NCO + NO */
    /*eqcon[280] *= 1;  */

    /*reaction 282: NCO + NO2 <=> N2O + CO2 */
    /*eqcon[281] *= 1;  */

    /*reaction 283: N + CO2 <=> NO + CO */
    /*eqcon[282] *= 1;  */

    /*reaction 284: O + CH3 => H + H2 + CO */
    eqcon[283] *= 1e-06; 

    /*reaction 285: O + C2H4 <=> H + CH2CHO */
    /*eqcon[284] *= 1;  */

    /*reaction 286: O + C2H5 <=> H + CH3CHO */
    /*eqcon[285] *= 1;  */

    /*reaction 287: OH + HO2 <=> O2 + H2O */
    /*eqcon[286] *= 1;  */

    /*reaction 288: OH + CH3 => H2 + CH2O */
    /*eqcon[287] *= 1;  */

    /*reaction 289: CH + H2 (+M) <=> CH3 (+M) */
    eqcon[288] *= 1e+06; 

    /*reaction 290: CH2 + O2 => 2 H + CO2 */
    eqcon[289] *= 1e-06; 

    /*reaction 291: CH2 + O2 <=> O + CH2O */
    /*eqcon[290] *= 1;  */

    /*reaction 292: CH2 + CH2 => 2 H + C2H2 */
    eqcon[291] *= 1e-06; 

    /*reaction 293: CH2(S) + H2O => H2 + CH2O */
    /*eqcon[292] *= 1;  */

    /*reaction 294: C2H3 + O2 <=> O + CH2CHO */
    /*eqcon[293] *= 1;  */

    /*reaction 295: C2H3 + O2 <=> HO2 + C2H2 */
    /*eqcon[294] *= 1;  */

    /*reaction 296: O + CH3CHO <=> OH + CH2CHO */
    /*eqcon[295] *= 1;  */

    /*reaction 297: O + CH3CHO => OH + CH3 + CO */
    eqcon[296] *= 1e-06; 

    /*reaction 298: O2 + CH3CHO => HO2 + CH3 + CO */
    eqcon[297] *= 1e-06; 

    /*reaction 299: H + CH3CHO <=> CH2CHO + H2 */
    /*eqcon[298] *= 1;  */

    /*reaction 300: H + CH3CHO => CH3 + H2 + CO */
    eqcon[299] *= 1e-06; 

    /*reaction 301: OH + CH3CHO => CH3 + H2O + CO */
    eqcon[300] *= 1e-06; 

    /*reaction 302: HO2 + CH3CHO => CH3 + H2O2 + CO */
    eqcon[301] *= 1e-06; 

    /*reaction 303: CH3 + CH3CHO => CH3 + CH4 + CO */
    eqcon[302] *= 1e-06; 

    /*reaction 304: H + CH2CO (+M) <=> CH2CHO (+M) */
    eqcon[303] *= 1e+06; 

    /*reaction 305: O + CH2CHO => H + CH2 + CO2 */
    eqcon[304] *= 1e-06; 

    /*reaction 306: O2 + CH2CHO => OH + CO + CH2O */
    eqcon[305] *= 1e-06; 

    /*reaction 307: O2 + CH2CHO => OH + 2 HCO */
    eqcon[306] *= 1e-06; 

    /*reaction 308: H + CH2CHO <=> CH3 + HCO */
    /*eqcon[307] *= 1;  */

    /*reaction 309: H + CH2CHO <=> CH2CO + H2 */
    /*eqcon[308] *= 1;  */

    /*reaction 310: OH + CH2CHO <=> H2O + CH2CO */
    /*eqcon[309] *= 1;  */

    /*reaction 311: OH + CH2CHO <=> HCO + CH2OH */
    /*eqcon[310] *= 1;  */

    /*reaction 312: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    eqcon[311] *= 1e+06; 

    /*reaction 313: O + C3H8 <=> OH + C3H7 */
    /*eqcon[312] *= 1;  */

    /*reaction 314: H + C3H8 <=> C3H7 + H2 */
    /*eqcon[313] *= 1;  */

    /*reaction 315: OH + C3H8 <=> C3H7 + H2O */
    /*eqcon[314] *= 1;  */

    /*reaction 316: C3H7 + H2O2 <=> HO2 + C3H8 */
    /*eqcon[315] *= 1;  */

    /*reaction 317: CH3 + C3H8 <=> C3H7 + CH4 */
    /*eqcon[316] *= 1;  */

    /*reaction 318: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    eqcon[317] *= 1e+06; 

    /*reaction 319: O + C3H7 <=> C2H5 + CH2O */
    /*eqcon[318] *= 1;  */

    /*reaction 320: H + C3H7 (+M) <=> C3H8 (+M) */
    eqcon[319] *= 1e+06; 

    /*reaction 321: H + C3H7 <=> CH3 + C2H5 */
    /*eqcon[320] *= 1;  */

    /*reaction 322: OH + C3H7 <=> C2H5 + CH2OH */
    /*eqcon[321] *= 1;  */

    /*reaction 323: HO2 + C3H7 <=> O2 + C3H8 */
    /*eqcon[322] *= 1;  */

    /*reaction 324: HO2 + C3H7 => OH + C2H5 + CH2O */
    eqcon[323] *= 1e-06; 

    /*reaction 325: CH3 + C3H7 <=> 2 C2H5 */
    /*eqcon[324] *= 1;  */
}


/*compute the production rate for each species */
void productionRate(double * wdot, double * sc, double T)
{
    double qdot;

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[53];                /*Gibbs free energy */
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
    double refC = 101325 / 8.31451 / T;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 53; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*zero out wdot */
    for (id = 0; id < 53; ++id) {
        wdot[id] = 0.0;
    }

    /*reaction 1: 2 O + M <=> O2 + M */
    phi_f = sc[2]*sc[2];
    alpha = mixture + 1.4*sc[0] + 14.4*sc[5] + sc[13] + 0.75*sc[14] + 2.6*sc[15] + 2*sc[26] + -0.17*sc[48];
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
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
    k_f = 1e-06 * 38700*exp(2.7*tc[0]-3150.48/tc[1]);
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
    k_f = 1e-06 * 5.06e+13;
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

    /*reaction 12: O + CO (+M) <=> CO2 (+M) */
    phi_f = sc[2]*sc[14];
    alpha = mixture + sc[0] + 5*sc[3] + 5*sc[5] + sc[13] + 0.5*sc[14] + 2.5*sc[15] + 2*sc[26] + -0.5*sc[48];
    k_f = 1e-06 * 1.8e+10*exp(-1200.3/tc[1]);
    redP = 1e-12 * alpha / k_f * 6.02e+14*exp(-1509.81/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
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
    k_f = 1e-06 * 1.35e+07*exp(2*tc[0]-956.215/tc[1]);
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
    k_f = 1e-06 * 6.94e+06*exp(2*tc[0]-956.215/tc[1]);
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
    k_f = 1e-06 * 1.25e+07*exp(1.83*tc[0]-110.72/tc[1]);
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
    k_f = 1e-06 * 2.24e+13;
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
    alpha = mixture + -1*sc[3] + -1*sc[5] + -0.25*sc[14] + 0.5*sc[15] + 0.5*sc[26] + -1*sc[47] + -1*sc[48];
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
    k_f = 1e-12 * 2.08e+19*exp(-1.24*tc[0]);
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
    k_f = 1e-12 * 1.126e+19*exp(-0.76*tc[0]);
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
    phi_f = sc[1]*sc[3]*sc[47];
    k_f = 1e-12 * 2.6e+19*exp(-1.24*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[47];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[47]) - (g_RT[6] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[47] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[47] += 1 * qdot;

    /*reaction 37: H + O2 + AR <=> HO2 + AR */
    phi_f = sc[1]*sc[3]*sc[48];
    k_f = 1e-12 * 7e+17*exp(-0.8*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[48];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[48]) - (g_RT[6] + g_RT[48]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[48] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[48] += 1 * qdot;

    /*reaction 38: H + O2 <=> O + OH */
    phi_f = sc[1]*sc[3];
    k_f = 1e-06 * 2.65e+16*exp(-0.6707*tc[0]-8576.25/tc[1]);
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
    alpha = mixture + -1*sc[0] + -1*sc[5] + sc[13] + -1*sc[15] + 2*sc[26] + -0.37*sc[48];
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
    alpha = mixture + -0.27*sc[0] + 2.65*sc[5] + sc[13] + 2*sc[26] + -0.62*sc[48];
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
    k_f = 1e-06 * 4.48e+13*exp(-537.494/tc[1]);
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
    k_f = 1e-06 * 8.4e+13*exp(-319.577/tc[1]);
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
    k_f = 1e-06 * 1.65e+14;
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 6e+14;
    redP = 1e-12 * alpha / k_f * 1.04e+26*exp(-2.76*tc[0]-805.234/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.438*exp(T/-91))+ (0.562*exp(T/-5836))+ (exp(-8552/T)));
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
    alpha = mixture + sc[0] + 5*sc[5] + 2*sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 1.39e+16*exp(-0.534*tc[0]-269.753/tc[1]);
    redP = 1e-12 * alpha / k_f * 2.62e+33*exp(-4.76*tc[0]-1227.98/tc[1]);
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 1.09e+12*exp(0.48*tc[0]+130.851/tc[1]);
    redP = 1e-12 * alpha / k_f * 2.47e+24*exp(-2.57*tc[0]-213.89/tc[1]);
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
    k_f = 1e-06 * 5.74e+07*exp(1.9*tc[0]-1379.97/tc[1]);
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
    k_f = 1e-06 * 1.055e+12*exp(0.5*tc[0]-43.2813/tc[1]);
    redP = 1e-12 * alpha / k_f * 4.36e+31*exp(-4.65*tc[0]-2556.62/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.4*exp(T/-100))+ (0.6*exp(T/-90000))+ (exp(-10000/T)));
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
    k_f = 1e-06 * 1.65e+11*exp(0.65*tc[0]+142.929/tc[1]);
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
    k_f = 1e-06 * 3.28e+13*exp(-0.09*tc[0]-306.995/tc[1]);
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
    k_f = 1e-06 * 2.43e+12*exp(0.515*tc[0]-25.1636/tc[1]);
    redP = 1e-12 * alpha / k_f * 4.66e+41*exp(-7.44*tc[0]-7086.06/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.3*exp(T/-100))+ (0.7*exp(T/-90000))+ (exp(-10000/T)));
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
    k_f = 1e-06 * 4.15e+07*exp(1.63*tc[0]-968.294/tc[1]);
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
    k_f = 1e-06 * 1.5e+12*exp(0.5*tc[0]+55.3598/tc[1]);
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
    k_f = 1e-06 * 2.62e+14*exp(-0.23*tc[0]-538.5/tc[1]);
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 5.4e+11*exp(0.454*tc[0]-915.954/tc[1]);
    redP = 1e-12 * alpha / k_f * 6e+41*exp(-7.62*tc[0]-3507.8/tc[1]);
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
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
    k_f = 1e-06 * 1.45e+13*exp(+251.636/tc[1]);
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
    k_f = 1e-06 * 2e+12*exp(-214.897/tc[1]);
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
    k_f = 1e-06 * 1.7e+18*exp(-14801.2/tc[1]);
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
    k_f = 1e-06 * 2.79e+18*exp(-1.43*tc[0]-669.351/tc[1]);
    redP = 1e-12 * alpha / k_f * 4e+36*exp(-5.92*tc[0]-1580.27/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.588*exp(T/-195))+ (0.412*exp(T/-5900))+ (exp(-6394/T)));
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
    k_f = 1e-06 * 6.44e+17*exp(-1.34*tc[0]-713.135/tc[1]);
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
    k_f = 1e-06 * 3.78e+13;
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
    k_f = 1e-06 * 5.6e+06*exp(2*tc[0]-6039.26/tc[1]);
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
    k_f = 1e-06 * 6.71e+13;
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
    k_f = 1e-06 * 1.08e+14*exp(-1565.17/tc[1]);
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
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
    k_f = 1e-06 * 1.9e+14*exp(-7947.66/tc[1]);
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

    /*reaction 135: CH2 + O2 => OH + H + CO */
    phi_f = sc[10]*sc[3];
    k_f = 1e-06 * 5e+12*exp(-754.907/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[14] += 1 * qdot;

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
    k_f = 1e-06 * 1.6e+15*exp(-6011.07/tc[1]);
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
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
    phi_f = sc[11]*sc[47];
    k_f = 1e-06 * 1.5e+13*exp(-301.963/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[47];
    Kc = exp((g_RT[11] + g_RT[47]) - (g_RT[10] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[47] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[47] += 1 * qdot;

    /*reaction 143: CH2(S) + AR <=> CH2 + AR */
    phi_f = sc[11]*sc[48];
    k_f = 1e-06 * 9e+12*exp(-301.963/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[48];
    Kc = exp((g_RT[11] + g_RT[48]) - (g_RT[10] + g_RT[48]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[48] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[48] += 1 * qdot;

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
    k_f = 1e-06 * 4.82e+17*exp(-1.16*tc[0]-576.246/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.88e+38*exp(-6.36*tc[0]-2536.49/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.3973*exp(T/-208))+ (0.6027*exp(T/-3922))+ (exp(-10180/T)));
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
    k_f = 1e-06 * 3.56e+13*exp(-15339.7/tc[1]);
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
    k_f = 1e-06 * 2.31e+12*exp(-10224/tc[1]);
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 6.77e+16*exp(-1.18*tc[0]-329.139/tc[1]);
    redP = 1e-12 * alpha / k_f * 3.4e+41*exp(-7.03*tc[0]-1390.04/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.381*exp(T/-73.2))+ (0.619*exp(T/-1180))+ (exp(-9999/T)));
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
    k_f = 1e-06 * 6.84e+12*exp(0.1*tc[0]-5334.68/tc[1]);
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
    k_f = 1e-06 * 1.5e+18*exp(-1*tc[0]-8555.61/tc[1]);
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
    k_f = 1e-06 * 1.345e+13*exp(-201.309/tc[1]);
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
    k_f = 1e-06 * 1e+13*exp(+379.97/tc[1]);
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
    k_f = 1e-06 * 5.68e+10*exp(0.9*tc[0]-1003.02/tc[1]);
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
    k_f = 1e-06 * 4.58e+16*exp(-1.39*tc[0]-510.82/tc[1]);
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1 * 8e+12*exp(0.44*tc[0]-43668.8/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.58e+51*exp(-9.3*tc[0]-49219.9/tc[1]);
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
    k_f = 1e-06 * 3.2e+12*exp(-429.794/tc[1]);
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

    /*reaction 178: N + NO <=> N2 + O */
    phi_f = sc[30]*sc[35];
    k_f = 1e-06 * 2.7e+13*exp(-178.661/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[2];
    Kc = exp((g_RT[30] + g_RT[35]) - (g_RT[47] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[30] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[47] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 179: N + O2 <=> NO + O */
    phi_f = sc[30]*sc[3];
    k_f = 1e-06 * 9e+09*exp(1*tc[0]-3271.26/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[2];
    Kc = exp((g_RT[30] + g_RT[3]) - (g_RT[35] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[30] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 180: N + OH <=> NO + H */
    phi_f = sc[30]*sc[4];
    k_f = 1e-06 * 3.36e+13*exp(-193.759/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[1];
    Kc = exp((g_RT[30] + g_RT[4]) - (g_RT[35] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[30] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 181: N2O + O <=> N2 + O2 */
    phi_f = sc[37]*sc[2];
    k_f = 1e-06 * 1.4e+12*exp(-5440.36/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[3];
    Kc = exp((g_RT[37] + g_RT[2]) - (g_RT[47] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[37] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[47] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 182: N2O + O <=> 2 NO */
    phi_f = sc[37]*sc[2];
    k_f = 1e-06 * 2.9e+13*exp(-11650.7/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[35];
    Kc = exp((g_RT[37] + g_RT[2]) - (2 * g_RT[35]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[37] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[35] += 2 * qdot;

    /*reaction 183: N2O + H <=> N2 + OH */
    phi_f = sc[37]*sc[1];
    k_f = 1e-06 * 3.87e+14*exp(-9501.76/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[4];
    Kc = exp((g_RT[37] + g_RT[1]) - (g_RT[47] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[37] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[47] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 184: N2O + OH <=> N2 + HO2 */
    phi_f = sc[37]*sc[4];
    k_f = 1e-06 * 2e+12*exp(-10598.9/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[6];
    Kc = exp((g_RT[37] + g_RT[4]) - (g_RT[47] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[37] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[47] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 185: N2O (+M) <=> N2 + O (+M) */
    phi_f = sc[37];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.375*sc[48];
    k_f = 1 * 7.91e+10*exp(-28193.3/tc[1]);
    redP = 1e-12 * alpha / k_f * 6.37e+14*exp(-28505.3/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[2];
    Kc = refC * exp((g_RT[37]) - (g_RT[47] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[37] -= 1 * qdot;
    wdot[47] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 186: HO2 + NO <=> NO2 + OH */
    phi_f = sc[6]*sc[35];
    k_f = 1e-06 * 2.11e+12*exp(+241.57/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[36]*sc[4];
    Kc = exp((g_RT[6] + g_RT[35]) - (g_RT[36] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[36] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 187: NO + O + M <=> NO2 + M */
    phi_f = sc[35]*sc[2];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-12 * alpha * 1.06e+20*exp(-1.41*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[36];
    Kc = 1.0 / (refC) * exp((g_RT[35] + g_RT[2]) - (g_RT[36]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[35] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[36] += 1 * qdot;

    /*reaction 188: NO2 + O <=> NO + O2 */
    phi_f = sc[36]*sc[2];
    k_f = 1e-06 * 3.9e+12*exp(+120.785/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[3];
    Kc = exp((g_RT[36] + g_RT[2]) - (g_RT[35] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 189: NO2 + H <=> NO + OH */
    phi_f = sc[36]*sc[1];
    k_f = 1e-06 * 1.32e+14*exp(-181.178/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[4];
    Kc = exp((g_RT[36] + g_RT[1]) - (g_RT[35] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 190: NH + O <=> NO + H */
    phi_f = sc[31]*sc[2];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[1];
    Kc = exp((g_RT[31] + g_RT[2]) - (g_RT[35] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 191: NH + H <=> N + H2 */
    phi_f = sc[31]*sc[1];
    k_f = 1e-06 * 3.2e+13*exp(-166.08/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[0];
    Kc = exp((g_RT[31] + g_RT[1]) - (g_RT[30] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[30] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 192: NH + OH <=> HNO + H */
    phi_f = sc[31]*sc[4];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[1];
    Kc = exp((g_RT[31] + g_RT[4]) - (g_RT[38] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[38] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 193: NH + OH <=> N + H2O */
    phi_f = sc[31]*sc[4];
    k_f = 1e-06 * 2e+09*exp(1.2*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[5];
    Kc = exp((g_RT[31] + g_RT[4]) - (g_RT[30] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[30] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 194: NH + O2 <=> HNO + O */
    phi_f = sc[31]*sc[3];
    k_f = 1e-06 * 461000*exp(2*tc[0]-3271.26/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[2];
    Kc = exp((g_RT[31] + g_RT[3]) - (g_RT[38] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[38] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 195: NH + O2 <=> NO + OH */
    phi_f = sc[31]*sc[3];
    k_f = 1e-06 * 1.28e+06*exp(1.5*tc[0]-50.3271/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[4];
    Kc = exp((g_RT[31] + g_RT[3]) - (g_RT[35] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 196: NH + N <=> N2 + H */
    phi_f = sc[31]*sc[30];
    k_f = 1e-06 * 1.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[1];
    Kc = exp((g_RT[31] + g_RT[30]) - (g_RT[47] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[30] -= 1 * qdot;
    wdot[47] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 197: NH + H2O <=> HNO + H2 */
    phi_f = sc[31]*sc[5];
    k_f = 1e-06 * 2e+13*exp(-6970.31/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[0];
    Kc = exp((g_RT[31] + g_RT[5]) - (g_RT[38] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[38] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 198: NH + NO <=> N2 + OH */
    phi_f = sc[31]*sc[35];
    k_f = 1e-06 * 2.16e+13*exp(-0.23*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[4];
    Kc = exp((g_RT[31] + g_RT[35]) - (g_RT[47] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[47] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 199: NH + NO <=> N2O + H */
    phi_f = sc[31]*sc[35];
    k_f = 1e-06 * 3.65e+14*exp(-0.45*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[37]*sc[1];
    Kc = exp((g_RT[31] + g_RT[35]) - (g_RT[37] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[37] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 200: NH2 + O <=> OH + NH */
    phi_f = sc[32]*sc[2];
    k_f = 1e-06 * 3e+12;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[31];
    Kc = exp((g_RT[32] + g_RT[2]) - (g_RT[4] + g_RT[31]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[32] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[31] += 1 * qdot;

    /*reaction 201: NH2 + O <=> H + HNO */
    phi_f = sc[32]*sc[2];
    k_f = 1e-06 * 3.9e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[38];
    Kc = exp((g_RT[32] + g_RT[2]) - (g_RT[1] + g_RT[38]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[32] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[38] += 1 * qdot;

    /*reaction 202: NH2 + H <=> NH + H2 */
    phi_f = sc[32]*sc[1];
    k_f = 1e-06 * 4e+13*exp(-1836.94/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[0];
    Kc = exp((g_RT[32] + g_RT[1]) - (g_RT[31] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[32] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[31] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 203: NH2 + OH <=> NH + H2O */
    phi_f = sc[32]*sc[4];
    k_f = 1e-06 * 9e+07*exp(1.5*tc[0]+231.505/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[5];
    Kc = exp((g_RT[32] + g_RT[4]) - (g_RT[31] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[32] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[31] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 204: NNH <=> N2 + H */
    phi_f = sc[34];
    k_f = 1 * 3.3e+08;
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[1];
    Kc = refC * exp((g_RT[34]) - (g_RT[47] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[34] -= 1 * qdot;
    wdot[47] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 205: NNH + M <=> N2 + H + M */
    phi_f = sc[34];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * alpha * 1.3e+14*exp(-0.11*tc[0]-2506.29/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[1];
    Kc = refC * exp((g_RT[34]) - (g_RT[47] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[34] -= 1 * qdot;
    wdot[47] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 206: NNH + O2 <=> HO2 + N2 */
    phi_f = sc[34]*sc[3];
    k_f = 1e-06 * 5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[47];
    Kc = exp((g_RT[34] + g_RT[3]) - (g_RT[6] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[34] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[47] += 1 * qdot;

    /*reaction 207: NNH + O <=> OH + N2 */
    phi_f = sc[34]*sc[2];
    k_f = 1e-06 * 2.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[47];
    Kc = exp((g_RT[34] + g_RT[2]) - (g_RT[4] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[34] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[47] += 1 * qdot;

    /*reaction 208: NNH + O <=> NH + NO */
    phi_f = sc[34]*sc[2];
    k_f = 1e-06 * 7e+13;
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[35];
    Kc = exp((g_RT[34] + g_RT[2]) - (g_RT[31] + g_RT[35]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[34] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[31] += 1 * qdot;
    wdot[35] += 1 * qdot;

    /*reaction 209: NNH + H <=> H2 + N2 */
    phi_f = sc[34]*sc[1];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[47];
    Kc = exp((g_RT[34] + g_RT[1]) - (g_RT[0] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[34] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[47] += 1 * qdot;

    /*reaction 210: NNH + OH <=> H2O + N2 */
    phi_f = sc[34]*sc[4];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[47];
    Kc = exp((g_RT[34] + g_RT[4]) - (g_RT[5] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[34] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[47] += 1 * qdot;

    /*reaction 211: NNH + CH3 <=> CH4 + N2 */
    phi_f = sc[34]*sc[12];
    k_f = 1e-06 * 2.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[47];
    Kc = exp((g_RT[34] + g_RT[12]) - (g_RT[13] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[34] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[47] += 1 * qdot;

    /*reaction 212: H + NO + M <=> HNO + M */
    phi_f = sc[1]*sc[35];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-12 * alpha * 4.48e+19*exp(-1.32*tc[0]-372.421/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[35]) - (g_RT[38]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[38] += 1 * qdot;

    /*reaction 213: HNO + O <=> NO + OH */
    phi_f = sc[38]*sc[2];
    k_f = 1e-06 * 2.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[4];
    Kc = exp((g_RT[38] + g_RT[2]) - (g_RT[35] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[38] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 214: HNO + H <=> H2 + NO */
    phi_f = sc[38]*sc[1];
    k_f = 1e-06 * 9e+11*exp(0.72*tc[0]-332.159/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[35];
    Kc = exp((g_RT[38] + g_RT[1]) - (g_RT[0] + g_RT[35]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[38] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[35] += 1 * qdot;

    /*reaction 215: HNO + OH <=> NO + H2O */
    phi_f = sc[38]*sc[4];
    k_f = 1e-06 * 1.3e+07*exp(1.9*tc[0]+478.108/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[5];
    Kc = exp((g_RT[38] + g_RT[4]) - (g_RT[35] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[38] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 216: HNO + O2 <=> HO2 + NO */
    phi_f = sc[38]*sc[3];
    k_f = 1e-06 * 1e+13*exp(-6542.53/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[35];
    Kc = exp((g_RT[38] + g_RT[3]) - (g_RT[6] + g_RT[35]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[38] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[35] += 1 * qdot;

    /*reaction 217: CN + O <=> CO + N */
    phi_f = sc[39]*sc[2];
    k_f = 1e-06 * 7.7e+13;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[30];
    Kc = exp((g_RT[39] + g_RT[2]) - (g_RT[14] + g_RT[30]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[39] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[30] += 1 * qdot;

    /*reaction 218: CN + OH <=> NCO + H */
    phi_f = sc[39]*sc[4];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[46]*sc[1];
    Kc = exp((g_RT[39] + g_RT[4]) - (g_RT[46] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[39] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[46] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 219: CN + H2O <=> HCN + OH */
    phi_f = sc[39]*sc[5];
    k_f = 1e-06 * 8e+12*exp(-3754.4/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[40]*sc[4];
    Kc = exp((g_RT[39] + g_RT[5]) - (g_RT[40] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[39] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[40] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 220: CN + O2 <=> NCO + O */
    phi_f = sc[39]*sc[3];
    k_f = 1e-06 * 6.14e+12*exp(+221.439/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[46]*sc[2];
    Kc = exp((g_RT[39] + g_RT[3]) - (g_RT[46] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[39] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[46] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 221: CN + H2 <=> HCN + H */
    phi_f = sc[39]*sc[0];
    k_f = 1e-06 * 295000*exp(2.45*tc[0]-1127.33/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[40]*sc[1];
    Kc = exp((g_RT[39] + g_RT[0]) - (g_RT[40] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[39] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[40] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 222: NCO + O <=> NO + CO */
    phi_f = sc[46]*sc[2];
    k_f = 1e-06 * 2.35e+13;
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[14];
    Kc = exp((g_RT[46] + g_RT[2]) - (g_RT[35] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[46] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 223: NCO + H <=> NH + CO */
    phi_f = sc[46]*sc[1];
    k_f = 1e-06 * 5.4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[14];
    Kc = exp((g_RT[46] + g_RT[1]) - (g_RT[31] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[46] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[31] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 224: NCO + OH <=> NO + H + CO */
    phi_f = sc[46]*sc[4];
    k_f = 1e-06 * 2.5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[1]*sc[14];
    Kc = refC * exp((g_RT[46] + g_RT[4]) - (g_RT[35] + g_RT[1] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[46] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 225: NCO + N <=> N2 + CO */
    phi_f = sc[46]*sc[30];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[14];
    Kc = exp((g_RT[46] + g_RT[30]) - (g_RT[47] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[46] -= 1 * qdot;
    wdot[30] -= 1 * qdot;
    wdot[47] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 226: NCO + O2 <=> NO + CO2 */
    phi_f = sc[46]*sc[3];
    k_f = 1e-06 * 2e+12*exp(-10065.4/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[15];
    Kc = exp((g_RT[46] + g_RT[3]) - (g_RT[35] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[46] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 227: NCO + M <=> N + CO + M */
    phi_f = sc[46];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * alpha * 3.1e+14*exp(-27201.8/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[14];
    Kc = refC * exp((g_RT[46]) - (g_RT[30] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[46] -= 1 * qdot;
    wdot[30] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 228: NCO + NO <=> N2O + CO */
    phi_f = sc[46]*sc[35];
    k_f = 1e-06 * 1.9e+17*exp(-1.52*tc[0]-372.421/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[37]*sc[14];
    Kc = exp((g_RT[46] + g_RT[35]) - (g_RT[37] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[46] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[37] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 229: NCO + NO <=> N2 + CO2 */
    phi_f = sc[46]*sc[35];
    k_f = 1e-06 * 3.8e+18*exp(-2*tc[0]-402.617/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[15];
    Kc = exp((g_RT[46] + g_RT[35]) - (g_RT[47] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[46] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[47] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 230: HCN + M <=> H + CN + M */
    phi_f = sc[40];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * alpha * 1.04e+29*exp(-3.3*tc[0]-63714.1/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[39];
    Kc = refC * exp((g_RT[40]) - (g_RT[1] + g_RT[39]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[40] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[39] += 1 * qdot;

    /*reaction 231: HCN + O <=> NCO + H */
    phi_f = sc[40]*sc[2];
    k_f = 1e-06 * 20300*exp(2.64*tc[0]-2506.29/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[46]*sc[1];
    Kc = exp((g_RT[40] + g_RT[2]) - (g_RT[46] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[40] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[46] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 232: HCN + O <=> NH + CO */
    phi_f = sc[40]*sc[2];
    k_f = 1e-06 * 5070*exp(2.64*tc[0]-2506.29/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[14];
    Kc = exp((g_RT[40] + g_RT[2]) - (g_RT[31] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[40] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[31] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 233: HCN + O <=> CN + OH */
    phi_f = sc[40]*sc[2];
    k_f = 1e-06 * 3.91e+09*exp(1.58*tc[0]-13387/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[39]*sc[4];
    Kc = exp((g_RT[40] + g_RT[2]) - (g_RT[39] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[40] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[39] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 234: HCN + OH <=> HOCN + H */
    phi_f = sc[40]*sc[4];
    k_f = 1e-06 * 1.1e+06*exp(2.03*tc[0]-6728.74/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[44]*sc[1];
    Kc = exp((g_RT[40] + g_RT[4]) - (g_RT[44] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[40] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[44] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 235: HCN + OH <=> HNCO + H */
    phi_f = sc[40]*sc[4];
    k_f = 1e-06 * 4400*exp(2.26*tc[0]-3220.94/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[45]*sc[1];
    Kc = exp((g_RT[40] + g_RT[4]) - (g_RT[45] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[40] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[45] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 236: HCN + OH <=> NH2 + CO */
    phi_f = sc[40]*sc[4];
    k_f = 1e-06 * 160*exp(2.56*tc[0]-4529.44/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[14];
    Kc = exp((g_RT[40] + g_RT[4]) - (g_RT[32] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[40] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 237: H + HCN (+M) <=> H2CN (+M) */
    phi_f = sc[1]*sc[40];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 3.3e+13;
    redP = 1e-12 * alpha / k_f * 1.4e+26*exp(-3.4*tc[0]-956.215/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[41];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[40]) - (g_RT[41]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[40] -= 1 * qdot;
    wdot[41] += 1 * qdot;

    /*reaction 238: H2CN + N <=> N2 + CH2 */
    phi_f = sc[41]*sc[30];
    k_f = 1e-06 * 6e+13*exp(-201.309/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[10];
    Kc = exp((g_RT[41] + g_RT[30]) - (g_RT[47] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[41] -= 1 * qdot;
    wdot[30] -= 1 * qdot;
    wdot[47] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 239: C + N2 <=> CN + N */
    phi_f = sc[8]*sc[47];
    k_f = 1e-06 * 6.3e+13*exp(-23160.5/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[39]*sc[30];
    Kc = exp((g_RT[8] + g_RT[47]) - (g_RT[39] + g_RT[30]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[47] -= 1 * qdot;
    wdot[39] += 1 * qdot;
    wdot[30] += 1 * qdot;

    /*reaction 240: CH + N2 <=> HCN + N */
    phi_f = sc[9]*sc[47];
    k_f = 1e-06 * 3.12e+09*exp(0.88*tc[0]-10130.9/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[40]*sc[30];
    Kc = exp((g_RT[9] + g_RT[47]) - (g_RT[40] + g_RT[30]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[47] -= 1 * qdot;
    wdot[40] += 1 * qdot;
    wdot[30] += 1 * qdot;

    /*reaction 241: CH + N2 (+M) <=> HCNN (+M) */
    phi_f = sc[9]*sc[47];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + 0*sc[48];
    k_f = 1e-06 * 3.1e+12*exp(0.15*tc[0]);
    redP = 1e-12 * alpha / k_f * 1.3e+25*exp(-3.16*tc[0]-372.421/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.333*exp(T/-235))+ (0.667*exp(T/-2117))+ (exp(-4536/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[42];
    Kc = 1.0 / (refC) * exp((g_RT[9] + g_RT[47]) - (g_RT[42]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[47] -= 1 * qdot;
    wdot[42] += 1 * qdot;

    /*reaction 242: CH2 + N2 <=> HCN + NH */
    phi_f = sc[10]*sc[47];
    k_f = 1e-06 * 1e+13*exp(-37242.1/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[40]*sc[31];
    Kc = exp((g_RT[10] + g_RT[47]) - (g_RT[40] + g_RT[31]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[47] -= 1 * qdot;
    wdot[40] += 1 * qdot;
    wdot[31] += 1 * qdot;

    /*reaction 243: CH2(S) + N2 <=> NH + HCN */
    phi_f = sc[11]*sc[47];
    k_f = 1e-06 * 1e+11*exp(-32712.6/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[40];
    Kc = exp((g_RT[11] + g_RT[47]) - (g_RT[31] + g_RT[40]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[47] -= 1 * qdot;
    wdot[31] += 1 * qdot;
    wdot[40] += 1 * qdot;

    /*reaction 244: C + NO <=> CN + O */
    phi_f = sc[8]*sc[35];
    k_f = 1e-06 * 1.9e+13;
    q_f = phi_f * k_f;
    phi_r = sc[39]*sc[2];
    Kc = exp((g_RT[8] + g_RT[35]) - (g_RT[39] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[39] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 245: C + NO <=> CO + N */
    phi_f = sc[8]*sc[35];
    k_f = 1e-06 * 2.9e+13;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[30];
    Kc = exp((g_RT[8] + g_RT[35]) - (g_RT[14] + g_RT[30]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[30] += 1 * qdot;

    /*reaction 246: CH + NO <=> HCN + O */
    phi_f = sc[9]*sc[35];
    k_f = 1e-06 * 4.1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[40]*sc[2];
    Kc = exp((g_RT[9] + g_RT[35]) - (g_RT[40] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[40] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 247: CH + NO <=> H + NCO */
    phi_f = sc[9]*sc[35];
    k_f = 1e-06 * 1.62e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[46];
    Kc = exp((g_RT[9] + g_RT[35]) - (g_RT[1] + g_RT[46]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[46] += 1 * qdot;

    /*reaction 248: CH + NO <=> N + HCO */
    phi_f = sc[9]*sc[35];
    k_f = 1e-06 * 2.46e+13;
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[16];
    Kc = exp((g_RT[9] + g_RT[35]) - (g_RT[30] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[30] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 249: CH2 + NO <=> H + HNCO */
    phi_f = sc[10]*sc[35];
    k_f = 1e-06 * 3.1e+17*exp(-1.38*tc[0]-639.155/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[45];
    Kc = exp((g_RT[10] + g_RT[35]) - (g_RT[1] + g_RT[45]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[45] += 1 * qdot;

    /*reaction 250: CH2 + NO <=> OH + HCN */
    phi_f = sc[10]*sc[35];
    k_f = 1e-06 * 2.9e+14*exp(-0.69*tc[0]-382.486/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[40];
    Kc = exp((g_RT[10] + g_RT[35]) - (g_RT[4] + g_RT[40]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[40] += 1 * qdot;

    /*reaction 251: CH2 + NO <=> H + HCNO */
    phi_f = sc[10]*sc[35];
    k_f = 1e-06 * 3.8e+13*exp(-0.36*tc[0]-291.897/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[43];
    Kc = exp((g_RT[10] + g_RT[35]) - (g_RT[1] + g_RT[43]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[43] += 1 * qdot;

    /*reaction 252: CH2(S) + NO <=> H + HNCO */
    phi_f = sc[11]*sc[35];
    k_f = 1e-06 * 3.1e+17*exp(-1.38*tc[0]-639.155/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[45];
    Kc = exp((g_RT[11] + g_RT[35]) - (g_RT[1] + g_RT[45]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[45] += 1 * qdot;

    /*reaction 253: CH2(S) + NO <=> OH + HCN */
    phi_f = sc[11]*sc[35];
    k_f = 1e-06 * 2.9e+14*exp(-0.69*tc[0]-382.486/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[40];
    Kc = exp((g_RT[11] + g_RT[35]) - (g_RT[4] + g_RT[40]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[40] += 1 * qdot;

    /*reaction 254: CH2(S) + NO <=> H + HCNO */
    phi_f = sc[11]*sc[35];
    k_f = 1e-06 * 3.8e+13*exp(-0.36*tc[0]-291.897/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[43];
    Kc = exp((g_RT[11] + g_RT[35]) - (g_RT[1] + g_RT[43]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[43] += 1 * qdot;

    /*reaction 255: CH3 + NO <=> HCN + H2O */
    phi_f = sc[12]*sc[35];
    k_f = 1e-06 * 9.6e+13*exp(-14494.2/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[40]*sc[5];
    Kc = exp((g_RT[12] + g_RT[35]) - (g_RT[40] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[40] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 256: CH3 + NO <=> H2CN + OH */
    phi_f = sc[12]*sc[35];
    k_f = 1e-06 * 1e+12*exp(-10946.1/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[41]*sc[4];
    Kc = exp((g_RT[12] + g_RT[35]) - (g_RT[41] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[41] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 257: HCNN + O <=> CO + H + N2 */
    phi_f = sc[42]*sc[2];
    k_f = 1e-06 * 2.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[1]*sc[47];
    Kc = refC * exp((g_RT[42] + g_RT[2]) - (g_RT[14] + g_RT[1] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[42] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[47] += 1 * qdot;

    /*reaction 258: HCNN + O <=> HCN + NO */
    phi_f = sc[42]*sc[2];
    k_f = 1e-06 * 2e+12;
    q_f = phi_f * k_f;
    phi_r = sc[40]*sc[35];
    Kc = exp((g_RT[42] + g_RT[2]) - (g_RT[40] + g_RT[35]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[42] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[40] += 1 * qdot;
    wdot[35] += 1 * qdot;

    /*reaction 259: HCNN + O2 <=> O + HCO + N2 */
    phi_f = sc[42]*sc[3];
    k_f = 1e-06 * 1.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[16]*sc[47];
    Kc = refC * exp((g_RT[42] + g_RT[3]) - (g_RT[2] + g_RT[16] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[42] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[47] += 1 * qdot;

    /*reaction 260: HCNN + OH <=> H + HCO + N2 */
    phi_f = sc[42]*sc[4];
    k_f = 1e-06 * 1.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16]*sc[47];
    Kc = refC * exp((g_RT[42] + g_RT[4]) - (g_RT[1] + g_RT[16] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[42] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[47] += 1 * qdot;

    /*reaction 261: HCNN + H <=> CH2 + N2 */
    phi_f = sc[42]*sc[1];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[47];
    Kc = exp((g_RT[42] + g_RT[1]) - (g_RT[10] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[42] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[47] += 1 * qdot;

    /*reaction 262: HNCO + O <=> NH + CO2 */
    phi_f = sc[45]*sc[2];
    k_f = 1e-06 * 9.8e+07*exp(1.41*tc[0]-4277.81/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[15];
    Kc = exp((g_RT[45] + g_RT[2]) - (g_RT[31] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[45] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[31] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 263: HNCO + O <=> HNO + CO */
    phi_f = sc[45]*sc[2];
    k_f = 1e-06 * 1.5e+08*exp(1.57*tc[0]-22143.9/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[14];
    Kc = exp((g_RT[45] + g_RT[2]) - (g_RT[38] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[45] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[38] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 264: HNCO + O <=> NCO + OH */
    phi_f = sc[45]*sc[2];
    k_f = 1e-06 * 2.2e+06*exp(2.11*tc[0]-5737.29/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[46]*sc[4];
    Kc = exp((g_RT[45] + g_RT[2]) - (g_RT[46] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[45] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[46] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 265: HNCO + H <=> NH2 + CO */
    phi_f = sc[45]*sc[1];
    k_f = 1e-06 * 2.25e+07*exp(1.7*tc[0]-1912.43/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[14];
    Kc = exp((g_RT[45] + g_RT[1]) - (g_RT[32] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[45] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 266: HNCO + H <=> H2 + NCO */
    phi_f = sc[45]*sc[1];
    k_f = 1e-06 * 105000*exp(2.5*tc[0]-6693.51/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[46];
    Kc = exp((g_RT[45] + g_RT[1]) - (g_RT[0] + g_RT[46]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[45] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[46] += 1 * qdot;

    /*reaction 267: HNCO + OH <=> NCO + H2O */
    phi_f = sc[45]*sc[4];
    k_f = 1e-06 * 3.3e+07*exp(1.5*tc[0]-1811.78/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[46]*sc[5];
    Kc = exp((g_RT[45] + g_RT[4]) - (g_RT[46] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[45] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[46] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 268: HNCO + OH <=> NH2 + CO2 */
    phi_f = sc[45]*sc[4];
    k_f = 1e-06 * 3.3e+06*exp(1.5*tc[0]-1811.78/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[15];
    Kc = exp((g_RT[45] + g_RT[4]) - (g_RT[32] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[45] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 269: HNCO + M <=> NH + CO + M */
    phi_f = sc[45];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * alpha * 1.18e+16*exp(-42637.1/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[14];
    Kc = refC * exp((g_RT[45]) - (g_RT[31] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[45] -= 1 * qdot;
    wdot[31] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 270: HCNO + H <=> H + HNCO */
    phi_f = sc[43]*sc[1];
    k_f = 1e-06 * 2.1e+15*exp(-0.69*tc[0]-1434.32/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[45];
    Kc = exp((g_RT[43] + g_RT[1]) - (g_RT[1] + g_RT[45]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[43] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[45] += 1 * qdot;

    /*reaction 271: HCNO + H <=> OH + HCN */
    phi_f = sc[43]*sc[1];
    k_f = 1e-06 * 2.7e+11*exp(0.18*tc[0]-1066.94/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[40];
    Kc = exp((g_RT[43] + g_RT[1]) - (g_RT[4] + g_RT[40]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[43] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[40] += 1 * qdot;

    /*reaction 272: HCNO + H <=> NH2 + CO */
    phi_f = sc[43]*sc[1];
    k_f = 1e-06 * 1.7e+14*exp(-0.75*tc[0]-1454.45/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[14];
    Kc = exp((g_RT[43] + g_RT[1]) - (g_RT[32] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[43] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 273: HOCN + H <=> H + HNCO */
    phi_f = sc[44]*sc[1];
    k_f = 1e-06 * 2e+07*exp(2*tc[0]-1006.54/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[45];
    Kc = exp((g_RT[44] + g_RT[1]) - (g_RT[1] + g_RT[45]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[44] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[45] += 1 * qdot;

    /*reaction 274: HCCO + NO <=> HCNO + CO */
    phi_f = sc[27]*sc[35];
    k_f = 1e-06 * 9e+12;
    q_f = phi_f * k_f;
    phi_r = sc[43]*sc[14];
    Kc = exp((g_RT[27] + g_RT[35]) - (g_RT[43] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[27] -= 1 * qdot;
    wdot[35] -= 1 * qdot;
    wdot[43] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 275: CH3 + N <=> H2CN + H */
    phi_f = sc[12]*sc[30];
    k_f = 1e-06 * 6.1e+14*exp(-0.31*tc[0]-145.949/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[41]*sc[1];
    Kc = exp((g_RT[12] + g_RT[30]) - (g_RT[41] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[30] -= 1 * qdot;
    wdot[41] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 276: CH3 + N <=> HCN + H2 */
    phi_f = sc[12]*sc[30];
    k_f = 1e-06 * 3.7e+12*exp(0.15*tc[0]+45.2944/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[40]*sc[0];
    Kc = exp((g_RT[12] + g_RT[30]) - (g_RT[40] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[30] -= 1 * qdot;
    wdot[40] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 277: NH3 + H <=> NH2 + H2 */
    phi_f = sc[33]*sc[1];
    k_f = 1e-06 * 540000*exp(2.4*tc[0]-4989.93/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[0];
    Kc = exp((g_RT[33] + g_RT[1]) - (g_RT[32] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 278: NH3 + OH <=> NH2 + H2O */
    phi_f = sc[33]*sc[4];
    k_f = 1e-06 * 5e+07*exp(1.6*tc[0]-480.624/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[5];
    Kc = exp((g_RT[33] + g_RT[4]) - (g_RT[32] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 279: NH3 + O <=> NH2 + OH */
    phi_f = sc[33]*sc[2];
    k_f = 1e-06 * 9.4e+06*exp(1.94*tc[0]-3251.13/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[4];
    Kc = exp((g_RT[33] + g_RT[2]) - (g_RT[32] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 280: NH + CO2 <=> HNO + CO */
    phi_f = sc[31]*sc[15];
    k_f = 1e-06 * 1e+13*exp(-7221.94/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[14];
    Kc = exp((g_RT[31] + g_RT[15]) - (g_RT[38] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[38] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 281: CN + NO2 <=> NCO + NO */
    phi_f = sc[39]*sc[36];
    k_f = 1e-06 * 6.16e+15*exp(-0.752*tc[0]-173.629/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[46]*sc[35];
    Kc = exp((g_RT[39] + g_RT[36]) - (g_RT[46] + g_RT[35]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[39] -= 1 * qdot;
    wdot[36] -= 1 * qdot;
    wdot[46] += 1 * qdot;
    wdot[35] += 1 * qdot;

    /*reaction 282: NCO + NO2 <=> N2O + CO2 */
    phi_f = sc[46]*sc[36];
    k_f = 1e-06 * 3.25e+12*exp(+354.806/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[37]*sc[15];
    Kc = exp((g_RT[46] + g_RT[36]) - (g_RT[37] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[46] -= 1 * qdot;
    wdot[36] -= 1 * qdot;
    wdot[37] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 283: N + CO2 <=> NO + CO */
    phi_f = sc[30]*sc[15];
    k_f = 1e-06 * 3e+12*exp(-5686.97/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[14];
    Kc = exp((g_RT[30] + g_RT[15]) - (g_RT[35] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[30] -= 1 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 284: O + CH3 => H + H2 + CO */
    phi_f = sc[2]*sc[12];
    k_f = 1e-06 * 3.37e+13;
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 285: O + C2H4 <=> H + CH2CHO */
    phi_f = sc[2]*sc[24];
    k_f = 1e-06 * 6.7e+06*exp(1.83*tc[0]-110.72/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[51];
    Kc = exp((g_RT[2] + g_RT[24]) - (g_RT[1] + g_RT[51]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[24] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[51] += 1 * qdot;

    /*reaction 286: O + C2H5 <=> H + CH3CHO */
    phi_f = sc[2]*sc[25];
    k_f = 1e-06 * 1.096e+14;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[52];
    Kc = exp((g_RT[2] + g_RT[25]) - (g_RT[1] + g_RT[52]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[25] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[52] += 1 * qdot;

    /*reaction 287: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = 1e-06 * 5e+15*exp(-8721.69/tc[1]);
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

    /*reaction 288: OH + CH3 => H2 + CH2O */
    phi_f = sc[4]*sc[12];
    k_f = 1e-06 * 8e+09*exp(0.5*tc[0]+883.241/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 289: CH + H2 (+M) <=> CH3 (+M) */
    phi_f = sc[9]*sc[0];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 1.97e+12*exp(0.43*tc[0]+186.21/tc[1]);
    redP = 1e-12 * alpha / k_f * 4.82e+25*exp(-2.8*tc[0]-296.93/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.422*exp(T/-122))+ (0.578*exp(T/-2535))+ (exp(-9365/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[12];
    Kc = 1.0 / (refC) * exp((g_RT[9] + g_RT[0]) - (g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 290: CH2 + O2 => 2 H + CO2 */
    phi_f = sc[10]*sc[3];
    k_f = 1e-06 * 5.8e+12*exp(-754.907/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[1] += 2 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 291: CH2 + O2 <=> O + CH2O */
    phi_f = sc[10]*sc[3];
    k_f = 1e-06 * 2.4e+12*exp(-754.907/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[17];
    Kc = exp((g_RT[10] + g_RT[3]) - (g_RT[2] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 292: CH2 + CH2 => 2 H + C2H2 */
    phi_f = sc[10]*sc[10];
    k_f = 1e-06 * 2e+14*exp(-5530.45/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[1] += 2 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 293: CH2(S) + H2O => H2 + CH2O */
    phi_f = sc[11]*sc[5];
    k_f = 1e-06 * 6.82e+10*exp(0.25*tc[0]+470.559/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 294: C2H3 + O2 <=> O + CH2CHO */
    phi_f = sc[23]*sc[3];
    k_f = 1e-06 * 3.03e+11*exp(0.29*tc[0]-5.53598/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[51];
    Kc = exp((g_RT[23] + g_RT[3]) - (g_RT[2] + g_RT[51]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[23] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[51] += 1 * qdot;

    /*reaction 295: C2H3 + O2 <=> HO2 + C2H2 */
    phi_f = sc[23]*sc[3];
    k_f = 1e-06 * 1.337e+06*exp(1.61*tc[0]+193.256/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[22];
    Kc = exp((g_RT[23] + g_RT[3]) - (g_RT[6] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[23] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 296: O + CH3CHO <=> OH + CH2CHO */
    phi_f = sc[2]*sc[52];
    k_f = 1e-06 * 5.84e+12*exp(-909.914/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[51];
    Kc = exp((g_RT[2] + g_RT[52]) - (g_RT[4] + g_RT[51]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[52] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[51] += 1 * qdot;

    /*reaction 297: O + CH3CHO => OH + CH3 + CO */
    phi_f = sc[2]*sc[52];
    k_f = 1e-06 * 5.84e+12*exp(-909.914/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[52] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 298: O2 + CH3CHO => HO2 + CH3 + CO */
    phi_f = sc[3]*sc[52];
    k_f = 1e-06 * 3.01e+13*exp(-19703.1/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[52] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 299: H + CH3CHO <=> CH2CHO + H2 */
    phi_f = sc[1]*sc[52];
    k_f = 1e-06 * 2.05e+09*exp(1.16*tc[0]-1210.37/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[51]*sc[0];
    Kc = exp((g_RT[1] + g_RT[52]) - (g_RT[51] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[52] -= 1 * qdot;
    wdot[51] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 300: H + CH3CHO => CH3 + H2 + CO */
    phi_f = sc[1]*sc[52];
    k_f = 1e-06 * 2.05e+09*exp(1.16*tc[0]-1210.37/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[52] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 301: OH + CH3CHO => CH3 + H2O + CO */
    phi_f = sc[4]*sc[52];
    k_f = 1e-06 * 2.343e+10*exp(0.73*tc[0]+560.141/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[52] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 302: HO2 + CH3CHO => CH3 + H2O2 + CO */
    phi_f = sc[6]*sc[52];
    k_f = 1e-06 * 3.01e+12*exp(-6000.5/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[52] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 303: CH3 + CH3CHO => CH3 + CH4 + CO */
    phi_f = sc[12]*sc[52];
    k_f = 1e-06 * 2.72e+06*exp(1.77*tc[0]-2979.37/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[52] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 304: H + CH2CO (+M) <=> CH2CHO (+M) */
    phi_f = sc[1]*sc[28];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 4.865e+11*exp(0.422*tc[0]+883.241/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.012e+42*exp(-7.63*tc[0]-1939.61/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.535*exp(T/-201))+ (0.465*exp(T/-1773))+ (exp(-5333/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[51];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[28]) - (g_RT[51]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[28] -= 1 * qdot;
    wdot[51] += 1 * qdot;

    /*reaction 305: O + CH2CHO => H + CH2 + CO2 */
    phi_f = sc[2]*sc[51];
    k_f = 1e-06 * 1.5e+14;
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[51] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 306: O2 + CH2CHO => OH + CO + CH2O */
    phi_f = sc[3]*sc[51];
    k_f = 1e-06 * 1.81e+10;
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[51] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 307: O2 + CH2CHO => OH + 2 HCO */
    phi_f = sc[3]*sc[51];
    k_f = 1e-06 * 2.35e+10;
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[51] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[16] += 2 * qdot;

    /*reaction 308: H + CH2CHO <=> CH3 + HCO */
    phi_f = sc[1]*sc[51];
    k_f = 1e-06 * 2.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[16];
    Kc = exp((g_RT[1] + g_RT[51]) - (g_RT[12] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[51] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 309: H + CH2CHO <=> CH2CO + H2 */
    phi_f = sc[1]*sc[51];
    k_f = 1e-06 * 1.1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[28]*sc[0];
    Kc = exp((g_RT[1] + g_RT[51]) - (g_RT[28] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[51] -= 1 * qdot;
    wdot[28] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 310: OH + CH2CHO <=> H2O + CH2CO */
    phi_f = sc[4]*sc[51];
    k_f = 1e-06 * 1.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[28];
    Kc = exp((g_RT[4] + g_RT[51]) - (g_RT[5] + g_RT[28]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[51] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[28] += 1 * qdot;

    /*reaction 311: OH + CH2CHO <=> HCO + CH2OH */
    phi_f = sc[4]*sc[51];
    k_f = 1e-06 * 3.01e+13;
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[18];
    Kc = exp((g_RT[4] + g_RT[51]) - (g_RT[16] + g_RT[18]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[51] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[18] += 1 * qdot;

    /*reaction 312: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    phi_f = sc[12]*sc[25];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 9.43e+12;
    redP = 1e-12 * alpha / k_f * 2.71e+74*exp(-16.82*tc[0]-6575.24/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.8473*exp(T/-291))+ (0.1527*exp(T/-2742))+ (exp(-7748/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[50];
    Kc = 1.0 / (refC) * exp((g_RT[12] + g_RT[25]) - (g_RT[50]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[25] -= 1 * qdot;
    wdot[50] += 1 * qdot;

    /*reaction 313: O + C3H8 <=> OH + C3H7 */
    phi_f = sc[2]*sc[50];
    k_f = 1e-06 * 193000*exp(2.68*tc[0]-1870.16/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[49];
    Kc = exp((g_RT[2] + g_RT[50]) - (g_RT[4] + g_RT[49]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[50] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[49] += 1 * qdot;

    /*reaction 314: H + C3H8 <=> C3H7 + H2 */
    phi_f = sc[1]*sc[50];
    k_f = 1e-06 * 1.32e+06*exp(2.54*tc[0]-3400.1/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[49]*sc[0];
    Kc = exp((g_RT[1] + g_RT[50]) - (g_RT[49] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[50] -= 1 * qdot;
    wdot[49] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 315: OH + C3H8 <=> C3H7 + H2O */
    phi_f = sc[4]*sc[50];
    k_f = 1e-06 * 3.16e+07*exp(1.8*tc[0]-470.055/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[49]*sc[5];
    Kc = exp((g_RT[4] + g_RT[50]) - (g_RT[49] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[50] -= 1 * qdot;
    wdot[49] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 316: C3H7 + H2O2 <=> HO2 + C3H8 */
    phi_f = sc[49]*sc[7];
    k_f = 1e-06 * 378*exp(2.72*tc[0]-754.907/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[50];
    Kc = exp((g_RT[49] + g_RT[7]) - (g_RT[6] + g_RT[50]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[49] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[50] += 1 * qdot;

    /*reaction 317: CH3 + C3H8 <=> C3H7 + CH4 */
    phi_f = sc[12]*sc[50];
    k_f = 1e-06 * 0.903*exp(3.65*tc[0]-3600.4/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[49]*sc[13];
    Kc = exp((g_RT[12] + g_RT[50]) - (g_RT[49] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[50] -= 1 * qdot;
    wdot[49] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 318: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    phi_f = sc[12]*sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 2.55e+06*exp(1.6*tc[0]-2868.65/tc[1]);
    redP = 1e-12 * alpha / k_f * 3e+63*exp(-14.6*tc[0]-9144.44/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.8106*exp(T/-277))+ (0.1894*exp(T/-8748))+ (exp(-7891/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[49];
    Kc = 1.0 / (refC) * exp((g_RT[12] + g_RT[24]) - (g_RT[49]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[24] -= 1 * qdot;
    wdot[49] += 1 * qdot;

    /*reaction 319: O + C3H7 <=> C2H5 + CH2O */
    phi_f = sc[2]*sc[49];
    k_f = 1e-06 * 9.64e+13;
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[17];
    Kc = exp((g_RT[2] + g_RT[49]) - (g_RT[25] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[49] -= 1 * qdot;
    wdot[25] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 320: H + C3H7 (+M) <=> C3H8 (+M) */
    phi_f = sc[1]*sc[49];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 3.613e+13;
    redP = 1e-12 * alpha / k_f * 4.42e+61*exp(-13.545*tc[0]-5715.65/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.685*exp(T/-369))+ (0.315*exp(T/-3285))+ (exp(-6667/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[50];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[49]) - (g_RT[50]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[49] -= 1 * qdot;
    wdot[50] += 1 * qdot;

    /*reaction 321: H + C3H7 <=> CH3 + C2H5 */
    phi_f = sc[1]*sc[49];
    k_f = 1e-06 * 4.06e+06*exp(2.19*tc[0]-447.911/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[25];
    Kc = exp((g_RT[1] + g_RT[49]) - (g_RT[12] + g_RT[25]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[49] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[25] += 1 * qdot;

    /*reaction 322: OH + C3H7 <=> C2H5 + CH2OH */
    phi_f = sc[4]*sc[49];
    k_f = 1e-06 * 2.41e+13;
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[18];
    Kc = exp((g_RT[4] + g_RT[49]) - (g_RT[25] + g_RT[18]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[49] -= 1 * qdot;
    wdot[25] += 1 * qdot;
    wdot[18] += 1 * qdot;

    /*reaction 323: HO2 + C3H7 <=> O2 + C3H8 */
    phi_f = sc[6]*sc[49];
    k_f = 1e-06 * 2.55e+10*exp(0.255*tc[0]+474.585/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[50];
    Kc = exp((g_RT[6] + g_RT[49]) - (g_RT[3] + g_RT[50]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[49] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[50] += 1 * qdot;

    /*reaction 324: HO2 + C3H7 => OH + C2H5 + CH2O */
    phi_f = sc[6]*sc[49];
    k_f = 1e-06 * 2.41e+13;
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[49] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[25] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 325: CH3 + C3H7 <=> 2 C2H5 */
    phi_f = sc[12]*sc[49];
    k_f = 1e-06 * 1.927e+13*exp(-0.32*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[25];
    Kc = exp((g_RT[12] + g_RT[49]) - (2 * g_RT[25]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[49] -= 1 * qdot;
    wdot[25] += 2 * qdot;

    return;
}


/*compute the progress rate for each reaction */
void progressRate(double * qdot, double * sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[53];                /*Gibbs free energy */
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
    double refC = 101325 / 8.31451 / T;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 53; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*reaction 1: 2 O + M <=> O2 + M */
    phi_f = sc[2]*sc[2];
    alpha = mixture + 1.4*sc[0] + 14.4*sc[5] + sc[13] + 0.75*sc[14] + 2.6*sc[15] + 2*sc[26] + -0.17*sc[48];
    k_f = 1e-12 * alpha * 1.2e+17*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[3];
    Kc = 1.0 / (refC) * exp((2 * g_RT[2]) - (g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[0] = q_f - q_r;

    /*reaction 2: O + H + M <=> OH + M */
    phi_f = sc[2]*sc[1];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-12 * alpha * 5e+17*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[4];
    Kc = 1.0 / (refC) * exp((g_RT[2] + g_RT[1]) - (g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[1] = q_f - q_r;

    /*reaction 3: O + H2 <=> H + OH */
    phi_f = sc[2]*sc[0];
    k_f = 1e-06 * 38700*exp(2.7*tc[0]-3150.48/tc[1]);
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
    k_f = 1e-06 * 5.06e+13;
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

    /*reaction 12: O + CO (+M) <=> CO2 (+M) */
    phi_f = sc[2]*sc[14];
    alpha = mixture + sc[0] + 5*sc[3] + 5*sc[5] + sc[13] + 0.5*sc[14] + 2.5*sc[15] + 2*sc[26] + -0.5*sc[48];
    k_f = 1e-06 * 1.8e+10*exp(-1200.3/tc[1]);
    redP = 1e-12 * alpha / k_f * 6.02e+14*exp(-1509.81/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
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
    k_f = 1e-06 * 1.35e+07*exp(2*tc[0]-956.215/tc[1]);
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
    k_f = 1e-06 * 6.94e+06*exp(2*tc[0]-956.215/tc[1]);
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
    k_f = 1e-06 * 1.25e+07*exp(1.83*tc[0]-110.72/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[16];
    Kc = exp((g_RT[2] + g_RT[24]) - (g_RT[12] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[24] = q_f - q_r;

    /*reaction 26: O + C2H5 <=> CH3 + CH2O */
    phi_f = sc[2]*sc[25];
    k_f = 1e-06 * 2.24e+13;
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
    alpha = mixture + -1*sc[3] + -1*sc[5] + -0.25*sc[14] + 0.5*sc[15] + 0.5*sc[26] + -1*sc[47] + -1*sc[48];
    k_f = 1e-12 * alpha * 2.8e+18*exp(-0.86*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[32] = q_f - q_r;

    /*reaction 34: H + 2 O2 <=> HO2 + O2 */
    phi_f = sc[1]*sc[3]*sc[3];
    k_f = 1e-12 * 2.08e+19*exp(-1.24*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = 1.0 / (refC) * exp((g_RT[1] + 2 * g_RT[3]) - (g_RT[6] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[33] = q_f - q_r;

    /*reaction 35: H + O2 + H2O <=> HO2 + H2O */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = 1e-12 * 1.126e+19*exp(-0.76*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[34] = q_f - q_r;

    /*reaction 36: H + O2 + N2 <=> HO2 + N2 */
    phi_f = sc[1]*sc[3]*sc[47];
    k_f = 1e-12 * 2.6e+19*exp(-1.24*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[47];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[47]) - (g_RT[6] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[35] = q_f - q_r;

    /*reaction 37: H + O2 + AR <=> HO2 + AR */
    phi_f = sc[1]*sc[3]*sc[48];
    k_f = 1e-12 * 7e+17*exp(-0.8*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[48];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[48]) - (g_RT[6] + g_RT[48]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[36] = q_f - q_r;

    /*reaction 38: H + O2 <=> O + OH */
    phi_f = sc[1]*sc[3];
    k_f = 1e-06 * 2.65e+16*exp(-0.6707*tc[0]-8576.25/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[37] = q_f - q_r;

    /*reaction 39: 2 H + M <=> H2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[0] + -1*sc[5] + sc[13] + -1*sc[15] + 2*sc[26] + -0.37*sc[48];
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
    alpha = mixture + -0.27*sc[0] + 2.65*sc[5] + sc[13] + 2*sc[26] + -0.62*sc[48];
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
    k_f = 1e-06 * 4.48e+13*exp(-537.494/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[44] = q_f - q_r;

    /*reaction 46: H + HO2 <=> 2 OH */
    phi_f = sc[1]*sc[6];
    k_f = 1e-06 * 8.4e+13*exp(-319.577/tc[1]);
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
    k_f = 1e-06 * 1.65e+14;
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = exp((g_RT[1] + g_RT[9]) - (g_RT[8] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[48] = q_f - q_r;

    /*reaction 50: H + CH2 (+M) <=> CH3 (+M) */
    phi_f = sc[1]*sc[10];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 6e+14;
    redP = 1e-12 * alpha / k_f * 1.04e+26*exp(-2.76*tc[0]-805.234/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.438*exp(T/-91))+ (0.562*exp(T/-5836))+ (exp(-8552/T)));
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
    alpha = mixture + sc[0] + 5*sc[5] + 2*sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 1.39e+16*exp(-0.534*tc[0]-269.753/tc[1]);
    redP = 1e-12 * alpha / k_f * 2.62e+33*exp(-4.76*tc[0]-1227.98/tc[1]);
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 1.09e+12*exp(0.48*tc[0]+130.851/tc[1]);
    redP = 1e-12 * alpha / k_f * 2.47e+24*exp(-2.57*tc[0]-213.89/tc[1]);
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
    k_f = 1e-06 * 5.74e+07*exp(1.9*tc[0]-1379.97/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[0];
    Kc = exp((g_RT[1] + g_RT[17]) - (g_RT[16] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[57] = q_f - q_r;

    /*reaction 59: H + CH2OH (+M) <=> CH3OH (+M) */
    phi_f = sc[1]*sc[18];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = 1e-06 * 1.055e+12*exp(0.5*tc[0]-43.2813/tc[1]);
    redP = 1e-12 * alpha / k_f * 4.36e+31*exp(-4.65*tc[0]-2556.62/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.4*exp(T/-100))+ (0.6*exp(T/-90000))+ (exp(-10000/T)));
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
    k_f = 1e-06 * 1.65e+11*exp(0.65*tc[0]+142.929/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = exp((g_RT[1] + g_RT[18]) - (g_RT[4] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[60] = q_f - q_r;

    /*reaction 62: H + CH2OH <=> CH2(S) + H2O */
    phi_f = sc[1]*sc[18];
    k_f = 1e-06 * 3.28e+13*exp(-0.09*tc[0]-306.995/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = exp((g_RT[1] + g_RT[18]) - (g_RT[11] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[61] = q_f - q_r;

    /*reaction 63: H + CH3O (+M) <=> CH3OH (+M) */
    phi_f = sc[1]*sc[19];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = 1e-06 * 2.43e+12*exp(0.515*tc[0]-25.1636/tc[1]);
    redP = 1e-12 * alpha / k_f * 4.66e+41*exp(-7.44*tc[0]-7086.06/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.3*exp(T/-100))+ (0.7*exp(T/-90000))+ (exp(-10000/T)));
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
    k_f = 1e-06 * 4.15e+07*exp(1.63*tc[0]-968.294/tc[1]);
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
    k_f = 1e-06 * 1.5e+12*exp(0.5*tc[0]+55.3598/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = exp((g_RT[1] + g_RT[19]) - (g_RT[4] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[65] = q_f - q_r;

    /*reaction 67: H + CH3O <=> CH2(S) + H2O */
    phi_f = sc[1]*sc[19];
    k_f = 1e-06 * 2.62e+14*exp(-0.23*tc[0]-538.5/tc[1]);
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 5.4e+11*exp(0.454*tc[0]-915.954/tc[1]);
    redP = 1e-12 * alpha / k_f * 6e+41*exp(-7.62*tc[0]-3507.8/tc[1]);
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
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
    k_f = 1e-06 * 1.45e+13*exp(+251.636/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[86] = q_f - q_r;

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    phi_f = sc[4]*sc[7];
    k_f = 1e-06 * 2e+12*exp(-214.897/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[87] = q_f - q_r;

    /*reaction 89: OH + H2O2 <=> HO2 + H2O */
    phi_f = sc[4]*sc[7];
    k_f = 1e-06 * 1.7e+18*exp(-14801.2/tc[1]);
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
    k_f = 1e-06 * 2.79e+18*exp(-1.43*tc[0]-669.351/tc[1]);
    redP = 1e-12 * alpha / k_f * 4e+36*exp(-5.92*tc[0]-1580.27/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.588*exp(T/-195))+ (0.412*exp(T/-5900))+ (exp(-6394/T)));
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
    k_f = 1e-06 * 6.44e+17*exp(-1.34*tc[0]-713.135/tc[1]);
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
    k_f = 1e-06 * 3.78e+13;
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
    k_f = 1e-06 * 5.6e+06*exp(2*tc[0]-6039.26/tc[1]);
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
    k_f = 1e-06 * 6.71e+13;
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[16];
    Kc = exp((g_RT[9] + g_RT[3]) - (g_RT[2] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[124] = q_f - q_r;

    /*reaction 126: CH + H2 <=> H + CH2 */
    phi_f = sc[9]*sc[0];
    k_f = 1e-06 * 1.08e+14*exp(-1565.17/tc[1]);
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
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
    k_f = 1e-06 * 1.9e+14*exp(-7947.66/tc[1]);
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

    /*reaction 135: CH2 + O2 => OH + H + CO */
    phi_f = sc[10]*sc[3];
    k_f = 1e-06 * 5e+12*exp(-754.907/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
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
    k_f = 1e-06 * 1.6e+15*exp(-6011.07/tc[1]);
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
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
    phi_f = sc[11]*sc[47];
    k_f = 1e-06 * 1.5e+13*exp(-301.963/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[47];
    Kc = exp((g_RT[11] + g_RT[47]) - (g_RT[10] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[141] = q_f - q_r;

    /*reaction 143: CH2(S) + AR <=> CH2 + AR */
    phi_f = sc[11]*sc[48];
    k_f = 1e-06 * 9e+12*exp(-301.963/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[48];
    Kc = exp((g_RT[11] + g_RT[48]) - (g_RT[10] + g_RT[48]));
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
    k_f = 1e-06 * 4.82e+17*exp(-1.16*tc[0]-576.246/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.88e+38*exp(-6.36*tc[0]-2536.49/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.3973*exp(T/-208))+ (0.6027*exp(T/-3922))+ (exp(-10180/T)));
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
    k_f = 1e-06 * 3.56e+13*exp(-15339.7/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[19];
    Kc = exp((g_RT[12] + g_RT[3]) - (g_RT[2] + g_RT[19]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[154] = q_f - q_r;

    /*reaction 156: CH3 + O2 <=> OH + CH2O */
    phi_f = sc[12]*sc[3];
    k_f = 1e-06 * 2.31e+12*exp(-10224/tc[1]);
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 6.77e+16*exp(-1.18*tc[0]-329.139/tc[1]);
    redP = 1e-12 * alpha / k_f * 3.4e+41*exp(-7.03*tc[0]-1390.04/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.381*exp(T/-73.2))+ (0.619*exp(T/-1180))+ (exp(-9999/T)));
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
    k_f = 1e-06 * 6.84e+12*exp(0.1*tc[0]-5334.68/tc[1]);
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
    k_f = 1e-06 * 1.5e+18*exp(-1*tc[0]-8555.61/tc[1]);
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
    k_f = 1e-06 * 1.345e+13*exp(-201.309/tc[1]);
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
    k_f = 1e-06 * 1e+13*exp(+379.97/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[14];
    Kc = exp((g_RT[21] + g_RT[3]) - (g_RT[16] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[170] = q_f - q_r;

    /*reaction 172: C2H + H2 <=> H + C2H2 */
    phi_f = sc[21]*sc[0];
    k_f = 1e-06 * 5.68e+10*exp(0.9*tc[0]-1003.02/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = exp((g_RT[21] + g_RT[0]) - (g_RT[1] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[171] = q_f - q_r;

    /*reaction 173: C2H3 + O2 <=> HCO + CH2O */
    phi_f = sc[23]*sc[3];
    k_f = 1e-06 * 4.58e+16*exp(-1.39*tc[0]-510.82/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[17];
    Kc = exp((g_RT[23] + g_RT[3]) - (g_RT[16] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[172] = q_f - q_r;

    /*reaction 174: C2H4 (+M) <=> H2 + C2H2 (+M) */
    phi_f = sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1 * 8e+12*exp(0.44*tc[0]-43668.8/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.58e+51*exp(-9.3*tc[0]-49219.9/tc[1]);
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
    k_f = 1e-06 * 3.2e+12*exp(-429.794/tc[1]);
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

    /*reaction 178: N + NO <=> N2 + O */
    phi_f = sc[30]*sc[35];
    k_f = 1e-06 * 2.7e+13*exp(-178.661/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[2];
    Kc = exp((g_RT[30] + g_RT[35]) - (g_RT[47] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[177] = q_f - q_r;

    /*reaction 179: N + O2 <=> NO + O */
    phi_f = sc[30]*sc[3];
    k_f = 1e-06 * 9e+09*exp(1*tc[0]-3271.26/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[2];
    Kc = exp((g_RT[30] + g_RT[3]) - (g_RT[35] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[178] = q_f - q_r;

    /*reaction 180: N + OH <=> NO + H */
    phi_f = sc[30]*sc[4];
    k_f = 1e-06 * 3.36e+13*exp(-193.759/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[1];
    Kc = exp((g_RT[30] + g_RT[4]) - (g_RT[35] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[179] = q_f - q_r;

    /*reaction 181: N2O + O <=> N2 + O2 */
    phi_f = sc[37]*sc[2];
    k_f = 1e-06 * 1.4e+12*exp(-5440.36/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[3];
    Kc = exp((g_RT[37] + g_RT[2]) - (g_RT[47] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[180] = q_f - q_r;

    /*reaction 182: N2O + O <=> 2 NO */
    phi_f = sc[37]*sc[2];
    k_f = 1e-06 * 2.9e+13*exp(-11650.7/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[35];
    Kc = exp((g_RT[37] + g_RT[2]) - (2 * g_RT[35]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[181] = q_f - q_r;

    /*reaction 183: N2O + H <=> N2 + OH */
    phi_f = sc[37]*sc[1];
    k_f = 1e-06 * 3.87e+14*exp(-9501.76/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[4];
    Kc = exp((g_RT[37] + g_RT[1]) - (g_RT[47] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[182] = q_f - q_r;

    /*reaction 184: N2O + OH <=> N2 + HO2 */
    phi_f = sc[37]*sc[4];
    k_f = 1e-06 * 2e+12*exp(-10598.9/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[6];
    Kc = exp((g_RT[37] + g_RT[4]) - (g_RT[47] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[183] = q_f - q_r;

    /*reaction 185: N2O (+M) <=> N2 + O (+M) */
    phi_f = sc[37];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.375*sc[48];
    k_f = 1 * 7.91e+10*exp(-28193.3/tc[1]);
    redP = 1e-12 * alpha / k_f * 6.37e+14*exp(-28505.3/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[2];
    Kc = refC * exp((g_RT[37]) - (g_RT[47] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[184] = q_f - q_r;

    /*reaction 186: HO2 + NO <=> NO2 + OH */
    phi_f = sc[6]*sc[35];
    k_f = 1e-06 * 2.11e+12*exp(+241.57/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[36]*sc[4];
    Kc = exp((g_RT[6] + g_RT[35]) - (g_RT[36] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[185] = q_f - q_r;

    /*reaction 187: NO + O + M <=> NO2 + M */
    phi_f = sc[35]*sc[2];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-12 * alpha * 1.06e+20*exp(-1.41*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[36];
    Kc = 1.0 / (refC) * exp((g_RT[35] + g_RT[2]) - (g_RT[36]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[186] = q_f - q_r;

    /*reaction 188: NO2 + O <=> NO + O2 */
    phi_f = sc[36]*sc[2];
    k_f = 1e-06 * 3.9e+12*exp(+120.785/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[3];
    Kc = exp((g_RT[36] + g_RT[2]) - (g_RT[35] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[187] = q_f - q_r;

    /*reaction 189: NO2 + H <=> NO + OH */
    phi_f = sc[36]*sc[1];
    k_f = 1e-06 * 1.32e+14*exp(-181.178/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[4];
    Kc = exp((g_RT[36] + g_RT[1]) - (g_RT[35] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[188] = q_f - q_r;

    /*reaction 190: NH + O <=> NO + H */
    phi_f = sc[31]*sc[2];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[1];
    Kc = exp((g_RT[31] + g_RT[2]) - (g_RT[35] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[189] = q_f - q_r;

    /*reaction 191: NH + H <=> N + H2 */
    phi_f = sc[31]*sc[1];
    k_f = 1e-06 * 3.2e+13*exp(-166.08/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[0];
    Kc = exp((g_RT[31] + g_RT[1]) - (g_RT[30] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[190] = q_f - q_r;

    /*reaction 192: NH + OH <=> HNO + H */
    phi_f = sc[31]*sc[4];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[1];
    Kc = exp((g_RT[31] + g_RT[4]) - (g_RT[38] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[191] = q_f - q_r;

    /*reaction 193: NH + OH <=> N + H2O */
    phi_f = sc[31]*sc[4];
    k_f = 1e-06 * 2e+09*exp(1.2*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[5];
    Kc = exp((g_RT[31] + g_RT[4]) - (g_RT[30] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[192] = q_f - q_r;

    /*reaction 194: NH + O2 <=> HNO + O */
    phi_f = sc[31]*sc[3];
    k_f = 1e-06 * 461000*exp(2*tc[0]-3271.26/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[2];
    Kc = exp((g_RT[31] + g_RT[3]) - (g_RT[38] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[193] = q_f - q_r;

    /*reaction 195: NH + O2 <=> NO + OH */
    phi_f = sc[31]*sc[3];
    k_f = 1e-06 * 1.28e+06*exp(1.5*tc[0]-50.3271/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[4];
    Kc = exp((g_RT[31] + g_RT[3]) - (g_RT[35] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[194] = q_f - q_r;

    /*reaction 196: NH + N <=> N2 + H */
    phi_f = sc[31]*sc[30];
    k_f = 1e-06 * 1.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[1];
    Kc = exp((g_RT[31] + g_RT[30]) - (g_RT[47] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[195] = q_f - q_r;

    /*reaction 197: NH + H2O <=> HNO + H2 */
    phi_f = sc[31]*sc[5];
    k_f = 1e-06 * 2e+13*exp(-6970.31/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[0];
    Kc = exp((g_RT[31] + g_RT[5]) - (g_RT[38] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[196] = q_f - q_r;

    /*reaction 198: NH + NO <=> N2 + OH */
    phi_f = sc[31]*sc[35];
    k_f = 1e-06 * 2.16e+13*exp(-0.23*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[4];
    Kc = exp((g_RT[31] + g_RT[35]) - (g_RT[47] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[197] = q_f - q_r;

    /*reaction 199: NH + NO <=> N2O + H */
    phi_f = sc[31]*sc[35];
    k_f = 1e-06 * 3.65e+14*exp(-0.45*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[37]*sc[1];
    Kc = exp((g_RT[31] + g_RT[35]) - (g_RT[37] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[198] = q_f - q_r;

    /*reaction 200: NH2 + O <=> OH + NH */
    phi_f = sc[32]*sc[2];
    k_f = 1e-06 * 3e+12;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[31];
    Kc = exp((g_RT[32] + g_RT[2]) - (g_RT[4] + g_RT[31]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[199] = q_f - q_r;

    /*reaction 201: NH2 + O <=> H + HNO */
    phi_f = sc[32]*sc[2];
    k_f = 1e-06 * 3.9e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[38];
    Kc = exp((g_RT[32] + g_RT[2]) - (g_RT[1] + g_RT[38]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[200] = q_f - q_r;

    /*reaction 202: NH2 + H <=> NH + H2 */
    phi_f = sc[32]*sc[1];
    k_f = 1e-06 * 4e+13*exp(-1836.94/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[0];
    Kc = exp((g_RT[32] + g_RT[1]) - (g_RT[31] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[201] = q_f - q_r;

    /*reaction 203: NH2 + OH <=> NH + H2O */
    phi_f = sc[32]*sc[4];
    k_f = 1e-06 * 9e+07*exp(1.5*tc[0]+231.505/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[5];
    Kc = exp((g_RT[32] + g_RT[4]) - (g_RT[31] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[202] = q_f - q_r;

    /*reaction 204: NNH <=> N2 + H */
    phi_f = sc[34];
    k_f = 1 * 3.3e+08;
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[1];
    Kc = refC * exp((g_RT[34]) - (g_RT[47] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[203] = q_f - q_r;

    /*reaction 205: NNH + M <=> N2 + H + M */
    phi_f = sc[34];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * alpha * 1.3e+14*exp(-0.11*tc[0]-2506.29/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[1];
    Kc = refC * exp((g_RT[34]) - (g_RT[47] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[204] = q_f - q_r;

    /*reaction 206: NNH + O2 <=> HO2 + N2 */
    phi_f = sc[34]*sc[3];
    k_f = 1e-06 * 5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[47];
    Kc = exp((g_RT[34] + g_RT[3]) - (g_RT[6] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[205] = q_f - q_r;

    /*reaction 207: NNH + O <=> OH + N2 */
    phi_f = sc[34]*sc[2];
    k_f = 1e-06 * 2.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[47];
    Kc = exp((g_RT[34] + g_RT[2]) - (g_RT[4] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[206] = q_f - q_r;

    /*reaction 208: NNH + O <=> NH + NO */
    phi_f = sc[34]*sc[2];
    k_f = 1e-06 * 7e+13;
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[35];
    Kc = exp((g_RT[34] + g_RT[2]) - (g_RT[31] + g_RT[35]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[207] = q_f - q_r;

    /*reaction 209: NNH + H <=> H2 + N2 */
    phi_f = sc[34]*sc[1];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[47];
    Kc = exp((g_RT[34] + g_RT[1]) - (g_RT[0] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[208] = q_f - q_r;

    /*reaction 210: NNH + OH <=> H2O + N2 */
    phi_f = sc[34]*sc[4];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[47];
    Kc = exp((g_RT[34] + g_RT[4]) - (g_RT[5] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[209] = q_f - q_r;

    /*reaction 211: NNH + CH3 <=> CH4 + N2 */
    phi_f = sc[34]*sc[12];
    k_f = 1e-06 * 2.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[47];
    Kc = exp((g_RT[34] + g_RT[12]) - (g_RT[13] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[210] = q_f - q_r;

    /*reaction 212: H + NO + M <=> HNO + M */
    phi_f = sc[1]*sc[35];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-12 * alpha * 4.48e+19*exp(-1.32*tc[0]-372.421/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[35]) - (g_RT[38]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[211] = q_f - q_r;

    /*reaction 213: HNO + O <=> NO + OH */
    phi_f = sc[38]*sc[2];
    k_f = 1e-06 * 2.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[4];
    Kc = exp((g_RT[38] + g_RT[2]) - (g_RT[35] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[212] = q_f - q_r;

    /*reaction 214: HNO + H <=> H2 + NO */
    phi_f = sc[38]*sc[1];
    k_f = 1e-06 * 9e+11*exp(0.72*tc[0]-332.159/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[35];
    Kc = exp((g_RT[38] + g_RT[1]) - (g_RT[0] + g_RT[35]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[213] = q_f - q_r;

    /*reaction 215: HNO + OH <=> NO + H2O */
    phi_f = sc[38]*sc[4];
    k_f = 1e-06 * 1.3e+07*exp(1.9*tc[0]+478.108/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[5];
    Kc = exp((g_RT[38] + g_RT[4]) - (g_RT[35] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[214] = q_f - q_r;

    /*reaction 216: HNO + O2 <=> HO2 + NO */
    phi_f = sc[38]*sc[3];
    k_f = 1e-06 * 1e+13*exp(-6542.53/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[35];
    Kc = exp((g_RT[38] + g_RT[3]) - (g_RT[6] + g_RT[35]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[215] = q_f - q_r;

    /*reaction 217: CN + O <=> CO + N */
    phi_f = sc[39]*sc[2];
    k_f = 1e-06 * 7.7e+13;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[30];
    Kc = exp((g_RT[39] + g_RT[2]) - (g_RT[14] + g_RT[30]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[216] = q_f - q_r;

    /*reaction 218: CN + OH <=> NCO + H */
    phi_f = sc[39]*sc[4];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[46]*sc[1];
    Kc = exp((g_RT[39] + g_RT[4]) - (g_RT[46] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[217] = q_f - q_r;

    /*reaction 219: CN + H2O <=> HCN + OH */
    phi_f = sc[39]*sc[5];
    k_f = 1e-06 * 8e+12*exp(-3754.4/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[40]*sc[4];
    Kc = exp((g_RT[39] + g_RT[5]) - (g_RT[40] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[218] = q_f - q_r;

    /*reaction 220: CN + O2 <=> NCO + O */
    phi_f = sc[39]*sc[3];
    k_f = 1e-06 * 6.14e+12*exp(+221.439/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[46]*sc[2];
    Kc = exp((g_RT[39] + g_RT[3]) - (g_RT[46] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[219] = q_f - q_r;

    /*reaction 221: CN + H2 <=> HCN + H */
    phi_f = sc[39]*sc[0];
    k_f = 1e-06 * 295000*exp(2.45*tc[0]-1127.33/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[40]*sc[1];
    Kc = exp((g_RT[39] + g_RT[0]) - (g_RT[40] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[220] = q_f - q_r;

    /*reaction 222: NCO + O <=> NO + CO */
    phi_f = sc[46]*sc[2];
    k_f = 1e-06 * 2.35e+13;
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[14];
    Kc = exp((g_RT[46] + g_RT[2]) - (g_RT[35] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[221] = q_f - q_r;

    /*reaction 223: NCO + H <=> NH + CO */
    phi_f = sc[46]*sc[1];
    k_f = 1e-06 * 5.4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[14];
    Kc = exp((g_RT[46] + g_RT[1]) - (g_RT[31] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[222] = q_f - q_r;

    /*reaction 224: NCO + OH <=> NO + H + CO */
    phi_f = sc[46]*sc[4];
    k_f = 1e-06 * 2.5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[1]*sc[14];
    Kc = refC * exp((g_RT[46] + g_RT[4]) - (g_RT[35] + g_RT[1] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[223] = q_f - q_r;

    /*reaction 225: NCO + N <=> N2 + CO */
    phi_f = sc[46]*sc[30];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[14];
    Kc = exp((g_RT[46] + g_RT[30]) - (g_RT[47] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[224] = q_f - q_r;

    /*reaction 226: NCO + O2 <=> NO + CO2 */
    phi_f = sc[46]*sc[3];
    k_f = 1e-06 * 2e+12*exp(-10065.4/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[15];
    Kc = exp((g_RT[46] + g_RT[3]) - (g_RT[35] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[225] = q_f - q_r;

    /*reaction 227: NCO + M <=> N + CO + M */
    phi_f = sc[46];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * alpha * 3.1e+14*exp(-27201.8/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[14];
    Kc = refC * exp((g_RT[46]) - (g_RT[30] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[226] = q_f - q_r;

    /*reaction 228: NCO + NO <=> N2O + CO */
    phi_f = sc[46]*sc[35];
    k_f = 1e-06 * 1.9e+17*exp(-1.52*tc[0]-372.421/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[37]*sc[14];
    Kc = exp((g_RT[46] + g_RT[35]) - (g_RT[37] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[227] = q_f - q_r;

    /*reaction 229: NCO + NO <=> N2 + CO2 */
    phi_f = sc[46]*sc[35];
    k_f = 1e-06 * 3.8e+18*exp(-2*tc[0]-402.617/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[15];
    Kc = exp((g_RT[46] + g_RT[35]) - (g_RT[47] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[228] = q_f - q_r;

    /*reaction 230: HCN + M <=> H + CN + M */
    phi_f = sc[40];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * alpha * 1.04e+29*exp(-3.3*tc[0]-63714.1/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[39];
    Kc = refC * exp((g_RT[40]) - (g_RT[1] + g_RT[39]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[229] = q_f - q_r;

    /*reaction 231: HCN + O <=> NCO + H */
    phi_f = sc[40]*sc[2];
    k_f = 1e-06 * 20300*exp(2.64*tc[0]-2506.29/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[46]*sc[1];
    Kc = exp((g_RT[40] + g_RT[2]) - (g_RT[46] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[230] = q_f - q_r;

    /*reaction 232: HCN + O <=> NH + CO */
    phi_f = sc[40]*sc[2];
    k_f = 1e-06 * 5070*exp(2.64*tc[0]-2506.29/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[14];
    Kc = exp((g_RT[40] + g_RT[2]) - (g_RT[31] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[231] = q_f - q_r;

    /*reaction 233: HCN + O <=> CN + OH */
    phi_f = sc[40]*sc[2];
    k_f = 1e-06 * 3.91e+09*exp(1.58*tc[0]-13387/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[39]*sc[4];
    Kc = exp((g_RT[40] + g_RT[2]) - (g_RT[39] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[232] = q_f - q_r;

    /*reaction 234: HCN + OH <=> HOCN + H */
    phi_f = sc[40]*sc[4];
    k_f = 1e-06 * 1.1e+06*exp(2.03*tc[0]-6728.74/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[44]*sc[1];
    Kc = exp((g_RT[40] + g_RT[4]) - (g_RT[44] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[233] = q_f - q_r;

    /*reaction 235: HCN + OH <=> HNCO + H */
    phi_f = sc[40]*sc[4];
    k_f = 1e-06 * 4400*exp(2.26*tc[0]-3220.94/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[45]*sc[1];
    Kc = exp((g_RT[40] + g_RT[4]) - (g_RT[45] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[234] = q_f - q_r;

    /*reaction 236: HCN + OH <=> NH2 + CO */
    phi_f = sc[40]*sc[4];
    k_f = 1e-06 * 160*exp(2.56*tc[0]-4529.44/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[14];
    Kc = exp((g_RT[40] + g_RT[4]) - (g_RT[32] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[235] = q_f - q_r;

    /*reaction 237: H + HCN (+M) <=> H2CN (+M) */
    phi_f = sc[1]*sc[40];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 3.3e+13;
    redP = 1e-12 * alpha / k_f * 1.4e+26*exp(-3.4*tc[0]-956.215/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[41];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[40]) - (g_RT[41]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[236] = q_f - q_r;

    /*reaction 238: H2CN + N <=> N2 + CH2 */
    phi_f = sc[41]*sc[30];
    k_f = 1e-06 * 6e+13*exp(-201.309/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[47]*sc[10];
    Kc = exp((g_RT[41] + g_RT[30]) - (g_RT[47] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[237] = q_f - q_r;

    /*reaction 239: C + N2 <=> CN + N */
    phi_f = sc[8]*sc[47];
    k_f = 1e-06 * 6.3e+13*exp(-23160.5/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[39]*sc[30];
    Kc = exp((g_RT[8] + g_RT[47]) - (g_RT[39] + g_RT[30]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[238] = q_f - q_r;

    /*reaction 240: CH + N2 <=> HCN + N */
    phi_f = sc[9]*sc[47];
    k_f = 1e-06 * 3.12e+09*exp(0.88*tc[0]-10130.9/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[40]*sc[30];
    Kc = exp((g_RT[9] + g_RT[47]) - (g_RT[40] + g_RT[30]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[239] = q_f - q_r;

    /*reaction 241: CH + N2 (+M) <=> HCNN (+M) */
    phi_f = sc[9]*sc[47];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + 0*sc[48];
    k_f = 1e-06 * 3.1e+12*exp(0.15*tc[0]);
    redP = 1e-12 * alpha / k_f * 1.3e+25*exp(-3.16*tc[0]-372.421/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.333*exp(T/-235))+ (0.667*exp(T/-2117))+ (exp(-4536/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[42];
    Kc = 1.0 / (refC) * exp((g_RT[9] + g_RT[47]) - (g_RT[42]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[240] = q_f - q_r;

    /*reaction 242: CH2 + N2 <=> HCN + NH */
    phi_f = sc[10]*sc[47];
    k_f = 1e-06 * 1e+13*exp(-37242.1/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[40]*sc[31];
    Kc = exp((g_RT[10] + g_RT[47]) - (g_RT[40] + g_RT[31]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[241] = q_f - q_r;

    /*reaction 243: CH2(S) + N2 <=> NH + HCN */
    phi_f = sc[11]*sc[47];
    k_f = 1e-06 * 1e+11*exp(-32712.6/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[40];
    Kc = exp((g_RT[11] + g_RT[47]) - (g_RT[31] + g_RT[40]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[242] = q_f - q_r;

    /*reaction 244: C + NO <=> CN + O */
    phi_f = sc[8]*sc[35];
    k_f = 1e-06 * 1.9e+13;
    q_f = phi_f * k_f;
    phi_r = sc[39]*sc[2];
    Kc = exp((g_RT[8] + g_RT[35]) - (g_RT[39] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[243] = q_f - q_r;

    /*reaction 245: C + NO <=> CO + N */
    phi_f = sc[8]*sc[35];
    k_f = 1e-06 * 2.9e+13;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[30];
    Kc = exp((g_RT[8] + g_RT[35]) - (g_RT[14] + g_RT[30]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[244] = q_f - q_r;

    /*reaction 246: CH + NO <=> HCN + O */
    phi_f = sc[9]*sc[35];
    k_f = 1e-06 * 4.1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[40]*sc[2];
    Kc = exp((g_RT[9] + g_RT[35]) - (g_RT[40] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[245] = q_f - q_r;

    /*reaction 247: CH + NO <=> H + NCO */
    phi_f = sc[9]*sc[35];
    k_f = 1e-06 * 1.62e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[46];
    Kc = exp((g_RT[9] + g_RT[35]) - (g_RT[1] + g_RT[46]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[246] = q_f - q_r;

    /*reaction 248: CH + NO <=> N + HCO */
    phi_f = sc[9]*sc[35];
    k_f = 1e-06 * 2.46e+13;
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[16];
    Kc = exp((g_RT[9] + g_RT[35]) - (g_RT[30] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[247] = q_f - q_r;

    /*reaction 249: CH2 + NO <=> H + HNCO */
    phi_f = sc[10]*sc[35];
    k_f = 1e-06 * 3.1e+17*exp(-1.38*tc[0]-639.155/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[45];
    Kc = exp((g_RT[10] + g_RT[35]) - (g_RT[1] + g_RT[45]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[248] = q_f - q_r;

    /*reaction 250: CH2 + NO <=> OH + HCN */
    phi_f = sc[10]*sc[35];
    k_f = 1e-06 * 2.9e+14*exp(-0.69*tc[0]-382.486/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[40];
    Kc = exp((g_RT[10] + g_RT[35]) - (g_RT[4] + g_RT[40]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[249] = q_f - q_r;

    /*reaction 251: CH2 + NO <=> H + HCNO */
    phi_f = sc[10]*sc[35];
    k_f = 1e-06 * 3.8e+13*exp(-0.36*tc[0]-291.897/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[43];
    Kc = exp((g_RT[10] + g_RT[35]) - (g_RT[1] + g_RT[43]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[250] = q_f - q_r;

    /*reaction 252: CH2(S) + NO <=> H + HNCO */
    phi_f = sc[11]*sc[35];
    k_f = 1e-06 * 3.1e+17*exp(-1.38*tc[0]-639.155/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[45];
    Kc = exp((g_RT[11] + g_RT[35]) - (g_RT[1] + g_RT[45]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[251] = q_f - q_r;

    /*reaction 253: CH2(S) + NO <=> OH + HCN */
    phi_f = sc[11]*sc[35];
    k_f = 1e-06 * 2.9e+14*exp(-0.69*tc[0]-382.486/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[40];
    Kc = exp((g_RT[11] + g_RT[35]) - (g_RT[4] + g_RT[40]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[252] = q_f - q_r;

    /*reaction 254: CH2(S) + NO <=> H + HCNO */
    phi_f = sc[11]*sc[35];
    k_f = 1e-06 * 3.8e+13*exp(-0.36*tc[0]-291.897/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[43];
    Kc = exp((g_RT[11] + g_RT[35]) - (g_RT[1] + g_RT[43]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[253] = q_f - q_r;

    /*reaction 255: CH3 + NO <=> HCN + H2O */
    phi_f = sc[12]*sc[35];
    k_f = 1e-06 * 9.6e+13*exp(-14494.2/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[40]*sc[5];
    Kc = exp((g_RT[12] + g_RT[35]) - (g_RT[40] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[254] = q_f - q_r;

    /*reaction 256: CH3 + NO <=> H2CN + OH */
    phi_f = sc[12]*sc[35];
    k_f = 1e-06 * 1e+12*exp(-10946.1/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[41]*sc[4];
    Kc = exp((g_RT[12] + g_RT[35]) - (g_RT[41] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[255] = q_f - q_r;

    /*reaction 257: HCNN + O <=> CO + H + N2 */
    phi_f = sc[42]*sc[2];
    k_f = 1e-06 * 2.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[1]*sc[47];
    Kc = refC * exp((g_RT[42] + g_RT[2]) - (g_RT[14] + g_RT[1] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[256] = q_f - q_r;

    /*reaction 258: HCNN + O <=> HCN + NO */
    phi_f = sc[42]*sc[2];
    k_f = 1e-06 * 2e+12;
    q_f = phi_f * k_f;
    phi_r = sc[40]*sc[35];
    Kc = exp((g_RT[42] + g_RT[2]) - (g_RT[40] + g_RT[35]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[257] = q_f - q_r;

    /*reaction 259: HCNN + O2 <=> O + HCO + N2 */
    phi_f = sc[42]*sc[3];
    k_f = 1e-06 * 1.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[16]*sc[47];
    Kc = refC * exp((g_RT[42] + g_RT[3]) - (g_RT[2] + g_RT[16] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[258] = q_f - q_r;

    /*reaction 260: HCNN + OH <=> H + HCO + N2 */
    phi_f = sc[42]*sc[4];
    k_f = 1e-06 * 1.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16]*sc[47];
    Kc = refC * exp((g_RT[42] + g_RT[4]) - (g_RT[1] + g_RT[16] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[259] = q_f - q_r;

    /*reaction 261: HCNN + H <=> CH2 + N2 */
    phi_f = sc[42]*sc[1];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[47];
    Kc = exp((g_RT[42] + g_RT[1]) - (g_RT[10] + g_RT[47]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[260] = q_f - q_r;

    /*reaction 262: HNCO + O <=> NH + CO2 */
    phi_f = sc[45]*sc[2];
    k_f = 1e-06 * 9.8e+07*exp(1.41*tc[0]-4277.81/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[15];
    Kc = exp((g_RT[45] + g_RT[2]) - (g_RT[31] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[261] = q_f - q_r;

    /*reaction 263: HNCO + O <=> HNO + CO */
    phi_f = sc[45]*sc[2];
    k_f = 1e-06 * 1.5e+08*exp(1.57*tc[0]-22143.9/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[14];
    Kc = exp((g_RT[45] + g_RT[2]) - (g_RT[38] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[262] = q_f - q_r;

    /*reaction 264: HNCO + O <=> NCO + OH */
    phi_f = sc[45]*sc[2];
    k_f = 1e-06 * 2.2e+06*exp(2.11*tc[0]-5737.29/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[46]*sc[4];
    Kc = exp((g_RT[45] + g_RT[2]) - (g_RT[46] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[263] = q_f - q_r;

    /*reaction 265: HNCO + H <=> NH2 + CO */
    phi_f = sc[45]*sc[1];
    k_f = 1e-06 * 2.25e+07*exp(1.7*tc[0]-1912.43/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[14];
    Kc = exp((g_RT[45] + g_RT[1]) - (g_RT[32] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[264] = q_f - q_r;

    /*reaction 266: HNCO + H <=> H2 + NCO */
    phi_f = sc[45]*sc[1];
    k_f = 1e-06 * 105000*exp(2.5*tc[0]-6693.51/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[46];
    Kc = exp((g_RT[45] + g_RT[1]) - (g_RT[0] + g_RT[46]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[265] = q_f - q_r;

    /*reaction 267: HNCO + OH <=> NCO + H2O */
    phi_f = sc[45]*sc[4];
    k_f = 1e-06 * 3.3e+07*exp(1.5*tc[0]-1811.78/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[46]*sc[5];
    Kc = exp((g_RT[45] + g_RT[4]) - (g_RT[46] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[266] = q_f - q_r;

    /*reaction 268: HNCO + OH <=> NH2 + CO2 */
    phi_f = sc[45]*sc[4];
    k_f = 1e-06 * 3.3e+06*exp(1.5*tc[0]-1811.78/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[15];
    Kc = exp((g_RT[45] + g_RT[4]) - (g_RT[32] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[267] = q_f - q_r;

    /*reaction 269: HNCO + M <=> NH + CO + M */
    phi_f = sc[45];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * alpha * 1.18e+16*exp(-42637.1/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[14];
    Kc = refC * exp((g_RT[45]) - (g_RT[31] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[268] = q_f - q_r;

    /*reaction 270: HCNO + H <=> H + HNCO */
    phi_f = sc[43]*sc[1];
    k_f = 1e-06 * 2.1e+15*exp(-0.69*tc[0]-1434.32/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[45];
    Kc = exp((g_RT[43] + g_RT[1]) - (g_RT[1] + g_RT[45]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[269] = q_f - q_r;

    /*reaction 271: HCNO + H <=> OH + HCN */
    phi_f = sc[43]*sc[1];
    k_f = 1e-06 * 2.7e+11*exp(0.18*tc[0]-1066.94/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[40];
    Kc = exp((g_RT[43] + g_RT[1]) - (g_RT[4] + g_RT[40]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[270] = q_f - q_r;

    /*reaction 272: HCNO + H <=> NH2 + CO */
    phi_f = sc[43]*sc[1];
    k_f = 1e-06 * 1.7e+14*exp(-0.75*tc[0]-1454.45/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[14];
    Kc = exp((g_RT[43] + g_RT[1]) - (g_RT[32] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[271] = q_f - q_r;

    /*reaction 273: HOCN + H <=> H + HNCO */
    phi_f = sc[44]*sc[1];
    k_f = 1e-06 * 2e+07*exp(2*tc[0]-1006.54/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[45];
    Kc = exp((g_RT[44] + g_RT[1]) - (g_RT[1] + g_RT[45]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[272] = q_f - q_r;

    /*reaction 274: HCCO + NO <=> HCNO + CO */
    phi_f = sc[27]*sc[35];
    k_f = 1e-06 * 9e+12;
    q_f = phi_f * k_f;
    phi_r = sc[43]*sc[14];
    Kc = exp((g_RT[27] + g_RT[35]) - (g_RT[43] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[273] = q_f - q_r;

    /*reaction 275: CH3 + N <=> H2CN + H */
    phi_f = sc[12]*sc[30];
    k_f = 1e-06 * 6.1e+14*exp(-0.31*tc[0]-145.949/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[41]*sc[1];
    Kc = exp((g_RT[12] + g_RT[30]) - (g_RT[41] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[274] = q_f - q_r;

    /*reaction 276: CH3 + N <=> HCN + H2 */
    phi_f = sc[12]*sc[30];
    k_f = 1e-06 * 3.7e+12*exp(0.15*tc[0]+45.2944/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[40]*sc[0];
    Kc = exp((g_RT[12] + g_RT[30]) - (g_RT[40] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[275] = q_f - q_r;

    /*reaction 277: NH3 + H <=> NH2 + H2 */
    phi_f = sc[33]*sc[1];
    k_f = 1e-06 * 540000*exp(2.4*tc[0]-4989.93/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[0];
    Kc = exp((g_RT[33] + g_RT[1]) - (g_RT[32] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[276] = q_f - q_r;

    /*reaction 278: NH3 + OH <=> NH2 + H2O */
    phi_f = sc[33]*sc[4];
    k_f = 1e-06 * 5e+07*exp(1.6*tc[0]-480.624/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[5];
    Kc = exp((g_RT[33] + g_RT[4]) - (g_RT[32] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[277] = q_f - q_r;

    /*reaction 279: NH3 + O <=> NH2 + OH */
    phi_f = sc[33]*sc[2];
    k_f = 1e-06 * 9.4e+06*exp(1.94*tc[0]-3251.13/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[4];
    Kc = exp((g_RT[33] + g_RT[2]) - (g_RT[32] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[278] = q_f - q_r;

    /*reaction 280: NH + CO2 <=> HNO + CO */
    phi_f = sc[31]*sc[15];
    k_f = 1e-06 * 1e+13*exp(-7221.94/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[38]*sc[14];
    Kc = exp((g_RT[31] + g_RT[15]) - (g_RT[38] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[279] = q_f - q_r;

    /*reaction 281: CN + NO2 <=> NCO + NO */
    phi_f = sc[39]*sc[36];
    k_f = 1e-06 * 6.16e+15*exp(-0.752*tc[0]-173.629/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[46]*sc[35];
    Kc = exp((g_RT[39] + g_RT[36]) - (g_RT[46] + g_RT[35]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[280] = q_f - q_r;

    /*reaction 282: NCO + NO2 <=> N2O + CO2 */
    phi_f = sc[46]*sc[36];
    k_f = 1e-06 * 3.25e+12*exp(+354.806/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[37]*sc[15];
    Kc = exp((g_RT[46] + g_RT[36]) - (g_RT[37] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[281] = q_f - q_r;

    /*reaction 283: N + CO2 <=> NO + CO */
    phi_f = sc[30]*sc[15];
    k_f = 1e-06 * 3e+12*exp(-5686.97/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[14];
    Kc = exp((g_RT[30] + g_RT[15]) - (g_RT[35] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[282] = q_f - q_r;

    /*reaction 284: O + CH3 => H + H2 + CO */
    phi_f = sc[2]*sc[12];
    k_f = 1e-06 * 3.37e+13;
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[283] = q_f - q_r;

    /*reaction 285: O + C2H4 <=> H + CH2CHO */
    phi_f = sc[2]*sc[24];
    k_f = 1e-06 * 6.7e+06*exp(1.83*tc[0]-110.72/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[51];
    Kc = exp((g_RT[2] + g_RT[24]) - (g_RT[1] + g_RT[51]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[284] = q_f - q_r;

    /*reaction 286: O + C2H5 <=> H + CH3CHO */
    phi_f = sc[2]*sc[25];
    k_f = 1e-06 * 1.096e+14;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[52];
    Kc = exp((g_RT[2] + g_RT[25]) - (g_RT[1] + g_RT[52]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[285] = q_f - q_r;

    /*reaction 287: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = 1e-06 * 5e+15*exp(-8721.69/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[286] = q_f - q_r;

    /*reaction 288: OH + CH3 => H2 + CH2O */
    phi_f = sc[4]*sc[12];
    k_f = 1e-06 * 8e+09*exp(0.5*tc[0]+883.241/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[287] = q_f - q_r;

    /*reaction 289: CH + H2 (+M) <=> CH3 (+M) */
    phi_f = sc[9]*sc[0];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 1.97e+12*exp(0.43*tc[0]+186.21/tc[1]);
    redP = 1e-12 * alpha / k_f * 4.82e+25*exp(-2.8*tc[0]-296.93/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.422*exp(T/-122))+ (0.578*exp(T/-2535))+ (exp(-9365/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[12];
    Kc = 1.0 / (refC) * exp((g_RT[9] + g_RT[0]) - (g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[288] = q_f - q_r;

    /*reaction 290: CH2 + O2 => 2 H + CO2 */
    phi_f = sc[10]*sc[3];
    k_f = 1e-06 * 5.8e+12*exp(-754.907/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[289] = q_f - q_r;

    /*reaction 291: CH2 + O2 <=> O + CH2O */
    phi_f = sc[10]*sc[3];
    k_f = 1e-06 * 2.4e+12*exp(-754.907/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[17];
    Kc = exp((g_RT[10] + g_RT[3]) - (g_RT[2] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[290] = q_f - q_r;

    /*reaction 292: CH2 + CH2 => 2 H + C2H2 */
    phi_f = sc[10]*sc[10];
    k_f = 1e-06 * 2e+14*exp(-5530.45/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[291] = q_f - q_r;

    /*reaction 293: CH2(S) + H2O => H2 + CH2O */
    phi_f = sc[11]*sc[5];
    k_f = 1e-06 * 6.82e+10*exp(0.25*tc[0]+470.559/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[292] = q_f - q_r;

    /*reaction 294: C2H3 + O2 <=> O + CH2CHO */
    phi_f = sc[23]*sc[3];
    k_f = 1e-06 * 3.03e+11*exp(0.29*tc[0]-5.53598/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[51];
    Kc = exp((g_RT[23] + g_RT[3]) - (g_RT[2] + g_RT[51]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[293] = q_f - q_r;

    /*reaction 295: C2H3 + O2 <=> HO2 + C2H2 */
    phi_f = sc[23]*sc[3];
    k_f = 1e-06 * 1.337e+06*exp(1.61*tc[0]+193.256/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[22];
    Kc = exp((g_RT[23] + g_RT[3]) - (g_RT[6] + g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[294] = q_f - q_r;

    /*reaction 296: O + CH3CHO <=> OH + CH2CHO */
    phi_f = sc[2]*sc[52];
    k_f = 1e-06 * 5.84e+12*exp(-909.914/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[51];
    Kc = exp((g_RT[2] + g_RT[52]) - (g_RT[4] + g_RT[51]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[295] = q_f - q_r;

    /*reaction 297: O + CH3CHO => OH + CH3 + CO */
    phi_f = sc[2]*sc[52];
    k_f = 1e-06 * 5.84e+12*exp(-909.914/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[296] = q_f - q_r;

    /*reaction 298: O2 + CH3CHO => HO2 + CH3 + CO */
    phi_f = sc[3]*sc[52];
    k_f = 1e-06 * 3.01e+13*exp(-19703.1/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[297] = q_f - q_r;

    /*reaction 299: H + CH3CHO <=> CH2CHO + H2 */
    phi_f = sc[1]*sc[52];
    k_f = 1e-06 * 2.05e+09*exp(1.16*tc[0]-1210.37/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[51]*sc[0];
    Kc = exp((g_RT[1] + g_RT[52]) - (g_RT[51] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[298] = q_f - q_r;

    /*reaction 300: H + CH3CHO => CH3 + H2 + CO */
    phi_f = sc[1]*sc[52];
    k_f = 1e-06 * 2.05e+09*exp(1.16*tc[0]-1210.37/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[299] = q_f - q_r;

    /*reaction 301: OH + CH3CHO => CH3 + H2O + CO */
    phi_f = sc[4]*sc[52];
    k_f = 1e-06 * 2.343e+10*exp(0.73*tc[0]+560.141/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[300] = q_f - q_r;

    /*reaction 302: HO2 + CH3CHO => CH3 + H2O2 + CO */
    phi_f = sc[6]*sc[52];
    k_f = 1e-06 * 3.01e+12*exp(-6000.5/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[301] = q_f - q_r;

    /*reaction 303: CH3 + CH3CHO => CH3 + CH4 + CO */
    phi_f = sc[12]*sc[52];
    k_f = 1e-06 * 2.72e+06*exp(1.77*tc[0]-2979.37/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[302] = q_f - q_r;

    /*reaction 304: H + CH2CO (+M) <=> CH2CHO (+M) */
    phi_f = sc[1]*sc[28];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 4.865e+11*exp(0.422*tc[0]+883.241/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.012e+42*exp(-7.63*tc[0]-1939.61/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.535*exp(T/-201))+ (0.465*exp(T/-1773))+ (exp(-5333/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[51];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[28]) - (g_RT[51]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[303] = q_f - q_r;

    /*reaction 305: O + CH2CHO => H + CH2 + CO2 */
    phi_f = sc[2]*sc[51];
    k_f = 1e-06 * 1.5e+14;
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[304] = q_f - q_r;

    /*reaction 306: O2 + CH2CHO => OH + CO + CH2O */
    phi_f = sc[3]*sc[51];
    k_f = 1e-06 * 1.81e+10;
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[305] = q_f - q_r;

    /*reaction 307: O2 + CH2CHO => OH + 2 HCO */
    phi_f = sc[3]*sc[51];
    k_f = 1e-06 * 2.35e+10;
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[306] = q_f - q_r;

    /*reaction 308: H + CH2CHO <=> CH3 + HCO */
    phi_f = sc[1]*sc[51];
    k_f = 1e-06 * 2.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[16];
    Kc = exp((g_RT[1] + g_RT[51]) - (g_RT[12] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[307] = q_f - q_r;

    /*reaction 309: H + CH2CHO <=> CH2CO + H2 */
    phi_f = sc[1]*sc[51];
    k_f = 1e-06 * 1.1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[28]*sc[0];
    Kc = exp((g_RT[1] + g_RT[51]) - (g_RT[28] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[308] = q_f - q_r;

    /*reaction 310: OH + CH2CHO <=> H2O + CH2CO */
    phi_f = sc[4]*sc[51];
    k_f = 1e-06 * 1.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[28];
    Kc = exp((g_RT[4] + g_RT[51]) - (g_RT[5] + g_RT[28]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[309] = q_f - q_r;

    /*reaction 311: OH + CH2CHO <=> HCO + CH2OH */
    phi_f = sc[4]*sc[51];
    k_f = 1e-06 * 3.01e+13;
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[18];
    Kc = exp((g_RT[4] + g_RT[51]) - (g_RT[16] + g_RT[18]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[310] = q_f - q_r;

    /*reaction 312: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    phi_f = sc[12]*sc[25];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 9.43e+12;
    redP = 1e-12 * alpha / k_f * 2.71e+74*exp(-16.82*tc[0]-6575.24/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.8473*exp(T/-291))+ (0.1527*exp(T/-2742))+ (exp(-7748/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[50];
    Kc = 1.0 / (refC) * exp((g_RT[12] + g_RT[25]) - (g_RT[50]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[311] = q_f - q_r;

    /*reaction 313: O + C3H8 <=> OH + C3H7 */
    phi_f = sc[2]*sc[50];
    k_f = 1e-06 * 193000*exp(2.68*tc[0]-1870.16/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[49];
    Kc = exp((g_RT[2] + g_RT[50]) - (g_RT[4] + g_RT[49]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[312] = q_f - q_r;

    /*reaction 314: H + C3H8 <=> C3H7 + H2 */
    phi_f = sc[1]*sc[50];
    k_f = 1e-06 * 1.32e+06*exp(2.54*tc[0]-3400.1/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[49]*sc[0];
    Kc = exp((g_RT[1] + g_RT[50]) - (g_RT[49] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[313] = q_f - q_r;

    /*reaction 315: OH + C3H8 <=> C3H7 + H2O */
    phi_f = sc[4]*sc[50];
    k_f = 1e-06 * 3.16e+07*exp(1.8*tc[0]-470.055/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[49]*sc[5];
    Kc = exp((g_RT[4] + g_RT[50]) - (g_RT[49] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[314] = q_f - q_r;

    /*reaction 316: C3H7 + H2O2 <=> HO2 + C3H8 */
    phi_f = sc[49]*sc[7];
    k_f = 1e-06 * 378*exp(2.72*tc[0]-754.907/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[50];
    Kc = exp((g_RT[49] + g_RT[7]) - (g_RT[6] + g_RT[50]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[315] = q_f - q_r;

    /*reaction 317: CH3 + C3H8 <=> C3H7 + CH4 */
    phi_f = sc[12]*sc[50];
    k_f = 1e-06 * 0.903*exp(3.65*tc[0]-3600.4/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[49]*sc[13];
    Kc = exp((g_RT[12] + g_RT[50]) - (g_RT[49] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[316] = q_f - q_r;

    /*reaction 318: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    phi_f = sc[12]*sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 2.55e+06*exp(1.6*tc[0]-2868.65/tc[1]);
    redP = 1e-12 * alpha / k_f * 3e+63*exp(-14.6*tc[0]-9144.44/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.8106*exp(T/-277))+ (0.1894*exp(T/-8748))+ (exp(-7891/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[49];
    Kc = 1.0 / (refC) * exp((g_RT[12] + g_RT[24]) - (g_RT[49]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[317] = q_f - q_r;

    /*reaction 319: O + C3H7 <=> C2H5 + CH2O */
    phi_f = sc[2]*sc[49];
    k_f = 1e-06 * 9.64e+13;
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[17];
    Kc = exp((g_RT[2] + g_RT[49]) - (g_RT[25] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[318] = q_f - q_r;

    /*reaction 320: H + C3H7 (+M) <=> C3H8 (+M) */
    phi_f = sc[1]*sc[49];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[48];
    k_f = 1e-06 * 3.613e+13;
    redP = 1e-12 * alpha / k_f * 4.42e+61*exp(-13.545*tc[0]-5715.65/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.685*exp(T/-369))+ (0.315*exp(T/-3285))+ (exp(-6667/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[50];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[49]) - (g_RT[50]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[319] = q_f - q_r;

    /*reaction 321: H + C3H7 <=> CH3 + C2H5 */
    phi_f = sc[1]*sc[49];
    k_f = 1e-06 * 4.06e+06*exp(2.19*tc[0]-447.911/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[25];
    Kc = exp((g_RT[1] + g_RT[49]) - (g_RT[12] + g_RT[25]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[320] = q_f - q_r;

    /*reaction 322: OH + C3H7 <=> C2H5 + CH2OH */
    phi_f = sc[4]*sc[49];
    k_f = 1e-06 * 2.41e+13;
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[18];
    Kc = exp((g_RT[4] + g_RT[49]) - (g_RT[25] + g_RT[18]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[321] = q_f - q_r;

    /*reaction 323: HO2 + C3H7 <=> O2 + C3H8 */
    phi_f = sc[6]*sc[49];
    k_f = 1e-06 * 2.55e+10*exp(0.255*tc[0]+474.585/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[50];
    Kc = exp((g_RT[6] + g_RT[49]) - (g_RT[3] + g_RT[50]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[322] = q_f - q_r;

    /*reaction 324: HO2 + C3H7 => OH + C2H5 + CH2O */
    phi_f = sc[6]*sc[49];
    k_f = 1e-06 * 2.41e+13;
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[323] = q_f - q_r;

    /*reaction 325: CH3 + C3H7 <=> 2 C2H5 */
    phi_f = sc[12]*sc[49];
    k_f = 1e-06 * 1.927e+13*exp(-0.32*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[25];
    Kc = exp((g_RT[12] + g_RT[49]) - (2 * g_RT[25]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[324] = q_f - q_r;

    return;
}


/*compute the equilibrium constants for each reaction */
void equilibriumConstants(double *kc, double * g_RT, double T)
{
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

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

    /*reaction 12: O + CO (+M) <=> CO2 (+M) */
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
    kc[35] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[47]) - (g_RT[6] + g_RT[47]));

    /*reaction 37: H + O2 + AR <=> HO2 + AR */
    kc[36] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[48]) - (g_RT[6] + g_RT[48]));

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

    /*reaction 135: CH2 + O2 => OH + H + CO */
    kc[134] = refC * exp((g_RT[10] + g_RT[3]) - (g_RT[4] + g_RT[1] + g_RT[14]));

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
    kc[141] = exp((g_RT[11] + g_RT[47]) - (g_RT[10] + g_RT[47]));

    /*reaction 143: CH2(S) + AR <=> CH2 + AR */
    kc[142] = exp((g_RT[11] + g_RT[48]) - (g_RT[10] + g_RT[48]));

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

    /*reaction 178: N + NO <=> N2 + O */
    kc[177] = exp((g_RT[30] + g_RT[35]) - (g_RT[47] + g_RT[2]));

    /*reaction 179: N + O2 <=> NO + O */
    kc[178] = exp((g_RT[30] + g_RT[3]) - (g_RT[35] + g_RT[2]));

    /*reaction 180: N + OH <=> NO + H */
    kc[179] = exp((g_RT[30] + g_RT[4]) - (g_RT[35] + g_RT[1]));

    /*reaction 181: N2O + O <=> N2 + O2 */
    kc[180] = exp((g_RT[37] + g_RT[2]) - (g_RT[47] + g_RT[3]));

    /*reaction 182: N2O + O <=> 2 NO */
    kc[181] = exp((g_RT[37] + g_RT[2]) - (2 * g_RT[35]));

    /*reaction 183: N2O + H <=> N2 + OH */
    kc[182] = exp((g_RT[37] + g_RT[1]) - (g_RT[47] + g_RT[4]));

    /*reaction 184: N2O + OH <=> N2 + HO2 */
    kc[183] = exp((g_RT[37] + g_RT[4]) - (g_RT[47] + g_RT[6]));

    /*reaction 185: N2O (+M) <=> N2 + O (+M) */
    kc[184] = refC * exp((g_RT[37]) - (g_RT[47] + g_RT[2]));

    /*reaction 186: HO2 + NO <=> NO2 + OH */
    kc[185] = exp((g_RT[6] + g_RT[35]) - (g_RT[36] + g_RT[4]));

    /*reaction 187: NO + O + M <=> NO2 + M */
    kc[186] = 1.0 / (refC) * exp((g_RT[35] + g_RT[2]) - (g_RT[36]));

    /*reaction 188: NO2 + O <=> NO + O2 */
    kc[187] = exp((g_RT[36] + g_RT[2]) - (g_RT[35] + g_RT[3]));

    /*reaction 189: NO2 + H <=> NO + OH */
    kc[188] = exp((g_RT[36] + g_RT[1]) - (g_RT[35] + g_RT[4]));

    /*reaction 190: NH + O <=> NO + H */
    kc[189] = exp((g_RT[31] + g_RT[2]) - (g_RT[35] + g_RT[1]));

    /*reaction 191: NH + H <=> N + H2 */
    kc[190] = exp((g_RT[31] + g_RT[1]) - (g_RT[30] + g_RT[0]));

    /*reaction 192: NH + OH <=> HNO + H */
    kc[191] = exp((g_RT[31] + g_RT[4]) - (g_RT[38] + g_RT[1]));

    /*reaction 193: NH + OH <=> N + H2O */
    kc[192] = exp((g_RT[31] + g_RT[4]) - (g_RT[30] + g_RT[5]));

    /*reaction 194: NH + O2 <=> HNO + O */
    kc[193] = exp((g_RT[31] + g_RT[3]) - (g_RT[38] + g_RT[2]));

    /*reaction 195: NH + O2 <=> NO + OH */
    kc[194] = exp((g_RT[31] + g_RT[3]) - (g_RT[35] + g_RT[4]));

    /*reaction 196: NH + N <=> N2 + H */
    kc[195] = exp((g_RT[31] + g_RT[30]) - (g_RT[47] + g_RT[1]));

    /*reaction 197: NH + H2O <=> HNO + H2 */
    kc[196] = exp((g_RT[31] + g_RT[5]) - (g_RT[38] + g_RT[0]));

    /*reaction 198: NH + NO <=> N2 + OH */
    kc[197] = exp((g_RT[31] + g_RT[35]) - (g_RT[47] + g_RT[4]));

    /*reaction 199: NH + NO <=> N2O + H */
    kc[198] = exp((g_RT[31] + g_RT[35]) - (g_RT[37] + g_RT[1]));

    /*reaction 200: NH2 + O <=> OH + NH */
    kc[199] = exp((g_RT[32] + g_RT[2]) - (g_RT[4] + g_RT[31]));

    /*reaction 201: NH2 + O <=> H + HNO */
    kc[200] = exp((g_RT[32] + g_RT[2]) - (g_RT[1] + g_RT[38]));

    /*reaction 202: NH2 + H <=> NH + H2 */
    kc[201] = exp((g_RT[32] + g_RT[1]) - (g_RT[31] + g_RT[0]));

    /*reaction 203: NH2 + OH <=> NH + H2O */
    kc[202] = exp((g_RT[32] + g_RT[4]) - (g_RT[31] + g_RT[5]));

    /*reaction 204: NNH <=> N2 + H */
    kc[203] = refC * exp((g_RT[34]) - (g_RT[47] + g_RT[1]));

    /*reaction 205: NNH + M <=> N2 + H + M */
    kc[204] = refC * exp((g_RT[34]) - (g_RT[47] + g_RT[1]));

    /*reaction 206: NNH + O2 <=> HO2 + N2 */
    kc[205] = exp((g_RT[34] + g_RT[3]) - (g_RT[6] + g_RT[47]));

    /*reaction 207: NNH + O <=> OH + N2 */
    kc[206] = exp((g_RT[34] + g_RT[2]) - (g_RT[4] + g_RT[47]));

    /*reaction 208: NNH + O <=> NH + NO */
    kc[207] = exp((g_RT[34] + g_RT[2]) - (g_RT[31] + g_RT[35]));

    /*reaction 209: NNH + H <=> H2 + N2 */
    kc[208] = exp((g_RT[34] + g_RT[1]) - (g_RT[0] + g_RT[47]));

    /*reaction 210: NNH + OH <=> H2O + N2 */
    kc[209] = exp((g_RT[34] + g_RT[4]) - (g_RT[5] + g_RT[47]));

    /*reaction 211: NNH + CH3 <=> CH4 + N2 */
    kc[210] = exp((g_RT[34] + g_RT[12]) - (g_RT[13] + g_RT[47]));

    /*reaction 212: H + NO + M <=> HNO + M */
    kc[211] = 1.0 / (refC) * exp((g_RT[1] + g_RT[35]) - (g_RT[38]));

    /*reaction 213: HNO + O <=> NO + OH */
    kc[212] = exp((g_RT[38] + g_RT[2]) - (g_RT[35] + g_RT[4]));

    /*reaction 214: HNO + H <=> H2 + NO */
    kc[213] = exp((g_RT[38] + g_RT[1]) - (g_RT[0] + g_RT[35]));

    /*reaction 215: HNO + OH <=> NO + H2O */
    kc[214] = exp((g_RT[38] + g_RT[4]) - (g_RT[35] + g_RT[5]));

    /*reaction 216: HNO + O2 <=> HO2 + NO */
    kc[215] = exp((g_RT[38] + g_RT[3]) - (g_RT[6] + g_RT[35]));

    /*reaction 217: CN + O <=> CO + N */
    kc[216] = exp((g_RT[39] + g_RT[2]) - (g_RT[14] + g_RT[30]));

    /*reaction 218: CN + OH <=> NCO + H */
    kc[217] = exp((g_RT[39] + g_RT[4]) - (g_RT[46] + g_RT[1]));

    /*reaction 219: CN + H2O <=> HCN + OH */
    kc[218] = exp((g_RT[39] + g_RT[5]) - (g_RT[40] + g_RT[4]));

    /*reaction 220: CN + O2 <=> NCO + O */
    kc[219] = exp((g_RT[39] + g_RT[3]) - (g_RT[46] + g_RT[2]));

    /*reaction 221: CN + H2 <=> HCN + H */
    kc[220] = exp((g_RT[39] + g_RT[0]) - (g_RT[40] + g_RT[1]));

    /*reaction 222: NCO + O <=> NO + CO */
    kc[221] = exp((g_RT[46] + g_RT[2]) - (g_RT[35] + g_RT[14]));

    /*reaction 223: NCO + H <=> NH + CO */
    kc[222] = exp((g_RT[46] + g_RT[1]) - (g_RT[31] + g_RT[14]));

    /*reaction 224: NCO + OH <=> NO + H + CO */
    kc[223] = refC * exp((g_RT[46] + g_RT[4]) - (g_RT[35] + g_RT[1] + g_RT[14]));

    /*reaction 225: NCO + N <=> N2 + CO */
    kc[224] = exp((g_RT[46] + g_RT[30]) - (g_RT[47] + g_RT[14]));

    /*reaction 226: NCO + O2 <=> NO + CO2 */
    kc[225] = exp((g_RT[46] + g_RT[3]) - (g_RT[35] + g_RT[15]));

    /*reaction 227: NCO + M <=> N + CO + M */
    kc[226] = refC * exp((g_RT[46]) - (g_RT[30] + g_RT[14]));

    /*reaction 228: NCO + NO <=> N2O + CO */
    kc[227] = exp((g_RT[46] + g_RT[35]) - (g_RT[37] + g_RT[14]));

    /*reaction 229: NCO + NO <=> N2 + CO2 */
    kc[228] = exp((g_RT[46] + g_RT[35]) - (g_RT[47] + g_RT[15]));

    /*reaction 230: HCN + M <=> H + CN + M */
    kc[229] = refC * exp((g_RT[40]) - (g_RT[1] + g_RT[39]));

    /*reaction 231: HCN + O <=> NCO + H */
    kc[230] = exp((g_RT[40] + g_RT[2]) - (g_RT[46] + g_RT[1]));

    /*reaction 232: HCN + O <=> NH + CO */
    kc[231] = exp((g_RT[40] + g_RT[2]) - (g_RT[31] + g_RT[14]));

    /*reaction 233: HCN + O <=> CN + OH */
    kc[232] = exp((g_RT[40] + g_RT[2]) - (g_RT[39] + g_RT[4]));

    /*reaction 234: HCN + OH <=> HOCN + H */
    kc[233] = exp((g_RT[40] + g_RT[4]) - (g_RT[44] + g_RT[1]));

    /*reaction 235: HCN + OH <=> HNCO + H */
    kc[234] = exp((g_RT[40] + g_RT[4]) - (g_RT[45] + g_RT[1]));

    /*reaction 236: HCN + OH <=> NH2 + CO */
    kc[235] = exp((g_RT[40] + g_RT[4]) - (g_RT[32] + g_RT[14]));

    /*reaction 237: H + HCN (+M) <=> H2CN (+M) */
    kc[236] = 1.0 / (refC) * exp((g_RT[1] + g_RT[40]) - (g_RT[41]));

    /*reaction 238: H2CN + N <=> N2 + CH2 */
    kc[237] = exp((g_RT[41] + g_RT[30]) - (g_RT[47] + g_RT[10]));

    /*reaction 239: C + N2 <=> CN + N */
    kc[238] = exp((g_RT[8] + g_RT[47]) - (g_RT[39] + g_RT[30]));

    /*reaction 240: CH + N2 <=> HCN + N */
    kc[239] = exp((g_RT[9] + g_RT[47]) - (g_RT[40] + g_RT[30]));

    /*reaction 241: CH + N2 (+M) <=> HCNN (+M) */
    kc[240] = 1.0 / (refC) * exp((g_RT[9] + g_RT[47]) - (g_RT[42]));

    /*reaction 242: CH2 + N2 <=> HCN + NH */
    kc[241] = exp((g_RT[10] + g_RT[47]) - (g_RT[40] + g_RT[31]));

    /*reaction 243: CH2(S) + N2 <=> NH + HCN */
    kc[242] = exp((g_RT[11] + g_RT[47]) - (g_RT[31] + g_RT[40]));

    /*reaction 244: C + NO <=> CN + O */
    kc[243] = exp((g_RT[8] + g_RT[35]) - (g_RT[39] + g_RT[2]));

    /*reaction 245: C + NO <=> CO + N */
    kc[244] = exp((g_RT[8] + g_RT[35]) - (g_RT[14] + g_RT[30]));

    /*reaction 246: CH + NO <=> HCN + O */
    kc[245] = exp((g_RT[9] + g_RT[35]) - (g_RT[40] + g_RT[2]));

    /*reaction 247: CH + NO <=> H + NCO */
    kc[246] = exp((g_RT[9] + g_RT[35]) - (g_RT[1] + g_RT[46]));

    /*reaction 248: CH + NO <=> N + HCO */
    kc[247] = exp((g_RT[9] + g_RT[35]) - (g_RT[30] + g_RT[16]));

    /*reaction 249: CH2 + NO <=> H + HNCO */
    kc[248] = exp((g_RT[10] + g_RT[35]) - (g_RT[1] + g_RT[45]));

    /*reaction 250: CH2 + NO <=> OH + HCN */
    kc[249] = exp((g_RT[10] + g_RT[35]) - (g_RT[4] + g_RT[40]));

    /*reaction 251: CH2 + NO <=> H + HCNO */
    kc[250] = exp((g_RT[10] + g_RT[35]) - (g_RT[1] + g_RT[43]));

    /*reaction 252: CH2(S) + NO <=> H + HNCO */
    kc[251] = exp((g_RT[11] + g_RT[35]) - (g_RT[1] + g_RT[45]));

    /*reaction 253: CH2(S) + NO <=> OH + HCN */
    kc[252] = exp((g_RT[11] + g_RT[35]) - (g_RT[4] + g_RT[40]));

    /*reaction 254: CH2(S) + NO <=> H + HCNO */
    kc[253] = exp((g_RT[11] + g_RT[35]) - (g_RT[1] + g_RT[43]));

    /*reaction 255: CH3 + NO <=> HCN + H2O */
    kc[254] = exp((g_RT[12] + g_RT[35]) - (g_RT[40] + g_RT[5]));

    /*reaction 256: CH3 + NO <=> H2CN + OH */
    kc[255] = exp((g_RT[12] + g_RT[35]) - (g_RT[41] + g_RT[4]));

    /*reaction 257: HCNN + O <=> CO + H + N2 */
    kc[256] = refC * exp((g_RT[42] + g_RT[2]) - (g_RT[14] + g_RT[1] + g_RT[47]));

    /*reaction 258: HCNN + O <=> HCN + NO */
    kc[257] = exp((g_RT[42] + g_RT[2]) - (g_RT[40] + g_RT[35]));

    /*reaction 259: HCNN + O2 <=> O + HCO + N2 */
    kc[258] = refC * exp((g_RT[42] + g_RT[3]) - (g_RT[2] + g_RT[16] + g_RT[47]));

    /*reaction 260: HCNN + OH <=> H + HCO + N2 */
    kc[259] = refC * exp((g_RT[42] + g_RT[4]) - (g_RT[1] + g_RT[16] + g_RT[47]));

    /*reaction 261: HCNN + H <=> CH2 + N2 */
    kc[260] = exp((g_RT[42] + g_RT[1]) - (g_RT[10] + g_RT[47]));

    /*reaction 262: HNCO + O <=> NH + CO2 */
    kc[261] = exp((g_RT[45] + g_RT[2]) - (g_RT[31] + g_RT[15]));

    /*reaction 263: HNCO + O <=> HNO + CO */
    kc[262] = exp((g_RT[45] + g_RT[2]) - (g_RT[38] + g_RT[14]));

    /*reaction 264: HNCO + O <=> NCO + OH */
    kc[263] = exp((g_RT[45] + g_RT[2]) - (g_RT[46] + g_RT[4]));

    /*reaction 265: HNCO + H <=> NH2 + CO */
    kc[264] = exp((g_RT[45] + g_RT[1]) - (g_RT[32] + g_RT[14]));

    /*reaction 266: HNCO + H <=> H2 + NCO */
    kc[265] = exp((g_RT[45] + g_RT[1]) - (g_RT[0] + g_RT[46]));

    /*reaction 267: HNCO + OH <=> NCO + H2O */
    kc[266] = exp((g_RT[45] + g_RT[4]) - (g_RT[46] + g_RT[5]));

    /*reaction 268: HNCO + OH <=> NH2 + CO2 */
    kc[267] = exp((g_RT[45] + g_RT[4]) - (g_RT[32] + g_RT[15]));

    /*reaction 269: HNCO + M <=> NH + CO + M */
    kc[268] = refC * exp((g_RT[45]) - (g_RT[31] + g_RT[14]));

    /*reaction 270: HCNO + H <=> H + HNCO */
    kc[269] = exp((g_RT[43] + g_RT[1]) - (g_RT[1] + g_RT[45]));

    /*reaction 271: HCNO + H <=> OH + HCN */
    kc[270] = exp((g_RT[43] + g_RT[1]) - (g_RT[4] + g_RT[40]));

    /*reaction 272: HCNO + H <=> NH2 + CO */
    kc[271] = exp((g_RT[43] + g_RT[1]) - (g_RT[32] + g_RT[14]));

    /*reaction 273: HOCN + H <=> H + HNCO */
    kc[272] = exp((g_RT[44] + g_RT[1]) - (g_RT[1] + g_RT[45]));

    /*reaction 274: HCCO + NO <=> HCNO + CO */
    kc[273] = exp((g_RT[27] + g_RT[35]) - (g_RT[43] + g_RT[14]));

    /*reaction 275: CH3 + N <=> H2CN + H */
    kc[274] = exp((g_RT[12] + g_RT[30]) - (g_RT[41] + g_RT[1]));

    /*reaction 276: CH3 + N <=> HCN + H2 */
    kc[275] = exp((g_RT[12] + g_RT[30]) - (g_RT[40] + g_RT[0]));

    /*reaction 277: NH3 + H <=> NH2 + H2 */
    kc[276] = exp((g_RT[33] + g_RT[1]) - (g_RT[32] + g_RT[0]));

    /*reaction 278: NH3 + OH <=> NH2 + H2O */
    kc[277] = exp((g_RT[33] + g_RT[4]) - (g_RT[32] + g_RT[5]));

    /*reaction 279: NH3 + O <=> NH2 + OH */
    kc[278] = exp((g_RT[33] + g_RT[2]) - (g_RT[32] + g_RT[4]));

    /*reaction 280: NH + CO2 <=> HNO + CO */
    kc[279] = exp((g_RT[31] + g_RT[15]) - (g_RT[38] + g_RT[14]));

    /*reaction 281: CN + NO2 <=> NCO + NO */
    kc[280] = exp((g_RT[39] + g_RT[36]) - (g_RT[46] + g_RT[35]));

    /*reaction 282: NCO + NO2 <=> N2O + CO2 */
    kc[281] = exp((g_RT[46] + g_RT[36]) - (g_RT[37] + g_RT[15]));

    /*reaction 283: N + CO2 <=> NO + CO */
    kc[282] = exp((g_RT[30] + g_RT[15]) - (g_RT[35] + g_RT[14]));

    /*reaction 284: O + CH3 => H + H2 + CO */
    kc[283] = refC * exp((g_RT[2] + g_RT[12]) - (g_RT[1] + g_RT[0] + g_RT[14]));

    /*reaction 285: O + C2H4 <=> H + CH2CHO */
    kc[284] = exp((g_RT[2] + g_RT[24]) - (g_RT[1] + g_RT[51]));

    /*reaction 286: O + C2H5 <=> H + CH3CHO */
    kc[285] = exp((g_RT[2] + g_RT[25]) - (g_RT[1] + g_RT[52]));

    /*reaction 287: OH + HO2 <=> O2 + H2O */
    kc[286] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));

    /*reaction 288: OH + CH3 => H2 + CH2O */
    kc[287] = exp((g_RT[4] + g_RT[12]) - (g_RT[0] + g_RT[17]));

    /*reaction 289: CH + H2 (+M) <=> CH3 (+M) */
    kc[288] = 1.0 / (refC) * exp((g_RT[9] + g_RT[0]) - (g_RT[12]));

    /*reaction 290: CH2 + O2 => 2 H + CO2 */
    kc[289] = refC * exp((g_RT[10] + g_RT[3]) - (2 * g_RT[1] + g_RT[15]));

    /*reaction 291: CH2 + O2 <=> O + CH2O */
    kc[290] = exp((g_RT[10] + g_RT[3]) - (g_RT[2] + g_RT[17]));

    /*reaction 292: CH2 + CH2 => 2 H + C2H2 */
    kc[291] = refC * exp((g_RT[10] + g_RT[10]) - (2 * g_RT[1] + g_RT[22]));

    /*reaction 293: CH2(S) + H2O => H2 + CH2O */
    kc[292] = exp((g_RT[11] + g_RT[5]) - (g_RT[0] + g_RT[17]));

    /*reaction 294: C2H3 + O2 <=> O + CH2CHO */
    kc[293] = exp((g_RT[23] + g_RT[3]) - (g_RT[2] + g_RT[51]));

    /*reaction 295: C2H3 + O2 <=> HO2 + C2H2 */
    kc[294] = exp((g_RT[23] + g_RT[3]) - (g_RT[6] + g_RT[22]));

    /*reaction 296: O + CH3CHO <=> OH + CH2CHO */
    kc[295] = exp((g_RT[2] + g_RT[52]) - (g_RT[4] + g_RT[51]));

    /*reaction 297: O + CH3CHO => OH + CH3 + CO */
    kc[296] = refC * exp((g_RT[2] + g_RT[52]) - (g_RT[4] + g_RT[12] + g_RT[14]));

    /*reaction 298: O2 + CH3CHO => HO2 + CH3 + CO */
    kc[297] = refC * exp((g_RT[3] + g_RT[52]) - (g_RT[6] + g_RT[12] + g_RT[14]));

    /*reaction 299: H + CH3CHO <=> CH2CHO + H2 */
    kc[298] = exp((g_RT[1] + g_RT[52]) - (g_RT[51] + g_RT[0]));

    /*reaction 300: H + CH3CHO => CH3 + H2 + CO */
    kc[299] = refC * exp((g_RT[1] + g_RT[52]) - (g_RT[12] + g_RT[0] + g_RT[14]));

    /*reaction 301: OH + CH3CHO => CH3 + H2O + CO */
    kc[300] = refC * exp((g_RT[4] + g_RT[52]) - (g_RT[12] + g_RT[5] + g_RT[14]));

    /*reaction 302: HO2 + CH3CHO => CH3 + H2O2 + CO */
    kc[301] = refC * exp((g_RT[6] + g_RT[52]) - (g_RT[12] + g_RT[7] + g_RT[14]));

    /*reaction 303: CH3 + CH3CHO => CH3 + CH4 + CO */
    kc[302] = refC * exp((g_RT[12] + g_RT[52]) - (g_RT[12] + g_RT[13] + g_RT[14]));

    /*reaction 304: H + CH2CO (+M) <=> CH2CHO (+M) */
    kc[303] = 1.0 / (refC) * exp((g_RT[1] + g_RT[28]) - (g_RT[51]));

    /*reaction 305: O + CH2CHO => H + CH2 + CO2 */
    kc[304] = refC * exp((g_RT[2] + g_RT[51]) - (g_RT[1] + g_RT[10] + g_RT[15]));

    /*reaction 306: O2 + CH2CHO => OH + CO + CH2O */
    kc[305] = refC * exp((g_RT[3] + g_RT[51]) - (g_RT[4] + g_RT[14] + g_RT[17]));

    /*reaction 307: O2 + CH2CHO => OH + 2 HCO */
    kc[306] = refC * exp((g_RT[3] + g_RT[51]) - (g_RT[4] + 2 * g_RT[16]));

    /*reaction 308: H + CH2CHO <=> CH3 + HCO */
    kc[307] = exp((g_RT[1] + g_RT[51]) - (g_RT[12] + g_RT[16]));

    /*reaction 309: H + CH2CHO <=> CH2CO + H2 */
    kc[308] = exp((g_RT[1] + g_RT[51]) - (g_RT[28] + g_RT[0]));

    /*reaction 310: OH + CH2CHO <=> H2O + CH2CO */
    kc[309] = exp((g_RT[4] + g_RT[51]) - (g_RT[5] + g_RT[28]));

    /*reaction 311: OH + CH2CHO <=> HCO + CH2OH */
    kc[310] = exp((g_RT[4] + g_RT[51]) - (g_RT[16] + g_RT[18]));

    /*reaction 312: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    kc[311] = 1.0 / (refC) * exp((g_RT[12] + g_RT[25]) - (g_RT[50]));

    /*reaction 313: O + C3H8 <=> OH + C3H7 */
    kc[312] = exp((g_RT[2] + g_RT[50]) - (g_RT[4] + g_RT[49]));

    /*reaction 314: H + C3H8 <=> C3H7 + H2 */
    kc[313] = exp((g_RT[1] + g_RT[50]) - (g_RT[49] + g_RT[0]));

    /*reaction 315: OH + C3H8 <=> C3H7 + H2O */
    kc[314] = exp((g_RT[4] + g_RT[50]) - (g_RT[49] + g_RT[5]));

    /*reaction 316: C3H7 + H2O2 <=> HO2 + C3H8 */
    kc[315] = exp((g_RT[49] + g_RT[7]) - (g_RT[6] + g_RT[50]));

    /*reaction 317: CH3 + C3H8 <=> C3H7 + CH4 */
    kc[316] = exp((g_RT[12] + g_RT[50]) - (g_RT[49] + g_RT[13]));

    /*reaction 318: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    kc[317] = 1.0 / (refC) * exp((g_RT[12] + g_RT[24]) - (g_RT[49]));

    /*reaction 319: O + C3H7 <=> C2H5 + CH2O */
    kc[318] = exp((g_RT[2] + g_RT[49]) - (g_RT[25] + g_RT[17]));

    /*reaction 320: H + C3H7 (+M) <=> C3H8 (+M) */
    kc[319] = 1.0 / (refC) * exp((g_RT[1] + g_RT[49]) - (g_RT[50]));

    /*reaction 321: H + C3H7 <=> CH3 + C2H5 */
    kc[320] = exp((g_RT[1] + g_RT[49]) - (g_RT[12] + g_RT[25]));

    /*reaction 322: OH + C3H7 <=> C2H5 + CH2OH */
    kc[321] = exp((g_RT[4] + g_RT[49]) - (g_RT[25] + g_RT[18]));

    /*reaction 323: HO2 + C3H7 <=> O2 + C3H8 */
    kc[322] = exp((g_RT[6] + g_RT[49]) - (g_RT[3] + g_RT[50]));

    /*reaction 324: HO2 + C3H7 => OH + C2H5 + CH2O */
    kc[323] = refC * exp((g_RT[6] + g_RT[49]) - (g_RT[4] + g_RT[25] + g_RT[17]));

    /*reaction 325: CH3 + C3H7 <=> 2 C2H5 */
    kc[324] = exp((g_RT[12] + g_RT[49]) - (2 * g_RT[25]));

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
        /*species 30: N */
        species[30] =
            +5.61046370e+04 / tc[1]
            -1.69390870e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 31: NH */
        species[31] =
            +4.18806290e+04 / tc[1]
            +1.64458070e+00
            -3.49290850e+00 * tc[0]
            -1.55895990e-04 * tc[1]
            +2.48174733e-07 * tc[2]
            -2.06803683e-10 * tc[3]
            +5.17848350e-14 * tc[4];
        /*species 32: NH2 */
        species[32] =
            +2.18859100e+04 / tc[1]
            +4.34584538e+00
            -4.20400290e+00 * tc[0]
            +1.05306925e-03 * tc[1]
            -1.18447247e-06 * tc[2]
            +4.67626642e-10 * tc[3]
            -8.22035850e-14 * tc[4];
        /*species 33: NH3 */
        species[33] =
            -6.74172850e+03 / tc[1]
            +4.91140017e+00
            -4.28602740e+00 * tc[0]
            +2.33026150e-03 * tc[1]
            -3.61975217e-06 * tc[2]
            +1.90074058e-09 * tc[3]
            -4.13190230e-13 * tc[4];
        /*species 34: NNH */
        species[34] =
            +2.87919730e+04 / tc[1]
            +1.36675170e+00
            -4.34469270e+00 * tc[0]
            +2.42485360e-03 * tc[1]
            -3.34324317e-06 * tc[2]
            +1.81053867e-09 * tc[3]
            -3.97347695e-13 * tc[4];
        /*species 35: NO */
        species[35] =
            +9.84462300e+03 / tc[1]
            +1.93762990e+00
            -4.21847630e+00 * tc[0]
            +2.31948800e-03 * tc[1]
            -1.84017033e-06 * tc[2]
            +7.78011283e-10 * tc[3]
            -1.40178850e-13 * tc[4];
        /*species 36: NO2 */
        species[36] =
            +2.89661790e+03 / tc[1]
            -2.36796050e+00
            -3.94403120e+00 * tc[0]
            +7.92714500e-04 * tc[1]
            -2.77630200e-06 * tc[2]
            +1.70628550e-09 * tc[3]
            -3.91752820e-13 * tc[4];
        /*species 37: N2O */
        species[37] =
            +8.74177440e+03 / tc[1]
            -8.50084180e+00
            -2.25715020e+00 * tc[0]
            -5.65236400e-03 * tc[1]
            +2.27855317e-06 * tc[2]
            -8.06831717e-10 * tc[3]
            +1.46535910e-13 * tc[4];
        /*species 38: HNO */
        species[38] =
            +1.15482970e+04 / tc[1]
            +2.78364990e+00
            -4.53349160e+00 * tc[0]
            +2.83480855e-03 * tc[1]
            -3.07886783e-06 * tc[2]
            +1.42809117e-09 * tc[3]
            -2.77272865e-13 * tc[4];
        /*species 39: CN */
        species[39] =
            +5.17083400e+04 / tc[1]
            -3.67564400e-01
            -3.61293510e+00 * tc[0]
            +4.77756635e-04 * tc[1]
            -3.57382950e-07 * tc[2]
            +2.62636025e-11 * tc[3]
            +2.32151780e-14 * tc[4];
        /*species 40: HCN */
        species[40] =
            +1.47126330e+04 / tc[1]
            -6.65745330e+00
            -2.25898860e+00 * tc[0]
            -5.02558500e-03 * tc[1]
            +2.22529383e-06 * tc[2]
            -8.41029083e-10 * tc[3]
            +1.50445140e-13 * tc[4];
        /*species 41: H2CN */
        species[41] =
            +2.86378200e+04 / tc[1]
            -6.14109010e+00
            -2.85166100e+00 * tc[0]
            -2.84761655e-03 * tc[1]
            -1.78523333e-07 * tc[2]
            +1.35217667e-10 * tc[3]
            +1.17555405e-14 * tc[4];
        /*species 42: HCNN */
        species[42] =
            +5.42619840e+04 / tc[1]
            -9.15155060e+00
            -2.52431940e+00 * tc[0]
            -7.98030950e-03 * tc[1]
            +3.13605900e-06 * tc[2]
            -1.01046167e-09 * tc[3]
            +1.61786890e-13 * tc[4];
        /*species 46: NCO */
        species[46] =
            +1.46824770e+04 / tc[1]
            -6.72353380e+00
            -2.82693080e+00 * tc[0]
            -4.40258440e-03 * tc[1]
            +1.39776890e-06 * tc[2]
            -4.00141367e-10 * tc[3]
            +6.65679750e-14 * tc[4];
        /*species 47: N2 */
        species[47] =
            -1.02089990e+03 / tc[1]
            -6.51695000e-01
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
        /*species 48: AR */
        species[48] =
            -7.45375000e+02 / tc[1]
            -1.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 49: C3H7 */
        species[49] =
            +1.06318630e+04 / tc[1]
            -2.00710072e+01
            -1.05155180e+00 * tc[0]
            -1.29959900e-02 * tc[1]
            -3.96675667e-07 * tc[2]
            +1.63413075e-09 * tc[3]
            -4.68662350e-13 * tc[4];
        /*species 50: C3H8 */
        species[50] =
            -1.39585200e+04 / tc[1]
            -1.82681372e+01
            -9.33553810e-01 * tc[0]
            -1.32122895e-02 * tc[1]
            -1.01766212e-06 * tc[2]
            +1.83145825e-09 * tc[3]
            -4.75746265e-13 * tc[4];
        /*species 51: CH2CHO */
        species[51] =
            +1.52147660e+03 / tc[1]
            -6.14922800e+00
            -3.40906200e+00 * tc[0]
            -5.36928700e-03 * tc[1]
            -3.15248667e-07 * tc[2]
            +5.96548583e-10 * tc[3]
            -1.43369250e-13 * tc[4];
        /*species 52: CH3CHO */
        species[52] =
            -2.15728780e+04 / tc[1]
            +6.26443600e-01
            -4.72945950e+00 * tc[0]
            +1.59664290e-03 * tc[1]
            -7.92248683e-06 * tc[2]
            +4.78821758e-09 * tc[3]
            -1.09655560e-12 * tc[4];
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
        /*species 30: N */
        species[30] =
            +5.61337730e+04 / tc[1]
            -2.23366670e+00
            -2.41594290e+00 * tc[0]
            -8.74453250e-05 * tc[1]
            +1.98372817e-08 * tc[2]
            -2.51885375e-12 * tc[3]
            +1.01804910e-16 * tc[4];
        /*species 31: NH */
        species[31] =
            +4.21208480e+04 / tc[1]
            -2.95708710e+00
            -2.78369280e+00 * tc[0]
            -6.64921500e-04 * tc[1]
            +7.07967450e-08 * tc[2]
            -6.52904175e-12 * tc[3]
            +2.75222350e-16 * tc[4];
        /*species 32: NH2 */
        species[32] =
            +2.21719570e+04 / tc[1]
            -3.68567420e+00
            -2.83474210e+00 * tc[0]
            -1.60365410e-03 * tc[1]
            +1.55651340e-07 * tc[2]
            -1.14191275e-11 * tc[3]
            +3.96030720e-16 * tc[4];
        /*species 33: NH3 */
        species[33] =
            -6.54469580e+03 / tc[1]
            -3.93184070e+00
            -2.63445210e+00 * tc[0]
            -2.83312800e-03 * tc[1]
            +2.87977933e-07 * tc[2]
            -1.98893008e-11 * tc[3]
            +6.28939300e-16 * tc[4];
        /*species 34: NNH */
        species[34] =
            +2.86506970e+04 / tc[1]
            -7.03752300e-01
            -3.76675440e+00 * tc[0]
            -1.44575410e-03 * tc[1]
            +1.73610333e-07 * tc[2]
            -1.40354950e-11 * tc[3]
            +5.04594800e-16 * tc[4];
        /*species 35: NO */
        species[35] =
            +9.92097460e+03 / tc[1]
            -3.10869710e+00
            -3.26060560e+00 * tc[0]
            -5.95552150e-04 * tc[1]
            +7.15284133e-08 * tc[2]
            -5.78813908e-12 * tc[3]
            +2.01680495e-16 * tc[4];
        /*species 36: NO2 */
        species[36] =
            +2.31649830e+03 / tc[1]
            +5.00217115e+00
            -4.88475420e+00 * tc[0]
            -1.08619780e-03 * tc[1]
            +1.38011510e-07 * tc[2]
            -1.31229250e-11 * tc[3]
            +5.25544750e-16 * tc[4];
        /*species 37: N2O */
        species[37] =
            +8.07340480e+03 / tc[1]
            +7.02479360e+00
            -4.82307290e+00 * tc[0]
            -1.31351255e-03 * tc[1]
            +1.59751457e-07 * tc[2]
            -1.33339267e-11 * tc[3]
            +4.88761515e-16 * tc[4];
        /*species 38: HNO */
        species[38] =
            +1.17505820e+04 / tc[1]
            -5.62712190e+00
            -2.97925090e+00 * tc[0]
            -1.74720295e-03 * tc[1]
            +1.30916297e-07 * tc[2]
            -4.78996617e-12 * tc[3]
            +9.66795800e-18 * tc[4];
        /*species 39: CN */
        species[39] =
            +5.15361880e+04 / tc[1]
            +9.59220400e-01
            -3.74598050e+00 * tc[0]
            -2.17253875e-05 * tc[1]
            -4.95099733e-08 * tc[2]
            +5.72098383e-12 * tc[3]
            -2.20670865e-16 * tc[4];
        /*species 40: HCN */
        species[40] =
            +1.44072920e+04 / tc[1]
            +2.22677910e+00
            -3.80223920e+00 * tc[0]
            -1.57321140e-03 * tc[1]
            +1.77203083e-07 * tc[2]
            -1.38497975e-11 * tc[3]
            +4.89987850e-16 * tc[4];
        /*species 41: H2CN */
        species[41] =
            +2.76771090e+04 / tc[1]
            +9.65418100e+00
            -5.20970300e+00 * tc[0]
            -1.48464555e-03 * tc[1]
            +4.75931517e-08 * tc[2]
            +1.36295833e-11 * tc[3]
            -1.52162945e-15 * tc[4];
        /*species 42: HCNN */
        species[42] =
            +5.34529410e+04 / tc[1]
            +1.09976864e+01
            -5.89463620e+00 * tc[0]
            -1.99479795e-03 * tc[1]
            +2.66373000e-07 * tc[2]
            -2.43744958e-11 * tc[3]
            +1.00473430e-15 * tc[4];
        /*species 46: NCO */
        species[46] =
            +1.40041230e+04 / tc[1]
            +7.69645050e+00
            -5.15218450e+00 * tc[0]
            -1.15258805e-03 * tc[1]
            +1.46721922e-07 * tc[2]
            -1.23242483e-11 * tc[3]
            +4.54889980e-16 * tc[4];
        /*species 47: N2 */
        species[47] =
            -9.22797700e+02 / tc[1]
            -3.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 48: AR */
        species[48] =
            -7.45375000e+02 / tc[1]
            -1.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 49: C3H7 */
        species[49] =
            +8.29843360e+03 / tc[1]
            +2.31828787e+01
            -7.70269870e+00 * tc[0]
            -8.02210150e-03 * tc[1]
            +8.80553667e-07 * tc[2]
            -6.35821583e-11 * tc[3]
            +1.96961420e-15 * tc[4];
        /*species 50: C3H8 */
        species[50] =
            -1.64675160e+04 / tc[1]
            +2.54264858e+01
            -7.53413680e+00 * tc[0]
            -9.43611950e-03 * tc[1]
            +1.04530818e-06 * tc[2]
            -7.62297075e-11 * tc[3]
            +2.39190345e-15 * tc[4];
        /*species 51: CH2CHO */
        species[51] =
            +4.90321800e+02 / tc[1]
            +1.10209210e+01
            -5.97567000e+00 * tc[0]
            -4.06529550e-03 * tc[1]
            +4.57270667e-07 * tc[2]
            -3.39192000e-11 * tc[3]
            +1.08800850e-15 * tc[4];
        /*species 52: CH3CHO */
        species[52] =
            -2.25931220e+04 / tc[1]
            +8.88490250e+00
            -5.40411080e+00 * tc[0]
            -5.86152950e-03 * tc[1]
            +7.04385617e-07 * tc[2]
            -5.69770425e-11 * tc[3]
            +2.04924315e-15 * tc[4];
    }

    /*species with midpoint at T=1368 kelvin */
    if (T < 1368) {
        /*species 44: HOCN */
        species[44] =
            -2.82698400e+03 / tc[1]
            -1.84687210e+00
            -3.78604952e+00 * tc[0]
            -3.44333961e-03 * tc[1]
            +5.35813107e-07 * tc[2]
            -4.30996473e-11 * tc[3]
            -5.96803940e-16 * tc[4];
    } else {
        /*species 44: HOCN */
        species[44] =
            -3.70653331e+03 / tc[1]
            +1.20795271e+01
            -5.89784885e+00 * tc[0]
            -1.58394696e-03 * tc[1]
            +1.86335107e-07 * tc[2]
            -1.47702620e-11 * tc[3]
            +5.21695885e-16 * tc[4];
    }

    /*species with midpoint at T=1478 kelvin */
    if (T < 1478) {
        /*species 45: HNCO */
        species[45] =
            -1.55873636e+04 / tc[1]
            -2.56361410e+00
            -3.63096317e+00 * tc[0]
            -3.65141179e-03 * tc[1]
            +3.80083338e-07 * tc[2]
            +5.51059415e-11 * tc[3]
            -1.81117876e-14 * tc[4];
    } else {
        /*species 45: HNCO */
        species[45] =
            -1.66599344e+04 / tc[1]
            +1.46061988e+01
            -6.22395134e+00 * tc[0]
            -1.58932002e-03 * tc[1]
            +1.82297925e-07 * tc[2]
            -1.42279302e-11 * tc[3]
            +4.97510977e-16 * tc[4];
    }

    /*species with midpoint at T=1382 kelvin */
    if (T < 1382) {
        /*species 43: HCNO */
        species[43] =
            +1.92990252e+04 / tc[1]
            -8.08601731e+00
            -2.64727989e+00 * tc[0]
            -6.37526710e-03 * tc[1]
            +1.74657060e-06 * tc[2]
            -3.67860697e-10 * tc[3]
            +3.78760733e-14 * tc[4];
    } else {
        /*species 43: HCNO */
        species[43] =
            +1.79661339e+04 / tc[1]
            +1.69292645e+01
            -6.59860456e+00 * tc[0]
            -1.51389313e-03 * tc[1]
            +1.79507243e-07 * tc[2]
            -1.43055440e-11 * tc[3]
            +5.07196955e-16 * tc[4];
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
        /*species 30: N */
        species[30] =
            +5.61046370e+04 / tc[1]
            -2.69390870e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 31: NH */
        species[31] =
            +4.18806290e+04 / tc[1]
            +6.44580700e-01
            -3.49290850e+00 * tc[0]
            -1.55895990e-04 * tc[1]
            +2.48174733e-07 * tc[2]
            -2.06803683e-10 * tc[3]
            +5.17848350e-14 * tc[4];
        /*species 32: NH2 */
        species[32] =
            +2.18859100e+04 / tc[1]
            +3.34584538e+00
            -4.20400290e+00 * tc[0]
            +1.05306925e-03 * tc[1]
            -1.18447247e-06 * tc[2]
            +4.67626642e-10 * tc[3]
            -8.22035850e-14 * tc[4];
        /*species 33: NH3 */
        species[33] =
            -6.74172850e+03 / tc[1]
            +3.91140017e+00
            -4.28602740e+00 * tc[0]
            +2.33026150e-03 * tc[1]
            -3.61975217e-06 * tc[2]
            +1.90074058e-09 * tc[3]
            -4.13190230e-13 * tc[4];
        /*species 34: NNH */
        species[34] =
            +2.87919730e+04 / tc[1]
            +3.66751700e-01
            -4.34469270e+00 * tc[0]
            +2.42485360e-03 * tc[1]
            -3.34324317e-06 * tc[2]
            +1.81053867e-09 * tc[3]
            -3.97347695e-13 * tc[4];
        /*species 35: NO */
        species[35] =
            +9.84462300e+03 / tc[1]
            +9.37629900e-01
            -4.21847630e+00 * tc[0]
            +2.31948800e-03 * tc[1]
            -1.84017033e-06 * tc[2]
            +7.78011283e-10 * tc[3]
            -1.40178850e-13 * tc[4];
        /*species 36: NO2 */
        species[36] =
            +2.89661790e+03 / tc[1]
            -3.36796050e+00
            -3.94403120e+00 * tc[0]
            +7.92714500e-04 * tc[1]
            -2.77630200e-06 * tc[2]
            +1.70628550e-09 * tc[3]
            -3.91752820e-13 * tc[4];
        /*species 37: N2O */
        species[37] =
            +8.74177440e+03 / tc[1]
            -9.50084180e+00
            -2.25715020e+00 * tc[0]
            -5.65236400e-03 * tc[1]
            +2.27855317e-06 * tc[2]
            -8.06831717e-10 * tc[3]
            +1.46535910e-13 * tc[4];
        /*species 38: HNO */
        species[38] =
            +1.15482970e+04 / tc[1]
            +1.78364990e+00
            -4.53349160e+00 * tc[0]
            +2.83480855e-03 * tc[1]
            -3.07886783e-06 * tc[2]
            +1.42809117e-09 * tc[3]
            -2.77272865e-13 * tc[4];
        /*species 39: CN */
        species[39] =
            +5.17083400e+04 / tc[1]
            -1.36756440e+00
            -3.61293510e+00 * tc[0]
            +4.77756635e-04 * tc[1]
            -3.57382950e-07 * tc[2]
            +2.62636025e-11 * tc[3]
            +2.32151780e-14 * tc[4];
        /*species 40: HCN */
        species[40] =
            +1.47126330e+04 / tc[1]
            -7.65745330e+00
            -2.25898860e+00 * tc[0]
            -5.02558500e-03 * tc[1]
            +2.22529383e-06 * tc[2]
            -8.41029083e-10 * tc[3]
            +1.50445140e-13 * tc[4];
        /*species 41: H2CN */
        species[41] =
            +2.86378200e+04 / tc[1]
            -7.14109010e+00
            -2.85166100e+00 * tc[0]
            -2.84761655e-03 * tc[1]
            -1.78523333e-07 * tc[2]
            +1.35217667e-10 * tc[3]
            +1.17555405e-14 * tc[4];
        /*species 42: HCNN */
        species[42] =
            +5.42619840e+04 / tc[1]
            -1.01515506e+01
            -2.52431940e+00 * tc[0]
            -7.98030950e-03 * tc[1]
            +3.13605900e-06 * tc[2]
            -1.01046167e-09 * tc[3]
            +1.61786890e-13 * tc[4];
        /*species 46: NCO */
        species[46] =
            +1.46824770e+04 / tc[1]
            -7.72353380e+00
            -2.82693080e+00 * tc[0]
            -4.40258440e-03 * tc[1]
            +1.39776890e-06 * tc[2]
            -4.00141367e-10 * tc[3]
            +6.65679750e-14 * tc[4];
        /*species 47: N2 */
        species[47] =
            -1.02089990e+03 / tc[1]
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
        /*species 48: AR */
        species[48] =
            -7.45375000e+02 / tc[1]
            -2.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 49: C3H7 */
        species[49] =
            +1.06318630e+04 / tc[1]
            -2.10710072e+01
            -1.05155180e+00 * tc[0]
            -1.29959900e-02 * tc[1]
            -3.96675667e-07 * tc[2]
            +1.63413075e-09 * tc[3]
            -4.68662350e-13 * tc[4];
        /*species 50: C3H8 */
        species[50] =
            -1.39585200e+04 / tc[1]
            -1.92681372e+01
            -9.33553810e-01 * tc[0]
            -1.32122895e-02 * tc[1]
            -1.01766212e-06 * tc[2]
            +1.83145825e-09 * tc[3]
            -4.75746265e-13 * tc[4];
        /*species 51: CH2CHO */
        species[51] =
            +1.52147660e+03 / tc[1]
            -7.14922800e+00
            -3.40906200e+00 * tc[0]
            -5.36928700e-03 * tc[1]
            -3.15248667e-07 * tc[2]
            +5.96548583e-10 * tc[3]
            -1.43369250e-13 * tc[4];
        /*species 52: CH3CHO */
        species[52] =
            -2.15728780e+04 / tc[1]
            -3.73556400e-01
            -4.72945950e+00 * tc[0]
            +1.59664290e-03 * tc[1]
            -7.92248683e-06 * tc[2]
            +4.78821758e-09 * tc[3]
            -1.09655560e-12 * tc[4];
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
        /*species 30: N */
        species[30] =
            +5.61337730e+04 / tc[1]
            -3.23366670e+00
            -2.41594290e+00 * tc[0]
            -8.74453250e-05 * tc[1]
            +1.98372817e-08 * tc[2]
            -2.51885375e-12 * tc[3]
            +1.01804910e-16 * tc[4];
        /*species 31: NH */
        species[31] =
            +4.21208480e+04 / tc[1]
            -3.95708710e+00
            -2.78369280e+00 * tc[0]
            -6.64921500e-04 * tc[1]
            +7.07967450e-08 * tc[2]
            -6.52904175e-12 * tc[3]
            +2.75222350e-16 * tc[4];
        /*species 32: NH2 */
        species[32] =
            +2.21719570e+04 / tc[1]
            -4.68567420e+00
            -2.83474210e+00 * tc[0]
            -1.60365410e-03 * tc[1]
            +1.55651340e-07 * tc[2]
            -1.14191275e-11 * tc[3]
            +3.96030720e-16 * tc[4];
        /*species 33: NH3 */
        species[33] =
            -6.54469580e+03 / tc[1]
            -4.93184070e+00
            -2.63445210e+00 * tc[0]
            -2.83312800e-03 * tc[1]
            +2.87977933e-07 * tc[2]
            -1.98893008e-11 * tc[3]
            +6.28939300e-16 * tc[4];
        /*species 34: NNH */
        species[34] =
            +2.86506970e+04 / tc[1]
            -1.70375230e+00
            -3.76675440e+00 * tc[0]
            -1.44575410e-03 * tc[1]
            +1.73610333e-07 * tc[2]
            -1.40354950e-11 * tc[3]
            +5.04594800e-16 * tc[4];
        /*species 35: NO */
        species[35] =
            +9.92097460e+03 / tc[1]
            -4.10869710e+00
            -3.26060560e+00 * tc[0]
            -5.95552150e-04 * tc[1]
            +7.15284133e-08 * tc[2]
            -5.78813908e-12 * tc[3]
            +2.01680495e-16 * tc[4];
        /*species 36: NO2 */
        species[36] =
            +2.31649830e+03 / tc[1]
            +4.00217115e+00
            -4.88475420e+00 * tc[0]
            -1.08619780e-03 * tc[1]
            +1.38011510e-07 * tc[2]
            -1.31229250e-11 * tc[3]
            +5.25544750e-16 * tc[4];
        /*species 37: N2O */
        species[37] =
            +8.07340480e+03 / tc[1]
            +6.02479360e+00
            -4.82307290e+00 * tc[0]
            -1.31351255e-03 * tc[1]
            +1.59751457e-07 * tc[2]
            -1.33339267e-11 * tc[3]
            +4.88761515e-16 * tc[4];
        /*species 38: HNO */
        species[38] =
            +1.17505820e+04 / tc[1]
            -6.62712190e+00
            -2.97925090e+00 * tc[0]
            -1.74720295e-03 * tc[1]
            +1.30916297e-07 * tc[2]
            -4.78996617e-12 * tc[3]
            +9.66795800e-18 * tc[4];
        /*species 39: CN */
        species[39] =
            +5.15361880e+04 / tc[1]
            -4.07796000e-02
            -3.74598050e+00 * tc[0]
            -2.17253875e-05 * tc[1]
            -4.95099733e-08 * tc[2]
            +5.72098383e-12 * tc[3]
            -2.20670865e-16 * tc[4];
        /*species 40: HCN */
        species[40] =
            +1.44072920e+04 / tc[1]
            +1.22677910e+00
            -3.80223920e+00 * tc[0]
            -1.57321140e-03 * tc[1]
            +1.77203083e-07 * tc[2]
            -1.38497975e-11 * tc[3]
            +4.89987850e-16 * tc[4];
        /*species 41: H2CN */
        species[41] =
            +2.76771090e+04 / tc[1]
            +8.65418100e+00
            -5.20970300e+00 * tc[0]
            -1.48464555e-03 * tc[1]
            +4.75931517e-08 * tc[2]
            +1.36295833e-11 * tc[3]
            -1.52162945e-15 * tc[4];
        /*species 42: HCNN */
        species[42] =
            +5.34529410e+04 / tc[1]
            +9.99768640e+00
            -5.89463620e+00 * tc[0]
            -1.99479795e-03 * tc[1]
            +2.66373000e-07 * tc[2]
            -2.43744958e-11 * tc[3]
            +1.00473430e-15 * tc[4];
        /*species 46: NCO */
        species[46] =
            +1.40041230e+04 / tc[1]
            +6.69645050e+00
            -5.15218450e+00 * tc[0]
            -1.15258805e-03 * tc[1]
            +1.46721922e-07 * tc[2]
            -1.23242483e-11 * tc[3]
            +4.54889980e-16 * tc[4];
        /*species 47: N2 */
        species[47] =
            -9.22797700e+02 / tc[1]
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 48: AR */
        species[48] =
            -7.45375000e+02 / tc[1]
            -2.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 49: C3H7 */
        species[49] =
            +8.29843360e+03 / tc[1]
            +2.21828787e+01
            -7.70269870e+00 * tc[0]
            -8.02210150e-03 * tc[1]
            +8.80553667e-07 * tc[2]
            -6.35821583e-11 * tc[3]
            +1.96961420e-15 * tc[4];
        /*species 50: C3H8 */
        species[50] =
            -1.64675160e+04 / tc[1]
            +2.44264858e+01
            -7.53413680e+00 * tc[0]
            -9.43611950e-03 * tc[1]
            +1.04530818e-06 * tc[2]
            -7.62297075e-11 * tc[3]
            +2.39190345e-15 * tc[4];
        /*species 51: CH2CHO */
        species[51] =
            +4.90321800e+02 / tc[1]
            +1.00209210e+01
            -5.97567000e+00 * tc[0]
            -4.06529550e-03 * tc[1]
            +4.57270667e-07 * tc[2]
            -3.39192000e-11 * tc[3]
            +1.08800850e-15 * tc[4];
        /*species 52: CH3CHO */
        species[52] =
            -2.25931220e+04 / tc[1]
            +7.88490250e+00
            -5.40411080e+00 * tc[0]
            -5.86152950e-03 * tc[1]
            +7.04385617e-07 * tc[2]
            -5.69770425e-11 * tc[3]
            +2.04924315e-15 * tc[4];
    }

    /*species with midpoint at T=1368 kelvin */
    if (T < 1368) {
        /*species 44: HOCN */
        species[44] =
            -2.82698400e+03 / tc[1]
            -2.84687210e+00
            -3.78604952e+00 * tc[0]
            -3.44333961e-03 * tc[1]
            +5.35813107e-07 * tc[2]
            -4.30996473e-11 * tc[3]
            -5.96803940e-16 * tc[4];
    } else {
        /*species 44: HOCN */
        species[44] =
            -3.70653331e+03 / tc[1]
            +1.10795271e+01
            -5.89784885e+00 * tc[0]
            -1.58394696e-03 * tc[1]
            +1.86335107e-07 * tc[2]
            -1.47702620e-11 * tc[3]
            +5.21695885e-16 * tc[4];
    }

    /*species with midpoint at T=1478 kelvin */
    if (T < 1478) {
        /*species 45: HNCO */
        species[45] =
            -1.55873636e+04 / tc[1]
            -3.56361410e+00
            -3.63096317e+00 * tc[0]
            -3.65141179e-03 * tc[1]
            +3.80083338e-07 * tc[2]
            +5.51059415e-11 * tc[3]
            -1.81117876e-14 * tc[4];
    } else {
        /*species 45: HNCO */
        species[45] =
            -1.66599344e+04 / tc[1]
            +1.36061988e+01
            -6.22395134e+00 * tc[0]
            -1.58932002e-03 * tc[1]
            +1.82297925e-07 * tc[2]
            -1.42279302e-11 * tc[3]
            +4.97510977e-16 * tc[4];
    }

    /*species with midpoint at T=1382 kelvin */
    if (T < 1382) {
        /*species 43: HCNO */
        species[43] =
            +1.92990252e+04 / tc[1]
            -9.08601731e+00
            -2.64727989e+00 * tc[0]
            -6.37526710e-03 * tc[1]
            +1.74657060e-06 * tc[2]
            -3.67860697e-10 * tc[3]
            +3.78760733e-14 * tc[4];
    } else {
        /*species 43: HCNO */
        species[43] =
            +1.79661339e+04 / tc[1]
            +1.59292645e+01
            -6.59860456e+00 * tc[0]
            -1.51389313e-03 * tc[1]
            +1.79507243e-07 * tc[2]
            -1.43055440e-11 * tc[3]
            +5.07196955e-16 * tc[4];
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
        /*species 30: N */
        species[30] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 31: NH */
        species[31] =
            +2.49290850e+00
            +3.11791980e-04 * tc[1]
            -1.48904840e-06 * tc[2]
            +2.48164420e-09 * tc[3]
            -1.03569670e-12 * tc[4];
        /*species 32: NH2 */
        species[32] =
            +3.20400290e+00
            -2.10613850e-03 * tc[1]
            +7.10683480e-06 * tc[2]
            -5.61151970e-09 * tc[3]
            +1.64407170e-12 * tc[4];
        /*species 33: NH3 */
        species[33] =
            +3.28602740e+00
            -4.66052300e-03 * tc[1]
            +2.17185130e-05 * tc[2]
            -2.28088870e-08 * tc[3]
            +8.26380460e-12 * tc[4];
        /*species 34: NNH */
        species[34] =
            +3.34469270e+00
            -4.84970720e-03 * tc[1]
            +2.00594590e-05 * tc[2]
            -2.17264640e-08 * tc[3]
            +7.94695390e-12 * tc[4];
        /*species 35: NO */
        species[35] =
            +3.21847630e+00
            -4.63897600e-03 * tc[1]
            +1.10410220e-05 * tc[2]
            -9.33613540e-09 * tc[3]
            +2.80357700e-12 * tc[4];
        /*species 36: NO2 */
        species[36] =
            +2.94403120e+00
            -1.58542900e-03 * tc[1]
            +1.66578120e-05 * tc[2]
            -2.04754260e-08 * tc[3]
            +7.83505640e-12 * tc[4];
        /*species 37: N2O */
        species[37] =
            +1.25715020e+00
            +1.13047280e-02 * tc[1]
            -1.36713190e-05 * tc[2]
            +9.68198060e-09 * tc[3]
            -2.93071820e-12 * tc[4];
        /*species 38: HNO */
        species[38] =
            +3.53349160e+00
            -5.66961710e-03 * tc[1]
            +1.84732070e-05 * tc[2]
            -1.71370940e-08 * tc[3]
            +5.54545730e-12 * tc[4];
        /*species 39: CN */
        species[39] =
            +2.61293510e+00
            -9.55513270e-04 * tc[1]
            +2.14429770e-06 * tc[2]
            -3.15163230e-10 * tc[3]
            -4.64303560e-13 * tc[4];
        /*species 40: HCN */
        species[40] =
            +1.25898860e+00
            +1.00511700e-02 * tc[1]
            -1.33517630e-05 * tc[2]
            +1.00923490e-08 * tc[3]
            -3.00890280e-12 * tc[4];
        /*species 41: H2CN */
        species[41] =
            +1.85166100e+00
            +5.69523310e-03 * tc[1]
            +1.07114000e-06 * tc[2]
            -1.62261200e-09 * tc[3]
            -2.35110810e-13 * tc[4];
        /*species 42: HCNN */
        species[42] =
            +1.52431940e+00
            +1.59606190e-02 * tc[1]
            -1.88163540e-05 * tc[2]
            +1.21255400e-08 * tc[3]
            -3.23573780e-12 * tc[4];
        /*species 46: NCO */
        species[46] =
            +1.82693080e+00
            +8.80516880e-03 * tc[1]
            -8.38661340e-06 * tc[2]
            +4.80169640e-09 * tc[3]
            -1.33135950e-12 * tc[4];
        /*species 47: N2 */
        species[47] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 48: AR */
        species[48] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 49: C3H7 */
        species[49] =
            +5.15518000e-02
            +2.59919800e-02 * tc[1]
            +2.38005400e-06 * tc[2]
            -1.96095690e-08 * tc[3]
            +9.37324700e-12 * tc[4];
        /*species 50: C3H8 */
        species[50] =
            -6.64461900e-02
            +2.64245790e-02 * tc[1]
            +6.10597270e-06 * tc[2]
            -2.19774990e-08 * tc[3]
            +9.51492530e-12 * tc[4];
        /*species 51: CH2CHO */
        species[51] =
            +2.40906200e+00
            +1.07385740e-02 * tc[1]
            +1.89149200e-06 * tc[2]
            -7.15858300e-09 * tc[3]
            +2.86738500e-12 * tc[4];
        /*species 52: CH3CHO */
        species[52] =
            +3.72945950e+00
            -3.19328580e-03 * tc[1]
            +4.75349210e-05 * tc[2]
            -5.74586110e-08 * tc[3]
            +2.19311120e-11 * tc[4];
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
        /*species 30: N */
        species[30] =
            +1.41594290e+00
            +1.74890650e-04 * tc[1]
            -1.19023690e-07 * tc[2]
            +3.02262450e-11 * tc[3]
            -2.03609820e-15 * tc[4];
        /*species 31: NH */
        species[31] =
            +1.78369280e+00
            +1.32984300e-03 * tc[1]
            -4.24780470e-07 * tc[2]
            +7.83485010e-11 * tc[3]
            -5.50444700e-15 * tc[4];
        /*species 32: NH2 */
        species[32] =
            +1.83474210e+00
            +3.20730820e-03 * tc[1]
            -9.33908040e-07 * tc[2]
            +1.37029530e-10 * tc[3]
            -7.92061440e-15 * tc[4];
        /*species 33: NH3 */
        species[33] =
            +1.63445210e+00
            +5.66625600e-03 * tc[1]
            -1.72786760e-06 * tc[2]
            +2.38671610e-10 * tc[3]
            -1.25787860e-14 * tc[4];
        /*species 34: NNH */
        species[34] =
            +2.76675440e+00
            +2.89150820e-03 * tc[1]
            -1.04166200e-06 * tc[2]
            +1.68425940e-10 * tc[3]
            -1.00918960e-14 * tc[4];
        /*species 35: NO */
        species[35] =
            +2.26060560e+00
            +1.19110430e-03 * tc[1]
            -4.29170480e-07 * tc[2]
            +6.94576690e-11 * tc[3]
            -4.03360990e-15 * tc[4];
        /*species 36: NO2 */
        species[36] =
            +3.88475420e+00
            +2.17239560e-03 * tc[1]
            -8.28069060e-07 * tc[2]
            +1.57475100e-10 * tc[3]
            -1.05108950e-14 * tc[4];
        /*species 37: N2O */
        species[37] =
            +3.82307290e+00
            +2.62702510e-03 * tc[1]
            -9.58508740e-07 * tc[2]
            +1.60007120e-10 * tc[3]
            -9.77523030e-15 * tc[4];
        /*species 38: HNO */
        species[38] =
            +1.97925090e+00
            +3.49440590e-03 * tc[1]
            -7.85497780e-07 * tc[2]
            +5.74795940e-11 * tc[3]
            -1.93359160e-16 * tc[4];
        /*species 39: CN */
        species[39] =
            +2.74598050e+00
            +4.34507750e-05 * tc[1]
            +2.97059840e-07 * tc[2]
            -6.86518060e-11 * tc[3]
            +4.41341730e-15 * tc[4];
        /*species 40: HCN */
        species[40] =
            +2.80223920e+00
            +3.14642280e-03 * tc[1]
            -1.06321850e-06 * tc[2]
            +1.66197570e-10 * tc[3]
            -9.79975700e-15 * tc[4];
        /*species 41: H2CN */
        species[41] =
            +4.20970300e+00
            +2.96929110e-03 * tc[1]
            -2.85558910e-07 * tc[2]
            -1.63555000e-10 * tc[3]
            +3.04325890e-14 * tc[4];
        /*species 42: HCNN */
        species[42] =
            +4.89463620e+00
            +3.98959590e-03 * tc[1]
            -1.59823800e-06 * tc[2]
            +2.92493950e-10 * tc[3]
            -2.00946860e-14 * tc[4];
        /*species 46: NCO */
        species[46] =
            +4.15218450e+00
            +2.30517610e-03 * tc[1]
            -8.80331530e-07 * tc[2]
            +1.47890980e-10 * tc[3]
            -9.09779960e-15 * tc[4];
        /*species 47: N2 */
        species[47] =
            +1.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 48: AR */
        species[48] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 49: C3H7 */
        species[49] =
            +6.70269870e+00
            +1.60442030e-02 * tc[1]
            -5.28332200e-06 * tc[2]
            +7.62985900e-10 * tc[3]
            -3.93922840e-14 * tc[4];
        /*species 50: C3H8 */
        species[50] =
            +6.53413680e+00
            +1.88722390e-02 * tc[1]
            -6.27184910e-06 * tc[2]
            +9.14756490e-10 * tc[3]
            -4.78380690e-14 * tc[4];
        /*species 51: CH2CHO */
        species[51] =
            +4.97567000e+00
            +8.13059100e-03 * tc[1]
            -2.74362400e-06 * tc[2]
            +4.07030400e-10 * tc[3]
            -2.17601700e-14 * tc[4];
        /*species 52: CH3CHO */
        species[52] =
            +4.40411080e+00
            +1.17230590e-02 * tc[1]
            -4.22631370e-06 * tc[2]
            +6.83724510e-10 * tc[3]
            -4.09848630e-14 * tc[4];
    }

    /*species with midpoint at T=1368 kelvin */
    if (T < 1368) {
        /*species 44: HOCN */
        species[44] =
            +2.78604952e+00
            +6.88667922e-03 * tc[1]
            -3.21487864e-06 * tc[2]
            +5.17195767e-10 * tc[3]
            +1.19360788e-14 * tc[4];
    } else {
        /*species 44: HOCN */
        species[44] =
            +4.89784885e+00
            +3.16789393e-03 * tc[1]
            -1.11801064e-06 * tc[2]
            +1.77243144e-10 * tc[3]
            -1.04339177e-14 * tc[4];
    }

    /*species with midpoint at T=1478 kelvin */
    if (T < 1478) {
        /*species 45: HNCO */
        species[45] =
            +2.63096317e+00
            +7.30282357e-03 * tc[1]
            -2.28050003e-06 * tc[2]
            -6.61271298e-10 * tc[3]
            +3.62235752e-13 * tc[4];
    } else {
        /*species 45: HNCO */
        species[45] =
            +5.22395134e+00
            +3.17864004e-03 * tc[1]
            -1.09378755e-06 * tc[2]
            +1.70735163e-10 * tc[3]
            -9.95021955e-15 * tc[4];
    }

    /*species with midpoint at T=1382 kelvin */
    if (T < 1382) {
        /*species 43: HCNO */
        species[43] =
            +1.64727989e+00
            +1.27505342e-02 * tc[1]
            -1.04794236e-05 * tc[2]
            +4.41432836e-09 * tc[3]
            -7.57521466e-13 * tc[4];
    } else {
        /*species 43: HCNO */
        species[43] =
            +5.59860456e+00
            +3.02778626e-03 * tc[1]
            -1.07704346e-06 * tc[2]
            +1.71666528e-10 * tc[3]
            -1.01439391e-14 * tc[4];
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
        /*species 30: N */
        species[30] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 31: NH */
        species[31] =
            +3.49290850e+00
            +3.11791980e-04 * tc[1]
            -1.48904840e-06 * tc[2]
            +2.48164420e-09 * tc[3]
            -1.03569670e-12 * tc[4];
        /*species 32: NH2 */
        species[32] =
            +4.20400290e+00
            -2.10613850e-03 * tc[1]
            +7.10683480e-06 * tc[2]
            -5.61151970e-09 * tc[3]
            +1.64407170e-12 * tc[4];
        /*species 33: NH3 */
        species[33] =
            +4.28602740e+00
            -4.66052300e-03 * tc[1]
            +2.17185130e-05 * tc[2]
            -2.28088870e-08 * tc[3]
            +8.26380460e-12 * tc[4];
        /*species 34: NNH */
        species[34] =
            +4.34469270e+00
            -4.84970720e-03 * tc[1]
            +2.00594590e-05 * tc[2]
            -2.17264640e-08 * tc[3]
            +7.94695390e-12 * tc[4];
        /*species 35: NO */
        species[35] =
            +4.21847630e+00
            -4.63897600e-03 * tc[1]
            +1.10410220e-05 * tc[2]
            -9.33613540e-09 * tc[3]
            +2.80357700e-12 * tc[4];
        /*species 36: NO2 */
        species[36] =
            +3.94403120e+00
            -1.58542900e-03 * tc[1]
            +1.66578120e-05 * tc[2]
            -2.04754260e-08 * tc[3]
            +7.83505640e-12 * tc[4];
        /*species 37: N2O */
        species[37] =
            +2.25715020e+00
            +1.13047280e-02 * tc[1]
            -1.36713190e-05 * tc[2]
            +9.68198060e-09 * tc[3]
            -2.93071820e-12 * tc[4];
        /*species 38: HNO */
        species[38] =
            +4.53349160e+00
            -5.66961710e-03 * tc[1]
            +1.84732070e-05 * tc[2]
            -1.71370940e-08 * tc[3]
            +5.54545730e-12 * tc[4];
        /*species 39: CN */
        species[39] =
            +3.61293510e+00
            -9.55513270e-04 * tc[1]
            +2.14429770e-06 * tc[2]
            -3.15163230e-10 * tc[3]
            -4.64303560e-13 * tc[4];
        /*species 40: HCN */
        species[40] =
            +2.25898860e+00
            +1.00511700e-02 * tc[1]
            -1.33517630e-05 * tc[2]
            +1.00923490e-08 * tc[3]
            -3.00890280e-12 * tc[4];
        /*species 41: H2CN */
        species[41] =
            +2.85166100e+00
            +5.69523310e-03 * tc[1]
            +1.07114000e-06 * tc[2]
            -1.62261200e-09 * tc[3]
            -2.35110810e-13 * tc[4];
        /*species 42: HCNN */
        species[42] =
            +2.52431940e+00
            +1.59606190e-02 * tc[1]
            -1.88163540e-05 * tc[2]
            +1.21255400e-08 * tc[3]
            -3.23573780e-12 * tc[4];
        /*species 46: NCO */
        species[46] =
            +2.82693080e+00
            +8.80516880e-03 * tc[1]
            -8.38661340e-06 * tc[2]
            +4.80169640e-09 * tc[3]
            -1.33135950e-12 * tc[4];
        /*species 47: N2 */
        species[47] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 48: AR */
        species[48] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 49: C3H7 */
        species[49] =
            +1.05155180e+00
            +2.59919800e-02 * tc[1]
            +2.38005400e-06 * tc[2]
            -1.96095690e-08 * tc[3]
            +9.37324700e-12 * tc[4];
        /*species 50: C3H8 */
        species[50] =
            +9.33553810e-01
            +2.64245790e-02 * tc[1]
            +6.10597270e-06 * tc[2]
            -2.19774990e-08 * tc[3]
            +9.51492530e-12 * tc[4];
        /*species 51: CH2CHO */
        species[51] =
            +3.40906200e+00
            +1.07385740e-02 * tc[1]
            +1.89149200e-06 * tc[2]
            -7.15858300e-09 * tc[3]
            +2.86738500e-12 * tc[4];
        /*species 52: CH3CHO */
        species[52] =
            +4.72945950e+00
            -3.19328580e-03 * tc[1]
            +4.75349210e-05 * tc[2]
            -5.74586110e-08 * tc[3]
            +2.19311120e-11 * tc[4];
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
        /*species 30: N */
        species[30] =
            +2.41594290e+00
            +1.74890650e-04 * tc[1]
            -1.19023690e-07 * tc[2]
            +3.02262450e-11 * tc[3]
            -2.03609820e-15 * tc[4];
        /*species 31: NH */
        species[31] =
            +2.78369280e+00
            +1.32984300e-03 * tc[1]
            -4.24780470e-07 * tc[2]
            +7.83485010e-11 * tc[3]
            -5.50444700e-15 * tc[4];
        /*species 32: NH2 */
        species[32] =
            +2.83474210e+00
            +3.20730820e-03 * tc[1]
            -9.33908040e-07 * tc[2]
            +1.37029530e-10 * tc[3]
            -7.92061440e-15 * tc[4];
        /*species 33: NH3 */
        species[33] =
            +2.63445210e+00
            +5.66625600e-03 * tc[1]
            -1.72786760e-06 * tc[2]
            +2.38671610e-10 * tc[3]
            -1.25787860e-14 * tc[4];
        /*species 34: NNH */
        species[34] =
            +3.76675440e+00
            +2.89150820e-03 * tc[1]
            -1.04166200e-06 * tc[2]
            +1.68425940e-10 * tc[3]
            -1.00918960e-14 * tc[4];
        /*species 35: NO */
        species[35] =
            +3.26060560e+00
            +1.19110430e-03 * tc[1]
            -4.29170480e-07 * tc[2]
            +6.94576690e-11 * tc[3]
            -4.03360990e-15 * tc[4];
        /*species 36: NO2 */
        species[36] =
            +4.88475420e+00
            +2.17239560e-03 * tc[1]
            -8.28069060e-07 * tc[2]
            +1.57475100e-10 * tc[3]
            -1.05108950e-14 * tc[4];
        /*species 37: N2O */
        species[37] =
            +4.82307290e+00
            +2.62702510e-03 * tc[1]
            -9.58508740e-07 * tc[2]
            +1.60007120e-10 * tc[3]
            -9.77523030e-15 * tc[4];
        /*species 38: HNO */
        species[38] =
            +2.97925090e+00
            +3.49440590e-03 * tc[1]
            -7.85497780e-07 * tc[2]
            +5.74795940e-11 * tc[3]
            -1.93359160e-16 * tc[4];
        /*species 39: CN */
        species[39] =
            +3.74598050e+00
            +4.34507750e-05 * tc[1]
            +2.97059840e-07 * tc[2]
            -6.86518060e-11 * tc[3]
            +4.41341730e-15 * tc[4];
        /*species 40: HCN */
        species[40] =
            +3.80223920e+00
            +3.14642280e-03 * tc[1]
            -1.06321850e-06 * tc[2]
            +1.66197570e-10 * tc[3]
            -9.79975700e-15 * tc[4];
        /*species 41: H2CN */
        species[41] =
            +5.20970300e+00
            +2.96929110e-03 * tc[1]
            -2.85558910e-07 * tc[2]
            -1.63555000e-10 * tc[3]
            +3.04325890e-14 * tc[4];
        /*species 42: HCNN */
        species[42] =
            +5.89463620e+00
            +3.98959590e-03 * tc[1]
            -1.59823800e-06 * tc[2]
            +2.92493950e-10 * tc[3]
            -2.00946860e-14 * tc[4];
        /*species 46: NCO */
        species[46] =
            +5.15218450e+00
            +2.30517610e-03 * tc[1]
            -8.80331530e-07 * tc[2]
            +1.47890980e-10 * tc[3]
            -9.09779960e-15 * tc[4];
        /*species 47: N2 */
        species[47] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 48: AR */
        species[48] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 49: C3H7 */
        species[49] =
            +7.70269870e+00
            +1.60442030e-02 * tc[1]
            -5.28332200e-06 * tc[2]
            +7.62985900e-10 * tc[3]
            -3.93922840e-14 * tc[4];
        /*species 50: C3H8 */
        species[50] =
            +7.53413680e+00
            +1.88722390e-02 * tc[1]
            -6.27184910e-06 * tc[2]
            +9.14756490e-10 * tc[3]
            -4.78380690e-14 * tc[4];
        /*species 51: CH2CHO */
        species[51] =
            +5.97567000e+00
            +8.13059100e-03 * tc[1]
            -2.74362400e-06 * tc[2]
            +4.07030400e-10 * tc[3]
            -2.17601700e-14 * tc[4];
        /*species 52: CH3CHO */
        species[52] =
            +5.40411080e+00
            +1.17230590e-02 * tc[1]
            -4.22631370e-06 * tc[2]
            +6.83724510e-10 * tc[3]
            -4.09848630e-14 * tc[4];
    }

    /*species with midpoint at T=1368 kelvin */
    if (T < 1368) {
        /*species 44: HOCN */
        species[44] =
            +3.78604952e+00
            +6.88667922e-03 * tc[1]
            -3.21487864e-06 * tc[2]
            +5.17195767e-10 * tc[3]
            +1.19360788e-14 * tc[4];
    } else {
        /*species 44: HOCN */
        species[44] =
            +5.89784885e+00
            +3.16789393e-03 * tc[1]
            -1.11801064e-06 * tc[2]
            +1.77243144e-10 * tc[3]
            -1.04339177e-14 * tc[4];
    }

    /*species with midpoint at T=1478 kelvin */
    if (T < 1478) {
        /*species 45: HNCO */
        species[45] =
            +3.63096317e+00
            +7.30282357e-03 * tc[1]
            -2.28050003e-06 * tc[2]
            -6.61271298e-10 * tc[3]
            +3.62235752e-13 * tc[4];
    } else {
        /*species 45: HNCO */
        species[45] =
            +6.22395134e+00
            +3.17864004e-03 * tc[1]
            -1.09378755e-06 * tc[2]
            +1.70735163e-10 * tc[3]
            -9.95021955e-15 * tc[4];
    }

    /*species with midpoint at T=1382 kelvin */
    if (T < 1382) {
        /*species 43: HCNO */
        species[43] =
            +2.64727989e+00
            +1.27505342e-02 * tc[1]
            -1.04794236e-05 * tc[2]
            +4.41432836e-09 * tc[3]
            -7.57521466e-13 * tc[4];
    } else {
        /*species 43: HCNO */
        species[43] =
            +6.59860456e+00
            +3.02778626e-03 * tc[1]
            -1.07704346e-06 * tc[2]
            +1.71666528e-10 * tc[3]
            -1.01439391e-14 * tc[4];
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
        /*species 30: N */
        species[30] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +5.61046370e+04 / tc[1];
        /*species 31: NH */
        species[31] =
            +2.49290850e+00
            +1.55895990e-04 * tc[1]
            -4.96349467e-07 * tc[2]
            +6.20411050e-10 * tc[3]
            -2.07139340e-13 * tc[4]
            +4.18806290e+04 / tc[1];
        /*species 32: NH2 */
        species[32] =
            +3.20400290e+00
            -1.05306925e-03 * tc[1]
            +2.36894493e-06 * tc[2]
            -1.40287992e-09 * tc[3]
            +3.28814340e-13 * tc[4]
            +2.18859100e+04 / tc[1];
        /*species 33: NH3 */
        species[33] =
            +3.28602740e+00
            -2.33026150e-03 * tc[1]
            +7.23950433e-06 * tc[2]
            -5.70222175e-09 * tc[3]
            +1.65276092e-12 * tc[4]
            -6.74172850e+03 / tc[1];
        /*species 34: NNH */
        species[34] =
            +3.34469270e+00
            -2.42485360e-03 * tc[1]
            +6.68648633e-06 * tc[2]
            -5.43161600e-09 * tc[3]
            +1.58939078e-12 * tc[4]
            +2.87919730e+04 / tc[1];
        /*species 35: NO */
        species[35] =
            +3.21847630e+00
            -2.31948800e-03 * tc[1]
            +3.68034067e-06 * tc[2]
            -2.33403385e-09 * tc[3]
            +5.60715400e-13 * tc[4]
            +9.84462300e+03 / tc[1];
        /*species 36: NO2 */
        species[36] =
            +2.94403120e+00
            -7.92714500e-04 * tc[1]
            +5.55260400e-06 * tc[2]
            -5.11885650e-09 * tc[3]
            +1.56701128e-12 * tc[4]
            +2.89661790e+03 / tc[1];
        /*species 37: N2O */
        species[37] =
            +1.25715020e+00
            +5.65236400e-03 * tc[1]
            -4.55710633e-06 * tc[2]
            +2.42049515e-09 * tc[3]
            -5.86143640e-13 * tc[4]
            +8.74177440e+03 / tc[1];
        /*species 38: HNO */
        species[38] =
            +3.53349160e+00
            -2.83480855e-03 * tc[1]
            +6.15773567e-06 * tc[2]
            -4.28427350e-09 * tc[3]
            +1.10909146e-12 * tc[4]
            +1.15482970e+04 / tc[1];
        /*species 39: CN */
        species[39] =
            +2.61293510e+00
            -4.77756635e-04 * tc[1]
            +7.14765900e-07 * tc[2]
            -7.87908075e-11 * tc[3]
            -9.28607120e-14 * tc[4]
            +5.17083400e+04 / tc[1];
        /*species 40: HCN */
        species[40] =
            +1.25898860e+00
            +5.02558500e-03 * tc[1]
            -4.45058767e-06 * tc[2]
            +2.52308725e-09 * tc[3]
            -6.01780560e-13 * tc[4]
            +1.47126330e+04 / tc[1];
        /*species 41: H2CN */
        species[41] =
            +1.85166100e+00
            +2.84761655e-03 * tc[1]
            +3.57046667e-07 * tc[2]
            -4.05653000e-10 * tc[3]
            -4.70221620e-14 * tc[4]
            +2.86378200e+04 / tc[1];
        /*species 42: HCNN */
        species[42] =
            +1.52431940e+00
            +7.98030950e-03 * tc[1]
            -6.27211800e-06 * tc[2]
            +3.03138500e-09 * tc[3]
            -6.47147560e-13 * tc[4]
            +5.42619840e+04 / tc[1];
        /*species 46: NCO */
        species[46] =
            +1.82693080e+00
            +4.40258440e-03 * tc[1]
            -2.79553780e-06 * tc[2]
            +1.20042410e-09 * tc[3]
            -2.66271900e-13 * tc[4]
            +1.46824770e+04 / tc[1];
        /*species 47: N2 */
        species[47] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 / tc[1];
        /*species 48: AR */
        species[48] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 / tc[1];
        /*species 49: C3H7 */
        species[49] =
            +5.15518000e-02
            +1.29959900e-02 * tc[1]
            +7.93351333e-07 * tc[2]
            -4.90239225e-09 * tc[3]
            +1.87464940e-12 * tc[4]
            +1.06318630e+04 / tc[1];
        /*species 50: C3H8 */
        species[50] =
            -6.64461900e-02
            +1.32122895e-02 * tc[1]
            +2.03532423e-06 * tc[2]
            -5.49437475e-09 * tc[3]
            +1.90298506e-12 * tc[4]
            -1.39585200e+04 / tc[1];
        /*species 51: CH2CHO */
        species[51] =
            +2.40906200e+00
            +5.36928700e-03 * tc[1]
            +6.30497333e-07 * tc[2]
            -1.78964575e-09 * tc[3]
            +5.73477000e-13 * tc[4]
            +1.52147660e+03 / tc[1];
        /*species 52: CH3CHO */
        species[52] =
            +3.72945950e+00
            -1.59664290e-03 * tc[1]
            +1.58449737e-05 * tc[2]
            -1.43646527e-08 * tc[3]
            +4.38622240e-12 * tc[4]
            -2.15728780e+04 / tc[1];
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
        /*species 30: N */
        species[30] =
            +1.41594290e+00
            +8.74453250e-05 * tc[1]
            -3.96745633e-08 * tc[2]
            +7.55656125e-12 * tc[3]
            -4.07219640e-16 * tc[4]
            +5.61337730e+04 / tc[1];
        /*species 31: NH */
        species[31] =
            +1.78369280e+00
            +6.64921500e-04 * tc[1]
            -1.41593490e-07 * tc[2]
            +1.95871253e-11 * tc[3]
            -1.10088940e-15 * tc[4]
            +4.21208480e+04 / tc[1];
        /*species 32: NH2 */
        species[32] =
            +1.83474210e+00
            +1.60365410e-03 * tc[1]
            -3.11302680e-07 * tc[2]
            +3.42573825e-11 * tc[3]
            -1.58412288e-15 * tc[4]
            +2.21719570e+04 / tc[1];
        /*species 33: NH3 */
        species[33] =
            +1.63445210e+00
            +2.83312800e-03 * tc[1]
            -5.75955867e-07 * tc[2]
            +5.96679025e-11 * tc[3]
            -2.51575720e-15 * tc[4]
            -6.54469580e+03 / tc[1];
        /*species 34: NNH */
        species[34] =
            +2.76675440e+00
            +1.44575410e-03 * tc[1]
            -3.47220667e-07 * tc[2]
            +4.21064850e-11 * tc[3]
            -2.01837920e-15 * tc[4]
            +2.86506970e+04 / tc[1];
        /*species 35: NO */
        species[35] =
            +2.26060560e+00
            +5.95552150e-04 * tc[1]
            -1.43056827e-07 * tc[2]
            +1.73644173e-11 * tc[3]
            -8.06721980e-16 * tc[4]
            +9.92097460e+03 / tc[1];
        /*species 36: NO2 */
        species[36] =
            +3.88475420e+00
            +1.08619780e-03 * tc[1]
            -2.76023020e-07 * tc[2]
            +3.93687750e-11 * tc[3]
            -2.10217900e-15 * tc[4]
            +2.31649830e+03 / tc[1];
        /*species 37: N2O */
        species[37] =
            +3.82307290e+00
            +1.31351255e-03 * tc[1]
            -3.19502913e-07 * tc[2]
            +4.00017800e-11 * tc[3]
            -1.95504606e-15 * tc[4]
            +8.07340480e+03 / tc[1];
        /*species 38: HNO */
        species[38] =
            +1.97925090e+00
            +1.74720295e-03 * tc[1]
            -2.61832593e-07 * tc[2]
            +1.43698985e-11 * tc[3]
            -3.86718320e-17 * tc[4]
            +1.17505820e+04 / tc[1];
        /*species 39: CN */
        species[39] =
            +2.74598050e+00
            +2.17253875e-05 * tc[1]
            +9.90199467e-08 * tc[2]
            -1.71629515e-11 * tc[3]
            +8.82683460e-16 * tc[4]
            +5.15361880e+04 / tc[1];
        /*species 40: HCN */
        species[40] =
            +2.80223920e+00
            +1.57321140e-03 * tc[1]
            -3.54406167e-07 * tc[2]
            +4.15493925e-11 * tc[3]
            -1.95995140e-15 * tc[4]
            +1.44072920e+04 / tc[1];
        /*species 41: H2CN */
        species[41] =
            +4.20970300e+00
            +1.48464555e-03 * tc[1]
            -9.51863033e-08 * tc[2]
            -4.08887500e-11 * tc[3]
            +6.08651780e-15 * tc[4]
            +2.76771090e+04 / tc[1];
        /*species 42: HCNN */
        species[42] =
            +4.89463620e+00
            +1.99479795e-03 * tc[1]
            -5.32746000e-07 * tc[2]
            +7.31234875e-11 * tc[3]
            -4.01893720e-15 * tc[4]
            +5.34529410e+04 / tc[1];
        /*species 46: NCO */
        species[46] =
            +4.15218450e+00
            +1.15258805e-03 * tc[1]
            -2.93443843e-07 * tc[2]
            +3.69727450e-11 * tc[3]
            -1.81955992e-15 * tc[4]
            +1.40041230e+04 / tc[1];
        /*species 47: N2 */
        species[47] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 / tc[1];
        /*species 48: AR */
        species[48] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 / tc[1];
        /*species 49: C3H7 */
        species[49] =
            +6.70269870e+00
            +8.02210150e-03 * tc[1]
            -1.76110733e-06 * tc[2]
            +1.90746475e-10 * tc[3]
            -7.87845680e-15 * tc[4]
            +8.29843360e+03 / tc[1];
        /*species 50: C3H8 */
        species[50] =
            +6.53413680e+00
            +9.43611950e-03 * tc[1]
            -2.09061637e-06 * tc[2]
            +2.28689123e-10 * tc[3]
            -9.56761380e-15 * tc[4]
            -1.64675160e+04 / tc[1];
        /*species 51: CH2CHO */
        species[51] =
            +4.97567000e+00
            +4.06529550e-03 * tc[1]
            -9.14541333e-07 * tc[2]
            +1.01757600e-10 * tc[3]
            -4.35203400e-15 * tc[4]
            +4.90321800e+02 / tc[1];
        /*species 52: CH3CHO */
        species[52] =
            +4.40411080e+00
            +5.86152950e-03 * tc[1]
            -1.40877123e-06 * tc[2]
            +1.70931128e-10 * tc[3]
            -8.19697260e-15 * tc[4]
            -2.25931220e+04 / tc[1];
    }

    /*species with midpoint at T=1368 kelvin */
    if (T < 1368) {
        /*species 44: HOCN */
        species[44] =
            +2.78604952e+00
            +3.44333961e-03 * tc[1]
            -1.07162621e-06 * tc[2]
            +1.29298942e-10 * tc[3]
            +2.38721576e-15 * tc[4]
            -2.82698400e+03 / tc[1];
    } else {
        /*species 44: HOCN */
        species[44] =
            +4.89784885e+00
            +1.58394696e-03 * tc[1]
            -3.72670213e-07 * tc[2]
            +4.43107860e-11 * tc[3]
            -2.08678354e-15 * tc[4]
            -3.70653331e+03 / tc[1];
    }

    /*species with midpoint at T=1478 kelvin */
    if (T < 1478) {
        /*species 45: HNCO */
        species[45] =
            +2.63096317e+00
            +3.65141179e-03 * tc[1]
            -7.60166677e-07 * tc[2]
            -1.65317825e-10 * tc[3]
            +7.24471504e-14 * tc[4]
            -1.55873636e+04 / tc[1];
    } else {
        /*species 45: HNCO */
        species[45] =
            +5.22395134e+00
            +1.58932002e-03 * tc[1]
            -3.64595850e-07 * tc[2]
            +4.26837908e-11 * tc[3]
            -1.99004391e-15 * tc[4]
            -1.66599344e+04 / tc[1];
    }

    /*species with midpoint at T=1382 kelvin */
    if (T < 1382) {
        /*species 43: HCNO */
        species[43] =
            +1.64727989e+00
            +6.37526710e-03 * tc[1]
            -3.49314120e-06 * tc[2]
            +1.10358209e-09 * tc[3]
            -1.51504293e-13 * tc[4]
            +1.92990252e+04 / tc[1];
    } else {
        /*species 43: HCNO */
        species[43] =
            +5.59860456e+00
            +1.51389313e-03 * tc[1]
            -3.59014487e-07 * tc[2]
            +4.29166320e-11 * tc[3]
            -2.02878782e-15 * tc[4]
            +1.79661339e+04 / tc[1];
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
        /*species 30: N */
        species[30] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +5.61046370e+04 / tc[1];
        /*species 31: NH */
        species[31] =
            +3.49290850e+00
            +1.55895990e-04 * tc[1]
            -4.96349467e-07 * tc[2]
            +6.20411050e-10 * tc[3]
            -2.07139340e-13 * tc[4]
            +4.18806290e+04 / tc[1];
        /*species 32: NH2 */
        species[32] =
            +4.20400290e+00
            -1.05306925e-03 * tc[1]
            +2.36894493e-06 * tc[2]
            -1.40287992e-09 * tc[3]
            +3.28814340e-13 * tc[4]
            +2.18859100e+04 / tc[1];
        /*species 33: NH3 */
        species[33] =
            +4.28602740e+00
            -2.33026150e-03 * tc[1]
            +7.23950433e-06 * tc[2]
            -5.70222175e-09 * tc[3]
            +1.65276092e-12 * tc[4]
            -6.74172850e+03 / tc[1];
        /*species 34: NNH */
        species[34] =
            +4.34469270e+00
            -2.42485360e-03 * tc[1]
            +6.68648633e-06 * tc[2]
            -5.43161600e-09 * tc[3]
            +1.58939078e-12 * tc[4]
            +2.87919730e+04 / tc[1];
        /*species 35: NO */
        species[35] =
            +4.21847630e+00
            -2.31948800e-03 * tc[1]
            +3.68034067e-06 * tc[2]
            -2.33403385e-09 * tc[3]
            +5.60715400e-13 * tc[4]
            +9.84462300e+03 / tc[1];
        /*species 36: NO2 */
        species[36] =
            +3.94403120e+00
            -7.92714500e-04 * tc[1]
            +5.55260400e-06 * tc[2]
            -5.11885650e-09 * tc[3]
            +1.56701128e-12 * tc[4]
            +2.89661790e+03 / tc[1];
        /*species 37: N2O */
        species[37] =
            +2.25715020e+00
            +5.65236400e-03 * tc[1]
            -4.55710633e-06 * tc[2]
            +2.42049515e-09 * tc[3]
            -5.86143640e-13 * tc[4]
            +8.74177440e+03 / tc[1];
        /*species 38: HNO */
        species[38] =
            +4.53349160e+00
            -2.83480855e-03 * tc[1]
            +6.15773567e-06 * tc[2]
            -4.28427350e-09 * tc[3]
            +1.10909146e-12 * tc[4]
            +1.15482970e+04 / tc[1];
        /*species 39: CN */
        species[39] =
            +3.61293510e+00
            -4.77756635e-04 * tc[1]
            +7.14765900e-07 * tc[2]
            -7.87908075e-11 * tc[3]
            -9.28607120e-14 * tc[4]
            +5.17083400e+04 / tc[1];
        /*species 40: HCN */
        species[40] =
            +2.25898860e+00
            +5.02558500e-03 * tc[1]
            -4.45058767e-06 * tc[2]
            +2.52308725e-09 * tc[3]
            -6.01780560e-13 * tc[4]
            +1.47126330e+04 / tc[1];
        /*species 41: H2CN */
        species[41] =
            +2.85166100e+00
            +2.84761655e-03 * tc[1]
            +3.57046667e-07 * tc[2]
            -4.05653000e-10 * tc[3]
            -4.70221620e-14 * tc[4]
            +2.86378200e+04 / tc[1];
        /*species 42: HCNN */
        species[42] =
            +2.52431940e+00
            +7.98030950e-03 * tc[1]
            -6.27211800e-06 * tc[2]
            +3.03138500e-09 * tc[3]
            -6.47147560e-13 * tc[4]
            +5.42619840e+04 / tc[1];
        /*species 46: NCO */
        species[46] =
            +2.82693080e+00
            +4.40258440e-03 * tc[1]
            -2.79553780e-06 * tc[2]
            +1.20042410e-09 * tc[3]
            -2.66271900e-13 * tc[4]
            +1.46824770e+04 / tc[1];
        /*species 47: N2 */
        species[47] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 / tc[1];
        /*species 48: AR */
        species[48] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 / tc[1];
        /*species 49: C3H7 */
        species[49] =
            +1.05155180e+00
            +1.29959900e-02 * tc[1]
            +7.93351333e-07 * tc[2]
            -4.90239225e-09 * tc[3]
            +1.87464940e-12 * tc[4]
            +1.06318630e+04 / tc[1];
        /*species 50: C3H8 */
        species[50] =
            +9.33553810e-01
            +1.32122895e-02 * tc[1]
            +2.03532423e-06 * tc[2]
            -5.49437475e-09 * tc[3]
            +1.90298506e-12 * tc[4]
            -1.39585200e+04 / tc[1];
        /*species 51: CH2CHO */
        species[51] =
            +3.40906200e+00
            +5.36928700e-03 * tc[1]
            +6.30497333e-07 * tc[2]
            -1.78964575e-09 * tc[3]
            +5.73477000e-13 * tc[4]
            +1.52147660e+03 / tc[1];
        /*species 52: CH3CHO */
        species[52] =
            +4.72945950e+00
            -1.59664290e-03 * tc[1]
            +1.58449737e-05 * tc[2]
            -1.43646527e-08 * tc[3]
            +4.38622240e-12 * tc[4]
            -2.15728780e+04 / tc[1];
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
        /*species 30: N */
        species[30] =
            +2.41594290e+00
            +8.74453250e-05 * tc[1]
            -3.96745633e-08 * tc[2]
            +7.55656125e-12 * tc[3]
            -4.07219640e-16 * tc[4]
            +5.61337730e+04 / tc[1];
        /*species 31: NH */
        species[31] =
            +2.78369280e+00
            +6.64921500e-04 * tc[1]
            -1.41593490e-07 * tc[2]
            +1.95871253e-11 * tc[3]
            -1.10088940e-15 * tc[4]
            +4.21208480e+04 / tc[1];
        /*species 32: NH2 */
        species[32] =
            +2.83474210e+00
            +1.60365410e-03 * tc[1]
            -3.11302680e-07 * tc[2]
            +3.42573825e-11 * tc[3]
            -1.58412288e-15 * tc[4]
            +2.21719570e+04 / tc[1];
        /*species 33: NH3 */
        species[33] =
            +2.63445210e+00
            +2.83312800e-03 * tc[1]
            -5.75955867e-07 * tc[2]
            +5.96679025e-11 * tc[3]
            -2.51575720e-15 * tc[4]
            -6.54469580e+03 / tc[1];
        /*species 34: NNH */
        species[34] =
            +3.76675440e+00
            +1.44575410e-03 * tc[1]
            -3.47220667e-07 * tc[2]
            +4.21064850e-11 * tc[3]
            -2.01837920e-15 * tc[4]
            +2.86506970e+04 / tc[1];
        /*species 35: NO */
        species[35] =
            +3.26060560e+00
            +5.95552150e-04 * tc[1]
            -1.43056827e-07 * tc[2]
            +1.73644173e-11 * tc[3]
            -8.06721980e-16 * tc[4]
            +9.92097460e+03 / tc[1];
        /*species 36: NO2 */
        species[36] =
            +4.88475420e+00
            +1.08619780e-03 * tc[1]
            -2.76023020e-07 * tc[2]
            +3.93687750e-11 * tc[3]
            -2.10217900e-15 * tc[4]
            +2.31649830e+03 / tc[1];
        /*species 37: N2O */
        species[37] =
            +4.82307290e+00
            +1.31351255e-03 * tc[1]
            -3.19502913e-07 * tc[2]
            +4.00017800e-11 * tc[3]
            -1.95504606e-15 * tc[4]
            +8.07340480e+03 / tc[1];
        /*species 38: HNO */
        species[38] =
            +2.97925090e+00
            +1.74720295e-03 * tc[1]
            -2.61832593e-07 * tc[2]
            +1.43698985e-11 * tc[3]
            -3.86718320e-17 * tc[4]
            +1.17505820e+04 / tc[1];
        /*species 39: CN */
        species[39] =
            +3.74598050e+00
            +2.17253875e-05 * tc[1]
            +9.90199467e-08 * tc[2]
            -1.71629515e-11 * tc[3]
            +8.82683460e-16 * tc[4]
            +5.15361880e+04 / tc[1];
        /*species 40: HCN */
        species[40] =
            +3.80223920e+00
            +1.57321140e-03 * tc[1]
            -3.54406167e-07 * tc[2]
            +4.15493925e-11 * tc[3]
            -1.95995140e-15 * tc[4]
            +1.44072920e+04 / tc[1];
        /*species 41: H2CN */
        species[41] =
            +5.20970300e+00
            +1.48464555e-03 * tc[1]
            -9.51863033e-08 * tc[2]
            -4.08887500e-11 * tc[3]
            +6.08651780e-15 * tc[4]
            +2.76771090e+04 / tc[1];
        /*species 42: HCNN */
        species[42] =
            +5.89463620e+00
            +1.99479795e-03 * tc[1]
            -5.32746000e-07 * tc[2]
            +7.31234875e-11 * tc[3]
            -4.01893720e-15 * tc[4]
            +5.34529410e+04 / tc[1];
        /*species 46: NCO */
        species[46] =
            +5.15218450e+00
            +1.15258805e-03 * tc[1]
            -2.93443843e-07 * tc[2]
            +3.69727450e-11 * tc[3]
            -1.81955992e-15 * tc[4]
            +1.40041230e+04 / tc[1];
        /*species 47: N2 */
        species[47] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 / tc[1];
        /*species 48: AR */
        species[48] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 / tc[1];
        /*species 49: C3H7 */
        species[49] =
            +7.70269870e+00
            +8.02210150e-03 * tc[1]
            -1.76110733e-06 * tc[2]
            +1.90746475e-10 * tc[3]
            -7.87845680e-15 * tc[4]
            +8.29843360e+03 / tc[1];
        /*species 50: C3H8 */
        species[50] =
            +7.53413680e+00
            +9.43611950e-03 * tc[1]
            -2.09061637e-06 * tc[2]
            +2.28689123e-10 * tc[3]
            -9.56761380e-15 * tc[4]
            -1.64675160e+04 / tc[1];
        /*species 51: CH2CHO */
        species[51] =
            +5.97567000e+00
            +4.06529550e-03 * tc[1]
            -9.14541333e-07 * tc[2]
            +1.01757600e-10 * tc[3]
            -4.35203400e-15 * tc[4]
            +4.90321800e+02 / tc[1];
        /*species 52: CH3CHO */
        species[52] =
            +5.40411080e+00
            +5.86152950e-03 * tc[1]
            -1.40877123e-06 * tc[2]
            +1.70931128e-10 * tc[3]
            -8.19697260e-15 * tc[4]
            -2.25931220e+04 / tc[1];
    }

    /*species with midpoint at T=1368 kelvin */
    if (T < 1368) {
        /*species 44: HOCN */
        species[44] =
            +3.78604952e+00
            +3.44333961e-03 * tc[1]
            -1.07162621e-06 * tc[2]
            +1.29298942e-10 * tc[3]
            +2.38721576e-15 * tc[4]
            -2.82698400e+03 / tc[1];
    } else {
        /*species 44: HOCN */
        species[44] =
            +5.89784885e+00
            +1.58394696e-03 * tc[1]
            -3.72670213e-07 * tc[2]
            +4.43107860e-11 * tc[3]
            -2.08678354e-15 * tc[4]
            -3.70653331e+03 / tc[1];
    }

    /*species with midpoint at T=1478 kelvin */
    if (T < 1478) {
        /*species 45: HNCO */
        species[45] =
            +3.63096317e+00
            +3.65141179e-03 * tc[1]
            -7.60166677e-07 * tc[2]
            -1.65317825e-10 * tc[3]
            +7.24471504e-14 * tc[4]
            -1.55873636e+04 / tc[1];
    } else {
        /*species 45: HNCO */
        species[45] =
            +6.22395134e+00
            +1.58932002e-03 * tc[1]
            -3.64595850e-07 * tc[2]
            +4.26837908e-11 * tc[3]
            -1.99004391e-15 * tc[4]
            -1.66599344e+04 / tc[1];
    }

    /*species with midpoint at T=1382 kelvin */
    if (T < 1382) {
        /*species 43: HCNO */
        species[43] =
            +2.64727989e+00
            +6.37526710e-03 * tc[1]
            -3.49314120e-06 * tc[2]
            +1.10358209e-09 * tc[3]
            -1.51504293e-13 * tc[4]
            +1.92990252e+04 / tc[1];
    } else {
        /*species 43: HCNO */
        species[43] =
            +6.59860456e+00
            +1.51389313e-03 * tc[1]
            -3.59014487e-07 * tc[2]
            +4.29166320e-11 * tc[3]
            -2.02878782e-15 * tc[4]
            +1.79661339e+04 / tc[1];
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
        /*species 30: N */
        species[30] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.19390870e+00 ;
        /*species 31: NH */
        species[31] =
            +3.49290850e+00 * tc[0]
            +3.11791980e-04 * tc[1]
            -7.44524200e-07 * tc[2]
            +8.27214733e-10 * tc[3]
            -2.58924175e-13 * tc[4]
            +1.84832780e+00 ;
        /*species 32: NH2 */
        species[32] =
            +4.20400290e+00 * tc[0]
            -2.10613850e-03 * tc[1]
            +3.55341740e-06 * tc[2]
            -1.87050657e-09 * tc[3]
            +4.11017925e-13 * tc[4]
            -1.41842480e-01 ;
        /*species 33: NH3 */
        species[33] =
            +4.28602740e+00 * tc[0]
            -4.66052300e-03 * tc[1]
            +1.08592565e-05 * tc[2]
            -7.60296233e-09 * tc[3]
            +2.06595115e-12 * tc[4]
            -6.25372770e-01 ;
        /*species 34: NNH */
        species[34] =
            +4.34469270e+00 * tc[0]
            -4.84970720e-03 * tc[1]
            +1.00297295e-05 * tc[2]
            -7.24215467e-09 * tc[3]
            +1.98673848e-12 * tc[4]
            +2.97794100e+00 ;
        /*species 35: NO */
        species[35] =
            +4.21847630e+00 * tc[0]
            -4.63897600e-03 * tc[1]
            +5.52051100e-06 * tc[2]
            -3.11204513e-09 * tc[3]
            +7.00894250e-13 * tc[4]
            +2.28084640e+00 ;
        /*species 36: NO2 */
        species[36] =
            +3.94403120e+00 * tc[0]
            -1.58542900e-03 * tc[1]
            +8.32890600e-06 * tc[2]
            -6.82514200e-09 * tc[3]
            +1.95876410e-12 * tc[4]
            +6.31199170e+00 ;
        /*species 37: N2O */
        species[37] =
            +2.25715020e+00 * tc[0]
            +1.13047280e-02 * tc[1]
            -6.83565950e-06 * tc[2]
            +3.22732687e-09 * tc[3]
            -7.32679550e-13 * tc[4]
            +1.07579920e+01 ;
        /*species 38: HNO */
        species[38] =
            +4.53349160e+00 * tc[0]
            -5.66961710e-03 * tc[1]
            +9.23660350e-06 * tc[2]
            -5.71236467e-09 * tc[3]
            +1.38636433e-12 * tc[4]
            +1.74984170e+00 ;
        /*species 39: CN */
        species[39] =
            +3.61293510e+00 * tc[0]
            -9.55513270e-04 * tc[1]
            +1.07214885e-06 * tc[2]
            -1.05054410e-10 * tc[3]
            -1.16075890e-13 * tc[4]
            +3.98049950e+00 ;
        /*species 40: HCN */
        species[40] =
            +2.25898860e+00 * tc[0]
            +1.00511700e-02 * tc[1]
            -6.67588150e-06 * tc[2]
            +3.36411633e-09 * tc[3]
            -7.52225700e-13 * tc[4]
            +8.91644190e+00 ;
        /*species 41: H2CN */
        species[41] =
            +2.85166100e+00 * tc[0]
            +5.69523310e-03 * tc[1]
            +5.35570000e-07 * tc[2]
            -5.40870667e-10 * tc[3]
            -5.87777025e-14 * tc[4]
            +8.99275110e+00 ;
        /*species 42: HCNN */
        species[42] =
            +2.52431940e+00 * tc[0]
            +1.59606190e-02 * tc[1]
            -9.40817700e-06 * tc[2]
            +4.04184667e-09 * tc[3]
            -8.08934450e-13 * tc[4]
            +1.16758700e+01 ;
        /*species 46: NCO */
        species[46] =
            +2.82693080e+00 * tc[0]
            +8.80516880e-03 * tc[1]
            -4.19330670e-06 * tc[2]
            +1.60056547e-09 * tc[3]
            -3.32839875e-13 * tc[4]
            +9.55046460e+00 ;
        /*species 47: N2 */
        species[47] =
            +3.29867700e+00 * tc[0]
            +1.40824040e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213500e-13 * tc[4]
            +3.95037200e+00 ;
        /*species 48: AR */
        species[48] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.36600000e+00 ;
        /*species 49: C3H7 */
        species[49] =
            +1.05155180e+00 * tc[0]
            +2.59919800e-02 * tc[1]
            +1.19002700e-06 * tc[2]
            -6.53652300e-09 * tc[3]
            +2.34331175e-12 * tc[4]
            +2.11225590e+01 ;
        /*species 50: C3H8 */
        species[50] =
            +9.33553810e-01 * tc[0]
            +2.64245790e-02 * tc[1]
            +3.05298635e-06 * tc[2]
            -7.32583300e-09 * tc[3]
            +2.37873132e-12 * tc[4]
            +1.92016910e+01 ;
        /*species 51: CH2CHO */
        species[51] =
            +3.40906200e+00 * tc[0]
            +1.07385740e-02 * tc[1]
            +9.45746000e-07 * tc[2]
            -2.38619433e-09 * tc[3]
            +7.16846250e-13 * tc[4]
            +9.55829000e+00 ;
        /*species 52: CH3CHO */
        species[52] =
            +4.72945950e+00 * tc[0]
            -3.19328580e-03 * tc[1]
            +2.37674605e-05 * tc[2]
            -1.91528703e-08 * tc[3]
            +5.48277800e-12 * tc[4]
            +4.10301590e+00 ;
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
        /*species 30: N */
        species[30] =
            +2.41594290e+00 * tc[0]
            +1.74890650e-04 * tc[1]
            -5.95118450e-08 * tc[2]
            +1.00754150e-11 * tc[3]
            -5.09024550e-16 * tc[4]
            +4.64960960e+00 ;
        /*species 31: NH */
        species[31] =
            +2.78369280e+00 * tc[0]
            +1.32984300e-03 * tc[1]
            -2.12390235e-07 * tc[2]
            +2.61161670e-11 * tc[3]
            -1.37611175e-15 * tc[4]
            +5.74077990e+00 ;
        /*species 32: NH2 */
        species[32] =
            +2.83474210e+00 * tc[0]
            +3.20730820e-03 * tc[1]
            -4.66954020e-07 * tc[2]
            +4.56765100e-11 * tc[3]
            -1.98015360e-15 * tc[4]
            +6.52041630e+00 ;
        /*species 33: NH3 */
        species[33] =
            +2.63445210e+00 * tc[0]
            +5.66625600e-03 * tc[1]
            -8.63933800e-07 * tc[2]
            +7.95572033e-11 * tc[3]
            -3.14469650e-15 * tc[4]
            +6.56629280e+00 ;
        /*species 34: NNH */
        species[34] =
            +3.76675440e+00 * tc[0]
            +2.89150820e-03 * tc[1]
            -5.20831000e-07 * tc[2]
            +5.61419800e-11 * tc[3]
            -2.52297400e-15 * tc[4]
            +4.47050670e+00 ;
        /*species 35: NO */
        species[35] =
            +3.26060560e+00 * tc[0]
            +1.19110430e-03 * tc[1]
            -2.14585240e-07 * tc[2]
            +2.31525563e-11 * tc[3]
            -1.00840247e-15 * tc[4]
            +6.36930270e+00 ;
        /*species 36: NO2 */
        species[36] =
            +4.88475420e+00 * tc[0]
            +2.17239560e-03 * tc[1]
            -4.14034530e-07 * tc[2]
            +5.24917000e-11 * tc[3]
            -2.62772375e-15 * tc[4]
            -1.17416950e-01 ;
        /*species 37: N2O */
        species[37] =
            +4.82307290e+00 * tc[0]
            +2.62702510e-03 * tc[1]
            -4.79254370e-07 * tc[2]
            +5.33357067e-11 * tc[3]
            -2.44380757e-15 * tc[4]
            -2.20172070e+00 ;
        /*species 38: HNO */
        species[38] =
            +2.97925090e+00 * tc[0]
            +3.49440590e-03 * tc[1]
            -3.92748890e-07 * tc[2]
            +1.91598647e-11 * tc[3]
            -4.83397900e-17 * tc[4]
            +8.60637280e+00 ;
        /*species 39: CN */
        species[39] =
            +3.74598050e+00 * tc[0]
            +4.34507750e-05 * tc[1]
            +1.48529920e-07 * tc[2]
            -2.28839353e-11 * tc[3]
            +1.10335433e-15 * tc[4]
            +2.78676010e+00 ;
        /*species 40: HCN */
        species[40] =
            +3.80223920e+00 * tc[0]
            +3.14642280e-03 * tc[1]
            -5.31609250e-07 * tc[2]
            +5.53991900e-11 * tc[3]
            -2.44993925e-15 * tc[4]
            +1.57546010e+00 ;
        /*species 41: H2CN */
        species[41] =
            +5.20970300e+00 * tc[0]
            +2.96929110e-03 * tc[1]
            -1.42779455e-07 * tc[2]
            -5.45183333e-11 * tc[3]
            +7.60814725e-15 * tc[4]
            -4.44447800e+00 ;
        /*species 42: HCNN */
        species[42] =
            +5.89463620e+00 * tc[0]
            +3.98959590e-03 * tc[1]
            -7.99119000e-07 * tc[2]
            +9.74979833e-11 * tc[3]
            -5.02367150e-15 * tc[4]
            -5.10305020e+00 ;
        /*species 46: NCO */
        species[46] =
            +5.15218450e+00 * tc[0]
            +2.30517610e-03 * tc[1]
            -4.40165765e-07 * tc[2]
            +4.92969933e-11 * tc[3]
            -2.27444990e-15 * tc[4]
            -2.54426600e+00 ;
        /*species 47: N2 */
        species[47] =
            +2.92664000e+00 * tc[0]
            +1.48797680e-03 * tc[1]
            -2.84238000e-07 * tc[2]
            +3.36567933e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
        /*species 48: AR */
        species[48] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.36600000e+00 ;
        /*species 49: C3H7 */
        species[49] =
            +7.70269870e+00 * tc[0]
            +1.60442030e-02 * tc[1]
            -2.64166100e-06 * tc[2]
            +2.54328633e-10 * tc[3]
            -9.84807100e-15 * tc[4]
            -1.54801800e+01 ;
        /*species 50: C3H8 */
        species[50] =
            +7.53413680e+00 * tc[0]
            +1.88722390e-02 * tc[1]
            -3.13592455e-06 * tc[2]
            +3.04918830e-10 * tc[3]
            -1.19595173e-14 * tc[4]
            -1.78923490e+01 ;
        /*species 51: CH2CHO */
        species[51] =
            +5.97567000e+00 * tc[0]
            +8.13059100e-03 * tc[1]
            -1.37181200e-06 * tc[2]
            +1.35676800e-10 * tc[3]
            -5.44004250e-15 * tc[4]
            -5.04525100e+00 ;
        /*species 52: CH3CHO */
        species[52] =
            +5.40411080e+00 * tc[0]
            +1.17230590e-02 * tc[1]
            -2.11315685e-06 * tc[2]
            +2.27908170e-10 * tc[3]
            -1.02462158e-14 * tc[4]
            -3.48079170e+00 ;
    }

    /*species with midpoint at T=1368 kelvin */
    if (T < 1368) {
        /*species 44: HOCN */
        species[44] =
            +3.78604952e+00 * tc[0]
            +6.88667922e-03 * tc[1]
            -1.60743932e-06 * tc[2]
            +1.72398589e-10 * tc[3]
            +2.98401970e-15 * tc[4]
            +5.63292162e+00 ;
    } else {
        /*species 44: HOCN */
        species[44] =
            +5.89784885e+00 * tc[0]
            +3.16789393e-03 * tc[1]
            -5.59005320e-07 * tc[2]
            +5.90810480e-11 * tc[3]
            -2.60847942e-15 * tc[4]
            -6.18167825e+00 ;
    }

    /*species with midpoint at T=1478 kelvin */
    if (T < 1478) {
        /*species 45: HNCO */
        species[45] =
            +3.63096317e+00 * tc[0]
            +7.30282357e-03 * tc[1]
            -1.14025001e-06 * tc[2]
            -2.20423766e-10 * tc[3]
            +9.05589380e-14 * tc[4]
            +6.19457727e+00 ;
    } else {
        /*species 45: HNCO */
        species[45] =
            +6.22395134e+00 * tc[0]
            +3.17864004e-03 * tc[1]
            -5.46893775e-07 * tc[2]
            +5.69117210e-11 * tc[3]
            -2.48755489e-15 * tc[4]
            -8.38224741e+00 ;
    }

    /*species with midpoint at T=1382 kelvin */
    if (T < 1382) {
        /*species 43: HCNO */
        species[43] =
            +2.64727989e+00 * tc[0]
            +1.27505342e-02 * tc[1]
            -5.23971180e-06 * tc[2]
            +1.47144279e-09 * tc[3]
            -1.89380367e-13 * tc[4]
            +1.07332972e+01 ;
    } else {
        /*species 43: HCNO */
        species[43] =
            +6.59860456e+00 * tc[0]
            +3.02778626e-03 * tc[1]
            -5.38521730e-07 * tc[2]
            +5.72221760e-11 * tc[3]
            -2.53598478e-15 * tc[4]
            -1.03306599e+01 ;
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
    wt[30] = 14.006700; /*N */
    wt[31] = 15.014670; /*NH */
    wt[32] = 16.022640; /*NH2 */
    wt[33] = 17.030610; /*NH3 */
    wt[34] = 29.021370; /*NNH */
    wt[35] = 30.006100; /*NO */
    wt[36] = 46.005500; /*NO2 */
    wt[37] = 44.012800; /*N2O */
    wt[38] = 31.014070; /*HNO */
    wt[39] = 26.017700; /*CN */
    wt[40] = 27.025670; /*HCN */
    wt[41] = 28.033640; /*H2CN */
    wt[42] = 41.032370; /*HCNN */
    wt[43] = 43.025070; /*HCNO */
    wt[44] = 43.025070; /*HOCN */
    wt[45] = 43.025070; /*HNCO */
    wt[46] = 42.017100; /*NCO */
    wt[47] = 28.013400; /*N2 */
    wt[48] = 39.948000; /*AR */
    wt[49] = 43.088790; /*C3H7 */
    wt[50] = 44.096760; /*C3H8 */
    wt[51] = 43.045310; /*CH2CHO */
    wt[52] = 44.053280; /*CH3CHO */

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
    y[30] = phi[30]*14.006700;   XW += y[30]; /*N */
    y[31] = phi[31]*15.014670;   XW += y[31]; /*NH */
    y[32] = phi[32]*16.022640;   XW += y[32]; /*NH2 */
    y[33] = phi[33]*17.030610;   XW += y[33]; /*NH3 */
    y[34] = phi[34]*29.021370;   XW += y[34]; /*NNH */
    y[35] = phi[35]*30.006100;   XW += y[35]; /*NO */
    y[36] = phi[36]*46.005500;   XW += y[36]; /*NO2 */
    y[37] = phi[37]*44.012800;   XW += y[37]; /*N2O */
    y[38] = phi[38]*31.014070;   XW += y[38]; /*HNO */
    y[39] = phi[39]*26.017700;   XW += y[39]; /*CN */
    y[40] = phi[40]*27.025670;   XW += y[40]; /*HCN */
    y[41] = phi[41]*28.033640;   XW += y[41]; /*H2CN */
    y[42] = phi[42]*41.032370;   XW += y[42]; /*HCNN */
    y[43] = phi[43]*43.025070;   XW += y[43]; /*HCNO */
    y[44] = phi[44]*43.025070;   XW += y[44]; /*HOCN */
    y[45] = phi[45]*43.025070;   XW += y[45]; /*HNCO */
    y[46] = phi[46]*42.017100;   XW += y[46]; /*NCO */
    y[47] = phi[47]*28.013400;   XW += y[47]; /*N2 */
    y[48] = phi[48]*39.948000;   XW += y[48]; /*AR */
    y[49] = phi[49]*43.088790;   XW += y[49]; /*C3H7 */
    y[50] = phi[50]*44.096760;   XW += y[50]; /*C3H8 */
    y[51] = phi[51]*43.045310;   XW += y[51]; /*CH2CHO */
    y[52] = phi[52]*44.053280;   XW += y[52]; /*CH3CHO */
    for (id = 0; id < 53; ++id) {
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
    phi[30] = y[30]/ 1.40067000e-02; /*N (wt in kg) */
    phi[31] = y[31]/ 1.50146700e-02; /*NH (wt in kg) */
    phi[32] = y[32]/ 1.60226400e-02; /*NH2 (wt in kg) */
    phi[33] = y[33]/ 1.70306100e-02; /*NH3 (wt in kg) */
    phi[34] = y[34]/ 2.90213700e-02; /*NNH (wt in kg) */
    phi[35] = y[35]/ 3.00061000e-02; /*NO (wt in kg) */
    phi[36] = y[36]/ 4.60055000e-02; /*NO2 (wt in kg) */
    phi[37] = y[37]/ 4.40128000e-02; /*N2O (wt in kg) */
    phi[38] = y[38]/ 3.10140700e-02; /*HNO (wt in kg) */
    phi[39] = y[39]/ 2.60177000e-02; /*CN (wt in kg) */
    phi[40] = y[40]/ 2.70256700e-02; /*HCN (wt in kg) */
    phi[41] = y[41]/ 2.80336400e-02; /*H2CN (wt in kg) */
    phi[42] = y[42]/ 4.10323700e-02; /*HCNN (wt in kg) */
    phi[43] = y[43]/ 4.30250700e-02; /*HCNO (wt in kg) */
    phi[44] = y[44]/ 4.30250700e-02; /*HOCN (wt in kg) */
    phi[45] = y[45]/ 4.30250700e-02; /*HNCO (wt in kg) */
    phi[46] = y[46]/ 4.20171000e-02; /*NCO (wt in kg) */
    phi[47] = y[47]/ 2.80134000e-02; /*N2 (wt in kg) */
    phi[48] = y[48]/ 3.99480000e-02; /*AR (wt in kg) */
    phi[49] = y[49]/ 4.30887900e-02; /*C3H7 (wt in kg) */
    phi[50] = y[50]/ 4.40967600e-02; /*C3H8 (wt in kg) */
    phi[51] = y[51]/ 4.30453100e-02; /*CH2CHO (wt in kg) */
    phi[52] = y[52]/ 4.40532800e-02; /*CH3CHO (wt in kg) */

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
    y[30] = c[30] * 14.006700 / (*rho); 
    y[31] = c[31] * 15.014670 / (*rho); 
    y[32] = c[32] * 16.022640 / (*rho); 
    y[33] = c[33] * 17.030610 / (*rho); 
    y[34] = c[34] * 29.021370 / (*rho); 
    y[35] = c[35] * 30.006100 / (*rho); 
    y[36] = c[36] * 46.005500 / (*rho); 
    y[37] = c[37] * 44.012800 / (*rho); 
    y[38] = c[38] * 31.014070 / (*rho); 
    y[39] = c[39] * 26.017700 / (*rho); 
    y[40] = c[40] * 27.025670 / (*rho); 
    y[41] = c[41] * 28.033640 / (*rho); 
    y[42] = c[42] * 41.032370 / (*rho); 
    y[43] = c[43] * 43.025070 / (*rho); 
    y[44] = c[44] * 43.025070 / (*rho); 
    y[45] = c[45] * 43.025070 / (*rho); 
    y[46] = c[46] * 42.017100 / (*rho); 
    y[47] = c[47] * 28.013400 / (*rho); 
    y[48] = c[48] * 39.948000 / (*rho); 
    y[49] = c[49] * 43.088790 / (*rho); 
    y[50] = c[50] * 44.096760 / (*rho); 
    y[51] = c[51] * 43.045310 / (*rho); 
    y[52] = c[52] * 44.053280 / (*rho); 

    return;
}


/*ddebdf compatible right hand side of CV burner */
/*rwrk[0] and rwrk[1] should contain rho and ene respectively */
/*working variable phi contains specific mole numbers */
void fecvrhs_(double * time, double * phi, double * phidot, double * rwrk, int * iwrk)
{
    double rho,ene; /*CV Parameters */
    double y[53], wdot[53]; /*temporary storage */
    int i; /*Loop counter */
    double temperature,pressure; /*temporary var */
    rho = rwrk[0];
    ene = rwrk[1];
    fephity_(phi, iwrk, rwrk, y);
    feeytt_(&ene, y, iwrk, rwrk, &temperature);
    CKPY(&rho, &temperature,  y, iwrk, rwrk, &pressure);
    CKWYP(&pressure, &temperature,  y, iwrk, rwrk, wdot);
    for (i=0; i<53; ++i) phidot[i] = wdot[i] / (rho/1000.0); 

    return;
}


/*returns the dimensionality of the cv burner (number of species) */
int fecvdim_()
{
    return 53;
}


/*ddebdf compatible right hand side of ZND solver */
/*rwrk[0] : scaling factor for pressure */
/*rwrk[1] : preshock density (g/cc)  */
/*rwrk[2] : detonation velocity (cm/s)  */
/*solution vector: [P; rho; y0 ... ylast]  */
void fezndrhs_(double * time, double * z, double * zdot, double * rwrk, int * iwrk)
{
    double psc,rho1,udet; /*ZND Parameters */
    double wt[53], hms[53], wdot[53]; /*temporary storage */
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
    for (i=0; i<53; ++i) {
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
    return 56;
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
    if (strcmp(s1, "N")==0) return 30; 
    if (strcmp(s1, "NH")==0) return 31; 
    if (strcmp(s1, "NH2")==0) return 32; 
    if (strcmp(s1, "NH3")==0) return 33; 
    if (strcmp(s1, "NNH")==0) return 34; 
    if (strcmp(s1, "NO")==0) return 35; 
    if (strcmp(s1, "NO2")==0) return 36; 
    if (strcmp(s1, "N2O")==0) return 37; 
    if (strcmp(s1, "HNO")==0) return 38; 
    if (strcmp(s1, "CN")==0) return 39; 
    if (strcmp(s1, "HCN")==0) return 40; 
    if (strcmp(s1, "H2CN")==0) return 41; 
    if (strcmp(s1, "HCNN")==0) return 42; 
    if (strcmp(s1, "HCNO")==0) return 43; 
    if (strcmp(s1, "HOCN")==0) return 44; 
    if (strcmp(s1, "HNCO")==0) return 45; 
    if (strcmp(s1, "NCO")==0) return 46; 
    if (strcmp(s1, "N2")==0) return 47; 
    if (strcmp(s1, "AR")==0) return 48; 
    if (strcmp(s1, "C3H7")==0) return 49; 
    if (strcmp(s1, "C3H8")==0) return 50; 
    if (strcmp(s1, "CH2CHO")==0) return 51; 
    if (strcmp(s1, "CH3CHO")==0) return 52; 
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
    if (sn==30) return "N"; 
    if (sn==31) return "NH"; 
    if (sn==32) return "NH2"; 
    if (sn==33) return "NH3"; 
    if (sn==34) return "NNH"; 
    if (sn==35) return "NO"; 
    if (sn==36) return "NO2"; 
    if (sn==37) return "N2O"; 
    if (sn==38) return "HNO"; 
    if (sn==39) return "CN"; 
    if (sn==40) return "HCN"; 
    if (sn==41) return "H2CN"; 
    if (sn==42) return "HCNN"; 
    if (sn==43) return "HCNO"; 
    if (sn==44) return "HOCN"; 
    if (sn==45) return "HNCO"; 
    if (sn==46) return "NCO"; 
    if (sn==47) return "N2"; 
    if (sn==48) return "AR"; 
    if (sn==49) return "C3H7"; 
    if (sn==50) return "C3H8"; 
    if (sn==51) return "CH2CHO"; 
    if (sn==52) return "CH3CHO"; 
    /*species name not found */
    return "NOTFOUND";
}

/* End of file  */
