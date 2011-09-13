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
    *mm = 5;
    *kk = 21;
    *ii = 84;
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
    kname[ 4*lenkname + 2 ] = ' ';

}


/* Returns the char strings of species names */
void CKSYMS(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*21; i++) {
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

    /* CH2  */
    kname[ 7*lenkname + 0 ] = 'C';
    kname[ 7*lenkname + 1 ] = 'H';
    kname[ 7*lenkname + 2 ] = '2';
    kname[ 7*lenkname + 3 ] = ' ';

    /* CH2(S)  */
    kname[ 8*lenkname + 0 ] = 'C';
    kname[ 8*lenkname + 1 ] = 'H';
    kname[ 8*lenkname + 2 ] = '2';
    kname[ 8*lenkname + 3 ] = '(';
    kname[ 8*lenkname + 4 ] = 'S';
    kname[ 8*lenkname + 5 ] = ')';
    kname[ 8*lenkname + 6 ] = ' ';

    /* CH3  */
    kname[ 9*lenkname + 0 ] = 'C';
    kname[ 9*lenkname + 1 ] = 'H';
    kname[ 9*lenkname + 2 ] = '3';
    kname[ 9*lenkname + 3 ] = ' ';

    /* CH4  */
    kname[ 10*lenkname + 0 ] = 'C';
    kname[ 10*lenkname + 1 ] = 'H';
    kname[ 10*lenkname + 2 ] = '4';
    kname[ 10*lenkname + 3 ] = ' ';

    /* CO  */
    kname[ 11*lenkname + 0 ] = 'C';
    kname[ 11*lenkname + 1 ] = 'O';
    kname[ 11*lenkname + 2 ] = ' ';

    /* CO2  */
    kname[ 12*lenkname + 0 ] = 'C';
    kname[ 12*lenkname + 1 ] = 'O';
    kname[ 12*lenkname + 2 ] = '2';
    kname[ 12*lenkname + 3 ] = ' ';

    /* HCO  */
    kname[ 13*lenkname + 0 ] = 'H';
    kname[ 13*lenkname + 1 ] = 'C';
    kname[ 13*lenkname + 2 ] = 'O';
    kname[ 13*lenkname + 3 ] = ' ';

    /* CH2O  */
    kname[ 14*lenkname + 0 ] = 'C';
    kname[ 14*lenkname + 1 ] = 'H';
    kname[ 14*lenkname + 2 ] = '2';
    kname[ 14*lenkname + 3 ] = 'O';
    kname[ 14*lenkname + 4 ] = ' ';

    /* CH3O  */
    kname[ 15*lenkname + 0 ] = 'C';
    kname[ 15*lenkname + 1 ] = 'H';
    kname[ 15*lenkname + 2 ] = '3';
    kname[ 15*lenkname + 3 ] = 'O';
    kname[ 15*lenkname + 4 ] = ' ';

    /* C2H4  */
    kname[ 16*lenkname + 0 ] = 'C';
    kname[ 16*lenkname + 1 ] = '2';
    kname[ 16*lenkname + 2 ] = 'H';
    kname[ 16*lenkname + 3 ] = '4';
    kname[ 16*lenkname + 4 ] = ' ';

    /* C2H5  */
    kname[ 17*lenkname + 0 ] = 'C';
    kname[ 17*lenkname + 1 ] = '2';
    kname[ 17*lenkname + 2 ] = 'H';
    kname[ 17*lenkname + 3 ] = '5';
    kname[ 17*lenkname + 4 ] = ' ';

    /* C2H6  */
    kname[ 18*lenkname + 0 ] = 'C';
    kname[ 18*lenkname + 1 ] = '2';
    kname[ 18*lenkname + 2 ] = 'H';
    kname[ 18*lenkname + 3 ] = '6';
    kname[ 18*lenkname + 4 ] = ' ';

    /* N2  */
    kname[ 19*lenkname + 0 ] = 'N';
    kname[ 19*lenkname + 1 ] = '2';
    kname[ 19*lenkname + 2 ] = ' ';

    /* AR  */
    kname[ 20*lenkname + 0 ] = 'A';
    kname[ 20*lenkname + 1 ] = 'R';
    kname[ 20*lenkname + 2 ] = ' ';

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
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*14.027090; /*CH2 */
    XW += x[8]*14.027090; /*CH2(S) */
    XW += x[9]*15.035060; /*CH3 */
    XW += x[10]*16.043030; /*CH4 */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*30.026490; /*CH2O */
    XW += x[15]*31.034460; /*CH3O */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*29.062150; /*C2H5 */
    XW += x[18]*30.070120; /*C2H6 */
    XW += x[19]*28.013400; /*N2 */
    XW += x[20]*39.948000; /*AR */
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
    YOW += y[7]/14.027090; /*CH2 */
    YOW += y[8]/14.027090; /*CH2(S) */
    YOW += y[9]/15.035060; /*CH3 */
    YOW += y[10]/16.043030; /*CH4 */
    YOW += y[11]/28.010550; /*CO */
    YOW += y[12]/44.009950; /*CO2 */
    YOW += y[13]/29.018520; /*HCO */
    YOW += y[14]/30.026490; /*CH2O */
    YOW += y[15]/31.034460; /*CH3O */
    YOW += y[16]/28.054180; /*C2H4 */
    YOW += y[17]/29.062150; /*C2H5 */
    YOW += y[18]/30.070120; /*C2H6 */
    YOW += y[19]/28.013400; /*N2 */
    YOW += y[20]/39.948000; /*AR */
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
    W += c[7]*14.027090; /*CH2 */
    W += c[8]*14.027090; /*CH2(S) */
    W += c[9]*15.035060; /*CH3 */
    W += c[10]*16.043030; /*CH4 */
    W += c[11]*28.010550; /*CO */
    W += c[12]*44.009950; /*CO2 */
    W += c[13]*29.018520; /*HCO */
    W += c[14]*30.026490; /*CH2O */
    W += c[15]*31.034460; /*CH3O */
    W += c[16]*28.054180; /*C2H4 */
    W += c[17]*29.062150; /*C2H5 */
    W += c[18]*30.070120; /*C2H6 */
    W += c[19]*28.013400; /*N2 */
    W += c[20]*39.948000; /*AR */

    for (id = 0; id < 21; ++id) {
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
    XW += x[7]*14.027090; /*CH2 */
    XW += x[8]*14.027090; /*CH2(S) */
    XW += x[9]*15.035060; /*CH3 */
    XW += x[10]*16.043030; /*CH4 */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*30.026490; /*CH2O */
    XW += x[15]*31.034460; /*CH3O */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*29.062150; /*C2H5 */
    XW += x[18]*30.070120; /*C2H6 */
    XW += x[19]*28.013400; /*N2 */
    XW += x[20]*39.948000; /*AR */
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
    YOW += y[7]/14.027090; /*CH2 */
    YOW += y[8]/14.027090; /*CH2(S) */
    YOW += y[9]/15.035060; /*CH3 */
    YOW += y[10]/16.043030; /*CH4 */
    YOW += y[11]/28.010550; /*CO */
    YOW += y[12]/44.009950; /*CO2 */
    YOW += y[13]/29.018520; /*HCO */
    YOW += y[14]/30.026490; /*CH2O */
    YOW += y[15]/31.034460; /*CH3O */
    YOW += y[16]/28.054180; /*C2H4 */
    YOW += y[17]/29.062150; /*C2H5 */
    YOW += y[18]/30.070120; /*C2H6 */
    YOW += y[19]/28.013400; /*N2 */
    YOW += y[20]/39.948000; /*AR */
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
    W += c[7]*14.027090; /*CH2 */
    W += c[8]*14.027090; /*CH2(S) */
    W += c[9]*15.035060; /*CH3 */
    W += c[10]*16.043030; /*CH4 */
    W += c[11]*28.010550; /*CO */
    W += c[12]*44.009950; /*CO2 */
    W += c[13]*29.018520; /*HCO */
    W += c[14]*30.026490; /*CH2O */
    W += c[15]*31.034460; /*CH3O */
    W += c[16]*28.054180; /*C2H4 */
    W += c[17]*29.062150; /*C2H5 */
    W += c[18]*30.070120; /*C2H6 */
    W += c[19]*28.013400; /*N2 */
    W += c[20]*39.948000; /*AR */

    for (id = 0; id < 21; ++id) {
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
    YOW += y[7]/14.027090; /*CH2 */
    YOW += y[8]/14.027090; /*CH2(S) */
    YOW += y[9]/15.035060; /*CH3 */
    YOW += y[10]/16.043030; /*CH4 */
    YOW += y[11]/28.010550; /*CO */
    YOW += y[12]/44.009950; /*CO2 */
    YOW += y[13]/29.018520; /*HCO */
    YOW += y[14]/30.026490; /*CH2O */
    YOW += y[15]/31.034460; /*CH3O */
    YOW += y[16]/28.054180; /*C2H4 */
    YOW += y[17]/29.062150; /*C2H5 */
    YOW += y[18]/30.070120; /*C2H6 */
    YOW += y[19]/28.013400; /*N2 */
    YOW += y[20]/39.948000; /*AR */
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
    XW += x[7]*14.027090; /*CH2 */
    XW += x[8]*14.027090; /*CH2(S) */
    XW += x[9]*15.035060; /*CH3 */
    XW += x[10]*16.043030; /*CH4 */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*30.026490; /*CH2O */
    XW += x[15]*31.034460; /*CH3O */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*29.062150; /*C2H5 */
    XW += x[18]*30.070120; /*C2H6 */
    XW += x[19]*28.013400; /*N2 */
    XW += x[20]*39.948000; /*AR */
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
    W += c[7]*14.027090; /*CH2 */
    W += c[8]*14.027090; /*CH2(S) */
    W += c[9]*15.035060; /*CH3 */
    W += c[10]*16.043030; /*CH4 */
    W += c[11]*28.010550; /*CO */
    W += c[12]*44.009950; /*CO2 */
    W += c[13]*29.018520; /*HCO */
    W += c[14]*30.026490; /*CH2O */
    W += c[15]*31.034460; /*CH3O */
    W += c[16]*28.054180; /*C2H4 */
    W += c[17]*29.062150; /*C2H5 */
    W += c[18]*30.070120; /*C2H6 */
    W += c[19]*28.013400; /*N2 */
    W += c[20]*39.948000; /*AR */

    for (id = 0; id < 21; ++id) {
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
    YOW += y[7]/14.027090; /*CH2 */
    YOW += y[8]/14.027090; /*CH2(S) */
    YOW += y[9]/15.035060; /*CH3 */
    YOW += y[10]/16.043030; /*CH4 */
    YOW += y[11]/28.010550; /*CO */
    YOW += y[12]/44.009950; /*CO2 */
    YOW += y[13]/29.018520; /*HCO */
    YOW += y[14]/30.026490; /*CH2O */
    YOW += y[15]/31.034460; /*CH3O */
    YOW += y[16]/28.054180; /*C2H4 */
    YOW += y[17]/29.062150; /*C2H5 */
    YOW += y[18]/30.070120; /*C2H6 */
    YOW += y[19]/28.013400; /*N2 */
    YOW += y[20]/39.948000; /*AR */
    /*Now compute conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(14.027090*YOW); 
    x[8] = y[8]/(14.027090*YOW); 
    x[9] = y[9]/(15.035060*YOW); 
    x[10] = y[10]/(16.043030*YOW); 
    x[11] = y[11]/(28.010550*YOW); 
    x[12] = y[12]/(44.009950*YOW); 
    x[13] = y[13]/(29.018520*YOW); 
    x[14] = y[14]/(30.026490*YOW); 
    x[15] = y[15]/(31.034460*YOW); 
    x[16] = y[16]/(28.054180*YOW); 
    x[17] = y[17]/(29.062150*YOW); 
    x[18] = y[18]/(30.070120*YOW); 
    x[19] = y[19]/(28.013400*YOW); 
    x[20] = y[20]/(39.948000*YOW); 

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
    YOW += y[7]/14.027090; /*CH2 */
    YOW += y[8]/14.027090; /*CH2(S) */
    YOW += y[9]/15.035060; /*CH3 */
    YOW += y[10]/16.043030; /*CH4 */
    YOW += y[11]/28.010550; /*CO */
    YOW += y[12]/44.009950; /*CO2 */
    YOW += y[13]/29.018520; /*HCO */
    YOW += y[14]/30.026490; /*CH2O */
    YOW += y[15]/31.034460; /*CH3O */
    YOW += y[16]/28.054180; /*C2H4 */
    YOW += y[17]/29.062150; /*C2H5 */
    YOW += y[18]/30.070120; /*C2H6 */
    YOW += y[19]/28.013400; /*N2 */
    YOW += y[20]/39.948000; /*AR */
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
    c[7] = PWORT * y[7]/14.027090; 
    c[8] = PWORT * y[8]/14.027090; 
    c[9] = PWORT * y[9]/15.035060; 
    c[10] = PWORT * y[10]/16.043030; 
    c[11] = PWORT * y[11]/28.010550; 
    c[12] = PWORT * y[12]/44.009950; 
    c[13] = PWORT * y[13]/29.018520; 
    c[14] = PWORT * y[14]/30.026490; 
    c[15] = PWORT * y[15]/31.034460; 
    c[16] = PWORT * y[16]/28.054180; 
    c[17] = PWORT * y[17]/29.062150; 
    c[18] = PWORT * y[18]/30.070120; 
    c[19] = PWORT * y[19]/28.013400; 
    c[20] = PWORT * y[20]/39.948000; 

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
    c[7] = (*rho) * y[7]/14.027090; 
    c[8] = (*rho) * y[8]/14.027090; 
    c[9] = (*rho) * y[9]/15.035060; 
    c[10] = (*rho) * y[10]/16.043030; 
    c[11] = (*rho) * y[11]/28.010550; 
    c[12] = (*rho) * y[12]/44.009950; 
    c[13] = (*rho) * y[13]/29.018520; 
    c[14] = (*rho) * y[14]/30.026490; 
    c[15] = (*rho) * y[15]/31.034460; 
    c[16] = (*rho) * y[16]/28.054180; 
    c[17] = (*rho) * y[17]/29.062150; 
    c[18] = (*rho) * y[18]/30.070120; 
    c[19] = (*rho) * y[19]/28.013400; 
    c[20] = (*rho) * y[20]/39.948000; 

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
    XW += x[7]*14.027090; /*CH2 */
    XW += x[8]*14.027090; /*CH2(S) */
    XW += x[9]*15.035060; /*CH3 */
    XW += x[10]*16.043030; /*CH4 */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*30.026490; /*CH2O */
    XW += x[15]*31.034460; /*CH3O */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*29.062150; /*C2H5 */
    XW += x[18]*30.070120; /*C2H6 */
    XW += x[19]*28.013400; /*N2 */
    XW += x[20]*39.948000; /*AR */
    /*Now compute conversion */
    y[0] = x[0]*2.015940/XW; 
    y[1] = x[1]*1.007970/XW; 
    y[2] = x[2]*15.999400/XW; 
    y[3] = x[3]*31.998800/XW; 
    y[4] = x[4]*17.007370/XW; 
    y[5] = x[5]*18.015340/XW; 
    y[6] = x[6]*33.006770/XW; 
    y[7] = x[7]*14.027090/XW; 
    y[8] = x[8]*14.027090/XW; 
    y[9] = x[9]*15.035060/XW; 
    y[10] = x[10]*16.043030/XW; 
    y[11] = x[11]*28.010550/XW; 
    y[12] = x[12]*44.009950/XW; 
    y[13] = x[13]*29.018520/XW; 
    y[14] = x[14]*30.026490/XW; 
    y[15] = x[15]*31.034460/XW; 
    y[16] = x[16]*28.054180/XW; 
    y[17] = x[17]*29.062150/XW; 
    y[18] = x[18]*30.070120/XW; 
    y[19] = x[19]*28.013400/XW; 
    y[20] = x[20]*39.948000/XW; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 21; ++id) {
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
    XW += x[7]*14.027090; /*CH2 */
    XW += x[8]*14.027090; /*CH2(S) */
    XW += x[9]*15.035060; /*CH3 */
    XW += x[10]*16.043030; /*CH4 */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*30.026490; /*CH2O */
    XW += x[15]*31.034460; /*CH3O */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*29.062150; /*C2H5 */
    XW += x[18]*30.070120; /*C2H6 */
    XW += x[19]*28.013400; /*N2 */
    XW += x[20]*39.948000; /*AR */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 21; ++id) {
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
    for (id = 0; id < 21; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    for (id = 0; id < 21; ++id) {
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
    CW += c[7]*14.027090; /*CH2 */
    CW += c[8]*14.027090; /*CH2(S) */
    CW += c[9]*15.035060; /*CH3 */
    CW += c[10]*16.043030; /*CH4 */
    CW += c[11]*28.010550; /*CO */
    CW += c[12]*44.009950; /*CO2 */
    CW += c[13]*29.018520; /*HCO */
    CW += c[14]*30.026490; /*CH2O */
    CW += c[15]*31.034460; /*CH3O */
    CW += c[16]*28.054180; /*C2H4 */
    CW += c[17]*29.062150; /*C2H5 */
    CW += c[18]*30.070120; /*C2H6 */
    CW += c[19]*28.013400; /*N2 */
    CW += c[20]*39.948000; /*AR */
    /*Now compute conversion */
    y[0] = c[0]*2.015940/CW; 
    y[1] = c[1]*1.007970/CW; 
    y[2] = c[2]*15.999400/CW; 
    y[3] = c[3]*31.998800/CW; 
    y[4] = c[4]*17.007370/CW; 
    y[5] = c[5]*18.015340/CW; 
    y[6] = c[6]*33.006770/CW; 
    y[7] = c[7]*14.027090/CW; 
    y[8] = c[8]*14.027090/CW; 
    y[9] = c[9]*15.035060/CW; 
    y[10] = c[10]*16.043030/CW; 
    y[11] = c[11]*28.010550/CW; 
    y[12] = c[12]*44.009950/CW; 
    y[13] = c[13]*29.018520/CW; 
    y[14] = c[14]*30.026490/CW; 
    y[15] = c[15]*31.034460/CW; 
    y[16] = c[16]*28.054180/CW; 
    y[17] = c[17]*29.062150/CW; 
    y[18] = c[18]*30.070120/CW; 
    y[19] = c[19]*28.013400/CW; 
    y[20] = c[20]*39.948000/CW; 

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
    for (id = 0; id < 21; ++id) {
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
    for (id = 0; id < 21; ++id) {
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
    for (id = 0; id < 21; ++id) {
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
    for (id = 0; id < 21; ++id) {
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
    for (id = 0; id < 21; ++id) {
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
    for (id = 0; id < 21; ++id) {
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
    for (id = 0; id < 21; ++id) {
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
    cvms[0] *= 4.124383662212169e+07; /*H2 */
    cvms[1] *= 8.248767324424338e+07; /*H */
    cvms[2] *= 5.196763628636074e+06; /*O */
    cvms[3] *= 2.598381814318037e+06; /*O2 */
    cvms[4] *= 4.888768810227566e+06; /*OH */
    cvms[5] *= 4.615239012974499e+06; /*H2O */
    cvms[6] *= 2.519031701678171e+06; /*HO2 */
    cvms[7] *= 5.927466067445207e+06; /*CH2 */
    cvms[8] *= 5.927466067445207e+06; /*CH2(S) */
    cvms[9] *= 5.530081023953346e+06; /*CH3 */
    cvms[10] *= 5.182630712527496e+06; /*CH4 */
    cvms[11] *= 2.968349425484326e+06; /*CO */
    cvms[12] *= 1.889234139098090e+06; /*CO2 */
    cvms[13] *= 2.865242610581105e+06; /*HCO */
    cvms[14] *= 2.769058254894261e+06; /*CH2O */
    cvms[15] *= 2.679121853578248e+06; /*CH3O */
    cvms[16] *= 2.963733033722604e+06; /*C2H4 */
    cvms[17] *= 2.860941121011349e+06; /*C2H5 */
    cvms[18] *= 2.765040511976673e+06; /*C2H6 */
    cvms[19] *= 2.968047434442088e+06; /*N2 */
    cvms[20] *= 2.081333233203164e+06; /*AR */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
void CKCPMS(double *T, int * iwrk, double * rwrk, double * cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 4.124383662212169e+07; /*H2 */
    cpms[1] *= 8.248767324424338e+07; /*H */
    cpms[2] *= 5.196763628636074e+06; /*O */
    cpms[3] *= 2.598381814318037e+06; /*O2 */
    cpms[4] *= 4.888768810227566e+06; /*OH */
    cpms[5] *= 4.615239012974499e+06; /*H2O */
    cpms[6] *= 2.519031701678171e+06; /*HO2 */
    cpms[7] *= 5.927466067445207e+06; /*CH2 */
    cpms[8] *= 5.927466067445207e+06; /*CH2(S) */
    cpms[9] *= 5.530081023953346e+06; /*CH3 */
    cpms[10] *= 5.182630712527496e+06; /*CH4 */
    cpms[11] *= 2.968349425484326e+06; /*CO */
    cpms[12] *= 1.889234139098090e+06; /*CO2 */
    cpms[13] *= 2.865242610581105e+06; /*HCO */
    cpms[14] *= 2.769058254894261e+06; /*CH2O */
    cpms[15] *= 2.679121853578248e+06; /*CH3O */
    cpms[16] *= 2.963733033722604e+06; /*C2H4 */
    cpms[17] *= 2.860941121011349e+06; /*C2H5 */
    cpms[18] *= 2.765040511976673e+06; /*C2H6 */
    cpms[19] *= 2.968047434442088e+06; /*N2 */
    cpms[20] *= 2.081333233203164e+06; /*AR */
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
    ums[7] *= RT/14.027090; /*CH2 */
    ums[8] *= RT/14.027090; /*CH2(S) */
    ums[9] *= RT/15.035060; /*CH3 */
    ums[10] *= RT/16.043030; /*CH4 */
    ums[11] *= RT/28.010550; /*CO */
    ums[12] *= RT/44.009950; /*CO2 */
    ums[13] *= RT/29.018520; /*HCO */
    ums[14] *= RT/30.026490; /*CH2O */
    ums[15] *= RT/31.034460; /*CH3O */
    ums[16] *= RT/28.054180; /*C2H4 */
    ums[17] *= RT/29.062150; /*C2H5 */
    ums[18] *= RT/30.070120; /*C2H6 */
    ums[19] *= RT/28.013400; /*N2 */
    ums[20] *= RT/39.948000; /*AR */
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
    hms[7] *= RT/14.027090; /*CH2 */
    hms[8] *= RT/14.027090; /*CH2(S) */
    hms[9] *= RT/15.035060; /*CH3 */
    hms[10] *= RT/16.043030; /*CH4 */
    hms[11] *= RT/28.010550; /*CO */
    hms[12] *= RT/44.009950; /*CO2 */
    hms[13] *= RT/29.018520; /*HCO */
    hms[14] *= RT/30.026490; /*CH2O */
    hms[15] *= RT/31.034460; /*CH3O */
    hms[16] *= RT/28.054180; /*C2H4 */
    hms[17] *= RT/29.062150; /*C2H5 */
    hms[18] *= RT/30.070120; /*C2H6 */
    hms[19] *= RT/28.013400; /*N2 */
    hms[20] *= RT/39.948000; /*AR */
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
    gms[7] *= RT/14.027090; /*CH2 */
    gms[8] *= RT/14.027090; /*CH2(S) */
    gms[9] *= RT/15.035060; /*CH3 */
    gms[10] *= RT/16.043030; /*CH4 */
    gms[11] *= RT/28.010550; /*CO */
    gms[12] *= RT/44.009950; /*CO2 */
    gms[13] *= RT/29.018520; /*HCO */
    gms[14] *= RT/30.026490; /*CH2O */
    gms[15] *= RT/31.034460; /*CH3O */
    gms[16] *= RT/28.054180; /*C2H4 */
    gms[17] *= RT/29.062150; /*C2H5 */
    gms[18] *= RT/30.070120; /*C2H6 */
    gms[19] *= RT/28.013400; /*N2 */
    gms[20] *= RT/39.948000; /*AR */
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
    ams[7] *= RT/14.027090; /*CH2 */
    ams[8] *= RT/14.027090; /*CH2(S) */
    ams[9] *= RT/15.035060; /*CH3 */
    ams[10] *= RT/16.043030; /*CH4 */
    ams[11] *= RT/28.010550; /*CO */
    ams[12] *= RT/44.009950; /*CO2 */
    ams[13] *= RT/29.018520; /*HCO */
    ams[14] *= RT/30.026490; /*CH2O */
    ams[15] *= RT/31.034460; /*CH3O */
    ams[16] *= RT/28.054180; /*C2H4 */
    ams[17] *= RT/29.062150; /*C2H5 */
    ams[18] *= RT/30.070120; /*C2H6 */
    ams[19] *= RT/28.013400; /*N2 */
    ams[20] *= RT/39.948000; /*AR */
}


/*Returns the entropies in mass units (Eq 28.) */
void CKSMS(double *T, int * iwrk, double * rwrk, double * sms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 4.124383662212169e+07; /*H2 */
    sms[1] *= 8.248767324424338e+07; /*H */
    sms[2] *= 5.196763628636074e+06; /*O */
    sms[3] *= 2.598381814318037e+06; /*O2 */
    sms[4] *= 4.888768810227566e+06; /*OH */
    sms[5] *= 4.615239012974499e+06; /*H2O */
    sms[6] *= 2.519031701678171e+06; /*HO2 */
    sms[7] *= 5.927466067445207e+06; /*CH2 */
    sms[8] *= 5.927466067445207e+06; /*CH2(S) */
    sms[9] *= 5.530081023953346e+06; /*CH3 */
    sms[10] *= 5.182630712527496e+06; /*CH4 */
    sms[11] *= 2.968349425484326e+06; /*CO */
    sms[12] *= 1.889234139098090e+06; /*CO2 */
    sms[13] *= 2.865242610581105e+06; /*HCO */
    sms[14] *= 2.769058254894261e+06; /*CH2O */
    sms[15] *= 2.679121853578248e+06; /*CH3O */
    sms[16] *= 2.963733033722604e+06; /*C2H4 */
    sms[17] *= 2.860941121011349e+06; /*C2H5 */
    sms[18] *= 2.765040511976673e+06; /*C2H6 */
    sms[19] *= 2.968047434442088e+06; /*N2 */
    sms[20] *= 2.081333233203164e+06; /*AR */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double *T, double *x, int * iwrk, double * rwrk, double * cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[21]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 21; ++id) {
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
    double cpor[21]; /* temporary storage */
    cp_R(cpor, tc);
    /*multiply by y/molecularweight */
    result += cpor[0]*y[0]/2.015940; /*H2 */
    result += cpor[1]*y[1]/1.007970; /*H */
    result += cpor[2]*y[2]/15.999400; /*O */
    result += cpor[3]*y[3]/31.998800; /*O2 */
    result += cpor[4]*y[4]/17.007370; /*OH */
    result += cpor[5]*y[5]/18.015340; /*H2O */
    result += cpor[6]*y[6]/33.006770; /*HO2 */
    result += cpor[7]*y[7]/14.027090; /*CH2 */
    result += cpor[8]*y[8]/14.027090; /*CH2(S) */
    result += cpor[9]*y[9]/15.035060; /*CH3 */
    result += cpor[10]*y[10]/16.043030; /*CH4 */
    result += cpor[11]*y[11]/28.010550; /*CO */
    result += cpor[12]*y[12]/44.009950; /*CO2 */
    result += cpor[13]*y[13]/29.018520; /*HCO */
    result += cpor[14]*y[14]/30.026490; /*CH2O */
    result += cpor[15]*y[15]/31.034460; /*CH3O */
    result += cpor[16]*y[16]/28.054180; /*C2H4 */
    result += cpor[17]*y[17]/29.062150; /*C2H5 */
    result += cpor[18]*y[18]/30.070120; /*C2H6 */
    result += cpor[19]*y[19]/28.013400; /*N2 */
    result += cpor[20]*y[20]/39.948000; /*AR */

    *cpbs = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double *T, double *x, int * iwrk, double * rwrk, double * cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[21]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 21; ++id) {
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
    double cvor[21]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]/2.015940; /*H2 */
    result += cvor[1]*y[1]/1.007970; /*H */
    result += cvor[2]*y[2]/15.999400; /*O */
    result += cvor[3]*y[3]/31.998800; /*O2 */
    result += cvor[4]*y[4]/17.007370; /*OH */
    result += cvor[5]*y[5]/18.015340; /*H2O */
    result += cvor[6]*y[6]/33.006770; /*HO2 */
    result += cvor[7]*y[7]/14.027090; /*CH2 */
    result += cvor[8]*y[8]/14.027090; /*CH2(S) */
    result += cvor[9]*y[9]/15.035060; /*CH3 */
    result += cvor[10]*y[10]/16.043030; /*CH4 */
    result += cvor[11]*y[11]/28.010550; /*CO */
    result += cvor[12]*y[12]/44.009950; /*CO2 */
    result += cvor[13]*y[13]/29.018520; /*HCO */
    result += cvor[14]*y[14]/30.026490; /*CH2O */
    result += cvor[15]*y[15]/31.034460; /*CH3O */
    result += cvor[16]*y[16]/28.054180; /*C2H4 */
    result += cvor[17]*y[17]/29.062150; /*C2H5 */
    result += cvor[18]*y[18]/30.070120; /*C2H6 */
    result += cvor[19]*y[19]/28.013400; /*N2 */
    result += cvor[20]*y[20]/39.948000; /*AR */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double *T, double *x, int * iwrk, double * rwrk, double * hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[21]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 21; ++id) {
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
    double hml[21]; /* temporary storage */
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
    result += y[7]*hml[7]/14.027090; /*CH2 */
    result += y[8]*hml[8]/14.027090; /*CH2(S) */
    result += y[9]*hml[9]/15.035060; /*CH3 */
    result += y[10]*hml[10]/16.043030; /*CH4 */
    result += y[11]*hml[11]/28.010550; /*CO */
    result += y[12]*hml[12]/44.009950; /*CO2 */
    result += y[13]*hml[13]/29.018520; /*HCO */
    result += y[14]*hml[14]/30.026490; /*CH2O */
    result += y[15]*hml[15]/31.034460; /*CH3O */
    result += y[16]*hml[16]/28.054180; /*C2H4 */
    result += y[17]*hml[17]/29.062150; /*C2H5 */
    result += y[18]*hml[18]/30.070120; /*C2H6 */
    result += y[19]*hml[19]/28.013400; /*N2 */
    result += y[20]*hml[20]/39.948000; /*AR */

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double *T, double *x, int * iwrk, double * rwrk, double * ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[21]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 21; ++id) {
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
    double ums[21]; /* temporary energy array */
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
    result += y[7]*ums[7]/14.027090; /*CH2 */
    result += y[8]*ums[8]/14.027090; /*CH2(S) */
    result += y[9]*ums[9]/15.035060; /*CH3 */
    result += y[10]*ums[10]/16.043030; /*CH4 */
    result += y[11]*ums[11]/28.010550; /*CO */
    result += y[12]*ums[12]/44.009950; /*CO2 */
    result += y[13]*ums[13]/29.018520; /*HCO */
    result += y[14]*ums[14]/30.026490; /*CH2O */
    result += y[15]*ums[15]/31.034460; /*CH3O */
    result += y[16]*ums[16]/28.054180; /*C2H4 */
    result += y[17]*ums[17]/29.062150; /*C2H5 */
    result += y[18]*ums[18]/30.070120; /*C2H6 */
    result += y[19]*ums[19]/28.013400; /*N2 */
    result += y[20]*ums[20]/39.948000; /*AR */

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
    double sor[21]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 21; ++id) {
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
    double sor[21]; /* temporary storage */
    double x[21]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/2.015940; /*H2 */
    YOW += y[1]/1.007970; /*H */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/18.015340; /*H2O */
    YOW += y[6]/33.006770; /*HO2 */
    YOW += y[7]/14.027090; /*CH2 */
    YOW += y[8]/14.027090; /*CH2(S) */
    YOW += y[9]/15.035060; /*CH3 */
    YOW += y[10]/16.043030; /*CH4 */
    YOW += y[11]/28.010550; /*CO */
    YOW += y[12]/44.009950; /*CO2 */
    YOW += y[13]/29.018520; /*HCO */
    YOW += y[14]/30.026490; /*CH2O */
    YOW += y[15]/31.034460; /*CH3O */
    YOW += y[16]/28.054180; /*C2H4 */
    YOW += y[17]/29.062150; /*C2H5 */
    YOW += y[18]/30.070120; /*C2H6 */
    YOW += y[19]/28.013400; /*N2 */
    YOW += y[20]/39.948000; /*AR */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(14.027090*YOW); 
    x[8] = y[8]/(14.027090*YOW); 
    x[9] = y[9]/(15.035060*YOW); 
    x[10] = y[10]/(16.043030*YOW); 
    x[11] = y[11]/(28.010550*YOW); 
    x[12] = y[12]/(44.009950*YOW); 
    x[13] = y[13]/(29.018520*YOW); 
    x[14] = y[14]/(30.026490*YOW); 
    x[15] = y[15]/(31.034460*YOW); 
    x[16] = y[16]/(28.054180*YOW); 
    x[17] = y[17]/(29.062150*YOW); 
    x[18] = y[18]/(30.070120*YOW); 
    x[19] = y[19]/(28.013400*YOW); 
    x[20] = y[20]/(39.948000*YOW); 
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
    double gort[21]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 21; ++id) {
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
    double gort[21]; /* temporary storage */
    double x[21]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/2.015940; /*H2 */
    YOW += y[1]/1.007970; /*H */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/18.015340; /*H2O */
    YOW += y[6]/33.006770; /*HO2 */
    YOW += y[7]/14.027090; /*CH2 */
    YOW += y[8]/14.027090; /*CH2(S) */
    YOW += y[9]/15.035060; /*CH3 */
    YOW += y[10]/16.043030; /*CH4 */
    YOW += y[11]/28.010550; /*CO */
    YOW += y[12]/44.009950; /*CO2 */
    YOW += y[13]/29.018520; /*HCO */
    YOW += y[14]/30.026490; /*CH2O */
    YOW += y[15]/31.034460; /*CH3O */
    YOW += y[16]/28.054180; /*C2H4 */
    YOW += y[17]/29.062150; /*C2H5 */
    YOW += y[18]/30.070120; /*C2H6 */
    YOW += y[19]/28.013400; /*N2 */
    YOW += y[20]/39.948000; /*AR */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(14.027090*YOW); 
    x[8] = y[8]/(14.027090*YOW); 
    x[9] = y[9]/(15.035060*YOW); 
    x[10] = y[10]/(16.043030*YOW); 
    x[11] = y[11]/(28.010550*YOW); 
    x[12] = y[12]/(44.009950*YOW); 
    x[13] = y[13]/(29.018520*YOW); 
    x[14] = y[14]/(30.026490*YOW); 
    x[15] = y[15]/(31.034460*YOW); 
    x[16] = y[16]/(28.054180*YOW); 
    x[17] = y[17]/(29.062150*YOW); 
    x[18] = y[18]/(30.070120*YOW); 
    x[19] = y[19]/(28.013400*YOW); 
    x[20] = y[20]/(39.948000*YOW); 
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
    double aort[21]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 21; ++id) {
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
    double aort[21]; /* temporary storage */
    double x[21]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/2.015940; /*H2 */
    YOW += y[1]/1.007970; /*H */
    YOW += y[2]/15.999400; /*O */
    YOW += y[3]/31.998800; /*O2 */
    YOW += y[4]/17.007370; /*OH */
    YOW += y[5]/18.015340; /*H2O */
    YOW += y[6]/33.006770; /*HO2 */
    YOW += y[7]/14.027090; /*CH2 */
    YOW += y[8]/14.027090; /*CH2(S) */
    YOW += y[9]/15.035060; /*CH3 */
    YOW += y[10]/16.043030; /*CH4 */
    YOW += y[11]/28.010550; /*CO */
    YOW += y[12]/44.009950; /*CO2 */
    YOW += y[13]/29.018520; /*HCO */
    YOW += y[14]/30.026490; /*CH2O */
    YOW += y[15]/31.034460; /*CH3O */
    YOW += y[16]/28.054180; /*C2H4 */
    YOW += y[17]/29.062150; /*C2H5 */
    YOW += y[18]/30.070120; /*C2H6 */
    YOW += y[19]/28.013400; /*N2 */
    YOW += y[20]/39.948000; /*AR */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(14.027090*YOW); 
    x[8] = y[8]/(14.027090*YOW); 
    x[9] = y[9]/(15.035060*YOW); 
    x[10] = y[10]/(16.043030*YOW); 
    x[11] = y[11]/(28.010550*YOW); 
    x[12] = y[12]/(44.009950*YOW); 
    x[13] = y[13]/(29.018520*YOW); 
    x[14] = y[14]/(30.026490*YOW); 
    x[15] = y[15]/(31.034460*YOW); 
    x[16] = y[16]/(28.054180*YOW); 
    x[17] = y[17]/(29.062150*YOW); 
    x[18] = y[18]/(30.070120*YOW); 
    x[19] = y[19]/(28.013400*YOW); 
    x[20] = y[20]/(39.948000*YOW); 
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
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double * T, double * C, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 21; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[21]; /*temporary storage */
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
    YOW += y[7]/14.027090; /*CH2 */
    YOW += y[8]/14.027090; /*CH2(S) */
    YOW += y[9]/15.035060; /*CH3 */
    YOW += y[10]/16.043030; /*CH4 */
    YOW += y[11]/28.010550; /*CO */
    YOW += y[12]/44.009950; /*CO2 */
    YOW += y[13]/29.018520; /*HCO */
    YOW += y[14]/30.026490; /*CH2O */
    YOW += y[15]/31.034460; /*CH3O */
    YOW += y[16]/28.054180; /*C2H4 */
    YOW += y[17]/29.062150; /*C2H5 */
    YOW += y[18]/30.070120; /*C2H6 */
    YOW += y[19]/28.013400; /*N2 */
    YOW += y[20]/39.948000; /*AR */
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
    c[7] = PWORT * y[7]/14.027090; 
    c[8] = PWORT * y[8]/14.027090; 
    c[9] = PWORT * y[9]/15.035060; 
    c[10] = PWORT * y[10]/16.043030; 
    c[11] = PWORT * y[11]/28.010550; 
    c[12] = PWORT * y[12]/44.009950; 
    c[13] = PWORT * y[13]/29.018520; 
    c[14] = PWORT * y[14]/30.026490; 
    c[15] = PWORT * y[15]/31.034460; 
    c[16] = PWORT * y[16]/28.054180; 
    c[17] = PWORT * y[17]/29.062150; 
    c[18] = PWORT * y[18]/30.070120; 
    c[19] = PWORT * y[19]/28.013400; 
    c[20] = PWORT * y[20]/39.948000; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[21]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 21; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[21]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]/2.015940; 
    c[1] = 1e6 * (*rho) * y[1]/1.007970; 
    c[2] = 1e6 * (*rho) * y[2]/15.999400; 
    c[3] = 1e6 * (*rho) * y[3]/31.998800; 
    c[4] = 1e6 * (*rho) * y[4]/17.007370; 
    c[5] = 1e6 * (*rho) * y[5]/18.015340; 
    c[6] = 1e6 * (*rho) * y[6]/33.006770; 
    c[7] = 1e6 * (*rho) * y[7]/14.027090; 
    c[8] = 1e6 * (*rho) * y[8]/14.027090; 
    c[9] = 1e6 * (*rho) * y[9]/15.035060; 
    c[10] = 1e6 * (*rho) * y[10]/16.043030; 
    c[11] = 1e6 * (*rho) * y[11]/28.010550; 
    c[12] = 1e6 * (*rho) * y[12]/44.009950; 
    c[13] = 1e6 * (*rho) * y[13]/29.018520; 
    c[14] = 1e6 * (*rho) * y[14]/30.026490; 
    c[15] = 1e6 * (*rho) * y[15]/31.034460; 
    c[16] = 1e6 * (*rho) * y[16]/28.054180; 
    c[17] = 1e6 * (*rho) * y[17]/29.062150; 
    c[18] = 1e6 * (*rho) * y[18]/30.070120; 
    c[19] = 1e6 * (*rho) * y[19]/28.013400; 
    c[20] = 1e6 * (*rho) * y[20]/39.948000; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[21]; /*temporary storage */
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
    XW += x[7]*14.027090; /*CH2 */
    XW += x[8]*14.027090; /*CH2(S) */
    XW += x[9]*15.035060; /*CH3 */
    XW += x[10]*16.043030; /*CH4 */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*30.026490; /*CH2O */
    XW += x[15]*31.034460; /*CH3O */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*29.062150; /*C2H5 */
    XW += x[18]*30.070120; /*C2H6 */
    XW += x[19]*28.013400; /*N2 */
    XW += x[20]*39.948000; /*AR */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 21; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double * T, double * C, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 21; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 84; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double * P, double * T, double * x, int * iwrk, double * rwrk, double * q_f, double * q_r)
{
    int id; /*loop counter */
    double c[21]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 21; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 84; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[21]; /*temporary storage */
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
    YOW += y[7]/14.027090; /*CH2 */
    YOW += y[8]/14.027090; /*CH2(S) */
    YOW += y[9]/15.035060; /*CH3 */
    YOW += y[10]/16.043030; /*CH4 */
    YOW += y[11]/28.010550; /*CO */
    YOW += y[12]/44.009950; /*CO2 */
    YOW += y[13]/29.018520; /*HCO */
    YOW += y[14]/30.026490; /*CH2O */
    YOW += y[15]/31.034460; /*CH3O */
    YOW += y[16]/28.054180; /*C2H4 */
    YOW += y[17]/29.062150; /*C2H5 */
    YOW += y[18]/30.070120; /*C2H6 */
    YOW += y[19]/28.013400; /*N2 */
    YOW += y[20]/39.948000; /*AR */
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
    c[7] = PWORT * y[7]/14.027090; 
    c[8] = PWORT * y[8]/14.027090; 
    c[9] = PWORT * y[9]/15.035060; 
    c[10] = PWORT * y[10]/16.043030; 
    c[11] = PWORT * y[11]/28.010550; 
    c[12] = PWORT * y[12]/44.009950; 
    c[13] = PWORT * y[13]/29.018520; 
    c[14] = PWORT * y[14]/30.026490; 
    c[15] = PWORT * y[15]/31.034460; 
    c[16] = PWORT * y[16]/28.054180; 
    c[17] = PWORT * y[17]/29.062150; 
    c[18] = PWORT * y[18]/30.070120; 
    c[19] = PWORT * y[19]/28.013400; 
    c[20] = PWORT * y[20]/39.948000; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 84; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[21]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 21; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 84; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[21]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]/2.015940; 
    c[1] = 1e6 * (*rho) * y[1]/1.007970; 
    c[2] = 1e6 * (*rho) * y[2]/15.999400; 
    c[3] = 1e6 * (*rho) * y[3]/31.998800; 
    c[4] = 1e6 * (*rho) * y[4]/17.007370; 
    c[5] = 1e6 * (*rho) * y[5]/18.015340; 
    c[6] = 1e6 * (*rho) * y[6]/33.006770; 
    c[7] = 1e6 * (*rho) * y[7]/14.027090; 
    c[8] = 1e6 * (*rho) * y[8]/14.027090; 
    c[9] = 1e6 * (*rho) * y[9]/15.035060; 
    c[10] = 1e6 * (*rho) * y[10]/16.043030; 
    c[11] = 1e6 * (*rho) * y[11]/28.010550; 
    c[12] = 1e6 * (*rho) * y[12]/44.009950; 
    c[13] = 1e6 * (*rho) * y[13]/29.018520; 
    c[14] = 1e6 * (*rho) * y[14]/30.026490; 
    c[15] = 1e6 * (*rho) * y[15]/31.034460; 
    c[16] = 1e6 * (*rho) * y[16]/28.054180; 
    c[17] = 1e6 * (*rho) * y[17]/29.062150; 
    c[18] = 1e6 * (*rho) * y[18]/30.070120; 
    c[19] = 1e6 * (*rho) * y[19]/28.013400; 
    c[20] = 1e6 * (*rho) * y[20]/39.948000; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 84; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[21]; /*temporary storage */
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
    XW += x[7]*14.027090; /*CH2 */
    XW += x[8]*14.027090; /*CH2(S) */
    XW += x[9]*15.035060; /*CH3 */
    XW += x[10]*16.043030; /*CH4 */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*30.026490; /*CH2O */
    XW += x[15]*31.034460; /*CH3O */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*29.062150; /*C2H5 */
    XW += x[18]*30.070120; /*C2H6 */
    XW += x[19]*28.013400; /*N2 */
    XW += x[20]*39.948000; /*AR */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 21; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 84; ++id) {
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
    for (id = 0; id < 21 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: O + H + M <=> OH + M */
    nuki[ 2 * kd + 0 ] += -1 ;
    nuki[ 1 * kd + 0 ] += -1 ;
    nuki[ 4 * kd + 0 ] += +1 ;

    /*reaction 2: O + H2 <=> H + OH */
    nuki[ 2 * kd + 1 ] += -1 ;
    nuki[ 0 * kd + 1 ] += -1 ;
    nuki[ 1 * kd + 1 ] += +1 ;
    nuki[ 4 * kd + 1 ] += +1 ;

    /*reaction 3: O + HO2 <=> OH + O2 */
    nuki[ 2 * kd + 2 ] += -1 ;
    nuki[ 6 * kd + 2 ] += -1 ;
    nuki[ 4 * kd + 2 ] += +1 ;
    nuki[ 3 * kd + 2 ] += +1 ;

    /*reaction 4: O + CH2 <=> H + HCO */
    nuki[ 2 * kd + 3 ] += -1 ;
    nuki[ 7 * kd + 3 ] += -1 ;
    nuki[ 1 * kd + 3 ] += +1 ;
    nuki[ 13 * kd + 3 ] += +1 ;

    /*reaction 5: O + CH2(S) <=> H + HCO */
    nuki[ 2 * kd + 4 ] += -1 ;
    nuki[ 8 * kd + 4 ] += -1 ;
    nuki[ 1 * kd + 4 ] += +1 ;
    nuki[ 13 * kd + 4 ] += +1 ;

    /*reaction 6: O + CH3 <=> H + CH2O */
    nuki[ 2 * kd + 5 ] += -1 ;
    nuki[ 9 * kd + 5 ] += -1 ;
    nuki[ 1 * kd + 5 ] += +1 ;
    nuki[ 14 * kd + 5 ] += +1 ;

    /*reaction 7: O + CH4 <=> OH + CH3 */
    nuki[ 2 * kd + 6 ] += -1 ;
    nuki[ 10 * kd + 6 ] += -1 ;
    nuki[ 4 * kd + 6 ] += +1 ;
    nuki[ 9 * kd + 6 ] += +1 ;

    /*reaction 8: O + CO + M <=> CO2 + M */
    nuki[ 2 * kd + 7 ] += -1 ;
    nuki[ 11 * kd + 7 ] += -1 ;
    nuki[ 12 * kd + 7 ] += +1 ;

    /*reaction 9: O + HCO <=> OH + CO */
    nuki[ 2 * kd + 8 ] += -1 ;
    nuki[ 13 * kd + 8 ] += -1 ;
    nuki[ 4 * kd + 8 ] += +1 ;
    nuki[ 11 * kd + 8 ] += +1 ;

    /*reaction 10: O + HCO <=> H + CO2 */
    nuki[ 2 * kd + 9 ] += -1 ;
    nuki[ 13 * kd + 9 ] += -1 ;
    nuki[ 1 * kd + 9 ] += +1 ;
    nuki[ 12 * kd + 9 ] += +1 ;

    /*reaction 11: O + CH2O <=> OH + HCO */
    nuki[ 2 * kd + 10 ] += -1 ;
    nuki[ 14 * kd + 10 ] += -1 ;
    nuki[ 4 * kd + 10 ] += +1 ;
    nuki[ 13 * kd + 10 ] += +1 ;

    /*reaction 12: O + C2H4 <=> CH3 + HCO */
    nuki[ 2 * kd + 11 ] += -1 ;
    nuki[ 16 * kd + 11 ] += -1 ;
    nuki[ 9 * kd + 11 ] += +1 ;
    nuki[ 13 * kd + 11 ] += +1 ;

    /*reaction 13: O + C2H5 <=> CH3 + CH2O */
    nuki[ 2 * kd + 12 ] += -1 ;
    nuki[ 17 * kd + 12 ] += -1 ;
    nuki[ 9 * kd + 12 ] += +1 ;
    nuki[ 14 * kd + 12 ] += +1 ;

    /*reaction 14: O + C2H6 <=> OH + C2H5 */
    nuki[ 2 * kd + 13 ] += -1 ;
    nuki[ 18 * kd + 13 ] += -1 ;
    nuki[ 4 * kd + 13 ] += +1 ;
    nuki[ 17 * kd + 13 ] += +1 ;

    /*reaction 15: O2 + CO <=> O + CO2 */
    nuki[ 3 * kd + 14 ] += -1 ;
    nuki[ 11 * kd + 14 ] += -1 ;
    nuki[ 2 * kd + 14 ] += +1 ;
    nuki[ 12 * kd + 14 ] += +1 ;

    /*reaction 16: O2 + CH2O <=> HO2 + HCO */
    nuki[ 3 * kd + 15 ] += -1 ;
    nuki[ 14 * kd + 15 ] += -1 ;
    nuki[ 6 * kd + 15 ] += +1 ;
    nuki[ 13 * kd + 15 ] += +1 ;

    /*reaction 17: H + O2 + M <=> HO2 + M */
    nuki[ 1 * kd + 16 ] += -1 ;
    nuki[ 3 * kd + 16 ] += -1 ;
    nuki[ 6 * kd + 16 ] += +1 ;

    /*reaction 18: H + 2 O2 <=> HO2 + O2 */
    nuki[ 1 * kd + 17 ] += -1 ;
    nuki[ 3 * kd + 17 ] += -2 ;
    nuki[ 6 * kd + 17 ] += +1 ;
    nuki[ 3 * kd + 17 ] += +1 ;

    /*reaction 19: H + O2 + H2O <=> HO2 + H2O */
    nuki[ 1 * kd + 18 ] += -1 ;
    nuki[ 3 * kd + 18 ] += -1 ;
    nuki[ 5 * kd + 18 ] += -1 ;
    nuki[ 6 * kd + 18 ] += +1 ;
    nuki[ 5 * kd + 18 ] += +1 ;

    /*reaction 20: H + O2 + N2 <=> HO2 + N2 */
    nuki[ 1 * kd + 19 ] += -1 ;
    nuki[ 3 * kd + 19 ] += -1 ;
    nuki[ 19 * kd + 19 ] += -1 ;
    nuki[ 6 * kd + 19 ] += +1 ;
    nuki[ 19 * kd + 19 ] += +1 ;

    /*reaction 21: H + O2 + AR <=> HO2 + AR */
    nuki[ 1 * kd + 20 ] += -1 ;
    nuki[ 3 * kd + 20 ] += -1 ;
    nuki[ 20 * kd + 20 ] += -1 ;
    nuki[ 6 * kd + 20 ] += +1 ;
    nuki[ 20 * kd + 20 ] += +1 ;

    /*reaction 22: H + O2 <=> O + OH */
    nuki[ 1 * kd + 21 ] += -1 ;
    nuki[ 3 * kd + 21 ] += -1 ;
    nuki[ 2 * kd + 21 ] += +1 ;
    nuki[ 4 * kd + 21 ] += +1 ;

    /*reaction 23: 2 H + M <=> H2 + M */
    nuki[ 1 * kd + 22 ] += -2 ;
    nuki[ 0 * kd + 22 ] += +1 ;

    /*reaction 24: 2 H + H2 <=> 2 H2 */
    nuki[ 1 * kd + 23 ] += -2 ;
    nuki[ 0 * kd + 23 ] += -1 ;
    nuki[ 0 * kd + 23 ] += +2 ;

    /*reaction 25: 2 H + H2O <=> H2 + H2O */
    nuki[ 1 * kd + 24 ] += -2 ;
    nuki[ 5 * kd + 24 ] += -1 ;
    nuki[ 0 * kd + 24 ] += +1 ;
    nuki[ 5 * kd + 24 ] += +1 ;

    /*reaction 26: 2 H + CO2 <=> H2 + CO2 */
    nuki[ 1 * kd + 25 ] += -2 ;
    nuki[ 12 * kd + 25 ] += -1 ;
    nuki[ 0 * kd + 25 ] += +1 ;
    nuki[ 12 * kd + 25 ] += +1 ;

    /*reaction 27: H + OH + M <=> H2O + M */
    nuki[ 1 * kd + 26 ] += -1 ;
    nuki[ 4 * kd + 26 ] += -1 ;
    nuki[ 5 * kd + 26 ] += +1 ;

    /*reaction 28: H + HO2 <=> O2 + H2 */
    nuki[ 1 * kd + 27 ] += -1 ;
    nuki[ 6 * kd + 27 ] += -1 ;
    nuki[ 3 * kd + 27 ] += +1 ;
    nuki[ 0 * kd + 27 ] += +1 ;

    /*reaction 29: H + HO2 <=> 2 OH */
    nuki[ 1 * kd + 28 ] += -1 ;
    nuki[ 6 * kd + 28 ] += -1 ;
    nuki[ 4 * kd + 28 ] += +2 ;

    /*reaction 30: H + CH2 (+M) <=> CH3 (+M) */
    nuki[ 1 * kd + 29 ] += -1 ;
    nuki[ 7 * kd + 29 ] += -1 ;
    nuki[ 9 * kd + 29 ] += +1 ;

    /*reaction 31: H + CH3 (+M) <=> CH4 (+M) */
    nuki[ 1 * kd + 30 ] += -1 ;
    nuki[ 9 * kd + 30 ] += -1 ;
    nuki[ 10 * kd + 30 ] += +1 ;

    /*reaction 32: H + CH4 <=> CH3 + H2 */
    nuki[ 1 * kd + 31 ] += -1 ;
    nuki[ 10 * kd + 31 ] += -1 ;
    nuki[ 9 * kd + 31 ] += +1 ;
    nuki[ 0 * kd + 31 ] += +1 ;

    /*reaction 33: H + HCO (+M) <=> CH2O (+M) */
    nuki[ 1 * kd + 32 ] += -1 ;
    nuki[ 13 * kd + 32 ] += -1 ;
    nuki[ 14 * kd + 32 ] += +1 ;

    /*reaction 34: H + HCO <=> H2 + CO */
    nuki[ 1 * kd + 33 ] += -1 ;
    nuki[ 13 * kd + 33 ] += -1 ;
    nuki[ 0 * kd + 33 ] += +1 ;
    nuki[ 11 * kd + 33 ] += +1 ;

    /*reaction 35: H + CH2O (+M) <=> CH3O (+M) */
    nuki[ 1 * kd + 34 ] += -1 ;
    nuki[ 14 * kd + 34 ] += -1 ;
    nuki[ 15 * kd + 34 ] += +1 ;

    /*reaction 36: H + CH2O <=> HCO + H2 */
    nuki[ 1 * kd + 35 ] += -1 ;
    nuki[ 14 * kd + 35 ] += -1 ;
    nuki[ 13 * kd + 35 ] += +1 ;
    nuki[ 0 * kd + 35 ] += +1 ;

    /*reaction 37: H + CH3O <=> OH + CH3 */
    nuki[ 1 * kd + 36 ] += -1 ;
    nuki[ 15 * kd + 36 ] += -1 ;
    nuki[ 4 * kd + 36 ] += +1 ;
    nuki[ 9 * kd + 36 ] += +1 ;

    /*reaction 38: H + C2H4 (+M) <=> C2H5 (+M) */
    nuki[ 1 * kd + 37 ] += -1 ;
    nuki[ 16 * kd + 37 ] += -1 ;
    nuki[ 17 * kd + 37 ] += +1 ;

    /*reaction 39: H + C2H5 (+M) <=> C2H6 (+M) */
    nuki[ 1 * kd + 38 ] += -1 ;
    nuki[ 17 * kd + 38 ] += -1 ;
    nuki[ 18 * kd + 38 ] += +1 ;

    /*reaction 40: H + C2H6 <=> C2H5 + H2 */
    nuki[ 1 * kd + 39 ] += -1 ;
    nuki[ 18 * kd + 39 ] += -1 ;
    nuki[ 17 * kd + 39 ] += +1 ;
    nuki[ 0 * kd + 39 ] += +1 ;

    /*reaction 41: H2 + CO (+M) <=> CH2O (+M) */
    nuki[ 0 * kd + 40 ] += -1 ;
    nuki[ 11 * kd + 40 ] += -1 ;
    nuki[ 14 * kd + 40 ] += +1 ;

    /*reaction 42: OH + H2 <=> H + H2O */
    nuki[ 4 * kd + 41 ] += -1 ;
    nuki[ 0 * kd + 41 ] += -1 ;
    nuki[ 1 * kd + 41 ] += +1 ;
    nuki[ 5 * kd + 41 ] += +1 ;

    /*reaction 43: 2 OH <=> O + H2O */
    nuki[ 4 * kd + 42 ] += -2 ;
    nuki[ 2 * kd + 42 ] += +1 ;
    nuki[ 5 * kd + 42 ] += +1 ;

    /*reaction 44: OH + HO2 <=> O2 + H2O */
    nuki[ 4 * kd + 43 ] += -1 ;
    nuki[ 6 * kd + 43 ] += -1 ;
    nuki[ 3 * kd + 43 ] += +1 ;
    nuki[ 5 * kd + 43 ] += +1 ;

    /*reaction 45: OH + CH2 <=> H + CH2O */
    nuki[ 4 * kd + 44 ] += -1 ;
    nuki[ 7 * kd + 44 ] += -1 ;
    nuki[ 1 * kd + 44 ] += +1 ;
    nuki[ 14 * kd + 44 ] += +1 ;

    /*reaction 46: OH + CH2(S) <=> H + CH2O */
    nuki[ 4 * kd + 45 ] += -1 ;
    nuki[ 8 * kd + 45 ] += -1 ;
    nuki[ 1 * kd + 45 ] += +1 ;
    nuki[ 14 * kd + 45 ] += +1 ;

    /*reaction 47: OH + CH3 <=> CH2 + H2O */
    nuki[ 4 * kd + 46 ] += -1 ;
    nuki[ 9 * kd + 46 ] += -1 ;
    nuki[ 7 * kd + 46 ] += +1 ;
    nuki[ 5 * kd + 46 ] += +1 ;

    /*reaction 48: OH + CH3 <=> CH2(S) + H2O */
    nuki[ 4 * kd + 47 ] += -1 ;
    nuki[ 9 * kd + 47 ] += -1 ;
    nuki[ 8 * kd + 47 ] += +1 ;
    nuki[ 5 * kd + 47 ] += +1 ;

    /*reaction 49: OH + CH4 <=> CH3 + H2O */
    nuki[ 4 * kd + 48 ] += -1 ;
    nuki[ 10 * kd + 48 ] += -1 ;
    nuki[ 9 * kd + 48 ] += +1 ;
    nuki[ 5 * kd + 48 ] += +1 ;

    /*reaction 50: OH + CO <=> H + CO2 */
    nuki[ 4 * kd + 49 ] += -1 ;
    nuki[ 11 * kd + 49 ] += -1 ;
    nuki[ 1 * kd + 49 ] += +1 ;
    nuki[ 12 * kd + 49 ] += +1 ;

    /*reaction 51: OH + HCO <=> H2O + CO */
    nuki[ 4 * kd + 50 ] += -1 ;
    nuki[ 13 * kd + 50 ] += -1 ;
    nuki[ 5 * kd + 50 ] += +1 ;
    nuki[ 11 * kd + 50 ] += +1 ;

    /*reaction 52: OH + CH2O <=> HCO + H2O */
    nuki[ 4 * kd + 51 ] += -1 ;
    nuki[ 14 * kd + 51 ] += -1 ;
    nuki[ 13 * kd + 51 ] += +1 ;
    nuki[ 5 * kd + 51 ] += +1 ;

    /*reaction 53: OH + C2H6 <=> C2H5 + H2O */
    nuki[ 4 * kd + 52 ] += -1 ;
    nuki[ 18 * kd + 52 ] += -1 ;
    nuki[ 17 * kd + 52 ] += +1 ;
    nuki[ 5 * kd + 52 ] += +1 ;

    /*reaction 54: HO2 + CH2 <=> OH + CH2O */
    nuki[ 6 * kd + 53 ] += -1 ;
    nuki[ 7 * kd + 53 ] += -1 ;
    nuki[ 4 * kd + 53 ] += +1 ;
    nuki[ 14 * kd + 53 ] += +1 ;

    /*reaction 55: HO2 + CH3 <=> O2 + CH4 */
    nuki[ 6 * kd + 54 ] += -1 ;
    nuki[ 9 * kd + 54 ] += -1 ;
    nuki[ 3 * kd + 54 ] += +1 ;
    nuki[ 10 * kd + 54 ] += +1 ;

    /*reaction 56: HO2 + CH3 <=> OH + CH3O */
    nuki[ 6 * kd + 55 ] += -1 ;
    nuki[ 9 * kd + 55 ] += -1 ;
    nuki[ 4 * kd + 55 ] += +1 ;
    nuki[ 15 * kd + 55 ] += +1 ;

    /*reaction 57: HO2 + CO <=> OH + CO2 */
    nuki[ 6 * kd + 56 ] += -1 ;
    nuki[ 11 * kd + 56 ] += -1 ;
    nuki[ 4 * kd + 56 ] += +1 ;
    nuki[ 12 * kd + 56 ] += +1 ;

    /*reaction 58: CH2 + O2 <=> OH + HCO */
    nuki[ 7 * kd + 57 ] += -1 ;
    nuki[ 3 * kd + 57 ] += -1 ;
    nuki[ 4 * kd + 57 ] += +1 ;
    nuki[ 13 * kd + 57 ] += +1 ;

    /*reaction 59: CH2 + H2 <=> H + CH3 */
    nuki[ 7 * kd + 58 ] += -1 ;
    nuki[ 0 * kd + 58 ] += -1 ;
    nuki[ 1 * kd + 58 ] += +1 ;
    nuki[ 9 * kd + 58 ] += +1 ;

    /*reaction 60: CH2 + CH3 <=> H + C2H4 */
    nuki[ 7 * kd + 59 ] += -1 ;
    nuki[ 9 * kd + 59 ] += -1 ;
    nuki[ 1 * kd + 59 ] += +1 ;
    nuki[ 16 * kd + 59 ] += +1 ;

    /*reaction 61: CH2 + CH4 <=> 2 CH3 */
    nuki[ 7 * kd + 60 ] += -1 ;
    nuki[ 10 * kd + 60 ] += -1 ;
    nuki[ 9 * kd + 60 ] += +2 ;

    /*reaction 62: CH2(S) + N2 <=> CH2 + N2 */
    nuki[ 8 * kd + 61 ] += -1 ;
    nuki[ 19 * kd + 61 ] += -1 ;
    nuki[ 7 * kd + 61 ] += +1 ;
    nuki[ 19 * kd + 61 ] += +1 ;

    /*reaction 63: CH2(S) + AR <=> CH2 + AR */
    nuki[ 8 * kd + 62 ] += -1 ;
    nuki[ 20 * kd + 62 ] += -1 ;
    nuki[ 7 * kd + 62 ] += +1 ;
    nuki[ 20 * kd + 62 ] += +1 ;

    /*reaction 64: CH2(S) + O2 <=> H + OH + CO */
    nuki[ 8 * kd + 63 ] += -1 ;
    nuki[ 3 * kd + 63 ] += -1 ;
    nuki[ 1 * kd + 63 ] += +1 ;
    nuki[ 4 * kd + 63 ] += +1 ;
    nuki[ 11 * kd + 63 ] += +1 ;

    /*reaction 65: CH2(S) + O2 <=> CO + H2O */
    nuki[ 8 * kd + 64 ] += -1 ;
    nuki[ 3 * kd + 64 ] += -1 ;
    nuki[ 11 * kd + 64 ] += +1 ;
    nuki[ 5 * kd + 64 ] += +1 ;

    /*reaction 66: CH2(S) + H2 <=> CH3 + H */
    nuki[ 8 * kd + 65 ] += -1 ;
    nuki[ 0 * kd + 65 ] += -1 ;
    nuki[ 9 * kd + 65 ] += +1 ;
    nuki[ 1 * kd + 65 ] += +1 ;

    /*reaction 67: CH2(S) + H2O <=> CH2 + H2O */
    nuki[ 8 * kd + 66 ] += -1 ;
    nuki[ 5 * kd + 66 ] += -1 ;
    nuki[ 7 * kd + 66 ] += +1 ;
    nuki[ 5 * kd + 66 ] += +1 ;

    /*reaction 68: CH2(S) + CH3 <=> H + C2H4 */
    nuki[ 8 * kd + 67 ] += -1 ;
    nuki[ 9 * kd + 67 ] += -1 ;
    nuki[ 1 * kd + 67 ] += +1 ;
    nuki[ 16 * kd + 67 ] += +1 ;

    /*reaction 69: CH2(S) + CH4 <=> 2 CH3 */
    nuki[ 8 * kd + 68 ] += -1 ;
    nuki[ 10 * kd + 68 ] += -1 ;
    nuki[ 9 * kd + 68 ] += +2 ;

    /*reaction 70: CH2(S) + CO <=> CH2 + CO */
    nuki[ 8 * kd + 69 ] += -1 ;
    nuki[ 11 * kd + 69 ] += -1 ;
    nuki[ 7 * kd + 69 ] += +1 ;
    nuki[ 11 * kd + 69 ] += +1 ;

    /*reaction 71: CH2(S) + CO2 <=> CH2 + CO2 */
    nuki[ 8 * kd + 70 ] += -1 ;
    nuki[ 12 * kd + 70 ] += -1 ;
    nuki[ 7 * kd + 70 ] += +1 ;
    nuki[ 12 * kd + 70 ] += +1 ;

    /*reaction 72: CH2(S) + CO2 <=> CO + CH2O */
    nuki[ 8 * kd + 71 ] += -1 ;
    nuki[ 12 * kd + 71 ] += -1 ;
    nuki[ 11 * kd + 71 ] += +1 ;
    nuki[ 14 * kd + 71 ] += +1 ;

    /*reaction 73: CH3 + O2 <=> O + CH3O */
    nuki[ 9 * kd + 72 ] += -1 ;
    nuki[ 3 * kd + 72 ] += -1 ;
    nuki[ 2 * kd + 72 ] += +1 ;
    nuki[ 15 * kd + 72 ] += +1 ;

    /*reaction 74: CH3 + O2 <=> OH + CH2O */
    nuki[ 9 * kd + 73 ] += -1 ;
    nuki[ 3 * kd + 73 ] += -1 ;
    nuki[ 4 * kd + 73 ] += +1 ;
    nuki[ 14 * kd + 73 ] += +1 ;

    /*reaction 75: 2 CH3 (+M) <=> C2H6 (+M) */
    nuki[ 9 * kd + 74 ] += -2 ;
    nuki[ 18 * kd + 74 ] += +1 ;

    /*reaction 76: 2 CH3 <=> H + C2H5 */
    nuki[ 9 * kd + 75 ] += -2 ;
    nuki[ 1 * kd + 75 ] += +1 ;
    nuki[ 17 * kd + 75 ] += +1 ;

    /*reaction 77: CH3 + HCO <=> CH4 + CO */
    nuki[ 9 * kd + 76 ] += -1 ;
    nuki[ 13 * kd + 76 ] += -1 ;
    nuki[ 10 * kd + 76 ] += +1 ;
    nuki[ 11 * kd + 76 ] += +1 ;

    /*reaction 78: CH3 + CH2O <=> HCO + CH4 */
    nuki[ 9 * kd + 77 ] += -1 ;
    nuki[ 14 * kd + 77 ] += -1 ;
    nuki[ 13 * kd + 77 ] += +1 ;
    nuki[ 10 * kd + 77 ] += +1 ;

    /*reaction 79: CH3 + C2H6 <=> C2H5 + CH4 */
    nuki[ 9 * kd + 78 ] += -1 ;
    nuki[ 18 * kd + 78 ] += -1 ;
    nuki[ 17 * kd + 78 ] += +1 ;
    nuki[ 10 * kd + 78 ] += +1 ;

    /*reaction 80: HCO + H2O <=> H + CO + H2O */
    nuki[ 13 * kd + 79 ] += -1 ;
    nuki[ 5 * kd + 79 ] += -1 ;
    nuki[ 1 * kd + 79 ] += +1 ;
    nuki[ 11 * kd + 79 ] += +1 ;
    nuki[ 5 * kd + 79 ] += +1 ;

    /*reaction 81: HCO + M <=> H + CO + M */
    nuki[ 13 * kd + 80 ] += -1 ;
    nuki[ 1 * kd + 80 ] += +1 ;
    nuki[ 11 * kd + 80 ] += +1 ;

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    nuki[ 13 * kd + 81 ] += -1 ;
    nuki[ 3 * kd + 81 ] += -1 ;
    nuki[ 6 * kd + 81 ] += +1 ;
    nuki[ 11 * kd + 81 ] += +1 ;

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    nuki[ 15 * kd + 82 ] += -1 ;
    nuki[ 3 * kd + 82 ] += -1 ;
    nuki[ 6 * kd + 82 ] += +1 ;
    nuki[ 14 * kd + 82 ] += +1 ;

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    nuki[ 17 * kd + 83 ] += -1 ;
    nuki[ 3 * kd + 83 ] += -1 ;
    nuki[ 6 * kd + 83 ] += +1 ;
    nuki[ 16 * kd + 83 ] += +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < 5 * 21; ++ id) {
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

    /*CH2 */
    ncf[ 7 * kd + 2 ] = 1; /*C */
    ncf[ 7 * kd + 1 ] = 2; /*H */

    /*CH2(S) */
    ncf[ 8 * kd + 2 ] = 1; /*C */
    ncf[ 8 * kd + 1 ] = 2; /*H */

    /*CH3 */
    ncf[ 9 * kd + 2 ] = 1; /*C */
    ncf[ 9 * kd + 1 ] = 3; /*H */

    /*CH4 */
    ncf[ 10 * kd + 2 ] = 1; /*C */
    ncf[ 10 * kd + 1 ] = 4; /*H */

    /*CO */
    ncf[ 11 * kd + 2 ] = 1; /*C */
    ncf[ 11 * kd + 0 ] = 1; /*O */

    /*CO2 */
    ncf[ 12 * kd + 2 ] = 1; /*C */
    ncf[ 12 * kd + 0 ] = 2; /*O */

    /*HCO */
    ncf[ 13 * kd + 1 ] = 1; /*H */
    ncf[ 13 * kd + 2 ] = 1; /*C */
    ncf[ 13 * kd + 0 ] = 1; /*O */

    /*CH2O */
    ncf[ 14 * kd + 1 ] = 2; /*H */
    ncf[ 14 * kd + 2 ] = 1; /*C */
    ncf[ 14 * kd + 0 ] = 1; /*O */

    /*CH3O */
    ncf[ 15 * kd + 2 ] = 1; /*C */
    ncf[ 15 * kd + 1 ] = 3; /*H */
    ncf[ 15 * kd + 0 ] = 1; /*O */

    /*C2H4 */
    ncf[ 16 * kd + 2 ] = 2; /*C */
    ncf[ 16 * kd + 1 ] = 4; /*H */

    /*C2H5 */
    ncf[ 17 * kd + 2 ] = 2; /*C */
    ncf[ 17 * kd + 1 ] = 5; /*H */

    /*C2H6 */
    ncf[ 18 * kd + 2 ] = 2; /*C */
    ncf[ 18 * kd + 1 ] = 6; /*H */

    /*N2 */
    ncf[ 19 * kd + 3 ] = 2; /*N */

    /*AR */
    ncf[ 20 * kd + 4 ] = 1; /*AR */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE(int * iwrk, double * rwrk, double * a, double * b, double * e)
{

    /*reaction 1: O + H + M <=> OH + M */
    a[0] = 5e+17;
    b[0] = -1;
    e[0] = 0;

    /*reaction 2: O + H2 <=> H + OH */
    a[1] = 50000;
    b[1] = 2.67;
    e[1] = 6290;

    /*reaction 3: O + HO2 <=> OH + O2 */
    a[2] = 2e+13;
    b[2] = 0;
    e[2] = 0;

    /*reaction 4: O + CH2 <=> H + HCO */
    a[3] = 8e+13;
    b[3] = 0;
    e[3] = 0;

    /*reaction 5: O + CH2(S) <=> H + HCO */
    a[4] = 1.5e+13;
    b[4] = 0;
    e[4] = 0;

    /*reaction 6: O + CH3 <=> H + CH2O */
    a[5] = 8.43e+13;
    b[5] = 0;
    e[5] = 0;

    /*reaction 7: O + CH4 <=> OH + CH3 */
    a[6] = 1.02e+09;
    b[6] = 1.5;
    e[6] = 8600;

    /*reaction 8: O + CO + M <=> CO2 + M */
    a[7] = 6.02e+14;
    b[7] = 0;
    e[7] = 3000;

    /*reaction 9: O + HCO <=> OH + CO */
    a[8] = 3e+13;
    b[8] = 0;
    e[8] = 0;

    /*reaction 10: O + HCO <=> H + CO2 */
    a[9] = 3e+13;
    b[9] = 0;
    e[9] = 0;

    /*reaction 11: O + CH2O <=> OH + HCO */
    a[10] = 3.9e+13;
    b[10] = 0;
    e[10] = 3540;

    /*reaction 12: O + C2H4 <=> CH3 + HCO */
    a[11] = 1.92e+07;
    b[11] = 1.83;
    e[11] = 220;

    /*reaction 13: O + C2H5 <=> CH3 + CH2O */
    a[12] = 1.32e+14;
    b[12] = 0;
    e[12] = 0;

    /*reaction 14: O + C2H6 <=> OH + C2H5 */
    a[13] = 8.98e+07;
    b[13] = 1.92;
    e[13] = 5690;

    /*reaction 15: O2 + CO <=> O + CO2 */
    a[14] = 2.5e+12;
    b[14] = 0;
    e[14] = 47800;

    /*reaction 16: O2 + CH2O <=> HO2 + HCO */
    a[15] = 1e+14;
    b[15] = 0;
    e[15] = 40000;

    /*reaction 17: H + O2 + M <=> HO2 + M */
    a[16] = 2.8e+18;
    b[16] = -0.86;
    e[16] = 0;

    /*reaction 18: H + 2 O2 <=> HO2 + O2 */
    a[17] = 3e+20;
    b[17] = -1.72;
    e[17] = 0;

    /*reaction 19: H + O2 + H2O <=> HO2 + H2O */
    a[18] = 9.38e+18;
    b[18] = -0.76;
    e[18] = 0;

    /*reaction 20: H + O2 + N2 <=> HO2 + N2 */
    a[19] = 3.75e+20;
    b[19] = -1.72;
    e[19] = 0;

    /*reaction 21: H + O2 + AR <=> HO2 + AR */
    a[20] = 7e+17;
    b[20] = -0.8;
    e[20] = 0;

    /*reaction 22: H + O2 <=> O + OH */
    a[21] = 8.3e+13;
    b[21] = 0;
    e[21] = 14413;

    /*reaction 23: 2 H + M <=> H2 + M */
    a[22] = 1e+18;
    b[22] = -1;
    e[22] = 0;

    /*reaction 24: 2 H + H2 <=> 2 H2 */
    a[23] = 9e+16;
    b[23] = -0.6;
    e[23] = 0;

    /*reaction 25: 2 H + H2O <=> H2 + H2O */
    a[24] = 6e+19;
    b[24] = -1.25;
    e[24] = 0;

    /*reaction 26: 2 H + CO2 <=> H2 + CO2 */
    a[25] = 5.5e+20;
    b[25] = -2;
    e[25] = 0;

    /*reaction 27: H + OH + M <=> H2O + M */
    a[26] = 2.2e+22;
    b[26] = -2;
    e[26] = 0;

    /*reaction 28: H + HO2 <=> O2 + H2 */
    a[27] = 2.8e+13;
    b[27] = 0;
    e[27] = 1068;

    /*reaction 29: H + HO2 <=> 2 OH */
    a[28] = 1.34e+14;
    b[28] = 0;
    e[28] = 635;

    /*reaction 30: H + CH2 (+M) <=> CH3 (+M) */
    a[29] = 2.5e+16;
    b[29] = -0.8;
    e[29] = 0;

    /*reaction 31: H + CH3 (+M) <=> CH4 (+M) */
    a[30] = 1.27e+16;
    b[30] = -0.63;
    e[30] = 383;

    /*reaction 32: H + CH4 <=> CH3 + H2 */
    a[31] = 6.6e+08;
    b[31] = 1.62;
    e[31] = 10840;

    /*reaction 33: H + HCO (+M) <=> CH2O (+M) */
    a[32] = 1.09e+12;
    b[32] = 0.48;
    e[32] = -260;

    /*reaction 34: H + HCO <=> H2 + CO */
    a[33] = 7.34e+13;
    b[33] = 0;
    e[33] = 0;

    /*reaction 35: H + CH2O (+M) <=> CH3O (+M) */
    a[34] = 5.4e+11;
    b[34] = 0.454;
    e[34] = 2600;

    /*reaction 36: H + CH2O <=> HCO + H2 */
    a[35] = 2.3e+10;
    b[35] = 1.05;
    e[35] = 3275;

    /*reaction 37: H + CH3O <=> OH + CH3 */
    a[36] = 3.2e+13;
    b[36] = 0;
    e[36] = 0;

    /*reaction 38: H + C2H4 (+M) <=> C2H5 (+M) */
    a[37] = 1.08e+12;
    b[37] = 0.454;
    e[37] = 1820;

    /*reaction 39: H + C2H5 (+M) <=> C2H6 (+M) */
    a[38] = 5.21e+17;
    b[38] = -0.99;
    e[38] = 1580;

    /*reaction 40: H + C2H6 <=> C2H5 + H2 */
    a[39] = 1.15e+08;
    b[39] = 1.9;
    e[39] = 7530;

    /*reaction 41: H2 + CO (+M) <=> CH2O (+M) */
    a[40] = 4.3e+07;
    b[40] = 1.5;
    e[40] = 79600;

    /*reaction 42: OH + H2 <=> H + H2O */
    a[41] = 2.16e+08;
    b[41] = 1.51;
    e[41] = 3430;

    /*reaction 43: 2 OH <=> O + H2O */
    a[42] = 35700;
    b[42] = 2.4;
    e[42] = -2110;

    /*reaction 44: OH + HO2 <=> O2 + H2O */
    a[43] = 2.9e+13;
    b[43] = 0;
    e[43] = -500;

    /*reaction 45: OH + CH2 <=> H + CH2O */
    a[44] = 2e+13;
    b[44] = 0;
    e[44] = 0;

    /*reaction 46: OH + CH2(S) <=> H + CH2O */
    a[45] = 3e+13;
    b[45] = 0;
    e[45] = 0;

    /*reaction 47: OH + CH3 <=> CH2 + H2O */
    a[46] = 5.6e+07;
    b[46] = 1.6;
    e[46] = 5420;

    /*reaction 48: OH + CH3 <=> CH2(S) + H2O */
    a[47] = 2.501e+13;
    b[47] = 0;
    e[47] = 0;

    /*reaction 49: OH + CH4 <=> CH3 + H2O */
    a[48] = 1e+08;
    b[48] = 1.6;
    e[48] = 3120;

    /*reaction 50: OH + CO <=> H + CO2 */
    a[49] = 4.76e+07;
    b[49] = 1.228;
    e[49] = 70;

    /*reaction 51: OH + HCO <=> H2O + CO */
    a[50] = 5e+13;
    b[50] = 0;
    e[50] = 0;

    /*reaction 52: OH + CH2O <=> HCO + H2O */
    a[51] = 3.43e+09;
    b[51] = 1.18;
    e[51] = -447;

    /*reaction 53: OH + C2H6 <=> C2H5 + H2O */
    a[52] = 3.54e+06;
    b[52] = 2.12;
    e[52] = 870;

    /*reaction 54: HO2 + CH2 <=> OH + CH2O */
    a[53] = 2e+13;
    b[53] = 0;
    e[53] = 0;

    /*reaction 55: HO2 + CH3 <=> O2 + CH4 */
    a[54] = 1e+12;
    b[54] = 0;
    e[54] = 0;

    /*reaction 56: HO2 + CH3 <=> OH + CH3O */
    a[55] = 2e+13;
    b[55] = 0;
    e[55] = 0;

    /*reaction 57: HO2 + CO <=> OH + CO2 */
    a[56] = 1.5e+14;
    b[56] = 0;
    e[56] = 23600;

    /*reaction 58: CH2 + O2 <=> OH + HCO */
    a[57] = 1.32e+13;
    b[57] = 0;
    e[57] = 1500;

    /*reaction 59: CH2 + H2 <=> H + CH3 */
    a[58] = 500000;
    b[58] = 2;
    e[58] = 7230;

    /*reaction 60: CH2 + CH3 <=> H + C2H4 */
    a[59] = 4e+13;
    b[59] = 0;
    e[59] = 0;

    /*reaction 61: CH2 + CH4 <=> 2 CH3 */
    a[60] = 2.46e+06;
    b[60] = 2;
    e[60] = 8270;

    /*reaction 62: CH2(S) + N2 <=> CH2 + N2 */
    a[61] = 1.5e+13;
    b[61] = 0;
    e[61] = 600;

    /*reaction 63: CH2(S) + AR <=> CH2 + AR */
    a[62] = 9e+12;
    b[62] = 0;
    e[62] = 600;

    /*reaction 64: CH2(S) + O2 <=> H + OH + CO */
    a[63] = 2.8e+13;
    b[63] = 0;
    e[63] = 0;

    /*reaction 65: CH2(S) + O2 <=> CO + H2O */
    a[64] = 1.2e+13;
    b[64] = 0;
    e[64] = 0;

    /*reaction 66: CH2(S) + H2 <=> CH3 + H */
    a[65] = 7e+13;
    b[65] = 0;
    e[65] = 0;

    /*reaction 67: CH2(S) + H2O <=> CH2 + H2O */
    a[66] = 3e+13;
    b[66] = 0;
    e[66] = 0;

    /*reaction 68: CH2(S) + CH3 <=> H + C2H4 */
    a[67] = 1.2e+13;
    b[67] = 0;
    e[67] = -570;

    /*reaction 69: CH2(S) + CH4 <=> 2 CH3 */
    a[68] = 1.6e+13;
    b[68] = 0;
    e[68] = -570;

    /*reaction 70: CH2(S) + CO <=> CH2 + CO */
    a[69] = 9e+12;
    b[69] = 0;
    e[69] = 0;

    /*reaction 71: CH2(S) + CO2 <=> CH2 + CO2 */
    a[70] = 7e+12;
    b[70] = 0;
    e[70] = 0;

    /*reaction 72: CH2(S) + CO2 <=> CO + CH2O */
    a[71] = 1.4e+13;
    b[71] = 0;
    e[71] = 0;

    /*reaction 73: CH3 + O2 <=> O + CH3O */
    a[72] = 2.675e+13;
    b[72] = 0;
    e[72] = 28800;

    /*reaction 74: CH3 + O2 <=> OH + CH2O */
    a[73] = 3.6e+10;
    b[73] = 0;
    e[73] = 8940;

    /*reaction 75: 2 CH3 (+M) <=> C2H6 (+M) */
    a[74] = 2.12e+16;
    b[74] = -0.97;
    e[74] = 620;

    /*reaction 76: 2 CH3 <=> H + C2H5 */
    a[75] = 4.99e+12;
    b[75] = 0.1;
    e[75] = 10600;

    /*reaction 77: CH3 + HCO <=> CH4 + CO */
    a[76] = 2.648e+13;
    b[76] = 0;
    e[76] = 0;

    /*reaction 78: CH3 + CH2O <=> HCO + CH4 */
    a[77] = 3320;
    b[77] = 2.81;
    e[77] = 5860;

    /*reaction 79: CH3 + C2H6 <=> C2H5 + CH4 */
    a[78] = 6.14e+06;
    b[78] = 1.74;
    e[78] = 10450;

    /*reaction 80: HCO + H2O <=> H + CO + H2O */
    a[79] = 2.244e+18;
    b[79] = -1;
    e[79] = 17000;

    /*reaction 81: HCO + M <=> H + CO + M */
    a[80] = 1.87e+17;
    b[80] = -1;
    e[80] = 17000;

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    a[81] = 7.6e+12;
    b[81] = 0;
    e[81] = 400;

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    a[82] = 4.28e-13;
    b[82] = 7.6;
    e[82] = -3530;

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    a[83] = 8.4e+11;
    b[83] = 0;
    e[83] = 3875;

    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double * T, double * C, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[21]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: O + H + M <=> OH + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: O + H2 <=> H + OH */
    /*eqcon[1] *= 1;  */

    /*reaction 3: O + HO2 <=> OH + O2 */
    /*eqcon[2] *= 1;  */

    /*reaction 4: O + CH2 <=> H + HCO */
    /*eqcon[3] *= 1;  */

    /*reaction 5: O + CH2(S) <=> H + HCO */
    /*eqcon[4] *= 1;  */

    /*reaction 6: O + CH3 <=> H + CH2O */
    /*eqcon[5] *= 1;  */

    /*reaction 7: O + CH4 <=> OH + CH3 */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + CO + M <=> CO2 + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: O + HCO <=> OH + CO */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + HCO <=> H + CO2 */
    /*eqcon[9] *= 1;  */

    /*reaction 11: O + CH2O <=> OH + HCO */
    /*eqcon[10] *= 1;  */

    /*reaction 12: O + C2H4 <=> CH3 + HCO */
    /*eqcon[11] *= 1;  */

    /*reaction 13: O + C2H5 <=> CH3 + CH2O */
    /*eqcon[12] *= 1;  */

    /*reaction 14: O + C2H6 <=> OH + C2H5 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O2 + CO <=> O + CO2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: O2 + CH2O <=> HO2 + HCO */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + O2 + M <=> HO2 + M */
    eqcon[16] *= 1e+06; 

    /*reaction 18: H + 2 O2 <=> HO2 + O2 */
    eqcon[17] *= 1e+06; 

    /*reaction 19: H + O2 + H2O <=> HO2 + H2O */
    eqcon[18] *= 1e+06; 

    /*reaction 20: H + O2 + N2 <=> HO2 + N2 */
    eqcon[19] *= 1e+06; 

    /*reaction 21: H + O2 + AR <=> HO2 + AR */
    eqcon[20] *= 1e+06; 

    /*reaction 22: H + O2 <=> O + OH */
    /*eqcon[21] *= 1;  */

    /*reaction 23: 2 H + M <=> H2 + M */
    eqcon[22] *= 1e+06; 

    /*reaction 24: 2 H + H2 <=> 2 H2 */
    eqcon[23] *= 1e+06; 

    /*reaction 25: 2 H + H2O <=> H2 + H2O */
    eqcon[24] *= 1e+06; 

    /*reaction 26: 2 H + CO2 <=> H2 + CO2 */
    eqcon[25] *= 1e+06; 

    /*reaction 27: H + OH + M <=> H2O + M */
    eqcon[26] *= 1e+06; 

    /*reaction 28: H + HO2 <=> O2 + H2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H + HO2 <=> 2 OH */
    /*eqcon[28] *= 1;  */

    /*reaction 30: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[29] *= 1e+06; 

    /*reaction 31: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[30] *= 1e+06; 

    /*reaction 32: H + CH4 <=> CH3 + H2 */
    /*eqcon[31] *= 1;  */

    /*reaction 33: H + HCO (+M) <=> CH2O (+M) */
    eqcon[32] *= 1e+06; 

    /*reaction 34: H + HCO <=> H2 + CO */
    /*eqcon[33] *= 1;  */

    /*reaction 35: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[34] *= 1e+06; 

    /*reaction 36: H + CH2O <=> HCO + H2 */
    /*eqcon[35] *= 1;  */

    /*reaction 37: H + CH3O <=> OH + CH3 */
    /*eqcon[36] *= 1;  */

    /*reaction 38: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[37] *= 1e+06; 

    /*reaction 39: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[38] *= 1e+06; 

    /*reaction 40: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[40] *= 1e+06; 

    /*reaction 42: OH + H2 <=> H + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: 2 OH <=> O + H2O */
    /*eqcon[42] *= 1;  */

    /*reaction 44: OH + HO2 <=> O2 + H2O */
    /*eqcon[43] *= 1;  */

    /*reaction 45: OH + CH2 <=> H + CH2O */
    /*eqcon[44] *= 1;  */

    /*reaction 46: OH + CH2(S) <=> H + CH2O */
    /*eqcon[45] *= 1;  */

    /*reaction 47: OH + CH3 <=> CH2 + H2O */
    /*eqcon[46] *= 1;  */

    /*reaction 48: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[47] *= 1;  */

    /*reaction 49: OH + CH4 <=> CH3 + H2O */
    /*eqcon[48] *= 1;  */

    /*reaction 50: OH + CO <=> H + CO2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: OH + HCO <=> H2O + CO */
    /*eqcon[50] *= 1;  */

    /*reaction 52: OH + CH2O <=> HCO + H2O */
    /*eqcon[51] *= 1;  */

    /*reaction 53: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[52] *= 1;  */

    /*reaction 54: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[53] *= 1;  */

    /*reaction 55: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[54] *= 1;  */

    /*reaction 56: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: HO2 + CO <=> OH + CO2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: CH2 + O2 <=> OH + HCO */
    /*eqcon[57] *= 1;  */

    /*reaction 59: CH2 + H2 <=> H + CH3 */
    /*eqcon[58] *= 1;  */

    /*reaction 60: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[60] *= 1;  */

    /*reaction 62: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[61] *= 1;  */

    /*reaction 63: CH2(S) + AR <=> CH2 + AR */
    /*eqcon[62] *= 1;  */

    /*reaction 64: CH2(S) + O2 <=> H + OH + CO */
    eqcon[63] *= 1e-06; 

    /*reaction 65: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[64] *= 1;  */

    /*reaction 66: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[65] *= 1;  */

    /*reaction 67: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[68] *= 1;  */

    /*reaction 70: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[69] *= 1;  */

    /*reaction 71: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[70] *= 1;  */

    /*reaction 72: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[71] *= 1;  */

    /*reaction 73: CH3 + O2 <=> O + CH3O */
    /*eqcon[72] *= 1;  */

    /*reaction 74: CH3 + O2 <=> OH + CH2O */
    /*eqcon[73] *= 1;  */

    /*reaction 75: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[74] *= 1e+06; 

    /*reaction 76: 2 CH3 <=> H + C2H5 */
    /*eqcon[75] *= 1;  */

    /*reaction 77: CH3 + HCO <=> CH4 + CO */
    /*eqcon[76] *= 1;  */

    /*reaction 78: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[77] *= 1;  */

    /*reaction 79: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: HCO + H2O <=> H + CO + H2O */
    eqcon[79] *= 1e-06; 

    /*reaction 81: HCO + M <=> H + CO + M */
    eqcon[80] *= 1e-06; 

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    /*eqcon[81] *= 1;  */

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[82] *= 1;  */

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[83] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[21]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: O + H + M <=> OH + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: O + H2 <=> H + OH */
    /*eqcon[1] *= 1;  */

    /*reaction 3: O + HO2 <=> OH + O2 */
    /*eqcon[2] *= 1;  */

    /*reaction 4: O + CH2 <=> H + HCO */
    /*eqcon[3] *= 1;  */

    /*reaction 5: O + CH2(S) <=> H + HCO */
    /*eqcon[4] *= 1;  */

    /*reaction 6: O + CH3 <=> H + CH2O */
    /*eqcon[5] *= 1;  */

    /*reaction 7: O + CH4 <=> OH + CH3 */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + CO + M <=> CO2 + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: O + HCO <=> OH + CO */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + HCO <=> H + CO2 */
    /*eqcon[9] *= 1;  */

    /*reaction 11: O + CH2O <=> OH + HCO */
    /*eqcon[10] *= 1;  */

    /*reaction 12: O + C2H4 <=> CH3 + HCO */
    /*eqcon[11] *= 1;  */

    /*reaction 13: O + C2H5 <=> CH3 + CH2O */
    /*eqcon[12] *= 1;  */

    /*reaction 14: O + C2H6 <=> OH + C2H5 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O2 + CO <=> O + CO2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: O2 + CH2O <=> HO2 + HCO */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + O2 + M <=> HO2 + M */
    eqcon[16] *= 1e+06; 

    /*reaction 18: H + 2 O2 <=> HO2 + O2 */
    eqcon[17] *= 1e+06; 

    /*reaction 19: H + O2 + H2O <=> HO2 + H2O */
    eqcon[18] *= 1e+06; 

    /*reaction 20: H + O2 + N2 <=> HO2 + N2 */
    eqcon[19] *= 1e+06; 

    /*reaction 21: H + O2 + AR <=> HO2 + AR */
    eqcon[20] *= 1e+06; 

    /*reaction 22: H + O2 <=> O + OH */
    /*eqcon[21] *= 1;  */

    /*reaction 23: 2 H + M <=> H2 + M */
    eqcon[22] *= 1e+06; 

    /*reaction 24: 2 H + H2 <=> 2 H2 */
    eqcon[23] *= 1e+06; 

    /*reaction 25: 2 H + H2O <=> H2 + H2O */
    eqcon[24] *= 1e+06; 

    /*reaction 26: 2 H + CO2 <=> H2 + CO2 */
    eqcon[25] *= 1e+06; 

    /*reaction 27: H + OH + M <=> H2O + M */
    eqcon[26] *= 1e+06; 

    /*reaction 28: H + HO2 <=> O2 + H2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H + HO2 <=> 2 OH */
    /*eqcon[28] *= 1;  */

    /*reaction 30: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[29] *= 1e+06; 

    /*reaction 31: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[30] *= 1e+06; 

    /*reaction 32: H + CH4 <=> CH3 + H2 */
    /*eqcon[31] *= 1;  */

    /*reaction 33: H + HCO (+M) <=> CH2O (+M) */
    eqcon[32] *= 1e+06; 

    /*reaction 34: H + HCO <=> H2 + CO */
    /*eqcon[33] *= 1;  */

    /*reaction 35: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[34] *= 1e+06; 

    /*reaction 36: H + CH2O <=> HCO + H2 */
    /*eqcon[35] *= 1;  */

    /*reaction 37: H + CH3O <=> OH + CH3 */
    /*eqcon[36] *= 1;  */

    /*reaction 38: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[37] *= 1e+06; 

    /*reaction 39: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[38] *= 1e+06; 

    /*reaction 40: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[40] *= 1e+06; 

    /*reaction 42: OH + H2 <=> H + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: 2 OH <=> O + H2O */
    /*eqcon[42] *= 1;  */

    /*reaction 44: OH + HO2 <=> O2 + H2O */
    /*eqcon[43] *= 1;  */

    /*reaction 45: OH + CH2 <=> H + CH2O */
    /*eqcon[44] *= 1;  */

    /*reaction 46: OH + CH2(S) <=> H + CH2O */
    /*eqcon[45] *= 1;  */

    /*reaction 47: OH + CH3 <=> CH2 + H2O */
    /*eqcon[46] *= 1;  */

    /*reaction 48: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[47] *= 1;  */

    /*reaction 49: OH + CH4 <=> CH3 + H2O */
    /*eqcon[48] *= 1;  */

    /*reaction 50: OH + CO <=> H + CO2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: OH + HCO <=> H2O + CO */
    /*eqcon[50] *= 1;  */

    /*reaction 52: OH + CH2O <=> HCO + H2O */
    /*eqcon[51] *= 1;  */

    /*reaction 53: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[52] *= 1;  */

    /*reaction 54: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[53] *= 1;  */

    /*reaction 55: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[54] *= 1;  */

    /*reaction 56: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: HO2 + CO <=> OH + CO2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: CH2 + O2 <=> OH + HCO */
    /*eqcon[57] *= 1;  */

    /*reaction 59: CH2 + H2 <=> H + CH3 */
    /*eqcon[58] *= 1;  */

    /*reaction 60: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[60] *= 1;  */

    /*reaction 62: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[61] *= 1;  */

    /*reaction 63: CH2(S) + AR <=> CH2 + AR */
    /*eqcon[62] *= 1;  */

    /*reaction 64: CH2(S) + O2 <=> H + OH + CO */
    eqcon[63] *= 1e-06; 

    /*reaction 65: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[64] *= 1;  */

    /*reaction 66: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[65] *= 1;  */

    /*reaction 67: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[68] *= 1;  */

    /*reaction 70: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[69] *= 1;  */

    /*reaction 71: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[70] *= 1;  */

    /*reaction 72: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[71] *= 1;  */

    /*reaction 73: CH3 + O2 <=> O + CH3O */
    /*eqcon[72] *= 1;  */

    /*reaction 74: CH3 + O2 <=> OH + CH2O */
    /*eqcon[73] *= 1;  */

    /*reaction 75: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[74] *= 1e+06; 

    /*reaction 76: 2 CH3 <=> H + C2H5 */
    /*eqcon[75] *= 1;  */

    /*reaction 77: CH3 + HCO <=> CH4 + CO */
    /*eqcon[76] *= 1;  */

    /*reaction 78: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[77] *= 1;  */

    /*reaction 79: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: HCO + H2O <=> H + CO + H2O */
    eqcon[79] *= 1e-06; 

    /*reaction 81: HCO + M <=> H + CO + M */
    eqcon[80] *= 1e-06; 

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    /*eqcon[81] *= 1;  */

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[82] *= 1;  */

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[83] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[21]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: O + H + M <=> OH + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: O + H2 <=> H + OH */
    /*eqcon[1] *= 1;  */

    /*reaction 3: O + HO2 <=> OH + O2 */
    /*eqcon[2] *= 1;  */

    /*reaction 4: O + CH2 <=> H + HCO */
    /*eqcon[3] *= 1;  */

    /*reaction 5: O + CH2(S) <=> H + HCO */
    /*eqcon[4] *= 1;  */

    /*reaction 6: O + CH3 <=> H + CH2O */
    /*eqcon[5] *= 1;  */

    /*reaction 7: O + CH4 <=> OH + CH3 */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + CO + M <=> CO2 + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: O + HCO <=> OH + CO */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + HCO <=> H + CO2 */
    /*eqcon[9] *= 1;  */

    /*reaction 11: O + CH2O <=> OH + HCO */
    /*eqcon[10] *= 1;  */

    /*reaction 12: O + C2H4 <=> CH3 + HCO */
    /*eqcon[11] *= 1;  */

    /*reaction 13: O + C2H5 <=> CH3 + CH2O */
    /*eqcon[12] *= 1;  */

    /*reaction 14: O + C2H6 <=> OH + C2H5 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O2 + CO <=> O + CO2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: O2 + CH2O <=> HO2 + HCO */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + O2 + M <=> HO2 + M */
    eqcon[16] *= 1e+06; 

    /*reaction 18: H + 2 O2 <=> HO2 + O2 */
    eqcon[17] *= 1e+06; 

    /*reaction 19: H + O2 + H2O <=> HO2 + H2O */
    eqcon[18] *= 1e+06; 

    /*reaction 20: H + O2 + N2 <=> HO2 + N2 */
    eqcon[19] *= 1e+06; 

    /*reaction 21: H + O2 + AR <=> HO2 + AR */
    eqcon[20] *= 1e+06; 

    /*reaction 22: H + O2 <=> O + OH */
    /*eqcon[21] *= 1;  */

    /*reaction 23: 2 H + M <=> H2 + M */
    eqcon[22] *= 1e+06; 

    /*reaction 24: 2 H + H2 <=> 2 H2 */
    eqcon[23] *= 1e+06; 

    /*reaction 25: 2 H + H2O <=> H2 + H2O */
    eqcon[24] *= 1e+06; 

    /*reaction 26: 2 H + CO2 <=> H2 + CO2 */
    eqcon[25] *= 1e+06; 

    /*reaction 27: H + OH + M <=> H2O + M */
    eqcon[26] *= 1e+06; 

    /*reaction 28: H + HO2 <=> O2 + H2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H + HO2 <=> 2 OH */
    /*eqcon[28] *= 1;  */

    /*reaction 30: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[29] *= 1e+06; 

    /*reaction 31: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[30] *= 1e+06; 

    /*reaction 32: H + CH4 <=> CH3 + H2 */
    /*eqcon[31] *= 1;  */

    /*reaction 33: H + HCO (+M) <=> CH2O (+M) */
    eqcon[32] *= 1e+06; 

    /*reaction 34: H + HCO <=> H2 + CO */
    /*eqcon[33] *= 1;  */

    /*reaction 35: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[34] *= 1e+06; 

    /*reaction 36: H + CH2O <=> HCO + H2 */
    /*eqcon[35] *= 1;  */

    /*reaction 37: H + CH3O <=> OH + CH3 */
    /*eqcon[36] *= 1;  */

    /*reaction 38: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[37] *= 1e+06; 

    /*reaction 39: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[38] *= 1e+06; 

    /*reaction 40: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[40] *= 1e+06; 

    /*reaction 42: OH + H2 <=> H + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: 2 OH <=> O + H2O */
    /*eqcon[42] *= 1;  */

    /*reaction 44: OH + HO2 <=> O2 + H2O */
    /*eqcon[43] *= 1;  */

    /*reaction 45: OH + CH2 <=> H + CH2O */
    /*eqcon[44] *= 1;  */

    /*reaction 46: OH + CH2(S) <=> H + CH2O */
    /*eqcon[45] *= 1;  */

    /*reaction 47: OH + CH3 <=> CH2 + H2O */
    /*eqcon[46] *= 1;  */

    /*reaction 48: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[47] *= 1;  */

    /*reaction 49: OH + CH4 <=> CH3 + H2O */
    /*eqcon[48] *= 1;  */

    /*reaction 50: OH + CO <=> H + CO2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: OH + HCO <=> H2O + CO */
    /*eqcon[50] *= 1;  */

    /*reaction 52: OH + CH2O <=> HCO + H2O */
    /*eqcon[51] *= 1;  */

    /*reaction 53: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[52] *= 1;  */

    /*reaction 54: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[53] *= 1;  */

    /*reaction 55: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[54] *= 1;  */

    /*reaction 56: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: HO2 + CO <=> OH + CO2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: CH2 + O2 <=> OH + HCO */
    /*eqcon[57] *= 1;  */

    /*reaction 59: CH2 + H2 <=> H + CH3 */
    /*eqcon[58] *= 1;  */

    /*reaction 60: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[60] *= 1;  */

    /*reaction 62: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[61] *= 1;  */

    /*reaction 63: CH2(S) + AR <=> CH2 + AR */
    /*eqcon[62] *= 1;  */

    /*reaction 64: CH2(S) + O2 <=> H + OH + CO */
    eqcon[63] *= 1e-06; 

    /*reaction 65: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[64] *= 1;  */

    /*reaction 66: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[65] *= 1;  */

    /*reaction 67: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[68] *= 1;  */

    /*reaction 70: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[69] *= 1;  */

    /*reaction 71: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[70] *= 1;  */

    /*reaction 72: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[71] *= 1;  */

    /*reaction 73: CH3 + O2 <=> O + CH3O */
    /*eqcon[72] *= 1;  */

    /*reaction 74: CH3 + O2 <=> OH + CH2O */
    /*eqcon[73] *= 1;  */

    /*reaction 75: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[74] *= 1e+06; 

    /*reaction 76: 2 CH3 <=> H + C2H5 */
    /*eqcon[75] *= 1;  */

    /*reaction 77: CH3 + HCO <=> CH4 + CO */
    /*eqcon[76] *= 1;  */

    /*reaction 78: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[77] *= 1;  */

    /*reaction 79: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: HCO + H2O <=> H + CO + H2O */
    eqcon[79] *= 1e-06; 

    /*reaction 81: HCO + M <=> H + CO + M */
    eqcon[80] *= 1e-06; 

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    /*eqcon[81] *= 1;  */

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[82] *= 1;  */

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[83] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[21]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: O + H + M <=> OH + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: O + H2 <=> H + OH */
    /*eqcon[1] *= 1;  */

    /*reaction 3: O + HO2 <=> OH + O2 */
    /*eqcon[2] *= 1;  */

    /*reaction 4: O + CH2 <=> H + HCO */
    /*eqcon[3] *= 1;  */

    /*reaction 5: O + CH2(S) <=> H + HCO */
    /*eqcon[4] *= 1;  */

    /*reaction 6: O + CH3 <=> H + CH2O */
    /*eqcon[5] *= 1;  */

    /*reaction 7: O + CH4 <=> OH + CH3 */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + CO + M <=> CO2 + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: O + HCO <=> OH + CO */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + HCO <=> H + CO2 */
    /*eqcon[9] *= 1;  */

    /*reaction 11: O + CH2O <=> OH + HCO */
    /*eqcon[10] *= 1;  */

    /*reaction 12: O + C2H4 <=> CH3 + HCO */
    /*eqcon[11] *= 1;  */

    /*reaction 13: O + C2H5 <=> CH3 + CH2O */
    /*eqcon[12] *= 1;  */

    /*reaction 14: O + C2H6 <=> OH + C2H5 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O2 + CO <=> O + CO2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: O2 + CH2O <=> HO2 + HCO */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + O2 + M <=> HO2 + M */
    eqcon[16] *= 1e+06; 

    /*reaction 18: H + 2 O2 <=> HO2 + O2 */
    eqcon[17] *= 1e+06; 

    /*reaction 19: H + O2 + H2O <=> HO2 + H2O */
    eqcon[18] *= 1e+06; 

    /*reaction 20: H + O2 + N2 <=> HO2 + N2 */
    eqcon[19] *= 1e+06; 

    /*reaction 21: H + O2 + AR <=> HO2 + AR */
    eqcon[20] *= 1e+06; 

    /*reaction 22: H + O2 <=> O + OH */
    /*eqcon[21] *= 1;  */

    /*reaction 23: 2 H + M <=> H2 + M */
    eqcon[22] *= 1e+06; 

    /*reaction 24: 2 H + H2 <=> 2 H2 */
    eqcon[23] *= 1e+06; 

    /*reaction 25: 2 H + H2O <=> H2 + H2O */
    eqcon[24] *= 1e+06; 

    /*reaction 26: 2 H + CO2 <=> H2 + CO2 */
    eqcon[25] *= 1e+06; 

    /*reaction 27: H + OH + M <=> H2O + M */
    eqcon[26] *= 1e+06; 

    /*reaction 28: H + HO2 <=> O2 + H2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H + HO2 <=> 2 OH */
    /*eqcon[28] *= 1;  */

    /*reaction 30: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[29] *= 1e+06; 

    /*reaction 31: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[30] *= 1e+06; 

    /*reaction 32: H + CH4 <=> CH3 + H2 */
    /*eqcon[31] *= 1;  */

    /*reaction 33: H + HCO (+M) <=> CH2O (+M) */
    eqcon[32] *= 1e+06; 

    /*reaction 34: H + HCO <=> H2 + CO */
    /*eqcon[33] *= 1;  */

    /*reaction 35: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[34] *= 1e+06; 

    /*reaction 36: H + CH2O <=> HCO + H2 */
    /*eqcon[35] *= 1;  */

    /*reaction 37: H + CH3O <=> OH + CH3 */
    /*eqcon[36] *= 1;  */

    /*reaction 38: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[37] *= 1e+06; 

    /*reaction 39: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[38] *= 1e+06; 

    /*reaction 40: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[40] *= 1e+06; 

    /*reaction 42: OH + H2 <=> H + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: 2 OH <=> O + H2O */
    /*eqcon[42] *= 1;  */

    /*reaction 44: OH + HO2 <=> O2 + H2O */
    /*eqcon[43] *= 1;  */

    /*reaction 45: OH + CH2 <=> H + CH2O */
    /*eqcon[44] *= 1;  */

    /*reaction 46: OH + CH2(S) <=> H + CH2O */
    /*eqcon[45] *= 1;  */

    /*reaction 47: OH + CH3 <=> CH2 + H2O */
    /*eqcon[46] *= 1;  */

    /*reaction 48: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[47] *= 1;  */

    /*reaction 49: OH + CH4 <=> CH3 + H2O */
    /*eqcon[48] *= 1;  */

    /*reaction 50: OH + CO <=> H + CO2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: OH + HCO <=> H2O + CO */
    /*eqcon[50] *= 1;  */

    /*reaction 52: OH + CH2O <=> HCO + H2O */
    /*eqcon[51] *= 1;  */

    /*reaction 53: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[52] *= 1;  */

    /*reaction 54: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[53] *= 1;  */

    /*reaction 55: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[54] *= 1;  */

    /*reaction 56: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: HO2 + CO <=> OH + CO2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: CH2 + O2 <=> OH + HCO */
    /*eqcon[57] *= 1;  */

    /*reaction 59: CH2 + H2 <=> H + CH3 */
    /*eqcon[58] *= 1;  */

    /*reaction 60: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[60] *= 1;  */

    /*reaction 62: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[61] *= 1;  */

    /*reaction 63: CH2(S) + AR <=> CH2 + AR */
    /*eqcon[62] *= 1;  */

    /*reaction 64: CH2(S) + O2 <=> H + OH + CO */
    eqcon[63] *= 1e-06; 

    /*reaction 65: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[64] *= 1;  */

    /*reaction 66: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[65] *= 1;  */

    /*reaction 67: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[68] *= 1;  */

    /*reaction 70: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[69] *= 1;  */

    /*reaction 71: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[70] *= 1;  */

    /*reaction 72: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[71] *= 1;  */

    /*reaction 73: CH3 + O2 <=> O + CH3O */
    /*eqcon[72] *= 1;  */

    /*reaction 74: CH3 + O2 <=> OH + CH2O */
    /*eqcon[73] *= 1;  */

    /*reaction 75: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[74] *= 1e+06; 

    /*reaction 76: 2 CH3 <=> H + C2H5 */
    /*eqcon[75] *= 1;  */

    /*reaction 77: CH3 + HCO <=> CH4 + CO */
    /*eqcon[76] *= 1;  */

    /*reaction 78: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[77] *= 1;  */

    /*reaction 79: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: HCO + H2O <=> H + CO + H2O */
    eqcon[79] *= 1e-06; 

    /*reaction 81: HCO + M <=> H + CO + M */
    eqcon[80] *= 1e-06; 

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    /*eqcon[81] *= 1;  */

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[82] *= 1;  */

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[83] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[21]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: O + H + M <=> OH + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: O + H2 <=> H + OH */
    /*eqcon[1] *= 1;  */

    /*reaction 3: O + HO2 <=> OH + O2 */
    /*eqcon[2] *= 1;  */

    /*reaction 4: O + CH2 <=> H + HCO */
    /*eqcon[3] *= 1;  */

    /*reaction 5: O + CH2(S) <=> H + HCO */
    /*eqcon[4] *= 1;  */

    /*reaction 6: O + CH3 <=> H + CH2O */
    /*eqcon[5] *= 1;  */

    /*reaction 7: O + CH4 <=> OH + CH3 */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + CO + M <=> CO2 + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: O + HCO <=> OH + CO */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + HCO <=> H + CO2 */
    /*eqcon[9] *= 1;  */

    /*reaction 11: O + CH2O <=> OH + HCO */
    /*eqcon[10] *= 1;  */

    /*reaction 12: O + C2H4 <=> CH3 + HCO */
    /*eqcon[11] *= 1;  */

    /*reaction 13: O + C2H5 <=> CH3 + CH2O */
    /*eqcon[12] *= 1;  */

    /*reaction 14: O + C2H6 <=> OH + C2H5 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O2 + CO <=> O + CO2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: O2 + CH2O <=> HO2 + HCO */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + O2 + M <=> HO2 + M */
    eqcon[16] *= 1e+06; 

    /*reaction 18: H + 2 O2 <=> HO2 + O2 */
    eqcon[17] *= 1e+06; 

    /*reaction 19: H + O2 + H2O <=> HO2 + H2O */
    eqcon[18] *= 1e+06; 

    /*reaction 20: H + O2 + N2 <=> HO2 + N2 */
    eqcon[19] *= 1e+06; 

    /*reaction 21: H + O2 + AR <=> HO2 + AR */
    eqcon[20] *= 1e+06; 

    /*reaction 22: H + O2 <=> O + OH */
    /*eqcon[21] *= 1;  */

    /*reaction 23: 2 H + M <=> H2 + M */
    eqcon[22] *= 1e+06; 

    /*reaction 24: 2 H + H2 <=> 2 H2 */
    eqcon[23] *= 1e+06; 

    /*reaction 25: 2 H + H2O <=> H2 + H2O */
    eqcon[24] *= 1e+06; 

    /*reaction 26: 2 H + CO2 <=> H2 + CO2 */
    eqcon[25] *= 1e+06; 

    /*reaction 27: H + OH + M <=> H2O + M */
    eqcon[26] *= 1e+06; 

    /*reaction 28: H + HO2 <=> O2 + H2 */
    /*eqcon[27] *= 1;  */

    /*reaction 29: H + HO2 <=> 2 OH */
    /*eqcon[28] *= 1;  */

    /*reaction 30: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[29] *= 1e+06; 

    /*reaction 31: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[30] *= 1e+06; 

    /*reaction 32: H + CH4 <=> CH3 + H2 */
    /*eqcon[31] *= 1;  */

    /*reaction 33: H + HCO (+M) <=> CH2O (+M) */
    eqcon[32] *= 1e+06; 

    /*reaction 34: H + HCO <=> H2 + CO */
    /*eqcon[33] *= 1;  */

    /*reaction 35: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[34] *= 1e+06; 

    /*reaction 36: H + CH2O <=> HCO + H2 */
    /*eqcon[35] *= 1;  */

    /*reaction 37: H + CH3O <=> OH + CH3 */
    /*eqcon[36] *= 1;  */

    /*reaction 38: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[37] *= 1e+06; 

    /*reaction 39: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[38] *= 1e+06; 

    /*reaction 40: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[39] *= 1;  */

    /*reaction 41: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[40] *= 1e+06; 

    /*reaction 42: OH + H2 <=> H + H2O */
    /*eqcon[41] *= 1;  */

    /*reaction 43: 2 OH <=> O + H2O */
    /*eqcon[42] *= 1;  */

    /*reaction 44: OH + HO2 <=> O2 + H2O */
    /*eqcon[43] *= 1;  */

    /*reaction 45: OH + CH2 <=> H + CH2O */
    /*eqcon[44] *= 1;  */

    /*reaction 46: OH + CH2(S) <=> H + CH2O */
    /*eqcon[45] *= 1;  */

    /*reaction 47: OH + CH3 <=> CH2 + H2O */
    /*eqcon[46] *= 1;  */

    /*reaction 48: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[47] *= 1;  */

    /*reaction 49: OH + CH4 <=> CH3 + H2O */
    /*eqcon[48] *= 1;  */

    /*reaction 50: OH + CO <=> H + CO2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: OH + HCO <=> H2O + CO */
    /*eqcon[50] *= 1;  */

    /*reaction 52: OH + CH2O <=> HCO + H2O */
    /*eqcon[51] *= 1;  */

    /*reaction 53: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[52] *= 1;  */

    /*reaction 54: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[53] *= 1;  */

    /*reaction 55: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[54] *= 1;  */

    /*reaction 56: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: HO2 + CO <=> OH + CO2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: CH2 + O2 <=> OH + HCO */
    /*eqcon[57] *= 1;  */

    /*reaction 59: CH2 + H2 <=> H + CH3 */
    /*eqcon[58] *= 1;  */

    /*reaction 60: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[60] *= 1;  */

    /*reaction 62: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[61] *= 1;  */

    /*reaction 63: CH2(S) + AR <=> CH2 + AR */
    /*eqcon[62] *= 1;  */

    /*reaction 64: CH2(S) + O2 <=> H + OH + CO */
    eqcon[63] *= 1e-06; 

    /*reaction 65: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[64] *= 1;  */

    /*reaction 66: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[65] *= 1;  */

    /*reaction 67: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[68] *= 1;  */

    /*reaction 70: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[69] *= 1;  */

    /*reaction 71: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[70] *= 1;  */

    /*reaction 72: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[71] *= 1;  */

    /*reaction 73: CH3 + O2 <=> O + CH3O */
    /*eqcon[72] *= 1;  */

    /*reaction 74: CH3 + O2 <=> OH + CH2O */
    /*eqcon[73] *= 1;  */

    /*reaction 75: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[74] *= 1e+06; 

    /*reaction 76: 2 CH3 <=> H + C2H5 */
    /*eqcon[75] *= 1;  */

    /*reaction 77: CH3 + HCO <=> CH4 + CO */
    /*eqcon[76] *= 1;  */

    /*reaction 78: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[77] *= 1;  */

    /*reaction 79: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: HCO + H2O <=> H + CO + H2O */
    eqcon[79] *= 1e-06; 

    /*reaction 81: HCO + M <=> H + CO + M */
    eqcon[80] *= 1e-06; 

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    /*eqcon[81] *= 1;  */

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[82] *= 1;  */

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[83] *= 1;  */
}

static double T_save = -1;
#ifdef BL_USE_OMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[84];
#ifdef BL_USE_OMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[84];
#ifdef BL_USE_OMP
#pragma omp threadprivate(Kc_save)
#endif

/*compute the production rate for each species */
void productionRate(double * wdot, double * sc, double T)
{
    double qdot;

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[21];                /*Gibbs free energy */
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
    for (id = 0; id < 21; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*zero out wdot */
    for (id = 0; id < 21; ++id) {
        wdot[id] = 0.0;
    }

    if (T != T_save)
    {
        T_save = T;

        k_f_save[0] = 1e-12 * 5e+17*exp(-1*tc[0]);
        k_f_save[1] = 1e-06 * 50000*exp(2.67*tc[0]-3165.2328279116868543*invT);
        k_f_save[2] = 1e-06 * 2e+13;
        k_f_save[3] = 1e-06 * 8e+13;
        k_f_save[4] = 1e-06 * 1.5e+13;
        k_f_save[5] = 1e-06 * 8.43e+13;
        k_f_save[6] = 1e-06 * 1.02e+09*exp(1.5*tc[0]-4327.6633259205891591*invT);
        k_f_save[7] = 1e-12 * 6.02e+14*exp(-1509.6499974141590883*invT);
        k_f_save[8] = 1e-06 * 3e+13;
        k_f_save[9] = 1e-06 * 3e+13;
        k_f_save[10] = 1e-06 * 3.9e+13*exp(-1781.3869969487077469*invT);
        k_f_save[11] = 1e-06 * 1.92e+07*exp(1.83*tc[0]-110.7076664770383303*invT);
        k_f_save[12] = 1e-06 * 1.32e+14;
        k_f_save[13] = 1e-06 * 8.98e+07*exp(1.92*tc[0]-2863.3028284288548093*invT);
        k_f_save[14] = 1e-06 * 2.5e+12*exp(-24053.756625465601246*invT);
        k_f_save[15] = 1e-06 * 1e+14*exp(-20128.666632188789663*invT);
        k_f_save[16] = 1e-12 * 2.8e+18*exp(-0.86*tc[0]);
        k_f_save[17] = 1e-12 * 3e+20*exp(-1.72*tc[0]);
        k_f_save[18] = 1e-12 * 9.38e+18*exp(-0.76*tc[0]);
        k_f_save[19] = 1e-12 * 3.75e+20*exp(-1.72*tc[0]);
        k_f_save[20] = 1e-12 * 7e+17*exp(-0.8*tc[0]);
        k_f_save[21] = 1e-06 * 8.3e+13*exp(-7252.8618042434254676*invT);
        k_f_save[22] = 1e-12 * 1e+18*exp(-1*tc[0]);
        k_f_save[23] = 1e-12 * 9e+16*exp(-0.6*tc[0]);
        k_f_save[24] = 1e-12 * 6e+19*exp(-1.25*tc[0]);
        k_f_save[25] = 1e-12 * 5.5e+20*exp(-2*tc[0]);
        k_f_save[26] = 1e-12 * 2.2e+22*exp(-2*tc[0]);
        k_f_save[27] = 1e-06 * 2.8e+13*exp(-537.43539907944057177*invT);
        k_f_save[28] = 1e-06 * 1.34e+14*exp(-319.54258278599701271*invT);
        k_f_save[29] = 1e-06 * 2.5e+16*exp(-0.8*tc[0]);
        k_f_save[30] = 1e-06 * 1.27e+16*exp(-0.63*tc[0]-192.73198300320765952*invT);
        k_f_save[31] = 1e-06 * 6.6e+08*exp(1.62*tc[0]-5454.8686573231616421*invT);
        k_f_save[32] = 1e-06 * 1.09e+12*exp(0.48*tc[0]+130.83633310922709825*invT);
        k_f_save[33] = 1e-06 * 7.34e+13;
        k_f_save[34] = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1308.3633310922712099*invT);
        k_f_save[35] = 1e-06 * 2.3e+10*exp(1.05*tc[0]-1648.0345805104568626*invT);
        k_f_save[36] = 1e-06 * 3.2e+13;
        k_f_save[37] = 1e-06 * 1.08e+12*exp(0.454*tc[0]-915.85433176458980142*invT);
        k_f_save[38] = 1e-06 * 5.21e+17*exp(-0.99*tc[0]-795.08233197145705162*invT);
        k_f_save[39] = 1e-06 * 1.15e+08*exp(1.9*tc[0]-3789.2214935095394139*invT);
        k_f_save[40] = 1e-06 * 4.3e+07*exp(1.5*tc[0]-40056.046598055690993*invT);
        k_f_save[41] = 1e-06 * 2.16e+08*exp(1.51*tc[0]-1726.0331637101885462*invT);
        k_f_save[42] = 1e-06 * 35700*exp(2.4*tc[0]+1061.7871648479585929*invT);
        k_f_save[43] = 1e-06 * 2.9e+13*exp(+251.60833290235981963*invT);
        k_f_save[44] = 1e-06 * 2e+13;
        k_f_save[45] = 1e-06 * 3e+13;
        k_f_save[46] = 1e-06 * 5.6e+07*exp(1.6*tc[0]-2727.4343286615808211*invT);
        k_f_save[47] = 1e-06 * 2.501e+13;
        k_f_save[48] = 1e-06 * 1e+08*exp(1.6*tc[0]-1570.0359973107254064*invT);
        k_f_save[49] = 1e-06 * 4.76e+07*exp(1.228*tc[0]-35.225166606330375885*invT);
        k_f_save[50] = 1e-06 * 5e+13;
        k_f_save[51] = 1e-06 * 3.43e+09*exp(1.18*tc[0]+224.93784961470970529*invT);
        k_f_save[52] = 1e-06 * 3.54e+06*exp(2.12*tc[0]-437.79849925010609013*invT);
        k_f_save[53] = 1e-06 * 2e+13;
        k_f_save[54] = 1e-06 * 1e+12;
        k_f_save[55] = 1e-06 * 2e+13;
        k_f_save[56] = 1e-06 * 1.5e+14*exp(-11875.913312991384373*invT);
        k_f_save[57] = 1e-06 * 1.32e+13*exp(-754.82499870707954415*invT);
        k_f_save[58] = 1e-06 * 500000*exp(2*tc[0]-3638.256493768123164*invT);
        k_f_save[59] = 1e-06 * 4e+13;
        k_f_save[60] = 1e-06 * 2.46e+06*exp(2*tc[0]-4161.6018262050320118*invT);
        k_f_save[61] = 1e-06 * 1.5e+13*exp(-301.92999948283181766*invT);
        k_f_save[62] = 1e-06 * 9e+12*exp(-301.92999948283181766*invT);
        k_f_save[63] = 1e-06 * 2.8e+13;
        k_f_save[64] = 1e-06 * 1.2e+13;
        k_f_save[65] = 1e-06 * 7e+13;
        k_f_save[66] = 1e-06 * 3e+13;
        k_f_save[67] = 1e-06 * 1.2e+13*exp(+286.83349950869023814*invT);
        k_f_save[68] = 1e-06 * 1.6e+13*exp(+286.83349950869023814*invT);
        k_f_save[69] = 1e-06 * 9e+12;
        k_f_save[70] = 1e-06 * 7e+12;
        k_f_save[71] = 1e-06 * 1.4e+13;
        k_f_save[72] = 1e-06 * 2.675e+13*exp(-14492.639975175927248*invT);
        k_f_save[73] = 1e-06 * 3.6e+10*exp(-4498.7569922941938785*invT);
        k_f_save[74] = 1e-06 * 2.12e+16*exp(-0.97*tc[0]-311.99433279892622295*invT);
        k_f_save[75] = 1e-06 * 4.99e+12*exp(0.1*tc[0]-5334.096657530029006*invT);
        k_f_save[76] = 1e-06 * 2.648e+13;
        k_f_save[77] = 1e-06 * 3320*exp(2.81*tc[0]-2948.8496616156576238*invT);
        k_f_save[78] = 1e-06 * 6.14e+06*exp(1.74*tc[0]-5258.6141576593208811*invT);
        k_f_save[79] = 1e-06 * 2.244e+18*exp(-1*tc[0]-8554.6833186802341515*invT);
        k_f_save[80] = 1e-06 * 1.87e+17*exp(-1*tc[0]-8554.6833186802341515*invT);
        k_f_save[81] = 1e-06 * 7.6e+12*exp(-201.28666632188787844*invT);
        k_f_save[82] = 1e-06 * 4.28e-13*exp(7.6*tc[0]+1776.3548302906606295*invT);
        k_f_save[83] = 1e-06 * 8.4e+11*exp(-1949.9645799932889076*invT);

        Kc_save[0] = 1.0 / (refC) * exp((g_RT[2] + g_RT[1]) - (g_RT[4]));
        Kc_save[1] = exp((g_RT[2] + g_RT[0]) - (g_RT[1] + g_RT[4]));
        Kc_save[2] = exp((g_RT[2] + g_RT[6]) - (g_RT[4] + g_RT[3]));
        Kc_save[3] = exp((g_RT[2] + g_RT[7]) - (g_RT[1] + g_RT[13]));
        Kc_save[4] = exp((g_RT[2] + g_RT[8]) - (g_RT[1] + g_RT[13]));
        Kc_save[5] = exp((g_RT[2] + g_RT[9]) - (g_RT[1] + g_RT[14]));
        Kc_save[6] = exp((g_RT[2] + g_RT[10]) - (g_RT[4] + g_RT[9]));
        Kc_save[7] = 1.0 / (refC) * exp((g_RT[2] + g_RT[11]) - (g_RT[12]));
        Kc_save[8] = exp((g_RT[2] + g_RT[13]) - (g_RT[4] + g_RT[11]));
        Kc_save[9] = exp((g_RT[2] + g_RT[13]) - (g_RT[1] + g_RT[12]));
        Kc_save[10] = exp((g_RT[2] + g_RT[14]) - (g_RT[4] + g_RT[13]));
        Kc_save[11] = exp((g_RT[2] + g_RT[16]) - (g_RT[9] + g_RT[13]));
        Kc_save[12] = exp((g_RT[2] + g_RT[17]) - (g_RT[9] + g_RT[14]));
        Kc_save[13] = exp((g_RT[2] + g_RT[18]) - (g_RT[4] + g_RT[17]));
        Kc_save[14] = exp((g_RT[3] + g_RT[11]) - (g_RT[2] + g_RT[12]));
        Kc_save[15] = exp((g_RT[3] + g_RT[14]) - (g_RT[6] + g_RT[13]));
        Kc_save[16] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));
        Kc_save[17] = 1.0 / (refC) * exp((g_RT[1] + 2 * g_RT[3]) - (g_RT[6] + g_RT[3]));
        Kc_save[18] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));
        Kc_save[19] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[19]) - (g_RT[6] + g_RT[19]));
        Kc_save[20] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[20]) - (g_RT[6] + g_RT[20]));
        Kc_save[21] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));
        Kc_save[22] = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));
        Kc_save[23] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));
        Kc_save[24] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));
        Kc_save[25] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[12]) - (g_RT[0] + g_RT[12]));
        Kc_save[26] = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));
        Kc_save[27] = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));
        Kc_save[28] = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));
        Kc_save[29] = 1.0 / (refC) * exp((g_RT[1] + g_RT[7]) - (g_RT[9]));
        Kc_save[30] = 1.0 / (refC) * exp((g_RT[1] + g_RT[9]) - (g_RT[10]));
        Kc_save[31] = exp((g_RT[1] + g_RT[10]) - (g_RT[9] + g_RT[0]));
        Kc_save[32] = 1.0 / (refC) * exp((g_RT[1] + g_RT[13]) - (g_RT[14]));
        Kc_save[33] = exp((g_RT[1] + g_RT[13]) - (g_RT[0] + g_RT[11]));
        Kc_save[34] = 1.0 / (refC) * exp((g_RT[1] + g_RT[14]) - (g_RT[15]));
        Kc_save[35] = exp((g_RT[1] + g_RT[14]) - (g_RT[13] + g_RT[0]));
        Kc_save[36] = exp((g_RT[1] + g_RT[15]) - (g_RT[4] + g_RT[9]));
        Kc_save[37] = 1.0 / (refC) * exp((g_RT[1] + g_RT[16]) - (g_RT[17]));
        Kc_save[38] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[18]));
        Kc_save[39] = exp((g_RT[1] + g_RT[18]) - (g_RT[17] + g_RT[0]));
        Kc_save[40] = 1.0 / (refC) * exp((g_RT[0] + g_RT[11]) - (g_RT[14]));
        Kc_save[41] = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));
        Kc_save[42] = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));
        Kc_save[43] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
        Kc_save[44] = exp((g_RT[4] + g_RT[7]) - (g_RT[1] + g_RT[14]));
        Kc_save[45] = exp((g_RT[4] + g_RT[8]) - (g_RT[1] + g_RT[14]));
        Kc_save[46] = exp((g_RT[4] + g_RT[9]) - (g_RT[7] + g_RT[5]));
        Kc_save[47] = exp((g_RT[4] + g_RT[9]) - (g_RT[8] + g_RT[5]));
        Kc_save[48] = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[5]));
        Kc_save[49] = exp((g_RT[4] + g_RT[11]) - (g_RT[1] + g_RT[12]));
        Kc_save[50] = exp((g_RT[4] + g_RT[13]) - (g_RT[5] + g_RT[11]));
        Kc_save[51] = exp((g_RT[4] + g_RT[14]) - (g_RT[13] + g_RT[5]));
        Kc_save[52] = exp((g_RT[4] + g_RT[18]) - (g_RT[17] + g_RT[5]));
        Kc_save[53] = exp((g_RT[6] + g_RT[7]) - (g_RT[4] + g_RT[14]));
        Kc_save[54] = exp((g_RT[6] + g_RT[9]) - (g_RT[3] + g_RT[10]));
        Kc_save[55] = exp((g_RT[6] + g_RT[9]) - (g_RT[4] + g_RT[15]));
        Kc_save[56] = exp((g_RT[6] + g_RT[11]) - (g_RT[4] + g_RT[12]));
        Kc_save[57] = exp((g_RT[7] + g_RT[3]) - (g_RT[4] + g_RT[13]));
        Kc_save[58] = exp((g_RT[7] + g_RT[0]) - (g_RT[1] + g_RT[9]));
        Kc_save[59] = exp((g_RT[7] + g_RT[9]) - (g_RT[1] + g_RT[16]));
        Kc_save[60] = exp((g_RT[7] + g_RT[10]) - (2 * g_RT[9]));
        Kc_save[61] = exp((g_RT[8] + g_RT[19]) - (g_RT[7] + g_RT[19]));
        Kc_save[62] = exp((g_RT[8] + g_RT[20]) - (g_RT[7] + g_RT[20]));
        Kc_save[63] = refC * exp((g_RT[8] + g_RT[3]) - (g_RT[1] + g_RT[4] + g_RT[11]));
        Kc_save[64] = exp((g_RT[8] + g_RT[3]) - (g_RT[11] + g_RT[5]));
        Kc_save[65] = exp((g_RT[8] + g_RT[0]) - (g_RT[9] + g_RT[1]));
        Kc_save[66] = exp((g_RT[8] + g_RT[5]) - (g_RT[7] + g_RT[5]));
        Kc_save[67] = exp((g_RT[8] + g_RT[9]) - (g_RT[1] + g_RT[16]));
        Kc_save[68] = exp((g_RT[8] + g_RT[10]) - (2 * g_RT[9]));
        Kc_save[69] = exp((g_RT[8] + g_RT[11]) - (g_RT[7] + g_RT[11]));
        Kc_save[70] = exp((g_RT[8] + g_RT[12]) - (g_RT[7] + g_RT[12]));
        Kc_save[71] = exp((g_RT[8] + g_RT[12]) - (g_RT[11] + g_RT[14]));
        Kc_save[72] = exp((g_RT[9] + g_RT[3]) - (g_RT[2] + g_RT[15]));
        Kc_save[73] = exp((g_RT[9] + g_RT[3]) - (g_RT[4] + g_RT[14]));
        Kc_save[74] = 1.0 / (refC) * exp((2 * g_RT[9]) - (g_RT[18]));
        Kc_save[75] = exp((2 * g_RT[9]) - (g_RT[1] + g_RT[17]));
        Kc_save[76] = exp((g_RT[9] + g_RT[13]) - (g_RT[10] + g_RT[11]));
        Kc_save[77] = exp((g_RT[9] + g_RT[14]) - (g_RT[13] + g_RT[10]));
        Kc_save[78] = exp((g_RT[9] + g_RT[18]) - (g_RT[17] + g_RT[10]));
        Kc_save[79] = refC * exp((g_RT[13] + g_RT[5]) - (g_RT[1] + g_RT[11] + g_RT[5]));
        Kc_save[80] = refC * exp((g_RT[13]) - (g_RT[1] + g_RT[11]));
        Kc_save[81] = exp((g_RT[13] + g_RT[3]) - (g_RT[6] + g_RT[11]));
        Kc_save[82] = exp((g_RT[15] + g_RT[3]) - (g_RT[6] + g_RT[14]));
        Kc_save[83] = exp((g_RT[17] + g_RT[3]) - (g_RT[6] + g_RT[16]));
    }

    /*reaction 1: O + H + M <=> OH + M */
    phi_f = sc[2]*sc[1];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = alpha * k_f_save[0];
    q_f = phi_f * k_f;
    phi_r = sc[4];
    Kc = Kc_save[0];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 2: O + H2 <=> H + OH */
    phi_f = sc[2]*sc[0];
    k_f = k_f_save[1];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[4];
    Kc = Kc_save[1];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 3: O + HO2 <=> OH + O2 */
    phi_f = sc[2]*sc[6];
    k_f = k_f_save[2];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = Kc_save[2];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 4: O + CH2 <=> H + HCO */
    phi_f = sc[2]*sc[7];
    k_f = k_f_save[3];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[13];
    Kc = Kc_save[3];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 5: O + CH2(S) <=> H + HCO */
    phi_f = sc[2]*sc[8];
    k_f = k_f_save[4];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[13];
    Kc = Kc_save[4];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 6: O + CH3 <=> H + CH2O */
    phi_f = sc[2]*sc[9];
    k_f = k_f_save[5];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[5];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[9] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 7: O + CH4 <=> OH + CH3 */
    phi_f = sc[2]*sc[10];
    k_f = k_f_save[6];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = Kc_save[6];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 8: O + CO + M <=> CO2 + M */
    phi_f = sc[2]*sc[11];
    alpha = mixture + sc[0] + 5*sc[3] + 5*sc[5] + sc[10] + 0.5*sc[11] + 2.5*sc[12] + 2*sc[18] + -0.5*sc[20];
    k_f = alpha * k_f_save[7];
    q_f = phi_f * k_f;
    phi_r = sc[12];
    Kc = Kc_save[7];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 9: O + HCO <=> OH + CO */
    phi_f = sc[2]*sc[13];
    k_f = k_f_save[8];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[11];
    Kc = Kc_save[8];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[11] += 1 * qdot;

    /*reaction 10: O + HCO <=> H + CO2 */
    phi_f = sc[2]*sc[13];
    k_f = k_f_save[9];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = Kc_save[9];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 11: O + CH2O <=> OH + HCO */
    phi_f = sc[2]*sc[14];
    k_f = k_f_save[10];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[13];
    Kc = Kc_save[10];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 12: O + C2H4 <=> CH3 + HCO */
    phi_f = sc[2]*sc[16];
    k_f = k_f_save[11];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[13];
    Kc = Kc_save[11];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[16] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 13: O + C2H5 <=> CH3 + CH2O */
    phi_f = sc[2]*sc[17];
    k_f = k_f_save[12];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[14];
    Kc = Kc_save[12];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 14: O + C2H6 <=> OH + C2H5 */
    phi_f = sc[2]*sc[18];
    k_f = k_f_save[13];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[13];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 15: O2 + CO <=> O + CO2 */
    phi_f = sc[3]*sc[11];
    k_f = k_f_save[14];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[12];
    Kc = Kc_save[14];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 16: O2 + CH2O <=> HO2 + HCO */
    phi_f = sc[3]*sc[14];
    k_f = k_f_save[15];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = Kc_save[15];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 17: H + O2 + M <=> HO2 + M */
    phi_f = sc[1]*sc[3];
    alpha = mixture + -1*sc[3] + -1*sc[5] + -0.25*sc[11] + 0.5*sc[12] + 0.5*sc[18] + -1*sc[19] + -1*sc[20];
    k_f = alpha * k_f_save[16];
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = Kc_save[16];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 18: H + 2 O2 <=> HO2 + O2 */
    phi_f = sc[1]*sc[3]*sc[3];
    k_f = k_f_save[17];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = Kc_save[17];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 2 * qdot;
    wdot[6] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 19: H + O2 + H2O <=> HO2 + H2O */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = k_f_save[18];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = Kc_save[18];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 20: H + O2 + N2 <=> HO2 + N2 */
    phi_f = sc[1]*sc[3]*sc[19];
    k_f = k_f_save[19];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[19];
    Kc = Kc_save[19];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[19] += 1 * qdot;

    /*reaction 21: H + O2 + AR <=> HO2 + AR */
    phi_f = sc[1]*sc[3]*sc[20];
    k_f = k_f_save[20];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[20];
    Kc = Kc_save[20];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[20] += 1 * qdot;

    /*reaction 22: H + O2 <=> O + OH */
    phi_f = sc[1]*sc[3];
    k_f = k_f_save[21];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = Kc_save[21];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 23: 2 H + M <=> H2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[0] + -1*sc[5] + sc[10] + -1*sc[12] + 2*sc[18] + -0.37*sc[20];
    k_f = alpha * k_f_save[22];
    q_f = phi_f * k_f;
    phi_r = sc[0];
    Kc = Kc_save[22];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 24: 2 H + H2 <=> 2 H2 */
    phi_f = sc[1]*sc[1]*sc[0];
    k_f = k_f_save[23];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[0];
    Kc = Kc_save[23];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[0] += 2 * qdot;

    /*reaction 25: 2 H + H2O <=> H2 + H2O */
    phi_f = sc[1]*sc[1]*sc[5];
    k_f = k_f_save[24];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[5];
    Kc = Kc_save[24];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 26: 2 H + CO2 <=> H2 + CO2 */
    phi_f = sc[1]*sc[1]*sc[12];
    k_f = k_f_save[25];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[12];
    Kc = Kc_save[25];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 27: H + OH + M <=> H2O + M */
    phi_f = sc[1]*sc[4];
    alpha = mixture + -0.27*sc[0] + 2.65*sc[5] + sc[10] + 2*sc[18] + -0.62*sc[20];
    k_f = alpha * k_f_save[26];
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = Kc_save[26];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 28: H + HO2 <=> O2 + H2 */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[27];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = Kc_save[27];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 29: H + HO2 <=> 2 OH */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[28];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = Kc_save[28];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[4] += 2 * qdot;

    /*reaction 30: H + CH2 (+M) <=> CH3 (+M) */
    phi_f = sc[1]*sc[7];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[29];
    redP = 1e-12 * alpha / k_f * 3.2e+27*exp(-3.14*tc[0]-618.95649893980521483*invT);
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
    phi_r = sc[9];
    Kc = Kc_save[29];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 31: H + CH3 (+M) <=> CH4 (+M) */
    phi_f = sc[1]*sc[9];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[30];
    redP = 1e-12 * alpha / k_f * 2.477e+33*exp(-4.76*tc[0]-1227.8486645635159675*invT);
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
    phi_r = sc[10];
    Kc = Kc_save[30];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[9] -= 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 32: H + CH4 <=> CH3 + H2 */
    phi_f = sc[1]*sc[10];
    k_f = k_f_save[31];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = Kc_save[31];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 33: H + HCO (+M) <=> CH2O (+M) */
    phi_f = sc[1]*sc[13];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[32];
    redP = 1e-12 * alpha / k_f * 1.35e+24*exp(-2.57*tc[0]-717.08374877172548167*invT);
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
    phi_r = sc[14];
    Kc = Kc_save[32];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 34: H + HCO <=> H2 + CO */
    phi_f = sc[1]*sc[13];
    k_f = k_f_save[33];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[11];
    Kc = Kc_save[33];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[11] += 1 * qdot;

    /*reaction 35: H + CH2O (+M) <=> CH3O (+M) */
    phi_f = sc[1]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18];
    k_f = k_f_save[34];
    redP = 1e-12 * alpha / k_f * 2.2e+30*exp(-4.8*tc[0]-2797.8846618742413739*invT);
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
    phi_r = sc[15];
    Kc = Kc_save[34];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 36: H + CH2O <=> HCO + H2 */
    phi_f = sc[1]*sc[14];
    k_f = k_f_save[35];
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[0];
    Kc = Kc_save[35];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 37: H + CH3O <=> OH + CH3 */
    phi_f = sc[1]*sc[15];
    k_f = k_f_save[36];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = Kc_save[36];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 38: H + C2H4 (+M) <=> C2H5 (+M) */
    phi_f = sc[1]*sc[16];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[37];
    redP = 1e-12 * alpha / k_f * 1.2e+42*exp(-7.62*tc[0]-3507.4201606588962932*invT);
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
    phi_r = sc[17];
    Kc = Kc_save[37];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[16] -= 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 39: H + C2H5 (+M) <=> C2H6 (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[38];
    redP = 1e-12 * alpha / k_f * 1.99e+41*exp(-7.08*tc[0]-3364.0034109045514015*invT);
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
    phi_r = sc[18];
    Kc = Kc_save[38];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[18] += 1 * qdot;

    /*reaction 40: H + C2H6 <=> C2H5 + H2 */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[39];
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[0];
    Kc = Kc_save[39];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 41: H2 + CO (+M) <=> CH2O (+M) */
    phi_f = sc[0]*sc[11];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[40];
    redP = 1e-12 * alpha / k_f * 5.07e+27*exp(-3.42*tc[0]-42446.325760628104035*invT);
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
    phi_r = sc[14];
    Kc = Kc_save[40];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 42: OH + H2 <=> H + H2O */
    phi_f = sc[4]*sc[0];
    k_f = k_f_save[41];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    Kc = Kc_save[41];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 43: 2 OH <=> O + H2O */
    phi_f = sc[4]*sc[4];
    k_f = k_f_save[42];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = Kc_save[42];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 2 * qdot;
    wdot[2] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 44: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = k_f_save[43];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = Kc_save[43];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 45: OH + CH2 <=> H + CH2O */
    phi_f = sc[4]*sc[7];
    k_f = k_f_save[44];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[44];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 46: OH + CH2(S) <=> H + CH2O */
    phi_f = sc[4]*sc[8];
    k_f = k_f_save[45];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[45];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 47: OH + CH3 <=> CH2 + H2O */
    phi_f = sc[4]*sc[9];
    k_f = k_f_save[46];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[5];
    Kc = Kc_save[46];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[9] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 48: OH + CH3 <=> CH2(S) + H2O */
    phi_f = sc[4]*sc[9];
    k_f = k_f_save[47];
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[5];
    Kc = Kc_save[47];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[9] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 49: OH + CH4 <=> CH3 + H2O */
    phi_f = sc[4]*sc[10];
    k_f = k_f_save[48];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[5];
    Kc = Kc_save[48];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 50: OH + CO <=> H + CO2 */
    phi_f = sc[4]*sc[11];
    k_f = k_f_save[49];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = Kc_save[49];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 51: OH + HCO <=> H2O + CO */
    phi_f = sc[4]*sc[13];
    k_f = k_f_save[50];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[11];
    Kc = Kc_save[50];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[11] += 1 * qdot;

    /*reaction 52: OH + CH2O <=> HCO + H2O */
    phi_f = sc[4]*sc[14];
    k_f = k_f_save[51];
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[5];
    Kc = Kc_save[51];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 53: OH + C2H6 <=> C2H5 + H2O */
    phi_f = sc[4]*sc[18];
    k_f = k_f_save[52];
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[5];
    Kc = Kc_save[52];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 54: HO2 + CH2 <=> OH + CH2O */
    phi_f = sc[6]*sc[7];
    k_f = k_f_save[53];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[14];
    Kc = Kc_save[53];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 55: HO2 + CH3 <=> O2 + CH4 */
    phi_f = sc[6]*sc[9];
    k_f = k_f_save[54];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[10];
    Kc = Kc_save[54];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[9] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 56: HO2 + CH3 <=> OH + CH3O */
    phi_f = sc[6]*sc[9];
    k_f = k_f_save[55];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[15];
    Kc = Kc_save[55];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[9] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 57: HO2 + CO <=> OH + CO2 */
    phi_f = sc[6]*sc[11];
    k_f = k_f_save[56];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = Kc_save[56];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 58: CH2 + O2 <=> OH + HCO */
    phi_f = sc[7]*sc[3];
    k_f = k_f_save[57];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[13];
    Kc = Kc_save[57];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 59: CH2 + H2 <=> H + CH3 */
    phi_f = sc[7]*sc[0];
    k_f = k_f_save[58];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[9];
    Kc = Kc_save[58];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 60: CH2 + CH3 <=> H + C2H4 */
    phi_f = sc[7]*sc[9];
    k_f = k_f_save[59];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_save[59];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[9] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 61: CH2 + CH4 <=> 2 CH3 */
    phi_f = sc[7]*sc[10];
    k_f = k_f_save[60];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[9];
    Kc = Kc_save[60];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[9] += 2 * qdot;

    /*reaction 62: CH2(S) + N2 <=> CH2 + N2 */
    phi_f = sc[8]*sc[19];
    k_f = k_f_save[61];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[19];
    Kc = Kc_save[61];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[19] += 1 * qdot;

    /*reaction 63: CH2(S) + AR <=> CH2 + AR */
    phi_f = sc[8]*sc[20];
    k_f = k_f_save[62];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[20];
    Kc = Kc_save[62];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[20] += 1 * qdot;

    /*reaction 64: CH2(S) + O2 <=> H + OH + CO */
    phi_f = sc[8]*sc[3];
    k_f = k_f_save[63];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[4]*sc[11];
    Kc = Kc_save[63];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[11] += 1 * qdot;

    /*reaction 65: CH2(S) + O2 <=> CO + H2O */
    phi_f = sc[8]*sc[3];
    k_f = k_f_save[64];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[64];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 66: CH2(S) + H2 <=> CH3 + H */
    phi_f = sc[8]*sc[0];
    k_f = k_f_save[65];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[1];
    Kc = Kc_save[65];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 67: CH2(S) + H2O <=> CH2 + H2O */
    phi_f = sc[8]*sc[5];
    k_f = k_f_save[66];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[5];
    Kc = Kc_save[66];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 68: CH2(S) + CH3 <=> H + C2H4 */
    phi_f = sc[8]*sc[9];
    k_f = k_f_save[67];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_save[67];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[9] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 69: CH2(S) + CH4 <=> 2 CH3 */
    phi_f = sc[8]*sc[10];
    k_f = k_f_save[68];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[9];
    Kc = Kc_save[68];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[9] += 2 * qdot;

    /*reaction 70: CH2(S) + CO <=> CH2 + CO */
    phi_f = sc[8]*sc[11];
    k_f = k_f_save[69];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[11];
    Kc = Kc_save[69];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[11] += 1 * qdot;

    /*reaction 71: CH2(S) + CO2 <=> CH2 + CO2 */
    phi_f = sc[8]*sc[12];
    k_f = k_f_save[70];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[12];
    Kc = Kc_save[70];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 72: CH2(S) + CO2 <=> CO + CH2O */
    phi_f = sc[8]*sc[12];
    k_f = k_f_save[71];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[14];
    Kc = Kc_save[71];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 73: CH3 + O2 <=> O + CH3O */
    phi_f = sc[9]*sc[3];
    k_f = k_f_save[72];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[15];
    Kc = Kc_save[72];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 74: CH3 + O2 <=> OH + CH2O */
    phi_f = sc[9]*sc[3];
    k_f = k_f_save[73];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[14];
    Kc = Kc_save[73];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 75: 2 CH3 (+M) <=> C2H6 (+M) */
    phi_f = sc[9]*sc[9];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[74];
    redP = 1e-12 * alpha / k_f * 1.77e+50*exp(-9.67*tc[0]-3130.0076613053565779*invT);
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
    phi_r = sc[18];
    Kc = Kc_save[74];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 2 * qdot;
    wdot[18] += 1 * qdot;

    /*reaction 76: 2 CH3 <=> H + C2H5 */
    phi_f = sc[9]*sc[9];
    k_f = k_f_save[75];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[75];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 2 * qdot;
    wdot[1] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 77: CH3 + HCO <=> CH4 + CO */
    phi_f = sc[9]*sc[13];
    k_f = k_f_save[76];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[11];
    Kc = Kc_save[76];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[11] += 1 * qdot;

    /*reaction 78: CH3 + CH2O <=> HCO + CH4 */
    phi_f = sc[9]*sc[14];
    k_f = k_f_save[77];
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[10];
    Kc = Kc_save[77];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 79: CH3 + C2H6 <=> C2H5 + CH4 */
    phi_f = sc[9]*sc[18];
    k_f = k_f_save[78];
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[10];
    Kc = Kc_save[78];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 80: HCO + H2O <=> H + CO + H2O */
    phi_f = sc[13]*sc[5];
    k_f = k_f_save[79];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[11]*sc[5];
    Kc = Kc_save[79];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 81: HCO + M <=> H + CO + M */
    phi_f = sc[13];
    alpha = mixture + sc[0] + -1*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18];
    k_f = alpha * k_f_save[80];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[11];
    Kc = Kc_save[80];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[11] += 1 * qdot;

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    phi_f = sc[13]*sc[3];
    k_f = k_f_save[81];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[11];
    Kc = Kc_save[81];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[11] += 1 * qdot;

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    phi_f = sc[15]*sc[3];
    k_f = k_f_save[82];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[14];
    Kc = Kc_save[82];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    phi_f = sc[17]*sc[3];
    k_f = k_f_save[83];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[16];
    Kc = Kc_save[83];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[17] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[16] += 1 * qdot;

    return;
}


/*compute the progress rate for each reaction */
void progressRate(double * qdot, double * sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[21];                /*Gibbs free energy */
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
    for (id = 0; id < 21; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    if (T != T_save)
    {
        T_save = T;

        k_f_save[0] = 1e-12 * 5e+17*exp(-1*tc[0]);
        k_f_save[1] = 1e-06 * 50000*exp(2.67*tc[0]-3165.2328279116868543*invT);
        k_f_save[2] = 1e-06 * 2e+13;
        k_f_save[3] = 1e-06 * 8e+13;
        k_f_save[4] = 1e-06 * 1.5e+13;
        k_f_save[5] = 1e-06 * 8.43e+13;
        k_f_save[6] = 1e-06 * 1.02e+09*exp(1.5*tc[0]-4327.6633259205891591*invT);
        k_f_save[7] = 1e-12 * 6.02e+14*exp(-1509.6499974141590883*invT);
        k_f_save[8] = 1e-06 * 3e+13;
        k_f_save[9] = 1e-06 * 3e+13;
        k_f_save[10] = 1e-06 * 3.9e+13*exp(-1781.3869969487077469*invT);
        k_f_save[11] = 1e-06 * 1.92e+07*exp(1.83*tc[0]-110.7076664770383303*invT);
        k_f_save[12] = 1e-06 * 1.32e+14;
        k_f_save[13] = 1e-06 * 8.98e+07*exp(1.92*tc[0]-2863.3028284288548093*invT);
        k_f_save[14] = 1e-06 * 2.5e+12*exp(-24053.756625465601246*invT);
        k_f_save[15] = 1e-06 * 1e+14*exp(-20128.666632188789663*invT);
        k_f_save[16] = 1e-12 * 2.8e+18*exp(-0.86*tc[0]);
        k_f_save[17] = 1e-12 * 3e+20*exp(-1.72*tc[0]);
        k_f_save[18] = 1e-12 * 9.38e+18*exp(-0.76*tc[0]);
        k_f_save[19] = 1e-12 * 3.75e+20*exp(-1.72*tc[0]);
        k_f_save[20] = 1e-12 * 7e+17*exp(-0.8*tc[0]);
        k_f_save[21] = 1e-06 * 8.3e+13*exp(-7252.8618042434254676*invT);
        k_f_save[22] = 1e-12 * 1e+18*exp(-1*tc[0]);
        k_f_save[23] = 1e-12 * 9e+16*exp(-0.6*tc[0]);
        k_f_save[24] = 1e-12 * 6e+19*exp(-1.25*tc[0]);
        k_f_save[25] = 1e-12 * 5.5e+20*exp(-2*tc[0]);
        k_f_save[26] = 1e-12 * 2.2e+22*exp(-2*tc[0]);
        k_f_save[27] = 1e-06 * 2.8e+13*exp(-537.43539907944057177*invT);
        k_f_save[28] = 1e-06 * 1.34e+14*exp(-319.54258278599701271*invT);
        k_f_save[29] = 1e-06 * 2.5e+16*exp(-0.8*tc[0]);
        k_f_save[30] = 1e-06 * 1.27e+16*exp(-0.63*tc[0]-192.73198300320765952*invT);
        k_f_save[31] = 1e-06 * 6.6e+08*exp(1.62*tc[0]-5454.8686573231616421*invT);
        k_f_save[32] = 1e-06 * 1.09e+12*exp(0.48*tc[0]+130.83633310922709825*invT);
        k_f_save[33] = 1e-06 * 7.34e+13;
        k_f_save[34] = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1308.3633310922712099*invT);
        k_f_save[35] = 1e-06 * 2.3e+10*exp(1.05*tc[0]-1648.0345805104568626*invT);
        k_f_save[36] = 1e-06 * 3.2e+13;
        k_f_save[37] = 1e-06 * 1.08e+12*exp(0.454*tc[0]-915.85433176458980142*invT);
        k_f_save[38] = 1e-06 * 5.21e+17*exp(-0.99*tc[0]-795.08233197145705162*invT);
        k_f_save[39] = 1e-06 * 1.15e+08*exp(1.9*tc[0]-3789.2214935095394139*invT);
        k_f_save[40] = 1e-06 * 4.3e+07*exp(1.5*tc[0]-40056.046598055690993*invT);
        k_f_save[41] = 1e-06 * 2.16e+08*exp(1.51*tc[0]-1726.0331637101885462*invT);
        k_f_save[42] = 1e-06 * 35700*exp(2.4*tc[0]+1061.7871648479585929*invT);
        k_f_save[43] = 1e-06 * 2.9e+13*exp(+251.60833290235981963*invT);
        k_f_save[44] = 1e-06 * 2e+13;
        k_f_save[45] = 1e-06 * 3e+13;
        k_f_save[46] = 1e-06 * 5.6e+07*exp(1.6*tc[0]-2727.4343286615808211*invT);
        k_f_save[47] = 1e-06 * 2.501e+13;
        k_f_save[48] = 1e-06 * 1e+08*exp(1.6*tc[0]-1570.0359973107254064*invT);
        k_f_save[49] = 1e-06 * 4.76e+07*exp(1.228*tc[0]-35.225166606330375885*invT);
        k_f_save[50] = 1e-06 * 5e+13;
        k_f_save[51] = 1e-06 * 3.43e+09*exp(1.18*tc[0]+224.93784961470970529*invT);
        k_f_save[52] = 1e-06 * 3.54e+06*exp(2.12*tc[0]-437.79849925010609013*invT);
        k_f_save[53] = 1e-06 * 2e+13;
        k_f_save[54] = 1e-06 * 1e+12;
        k_f_save[55] = 1e-06 * 2e+13;
        k_f_save[56] = 1e-06 * 1.5e+14*exp(-11875.913312991384373*invT);
        k_f_save[57] = 1e-06 * 1.32e+13*exp(-754.82499870707954415*invT);
        k_f_save[58] = 1e-06 * 500000*exp(2*tc[0]-3638.256493768123164*invT);
        k_f_save[59] = 1e-06 * 4e+13;
        k_f_save[60] = 1e-06 * 2.46e+06*exp(2*tc[0]-4161.6018262050320118*invT);
        k_f_save[61] = 1e-06 * 1.5e+13*exp(-301.92999948283181766*invT);
        k_f_save[62] = 1e-06 * 9e+12*exp(-301.92999948283181766*invT);
        k_f_save[63] = 1e-06 * 2.8e+13;
        k_f_save[64] = 1e-06 * 1.2e+13;
        k_f_save[65] = 1e-06 * 7e+13;
        k_f_save[66] = 1e-06 * 3e+13;
        k_f_save[67] = 1e-06 * 1.2e+13*exp(+286.83349950869023814*invT);
        k_f_save[68] = 1e-06 * 1.6e+13*exp(+286.83349950869023814*invT);
        k_f_save[69] = 1e-06 * 9e+12;
        k_f_save[70] = 1e-06 * 7e+12;
        k_f_save[71] = 1e-06 * 1.4e+13;
        k_f_save[72] = 1e-06 * 2.675e+13*exp(-14492.639975175927248*invT);
        k_f_save[73] = 1e-06 * 3.6e+10*exp(-4498.7569922941938785*invT);
        k_f_save[74] = 1e-06 * 2.12e+16*exp(-0.97*tc[0]-311.99433279892622295*invT);
        k_f_save[75] = 1e-06 * 4.99e+12*exp(0.1*tc[0]-5334.096657530029006*invT);
        k_f_save[76] = 1e-06 * 2.648e+13;
        k_f_save[77] = 1e-06 * 3320*exp(2.81*tc[0]-2948.8496616156576238*invT);
        k_f_save[78] = 1e-06 * 6.14e+06*exp(1.74*tc[0]-5258.6141576593208811*invT);
        k_f_save[79] = 1e-06 * 2.244e+18*exp(-1*tc[0]-8554.6833186802341515*invT);
        k_f_save[80] = 1e-06 * 1.87e+17*exp(-1*tc[0]-8554.6833186802341515*invT);
        k_f_save[81] = 1e-06 * 7.6e+12*exp(-201.28666632188787844*invT);
        k_f_save[82] = 1e-06 * 4.28e-13*exp(7.6*tc[0]+1776.3548302906606295*invT);
        k_f_save[83] = 1e-06 * 8.4e+11*exp(-1949.9645799932889076*invT);

        Kc_save[0] = 1.0 / (refC) * exp((g_RT[2] + g_RT[1]) - (g_RT[4]));
        Kc_save[1] = exp((g_RT[2] + g_RT[0]) - (g_RT[1] + g_RT[4]));
        Kc_save[2] = exp((g_RT[2] + g_RT[6]) - (g_RT[4] + g_RT[3]));
        Kc_save[3] = exp((g_RT[2] + g_RT[7]) - (g_RT[1] + g_RT[13]));
        Kc_save[4] = exp((g_RT[2] + g_RT[8]) - (g_RT[1] + g_RT[13]));
        Kc_save[5] = exp((g_RT[2] + g_RT[9]) - (g_RT[1] + g_RT[14]));
        Kc_save[6] = exp((g_RT[2] + g_RT[10]) - (g_RT[4] + g_RT[9]));
        Kc_save[7] = 1.0 / (refC) * exp((g_RT[2] + g_RT[11]) - (g_RT[12]));
        Kc_save[8] = exp((g_RT[2] + g_RT[13]) - (g_RT[4] + g_RT[11]));
        Kc_save[9] = exp((g_RT[2] + g_RT[13]) - (g_RT[1] + g_RT[12]));
        Kc_save[10] = exp((g_RT[2] + g_RT[14]) - (g_RT[4] + g_RT[13]));
        Kc_save[11] = exp((g_RT[2] + g_RT[16]) - (g_RT[9] + g_RT[13]));
        Kc_save[12] = exp((g_RT[2] + g_RT[17]) - (g_RT[9] + g_RT[14]));
        Kc_save[13] = exp((g_RT[2] + g_RT[18]) - (g_RT[4] + g_RT[17]));
        Kc_save[14] = exp((g_RT[3] + g_RT[11]) - (g_RT[2] + g_RT[12]));
        Kc_save[15] = exp((g_RT[3] + g_RT[14]) - (g_RT[6] + g_RT[13]));
        Kc_save[16] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));
        Kc_save[17] = 1.0 / (refC) * exp((g_RT[1] + 2 * g_RT[3]) - (g_RT[6] + g_RT[3]));
        Kc_save[18] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));
        Kc_save[19] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[19]) - (g_RT[6] + g_RT[19]));
        Kc_save[20] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[20]) - (g_RT[6] + g_RT[20]));
        Kc_save[21] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));
        Kc_save[22] = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));
        Kc_save[23] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));
        Kc_save[24] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));
        Kc_save[25] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[12]) - (g_RT[0] + g_RT[12]));
        Kc_save[26] = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));
        Kc_save[27] = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));
        Kc_save[28] = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));
        Kc_save[29] = 1.0 / (refC) * exp((g_RT[1] + g_RT[7]) - (g_RT[9]));
        Kc_save[30] = 1.0 / (refC) * exp((g_RT[1] + g_RT[9]) - (g_RT[10]));
        Kc_save[31] = exp((g_RT[1] + g_RT[10]) - (g_RT[9] + g_RT[0]));
        Kc_save[32] = 1.0 / (refC) * exp((g_RT[1] + g_RT[13]) - (g_RT[14]));
        Kc_save[33] = exp((g_RT[1] + g_RT[13]) - (g_RT[0] + g_RT[11]));
        Kc_save[34] = 1.0 / (refC) * exp((g_RT[1] + g_RT[14]) - (g_RT[15]));
        Kc_save[35] = exp((g_RT[1] + g_RT[14]) - (g_RT[13] + g_RT[0]));
        Kc_save[36] = exp((g_RT[1] + g_RT[15]) - (g_RT[4] + g_RT[9]));
        Kc_save[37] = 1.0 / (refC) * exp((g_RT[1] + g_RT[16]) - (g_RT[17]));
        Kc_save[38] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[18]));
        Kc_save[39] = exp((g_RT[1] + g_RT[18]) - (g_RT[17] + g_RT[0]));
        Kc_save[40] = 1.0 / (refC) * exp((g_RT[0] + g_RT[11]) - (g_RT[14]));
        Kc_save[41] = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));
        Kc_save[42] = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));
        Kc_save[43] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
        Kc_save[44] = exp((g_RT[4] + g_RT[7]) - (g_RT[1] + g_RT[14]));
        Kc_save[45] = exp((g_RT[4] + g_RT[8]) - (g_RT[1] + g_RT[14]));
        Kc_save[46] = exp((g_RT[4] + g_RT[9]) - (g_RT[7] + g_RT[5]));
        Kc_save[47] = exp((g_RT[4] + g_RT[9]) - (g_RT[8] + g_RT[5]));
        Kc_save[48] = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[5]));
        Kc_save[49] = exp((g_RT[4] + g_RT[11]) - (g_RT[1] + g_RT[12]));
        Kc_save[50] = exp((g_RT[4] + g_RT[13]) - (g_RT[5] + g_RT[11]));
        Kc_save[51] = exp((g_RT[4] + g_RT[14]) - (g_RT[13] + g_RT[5]));
        Kc_save[52] = exp((g_RT[4] + g_RT[18]) - (g_RT[17] + g_RT[5]));
        Kc_save[53] = exp((g_RT[6] + g_RT[7]) - (g_RT[4] + g_RT[14]));
        Kc_save[54] = exp((g_RT[6] + g_RT[9]) - (g_RT[3] + g_RT[10]));
        Kc_save[55] = exp((g_RT[6] + g_RT[9]) - (g_RT[4] + g_RT[15]));
        Kc_save[56] = exp((g_RT[6] + g_RT[11]) - (g_RT[4] + g_RT[12]));
        Kc_save[57] = exp((g_RT[7] + g_RT[3]) - (g_RT[4] + g_RT[13]));
        Kc_save[58] = exp((g_RT[7] + g_RT[0]) - (g_RT[1] + g_RT[9]));
        Kc_save[59] = exp((g_RT[7] + g_RT[9]) - (g_RT[1] + g_RT[16]));
        Kc_save[60] = exp((g_RT[7] + g_RT[10]) - (2 * g_RT[9]));
        Kc_save[61] = exp((g_RT[8] + g_RT[19]) - (g_RT[7] + g_RT[19]));
        Kc_save[62] = exp((g_RT[8] + g_RT[20]) - (g_RT[7] + g_RT[20]));
        Kc_save[63] = refC * exp((g_RT[8] + g_RT[3]) - (g_RT[1] + g_RT[4] + g_RT[11]));
        Kc_save[64] = exp((g_RT[8] + g_RT[3]) - (g_RT[11] + g_RT[5]));
        Kc_save[65] = exp((g_RT[8] + g_RT[0]) - (g_RT[9] + g_RT[1]));
        Kc_save[66] = exp((g_RT[8] + g_RT[5]) - (g_RT[7] + g_RT[5]));
        Kc_save[67] = exp((g_RT[8] + g_RT[9]) - (g_RT[1] + g_RT[16]));
        Kc_save[68] = exp((g_RT[8] + g_RT[10]) - (2 * g_RT[9]));
        Kc_save[69] = exp((g_RT[8] + g_RT[11]) - (g_RT[7] + g_RT[11]));
        Kc_save[70] = exp((g_RT[8] + g_RT[12]) - (g_RT[7] + g_RT[12]));
        Kc_save[71] = exp((g_RT[8] + g_RT[12]) - (g_RT[11] + g_RT[14]));
        Kc_save[72] = exp((g_RT[9] + g_RT[3]) - (g_RT[2] + g_RT[15]));
        Kc_save[73] = exp((g_RT[9] + g_RT[3]) - (g_RT[4] + g_RT[14]));
        Kc_save[74] = 1.0 / (refC) * exp((2 * g_RT[9]) - (g_RT[18]));
        Kc_save[75] = exp((2 * g_RT[9]) - (g_RT[1] + g_RT[17]));
        Kc_save[76] = exp((g_RT[9] + g_RT[13]) - (g_RT[10] + g_RT[11]));
        Kc_save[77] = exp((g_RT[9] + g_RT[14]) - (g_RT[13] + g_RT[10]));
        Kc_save[78] = exp((g_RT[9] + g_RT[18]) - (g_RT[17] + g_RT[10]));
        Kc_save[79] = refC * exp((g_RT[13] + g_RT[5]) - (g_RT[1] + g_RT[11] + g_RT[5]));
        Kc_save[80] = refC * exp((g_RT[13]) - (g_RT[1] + g_RT[11]));
        Kc_save[81] = exp((g_RT[13] + g_RT[3]) - (g_RT[6] + g_RT[11]));
        Kc_save[82] = exp((g_RT[15] + g_RT[3]) - (g_RT[6] + g_RT[14]));
        Kc_save[83] = exp((g_RT[17] + g_RT[3]) - (g_RT[6] + g_RT[16]));
    }

    /*reaction 1: O + H + M <=> OH + M */
    phi_f = sc[2]*sc[1];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = alpha * k_f_save[0];
    q_f = phi_f * k_f;
    phi_r = sc[4];
    Kc = Kc_save[0];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[0] = q_f - q_r;

    /*reaction 2: O + H2 <=> H + OH */
    phi_f = sc[2]*sc[0];
    k_f = k_f_save[1];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[4];
    Kc = Kc_save[1];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[1] = q_f - q_r;

    /*reaction 3: O + HO2 <=> OH + O2 */
    phi_f = sc[2]*sc[6];
    k_f = k_f_save[2];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = Kc_save[2];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[2] = q_f - q_r;

    /*reaction 4: O + CH2 <=> H + HCO */
    phi_f = sc[2]*sc[7];
    k_f = k_f_save[3];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[13];
    Kc = Kc_save[3];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[3] = q_f - q_r;

    /*reaction 5: O + CH2(S) <=> H + HCO */
    phi_f = sc[2]*sc[8];
    k_f = k_f_save[4];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[13];
    Kc = Kc_save[4];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[4] = q_f - q_r;

    /*reaction 6: O + CH3 <=> H + CH2O */
    phi_f = sc[2]*sc[9];
    k_f = k_f_save[5];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[5];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[5] = q_f - q_r;

    /*reaction 7: O + CH4 <=> OH + CH3 */
    phi_f = sc[2]*sc[10];
    k_f = k_f_save[6];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = Kc_save[6];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[6] = q_f - q_r;

    /*reaction 8: O + CO + M <=> CO2 + M */
    phi_f = sc[2]*sc[11];
    alpha = mixture + sc[0] + 5*sc[3] + 5*sc[5] + sc[10] + 0.5*sc[11] + 2.5*sc[12] + 2*sc[18] + -0.5*sc[20];
    k_f = alpha * k_f_save[7];
    q_f = phi_f * k_f;
    phi_r = sc[12];
    Kc = Kc_save[7];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[7] = q_f - q_r;

    /*reaction 9: O + HCO <=> OH + CO */
    phi_f = sc[2]*sc[13];
    k_f = k_f_save[8];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[11];
    Kc = Kc_save[8];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[8] = q_f - q_r;

    /*reaction 10: O + HCO <=> H + CO2 */
    phi_f = sc[2]*sc[13];
    k_f = k_f_save[9];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = Kc_save[9];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[9] = q_f - q_r;

    /*reaction 11: O + CH2O <=> OH + HCO */
    phi_f = sc[2]*sc[14];
    k_f = k_f_save[10];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[13];
    Kc = Kc_save[10];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[10] = q_f - q_r;

    /*reaction 12: O + C2H4 <=> CH3 + HCO */
    phi_f = sc[2]*sc[16];
    k_f = k_f_save[11];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[13];
    Kc = Kc_save[11];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[11] = q_f - q_r;

    /*reaction 13: O + C2H5 <=> CH3 + CH2O */
    phi_f = sc[2]*sc[17];
    k_f = k_f_save[12];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[14];
    Kc = Kc_save[12];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[12] = q_f - q_r;

    /*reaction 14: O + C2H6 <=> OH + C2H5 */
    phi_f = sc[2]*sc[18];
    k_f = k_f_save[13];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[13];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[13] = q_f - q_r;

    /*reaction 15: O2 + CO <=> O + CO2 */
    phi_f = sc[3]*sc[11];
    k_f = k_f_save[14];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[12];
    Kc = Kc_save[14];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[14] = q_f - q_r;

    /*reaction 16: O2 + CH2O <=> HO2 + HCO */
    phi_f = sc[3]*sc[14];
    k_f = k_f_save[15];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = Kc_save[15];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[15] = q_f - q_r;

    /*reaction 17: H + O2 + M <=> HO2 + M */
    phi_f = sc[1]*sc[3];
    alpha = mixture + -1*sc[3] + -1*sc[5] + -0.25*sc[11] + 0.5*sc[12] + 0.5*sc[18] + -1*sc[19] + -1*sc[20];
    k_f = alpha * k_f_save[16];
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = Kc_save[16];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[16] = q_f - q_r;

    /*reaction 18: H + 2 O2 <=> HO2 + O2 */
    phi_f = sc[1]*sc[3]*sc[3];
    k_f = k_f_save[17];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = Kc_save[17];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[17] = q_f - q_r;

    /*reaction 19: H + O2 + H2O <=> HO2 + H2O */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = k_f_save[18];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = Kc_save[18];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[18] = q_f - q_r;

    /*reaction 20: H + O2 + N2 <=> HO2 + N2 */
    phi_f = sc[1]*sc[3]*sc[19];
    k_f = k_f_save[19];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[19];
    Kc = Kc_save[19];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[19] = q_f - q_r;

    /*reaction 21: H + O2 + AR <=> HO2 + AR */
    phi_f = sc[1]*sc[3]*sc[20];
    k_f = k_f_save[20];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[20];
    Kc = Kc_save[20];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[20] = q_f - q_r;

    /*reaction 22: H + O2 <=> O + OH */
    phi_f = sc[1]*sc[3];
    k_f = k_f_save[21];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = Kc_save[21];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[21] = q_f - q_r;

    /*reaction 23: 2 H + M <=> H2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[0] + -1*sc[5] + sc[10] + -1*sc[12] + 2*sc[18] + -0.37*sc[20];
    k_f = alpha * k_f_save[22];
    q_f = phi_f * k_f;
    phi_r = sc[0];
    Kc = Kc_save[22];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[22] = q_f - q_r;

    /*reaction 24: 2 H + H2 <=> 2 H2 */
    phi_f = sc[1]*sc[1]*sc[0];
    k_f = k_f_save[23];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[0];
    Kc = Kc_save[23];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[23] = q_f - q_r;

    /*reaction 25: 2 H + H2O <=> H2 + H2O */
    phi_f = sc[1]*sc[1]*sc[5];
    k_f = k_f_save[24];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[5];
    Kc = Kc_save[24];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[24] = q_f - q_r;

    /*reaction 26: 2 H + CO2 <=> H2 + CO2 */
    phi_f = sc[1]*sc[1]*sc[12];
    k_f = k_f_save[25];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[12];
    Kc = Kc_save[25];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[25] = q_f - q_r;

    /*reaction 27: H + OH + M <=> H2O + M */
    phi_f = sc[1]*sc[4];
    alpha = mixture + -0.27*sc[0] + 2.65*sc[5] + sc[10] + 2*sc[18] + -0.62*sc[20];
    k_f = alpha * k_f_save[26];
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = Kc_save[26];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[26] = q_f - q_r;

    /*reaction 28: H + HO2 <=> O2 + H2 */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[27];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = Kc_save[27];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[27] = q_f - q_r;

    /*reaction 29: H + HO2 <=> 2 OH */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[28];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = Kc_save[28];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[28] = q_f - q_r;

    /*reaction 30: H + CH2 (+M) <=> CH3 (+M) */
    phi_f = sc[1]*sc[7];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[29];
    redP = 1e-12 * alpha / k_f * 3.2e+27*exp(-3.14*tc[0]-618.95649893980521483*invT);
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
    phi_r = sc[9];
    Kc = Kc_save[29];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[29] = q_f - q_r;

    /*reaction 31: H + CH3 (+M) <=> CH4 (+M) */
    phi_f = sc[1]*sc[9];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[30];
    redP = 1e-12 * alpha / k_f * 2.477e+33*exp(-4.76*tc[0]-1227.8486645635159675*invT);
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
    phi_r = sc[10];
    Kc = Kc_save[30];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[30] = q_f - q_r;

    /*reaction 32: H + CH4 <=> CH3 + H2 */
    phi_f = sc[1]*sc[10];
    k_f = k_f_save[31];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = Kc_save[31];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[31] = q_f - q_r;

    /*reaction 33: H + HCO (+M) <=> CH2O (+M) */
    phi_f = sc[1]*sc[13];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[32];
    redP = 1e-12 * alpha / k_f * 1.35e+24*exp(-2.57*tc[0]-717.08374877172548167*invT);
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
    phi_r = sc[14];
    Kc = Kc_save[32];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[32] = q_f - q_r;

    /*reaction 34: H + HCO <=> H2 + CO */
    phi_f = sc[1]*sc[13];
    k_f = k_f_save[33];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[11];
    Kc = Kc_save[33];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[33] = q_f - q_r;

    /*reaction 35: H + CH2O (+M) <=> CH3O (+M) */
    phi_f = sc[1]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18];
    k_f = k_f_save[34];
    redP = 1e-12 * alpha / k_f * 2.2e+30*exp(-4.8*tc[0]-2797.8846618742413739*invT);
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
    phi_r = sc[15];
    Kc = Kc_save[34];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[34] = q_f - q_r;

    /*reaction 36: H + CH2O <=> HCO + H2 */
    phi_f = sc[1]*sc[14];
    k_f = k_f_save[35];
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[0];
    Kc = Kc_save[35];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[35] = q_f - q_r;

    /*reaction 37: H + CH3O <=> OH + CH3 */
    phi_f = sc[1]*sc[15];
    k_f = k_f_save[36];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = Kc_save[36];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[36] = q_f - q_r;

    /*reaction 38: H + C2H4 (+M) <=> C2H5 (+M) */
    phi_f = sc[1]*sc[16];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[37];
    redP = 1e-12 * alpha / k_f * 1.2e+42*exp(-7.62*tc[0]-3507.4201606588962932*invT);
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
    phi_r = sc[17];
    Kc = Kc_save[37];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[37] = q_f - q_r;

    /*reaction 39: H + C2H5 (+M) <=> C2H6 (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[38];
    redP = 1e-12 * alpha / k_f * 1.99e+41*exp(-7.08*tc[0]-3364.0034109045514015*invT);
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
    phi_r = sc[18];
    Kc = Kc_save[38];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[38] = q_f - q_r;

    /*reaction 40: H + C2H6 <=> C2H5 + H2 */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[39];
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[0];
    Kc = Kc_save[39];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[39] = q_f - q_r;

    /*reaction 41: H2 + CO (+M) <=> CH2O (+M) */
    phi_f = sc[0]*sc[11];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[40];
    redP = 1e-12 * alpha / k_f * 5.07e+27*exp(-3.42*tc[0]-42446.325760628104035*invT);
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
    phi_r = sc[14];
    Kc = Kc_save[40];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[40] = q_f - q_r;

    /*reaction 42: OH + H2 <=> H + H2O */
    phi_f = sc[4]*sc[0];
    k_f = k_f_save[41];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    Kc = Kc_save[41];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[41] = q_f - q_r;

    /*reaction 43: 2 OH <=> O + H2O */
    phi_f = sc[4]*sc[4];
    k_f = k_f_save[42];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = Kc_save[42];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[42] = q_f - q_r;

    /*reaction 44: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = k_f_save[43];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = Kc_save[43];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[43] = q_f - q_r;

    /*reaction 45: OH + CH2 <=> H + CH2O */
    phi_f = sc[4]*sc[7];
    k_f = k_f_save[44];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[44];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[44] = q_f - q_r;

    /*reaction 46: OH + CH2(S) <=> H + CH2O */
    phi_f = sc[4]*sc[8];
    k_f = k_f_save[45];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[45];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[45] = q_f - q_r;

    /*reaction 47: OH + CH3 <=> CH2 + H2O */
    phi_f = sc[4]*sc[9];
    k_f = k_f_save[46];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[5];
    Kc = Kc_save[46];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[46] = q_f - q_r;

    /*reaction 48: OH + CH3 <=> CH2(S) + H2O */
    phi_f = sc[4]*sc[9];
    k_f = k_f_save[47];
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[5];
    Kc = Kc_save[47];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[47] = q_f - q_r;

    /*reaction 49: OH + CH4 <=> CH3 + H2O */
    phi_f = sc[4]*sc[10];
    k_f = k_f_save[48];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[5];
    Kc = Kc_save[48];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[48] = q_f - q_r;

    /*reaction 50: OH + CO <=> H + CO2 */
    phi_f = sc[4]*sc[11];
    k_f = k_f_save[49];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = Kc_save[49];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[49] = q_f - q_r;

    /*reaction 51: OH + HCO <=> H2O + CO */
    phi_f = sc[4]*sc[13];
    k_f = k_f_save[50];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[11];
    Kc = Kc_save[50];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[50] = q_f - q_r;

    /*reaction 52: OH + CH2O <=> HCO + H2O */
    phi_f = sc[4]*sc[14];
    k_f = k_f_save[51];
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[5];
    Kc = Kc_save[51];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[51] = q_f - q_r;

    /*reaction 53: OH + C2H6 <=> C2H5 + H2O */
    phi_f = sc[4]*sc[18];
    k_f = k_f_save[52];
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[5];
    Kc = Kc_save[52];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[52] = q_f - q_r;

    /*reaction 54: HO2 + CH2 <=> OH + CH2O */
    phi_f = sc[6]*sc[7];
    k_f = k_f_save[53];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[14];
    Kc = Kc_save[53];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[53] = q_f - q_r;

    /*reaction 55: HO2 + CH3 <=> O2 + CH4 */
    phi_f = sc[6]*sc[9];
    k_f = k_f_save[54];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[10];
    Kc = Kc_save[54];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[54] = q_f - q_r;

    /*reaction 56: HO2 + CH3 <=> OH + CH3O */
    phi_f = sc[6]*sc[9];
    k_f = k_f_save[55];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[15];
    Kc = Kc_save[55];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[55] = q_f - q_r;

    /*reaction 57: HO2 + CO <=> OH + CO2 */
    phi_f = sc[6]*sc[11];
    k_f = k_f_save[56];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = Kc_save[56];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[56] = q_f - q_r;

    /*reaction 58: CH2 + O2 <=> OH + HCO */
    phi_f = sc[7]*sc[3];
    k_f = k_f_save[57];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[13];
    Kc = Kc_save[57];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[57] = q_f - q_r;

    /*reaction 59: CH2 + H2 <=> H + CH3 */
    phi_f = sc[7]*sc[0];
    k_f = k_f_save[58];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[9];
    Kc = Kc_save[58];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[58] = q_f - q_r;

    /*reaction 60: CH2 + CH3 <=> H + C2H4 */
    phi_f = sc[7]*sc[9];
    k_f = k_f_save[59];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_save[59];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[59] = q_f - q_r;

    /*reaction 61: CH2 + CH4 <=> 2 CH3 */
    phi_f = sc[7]*sc[10];
    k_f = k_f_save[60];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[9];
    Kc = Kc_save[60];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[60] = q_f - q_r;

    /*reaction 62: CH2(S) + N2 <=> CH2 + N2 */
    phi_f = sc[8]*sc[19];
    k_f = k_f_save[61];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[19];
    Kc = Kc_save[61];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[61] = q_f - q_r;

    /*reaction 63: CH2(S) + AR <=> CH2 + AR */
    phi_f = sc[8]*sc[20];
    k_f = k_f_save[62];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[20];
    Kc = Kc_save[62];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[62] = q_f - q_r;

    /*reaction 64: CH2(S) + O2 <=> H + OH + CO */
    phi_f = sc[8]*sc[3];
    k_f = k_f_save[63];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[4]*sc[11];
    Kc = Kc_save[63];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[63] = q_f - q_r;

    /*reaction 65: CH2(S) + O2 <=> CO + H2O */
    phi_f = sc[8]*sc[3];
    k_f = k_f_save[64];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[64];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[64] = q_f - q_r;

    /*reaction 66: CH2(S) + H2 <=> CH3 + H */
    phi_f = sc[8]*sc[0];
    k_f = k_f_save[65];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[1];
    Kc = Kc_save[65];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[65] = q_f - q_r;

    /*reaction 67: CH2(S) + H2O <=> CH2 + H2O */
    phi_f = sc[8]*sc[5];
    k_f = k_f_save[66];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[5];
    Kc = Kc_save[66];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[66] = q_f - q_r;

    /*reaction 68: CH2(S) + CH3 <=> H + C2H4 */
    phi_f = sc[8]*sc[9];
    k_f = k_f_save[67];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_save[67];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[67] = q_f - q_r;

    /*reaction 69: CH2(S) + CH4 <=> 2 CH3 */
    phi_f = sc[8]*sc[10];
    k_f = k_f_save[68];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[9];
    Kc = Kc_save[68];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[68] = q_f - q_r;

    /*reaction 70: CH2(S) + CO <=> CH2 + CO */
    phi_f = sc[8]*sc[11];
    k_f = k_f_save[69];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[11];
    Kc = Kc_save[69];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[69] = q_f - q_r;

    /*reaction 71: CH2(S) + CO2 <=> CH2 + CO2 */
    phi_f = sc[8]*sc[12];
    k_f = k_f_save[70];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[12];
    Kc = Kc_save[70];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[70] = q_f - q_r;

    /*reaction 72: CH2(S) + CO2 <=> CO + CH2O */
    phi_f = sc[8]*sc[12];
    k_f = k_f_save[71];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[14];
    Kc = Kc_save[71];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[71] = q_f - q_r;

    /*reaction 73: CH3 + O2 <=> O + CH3O */
    phi_f = sc[9]*sc[3];
    k_f = k_f_save[72];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[15];
    Kc = Kc_save[72];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[72] = q_f - q_r;

    /*reaction 74: CH3 + O2 <=> OH + CH2O */
    phi_f = sc[9]*sc[3];
    k_f = k_f_save[73];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[14];
    Kc = Kc_save[73];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[73] = q_f - q_r;

    /*reaction 75: 2 CH3 (+M) <=> C2H6 (+M) */
    phi_f = sc[9]*sc[9];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[74];
    redP = 1e-12 * alpha / k_f * 1.77e+50*exp(-9.67*tc[0]-3130.0076613053565779*invT);
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
    phi_r = sc[18];
    Kc = Kc_save[74];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[74] = q_f - q_r;

    /*reaction 76: 2 CH3 <=> H + C2H5 */
    phi_f = sc[9]*sc[9];
    k_f = k_f_save[75];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[75];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[75] = q_f - q_r;

    /*reaction 77: CH3 + HCO <=> CH4 + CO */
    phi_f = sc[9]*sc[13];
    k_f = k_f_save[76];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[11];
    Kc = Kc_save[76];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[76] = q_f - q_r;

    /*reaction 78: CH3 + CH2O <=> HCO + CH4 */
    phi_f = sc[9]*sc[14];
    k_f = k_f_save[77];
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[10];
    Kc = Kc_save[77];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[77] = q_f - q_r;

    /*reaction 79: CH3 + C2H6 <=> C2H5 + CH4 */
    phi_f = sc[9]*sc[18];
    k_f = k_f_save[78];
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[10];
    Kc = Kc_save[78];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[78] = q_f - q_r;

    /*reaction 80: HCO + H2O <=> H + CO + H2O */
    phi_f = sc[13]*sc[5];
    k_f = k_f_save[79];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[11]*sc[5];
    Kc = Kc_save[79];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[79] = q_f - q_r;

    /*reaction 81: HCO + M <=> H + CO + M */
    phi_f = sc[13];
    alpha = mixture + sc[0] + -1*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18];
    k_f = alpha * k_f_save[80];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[11];
    Kc = Kc_save[80];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[80] = q_f - q_r;

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    phi_f = sc[13]*sc[3];
    k_f = k_f_save[81];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[11];
    Kc = Kc_save[81];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[81] = q_f - q_r;

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    phi_f = sc[15]*sc[3];
    k_f = k_f_save[82];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[14];
    Kc = Kc_save[82];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[82] = q_f - q_r;

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    phi_f = sc[17]*sc[3];
    k_f = k_f_save[83];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[16];
    Kc = Kc_save[83];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[83] = q_f - q_r;

    return;
}


/*compute the progress rate for each reaction */
void progressRateFR(double * q_f, double * q_r, double * sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[21];                /*Gibbs free energy */
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
    for (id = 0; id < 21; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    if (T != T_save)
    {
        T_save = T;

        k_f_save[0] = 1e-12 * 5e+17*exp(-1*tc[0]);
        k_f_save[1] = 1e-06 * 50000*exp(2.67*tc[0]-3165.2328279116868543*invT);
        k_f_save[2] = 1e-06 * 2e+13;
        k_f_save[3] = 1e-06 * 8e+13;
        k_f_save[4] = 1e-06 * 1.5e+13;
        k_f_save[5] = 1e-06 * 8.43e+13;
        k_f_save[6] = 1e-06 * 1.02e+09*exp(1.5*tc[0]-4327.6633259205891591*invT);
        k_f_save[7] = 1e-12 * 6.02e+14*exp(-1509.6499974141590883*invT);
        k_f_save[8] = 1e-06 * 3e+13;
        k_f_save[9] = 1e-06 * 3e+13;
        k_f_save[10] = 1e-06 * 3.9e+13*exp(-1781.3869969487077469*invT);
        k_f_save[11] = 1e-06 * 1.92e+07*exp(1.83*tc[0]-110.7076664770383303*invT);
        k_f_save[12] = 1e-06 * 1.32e+14;
        k_f_save[13] = 1e-06 * 8.98e+07*exp(1.92*tc[0]-2863.3028284288548093*invT);
        k_f_save[14] = 1e-06 * 2.5e+12*exp(-24053.756625465601246*invT);
        k_f_save[15] = 1e-06 * 1e+14*exp(-20128.666632188789663*invT);
        k_f_save[16] = 1e-12 * 2.8e+18*exp(-0.86*tc[0]);
        k_f_save[17] = 1e-12 * 3e+20*exp(-1.72*tc[0]);
        k_f_save[18] = 1e-12 * 9.38e+18*exp(-0.76*tc[0]);
        k_f_save[19] = 1e-12 * 3.75e+20*exp(-1.72*tc[0]);
        k_f_save[20] = 1e-12 * 7e+17*exp(-0.8*tc[0]);
        k_f_save[21] = 1e-06 * 8.3e+13*exp(-7252.8618042434254676*invT);
        k_f_save[22] = 1e-12 * 1e+18*exp(-1*tc[0]);
        k_f_save[23] = 1e-12 * 9e+16*exp(-0.6*tc[0]);
        k_f_save[24] = 1e-12 * 6e+19*exp(-1.25*tc[0]);
        k_f_save[25] = 1e-12 * 5.5e+20*exp(-2*tc[0]);
        k_f_save[26] = 1e-12 * 2.2e+22*exp(-2*tc[0]);
        k_f_save[27] = 1e-06 * 2.8e+13*exp(-537.43539907944057177*invT);
        k_f_save[28] = 1e-06 * 1.34e+14*exp(-319.54258278599701271*invT);
        k_f_save[29] = 1e-06 * 2.5e+16*exp(-0.8*tc[0]);
        k_f_save[30] = 1e-06 * 1.27e+16*exp(-0.63*tc[0]-192.73198300320765952*invT);
        k_f_save[31] = 1e-06 * 6.6e+08*exp(1.62*tc[0]-5454.8686573231616421*invT);
        k_f_save[32] = 1e-06 * 1.09e+12*exp(0.48*tc[0]+130.83633310922709825*invT);
        k_f_save[33] = 1e-06 * 7.34e+13;
        k_f_save[34] = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1308.3633310922712099*invT);
        k_f_save[35] = 1e-06 * 2.3e+10*exp(1.05*tc[0]-1648.0345805104568626*invT);
        k_f_save[36] = 1e-06 * 3.2e+13;
        k_f_save[37] = 1e-06 * 1.08e+12*exp(0.454*tc[0]-915.85433176458980142*invT);
        k_f_save[38] = 1e-06 * 5.21e+17*exp(-0.99*tc[0]-795.08233197145705162*invT);
        k_f_save[39] = 1e-06 * 1.15e+08*exp(1.9*tc[0]-3789.2214935095394139*invT);
        k_f_save[40] = 1e-06 * 4.3e+07*exp(1.5*tc[0]-40056.046598055690993*invT);
        k_f_save[41] = 1e-06 * 2.16e+08*exp(1.51*tc[0]-1726.0331637101885462*invT);
        k_f_save[42] = 1e-06 * 35700*exp(2.4*tc[0]+1061.7871648479585929*invT);
        k_f_save[43] = 1e-06 * 2.9e+13*exp(+251.60833290235981963*invT);
        k_f_save[44] = 1e-06 * 2e+13;
        k_f_save[45] = 1e-06 * 3e+13;
        k_f_save[46] = 1e-06 * 5.6e+07*exp(1.6*tc[0]-2727.4343286615808211*invT);
        k_f_save[47] = 1e-06 * 2.501e+13;
        k_f_save[48] = 1e-06 * 1e+08*exp(1.6*tc[0]-1570.0359973107254064*invT);
        k_f_save[49] = 1e-06 * 4.76e+07*exp(1.228*tc[0]-35.225166606330375885*invT);
        k_f_save[50] = 1e-06 * 5e+13;
        k_f_save[51] = 1e-06 * 3.43e+09*exp(1.18*tc[0]+224.93784961470970529*invT);
        k_f_save[52] = 1e-06 * 3.54e+06*exp(2.12*tc[0]-437.79849925010609013*invT);
        k_f_save[53] = 1e-06 * 2e+13;
        k_f_save[54] = 1e-06 * 1e+12;
        k_f_save[55] = 1e-06 * 2e+13;
        k_f_save[56] = 1e-06 * 1.5e+14*exp(-11875.913312991384373*invT);
        k_f_save[57] = 1e-06 * 1.32e+13*exp(-754.82499870707954415*invT);
        k_f_save[58] = 1e-06 * 500000*exp(2*tc[0]-3638.256493768123164*invT);
        k_f_save[59] = 1e-06 * 4e+13;
        k_f_save[60] = 1e-06 * 2.46e+06*exp(2*tc[0]-4161.6018262050320118*invT);
        k_f_save[61] = 1e-06 * 1.5e+13*exp(-301.92999948283181766*invT);
        k_f_save[62] = 1e-06 * 9e+12*exp(-301.92999948283181766*invT);
        k_f_save[63] = 1e-06 * 2.8e+13;
        k_f_save[64] = 1e-06 * 1.2e+13;
        k_f_save[65] = 1e-06 * 7e+13;
        k_f_save[66] = 1e-06 * 3e+13;
        k_f_save[67] = 1e-06 * 1.2e+13*exp(+286.83349950869023814*invT);
        k_f_save[68] = 1e-06 * 1.6e+13*exp(+286.83349950869023814*invT);
        k_f_save[69] = 1e-06 * 9e+12;
        k_f_save[70] = 1e-06 * 7e+12;
        k_f_save[71] = 1e-06 * 1.4e+13;
        k_f_save[72] = 1e-06 * 2.675e+13*exp(-14492.639975175927248*invT);
        k_f_save[73] = 1e-06 * 3.6e+10*exp(-4498.7569922941938785*invT);
        k_f_save[74] = 1e-06 * 2.12e+16*exp(-0.97*tc[0]-311.99433279892622295*invT);
        k_f_save[75] = 1e-06 * 4.99e+12*exp(0.1*tc[0]-5334.096657530029006*invT);
        k_f_save[76] = 1e-06 * 2.648e+13;
        k_f_save[77] = 1e-06 * 3320*exp(2.81*tc[0]-2948.8496616156576238*invT);
        k_f_save[78] = 1e-06 * 6.14e+06*exp(1.74*tc[0]-5258.6141576593208811*invT);
        k_f_save[79] = 1e-06 * 2.244e+18*exp(-1*tc[0]-8554.6833186802341515*invT);
        k_f_save[80] = 1e-06 * 1.87e+17*exp(-1*tc[0]-8554.6833186802341515*invT);
        k_f_save[81] = 1e-06 * 7.6e+12*exp(-201.28666632188787844*invT);
        k_f_save[82] = 1e-06 * 4.28e-13*exp(7.6*tc[0]+1776.3548302906606295*invT);
        k_f_save[83] = 1e-06 * 8.4e+11*exp(-1949.9645799932889076*invT);

        Kc_save[0] = 1.0 / (refC) * exp((g_RT[2] + g_RT[1]) - (g_RT[4]));
        Kc_save[1] = exp((g_RT[2] + g_RT[0]) - (g_RT[1] + g_RT[4]));
        Kc_save[2] = exp((g_RT[2] + g_RT[6]) - (g_RT[4] + g_RT[3]));
        Kc_save[3] = exp((g_RT[2] + g_RT[7]) - (g_RT[1] + g_RT[13]));
        Kc_save[4] = exp((g_RT[2] + g_RT[8]) - (g_RT[1] + g_RT[13]));
        Kc_save[5] = exp((g_RT[2] + g_RT[9]) - (g_RT[1] + g_RT[14]));
        Kc_save[6] = exp((g_RT[2] + g_RT[10]) - (g_RT[4] + g_RT[9]));
        Kc_save[7] = 1.0 / (refC) * exp((g_RT[2] + g_RT[11]) - (g_RT[12]));
        Kc_save[8] = exp((g_RT[2] + g_RT[13]) - (g_RT[4] + g_RT[11]));
        Kc_save[9] = exp((g_RT[2] + g_RT[13]) - (g_RT[1] + g_RT[12]));
        Kc_save[10] = exp((g_RT[2] + g_RT[14]) - (g_RT[4] + g_RT[13]));
        Kc_save[11] = exp((g_RT[2] + g_RT[16]) - (g_RT[9] + g_RT[13]));
        Kc_save[12] = exp((g_RT[2] + g_RT[17]) - (g_RT[9] + g_RT[14]));
        Kc_save[13] = exp((g_RT[2] + g_RT[18]) - (g_RT[4] + g_RT[17]));
        Kc_save[14] = exp((g_RT[3] + g_RT[11]) - (g_RT[2] + g_RT[12]));
        Kc_save[15] = exp((g_RT[3] + g_RT[14]) - (g_RT[6] + g_RT[13]));
        Kc_save[16] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));
        Kc_save[17] = 1.0 / (refC) * exp((g_RT[1] + 2 * g_RT[3]) - (g_RT[6] + g_RT[3]));
        Kc_save[18] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));
        Kc_save[19] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[19]) - (g_RT[6] + g_RT[19]));
        Kc_save[20] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[20]) - (g_RT[6] + g_RT[20]));
        Kc_save[21] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));
        Kc_save[22] = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));
        Kc_save[23] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));
        Kc_save[24] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));
        Kc_save[25] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[12]) - (g_RT[0] + g_RT[12]));
        Kc_save[26] = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));
        Kc_save[27] = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));
        Kc_save[28] = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));
        Kc_save[29] = 1.0 / (refC) * exp((g_RT[1] + g_RT[7]) - (g_RT[9]));
        Kc_save[30] = 1.0 / (refC) * exp((g_RT[1] + g_RT[9]) - (g_RT[10]));
        Kc_save[31] = exp((g_RT[1] + g_RT[10]) - (g_RT[9] + g_RT[0]));
        Kc_save[32] = 1.0 / (refC) * exp((g_RT[1] + g_RT[13]) - (g_RT[14]));
        Kc_save[33] = exp((g_RT[1] + g_RT[13]) - (g_RT[0] + g_RT[11]));
        Kc_save[34] = 1.0 / (refC) * exp((g_RT[1] + g_RT[14]) - (g_RT[15]));
        Kc_save[35] = exp((g_RT[1] + g_RT[14]) - (g_RT[13] + g_RT[0]));
        Kc_save[36] = exp((g_RT[1] + g_RT[15]) - (g_RT[4] + g_RT[9]));
        Kc_save[37] = 1.0 / (refC) * exp((g_RT[1] + g_RT[16]) - (g_RT[17]));
        Kc_save[38] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[18]));
        Kc_save[39] = exp((g_RT[1] + g_RT[18]) - (g_RT[17] + g_RT[0]));
        Kc_save[40] = 1.0 / (refC) * exp((g_RT[0] + g_RT[11]) - (g_RT[14]));
        Kc_save[41] = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));
        Kc_save[42] = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));
        Kc_save[43] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
        Kc_save[44] = exp((g_RT[4] + g_RT[7]) - (g_RT[1] + g_RT[14]));
        Kc_save[45] = exp((g_RT[4] + g_RT[8]) - (g_RT[1] + g_RT[14]));
        Kc_save[46] = exp((g_RT[4] + g_RT[9]) - (g_RT[7] + g_RT[5]));
        Kc_save[47] = exp((g_RT[4] + g_RT[9]) - (g_RT[8] + g_RT[5]));
        Kc_save[48] = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[5]));
        Kc_save[49] = exp((g_RT[4] + g_RT[11]) - (g_RT[1] + g_RT[12]));
        Kc_save[50] = exp((g_RT[4] + g_RT[13]) - (g_RT[5] + g_RT[11]));
        Kc_save[51] = exp((g_RT[4] + g_RT[14]) - (g_RT[13] + g_RT[5]));
        Kc_save[52] = exp((g_RT[4] + g_RT[18]) - (g_RT[17] + g_RT[5]));
        Kc_save[53] = exp((g_RT[6] + g_RT[7]) - (g_RT[4] + g_RT[14]));
        Kc_save[54] = exp((g_RT[6] + g_RT[9]) - (g_RT[3] + g_RT[10]));
        Kc_save[55] = exp((g_RT[6] + g_RT[9]) - (g_RT[4] + g_RT[15]));
        Kc_save[56] = exp((g_RT[6] + g_RT[11]) - (g_RT[4] + g_RT[12]));
        Kc_save[57] = exp((g_RT[7] + g_RT[3]) - (g_RT[4] + g_RT[13]));
        Kc_save[58] = exp((g_RT[7] + g_RT[0]) - (g_RT[1] + g_RT[9]));
        Kc_save[59] = exp((g_RT[7] + g_RT[9]) - (g_RT[1] + g_RT[16]));
        Kc_save[60] = exp((g_RT[7] + g_RT[10]) - (2 * g_RT[9]));
        Kc_save[61] = exp((g_RT[8] + g_RT[19]) - (g_RT[7] + g_RT[19]));
        Kc_save[62] = exp((g_RT[8] + g_RT[20]) - (g_RT[7] + g_RT[20]));
        Kc_save[63] = refC * exp((g_RT[8] + g_RT[3]) - (g_RT[1] + g_RT[4] + g_RT[11]));
        Kc_save[64] = exp((g_RT[8] + g_RT[3]) - (g_RT[11] + g_RT[5]));
        Kc_save[65] = exp((g_RT[8] + g_RT[0]) - (g_RT[9] + g_RT[1]));
        Kc_save[66] = exp((g_RT[8] + g_RT[5]) - (g_RT[7] + g_RT[5]));
        Kc_save[67] = exp((g_RT[8] + g_RT[9]) - (g_RT[1] + g_RT[16]));
        Kc_save[68] = exp((g_RT[8] + g_RT[10]) - (2 * g_RT[9]));
        Kc_save[69] = exp((g_RT[8] + g_RT[11]) - (g_RT[7] + g_RT[11]));
        Kc_save[70] = exp((g_RT[8] + g_RT[12]) - (g_RT[7] + g_RT[12]));
        Kc_save[71] = exp((g_RT[8] + g_RT[12]) - (g_RT[11] + g_RT[14]));
        Kc_save[72] = exp((g_RT[9] + g_RT[3]) - (g_RT[2] + g_RT[15]));
        Kc_save[73] = exp((g_RT[9] + g_RT[3]) - (g_RT[4] + g_RT[14]));
        Kc_save[74] = 1.0 / (refC) * exp((2 * g_RT[9]) - (g_RT[18]));
        Kc_save[75] = exp((2 * g_RT[9]) - (g_RT[1] + g_RT[17]));
        Kc_save[76] = exp((g_RT[9] + g_RT[13]) - (g_RT[10] + g_RT[11]));
        Kc_save[77] = exp((g_RT[9] + g_RT[14]) - (g_RT[13] + g_RT[10]));
        Kc_save[78] = exp((g_RT[9] + g_RT[18]) - (g_RT[17] + g_RT[10]));
        Kc_save[79] = refC * exp((g_RT[13] + g_RT[5]) - (g_RT[1] + g_RT[11] + g_RT[5]));
        Kc_save[80] = refC * exp((g_RT[13]) - (g_RT[1] + g_RT[11]));
        Kc_save[81] = exp((g_RT[13] + g_RT[3]) - (g_RT[6] + g_RT[11]));
        Kc_save[82] = exp((g_RT[15] + g_RT[3]) - (g_RT[6] + g_RT[14]));
        Kc_save[83] = exp((g_RT[17] + g_RT[3]) - (g_RT[6] + g_RT[16]));
    }

    /*reaction 1: O + H + M <=> OH + M */
    phi_f = sc[2]*sc[1];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = alpha * k_f_save[0];
    q_f[0] = phi_f * k_f;
    phi_r = sc[4];
    Kc = Kc_save[0];
    k_r = k_f / Kc;
    q_r[0] = phi_r * k_r;

    /*reaction 2: O + H2 <=> H + OH */
    phi_f = sc[2]*sc[0];
    k_f = k_f_save[1];
    q_f[1] = phi_f * k_f;
    phi_r = sc[1]*sc[4];
    Kc = Kc_save[1];
    k_r = k_f / Kc;
    q_r[1] = phi_r * k_r;

    /*reaction 3: O + HO2 <=> OH + O2 */
    phi_f = sc[2]*sc[6];
    k_f = k_f_save[2];
    q_f[2] = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = Kc_save[2];
    k_r = k_f / Kc;
    q_r[2] = phi_r * k_r;

    /*reaction 4: O + CH2 <=> H + HCO */
    phi_f = sc[2]*sc[7];
    k_f = k_f_save[3];
    q_f[3] = phi_f * k_f;
    phi_r = sc[1]*sc[13];
    Kc = Kc_save[3];
    k_r = k_f / Kc;
    q_r[3] = phi_r * k_r;

    /*reaction 5: O + CH2(S) <=> H + HCO */
    phi_f = sc[2]*sc[8];
    k_f = k_f_save[4];
    q_f[4] = phi_f * k_f;
    phi_r = sc[1]*sc[13];
    Kc = Kc_save[4];
    k_r = k_f / Kc;
    q_r[4] = phi_r * k_r;

    /*reaction 6: O + CH3 <=> H + CH2O */
    phi_f = sc[2]*sc[9];
    k_f = k_f_save[5];
    q_f[5] = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[5];
    k_r = k_f / Kc;
    q_r[5] = phi_r * k_r;

    /*reaction 7: O + CH4 <=> OH + CH3 */
    phi_f = sc[2]*sc[10];
    k_f = k_f_save[6];
    q_f[6] = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = Kc_save[6];
    k_r = k_f / Kc;
    q_r[6] = phi_r * k_r;

    /*reaction 8: O + CO + M <=> CO2 + M */
    phi_f = sc[2]*sc[11];
    alpha = mixture + sc[0] + 5*sc[3] + 5*sc[5] + sc[10] + 0.5*sc[11] + 2.5*sc[12] + 2*sc[18] + -0.5*sc[20];
    k_f = alpha * k_f_save[7];
    q_f[7] = phi_f * k_f;
    phi_r = sc[12];
    Kc = Kc_save[7];
    k_r = k_f / Kc;
    q_r[7] = phi_r * k_r;

    /*reaction 9: O + HCO <=> OH + CO */
    phi_f = sc[2]*sc[13];
    k_f = k_f_save[8];
    q_f[8] = phi_f * k_f;
    phi_r = sc[4]*sc[11];
    Kc = Kc_save[8];
    k_r = k_f / Kc;
    q_r[8] = phi_r * k_r;

    /*reaction 10: O + HCO <=> H + CO2 */
    phi_f = sc[2]*sc[13];
    k_f = k_f_save[9];
    q_f[9] = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = Kc_save[9];
    k_r = k_f / Kc;
    q_r[9] = phi_r * k_r;

    /*reaction 11: O + CH2O <=> OH + HCO */
    phi_f = sc[2]*sc[14];
    k_f = k_f_save[10];
    q_f[10] = phi_f * k_f;
    phi_r = sc[4]*sc[13];
    Kc = Kc_save[10];
    k_r = k_f / Kc;
    q_r[10] = phi_r * k_r;

    /*reaction 12: O + C2H4 <=> CH3 + HCO */
    phi_f = sc[2]*sc[16];
    k_f = k_f_save[11];
    q_f[11] = phi_f * k_f;
    phi_r = sc[9]*sc[13];
    Kc = Kc_save[11];
    k_r = k_f / Kc;
    q_r[11] = phi_r * k_r;

    /*reaction 13: O + C2H5 <=> CH3 + CH2O */
    phi_f = sc[2]*sc[17];
    k_f = k_f_save[12];
    q_f[12] = phi_f * k_f;
    phi_r = sc[9]*sc[14];
    Kc = Kc_save[12];
    k_r = k_f / Kc;
    q_r[12] = phi_r * k_r;

    /*reaction 14: O + C2H6 <=> OH + C2H5 */
    phi_f = sc[2]*sc[18];
    k_f = k_f_save[13];
    q_f[13] = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[13];
    k_r = k_f / Kc;
    q_r[13] = phi_r * k_r;

    /*reaction 15: O2 + CO <=> O + CO2 */
    phi_f = sc[3]*sc[11];
    k_f = k_f_save[14];
    q_f[14] = phi_f * k_f;
    phi_r = sc[2]*sc[12];
    Kc = Kc_save[14];
    k_r = k_f / Kc;
    q_r[14] = phi_r * k_r;

    /*reaction 16: O2 + CH2O <=> HO2 + HCO */
    phi_f = sc[3]*sc[14];
    k_f = k_f_save[15];
    q_f[15] = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = Kc_save[15];
    k_r = k_f / Kc;
    q_r[15] = phi_r * k_r;

    /*reaction 17: H + O2 + M <=> HO2 + M */
    phi_f = sc[1]*sc[3];
    alpha = mixture + -1*sc[3] + -1*sc[5] + -0.25*sc[11] + 0.5*sc[12] + 0.5*sc[18] + -1*sc[19] + -1*sc[20];
    k_f = alpha * k_f_save[16];
    q_f[16] = phi_f * k_f;
    phi_r = sc[6];
    Kc = Kc_save[16];
    k_r = k_f / Kc;
    q_r[16] = phi_r * k_r;

    /*reaction 18: H + 2 O2 <=> HO2 + O2 */
    phi_f = sc[1]*sc[3]*sc[3];
    k_f = k_f_save[17];
    q_f[17] = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = Kc_save[17];
    k_r = k_f / Kc;
    q_r[17] = phi_r * k_r;

    /*reaction 19: H + O2 + H2O <=> HO2 + H2O */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = k_f_save[18];
    q_f[18] = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = Kc_save[18];
    k_r = k_f / Kc;
    q_r[18] = phi_r * k_r;

    /*reaction 20: H + O2 + N2 <=> HO2 + N2 */
    phi_f = sc[1]*sc[3]*sc[19];
    k_f = k_f_save[19];
    q_f[19] = phi_f * k_f;
    phi_r = sc[6]*sc[19];
    Kc = Kc_save[19];
    k_r = k_f / Kc;
    q_r[19] = phi_r * k_r;

    /*reaction 21: H + O2 + AR <=> HO2 + AR */
    phi_f = sc[1]*sc[3]*sc[20];
    k_f = k_f_save[20];
    q_f[20] = phi_f * k_f;
    phi_r = sc[6]*sc[20];
    Kc = Kc_save[20];
    k_r = k_f / Kc;
    q_r[20] = phi_r * k_r;

    /*reaction 22: H + O2 <=> O + OH */
    phi_f = sc[1]*sc[3];
    k_f = k_f_save[21];
    q_f[21] = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = Kc_save[21];
    k_r = k_f / Kc;
    q_r[21] = phi_r * k_r;

    /*reaction 23: 2 H + M <=> H2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[0] + -1*sc[5] + sc[10] + -1*sc[12] + 2*sc[18] + -0.37*sc[20];
    k_f = alpha * k_f_save[22];
    q_f[22] = phi_f * k_f;
    phi_r = sc[0];
    Kc = Kc_save[22];
    k_r = k_f / Kc;
    q_r[22] = phi_r * k_r;

    /*reaction 24: 2 H + H2 <=> 2 H2 */
    phi_f = sc[1]*sc[1]*sc[0];
    k_f = k_f_save[23];
    q_f[23] = phi_f * k_f;
    phi_r = sc[0]*sc[0];
    Kc = Kc_save[23];
    k_r = k_f / Kc;
    q_r[23] = phi_r * k_r;

    /*reaction 25: 2 H + H2O <=> H2 + H2O */
    phi_f = sc[1]*sc[1]*sc[5];
    k_f = k_f_save[24];
    q_f[24] = phi_f * k_f;
    phi_r = sc[0]*sc[5];
    Kc = Kc_save[24];
    k_r = k_f / Kc;
    q_r[24] = phi_r * k_r;

    /*reaction 26: 2 H + CO2 <=> H2 + CO2 */
    phi_f = sc[1]*sc[1]*sc[12];
    k_f = k_f_save[25];
    q_f[25] = phi_f * k_f;
    phi_r = sc[0]*sc[12];
    Kc = Kc_save[25];
    k_r = k_f / Kc;
    q_r[25] = phi_r * k_r;

    /*reaction 27: H + OH + M <=> H2O + M */
    phi_f = sc[1]*sc[4];
    alpha = mixture + -0.27*sc[0] + 2.65*sc[5] + sc[10] + 2*sc[18] + -0.62*sc[20];
    k_f = alpha * k_f_save[26];
    q_f[26] = phi_f * k_f;
    phi_r = sc[5];
    Kc = Kc_save[26];
    k_r = k_f / Kc;
    q_r[26] = phi_r * k_r;

    /*reaction 28: H + HO2 <=> O2 + H2 */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[27];
    q_f[27] = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = Kc_save[27];
    k_r = k_f / Kc;
    q_r[27] = phi_r * k_r;

    /*reaction 29: H + HO2 <=> 2 OH */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[28];
    q_f[28] = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = Kc_save[28];
    k_r = k_f / Kc;
    q_r[28] = phi_r * k_r;

    /*reaction 30: H + CH2 (+M) <=> CH3 (+M) */
    phi_f = sc[1]*sc[7];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[29];
    redP = 1e-12 * alpha / k_f * 3.2e+27*exp(-3.14*tc[0]-618.95649893980521483*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.32*exp(T/-78))+ (0.68*exp(T/-1995))+ (exp(-5590/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[29] = phi_f * k_f;
    phi_r = sc[9];
    Kc = Kc_save[29];
    k_r = k_f / Kc;
    q_r[29] = phi_r * k_r;

    /*reaction 31: H + CH3 (+M) <=> CH4 (+M) */
    phi_f = sc[1]*sc[9];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[30];
    redP = 1e-12 * alpha / k_f * 2.477e+33*exp(-4.76*tc[0]-1227.8486645635159675*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.217*exp(T/-74))+ (0.783*exp(T/-2941))+ (exp(-6964/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[30] = phi_f * k_f;
    phi_r = sc[10];
    Kc = Kc_save[30];
    k_r = k_f / Kc;
    q_r[30] = phi_r * k_r;

    /*reaction 32: H + CH4 <=> CH3 + H2 */
    phi_f = sc[1]*sc[10];
    k_f = k_f_save[31];
    q_f[31] = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = Kc_save[31];
    k_r = k_f / Kc;
    q_r[31] = phi_r * k_r;

    /*reaction 33: H + HCO (+M) <=> CH2O (+M) */
    phi_f = sc[1]*sc[13];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[32];
    redP = 1e-12 * alpha / k_f * 1.35e+24*exp(-2.57*tc[0]-717.08374877172548167*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2176*exp(T/-271))+ (0.7824*exp(T/-2755))+ (exp(-6570/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[32] = phi_f * k_f;
    phi_r = sc[14];
    Kc = Kc_save[32];
    k_r = k_f / Kc;
    q_r[32] = phi_r * k_r;

    /*reaction 34: H + HCO <=> H2 + CO */
    phi_f = sc[1]*sc[13];
    k_f = k_f_save[33];
    q_f[33] = phi_f * k_f;
    phi_r = sc[0]*sc[11];
    Kc = Kc_save[33];
    k_r = k_f / Kc;
    q_r[33] = phi_r * k_r;

    /*reaction 35: H + CH2O (+M) <=> CH3O (+M) */
    phi_f = sc[1]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18];
    k_f = k_f_save[34];
    redP = 1e-12 * alpha / k_f * 2.2e+30*exp(-4.8*tc[0]-2797.8846618742413739*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.242*exp(T/-94))+ (0.758*exp(T/-1555))+ (exp(-4200/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[34] = phi_f * k_f;
    phi_r = sc[15];
    Kc = Kc_save[34];
    k_r = k_f / Kc;
    q_r[34] = phi_r * k_r;

    /*reaction 36: H + CH2O <=> HCO + H2 */
    phi_f = sc[1]*sc[14];
    k_f = k_f_save[35];
    q_f[35] = phi_f * k_f;
    phi_r = sc[13]*sc[0];
    Kc = Kc_save[35];
    k_r = k_f / Kc;
    q_r[35] = phi_r * k_r;

    /*reaction 37: H + CH3O <=> OH + CH3 */
    phi_f = sc[1]*sc[15];
    k_f = k_f_save[36];
    q_f[36] = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = Kc_save[36];
    k_r = k_f / Kc;
    q_r[36] = phi_r * k_r;

    /*reaction 38: H + C2H4 (+M) <=> C2H5 (+M) */
    phi_f = sc[1]*sc[16];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[37];
    redP = 1e-12 * alpha / k_f * 1.2e+42*exp(-7.62*tc[0]-3507.4201606588962932*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.0247*exp(T/-210))+ (0.9753*exp(T/-984))+ (exp(-4374/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[37] = phi_f * k_f;
    phi_r = sc[17];
    Kc = Kc_save[37];
    k_r = k_f / Kc;
    q_r[37] = phi_r * k_r;

    /*reaction 39: H + C2H5 (+M) <=> C2H6 (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[38];
    redP = 1e-12 * alpha / k_f * 1.99e+41*exp(-7.08*tc[0]-3364.0034109045514015*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.1578*exp(T/-125))+ (0.8422*exp(T/-2219))+ (exp(-6882/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[38] = phi_f * k_f;
    phi_r = sc[18];
    Kc = Kc_save[38];
    k_r = k_f / Kc;
    q_r[38] = phi_r * k_r;

    /*reaction 40: H + C2H6 <=> C2H5 + H2 */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[39];
    q_f[39] = phi_f * k_f;
    phi_r = sc[17]*sc[0];
    Kc = Kc_save[39];
    k_r = k_f / Kc;
    q_r[39] = phi_r * k_r;

    /*reaction 41: H2 + CO (+M) <=> CH2O (+M) */
    phi_f = sc[0]*sc[11];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[40];
    redP = 1e-12 * alpha / k_f * 5.07e+27*exp(-3.42*tc[0]-42446.325760628104035*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.068*exp(T/-197))+ (0.932*exp(T/-1540))+ (exp(-10300/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[40] = phi_f * k_f;
    phi_r = sc[14];
    Kc = Kc_save[40];
    k_r = k_f / Kc;
    q_r[40] = phi_r * k_r;

    /*reaction 42: OH + H2 <=> H + H2O */
    phi_f = sc[4]*sc[0];
    k_f = k_f_save[41];
    q_f[41] = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    Kc = Kc_save[41];
    k_r = k_f / Kc;
    q_r[41] = phi_r * k_r;

    /*reaction 43: 2 OH <=> O + H2O */
    phi_f = sc[4]*sc[4];
    k_f = k_f_save[42];
    q_f[42] = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = Kc_save[42];
    k_r = k_f / Kc;
    q_r[42] = phi_r * k_r;

    /*reaction 44: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = k_f_save[43];
    q_f[43] = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = Kc_save[43];
    k_r = k_f / Kc;
    q_r[43] = phi_r * k_r;

    /*reaction 45: OH + CH2 <=> H + CH2O */
    phi_f = sc[4]*sc[7];
    k_f = k_f_save[44];
    q_f[44] = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[44];
    k_r = k_f / Kc;
    q_r[44] = phi_r * k_r;

    /*reaction 46: OH + CH2(S) <=> H + CH2O */
    phi_f = sc[4]*sc[8];
    k_f = k_f_save[45];
    q_f[45] = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[45];
    k_r = k_f / Kc;
    q_r[45] = phi_r * k_r;

    /*reaction 47: OH + CH3 <=> CH2 + H2O */
    phi_f = sc[4]*sc[9];
    k_f = k_f_save[46];
    q_f[46] = phi_f * k_f;
    phi_r = sc[7]*sc[5];
    Kc = Kc_save[46];
    k_r = k_f / Kc;
    q_r[46] = phi_r * k_r;

    /*reaction 48: OH + CH3 <=> CH2(S) + H2O */
    phi_f = sc[4]*sc[9];
    k_f = k_f_save[47];
    q_f[47] = phi_f * k_f;
    phi_r = sc[8]*sc[5];
    Kc = Kc_save[47];
    k_r = k_f / Kc;
    q_r[47] = phi_r * k_r;

    /*reaction 49: OH + CH4 <=> CH3 + H2O */
    phi_f = sc[4]*sc[10];
    k_f = k_f_save[48];
    q_f[48] = phi_f * k_f;
    phi_r = sc[9]*sc[5];
    Kc = Kc_save[48];
    k_r = k_f / Kc;
    q_r[48] = phi_r * k_r;

    /*reaction 50: OH + CO <=> H + CO2 */
    phi_f = sc[4]*sc[11];
    k_f = k_f_save[49];
    q_f[49] = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = Kc_save[49];
    k_r = k_f / Kc;
    q_r[49] = phi_r * k_r;

    /*reaction 51: OH + HCO <=> H2O + CO */
    phi_f = sc[4]*sc[13];
    k_f = k_f_save[50];
    q_f[50] = phi_f * k_f;
    phi_r = sc[5]*sc[11];
    Kc = Kc_save[50];
    k_r = k_f / Kc;
    q_r[50] = phi_r * k_r;

    /*reaction 52: OH + CH2O <=> HCO + H2O */
    phi_f = sc[4]*sc[14];
    k_f = k_f_save[51];
    q_f[51] = phi_f * k_f;
    phi_r = sc[13]*sc[5];
    Kc = Kc_save[51];
    k_r = k_f / Kc;
    q_r[51] = phi_r * k_r;

    /*reaction 53: OH + C2H6 <=> C2H5 + H2O */
    phi_f = sc[4]*sc[18];
    k_f = k_f_save[52];
    q_f[52] = phi_f * k_f;
    phi_r = sc[17]*sc[5];
    Kc = Kc_save[52];
    k_r = k_f / Kc;
    q_r[52] = phi_r * k_r;

    /*reaction 54: HO2 + CH2 <=> OH + CH2O */
    phi_f = sc[6]*sc[7];
    k_f = k_f_save[53];
    q_f[53] = phi_f * k_f;
    phi_r = sc[4]*sc[14];
    Kc = Kc_save[53];
    k_r = k_f / Kc;
    q_r[53] = phi_r * k_r;

    /*reaction 55: HO2 + CH3 <=> O2 + CH4 */
    phi_f = sc[6]*sc[9];
    k_f = k_f_save[54];
    q_f[54] = phi_f * k_f;
    phi_r = sc[3]*sc[10];
    Kc = Kc_save[54];
    k_r = k_f / Kc;
    q_r[54] = phi_r * k_r;

    /*reaction 56: HO2 + CH3 <=> OH + CH3O */
    phi_f = sc[6]*sc[9];
    k_f = k_f_save[55];
    q_f[55] = phi_f * k_f;
    phi_r = sc[4]*sc[15];
    Kc = Kc_save[55];
    k_r = k_f / Kc;
    q_r[55] = phi_r * k_r;

    /*reaction 57: HO2 + CO <=> OH + CO2 */
    phi_f = sc[6]*sc[11];
    k_f = k_f_save[56];
    q_f[56] = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = Kc_save[56];
    k_r = k_f / Kc;
    q_r[56] = phi_r * k_r;

    /*reaction 58: CH2 + O2 <=> OH + HCO */
    phi_f = sc[7]*sc[3];
    k_f = k_f_save[57];
    q_f[57] = phi_f * k_f;
    phi_r = sc[4]*sc[13];
    Kc = Kc_save[57];
    k_r = k_f / Kc;
    q_r[57] = phi_r * k_r;

    /*reaction 59: CH2 + H2 <=> H + CH3 */
    phi_f = sc[7]*sc[0];
    k_f = k_f_save[58];
    q_f[58] = phi_f * k_f;
    phi_r = sc[1]*sc[9];
    Kc = Kc_save[58];
    k_r = k_f / Kc;
    q_r[58] = phi_r * k_r;

    /*reaction 60: CH2 + CH3 <=> H + C2H4 */
    phi_f = sc[7]*sc[9];
    k_f = k_f_save[59];
    q_f[59] = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_save[59];
    k_r = k_f / Kc;
    q_r[59] = phi_r * k_r;

    /*reaction 61: CH2 + CH4 <=> 2 CH3 */
    phi_f = sc[7]*sc[10];
    k_f = k_f_save[60];
    q_f[60] = phi_f * k_f;
    phi_r = sc[9]*sc[9];
    Kc = Kc_save[60];
    k_r = k_f / Kc;
    q_r[60] = phi_r * k_r;

    /*reaction 62: CH2(S) + N2 <=> CH2 + N2 */
    phi_f = sc[8]*sc[19];
    k_f = k_f_save[61];
    q_f[61] = phi_f * k_f;
    phi_r = sc[7]*sc[19];
    Kc = Kc_save[61];
    k_r = k_f / Kc;
    q_r[61] = phi_r * k_r;

    /*reaction 63: CH2(S) + AR <=> CH2 + AR */
    phi_f = sc[8]*sc[20];
    k_f = k_f_save[62];
    q_f[62] = phi_f * k_f;
    phi_r = sc[7]*sc[20];
    Kc = Kc_save[62];
    k_r = k_f / Kc;
    q_r[62] = phi_r * k_r;

    /*reaction 64: CH2(S) + O2 <=> H + OH + CO */
    phi_f = sc[8]*sc[3];
    k_f = k_f_save[63];
    q_f[63] = phi_f * k_f;
    phi_r = sc[1]*sc[4]*sc[11];
    Kc = Kc_save[63];
    k_r = k_f / Kc;
    q_r[63] = phi_r * k_r;

    /*reaction 65: CH2(S) + O2 <=> CO + H2O */
    phi_f = sc[8]*sc[3];
    k_f = k_f_save[64];
    q_f[64] = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[64];
    k_r = k_f / Kc;
    q_r[64] = phi_r * k_r;

    /*reaction 66: CH2(S) + H2 <=> CH3 + H */
    phi_f = sc[8]*sc[0];
    k_f = k_f_save[65];
    q_f[65] = phi_f * k_f;
    phi_r = sc[9]*sc[1];
    Kc = Kc_save[65];
    k_r = k_f / Kc;
    q_r[65] = phi_r * k_r;

    /*reaction 67: CH2(S) + H2O <=> CH2 + H2O */
    phi_f = sc[8]*sc[5];
    k_f = k_f_save[66];
    q_f[66] = phi_f * k_f;
    phi_r = sc[7]*sc[5];
    Kc = Kc_save[66];
    k_r = k_f / Kc;
    q_r[66] = phi_r * k_r;

    /*reaction 68: CH2(S) + CH3 <=> H + C2H4 */
    phi_f = sc[8]*sc[9];
    k_f = k_f_save[67];
    q_f[67] = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_save[67];
    k_r = k_f / Kc;
    q_r[67] = phi_r * k_r;

    /*reaction 69: CH2(S) + CH4 <=> 2 CH3 */
    phi_f = sc[8]*sc[10];
    k_f = k_f_save[68];
    q_f[68] = phi_f * k_f;
    phi_r = sc[9]*sc[9];
    Kc = Kc_save[68];
    k_r = k_f / Kc;
    q_r[68] = phi_r * k_r;

    /*reaction 70: CH2(S) + CO <=> CH2 + CO */
    phi_f = sc[8]*sc[11];
    k_f = k_f_save[69];
    q_f[69] = phi_f * k_f;
    phi_r = sc[7]*sc[11];
    Kc = Kc_save[69];
    k_r = k_f / Kc;
    q_r[69] = phi_r * k_r;

    /*reaction 71: CH2(S) + CO2 <=> CH2 + CO2 */
    phi_f = sc[8]*sc[12];
    k_f = k_f_save[70];
    q_f[70] = phi_f * k_f;
    phi_r = sc[7]*sc[12];
    Kc = Kc_save[70];
    k_r = k_f / Kc;
    q_r[70] = phi_r * k_r;

    /*reaction 72: CH2(S) + CO2 <=> CO + CH2O */
    phi_f = sc[8]*sc[12];
    k_f = k_f_save[71];
    q_f[71] = phi_f * k_f;
    phi_r = sc[11]*sc[14];
    Kc = Kc_save[71];
    k_r = k_f / Kc;
    q_r[71] = phi_r * k_r;

    /*reaction 73: CH3 + O2 <=> O + CH3O */
    phi_f = sc[9]*sc[3];
    k_f = k_f_save[72];
    q_f[72] = phi_f * k_f;
    phi_r = sc[2]*sc[15];
    Kc = Kc_save[72];
    k_r = k_f / Kc;
    q_r[72] = phi_r * k_r;

    /*reaction 74: CH3 + O2 <=> OH + CH2O */
    phi_f = sc[9]*sc[3];
    k_f = k_f_save[73];
    q_f[73] = phi_f * k_f;
    phi_r = sc[4]*sc[14];
    Kc = Kc_save[73];
    k_r = k_f / Kc;
    q_r[73] = phi_r * k_r;

    /*reaction 75: 2 CH3 (+M) <=> C2H6 (+M) */
    phi_f = sc[9]*sc[9];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_save[74];
    redP = 1e-12 * alpha / k_f * 1.77e+50*exp(-9.67*tc[0]-3130.0076613053565779*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.4675*exp(T/-151))+ (0.5325*exp(T/-1038))+ (exp(-4970/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[74] = phi_f * k_f;
    phi_r = sc[18];
    Kc = Kc_save[74];
    k_r = k_f / Kc;
    q_r[74] = phi_r * k_r;

    /*reaction 76: 2 CH3 <=> H + C2H5 */
    phi_f = sc[9]*sc[9];
    k_f = k_f_save[75];
    q_f[75] = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[75];
    k_r = k_f / Kc;
    q_r[75] = phi_r * k_r;

    /*reaction 77: CH3 + HCO <=> CH4 + CO */
    phi_f = sc[9]*sc[13];
    k_f = k_f_save[76];
    q_f[76] = phi_f * k_f;
    phi_r = sc[10]*sc[11];
    Kc = Kc_save[76];
    k_r = k_f / Kc;
    q_r[76] = phi_r * k_r;

    /*reaction 78: CH3 + CH2O <=> HCO + CH4 */
    phi_f = sc[9]*sc[14];
    k_f = k_f_save[77];
    q_f[77] = phi_f * k_f;
    phi_r = sc[13]*sc[10];
    Kc = Kc_save[77];
    k_r = k_f / Kc;
    q_r[77] = phi_r * k_r;

    /*reaction 79: CH3 + C2H6 <=> C2H5 + CH4 */
    phi_f = sc[9]*sc[18];
    k_f = k_f_save[78];
    q_f[78] = phi_f * k_f;
    phi_r = sc[17]*sc[10];
    Kc = Kc_save[78];
    k_r = k_f / Kc;
    q_r[78] = phi_r * k_r;

    /*reaction 80: HCO + H2O <=> H + CO + H2O */
    phi_f = sc[13]*sc[5];
    k_f = k_f_save[79];
    q_f[79] = phi_f * k_f;
    phi_r = sc[1]*sc[11]*sc[5];
    Kc = Kc_save[79];
    k_r = k_f / Kc;
    q_r[79] = phi_r * k_r;

    /*reaction 81: HCO + M <=> H + CO + M */
    phi_f = sc[13];
    alpha = mixture + sc[0] + -1*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18];
    k_f = alpha * k_f_save[80];
    q_f[80] = phi_f * k_f;
    phi_r = sc[1]*sc[11];
    Kc = Kc_save[80];
    k_r = k_f / Kc;
    q_r[80] = phi_r * k_r;

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    phi_f = sc[13]*sc[3];
    k_f = k_f_save[81];
    q_f[81] = phi_f * k_f;
    phi_r = sc[6]*sc[11];
    Kc = Kc_save[81];
    k_r = k_f / Kc;
    q_r[81] = phi_r * k_r;

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    phi_f = sc[15]*sc[3];
    k_f = k_f_save[82];
    q_f[82] = phi_f * k_f;
    phi_r = sc[6]*sc[14];
    Kc = Kc_save[82];
    k_r = k_f / Kc;
    q_r[82] = phi_r * k_r;

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    phi_f = sc[17]*sc[3];
    k_f = k_f_save[83];
    q_f[83] = phi_f * k_f;
    phi_r = sc[6]*sc[16];
    Kc = Kc_save[83];
    k_r = k_f / Kc;
    q_r[83] = phi_r * k_r;

    return;
}


/*compute the equilibrium constants for each reaction */
void equilibriumConstants(double *kc, double * g_RT, double T)
{
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    /*reaction 1: O + H + M <=> OH + M */
    kc[0] = 1.0 / (refC) * exp((g_RT[2] + g_RT[1]) - (g_RT[4]));

    /*reaction 2: O + H2 <=> H + OH */
    kc[1] = exp((g_RT[2] + g_RT[0]) - (g_RT[1] + g_RT[4]));

    /*reaction 3: O + HO2 <=> OH + O2 */
    kc[2] = exp((g_RT[2] + g_RT[6]) - (g_RT[4] + g_RT[3]));

    /*reaction 4: O + CH2 <=> H + HCO */
    kc[3] = exp((g_RT[2] + g_RT[7]) - (g_RT[1] + g_RT[13]));

    /*reaction 5: O + CH2(S) <=> H + HCO */
    kc[4] = exp((g_RT[2] + g_RT[8]) - (g_RT[1] + g_RT[13]));

    /*reaction 6: O + CH3 <=> H + CH2O */
    kc[5] = exp((g_RT[2] + g_RT[9]) - (g_RT[1] + g_RT[14]));

    /*reaction 7: O + CH4 <=> OH + CH3 */
    kc[6] = exp((g_RT[2] + g_RT[10]) - (g_RT[4] + g_RT[9]));

    /*reaction 8: O + CO + M <=> CO2 + M */
    kc[7] = 1.0 / (refC) * exp((g_RT[2] + g_RT[11]) - (g_RT[12]));

    /*reaction 9: O + HCO <=> OH + CO */
    kc[8] = exp((g_RT[2] + g_RT[13]) - (g_RT[4] + g_RT[11]));

    /*reaction 10: O + HCO <=> H + CO2 */
    kc[9] = exp((g_RT[2] + g_RT[13]) - (g_RT[1] + g_RT[12]));

    /*reaction 11: O + CH2O <=> OH + HCO */
    kc[10] = exp((g_RT[2] + g_RT[14]) - (g_RT[4] + g_RT[13]));

    /*reaction 12: O + C2H4 <=> CH3 + HCO */
    kc[11] = exp((g_RT[2] + g_RT[16]) - (g_RT[9] + g_RT[13]));

    /*reaction 13: O + C2H5 <=> CH3 + CH2O */
    kc[12] = exp((g_RT[2] + g_RT[17]) - (g_RT[9] + g_RT[14]));

    /*reaction 14: O + C2H6 <=> OH + C2H5 */
    kc[13] = exp((g_RT[2] + g_RT[18]) - (g_RT[4] + g_RT[17]));

    /*reaction 15: O2 + CO <=> O + CO2 */
    kc[14] = exp((g_RT[3] + g_RT[11]) - (g_RT[2] + g_RT[12]));

    /*reaction 16: O2 + CH2O <=> HO2 + HCO */
    kc[15] = exp((g_RT[3] + g_RT[14]) - (g_RT[6] + g_RT[13]));

    /*reaction 17: H + O2 + M <=> HO2 + M */
    kc[16] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));

    /*reaction 18: H + 2 O2 <=> HO2 + O2 */
    kc[17] = 1.0 / (refC) * exp((g_RT[1] + 2 * g_RT[3]) - (g_RT[6] + g_RT[3]));

    /*reaction 19: H + O2 + H2O <=> HO2 + H2O */
    kc[18] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));

    /*reaction 20: H + O2 + N2 <=> HO2 + N2 */
    kc[19] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[19]) - (g_RT[6] + g_RT[19]));

    /*reaction 21: H + O2 + AR <=> HO2 + AR */
    kc[20] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[20]) - (g_RT[6] + g_RT[20]));

    /*reaction 22: H + O2 <=> O + OH */
    kc[21] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));

    /*reaction 23: 2 H + M <=> H2 + M */
    kc[22] = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));

    /*reaction 24: 2 H + H2 <=> 2 H2 */
    kc[23] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));

    /*reaction 25: 2 H + H2O <=> H2 + H2O */
    kc[24] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));

    /*reaction 26: 2 H + CO2 <=> H2 + CO2 */
    kc[25] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[12]) - (g_RT[0] + g_RT[12]));

    /*reaction 27: H + OH + M <=> H2O + M */
    kc[26] = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));

    /*reaction 28: H + HO2 <=> O2 + H2 */
    kc[27] = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));

    /*reaction 29: H + HO2 <=> 2 OH */
    kc[28] = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));

    /*reaction 30: H + CH2 (+M) <=> CH3 (+M) */
    kc[29] = 1.0 / (refC) * exp((g_RT[1] + g_RT[7]) - (g_RT[9]));

    /*reaction 31: H + CH3 (+M) <=> CH4 (+M) */
    kc[30] = 1.0 / (refC) * exp((g_RT[1] + g_RT[9]) - (g_RT[10]));

    /*reaction 32: H + CH4 <=> CH3 + H2 */
    kc[31] = exp((g_RT[1] + g_RT[10]) - (g_RT[9] + g_RT[0]));

    /*reaction 33: H + HCO (+M) <=> CH2O (+M) */
    kc[32] = 1.0 / (refC) * exp((g_RT[1] + g_RT[13]) - (g_RT[14]));

    /*reaction 34: H + HCO <=> H2 + CO */
    kc[33] = exp((g_RT[1] + g_RT[13]) - (g_RT[0] + g_RT[11]));

    /*reaction 35: H + CH2O (+M) <=> CH3O (+M) */
    kc[34] = 1.0 / (refC) * exp((g_RT[1] + g_RT[14]) - (g_RT[15]));

    /*reaction 36: H + CH2O <=> HCO + H2 */
    kc[35] = exp((g_RT[1] + g_RT[14]) - (g_RT[13] + g_RT[0]));

    /*reaction 37: H + CH3O <=> OH + CH3 */
    kc[36] = exp((g_RT[1] + g_RT[15]) - (g_RT[4] + g_RT[9]));

    /*reaction 38: H + C2H4 (+M) <=> C2H5 (+M) */
    kc[37] = 1.0 / (refC) * exp((g_RT[1] + g_RT[16]) - (g_RT[17]));

    /*reaction 39: H + C2H5 (+M) <=> C2H6 (+M) */
    kc[38] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[18]));

    /*reaction 40: H + C2H6 <=> C2H5 + H2 */
    kc[39] = exp((g_RT[1] + g_RT[18]) - (g_RT[17] + g_RT[0]));

    /*reaction 41: H2 + CO (+M) <=> CH2O (+M) */
    kc[40] = 1.0 / (refC) * exp((g_RT[0] + g_RT[11]) - (g_RT[14]));

    /*reaction 42: OH + H2 <=> H + H2O */
    kc[41] = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));

    /*reaction 43: 2 OH <=> O + H2O */
    kc[42] = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));

    /*reaction 44: OH + HO2 <=> O2 + H2O */
    kc[43] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));

    /*reaction 45: OH + CH2 <=> H + CH2O */
    kc[44] = exp((g_RT[4] + g_RT[7]) - (g_RT[1] + g_RT[14]));

    /*reaction 46: OH + CH2(S) <=> H + CH2O */
    kc[45] = exp((g_RT[4] + g_RT[8]) - (g_RT[1] + g_RT[14]));

    /*reaction 47: OH + CH3 <=> CH2 + H2O */
    kc[46] = exp((g_RT[4] + g_RT[9]) - (g_RT[7] + g_RT[5]));

    /*reaction 48: OH + CH3 <=> CH2(S) + H2O */
    kc[47] = exp((g_RT[4] + g_RT[9]) - (g_RT[8] + g_RT[5]));

    /*reaction 49: OH + CH4 <=> CH3 + H2O */
    kc[48] = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[5]));

    /*reaction 50: OH + CO <=> H + CO2 */
    kc[49] = exp((g_RT[4] + g_RT[11]) - (g_RT[1] + g_RT[12]));

    /*reaction 51: OH + HCO <=> H2O + CO */
    kc[50] = exp((g_RT[4] + g_RT[13]) - (g_RT[5] + g_RT[11]));

    /*reaction 52: OH + CH2O <=> HCO + H2O */
    kc[51] = exp((g_RT[4] + g_RT[14]) - (g_RT[13] + g_RT[5]));

    /*reaction 53: OH + C2H6 <=> C2H5 + H2O */
    kc[52] = exp((g_RT[4] + g_RT[18]) - (g_RT[17] + g_RT[5]));

    /*reaction 54: HO2 + CH2 <=> OH + CH2O */
    kc[53] = exp((g_RT[6] + g_RT[7]) - (g_RT[4] + g_RT[14]));

    /*reaction 55: HO2 + CH3 <=> O2 + CH4 */
    kc[54] = exp((g_RT[6] + g_RT[9]) - (g_RT[3] + g_RT[10]));

    /*reaction 56: HO2 + CH3 <=> OH + CH3O */
    kc[55] = exp((g_RT[6] + g_RT[9]) - (g_RT[4] + g_RT[15]));

    /*reaction 57: HO2 + CO <=> OH + CO2 */
    kc[56] = exp((g_RT[6] + g_RT[11]) - (g_RT[4] + g_RT[12]));

    /*reaction 58: CH2 + O2 <=> OH + HCO */
    kc[57] = exp((g_RT[7] + g_RT[3]) - (g_RT[4] + g_RT[13]));

    /*reaction 59: CH2 + H2 <=> H + CH3 */
    kc[58] = exp((g_RT[7] + g_RT[0]) - (g_RT[1] + g_RT[9]));

    /*reaction 60: CH2 + CH3 <=> H + C2H4 */
    kc[59] = exp((g_RT[7] + g_RT[9]) - (g_RT[1] + g_RT[16]));

    /*reaction 61: CH2 + CH4 <=> 2 CH3 */
    kc[60] = exp((g_RT[7] + g_RT[10]) - (2 * g_RT[9]));

    /*reaction 62: CH2(S) + N2 <=> CH2 + N2 */
    kc[61] = exp((g_RT[8] + g_RT[19]) - (g_RT[7] + g_RT[19]));

    /*reaction 63: CH2(S) + AR <=> CH2 + AR */
    kc[62] = exp((g_RT[8] + g_RT[20]) - (g_RT[7] + g_RT[20]));

    /*reaction 64: CH2(S) + O2 <=> H + OH + CO */
    kc[63] = refC * exp((g_RT[8] + g_RT[3]) - (g_RT[1] + g_RT[4] + g_RT[11]));

    /*reaction 65: CH2(S) + O2 <=> CO + H2O */
    kc[64] = exp((g_RT[8] + g_RT[3]) - (g_RT[11] + g_RT[5]));

    /*reaction 66: CH2(S) + H2 <=> CH3 + H */
    kc[65] = exp((g_RT[8] + g_RT[0]) - (g_RT[9] + g_RT[1]));

    /*reaction 67: CH2(S) + H2O <=> CH2 + H2O */
    kc[66] = exp((g_RT[8] + g_RT[5]) - (g_RT[7] + g_RT[5]));

    /*reaction 68: CH2(S) + CH3 <=> H + C2H4 */
    kc[67] = exp((g_RT[8] + g_RT[9]) - (g_RT[1] + g_RT[16]));

    /*reaction 69: CH2(S) + CH4 <=> 2 CH3 */
    kc[68] = exp((g_RT[8] + g_RT[10]) - (2 * g_RT[9]));

    /*reaction 70: CH2(S) + CO <=> CH2 + CO */
    kc[69] = exp((g_RT[8] + g_RT[11]) - (g_RT[7] + g_RT[11]));

    /*reaction 71: CH2(S) + CO2 <=> CH2 + CO2 */
    kc[70] = exp((g_RT[8] + g_RT[12]) - (g_RT[7] + g_RT[12]));

    /*reaction 72: CH2(S) + CO2 <=> CO + CH2O */
    kc[71] = exp((g_RT[8] + g_RT[12]) - (g_RT[11] + g_RT[14]));

    /*reaction 73: CH3 + O2 <=> O + CH3O */
    kc[72] = exp((g_RT[9] + g_RT[3]) - (g_RT[2] + g_RT[15]));

    /*reaction 74: CH3 + O2 <=> OH + CH2O */
    kc[73] = exp((g_RT[9] + g_RT[3]) - (g_RT[4] + g_RT[14]));

    /*reaction 75: 2 CH3 (+M) <=> C2H6 (+M) */
    kc[74] = 1.0 / (refC) * exp((2 * g_RT[9]) - (g_RT[18]));

    /*reaction 76: 2 CH3 <=> H + C2H5 */
    kc[75] = exp((2 * g_RT[9]) - (g_RT[1] + g_RT[17]));

    /*reaction 77: CH3 + HCO <=> CH4 + CO */
    kc[76] = exp((g_RT[9] + g_RT[13]) - (g_RT[10] + g_RT[11]));

    /*reaction 78: CH3 + CH2O <=> HCO + CH4 */
    kc[77] = exp((g_RT[9] + g_RT[14]) - (g_RT[13] + g_RT[10]));

    /*reaction 79: CH3 + C2H6 <=> C2H5 + CH4 */
    kc[78] = exp((g_RT[9] + g_RT[18]) - (g_RT[17] + g_RT[10]));

    /*reaction 80: HCO + H2O <=> H + CO + H2O */
    kc[79] = refC * exp((g_RT[13] + g_RT[5]) - (g_RT[1] + g_RT[11] + g_RT[5]));

    /*reaction 81: HCO + M <=> H + CO + M */
    kc[80] = refC * exp((g_RT[13]) - (g_RT[1] + g_RT[11]));

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    kc[81] = exp((g_RT[13] + g_RT[3]) - (g_RT[6] + g_RT[11]));

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    kc[82] = exp((g_RT[15] + g_RT[3]) - (g_RT[6] + g_RT[14]));

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    kc[83] = exp((g_RT[17] + g_RT[3]) - (g_RT[6] + g_RT[16]));

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
        /*species 0: H2 */
        species[0] =
            -9.179351730000000e+02 * invT
            +1.661320882000000e+00
            -2.344331120000000e+00 * tc[0]
            -3.990260375000000e-03 * tc[1]
            +3.246358500000000e-06 * tc[2]
            -1.679767450000000e-09 * tc[3]
            +3.688058805000000e-13 * tc[4];
        /*species 1: H */
        species[1] =
            +2.547365990000000e+04 * invT
            +2.946682853000000e+00
            -2.500000000000000e+00 * tc[0]
            -3.526664095000000e-13 * tc[1]
            +3.326532733333333e-16 * tc[2]
            -1.917346933333333e-19 * tc[3]
            +4.638661660000000e-23 * tc[4];
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
            +1.498367080000000e-03 * tc[1]
            -1.641217001666667e-06 * tc[2]
            +8.067745908333334e-10 * tc[3]
            -1.621864185000000e-13 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.615080560000000e+03 * invT
            +4.095940888000000e+00
            -3.992015430000000e+00 * tc[0]
            +1.200658760000000e-03 * tc[1]
            -7.696564016666666e-07 * tc[2]
            +3.234277775000000e-10 * tc[3]
            -6.820573500000000e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.029372670000000e+04 * invT
            +5.047672768000000e+00
            -4.198640560000000e+00 * tc[0]
            +1.018217050000000e-03 * tc[1]
            -1.086733685000000e-06 * tc[2]
            +4.573308850000000e-10 * tc[3]
            -8.859890850000000e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +2.948080400000000e+02 * invT
            +5.851355599999999e-01
            -4.301798010000000e+00 * tc[0]
            +2.374560255000000e-03 * tc[1]
            -3.526381516666666e-06 * tc[2]
            +2.023032450000000e-09 * tc[3]
            -4.646125620000001e-13 * tc[4];
        /*species 7: CH2 */
        species[7] =
            +4.600404010000000e+04 * invT
            +2.200146820000000e+00
            -3.762678670000000e+00 * tc[0]
            -4.844360715000000e-04 * tc[1]
            -4.658164016666667e-07 * tc[2]
            +3.209092941666667e-10 * tc[3]
            -8.437085950000000e-14 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +5.049681630000000e+04 * invT
            +4.967723077000000e+00
            -4.198604110000000e+00 * tc[0]
            +1.183307095000000e-03 * tc[1]
            -1.372160366666667e-06 * tc[2]
            +5.573466508333334e-10 * tc[3]
            -9.715736850000000e-14 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +1.644499880000000e+04 * invT
            +2.069026070000000e+00
            -3.673590400000000e+00 * tc[0]
            -1.005475875000000e-03 * tc[1]
            -9.550364266666668e-07 * tc[2]
            +5.725978541666666e-10 * tc[3]
            -1.271928670000000e-13 * tc[4];
        /*species 10: CH4 */
        species[10] =
            -1.024664760000000e+04 * invT
            +9.791179889999999e+00
            -5.149876130000000e+00 * tc[0]
            +6.835489400000000e-03 * tc[1]
            -8.196676650000000e-06 * tc[2]
            +4.039525216666667e-09 * tc[3]
            -8.334697800000000e-13 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.434408600000000e+04 * invT
            +7.112418999999992e-02
            -3.579533470000000e+00 * tc[0]
            +3.051768400000000e-04 * tc[1]
            -1.694690550000000e-07 * tc[2]
            -7.558382366666667e-11 * tc[3]
            +4.522122495000000e-14 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.837196970000000e+04 * invT
            -7.544278700000000e+00
            -2.356773520000000e+00 * tc[0]
            -4.492298385000000e-03 * tc[1]
            +1.187260448333333e-06 * tc[2]
            -2.049325183333333e-10 * tc[3]
            +7.184977399999999e-15 * tc[4];
        /*species 13: HCO */
        species[13] =
            +3.839564960000000e+03 * invT
            +8.268134100000002e-01
            -4.221185840000000e+00 * tc[0]
            +1.621962660000000e-03 * tc[1]
            -2.296657433333333e-06 * tc[2]
            +1.109534108333333e-09 * tc[3]
            -2.168844325000000e-13 * tc[4];
        /*species 14: CH2O */
        species[14] =
            -1.430895670000000e+04 * invT
            +4.190910250000000e+00
            -4.793723150000000e+00 * tc[0]
            +4.954166845000000e-03 * tc[1]
            -6.220333466666666e-06 * tc[2]
            +3.160710508333333e-09 * tc[3]
            -6.588632600000000e-13 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +9.786011000000000e+02 * invT
            -1.104597300000000e+01
            -2.106204000000000e+00 * tc[0]
            -3.608297500000000e-03 * tc[1]
            -8.897453333333333e-07 * tc[2]
            +6.148030000000000e-10 * tc[3]
            -1.037805000000000e-13 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +5.089775930000000e+03 * invT
            -1.381294799999999e-01
            -3.959201480000000e+00 * tc[0]
            +3.785261235000000e-03 * tc[1]
            -9.516504866666667e-06 * tc[2]
            +5.763239608333333e-09 * tc[3]
            -1.349421865000000e-12 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +1.284162650000000e+04 * invT
            -4.007435600000004e-01
            -4.306465680000000e+00 * tc[0]
            +2.093294460000000e-03 * tc[1]
            -8.285713450000000e-06 * tc[2]
            +4.992721716666666e-09 * tc[3]
            -1.152545020000000e-12 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            -1.152220550000000e+04 * invT
            +1.624601760000000e+00
            -4.291424920000000e+00 * tc[0]
            +2.750771350000000e-03 * tc[1]
            -9.990638133333334e-06 * tc[2]
            +5.903885708333334e-09 * tc[3]
            -1.343428855000000e-12 * tc[4];
        /*species 19: N2 */
        species[19] =
            -1.020899900000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041202000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427000000000e-13 * tc[4];
        /*species 20: AR */
        species[20] =
            -7.453750000000000e+02 * invT
            -1.866000000000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -9.501589220000000e+02 * invT
            +6.542302510000000e+00
            -3.337279200000000e+00 * tc[0]
            +2.470123655000000e-05 * tc[1]
            -8.324279633333333e-08 * tc[2]
            +1.496386616666667e-11 * tc[3]
            -1.001276880000000e-15 * tc[4];
        /*species 1: H */
        species[1] =
            +2.547365990000000e+04 * invT
            +2.946682924000000e+00
            -2.500000010000000e+00 * tc[0]
            +1.154214865000000e-11 * tc[1]
            -2.692699133333334e-15 * tc[2]
            +3.945960291666667e-19 * tc[3]
            -2.490986785000000e-23 * tc[4];
        /*species 2: O */
        species[2] =
            +2.921757910000000e+04 * invT
            -2.214917859999999e+00
            -2.569420780000000e+00 * tc[0]
            +4.298705685000000e-05 * tc[1]
            -6.991409816666667e-09 * tc[2]
            +8.348149916666666e-13 * tc[3]
            -6.141684549999999e-17 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.088457720000000e+03 * invT
            -2.170693450000000e+00
            -3.282537840000000e+00 * tc[0]
            -7.415437700000000e-04 * tc[1]
            +1.263277781666667e-07 * tc[2]
            -1.745587958333333e-11 * tc[3]
            +1.083588970000000e-15 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.858657000000000e+03 * invT
            -1.383808430000000e+00
            -3.092887670000000e+00 * tc[0]
            -2.742148580000000e-04 * tc[1]
            -2.108420466666667e-08 * tc[2]
            +7.328846300000000e-12 * tc[3]
            -5.870618800000000e-16 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.000429710000000e+04 * invT
            -1.932777610000000e+00
            -3.033992490000000e+00 * tc[0]
            -1.088459020000000e-03 * tc[1]
            +2.734541966666666e-08 * tc[2]
            +8.086832250000000e-12 * tc[3]
            -8.410049600000000e-16 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.118567130000000e+02 * invT
            +2.321087500000001e-01
            -4.017210900000000e+00 * tc[0]
            -1.119910065000000e-03 * tc[1]
            +1.056096916666667e-07 * tc[2]
            -9.520530833333334e-12 * tc[3]
            +5.395426750000000e-16 * tc[4];
        /*species 7: CH2 */
        species[7] =
            +4.626360400000000e+04 * invT
            -3.297092110000000e+00
            -2.874101130000000e+00 * tc[0]
            -1.828196460000000e-03 * tc[1]
            +2.348243283333333e-07 * tc[2]
            -2.168162908333333e-11 * tc[3]
            +9.386378350000000e-16 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +5.092599970000000e+04 * invT
            -6.334463270000000e+00
            -2.292038420000000e+00 * tc[0]
            -2.327943185000000e-03 * tc[1]
            +3.353199116666667e-07 * tc[2]
            -3.482550000000000e-11 * tc[3]
            +1.698581825000000e-15 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +1.677558430000000e+04 * invT
            -6.194354070000000e+00
            -2.285717720000000e+00 * tc[0]
            -3.619950185000000e-03 * tc[1]
            +4.978572466666667e-07 * tc[2]
            -4.964038700000000e-11 * tc[3]
            +2.335771970000000e-15 * tc[4];
        /*species 10: CH4 */
        species[10] =
            -9.468344590000001e+03 * invT
            -1.836246650500000e+01
            -7.485149500000000e-02 * tc[0]
            -6.695473350000000e-03 * tc[1]
            +9.554763483333333e-07 * tc[2]
            -1.019104458333333e-10 * tc[3]
            +5.090761500000000e-15 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.415187240000000e+04 * invT
            -5.103502110000000e+00
            -2.715185610000000e+00 * tc[0]
            -1.031263715000000e-03 * tc[1]
            +1.664709618333334e-07 * tc[2]
            -1.917108400000000e-11 * tc[3]
            +1.018238580000000e-15 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.875916600000000e+04 * invT
            +1.585822230000000e+00
            -3.857460290000000e+00 * tc[0]
            -2.207185130000000e-03 * tc[1]
            +3.691356733333334e-07 * tc[2]
            -4.362418233333334e-11 * tc[3]
            +2.360420820000000e-15 * tc[4];
        /*species 13: HCO */
        species[13] =
            +4.011918150000000e+03 * invT
            -7.026170540000000e+00
            -2.772174380000000e+00 * tc[0]
            -2.478477630000000e-03 * tc[1]
            +4.140760216666667e-07 * tc[2]
            -4.909681483333334e-11 * tc[3]
            +2.667543555000000e-15 * tc[4];
        /*species 14: CH2O */
        species[14] =
            -1.399583230000000e+04 * invT
            -1.189563292000000e+01
            -1.760690080000000e+00 * tc[0]
            -4.600000410000000e-03 * tc[1]
            +7.370980216666666e-07 * tc[2]
            -8.386767666666666e-11 * tc[3]
            +4.419278200000001e-15 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +1.278325200000000e+02 * invT
            +8.412240000000000e-01
            -3.770799000000000e+00 * tc[0]
            -3.935748500000000e-03 * tc[1]
            +4.427306666666667e-07 * tc[2]
            -3.287025833333333e-11 * tc[3]
            +1.056308000000000e-15 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +4.939886140000000e+03 * invT
            -8.269258140000002e+00
            -2.036111160000000e+00 * tc[0]
            -7.322707550000000e-03 * tc[1]
            +1.118463191666667e-06 * tc[2]
            -1.226857691666667e-10 * tc[3]
            +6.285303050000000e-15 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +1.285752000000000e+04 * invT
            -1.150777788000000e+01
            -1.954656420000000e+00 * tc[0]
            -8.698636100000001e-03 * tc[1]
            +1.330344446666667e-06 * tc[2]
            -1.460147408333333e-10 * tc[3]
            +7.482078800000000e-15 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            -1.142639320000000e+04 * invT
            -1.404372920000000e+01
            -1.071881500000000e+00 * tc[0]
            -1.084263385000000e-02 * tc[1]
            +1.670934450000000e-06 * tc[2]
            -1.845100008333333e-10 * tc[3]
            +9.500144500000000e-15 * tc[4];
        /*species 19: N2 */
        species[19] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439884000000000e-04 * tc[1]
            +9.474600000000001e-08 * tc[2]
            -8.414198333333333e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
        /*species 20: AR */
        species[20] =
            -7.453750000000000e+02 * invT
            -1.866000000000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
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
        /*species 0: H2 */
        species[0] =
            -9.17935173e+02 * invT
            +6.61320882e-01
            -2.34433112e+00 * tc[0]
            -3.99026037e-03 * tc[1]
            +3.24635850e-06 * tc[2]
            -1.67976745e-09 * tc[3]
            +3.68805881e-13 * tc[4];
        /*species 1: H */
        species[1] =
            +2.54736599e+04 * invT
            +1.94668285e+00
            -2.50000000e+00 * tc[0]
            -3.52666409e-13 * tc[1]
            +3.32653273e-16 * tc[2]
            -1.91734693e-19 * tc[3]
            +4.63866166e-23 * tc[4];
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
            +1.49836708e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774591e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.61508056e+03 * invT
            +3.09594089e+00
            -3.99201543e+00 * tc[0]
            +1.20065876e-03 * tc[1]
            -7.69656402e-07 * tc[2]
            +3.23427778e-10 * tc[3]
            -6.82057350e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.02937267e+04 * invT
            +4.04767277e+00
            -4.19864056e+00 * tc[0]
            +1.01821705e-03 * tc[1]
            -1.08673369e-06 * tc[2]
            +4.57330885e-10 * tc[3]
            -8.85989085e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +2.94808040e+02 * invT
            -4.14864440e-01
            -4.30179801e+00 * tc[0]
            +2.37456025e-03 * tc[1]
            -3.52638152e-06 * tc[2]
            +2.02303245e-09 * tc[3]
            -4.64612562e-13 * tc[4];
        /*species 7: CH2 */
        species[7] =
            +4.60040401e+04 * invT
            +1.20014682e+00
            -3.76267867e+00 * tc[0]
            -4.84436072e-04 * tc[1]
            -4.65816402e-07 * tc[2]
            +3.20909294e-10 * tc[3]
            -8.43708595e-14 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +5.04968163e+04 * invT
            +3.96772308e+00
            -4.19860411e+00 * tc[0]
            +1.18330710e-03 * tc[1]
            -1.37216037e-06 * tc[2]
            +5.57346651e-10 * tc[3]
            -9.71573685e-14 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +1.64449988e+04 * invT
            +1.06902607e+00
            -3.67359040e+00 * tc[0]
            -1.00547588e-03 * tc[1]
            -9.55036427e-07 * tc[2]
            +5.72597854e-10 * tc[3]
            -1.27192867e-13 * tc[4];
        /*species 10: CH4 */
        species[10] =
            -1.02466476e+04 * invT
            +8.79117989e+00
            -5.14987613e+00 * tc[0]
            +6.83548940e-03 * tc[1]
            -8.19667665e-06 * tc[2]
            +4.03952522e-09 * tc[3]
            -8.33469780e-13 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.43440860e+04 * invT
            -9.28875810e-01
            -3.57953347e+00 * tc[0]
            +3.05176840e-04 * tc[1]
            -1.69469055e-07 * tc[2]
            -7.55838237e-11 * tc[3]
            +4.52212249e-14 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.83719697e+04 * invT
            -8.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
        /*species 13: HCO */
        species[13] =
            +3.83956496e+03 * invT
            -1.73186590e-01
            -4.22118584e+00 * tc[0]
            +1.62196266e-03 * tc[1]
            -2.29665743e-06 * tc[2]
            +1.10953411e-09 * tc[3]
            -2.16884432e-13 * tc[4];
        /*species 14: CH2O */
        species[14] =
            -1.43089567e+04 * invT
            +3.19091025e+00
            -4.79372315e+00 * tc[0]
            +4.95416684e-03 * tc[1]
            -6.22033347e-06 * tc[2]
            +3.16071051e-09 * tc[3]
            -6.58863260e-13 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +9.78601100e+02 * invT
            -1.20459730e+01
            -2.10620400e+00 * tc[0]
            -3.60829750e-03 * tc[1]
            -8.89745333e-07 * tc[2]
            +6.14803000e-10 * tc[3]
            -1.03780500e-13 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +5.08977593e+03 * invT
            -1.13812948e+00
            -3.95920148e+00 * tc[0]
            +3.78526124e-03 * tc[1]
            -9.51650487e-06 * tc[2]
            +5.76323961e-09 * tc[3]
            -1.34942187e-12 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +1.28416265e+04 * invT
            -1.40074356e+00
            -4.30646568e+00 * tc[0]
            +2.09329446e-03 * tc[1]
            -8.28571345e-06 * tc[2]
            +4.99272172e-09 * tc[3]
            -1.15254502e-12 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            -1.15222055e+04 * invT
            +6.24601760e-01
            -4.29142492e+00 * tc[0]
            +2.75077135e-03 * tc[1]
            -9.99063813e-06 * tc[2]
            +5.90388571e-09 * tc[3]
            -1.34342886e-12 * tc[4];
        /*species 19: N2 */
        species[19] =
            -1.02089990e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
        /*species 20: AR */
        species[20] =
            -7.45375000e+02 * invT
            -2.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -9.50158922e+02 * invT
            +5.54230251e+00
            -3.33727920e+00 * tc[0]
            +2.47012365e-05 * tc[1]
            -8.32427963e-08 * tc[2]
            +1.49638662e-11 * tc[3]
            -1.00127688e-15 * tc[4];
        /*species 1: H */
        species[1] =
            +2.54736599e+04 * invT
            +1.94668292e+00
            -2.50000001e+00 * tc[0]
            +1.15421486e-11 * tc[1]
            -2.69269913e-15 * tc[2]
            +3.94596029e-19 * tc[3]
            -2.49098679e-23 * tc[4];
        /*species 2: O */
        species[2] =
            +2.92175791e+04 * invT
            -3.21491786e+00
            -2.56942078e+00 * tc[0]
            +4.29870569e-05 * tc[1]
            -6.99140982e-09 * tc[2]
            +8.34814992e-13 * tc[3]
            -6.14168455e-17 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.08845772e+03 * invT
            -3.17069345e+00
            -3.28253784e+00 * tc[0]
            -7.41543770e-04 * tc[1]
            +1.26327778e-07 * tc[2]
            -1.74558796e-11 * tc[3]
            +1.08358897e-15 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.85865700e+03 * invT
            -2.38380843e+00
            -3.09288767e+00 * tc[0]
            -2.74214858e-04 * tc[1]
            -2.10842047e-08 * tc[2]
            +7.32884630e-12 * tc[3]
            -5.87061880e-16 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.00042971e+04 * invT
            -2.93277761e+00
            -3.03399249e+00 * tc[0]
            -1.08845902e-03 * tc[1]
            +2.73454197e-08 * tc[2]
            +8.08683225e-12 * tc[3]
            -8.41004960e-16 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.11856713e+02 * invT
            -7.67891250e-01
            -4.01721090e+00 * tc[0]
            -1.11991006e-03 * tc[1]
            +1.05609692e-07 * tc[2]
            -9.52053083e-12 * tc[3]
            +5.39542675e-16 * tc[4];
        /*species 7: CH2 */
        species[7] =
            +4.62636040e+04 * invT
            -4.29709211e+00
            -2.87410113e+00 * tc[0]
            -1.82819646e-03 * tc[1]
            +2.34824328e-07 * tc[2]
            -2.16816291e-11 * tc[3]
            +9.38637835e-16 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +5.09259997e+04 * invT
            -7.33446327e+00
            -2.29203842e+00 * tc[0]
            -2.32794318e-03 * tc[1]
            +3.35319912e-07 * tc[2]
            -3.48255000e-11 * tc[3]
            +1.69858182e-15 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +1.67755843e+04 * invT
            -7.19435407e+00
            -2.28571772e+00 * tc[0]
            -3.61995018e-03 * tc[1]
            +4.97857247e-07 * tc[2]
            -4.96403870e-11 * tc[3]
            +2.33577197e-15 * tc[4];
        /*species 10: CH4 */
        species[10] =
            -9.46834459e+03 * invT
            -1.93624665e+01
            -7.48514950e-02 * tc[0]
            -6.69547335e-03 * tc[1]
            +9.55476348e-07 * tc[2]
            -1.01910446e-10 * tc[3]
            +5.09076150e-15 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.41518724e+04 * invT
            -6.10350211e+00
            -2.71518561e+00 * tc[0]
            -1.03126372e-03 * tc[1]
            +1.66470962e-07 * tc[2]
            -1.91710840e-11 * tc[3]
            +1.01823858e-15 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.87591660e+04 * invT
            +5.85822230e-01
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
        /*species 13: HCO */
        species[13] =
            +4.01191815e+03 * invT
            -8.02617054e+00
            -2.77217438e+00 * tc[0]
            -2.47847763e-03 * tc[1]
            +4.14076022e-07 * tc[2]
            -4.90968148e-11 * tc[3]
            +2.66754356e-15 * tc[4];
        /*species 14: CH2O */
        species[14] =
            -1.39958323e+04 * invT
            -1.28956329e+01
            -1.76069008e+00 * tc[0]
            -4.60000041e-03 * tc[1]
            +7.37098022e-07 * tc[2]
            -8.38676767e-11 * tc[3]
            +4.41927820e-15 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +1.27832520e+02 * invT
            -1.58776000e-01
            -3.77079900e+00 * tc[0]
            -3.93574850e-03 * tc[1]
            +4.42730667e-07 * tc[2]
            -3.28702583e-11 * tc[3]
            +1.05630800e-15 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +4.93988614e+03 * invT
            -9.26925814e+00
            -2.03611116e+00 * tc[0]
            -7.32270755e-03 * tc[1]
            +1.11846319e-06 * tc[2]
            -1.22685769e-10 * tc[3]
            +6.28530305e-15 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +1.28575200e+04 * invT
            -1.25077779e+01
            -1.95465642e+00 * tc[0]
            -8.69863610e-03 * tc[1]
            +1.33034445e-06 * tc[2]
            -1.46014741e-10 * tc[3]
            +7.48207880e-15 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            -1.14263932e+04 * invT
            -1.50437292e+01
            -1.07188150e+00 * tc[0]
            -1.08426339e-02 * tc[1]
            +1.67093445e-06 * tc[2]
            -1.84510001e-10 * tc[3]
            +9.50014450e-15 * tc[4];
        /*species 19: N2 */
        species[19] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 20: AR */
        species[20] =
            -7.45375000e+02 * invT
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
    double T = tc[1], invT = 1.0 / T;

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
        /*species 7: CH2 */
        species[7] =
            +2.76267867e+00
            +9.68872143e-04 * tc[1]
            +2.79489841e-06 * tc[2]
            -3.85091153e-09 * tc[3]
            +1.68741719e-12 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +3.19860411e+00
            -2.36661419e-03 * tc[1]
            +8.23296220e-06 * tc[2]
            -6.68815981e-09 * tc[3]
            +1.94314737e-12 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +2.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 10: CH4 */
        species[10] =
            +4.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 11: CO */
        species[11] =
            +2.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +1.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 13: HCO */
        species[13] =
            +3.22118584e+00
            -3.24392532e-03 * tc[1]
            +1.37799446e-05 * tc[2]
            -1.33144093e-08 * tc[3]
            +4.33768865e-12 * tc[4];
        /*species 14: CH2O */
        species[14] =
            +3.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +1.10620400e+00
            +7.21659500e-03 * tc[1]
            +5.33847200e-06 * tc[2]
            -7.37763600e-09 * tc[3]
            +2.07561000e-12 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +2.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +3.30646568e+00
            -4.18658892e-03 * tc[1]
            +4.97142807e-05 * tc[2]
            -5.99126606e-08 * tc[3]
            +2.30509004e-11 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            +3.29142492e+00
            -5.50154270e-03 * tc[1]
            +5.99438288e-05 * tc[2]
            -7.08466285e-08 * tc[3]
            +2.68685771e-11 * tc[4];
        /*species 19: N2 */
        species[19] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 20: AR */
        species[20] =
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
        /*species 7: CH2 */
        species[7] =
            +1.87410113e+00
            +3.65639292e-03 * tc[1]
            -1.40894597e-06 * tc[2]
            +2.60179549e-10 * tc[3]
            -1.87727567e-14 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +1.29203842e+00
            +4.65588637e-03 * tc[1]
            -2.01191947e-06 * tc[2]
            +4.17906000e-10 * tc[3]
            -3.39716365e-14 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +1.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 10: CH4 */
        species[10] =
            -9.25148505e-01
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 11: CO */
        species[11] =
            +1.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +2.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 13: HCO */
        species[13] =
            +1.77217438e+00
            +4.95695526e-03 * tc[1]
            -2.48445613e-06 * tc[2]
            +5.89161778e-10 * tc[3]
            -5.33508711e-14 * tc[4];
        /*species 14: CH2O */
        species[14] =
            +7.60690080e-01
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +2.77079900e+00
            +7.87149700e-03 * tc[1]
            -2.65638400e-06 * tc[2]
            +3.94443100e-10 * tc[3]
            -2.11261600e-14 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +1.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +9.54656420e-01
            +1.73972722e-02 * tc[1]
            -7.98206668e-06 * tc[2]
            +1.75217689e-09 * tc[3]
            -1.49641576e-13 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            +7.18815000e-02
            +2.16852677e-02 * tc[1]
            -1.00256067e-05 * tc[2]
            +2.21412001e-09 * tc[3]
            -1.90002890e-13 * tc[4];
        /*species 19: N2 */
        species[19] =
            +1.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 20: AR */
        species[20] =
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
    double T = tc[1], invT = 1.0 / T;

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
        /*species 7: CH2 */
        species[7] =
            +3.76267867e+00
            +9.68872143e-04 * tc[1]
            +2.79489841e-06 * tc[2]
            -3.85091153e-09 * tc[3]
            +1.68741719e-12 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +4.19860411e+00
            -2.36661419e-03 * tc[1]
            +8.23296220e-06 * tc[2]
            -6.68815981e-09 * tc[3]
            +1.94314737e-12 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +3.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 10: CH4 */
        species[10] =
            +5.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 11: CO */
        species[11] =
            +3.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +2.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 13: HCO */
        species[13] =
            +4.22118584e+00
            -3.24392532e-03 * tc[1]
            +1.37799446e-05 * tc[2]
            -1.33144093e-08 * tc[3]
            +4.33768865e-12 * tc[4];
        /*species 14: CH2O */
        species[14] =
            +4.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +2.10620400e+00
            +7.21659500e-03 * tc[1]
            +5.33847200e-06 * tc[2]
            -7.37763600e-09 * tc[3]
            +2.07561000e-12 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +3.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +4.30646568e+00
            -4.18658892e-03 * tc[1]
            +4.97142807e-05 * tc[2]
            -5.99126606e-08 * tc[3]
            +2.30509004e-11 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            +4.29142492e+00
            -5.50154270e-03 * tc[1]
            +5.99438288e-05 * tc[2]
            -7.08466285e-08 * tc[3]
            +2.68685771e-11 * tc[4];
        /*species 19: N2 */
        species[19] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 20: AR */
        species[20] =
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
        /*species 7: CH2 */
        species[7] =
            +2.87410113e+00
            +3.65639292e-03 * tc[1]
            -1.40894597e-06 * tc[2]
            +2.60179549e-10 * tc[3]
            -1.87727567e-14 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +2.29203842e+00
            +4.65588637e-03 * tc[1]
            -2.01191947e-06 * tc[2]
            +4.17906000e-10 * tc[3]
            -3.39716365e-14 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +2.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 10: CH4 */
        species[10] =
            +7.48514950e-02
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 11: CO */
        species[11] =
            +2.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +3.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 13: HCO */
        species[13] =
            +2.77217438e+00
            +4.95695526e-03 * tc[1]
            -2.48445613e-06 * tc[2]
            +5.89161778e-10 * tc[3]
            -5.33508711e-14 * tc[4];
        /*species 14: CH2O */
        species[14] =
            +1.76069008e+00
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +3.77079900e+00
            +7.87149700e-03 * tc[1]
            -2.65638400e-06 * tc[2]
            +3.94443100e-10 * tc[3]
            -2.11261600e-14 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +2.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +1.95465642e+00
            +1.73972722e-02 * tc[1]
            -7.98206668e-06 * tc[2]
            +1.75217689e-09 * tc[3]
            -1.49641576e-13 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            +1.07188150e+00
            +2.16852677e-02 * tc[1]
            -1.00256067e-05 * tc[2]
            +2.21412001e-09 * tc[3]
            -1.90002890e-13 * tc[4];
        /*species 19: N2 */
        species[19] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 20: AR */
        species[20] =
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
    double T = tc[1], invT = 1.0 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +1.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 * invT;
        /*species 1: H */
        species[1] =
            +1.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 * invT;
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
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 4: OH */
        species[4] =
            +2.99201543e+00
            -1.20065876e-03 * tc[1]
            +1.53931280e-06 * tc[2]
            -9.70283332e-10 * tc[3]
            +2.72822940e-13 * tc[4]
            +3.61508056e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +3.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 7: CH2 */
        species[7] =
            +2.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 * invT;
        /*species 8: CH2(S) */
        species[8] =
            +3.19860411e+00
            -1.18330710e-03 * tc[1]
            +2.74432073e-06 * tc[2]
            -1.67203995e-09 * tc[3]
            +3.88629474e-13 * tc[4]
            +5.04968163e+04 * invT;
        /*species 9: CH3 */
        species[9] =
            +2.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 * invT;
        /*species 10: CH4 */
        species[10] =
            +4.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 11: CO */
        species[11] =
            +2.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +1.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 13: HCO */
        species[13] =
            +3.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 * invT;
        /*species 14: CH2O */
        species[14] =
            +3.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 * invT;
        /*species 15: CH3O */
        species[15] =
            +1.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122000e-13 * tc[4]
            +9.78601100e+02 * invT;
        /*species 16: C2H4 */
        species[16] =
            +2.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 * invT;
        /*species 17: C2H5 */
        species[17] =
            +3.30646568e+00
            -2.09329446e-03 * tc[1]
            +1.65714269e-05 * tc[2]
            -1.49781651e-08 * tc[3]
            +4.61018008e-12 * tc[4]
            +1.28416265e+04 * invT;
        /*species 18: C2H6 */
        species[18] =
            +3.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 * invT;
        /*species 19: N2 */
        species[19] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
        /*species 20: AR */
        species[20] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +2.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 * invT;
        /*species 1: H */
        species[1] =
            +1.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 * invT;
        /*species 2: O */
        species[2] =
            +1.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +2.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 4: OH */
        species[4] =
            +2.09288767e+00
            +2.74214858e-04 * tc[1]
            +4.21684093e-08 * tc[2]
            -2.19865389e-11 * tc[3]
            +2.34824752e-15 * tc[4]
            +3.85865700e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +2.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 7: CH2 */
        species[7] =
            +1.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 * invT;
        /*species 8: CH2(S) */
        species[8] =
            +1.29203842e+00
            +2.32794318e-03 * tc[1]
            -6.70639823e-07 * tc[2]
            +1.04476500e-10 * tc[3]
            -6.79432730e-15 * tc[4]
            +5.09259997e+04 * invT;
        /*species 9: CH3 */
        species[9] =
            +1.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 * invT;
        /*species 10: CH4 */
        species[10] =
            -9.25148505e-01
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 11: CO */
        species[11] =
            +1.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +2.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 13: HCO */
        species[13] =
            +1.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 * invT;
        /*species 14: CH2O */
        species[14] =
            +7.60690080e-01
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 * invT;
        /*species 15: CH3O */
        species[15] =
            +2.77079900e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832520e+02 * invT;
        /*species 16: C2H4 */
        species[16] =
            +1.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 * invT;
        /*species 17: C2H5 */
        species[17] =
            +9.54656420e-01
            +8.69863610e-03 * tc[1]
            -2.66068889e-06 * tc[2]
            +4.38044223e-10 * tc[3]
            -2.99283152e-14 * tc[4]
            +1.28575200e+04 * invT;
        /*species 18: C2H6 */
        species[18] =
            +7.18815000e-02
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 * invT;
        /*species 19: N2 */
        species[19] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 20: AR */
        species[20] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
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
        /*species 0: H2 */
        species[0] =
            +2.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 * invT;
        /*species 1: H */
        species[1] =
            +2.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 * invT;
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
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 4: OH */
        species[4] =
            +3.99201543e+00
            -1.20065876e-03 * tc[1]
            +1.53931280e-06 * tc[2]
            -9.70283332e-10 * tc[3]
            +2.72822940e-13 * tc[4]
            +3.61508056e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +4.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 7: CH2 */
        species[7] =
            +3.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 * invT;
        /*species 8: CH2(S) */
        species[8] =
            +4.19860411e+00
            -1.18330710e-03 * tc[1]
            +2.74432073e-06 * tc[2]
            -1.67203995e-09 * tc[3]
            +3.88629474e-13 * tc[4]
            +5.04968163e+04 * invT;
        /*species 9: CH3 */
        species[9] =
            +3.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 * invT;
        /*species 10: CH4 */
        species[10] =
            +5.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 11: CO */
        species[11] =
            +3.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +2.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 13: HCO */
        species[13] =
            +4.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 * invT;
        /*species 14: CH2O */
        species[14] =
            +4.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 * invT;
        /*species 15: CH3O */
        species[15] =
            +2.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122000e-13 * tc[4]
            +9.78601100e+02 * invT;
        /*species 16: C2H4 */
        species[16] =
            +3.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 * invT;
        /*species 17: C2H5 */
        species[17] =
            +4.30646568e+00
            -2.09329446e-03 * tc[1]
            +1.65714269e-05 * tc[2]
            -1.49781651e-08 * tc[3]
            +4.61018008e-12 * tc[4]
            +1.28416265e+04 * invT;
        /*species 18: C2H6 */
        species[18] =
            +4.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 * invT;
        /*species 19: N2 */
        species[19] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
        /*species 20: AR */
        species[20] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +3.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 * invT;
        /*species 1: H */
        species[1] =
            +2.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 * invT;
        /*species 2: O */
        species[2] =
            +2.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 4: OH */
        species[4] =
            +3.09288767e+00
            +2.74214858e-04 * tc[1]
            +4.21684093e-08 * tc[2]
            -2.19865389e-11 * tc[3]
            +2.34824752e-15 * tc[4]
            +3.85865700e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +3.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 7: CH2 */
        species[7] =
            +2.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 * invT;
        /*species 8: CH2(S) */
        species[8] =
            +2.29203842e+00
            +2.32794318e-03 * tc[1]
            -6.70639823e-07 * tc[2]
            +1.04476500e-10 * tc[3]
            -6.79432730e-15 * tc[4]
            +5.09259997e+04 * invT;
        /*species 9: CH3 */
        species[9] =
            +2.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 * invT;
        /*species 10: CH4 */
        species[10] =
            +7.48514950e-02
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 11: CO */
        species[11] =
            +2.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +3.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 13: HCO */
        species[13] =
            +2.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 * invT;
        /*species 14: CH2O */
        species[14] =
            +1.76069008e+00
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 * invT;
        /*species 15: CH3O */
        species[15] =
            +3.77079900e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832520e+02 * invT;
        /*species 16: C2H4 */
        species[16] =
            +2.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 * invT;
        /*species 17: C2H5 */
        species[17] =
            +1.95465642e+00
            +8.69863610e-03 * tc[1]
            -2.66068889e-06 * tc[2]
            +4.38044223e-10 * tc[3]
            -2.99283152e-14 * tc[4]
            +1.28575200e+04 * invT;
        /*species 18: C2H6 */
        species[18] =
            +1.07188150e+00
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 * invT;
        /*species 19: N2 */
        species[19] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 20: AR */
        species[20] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
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
        /*species 7: CH2 */
        species[7] =
            +3.76267867e+00 * tc[0]
            +9.68872143e-04 * tc[1]
            +1.39744921e-06 * tc[2]
            -1.28363718e-09 * tc[3]
            +4.21854298e-13 * tc[4]
            +1.56253185e+00 ;
        /*species 8: CH2(S) */
        species[8] =
            +4.19860411e+00 * tc[0]
            -2.36661419e-03 * tc[1]
            +4.11648110e-06 * tc[2]
            -2.22938660e-09 * tc[3]
            +4.85786843e-13 * tc[4]
            -7.69118967e-01 ;
        /*species 9: CH3 */
        species[9] =
            +3.67359040e+00 * tc[0]
            +2.01095175e-03 * tc[1]
            +2.86510928e-06 * tc[2]
            -2.29039142e-09 * tc[3]
            +6.35964335e-13 * tc[4]
            +1.60456433e+00 ;
        /*species 10: CH4 */
        species[10] =
            +5.14987613e+00 * tc[0]
            -1.36709788e-02 * tc[1]
            +2.45900299e-05 * tc[2]
            -1.61581009e-08 * tc[3]
            +4.16734890e-12 * tc[4]
            -4.64130376e+00 ;
        /*species 11: CO */
        species[11] =
            +3.57953347e+00 * tc[0]
            -6.10353680e-04 * tc[1]
            +5.08407165e-07 * tc[2]
            +3.02335295e-10 * tc[3]
            -2.26106125e-13 * tc[4]
            +3.50840928e+00 ;
        /*species 12: CO2 */
        species[12] =
            +2.35677352e+00 * tc[0]
            +8.98459677e-03 * tc[1]
            -3.56178134e-06 * tc[2]
            +8.19730073e-10 * tc[3]
            -3.59248870e-14 * tc[4]
            +9.90105222e+00 ;
        /*species 13: HCO */
        species[13] =
            +4.22118584e+00 * tc[0]
            -3.24392532e-03 * tc[1]
            +6.88997230e-06 * tc[2]
            -4.43813643e-09 * tc[3]
            +1.08442216e-12 * tc[4]
            +3.39437243e+00 ;
        /*species 14: CH2O */
        species[14] =
            +4.79372315e+00 * tc[0]
            -9.90833369e-03 * tc[1]
            +1.86610004e-05 * tc[2]
            -1.26428420e-08 * tc[3]
            +3.29431630e-12 * tc[4]
            +6.02812900e-01 ;
        /*species 15: CH3O */
        species[15] =
            +2.10620400e+00 * tc[0]
            +7.21659500e-03 * tc[1]
            +2.66923600e-06 * tc[2]
            -2.45921200e-09 * tc[3]
            +5.18902500e-13 * tc[4]
            +1.31521770e+01 ;
        /*species 16: C2H4 */
        species[16] =
            +3.95920148e+00 * tc[0]
            -7.57052247e-03 * tc[1]
            +2.85495146e-05 * tc[2]
            -2.30529584e-08 * tc[3]
            +6.74710933e-12 * tc[4]
            +4.09733096e+00 ;
        /*species 17: C2H5 */
        species[17] =
            +4.30646568e+00 * tc[0]
            -4.18658892e-03 * tc[1]
            +2.48571403e-05 * tc[2]
            -1.99708869e-08 * tc[3]
            +5.76272510e-12 * tc[4]
            +4.70720924e+00 ;
        /*species 18: C2H6 */
        species[18] =
            +4.29142492e+00 * tc[0]
            -5.50154270e-03 * tc[1]
            +2.99719144e-05 * tc[2]
            -2.36155428e-08 * tc[3]
            +6.71714427e-12 * tc[4]
            +2.66682316e+00 ;
        /*species 19: N2 */
        species[19] =
            +3.29867700e+00 * tc[0]
            +1.40824040e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213500e-13 * tc[4]
            +3.95037200e+00 ;
        /*species 20: AR */
        species[20] =
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
        /*species 7: CH2 */
        species[7] =
            +2.87410113e+00 * tc[0]
            +3.65639292e-03 * tc[1]
            -7.04472985e-07 * tc[2]
            +8.67265163e-11 * tc[3]
            -4.69318918e-15 * tc[4]
            +6.17119324e+00 ;
        /*species 8: CH2(S) */
        species[8] =
            +2.29203842e+00 * tc[0]
            +4.65588637e-03 * tc[1]
            -1.00595973e-06 * tc[2]
            +1.39302000e-10 * tc[3]
            -8.49290912e-15 * tc[4]
            +8.62650169e+00 ;
        /*species 9: CH3 */
        species[9] =
            +2.28571772e+00 * tc[0]
            +7.23990037e-03 * tc[1]
            -1.49357174e-06 * tc[2]
            +1.98561548e-10 * tc[3]
            -1.16788599e-14 * tc[4]
            +8.48007179e+00 ;
        /*species 10: CH4 */
        species[10] =
            +7.48514950e-02 * tc[0]
            +1.33909467e-02 * tc[1]
            -2.86642905e-06 * tc[2]
            +4.07641783e-10 * tc[3]
            -2.54538075e-14 * tc[4]
            +1.84373180e+01 ;
        /*species 11: CO */
        species[11] =
            +2.71518561e+00 * tc[0]
            +2.06252743e-03 * tc[1]
            -4.99412886e-07 * tc[2]
            +7.66843360e-11 * tc[3]
            -5.09119290e-15 * tc[4]
            +7.81868772e+00 ;
        /*species 12: CO2 */
        species[12] =
            +3.85746029e+00 * tc[0]
            +4.41437026e-03 * tc[1]
            -1.10740702e-06 * tc[2]
            +1.74496729e-10 * tc[3]
            -1.18021041e-14 * tc[4]
            +2.27163806e+00 ;
        /*species 13: HCO */
        species[13] =
            +2.77217438e+00 * tc[0]
            +4.95695526e-03 * tc[1]
            -1.24222806e-06 * tc[2]
            +1.96387259e-10 * tc[3]
            -1.33377178e-14 * tc[4]
            +9.79834492e+00 ;
        /*species 14: CH2O */
        species[14] =
            +1.76069008e+00 * tc[0]
            +9.20000082e-03 * tc[1]
            -2.21129406e-06 * tc[2]
            +3.35470707e-10 * tc[3]
            -2.20963910e-14 * tc[4]
            +1.36563230e+01 ;
        /*species 15: CH3O */
        species[15] =
            +3.77079900e+00 * tc[0]
            +7.87149700e-03 * tc[1]
            -1.32819200e-06 * tc[2]
            +1.31481033e-10 * tc[3]
            -5.28154000e-15 * tc[4]
            +2.92957500e+00 ;
        /*species 16: C2H4 */
        species[16] =
            +2.03611116e+00 * tc[0]
            +1.46454151e-02 * tc[1]
            -3.35538958e-06 * tc[2]
            +4.90743077e-10 * tc[3]
            -3.14265152e-14 * tc[4]
            +1.03053693e+01 ;
        /*species 17: C2H5 */
        species[17] =
            +1.95465642e+00 * tc[0]
            +1.73972722e-02 * tc[1]
            -3.99103334e-06 * tc[2]
            +5.84058963e-10 * tc[3]
            -3.74103940e-14 * tc[4]
            +1.34624343e+01 ;
        /*species 18: C2H6 */
        species[18] =
            +1.07188150e+00 * tc[0]
            +2.16852677e-02 * tc[1]
            -5.01280335e-06 * tc[2]
            +7.38040003e-10 * tc[3]
            -4.75007225e-14 * tc[4]
            +1.51156107e+01 ;
        /*species 19: N2 */
        species[19] =
            +2.92664000e+00 * tc[0]
            +1.48797680e-03 * tc[1]
            -2.84238000e-07 * tc[2]
            +3.36567933e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
        /*species 20: AR */
        species[20] =
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
    wt[7] = 14.027090; /*CH2 */
    wt[8] = 14.027090; /*CH2(S) */
    wt[9] = 15.035060; /*CH3 */
    wt[10] = 16.043030; /*CH4 */
    wt[11] = 28.010550; /*CO */
    wt[12] = 44.009950; /*CO2 */
    wt[13] = 29.018520; /*HCO */
    wt[14] = 30.026490; /*CH2O */
    wt[15] = 31.034460; /*CH3O */
    wt[16] = 28.054180; /*C2H4 */
    wt[17] = 29.062150; /*C2H5 */
    wt[18] = 30.070120; /*C2H6 */
    wt[19] = 28.013400; /*N2 */
    wt[20] = 39.948000; /*AR */

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
    y[7] = phi[7]*14.027090;   XW += y[7]; /*CH2 */
    y[8] = phi[8]*14.027090;   XW += y[8]; /*CH2(S) */
    y[9] = phi[9]*15.035060;   XW += y[9]; /*CH3 */
    y[10] = phi[10]*16.043030;   XW += y[10]; /*CH4 */
    y[11] = phi[11]*28.010550;   XW += y[11]; /*CO */
    y[12] = phi[12]*44.009950;   XW += y[12]; /*CO2 */
    y[13] = phi[13]*29.018520;   XW += y[13]; /*HCO */
    y[14] = phi[14]*30.026490;   XW += y[14]; /*CH2O */
    y[15] = phi[15]*31.034460;   XW += y[15]; /*CH3O */
    y[16] = phi[16]*28.054180;   XW += y[16]; /*C2H4 */
    y[17] = phi[17]*29.062150;   XW += y[17]; /*C2H5 */
    y[18] = phi[18]*30.070120;   XW += y[18]; /*C2H6 */
    y[19] = phi[19]*28.013400;   XW += y[19]; /*N2 */
    y[20] = phi[20]*39.948000;   XW += y[20]; /*AR */
    for (id = 0; id < 21; ++id) {
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
    phi[7] = y[7]/ 1.40270900e-02; /*CH2 (wt in kg) */
    phi[8] = y[8]/ 1.40270900e-02; /*CH2(S) (wt in kg) */
    phi[9] = y[9]/ 1.50350600e-02; /*CH3 (wt in kg) */
    phi[10] = y[10]/ 1.60430300e-02; /*CH4 (wt in kg) */
    phi[11] = y[11]/ 2.80105500e-02; /*CO (wt in kg) */
    phi[12] = y[12]/ 4.40099500e-02; /*CO2 (wt in kg) */
    phi[13] = y[13]/ 2.90185200e-02; /*HCO (wt in kg) */
    phi[14] = y[14]/ 3.00264900e-02; /*CH2O (wt in kg) */
    phi[15] = y[15]/ 3.10344600e-02; /*CH3O (wt in kg) */
    phi[16] = y[16]/ 2.80541800e-02; /*C2H4 (wt in kg) */
    phi[17] = y[17]/ 2.90621500e-02; /*C2H5 (wt in kg) */
    phi[18] = y[18]/ 3.00701200e-02; /*C2H6 (wt in kg) */
    phi[19] = y[19]/ 2.80134000e-02; /*N2 (wt in kg) */
    phi[20] = y[20]/ 3.99480000e-02; /*AR (wt in kg) */

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
    y[7] = c[7] * 14.027090 / (*rho); 
    y[8] = c[8] * 14.027090 / (*rho); 
    y[9] = c[9] * 15.035060 / (*rho); 
    y[10] = c[10] * 16.043030 / (*rho); 
    y[11] = c[11] * 28.010550 / (*rho); 
    y[12] = c[12] * 44.009950 / (*rho); 
    y[13] = c[13] * 29.018520 / (*rho); 
    y[14] = c[14] * 30.026490 / (*rho); 
    y[15] = c[15] * 31.034460 / (*rho); 
    y[16] = c[16] * 28.054180 / (*rho); 
    y[17] = c[17] * 29.062150 / (*rho); 
    y[18] = c[18] * 30.070120 / (*rho); 
    y[19] = c[19] * 28.013400 / (*rho); 
    y[20] = c[20] * 39.948000 / (*rho); 

    return;
}


/*ddebdf compatible right hand side of CV burner */
/*rwrk[0] and rwrk[1] should contain rho and ene respectively */
/*working variable phi contains specific mole numbers */
void fecvrhs_(double * time, double * phi, double * phidot, double * rwrk, int * iwrk)
{
    double rho,ene; /*CV Parameters */
    double y[21], wdot[21]; /*temporary storage */
    int i; /*Loop counter */
    double temperature,pressure; /*temporary var */
    rho = rwrk[0];
    ene = rwrk[1];
    fephity_(phi, iwrk, rwrk, y);
    feeytt_(&ene, y, iwrk, rwrk, &temperature);
    CKPY(&rho, &temperature,  y, iwrk, rwrk, &pressure);
    CKWYP(&pressure, &temperature,  y, iwrk, rwrk, wdot);
    for (i=0; i<21; ++i) phidot[i] = wdot[i] / (rho/1000.0); 

    return;
}


/*returns the dimensionality of the cv burner (number of species) */
int fecvdim_()
{
    return 21;
}


/*ddebdf compatible right hand side of ZND solver */
/*rwrk[0] : scaling factor for pressure */
/*rwrk[1] : preshock density (g/cc)  */
/*rwrk[2] : detonation velocity (cm/s)  */
/*solution vector: [P; rho; y0 ... ylast]  */
void fezndrhs_(double * time, double * z, double * zdot, double * rwrk, int * iwrk)
{
    double psc,rho1,udet; /*ZND Parameters */
    double wt[21], hms[21], wdot[21]; /*temporary storage */
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
    for (i=0; i<21; ++i) {
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
    return 24;
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
    if (strcmp(s1, "CH2")==0) return 7; 
    if (strcmp(s1, "CH2(S)")==0) return 8; 
    if (strcmp(s1, "CH3")==0) return 9; 
    if (strcmp(s1, "CH4")==0) return 10; 
    if (strcmp(s1, "CO")==0) return 11; 
    if (strcmp(s1, "CO2")==0) return 12; 
    if (strcmp(s1, "HCO")==0) return 13; 
    if (strcmp(s1, "CH2O")==0) return 14; 
    if (strcmp(s1, "CH3O")==0) return 15; 
    if (strcmp(s1, "C2H4")==0) return 16; 
    if (strcmp(s1, "C2H5")==0) return 17; 
    if (strcmp(s1, "C2H6")==0) return 18; 
    if (strcmp(s1, "N2")==0) return 19; 
    if (strcmp(s1, "AR")==0) return 20; 
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
    if (sn==7) return "CH2"; 
    if (sn==8) return "CH2(S)"; 
    if (sn==9) return "CH3"; 
    if (sn==10) return "CH4"; 
    if (sn==11) return "CO"; 
    if (sn==12) return "CO2"; 
    if (sn==13) return "HCO"; 
    if (sn==14) return "CH2O"; 
    if (sn==15) return "CH3O"; 
    if (sn==16) return "C2H4"; 
    if (sn==17) return "C2H5"; 
    if (sn==18) return "C2H6"; 
    if (sn==19) return "N2"; 
    if (sn==20) return "AR"; 
    /*species name not found */
    return "NOTFOUND";
}

/* End of file  */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENIMC EGTRANSETLENIMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENIMC egtransetlenimc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENIMC egtransetlenimc_
#endif
extern "C" { void egtransetLENIMC(int* LENIMC) {
  *LENIMC =           86;}}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENRMC EGTRANSETLENRMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENRMC egtransetlenrmc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENRMC egtransetlenrmc_
#endif
extern "C" { void egtransetLENRMC(int* LENRMC) {
  *LENRMC =         9114;}}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNO EGTRANSETNO
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNO egtransetno
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNO egtransetno_
#endif
extern "C" { void egtransetNO(int* NO) {
  *NO =            4;}}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetKK EGTRANSETKK
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetKK egtransetkk
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetKK egtransetkk_
#endif
extern "C" { void egtransetKK(int* KK) {
  *KK =           21;}}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLITE EGTRANSETNLITE
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLITE egtransetnlite
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLITE egtransetnlite_
#endif
extern "C" { void egtransetNLITE(int* NLITE) {
  *NLITE =            2;}}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetPATM EGTRANSETPATM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetPATM egtransetpatm
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetPATM egtransetpatm_
#endif
extern "C" { void egtransetPATM(double* PATM) {
  *PATM =   0.1013250000000000E+07;}}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetWT EGTRANSETWT
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetWT egtransetwt
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetWT egtransetwt_
#endif
extern "C" { void egtransetWT(double* WT) {
  WT[           0] =   0.2015939950942993E+01;
  WT[           1] =   0.1007969975471497E+01;
  WT[           2] =   0.1599940013885498E+02;
  WT[           3] =   0.3199880027770996E+02;
  WT[           4] =   0.1700737011432648E+02;
  WT[           5] =   0.1801534008979797E+02;
  WT[           6] =   0.3300677025318146E+02;
  WT[           7] =   0.1402709031105042E+02;
  WT[           8] =   0.1402709031105042E+02;
  WT[           9] =   0.1503506028652191E+02;
  WT[          10] =   0.1604303026199341E+02;
  WT[          11] =   0.2801055049896240E+02;
  WT[          12] =   0.4400995063781738E+02;
  WT[          13] =   0.2901852047443390E+02;
  WT[          14] =   0.3002649044990540E+02;
  WT[          15] =   0.3103446042537689E+02;
  WT[          16] =   0.2805418062210083E+02;
  WT[          17] =   0.2906215059757233E+02;
  WT[          18] =   0.3007012057304382E+02;
  WT[          19] =   0.2801339912414551E+02;
  WT[          20] =   0.3994800186157227E+02;
};  }
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetEPS EGTRANSETEPS
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetEPS egtranseteps
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetEPS egtranseteps_
#endif
extern "C" { void egtransetEPS(double* EPS) {
  EPS[           0] =   0.3800000000000000E+02;
  EPS[           1] =   0.1450000000000000E+03;
  EPS[           2] =   0.8000000000000000E+02;
  EPS[           3] =   0.1074000000000000E+03;
  EPS[           4] =   0.8000000000000000E+02;
  EPS[           5] =   0.5724000000000000E+03;
  EPS[           6] =   0.1074000000000000E+03;
  EPS[           7] =   0.1440000000000000E+03;
  EPS[           8] =   0.1440000000000000E+03;
  EPS[           9] =   0.1440000000000000E+03;
  EPS[          10] =   0.1414000000000000E+03;
  EPS[          11] =   0.9809999999999999E+02;
  EPS[          12] =   0.2440000000000000E+03;
  EPS[          13] =   0.4980000000000000E+03;
  EPS[          14] =   0.4980000000000000E+03;
  EPS[          15] =   0.4170000000000000E+03;
  EPS[          16] =   0.2808000000000000E+03;
  EPS[          17] =   0.2523000000000000E+03;
  EPS[          18] =   0.2523000000000000E+03;
  EPS[          19] =   0.9753000000000000E+02;
  EPS[          20] =   0.1365000000000000E+03;
};  }
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetSIG EGTRANSETSIG
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetSIG egtransetsig
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetSIG egtransetsig_
#endif
extern "C" { void egtransetSIG(double* SIG) {
  SIG[           0] =   0.2920000000000000E+01;
  SIG[           1] =   0.2050000000000000E+01;
  SIG[           2] =   0.2750000000000000E+01;
  SIG[           3] =   0.3458000000000000E+01;
  SIG[           4] =   0.2750000000000000E+01;
  SIG[           5] =   0.2605000000000000E+01;
  SIG[           6] =   0.3458000000000000E+01;
  SIG[           7] =   0.3800000000000000E+01;
  SIG[           8] =   0.3800000000000000E+01;
  SIG[           9] =   0.3800000000000000E+01;
  SIG[          10] =   0.3746000000000000E+01;
  SIG[          11] =   0.3650000000000000E+01;
  SIG[          12] =   0.3763000000000000E+01;
  SIG[          13] =   0.3590000000000000E+01;
  SIG[          14] =   0.3590000000000000E+01;
  SIG[          15] =   0.3690000000000000E+01;
  SIG[          16] =   0.3971000000000000E+01;
  SIG[          17] =   0.4302000000000000E+01;
  SIG[          18] =   0.4302000000000000E+01;
  SIG[          19] =   0.3621000000000000E+01;
  SIG[          20] =   0.3330000000000000E+01;
};  }
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetDIP EGTRANSETDIP
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetDIP egtransetdip
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetDIP egtransetdip_
#endif
extern "C" { void egtransetDIP(double* DIP) {
  DIP[           0] =   0.0000000000000000E+00;
  DIP[           1] =   0.0000000000000000E+00;
  DIP[           2] =   0.0000000000000000E+00;
  DIP[           3] =   0.0000000000000000E+00;
  DIP[           4] =   0.0000000000000000E+00;
  DIP[           5] =   0.1844000000000000E+01;
  DIP[           6] =   0.0000000000000000E+00;
  DIP[           7] =   0.0000000000000000E+00;
  DIP[           8] =   0.0000000000000000E+00;
  DIP[           9] =   0.0000000000000000E+00;
  DIP[          10] =   0.0000000000000000E+00;
  DIP[          11] =   0.0000000000000000E+00;
  DIP[          12] =   0.0000000000000000E+00;
  DIP[          13] =   0.0000000000000000E+00;
  DIP[          14] =   0.0000000000000000E+00;
  DIP[          15] =   0.1700000000000000E+01;
  DIP[          16] =   0.0000000000000000E+00;
  DIP[          17] =   0.0000000000000000E+00;
  DIP[          18] =   0.0000000000000000E+00;
  DIP[          19] =   0.0000000000000000E+00;
  DIP[          20] =   0.0000000000000000E+00;
};  }
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetPOL EGTRANSETPOL
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetPOL egtransetpol
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetPOL egtransetpol_
#endif
extern "C" { void egtransetPOL(double* POL) {
  POL[           0] =   0.7900000000000000E+00;
  POL[           1] =   0.0000000000000000E+00;
  POL[           2] =   0.0000000000000000E+00;
  POL[           3] =   0.1600000000000000E+01;
  POL[           4] =   0.0000000000000000E+00;
  POL[           5] =   0.0000000000000000E+00;
  POL[           6] =   0.0000000000000000E+00;
  POL[           7] =   0.0000000000000000E+00;
  POL[           8] =   0.0000000000000000E+00;
  POL[           9] =   0.0000000000000000E+00;
  POL[          10] =   0.2600000000000000E+01;
  POL[          11] =   0.1950000000000000E+01;
  POL[          12] =   0.2650000000000000E+01;
  POL[          13] =   0.0000000000000000E+00;
  POL[          14] =   0.0000000000000000E+00;
  POL[          15] =   0.0000000000000000E+00;
  POL[          16] =   0.0000000000000000E+00;
  POL[          17] =   0.0000000000000000E+00;
  POL[          18] =   0.0000000000000000E+00;
  POL[          19] =   0.1760000000000000E+01;
  POL[          20] =   0.0000000000000000E+00;
};  }
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetZROT EGTRANSETZROT
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetZROT egtransetzrot
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetZROT egtransetzrot_
#endif
extern "C" { void egtransetZROT(double* ZROT) {
  ZROT[           0] =   0.2800000000000000E+03;
  ZROT[           1] =   0.0000000000000000E+00;
  ZROT[           2] =   0.0000000000000000E+00;
  ZROT[           3] =   0.3800000000000000E+01;
  ZROT[           4] =   0.0000000000000000E+00;
  ZROT[           5] =   0.4000000000000000E+01;
  ZROT[           6] =   0.1000000000000000E+01;
  ZROT[           7] =   0.0000000000000000E+00;
  ZROT[           8] =   0.0000000000000000E+00;
  ZROT[           9] =   0.0000000000000000E+00;
  ZROT[          10] =   0.1300000000000000E+02;
  ZROT[          11] =   0.1800000000000000E+01;
  ZROT[          12] =   0.2100000000000000E+01;
  ZROT[          13] =   0.0000000000000000E+00;
  ZROT[          14] =   0.2000000000000000E+01;
  ZROT[          15] =   0.2000000000000000E+01;
  ZROT[          16] =   0.1500000000000000E+01;
  ZROT[          17] =   0.1500000000000000E+01;
  ZROT[          18] =   0.1500000000000000E+01;
  ZROT[          19] =   0.4000000000000000E+01;
  ZROT[          20] =   0.0000000000000000E+00;
};  }
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLIN EGTRANSETNLIN
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLIN egtransetnlin
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLIN egtransetnlin_
#endif
extern "C" { void egtransetNLIN(int* NLIN) {
  NLIN[           0] =            1;
  NLIN[           1] =            0;
  NLIN[           2] =            0;
  NLIN[           3] =            1;
  NLIN[           4] =            1;
  NLIN[           5] =            2;
  NLIN[           6] =            2;
  NLIN[           7] =            1;
  NLIN[           8] =            1;
  NLIN[           9] =            1;
  NLIN[          10] =            2;
  NLIN[          11] =            1;
  NLIN[          12] =            1;
  NLIN[          13] =            2;
  NLIN[          14] =            2;
  NLIN[          15] =            2;
  NLIN[          16] =            2;
  NLIN[          17] =            2;
  NLIN[          18] =            2;
  NLIN[          19] =            1;
  NLIN[          20] =            0;
};  }
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFLAM EGTRANSETCOFLAM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFLAM egtransetcoflam
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFLAM egtransetcoflam_
#endif
extern "C" { void egtransetCOFLAM(double* COFLAM) {
  COFLAM[           0] =   0.9235828509713865E+01;
  COFLAM[           1] =  -0.4674267763951823E+00;
  COFLAM[           2] =   0.1156734789350108E+00;
  COFLAM[           3] =  -0.2596025563286716E-02;
  COFLAM[           4] =  -0.8554281253040258E+00;
  COFLAM[           5] =   0.3652691486294877E+01;
  COFLAM[           6] =  -0.3980303020607808E+00;
  COFLAM[           7] =   0.1757072886401263E-01;
  COFLAM[           8] =   0.1685025383414630E+01;
  COFLAM[           9] =   0.1929024677614472E+01;
  COFLAM[          10] =  -0.1738657445088021E+00;
  COFLAM[          11] =   0.7841476915125366E-02;
  COFLAM[          12] =  -0.1936698464483141E+01;
  COFLAM[          13] =   0.2890477542404777E+01;
  COFLAM[          14] =  -0.2709591162375214E+00;
  COFLAM[          15] =   0.1152570281409322E-01;
  COFLAM[          16] =   0.1416223092993048E+02;
  COFLAM[          17] =  -0.3244626711115035E+01;
  COFLAM[          18] =   0.5336588172885086E+00;
  COFLAM[          19] =  -0.2328116832308653E-01;
  COFLAM[          20] =   0.2336546540925643E+02;
  COFLAM[          21] =  -0.8965822807102029E+01;
  COFLAM[          22] =   0.1528828068043276E+01;
  COFLAM[          23] =  -0.7590175979167682E-01;
  COFLAM[          24] =  -0.1130096296226325E+01;
  COFLAM[          25] =   0.2340066562944334E+01;
  COFLAM[          26] =  -0.1632055932993225E+00;
  COFLAM[          27] =   0.5799980518407136E-02;
  COFLAM[          28] =   0.1291635112732534E+02;
  COFLAM[          29] =  -0.3737065816673101E+01;
  COFLAM[          30] =   0.7157937477351572E+00;
  COFLAM[          31] =  -0.3638374043162176E-01;
  COFLAM[          32] =   0.1893642052511076E+02;
  COFLAM[          33] =  -0.6509584574818159E+01;
  COFLAM[          34] =   0.1132853049549838E+01;
  COFLAM[          35] =  -0.5695791215151112E-01;
  COFLAM[          36] =   0.1399177056183472E+02;
  COFLAM[          37] =  -0.4641882989478738E+01;
  COFLAM[          38] =   0.9076435786312527E+00;
  COFLAM[          39] =  -0.4772397826828819E-01;
  COFLAM[          40] =   0.1330618431273573E+02;
  COFLAM[          41] =  -0.4960294456747341E+01;
  COFLAM[          42] =   0.1032808842591664E+01;
  COFLAM[          43] =  -0.5633567903356620E-01;
  COFLAM[          44] =   0.1187710072654253E+02;
  COFLAM[          45] =  -0.3154801252509626E+01;
  COFLAM[          46] =   0.6020483454882906E+00;
  COFLAM[          47] =  -0.3032714732778013E-01;
  COFLAM[          48] =  -0.1135070822782991E+02;
  COFLAM[          49] =   0.5875667873528579E+01;
  COFLAM[          50] =  -0.5677982250304914E+00;
  COFLAM[          51] =   0.2031670238952729E-01;
  COFLAM[          52] =   0.6296028277200143E+01;
  COFLAM[          53] =  -0.2225281214264024E+01;
  COFLAM[          54] =   0.6369149369376514E+00;
  COFLAM[          55] =  -0.3808451133860326E-01;
  COFLAM[          56] =   0.5384151635429418E+01;
  COFLAM[          57] =  -0.2389146833355944E+01;
  COFLAM[          58] =   0.7389861798813848E+00;
  COFLAM[          59] =  -0.4581386404459022E-01;
  COFLAM[          60] =  -0.6138056391436082E+01;
  COFLAM[          61] =   0.2471262972505274E+01;
  COFLAM[          62] =   0.6476625884336293E-01;
  COFLAM[          63] =  -0.1455104151494988E-01;
  COFLAM[          64] =  -0.1460870939008838E+02;
  COFLAM[          65] =   0.6359722173433303E+01;
  COFLAM[          66] =  -0.5034536870651064E+00;
  COFLAM[          67] =   0.1259513381293373E-01;
  COFLAM[          68] =  -0.8941827721234478E+01;
  COFLAM[          69] =   0.4021583693700529E+01;
  COFLAM[          70] =  -0.1835686948063458E+00;
  COFLAM[          71] =  -0.1963330023789806E-02;
  COFLAM[          72] =  -0.1098227802311090E+02;
  COFLAM[          73] =   0.4703047406066884E+01;
  COFLAM[          74] =  -0.2517963879729279E+00;
  COFLAM[          75] =   0.1532853967214294E-03;
  COFLAM[          76] =   0.1293004274651541E+02;
  COFLAM[          77] =  -0.3528374680486067E+01;
  COFLAM[          78] =   0.6455829015382131E+00;
  COFLAM[          79] =  -0.3194413600157287E-01;
  COFLAM[          80] =  -0.3166636127434923E+01;
  COFLAM[          81] =   0.3467381630049378E+01;
  COFLAM[          82] =  -0.3746257298439007E+00;
  COFLAM[          83] =   0.1658331946651994E-01;
};  }
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFETA EGTRANSETCOFETA
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFETA egtransetcofeta
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFETA egtransetcofeta_
#endif
extern "C" { void egtransetCOFETA(double* COFETA) {
  COFETA[           0] =  -0.1384035451131115E+02;
  COFETA[           1] =   0.1003491325708292E+01;
  COFETA[           2] =  -0.5016044554561987E-01;
  COFETA[           3] =   0.2330995224065489E-02;
  COFETA[           4] =  -0.2040534341454124E+02;
  COFETA[           5] =   0.3652691486295009E+01;
  COFETA[           6] =  -0.3980303020608006E+00;
  COFETA[           7] =   0.1757072886401361E-01;
  COFETA[           8] =  -0.1510027705857393E+02;
  COFETA[           9] =   0.1929024677614498E+01;
  COFETA[          10] =  -0.1738657445088048E+00;
  COFETA[          11] =   0.7841476915125451E-02;
  COFETA[          12] =  -0.1715809053469514E+02;
  COFETA[          13] =   0.2678088349030608E+01;
  COFETA[          14] =  -0.2721592407921913E+00;
  COFETA[          15] =   0.1214173232611473E-01;
  COFETA[          16] =  -0.1506972928055982E+02;
  COFETA[          17] =   0.1929024677614455E+01;
  COFETA[          18] =  -0.1738657445087985E+00;
  COFETA[          19] =   0.7841476915125141E-02;
  COFETA[          20] =  -0.1055754979620539E+02;
  COFETA[          21] =  -0.1377850378695579E+01;
  COFETA[          22] =   0.4213981638352425E+00;
  COFETA[          23] =  -0.2414423055966423E-01;
  COFETA[          24] =  -0.1714258339027710E+02;
  COFETA[          25] =   0.2678088349030630E+01;
  COFETA[          26] =  -0.2721592407921950E+00;
  COFETA[          27] =   0.1214173232611492E-01;
  COFETA[          28] =  -0.2026303235679370E+02;
  COFETA[          29] =   0.3630408718312100E+01;
  COFETA[          30] =  -0.3952256976799819E+00;
  COFETA[          31] =   0.1745288966941290E-01;
  COFETA[          32] =  -0.2026303235679370E+02;
  COFETA[          33] =   0.3630408718312100E+01;
  COFETA[          34] =  -0.3952256976799819E+00;
  COFETA[          35] =   0.1745288966941290E-01;
  COFETA[          36] =  -0.2022833518474935E+02;
  COFETA[          37] =   0.3630408718312077E+01;
  COFETA[          38] =  -0.3952256976799788E+00;
  COFETA[          39] =   0.1745288966941275E-01;
  COFETA[          40] =  -0.2000457400541702E+02;
  COFETA[          41] =   0.3569542093003574E+01;
  COFETA[          42] =  -0.3874920392785122E+00;
  COFETA[          43] =   0.1712461411247914E-01;
  COFETA[          44] =  -0.1661561262594180E+02;
  COFETA[          45] =   0.2400975158113569E+01;
  COFETA[          46] =  -0.2357717790312281E+00;
  COFETA[          47] =   0.1054820948438182E-01;
  COFETA[          48] =  -0.2397057295682966E+02;
  COFETA[          49] =   0.5130426196036950E+01;
  COFETA[          50] =  -0.5724284704186094E+00;
  COFETA[          51] =   0.2440888721969576E-01;
  COFETA[          52] =  -0.1987003264957855E+02;
  COFETA[          53] =   0.2703811616927088E+01;
  COFETA[          54] =  -0.1672355189642018E+00;
  COFETA[          55] =   0.3212257118450056E-02;
  COFETA[          56] =  -0.1985295977766010E+02;
  COFETA[          57] =   0.2703811616927011E+01;
  COFETA[          58] =  -0.1672355189641908E+00;
  COFETA[          59] =   0.3212257118449540E-02;
  COFETA[          60] =  -0.1997384799268116E+02;
  COFETA[          61] =   0.2861750522819142E+01;
  COFETA[          62] =  -0.2024899143538781E+00;
  COFETA[          63] =   0.5362464436156140E-02;
  COFETA[          64] =  -0.2504838169296130E+02;
  COFETA[          65] =   0.5332796564390715E+01;
  COFETA[          66] =  -0.5890465746677508E+00;
  COFETA[          67] =   0.2473730756384840E-01;
  COFETA[          68] =  -0.2463014790579715E+02;
  COFETA[          69] =   0.5183142307718661E+01;
  COFETA[          70] =  -0.5771878965089483E+00;
  COFETA[          71] =   0.2453023678847775E-01;
  COFETA[          72] =  -0.2461310023284140E+02;
  COFETA[          73] =   0.5183142307718687E+01;
  COFETA[          74] =  -0.5771878965089522E+00;
  COFETA[          75] =   0.2453023678847794E-01;
  COFETA[          76] =  -0.1656563666406150E+02;
  COFETA[          77] =   0.2388167035581858E+01;
  COFETA[          78] =  -0.2341208182867099E+00;
  COFETA[          79] =   0.1047727172770704E-01;
  COFETA[          80] =  -0.1903691114465788E+02;
  COFETA[          81] =   0.3467381630049308E+01;
  COFETA[          82] =  -0.3746257298438905E+00;
  COFETA[          83] =   0.1658331946651944E-01;
};  }
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFD EGTRANSETCOFD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFD egtransetcofd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFD egtransetcofd_
#endif
extern "C" { void egtransetCOFD(double* COFD) {
  COFD[           0] =  -0.1031825736753964E+02;
  COFD[           1] =   0.2192412640005562E+01;
  COFD[           2] =  -0.7539036147594932E-01;
  COFD[           3] =   0.3511225190832653E-02;
  COFD[           4] =  -0.1201910564637032E+02;
  COFD[           5] =   0.3033139269467014E+01;
  COFD[           6] =  -0.1860150369684279E+00;
  COFD[           7] =   0.8361628520974778E-02;
  COFD[           8] =  -0.1093372460559999E+02;
  COFD[           9] =   0.2306133612968357E+01;
  COFD[          10] =  -0.8736427899565385E-01;
  COFD[          11] =   0.3897526119633928E-02;
  COFD[          12] =  -0.1270707915711888E+02;
  COFD[          13] =   0.2929619774221139E+01;
  COFD[          14] =  -0.1739670601708100E+00;
  COFD[          15] =   0.7901198812704345E-02;
  COFD[          16] =  -0.1093673058791694E+02;
  COFD[          17] =   0.2306118133758978E+01;
  COFD[          18] =  -0.8736238957367691E-01;
  COFD[          19] =   0.3897448280154413E-02;
  COFD[          20] =  -0.1818497883283899E+02;
  COFD[          21] =   0.5043624597310085E+01;
  COFD[          22] =  -0.4377066809210735E+00;
  COFD[          23] =   0.1887643866477002E-01;
  COFD[          24] =  -0.1271124360389947E+02;
  COFD[          25] =   0.2931047949505553E+01;
  COFD[          26] =  -0.1741679825526774E+00;
  COFD[          27] =   0.7910585398301936E-02;
  COFD[          28] =  -0.1355980676188247E+02;
  COFD[          29] =   0.3231182295054575E+01;
  COFD[          30] =  -0.2138406992439825E+00;
  COFD[          31] =   0.9659145053095021E-02;
  COFD[          32] =  -0.1355980676188247E+02;
  COFD[          33] =   0.3231182295054575E+01;
  COFD[          34] =  -0.2138406992439825E+00;
  COFD[          35] =   0.9659145053095021E-02;
  COFD[          36] =  -0.1358415095584097E+02;
  COFD[          37] =   0.3239845970896248E+01;
  COFD[          38] =  -0.2150500406377150E+00;
  COFD[          39] =   0.9715225271588263E-02;
  COFD[          40] =  -0.1353824314289343E+02;
  COFD[          41] =   0.3227352939332374E+01;
  COFD[          42] =  -0.2134245159235849E+00;
  COFD[          43] =   0.9644648995332943E-02;
  COFD[          44] =  -0.1268296956267843E+02;
  COFD[          45] =   0.2903428401599904E+01;
  COFD[          46] =  -0.1707667302851327E+00;
  COFD[          47] =   0.7771111802295727E-02;
  COFD[          48] =  -0.1564739122681956E+02;
  COFD[          49] =   0.4025825276486324E+01;
  COFD[          50] =  -0.3181195147607377E+00;
  COFD[          51] =   0.1422136468533644E-01;
  COFD[          52] =  -0.1746999275651299E+02;
  COFD[          53] =   0.4670556529920976E+01;
  COFD[          54] =  -0.3931788949833612E+00;
  COFD[          55] =   0.1709979614643272E-01;
  COFD[          56] =  -0.1747771970561251E+02;
  COFD[          57] =   0.4673159679696371E+01;
  COFD[          58] =  -0.3935109715079382E+00;
  COFD[          59] =   0.1711386921544412E-01;
  COFD[          60] =  -0.1797854335218850E+02;
  COFD[          61] =   0.4905020138101169E+01;
  COFD[          62] =  -0.4293084201751438E+00;
  COFD[          63] =   0.1891302999394448E-01;
  COFD[          64] =  -0.1646041046594534E+02;
  COFD[          65] =   0.4329097783536894E+01;
  COFD[          66] =  -0.3588538378541796E+00;
  COFD[          67] =   0.1604603061706265E-01;
  COFD[          68] =  -0.1615267929807825E+02;
  COFD[          69] =   0.4174121490551778E+01;
  COFD[          70] =  -0.3383475624805202E+00;
  COFD[          71] =   0.1513993577326728E-01;
  COFD[          72] =  -0.1616515483252968E+02;
  COFD[          73] =   0.4178910660883203E+01;
  COFD[          74] =  -0.3390062877759687E+00;
  COFD[          75] =   0.1517006334795817E-01;
  COFD[          76] =  -0.1266066516595179E+02;
  COFD[          77] =   0.2898076057146818E+01;
  COFD[          78] =  -0.1700453759318706E+00;
  COFD[          79] =   0.7738690296212197E-02;
  COFD[          80] =  -0.1342264266757377E+02;
  COFD[          81] =   0.3221777556990338E+01;
  COFD[          82] =  -0.2128673394230190E+00;
  COFD[          83] =   0.9627787551744238E-02;
  COFD[          84] =  -0.1201910564637032E+02;
  COFD[          85] =   0.3033139269467014E+01;
  COFD[          86] =  -0.1860150369684279E+00;
  COFD[          87] =   0.8361628520974778E-02;
  COFD[          88] =  -0.1519084593013382E+02;
  COFD[          89] =   0.4384306309333000E+01;
  COFD[          90] =  -0.3557550068895214E+00;
  COFD[          91] =   0.1548014912218549E-01;
  COFD[          92] =  -0.1531287529990816E+02;
  COFD[          93] =   0.4270888061293578E+01;
  COFD[          94] =  -0.3482328696611331E+00;
  COFD[          95] =   0.1545097630860032E-01;
  COFD[          96] =  -0.1761107854705398E+02;
  COFD[          97] =   0.5048254268555682E+01;
  COFD[          98] =  -0.4490719655111177E+00;
  COFD[          99] =   0.1981332179135272E-01;
  COFD[         100] =  -0.1533072430565730E+02;
  COFD[         101] =   0.4277567720643019E+01;
  COFD[         102] =  -0.3491407301790579E+00;
  COFD[         103] =   0.1549202933974097E-01;
  COFD[         104] =  -0.1651315089501800E+02;
  COFD[         105] =   0.4220701654162998E+01;
  COFD[         106] =  -0.2825947916575933E+00;
  COFD[         107] =   0.1018065117940972E-01;
  COFD[         108] =  -0.1762111333528786E+02;
  COFD[         109] =   0.5052180357519856E+01;
  COFD[         110] =  -0.4496001366665920E+00;
  COFD[         111] =   0.1983696581079775E-01;
  COFD[         112] =  -0.1798724510881690E+02;
  COFD[         113] =   0.5079566464178115E+01;
  COFD[         114] =  -0.4444140659462563E+00;
  COFD[         115] =   0.1923462920552580E-01;
  COFD[         116] =  -0.1798724510881690E+02;
  COFD[         117] =   0.5079566464178115E+01;
  COFD[         118] =  -0.4444140659462563E+00;
  COFD[         119] =   0.1923462920552580E-01;
  COFD[         120] =  -0.1801636950310948E+02;
  COFD[         121] =   0.5090088857460223E+01;
  COFD[         122] =  -0.4457552498826564E+00;
  COFD[         123] =   0.1929141580461785E-01;
  COFD[         124] =  -0.1790236652645725E+02;
  COFD[         125] =   0.5054241512059871E+01;
  COFD[         126] =  -0.4413163862361035E+00;
  COFD[         127] =   0.1910841720545356E-01;
  COFD[         128] =  -0.1719132602269712E+02;
  COFD[         129] =   0.4859367175853364E+01;
  COFD[         130] =  -0.4240580935572132E+00;
  COFD[         131] =   0.1870696403708276E-01;
  COFD[         132] =  -0.1725034548774479E+02;
  COFD[         133] =   0.4566921339972978E+01;
  COFD[         134] =  -0.3509192271095338E+00;
  COFD[         135] =   0.1402685086906688E-01;
  COFD[         136] =  -0.1540645026491837E+02;
  COFD[         137] =   0.3572427121876953E+01;
  COFD[         138] =  -0.1878763535992870E+00;
  COFD[         139] =   0.5591984032298197E-02;
  COFD[         140] =  -0.1539084008249366E+02;
  COFD[         141] =   0.3565006575777335E+01;
  COFD[         142] =  -0.1867580729432450E+00;
  COFD[         143] =   0.5536884066371989E-02;
  COFD[         144] =  -0.1561125740910447E+02;
  COFD[         145] =   0.3698144089263720E+01;
  COFD[         146] =  -0.2103080207621357E+00;
  COFD[         147] =   0.6792204226221277E-02;
  COFD[         148] =  -0.1694497223361135E+02;
  COFD[         149] =   0.4353290721870846E+01;
  COFD[         150] =  -0.3147441391485086E+00;
  COFD[         151] =   0.1210201943729390E-01;
  COFD[         152] =  -0.1738386031480564E+02;
  COFD[         153] =   0.4515894896789344E+01;
  COFD[         154] =  -0.3391856065220307E+00;
  COFD[         155] =   0.1329276720561771E-01;
  COFD[         156] =  -0.1738223879329547E+02;
  COFD[         157] =   0.4514299549461490E+01;
  COFD[         158] =  -0.3388651702152571E+00;
  COFD[         159] =   0.1327383213284630E-01;
  COFD[         160] =  -0.1711567017943538E+02;
  COFD[         161] =   0.4833264586465570E+01;
  COFD[         162] =  -0.4205908440542898E+00;
  COFD[         163] =   0.1855345683322859E-01;
  COFD[         164] =  -0.1752038279713604E+02;
  COFD[         165] =   0.4959639086677385E+01;
  COFD[         166] =  -0.4293974875408633E+00;
  COFD[         167] =   0.1860806087249651E-01;
  COFD[         168] =  -0.1093372460559999E+02;
  COFD[         169] =   0.2306133612968357E+01;
  COFD[         170] =  -0.8736427899565385E-01;
  COFD[         171] =   0.3897526119633928E-02;
  COFD[         172] =  -0.1531287529990816E+02;
  COFD[         173] =   0.4270888061293578E+01;
  COFD[         174] =  -0.3482328696611331E+00;
  COFD[         175] =   0.1545097630860032E-01;
  COFD[         176] =  -0.1356847565002614E+02;
  COFD[         177] =   0.3062175626208075E+01;
  COFD[         178] =  -0.1889691489381823E+00;
  COFD[         179] =   0.8454139853458284E-02;
  COFD[         180] =  -0.1511165027002312E+02;
  COFD[         181] =   0.3515441252015493E+01;
  COFD[         182] =  -0.2491198761221894E+00;
  COFD[         183] =   0.1111625360314578E-01;
  COFD[         184] =  -0.1358466752355666E+02;
  COFD[         185] =   0.3062671492876047E+01;
  COFD[         186] =  -0.1890385291581307E+00;
  COFD[         187] =   0.8457362732499073E-02;
  COFD[         188] =  -0.1906508600184398E+02;
  COFD[         189] =   0.4989055420187643E+01;
  COFD[         190] =  -0.4192832417141421E+00;
  COFD[         191] =   0.1761870656062123E-01;
  COFD[         192] =  -0.1513524421958843E+02;
  COFD[         193] =   0.3523297948981071E+01;
  COFD[         194] =  -0.2502118585975575E+00;
  COFD[         195] =   0.1116667641067510E-01;
  COFD[         196] =  -0.1558834457822415E+02;
  COFD[         197] =   0.3699567664433178E+01;
  COFD[         198] =  -0.2706504848826866E+00;
  COFD[         199] =   0.1194595865045622E-01;
  COFD[         200] =  -0.1558834457822415E+02;
  COFD[         201] =   0.3699567664433178E+01;
  COFD[         202] =  -0.2706504848826866E+00;
  COFD[         203] =   0.1194595865045622E-01;
  COFD[         204] =  -0.1561687122378518E+02;
  COFD[         205] =   0.3703886367486989E+01;
  COFD[         206] =  -0.2712410029515650E+00;
  COFD[         207] =   0.1197280722366742E-01;
  COFD[         208] =  -0.1556217965348612E+02;
  COFD[         209] =   0.3683628994574988E+01;
  COFD[         210] =  -0.2686556693855442E+00;
  COFD[         211] =   0.1186253607235680E-01;
  COFD[         212] =  -0.1478198899238034E+02;
  COFD[         213] =   0.3371092581023257E+01;
  COFD[         214] =  -0.2298567856863946E+00;
  COFD[         215] =   0.1025846107153025E-01;
  COFD[         216] =  -0.1793673861488134E+02;
  COFD[         217] =   0.4489875914439151E+01;
  COFD[         218] =  -0.3697556337020395E+00;
  COFD[         219] =   0.1609253879718407E-01;
  COFD[         220] =  -0.1915075646861156E+02;
  COFD[         221] =   0.4868838487330766E+01;
  COFD[         222] =  -0.4046210142562804E+00;
  COFD[         223] =   0.1701122981099066E-01;
  COFD[         224] =  -0.1914986770399175E+02;
  COFD[         225] =   0.4864899621234062E+01;
  COFD[         226] =  -0.4039220297382060E+00;
  COFD[         227] =   0.1697258735720344E-01;
  COFD[         228] =  -0.1878355293596396E+02;
  COFD[         229] =   0.4745970603113269E+01;
  COFD[         230] =  -0.3924372163577848E+00;
  COFD[         231] =   0.1663659729501466E-01;
  COFD[         232] =  -0.1800091724620567E+02;
  COFD[         233] =   0.4493072014688184E+01;
  COFD[         234] =  -0.3678451189139861E+00;
  COFD[         235] =   0.1591448159607913E-01;
  COFD[         236] =  -0.1802433609839966E+02;
  COFD[         237] =   0.4483646303442399E+01;
  COFD[         238] =  -0.3688151334706006E+00;
  COFD[         239] =   0.1604826293273014E-01;
  COFD[         240] =  -0.1805221831789452E+02;
  COFD[         241] =   0.4492279525932131E+01;
  COFD[         242] =  -0.3699242263001061E+00;
  COFD[         243] =   0.1609559465701492E-01;
  COFD[         244] =  -0.1474819919893963E+02;
  COFD[         245] =   0.3361311502126538E+01;
  COFD[         246] =  -0.2285465363406641E+00;
  COFD[         247] =   0.1019990809984005E-01;
  COFD[         248] =  -0.1581858256688188E+02;
  COFD[         249] =   0.3775899039344615E+01;
  COFD[         250] =  -0.2814687920927101E+00;
  COFD[         251] =   0.1245218573218981E-01;
  COFD[         252] =  -0.1270707915711888E+02;
  COFD[         253] =   0.2929619774221139E+01;
  COFD[         254] =  -0.1739670601708100E+00;
  COFD[         255] =   0.7901198812704345E-02;
  COFD[         256] =  -0.1761107854705398E+02;
  COFD[         257] =   0.5048254268555682E+01;
  COFD[         258] =  -0.4490719655111177E+00;
  COFD[         259] =   0.1981332179135272E-01;
  COFD[         260] =  -0.1511165027002312E+02;
  COFD[         261] =   0.3515441252015493E+01;
  COFD[         262] =  -0.2491198761221894E+00;
  COFD[         263] =   0.1111625360314578E-01;
  COFD[         264] =  -0.1605076758771229E+02;
  COFD[         265] =   0.3688226344868140E+01;
  COFD[         266] =  -0.2690947493495641E+00;
  COFD[         267] =   0.1187504521154671E-01;
  COFD[         268] =  -0.1509620843583920E+02;
  COFD[         269] =   0.3500304707531062E+01;
  COFD[         270] =  -0.2470160084537919E+00;
  COFD[         271] =   0.1101910392929726E-01;
  COFD[         272] =  -0.2029800570223752E+02;
  COFD[         273] =   0.5167991782317953E+01;
  COFD[         274] =  -0.4260628585538284E+00;
  COFD[         275] =   0.1725358239585169E-01;
  COFD[         276] =  -0.1605887758547131E+02;
  COFD[         277] =   0.3688399632843664E+01;
  COFD[         278] =  -0.2691183502617316E+00;
  COFD[         279] =   0.1187611401222301E-01;
  COFD[         280] =  -0.1707915522153321E+02;
  COFD[         281] =   0.4132970234506635E+01;
  COFD[         282] =  -0.3260547355537204E+00;
  COFD[         283] =   0.1431106723651207E-01;
  COFD[         284] =  -0.1707915522153321E+02;
  COFD[         285] =   0.4132970234506635E+01;
  COFD[         286] =  -0.3260547355537204E+00;
  COFD[         287] =   0.1431106723651207E-01;
  COFD[         288] =  -0.1707094794267795E+02;
  COFD[         289] =   0.4119810624341889E+01;
  COFD[         290] =  -0.3242844616145015E+00;
  COFD[         291] =   0.1423184749658152E-01;
  COFD[         292] =  -0.1696131964138645E+02;
  COFD[         293] =   0.4074167708212651E+01;
  COFD[         294] =  -0.3182576642303228E+00;
  COFD[         295] =   0.1396612551510442E-01;
  COFD[         296] =  -0.1582244348422915E+02;
  COFD[         297] =   0.3599211636009877E+01;
  COFD[         298] =  -0.2581898369286856E+00;
  COFD[         299] =   0.1143096575439008E-01;
  COFD[         300] =  -0.1868514242925209E+02;
  COFD[         301] =   0.4565512377700616E+01;
  COFD[         302] =  -0.3754538805580755E+00;
  COFD[         303] =   0.1617305773699574E-01;
  COFD[         304] =  -0.2031657960333580E+02;
  COFD[         305] =   0.5113614640951361E+01;
  COFD[         306] =  -0.4319293017089550E+00;
  COFD[         307] =   0.1802805482147992E-01;
  COFD[         308] =  -0.2032689236579443E+02;
  COFD[         309] =   0.5114221260960875E+01;
  COFD[         310] =  -0.4320136109530544E+00;
  COFD[         311] =   0.1803190494311217E-01;
  COFD[         312] =  -0.2016013281728147E+02;
  COFD[         313] =   0.5063566265537697E+01;
  COFD[         314] =  -0.4281515792120217E+00;
  COFD[         315] =   0.1797180686456706E-01;
  COFD[         316] =  -0.1901929468754062E+02;
  COFD[         317] =   0.4693839754776671E+01;
  COFD[         318] =  -0.3901000075815114E+00;
  COFD[         319] =   0.1672631492507185E-01;
  COFD[         320] =  -0.1881741117557303E+02;
  COFD[         321] =   0.4595519619722809E+01;
  COFD[         322] =  -0.3791053358880415E+00;
  COFD[         323] =   0.1632232422675412E-01;
  COFD[         324] =  -0.1882824574536174E+02;
  COFD[         325] =   0.4596254296231404E+01;
  COFD[         326] =  -0.3791928797888692E+00;
  COFD[         327] =   0.1632573982847240E-01;
  COFD[         328] =  -0.1578873902175925E+02;
  COFD[         329] =   0.3589258101555060E+01;
  COFD[         330] =  -0.2568718883472193E+00;
  COFD[         331] =   0.1137271971440306E-01;
  COFD[         332] =  -0.1677455230511400E+02;
  COFD[         333] =   0.3946214344521416E+01;
  COFD[         334] =  -0.3012613112287137E+00;
  COFD[         335] =   0.1321325422150532E-01;
  COFD[         336] =  -0.1093673058791694E+02;
  COFD[         337] =   0.2306118133758978E+01;
  COFD[         338] =  -0.8736238957367691E-01;
  COFD[         339] =   0.3897448280154413E-02;
  COFD[         340] =  -0.1533072430565730E+02;
  COFD[         341] =   0.4277567720643019E+01;
  COFD[         342] =  -0.3491407301790579E+00;
  COFD[         343] =   0.1549202933974097E-01;
  COFD[         344] =  -0.1358466752355666E+02;
  COFD[         345] =   0.3062671492876047E+01;
  COFD[         346] =  -0.1890385291581307E+00;
  COFD[         347] =   0.8457362732499073E-02;
  COFD[         348] =  -0.1509620843583920E+02;
  COFD[         349] =   0.3500304707531062E+01;
  COFD[         350] =  -0.2470160084537919E+00;
  COFD[         351] =   0.1101910392929726E-01;
  COFD[         352] =  -0.1359902342804023E+02;
  COFD[         353] =   0.3062175626208111E+01;
  COFD[         354] =  -0.1889691489381874E+00;
  COFD[         355] =   0.8454139853458534E-02;
  COFD[         356] =  -0.1908171031786331E+02;
  COFD[         357] =   0.4989633286602976E+01;
  COFD[         358] =  -0.4194071192530971E+00;
  COFD[         359] =   0.1762632910887997E-01;
  COFD[         360] =  -0.1511946285859442E+02;
  COFD[         361] =   0.3507927470352835E+01;
  COFD[         362] =  -0.2480755290441053E+00;
  COFD[         363] =   0.1106802953504252E-01;
  COFD[         364] =  -0.1559686009787643E+02;
  COFD[         365] =   0.3697254029246602E+01;
  COFD[         366] =  -0.2703332731709061E+00;
  COFD[         367] =   0.1193149792533076E-01;
  COFD[         368] =  -0.1559686009787643E+02;
  COFD[         369] =   0.3697254029246602E+01;
  COFD[         370] =  -0.2703332731709061E+00;
  COFD[         371] =   0.1193149792533076E-01;
  COFD[         372] =  -0.1562210625534589E+02;
  COFD[         373] =   0.3699988818188416E+01;
  COFD[         374] =  -0.2707081263160739E+00;
  COFD[         375] =   0.1194858184298051E-01;
  COFD[         376] =  -0.1556532980903788E+02;
  COFD[         377] =   0.3678661452138608E+01;
  COFD[         378] =  -0.2679777825820168E+00;
  COFD[         379] =   0.1183177506717497E-01;
  COFD[         380] =  -0.1476985219999908E+02;
  COFD[         381] =   0.3357751196108267E+01;
  COFD[         382] =  -0.2280046831004087E+00;
  COFD[         383] =   0.1017303667696902E-01;
  COFD[         384] =  -0.1792454027766008E+02;
  COFD[         385] =   0.4476326984714260E+01;
  COFD[         386] =  -0.3680191278135618E+00;
  COFD[         387] =   0.1601860672433683E-01;
  COFD[         388] =  -0.1918204714711579E+02;
  COFD[         389] =   0.4875615818554563E+01;
  COFD[         390] =  -0.4058262469897604E+00;
  COFD[         391] =   0.1707794755285167E-01;
  COFD[         392] =  -0.1918177070837826E+02;
  COFD[         393] =   0.4871873199573439E+01;
  COFD[         394] =  -0.4051602536432213E+00;
  COFD[         395] =   0.1704106539271711E-01;
  COFD[         396] =  -0.1880409131690931E+02;
  COFD[         397] =   0.4747896766189557E+01;
  COFD[         398] =  -0.3929513558466455E+00;
  COFD[         399] =   0.1667071927589917E-01;
  COFD[         400] =  -0.1799417924090561E+02;
  COFD[         401] =   0.4483237269191199E+01;
  COFD[         402] =  -0.3666457816700467E+00;
  COFD[         403] =   0.1586624138953086E-01;
  COFD[         404] =  -0.1800542957417654E+02;
  COFD[         405] =   0.4468504271448043E+01;
  COFD[         406] =  -0.3668692513471102E+00;
  COFD[         407] =   0.1596519376748332E-01;
  COFD[         408] =  -0.1803293840758895E+02;
  COFD[         409] =   0.4476898762718537E+01;
  COFD[         410] =  -0.3679481179513786E+00;
  COFD[         411] =   0.1601125468656841E-01;
  COFD[         412] =  -0.1473661449825940E+02;
  COFD[         413] =   0.3348204558833826E+01;
  COFD[         414] =  -0.2267271723233180E+00;
  COFD[         415] =   0.1011600240359858E-01;
  COFD[         416] =  -0.1580533640016843E+02;
  COFD[         417] =   0.3761415063014861E+01;
  COFD[         418] =  -0.2795015971789234E+00;
  COFD[         419] =   0.1236332420639666E-01;
  COFD[         420] =  -0.1818497883283899E+02;
  COFD[         421] =   0.5043624597310085E+01;
  COFD[         422] =  -0.4377066809210735E+00;
  COFD[         423] =   0.1887643866477002E-01;
  COFD[         424] =  -0.1651315089501800E+02;
  COFD[         425] =   0.4220701654162998E+01;
  COFD[         426] =  -0.2825947916575933E+00;
  COFD[         427] =   0.1018065117940972E-01;
  COFD[         428] =  -0.1906508600184398E+02;
  COFD[         429] =   0.4989055420187643E+01;
  COFD[         430] =  -0.4192832417141421E+00;
  COFD[         431] =   0.1761870656062123E-01;
  COFD[         432] =  -0.2029800570223752E+02;
  COFD[         433] =   0.5167991782317953E+01;
  COFD[         434] =  -0.4260628585538284E+00;
  COFD[         435] =   0.1725358239585169E-01;
  COFD[         436] =  -0.1908171031786331E+02;
  COFD[         437] =   0.4989633286602976E+01;
  COFD[         438] =  -0.4194071192530971E+00;
  COFD[         439] =   0.1762632910887997E-01;
  COFD[         440] =  -0.1180242969654581E+02;
  COFD[         441] =   0.8887910786873973E+00;
  COFD[         442] =   0.2461702200576275E+00;
  COFD[         443] =  -0.1607058638443385E-01;
  COFD[         444] =  -0.1969226895876598E+02;
  COFD[         445] =   0.4975023535625081E+01;
  COFD[         446] =  -0.4062730945223436E+00;
  COFD[         447] =   0.1660121367605663E-01;
  COFD[         448] =  -0.2040913503195636E+02;
  COFD[         449] =   0.5224092639655129E+01;
  COFD[         450] =  -0.4327353018434304E+00;
  COFD[         451] =   0.1753064056853202E-01;
  COFD[         452] =  -0.2040913503195636E+02;
  COFD[         453] =   0.5224092639655129E+01;
  COFD[         454] =  -0.4327353018434304E+00;
  COFD[         455] =   0.1753064056853202E-01;
  COFD[         456] =  -0.2042453573019913E+02;
  COFD[         457] =   0.5222435968901098E+01;
  COFD[         458] =  -0.4324987976330752E+00;
  COFD[         459] =   0.1751958312757915E-01;
  COFD[         460] =  -0.2052802769939090E+02;
  COFD[         461] =   0.5190809253464677E+01;
  COFD[         462] =  -0.4197156907449909E+00;
  COFD[         463] =   0.1661902490158712E-01;
  COFD[         464] =  -0.2031234004461199E+02;
  COFD[         465] =   0.5186312454417227E+01;
  COFD[         466] =  -0.4313742876779729E+00;
  COFD[         467] =   0.1759974844922337E-01;
  COFD[         468] =  -0.1924815690645235E+02;
  COFD[         469] =   0.4356731855020855E+01;
  COFD[         470] =  -0.2835834611821879E+00;
  COFD[         471] =   0.9647123135613308E-02;
  COFD[         472] =  -0.1757312087019186E+02;
  COFD[         473] =   0.3511204681563826E+01;
  COFD[         474] =  -0.1497966127823777E+00;
  COFD[         475] =   0.2966580974304371E-02;
  COFD[         476] =  -0.1757470287826831E+02;
  COFD[         477] =   0.3509436887842453E+01;
  COFD[         478] =  -0.1496102927025162E+00;
  COFD[         479] =   0.2961974572337444E-02;
  COFD[         480] =  -0.1619587514529141E+02;
  COFD[         481] =   0.2843235032975332E+01;
  COFD[         482] =  -0.5104314541884596E-01;
  COFD[         483] =  -0.1707937587486378E-02;
  COFD[         484] =  -0.1988688758541057E+02;
  COFD[         485] =   0.4655348588622681E+01;
  COFD[         486] =  -0.3287200044949488E+00;
  COFD[         487] =   0.1186464366263145E-01;
  COFD[         488] =  -0.2010033725040245E+02;
  COFD[         489] =   0.4740608584978043E+01;
  COFD[         490] =  -0.3442903017689571E+00;
  COFD[         491] =   0.1271241178230842E-01;
  COFD[         492] =  -0.2007381993778366E+02;
  COFD[         493] =   0.4726146050811121E+01;
  COFD[         494] =  -0.3422109775805109E+00;
  COFD[         495] =   0.1261488904944177E-01;
  COFD[         496] =  -0.2024331300258610E+02;
  COFD[         497] =   0.5167058816761965E+01;
  COFD[         498] =  -0.4292618020154789E+00;
  COFD[         499] =   0.1752039422653532E-01;
  COFD[         500] =  -0.1991427613539219E+02;
  COFD[         501] =   0.4985074718206674E+01;
  COFD[         502] =  -0.3988046631693363E+00;
  COFD[         503] =   0.1592993698509921E-01;
  COFD[         504] =  -0.1271124360389947E+02;
  COFD[         505] =   0.2931047949505553E+01;
  COFD[         506] =  -0.1741679825526774E+00;
  COFD[         507] =   0.7910585398301936E-02;
  COFD[         508] =  -0.1762111333528786E+02;
  COFD[         509] =   0.5052180357519856E+01;
  COFD[         510] =  -0.4496001366665920E+00;
  COFD[         511] =   0.1983696581079775E-01;
  COFD[         512] =  -0.1513524421958843E+02;
  COFD[         513] =   0.3523297948981071E+01;
  COFD[         514] =  -0.2502118585975575E+00;
  COFD[         515] =   0.1116667641067510E-01;
  COFD[         516] =  -0.1605887758547131E+02;
  COFD[         517] =   0.3688399632843664E+01;
  COFD[         518] =  -0.2691183502617316E+00;
  COFD[         519] =   0.1187611401222301E-01;
  COFD[         520] =  -0.1511946285859442E+02;
  COFD[         521] =   0.3507927470352835E+01;
  COFD[         522] =  -0.2480755290441053E+00;
  COFD[         523] =   0.1106802953504252E-01;
  COFD[         524] =  -0.1969226895876598E+02;
  COFD[         525] =   0.4975023535625081E+01;
  COFD[         526] =  -0.4062730945223436E+00;
  COFD[         527] =   0.1660121367605663E-01;
  COFD[         528] =  -0.1606627473213039E+02;
  COFD[         529] =   0.3688226344868143E+01;
  COFD[         530] =  -0.2690947493495643E+00;
  COFD[         531] =   0.1187504521154673E-01;
  COFD[         532] =  -0.1709853282395253E+02;
  COFD[         533] =   0.4139002691011624E+01;
  COFD[         534] =  -0.3268662406958736E+00;
  COFD[         535] =   0.1434738242897968E-01;
  COFD[         536] =  -0.1709853282395253E+02;
  COFD[         537] =   0.4139002691011624E+01;
  COFD[         538] =  -0.3268662406958736E+00;
  COFD[         539] =   0.1434738242897968E-01;
  COFD[         540] =  -0.1709003104221693E+02;
  COFD[         541] =   0.4125627998285617E+01;
  COFD[         542] =  -0.3250670331906144E+00;
  COFD[         543] =   0.1426686739847163E-01;
  COFD[         544] =  -0.1698006191382480E+02;
  COFD[         545] =   0.4079755891690070E+01;
  COFD[         546] =  -0.3190093940831370E+00;
  COFD[         547] =   0.1399976445173436E-01;
  COFD[         548] =  -0.1583284545979817E+02;
  COFD[         549] =   0.3600569993944359E+01;
  COFD[         550] =  -0.2583773478248265E+00;
  COFD[         551] =   0.1143956672273590E-01;
  COFD[         552] =  -0.1868972720602865E+02;
  COFD[         553] =   0.4564003101584411E+01;
  COFD[         554] =  -0.3752914154261894E+00;
  COFD[         555] =   0.1616756912994680E-01;
  COFD[         556] =  -0.2032173485009814E+02;
  COFD[         557] =   0.5112656067238100E+01;
  COFD[         558] =  -0.4317894593025599E+00;
  COFD[         559] =   0.1802136540679988E-01;
  COFD[         560] =  -0.2033308000404864E+02;
  COFD[         561] =   0.5113689932626547E+01;
  COFD[         562] =  -0.4319400850127593E+00;
  COFD[         563] =   0.1802856188319298E-01;
  COFD[         564] =  -0.1994316533928780E+02;
  COFD[         565] =   0.4989293390342089E+01;
  COFD[         566] =  -0.4200959859596722E+00;
  COFD[         567] =   0.1768643162186107E-01;
  COFD[         568] =  -0.1902570873205879E+02;
  COFD[         569] =   0.4693487292163261E+01;
  COFD[         570] =  -0.3900471625097403E+00;
  COFD[         571] =   0.1672372332121293E-01;
  COFD[         572] =  -0.1882326100183751E+02;
  COFD[         573] =   0.4594957876626870E+01;
  COFD[         574] =  -0.3790360690538448E+00;
  COFD[         575] =   0.1631950746180167E-01;
  COFD[         576] =  -0.1883388798798641E+02;
  COFD[         577] =   0.4595580957493091E+01;
  COFD[         578] =  -0.3791127554800276E+00;
  COFD[         579] =   0.1632261913288707E-01;
  COFD[         580] =  -0.1579929195048582E+02;
  COFD[         581] =   0.3590679879226839E+01;
  COFD[         582] =  -0.2570681340587627E+00;
  COFD[         583] =   0.1138172060899465E-01;
  COFD[         584] =  -0.1677798202265061E+02;
  COFD[         585] =   0.3944105793182093E+01;
  COFD[         586] =  -0.3009768138615855E+00;
  COFD[         587] =   0.1320048849753941E-01;
  COFD[         588] =  -0.1355980676188247E+02;
  COFD[         589] =   0.3231182295054575E+01;
  COFD[         590] =  -0.2138406992439825E+00;
  COFD[         591] =   0.9659145053095021E-02;
  COFD[         592] =  -0.1798724510881690E+02;
  COFD[         593] =   0.5079566464178115E+01;
  COFD[         594] =  -0.4444140659462563E+00;
  COFD[         595] =   0.1923462920552580E-01;
  COFD[         596] =  -0.1558834457822415E+02;
  COFD[         597] =   0.3699567664433178E+01;
  COFD[         598] =  -0.2706504848826866E+00;
  COFD[         599] =   0.1194595865045622E-01;
  COFD[         600] =  -0.1707915522153321E+02;
  COFD[         601] =   0.4132970234506635E+01;
  COFD[         602] =  -0.3260547355537204E+00;
  COFD[         603] =   0.1431106723651207E-01;
  COFD[         604] =  -0.1559686009787643E+02;
  COFD[         605] =   0.3697254029246602E+01;
  COFD[         606] =  -0.2703332731709061E+00;
  COFD[         607] =   0.1193149792533076E-01;
  COFD[         608] =  -0.2040913503195636E+02;
  COFD[         609] =   0.5224092639655129E+01;
  COFD[         610] =  -0.4327353018434304E+00;
  COFD[         611] =   0.1753064056853202E-01;
  COFD[         612] =  -0.1709853282395253E+02;
  COFD[         613] =   0.4139002691011624E+01;
  COFD[         614] =  -0.3268662406958736E+00;
  COFD[         615] =   0.1434738242897968E-01;
  COFD[         616] =  -0.1769521074851900E+02;
  COFD[         617] =   0.4367699379901299E+01;
  COFD[         618] =  -0.3537277912984029E+00;
  COFD[         619] =   0.1539755051214757E-01;
  COFD[         620] =  -0.1769521074851900E+02;
  COFD[         621] =   0.4367699379901299E+01;
  COFD[         622] =  -0.3537277912984029E+00;
  COFD[         623] =   0.1539755051214757E-01;
  COFD[         624] =  -0.1771391298979426E+02;
  COFD[         625] =   0.4368348686795481E+01;
  COFD[         626] =  -0.3538107613054961E+00;
  COFD[         627] =   0.1540107028020538E-01;
  COFD[         628] =  -0.1765053635794619E+02;
  COFD[         629] =   0.4344873087030904E+01;
  COFD[         630] =  -0.3508852749505070E+00;
  COFD[         631] =   0.1527908501173896E-01;
  COFD[         632] =  -0.1661775419297465E+02;
  COFD[         633] =   0.3942050022877351E+01;
  COFD[         634] =  -0.3009976504232300E+00;
  COFD[         635] =   0.1321284000123418E-01;
  COFD[         636] =  -0.1920160126704923E+02;
  COFD[         637] =   0.4768586004236171E+01;
  COFD[         638] =  -0.3938478869794211E+00;
  COFD[         639] =   0.1663350528215831E-01;
  COFD[         640] =  -0.2039806228291693E+02;
  COFD[         641] =   0.5097933786196286E+01;
  COFD[         642] =  -0.4183216752667500E+00;
  COFD[         643] =   0.1696948938817924E-01;
  COFD[         644] =  -0.2038049892591359E+02;
  COFD[         645] =   0.5087331121491625E+01;
  COFD[         646] =  -0.4167206268593121E+00;
  COFD[         647] =   0.1689044439564830E-01;
  COFD[         648] =  -0.2006425198525411E+02;
  COFD[         649] =   0.4999965267895282E+01;
  COFD[         650] =  -0.4103327875262441E+00;
  COFD[         651] =   0.1680582455088142E-01;
  COFD[         652] =  -0.1956675647234025E+02;
  COFD[         653] =   0.4893594117189979E+01;
  COFD[         654] =  -0.4083224878273328E+00;
  COFD[         655] =   0.1719480931404233E-01;
  COFD[         656] =  -0.1941109254354329E+02;
  COFD[         657] =   0.4819578196751006E+01;
  COFD[         658] =  -0.4008189828944452E+00;
  COFD[         659] =   0.1695355424415370E-01;
  COFD[         660] =  -0.1941401903486105E+02;
  COFD[         661] =   0.4817662369355102E+01;
  COFD[         662] =  -0.4004301095900044E+00;
  COFD[         663] =   0.1693043453375954E-01;
  COFD[         664] =  -0.1658740614019908E+02;
  COFD[         665] =   0.3933646317245189E+01;
  COFD[         666] =  -0.2999344650396615E+00;
  COFD[         667] =   0.1316797448289397E-01;
  COFD[         668] =  -0.1779539137381903E+02;
  COFD[         669] =   0.4385582017594726E+01;
  COFD[         670] =  -0.3562391667708982E+00;
  COFD[         671] =   0.1551072662807571E-01;
  COFD[         672] =  -0.1355980676188247E+02;
  COFD[         673] =   0.3231182295054575E+01;
  COFD[         674] =  -0.2138406992439825E+00;
  COFD[         675] =   0.9659145053095021E-02;
  COFD[         676] =  -0.1798724510881690E+02;
  COFD[         677] =   0.5079566464178115E+01;
  COFD[         678] =  -0.4444140659462563E+00;
  COFD[         679] =   0.1923462920552580E-01;
  COFD[         680] =  -0.1558834457822415E+02;
  COFD[         681] =   0.3699567664433178E+01;
  COFD[         682] =  -0.2706504848826866E+00;
  COFD[         683] =   0.1194595865045622E-01;
  COFD[         684] =  -0.1707915522153321E+02;
  COFD[         685] =   0.4132970234506635E+01;
  COFD[         686] =  -0.3260547355537204E+00;
  COFD[         687] =   0.1431106723651207E-01;
  COFD[         688] =  -0.1559686009787643E+02;
  COFD[         689] =   0.3697254029246602E+01;
  COFD[         690] =  -0.2703332731709061E+00;
  COFD[         691] =   0.1193149792533076E-01;
  COFD[         692] =  -0.2040913503195636E+02;
  COFD[         693] =   0.5224092639655129E+01;
  COFD[         694] =  -0.4327353018434304E+00;
  COFD[         695] =   0.1753064056853202E-01;
  COFD[         696] =  -0.1709853282395253E+02;
  COFD[         697] =   0.4139002691011624E+01;
  COFD[         698] =  -0.3268662406958736E+00;
  COFD[         699] =   0.1434738242897968E-01;
  COFD[         700] =  -0.1769521074851900E+02;
  COFD[         701] =   0.4367699379901299E+01;
  COFD[         702] =  -0.3537277912984029E+00;
  COFD[         703] =   0.1539755051214757E-01;
  COFD[         704] =  -0.1769521074851900E+02;
  COFD[         705] =   0.4367699379901299E+01;
  COFD[         706] =  -0.3537277912984029E+00;
  COFD[         707] =   0.1539755051214757E-01;
  COFD[         708] =  -0.1771391298979426E+02;
  COFD[         709] =   0.4368348686795481E+01;
  COFD[         710] =  -0.3538107613054961E+00;
  COFD[         711] =   0.1540107028020538E-01;
  COFD[         712] =  -0.1765053635794619E+02;
  COFD[         713] =   0.4344873087030904E+01;
  COFD[         714] =  -0.3508852749505070E+00;
  COFD[         715] =   0.1527908501173896E-01;
  COFD[         716] =  -0.1661775419297465E+02;
  COFD[         717] =   0.3942050022877351E+01;
  COFD[         718] =  -0.3009976504232300E+00;
  COFD[         719] =   0.1321284000123418E-01;
  COFD[         720] =  -0.1920160126704923E+02;
  COFD[         721] =   0.4768586004236171E+01;
  COFD[         722] =  -0.3938478869794211E+00;
  COFD[         723] =   0.1663350528215831E-01;
  COFD[         724] =  -0.2039806228291693E+02;
  COFD[         725] =   0.5097933786196286E+01;
  COFD[         726] =  -0.4183216752667500E+00;
  COFD[         727] =   0.1696948938817924E-01;
  COFD[         728] =  -0.2038049892591359E+02;
  COFD[         729] =   0.5087331121491625E+01;
  COFD[         730] =  -0.4167206268593121E+00;
  COFD[         731] =   0.1689044439564830E-01;
  COFD[         732] =  -0.2006425198525411E+02;
  COFD[         733] =   0.4999965267895282E+01;
  COFD[         734] =  -0.4103327875262441E+00;
  COFD[         735] =   0.1680582455088142E-01;
  COFD[         736] =  -0.1956675647234025E+02;
  COFD[         737] =   0.4893594117189979E+01;
  COFD[         738] =  -0.4083224878273328E+00;
  COFD[         739] =   0.1719480931404233E-01;
  COFD[         740] =  -0.1941109254354329E+02;
  COFD[         741] =   0.4819578196751006E+01;
  COFD[         742] =  -0.4008189828944452E+00;
  COFD[         743] =   0.1695355424415370E-01;
  COFD[         744] =  -0.1941401903486105E+02;
  COFD[         745] =   0.4817662369355102E+01;
  COFD[         746] =  -0.4004301095900044E+00;
  COFD[         747] =   0.1693043453375954E-01;
  COFD[         748] =  -0.1658740614019908E+02;
  COFD[         749] =   0.3933646317245189E+01;
  COFD[         750] =  -0.2999344650396615E+00;
  COFD[         751] =   0.1316797448289397E-01;
  COFD[         752] =  -0.1779539137381903E+02;
  COFD[         753] =   0.4385582017594726E+01;
  COFD[         754] =  -0.3562391667708982E+00;
  COFD[         755] =   0.1551072662807571E-01;
  COFD[         756] =  -0.1358415095584097E+02;
  COFD[         757] =   0.3239845970896248E+01;
  COFD[         758] =  -0.2150500406377150E+00;
  COFD[         759] =   0.9715225271588263E-02;
  COFD[         760] =  -0.1801636950310948E+02;
  COFD[         761] =   0.5090088857460223E+01;
  COFD[         762] =  -0.4457552498826564E+00;
  COFD[         763] =   0.1929141580461785E-01;
  COFD[         764] =  -0.1561687122378518E+02;
  COFD[         765] =   0.3703886367486989E+01;
  COFD[         766] =  -0.2712410029515650E+00;
  COFD[         767] =   0.1197280722366742E-01;
  COFD[         768] =  -0.1707094794267795E+02;
  COFD[         769] =   0.4119810624341889E+01;
  COFD[         770] =  -0.3242844616145015E+00;
  COFD[         771] =   0.1423184749658152E-01;
  COFD[         772] =  -0.1562210625534589E+02;
  COFD[         773] =   0.3699988818188416E+01;
  COFD[         774] =  -0.2707081263160739E+00;
  COFD[         775] =   0.1194858184298051E-01;
  COFD[         776] =  -0.2042453573019913E+02;
  COFD[         777] =   0.5222435968901098E+01;
  COFD[         778] =  -0.4324987976330752E+00;
  COFD[         779] =   0.1751958312757915E-01;
  COFD[         780] =  -0.1709003104221693E+02;
  COFD[         781] =   0.4125627998285617E+01;
  COFD[         782] =  -0.3250670331906144E+00;
  COFD[         783] =   0.1426686739847163E-01;
  COFD[         784] =  -0.1771391298979426E+02;
  COFD[         785] =   0.4368348686795481E+01;
  COFD[         786] =  -0.3538107613054961E+00;
  COFD[         787] =   0.1540107028020538E-01;
  COFD[         788] =  -0.1771391298979426E+02;
  COFD[         789] =   0.4368348686795481E+01;
  COFD[         790] =  -0.3538107613054961E+00;
  COFD[         791] =   0.1540107028020538E-01;
  COFD[         792] =  -0.1772990792056323E+02;
  COFD[         793] =   0.4367699379901269E+01;
  COFD[         794] =  -0.3537277912983985E+00;
  COFD[         795] =   0.1539755051214736E-01;
  COFD[         796] =  -0.1766459414605863E+02;
  COFD[         797] =   0.4343242779501765E+01;
  COFD[         798] =  -0.3506769601296933E+00;
  COFD[         799] =   0.1527024825917116E-01;
  COFD[         800] =  -0.1661214673917502E+02;
  COFD[         801] =   0.3930381940303202E+01;
  COFD[         802] =  -0.2994310262293675E+00;
  COFD[         803] =   0.1314286546846328E-01;
  COFD[         804] =  -0.1923136095279614E+02;
  COFD[         805] =   0.4771652981887103E+01;
  COFD[         806] =  -0.3944968044729790E+00;
  COFD[         807] =   0.1667283377652949E-01;
  COFD[         808] =  -0.2046646394524606E+02;
  COFD[         809] =   0.5118752767847491E+01;
  COFD[         810] =  -0.4214659652131379E+00;
  COFD[         811] =   0.1712474984073264E-01;
  COFD[         812] =  -0.2045023591963550E+02;
  COFD[         813] =   0.5108641161759671E+01;
  COFD[         814] =  -0.4199387136777197E+00;
  COFD[         815] =   0.1704933199723114E-01;
  COFD[         816] =  -0.2012688925964889E+02;
  COFD[         817] =   0.5018247199767933E+01;
  COFD[         818] =  -0.4131414445091446E+00;
  COFD[         819] =   0.1694646037906594E-01;
  COFD[         820] =  -0.1960184896363511E+02;
  COFD[         821] =   0.4900346283733346E+01;
  COFD[         822] =  -0.4094945328875742E+00;
  COFD[         823] =   0.1725872606725634E-01;
  COFD[         824] =  -0.1943925299002717E+02;
  COFD[         825] =   0.4823356728877616E+01;
  COFD[         826] =  -0.4015874038077185E+00;
  COFD[         827] =   0.1699928308073373E-01;
  COFD[         828] =  -0.1944256474839717E+02;
  COFD[         829] =   0.4821522129341849E+01;
  COFD[         830] =  -0.4012140535427826E+00;
  COFD[         831] =   0.1697705724216803E-01;
  COFD[         832] =  -0.1658229952875101E+02;
  COFD[         833] =   0.3922184822536529E+01;
  COFD[         834] =  -0.2983959485115160E+00;
  COFD[         835] =   0.1309927190981370E-01;
  COFD[         836] =  -0.1779710613097560E+02;
  COFD[         837] =   0.4376269742117769E+01;
  COFD[         838] =  -0.3550500828076509E+00;
  COFD[         839] =   0.1546030915508793E-01;
  COFD[         840] =  -0.1353824314289343E+02;
  COFD[         841] =   0.3227352939332374E+01;
  COFD[         842] =  -0.2134245159235849E+00;
  COFD[         843] =   0.9644648995332943E-02;
  COFD[         844] =  -0.1790236652645725E+02;
  COFD[         845] =   0.5054241512059871E+01;
  COFD[         846] =  -0.4413163862361035E+00;
  COFD[         847] =   0.1910841720545356E-01;
  COFD[         848] =  -0.1556217965348612E+02;
  COFD[         849] =   0.3683628994574988E+01;
  COFD[         850] =  -0.2686556693855442E+00;
  COFD[         851] =   0.1186253607235680E-01;
  COFD[         852] =  -0.1696131964138645E+02;
  COFD[         853] =   0.4074167708212651E+01;
  COFD[         854] =  -0.3182576642303228E+00;
  COFD[         855] =   0.1396612551510442E-01;
  COFD[         856] =  -0.1556532980903788E+02;
  COFD[         857] =   0.3678661452138608E+01;
  COFD[         858] =  -0.2679777825820168E+00;
  COFD[         859] =   0.1183177506717497E-01;
  COFD[         860] =  -0.2052802769939090E+02;
  COFD[         861] =   0.5190809253464677E+01;
  COFD[         862] =  -0.4197156907449909E+00;
  COFD[         863] =   0.1661902490158712E-01;
  COFD[         864] =  -0.1698006191382480E+02;
  COFD[         865] =   0.4079755891690070E+01;
  COFD[         866] =  -0.3190093940831370E+00;
  COFD[         867] =   0.1399976445173436E-01;
  COFD[         868] =  -0.1765053635794619E+02;
  COFD[         869] =   0.4344873087030904E+01;
  COFD[         870] =  -0.3508852749505070E+00;
  COFD[         871] =   0.1527908501173896E-01;
  COFD[         872] =  -0.1765053635794619E+02;
  COFD[         873] =   0.4344873087030904E+01;
  COFD[         874] =  -0.3508852749505070E+00;
  COFD[         875] =   0.1527908501173896E-01;
  COFD[         876] =  -0.1766459414605863E+02;
  COFD[         877] =   0.4343242779501765E+01;
  COFD[         878] =  -0.3506769601296933E+00;
  COFD[         879] =   0.1527024825917116E-01;
  COFD[         880] =  -0.1760280207090635E+02;
  COFD[         881] =   0.4320037087778842E+01;
  COFD[         882] =  -0.3477962090917468E+00;
  COFD[         883] =   0.1515062790890059E-01;
  COFD[         884] =  -0.1654504240089281E+02;
  COFD[         885] =   0.3902990496709335E+01;
  COFD[         886] =  -0.2959726841289695E+00;
  COFD[         887] =   0.1299731318228341E-01;
  COFD[         888] =  -0.1920918070133062E+02;
  COFD[         889] =   0.4764713252836613E+01;
  COFD[         890] =  -0.3942405534249802E+00;
  COFD[         891] =   0.1668923743050566E-01;
  COFD[         892] =  -0.2048733868393580E+02;
  COFD[         893] =   0.5132889365248658E+01;
  COFD[         894] =  -0.4243039026916984E+00;
  COFD[         895] =   0.1728905581961421E-01;
  COFD[         896] =  -0.2047278277859408E+02;
  COFD[         897] =   0.5123414930187060E+01;
  COFD[         898] =  -0.4228706903648477E+00;
  COFD[         899] =   0.1721819830097161E-01;
  COFD[         900] =  -0.2034321636694622E+02;
  COFD[         901] =   0.5089562300279502E+01;
  COFD[         902] =  -0.4211891334428361E+00;
  COFD[         903] =   0.1724909778947470E-01;
  COFD[         904] =  -0.1957015754946189E+02;
  COFD[         905] =   0.4890650507507935E+01;
  COFD[         906] =  -0.4088565156459115E+00;
  COFD[         907] =   0.1725720932243744E-01;
  COFD[         908] =  -0.1938496483487797E+02;
  COFD[         909] =   0.4803651837482498E+01;
  COFD[         910] =  -0.3995406103851766E+00;
  COFD[         911] =   0.1693221082596063E-01;
  COFD[         912] =  -0.1938900888397246E+02;
  COFD[         913] =   0.4802051950424429E+01;
  COFD[         914] =  -0.3992044152028926E+00;
  COFD[         915] =   0.1691189193038190E-01;
  COFD[         916] =  -0.1651845590523054E+02;
  COFD[         917] =   0.3896132878391479E+01;
  COFD[         918] =  -0.2951170694006329E+00;
  COFD[         919] =   0.1296170338310411E-01;
  COFD[         920] =  -0.1772757250170383E+02;
  COFD[         921] =   0.4347725932428951E+01;
  COFD[         922] =  -0.3515015798062593E+00;
  COFD[         923] =   0.1531308129463464E-01;
  COFD[         924] =  -0.1268296956267843E+02;
  COFD[         925] =   0.2903428401599904E+01;
  COFD[         926] =  -0.1707667302851327E+00;
  COFD[         927] =   0.7771111802295727E-02;
  COFD[         928] =  -0.1719132602269712E+02;
  COFD[         929] =   0.4859367175853364E+01;
  COFD[         930] =  -0.4240580935572132E+00;
  COFD[         931] =   0.1870696403708276E-01;
  COFD[         932] =  -0.1478198899238034E+02;
  COFD[         933] =   0.3371092581023257E+01;
  COFD[         934] =  -0.2298567856863946E+00;
  COFD[         935] =   0.1025846107153025E-01;
  COFD[         936] =  -0.1582244348422915E+02;
  COFD[         937] =   0.3599211636009877E+01;
  COFD[         938] =  -0.2581898369286856E+00;
  COFD[         939] =   0.1143096575439008E-01;
  COFD[         940] =  -0.1476985219999908E+02;
  COFD[         941] =   0.3357751196108267E+01;
  COFD[         942] =  -0.2280046831004087E+00;
  COFD[         943] =   0.1017303667696902E-01;
  COFD[         944] =  -0.2031234004461199E+02;
  COFD[         945] =   0.5186312454417227E+01;
  COFD[         946] =  -0.4313742876779729E+00;
  COFD[         947] =   0.1759974844922337E-01;
  COFD[         948] =  -0.1583284545979817E+02;
  COFD[         949] =   0.3600569993944359E+01;
  COFD[         950] =  -0.2583773478248265E+00;
  COFD[         951] =   0.1143956672273590E-01;
  COFD[         952] =  -0.1661775419297465E+02;
  COFD[         953] =   0.3942050022877351E+01;
  COFD[         954] =  -0.3009976504232300E+00;
  COFD[         955] =   0.1321284000123418E-01;
  COFD[         956] =  -0.1661775419297465E+02;
  COFD[         957] =   0.3942050022877351E+01;
  COFD[         958] =  -0.3009976504232300E+00;
  COFD[         959] =   0.1321284000123418E-01;
  COFD[         960] =  -0.1661214673917502E+02;
  COFD[         961] =   0.3930381940303202E+01;
  COFD[         962] =  -0.2994310262293675E+00;
  COFD[         963] =   0.1314286546846328E-01;
  COFD[         964] =  -0.1654504240089281E+02;
  COFD[         965] =   0.3902990496709335E+01;
  COFD[         966] =  -0.2959726841289695E+00;
  COFD[         967] =   0.1299731318228341E-01;
  COFD[         968] =  -0.1551400389810888E+02;
  COFD[         969] =   0.3472777644794381E+01;
  COFD[         970] =  -0.2416432013010369E+00;
  COFD[         971] =   0.1070799930256162E-01;
  COFD[         972] =  -0.1853079633836986E+02;
  COFD[         973] =   0.4516962151189352E+01;
  COFD[         974] =  -0.3708427592293742E+00;
  COFD[         975] =   0.1604368344993687E-01;
  COFD[         976] =  -0.2009021453031469E+02;
  COFD[         977] =   0.5037634138971256E+01;
  COFD[         978] =  -0.4241196608853510E+00;
  COFD[         979] =   0.1777020759153531E-01;
  COFD[         980] =  -0.2009778046496907E+02;
  COFD[         981] =   0.5037181430781024E+01;
  COFD[         982] =  -0.4240362270814515E+00;
  COFD[         983] =   0.1776546650877937E-01;
  COFD[         984] =  -0.2003301781547066E+02;
  COFD[         985] =   0.5028625403684593E+01;
  COFD[         986] =  -0.4259384280972174E+00;
  COFD[         987] =   0.1797007932483027E-01;
  COFD[         988] =  -0.1879103964706305E+02;
  COFD[         989] =   0.4613167070471964E+01;
  COFD[         990] =  -0.3811821259847005E+00;
  COFD[         991] =   0.1640353922053438E-01;
  COFD[         992] =  -0.1864571977877140E+02;
  COFD[         993] =   0.4539069966095800E+01;
  COFD[         994] =  -0.3735011850255855E+00;
  COFD[         995] =   0.1615133908644331E-01;
  COFD[         996] =  -0.1865835471859333E+02;
  COFD[         997] =   0.4540734987509840E+01;
  COFD[         998] =  -0.3737078963922989E+00;
  COFD[         999] =   0.1615982140755382E-01;
  COFD[        1000] =  -0.1549017806948574E+02;
  COFD[        1001] =   0.3466859355125401E+01;
  COFD[        1002] =  -0.2408856054712291E+00;
  COFD[        1003] =   0.1067561900190753E-01;
  COFD[        1004] =  -0.1650781936547825E+02;
  COFD[        1005] =   0.3842122413055806E+01;
  COFD[        1006] =  -0.2882098997185555E+00;
  COFD[        1007] =   0.1266705263344362E-01;
  COFD[        1008] =  -0.1564739122681956E+02;
  COFD[        1009] =   0.4025825276486324E+01;
  COFD[        1010] =  -0.3181195147607377E+00;
  COFD[        1011] =   0.1422136468533644E-01;
  COFD[        1012] =  -0.1725034548774479E+02;
  COFD[        1013] =   0.4566921339972978E+01;
  COFD[        1014] =  -0.3509192271095338E+00;
  COFD[        1015] =   0.1402685086906688E-01;
  COFD[        1016] =  -0.1793673861488134E+02;
  COFD[        1017] =   0.4489875914439151E+01;
  COFD[        1018] =  -0.3697556337020395E+00;
  COFD[        1019] =   0.1609253879718407E-01;
  COFD[        1020] =  -0.1868514242925209E+02;
  COFD[        1021] =   0.4565512377700616E+01;
  COFD[        1022] =  -0.3754538805580755E+00;
  COFD[        1023] =   0.1617305773699574E-01;
  COFD[        1024] =  -0.1792454027766008E+02;
  COFD[        1025] =   0.4476326984714260E+01;
  COFD[        1026] =  -0.3680191278135618E+00;
  COFD[        1027] =   0.1601860672433683E-01;
  COFD[        1028] =  -0.1924815690645235E+02;
  COFD[        1029] =   0.4356731855020855E+01;
  COFD[        1030] =  -0.2835834611821879E+00;
  COFD[        1031] =   0.9647123135613308E-02;
  COFD[        1032] =  -0.1868972720602865E+02;
  COFD[        1033] =   0.4564003101584411E+01;
  COFD[        1034] =  -0.3752914154261894E+00;
  COFD[        1035] =   0.1616756912994680E-01;
  COFD[        1036] =  -0.1920160126704923E+02;
  COFD[        1037] =   0.4768586004236171E+01;
  COFD[        1038] =  -0.3938478869794211E+00;
  COFD[        1039] =   0.1663350528215831E-01;
  COFD[        1040] =  -0.1920160126704923E+02;
  COFD[        1041] =   0.4768586004236171E+01;
  COFD[        1042] =  -0.3938478869794211E+00;
  COFD[        1043] =   0.1663350528215831E-01;
  COFD[        1044] =  -0.1923136095279614E+02;
  COFD[        1045] =   0.4771652981887103E+01;
  COFD[        1046] =  -0.3944968044729790E+00;
  COFD[        1047] =   0.1667283377652949E-01;
  COFD[        1048] =  -0.1920918070133062E+02;
  COFD[        1049] =   0.4764713252836613E+01;
  COFD[        1050] =  -0.3942405534249802E+00;
  COFD[        1051] =   0.1668923743050566E-01;
  COFD[        1052] =  -0.1853079633836986E+02;
  COFD[        1053] =   0.4516962151189352E+01;
  COFD[        1054] =  -0.3708427592293742E+00;
  COFD[        1055] =   0.1604368344993687E-01;
  COFD[        1056] =  -0.2074682089064003E+02;
  COFD[        1057] =   0.5125913401864712E+01;
  COFD[        1058] =  -0.4300036496008145E+00;
  COFD[        1059] =   0.1780187794631543E-01;
  COFD[        1060] =  -0.2103480542370173E+02;
  COFD[        1061] =   0.5058663047911905E+01;
  COFD[        1062] =  -0.3955876986371539E+00;
  COFD[        1063] =   0.1531365269600822E-01;
  COFD[        1064] =  -0.2106203504941583E+02;
  COFD[        1065] =   0.5066201505062354E+01;
  COFD[        1066] =  -0.3966848046994042E+00;
  COFD[        1067] =   0.1536587788350550E-01;
  COFD[        1068] =  -0.2111193727907456E+02;
  COFD[        1069] =   0.5115004803331022E+01;
  COFD[        1070] =  -0.4073081193984433E+00;
  COFD[        1071] =   0.1598596435486438E-01;
  COFD[        1072] =  -0.2089328018920450E+02;
  COFD[        1073] =   0.5166930585415199E+01;
  COFD[        1074] =  -0.4306289589876920E+00;
  COFD[        1075] =   0.1764267219152353E-01;
  COFD[        1076] =  -0.2075862463488640E+02;
  COFD[        1077] =   0.5108417296497037E+01;
  COFD[        1078] =  -0.4262851069696628E+00;
  COFD[        1079] =   0.1758280414396813E-01;
  COFD[        1080] =  -0.2077840705218999E+02;
  COFD[        1081] =   0.5112879719806333E+01;
  COFD[        1082] =  -0.4269664397363305E+00;
  COFD[        1083] =   0.1761673548420319E-01;
  COFD[        1084] =  -0.1850244542172855E+02;
  COFD[        1085] =   0.4509546506902997E+01;
  COFD[        1086] =  -0.3699250247632380E+00;
  COFD[        1087] =   0.1600569104209367E-01;
  COFD[        1088] =  -0.1931016351705685E+02;
  COFD[        1089] =   0.4757139666603682E+01;
  COFD[        1090] =  -0.3961700036764624E+00;
  COFD[        1091] =   0.1690081368592412E-01;
  COFD[        1092] =  -0.1746999275651299E+02;
  COFD[        1093] =   0.4670556529920976E+01;
  COFD[        1094] =  -0.3931788949833612E+00;
  COFD[        1095] =   0.1709979614643272E-01;
  COFD[        1096] =  -0.1540645026491837E+02;
  COFD[        1097] =   0.3572427121876953E+01;
  COFD[        1098] =  -0.1878763535992870E+00;
  COFD[        1099] =   0.5591984032298197E-02;
  COFD[        1100] =  -0.1915075646861156E+02;
  COFD[        1101] =   0.4868838487330766E+01;
  COFD[        1102] =  -0.4046210142562804E+00;
  COFD[        1103] =   0.1701122981099066E-01;
  COFD[        1104] =  -0.2031657960333580E+02;
  COFD[        1105] =   0.5113614640951361E+01;
  COFD[        1106] =  -0.4319293017089550E+00;
  COFD[        1107] =   0.1802805482147992E-01;
  COFD[        1108] =  -0.1918204714711579E+02;
  COFD[        1109] =   0.4875615818554563E+01;
  COFD[        1110] =  -0.4058262469897604E+00;
  COFD[        1111] =   0.1707794755285167E-01;
  COFD[        1112] =  -0.1757312087019186E+02;
  COFD[        1113] =   0.3511204681563826E+01;
  COFD[        1114] =  -0.1497966127823777E+00;
  COFD[        1115] =   0.2966580974304371E-02;
  COFD[        1116] =  -0.2032173485009814E+02;
  COFD[        1117] =   0.5112656067238100E+01;
  COFD[        1118] =  -0.4317894593025599E+00;
  COFD[        1119] =   0.1802136540679988E-01;
  COFD[        1120] =  -0.2039806228291693E+02;
  COFD[        1121] =   0.5097933786196286E+01;
  COFD[        1122] =  -0.4183216752667500E+00;
  COFD[        1123] =   0.1696948938817924E-01;
  COFD[        1124] =  -0.2039806228291693E+02;
  COFD[        1125] =   0.5097933786196286E+01;
  COFD[        1126] =  -0.4183216752667500E+00;
  COFD[        1127] =   0.1696948938817924E-01;
  COFD[        1128] =  -0.2046646394524606E+02;
  COFD[        1129] =   0.5118752767847491E+01;
  COFD[        1130] =  -0.4214659652131379E+00;
  COFD[        1131] =   0.1712474984073264E-01;
  COFD[        1132] =  -0.2048733868393580E+02;
  COFD[        1133] =   0.5132889365248658E+01;
  COFD[        1134] =  -0.4243039026916984E+00;
  COFD[        1135] =   0.1728905581961421E-01;
  COFD[        1136] =  -0.2009021453031469E+02;
  COFD[        1137] =   0.5037634138971256E+01;
  COFD[        1138] =  -0.4241196608853510E+00;
  COFD[        1139] =   0.1777020759153531E-01;
  COFD[        1140] =  -0.2103480542370173E+02;
  COFD[        1141] =   0.5058663047911905E+01;
  COFD[        1142] =  -0.3955876986371539E+00;
  COFD[        1143] =   0.1531365269600822E-01;
  COFD[        1144] =  -0.1885553227873480E+02;
  COFD[        1145] =   0.3923352162852614E+01;
  COFD[        1146] =  -0.2115865709699218E+00;
  COFD[        1147] =   0.5939285843445688E-02;
  COFD[        1148] =  -0.1886381013480732E+02;
  COFD[        1149] =   0.3923279625734948E+01;
  COFD[        1150] =  -0.2115777534001486E+00;
  COFD[        1151] =   0.5938970665371634E-02;
  COFD[        1152] =  -0.1981646458797120E+02;
  COFD[        1153] =   0.4378886175261416E+01;
  COFD[        1154] =  -0.2813828834600260E+00;
  COFD[        1155] =   0.9370857616443682E-02;
  COFD[        1156] =  -0.2090933405715619E+02;
  COFD[        1157] =   0.4972999527197874E+01;
  COFD[        1158] =  -0.3786889263628220E+00;
  COFD[        1159] =   0.1435806171655258E-01;
  COFD[        1160] =  -0.2110606940488274E+02;
  COFD[        1161] =   0.5056386072237322E+01;
  COFD[        1162] =  -0.3941314241059736E+00;
  COFD[        1163] =   0.1520386988001248E-01;
  COFD[        1164] =  -0.2110836259086301E+02;
  COFD[        1165] =   0.5053674020415143E+01;
  COFD[        1166] =  -0.3937388243323189E+00;
  COFD[        1167] =   0.1518528106281032E-01;
  COFD[        1168] =  -0.2006577079410498E+02;
  COFD[        1169] =   0.5032536447034316E+01;
  COFD[        1170] =  -0.4235925902580250E+00;
  COFD[        1171] =   0.1775279409683131E-01;
  COFD[        1172] =  -0.2070496369832456E+02;
  COFD[        1173] =   0.5191072368701859E+01;
  COFD[        1174] =  -0.4346910627906122E+00;
  COFD[        1175] =   0.1785842811134005E-01;
  COFD[        1176] =  -0.1747771970561251E+02;
  COFD[        1177] =   0.4673159679696371E+01;
  COFD[        1178] =  -0.3935109715079382E+00;
  COFD[        1179] =   0.1711386921544412E-01;
  COFD[        1180] =  -0.1539084008249366E+02;
  COFD[        1181] =   0.3565006575777335E+01;
  COFD[        1182] =  -0.1867580729432450E+00;
  COFD[        1183] =   0.5536884066371989E-02;
  COFD[        1184] =  -0.1914986770399175E+02;
  COFD[        1185] =   0.4864899621234062E+01;
  COFD[        1186] =  -0.4039220297382060E+00;
  COFD[        1187] =   0.1697258735720344E-01;
  COFD[        1188] =  -0.2032689236579443E+02;
  COFD[        1189] =   0.5114221260960875E+01;
  COFD[        1190] =  -0.4320136109530544E+00;
  COFD[        1191] =   0.1803190494311217E-01;
  COFD[        1192] =  -0.1918177070837826E+02;
  COFD[        1193] =   0.4871873199573439E+01;
  COFD[        1194] =  -0.4051602536432213E+00;
  COFD[        1195] =   0.1704106539271711E-01;
  COFD[        1196] =  -0.1757470287826831E+02;
  COFD[        1197] =   0.3509436887842453E+01;
  COFD[        1198] =  -0.1496102927025162E+00;
  COFD[        1199] =   0.2961974572337444E-02;
  COFD[        1200] =  -0.2033308000404864E+02;
  COFD[        1201] =   0.5113689932626547E+01;
  COFD[        1202] =  -0.4319400850127593E+00;
  COFD[        1203] =   0.1802856188319298E-01;
  COFD[        1204] =  -0.2038049892591359E+02;
  COFD[        1205] =   0.5087331121491625E+01;
  COFD[        1206] =  -0.4167206268593121E+00;
  COFD[        1207] =   0.1689044439564830E-01;
  COFD[        1208] =  -0.2038049892591359E+02;
  COFD[        1209] =   0.5087331121491625E+01;
  COFD[        1210] =  -0.4167206268593121E+00;
  COFD[        1211] =   0.1689044439564830E-01;
  COFD[        1212] =  -0.2045023591963550E+02;
  COFD[        1213] =   0.5108641161759671E+01;
  COFD[        1214] =  -0.4199387136777197E+00;
  COFD[        1215] =   0.1704933199723114E-01;
  COFD[        1216] =  -0.2047278277859408E+02;
  COFD[        1217] =   0.5123414930187060E+01;
  COFD[        1218] =  -0.4228706903648477E+00;
  COFD[        1219] =   0.1721819830097161E-01;
  COFD[        1220] =  -0.2009778046496907E+02;
  COFD[        1221] =   0.5037181430781024E+01;
  COFD[        1222] =  -0.4240362270814515E+00;
  COFD[        1223] =   0.1776546650877937E-01;
  COFD[        1224] =  -0.2106203504941583E+02;
  COFD[        1225] =   0.5066201505062354E+01;
  COFD[        1226] =  -0.3966848046994042E+00;
  COFD[        1227] =   0.1536587788350550E-01;
  COFD[        1228] =  -0.1886381013480732E+02;
  COFD[        1229] =   0.3923279625734948E+01;
  COFD[        1230] =  -0.2115777534001486E+00;
  COFD[        1231] =   0.5938970665371634E-02;
  COFD[        1232] =  -0.1887260515065304E+02;
  COFD[        1233] =   0.3923352162852601E+01;
  COFD[        1234] =  -0.2115865709699200E+00;
  COFD[        1235] =   0.5939285843445600E-02;
  COFD[        1236] =  -0.1982650591715289E+02;
  COFD[        1237] =   0.4379424471801094E+01;
  COFD[        1238] =  -0.2814553408580349E+00;
  COFD[        1239] =   0.9373962429419786E-02;
  COFD[        1240] =  -0.2091775085800750E+02;
  COFD[        1241] =   0.4973029974129495E+01;
  COFD[        1242] =  -0.3786913370630364E+00;
  COFD[        1243] =   0.1435807083846824E-01;
  COFD[        1244] =  -0.2111907495038380E+02;
  COFD[        1245] =   0.5058381671032987E+01;
  COFD[        1246] =  -0.3944196800959486E+00;
  COFD[        1247] =   0.1521748580862574E-01;
  COFD[        1248] =  -0.2112313758034541E+02;
  COFD[        1249] =   0.5056389549388699E+01;
  COFD[        1250] =  -0.3941319269977342E+00;
  COFD[        1251] =   0.1520389366699413E-01;
  COFD[        1252] =  -0.2007323696432831E+02;
  COFD[        1253] =   0.5032033200707714E+01;
  COFD[        1254] =  -0.4235009227577490E+00;
  COFD[        1255] =   0.1774762291235686E-01;
  COFD[        1256] =  -0.2072332061257778E+02;
  COFD[        1257] =   0.5194960177957292E+01;
  COFD[        1258] =  -0.4352726149968372E+00;
  COFD[        1259] =   0.1788688372961993E-01;
  COFD[        1260] =  -0.1797854335218850E+02;
  COFD[        1261] =   0.4905020138101169E+01;
  COFD[        1262] =  -0.4293084201751438E+00;
  COFD[        1263] =   0.1891302999394448E-01;
  COFD[        1264] =  -0.1561125740910447E+02;
  COFD[        1265] =   0.3698144089263720E+01;
  COFD[        1266] =  -0.2103080207621357E+00;
  COFD[        1267] =   0.6792204226221277E-02;
  COFD[        1268] =  -0.1878355293596396E+02;
  COFD[        1269] =   0.4745970603113269E+01;
  COFD[        1270] =  -0.3924372163577848E+00;
  COFD[        1271] =   0.1663659729501466E-01;
  COFD[        1272] =  -0.2016013281728147E+02;
  COFD[        1273] =   0.5063566265537697E+01;
  COFD[        1274] =  -0.4281515792120217E+00;
  COFD[        1275] =   0.1797180686456706E-01;
  COFD[        1276] =  -0.1880409131690931E+02;
  COFD[        1277] =   0.4747896766189557E+01;
  COFD[        1278] =  -0.3929513558466455E+00;
  COFD[        1279] =   0.1667071927589917E-01;
  COFD[        1280] =  -0.1619587514529141E+02;
  COFD[        1281] =   0.2843235032975332E+01;
  COFD[        1282] =  -0.5104314541884596E-01;
  COFD[        1283] =  -0.1707937587486378E-02;
  COFD[        1284] =  -0.1994316533928780E+02;
  COFD[        1285] =   0.4989293390342089E+01;
  COFD[        1286] =  -0.4200959859596722E+00;
  COFD[        1287] =   0.1768643162186107E-01;
  COFD[        1288] =  -0.2006425198525411E+02;
  COFD[        1289] =   0.4999965267895282E+01;
  COFD[        1290] =  -0.4103327875262441E+00;
  COFD[        1291] =   0.1680582455088142E-01;
  COFD[        1292] =  -0.2006425198525411E+02;
  COFD[        1293] =   0.4999965267895282E+01;
  COFD[        1294] =  -0.4103327875262441E+00;
  COFD[        1295] =   0.1680582455088142E-01;
  COFD[        1296] =  -0.2012688925964889E+02;
  COFD[        1297] =   0.5018247199767933E+01;
  COFD[        1298] =  -0.4131414445091446E+00;
  COFD[        1299] =   0.1694646037906594E-01;
  COFD[        1300] =  -0.2034321636694622E+02;
  COFD[        1301] =   0.5089562300279502E+01;
  COFD[        1302] =  -0.4211891334428361E+00;
  COFD[        1303] =   0.1724909778947470E-01;
  COFD[        1304] =  -0.2003301781547066E+02;
  COFD[        1305] =   0.5028625403684593E+01;
  COFD[        1306] =  -0.4259384280972174E+00;
  COFD[        1307] =   0.1797007932483027E-01;
  COFD[        1308] =  -0.2111193727907456E+02;
  COFD[        1309] =   0.5115004803331022E+01;
  COFD[        1310] =  -0.4073081193984433E+00;
  COFD[        1311] =   0.1598596435486438E-01;
  COFD[        1312] =  -0.1981646458797120E+02;
  COFD[        1313] =   0.4378886175261416E+01;
  COFD[        1314] =  -0.2813828834600260E+00;
  COFD[        1315] =   0.9370857616443682E-02;
  COFD[        1316] =  -0.1982650591715289E+02;
  COFD[        1317] =   0.4379424471801094E+01;
  COFD[        1318] =  -0.2814553408580349E+00;
  COFD[        1319] =   0.9373962429419786E-02;
  COFD[        1320] =  -0.1906244143057827E+02;
  COFD[        1321] =   0.4049325847301449E+01;
  COFD[        1322] =  -0.2359825286728792E+00;
  COFD[        1323] =   0.7307807686876559E-02;
  COFD[        1324] =  -0.2114547150354106E+02;
  COFD[        1325] =   0.5122846710851252E+01;
  COFD[        1326] =  -0.4061105011633706E+00;
  COFD[        1327] =   0.1585327078701023E-01;
  COFD[        1328] =  -0.2123091991916349E+02;
  COFD[        1329] =   0.5158619088859661E+01;
  COFD[        1330] =  -0.4148430725527501E+00;
  COFD[        1331] =   0.1638779761015646E-01;
  COFD[        1332] =  -0.2123706489334075E+02;
  COFD[        1333] =   0.5157473227475789E+01;
  COFD[        1334] =  -0.4146762094312990E+00;
  COFD[        1335] =   0.1637983966492620E-01;
  COFD[        1336] =  -0.2000088490858467E+02;
  COFD[        1337] =   0.5021789747731799E+01;
  COFD[        1338] =  -0.4253480595637520E+00;
  COFD[        1339] =   0.1795667045275812E-01;
  COFD[        1340] =  -0.2043431836412989E+02;
  COFD[        1341] =   0.5116486808224121E+01;
  COFD[        1342] =  -0.4301998553589693E+00;
  COFD[        1343] =   0.1786838693497468E-01;
  COFD[        1344] =  -0.1646041046594534E+02;
  COFD[        1345] =   0.4329097783536894E+01;
  COFD[        1346] =  -0.3588538378541796E+00;
  COFD[        1347] =   0.1604603061706265E-01;
  COFD[        1348] =  -0.1694497223361135E+02;
  COFD[        1349] =   0.4353290721870846E+01;
  COFD[        1350] =  -0.3147441391485086E+00;
  COFD[        1351] =   0.1210201943729390E-01;
  COFD[        1352] =  -0.1800091724620567E+02;
  COFD[        1353] =   0.4493072014688184E+01;
  COFD[        1354] =  -0.3678451189139861E+00;
  COFD[        1355] =   0.1591448159607913E-01;
  COFD[        1356] =  -0.1901929468754062E+02;
  COFD[        1357] =   0.4693839754776671E+01;
  COFD[        1358] =  -0.3901000075815114E+00;
  COFD[        1359] =   0.1672631492507185E-01;
  COFD[        1360] =  -0.1799417924090561E+02;
  COFD[        1361] =   0.4483237269191199E+01;
  COFD[        1362] =  -0.3666457816700467E+00;
  COFD[        1363] =   0.1586624138953086E-01;
  COFD[        1364] =  -0.1988688758541057E+02;
  COFD[        1365] =   0.4655348588622681E+01;
  COFD[        1366] =  -0.3287200044949488E+00;
  COFD[        1367] =   0.1186464366263145E-01;
  COFD[        1368] =  -0.1902570873205879E+02;
  COFD[        1369] =   0.4693487292163261E+01;
  COFD[        1370] =  -0.3900471625097403E+00;
  COFD[        1371] =   0.1672372332121293E-01;
  COFD[        1372] =  -0.1956675647234025E+02;
  COFD[        1373] =   0.4893594117189979E+01;
  COFD[        1374] =  -0.4083224878273328E+00;
  COFD[        1375] =   0.1719480931404233E-01;
  COFD[        1376] =  -0.1956675647234025E+02;
  COFD[        1377] =   0.4893594117189979E+01;
  COFD[        1378] =  -0.4083224878273328E+00;
  COFD[        1379] =   0.1719480931404233E-01;
  COFD[        1380] =  -0.1960184896363511E+02;
  COFD[        1381] =   0.4900346283733346E+01;
  COFD[        1382] =  -0.4094945328875742E+00;
  COFD[        1383] =   0.1725872606725634E-01;
  COFD[        1384] =  -0.1957015754946189E+02;
  COFD[        1385] =   0.4890650507507935E+01;
  COFD[        1386] =  -0.4088565156459115E+00;
  COFD[        1387] =   0.1725720932243744E-01;
  COFD[        1388] =  -0.1879103964706305E+02;
  COFD[        1389] =   0.4613167070471964E+01;
  COFD[        1390] =  -0.3811821259847005E+00;
  COFD[        1391] =   0.1640353922053438E-01;
  COFD[        1392] =  -0.2089328018920450E+02;
  COFD[        1393] =   0.5166930585415199E+01;
  COFD[        1394] =  -0.4306289589876920E+00;
  COFD[        1395] =   0.1764267219152353E-01;
  COFD[        1396] =  -0.2090933405715619E+02;
  COFD[        1397] =   0.4972999527197874E+01;
  COFD[        1398] =  -0.3786889263628220E+00;
  COFD[        1399] =   0.1435806171655258E-01;
  COFD[        1400] =  -0.2091775085800750E+02;
  COFD[        1401] =   0.4973029974129495E+01;
  COFD[        1402] =  -0.3786913370630364E+00;
  COFD[        1403] =   0.1435807083846824E-01;
  COFD[        1404] =  -0.2114547150354106E+02;
  COFD[        1405] =   0.5122846710851252E+01;
  COFD[        1406] =  -0.4061105011633706E+00;
  COFD[        1407] =   0.1585327078701023E-01;
  COFD[        1408] =  -0.2117733141272616E+02;
  COFD[        1409] =   0.5262507718327111E+01;
  COFD[        1410] =  -0.4398579014089885E+00;
  COFD[        1411] =   0.1792489830320274E-01;
  COFD[        1412] =  -0.2108189705348579E+02;
  COFD[        1413] =   0.5219096733932689E+01;
  COFD[        1414] =  -0.4371166395377749E+00;
  COFD[        1415] =   0.1791475976287751E-01;
  COFD[        1416] =  -0.2108550712837031E+02;
  COFD[        1417] =   0.5216945829110916E+01;
  COFD[        1418] =  -0.4367926562187473E+00;
  COFD[        1419] =   0.1789880186533132E-01;
  COFD[        1420] =  -0.1876345278233955E+02;
  COFD[        1421] =   0.4606235123260915E+01;
  COFD[        1422] =  -0.3803641770348614E+00;
  COFD[        1423] =   0.1637150008428244E-01;
  COFD[        1424] =  -0.1962450162993201E+02;
  COFD[        1425] =   0.4875751864229406E+01;
  COFD[        1426] =  -0.4087050411778894E+00;
  COFD[        1427] =   0.1732743672805724E-01;
  COFD[        1428] =  -0.1615267929807825E+02;
  COFD[        1429] =   0.4174121490551778E+01;
  COFD[        1430] =  -0.3383475624805202E+00;
  COFD[        1431] =   0.1513993577326728E-01;
  COFD[        1432] =  -0.1738386031480564E+02;
  COFD[        1433] =   0.4515894896789344E+01;
  COFD[        1434] =  -0.3391856065220307E+00;
  COFD[        1435] =   0.1329276720561771E-01;
  COFD[        1436] =  -0.1802433609839966E+02;
  COFD[        1437] =   0.4483646303442399E+01;
  COFD[        1438] =  -0.3688151334706006E+00;
  COFD[        1439] =   0.1604826293273014E-01;
  COFD[        1440] =  -0.1881741117557303E+02;
  COFD[        1441] =   0.4595519619722809E+01;
  COFD[        1442] =  -0.3791053358880415E+00;
  COFD[        1443] =   0.1632232422675412E-01;
  COFD[        1444] =  -0.1800542957417654E+02;
  COFD[        1445] =   0.4468504271448043E+01;
  COFD[        1446] =  -0.3668692513471102E+00;
  COFD[        1447] =   0.1596519376748332E-01;
  COFD[        1448] =  -0.2010033725040245E+02;
  COFD[        1449] =   0.4740608584978043E+01;
  COFD[        1450] =  -0.3442903017689571E+00;
  COFD[        1451] =   0.1271241178230842E-01;
  COFD[        1452] =  -0.1882326100183751E+02;
  COFD[        1453] =   0.4594957876626870E+01;
  COFD[        1454] =  -0.3790360690538448E+00;
  COFD[        1455] =   0.1631950746180167E-01;
  COFD[        1456] =  -0.1941109254354329E+02;
  COFD[        1457] =   0.4819578196751006E+01;
  COFD[        1458] =  -0.4008189828944452E+00;
  COFD[        1459] =   0.1695355424415370E-01;
  COFD[        1460] =  -0.1941109254354329E+02;
  COFD[        1461] =   0.4819578196751006E+01;
  COFD[        1462] =  -0.4008189828944452E+00;
  COFD[        1463] =   0.1695355424415370E-01;
  COFD[        1464] =  -0.1943925299002717E+02;
  COFD[        1465] =   0.4823356728877616E+01;
  COFD[        1466] =  -0.4015874038077185E+00;
  COFD[        1467] =   0.1699928308073373E-01;
  COFD[        1468] =  -0.1938496483487797E+02;
  COFD[        1469] =   0.4803651837482498E+01;
  COFD[        1470] =  -0.3995406103851766E+00;
  COFD[        1471] =   0.1693221082596063E-01;
  COFD[        1472] =  -0.1864571977877140E+02;
  COFD[        1473] =   0.4539069966095800E+01;
  COFD[        1474] =  -0.3735011850255855E+00;
  COFD[        1475] =   0.1615133908644331E-01;
  COFD[        1476] =  -0.2075862463488640E+02;
  COFD[        1477] =   0.5108417296497037E+01;
  COFD[        1478] =  -0.4262851069696628E+00;
  COFD[        1479] =   0.1758280414396813E-01;
  COFD[        1480] =  -0.2110606940488274E+02;
  COFD[        1481] =   0.5056386072237322E+01;
  COFD[        1482] =  -0.3941314241059736E+00;
  COFD[        1483] =   0.1520386988001248E-01;
  COFD[        1484] =  -0.2111907495038380E+02;
  COFD[        1485] =   0.5058381671032987E+01;
  COFD[        1486] =  -0.3944196800959486E+00;
  COFD[        1487] =   0.1521748580862574E-01;
  COFD[        1488] =  -0.2123091991916349E+02;
  COFD[        1489] =   0.5158619088859661E+01;
  COFD[        1490] =  -0.4148430725527501E+00;
  COFD[        1491] =   0.1638779761015646E-01;
  COFD[        1492] =  -0.2108189705348579E+02;
  COFD[        1493] =   0.5219096733932689E+01;
  COFD[        1494] =  -0.4371166395377749E+00;
  COFD[        1495] =   0.1791475976287751E-01;
  COFD[        1496] =  -0.2091185026658962E+02;
  COFD[        1497] =   0.5148560230790555E+01;
  COFD[        1498] =  -0.4310562392852861E+00;
  COFD[        1499] =   0.1777358978417338E-01;
  COFD[        1500] =  -0.2091965964592144E+02;
  COFD[        1501] =   0.5148262389934279E+01;
  COFD[        1502] =  -0.4310108894089904E+00;
  COFD[        1503] =   0.1777133583281057E-01;
  COFD[        1504] =  -0.1861866210245701E+02;
  COFD[        1505] =   0.4532060302179115E+01;
  COFD[        1506] =  -0.3726516669477872E+00;
  COFD[        1507] =   0.1611699967368039E-01;
  COFD[        1508] =  -0.1943393371703941E+02;
  COFD[        1509] =   0.4784664368892784E+01;
  COFD[        1510] =  -0.3991077435475853E+00;
  COFD[        1511] =   0.1700210403998681E-01;
  COFD[        1512] =  -0.1616515483252968E+02;
  COFD[        1513] =   0.4178910660883203E+01;
  COFD[        1514] =  -0.3390062877759687E+00;
  COFD[        1515] =   0.1517006334795817E-01;
  COFD[        1516] =  -0.1738223879329547E+02;
  COFD[        1517] =   0.4514299549461490E+01;
  COFD[        1518] =  -0.3388651702152571E+00;
  COFD[        1519] =   0.1327383213284630E-01;
  COFD[        1520] =  -0.1805221831789452E+02;
  COFD[        1521] =   0.4492279525932131E+01;
  COFD[        1522] =  -0.3699242263001061E+00;
  COFD[        1523] =   0.1609559465701492E-01;
  COFD[        1524] =  -0.1882824574536174E+02;
  COFD[        1525] =   0.4596254296231404E+01;
  COFD[        1526] =  -0.3791928797888692E+00;
  COFD[        1527] =   0.1632573982847240E-01;
  COFD[        1528] =  -0.1803293840758895E+02;
  COFD[        1529] =   0.4476898762718537E+01;
  COFD[        1530] =  -0.3679481179513786E+00;
  COFD[        1531] =   0.1601125468656841E-01;
  COFD[        1532] =  -0.2007381993778366E+02;
  COFD[        1533] =   0.4726146050811121E+01;
  COFD[        1534] =  -0.3422109775805109E+00;
  COFD[        1535] =   0.1261488904944177E-01;
  COFD[        1536] =  -0.1883388798798641E+02;
  COFD[        1537] =   0.4595580957493091E+01;
  COFD[        1538] =  -0.3791127554800276E+00;
  COFD[        1539] =   0.1632261913288707E-01;
  COFD[        1540] =  -0.1941401903486105E+02;
  COFD[        1541] =   0.4817662369355102E+01;
  COFD[        1542] =  -0.4004301095900044E+00;
  COFD[        1543] =   0.1693043453375954E-01;
  COFD[        1544] =  -0.1941401903486105E+02;
  COFD[        1545] =   0.4817662369355102E+01;
  COFD[        1546] =  -0.4004301095900044E+00;
  COFD[        1547] =   0.1693043453375954E-01;
  COFD[        1548] =  -0.1944256474839717E+02;
  COFD[        1549] =   0.4821522129341849E+01;
  COFD[        1550] =  -0.4012140535427826E+00;
  COFD[        1551] =   0.1697705724216803E-01;
  COFD[        1552] =  -0.1938900888397246E+02;
  COFD[        1553] =   0.4802051950424429E+01;
  COFD[        1554] =  -0.3992044152028926E+00;
  COFD[        1555] =   0.1691189193038190E-01;
  COFD[        1556] =  -0.1865835471859333E+02;
  COFD[        1557] =   0.4540734987509840E+01;
  COFD[        1558] =  -0.3737078963922989E+00;
  COFD[        1559] =   0.1615982140755382E-01;
  COFD[        1560] =  -0.2077840705218999E+02;
  COFD[        1561] =   0.5112879719806333E+01;
  COFD[        1562] =  -0.4269664397363305E+00;
  COFD[        1563] =   0.1761673548420319E-01;
  COFD[        1564] =  -0.2110836259086301E+02;
  COFD[        1565] =   0.5053674020415143E+01;
  COFD[        1566] =  -0.3937388243323189E+00;
  COFD[        1567] =   0.1518528106281032E-01;
  COFD[        1568] =  -0.2112313758034541E+02;
  COFD[        1569] =   0.5056389549388699E+01;
  COFD[        1570] =  -0.3941319269977342E+00;
  COFD[        1571] =   0.1520389366699413E-01;
  COFD[        1572] =  -0.2123706489334075E+02;
  COFD[        1573] =   0.5157473227475789E+01;
  COFD[        1574] =  -0.4146762094312990E+00;
  COFD[        1575] =   0.1637983966492620E-01;
  COFD[        1576] =  -0.2108550712837031E+02;
  COFD[        1577] =   0.5216945829110916E+01;
  COFD[        1578] =  -0.4367926562187473E+00;
  COFD[        1579] =   0.1789880186533132E-01;
  COFD[        1580] =  -0.2091965964592144E+02;
  COFD[        1581] =   0.5148262389934279E+01;
  COFD[        1582] =  -0.4310108894089904E+00;
  COFD[        1583] =   0.1777133583281057E-01;
  COFD[        1584] =  -0.2092889793954542E+02;
  COFD[        1585] =   0.5148560230790550E+01;
  COFD[        1586] =  -0.4310562392852854E+00;
  COFD[        1587] =   0.1777358978417334E-01;
  COFD[        1588] =  -0.1863141958605046E+02;
  COFD[        1589] =   0.4533772315478532E+01;
  COFD[        1590] =  -0.3728642283605221E+00;
  COFD[        1591] =   0.1612572322702682E-01;
  COFD[        1592] =  -0.1944451603679417E+02;
  COFD[        1593] =   0.4785108067913721E+01;
  COFD[        1594] =  -0.3991858464067674E+00;
  COFD[        1595] =   0.1700639619802289E-01;
  COFD[        1596] =  -0.1266066516595179E+02;
  COFD[        1597] =   0.2898076057146818E+01;
  COFD[        1598] =  -0.1700453759318706E+00;
  COFD[        1599] =   0.7738690296212197E-02;
  COFD[        1600] =  -0.1711567017943538E+02;
  COFD[        1601] =   0.4833264586465570E+01;
  COFD[        1602] =  -0.4205908440542898E+00;
  COFD[        1603] =   0.1855345683322859E-01;
  COFD[        1604] =  -0.1474819919893963E+02;
  COFD[        1605] =   0.3361311502126538E+01;
  COFD[        1606] =  -0.2285465363406641E+00;
  COFD[        1607] =   0.1019990809984005E-01;
  COFD[        1608] =  -0.1578873902175925E+02;
  COFD[        1609] =   0.3589258101555060E+01;
  COFD[        1610] =  -0.2568718883472193E+00;
  COFD[        1611] =   0.1137271971440306E-01;
  COFD[        1612] =  -0.1473661449825940E+02;
  COFD[        1613] =   0.3348204558833826E+01;
  COFD[        1614] =  -0.2267271723233180E+00;
  COFD[        1615] =   0.1011600240359858E-01;
  COFD[        1616] =  -0.2024331300258610E+02;
  COFD[        1617] =   0.5167058816761965E+01;
  COFD[        1618] =  -0.4292618020154789E+00;
  COFD[        1619] =   0.1752039422653532E-01;
  COFD[        1620] =  -0.1579929195048582E+02;
  COFD[        1621] =   0.3590679879226839E+01;
  COFD[        1622] =  -0.2570681340587627E+00;
  COFD[        1623] =   0.1138172060899465E-01;
  COFD[        1624] =  -0.1658740614019908E+02;
  COFD[        1625] =   0.3933646317245189E+01;
  COFD[        1626] =  -0.2999344650396615E+00;
  COFD[        1627] =   0.1316797448289397E-01;
  COFD[        1628] =  -0.1658740614019908E+02;
  COFD[        1629] =   0.3933646317245189E+01;
  COFD[        1630] =  -0.2999344650396615E+00;
  COFD[        1631] =   0.1316797448289397E-01;
  COFD[        1632] =  -0.1658229952875101E+02;
  COFD[        1633] =   0.3922184822536529E+01;
  COFD[        1634] =  -0.2983959485115160E+00;
  COFD[        1635] =   0.1309927190981370E-01;
  COFD[        1636] =  -0.1651845590523054E+02;
  COFD[        1637] =   0.3896132878391479E+01;
  COFD[        1638] =  -0.2951170694006329E+00;
  COFD[        1639] =   0.1296170338310411E-01;
  COFD[        1640] =  -0.1549017806948574E+02;
  COFD[        1641] =   0.3466859355125401E+01;
  COFD[        1642] =  -0.2408856054712291E+00;
  COFD[        1643] =   0.1067561900190753E-01;
  COFD[        1644] =  -0.1850244542172855E+02;
  COFD[        1645] =   0.4509546506902997E+01;
  COFD[        1646] =  -0.3699250247632380E+00;
  COFD[        1647] =   0.1600569104209367E-01;
  COFD[        1648] =  -0.2006577079410498E+02;
  COFD[        1649] =   0.5032536447034316E+01;
  COFD[        1650] =  -0.4235925902580250E+00;
  COFD[        1651] =   0.1775279409683131E-01;
  COFD[        1652] =  -0.2007323696432831E+02;
  COFD[        1653] =   0.5032033200707714E+01;
  COFD[        1654] =  -0.4235009227577490E+00;
  COFD[        1655] =   0.1774762291235686E-01;
  COFD[        1656] =  -0.2000088490858467E+02;
  COFD[        1657] =   0.5021789747731799E+01;
  COFD[        1658] =  -0.4253480595637520E+00;
  COFD[        1659] =   0.1795667045275812E-01;
  COFD[        1660] =  -0.1876345278233955E+02;
  COFD[        1661] =   0.4606235123260915E+01;
  COFD[        1662] =  -0.3803641770348614E+00;
  COFD[        1663] =   0.1637150008428244E-01;
  COFD[        1664] =  -0.1861866210245701E+02;
  COFD[        1665] =   0.4532060302179115E+01;
  COFD[        1666] =  -0.3726516669477872E+00;
  COFD[        1667] =   0.1611699967368039E-01;
  COFD[        1668] =  -0.1863141958605046E+02;
  COFD[        1669] =   0.4533772315478532E+01;
  COFD[        1670] =  -0.3728642283605221E+00;
  COFD[        1671] =   0.1612572322702682E-01;
  COFD[        1672] =  -0.1546612761154138E+02;
  COFD[        1673] =   0.3460858520880728E+01;
  COFD[        1674] =  -0.2401164641793465E+00;
  COFD[        1675] =   0.1064270570979806E-01;
  COFD[        1676] =  -0.1648483159780408E+02;
  COFD[        1677] =   0.3836871316958563E+01;
  COFD[        1678] =  -0.2875628541513582E+00;
  COFD[        1679] =   0.1264044562964702E-01;
  COFD[        1680] =  -0.1342264266757377E+02;
  COFD[        1681] =   0.3221777556990338E+01;
  COFD[        1682] =  -0.2128673394230190E+00;
  COFD[        1683] =   0.9627787551744238E-02;
  COFD[        1684] =  -0.1752038279713604E+02;
  COFD[        1685] =   0.4959639086677385E+01;
  COFD[        1686] =  -0.4293974875408633E+00;
  COFD[        1687] =   0.1860806087249651E-01;
  COFD[        1688] =  -0.1581858256688188E+02;
  COFD[        1689] =   0.3775899039344615E+01;
  COFD[        1690] =  -0.2814687920927101E+00;
  COFD[        1691] =   0.1245218573218981E-01;
  COFD[        1692] =  -0.1677455230511400E+02;
  COFD[        1693] =   0.3946214344521416E+01;
  COFD[        1694] =  -0.3012613112287137E+00;
  COFD[        1695] =   0.1321325422150532E-01;
  COFD[        1696] =  -0.1580533640016843E+02;
  COFD[        1697] =   0.3761415063014861E+01;
  COFD[        1698] =  -0.2795015971789234E+00;
  COFD[        1699] =   0.1236332420639666E-01;
  COFD[        1700] =  -0.1991427613539219E+02;
  COFD[        1701] =   0.4985074718206674E+01;
  COFD[        1702] =  -0.3988046631693363E+00;
  COFD[        1703] =   0.1592993698509921E-01;
  COFD[        1704] =  -0.1677798202265061E+02;
  COFD[        1705] =   0.3944105793182093E+01;
  COFD[        1706] =  -0.3009768138615855E+00;
  COFD[        1707] =   0.1320048849753941E-01;
  COFD[        1708] =  -0.1779539137381903E+02;
  COFD[        1709] =   0.4385582017594726E+01;
  COFD[        1710] =  -0.3562391667708982E+00;
  COFD[        1711] =   0.1551072662807571E-01;
  COFD[        1712] =  -0.1779539137381903E+02;
  COFD[        1713] =   0.4385582017594726E+01;
  COFD[        1714] =  -0.3562391667708982E+00;
  COFD[        1715] =   0.1551072662807571E-01;
  COFD[        1716] =  -0.1779710613097560E+02;
  COFD[        1717] =   0.4376269742117769E+01;
  COFD[        1718] =  -0.3550500828076509E+00;
  COFD[        1719] =   0.1546030915508793E-01;
  COFD[        1720] =  -0.1772757250170383E+02;
  COFD[        1721] =   0.4347725932428951E+01;
  COFD[        1722] =  -0.3515015798062593E+00;
  COFD[        1723] =   0.1531308129463464E-01;
  COFD[        1724] =  -0.1650781936547825E+02;
  COFD[        1725] =   0.3842122413055806E+01;
  COFD[        1726] =  -0.2882098997185555E+00;
  COFD[        1727] =   0.1266705263344362E-01;
  COFD[        1728] =  -0.1931016351705685E+02;
  COFD[        1729] =   0.4757139666603682E+01;
  COFD[        1730] =  -0.3961700036764624E+00;
  COFD[        1731] =   0.1690081368592412E-01;
  COFD[        1732] =  -0.2070496369832456E+02;
  COFD[        1733] =   0.5191072368701859E+01;
  COFD[        1734] =  -0.4346910627906122E+00;
  COFD[        1735] =   0.1785842811134005E-01;
  COFD[        1736] =  -0.2072332061257778E+02;
  COFD[        1737] =   0.5194960177957292E+01;
  COFD[        1738] =  -0.4352726149968372E+00;
  COFD[        1739] =   0.1788688372961993E-01;
  COFD[        1740] =  -0.2043431836412989E+02;
  COFD[        1741] =   0.5116486808224121E+01;
  COFD[        1742] =  -0.4301998553589693E+00;
  COFD[        1743] =   0.1786838693497468E-01;
  COFD[        1744] =  -0.1962450162993201E+02;
  COFD[        1745] =   0.4875751864229406E+01;
  COFD[        1746] =  -0.4087050411778894E+00;
  COFD[        1747] =   0.1732743672805724E-01;
  COFD[        1748] =  -0.1943393371703941E+02;
  COFD[        1749] =   0.4784664368892784E+01;
  COFD[        1750] =  -0.3991077435475853E+00;
  COFD[        1751] =   0.1700210403998681E-01;
  COFD[        1752] =  -0.1944451603679417E+02;
  COFD[        1753] =   0.4785108067913721E+01;
  COFD[        1754] =  -0.3991858464067674E+00;
  COFD[        1755] =   0.1700639619802289E-01;
  COFD[        1756] =  -0.1648483159780408E+02;
  COFD[        1757] =   0.3836871316958563E+01;
  COFD[        1758] =  -0.2875628541513582E+00;
  COFD[        1759] =   0.1264044562964702E-01;
  COFD[        1760] =  -0.1771336337770503E+02;
  COFD[        1761] =   0.4289097727284545E+01;
  COFD[        1762] =  -0.3450079724750131E+00;
  COFD[        1763] =   0.1508113614518283E-01;
};  }




#if 0




\\
\\
\\  This is the mechanism file
\\
\\
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!     
!  Reduced version of GRI-MECH 1.2. 19 species ( + N2, AR); 84 reactions. !     
!                                 PennState Dec, 1994                     !     
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!     
ELEMENTS                                                                        
O  H  C  N  AR                                                                  
END                                                                             
SPECIES                                                                         
H2      H       O       O2      OH      H2O     HO2                             
CH2     CH2(S)  CH3     CH4     CO      CO2     HCO                             
CH2O    CH3O    C2H4    C2H5    C2H6                                            
N2      AR                                                                      
END                                                                             
REACTIONS                                                                       
O+H+M<=>OH+M                             5.000E+17   -1.000      0.00           
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
O+H2<=>H+OH                              5.000E+04    2.670   6290.00           
O+HO2<=>OH+O2                            2.000E+13    0.000      0.00           
O+CH2<=>H+HCO                            8.000E+13    0.000      0.00           
O+CH2(S)<=>H+HCO                         1.500E+13    0.000      0.00           
O+CH3<=>H+CH2O                           8.430E+13    0.000      0.00           
O+CH4<=>OH+CH3                           1.020E+09    1.500   8600.00           
O+CO+M<=>CO2+M                           6.020E+14    0.000   3000.00           
H2/2.00/ O2/6.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/3.50/ C2H6/3.00/ AR/0.50/    
O+HCO<=>OH+CO                            3.000E+13    0.000      0.00           
O+HCO<=>H+CO2                            3.000E+13    0.000      0.00           
O+CH2O<=>OH+HCO                          3.900E+13    0.000   3540.00           
O+C2H4<=>CH3+HCO                         1.920E+07    1.830    220.00           
O+C2H5<=>CH3+CH2O                        1.320E+14    0.000      0.00           
O+C2H6<=>OH+C2H5                         8.980E+07    1.920   5690.00           
O2+CO<=>O+CO2                            2.500E+12    0.000  47800.00           
O2+CH2O<=>HO2+HCO                        1.000E+14    0.000  40000.00           
H+O2+M<=>HO2+M                           2.800E+18   -0.860      0.00           
O2/0.00/ H2O/0.00/ CO/0.75/ CO2/1.50/ C2H6/1.50/ N2/0.00/ AR/0.00/              
H+2O2<=>HO2+O2                           3.000E+20   -1.720      0.00           
H+O2+H2O<=>HO2+H2O                       9.380E+18   -0.760      0.00           
H+O2+N2<=>HO2+N2                         3.750E+20   -1.720      0.00           
H+O2+AR<=>HO2+AR                         7.000E+17   -0.800      0.00           
H+O2<=>O+OH                              8.300E+13    0.000  14413.00           
2H+M<=>H2+M                              1.000E+18   -1.000      0.00           
H2/0.00/ H2O/0.00/ CH4/2.00/ CO2/0.00/ C2H6/3.00/ AR/0.63/                      
2H+H2<=>2H2                              9.000E+16   -0.600      0.00           
2H+H2O<=>H2+H2O                          6.000E+19   -1.250      0.00           
2H+CO2<=>H2+CO2                          5.500E+20   -2.000      0.00           
H+OH+M<=>H2O+M                           2.200E+22   -2.000      0.00           
H2/0.73/ H2O/3.65/ CH4/2.00/ C2H6/3.00/ AR/0.38/                                
H+HO2<=>O2+H2                            2.800E+13    0.000   1068.00           
H+HO2<=>2OH                              1.340E+14    0.000    635.00           
H+CH2(+M)<=>CH3(+M)                      2.500E+16   -0.800      0.00           
     LOW  /  3.200E+27   -3.140   1230.00/                                      
     TROE/  0.6800   78.00  1995.00  5590.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
H+CH3(+M)<=>CH4(+M)                      1.270E+16   -0.630    383.00           
     LOW  /  2.477E+33   -4.760   2440.00/                                      
     TROE/  0.7830   74.00  2941.00  6964.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
H+CH4<=>CH3+H2                           6.600E+08    1.620  10840.00           
H+HCO(+M)<=>CH2O(+M)                     1.090E+12    0.480   -260.00           
     LOW  /  1.350E+24   -2.570   1425.00/                                      
     TROE/  0.7824  271.00  2755.00  6570.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
H+HCO<=>H2+CO                            7.340E+13    0.000      0.00           
H+CH2O(+M)<=>CH3O(+M)                    5.400E+11    0.454   2600.00           
     LOW  /  2.200E+30   -4.800   5560.00/                                      
     TROE/  0.7580   94.00  1555.00  4200.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/                      
H+CH2O<=>HCO+H2                          2.300E+10    1.050   3275.00           
H+CH3O<=>OH+CH3                          3.200E+13    0.000      0.00           
H+C2H4(+M)<=>C2H5(+M)                    1.080E+12    0.454   1820.00           
     LOW  /  1.200E+42   -7.620   6970.00/                                      
     TROE/  0.9753  210.00   984.00  4374.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
H+C2H5(+M)<=>C2H6(+M)                    5.210E+17   -0.990   1580.00           
     LOW  /  1.990E+41   -7.080   6685.00/                                      
     TROE/  0.8422  125.00  2219.00  6882.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
H+C2H6<=>C2H5+H2                         1.150E+08    1.900   7530.00           
H2+CO(+M)<=>CH2O(+M)                     4.300E+07    1.500  79600.00           
     LOW  /  5.070E+27   -3.420  84350.00/                                      
     TROE/  0.9320  197.00  1540.00 10300.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
OH+H2<=>H+H2O                            2.160E+08    1.510   3430.00           
2OH<=>O+H2O                              3.570E+04    2.400  -2110.00           
OH+HO2<=>O2+H2O                          2.900E+13    0.000   -500.00           
OH+CH2<=>H+CH2O                          2.000E+13    0.000      0.00           
OH+CH2(S)<=>H+CH2O                       3.000E+13    0.000      0.00           
OH+CH3<=>CH2+H2O                         5.600E+07    1.600   5420.00           
OH+CH3<=>CH2(S)+H2O                      2.501E+13    0.000      0.00           
OH+CH4<=>CH3+H2O                         1.000E+08    1.600   3120.00           
OH+CO<=>H+CO2                            4.760E+07    1.228     70.00           
OH+HCO<=>H2O+CO                          5.000E+13    0.000      0.00           
OH+CH2O<=>HCO+H2O                        3.430E+09    1.180   -447.00           
OH+C2H6<=>C2H5+H2O                       3.540E+06    2.120    870.00           
HO2+CH2<=>OH+CH2O                        2.000E+13    0.000      0.00           
HO2+CH3<=>O2+CH4                         1.000E+12    0.000      0.00           
HO2+CH3<=>OH+CH3O                        2.000E+13    0.000      0.00           
HO2+CO<=>OH+CO2                          1.500E+14    0.000  23600.00           
CH2+O2<=>OH+HCO                          1.320E+13    0.000   1500.00           
CH2+H2<=>H+CH3                           5.000E+05    2.000   7230.00           
CH2+CH3<=>H+C2H4                         4.000E+13    0.000      0.00           
CH2+CH4<=>2CH3                           2.460E+06    2.000   8270.00           
CH2(S)+N2<=>CH2+N2                       1.500E+13    0.000    600.00           
CH2(S)+AR<=>CH2+AR                       9.000E+12    0.000    600.00           
CH2(S)+O2<=>H+OH+CO                      2.800E+13    0.000      0.00           
CH2(S)+O2<=>CO+H2O                       1.200E+13    0.000      0.00           
CH2(S)+H2<=>CH3+H                        7.000E+13    0.000      0.00           
CH2(S)+H2O<=>CH2+H2O                     3.000E+13    0.000      0.00           
CH2(S)+CH3<=>H+C2H4                      1.200E+13    0.000   -570.00           
CH2(S)+CH4<=>2CH3                        1.600E+13    0.000   -570.00           
CH2(S)+CO<=>CH2+CO                       9.000E+12    0.000      0.00           
CH2(S)+CO2<=>CH2+CO2                     7.000E+12    0.000      0.00           
CH2(S)+CO2<=>CO+CH2O                     1.400E+13    0.000      0.00           
CH3+O2<=>O+CH3O                          2.675E+13    0.000  28800.00           
CH3+O2<=>OH+CH2O                         3.600E+10    0.000   8940.00           
2CH3(+M)<=>C2H6(+M)                      2.120E+16   -0.970    620.00           
     LOW  /  1.770E+50   -9.670   6220.00/                                      
     TROE/  0.5325  151.00  1038.00  4970.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
2CH3<=>H+C2H5                            4.990E+12    0.100  10600.00           
CH3+HCO<=>CH4+CO                         2.648E+13    0.000      0.00           
CH3+CH2O<=>HCO+CH4                       3.320E+03    2.810   5860.00           
CH3+C2H6<=>C2H5+CH4                      6.140E+06    1.740  10450.00           
HCO+H2O<=>H+CO+H2O                       2.244E+18   -1.000  17000.00           
HCO+M<=>H+CO+M                           1.870E+17   -1.000  17000.00           
H2/2.00/ H2O/0.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/                      
HCO+O2<=>HO2+CO                          7.600E+12    0.000    400.00           
CH3O+O2<=>HO2+CH2O                       4.280E-13    7.600  -3530.00           
C2H5+O2<=>HO2+C2H4                       8.400E+11    0.000   3875.00           
END                                                                             

\\
\\
\\  This is the therm file
\\
\\
! GRI-MECH version 1.2 Thermodynamics released 11/16/94
! NASA Polynomial format for CHEMKIN-II
! CH2* symbol changed to CH2(S); only change from version 1.1
! see README file for disclaimer
THERMO ALL
300.      1000.     5000.
O                 L 1/90O   1   00   00   00G   200.000  3500.000  1000.000    1
 2.56942078E+00-8.59741137E-05 4.19484589E-08-1.00177799E-11 1.22833691E-15    2
 2.92175791E+04 4.78433864E+00 3.16826710E+00-3.27931884E-03 6.64306396E-06    3
-6.12806624E-09 2.11265971E-12 2.91222592E+04 2.05193346E+00 6.72540300E+03    4
O2                TPIS89O   2   00   00   00G   200.000  3500.000  1000.000    1
 3.28253784E+00 1.48308754E-03-7.57966669E-07 2.09470555E-10-2.16717794E-14    2
-1.08845772E+03 5.45323129E+00 3.78245636E+00-2.99673416E-03 9.84730201E-06    3
-9.68129509E-09 3.24372837E-12-1.06394356E+03 3.65767573E+00 8.68010400E+03    4
H                 L 7/88H   1   00   00   00G   200.000  3500.000   1000.00    1
 2.50000001E+00-2.30842973E-11 1.61561948E-14-4.73515235E-18 4.98197357E-22    2
 2.54736599E+04-4.46682914E-01 2.50000000E+00 7.05332819E-13-1.99591964E-15    3
 2.30081632E-18-9.27732332E-22 2.54736599E+04-4.46682853E-01 6.19742800E+03    4
H2                TPIS78H   2   00   00   00G   200.000  3500.000   1000.00    1
 3.33727920E+00-4.94024731E-05 4.99456778E-07-1.79566394E-10 2.00255376E-14    2
-9.50158922E+02-3.20502331E+00 2.34433112E+00 7.98052075E-03-1.94781510E-05    3
 2.01572094E-08-7.37611761E-12-9.17935173E+02 6.83010238E-01 8.46810200E+03    4
OH                RUS 78O   1H   1   00   00G   200.000  3500.000  1000.000    1
 3.09288767E+00 5.48429716E-04 1.26505228E-07-8.79461556E-11 1.17412376E-14    2
 3.85865700E+03 4.47669610E+00 3.99201543E+00-2.40131752E-03 4.61793841E-06    3
-3.88113333E-09 1.36411470E-12 3.61508056E+03-1.03925458E-01 8.81310600E+03    4
H2O               L 8/89H   2O   1   00   00G   200.000  3500.000  1000.000    1
 3.03399249E+00 2.17691804E-03-1.64072518E-07-9.70419870E-11 1.68200992E-14    2
-3.00042971E+04 4.96677010E+00 4.19864056E+00-2.03643410E-03 6.52040211E-06    3
-5.48797062E-09 1.77197817E-12-3.02937267E+04-8.49032208E-01 9.90409200E+03    4
HO2               L 5/89H   1O   2   00   00G   200.000  3500.000  1000.000    1
 4.01721090E+00 2.23982013E-03-6.33658150E-07 1.14246370E-10-1.07908535E-14    2
 1.11856713E+02 3.78510215E+00 4.30179801E+00-4.74912051E-03 2.11582891E-05    3
-2.42763894E-08 9.29225124E-12 2.94808040E+02 3.71666245E+00 1.00021620E+04    4
H2O2              L 7/88H   2O   2   00   00G   200.000  3500.000  1000.000    1
 4.16500285E+00 4.90831694E-03-1.90139225E-06 3.71185986E-10-2.87908305E-14    2
-1.78617877E+04 2.91615662E+00 4.27611269E+00-5.42822417E-04 1.67335701E-05    3
-2.15770813E-08 8.62454363E-12-1.77025821E+04 3.43505074E+00 1.11588350E+04    4
C                 L11/88C   1   00   00   00G   200.000  3500.000  1000.000    1
 2.49266888E+00 4.79889284E-05-7.24335020E-08 3.74291029E-11-4.87277893E-15    2
 8.54512953E+04 4.80150373E+00 2.55423955E+00-3.21537724E-04 7.33792245E-07    3
-7.32234889E-10 2.66521446E-13 8.54438832E+04 4.53130848E+00 6.53589500E+03    4
CH                TPIS79C   1H   1   00   00G   200.000  3500.000  1000.000    1
 2.87846473E+00 9.70913681E-04 1.44445655E-07-1.30687849E-10 1.76079383E-14    2
 7.10124364E+04 5.48497999E+00 3.48981665E+00 3.23835541E-04-1.68899065E-06    3
 3.16217327E-09-1.40609067E-12 7.07972934E+04 2.08401108E+00 8.62500000E+03    4
CH2               L S/93C   1H   2   00   00G   200.000  3500.000  1000.000    1
 2.87410113E+00 3.65639292E-03-1.40894597E-06 2.60179549E-10-1.87727567E-14    2
 4.62636040E+04 6.17119324E+00 3.76267867E+00 9.68872143E-04 2.79489841E-06    3
-3.85091153E-09 1.68741719E-12 4.60040401E+04 1.56253185E+00 1.00274170E+04    4
CH2(S)            L S/93C   1H   2   00   00G   200.000  3500.000  1000.000    1
 2.29203842E+00 4.65588637E-03-2.01191947E-06 4.17906000E-10-3.39716365E-14    2
 5.09259997E+04 8.62650169E+00 4.19860411E+00-2.36661419E-03 8.23296220E-06    3
-6.68815981E-09 1.94314737E-12 5.04968163E+04-7.69118967E-01 9.93967200E+03    4
CH3               L11/89C   1H   3   00   00G   200.000  3500.000  1000.000    1
 2.28571772E+00 7.23990037E-03-2.98714348E-06 5.95684644E-10-4.67154394E-14    2
 1.67755843E+04 8.48007179E+00 3.67359040E+00 2.01095175E-03 5.73021856E-06    3
-6.87117425E-09 2.54385734E-12 1.64449988E+04 1.60456433E+00 1.03663400E+04    4
CH4               L 8/88C   1H   4   00   00G   200.000  3500.000  1000.000    1
 7.48514950E-02 1.33909467E-02-5.73285809E-06 1.22292535E-09-1.01815230E-13    2
-9.46834459E+03 1.84373180E+01 5.14987613E+00-1.36709788E-02 4.91800599E-05    3
-4.84743026E-08 1.66693956E-11-1.02466476E+04-4.64130376E+00 1.00161980E+04    4
CO                TPIS79C   1O   1   00   00G   200.000  3500.000  1000.000    1
 2.71518561E+00 2.06252743E-03-9.98825771E-07 2.30053008E-10-2.03647716E-14    2
-1.41518724E+04 7.81868772E+00 3.57953347E+00-6.10353680E-04 1.01681433E-06    3
 9.07005884E-10-9.04424499E-13-1.43440860E+04 3.50840928E+00 8.67100000E+03    4
CO2               L 7/88C   1O   2   00   00G   200.000  3500.000  1000.000    1
 3.85746029E+00 4.41437026E-03-2.21481404E-06 5.23490188E-10-4.72084164E-14    2
-4.87591660E+04 2.27163806E+00 2.35677352E+00 8.98459677E-03-7.12356269E-06    3
 2.45919022E-09-1.43699548E-13-4.83719697E+04 9.90105222E+00 9.36546900E+03    4
HCO               L12/89H   1C   1O   1   00G   200.000  3500.000  1000.000    1
 2.77217438E+00 4.95695526E-03-2.48445613E-06 5.89161778E-10-5.33508711E-14    2
 4.01191815E+03 9.79834492E+00 4.22118584E+00-3.24392532E-03 1.37799446E-05    3
-1.33144093E-08 4.33768865E-12 3.83956496E+03 3.39437243E+00 9.98945000E+03    4
CH2O              L 8/88H   2C   1O   1   00G   200.000  3500.000  1000.000    1
 1.76069008E+00 9.20000082E-03-4.42258813E-06 1.00641212E-09-8.83855640E-14    2
-1.39958323E+04 1.36563230E+01 4.79372315E+00-9.90833369E-03 3.73220008E-05    3
-3.79285261E-08 1.31772652E-11-1.43089567E+04 6.02812900E-01 1.00197170E+04    4
CH2OH             GUNL93C   1H   3O   1   00G   200.000  3500.000 1000.0       1
 3.69266569E+00 8.64576797E-03-3.75101120E-06 7.87234636E-10-6.48554201E-14    2
-3.24250627E+03 5.81043215E+00 3.86388918E+00 5.59672304E-03 5.93271791E-06    3
-1.04532012E-08 4.36967278E-12-3.19391367E+03 5.47302243E+00 1.18339080E+04    4
CH3O              121686C   1H   3O   1     G  0300.00   3000.00  1000.00      1
 0.03770799E+02 0.07871497E-01-0.02656384E-04 0.03944431E-08-0.02112616E-12    2
 0.12783252E+03 0.02929575E+02 0.02106204E+02 0.07216595E-01 0.05338472E-04    3
-0.07377636E-07 0.02075610E-10 0.09786011E+04 0.13152177E+02                   4
CH3OH             L 8/88C   1H   4O   1   00G   200.000  3500.000  1000.000    1
 1.78970791E+00 1.40938292E-02-6.36500835E-06 1.38171085E-09-1.17060220E-13    2
-2.53748747E+04 1.45023623E+01 5.71539582E+00-1.52309129E-02 6.52441155E-05    3
-7.10806889E-08 2.61352698E-11-2.56427656E+04-1.50409823E+00 1.14352770E+04    4
C2H               L 1/91C   2H   1   00   00G   200.000  3500.000  1000.000    1
 3.16780652E+00 4.75221902E-03-1.83787077E-06 3.04190252E-10-1.77232770E-14    2
 6.71210650E+04 6.63589475E+00 2.88965733E+00 1.34099611E-02-2.84769501E-05    3
 2.94791045E-08-1.09331511E-11 6.68393932E+04 6.22296438E+00 1.04544720E+04    4
C2H2              L 1/91C   2H   2   00   00G   200.000  3500.000  1000.000    1
 4.14756964E+00 5.96166664E-03-2.37294852E-06 4.67412171E-10-3.61235213E-14    2
 2.59359992E+04-1.23028121E+00 8.08681094E-01 2.33615629E-02-3.55171815E-05    3
 2.80152437E-08-8.50072974E-12 2.64289807E+04 1.39397051E+01 1.00058390E+04    4
C2H3              L 2/92C   2H   3   00   00G   200.000  3500.000  1000.000    1
 3.01672400E+00 1.03302292E-02-4.68082349E-06 1.01763288E-09-8.62607041E-14    2
 3.46128739E+04 7.78732378E+00 3.21246645E+00 1.51479162E-03 2.59209412E-05    3
-3.57657847E-08 1.47150873E-11 3.48598468E+04 8.51054025E+00 1.05750490E+04    4
C2H4              L 1/91C   2H   4   00   00G   200.000  3500.000  1000.000    1
 2.03611116E+00 1.46454151E-02-6.71077915E-06 1.47222923E-09-1.25706061E-13    2
 4.93988614E+03 1.03053693E+01 3.95920148E+00-7.57052247E-03 5.70990292E-05    3
-6.91588753E-08 2.69884373E-11 5.08977593E+03 4.09733096E+00 1.05186890E+04    4
C2H5              L12/92C   2H   5   00   00G   200.000  3500.000  1000.000    1
 1.95465642E+00 1.73972722E-02-7.98206668E-06 1.75217689E-09-1.49641576E-13    2
 1.28575200E+04 1.34624343E+01 4.30646568E+00-4.18658892E-03 4.97142807E-05    3
-5.99126606E-08 2.30509004E-11 1.28416265E+04 4.70720924E+00 1.21852440E+04    4
C2H6              L 8/88C   2H   6   00   00G   200.000  3500.000  1000.000    1
 1.07188150E+00 2.16852677E-02-1.00256067E-05 2.21412001E-09-1.90002890E-13    2
-1.14263932E+04 1.51156107E+01 4.29142492E+00-5.50154270E-03 5.99438288E-05    3
-7.08466285E-08 2.68685771E-11-1.15222055E+04 2.66682316E+00 1.18915940E+04    4
CH2CO             L 5/90C   2H   2O   1   00G   200.000  3500.000  1000.000    1
 4.51129732E+00 9.00359745E-03-4.16939635E-06 9.23345882E-10-7.94838201E-14    2
-7.55105311E+03 6.32247205E-01 2.13583630E+00 1.81188721E-02-1.73947474E-05    3
 9.34397568E-09-2.01457615E-12-7.04291804E+03 1.22156480E+01 1.17977430E+04    4
HCCO              SRIC91H   1C   2O   1     G  0300.00   4000.00  1000.00      1
 0.56282058E+01 0.40853401E-02-0.15934547E-05 0.28626052E-09-0.19407832E-13    2
 0.19327215E+05-0.39302595E+01 0.22517214E+01 0.17655021E-01-0.23729101E-04    3
 0.17275759E-07-0.50664811E-11 0.20059449E+05 0.12490417E+02                   4
HCCOH              SRI91C   2O   1H   20   0G   300.000  5000.000   1000.      1
 0.59238291E+01 0.67923600E-02-0.25658564E-05 0.44987841E-09-0.29940101E-13    2
 0.72646260E+04-0.76017742E+01 0.12423733E+01 0.31072201E-01-0.50866864E-04    3
 0.43137131E-07-0.14014594E-10 0.80316143E+04 0.13874319E+02                   4
N2                121286N   2               G  0300.00   5000.00  1000.00      1
 0.02926640E+02 0.14879768E-02-0.05684760E-05 0.10097038E-09-0.06753351E-13    2
-0.09227977E+04 0.05980528E+02 0.03298677E+02 0.14082404E-02-0.03963222E-04    3
 0.05641515E-07-0.02444854E-10-0.10208999E+04 0.03950372E+02                   4
AR                120186AR  1               G  0300.00   5000.00  1000.00      1
 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-0.07453750E+04 0.04366000E+02 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-0.07453750E+04 0.04366000E+02                   4
END

\\
\\
\\  This is the tran file
\\
\\
AR                 0   136.500     3.330     0.000     0.000     0.000
C                  0    71.400     3.298     0.000     0.000     0.000 ! *
C2                 1    97.530     3.621     0.000     1.760     4.000
C2O                1   232.400     3.828     0.000     0.000     1.000 ! *
CN2                1   232.400     3.828     0.000     0.000     1.000 ! OIS
C2H                1   209.000     4.100     0.000     0.000     2.500
C2H2               1   209.000     4.100     0.000     0.000     2.500
C2H2OH             2   224.700     4.162     0.000     0.000     1.000 ! *
C2H3               2   209.000     4.100     0.000     0.000     1.000 ! *
C2H4               2   280.800     3.971     0.000     0.000     1.500
C2H5               2   252.300     4.302     0.000     0.000     1.500
C2H6               2   252.300     4.302     0.000     0.000     1.500
C2N                1   232.400     3.828     0.000     0.000     1.000 ! OIS
C2N2               1   349.000     4.361     0.000     0.000     1.000 ! OIS
C3H2               2   209.000     4.100     0.000     0.000     1.000 ! *
C3H4               1   252.000     4.760     0.000     0.000     1.000
C3H6               2   266.800     4.982     0.000     0.000     1.000
C3H7               2   266.800     4.982     0.000     0.000     1.000
C4H6               2   357.000     5.180     0.000     0.000     1.000
I*C3H7             2   266.800     4.982     0.000     0.000     1.000
N*C3H7             2   266.800     4.982     0.000     0.000     1.000
C3H8               2   266.800     4.982     0.000     0.000     1.000
C4H                1   357.000     5.180     0.000     0.000     1.000
C4H2               1   357.000     5.180     0.000     0.000     1.000
C4H2OH             2   224.700     4.162     0.000     0.000     1.000 ! *
C4H8               2   357.000     5.176     0.000     0.000     1.000
C4H9               2   357.000     5.176     0.000     0.000     1.000
I*C4H9             2   357.000     5.176     0.000     0.000     1.000
C5H2               1   357.000     5.180     0.000     0.000     1.000
C5H3               1   357.000     5.180     0.000     0.000     1.000
C6H2               1   357.000     5.180     0.000     0.000     1.000
C6H5               2   412.300     5.349     0.000     0.000     1.000 ! JAM
C6H5O              2   450.000     5.500     0.000     0.000     1.000 ! JAM
C5H5OH             2   450.000     5.500     0.000     0.000     1.000 ! JAM
C6H6               2   412.300     5.349     0.000     0.000     1.000 ! SVE
C6H7               2   412.300     5.349     0.000     0.000     1.000 ! JAM
CH                 1    80.000     2.750     0.000     0.000     0.000
CH2                1   144.000     3.800     0.000     0.000     0.000
CH2(S)             1   144.000     3.800     0.000     0.000     0.000
CH2*               1   144.000     3.800     0.000     0.000     0.000
CH2CHCCH           2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH2CHCCH2          2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH2CHCH2           2   260.000     4.850     0.000     0.000     1.000 ! JAM
CH2CHCHCH          2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH2CHCHCH2         2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH2CO              2   436.000     3.970     0.000     0.000     2.000
CH2O               2   498.000     3.590     0.000     0.000     2.000
CH2OH              2   417.000     3.690     1.700     0.000     2.000
CH3                1   144.000     3.800     0.000     0.000     0.000
CH3CC              2   252.000     4.760     0.000     0.000     1.000 ! JAM
CH3CCCH2           2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH3CCCH3           2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH3CCH2            2   260.000     4.850     0.000     0.000     1.000 ! JAM
CH3CHCH            2   260.000     4.850     0.000     0.000     1.000 ! JAM
CH3CH2CCH          2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH3CHO             2   436.000     3.970     0.000     0.000     2.000
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
H                  0   145.000     2.050     0.000     0.000     0.000
H2C4O              2   357.000     5.180     0.000     0.000     1.000 ! JAM
H2                 1    38.000     2.920     0.000     0.790   280.000
H2CCCCH            2   357.000     5.180     0.000     0.000     1.000 ! JAM
H2CCCCH2           2   357.000     5.180     0.000     0.000     1.000 ! JAM
H2CCCH             2   252.000     4.760     0.000     0.000     1.000 ! JAM
H2CN               1   569.000     3.630     0.000     0.000     1.000 ! os/jm
H2NO               2   116.700     3.492     0.000     0.000     1.000 ! JAM
H2O                2   572.400     2.605     1.844     0.000     4.000
H2O2               2   107.400     3.458     0.000     0.000     3.800
HC2N2              1   349.000     4.361     0.000     0.000     1.000 ! OIS
HCCHCCH            2   357.000     5.180     0.000     0.000     1.000 ! JAM
HCCO               2   150.000     2.500     0.000     0.000     1.000 ! *
HCNN               2   150.000     2.500     0.000     0.000     1.000 ! *
HCCOH              2   436.000     3.970     0.000     0.000     2.000
HCN                1   569.000     3.630     0.000     0.000     1.000 ! OIS
HCO                2   498.000     3.590     0.000     0.000     0.000
HE                 0    10.200     2.576     0.000     0.000     0.000 ! *
HCNO               2   232.400     3.828     0.000     0.000     1.000 ! JAM
HOCN               2   232.400     3.828     0.000     0.000     1.000 ! JAM
HNCO               2   232.400     3.828     0.000     0.000     1.000 ! OIS
HNNO               2   232.400     3.828     0.000     0.000     1.000 ! *
HNO                2   116.700     3.492     0.000     0.000     1.000 ! *
HNOH               2   116.700     3.492     0.000     0.000     1.000 ! JAM
HO2                2   107.400     3.458     0.000     0.000     1.000 ! *
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
OH                 1    80.000     2.750     0.000     0.000     0.000

#endif
