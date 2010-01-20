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
    *kk = 38;
    *ii = 228;
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

    /* h  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = ' ';

    /* o  */
    kname[ 1*lenkname + 0 ] = 'O';
    kname[ 1*lenkname + 1 ] = ' ';

    /* c  */
    kname[ 2*lenkname + 0 ] = 'C';
    kname[ 2*lenkname + 1 ] = ' ';

    /* n  */
    kname[ 3*lenkname + 0 ] = 'N';
    kname[ 3*lenkname + 1 ] = ' ';

}


/* Returns the char strings of species names */
void CKSYMS(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*38; i++) {
        kname[i] = ' ';
    }

    /* h2  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = '2';
    kname[ 0*lenkname + 2 ] = ' ';

    /* h  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = ' ';

    /* ch4  */
    kname[ 2*lenkname + 0 ] = 'C';
    kname[ 2*lenkname + 1 ] = 'H';
    kname[ 2*lenkname + 2 ] = '4';
    kname[ 2*lenkname + 3 ] = ' ';

    /* ch3  */
    kname[ 3*lenkname + 0 ] = 'C';
    kname[ 3*lenkname + 1 ] = 'H';
    kname[ 3*lenkname + 2 ] = '3';
    kname[ 3*lenkname + 3 ] = ' ';

    /* ch2  */
    kname[ 4*lenkname + 0 ] = 'C';
    kname[ 4*lenkname + 1 ] = 'H';
    kname[ 4*lenkname + 2 ] = '2';
    kname[ 4*lenkname + 3 ] = ' ';

    /* ch  */
    kname[ 5*lenkname + 0 ] = 'C';
    kname[ 5*lenkname + 1 ] = 'H';
    kname[ 5*lenkname + 2 ] = ' ';

    /* ch2o  */
    kname[ 6*lenkname + 0 ] = 'C';
    kname[ 6*lenkname + 1 ] = 'H';
    kname[ 6*lenkname + 2 ] = '2';
    kname[ 6*lenkname + 3 ] = 'O';
    kname[ 6*lenkname + 4 ] = ' ';

    /* hco  */
    kname[ 7*lenkname + 0 ] = 'H';
    kname[ 7*lenkname + 1 ] = 'C';
    kname[ 7*lenkname + 2 ] = 'O';
    kname[ 7*lenkname + 3 ] = ' ';

    /* co2  */
    kname[ 8*lenkname + 0 ] = 'C';
    kname[ 8*lenkname + 1 ] = 'O';
    kname[ 8*lenkname + 2 ] = '2';
    kname[ 8*lenkname + 3 ] = ' ';

    /* co  */
    kname[ 9*lenkname + 0 ] = 'C';
    kname[ 9*lenkname + 1 ] = 'O';
    kname[ 9*lenkname + 2 ] = ' ';

    /* o2  */
    kname[ 10*lenkname + 0 ] = 'O';
    kname[ 10*lenkname + 1 ] = '2';
    kname[ 10*lenkname + 2 ] = ' ';

    /* o  */
    kname[ 11*lenkname + 0 ] = 'O';
    kname[ 11*lenkname + 1 ] = ' ';

    /* oh  */
    kname[ 12*lenkname + 0 ] = 'O';
    kname[ 12*lenkname + 1 ] = 'H';
    kname[ 12*lenkname + 2 ] = ' ';

    /* ho2  */
    kname[ 13*lenkname + 0 ] = 'H';
    kname[ 13*lenkname + 1 ] = 'O';
    kname[ 13*lenkname + 2 ] = '2';
    kname[ 13*lenkname + 3 ] = ' ';

    /* h2o2  */
    kname[ 14*lenkname + 0 ] = 'H';
    kname[ 14*lenkname + 1 ] = '2';
    kname[ 14*lenkname + 2 ] = 'O';
    kname[ 14*lenkname + 3 ] = '2';
    kname[ 14*lenkname + 4 ] = ' ';

    /* h2o  */
    kname[ 15*lenkname + 0 ] = 'H';
    kname[ 15*lenkname + 1 ] = '2';
    kname[ 15*lenkname + 2 ] = 'O';
    kname[ 15*lenkname + 3 ] = ' ';

    /* c2h  */
    kname[ 16*lenkname + 0 ] = 'C';
    kname[ 16*lenkname + 1 ] = '2';
    kname[ 16*lenkname + 2 ] = 'H';
    kname[ 16*lenkname + 3 ] = ' ';

    /* hcco  */
    kname[ 17*lenkname + 0 ] = 'H';
    kname[ 17*lenkname + 1 ] = 'C';
    kname[ 17*lenkname + 2 ] = 'C';
    kname[ 17*lenkname + 3 ] = 'O';
    kname[ 17*lenkname + 4 ] = ' ';

    /* c2h2  */
    kname[ 18*lenkname + 0 ] = 'C';
    kname[ 18*lenkname + 1 ] = '2';
    kname[ 18*lenkname + 2 ] = 'H';
    kname[ 18*lenkname + 3 ] = '2';
    kname[ 18*lenkname + 4 ] = ' ';

    /* c2h3  */
    kname[ 19*lenkname + 0 ] = 'C';
    kname[ 19*lenkname + 1 ] = '2';
    kname[ 19*lenkname + 2 ] = 'H';
    kname[ 19*lenkname + 3 ] = '3';
    kname[ 19*lenkname + 4 ] = ' ';

    /* c2h4  */
    kname[ 20*lenkname + 0 ] = 'C';
    kname[ 20*lenkname + 1 ] = '2';
    kname[ 20*lenkname + 2 ] = 'H';
    kname[ 20*lenkname + 3 ] = '4';
    kname[ 20*lenkname + 4 ] = ' ';

    /* c2h5  */
    kname[ 21*lenkname + 0 ] = 'C';
    kname[ 21*lenkname + 1 ] = '2';
    kname[ 21*lenkname + 2 ] = 'H';
    kname[ 21*lenkname + 3 ] = '5';
    kname[ 21*lenkname + 4 ] = ' ';

    /* c2h6  */
    kname[ 22*lenkname + 0 ] = 'C';
    kname[ 22*lenkname + 1 ] = '2';
    kname[ 22*lenkname + 2 ] = 'H';
    kname[ 22*lenkname + 3 ] = '6';
    kname[ 22*lenkname + 4 ] = ' ';

    /* ch2oh  */
    kname[ 23*lenkname + 0 ] = 'C';
    kname[ 23*lenkname + 1 ] = 'H';
    kname[ 23*lenkname + 2 ] = '2';
    kname[ 23*lenkname + 3 ] = 'O';
    kname[ 23*lenkname + 4 ] = 'H';
    kname[ 23*lenkname + 5 ] = ' ';

    /* ch3o  */
    kname[ 24*lenkname + 0 ] = 'C';
    kname[ 24*lenkname + 1 ] = 'H';
    kname[ 24*lenkname + 2 ] = '3';
    kname[ 24*lenkname + 3 ] = 'O';
    kname[ 24*lenkname + 4 ] = ' ';

    /* hccoh  */
    kname[ 25*lenkname + 0 ] = 'H';
    kname[ 25*lenkname + 1 ] = 'C';
    kname[ 25*lenkname + 2 ] = 'C';
    kname[ 25*lenkname + 3 ] = 'O';
    kname[ 25*lenkname + 4 ] = 'H';
    kname[ 25*lenkname + 5 ] = ' ';

    /* h2ccch  */
    kname[ 26*lenkname + 0 ] = 'H';
    kname[ 26*lenkname + 1 ] = '2';
    kname[ 26*lenkname + 2 ] = 'C';
    kname[ 26*lenkname + 3 ] = 'C';
    kname[ 26*lenkname + 4 ] = 'C';
    kname[ 26*lenkname + 5 ] = 'H';
    kname[ 26*lenkname + 6 ] = ' ';

    /* c3h2  */
    kname[ 27*lenkname + 0 ] = 'C';
    kname[ 27*lenkname + 1 ] = '3';
    kname[ 27*lenkname + 2 ] = 'H';
    kname[ 27*lenkname + 3 ] = '2';
    kname[ 27*lenkname + 4 ] = ' ';

    /* ch2s  */
    kname[ 28*lenkname + 0 ] = 'C';
    kname[ 28*lenkname + 1 ] = 'H';
    kname[ 28*lenkname + 2 ] = '2';
    kname[ 28*lenkname + 3 ] = 'S';
    kname[ 28*lenkname + 4 ] = ' ';

    /* ch2co  */
    kname[ 29*lenkname + 0 ] = 'C';
    kname[ 29*lenkname + 1 ] = 'H';
    kname[ 29*lenkname + 2 ] = '2';
    kname[ 29*lenkname + 3 ] = 'C';
    kname[ 29*lenkname + 4 ] = 'O';
    kname[ 29*lenkname + 5 ] = ' ';

    /* ch2hco  */
    kname[ 30*lenkname + 0 ] = 'C';
    kname[ 30*lenkname + 1 ] = 'H';
    kname[ 30*lenkname + 2 ] = '2';
    kname[ 30*lenkname + 3 ] = 'H';
    kname[ 30*lenkname + 4 ] = 'C';
    kname[ 30*lenkname + 5 ] = 'O';
    kname[ 30*lenkname + 6 ] = ' ';

    /* ch3co  */
    kname[ 31*lenkname + 0 ] = 'C';
    kname[ 31*lenkname + 1 ] = 'H';
    kname[ 31*lenkname + 2 ] = '3';
    kname[ 31*lenkname + 3 ] = 'C';
    kname[ 31*lenkname + 4 ] = 'O';
    kname[ 31*lenkname + 5 ] = ' ';

    /* ch3hco  */
    kname[ 32*lenkname + 0 ] = 'C';
    kname[ 32*lenkname + 1 ] = 'H';
    kname[ 32*lenkname + 2 ] = '3';
    kname[ 32*lenkname + 3 ] = 'H';
    kname[ 32*lenkname + 4 ] = 'C';
    kname[ 32*lenkname + 5 ] = 'O';
    kname[ 32*lenkname + 6 ] = ' ';

    /* c2h5oh  */
    kname[ 33*lenkname + 0 ] = 'C';
    kname[ 33*lenkname + 1 ] = '2';
    kname[ 33*lenkname + 2 ] = 'H';
    kname[ 33*lenkname + 3 ] = '5';
    kname[ 33*lenkname + 4 ] = 'O';
    kname[ 33*lenkname + 5 ] = 'H';
    kname[ 33*lenkname + 6 ] = ' ';

    /* c2h4oh  */
    kname[ 34*lenkname + 0 ] = 'C';
    kname[ 34*lenkname + 1 ] = '2';
    kname[ 34*lenkname + 2 ] = 'H';
    kname[ 34*lenkname + 3 ] = '4';
    kname[ 34*lenkname + 4 ] = 'O';
    kname[ 34*lenkname + 5 ] = 'H';
    kname[ 34*lenkname + 6 ] = ' ';

    /* ch3choh  */
    kname[ 35*lenkname + 0 ] = 'C';
    kname[ 35*lenkname + 1 ] = 'H';
    kname[ 35*lenkname + 2 ] = '3';
    kname[ 35*lenkname + 3 ] = 'C';
    kname[ 35*lenkname + 4 ] = 'H';
    kname[ 35*lenkname + 5 ] = 'O';
    kname[ 35*lenkname + 6 ] = 'H';
    kname[ 35*lenkname + 7 ] = ' ';

    /* ch3ch2o  */
    kname[ 36*lenkname + 0 ] = 'C';
    kname[ 36*lenkname + 1 ] = 'H';
    kname[ 36*lenkname + 2 ] = '3';
    kname[ 36*lenkname + 3 ] = 'C';
    kname[ 36*lenkname + 4 ] = 'H';
    kname[ 36*lenkname + 5 ] = '2';
    kname[ 36*lenkname + 6 ] = 'O';
    kname[ 36*lenkname + 7 ] = ' ';

    /* n2  */
    kname[ 37*lenkname + 0 ] = 'N';
    kname[ 37*lenkname + 1 ] = '2';
    kname[ 37*lenkname + 2 ] = ' ';

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
    XW += x[0]*2.015940; /*h2 */
    XW += x[1]*1.007970; /*h */
    XW += x[2]*16.043030; /*ch4 */
    XW += x[3]*15.035060; /*ch3 */
    XW += x[4]*14.027090; /*ch2 */
    XW += x[5]*13.019120; /*ch */
    XW += x[6]*30.026490; /*ch2o */
    XW += x[7]*29.018520; /*hco */
    XW += x[8]*44.009950; /*co2 */
    XW += x[9]*28.010550; /*co */
    XW += x[10]*31.998800; /*o2 */
    XW += x[11]*15.999400; /*o */
    XW += x[12]*17.007370; /*oh */
    XW += x[13]*33.006770; /*ho2 */
    XW += x[14]*34.014740; /*h2o2 */
    XW += x[15]*18.015340; /*h2o */
    XW += x[16]*25.030270; /*c2h */
    XW += x[17]*41.029670; /*hcco */
    XW += x[18]*26.038240; /*c2h2 */
    XW += x[19]*27.046210; /*c2h3 */
    XW += x[20]*28.054180; /*c2h4 */
    XW += x[21]*29.062150; /*c2h5 */
    XW += x[22]*30.070120; /*c2h6 */
    XW += x[23]*31.034460; /*ch2oh */
    XW += x[24]*31.034460; /*ch3o */
    XW += x[25]*42.037640; /*hccoh */
    XW += x[26]*39.057360; /*h2ccch */
    XW += x[27]*38.049390; /*c3h2 */
    XW += x[28]*14.027090; /*ch2s */
    XW += x[29]*42.037640; /*ch2co */
    XW += x[30]*43.045610; /*ch2hco */
    XW += x[31]*43.045610; /*ch3co */
    XW += x[32]*44.053580; /*ch3hco */
    XW += x[33]*46.069520; /*c2h5oh */
    XW += x[34]*45.061550; /*c2h4oh */
    XW += x[35]*45.061550; /*ch3choh */
    XW += x[36]*45.061550; /*ch3ch2o */
    XW += x[37]*28.013400; /*n2 */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void CKPY(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]/2.015940; /*h2 */
    YOW += y[1]/1.007970; /*h */
    YOW += y[2]/16.043030; /*ch4 */
    YOW += y[3]/15.035060; /*ch3 */
    YOW += y[4]/14.027090; /*ch2 */
    YOW += y[5]/13.019120; /*ch */
    YOW += y[6]/30.026490; /*ch2o */
    YOW += y[7]/29.018520; /*hco */
    YOW += y[8]/44.009950; /*co2 */
    YOW += y[9]/28.010550; /*co */
    YOW += y[10]/31.998800; /*o2 */
    YOW += y[11]/15.999400; /*o */
    YOW += y[12]/17.007370; /*oh */
    YOW += y[13]/33.006770; /*ho2 */
    YOW += y[14]/34.014740; /*h2o2 */
    YOW += y[15]/18.015340; /*h2o */
    YOW += y[16]/25.030270; /*c2h */
    YOW += y[17]/41.029670; /*hcco */
    YOW += y[18]/26.038240; /*c2h2 */
    YOW += y[19]/27.046210; /*c2h3 */
    YOW += y[20]/28.054180; /*c2h4 */
    YOW += y[21]/29.062150; /*c2h5 */
    YOW += y[22]/30.070120; /*c2h6 */
    YOW += y[23]/31.034460; /*ch2oh */
    YOW += y[24]/31.034460; /*ch3o */
    YOW += y[25]/42.037640; /*hccoh */
    YOW += y[26]/39.057360; /*h2ccch */
    YOW += y[27]/38.049390; /*c3h2 */
    YOW += y[28]/14.027090; /*ch2s */
    YOW += y[29]/42.037640; /*ch2co */
    YOW += y[30]/43.045610; /*ch2hco */
    YOW += y[31]/43.045610; /*ch3co */
    YOW += y[32]/44.053580; /*ch3hco */
    YOW += y[33]/46.069520; /*c2h5oh */
    YOW += y[34]/45.061550; /*c2h4oh */
    YOW += y[35]/45.061550; /*ch3choh */
    YOW += y[36]/45.061550; /*ch3ch2o */
    YOW += y[37]/28.013400; /*n2 */
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
    W += c[0]*2.015940; /*h2 */
    W += c[1]*1.007970; /*h */
    W += c[2]*16.043030; /*ch4 */
    W += c[3]*15.035060; /*ch3 */
    W += c[4]*14.027090; /*ch2 */
    W += c[5]*13.019120; /*ch */
    W += c[6]*30.026490; /*ch2o */
    W += c[7]*29.018520; /*hco */
    W += c[8]*44.009950; /*co2 */
    W += c[9]*28.010550; /*co */
    W += c[10]*31.998800; /*o2 */
    W += c[11]*15.999400; /*o */
    W += c[12]*17.007370; /*oh */
    W += c[13]*33.006770; /*ho2 */
    W += c[14]*34.014740; /*h2o2 */
    W += c[15]*18.015340; /*h2o */
    W += c[16]*25.030270; /*c2h */
    W += c[17]*41.029670; /*hcco */
    W += c[18]*26.038240; /*c2h2 */
    W += c[19]*27.046210; /*c2h3 */
    W += c[20]*28.054180; /*c2h4 */
    W += c[21]*29.062150; /*c2h5 */
    W += c[22]*30.070120; /*c2h6 */
    W += c[23]*31.034460; /*ch2oh */
    W += c[24]*31.034460; /*ch3o */
    W += c[25]*42.037640; /*hccoh */
    W += c[26]*39.057360; /*h2ccch */
    W += c[27]*38.049390; /*c3h2 */
    W += c[28]*14.027090; /*ch2s */
    W += c[29]*42.037640; /*ch2co */
    W += c[30]*43.045610; /*ch2hco */
    W += c[31]*43.045610; /*ch3co */
    W += c[32]*44.053580; /*ch3hco */
    W += c[33]*46.069520; /*c2h5oh */
    W += c[34]*45.061550; /*c2h4oh */
    W += c[35]*45.061550; /*ch3choh */
    W += c[36]*45.061550; /*ch3ch2o */
    W += c[37]*28.013400; /*n2 */

    for (id = 0; id < 38; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31451e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double * P, double * T, double * x, int * iwrk, double * rwrk, double * rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*h2 */
    XW += x[1]*1.007970; /*h */
    XW += x[2]*16.043030; /*ch4 */
    XW += x[3]*15.035060; /*ch3 */
    XW += x[4]*14.027090; /*ch2 */
    XW += x[5]*13.019120; /*ch */
    XW += x[6]*30.026490; /*ch2o */
    XW += x[7]*29.018520; /*hco */
    XW += x[8]*44.009950; /*co2 */
    XW += x[9]*28.010550; /*co */
    XW += x[10]*31.998800; /*o2 */
    XW += x[11]*15.999400; /*o */
    XW += x[12]*17.007370; /*oh */
    XW += x[13]*33.006770; /*ho2 */
    XW += x[14]*34.014740; /*h2o2 */
    XW += x[15]*18.015340; /*h2o */
    XW += x[16]*25.030270; /*c2h */
    XW += x[17]*41.029670; /*hcco */
    XW += x[18]*26.038240; /*c2h2 */
    XW += x[19]*27.046210; /*c2h3 */
    XW += x[20]*28.054180; /*c2h4 */
    XW += x[21]*29.062150; /*c2h5 */
    XW += x[22]*30.070120; /*c2h6 */
    XW += x[23]*31.034460; /*ch2oh */
    XW += x[24]*31.034460; /*ch3o */
    XW += x[25]*42.037640; /*hccoh */
    XW += x[26]*39.057360; /*h2ccch */
    XW += x[27]*38.049390; /*c3h2 */
    XW += x[28]*14.027090; /*ch2s */
    XW += x[29]*42.037640; /*ch2co */
    XW += x[30]*43.045610; /*ch2hco */
    XW += x[31]*43.045610; /*ch3co */
    XW += x[32]*44.053580; /*ch3hco */
    XW += x[33]*46.069520; /*c2h5oh */
    XW += x[34]*45.061550; /*c2h4oh */
    XW += x[35]*45.061550; /*ch3choh */
    XW += x[36]*45.061550; /*ch3ch2o */
    XW += x[37]*28.013400; /*n2 */
    *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double * P, double * T, double * y, int * iwrk, double * rwrk, double * rho)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]/2.015940; /*h2 */
    YOW += y[1]/1.007970; /*h */
    YOW += y[2]/16.043030; /*ch4 */
    YOW += y[3]/15.035060; /*ch3 */
    YOW += y[4]/14.027090; /*ch2 */
    YOW += y[5]/13.019120; /*ch */
    YOW += y[6]/30.026490; /*ch2o */
    YOW += y[7]/29.018520; /*hco */
    YOW += y[8]/44.009950; /*co2 */
    YOW += y[9]/28.010550; /*co */
    YOW += y[10]/31.998800; /*o2 */
    YOW += y[11]/15.999400; /*o */
    YOW += y[12]/17.007370; /*oh */
    YOW += y[13]/33.006770; /*ho2 */
    YOW += y[14]/34.014740; /*h2o2 */
    YOW += y[15]/18.015340; /*h2o */
    YOW += y[16]/25.030270; /*c2h */
    YOW += y[17]/41.029670; /*hcco */
    YOW += y[18]/26.038240; /*c2h2 */
    YOW += y[19]/27.046210; /*c2h3 */
    YOW += y[20]/28.054180; /*c2h4 */
    YOW += y[21]/29.062150; /*c2h5 */
    YOW += y[22]/30.070120; /*c2h6 */
    YOW += y[23]/31.034460; /*ch2oh */
    YOW += y[24]/31.034460; /*ch3o */
    YOW += y[25]/42.037640; /*hccoh */
    YOW += y[26]/39.057360; /*h2ccch */
    YOW += y[27]/38.049390; /*c3h2 */
    YOW += y[28]/14.027090; /*ch2s */
    YOW += y[29]/42.037640; /*ch2co */
    YOW += y[30]/43.045610; /*ch2hco */
    YOW += y[31]/43.045610; /*ch3co */
    YOW += y[32]/44.053580; /*ch3hco */
    YOW += y[33]/46.069520; /*c2h5oh */
    YOW += y[34]/45.061550; /*c2h4oh */
    YOW += y[35]/45.061550; /*ch3choh */
    YOW += y[36]/45.061550; /*ch3ch2o */
    YOW += y[37]/28.013400; /*n2 */
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
    W += c[0]*2.015940; /*h2 */
    W += c[1]*1.007970; /*h */
    W += c[2]*16.043030; /*ch4 */
    W += c[3]*15.035060; /*ch3 */
    W += c[4]*14.027090; /*ch2 */
    W += c[5]*13.019120; /*ch */
    W += c[6]*30.026490; /*ch2o */
    W += c[7]*29.018520; /*hco */
    W += c[8]*44.009950; /*co2 */
    W += c[9]*28.010550; /*co */
    W += c[10]*31.998800; /*o2 */
    W += c[11]*15.999400; /*o */
    W += c[12]*17.007370; /*oh */
    W += c[13]*33.006770; /*ho2 */
    W += c[14]*34.014740; /*h2o2 */
    W += c[15]*18.015340; /*h2o */
    W += c[16]*25.030270; /*c2h */
    W += c[17]*41.029670; /*hcco */
    W += c[18]*26.038240; /*c2h2 */
    W += c[19]*27.046210; /*c2h3 */
    W += c[20]*28.054180; /*c2h4 */
    W += c[21]*29.062150; /*c2h5 */
    W += c[22]*30.070120; /*c2h6 */
    W += c[23]*31.034460; /*ch2oh */
    W += c[24]*31.034460; /*ch3o */
    W += c[25]*42.037640; /*hccoh */
    W += c[26]*39.057360; /*h2ccch */
    W += c[27]*38.049390; /*c3h2 */
    W += c[28]*14.027090; /*ch2s */
    W += c[29]*42.037640; /*ch2co */
    W += c[30]*43.045610; /*ch2hco */
    W += c[31]*43.045610; /*ch3co */
    W += c[32]*44.053580; /*ch3hco */
    W += c[33]*46.069520; /*c2h5oh */
    W += c[34]*45.061550; /*c2h4oh */
    W += c[35]*45.061550; /*ch3choh */
    W += c[36]*45.061550; /*ch3ch2o */
    W += c[37]*28.013400; /*n2 */

    for (id = 0; id < 38; ++id) {
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
    YOW += y[0]/2.015940; /*h2 */
    YOW += y[1]/1.007970; /*h */
    YOW += y[2]/16.043030; /*ch4 */
    YOW += y[3]/15.035060; /*ch3 */
    YOW += y[4]/14.027090; /*ch2 */
    YOW += y[5]/13.019120; /*ch */
    YOW += y[6]/30.026490; /*ch2o */
    YOW += y[7]/29.018520; /*hco */
    YOW += y[8]/44.009950; /*co2 */
    YOW += y[9]/28.010550; /*co */
    YOW += y[10]/31.998800; /*o2 */
    YOW += y[11]/15.999400; /*o */
    YOW += y[12]/17.007370; /*oh */
    YOW += y[13]/33.006770; /*ho2 */
    YOW += y[14]/34.014740; /*h2o2 */
    YOW += y[15]/18.015340; /*h2o */
    YOW += y[16]/25.030270; /*c2h */
    YOW += y[17]/41.029670; /*hcco */
    YOW += y[18]/26.038240; /*c2h2 */
    YOW += y[19]/27.046210; /*c2h3 */
    YOW += y[20]/28.054180; /*c2h4 */
    YOW += y[21]/29.062150; /*c2h5 */
    YOW += y[22]/30.070120; /*c2h6 */
    YOW += y[23]/31.034460; /*ch2oh */
    YOW += y[24]/31.034460; /*ch3o */
    YOW += y[25]/42.037640; /*hccoh */
    YOW += y[26]/39.057360; /*h2ccch */
    YOW += y[27]/38.049390; /*c3h2 */
    YOW += y[28]/14.027090; /*ch2s */
    YOW += y[29]/42.037640; /*ch2co */
    YOW += y[30]/43.045610; /*ch2hco */
    YOW += y[31]/43.045610; /*ch3co */
    YOW += y[32]/44.053580; /*ch3hco */
    YOW += y[33]/46.069520; /*c2h5oh */
    YOW += y[34]/45.061550; /*c2h4oh */
    YOW += y[35]/45.061550; /*ch3choh */
    YOW += y[36]/45.061550; /*ch3ch2o */
    YOW += y[37]/28.013400; /*n2 */
    *wtm = 1.0 / YOW;

    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWX(double *x, int * iwrk, double * rwrk, double * wtm)
{
    double XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*2.015940; /*h2 */
    XW += x[1]*1.007970; /*h */
    XW += x[2]*16.043030; /*ch4 */
    XW += x[3]*15.035060; /*ch3 */
    XW += x[4]*14.027090; /*ch2 */
    XW += x[5]*13.019120; /*ch */
    XW += x[6]*30.026490; /*ch2o */
    XW += x[7]*29.018520; /*hco */
    XW += x[8]*44.009950; /*co2 */
    XW += x[9]*28.010550; /*co */
    XW += x[10]*31.998800; /*o2 */
    XW += x[11]*15.999400; /*o */
    XW += x[12]*17.007370; /*oh */
    XW += x[13]*33.006770; /*ho2 */
    XW += x[14]*34.014740; /*h2o2 */
    XW += x[15]*18.015340; /*h2o */
    XW += x[16]*25.030270; /*c2h */
    XW += x[17]*41.029670; /*hcco */
    XW += x[18]*26.038240; /*c2h2 */
    XW += x[19]*27.046210; /*c2h3 */
    XW += x[20]*28.054180; /*c2h4 */
    XW += x[21]*29.062150; /*c2h5 */
    XW += x[22]*30.070120; /*c2h6 */
    XW += x[23]*31.034460; /*ch2oh */
    XW += x[24]*31.034460; /*ch3o */
    XW += x[25]*42.037640; /*hccoh */
    XW += x[26]*39.057360; /*h2ccch */
    XW += x[27]*38.049390; /*c3h2 */
    XW += x[28]*14.027090; /*ch2s */
    XW += x[29]*42.037640; /*ch2co */
    XW += x[30]*43.045610; /*ch2hco */
    XW += x[31]*43.045610; /*ch3co */
    XW += x[32]*44.053580; /*ch3hco */
    XW += x[33]*46.069520; /*c2h5oh */
    XW += x[34]*45.061550; /*c2h4oh */
    XW += x[35]*45.061550; /*ch3choh */
    XW += x[36]*45.061550; /*ch3ch2o */
    XW += x[37]*28.013400; /*n2 */
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
    W += c[0]*2.015940; /*h2 */
    W += c[1]*1.007970; /*h */
    W += c[2]*16.043030; /*ch4 */
    W += c[3]*15.035060; /*ch3 */
    W += c[4]*14.027090; /*ch2 */
    W += c[5]*13.019120; /*ch */
    W += c[6]*30.026490; /*ch2o */
    W += c[7]*29.018520; /*hco */
    W += c[8]*44.009950; /*co2 */
    W += c[9]*28.010550; /*co */
    W += c[10]*31.998800; /*o2 */
    W += c[11]*15.999400; /*o */
    W += c[12]*17.007370; /*oh */
    W += c[13]*33.006770; /*ho2 */
    W += c[14]*34.014740; /*h2o2 */
    W += c[15]*18.015340; /*h2o */
    W += c[16]*25.030270; /*c2h */
    W += c[17]*41.029670; /*hcco */
    W += c[18]*26.038240; /*c2h2 */
    W += c[19]*27.046210; /*c2h3 */
    W += c[20]*28.054180; /*c2h4 */
    W += c[21]*29.062150; /*c2h5 */
    W += c[22]*30.070120; /*c2h6 */
    W += c[23]*31.034460; /*ch2oh */
    W += c[24]*31.034460; /*ch3o */
    W += c[25]*42.037640; /*hccoh */
    W += c[26]*39.057360; /*h2ccch */
    W += c[27]*38.049390; /*c3h2 */
    W += c[28]*14.027090; /*ch2s */
    W += c[29]*42.037640; /*ch2co */
    W += c[30]*43.045610; /*ch2hco */
    W += c[31]*43.045610; /*ch3co */
    W += c[32]*44.053580; /*ch3hco */
    W += c[33]*46.069520; /*c2h5oh */
    W += c[34]*45.061550; /*c2h4oh */
    W += c[35]*45.061550; /*ch3choh */
    W += c[36]*45.061550; /*ch3ch2o */
    W += c[37]*28.013400; /*n2 */

    for (id = 0; id < 38; ++id) {
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
    YOW += y[0]/2.015940; /*h2 */
    YOW += y[1]/1.007970; /*h */
    YOW += y[2]/16.043030; /*ch4 */
    YOW += y[3]/15.035060; /*ch3 */
    YOW += y[4]/14.027090; /*ch2 */
    YOW += y[5]/13.019120; /*ch */
    YOW += y[6]/30.026490; /*ch2o */
    YOW += y[7]/29.018520; /*hco */
    YOW += y[8]/44.009950; /*co2 */
    YOW += y[9]/28.010550; /*co */
    YOW += y[10]/31.998800; /*o2 */
    YOW += y[11]/15.999400; /*o */
    YOW += y[12]/17.007370; /*oh */
    YOW += y[13]/33.006770; /*ho2 */
    YOW += y[14]/34.014740; /*h2o2 */
    YOW += y[15]/18.015340; /*h2o */
    YOW += y[16]/25.030270; /*c2h */
    YOW += y[17]/41.029670; /*hcco */
    YOW += y[18]/26.038240; /*c2h2 */
    YOW += y[19]/27.046210; /*c2h3 */
    YOW += y[20]/28.054180; /*c2h4 */
    YOW += y[21]/29.062150; /*c2h5 */
    YOW += y[22]/30.070120; /*c2h6 */
    YOW += y[23]/31.034460; /*ch2oh */
    YOW += y[24]/31.034460; /*ch3o */
    YOW += y[25]/42.037640; /*hccoh */
    YOW += y[26]/39.057360; /*h2ccch */
    YOW += y[27]/38.049390; /*c3h2 */
    YOW += y[28]/14.027090; /*ch2s */
    YOW += y[29]/42.037640; /*ch2co */
    YOW += y[30]/43.045610; /*ch2hco */
    YOW += y[31]/43.045610; /*ch3co */
    YOW += y[32]/44.053580; /*ch3hco */
    YOW += y[33]/46.069520; /*c2h5oh */
    YOW += y[34]/45.061550; /*c2h4oh */
    YOW += y[35]/45.061550; /*ch3choh */
    YOW += y[36]/45.061550; /*ch3ch2o */
    YOW += y[37]/28.013400; /*n2 */
    /*Now compute conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(16.043030*YOW); 
    x[3] = y[3]/(15.035060*YOW); 
    x[4] = y[4]/(14.027090*YOW); 
    x[5] = y[5]/(13.019120*YOW); 
    x[6] = y[6]/(30.026490*YOW); 
    x[7] = y[7]/(29.018520*YOW); 
    x[8] = y[8]/(44.009950*YOW); 
    x[9] = y[9]/(28.010550*YOW); 
    x[10] = y[10]/(31.998800*YOW); 
    x[11] = y[11]/(15.999400*YOW); 
    x[12] = y[12]/(17.007370*YOW); 
    x[13] = y[13]/(33.006770*YOW); 
    x[14] = y[14]/(34.014740*YOW); 
    x[15] = y[15]/(18.015340*YOW); 
    x[16] = y[16]/(25.030270*YOW); 
    x[17] = y[17]/(41.029670*YOW); 
    x[18] = y[18]/(26.038240*YOW); 
    x[19] = y[19]/(27.046210*YOW); 
    x[20] = y[20]/(28.054180*YOW); 
    x[21] = y[21]/(29.062150*YOW); 
    x[22] = y[22]/(30.070120*YOW); 
    x[23] = y[23]/(31.034460*YOW); 
    x[24] = y[24]/(31.034460*YOW); 
    x[25] = y[25]/(42.037640*YOW); 
    x[26] = y[26]/(39.057360*YOW); 
    x[27] = y[27]/(38.049390*YOW); 
    x[28] = y[28]/(14.027090*YOW); 
    x[29] = y[29]/(42.037640*YOW); 
    x[30] = y[30]/(43.045610*YOW); 
    x[31] = y[31]/(43.045610*YOW); 
    x[32] = y[32]/(44.053580*YOW); 
    x[33] = y[33]/(46.069520*YOW); 
    x[34] = y[34]/(45.061550*YOW); 
    x[35] = y[35]/(45.061550*YOW); 
    x[36] = y[36]/(45.061550*YOW); 
    x[37] = y[37]/(28.013400*YOW); 

    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * c)
{
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/2.015940; /*h2 */
    YOW += y[1]/1.007970; /*h */
    YOW += y[2]/16.043030; /*ch4 */
    YOW += y[3]/15.035060; /*ch3 */
    YOW += y[4]/14.027090; /*ch2 */
    YOW += y[5]/13.019120; /*ch */
    YOW += y[6]/30.026490; /*ch2o */
    YOW += y[7]/29.018520; /*hco */
    YOW += y[8]/44.009950; /*co2 */
    YOW += y[9]/28.010550; /*co */
    YOW += y[10]/31.998800; /*o2 */
    YOW += y[11]/15.999400; /*o */
    YOW += y[12]/17.007370; /*oh */
    YOW += y[13]/33.006770; /*ho2 */
    YOW += y[14]/34.014740; /*h2o2 */
    YOW += y[15]/18.015340; /*h2o */
    YOW += y[16]/25.030270; /*c2h */
    YOW += y[17]/41.029670; /*hcco */
    YOW += y[18]/26.038240; /*c2h2 */
    YOW += y[19]/27.046210; /*c2h3 */
    YOW += y[20]/28.054180; /*c2h4 */
    YOW += y[21]/29.062150; /*c2h5 */
    YOW += y[22]/30.070120; /*c2h6 */
    YOW += y[23]/31.034460; /*ch2oh */
    YOW += y[24]/31.034460; /*ch3o */
    YOW += y[25]/42.037640; /*hccoh */
    YOW += y[26]/39.057360; /*h2ccch */
    YOW += y[27]/38.049390; /*c3h2 */
    YOW += y[28]/14.027090; /*ch2s */
    YOW += y[29]/42.037640; /*ch2co */
    YOW += y[30]/43.045610; /*ch2hco */
    YOW += y[31]/43.045610; /*ch3co */
    YOW += y[32]/44.053580; /*ch3hco */
    YOW += y[33]/46.069520; /*c2h5oh */
    YOW += y[34]/45.061550; /*c2h4oh */
    YOW += y[35]/45.061550; /*ch3choh */
    YOW += y[36]/45.061550; /*ch3ch2o */
    YOW += y[37]/28.013400; /*n2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*Now compute conversion */
    c[0] = PWORT * y[0]/2.015940; 
    c[1] = PWORT * y[1]/1.007970; 
    c[2] = PWORT * y[2]/16.043030; 
    c[3] = PWORT * y[3]/15.035060; 
    c[4] = PWORT * y[4]/14.027090; 
    c[5] = PWORT * y[5]/13.019120; 
    c[6] = PWORT * y[6]/30.026490; 
    c[7] = PWORT * y[7]/29.018520; 
    c[8] = PWORT * y[8]/44.009950; 
    c[9] = PWORT * y[9]/28.010550; 
    c[10] = PWORT * y[10]/31.998800; 
    c[11] = PWORT * y[11]/15.999400; 
    c[12] = PWORT * y[12]/17.007370; 
    c[13] = PWORT * y[13]/33.006770; 
    c[14] = PWORT * y[14]/34.014740; 
    c[15] = PWORT * y[15]/18.015340; 
    c[16] = PWORT * y[16]/25.030270; 
    c[17] = PWORT * y[17]/41.029670; 
    c[18] = PWORT * y[18]/26.038240; 
    c[19] = PWORT * y[19]/27.046210; 
    c[20] = PWORT * y[20]/28.054180; 
    c[21] = PWORT * y[21]/29.062150; 
    c[22] = PWORT * y[22]/30.070120; 
    c[23] = PWORT * y[23]/31.034460; 
    c[24] = PWORT * y[24]/31.034460; 
    c[25] = PWORT * y[25]/42.037640; 
    c[26] = PWORT * y[26]/39.057360; 
    c[27] = PWORT * y[27]/38.049390; 
    c[28] = PWORT * y[28]/14.027090; 
    c[29] = PWORT * y[29]/42.037640; 
    c[30] = PWORT * y[30]/43.045610; 
    c[31] = PWORT * y[31]/43.045610; 
    c[32] = PWORT * y[32]/44.053580; 
    c[33] = PWORT * y[33]/46.069520; 
    c[34] = PWORT * y[34]/45.061550; 
    c[35] = PWORT * y[35]/45.061550; 
    c[36] = PWORT * y[36]/45.061550; 
    c[37] = PWORT * y[37]/28.013400; 

    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * c)
{
    /*See Eq 8 (Temperature not used) */
    c[0] = (*rho) * y[0]/2.015940; 
    c[1] = (*rho) * y[1]/1.007970; 
    c[2] = (*rho) * y[2]/16.043030; 
    c[3] = (*rho) * y[3]/15.035060; 
    c[4] = (*rho) * y[4]/14.027090; 
    c[5] = (*rho) * y[5]/13.019120; 
    c[6] = (*rho) * y[6]/30.026490; 
    c[7] = (*rho) * y[7]/29.018520; 
    c[8] = (*rho) * y[8]/44.009950; 
    c[9] = (*rho) * y[9]/28.010550; 
    c[10] = (*rho) * y[10]/31.998800; 
    c[11] = (*rho) * y[11]/15.999400; 
    c[12] = (*rho) * y[12]/17.007370; 
    c[13] = (*rho) * y[13]/33.006770; 
    c[14] = (*rho) * y[14]/34.014740; 
    c[15] = (*rho) * y[15]/18.015340; 
    c[16] = (*rho) * y[16]/25.030270; 
    c[17] = (*rho) * y[17]/41.029670; 
    c[18] = (*rho) * y[18]/26.038240; 
    c[19] = (*rho) * y[19]/27.046210; 
    c[20] = (*rho) * y[20]/28.054180; 
    c[21] = (*rho) * y[21]/29.062150; 
    c[22] = (*rho) * y[22]/30.070120; 
    c[23] = (*rho) * y[23]/31.034460; 
    c[24] = (*rho) * y[24]/31.034460; 
    c[25] = (*rho) * y[25]/42.037640; 
    c[26] = (*rho) * y[26]/39.057360; 
    c[27] = (*rho) * y[27]/38.049390; 
    c[28] = (*rho) * y[28]/14.027090; 
    c[29] = (*rho) * y[29]/42.037640; 
    c[30] = (*rho) * y[30]/43.045610; 
    c[31] = (*rho) * y[31]/43.045610; 
    c[32] = (*rho) * y[32]/44.053580; 
    c[33] = (*rho) * y[33]/46.069520; 
    c[34] = (*rho) * y[34]/45.061550; 
    c[35] = (*rho) * y[35]/45.061550; 
    c[36] = (*rho) * y[36]/45.061550; 
    c[37] = (*rho) * y[37]/28.013400; 

    return;
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
void CKXTY(double * x, int * iwrk, double * rwrk, double * y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*h2 */
    XW += x[1]*1.007970; /*h */
    XW += x[2]*16.043030; /*ch4 */
    XW += x[3]*15.035060; /*ch3 */
    XW += x[4]*14.027090; /*ch2 */
    XW += x[5]*13.019120; /*ch */
    XW += x[6]*30.026490; /*ch2o */
    XW += x[7]*29.018520; /*hco */
    XW += x[8]*44.009950; /*co2 */
    XW += x[9]*28.010550; /*co */
    XW += x[10]*31.998800; /*o2 */
    XW += x[11]*15.999400; /*o */
    XW += x[12]*17.007370; /*oh */
    XW += x[13]*33.006770; /*ho2 */
    XW += x[14]*34.014740; /*h2o2 */
    XW += x[15]*18.015340; /*h2o */
    XW += x[16]*25.030270; /*c2h */
    XW += x[17]*41.029670; /*hcco */
    XW += x[18]*26.038240; /*c2h2 */
    XW += x[19]*27.046210; /*c2h3 */
    XW += x[20]*28.054180; /*c2h4 */
    XW += x[21]*29.062150; /*c2h5 */
    XW += x[22]*30.070120; /*c2h6 */
    XW += x[23]*31.034460; /*ch2oh */
    XW += x[24]*31.034460; /*ch3o */
    XW += x[25]*42.037640; /*hccoh */
    XW += x[26]*39.057360; /*h2ccch */
    XW += x[27]*38.049390; /*c3h2 */
    XW += x[28]*14.027090; /*ch2s */
    XW += x[29]*42.037640; /*ch2co */
    XW += x[30]*43.045610; /*ch2hco */
    XW += x[31]*43.045610; /*ch3co */
    XW += x[32]*44.053580; /*ch3hco */
    XW += x[33]*46.069520; /*c2h5oh */
    XW += x[34]*45.061550; /*c2h4oh */
    XW += x[35]*45.061550; /*ch3choh */
    XW += x[36]*45.061550; /*ch3ch2o */
    XW += x[37]*28.013400; /*n2 */
    /*Now compute conversion */
    y[0] = x[0]*2.015940/XW; 
    y[1] = x[1]*1.007970/XW; 
    y[2] = x[2]*16.043030/XW; 
    y[3] = x[3]*15.035060/XW; 
    y[4] = x[4]*14.027090/XW; 
    y[5] = x[5]*13.019120/XW; 
    y[6] = x[6]*30.026490/XW; 
    y[7] = x[7]*29.018520/XW; 
    y[8] = x[8]*44.009950/XW; 
    y[9] = x[9]*28.010550/XW; 
    y[10] = x[10]*31.998800/XW; 
    y[11] = x[11]*15.999400/XW; 
    y[12] = x[12]*17.007370/XW; 
    y[13] = x[13]*33.006770/XW; 
    y[14] = x[14]*34.014740/XW; 
    y[15] = x[15]*18.015340/XW; 
    y[16] = x[16]*25.030270/XW; 
    y[17] = x[17]*41.029670/XW; 
    y[18] = x[18]*26.038240/XW; 
    y[19] = x[19]*27.046210/XW; 
    y[20] = x[20]*28.054180/XW; 
    y[21] = x[21]*29.062150/XW; 
    y[22] = x[22]*30.070120/XW; 
    y[23] = x[23]*31.034460/XW; 
    y[24] = x[24]*31.034460/XW; 
    y[25] = x[25]*42.037640/XW; 
    y[26] = x[26]*39.057360/XW; 
    y[27] = x[27]*38.049390/XW; 
    y[28] = x[28]*14.027090/XW; 
    y[29] = x[29]*42.037640/XW; 
    y[30] = x[30]*43.045610/XW; 
    y[31] = x[31]*43.045610/XW; 
    y[32] = x[32]*44.053580/XW; 
    y[33] = x[33]*46.069520/XW; 
    y[34] = x[34]*45.061550/XW; 
    y[35] = x[35]*45.061550/XW; 
    y[36] = x[36]*45.061550/XW; 
    y[37] = x[37]*28.013400/XW; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 38; ++id) {
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
    XW += x[0]*2.015940; /*h2 */
    XW += x[1]*1.007970; /*h */
    XW += x[2]*16.043030; /*ch4 */
    XW += x[3]*15.035060; /*ch3 */
    XW += x[4]*14.027090; /*ch2 */
    XW += x[5]*13.019120; /*ch */
    XW += x[6]*30.026490; /*ch2o */
    XW += x[7]*29.018520; /*hco */
    XW += x[8]*44.009950; /*co2 */
    XW += x[9]*28.010550; /*co */
    XW += x[10]*31.998800; /*o2 */
    XW += x[11]*15.999400; /*o */
    XW += x[12]*17.007370; /*oh */
    XW += x[13]*33.006770; /*ho2 */
    XW += x[14]*34.014740; /*h2o2 */
    XW += x[15]*18.015340; /*h2o */
    XW += x[16]*25.030270; /*c2h */
    XW += x[17]*41.029670; /*hcco */
    XW += x[18]*26.038240; /*c2h2 */
    XW += x[19]*27.046210; /*c2h3 */
    XW += x[20]*28.054180; /*c2h4 */
    XW += x[21]*29.062150; /*c2h5 */
    XW += x[22]*30.070120; /*c2h6 */
    XW += x[23]*31.034460; /*ch2oh */
    XW += x[24]*31.034460; /*ch3o */
    XW += x[25]*42.037640; /*hccoh */
    XW += x[26]*39.057360; /*h2ccch */
    XW += x[27]*38.049390; /*c3h2 */
    XW += x[28]*14.027090; /*ch2s */
    XW += x[29]*42.037640; /*ch2co */
    XW += x[30]*43.045610; /*ch2hco */
    XW += x[31]*43.045610; /*ch3co */
    XW += x[32]*44.053580; /*ch3hco */
    XW += x[33]*46.069520; /*c2h5oh */
    XW += x[34]*45.061550; /*c2h4oh */
    XW += x[35]*45.061550; /*ch3choh */
    XW += x[36]*45.061550; /*ch3ch2o */
    XW += x[37]*28.013400; /*n2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 38; ++id) {
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
    for (id = 0; id < 38; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    for (id = 0; id < 38; ++id) {
        x[id] = c[id]/sumC;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
void CKCTY(double * c, int * iwrk, double * rwrk, double * y)
{
    double CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*2.015940; /*h2 */
    CW += c[1]*1.007970; /*h */
    CW += c[2]*16.043030; /*ch4 */
    CW += c[3]*15.035060; /*ch3 */
    CW += c[4]*14.027090; /*ch2 */
    CW += c[5]*13.019120; /*ch */
    CW += c[6]*30.026490; /*ch2o */
    CW += c[7]*29.018520; /*hco */
    CW += c[8]*44.009950; /*co2 */
    CW += c[9]*28.010550; /*co */
    CW += c[10]*31.998800; /*o2 */
    CW += c[11]*15.999400; /*o */
    CW += c[12]*17.007370; /*oh */
    CW += c[13]*33.006770; /*ho2 */
    CW += c[14]*34.014740; /*h2o2 */
    CW += c[15]*18.015340; /*h2o */
    CW += c[16]*25.030270; /*c2h */
    CW += c[17]*41.029670; /*hcco */
    CW += c[18]*26.038240; /*c2h2 */
    CW += c[19]*27.046210; /*c2h3 */
    CW += c[20]*28.054180; /*c2h4 */
    CW += c[21]*29.062150; /*c2h5 */
    CW += c[22]*30.070120; /*c2h6 */
    CW += c[23]*31.034460; /*ch2oh */
    CW += c[24]*31.034460; /*ch3o */
    CW += c[25]*42.037640; /*hccoh */
    CW += c[26]*39.057360; /*h2ccch */
    CW += c[27]*38.049390; /*c3h2 */
    CW += c[28]*14.027090; /*ch2s */
    CW += c[29]*42.037640; /*ch2co */
    CW += c[30]*43.045610; /*ch2hco */
    CW += c[31]*43.045610; /*ch3co */
    CW += c[32]*44.053580; /*ch3hco */
    CW += c[33]*46.069520; /*c2h5oh */
    CW += c[34]*45.061550; /*c2h4oh */
    CW += c[35]*45.061550; /*ch3choh */
    CW += c[36]*45.061550; /*ch3ch2o */
    CW += c[37]*28.013400; /*n2 */
    /*Now compute conversion */
    y[0] = c[0]*2.015940/CW; 
    y[1] = c[1]*1.007970/CW; 
    y[2] = c[2]*16.043030/CW; 
    y[3] = c[3]*15.035060/CW; 
    y[4] = c[4]*14.027090/CW; 
    y[5] = c[5]*13.019120/CW; 
    y[6] = c[6]*30.026490/CW; 
    y[7] = c[7]*29.018520/CW; 
    y[8] = c[8]*44.009950/CW; 
    y[9] = c[9]*28.010550/CW; 
    y[10] = c[10]*31.998800/CW; 
    y[11] = c[11]*15.999400/CW; 
    y[12] = c[12]*17.007370/CW; 
    y[13] = c[13]*33.006770/CW; 
    y[14] = c[14]*34.014740/CW; 
    y[15] = c[15]*18.015340/CW; 
    y[16] = c[16]*25.030270/CW; 
    y[17] = c[17]*41.029670/CW; 
    y[18] = c[18]*26.038240/CW; 
    y[19] = c[19]*27.046210/CW; 
    y[20] = c[20]*28.054180/CW; 
    y[21] = c[21]*29.062150/CW; 
    y[22] = c[22]*30.070120/CW; 
    y[23] = c[23]*31.034460/CW; 
    y[24] = c[24]*31.034460/CW; 
    y[25] = c[25]*42.037640/CW; 
    y[26] = c[26]*39.057360/CW; 
    y[27] = c[27]*38.049390/CW; 
    y[28] = c[28]*14.027090/CW; 
    y[29] = c[29]*42.037640/CW; 
    y[30] = c[30]*43.045610/CW; 
    y[31] = c[31]*43.045610/CW; 
    y[32] = c[32]*44.053580/CW; 
    y[33] = c[33]*46.069520/CW; 
    y[34] = c[34]*45.061550/CW; 
    y[35] = c[35]*45.061550/CW; 
    y[36] = c[36]*45.061550/CW; 
    y[37] = c[37]*28.013400/CW; 

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
    for (id = 0; id < 38; ++id) {
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
    for (id = 0; id < 38; ++id) {
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
    for (id = 0; id < 38; ++id) {
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
    for (id = 0; id < 38; ++id) {
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
    for (id = 0; id < 38; ++id) {
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
    for (id = 0; id < 38; ++id) {
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
    for (id = 0; id < 38; ++id) {
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
    cvms[0] *= 4.124383662212169e+07; /*h2 */
    cvms[1] *= 8.248767324424338e+07; /*h */
    cvms[2] *= 5.182630712527496e+06; /*ch4 */
    cvms[3] *= 5.530081023953346e+06; /*ch3 */
    cvms[4] *= 5.927466067445207e+06; /*ch2 */
    cvms[5] *= 6.386384025955671e+06; /*ch */
    cvms[6] *= 2.769058254894261e+06; /*ch2o */
    cvms[7] *= 2.865242610581105e+06; /*hco */
    cvms[8] *= 1.889234139098090e+06; /*co2 */
    cvms[9] *= 2.968349425484326e+06; /*co */
    cvms[10] *= 2.598381814318037e+06; /*o2 */
    cvms[11] *= 5.196763628636074e+06; /*o */
    cvms[12] *= 4.888768810227566e+06; /*oh */
    cvms[13] *= 2.519031701678171e+06; /*ho2 */
    cvms[14] *= 2.444384405113783e+06; /*h2o2 */
    cvms[15] *= 4.615239012974499e+06; /*h2o */
    cvms[16] *= 3.321781986370902e+06; /*c2h */
    cvms[17] *= 2.026462801187531e+06; /*hcco */
    cvms[18] *= 3.193192012977835e+06; /*c2h2 */
    cvms[19] *= 3.074186734481467e+06; /*c2h3 */
    cvms[20] *= 2.963733033722604e+06; /*c2h4 */
    cvms[21] *= 2.860941121011349e+06; /*c2h5 */
    cvms[22] *= 2.765040511976673e+06; /*c2h6 */
    cvms[23] *= 2.679121853578248e+06; /*ch2oh */
    cvms[24] *= 2.679121853578248e+06; /*ch3o */
    cvms[25] *= 1.977872687429646e+06; /*hccoh */
    cvms[26] *= 2.128794675318557e+06; /*h2ccch */
    cvms[27] *= 2.185188777007989e+06; /*c3h2 */
    cvms[28] *= 5.927466067445207e+06; /*ch2s */
    cvms[29] *= 1.977872687429646e+06; /*ch2co */
    cvms[30] *= 1.931558177477332e+06; /*ch2hco */
    cvms[31] *= 1.931558177477331e+06; /*ch3co */
    cvms[32] *= 1.887363070152301e+06; /*ch3hco */
    cvms[33] *= 1.804774610197805e+06; /*c2h5oh */
    cvms[34] *= 1.845145140369117e+06; /*c2h4oh */
    cvms[35] *= 1.845145140369117e+06; /*ch3choh */
    cvms[36] *= 1.845145140369117e+06; /*ch3ch2o */
    cvms[37] *= 2.968047434442088e+06; /*n2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
void CKCPMS(double *T, int * iwrk, double * rwrk, double * cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 4.124383662212169e+07; /*h2 */
    cpms[1] *= 8.248767324424338e+07; /*h */
    cpms[2] *= 5.182630712527496e+06; /*ch4 */
    cpms[3] *= 5.530081023953346e+06; /*ch3 */
    cpms[4] *= 5.927466067445207e+06; /*ch2 */
    cpms[5] *= 6.386384025955671e+06; /*ch */
    cpms[6] *= 2.769058254894261e+06; /*ch2o */
    cpms[7] *= 2.865242610581105e+06; /*hco */
    cpms[8] *= 1.889234139098090e+06; /*co2 */
    cpms[9] *= 2.968349425484326e+06; /*co */
    cpms[10] *= 2.598381814318037e+06; /*o2 */
    cpms[11] *= 5.196763628636074e+06; /*o */
    cpms[12] *= 4.888768810227566e+06; /*oh */
    cpms[13] *= 2.519031701678171e+06; /*ho2 */
    cpms[14] *= 2.444384405113783e+06; /*h2o2 */
    cpms[15] *= 4.615239012974499e+06; /*h2o */
    cpms[16] *= 3.321781986370902e+06; /*c2h */
    cpms[17] *= 2.026462801187531e+06; /*hcco */
    cpms[18] *= 3.193192012977835e+06; /*c2h2 */
    cpms[19] *= 3.074186734481467e+06; /*c2h3 */
    cpms[20] *= 2.963733033722604e+06; /*c2h4 */
    cpms[21] *= 2.860941121011349e+06; /*c2h5 */
    cpms[22] *= 2.765040511976673e+06; /*c2h6 */
    cpms[23] *= 2.679121853578248e+06; /*ch2oh */
    cpms[24] *= 2.679121853578248e+06; /*ch3o */
    cpms[25] *= 1.977872687429646e+06; /*hccoh */
    cpms[26] *= 2.128794675318557e+06; /*h2ccch */
    cpms[27] *= 2.185188777007989e+06; /*c3h2 */
    cpms[28] *= 5.927466067445207e+06; /*ch2s */
    cpms[29] *= 1.977872687429646e+06; /*ch2co */
    cpms[30] *= 1.931558177477332e+06; /*ch2hco */
    cpms[31] *= 1.931558177477331e+06; /*ch3co */
    cpms[32] *= 1.887363070152301e+06; /*ch3hco */
    cpms[33] *= 1.804774610197805e+06; /*c2h5oh */
    cpms[34] *= 1.845145140369117e+06; /*c2h4oh */
    cpms[35] *= 1.845145140369117e+06; /*ch3choh */
    cpms[36] *= 1.845145140369117e+06; /*ch3ch2o */
    cpms[37] *= 2.968047434442088e+06; /*n2 */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double *T, int * iwrk, double * rwrk, double * ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    ums[0] *= RT/2.015940; /*h2 */
    ums[1] *= RT/1.007970; /*h */
    ums[2] *= RT/16.043030; /*ch4 */
    ums[3] *= RT/15.035060; /*ch3 */
    ums[4] *= RT/14.027090; /*ch2 */
    ums[5] *= RT/13.019120; /*ch */
    ums[6] *= RT/30.026490; /*ch2o */
    ums[7] *= RT/29.018520; /*hco */
    ums[8] *= RT/44.009950; /*co2 */
    ums[9] *= RT/28.010550; /*co */
    ums[10] *= RT/31.998800; /*o2 */
    ums[11] *= RT/15.999400; /*o */
    ums[12] *= RT/17.007370; /*oh */
    ums[13] *= RT/33.006770; /*ho2 */
    ums[14] *= RT/34.014740; /*h2o2 */
    ums[15] *= RT/18.015340; /*h2o */
    ums[16] *= RT/25.030270; /*c2h */
    ums[17] *= RT/41.029670; /*hcco */
    ums[18] *= RT/26.038240; /*c2h2 */
    ums[19] *= RT/27.046210; /*c2h3 */
    ums[20] *= RT/28.054180; /*c2h4 */
    ums[21] *= RT/29.062150; /*c2h5 */
    ums[22] *= RT/30.070120; /*c2h6 */
    ums[23] *= RT/31.034460; /*ch2oh */
    ums[24] *= RT/31.034460; /*ch3o */
    ums[25] *= RT/42.037640; /*hccoh */
    ums[26] *= RT/39.057360; /*h2ccch */
    ums[27] *= RT/38.049390; /*c3h2 */
    ums[28] *= RT/14.027090; /*ch2s */
    ums[29] *= RT/42.037640; /*ch2co */
    ums[30] *= RT/43.045610; /*ch2hco */
    ums[31] *= RT/43.045610; /*ch3co */
    ums[32] *= RT/44.053580; /*ch3hco */
    ums[33] *= RT/46.069520; /*c2h5oh */
    ums[34] *= RT/45.061550; /*c2h4oh */
    ums[35] *= RT/45.061550; /*ch3choh */
    ums[36] *= RT/45.061550; /*ch3ch2o */
    ums[37] *= RT/28.013400; /*n2 */
}


/*Returns enthalpy in mass units (Eq 27.) */
void CKHMS(double *T, int * iwrk, double * rwrk, double * hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hms, tc);
    hms[0] *= RT/2.015940; /*h2 */
    hms[1] *= RT/1.007970; /*h */
    hms[2] *= RT/16.043030; /*ch4 */
    hms[3] *= RT/15.035060; /*ch3 */
    hms[4] *= RT/14.027090; /*ch2 */
    hms[5] *= RT/13.019120; /*ch */
    hms[6] *= RT/30.026490; /*ch2o */
    hms[7] *= RT/29.018520; /*hco */
    hms[8] *= RT/44.009950; /*co2 */
    hms[9] *= RT/28.010550; /*co */
    hms[10] *= RT/31.998800; /*o2 */
    hms[11] *= RT/15.999400; /*o */
    hms[12] *= RT/17.007370; /*oh */
    hms[13] *= RT/33.006770; /*ho2 */
    hms[14] *= RT/34.014740; /*h2o2 */
    hms[15] *= RT/18.015340; /*h2o */
    hms[16] *= RT/25.030270; /*c2h */
    hms[17] *= RT/41.029670; /*hcco */
    hms[18] *= RT/26.038240; /*c2h2 */
    hms[19] *= RT/27.046210; /*c2h3 */
    hms[20] *= RT/28.054180; /*c2h4 */
    hms[21] *= RT/29.062150; /*c2h5 */
    hms[22] *= RT/30.070120; /*c2h6 */
    hms[23] *= RT/31.034460; /*ch2oh */
    hms[24] *= RT/31.034460; /*ch3o */
    hms[25] *= RT/42.037640; /*hccoh */
    hms[26] *= RT/39.057360; /*h2ccch */
    hms[27] *= RT/38.049390; /*c3h2 */
    hms[28] *= RT/14.027090; /*ch2s */
    hms[29] *= RT/42.037640; /*ch2co */
    hms[30] *= RT/43.045610; /*ch2hco */
    hms[31] *= RT/43.045610; /*ch3co */
    hms[32] *= RT/44.053580; /*ch3hco */
    hms[33] *= RT/46.069520; /*c2h5oh */
    hms[34] *= RT/45.061550; /*c2h4oh */
    hms[35] *= RT/45.061550; /*ch3choh */
    hms[36] *= RT/45.061550; /*ch3ch2o */
    hms[37] *= RT/28.013400; /*n2 */
}


/*Returns gibbs in mass units (Eq 31.) */
void CKGMS(double *T, int * iwrk, double * rwrk, double * gms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gms, tc);
    gms[0] *= RT/2.015940; /*h2 */
    gms[1] *= RT/1.007970; /*h */
    gms[2] *= RT/16.043030; /*ch4 */
    gms[3] *= RT/15.035060; /*ch3 */
    gms[4] *= RT/14.027090; /*ch2 */
    gms[5] *= RT/13.019120; /*ch */
    gms[6] *= RT/30.026490; /*ch2o */
    gms[7] *= RT/29.018520; /*hco */
    gms[8] *= RT/44.009950; /*co2 */
    gms[9] *= RT/28.010550; /*co */
    gms[10] *= RT/31.998800; /*o2 */
    gms[11] *= RT/15.999400; /*o */
    gms[12] *= RT/17.007370; /*oh */
    gms[13] *= RT/33.006770; /*ho2 */
    gms[14] *= RT/34.014740; /*h2o2 */
    gms[15] *= RT/18.015340; /*h2o */
    gms[16] *= RT/25.030270; /*c2h */
    gms[17] *= RT/41.029670; /*hcco */
    gms[18] *= RT/26.038240; /*c2h2 */
    gms[19] *= RT/27.046210; /*c2h3 */
    gms[20] *= RT/28.054180; /*c2h4 */
    gms[21] *= RT/29.062150; /*c2h5 */
    gms[22] *= RT/30.070120; /*c2h6 */
    gms[23] *= RT/31.034460; /*ch2oh */
    gms[24] *= RT/31.034460; /*ch3o */
    gms[25] *= RT/42.037640; /*hccoh */
    gms[26] *= RT/39.057360; /*h2ccch */
    gms[27] *= RT/38.049390; /*c3h2 */
    gms[28] *= RT/14.027090; /*ch2s */
    gms[29] *= RT/42.037640; /*ch2co */
    gms[30] *= RT/43.045610; /*ch2hco */
    gms[31] *= RT/43.045610; /*ch3co */
    gms[32] *= RT/44.053580; /*ch3hco */
    gms[33] *= RT/46.069520; /*c2h5oh */
    gms[34] *= RT/45.061550; /*c2h4oh */
    gms[35] *= RT/45.061550; /*ch3choh */
    gms[36] *= RT/45.061550; /*ch3ch2o */
    gms[37] *= RT/28.013400; /*n2 */
}


/*Returns helmholtz in mass units (Eq 32.) */
void CKAMS(double *T, int * iwrk, double * rwrk, double * ams)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(ams, tc);
    ams[0] *= RT/2.015940; /*h2 */
    ams[1] *= RT/1.007970; /*h */
    ams[2] *= RT/16.043030; /*ch4 */
    ams[3] *= RT/15.035060; /*ch3 */
    ams[4] *= RT/14.027090; /*ch2 */
    ams[5] *= RT/13.019120; /*ch */
    ams[6] *= RT/30.026490; /*ch2o */
    ams[7] *= RT/29.018520; /*hco */
    ams[8] *= RT/44.009950; /*co2 */
    ams[9] *= RT/28.010550; /*co */
    ams[10] *= RT/31.998800; /*o2 */
    ams[11] *= RT/15.999400; /*o */
    ams[12] *= RT/17.007370; /*oh */
    ams[13] *= RT/33.006770; /*ho2 */
    ams[14] *= RT/34.014740; /*h2o2 */
    ams[15] *= RT/18.015340; /*h2o */
    ams[16] *= RT/25.030270; /*c2h */
    ams[17] *= RT/41.029670; /*hcco */
    ams[18] *= RT/26.038240; /*c2h2 */
    ams[19] *= RT/27.046210; /*c2h3 */
    ams[20] *= RT/28.054180; /*c2h4 */
    ams[21] *= RT/29.062150; /*c2h5 */
    ams[22] *= RT/30.070120; /*c2h6 */
    ams[23] *= RT/31.034460; /*ch2oh */
    ams[24] *= RT/31.034460; /*ch3o */
    ams[25] *= RT/42.037640; /*hccoh */
    ams[26] *= RT/39.057360; /*h2ccch */
    ams[27] *= RT/38.049390; /*c3h2 */
    ams[28] *= RT/14.027090; /*ch2s */
    ams[29] *= RT/42.037640; /*ch2co */
    ams[30] *= RT/43.045610; /*ch2hco */
    ams[31] *= RT/43.045610; /*ch3co */
    ams[32] *= RT/44.053580; /*ch3hco */
    ams[33] *= RT/46.069520; /*c2h5oh */
    ams[34] *= RT/45.061550; /*c2h4oh */
    ams[35] *= RT/45.061550; /*ch3choh */
    ams[36] *= RT/45.061550; /*ch3ch2o */
    ams[37] *= RT/28.013400; /*n2 */
}


/*Returns the entropies in mass units (Eq 28.) */
void CKSMS(double *T, int * iwrk, double * rwrk, double * sms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 4.124383662212169e+07; /*h2 */
    sms[1] *= 8.248767324424338e+07; /*h */
    sms[2] *= 5.182630712527496e+06; /*ch4 */
    sms[3] *= 5.530081023953346e+06; /*ch3 */
    sms[4] *= 5.927466067445207e+06; /*ch2 */
    sms[5] *= 6.386384025955671e+06; /*ch */
    sms[6] *= 2.769058254894261e+06; /*ch2o */
    sms[7] *= 2.865242610581105e+06; /*hco */
    sms[8] *= 1.889234139098090e+06; /*co2 */
    sms[9] *= 2.968349425484326e+06; /*co */
    sms[10] *= 2.598381814318037e+06; /*o2 */
    sms[11] *= 5.196763628636074e+06; /*o */
    sms[12] *= 4.888768810227566e+06; /*oh */
    sms[13] *= 2.519031701678171e+06; /*ho2 */
    sms[14] *= 2.444384405113783e+06; /*h2o2 */
    sms[15] *= 4.615239012974499e+06; /*h2o */
    sms[16] *= 3.321781986370902e+06; /*c2h */
    sms[17] *= 2.026462801187531e+06; /*hcco */
    sms[18] *= 3.193192012977835e+06; /*c2h2 */
    sms[19] *= 3.074186734481467e+06; /*c2h3 */
    sms[20] *= 2.963733033722604e+06; /*c2h4 */
    sms[21] *= 2.860941121011349e+06; /*c2h5 */
    sms[22] *= 2.765040511976673e+06; /*c2h6 */
    sms[23] *= 2.679121853578248e+06; /*ch2oh */
    sms[24] *= 2.679121853578248e+06; /*ch3o */
    sms[25] *= 1.977872687429646e+06; /*hccoh */
    sms[26] *= 2.128794675318557e+06; /*h2ccch */
    sms[27] *= 2.185188777007989e+06; /*c3h2 */
    sms[28] *= 5.927466067445207e+06; /*ch2s */
    sms[29] *= 1.977872687429646e+06; /*ch2co */
    sms[30] *= 1.931558177477332e+06; /*ch2hco */
    sms[31] *= 1.931558177477331e+06; /*ch3co */
    sms[32] *= 1.887363070152301e+06; /*ch3hco */
    sms[33] *= 1.804774610197805e+06; /*c2h5oh */
    sms[34] *= 1.845145140369117e+06; /*c2h4oh */
    sms[35] *= 1.845145140369117e+06; /*ch3choh */
    sms[36] *= 1.845145140369117e+06; /*ch3ch2o */
    sms[37] *= 2.968047434442088e+06; /*n2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double *T, double *x, int * iwrk, double * rwrk, double * cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[38]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 38; ++id) {
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
    double cpor[38]; /* temporary storage */
    cp_R(cpor, tc);
    /*multiply by y/molecularweight */
    result += cpor[0]*y[0]/2.015940; /*h2 */
    result += cpor[1]*y[1]/1.007970; /*h */
    result += cpor[2]*y[2]/16.043030; /*ch4 */
    result += cpor[3]*y[3]/15.035060; /*ch3 */
    result += cpor[4]*y[4]/14.027090; /*ch2 */
    result += cpor[5]*y[5]/13.019120; /*ch */
    result += cpor[6]*y[6]/30.026490; /*ch2o */
    result += cpor[7]*y[7]/29.018520; /*hco */
    result += cpor[8]*y[8]/44.009950; /*co2 */
    result += cpor[9]*y[9]/28.010550; /*co */
    result += cpor[10]*y[10]/31.998800; /*o2 */
    result += cpor[11]*y[11]/15.999400; /*o */
    result += cpor[12]*y[12]/17.007370; /*oh */
    result += cpor[13]*y[13]/33.006770; /*ho2 */
    result += cpor[14]*y[14]/34.014740; /*h2o2 */
    result += cpor[15]*y[15]/18.015340; /*h2o */
    result += cpor[16]*y[16]/25.030270; /*c2h */
    result += cpor[17]*y[17]/41.029670; /*hcco */
    result += cpor[18]*y[18]/26.038240; /*c2h2 */
    result += cpor[19]*y[19]/27.046210; /*c2h3 */
    result += cpor[20]*y[20]/28.054180; /*c2h4 */
    result += cpor[21]*y[21]/29.062150; /*c2h5 */
    result += cpor[22]*y[22]/30.070120; /*c2h6 */
    result += cpor[23]*y[23]/31.034460; /*ch2oh */
    result += cpor[24]*y[24]/31.034460; /*ch3o */
    result += cpor[25]*y[25]/42.037640; /*hccoh */
    result += cpor[26]*y[26]/39.057360; /*h2ccch */
    result += cpor[27]*y[27]/38.049390; /*c3h2 */
    result += cpor[28]*y[28]/14.027090; /*ch2s */
    result += cpor[29]*y[29]/42.037640; /*ch2co */
    result += cpor[30]*y[30]/43.045610; /*ch2hco */
    result += cpor[31]*y[31]/43.045610; /*ch3co */
    result += cpor[32]*y[32]/44.053580; /*ch3hco */
    result += cpor[33]*y[33]/46.069520; /*c2h5oh */
    result += cpor[34]*y[34]/45.061550; /*c2h4oh */
    result += cpor[35]*y[35]/45.061550; /*ch3choh */
    result += cpor[36]*y[36]/45.061550; /*ch3ch2o */
    result += cpor[37]*y[37]/28.013400; /*n2 */

    *cpbs = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double *T, double *x, int * iwrk, double * rwrk, double * cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[38]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 38; ++id) {
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
    double cvor[38]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]/2.015940; /*h2 */
    result += cvor[1]*y[1]/1.007970; /*h */
    result += cvor[2]*y[2]/16.043030; /*ch4 */
    result += cvor[3]*y[3]/15.035060; /*ch3 */
    result += cvor[4]*y[4]/14.027090; /*ch2 */
    result += cvor[5]*y[5]/13.019120; /*ch */
    result += cvor[6]*y[6]/30.026490; /*ch2o */
    result += cvor[7]*y[7]/29.018520; /*hco */
    result += cvor[8]*y[8]/44.009950; /*co2 */
    result += cvor[9]*y[9]/28.010550; /*co */
    result += cvor[10]*y[10]/31.998800; /*o2 */
    result += cvor[11]*y[11]/15.999400; /*o */
    result += cvor[12]*y[12]/17.007370; /*oh */
    result += cvor[13]*y[13]/33.006770; /*ho2 */
    result += cvor[14]*y[14]/34.014740; /*h2o2 */
    result += cvor[15]*y[15]/18.015340; /*h2o */
    result += cvor[16]*y[16]/25.030270; /*c2h */
    result += cvor[17]*y[17]/41.029670; /*hcco */
    result += cvor[18]*y[18]/26.038240; /*c2h2 */
    result += cvor[19]*y[19]/27.046210; /*c2h3 */
    result += cvor[20]*y[20]/28.054180; /*c2h4 */
    result += cvor[21]*y[21]/29.062150; /*c2h5 */
    result += cvor[22]*y[22]/30.070120; /*c2h6 */
    result += cvor[23]*y[23]/31.034460; /*ch2oh */
    result += cvor[24]*y[24]/31.034460; /*ch3o */
    result += cvor[25]*y[25]/42.037640; /*hccoh */
    result += cvor[26]*y[26]/39.057360; /*h2ccch */
    result += cvor[27]*y[27]/38.049390; /*c3h2 */
    result += cvor[28]*y[28]/14.027090; /*ch2s */
    result += cvor[29]*y[29]/42.037640; /*ch2co */
    result += cvor[30]*y[30]/43.045610; /*ch2hco */
    result += cvor[31]*y[31]/43.045610; /*ch3co */
    result += cvor[32]*y[32]/44.053580; /*ch3hco */
    result += cvor[33]*y[33]/46.069520; /*c2h5oh */
    result += cvor[34]*y[34]/45.061550; /*c2h4oh */
    result += cvor[35]*y[35]/45.061550; /*ch3choh */
    result += cvor[36]*y[36]/45.061550; /*ch3ch2o */
    result += cvor[37]*y[37]/28.013400; /*n2 */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double *T, double *x, int * iwrk, double * rwrk, double * hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[38]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 38; ++id) {
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
    double hml[38]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*hml[0]/2.015940; /*h2 */
    result += y[1]*hml[1]/1.007970; /*h */
    result += y[2]*hml[2]/16.043030; /*ch4 */
    result += y[3]*hml[3]/15.035060; /*ch3 */
    result += y[4]*hml[4]/14.027090; /*ch2 */
    result += y[5]*hml[5]/13.019120; /*ch */
    result += y[6]*hml[6]/30.026490; /*ch2o */
    result += y[7]*hml[7]/29.018520; /*hco */
    result += y[8]*hml[8]/44.009950; /*co2 */
    result += y[9]*hml[9]/28.010550; /*co */
    result += y[10]*hml[10]/31.998800; /*o2 */
    result += y[11]*hml[11]/15.999400; /*o */
    result += y[12]*hml[12]/17.007370; /*oh */
    result += y[13]*hml[13]/33.006770; /*ho2 */
    result += y[14]*hml[14]/34.014740; /*h2o2 */
    result += y[15]*hml[15]/18.015340; /*h2o */
    result += y[16]*hml[16]/25.030270; /*c2h */
    result += y[17]*hml[17]/41.029670; /*hcco */
    result += y[18]*hml[18]/26.038240; /*c2h2 */
    result += y[19]*hml[19]/27.046210; /*c2h3 */
    result += y[20]*hml[20]/28.054180; /*c2h4 */
    result += y[21]*hml[21]/29.062150; /*c2h5 */
    result += y[22]*hml[22]/30.070120; /*c2h6 */
    result += y[23]*hml[23]/31.034460; /*ch2oh */
    result += y[24]*hml[24]/31.034460; /*ch3o */
    result += y[25]*hml[25]/42.037640; /*hccoh */
    result += y[26]*hml[26]/39.057360; /*h2ccch */
    result += y[27]*hml[27]/38.049390; /*c3h2 */
    result += y[28]*hml[28]/14.027090; /*ch2s */
    result += y[29]*hml[29]/42.037640; /*ch2co */
    result += y[30]*hml[30]/43.045610; /*ch2hco */
    result += y[31]*hml[31]/43.045610; /*ch3co */
    result += y[32]*hml[32]/44.053580; /*ch3hco */
    result += y[33]*hml[33]/46.069520; /*c2h5oh */
    result += y[34]*hml[34]/45.061550; /*c2h4oh */
    result += y[35]*hml[35]/45.061550; /*ch3choh */
    result += y[36]*hml[36]/45.061550; /*ch3ch2o */
    result += y[37]*hml[37]/28.013400; /*n2 */

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double *T, double *x, int * iwrk, double * rwrk, double * ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[38]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 38; ++id) {
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
    double ums[38]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]/2.015940; /*h2 */
    result += y[1]*ums[1]/1.007970; /*h */
    result += y[2]*ums[2]/16.043030; /*ch4 */
    result += y[3]*ums[3]/15.035060; /*ch3 */
    result += y[4]*ums[4]/14.027090; /*ch2 */
    result += y[5]*ums[5]/13.019120; /*ch */
    result += y[6]*ums[6]/30.026490; /*ch2o */
    result += y[7]*ums[7]/29.018520; /*hco */
    result += y[8]*ums[8]/44.009950; /*co2 */
    result += y[9]*ums[9]/28.010550; /*co */
    result += y[10]*ums[10]/31.998800; /*o2 */
    result += y[11]*ums[11]/15.999400; /*o */
    result += y[12]*ums[12]/17.007370; /*oh */
    result += y[13]*ums[13]/33.006770; /*ho2 */
    result += y[14]*ums[14]/34.014740; /*h2o2 */
    result += y[15]*ums[15]/18.015340; /*h2o */
    result += y[16]*ums[16]/25.030270; /*c2h */
    result += y[17]*ums[17]/41.029670; /*hcco */
    result += y[18]*ums[18]/26.038240; /*c2h2 */
    result += y[19]*ums[19]/27.046210; /*c2h3 */
    result += y[20]*ums[20]/28.054180; /*c2h4 */
    result += y[21]*ums[21]/29.062150; /*c2h5 */
    result += y[22]*ums[22]/30.070120; /*c2h6 */
    result += y[23]*ums[23]/31.034460; /*ch2oh */
    result += y[24]*ums[24]/31.034460; /*ch3o */
    result += y[25]*ums[25]/42.037640; /*hccoh */
    result += y[26]*ums[26]/39.057360; /*h2ccch */
    result += y[27]*ums[27]/38.049390; /*c3h2 */
    result += y[28]*ums[28]/14.027090; /*ch2s */
    result += y[29]*ums[29]/42.037640; /*ch2co */
    result += y[30]*ums[30]/43.045610; /*ch2hco */
    result += y[31]*ums[31]/43.045610; /*ch3co */
    result += y[32]*ums[32]/44.053580; /*ch3hco */
    result += y[33]*ums[33]/46.069520; /*c2h5oh */
    result += y[34]*ums[34]/45.061550; /*c2h4oh */
    result += y[35]*ums[35]/45.061550; /*ch3choh */
    result += y[36]*ums[36]/45.061550; /*ch3ch2o */
    result += y[37]*ums[37]/28.013400; /*n2 */

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
    double sor[38]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 38; ++id) {
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
    double sor[38]; /* temporary storage */
    double x[38]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/2.015940; /*h2 */
    YOW += y[1]/1.007970; /*h */
    YOW += y[2]/16.043030; /*ch4 */
    YOW += y[3]/15.035060; /*ch3 */
    YOW += y[4]/14.027090; /*ch2 */
    YOW += y[5]/13.019120; /*ch */
    YOW += y[6]/30.026490; /*ch2o */
    YOW += y[7]/29.018520; /*hco */
    YOW += y[8]/44.009950; /*co2 */
    YOW += y[9]/28.010550; /*co */
    YOW += y[10]/31.998800; /*o2 */
    YOW += y[11]/15.999400; /*o */
    YOW += y[12]/17.007370; /*oh */
    YOW += y[13]/33.006770; /*ho2 */
    YOW += y[14]/34.014740; /*h2o2 */
    YOW += y[15]/18.015340; /*h2o */
    YOW += y[16]/25.030270; /*c2h */
    YOW += y[17]/41.029670; /*hcco */
    YOW += y[18]/26.038240; /*c2h2 */
    YOW += y[19]/27.046210; /*c2h3 */
    YOW += y[20]/28.054180; /*c2h4 */
    YOW += y[21]/29.062150; /*c2h5 */
    YOW += y[22]/30.070120; /*c2h6 */
    YOW += y[23]/31.034460; /*ch2oh */
    YOW += y[24]/31.034460; /*ch3o */
    YOW += y[25]/42.037640; /*hccoh */
    YOW += y[26]/39.057360; /*h2ccch */
    YOW += y[27]/38.049390; /*c3h2 */
    YOW += y[28]/14.027090; /*ch2s */
    YOW += y[29]/42.037640; /*ch2co */
    YOW += y[30]/43.045610; /*ch2hco */
    YOW += y[31]/43.045610; /*ch3co */
    YOW += y[32]/44.053580; /*ch3hco */
    YOW += y[33]/46.069520; /*c2h5oh */
    YOW += y[34]/45.061550; /*c2h4oh */
    YOW += y[35]/45.061550; /*ch3choh */
    YOW += y[36]/45.061550; /*ch3ch2o */
    YOW += y[37]/28.013400; /*n2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(16.043030*YOW); 
    x[3] = y[3]/(15.035060*YOW); 
    x[4] = y[4]/(14.027090*YOW); 
    x[5] = y[5]/(13.019120*YOW); 
    x[6] = y[6]/(30.026490*YOW); 
    x[7] = y[7]/(29.018520*YOW); 
    x[8] = y[8]/(44.009950*YOW); 
    x[9] = y[9]/(28.010550*YOW); 
    x[10] = y[10]/(31.998800*YOW); 
    x[11] = y[11]/(15.999400*YOW); 
    x[12] = y[12]/(17.007370*YOW); 
    x[13] = y[13]/(33.006770*YOW); 
    x[14] = y[14]/(34.014740*YOW); 
    x[15] = y[15]/(18.015340*YOW); 
    x[16] = y[16]/(25.030270*YOW); 
    x[17] = y[17]/(41.029670*YOW); 
    x[18] = y[18]/(26.038240*YOW); 
    x[19] = y[19]/(27.046210*YOW); 
    x[20] = y[20]/(28.054180*YOW); 
    x[21] = y[21]/(29.062150*YOW); 
    x[22] = y[22]/(30.070120*YOW); 
    x[23] = y[23]/(31.034460*YOW); 
    x[24] = y[24]/(31.034460*YOW); 
    x[25] = y[25]/(42.037640*YOW); 
    x[26] = y[26]/(39.057360*YOW); 
    x[27] = y[27]/(38.049390*YOW); 
    x[28] = y[28]/(14.027090*YOW); 
    x[29] = y[29]/(42.037640*YOW); 
    x[30] = y[30]/(43.045610*YOW); 
    x[31] = y[31]/(43.045610*YOW); 
    x[32] = y[32]/(44.053580*YOW); 
    x[33] = y[33]/(46.069520*YOW); 
    x[34] = y[34]/(45.061550*YOW); 
    x[35] = y[35]/(45.061550*YOW); 
    x[36] = y[36]/(45.061550*YOW); 
    x[37] = y[37]/(28.013400*YOW); 
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
    double gort[38]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 38; ++id) {
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
    double gort[38]; /* temporary storage */
    double x[38]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/2.015940; /*h2 */
    YOW += y[1]/1.007970; /*h */
    YOW += y[2]/16.043030; /*ch4 */
    YOW += y[3]/15.035060; /*ch3 */
    YOW += y[4]/14.027090; /*ch2 */
    YOW += y[5]/13.019120; /*ch */
    YOW += y[6]/30.026490; /*ch2o */
    YOW += y[7]/29.018520; /*hco */
    YOW += y[8]/44.009950; /*co2 */
    YOW += y[9]/28.010550; /*co */
    YOW += y[10]/31.998800; /*o2 */
    YOW += y[11]/15.999400; /*o */
    YOW += y[12]/17.007370; /*oh */
    YOW += y[13]/33.006770; /*ho2 */
    YOW += y[14]/34.014740; /*h2o2 */
    YOW += y[15]/18.015340; /*h2o */
    YOW += y[16]/25.030270; /*c2h */
    YOW += y[17]/41.029670; /*hcco */
    YOW += y[18]/26.038240; /*c2h2 */
    YOW += y[19]/27.046210; /*c2h3 */
    YOW += y[20]/28.054180; /*c2h4 */
    YOW += y[21]/29.062150; /*c2h5 */
    YOW += y[22]/30.070120; /*c2h6 */
    YOW += y[23]/31.034460; /*ch2oh */
    YOW += y[24]/31.034460; /*ch3o */
    YOW += y[25]/42.037640; /*hccoh */
    YOW += y[26]/39.057360; /*h2ccch */
    YOW += y[27]/38.049390; /*c3h2 */
    YOW += y[28]/14.027090; /*ch2s */
    YOW += y[29]/42.037640; /*ch2co */
    YOW += y[30]/43.045610; /*ch2hco */
    YOW += y[31]/43.045610; /*ch3co */
    YOW += y[32]/44.053580; /*ch3hco */
    YOW += y[33]/46.069520; /*c2h5oh */
    YOW += y[34]/45.061550; /*c2h4oh */
    YOW += y[35]/45.061550; /*ch3choh */
    YOW += y[36]/45.061550; /*ch3ch2o */
    YOW += y[37]/28.013400; /*n2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(16.043030*YOW); 
    x[3] = y[3]/(15.035060*YOW); 
    x[4] = y[4]/(14.027090*YOW); 
    x[5] = y[5]/(13.019120*YOW); 
    x[6] = y[6]/(30.026490*YOW); 
    x[7] = y[7]/(29.018520*YOW); 
    x[8] = y[8]/(44.009950*YOW); 
    x[9] = y[9]/(28.010550*YOW); 
    x[10] = y[10]/(31.998800*YOW); 
    x[11] = y[11]/(15.999400*YOW); 
    x[12] = y[12]/(17.007370*YOW); 
    x[13] = y[13]/(33.006770*YOW); 
    x[14] = y[14]/(34.014740*YOW); 
    x[15] = y[15]/(18.015340*YOW); 
    x[16] = y[16]/(25.030270*YOW); 
    x[17] = y[17]/(41.029670*YOW); 
    x[18] = y[18]/(26.038240*YOW); 
    x[19] = y[19]/(27.046210*YOW); 
    x[20] = y[20]/(28.054180*YOW); 
    x[21] = y[21]/(29.062150*YOW); 
    x[22] = y[22]/(30.070120*YOW); 
    x[23] = y[23]/(31.034460*YOW); 
    x[24] = y[24]/(31.034460*YOW); 
    x[25] = y[25]/(42.037640*YOW); 
    x[26] = y[26]/(39.057360*YOW); 
    x[27] = y[27]/(38.049390*YOW); 
    x[28] = y[28]/(14.027090*YOW); 
    x[29] = y[29]/(42.037640*YOW); 
    x[30] = y[30]/(43.045610*YOW); 
    x[31] = y[31]/(43.045610*YOW); 
    x[32] = y[32]/(44.053580*YOW); 
    x[33] = y[33]/(46.069520*YOW); 
    x[34] = y[34]/(45.061550*YOW); 
    x[35] = y[35]/(45.061550*YOW); 
    x[36] = y[36]/(45.061550*YOW); 
    x[37] = y[37]/(28.013400*YOW); 
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
    double aort[38]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 38; ++id) {
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
    double aort[38]; /* temporary storage */
    double x[38]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/2.015940; /*h2 */
    YOW += y[1]/1.007970; /*h */
    YOW += y[2]/16.043030; /*ch4 */
    YOW += y[3]/15.035060; /*ch3 */
    YOW += y[4]/14.027090; /*ch2 */
    YOW += y[5]/13.019120; /*ch */
    YOW += y[6]/30.026490; /*ch2o */
    YOW += y[7]/29.018520; /*hco */
    YOW += y[8]/44.009950; /*co2 */
    YOW += y[9]/28.010550; /*co */
    YOW += y[10]/31.998800; /*o2 */
    YOW += y[11]/15.999400; /*o */
    YOW += y[12]/17.007370; /*oh */
    YOW += y[13]/33.006770; /*ho2 */
    YOW += y[14]/34.014740; /*h2o2 */
    YOW += y[15]/18.015340; /*h2o */
    YOW += y[16]/25.030270; /*c2h */
    YOW += y[17]/41.029670; /*hcco */
    YOW += y[18]/26.038240; /*c2h2 */
    YOW += y[19]/27.046210; /*c2h3 */
    YOW += y[20]/28.054180; /*c2h4 */
    YOW += y[21]/29.062150; /*c2h5 */
    YOW += y[22]/30.070120; /*c2h6 */
    YOW += y[23]/31.034460; /*ch2oh */
    YOW += y[24]/31.034460; /*ch3o */
    YOW += y[25]/42.037640; /*hccoh */
    YOW += y[26]/39.057360; /*h2ccch */
    YOW += y[27]/38.049390; /*c3h2 */
    YOW += y[28]/14.027090; /*ch2s */
    YOW += y[29]/42.037640; /*ch2co */
    YOW += y[30]/43.045610; /*ch2hco */
    YOW += y[31]/43.045610; /*ch3co */
    YOW += y[32]/44.053580; /*ch3hco */
    YOW += y[33]/46.069520; /*c2h5oh */
    YOW += y[34]/45.061550; /*c2h4oh */
    YOW += y[35]/45.061550; /*ch3choh */
    YOW += y[36]/45.061550; /*ch3ch2o */
    YOW += y[37]/28.013400; /*n2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(16.043030*YOW); 
    x[3] = y[3]/(15.035060*YOW); 
    x[4] = y[4]/(14.027090*YOW); 
    x[5] = y[5]/(13.019120*YOW); 
    x[6] = y[6]/(30.026490*YOW); 
    x[7] = y[7]/(29.018520*YOW); 
    x[8] = y[8]/(44.009950*YOW); 
    x[9] = y[9]/(28.010550*YOW); 
    x[10] = y[10]/(31.998800*YOW); 
    x[11] = y[11]/(15.999400*YOW); 
    x[12] = y[12]/(17.007370*YOW); 
    x[13] = y[13]/(33.006770*YOW); 
    x[14] = y[14]/(34.014740*YOW); 
    x[15] = y[15]/(18.015340*YOW); 
    x[16] = y[16]/(25.030270*YOW); 
    x[17] = y[17]/(41.029670*YOW); 
    x[18] = y[18]/(26.038240*YOW); 
    x[19] = y[19]/(27.046210*YOW); 
    x[20] = y[20]/(28.054180*YOW); 
    x[21] = y[21]/(29.062150*YOW); 
    x[22] = y[22]/(30.070120*YOW); 
    x[23] = y[23]/(31.034460*YOW); 
    x[24] = y[24]/(31.034460*YOW); 
    x[25] = y[25]/(42.037640*YOW); 
    x[26] = y[26]/(39.057360*YOW); 
    x[27] = y[27]/(38.049390*YOW); 
    x[28] = y[28]/(14.027090*YOW); 
    x[29] = y[29]/(42.037640*YOW); 
    x[30] = y[30]/(43.045610*YOW); 
    x[31] = y[31]/(43.045610*YOW); 
    x[32] = y[32]/(44.053580*YOW); 
    x[33] = y[33]/(46.069520*YOW); 
    x[34] = y[34]/(45.061550*YOW); 
    x[35] = y[35]/(45.061550*YOW); 
    x[36] = y[36]/(45.061550*YOW); 
    x[37] = y[37]/(28.013400*YOW); 
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
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double * T, double * C, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 38; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 38; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[38]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/2.015940; /*h2 */
    YOW += y[1]/1.007970; /*h */
    YOW += y[2]/16.043030; /*ch4 */
    YOW += y[3]/15.035060; /*ch3 */
    YOW += y[4]/14.027090; /*ch2 */
    YOW += y[5]/13.019120; /*ch */
    YOW += y[6]/30.026490; /*ch2o */
    YOW += y[7]/29.018520; /*hco */
    YOW += y[8]/44.009950; /*co2 */
    YOW += y[9]/28.010550; /*co */
    YOW += y[10]/31.998800; /*o2 */
    YOW += y[11]/15.999400; /*o */
    YOW += y[12]/17.007370; /*oh */
    YOW += y[13]/33.006770; /*ho2 */
    YOW += y[14]/34.014740; /*h2o2 */
    YOW += y[15]/18.015340; /*h2o */
    YOW += y[16]/25.030270; /*c2h */
    YOW += y[17]/41.029670; /*hcco */
    YOW += y[18]/26.038240; /*c2h2 */
    YOW += y[19]/27.046210; /*c2h3 */
    YOW += y[20]/28.054180; /*c2h4 */
    YOW += y[21]/29.062150; /*c2h5 */
    YOW += y[22]/30.070120; /*c2h6 */
    YOW += y[23]/31.034460; /*ch2oh */
    YOW += y[24]/31.034460; /*ch3o */
    YOW += y[25]/42.037640; /*hccoh */
    YOW += y[26]/39.057360; /*h2ccch */
    YOW += y[27]/38.049390; /*c3h2 */
    YOW += y[28]/14.027090; /*ch2s */
    YOW += y[29]/42.037640; /*ch2co */
    YOW += y[30]/43.045610; /*ch2hco */
    YOW += y[31]/43.045610; /*ch3co */
    YOW += y[32]/44.053580; /*ch3hco */
    YOW += y[33]/46.069520; /*c2h5oh */
    YOW += y[34]/45.061550; /*c2h4oh */
    YOW += y[35]/45.061550; /*ch3choh */
    YOW += y[36]/45.061550; /*ch3ch2o */
    YOW += y[37]/28.013400; /*n2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]/2.015940; 
    c[1] = PWORT * y[1]/1.007970; 
    c[2] = PWORT * y[2]/16.043030; 
    c[3] = PWORT * y[3]/15.035060; 
    c[4] = PWORT * y[4]/14.027090; 
    c[5] = PWORT * y[5]/13.019120; 
    c[6] = PWORT * y[6]/30.026490; 
    c[7] = PWORT * y[7]/29.018520; 
    c[8] = PWORT * y[8]/44.009950; 
    c[9] = PWORT * y[9]/28.010550; 
    c[10] = PWORT * y[10]/31.998800; 
    c[11] = PWORT * y[11]/15.999400; 
    c[12] = PWORT * y[12]/17.007370; 
    c[13] = PWORT * y[13]/33.006770; 
    c[14] = PWORT * y[14]/34.014740; 
    c[15] = PWORT * y[15]/18.015340; 
    c[16] = PWORT * y[16]/25.030270; 
    c[17] = PWORT * y[17]/41.029670; 
    c[18] = PWORT * y[18]/26.038240; 
    c[19] = PWORT * y[19]/27.046210; 
    c[20] = PWORT * y[20]/28.054180; 
    c[21] = PWORT * y[21]/29.062150; 
    c[22] = PWORT * y[22]/30.070120; 
    c[23] = PWORT * y[23]/31.034460; 
    c[24] = PWORT * y[24]/31.034460; 
    c[25] = PWORT * y[25]/42.037640; 
    c[26] = PWORT * y[26]/39.057360; 
    c[27] = PWORT * y[27]/38.049390; 
    c[28] = PWORT * y[28]/14.027090; 
    c[29] = PWORT * y[29]/42.037640; 
    c[30] = PWORT * y[30]/43.045610; 
    c[31] = PWORT * y[31]/43.045610; 
    c[32] = PWORT * y[32]/44.053580; 
    c[33] = PWORT * y[33]/46.069520; 
    c[34] = PWORT * y[34]/45.061550; 
    c[35] = PWORT * y[35]/45.061550; 
    c[36] = PWORT * y[36]/45.061550; 
    c[37] = PWORT * y[37]/28.013400; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 38; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[38]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 38; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 38; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[38]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]/2.015940; 
    c[1] = 1e6 * (*rho) * y[1]/1.007970; 
    c[2] = 1e6 * (*rho) * y[2]/16.043030; 
    c[3] = 1e6 * (*rho) * y[3]/15.035060; 
    c[4] = 1e6 * (*rho) * y[4]/14.027090; 
    c[5] = 1e6 * (*rho) * y[5]/13.019120; 
    c[6] = 1e6 * (*rho) * y[6]/30.026490; 
    c[7] = 1e6 * (*rho) * y[7]/29.018520; 
    c[8] = 1e6 * (*rho) * y[8]/44.009950; 
    c[9] = 1e6 * (*rho) * y[9]/28.010550; 
    c[10] = 1e6 * (*rho) * y[10]/31.998800; 
    c[11] = 1e6 * (*rho) * y[11]/15.999400; 
    c[12] = 1e6 * (*rho) * y[12]/17.007370; 
    c[13] = 1e6 * (*rho) * y[13]/33.006770; 
    c[14] = 1e6 * (*rho) * y[14]/34.014740; 
    c[15] = 1e6 * (*rho) * y[15]/18.015340; 
    c[16] = 1e6 * (*rho) * y[16]/25.030270; 
    c[17] = 1e6 * (*rho) * y[17]/41.029670; 
    c[18] = 1e6 * (*rho) * y[18]/26.038240; 
    c[19] = 1e6 * (*rho) * y[19]/27.046210; 
    c[20] = 1e6 * (*rho) * y[20]/28.054180; 
    c[21] = 1e6 * (*rho) * y[21]/29.062150; 
    c[22] = 1e6 * (*rho) * y[22]/30.070120; 
    c[23] = 1e6 * (*rho) * y[23]/31.034460; 
    c[24] = 1e6 * (*rho) * y[24]/31.034460; 
    c[25] = 1e6 * (*rho) * y[25]/42.037640; 
    c[26] = 1e6 * (*rho) * y[26]/39.057360; 
    c[27] = 1e6 * (*rho) * y[27]/38.049390; 
    c[28] = 1e6 * (*rho) * y[28]/14.027090; 
    c[29] = 1e6 * (*rho) * y[29]/42.037640; 
    c[30] = 1e6 * (*rho) * y[30]/43.045610; 
    c[31] = 1e6 * (*rho) * y[31]/43.045610; 
    c[32] = 1e6 * (*rho) * y[32]/44.053580; 
    c[33] = 1e6 * (*rho) * y[33]/46.069520; 
    c[34] = 1e6 * (*rho) * y[34]/45.061550; 
    c[35] = 1e6 * (*rho) * y[35]/45.061550; 
    c[36] = 1e6 * (*rho) * y[36]/45.061550; 
    c[37] = 1e6 * (*rho) * y[37]/28.013400; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 38; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[38]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*h2 */
    XW += x[1]*1.007970; /*h */
    XW += x[2]*16.043030; /*ch4 */
    XW += x[3]*15.035060; /*ch3 */
    XW += x[4]*14.027090; /*ch2 */
    XW += x[5]*13.019120; /*ch */
    XW += x[6]*30.026490; /*ch2o */
    XW += x[7]*29.018520; /*hco */
    XW += x[8]*44.009950; /*co2 */
    XW += x[9]*28.010550; /*co */
    XW += x[10]*31.998800; /*o2 */
    XW += x[11]*15.999400; /*o */
    XW += x[12]*17.007370; /*oh */
    XW += x[13]*33.006770; /*ho2 */
    XW += x[14]*34.014740; /*h2o2 */
    XW += x[15]*18.015340; /*h2o */
    XW += x[16]*25.030270; /*c2h */
    XW += x[17]*41.029670; /*hcco */
    XW += x[18]*26.038240; /*c2h2 */
    XW += x[19]*27.046210; /*c2h3 */
    XW += x[20]*28.054180; /*c2h4 */
    XW += x[21]*29.062150; /*c2h5 */
    XW += x[22]*30.070120; /*c2h6 */
    XW += x[23]*31.034460; /*ch2oh */
    XW += x[24]*31.034460; /*ch3o */
    XW += x[25]*42.037640; /*hccoh */
    XW += x[26]*39.057360; /*h2ccch */
    XW += x[27]*38.049390; /*c3h2 */
    XW += x[28]*14.027090; /*ch2s */
    XW += x[29]*42.037640; /*ch2co */
    XW += x[30]*43.045610; /*ch2hco */
    XW += x[31]*43.045610; /*ch3co */
    XW += x[32]*44.053580; /*ch3hco */
    XW += x[33]*46.069520; /*c2h5oh */
    XW += x[34]*45.061550; /*c2h4oh */
    XW += x[35]*45.061550; /*ch3choh */
    XW += x[36]*45.061550; /*ch3ch2o */
    XW += x[37]*28.013400; /*n2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 38; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 38; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double * T, double * C, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 38; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 38; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 228; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double * P, double * T, double * x, int * iwrk, double * rwrk, double * q_f, double * q_r)
{
    int id; /*loop counter */
    double c[38]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 38; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 228; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[38]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/2.015940; /*h2 */
    YOW += y[1]/1.007970; /*h */
    YOW += y[2]/16.043030; /*ch4 */
    YOW += y[3]/15.035060; /*ch3 */
    YOW += y[4]/14.027090; /*ch2 */
    YOW += y[5]/13.019120; /*ch */
    YOW += y[6]/30.026490; /*ch2o */
    YOW += y[7]/29.018520; /*hco */
    YOW += y[8]/44.009950; /*co2 */
    YOW += y[9]/28.010550; /*co */
    YOW += y[10]/31.998800; /*o2 */
    YOW += y[11]/15.999400; /*o */
    YOW += y[12]/17.007370; /*oh */
    YOW += y[13]/33.006770; /*ho2 */
    YOW += y[14]/34.014740; /*h2o2 */
    YOW += y[15]/18.015340; /*h2o */
    YOW += y[16]/25.030270; /*c2h */
    YOW += y[17]/41.029670; /*hcco */
    YOW += y[18]/26.038240; /*c2h2 */
    YOW += y[19]/27.046210; /*c2h3 */
    YOW += y[20]/28.054180; /*c2h4 */
    YOW += y[21]/29.062150; /*c2h5 */
    YOW += y[22]/30.070120; /*c2h6 */
    YOW += y[23]/31.034460; /*ch2oh */
    YOW += y[24]/31.034460; /*ch3o */
    YOW += y[25]/42.037640; /*hccoh */
    YOW += y[26]/39.057360; /*h2ccch */
    YOW += y[27]/38.049390; /*c3h2 */
    YOW += y[28]/14.027090; /*ch2s */
    YOW += y[29]/42.037640; /*ch2co */
    YOW += y[30]/43.045610; /*ch2hco */
    YOW += y[31]/43.045610; /*ch3co */
    YOW += y[32]/44.053580; /*ch3hco */
    YOW += y[33]/46.069520; /*c2h5oh */
    YOW += y[34]/45.061550; /*c2h4oh */
    YOW += y[35]/45.061550; /*ch3choh */
    YOW += y[36]/45.061550; /*ch3ch2o */
    YOW += y[37]/28.013400; /*n2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]/2.015940; 
    c[1] = PWORT * y[1]/1.007970; 
    c[2] = PWORT * y[2]/16.043030; 
    c[3] = PWORT * y[3]/15.035060; 
    c[4] = PWORT * y[4]/14.027090; 
    c[5] = PWORT * y[5]/13.019120; 
    c[6] = PWORT * y[6]/30.026490; 
    c[7] = PWORT * y[7]/29.018520; 
    c[8] = PWORT * y[8]/44.009950; 
    c[9] = PWORT * y[9]/28.010550; 
    c[10] = PWORT * y[10]/31.998800; 
    c[11] = PWORT * y[11]/15.999400; 
    c[12] = PWORT * y[12]/17.007370; 
    c[13] = PWORT * y[13]/33.006770; 
    c[14] = PWORT * y[14]/34.014740; 
    c[15] = PWORT * y[15]/18.015340; 
    c[16] = PWORT * y[16]/25.030270; 
    c[17] = PWORT * y[17]/41.029670; 
    c[18] = PWORT * y[18]/26.038240; 
    c[19] = PWORT * y[19]/27.046210; 
    c[20] = PWORT * y[20]/28.054180; 
    c[21] = PWORT * y[21]/29.062150; 
    c[22] = PWORT * y[22]/30.070120; 
    c[23] = PWORT * y[23]/31.034460; 
    c[24] = PWORT * y[24]/31.034460; 
    c[25] = PWORT * y[25]/42.037640; 
    c[26] = PWORT * y[26]/39.057360; 
    c[27] = PWORT * y[27]/38.049390; 
    c[28] = PWORT * y[28]/14.027090; 
    c[29] = PWORT * y[29]/42.037640; 
    c[30] = PWORT * y[30]/43.045610; 
    c[31] = PWORT * y[31]/43.045610; 
    c[32] = PWORT * y[32]/44.053580; 
    c[33] = PWORT * y[33]/46.069520; 
    c[34] = PWORT * y[34]/45.061550; 
    c[35] = PWORT * y[35]/45.061550; 
    c[36] = PWORT * y[36]/45.061550; 
    c[37] = PWORT * y[37]/28.013400; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 228; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[38]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 38; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 228; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[38]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]/2.015940; 
    c[1] = 1e6 * (*rho) * y[1]/1.007970; 
    c[2] = 1e6 * (*rho) * y[2]/16.043030; 
    c[3] = 1e6 * (*rho) * y[3]/15.035060; 
    c[4] = 1e6 * (*rho) * y[4]/14.027090; 
    c[5] = 1e6 * (*rho) * y[5]/13.019120; 
    c[6] = 1e6 * (*rho) * y[6]/30.026490; 
    c[7] = 1e6 * (*rho) * y[7]/29.018520; 
    c[8] = 1e6 * (*rho) * y[8]/44.009950; 
    c[9] = 1e6 * (*rho) * y[9]/28.010550; 
    c[10] = 1e6 * (*rho) * y[10]/31.998800; 
    c[11] = 1e6 * (*rho) * y[11]/15.999400; 
    c[12] = 1e6 * (*rho) * y[12]/17.007370; 
    c[13] = 1e6 * (*rho) * y[13]/33.006770; 
    c[14] = 1e6 * (*rho) * y[14]/34.014740; 
    c[15] = 1e6 * (*rho) * y[15]/18.015340; 
    c[16] = 1e6 * (*rho) * y[16]/25.030270; 
    c[17] = 1e6 * (*rho) * y[17]/41.029670; 
    c[18] = 1e6 * (*rho) * y[18]/26.038240; 
    c[19] = 1e6 * (*rho) * y[19]/27.046210; 
    c[20] = 1e6 * (*rho) * y[20]/28.054180; 
    c[21] = 1e6 * (*rho) * y[21]/29.062150; 
    c[22] = 1e6 * (*rho) * y[22]/30.070120; 
    c[23] = 1e6 * (*rho) * y[23]/31.034460; 
    c[24] = 1e6 * (*rho) * y[24]/31.034460; 
    c[25] = 1e6 * (*rho) * y[25]/42.037640; 
    c[26] = 1e6 * (*rho) * y[26]/39.057360; 
    c[27] = 1e6 * (*rho) * y[27]/38.049390; 
    c[28] = 1e6 * (*rho) * y[28]/14.027090; 
    c[29] = 1e6 * (*rho) * y[29]/42.037640; 
    c[30] = 1e6 * (*rho) * y[30]/43.045610; 
    c[31] = 1e6 * (*rho) * y[31]/43.045610; 
    c[32] = 1e6 * (*rho) * y[32]/44.053580; 
    c[33] = 1e6 * (*rho) * y[33]/46.069520; 
    c[34] = 1e6 * (*rho) * y[34]/45.061550; 
    c[35] = 1e6 * (*rho) * y[35]/45.061550; 
    c[36] = 1e6 * (*rho) * y[36]/45.061550; 
    c[37] = 1e6 * (*rho) * y[37]/28.013400; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 228; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[38]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*h2 */
    XW += x[1]*1.007970; /*h */
    XW += x[2]*16.043030; /*ch4 */
    XW += x[3]*15.035060; /*ch3 */
    XW += x[4]*14.027090; /*ch2 */
    XW += x[5]*13.019120; /*ch */
    XW += x[6]*30.026490; /*ch2o */
    XW += x[7]*29.018520; /*hco */
    XW += x[8]*44.009950; /*co2 */
    XW += x[9]*28.010550; /*co */
    XW += x[10]*31.998800; /*o2 */
    XW += x[11]*15.999400; /*o */
    XW += x[12]*17.007370; /*oh */
    XW += x[13]*33.006770; /*ho2 */
    XW += x[14]*34.014740; /*h2o2 */
    XW += x[15]*18.015340; /*h2o */
    XW += x[16]*25.030270; /*c2h */
    XW += x[17]*41.029670; /*hcco */
    XW += x[18]*26.038240; /*c2h2 */
    XW += x[19]*27.046210; /*c2h3 */
    XW += x[20]*28.054180; /*c2h4 */
    XW += x[21]*29.062150; /*c2h5 */
    XW += x[22]*30.070120; /*c2h6 */
    XW += x[23]*31.034460; /*ch2oh */
    XW += x[24]*31.034460; /*ch3o */
    XW += x[25]*42.037640; /*hccoh */
    XW += x[26]*39.057360; /*h2ccch */
    XW += x[27]*38.049390; /*c3h2 */
    XW += x[28]*14.027090; /*ch2s */
    XW += x[29]*42.037640; /*ch2co */
    XW += x[30]*43.045610; /*ch2hco */
    XW += x[31]*43.045610; /*ch3co */
    XW += x[32]*44.053580; /*ch3hco */
    XW += x[33]*46.069520; /*c2h5oh */
    XW += x[34]*45.061550; /*c2h4oh */
    XW += x[35]*45.061550; /*ch3choh */
    XW += x[36]*45.061550; /*ch3ch2o */
    XW += x[37]*28.013400; /*n2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 38; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 228; ++id) {
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
    for (id = 0; id < 38 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: oh + h2 <=> h + h2o */
    nuki[ 12 * kd + 0 ] += -1 ;
    nuki[ 0 * kd + 0 ] += -1 ;
    nuki[ 1 * kd + 0 ] += +1 ;
    nuki[ 15 * kd + 0 ] += +1 ;

    /*reaction 2: o + oh <=> o2 + h */
    nuki[ 11 * kd + 1 ] += -1 ;
    nuki[ 12 * kd + 1 ] += -1 ;
    nuki[ 10 * kd + 1 ] += +1 ;
    nuki[ 1 * kd + 1 ] += +1 ;

    /*reaction 3: o + h2 <=> oh + h */
    nuki[ 11 * kd + 2 ] += -1 ;
    nuki[ 0 * kd + 2 ] += -1 ;
    nuki[ 12 * kd + 2 ] += +1 ;
    nuki[ 1 * kd + 2 ] += +1 ;

    /*reaction 4: h + o2 (+M) <=> ho2 (+M) */
    nuki[ 1 * kd + 3 ] += -1 ;
    nuki[ 10 * kd + 3 ] += -1 ;
    nuki[ 13 * kd + 3 ] += +1 ;

    /*reaction 5: h + o2 (+n2) <=> ho2 (+n2) */
    nuki[ 1 * kd + 4 ] += -1 ;
    nuki[ 10 * kd + 4 ] += -1 ;
    nuki[ 13 * kd + 4 ] += +1 ;

    /*reaction 6: h + o2 (+h2) <=> ho2 (+h2) */
    nuki[ 1 * kd + 5 ] += -1 ;
    nuki[ 10 * kd + 5 ] += -1 ;
    nuki[ 13 * kd + 5 ] += +1 ;

    /*reaction 7: h + o2 (+h2o) <=> ho2 (+h2o) */
    nuki[ 1 * kd + 6 ] += -1 ;
    nuki[ 10 * kd + 6 ] += -1 ;
    nuki[ 13 * kd + 6 ] += +1 ;

    /*reaction 8: oh + ho2 <=> h2o + o2 */
    nuki[ 12 * kd + 7 ] += -1 ;
    nuki[ 13 * kd + 7 ] += -1 ;
    nuki[ 15 * kd + 7 ] += +1 ;
    nuki[ 10 * kd + 7 ] += +1 ;

    /*reaction 9: oh + ho2 <=> h2o + o2 */
    nuki[ 12 * kd + 8 ] += -1 ;
    nuki[ 13 * kd + 8 ] += -1 ;
    nuki[ 15 * kd + 8 ] += +1 ;
    nuki[ 10 * kd + 8 ] += +1 ;

    /*reaction 10: h + ho2 <=> oh + oh */
    nuki[ 1 * kd + 9 ] += -1 ;
    nuki[ 13 * kd + 9 ] += -1 ;
    nuki[ 12 * kd + 9 ] += +1 ;
    nuki[ 12 * kd + 9 ] += +1 ;

    /*reaction 11: h + ho2 <=> h2 + o2 */
    nuki[ 1 * kd + 10 ] += -1 ;
    nuki[ 13 * kd + 10 ] += -1 ;
    nuki[ 0 * kd + 10 ] += +1 ;
    nuki[ 10 * kd + 10 ] += +1 ;

    /*reaction 12: h + ho2 <=> o + h2o */
    nuki[ 1 * kd + 11 ] += -1 ;
    nuki[ 13 * kd + 11 ] += -1 ;
    nuki[ 11 * kd + 11 ] += +1 ;
    nuki[ 15 * kd + 11 ] += +1 ;

    /*reaction 13: o + ho2 <=> o2 + oh */
    nuki[ 11 * kd + 12 ] += -1 ;
    nuki[ 13 * kd + 12 ] += -1 ;
    nuki[ 10 * kd + 12 ] += +1 ;
    nuki[ 12 * kd + 12 ] += +1 ;

    /*reaction 14: oh + oh <=> o + h2o */
    nuki[ 12 * kd + 13 ] += -1 ;
    nuki[ 12 * kd + 13 ] += -1 ;
    nuki[ 11 * kd + 13 ] += +1 ;
    nuki[ 15 * kd + 13 ] += +1 ;

    /*reaction 15: h + h + M <=> h2 + M */
    nuki[ 1 * kd + 14 ] += -1 ;
    nuki[ 1 * kd + 14 ] += -1 ;
    nuki[ 0 * kd + 14 ] += +1 ;

    /*reaction 16: h + h + h2 <=> h2 + h2 */
    nuki[ 1 * kd + 15 ] += -1 ;
    nuki[ 1 * kd + 15 ] += -1 ;
    nuki[ 0 * kd + 15 ] += -1 ;
    nuki[ 0 * kd + 15 ] += +1 ;
    nuki[ 0 * kd + 15 ] += +1 ;

    /*reaction 17: h + h + h2o <=> h2 + h2o */
    nuki[ 1 * kd + 16 ] += -1 ;
    nuki[ 1 * kd + 16 ] += -1 ;
    nuki[ 15 * kd + 16 ] += -1 ;
    nuki[ 0 * kd + 16 ] += +1 ;
    nuki[ 15 * kd + 16 ] += +1 ;

    /*reaction 18: h + oh + M <=> h2o + M */
    nuki[ 1 * kd + 17 ] += -1 ;
    nuki[ 12 * kd + 17 ] += -1 ;
    nuki[ 15 * kd + 17 ] += +1 ;

    /*reaction 19: h + o + M <=> oh + M */
    nuki[ 1 * kd + 18 ] += -1 ;
    nuki[ 11 * kd + 18 ] += -1 ;
    nuki[ 12 * kd + 18 ] += +1 ;

    /*reaction 20: o + o + M <=> o2 + M */
    nuki[ 11 * kd + 19 ] += -1 ;
    nuki[ 11 * kd + 19 ] += -1 ;
    nuki[ 10 * kd + 19 ] += +1 ;

    /*reaction 21: ho2 + ho2 <=> h2o2 + o2 */
    nuki[ 13 * kd + 20 ] += -1 ;
    nuki[ 13 * kd + 20 ] += -1 ;
    nuki[ 14 * kd + 20 ] += +1 ;
    nuki[ 10 * kd + 20 ] += +1 ;

    /*reaction 22: ho2 + ho2 <=> h2o2 + o2 */
    nuki[ 13 * kd + 21 ] += -1 ;
    nuki[ 13 * kd + 21 ] += -1 ;
    nuki[ 14 * kd + 21 ] += +1 ;
    nuki[ 10 * kd + 21 ] += +1 ;

    /*reaction 23: oh + oh (+M) <=> h2o2 (+M) */
    nuki[ 12 * kd + 22 ] += -1 ;
    nuki[ 12 * kd + 22 ] += -1 ;
    nuki[ 14 * kd + 22 ] += +1 ;

    /*reaction 24: h2o2 + h <=> ho2 + h2 */
    nuki[ 14 * kd + 23 ] += -1 ;
    nuki[ 1 * kd + 23 ] += -1 ;
    nuki[ 13 * kd + 23 ] += +1 ;
    nuki[ 0 * kd + 23 ] += +1 ;

    /*reaction 25: h2o2 + h <=> oh + h2o */
    nuki[ 14 * kd + 24 ] += -1 ;
    nuki[ 1 * kd + 24 ] += -1 ;
    nuki[ 12 * kd + 24 ] += +1 ;
    nuki[ 15 * kd + 24 ] += +1 ;

    /*reaction 26: h2o2 + o <=> oh + ho2 */
    nuki[ 14 * kd + 25 ] += -1 ;
    nuki[ 11 * kd + 25 ] += -1 ;
    nuki[ 12 * kd + 25 ] += +1 ;
    nuki[ 13 * kd + 25 ] += +1 ;

    /*reaction 27: h2o2 + oh <=> h2o + ho2 */
    nuki[ 14 * kd + 26 ] += -1 ;
    nuki[ 12 * kd + 26 ] += -1 ;
    nuki[ 15 * kd + 26 ] += +1 ;
    nuki[ 13 * kd + 26 ] += +1 ;

    /*reaction 28: ch3 + ch3 (+M) <=> c2h6 (+M) */
    nuki[ 3 * kd + 27 ] += -1 ;
    nuki[ 3 * kd + 27 ] += -1 ;
    nuki[ 22 * kd + 27 ] += +1 ;

    /*reaction 29: ch3 + h (+M) <=> ch4 (+M) */
    nuki[ 3 * kd + 28 ] += -1 ;
    nuki[ 1 * kd + 28 ] += -1 ;
    nuki[ 2 * kd + 28 ] += +1 ;

    /*reaction 30: ch4 + h <=> ch3 + h2 */
    nuki[ 2 * kd + 29 ] += -1 ;
    nuki[ 1 * kd + 29 ] += -1 ;
    nuki[ 3 * kd + 29 ] += +1 ;
    nuki[ 0 * kd + 29 ] += +1 ;

    /*reaction 31: ch4 + oh <=> ch3 + h2o */
    nuki[ 2 * kd + 30 ] += -1 ;
    nuki[ 12 * kd + 30 ] += -1 ;
    nuki[ 3 * kd + 30 ] += +1 ;
    nuki[ 15 * kd + 30 ] += +1 ;

    /*reaction 32: ch4 + o <=> ch3 + oh */
    nuki[ 2 * kd + 31 ] += -1 ;
    nuki[ 11 * kd + 31 ] += -1 ;
    nuki[ 3 * kd + 31 ] += +1 ;
    nuki[ 12 * kd + 31 ] += +1 ;

    /*reaction 33: ch4 + ho2 <=> ch3 + h2o2 */
    nuki[ 2 * kd + 32 ] += -1 ;
    nuki[ 13 * kd + 32 ] += -1 ;
    nuki[ 3 * kd + 32 ] += +1 ;
    nuki[ 14 * kd + 32 ] += +1 ;

    /*reaction 34: ch3 + ho2 <=> ch3o + oh */
    nuki[ 3 * kd + 33 ] += -1 ;
    nuki[ 13 * kd + 33 ] += -1 ;
    nuki[ 24 * kd + 33 ] += +1 ;
    nuki[ 12 * kd + 33 ] += +1 ;

    /*reaction 35: ch3 + ho2 <=> ch4 + o2 */
    nuki[ 3 * kd + 34 ] += -1 ;
    nuki[ 13 * kd + 34 ] += -1 ;
    nuki[ 2 * kd + 34 ] += +1 ;
    nuki[ 10 * kd + 34 ] += +1 ;

    /*reaction 36: ch3 + o <=> ch2o + h */
    nuki[ 3 * kd + 35 ] += -1 ;
    nuki[ 11 * kd + 35 ] += -1 ;
    nuki[ 6 * kd + 35 ] += +1 ;
    nuki[ 1 * kd + 35 ] += +1 ;

    /*reaction 37: ch3 + o2 <=> ch3o + o */
    nuki[ 3 * kd + 36 ] += -1 ;
    nuki[ 10 * kd + 36 ] += -1 ;
    nuki[ 24 * kd + 36 ] += +1 ;
    nuki[ 11 * kd + 36 ] += +1 ;

    /*reaction 38: ch3 + o2 <=> ch2o + oh */
    nuki[ 3 * kd + 37 ] += -1 ;
    nuki[ 10 * kd + 37 ] += -1 ;
    nuki[ 6 * kd + 37 ] += +1 ;
    nuki[ 12 * kd + 37 ] += +1 ;

    /*reaction 39: ch3o + h <=> ch3 + oh */
    nuki[ 24 * kd + 38 ] += -1 ;
    nuki[ 1 * kd + 38 ] += -1 ;
    nuki[ 3 * kd + 38 ] += +1 ;
    nuki[ 12 * kd + 38 ] += +1 ;

    /*reaction 40: ch2oh + h <=> ch3 + oh */
    nuki[ 23 * kd + 39 ] += -1 ;
    nuki[ 1 * kd + 39 ] += -1 ;
    nuki[ 3 * kd + 39 ] += +1 ;
    nuki[ 12 * kd + 39 ] += +1 ;

    /*reaction 41: ch3 + oh <=> ch2s + h2o */
    nuki[ 3 * kd + 40 ] += -1 ;
    nuki[ 12 * kd + 40 ] += -1 ;
    nuki[ 28 * kd + 40 ] += +1 ;
    nuki[ 15 * kd + 40 ] += +1 ;

    /*reaction 42: ch3 + oh <=> ch2 + h2o */
    nuki[ 3 * kd + 41 ] += -1 ;
    nuki[ 12 * kd + 41 ] += -1 ;
    nuki[ 4 * kd + 41 ] += +1 ;
    nuki[ 15 * kd + 41 ] += +1 ;

    /*reaction 43: ch3 + h <=> ch2 + h2 */
    nuki[ 3 * kd + 42 ] += -1 ;
    nuki[ 1 * kd + 42 ] += -1 ;
    nuki[ 4 * kd + 42 ] += +1 ;
    nuki[ 0 * kd + 42 ] += +1 ;

    /*reaction 44: ch3 + M <=> ch + h2 + M */
    nuki[ 3 * kd + 43 ] += -1 ;
    nuki[ 5 * kd + 43 ] += +1 ;
    nuki[ 0 * kd + 43 ] += +1 ;

    /*reaction 45: ch3 + M <=> ch2 + h + M */
    nuki[ 3 * kd + 44 ] += -1 ;
    nuki[ 4 * kd + 44 ] += +1 ;
    nuki[ 1 * kd + 44 ] += +1 ;

    /*reaction 46: ch2o + h (+M) <=> ch3o (+M) */
    nuki[ 6 * kd + 45 ] += -1 ;
    nuki[ 1 * kd + 45 ] += -1 ;
    nuki[ 24 * kd + 45 ] += +1 ;

    /*reaction 47: ch2o + h (+M) <=> ch2oh (+M) */
    nuki[ 6 * kd + 46 ] += -1 ;
    nuki[ 1 * kd + 46 ] += -1 ;
    nuki[ 23 * kd + 46 ] += +1 ;

    /*reaction 48: ch3o + ch3 <=> ch2o + ch4 */
    nuki[ 24 * kd + 47 ] += -1 ;
    nuki[ 3 * kd + 47 ] += -1 ;
    nuki[ 6 * kd + 47 ] += +1 ;
    nuki[ 2 * kd + 47 ] += +1 ;

    /*reaction 49: ch3o + h <=> ch2o + h2 */
    nuki[ 24 * kd + 48 ] += -1 ;
    nuki[ 1 * kd + 48 ] += -1 ;
    nuki[ 6 * kd + 48 ] += +1 ;
    nuki[ 0 * kd + 48 ] += +1 ;

    /*reaction 50: ch2oh + h <=> ch2o + h2 */
    nuki[ 23 * kd + 49 ] += -1 ;
    nuki[ 1 * kd + 49 ] += -1 ;
    nuki[ 6 * kd + 49 ] += +1 ;
    nuki[ 0 * kd + 49 ] += +1 ;

    /*reaction 51: ch3o + oh <=> ch2o + h2o */
    nuki[ 24 * kd + 50 ] += -1 ;
    nuki[ 12 * kd + 50 ] += -1 ;
    nuki[ 6 * kd + 50 ] += +1 ;
    nuki[ 15 * kd + 50 ] += +1 ;

    /*reaction 52: ch2oh + oh <=> ch2o + h2o */
    nuki[ 23 * kd + 51 ] += -1 ;
    nuki[ 12 * kd + 51 ] += -1 ;
    nuki[ 6 * kd + 51 ] += +1 ;
    nuki[ 15 * kd + 51 ] += +1 ;

    /*reaction 53: ch3o + o <=> ch2o + oh */
    nuki[ 24 * kd + 52 ] += -1 ;
    nuki[ 11 * kd + 52 ] += -1 ;
    nuki[ 6 * kd + 52 ] += +1 ;
    nuki[ 12 * kd + 52 ] += +1 ;

    /*reaction 54: ch2oh + o <=> ch2o + oh */
    nuki[ 23 * kd + 53 ] += -1 ;
    nuki[ 11 * kd + 53 ] += -1 ;
    nuki[ 6 * kd + 53 ] += +1 ;
    nuki[ 12 * kd + 53 ] += +1 ;

    /*reaction 55: ch3o + o2 <=> ch2o + ho2 */
    nuki[ 24 * kd + 54 ] += -1 ;
    nuki[ 10 * kd + 54 ] += -1 ;
    nuki[ 6 * kd + 54 ] += +1 ;
    nuki[ 13 * kd + 54 ] += +1 ;

    /*reaction 56: ch3o + co <=> ch3 + co2 */
    nuki[ 24 * kd + 55 ] += -1 ;
    nuki[ 9 * kd + 55 ] += -1 ;
    nuki[ 3 * kd + 55 ] += +1 ;
    nuki[ 8 * kd + 55 ] += +1 ;

    /*reaction 57: ch2oh + o2 <=> ch2o + ho2 */
    nuki[ 23 * kd + 56 ] += -1 ;
    nuki[ 10 * kd + 56 ] += -1 ;
    nuki[ 6 * kd + 56 ] += +1 ;
    nuki[ 13 * kd + 56 ] += +1 ;

    /*reaction 58: ch2oh + o2 <=> ch2o + ho2 */
    nuki[ 23 * kd + 57 ] += -1 ;
    nuki[ 10 * kd + 57 ] += -1 ;
    nuki[ 6 * kd + 57 ] += +1 ;
    nuki[ 13 * kd + 57 ] += +1 ;

    /*reaction 59: ch2 + h <=> ch + h2 */
    nuki[ 4 * kd + 58 ] += -1 ;
    nuki[ 1 * kd + 58 ] += -1 ;
    nuki[ 5 * kd + 58 ] += +1 ;
    nuki[ 0 * kd + 58 ] += +1 ;

    /*reaction 60: ch2 + oh <=> ch + h2o */
    nuki[ 4 * kd + 59 ] += -1 ;
    nuki[ 12 * kd + 59 ] += -1 ;
    nuki[ 5 * kd + 59 ] += +1 ;
    nuki[ 15 * kd + 59 ] += +1 ;

    /*reaction 61: ch2 + oh <=> ch2o + h */
    nuki[ 4 * kd + 60 ] += -1 ;
    nuki[ 12 * kd + 60 ] += -1 ;
    nuki[ 6 * kd + 60 ] += +1 ;
    nuki[ 1 * kd + 60 ] += +1 ;

    /*reaction 62: ch2 + co2 <=> ch2o + co */
    nuki[ 4 * kd + 61 ] += -1 ;
    nuki[ 8 * kd + 61 ] += -1 ;
    nuki[ 6 * kd + 61 ] += +1 ;
    nuki[ 9 * kd + 61 ] += +1 ;

    /*reaction 63: ch2 + o <=> co + h + h */
    nuki[ 4 * kd + 62 ] += -1 ;
    nuki[ 11 * kd + 62 ] += -1 ;
    nuki[ 9 * kd + 62 ] += +1 ;
    nuki[ 1 * kd + 62 ] += +1 ;
    nuki[ 1 * kd + 62 ] += +1 ;

    /*reaction 64: ch2 + o <=> co + h2 */
    nuki[ 4 * kd + 63 ] += -1 ;
    nuki[ 11 * kd + 63 ] += -1 ;
    nuki[ 9 * kd + 63 ] += +1 ;
    nuki[ 0 * kd + 63 ] += +1 ;

    /*reaction 65: ch2 + o2 <=> ch2o + o */
    nuki[ 4 * kd + 64 ] += -1 ;
    nuki[ 10 * kd + 64 ] += -1 ;
    nuki[ 6 * kd + 64 ] += +1 ;
    nuki[ 11 * kd + 64 ] += +1 ;

    /*reaction 66: ch2 + o2 <=> co2 + h + h */
    nuki[ 4 * kd + 65 ] += -1 ;
    nuki[ 10 * kd + 65 ] += -1 ;
    nuki[ 8 * kd + 65 ] += +1 ;
    nuki[ 1 * kd + 65 ] += +1 ;
    nuki[ 1 * kd + 65 ] += +1 ;

    /*reaction 67: ch2 + o2 <=> co2 + h2 */
    nuki[ 4 * kd + 66 ] += -1 ;
    nuki[ 10 * kd + 66 ] += -1 ;
    nuki[ 8 * kd + 66 ] += +1 ;
    nuki[ 0 * kd + 66 ] += +1 ;

    /*reaction 68: ch2 + o2 <=> co + h2o */
    nuki[ 4 * kd + 67 ] += -1 ;
    nuki[ 10 * kd + 67 ] += -1 ;
    nuki[ 9 * kd + 67 ] += +1 ;
    nuki[ 15 * kd + 67 ] += +1 ;

    /*reaction 69: ch2 + o2 <=> hco + oh */
    nuki[ 4 * kd + 68 ] += -1 ;
    nuki[ 10 * kd + 68 ] += -1 ;
    nuki[ 7 * kd + 68 ] += +1 ;
    nuki[ 12 * kd + 68 ] += +1 ;

    /*reaction 70: ch2 + ch3 <=> c2h4 + h */
    nuki[ 4 * kd + 69 ] += -1 ;
    nuki[ 3 * kd + 69 ] += -1 ;
    nuki[ 20 * kd + 69 ] += +1 ;
    nuki[ 1 * kd + 69 ] += +1 ;

    /*reaction 71: ch2 + ch2 <=> c2h2 + h + h */
    nuki[ 4 * kd + 70 ] += -1 ;
    nuki[ 4 * kd + 70 ] += -1 ;
    nuki[ 18 * kd + 70 ] += +1 ;
    nuki[ 1 * kd + 70 ] += +1 ;
    nuki[ 1 * kd + 70 ] += +1 ;

    /*reaction 72: ch2 + hcco <=> c2h3 + co */
    nuki[ 4 * kd + 71 ] += -1 ;
    nuki[ 17 * kd + 71 ] += -1 ;
    nuki[ 19 * kd + 71 ] += +1 ;
    nuki[ 9 * kd + 71 ] += +1 ;

    /*reaction 73: ch2 + c2h2 <=> h2ccch + h */
    nuki[ 4 * kd + 72 ] += -1 ;
    nuki[ 18 * kd + 72 ] += -1 ;
    nuki[ 26 * kd + 72 ] += +1 ;
    nuki[ 1 * kd + 72 ] += +1 ;

    /*reaction 74: ch2s + M <=> ch2 + M */
    nuki[ 28 * kd + 73 ] += -1 ;
    nuki[ 4 * kd + 73 ] += +1 ;

    /*reaction 75: ch2s + ch4 <=> ch3 + ch3 */
    nuki[ 28 * kd + 74 ] += -1 ;
    nuki[ 2 * kd + 74 ] += -1 ;
    nuki[ 3 * kd + 74 ] += +1 ;
    nuki[ 3 * kd + 74 ] += +1 ;

    /*reaction 76: ch2s + c2h6 <=> ch3 + c2h5 */
    nuki[ 28 * kd + 75 ] += -1 ;
    nuki[ 22 * kd + 75 ] += -1 ;
    nuki[ 3 * kd + 75 ] += +1 ;
    nuki[ 21 * kd + 75 ] += +1 ;

    /*reaction 77: ch2s + o2 <=> co + oh + h */
    nuki[ 28 * kd + 76 ] += -1 ;
    nuki[ 10 * kd + 76 ] += -1 ;
    nuki[ 9 * kd + 76 ] += +1 ;
    nuki[ 12 * kd + 76 ] += +1 ;
    nuki[ 1 * kd + 76 ] += +1 ;

    /*reaction 78: ch2s + h2 <=> ch3 + h */
    nuki[ 28 * kd + 77 ] += -1 ;
    nuki[ 0 * kd + 77 ] += -1 ;
    nuki[ 3 * kd + 77 ] += +1 ;
    nuki[ 1 * kd + 77 ] += +1 ;

    /*reaction 79: ch2s + c2h2 <=> h2ccch + h */
    nuki[ 28 * kd + 78 ] += -1 ;
    nuki[ 18 * kd + 78 ] += -1 ;
    nuki[ 26 * kd + 78 ] += +1 ;
    nuki[ 1 * kd + 78 ] += +1 ;

    /*reaction 80: ch2s + o <=> co + h + h */
    nuki[ 28 * kd + 79 ] += -1 ;
    nuki[ 11 * kd + 79 ] += -1 ;
    nuki[ 9 * kd + 79 ] += +1 ;
    nuki[ 1 * kd + 79 ] += +1 ;
    nuki[ 1 * kd + 79 ] += +1 ;

    /*reaction 81: ch2s + oh <=> ch2o + h */
    nuki[ 28 * kd + 80 ] += -1 ;
    nuki[ 12 * kd + 80 ] += -1 ;
    nuki[ 6 * kd + 80 ] += +1 ;
    nuki[ 1 * kd + 80 ] += +1 ;

    /*reaction 82: ch2s + h <=> ch + h2 */
    nuki[ 28 * kd + 81 ] += -1 ;
    nuki[ 1 * kd + 81 ] += -1 ;
    nuki[ 5 * kd + 81 ] += +1 ;
    nuki[ 0 * kd + 81 ] += +1 ;

    /*reaction 83: ch2s + co2 <=> ch2o + co */
    nuki[ 28 * kd + 82 ] += -1 ;
    nuki[ 8 * kd + 82 ] += -1 ;
    nuki[ 6 * kd + 82 ] += +1 ;
    nuki[ 9 * kd + 82 ] += +1 ;

    /*reaction 84: ch2s + ch3 <=> c2h4 + h */
    nuki[ 28 * kd + 83 ] += -1 ;
    nuki[ 3 * kd + 83 ] += -1 ;
    nuki[ 20 * kd + 83 ] += +1 ;
    nuki[ 1 * kd + 83 ] += +1 ;

    /*reaction 85: ch2s + ch2co <=> c2h4 + co */
    nuki[ 28 * kd + 84 ] += -1 ;
    nuki[ 29 * kd + 84 ] += -1 ;
    nuki[ 20 * kd + 84 ] += +1 ;
    nuki[ 9 * kd + 84 ] += +1 ;

    /*reaction 86: ch + o2 <=> hco + o */
    nuki[ 5 * kd + 85 ] += -1 ;
    nuki[ 10 * kd + 85 ] += -1 ;
    nuki[ 7 * kd + 85 ] += +1 ;
    nuki[ 11 * kd + 85 ] += +1 ;

    /*reaction 87: ch + o <=> co + h */
    nuki[ 5 * kd + 86 ] += -1 ;
    nuki[ 11 * kd + 86 ] += -1 ;
    nuki[ 9 * kd + 86 ] += +1 ;
    nuki[ 1 * kd + 86 ] += +1 ;

    /*reaction 88: ch + oh <=> hco + h */
    nuki[ 5 * kd + 87 ] += -1 ;
    nuki[ 12 * kd + 87 ] += -1 ;
    nuki[ 7 * kd + 87 ] += +1 ;
    nuki[ 1 * kd + 87 ] += +1 ;

    /*reaction 89: ch + co2 <=> hco + co */
    nuki[ 5 * kd + 88 ] += -1 ;
    nuki[ 8 * kd + 88 ] += -1 ;
    nuki[ 7 * kd + 88 ] += +1 ;
    nuki[ 9 * kd + 88 ] += +1 ;

    /*reaction 90: ch + h2o <=> ch2o + h */
    nuki[ 5 * kd + 89 ] += -1 ;
    nuki[ 15 * kd + 89 ] += -1 ;
    nuki[ 6 * kd + 89 ] += +1 ;
    nuki[ 1 * kd + 89 ] += +1 ;

    /*reaction 91: ch + ch2o <=> ch2co + h */
    nuki[ 5 * kd + 90 ] += -1 ;
    nuki[ 6 * kd + 90 ] += -1 ;
    nuki[ 29 * kd + 90 ] += +1 ;
    nuki[ 1 * kd + 90 ] += +1 ;

    /*reaction 92: ch + c2h2 <=> c3h2 + h */
    nuki[ 5 * kd + 91 ] += -1 ;
    nuki[ 18 * kd + 91 ] += -1 ;
    nuki[ 27 * kd + 91 ] += +1 ;
    nuki[ 1 * kd + 91 ] += +1 ;

    /*reaction 93: ch + ch2 <=> c2h2 + h */
    nuki[ 5 * kd + 92 ] += -1 ;
    nuki[ 4 * kd + 92 ] += -1 ;
    nuki[ 18 * kd + 92 ] += +1 ;
    nuki[ 1 * kd + 92 ] += +1 ;

    /*reaction 94: ch + ch3 <=> c2h3 + h */
    nuki[ 5 * kd + 93 ] += -1 ;
    nuki[ 3 * kd + 93 ] += -1 ;
    nuki[ 19 * kd + 93 ] += +1 ;
    nuki[ 1 * kd + 93 ] += +1 ;

    /*reaction 95: ch + ch4 <=> c2h4 + h */
    nuki[ 5 * kd + 94 ] += -1 ;
    nuki[ 2 * kd + 94 ] += -1 ;
    nuki[ 20 * kd + 94 ] += +1 ;
    nuki[ 1 * kd + 94 ] += +1 ;

    /*reaction 96: ch2o + oh <=> hco + h2o */
    nuki[ 6 * kd + 95 ] += -1 ;
    nuki[ 12 * kd + 95 ] += -1 ;
    nuki[ 7 * kd + 95 ] += +1 ;
    nuki[ 15 * kd + 95 ] += +1 ;

    /*reaction 97: ch2o + h <=> hco + h2 */
    nuki[ 6 * kd + 96 ] += -1 ;
    nuki[ 1 * kd + 96 ] += -1 ;
    nuki[ 7 * kd + 96 ] += +1 ;
    nuki[ 0 * kd + 96 ] += +1 ;

    /*reaction 98: ch2o + M <=> hco + h + M */
    nuki[ 6 * kd + 97 ] += -1 ;
    nuki[ 7 * kd + 97 ] += +1 ;
    nuki[ 1 * kd + 97 ] += +1 ;

    /*reaction 99: ch2o + o <=> hco + oh */
    nuki[ 6 * kd + 98 ] += -1 ;
    nuki[ 11 * kd + 98 ] += -1 ;
    nuki[ 7 * kd + 98 ] += +1 ;
    nuki[ 12 * kd + 98 ] += +1 ;

    /*reaction 100: hco + o2 <=> co + ho2 */
    nuki[ 7 * kd + 99 ] += -1 ;
    nuki[ 10 * kd + 99 ] += -1 ;
    nuki[ 9 * kd + 99 ] += +1 ;
    nuki[ 13 * kd + 99 ] += +1 ;

    /*reaction 101: hco + M <=> h + co + M */
    nuki[ 7 * kd + 100 ] += -1 ;
    nuki[ 1 * kd + 100 ] += +1 ;
    nuki[ 9 * kd + 100 ] += +1 ;

    /*reaction 102: hco + oh <=> h2o + co */
    nuki[ 7 * kd + 101 ] += -1 ;
    nuki[ 12 * kd + 101 ] += -1 ;
    nuki[ 15 * kd + 101 ] += +1 ;
    nuki[ 9 * kd + 101 ] += +1 ;

    /*reaction 103: hco + h <=> co + h2 */
    nuki[ 7 * kd + 102 ] += -1 ;
    nuki[ 1 * kd + 102 ] += -1 ;
    nuki[ 9 * kd + 102 ] += +1 ;
    nuki[ 0 * kd + 102 ] += +1 ;

    /*reaction 104: hco + o <=> co + oh */
    nuki[ 7 * kd + 103 ] += -1 ;
    nuki[ 11 * kd + 103 ] += -1 ;
    nuki[ 9 * kd + 103 ] += +1 ;
    nuki[ 12 * kd + 103 ] += +1 ;

    /*reaction 105: hco + o <=> co2 + h */
    nuki[ 7 * kd + 104 ] += -1 ;
    nuki[ 11 * kd + 104 ] += -1 ;
    nuki[ 8 * kd + 104 ] += +1 ;
    nuki[ 1 * kd + 104 ] += +1 ;

    /*reaction 106: co + oh <=> co2 + h */
    nuki[ 9 * kd + 105 ] += -1 ;
    nuki[ 12 * kd + 105 ] += -1 ;
    nuki[ 8 * kd + 105 ] += +1 ;
    nuki[ 1 * kd + 105 ] += +1 ;

    /*reaction 107: co + o + M <=> co2 + M */
    nuki[ 9 * kd + 106 ] += -1 ;
    nuki[ 11 * kd + 106 ] += -1 ;
    nuki[ 8 * kd + 106 ] += +1 ;

    /*reaction 108: co + o2 <=> co2 + o */
    nuki[ 9 * kd + 107 ] += -1 ;
    nuki[ 10 * kd + 107 ] += -1 ;
    nuki[ 8 * kd + 107 ] += +1 ;
    nuki[ 11 * kd + 107 ] += +1 ;

    /*reaction 109: co + ho2 <=> co2 + oh */
    nuki[ 9 * kd + 108 ] += -1 ;
    nuki[ 13 * kd + 108 ] += -1 ;
    nuki[ 8 * kd + 108 ] += +1 ;
    nuki[ 12 * kd + 108 ] += +1 ;

    /*reaction 110: c2h5oh (+M) <=> ch3 + ch2oh (+M) */
    nuki[ 33 * kd + 109 ] += -1 ;
    nuki[ 3 * kd + 109 ] += +1 ;
    nuki[ 23 * kd + 109 ] += +1 ;

    /*reaction 111: c2h5oh (+M) <=> c2h5 + oh (+M) */
    nuki[ 33 * kd + 110 ] += -1 ;
    nuki[ 21 * kd + 110 ] += +1 ;
    nuki[ 12 * kd + 110 ] += +1 ;

    /*reaction 112: c2h5oh (+M) <=> c2h4 + h2o (+M) */
    nuki[ 33 * kd + 111 ] += -1 ;
    nuki[ 20 * kd + 111 ] += +1 ;
    nuki[ 15 * kd + 111 ] += +1 ;

    /*reaction 113: c2h5oh (+M) <=> ch3hco + h2 (+M) */
    nuki[ 33 * kd + 112 ] += -1 ;
    nuki[ 32 * kd + 112 ] += +1 ;
    nuki[ 0 * kd + 112 ] += +1 ;

    /*reaction 114: c2h5oh + oh <=> c2h4oh + h2o */
    nuki[ 33 * kd + 113 ] += -1 ;
    nuki[ 12 * kd + 113 ] += -1 ;
    nuki[ 34 * kd + 113 ] += +1 ;
    nuki[ 15 * kd + 113 ] += +1 ;

    /*reaction 115: c2h5oh + oh <=> ch3choh + h2o */
    nuki[ 33 * kd + 114 ] += -1 ;
    nuki[ 12 * kd + 114 ] += -1 ;
    nuki[ 35 * kd + 114 ] += +1 ;
    nuki[ 15 * kd + 114 ] += +1 ;

    /*reaction 116: c2h5oh + oh <=> ch3ch2o + h2o */
    nuki[ 33 * kd + 115 ] += -1 ;
    nuki[ 12 * kd + 115 ] += -1 ;
    nuki[ 36 * kd + 115 ] += +1 ;
    nuki[ 15 * kd + 115 ] += +1 ;

    /*reaction 117: c2h5oh + h <=> c2h4oh + h2 */
    nuki[ 33 * kd + 116 ] += -1 ;
    nuki[ 1 * kd + 116 ] += -1 ;
    nuki[ 34 * kd + 116 ] += +1 ;
    nuki[ 0 * kd + 116 ] += +1 ;

    /*reaction 118: c2h5oh + h <=> ch3choh + h2 */
    nuki[ 33 * kd + 117 ] += -1 ;
    nuki[ 1 * kd + 117 ] += -1 ;
    nuki[ 35 * kd + 117 ] += +1 ;
    nuki[ 0 * kd + 117 ] += +1 ;

    /*reaction 119: c2h5oh + h <=> ch3ch2o + h2 */
    nuki[ 33 * kd + 118 ] += -1 ;
    nuki[ 1 * kd + 118 ] += -1 ;
    nuki[ 36 * kd + 118 ] += +1 ;
    nuki[ 0 * kd + 118 ] += +1 ;

    /*reaction 120: c2h5oh + o <=> c2h4oh + oh */
    nuki[ 33 * kd + 119 ] += -1 ;
    nuki[ 11 * kd + 119 ] += -1 ;
    nuki[ 34 * kd + 119 ] += +1 ;
    nuki[ 12 * kd + 119 ] += +1 ;

    /*reaction 121: c2h5oh + o <=> ch3choh + oh */
    nuki[ 33 * kd + 120 ] += -1 ;
    nuki[ 11 * kd + 120 ] += -1 ;
    nuki[ 35 * kd + 120 ] += +1 ;
    nuki[ 12 * kd + 120 ] += +1 ;

    /*reaction 122: c2h5oh + o <=> ch3ch2o + oh */
    nuki[ 33 * kd + 121 ] += -1 ;
    nuki[ 11 * kd + 121 ] += -1 ;
    nuki[ 36 * kd + 121 ] += +1 ;
    nuki[ 12 * kd + 121 ] += +1 ;

    /*reaction 123: c2h5oh + ch3 <=> c2h4oh + ch4 */
    nuki[ 33 * kd + 122 ] += -1 ;
    nuki[ 3 * kd + 122 ] += -1 ;
    nuki[ 34 * kd + 122 ] += +1 ;
    nuki[ 2 * kd + 122 ] += +1 ;

    /*reaction 124: c2h5oh + ch3 <=> ch3choh + ch4 */
    nuki[ 33 * kd + 123 ] += -1 ;
    nuki[ 3 * kd + 123 ] += -1 ;
    nuki[ 35 * kd + 123 ] += +1 ;
    nuki[ 2 * kd + 123 ] += +1 ;

    /*reaction 125: c2h5oh + ch3 <=> ch3ch2o + ch4 */
    nuki[ 33 * kd + 124 ] += -1 ;
    nuki[ 3 * kd + 124 ] += -1 ;
    nuki[ 36 * kd + 124 ] += +1 ;
    nuki[ 2 * kd + 124 ] += +1 ;

    /*reaction 126: c2h5oh + ho2 <=> ch3choh + h2o2 */
    nuki[ 33 * kd + 125 ] += -1 ;
    nuki[ 13 * kd + 125 ] += -1 ;
    nuki[ 35 * kd + 125 ] += +1 ;
    nuki[ 14 * kd + 125 ] += +1 ;

    /*reaction 127: c2h5oh + ho2 <=> c2h4oh + h2o2 */
    nuki[ 33 * kd + 126 ] += -1 ;
    nuki[ 13 * kd + 126 ] += -1 ;
    nuki[ 34 * kd + 126 ] += +1 ;
    nuki[ 14 * kd + 126 ] += +1 ;

    /*reaction 128: c2h5oh + ho2 <=> ch3ch2o + h2o2 */
    nuki[ 33 * kd + 127 ] += -1 ;
    nuki[ 13 * kd + 127 ] += -1 ;
    nuki[ 36 * kd + 127 ] += +1 ;
    nuki[ 14 * kd + 127 ] += +1 ;

    /*reaction 129: ch3ch2o + M <=> ch3hco + h + M */
    nuki[ 36 * kd + 128 ] += -1 ;
    nuki[ 32 * kd + 128 ] += +1 ;
    nuki[ 1 * kd + 128 ] += +1 ;

    /*reaction 130: ch3ch2o + M <=> ch3 + ch2o + M */
    nuki[ 36 * kd + 129 ] += -1 ;
    nuki[ 3 * kd + 129 ] += +1 ;
    nuki[ 6 * kd + 129 ] += +1 ;

    /*reaction 131: ch3ch2o + o2 <=> ch3hco + ho2 */
    nuki[ 36 * kd + 130 ] += -1 ;
    nuki[ 10 * kd + 130 ] += -1 ;
    nuki[ 32 * kd + 130 ] += +1 ;
    nuki[ 13 * kd + 130 ] += +1 ;

    /*reaction 132: ch3ch2o + co <=> c2h5 + co2 */
    nuki[ 36 * kd + 131 ] += -1 ;
    nuki[ 9 * kd + 131 ] += -1 ;
    nuki[ 21 * kd + 131 ] += +1 ;
    nuki[ 8 * kd + 131 ] += +1 ;

    /*reaction 133: ch3ch2o + h <=> ch3 + ch2oh */
    nuki[ 36 * kd + 132 ] += -1 ;
    nuki[ 1 * kd + 132 ] += -1 ;
    nuki[ 3 * kd + 132 ] += +1 ;
    nuki[ 23 * kd + 132 ] += +1 ;

    /*reaction 134: ch3ch2o + h <=> c2h4 + h2o */
    nuki[ 36 * kd + 133 ] += -1 ;
    nuki[ 1 * kd + 133 ] += -1 ;
    nuki[ 20 * kd + 133 ] += +1 ;
    nuki[ 15 * kd + 133 ] += +1 ;

    /*reaction 135: ch3ch2o + oh <=> ch3hco + h2o */
    nuki[ 36 * kd + 134 ] += -1 ;
    nuki[ 12 * kd + 134 ] += -1 ;
    nuki[ 32 * kd + 134 ] += +1 ;
    nuki[ 15 * kd + 134 ] += +1 ;

    /*reaction 136: ch3choh + o2 <=> ch3hco + ho2 */
    nuki[ 35 * kd + 135 ] += -1 ;
    nuki[ 10 * kd + 135 ] += -1 ;
    nuki[ 32 * kd + 135 ] += +1 ;
    nuki[ 13 * kd + 135 ] += +1 ;

    /*reaction 137: ch3choh + o2 <=> ch3hco + ho2 */
    nuki[ 35 * kd + 136 ] += -1 ;
    nuki[ 10 * kd + 136 ] += -1 ;
    nuki[ 32 * kd + 136 ] += +1 ;
    nuki[ 13 * kd + 136 ] += +1 ;

    /*reaction 138: ch3choh + o <=> ch3hco + oh */
    nuki[ 35 * kd + 137 ] += -1 ;
    nuki[ 11 * kd + 137 ] += -1 ;
    nuki[ 32 * kd + 137 ] += +1 ;
    nuki[ 12 * kd + 137 ] += +1 ;

    /*reaction 139: ch3choh + h <=> c2h4 + h2o */
    nuki[ 35 * kd + 138 ] += -1 ;
    nuki[ 1 * kd + 138 ] += -1 ;
    nuki[ 20 * kd + 138 ] += +1 ;
    nuki[ 15 * kd + 138 ] += +1 ;

    /*reaction 140: ch3choh + h <=> ch3 + ch2oh */
    nuki[ 35 * kd + 139 ] += -1 ;
    nuki[ 1 * kd + 139 ] += -1 ;
    nuki[ 3 * kd + 139 ] += +1 ;
    nuki[ 23 * kd + 139 ] += +1 ;

    /*reaction 141: ch3choh + ho2 <=> ch3hco + oh + oh */
    nuki[ 35 * kd + 140 ] += -1 ;
    nuki[ 13 * kd + 140 ] += -1 ;
    nuki[ 32 * kd + 140 ] += +1 ;
    nuki[ 12 * kd + 140 ] += +1 ;
    nuki[ 12 * kd + 140 ] += +1 ;

    /*reaction 142: ch3choh + oh <=> ch3hco + h2o */
    nuki[ 35 * kd + 141 ] += -1 ;
    nuki[ 12 * kd + 141 ] += -1 ;
    nuki[ 32 * kd + 141 ] += +1 ;
    nuki[ 15 * kd + 141 ] += +1 ;

    /*reaction 143: ch3choh + M <=> ch3hco + h + M */
    nuki[ 35 * kd + 142 ] += -1 ;
    nuki[ 32 * kd + 142 ] += +1 ;
    nuki[ 1 * kd + 142 ] += +1 ;

    /*reaction 144: ch3hco + oh <=> ch3co + h2o */
    nuki[ 32 * kd + 143 ] += -1 ;
    nuki[ 12 * kd + 143 ] += -1 ;
    nuki[ 31 * kd + 143 ] += +1 ;
    nuki[ 15 * kd + 143 ] += +1 ;

    /*reaction 145: ch3hco + oh <=> ch2hco + h2o */
    nuki[ 32 * kd + 144 ] += -1 ;
    nuki[ 12 * kd + 144 ] += -1 ;
    nuki[ 30 * kd + 144 ] += +1 ;
    nuki[ 15 * kd + 144 ] += +1 ;

    /*reaction 146: ch3hco + o <=> ch3co + oh */
    nuki[ 32 * kd + 145 ] += -1 ;
    nuki[ 11 * kd + 145 ] += -1 ;
    nuki[ 31 * kd + 145 ] += +1 ;
    nuki[ 12 * kd + 145 ] += +1 ;

    /*reaction 147: ch3hco + o <=> ch2hco + oh */
    nuki[ 32 * kd + 146 ] += -1 ;
    nuki[ 11 * kd + 146 ] += -1 ;
    nuki[ 30 * kd + 146 ] += +1 ;
    nuki[ 12 * kd + 146 ] += +1 ;

    /*reaction 148: ch3hco + h <=> ch3co + h2 */
    nuki[ 32 * kd + 147 ] += -1 ;
    nuki[ 1 * kd + 147 ] += -1 ;
    nuki[ 31 * kd + 147 ] += +1 ;
    nuki[ 0 * kd + 147 ] += +1 ;

    /*reaction 149: ch3hco + h <=> ch2hco + h2 */
    nuki[ 32 * kd + 148 ] += -1 ;
    nuki[ 1 * kd + 148 ] += -1 ;
    nuki[ 30 * kd + 148 ] += +1 ;
    nuki[ 0 * kd + 148 ] += +1 ;

    /*reaction 150: ch3hco + ch3 <=> ch3co + ch4 */
    nuki[ 32 * kd + 149 ] += -1 ;
    nuki[ 3 * kd + 149 ] += -1 ;
    nuki[ 31 * kd + 149 ] += +1 ;
    nuki[ 2 * kd + 149 ] += +1 ;

    /*reaction 151: ch3hco + ch3 <=> ch2hco + ch4 */
    nuki[ 32 * kd + 150 ] += -1 ;
    nuki[ 3 * kd + 150 ] += -1 ;
    nuki[ 30 * kd + 150 ] += +1 ;
    nuki[ 2 * kd + 150 ] += +1 ;

    /*reaction 152: ch3hco + ho2 <=> ch3co + h2o2 */
    nuki[ 32 * kd + 151 ] += -1 ;
    nuki[ 13 * kd + 151 ] += -1 ;
    nuki[ 31 * kd + 151 ] += +1 ;
    nuki[ 14 * kd + 151 ] += +1 ;

    /*reaction 153: ch3hco + ho2 <=> ch2hco + h2o2 */
    nuki[ 32 * kd + 152 ] += -1 ;
    nuki[ 13 * kd + 152 ] += -1 ;
    nuki[ 30 * kd + 152 ] += +1 ;
    nuki[ 14 * kd + 152 ] += +1 ;

    /*reaction 154: ch3hco + o2 <=> ch3co + ho2 */
    nuki[ 32 * kd + 153 ] += -1 ;
    nuki[ 10 * kd + 153 ] += -1 ;
    nuki[ 31 * kd + 153 ] += +1 ;
    nuki[ 13 * kd + 153 ] += +1 ;

    /*reaction 155: c2h6 + ch3 <=> c2h5 + ch4 */
    nuki[ 22 * kd + 154 ] += -1 ;
    nuki[ 3 * kd + 154 ] += -1 ;
    nuki[ 21 * kd + 154 ] += +1 ;
    nuki[ 2 * kd + 154 ] += +1 ;

    /*reaction 156: c2h6 + h <=> c2h5 + h2 */
    nuki[ 22 * kd + 155 ] += -1 ;
    nuki[ 1 * kd + 155 ] += -1 ;
    nuki[ 21 * kd + 155 ] += +1 ;
    nuki[ 0 * kd + 155 ] += +1 ;

    /*reaction 157: c2h6 + o <=> c2h5 + oh */
    nuki[ 22 * kd + 156 ] += -1 ;
    nuki[ 11 * kd + 156 ] += -1 ;
    nuki[ 21 * kd + 156 ] += +1 ;
    nuki[ 12 * kd + 156 ] += +1 ;

    /*reaction 158: c2h6 + oh <=> c2h5 + h2o */
    nuki[ 22 * kd + 157 ] += -1 ;
    nuki[ 12 * kd + 157 ] += -1 ;
    nuki[ 21 * kd + 157 ] += +1 ;
    nuki[ 15 * kd + 157 ] += +1 ;

    /*reaction 159: c2h5 + h <=> c2h4 + h2 */
    nuki[ 21 * kd + 158 ] += -1 ;
    nuki[ 1 * kd + 158 ] += -1 ;
    nuki[ 20 * kd + 158 ] += +1 ;
    nuki[ 0 * kd + 158 ] += +1 ;

    /*reaction 160: c2h5 + h <=> ch3 + ch3 */
    nuki[ 21 * kd + 159 ] += -1 ;
    nuki[ 1 * kd + 159 ] += -1 ;
    nuki[ 3 * kd + 159 ] += +1 ;
    nuki[ 3 * kd + 159 ] += +1 ;

    /*reaction 161: c2h5 + h <=> c2h6 */
    nuki[ 21 * kd + 160 ] += -1 ;
    nuki[ 1 * kd + 160 ] += -1 ;
    nuki[ 22 * kd + 160 ] += +1 ;

    /*reaction 162: c2h5 + oh <=> c2h4 + h2o */
    nuki[ 21 * kd + 161 ] += -1 ;
    nuki[ 12 * kd + 161 ] += -1 ;
    nuki[ 20 * kd + 161 ] += +1 ;
    nuki[ 15 * kd + 161 ] += +1 ;

    /*reaction 163: c2h5 + o <=> ch3 + ch2o */
    nuki[ 21 * kd + 162 ] += -1 ;
    nuki[ 11 * kd + 162 ] += -1 ;
    nuki[ 3 * kd + 162 ] += +1 ;
    nuki[ 6 * kd + 162 ] += +1 ;

    /*reaction 164: c2h5 + ho2 <=> c2h6 + o2 */
    nuki[ 21 * kd + 163 ] += -1 ;
    nuki[ 13 * kd + 163 ] += -1 ;
    nuki[ 22 * kd + 163 ] += +1 ;
    nuki[ 10 * kd + 163 ] += +1 ;

    /*reaction 165: c2h5 + ho2 <=> ch3ch2o + oh */
    nuki[ 21 * kd + 164 ] += -1 ;
    nuki[ 13 * kd + 164 ] += -1 ;
    nuki[ 36 * kd + 164 ] += +1 ;
    nuki[ 12 * kd + 164 ] += +1 ;

    /*reaction 166: c2h5 + o2 <=> c2h4 + ho2 */
    nuki[ 21 * kd + 165 ] += -1 ;
    nuki[ 10 * kd + 165 ] += -1 ;
    nuki[ 20 * kd + 165 ] += +1 ;
    nuki[ 13 * kd + 165 ] += +1 ;

    /*reaction 167: c2h5 + o2 <=> ch3hco + oh */
    nuki[ 21 * kd + 166 ] += -1 ;
    nuki[ 10 * kd + 166 ] += -1 ;
    nuki[ 32 * kd + 166 ] += +1 ;
    nuki[ 12 * kd + 166 ] += +1 ;

    /*reaction 168: c2h4 + oh <=> c2h4oh */
    nuki[ 20 * kd + 167 ] += -1 ;
    nuki[ 12 * kd + 167 ] += -1 ;
    nuki[ 34 * kd + 167 ] += +1 ;

    /*reaction 169: c2h4 + oh <=> c2h3 + h2o */
    nuki[ 20 * kd + 168 ] += -1 ;
    nuki[ 12 * kd + 168 ] += -1 ;
    nuki[ 19 * kd + 168 ] += +1 ;
    nuki[ 15 * kd + 168 ] += +1 ;

    /*reaction 170: c2h4 + o <=> ch3 + hco */
    nuki[ 20 * kd + 169 ] += -1 ;
    nuki[ 11 * kd + 169 ] += -1 ;
    nuki[ 3 * kd + 169 ] += +1 ;
    nuki[ 7 * kd + 169 ] += +1 ;

    /*reaction 171: c2h4 + o <=> ch2hco + h */
    nuki[ 20 * kd + 170 ] += -1 ;
    nuki[ 11 * kd + 170 ] += -1 ;
    nuki[ 30 * kd + 170 ] += +1 ;
    nuki[ 1 * kd + 170 ] += +1 ;

    /*reaction 172: c2h4 + ch3 <=> c2h3 + ch4 */
    nuki[ 20 * kd + 171 ] += -1 ;
    nuki[ 3 * kd + 171 ] += -1 ;
    nuki[ 19 * kd + 171 ] += +1 ;
    nuki[ 2 * kd + 171 ] += +1 ;

    /*reaction 173: c2h4 + h <=> c2h3 + h2 */
    nuki[ 20 * kd + 172 ] += -1 ;
    nuki[ 1 * kd + 172 ] += -1 ;
    nuki[ 19 * kd + 172 ] += +1 ;
    nuki[ 0 * kd + 172 ] += +1 ;

    /*reaction 174: c2h4 + h (+M) <=> c2h5 (+M) */
    nuki[ 20 * kd + 173 ] += -1 ;
    nuki[ 1 * kd + 173 ] += -1 ;
    nuki[ 21 * kd + 173 ] += +1 ;

    /*reaction 175: c2h4 (+M) <=> c2h2 + h2 (+M) */
    nuki[ 20 * kd + 174 ] += -1 ;
    nuki[ 18 * kd + 174 ] += +1 ;
    nuki[ 0 * kd + 174 ] += +1 ;

    /*reaction 176: c2h3 + h (+M) <=> c2h4 (+M) */
    nuki[ 19 * kd + 175 ] += -1 ;
    nuki[ 1 * kd + 175 ] += -1 ;
    nuki[ 20 * kd + 175 ] += +1 ;

    /*reaction 177: c2h3 + h <=> c2h2 + h2 */
    nuki[ 19 * kd + 176 ] += -1 ;
    nuki[ 1 * kd + 176 ] += -1 ;
    nuki[ 18 * kd + 176 ] += +1 ;
    nuki[ 0 * kd + 176 ] += +1 ;

    /*reaction 178: c2h3 + o <=> ch2co + h */
    nuki[ 19 * kd + 177 ] += -1 ;
    nuki[ 11 * kd + 177 ] += -1 ;
    nuki[ 29 * kd + 177 ] += +1 ;
    nuki[ 1 * kd + 177 ] += +1 ;

    /*reaction 179: c2h3 + o2 <=> ch2o + hco */
    nuki[ 19 * kd + 178 ] += -1 ;
    nuki[ 10 * kd + 178 ] += -1 ;
    nuki[ 6 * kd + 178 ] += +1 ;
    nuki[ 7 * kd + 178 ] += +1 ;

    /*reaction 180: c2h3 + o2 <=> ch2hco + o */
    nuki[ 19 * kd + 179 ] += -1 ;
    nuki[ 10 * kd + 179 ] += -1 ;
    nuki[ 30 * kd + 179 ] += +1 ;
    nuki[ 11 * kd + 179 ] += +1 ;

    /*reaction 181: c2h3 + o2 <=> c2h2 + ho2 */
    nuki[ 19 * kd + 180 ] += -1 ;
    nuki[ 10 * kd + 180 ] += -1 ;
    nuki[ 18 * kd + 180 ] += +1 ;
    nuki[ 13 * kd + 180 ] += +1 ;

    /*reaction 182: c2h3 + oh <=> c2h2 + h2o */
    nuki[ 19 * kd + 181 ] += -1 ;
    nuki[ 12 * kd + 181 ] += -1 ;
    nuki[ 18 * kd + 181 ] += +1 ;
    nuki[ 15 * kd + 181 ] += +1 ;

    /*reaction 183: c2h3 + c2h <=> c2h2 + c2h2 */
    nuki[ 19 * kd + 182 ] += -1 ;
    nuki[ 16 * kd + 182 ] += -1 ;
    nuki[ 18 * kd + 182 ] += +1 ;
    nuki[ 18 * kd + 182 ] += +1 ;

    /*reaction 184: c2h3 + ch <=> ch2 + c2h2 */
    nuki[ 19 * kd + 183 ] += -1 ;
    nuki[ 5 * kd + 183 ] += -1 ;
    nuki[ 4 * kd + 183 ] += +1 ;
    nuki[ 18 * kd + 183 ] += +1 ;

    /*reaction 185: c2h3 + ch3 <=> c2h2 + ch4 */
    nuki[ 19 * kd + 184 ] += -1 ;
    nuki[ 3 * kd + 184 ] += -1 ;
    nuki[ 18 * kd + 184 ] += +1 ;
    nuki[ 2 * kd + 184 ] += +1 ;

    /*reaction 186: c2h2 + oh <=> c2h + h2o */
    nuki[ 18 * kd + 185 ] += -1 ;
    nuki[ 12 * kd + 185 ] += -1 ;
    nuki[ 16 * kd + 185 ] += +1 ;
    nuki[ 15 * kd + 185 ] += +1 ;

    /*reaction 187: c2h2 + oh <=> hccoh + h */
    nuki[ 18 * kd + 186 ] += -1 ;
    nuki[ 12 * kd + 186 ] += -1 ;
    nuki[ 25 * kd + 186 ] += +1 ;
    nuki[ 1 * kd + 186 ] += +1 ;

    /*reaction 188: c2h2 + oh <=> ch2co + h */
    nuki[ 18 * kd + 187 ] += -1 ;
    nuki[ 12 * kd + 187 ] += -1 ;
    nuki[ 29 * kd + 187 ] += +1 ;
    nuki[ 1 * kd + 187 ] += +1 ;

    /*reaction 189: c2h2 + oh <=> ch2co + h */
    nuki[ 18 * kd + 188 ] += -1 ;
    nuki[ 12 * kd + 188 ] += -1 ;
    nuki[ 29 * kd + 188 ] += +1 ;
    nuki[ 1 * kd + 188 ] += +1 ;

    /*reaction 190: c2h2 + oh <=> ch3 + co */
    nuki[ 18 * kd + 189 ] += -1 ;
    nuki[ 12 * kd + 189 ] += -1 ;
    nuki[ 3 * kd + 189 ] += +1 ;
    nuki[ 9 * kd + 189 ] += +1 ;

    /*reaction 191: hccoh + h <=> ch2co + h */
    nuki[ 25 * kd + 190 ] += -1 ;
    nuki[ 1 * kd + 190 ] += -1 ;
    nuki[ 29 * kd + 190 ] += +1 ;
    nuki[ 1 * kd + 190 ] += +1 ;

    /*reaction 192: c2h2 + o <=> ch2 + co */
    nuki[ 18 * kd + 191 ] += -1 ;
    nuki[ 11 * kd + 191 ] += -1 ;
    nuki[ 4 * kd + 191 ] += +1 ;
    nuki[ 9 * kd + 191 ] += +1 ;

    /*reaction 193: c2h2 + o <=> hcco + h */
    nuki[ 18 * kd + 192 ] += -1 ;
    nuki[ 11 * kd + 192 ] += -1 ;
    nuki[ 17 * kd + 192 ] += +1 ;
    nuki[ 1 * kd + 192 ] += +1 ;

    /*reaction 194: c2h2 + o <=> c2h + oh */
    nuki[ 18 * kd + 193 ] += -1 ;
    nuki[ 11 * kd + 193 ] += -1 ;
    nuki[ 16 * kd + 193 ] += +1 ;
    nuki[ 12 * kd + 193 ] += +1 ;

    /*reaction 195: c2h2 + ch3 <=> c2h + ch4 */
    nuki[ 18 * kd + 194 ] += -1 ;
    nuki[ 3 * kd + 194 ] += -1 ;
    nuki[ 16 * kd + 194 ] += +1 ;
    nuki[ 2 * kd + 194 ] += +1 ;

    /*reaction 196: c2h2 + o2 <=> hcco + oh */
    nuki[ 18 * kd + 195 ] += -1 ;
    nuki[ 10 * kd + 195 ] += -1 ;
    nuki[ 17 * kd + 195 ] += +1 ;
    nuki[ 12 * kd + 195 ] += +1 ;

    /*reaction 197: c2h2 + M <=> c2h + h + M */
    nuki[ 18 * kd + 196 ] += -1 ;
    nuki[ 16 * kd + 196 ] += +1 ;
    nuki[ 1 * kd + 196 ] += +1 ;

    /*reaction 198: c2h2 + h (+M) <=> c2h3 (+M) */
    nuki[ 18 * kd + 197 ] += -1 ;
    nuki[ 1 * kd + 197 ] += -1 ;
    nuki[ 19 * kd + 197 ] += +1 ;

    /*reaction 199: ch2hco + h <=> ch3 + hco */
    nuki[ 30 * kd + 198 ] += -1 ;
    nuki[ 1 * kd + 198 ] += -1 ;
    nuki[ 3 * kd + 198 ] += +1 ;
    nuki[ 7 * kd + 198 ] += +1 ;

    /*reaction 200: ch2hco + h <=> ch2co + h2 */
    nuki[ 30 * kd + 199 ] += -1 ;
    nuki[ 1 * kd + 199 ] += -1 ;
    nuki[ 29 * kd + 199 ] += +1 ;
    nuki[ 0 * kd + 199 ] += +1 ;

    /*reaction 201: ch2hco + o <=> ch2o + hco */
    nuki[ 30 * kd + 200 ] += -1 ;
    nuki[ 11 * kd + 200 ] += -1 ;
    nuki[ 6 * kd + 200 ] += +1 ;
    nuki[ 7 * kd + 200 ] += +1 ;

    /*reaction 202: ch2hco + oh <=> ch2co + h2o */
    nuki[ 30 * kd + 201 ] += -1 ;
    nuki[ 12 * kd + 201 ] += -1 ;
    nuki[ 29 * kd + 201 ] += +1 ;
    nuki[ 15 * kd + 201 ] += +1 ;

    /*reaction 203: ch2hco + o2 <=> ch2o + co + oh */
    nuki[ 30 * kd + 202 ] += -1 ;
    nuki[ 10 * kd + 202 ] += -1 ;
    nuki[ 6 * kd + 202 ] += +1 ;
    nuki[ 9 * kd + 202 ] += +1 ;
    nuki[ 12 * kd + 202 ] += +1 ;

    /*reaction 204: ch2hco + ch3 <=> c2h5 + co + h */
    nuki[ 30 * kd + 203 ] += -1 ;
    nuki[ 3 * kd + 203 ] += -1 ;
    nuki[ 21 * kd + 203 ] += +1 ;
    nuki[ 9 * kd + 203 ] += +1 ;
    nuki[ 1 * kd + 203 ] += +1 ;

    /*reaction 205: ch2hco + ho2 <=> ch2o + hco + oh */
    nuki[ 30 * kd + 204 ] += -1 ;
    nuki[ 13 * kd + 204 ] += -1 ;
    nuki[ 6 * kd + 204 ] += +1 ;
    nuki[ 7 * kd + 204 ] += +1 ;
    nuki[ 12 * kd + 204 ] += +1 ;

    /*reaction 206: ch2hco + ho2 <=> ch3hco + o2 */
    nuki[ 30 * kd + 205 ] += -1 ;
    nuki[ 13 * kd + 205 ] += -1 ;
    nuki[ 32 * kd + 205 ] += +1 ;
    nuki[ 10 * kd + 205 ] += +1 ;

    /*reaction 207: ch2hco <=> ch3 + co */
    nuki[ 30 * kd + 206 ] += -1 ;
    nuki[ 3 * kd + 206 ] += +1 ;
    nuki[ 9 * kd + 206 ] += +1 ;

    /*reaction 208: ch2hco <=> ch2co + h */
    nuki[ 30 * kd + 207 ] += -1 ;
    nuki[ 29 * kd + 207 ] += +1 ;
    nuki[ 1 * kd + 207 ] += +1 ;

    /*reaction 209: ch3co (+M) <=> ch3 + co (+M) */
    nuki[ 31 * kd + 208 ] += -1 ;
    nuki[ 3 * kd + 208 ] += +1 ;
    nuki[ 9 * kd + 208 ] += +1 ;

    /*reaction 210: ch2co + o <=> co2 + ch2 */
    nuki[ 29 * kd + 209 ] += -1 ;
    nuki[ 11 * kd + 209 ] += -1 ;
    nuki[ 8 * kd + 209 ] += +1 ;
    nuki[ 4 * kd + 209 ] += +1 ;

    /*reaction 211: ch2co + h <=> ch3 + co */
    nuki[ 29 * kd + 210 ] += -1 ;
    nuki[ 1 * kd + 210 ] += -1 ;
    nuki[ 3 * kd + 210 ] += +1 ;
    nuki[ 9 * kd + 210 ] += +1 ;

    /*reaction 212: ch2co + h <=> hcco + h2 */
    nuki[ 29 * kd + 211 ] += -1 ;
    nuki[ 1 * kd + 211 ] += -1 ;
    nuki[ 17 * kd + 211 ] += +1 ;
    nuki[ 0 * kd + 211 ] += +1 ;

    /*reaction 213: ch2co + o <=> hcco + oh */
    nuki[ 29 * kd + 212 ] += -1 ;
    nuki[ 11 * kd + 212 ] += -1 ;
    nuki[ 17 * kd + 212 ] += +1 ;
    nuki[ 12 * kd + 212 ] += +1 ;

    /*reaction 214: ch2co + oh <=> hcco + h2o */
    nuki[ 29 * kd + 213 ] += -1 ;
    nuki[ 12 * kd + 213 ] += -1 ;
    nuki[ 17 * kd + 213 ] += +1 ;
    nuki[ 15 * kd + 213 ] += +1 ;

    /*reaction 215: ch2co + oh <=> ch2oh + co */
    nuki[ 29 * kd + 214 ] += -1 ;
    nuki[ 12 * kd + 214 ] += -1 ;
    nuki[ 23 * kd + 214 ] += +1 ;
    nuki[ 9 * kd + 214 ] += +1 ;

    /*reaction 216: ch2co (+M) <=> ch2 + co (+M) */
    nuki[ 29 * kd + 215 ] += -1 ;
    nuki[ 4 * kd + 215 ] += +1 ;
    nuki[ 9 * kd + 215 ] += +1 ;

    /*reaction 217: c2h + h2 <=> c2h2 + h */
    nuki[ 16 * kd + 216 ] += -1 ;
    nuki[ 0 * kd + 216 ] += -1 ;
    nuki[ 18 * kd + 216 ] += +1 ;
    nuki[ 1 * kd + 216 ] += +1 ;

    /*reaction 218: c2h + o <=> ch + co */
    nuki[ 16 * kd + 217 ] += -1 ;
    nuki[ 11 * kd + 217 ] += -1 ;
    nuki[ 5 * kd + 217 ] += +1 ;
    nuki[ 9 * kd + 217 ] += +1 ;

    /*reaction 219: c2h + oh <=> hcco + h */
    nuki[ 16 * kd + 218 ] += -1 ;
    nuki[ 12 * kd + 218 ] += -1 ;
    nuki[ 17 * kd + 218 ] += +1 ;
    nuki[ 1 * kd + 218 ] += +1 ;

    /*reaction 220: c2h + o2 <=> co + co + h */
    nuki[ 16 * kd + 219 ] += -1 ;
    nuki[ 10 * kd + 219 ] += -1 ;
    nuki[ 9 * kd + 219 ] += +1 ;
    nuki[ 9 * kd + 219 ] += +1 ;
    nuki[ 1 * kd + 219 ] += +1 ;

    /*reaction 221: hcco + c2h2 <=> h2ccch + co */
    nuki[ 17 * kd + 220 ] += -1 ;
    nuki[ 18 * kd + 220 ] += -1 ;
    nuki[ 26 * kd + 220 ] += +1 ;
    nuki[ 9 * kd + 220 ] += +1 ;

    /*reaction 222: hcco + h <=> ch2s + co */
    nuki[ 17 * kd + 221 ] += -1 ;
    nuki[ 1 * kd + 221 ] += -1 ;
    nuki[ 28 * kd + 221 ] += +1 ;
    nuki[ 9 * kd + 221 ] += +1 ;

    /*reaction 223: hcco + o <=> h + co + co */
    nuki[ 17 * kd + 222 ] += -1 ;
    nuki[ 11 * kd + 222 ] += -1 ;
    nuki[ 1 * kd + 222 ] += +1 ;
    nuki[ 9 * kd + 222 ] += +1 ;
    nuki[ 9 * kd + 222 ] += +1 ;

    /*reaction 224: hcco + o <=> ch + co2 */
    nuki[ 17 * kd + 223 ] += -1 ;
    nuki[ 11 * kd + 223 ] += -1 ;
    nuki[ 5 * kd + 223 ] += +1 ;
    nuki[ 8 * kd + 223 ] += +1 ;

    /*reaction 225: hcco + o2 <=> hco + co + o */
    nuki[ 17 * kd + 224 ] += -1 ;
    nuki[ 10 * kd + 224 ] += -1 ;
    nuki[ 7 * kd + 224 ] += +1 ;
    nuki[ 9 * kd + 224 ] += +1 ;
    nuki[ 11 * kd + 224 ] += +1 ;

    /*reaction 226: hcco + o2 <=> co2 + hco */
    nuki[ 17 * kd + 225 ] += -1 ;
    nuki[ 10 * kd + 225 ] += -1 ;
    nuki[ 8 * kd + 225 ] += +1 ;
    nuki[ 7 * kd + 225 ] += +1 ;

    /*reaction 227: hcco + ch <=> c2h2 + co */
    nuki[ 17 * kd + 226 ] += -1 ;
    nuki[ 5 * kd + 226 ] += -1 ;
    nuki[ 18 * kd + 226 ] += +1 ;
    nuki[ 9 * kd + 226 ] += +1 ;

    /*reaction 228: hcco + hcco <=> c2h2 + co + co */
    nuki[ 17 * kd + 227 ] += -1 ;
    nuki[ 17 * kd + 227 ] += -1 ;
    nuki[ 18 * kd + 227 ] += +1 ;
    nuki[ 9 * kd + 227 ] += +1 ;
    nuki[ 9 * kd + 227 ] += +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < 4 * 38; ++ id) {
         ncf[id] = 0; 
    }

    /*h2 */
    ncf[ 0 * kd + 0 ] = 2; /*h */

    /*h */
    ncf[ 1 * kd + 0 ] = 1; /*h */

    /*ch4 */
    ncf[ 2 * kd + 2 ] = 1; /*c */
    ncf[ 2 * kd + 0 ] = 4; /*h */

    /*ch3 */
    ncf[ 3 * kd + 2 ] = 1; /*c */
    ncf[ 3 * kd + 0 ] = 3; /*h */

    /*ch2 */
    ncf[ 4 * kd + 2 ] = 1; /*c */
    ncf[ 4 * kd + 0 ] = 2; /*h */

    /*ch */
    ncf[ 5 * kd + 2 ] = 1; /*c */
    ncf[ 5 * kd + 0 ] = 1; /*h */

    /*ch2o */
    ncf[ 6 * kd + 2 ] = 1; /*c */
    ncf[ 6 * kd + 0 ] = 2; /*h */
    ncf[ 6 * kd + 1 ] = 1; /*o */

    /*hco */
    ncf[ 7 * kd + 0 ] = 1; /*h */
    ncf[ 7 * kd + 2 ] = 1; /*c */
    ncf[ 7 * kd + 1 ] = 1; /*o */

    /*co2 */
    ncf[ 8 * kd + 2 ] = 1; /*c */
    ncf[ 8 * kd + 1 ] = 2; /*o */

    /*co */
    ncf[ 9 * kd + 2 ] = 1; /*c */
    ncf[ 9 * kd + 1 ] = 1; /*o */

    /*o2 */
    ncf[ 10 * kd + 1 ] = 2; /*o */

    /*o */
    ncf[ 11 * kd + 1 ] = 1; /*o */

    /*oh */
    ncf[ 12 * kd + 1 ] = 1; /*o */
    ncf[ 12 * kd + 0 ] = 1; /*h */

    /*ho2 */
    ncf[ 13 * kd + 0 ] = 1; /*h */
    ncf[ 13 * kd + 1 ] = 2; /*o */

    /*h2o2 */
    ncf[ 14 * kd + 0 ] = 2; /*h */
    ncf[ 14 * kd + 1 ] = 2; /*o */

    /*h2o */
    ncf[ 15 * kd + 0 ] = 2; /*h */
    ncf[ 15 * kd + 1 ] = 1; /*o */

    /*c2h */
    ncf[ 16 * kd + 2 ] = 2; /*c */
    ncf[ 16 * kd + 0 ] = 1; /*h */

    /*hcco */
    ncf[ 17 * kd + 0 ] = 1; /*h */
    ncf[ 17 * kd + 2 ] = 2; /*c */
    ncf[ 17 * kd + 1 ] = 1; /*o */

    /*c2h2 */
    ncf[ 18 * kd + 2 ] = 2; /*c */
    ncf[ 18 * kd + 0 ] = 2; /*h */

    /*c2h3 */
    ncf[ 19 * kd + 2 ] = 2; /*c */
    ncf[ 19 * kd + 0 ] = 3; /*h */

    /*c2h4 */
    ncf[ 20 * kd + 2 ] = 2; /*c */
    ncf[ 20 * kd + 0 ] = 4; /*h */

    /*c2h5 */
    ncf[ 21 * kd + 2 ] = 2; /*c */
    ncf[ 21 * kd + 0 ] = 5; /*h */

    /*c2h6 */
    ncf[ 22 * kd + 2 ] = 2; /*c */
    ncf[ 22 * kd + 0 ] = 6; /*h */

    /*ch2oh */
    ncf[ 23 * kd + 2 ] = 1; /*c */
    ncf[ 23 * kd + 0 ] = 3; /*h */
    ncf[ 23 * kd + 1 ] = 1; /*o */

    /*ch3o */
    ncf[ 24 * kd + 2 ] = 1; /*c */
    ncf[ 24 * kd + 0 ] = 3; /*h */
    ncf[ 24 * kd + 1 ] = 1; /*o */

    /*hccoh */
    ncf[ 25 * kd + 0 ] = 2; /*h */
    ncf[ 25 * kd + 2 ] = 2; /*c */
    ncf[ 25 * kd + 1 ] = 1; /*o */

    /*h2ccch */
    ncf[ 26 * kd + 2 ] = 3; /*c */
    ncf[ 26 * kd + 0 ] = 3; /*h */

    /*c3h2 */
    ncf[ 27 * kd + 0 ] = 2; /*h */
    ncf[ 27 * kd + 2 ] = 3; /*c */

    /*ch2s */
    ncf[ 28 * kd + 2 ] = 1; /*c */
    ncf[ 28 * kd + 0 ] = 2; /*h */

    /*ch2co */
    ncf[ 29 * kd + 2 ] = 2; /*c */
    ncf[ 29 * kd + 0 ] = 2; /*h */
    ncf[ 29 * kd + 1 ] = 1; /*o */

    /*ch2hco */
    ncf[ 30 * kd + 1 ] = 1; /*o */
    ncf[ 30 * kd + 0 ] = 3; /*h */
    ncf[ 30 * kd + 2 ] = 2; /*c */

    /*ch3co */
    ncf[ 31 * kd + 2 ] = 2; /*c */
    ncf[ 31 * kd + 0 ] = 3; /*h */
    ncf[ 31 * kd + 1 ] = 1; /*o */

    /*ch3hco */
    ncf[ 32 * kd + 2 ] = 2; /*c */
    ncf[ 32 * kd + 1 ] = 1; /*o */
    ncf[ 32 * kd + 0 ] = 4; /*h */

    /*c2h5oh */
    ncf[ 33 * kd + 2 ] = 2; /*c */
    ncf[ 33 * kd + 0 ] = 6; /*h */
    ncf[ 33 * kd + 1 ] = 1; /*o */

    /*c2h4oh */
    ncf[ 34 * kd + 2 ] = 2; /*c */
    ncf[ 34 * kd + 0 ] = 5; /*h */
    ncf[ 34 * kd + 1 ] = 1; /*o */

    /*ch3choh */
    ncf[ 35 * kd + 2 ] = 2; /*c */
    ncf[ 35 * kd + 0 ] = 5; /*h */
    ncf[ 35 * kd + 1 ] = 1; /*o */

    /*ch3ch2o */
    ncf[ 36 * kd + 2 ] = 2; /*c */
    ncf[ 36 * kd + 0 ] = 5; /*h */
    ncf[ 36 * kd + 1 ] = 1; /*o */

    /*n2 */
    ncf[ 37 * kd + 3 ] = 2; /*n */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE(int * iwrk, double * rwrk, double * a, double * b, double * e)
{

    /*reaction 1: oh + h2 <=> h + h2o */
    a[0] = 2.14e+08;
    b[0] = 1.52;
    e[0] = 3449;

    /*reaction 2: o + oh <=> o2 + h */
    a[1] = 2.02e+14;
    b[1] = -0.4;
    e[1] = 0;

    /*reaction 3: o + h2 <=> oh + h */
    a[2] = 50600;
    b[2] = 2.67;
    e[2] = 6290;

    /*reaction 4: h + o2 (+M) <=> ho2 (+M) */
    a[3] = 4.52e+13;
    b[3] = 0;
    e[3] = 0;

    /*reaction 5: h + o2 (+n2) <=> ho2 (+n2) */
    a[4] = 4.52e+13;
    b[4] = 0;
    e[4] = 0;

    /*reaction 6: h + o2 (+h2) <=> ho2 (+h2) */
    a[5] = 4.52e+13;
    b[5] = 0;
    e[5] = 0;

    /*reaction 7: h + o2 (+h2o) <=> ho2 (+h2o) */
    a[6] = 4.52e+13;
    b[6] = 0;
    e[6] = 0;

    /*reaction 8: oh + ho2 <=> h2o + o2 */
    a[7] = 2.13e+28;
    b[7] = -4.827;
    e[7] = 3500;

    /*reaction 9: oh + ho2 <=> h2o + o2 */
    a[8] = 9.1e+14;
    b[8] = 0;
    e[8] = 10964;

    /*reaction 10: h + ho2 <=> oh + oh */
    a[9] = 1.5e+14;
    b[9] = 0;
    e[9] = 1000;

    /*reaction 11: h + ho2 <=> h2 + o2 */
    a[10] = 6.63e+13;
    b[10] = 0;
    e[10] = 2126;

    /*reaction 12: h + ho2 <=> o + h2o */
    a[11] = 3.01e+13;
    b[11] = 0;
    e[11] = 1721;

    /*reaction 13: o + ho2 <=> o2 + oh */
    a[12] = 3.25e+13;
    b[12] = 0;
    e[12] = 0;

    /*reaction 14: oh + oh <=> o + h2o */
    a[13] = 35700;
    b[13] = 2.4;
    e[13] = -2112;

    /*reaction 15: h + h + M <=> h2 + M */
    a[14] = 1e+18;
    b[14] = -1;
    e[14] = 0;

    /*reaction 16: h + h + h2 <=> h2 + h2 */
    a[15] = 9.2e+16;
    b[15] = -0.6;
    e[15] = 0;

    /*reaction 17: h + h + h2o <=> h2 + h2o */
    a[16] = 6e+19;
    b[16] = -1.25;
    e[16] = 0;

    /*reaction 18: h + oh + M <=> h2o + M */
    a[17] = 2.21e+22;
    b[17] = -2;
    e[17] = 0;

    /*reaction 19: h + o + M <=> oh + M */
    a[18] = 4.71e+18;
    b[18] = -1;
    e[18] = 0;

    /*reaction 20: o + o + M <=> o2 + M */
    a[19] = 1.89e+13;
    b[19] = 0;
    e[19] = -1788;

    /*reaction 21: ho2 + ho2 <=> h2o2 + o2 */
    a[20] = 4.2e+14;
    b[20] = 0;
    e[20] = 11982;

    /*reaction 22: ho2 + ho2 <=> h2o2 + o2 */
    a[21] = 1.3e+11;
    b[21] = 0;
    e[21] = -1629;

    /*reaction 23: oh + oh (+M) <=> h2o2 (+M) */
    a[22] = 1.24e+14;
    b[22] = -0.37;
    e[22] = 0;

    /*reaction 24: h2o2 + h <=> ho2 + h2 */
    a[23] = 1.98e+06;
    b[23] = 2;
    e[23] = 2435;

    /*reaction 25: h2o2 + h <=> oh + h2o */
    a[24] = 3.07e+13;
    b[24] = 0;
    e[24] = 4217;

    /*reaction 26: h2o2 + o <=> oh + ho2 */
    a[25] = 9.55e+06;
    b[25] = 2;
    e[25] = 3970;

    /*reaction 27: h2o2 + oh <=> h2o + ho2 */
    a[26] = 2.4;
    b[26] = 4.042;
    e[26] = -2162;

    /*reaction 28: ch3 + ch3 (+M) <=> c2h6 (+M) */
    a[27] = 9.22e+16;
    b[27] = -1.174;
    e[27] = 636;

    /*reaction 29: ch3 + h (+M) <=> ch4 (+M) */
    a[28] = 2.14e+15;
    b[28] = -0.4;
    e[28] = 0;

    /*reaction 30: ch4 + h <=> ch3 + h2 */
    a[29] = 22000;
    b[29] = 3;
    e[29] = 8750;

    /*reaction 31: ch4 + oh <=> ch3 + h2o */
    a[30] = 4.19e+06;
    b[30] = 2;
    e[30] = 2547;

    /*reaction 32: ch4 + o <=> ch3 + oh */
    a[31] = 6.92e+08;
    b[31] = 1.56;
    e[31] = 8485;

    /*reaction 33: ch4 + ho2 <=> ch3 + h2o2 */
    a[32] = 1.12e+13;
    b[32] = 0;
    e[32] = 24640;

    /*reaction 34: ch3 + ho2 <=> ch3o + oh */
    a[33] = 7e+12;
    b[33] = 0;
    e[33] = 0;

    /*reaction 35: ch3 + ho2 <=> ch4 + o2 */
    a[34] = 3e+12;
    b[34] = 0;
    e[34] = 0;

    /*reaction 36: ch3 + o <=> ch2o + h */
    a[35] = 8e+13;
    b[35] = 0;
    e[35] = 0;

    /*reaction 37: ch3 + o2 <=> ch3o + o */
    a[36] = 1.45e+13;
    b[36] = 0;
    e[36] = 29209;

    /*reaction 38: ch3 + o2 <=> ch2o + oh */
    a[37] = 2.51e+11;
    b[37] = 0;
    e[37] = 14640;

    /*reaction 39: ch3o + h <=> ch3 + oh */
    a[38] = 1e+13;
    b[38] = 0;
    e[38] = 0;

    /*reaction 40: ch2oh + h <=> ch3 + oh */
    a[39] = 1e+13;
    b[39] = 0;
    e[39] = 0;

    /*reaction 41: ch3 + oh <=> ch2s + h2o */
    a[40] = 2e+13;
    b[40] = 0;
    e[40] = 550;

    /*reaction 42: ch3 + oh <=> ch2 + h2o */
    a[41] = 3e+06;
    b[41] = 2;
    e[41] = 2500;

    /*reaction 43: ch3 + h <=> ch2 + h2 */
    a[42] = 9e+13;
    b[42] = 0;
    e[42] = 15100;

    /*reaction 44: ch3 + M <=> ch + h2 + M */
    a[43] = 6.9e+14;
    b[43] = 0;
    e[43] = 82469;

    /*reaction 45: ch3 + M <=> ch2 + h + M */
    a[44] = 1.9e+16;
    b[44] = 0;
    e[44] = 91411;

    /*reaction 46: ch2o + h (+M) <=> ch3o (+M) */
    a[45] = 5.4e+11;
    b[45] = 0.454;
    e[45] = 2600;

    /*reaction 47: ch2o + h (+M) <=> ch2oh (+M) */
    a[46] = 5.4e+11;
    b[46] = 0.454;
    e[46] = 3600;

    /*reaction 48: ch3o + ch3 <=> ch2o + ch4 */
    a[47] = 1.2e+13;
    b[47] = 0;
    e[47] = 0;

    /*reaction 49: ch3o + h <=> ch2o + h2 */
    a[48] = 2e+13;
    b[48] = 0;
    e[48] = 0;

    /*reaction 50: ch2oh + h <=> ch2o + h2 */
    a[49] = 2e+13;
    b[49] = 0;
    e[49] = 0;

    /*reaction 51: ch3o + oh <=> ch2o + h2o */
    a[50] = 1e+13;
    b[50] = 0;
    e[50] = 0;

    /*reaction 52: ch2oh + oh <=> ch2o + h2o */
    a[51] = 1e+13;
    b[51] = 0;
    e[51] = 0;

    /*reaction 53: ch3o + o <=> ch2o + oh */
    a[52] = 1e+13;
    b[52] = 0;
    e[52] = 0;

    /*reaction 54: ch2oh + o <=> ch2o + oh */
    a[53] = 1e+13;
    b[53] = 0;
    e[53] = 0;

    /*reaction 55: ch3o + o2 <=> ch2o + ho2 */
    a[54] = 6.3e+10;
    b[54] = 0;
    e[54] = 2600;

    /*reaction 56: ch3o + co <=> ch3 + co2 */
    a[55] = 468;
    b[55] = 3.16;
    e[55] = 5380;

    /*reaction 57: ch2oh + o2 <=> ch2o + ho2 */
    a[56] = 1.57e+15;
    b[56] = -1;
    e[56] = 0;

    /*reaction 58: ch2oh + o2 <=> ch2o + ho2 */
    a[57] = 7.23e+13;
    b[57] = 0;
    e[57] = 3577;

    /*reaction 59: ch2 + h <=> ch + h2 */
    a[58] = 1e+18;
    b[58] = -1.56;
    e[58] = 0;

    /*reaction 60: ch2 + oh <=> ch + h2o */
    a[59] = 1.13e+07;
    b[59] = 2;
    e[59] = 3000;

    /*reaction 61: ch2 + oh <=> ch2o + h */
    a[60] = 2.5e+13;
    b[60] = 0;
    e[60] = 0;

    /*reaction 62: ch2 + co2 <=> ch2o + co */
    a[61] = 1.1e+11;
    b[61] = 0;
    e[61] = 1000;

    /*reaction 63: ch2 + o <=> co + h + h */
    a[62] = 5e+13;
    b[62] = 0;
    e[62] = 0;

    /*reaction 64: ch2 + o <=> co + h2 */
    a[63] = 3e+13;
    b[63] = 0;
    e[63] = 0;

    /*reaction 65: ch2 + o2 <=> ch2o + o */
    a[64] = 3.29e+21;
    b[64] = -3.3;
    e[64] = 2868;

    /*reaction 66: ch2 + o2 <=> co2 + h + h */
    a[65] = 3.29e+21;
    b[65] = -3.3;
    e[65] = 2868;

    /*reaction 67: ch2 + o2 <=> co2 + h2 */
    a[66] = 1.01e+21;
    b[66] = -3.3;
    e[66] = 1508;

    /*reaction 68: ch2 + o2 <=> co + h2o */
    a[67] = 7.28e+19;
    b[67] = -2.54;
    e[67] = 1809;

    /*reaction 69: ch2 + o2 <=> hco + oh */
    a[68] = 1.29e+20;
    b[68] = -3.3;
    e[68] = 284;

    /*reaction 70: ch2 + ch3 <=> c2h4 + h */
    a[69] = 4e+13;
    b[69] = 0;
    e[69] = 0;

    /*reaction 71: ch2 + ch2 <=> c2h2 + h + h */
    a[70] = 4e+13;
    b[70] = 0;
    e[70] = 0;

    /*reaction 72: ch2 + hcco <=> c2h3 + co */
    a[71] = 3e+13;
    b[71] = 0;
    e[71] = 0;

    /*reaction 73: ch2 + c2h2 <=> h2ccch + h */
    a[72] = 1.2e+13;
    b[72] = 0;
    e[72] = 6600;

    /*reaction 74: ch2s + M <=> ch2 + M */
    a[73] = 1e+13;
    b[73] = 0;
    e[73] = 0;

    /*reaction 75: ch2s + ch4 <=> ch3 + ch3 */
    a[74] = 4e+13;
    b[74] = 0;
    e[74] = 0;

    /*reaction 76: ch2s + c2h6 <=> ch3 + c2h5 */
    a[75] = 1.2e+14;
    b[75] = 0;
    e[75] = 0;

    /*reaction 77: ch2s + o2 <=> co + oh + h */
    a[76] = 7e+13;
    b[76] = 0;
    e[76] = 0;

    /*reaction 78: ch2s + h2 <=> ch3 + h */
    a[77] = 7e+13;
    b[77] = 0;
    e[77] = 0;

    /*reaction 79: ch2s + c2h2 <=> h2ccch + h */
    a[78] = 1.5e+14;
    b[78] = 0;
    e[78] = 0;

    /*reaction 80: ch2s + o <=> co + h + h */
    a[79] = 3e+13;
    b[79] = 0;
    e[79] = 0;

    /*reaction 81: ch2s + oh <=> ch2o + h */
    a[80] = 3e+13;
    b[80] = 0;
    e[80] = 0;

    /*reaction 82: ch2s + h <=> ch + h2 */
    a[81] = 3e+13;
    b[81] = 0;
    e[81] = 0;

    /*reaction 83: ch2s + co2 <=> ch2o + co */
    a[82] = 3e+12;
    b[82] = 0;
    e[82] = 0;

    /*reaction 84: ch2s + ch3 <=> c2h4 + h */
    a[83] = 2e+13;
    b[83] = 0;
    e[83] = 0;

    /*reaction 85: ch2s + ch2co <=> c2h4 + co */
    a[84] = 1.6e+14;
    b[84] = 0;
    e[84] = 0;

    /*reaction 86: ch + o2 <=> hco + o */
    a[85] = 3.3e+13;
    b[85] = 0;
    e[85] = 0;

    /*reaction 87: ch + o <=> co + h */
    a[86] = 5.7e+13;
    b[86] = 0;
    e[86] = 0;

    /*reaction 88: ch + oh <=> hco + h */
    a[87] = 3e+13;
    b[87] = 0;
    e[87] = 0;

    /*reaction 89: ch + co2 <=> hco + co */
    a[88] = 3.4e+12;
    b[88] = 0;
    e[88] = 690;

    /*reaction 90: ch + h2o <=> ch2o + h */
    a[89] = 1.17e+15;
    b[89] = -0.75;
    e[89] = 0;

    /*reaction 91: ch + ch2o <=> ch2co + h */
    a[90] = 9.46e+13;
    b[90] = 0;
    e[90] = -515;

    /*reaction 92: ch + c2h2 <=> c3h2 + h */
    a[91] = 1e+14;
    b[91] = 0;
    e[91] = 0;

    /*reaction 93: ch + ch2 <=> c2h2 + h */
    a[92] = 4e+13;
    b[92] = 0;
    e[92] = 0;

    /*reaction 94: ch + ch3 <=> c2h3 + h */
    a[93] = 3e+13;
    b[93] = 0;
    e[93] = 0;

    /*reaction 95: ch + ch4 <=> c2h4 + h */
    a[94] = 6e+13;
    b[94] = 0;
    e[94] = 0;

    /*reaction 96: ch2o + oh <=> hco + h2o */
    a[95] = 3.43e+09;
    b[95] = 1.18;
    e[95] = -447;

    /*reaction 97: ch2o + h <=> hco + h2 */
    a[96] = 2.19e+08;
    b[96] = 1.77;
    e[96] = 3000;

    /*reaction 98: ch2o + M <=> hco + h + M */
    a[97] = 3.31e+16;
    b[97] = 0;
    e[97] = 81000;

    /*reaction 99: ch2o + o <=> hco + oh */
    a[98] = 1.8e+13;
    b[98] = 0;
    e[98] = 3080;

    /*reaction 100: hco + o2 <=> co + ho2 */
    a[99] = 7.58e+12;
    b[99] = 0;
    e[99] = 410;

    /*reaction 101: hco + M <=> h + co + M */
    a[100] = 1.86e+17;
    b[100] = -1;
    e[100] = 17000;

    /*reaction 102: hco + oh <=> h2o + co */
    a[101] = 1e+14;
    b[101] = 0;
    e[101] = 0;

    /*reaction 103: hco + h <=> co + h2 */
    a[102] = 1.19e+13;
    b[102] = 0.25;
    e[102] = 0;

    /*reaction 104: hco + o <=> co + oh */
    a[103] = 3e+13;
    b[103] = 0;
    e[103] = 0;

    /*reaction 105: hco + o <=> co2 + h */
    a[104] = 3e+13;
    b[104] = 0;
    e[104] = 0;

    /*reaction 106: co + oh <=> co2 + h */
    a[105] = 9420;
    b[105] = 2.25;
    e[105] = -2351;

    /*reaction 107: co + o + M <=> co2 + M */
    a[106] = 6.17e+14;
    b[106] = 0;
    e[106] = 3000;

    /*reaction 108: co + o2 <=> co2 + o */
    a[107] = 2.53e+12;
    b[107] = 0;
    e[107] = 47688;

    /*reaction 109: co + ho2 <=> co2 + oh */
    a[108] = 5.8e+13;
    b[108] = 0;
    e[108] = 22934;

    /*reaction 110: c2h5oh (+M) <=> ch3 + ch2oh (+M) */
    a[109] = 5.94e+23;
    b[109] = -1.68;
    e[109] = 91163;

    /*reaction 111: c2h5oh (+M) <=> c2h5 + oh (+M) */
    a[110] = 1.25e+23;
    b[110] = -1.54;
    e[110] = 96005;

    /*reaction 112: c2h5oh (+M) <=> c2h4 + h2o (+M) */
    a[111] = 2.79e+13;
    b[111] = 0.09;
    e[111] = 66136;

    /*reaction 113: c2h5oh (+M) <=> ch3hco + h2 (+M) */
    a[112] = 7.24e+11;
    b[112] = 0.095;
    e[112] = 91007;

    /*reaction 114: c2h5oh + oh <=> c2h4oh + h2o */
    a[113] = 1.74e+11;
    b[113] = 0.27;
    e[113] = 600;

    /*reaction 115: c2h5oh + oh <=> ch3choh + h2o */
    a[114] = 4.64e+11;
    b[114] = 0.15;
    e[114] = 0;

    /*reaction 116: c2h5oh + oh <=> ch3ch2o + h2o */
    a[115] = 7.46e+11;
    b[115] = 0.3;
    e[115] = 1634;

    /*reaction 117: c2h5oh + h <=> c2h4oh + h2 */
    a[116] = 1.23e+07;
    b[116] = 1.8;
    e[116] = 5098;

    /*reaction 118: c2h5oh + h <=> ch3choh + h2 */
    a[117] = 2.58e+07;
    b[117] = 1.65;
    e[117] = 2827;

    /*reaction 119: c2h5oh + h <=> ch3ch2o + h2 */
    a[118] = 1.5e+07;
    b[118] = 1.6;
    e[118] = 3038;

    /*reaction 120: c2h5oh + o <=> c2h4oh + oh */
    a[119] = 9.41e+07;
    b[119] = 1.7;
    e[119] = 5459;

    /*reaction 121: c2h5oh + o <=> ch3choh + oh */
    a[120] = 1.88e+07;
    b[120] = 1.85;
    e[120] = 1824;

    /*reaction 122: c2h5oh + o <=> ch3ch2o + oh */
    a[121] = 1.58e+07;
    b[121] = 2;
    e[121] = 4448;

    /*reaction 123: c2h5oh + ch3 <=> c2h4oh + ch4 */
    a[122] = 219;
    b[122] = 3.18;
    e[122] = 9622;

    /*reaction 124: c2h5oh + ch3 <=> ch3choh + ch4 */
    a[123] = 728;
    b[123] = 2.99;
    e[123] = 7948;

    /*reaction 125: c2h5oh + ch3 <=> ch3ch2o + ch4 */
    a[124] = 145;
    b[124] = 2.99;
    e[124] = 7649;

    /*reaction 126: c2h5oh + ho2 <=> ch3choh + h2o2 */
    a[125] = 8200;
    b[125] = 2.55;
    e[125] = 10750;

    /*reaction 127: c2h5oh + ho2 <=> c2h4oh + h2o2 */
    a[126] = 12300;
    b[126] = 2.55;
    e[126] = 15750;

    /*reaction 128: c2h5oh + ho2 <=> ch3ch2o + h2o2 */
    a[127] = 2.5e+12;
    b[127] = 0;
    e[127] = 24000;

    /*reaction 129: ch3ch2o + M <=> ch3hco + h + M */
    a[128] = 1.16e+35;
    b[128] = -5.89;
    e[128] = 25274;

    /*reaction 130: ch3ch2o + M <=> ch3 + ch2o + M */
    a[129] = 1.35e+38;
    b[129] = -6.96;
    e[129] = 23800;

    /*reaction 131: ch3ch2o + o2 <=> ch3hco + ho2 */
    a[130] = 4e+10;
    b[130] = 0;
    e[130] = 1100;

    /*reaction 132: ch3ch2o + co <=> c2h5 + co2 */
    a[131] = 468;
    b[131] = 3.16;
    e[131] = 5380;

    /*reaction 133: ch3ch2o + h <=> ch3 + ch2oh */
    a[132] = 3e+13;
    b[132] = 0;
    e[132] = 0;

    /*reaction 134: ch3ch2o + h <=> c2h4 + h2o */
    a[133] = 3e+13;
    b[133] = 0;
    e[133] = 0;

    /*reaction 135: ch3ch2o + oh <=> ch3hco + h2o */
    a[134] = 1e+13;
    b[134] = 0;
    e[134] = 0;

    /*reaction 136: ch3choh + o2 <=> ch3hco + ho2 */
    a[135] = 4.82e+14;
    b[135] = 0;
    e[135] = 5017;

    /*reaction 137: ch3choh + o2 <=> ch3hco + ho2 */
    a[136] = 8.43e+15;
    b[136] = -1.2;
    e[136] = 0;

    /*reaction 138: ch3choh + o <=> ch3hco + oh */
    a[137] = 1e+14;
    b[137] = 0;
    e[137] = 0;

    /*reaction 139: ch3choh + h <=> c2h4 + h2o */
    a[138] = 3e+13;
    b[138] = 0;
    e[138] = 0;

    /*reaction 140: ch3choh + h <=> ch3 + ch2oh */
    a[139] = 3e+13;
    b[139] = 0;
    e[139] = 0;

    /*reaction 141: ch3choh + ho2 <=> ch3hco + oh + oh */
    a[140] = 4e+13;
    b[140] = 0;
    e[140] = 0;

    /*reaction 142: ch3choh + oh <=> ch3hco + h2o */
    a[141] = 5e+12;
    b[141] = 0;
    e[141] = 0;

    /*reaction 143: ch3choh + M <=> ch3hco + h + M */
    a[142] = 1e+14;
    b[142] = 0;
    e[142] = 25000;

    /*reaction 144: ch3hco + oh <=> ch3co + h2o */
    a[143] = 9.24e+06;
    b[143] = 1.5;
    e[143] = -962;

    /*reaction 145: ch3hco + oh <=> ch2hco + h2o */
    a[144] = 172000;
    b[144] = 2.4;
    e[144] = 815;

    /*reaction 146: ch3hco + o <=> ch3co + oh */
    a[145] = 1.77e+18;
    b[145] = -1.9;
    e[145] = 2975;

    /*reaction 147: ch3hco + o <=> ch2hco + oh */
    a[146] = 3.72e+13;
    b[146] = -0.2;
    e[146] = 3556;

    /*reaction 148: ch3hco + h <=> ch3co + h2 */
    a[147] = 4.66e+13;
    b[147] = -0.35;
    e[147] = 2988;

    /*reaction 149: ch3hco + h <=> ch2hco + h2 */
    a[148] = 1.85e+12;
    b[148] = 0.4;
    e[148] = 5359;

    /*reaction 150: ch3hco + ch3 <=> ch3co + ch4 */
    a[149] = 3.9e-07;
    b[149] = 5.8;
    e[149] = 2200;

    /*reaction 151: ch3hco + ch3 <=> ch2hco + ch4 */
    a[150] = 24.5;
    b[150] = 3.15;
    e[150] = 5727;

    /*reaction 152: ch3hco + ho2 <=> ch3co + h2o2 */
    a[151] = 2.4e+19;
    b[151] = -2.2;
    e[151] = 14030;

    /*reaction 153: ch3hco + ho2 <=> ch2hco + h2o2 */
    a[152] = 2.32e+11;
    b[152] = 0.4;
    e[152] = 14864;

    /*reaction 154: ch3hco + o2 <=> ch3co + ho2 */
    a[153] = 1e+14;
    b[153] = 0;
    e[153] = 42200;

    /*reaction 155: c2h6 + ch3 <=> c2h5 + ch4 */
    a[154] = 0.55;
    b[154] = 4;
    e[154] = 8300;

    /*reaction 156: c2h6 + h <=> c2h5 + h2 */
    a[155] = 540;
    b[155] = 3.5;
    e[155] = 5210;

    /*reaction 157: c2h6 + o <=> c2h5 + oh */
    a[156] = 3e+07;
    b[156] = 2;
    e[156] = 5115;

    /*reaction 158: c2h6 + oh <=> c2h5 + h2o */
    a[157] = 7.23e+06;
    b[157] = 2;
    e[157] = 864;

    /*reaction 159: c2h5 + h <=> c2h4 + h2 */
    a[158] = 1.25e+14;
    b[158] = 0;
    e[158] = 8000;

    /*reaction 160: c2h5 + h <=> ch3 + ch3 */
    a[159] = 3e+13;
    b[159] = 0;
    e[159] = 0;

    /*reaction 161: c2h5 + h <=> c2h6 */
    a[160] = 3e+13;
    b[160] = 0;
    e[160] = 0;

    /*reaction 162: c2h5 + oh <=> c2h4 + h2o */
    a[161] = 4e+13;
    b[161] = 0;
    e[161] = 0;

    /*reaction 163: c2h5 + o <=> ch3 + ch2o */
    a[162] = 1e+14;
    b[162] = 0;
    e[162] = 0;

    /*reaction 164: c2h5 + ho2 <=> c2h6 + o2 */
    a[163] = 3e+12;
    b[163] = 0;
    e[163] = 0;

    /*reaction 165: c2h5 + ho2 <=> ch3ch2o + oh */
    a[164] = 3e+13;
    b[164] = 0;
    e[164] = 0;

    /*reaction 166: c2h5 + o2 <=> c2h4 + ho2 */
    a[165] = 2.89e+28;
    b[165] = -5.4;
    e[165] = 7585;

    /*reaction 167: c2h5 + o2 <=> ch3hco + oh */
    a[166] = 4.9e+11;
    b[166] = -0.48;
    e[166] = 8357;

    /*reaction 168: c2h4 + oh <=> c2h4oh */
    a[167] = 1.29e+12;
    b[167] = 0;
    e[167] = -817;

    /*reaction 169: c2h4 + oh <=> c2h3 + h2o */
    a[168] = 2.02e+13;
    b[168] = 0;
    e[168] = 5936;

    /*reaction 170: c2h4 + o <=> ch3 + hco */
    a[169] = 1.02e+07;
    b[169] = 1.88;
    e[169] = 179;

    /*reaction 171: c2h4 + o <=> ch2hco + h */
    a[170] = 3.39e+06;
    b[170] = 1.88;
    e[170] = 179;

    /*reaction 172: c2h4 + ch3 <=> c2h3 + ch4 */
    a[171] = 6.62;
    b[171] = 3.7;
    e[171] = 9500;

    /*reaction 173: c2h4 + h <=> c2h3 + h2 */
    a[172] = 3.36e-07;
    b[172] = 6;
    e[172] = 1692;

    /*reaction 174: c2h4 + h (+M) <=> c2h5 (+M) */
    a[173] = 1.08e+12;
    b[173] = 0.454;
    e[173] = 1822;

    /*reaction 175: c2h4 (+M) <=> c2h2 + h2 (+M) */
    a[174] = 1.8e+14;
    b[174] = 0;
    e[174] = 87000;

    /*reaction 176: c2h3 + h (+M) <=> c2h4 (+M) */
    a[175] = 6.1e+12;
    b[175] = 0.27;
    e[175] = 280;

    /*reaction 177: c2h3 + h <=> c2h2 + h2 */
    a[176] = 9e+13;
    b[176] = 0;
    e[176] = 0;

    /*reaction 178: c2h3 + o <=> ch2co + h */
    a[177] = 3e+13;
    b[177] = 0;
    e[177] = 0;

    /*reaction 179: c2h3 + o2 <=> ch2o + hco */
    a[178] = 1.7e+29;
    b[178] = -5.312;
    e[178] = 6500;

    /*reaction 180: c2h3 + o2 <=> ch2hco + o */
    a[179] = 5.5e+14;
    b[179] = -0.611;
    e[179] = 5260;

    /*reaction 181: c2h3 + o2 <=> c2h2 + ho2 */
    a[180] = 2.12e-06;
    b[180] = 6;
    e[180] = 9484;

    /*reaction 182: c2h3 + oh <=> c2h2 + h2o */
    a[181] = 2e+13;
    b[181] = 0;
    e[181] = 0;

    /*reaction 183: c2h3 + c2h <=> c2h2 + c2h2 */
    a[182] = 3e+13;
    b[182] = 0;
    e[182] = 0;

    /*reaction 184: c2h3 + ch <=> ch2 + c2h2 */
    a[183] = 5e+13;
    b[183] = 0;
    e[183] = 0;

    /*reaction 185: c2h3 + ch3 <=> c2h2 + ch4 */
    a[184] = 2e+13;
    b[184] = 0;
    e[184] = 0;

    /*reaction 186: c2h2 + oh <=> c2h + h2o */
    a[185] = 3.37e+07;
    b[185] = 2;
    e[185] = 14000;

    /*reaction 187: c2h2 + oh <=> hccoh + h */
    a[186] = 504000;
    b[186] = 2.3;
    e[186] = 13500;

    /*reaction 188: c2h2 + oh <=> ch2co + h */
    a[187] = 0.000218;
    b[187] = 4.5;
    e[187] = -1000;

    /*reaction 189: c2h2 + oh <=> ch2co + h */
    a[188] = 2e+11;
    b[188] = 0;
    e[188] = 0;

    /*reaction 190: c2h2 + oh <=> ch3 + co */
    a[189] = 0.000483;
    b[189] = 4;
    e[189] = -2000;

    /*reaction 191: hccoh + h <=> ch2co + h */
    a[190] = 1e+13;
    b[190] = 0;
    e[190] = 0;

    /*reaction 192: c2h2 + o <=> ch2 + co */
    a[191] = 6.12e+06;
    b[191] = 2;
    e[191] = 1900;

    /*reaction 193: c2h2 + o <=> hcco + h */
    a[192] = 1.43e+07;
    b[192] = 2;
    e[192] = 1900;

    /*reaction 194: c2h2 + o <=> c2h + oh */
    a[193] = 3.16e+15;
    b[193] = -0.6;
    e[193] = 15000;

    /*reaction 195: c2h2 + ch3 <=> c2h + ch4 */
    a[194] = 1.81e+11;
    b[194] = 0;
    e[194] = 17289;

    /*reaction 196: c2h2 + o2 <=> hcco + oh */
    a[195] = 4e+07;
    b[195] = 1.5;
    e[195] = 30100;

    /*reaction 197: c2h2 + M <=> c2h + h + M */
    a[196] = 4.2e+16;
    b[196] = 0;
    e[196] = 107000;

    /*reaction 198: c2h2 + h (+M) <=> c2h3 (+M) */
    a[197] = 3.11e+11;
    b[197] = 0.58;
    e[197] = 2589;

    /*reaction 199: ch2hco + h <=> ch3 + hco */
    a[198] = 5e+13;
    b[198] = 0;
    e[198] = 0;

    /*reaction 200: ch2hco + h <=> ch2co + h2 */
    a[199] = 2e+13;
    b[199] = 0;
    e[199] = 0;

    /*reaction 201: ch2hco + o <=> ch2o + hco */
    a[200] = 1e+14;
    b[200] = 0;
    e[200] = 0;

    /*reaction 202: ch2hco + oh <=> ch2co + h2o */
    a[201] = 3e+13;
    b[201] = 0;
    e[201] = 0;

    /*reaction 203: ch2hco + o2 <=> ch2o + co + oh */
    a[202] = 3e+10;
    b[202] = 0;
    e[202] = 0;

    /*reaction 204: ch2hco + ch3 <=> c2h5 + co + h */
    a[203] = 4.9e+14;
    b[203] = -0.5;
    e[203] = 0;

    /*reaction 205: ch2hco + ho2 <=> ch2o + hco + oh */
    a[204] = 7e+12;
    b[204] = 0;
    e[204] = 0;

    /*reaction 206: ch2hco + ho2 <=> ch3hco + o2 */
    a[205] = 3e+12;
    b[205] = 0;
    e[205] = 0;

    /*reaction 207: ch2hco <=> ch3 + co */
    a[206] = 1.17e+43;
    b[206] = -9.83;
    e[206] = 43756;

    /*reaction 208: ch2hco <=> ch2co + h */
    a[207] = 1.81e+43;
    b[207] = -9.61;
    e[207] = 45868;

    /*reaction 209: ch3co (+M) <=> ch3 + co (+M) */
    a[208] = 3e+12;
    b[208] = 0;
    e[208] = 16722;

    /*reaction 210: ch2co + o <=> co2 + ch2 */
    a[209] = 1.75e+12;
    b[209] = 0;
    e[209] = 1350;

    /*reaction 211: ch2co + h <=> ch3 + co */
    a[210] = 27100;
    b[210] = 2.75;
    e[210] = 714;

    /*reaction 212: ch2co + h <=> hcco + h2 */
    a[211] = 2e+14;
    b[211] = 0;
    e[211] = 8000;

    /*reaction 213: ch2co + o <=> hcco + oh */
    a[212] = 1e+13;
    b[212] = 0;
    e[212] = 8000;

    /*reaction 214: ch2co + oh <=> hcco + h2o */
    a[213] = 1e+13;
    b[213] = 0;
    e[213] = 2000;

    /*reaction 215: ch2co + oh <=> ch2oh + co */
    a[214] = 3.73e+12;
    b[214] = 0;
    e[214] = -1013;

    /*reaction 216: ch2co (+M) <=> ch2 + co (+M) */
    a[215] = 3e+14;
    b[215] = 0;
    e[215] = 70980;

    /*reaction 217: c2h + h2 <=> c2h2 + h */
    a[216] = 409000;
    b[216] = 2.39;
    e[216] = 864.3;

    /*reaction 218: c2h + o <=> ch + co */
    a[217] = 5e+13;
    b[217] = 0;
    e[217] = 0;

    /*reaction 219: c2h + oh <=> hcco + h */
    a[218] = 2e+13;
    b[218] = 0;
    e[218] = 0;

    /*reaction 220: c2h + o2 <=> co + co + h */
    a[219] = 9.04e+12;
    b[219] = 0;
    e[219] = -457;

    /*reaction 221: hcco + c2h2 <=> h2ccch + co */
    a[220] = 1e+11;
    b[220] = 0;
    e[220] = 3000;

    /*reaction 222: hcco + h <=> ch2s + co */
    a[221] = 1e+14;
    b[221] = 0;
    e[221] = 0;

    /*reaction 223: hcco + o <=> h + co + co */
    a[222] = 8e+13;
    b[222] = 0;
    e[222] = 0;

    /*reaction 224: hcco + o <=> ch + co2 */
    a[223] = 2.95e+13;
    b[223] = 0;
    e[223] = 1113;

    /*reaction 225: hcco + o2 <=> hco + co + o */
    a[224] = 2.5e+08;
    b[224] = 1;
    e[224] = 0;

    /*reaction 226: hcco + o2 <=> co2 + hco */
    a[225] = 2.4e+11;
    b[225] = 0;
    e[225] = -854;

    /*reaction 227: hcco + ch <=> c2h2 + co */
    a[226] = 5e+13;
    b[226] = 0;
    e[226] = 0;

    /*reaction 228: hcco + hcco <=> c2h2 + co + co */
    a[227] = 1e+13;
    b[227] = 0;
    e[227] = 0;

    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double * T, double * C, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[38]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: oh + h2 <=> h + h2o */
    /*eqcon[0] *= 1;  */

    /*reaction 2: o + oh <=> o2 + h */
    /*eqcon[1] *= 1;  */

    /*reaction 3: o + h2 <=> oh + h */
    /*eqcon[2] *= 1;  */

    /*reaction 4: h + o2 (+M) <=> ho2 (+M) */
    eqcon[3] *= 1e+06; 

    /*reaction 5: h + o2 (+n2) <=> ho2 (+n2) */
    eqcon[4] *= 1e+06; 

    /*reaction 6: h + o2 (+h2) <=> ho2 (+h2) */
    eqcon[5] *= 1e+06; 

    /*reaction 7: h + o2 (+h2o) <=> ho2 (+h2o) */
    eqcon[6] *= 1e+06; 

    /*reaction 8: oh + ho2 <=> h2o + o2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: oh + ho2 <=> h2o + o2 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: h + ho2 <=> oh + oh */
    /*eqcon[9] *= 1;  */

    /*reaction 11: h + ho2 <=> h2 + o2 */
    /*eqcon[10] *= 1;  */

    /*reaction 12: h + ho2 <=> o + h2o */
    /*eqcon[11] *= 1;  */

    /*reaction 13: o + ho2 <=> o2 + oh */
    /*eqcon[12] *= 1;  */

    /*reaction 14: oh + oh <=> o + h2o */
    /*eqcon[13] *= 1;  */

    /*reaction 15: h + h + M <=> h2 + M */
    eqcon[14] *= 1e+06; 

    /*reaction 16: h + h + h2 <=> h2 + h2 */
    eqcon[15] *= 1e+06; 

    /*reaction 17: h + h + h2o <=> h2 + h2o */
    eqcon[16] *= 1e+06; 

    /*reaction 18: h + oh + M <=> h2o + M */
    eqcon[17] *= 1e+06; 

    /*reaction 19: h + o + M <=> oh + M */
    eqcon[18] *= 1e+06; 

    /*reaction 20: o + o + M <=> o2 + M */
    eqcon[19] *= 1e+06; 

    /*reaction 21: ho2 + ho2 <=> h2o2 + o2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: ho2 + ho2 <=> h2o2 + o2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: oh + oh (+M) <=> h2o2 (+M) */
    eqcon[22] *= 1e+06; 

    /*reaction 24: h2o2 + h <=> ho2 + h2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: h2o2 + h <=> oh + h2o */
    /*eqcon[24] *= 1;  */

    /*reaction 26: h2o2 + o <=> oh + ho2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: h2o2 + oh <=> h2o + ho2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: ch3 + ch3 (+M) <=> c2h6 (+M) */
    eqcon[27] *= 1e+06; 

    /*reaction 29: ch3 + h (+M) <=> ch4 (+M) */
    eqcon[28] *= 1e+06; 

    /*reaction 30: ch4 + h <=> ch3 + h2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: ch4 + oh <=> ch3 + h2o */
    /*eqcon[30] *= 1;  */

    /*reaction 32: ch4 + o <=> ch3 + oh */
    /*eqcon[31] *= 1;  */

    /*reaction 33: ch4 + ho2 <=> ch3 + h2o2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: ch3 + ho2 <=> ch3o + oh */
    /*eqcon[33] *= 1;  */

    /*reaction 35: ch3 + ho2 <=> ch4 + o2 */
    /*eqcon[34] *= 1;  */

    /*reaction 36: ch3 + o <=> ch2o + h */
    /*eqcon[35] *= 1;  */

    /*reaction 37: ch3 + o2 <=> ch3o + o */
    /*eqcon[36] *= 1;  */

    /*reaction 38: ch3 + o2 <=> ch2o + oh */
    /*eqcon[37] *= 1;  */

    /*reaction 39: ch3o + h <=> ch3 + oh */
    /*eqcon[38] *= 1;  */

    /*reaction 40: ch2oh + h <=> ch3 + oh */
    /*eqcon[39] *= 1;  */

    /*reaction 41: ch3 + oh <=> ch2s + h2o */
    /*eqcon[40] *= 1;  */

    /*reaction 42: ch3 + oh <=> ch2 + h2o */
    /*eqcon[41] *= 1;  */

    /*reaction 43: ch3 + h <=> ch2 + h2 */
    /*eqcon[42] *= 1;  */

    /*reaction 44: ch3 + M <=> ch + h2 + M */
    eqcon[43] *= 1e-06; 

    /*reaction 45: ch3 + M <=> ch2 + h + M */
    eqcon[44] *= 1e-06; 

    /*reaction 46: ch2o + h (+M) <=> ch3o (+M) */
    eqcon[45] *= 1e+06; 

    /*reaction 47: ch2o + h (+M) <=> ch2oh (+M) */
    eqcon[46] *= 1e+06; 

    /*reaction 48: ch3o + ch3 <=> ch2o + ch4 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: ch3o + h <=> ch2o + h2 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: ch2oh + h <=> ch2o + h2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: ch3o + oh <=> ch2o + h2o */
    /*eqcon[50] *= 1;  */

    /*reaction 52: ch2oh + oh <=> ch2o + h2o */
    /*eqcon[51] *= 1;  */

    /*reaction 53: ch3o + o <=> ch2o + oh */
    /*eqcon[52] *= 1;  */

    /*reaction 54: ch2oh + o <=> ch2o + oh */
    /*eqcon[53] *= 1;  */

    /*reaction 55: ch3o + o2 <=> ch2o + ho2 */
    /*eqcon[54] *= 1;  */

    /*reaction 56: ch3o + co <=> ch3 + co2 */
    /*eqcon[55] *= 1;  */

    /*reaction 57: ch2oh + o2 <=> ch2o + ho2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: ch2oh + o2 <=> ch2o + ho2 */
    /*eqcon[57] *= 1;  */

    /*reaction 59: ch2 + h <=> ch + h2 */
    /*eqcon[58] *= 1;  */

    /*reaction 60: ch2 + oh <=> ch + h2o */
    /*eqcon[59] *= 1;  */

    /*reaction 61: ch2 + oh <=> ch2o + h */
    /*eqcon[60] *= 1;  */

    /*reaction 62: ch2 + co2 <=> ch2o + co */
    /*eqcon[61] *= 1;  */

    /*reaction 63: ch2 + o <=> co + h + h */
    eqcon[62] *= 1e-06; 

    /*reaction 64: ch2 + o <=> co + h2 */
    /*eqcon[63] *= 1;  */

    /*reaction 65: ch2 + o2 <=> ch2o + o */
    /*eqcon[64] *= 1;  */

    /*reaction 66: ch2 + o2 <=> co2 + h + h */
    eqcon[65] *= 1e-06; 

    /*reaction 67: ch2 + o2 <=> co2 + h2 */
    /*eqcon[66] *= 1;  */

    /*reaction 68: ch2 + o2 <=> co + h2o */
    /*eqcon[67] *= 1;  */

    /*reaction 69: ch2 + o2 <=> hco + oh */
    /*eqcon[68] *= 1;  */

    /*reaction 70: ch2 + ch3 <=> c2h4 + h */
    /*eqcon[69] *= 1;  */

    /*reaction 71: ch2 + ch2 <=> c2h2 + h + h */
    eqcon[70] *= 1e-06; 

    /*reaction 72: ch2 + hcco <=> c2h3 + co */
    /*eqcon[71] *= 1;  */

    /*reaction 73: ch2 + c2h2 <=> h2ccch + h */
    /*eqcon[72] *= 1;  */

    /*reaction 74: ch2s + M <=> ch2 + M */
    /*eqcon[73] *= 1;  */

    /*reaction 75: ch2s + ch4 <=> ch3 + ch3 */
    /*eqcon[74] *= 1;  */

    /*reaction 76: ch2s + c2h6 <=> ch3 + c2h5 */
    /*eqcon[75] *= 1;  */

    /*reaction 77: ch2s + o2 <=> co + oh + h */
    eqcon[76] *= 1e-06; 

    /*reaction 78: ch2s + h2 <=> ch3 + h */
    /*eqcon[77] *= 1;  */

    /*reaction 79: ch2s + c2h2 <=> h2ccch + h */
    /*eqcon[78] *= 1;  */

    /*reaction 80: ch2s + o <=> co + h + h */
    eqcon[79] *= 1e-06; 

    /*reaction 81: ch2s + oh <=> ch2o + h */
    /*eqcon[80] *= 1;  */

    /*reaction 82: ch2s + h <=> ch + h2 */
    /*eqcon[81] *= 1;  */

    /*reaction 83: ch2s + co2 <=> ch2o + co */
    /*eqcon[82] *= 1;  */

    /*reaction 84: ch2s + ch3 <=> c2h4 + h */
    /*eqcon[83] *= 1;  */

    /*reaction 85: ch2s + ch2co <=> c2h4 + co */
    /*eqcon[84] *= 1;  */

    /*reaction 86: ch + o2 <=> hco + o */
    /*eqcon[85] *= 1;  */

    /*reaction 87: ch + o <=> co + h */
    /*eqcon[86] *= 1;  */

    /*reaction 88: ch + oh <=> hco + h */
    /*eqcon[87] *= 1;  */

    /*reaction 89: ch + co2 <=> hco + co */
    /*eqcon[88] *= 1;  */

    /*reaction 90: ch + h2o <=> ch2o + h */
    /*eqcon[89] *= 1;  */

    /*reaction 91: ch + ch2o <=> ch2co + h */
    /*eqcon[90] *= 1;  */

    /*reaction 92: ch + c2h2 <=> c3h2 + h */
    /*eqcon[91] *= 1;  */

    /*reaction 93: ch + ch2 <=> c2h2 + h */
    /*eqcon[92] *= 1;  */

    /*reaction 94: ch + ch3 <=> c2h3 + h */
    /*eqcon[93] *= 1;  */

    /*reaction 95: ch + ch4 <=> c2h4 + h */
    /*eqcon[94] *= 1;  */

    /*reaction 96: ch2o + oh <=> hco + h2o */
    /*eqcon[95] *= 1;  */

    /*reaction 97: ch2o + h <=> hco + h2 */
    /*eqcon[96] *= 1;  */

    /*reaction 98: ch2o + M <=> hco + h + M */
    eqcon[97] *= 1e-06; 

    /*reaction 99: ch2o + o <=> hco + oh */
    /*eqcon[98] *= 1;  */

    /*reaction 100: hco + o2 <=> co + ho2 */
    /*eqcon[99] *= 1;  */

    /*reaction 101: hco + M <=> h + co + M */
    eqcon[100] *= 1e-06; 

    /*reaction 102: hco + oh <=> h2o + co */
    /*eqcon[101] *= 1;  */

    /*reaction 103: hco + h <=> co + h2 */
    /*eqcon[102] *= 1;  */

    /*reaction 104: hco + o <=> co + oh */
    /*eqcon[103] *= 1;  */

    /*reaction 105: hco + o <=> co2 + h */
    /*eqcon[104] *= 1;  */

    /*reaction 106: co + oh <=> co2 + h */
    /*eqcon[105] *= 1;  */

    /*reaction 107: co + o + M <=> co2 + M */
    eqcon[106] *= 1e+06; 

    /*reaction 108: co + o2 <=> co2 + o */
    /*eqcon[107] *= 1;  */

    /*reaction 109: co + ho2 <=> co2 + oh */
    /*eqcon[108] *= 1;  */

    /*reaction 110: c2h5oh (+M) <=> ch3 + ch2oh (+M) */
    eqcon[109] *= 1e-06; 

    /*reaction 111: c2h5oh (+M) <=> c2h5 + oh (+M) */
    eqcon[110] *= 1e-06; 

    /*reaction 112: c2h5oh (+M) <=> c2h4 + h2o (+M) */
    eqcon[111] *= 1e-06; 

    /*reaction 113: c2h5oh (+M) <=> ch3hco + h2 (+M) */
    eqcon[112] *= 1e-06; 

    /*reaction 114: c2h5oh + oh <=> c2h4oh + h2o */
    /*eqcon[113] *= 1;  */

    /*reaction 115: c2h5oh + oh <=> ch3choh + h2o */
    /*eqcon[114] *= 1;  */

    /*reaction 116: c2h5oh + oh <=> ch3ch2o + h2o */
    /*eqcon[115] *= 1;  */

    /*reaction 117: c2h5oh + h <=> c2h4oh + h2 */
    /*eqcon[116] *= 1;  */

    /*reaction 118: c2h5oh + h <=> ch3choh + h2 */
    /*eqcon[117] *= 1;  */

    /*reaction 119: c2h5oh + h <=> ch3ch2o + h2 */
    /*eqcon[118] *= 1;  */

    /*reaction 120: c2h5oh + o <=> c2h4oh + oh */
    /*eqcon[119] *= 1;  */

    /*reaction 121: c2h5oh + o <=> ch3choh + oh */
    /*eqcon[120] *= 1;  */

    /*reaction 122: c2h5oh + o <=> ch3ch2o + oh */
    /*eqcon[121] *= 1;  */

    /*reaction 123: c2h5oh + ch3 <=> c2h4oh + ch4 */
    /*eqcon[122] *= 1;  */

    /*reaction 124: c2h5oh + ch3 <=> ch3choh + ch4 */
    /*eqcon[123] *= 1;  */

    /*reaction 125: c2h5oh + ch3 <=> ch3ch2o + ch4 */
    /*eqcon[124] *= 1;  */

    /*reaction 126: c2h5oh + ho2 <=> ch3choh + h2o2 */
    /*eqcon[125] *= 1;  */

    /*reaction 127: c2h5oh + ho2 <=> c2h4oh + h2o2 */
    /*eqcon[126] *= 1;  */

    /*reaction 128: c2h5oh + ho2 <=> ch3ch2o + h2o2 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: ch3ch2o + M <=> ch3hco + h + M */
    eqcon[128] *= 1e-06; 

    /*reaction 130: ch3ch2o + M <=> ch3 + ch2o + M */
    eqcon[129] *= 1e-06; 

    /*reaction 131: ch3ch2o + o2 <=> ch3hco + ho2 */
    /*eqcon[130] *= 1;  */

    /*reaction 132: ch3ch2o + co <=> c2h5 + co2 */
    /*eqcon[131] *= 1;  */

    /*reaction 133: ch3ch2o + h <=> ch3 + ch2oh */
    /*eqcon[132] *= 1;  */

    /*reaction 134: ch3ch2o + h <=> c2h4 + h2o */
    /*eqcon[133] *= 1;  */

    /*reaction 135: ch3ch2o + oh <=> ch3hco + h2o */
    /*eqcon[134] *= 1;  */

    /*reaction 136: ch3choh + o2 <=> ch3hco + ho2 */
    /*eqcon[135] *= 1;  */

    /*reaction 137: ch3choh + o2 <=> ch3hco + ho2 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: ch3choh + o <=> ch3hco + oh */
    /*eqcon[137] *= 1;  */

    /*reaction 139: ch3choh + h <=> c2h4 + h2o */
    /*eqcon[138] *= 1;  */

    /*reaction 140: ch3choh + h <=> ch3 + ch2oh */
    /*eqcon[139] *= 1;  */

    /*reaction 141: ch3choh + ho2 <=> ch3hco + oh + oh */
    eqcon[140] *= 1e-06; 

    /*reaction 142: ch3choh + oh <=> ch3hco + h2o */
    /*eqcon[141] *= 1;  */

    /*reaction 143: ch3choh + M <=> ch3hco + h + M */
    eqcon[142] *= 1e-06; 

    /*reaction 144: ch3hco + oh <=> ch3co + h2o */
    /*eqcon[143] *= 1;  */

    /*reaction 145: ch3hco + oh <=> ch2hco + h2o */
    /*eqcon[144] *= 1;  */

    /*reaction 146: ch3hco + o <=> ch3co + oh */
    /*eqcon[145] *= 1;  */

    /*reaction 147: ch3hco + o <=> ch2hco + oh */
    /*eqcon[146] *= 1;  */

    /*reaction 148: ch3hco + h <=> ch3co + h2 */
    /*eqcon[147] *= 1;  */

    /*reaction 149: ch3hco + h <=> ch2hco + h2 */
    /*eqcon[148] *= 1;  */

    /*reaction 150: ch3hco + ch3 <=> ch3co + ch4 */
    /*eqcon[149] *= 1;  */

    /*reaction 151: ch3hco + ch3 <=> ch2hco + ch4 */
    /*eqcon[150] *= 1;  */

    /*reaction 152: ch3hco + ho2 <=> ch3co + h2o2 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: ch3hco + ho2 <=> ch2hco + h2o2 */
    /*eqcon[152] *= 1;  */

    /*reaction 154: ch3hco + o2 <=> ch3co + ho2 */
    /*eqcon[153] *= 1;  */

    /*reaction 155: c2h6 + ch3 <=> c2h5 + ch4 */
    /*eqcon[154] *= 1;  */

    /*reaction 156: c2h6 + h <=> c2h5 + h2 */
    /*eqcon[155] *= 1;  */

    /*reaction 157: c2h6 + o <=> c2h5 + oh */
    /*eqcon[156] *= 1;  */

    /*reaction 158: c2h6 + oh <=> c2h5 + h2o */
    /*eqcon[157] *= 1;  */

    /*reaction 159: c2h5 + h <=> c2h4 + h2 */
    /*eqcon[158] *= 1;  */

    /*reaction 160: c2h5 + h <=> ch3 + ch3 */
    /*eqcon[159] *= 1;  */

    /*reaction 161: c2h5 + h <=> c2h6 */
    eqcon[160] *= 1e+06; 

    /*reaction 162: c2h5 + oh <=> c2h4 + h2o */
    /*eqcon[161] *= 1;  */

    /*reaction 163: c2h5 + o <=> ch3 + ch2o */
    /*eqcon[162] *= 1;  */

    /*reaction 164: c2h5 + ho2 <=> c2h6 + o2 */
    /*eqcon[163] *= 1;  */

    /*reaction 165: c2h5 + ho2 <=> ch3ch2o + oh */
    /*eqcon[164] *= 1;  */

    /*reaction 166: c2h5 + o2 <=> c2h4 + ho2 */
    /*eqcon[165] *= 1;  */

    /*reaction 167: c2h5 + o2 <=> ch3hco + oh */
    /*eqcon[166] *= 1;  */

    /*reaction 168: c2h4 + oh <=> c2h4oh */
    eqcon[167] *= 1e+06; 

    /*reaction 169: c2h4 + oh <=> c2h3 + h2o */
    /*eqcon[168] *= 1;  */

    /*reaction 170: c2h4 + o <=> ch3 + hco */
    /*eqcon[169] *= 1;  */

    /*reaction 171: c2h4 + o <=> ch2hco + h */
    /*eqcon[170] *= 1;  */

    /*reaction 172: c2h4 + ch3 <=> c2h3 + ch4 */
    /*eqcon[171] *= 1;  */

    /*reaction 173: c2h4 + h <=> c2h3 + h2 */
    /*eqcon[172] *= 1;  */

    /*reaction 174: c2h4 + h (+M) <=> c2h5 (+M) */
    eqcon[173] *= 1e+06; 

    /*reaction 175: c2h4 (+M) <=> c2h2 + h2 (+M) */
    eqcon[174] *= 1e-06; 

    /*reaction 176: c2h3 + h (+M) <=> c2h4 (+M) */
    eqcon[175] *= 1e+06; 

    /*reaction 177: c2h3 + h <=> c2h2 + h2 */
    /*eqcon[176] *= 1;  */

    /*reaction 178: c2h3 + o <=> ch2co + h */
    /*eqcon[177] *= 1;  */

    /*reaction 179: c2h3 + o2 <=> ch2o + hco */
    /*eqcon[178] *= 1;  */

    /*reaction 180: c2h3 + o2 <=> ch2hco + o */
    /*eqcon[179] *= 1;  */

    /*reaction 181: c2h3 + o2 <=> c2h2 + ho2 */
    /*eqcon[180] *= 1;  */

    /*reaction 182: c2h3 + oh <=> c2h2 + h2o */
    /*eqcon[181] *= 1;  */

    /*reaction 183: c2h3 + c2h <=> c2h2 + c2h2 */
    /*eqcon[182] *= 1;  */

    /*reaction 184: c2h3 + ch <=> ch2 + c2h2 */
    /*eqcon[183] *= 1;  */

    /*reaction 185: c2h3 + ch3 <=> c2h2 + ch4 */
    /*eqcon[184] *= 1;  */

    /*reaction 186: c2h2 + oh <=> c2h + h2o */
    /*eqcon[185] *= 1;  */

    /*reaction 187: c2h2 + oh <=> hccoh + h */
    /*eqcon[186] *= 1;  */

    /*reaction 188: c2h2 + oh <=> ch2co + h */
    /*eqcon[187] *= 1;  */

    /*reaction 189: c2h2 + oh <=> ch2co + h */
    /*eqcon[188] *= 1;  */

    /*reaction 190: c2h2 + oh <=> ch3 + co */
    /*eqcon[189] *= 1;  */

    /*reaction 191: hccoh + h <=> ch2co + h */
    /*eqcon[190] *= 1;  */

    /*reaction 192: c2h2 + o <=> ch2 + co */
    /*eqcon[191] *= 1;  */

    /*reaction 193: c2h2 + o <=> hcco + h */
    /*eqcon[192] *= 1;  */

    /*reaction 194: c2h2 + o <=> c2h + oh */
    /*eqcon[193] *= 1;  */

    /*reaction 195: c2h2 + ch3 <=> c2h + ch4 */
    /*eqcon[194] *= 1;  */

    /*reaction 196: c2h2 + o2 <=> hcco + oh */
    /*eqcon[195] *= 1;  */

    /*reaction 197: c2h2 + M <=> c2h + h + M */
    eqcon[196] *= 1e-06; 

    /*reaction 198: c2h2 + h (+M) <=> c2h3 (+M) */
    eqcon[197] *= 1e+06; 

    /*reaction 199: ch2hco + h <=> ch3 + hco */
    /*eqcon[198] *= 1;  */

    /*reaction 200: ch2hco + h <=> ch2co + h2 */
    /*eqcon[199] *= 1;  */

    /*reaction 201: ch2hco + o <=> ch2o + hco */
    /*eqcon[200] *= 1;  */

    /*reaction 202: ch2hco + oh <=> ch2co + h2o */
    /*eqcon[201] *= 1;  */

    /*reaction 203: ch2hco + o2 <=> ch2o + co + oh */
    eqcon[202] *= 1e-06; 

    /*reaction 204: ch2hco + ch3 <=> c2h5 + co + h */
    eqcon[203] *= 1e-06; 

    /*reaction 205: ch2hco + ho2 <=> ch2o + hco + oh */
    eqcon[204] *= 1e-06; 

    /*reaction 206: ch2hco + ho2 <=> ch3hco + o2 */
    /*eqcon[205] *= 1;  */

    /*reaction 207: ch2hco <=> ch3 + co */
    eqcon[206] *= 1e-06; 

    /*reaction 208: ch2hco <=> ch2co + h */
    eqcon[207] *= 1e-06; 

    /*reaction 209: ch3co (+M) <=> ch3 + co (+M) */
    eqcon[208] *= 1e-06; 

    /*reaction 210: ch2co + o <=> co2 + ch2 */
    /*eqcon[209] *= 1;  */

    /*reaction 211: ch2co + h <=> ch3 + co */
    /*eqcon[210] *= 1;  */

    /*reaction 212: ch2co + h <=> hcco + h2 */
    /*eqcon[211] *= 1;  */

    /*reaction 213: ch2co + o <=> hcco + oh */
    /*eqcon[212] *= 1;  */

    /*reaction 214: ch2co + oh <=> hcco + h2o */
    /*eqcon[213] *= 1;  */

    /*reaction 215: ch2co + oh <=> ch2oh + co */
    /*eqcon[214] *= 1;  */

    /*reaction 216: ch2co (+M) <=> ch2 + co (+M) */
    eqcon[215] *= 1e-06; 

    /*reaction 217: c2h + h2 <=> c2h2 + h */
    /*eqcon[216] *= 1;  */

    /*reaction 218: c2h + o <=> ch + co */
    /*eqcon[217] *= 1;  */

    /*reaction 219: c2h + oh <=> hcco + h */
    /*eqcon[218] *= 1;  */

    /*reaction 220: c2h + o2 <=> co + co + h */
    eqcon[219] *= 1e-06; 

    /*reaction 221: hcco + c2h2 <=> h2ccch + co */
    /*eqcon[220] *= 1;  */

    /*reaction 222: hcco + h <=> ch2s + co */
    /*eqcon[221] *= 1;  */

    /*reaction 223: hcco + o <=> h + co + co */
    eqcon[222] *= 1e-06; 

    /*reaction 224: hcco + o <=> ch + co2 */
    /*eqcon[223] *= 1;  */

    /*reaction 225: hcco + o2 <=> hco + co + o */
    eqcon[224] *= 1e-06; 

    /*reaction 226: hcco + o2 <=> co2 + hco */
    /*eqcon[225] *= 1;  */

    /*reaction 227: hcco + ch <=> c2h2 + co */
    /*eqcon[226] *= 1;  */

    /*reaction 228: hcco + hcco <=> c2h2 + co + co */
    eqcon[227] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[38]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: oh + h2 <=> h + h2o */
    /*eqcon[0] *= 1;  */

    /*reaction 2: o + oh <=> o2 + h */
    /*eqcon[1] *= 1;  */

    /*reaction 3: o + h2 <=> oh + h */
    /*eqcon[2] *= 1;  */

    /*reaction 4: h + o2 (+M) <=> ho2 (+M) */
    eqcon[3] *= 1e+06; 

    /*reaction 5: h + o2 (+n2) <=> ho2 (+n2) */
    eqcon[4] *= 1e+06; 

    /*reaction 6: h + o2 (+h2) <=> ho2 (+h2) */
    eqcon[5] *= 1e+06; 

    /*reaction 7: h + o2 (+h2o) <=> ho2 (+h2o) */
    eqcon[6] *= 1e+06; 

    /*reaction 8: oh + ho2 <=> h2o + o2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: oh + ho2 <=> h2o + o2 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: h + ho2 <=> oh + oh */
    /*eqcon[9] *= 1;  */

    /*reaction 11: h + ho2 <=> h2 + o2 */
    /*eqcon[10] *= 1;  */

    /*reaction 12: h + ho2 <=> o + h2o */
    /*eqcon[11] *= 1;  */

    /*reaction 13: o + ho2 <=> o2 + oh */
    /*eqcon[12] *= 1;  */

    /*reaction 14: oh + oh <=> o + h2o */
    /*eqcon[13] *= 1;  */

    /*reaction 15: h + h + M <=> h2 + M */
    eqcon[14] *= 1e+06; 

    /*reaction 16: h + h + h2 <=> h2 + h2 */
    eqcon[15] *= 1e+06; 

    /*reaction 17: h + h + h2o <=> h2 + h2o */
    eqcon[16] *= 1e+06; 

    /*reaction 18: h + oh + M <=> h2o + M */
    eqcon[17] *= 1e+06; 

    /*reaction 19: h + o + M <=> oh + M */
    eqcon[18] *= 1e+06; 

    /*reaction 20: o + o + M <=> o2 + M */
    eqcon[19] *= 1e+06; 

    /*reaction 21: ho2 + ho2 <=> h2o2 + o2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: ho2 + ho2 <=> h2o2 + o2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: oh + oh (+M) <=> h2o2 (+M) */
    eqcon[22] *= 1e+06; 

    /*reaction 24: h2o2 + h <=> ho2 + h2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: h2o2 + h <=> oh + h2o */
    /*eqcon[24] *= 1;  */

    /*reaction 26: h2o2 + o <=> oh + ho2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: h2o2 + oh <=> h2o + ho2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: ch3 + ch3 (+M) <=> c2h6 (+M) */
    eqcon[27] *= 1e+06; 

    /*reaction 29: ch3 + h (+M) <=> ch4 (+M) */
    eqcon[28] *= 1e+06; 

    /*reaction 30: ch4 + h <=> ch3 + h2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: ch4 + oh <=> ch3 + h2o */
    /*eqcon[30] *= 1;  */

    /*reaction 32: ch4 + o <=> ch3 + oh */
    /*eqcon[31] *= 1;  */

    /*reaction 33: ch4 + ho2 <=> ch3 + h2o2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: ch3 + ho2 <=> ch3o + oh */
    /*eqcon[33] *= 1;  */

    /*reaction 35: ch3 + ho2 <=> ch4 + o2 */
    /*eqcon[34] *= 1;  */

    /*reaction 36: ch3 + o <=> ch2o + h */
    /*eqcon[35] *= 1;  */

    /*reaction 37: ch3 + o2 <=> ch3o + o */
    /*eqcon[36] *= 1;  */

    /*reaction 38: ch3 + o2 <=> ch2o + oh */
    /*eqcon[37] *= 1;  */

    /*reaction 39: ch3o + h <=> ch3 + oh */
    /*eqcon[38] *= 1;  */

    /*reaction 40: ch2oh + h <=> ch3 + oh */
    /*eqcon[39] *= 1;  */

    /*reaction 41: ch3 + oh <=> ch2s + h2o */
    /*eqcon[40] *= 1;  */

    /*reaction 42: ch3 + oh <=> ch2 + h2o */
    /*eqcon[41] *= 1;  */

    /*reaction 43: ch3 + h <=> ch2 + h2 */
    /*eqcon[42] *= 1;  */

    /*reaction 44: ch3 + M <=> ch + h2 + M */
    eqcon[43] *= 1e-06; 

    /*reaction 45: ch3 + M <=> ch2 + h + M */
    eqcon[44] *= 1e-06; 

    /*reaction 46: ch2o + h (+M) <=> ch3o (+M) */
    eqcon[45] *= 1e+06; 

    /*reaction 47: ch2o + h (+M) <=> ch2oh (+M) */
    eqcon[46] *= 1e+06; 

    /*reaction 48: ch3o + ch3 <=> ch2o + ch4 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: ch3o + h <=> ch2o + h2 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: ch2oh + h <=> ch2o + h2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: ch3o + oh <=> ch2o + h2o */
    /*eqcon[50] *= 1;  */

    /*reaction 52: ch2oh + oh <=> ch2o + h2o */
    /*eqcon[51] *= 1;  */

    /*reaction 53: ch3o + o <=> ch2o + oh */
    /*eqcon[52] *= 1;  */

    /*reaction 54: ch2oh + o <=> ch2o + oh */
    /*eqcon[53] *= 1;  */

    /*reaction 55: ch3o + o2 <=> ch2o + ho2 */
    /*eqcon[54] *= 1;  */

    /*reaction 56: ch3o + co <=> ch3 + co2 */
    /*eqcon[55] *= 1;  */

    /*reaction 57: ch2oh + o2 <=> ch2o + ho2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: ch2oh + o2 <=> ch2o + ho2 */
    /*eqcon[57] *= 1;  */

    /*reaction 59: ch2 + h <=> ch + h2 */
    /*eqcon[58] *= 1;  */

    /*reaction 60: ch2 + oh <=> ch + h2o */
    /*eqcon[59] *= 1;  */

    /*reaction 61: ch2 + oh <=> ch2o + h */
    /*eqcon[60] *= 1;  */

    /*reaction 62: ch2 + co2 <=> ch2o + co */
    /*eqcon[61] *= 1;  */

    /*reaction 63: ch2 + o <=> co + h + h */
    eqcon[62] *= 1e-06; 

    /*reaction 64: ch2 + o <=> co + h2 */
    /*eqcon[63] *= 1;  */

    /*reaction 65: ch2 + o2 <=> ch2o + o */
    /*eqcon[64] *= 1;  */

    /*reaction 66: ch2 + o2 <=> co2 + h + h */
    eqcon[65] *= 1e-06; 

    /*reaction 67: ch2 + o2 <=> co2 + h2 */
    /*eqcon[66] *= 1;  */

    /*reaction 68: ch2 + o2 <=> co + h2o */
    /*eqcon[67] *= 1;  */

    /*reaction 69: ch2 + o2 <=> hco + oh */
    /*eqcon[68] *= 1;  */

    /*reaction 70: ch2 + ch3 <=> c2h4 + h */
    /*eqcon[69] *= 1;  */

    /*reaction 71: ch2 + ch2 <=> c2h2 + h + h */
    eqcon[70] *= 1e-06; 

    /*reaction 72: ch2 + hcco <=> c2h3 + co */
    /*eqcon[71] *= 1;  */

    /*reaction 73: ch2 + c2h2 <=> h2ccch + h */
    /*eqcon[72] *= 1;  */

    /*reaction 74: ch2s + M <=> ch2 + M */
    /*eqcon[73] *= 1;  */

    /*reaction 75: ch2s + ch4 <=> ch3 + ch3 */
    /*eqcon[74] *= 1;  */

    /*reaction 76: ch2s + c2h6 <=> ch3 + c2h5 */
    /*eqcon[75] *= 1;  */

    /*reaction 77: ch2s + o2 <=> co + oh + h */
    eqcon[76] *= 1e-06; 

    /*reaction 78: ch2s + h2 <=> ch3 + h */
    /*eqcon[77] *= 1;  */

    /*reaction 79: ch2s + c2h2 <=> h2ccch + h */
    /*eqcon[78] *= 1;  */

    /*reaction 80: ch2s + o <=> co + h + h */
    eqcon[79] *= 1e-06; 

    /*reaction 81: ch2s + oh <=> ch2o + h */
    /*eqcon[80] *= 1;  */

    /*reaction 82: ch2s + h <=> ch + h2 */
    /*eqcon[81] *= 1;  */

    /*reaction 83: ch2s + co2 <=> ch2o + co */
    /*eqcon[82] *= 1;  */

    /*reaction 84: ch2s + ch3 <=> c2h4 + h */
    /*eqcon[83] *= 1;  */

    /*reaction 85: ch2s + ch2co <=> c2h4 + co */
    /*eqcon[84] *= 1;  */

    /*reaction 86: ch + o2 <=> hco + o */
    /*eqcon[85] *= 1;  */

    /*reaction 87: ch + o <=> co + h */
    /*eqcon[86] *= 1;  */

    /*reaction 88: ch + oh <=> hco + h */
    /*eqcon[87] *= 1;  */

    /*reaction 89: ch + co2 <=> hco + co */
    /*eqcon[88] *= 1;  */

    /*reaction 90: ch + h2o <=> ch2o + h */
    /*eqcon[89] *= 1;  */

    /*reaction 91: ch + ch2o <=> ch2co + h */
    /*eqcon[90] *= 1;  */

    /*reaction 92: ch + c2h2 <=> c3h2 + h */
    /*eqcon[91] *= 1;  */

    /*reaction 93: ch + ch2 <=> c2h2 + h */
    /*eqcon[92] *= 1;  */

    /*reaction 94: ch + ch3 <=> c2h3 + h */
    /*eqcon[93] *= 1;  */

    /*reaction 95: ch + ch4 <=> c2h4 + h */
    /*eqcon[94] *= 1;  */

    /*reaction 96: ch2o + oh <=> hco + h2o */
    /*eqcon[95] *= 1;  */

    /*reaction 97: ch2o + h <=> hco + h2 */
    /*eqcon[96] *= 1;  */

    /*reaction 98: ch2o + M <=> hco + h + M */
    eqcon[97] *= 1e-06; 

    /*reaction 99: ch2o + o <=> hco + oh */
    /*eqcon[98] *= 1;  */

    /*reaction 100: hco + o2 <=> co + ho2 */
    /*eqcon[99] *= 1;  */

    /*reaction 101: hco + M <=> h + co + M */
    eqcon[100] *= 1e-06; 

    /*reaction 102: hco + oh <=> h2o + co */
    /*eqcon[101] *= 1;  */

    /*reaction 103: hco + h <=> co + h2 */
    /*eqcon[102] *= 1;  */

    /*reaction 104: hco + o <=> co + oh */
    /*eqcon[103] *= 1;  */

    /*reaction 105: hco + o <=> co2 + h */
    /*eqcon[104] *= 1;  */

    /*reaction 106: co + oh <=> co2 + h */
    /*eqcon[105] *= 1;  */

    /*reaction 107: co + o + M <=> co2 + M */
    eqcon[106] *= 1e+06; 

    /*reaction 108: co + o2 <=> co2 + o */
    /*eqcon[107] *= 1;  */

    /*reaction 109: co + ho2 <=> co2 + oh */
    /*eqcon[108] *= 1;  */

    /*reaction 110: c2h5oh (+M) <=> ch3 + ch2oh (+M) */
    eqcon[109] *= 1e-06; 

    /*reaction 111: c2h5oh (+M) <=> c2h5 + oh (+M) */
    eqcon[110] *= 1e-06; 

    /*reaction 112: c2h5oh (+M) <=> c2h4 + h2o (+M) */
    eqcon[111] *= 1e-06; 

    /*reaction 113: c2h5oh (+M) <=> ch3hco + h2 (+M) */
    eqcon[112] *= 1e-06; 

    /*reaction 114: c2h5oh + oh <=> c2h4oh + h2o */
    /*eqcon[113] *= 1;  */

    /*reaction 115: c2h5oh + oh <=> ch3choh + h2o */
    /*eqcon[114] *= 1;  */

    /*reaction 116: c2h5oh + oh <=> ch3ch2o + h2o */
    /*eqcon[115] *= 1;  */

    /*reaction 117: c2h5oh + h <=> c2h4oh + h2 */
    /*eqcon[116] *= 1;  */

    /*reaction 118: c2h5oh + h <=> ch3choh + h2 */
    /*eqcon[117] *= 1;  */

    /*reaction 119: c2h5oh + h <=> ch3ch2o + h2 */
    /*eqcon[118] *= 1;  */

    /*reaction 120: c2h5oh + o <=> c2h4oh + oh */
    /*eqcon[119] *= 1;  */

    /*reaction 121: c2h5oh + o <=> ch3choh + oh */
    /*eqcon[120] *= 1;  */

    /*reaction 122: c2h5oh + o <=> ch3ch2o + oh */
    /*eqcon[121] *= 1;  */

    /*reaction 123: c2h5oh + ch3 <=> c2h4oh + ch4 */
    /*eqcon[122] *= 1;  */

    /*reaction 124: c2h5oh + ch3 <=> ch3choh + ch4 */
    /*eqcon[123] *= 1;  */

    /*reaction 125: c2h5oh + ch3 <=> ch3ch2o + ch4 */
    /*eqcon[124] *= 1;  */

    /*reaction 126: c2h5oh + ho2 <=> ch3choh + h2o2 */
    /*eqcon[125] *= 1;  */

    /*reaction 127: c2h5oh + ho2 <=> c2h4oh + h2o2 */
    /*eqcon[126] *= 1;  */

    /*reaction 128: c2h5oh + ho2 <=> ch3ch2o + h2o2 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: ch3ch2o + M <=> ch3hco + h + M */
    eqcon[128] *= 1e-06; 

    /*reaction 130: ch3ch2o + M <=> ch3 + ch2o + M */
    eqcon[129] *= 1e-06; 

    /*reaction 131: ch3ch2o + o2 <=> ch3hco + ho2 */
    /*eqcon[130] *= 1;  */

    /*reaction 132: ch3ch2o + co <=> c2h5 + co2 */
    /*eqcon[131] *= 1;  */

    /*reaction 133: ch3ch2o + h <=> ch3 + ch2oh */
    /*eqcon[132] *= 1;  */

    /*reaction 134: ch3ch2o + h <=> c2h4 + h2o */
    /*eqcon[133] *= 1;  */

    /*reaction 135: ch3ch2o + oh <=> ch3hco + h2o */
    /*eqcon[134] *= 1;  */

    /*reaction 136: ch3choh + o2 <=> ch3hco + ho2 */
    /*eqcon[135] *= 1;  */

    /*reaction 137: ch3choh + o2 <=> ch3hco + ho2 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: ch3choh + o <=> ch3hco + oh */
    /*eqcon[137] *= 1;  */

    /*reaction 139: ch3choh + h <=> c2h4 + h2o */
    /*eqcon[138] *= 1;  */

    /*reaction 140: ch3choh + h <=> ch3 + ch2oh */
    /*eqcon[139] *= 1;  */

    /*reaction 141: ch3choh + ho2 <=> ch3hco + oh + oh */
    eqcon[140] *= 1e-06; 

    /*reaction 142: ch3choh + oh <=> ch3hco + h2o */
    /*eqcon[141] *= 1;  */

    /*reaction 143: ch3choh + M <=> ch3hco + h + M */
    eqcon[142] *= 1e-06; 

    /*reaction 144: ch3hco + oh <=> ch3co + h2o */
    /*eqcon[143] *= 1;  */

    /*reaction 145: ch3hco + oh <=> ch2hco + h2o */
    /*eqcon[144] *= 1;  */

    /*reaction 146: ch3hco + o <=> ch3co + oh */
    /*eqcon[145] *= 1;  */

    /*reaction 147: ch3hco + o <=> ch2hco + oh */
    /*eqcon[146] *= 1;  */

    /*reaction 148: ch3hco + h <=> ch3co + h2 */
    /*eqcon[147] *= 1;  */

    /*reaction 149: ch3hco + h <=> ch2hco + h2 */
    /*eqcon[148] *= 1;  */

    /*reaction 150: ch3hco + ch3 <=> ch3co + ch4 */
    /*eqcon[149] *= 1;  */

    /*reaction 151: ch3hco + ch3 <=> ch2hco + ch4 */
    /*eqcon[150] *= 1;  */

    /*reaction 152: ch3hco + ho2 <=> ch3co + h2o2 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: ch3hco + ho2 <=> ch2hco + h2o2 */
    /*eqcon[152] *= 1;  */

    /*reaction 154: ch3hco + o2 <=> ch3co + ho2 */
    /*eqcon[153] *= 1;  */

    /*reaction 155: c2h6 + ch3 <=> c2h5 + ch4 */
    /*eqcon[154] *= 1;  */

    /*reaction 156: c2h6 + h <=> c2h5 + h2 */
    /*eqcon[155] *= 1;  */

    /*reaction 157: c2h6 + o <=> c2h5 + oh */
    /*eqcon[156] *= 1;  */

    /*reaction 158: c2h6 + oh <=> c2h5 + h2o */
    /*eqcon[157] *= 1;  */

    /*reaction 159: c2h5 + h <=> c2h4 + h2 */
    /*eqcon[158] *= 1;  */

    /*reaction 160: c2h5 + h <=> ch3 + ch3 */
    /*eqcon[159] *= 1;  */

    /*reaction 161: c2h5 + h <=> c2h6 */
    eqcon[160] *= 1e+06; 

    /*reaction 162: c2h5 + oh <=> c2h4 + h2o */
    /*eqcon[161] *= 1;  */

    /*reaction 163: c2h5 + o <=> ch3 + ch2o */
    /*eqcon[162] *= 1;  */

    /*reaction 164: c2h5 + ho2 <=> c2h6 + o2 */
    /*eqcon[163] *= 1;  */

    /*reaction 165: c2h5 + ho2 <=> ch3ch2o + oh */
    /*eqcon[164] *= 1;  */

    /*reaction 166: c2h5 + o2 <=> c2h4 + ho2 */
    /*eqcon[165] *= 1;  */

    /*reaction 167: c2h5 + o2 <=> ch3hco + oh */
    /*eqcon[166] *= 1;  */

    /*reaction 168: c2h4 + oh <=> c2h4oh */
    eqcon[167] *= 1e+06; 

    /*reaction 169: c2h4 + oh <=> c2h3 + h2o */
    /*eqcon[168] *= 1;  */

    /*reaction 170: c2h4 + o <=> ch3 + hco */
    /*eqcon[169] *= 1;  */

    /*reaction 171: c2h4 + o <=> ch2hco + h */
    /*eqcon[170] *= 1;  */

    /*reaction 172: c2h4 + ch3 <=> c2h3 + ch4 */
    /*eqcon[171] *= 1;  */

    /*reaction 173: c2h4 + h <=> c2h3 + h2 */
    /*eqcon[172] *= 1;  */

    /*reaction 174: c2h4 + h (+M) <=> c2h5 (+M) */
    eqcon[173] *= 1e+06; 

    /*reaction 175: c2h4 (+M) <=> c2h2 + h2 (+M) */
    eqcon[174] *= 1e-06; 

    /*reaction 176: c2h3 + h (+M) <=> c2h4 (+M) */
    eqcon[175] *= 1e+06; 

    /*reaction 177: c2h3 + h <=> c2h2 + h2 */
    /*eqcon[176] *= 1;  */

    /*reaction 178: c2h3 + o <=> ch2co + h */
    /*eqcon[177] *= 1;  */

    /*reaction 179: c2h3 + o2 <=> ch2o + hco */
    /*eqcon[178] *= 1;  */

    /*reaction 180: c2h3 + o2 <=> ch2hco + o */
    /*eqcon[179] *= 1;  */

    /*reaction 181: c2h3 + o2 <=> c2h2 + ho2 */
    /*eqcon[180] *= 1;  */

    /*reaction 182: c2h3 + oh <=> c2h2 + h2o */
    /*eqcon[181] *= 1;  */

    /*reaction 183: c2h3 + c2h <=> c2h2 + c2h2 */
    /*eqcon[182] *= 1;  */

    /*reaction 184: c2h3 + ch <=> ch2 + c2h2 */
    /*eqcon[183] *= 1;  */

    /*reaction 185: c2h3 + ch3 <=> c2h2 + ch4 */
    /*eqcon[184] *= 1;  */

    /*reaction 186: c2h2 + oh <=> c2h + h2o */
    /*eqcon[185] *= 1;  */

    /*reaction 187: c2h2 + oh <=> hccoh + h */
    /*eqcon[186] *= 1;  */

    /*reaction 188: c2h2 + oh <=> ch2co + h */
    /*eqcon[187] *= 1;  */

    /*reaction 189: c2h2 + oh <=> ch2co + h */
    /*eqcon[188] *= 1;  */

    /*reaction 190: c2h2 + oh <=> ch3 + co */
    /*eqcon[189] *= 1;  */

    /*reaction 191: hccoh + h <=> ch2co + h */
    /*eqcon[190] *= 1;  */

    /*reaction 192: c2h2 + o <=> ch2 + co */
    /*eqcon[191] *= 1;  */

    /*reaction 193: c2h2 + o <=> hcco + h */
    /*eqcon[192] *= 1;  */

    /*reaction 194: c2h2 + o <=> c2h + oh */
    /*eqcon[193] *= 1;  */

    /*reaction 195: c2h2 + ch3 <=> c2h + ch4 */
    /*eqcon[194] *= 1;  */

    /*reaction 196: c2h2 + o2 <=> hcco + oh */
    /*eqcon[195] *= 1;  */

    /*reaction 197: c2h2 + M <=> c2h + h + M */
    eqcon[196] *= 1e-06; 

    /*reaction 198: c2h2 + h (+M) <=> c2h3 (+M) */
    eqcon[197] *= 1e+06; 

    /*reaction 199: ch2hco + h <=> ch3 + hco */
    /*eqcon[198] *= 1;  */

    /*reaction 200: ch2hco + h <=> ch2co + h2 */
    /*eqcon[199] *= 1;  */

    /*reaction 201: ch2hco + o <=> ch2o + hco */
    /*eqcon[200] *= 1;  */

    /*reaction 202: ch2hco + oh <=> ch2co + h2o */
    /*eqcon[201] *= 1;  */

    /*reaction 203: ch2hco + o2 <=> ch2o + co + oh */
    eqcon[202] *= 1e-06; 

    /*reaction 204: ch2hco + ch3 <=> c2h5 + co + h */
    eqcon[203] *= 1e-06; 

    /*reaction 205: ch2hco + ho2 <=> ch2o + hco + oh */
    eqcon[204] *= 1e-06; 

    /*reaction 206: ch2hco + ho2 <=> ch3hco + o2 */
    /*eqcon[205] *= 1;  */

    /*reaction 207: ch2hco <=> ch3 + co */
    eqcon[206] *= 1e-06; 

    /*reaction 208: ch2hco <=> ch2co + h */
    eqcon[207] *= 1e-06; 

    /*reaction 209: ch3co (+M) <=> ch3 + co (+M) */
    eqcon[208] *= 1e-06; 

    /*reaction 210: ch2co + o <=> co2 + ch2 */
    /*eqcon[209] *= 1;  */

    /*reaction 211: ch2co + h <=> ch3 + co */
    /*eqcon[210] *= 1;  */

    /*reaction 212: ch2co + h <=> hcco + h2 */
    /*eqcon[211] *= 1;  */

    /*reaction 213: ch2co + o <=> hcco + oh */
    /*eqcon[212] *= 1;  */

    /*reaction 214: ch2co + oh <=> hcco + h2o */
    /*eqcon[213] *= 1;  */

    /*reaction 215: ch2co + oh <=> ch2oh + co */
    /*eqcon[214] *= 1;  */

    /*reaction 216: ch2co (+M) <=> ch2 + co (+M) */
    eqcon[215] *= 1e-06; 

    /*reaction 217: c2h + h2 <=> c2h2 + h */
    /*eqcon[216] *= 1;  */

    /*reaction 218: c2h + o <=> ch + co */
    /*eqcon[217] *= 1;  */

    /*reaction 219: c2h + oh <=> hcco + h */
    /*eqcon[218] *= 1;  */

    /*reaction 220: c2h + o2 <=> co + co + h */
    eqcon[219] *= 1e-06; 

    /*reaction 221: hcco + c2h2 <=> h2ccch + co */
    /*eqcon[220] *= 1;  */

    /*reaction 222: hcco + h <=> ch2s + co */
    /*eqcon[221] *= 1;  */

    /*reaction 223: hcco + o <=> h + co + co */
    eqcon[222] *= 1e-06; 

    /*reaction 224: hcco + o <=> ch + co2 */
    /*eqcon[223] *= 1;  */

    /*reaction 225: hcco + o2 <=> hco + co + o */
    eqcon[224] *= 1e-06; 

    /*reaction 226: hcco + o2 <=> co2 + hco */
    /*eqcon[225] *= 1;  */

    /*reaction 227: hcco + ch <=> c2h2 + co */
    /*eqcon[226] *= 1;  */

    /*reaction 228: hcco + hcco <=> c2h2 + co + co */
    eqcon[227] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[38]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: oh + h2 <=> h + h2o */
    /*eqcon[0] *= 1;  */

    /*reaction 2: o + oh <=> o2 + h */
    /*eqcon[1] *= 1;  */

    /*reaction 3: o + h2 <=> oh + h */
    /*eqcon[2] *= 1;  */

    /*reaction 4: h + o2 (+M) <=> ho2 (+M) */
    eqcon[3] *= 1e+06; 

    /*reaction 5: h + o2 (+n2) <=> ho2 (+n2) */
    eqcon[4] *= 1e+06; 

    /*reaction 6: h + o2 (+h2) <=> ho2 (+h2) */
    eqcon[5] *= 1e+06; 

    /*reaction 7: h + o2 (+h2o) <=> ho2 (+h2o) */
    eqcon[6] *= 1e+06; 

    /*reaction 8: oh + ho2 <=> h2o + o2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: oh + ho2 <=> h2o + o2 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: h + ho2 <=> oh + oh */
    /*eqcon[9] *= 1;  */

    /*reaction 11: h + ho2 <=> h2 + o2 */
    /*eqcon[10] *= 1;  */

    /*reaction 12: h + ho2 <=> o + h2o */
    /*eqcon[11] *= 1;  */

    /*reaction 13: o + ho2 <=> o2 + oh */
    /*eqcon[12] *= 1;  */

    /*reaction 14: oh + oh <=> o + h2o */
    /*eqcon[13] *= 1;  */

    /*reaction 15: h + h + M <=> h2 + M */
    eqcon[14] *= 1e+06; 

    /*reaction 16: h + h + h2 <=> h2 + h2 */
    eqcon[15] *= 1e+06; 

    /*reaction 17: h + h + h2o <=> h2 + h2o */
    eqcon[16] *= 1e+06; 

    /*reaction 18: h + oh + M <=> h2o + M */
    eqcon[17] *= 1e+06; 

    /*reaction 19: h + o + M <=> oh + M */
    eqcon[18] *= 1e+06; 

    /*reaction 20: o + o + M <=> o2 + M */
    eqcon[19] *= 1e+06; 

    /*reaction 21: ho2 + ho2 <=> h2o2 + o2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: ho2 + ho2 <=> h2o2 + o2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: oh + oh (+M) <=> h2o2 (+M) */
    eqcon[22] *= 1e+06; 

    /*reaction 24: h2o2 + h <=> ho2 + h2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: h2o2 + h <=> oh + h2o */
    /*eqcon[24] *= 1;  */

    /*reaction 26: h2o2 + o <=> oh + ho2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: h2o2 + oh <=> h2o + ho2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: ch3 + ch3 (+M) <=> c2h6 (+M) */
    eqcon[27] *= 1e+06; 

    /*reaction 29: ch3 + h (+M) <=> ch4 (+M) */
    eqcon[28] *= 1e+06; 

    /*reaction 30: ch4 + h <=> ch3 + h2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: ch4 + oh <=> ch3 + h2o */
    /*eqcon[30] *= 1;  */

    /*reaction 32: ch4 + o <=> ch3 + oh */
    /*eqcon[31] *= 1;  */

    /*reaction 33: ch4 + ho2 <=> ch3 + h2o2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: ch3 + ho2 <=> ch3o + oh */
    /*eqcon[33] *= 1;  */

    /*reaction 35: ch3 + ho2 <=> ch4 + o2 */
    /*eqcon[34] *= 1;  */

    /*reaction 36: ch3 + o <=> ch2o + h */
    /*eqcon[35] *= 1;  */

    /*reaction 37: ch3 + o2 <=> ch3o + o */
    /*eqcon[36] *= 1;  */

    /*reaction 38: ch3 + o2 <=> ch2o + oh */
    /*eqcon[37] *= 1;  */

    /*reaction 39: ch3o + h <=> ch3 + oh */
    /*eqcon[38] *= 1;  */

    /*reaction 40: ch2oh + h <=> ch3 + oh */
    /*eqcon[39] *= 1;  */

    /*reaction 41: ch3 + oh <=> ch2s + h2o */
    /*eqcon[40] *= 1;  */

    /*reaction 42: ch3 + oh <=> ch2 + h2o */
    /*eqcon[41] *= 1;  */

    /*reaction 43: ch3 + h <=> ch2 + h2 */
    /*eqcon[42] *= 1;  */

    /*reaction 44: ch3 + M <=> ch + h2 + M */
    eqcon[43] *= 1e-06; 

    /*reaction 45: ch3 + M <=> ch2 + h + M */
    eqcon[44] *= 1e-06; 

    /*reaction 46: ch2o + h (+M) <=> ch3o (+M) */
    eqcon[45] *= 1e+06; 

    /*reaction 47: ch2o + h (+M) <=> ch2oh (+M) */
    eqcon[46] *= 1e+06; 

    /*reaction 48: ch3o + ch3 <=> ch2o + ch4 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: ch3o + h <=> ch2o + h2 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: ch2oh + h <=> ch2o + h2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: ch3o + oh <=> ch2o + h2o */
    /*eqcon[50] *= 1;  */

    /*reaction 52: ch2oh + oh <=> ch2o + h2o */
    /*eqcon[51] *= 1;  */

    /*reaction 53: ch3o + o <=> ch2o + oh */
    /*eqcon[52] *= 1;  */

    /*reaction 54: ch2oh + o <=> ch2o + oh */
    /*eqcon[53] *= 1;  */

    /*reaction 55: ch3o + o2 <=> ch2o + ho2 */
    /*eqcon[54] *= 1;  */

    /*reaction 56: ch3o + co <=> ch3 + co2 */
    /*eqcon[55] *= 1;  */

    /*reaction 57: ch2oh + o2 <=> ch2o + ho2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: ch2oh + o2 <=> ch2o + ho2 */
    /*eqcon[57] *= 1;  */

    /*reaction 59: ch2 + h <=> ch + h2 */
    /*eqcon[58] *= 1;  */

    /*reaction 60: ch2 + oh <=> ch + h2o */
    /*eqcon[59] *= 1;  */

    /*reaction 61: ch2 + oh <=> ch2o + h */
    /*eqcon[60] *= 1;  */

    /*reaction 62: ch2 + co2 <=> ch2o + co */
    /*eqcon[61] *= 1;  */

    /*reaction 63: ch2 + o <=> co + h + h */
    eqcon[62] *= 1e-06; 

    /*reaction 64: ch2 + o <=> co + h2 */
    /*eqcon[63] *= 1;  */

    /*reaction 65: ch2 + o2 <=> ch2o + o */
    /*eqcon[64] *= 1;  */

    /*reaction 66: ch2 + o2 <=> co2 + h + h */
    eqcon[65] *= 1e-06; 

    /*reaction 67: ch2 + o2 <=> co2 + h2 */
    /*eqcon[66] *= 1;  */

    /*reaction 68: ch2 + o2 <=> co + h2o */
    /*eqcon[67] *= 1;  */

    /*reaction 69: ch2 + o2 <=> hco + oh */
    /*eqcon[68] *= 1;  */

    /*reaction 70: ch2 + ch3 <=> c2h4 + h */
    /*eqcon[69] *= 1;  */

    /*reaction 71: ch2 + ch2 <=> c2h2 + h + h */
    eqcon[70] *= 1e-06; 

    /*reaction 72: ch2 + hcco <=> c2h3 + co */
    /*eqcon[71] *= 1;  */

    /*reaction 73: ch2 + c2h2 <=> h2ccch + h */
    /*eqcon[72] *= 1;  */

    /*reaction 74: ch2s + M <=> ch2 + M */
    /*eqcon[73] *= 1;  */

    /*reaction 75: ch2s + ch4 <=> ch3 + ch3 */
    /*eqcon[74] *= 1;  */

    /*reaction 76: ch2s + c2h6 <=> ch3 + c2h5 */
    /*eqcon[75] *= 1;  */

    /*reaction 77: ch2s + o2 <=> co + oh + h */
    eqcon[76] *= 1e-06; 

    /*reaction 78: ch2s + h2 <=> ch3 + h */
    /*eqcon[77] *= 1;  */

    /*reaction 79: ch2s + c2h2 <=> h2ccch + h */
    /*eqcon[78] *= 1;  */

    /*reaction 80: ch2s + o <=> co + h + h */
    eqcon[79] *= 1e-06; 

    /*reaction 81: ch2s + oh <=> ch2o + h */
    /*eqcon[80] *= 1;  */

    /*reaction 82: ch2s + h <=> ch + h2 */
    /*eqcon[81] *= 1;  */

    /*reaction 83: ch2s + co2 <=> ch2o + co */
    /*eqcon[82] *= 1;  */

    /*reaction 84: ch2s + ch3 <=> c2h4 + h */
    /*eqcon[83] *= 1;  */

    /*reaction 85: ch2s + ch2co <=> c2h4 + co */
    /*eqcon[84] *= 1;  */

    /*reaction 86: ch + o2 <=> hco + o */
    /*eqcon[85] *= 1;  */

    /*reaction 87: ch + o <=> co + h */
    /*eqcon[86] *= 1;  */

    /*reaction 88: ch + oh <=> hco + h */
    /*eqcon[87] *= 1;  */

    /*reaction 89: ch + co2 <=> hco + co */
    /*eqcon[88] *= 1;  */

    /*reaction 90: ch + h2o <=> ch2o + h */
    /*eqcon[89] *= 1;  */

    /*reaction 91: ch + ch2o <=> ch2co + h */
    /*eqcon[90] *= 1;  */

    /*reaction 92: ch + c2h2 <=> c3h2 + h */
    /*eqcon[91] *= 1;  */

    /*reaction 93: ch + ch2 <=> c2h2 + h */
    /*eqcon[92] *= 1;  */

    /*reaction 94: ch + ch3 <=> c2h3 + h */
    /*eqcon[93] *= 1;  */

    /*reaction 95: ch + ch4 <=> c2h4 + h */
    /*eqcon[94] *= 1;  */

    /*reaction 96: ch2o + oh <=> hco + h2o */
    /*eqcon[95] *= 1;  */

    /*reaction 97: ch2o + h <=> hco + h2 */
    /*eqcon[96] *= 1;  */

    /*reaction 98: ch2o + M <=> hco + h + M */
    eqcon[97] *= 1e-06; 

    /*reaction 99: ch2o + o <=> hco + oh */
    /*eqcon[98] *= 1;  */

    /*reaction 100: hco + o2 <=> co + ho2 */
    /*eqcon[99] *= 1;  */

    /*reaction 101: hco + M <=> h + co + M */
    eqcon[100] *= 1e-06; 

    /*reaction 102: hco + oh <=> h2o + co */
    /*eqcon[101] *= 1;  */

    /*reaction 103: hco + h <=> co + h2 */
    /*eqcon[102] *= 1;  */

    /*reaction 104: hco + o <=> co + oh */
    /*eqcon[103] *= 1;  */

    /*reaction 105: hco + o <=> co2 + h */
    /*eqcon[104] *= 1;  */

    /*reaction 106: co + oh <=> co2 + h */
    /*eqcon[105] *= 1;  */

    /*reaction 107: co + o + M <=> co2 + M */
    eqcon[106] *= 1e+06; 

    /*reaction 108: co + o2 <=> co2 + o */
    /*eqcon[107] *= 1;  */

    /*reaction 109: co + ho2 <=> co2 + oh */
    /*eqcon[108] *= 1;  */

    /*reaction 110: c2h5oh (+M) <=> ch3 + ch2oh (+M) */
    eqcon[109] *= 1e-06; 

    /*reaction 111: c2h5oh (+M) <=> c2h5 + oh (+M) */
    eqcon[110] *= 1e-06; 

    /*reaction 112: c2h5oh (+M) <=> c2h4 + h2o (+M) */
    eqcon[111] *= 1e-06; 

    /*reaction 113: c2h5oh (+M) <=> ch3hco + h2 (+M) */
    eqcon[112] *= 1e-06; 

    /*reaction 114: c2h5oh + oh <=> c2h4oh + h2o */
    /*eqcon[113] *= 1;  */

    /*reaction 115: c2h5oh + oh <=> ch3choh + h2o */
    /*eqcon[114] *= 1;  */

    /*reaction 116: c2h5oh + oh <=> ch3ch2o + h2o */
    /*eqcon[115] *= 1;  */

    /*reaction 117: c2h5oh + h <=> c2h4oh + h2 */
    /*eqcon[116] *= 1;  */

    /*reaction 118: c2h5oh + h <=> ch3choh + h2 */
    /*eqcon[117] *= 1;  */

    /*reaction 119: c2h5oh + h <=> ch3ch2o + h2 */
    /*eqcon[118] *= 1;  */

    /*reaction 120: c2h5oh + o <=> c2h4oh + oh */
    /*eqcon[119] *= 1;  */

    /*reaction 121: c2h5oh + o <=> ch3choh + oh */
    /*eqcon[120] *= 1;  */

    /*reaction 122: c2h5oh + o <=> ch3ch2o + oh */
    /*eqcon[121] *= 1;  */

    /*reaction 123: c2h5oh + ch3 <=> c2h4oh + ch4 */
    /*eqcon[122] *= 1;  */

    /*reaction 124: c2h5oh + ch3 <=> ch3choh + ch4 */
    /*eqcon[123] *= 1;  */

    /*reaction 125: c2h5oh + ch3 <=> ch3ch2o + ch4 */
    /*eqcon[124] *= 1;  */

    /*reaction 126: c2h5oh + ho2 <=> ch3choh + h2o2 */
    /*eqcon[125] *= 1;  */

    /*reaction 127: c2h5oh + ho2 <=> c2h4oh + h2o2 */
    /*eqcon[126] *= 1;  */

    /*reaction 128: c2h5oh + ho2 <=> ch3ch2o + h2o2 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: ch3ch2o + M <=> ch3hco + h + M */
    eqcon[128] *= 1e-06; 

    /*reaction 130: ch3ch2o + M <=> ch3 + ch2o + M */
    eqcon[129] *= 1e-06; 

    /*reaction 131: ch3ch2o + o2 <=> ch3hco + ho2 */
    /*eqcon[130] *= 1;  */

    /*reaction 132: ch3ch2o + co <=> c2h5 + co2 */
    /*eqcon[131] *= 1;  */

    /*reaction 133: ch3ch2o + h <=> ch3 + ch2oh */
    /*eqcon[132] *= 1;  */

    /*reaction 134: ch3ch2o + h <=> c2h4 + h2o */
    /*eqcon[133] *= 1;  */

    /*reaction 135: ch3ch2o + oh <=> ch3hco + h2o */
    /*eqcon[134] *= 1;  */

    /*reaction 136: ch3choh + o2 <=> ch3hco + ho2 */
    /*eqcon[135] *= 1;  */

    /*reaction 137: ch3choh + o2 <=> ch3hco + ho2 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: ch3choh + o <=> ch3hco + oh */
    /*eqcon[137] *= 1;  */

    /*reaction 139: ch3choh + h <=> c2h4 + h2o */
    /*eqcon[138] *= 1;  */

    /*reaction 140: ch3choh + h <=> ch3 + ch2oh */
    /*eqcon[139] *= 1;  */

    /*reaction 141: ch3choh + ho2 <=> ch3hco + oh + oh */
    eqcon[140] *= 1e-06; 

    /*reaction 142: ch3choh + oh <=> ch3hco + h2o */
    /*eqcon[141] *= 1;  */

    /*reaction 143: ch3choh + M <=> ch3hco + h + M */
    eqcon[142] *= 1e-06; 

    /*reaction 144: ch3hco + oh <=> ch3co + h2o */
    /*eqcon[143] *= 1;  */

    /*reaction 145: ch3hco + oh <=> ch2hco + h2o */
    /*eqcon[144] *= 1;  */

    /*reaction 146: ch3hco + o <=> ch3co + oh */
    /*eqcon[145] *= 1;  */

    /*reaction 147: ch3hco + o <=> ch2hco + oh */
    /*eqcon[146] *= 1;  */

    /*reaction 148: ch3hco + h <=> ch3co + h2 */
    /*eqcon[147] *= 1;  */

    /*reaction 149: ch3hco + h <=> ch2hco + h2 */
    /*eqcon[148] *= 1;  */

    /*reaction 150: ch3hco + ch3 <=> ch3co + ch4 */
    /*eqcon[149] *= 1;  */

    /*reaction 151: ch3hco + ch3 <=> ch2hco + ch4 */
    /*eqcon[150] *= 1;  */

    /*reaction 152: ch3hco + ho2 <=> ch3co + h2o2 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: ch3hco + ho2 <=> ch2hco + h2o2 */
    /*eqcon[152] *= 1;  */

    /*reaction 154: ch3hco + o2 <=> ch3co + ho2 */
    /*eqcon[153] *= 1;  */

    /*reaction 155: c2h6 + ch3 <=> c2h5 + ch4 */
    /*eqcon[154] *= 1;  */

    /*reaction 156: c2h6 + h <=> c2h5 + h2 */
    /*eqcon[155] *= 1;  */

    /*reaction 157: c2h6 + o <=> c2h5 + oh */
    /*eqcon[156] *= 1;  */

    /*reaction 158: c2h6 + oh <=> c2h5 + h2o */
    /*eqcon[157] *= 1;  */

    /*reaction 159: c2h5 + h <=> c2h4 + h2 */
    /*eqcon[158] *= 1;  */

    /*reaction 160: c2h5 + h <=> ch3 + ch3 */
    /*eqcon[159] *= 1;  */

    /*reaction 161: c2h5 + h <=> c2h6 */
    eqcon[160] *= 1e+06; 

    /*reaction 162: c2h5 + oh <=> c2h4 + h2o */
    /*eqcon[161] *= 1;  */

    /*reaction 163: c2h5 + o <=> ch3 + ch2o */
    /*eqcon[162] *= 1;  */

    /*reaction 164: c2h5 + ho2 <=> c2h6 + o2 */
    /*eqcon[163] *= 1;  */

    /*reaction 165: c2h5 + ho2 <=> ch3ch2o + oh */
    /*eqcon[164] *= 1;  */

    /*reaction 166: c2h5 + o2 <=> c2h4 + ho2 */
    /*eqcon[165] *= 1;  */

    /*reaction 167: c2h5 + o2 <=> ch3hco + oh */
    /*eqcon[166] *= 1;  */

    /*reaction 168: c2h4 + oh <=> c2h4oh */
    eqcon[167] *= 1e+06; 

    /*reaction 169: c2h4 + oh <=> c2h3 + h2o */
    /*eqcon[168] *= 1;  */

    /*reaction 170: c2h4 + o <=> ch3 + hco */
    /*eqcon[169] *= 1;  */

    /*reaction 171: c2h4 + o <=> ch2hco + h */
    /*eqcon[170] *= 1;  */

    /*reaction 172: c2h4 + ch3 <=> c2h3 + ch4 */
    /*eqcon[171] *= 1;  */

    /*reaction 173: c2h4 + h <=> c2h3 + h2 */
    /*eqcon[172] *= 1;  */

    /*reaction 174: c2h4 + h (+M) <=> c2h5 (+M) */
    eqcon[173] *= 1e+06; 

    /*reaction 175: c2h4 (+M) <=> c2h2 + h2 (+M) */
    eqcon[174] *= 1e-06; 

    /*reaction 176: c2h3 + h (+M) <=> c2h4 (+M) */
    eqcon[175] *= 1e+06; 

    /*reaction 177: c2h3 + h <=> c2h2 + h2 */
    /*eqcon[176] *= 1;  */

    /*reaction 178: c2h3 + o <=> ch2co + h */
    /*eqcon[177] *= 1;  */

    /*reaction 179: c2h3 + o2 <=> ch2o + hco */
    /*eqcon[178] *= 1;  */

    /*reaction 180: c2h3 + o2 <=> ch2hco + o */
    /*eqcon[179] *= 1;  */

    /*reaction 181: c2h3 + o2 <=> c2h2 + ho2 */
    /*eqcon[180] *= 1;  */

    /*reaction 182: c2h3 + oh <=> c2h2 + h2o */
    /*eqcon[181] *= 1;  */

    /*reaction 183: c2h3 + c2h <=> c2h2 + c2h2 */
    /*eqcon[182] *= 1;  */

    /*reaction 184: c2h3 + ch <=> ch2 + c2h2 */
    /*eqcon[183] *= 1;  */

    /*reaction 185: c2h3 + ch3 <=> c2h2 + ch4 */
    /*eqcon[184] *= 1;  */

    /*reaction 186: c2h2 + oh <=> c2h + h2o */
    /*eqcon[185] *= 1;  */

    /*reaction 187: c2h2 + oh <=> hccoh + h */
    /*eqcon[186] *= 1;  */

    /*reaction 188: c2h2 + oh <=> ch2co + h */
    /*eqcon[187] *= 1;  */

    /*reaction 189: c2h2 + oh <=> ch2co + h */
    /*eqcon[188] *= 1;  */

    /*reaction 190: c2h2 + oh <=> ch3 + co */
    /*eqcon[189] *= 1;  */

    /*reaction 191: hccoh + h <=> ch2co + h */
    /*eqcon[190] *= 1;  */

    /*reaction 192: c2h2 + o <=> ch2 + co */
    /*eqcon[191] *= 1;  */

    /*reaction 193: c2h2 + o <=> hcco + h */
    /*eqcon[192] *= 1;  */

    /*reaction 194: c2h2 + o <=> c2h + oh */
    /*eqcon[193] *= 1;  */

    /*reaction 195: c2h2 + ch3 <=> c2h + ch4 */
    /*eqcon[194] *= 1;  */

    /*reaction 196: c2h2 + o2 <=> hcco + oh */
    /*eqcon[195] *= 1;  */

    /*reaction 197: c2h2 + M <=> c2h + h + M */
    eqcon[196] *= 1e-06; 

    /*reaction 198: c2h2 + h (+M) <=> c2h3 (+M) */
    eqcon[197] *= 1e+06; 

    /*reaction 199: ch2hco + h <=> ch3 + hco */
    /*eqcon[198] *= 1;  */

    /*reaction 200: ch2hco + h <=> ch2co + h2 */
    /*eqcon[199] *= 1;  */

    /*reaction 201: ch2hco + o <=> ch2o + hco */
    /*eqcon[200] *= 1;  */

    /*reaction 202: ch2hco + oh <=> ch2co + h2o */
    /*eqcon[201] *= 1;  */

    /*reaction 203: ch2hco + o2 <=> ch2o + co + oh */
    eqcon[202] *= 1e-06; 

    /*reaction 204: ch2hco + ch3 <=> c2h5 + co + h */
    eqcon[203] *= 1e-06; 

    /*reaction 205: ch2hco + ho2 <=> ch2o + hco + oh */
    eqcon[204] *= 1e-06; 

    /*reaction 206: ch2hco + ho2 <=> ch3hco + o2 */
    /*eqcon[205] *= 1;  */

    /*reaction 207: ch2hco <=> ch3 + co */
    eqcon[206] *= 1e-06; 

    /*reaction 208: ch2hco <=> ch2co + h */
    eqcon[207] *= 1e-06; 

    /*reaction 209: ch3co (+M) <=> ch3 + co (+M) */
    eqcon[208] *= 1e-06; 

    /*reaction 210: ch2co + o <=> co2 + ch2 */
    /*eqcon[209] *= 1;  */

    /*reaction 211: ch2co + h <=> ch3 + co */
    /*eqcon[210] *= 1;  */

    /*reaction 212: ch2co + h <=> hcco + h2 */
    /*eqcon[211] *= 1;  */

    /*reaction 213: ch2co + o <=> hcco + oh */
    /*eqcon[212] *= 1;  */

    /*reaction 214: ch2co + oh <=> hcco + h2o */
    /*eqcon[213] *= 1;  */

    /*reaction 215: ch2co + oh <=> ch2oh + co */
    /*eqcon[214] *= 1;  */

    /*reaction 216: ch2co (+M) <=> ch2 + co (+M) */
    eqcon[215] *= 1e-06; 

    /*reaction 217: c2h + h2 <=> c2h2 + h */
    /*eqcon[216] *= 1;  */

    /*reaction 218: c2h + o <=> ch + co */
    /*eqcon[217] *= 1;  */

    /*reaction 219: c2h + oh <=> hcco + h */
    /*eqcon[218] *= 1;  */

    /*reaction 220: c2h + o2 <=> co + co + h */
    eqcon[219] *= 1e-06; 

    /*reaction 221: hcco + c2h2 <=> h2ccch + co */
    /*eqcon[220] *= 1;  */

    /*reaction 222: hcco + h <=> ch2s + co */
    /*eqcon[221] *= 1;  */

    /*reaction 223: hcco + o <=> h + co + co */
    eqcon[222] *= 1e-06; 

    /*reaction 224: hcco + o <=> ch + co2 */
    /*eqcon[223] *= 1;  */

    /*reaction 225: hcco + o2 <=> hco + co + o */
    eqcon[224] *= 1e-06; 

    /*reaction 226: hcco + o2 <=> co2 + hco */
    /*eqcon[225] *= 1;  */

    /*reaction 227: hcco + ch <=> c2h2 + co */
    /*eqcon[226] *= 1;  */

    /*reaction 228: hcco + hcco <=> c2h2 + co + co */
    eqcon[227] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[38]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: oh + h2 <=> h + h2o */
    /*eqcon[0] *= 1;  */

    /*reaction 2: o + oh <=> o2 + h */
    /*eqcon[1] *= 1;  */

    /*reaction 3: o + h2 <=> oh + h */
    /*eqcon[2] *= 1;  */

    /*reaction 4: h + o2 (+M) <=> ho2 (+M) */
    eqcon[3] *= 1e+06; 

    /*reaction 5: h + o2 (+n2) <=> ho2 (+n2) */
    eqcon[4] *= 1e+06; 

    /*reaction 6: h + o2 (+h2) <=> ho2 (+h2) */
    eqcon[5] *= 1e+06; 

    /*reaction 7: h + o2 (+h2o) <=> ho2 (+h2o) */
    eqcon[6] *= 1e+06; 

    /*reaction 8: oh + ho2 <=> h2o + o2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: oh + ho2 <=> h2o + o2 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: h + ho2 <=> oh + oh */
    /*eqcon[9] *= 1;  */

    /*reaction 11: h + ho2 <=> h2 + o2 */
    /*eqcon[10] *= 1;  */

    /*reaction 12: h + ho2 <=> o + h2o */
    /*eqcon[11] *= 1;  */

    /*reaction 13: o + ho2 <=> o2 + oh */
    /*eqcon[12] *= 1;  */

    /*reaction 14: oh + oh <=> o + h2o */
    /*eqcon[13] *= 1;  */

    /*reaction 15: h + h + M <=> h2 + M */
    eqcon[14] *= 1e+06; 

    /*reaction 16: h + h + h2 <=> h2 + h2 */
    eqcon[15] *= 1e+06; 

    /*reaction 17: h + h + h2o <=> h2 + h2o */
    eqcon[16] *= 1e+06; 

    /*reaction 18: h + oh + M <=> h2o + M */
    eqcon[17] *= 1e+06; 

    /*reaction 19: h + o + M <=> oh + M */
    eqcon[18] *= 1e+06; 

    /*reaction 20: o + o + M <=> o2 + M */
    eqcon[19] *= 1e+06; 

    /*reaction 21: ho2 + ho2 <=> h2o2 + o2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: ho2 + ho2 <=> h2o2 + o2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: oh + oh (+M) <=> h2o2 (+M) */
    eqcon[22] *= 1e+06; 

    /*reaction 24: h2o2 + h <=> ho2 + h2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: h2o2 + h <=> oh + h2o */
    /*eqcon[24] *= 1;  */

    /*reaction 26: h2o2 + o <=> oh + ho2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: h2o2 + oh <=> h2o + ho2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: ch3 + ch3 (+M) <=> c2h6 (+M) */
    eqcon[27] *= 1e+06; 

    /*reaction 29: ch3 + h (+M) <=> ch4 (+M) */
    eqcon[28] *= 1e+06; 

    /*reaction 30: ch4 + h <=> ch3 + h2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: ch4 + oh <=> ch3 + h2o */
    /*eqcon[30] *= 1;  */

    /*reaction 32: ch4 + o <=> ch3 + oh */
    /*eqcon[31] *= 1;  */

    /*reaction 33: ch4 + ho2 <=> ch3 + h2o2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: ch3 + ho2 <=> ch3o + oh */
    /*eqcon[33] *= 1;  */

    /*reaction 35: ch3 + ho2 <=> ch4 + o2 */
    /*eqcon[34] *= 1;  */

    /*reaction 36: ch3 + o <=> ch2o + h */
    /*eqcon[35] *= 1;  */

    /*reaction 37: ch3 + o2 <=> ch3o + o */
    /*eqcon[36] *= 1;  */

    /*reaction 38: ch3 + o2 <=> ch2o + oh */
    /*eqcon[37] *= 1;  */

    /*reaction 39: ch3o + h <=> ch3 + oh */
    /*eqcon[38] *= 1;  */

    /*reaction 40: ch2oh + h <=> ch3 + oh */
    /*eqcon[39] *= 1;  */

    /*reaction 41: ch3 + oh <=> ch2s + h2o */
    /*eqcon[40] *= 1;  */

    /*reaction 42: ch3 + oh <=> ch2 + h2o */
    /*eqcon[41] *= 1;  */

    /*reaction 43: ch3 + h <=> ch2 + h2 */
    /*eqcon[42] *= 1;  */

    /*reaction 44: ch3 + M <=> ch + h2 + M */
    eqcon[43] *= 1e-06; 

    /*reaction 45: ch3 + M <=> ch2 + h + M */
    eqcon[44] *= 1e-06; 

    /*reaction 46: ch2o + h (+M) <=> ch3o (+M) */
    eqcon[45] *= 1e+06; 

    /*reaction 47: ch2o + h (+M) <=> ch2oh (+M) */
    eqcon[46] *= 1e+06; 

    /*reaction 48: ch3o + ch3 <=> ch2o + ch4 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: ch3o + h <=> ch2o + h2 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: ch2oh + h <=> ch2o + h2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: ch3o + oh <=> ch2o + h2o */
    /*eqcon[50] *= 1;  */

    /*reaction 52: ch2oh + oh <=> ch2o + h2o */
    /*eqcon[51] *= 1;  */

    /*reaction 53: ch3o + o <=> ch2o + oh */
    /*eqcon[52] *= 1;  */

    /*reaction 54: ch2oh + o <=> ch2o + oh */
    /*eqcon[53] *= 1;  */

    /*reaction 55: ch3o + o2 <=> ch2o + ho2 */
    /*eqcon[54] *= 1;  */

    /*reaction 56: ch3o + co <=> ch3 + co2 */
    /*eqcon[55] *= 1;  */

    /*reaction 57: ch2oh + o2 <=> ch2o + ho2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: ch2oh + o2 <=> ch2o + ho2 */
    /*eqcon[57] *= 1;  */

    /*reaction 59: ch2 + h <=> ch + h2 */
    /*eqcon[58] *= 1;  */

    /*reaction 60: ch2 + oh <=> ch + h2o */
    /*eqcon[59] *= 1;  */

    /*reaction 61: ch2 + oh <=> ch2o + h */
    /*eqcon[60] *= 1;  */

    /*reaction 62: ch2 + co2 <=> ch2o + co */
    /*eqcon[61] *= 1;  */

    /*reaction 63: ch2 + o <=> co + h + h */
    eqcon[62] *= 1e-06; 

    /*reaction 64: ch2 + o <=> co + h2 */
    /*eqcon[63] *= 1;  */

    /*reaction 65: ch2 + o2 <=> ch2o + o */
    /*eqcon[64] *= 1;  */

    /*reaction 66: ch2 + o2 <=> co2 + h + h */
    eqcon[65] *= 1e-06; 

    /*reaction 67: ch2 + o2 <=> co2 + h2 */
    /*eqcon[66] *= 1;  */

    /*reaction 68: ch2 + o2 <=> co + h2o */
    /*eqcon[67] *= 1;  */

    /*reaction 69: ch2 + o2 <=> hco + oh */
    /*eqcon[68] *= 1;  */

    /*reaction 70: ch2 + ch3 <=> c2h4 + h */
    /*eqcon[69] *= 1;  */

    /*reaction 71: ch2 + ch2 <=> c2h2 + h + h */
    eqcon[70] *= 1e-06; 

    /*reaction 72: ch2 + hcco <=> c2h3 + co */
    /*eqcon[71] *= 1;  */

    /*reaction 73: ch2 + c2h2 <=> h2ccch + h */
    /*eqcon[72] *= 1;  */

    /*reaction 74: ch2s + M <=> ch2 + M */
    /*eqcon[73] *= 1;  */

    /*reaction 75: ch2s + ch4 <=> ch3 + ch3 */
    /*eqcon[74] *= 1;  */

    /*reaction 76: ch2s + c2h6 <=> ch3 + c2h5 */
    /*eqcon[75] *= 1;  */

    /*reaction 77: ch2s + o2 <=> co + oh + h */
    eqcon[76] *= 1e-06; 

    /*reaction 78: ch2s + h2 <=> ch3 + h */
    /*eqcon[77] *= 1;  */

    /*reaction 79: ch2s + c2h2 <=> h2ccch + h */
    /*eqcon[78] *= 1;  */

    /*reaction 80: ch2s + o <=> co + h + h */
    eqcon[79] *= 1e-06; 

    /*reaction 81: ch2s + oh <=> ch2o + h */
    /*eqcon[80] *= 1;  */

    /*reaction 82: ch2s + h <=> ch + h2 */
    /*eqcon[81] *= 1;  */

    /*reaction 83: ch2s + co2 <=> ch2o + co */
    /*eqcon[82] *= 1;  */

    /*reaction 84: ch2s + ch3 <=> c2h4 + h */
    /*eqcon[83] *= 1;  */

    /*reaction 85: ch2s + ch2co <=> c2h4 + co */
    /*eqcon[84] *= 1;  */

    /*reaction 86: ch + o2 <=> hco + o */
    /*eqcon[85] *= 1;  */

    /*reaction 87: ch + o <=> co + h */
    /*eqcon[86] *= 1;  */

    /*reaction 88: ch + oh <=> hco + h */
    /*eqcon[87] *= 1;  */

    /*reaction 89: ch + co2 <=> hco + co */
    /*eqcon[88] *= 1;  */

    /*reaction 90: ch + h2o <=> ch2o + h */
    /*eqcon[89] *= 1;  */

    /*reaction 91: ch + ch2o <=> ch2co + h */
    /*eqcon[90] *= 1;  */

    /*reaction 92: ch + c2h2 <=> c3h2 + h */
    /*eqcon[91] *= 1;  */

    /*reaction 93: ch + ch2 <=> c2h2 + h */
    /*eqcon[92] *= 1;  */

    /*reaction 94: ch + ch3 <=> c2h3 + h */
    /*eqcon[93] *= 1;  */

    /*reaction 95: ch + ch4 <=> c2h4 + h */
    /*eqcon[94] *= 1;  */

    /*reaction 96: ch2o + oh <=> hco + h2o */
    /*eqcon[95] *= 1;  */

    /*reaction 97: ch2o + h <=> hco + h2 */
    /*eqcon[96] *= 1;  */

    /*reaction 98: ch2o + M <=> hco + h + M */
    eqcon[97] *= 1e-06; 

    /*reaction 99: ch2o + o <=> hco + oh */
    /*eqcon[98] *= 1;  */

    /*reaction 100: hco + o2 <=> co + ho2 */
    /*eqcon[99] *= 1;  */

    /*reaction 101: hco + M <=> h + co + M */
    eqcon[100] *= 1e-06; 

    /*reaction 102: hco + oh <=> h2o + co */
    /*eqcon[101] *= 1;  */

    /*reaction 103: hco + h <=> co + h2 */
    /*eqcon[102] *= 1;  */

    /*reaction 104: hco + o <=> co + oh */
    /*eqcon[103] *= 1;  */

    /*reaction 105: hco + o <=> co2 + h */
    /*eqcon[104] *= 1;  */

    /*reaction 106: co + oh <=> co2 + h */
    /*eqcon[105] *= 1;  */

    /*reaction 107: co + o + M <=> co2 + M */
    eqcon[106] *= 1e+06; 

    /*reaction 108: co + o2 <=> co2 + o */
    /*eqcon[107] *= 1;  */

    /*reaction 109: co + ho2 <=> co2 + oh */
    /*eqcon[108] *= 1;  */

    /*reaction 110: c2h5oh (+M) <=> ch3 + ch2oh (+M) */
    eqcon[109] *= 1e-06; 

    /*reaction 111: c2h5oh (+M) <=> c2h5 + oh (+M) */
    eqcon[110] *= 1e-06; 

    /*reaction 112: c2h5oh (+M) <=> c2h4 + h2o (+M) */
    eqcon[111] *= 1e-06; 

    /*reaction 113: c2h5oh (+M) <=> ch3hco + h2 (+M) */
    eqcon[112] *= 1e-06; 

    /*reaction 114: c2h5oh + oh <=> c2h4oh + h2o */
    /*eqcon[113] *= 1;  */

    /*reaction 115: c2h5oh + oh <=> ch3choh + h2o */
    /*eqcon[114] *= 1;  */

    /*reaction 116: c2h5oh + oh <=> ch3ch2o + h2o */
    /*eqcon[115] *= 1;  */

    /*reaction 117: c2h5oh + h <=> c2h4oh + h2 */
    /*eqcon[116] *= 1;  */

    /*reaction 118: c2h5oh + h <=> ch3choh + h2 */
    /*eqcon[117] *= 1;  */

    /*reaction 119: c2h5oh + h <=> ch3ch2o + h2 */
    /*eqcon[118] *= 1;  */

    /*reaction 120: c2h5oh + o <=> c2h4oh + oh */
    /*eqcon[119] *= 1;  */

    /*reaction 121: c2h5oh + o <=> ch3choh + oh */
    /*eqcon[120] *= 1;  */

    /*reaction 122: c2h5oh + o <=> ch3ch2o + oh */
    /*eqcon[121] *= 1;  */

    /*reaction 123: c2h5oh + ch3 <=> c2h4oh + ch4 */
    /*eqcon[122] *= 1;  */

    /*reaction 124: c2h5oh + ch3 <=> ch3choh + ch4 */
    /*eqcon[123] *= 1;  */

    /*reaction 125: c2h5oh + ch3 <=> ch3ch2o + ch4 */
    /*eqcon[124] *= 1;  */

    /*reaction 126: c2h5oh + ho2 <=> ch3choh + h2o2 */
    /*eqcon[125] *= 1;  */

    /*reaction 127: c2h5oh + ho2 <=> c2h4oh + h2o2 */
    /*eqcon[126] *= 1;  */

    /*reaction 128: c2h5oh + ho2 <=> ch3ch2o + h2o2 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: ch3ch2o + M <=> ch3hco + h + M */
    eqcon[128] *= 1e-06; 

    /*reaction 130: ch3ch2o + M <=> ch3 + ch2o + M */
    eqcon[129] *= 1e-06; 

    /*reaction 131: ch3ch2o + o2 <=> ch3hco + ho2 */
    /*eqcon[130] *= 1;  */

    /*reaction 132: ch3ch2o + co <=> c2h5 + co2 */
    /*eqcon[131] *= 1;  */

    /*reaction 133: ch3ch2o + h <=> ch3 + ch2oh */
    /*eqcon[132] *= 1;  */

    /*reaction 134: ch3ch2o + h <=> c2h4 + h2o */
    /*eqcon[133] *= 1;  */

    /*reaction 135: ch3ch2o + oh <=> ch3hco + h2o */
    /*eqcon[134] *= 1;  */

    /*reaction 136: ch3choh + o2 <=> ch3hco + ho2 */
    /*eqcon[135] *= 1;  */

    /*reaction 137: ch3choh + o2 <=> ch3hco + ho2 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: ch3choh + o <=> ch3hco + oh */
    /*eqcon[137] *= 1;  */

    /*reaction 139: ch3choh + h <=> c2h4 + h2o */
    /*eqcon[138] *= 1;  */

    /*reaction 140: ch3choh + h <=> ch3 + ch2oh */
    /*eqcon[139] *= 1;  */

    /*reaction 141: ch3choh + ho2 <=> ch3hco + oh + oh */
    eqcon[140] *= 1e-06; 

    /*reaction 142: ch3choh + oh <=> ch3hco + h2o */
    /*eqcon[141] *= 1;  */

    /*reaction 143: ch3choh + M <=> ch3hco + h + M */
    eqcon[142] *= 1e-06; 

    /*reaction 144: ch3hco + oh <=> ch3co + h2o */
    /*eqcon[143] *= 1;  */

    /*reaction 145: ch3hco + oh <=> ch2hco + h2o */
    /*eqcon[144] *= 1;  */

    /*reaction 146: ch3hco + o <=> ch3co + oh */
    /*eqcon[145] *= 1;  */

    /*reaction 147: ch3hco + o <=> ch2hco + oh */
    /*eqcon[146] *= 1;  */

    /*reaction 148: ch3hco + h <=> ch3co + h2 */
    /*eqcon[147] *= 1;  */

    /*reaction 149: ch3hco + h <=> ch2hco + h2 */
    /*eqcon[148] *= 1;  */

    /*reaction 150: ch3hco + ch3 <=> ch3co + ch4 */
    /*eqcon[149] *= 1;  */

    /*reaction 151: ch3hco + ch3 <=> ch2hco + ch4 */
    /*eqcon[150] *= 1;  */

    /*reaction 152: ch3hco + ho2 <=> ch3co + h2o2 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: ch3hco + ho2 <=> ch2hco + h2o2 */
    /*eqcon[152] *= 1;  */

    /*reaction 154: ch3hco + o2 <=> ch3co + ho2 */
    /*eqcon[153] *= 1;  */

    /*reaction 155: c2h6 + ch3 <=> c2h5 + ch4 */
    /*eqcon[154] *= 1;  */

    /*reaction 156: c2h6 + h <=> c2h5 + h2 */
    /*eqcon[155] *= 1;  */

    /*reaction 157: c2h6 + o <=> c2h5 + oh */
    /*eqcon[156] *= 1;  */

    /*reaction 158: c2h6 + oh <=> c2h5 + h2o */
    /*eqcon[157] *= 1;  */

    /*reaction 159: c2h5 + h <=> c2h4 + h2 */
    /*eqcon[158] *= 1;  */

    /*reaction 160: c2h5 + h <=> ch3 + ch3 */
    /*eqcon[159] *= 1;  */

    /*reaction 161: c2h5 + h <=> c2h6 */
    eqcon[160] *= 1e+06; 

    /*reaction 162: c2h5 + oh <=> c2h4 + h2o */
    /*eqcon[161] *= 1;  */

    /*reaction 163: c2h5 + o <=> ch3 + ch2o */
    /*eqcon[162] *= 1;  */

    /*reaction 164: c2h5 + ho2 <=> c2h6 + o2 */
    /*eqcon[163] *= 1;  */

    /*reaction 165: c2h5 + ho2 <=> ch3ch2o + oh */
    /*eqcon[164] *= 1;  */

    /*reaction 166: c2h5 + o2 <=> c2h4 + ho2 */
    /*eqcon[165] *= 1;  */

    /*reaction 167: c2h5 + o2 <=> ch3hco + oh */
    /*eqcon[166] *= 1;  */

    /*reaction 168: c2h4 + oh <=> c2h4oh */
    eqcon[167] *= 1e+06; 

    /*reaction 169: c2h4 + oh <=> c2h3 + h2o */
    /*eqcon[168] *= 1;  */

    /*reaction 170: c2h4 + o <=> ch3 + hco */
    /*eqcon[169] *= 1;  */

    /*reaction 171: c2h4 + o <=> ch2hco + h */
    /*eqcon[170] *= 1;  */

    /*reaction 172: c2h4 + ch3 <=> c2h3 + ch4 */
    /*eqcon[171] *= 1;  */

    /*reaction 173: c2h4 + h <=> c2h3 + h2 */
    /*eqcon[172] *= 1;  */

    /*reaction 174: c2h4 + h (+M) <=> c2h5 (+M) */
    eqcon[173] *= 1e+06; 

    /*reaction 175: c2h4 (+M) <=> c2h2 + h2 (+M) */
    eqcon[174] *= 1e-06; 

    /*reaction 176: c2h3 + h (+M) <=> c2h4 (+M) */
    eqcon[175] *= 1e+06; 

    /*reaction 177: c2h3 + h <=> c2h2 + h2 */
    /*eqcon[176] *= 1;  */

    /*reaction 178: c2h3 + o <=> ch2co + h */
    /*eqcon[177] *= 1;  */

    /*reaction 179: c2h3 + o2 <=> ch2o + hco */
    /*eqcon[178] *= 1;  */

    /*reaction 180: c2h3 + o2 <=> ch2hco + o */
    /*eqcon[179] *= 1;  */

    /*reaction 181: c2h3 + o2 <=> c2h2 + ho2 */
    /*eqcon[180] *= 1;  */

    /*reaction 182: c2h3 + oh <=> c2h2 + h2o */
    /*eqcon[181] *= 1;  */

    /*reaction 183: c2h3 + c2h <=> c2h2 + c2h2 */
    /*eqcon[182] *= 1;  */

    /*reaction 184: c2h3 + ch <=> ch2 + c2h2 */
    /*eqcon[183] *= 1;  */

    /*reaction 185: c2h3 + ch3 <=> c2h2 + ch4 */
    /*eqcon[184] *= 1;  */

    /*reaction 186: c2h2 + oh <=> c2h + h2o */
    /*eqcon[185] *= 1;  */

    /*reaction 187: c2h2 + oh <=> hccoh + h */
    /*eqcon[186] *= 1;  */

    /*reaction 188: c2h2 + oh <=> ch2co + h */
    /*eqcon[187] *= 1;  */

    /*reaction 189: c2h2 + oh <=> ch2co + h */
    /*eqcon[188] *= 1;  */

    /*reaction 190: c2h2 + oh <=> ch3 + co */
    /*eqcon[189] *= 1;  */

    /*reaction 191: hccoh + h <=> ch2co + h */
    /*eqcon[190] *= 1;  */

    /*reaction 192: c2h2 + o <=> ch2 + co */
    /*eqcon[191] *= 1;  */

    /*reaction 193: c2h2 + o <=> hcco + h */
    /*eqcon[192] *= 1;  */

    /*reaction 194: c2h2 + o <=> c2h + oh */
    /*eqcon[193] *= 1;  */

    /*reaction 195: c2h2 + ch3 <=> c2h + ch4 */
    /*eqcon[194] *= 1;  */

    /*reaction 196: c2h2 + o2 <=> hcco + oh */
    /*eqcon[195] *= 1;  */

    /*reaction 197: c2h2 + M <=> c2h + h + M */
    eqcon[196] *= 1e-06; 

    /*reaction 198: c2h2 + h (+M) <=> c2h3 (+M) */
    eqcon[197] *= 1e+06; 

    /*reaction 199: ch2hco + h <=> ch3 + hco */
    /*eqcon[198] *= 1;  */

    /*reaction 200: ch2hco + h <=> ch2co + h2 */
    /*eqcon[199] *= 1;  */

    /*reaction 201: ch2hco + o <=> ch2o + hco */
    /*eqcon[200] *= 1;  */

    /*reaction 202: ch2hco + oh <=> ch2co + h2o */
    /*eqcon[201] *= 1;  */

    /*reaction 203: ch2hco + o2 <=> ch2o + co + oh */
    eqcon[202] *= 1e-06; 

    /*reaction 204: ch2hco + ch3 <=> c2h5 + co + h */
    eqcon[203] *= 1e-06; 

    /*reaction 205: ch2hco + ho2 <=> ch2o + hco + oh */
    eqcon[204] *= 1e-06; 

    /*reaction 206: ch2hco + ho2 <=> ch3hco + o2 */
    /*eqcon[205] *= 1;  */

    /*reaction 207: ch2hco <=> ch3 + co */
    eqcon[206] *= 1e-06; 

    /*reaction 208: ch2hco <=> ch2co + h */
    eqcon[207] *= 1e-06; 

    /*reaction 209: ch3co (+M) <=> ch3 + co (+M) */
    eqcon[208] *= 1e-06; 

    /*reaction 210: ch2co + o <=> co2 + ch2 */
    /*eqcon[209] *= 1;  */

    /*reaction 211: ch2co + h <=> ch3 + co */
    /*eqcon[210] *= 1;  */

    /*reaction 212: ch2co + h <=> hcco + h2 */
    /*eqcon[211] *= 1;  */

    /*reaction 213: ch2co + o <=> hcco + oh */
    /*eqcon[212] *= 1;  */

    /*reaction 214: ch2co + oh <=> hcco + h2o */
    /*eqcon[213] *= 1;  */

    /*reaction 215: ch2co + oh <=> ch2oh + co */
    /*eqcon[214] *= 1;  */

    /*reaction 216: ch2co (+M) <=> ch2 + co (+M) */
    eqcon[215] *= 1e-06; 

    /*reaction 217: c2h + h2 <=> c2h2 + h */
    /*eqcon[216] *= 1;  */

    /*reaction 218: c2h + o <=> ch + co */
    /*eqcon[217] *= 1;  */

    /*reaction 219: c2h + oh <=> hcco + h */
    /*eqcon[218] *= 1;  */

    /*reaction 220: c2h + o2 <=> co + co + h */
    eqcon[219] *= 1e-06; 

    /*reaction 221: hcco + c2h2 <=> h2ccch + co */
    /*eqcon[220] *= 1;  */

    /*reaction 222: hcco + h <=> ch2s + co */
    /*eqcon[221] *= 1;  */

    /*reaction 223: hcco + o <=> h + co + co */
    eqcon[222] *= 1e-06; 

    /*reaction 224: hcco + o <=> ch + co2 */
    /*eqcon[223] *= 1;  */

    /*reaction 225: hcco + o2 <=> hco + co + o */
    eqcon[224] *= 1e-06; 

    /*reaction 226: hcco + o2 <=> co2 + hco */
    /*eqcon[225] *= 1;  */

    /*reaction 227: hcco + ch <=> c2h2 + co */
    /*eqcon[226] *= 1;  */

    /*reaction 228: hcco + hcco <=> c2h2 + co + co */
    eqcon[227] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[38]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: oh + h2 <=> h + h2o */
    /*eqcon[0] *= 1;  */

    /*reaction 2: o + oh <=> o2 + h */
    /*eqcon[1] *= 1;  */

    /*reaction 3: o + h2 <=> oh + h */
    /*eqcon[2] *= 1;  */

    /*reaction 4: h + o2 (+M) <=> ho2 (+M) */
    eqcon[3] *= 1e+06; 

    /*reaction 5: h + o2 (+n2) <=> ho2 (+n2) */
    eqcon[4] *= 1e+06; 

    /*reaction 6: h + o2 (+h2) <=> ho2 (+h2) */
    eqcon[5] *= 1e+06; 

    /*reaction 7: h + o2 (+h2o) <=> ho2 (+h2o) */
    eqcon[6] *= 1e+06; 

    /*reaction 8: oh + ho2 <=> h2o + o2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: oh + ho2 <=> h2o + o2 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: h + ho2 <=> oh + oh */
    /*eqcon[9] *= 1;  */

    /*reaction 11: h + ho2 <=> h2 + o2 */
    /*eqcon[10] *= 1;  */

    /*reaction 12: h + ho2 <=> o + h2o */
    /*eqcon[11] *= 1;  */

    /*reaction 13: o + ho2 <=> o2 + oh */
    /*eqcon[12] *= 1;  */

    /*reaction 14: oh + oh <=> o + h2o */
    /*eqcon[13] *= 1;  */

    /*reaction 15: h + h + M <=> h2 + M */
    eqcon[14] *= 1e+06; 

    /*reaction 16: h + h + h2 <=> h2 + h2 */
    eqcon[15] *= 1e+06; 

    /*reaction 17: h + h + h2o <=> h2 + h2o */
    eqcon[16] *= 1e+06; 

    /*reaction 18: h + oh + M <=> h2o + M */
    eqcon[17] *= 1e+06; 

    /*reaction 19: h + o + M <=> oh + M */
    eqcon[18] *= 1e+06; 

    /*reaction 20: o + o + M <=> o2 + M */
    eqcon[19] *= 1e+06; 

    /*reaction 21: ho2 + ho2 <=> h2o2 + o2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: ho2 + ho2 <=> h2o2 + o2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: oh + oh (+M) <=> h2o2 (+M) */
    eqcon[22] *= 1e+06; 

    /*reaction 24: h2o2 + h <=> ho2 + h2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: h2o2 + h <=> oh + h2o */
    /*eqcon[24] *= 1;  */

    /*reaction 26: h2o2 + o <=> oh + ho2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: h2o2 + oh <=> h2o + ho2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: ch3 + ch3 (+M) <=> c2h6 (+M) */
    eqcon[27] *= 1e+06; 

    /*reaction 29: ch3 + h (+M) <=> ch4 (+M) */
    eqcon[28] *= 1e+06; 

    /*reaction 30: ch4 + h <=> ch3 + h2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: ch4 + oh <=> ch3 + h2o */
    /*eqcon[30] *= 1;  */

    /*reaction 32: ch4 + o <=> ch3 + oh */
    /*eqcon[31] *= 1;  */

    /*reaction 33: ch4 + ho2 <=> ch3 + h2o2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: ch3 + ho2 <=> ch3o + oh */
    /*eqcon[33] *= 1;  */

    /*reaction 35: ch3 + ho2 <=> ch4 + o2 */
    /*eqcon[34] *= 1;  */

    /*reaction 36: ch3 + o <=> ch2o + h */
    /*eqcon[35] *= 1;  */

    /*reaction 37: ch3 + o2 <=> ch3o + o */
    /*eqcon[36] *= 1;  */

    /*reaction 38: ch3 + o2 <=> ch2o + oh */
    /*eqcon[37] *= 1;  */

    /*reaction 39: ch3o + h <=> ch3 + oh */
    /*eqcon[38] *= 1;  */

    /*reaction 40: ch2oh + h <=> ch3 + oh */
    /*eqcon[39] *= 1;  */

    /*reaction 41: ch3 + oh <=> ch2s + h2o */
    /*eqcon[40] *= 1;  */

    /*reaction 42: ch3 + oh <=> ch2 + h2o */
    /*eqcon[41] *= 1;  */

    /*reaction 43: ch3 + h <=> ch2 + h2 */
    /*eqcon[42] *= 1;  */

    /*reaction 44: ch3 + M <=> ch + h2 + M */
    eqcon[43] *= 1e-06; 

    /*reaction 45: ch3 + M <=> ch2 + h + M */
    eqcon[44] *= 1e-06; 

    /*reaction 46: ch2o + h (+M) <=> ch3o (+M) */
    eqcon[45] *= 1e+06; 

    /*reaction 47: ch2o + h (+M) <=> ch2oh (+M) */
    eqcon[46] *= 1e+06; 

    /*reaction 48: ch3o + ch3 <=> ch2o + ch4 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: ch3o + h <=> ch2o + h2 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: ch2oh + h <=> ch2o + h2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: ch3o + oh <=> ch2o + h2o */
    /*eqcon[50] *= 1;  */

    /*reaction 52: ch2oh + oh <=> ch2o + h2o */
    /*eqcon[51] *= 1;  */

    /*reaction 53: ch3o + o <=> ch2o + oh */
    /*eqcon[52] *= 1;  */

    /*reaction 54: ch2oh + o <=> ch2o + oh */
    /*eqcon[53] *= 1;  */

    /*reaction 55: ch3o + o2 <=> ch2o + ho2 */
    /*eqcon[54] *= 1;  */

    /*reaction 56: ch3o + co <=> ch3 + co2 */
    /*eqcon[55] *= 1;  */

    /*reaction 57: ch2oh + o2 <=> ch2o + ho2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: ch2oh + o2 <=> ch2o + ho2 */
    /*eqcon[57] *= 1;  */

    /*reaction 59: ch2 + h <=> ch + h2 */
    /*eqcon[58] *= 1;  */

    /*reaction 60: ch2 + oh <=> ch + h2o */
    /*eqcon[59] *= 1;  */

    /*reaction 61: ch2 + oh <=> ch2o + h */
    /*eqcon[60] *= 1;  */

    /*reaction 62: ch2 + co2 <=> ch2o + co */
    /*eqcon[61] *= 1;  */

    /*reaction 63: ch2 + o <=> co + h + h */
    eqcon[62] *= 1e-06; 

    /*reaction 64: ch2 + o <=> co + h2 */
    /*eqcon[63] *= 1;  */

    /*reaction 65: ch2 + o2 <=> ch2o + o */
    /*eqcon[64] *= 1;  */

    /*reaction 66: ch2 + o2 <=> co2 + h + h */
    eqcon[65] *= 1e-06; 

    /*reaction 67: ch2 + o2 <=> co2 + h2 */
    /*eqcon[66] *= 1;  */

    /*reaction 68: ch2 + o2 <=> co + h2o */
    /*eqcon[67] *= 1;  */

    /*reaction 69: ch2 + o2 <=> hco + oh */
    /*eqcon[68] *= 1;  */

    /*reaction 70: ch2 + ch3 <=> c2h4 + h */
    /*eqcon[69] *= 1;  */

    /*reaction 71: ch2 + ch2 <=> c2h2 + h + h */
    eqcon[70] *= 1e-06; 

    /*reaction 72: ch2 + hcco <=> c2h3 + co */
    /*eqcon[71] *= 1;  */

    /*reaction 73: ch2 + c2h2 <=> h2ccch + h */
    /*eqcon[72] *= 1;  */

    /*reaction 74: ch2s + M <=> ch2 + M */
    /*eqcon[73] *= 1;  */

    /*reaction 75: ch2s + ch4 <=> ch3 + ch3 */
    /*eqcon[74] *= 1;  */

    /*reaction 76: ch2s + c2h6 <=> ch3 + c2h5 */
    /*eqcon[75] *= 1;  */

    /*reaction 77: ch2s + o2 <=> co + oh + h */
    eqcon[76] *= 1e-06; 

    /*reaction 78: ch2s + h2 <=> ch3 + h */
    /*eqcon[77] *= 1;  */

    /*reaction 79: ch2s + c2h2 <=> h2ccch + h */
    /*eqcon[78] *= 1;  */

    /*reaction 80: ch2s + o <=> co + h + h */
    eqcon[79] *= 1e-06; 

    /*reaction 81: ch2s + oh <=> ch2o + h */
    /*eqcon[80] *= 1;  */

    /*reaction 82: ch2s + h <=> ch + h2 */
    /*eqcon[81] *= 1;  */

    /*reaction 83: ch2s + co2 <=> ch2o + co */
    /*eqcon[82] *= 1;  */

    /*reaction 84: ch2s + ch3 <=> c2h4 + h */
    /*eqcon[83] *= 1;  */

    /*reaction 85: ch2s + ch2co <=> c2h4 + co */
    /*eqcon[84] *= 1;  */

    /*reaction 86: ch + o2 <=> hco + o */
    /*eqcon[85] *= 1;  */

    /*reaction 87: ch + o <=> co + h */
    /*eqcon[86] *= 1;  */

    /*reaction 88: ch + oh <=> hco + h */
    /*eqcon[87] *= 1;  */

    /*reaction 89: ch + co2 <=> hco + co */
    /*eqcon[88] *= 1;  */

    /*reaction 90: ch + h2o <=> ch2o + h */
    /*eqcon[89] *= 1;  */

    /*reaction 91: ch + ch2o <=> ch2co + h */
    /*eqcon[90] *= 1;  */

    /*reaction 92: ch + c2h2 <=> c3h2 + h */
    /*eqcon[91] *= 1;  */

    /*reaction 93: ch + ch2 <=> c2h2 + h */
    /*eqcon[92] *= 1;  */

    /*reaction 94: ch + ch3 <=> c2h3 + h */
    /*eqcon[93] *= 1;  */

    /*reaction 95: ch + ch4 <=> c2h4 + h */
    /*eqcon[94] *= 1;  */

    /*reaction 96: ch2o + oh <=> hco + h2o */
    /*eqcon[95] *= 1;  */

    /*reaction 97: ch2o + h <=> hco + h2 */
    /*eqcon[96] *= 1;  */

    /*reaction 98: ch2o + M <=> hco + h + M */
    eqcon[97] *= 1e-06; 

    /*reaction 99: ch2o + o <=> hco + oh */
    /*eqcon[98] *= 1;  */

    /*reaction 100: hco + o2 <=> co + ho2 */
    /*eqcon[99] *= 1;  */

    /*reaction 101: hco + M <=> h + co + M */
    eqcon[100] *= 1e-06; 

    /*reaction 102: hco + oh <=> h2o + co */
    /*eqcon[101] *= 1;  */

    /*reaction 103: hco + h <=> co + h2 */
    /*eqcon[102] *= 1;  */

    /*reaction 104: hco + o <=> co + oh */
    /*eqcon[103] *= 1;  */

    /*reaction 105: hco + o <=> co2 + h */
    /*eqcon[104] *= 1;  */

    /*reaction 106: co + oh <=> co2 + h */
    /*eqcon[105] *= 1;  */

    /*reaction 107: co + o + M <=> co2 + M */
    eqcon[106] *= 1e+06; 

    /*reaction 108: co + o2 <=> co2 + o */
    /*eqcon[107] *= 1;  */

    /*reaction 109: co + ho2 <=> co2 + oh */
    /*eqcon[108] *= 1;  */

    /*reaction 110: c2h5oh (+M) <=> ch3 + ch2oh (+M) */
    eqcon[109] *= 1e-06; 

    /*reaction 111: c2h5oh (+M) <=> c2h5 + oh (+M) */
    eqcon[110] *= 1e-06; 

    /*reaction 112: c2h5oh (+M) <=> c2h4 + h2o (+M) */
    eqcon[111] *= 1e-06; 

    /*reaction 113: c2h5oh (+M) <=> ch3hco + h2 (+M) */
    eqcon[112] *= 1e-06; 

    /*reaction 114: c2h5oh + oh <=> c2h4oh + h2o */
    /*eqcon[113] *= 1;  */

    /*reaction 115: c2h5oh + oh <=> ch3choh + h2o */
    /*eqcon[114] *= 1;  */

    /*reaction 116: c2h5oh + oh <=> ch3ch2o + h2o */
    /*eqcon[115] *= 1;  */

    /*reaction 117: c2h5oh + h <=> c2h4oh + h2 */
    /*eqcon[116] *= 1;  */

    /*reaction 118: c2h5oh + h <=> ch3choh + h2 */
    /*eqcon[117] *= 1;  */

    /*reaction 119: c2h5oh + h <=> ch3ch2o + h2 */
    /*eqcon[118] *= 1;  */

    /*reaction 120: c2h5oh + o <=> c2h4oh + oh */
    /*eqcon[119] *= 1;  */

    /*reaction 121: c2h5oh + o <=> ch3choh + oh */
    /*eqcon[120] *= 1;  */

    /*reaction 122: c2h5oh + o <=> ch3ch2o + oh */
    /*eqcon[121] *= 1;  */

    /*reaction 123: c2h5oh + ch3 <=> c2h4oh + ch4 */
    /*eqcon[122] *= 1;  */

    /*reaction 124: c2h5oh + ch3 <=> ch3choh + ch4 */
    /*eqcon[123] *= 1;  */

    /*reaction 125: c2h5oh + ch3 <=> ch3ch2o + ch4 */
    /*eqcon[124] *= 1;  */

    /*reaction 126: c2h5oh + ho2 <=> ch3choh + h2o2 */
    /*eqcon[125] *= 1;  */

    /*reaction 127: c2h5oh + ho2 <=> c2h4oh + h2o2 */
    /*eqcon[126] *= 1;  */

    /*reaction 128: c2h5oh + ho2 <=> ch3ch2o + h2o2 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: ch3ch2o + M <=> ch3hco + h + M */
    eqcon[128] *= 1e-06; 

    /*reaction 130: ch3ch2o + M <=> ch3 + ch2o + M */
    eqcon[129] *= 1e-06; 

    /*reaction 131: ch3ch2o + o2 <=> ch3hco + ho2 */
    /*eqcon[130] *= 1;  */

    /*reaction 132: ch3ch2o + co <=> c2h5 + co2 */
    /*eqcon[131] *= 1;  */

    /*reaction 133: ch3ch2o + h <=> ch3 + ch2oh */
    /*eqcon[132] *= 1;  */

    /*reaction 134: ch3ch2o + h <=> c2h4 + h2o */
    /*eqcon[133] *= 1;  */

    /*reaction 135: ch3ch2o + oh <=> ch3hco + h2o */
    /*eqcon[134] *= 1;  */

    /*reaction 136: ch3choh + o2 <=> ch3hco + ho2 */
    /*eqcon[135] *= 1;  */

    /*reaction 137: ch3choh + o2 <=> ch3hco + ho2 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: ch3choh + o <=> ch3hco + oh */
    /*eqcon[137] *= 1;  */

    /*reaction 139: ch3choh + h <=> c2h4 + h2o */
    /*eqcon[138] *= 1;  */

    /*reaction 140: ch3choh + h <=> ch3 + ch2oh */
    /*eqcon[139] *= 1;  */

    /*reaction 141: ch3choh + ho2 <=> ch3hco + oh + oh */
    eqcon[140] *= 1e-06; 

    /*reaction 142: ch3choh + oh <=> ch3hco + h2o */
    /*eqcon[141] *= 1;  */

    /*reaction 143: ch3choh + M <=> ch3hco + h + M */
    eqcon[142] *= 1e-06; 

    /*reaction 144: ch3hco + oh <=> ch3co + h2o */
    /*eqcon[143] *= 1;  */

    /*reaction 145: ch3hco + oh <=> ch2hco + h2o */
    /*eqcon[144] *= 1;  */

    /*reaction 146: ch3hco + o <=> ch3co + oh */
    /*eqcon[145] *= 1;  */

    /*reaction 147: ch3hco + o <=> ch2hco + oh */
    /*eqcon[146] *= 1;  */

    /*reaction 148: ch3hco + h <=> ch3co + h2 */
    /*eqcon[147] *= 1;  */

    /*reaction 149: ch3hco + h <=> ch2hco + h2 */
    /*eqcon[148] *= 1;  */

    /*reaction 150: ch3hco + ch3 <=> ch3co + ch4 */
    /*eqcon[149] *= 1;  */

    /*reaction 151: ch3hco + ch3 <=> ch2hco + ch4 */
    /*eqcon[150] *= 1;  */

    /*reaction 152: ch3hco + ho2 <=> ch3co + h2o2 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: ch3hco + ho2 <=> ch2hco + h2o2 */
    /*eqcon[152] *= 1;  */

    /*reaction 154: ch3hco + o2 <=> ch3co + ho2 */
    /*eqcon[153] *= 1;  */

    /*reaction 155: c2h6 + ch3 <=> c2h5 + ch4 */
    /*eqcon[154] *= 1;  */

    /*reaction 156: c2h6 + h <=> c2h5 + h2 */
    /*eqcon[155] *= 1;  */

    /*reaction 157: c2h6 + o <=> c2h5 + oh */
    /*eqcon[156] *= 1;  */

    /*reaction 158: c2h6 + oh <=> c2h5 + h2o */
    /*eqcon[157] *= 1;  */

    /*reaction 159: c2h5 + h <=> c2h4 + h2 */
    /*eqcon[158] *= 1;  */

    /*reaction 160: c2h5 + h <=> ch3 + ch3 */
    /*eqcon[159] *= 1;  */

    /*reaction 161: c2h5 + h <=> c2h6 */
    eqcon[160] *= 1e+06; 

    /*reaction 162: c2h5 + oh <=> c2h4 + h2o */
    /*eqcon[161] *= 1;  */

    /*reaction 163: c2h5 + o <=> ch3 + ch2o */
    /*eqcon[162] *= 1;  */

    /*reaction 164: c2h5 + ho2 <=> c2h6 + o2 */
    /*eqcon[163] *= 1;  */

    /*reaction 165: c2h5 + ho2 <=> ch3ch2o + oh */
    /*eqcon[164] *= 1;  */

    /*reaction 166: c2h5 + o2 <=> c2h4 + ho2 */
    /*eqcon[165] *= 1;  */

    /*reaction 167: c2h5 + o2 <=> ch3hco + oh */
    /*eqcon[166] *= 1;  */

    /*reaction 168: c2h4 + oh <=> c2h4oh */
    eqcon[167] *= 1e+06; 

    /*reaction 169: c2h4 + oh <=> c2h3 + h2o */
    /*eqcon[168] *= 1;  */

    /*reaction 170: c2h4 + o <=> ch3 + hco */
    /*eqcon[169] *= 1;  */

    /*reaction 171: c2h4 + o <=> ch2hco + h */
    /*eqcon[170] *= 1;  */

    /*reaction 172: c2h4 + ch3 <=> c2h3 + ch4 */
    /*eqcon[171] *= 1;  */

    /*reaction 173: c2h4 + h <=> c2h3 + h2 */
    /*eqcon[172] *= 1;  */

    /*reaction 174: c2h4 + h (+M) <=> c2h5 (+M) */
    eqcon[173] *= 1e+06; 

    /*reaction 175: c2h4 (+M) <=> c2h2 + h2 (+M) */
    eqcon[174] *= 1e-06; 

    /*reaction 176: c2h3 + h (+M) <=> c2h4 (+M) */
    eqcon[175] *= 1e+06; 

    /*reaction 177: c2h3 + h <=> c2h2 + h2 */
    /*eqcon[176] *= 1;  */

    /*reaction 178: c2h3 + o <=> ch2co + h */
    /*eqcon[177] *= 1;  */

    /*reaction 179: c2h3 + o2 <=> ch2o + hco */
    /*eqcon[178] *= 1;  */

    /*reaction 180: c2h3 + o2 <=> ch2hco + o */
    /*eqcon[179] *= 1;  */

    /*reaction 181: c2h3 + o2 <=> c2h2 + ho2 */
    /*eqcon[180] *= 1;  */

    /*reaction 182: c2h3 + oh <=> c2h2 + h2o */
    /*eqcon[181] *= 1;  */

    /*reaction 183: c2h3 + c2h <=> c2h2 + c2h2 */
    /*eqcon[182] *= 1;  */

    /*reaction 184: c2h3 + ch <=> ch2 + c2h2 */
    /*eqcon[183] *= 1;  */

    /*reaction 185: c2h3 + ch3 <=> c2h2 + ch4 */
    /*eqcon[184] *= 1;  */

    /*reaction 186: c2h2 + oh <=> c2h + h2o */
    /*eqcon[185] *= 1;  */

    /*reaction 187: c2h2 + oh <=> hccoh + h */
    /*eqcon[186] *= 1;  */

    /*reaction 188: c2h2 + oh <=> ch2co + h */
    /*eqcon[187] *= 1;  */

    /*reaction 189: c2h2 + oh <=> ch2co + h */
    /*eqcon[188] *= 1;  */

    /*reaction 190: c2h2 + oh <=> ch3 + co */
    /*eqcon[189] *= 1;  */

    /*reaction 191: hccoh + h <=> ch2co + h */
    /*eqcon[190] *= 1;  */

    /*reaction 192: c2h2 + o <=> ch2 + co */
    /*eqcon[191] *= 1;  */

    /*reaction 193: c2h2 + o <=> hcco + h */
    /*eqcon[192] *= 1;  */

    /*reaction 194: c2h2 + o <=> c2h + oh */
    /*eqcon[193] *= 1;  */

    /*reaction 195: c2h2 + ch3 <=> c2h + ch4 */
    /*eqcon[194] *= 1;  */

    /*reaction 196: c2h2 + o2 <=> hcco + oh */
    /*eqcon[195] *= 1;  */

    /*reaction 197: c2h2 + M <=> c2h + h + M */
    eqcon[196] *= 1e-06; 

    /*reaction 198: c2h2 + h (+M) <=> c2h3 (+M) */
    eqcon[197] *= 1e+06; 

    /*reaction 199: ch2hco + h <=> ch3 + hco */
    /*eqcon[198] *= 1;  */

    /*reaction 200: ch2hco + h <=> ch2co + h2 */
    /*eqcon[199] *= 1;  */

    /*reaction 201: ch2hco + o <=> ch2o + hco */
    /*eqcon[200] *= 1;  */

    /*reaction 202: ch2hco + oh <=> ch2co + h2o */
    /*eqcon[201] *= 1;  */

    /*reaction 203: ch2hco + o2 <=> ch2o + co + oh */
    eqcon[202] *= 1e-06; 

    /*reaction 204: ch2hco + ch3 <=> c2h5 + co + h */
    eqcon[203] *= 1e-06; 

    /*reaction 205: ch2hco + ho2 <=> ch2o + hco + oh */
    eqcon[204] *= 1e-06; 

    /*reaction 206: ch2hco + ho2 <=> ch3hco + o2 */
    /*eqcon[205] *= 1;  */

    /*reaction 207: ch2hco <=> ch3 + co */
    eqcon[206] *= 1e-06; 

    /*reaction 208: ch2hco <=> ch2co + h */
    eqcon[207] *= 1e-06; 

    /*reaction 209: ch3co (+M) <=> ch3 + co (+M) */
    eqcon[208] *= 1e-06; 

    /*reaction 210: ch2co + o <=> co2 + ch2 */
    /*eqcon[209] *= 1;  */

    /*reaction 211: ch2co + h <=> ch3 + co */
    /*eqcon[210] *= 1;  */

    /*reaction 212: ch2co + h <=> hcco + h2 */
    /*eqcon[211] *= 1;  */

    /*reaction 213: ch2co + o <=> hcco + oh */
    /*eqcon[212] *= 1;  */

    /*reaction 214: ch2co + oh <=> hcco + h2o */
    /*eqcon[213] *= 1;  */

    /*reaction 215: ch2co + oh <=> ch2oh + co */
    /*eqcon[214] *= 1;  */

    /*reaction 216: ch2co (+M) <=> ch2 + co (+M) */
    eqcon[215] *= 1e-06; 

    /*reaction 217: c2h + h2 <=> c2h2 + h */
    /*eqcon[216] *= 1;  */

    /*reaction 218: c2h + o <=> ch + co */
    /*eqcon[217] *= 1;  */

    /*reaction 219: c2h + oh <=> hcco + h */
    /*eqcon[218] *= 1;  */

    /*reaction 220: c2h + o2 <=> co + co + h */
    eqcon[219] *= 1e-06; 

    /*reaction 221: hcco + c2h2 <=> h2ccch + co */
    /*eqcon[220] *= 1;  */

    /*reaction 222: hcco + h <=> ch2s + co */
    /*eqcon[221] *= 1;  */

    /*reaction 223: hcco + o <=> h + co + co */
    eqcon[222] *= 1e-06; 

    /*reaction 224: hcco + o <=> ch + co2 */
    /*eqcon[223] *= 1;  */

    /*reaction 225: hcco + o2 <=> hco + co + o */
    eqcon[224] *= 1e-06; 

    /*reaction 226: hcco + o2 <=> co2 + hco */
    /*eqcon[225] *= 1;  */

    /*reaction 227: hcco + ch <=> c2h2 + co */
    /*eqcon[226] *= 1;  */

    /*reaction 228: hcco + hcco <=> c2h2 + co + co */
    eqcon[227] *= 1e-06; 
}


/*compute the production rate for each species */
void productionRate(double * wdot, double * sc, double T)
{
    double qdot;

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[38];                /*Gibbs free energy */
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
    for (id = 0; id < 38; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*zero out wdot */
    for (id = 0; id < 38; ++id) {
        wdot[id] = 0.0;
    }

    /*reaction 1: oh + h2 <=> h + h2o */
    phi_f = sc[12]*sc[0];
    k_f = 1e-06 * 2.14e+08*exp(1.52*tc[0]-1735.5942803604782512/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[15];
    Kc = exp((g_RT[12] + g_RT[0]) - (g_RT[1] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 2: o + oh <=> o2 + h */
    phi_f = sc[11]*sc[12];
    k_f = 1e-06 * 2.02e+14*exp(-0.4*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[1];
    Kc = exp((g_RT[11] + g_RT[12]) - (g_RT[10] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 3: o + h2 <=> oh + h */
    phi_f = sc[11]*sc[0];
    k_f = 1e-06 * 50600*exp(2.67*tc[0]-3165.2328279116868543/tc[1]);
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

    /*reaction 4: h + o2 (+M) <=> ho2 (+M) */
    phi_f = sc[1]*sc[10];
    alpha = mixture + -1*sc[15] + -1*sc[0] + -1*sc[37] + 9*sc[2] + 2.8*sc[8] + 0.9*sc[9];
    k_f = 1e-06 * 4.52e+13;
    redP = 1e-12 * alpha / k_f * 1.05e+19*exp(-1.257*tc[0]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[13];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 5: h + o2 (+n2) <=> ho2 (+n2) */
    phi_f = sc[1]*sc[10];
    alpha = sc[37];
    k_f = 1e-06 * 4.52e+13;
    redP = 1e-12 * alpha / k_f * 2.03e+20*exp(-1.59*tc[0]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[13];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 6: h + o2 (+h2) <=> ho2 (+h2) */
    phi_f = sc[1]*sc[10];
    alpha = sc[0];
    k_f = 1e-06 * 4.52e+13;
    redP = 1e-12 * alpha / k_f * 1.52e+19*exp(-1.133*tc[0]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[13];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 7: h + o2 (+h2o) <=> ho2 (+h2o) */
    phi_f = sc[1]*sc[10];
    alpha = sc[15];
    k_f = 1e-06 * 4.52e+13;
    redP = 1e-12 * alpha / k_f * 2.1e+23*exp(-2.437*tc[0]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[13];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 8: oh + ho2 <=> h2o + o2 */
    phi_f = sc[12]*sc[13];
    k_f = 1e-06 * 2.13e+28*exp(-4.827*tc[0]-1761.2583303165188227/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[10];
    Kc = exp((g_RT[12] + g_RT[13]) - (g_RT[15] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[15] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 9: oh + ho2 <=> h2o + o2 */
    phi_f = sc[12]*sc[13];
    k_f = 1e-06 * 9.1e+14*exp(-5517.267523882946989/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[10];
    Kc = exp((g_RT[12] + g_RT[13]) - (g_RT[15] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[15] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 10: h + ho2 <=> oh + oh */
    phi_f = sc[1]*sc[13];
    k_f = 1e-06 * 1.5e+14*exp(-503.21666580471963925/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[12];
    Kc = exp((g_RT[1] + g_RT[13]) - (g_RT[12] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 11: h + ho2 <=> h2 + o2 */
    phi_f = sc[1]*sc[13];
    k_f = 1e-06 * 6.63e+13*exp(-1069.8386315008340262/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[10];
    Kc = exp((g_RT[1] + g_RT[13]) - (g_RT[0] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 12: h + ho2 <=> o + h2o */
    phi_f = sc[1]*sc[13];
    k_f = 1e-06 * 3.01e+13*exp(-866.03588184992258903/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[15];
    Kc = exp((g_RT[1] + g_RT[13]) - (g_RT[11] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 13: o + ho2 <=> o2 + oh */
    phi_f = sc[11]*sc[13];
    k_f = 1e-06 * 3.25e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[12];
    Kc = exp((g_RT[11] + g_RT[13]) - (g_RT[10] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 14: oh + oh <=> o + h2o */
    phi_f = sc[12]*sc[12];
    k_f = 1e-06 * 35700*exp(2.4*tc[0]+1062.7935981795678799/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[15];
    Kc = exp((g_RT[12] + g_RT[12]) - (g_RT[11] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 15: h + h + M <=> h2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[15] + -1*sc[0];
    k_f = 1e-12 * alpha * 1e+18*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[1]) - (g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 16: h + h + h2 <=> h2 + h2 */
    phi_f = sc[1]*sc[1]*sc[0];
    k_f = 1e-12 * 9.2e+16*exp(-0.6*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[0];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[1] + g_RT[0]) - (g_RT[0] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 17: h + h + h2o <=> h2 + h2o */
    phi_f = sc[1]*sc[1]*sc[15];
    k_f = 1e-12 * 6e+19*exp(-1.25*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[15];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[1] + g_RT[15]) - (g_RT[0] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 18: h + oh + M <=> h2o + M */
    phi_f = sc[1]*sc[12];
    alpha = mixture + 5.4*sc[15];
    k_f = 1e-12 * alpha * 2.21e+22*exp(-2*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[15];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[12]) - (g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 19: h + o + M <=> oh + M */
    phi_f = sc[1]*sc[11];
    alpha = mixture + 5.4*sc[15];
    k_f = 1e-12 * alpha * 4.71e+18*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[12];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[11]) - (g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 20: o + o + M <=> o2 + M */
    phi_f = sc[11]*sc[11];
    alpha = mixture;
    k_f = 1e-12 * alpha * 1.89e+13*exp(+899.75139845883882117/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10];
    Kc = 1.0 / (refC) * exp((g_RT[11] + g_RT[11]) - (g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 21: ho2 + ho2 <=> h2o2 + o2 */
    phi_f = sc[13]*sc[13];
    k_f = 1e-06 * 4.2e+14*exp(-6029.5420896721516328/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[10];
    Kc = exp((g_RT[13] + g_RT[13]) - (g_RT[14] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 22: ho2 + ho2 <=> h2o2 + o2 */
    phi_f = sc[13]*sc[13];
    k_f = 1e-06 * 1.3e+11*exp(+819.73994859588844974/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[10];
    Kc = exp((g_RT[13] + g_RT[13]) - (g_RT[14] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 23: oh + oh (+M) <=> h2o2 (+M) */
    phi_f = sc[12]*sc[12];
    alpha = mixture;
    k_f = 1e-06 * 1.24e+14*exp(-0.37*tc[0]);
    redP = 1e-12 * alpha / k_f * 3.04e+30*exp(-4.63*tc[0]-1031.0909482338706766/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.53*exp(T/-100))+ (0.47*exp(T/-2000))+ (exp(-1e+15/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[14];
    Kc = 1.0 / (refC) * exp((g_RT[12] + g_RT[12]) - (g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 24: h2o2 + h <=> ho2 + h2 */
    phi_f = sc[14]*sc[1];
    k_f = 1e-06 * 1.98e+06*exp(2*tc[0]-1225.3325812344924088/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[0];
    Kc = exp((g_RT[14] + g_RT[1]) - (g_RT[13] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[14] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 25: h2o2 + h <=> oh + h2o */
    phi_f = sc[14]*sc[1];
    k_f = 1e-06 * 3.07e+13*exp(-2122.0646796985029141/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[15];
    Kc = exp((g_RT[14] + g_RT[1]) - (g_RT[12] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[14] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 26: h2o2 + o <=> oh + ho2 */
    phi_f = sc[14]*sc[11];
    k_f = 1e-06 * 9.55e+06*exp(2*tc[0]-1997.7701632447369775/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[13];
    Kc = exp((g_RT[14] + g_RT[11]) - (g_RT[12] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[14] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 27: h2o2 + oh <=> h2o + ho2 */
    phi_f = sc[14]*sc[12];
    k_f = 1e-06 * 2.4*exp(4.042*tc[0]+1087.9544314698041489/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[13];
    Kc = exp((g_RT[14] + g_RT[12]) - (g_RT[15] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[14] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[15] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 28: ch3 + ch3 (+M) <=> c2h6 (+M) */
    phi_f = sc[3]*sc[3];
    alpha = mixture + 4*sc[15] + sc[0] + 2*sc[8] + sc[9];
    k_f = 1e-06 * 9.22e+16*exp(-1.174*tc[0]-320.04579945180171308/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.14e+36*exp(-5.246*tc[0]-857.98441519704704206/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.595*exp(T/-1120))+ (0.405*exp(T/-69.6))+ (exp(-1e+15/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[22];
    Kc = 1.0 / (refC) * exp((g_RT[3] + g_RT[3]) - (g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 29: ch3 + h (+M) <=> ch4 (+M) */
    phi_f = sc[3]*sc[1];
    alpha = mixture + 4*sc[15] + sc[0] + 2*sc[8] + sc[9];
    k_f = 1e-06 * 2.14e+15*exp(-0.4*tc[0]);
    redP = 1e-12 * alpha / k_f * 3.31e+30*exp(-4*tc[0]-1060.7807315163490784/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((1*exp(T/-1e-15))+ (0*exp(T/-1e-15))+ (exp(-40/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[2];
    Kc = 1.0 / (refC) * exp((g_RT[3] + g_RT[1]) - (g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 30: ch4 + h <=> ch3 + h2 */
    phi_f = sc[2]*sc[1];
    k_f = 1e-06 * 22000*exp(3*tc[0]-4403.145825791297284/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = exp((g_RT[2] + g_RT[1]) - (g_RT[3] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 31: ch4 + oh <=> ch3 + h2o */
    phi_f = sc[2]*sc[12];
    k_f = 1e-06 * 4.19e+06*exp(2*tc[0]-1281.6928478046211239/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[15];
    Kc = exp((g_RT[2] + g_RT[12]) - (g_RT[3] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 32: ch4 + o <=> ch3 + oh */
    phi_f = sc[2]*sc[11];
    k_f = 1e-06 * 6.92e+08*exp(1.56*tc[0]-4269.7934093530466271/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[12];
    Kc = exp((g_RT[2] + g_RT[11]) - (g_RT[3] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 33: ch4 + ho2 <=> ch3 + h2o2 */
    phi_f = sc[2]*sc[13];
    k_f = 1e-06 * 1.12e+13*exp(-12399.258645428293676/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[14];
    Kc = exp((g_RT[2] + g_RT[13]) - (g_RT[3] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 34: ch3 + ho2 <=> ch3o + oh */
    phi_f = sc[3]*sc[13];
    k_f = 1e-06 * 7e+12;
    q_f = phi_f * k_f;
    phi_r = sc[24]*sc[12];
    Kc = exp((g_RT[3] + g_RT[13]) - (g_RT[24] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[24] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 35: ch3 + ho2 <=> ch4 + o2 */
    phi_f = sc[3]*sc[13];
    k_f = 1e-06 * 3e+12;
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[10];
    Kc = exp((g_RT[3] + g_RT[13]) - (g_RT[2] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 36: ch3 + o <=> ch2o + h */
    phi_f = sc[3]*sc[11];
    k_f = 1e-06 * 8e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[1];
    Kc = exp((g_RT[3] + g_RT[11]) - (g_RT[6] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 37: ch3 + o2 <=> ch3o + o */
    phi_f = sc[3]*sc[10];
    k_f = 1e-06 * 1.45e+13*exp(-14698.455591490057486/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[24]*sc[11];
    Kc = exp((g_RT[3] + g_RT[10]) - (g_RT[24] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[24] += 1 * qdot;
    wdot[11] += 1 * qdot;

    /*reaction 38: ch3 + o2 <=> ch2o + oh */
    phi_f = sc[3]*sc[10];
    k_f = 1e-06 * 2.51e+11*exp(-7367.0919873810962599/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[12];
    Kc = exp((g_RT[3] + g_RT[10]) - (g_RT[6] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 39: ch3o + h <=> ch3 + oh */
    phi_f = sc[24]*sc[1];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[12];
    Kc = exp((g_RT[24] + g_RT[1]) - (g_RT[3] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[24] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 40: ch2oh + h <=> ch3 + oh */
    phi_f = sc[23]*sc[1];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[12];
    Kc = exp((g_RT[23] + g_RT[1]) - (g_RT[3] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[23] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 41: ch3 + oh <=> ch2s + h2o */
    phi_f = sc[3]*sc[12];
    k_f = 1e-06 * 2e+13*exp(-276.76916619259583285/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[28]*sc[15];
    Kc = exp((g_RT[3] + g_RT[12]) - (g_RT[28] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[28] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 42: ch3 + oh <=> ch2 + h2o */
    phi_f = sc[3]*sc[12];
    k_f = 1e-06 * 3e+06*exp(2*tc[0]-1258.0416645117993539/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[15];
    Kc = exp((g_RT[3] + g_RT[12]) - (g_RT[4] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 43: ch3 + h <=> ch2 + h2 */
    phi_f = sc[3]*sc[1];
    k_f = 1e-06 * 9e+13*exp(-7598.5716536512672974/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[0];
    Kc = exp((g_RT[3] + g_RT[1]) - (g_RT[4] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 44: ch3 + M <=> ch + h2 + M */
    phi_f = sc[3];
    alpha = mixture;
    k_f = 1e-06 * alpha * 6.9e+14*exp(-41499.775212249434844/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[0];
    Kc = refC * exp((g_RT[3]) - (g_RT[5] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 45: ch3 + M <=> ch2 + h + M */
    phi_f = sc[3];
    alpha = mixture;
    k_f = 1e-06 * alpha * 1.9e+16*exp(-45999.538637875230052/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[1];
    Kc = refC * exp((g_RT[3]) - (g_RT[4] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 46: ch2o + h (+M) <=> ch3o (+M) */
    phi_f = sc[6]*sc[1];
    alpha = mixture + 4*sc[15];
    k_f = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1308.3633310922712099/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.5e+30*exp(-4.8*tc[0]-2797.8846618742413739/tc[1]);
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
    phi_r = sc[24];
    Kc = 1.0 / (refC) * exp((g_RT[6] + g_RT[1]) - (g_RT[24]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[24] += 1 * qdot;

    /*reaction 47: ch2o + h (+M) <=> ch2oh (+M) */
    phi_f = sc[6]*sc[1];
    alpha = mixture + 4*sc[15];
    k_f = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1811.579996896990906/tc[1]);
    redP = 1e-12 * alpha / k_f * 9.1e+31*exp(-4.82*tc[0]-3286.0048277048194905/tc[1]);
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
    phi_r = sc[23];
    Kc = 1.0 / (refC) * exp((g_RT[6] + g_RT[1]) - (g_RT[23]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[23] += 1 * qdot;

    /*reaction 48: ch3o + ch3 <=> ch2o + ch4 */
    phi_f = sc[24]*sc[3];
    k_f = 1e-06 * 1.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[2];
    Kc = exp((g_RT[24] + g_RT[3]) - (g_RT[6] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[24] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 49: ch3o + h <=> ch2o + h2 */
    phi_f = sc[24]*sc[1];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = exp((g_RT[24] + g_RT[1]) - (g_RT[6] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[24] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 50: ch2oh + h <=> ch2o + h2 */
    phi_f = sc[23]*sc[1];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = exp((g_RT[23] + g_RT[1]) - (g_RT[6] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[23] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 51: ch3o + oh <=> ch2o + h2o */
    phi_f = sc[24]*sc[12];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[15];
    Kc = exp((g_RT[24] + g_RT[12]) - (g_RT[6] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[24] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 52: ch2oh + oh <=> ch2o + h2o */
    phi_f = sc[23]*sc[12];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[15];
    Kc = exp((g_RT[23] + g_RT[12]) - (g_RT[6] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[23] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 53: ch3o + o <=> ch2o + oh */
    phi_f = sc[24]*sc[11];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[12];
    Kc = exp((g_RT[24] + g_RT[11]) - (g_RT[6] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[24] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 54: ch2oh + o <=> ch2o + oh */
    phi_f = sc[23]*sc[11];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[12];
    Kc = exp((g_RT[23] + g_RT[11]) - (g_RT[6] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[23] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 55: ch3o + o2 <=> ch2o + ho2 */
    phi_f = sc[24]*sc[10];
    k_f = 1e-06 * 6.3e+10*exp(-1308.3633310922712099/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = exp((g_RT[24] + g_RT[10]) - (g_RT[6] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[24] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 56: ch3o + co <=> ch3 + co2 */
    phi_f = sc[24]*sc[9];
    k_f = 1e-06 * 468*exp(3.16*tc[0]-2707.3056620293918968/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[8];
    Kc = exp((g_RT[24] + g_RT[9]) - (g_RT[3] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[24] -= 1 * qdot;
    wdot[9] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 57: ch2oh + o2 <=> ch2o + ho2 */
    phi_f = sc[23]*sc[10];
    k_f = 1e-06 * 1.57e+15*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = exp((g_RT[23] + g_RT[10]) - (g_RT[6] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[23] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 58: ch2oh + o2 <=> ch2o + ho2 */
    phi_f = sc[23]*sc[10];
    k_f = 1e-06 * 7.23e+13*exp(-1800.0060135834821722/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = exp((g_RT[23] + g_RT[10]) - (g_RT[6] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[23] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 59: ch2 + h <=> ch + h2 */
    phi_f = sc[4]*sc[1];
    k_f = 1e-06 * 1e+18*exp(-1.56*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[0];
    Kc = exp((g_RT[4] + g_RT[1]) - (g_RT[5] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 60: ch2 + oh <=> ch + h2o */
    phi_f = sc[4]*sc[12];
    k_f = 1e-06 * 1.13e+07*exp(2*tc[0]-1509.6499974141590883/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[15];
    Kc = exp((g_RT[4] + g_RT[12]) - (g_RT[5] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 61: ch2 + oh <=> ch2o + h */
    phi_f = sc[4]*sc[12];
    k_f = 1e-06 * 2.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[1];
    Kc = exp((g_RT[4] + g_RT[12]) - (g_RT[6] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 62: ch2 + co2 <=> ch2o + co */
    phi_f = sc[4]*sc[8];
    k_f = 1e-06 * 1.1e+11*exp(-503.21666580471963925/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[9];
    Kc = exp((g_RT[4] + g_RT[8]) - (g_RT[6] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 63: ch2 + o <=> co + h + h */
    phi_f = sc[4]*sc[11];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[1]*sc[1];
    Kc = refC * exp((g_RT[4] + g_RT[11]) - (g_RT[9] + g_RT[1] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 64: ch2 + o <=> co + h2 */
    phi_f = sc[4]*sc[11];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = exp((g_RT[4] + g_RT[11]) - (g_RT[9] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 65: ch2 + o2 <=> ch2o + o */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 3.29e+21*exp(-3.3*tc[0]-1443.2253975279359111/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[11];
    Kc = exp((g_RT[4] + g_RT[10]) - (g_RT[6] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[11] += 1 * qdot;

    /*reaction 66: ch2 + o2 <=> co2 + h + h */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 3.29e+21*exp(-3.3*tc[0]-1443.2253975279359111/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[1]*sc[1];
    Kc = refC * exp((g_RT[4] + g_RT[10]) - (g_RT[8] + g_RT[1] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 67: ch2 + o2 <=> co2 + h2 */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 1.01e+21*exp(-3.3*tc[0]-758.85073203351726079/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = exp((g_RT[4] + g_RT[10]) - (g_RT[8] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 68: ch2 + o2 <=> co + h2o */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 7.28e+19*exp(-2.54*tc[0]-910.31894844073792683/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[15];
    Kc = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 69: ch2 + o2 <=> hco + oh */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 1.29e+20*exp(-3.3*tc[0]-142.9135330885404187/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[12];
    Kc = exp((g_RT[4] + g_RT[10]) - (g_RT[7] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 70: ch2 + ch3 <=> c2h4 + h */
    phi_f = sc[4]*sc[3];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[1];
    Kc = exp((g_RT[4] + g_RT[3]) - (g_RT[20] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[20] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 71: ch2 + ch2 <=> c2h2 + h + h */
    phi_f = sc[4]*sc[4];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[1]*sc[1];
    Kc = refC * exp((g_RT[4] + g_RT[4]) - (g_RT[18] + g_RT[1] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 72: ch2 + hcco <=> c2h3 + co */
    phi_f = sc[4]*sc[17];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[9];
    Kc = exp((g_RT[4] + g_RT[17]) - (g_RT[19] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[19] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 73: ch2 + c2h2 <=> h2ccch + h */
    phi_f = sc[4]*sc[18];
    k_f = 1e-06 * 1.2e+13*exp(-3321.2299943111497669/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[26]*sc[1];
    Kc = exp((g_RT[4] + g_RT[18]) - (g_RT[26] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[26] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 74: ch2s + M <=> ch2 + M */
    phi_f = sc[28];
    alpha = mixture + 11*sc[1] + 3*sc[18] + 2*sc[15];
    k_f = 1e-06 * alpha * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4];
    Kc = exp((g_RT[28]) - (g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[28] -= 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 75: ch2s + ch4 <=> ch3 + ch3 */
    phi_f = sc[28]*sc[2];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[3];
    Kc = exp((g_RT[28] + g_RT[2]) - (g_RT[3] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[28] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 76: ch2s + c2h6 <=> ch3 + c2h5 */
    phi_f = sc[28]*sc[22];
    k_f = 1e-06 * 1.2e+14;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[21];
    Kc = exp((g_RT[28] + g_RT[22]) - (g_RT[3] + g_RT[21]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[28] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[21] += 1 * qdot;

    /*reaction 77: ch2s + o2 <=> co + oh + h */
    phi_f = sc[28]*sc[10];
    k_f = 1e-06 * 7e+13;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[12]*sc[1];
    Kc = refC * exp((g_RT[28] + g_RT[10]) - (g_RT[9] + g_RT[12] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[28] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 78: ch2s + h2 <=> ch3 + h */
    phi_f = sc[28]*sc[0];
    k_f = 1e-06 * 7e+13;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[1];
    Kc = exp((g_RT[28] + g_RT[0]) - (g_RT[3] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[28] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 79: ch2s + c2h2 <=> h2ccch + h */
    phi_f = sc[28]*sc[18];
    k_f = 1e-06 * 1.5e+14;
    q_f = phi_f * k_f;
    phi_r = sc[26]*sc[1];
    Kc = exp((g_RT[28] + g_RT[18]) - (g_RT[26] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[28] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[26] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 80: ch2s + o <=> co + h + h */
    phi_f = sc[28]*sc[11];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[1]*sc[1];
    Kc = refC * exp((g_RT[28] + g_RT[11]) - (g_RT[9] + g_RT[1] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[28] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 81: ch2s + oh <=> ch2o + h */
    phi_f = sc[28]*sc[12];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[1];
    Kc = exp((g_RT[28] + g_RT[12]) - (g_RT[6] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[28] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 82: ch2s + h <=> ch + h2 */
    phi_f = sc[28]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[0];
    Kc = exp((g_RT[28] + g_RT[1]) - (g_RT[5] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[28] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 83: ch2s + co2 <=> ch2o + co */
    phi_f = sc[28]*sc[8];
    k_f = 1e-06 * 3e+12;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[9];
    Kc = exp((g_RT[28] + g_RT[8]) - (g_RT[6] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[28] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 84: ch2s + ch3 <=> c2h4 + h */
    phi_f = sc[28]*sc[3];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[1];
    Kc = exp((g_RT[28] + g_RT[3]) - (g_RT[20] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[28] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[20] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 85: ch2s + ch2co <=> c2h4 + co */
    phi_f = sc[28]*sc[29];
    k_f = 1e-06 * 1.6e+14;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[9];
    Kc = exp((g_RT[28] + g_RT[29]) - (g_RT[20] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[28] -= 1 * qdot;
    wdot[29] -= 1 * qdot;
    wdot[20] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 86: ch + o2 <=> hco + o */
    phi_f = sc[5]*sc[10];
    k_f = 1e-06 * 3.3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[11];
    Kc = exp((g_RT[5] + g_RT[10]) - (g_RT[7] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[11] += 1 * qdot;

    /*reaction 87: ch + o <=> co + h */
    phi_f = sc[5]*sc[11];
    k_f = 1e-06 * 5.7e+13;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[1];
    Kc = exp((g_RT[5] + g_RT[11]) - (g_RT[9] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 88: ch + oh <=> hco + h */
    phi_f = sc[5]*sc[12];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[1];
    Kc = exp((g_RT[5] + g_RT[12]) - (g_RT[7] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 89: ch + co2 <=> hco + co */
    phi_f = sc[5]*sc[8];
    k_f = 1e-06 * 3.4e+12*exp(-347.2194994052565562/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[9];
    Kc = exp((g_RT[5] + g_RT[8]) - (g_RT[7] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 90: ch + h2o <=> ch2o + h */
    phi_f = sc[5]*sc[15];
    k_f = 1e-06 * 1.17e+15*exp(-0.75*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[1];
    Kc = exp((g_RT[5] + g_RT[15]) - (g_RT[6] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 91: ch + ch2o <=> ch2co + h */
    phi_f = sc[5]*sc[6];
    k_f = 1e-06 * 9.46e+13*exp(+259.15658288943063781/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[29]*sc[1];
    Kc = exp((g_RT[5] + g_RT[6]) - (g_RT[29] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[29] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 92: ch + c2h2 <=> c3h2 + h */
    phi_f = sc[5]*sc[18];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[1];
    Kc = exp((g_RT[5] + g_RT[18]) - (g_RT[27] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[27] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 93: ch + ch2 <=> c2h2 + h */
    phi_f = sc[5]*sc[4];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[1];
    Kc = exp((g_RT[5] + g_RT[4]) - (g_RT[18] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 94: ch + ch3 <=> c2h3 + h */
    phi_f = sc[5]*sc[3];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[1];
    Kc = exp((g_RT[5] + g_RT[3]) - (g_RT[19] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[19] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 95: ch + ch4 <=> c2h4 + h */
    phi_f = sc[5]*sc[2];
    k_f = 1e-06 * 6e+13;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[1];
    Kc = exp((g_RT[5] + g_RT[2]) - (g_RT[20] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[5] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[20] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 96: ch2o + oh <=> hco + h2o */
    phi_f = sc[6]*sc[12];
    k_f = 1e-06 * 3.43e+09*exp(1.18*tc[0]+224.93784961470970529/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[15];
    Kc = exp((g_RT[6] + g_RT[12]) - (g_RT[7] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 97: ch2o + h <=> hco + h2 */
    phi_f = sc[6]*sc[1];
    k_f = 1e-06 * 2.19e+08*exp(1.77*tc[0]-1509.6499974141590883/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[0];
    Kc = exp((g_RT[6] + g_RT[1]) - (g_RT[7] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 98: ch2o + M <=> hco + h + M */
    phi_f = sc[6];
    alpha = mixture;
    k_f = 1e-06 * alpha * 3.31e+16*exp(-40760.549930182300159/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[1];
    Kc = refC * exp((g_RT[6]) - (g_RT[7] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 99: ch2o + o <=> hco + oh */
    phi_f = sc[6]*sc[11];
    k_f = 1e-06 * 1.8e+13*exp(-1549.9073306785367095/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[12];
    Kc = exp((g_RT[6] + g_RT[11]) - (g_RT[7] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 100: hco + o2 <=> co + ho2 */
    phi_f = sc[7]*sc[10];
    k_f = 1e-06 * 7.58e+12*exp(-206.31883297993508108/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[13];
    Kc = exp((g_RT[7] + g_RT[10]) - (g_RT[9] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 101: hco + M <=> h + co + M */
    phi_f = sc[7];
    alpha = mixture + 4*sc[15] + 0.87*sc[0] + 2*sc[8] + 0.87*sc[9] + 1.81*sc[2];
    k_f = 1e-06 * alpha * 1.86e+17*exp(-1*tc[0]-8554.6833186802341515/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[9];
    Kc = refC * exp((g_RT[7]) - (g_RT[1] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 102: hco + oh <=> h2o + co */
    phi_f = sc[7]*sc[12];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[9];
    Kc = exp((g_RT[7] + g_RT[12]) - (g_RT[15] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[15] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 103: hco + h <=> co + h2 */
    phi_f = sc[7]*sc[1];
    k_f = 1e-06 * 1.19e+13*exp(0.25*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = exp((g_RT[7] + g_RT[1]) - (g_RT[9] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 104: hco + o <=> co + oh */
    phi_f = sc[7]*sc[11];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[12];
    Kc = exp((g_RT[7] + g_RT[11]) - (g_RT[9] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 105: hco + o <=> co2 + h */
    phi_f = sc[7]*sc[11];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[1];
    Kc = exp((g_RT[7] + g_RT[11]) - (g_RT[8] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 106: co + oh <=> co2 + h */
    phi_f = sc[9]*sc[12];
    k_f = 1e-06 * 9420*exp(2.25*tc[0]+1183.0623813068959862/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[1];
    Kc = exp((g_RT[9] + g_RT[12]) - (g_RT[8] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 107: co + o + M <=> co2 + M */
    phi_f = sc[9]*sc[11];
    alpha = mixture;
    k_f = 1e-12 * alpha * 6.17e+14*exp(-1509.6499974141590883/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8];
    Kc = 1.0 / (refC) * exp((g_RT[9] + g_RT[11]) - (g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 108: co + o2 <=> co2 + o */
    phi_f = sc[9]*sc[10];
    k_f = 1e-06 * 2.53e+12*exp(-23997.396358895472076/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[11];
    Kc = exp((g_RT[9] + g_RT[10]) - (g_RT[8] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[11] += 1 * qdot;

    /*reaction 109: co + ho2 <=> co2 + oh */
    phi_f = sc[9]*sc[13];
    k_f = 1e-06 * 5.8e+13*exp(-11540.771013565441535/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[12];
    Kc = exp((g_RT[9] + g_RT[13]) - (g_RT[8] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 110: c2h5oh (+M) <=> ch3 + ch2oh (+M) */
    phi_f = sc[33];
    alpha = mixture + 4*sc[15] + sc[0] + 2*sc[8] + sc[9];
    k_f = 1 * 5.94e+23*exp(-1.68*tc[0]-45874.740904755664815/tc[1]);
    redP = 1e-6 * alpha / k_f * 2.88e+85*exp(-18.9*tc[0]-55310.556605259960634/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-200))+ (0.5*exp(T/-890))+ (exp(-4600/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[23];
    Kc = refC * exp((g_RT[33]) - (g_RT[3] + g_RT[23]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[23] += 1 * qdot;

    /*reaction 111: c2h5oh (+M) <=> c2h5 + oh (+M) */
    phi_f = sc[33];
    alpha = mixture + 4*sc[15] + sc[0] + 2*sc[8] + sc[9];
    k_f = 1 * 1.25e+23*exp(-1.54*tc[0]-48311.316000582111883/tc[1]);
    redP = 1e-6 * alpha / k_f * 3.252e+85*exp(-18.81*tc[0]-57834.691400936433638/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-300))+ (0.5*exp(T/-900))+ (exp(-5000/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[12];
    Kc = refC * exp((g_RT[33]) - (g_RT[21] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[21] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 112: c2h5oh (+M) <=> c2h4 + h2o (+M) */
    phi_f = sc[33];
    alpha = mixture + 4*sc[15];
    k_f = 1 * 2.79e+13*exp(0.09*tc[0]-33280.737409660941921/tc[1]);
    redP = 1e-6 * alpha / k_f * 2.57e+83*exp(-18.85*tc[0]-43504.087192149629118/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.3*exp(T/-350))+ (0.7*exp(T/-800))+ (exp(-3800/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[15];
    Kc = refC * exp((g_RT[33]) - (g_RT[20] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[20] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 113: c2h5oh (+M) <=> ch3hco + h2 (+M) */
    phi_f = sc[33];
    alpha = mixture + 4*sc[15];
    k_f = 1 * 7.24e+11*exp(0.095*tc[0]-45796.239104890126328/tc[1]);
    redP = 1e-6 * alpha / k_f * 4.46e+87*exp(-19.42*tc[0]-58164.801533704332542/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.1*exp(T/-900))+ (0.9*exp(T/-1100))+ (exp(-3500/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[0];
    Kc = refC * exp((g_RT[33]) - (g_RT[32] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 114: c2h5oh + oh <=> c2h4oh + h2o */
    phi_f = sc[33]*sc[12];
    k_f = 1e-06 * 1.74e+11*exp(0.27*tc[0]-301.92999948283181766/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[34]*sc[15];
    Kc = exp((g_RT[33] + g_RT[12]) - (g_RT[34] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[34] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 115: c2h5oh + oh <=> ch3choh + h2o */
    phi_f = sc[33]*sc[12];
    k_f = 1e-06 * 4.64e+11*exp(0.15*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[15];
    Kc = exp((g_RT[33] + g_RT[12]) - (g_RT[35] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 116: c2h5oh + oh <=> ch3ch2o + h2o */
    phi_f = sc[33]*sc[12];
    k_f = 1e-06 * 7.46e+11*exp(0.3*tc[0]-822.25603192491200844/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[36]*sc[15];
    Kc = exp((g_RT[33] + g_RT[12]) - (g_RT[36] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[36] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 117: c2h5oh + h <=> c2h4oh + h2 */
    phi_f = sc[33]*sc[1];
    k_f = 1e-06 * 1.23e+07*exp(1.8*tc[0]-2565.3985622724612767/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[34]*sc[0];
    Kc = exp((g_RT[33] + g_RT[1]) - (g_RT[34] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[34] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 118: c2h5oh + h <=> ch3choh + h2 */
    phi_f = sc[33]*sc[1];
    k_f = 1e-06 * 2.58e+07*exp(1.65*tc[0]-1422.593514229942457/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[0];
    Kc = exp((g_RT[33] + g_RT[1]) - (g_RT[35] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 119: c2h5oh + h <=> ch3ch2o + h2 */
    phi_f = sc[33]*sc[1];
    k_f = 1e-06 * 1.5e+07*exp(1.6*tc[0]-1528.7722307147384981/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[36]*sc[0];
    Kc = exp((g_RT[33] + g_RT[1]) - (g_RT[36] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[36] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 120: c2h5oh + o <=> c2h4oh + oh */
    phi_f = sc[33]*sc[11];
    k_f = 1e-06 * 9.41e+07*exp(1.7*tc[0]-2747.0597786279645334/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[34]*sc[12];
    Kc = exp((g_RT[33] + g_RT[11]) - (g_RT[34] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[34] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 121: c2h5oh + o <=> ch3choh + oh */
    phi_f = sc[33]*sc[11];
    k_f = 1e-06 * 1.88e+07*exp(1.85*tc[0]-917.86719842780871659/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[12];
    Kc = exp((g_RT[33] + g_RT[11]) - (g_RT[35] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 122: c2h5oh + o <=> ch3ch2o + oh */
    phi_f = sc[33]*sc[11];
    k_f = 1e-06 * 1.58e+07*exp(2*tc[0]-2238.3077294993931901/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[36]*sc[12];
    Kc = exp((g_RT[33] + g_RT[11]) - (g_RT[36] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[36] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 123: c2h5oh + ch3 <=> c2h4oh + ch4 */
    phi_f = sc[33]*sc[3];
    k_f = 1e-06 * 219*exp(3.18*tc[0]-4841.9507583730128317/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[34]*sc[2];
    Kc = exp((g_RT[33] + g_RT[3]) - (g_RT[34] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[34] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 124: c2h5oh + ch3 <=> ch3choh + ch4 */
    phi_f = sc[33]*sc[3];
    k_f = 1e-06 * 728*exp(2.99*tc[0]-3999.5660598159120127/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[2];
    Kc = exp((g_RT[33] + g_RT[3]) - (g_RT[35] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 125: c2h5oh + ch3 <=> ch3ch2o + ch4 */
    phi_f = sc[33]*sc[3];
    k_f = 1e-06 * 145*exp(2.99*tc[0]-3849.1042767403009748/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[36]*sc[2];
    Kc = exp((g_RT[33] + g_RT[3]) - (g_RT[36] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[36] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 126: c2h5oh + ho2 <=> ch3choh + h2o2 */
    phi_f = sc[33]*sc[13];
    k_f = 1e-06 * 8200*exp(2.55*tc[0]-5409.579157400737131/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[14];
    Kc = exp((g_RT[33] + g_RT[13]) - (g_RT[35] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[35] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 127: c2h5oh + ho2 <=> c2h4oh + h2o2 */
    phi_f = sc[33]*sc[13];
    k_f = 1e-06 * 12300*exp(2.55*tc[0]-7925.6624864243349293/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[34]*sc[14];
    Kc = exp((g_RT[33] + g_RT[13]) - (g_RT[34] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[34] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 128: c2h5oh + ho2 <=> ch3ch2o + h2o2 */
    phi_f = sc[33]*sc[13];
    k_f = 1e-06 * 2.5e+12*exp(-12077.199979313272706/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[36]*sc[14];
    Kc = exp((g_RT[33] + g_RT[13]) - (g_RT[36] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[33] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[36] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 129: ch3ch2o + M <=> ch3hco + h + M */
    phi_f = sc[36];
    alpha = mixture;
    k_f = 1e-06 * alpha * 1.16e+35*exp(-5.89*tc[0]-12718.298011548486102/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[1];
    Kc = refC * exp((g_RT[36]) - (g_RT[32] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 130: ch3ch2o + M <=> ch3 + ch2o + M */
    phi_f = sc[36];
    alpha = mixture;
    k_f = 1e-06 * alpha * 1.35e+38*exp(-6.96*tc[0]-11976.55664615232854/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[6];
    Kc = refC * exp((g_RT[36]) - (g_RT[3] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 131: ch3ch2o + o2 <=> ch3hco + ho2 */
    phi_f = sc[36]*sc[10];
    k_f = 1e-06 * 4e+10*exp(-553.53833238519166571/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[13];
    Kc = exp((g_RT[36] + g_RT[10]) - (g_RT[32] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 132: ch3ch2o + co <=> c2h5 + co2 */
    phi_f = sc[36]*sc[9];
    k_f = 1e-06 * 468*exp(3.16*tc[0]-2707.3056620293918968/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[8];
    Kc = exp((g_RT[36] + g_RT[9]) - (g_RT[21] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[9] -= 1 * qdot;
    wdot[21] += 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 133: ch3ch2o + h <=> ch3 + ch2oh */
    phi_f = sc[36]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[23];
    Kc = exp((g_RT[36] + g_RT[1]) - (g_RT[3] + g_RT[23]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[23] += 1 * qdot;

    /*reaction 134: ch3ch2o + h <=> c2h4 + h2o */
    phi_f = sc[36]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[15];
    Kc = exp((g_RT[36] + g_RT[1]) - (g_RT[20] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[20] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 135: ch3ch2o + oh <=> ch3hco + h2o */
    phi_f = sc[36]*sc[12];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[15];
    Kc = exp((g_RT[36] + g_RT[12]) - (g_RT[32] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[36] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 136: ch3choh + o2 <=> ch3hco + ho2 */
    phi_f = sc[35]*sc[10];
    k_f = 1e-06 * 4.82e+14*exp(-2524.638012342278671/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[13];
    Kc = exp((g_RT[35] + g_RT[10]) - (g_RT[32] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[35] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 137: ch3choh + o2 <=> ch3hco + ho2 */
    phi_f = sc[35]*sc[10];
    k_f = 1e-06 * 8.43e+15*exp(-1.2*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[13];
    Kc = exp((g_RT[35] + g_RT[10]) - (g_RT[32] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[35] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 138: ch3choh + o <=> ch3hco + oh */
    phi_f = sc[35]*sc[11];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[12];
    Kc = exp((g_RT[35] + g_RT[11]) - (g_RT[32] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[35] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 139: ch3choh + h <=> c2h4 + h2o */
    phi_f = sc[35]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[15];
    Kc = exp((g_RT[35] + g_RT[1]) - (g_RT[20] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[35] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[20] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 140: ch3choh + h <=> ch3 + ch2oh */
    phi_f = sc[35]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[23];
    Kc = exp((g_RT[35] + g_RT[1]) - (g_RT[3] + g_RT[23]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[35] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[23] += 1 * qdot;

    /*reaction 141: ch3choh + ho2 <=> ch3hco + oh + oh */
    phi_f = sc[35]*sc[13];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[12]*sc[12];
    Kc = refC * exp((g_RT[35] + g_RT[13]) - (g_RT[32] + g_RT[12] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[35] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 142: ch3choh + oh <=> ch3hco + h2o */
    phi_f = sc[35]*sc[12];
    k_f = 1e-06 * 5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[15];
    Kc = exp((g_RT[35] + g_RT[12]) - (g_RT[32] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[35] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 143: ch3choh + M <=> ch3hco + h + M */
    phi_f = sc[35];
    alpha = mixture;
    k_f = 1e-06 * alpha * 1e+14*exp(-12580.416645117993539/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[1];
    Kc = refC * exp((g_RT[35]) - (g_RT[32] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[35] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 144: ch3hco + oh <=> ch3co + h2o */
    phi_f = sc[32]*sc[12];
    k_f = 1e-06 * 9.24e+06*exp(1.5*tc[0]+484.0944325041403431/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[15];
    Kc = exp((g_RT[32] + g_RT[12]) - (g_RT[31] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[32] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[31] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 145: ch3hco + oh <=> ch2hco + h2o */
    phi_f = sc[32]*sc[12];
    k_f = 1e-06 * 172000*exp(2.4*tc[0]-410.12158263084654664/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[15];
    Kc = exp((g_RT[32] + g_RT[12]) - (g_RT[30] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[32] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[30] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 146: ch3hco + o <=> ch3co + oh */
    phi_f = sc[32]*sc[11];
    k_f = 1e-06 * 1.77e+18*exp(-1.9*tc[0]-1497.0695807690410675/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[12];
    Kc = exp((g_RT[32] + g_RT[11]) - (g_RT[31] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[32] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[31] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 147: ch3hco + o <=> ch2hco + oh */
    phi_f = sc[32]*sc[11];
    k_f = 1e-06 * 3.72e+13*exp(-0.2*tc[0]-1789.4384636015831802/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[12];
    Kc = exp((g_RT[32] + g_RT[11]) - (g_RT[30] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[32] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[30] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 148: ch3hco + h <=> ch3co + h2 */
    phi_f = sc[32]*sc[1];
    k_f = 1e-06 * 4.66e+13*exp(-0.35*tc[0]-1503.6113974245024565/tc[1]);
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

    /*reaction 149: ch3hco + h <=> ch2hco + h2 */
    phi_f = sc[32]*sc[1];
    k_f = 1e-06 * 1.85e+12*exp(0.4*tc[0]-2696.7381120474929048/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[0];
    Kc = exp((g_RT[32] + g_RT[1]) - (g_RT[30] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[32] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[30] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 150: ch3hco + ch3 <=> ch3co + ch4 */
    phi_f = sc[32]*sc[3];
    k_f = 1e-06 * 3.9e-07*exp(5.8*tc[0]-1107.0766647703833314/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[2];
    Kc = exp((g_RT[32] + g_RT[3]) - (g_RT[31] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[32] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[31] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 151: ch3hco + ch3 <=> ch2hco + ch4 */
    phi_f = sc[32]*sc[3];
    k_f = 1e-06 * 24.5*exp(3.15*tc[0]-2881.9218450636299167/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[2];
    Kc = exp((g_RT[32] + g_RT[3]) - (g_RT[30] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[32] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[30] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 152: ch3hco + ho2 <=> ch3co + h2o2 */
    phi_f = sc[32]*sc[13];
    k_f = 1e-06 * 2.4e+19*exp(-2.2*tc[0]-7060.1298212402170975/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[14];
    Kc = exp((g_RT[32] + g_RT[13]) - (g_RT[31] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[32] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[31] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 153: ch3hco + ho2 <=> ch2hco + h2o2 */
    phi_f = sc[32]*sc[13];
    k_f = 1e-06 * 2.32e+11*exp(0.4*tc[0]-7479.8125205213536901/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[14];
    Kc = exp((g_RT[32] + g_RT[13]) - (g_RT[30] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[32] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[30] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 154: ch3hco + o2 <=> ch3co + ho2 */
    phi_f = sc[32]*sc[10];
    k_f = 1e-06 * 1e+14*exp(-21235.743296959171857/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[13];
    Kc = exp((g_RT[32] + g_RT[10]) - (g_RT[31] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[32] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[31] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 155: c2h6 + ch3 <=> c2h5 + ch4 */
    phi_f = sc[22]*sc[3];
    k_f = 1e-06 * 0.55*exp(4*tc[0]-4176.6983261791738187/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[2];
    Kc = exp((g_RT[22] + g_RT[3]) - (g_RT[21] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[22] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[21] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 156: c2h6 + h <=> c2h5 + h2 */
    phi_f = sc[22]*sc[1];
    k_f = 1e-06 * 540*exp(3.5*tc[0]-2621.7588288425899918/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[0];
    Kc = exp((g_RT[22] + g_RT[1]) - (g_RT[21] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[22] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[21] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 157: c2h6 + o <=> c2h5 + oh */
    phi_f = sc[22]*sc[11];
    k_f = 1e-06 * 3e+07*exp(2*tc[0]-2573.9532455911412399/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[12];
    Kc = exp((g_RT[22] + g_RT[11]) - (g_RT[21] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[22] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[21] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 158: c2h6 + oh <=> c2h5 + h2o */
    phi_f = sc[22]*sc[12];
    k_f = 1e-06 * 7.23e+06*exp(2*tc[0]-434.77919925527783107/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[15];
    Kc = exp((g_RT[22] + g_RT[12]) - (g_RT[21] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[22] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[21] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 159: c2h5 + h <=> c2h4 + h2 */
    phi_f = sc[21]*sc[1];
    k_f = 1e-06 * 1.25e+14*exp(-4025.733326437757114/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[0];
    Kc = exp((g_RT[21] + g_RT[1]) - (g_RT[20] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[20] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 160: c2h5 + h <=> ch3 + ch3 */
    phi_f = sc[21]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[3];
    Kc = exp((g_RT[21] + g_RT[1]) - (g_RT[3] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 161: c2h5 + h <=> c2h6 */
    phi_f = sc[21]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[22];
    Kc = 1.0 / (refC) * exp((g_RT[21] + g_RT[1]) - (g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 162: c2h5 + oh <=> c2h4 + h2o */
    phi_f = sc[21]*sc[12];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[15];
    Kc = exp((g_RT[21] + g_RT[12]) - (g_RT[20] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[20] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 163: c2h5 + o <=> ch3 + ch2o */
    phi_f = sc[21]*sc[11];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[6];
    Kc = exp((g_RT[21] + g_RT[11]) - (g_RT[3] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 164: c2h5 + ho2 <=> c2h6 + o2 */
    phi_f = sc[21]*sc[13];
    k_f = 1e-06 * 3e+12;
    q_f = phi_f * k_f;
    phi_r = sc[22]*sc[10];
    Kc = exp((g_RT[21] + g_RT[13]) - (g_RT[22] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[22] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 165: c2h5 + ho2 <=> ch3ch2o + oh */
    phi_f = sc[21]*sc[13];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[36]*sc[12];
    Kc = exp((g_RT[21] + g_RT[13]) - (g_RT[36] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[36] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 166: c2h5 + o2 <=> c2h4 + ho2 */
    phi_f = sc[21]*sc[10];
    k_f = 1e-06 * 2.89e+28*exp(-5.4*tc[0]-3816.8984101287987869/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[13];
    Kc = exp((g_RT[21] + g_RT[10]) - (g_RT[20] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[20] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 167: c2h5 + o2 <=> ch3hco + oh */
    phi_f = sc[21]*sc[10];
    k_f = 1e-06 * 4.9e+11*exp(-0.48*tc[0]-4205.3816761300422513/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[12];
    Kc = exp((g_RT[21] + g_RT[10]) - (g_RT[32] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 168: c2h4 + oh <=> c2h4oh */
    phi_f = sc[20]*sc[12];
    k_f = 1e-06 * 1.29e+12*exp(+411.12801596245600422/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[34];
    Kc = 1.0 / (refC) * exp((g_RT[20] + g_RT[12]) - (g_RT[34]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[20] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[34] += 1 * qdot;

    /*reaction 169: c2h4 + oh <=> c2h3 + h2o */
    phi_f = sc[20]*sc[12];
    k_f = 1e-06 * 2.02e+13*exp(-2987.0941282168159887/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[15];
    Kc = exp((g_RT[20] + g_RT[12]) - (g_RT[19] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[20] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[19] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 170: c2h4 + o <=> ch3 + hco */
    phi_f = sc[20]*sc[11];
    k_f = 1e-06 * 1.02e+07*exp(1.88*tc[0]-90.07578317904483356/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = exp((g_RT[20] + g_RT[11]) - (g_RT[3] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[20] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 171: c2h4 + o <=> ch2hco + h */
    phi_f = sc[20]*sc[11];
    k_f = 1e-06 * 3.39e+06*exp(1.88*tc[0]-90.07578317904483356/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[1];
    Kc = exp((g_RT[20] + g_RT[11]) - (g_RT[30] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[20] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[30] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 172: c2h4 + ch3 <=> c2h3 + ch4 */
    phi_f = sc[20]*sc[3];
    k_f = 1e-06 * 6.62*exp(3.7*tc[0]-4780.5583251448369992/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[2];
    Kc = exp((g_RT[20] + g_RT[3]) - (g_RT[19] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[20] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[19] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 173: c2h4 + h <=> c2h3 + h2 */
    phi_f = sc[20]*sc[1];
    k_f = 1e-06 * 3.36e-07*exp(6*tc[0]-851.44259854158576672/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[0];
    Kc = exp((g_RT[20] + g_RT[1]) - (g_RT[19] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[20] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[19] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 174: c2h4 + h (+M) <=> c2h5 (+M) */
    phi_f = sc[20]*sc[1];
    alpha = mixture + 4*sc[15] + sc[0] + 2*sc[8] + sc[9];
    k_f = 1e-06 * 1.08e+12*exp(0.454*tc[0]-916.86076509619920216/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.112e+34*exp(-5*tc[0]-2238.3077294993931901/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0*exp(T/-1e-15))+ (1*exp(T/-95))+ (exp(-200/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[21];
    Kc = 1.0 / (refC) * exp((g_RT[20] + g_RT[1]) - (g_RT[21]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[20] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[21] += 1 * qdot;

    /*reaction 175: c2h4 (+M) <=> c2h2 + h2 (+M) */
    phi_f = sc[20];
    alpha = mixture;
    k_f = 1 * 1.8e+14*exp(-43779.849925010610605/tc[1]);
    redP = 1e-6 * alpha / k_f * 1.5e+15*exp(-27899.841602211075951/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[0];
    Kc = refC * exp((g_RT[20]) - (g_RT[18] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[20] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 176: c2h3 + h (+M) <=> c2h4 (+M) */
    phi_f = sc[19]*sc[1];
    alpha = mixture + 4*sc[15];
    k_f = 1e-06 * 6.1e+12*exp(0.27*tc[0]-140.90066642532150354/tc[1]);
    redP = 1e-12 * alpha / k_f * 9.8e+29*exp(-3.86*tc[0]-1670.6793304716693456/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.218*exp(T/-208))+ (0.782*exp(T/-2663))+ (exp(-6095/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[20];
    Kc = 1.0 / (refC) * exp((g_RT[19] + g_RT[1]) - (g_RT[20]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[20] += 1 * qdot;

    /*reaction 177: c2h3 + h <=> c2h2 + h2 */
    phi_f = sc[19]*sc[1];
    k_f = 1e-06 * 9e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[0];
    Kc = exp((g_RT[19] + g_RT[1]) - (g_RT[18] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 178: c2h3 + o <=> ch2co + h */
    phi_f = sc[19]*sc[11];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[29]*sc[1];
    Kc = exp((g_RT[19] + g_RT[11]) - (g_RT[29] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[29] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 179: c2h3 + o2 <=> ch2o + hco */
    phi_f = sc[19]*sc[10];
    k_f = 1e-06 * 1.7e+29*exp(-5.312*tc[0]-3270.9083277306781383/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[7];
    Kc = exp((g_RT[19] + g_RT[10]) - (g_RT[6] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 180: c2h3 + o2 <=> ch2hco + o */
    phi_f = sc[19]*sc[10];
    k_f = 1e-06 * 5.5e+14*exp(-0.611*tc[0]-2646.9196621328255787/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[11];
    Kc = exp((g_RT[19] + g_RT[10]) - (g_RT[30] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[30] += 1 * qdot;
    wdot[11] += 1 * qdot;

    /*reaction 181: c2h3 + o2 <=> c2h2 + ho2 */
    phi_f = sc[19]*sc[10];
    k_f = 1e-06 * 2.12e-06*exp(6*tc[0]-4772.5068584919617933/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[13];
    Kc = exp((g_RT[19] + g_RT[10]) - (g_RT[18] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 182: c2h3 + oh <=> c2h2 + h2o */
    phi_f = sc[19]*sc[12];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[15];
    Kc = exp((g_RT[19] + g_RT[12]) - (g_RT[18] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 183: c2h3 + c2h <=> c2h2 + c2h2 */
    phi_f = sc[19]*sc[16];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[18];
    Kc = exp((g_RT[19] + g_RT[16]) - (g_RT[18] + g_RT[18]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[16] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[18] += 1 * qdot;

    /*reaction 184: c2h3 + ch <=> ch2 + c2h2 */
    phi_f = sc[19]*sc[5];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[18];
    Kc = exp((g_RT[19] + g_RT[5]) - (g_RT[4] + g_RT[18]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[18] += 1 * qdot;

    /*reaction 185: c2h3 + ch3 <=> c2h2 + ch4 */
    phi_f = sc[19]*sc[3];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[2];
    Kc = exp((g_RT[19] + g_RT[3]) - (g_RT[18] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 186: c2h2 + oh <=> c2h + h2o */
    phi_f = sc[18]*sc[12];
    k_f = 1e-06 * 3.37e+07*exp(2*tc[0]-7045.0333212660752906/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[15];
    Kc = exp((g_RT[18] + g_RT[12]) - (g_RT[16] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 187: c2h2 + oh <=> hccoh + h */
    phi_f = sc[18]*sc[12];
    k_f = 1e-06 * 504000*exp(2.3*tc[0]-6793.4249883637157836/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[1];
    Kc = exp((g_RT[18] + g_RT[12]) - (g_RT[25] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[25] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 188: c2h2 + oh <=> ch2co + h */
    phi_f = sc[18]*sc[12];
    k_f = 1e-06 * 0.000218*exp(4.5*tc[0]+503.21666580471963925/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[29]*sc[1];
    Kc = exp((g_RT[18] + g_RT[12]) - (g_RT[29] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[29] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 189: c2h2 + oh <=> ch2co + h */
    phi_f = sc[18]*sc[12];
    k_f = 1e-06 * 2e+11;
    q_f = phi_f * k_f;
    phi_r = sc[29]*sc[1];
    Kc = exp((g_RT[18] + g_RT[12]) - (g_RT[29] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[29] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 190: c2h2 + oh <=> ch3 + co */
    phi_f = sc[18]*sc[12];
    k_f = 1e-06 * 0.000483*exp(4*tc[0]+1006.4333316094392785/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[9];
    Kc = exp((g_RT[18] + g_RT[12]) - (g_RT[3] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 191: hccoh + h <=> ch2co + h */
    phi_f = sc[25]*sc[1];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[29]*sc[1];
    Kc = exp((g_RT[25] + g_RT[1]) - (g_RT[29] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[25] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[29] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 192: c2h2 + o <=> ch2 + co */
    phi_f = sc[18]*sc[11];
    k_f = 1e-06 * 6.12e+06*exp(2*tc[0]-956.11166502896742259/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = exp((g_RT[18] + g_RT[11]) - (g_RT[4] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 193: c2h2 + o <=> hcco + h */
    phi_f = sc[18]*sc[11];
    k_f = 1e-06 * 1.43e+07*exp(2*tc[0]-956.11166502896742259/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[1];
    Kc = exp((g_RT[18] + g_RT[11]) - (g_RT[17] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 194: c2h2 + o <=> c2h + oh */
    phi_f = sc[18]*sc[11];
    k_f = 1e-06 * 3.16e+15*exp(-0.6*tc[0]-7548.2499870707952141/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[12];
    Kc = exp((g_RT[18] + g_RT[11]) - (g_RT[16] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 195: c2h2 + ch3 <=> c2h + ch4 */
    phi_f = sc[18]*sc[3];
    k_f = 1e-06 * 1.81e+11*exp(-8700.1129350977989816/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[2];
    Kc = exp((g_RT[18] + g_RT[3]) - (g_RT[16] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 196: c2h2 + o2 <=> hcco + oh */
    phi_f = sc[18]*sc[10];
    k_f = 1e-06 * 4e+07*exp(1.5*tc[0]-15146.821640722062511/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[12];
    Kc = exp((g_RT[18] + g_RT[10]) - (g_RT[17] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 197: c2h2 + M <=> c2h + h + M */
    phi_f = sc[18];
    alpha = mixture;
    k_f = 1e-06 * alpha * 4.2e+16*exp(-53844.183241105005436/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[1];
    Kc = refC * exp((g_RT[18]) - (g_RT[16] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 198: c2h2 + h (+M) <=> c2h3 (+M) */
    phi_f = sc[18]*sc[1];
    alpha = mixture + 4*sc[15] + sc[0] + 2*sc[8] + sc[9];
    k_f = 1e-06 * 3.11e+11*exp(0.58*tc[0]-1302.8279477684191079/tc[1]);
    redP = 1e-12 * alpha / k_f * 2.25e+40*exp(-7.269*tc[0]-3309.6560109976417152/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0*exp(T/-1e-15))+ (1*exp(T/-675))+ (exp(-1e+15/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[19];
    Kc = 1.0 / (refC) * exp((g_RT[18] + g_RT[1]) - (g_RT[19]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[19] += 1 * qdot;

    /*reaction 199: ch2hco + h <=> ch3 + hco */
    phi_f = sc[30]*sc[1];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = exp((g_RT[30] + g_RT[1]) - (g_RT[3] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[30] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 200: ch2hco + h <=> ch2co + h2 */
    phi_f = sc[30]*sc[1];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[29]*sc[0];
    Kc = exp((g_RT[30] + g_RT[1]) - (g_RT[29] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[30] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[29] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 201: ch2hco + o <=> ch2o + hco */
    phi_f = sc[30]*sc[11];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[7];
    Kc = exp((g_RT[30] + g_RT[11]) - (g_RT[6] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[30] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 202: ch2hco + oh <=> ch2co + h2o */
    phi_f = sc[30]*sc[12];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[29]*sc[15];
    Kc = exp((g_RT[30] + g_RT[12]) - (g_RT[29] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[30] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[29] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 203: ch2hco + o2 <=> ch2o + co + oh */
    phi_f = sc[30]*sc[10];
    k_f = 1e-06 * 3e+10;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[9]*sc[12];
    Kc = refC * exp((g_RT[30] + g_RT[10]) - (g_RT[6] + g_RT[9] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[30] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 204: ch2hco + ch3 <=> c2h5 + co + h */
    phi_f = sc[30]*sc[3];
    k_f = 1e-06 * 4.9e+14*exp(-0.5*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[9]*sc[1];
    Kc = refC * exp((g_RT[30] + g_RT[3]) - (g_RT[21] + g_RT[9] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[30] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[21] += 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 205: ch2hco + ho2 <=> ch2o + hco + oh */
    phi_f = sc[30]*sc[13];
    k_f = 1e-06 * 7e+12;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[7]*sc[12];
    Kc = refC * exp((g_RT[30] + g_RT[13]) - (g_RT[6] + g_RT[7] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[30] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 206: ch2hco + ho2 <=> ch3hco + o2 */
    phi_f = sc[30]*sc[13];
    k_f = 1e-06 * 3e+12;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[10];
    Kc = exp((g_RT[30] + g_RT[13]) - (g_RT[32] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[30] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[32] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 207: ch2hco <=> ch3 + co */
    phi_f = sc[30];
    k_f = 1 * 1.17e+43*exp(-9.83*tc[0]-22018.748428951312235/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[9];
    Kc = refC * exp((g_RT[30]) - (g_RT[3] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[30] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 208: ch2hco <=> ch2co + h */
    phi_f = sc[30];
    k_f = 1 * 1.81e+43*exp(-9.61*tc[0]-23081.542027130883071/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[29]*sc[1];
    Kc = refC * exp((g_RT[30]) - (g_RT[29] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[30] -= 1 * qdot;
    wdot[29] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 209: ch3co (+M) <=> ch3 + co (+M) */
    phi_f = sc[31];
    alpha = mixture;
    k_f = 1 * 3e+12*exp(-8414.7890855865243793/tc[1]);
    redP = 1e-6 * alpha / k_f * 1.2e+15*exp(-6299.2662225434805805/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[9];
    Kc = refC * exp((g_RT[31]) - (g_RT[3] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 210: ch2co + o <=> co2 + ch2 */
    phi_f = sc[29]*sc[11];
    k_f = 1e-06 * 1.75e+12*exp(-679.34249883637164658/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[4];
    Kc = exp((g_RT[29] + g_RT[11]) - (g_RT[8] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[29] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 211: ch2co + h <=> ch3 + co */
    phi_f = sc[29]*sc[1];
    k_f = 1e-06 * 27100*exp(2.75*tc[0]-359.29669938456987666/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[9];
    Kc = exp((g_RT[29] + g_RT[1]) - (g_RT[3] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[29] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 212: ch2co + h <=> hcco + h2 */
    phi_f = sc[29]*sc[1];
    k_f = 1e-06 * 2e+14*exp(-4025.733326437757114/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[0];
    Kc = exp((g_RT[29] + g_RT[1]) - (g_RT[17] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[29] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 213: ch2co + o <=> hcco + oh */
    phi_f = sc[29]*sc[11];
    k_f = 1e-06 * 1e+13*exp(-4025.733326437757114/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[12];
    Kc = exp((g_RT[29] + g_RT[11]) - (g_RT[17] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[29] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 214: ch2co + oh <=> hcco + h2o */
    phi_f = sc[29]*sc[12];
    k_f = 1e-06 * 1e+13*exp(-1006.4333316094392785/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[15];
    Kc = exp((g_RT[29] + g_RT[12]) - (g_RT[17] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[29] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 215: ch2co + oh <=> ch2oh + co */
    phi_f = sc[29]*sc[12];
    k_f = 1e-06 * 3.73e+12*exp(+509.75848246018102827/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[9];
    Kc = exp((g_RT[29] + g_RT[12]) - (g_RT[23] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[29] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[23] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 216: ch2co (+M) <=> ch2 + co (+M) */
    phi_f = sc[29];
    alpha = mixture;
    k_f = 1 * 3e+14*exp(-35718.318938819007599/tc[1]);
    redP = 1e-6 * alpha / k_f * 3.6e+15*exp(-29825.651782245739014/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = refC * exp((g_RT[29]) - (g_RT[4] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[29] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 217: c2h + h2 <=> c2h2 + h */
    phi_f = sc[16]*sc[0];
    k_f = 1e-06 * 409000*exp(2.39*tc[0]-434.93016425501917865/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[1];
    Kc = exp((g_RT[16] + g_RT[0]) - (g_RT[18] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 218: c2h + o <=> ch + co */
    phi_f = sc[16]*sc[11];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[9];
    Kc = exp((g_RT[16] + g_RT[11]) - (g_RT[5] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 219: c2h + oh <=> hcco + h */
    phi_f = sc[16]*sc[12];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[1];
    Kc = exp((g_RT[16] + g_RT[12]) - (g_RT[17] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 220: c2h + o2 <=> co + co + h */
    phi_f = sc[16]*sc[10];
    k_f = 1e-06 * 9.04e+12*exp(+229.97001627275690794/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[9]*sc[1];
    Kc = refC * exp((g_RT[16] + g_RT[10]) - (g_RT[9] + g_RT[9] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 221: hcco + c2h2 <=> h2ccch + co */
    phi_f = sc[17]*sc[18];
    k_f = 1e-06 * 1e+11*exp(-1509.6499974141590883/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[26]*sc[9];
    Kc = exp((g_RT[17] + g_RT[18]) - (g_RT[26] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[17] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[26] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 222: hcco + h <=> ch2s + co */
    phi_f = sc[17]*sc[1];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[28]*sc[9];
    Kc = exp((g_RT[17] + g_RT[1]) - (g_RT[28] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[17] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[28] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 223: hcco + o <=> h + co + co */
    phi_f = sc[17]*sc[11];
    k_f = 1e-06 * 8e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[9]*sc[9];
    Kc = refC * exp((g_RT[17] + g_RT[11]) - (g_RT[1] + g_RT[9] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[17] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 224: hcco + o <=> ch + co2 */
    phi_f = sc[17]*sc[11];
    k_f = 1e-06 * 2.95e+13*exp(-560.08014904065305473/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[8];
    Kc = exp((g_RT[17] + g_RT[11]) - (g_RT[5] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[17] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[8] += 1 * qdot;

    /*reaction 225: hcco + o2 <=> hco + co + o */
    phi_f = sc[17]*sc[10];
    k_f = 1e-06 * 2.5e+08*exp(1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[9]*sc[11];
    Kc = refC * exp((g_RT[17] + g_RT[10]) - (g_RT[7] + g_RT[9] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[17] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[11] += 1 * qdot;

    /*reaction 226: hcco + o2 <=> co2 + hco */
    phi_f = sc[17]*sc[10];
    k_f = 1e-06 * 2.4e+11*exp(+429.7470325972306/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[7];
    Kc = exp((g_RT[17] + g_RT[10]) - (g_RT[8] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[17] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 227: hcco + ch <=> c2h2 + co */
    phi_f = sc[17]*sc[5];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[9];
    Kc = exp((g_RT[17] + g_RT[5]) - (g_RT[18] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[17] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 228: hcco + hcco <=> c2h2 + co + co */
    phi_f = sc[17]*sc[17];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[9]*sc[9];
    Kc = refC * exp((g_RT[17] + g_RT[17]) - (g_RT[18] + g_RT[9] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[17] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[9] += 1 * qdot;

    return;
}


/*compute the progress rate for each reaction */
void progressRate(double * qdot, double * sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[38];                /*Gibbs free energy */
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
    for (id = 0; id < 38; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*reaction 1: oh + h2 <=> h + h2o */
    phi_f = sc[12]*sc[0];
    k_f = 1e-06 * 2.14e+08*exp(1.52*tc[0]-1735.5942803604782512/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[15];
    Kc = exp((g_RT[12] + g_RT[0]) - (g_RT[1] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[0] = q_f - q_r;

    /*reaction 2: o + oh <=> o2 + h */
    phi_f = sc[11]*sc[12];
    k_f = 1e-06 * 2.02e+14*exp(-0.4*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[1];
    Kc = exp((g_RT[11] + g_RT[12]) - (g_RT[10] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[1] = q_f - q_r;

    /*reaction 3: o + h2 <=> oh + h */
    phi_f = sc[11]*sc[0];
    k_f = 1e-06 * 50600*exp(2.67*tc[0]-3165.2328279116868543/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[1];
    Kc = exp((g_RT[11] + g_RT[0]) - (g_RT[12] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[2] = q_f - q_r;

    /*reaction 4: h + o2 (+M) <=> ho2 (+M) */
    phi_f = sc[1]*sc[10];
    alpha = mixture + -1*sc[15] + -1*sc[0] + -1*sc[37] + 9*sc[2] + 2.8*sc[8] + 0.9*sc[9];
    k_f = 1e-06 * 4.52e+13;
    redP = 1e-12 * alpha / k_f * 1.05e+19*exp(-1.257*tc[0]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[13];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[3] = q_f - q_r;

    /*reaction 5: h + o2 (+n2) <=> ho2 (+n2) */
    phi_f = sc[1]*sc[10];
    alpha = sc[37];
    k_f = 1e-06 * 4.52e+13;
    redP = 1e-12 * alpha / k_f * 2.03e+20*exp(-1.59*tc[0]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[13];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[4] = q_f - q_r;

    /*reaction 6: h + o2 (+h2) <=> ho2 (+h2) */
    phi_f = sc[1]*sc[10];
    alpha = sc[0];
    k_f = 1e-06 * 4.52e+13;
    redP = 1e-12 * alpha / k_f * 1.52e+19*exp(-1.133*tc[0]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[13];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[5] = q_f - q_r;

    /*reaction 7: h + o2 (+h2o) <=> ho2 (+h2o) */
    phi_f = sc[1]*sc[10];
    alpha = sc[15];
    k_f = 1e-06 * 4.52e+13;
    redP = 1e-12 * alpha / k_f * 2.1e+23*exp(-2.437*tc[0]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[13];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[6] = q_f - q_r;

    /*reaction 8: oh + ho2 <=> h2o + o2 */
    phi_f = sc[12]*sc[13];
    k_f = 1e-06 * 2.13e+28*exp(-4.827*tc[0]-1761.2583303165188227/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[10];
    Kc = exp((g_RT[12] + g_RT[13]) - (g_RT[15] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[7] = q_f - q_r;

    /*reaction 9: oh + ho2 <=> h2o + o2 */
    phi_f = sc[12]*sc[13];
    k_f = 1e-06 * 9.1e+14*exp(-5517.267523882946989/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[10];
    Kc = exp((g_RT[12] + g_RT[13]) - (g_RT[15] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[8] = q_f - q_r;

    /*reaction 10: h + ho2 <=> oh + oh */
    phi_f = sc[1]*sc[13];
    k_f = 1e-06 * 1.5e+14*exp(-503.21666580471963925/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[12];
    Kc = exp((g_RT[1] + g_RT[13]) - (g_RT[12] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[9] = q_f - q_r;

    /*reaction 11: h + ho2 <=> h2 + o2 */
    phi_f = sc[1]*sc[13];
    k_f = 1e-06 * 6.63e+13*exp(-1069.8386315008340262/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[10];
    Kc = exp((g_RT[1] + g_RT[13]) - (g_RT[0] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[10] = q_f - q_r;

    /*reaction 12: h + ho2 <=> o + h2o */
    phi_f = sc[1]*sc[13];
    k_f = 1e-06 * 3.01e+13*exp(-866.03588184992258903/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[15];
    Kc = exp((g_RT[1] + g_RT[13]) - (g_RT[11] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[11] = q_f - q_r;

    /*reaction 13: o + ho2 <=> o2 + oh */
    phi_f = sc[11]*sc[13];
    k_f = 1e-06 * 3.25e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[12];
    Kc = exp((g_RT[11] + g_RT[13]) - (g_RT[10] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[12] = q_f - q_r;

    /*reaction 14: oh + oh <=> o + h2o */
    phi_f = sc[12]*sc[12];
    k_f = 1e-06 * 35700*exp(2.4*tc[0]+1062.7935981795678799/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[15];
    Kc = exp((g_RT[12] + g_RT[12]) - (g_RT[11] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[13] = q_f - q_r;

    /*reaction 15: h + h + M <=> h2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[15] + -1*sc[0];
    k_f = 1e-12 * alpha * 1e+18*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[1]) - (g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[14] = q_f - q_r;

    /*reaction 16: h + h + h2 <=> h2 + h2 */
    phi_f = sc[1]*sc[1]*sc[0];
    k_f = 1e-12 * 9.2e+16*exp(-0.6*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[0];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[1] + g_RT[0]) - (g_RT[0] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[15] = q_f - q_r;

    /*reaction 17: h + h + h2o <=> h2 + h2o */
    phi_f = sc[1]*sc[1]*sc[15];
    k_f = 1e-12 * 6e+19*exp(-1.25*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[15];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[1] + g_RT[15]) - (g_RT[0] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[16] = q_f - q_r;

    /*reaction 18: h + oh + M <=> h2o + M */
    phi_f = sc[1]*sc[12];
    alpha = mixture + 5.4*sc[15];
    k_f = 1e-12 * alpha * 2.21e+22*exp(-2*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[15];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[12]) - (g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[17] = q_f - q_r;

    /*reaction 19: h + o + M <=> oh + M */
    phi_f = sc[1]*sc[11];
    alpha = mixture + 5.4*sc[15];
    k_f = 1e-12 * alpha * 4.71e+18*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[12];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[11]) - (g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[18] = q_f - q_r;

    /*reaction 20: o + o + M <=> o2 + M */
    phi_f = sc[11]*sc[11];
    alpha = mixture;
    k_f = 1e-12 * alpha * 1.89e+13*exp(+899.75139845883882117/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10];
    Kc = 1.0 / (refC) * exp((g_RT[11] + g_RT[11]) - (g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[19] = q_f - q_r;

    /*reaction 21: ho2 + ho2 <=> h2o2 + o2 */
    phi_f = sc[13]*sc[13];
    k_f = 1e-06 * 4.2e+14*exp(-6029.5420896721516328/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[10];
    Kc = exp((g_RT[13] + g_RT[13]) - (g_RT[14] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[20] = q_f - q_r;

    /*reaction 22: ho2 + ho2 <=> h2o2 + o2 */
    phi_f = sc[13]*sc[13];
    k_f = 1e-06 * 1.3e+11*exp(+819.73994859588844974/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[10];
    Kc = exp((g_RT[13] + g_RT[13]) - (g_RT[14] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[21] = q_f - q_r;

    /*reaction 23: oh + oh (+M) <=> h2o2 (+M) */
    phi_f = sc[12]*sc[12];
    alpha = mixture;
    k_f = 1e-06 * 1.24e+14*exp(-0.37*tc[0]);
    redP = 1e-12 * alpha / k_f * 3.04e+30*exp(-4.63*tc[0]-1031.0909482338706766/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.53*exp(T/-100))+ (0.47*exp(T/-2000))+ (exp(-1e+15/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[14];
    Kc = 1.0 / (refC) * exp((g_RT[12] + g_RT[12]) - (g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[22] = q_f - q_r;

    /*reaction 24: h2o2 + h <=> ho2 + h2 */
    phi_f = sc[14]*sc[1];
    k_f = 1e-06 * 1.98e+06*exp(2*tc[0]-1225.3325812344924088/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[0];
    Kc = exp((g_RT[14] + g_RT[1]) - (g_RT[13] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[23] = q_f - q_r;

    /*reaction 25: h2o2 + h <=> oh + h2o */
    phi_f = sc[14]*sc[1];
    k_f = 1e-06 * 3.07e+13*exp(-2122.0646796985029141/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[15];
    Kc = exp((g_RT[14] + g_RT[1]) - (g_RT[12] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[24] = q_f - q_r;

    /*reaction 26: h2o2 + o <=> oh + ho2 */
    phi_f = sc[14]*sc[11];
    k_f = 1e-06 * 9.55e+06*exp(2*tc[0]-1997.7701632447369775/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[13];
    Kc = exp((g_RT[14] + g_RT[11]) - (g_RT[12] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[25] = q_f - q_r;

    /*reaction 27: h2o2 + oh <=> h2o + ho2 */
    phi_f = sc[14]*sc[12];
    k_f = 1e-06 * 2.4*exp(4.042*tc[0]+1087.9544314698041489/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[13];
    Kc = exp((g_RT[14] + g_RT[12]) - (g_RT[15] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[26] = q_f - q_r;

    /*reaction 28: ch3 + ch3 (+M) <=> c2h6 (+M) */
    phi_f = sc[3]*sc[3];
    alpha = mixture + 4*sc[15] + sc[0] + 2*sc[8] + sc[9];
    k_f = 1e-06 * 9.22e+16*exp(-1.174*tc[0]-320.04579945180171308/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.14e+36*exp(-5.246*tc[0]-857.98441519704704206/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.595*exp(T/-1120))+ (0.405*exp(T/-69.6))+ (exp(-1e+15/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[22];
    Kc = 1.0 / (refC) * exp((g_RT[3] + g_RT[3]) - (g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[27] = q_f - q_r;

    /*reaction 29: ch3 + h (+M) <=> ch4 (+M) */
    phi_f = sc[3]*sc[1];
    alpha = mixture + 4*sc[15] + sc[0] + 2*sc[8] + sc[9];
    k_f = 1e-06 * 2.14e+15*exp(-0.4*tc[0]);
    redP = 1e-12 * alpha / k_f * 3.31e+30*exp(-4*tc[0]-1060.7807315163490784/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((1*exp(T/-1e-15))+ (0*exp(T/-1e-15))+ (exp(-40/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[2];
    Kc = 1.0 / (refC) * exp((g_RT[3] + g_RT[1]) - (g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[28] = q_f - q_r;

    /*reaction 30: ch4 + h <=> ch3 + h2 */
    phi_f = sc[2]*sc[1];
    k_f = 1e-06 * 22000*exp(3*tc[0]-4403.145825791297284/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = exp((g_RT[2] + g_RT[1]) - (g_RT[3] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[29] = q_f - q_r;

    /*reaction 31: ch4 + oh <=> ch3 + h2o */
    phi_f = sc[2]*sc[12];
    k_f = 1e-06 * 4.19e+06*exp(2*tc[0]-1281.6928478046211239/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[15];
    Kc = exp((g_RT[2] + g_RT[12]) - (g_RT[3] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[30] = q_f - q_r;

    /*reaction 32: ch4 + o <=> ch3 + oh */
    phi_f = sc[2]*sc[11];
    k_f = 1e-06 * 6.92e+08*exp(1.56*tc[0]-4269.7934093530466271/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[12];
    Kc = exp((g_RT[2] + g_RT[11]) - (g_RT[3] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[31] = q_f - q_r;

    /*reaction 33: ch4 + ho2 <=> ch3 + h2o2 */
    phi_f = sc[2]*sc[13];
    k_f = 1e-06 * 1.12e+13*exp(-12399.258645428293676/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[14];
    Kc = exp((g_RT[2] + g_RT[13]) - (g_RT[3] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[32] = q_f - q_r;

    /*reaction 34: ch3 + ho2 <=> ch3o + oh */
    phi_f = sc[3]*sc[13];
    k_f = 1e-06 * 7e+12;
    q_f = phi_f * k_f;
    phi_r = sc[24]*sc[12];
    Kc = exp((g_RT[3] + g_RT[13]) - (g_RT[24] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[33] = q_f - q_r;

    /*reaction 35: ch3 + ho2 <=> ch4 + o2 */
    phi_f = sc[3]*sc[13];
    k_f = 1e-06 * 3e+12;
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[10];
    Kc = exp((g_RT[3] + g_RT[13]) - (g_RT[2] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[34] = q_f - q_r;

    /*reaction 36: ch3 + o <=> ch2o + h */
    phi_f = sc[3]*sc[11];
    k_f = 1e-06 * 8e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[1];
    Kc = exp((g_RT[3] + g_RT[11]) - (g_RT[6] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[35] = q_f - q_r;

    /*reaction 37: ch3 + o2 <=> ch3o + o */
    phi_f = sc[3]*sc[10];
    k_f = 1e-06 * 1.45e+13*exp(-14698.455591490057486/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[24]*sc[11];
    Kc = exp((g_RT[3] + g_RT[10]) - (g_RT[24] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[36] = q_f - q_r;

    /*reaction 38: ch3 + o2 <=> ch2o + oh */
    phi_f = sc[3]*sc[10];
    k_f = 1e-06 * 2.51e+11*exp(-7367.0919873810962599/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[12];
    Kc = exp((g_RT[3] + g_RT[10]) - (g_RT[6] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[37] = q_f - q_r;

    /*reaction 39: ch3o + h <=> ch3 + oh */
    phi_f = sc[24]*sc[1];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[12];
    Kc = exp((g_RT[24] + g_RT[1]) - (g_RT[3] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[38] = q_f - q_r;

    /*reaction 40: ch2oh + h <=> ch3 + oh */
    phi_f = sc[23]*sc[1];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[12];
    Kc = exp((g_RT[23] + g_RT[1]) - (g_RT[3] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[39] = q_f - q_r;

    /*reaction 41: ch3 + oh <=> ch2s + h2o */
    phi_f = sc[3]*sc[12];
    k_f = 1e-06 * 2e+13*exp(-276.76916619259583285/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[28]*sc[15];
    Kc = exp((g_RT[3] + g_RT[12]) - (g_RT[28] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[40] = q_f - q_r;

    /*reaction 42: ch3 + oh <=> ch2 + h2o */
    phi_f = sc[3]*sc[12];
    k_f = 1e-06 * 3e+06*exp(2*tc[0]-1258.0416645117993539/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[15];
    Kc = exp((g_RT[3] + g_RT[12]) - (g_RT[4] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[41] = q_f - q_r;

    /*reaction 43: ch3 + h <=> ch2 + h2 */
    phi_f = sc[3]*sc[1];
    k_f = 1e-06 * 9e+13*exp(-7598.5716536512672974/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[0];
    Kc = exp((g_RT[3] + g_RT[1]) - (g_RT[4] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[42] = q_f - q_r;

    /*reaction 44: ch3 + M <=> ch + h2 + M */
    phi_f = sc[3];
    alpha = mixture;
    k_f = 1e-06 * alpha * 6.9e+14*exp(-41499.775212249434844/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[0];
    Kc = refC * exp((g_RT[3]) - (g_RT[5] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[43] = q_f - q_r;

    /*reaction 45: ch3 + M <=> ch2 + h + M */
    phi_f = sc[3];
    alpha = mixture;
    k_f = 1e-06 * alpha * 1.9e+16*exp(-45999.538637875230052/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[1];
    Kc = refC * exp((g_RT[3]) - (g_RT[4] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[44] = q_f - q_r;

    /*reaction 46: ch2o + h (+M) <=> ch3o (+M) */
    phi_f = sc[6]*sc[1];
    alpha = mixture + 4*sc[15];
    k_f = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1308.3633310922712099/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.5e+30*exp(-4.8*tc[0]-2797.8846618742413739/tc[1]);
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
    phi_r = sc[24];
    Kc = 1.0 / (refC) * exp((g_RT[6] + g_RT[1]) - (g_RT[24]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[45] = q_f - q_r;

    /*reaction 47: ch2o + h (+M) <=> ch2oh (+M) */
    phi_f = sc[6]*sc[1];
    alpha = mixture + 4*sc[15];
    k_f = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1811.579996896990906/tc[1]);
    redP = 1e-12 * alpha / k_f * 9.1e+31*exp(-4.82*tc[0]-3286.0048277048194905/tc[1]);
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
    phi_r = sc[23];
    Kc = 1.0 / (refC) * exp((g_RT[6] + g_RT[1]) - (g_RT[23]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[46] = q_f - q_r;

    /*reaction 48: ch3o + ch3 <=> ch2o + ch4 */
    phi_f = sc[24]*sc[3];
    k_f = 1e-06 * 1.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[2];
    Kc = exp((g_RT[24] + g_RT[3]) - (g_RT[6] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[47] = q_f - q_r;

    /*reaction 49: ch3o + h <=> ch2o + h2 */
    phi_f = sc[24]*sc[1];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = exp((g_RT[24] + g_RT[1]) - (g_RT[6] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[48] = q_f - q_r;

    /*reaction 50: ch2oh + h <=> ch2o + h2 */
    phi_f = sc[23]*sc[1];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = exp((g_RT[23] + g_RT[1]) - (g_RT[6] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[49] = q_f - q_r;

    /*reaction 51: ch3o + oh <=> ch2o + h2o */
    phi_f = sc[24]*sc[12];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[15];
    Kc = exp((g_RT[24] + g_RT[12]) - (g_RT[6] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[50] = q_f - q_r;

    /*reaction 52: ch2oh + oh <=> ch2o + h2o */
    phi_f = sc[23]*sc[12];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[15];
    Kc = exp((g_RT[23] + g_RT[12]) - (g_RT[6] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[51] = q_f - q_r;

    /*reaction 53: ch3o + o <=> ch2o + oh */
    phi_f = sc[24]*sc[11];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[12];
    Kc = exp((g_RT[24] + g_RT[11]) - (g_RT[6] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[52] = q_f - q_r;

    /*reaction 54: ch2oh + o <=> ch2o + oh */
    phi_f = sc[23]*sc[11];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[12];
    Kc = exp((g_RT[23] + g_RT[11]) - (g_RT[6] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[53] = q_f - q_r;

    /*reaction 55: ch3o + o2 <=> ch2o + ho2 */
    phi_f = sc[24]*sc[10];
    k_f = 1e-06 * 6.3e+10*exp(-1308.3633310922712099/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = exp((g_RT[24] + g_RT[10]) - (g_RT[6] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[54] = q_f - q_r;

    /*reaction 56: ch3o + co <=> ch3 + co2 */
    phi_f = sc[24]*sc[9];
    k_f = 1e-06 * 468*exp(3.16*tc[0]-2707.3056620293918968/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[8];
    Kc = exp((g_RT[24] + g_RT[9]) - (g_RT[3] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[55] = q_f - q_r;

    /*reaction 57: ch2oh + o2 <=> ch2o + ho2 */
    phi_f = sc[23]*sc[10];
    k_f = 1e-06 * 1.57e+15*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = exp((g_RT[23] + g_RT[10]) - (g_RT[6] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[56] = q_f - q_r;

    /*reaction 58: ch2oh + o2 <=> ch2o + ho2 */
    phi_f = sc[23]*sc[10];
    k_f = 1e-06 * 7.23e+13*exp(-1800.0060135834821722/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = exp((g_RT[23] + g_RT[10]) - (g_RT[6] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[57] = q_f - q_r;

    /*reaction 59: ch2 + h <=> ch + h2 */
    phi_f = sc[4]*sc[1];
    k_f = 1e-06 * 1e+18*exp(-1.56*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[0];
    Kc = exp((g_RT[4] + g_RT[1]) - (g_RT[5] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[58] = q_f - q_r;

    /*reaction 60: ch2 + oh <=> ch + h2o */
    phi_f = sc[4]*sc[12];
    k_f = 1e-06 * 1.13e+07*exp(2*tc[0]-1509.6499974141590883/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[15];
    Kc = exp((g_RT[4] + g_RT[12]) - (g_RT[5] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[59] = q_f - q_r;

    /*reaction 61: ch2 + oh <=> ch2o + h */
    phi_f = sc[4]*sc[12];
    k_f = 1e-06 * 2.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[1];
    Kc = exp((g_RT[4] + g_RT[12]) - (g_RT[6] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[60] = q_f - q_r;

    /*reaction 62: ch2 + co2 <=> ch2o + co */
    phi_f = sc[4]*sc[8];
    k_f = 1e-06 * 1.1e+11*exp(-503.21666580471963925/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[9];
    Kc = exp((g_RT[4] + g_RT[8]) - (g_RT[6] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[61] = q_f - q_r;

    /*reaction 63: ch2 + o <=> co + h + h */
    phi_f = sc[4]*sc[11];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[1]*sc[1];
    Kc = refC * exp((g_RT[4] + g_RT[11]) - (g_RT[9] + g_RT[1] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[62] = q_f - q_r;

    /*reaction 64: ch2 + o <=> co + h2 */
    phi_f = sc[4]*sc[11];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = exp((g_RT[4] + g_RT[11]) - (g_RT[9] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[63] = q_f - q_r;

    /*reaction 65: ch2 + o2 <=> ch2o + o */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 3.29e+21*exp(-3.3*tc[0]-1443.2253975279359111/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[11];
    Kc = exp((g_RT[4] + g_RT[10]) - (g_RT[6] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[64] = q_f - q_r;

    /*reaction 66: ch2 + o2 <=> co2 + h + h */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 3.29e+21*exp(-3.3*tc[0]-1443.2253975279359111/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[1]*sc[1];
    Kc = refC * exp((g_RT[4] + g_RT[10]) - (g_RT[8] + g_RT[1] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[65] = q_f - q_r;

    /*reaction 67: ch2 + o2 <=> co2 + h2 */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 1.01e+21*exp(-3.3*tc[0]-758.85073203351726079/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = exp((g_RT[4] + g_RT[10]) - (g_RT[8] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[66] = q_f - q_r;

    /*reaction 68: ch2 + o2 <=> co + h2o */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 7.28e+19*exp(-2.54*tc[0]-910.31894844073792683/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[15];
    Kc = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[67] = q_f - q_r;

    /*reaction 69: ch2 + o2 <=> hco + oh */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 1.29e+20*exp(-3.3*tc[0]-142.9135330885404187/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[12];
    Kc = exp((g_RT[4] + g_RT[10]) - (g_RT[7] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[68] = q_f - q_r;

    /*reaction 70: ch2 + ch3 <=> c2h4 + h */
    phi_f = sc[4]*sc[3];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[1];
    Kc = exp((g_RT[4] + g_RT[3]) - (g_RT[20] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[69] = q_f - q_r;

    /*reaction 71: ch2 + ch2 <=> c2h2 + h + h */
    phi_f = sc[4]*sc[4];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[1]*sc[1];
    Kc = refC * exp((g_RT[4] + g_RT[4]) - (g_RT[18] + g_RT[1] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[70] = q_f - q_r;

    /*reaction 72: ch2 + hcco <=> c2h3 + co */
    phi_f = sc[4]*sc[17];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[9];
    Kc = exp((g_RT[4] + g_RT[17]) - (g_RT[19] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[71] = q_f - q_r;

    /*reaction 73: ch2 + c2h2 <=> h2ccch + h */
    phi_f = sc[4]*sc[18];
    k_f = 1e-06 * 1.2e+13*exp(-3321.2299943111497669/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[26]*sc[1];
    Kc = exp((g_RT[4] + g_RT[18]) - (g_RT[26] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[72] = q_f - q_r;

    /*reaction 74: ch2s + M <=> ch2 + M */
    phi_f = sc[28];
    alpha = mixture + 11*sc[1] + 3*sc[18] + 2*sc[15];
    k_f = 1e-06 * alpha * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4];
    Kc = exp((g_RT[28]) - (g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[73] = q_f - q_r;

    /*reaction 75: ch2s + ch4 <=> ch3 + ch3 */
    phi_f = sc[28]*sc[2];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[3];
    Kc = exp((g_RT[28] + g_RT[2]) - (g_RT[3] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[74] = q_f - q_r;

    /*reaction 76: ch2s + c2h6 <=> ch3 + c2h5 */
    phi_f = sc[28]*sc[22];
    k_f = 1e-06 * 1.2e+14;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[21];
    Kc = exp((g_RT[28] + g_RT[22]) - (g_RT[3] + g_RT[21]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[75] = q_f - q_r;

    /*reaction 77: ch2s + o2 <=> co + oh + h */
    phi_f = sc[28]*sc[10];
    k_f = 1e-06 * 7e+13;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[12]*sc[1];
    Kc = refC * exp((g_RT[28] + g_RT[10]) - (g_RT[9] + g_RT[12] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[76] = q_f - q_r;

    /*reaction 78: ch2s + h2 <=> ch3 + h */
    phi_f = sc[28]*sc[0];
    k_f = 1e-06 * 7e+13;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[1];
    Kc = exp((g_RT[28] + g_RT[0]) - (g_RT[3] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[77] = q_f - q_r;

    /*reaction 79: ch2s + c2h2 <=> h2ccch + h */
    phi_f = sc[28]*sc[18];
    k_f = 1e-06 * 1.5e+14;
    q_f = phi_f * k_f;
    phi_r = sc[26]*sc[1];
    Kc = exp((g_RT[28] + g_RT[18]) - (g_RT[26] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[78] = q_f - q_r;

    /*reaction 80: ch2s + o <=> co + h + h */
    phi_f = sc[28]*sc[11];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[1]*sc[1];
    Kc = refC * exp((g_RT[28] + g_RT[11]) - (g_RT[9] + g_RT[1] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[79] = q_f - q_r;

    /*reaction 81: ch2s + oh <=> ch2o + h */
    phi_f = sc[28]*sc[12];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[1];
    Kc = exp((g_RT[28] + g_RT[12]) - (g_RT[6] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[80] = q_f - q_r;

    /*reaction 82: ch2s + h <=> ch + h2 */
    phi_f = sc[28]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[0];
    Kc = exp((g_RT[28] + g_RT[1]) - (g_RT[5] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[81] = q_f - q_r;

    /*reaction 83: ch2s + co2 <=> ch2o + co */
    phi_f = sc[28]*sc[8];
    k_f = 1e-06 * 3e+12;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[9];
    Kc = exp((g_RT[28] + g_RT[8]) - (g_RT[6] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[82] = q_f - q_r;

    /*reaction 84: ch2s + ch3 <=> c2h4 + h */
    phi_f = sc[28]*sc[3];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[1];
    Kc = exp((g_RT[28] + g_RT[3]) - (g_RT[20] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[83] = q_f - q_r;

    /*reaction 85: ch2s + ch2co <=> c2h4 + co */
    phi_f = sc[28]*sc[29];
    k_f = 1e-06 * 1.6e+14;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[9];
    Kc = exp((g_RT[28] + g_RT[29]) - (g_RT[20] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[84] = q_f - q_r;

    /*reaction 86: ch + o2 <=> hco + o */
    phi_f = sc[5]*sc[10];
    k_f = 1e-06 * 3.3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[11];
    Kc = exp((g_RT[5] + g_RT[10]) - (g_RT[7] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[85] = q_f - q_r;

    /*reaction 87: ch + o <=> co + h */
    phi_f = sc[5]*sc[11];
    k_f = 1e-06 * 5.7e+13;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[1];
    Kc = exp((g_RT[5] + g_RT[11]) - (g_RT[9] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[86] = q_f - q_r;

    /*reaction 88: ch + oh <=> hco + h */
    phi_f = sc[5]*sc[12];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[1];
    Kc = exp((g_RT[5] + g_RT[12]) - (g_RT[7] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[87] = q_f - q_r;

    /*reaction 89: ch + co2 <=> hco + co */
    phi_f = sc[5]*sc[8];
    k_f = 1e-06 * 3.4e+12*exp(-347.2194994052565562/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[9];
    Kc = exp((g_RT[5] + g_RT[8]) - (g_RT[7] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[88] = q_f - q_r;

    /*reaction 90: ch + h2o <=> ch2o + h */
    phi_f = sc[5]*sc[15];
    k_f = 1e-06 * 1.17e+15*exp(-0.75*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[1];
    Kc = exp((g_RT[5] + g_RT[15]) - (g_RT[6] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[89] = q_f - q_r;

    /*reaction 91: ch + ch2o <=> ch2co + h */
    phi_f = sc[5]*sc[6];
    k_f = 1e-06 * 9.46e+13*exp(+259.15658288943063781/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[29]*sc[1];
    Kc = exp((g_RT[5] + g_RT[6]) - (g_RT[29] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[90] = q_f - q_r;

    /*reaction 92: ch + c2h2 <=> c3h2 + h */
    phi_f = sc[5]*sc[18];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[1];
    Kc = exp((g_RT[5] + g_RT[18]) - (g_RT[27] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[91] = q_f - q_r;

    /*reaction 93: ch + ch2 <=> c2h2 + h */
    phi_f = sc[5]*sc[4];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[1];
    Kc = exp((g_RT[5] + g_RT[4]) - (g_RT[18] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[92] = q_f - q_r;

    /*reaction 94: ch + ch3 <=> c2h3 + h */
    phi_f = sc[5]*sc[3];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[1];
    Kc = exp((g_RT[5] + g_RT[3]) - (g_RT[19] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[93] = q_f - q_r;

    /*reaction 95: ch + ch4 <=> c2h4 + h */
    phi_f = sc[5]*sc[2];
    k_f = 1e-06 * 6e+13;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[1];
    Kc = exp((g_RT[5] + g_RT[2]) - (g_RT[20] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[94] = q_f - q_r;

    /*reaction 96: ch2o + oh <=> hco + h2o */
    phi_f = sc[6]*sc[12];
    k_f = 1e-06 * 3.43e+09*exp(1.18*tc[0]+224.93784961470970529/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[15];
    Kc = exp((g_RT[6] + g_RT[12]) - (g_RT[7] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[95] = q_f - q_r;

    /*reaction 97: ch2o + h <=> hco + h2 */
    phi_f = sc[6]*sc[1];
    k_f = 1e-06 * 2.19e+08*exp(1.77*tc[0]-1509.6499974141590883/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[0];
    Kc = exp((g_RT[6] + g_RT[1]) - (g_RT[7] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[96] = q_f - q_r;

    /*reaction 98: ch2o + M <=> hco + h + M */
    phi_f = sc[6];
    alpha = mixture;
    k_f = 1e-06 * alpha * 3.31e+16*exp(-40760.549930182300159/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[1];
    Kc = refC * exp((g_RT[6]) - (g_RT[7] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[97] = q_f - q_r;

    /*reaction 99: ch2o + o <=> hco + oh */
    phi_f = sc[6]*sc[11];
    k_f = 1e-06 * 1.8e+13*exp(-1549.9073306785367095/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[12];
    Kc = exp((g_RT[6] + g_RT[11]) - (g_RT[7] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[98] = q_f - q_r;

    /*reaction 100: hco + o2 <=> co + ho2 */
    phi_f = sc[7]*sc[10];
    k_f = 1e-06 * 7.58e+12*exp(-206.31883297993508108/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[13];
    Kc = exp((g_RT[7] + g_RT[10]) - (g_RT[9] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[99] = q_f - q_r;

    /*reaction 101: hco + M <=> h + co + M */
    phi_f = sc[7];
    alpha = mixture + 4*sc[15] + 0.87*sc[0] + 2*sc[8] + 0.87*sc[9] + 1.81*sc[2];
    k_f = 1e-06 * alpha * 1.86e+17*exp(-1*tc[0]-8554.6833186802341515/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[9];
    Kc = refC * exp((g_RT[7]) - (g_RT[1] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[100] = q_f - q_r;

    /*reaction 102: hco + oh <=> h2o + co */
    phi_f = sc[7]*sc[12];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[9];
    Kc = exp((g_RT[7] + g_RT[12]) - (g_RT[15] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[101] = q_f - q_r;

    /*reaction 103: hco + h <=> co + h2 */
    phi_f = sc[7]*sc[1];
    k_f = 1e-06 * 1.19e+13*exp(0.25*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = exp((g_RT[7] + g_RT[1]) - (g_RT[9] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[102] = q_f - q_r;

    /*reaction 104: hco + o <=> co + oh */
    phi_f = sc[7]*sc[11];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[12];
    Kc = exp((g_RT[7] + g_RT[11]) - (g_RT[9] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[103] = q_f - q_r;

    /*reaction 105: hco + o <=> co2 + h */
    phi_f = sc[7]*sc[11];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[1];
    Kc = exp((g_RT[7] + g_RT[11]) - (g_RT[8] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[104] = q_f - q_r;

    /*reaction 106: co + oh <=> co2 + h */
    phi_f = sc[9]*sc[12];
    k_f = 1e-06 * 9420*exp(2.25*tc[0]+1183.0623813068959862/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[1];
    Kc = exp((g_RT[9] + g_RT[12]) - (g_RT[8] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[105] = q_f - q_r;

    /*reaction 107: co + o + M <=> co2 + M */
    phi_f = sc[9]*sc[11];
    alpha = mixture;
    k_f = 1e-12 * alpha * 6.17e+14*exp(-1509.6499974141590883/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8];
    Kc = 1.0 / (refC) * exp((g_RT[9] + g_RT[11]) - (g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[106] = q_f - q_r;

    /*reaction 108: co + o2 <=> co2 + o */
    phi_f = sc[9]*sc[10];
    k_f = 1e-06 * 2.53e+12*exp(-23997.396358895472076/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[11];
    Kc = exp((g_RT[9] + g_RT[10]) - (g_RT[8] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[107] = q_f - q_r;

    /*reaction 109: co + ho2 <=> co2 + oh */
    phi_f = sc[9]*sc[13];
    k_f = 1e-06 * 5.8e+13*exp(-11540.771013565441535/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[12];
    Kc = exp((g_RT[9] + g_RT[13]) - (g_RT[8] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[108] = q_f - q_r;

    /*reaction 110: c2h5oh (+M) <=> ch3 + ch2oh (+M) */
    phi_f = sc[33];
    alpha = mixture + 4*sc[15] + sc[0] + 2*sc[8] + sc[9];
    k_f = 1 * 5.94e+23*exp(-1.68*tc[0]-45874.740904755664815/tc[1]);
    redP = 1e-6 * alpha / k_f * 2.88e+85*exp(-18.9*tc[0]-55310.556605259960634/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-200))+ (0.5*exp(T/-890))+ (exp(-4600/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[23];
    Kc = refC * exp((g_RT[33]) - (g_RT[3] + g_RT[23]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[109] = q_f - q_r;

    /*reaction 111: c2h5oh (+M) <=> c2h5 + oh (+M) */
    phi_f = sc[33];
    alpha = mixture + 4*sc[15] + sc[0] + 2*sc[8] + sc[9];
    k_f = 1 * 1.25e+23*exp(-1.54*tc[0]-48311.316000582111883/tc[1]);
    redP = 1e-6 * alpha / k_f * 3.252e+85*exp(-18.81*tc[0]-57834.691400936433638/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-300))+ (0.5*exp(T/-900))+ (exp(-5000/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[12];
    Kc = refC * exp((g_RT[33]) - (g_RT[21] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[110] = q_f - q_r;

    /*reaction 112: c2h5oh (+M) <=> c2h4 + h2o (+M) */
    phi_f = sc[33];
    alpha = mixture + 4*sc[15];
    k_f = 1 * 2.79e+13*exp(0.09*tc[0]-33280.737409660941921/tc[1]);
    redP = 1e-6 * alpha / k_f * 2.57e+83*exp(-18.85*tc[0]-43504.087192149629118/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.3*exp(T/-350))+ (0.7*exp(T/-800))+ (exp(-3800/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[15];
    Kc = refC * exp((g_RT[33]) - (g_RT[20] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[111] = q_f - q_r;

    /*reaction 113: c2h5oh (+M) <=> ch3hco + h2 (+M) */
    phi_f = sc[33];
    alpha = mixture + 4*sc[15];
    k_f = 1 * 7.24e+11*exp(0.095*tc[0]-45796.239104890126328/tc[1]);
    redP = 1e-6 * alpha / k_f * 4.46e+87*exp(-19.42*tc[0]-58164.801533704332542/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.1*exp(T/-900))+ (0.9*exp(T/-1100))+ (exp(-3500/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[0];
    Kc = refC * exp((g_RT[33]) - (g_RT[32] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[112] = q_f - q_r;

    /*reaction 114: c2h5oh + oh <=> c2h4oh + h2o */
    phi_f = sc[33]*sc[12];
    k_f = 1e-06 * 1.74e+11*exp(0.27*tc[0]-301.92999948283181766/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[34]*sc[15];
    Kc = exp((g_RT[33] + g_RT[12]) - (g_RT[34] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[113] = q_f - q_r;

    /*reaction 115: c2h5oh + oh <=> ch3choh + h2o */
    phi_f = sc[33]*sc[12];
    k_f = 1e-06 * 4.64e+11*exp(0.15*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[15];
    Kc = exp((g_RT[33] + g_RT[12]) - (g_RT[35] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[114] = q_f - q_r;

    /*reaction 116: c2h5oh + oh <=> ch3ch2o + h2o */
    phi_f = sc[33]*sc[12];
    k_f = 1e-06 * 7.46e+11*exp(0.3*tc[0]-822.25603192491200844/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[36]*sc[15];
    Kc = exp((g_RT[33] + g_RT[12]) - (g_RT[36] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[115] = q_f - q_r;

    /*reaction 117: c2h5oh + h <=> c2h4oh + h2 */
    phi_f = sc[33]*sc[1];
    k_f = 1e-06 * 1.23e+07*exp(1.8*tc[0]-2565.3985622724612767/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[34]*sc[0];
    Kc = exp((g_RT[33] + g_RT[1]) - (g_RT[34] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[116] = q_f - q_r;

    /*reaction 118: c2h5oh + h <=> ch3choh + h2 */
    phi_f = sc[33]*sc[1];
    k_f = 1e-06 * 2.58e+07*exp(1.65*tc[0]-1422.593514229942457/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[0];
    Kc = exp((g_RT[33] + g_RT[1]) - (g_RT[35] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[117] = q_f - q_r;

    /*reaction 119: c2h5oh + h <=> ch3ch2o + h2 */
    phi_f = sc[33]*sc[1];
    k_f = 1e-06 * 1.5e+07*exp(1.6*tc[0]-1528.7722307147384981/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[36]*sc[0];
    Kc = exp((g_RT[33] + g_RT[1]) - (g_RT[36] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[118] = q_f - q_r;

    /*reaction 120: c2h5oh + o <=> c2h4oh + oh */
    phi_f = sc[33]*sc[11];
    k_f = 1e-06 * 9.41e+07*exp(1.7*tc[0]-2747.0597786279645334/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[34]*sc[12];
    Kc = exp((g_RT[33] + g_RT[11]) - (g_RT[34] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[119] = q_f - q_r;

    /*reaction 121: c2h5oh + o <=> ch3choh + oh */
    phi_f = sc[33]*sc[11];
    k_f = 1e-06 * 1.88e+07*exp(1.85*tc[0]-917.86719842780871659/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[12];
    Kc = exp((g_RT[33] + g_RT[11]) - (g_RT[35] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[120] = q_f - q_r;

    /*reaction 122: c2h5oh + o <=> ch3ch2o + oh */
    phi_f = sc[33]*sc[11];
    k_f = 1e-06 * 1.58e+07*exp(2*tc[0]-2238.3077294993931901/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[36]*sc[12];
    Kc = exp((g_RT[33] + g_RT[11]) - (g_RT[36] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[121] = q_f - q_r;

    /*reaction 123: c2h5oh + ch3 <=> c2h4oh + ch4 */
    phi_f = sc[33]*sc[3];
    k_f = 1e-06 * 219*exp(3.18*tc[0]-4841.9507583730128317/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[34]*sc[2];
    Kc = exp((g_RT[33] + g_RT[3]) - (g_RT[34] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[122] = q_f - q_r;

    /*reaction 124: c2h5oh + ch3 <=> ch3choh + ch4 */
    phi_f = sc[33]*sc[3];
    k_f = 1e-06 * 728*exp(2.99*tc[0]-3999.5660598159120127/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[2];
    Kc = exp((g_RT[33] + g_RT[3]) - (g_RT[35] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[123] = q_f - q_r;

    /*reaction 125: c2h5oh + ch3 <=> ch3ch2o + ch4 */
    phi_f = sc[33]*sc[3];
    k_f = 1e-06 * 145*exp(2.99*tc[0]-3849.1042767403009748/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[36]*sc[2];
    Kc = exp((g_RT[33] + g_RT[3]) - (g_RT[36] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[124] = q_f - q_r;

    /*reaction 126: c2h5oh + ho2 <=> ch3choh + h2o2 */
    phi_f = sc[33]*sc[13];
    k_f = 1e-06 * 8200*exp(2.55*tc[0]-5409.579157400737131/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[35]*sc[14];
    Kc = exp((g_RT[33] + g_RT[13]) - (g_RT[35] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[125] = q_f - q_r;

    /*reaction 127: c2h5oh + ho2 <=> c2h4oh + h2o2 */
    phi_f = sc[33]*sc[13];
    k_f = 1e-06 * 12300*exp(2.55*tc[0]-7925.6624864243349293/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[34]*sc[14];
    Kc = exp((g_RT[33] + g_RT[13]) - (g_RT[34] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[126] = q_f - q_r;

    /*reaction 128: c2h5oh + ho2 <=> ch3ch2o + h2o2 */
    phi_f = sc[33]*sc[13];
    k_f = 1e-06 * 2.5e+12*exp(-12077.199979313272706/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[36]*sc[14];
    Kc = exp((g_RT[33] + g_RT[13]) - (g_RT[36] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[127] = q_f - q_r;

    /*reaction 129: ch3ch2o + M <=> ch3hco + h + M */
    phi_f = sc[36];
    alpha = mixture;
    k_f = 1e-06 * alpha * 1.16e+35*exp(-5.89*tc[0]-12718.298011548486102/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[1];
    Kc = refC * exp((g_RT[36]) - (g_RT[32] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[128] = q_f - q_r;

    /*reaction 130: ch3ch2o + M <=> ch3 + ch2o + M */
    phi_f = sc[36];
    alpha = mixture;
    k_f = 1e-06 * alpha * 1.35e+38*exp(-6.96*tc[0]-11976.55664615232854/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[6];
    Kc = refC * exp((g_RT[36]) - (g_RT[3] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[129] = q_f - q_r;

    /*reaction 131: ch3ch2o + o2 <=> ch3hco + ho2 */
    phi_f = sc[36]*sc[10];
    k_f = 1e-06 * 4e+10*exp(-553.53833238519166571/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[13];
    Kc = exp((g_RT[36] + g_RT[10]) - (g_RT[32] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[130] = q_f - q_r;

    /*reaction 132: ch3ch2o + co <=> c2h5 + co2 */
    phi_f = sc[36]*sc[9];
    k_f = 1e-06 * 468*exp(3.16*tc[0]-2707.3056620293918968/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[8];
    Kc = exp((g_RT[36] + g_RT[9]) - (g_RT[21] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[131] = q_f - q_r;

    /*reaction 133: ch3ch2o + h <=> ch3 + ch2oh */
    phi_f = sc[36]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[23];
    Kc = exp((g_RT[36] + g_RT[1]) - (g_RT[3] + g_RT[23]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[132] = q_f - q_r;

    /*reaction 134: ch3ch2o + h <=> c2h4 + h2o */
    phi_f = sc[36]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[15];
    Kc = exp((g_RT[36] + g_RT[1]) - (g_RT[20] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[133] = q_f - q_r;

    /*reaction 135: ch3ch2o + oh <=> ch3hco + h2o */
    phi_f = sc[36]*sc[12];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[15];
    Kc = exp((g_RT[36] + g_RT[12]) - (g_RT[32] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[134] = q_f - q_r;

    /*reaction 136: ch3choh + o2 <=> ch3hco + ho2 */
    phi_f = sc[35]*sc[10];
    k_f = 1e-06 * 4.82e+14*exp(-2524.638012342278671/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[13];
    Kc = exp((g_RT[35] + g_RT[10]) - (g_RT[32] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[135] = q_f - q_r;

    /*reaction 137: ch3choh + o2 <=> ch3hco + ho2 */
    phi_f = sc[35]*sc[10];
    k_f = 1e-06 * 8.43e+15*exp(-1.2*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[13];
    Kc = exp((g_RT[35] + g_RT[10]) - (g_RT[32] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[136] = q_f - q_r;

    /*reaction 138: ch3choh + o <=> ch3hco + oh */
    phi_f = sc[35]*sc[11];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[12];
    Kc = exp((g_RT[35] + g_RT[11]) - (g_RT[32] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[137] = q_f - q_r;

    /*reaction 139: ch3choh + h <=> c2h4 + h2o */
    phi_f = sc[35]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[15];
    Kc = exp((g_RT[35] + g_RT[1]) - (g_RT[20] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[138] = q_f - q_r;

    /*reaction 140: ch3choh + h <=> ch3 + ch2oh */
    phi_f = sc[35]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[23];
    Kc = exp((g_RT[35] + g_RT[1]) - (g_RT[3] + g_RT[23]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[139] = q_f - q_r;

    /*reaction 141: ch3choh + ho2 <=> ch3hco + oh + oh */
    phi_f = sc[35]*sc[13];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[12]*sc[12];
    Kc = refC * exp((g_RT[35] + g_RT[13]) - (g_RT[32] + g_RT[12] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[140] = q_f - q_r;

    /*reaction 142: ch3choh + oh <=> ch3hco + h2o */
    phi_f = sc[35]*sc[12];
    k_f = 1e-06 * 5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[15];
    Kc = exp((g_RT[35] + g_RT[12]) - (g_RT[32] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[141] = q_f - q_r;

    /*reaction 143: ch3choh + M <=> ch3hco + h + M */
    phi_f = sc[35];
    alpha = mixture;
    k_f = 1e-06 * alpha * 1e+14*exp(-12580.416645117993539/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[1];
    Kc = refC * exp((g_RT[35]) - (g_RT[32] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[142] = q_f - q_r;

    /*reaction 144: ch3hco + oh <=> ch3co + h2o */
    phi_f = sc[32]*sc[12];
    k_f = 1e-06 * 9.24e+06*exp(1.5*tc[0]+484.0944325041403431/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[15];
    Kc = exp((g_RT[32] + g_RT[12]) - (g_RT[31] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[143] = q_f - q_r;

    /*reaction 145: ch3hco + oh <=> ch2hco + h2o */
    phi_f = sc[32]*sc[12];
    k_f = 1e-06 * 172000*exp(2.4*tc[0]-410.12158263084654664/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[15];
    Kc = exp((g_RT[32] + g_RT[12]) - (g_RT[30] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[144] = q_f - q_r;

    /*reaction 146: ch3hco + o <=> ch3co + oh */
    phi_f = sc[32]*sc[11];
    k_f = 1e-06 * 1.77e+18*exp(-1.9*tc[0]-1497.0695807690410675/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[12];
    Kc = exp((g_RT[32] + g_RT[11]) - (g_RT[31] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[145] = q_f - q_r;

    /*reaction 147: ch3hco + o <=> ch2hco + oh */
    phi_f = sc[32]*sc[11];
    k_f = 1e-06 * 3.72e+13*exp(-0.2*tc[0]-1789.4384636015831802/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[12];
    Kc = exp((g_RT[32] + g_RT[11]) - (g_RT[30] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[146] = q_f - q_r;

    /*reaction 148: ch3hco + h <=> ch3co + h2 */
    phi_f = sc[32]*sc[1];
    k_f = 1e-06 * 4.66e+13*exp(-0.35*tc[0]-1503.6113974245024565/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[0];
    Kc = exp((g_RT[32] + g_RT[1]) - (g_RT[31] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[147] = q_f - q_r;

    /*reaction 149: ch3hco + h <=> ch2hco + h2 */
    phi_f = sc[32]*sc[1];
    k_f = 1e-06 * 1.85e+12*exp(0.4*tc[0]-2696.7381120474929048/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[0];
    Kc = exp((g_RT[32] + g_RT[1]) - (g_RT[30] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[148] = q_f - q_r;

    /*reaction 150: ch3hco + ch3 <=> ch3co + ch4 */
    phi_f = sc[32]*sc[3];
    k_f = 1e-06 * 3.9e-07*exp(5.8*tc[0]-1107.0766647703833314/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[2];
    Kc = exp((g_RT[32] + g_RT[3]) - (g_RT[31] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[149] = q_f - q_r;

    /*reaction 151: ch3hco + ch3 <=> ch2hco + ch4 */
    phi_f = sc[32]*sc[3];
    k_f = 1e-06 * 24.5*exp(3.15*tc[0]-2881.9218450636299167/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[2];
    Kc = exp((g_RT[32] + g_RT[3]) - (g_RT[30] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[150] = q_f - q_r;

    /*reaction 152: ch3hco + ho2 <=> ch3co + h2o2 */
    phi_f = sc[32]*sc[13];
    k_f = 1e-06 * 2.4e+19*exp(-2.2*tc[0]-7060.1298212402170975/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[14];
    Kc = exp((g_RT[32] + g_RT[13]) - (g_RT[31] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[151] = q_f - q_r;

    /*reaction 153: ch3hco + ho2 <=> ch2hco + h2o2 */
    phi_f = sc[32]*sc[13];
    k_f = 1e-06 * 2.32e+11*exp(0.4*tc[0]-7479.8125205213536901/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[14];
    Kc = exp((g_RT[32] + g_RT[13]) - (g_RT[30] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[152] = q_f - q_r;

    /*reaction 154: ch3hco + o2 <=> ch3co + ho2 */
    phi_f = sc[32]*sc[10];
    k_f = 1e-06 * 1e+14*exp(-21235.743296959171857/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[13];
    Kc = exp((g_RT[32] + g_RT[10]) - (g_RT[31] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[153] = q_f - q_r;

    /*reaction 155: c2h6 + ch3 <=> c2h5 + ch4 */
    phi_f = sc[22]*sc[3];
    k_f = 1e-06 * 0.55*exp(4*tc[0]-4176.6983261791738187/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[2];
    Kc = exp((g_RT[22] + g_RT[3]) - (g_RT[21] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[154] = q_f - q_r;

    /*reaction 156: c2h6 + h <=> c2h5 + h2 */
    phi_f = sc[22]*sc[1];
    k_f = 1e-06 * 540*exp(3.5*tc[0]-2621.7588288425899918/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[0];
    Kc = exp((g_RT[22] + g_RT[1]) - (g_RT[21] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[155] = q_f - q_r;

    /*reaction 157: c2h6 + o <=> c2h5 + oh */
    phi_f = sc[22]*sc[11];
    k_f = 1e-06 * 3e+07*exp(2*tc[0]-2573.9532455911412399/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[12];
    Kc = exp((g_RT[22] + g_RT[11]) - (g_RT[21] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[156] = q_f - q_r;

    /*reaction 158: c2h6 + oh <=> c2h5 + h2o */
    phi_f = sc[22]*sc[12];
    k_f = 1e-06 * 7.23e+06*exp(2*tc[0]-434.77919925527783107/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[15];
    Kc = exp((g_RT[22] + g_RT[12]) - (g_RT[21] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[157] = q_f - q_r;

    /*reaction 159: c2h5 + h <=> c2h4 + h2 */
    phi_f = sc[21]*sc[1];
    k_f = 1e-06 * 1.25e+14*exp(-4025.733326437757114/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[0];
    Kc = exp((g_RT[21] + g_RT[1]) - (g_RT[20] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[158] = q_f - q_r;

    /*reaction 160: c2h5 + h <=> ch3 + ch3 */
    phi_f = sc[21]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[3];
    Kc = exp((g_RT[21] + g_RT[1]) - (g_RT[3] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[159] = q_f - q_r;

    /*reaction 161: c2h5 + h <=> c2h6 */
    phi_f = sc[21]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[22];
    Kc = 1.0 / (refC) * exp((g_RT[21] + g_RT[1]) - (g_RT[22]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[160] = q_f - q_r;

    /*reaction 162: c2h5 + oh <=> c2h4 + h2o */
    phi_f = sc[21]*sc[12];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[15];
    Kc = exp((g_RT[21] + g_RT[12]) - (g_RT[20] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[161] = q_f - q_r;

    /*reaction 163: c2h5 + o <=> ch3 + ch2o */
    phi_f = sc[21]*sc[11];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[6];
    Kc = exp((g_RT[21] + g_RT[11]) - (g_RT[3] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[162] = q_f - q_r;

    /*reaction 164: c2h5 + ho2 <=> c2h6 + o2 */
    phi_f = sc[21]*sc[13];
    k_f = 1e-06 * 3e+12;
    q_f = phi_f * k_f;
    phi_r = sc[22]*sc[10];
    Kc = exp((g_RT[21] + g_RT[13]) - (g_RT[22] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[163] = q_f - q_r;

    /*reaction 165: c2h5 + ho2 <=> ch3ch2o + oh */
    phi_f = sc[21]*sc[13];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[36]*sc[12];
    Kc = exp((g_RT[21] + g_RT[13]) - (g_RT[36] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[164] = q_f - q_r;

    /*reaction 166: c2h5 + o2 <=> c2h4 + ho2 */
    phi_f = sc[21]*sc[10];
    k_f = 1e-06 * 2.89e+28*exp(-5.4*tc[0]-3816.8984101287987869/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[20]*sc[13];
    Kc = exp((g_RT[21] + g_RT[10]) - (g_RT[20] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[165] = q_f - q_r;

    /*reaction 167: c2h5 + o2 <=> ch3hco + oh */
    phi_f = sc[21]*sc[10];
    k_f = 1e-06 * 4.9e+11*exp(-0.48*tc[0]-4205.3816761300422513/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[12];
    Kc = exp((g_RT[21] + g_RT[10]) - (g_RT[32] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[166] = q_f - q_r;

    /*reaction 168: c2h4 + oh <=> c2h4oh */
    phi_f = sc[20]*sc[12];
    k_f = 1e-06 * 1.29e+12*exp(+411.12801596245600422/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[34];
    Kc = 1.0 / (refC) * exp((g_RT[20] + g_RT[12]) - (g_RT[34]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[167] = q_f - q_r;

    /*reaction 169: c2h4 + oh <=> c2h3 + h2o */
    phi_f = sc[20]*sc[12];
    k_f = 1e-06 * 2.02e+13*exp(-2987.0941282168159887/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[15];
    Kc = exp((g_RT[20] + g_RT[12]) - (g_RT[19] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[168] = q_f - q_r;

    /*reaction 170: c2h4 + o <=> ch3 + hco */
    phi_f = sc[20]*sc[11];
    k_f = 1e-06 * 1.02e+07*exp(1.88*tc[0]-90.07578317904483356/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = exp((g_RT[20] + g_RT[11]) - (g_RT[3] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[169] = q_f - q_r;

    /*reaction 171: c2h4 + o <=> ch2hco + h */
    phi_f = sc[20]*sc[11];
    k_f = 1e-06 * 3.39e+06*exp(1.88*tc[0]-90.07578317904483356/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[1];
    Kc = exp((g_RT[20] + g_RT[11]) - (g_RT[30] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[170] = q_f - q_r;

    /*reaction 172: c2h4 + ch3 <=> c2h3 + ch4 */
    phi_f = sc[20]*sc[3];
    k_f = 1e-06 * 6.62*exp(3.7*tc[0]-4780.5583251448369992/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[2];
    Kc = exp((g_RT[20] + g_RT[3]) - (g_RT[19] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[171] = q_f - q_r;

    /*reaction 173: c2h4 + h <=> c2h3 + h2 */
    phi_f = sc[20]*sc[1];
    k_f = 1e-06 * 3.36e-07*exp(6*tc[0]-851.44259854158576672/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[0];
    Kc = exp((g_RT[20] + g_RT[1]) - (g_RT[19] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[172] = q_f - q_r;

    /*reaction 174: c2h4 + h (+M) <=> c2h5 (+M) */
    phi_f = sc[20]*sc[1];
    alpha = mixture + 4*sc[15] + sc[0] + 2*sc[8] + sc[9];
    k_f = 1e-06 * 1.08e+12*exp(0.454*tc[0]-916.86076509619920216/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.112e+34*exp(-5*tc[0]-2238.3077294993931901/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0*exp(T/-1e-15))+ (1*exp(T/-95))+ (exp(-200/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[21];
    Kc = 1.0 / (refC) * exp((g_RT[20] + g_RT[1]) - (g_RT[21]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[173] = q_f - q_r;

    /*reaction 175: c2h4 (+M) <=> c2h2 + h2 (+M) */
    phi_f = sc[20];
    alpha = mixture;
    k_f = 1 * 1.8e+14*exp(-43779.849925010610605/tc[1]);
    redP = 1e-6 * alpha / k_f * 1.5e+15*exp(-27899.841602211075951/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[0];
    Kc = refC * exp((g_RT[20]) - (g_RT[18] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[174] = q_f - q_r;

    /*reaction 176: c2h3 + h (+M) <=> c2h4 (+M) */
    phi_f = sc[19]*sc[1];
    alpha = mixture + 4*sc[15];
    k_f = 1e-06 * 6.1e+12*exp(0.27*tc[0]-140.90066642532150354/tc[1]);
    redP = 1e-12 * alpha / k_f * 9.8e+29*exp(-3.86*tc[0]-1670.6793304716693456/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.218*exp(T/-208))+ (0.782*exp(T/-2663))+ (exp(-6095/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[20];
    Kc = 1.0 / (refC) * exp((g_RT[19] + g_RT[1]) - (g_RT[20]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[175] = q_f - q_r;

    /*reaction 177: c2h3 + h <=> c2h2 + h2 */
    phi_f = sc[19]*sc[1];
    k_f = 1e-06 * 9e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[0];
    Kc = exp((g_RT[19] + g_RT[1]) - (g_RT[18] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[176] = q_f - q_r;

    /*reaction 178: c2h3 + o <=> ch2co + h */
    phi_f = sc[19]*sc[11];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[29]*sc[1];
    Kc = exp((g_RT[19] + g_RT[11]) - (g_RT[29] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[177] = q_f - q_r;

    /*reaction 179: c2h3 + o2 <=> ch2o + hco */
    phi_f = sc[19]*sc[10];
    k_f = 1e-06 * 1.7e+29*exp(-5.312*tc[0]-3270.9083277306781383/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[7];
    Kc = exp((g_RT[19] + g_RT[10]) - (g_RT[6] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[178] = q_f - q_r;

    /*reaction 180: c2h3 + o2 <=> ch2hco + o */
    phi_f = sc[19]*sc[10];
    k_f = 1e-06 * 5.5e+14*exp(-0.611*tc[0]-2646.9196621328255787/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[30]*sc[11];
    Kc = exp((g_RT[19] + g_RT[10]) - (g_RT[30] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[179] = q_f - q_r;

    /*reaction 181: c2h3 + o2 <=> c2h2 + ho2 */
    phi_f = sc[19]*sc[10];
    k_f = 1e-06 * 2.12e-06*exp(6*tc[0]-4772.5068584919617933/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[13];
    Kc = exp((g_RT[19] + g_RT[10]) - (g_RT[18] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[180] = q_f - q_r;

    /*reaction 182: c2h3 + oh <=> c2h2 + h2o */
    phi_f = sc[19]*sc[12];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[15];
    Kc = exp((g_RT[19] + g_RT[12]) - (g_RT[18] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[181] = q_f - q_r;

    /*reaction 183: c2h3 + c2h <=> c2h2 + c2h2 */
    phi_f = sc[19]*sc[16];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[18];
    Kc = exp((g_RT[19] + g_RT[16]) - (g_RT[18] + g_RT[18]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[182] = q_f - q_r;

    /*reaction 184: c2h3 + ch <=> ch2 + c2h2 */
    phi_f = sc[19]*sc[5];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[18];
    Kc = exp((g_RT[19] + g_RT[5]) - (g_RT[4] + g_RT[18]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[183] = q_f - q_r;

    /*reaction 185: c2h3 + ch3 <=> c2h2 + ch4 */
    phi_f = sc[19]*sc[3];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[2];
    Kc = exp((g_RT[19] + g_RT[3]) - (g_RT[18] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[184] = q_f - q_r;

    /*reaction 186: c2h2 + oh <=> c2h + h2o */
    phi_f = sc[18]*sc[12];
    k_f = 1e-06 * 3.37e+07*exp(2*tc[0]-7045.0333212660752906/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[15];
    Kc = exp((g_RT[18] + g_RT[12]) - (g_RT[16] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[185] = q_f - q_r;

    /*reaction 187: c2h2 + oh <=> hccoh + h */
    phi_f = sc[18]*sc[12];
    k_f = 1e-06 * 504000*exp(2.3*tc[0]-6793.4249883637157836/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[1];
    Kc = exp((g_RT[18] + g_RT[12]) - (g_RT[25] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[186] = q_f - q_r;

    /*reaction 188: c2h2 + oh <=> ch2co + h */
    phi_f = sc[18]*sc[12];
    k_f = 1e-06 * 0.000218*exp(4.5*tc[0]+503.21666580471963925/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[29]*sc[1];
    Kc = exp((g_RT[18] + g_RT[12]) - (g_RT[29] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[187] = q_f - q_r;

    /*reaction 189: c2h2 + oh <=> ch2co + h */
    phi_f = sc[18]*sc[12];
    k_f = 1e-06 * 2e+11;
    q_f = phi_f * k_f;
    phi_r = sc[29]*sc[1];
    Kc = exp((g_RT[18] + g_RT[12]) - (g_RT[29] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[188] = q_f - q_r;

    /*reaction 190: c2h2 + oh <=> ch3 + co */
    phi_f = sc[18]*sc[12];
    k_f = 1e-06 * 0.000483*exp(4*tc[0]+1006.4333316094392785/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[9];
    Kc = exp((g_RT[18] + g_RT[12]) - (g_RT[3] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[189] = q_f - q_r;

    /*reaction 191: hccoh + h <=> ch2co + h */
    phi_f = sc[25]*sc[1];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[29]*sc[1];
    Kc = exp((g_RT[25] + g_RT[1]) - (g_RT[29] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[190] = q_f - q_r;

    /*reaction 192: c2h2 + o <=> ch2 + co */
    phi_f = sc[18]*sc[11];
    k_f = 1e-06 * 6.12e+06*exp(2*tc[0]-956.11166502896742259/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = exp((g_RT[18] + g_RT[11]) - (g_RT[4] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[191] = q_f - q_r;

    /*reaction 193: c2h2 + o <=> hcco + h */
    phi_f = sc[18]*sc[11];
    k_f = 1e-06 * 1.43e+07*exp(2*tc[0]-956.11166502896742259/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[1];
    Kc = exp((g_RT[18] + g_RT[11]) - (g_RT[17] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[192] = q_f - q_r;

    /*reaction 194: c2h2 + o <=> c2h + oh */
    phi_f = sc[18]*sc[11];
    k_f = 1e-06 * 3.16e+15*exp(-0.6*tc[0]-7548.2499870707952141/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[12];
    Kc = exp((g_RT[18] + g_RT[11]) - (g_RT[16] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[193] = q_f - q_r;

    /*reaction 195: c2h2 + ch3 <=> c2h + ch4 */
    phi_f = sc[18]*sc[3];
    k_f = 1e-06 * 1.81e+11*exp(-8700.1129350977989816/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[2];
    Kc = exp((g_RT[18] + g_RT[3]) - (g_RT[16] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[194] = q_f - q_r;

    /*reaction 196: c2h2 + o2 <=> hcco + oh */
    phi_f = sc[18]*sc[10];
    k_f = 1e-06 * 4e+07*exp(1.5*tc[0]-15146.821640722062511/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[12];
    Kc = exp((g_RT[18] + g_RT[10]) - (g_RT[17] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[195] = q_f - q_r;

    /*reaction 197: c2h2 + M <=> c2h + h + M */
    phi_f = sc[18];
    alpha = mixture;
    k_f = 1e-06 * alpha * 4.2e+16*exp(-53844.183241105005436/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[1];
    Kc = refC * exp((g_RT[18]) - (g_RT[16] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[196] = q_f - q_r;

    /*reaction 198: c2h2 + h (+M) <=> c2h3 (+M) */
    phi_f = sc[18]*sc[1];
    alpha = mixture + 4*sc[15] + sc[0] + 2*sc[8] + sc[9];
    k_f = 1e-06 * 3.11e+11*exp(0.58*tc[0]-1302.8279477684191079/tc[1]);
    redP = 1e-12 * alpha / k_f * 2.25e+40*exp(-7.269*tc[0]-3309.6560109976417152/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0*exp(T/-1e-15))+ (1*exp(T/-675))+ (exp(-1e+15/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[19];
    Kc = 1.0 / (refC) * exp((g_RT[18] + g_RT[1]) - (g_RT[19]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[197] = q_f - q_r;

    /*reaction 199: ch2hco + h <=> ch3 + hco */
    phi_f = sc[30]*sc[1];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = exp((g_RT[30] + g_RT[1]) - (g_RT[3] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[198] = q_f - q_r;

    /*reaction 200: ch2hco + h <=> ch2co + h2 */
    phi_f = sc[30]*sc[1];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[29]*sc[0];
    Kc = exp((g_RT[30] + g_RT[1]) - (g_RT[29] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[199] = q_f - q_r;

    /*reaction 201: ch2hco + o <=> ch2o + hco */
    phi_f = sc[30]*sc[11];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[7];
    Kc = exp((g_RT[30] + g_RT[11]) - (g_RT[6] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[200] = q_f - q_r;

    /*reaction 202: ch2hco + oh <=> ch2co + h2o */
    phi_f = sc[30]*sc[12];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[29]*sc[15];
    Kc = exp((g_RT[30] + g_RT[12]) - (g_RT[29] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[201] = q_f - q_r;

    /*reaction 203: ch2hco + o2 <=> ch2o + co + oh */
    phi_f = sc[30]*sc[10];
    k_f = 1e-06 * 3e+10;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[9]*sc[12];
    Kc = refC * exp((g_RT[30] + g_RT[10]) - (g_RT[6] + g_RT[9] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[202] = q_f - q_r;

    /*reaction 204: ch2hco + ch3 <=> c2h5 + co + h */
    phi_f = sc[30]*sc[3];
    k_f = 1e-06 * 4.9e+14*exp(-0.5*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[9]*sc[1];
    Kc = refC * exp((g_RT[30] + g_RT[3]) - (g_RT[21] + g_RT[9] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[203] = q_f - q_r;

    /*reaction 205: ch2hco + ho2 <=> ch2o + hco + oh */
    phi_f = sc[30]*sc[13];
    k_f = 1e-06 * 7e+12;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[7]*sc[12];
    Kc = refC * exp((g_RT[30] + g_RT[13]) - (g_RT[6] + g_RT[7] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[204] = q_f - q_r;

    /*reaction 206: ch2hco + ho2 <=> ch3hco + o2 */
    phi_f = sc[30]*sc[13];
    k_f = 1e-06 * 3e+12;
    q_f = phi_f * k_f;
    phi_r = sc[32]*sc[10];
    Kc = exp((g_RT[30] + g_RT[13]) - (g_RT[32] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[205] = q_f - q_r;

    /*reaction 207: ch2hco <=> ch3 + co */
    phi_f = sc[30];
    k_f = 1 * 1.17e+43*exp(-9.83*tc[0]-22018.748428951312235/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[9];
    Kc = refC * exp((g_RT[30]) - (g_RT[3] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[206] = q_f - q_r;

    /*reaction 208: ch2hco <=> ch2co + h */
    phi_f = sc[30];
    k_f = 1 * 1.81e+43*exp(-9.61*tc[0]-23081.542027130883071/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[29]*sc[1];
    Kc = refC * exp((g_RT[30]) - (g_RT[29] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[207] = q_f - q_r;

    /*reaction 209: ch3co (+M) <=> ch3 + co (+M) */
    phi_f = sc[31];
    alpha = mixture;
    k_f = 1 * 3e+12*exp(-8414.7890855865243793/tc[1]);
    redP = 1e-6 * alpha / k_f * 1.2e+15*exp(-6299.2662225434805805/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[9];
    Kc = refC * exp((g_RT[31]) - (g_RT[3] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[208] = q_f - q_r;

    /*reaction 210: ch2co + o <=> co2 + ch2 */
    phi_f = sc[29]*sc[11];
    k_f = 1e-06 * 1.75e+12*exp(-679.34249883637164658/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[4];
    Kc = exp((g_RT[29] + g_RT[11]) - (g_RT[8] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[209] = q_f - q_r;

    /*reaction 211: ch2co + h <=> ch3 + co */
    phi_f = sc[29]*sc[1];
    k_f = 1e-06 * 27100*exp(2.75*tc[0]-359.29669938456987666/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[9];
    Kc = exp((g_RT[29] + g_RT[1]) - (g_RT[3] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[210] = q_f - q_r;

    /*reaction 212: ch2co + h <=> hcco + h2 */
    phi_f = sc[29]*sc[1];
    k_f = 1e-06 * 2e+14*exp(-4025.733326437757114/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[0];
    Kc = exp((g_RT[29] + g_RT[1]) - (g_RT[17] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[211] = q_f - q_r;

    /*reaction 213: ch2co + o <=> hcco + oh */
    phi_f = sc[29]*sc[11];
    k_f = 1e-06 * 1e+13*exp(-4025.733326437757114/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[12];
    Kc = exp((g_RT[29] + g_RT[11]) - (g_RT[17] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[212] = q_f - q_r;

    /*reaction 214: ch2co + oh <=> hcco + h2o */
    phi_f = sc[29]*sc[12];
    k_f = 1e-06 * 1e+13*exp(-1006.4333316094392785/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[15];
    Kc = exp((g_RT[29] + g_RT[12]) - (g_RT[17] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[213] = q_f - q_r;

    /*reaction 215: ch2co + oh <=> ch2oh + co */
    phi_f = sc[29]*sc[12];
    k_f = 1e-06 * 3.73e+12*exp(+509.75848246018102827/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[9];
    Kc = exp((g_RT[29] + g_RT[12]) - (g_RT[23] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[214] = q_f - q_r;

    /*reaction 216: ch2co (+M) <=> ch2 + co (+M) */
    phi_f = sc[29];
    alpha = mixture;
    k_f = 1 * 3e+14*exp(-35718.318938819007599/tc[1]);
    redP = 1e-6 * alpha / k_f * 3.6e+15*exp(-29825.651782245739014/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = refC * exp((g_RT[29]) - (g_RT[4] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[215] = q_f - q_r;

    /*reaction 217: c2h + h2 <=> c2h2 + h */
    phi_f = sc[16]*sc[0];
    k_f = 1e-06 * 409000*exp(2.39*tc[0]-434.93016425501917865/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[1];
    Kc = exp((g_RT[16] + g_RT[0]) - (g_RT[18] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[216] = q_f - q_r;

    /*reaction 218: c2h + o <=> ch + co */
    phi_f = sc[16]*sc[11];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[9];
    Kc = exp((g_RT[16] + g_RT[11]) - (g_RT[5] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[217] = q_f - q_r;

    /*reaction 219: c2h + oh <=> hcco + h */
    phi_f = sc[16]*sc[12];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[1];
    Kc = exp((g_RT[16] + g_RT[12]) - (g_RT[17] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[218] = q_f - q_r;

    /*reaction 220: c2h + o2 <=> co + co + h */
    phi_f = sc[16]*sc[10];
    k_f = 1e-06 * 9.04e+12*exp(+229.97001627275690794/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[9]*sc[1];
    Kc = refC * exp((g_RT[16] + g_RT[10]) - (g_RT[9] + g_RT[9] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[219] = q_f - q_r;

    /*reaction 221: hcco + c2h2 <=> h2ccch + co */
    phi_f = sc[17]*sc[18];
    k_f = 1e-06 * 1e+11*exp(-1509.6499974141590883/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[26]*sc[9];
    Kc = exp((g_RT[17] + g_RT[18]) - (g_RT[26] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[220] = q_f - q_r;

    /*reaction 222: hcco + h <=> ch2s + co */
    phi_f = sc[17]*sc[1];
    k_f = 1e-06 * 1e+14;
    q_f = phi_f * k_f;
    phi_r = sc[28]*sc[9];
    Kc = exp((g_RT[17] + g_RT[1]) - (g_RT[28] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[221] = q_f - q_r;

    /*reaction 223: hcco + o <=> h + co + co */
    phi_f = sc[17]*sc[11];
    k_f = 1e-06 * 8e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[9]*sc[9];
    Kc = refC * exp((g_RT[17] + g_RT[11]) - (g_RT[1] + g_RT[9] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[222] = q_f - q_r;

    /*reaction 224: hcco + o <=> ch + co2 */
    phi_f = sc[17]*sc[11];
    k_f = 1e-06 * 2.95e+13*exp(-560.08014904065305473/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[8];
    Kc = exp((g_RT[17] + g_RT[11]) - (g_RT[5] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[223] = q_f - q_r;

    /*reaction 225: hcco + o2 <=> hco + co + o */
    phi_f = sc[17]*sc[10];
    k_f = 1e-06 * 2.5e+08*exp(1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[9]*sc[11];
    Kc = refC * exp((g_RT[17] + g_RT[10]) - (g_RT[7] + g_RT[9] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[224] = q_f - q_r;

    /*reaction 226: hcco + o2 <=> co2 + hco */
    phi_f = sc[17]*sc[10];
    k_f = 1e-06 * 2.4e+11*exp(+429.7470325972306/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[7];
    Kc = exp((g_RT[17] + g_RT[10]) - (g_RT[8] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[225] = q_f - q_r;

    /*reaction 227: hcco + ch <=> c2h2 + co */
    phi_f = sc[17]*sc[5];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[9];
    Kc = exp((g_RT[17] + g_RT[5]) - (g_RT[18] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[226] = q_f - q_r;

    /*reaction 228: hcco + hcco <=> c2h2 + co + co */
    phi_f = sc[17]*sc[17];
    k_f = 1e-06 * 1e+13;
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[9]*sc[9];
    Kc = refC * exp((g_RT[17] + g_RT[17]) - (g_RT[18] + g_RT[9] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[227] = q_f - q_r;

    return;
}


/*compute the progress rate for each reaction */
void progressRateFR(double * q_f, double * q_r, double * sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[38];                /*Gibbs free energy */
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
    double refC = 101325 / 8.31451 / T;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 38; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*reaction 1: oh + h2 <=> h + h2o */
    phi_f = sc[12]*sc[0];
    k_f = 1e-06 * 2.14e+08*exp(1.52*tc[0]-1735.5942803604782512/tc[1]);
    q_f[0] = phi_f * k_f;
    phi_r = sc[1]*sc[15];
    Kc = exp((g_RT[12] + g_RT[0]) - (g_RT[1] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[0] = phi_r * k_r;

    /*reaction 2: o + oh <=> o2 + h */
    phi_f = sc[11]*sc[12];
    k_f = 1e-06 * 2.02e+14*exp(-0.4*tc[0]);
    q_f[1] = phi_f * k_f;
    phi_r = sc[10]*sc[1];
    Kc = exp((g_RT[11] + g_RT[12]) - (g_RT[10] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[1] = phi_r * k_r;

    /*reaction 3: o + h2 <=> oh + h */
    phi_f = sc[11]*sc[0];
    k_f = 1e-06 * 50600*exp(2.67*tc[0]-3165.2328279116868543/tc[1]);
    q_f[2] = phi_f * k_f;
    phi_r = sc[12]*sc[1];
    Kc = exp((g_RT[11] + g_RT[0]) - (g_RT[12] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[2] = phi_r * k_r;

    /*reaction 4: h + o2 (+M) <=> ho2 (+M) */
    phi_f = sc[1]*sc[10];
    alpha = mixture + -1*sc[15] + -1*sc[0] + -1*sc[37] + 9*sc[2] + 2.8*sc[8] + 0.9*sc[9];
    k_f = 1e-06 * 4.52e+13;
    redP = 1e-12 * alpha / k_f * 1.05e+19*exp(-1.257*tc[0]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f[3] = phi_f * k_f;
    phi_r = sc[13];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[13]));
    k_r = k_f / Kc;
    q_r[3] = phi_r * k_r;

    /*reaction 5: h + o2 (+n2) <=> ho2 (+n2) */
    phi_f = sc[1]*sc[10];
    alpha = sc[37];
    k_f = 1e-06 * 4.52e+13;
    redP = 1e-12 * alpha / k_f * 2.03e+20*exp(-1.59*tc[0]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f[4] = phi_f * k_f;
    phi_r = sc[13];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[13]));
    k_r = k_f / Kc;
    q_r[4] = phi_r * k_r;

    /*reaction 6: h + o2 (+h2) <=> ho2 (+h2) */
    phi_f = sc[1]*sc[10];
    alpha = sc[0];
    k_f = 1e-06 * 4.52e+13;
    redP = 1e-12 * alpha / k_f * 1.52e+19*exp(-1.133*tc[0]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f[5] = phi_f * k_f;
    phi_r = sc[13];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[13]));
    k_r = k_f / Kc;
    q_r[5] = phi_r * k_r;

    /*reaction 7: h + o2 (+h2o) <=> ho2 (+h2o) */
    phi_f = sc[1]*sc[10];
    alpha = sc[15];
    k_f = 1e-06 * 4.52e+13;
    redP = 1e-12 * alpha / k_f * 2.1e+23*exp(-2.437*tc[0]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f[6] = phi_f * k_f;
    phi_r = sc[13];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[13]));
    k_r = k_f / Kc;
    q_r[6] = phi_r * k_r;

    /*reaction 8: oh + ho2 <=> h2o + o2 */
    phi_f = sc[12]*sc[13];
    k_f = 1e-06 * 2.13e+28*exp(-4.827*tc[0]-1761.2583303165188227/tc[1]);
    q_f[7] = phi_f * k_f;
    phi_r = sc[15]*sc[10];
    Kc = exp((g_RT[12] + g_RT[13]) - (g_RT[15] + g_RT[10]));
    k_r = k_f / Kc;
    q_r[7] = phi_r * k_r;

    /*reaction 9: oh + ho2 <=> h2o + o2 */
    phi_f = sc[12]*sc[13];
    k_f = 1e-06 * 9.1e+14*exp(-5517.267523882946989/tc[1]);
    q_f[8] = phi_f * k_f;
    phi_r = sc[15]*sc[10];
    Kc = exp((g_RT[12] + g_RT[13]) - (g_RT[15] + g_RT[10]));
    k_r = k_f / Kc;
    q_r[8] = phi_r * k_r;

    /*reaction 10: h + ho2 <=> oh + oh */
    phi_f = sc[1]*sc[13];
    k_f = 1e-06 * 1.5e+14*exp(-503.21666580471963925/tc[1]);
    q_f[9] = phi_f * k_f;
    phi_r = sc[12]*sc[12];
    Kc = exp((g_RT[1] + g_RT[13]) - (g_RT[12] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[9] = phi_r * k_r;

    /*reaction 11: h + ho2 <=> h2 + o2 */
    phi_f = sc[1]*sc[13];
    k_f = 1e-06 * 6.63e+13*exp(-1069.8386315008340262/tc[1]);
    q_f[10] = phi_f * k_f;
    phi_r = sc[0]*sc[10];
    Kc = exp((g_RT[1] + g_RT[13]) - (g_RT[0] + g_RT[10]));
    k_r = k_f / Kc;
    q_r[10] = phi_r * k_r;

    /*reaction 12: h + ho2 <=> o + h2o */
    phi_f = sc[1]*sc[13];
    k_f = 1e-06 * 3.01e+13*exp(-866.03588184992258903/tc[1]);
    q_f[11] = phi_f * k_f;
    phi_r = sc[11]*sc[15];
    Kc = exp((g_RT[1] + g_RT[13]) - (g_RT[11] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[11] = phi_r * k_r;

    /*reaction 13: o + ho2 <=> o2 + oh */
    phi_f = sc[11]*sc[13];
    k_f = 1e-06 * 3.25e+13;
    q_f[12] = phi_f * k_f;
    phi_r = sc[10]*sc[12];
    Kc = exp((g_RT[11] + g_RT[13]) - (g_RT[10] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[12] = phi_r * k_r;

    /*reaction 14: oh + oh <=> o + h2o */
    phi_f = sc[12]*sc[12];
    k_f = 1e-06 * 35700*exp(2.4*tc[0]+1062.7935981795678799/tc[1]);
    q_f[13] = phi_f * k_f;
    phi_r = sc[11]*sc[15];
    Kc = exp((g_RT[12] + g_RT[12]) - (g_RT[11] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[13] = phi_r * k_r;

    /*reaction 15: h + h + M <=> h2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[15] + -1*sc[0];
    k_f = 1e-12 * alpha * 1e+18*exp(-1*tc[0]);
    q_f[14] = phi_f * k_f;
    phi_r = sc[0];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[1]) - (g_RT[0]));
    k_r = k_f / Kc;
    q_r[14] = phi_r * k_r;

    /*reaction 16: h + h + h2 <=> h2 + h2 */
    phi_f = sc[1]*sc[1]*sc[0];
    k_f = 1e-12 * 9.2e+16*exp(-0.6*tc[0]);
    q_f[15] = phi_f * k_f;
    phi_r = sc[0]*sc[0];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[1] + g_RT[0]) - (g_RT[0] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[15] = phi_r * k_r;

    /*reaction 17: h + h + h2o <=> h2 + h2o */
    phi_f = sc[1]*sc[1]*sc[15];
    k_f = 1e-12 * 6e+19*exp(-1.25*tc[0]);
    q_f[16] = phi_f * k_f;
    phi_r = sc[0]*sc[15];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[1] + g_RT[15]) - (g_RT[0] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[16] = phi_r * k_r;

    /*reaction 18: h + oh + M <=> h2o + M */
    phi_f = sc[1]*sc[12];
    alpha = mixture + 5.4*sc[15];
    k_f = 1e-12 * alpha * 2.21e+22*exp(-2*tc[0]);
    q_f[17] = phi_f * k_f;
    phi_r = sc[15];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[12]) - (g_RT[15]));
    k_r = k_f / Kc;
    q_r[17] = phi_r * k_r;

    /*reaction 19: h + o + M <=> oh + M */
    phi_f = sc[1]*sc[11];
    alpha = mixture + 5.4*sc[15];
    k_f = 1e-12 * alpha * 4.71e+18*exp(-1*tc[0]);
    q_f[18] = phi_f * k_f;
    phi_r = sc[12];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[11]) - (g_RT[12]));
    k_r = k_f / Kc;
    q_r[18] = phi_r * k_r;

    /*reaction 20: o + o + M <=> o2 + M */
    phi_f = sc[11]*sc[11];
    alpha = mixture;
    k_f = 1e-12 * alpha * 1.89e+13*exp(+899.75139845883882117/tc[1]);
    q_f[19] = phi_f * k_f;
    phi_r = sc[10];
    Kc = 1.0 / (refC) * exp((g_RT[11] + g_RT[11]) - (g_RT[10]));
    k_r = k_f / Kc;
    q_r[19] = phi_r * k_r;

    /*reaction 21: ho2 + ho2 <=> h2o2 + o2 */
    phi_f = sc[13]*sc[13];
    k_f = 1e-06 * 4.2e+14*exp(-6029.5420896721516328/tc[1]);
    q_f[20] = phi_f * k_f;
    phi_r = sc[14]*sc[10];
    Kc = exp((g_RT[13] + g_RT[13]) - (g_RT[14] + g_RT[10]));
    k_r = k_f / Kc;
    q_r[20] = phi_r * k_r;

    /*reaction 22: ho2 + ho2 <=> h2o2 + o2 */
    phi_f = sc[13]*sc[13];
    k_f = 1e-06 * 1.3e+11*exp(+819.73994859588844974/tc[1]);
    q_f[21] = phi_f * k_f;
    phi_r = sc[14]*sc[10];
    Kc = exp((g_RT[13] + g_RT[13]) - (g_RT[14] + g_RT[10]));
    k_r = k_f / Kc;
    q_r[21] = phi_r * k_r;

    /*reaction 23: oh + oh (+M) <=> h2o2 (+M) */
    phi_f = sc[12]*sc[12];
    alpha = mixture;
    k_f = 1e-06 * 1.24e+14*exp(-0.37*tc[0]);
    redP = 1e-12 * alpha / k_f * 3.04e+30*exp(-4.63*tc[0]-1031.0909482338706766/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.53*exp(T/-100))+ (0.47*exp(T/-2000))+ (exp(-1e+15/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[22] = phi_f * k_f;
    phi_r = sc[14];
    Kc = 1.0 / (refC) * exp((g_RT[12] + g_RT[12]) - (g_RT[14]));
    k_r = k_f / Kc;
    q_r[22] = phi_r * k_r;

    /*reaction 24: h2o2 + h <=> ho2 + h2 */
    phi_f = sc[14]*sc[1];
    k_f = 1e-06 * 1.98e+06*exp(2*tc[0]-1225.3325812344924088/tc[1]);
    q_f[23] = phi_f * k_f;
    phi_r = sc[13]*sc[0];
    Kc = exp((g_RT[14] + g_RT[1]) - (g_RT[13] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[23] = phi_r * k_r;

    /*reaction 25: h2o2 + h <=> oh + h2o */
    phi_f = sc[14]*sc[1];
    k_f = 1e-06 * 3.07e+13*exp(-2122.0646796985029141/tc[1]);
    q_f[24] = phi_f * k_f;
    phi_r = sc[12]*sc[15];
    Kc = exp((g_RT[14] + g_RT[1]) - (g_RT[12] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[24] = phi_r * k_r;

    /*reaction 26: h2o2 + o <=> oh + ho2 */
    phi_f = sc[14]*sc[11];
    k_f = 1e-06 * 9.55e+06*exp(2*tc[0]-1997.7701632447369775/tc[1]);
    q_f[25] = phi_f * k_f;
    phi_r = sc[12]*sc[13];
    Kc = exp((g_RT[14] + g_RT[11]) - (g_RT[12] + g_RT[13]));
    k_r = k_f / Kc;
    q_r[25] = phi_r * k_r;

    /*reaction 27: h2o2 + oh <=> h2o + ho2 */
    phi_f = sc[14]*sc[12];
    k_f = 1e-06 * 2.4*exp(4.042*tc[0]+1087.9544314698041489/tc[1]);
    q_f[26] = phi_f * k_f;
    phi_r = sc[15]*sc[13];
    Kc = exp((g_RT[14] + g_RT[12]) - (g_RT[15] + g_RT[13]));
    k_r = k_f / Kc;
    q_r[26] = phi_r * k_r;

    /*reaction 28: ch3 + ch3 (+M) <=> c2h6 (+M) */
    phi_f = sc[3]*sc[3];
    alpha = mixture + 4*sc[15] + sc[0] + 2*sc[8] + sc[9];
    k_f = 1e-06 * 9.22e+16*exp(-1.174*tc[0]-320.04579945180171308/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.14e+36*exp(-5.246*tc[0]-857.98441519704704206/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.595*exp(T/-1120))+ (0.405*exp(T/-69.6))+ (exp(-1e+15/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[27] = phi_f * k_f;
    phi_r = sc[22];
    Kc = 1.0 / (refC) * exp((g_RT[3] + g_RT[3]) - (g_RT[22]));
    k_r = k_f / Kc;
    q_r[27] = phi_r * k_r;

    /*reaction 29: ch3 + h (+M) <=> ch4 (+M) */
    phi_f = sc[3]*sc[1];
    alpha = mixture + 4*sc[15] + sc[0] + 2*sc[8] + sc[9];
    k_f = 1e-06 * 2.14e+15*exp(-0.4*tc[0]);
    redP = 1e-12 * alpha / k_f * 3.31e+30*exp(-4*tc[0]-1060.7807315163490784/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((1*exp(T/-1e-15))+ (0*exp(T/-1e-15))+ (exp(-40/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[28] = phi_f * k_f;
    phi_r = sc[2];
    Kc = 1.0 / (refC) * exp((g_RT[3] + g_RT[1]) - (g_RT[2]));
    k_r = k_f / Kc;
    q_r[28] = phi_r * k_r;

    /*reaction 30: ch4 + h <=> ch3 + h2 */
    phi_f = sc[2]*sc[1];
    k_f = 1e-06 * 22000*exp(3*tc[0]-4403.145825791297284/tc[1]);
    q_f[29] = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = exp((g_RT[2] + g_RT[1]) - (g_RT[3] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[29] = phi_r * k_r;

    /*reaction 31: ch4 + oh <=> ch3 + h2o */
    phi_f = sc[2]*sc[12];
    k_f = 1e-06 * 4.19e+06*exp(2*tc[0]-1281.6928478046211239/tc[1]);
    q_f[30] = phi_f * k_f;
    phi_r = sc[3]*sc[15];
    Kc = exp((g_RT[2] + g_RT[12]) - (g_RT[3] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[30] = phi_r * k_r;

    /*reaction 32: ch4 + o <=> ch3 + oh */
    phi_f = sc[2]*sc[11];
    k_f = 1e-06 * 6.92e+08*exp(1.56*tc[0]-4269.7934093530466271/tc[1]);
    q_f[31] = phi_f * k_f;
    phi_r = sc[3]*sc[12];
    Kc = exp((g_RT[2] + g_RT[11]) - (g_RT[3] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[31] = phi_r * k_r;

    /*reaction 33: ch4 + ho2 <=> ch3 + h2o2 */
    phi_f = sc[2]*sc[13];
    k_f = 1e-06 * 1.12e+13*exp(-12399.258645428293676/tc[1]);
    q_f[32] = phi_f * k_f;
    phi_r = sc[3]*sc[14];
    Kc = exp((g_RT[2] + g_RT[13]) - (g_RT[3] + g_RT[14]));
    k_r = k_f / Kc;
    q_r[32] = phi_r * k_r;

    /*reaction 34: ch3 + ho2 <=> ch3o + oh */
    phi_f = sc[3]*sc[13];
    k_f = 1e-06 * 7e+12;
    q_f[33] = phi_f * k_f;
    phi_r = sc[24]*sc[12];
    Kc = exp((g_RT[3] + g_RT[13]) - (g_RT[24] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[33] = phi_r * k_r;

    /*reaction 35: ch3 + ho2 <=> ch4 + o2 */
    phi_f = sc[3]*sc[13];
    k_f = 1e-06 * 3e+12;
    q_f[34] = phi_f * k_f;
    phi_r = sc[2]*sc[10];
    Kc = exp((g_RT[3] + g_RT[13]) - (g_RT[2] + g_RT[10]));
    k_r = k_f / Kc;
    q_r[34] = phi_r * k_r;

    /*reaction 36: ch3 + o <=> ch2o + h */
    phi_f = sc[3]*sc[11];
    k_f = 1e-06 * 8e+13;
    q_f[35] = phi_f * k_f;
    phi_r = sc[6]*sc[1];
    Kc = exp((g_RT[3] + g_RT[11]) - (g_RT[6] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[35] = phi_r * k_r;

    /*reaction 37: ch3 + o2 <=> ch3o + o */
    phi_f = sc[3]*sc[10];
    k_f = 1e-06 * 1.45e+13*exp(-14698.455591490057486/tc[1]);
    q_f[36] = phi_f * k_f;
    phi_r = sc[24]*sc[11];
    Kc = exp((g_RT[3] + g_RT[10]) - (g_RT[24] + g_RT[11]));
    k_r = k_f / Kc;
    q_r[36] = phi_r * k_r;

    /*reaction 38: ch3 + o2 <=> ch2o + oh */
    phi_f = sc[3]*sc[10];
    k_f = 1e-06 * 2.51e+11*exp(-7367.0919873810962599/tc[1]);
    q_f[37] = phi_f * k_f;
    phi_r = sc[6]*sc[12];
    Kc = exp((g_RT[3] + g_RT[10]) - (g_RT[6] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[37] = phi_r * k_r;

    /*reaction 39: ch3o + h <=> ch3 + oh */
    phi_f = sc[24]*sc[1];
    k_f = 1e-06 * 1e+13;
    q_f[38] = phi_f * k_f;
    phi_r = sc[3]*sc[12];
    Kc = exp((g_RT[24] + g_RT[1]) - (g_RT[3] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[38] = phi_r * k_r;

    /*reaction 40: ch2oh + h <=> ch3 + oh */
    phi_f = sc[23]*sc[1];
    k_f = 1e-06 * 1e+13;
    q_f[39] = phi_f * k_f;
    phi_r = sc[3]*sc[12];
    Kc = exp((g_RT[23] + g_RT[1]) - (g_RT[3] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[39] = phi_r * k_r;

    /*reaction 41: ch3 + oh <=> ch2s + h2o */
    phi_f = sc[3]*sc[12];
    k_f = 1e-06 * 2e+13*exp(-276.76916619259583285/tc[1]);
    q_f[40] = phi_f * k_f;
    phi_r = sc[28]*sc[15];
    Kc = exp((g_RT[3] + g_RT[12]) - (g_RT[28] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[40] = phi_r * k_r;

    /*reaction 42: ch3 + oh <=> ch2 + h2o */
    phi_f = sc[3]*sc[12];
    k_f = 1e-06 * 3e+06*exp(2*tc[0]-1258.0416645117993539/tc[1]);
    q_f[41] = phi_f * k_f;
    phi_r = sc[4]*sc[15];
    Kc = exp((g_RT[3] + g_RT[12]) - (g_RT[4] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[41] = phi_r * k_r;

    /*reaction 43: ch3 + h <=> ch2 + h2 */
    phi_f = sc[3]*sc[1];
    k_f = 1e-06 * 9e+13*exp(-7598.5716536512672974/tc[1]);
    q_f[42] = phi_f * k_f;
    phi_r = sc[4]*sc[0];
    Kc = exp((g_RT[3] + g_RT[1]) - (g_RT[4] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[42] = phi_r * k_r;

    /*reaction 44: ch3 + M <=> ch + h2 + M */
    phi_f = sc[3];
    alpha = mixture;
    k_f = 1e-06 * alpha * 6.9e+14*exp(-41499.775212249434844/tc[1]);
    q_f[43] = phi_f * k_f;
    phi_r = sc[5]*sc[0];
    Kc = refC * exp((g_RT[3]) - (g_RT[5] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[43] = phi_r * k_r;

    /*reaction 45: ch3 + M <=> ch2 + h + M */
    phi_f = sc[3];
    alpha = mixture;
    k_f = 1e-06 * alpha * 1.9e+16*exp(-45999.538637875230052/tc[1]);
    q_f[44] = phi_f * k_f;
    phi_r = sc[4]*sc[1];
    Kc = refC * exp((g_RT[3]) - (g_RT[4] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[44] = phi_r * k_r;

    /*reaction 46: ch2o + h (+M) <=> ch3o (+M) */
    phi_f = sc[6]*sc[1];
    alpha = mixture + 4*sc[15];
    k_f = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1308.3633310922712099/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.5e+30*exp(-4.8*tc[0]-2797.8846618742413739/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.242*exp(T/-94))+ (0.758*exp(T/-1555))+ (exp(-4200/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[45] = phi_f * k_f;
    phi_r = sc[24];
    Kc = 1.0 / (refC) * exp((g_RT[6] + g_RT[1]) - (g_RT[24]));
    k_r = k_f / Kc;
    q_r[45] = phi_r * k_r;

    /*reaction 47: ch2o + h (+M) <=> ch2oh (+M) */
    phi_f = sc[6]*sc[1];
    alpha = mixture + 4*sc[15];
    k_f = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1811.579996896990906/tc[1]);
    redP = 1e-12 * alpha / k_f * 9.1e+31*exp(-4.82*tc[0]-3286.0048277048194905/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2813*exp(T/-103))+ (0.7187*exp(T/-1291))+ (exp(-4160/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[46] = phi_f * k_f;
    phi_r = sc[23];
    Kc = 1.0 / (refC) * exp((g_RT[6] + g_RT[1]) - (g_RT[23]));
    k_r = k_f / Kc;
    q_r[46] = phi_r * k_r;

    /*reaction 48: ch3o + ch3 <=> ch2o + ch4 */
    phi_f = sc[24]*sc[3];
    k_f = 1e-06 * 1.2e+13;
    q_f[47] = phi_f * k_f;
    phi_r = sc[6]*sc[2];
    Kc = exp((g_RT[24] + g_RT[3]) - (g_RT[6] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[47] = phi_r * k_r;

    /*reaction 49: ch3o + h <=> ch2o + h2 */
    phi_f = sc[24]*sc[1];
    k_f = 1e-06 * 2e+13;
    q_f[48] = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = exp((g_RT[24] + g_RT[1]) - (g_RT[6] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[48] = phi_r * k_r;

    /*reaction 50: ch2oh + h <=> ch2o + h2 */
    phi_f = sc[23]*sc[1];
    k_f = 1e-06 * 2e+13;
    q_f[49] = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = exp((g_RT[23] + g_RT[1]) - (g_RT[6] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[49] = phi_r * k_r;

    /*reaction 51: ch3o + oh <=> ch2o + h2o */
    phi_f = sc[24]*sc[12];
    k_f = 1e-06 * 1e+13;
    q_f[50] = phi_f * k_f;
    phi_r = sc[6]*sc[15];
    Kc = exp((g_RT[24] + g_RT[12]) - (g_RT[6] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[50] = phi_r * k_r;

    /*reaction 52: ch2oh + oh <=> ch2o + h2o */
    phi_f = sc[23]*sc[12];
    k_f = 1e-06 * 1e+13;
    q_f[51] = phi_f * k_f;
    phi_r = sc[6]*sc[15];
    Kc = exp((g_RT[23] + g_RT[12]) - (g_RT[6] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[51] = phi_r * k_r;

    /*reaction 53: ch3o + o <=> ch2o + oh */
    phi_f = sc[24]*sc[11];
    k_f = 1e-06 * 1e+13;
    q_f[52] = phi_f * k_f;
    phi_r = sc[6]*sc[12];
    Kc = exp((g_RT[24] + g_RT[11]) - (g_RT[6] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[52] = phi_r * k_r;

    /*reaction 54: ch2oh + o <=> ch2o + oh */
    phi_f = sc[23]*sc[11];
    k_f = 1e-06 * 1e+13;
    q_f[53] = phi_f * k_f;
    phi_r = sc[6]*sc[12];
    Kc = exp((g_RT[23] + g_RT[11]) - (g_RT[6] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[53] = phi_r * k_r;

    /*reaction 55: ch3o + o2 <=> ch2o + ho2 */
    phi_f = sc[24]*sc[10];
    k_f = 1e-06 * 6.3e+10*exp(-1308.3633310922712099/tc[1]);
    q_f[54] = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = exp((g_RT[24] + g_RT[10]) - (g_RT[6] + g_RT[13]));
    k_r = k_f / Kc;
    q_r[54] = phi_r * k_r;

    /*reaction 56: ch3o + co <=> ch3 + co2 */
    phi_f = sc[24]*sc[9];
    k_f = 1e-06 * 468*exp(3.16*tc[0]-2707.3056620293918968/tc[1]);
    q_f[55] = phi_f * k_f;
    phi_r = sc[3]*sc[8];
    Kc = exp((g_RT[24] + g_RT[9]) - (g_RT[3] + g_RT[8]));
    k_r = k_f / Kc;
    q_r[55] = phi_r * k_r;

    /*reaction 57: ch2oh + o2 <=> ch2o + ho2 */
    phi_f = sc[23]*sc[10];
    k_f = 1e-06 * 1.57e+15*exp(-1*tc[0]);
    q_f[56] = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = exp((g_RT[23] + g_RT[10]) - (g_RT[6] + g_RT[13]));
    k_r = k_f / Kc;
    q_r[56] = phi_r * k_r;

    /*reaction 58: ch2oh + o2 <=> ch2o + ho2 */
    phi_f = sc[23]*sc[10];
    k_f = 1e-06 * 7.23e+13*exp(-1800.0060135834821722/tc[1]);
    q_f[57] = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = exp((g_RT[23] + g_RT[10]) - (g_RT[6] + g_RT[13]));
    k_r = k_f / Kc;
    q_r[57] = phi_r * k_r;

    /*reaction 59: ch2 + h <=> ch + h2 */
    phi_f = sc[4]*sc[1];
    k_f = 1e-06 * 1e+18*exp(-1.56*tc[0]);
    q_f[58] = phi_f * k_f;
    phi_r = sc[5]*sc[0];
    Kc = exp((g_RT[4] + g_RT[1]) - (g_RT[5] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[58] = phi_r * k_r;

    /*reaction 60: ch2 + oh <=> ch + h2o */
    phi_f = sc[4]*sc[12];
    k_f = 1e-06 * 1.13e+07*exp(2*tc[0]-1509.6499974141590883/tc[1]);
    q_f[59] = phi_f * k_f;
    phi_r = sc[5]*sc[15];
    Kc = exp((g_RT[4] + g_RT[12]) - (g_RT[5] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[59] = phi_r * k_r;

    /*reaction 61: ch2 + oh <=> ch2o + h */
    phi_f = sc[4]*sc[12];
    k_f = 1e-06 * 2.5e+13;
    q_f[60] = phi_f * k_f;
    phi_r = sc[6]*sc[1];
    Kc = exp((g_RT[4] + g_RT[12]) - (g_RT[6] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[60] = phi_r * k_r;

    /*reaction 62: ch2 + co2 <=> ch2o + co */
    phi_f = sc[4]*sc[8];
    k_f = 1e-06 * 1.1e+11*exp(-503.21666580471963925/tc[1]);
    q_f[61] = phi_f * k_f;
    phi_r = sc[6]*sc[9];
    Kc = exp((g_RT[4] + g_RT[8]) - (g_RT[6] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[61] = phi_r * k_r;

    /*reaction 63: ch2 + o <=> co + h + h */
    phi_f = sc[4]*sc[11];
    k_f = 1e-06 * 5e+13;
    q_f[62] = phi_f * k_f;
    phi_r = sc[9]*sc[1]*sc[1];
    Kc = refC * exp((g_RT[4] + g_RT[11]) - (g_RT[9] + g_RT[1] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[62] = phi_r * k_r;

    /*reaction 64: ch2 + o <=> co + h2 */
    phi_f = sc[4]*sc[11];
    k_f = 1e-06 * 3e+13;
    q_f[63] = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = exp((g_RT[4] + g_RT[11]) - (g_RT[9] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[63] = phi_r * k_r;

    /*reaction 65: ch2 + o2 <=> ch2o + o */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 3.29e+21*exp(-3.3*tc[0]-1443.2253975279359111/tc[1]);
    q_f[64] = phi_f * k_f;
    phi_r = sc[6]*sc[11];
    Kc = exp((g_RT[4] + g_RT[10]) - (g_RT[6] + g_RT[11]));
    k_r = k_f / Kc;
    q_r[64] = phi_r * k_r;

    /*reaction 66: ch2 + o2 <=> co2 + h + h */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 3.29e+21*exp(-3.3*tc[0]-1443.2253975279359111/tc[1]);
    q_f[65] = phi_f * k_f;
    phi_r = sc[8]*sc[1]*sc[1];
    Kc = refC * exp((g_RT[4] + g_RT[10]) - (g_RT[8] + g_RT[1] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[65] = phi_r * k_r;

    /*reaction 67: ch2 + o2 <=> co2 + h2 */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 1.01e+21*exp(-3.3*tc[0]-758.85073203351726079/tc[1]);
    q_f[66] = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = exp((g_RT[4] + g_RT[10]) - (g_RT[8] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[66] = phi_r * k_r;

    /*reaction 68: ch2 + o2 <=> co + h2o */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 7.28e+19*exp(-2.54*tc[0]-910.31894844073792683/tc[1]);
    q_f[67] = phi_f * k_f;
    phi_r = sc[9]*sc[15];
    Kc = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[67] = phi_r * k_r;

    /*reaction 69: ch2 + o2 <=> hco + oh */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 1.29e+20*exp(-3.3*tc[0]-142.9135330885404187/tc[1]);
    q_f[68] = phi_f * k_f;
    phi_r = sc[7]*sc[12];
    Kc = exp((g_RT[4] + g_RT[10]) - (g_RT[7] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[68] = phi_r * k_r;

    /*reaction 70: ch2 + ch3 <=> c2h4 + h */
    phi_f = sc[4]*sc[3];
    k_f = 1e-06 * 4e+13;
    q_f[69] = phi_f * k_f;
    phi_r = sc[20]*sc[1];
    Kc = exp((g_RT[4] + g_RT[3]) - (g_RT[20] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[69] = phi_r * k_r;

    /*reaction 71: ch2 + ch2 <=> c2h2 + h + h */
    phi_f = sc[4]*sc[4];
    k_f = 1e-06 * 4e+13;
    q_f[70] = phi_f * k_f;
    phi_r = sc[18]*sc[1]*sc[1];
    Kc = refC * exp((g_RT[4] + g_RT[4]) - (g_RT[18] + g_RT[1] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[70] = phi_r * k_r;

    /*reaction 72: ch2 + hcco <=> c2h3 + co */
    phi_f = sc[4]*sc[17];
    k_f = 1e-06 * 3e+13;
    q_f[71] = phi_f * k_f;
    phi_r = sc[19]*sc[9];
    Kc = exp((g_RT[4] + g_RT[17]) - (g_RT[19] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[71] = phi_r * k_r;

    /*reaction 73: ch2 + c2h2 <=> h2ccch + h */
    phi_f = sc[4]*sc[18];
    k_f = 1e-06 * 1.2e+13*exp(-3321.2299943111497669/tc[1]);
    q_f[72] = phi_f * k_f;
    phi_r = sc[26]*sc[1];
    Kc = exp((g_RT[4] + g_RT[18]) - (g_RT[26] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[72] = phi_r * k_r;

    /*reaction 74: ch2s + M <=> ch2 + M */
    phi_f = sc[28];
    alpha = mixture + 11*sc[1] + 3*sc[18] + 2*sc[15];
    k_f = 1e-06 * alpha * 1e+13;
    q_f[73] = phi_f * k_f;
    phi_r = sc[4];
    Kc = exp((g_RT[28]) - (g_RT[4]));
    k_r = k_f / Kc;
    q_r[73] = phi_r * k_r;

    /*reaction 75: ch2s + ch4 <=> ch3 + ch3 */
    phi_f = sc[28]*sc[2];
    k_f = 1e-06 * 4e+13;
    q_f[74] = phi_f * k_f;
    phi_r = sc[3]*sc[3];
    Kc = exp((g_RT[28] + g_RT[2]) - (g_RT[3] + g_RT[3]));
    k_r = k_f / Kc;
    q_r[74] = phi_r * k_r;

    /*reaction 76: ch2s + c2h6 <=> ch3 + c2h5 */
    phi_f = sc[28]*sc[22];
    k_f = 1e-06 * 1.2e+14;
    q_f[75] = phi_f * k_f;
    phi_r = sc[3]*sc[21];
    Kc = exp((g_RT[28] + g_RT[22]) - (g_RT[3] + g_RT[21]));
    k_r = k_f / Kc;
    q_r[75] = phi_r * k_r;

    /*reaction 77: ch2s + o2 <=> co + oh + h */
    phi_f = sc[28]*sc[10];
    k_f = 1e-06 * 7e+13;
    q_f[76] = phi_f * k_f;
    phi_r = sc[9]*sc[12]*sc[1];
    Kc = refC * exp((g_RT[28] + g_RT[10]) - (g_RT[9] + g_RT[12] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[76] = phi_r * k_r;

    /*reaction 78: ch2s + h2 <=> ch3 + h */
    phi_f = sc[28]*sc[0];
    k_f = 1e-06 * 7e+13;
    q_f[77] = phi_f * k_f;
    phi_r = sc[3]*sc[1];
    Kc = exp((g_RT[28] + g_RT[0]) - (g_RT[3] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[77] = phi_r * k_r;

    /*reaction 79: ch2s + c2h2 <=> h2ccch + h */
    phi_f = sc[28]*sc[18];
    k_f = 1e-06 * 1.5e+14;
    q_f[78] = phi_f * k_f;
    phi_r = sc[26]*sc[1];
    Kc = exp((g_RT[28] + g_RT[18]) - (g_RT[26] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[78] = phi_r * k_r;

    /*reaction 80: ch2s + o <=> co + h + h */
    phi_f = sc[28]*sc[11];
    k_f = 1e-06 * 3e+13;
    q_f[79] = phi_f * k_f;
    phi_r = sc[9]*sc[1]*sc[1];
    Kc = refC * exp((g_RT[28] + g_RT[11]) - (g_RT[9] + g_RT[1] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[79] = phi_r * k_r;

    /*reaction 81: ch2s + oh <=> ch2o + h */
    phi_f = sc[28]*sc[12];
    k_f = 1e-06 * 3e+13;
    q_f[80] = phi_f * k_f;
    phi_r = sc[6]*sc[1];
    Kc = exp((g_RT[28] + g_RT[12]) - (g_RT[6] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[80] = phi_r * k_r;

    /*reaction 82: ch2s + h <=> ch + h2 */
    phi_f = sc[28]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f[81] = phi_f * k_f;
    phi_r = sc[5]*sc[0];
    Kc = exp((g_RT[28] + g_RT[1]) - (g_RT[5] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[81] = phi_r * k_r;

    /*reaction 83: ch2s + co2 <=> ch2o + co */
    phi_f = sc[28]*sc[8];
    k_f = 1e-06 * 3e+12;
    q_f[82] = phi_f * k_f;
    phi_r = sc[6]*sc[9];
    Kc = exp((g_RT[28] + g_RT[8]) - (g_RT[6] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[82] = phi_r * k_r;

    /*reaction 84: ch2s + ch3 <=> c2h4 + h */
    phi_f = sc[28]*sc[3];
    k_f = 1e-06 * 2e+13;
    q_f[83] = phi_f * k_f;
    phi_r = sc[20]*sc[1];
    Kc = exp((g_RT[28] + g_RT[3]) - (g_RT[20] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[83] = phi_r * k_r;

    /*reaction 85: ch2s + ch2co <=> c2h4 + co */
    phi_f = sc[28]*sc[29];
    k_f = 1e-06 * 1.6e+14;
    q_f[84] = phi_f * k_f;
    phi_r = sc[20]*sc[9];
    Kc = exp((g_RT[28] + g_RT[29]) - (g_RT[20] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[84] = phi_r * k_r;

    /*reaction 86: ch + o2 <=> hco + o */
    phi_f = sc[5]*sc[10];
    k_f = 1e-06 * 3.3e+13;
    q_f[85] = phi_f * k_f;
    phi_r = sc[7]*sc[11];
    Kc = exp((g_RT[5] + g_RT[10]) - (g_RT[7] + g_RT[11]));
    k_r = k_f / Kc;
    q_r[85] = phi_r * k_r;

    /*reaction 87: ch + o <=> co + h */
    phi_f = sc[5]*sc[11];
    k_f = 1e-06 * 5.7e+13;
    q_f[86] = phi_f * k_f;
    phi_r = sc[9]*sc[1];
    Kc = exp((g_RT[5] + g_RT[11]) - (g_RT[9] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[86] = phi_r * k_r;

    /*reaction 88: ch + oh <=> hco + h */
    phi_f = sc[5]*sc[12];
    k_f = 1e-06 * 3e+13;
    q_f[87] = phi_f * k_f;
    phi_r = sc[7]*sc[1];
    Kc = exp((g_RT[5] + g_RT[12]) - (g_RT[7] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[87] = phi_r * k_r;

    /*reaction 89: ch + co2 <=> hco + co */
    phi_f = sc[5]*sc[8];
    k_f = 1e-06 * 3.4e+12*exp(-347.2194994052565562/tc[1]);
    q_f[88] = phi_f * k_f;
    phi_r = sc[7]*sc[9];
    Kc = exp((g_RT[5] + g_RT[8]) - (g_RT[7] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[88] = phi_r * k_r;

    /*reaction 90: ch + h2o <=> ch2o + h */
    phi_f = sc[5]*sc[15];
    k_f = 1e-06 * 1.17e+15*exp(-0.75*tc[0]);
    q_f[89] = phi_f * k_f;
    phi_r = sc[6]*sc[1];
    Kc = exp((g_RT[5] + g_RT[15]) - (g_RT[6] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[89] = phi_r * k_r;

    /*reaction 91: ch + ch2o <=> ch2co + h */
    phi_f = sc[5]*sc[6];
    k_f = 1e-06 * 9.46e+13*exp(+259.15658288943063781/tc[1]);
    q_f[90] = phi_f * k_f;
    phi_r = sc[29]*sc[1];
    Kc = exp((g_RT[5] + g_RT[6]) - (g_RT[29] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[90] = phi_r * k_r;

    /*reaction 92: ch + c2h2 <=> c3h2 + h */
    phi_f = sc[5]*sc[18];
    k_f = 1e-06 * 1e+14;
    q_f[91] = phi_f * k_f;
    phi_r = sc[27]*sc[1];
    Kc = exp((g_RT[5] + g_RT[18]) - (g_RT[27] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[91] = phi_r * k_r;

    /*reaction 93: ch + ch2 <=> c2h2 + h */
    phi_f = sc[5]*sc[4];
    k_f = 1e-06 * 4e+13;
    q_f[92] = phi_f * k_f;
    phi_r = sc[18]*sc[1];
    Kc = exp((g_RT[5] + g_RT[4]) - (g_RT[18] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[92] = phi_r * k_r;

    /*reaction 94: ch + ch3 <=> c2h3 + h */
    phi_f = sc[5]*sc[3];
    k_f = 1e-06 * 3e+13;
    q_f[93] = phi_f * k_f;
    phi_r = sc[19]*sc[1];
    Kc = exp((g_RT[5] + g_RT[3]) - (g_RT[19] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[93] = phi_r * k_r;

    /*reaction 95: ch + ch4 <=> c2h4 + h */
    phi_f = sc[5]*sc[2];
    k_f = 1e-06 * 6e+13;
    q_f[94] = phi_f * k_f;
    phi_r = sc[20]*sc[1];
    Kc = exp((g_RT[5] + g_RT[2]) - (g_RT[20] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[94] = phi_r * k_r;

    /*reaction 96: ch2o + oh <=> hco + h2o */
    phi_f = sc[6]*sc[12];
    k_f = 1e-06 * 3.43e+09*exp(1.18*tc[0]+224.93784961470970529/tc[1]);
    q_f[95] = phi_f * k_f;
    phi_r = sc[7]*sc[15];
    Kc = exp((g_RT[6] + g_RT[12]) - (g_RT[7] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[95] = phi_r * k_r;

    /*reaction 97: ch2o + h <=> hco + h2 */
    phi_f = sc[6]*sc[1];
    k_f = 1e-06 * 2.19e+08*exp(1.77*tc[0]-1509.6499974141590883/tc[1]);
    q_f[96] = phi_f * k_f;
    phi_r = sc[7]*sc[0];
    Kc = exp((g_RT[6] + g_RT[1]) - (g_RT[7] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[96] = phi_r * k_r;

    /*reaction 98: ch2o + M <=> hco + h + M */
    phi_f = sc[6];
    alpha = mixture;
    k_f = 1e-06 * alpha * 3.31e+16*exp(-40760.549930182300159/tc[1]);
    q_f[97] = phi_f * k_f;
    phi_r = sc[7]*sc[1];
    Kc = refC * exp((g_RT[6]) - (g_RT[7] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[97] = phi_r * k_r;

    /*reaction 99: ch2o + o <=> hco + oh */
    phi_f = sc[6]*sc[11];
    k_f = 1e-06 * 1.8e+13*exp(-1549.9073306785367095/tc[1]);
    q_f[98] = phi_f * k_f;
    phi_r = sc[7]*sc[12];
    Kc = exp((g_RT[6] + g_RT[11]) - (g_RT[7] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[98] = phi_r * k_r;

    /*reaction 100: hco + o2 <=> co + ho2 */
    phi_f = sc[7]*sc[10];
    k_f = 1e-06 * 7.58e+12*exp(-206.31883297993508108/tc[1]);
    q_f[99] = phi_f * k_f;
    phi_r = sc[9]*sc[13];
    Kc = exp((g_RT[7] + g_RT[10]) - (g_RT[9] + g_RT[13]));
    k_r = k_f / Kc;
    q_r[99] = phi_r * k_r;

    /*reaction 101: hco + M <=> h + co + M */
    phi_f = sc[7];
    alpha = mixture + 4*sc[15] + 0.87*sc[0] + 2*sc[8] + 0.87*sc[9] + 1.81*sc[2];
    k_f = 1e-06 * alpha * 1.86e+17*exp(-1*tc[0]-8554.6833186802341515/tc[1]);
    q_f[100] = phi_f * k_f;
    phi_r = sc[1]*sc[9];
    Kc = refC * exp((g_RT[7]) - (g_RT[1] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[100] = phi_r * k_r;

    /*reaction 102: hco + oh <=> h2o + co */
    phi_f = sc[7]*sc[12];
    k_f = 1e-06 * 1e+14;
    q_f[101] = phi_f * k_f;
    phi_r = sc[15]*sc[9];
    Kc = exp((g_RT[7] + g_RT[12]) - (g_RT[15] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[101] = phi_r * k_r;

    /*reaction 103: hco + h <=> co + h2 */
    phi_f = sc[7]*sc[1];
    k_f = 1e-06 * 1.19e+13*exp(0.25*tc[0]);
    q_f[102] = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = exp((g_RT[7] + g_RT[1]) - (g_RT[9] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[102] = phi_r * k_r;

    /*reaction 104: hco + o <=> co + oh */
    phi_f = sc[7]*sc[11];
    k_f = 1e-06 * 3e+13;
    q_f[103] = phi_f * k_f;
    phi_r = sc[9]*sc[12];
    Kc = exp((g_RT[7] + g_RT[11]) - (g_RT[9] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[103] = phi_r * k_r;

    /*reaction 105: hco + o <=> co2 + h */
    phi_f = sc[7]*sc[11];
    k_f = 1e-06 * 3e+13;
    q_f[104] = phi_f * k_f;
    phi_r = sc[8]*sc[1];
    Kc = exp((g_RT[7] + g_RT[11]) - (g_RT[8] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[104] = phi_r * k_r;

    /*reaction 106: co + oh <=> co2 + h */
    phi_f = sc[9]*sc[12];
    k_f = 1e-06 * 9420*exp(2.25*tc[0]+1183.0623813068959862/tc[1]);
    q_f[105] = phi_f * k_f;
    phi_r = sc[8]*sc[1];
    Kc = exp((g_RT[9] + g_RT[12]) - (g_RT[8] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[105] = phi_r * k_r;

    /*reaction 107: co + o + M <=> co2 + M */
    phi_f = sc[9]*sc[11];
    alpha = mixture;
    k_f = 1e-12 * alpha * 6.17e+14*exp(-1509.6499974141590883/tc[1]);
    q_f[106] = phi_f * k_f;
    phi_r = sc[8];
    Kc = 1.0 / (refC) * exp((g_RT[9] + g_RT[11]) - (g_RT[8]));
    k_r = k_f / Kc;
    q_r[106] = phi_r * k_r;

    /*reaction 108: co + o2 <=> co2 + o */
    phi_f = sc[9]*sc[10];
    k_f = 1e-06 * 2.53e+12*exp(-23997.396358895472076/tc[1]);
    q_f[107] = phi_f * k_f;
    phi_r = sc[8]*sc[11];
    Kc = exp((g_RT[9] + g_RT[10]) - (g_RT[8] + g_RT[11]));
    k_r = k_f / Kc;
    q_r[107] = phi_r * k_r;

    /*reaction 109: co + ho2 <=> co2 + oh */
    phi_f = sc[9]*sc[13];
    k_f = 1e-06 * 5.8e+13*exp(-11540.771013565441535/tc[1]);
    q_f[108] = phi_f * k_f;
    phi_r = sc[8]*sc[12];
    Kc = exp((g_RT[9] + g_RT[13]) - (g_RT[8] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[108] = phi_r * k_r;

    /*reaction 110: c2h5oh (+M) <=> ch3 + ch2oh (+M) */
    phi_f = sc[33];
    alpha = mixture + 4*sc[15] + sc[0] + 2*sc[8] + sc[9];
    k_f = 1 * 5.94e+23*exp(-1.68*tc[0]-45874.740904755664815/tc[1]);
    redP = 1e-6 * alpha / k_f * 2.88e+85*exp(-18.9*tc[0]-55310.556605259960634/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-200))+ (0.5*exp(T/-890))+ (exp(-4600/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[109] = phi_f * k_f;
    phi_r = sc[3]*sc[23];
    Kc = refC * exp((g_RT[33]) - (g_RT[3] + g_RT[23]));
    k_r = k_f / Kc;
    q_r[109] = phi_r * k_r;

    /*reaction 111: c2h5oh (+M) <=> c2h5 + oh (+M) */
    phi_f = sc[33];
    alpha = mixture + 4*sc[15] + sc[0] + 2*sc[8] + sc[9];
    k_f = 1 * 1.25e+23*exp(-1.54*tc[0]-48311.316000582111883/tc[1]);
    redP = 1e-6 * alpha / k_f * 3.252e+85*exp(-18.81*tc[0]-57834.691400936433638/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.5*exp(T/-300))+ (0.5*exp(T/-900))+ (exp(-5000/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[110] = phi_f * k_f;
    phi_r = sc[21]*sc[12];
    Kc = refC * exp((g_RT[33]) - (g_RT[21] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[110] = phi_r * k_r;

    /*reaction 112: c2h5oh (+M) <=> c2h4 + h2o (+M) */
    phi_f = sc[33];
    alpha = mixture + 4*sc[15];
    k_f = 1 * 2.79e+13*exp(0.09*tc[0]-33280.737409660941921/tc[1]);
    redP = 1e-6 * alpha / k_f * 2.57e+83*exp(-18.85*tc[0]-43504.087192149629118/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.3*exp(T/-350))+ (0.7*exp(T/-800))+ (exp(-3800/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[111] = phi_f * k_f;
    phi_r = sc[20]*sc[15];
    Kc = refC * exp((g_RT[33]) - (g_RT[20] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[111] = phi_r * k_r;

    /*reaction 113: c2h5oh (+M) <=> ch3hco + h2 (+M) */
    phi_f = sc[33];
    alpha = mixture + 4*sc[15];
    k_f = 1 * 7.24e+11*exp(0.095*tc[0]-45796.239104890126328/tc[1]);
    redP = 1e-6 * alpha / k_f * 4.46e+87*exp(-19.42*tc[0]-58164.801533704332542/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.1*exp(T/-900))+ (0.9*exp(T/-1100))+ (exp(-3500/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[112] = phi_f * k_f;
    phi_r = sc[32]*sc[0];
    Kc = refC * exp((g_RT[33]) - (g_RT[32] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[112] = phi_r * k_r;

    /*reaction 114: c2h5oh + oh <=> c2h4oh + h2o */
    phi_f = sc[33]*sc[12];
    k_f = 1e-06 * 1.74e+11*exp(0.27*tc[0]-301.92999948283181766/tc[1]);
    q_f[113] = phi_f * k_f;
    phi_r = sc[34]*sc[15];
    Kc = exp((g_RT[33] + g_RT[12]) - (g_RT[34] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[113] = phi_r * k_r;

    /*reaction 115: c2h5oh + oh <=> ch3choh + h2o */
    phi_f = sc[33]*sc[12];
    k_f = 1e-06 * 4.64e+11*exp(0.15*tc[0]);
    q_f[114] = phi_f * k_f;
    phi_r = sc[35]*sc[15];
    Kc = exp((g_RT[33] + g_RT[12]) - (g_RT[35] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[114] = phi_r * k_r;

    /*reaction 116: c2h5oh + oh <=> ch3ch2o + h2o */
    phi_f = sc[33]*sc[12];
    k_f = 1e-06 * 7.46e+11*exp(0.3*tc[0]-822.25603192491200844/tc[1]);
    q_f[115] = phi_f * k_f;
    phi_r = sc[36]*sc[15];
    Kc = exp((g_RT[33] + g_RT[12]) - (g_RT[36] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[115] = phi_r * k_r;

    /*reaction 117: c2h5oh + h <=> c2h4oh + h2 */
    phi_f = sc[33]*sc[1];
    k_f = 1e-06 * 1.23e+07*exp(1.8*tc[0]-2565.3985622724612767/tc[1]);
    q_f[116] = phi_f * k_f;
    phi_r = sc[34]*sc[0];
    Kc = exp((g_RT[33] + g_RT[1]) - (g_RT[34] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[116] = phi_r * k_r;

    /*reaction 118: c2h5oh + h <=> ch3choh + h2 */
    phi_f = sc[33]*sc[1];
    k_f = 1e-06 * 2.58e+07*exp(1.65*tc[0]-1422.593514229942457/tc[1]);
    q_f[117] = phi_f * k_f;
    phi_r = sc[35]*sc[0];
    Kc = exp((g_RT[33] + g_RT[1]) - (g_RT[35] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[117] = phi_r * k_r;

    /*reaction 119: c2h5oh + h <=> ch3ch2o + h2 */
    phi_f = sc[33]*sc[1];
    k_f = 1e-06 * 1.5e+07*exp(1.6*tc[0]-1528.7722307147384981/tc[1]);
    q_f[118] = phi_f * k_f;
    phi_r = sc[36]*sc[0];
    Kc = exp((g_RT[33] + g_RT[1]) - (g_RT[36] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[118] = phi_r * k_r;

    /*reaction 120: c2h5oh + o <=> c2h4oh + oh */
    phi_f = sc[33]*sc[11];
    k_f = 1e-06 * 9.41e+07*exp(1.7*tc[0]-2747.0597786279645334/tc[1]);
    q_f[119] = phi_f * k_f;
    phi_r = sc[34]*sc[12];
    Kc = exp((g_RT[33] + g_RT[11]) - (g_RT[34] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[119] = phi_r * k_r;

    /*reaction 121: c2h5oh + o <=> ch3choh + oh */
    phi_f = sc[33]*sc[11];
    k_f = 1e-06 * 1.88e+07*exp(1.85*tc[0]-917.86719842780871659/tc[1]);
    q_f[120] = phi_f * k_f;
    phi_r = sc[35]*sc[12];
    Kc = exp((g_RT[33] + g_RT[11]) - (g_RT[35] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[120] = phi_r * k_r;

    /*reaction 122: c2h5oh + o <=> ch3ch2o + oh */
    phi_f = sc[33]*sc[11];
    k_f = 1e-06 * 1.58e+07*exp(2*tc[0]-2238.3077294993931901/tc[1]);
    q_f[121] = phi_f * k_f;
    phi_r = sc[36]*sc[12];
    Kc = exp((g_RT[33] + g_RT[11]) - (g_RT[36] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[121] = phi_r * k_r;

    /*reaction 123: c2h5oh + ch3 <=> c2h4oh + ch4 */
    phi_f = sc[33]*sc[3];
    k_f = 1e-06 * 219*exp(3.18*tc[0]-4841.9507583730128317/tc[1]);
    q_f[122] = phi_f * k_f;
    phi_r = sc[34]*sc[2];
    Kc = exp((g_RT[33] + g_RT[3]) - (g_RT[34] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[122] = phi_r * k_r;

    /*reaction 124: c2h5oh + ch3 <=> ch3choh + ch4 */
    phi_f = sc[33]*sc[3];
    k_f = 1e-06 * 728*exp(2.99*tc[0]-3999.5660598159120127/tc[1]);
    q_f[123] = phi_f * k_f;
    phi_r = sc[35]*sc[2];
    Kc = exp((g_RT[33] + g_RT[3]) - (g_RT[35] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[123] = phi_r * k_r;

    /*reaction 125: c2h5oh + ch3 <=> ch3ch2o + ch4 */
    phi_f = sc[33]*sc[3];
    k_f = 1e-06 * 145*exp(2.99*tc[0]-3849.1042767403009748/tc[1]);
    q_f[124] = phi_f * k_f;
    phi_r = sc[36]*sc[2];
    Kc = exp((g_RT[33] + g_RT[3]) - (g_RT[36] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[124] = phi_r * k_r;

    /*reaction 126: c2h5oh + ho2 <=> ch3choh + h2o2 */
    phi_f = sc[33]*sc[13];
    k_f = 1e-06 * 8200*exp(2.55*tc[0]-5409.579157400737131/tc[1]);
    q_f[125] = phi_f * k_f;
    phi_r = sc[35]*sc[14];
    Kc = exp((g_RT[33] + g_RT[13]) - (g_RT[35] + g_RT[14]));
    k_r = k_f / Kc;
    q_r[125] = phi_r * k_r;

    /*reaction 127: c2h5oh + ho2 <=> c2h4oh + h2o2 */
    phi_f = sc[33]*sc[13];
    k_f = 1e-06 * 12300*exp(2.55*tc[0]-7925.6624864243349293/tc[1]);
    q_f[126] = phi_f * k_f;
    phi_r = sc[34]*sc[14];
    Kc = exp((g_RT[33] + g_RT[13]) - (g_RT[34] + g_RT[14]));
    k_r = k_f / Kc;
    q_r[126] = phi_r * k_r;

    /*reaction 128: c2h5oh + ho2 <=> ch3ch2o + h2o2 */
    phi_f = sc[33]*sc[13];
    k_f = 1e-06 * 2.5e+12*exp(-12077.199979313272706/tc[1]);
    q_f[127] = phi_f * k_f;
    phi_r = sc[36]*sc[14];
    Kc = exp((g_RT[33] + g_RT[13]) - (g_RT[36] + g_RT[14]));
    k_r = k_f / Kc;
    q_r[127] = phi_r * k_r;

    /*reaction 129: ch3ch2o + M <=> ch3hco + h + M */
    phi_f = sc[36];
    alpha = mixture;
    k_f = 1e-06 * alpha * 1.16e+35*exp(-5.89*tc[0]-12718.298011548486102/tc[1]);
    q_f[128] = phi_f * k_f;
    phi_r = sc[32]*sc[1];
    Kc = refC * exp((g_RT[36]) - (g_RT[32] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[128] = phi_r * k_r;

    /*reaction 130: ch3ch2o + M <=> ch3 + ch2o + M */
    phi_f = sc[36];
    alpha = mixture;
    k_f = 1e-06 * alpha * 1.35e+38*exp(-6.96*tc[0]-11976.55664615232854/tc[1]);
    q_f[129] = phi_f * k_f;
    phi_r = sc[3]*sc[6];
    Kc = refC * exp((g_RT[36]) - (g_RT[3] + g_RT[6]));
    k_r = k_f / Kc;
    q_r[129] = phi_r * k_r;

    /*reaction 131: ch3ch2o + o2 <=> ch3hco + ho2 */
    phi_f = sc[36]*sc[10];
    k_f = 1e-06 * 4e+10*exp(-553.53833238519166571/tc[1]);
    q_f[130] = phi_f * k_f;
    phi_r = sc[32]*sc[13];
    Kc = exp((g_RT[36] + g_RT[10]) - (g_RT[32] + g_RT[13]));
    k_r = k_f / Kc;
    q_r[130] = phi_r * k_r;

    /*reaction 132: ch3ch2o + co <=> c2h5 + co2 */
    phi_f = sc[36]*sc[9];
    k_f = 1e-06 * 468*exp(3.16*tc[0]-2707.3056620293918968/tc[1]);
    q_f[131] = phi_f * k_f;
    phi_r = sc[21]*sc[8];
    Kc = exp((g_RT[36] + g_RT[9]) - (g_RT[21] + g_RT[8]));
    k_r = k_f / Kc;
    q_r[131] = phi_r * k_r;

    /*reaction 133: ch3ch2o + h <=> ch3 + ch2oh */
    phi_f = sc[36]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f[132] = phi_f * k_f;
    phi_r = sc[3]*sc[23];
    Kc = exp((g_RT[36] + g_RT[1]) - (g_RT[3] + g_RT[23]));
    k_r = k_f / Kc;
    q_r[132] = phi_r * k_r;

    /*reaction 134: ch3ch2o + h <=> c2h4 + h2o */
    phi_f = sc[36]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f[133] = phi_f * k_f;
    phi_r = sc[20]*sc[15];
    Kc = exp((g_RT[36] + g_RT[1]) - (g_RT[20] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[133] = phi_r * k_r;

    /*reaction 135: ch3ch2o + oh <=> ch3hco + h2o */
    phi_f = sc[36]*sc[12];
    k_f = 1e-06 * 1e+13;
    q_f[134] = phi_f * k_f;
    phi_r = sc[32]*sc[15];
    Kc = exp((g_RT[36] + g_RT[12]) - (g_RT[32] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[134] = phi_r * k_r;

    /*reaction 136: ch3choh + o2 <=> ch3hco + ho2 */
    phi_f = sc[35]*sc[10];
    k_f = 1e-06 * 4.82e+14*exp(-2524.638012342278671/tc[1]);
    q_f[135] = phi_f * k_f;
    phi_r = sc[32]*sc[13];
    Kc = exp((g_RT[35] + g_RT[10]) - (g_RT[32] + g_RT[13]));
    k_r = k_f / Kc;
    q_r[135] = phi_r * k_r;

    /*reaction 137: ch3choh + o2 <=> ch3hco + ho2 */
    phi_f = sc[35]*sc[10];
    k_f = 1e-06 * 8.43e+15*exp(-1.2*tc[0]);
    q_f[136] = phi_f * k_f;
    phi_r = sc[32]*sc[13];
    Kc = exp((g_RT[35] + g_RT[10]) - (g_RT[32] + g_RT[13]));
    k_r = k_f / Kc;
    q_r[136] = phi_r * k_r;

    /*reaction 138: ch3choh + o <=> ch3hco + oh */
    phi_f = sc[35]*sc[11];
    k_f = 1e-06 * 1e+14;
    q_f[137] = phi_f * k_f;
    phi_r = sc[32]*sc[12];
    Kc = exp((g_RT[35] + g_RT[11]) - (g_RT[32] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[137] = phi_r * k_r;

    /*reaction 139: ch3choh + h <=> c2h4 + h2o */
    phi_f = sc[35]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f[138] = phi_f * k_f;
    phi_r = sc[20]*sc[15];
    Kc = exp((g_RT[35] + g_RT[1]) - (g_RT[20] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[138] = phi_r * k_r;

    /*reaction 140: ch3choh + h <=> ch3 + ch2oh */
    phi_f = sc[35]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f[139] = phi_f * k_f;
    phi_r = sc[3]*sc[23];
    Kc = exp((g_RT[35] + g_RT[1]) - (g_RT[3] + g_RT[23]));
    k_r = k_f / Kc;
    q_r[139] = phi_r * k_r;

    /*reaction 141: ch3choh + ho2 <=> ch3hco + oh + oh */
    phi_f = sc[35]*sc[13];
    k_f = 1e-06 * 4e+13;
    q_f[140] = phi_f * k_f;
    phi_r = sc[32]*sc[12]*sc[12];
    Kc = refC * exp((g_RT[35] + g_RT[13]) - (g_RT[32] + g_RT[12] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[140] = phi_r * k_r;

    /*reaction 142: ch3choh + oh <=> ch3hco + h2o */
    phi_f = sc[35]*sc[12];
    k_f = 1e-06 * 5e+12;
    q_f[141] = phi_f * k_f;
    phi_r = sc[32]*sc[15];
    Kc = exp((g_RT[35] + g_RT[12]) - (g_RT[32] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[141] = phi_r * k_r;

    /*reaction 143: ch3choh + M <=> ch3hco + h + M */
    phi_f = sc[35];
    alpha = mixture;
    k_f = 1e-06 * alpha * 1e+14*exp(-12580.416645117993539/tc[1]);
    q_f[142] = phi_f * k_f;
    phi_r = sc[32]*sc[1];
    Kc = refC * exp((g_RT[35]) - (g_RT[32] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[142] = phi_r * k_r;

    /*reaction 144: ch3hco + oh <=> ch3co + h2o */
    phi_f = sc[32]*sc[12];
    k_f = 1e-06 * 9.24e+06*exp(1.5*tc[0]+484.0944325041403431/tc[1]);
    q_f[143] = phi_f * k_f;
    phi_r = sc[31]*sc[15];
    Kc = exp((g_RT[32] + g_RT[12]) - (g_RT[31] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[143] = phi_r * k_r;

    /*reaction 145: ch3hco + oh <=> ch2hco + h2o */
    phi_f = sc[32]*sc[12];
    k_f = 1e-06 * 172000*exp(2.4*tc[0]-410.12158263084654664/tc[1]);
    q_f[144] = phi_f * k_f;
    phi_r = sc[30]*sc[15];
    Kc = exp((g_RT[32] + g_RT[12]) - (g_RT[30] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[144] = phi_r * k_r;

    /*reaction 146: ch3hco + o <=> ch3co + oh */
    phi_f = sc[32]*sc[11];
    k_f = 1e-06 * 1.77e+18*exp(-1.9*tc[0]-1497.0695807690410675/tc[1]);
    q_f[145] = phi_f * k_f;
    phi_r = sc[31]*sc[12];
    Kc = exp((g_RT[32] + g_RT[11]) - (g_RT[31] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[145] = phi_r * k_r;

    /*reaction 147: ch3hco + o <=> ch2hco + oh */
    phi_f = sc[32]*sc[11];
    k_f = 1e-06 * 3.72e+13*exp(-0.2*tc[0]-1789.4384636015831802/tc[1]);
    q_f[146] = phi_f * k_f;
    phi_r = sc[30]*sc[12];
    Kc = exp((g_RT[32] + g_RT[11]) - (g_RT[30] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[146] = phi_r * k_r;

    /*reaction 148: ch3hco + h <=> ch3co + h2 */
    phi_f = sc[32]*sc[1];
    k_f = 1e-06 * 4.66e+13*exp(-0.35*tc[0]-1503.6113974245024565/tc[1]);
    q_f[147] = phi_f * k_f;
    phi_r = sc[31]*sc[0];
    Kc = exp((g_RT[32] + g_RT[1]) - (g_RT[31] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[147] = phi_r * k_r;

    /*reaction 149: ch3hco + h <=> ch2hco + h2 */
    phi_f = sc[32]*sc[1];
    k_f = 1e-06 * 1.85e+12*exp(0.4*tc[0]-2696.7381120474929048/tc[1]);
    q_f[148] = phi_f * k_f;
    phi_r = sc[30]*sc[0];
    Kc = exp((g_RT[32] + g_RT[1]) - (g_RT[30] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[148] = phi_r * k_r;

    /*reaction 150: ch3hco + ch3 <=> ch3co + ch4 */
    phi_f = sc[32]*sc[3];
    k_f = 1e-06 * 3.9e-07*exp(5.8*tc[0]-1107.0766647703833314/tc[1]);
    q_f[149] = phi_f * k_f;
    phi_r = sc[31]*sc[2];
    Kc = exp((g_RT[32] + g_RT[3]) - (g_RT[31] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[149] = phi_r * k_r;

    /*reaction 151: ch3hco + ch3 <=> ch2hco + ch4 */
    phi_f = sc[32]*sc[3];
    k_f = 1e-06 * 24.5*exp(3.15*tc[0]-2881.9218450636299167/tc[1]);
    q_f[150] = phi_f * k_f;
    phi_r = sc[30]*sc[2];
    Kc = exp((g_RT[32] + g_RT[3]) - (g_RT[30] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[150] = phi_r * k_r;

    /*reaction 152: ch3hco + ho2 <=> ch3co + h2o2 */
    phi_f = sc[32]*sc[13];
    k_f = 1e-06 * 2.4e+19*exp(-2.2*tc[0]-7060.1298212402170975/tc[1]);
    q_f[151] = phi_f * k_f;
    phi_r = sc[31]*sc[14];
    Kc = exp((g_RT[32] + g_RT[13]) - (g_RT[31] + g_RT[14]));
    k_r = k_f / Kc;
    q_r[151] = phi_r * k_r;

    /*reaction 153: ch3hco + ho2 <=> ch2hco + h2o2 */
    phi_f = sc[32]*sc[13];
    k_f = 1e-06 * 2.32e+11*exp(0.4*tc[0]-7479.8125205213536901/tc[1]);
    q_f[152] = phi_f * k_f;
    phi_r = sc[30]*sc[14];
    Kc = exp((g_RT[32] + g_RT[13]) - (g_RT[30] + g_RT[14]));
    k_r = k_f / Kc;
    q_r[152] = phi_r * k_r;

    /*reaction 154: ch3hco + o2 <=> ch3co + ho2 */
    phi_f = sc[32]*sc[10];
    k_f = 1e-06 * 1e+14*exp(-21235.743296959171857/tc[1]);
    q_f[153] = phi_f * k_f;
    phi_r = sc[31]*sc[13];
    Kc = exp((g_RT[32] + g_RT[10]) - (g_RT[31] + g_RT[13]));
    k_r = k_f / Kc;
    q_r[153] = phi_r * k_r;

    /*reaction 155: c2h6 + ch3 <=> c2h5 + ch4 */
    phi_f = sc[22]*sc[3];
    k_f = 1e-06 * 0.55*exp(4*tc[0]-4176.6983261791738187/tc[1]);
    q_f[154] = phi_f * k_f;
    phi_r = sc[21]*sc[2];
    Kc = exp((g_RT[22] + g_RT[3]) - (g_RT[21] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[154] = phi_r * k_r;

    /*reaction 156: c2h6 + h <=> c2h5 + h2 */
    phi_f = sc[22]*sc[1];
    k_f = 1e-06 * 540*exp(3.5*tc[0]-2621.7588288425899918/tc[1]);
    q_f[155] = phi_f * k_f;
    phi_r = sc[21]*sc[0];
    Kc = exp((g_RT[22] + g_RT[1]) - (g_RT[21] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[155] = phi_r * k_r;

    /*reaction 157: c2h6 + o <=> c2h5 + oh */
    phi_f = sc[22]*sc[11];
    k_f = 1e-06 * 3e+07*exp(2*tc[0]-2573.9532455911412399/tc[1]);
    q_f[156] = phi_f * k_f;
    phi_r = sc[21]*sc[12];
    Kc = exp((g_RT[22] + g_RT[11]) - (g_RT[21] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[156] = phi_r * k_r;

    /*reaction 158: c2h6 + oh <=> c2h5 + h2o */
    phi_f = sc[22]*sc[12];
    k_f = 1e-06 * 7.23e+06*exp(2*tc[0]-434.77919925527783107/tc[1]);
    q_f[157] = phi_f * k_f;
    phi_r = sc[21]*sc[15];
    Kc = exp((g_RT[22] + g_RT[12]) - (g_RT[21] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[157] = phi_r * k_r;

    /*reaction 159: c2h5 + h <=> c2h4 + h2 */
    phi_f = sc[21]*sc[1];
    k_f = 1e-06 * 1.25e+14*exp(-4025.733326437757114/tc[1]);
    q_f[158] = phi_f * k_f;
    phi_r = sc[20]*sc[0];
    Kc = exp((g_RT[21] + g_RT[1]) - (g_RT[20] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[158] = phi_r * k_r;

    /*reaction 160: c2h5 + h <=> ch3 + ch3 */
    phi_f = sc[21]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f[159] = phi_f * k_f;
    phi_r = sc[3]*sc[3];
    Kc = exp((g_RT[21] + g_RT[1]) - (g_RT[3] + g_RT[3]));
    k_r = k_f / Kc;
    q_r[159] = phi_r * k_r;

    /*reaction 161: c2h5 + h <=> c2h6 */
    phi_f = sc[21]*sc[1];
    k_f = 1e-06 * 3e+13;
    q_f[160] = phi_f * k_f;
    phi_r = sc[22];
    Kc = 1.0 / (refC) * exp((g_RT[21] + g_RT[1]) - (g_RT[22]));
    k_r = k_f / Kc;
    q_r[160] = phi_r * k_r;

    /*reaction 162: c2h5 + oh <=> c2h4 + h2o */
    phi_f = sc[21]*sc[12];
    k_f = 1e-06 * 4e+13;
    q_f[161] = phi_f * k_f;
    phi_r = sc[20]*sc[15];
    Kc = exp((g_RT[21] + g_RT[12]) - (g_RT[20] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[161] = phi_r * k_r;

    /*reaction 163: c2h5 + o <=> ch3 + ch2o */
    phi_f = sc[21]*sc[11];
    k_f = 1e-06 * 1e+14;
    q_f[162] = phi_f * k_f;
    phi_r = sc[3]*sc[6];
    Kc = exp((g_RT[21] + g_RT[11]) - (g_RT[3] + g_RT[6]));
    k_r = k_f / Kc;
    q_r[162] = phi_r * k_r;

    /*reaction 164: c2h5 + ho2 <=> c2h6 + o2 */
    phi_f = sc[21]*sc[13];
    k_f = 1e-06 * 3e+12;
    q_f[163] = phi_f * k_f;
    phi_r = sc[22]*sc[10];
    Kc = exp((g_RT[21] + g_RT[13]) - (g_RT[22] + g_RT[10]));
    k_r = k_f / Kc;
    q_r[163] = phi_r * k_r;

    /*reaction 165: c2h5 + ho2 <=> ch3ch2o + oh */
    phi_f = sc[21]*sc[13];
    k_f = 1e-06 * 3e+13;
    q_f[164] = phi_f * k_f;
    phi_r = sc[36]*sc[12];
    Kc = exp((g_RT[21] + g_RT[13]) - (g_RT[36] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[164] = phi_r * k_r;

    /*reaction 166: c2h5 + o2 <=> c2h4 + ho2 */
    phi_f = sc[21]*sc[10];
    k_f = 1e-06 * 2.89e+28*exp(-5.4*tc[0]-3816.8984101287987869/tc[1]);
    q_f[165] = phi_f * k_f;
    phi_r = sc[20]*sc[13];
    Kc = exp((g_RT[21] + g_RT[10]) - (g_RT[20] + g_RT[13]));
    k_r = k_f / Kc;
    q_r[165] = phi_r * k_r;

    /*reaction 167: c2h5 + o2 <=> ch3hco + oh */
    phi_f = sc[21]*sc[10];
    k_f = 1e-06 * 4.9e+11*exp(-0.48*tc[0]-4205.3816761300422513/tc[1]);
    q_f[166] = phi_f * k_f;
    phi_r = sc[32]*sc[12];
    Kc = exp((g_RT[21] + g_RT[10]) - (g_RT[32] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[166] = phi_r * k_r;

    /*reaction 168: c2h4 + oh <=> c2h4oh */
    phi_f = sc[20]*sc[12];
    k_f = 1e-06 * 1.29e+12*exp(+411.12801596245600422/tc[1]);
    q_f[167] = phi_f * k_f;
    phi_r = sc[34];
    Kc = 1.0 / (refC) * exp((g_RT[20] + g_RT[12]) - (g_RT[34]));
    k_r = k_f / Kc;
    q_r[167] = phi_r * k_r;

    /*reaction 169: c2h4 + oh <=> c2h3 + h2o */
    phi_f = sc[20]*sc[12];
    k_f = 1e-06 * 2.02e+13*exp(-2987.0941282168159887/tc[1]);
    q_f[168] = phi_f * k_f;
    phi_r = sc[19]*sc[15];
    Kc = exp((g_RT[20] + g_RT[12]) - (g_RT[19] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[168] = phi_r * k_r;

    /*reaction 170: c2h4 + o <=> ch3 + hco */
    phi_f = sc[20]*sc[11];
    k_f = 1e-06 * 1.02e+07*exp(1.88*tc[0]-90.07578317904483356/tc[1]);
    q_f[169] = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = exp((g_RT[20] + g_RT[11]) - (g_RT[3] + g_RT[7]));
    k_r = k_f / Kc;
    q_r[169] = phi_r * k_r;

    /*reaction 171: c2h4 + o <=> ch2hco + h */
    phi_f = sc[20]*sc[11];
    k_f = 1e-06 * 3.39e+06*exp(1.88*tc[0]-90.07578317904483356/tc[1]);
    q_f[170] = phi_f * k_f;
    phi_r = sc[30]*sc[1];
    Kc = exp((g_RT[20] + g_RT[11]) - (g_RT[30] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[170] = phi_r * k_r;

    /*reaction 172: c2h4 + ch3 <=> c2h3 + ch4 */
    phi_f = sc[20]*sc[3];
    k_f = 1e-06 * 6.62*exp(3.7*tc[0]-4780.5583251448369992/tc[1]);
    q_f[171] = phi_f * k_f;
    phi_r = sc[19]*sc[2];
    Kc = exp((g_RT[20] + g_RT[3]) - (g_RT[19] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[171] = phi_r * k_r;

    /*reaction 173: c2h4 + h <=> c2h3 + h2 */
    phi_f = sc[20]*sc[1];
    k_f = 1e-06 * 3.36e-07*exp(6*tc[0]-851.44259854158576672/tc[1]);
    q_f[172] = phi_f * k_f;
    phi_r = sc[19]*sc[0];
    Kc = exp((g_RT[20] + g_RT[1]) - (g_RT[19] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[172] = phi_r * k_r;

    /*reaction 174: c2h4 + h (+M) <=> c2h5 (+M) */
    phi_f = sc[20]*sc[1];
    alpha = mixture + 4*sc[15] + sc[0] + 2*sc[8] + sc[9];
    k_f = 1e-06 * 1.08e+12*exp(0.454*tc[0]-916.86076509619920216/tc[1]);
    redP = 1e-12 * alpha / k_f * 1.112e+34*exp(-5*tc[0]-2238.3077294993931901/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0*exp(T/-1e-15))+ (1*exp(T/-95))+ (exp(-200/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[173] = phi_f * k_f;
    phi_r = sc[21];
    Kc = 1.0 / (refC) * exp((g_RT[20] + g_RT[1]) - (g_RT[21]));
    k_r = k_f / Kc;
    q_r[173] = phi_r * k_r;

    /*reaction 175: c2h4 (+M) <=> c2h2 + h2 (+M) */
    phi_f = sc[20];
    alpha = mixture;
    k_f = 1 * 1.8e+14*exp(-43779.849925010610605/tc[1]);
    redP = 1e-6 * alpha / k_f * 1.5e+15*exp(-27899.841602211075951/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f[174] = phi_f * k_f;
    phi_r = sc[18]*sc[0];
    Kc = refC * exp((g_RT[20]) - (g_RT[18] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[174] = phi_r * k_r;

    /*reaction 176: c2h3 + h (+M) <=> c2h4 (+M) */
    phi_f = sc[19]*sc[1];
    alpha = mixture + 4*sc[15];
    k_f = 1e-06 * 6.1e+12*exp(0.27*tc[0]-140.90066642532150354/tc[1]);
    redP = 1e-12 * alpha / k_f * 9.8e+29*exp(-3.86*tc[0]-1670.6793304716693456/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.218*exp(T/-208))+ (0.782*exp(T/-2663))+ (exp(-6095/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[175] = phi_f * k_f;
    phi_r = sc[20];
    Kc = 1.0 / (refC) * exp((g_RT[19] + g_RT[1]) - (g_RT[20]));
    k_r = k_f / Kc;
    q_r[175] = phi_r * k_r;

    /*reaction 177: c2h3 + h <=> c2h2 + h2 */
    phi_f = sc[19]*sc[1];
    k_f = 1e-06 * 9e+13;
    q_f[176] = phi_f * k_f;
    phi_r = sc[18]*sc[0];
    Kc = exp((g_RT[19] + g_RT[1]) - (g_RT[18] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[176] = phi_r * k_r;

    /*reaction 178: c2h3 + o <=> ch2co + h */
    phi_f = sc[19]*sc[11];
    k_f = 1e-06 * 3e+13;
    q_f[177] = phi_f * k_f;
    phi_r = sc[29]*sc[1];
    Kc = exp((g_RT[19] + g_RT[11]) - (g_RT[29] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[177] = phi_r * k_r;

    /*reaction 179: c2h3 + o2 <=> ch2o + hco */
    phi_f = sc[19]*sc[10];
    k_f = 1e-06 * 1.7e+29*exp(-5.312*tc[0]-3270.9083277306781383/tc[1]);
    q_f[178] = phi_f * k_f;
    phi_r = sc[6]*sc[7];
    Kc = exp((g_RT[19] + g_RT[10]) - (g_RT[6] + g_RT[7]));
    k_r = k_f / Kc;
    q_r[178] = phi_r * k_r;

    /*reaction 180: c2h3 + o2 <=> ch2hco + o */
    phi_f = sc[19]*sc[10];
    k_f = 1e-06 * 5.5e+14*exp(-0.611*tc[0]-2646.9196621328255787/tc[1]);
    q_f[179] = phi_f * k_f;
    phi_r = sc[30]*sc[11];
    Kc = exp((g_RT[19] + g_RT[10]) - (g_RT[30] + g_RT[11]));
    k_r = k_f / Kc;
    q_r[179] = phi_r * k_r;

    /*reaction 181: c2h3 + o2 <=> c2h2 + ho2 */
    phi_f = sc[19]*sc[10];
    k_f = 1e-06 * 2.12e-06*exp(6*tc[0]-4772.5068584919617933/tc[1]);
    q_f[180] = phi_f * k_f;
    phi_r = sc[18]*sc[13];
    Kc = exp((g_RT[19] + g_RT[10]) - (g_RT[18] + g_RT[13]));
    k_r = k_f / Kc;
    q_r[180] = phi_r * k_r;

    /*reaction 182: c2h3 + oh <=> c2h2 + h2o */
    phi_f = sc[19]*sc[12];
    k_f = 1e-06 * 2e+13;
    q_f[181] = phi_f * k_f;
    phi_r = sc[18]*sc[15];
    Kc = exp((g_RT[19] + g_RT[12]) - (g_RT[18] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[181] = phi_r * k_r;

    /*reaction 183: c2h3 + c2h <=> c2h2 + c2h2 */
    phi_f = sc[19]*sc[16];
    k_f = 1e-06 * 3e+13;
    q_f[182] = phi_f * k_f;
    phi_r = sc[18]*sc[18];
    Kc = exp((g_RT[19] + g_RT[16]) - (g_RT[18] + g_RT[18]));
    k_r = k_f / Kc;
    q_r[182] = phi_r * k_r;

    /*reaction 184: c2h3 + ch <=> ch2 + c2h2 */
    phi_f = sc[19]*sc[5];
    k_f = 1e-06 * 5e+13;
    q_f[183] = phi_f * k_f;
    phi_r = sc[4]*sc[18];
    Kc = exp((g_RT[19] + g_RT[5]) - (g_RT[4] + g_RT[18]));
    k_r = k_f / Kc;
    q_r[183] = phi_r * k_r;

    /*reaction 185: c2h3 + ch3 <=> c2h2 + ch4 */
    phi_f = sc[19]*sc[3];
    k_f = 1e-06 * 2e+13;
    q_f[184] = phi_f * k_f;
    phi_r = sc[18]*sc[2];
    Kc = exp((g_RT[19] + g_RT[3]) - (g_RT[18] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[184] = phi_r * k_r;

    /*reaction 186: c2h2 + oh <=> c2h + h2o */
    phi_f = sc[18]*sc[12];
    k_f = 1e-06 * 3.37e+07*exp(2*tc[0]-7045.0333212660752906/tc[1]);
    q_f[185] = phi_f * k_f;
    phi_r = sc[16]*sc[15];
    Kc = exp((g_RT[18] + g_RT[12]) - (g_RT[16] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[185] = phi_r * k_r;

    /*reaction 187: c2h2 + oh <=> hccoh + h */
    phi_f = sc[18]*sc[12];
    k_f = 1e-06 * 504000*exp(2.3*tc[0]-6793.4249883637157836/tc[1]);
    q_f[186] = phi_f * k_f;
    phi_r = sc[25]*sc[1];
    Kc = exp((g_RT[18] + g_RT[12]) - (g_RT[25] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[186] = phi_r * k_r;

    /*reaction 188: c2h2 + oh <=> ch2co + h */
    phi_f = sc[18]*sc[12];
    k_f = 1e-06 * 0.000218*exp(4.5*tc[0]+503.21666580471963925/tc[1]);
    q_f[187] = phi_f * k_f;
    phi_r = sc[29]*sc[1];
    Kc = exp((g_RT[18] + g_RT[12]) - (g_RT[29] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[187] = phi_r * k_r;

    /*reaction 189: c2h2 + oh <=> ch2co + h */
    phi_f = sc[18]*sc[12];
    k_f = 1e-06 * 2e+11;
    q_f[188] = phi_f * k_f;
    phi_r = sc[29]*sc[1];
    Kc = exp((g_RT[18] + g_RT[12]) - (g_RT[29] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[188] = phi_r * k_r;

    /*reaction 190: c2h2 + oh <=> ch3 + co */
    phi_f = sc[18]*sc[12];
    k_f = 1e-06 * 0.000483*exp(4*tc[0]+1006.4333316094392785/tc[1]);
    q_f[189] = phi_f * k_f;
    phi_r = sc[3]*sc[9];
    Kc = exp((g_RT[18] + g_RT[12]) - (g_RT[3] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[189] = phi_r * k_r;

    /*reaction 191: hccoh + h <=> ch2co + h */
    phi_f = sc[25]*sc[1];
    k_f = 1e-06 * 1e+13;
    q_f[190] = phi_f * k_f;
    phi_r = sc[29]*sc[1];
    Kc = exp((g_RT[25] + g_RT[1]) - (g_RT[29] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[190] = phi_r * k_r;

    /*reaction 192: c2h2 + o <=> ch2 + co */
    phi_f = sc[18]*sc[11];
    k_f = 1e-06 * 6.12e+06*exp(2*tc[0]-956.11166502896742259/tc[1]);
    q_f[191] = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = exp((g_RT[18] + g_RT[11]) - (g_RT[4] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[191] = phi_r * k_r;

    /*reaction 193: c2h2 + o <=> hcco + h */
    phi_f = sc[18]*sc[11];
    k_f = 1e-06 * 1.43e+07*exp(2*tc[0]-956.11166502896742259/tc[1]);
    q_f[192] = phi_f * k_f;
    phi_r = sc[17]*sc[1];
    Kc = exp((g_RT[18] + g_RT[11]) - (g_RT[17] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[192] = phi_r * k_r;

    /*reaction 194: c2h2 + o <=> c2h + oh */
    phi_f = sc[18]*sc[11];
    k_f = 1e-06 * 3.16e+15*exp(-0.6*tc[0]-7548.2499870707952141/tc[1]);
    q_f[193] = phi_f * k_f;
    phi_r = sc[16]*sc[12];
    Kc = exp((g_RT[18] + g_RT[11]) - (g_RT[16] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[193] = phi_r * k_r;

    /*reaction 195: c2h2 + ch3 <=> c2h + ch4 */
    phi_f = sc[18]*sc[3];
    k_f = 1e-06 * 1.81e+11*exp(-8700.1129350977989816/tc[1]);
    q_f[194] = phi_f * k_f;
    phi_r = sc[16]*sc[2];
    Kc = exp((g_RT[18] + g_RT[3]) - (g_RT[16] + g_RT[2]));
    k_r = k_f / Kc;
    q_r[194] = phi_r * k_r;

    /*reaction 196: c2h2 + o2 <=> hcco + oh */
    phi_f = sc[18]*sc[10];
    k_f = 1e-06 * 4e+07*exp(1.5*tc[0]-15146.821640722062511/tc[1]);
    q_f[195] = phi_f * k_f;
    phi_r = sc[17]*sc[12];
    Kc = exp((g_RT[18] + g_RT[10]) - (g_RT[17] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[195] = phi_r * k_r;

    /*reaction 197: c2h2 + M <=> c2h + h + M */
    phi_f = sc[18];
    alpha = mixture;
    k_f = 1e-06 * alpha * 4.2e+16*exp(-53844.183241105005436/tc[1]);
    q_f[196] = phi_f * k_f;
    phi_r = sc[16]*sc[1];
    Kc = refC * exp((g_RT[18]) - (g_RT[16] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[196] = phi_r * k_r;

    /*reaction 198: c2h2 + h (+M) <=> c2h3 (+M) */
    phi_f = sc[18]*sc[1];
    alpha = mixture + 4*sc[15] + sc[0] + 2*sc[8] + sc[9];
    k_f = 1e-06 * 3.11e+11*exp(0.58*tc[0]-1302.8279477684191079/tc[1]);
    redP = 1e-12 * alpha / k_f * 2.25e+40*exp(-7.269*tc[0]-3309.6560109976417152/tc[1]);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0*exp(T/-1e-15))+ (1*exp(T/-675))+ (exp(-1e+15/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[197] = phi_f * k_f;
    phi_r = sc[19];
    Kc = 1.0 / (refC) * exp((g_RT[18] + g_RT[1]) - (g_RT[19]));
    k_r = k_f / Kc;
    q_r[197] = phi_r * k_r;

    /*reaction 199: ch2hco + h <=> ch3 + hco */
    phi_f = sc[30]*sc[1];
    k_f = 1e-06 * 5e+13;
    q_f[198] = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = exp((g_RT[30] + g_RT[1]) - (g_RT[3] + g_RT[7]));
    k_r = k_f / Kc;
    q_r[198] = phi_r * k_r;

    /*reaction 200: ch2hco + h <=> ch2co + h2 */
    phi_f = sc[30]*sc[1];
    k_f = 1e-06 * 2e+13;
    q_f[199] = phi_f * k_f;
    phi_r = sc[29]*sc[0];
    Kc = exp((g_RT[30] + g_RT[1]) - (g_RT[29] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[199] = phi_r * k_r;

    /*reaction 201: ch2hco + o <=> ch2o + hco */
    phi_f = sc[30]*sc[11];
    k_f = 1e-06 * 1e+14;
    q_f[200] = phi_f * k_f;
    phi_r = sc[6]*sc[7];
    Kc = exp((g_RT[30] + g_RT[11]) - (g_RT[6] + g_RT[7]));
    k_r = k_f / Kc;
    q_r[200] = phi_r * k_r;

    /*reaction 202: ch2hco + oh <=> ch2co + h2o */
    phi_f = sc[30]*sc[12];
    k_f = 1e-06 * 3e+13;
    q_f[201] = phi_f * k_f;
    phi_r = sc[29]*sc[15];
    Kc = exp((g_RT[30] + g_RT[12]) - (g_RT[29] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[201] = phi_r * k_r;

    /*reaction 203: ch2hco + o2 <=> ch2o + co + oh */
    phi_f = sc[30]*sc[10];
    k_f = 1e-06 * 3e+10;
    q_f[202] = phi_f * k_f;
    phi_r = sc[6]*sc[9]*sc[12];
    Kc = refC * exp((g_RT[30] + g_RT[10]) - (g_RT[6] + g_RT[9] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[202] = phi_r * k_r;

    /*reaction 204: ch2hco + ch3 <=> c2h5 + co + h */
    phi_f = sc[30]*sc[3];
    k_f = 1e-06 * 4.9e+14*exp(-0.5*tc[0]);
    q_f[203] = phi_f * k_f;
    phi_r = sc[21]*sc[9]*sc[1];
    Kc = refC * exp((g_RT[30] + g_RT[3]) - (g_RT[21] + g_RT[9] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[203] = phi_r * k_r;

    /*reaction 205: ch2hco + ho2 <=> ch2o + hco + oh */
    phi_f = sc[30]*sc[13];
    k_f = 1e-06 * 7e+12;
    q_f[204] = phi_f * k_f;
    phi_r = sc[6]*sc[7]*sc[12];
    Kc = refC * exp((g_RT[30] + g_RT[13]) - (g_RT[6] + g_RT[7] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[204] = phi_r * k_r;

    /*reaction 206: ch2hco + ho2 <=> ch3hco + o2 */
    phi_f = sc[30]*sc[13];
    k_f = 1e-06 * 3e+12;
    q_f[205] = phi_f * k_f;
    phi_r = sc[32]*sc[10];
    Kc = exp((g_RT[30] + g_RT[13]) - (g_RT[32] + g_RT[10]));
    k_r = k_f / Kc;
    q_r[205] = phi_r * k_r;

    /*reaction 207: ch2hco <=> ch3 + co */
    phi_f = sc[30];
    k_f = 1 * 1.17e+43*exp(-9.83*tc[0]-22018.748428951312235/tc[1]);
    q_f[206] = phi_f * k_f;
    phi_r = sc[3]*sc[9];
    Kc = refC * exp((g_RT[30]) - (g_RT[3] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[206] = phi_r * k_r;

    /*reaction 208: ch2hco <=> ch2co + h */
    phi_f = sc[30];
    k_f = 1 * 1.81e+43*exp(-9.61*tc[0]-23081.542027130883071/tc[1]);
    q_f[207] = phi_f * k_f;
    phi_r = sc[29]*sc[1];
    Kc = refC * exp((g_RT[30]) - (g_RT[29] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[207] = phi_r * k_r;

    /*reaction 209: ch3co (+M) <=> ch3 + co (+M) */
    phi_f = sc[31];
    alpha = mixture;
    k_f = 1 * 3e+12*exp(-8414.7890855865243793/tc[1]);
    redP = 1e-6 * alpha / k_f * 1.2e+15*exp(-6299.2662225434805805/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f[208] = phi_f * k_f;
    phi_r = sc[3]*sc[9];
    Kc = refC * exp((g_RT[31]) - (g_RT[3] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[208] = phi_r * k_r;

    /*reaction 210: ch2co + o <=> co2 + ch2 */
    phi_f = sc[29]*sc[11];
    k_f = 1e-06 * 1.75e+12*exp(-679.34249883637164658/tc[1]);
    q_f[209] = phi_f * k_f;
    phi_r = sc[8]*sc[4];
    Kc = exp((g_RT[29] + g_RT[11]) - (g_RT[8] + g_RT[4]));
    k_r = k_f / Kc;
    q_r[209] = phi_r * k_r;

    /*reaction 211: ch2co + h <=> ch3 + co */
    phi_f = sc[29]*sc[1];
    k_f = 1e-06 * 27100*exp(2.75*tc[0]-359.29669938456987666/tc[1]);
    q_f[210] = phi_f * k_f;
    phi_r = sc[3]*sc[9];
    Kc = exp((g_RT[29] + g_RT[1]) - (g_RT[3] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[210] = phi_r * k_r;

    /*reaction 212: ch2co + h <=> hcco + h2 */
    phi_f = sc[29]*sc[1];
    k_f = 1e-06 * 2e+14*exp(-4025.733326437757114/tc[1]);
    q_f[211] = phi_f * k_f;
    phi_r = sc[17]*sc[0];
    Kc = exp((g_RT[29] + g_RT[1]) - (g_RT[17] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[211] = phi_r * k_r;

    /*reaction 213: ch2co + o <=> hcco + oh */
    phi_f = sc[29]*sc[11];
    k_f = 1e-06 * 1e+13*exp(-4025.733326437757114/tc[1]);
    q_f[212] = phi_f * k_f;
    phi_r = sc[17]*sc[12];
    Kc = exp((g_RT[29] + g_RT[11]) - (g_RT[17] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[212] = phi_r * k_r;

    /*reaction 214: ch2co + oh <=> hcco + h2o */
    phi_f = sc[29]*sc[12];
    k_f = 1e-06 * 1e+13*exp(-1006.4333316094392785/tc[1]);
    q_f[213] = phi_f * k_f;
    phi_r = sc[17]*sc[15];
    Kc = exp((g_RT[29] + g_RT[12]) - (g_RT[17] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[213] = phi_r * k_r;

    /*reaction 215: ch2co + oh <=> ch2oh + co */
    phi_f = sc[29]*sc[12];
    k_f = 1e-06 * 3.73e+12*exp(+509.75848246018102827/tc[1]);
    q_f[214] = phi_f * k_f;
    phi_r = sc[23]*sc[9];
    Kc = exp((g_RT[29] + g_RT[12]) - (g_RT[23] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[214] = phi_r * k_r;

    /*reaction 216: ch2co (+M) <=> ch2 + co (+M) */
    phi_f = sc[29];
    alpha = mixture;
    k_f = 1 * 3e+14*exp(-35718.318938819007599/tc[1]);
    redP = 1e-6 * alpha / k_f * 3.6e+15*exp(-29825.651782245739014/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f[215] = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = refC * exp((g_RT[29]) - (g_RT[4] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[215] = phi_r * k_r;

    /*reaction 217: c2h + h2 <=> c2h2 + h */
    phi_f = sc[16]*sc[0];
    k_f = 1e-06 * 409000*exp(2.39*tc[0]-434.93016425501917865/tc[1]);
    q_f[216] = phi_f * k_f;
    phi_r = sc[18]*sc[1];
    Kc = exp((g_RT[16] + g_RT[0]) - (g_RT[18] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[216] = phi_r * k_r;

    /*reaction 218: c2h + o <=> ch + co */
    phi_f = sc[16]*sc[11];
    k_f = 1e-06 * 5e+13;
    q_f[217] = phi_f * k_f;
    phi_r = sc[5]*sc[9];
    Kc = exp((g_RT[16] + g_RT[11]) - (g_RT[5] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[217] = phi_r * k_r;

    /*reaction 219: c2h + oh <=> hcco + h */
    phi_f = sc[16]*sc[12];
    k_f = 1e-06 * 2e+13;
    q_f[218] = phi_f * k_f;
    phi_r = sc[17]*sc[1];
    Kc = exp((g_RT[16] + g_RT[12]) - (g_RT[17] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[218] = phi_r * k_r;

    /*reaction 220: c2h + o2 <=> co + co + h */
    phi_f = sc[16]*sc[10];
    k_f = 1e-06 * 9.04e+12*exp(+229.97001627275690794/tc[1]);
    q_f[219] = phi_f * k_f;
    phi_r = sc[9]*sc[9]*sc[1];
    Kc = refC * exp((g_RT[16] + g_RT[10]) - (g_RT[9] + g_RT[9] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[219] = phi_r * k_r;

    /*reaction 221: hcco + c2h2 <=> h2ccch + co */
    phi_f = sc[17]*sc[18];
    k_f = 1e-06 * 1e+11*exp(-1509.6499974141590883/tc[1]);
    q_f[220] = phi_f * k_f;
    phi_r = sc[26]*sc[9];
    Kc = exp((g_RT[17] + g_RT[18]) - (g_RT[26] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[220] = phi_r * k_r;

    /*reaction 222: hcco + h <=> ch2s + co */
    phi_f = sc[17]*sc[1];
    k_f = 1e-06 * 1e+14;
    q_f[221] = phi_f * k_f;
    phi_r = sc[28]*sc[9];
    Kc = exp((g_RT[17] + g_RT[1]) - (g_RT[28] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[221] = phi_r * k_r;

    /*reaction 223: hcco + o <=> h + co + co */
    phi_f = sc[17]*sc[11];
    k_f = 1e-06 * 8e+13;
    q_f[222] = phi_f * k_f;
    phi_r = sc[1]*sc[9]*sc[9];
    Kc = refC * exp((g_RT[17] + g_RT[11]) - (g_RT[1] + g_RT[9] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[222] = phi_r * k_r;

    /*reaction 224: hcco + o <=> ch + co2 */
    phi_f = sc[17]*sc[11];
    k_f = 1e-06 * 2.95e+13*exp(-560.08014904065305473/tc[1]);
    q_f[223] = phi_f * k_f;
    phi_r = sc[5]*sc[8];
    Kc = exp((g_RT[17] + g_RT[11]) - (g_RT[5] + g_RT[8]));
    k_r = k_f / Kc;
    q_r[223] = phi_r * k_r;

    /*reaction 225: hcco + o2 <=> hco + co + o */
    phi_f = sc[17]*sc[10];
    k_f = 1e-06 * 2.5e+08*exp(1*tc[0]);
    q_f[224] = phi_f * k_f;
    phi_r = sc[7]*sc[9]*sc[11];
    Kc = refC * exp((g_RT[17] + g_RT[10]) - (g_RT[7] + g_RT[9] + g_RT[11]));
    k_r = k_f / Kc;
    q_r[224] = phi_r * k_r;

    /*reaction 226: hcco + o2 <=> co2 + hco */
    phi_f = sc[17]*sc[10];
    k_f = 1e-06 * 2.4e+11*exp(+429.7470325972306/tc[1]);
    q_f[225] = phi_f * k_f;
    phi_r = sc[8]*sc[7];
    Kc = exp((g_RT[17] + g_RT[10]) - (g_RT[8] + g_RT[7]));
    k_r = k_f / Kc;
    q_r[225] = phi_r * k_r;

    /*reaction 227: hcco + ch <=> c2h2 + co */
    phi_f = sc[17]*sc[5];
    k_f = 1e-06 * 5e+13;
    q_f[226] = phi_f * k_f;
    phi_r = sc[18]*sc[9];
    Kc = exp((g_RT[17] + g_RT[5]) - (g_RT[18] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[226] = phi_r * k_r;

    /*reaction 228: hcco + hcco <=> c2h2 + co + co */
    phi_f = sc[17]*sc[17];
    k_f = 1e-06 * 1e+13;
    q_f[227] = phi_f * k_f;
    phi_r = sc[18]*sc[9]*sc[9];
    Kc = refC * exp((g_RT[17] + g_RT[17]) - (g_RT[18] + g_RT[9] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[227] = phi_r * k_r;

    return;
}


/*compute the equilibrium constants for each reaction */
void equilibriumConstants(double *kc, double * g_RT, double T)
{
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    /*reaction 1: oh + h2 <=> h + h2o */
    kc[0] = exp((g_RT[12] + g_RT[0]) - (g_RT[1] + g_RT[15]));

    /*reaction 2: o + oh <=> o2 + h */
    kc[1] = exp((g_RT[11] + g_RT[12]) - (g_RT[10] + g_RT[1]));

    /*reaction 3: o + h2 <=> oh + h */
    kc[2] = exp((g_RT[11] + g_RT[0]) - (g_RT[12] + g_RT[1]));

    /*reaction 4: h + o2 (+M) <=> ho2 (+M) */
    kc[3] = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[13]));

    /*reaction 5: h + o2 (+n2) <=> ho2 (+n2) */
    kc[4] = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[13]));

    /*reaction 6: h + o2 (+h2) <=> ho2 (+h2) */
    kc[5] = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[13]));

    /*reaction 7: h + o2 (+h2o) <=> ho2 (+h2o) */
    kc[6] = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[13]));

    /*reaction 8: oh + ho2 <=> h2o + o2 */
    kc[7] = exp((g_RT[12] + g_RT[13]) - (g_RT[15] + g_RT[10]));

    /*reaction 9: oh + ho2 <=> h2o + o2 */
    kc[8] = exp((g_RT[12] + g_RT[13]) - (g_RT[15] + g_RT[10]));

    /*reaction 10: h + ho2 <=> oh + oh */
    kc[9] = exp((g_RT[1] + g_RT[13]) - (g_RT[12] + g_RT[12]));

    /*reaction 11: h + ho2 <=> h2 + o2 */
    kc[10] = exp((g_RT[1] + g_RT[13]) - (g_RT[0] + g_RT[10]));

    /*reaction 12: h + ho2 <=> o + h2o */
    kc[11] = exp((g_RT[1] + g_RT[13]) - (g_RT[11] + g_RT[15]));

    /*reaction 13: o + ho2 <=> o2 + oh */
    kc[12] = exp((g_RT[11] + g_RT[13]) - (g_RT[10] + g_RT[12]));

    /*reaction 14: oh + oh <=> o + h2o */
    kc[13] = exp((g_RT[12] + g_RT[12]) - (g_RT[11] + g_RT[15]));

    /*reaction 15: h + h + M <=> h2 + M */
    kc[14] = 1.0 / (refC) * exp((g_RT[1] + g_RT[1]) - (g_RT[0]));

    /*reaction 16: h + h + h2 <=> h2 + h2 */
    kc[15] = 1.0 / (refC) * exp((g_RT[1] + g_RT[1] + g_RT[0]) - (g_RT[0] + g_RT[0]));

    /*reaction 17: h + h + h2o <=> h2 + h2o */
    kc[16] = 1.0 / (refC) * exp((g_RT[1] + g_RT[1] + g_RT[15]) - (g_RT[0] + g_RT[15]));

    /*reaction 18: h + oh + M <=> h2o + M */
    kc[17] = 1.0 / (refC) * exp((g_RT[1] + g_RT[12]) - (g_RT[15]));

    /*reaction 19: h + o + M <=> oh + M */
    kc[18] = 1.0 / (refC) * exp((g_RT[1] + g_RT[11]) - (g_RT[12]));

    /*reaction 20: o + o + M <=> o2 + M */
    kc[19] = 1.0 / (refC) * exp((g_RT[11] + g_RT[11]) - (g_RT[10]));

    /*reaction 21: ho2 + ho2 <=> h2o2 + o2 */
    kc[20] = exp((g_RT[13] + g_RT[13]) - (g_RT[14] + g_RT[10]));

    /*reaction 22: ho2 + ho2 <=> h2o2 + o2 */
    kc[21] = exp((g_RT[13] + g_RT[13]) - (g_RT[14] + g_RT[10]));

    /*reaction 23: oh + oh (+M) <=> h2o2 (+M) */
    kc[22] = 1.0 / (refC) * exp((g_RT[12] + g_RT[12]) - (g_RT[14]));

    /*reaction 24: h2o2 + h <=> ho2 + h2 */
    kc[23] = exp((g_RT[14] + g_RT[1]) - (g_RT[13] + g_RT[0]));

    /*reaction 25: h2o2 + h <=> oh + h2o */
    kc[24] = exp((g_RT[14] + g_RT[1]) - (g_RT[12] + g_RT[15]));

    /*reaction 26: h2o2 + o <=> oh + ho2 */
    kc[25] = exp((g_RT[14] + g_RT[11]) - (g_RT[12] + g_RT[13]));

    /*reaction 27: h2o2 + oh <=> h2o + ho2 */
    kc[26] = exp((g_RT[14] + g_RT[12]) - (g_RT[15] + g_RT[13]));

    /*reaction 28: ch3 + ch3 (+M) <=> c2h6 (+M) */
    kc[27] = 1.0 / (refC) * exp((g_RT[3] + g_RT[3]) - (g_RT[22]));

    /*reaction 29: ch3 + h (+M) <=> ch4 (+M) */
    kc[28] = 1.0 / (refC) * exp((g_RT[3] + g_RT[1]) - (g_RT[2]));

    /*reaction 30: ch4 + h <=> ch3 + h2 */
    kc[29] = exp((g_RT[2] + g_RT[1]) - (g_RT[3] + g_RT[0]));

    /*reaction 31: ch4 + oh <=> ch3 + h2o */
    kc[30] = exp((g_RT[2] + g_RT[12]) - (g_RT[3] + g_RT[15]));

    /*reaction 32: ch4 + o <=> ch3 + oh */
    kc[31] = exp((g_RT[2] + g_RT[11]) - (g_RT[3] + g_RT[12]));

    /*reaction 33: ch4 + ho2 <=> ch3 + h2o2 */
    kc[32] = exp((g_RT[2] + g_RT[13]) - (g_RT[3] + g_RT[14]));

    /*reaction 34: ch3 + ho2 <=> ch3o + oh */
    kc[33] = exp((g_RT[3] + g_RT[13]) - (g_RT[24] + g_RT[12]));

    /*reaction 35: ch3 + ho2 <=> ch4 + o2 */
    kc[34] = exp((g_RT[3] + g_RT[13]) - (g_RT[2] + g_RT[10]));

    /*reaction 36: ch3 + o <=> ch2o + h */
    kc[35] = exp((g_RT[3] + g_RT[11]) - (g_RT[6] + g_RT[1]));

    /*reaction 37: ch3 + o2 <=> ch3o + o */
    kc[36] = exp((g_RT[3] + g_RT[10]) - (g_RT[24] + g_RT[11]));

    /*reaction 38: ch3 + o2 <=> ch2o + oh */
    kc[37] = exp((g_RT[3] + g_RT[10]) - (g_RT[6] + g_RT[12]));

    /*reaction 39: ch3o + h <=> ch3 + oh */
    kc[38] = exp((g_RT[24] + g_RT[1]) - (g_RT[3] + g_RT[12]));

    /*reaction 40: ch2oh + h <=> ch3 + oh */
    kc[39] = exp((g_RT[23] + g_RT[1]) - (g_RT[3] + g_RT[12]));

    /*reaction 41: ch3 + oh <=> ch2s + h2o */
    kc[40] = exp((g_RT[3] + g_RT[12]) - (g_RT[28] + g_RT[15]));

    /*reaction 42: ch3 + oh <=> ch2 + h2o */
    kc[41] = exp((g_RT[3] + g_RT[12]) - (g_RT[4] + g_RT[15]));

    /*reaction 43: ch3 + h <=> ch2 + h2 */
    kc[42] = exp((g_RT[3] + g_RT[1]) - (g_RT[4] + g_RT[0]));

    /*reaction 44: ch3 + M <=> ch + h2 + M */
    kc[43] = refC * exp((g_RT[3]) - (g_RT[5] + g_RT[0]));

    /*reaction 45: ch3 + M <=> ch2 + h + M */
    kc[44] = refC * exp((g_RT[3]) - (g_RT[4] + g_RT[1]));

    /*reaction 46: ch2o + h (+M) <=> ch3o (+M) */
    kc[45] = 1.0 / (refC) * exp((g_RT[6] + g_RT[1]) - (g_RT[24]));

    /*reaction 47: ch2o + h (+M) <=> ch2oh (+M) */
    kc[46] = 1.0 / (refC) * exp((g_RT[6] + g_RT[1]) - (g_RT[23]));

    /*reaction 48: ch3o + ch3 <=> ch2o + ch4 */
    kc[47] = exp((g_RT[24] + g_RT[3]) - (g_RT[6] + g_RT[2]));

    /*reaction 49: ch3o + h <=> ch2o + h2 */
    kc[48] = exp((g_RT[24] + g_RT[1]) - (g_RT[6] + g_RT[0]));

    /*reaction 50: ch2oh + h <=> ch2o + h2 */
    kc[49] = exp((g_RT[23] + g_RT[1]) - (g_RT[6] + g_RT[0]));

    /*reaction 51: ch3o + oh <=> ch2o + h2o */
    kc[50] = exp((g_RT[24] + g_RT[12]) - (g_RT[6] + g_RT[15]));

    /*reaction 52: ch2oh + oh <=> ch2o + h2o */
    kc[51] = exp((g_RT[23] + g_RT[12]) - (g_RT[6] + g_RT[15]));

    /*reaction 53: ch3o + o <=> ch2o + oh */
    kc[52] = exp((g_RT[24] + g_RT[11]) - (g_RT[6] + g_RT[12]));

    /*reaction 54: ch2oh + o <=> ch2o + oh */
    kc[53] = exp((g_RT[23] + g_RT[11]) - (g_RT[6] + g_RT[12]));

    /*reaction 55: ch3o + o2 <=> ch2o + ho2 */
    kc[54] = exp((g_RT[24] + g_RT[10]) - (g_RT[6] + g_RT[13]));

    /*reaction 56: ch3o + co <=> ch3 + co2 */
    kc[55] = exp((g_RT[24] + g_RT[9]) - (g_RT[3] + g_RT[8]));

    /*reaction 57: ch2oh + o2 <=> ch2o + ho2 */
    kc[56] = exp((g_RT[23] + g_RT[10]) - (g_RT[6] + g_RT[13]));

    /*reaction 58: ch2oh + o2 <=> ch2o + ho2 */
    kc[57] = exp((g_RT[23] + g_RT[10]) - (g_RT[6] + g_RT[13]));

    /*reaction 59: ch2 + h <=> ch + h2 */
    kc[58] = exp((g_RT[4] + g_RT[1]) - (g_RT[5] + g_RT[0]));

    /*reaction 60: ch2 + oh <=> ch + h2o */
    kc[59] = exp((g_RT[4] + g_RT[12]) - (g_RT[5] + g_RT[15]));

    /*reaction 61: ch2 + oh <=> ch2o + h */
    kc[60] = exp((g_RT[4] + g_RT[12]) - (g_RT[6] + g_RT[1]));

    /*reaction 62: ch2 + co2 <=> ch2o + co */
    kc[61] = exp((g_RT[4] + g_RT[8]) - (g_RT[6] + g_RT[9]));

    /*reaction 63: ch2 + o <=> co + h + h */
    kc[62] = refC * exp((g_RT[4] + g_RT[11]) - (g_RT[9] + g_RT[1] + g_RT[1]));

    /*reaction 64: ch2 + o <=> co + h2 */
    kc[63] = exp((g_RT[4] + g_RT[11]) - (g_RT[9] + g_RT[0]));

    /*reaction 65: ch2 + o2 <=> ch2o + o */
    kc[64] = exp((g_RT[4] + g_RT[10]) - (g_RT[6] + g_RT[11]));

    /*reaction 66: ch2 + o2 <=> co2 + h + h */
    kc[65] = refC * exp((g_RT[4] + g_RT[10]) - (g_RT[8] + g_RT[1] + g_RT[1]));

    /*reaction 67: ch2 + o2 <=> co2 + h2 */
    kc[66] = exp((g_RT[4] + g_RT[10]) - (g_RT[8] + g_RT[0]));

    /*reaction 68: ch2 + o2 <=> co + h2o */
    kc[67] = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[15]));

    /*reaction 69: ch2 + o2 <=> hco + oh */
    kc[68] = exp((g_RT[4] + g_RT[10]) - (g_RT[7] + g_RT[12]));

    /*reaction 70: ch2 + ch3 <=> c2h4 + h */
    kc[69] = exp((g_RT[4] + g_RT[3]) - (g_RT[20] + g_RT[1]));

    /*reaction 71: ch2 + ch2 <=> c2h2 + h + h */
    kc[70] = refC * exp((g_RT[4] + g_RT[4]) - (g_RT[18] + g_RT[1] + g_RT[1]));

    /*reaction 72: ch2 + hcco <=> c2h3 + co */
    kc[71] = exp((g_RT[4] + g_RT[17]) - (g_RT[19] + g_RT[9]));

    /*reaction 73: ch2 + c2h2 <=> h2ccch + h */
    kc[72] = exp((g_RT[4] + g_RT[18]) - (g_RT[26] + g_RT[1]));

    /*reaction 74: ch2s + M <=> ch2 + M */
    kc[73] = exp((g_RT[28]) - (g_RT[4]));

    /*reaction 75: ch2s + ch4 <=> ch3 + ch3 */
    kc[74] = exp((g_RT[28] + g_RT[2]) - (g_RT[3] + g_RT[3]));

    /*reaction 76: ch2s + c2h6 <=> ch3 + c2h5 */
    kc[75] = exp((g_RT[28] + g_RT[22]) - (g_RT[3] + g_RT[21]));

    /*reaction 77: ch2s + o2 <=> co + oh + h */
    kc[76] = refC * exp((g_RT[28] + g_RT[10]) - (g_RT[9] + g_RT[12] + g_RT[1]));

    /*reaction 78: ch2s + h2 <=> ch3 + h */
    kc[77] = exp((g_RT[28] + g_RT[0]) - (g_RT[3] + g_RT[1]));

    /*reaction 79: ch2s + c2h2 <=> h2ccch + h */
    kc[78] = exp((g_RT[28] + g_RT[18]) - (g_RT[26] + g_RT[1]));

    /*reaction 80: ch2s + o <=> co + h + h */
    kc[79] = refC * exp((g_RT[28] + g_RT[11]) - (g_RT[9] + g_RT[1] + g_RT[1]));

    /*reaction 81: ch2s + oh <=> ch2o + h */
    kc[80] = exp((g_RT[28] + g_RT[12]) - (g_RT[6] + g_RT[1]));

    /*reaction 82: ch2s + h <=> ch + h2 */
    kc[81] = exp((g_RT[28] + g_RT[1]) - (g_RT[5] + g_RT[0]));

    /*reaction 83: ch2s + co2 <=> ch2o + co */
    kc[82] = exp((g_RT[28] + g_RT[8]) - (g_RT[6] + g_RT[9]));

    /*reaction 84: ch2s + ch3 <=> c2h4 + h */
    kc[83] = exp((g_RT[28] + g_RT[3]) - (g_RT[20] + g_RT[1]));

    /*reaction 85: ch2s + ch2co <=> c2h4 + co */
    kc[84] = exp((g_RT[28] + g_RT[29]) - (g_RT[20] + g_RT[9]));

    /*reaction 86: ch + o2 <=> hco + o */
    kc[85] = exp((g_RT[5] + g_RT[10]) - (g_RT[7] + g_RT[11]));

    /*reaction 87: ch + o <=> co + h */
    kc[86] = exp((g_RT[5] + g_RT[11]) - (g_RT[9] + g_RT[1]));

    /*reaction 88: ch + oh <=> hco + h */
    kc[87] = exp((g_RT[5] + g_RT[12]) - (g_RT[7] + g_RT[1]));

    /*reaction 89: ch + co2 <=> hco + co */
    kc[88] = exp((g_RT[5] + g_RT[8]) - (g_RT[7] + g_RT[9]));

    /*reaction 90: ch + h2o <=> ch2o + h */
    kc[89] = exp((g_RT[5] + g_RT[15]) - (g_RT[6] + g_RT[1]));

    /*reaction 91: ch + ch2o <=> ch2co + h */
    kc[90] = exp((g_RT[5] + g_RT[6]) - (g_RT[29] + g_RT[1]));

    /*reaction 92: ch + c2h2 <=> c3h2 + h */
    kc[91] = exp((g_RT[5] + g_RT[18]) - (g_RT[27] + g_RT[1]));

    /*reaction 93: ch + ch2 <=> c2h2 + h */
    kc[92] = exp((g_RT[5] + g_RT[4]) - (g_RT[18] + g_RT[1]));

    /*reaction 94: ch + ch3 <=> c2h3 + h */
    kc[93] = exp((g_RT[5] + g_RT[3]) - (g_RT[19] + g_RT[1]));

    /*reaction 95: ch + ch4 <=> c2h4 + h */
    kc[94] = exp((g_RT[5] + g_RT[2]) - (g_RT[20] + g_RT[1]));

    /*reaction 96: ch2o + oh <=> hco + h2o */
    kc[95] = exp((g_RT[6] + g_RT[12]) - (g_RT[7] + g_RT[15]));

    /*reaction 97: ch2o + h <=> hco + h2 */
    kc[96] = exp((g_RT[6] + g_RT[1]) - (g_RT[7] + g_RT[0]));

    /*reaction 98: ch2o + M <=> hco + h + M */
    kc[97] = refC * exp((g_RT[6]) - (g_RT[7] + g_RT[1]));

    /*reaction 99: ch2o + o <=> hco + oh */
    kc[98] = exp((g_RT[6] + g_RT[11]) - (g_RT[7] + g_RT[12]));

    /*reaction 100: hco + o2 <=> co + ho2 */
    kc[99] = exp((g_RT[7] + g_RT[10]) - (g_RT[9] + g_RT[13]));

    /*reaction 101: hco + M <=> h + co + M */
    kc[100] = refC * exp((g_RT[7]) - (g_RT[1] + g_RT[9]));

    /*reaction 102: hco + oh <=> h2o + co */
    kc[101] = exp((g_RT[7] + g_RT[12]) - (g_RT[15] + g_RT[9]));

    /*reaction 103: hco + h <=> co + h2 */
    kc[102] = exp((g_RT[7] + g_RT[1]) - (g_RT[9] + g_RT[0]));

    /*reaction 104: hco + o <=> co + oh */
    kc[103] = exp((g_RT[7] + g_RT[11]) - (g_RT[9] + g_RT[12]));

    /*reaction 105: hco + o <=> co2 + h */
    kc[104] = exp((g_RT[7] + g_RT[11]) - (g_RT[8] + g_RT[1]));

    /*reaction 106: co + oh <=> co2 + h */
    kc[105] = exp((g_RT[9] + g_RT[12]) - (g_RT[8] + g_RT[1]));

    /*reaction 107: co + o + M <=> co2 + M */
    kc[106] = 1.0 / (refC) * exp((g_RT[9] + g_RT[11]) - (g_RT[8]));

    /*reaction 108: co + o2 <=> co2 + o */
    kc[107] = exp((g_RT[9] + g_RT[10]) - (g_RT[8] + g_RT[11]));

    /*reaction 109: co + ho2 <=> co2 + oh */
    kc[108] = exp((g_RT[9] + g_RT[13]) - (g_RT[8] + g_RT[12]));

    /*reaction 110: c2h5oh (+M) <=> ch3 + ch2oh (+M) */
    kc[109] = refC * exp((g_RT[33]) - (g_RT[3] + g_RT[23]));

    /*reaction 111: c2h5oh (+M) <=> c2h5 + oh (+M) */
    kc[110] = refC * exp((g_RT[33]) - (g_RT[21] + g_RT[12]));

    /*reaction 112: c2h5oh (+M) <=> c2h4 + h2o (+M) */
    kc[111] = refC * exp((g_RT[33]) - (g_RT[20] + g_RT[15]));

    /*reaction 113: c2h5oh (+M) <=> ch3hco + h2 (+M) */
    kc[112] = refC * exp((g_RT[33]) - (g_RT[32] + g_RT[0]));

    /*reaction 114: c2h5oh + oh <=> c2h4oh + h2o */
    kc[113] = exp((g_RT[33] + g_RT[12]) - (g_RT[34] + g_RT[15]));

    /*reaction 115: c2h5oh + oh <=> ch3choh + h2o */
    kc[114] = exp((g_RT[33] + g_RT[12]) - (g_RT[35] + g_RT[15]));

    /*reaction 116: c2h5oh + oh <=> ch3ch2o + h2o */
    kc[115] = exp((g_RT[33] + g_RT[12]) - (g_RT[36] + g_RT[15]));

    /*reaction 117: c2h5oh + h <=> c2h4oh + h2 */
    kc[116] = exp((g_RT[33] + g_RT[1]) - (g_RT[34] + g_RT[0]));

    /*reaction 118: c2h5oh + h <=> ch3choh + h2 */
    kc[117] = exp((g_RT[33] + g_RT[1]) - (g_RT[35] + g_RT[0]));

    /*reaction 119: c2h5oh + h <=> ch3ch2o + h2 */
    kc[118] = exp((g_RT[33] + g_RT[1]) - (g_RT[36] + g_RT[0]));

    /*reaction 120: c2h5oh + o <=> c2h4oh + oh */
    kc[119] = exp((g_RT[33] + g_RT[11]) - (g_RT[34] + g_RT[12]));

    /*reaction 121: c2h5oh + o <=> ch3choh + oh */
    kc[120] = exp((g_RT[33] + g_RT[11]) - (g_RT[35] + g_RT[12]));

    /*reaction 122: c2h5oh + o <=> ch3ch2o + oh */
    kc[121] = exp((g_RT[33] + g_RT[11]) - (g_RT[36] + g_RT[12]));

    /*reaction 123: c2h5oh + ch3 <=> c2h4oh + ch4 */
    kc[122] = exp((g_RT[33] + g_RT[3]) - (g_RT[34] + g_RT[2]));

    /*reaction 124: c2h5oh + ch3 <=> ch3choh + ch4 */
    kc[123] = exp((g_RT[33] + g_RT[3]) - (g_RT[35] + g_RT[2]));

    /*reaction 125: c2h5oh + ch3 <=> ch3ch2o + ch4 */
    kc[124] = exp((g_RT[33] + g_RT[3]) - (g_RT[36] + g_RT[2]));

    /*reaction 126: c2h5oh + ho2 <=> ch3choh + h2o2 */
    kc[125] = exp((g_RT[33] + g_RT[13]) - (g_RT[35] + g_RT[14]));

    /*reaction 127: c2h5oh + ho2 <=> c2h4oh + h2o2 */
    kc[126] = exp((g_RT[33] + g_RT[13]) - (g_RT[34] + g_RT[14]));

    /*reaction 128: c2h5oh + ho2 <=> ch3ch2o + h2o2 */
    kc[127] = exp((g_RT[33] + g_RT[13]) - (g_RT[36] + g_RT[14]));

    /*reaction 129: ch3ch2o + M <=> ch3hco + h + M */
    kc[128] = refC * exp((g_RT[36]) - (g_RT[32] + g_RT[1]));

    /*reaction 130: ch3ch2o + M <=> ch3 + ch2o + M */
    kc[129] = refC * exp((g_RT[36]) - (g_RT[3] + g_RT[6]));

    /*reaction 131: ch3ch2o + o2 <=> ch3hco + ho2 */
    kc[130] = exp((g_RT[36] + g_RT[10]) - (g_RT[32] + g_RT[13]));

    /*reaction 132: ch3ch2o + co <=> c2h5 + co2 */
    kc[131] = exp((g_RT[36] + g_RT[9]) - (g_RT[21] + g_RT[8]));

    /*reaction 133: ch3ch2o + h <=> ch3 + ch2oh */
    kc[132] = exp((g_RT[36] + g_RT[1]) - (g_RT[3] + g_RT[23]));

    /*reaction 134: ch3ch2o + h <=> c2h4 + h2o */
    kc[133] = exp((g_RT[36] + g_RT[1]) - (g_RT[20] + g_RT[15]));

    /*reaction 135: ch3ch2o + oh <=> ch3hco + h2o */
    kc[134] = exp((g_RT[36] + g_RT[12]) - (g_RT[32] + g_RT[15]));

    /*reaction 136: ch3choh + o2 <=> ch3hco + ho2 */
    kc[135] = exp((g_RT[35] + g_RT[10]) - (g_RT[32] + g_RT[13]));

    /*reaction 137: ch3choh + o2 <=> ch3hco + ho2 */
    kc[136] = exp((g_RT[35] + g_RT[10]) - (g_RT[32] + g_RT[13]));

    /*reaction 138: ch3choh + o <=> ch3hco + oh */
    kc[137] = exp((g_RT[35] + g_RT[11]) - (g_RT[32] + g_RT[12]));

    /*reaction 139: ch3choh + h <=> c2h4 + h2o */
    kc[138] = exp((g_RT[35] + g_RT[1]) - (g_RT[20] + g_RT[15]));

    /*reaction 140: ch3choh + h <=> ch3 + ch2oh */
    kc[139] = exp((g_RT[35] + g_RT[1]) - (g_RT[3] + g_RT[23]));

    /*reaction 141: ch3choh + ho2 <=> ch3hco + oh + oh */
    kc[140] = refC * exp((g_RT[35] + g_RT[13]) - (g_RT[32] + g_RT[12] + g_RT[12]));

    /*reaction 142: ch3choh + oh <=> ch3hco + h2o */
    kc[141] = exp((g_RT[35] + g_RT[12]) - (g_RT[32] + g_RT[15]));

    /*reaction 143: ch3choh + M <=> ch3hco + h + M */
    kc[142] = refC * exp((g_RT[35]) - (g_RT[32] + g_RT[1]));

    /*reaction 144: ch3hco + oh <=> ch3co + h2o */
    kc[143] = exp((g_RT[32] + g_RT[12]) - (g_RT[31] + g_RT[15]));

    /*reaction 145: ch3hco + oh <=> ch2hco + h2o */
    kc[144] = exp((g_RT[32] + g_RT[12]) - (g_RT[30] + g_RT[15]));

    /*reaction 146: ch3hco + o <=> ch3co + oh */
    kc[145] = exp((g_RT[32] + g_RT[11]) - (g_RT[31] + g_RT[12]));

    /*reaction 147: ch3hco + o <=> ch2hco + oh */
    kc[146] = exp((g_RT[32] + g_RT[11]) - (g_RT[30] + g_RT[12]));

    /*reaction 148: ch3hco + h <=> ch3co + h2 */
    kc[147] = exp((g_RT[32] + g_RT[1]) - (g_RT[31] + g_RT[0]));

    /*reaction 149: ch3hco + h <=> ch2hco + h2 */
    kc[148] = exp((g_RT[32] + g_RT[1]) - (g_RT[30] + g_RT[0]));

    /*reaction 150: ch3hco + ch3 <=> ch3co + ch4 */
    kc[149] = exp((g_RT[32] + g_RT[3]) - (g_RT[31] + g_RT[2]));

    /*reaction 151: ch3hco + ch3 <=> ch2hco + ch4 */
    kc[150] = exp((g_RT[32] + g_RT[3]) - (g_RT[30] + g_RT[2]));

    /*reaction 152: ch3hco + ho2 <=> ch3co + h2o2 */
    kc[151] = exp((g_RT[32] + g_RT[13]) - (g_RT[31] + g_RT[14]));

    /*reaction 153: ch3hco + ho2 <=> ch2hco + h2o2 */
    kc[152] = exp((g_RT[32] + g_RT[13]) - (g_RT[30] + g_RT[14]));

    /*reaction 154: ch3hco + o2 <=> ch3co + ho2 */
    kc[153] = exp((g_RT[32] + g_RT[10]) - (g_RT[31] + g_RT[13]));

    /*reaction 155: c2h6 + ch3 <=> c2h5 + ch4 */
    kc[154] = exp((g_RT[22] + g_RT[3]) - (g_RT[21] + g_RT[2]));

    /*reaction 156: c2h6 + h <=> c2h5 + h2 */
    kc[155] = exp((g_RT[22] + g_RT[1]) - (g_RT[21] + g_RT[0]));

    /*reaction 157: c2h6 + o <=> c2h5 + oh */
    kc[156] = exp((g_RT[22] + g_RT[11]) - (g_RT[21] + g_RT[12]));

    /*reaction 158: c2h6 + oh <=> c2h5 + h2o */
    kc[157] = exp((g_RT[22] + g_RT[12]) - (g_RT[21] + g_RT[15]));

    /*reaction 159: c2h5 + h <=> c2h4 + h2 */
    kc[158] = exp((g_RT[21] + g_RT[1]) - (g_RT[20] + g_RT[0]));

    /*reaction 160: c2h5 + h <=> ch3 + ch3 */
    kc[159] = exp((g_RT[21] + g_RT[1]) - (g_RT[3] + g_RT[3]));

    /*reaction 161: c2h5 + h <=> c2h6 */
    kc[160] = 1.0 / (refC) * exp((g_RT[21] + g_RT[1]) - (g_RT[22]));

    /*reaction 162: c2h5 + oh <=> c2h4 + h2o */
    kc[161] = exp((g_RT[21] + g_RT[12]) - (g_RT[20] + g_RT[15]));

    /*reaction 163: c2h5 + o <=> ch3 + ch2o */
    kc[162] = exp((g_RT[21] + g_RT[11]) - (g_RT[3] + g_RT[6]));

    /*reaction 164: c2h5 + ho2 <=> c2h6 + o2 */
    kc[163] = exp((g_RT[21] + g_RT[13]) - (g_RT[22] + g_RT[10]));

    /*reaction 165: c2h5 + ho2 <=> ch3ch2o + oh */
    kc[164] = exp((g_RT[21] + g_RT[13]) - (g_RT[36] + g_RT[12]));

    /*reaction 166: c2h5 + o2 <=> c2h4 + ho2 */
    kc[165] = exp((g_RT[21] + g_RT[10]) - (g_RT[20] + g_RT[13]));

    /*reaction 167: c2h5 + o2 <=> ch3hco + oh */
    kc[166] = exp((g_RT[21] + g_RT[10]) - (g_RT[32] + g_RT[12]));

    /*reaction 168: c2h4 + oh <=> c2h4oh */
    kc[167] = 1.0 / (refC) * exp((g_RT[20] + g_RT[12]) - (g_RT[34]));

    /*reaction 169: c2h4 + oh <=> c2h3 + h2o */
    kc[168] = exp((g_RT[20] + g_RT[12]) - (g_RT[19] + g_RT[15]));

    /*reaction 170: c2h4 + o <=> ch3 + hco */
    kc[169] = exp((g_RT[20] + g_RT[11]) - (g_RT[3] + g_RT[7]));

    /*reaction 171: c2h4 + o <=> ch2hco + h */
    kc[170] = exp((g_RT[20] + g_RT[11]) - (g_RT[30] + g_RT[1]));

    /*reaction 172: c2h4 + ch3 <=> c2h3 + ch4 */
    kc[171] = exp((g_RT[20] + g_RT[3]) - (g_RT[19] + g_RT[2]));

    /*reaction 173: c2h4 + h <=> c2h3 + h2 */
    kc[172] = exp((g_RT[20] + g_RT[1]) - (g_RT[19] + g_RT[0]));

    /*reaction 174: c2h4 + h (+M) <=> c2h5 (+M) */
    kc[173] = 1.0 / (refC) * exp((g_RT[20] + g_RT[1]) - (g_RT[21]));

    /*reaction 175: c2h4 (+M) <=> c2h2 + h2 (+M) */
    kc[174] = refC * exp((g_RT[20]) - (g_RT[18] + g_RT[0]));

    /*reaction 176: c2h3 + h (+M) <=> c2h4 (+M) */
    kc[175] = 1.0 / (refC) * exp((g_RT[19] + g_RT[1]) - (g_RT[20]));

    /*reaction 177: c2h3 + h <=> c2h2 + h2 */
    kc[176] = exp((g_RT[19] + g_RT[1]) - (g_RT[18] + g_RT[0]));

    /*reaction 178: c2h3 + o <=> ch2co + h */
    kc[177] = exp((g_RT[19] + g_RT[11]) - (g_RT[29] + g_RT[1]));

    /*reaction 179: c2h3 + o2 <=> ch2o + hco */
    kc[178] = exp((g_RT[19] + g_RT[10]) - (g_RT[6] + g_RT[7]));

    /*reaction 180: c2h3 + o2 <=> ch2hco + o */
    kc[179] = exp((g_RT[19] + g_RT[10]) - (g_RT[30] + g_RT[11]));

    /*reaction 181: c2h3 + o2 <=> c2h2 + ho2 */
    kc[180] = exp((g_RT[19] + g_RT[10]) - (g_RT[18] + g_RT[13]));

    /*reaction 182: c2h3 + oh <=> c2h2 + h2o */
    kc[181] = exp((g_RT[19] + g_RT[12]) - (g_RT[18] + g_RT[15]));

    /*reaction 183: c2h3 + c2h <=> c2h2 + c2h2 */
    kc[182] = exp((g_RT[19] + g_RT[16]) - (g_RT[18] + g_RT[18]));

    /*reaction 184: c2h3 + ch <=> ch2 + c2h2 */
    kc[183] = exp((g_RT[19] + g_RT[5]) - (g_RT[4] + g_RT[18]));

    /*reaction 185: c2h3 + ch3 <=> c2h2 + ch4 */
    kc[184] = exp((g_RT[19] + g_RT[3]) - (g_RT[18] + g_RT[2]));

    /*reaction 186: c2h2 + oh <=> c2h + h2o */
    kc[185] = exp((g_RT[18] + g_RT[12]) - (g_RT[16] + g_RT[15]));

    /*reaction 187: c2h2 + oh <=> hccoh + h */
    kc[186] = exp((g_RT[18] + g_RT[12]) - (g_RT[25] + g_RT[1]));

    /*reaction 188: c2h2 + oh <=> ch2co + h */
    kc[187] = exp((g_RT[18] + g_RT[12]) - (g_RT[29] + g_RT[1]));

    /*reaction 189: c2h2 + oh <=> ch2co + h */
    kc[188] = exp((g_RT[18] + g_RT[12]) - (g_RT[29] + g_RT[1]));

    /*reaction 190: c2h2 + oh <=> ch3 + co */
    kc[189] = exp((g_RT[18] + g_RT[12]) - (g_RT[3] + g_RT[9]));

    /*reaction 191: hccoh + h <=> ch2co + h */
    kc[190] = exp((g_RT[25] + g_RT[1]) - (g_RT[29] + g_RT[1]));

    /*reaction 192: c2h2 + o <=> ch2 + co */
    kc[191] = exp((g_RT[18] + g_RT[11]) - (g_RT[4] + g_RT[9]));

    /*reaction 193: c2h2 + o <=> hcco + h */
    kc[192] = exp((g_RT[18] + g_RT[11]) - (g_RT[17] + g_RT[1]));

    /*reaction 194: c2h2 + o <=> c2h + oh */
    kc[193] = exp((g_RT[18] + g_RT[11]) - (g_RT[16] + g_RT[12]));

    /*reaction 195: c2h2 + ch3 <=> c2h + ch4 */
    kc[194] = exp((g_RT[18] + g_RT[3]) - (g_RT[16] + g_RT[2]));

    /*reaction 196: c2h2 + o2 <=> hcco + oh */
    kc[195] = exp((g_RT[18] + g_RT[10]) - (g_RT[17] + g_RT[12]));

    /*reaction 197: c2h2 + M <=> c2h + h + M */
    kc[196] = refC * exp((g_RT[18]) - (g_RT[16] + g_RT[1]));

    /*reaction 198: c2h2 + h (+M) <=> c2h3 (+M) */
    kc[197] = 1.0 / (refC) * exp((g_RT[18] + g_RT[1]) - (g_RT[19]));

    /*reaction 199: ch2hco + h <=> ch3 + hco */
    kc[198] = exp((g_RT[30] + g_RT[1]) - (g_RT[3] + g_RT[7]));

    /*reaction 200: ch2hco + h <=> ch2co + h2 */
    kc[199] = exp((g_RT[30] + g_RT[1]) - (g_RT[29] + g_RT[0]));

    /*reaction 201: ch2hco + o <=> ch2o + hco */
    kc[200] = exp((g_RT[30] + g_RT[11]) - (g_RT[6] + g_RT[7]));

    /*reaction 202: ch2hco + oh <=> ch2co + h2o */
    kc[201] = exp((g_RT[30] + g_RT[12]) - (g_RT[29] + g_RT[15]));

    /*reaction 203: ch2hco + o2 <=> ch2o + co + oh */
    kc[202] = refC * exp((g_RT[30] + g_RT[10]) - (g_RT[6] + g_RT[9] + g_RT[12]));

    /*reaction 204: ch2hco + ch3 <=> c2h5 + co + h */
    kc[203] = refC * exp((g_RT[30] + g_RT[3]) - (g_RT[21] + g_RT[9] + g_RT[1]));

    /*reaction 205: ch2hco + ho2 <=> ch2o + hco + oh */
    kc[204] = refC * exp((g_RT[30] + g_RT[13]) - (g_RT[6] + g_RT[7] + g_RT[12]));

    /*reaction 206: ch2hco + ho2 <=> ch3hco + o2 */
    kc[205] = exp((g_RT[30] + g_RT[13]) - (g_RT[32] + g_RT[10]));

    /*reaction 207: ch2hco <=> ch3 + co */
    kc[206] = refC * exp((g_RT[30]) - (g_RT[3] + g_RT[9]));

    /*reaction 208: ch2hco <=> ch2co + h */
    kc[207] = refC * exp((g_RT[30]) - (g_RT[29] + g_RT[1]));

    /*reaction 209: ch3co (+M) <=> ch3 + co (+M) */
    kc[208] = refC * exp((g_RT[31]) - (g_RT[3] + g_RT[9]));

    /*reaction 210: ch2co + o <=> co2 + ch2 */
    kc[209] = exp((g_RT[29] + g_RT[11]) - (g_RT[8] + g_RT[4]));

    /*reaction 211: ch2co + h <=> ch3 + co */
    kc[210] = exp((g_RT[29] + g_RT[1]) - (g_RT[3] + g_RT[9]));

    /*reaction 212: ch2co + h <=> hcco + h2 */
    kc[211] = exp((g_RT[29] + g_RT[1]) - (g_RT[17] + g_RT[0]));

    /*reaction 213: ch2co + o <=> hcco + oh */
    kc[212] = exp((g_RT[29] + g_RT[11]) - (g_RT[17] + g_RT[12]));

    /*reaction 214: ch2co + oh <=> hcco + h2o */
    kc[213] = exp((g_RT[29] + g_RT[12]) - (g_RT[17] + g_RT[15]));

    /*reaction 215: ch2co + oh <=> ch2oh + co */
    kc[214] = exp((g_RT[29] + g_RT[12]) - (g_RT[23] + g_RT[9]));

    /*reaction 216: ch2co (+M) <=> ch2 + co (+M) */
    kc[215] = refC * exp((g_RT[29]) - (g_RT[4] + g_RT[9]));

    /*reaction 217: c2h + h2 <=> c2h2 + h */
    kc[216] = exp((g_RT[16] + g_RT[0]) - (g_RT[18] + g_RT[1]));

    /*reaction 218: c2h + o <=> ch + co */
    kc[217] = exp((g_RT[16] + g_RT[11]) - (g_RT[5] + g_RT[9]));

    /*reaction 219: c2h + oh <=> hcco + h */
    kc[218] = exp((g_RT[16] + g_RT[12]) - (g_RT[17] + g_RT[1]));

    /*reaction 220: c2h + o2 <=> co + co + h */
    kc[219] = refC * exp((g_RT[16] + g_RT[10]) - (g_RT[9] + g_RT[9] + g_RT[1]));

    /*reaction 221: hcco + c2h2 <=> h2ccch + co */
    kc[220] = exp((g_RT[17] + g_RT[18]) - (g_RT[26] + g_RT[9]));

    /*reaction 222: hcco + h <=> ch2s + co */
    kc[221] = exp((g_RT[17] + g_RT[1]) - (g_RT[28] + g_RT[9]));

    /*reaction 223: hcco + o <=> h + co + co */
    kc[222] = refC * exp((g_RT[17] + g_RT[11]) - (g_RT[1] + g_RT[9] + g_RT[9]));

    /*reaction 224: hcco + o <=> ch + co2 */
    kc[223] = exp((g_RT[17] + g_RT[11]) - (g_RT[5] + g_RT[8]));

    /*reaction 225: hcco + o2 <=> hco + co + o */
    kc[224] = refC * exp((g_RT[17] + g_RT[10]) - (g_RT[7] + g_RT[9] + g_RT[11]));

    /*reaction 226: hcco + o2 <=> co2 + hco */
    kc[225] = exp((g_RT[17] + g_RT[10]) - (g_RT[8] + g_RT[7]));

    /*reaction 227: hcco + ch <=> c2h2 + co */
    kc[226] = exp((g_RT[17] + g_RT[5]) - (g_RT[18] + g_RT[9]));

    /*reaction 228: hcco + hcco <=> c2h2 + co + co */
    kc[227] = refC * exp((g_RT[17] + g_RT[17]) - (g_RT[18] + g_RT[9] + g_RT[9]));

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void gibbs(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1410 kelvin */
    if (T < 1410) {
        /*species 23: ch2oh */
        species[23] =
            -2.334787240000000e+03 / tc[1]
            -8.726552209999999e+00
            -2.600678490000000e+00 * tc[0]
            -6.422392700000000e-03 * tc[1]
            +1.389660308333334e-06 * tc[2]
            -2.297730050000000e-10 * tc[3]
            +1.785205530000000e-14 * tc[4];
    } else {
        /*species 23: ch2oh */
        species[23] =
            -3.501570980000000e+03 / tc[1]
            +1.292964647000000e+01
            -6.001278030000000e+00 * tc[0]
            -2.493607840000000e-03 * tc[1]
            +2.682558583333334e-07 * tc[2]
            -2.001892791666667e-11 * tc[3]
            +6.779135000000000e-16 * tc[4];
    }

    /*species with midpoint at T=1350 kelvin */
    if (T < 1350) {
        /*species 4: ch2 */
        species[4] =
            +4.591927710000000e+04 / tc[1]
            +5.465963500000002e-01
            -3.496239890000000e+00 * tc[0]
            -1.101293225000000e-03 * tc[1]
            +7.063813316666667e-08 * tc[2]
            -1.137372550000000e-12 * tc[3]
            -2.159686860000000e-16 * tc[4];
    } else {
        /*species 4: ch2 */
        species[4] =
            +4.593757430000000e+04 / tc[1]
            +2.082575200000001e-01
            -3.443104660000000e+00 * tc[0]
            -1.163190055000000e-03 * tc[1]
            +8.885799966666666e-08 * tc[2]
            -4.783699733333333e-12 * tc[3]
            +1.213474245000000e-16 * tc[4];
    }

    /*species with midpoint at T=1671 kelvin */
    if (T < 1671) {
        /*species 19: c2h3 */
        species[19] =
            +3.428689790000000e+04 / tc[1]
            -7.413894080000000e+00
            -2.739259420000000e+00 * tc[0]
            -3.515057955000000e-03 * tc[1]
            -3.944121650000000e-07 * tc[2]
            +2.996414475000000e-10 * tc[3]
            -4.458782455000000e-14 * tc[4];
    } else {
        /*species 19: c2h3 */
        species[19] =
            +3.351535440000000e+04 / tc[1]
            +1.703842990000000e+00
            -3.960477130000000e+00 * tc[0]
            -3.997130065000000e-03 * tc[1]
            +4.760134783333333e-07 * tc[2]
            -3.819590108333334e-11 * tc[3]
            +1.362850700000000e-15 * tc[4];
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: h2 */
        species[0] =
            -1.012521000000000e+03 / tc[1]
            +6.592218000000000e+00
            -3.298124000000000e+00 * tc[0]
            -4.124721000000000e-04 * tc[1]
            +1.357169166666667e-07 * tc[2]
            +7.896194999999999e-12 * tc[3]
            -2.067436000000000e-14 * tc[4];
        /*species 1: h */
        species[1] =
            +2.547163000000000e+04 / tc[1]
            +2.960117600000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 2: ch4 */
        species[2] =
            -9.825228999999999e+03 / tc[1]
            -1.294344850000000e+01
            -7.787415000000000e-01 * tc[0]
            -8.738340000000001e-03 * tc[1]
            +4.639015000000000e-06 * tc[2]
            -2.541423333333333e-09 * tc[3]
            +6.119655000000000e-13 * tc[4];
        /*species 3: ch3 */
        species[3] =
            +1.642378000000000e+04 / tc[1]
            -4.359351000000000e+00
            -2.430443000000000e+00 * tc[0]
            -5.562050000000000e-03 * tc[1]
            +2.800366666666666e-06 * tc[2]
            -1.351524166666667e-09 * tc[3]
            +2.932476500000000e-13 * tc[4];
        /*species 5: ch */
        species[5] =
            +7.045259000000000e+04 / tc[1]
            -1.313860000000000e-01
            -3.200202000000000e+00 * tc[0]
            -1.036438000000000e-03 * tc[1]
            +8.557385000000000e-07 * tc[2]
            -4.778241666666666e-10 * tc[3]
            +9.777665000000000e-14 * tc[4];
        /*species 6: ch2o */
        species[6] =
            -1.486540000000000e+04 / tc[1]
            -1.213208900000000e+01
            -1.652731000000000e+00 * tc[0]
            -6.315720000000000e-03 * tc[1]
            +3.146946666666667e-06 * tc[2]
            -1.708359166666667e-09 * tc[3]
            +4.206618500000000e-13 * tc[4];
        /*species 7: hco */
        species[7] =
            +4.159922000000000e+03 / tc[1]
            -6.085284000000000e+00
            -2.898330000000000e+00 * tc[0]
            -3.099573500000000e-03 * tc[1]
            +1.603847333333333e-06 * tc[2]
            -9.081875000000000e-10 * tc[3]
            +2.287442500000000e-13 * tc[4];
        /*species 8: co2 */
        species[8] =
            -4.837314000000000e+04 / tc[1]
            -7.912765000000000e+00
            -2.275725000000000e+00 * tc[0]
            -4.961036000000000e-03 * tc[1]
            +1.734851666666667e-06 * tc[2]
            -5.722239166666667e-10 * tc[3]
            +1.058640000000000e-13 * tc[4];
        /*species 9: co */
        species[9] =
            -1.431054000000000e+04 / tc[1]
            -1.586445000000000e+00
            -3.262452000000000e+00 * tc[0]
            -7.559705000000000e-04 * tc[1]
            +6.469591666666667e-07 * tc[2]
            -4.651620000000000e-10 * tc[3]
            +1.237475500000000e-13 * tc[4];
        /*species 10: o2 */
        species[10] =
            -1.005249000000000e+03 / tc[1]
            -2.821802000000000e+00
            -3.212936000000000e+00 * tc[0]
            -5.637430000000000e-04 * tc[1]
            +9.593583333333333e-08 * tc[2]
            -1.094897500000000e-10 * tc[3]
            +4.384277000000000e-14 * tc[4];
        /*species 11: o */
        species[11] =
            +2.914764000000000e+04 / tc[1]
            -1.756599999999997e-02
            -2.946429000000000e+00 * tc[0]
            +8.190830000000000e-04 * tc[1]
            -4.035053333333333e-07 * tc[2]
            +1.335702500000000e-10 * tc[3]
            -1.945348000000000e-14 * tc[4];
        /*species 12: oh */
        species[12] =
            +3.606782000000000e+03 / tc[1]
            +2.278406000000000e+00
            -3.637266000000000e+00 * tc[0]
            -9.254550000000000e-05 * tc[1]
            +2.793608333333333e-07 * tc[2]
            -1.989335833333333e-10 * tc[3]
            +4.215721000000000e-14 * tc[4];
        /*species 14: h2o2 */
        species[14] =
            -1.766315000000000e+04 / tc[1]
            -3.396609000000000e+00
            -3.388754000000000e+00 * tc[0]
            -3.284613000000000e-03 * tc[1]
            +2.475021666666666e-08 * tc[2]
            +3.854838333333333e-10 * tc[3]
            -1.235757500000000e-13 * tc[4];
        /*species 15: h2o */
        species[15] =
            -3.020811000000000e+04 / tc[1]
            +7.966090000000001e-01
            -3.386842000000000e+00 * tc[0]
            -1.737491000000000e-03 * tc[1]
            +1.059116000000000e-06 * tc[2]
            -5.807150833333333e-10 * tc[3]
            +1.253294000000000e-13 * tc[4];
        /*species 16: c2h */
        species[16] =
            +6.683813000000000e+04 / tc[1]
            -4.562516000000000e+00
            -2.737704000000000e+00 * tc[0]
            -4.024223000000000e-03 * tc[1]
            +1.540718333333333e-06 * tc[2]
            -5.437715833333333e-10 * tc[3]
            +9.697900000000001e-14 * tc[4];
        /*species 17: hcco */
        species[17] =
            +1.965892000000000e+04 / tc[1]
            +4.566121099999999e+00
            -5.047965000000000e+00 * tc[0]
            -2.226739000000000e-03 * tc[1]
            -3.780471666666667e-08 * tc[2]
            +1.235079166666667e-10 * tc[3]
            -1.125371000000000e-14 * tc[4];
        /*species 18: c2h2 */
        species[18] =
            +2.612444000000000e+04 / tc[1]
            -6.791815999999999e+00
            -2.013562000000000e+00 * tc[0]
            -7.595225000000000e-03 * tc[1]
            +2.693865000000000e-06 * tc[2]
            -7.565826666666667e-10 * tc[3]
            +9.563730000000000e-14 * tc[4];
        /*species 20: c2h4 */
        species[20] =
            +5.573046000000000e+03 / tc[1]
            -2.507297800000000e+01
            +8.614880000000000e-01 * tc[0]
            -1.398081500000000e-02 * tc[1]
            +5.647795000000000e-06 * tc[2]
            -2.320960000000000e-09 * tc[3]
            +4.868939500000000e-13 * tc[4];
        /*species 22: c2h6 */
        species[22] =
            -1.123918000000000e+04 / tc[1]
            -1.296975100000000e+01
            -1.462539000000000e+00 * tc[0]
            -7.747335000000000e-03 * tc[1]
            -9.634178333333333e-07 * tc[2]
            +1.048193333333333e-09 * tc[3]
            -2.293133500000000e-13 * tc[4];
        /*species 24: ch3o */
        species[24] =
            +9.786011000000000e+02 / tc[1]
            -1.104597600000000e+01
            -2.106204000000000e+00 * tc[0]
            -3.608297500000000e-03 * tc[1]
            -8.897453333333333e-07 * tc[2]
            +6.148030000000000e-10 * tc[3]
            -1.037805500000000e-13 * tc[4];
        /*species 25: hccoh */
        species[25] =
            +8.701190000000001e+03 / tc[1]
            -5.924100000000001e-01
            -3.899465000000000e+00 * tc[0]
            -4.850537500000000e-03 * tc[1]
            +5.198848333333333e-08 * tc[2]
            +4.614776666666667e-10 * tc[3]
            -1.232866000000000e-13 * tc[4];
        /*species 26: h2ccch */
        species[26] =
            +3.988883000000000e+04 / tc[1]
            +4.168745100000000e+00
            -4.754200000000000e+00 * tc[0]
            -5.540140000000000e-03 * tc[1]
            -4.655538333333333e-08 * tc[2]
            +4.566010000000000e-10 * tc[3]
            -9.748145000000000e-14 * tc[4];
        /*species 27: c3h2 */
        species[27] =
            +6.350421000000000e+04 / tc[1]
            -5.702732000000000e+00
            -3.166714000000000e+00 * tc[0]
            -1.241286000000000e-02 * tc[1]
            +7.652728333333333e-06 * tc[2]
            -3.556682500000000e-09 * tc[3]
            +7.410759999999999e-13 * tc[4];
        /*species 29: ch2co */
        species[29] =
            -7.632637000000000e+03 / tc[1]
            -5.698582000000000e+00
            -2.974971000000000e+00 * tc[0]
            -6.059355000000000e-03 * tc[1]
            +3.908410000000000e-07 * tc[2]
            +5.388904166666666e-10 * tc[3]
            -1.952824500000000e-13 * tc[4];
        /*species 30: ch2hco */
        species[30] =
            +1.521477000000000e+03 / tc[1]
            -6.149227999999999e+00
            -3.409062000000000e+00 * tc[0]
            -5.369285000000000e-03 * tc[1]
            -3.152486666666667e-07 * tc[2]
            +5.965485833333333e-10 * tc[3]
            -1.433692500000000e-13 * tc[4];
        /*species 31: ch3co */
        species[31] =
            -4.108508000000000e+03 / tc[1]
            -8.103572000000000e+00
            -3.125278000000000e+00 * tc[0]
            -4.889110000000000e-03 * tc[1]
            -7.535746666666667e-07 * tc[2]
            +7.507885000000000e-10 * tc[3]
            -1.596859000000000e-13 * tc[4];
        /*species 32: ch3hco */
        species[32] =
            -2.124589000000000e+04 / tc[1]
            -1.084519500000000e+01
            -2.505695000000000e+00 * tc[0]
            -6.684955000000000e-03 * tc[1]
            -7.786588333333333e-07 * tc[2]
            +9.401166666666668e-10 * tc[3]
            -2.131783000000000e-13 * tc[4];
        /*species 33: c2h5oh */
        species[33] =
            -2.999613200000000e+04 / tc[1]
            +5.684120000000004e-02
            -4.858695700000000e+00 * tc[0]
            +1.870086300000000e-03 * tc[1]
            -1.159256300000000e-05 * tc[2]
            +7.387899666666666e-09 * tc[3]
            -1.758441750000000e-12 * tc[4];
        /*species 37: n2 */
        species[37] =
            -1.020900000000000e+03 / tc[1]
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041200000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427500000000e-13 * tc[4];
    } else {
        /*species 0: h2 */
        species[0] =
            -8.350340000000000e+02 / tc[1]
            +4.346533000000000e+00
            -2.991423000000000e+00 * tc[0]
            -3.500322000000000e-04 * tc[1]
            +9.389715000000000e-09 * tc[2]
            +7.692981666666667e-13 * tc[3]
            -7.913760000000000e-17 * tc[4];
        /*species 1: h */
        species[1] =
            +2.547163000000000e+04 / tc[1]
            +2.960117600000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 2: ch4 */
        species[2] =
            -1.008079000000000e+04 / tc[1]
            -7.939916000000000e+00
            -1.683479000000000e+00 * tc[0]
            -5.118620000000000e-03 * tc[1]
            +6.458548333333333e-07 * tc[2]
            -5.654654166666667e-11 * tc[3]
            +2.251711500000000e-15 * tc[4];
        /*species 3: ch3 */
        species[3] =
            +1.643781000000000e+04 / tc[1]
            -2.608645000000000e+00
            -2.844052000000000e+00 * tc[0]
            -3.068987000000000e-03 * tc[1]
            +3.717241666666666e-07 * tc[2]
            -3.154300833333333e-11 * tc[3]
            +1.226079500000000e-15 * tc[4];
        /*species 5: ch */
        species[5] =
            +7.086723000000000e+04 / tc[1]
            -6.982150000000001e+00
            -2.196223000000000e+00 * tc[0]
            -1.170190500000000e-03 * tc[1]
            +1.176366833333333e-07 * tc[2]
            -7.506318333333334e-12 * tc[3]
            +1.927520000000000e-16 * tc[4];
        /*species 6: ch2o */
        species[6] =
            -1.532037000000000e+04 / tc[1]
            -3.916966000000000e+00
            -2.995606000000000e+00 * tc[0]
            -3.340660500000000e-03 * tc[1]
            +4.381591666666666e-07 * tc[2]
            -3.947627500000000e-11 * tc[3]
            +1.606258500000000e-15 * tc[4];
        /*species 7: hco */
        species[7] =
            +3.916324000000000e+03 / tc[1]
            -1.995028000000000e+00
            -3.557271000000000e+00 * tc[0]
            -1.672786500000000e-03 * tc[1]
            +2.225010000000000e-07 * tc[2]
            -2.058810833333333e-11 * tc[3]
            +8.569255000000000e-16 * tc[4];
        /*species 8: co2 */
        species[8] =
            -4.896696000000000e+04 / tc[1]
            +5.409018900000000e+00
            -4.453623000000000e+00 * tc[0]
            -1.570084500000000e-03 * tc[1]
            +2.130685000000000e-07 * tc[2]
            -1.994997500000000e-11 * tc[3]
            +8.345165000000000e-16 * tc[4];
        /*species 9: co */
        species[9] =
            -1.426835000000000e+04 / tc[1]
            -3.083140000000000e+00
            -3.025078000000000e+00 * tc[0]
            -7.213445000000000e-04 * tc[1]
            +9.384713333333334e-08 * tc[2]
            -8.488174999999999e-12 * tc[3]
            +3.455476000000000e-16 * tc[4];
        /*species 10: o2 */
        species[10] =
            -1.233930000000000e+03 / tc[1]
            +5.084119999999999e-01
            -3.697578000000000e+00 * tc[0]
            -3.067598500000000e-04 * tc[1]
            +2.098070000000000e-08 * tc[2]
            -1.479400833333333e-12 * tc[3]
            +5.682175000000001e-17 * tc[4];
        /*species 11: o */
        species[11] =
            +2.923080000000000e+04 / tc[1]
            -2.378248000000000e+00
            -2.542060000000000e+00 * tc[0]
            +1.377531000000000e-05 * tc[1]
            +5.171338333333333e-10 * tc[2]
            -3.792555833333334e-13 * tc[3]
            +2.184026000000000e-17 * tc[4];
        /*species 12: oh */
        species[12] =
            +3.886888000000000e+03 / tc[1]
            -2.712982000000000e+00
            -2.882730000000000e+00 * tc[0]
            -5.069870000000000e-04 * tc[1]
            +3.794795000000000e-08 * tc[2]
            -1.812236666666667e-12 * tc[3]
            +2.563152500000000e-17 * tc[4];
        /*species 14: h2o2 */
        species[14] =
            -1.800696000000000e+04 / tc[1]
            +4.072030000000000e+00
            -4.573167000000000e+00 * tc[0]
            -2.168068000000000e-03 * tc[1]
            +2.457815000000000e-07 * tc[2]
            -1.957420000000000e-11 * tc[3]
            +7.158270000000000e-16 * tc[4];
        /*species 15: h2o */
        species[15] =
            -2.989921000000000e+04 / tc[1]
            -4.190671000000000e+00
            -2.672146000000000e+00 * tc[0]
            -1.528146500000000e-03 * tc[1]
            +1.455043333333333e-07 * tc[2]
            -1.000830000000000e-11 * tc[3]
            +3.195809000000000e-16 * tc[4];
        /*species 16: c2h */
        species[16] =
            +6.655884000000000e+04 / tc[1]
            +2.795304000000000e+00
            -3.986367000000000e+00 * tc[0]
            -1.571561500000000e-03 * tc[1]
            +2.112071666666667e-07 * tc[2]
            -2.436969166666667e-11 * tc[3]
            +1.358160000000000e-15 * tc[4];
        /*species 17: hcco */
        species[17] =
            +1.901513000000000e+04 / tc[1]
            +1.582933500000000e+01
            -6.758073000000000e+00 * tc[0]
            -1.000200000000000e-03 * tc[1]
            +3.379345000000000e-08 * tc[2]
            +8.676100000000000e-12 * tc[3]
            -9.825825000000000e-16 * tc[4];
        /*species 18: c2h2 */
        species[18] =
            +2.566766000000000e+04 / tc[1]
            +7.237108000000000e+00
            -4.436770000000000e+00 * tc[0]
            -2.688019500000000e-03 * tc[1]
            +3.188028333333333e-07 * tc[2]
            -2.738649166666667e-11 * tc[3]
            +1.078355000000000e-15 * tc[4];
        /*species 20: c2h4 */
        species[20] =
            +4.428289000000000e+03 / tc[1]
            +1.298030000000000e+00
            -3.528419000000000e+00 * tc[0]
            -5.742590000000000e-03 * tc[1]
            +7.363975000000000e-07 * tc[2]
            -6.537167500000001e-11 * tc[3]
            +2.633424000000000e-15 * tc[4];
        /*species 22: c2h6 */
        species[22] =
            -1.271779000000000e+04 / tc[1]
            +1.006544500000000e+01
            -4.825938000000000e+00 * tc[0]
            -6.920215000000000e-03 * tc[1]
            +7.595431666666667e-07 * tc[2]
            -5.604139166666666e-11 * tc[3]
            +1.799080500000000e-15 * tc[4];
        /*species 24: ch3o */
        species[24] =
            +1.278325000000000e+02 / tc[1]
            +8.412250000000001e-01
            -3.770800000000000e+00 * tc[0]
            -3.935748500000000e-03 * tc[1]
            +4.427306666666667e-07 * tc[2]
            -3.287025833333333e-11 * tc[3]
            +1.056308000000000e-15 * tc[4];
        /*species 25: hccoh */
        species[25] =
            +7.598258000000000e+03 / tc[1]
            +2.134046400000000e+01
            -7.328324000000000e+00 * tc[0]
            -1.668208000000000e-03 * tc[1]
            +5.041175000000000e-08 * tc[2]
            +1.484255000000000e-11 * tc[3]
            -1.622584000000000e-15 * tc[4];
        /*species 26: h2ccch */
        species[26] =
            +3.847420000000000e+04 / tc[1]
            +3.061023700000000e+01
            -8.831047000000000e+00 * tc[0]
            -2.178597500000000e-03 * tc[1]
            +6.848445000000000e-08 * tc[2]
            +1.973935833333333e-11 * tc[3]
            -2.188260000000000e-15 * tc[4];
        /*species 27: c3h2 */
        species[27] =
            +6.259722000000000e+04 / tc[1]
            +2.003988100000000e+01
            -7.670981000000000e+00 * tc[0]
            -1.374374500000000e-03 * tc[1]
            +7.284905000000000e-08 * tc[2]
            +5.379665833333334e-12 * tc[3]
            -8.319435000000000e-16 * tc[4];
        /*species 29: ch2co */
        species[29] =
            -8.583402000000000e+03 / tc[1]
            +1.369639800000000e+01
            -6.038817000000000e+00 * tc[0]
            -2.902420000000000e-03 * tc[1]
            +3.201590000000000e-07 * tc[2]
            -2.328737500000000e-11 * tc[3]
            +7.294340000000000e-16 * tc[4];
        /*species 30: ch2hco */
        species[30] =
            +4.903218000000000e+02 / tc[1]
            +1.102092100000000e+01
            -5.975670000000000e+00 * tc[0]
            -4.065295500000000e-03 * tc[1]
            +4.572706666666667e-07 * tc[2]
            -3.391920000000000e-11 * tc[3]
            +1.088008500000000e-15 * tc[4];
        /*species 31: ch3co */
        species[31] =
            -5.187863000000000e+03 / tc[1]
            +8.887228000000000e+00
            -5.612279000000000e+00 * tc[0]
            -4.224943000000000e-03 * tc[1]
            +4.756911666666667e-07 * tc[2]
            -3.531980000000000e-11 * tc[3]
            +1.134202000000000e-15 * tc[4];
        /*species 32: ch3hco */
        species[32] =
            -2.264569000000000e+04 / tc[1]
            +1.188159600000000e+01
            -5.868650000000000e+00 * tc[0]
            -5.397120000000000e-03 * tc[1]
            +6.075883333333333e-07 * tc[2]
            -4.510760000000000e-11 * tc[3]
            +1.448422000000000e-15 * tc[4];
        /*species 33: c2h5oh */
        species[33] =
            -3.152562100000000e+04 / tc[1]
            +1.603545670000000e+01
            -6.562436500000000e+00 * tc[0]
            -7.602111000000000e-03 * tc[1]
            +8.982799166666666e-07 * tc[2]
            -7.185417583333333e-11 * tc[3]
            +2.564489350000000e-15 * tc[4];
        /*species 37: n2 */
        species[37] =
            -9.227977000000000e+02 / tc[1]
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439885000000000e-04 * tc[1]
            +9.474601666666666e-08 * tc[2]
            -8.414199999999999e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
    }

    /*species with midpoint at T=1451 kelvin */
    if (T < 1451) {
        /*species 13: ho2 */
        species[13] =
            +6.170734460000000e+02 / tc[1]
            -3.546790170000000e+00
            -3.476294990000000e+00 * tc[0]
            -1.102338725000000e-03 * tc[1]
            -2.614031266666667e-07 * tc[2]
            +1.772965016666667e-10 * tc[3]
            -2.915703590000000e-14 * tc[4];
    } else {
        /*species 13: ho2 */
        species[13] =
            +3.133658590000000e+01 / tc[1]
            +4.236915604000000e+00
            -4.583114990000000e+00 * tc[0]
            -8.636516950000000e-04 * tc[1]
            +1.031971893333333e-07 * tc[2]
            -8.266159941666667e-12 * tc[3]
            +2.941107835000000e-16 * tc[4];
    }

    /*species with midpoint at T=1359 kelvin */
    if (T < 1359) {
        /*species 28: ch2s */
        species[28] =
            +5.049544690000000e+04 / tc[1]
            +2.293920400000000e-01
            -3.323407010000000e+00 * tc[0]
            -1.140312895000000e-03 * tc[1]
            +4.135933600000000e-08 * tc[2]
            +8.505448833333333e-12 * tc[3]
            -1.244133095000000e-15 * tc[4];
    } else {
        /*species 28: ch2s */
        species[28] =
            +5.057396970000000e+04 / tc[1]
            -1.217650890000000e+00
            -3.096271800000000e+00 * tc[0]
            -1.404195925000000e-03 * tc[1]
            +1.189357543333333e-07 * tc[2]
            -7.017802483333333e-12 * tc[3]
            +1.923204140000000e-16 * tc[4];
    }

    /*species with midpoint at T=1524 kelvin */
    if (T < 1524) {
        /*species 35: ch3choh */
        species[35] =
            -6.128730360000000e+03 / tc[1]
            -1.199066214000000e+01
            -2.100268060000000e+00 * tc[0]
            -9.631333949999999e-03 * tc[1]
            +8.623957266666667e-07 * tc[2]
            +1.634921516666667e-10 * tc[3]
            -4.681985300000000e-14 * tc[4];
    } else {
        /*species 35: ch3choh */
        species[35] =
            -8.621997680000000e+03 / tc[1]
            +2.699544037000000e+01
            -7.954327470000000e+00 * tc[0]
            -5.250408500000000e-03 * tc[1]
            +5.843859833333333e-07 * tc[2]
            -4.468121733333333e-11 * tc[3]
            +1.539892295000000e-15 * tc[4];
    }

    /*species with midpoint at T=1405 kelvin */
    if (T < 1405) {
        /*species 34: c2h4oh */
        species[34] =
            -5.468083960000000e+03 / tc[1]
            -2.132986089400000e+01
            -7.761576060000001e-01 * tc[0]
            -1.455012570000000e-02 * tc[1]
            +3.566342816666667e-06 * tc[2]
            -6.806949116666667e-10 * tc[3]
            +6.279580000000001e-14 * tc[4];
        /*species 36: ch3ch2o */
        species[36] =
            -3.163971960000000e+03 / tc[1]
            -2.504189667800000e+01
            +2.712963780000000e-01 * tc[0]
            -1.494199060000000e-02 * tc[1]
            +3.284842466666667e-06 * tc[2]
            -5.311165775000000e-10 * tc[3]
            +3.889825270000000e-14 * tc[4];
    } else {
        /*species 34: c2h4oh */
        species[34] =
            -8.096020470000000e+03 / tc[1]
            +2.908263786000000e+01
            -8.755444960000000e+00 * tc[0]
            -4.747265735000000e-03 * tc[1]
            +5.146726250000000e-07 * tc[2]
            -3.862813508333334e-11 * tc[3]
            +1.313650220000000e-15 * tc[4];
        /*species 36: ch3ch2o */
        species[36] =
            -6.130979540000000e+03 / tc[1]
            +2.971038202000000e+01
            -8.311823920000000e+00 * tc[0]
            -5.171315950000000e-03 * tc[1]
            +5.653101483333333e-07 * tc[2]
            -4.268438475000000e-11 * tc[3]
            +1.458008565000000e-15 * tc[4];
    }

    /*species with midpoint at T=1375 kelvin */
    if (T < 1375) {
        /*species 21: c2h5 */
        species[21] =
            +1.333269920000000e+04 / tc[1]
            -1.516123668000000e+01
            -1.473748520000000e+00 * tc[0]
            -8.180330550000000e-03 * tc[1]
            +7.214947100000000e-07 * tc[2]
            +9.772555416666666e-11 * tc[3]
            -2.861629415000000e-14 * tc[4];
    } else {
        /*species 21: c2h5 */
        species[21] =
            +1.145398450000000e+04 / tc[1]
            +1.262368499000000e+01
            -5.601160910000000e+00 * tc[0]
            -5.348854050000000e-03 * tc[1]
            +6.058413000000000e-07 * tc[2]
            -4.681805450000000e-11 * tc[3]
            +1.624573570000000e-15 * tc[4];
    }
    return;
}


/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void helmholtz(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1410 kelvin */
    if (T < 1410) {
        /*species 23: ch2oh */
        species[23] =
            -2.33478724e+03 / tc[1]
            -9.72655221e+00
            -2.60067849e+00 * tc[0]
            -6.42239270e-03 * tc[1]
            +1.38966031e-06 * tc[2]
            -2.29773005e-10 * tc[3]
            +1.78520553e-14 * tc[4];
    } else {
        /*species 23: ch2oh */
        species[23] =
            -3.50157098e+03 / tc[1]
            +1.19296465e+01
            -6.00127803e+00 * tc[0]
            -2.49360784e-03 * tc[1]
            +2.68255858e-07 * tc[2]
            -2.00189279e-11 * tc[3]
            +6.77913500e-16 * tc[4];
    }

    /*species with midpoint at T=1350 kelvin */
    if (T < 1350) {
        /*species 4: ch2 */
        species[4] =
            +4.59192771e+04 / tc[1]
            -4.53403650e-01
            -3.49623989e+00 * tc[0]
            -1.10129322e-03 * tc[1]
            +7.06381332e-08 * tc[2]
            -1.13737255e-12 * tc[3]
            -2.15968686e-16 * tc[4];
    } else {
        /*species 4: ch2 */
        species[4] =
            +4.59375743e+04 / tc[1]
            -7.91742480e-01
            -3.44310466e+00 * tc[0]
            -1.16319005e-03 * tc[1]
            +8.88579997e-08 * tc[2]
            -4.78369973e-12 * tc[3]
            +1.21347424e-16 * tc[4];
    }

    /*species with midpoint at T=1671 kelvin */
    if (T < 1671) {
        /*species 19: c2h3 */
        species[19] =
            +3.42868979e+04 / tc[1]
            -8.41389408e+00
            -2.73925942e+00 * tc[0]
            -3.51505795e-03 * tc[1]
            -3.94412165e-07 * tc[2]
            +2.99641447e-10 * tc[3]
            -4.45878245e-14 * tc[4];
    } else {
        /*species 19: c2h3 */
        species[19] =
            +3.35153544e+04 / tc[1]
            +7.03842990e-01
            -3.96047713e+00 * tc[0]
            -3.99713006e-03 * tc[1]
            +4.76013478e-07 * tc[2]
            -3.81959011e-11 * tc[3]
            +1.36285070e-15 * tc[4];
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: h2 */
        species[0] =
            -1.01252100e+03 / tc[1]
            +5.59221800e+00
            -3.29812400e+00 * tc[0]
            -4.12472100e-04 * tc[1]
            +1.35716917e-07 * tc[2]
            +7.89619500e-12 * tc[3]
            -2.06743600e-14 * tc[4];
        /*species 1: h */
        species[1] =
            +2.54716300e+04 / tc[1]
            +1.96011760e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 2: ch4 */
        species[2] =
            -9.82522900e+03 / tc[1]
            -1.39434485e+01
            -7.78741500e-01 * tc[0]
            -8.73834000e-03 * tc[1]
            +4.63901500e-06 * tc[2]
            -2.54142333e-09 * tc[3]
            +6.11965500e-13 * tc[4];
        /*species 3: ch3 */
        species[3] =
            +1.64237800e+04 / tc[1]
            -5.35935100e+00
            -2.43044300e+00 * tc[0]
            -5.56205000e-03 * tc[1]
            +2.80036667e-06 * tc[2]
            -1.35152417e-09 * tc[3]
            +2.93247650e-13 * tc[4];
        /*species 5: ch */
        species[5] =
            +7.04525900e+04 / tc[1]
            -1.13138600e+00
            -3.20020200e+00 * tc[0]
            -1.03643800e-03 * tc[1]
            +8.55738500e-07 * tc[2]
            -4.77824167e-10 * tc[3]
            +9.77766500e-14 * tc[4];
        /*species 6: ch2o */
        species[6] =
            -1.48654000e+04 / tc[1]
            -1.31320890e+01
            -1.65273100e+00 * tc[0]
            -6.31572000e-03 * tc[1]
            +3.14694667e-06 * tc[2]
            -1.70835917e-09 * tc[3]
            +4.20661850e-13 * tc[4];
        /*species 7: hco */
        species[7] =
            +4.15992200e+03 / tc[1]
            -7.08528400e+00
            -2.89833000e+00 * tc[0]
            -3.09957350e-03 * tc[1]
            +1.60384733e-06 * tc[2]
            -9.08187500e-10 * tc[3]
            +2.28744250e-13 * tc[4];
        /*species 8: co2 */
        species[8] =
            -4.83731400e+04 / tc[1]
            -8.91276500e+00
            -2.27572500e+00 * tc[0]
            -4.96103600e-03 * tc[1]
            +1.73485167e-06 * tc[2]
            -5.72223917e-10 * tc[3]
            +1.05864000e-13 * tc[4];
        /*species 9: co */
        species[9] =
            -1.43105400e+04 / tc[1]
            -2.58644500e+00
            -3.26245200e+00 * tc[0]
            -7.55970500e-04 * tc[1]
            +6.46959167e-07 * tc[2]
            -4.65162000e-10 * tc[3]
            +1.23747550e-13 * tc[4];
        /*species 10: o2 */
        species[10] =
            -1.00524900e+03 / tc[1]
            -3.82180200e+00
            -3.21293600e+00 * tc[0]
            -5.63743000e-04 * tc[1]
            +9.59358333e-08 * tc[2]
            -1.09489750e-10 * tc[3]
            +4.38427700e-14 * tc[4];
        /*species 11: o */
        species[11] =
            +2.91476400e+04 / tc[1]
            -1.01756600e+00
            -2.94642900e+00 * tc[0]
            +8.19083000e-04 * tc[1]
            -4.03505333e-07 * tc[2]
            +1.33570250e-10 * tc[3]
            -1.94534800e-14 * tc[4];
        /*species 12: oh */
        species[12] =
            +3.60678200e+03 / tc[1]
            +1.27840600e+00
            -3.63726600e+00 * tc[0]
            -9.25455000e-05 * tc[1]
            +2.79360833e-07 * tc[2]
            -1.98933583e-10 * tc[3]
            +4.21572100e-14 * tc[4];
        /*species 14: h2o2 */
        species[14] =
            -1.76631500e+04 / tc[1]
            -4.39660900e+00
            -3.38875400e+00 * tc[0]
            -3.28461300e-03 * tc[1]
            +2.47502167e-08 * tc[2]
            +3.85483833e-10 * tc[3]
            -1.23575750e-13 * tc[4];
        /*species 15: h2o */
        species[15] =
            -3.02081100e+04 / tc[1]
            -2.03391000e-01
            -3.38684200e+00 * tc[0]
            -1.73749100e-03 * tc[1]
            +1.05911600e-06 * tc[2]
            -5.80715083e-10 * tc[3]
            +1.25329400e-13 * tc[4];
        /*species 16: c2h */
        species[16] =
            +6.68381300e+04 / tc[1]
            -5.56251600e+00
            -2.73770400e+00 * tc[0]
            -4.02422300e-03 * tc[1]
            +1.54071833e-06 * tc[2]
            -5.43771583e-10 * tc[3]
            +9.69790000e-14 * tc[4];
        /*species 17: hcco */
        species[17] =
            +1.96589200e+04 / tc[1]
            +3.56612110e+00
            -5.04796500e+00 * tc[0]
            -2.22673900e-03 * tc[1]
            -3.78047167e-08 * tc[2]
            +1.23507917e-10 * tc[3]
            -1.12537100e-14 * tc[4];
        /*species 18: c2h2 */
        species[18] =
            +2.61244400e+04 / tc[1]
            -7.79181600e+00
            -2.01356200e+00 * tc[0]
            -7.59522500e-03 * tc[1]
            +2.69386500e-06 * tc[2]
            -7.56582667e-10 * tc[3]
            +9.56373000e-14 * tc[4];
        /*species 20: c2h4 */
        species[20] =
            +5.57304600e+03 / tc[1]
            -2.60729780e+01
            +8.61488000e-01 * tc[0]
            -1.39808150e-02 * tc[1]
            +5.64779500e-06 * tc[2]
            -2.32096000e-09 * tc[3]
            +4.86893950e-13 * tc[4];
        /*species 22: c2h6 */
        species[22] =
            -1.12391800e+04 / tc[1]
            -1.39697510e+01
            -1.46253900e+00 * tc[0]
            -7.74733500e-03 * tc[1]
            -9.63417833e-07 * tc[2]
            +1.04819333e-09 * tc[3]
            -2.29313350e-13 * tc[4];
        /*species 24: ch3o */
        species[24] =
            +9.78601100e+02 / tc[1]
            -1.20459760e+01
            -2.10620400e+00 * tc[0]
            -3.60829750e-03 * tc[1]
            -8.89745333e-07 * tc[2]
            +6.14803000e-10 * tc[3]
            -1.03780550e-13 * tc[4];
        /*species 25: hccoh */
        species[25] =
            +8.70119000e+03 / tc[1]
            -1.59241000e+00
            -3.89946500e+00 * tc[0]
            -4.85053750e-03 * tc[1]
            +5.19884833e-08 * tc[2]
            +4.61477667e-10 * tc[3]
            -1.23286600e-13 * tc[4];
        /*species 26: h2ccch */
        species[26] =
            +3.98888300e+04 / tc[1]
            +3.16874510e+00
            -4.75420000e+00 * tc[0]
            -5.54014000e-03 * tc[1]
            -4.65553833e-08 * tc[2]
            +4.56601000e-10 * tc[3]
            -9.74814500e-14 * tc[4];
        /*species 27: c3h2 */
        species[27] =
            +6.35042100e+04 / tc[1]
            -6.70273200e+00
            -3.16671400e+00 * tc[0]
            -1.24128600e-02 * tc[1]
            +7.65272833e-06 * tc[2]
            -3.55668250e-09 * tc[3]
            +7.41076000e-13 * tc[4];
        /*species 29: ch2co */
        species[29] =
            -7.63263700e+03 / tc[1]
            -6.69858200e+00
            -2.97497100e+00 * tc[0]
            -6.05935500e-03 * tc[1]
            +3.90841000e-07 * tc[2]
            +5.38890417e-10 * tc[3]
            -1.95282450e-13 * tc[4];
        /*species 30: ch2hco */
        species[30] =
            +1.52147700e+03 / tc[1]
            -7.14922800e+00
            -3.40906200e+00 * tc[0]
            -5.36928500e-03 * tc[1]
            -3.15248667e-07 * tc[2]
            +5.96548583e-10 * tc[3]
            -1.43369250e-13 * tc[4];
        /*species 31: ch3co */
        species[31] =
            -4.10850800e+03 / tc[1]
            -9.10357200e+00
            -3.12527800e+00 * tc[0]
            -4.88911000e-03 * tc[1]
            -7.53574667e-07 * tc[2]
            +7.50788500e-10 * tc[3]
            -1.59685900e-13 * tc[4];
        /*species 32: ch3hco */
        species[32] =
            -2.12458900e+04 / tc[1]
            -1.18451950e+01
            -2.50569500e+00 * tc[0]
            -6.68495500e-03 * tc[1]
            -7.78658833e-07 * tc[2]
            +9.40116667e-10 * tc[3]
            -2.13178300e-13 * tc[4];
        /*species 33: c2h5oh */
        species[33] =
            -2.99961320e+04 / tc[1]
            -9.43158800e-01
            -4.85869570e+00 * tc[0]
            +1.87008630e-03 * tc[1]
            -1.15925630e-05 * tc[2]
            +7.38789967e-09 * tc[3]
            -1.75844175e-12 * tc[4];
        /*species 37: n2 */
        species[37] =
            -1.02090000e+03 / tc[1]
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120000e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242750e-13 * tc[4];
    } else {
        /*species 0: h2 */
        species[0] =
            -8.35034000e+02 / tc[1]
            +3.34653300e+00
            -2.99142300e+00 * tc[0]
            -3.50032200e-04 * tc[1]
            +9.38971500e-09 * tc[2]
            +7.69298167e-13 * tc[3]
            -7.91376000e-17 * tc[4];
        /*species 1: h */
        species[1] =
            +2.54716300e+04 / tc[1]
            +1.96011760e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 2: ch4 */
        species[2] =
            -1.00807900e+04 / tc[1]
            -8.93991600e+00
            -1.68347900e+00 * tc[0]
            -5.11862000e-03 * tc[1]
            +6.45854833e-07 * tc[2]
            -5.65465417e-11 * tc[3]
            +2.25171150e-15 * tc[4];
        /*species 3: ch3 */
        species[3] =
            +1.64378100e+04 / tc[1]
            -3.60864500e+00
            -2.84405200e+00 * tc[0]
            -3.06898700e-03 * tc[1]
            +3.71724167e-07 * tc[2]
            -3.15430083e-11 * tc[3]
            +1.22607950e-15 * tc[4];
        /*species 5: ch */
        species[5] =
            +7.08672300e+04 / tc[1]
            -7.98215000e+00
            -2.19622300e+00 * tc[0]
            -1.17019050e-03 * tc[1]
            +1.17636683e-07 * tc[2]
            -7.50631833e-12 * tc[3]
            +1.92752000e-16 * tc[4];
        /*species 6: ch2o */
        species[6] =
            -1.53203700e+04 / tc[1]
            -4.91696600e+00
            -2.99560600e+00 * tc[0]
            -3.34066050e-03 * tc[1]
            +4.38159167e-07 * tc[2]
            -3.94762750e-11 * tc[3]
            +1.60625850e-15 * tc[4];
        /*species 7: hco */
        species[7] =
            +3.91632400e+03 / tc[1]
            -2.99502800e+00
            -3.55727100e+00 * tc[0]
            -1.67278650e-03 * tc[1]
            +2.22501000e-07 * tc[2]
            -2.05881083e-11 * tc[3]
            +8.56925500e-16 * tc[4];
        /*species 8: co2 */
        species[8] =
            -4.89669600e+04 / tc[1]
            +4.40901890e+00
            -4.45362300e+00 * tc[0]
            -1.57008450e-03 * tc[1]
            +2.13068500e-07 * tc[2]
            -1.99499750e-11 * tc[3]
            +8.34516500e-16 * tc[4];
        /*species 9: co */
        species[9] =
            -1.42683500e+04 / tc[1]
            -4.08314000e+00
            -3.02507800e+00 * tc[0]
            -7.21344500e-04 * tc[1]
            +9.38471333e-08 * tc[2]
            -8.48817500e-12 * tc[3]
            +3.45547600e-16 * tc[4];
        /*species 10: o2 */
        species[10] =
            -1.23393000e+03 / tc[1]
            -4.91588000e-01
            -3.69757800e+00 * tc[0]
            -3.06759850e-04 * tc[1]
            +2.09807000e-08 * tc[2]
            -1.47940083e-12 * tc[3]
            +5.68217500e-17 * tc[4];
        /*species 11: o */
        species[11] =
            +2.92308000e+04 / tc[1]
            -3.37824800e+00
            -2.54206000e+00 * tc[0]
            +1.37753100e-05 * tc[1]
            +5.17133833e-10 * tc[2]
            -3.79255583e-13 * tc[3]
            +2.18402600e-17 * tc[4];
        /*species 12: oh */
        species[12] =
            +3.88688800e+03 / tc[1]
            -3.71298200e+00
            -2.88273000e+00 * tc[0]
            -5.06987000e-04 * tc[1]
            +3.79479500e-08 * tc[2]
            -1.81223667e-12 * tc[3]
            +2.56315250e-17 * tc[4];
        /*species 14: h2o2 */
        species[14] =
            -1.80069600e+04 / tc[1]
            +3.07203000e+00
            -4.57316700e+00 * tc[0]
            -2.16806800e-03 * tc[1]
            +2.45781500e-07 * tc[2]
            -1.95742000e-11 * tc[3]
            +7.15827000e-16 * tc[4];
        /*species 15: h2o */
        species[15] =
            -2.98992100e+04 / tc[1]
            -5.19067100e+00
            -2.67214600e+00 * tc[0]
            -1.52814650e-03 * tc[1]
            +1.45504333e-07 * tc[2]
            -1.00083000e-11 * tc[3]
            +3.19580900e-16 * tc[4];
        /*species 16: c2h */
        species[16] =
            +6.65588400e+04 / tc[1]
            +1.79530400e+00
            -3.98636700e+00 * tc[0]
            -1.57156150e-03 * tc[1]
            +2.11207167e-07 * tc[2]
            -2.43696917e-11 * tc[3]
            +1.35816000e-15 * tc[4];
        /*species 17: hcco */
        species[17] =
            +1.90151300e+04 / tc[1]
            +1.48293350e+01
            -6.75807300e+00 * tc[0]
            -1.00020000e-03 * tc[1]
            +3.37934500e-08 * tc[2]
            +8.67610000e-12 * tc[3]
            -9.82582500e-16 * tc[4];
        /*species 18: c2h2 */
        species[18] =
            +2.56676600e+04 / tc[1]
            +6.23710800e+00
            -4.43677000e+00 * tc[0]
            -2.68801950e-03 * tc[1]
            +3.18802833e-07 * tc[2]
            -2.73864917e-11 * tc[3]
            +1.07835500e-15 * tc[4];
        /*species 20: c2h4 */
        species[20] =
            +4.42828900e+03 / tc[1]
            +2.98030000e-01
            -3.52841900e+00 * tc[0]
            -5.74259000e-03 * tc[1]
            +7.36397500e-07 * tc[2]
            -6.53716750e-11 * tc[3]
            +2.63342400e-15 * tc[4];
        /*species 22: c2h6 */
        species[22] =
            -1.27177900e+04 / tc[1]
            +9.06544500e+00
            -4.82593800e+00 * tc[0]
            -6.92021500e-03 * tc[1]
            +7.59543167e-07 * tc[2]
            -5.60413917e-11 * tc[3]
            +1.79908050e-15 * tc[4];
        /*species 24: ch3o */
        species[24] =
            +1.27832500e+02 / tc[1]
            -1.58775000e-01
            -3.77080000e+00 * tc[0]
            -3.93574850e-03 * tc[1]
            +4.42730667e-07 * tc[2]
            -3.28702583e-11 * tc[3]
            +1.05630800e-15 * tc[4];
        /*species 25: hccoh */
        species[25] =
            +7.59825800e+03 / tc[1]
            +2.03404640e+01
            -7.32832400e+00 * tc[0]
            -1.66820800e-03 * tc[1]
            +5.04117500e-08 * tc[2]
            +1.48425500e-11 * tc[3]
            -1.62258400e-15 * tc[4];
        /*species 26: h2ccch */
        species[26] =
            +3.84742000e+04 / tc[1]
            +2.96102370e+01
            -8.83104700e+00 * tc[0]
            -2.17859750e-03 * tc[1]
            +6.84844500e-08 * tc[2]
            +1.97393583e-11 * tc[3]
            -2.18826000e-15 * tc[4];
        /*species 27: c3h2 */
        species[27] =
            +6.25972200e+04 / tc[1]
            +1.90398810e+01
            -7.67098100e+00 * tc[0]
            -1.37437450e-03 * tc[1]
            +7.28490500e-08 * tc[2]
            +5.37966583e-12 * tc[3]
            -8.31943500e-16 * tc[4];
        /*species 29: ch2co */
        species[29] =
            -8.58340200e+03 / tc[1]
            +1.26963980e+01
            -6.03881700e+00 * tc[0]
            -2.90242000e-03 * tc[1]
            +3.20159000e-07 * tc[2]
            -2.32873750e-11 * tc[3]
            +7.29434000e-16 * tc[4];
        /*species 30: ch2hco */
        species[30] =
            +4.90321800e+02 / tc[1]
            +1.00209210e+01
            -5.97567000e+00 * tc[0]
            -4.06529550e-03 * tc[1]
            +4.57270667e-07 * tc[2]
            -3.39192000e-11 * tc[3]
            +1.08800850e-15 * tc[4];
        /*species 31: ch3co */
        species[31] =
            -5.18786300e+03 / tc[1]
            +7.88722800e+00
            -5.61227900e+00 * tc[0]
            -4.22494300e-03 * tc[1]
            +4.75691167e-07 * tc[2]
            -3.53198000e-11 * tc[3]
            +1.13420200e-15 * tc[4];
        /*species 32: ch3hco */
        species[32] =
            -2.26456900e+04 / tc[1]
            +1.08815960e+01
            -5.86865000e+00 * tc[0]
            -5.39712000e-03 * tc[1]
            +6.07588333e-07 * tc[2]
            -4.51076000e-11 * tc[3]
            +1.44842200e-15 * tc[4];
        /*species 33: c2h5oh */
        species[33] =
            -3.15256210e+04 / tc[1]
            +1.50354567e+01
            -6.56243650e+00 * tc[0]
            -7.60211100e-03 * tc[1]
            +8.98279917e-07 * tc[2]
            -7.18541758e-11 * tc[3]
            +2.56448935e-15 * tc[4];
        /*species 37: n2 */
        species[37] =
            -9.22797700e+02 / tc[1]
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988500e-04 * tc[1]
            +9.47460167e-08 * tc[2]
            -8.41420000e-12 * tc[3]
            +3.37667550e-16 * tc[4];
    }

    /*species with midpoint at T=1451 kelvin */
    if (T < 1451) {
        /*species 13: ho2 */
        species[13] =
            +6.17073446e+02 / tc[1]
            -4.54679017e+00
            -3.47629499e+00 * tc[0]
            -1.10233872e-03 * tc[1]
            -2.61403127e-07 * tc[2]
            +1.77296502e-10 * tc[3]
            -2.91570359e-14 * tc[4];
    } else {
        /*species 13: ho2 */
        species[13] =
            +3.13365859e+01 / tc[1]
            +3.23691560e+00
            -4.58311499e+00 * tc[0]
            -8.63651695e-04 * tc[1]
            +1.03197189e-07 * tc[2]
            -8.26615994e-12 * tc[3]
            +2.94110784e-16 * tc[4];
    }

    /*species with midpoint at T=1359 kelvin */
    if (T < 1359) {
        /*species 28: ch2s */
        species[28] =
            +5.04954469e+04 / tc[1]
            -7.70607960e-01
            -3.32340701e+00 * tc[0]
            -1.14031289e-03 * tc[1]
            +4.13593360e-08 * tc[2]
            +8.50544883e-12 * tc[3]
            -1.24413309e-15 * tc[4];
    } else {
        /*species 28: ch2s */
        species[28] =
            +5.05739697e+04 / tc[1]
            -2.21765089e+00
            -3.09627180e+00 * tc[0]
            -1.40419592e-03 * tc[1]
            +1.18935754e-07 * tc[2]
            -7.01780248e-12 * tc[3]
            +1.92320414e-16 * tc[4];
    }

    /*species with midpoint at T=1524 kelvin */
    if (T < 1524) {
        /*species 35: ch3choh */
        species[35] =
            -6.12873036e+03 / tc[1]
            -1.29906621e+01
            -2.10026806e+00 * tc[0]
            -9.63133395e-03 * tc[1]
            +8.62395727e-07 * tc[2]
            +1.63492152e-10 * tc[3]
            -4.68198530e-14 * tc[4];
    } else {
        /*species 35: ch3choh */
        species[35] =
            -8.62199768e+03 / tc[1]
            +2.59954404e+01
            -7.95432747e+00 * tc[0]
            -5.25040850e-03 * tc[1]
            +5.84385983e-07 * tc[2]
            -4.46812173e-11 * tc[3]
            +1.53989229e-15 * tc[4];
    }

    /*species with midpoint at T=1405 kelvin */
    if (T < 1405) {
        /*species 34: c2h4oh */
        species[34] =
            -5.46808396e+03 / tc[1]
            -2.23298609e+01
            -7.76157606e-01 * tc[0]
            -1.45501257e-02 * tc[1]
            +3.56634282e-06 * tc[2]
            -6.80694912e-10 * tc[3]
            +6.27958000e-14 * tc[4];
        /*species 36: ch3ch2o */
        species[36] =
            -3.16397196e+03 / tc[1]
            -2.60418967e+01
            +2.71296378e-01 * tc[0]
            -1.49419906e-02 * tc[1]
            +3.28484247e-06 * tc[2]
            -5.31116578e-10 * tc[3]
            +3.88982527e-14 * tc[4];
    } else {
        /*species 34: c2h4oh */
        species[34] =
            -8.09602047e+03 / tc[1]
            +2.80826379e+01
            -8.75544496e+00 * tc[0]
            -4.74726573e-03 * tc[1]
            +5.14672625e-07 * tc[2]
            -3.86281351e-11 * tc[3]
            +1.31365022e-15 * tc[4];
        /*species 36: ch3ch2o */
        species[36] =
            -6.13097954e+03 / tc[1]
            +2.87103820e+01
            -8.31182392e+00 * tc[0]
            -5.17131595e-03 * tc[1]
            +5.65310148e-07 * tc[2]
            -4.26843848e-11 * tc[3]
            +1.45800857e-15 * tc[4];
    }

    /*species with midpoint at T=1375 kelvin */
    if (T < 1375) {
        /*species 21: c2h5 */
        species[21] =
            +1.33326992e+04 / tc[1]
            -1.61612367e+01
            -1.47374852e+00 * tc[0]
            -8.18033055e-03 * tc[1]
            +7.21494710e-07 * tc[2]
            +9.77255542e-11 * tc[3]
            -2.86162942e-14 * tc[4];
    } else {
        /*species 21: c2h5 */
        species[21] =
            +1.14539845e+04 / tc[1]
            +1.16236850e+01
            -5.60116091e+00 * tc[0]
            -5.34885405e-03 * tc[1]
            +6.05841300e-07 * tc[2]
            -4.68180545e-11 * tc[3]
            +1.62457357e-15 * tc[4];
    }
    return;
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cv_R(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1410 kelvin */
    if (T < 1410) {
        /*species 23: ch2oh */
        species[23] =
            +1.60067849e+00
            +1.28447854e-02 * tc[1]
            -8.33796185e-06 * tc[2]
            +2.75727606e-09 * tc[3]
            -3.57041106e-13 * tc[4];
    } else {
        /*species 23: ch2oh */
        species[23] =
            +5.00127803e+00
            +4.98721568e-03 * tc[1]
            -1.60953515e-06 * tc[2]
            +2.40227135e-10 * tc[3]
            -1.35582700e-14 * tc[4];
    }

    /*species with midpoint at T=1350 kelvin */
    if (T < 1350) {
        /*species 4: ch2 */
        species[4] =
            +2.49623989e+00
            +2.20258645e-03 * tc[1]
            -4.23828799e-07 * tc[2]
            +1.36484706e-11 * tc[3]
            +4.31937372e-15 * tc[4];
    } else {
        /*species 4: ch2 */
        species[4] =
            +2.44310466e+00
            +2.32638011e-03 * tc[1]
            -5.33147998e-07 * tc[2]
            +5.74043968e-11 * tc[3]
            -2.42694849e-15 * tc[4];
    }

    /*species with midpoint at T=1671 kelvin */
    if (T < 1671) {
        /*species 19: c2h3 */
        species[19] =
            +1.73925942e+00
            +7.03011591e-03 * tc[1]
            +2.36647299e-06 * tc[2]
            -3.59569737e-09 * tc[3]
            +8.91756491e-13 * tc[4];
    } else {
        /*species 19: c2h3 */
        species[19] =
            +2.96047713e+00
            +7.99426013e-03 * tc[1]
            -2.85608087e-06 * tc[2]
            +4.58350813e-10 * tc[3]
            -2.72570140e-14 * tc[4];
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: h2 */
        species[0] =
            +2.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 1: h */
        species[1] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 2: ch4 */
        species[2] =
            -2.21258500e-01
            +1.74766800e-02 * tc[1]
            -2.78340900e-05 * tc[2]
            +3.04970800e-08 * tc[3]
            -1.22393100e-11 * tc[4];
        /*species 3: ch3 */
        species[3] =
            +1.43044300e+00
            +1.11241000e-02 * tc[1]
            -1.68022000e-05 * tc[2]
            +1.62182900e-08 * tc[3]
            -5.86495300e-12 * tc[4];
        /*species 5: ch */
        species[5] =
            +2.20020200e+00
            +2.07287600e-03 * tc[1]
            -5.13443100e-06 * tc[2]
            +5.73389000e-09 * tc[3]
            -1.95553300e-12 * tc[4];
        /*species 6: ch2o */
        species[6] =
            +6.52731000e-01
            +1.26314400e-02 * tc[1]
            -1.88816800e-05 * tc[2]
            +2.05003100e-08 * tc[3]
            -8.41323700e-12 * tc[4];
        /*species 7: hco */
        species[7] =
            +1.89833000e+00
            +6.19914700e-03 * tc[1]
            -9.62308400e-06 * tc[2]
            +1.08982500e-08 * tc[3]
            -4.57488500e-12 * tc[4];
        /*species 8: co2 */
        species[8] =
            +1.27572500e+00
            +9.92207200e-03 * tc[1]
            -1.04091100e-05 * tc[2]
            +6.86668700e-09 * tc[3]
            -2.11728000e-12 * tc[4];
        /*species 9: co */
        species[9] =
            +2.26245200e+00
            +1.51194100e-03 * tc[1]
            -3.88175500e-06 * tc[2]
            +5.58194400e-09 * tc[3]
            -2.47495100e-12 * tc[4];
        /*species 10: o2 */
        species[10] =
            +2.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
        /*species 11: o */
        species[11] =
            +1.94642900e+00
            -1.63816600e-03 * tc[1]
            +2.42103200e-06 * tc[2]
            -1.60284300e-09 * tc[3]
            +3.89069600e-13 * tc[4];
        /*species 12: oh */
        species[12] =
            +2.63726600e+00
            +1.85091000e-04 * tc[1]
            -1.67616500e-06 * tc[2]
            +2.38720300e-09 * tc[3]
            -8.43144200e-13 * tc[4];
        /*species 14: h2o2 */
        species[14] =
            +2.38875400e+00
            +6.56922600e-03 * tc[1]
            -1.48501300e-07 * tc[2]
            -4.62580600e-09 * tc[3]
            +2.47151500e-12 * tc[4];
        /*species 15: h2o */
        species[15] =
            +2.38684200e+00
            +3.47498200e-03 * tc[1]
            -6.35469600e-06 * tc[2]
            +6.96858100e-09 * tc[3]
            -2.50658800e-12 * tc[4];
        /*species 16: c2h */
        species[16] =
            +1.73770400e+00
            +8.04844600e-03 * tc[1]
            -9.24431000e-06 * tc[2]
            +6.52525900e-09 * tc[3]
            -1.93958000e-12 * tc[4];
        /*species 17: hcco */
        species[17] =
            +4.04796500e+00
            +4.45347800e-03 * tc[1]
            +2.26828300e-07 * tc[2]
            -1.48209500e-09 * tc[3]
            +2.25074200e-13 * tc[4];
        /*species 18: c2h2 */
        species[18] =
            +1.01356200e+00
            +1.51904500e-02 * tc[1]
            -1.61631900e-05 * tc[2]
            +9.07899200e-09 * tc[3]
            -1.91274600e-12 * tc[4];
        /*species 20: c2h4 */
        species[20] =
            -1.86148800e+00
            +2.79616300e-02 * tc[1]
            -3.38867700e-05 * tc[2]
            +2.78515200e-08 * tc[3]
            -9.73787900e-12 * tc[4];
        /*species 22: c2h6 */
        species[22] =
            +4.62539000e-01
            +1.54946700e-02 * tc[1]
            +5.78050700e-06 * tc[2]
            -1.25783200e-08 * tc[3]
            +4.58626700e-12 * tc[4];
        /*species 24: ch3o */
        species[24] =
            +1.10620400e+00
            +7.21659500e-03 * tc[1]
            +5.33847200e-06 * tc[2]
            -7.37763600e-09 * tc[3]
            +2.07561100e-12 * tc[4];
        /*species 25: hccoh */
        species[25] =
            +2.89946500e+00
            +9.70107500e-03 * tc[1]
            -3.11930900e-07 * tc[2]
            -5.53773200e-09 * tc[3]
            +2.46573200e-12 * tc[4];
        /*species 26: h2ccch */
        species[26] =
            +3.75420000e+00
            +1.10802800e-02 * tc[1]
            +2.79332300e-07 * tc[2]
            -5.47921200e-09 * tc[3]
            +1.94962900e-12 * tc[4];
        /*species 27: c3h2 */
        species[27] =
            +2.16671400e+00
            +2.48257200e-02 * tc[1]
            -4.59163700e-05 * tc[2]
            +4.26801900e-08 * tc[3]
            -1.48215200e-11 * tc[4];
        /*species 29: ch2co */
        species[29] =
            +1.97497100e+00
            +1.21187100e-02 * tc[1]
            -2.34504600e-06 * tc[2]
            -6.46668500e-09 * tc[3]
            +3.90564900e-12 * tc[4];
        /*species 30: ch2hco */
        species[30] =
            +2.40906200e+00
            +1.07385700e-02 * tc[1]
            +1.89149200e-06 * tc[2]
            -7.15858300e-09 * tc[3]
            +2.86738500e-12 * tc[4];
        /*species 31: ch3co */
        species[31] =
            +2.12527800e+00
            +9.77822000e-03 * tc[1]
            +4.52144800e-06 * tc[2]
            -9.00946200e-09 * tc[3]
            +3.19371800e-12 * tc[4];
        /*species 32: ch3hco */
        species[32] =
            +1.50569500e+00
            +1.33699100e-02 * tc[1]
            +4.67195300e-06 * tc[2]
            -1.12814000e-08 * tc[3]
            +4.26356600e-12 * tc[4];
        /*species 33: c2h5oh */
        species[33] =
            +3.85869570e+00
            -3.74017260e-03 * tc[1]
            +6.95553780e-05 * tc[2]
            -8.86547960e-08 * tc[3]
            +3.51688350e-11 * tc[4];
        /*species 37: n2 */
        species[37] =
            +2.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: h2 */
        species[0] =
            +1.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 1: h */
        species[1] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 2: ch4 */
        species[2] =
            +6.83479000e-01
            +1.02372400e-02 * tc[1]
            -3.87512900e-06 * tc[2]
            +6.78558500e-10 * tc[3]
            -4.50342300e-14 * tc[4];
        /*species 3: ch3 */
        species[3] =
            +1.84405200e+00
            +6.13797400e-03 * tc[1]
            -2.23034500e-06 * tc[2]
            +3.78516100e-10 * tc[3]
            -2.45215900e-14 * tc[4];
        /*species 5: ch */
        species[5] =
            +1.19622300e+00
            +2.34038100e-03 * tc[1]
            -7.05820100e-07 * tc[2]
            +9.00758200e-11 * tc[3]
            -3.85504000e-15 * tc[4];
        /*species 6: ch2o */
        species[6] =
            +1.99560600e+00
            +6.68132100e-03 * tc[1]
            -2.62895500e-06 * tc[2]
            +4.73715300e-10 * tc[3]
            -3.21251700e-14 * tc[4];
        /*species 7: hco */
        species[7] =
            +2.55727100e+00
            +3.34557300e-03 * tc[1]
            -1.33500600e-06 * tc[2]
            +2.47057300e-10 * tc[3]
            -1.71385100e-14 * tc[4];
        /*species 8: co2 */
        species[8] =
            +3.45362300e+00
            +3.14016900e-03 * tc[1]
            -1.27841100e-06 * tc[2]
            +2.39399700e-10 * tc[3]
            -1.66903300e-14 * tc[4];
        /*species 9: co */
        species[9] =
            +2.02507800e+00
            +1.44268900e-03 * tc[1]
            -5.63082800e-07 * tc[2]
            +1.01858100e-10 * tc[3]
            -6.91095200e-15 * tc[4];
        /*species 10: o2 */
        species[10] =
            +2.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
        /*species 11: o */
        species[11] =
            +1.54206000e+00
            -2.75506200e-05 * tc[1]
            -3.10280300e-09 * tc[2]
            +4.55106700e-12 * tc[3]
            -4.36805200e-16 * tc[4];
        /*species 12: oh */
        species[12] =
            +1.88273000e+00
            +1.01397400e-03 * tc[1]
            -2.27687700e-07 * tc[2]
            +2.17468400e-11 * tc[3]
            -5.12630500e-16 * tc[4];
        /*species 14: h2o2 */
        species[14] =
            +3.57316700e+00
            +4.33613600e-03 * tc[1]
            -1.47468900e-06 * tc[2]
            +2.34890400e-10 * tc[3]
            -1.43165400e-14 * tc[4];
        /*species 15: h2o */
        species[15] =
            +1.67214600e+00
            +3.05629300e-03 * tc[1]
            -8.73026000e-07 * tc[2]
            +1.20099600e-10 * tc[3]
            -6.39161800e-15 * tc[4];
        /*species 16: c2h */
        species[16] =
            +2.98636700e+00
            +3.14312300e-03 * tc[1]
            -1.26724300e-06 * tc[2]
            +2.92436300e-10 * tc[3]
            -2.71632000e-14 * tc[4];
        /*species 17: hcco */
        species[17] =
            +5.75807300e+00
            +2.00040000e-03 * tc[1]
            -2.02760700e-07 * tc[2]
            -1.04113200e-10 * tc[3]
            +1.96516500e-14 * tc[4];
        /*species 18: c2h2 */
        species[18] =
            +3.43677000e+00
            +5.37603900e-03 * tc[1]
            -1.91281700e-06 * tc[2]
            +3.28637900e-10 * tc[3]
            -2.15671000e-14 * tc[4];
        /*species 20: c2h4 */
        species[20] =
            +2.52841900e+00
            +1.14851800e-02 * tc[1]
            -4.41838500e-06 * tc[2]
            +7.84460100e-10 * tc[3]
            -5.26684800e-14 * tc[4];
        /*species 22: c2h6 */
        species[22] =
            +3.82593800e+00
            +1.38404300e-02 * tc[1]
            -4.55725900e-06 * tc[2]
            +6.72496700e-10 * tc[3]
            -3.59816100e-14 * tc[4];
        /*species 24: ch3o */
        species[24] =
            +2.77080000e+00
            +7.87149700e-03 * tc[1]
            -2.65638400e-06 * tc[2]
            +3.94443100e-10 * tc[3]
            -2.11261600e-14 * tc[4];
        /*species 25: hccoh */
        species[25] =
            +6.32832400e+00
            +3.33641600e-03 * tc[1]
            -3.02470500e-07 * tc[2]
            -1.78110600e-10 * tc[3]
            +3.24516800e-14 * tc[4];
        /*species 26: h2ccch */
        species[26] =
            +7.83104700e+00
            +4.35719500e-03 * tc[1]
            -4.10906700e-07 * tc[2]
            -2.36872300e-10 * tc[3]
            +4.37652000e-14 * tc[4];
        /*species 27: c3h2 */
        species[27] =
            +6.67098100e+00
            +2.74874900e-03 * tc[1]
            -4.37094300e-07 * tc[2]
            -6.45559900e-11 * tc[3]
            +1.66388700e-14 * tc[4];
        /*species 29: ch2co */
        species[29] =
            +5.03881700e+00
            +5.80484000e-03 * tc[1]
            -1.92095400e-06 * tc[2]
            +2.79448500e-10 * tc[3]
            -1.45886800e-14 * tc[4];
        /*species 30: ch2hco */
        species[30] =
            +4.97567000e+00
            +8.13059100e-03 * tc[1]
            -2.74362400e-06 * tc[2]
            +4.07030400e-10 * tc[3]
            -2.17601700e-14 * tc[4];
        /*species 31: ch3co */
        species[31] =
            +4.61227900e+00
            +8.44988600e-03 * tc[1]
            -2.85414700e-06 * tc[2]
            +4.23837600e-10 * tc[3]
            -2.26840400e-14 * tc[4];
        /*species 32: ch3hco */
        species[32] =
            +4.86865000e+00
            +1.07942400e-02 * tc[1]
            -3.64553000e-06 * tc[2]
            +5.41291200e-10 * tc[3]
            -2.89684400e-14 * tc[4];
        /*species 33: c2h5oh */
        species[33] =
            +5.56243650e+00
            +1.52042220e-02 * tc[1]
            -5.38967950e-06 * tc[2]
            +8.62250110e-10 * tc[3]
            -5.12897870e-14 * tc[4];
        /*species 37: n2 */
        species[37] =
            +1.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }

    /*species with midpoint at T=1451 kelvin */
    if (T < 1451) {
        /*species 13: ho2 */
        species[13] =
            +2.47629499e+00
            +2.20467745e-03 * tc[1]
            +1.56841876e-06 * tc[2]
            -2.12755802e-09 * tc[3]
            +5.83140718e-13 * tc[4];
    } else {
        /*species 13: ho2 */
        species[13] =
            +3.58311499e+00
            +1.72730339e-03 * tc[1]
            -6.19183136e-07 * tc[2]
            +9.91939193e-11 * tc[3]
            -5.88221567e-15 * tc[4];
    }

    /*species with midpoint at T=1359 kelvin */
    if (T < 1359) {
        /*species 28: ch2s */
        species[28] =
            +2.32340701e+00
            +2.28062579e-03 * tc[1]
            -2.48156016e-07 * tc[2]
            -1.02065386e-10 * tc[3]
            +2.48826619e-14 * tc[4];
    } else {
        /*species 28: ch2s */
        species[28] =
            +2.09627180e+00
            +2.80839185e-03 * tc[1]
            -7.13614526e-07 * tc[2]
            +8.42136298e-11 * tc[3]
            -3.84640828e-15 * tc[4];
    }

    /*species with midpoint at T=1524 kelvin */
    if (T < 1524) {
        /*species 35: ch3choh */
        species[35] =
            +1.10026806e+00
            +1.92626679e-02 * tc[1]
            -5.17437436e-06 * tc[2]
            -1.96190582e-09 * tc[3]
            +9.36397060e-13 * tc[4];
    } else {
        /*species 35: ch3choh */
        species[35] =
            +6.95432747e+00
            +1.05008170e-02 * tc[1]
            -3.50631590e-06 * tc[2]
            +5.36174608e-10 * tc[3]
            -3.07978459e-14 * tc[4];
    }

    /*species with midpoint at T=1405 kelvin */
    if (T < 1405) {
        /*species 34: c2h4oh */
        species[34] =
            -2.23842394e-01
            +2.91002514e-02 * tc[1]
            -2.13980569e-05 * tc[2]
            +8.16833894e-09 * tc[3]
            -1.25591600e-12 * tc[4];
        /*species 36: ch3ch2o */
        species[36] =
            -1.27129638e+00
            +2.98839812e-02 * tc[1]
            -1.97090548e-05 * tc[2]
            +6.37339893e-09 * tc[3]
            -7.77965054e-13 * tc[4];
    } else {
        /*species 34: c2h4oh */
        species[34] =
            +7.75544496e+00
            +9.49453147e-03 * tc[1]
            -3.08803575e-06 * tc[2]
            +4.63537621e-10 * tc[3]
            -2.62730044e-14 * tc[4];
        /*species 36: ch3ch2o */
        species[36] =
            +7.31182392e+00
            +1.03426319e-02 * tc[1]
            -3.39186089e-06 * tc[2]
            +5.12212617e-10 * tc[3]
            -2.91601713e-14 * tc[4];
    }

    /*species with midpoint at T=1375 kelvin */
    if (T < 1375) {
        /*species 21: c2h5 */
        species[21] =
            +4.73748520e-01
            +1.63606611e-02 * tc[1]
            -4.32896826e-06 * tc[2]
            -1.17270665e-09 * tc[3]
            +5.72325883e-13 * tc[4];
    } else {
        /*species 21: c2h5 */
        species[21] =
            +4.60116091e+00
            +1.06977081e-02 * tc[1]
            -3.63504780e-06 * tc[2]
            +5.61816654e-10 * tc[3]
            -3.24914714e-14 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cp_R(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1410 kelvin */
    if (T < 1410) {
        /*species 23: ch2oh */
        species[23] =
            +2.60067849e+00
            +1.28447854e-02 * tc[1]
            -8.33796185e-06 * tc[2]
            +2.75727606e-09 * tc[3]
            -3.57041106e-13 * tc[4];
    } else {
        /*species 23: ch2oh */
        species[23] =
            +6.00127803e+00
            +4.98721568e-03 * tc[1]
            -1.60953515e-06 * tc[2]
            +2.40227135e-10 * tc[3]
            -1.35582700e-14 * tc[4];
    }

    /*species with midpoint at T=1350 kelvin */
    if (T < 1350) {
        /*species 4: ch2 */
        species[4] =
            +3.49623989e+00
            +2.20258645e-03 * tc[1]
            -4.23828799e-07 * tc[2]
            +1.36484706e-11 * tc[3]
            +4.31937372e-15 * tc[4];
    } else {
        /*species 4: ch2 */
        species[4] =
            +3.44310466e+00
            +2.32638011e-03 * tc[1]
            -5.33147998e-07 * tc[2]
            +5.74043968e-11 * tc[3]
            -2.42694849e-15 * tc[4];
    }

    /*species with midpoint at T=1671 kelvin */
    if (T < 1671) {
        /*species 19: c2h3 */
        species[19] =
            +2.73925942e+00
            +7.03011591e-03 * tc[1]
            +2.36647299e-06 * tc[2]
            -3.59569737e-09 * tc[3]
            +8.91756491e-13 * tc[4];
    } else {
        /*species 19: c2h3 */
        species[19] =
            +3.96047713e+00
            +7.99426013e-03 * tc[1]
            -2.85608087e-06 * tc[2]
            +4.58350813e-10 * tc[3]
            -2.72570140e-14 * tc[4];
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: h2 */
        species[0] =
            +3.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 1: h */
        species[1] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 2: ch4 */
        species[2] =
            +7.78741500e-01
            +1.74766800e-02 * tc[1]
            -2.78340900e-05 * tc[2]
            +3.04970800e-08 * tc[3]
            -1.22393100e-11 * tc[4];
        /*species 3: ch3 */
        species[3] =
            +2.43044300e+00
            +1.11241000e-02 * tc[1]
            -1.68022000e-05 * tc[2]
            +1.62182900e-08 * tc[3]
            -5.86495300e-12 * tc[4];
        /*species 5: ch */
        species[5] =
            +3.20020200e+00
            +2.07287600e-03 * tc[1]
            -5.13443100e-06 * tc[2]
            +5.73389000e-09 * tc[3]
            -1.95553300e-12 * tc[4];
        /*species 6: ch2o */
        species[6] =
            +1.65273100e+00
            +1.26314400e-02 * tc[1]
            -1.88816800e-05 * tc[2]
            +2.05003100e-08 * tc[3]
            -8.41323700e-12 * tc[4];
        /*species 7: hco */
        species[7] =
            +2.89833000e+00
            +6.19914700e-03 * tc[1]
            -9.62308400e-06 * tc[2]
            +1.08982500e-08 * tc[3]
            -4.57488500e-12 * tc[4];
        /*species 8: co2 */
        species[8] =
            +2.27572500e+00
            +9.92207200e-03 * tc[1]
            -1.04091100e-05 * tc[2]
            +6.86668700e-09 * tc[3]
            -2.11728000e-12 * tc[4];
        /*species 9: co */
        species[9] =
            +3.26245200e+00
            +1.51194100e-03 * tc[1]
            -3.88175500e-06 * tc[2]
            +5.58194400e-09 * tc[3]
            -2.47495100e-12 * tc[4];
        /*species 10: o2 */
        species[10] =
            +3.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
        /*species 11: o */
        species[11] =
            +2.94642900e+00
            -1.63816600e-03 * tc[1]
            +2.42103200e-06 * tc[2]
            -1.60284300e-09 * tc[3]
            +3.89069600e-13 * tc[4];
        /*species 12: oh */
        species[12] =
            +3.63726600e+00
            +1.85091000e-04 * tc[1]
            -1.67616500e-06 * tc[2]
            +2.38720300e-09 * tc[3]
            -8.43144200e-13 * tc[4];
        /*species 14: h2o2 */
        species[14] =
            +3.38875400e+00
            +6.56922600e-03 * tc[1]
            -1.48501300e-07 * tc[2]
            -4.62580600e-09 * tc[3]
            +2.47151500e-12 * tc[4];
        /*species 15: h2o */
        species[15] =
            +3.38684200e+00
            +3.47498200e-03 * tc[1]
            -6.35469600e-06 * tc[2]
            +6.96858100e-09 * tc[3]
            -2.50658800e-12 * tc[4];
        /*species 16: c2h */
        species[16] =
            +2.73770400e+00
            +8.04844600e-03 * tc[1]
            -9.24431000e-06 * tc[2]
            +6.52525900e-09 * tc[3]
            -1.93958000e-12 * tc[4];
        /*species 17: hcco */
        species[17] =
            +5.04796500e+00
            +4.45347800e-03 * tc[1]
            +2.26828300e-07 * tc[2]
            -1.48209500e-09 * tc[3]
            +2.25074200e-13 * tc[4];
        /*species 18: c2h2 */
        species[18] =
            +2.01356200e+00
            +1.51904500e-02 * tc[1]
            -1.61631900e-05 * tc[2]
            +9.07899200e-09 * tc[3]
            -1.91274600e-12 * tc[4];
        /*species 20: c2h4 */
        species[20] =
            -8.61488000e-01
            +2.79616300e-02 * tc[1]
            -3.38867700e-05 * tc[2]
            +2.78515200e-08 * tc[3]
            -9.73787900e-12 * tc[4];
        /*species 22: c2h6 */
        species[22] =
            +1.46253900e+00
            +1.54946700e-02 * tc[1]
            +5.78050700e-06 * tc[2]
            -1.25783200e-08 * tc[3]
            +4.58626700e-12 * tc[4];
        /*species 24: ch3o */
        species[24] =
            +2.10620400e+00
            +7.21659500e-03 * tc[1]
            +5.33847200e-06 * tc[2]
            -7.37763600e-09 * tc[3]
            +2.07561100e-12 * tc[4];
        /*species 25: hccoh */
        species[25] =
            +3.89946500e+00
            +9.70107500e-03 * tc[1]
            -3.11930900e-07 * tc[2]
            -5.53773200e-09 * tc[3]
            +2.46573200e-12 * tc[4];
        /*species 26: h2ccch */
        species[26] =
            +4.75420000e+00
            +1.10802800e-02 * tc[1]
            +2.79332300e-07 * tc[2]
            -5.47921200e-09 * tc[3]
            +1.94962900e-12 * tc[4];
        /*species 27: c3h2 */
        species[27] =
            +3.16671400e+00
            +2.48257200e-02 * tc[1]
            -4.59163700e-05 * tc[2]
            +4.26801900e-08 * tc[3]
            -1.48215200e-11 * tc[4];
        /*species 29: ch2co */
        species[29] =
            +2.97497100e+00
            +1.21187100e-02 * tc[1]
            -2.34504600e-06 * tc[2]
            -6.46668500e-09 * tc[3]
            +3.90564900e-12 * tc[4];
        /*species 30: ch2hco */
        species[30] =
            +3.40906200e+00
            +1.07385700e-02 * tc[1]
            +1.89149200e-06 * tc[2]
            -7.15858300e-09 * tc[3]
            +2.86738500e-12 * tc[4];
        /*species 31: ch3co */
        species[31] =
            +3.12527800e+00
            +9.77822000e-03 * tc[1]
            +4.52144800e-06 * tc[2]
            -9.00946200e-09 * tc[3]
            +3.19371800e-12 * tc[4];
        /*species 32: ch3hco */
        species[32] =
            +2.50569500e+00
            +1.33699100e-02 * tc[1]
            +4.67195300e-06 * tc[2]
            -1.12814000e-08 * tc[3]
            +4.26356600e-12 * tc[4];
        /*species 33: c2h5oh */
        species[33] =
            +4.85869570e+00
            -3.74017260e-03 * tc[1]
            +6.95553780e-05 * tc[2]
            -8.86547960e-08 * tc[3]
            +3.51688350e-11 * tc[4];
        /*species 37: n2 */
        species[37] =
            +3.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: h2 */
        species[0] =
            +2.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 1: h */
        species[1] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 2: ch4 */
        species[2] =
            +1.68347900e+00
            +1.02372400e-02 * tc[1]
            -3.87512900e-06 * tc[2]
            +6.78558500e-10 * tc[3]
            -4.50342300e-14 * tc[4];
        /*species 3: ch3 */
        species[3] =
            +2.84405200e+00
            +6.13797400e-03 * tc[1]
            -2.23034500e-06 * tc[2]
            +3.78516100e-10 * tc[3]
            -2.45215900e-14 * tc[4];
        /*species 5: ch */
        species[5] =
            +2.19622300e+00
            +2.34038100e-03 * tc[1]
            -7.05820100e-07 * tc[2]
            +9.00758200e-11 * tc[3]
            -3.85504000e-15 * tc[4];
        /*species 6: ch2o */
        species[6] =
            +2.99560600e+00
            +6.68132100e-03 * tc[1]
            -2.62895500e-06 * tc[2]
            +4.73715300e-10 * tc[3]
            -3.21251700e-14 * tc[4];
        /*species 7: hco */
        species[7] =
            +3.55727100e+00
            +3.34557300e-03 * tc[1]
            -1.33500600e-06 * tc[2]
            +2.47057300e-10 * tc[3]
            -1.71385100e-14 * tc[4];
        /*species 8: co2 */
        species[8] =
            +4.45362300e+00
            +3.14016900e-03 * tc[1]
            -1.27841100e-06 * tc[2]
            +2.39399700e-10 * tc[3]
            -1.66903300e-14 * tc[4];
        /*species 9: co */
        species[9] =
            +3.02507800e+00
            +1.44268900e-03 * tc[1]
            -5.63082800e-07 * tc[2]
            +1.01858100e-10 * tc[3]
            -6.91095200e-15 * tc[4];
        /*species 10: o2 */
        species[10] =
            +3.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
        /*species 11: o */
        species[11] =
            +2.54206000e+00
            -2.75506200e-05 * tc[1]
            -3.10280300e-09 * tc[2]
            +4.55106700e-12 * tc[3]
            -4.36805200e-16 * tc[4];
        /*species 12: oh */
        species[12] =
            +2.88273000e+00
            +1.01397400e-03 * tc[1]
            -2.27687700e-07 * tc[2]
            +2.17468400e-11 * tc[3]
            -5.12630500e-16 * tc[4];
        /*species 14: h2o2 */
        species[14] =
            +4.57316700e+00
            +4.33613600e-03 * tc[1]
            -1.47468900e-06 * tc[2]
            +2.34890400e-10 * tc[3]
            -1.43165400e-14 * tc[4];
        /*species 15: h2o */
        species[15] =
            +2.67214600e+00
            +3.05629300e-03 * tc[1]
            -8.73026000e-07 * tc[2]
            +1.20099600e-10 * tc[3]
            -6.39161800e-15 * tc[4];
        /*species 16: c2h */
        species[16] =
            +3.98636700e+00
            +3.14312300e-03 * tc[1]
            -1.26724300e-06 * tc[2]
            +2.92436300e-10 * tc[3]
            -2.71632000e-14 * tc[4];
        /*species 17: hcco */
        species[17] =
            +6.75807300e+00
            +2.00040000e-03 * tc[1]
            -2.02760700e-07 * tc[2]
            -1.04113200e-10 * tc[3]
            +1.96516500e-14 * tc[4];
        /*species 18: c2h2 */
        species[18] =
            +4.43677000e+00
            +5.37603900e-03 * tc[1]
            -1.91281700e-06 * tc[2]
            +3.28637900e-10 * tc[3]
            -2.15671000e-14 * tc[4];
        /*species 20: c2h4 */
        species[20] =
            +3.52841900e+00
            +1.14851800e-02 * tc[1]
            -4.41838500e-06 * tc[2]
            +7.84460100e-10 * tc[3]
            -5.26684800e-14 * tc[4];
        /*species 22: c2h6 */
        species[22] =
            +4.82593800e+00
            +1.38404300e-02 * tc[1]
            -4.55725900e-06 * tc[2]
            +6.72496700e-10 * tc[3]
            -3.59816100e-14 * tc[4];
        /*species 24: ch3o */
        species[24] =
            +3.77080000e+00
            +7.87149700e-03 * tc[1]
            -2.65638400e-06 * tc[2]
            +3.94443100e-10 * tc[3]
            -2.11261600e-14 * tc[4];
        /*species 25: hccoh */
        species[25] =
            +7.32832400e+00
            +3.33641600e-03 * tc[1]
            -3.02470500e-07 * tc[2]
            -1.78110600e-10 * tc[3]
            +3.24516800e-14 * tc[4];
        /*species 26: h2ccch */
        species[26] =
            +8.83104700e+00
            +4.35719500e-03 * tc[1]
            -4.10906700e-07 * tc[2]
            -2.36872300e-10 * tc[3]
            +4.37652000e-14 * tc[4];
        /*species 27: c3h2 */
        species[27] =
            +7.67098100e+00
            +2.74874900e-03 * tc[1]
            -4.37094300e-07 * tc[2]
            -6.45559900e-11 * tc[3]
            +1.66388700e-14 * tc[4];
        /*species 29: ch2co */
        species[29] =
            +6.03881700e+00
            +5.80484000e-03 * tc[1]
            -1.92095400e-06 * tc[2]
            +2.79448500e-10 * tc[3]
            -1.45886800e-14 * tc[4];
        /*species 30: ch2hco */
        species[30] =
            +5.97567000e+00
            +8.13059100e-03 * tc[1]
            -2.74362400e-06 * tc[2]
            +4.07030400e-10 * tc[3]
            -2.17601700e-14 * tc[4];
        /*species 31: ch3co */
        species[31] =
            +5.61227900e+00
            +8.44988600e-03 * tc[1]
            -2.85414700e-06 * tc[2]
            +4.23837600e-10 * tc[3]
            -2.26840400e-14 * tc[4];
        /*species 32: ch3hco */
        species[32] =
            +5.86865000e+00
            +1.07942400e-02 * tc[1]
            -3.64553000e-06 * tc[2]
            +5.41291200e-10 * tc[3]
            -2.89684400e-14 * tc[4];
        /*species 33: c2h5oh */
        species[33] =
            +6.56243650e+00
            +1.52042220e-02 * tc[1]
            -5.38967950e-06 * tc[2]
            +8.62250110e-10 * tc[3]
            -5.12897870e-14 * tc[4];
        /*species 37: n2 */
        species[37] =
            +2.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }

    /*species with midpoint at T=1451 kelvin */
    if (T < 1451) {
        /*species 13: ho2 */
        species[13] =
            +3.47629499e+00
            +2.20467745e-03 * tc[1]
            +1.56841876e-06 * tc[2]
            -2.12755802e-09 * tc[3]
            +5.83140718e-13 * tc[4];
    } else {
        /*species 13: ho2 */
        species[13] =
            +4.58311499e+00
            +1.72730339e-03 * tc[1]
            -6.19183136e-07 * tc[2]
            +9.91939193e-11 * tc[3]
            -5.88221567e-15 * tc[4];
    }

    /*species with midpoint at T=1359 kelvin */
    if (T < 1359) {
        /*species 28: ch2s */
        species[28] =
            +3.32340701e+00
            +2.28062579e-03 * tc[1]
            -2.48156016e-07 * tc[2]
            -1.02065386e-10 * tc[3]
            +2.48826619e-14 * tc[4];
    } else {
        /*species 28: ch2s */
        species[28] =
            +3.09627180e+00
            +2.80839185e-03 * tc[1]
            -7.13614526e-07 * tc[2]
            +8.42136298e-11 * tc[3]
            -3.84640828e-15 * tc[4];
    }

    /*species with midpoint at T=1524 kelvin */
    if (T < 1524) {
        /*species 35: ch3choh */
        species[35] =
            +2.10026806e+00
            +1.92626679e-02 * tc[1]
            -5.17437436e-06 * tc[2]
            -1.96190582e-09 * tc[3]
            +9.36397060e-13 * tc[4];
    } else {
        /*species 35: ch3choh */
        species[35] =
            +7.95432747e+00
            +1.05008170e-02 * tc[1]
            -3.50631590e-06 * tc[2]
            +5.36174608e-10 * tc[3]
            -3.07978459e-14 * tc[4];
    }

    /*species with midpoint at T=1405 kelvin */
    if (T < 1405) {
        /*species 34: c2h4oh */
        species[34] =
            +7.76157606e-01
            +2.91002514e-02 * tc[1]
            -2.13980569e-05 * tc[2]
            +8.16833894e-09 * tc[3]
            -1.25591600e-12 * tc[4];
        /*species 36: ch3ch2o */
        species[36] =
            -2.71296378e-01
            +2.98839812e-02 * tc[1]
            -1.97090548e-05 * tc[2]
            +6.37339893e-09 * tc[3]
            -7.77965054e-13 * tc[4];
    } else {
        /*species 34: c2h4oh */
        species[34] =
            +8.75544496e+00
            +9.49453147e-03 * tc[1]
            -3.08803575e-06 * tc[2]
            +4.63537621e-10 * tc[3]
            -2.62730044e-14 * tc[4];
        /*species 36: ch3ch2o */
        species[36] =
            +8.31182392e+00
            +1.03426319e-02 * tc[1]
            -3.39186089e-06 * tc[2]
            +5.12212617e-10 * tc[3]
            -2.91601713e-14 * tc[4];
    }

    /*species with midpoint at T=1375 kelvin */
    if (T < 1375) {
        /*species 21: c2h5 */
        species[21] =
            +1.47374852e+00
            +1.63606611e-02 * tc[1]
            -4.32896826e-06 * tc[2]
            -1.17270665e-09 * tc[3]
            +5.72325883e-13 * tc[4];
    } else {
        /*species 21: c2h5 */
        species[21] =
            +5.60116091e+00
            +1.06977081e-02 * tc[1]
            -3.63504780e-06 * tc[2]
            +5.61816654e-10 * tc[3]
            -3.24914714e-14 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesInternalEnergy(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1410 kelvin */
    if (T < 1410) {
        /*species 23: ch2oh */
        species[23] =
            +1.60067849e+00
            +6.42239270e-03 * tc[1]
            -2.77932062e-06 * tc[2]
            +6.89319015e-10 * tc[3]
            -7.14082212e-14 * tc[4]
            -2.33478724e+03 / tc[1];
    } else {
        /*species 23: ch2oh */
        species[23] =
            +5.00127803e+00
            +2.49360784e-03 * tc[1]
            -5.36511717e-07 * tc[2]
            +6.00567838e-11 * tc[3]
            -2.71165400e-15 * tc[4]
            -3.50157098e+03 / tc[1];
    }

    /*species with midpoint at T=1350 kelvin */
    if (T < 1350) {
        /*species 4: ch2 */
        species[4] =
            +2.49623989e+00
            +1.10129322e-03 * tc[1]
            -1.41276266e-07 * tc[2]
            +3.41211765e-12 * tc[3]
            +8.63874744e-16 * tc[4]
            +4.59192771e+04 / tc[1];
    } else {
        /*species 4: ch2 */
        species[4] =
            +2.44310466e+00
            +1.16319005e-03 * tc[1]
            -1.77715999e-07 * tc[2]
            +1.43510992e-11 * tc[3]
            -4.85389698e-16 * tc[4]
            +4.59375743e+04 / tc[1];
    }

    /*species with midpoint at T=1671 kelvin */
    if (T < 1671) {
        /*species 19: c2h3 */
        species[19] =
            +1.73925942e+00
            +3.51505795e-03 * tc[1]
            +7.88824330e-07 * tc[2]
            -8.98924342e-10 * tc[3]
            +1.78351298e-13 * tc[4]
            +3.42868979e+04 / tc[1];
    } else {
        /*species 19: c2h3 */
        species[19] =
            +2.96047713e+00
            +3.99713006e-03 * tc[1]
            -9.52026957e-07 * tc[2]
            +1.14587703e-10 * tc[3]
            -5.45140280e-15 * tc[4]
            +3.35153544e+04 / tc[1];
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: h2 */
        species[0] =
            +2.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 / tc[1];
        /*species 1: h */
        species[1] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 / tc[1];
        /*species 2: ch4 */
        species[2] =
            -2.21258500e-01
            +8.73834000e-03 * tc[1]
            -9.27803000e-06 * tc[2]
            +7.62427000e-09 * tc[3]
            -2.44786200e-12 * tc[4]
            -9.82522900e+03 / tc[1];
        /*species 3: ch3 */
        species[3] =
            +1.43044300e+00
            +5.56205000e-03 * tc[1]
            -5.60073333e-06 * tc[2]
            +4.05457250e-09 * tc[3]
            -1.17299060e-12 * tc[4]
            +1.64237800e+04 / tc[1];
        /*species 5: ch */
        species[5] =
            +2.20020200e+00
            +1.03643800e-03 * tc[1]
            -1.71147700e-06 * tc[2]
            +1.43347250e-09 * tc[3]
            -3.91106600e-13 * tc[4]
            +7.04525900e+04 / tc[1];
        /*species 6: ch2o */
        species[6] =
            +6.52731000e-01
            +6.31572000e-03 * tc[1]
            -6.29389333e-06 * tc[2]
            +5.12507750e-09 * tc[3]
            -1.68264740e-12 * tc[4]
            -1.48654000e+04 / tc[1];
        /*species 7: hco */
        species[7] =
            +1.89833000e+00
            +3.09957350e-03 * tc[1]
            -3.20769467e-06 * tc[2]
            +2.72456250e-09 * tc[3]
            -9.14977000e-13 * tc[4]
            +4.15992200e+03 / tc[1];
        /*species 8: co2 */
        species[8] =
            +1.27572500e+00
            +4.96103600e-03 * tc[1]
            -3.46970333e-06 * tc[2]
            +1.71667175e-09 * tc[3]
            -4.23456000e-13 * tc[4]
            -4.83731400e+04 / tc[1];
        /*species 9: co */
        species[9] =
            +2.26245200e+00
            +7.55970500e-04 * tc[1]
            -1.29391833e-06 * tc[2]
            +1.39548600e-09 * tc[3]
            -4.94990200e-13 * tc[4]
            -1.43105400e+04 / tc[1];
        /*species 10: o2 */
        species[10] =
            +2.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 / tc[1];
        /*species 11: o */
        species[11] =
            +1.94642900e+00
            -8.19083000e-04 * tc[1]
            +8.07010667e-07 * tc[2]
            -4.00710750e-10 * tc[3]
            +7.78139200e-14 * tc[4]
            +2.91476400e+04 / tc[1];
        /*species 12: oh */
        species[12] =
            +2.63726600e+00
            +9.25455000e-05 * tc[1]
            -5.58721667e-07 * tc[2]
            +5.96800750e-10 * tc[3]
            -1.68628840e-13 * tc[4]
            +3.60678200e+03 / tc[1];
        /*species 14: h2o2 */
        species[14] =
            +2.38875400e+00
            +3.28461300e-03 * tc[1]
            -4.95004333e-08 * tc[2]
            -1.15645150e-09 * tc[3]
            +4.94303000e-13 * tc[4]
            -1.76631500e+04 / tc[1];
        /*species 15: h2o */
        species[15] =
            +2.38684200e+00
            +1.73749100e-03 * tc[1]
            -2.11823200e-06 * tc[2]
            +1.74214525e-09 * tc[3]
            -5.01317600e-13 * tc[4]
            -3.02081100e+04 / tc[1];
        /*species 16: c2h */
        species[16] =
            +1.73770400e+00
            +4.02422300e-03 * tc[1]
            -3.08143667e-06 * tc[2]
            +1.63131475e-09 * tc[3]
            -3.87916000e-13 * tc[4]
            +6.68381300e+04 / tc[1];
        /*species 17: hcco */
        species[17] =
            +4.04796500e+00
            +2.22673900e-03 * tc[1]
            +7.56094333e-08 * tc[2]
            -3.70523750e-10 * tc[3]
            +4.50148400e-14 * tc[4]
            +1.96589200e+04 / tc[1];
        /*species 18: c2h2 */
        species[18] =
            +1.01356200e+00
            +7.59522500e-03 * tc[1]
            -5.38773000e-06 * tc[2]
            +2.26974800e-09 * tc[3]
            -3.82549200e-13 * tc[4]
            +2.61244400e+04 / tc[1];
        /*species 20: c2h4 */
        species[20] =
            -1.86148800e+00
            +1.39808150e-02 * tc[1]
            -1.12955900e-05 * tc[2]
            +6.96288000e-09 * tc[3]
            -1.94757580e-12 * tc[4]
            +5.57304600e+03 / tc[1];
        /*species 22: c2h6 */
        species[22] =
            +4.62539000e-01
            +7.74733500e-03 * tc[1]
            +1.92683567e-06 * tc[2]
            -3.14458000e-09 * tc[3]
            +9.17253400e-13 * tc[4]
            -1.12391800e+04 / tc[1];
        /*species 24: ch3o */
        species[24] =
            +1.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122200e-13 * tc[4]
            +9.78601100e+02 / tc[1];
        /*species 25: hccoh */
        species[25] =
            +2.89946500e+00
            +4.85053750e-03 * tc[1]
            -1.03976967e-07 * tc[2]
            -1.38443300e-09 * tc[3]
            +4.93146400e-13 * tc[4]
            +8.70119000e+03 / tc[1];
        /*species 26: h2ccch */
        species[26] =
            +3.75420000e+00
            +5.54014000e-03 * tc[1]
            +9.31107667e-08 * tc[2]
            -1.36980300e-09 * tc[3]
            +3.89925800e-13 * tc[4]
            +3.98888300e+04 / tc[1];
        /*species 27: c3h2 */
        species[27] =
            +2.16671400e+00
            +1.24128600e-02 * tc[1]
            -1.53054567e-05 * tc[2]
            +1.06700475e-08 * tc[3]
            -2.96430400e-12 * tc[4]
            +6.35042100e+04 / tc[1];
        /*species 29: ch2co */
        species[29] =
            +1.97497100e+00
            +6.05935500e-03 * tc[1]
            -7.81682000e-07 * tc[2]
            -1.61667125e-09 * tc[3]
            +7.81129800e-13 * tc[4]
            -7.63263700e+03 / tc[1];
        /*species 30: ch2hco */
        species[30] =
            +2.40906200e+00
            +5.36928500e-03 * tc[1]
            +6.30497333e-07 * tc[2]
            -1.78964575e-09 * tc[3]
            +5.73477000e-13 * tc[4]
            +1.52147700e+03 / tc[1];
        /*species 31: ch3co */
        species[31] =
            +2.12527800e+00
            +4.88911000e-03 * tc[1]
            +1.50714933e-06 * tc[2]
            -2.25236550e-09 * tc[3]
            +6.38743600e-13 * tc[4]
            -4.10850800e+03 / tc[1];
        /*species 32: ch3hco */
        species[32] =
            +1.50569500e+00
            +6.68495500e-03 * tc[1]
            +1.55731767e-06 * tc[2]
            -2.82035000e-09 * tc[3]
            +8.52713200e-13 * tc[4]
            -2.12458900e+04 / tc[1];
        /*species 33: c2h5oh */
        species[33] =
            +3.85869570e+00
            -1.87008630e-03 * tc[1]
            +2.31851260e-05 * tc[2]
            -2.21636990e-08 * tc[3]
            +7.03376700e-12 * tc[4]
            -2.99961320e+04 / tc[1];
        /*species 37: n2 */
        species[37] =
            +2.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 / tc[1];
    } else {
        /*species 0: h2 */
        species[0] =
            +1.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 / tc[1];
        /*species 1: h */
        species[1] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 / tc[1];
        /*species 2: ch4 */
        species[2] =
            +6.83479000e-01
            +5.11862000e-03 * tc[1]
            -1.29170967e-06 * tc[2]
            +1.69639625e-10 * tc[3]
            -9.00684600e-15 * tc[4]
            -1.00807900e+04 / tc[1];
        /*species 3: ch3 */
        species[3] =
            +1.84405200e+00
            +3.06898700e-03 * tc[1]
            -7.43448333e-07 * tc[2]
            +9.46290250e-11 * tc[3]
            -4.90431800e-15 * tc[4]
            +1.64378100e+04 / tc[1];
        /*species 5: ch */
        species[5] =
            +1.19622300e+00
            +1.17019050e-03 * tc[1]
            -2.35273367e-07 * tc[2]
            +2.25189550e-11 * tc[3]
            -7.71008000e-16 * tc[4]
            +7.08672300e+04 / tc[1];
        /*species 6: ch2o */
        species[6] =
            +1.99560600e+00
            +3.34066050e-03 * tc[1]
            -8.76318333e-07 * tc[2]
            +1.18428825e-10 * tc[3]
            -6.42503400e-15 * tc[4]
            -1.53203700e+04 / tc[1];
        /*species 7: hco */
        species[7] =
            +2.55727100e+00
            +1.67278650e-03 * tc[1]
            -4.45002000e-07 * tc[2]
            +6.17643250e-11 * tc[3]
            -3.42770200e-15 * tc[4]
            +3.91632400e+03 / tc[1];
        /*species 8: co2 */
        species[8] =
            +3.45362300e+00
            +1.57008450e-03 * tc[1]
            -4.26137000e-07 * tc[2]
            +5.98499250e-11 * tc[3]
            -3.33806600e-15 * tc[4]
            -4.89669600e+04 / tc[1];
        /*species 9: co */
        species[9] =
            +2.02507800e+00
            +7.21344500e-04 * tc[1]
            -1.87694267e-07 * tc[2]
            +2.54645250e-11 * tc[3]
            -1.38219040e-15 * tc[4]
            -1.42683500e+04 / tc[1];
        /*species 10: o2 */
        species[10] =
            +2.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 / tc[1];
        /*species 11: o */
        species[11] =
            +1.54206000e+00
            -1.37753100e-05 * tc[1]
            -1.03426767e-09 * tc[2]
            +1.13776675e-12 * tc[3]
            -8.73610400e-17 * tc[4]
            +2.92308000e+04 / tc[1];
        /*species 12: oh */
        species[12] =
            +1.88273000e+00
            +5.06987000e-04 * tc[1]
            -7.58959000e-08 * tc[2]
            +5.43671000e-12 * tc[3]
            -1.02526100e-16 * tc[4]
            +3.88688800e+03 / tc[1];
        /*species 14: h2o2 */
        species[14] =
            +3.57316700e+00
            +2.16806800e-03 * tc[1]
            -4.91563000e-07 * tc[2]
            +5.87226000e-11 * tc[3]
            -2.86330800e-15 * tc[4]
            -1.80069600e+04 / tc[1];
        /*species 15: h2o */
        species[15] =
            +1.67214600e+00
            +1.52814650e-03 * tc[1]
            -2.91008667e-07 * tc[2]
            +3.00249000e-11 * tc[3]
            -1.27832360e-15 * tc[4]
            -2.98992100e+04 / tc[1];
        /*species 16: c2h */
        species[16] =
            +2.98636700e+00
            +1.57156150e-03 * tc[1]
            -4.22414333e-07 * tc[2]
            +7.31090750e-11 * tc[3]
            -5.43264000e-15 * tc[4]
            +6.65588400e+04 / tc[1];
        /*species 17: hcco */
        species[17] =
            +5.75807300e+00
            +1.00020000e-03 * tc[1]
            -6.75869000e-08 * tc[2]
            -2.60283000e-11 * tc[3]
            +3.93033000e-15 * tc[4]
            +1.90151300e+04 / tc[1];
        /*species 18: c2h2 */
        species[18] =
            +3.43677000e+00
            +2.68801950e-03 * tc[1]
            -6.37605667e-07 * tc[2]
            +8.21594750e-11 * tc[3]
            -4.31342000e-15 * tc[4]
            +2.56676600e+04 / tc[1];
        /*species 20: c2h4 */
        species[20] =
            +2.52841900e+00
            +5.74259000e-03 * tc[1]
            -1.47279500e-06 * tc[2]
            +1.96115025e-10 * tc[3]
            -1.05336960e-14 * tc[4]
            +4.42828900e+03 / tc[1];
        /*species 22: c2h6 */
        species[22] =
            +3.82593800e+00
            +6.92021500e-03 * tc[1]
            -1.51908633e-06 * tc[2]
            +1.68124175e-10 * tc[3]
            -7.19632200e-15 * tc[4]
            -1.27177900e+04 / tc[1];
        /*species 24: ch3o */
        species[24] =
            +2.77080000e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832500e+02 / tc[1];
        /*species 25: hccoh */
        species[25] =
            +6.32832400e+00
            +1.66820800e-03 * tc[1]
            -1.00823500e-07 * tc[2]
            -4.45276500e-11 * tc[3]
            +6.49033600e-15 * tc[4]
            +7.59825800e+03 / tc[1];
        /*species 26: h2ccch */
        species[26] =
            +7.83104700e+00
            +2.17859750e-03 * tc[1]
            -1.36968900e-07 * tc[2]
            -5.92180750e-11 * tc[3]
            +8.75304000e-15 * tc[4]
            +3.84742000e+04 / tc[1];
        /*species 27: c3h2 */
        species[27] =
            +6.67098100e+00
            +1.37437450e-03 * tc[1]
            -1.45698100e-07 * tc[2]
            -1.61389975e-11 * tc[3]
            +3.32777400e-15 * tc[4]
            +6.25972200e+04 / tc[1];
        /*species 29: ch2co */
        species[29] =
            +5.03881700e+00
            +2.90242000e-03 * tc[1]
            -6.40318000e-07 * tc[2]
            +6.98621250e-11 * tc[3]
            -2.91773600e-15 * tc[4]
            -8.58340200e+03 / tc[1];
        /*species 30: ch2hco */
        species[30] =
            +4.97567000e+00
            +4.06529550e-03 * tc[1]
            -9.14541333e-07 * tc[2]
            +1.01757600e-10 * tc[3]
            -4.35203400e-15 * tc[4]
            +4.90321800e+02 / tc[1];
        /*species 31: ch3co */
        species[31] =
            +4.61227900e+00
            +4.22494300e-03 * tc[1]
            -9.51382333e-07 * tc[2]
            +1.05959400e-10 * tc[3]
            -4.53680800e-15 * tc[4]
            -5.18786300e+03 / tc[1];
        /*species 32: ch3hco */
        species[32] =
            +4.86865000e+00
            +5.39712000e-03 * tc[1]
            -1.21517667e-06 * tc[2]
            +1.35322800e-10 * tc[3]
            -5.79368800e-15 * tc[4]
            -2.26456900e+04 / tc[1];
        /*species 33: c2h5oh */
        species[33] =
            +5.56243650e+00
            +7.60211100e-03 * tc[1]
            -1.79655983e-06 * tc[2]
            +2.15562527e-10 * tc[3]
            -1.02579574e-14 * tc[4]
            -3.15256210e+04 / tc[1];
        /*species 37: n2 */
        species[37] =
            +1.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 / tc[1];
    }

    /*species with midpoint at T=1451 kelvin */
    if (T < 1451) {
        /*species 13: ho2 */
        species[13] =
            +2.47629499e+00
            +1.10233872e-03 * tc[1]
            +5.22806253e-07 * tc[2]
            -5.31889505e-10 * tc[3]
            +1.16628144e-13 * tc[4]
            +6.17073446e+02 / tc[1];
    } else {
        /*species 13: ho2 */
        species[13] =
            +3.58311499e+00
            +8.63651695e-04 * tc[1]
            -2.06394379e-07 * tc[2]
            +2.47984798e-11 * tc[3]
            -1.17644313e-15 * tc[4]
            +3.13365859e+01 / tc[1];
    }

    /*species with midpoint at T=1359 kelvin */
    if (T < 1359) {
        /*species 28: ch2s */
        species[28] =
            +2.32340701e+00
            +1.14031289e-03 * tc[1]
            -8.27186720e-08 * tc[2]
            -2.55163465e-11 * tc[3]
            +4.97653238e-15 * tc[4]
            +5.04954469e+04 / tc[1];
    } else {
        /*species 28: ch2s */
        species[28] =
            +2.09627180e+00
            +1.40419592e-03 * tc[1]
            -2.37871509e-07 * tc[2]
            +2.10534074e-11 * tc[3]
            -7.69281656e-16 * tc[4]
            +5.05739697e+04 / tc[1];
    }

    /*species with midpoint at T=1524 kelvin */
    if (T < 1524) {
        /*species 35: ch3choh */
        species[35] =
            +1.10026806e+00
            +9.63133395e-03 * tc[1]
            -1.72479145e-06 * tc[2]
            -4.90476455e-10 * tc[3]
            +1.87279412e-13 * tc[4]
            -6.12873036e+03 / tc[1];
    } else {
        /*species 35: ch3choh */
        species[35] =
            +6.95432747e+00
            +5.25040850e-03 * tc[1]
            -1.16877197e-06 * tc[2]
            +1.34043652e-10 * tc[3]
            -6.15956918e-15 * tc[4]
            -8.62199768e+03 / tc[1];
    }

    /*species with midpoint at T=1405 kelvin */
    if (T < 1405) {
        /*species 34: c2h4oh */
        species[34] =
            -2.23842394e-01
            +1.45501257e-02 * tc[1]
            -7.13268563e-06 * tc[2]
            +2.04208474e-09 * tc[3]
            -2.51183200e-13 * tc[4]
            -5.46808396e+03 / tc[1];
        /*species 36: ch3ch2o */
        species[36] =
            -1.27129638e+00
            +1.49419906e-02 * tc[1]
            -6.56968493e-06 * tc[2]
            +1.59334973e-09 * tc[3]
            -1.55593011e-13 * tc[4]
            -3.16397196e+03 / tc[1];
    } else {
        /*species 34: c2h4oh */
        species[34] =
            +7.75544496e+00
            +4.74726573e-03 * tc[1]
            -1.02934525e-06 * tc[2]
            +1.15884405e-10 * tc[3]
            -5.25460088e-15 * tc[4]
            -8.09602047e+03 / tc[1];
        /*species 36: ch3ch2o */
        species[36] =
            +7.31182392e+00
            +5.17131595e-03 * tc[1]
            -1.13062030e-06 * tc[2]
            +1.28053154e-10 * tc[3]
            -5.83203426e-15 * tc[4]
            -6.13097954e+03 / tc[1];
    }

    /*species with midpoint at T=1375 kelvin */
    if (T < 1375) {
        /*species 21: c2h5 */
        species[21] =
            +4.73748520e-01
            +8.18033055e-03 * tc[1]
            -1.44298942e-06 * tc[2]
            -2.93176662e-10 * tc[3]
            +1.14465177e-13 * tc[4]
            +1.33326992e+04 / tc[1];
    } else {
        /*species 21: c2h5 */
        species[21] =
            +4.60116091e+00
            +5.34885405e-03 * tc[1]
            -1.21168260e-06 * tc[2]
            +1.40454163e-10 * tc[3]
            -6.49829428e-15 * tc[4]
            +1.14539845e+04 / tc[1];
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEnthalpy(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1410 kelvin */
    if (T < 1410) {
        /*species 23: ch2oh */
        species[23] =
            +2.60067849e+00
            +6.42239270e-03 * tc[1]
            -2.77932062e-06 * tc[2]
            +6.89319015e-10 * tc[3]
            -7.14082212e-14 * tc[4]
            -2.33478724e+03 / tc[1];
    } else {
        /*species 23: ch2oh */
        species[23] =
            +6.00127803e+00
            +2.49360784e-03 * tc[1]
            -5.36511717e-07 * tc[2]
            +6.00567838e-11 * tc[3]
            -2.71165400e-15 * tc[4]
            -3.50157098e+03 / tc[1];
    }

    /*species with midpoint at T=1350 kelvin */
    if (T < 1350) {
        /*species 4: ch2 */
        species[4] =
            +3.49623989e+00
            +1.10129322e-03 * tc[1]
            -1.41276266e-07 * tc[2]
            +3.41211765e-12 * tc[3]
            +8.63874744e-16 * tc[4]
            +4.59192771e+04 / tc[1];
    } else {
        /*species 4: ch2 */
        species[4] =
            +3.44310466e+00
            +1.16319005e-03 * tc[1]
            -1.77715999e-07 * tc[2]
            +1.43510992e-11 * tc[3]
            -4.85389698e-16 * tc[4]
            +4.59375743e+04 / tc[1];
    }

    /*species with midpoint at T=1671 kelvin */
    if (T < 1671) {
        /*species 19: c2h3 */
        species[19] =
            +2.73925942e+00
            +3.51505795e-03 * tc[1]
            +7.88824330e-07 * tc[2]
            -8.98924342e-10 * tc[3]
            +1.78351298e-13 * tc[4]
            +3.42868979e+04 / tc[1];
    } else {
        /*species 19: c2h3 */
        species[19] =
            +3.96047713e+00
            +3.99713006e-03 * tc[1]
            -9.52026957e-07 * tc[2]
            +1.14587703e-10 * tc[3]
            -5.45140280e-15 * tc[4]
            +3.35153544e+04 / tc[1];
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: h2 */
        species[0] =
            +3.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 / tc[1];
        /*species 1: h */
        species[1] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 / tc[1];
        /*species 2: ch4 */
        species[2] =
            +7.78741500e-01
            +8.73834000e-03 * tc[1]
            -9.27803000e-06 * tc[2]
            +7.62427000e-09 * tc[3]
            -2.44786200e-12 * tc[4]
            -9.82522900e+03 / tc[1];
        /*species 3: ch3 */
        species[3] =
            +2.43044300e+00
            +5.56205000e-03 * tc[1]
            -5.60073333e-06 * tc[2]
            +4.05457250e-09 * tc[3]
            -1.17299060e-12 * tc[4]
            +1.64237800e+04 / tc[1];
        /*species 5: ch */
        species[5] =
            +3.20020200e+00
            +1.03643800e-03 * tc[1]
            -1.71147700e-06 * tc[2]
            +1.43347250e-09 * tc[3]
            -3.91106600e-13 * tc[4]
            +7.04525900e+04 / tc[1];
        /*species 6: ch2o */
        species[6] =
            +1.65273100e+00
            +6.31572000e-03 * tc[1]
            -6.29389333e-06 * tc[2]
            +5.12507750e-09 * tc[3]
            -1.68264740e-12 * tc[4]
            -1.48654000e+04 / tc[1];
        /*species 7: hco */
        species[7] =
            +2.89833000e+00
            +3.09957350e-03 * tc[1]
            -3.20769467e-06 * tc[2]
            +2.72456250e-09 * tc[3]
            -9.14977000e-13 * tc[4]
            +4.15992200e+03 / tc[1];
        /*species 8: co2 */
        species[8] =
            +2.27572500e+00
            +4.96103600e-03 * tc[1]
            -3.46970333e-06 * tc[2]
            +1.71667175e-09 * tc[3]
            -4.23456000e-13 * tc[4]
            -4.83731400e+04 / tc[1];
        /*species 9: co */
        species[9] =
            +3.26245200e+00
            +7.55970500e-04 * tc[1]
            -1.29391833e-06 * tc[2]
            +1.39548600e-09 * tc[3]
            -4.94990200e-13 * tc[4]
            -1.43105400e+04 / tc[1];
        /*species 10: o2 */
        species[10] =
            +3.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 / tc[1];
        /*species 11: o */
        species[11] =
            +2.94642900e+00
            -8.19083000e-04 * tc[1]
            +8.07010667e-07 * tc[2]
            -4.00710750e-10 * tc[3]
            +7.78139200e-14 * tc[4]
            +2.91476400e+04 / tc[1];
        /*species 12: oh */
        species[12] =
            +3.63726600e+00
            +9.25455000e-05 * tc[1]
            -5.58721667e-07 * tc[2]
            +5.96800750e-10 * tc[3]
            -1.68628840e-13 * tc[4]
            +3.60678200e+03 / tc[1];
        /*species 14: h2o2 */
        species[14] =
            +3.38875400e+00
            +3.28461300e-03 * tc[1]
            -4.95004333e-08 * tc[2]
            -1.15645150e-09 * tc[3]
            +4.94303000e-13 * tc[4]
            -1.76631500e+04 / tc[1];
        /*species 15: h2o */
        species[15] =
            +3.38684200e+00
            +1.73749100e-03 * tc[1]
            -2.11823200e-06 * tc[2]
            +1.74214525e-09 * tc[3]
            -5.01317600e-13 * tc[4]
            -3.02081100e+04 / tc[1];
        /*species 16: c2h */
        species[16] =
            +2.73770400e+00
            +4.02422300e-03 * tc[1]
            -3.08143667e-06 * tc[2]
            +1.63131475e-09 * tc[3]
            -3.87916000e-13 * tc[4]
            +6.68381300e+04 / tc[1];
        /*species 17: hcco */
        species[17] =
            +5.04796500e+00
            +2.22673900e-03 * tc[1]
            +7.56094333e-08 * tc[2]
            -3.70523750e-10 * tc[3]
            +4.50148400e-14 * tc[4]
            +1.96589200e+04 / tc[1];
        /*species 18: c2h2 */
        species[18] =
            +2.01356200e+00
            +7.59522500e-03 * tc[1]
            -5.38773000e-06 * tc[2]
            +2.26974800e-09 * tc[3]
            -3.82549200e-13 * tc[4]
            +2.61244400e+04 / tc[1];
        /*species 20: c2h4 */
        species[20] =
            -8.61488000e-01
            +1.39808150e-02 * tc[1]
            -1.12955900e-05 * tc[2]
            +6.96288000e-09 * tc[3]
            -1.94757580e-12 * tc[4]
            +5.57304600e+03 / tc[1];
        /*species 22: c2h6 */
        species[22] =
            +1.46253900e+00
            +7.74733500e-03 * tc[1]
            +1.92683567e-06 * tc[2]
            -3.14458000e-09 * tc[3]
            +9.17253400e-13 * tc[4]
            -1.12391800e+04 / tc[1];
        /*species 24: ch3o */
        species[24] =
            +2.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122200e-13 * tc[4]
            +9.78601100e+02 / tc[1];
        /*species 25: hccoh */
        species[25] =
            +3.89946500e+00
            +4.85053750e-03 * tc[1]
            -1.03976967e-07 * tc[2]
            -1.38443300e-09 * tc[3]
            +4.93146400e-13 * tc[4]
            +8.70119000e+03 / tc[1];
        /*species 26: h2ccch */
        species[26] =
            +4.75420000e+00
            +5.54014000e-03 * tc[1]
            +9.31107667e-08 * tc[2]
            -1.36980300e-09 * tc[3]
            +3.89925800e-13 * tc[4]
            +3.98888300e+04 / tc[1];
        /*species 27: c3h2 */
        species[27] =
            +3.16671400e+00
            +1.24128600e-02 * tc[1]
            -1.53054567e-05 * tc[2]
            +1.06700475e-08 * tc[3]
            -2.96430400e-12 * tc[4]
            +6.35042100e+04 / tc[1];
        /*species 29: ch2co */
        species[29] =
            +2.97497100e+00
            +6.05935500e-03 * tc[1]
            -7.81682000e-07 * tc[2]
            -1.61667125e-09 * tc[3]
            +7.81129800e-13 * tc[4]
            -7.63263700e+03 / tc[1];
        /*species 30: ch2hco */
        species[30] =
            +3.40906200e+00
            +5.36928500e-03 * tc[1]
            +6.30497333e-07 * tc[2]
            -1.78964575e-09 * tc[3]
            +5.73477000e-13 * tc[4]
            +1.52147700e+03 / tc[1];
        /*species 31: ch3co */
        species[31] =
            +3.12527800e+00
            +4.88911000e-03 * tc[1]
            +1.50714933e-06 * tc[2]
            -2.25236550e-09 * tc[3]
            +6.38743600e-13 * tc[4]
            -4.10850800e+03 / tc[1];
        /*species 32: ch3hco */
        species[32] =
            +2.50569500e+00
            +6.68495500e-03 * tc[1]
            +1.55731767e-06 * tc[2]
            -2.82035000e-09 * tc[3]
            +8.52713200e-13 * tc[4]
            -2.12458900e+04 / tc[1];
        /*species 33: c2h5oh */
        species[33] =
            +4.85869570e+00
            -1.87008630e-03 * tc[1]
            +2.31851260e-05 * tc[2]
            -2.21636990e-08 * tc[3]
            +7.03376700e-12 * tc[4]
            -2.99961320e+04 / tc[1];
        /*species 37: n2 */
        species[37] =
            +3.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 / tc[1];
    } else {
        /*species 0: h2 */
        species[0] =
            +2.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 / tc[1];
        /*species 1: h */
        species[1] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 / tc[1];
        /*species 2: ch4 */
        species[2] =
            +1.68347900e+00
            +5.11862000e-03 * tc[1]
            -1.29170967e-06 * tc[2]
            +1.69639625e-10 * tc[3]
            -9.00684600e-15 * tc[4]
            -1.00807900e+04 / tc[1];
        /*species 3: ch3 */
        species[3] =
            +2.84405200e+00
            +3.06898700e-03 * tc[1]
            -7.43448333e-07 * tc[2]
            +9.46290250e-11 * tc[3]
            -4.90431800e-15 * tc[4]
            +1.64378100e+04 / tc[1];
        /*species 5: ch */
        species[5] =
            +2.19622300e+00
            +1.17019050e-03 * tc[1]
            -2.35273367e-07 * tc[2]
            +2.25189550e-11 * tc[3]
            -7.71008000e-16 * tc[4]
            +7.08672300e+04 / tc[1];
        /*species 6: ch2o */
        species[6] =
            +2.99560600e+00
            +3.34066050e-03 * tc[1]
            -8.76318333e-07 * tc[2]
            +1.18428825e-10 * tc[3]
            -6.42503400e-15 * tc[4]
            -1.53203700e+04 / tc[1];
        /*species 7: hco */
        species[7] =
            +3.55727100e+00
            +1.67278650e-03 * tc[1]
            -4.45002000e-07 * tc[2]
            +6.17643250e-11 * tc[3]
            -3.42770200e-15 * tc[4]
            +3.91632400e+03 / tc[1];
        /*species 8: co2 */
        species[8] =
            +4.45362300e+00
            +1.57008450e-03 * tc[1]
            -4.26137000e-07 * tc[2]
            +5.98499250e-11 * tc[3]
            -3.33806600e-15 * tc[4]
            -4.89669600e+04 / tc[1];
        /*species 9: co */
        species[9] =
            +3.02507800e+00
            +7.21344500e-04 * tc[1]
            -1.87694267e-07 * tc[2]
            +2.54645250e-11 * tc[3]
            -1.38219040e-15 * tc[4]
            -1.42683500e+04 / tc[1];
        /*species 10: o2 */
        species[10] =
            +3.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 / tc[1];
        /*species 11: o */
        species[11] =
            +2.54206000e+00
            -1.37753100e-05 * tc[1]
            -1.03426767e-09 * tc[2]
            +1.13776675e-12 * tc[3]
            -8.73610400e-17 * tc[4]
            +2.92308000e+04 / tc[1];
        /*species 12: oh */
        species[12] =
            +2.88273000e+00
            +5.06987000e-04 * tc[1]
            -7.58959000e-08 * tc[2]
            +5.43671000e-12 * tc[3]
            -1.02526100e-16 * tc[4]
            +3.88688800e+03 / tc[1];
        /*species 14: h2o2 */
        species[14] =
            +4.57316700e+00
            +2.16806800e-03 * tc[1]
            -4.91563000e-07 * tc[2]
            +5.87226000e-11 * tc[3]
            -2.86330800e-15 * tc[4]
            -1.80069600e+04 / tc[1];
        /*species 15: h2o */
        species[15] =
            +2.67214600e+00
            +1.52814650e-03 * tc[1]
            -2.91008667e-07 * tc[2]
            +3.00249000e-11 * tc[3]
            -1.27832360e-15 * tc[4]
            -2.98992100e+04 / tc[1];
        /*species 16: c2h */
        species[16] =
            +3.98636700e+00
            +1.57156150e-03 * tc[1]
            -4.22414333e-07 * tc[2]
            +7.31090750e-11 * tc[3]
            -5.43264000e-15 * tc[4]
            +6.65588400e+04 / tc[1];
        /*species 17: hcco */
        species[17] =
            +6.75807300e+00
            +1.00020000e-03 * tc[1]
            -6.75869000e-08 * tc[2]
            -2.60283000e-11 * tc[3]
            +3.93033000e-15 * tc[4]
            +1.90151300e+04 / tc[1];
        /*species 18: c2h2 */
        species[18] =
            +4.43677000e+00
            +2.68801950e-03 * tc[1]
            -6.37605667e-07 * tc[2]
            +8.21594750e-11 * tc[3]
            -4.31342000e-15 * tc[4]
            +2.56676600e+04 / tc[1];
        /*species 20: c2h4 */
        species[20] =
            +3.52841900e+00
            +5.74259000e-03 * tc[1]
            -1.47279500e-06 * tc[2]
            +1.96115025e-10 * tc[3]
            -1.05336960e-14 * tc[4]
            +4.42828900e+03 / tc[1];
        /*species 22: c2h6 */
        species[22] =
            +4.82593800e+00
            +6.92021500e-03 * tc[1]
            -1.51908633e-06 * tc[2]
            +1.68124175e-10 * tc[3]
            -7.19632200e-15 * tc[4]
            -1.27177900e+04 / tc[1];
        /*species 24: ch3o */
        species[24] =
            +3.77080000e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832500e+02 / tc[1];
        /*species 25: hccoh */
        species[25] =
            +7.32832400e+00
            +1.66820800e-03 * tc[1]
            -1.00823500e-07 * tc[2]
            -4.45276500e-11 * tc[3]
            +6.49033600e-15 * tc[4]
            +7.59825800e+03 / tc[1];
        /*species 26: h2ccch */
        species[26] =
            +8.83104700e+00
            +2.17859750e-03 * tc[1]
            -1.36968900e-07 * tc[2]
            -5.92180750e-11 * tc[3]
            +8.75304000e-15 * tc[4]
            +3.84742000e+04 / tc[1];
        /*species 27: c3h2 */
        species[27] =
            +7.67098100e+00
            +1.37437450e-03 * tc[1]
            -1.45698100e-07 * tc[2]
            -1.61389975e-11 * tc[3]
            +3.32777400e-15 * tc[4]
            +6.25972200e+04 / tc[1];
        /*species 29: ch2co */
        species[29] =
            +6.03881700e+00
            +2.90242000e-03 * tc[1]
            -6.40318000e-07 * tc[2]
            +6.98621250e-11 * tc[3]
            -2.91773600e-15 * tc[4]
            -8.58340200e+03 / tc[1];
        /*species 30: ch2hco */
        species[30] =
            +5.97567000e+00
            +4.06529550e-03 * tc[1]
            -9.14541333e-07 * tc[2]
            +1.01757600e-10 * tc[3]
            -4.35203400e-15 * tc[4]
            +4.90321800e+02 / tc[1];
        /*species 31: ch3co */
        species[31] =
            +5.61227900e+00
            +4.22494300e-03 * tc[1]
            -9.51382333e-07 * tc[2]
            +1.05959400e-10 * tc[3]
            -4.53680800e-15 * tc[4]
            -5.18786300e+03 / tc[1];
        /*species 32: ch3hco */
        species[32] =
            +5.86865000e+00
            +5.39712000e-03 * tc[1]
            -1.21517667e-06 * tc[2]
            +1.35322800e-10 * tc[3]
            -5.79368800e-15 * tc[4]
            -2.26456900e+04 / tc[1];
        /*species 33: c2h5oh */
        species[33] =
            +6.56243650e+00
            +7.60211100e-03 * tc[1]
            -1.79655983e-06 * tc[2]
            +2.15562527e-10 * tc[3]
            -1.02579574e-14 * tc[4]
            -3.15256210e+04 / tc[1];
        /*species 37: n2 */
        species[37] =
            +2.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 / tc[1];
    }

    /*species with midpoint at T=1451 kelvin */
    if (T < 1451) {
        /*species 13: ho2 */
        species[13] =
            +3.47629499e+00
            +1.10233872e-03 * tc[1]
            +5.22806253e-07 * tc[2]
            -5.31889505e-10 * tc[3]
            +1.16628144e-13 * tc[4]
            +6.17073446e+02 / tc[1];
    } else {
        /*species 13: ho2 */
        species[13] =
            +4.58311499e+00
            +8.63651695e-04 * tc[1]
            -2.06394379e-07 * tc[2]
            +2.47984798e-11 * tc[3]
            -1.17644313e-15 * tc[4]
            +3.13365859e+01 / tc[1];
    }

    /*species with midpoint at T=1359 kelvin */
    if (T < 1359) {
        /*species 28: ch2s */
        species[28] =
            +3.32340701e+00
            +1.14031289e-03 * tc[1]
            -8.27186720e-08 * tc[2]
            -2.55163465e-11 * tc[3]
            +4.97653238e-15 * tc[4]
            +5.04954469e+04 / tc[1];
    } else {
        /*species 28: ch2s */
        species[28] =
            +3.09627180e+00
            +1.40419592e-03 * tc[1]
            -2.37871509e-07 * tc[2]
            +2.10534074e-11 * tc[3]
            -7.69281656e-16 * tc[4]
            +5.05739697e+04 / tc[1];
    }

    /*species with midpoint at T=1524 kelvin */
    if (T < 1524) {
        /*species 35: ch3choh */
        species[35] =
            +2.10026806e+00
            +9.63133395e-03 * tc[1]
            -1.72479145e-06 * tc[2]
            -4.90476455e-10 * tc[3]
            +1.87279412e-13 * tc[4]
            -6.12873036e+03 / tc[1];
    } else {
        /*species 35: ch3choh */
        species[35] =
            +7.95432747e+00
            +5.25040850e-03 * tc[1]
            -1.16877197e-06 * tc[2]
            +1.34043652e-10 * tc[3]
            -6.15956918e-15 * tc[4]
            -8.62199768e+03 / tc[1];
    }

    /*species with midpoint at T=1405 kelvin */
    if (T < 1405) {
        /*species 34: c2h4oh */
        species[34] =
            +7.76157606e-01
            +1.45501257e-02 * tc[1]
            -7.13268563e-06 * tc[2]
            +2.04208474e-09 * tc[3]
            -2.51183200e-13 * tc[4]
            -5.46808396e+03 / tc[1];
        /*species 36: ch3ch2o */
        species[36] =
            -2.71296378e-01
            +1.49419906e-02 * tc[1]
            -6.56968493e-06 * tc[2]
            +1.59334973e-09 * tc[3]
            -1.55593011e-13 * tc[4]
            -3.16397196e+03 / tc[1];
    } else {
        /*species 34: c2h4oh */
        species[34] =
            +8.75544496e+00
            +4.74726573e-03 * tc[1]
            -1.02934525e-06 * tc[2]
            +1.15884405e-10 * tc[3]
            -5.25460088e-15 * tc[4]
            -8.09602047e+03 / tc[1];
        /*species 36: ch3ch2o */
        species[36] =
            +8.31182392e+00
            +5.17131595e-03 * tc[1]
            -1.13062030e-06 * tc[2]
            +1.28053154e-10 * tc[3]
            -5.83203426e-15 * tc[4]
            -6.13097954e+03 / tc[1];
    }

    /*species with midpoint at T=1375 kelvin */
    if (T < 1375) {
        /*species 21: c2h5 */
        species[21] =
            +1.47374852e+00
            +8.18033055e-03 * tc[1]
            -1.44298942e-06 * tc[2]
            -2.93176662e-10 * tc[3]
            +1.14465177e-13 * tc[4]
            +1.33326992e+04 / tc[1];
    } else {
        /*species 21: c2h5 */
        species[21] =
            +5.60116091e+00
            +5.34885405e-03 * tc[1]
            -1.21168260e-06 * tc[2]
            +1.40454163e-10 * tc[3]
            -6.49829428e-15 * tc[4]
            +1.14539845e+04 / tc[1];
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEntropy(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1410 kelvin */
    if (T < 1410) {
        /*species 23: ch2oh */
        species[23] =
            +2.60067849e+00 * tc[0]
            +1.28447854e-02 * tc[1]
            -4.16898093e-06 * tc[2]
            +9.19092020e-10 * tc[3]
            -8.92602765e-14 * tc[4]
            +1.13272307e+01 ;
    } else {
        /*species 23: ch2oh */
        species[23] =
            +6.00127803e+00 * tc[0]
            +4.98721568e-03 * tc[1]
            -8.04767575e-07 * tc[2]
            +8.00757117e-11 * tc[3]
            -3.38956750e-15 * tc[4]
            -6.92836844e+00 ;
    }

    /*species with midpoint at T=1350 kelvin */
    if (T < 1350) {
        /*species 4: ch2 */
        species[4] =
            +3.49623989e+00 * tc[0]
            +2.20258645e-03 * tc[1]
            -2.11914399e-07 * tc[2]
            +4.54949020e-12 * tc[3]
            +1.07984343e-15 * tc[4]
            +2.94964354e+00 ;
    } else {
        /*species 4: ch2 */
        species[4] =
            +3.44310466e+00 * tc[0]
            +2.32638011e-03 * tc[1]
            -2.66573999e-07 * tc[2]
            +1.91347989e-11 * tc[3]
            -6.06737122e-16 * tc[4]
            +3.23484714e+00 ;
    }

    /*species with midpoint at T=1671 kelvin */
    if (T < 1671) {
        /*species 19: c2h3 */
        species[19] =
            +2.73925942e+00 * tc[0]
            +7.03011591e-03 * tc[1]
            +1.18323650e-06 * tc[2]
            -1.19856579e-09 * tc[3]
            +2.22939123e-13 * tc[4]
            +1.01531535e+01 ;
    } else {
        /*species 19: c2h3 */
        species[19] =
            +3.96047713e+00 * tc[0]
            +7.99426013e-03 * tc[1]
            -1.42804043e-06 * tc[2]
            +1.52783604e-10 * tc[3]
            -6.81425350e-15 * tc[4]
            +2.25663414e+00 ;
    }

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: h2 */
        species[0] =
            +3.29812400e+00 * tc[0]
            +8.24944200e-04 * tc[1]
            -4.07150750e-07 * tc[2]
            -3.15847800e-11 * tc[3]
            +1.03371800e-13 * tc[4]
            -3.29409400e+00 ;
        /*species 1: h */
        species[1] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117600e-01 ;
        /*species 2: ch4 */
        species[2] =
            +7.78741500e-01 * tc[0]
            +1.74766800e-02 * tc[1]
            -1.39170450e-05 * tc[2]
            +1.01656933e-08 * tc[3]
            -3.05982750e-12 * tc[4]
            +1.37221900e+01 ;
        /*species 3: ch3 */
        species[3] =
            +2.43044300e+00 * tc[0]
            +1.11241000e-02 * tc[1]
            -8.40110000e-06 * tc[2]
            +5.40609667e-09 * tc[3]
            -1.46623825e-12 * tc[4]
            +6.78979400e+00 ;
        /*species 5: ch */
        species[5] =
            +3.20020200e+00 * tc[0]
            +2.07287600e-03 * tc[1]
            -2.56721550e-06 * tc[2]
            +1.91129667e-09 * tc[3]
            -4.88883250e-13 * tc[4]
            +3.33158800e+00 ;
        /*species 6: ch2o */
        species[6] =
            +1.65273100e+00 * tc[0]
            +1.26314400e-02 * tc[1]
            -9.44084000e-06 * tc[2]
            +6.83343667e-09 * tc[3]
            -2.10330925e-12 * tc[4]
            +1.37848200e+01 ;
        /*species 7: hco */
        species[7] =
            +2.89833000e+00 * tc[0]
            +6.19914700e-03 * tc[1]
            -4.81154200e-06 * tc[2]
            +3.63275000e-09 * tc[3]
            -1.14372125e-12 * tc[4]
            +8.98361400e+00 ;
        /*species 8: co2 */
        species[8] =
            +2.27572500e+00 * tc[0]
            +9.92207200e-03 * tc[1]
            -5.20455500e-06 * tc[2]
            +2.28889567e-09 * tc[3]
            -5.29320000e-13 * tc[4]
            +1.01884900e+01 ;
        /*species 9: co */
        species[9] =
            +3.26245200e+00 * tc[0]
            +1.51194100e-03 * tc[1]
            -1.94087750e-06 * tc[2]
            +1.86064800e-09 * tc[3]
            -6.18737750e-13 * tc[4]
            +4.84889700e+00 ;
        /*species 10: o2 */
        species[10] =
            +3.21293600e+00 * tc[0]
            +1.12748600e-03 * tc[1]
            -2.87807500e-07 * tc[2]
            +4.37959000e-10 * tc[3]
            -2.19213850e-13 * tc[4]
            +6.03473800e+00 ;
        /*species 11: o */
        species[11] =
            +2.94642900e+00 * tc[0]
            -1.63816600e-03 * tc[1]
            +1.21051600e-06 * tc[2]
            -5.34281000e-10 * tc[3]
            +9.72674000e-14 * tc[4]
            +2.96399500e+00 ;
        /*species 12: oh */
        species[12] =
            +3.63726600e+00 * tc[0]
            +1.85091000e-04 * tc[1]
            -8.38082500e-07 * tc[2]
            +7.95734333e-10 * tc[3]
            -2.10786050e-13 * tc[4]
            +1.35886000e+00 ;
        /*species 14: h2o2 */
        species[14] =
            +3.38875400e+00 * tc[0]
            +6.56922600e-03 * tc[1]
            -7.42506500e-08 * tc[2]
            -1.54193533e-09 * tc[3]
            +6.17878750e-13 * tc[4]
            +6.78536300e+00 ;
        /*species 15: h2o */
        species[15] =
            +3.38684200e+00 * tc[0]
            +3.47498200e-03 * tc[1]
            -3.17734800e-06 * tc[2]
            +2.32286033e-09 * tc[3]
            -6.26647000e-13 * tc[4]
            +2.59023300e+00 ;
        /*species 16: c2h */
        species[16] =
            +2.73770400e+00 * tc[0]
            +8.04844600e-03 * tc[1]
            -4.62215500e-06 * tc[2]
            +2.17508633e-09 * tc[3]
            -4.84895000e-13 * tc[4]
            +7.30022000e+00 ;
        /*species 17: hcco */
        species[17] =
            +5.04796500e+00 * tc[0]
            +4.45347800e-03 * tc[1]
            +1.13414150e-07 * tc[2]
            -4.94031667e-10 * tc[3]
            +5.62685500e-14 * tc[4]
            +4.81843900e-01 ;
        /*species 18: c2h2 */
        species[18] =
            +2.01356200e+00 * tc[0]
            +1.51904500e-02 * tc[1]
            -8.08159500e-06 * tc[2]
            +3.02633067e-09 * tc[3]
            -4.78186500e-13 * tc[4]
            +8.80537800e+00 ;
        /*species 20: c2h4 */
        species[20] =
            -8.61488000e-01 * tc[0]
            +2.79616300e-02 * tc[1]
            -1.69433850e-05 * tc[2]
            +9.28384000e-09 * tc[3]
            -2.43446975e-12 * tc[4]
            +2.42114900e+01 ;
        /*species 22: c2h6 */
        species[22] =
            +1.46253900e+00 * tc[0]
            +1.54946700e-02 * tc[1]
            +2.89025350e-06 * tc[2]
            -4.19277333e-09 * tc[3]
            +1.14656675e-12 * tc[4]
            +1.44322900e+01 ;
        /*species 24: ch3o */
        species[24] =
            +2.10620400e+00 * tc[0]
            +7.21659500e-03 * tc[1]
            +2.66923600e-06 * tc[2]
            -2.45921200e-09 * tc[3]
            +5.18902750e-13 * tc[4]
            +1.31521800e+01 ;
        /*species 25: hccoh */
        species[25] =
            +3.89946500e+00 * tc[0]
            +9.70107500e-03 * tc[1]
            -1.55965450e-07 * tc[2]
            -1.84591067e-09 * tc[3]
            +6.16433000e-13 * tc[4]
            +4.49187500e+00 ;
        /*species 26: h2ccch */
        species[26] =
            +4.75420000e+00 * tc[0]
            +1.10802800e-02 * tc[1]
            +1.39666150e-07 * tc[2]
            -1.82640400e-09 * tc[3]
            +4.87407250e-13 * tc[4]
            +5.85454900e-01 ;
        /*species 27: c3h2 */
        species[27] =
            +3.16671400e+00 * tc[0]
            +2.48257200e-02 * tc[1]
            -2.29581850e-05 * tc[2]
            +1.42267300e-08 * tc[3]
            -3.70538000e-12 * tc[4]
            +8.86944600e+00 ;
        /*species 29: ch2co */
        species[29] =
            +2.97497100e+00 * tc[0]
            +1.21187100e-02 * tc[1]
            -1.17252300e-06 * tc[2]
            -2.15556167e-09 * tc[3]
            +9.76412250e-13 * tc[4]
            +8.67355300e+00 ;
        /*species 30: ch2hco */
        species[30] =
            +3.40906200e+00 * tc[0]
            +1.07385700e-02 * tc[1]
            +9.45746000e-07 * tc[2]
            -2.38619433e-09 * tc[3]
            +7.16846250e-13 * tc[4]
            +9.55829000e+00 ;
        /*species 31: ch3co */
        species[31] =
            +3.12527800e+00 * tc[0]
            +9.77822000e-03 * tc[1]
            +2.26072400e-06 * tc[2]
            -3.00315400e-09 * tc[3]
            +7.98429500e-13 * tc[4]
            +1.12288500e+01 ;
        /*species 32: ch3hco */
        species[32] =
            +2.50569500e+00 * tc[0]
            +1.33699100e-02 * tc[1]
            +2.33597650e-06 * tc[2]
            -3.76046667e-09 * tc[3]
            +1.06589150e-12 * tc[4]
            +1.33508900e+01 ;
        /*species 33: c2h5oh */
        species[33] =
            +4.85869570e+00 * tc[0]
            -3.74017260e-03 * tc[1]
            +3.47776890e-05 * tc[2]
            -2.95515987e-08 * tc[3]
            +8.79220875e-12 * tc[4]
            +4.80185450e+00 ;
        /*species 37: n2 */
        species[37] =
            +3.29867700e+00 * tc[0]
            +1.40824000e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213750e-13 * tc[4]
            +3.95037200e+00 ;
    } else {
        /*species 0: h2 */
        species[0] =
            +2.99142300e+00 * tc[0]
            +7.00064400e-04 * tc[1]
            -2.81691450e-08 * tc[2]
            -3.07719267e-12 * tc[3]
            +3.95688000e-16 * tc[4]
            -1.35511000e+00 ;
        /*species 1: h */
        species[1] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117600e-01 ;
        /*species 2: ch4 */
        species[2] =
            +1.68347900e+00 * tc[0]
            +1.02372400e-02 * tc[1]
            -1.93756450e-06 * tc[2]
            +2.26186167e-10 * tc[3]
            -1.12585575e-14 * tc[4]
            +9.62339500e+00 ;
        /*species 3: ch3 */
        species[3] =
            +2.84405200e+00 * tc[0]
            +6.13797400e-03 * tc[1]
            -1.11517250e-06 * tc[2]
            +1.26172033e-10 * tc[3]
            -6.13039750e-15 * tc[4]
            +5.45269700e+00 ;
        /*species 5: ch */
        species[5] =
            +2.19622300e+00 * tc[0]
            +2.34038100e-03 * tc[1]
            -3.52910050e-07 * tc[2]
            +3.00252733e-11 * tc[3]
            -9.63760000e-16 * tc[4]
            +9.17837300e+00 ;
        /*species 6: ch2o */
        species[6] =
            +2.99560600e+00 * tc[0]
            +6.68132100e-03 * tc[1]
            -1.31447750e-06 * tc[2]
            +1.57905100e-10 * tc[3]
            -8.03129250e-15 * tc[4]
            +6.91257200e+00 ;
        /*species 7: hco */
        species[7] =
            +3.55727100e+00 * tc[0]
            +3.34557300e-03 * tc[1]
            -6.67503000e-07 * tc[2]
            +8.23524333e-11 * tc[3]
            -4.28462750e-15 * tc[4]
            +5.55229900e+00 ;
        /*species 8: co2 */
        species[8] =
            +4.45362300e+00 * tc[0]
            +3.14016900e-03 * tc[1]
            -6.39205500e-07 * tc[2]
            +7.97999000e-11 * tc[3]
            -4.17258250e-15 * tc[4]
            -9.55395900e-01 ;
        /*species 9: co */
        species[9] =
            +3.02507800e+00 * tc[0]
            +1.44268900e-03 * tc[1]
            -2.81541400e-07 * tc[2]
            +3.39527000e-11 * tc[3]
            -1.72773800e-15 * tc[4]
            +6.10821800e+00 ;
        /*species 10: o2 */
        species[10] =
            +3.69757800e+00 * tc[0]
            +6.13519700e-04 * tc[1]
            -6.29421000e-08 * tc[2]
            +5.91760333e-12 * tc[3]
            -2.84108750e-16 * tc[4]
            +3.18916600e+00 ;
        /*species 11: o */
        species[11] =
            +2.54206000e+00 * tc[0]
            -2.75506200e-05 * tc[1]
            -1.55140150e-09 * tc[2]
            +1.51702233e-12 * tc[3]
            -1.09201300e-16 * tc[4]
            +4.92030800e+00 ;
        /*species 12: oh */
        species[12] =
            +2.88273000e+00 * tc[0]
            +1.01397400e-03 * tc[1]
            -1.13843850e-07 * tc[2]
            +7.24894667e-12 * tc[3]
            -1.28157625e-16 * tc[4]
            +5.59571200e+00 ;
        /*species 14: h2o2 */
        species[14] =
            +4.57316700e+00 * tc[0]
            +4.33613600e-03 * tc[1]
            -7.37344500e-07 * tc[2]
            +7.82968000e-11 * tc[3]
            -3.57913500e-15 * tc[4]
            +5.01137000e-01 ;
        /*species 15: h2o */
        species[15] =
            +2.67214600e+00 * tc[0]
            +3.05629300e-03 * tc[1]
            -4.36513000e-07 * tc[2]
            +4.00332000e-11 * tc[3]
            -1.59790450e-15 * tc[4]
            +6.86281700e+00 ;
        /*species 16: c2h */
        species[16] =
            +3.98636700e+00 * tc[0]
            +3.14312300e-03 * tc[1]
            -6.33621500e-07 * tc[2]
            +9.74787667e-11 * tc[3]
            -6.79080000e-15 * tc[4]
            +1.19106300e+00 ;
        /*species 17: hcco */
        species[17] =
            +6.75807300e+00 * tc[0]
            +2.00040000e-03 * tc[1]
            -1.01380350e-07 * tc[2]
            -3.47044000e-11 * tc[3]
            +4.91291250e-15 * tc[4]
            -9.07126200e+00 ;
        /*species 18: c2h2 */
        species[18] =
            +4.43677000e+00 * tc[0]
            +5.37603900e-03 * tc[1]
            -9.56408500e-07 * tc[2]
            +1.09545967e-10 * tc[3]
            -5.39177500e-15 * tc[4]
            -2.80033800e+00 ;
        /*species 20: c2h4 */
        species[20] =
            +3.52841900e+00 * tc[0]
            +1.14851800e-02 * tc[1]
            -2.20919250e-06 * tc[2]
            +2.61486700e-10 * tc[3]
            -1.31671200e-14 * tc[4]
            +2.23038900e+00 ;
        /*species 22: c2h6 */
        species[22] =
            +4.82593800e+00 * tc[0]
            +1.38404300e-02 * tc[1]
            -2.27862950e-06 * tc[2]
            +2.24165567e-10 * tc[3]
            -8.99540250e-15 * tc[4]
            -5.23950700e+00 ;
        /*species 24: ch3o */
        species[24] =
            +3.77080000e+00 * tc[0]
            +7.87149700e-03 * tc[1]
            -1.32819200e-06 * tc[2]
            +1.31481033e-10 * tc[3]
            -5.28154000e-15 * tc[4]
            +2.92957500e+00 ;
        /*species 25: hccoh */
        species[25] =
            +7.32832400e+00 * tc[0]
            +3.33641600e-03 * tc[1]
            -1.51235250e-07 * tc[2]
            -5.93702000e-11 * tc[3]
            +8.11292000e-15 * tc[4]
            -1.40121400e+01 ;
        /*species 26: h2ccch */
        species[26] =
            +8.83104700e+00 * tc[0]
            +4.35719500e-03 * tc[1]
            -2.05453350e-07 * tc[2]
            -7.89574333e-11 * tc[3]
            +1.09413000e-14 * tc[4]
            -2.17791900e+01 ;
        /*species 27: c3h2 */
        species[27] =
            +7.67098100e+00 * tc[0]
            +2.74874900e-03 * tc[1]
            -2.18547150e-07 * tc[2]
            -2.15186633e-11 * tc[3]
            +4.15971750e-15 * tc[4]
            -1.23689000e+01 ;
        /*species 29: ch2co */
        species[29] =
            +6.03881700e+00 * tc[0]
            +5.80484000e-03 * tc[1]
            -9.60477000e-07 * tc[2]
            +9.31495000e-11 * tc[3]
            -3.64717000e-15 * tc[4]
            -7.65758100e+00 ;
        /*species 30: ch2hco */
        species[30] =
            +5.97567000e+00 * tc[0]
            +8.13059100e-03 * tc[1]
            -1.37181200e-06 * tc[2]
            +1.35676800e-10 * tc[3]
            -5.44004250e-15 * tc[4]
            -5.04525100e+00 ;
        /*species 31: ch3co */
        species[31] =
            +5.61227900e+00 * tc[0]
            +8.44988600e-03 * tc[1]
            -1.42707350e-06 * tc[2]
            +1.41279200e-10 * tc[3]
            -5.67101000e-15 * tc[4]
            -3.27494900e+00 ;
        /*species 32: ch3hco */
        species[32] =
            +5.86865000e+00 * tc[0]
            +1.07942400e-02 * tc[1]
            -1.82276500e-06 * tc[2]
            +1.80430400e-10 * tc[3]
            -7.24211000e-15 * tc[4]
            -6.01294600e+00 ;
        /*species 33: c2h5oh */
        species[33] =
            +6.56243650e+00 * tc[0]
            +1.52042220e-02 * tc[1]
            -2.69483975e-06 * tc[2]
            +2.87416703e-10 * tc[3]
            -1.28224468e-14 * tc[4]
            -9.47302020e+00 ;
        /*species 37: n2 */
        species[37] =
            +2.92664000e+00 * tc[0]
            +1.48797700e-03 * tc[1]
            -2.84238050e-07 * tc[2]
            +3.36568000e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
    }

    /*species with midpoint at T=1451 kelvin */
    if (T < 1451) {
        /*species 13: ho2 */
        species[13] =
            +3.47629499e+00 * tc[0]
            +2.20467745e-03 * tc[1]
            +7.84209380e-07 * tc[2]
            -7.09186007e-10 * tc[3]
            +1.45785180e-13 * tc[4]
            +7.02308516e+00 ;
    } else {
        /*species 13: ho2 */
        species[13] =
            +4.58311499e+00 * tc[0]
            +1.72730339e-03 * tc[1]
            -3.09591568e-07 * tc[2]
            +3.30646398e-11 * tc[3]
            -1.47055392e-15 * tc[4]
            +3.46199386e-01 ;
    }

    /*species with midpoint at T=1359 kelvin */
    if (T < 1359) {
        /*species 28: ch2s */
        species[28] =
            +3.32340701e+00 * tc[0]
            +2.28062579e-03 * tc[1]
            -1.24078008e-07 * tc[2]
            -3.40217953e-11 * tc[3]
            +6.22066548e-15 * tc[4]
            +3.09401497e+00 ;
    } else {
        /*species 28: ch2s */
        species[28] =
            +3.09627180e+00 * tc[0]
            +2.80839185e-03 * tc[1]
            -3.56807263e-07 * tc[2]
            +2.80712099e-11 * tc[3]
            -9.61602070e-16 * tc[4]
            +4.31392269e+00 ;
    }

    /*species with midpoint at T=1524 kelvin */
    if (T < 1524) {
        /*species 35: ch3choh */
        species[35] =
            +2.10026806e+00 * tc[0]
            +1.92626679e-02 * tc[1]
            -2.58718718e-06 * tc[2]
            -6.53968607e-10 * tc[3]
            +2.34099265e-13 * tc[4]
            +1.40909302e+01 ;
    } else {
        /*species 35: ch3choh */
        species[35] =
            +7.95432747e+00 * tc[0]
            +1.05008170e-02 * tc[1]
            -1.75315795e-06 * tc[2]
            +1.78724869e-10 * tc[3]
            -7.69946147e-15 * tc[4]
            -1.90411129e+01 ;
    }

    /*species with midpoint at T=1405 kelvin */
    if (T < 1405) {
        /*species 34: c2h4oh */
        species[34] =
            +7.76157606e-01 * tc[0]
            +2.91002514e-02 * tc[1]
            -1.06990284e-05 * tc[2]
            +2.72277965e-09 * tc[3]
            -3.13979000e-13 * tc[4]
            +2.21060185e+01 ;
        /*species 36: ch3ch2o */
        species[36] =
            -2.71296378e-01 * tc[0]
            +2.98839812e-02 * tc[1]
            -9.85452740e-06 * tc[2]
            +2.12446631e-09 * tc[3]
            -1.94491264e-13 * tc[4]
            +2.47706003e+01 ;
    } else {
        /*species 34: c2h4oh */
        species[34] =
            +8.75544496e+00 * tc[0]
            +9.49453147e-03 * tc[1]
            -1.54401788e-06 * tc[2]
            +1.54512540e-10 * tc[3]
            -6.56825110e-15 * tc[4]
            -2.03271929e+01 ;
        /*species 36: ch3ch2o */
        species[36] =
            +8.31182392e+00 * tc[0]
            +1.03426319e-02 * tc[1]
            -1.69593044e-06 * tc[2]
            +1.70737539e-10 * tc[3]
            -7.29004283e-15 * tc[4]
            -2.13985581e+01 ;
    }

    /*species with midpoint at T=1375 kelvin */
    if (T < 1375) {
        /*species 21: c2h5 */
        species[21] =
            +1.47374852e+00 * tc[0]
            +1.63606611e-02 * tc[1]
            -2.16448413e-06 * tc[2]
            -3.90902217e-10 * tc[3]
            +1.43081471e-13 * tc[4]
            +1.66349852e+01 ;
    } else {
        /*species 21: c2h5 */
        species[21] =
            +5.60116091e+00 * tc[0]
            +1.06977081e-02 * tc[1]
            -1.81752390e-06 * tc[2]
            +1.87272218e-10 * tc[3]
            -8.12286785e-15 * tc[4]
            -7.02252408e+00 ;
    }
    return;
}


/*save molecular weights into array */
void molecularWeight(double * wt)
{
    wt[0] = 2.015940; /*h2 */
    wt[1] = 1.007970; /*h */
    wt[2] = 16.043030; /*ch4 */
    wt[3] = 15.035060; /*ch3 */
    wt[4] = 14.027090; /*ch2 */
    wt[5] = 13.019120; /*ch */
    wt[6] = 30.026490; /*ch2o */
    wt[7] = 29.018520; /*hco */
    wt[8] = 44.009950; /*co2 */
    wt[9] = 28.010550; /*co */
    wt[10] = 31.998800; /*o2 */
    wt[11] = 15.999400; /*o */
    wt[12] = 17.007370; /*oh */
    wt[13] = 33.006770; /*ho2 */
    wt[14] = 34.014740; /*h2o2 */
    wt[15] = 18.015340; /*h2o */
    wt[16] = 25.030270; /*c2h */
    wt[17] = 41.029670; /*hcco */
    wt[18] = 26.038240; /*c2h2 */
    wt[19] = 27.046210; /*c2h3 */
    wt[20] = 28.054180; /*c2h4 */
    wt[21] = 29.062150; /*c2h5 */
    wt[22] = 30.070120; /*c2h6 */
    wt[23] = 31.034460; /*ch2oh */
    wt[24] = 31.034460; /*ch3o */
    wt[25] = 42.037640; /*hccoh */
    wt[26] = 39.057360; /*h2ccch */
    wt[27] = 38.049390; /*c3h2 */
    wt[28] = 14.027090; /*ch2s */
    wt[29] = 42.037640; /*ch2co */
    wt[30] = 43.045610; /*ch2hco */
    wt[31] = 43.045610; /*ch3co */
    wt[32] = 44.053580; /*ch3hco */
    wt[33] = 46.069520; /*c2h5oh */
    wt[34] = 45.061550; /*c2h4oh */
    wt[35] = 45.061550; /*ch3choh */
    wt[36] = 45.061550; /*ch3ch2o */
    wt[37] = 28.013400; /*n2 */

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
    y[0] = phi[0]*2.015940;   XW += y[0]; /*h2 */
    y[1] = phi[1]*1.007970;   XW += y[1]; /*h */
    y[2] = phi[2]*16.043030;   XW += y[2]; /*ch4 */
    y[3] = phi[3]*15.035060;   XW += y[3]; /*ch3 */
    y[4] = phi[4]*14.027090;   XW += y[4]; /*ch2 */
    y[5] = phi[5]*13.019120;   XW += y[5]; /*ch */
    y[6] = phi[6]*30.026490;   XW += y[6]; /*ch2o */
    y[7] = phi[7]*29.018520;   XW += y[7]; /*hco */
    y[8] = phi[8]*44.009950;   XW += y[8]; /*co2 */
    y[9] = phi[9]*28.010550;   XW += y[9]; /*co */
    y[10] = phi[10]*31.998800;   XW += y[10]; /*o2 */
    y[11] = phi[11]*15.999400;   XW += y[11]; /*o */
    y[12] = phi[12]*17.007370;   XW += y[12]; /*oh */
    y[13] = phi[13]*33.006770;   XW += y[13]; /*ho2 */
    y[14] = phi[14]*34.014740;   XW += y[14]; /*h2o2 */
    y[15] = phi[15]*18.015340;   XW += y[15]; /*h2o */
    y[16] = phi[16]*25.030270;   XW += y[16]; /*c2h */
    y[17] = phi[17]*41.029670;   XW += y[17]; /*hcco */
    y[18] = phi[18]*26.038240;   XW += y[18]; /*c2h2 */
    y[19] = phi[19]*27.046210;   XW += y[19]; /*c2h3 */
    y[20] = phi[20]*28.054180;   XW += y[20]; /*c2h4 */
    y[21] = phi[21]*29.062150;   XW += y[21]; /*c2h5 */
    y[22] = phi[22]*30.070120;   XW += y[22]; /*c2h6 */
    y[23] = phi[23]*31.034460;   XW += y[23]; /*ch2oh */
    y[24] = phi[24]*31.034460;   XW += y[24]; /*ch3o */
    y[25] = phi[25]*42.037640;   XW += y[25]; /*hccoh */
    y[26] = phi[26]*39.057360;   XW += y[26]; /*h2ccch */
    y[27] = phi[27]*38.049390;   XW += y[27]; /*c3h2 */
    y[28] = phi[28]*14.027090;   XW += y[28]; /*ch2s */
    y[29] = phi[29]*42.037640;   XW += y[29]; /*ch2co */
    y[30] = phi[30]*43.045610;   XW += y[30]; /*ch2hco */
    y[31] = phi[31]*43.045610;   XW += y[31]; /*ch3co */
    y[32] = phi[32]*44.053580;   XW += y[32]; /*ch3hco */
    y[33] = phi[33]*46.069520;   XW += y[33]; /*c2h5oh */
    y[34] = phi[34]*45.061550;   XW += y[34]; /*c2h4oh */
    y[35] = phi[35]*45.061550;   XW += y[35]; /*ch3choh */
    y[36] = phi[36]*45.061550;   XW += y[36]; /*ch3ch2o */
    y[37] = phi[37]*28.013400;   XW += y[37]; /*n2 */
    for (id = 0; id < 38; ++id) {
        y[id] = y[id]/XW;
    }

    return;
}


/*convert y[species] (mass fracs) to phi[species] (specific mole num) */
void feytphi_(double * y, int * iwrk, double * rwrk, double * phi)
{
    phi[0] = y[0]/ 2.01594000e-03; /*h2 (wt in kg) */
    phi[1] = y[1]/ 1.00797000e-03; /*h (wt in kg) */
    phi[2] = y[2]/ 1.60430300e-02; /*ch4 (wt in kg) */
    phi[3] = y[3]/ 1.50350600e-02; /*ch3 (wt in kg) */
    phi[4] = y[4]/ 1.40270900e-02; /*ch2 (wt in kg) */
    phi[5] = y[5]/ 1.30191200e-02; /*ch (wt in kg) */
    phi[6] = y[6]/ 3.00264900e-02; /*ch2o (wt in kg) */
    phi[7] = y[7]/ 2.90185200e-02; /*hco (wt in kg) */
    phi[8] = y[8]/ 4.40099500e-02; /*co2 (wt in kg) */
    phi[9] = y[9]/ 2.80105500e-02; /*co (wt in kg) */
    phi[10] = y[10]/ 3.19988000e-02; /*o2 (wt in kg) */
    phi[11] = y[11]/ 1.59994000e-02; /*o (wt in kg) */
    phi[12] = y[12]/ 1.70073700e-02; /*oh (wt in kg) */
    phi[13] = y[13]/ 3.30067700e-02; /*ho2 (wt in kg) */
    phi[14] = y[14]/ 3.40147400e-02; /*h2o2 (wt in kg) */
    phi[15] = y[15]/ 1.80153400e-02; /*h2o (wt in kg) */
    phi[16] = y[16]/ 2.50302700e-02; /*c2h (wt in kg) */
    phi[17] = y[17]/ 4.10296700e-02; /*hcco (wt in kg) */
    phi[18] = y[18]/ 2.60382400e-02; /*c2h2 (wt in kg) */
    phi[19] = y[19]/ 2.70462100e-02; /*c2h3 (wt in kg) */
    phi[20] = y[20]/ 2.80541800e-02; /*c2h4 (wt in kg) */
    phi[21] = y[21]/ 2.90621500e-02; /*c2h5 (wt in kg) */
    phi[22] = y[22]/ 3.00701200e-02; /*c2h6 (wt in kg) */
    phi[23] = y[23]/ 3.10344600e-02; /*ch2oh (wt in kg) */
    phi[24] = y[24]/ 3.10344600e-02; /*ch3o (wt in kg) */
    phi[25] = y[25]/ 4.20376400e-02; /*hccoh (wt in kg) */
    phi[26] = y[26]/ 3.90573600e-02; /*h2ccch (wt in kg) */
    phi[27] = y[27]/ 3.80493900e-02; /*c3h2 (wt in kg) */
    phi[28] = y[28]/ 1.40270900e-02; /*ch2s (wt in kg) */
    phi[29] = y[29]/ 4.20376400e-02; /*ch2co (wt in kg) */
    phi[30] = y[30]/ 4.30456100e-02; /*ch2hco (wt in kg) */
    phi[31] = y[31]/ 4.30456100e-02; /*ch3co (wt in kg) */
    phi[32] = y[32]/ 4.40535800e-02; /*ch3hco (wt in kg) */
    phi[33] = y[33]/ 4.60695200e-02; /*c2h5oh (wt in kg) */
    phi[34] = y[34]/ 4.50615500e-02; /*c2h4oh (wt in kg) */
    phi[35] = y[35]/ 4.50615500e-02; /*ch3choh (wt in kg) */
    phi[36] = y[36]/ 4.50615500e-02; /*ch3ch2o (wt in kg) */
    phi[37] = y[37]/ 2.80134000e-02; /*n2 (wt in kg) */

    return;
}


/*reverse of ytcr, useful for rate computations */
void fectyr_(double * c, double * rho, int * iwrk, double * rwrk, double * y)
{
    y[0] = c[0] * 2.015940 / (*rho); 
    y[1] = c[1] * 1.007970 / (*rho); 
    y[2] = c[2] * 16.043030 / (*rho); 
    y[3] = c[3] * 15.035060 / (*rho); 
    y[4] = c[4] * 14.027090 / (*rho); 
    y[5] = c[5] * 13.019120 / (*rho); 
    y[6] = c[6] * 30.026490 / (*rho); 
    y[7] = c[7] * 29.018520 / (*rho); 
    y[8] = c[8] * 44.009950 / (*rho); 
    y[9] = c[9] * 28.010550 / (*rho); 
    y[10] = c[10] * 31.998800 / (*rho); 
    y[11] = c[11] * 15.999400 / (*rho); 
    y[12] = c[12] * 17.007370 / (*rho); 
    y[13] = c[13] * 33.006770 / (*rho); 
    y[14] = c[14] * 34.014740 / (*rho); 
    y[15] = c[15] * 18.015340 / (*rho); 
    y[16] = c[16] * 25.030270 / (*rho); 
    y[17] = c[17] * 41.029670 / (*rho); 
    y[18] = c[18] * 26.038240 / (*rho); 
    y[19] = c[19] * 27.046210 / (*rho); 
    y[20] = c[20] * 28.054180 / (*rho); 
    y[21] = c[21] * 29.062150 / (*rho); 
    y[22] = c[22] * 30.070120 / (*rho); 
    y[23] = c[23] * 31.034460 / (*rho); 
    y[24] = c[24] * 31.034460 / (*rho); 
    y[25] = c[25] * 42.037640 / (*rho); 
    y[26] = c[26] * 39.057360 / (*rho); 
    y[27] = c[27] * 38.049390 / (*rho); 
    y[28] = c[28] * 14.027090 / (*rho); 
    y[29] = c[29] * 42.037640 / (*rho); 
    y[30] = c[30] * 43.045610 / (*rho); 
    y[31] = c[31] * 43.045610 / (*rho); 
    y[32] = c[32] * 44.053580 / (*rho); 
    y[33] = c[33] * 46.069520 / (*rho); 
    y[34] = c[34] * 45.061550 / (*rho); 
    y[35] = c[35] * 45.061550 / (*rho); 
    y[36] = c[36] * 45.061550 / (*rho); 
    y[37] = c[37] * 28.013400 / (*rho); 

    return;
}


/*ddebdf compatible right hand side of CV burner */
/*rwrk[0] and rwrk[1] should contain rho and ene respectively */
/*working variable phi contains specific mole numbers */
void fecvrhs_(double * time, double * phi, double * phidot, double * rwrk, int * iwrk)
{
    double rho,ene; /*CV Parameters */
    double y[38], wdot[38]; /*temporary storage */
    int i; /*Loop counter */
    double temperature,pressure; /*temporary var */
    rho = rwrk[0];
    ene = rwrk[1];
    fephity_(phi, iwrk, rwrk, y);
    feeytt_(&ene, y, iwrk, rwrk, &temperature);
    CKPY(&rho, &temperature,  y, iwrk, rwrk, &pressure);
    CKWYP(&pressure, &temperature,  y, iwrk, rwrk, wdot);
    for (i=0; i<38; ++i) phidot[i] = wdot[i] / (rho/1000.0); 

    return;
}


/*returns the dimensionality of the cv burner (number of species) */
int fecvdim_()
{
    return 38;
}


/*ddebdf compatible right hand side of ZND solver */
/*rwrk[0] : scaling factor for pressure */
/*rwrk[1] : preshock density (g/cc)  */
/*rwrk[2] : detonation velocity (cm/s)  */
/*solution vector: [P; rho; y0 ... ylast]  */
void fezndrhs_(double * time, double * z, double * zdot, double * rwrk, int * iwrk)
{
    double psc,rho1,udet; /*ZND Parameters */
    double wt[38], hms[38], wdot[38]; /*temporary storage */
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
    for (i=0; i<38; ++i) {
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
    return 41;
}


/*returns the name of the source mechanism file  */
char* femechfile_()
{
    return "";
}


/*returns the species number */
int fesymnum_(const char* s1)
{
    if (strcmp(s1, "h2")==0) return 0; 
    if (strcmp(s1, "h")==0) return 1; 
    if (strcmp(s1, "ch4")==0) return 2; 
    if (strcmp(s1, "ch3")==0) return 3; 
    if (strcmp(s1, "ch2")==0) return 4; 
    if (strcmp(s1, "ch")==0) return 5; 
    if (strcmp(s1, "ch2o")==0) return 6; 
    if (strcmp(s1, "hco")==0) return 7; 
    if (strcmp(s1, "co2")==0) return 8; 
    if (strcmp(s1, "co")==0) return 9; 
    if (strcmp(s1, "o2")==0) return 10; 
    if (strcmp(s1, "o")==0) return 11; 
    if (strcmp(s1, "oh")==0) return 12; 
    if (strcmp(s1, "ho2")==0) return 13; 
    if (strcmp(s1, "h2o2")==0) return 14; 
    if (strcmp(s1, "h2o")==0) return 15; 
    if (strcmp(s1, "c2h")==0) return 16; 
    if (strcmp(s1, "hcco")==0) return 17; 
    if (strcmp(s1, "c2h2")==0) return 18; 
    if (strcmp(s1, "c2h3")==0) return 19; 
    if (strcmp(s1, "c2h4")==0) return 20; 
    if (strcmp(s1, "c2h5")==0) return 21; 
    if (strcmp(s1, "c2h6")==0) return 22; 
    if (strcmp(s1, "ch2oh")==0) return 23; 
    if (strcmp(s1, "ch3o")==0) return 24; 
    if (strcmp(s1, "hccoh")==0) return 25; 
    if (strcmp(s1, "h2ccch")==0) return 26; 
    if (strcmp(s1, "c3h2")==0) return 27; 
    if (strcmp(s1, "ch2s")==0) return 28; 
    if (strcmp(s1, "ch2co")==0) return 29; 
    if (strcmp(s1, "ch2hco")==0) return 30; 
    if (strcmp(s1, "ch3co")==0) return 31; 
    if (strcmp(s1, "ch3hco")==0) return 32; 
    if (strcmp(s1, "c2h5oh")==0) return 33; 
    if (strcmp(s1, "c2h4oh")==0) return 34; 
    if (strcmp(s1, "ch3choh")==0) return 35; 
    if (strcmp(s1, "ch3ch2o")==0) return 36; 
    if (strcmp(s1, "n2")==0) return 37; 
    /*species name not found */
    return -1;
}


/*returns the species name */
char* fesymname_(int sn)
{
    if (sn==0) return "h2"; 
    if (sn==1) return "h"; 
    if (sn==2) return "ch4"; 
    if (sn==3) return "ch3"; 
    if (sn==4) return "ch2"; 
    if (sn==5) return "ch"; 
    if (sn==6) return "ch2o"; 
    if (sn==7) return "hco"; 
    if (sn==8) return "co2"; 
    if (sn==9) return "co"; 
    if (sn==10) return "o2"; 
    if (sn==11) return "o"; 
    if (sn==12) return "oh"; 
    if (sn==13) return "ho2"; 
    if (sn==14) return "h2o2"; 
    if (sn==15) return "h2o"; 
    if (sn==16) return "c2h"; 
    if (sn==17) return "hcco"; 
    if (sn==18) return "c2h2"; 
    if (sn==19) return "c2h3"; 
    if (sn==20) return "c2h4"; 
    if (sn==21) return "c2h5"; 
    if (sn==22) return "c2h6"; 
    if (sn==23) return "ch2oh"; 
    if (sn==24) return "ch3o"; 
    if (sn==25) return "hccoh"; 
    if (sn==26) return "h2ccch"; 
    if (sn==27) return "c3h2"; 
    if (sn==28) return "ch2s"; 
    if (sn==29) return "ch2co"; 
    if (sn==30) return "ch2hco"; 
    if (sn==31) return "ch3co"; 
    if (sn==32) return "ch3hco"; 
    if (sn==33) return "c2h5oh"; 
    if (sn==34) return "c2h4oh"; 
    if (sn==35) return "ch3choh"; 
    if (sn==36) return "ch3ch2o"; 
    if (sn==37) return "n2"; 
    /*species name not found */
    return "NOTFOUND";
}

/* End of file  */
