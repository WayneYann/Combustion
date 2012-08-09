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
    kname[ 4*lenkname + 2 ] = ' ';

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
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*12.011150; /*C */
    XW += x[9]*13.019120; /*CH */
    XW += x[10]*14.027090; /*CH2 */
    XW += x[11]*14.027090; /*CH2(S) */
    XW += x[12]*15.035060; /*CH3 */
    XW += x[13]*16.043030; /*CH4 */
    XW += x[14]*28.010550; /*CO */
    XW += x[15]*44.009950; /*CO2 */
    XW += x[16]*29.018520; /*HCO */
    XW += x[17]*30.026490; /*CH2O */
    XW += x[18]*31.034460; /*CH2OH */
    XW += x[19]*31.034460; /*CH3O */
    XW += x[20]*32.042430; /*CH3OH */
    XW += x[21]*25.030270; /*C2H */
    XW += x[22]*26.038240; /*C2H2 */
    XW += x[23]*27.046210; /*C2H3 */
    XW += x[24]*28.054180; /*C2H4 */
    XW += x[25]*29.062150; /*C2H5 */
    XW += x[26]*30.070120; /*C2H6 */
    XW += x[27]*41.029670; /*HCCO */
    XW += x[28]*42.037640; /*CH2CO */
    XW += x[29]*42.037640; /*HCCOH */
    XW += x[30]*28.013400; /*N2 */
    XW += x[31]*39.948000; /*AR */
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
    YOW += y[8]/12.011150; /*C */
    YOW += y[9]/13.019120; /*CH */
    YOW += y[10]/14.027090; /*CH2 */
    YOW += y[11]/14.027090; /*CH2(S) */
    YOW += y[12]/15.035060; /*CH3 */
    YOW += y[13]/16.043030; /*CH4 */
    YOW += y[14]/28.010550; /*CO */
    YOW += y[15]/44.009950; /*CO2 */
    YOW += y[16]/29.018520; /*HCO */
    YOW += y[17]/30.026490; /*CH2O */
    YOW += y[18]/31.034460; /*CH2OH */
    YOW += y[19]/31.034460; /*CH3O */
    YOW += y[20]/32.042430; /*CH3OH */
    YOW += y[21]/25.030270; /*C2H */
    YOW += y[22]/26.038240; /*C2H2 */
    YOW += y[23]/27.046210; /*C2H3 */
    YOW += y[24]/28.054180; /*C2H4 */
    YOW += y[25]/29.062150; /*C2H5 */
    YOW += y[26]/30.070120; /*C2H6 */
    YOW += y[27]/41.029670; /*HCCO */
    YOW += y[28]/42.037640; /*CH2CO */
    YOW += y[29]/42.037640; /*HCCOH */
    YOW += y[30]/28.013400; /*N2 */
    YOW += y[31]/39.948000; /*AR */
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
    W += c[8]*12.011150; /*C */
    W += c[9]*13.019120; /*CH */
    W += c[10]*14.027090; /*CH2 */
    W += c[11]*14.027090; /*CH2(S) */
    W += c[12]*15.035060; /*CH3 */
    W += c[13]*16.043030; /*CH4 */
    W += c[14]*28.010550; /*CO */
    W += c[15]*44.009950; /*CO2 */
    W += c[16]*29.018520; /*HCO */
    W += c[17]*30.026490; /*CH2O */
    W += c[18]*31.034460; /*CH2OH */
    W += c[19]*31.034460; /*CH3O */
    W += c[20]*32.042430; /*CH3OH */
    W += c[21]*25.030270; /*C2H */
    W += c[22]*26.038240; /*C2H2 */
    W += c[23]*27.046210; /*C2H3 */
    W += c[24]*28.054180; /*C2H4 */
    W += c[25]*29.062150; /*C2H5 */
    W += c[26]*30.070120; /*C2H6 */
    W += c[27]*41.029670; /*HCCO */
    W += c[28]*42.037640; /*CH2CO */
    W += c[29]*42.037640; /*HCCOH */
    W += c[30]*28.013400; /*N2 */
    W += c[31]*39.948000; /*AR */

    for (id = 0; id < 32; ++id) {
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
    XW += x[8]*12.011150; /*C */
    XW += x[9]*13.019120; /*CH */
    XW += x[10]*14.027090; /*CH2 */
    XW += x[11]*14.027090; /*CH2(S) */
    XW += x[12]*15.035060; /*CH3 */
    XW += x[13]*16.043030; /*CH4 */
    XW += x[14]*28.010550; /*CO */
    XW += x[15]*44.009950; /*CO2 */
    XW += x[16]*29.018520; /*HCO */
    XW += x[17]*30.026490; /*CH2O */
    XW += x[18]*31.034460; /*CH2OH */
    XW += x[19]*31.034460; /*CH3O */
    XW += x[20]*32.042430; /*CH3OH */
    XW += x[21]*25.030270; /*C2H */
    XW += x[22]*26.038240; /*C2H2 */
    XW += x[23]*27.046210; /*C2H3 */
    XW += x[24]*28.054180; /*C2H4 */
    XW += x[25]*29.062150; /*C2H5 */
    XW += x[26]*30.070120; /*C2H6 */
    XW += x[27]*41.029670; /*HCCO */
    XW += x[28]*42.037640; /*CH2CO */
    XW += x[29]*42.037640; /*HCCOH */
    XW += x[30]*28.013400; /*N2 */
    XW += x[31]*39.948000; /*AR */
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
    YOW += y[8]/12.011150; /*C */
    YOW += y[9]/13.019120; /*CH */
    YOW += y[10]/14.027090; /*CH2 */
    YOW += y[11]/14.027090; /*CH2(S) */
    YOW += y[12]/15.035060; /*CH3 */
    YOW += y[13]/16.043030; /*CH4 */
    YOW += y[14]/28.010550; /*CO */
    YOW += y[15]/44.009950; /*CO2 */
    YOW += y[16]/29.018520; /*HCO */
    YOW += y[17]/30.026490; /*CH2O */
    YOW += y[18]/31.034460; /*CH2OH */
    YOW += y[19]/31.034460; /*CH3O */
    YOW += y[20]/32.042430; /*CH3OH */
    YOW += y[21]/25.030270; /*C2H */
    YOW += y[22]/26.038240; /*C2H2 */
    YOW += y[23]/27.046210; /*C2H3 */
    YOW += y[24]/28.054180; /*C2H4 */
    YOW += y[25]/29.062150; /*C2H5 */
    YOW += y[26]/30.070120; /*C2H6 */
    YOW += y[27]/41.029670; /*HCCO */
    YOW += y[28]/42.037640; /*CH2CO */
    YOW += y[29]/42.037640; /*HCCOH */
    YOW += y[30]/28.013400; /*N2 */
    YOW += y[31]/39.948000; /*AR */
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
    W += c[8]*12.011150; /*C */
    W += c[9]*13.019120; /*CH */
    W += c[10]*14.027090; /*CH2 */
    W += c[11]*14.027090; /*CH2(S) */
    W += c[12]*15.035060; /*CH3 */
    W += c[13]*16.043030; /*CH4 */
    W += c[14]*28.010550; /*CO */
    W += c[15]*44.009950; /*CO2 */
    W += c[16]*29.018520; /*HCO */
    W += c[17]*30.026490; /*CH2O */
    W += c[18]*31.034460; /*CH2OH */
    W += c[19]*31.034460; /*CH3O */
    W += c[20]*32.042430; /*CH3OH */
    W += c[21]*25.030270; /*C2H */
    W += c[22]*26.038240; /*C2H2 */
    W += c[23]*27.046210; /*C2H3 */
    W += c[24]*28.054180; /*C2H4 */
    W += c[25]*29.062150; /*C2H5 */
    W += c[26]*30.070120; /*C2H6 */
    W += c[27]*41.029670; /*HCCO */
    W += c[28]*42.037640; /*CH2CO */
    W += c[29]*42.037640; /*HCCOH */
    W += c[30]*28.013400; /*N2 */
    W += c[31]*39.948000; /*AR */

    for (id = 0; id < 32; ++id) {
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
    YOW += y[8]/12.011150; /*C */
    YOW += y[9]/13.019120; /*CH */
    YOW += y[10]/14.027090; /*CH2 */
    YOW += y[11]/14.027090; /*CH2(S) */
    YOW += y[12]/15.035060; /*CH3 */
    YOW += y[13]/16.043030; /*CH4 */
    YOW += y[14]/28.010550; /*CO */
    YOW += y[15]/44.009950; /*CO2 */
    YOW += y[16]/29.018520; /*HCO */
    YOW += y[17]/30.026490; /*CH2O */
    YOW += y[18]/31.034460; /*CH2OH */
    YOW += y[19]/31.034460; /*CH3O */
    YOW += y[20]/32.042430; /*CH3OH */
    YOW += y[21]/25.030270; /*C2H */
    YOW += y[22]/26.038240; /*C2H2 */
    YOW += y[23]/27.046210; /*C2H3 */
    YOW += y[24]/28.054180; /*C2H4 */
    YOW += y[25]/29.062150; /*C2H5 */
    YOW += y[26]/30.070120; /*C2H6 */
    YOW += y[27]/41.029670; /*HCCO */
    YOW += y[28]/42.037640; /*CH2CO */
    YOW += y[29]/42.037640; /*HCCOH */
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
    XW += x[8]*12.011150; /*C */
    XW += x[9]*13.019120; /*CH */
    XW += x[10]*14.027090; /*CH2 */
    XW += x[11]*14.027090; /*CH2(S) */
    XW += x[12]*15.035060; /*CH3 */
    XW += x[13]*16.043030; /*CH4 */
    XW += x[14]*28.010550; /*CO */
    XW += x[15]*44.009950; /*CO2 */
    XW += x[16]*29.018520; /*HCO */
    XW += x[17]*30.026490; /*CH2O */
    XW += x[18]*31.034460; /*CH2OH */
    XW += x[19]*31.034460; /*CH3O */
    XW += x[20]*32.042430; /*CH3OH */
    XW += x[21]*25.030270; /*C2H */
    XW += x[22]*26.038240; /*C2H2 */
    XW += x[23]*27.046210; /*C2H3 */
    XW += x[24]*28.054180; /*C2H4 */
    XW += x[25]*29.062150; /*C2H5 */
    XW += x[26]*30.070120; /*C2H6 */
    XW += x[27]*41.029670; /*HCCO */
    XW += x[28]*42.037640; /*CH2CO */
    XW += x[29]*42.037640; /*HCCOH */
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
    W += c[8]*12.011150; /*C */
    W += c[9]*13.019120; /*CH */
    W += c[10]*14.027090; /*CH2 */
    W += c[11]*14.027090; /*CH2(S) */
    W += c[12]*15.035060; /*CH3 */
    W += c[13]*16.043030; /*CH4 */
    W += c[14]*28.010550; /*CO */
    W += c[15]*44.009950; /*CO2 */
    W += c[16]*29.018520; /*HCO */
    W += c[17]*30.026490; /*CH2O */
    W += c[18]*31.034460; /*CH2OH */
    W += c[19]*31.034460; /*CH3O */
    W += c[20]*32.042430; /*CH3OH */
    W += c[21]*25.030270; /*C2H */
    W += c[22]*26.038240; /*C2H2 */
    W += c[23]*27.046210; /*C2H3 */
    W += c[24]*28.054180; /*C2H4 */
    W += c[25]*29.062150; /*C2H5 */
    W += c[26]*30.070120; /*C2H6 */
    W += c[27]*41.029670; /*HCCO */
    W += c[28]*42.037640; /*CH2CO */
    W += c[29]*42.037640; /*HCCOH */
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
    YOW += y[8]/12.011150; /*C */
    YOW += y[9]/13.019120; /*CH */
    YOW += y[10]/14.027090; /*CH2 */
    YOW += y[11]/14.027090; /*CH2(S) */
    YOW += y[12]/15.035060; /*CH3 */
    YOW += y[13]/16.043030; /*CH4 */
    YOW += y[14]/28.010550; /*CO */
    YOW += y[15]/44.009950; /*CO2 */
    YOW += y[16]/29.018520; /*HCO */
    YOW += y[17]/30.026490; /*CH2O */
    YOW += y[18]/31.034460; /*CH2OH */
    YOW += y[19]/31.034460; /*CH3O */
    YOW += y[20]/32.042430; /*CH3OH */
    YOW += y[21]/25.030270; /*C2H */
    YOW += y[22]/26.038240; /*C2H2 */
    YOW += y[23]/27.046210; /*C2H3 */
    YOW += y[24]/28.054180; /*C2H4 */
    YOW += y[25]/29.062150; /*C2H5 */
    YOW += y[26]/30.070120; /*C2H6 */
    YOW += y[27]/41.029670; /*HCCO */
    YOW += y[28]/42.037640; /*CH2CO */
    YOW += y[29]/42.037640; /*HCCOH */
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
    x[8] = y[8]/(12.011150*YOW); 
    x[9] = y[9]/(13.019120*YOW); 
    x[10] = y[10]/(14.027090*YOW); 
    x[11] = y[11]/(14.027090*YOW); 
    x[12] = y[12]/(15.035060*YOW); 
    x[13] = y[13]/(16.043030*YOW); 
    x[14] = y[14]/(28.010550*YOW); 
    x[15] = y[15]/(44.009950*YOW); 
    x[16] = y[16]/(29.018520*YOW); 
    x[17] = y[17]/(30.026490*YOW); 
    x[18] = y[18]/(31.034460*YOW); 
    x[19] = y[19]/(31.034460*YOW); 
    x[20] = y[20]/(32.042430*YOW); 
    x[21] = y[21]/(25.030270*YOW); 
    x[22] = y[22]/(26.038240*YOW); 
    x[23] = y[23]/(27.046210*YOW); 
    x[24] = y[24]/(28.054180*YOW); 
    x[25] = y[25]/(29.062150*YOW); 
    x[26] = y[26]/(30.070120*YOW); 
    x[27] = y[27]/(41.029670*YOW); 
    x[28] = y[28]/(42.037640*YOW); 
    x[29] = y[29]/(42.037640*YOW); 
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
    YOW += y[8]/12.011150; /*C */
    YOW += y[9]/13.019120; /*CH */
    YOW += y[10]/14.027090; /*CH2 */
    YOW += y[11]/14.027090; /*CH2(S) */
    YOW += y[12]/15.035060; /*CH3 */
    YOW += y[13]/16.043030; /*CH4 */
    YOW += y[14]/28.010550; /*CO */
    YOW += y[15]/44.009950; /*CO2 */
    YOW += y[16]/29.018520; /*HCO */
    YOW += y[17]/30.026490; /*CH2O */
    YOW += y[18]/31.034460; /*CH2OH */
    YOW += y[19]/31.034460; /*CH3O */
    YOW += y[20]/32.042430; /*CH3OH */
    YOW += y[21]/25.030270; /*C2H */
    YOW += y[22]/26.038240; /*C2H2 */
    YOW += y[23]/27.046210; /*C2H3 */
    YOW += y[24]/28.054180; /*C2H4 */
    YOW += y[25]/29.062150; /*C2H5 */
    YOW += y[26]/30.070120; /*C2H6 */
    YOW += y[27]/41.029670; /*HCCO */
    YOW += y[28]/42.037640; /*CH2CO */
    YOW += y[29]/42.037640; /*HCCOH */
    YOW += y[30]/28.013400; /*N2 */
    YOW += y[31]/39.948000; /*AR */
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
    c[8] = PWORT * y[8]/12.011150; 
    c[9] = PWORT * y[9]/13.019120; 
    c[10] = PWORT * y[10]/14.027090; 
    c[11] = PWORT * y[11]/14.027090; 
    c[12] = PWORT * y[12]/15.035060; 
    c[13] = PWORT * y[13]/16.043030; 
    c[14] = PWORT * y[14]/28.010550; 
    c[15] = PWORT * y[15]/44.009950; 
    c[16] = PWORT * y[16]/29.018520; 
    c[17] = PWORT * y[17]/30.026490; 
    c[18] = PWORT * y[18]/31.034460; 
    c[19] = PWORT * y[19]/31.034460; 
    c[20] = PWORT * y[20]/32.042430; 
    c[21] = PWORT * y[21]/25.030270; 
    c[22] = PWORT * y[22]/26.038240; 
    c[23] = PWORT * y[23]/27.046210; 
    c[24] = PWORT * y[24]/28.054180; 
    c[25] = PWORT * y[25]/29.062150; 
    c[26] = PWORT * y[26]/30.070120; 
    c[27] = PWORT * y[27]/41.029670; 
    c[28] = PWORT * y[28]/42.037640; 
    c[29] = PWORT * y[29]/42.037640; 
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
    c[8] = (*rho) * y[8]/12.011150; 
    c[9] = (*rho) * y[9]/13.019120; 
    c[10] = (*rho) * y[10]/14.027090; 
    c[11] = (*rho) * y[11]/14.027090; 
    c[12] = (*rho) * y[12]/15.035060; 
    c[13] = (*rho) * y[13]/16.043030; 
    c[14] = (*rho) * y[14]/28.010550; 
    c[15] = (*rho) * y[15]/44.009950; 
    c[16] = (*rho) * y[16]/29.018520; 
    c[17] = (*rho) * y[17]/30.026490; 
    c[18] = (*rho) * y[18]/31.034460; 
    c[19] = (*rho) * y[19]/31.034460; 
    c[20] = (*rho) * y[20]/32.042430; 
    c[21] = (*rho) * y[21]/25.030270; 
    c[22] = (*rho) * y[22]/26.038240; 
    c[23] = (*rho) * y[23]/27.046210; 
    c[24] = (*rho) * y[24]/28.054180; 
    c[25] = (*rho) * y[25]/29.062150; 
    c[26] = (*rho) * y[26]/30.070120; 
    c[27] = (*rho) * y[27]/41.029670; 
    c[28] = (*rho) * y[28]/42.037640; 
    c[29] = (*rho) * y[29]/42.037640; 
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
    XW += x[8]*12.011150; /*C */
    XW += x[9]*13.019120; /*CH */
    XW += x[10]*14.027090; /*CH2 */
    XW += x[11]*14.027090; /*CH2(S) */
    XW += x[12]*15.035060; /*CH3 */
    XW += x[13]*16.043030; /*CH4 */
    XW += x[14]*28.010550; /*CO */
    XW += x[15]*44.009950; /*CO2 */
    XW += x[16]*29.018520; /*HCO */
    XW += x[17]*30.026490; /*CH2O */
    XW += x[18]*31.034460; /*CH2OH */
    XW += x[19]*31.034460; /*CH3O */
    XW += x[20]*32.042430; /*CH3OH */
    XW += x[21]*25.030270; /*C2H */
    XW += x[22]*26.038240; /*C2H2 */
    XW += x[23]*27.046210; /*C2H3 */
    XW += x[24]*28.054180; /*C2H4 */
    XW += x[25]*29.062150; /*C2H5 */
    XW += x[26]*30.070120; /*C2H6 */
    XW += x[27]*41.029670; /*HCCO */
    XW += x[28]*42.037640; /*CH2CO */
    XW += x[29]*42.037640; /*HCCOH */
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
    y[8] = x[8]*12.011150/XW; 
    y[9] = x[9]*13.019120/XW; 
    y[10] = x[10]*14.027090/XW; 
    y[11] = x[11]*14.027090/XW; 
    y[12] = x[12]*15.035060/XW; 
    y[13] = x[13]*16.043030/XW; 
    y[14] = x[14]*28.010550/XW; 
    y[15] = x[15]*44.009950/XW; 
    y[16] = x[16]*29.018520/XW; 
    y[17] = x[17]*30.026490/XW; 
    y[18] = x[18]*31.034460/XW; 
    y[19] = x[19]*31.034460/XW; 
    y[20] = x[20]*32.042430/XW; 
    y[21] = x[21]*25.030270/XW; 
    y[22] = x[22]*26.038240/XW; 
    y[23] = x[23]*27.046210/XW; 
    y[24] = x[24]*28.054180/XW; 
    y[25] = x[25]*29.062150/XW; 
    y[26] = x[26]*30.070120/XW; 
    y[27] = x[27]*41.029670/XW; 
    y[28] = x[28]*42.037640/XW; 
    y[29] = x[29]*42.037640/XW; 
    y[30] = x[30]*28.013400/XW; 
    y[31] = x[31]*39.948000/XW; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

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
    XW += x[8]*12.011150; /*C */
    XW += x[9]*13.019120; /*CH */
    XW += x[10]*14.027090; /*CH2 */
    XW += x[11]*14.027090; /*CH2(S) */
    XW += x[12]*15.035060; /*CH3 */
    XW += x[13]*16.043030; /*CH4 */
    XW += x[14]*28.010550; /*CO */
    XW += x[15]*44.009950; /*CO2 */
    XW += x[16]*29.018520; /*HCO */
    XW += x[17]*30.026490; /*CH2O */
    XW += x[18]*31.034460; /*CH2OH */
    XW += x[19]*31.034460; /*CH3O */
    XW += x[20]*32.042430; /*CH3OH */
    XW += x[21]*25.030270; /*C2H */
    XW += x[22]*26.038240; /*C2H2 */
    XW += x[23]*27.046210; /*C2H3 */
    XW += x[24]*28.054180; /*C2H4 */
    XW += x[25]*29.062150; /*C2H5 */
    XW += x[26]*30.070120; /*C2H6 */
    XW += x[27]*41.029670; /*HCCO */
    XW += x[28]*42.037640; /*CH2CO */
    XW += x[29]*42.037640; /*HCCOH */
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
    CW += c[8]*12.011150; /*C */
    CW += c[9]*13.019120; /*CH */
    CW += c[10]*14.027090; /*CH2 */
    CW += c[11]*14.027090; /*CH2(S) */
    CW += c[12]*15.035060; /*CH3 */
    CW += c[13]*16.043030; /*CH4 */
    CW += c[14]*28.010550; /*CO */
    CW += c[15]*44.009950; /*CO2 */
    CW += c[16]*29.018520; /*HCO */
    CW += c[17]*30.026490; /*CH2O */
    CW += c[18]*31.034460; /*CH2OH */
    CW += c[19]*31.034460; /*CH3O */
    CW += c[20]*32.042430; /*CH3OH */
    CW += c[21]*25.030270; /*C2H */
    CW += c[22]*26.038240; /*C2H2 */
    CW += c[23]*27.046210; /*C2H3 */
    CW += c[24]*28.054180; /*C2H4 */
    CW += c[25]*29.062150; /*C2H5 */
    CW += c[26]*30.070120; /*C2H6 */
    CW += c[27]*41.029670; /*HCCO */
    CW += c[28]*42.037640; /*CH2CO */
    CW += c[29]*42.037640; /*HCCOH */
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
    y[8] = c[8]*12.011150/CW; 
    y[9] = c[9]*13.019120/CW; 
    y[10] = c[10]*14.027090/CW; 
    y[11] = c[11]*14.027090/CW; 
    y[12] = c[12]*15.035060/CW; 
    y[13] = c[13]*16.043030/CW; 
    y[14] = c[14]*28.010550/CW; 
    y[15] = c[15]*44.009950/CW; 
    y[16] = c[16]*29.018520/CW; 
    y[17] = c[17]*30.026490/CW; 
    y[18] = c[18]*31.034460/CW; 
    y[19] = c[19]*31.034460/CW; 
    y[20] = c[20]*32.042430/CW; 
    y[21] = c[21]*25.030270/CW; 
    y[22] = c[22]*26.038240/CW; 
    y[23] = c[23]*27.046210/CW; 
    y[24] = c[24]*28.054180/CW; 
    y[25] = c[25]*29.062150/CW; 
    y[26] = c[26]*30.070120/CW; 
    y[27] = c[27]*41.029670/CW; 
    y[28] = c[28]*42.037640/CW; 
    y[29] = c[29]*42.037640/CW; 
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
    for (id = 0; id < 32; ++id) {
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
    double RT = 8.31451e+07*tT; /*R*T */
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
    double RT = 8.31451e+07*tT; /*R*T */
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
    double RT = 8.31451e+07*tT; /*R*T */
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
    cvms[7] *= 2.444384405113783e+06; /*H2O2 */
    cvms[8] *= 6.922326338443862e+06; /*C */
    cvms[9] *= 6.386384025955671e+06; /*CH */
    cvms[10] *= 5.927466067445207e+06; /*CH2 */
    cvms[11] *= 5.927466067445207e+06; /*CH2(S) */
    cvms[12] *= 5.530081023953346e+06; /*CH3 */
    cvms[13] *= 5.182630712527496e+06; /*CH4 */
    cvms[14] *= 2.968349425484326e+06; /*CO */
    cvms[15] *= 1.889234139098090e+06; /*CO2 */
    cvms[16] *= 2.865242610581105e+06; /*HCO */
    cvms[17] *= 2.769058254894261e+06; /*CH2O */
    cvms[18] *= 2.679121853578248e+06; /*CH2OH */
    cvms[19] *= 2.679121853578248e+06; /*CH3O */
    cvms[20] *= 2.594843774332970e+06; /*CH3OH */
    cvms[21] *= 3.321781986370902e+06; /*C2H */
    cvms[22] *= 3.193192012977835e+06; /*C2H2 */
    cvms[23] *= 3.074186734481467e+06; /*C2H3 */
    cvms[24] *= 2.963733033722604e+06; /*C2H4 */
    cvms[25] *= 2.860941121011349e+06; /*C2H5 */
    cvms[26] *= 2.765040511976673e+06; /*C2H6 */
    cvms[27] *= 2.026462801187531e+06; /*HCCO */
    cvms[28] *= 1.977872687429646e+06; /*CH2CO */
    cvms[29] *= 1.977872687429646e+06; /*HCCOH */
    cvms[30] *= 2.968047434442088e+06; /*N2 */
    cvms[31] *= 2.081333233203164e+06; /*AR */
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
    cpms[7] *= 2.444384405113783e+06; /*H2O2 */
    cpms[8] *= 6.922326338443862e+06; /*C */
    cpms[9] *= 6.386384025955671e+06; /*CH */
    cpms[10] *= 5.927466067445207e+06; /*CH2 */
    cpms[11] *= 5.927466067445207e+06; /*CH2(S) */
    cpms[12] *= 5.530081023953346e+06; /*CH3 */
    cpms[13] *= 5.182630712527496e+06; /*CH4 */
    cpms[14] *= 2.968349425484326e+06; /*CO */
    cpms[15] *= 1.889234139098090e+06; /*CO2 */
    cpms[16] *= 2.865242610581105e+06; /*HCO */
    cpms[17] *= 2.769058254894261e+06; /*CH2O */
    cpms[18] *= 2.679121853578248e+06; /*CH2OH */
    cpms[19] *= 2.679121853578248e+06; /*CH3O */
    cpms[20] *= 2.594843774332970e+06; /*CH3OH */
    cpms[21] *= 3.321781986370902e+06; /*C2H */
    cpms[22] *= 3.193192012977835e+06; /*C2H2 */
    cpms[23] *= 3.074186734481467e+06; /*C2H3 */
    cpms[24] *= 2.963733033722604e+06; /*C2H4 */
    cpms[25] *= 2.860941121011349e+06; /*C2H5 */
    cpms[26] *= 2.765040511976673e+06; /*C2H6 */
    cpms[27] *= 2.026462801187531e+06; /*HCCO */
    cpms[28] *= 1.977872687429646e+06; /*CH2CO */
    cpms[29] *= 1.977872687429646e+06; /*HCCOH */
    cpms[30] *= 2.968047434442088e+06; /*N2 */
    cpms[31] *= 2.081333233203164e+06; /*AR */
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
    ums[8] *= RT/12.011150; /*C */
    ums[9] *= RT/13.019120; /*CH */
    ums[10] *= RT/14.027090; /*CH2 */
    ums[11] *= RT/14.027090; /*CH2(S) */
    ums[12] *= RT/15.035060; /*CH3 */
    ums[13] *= RT/16.043030; /*CH4 */
    ums[14] *= RT/28.010550; /*CO */
    ums[15] *= RT/44.009950; /*CO2 */
    ums[16] *= RT/29.018520; /*HCO */
    ums[17] *= RT/30.026490; /*CH2O */
    ums[18] *= RT/31.034460; /*CH2OH */
    ums[19] *= RT/31.034460; /*CH3O */
    ums[20] *= RT/32.042430; /*CH3OH */
    ums[21] *= RT/25.030270; /*C2H */
    ums[22] *= RT/26.038240; /*C2H2 */
    ums[23] *= RT/27.046210; /*C2H3 */
    ums[24] *= RT/28.054180; /*C2H4 */
    ums[25] *= RT/29.062150; /*C2H5 */
    ums[26] *= RT/30.070120; /*C2H6 */
    ums[27] *= RT/41.029670; /*HCCO */
    ums[28] *= RT/42.037640; /*CH2CO */
    ums[29] *= RT/42.037640; /*HCCOH */
    ums[30] *= RT/28.013400; /*N2 */
    ums[31] *= RT/39.948000; /*AR */
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
    hms[8] *= RT/12.011150; /*C */
    hms[9] *= RT/13.019120; /*CH */
    hms[10] *= RT/14.027090; /*CH2 */
    hms[11] *= RT/14.027090; /*CH2(S) */
    hms[12] *= RT/15.035060; /*CH3 */
    hms[13] *= RT/16.043030; /*CH4 */
    hms[14] *= RT/28.010550; /*CO */
    hms[15] *= RT/44.009950; /*CO2 */
    hms[16] *= RT/29.018520; /*HCO */
    hms[17] *= RT/30.026490; /*CH2O */
    hms[18] *= RT/31.034460; /*CH2OH */
    hms[19] *= RT/31.034460; /*CH3O */
    hms[20] *= RT/32.042430; /*CH3OH */
    hms[21] *= RT/25.030270; /*C2H */
    hms[22] *= RT/26.038240; /*C2H2 */
    hms[23] *= RT/27.046210; /*C2H3 */
    hms[24] *= RT/28.054180; /*C2H4 */
    hms[25] *= RT/29.062150; /*C2H5 */
    hms[26] *= RT/30.070120; /*C2H6 */
    hms[27] *= RT/41.029670; /*HCCO */
    hms[28] *= RT/42.037640; /*CH2CO */
    hms[29] *= RT/42.037640; /*HCCOH */
    hms[30] *= RT/28.013400; /*N2 */
    hms[31] *= RT/39.948000; /*AR */
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
    gms[8] *= RT/12.011150; /*C */
    gms[9] *= RT/13.019120; /*CH */
    gms[10] *= RT/14.027090; /*CH2 */
    gms[11] *= RT/14.027090; /*CH2(S) */
    gms[12] *= RT/15.035060; /*CH3 */
    gms[13] *= RT/16.043030; /*CH4 */
    gms[14] *= RT/28.010550; /*CO */
    gms[15] *= RT/44.009950; /*CO2 */
    gms[16] *= RT/29.018520; /*HCO */
    gms[17] *= RT/30.026490; /*CH2O */
    gms[18] *= RT/31.034460; /*CH2OH */
    gms[19] *= RT/31.034460; /*CH3O */
    gms[20] *= RT/32.042430; /*CH3OH */
    gms[21] *= RT/25.030270; /*C2H */
    gms[22] *= RT/26.038240; /*C2H2 */
    gms[23] *= RT/27.046210; /*C2H3 */
    gms[24] *= RT/28.054180; /*C2H4 */
    gms[25] *= RT/29.062150; /*C2H5 */
    gms[26] *= RT/30.070120; /*C2H6 */
    gms[27] *= RT/41.029670; /*HCCO */
    gms[28] *= RT/42.037640; /*CH2CO */
    gms[29] *= RT/42.037640; /*HCCOH */
    gms[30] *= RT/28.013400; /*N2 */
    gms[31] *= RT/39.948000; /*AR */
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
    ams[8] *= RT/12.011150; /*C */
    ams[9] *= RT/13.019120; /*CH */
    ams[10] *= RT/14.027090; /*CH2 */
    ams[11] *= RT/14.027090; /*CH2(S) */
    ams[12] *= RT/15.035060; /*CH3 */
    ams[13] *= RT/16.043030; /*CH4 */
    ams[14] *= RT/28.010550; /*CO */
    ams[15] *= RT/44.009950; /*CO2 */
    ams[16] *= RT/29.018520; /*HCO */
    ams[17] *= RT/30.026490; /*CH2O */
    ams[18] *= RT/31.034460; /*CH2OH */
    ams[19] *= RT/31.034460; /*CH3O */
    ams[20] *= RT/32.042430; /*CH3OH */
    ams[21] *= RT/25.030270; /*C2H */
    ams[22] *= RT/26.038240; /*C2H2 */
    ams[23] *= RT/27.046210; /*C2H3 */
    ams[24] *= RT/28.054180; /*C2H4 */
    ams[25] *= RT/29.062150; /*C2H5 */
    ams[26] *= RT/30.070120; /*C2H6 */
    ams[27] *= RT/41.029670; /*HCCO */
    ams[28] *= RT/42.037640; /*CH2CO */
    ams[29] *= RT/42.037640; /*HCCOH */
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
    sms[0] *= 4.124383662212169e+07; /*H2 */
    sms[1] *= 8.248767324424338e+07; /*H */
    sms[2] *= 5.196763628636074e+06; /*O */
    sms[3] *= 2.598381814318037e+06; /*O2 */
    sms[4] *= 4.888768810227566e+06; /*OH */
    sms[5] *= 4.615239012974499e+06; /*H2O */
    sms[6] *= 2.519031701678171e+06; /*HO2 */
    sms[7] *= 2.444384405113783e+06; /*H2O2 */
    sms[8] *= 6.922326338443862e+06; /*C */
    sms[9] *= 6.386384025955671e+06; /*CH */
    sms[10] *= 5.927466067445207e+06; /*CH2 */
    sms[11] *= 5.927466067445207e+06; /*CH2(S) */
    sms[12] *= 5.530081023953346e+06; /*CH3 */
    sms[13] *= 5.182630712527496e+06; /*CH4 */
    sms[14] *= 2.968349425484326e+06; /*CO */
    sms[15] *= 1.889234139098090e+06; /*CO2 */
    sms[16] *= 2.865242610581105e+06; /*HCO */
    sms[17] *= 2.769058254894261e+06; /*CH2O */
    sms[18] *= 2.679121853578248e+06; /*CH2OH */
    sms[19] *= 2.679121853578248e+06; /*CH3O */
    sms[20] *= 2.594843774332970e+06; /*CH3OH */
    sms[21] *= 3.321781986370902e+06; /*C2H */
    sms[22] *= 3.193192012977835e+06; /*C2H2 */
    sms[23] *= 3.074186734481467e+06; /*C2H3 */
    sms[24] *= 2.963733033722604e+06; /*C2H4 */
    sms[25] *= 2.860941121011349e+06; /*C2H5 */
    sms[26] *= 2.765040511976673e+06; /*C2H6 */
    sms[27] *= 2.026462801187531e+06; /*HCCO */
    sms[28] *= 1.977872687429646e+06; /*CH2CO */
    sms[29] *= 1.977872687429646e+06; /*HCCOH */
    sms[30] *= 2.968047434442088e+06; /*N2 */
    sms[31] *= 2.081333233203164e+06; /*AR */
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

    *cpbl = result * 8.31451e+07;
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
    result += cpor[0]*y[0]/2.015940; /*H2 */
    result += cpor[1]*y[1]/1.007970; /*H */
    result += cpor[2]*y[2]/15.999400; /*O */
    result += cpor[3]*y[3]/31.998800; /*O2 */
    result += cpor[4]*y[4]/17.007370; /*OH */
    result += cpor[5]*y[5]/18.015340; /*H2O */
    result += cpor[6]*y[6]/33.006770; /*HO2 */
    result += cpor[7]*y[7]/34.014740; /*H2O2 */
    result += cpor[8]*y[8]/12.011150; /*C */
    result += cpor[9]*y[9]/13.019120; /*CH */
    result += cpor[10]*y[10]/14.027090; /*CH2 */
    result += cpor[11]*y[11]/14.027090; /*CH2(S) */
    result += cpor[12]*y[12]/15.035060; /*CH3 */
    result += cpor[13]*y[13]/16.043030; /*CH4 */
    result += cpor[14]*y[14]/28.010550; /*CO */
    result += cpor[15]*y[15]/44.009950; /*CO2 */
    result += cpor[16]*y[16]/29.018520; /*HCO */
    result += cpor[17]*y[17]/30.026490; /*CH2O */
    result += cpor[18]*y[18]/31.034460; /*CH2OH */
    result += cpor[19]*y[19]/31.034460; /*CH3O */
    result += cpor[20]*y[20]/32.042430; /*CH3OH */
    result += cpor[21]*y[21]/25.030270; /*C2H */
    result += cpor[22]*y[22]/26.038240; /*C2H2 */
    result += cpor[23]*y[23]/27.046210; /*C2H3 */
    result += cpor[24]*y[24]/28.054180; /*C2H4 */
    result += cpor[25]*y[25]/29.062150; /*C2H5 */
    result += cpor[26]*y[26]/30.070120; /*C2H6 */
    result += cpor[27]*y[27]/41.029670; /*HCCO */
    result += cpor[28]*y[28]/42.037640; /*CH2CO */
    result += cpor[29]*y[29]/42.037640; /*HCCOH */
    result += cpor[30]*y[30]/28.013400; /*N2 */
    result += cpor[31]*y[31]/39.948000; /*AR */

    *cpbs = result * 8.31451e+07;
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

    *cvbl = result * 8.31451e+07;
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
    result += cvor[0]*y[0]/2.015940; /*H2 */
    result += cvor[1]*y[1]/1.007970; /*H */
    result += cvor[2]*y[2]/15.999400; /*O */
    result += cvor[3]*y[3]/31.998800; /*O2 */
    result += cvor[4]*y[4]/17.007370; /*OH */
    result += cvor[5]*y[5]/18.015340; /*H2O */
    result += cvor[6]*y[6]/33.006770; /*HO2 */
    result += cvor[7]*y[7]/34.014740; /*H2O2 */
    result += cvor[8]*y[8]/12.011150; /*C */
    result += cvor[9]*y[9]/13.019120; /*CH */
    result += cvor[10]*y[10]/14.027090; /*CH2 */
    result += cvor[11]*y[11]/14.027090; /*CH2(S) */
    result += cvor[12]*y[12]/15.035060; /*CH3 */
    result += cvor[13]*y[13]/16.043030; /*CH4 */
    result += cvor[14]*y[14]/28.010550; /*CO */
    result += cvor[15]*y[15]/44.009950; /*CO2 */
    result += cvor[16]*y[16]/29.018520; /*HCO */
    result += cvor[17]*y[17]/30.026490; /*CH2O */
    result += cvor[18]*y[18]/31.034460; /*CH2OH */
    result += cvor[19]*y[19]/31.034460; /*CH3O */
    result += cvor[20]*y[20]/32.042430; /*CH3OH */
    result += cvor[21]*y[21]/25.030270; /*C2H */
    result += cvor[22]*y[22]/26.038240; /*C2H2 */
    result += cvor[23]*y[23]/27.046210; /*C2H3 */
    result += cvor[24]*y[24]/28.054180; /*C2H4 */
    result += cvor[25]*y[25]/29.062150; /*C2H5 */
    result += cvor[26]*y[26]/30.070120; /*C2H6 */
    result += cvor[27]*y[27]/41.029670; /*HCCO */
    result += cvor[28]*y[28]/42.037640; /*CH2CO */
    result += cvor[29]*y[29]/42.037640; /*HCCOH */
    result += cvor[30]*y[30]/28.013400; /*N2 */
    result += cvor[31]*y[31]/39.948000; /*AR */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double *T, double *x, int * iwrk, double * rwrk, double * hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[32]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
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
    result += y[8]*hml[8]/12.011150; /*C */
    result += y[9]*hml[9]/13.019120; /*CH */
    result += y[10]*hml[10]/14.027090; /*CH2 */
    result += y[11]*hml[11]/14.027090; /*CH2(S) */
    result += y[12]*hml[12]/15.035060; /*CH3 */
    result += y[13]*hml[13]/16.043030; /*CH4 */
    result += y[14]*hml[14]/28.010550; /*CO */
    result += y[15]*hml[15]/44.009950; /*CO2 */
    result += y[16]*hml[16]/29.018520; /*HCO */
    result += y[17]*hml[17]/30.026490; /*CH2O */
    result += y[18]*hml[18]/31.034460; /*CH2OH */
    result += y[19]*hml[19]/31.034460; /*CH3O */
    result += y[20]*hml[20]/32.042430; /*CH3OH */
    result += y[21]*hml[21]/25.030270; /*C2H */
    result += y[22]*hml[22]/26.038240; /*C2H2 */
    result += y[23]*hml[23]/27.046210; /*C2H3 */
    result += y[24]*hml[24]/28.054180; /*C2H4 */
    result += y[25]*hml[25]/29.062150; /*C2H5 */
    result += y[26]*hml[26]/30.070120; /*C2H6 */
    result += y[27]*hml[27]/41.029670; /*HCCO */
    result += y[28]*hml[28]/42.037640; /*CH2CO */
    result += y[29]*hml[29]/42.037640; /*HCCOH */
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
    double RT = 8.31451e+07*tT; /*R*T */
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
    result += y[8]*ums[8]/12.011150; /*C */
    result += y[9]*ums[9]/13.019120; /*CH */
    result += y[10]*ums[10]/14.027090; /*CH2 */
    result += y[11]*ums[11]/14.027090; /*CH2(S) */
    result += y[12]*ums[12]/15.035060; /*CH3 */
    result += y[13]*ums[13]/16.043030; /*CH4 */
    result += y[14]*ums[14]/28.010550; /*CO */
    result += y[15]*ums[15]/44.009950; /*CO2 */
    result += y[16]*ums[16]/29.018520; /*HCO */
    result += y[17]*ums[17]/30.026490; /*CH2O */
    result += y[18]*ums[18]/31.034460; /*CH2OH */
    result += y[19]*ums[19]/31.034460; /*CH3O */
    result += y[20]*ums[20]/32.042430; /*CH3OH */
    result += y[21]*ums[21]/25.030270; /*C2H */
    result += y[22]*ums[22]/26.038240; /*C2H2 */
    result += y[23]*ums[23]/27.046210; /*C2H3 */
    result += y[24]*ums[24]/28.054180; /*C2H4 */
    result += y[25]*ums[25]/29.062150; /*C2H5 */
    result += y[26]*ums[26]/30.070120; /*C2H6 */
    result += y[27]*ums[27]/41.029670; /*HCCO */
    result += y[28]*ums[28]/42.037640; /*CH2CO */
    result += y[29]*ums[29]/42.037640; /*HCCOH */
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
    YOW += y[8]/12.011150; /*C */
    YOW += y[9]/13.019120; /*CH */
    YOW += y[10]/14.027090; /*CH2 */
    YOW += y[11]/14.027090; /*CH2(S) */
    YOW += y[12]/15.035060; /*CH3 */
    YOW += y[13]/16.043030; /*CH4 */
    YOW += y[14]/28.010550; /*CO */
    YOW += y[15]/44.009950; /*CO2 */
    YOW += y[16]/29.018520; /*HCO */
    YOW += y[17]/30.026490; /*CH2O */
    YOW += y[18]/31.034460; /*CH2OH */
    YOW += y[19]/31.034460; /*CH3O */
    YOW += y[20]/32.042430; /*CH3OH */
    YOW += y[21]/25.030270; /*C2H */
    YOW += y[22]/26.038240; /*C2H2 */
    YOW += y[23]/27.046210; /*C2H3 */
    YOW += y[24]/28.054180; /*C2H4 */
    YOW += y[25]/29.062150; /*C2H5 */
    YOW += y[26]/30.070120; /*C2H6 */
    YOW += y[27]/41.029670; /*HCCO */
    YOW += y[28]/42.037640; /*CH2CO */
    YOW += y[29]/42.037640; /*HCCOH */
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
    x[8] = y[8]/(12.011150*YOW); 
    x[9] = y[9]/(13.019120*YOW); 
    x[10] = y[10]/(14.027090*YOW); 
    x[11] = y[11]/(14.027090*YOW); 
    x[12] = y[12]/(15.035060*YOW); 
    x[13] = y[13]/(16.043030*YOW); 
    x[14] = y[14]/(28.010550*YOW); 
    x[15] = y[15]/(44.009950*YOW); 
    x[16] = y[16]/(29.018520*YOW); 
    x[17] = y[17]/(30.026490*YOW); 
    x[18] = y[18]/(31.034460*YOW); 
    x[19] = y[19]/(31.034460*YOW); 
    x[20] = y[20]/(32.042430*YOW); 
    x[21] = y[21]/(25.030270*YOW); 
    x[22] = y[22]/(26.038240*YOW); 
    x[23] = y[23]/(27.046210*YOW); 
    x[24] = y[24]/(28.054180*YOW); 
    x[25] = y[25]/(29.062150*YOW); 
    x[26] = y[26]/(30.070120*YOW); 
    x[27] = y[27]/(41.029670*YOW); 
    x[28] = y[28]/(42.037640*YOW); 
    x[29] = y[29]/(42.037640*YOW); 
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
    double RT = 8.31451e+07*tT; /*R*T */
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
    YOW += y[8]/12.011150; /*C */
    YOW += y[9]/13.019120; /*CH */
    YOW += y[10]/14.027090; /*CH2 */
    YOW += y[11]/14.027090; /*CH2(S) */
    YOW += y[12]/15.035060; /*CH3 */
    YOW += y[13]/16.043030; /*CH4 */
    YOW += y[14]/28.010550; /*CO */
    YOW += y[15]/44.009950; /*CO2 */
    YOW += y[16]/29.018520; /*HCO */
    YOW += y[17]/30.026490; /*CH2O */
    YOW += y[18]/31.034460; /*CH2OH */
    YOW += y[19]/31.034460; /*CH3O */
    YOW += y[20]/32.042430; /*CH3OH */
    YOW += y[21]/25.030270; /*C2H */
    YOW += y[22]/26.038240; /*C2H2 */
    YOW += y[23]/27.046210; /*C2H3 */
    YOW += y[24]/28.054180; /*C2H4 */
    YOW += y[25]/29.062150; /*C2H5 */
    YOW += y[26]/30.070120; /*C2H6 */
    YOW += y[27]/41.029670; /*HCCO */
    YOW += y[28]/42.037640; /*CH2CO */
    YOW += y[29]/42.037640; /*HCCOH */
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
    x[8] = y[8]/(12.011150*YOW); 
    x[9] = y[9]/(13.019120*YOW); 
    x[10] = y[10]/(14.027090*YOW); 
    x[11] = y[11]/(14.027090*YOW); 
    x[12] = y[12]/(15.035060*YOW); 
    x[13] = y[13]/(16.043030*YOW); 
    x[14] = y[14]/(28.010550*YOW); 
    x[15] = y[15]/(44.009950*YOW); 
    x[16] = y[16]/(29.018520*YOW); 
    x[17] = y[17]/(30.026490*YOW); 
    x[18] = y[18]/(31.034460*YOW); 
    x[19] = y[19]/(31.034460*YOW); 
    x[20] = y[20]/(32.042430*YOW); 
    x[21] = y[21]/(25.030270*YOW); 
    x[22] = y[22]/(26.038240*YOW); 
    x[23] = y[23]/(27.046210*YOW); 
    x[24] = y[24]/(28.054180*YOW); 
    x[25] = y[25]/(29.062150*YOW); 
    x[26] = y[26]/(30.070120*YOW); 
    x[27] = y[27]/(41.029670*YOW); 
    x[28] = y[28]/(42.037640*YOW); 
    x[29] = y[29]/(42.037640*YOW); 
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
    double RT = 8.31451e+07*tT; /*R*T */
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
    double RT = 8.31451e+07*tT; /*R*T */
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
    YOW += y[8]/12.011150; /*C */
    YOW += y[9]/13.019120; /*CH */
    YOW += y[10]/14.027090; /*CH2 */
    YOW += y[11]/14.027090; /*CH2(S) */
    YOW += y[12]/15.035060; /*CH3 */
    YOW += y[13]/16.043030; /*CH4 */
    YOW += y[14]/28.010550; /*CO */
    YOW += y[15]/44.009950; /*CO2 */
    YOW += y[16]/29.018520; /*HCO */
    YOW += y[17]/30.026490; /*CH2O */
    YOW += y[18]/31.034460; /*CH2OH */
    YOW += y[19]/31.034460; /*CH3O */
    YOW += y[20]/32.042430; /*CH3OH */
    YOW += y[21]/25.030270; /*C2H */
    YOW += y[22]/26.038240; /*C2H2 */
    YOW += y[23]/27.046210; /*C2H3 */
    YOW += y[24]/28.054180; /*C2H4 */
    YOW += y[25]/29.062150; /*C2H5 */
    YOW += y[26]/30.070120; /*C2H6 */
    YOW += y[27]/41.029670; /*HCCO */
    YOW += y[28]/42.037640; /*CH2CO */
    YOW += y[29]/42.037640; /*HCCOH */
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
    x[8] = y[8]/(12.011150*YOW); 
    x[9] = y[9]/(13.019120*YOW); 
    x[10] = y[10]/(14.027090*YOW); 
    x[11] = y[11]/(14.027090*YOW); 
    x[12] = y[12]/(15.035060*YOW); 
    x[13] = y[13]/(16.043030*YOW); 
    x[14] = y[14]/(28.010550*YOW); 
    x[15] = y[15]/(44.009950*YOW); 
    x[16] = y[16]/(29.018520*YOW); 
    x[17] = y[17]/(30.026490*YOW); 
    x[18] = y[18]/(31.034460*YOW); 
    x[19] = y[19]/(31.034460*YOW); 
    x[20] = y[20]/(32.042430*YOW); 
    x[21] = y[21]/(25.030270*YOW); 
    x[22] = y[22]/(26.038240*YOW); 
    x[23] = y[23]/(27.046210*YOW); 
    x[24] = y[24]/(28.054180*YOW); 
    x[25] = y[25]/(29.062150*YOW); 
    x[26] = y[26]/(30.070120*YOW); 
    x[27] = y[27]/(41.029670*YOW); 
    x[28] = y[28]/(42.037640*YOW); 
    x[29] = y[29]/(42.037640*YOW); 
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
    YOW += y[8]/12.011150; /*C */
    YOW += y[9]/13.019120; /*CH */
    YOW += y[10]/14.027090; /*CH2 */
    YOW += y[11]/14.027090; /*CH2(S) */
    YOW += y[12]/15.035060; /*CH3 */
    YOW += y[13]/16.043030; /*CH4 */
    YOW += y[14]/28.010550; /*CO */
    YOW += y[15]/44.009950; /*CO2 */
    YOW += y[16]/29.018520; /*HCO */
    YOW += y[17]/30.026490; /*CH2O */
    YOW += y[18]/31.034460; /*CH2OH */
    YOW += y[19]/31.034460; /*CH3O */
    YOW += y[20]/32.042430; /*CH3OH */
    YOW += y[21]/25.030270; /*C2H */
    YOW += y[22]/26.038240; /*C2H2 */
    YOW += y[23]/27.046210; /*C2H3 */
    YOW += y[24]/28.054180; /*C2H4 */
    YOW += y[25]/29.062150; /*C2H5 */
    YOW += y[26]/30.070120; /*C2H6 */
    YOW += y[27]/41.029670; /*HCCO */
    YOW += y[28]/42.037640; /*CH2CO */
    YOW += y[29]/42.037640; /*HCCOH */
    YOW += y[30]/28.013400; /*N2 */
    YOW += y[31]/39.948000; /*AR */
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
    c[8] = PWORT * y[8]/12.011150; 
    c[9] = PWORT * y[9]/13.019120; 
    c[10] = PWORT * y[10]/14.027090; 
    c[11] = PWORT * y[11]/14.027090; 
    c[12] = PWORT * y[12]/15.035060; 
    c[13] = PWORT * y[13]/16.043030; 
    c[14] = PWORT * y[14]/28.010550; 
    c[15] = PWORT * y[15]/44.009950; 
    c[16] = PWORT * y[16]/29.018520; 
    c[17] = PWORT * y[17]/30.026490; 
    c[18] = PWORT * y[18]/31.034460; 
    c[19] = PWORT * y[19]/31.034460; 
    c[20] = PWORT * y[20]/32.042430; 
    c[21] = PWORT * y[21]/25.030270; 
    c[22] = PWORT * y[22]/26.038240; 
    c[23] = PWORT * y[23]/27.046210; 
    c[24] = PWORT * y[24]/28.054180; 
    c[25] = PWORT * y[25]/29.062150; 
    c[26] = PWORT * y[26]/30.070120; 
    c[27] = PWORT * y[27]/41.029670; 
    c[28] = PWORT * y[28]/42.037640; 
    c[29] = PWORT * y[29]/42.037640; 
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
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

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
    c[8] = 1e6 * (*rho) * y[8]/12.011150; 
    c[9] = 1e6 * (*rho) * y[9]/13.019120; 
    c[10] = 1e6 * (*rho) * y[10]/14.027090; 
    c[11] = 1e6 * (*rho) * y[11]/14.027090; 
    c[12] = 1e6 * (*rho) * y[12]/15.035060; 
    c[13] = 1e6 * (*rho) * y[13]/16.043030; 
    c[14] = 1e6 * (*rho) * y[14]/28.010550; 
    c[15] = 1e6 * (*rho) * y[15]/44.009950; 
    c[16] = 1e6 * (*rho) * y[16]/29.018520; 
    c[17] = 1e6 * (*rho) * y[17]/30.026490; 
    c[18] = 1e6 * (*rho) * y[18]/31.034460; 
    c[19] = 1e6 * (*rho) * y[19]/31.034460; 
    c[20] = 1e6 * (*rho) * y[20]/32.042430; 
    c[21] = 1e6 * (*rho) * y[21]/25.030270; 
    c[22] = 1e6 * (*rho) * y[22]/26.038240; 
    c[23] = 1e6 * (*rho) * y[23]/27.046210; 
    c[24] = 1e6 * (*rho) * y[24]/28.054180; 
    c[25] = 1e6 * (*rho) * y[25]/29.062150; 
    c[26] = 1e6 * (*rho) * y[26]/30.070120; 
    c[27] = 1e6 * (*rho) * y[27]/41.029670; 
    c[28] = 1e6 * (*rho) * y[28]/42.037640; 
    c[29] = 1e6 * (*rho) * y[29]/42.037640; 
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
    XW += x[8]*12.011150; /*C */
    XW += x[9]*13.019120; /*CH */
    XW += x[10]*14.027090; /*CH2 */
    XW += x[11]*14.027090; /*CH2(S) */
    XW += x[12]*15.035060; /*CH3 */
    XW += x[13]*16.043030; /*CH4 */
    XW += x[14]*28.010550; /*CO */
    XW += x[15]*44.009950; /*CO2 */
    XW += x[16]*29.018520; /*HCO */
    XW += x[17]*30.026490; /*CH2O */
    XW += x[18]*31.034460; /*CH2OH */
    XW += x[19]*31.034460; /*CH3O */
    XW += x[20]*32.042430; /*CH3OH */
    XW += x[21]*25.030270; /*C2H */
    XW += x[22]*26.038240; /*C2H2 */
    XW += x[23]*27.046210; /*C2H3 */
    XW += x[24]*28.054180; /*C2H4 */
    XW += x[25]*29.062150; /*C2H5 */
    XW += x[26]*30.070120; /*C2H6 */
    XW += x[27]*41.029670; /*HCCO */
    XW += x[28]*42.037640; /*CH2CO */
    XW += x[29]*42.037640; /*HCCOH */
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
/*Given P, T, and mole fractions */
void CKKFKR(double * P, double * T, double * x, int * iwrk, double * rwrk, double * q_f, double * q_r)
{
    int id; /*loop counter */
    double c[32]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 32; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 177; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
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
    YOW += y[8]/12.011150; /*C */
    YOW += y[9]/13.019120; /*CH */
    YOW += y[10]/14.027090; /*CH2 */
    YOW += y[11]/14.027090; /*CH2(S) */
    YOW += y[12]/15.035060; /*CH3 */
    YOW += y[13]/16.043030; /*CH4 */
    YOW += y[14]/28.010550; /*CO */
    YOW += y[15]/44.009950; /*CO2 */
    YOW += y[16]/29.018520; /*HCO */
    YOW += y[17]/30.026490; /*CH2O */
    YOW += y[18]/31.034460; /*CH2OH */
    YOW += y[19]/31.034460; /*CH3O */
    YOW += y[20]/32.042430; /*CH3OH */
    YOW += y[21]/25.030270; /*C2H */
    YOW += y[22]/26.038240; /*C2H2 */
    YOW += y[23]/27.046210; /*C2H3 */
    YOW += y[24]/28.054180; /*C2H4 */
    YOW += y[25]/29.062150; /*C2H5 */
    YOW += y[26]/30.070120; /*C2H6 */
    YOW += y[27]/41.029670; /*HCCO */
    YOW += y[28]/42.037640; /*CH2CO */
    YOW += y[29]/42.037640; /*HCCOH */
    YOW += y[30]/28.013400; /*N2 */
    YOW += y[31]/39.948000; /*AR */
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
    c[8] = PWORT * y[8]/12.011150; 
    c[9] = PWORT * y[9]/13.019120; 
    c[10] = PWORT * y[10]/14.027090; 
    c[11] = PWORT * y[11]/14.027090; 
    c[12] = PWORT * y[12]/15.035060; 
    c[13] = PWORT * y[13]/16.043030; 
    c[14] = PWORT * y[14]/28.010550; 
    c[15] = PWORT * y[15]/44.009950; 
    c[16] = PWORT * y[16]/29.018520; 
    c[17] = PWORT * y[17]/30.026490; 
    c[18] = PWORT * y[18]/31.034460; 
    c[19] = PWORT * y[19]/31.034460; 
    c[20] = PWORT * y[20]/32.042430; 
    c[21] = PWORT * y[21]/25.030270; 
    c[22] = PWORT * y[22]/26.038240; 
    c[23] = PWORT * y[23]/27.046210; 
    c[24] = PWORT * y[24]/28.054180; 
    c[25] = PWORT * y[25]/29.062150; 
    c[26] = PWORT * y[26]/30.070120; 
    c[27] = PWORT * y[27]/41.029670; 
    c[28] = PWORT * y[28]/42.037640; 
    c[29] = PWORT * y[29]/42.037640; 
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
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

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
    c[8] = 1e6 * (*rho) * y[8]/12.011150; 
    c[9] = 1e6 * (*rho) * y[9]/13.019120; 
    c[10] = 1e6 * (*rho) * y[10]/14.027090; 
    c[11] = 1e6 * (*rho) * y[11]/14.027090; 
    c[12] = 1e6 * (*rho) * y[12]/15.035060; 
    c[13] = 1e6 * (*rho) * y[13]/16.043030; 
    c[14] = 1e6 * (*rho) * y[14]/28.010550; 
    c[15] = 1e6 * (*rho) * y[15]/44.009950; 
    c[16] = 1e6 * (*rho) * y[16]/29.018520; 
    c[17] = 1e6 * (*rho) * y[17]/30.026490; 
    c[18] = 1e6 * (*rho) * y[18]/31.034460; 
    c[19] = 1e6 * (*rho) * y[19]/31.034460; 
    c[20] = 1e6 * (*rho) * y[20]/32.042430; 
    c[21] = 1e6 * (*rho) * y[21]/25.030270; 
    c[22] = 1e6 * (*rho) * y[22]/26.038240; 
    c[23] = 1e6 * (*rho) * y[23]/27.046210; 
    c[24] = 1e6 * (*rho) * y[24]/28.054180; 
    c[25] = 1e6 * (*rho) * y[25]/29.062150; 
    c[26] = 1e6 * (*rho) * y[26]/30.070120; 
    c[27] = 1e6 * (*rho) * y[27]/41.029670; 
    c[28] = 1e6 * (*rho) * y[28]/42.037640; 
    c[29] = 1e6 * (*rho) * y[29]/42.037640; 
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
    XW += x[8]*12.011150; /*C */
    XW += x[9]*13.019120; /*CH */
    XW += x[10]*14.027090; /*CH2 */
    XW += x[11]*14.027090; /*CH2(S) */
    XW += x[12]*15.035060; /*CH3 */
    XW += x[13]*16.043030; /*CH4 */
    XW += x[14]*28.010550; /*CO */
    XW += x[15]*44.009950; /*CO2 */
    XW += x[16]*29.018520; /*HCO */
    XW += x[17]*30.026490; /*CH2O */
    XW += x[18]*31.034460; /*CH2OH */
    XW += x[19]*31.034460; /*CH3O */
    XW += x[20]*32.042430; /*CH3OH */
    XW += x[21]*25.030270; /*C2H */
    XW += x[22]*26.038240; /*C2H2 */
    XW += x[23]*27.046210; /*C2H3 */
    XW += x[24]*28.054180; /*C2H4 */
    XW += x[25]*29.062150; /*C2H5 */
    XW += x[26]*30.070120; /*C2H6 */
    XW += x[27]*41.029670; /*HCCO */
    XW += x[28]*42.037640; /*CH2CO */
    XW += x[29]*42.037640; /*HCCOH */
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
    for (id = 0; id < 32 * kd; ++ id) {
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

static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[177];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[177];
#ifdef _OPENMP
#pragma omp threadprivate(Kc_save)
#endif

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

    double invT = 1.0 / tc[1];

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

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

    if (T != T_save)
    {
        T_save = T;

        k_f_save[0] = 1e-12 * 1.2e+17*exp(-1*tc[0]);
        k_f_save[1] = 1e-12 * 5e+17*exp(-1*tc[0]);
        k_f_save[2] = 1e-06 * 50000*exp(2.67*tc[0]-3165.2328279116868543*invT);
        k_f_save[3] = 1e-06 * 2e+13;
        k_f_save[4] = 1e-06 * 9.63e+06*exp(2*tc[0]-2012.866663218878557*invT);
        k_f_save[5] = 1e-06 * 5.7e+13;
        k_f_save[6] = 1e-06 * 8e+13;
        k_f_save[7] = 1e-06 * 1.5e+13;
        k_f_save[8] = 1e-06 * 1.5e+13;
        k_f_save[9] = 1e-06 * 8.43e+13;
        k_f_save[10] = 1e-06 * 1.02e+09*exp(1.5*tc[0]-4327.6633259205891591*invT);
        k_f_save[11] = 1e-12 * 6.02e+14*exp(-1509.6499974141590883*invT);
        k_f_save[12] = 1e-06 * 3e+13;
        k_f_save[13] = 1e-06 * 3e+13;
        k_f_save[14] = 1e-06 * 3.9e+13*exp(-1781.3869969487077469*invT);
        k_f_save[15] = 1e-06 * 1e+13;
        k_f_save[16] = 1e-06 * 1e+13;
        k_f_save[17] = 1e-06 * 388000*exp(2.5*tc[0]-1559.9716639946311716*invT);
        k_f_save[18] = 1e-06 * 130000*exp(2.5*tc[0]-2516.0833290235987079*invT);
        k_f_save[19] = 1e-06 * 5e+13;
        k_f_save[20] = 1e-06 * 1.02e+07*exp(2*tc[0]-956.11166502896742259*invT);
        k_f_save[21] = 1e-06 * 4.6e+19*exp(-1.41*tc[0]-14568.122475046635373*invT);
        k_f_save[22] = 1e-06 * 1.02e+07*exp(2*tc[0]-956.11166502896742259*invT);
        k_f_save[23] = 1e-06 * 3e+13;
        k_f_save[24] = 1e-06 * 1.92e+07*exp(1.83*tc[0]-110.7076664770383303*invT);
        k_f_save[25] = 1e-06 * 1.32e+14;
        k_f_save[26] = 1e-06 * 8.98e+07*exp(1.92*tc[0]-2863.3028284288548093*invT);
        k_f_save[27] = 1e-06 * 1e+14;
        k_f_save[28] = 1e-06 * 1e+13*exp(-4025.733326437757114*invT);
        k_f_save[29] = 1e-06 * 1.75e+12*exp(-679.34249883637164658*invT);
        k_f_save[30] = 1e-06 * 2.5e+12*exp(-24053.756625465601246*invT);
        k_f_save[31] = 1e-06 * 1e+14*exp(-20128.666632188789663*invT);
        k_f_save[32] = 1e-12 * 2.8e+18*exp(-0.86*tc[0]);
        k_f_save[33] = 1e-12 * 3e+20*exp(-1.72*tc[0]);
        k_f_save[34] = 1e-12 * 9.38e+18*exp(-0.76*tc[0]);
        k_f_save[35] = 1e-12 * 3.75e+20*exp(-1.72*tc[0]);
        k_f_save[36] = 1e-12 * 7e+17*exp(-0.8*tc[0]);
        k_f_save[37] = 1e-06 * 8.3e+13*exp(-7252.8618042434254676*invT);
        k_f_save[38] = 1e-12 * 1e+18*exp(-1*tc[0]);
        k_f_save[39] = 1e-12 * 9e+16*exp(-0.6*tc[0]);
        k_f_save[40] = 1e-12 * 6e+19*exp(-1.25*tc[0]);
        k_f_save[41] = 1e-12 * 5.5e+20*exp(-2*tc[0]);
        k_f_save[42] = 1e-12 * 2.2e+22*exp(-2*tc[0]);
        k_f_save[43] = 1e-06 * 3.97e+12*exp(-337.65838275496690812*invT);
        k_f_save[44] = 1e-06 * 2.8e+13*exp(-537.43539907944057177*invT);
        k_f_save[45] = 1e-06 * 1.34e+14*exp(-319.54258278599701271*invT);
        k_f_save[46] = 1e-06 * 1.21e+07*exp(2*tc[0]-2616.7266621845424197*invT);
        k_f_save[47] = 1e-06 * 1e+13*exp(-1811.579996896990906*invT);
        k_f_save[48] = 1e-06 * 1.1e+14;
        k_f_save[49] = 1e-06 * 2.5e+16*exp(-0.8*tc[0]);
        k_f_save[50] = 1e-06 * 3e+13;
        k_f_save[51] = 1e-06 * 1.27e+16*exp(-0.63*tc[0]-192.73198300320765952*invT);
        k_f_save[52] = 1e-06 * 6.6e+08*exp(1.62*tc[0]-5454.8686573231616421*invT);
        k_f_save[53] = 1e-06 * 1.09e+12*exp(0.48*tc[0]+130.83633310922709825*invT);
        k_f_save[54] = 1e-06 * 7.34e+13;
        k_f_save[55] = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1811.579996896990906*invT);
        k_f_save[56] = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1308.3633310922712099*invT);
        k_f_save[57] = 1e-06 * 2.3e+10*exp(1.05*tc[0]-1648.0345805104568626*invT);
        k_f_save[58] = 1e-06 * 1.8e+13;
        k_f_save[59] = 1e-06 * 2e+13;
        k_f_save[60] = 1e-06 * 1.2e+13;
        k_f_save[61] = 1e-06 * 6e+12;
        k_f_save[62] = 1e-06 * 5e+13;
        k_f_save[63] = 1e-06 * 3.4e+06*exp(1.6*tc[0]);
        k_f_save[64] = 1e-06 * 2e+13;
        k_f_save[65] = 1e-06 * 3.2e+13;
        k_f_save[66] = 1e-06 * 1.6e+13;
        k_f_save[67] = 1e-06 * 1.7e+07*exp(2.1*tc[0]-2450.6651624689848177*invT);
        k_f_save[68] = 1e-06 * 4.2e+06*exp(2.1*tc[0]-2450.6651624689848177*invT);
        k_f_save[69] = 1e-06 * 1e+17*exp(-1*tc[0]);
        k_f_save[70] = 1e-06 * 5.6e+12*exp(-1207.7199979313272706*invT);
        k_f_save[71] = 1e-06 * 6.08e+12*exp(0.27*tc[0]-140.90066642532150354*invT);
        k_f_save[72] = 1e-06 * 3e+13;
        k_f_save[73] = 1e-06 * 1.08e+12*exp(0.454*tc[0]-915.85433176458980142*invT);
        k_f_save[74] = 1e-06 * 1.325e+06*exp(2.53*tc[0]-6159.3719894497689893*invT);
        k_f_save[75] = 1e-06 * 5.21e+17*exp(-0.99*tc[0]-795.08233197145705162*invT);
        k_f_save[76] = 1e-06 * 2e+12;
        k_f_save[77] = 1e-06 * 1.15e+08*exp(1.9*tc[0]-3789.2214935095394139*invT);
        k_f_save[78] = 1e-06 * 1e+14;
        k_f_save[79] = 1e-06 * 5e+13*exp(-4025.733326437757114*invT);
        k_f_save[80] = 1e-06 * 1.13e+13*exp(-1725.0267303785790318*invT);
        k_f_save[81] = 1e-06 * 1e+13;
        k_f_save[82] = 1e-06 * 4.3e+07*exp(1.5*tc[0]-40056.046598055690993*invT);
        k_f_save[83] = 1e-06 * 2.16e+08*exp(1.51*tc[0]-1726.0331637101885462*invT);
        k_f_save[84] = 1e-06 * 7.4e+13*exp(-0.37*tc[0]);
        k_f_save[85] = 1e-06 * 35700*exp(2.4*tc[0]+1061.7871648479585929*invT);
        k_f_save[86] = 1e-06 * 2.9e+13*exp(+251.60833290235981963*invT);
        k_f_save[87] = 1e-06 * 1.75e+12*exp(-161.02933305751031412*invT);
        k_f_save[88] = 1e-06 * 5.8e+14*exp(-4810.751325093120613*invT);
        k_f_save[89] = 1e-06 * 5e+13;
        k_f_save[90] = 1e-06 * 3e+13;
        k_f_save[91] = 1e-06 * 2e+13;
        k_f_save[92] = 1e-06 * 1.13e+07*exp(2*tc[0]-1509.6499974141590883*invT);
        k_f_save[93] = 1e-06 * 3e+13;
        k_f_save[94] = 1e-06 * 6.3e+13;
        k_f_save[95] = 1e-06 * 5.6e+07*exp(1.6*tc[0]-2727.4343286615808211*invT);
        k_f_save[96] = 1e-06 * 2.501e+13;
        k_f_save[97] = 1e-06 * 1e+08*exp(1.6*tc[0]-1570.0359973107254064*invT);
        k_f_save[98] = 1e-06 * 4.76e+07*exp(1.228*tc[0]-35.225166606330375885*invT);
        k_f_save[99] = 1e-06 * 5e+13;
        k_f_save[100] = 1e-06 * 3.43e+09*exp(1.18*tc[0]+224.93784961470970529*invT);
        k_f_save[101] = 1e-06 * 5e+12;
        k_f_save[102] = 1e-06 * 5e+12;
        k_f_save[103] = 1e-06 * 1.44e+06*exp(2*tc[0]+422.70199927596451062*invT);
        k_f_save[104] = 1e-06 * 6.3e+06*exp(2*tc[0]-754.82499870707954415*invT);
        k_f_save[105] = 1e-06 * 2e+13;
        k_f_save[106] = 1e-06 * 0.000218*exp(4.5*tc[0]+503.21666580471963925*invT);
        k_f_save[107] = 1e-06 * 504000*exp(2.3*tc[0]-6793.4249883637157836*invT);
        k_f_save[108] = 1e-06 * 3.37e+07*exp(2*tc[0]-7045.0333212660752906*invT);
        k_f_save[109] = 1e-06 * 0.000483*exp(4*tc[0]+1006.4333316094392785*invT);
        k_f_save[110] = 1e-06 * 5e+12;
        k_f_save[111] = 1e-06 * 3.6e+06*exp(2*tc[0]-1258.0416645117993539*invT);
        k_f_save[112] = 1e-06 * 3.54e+06*exp(2.12*tc[0]-437.79849925010609013*invT);
        k_f_save[113] = 1e-06 * 7.5e+12*exp(-1006.4333316094392785*invT);
        k_f_save[114] = 1e-06 * 1.3e+11*exp(+820.24316526169309327*invT);
        k_f_save[115] = 1e-06 * 4.2e+14*exp(-6038.5999896566363532*invT);
        k_f_save[116] = 1e-06 * 2e+13;
        k_f_save[117] = 1e-06 * 1e+12;
        k_f_save[118] = 1e-06 * 2e+13;
        k_f_save[119] = 1e-06 * 1.5e+14*exp(-11875.913312991384373*invT);
        k_f_save[120] = 1e-06 * 1e+12*exp(-4025.733326437757114*invT);
        k_f_save[121] = 1e-06 * 5.8e+13*exp(-289.85279950351855405*invT);
        k_f_save[122] = 1e-06 * 5e+13;
        k_f_save[123] = 1e-06 * 5e+13;
        k_f_save[124] = 1e-06 * 3.3e+13;
        k_f_save[125] = 1e-06 * 1.107e+08*exp(1.79*tc[0]-840.37183189388190385*invT);
        k_f_save[126] = 1e-06 * 5.71e+12*exp(+379.92858268256338761*invT);
        k_f_save[127] = 1e-06 * 4e+13;
        k_f_save[128] = 1e-06 * 3e+13;
        k_f_save[129] = 1e-06 * 6e+13;
        k_f_save[130] = 1e-06 * 5e+13;
        k_f_save[131] = 1e-06 * 3.4e+12*exp(-347.2194994052565562*invT);
        k_f_save[132] = 1e-06 * 9.46e+13*exp(+259.15658288943063781*invT);
        k_f_save[133] = 1e-06 * 5e+13;
        k_f_save[134] = 1e-06 * 1.32e+13*exp(-754.82499870707954415*invT);
        k_f_save[135] = 1e-06 * 500000*exp(2*tc[0]-3638.256493768123164*invT);
        k_f_save[136] = 1e-06 * 3.2e+13;
        k_f_save[137] = 1e-06 * 4e+13;
        k_f_save[138] = 1e-06 * 2.46e+06*exp(2*tc[0]-4161.6018262050320118*invT);
        k_f_save[139] = 1e-06 * 8.1e+11*exp(0.5*tc[0]-2269.5071627792858635*invT);
        k_f_save[140] = 1e-06 * 3e+13;
        k_f_save[141] = 1e-06 * 1.5e+13*exp(-301.92999948283181766*invT);
        k_f_save[142] = 1e-06 * 9e+12*exp(-301.92999948283181766*invT);
        k_f_save[143] = 1e-06 * 2.8e+13;
        k_f_save[144] = 1e-06 * 1.2e+13;
        k_f_save[145] = 1e-06 * 7e+13;
        k_f_save[146] = 1e-06 * 2e+13;
        k_f_save[147] = 1e-06 * 3e+13;
        k_f_save[148] = 1e-06 * 1.2e+13*exp(+286.83349950869023814*invT);
        k_f_save[149] = 1e-06 * 1.6e+13*exp(+286.83349950869023814*invT);
        k_f_save[150] = 1e-06 * 9e+12;
        k_f_save[151] = 1e-06 * 7e+12;
        k_f_save[152] = 1e-06 * 1.4e+13;
        k_f_save[153] = 1e-06 * 4e+13*exp(+276.76916619259583285*invT);
        k_f_save[154] = 1e-06 * 2.675e+13*exp(-14492.639975175927248*invT);
        k_f_save[155] = 1e-06 * 3.6e+10*exp(-4498.7569922941938785*invT);
        k_f_save[156] = 1e-06 * 24500*exp(2.47*tc[0]-2606.6623288684481849*invT);
        k_f_save[157] = 1e-06 * 2.12e+16*exp(-0.97*tc[0]-311.99433279892622295*invT);
        k_f_save[158] = 1e-06 * 4.99e+12*exp(0.1*tc[0]-5334.096657530029006*invT);
        k_f_save[159] = 1e-06 * 2.648e+13;
        k_f_save[160] = 1e-06 * 3320*exp(2.81*tc[0]-2948.8496616156576238*invT);
        k_f_save[161] = 1e-06 * 3e+07*exp(1.5*tc[0]-5001.973658098913802*invT);
        k_f_save[162] = 1e-06 * 1e+07*exp(1.5*tc[0]-5001.973658098913802*invT);
        k_f_save[163] = 1e-06 * 227000*exp(2*tc[0]-4629.5933254034207494*invT);
        k_f_save[164] = 1e-06 * 6.14e+06*exp(1.74*tc[0]-5258.6141576593208811*invT);
        k_f_save[165] = 1e-06 * 2.244e+18*exp(-1*tc[0]-8554.6833186802341515*invT);
        k_f_save[166] = 1e-06 * 1.87e+17*exp(-1*tc[0]-8554.6833186802341515*invT);
        k_f_save[167] = 1e-06 * 7.6e+12*exp(-201.28666632188787844*invT);
        k_f_save[168] = 1e-06 * 1.8e+13*exp(-452.89499922424772649*invT);
        k_f_save[169] = 1e-06 * 4.28e-13*exp(7.6*tc[0]+1776.3548302906606295*invT);
        k_f_save[170] = 1e-06 * 5e+13*exp(-754.82499870707954415*invT);
        k_f_save[171] = 1e-06 * 407000*exp(2.4*tc[0]-100.64333316094393922*invT);
        k_f_save[172] = 1e-06 * 3.98e+12*exp(+120.77199979313272138*invT);
        k_f_save[173] = 1 * 8e+12*exp(0.44*tc[0]-44670.543423484967207*invT);
        k_f_save[174] = 1e-06 * 8.4e+11*exp(-1949.9645799932889076*invT);
        k_f_save[175] = 1e-06 * 1.6e+12*exp(-429.7470325972306*invT);
        k_f_save[176] = 1e-06 * 1e+13;

        Kc_save[0] = 1.0 / (refC) * exp((2 * g_RT[2]) - (g_RT[3]));
        Kc_save[1] = 1.0 / (refC) * exp((g_RT[2] + g_RT[1]) - (g_RT[4]));
        Kc_save[2] = exp((g_RT[2] + g_RT[0]) - (g_RT[1] + g_RT[4]));
        Kc_save[3] = exp((g_RT[2] + g_RT[6]) - (g_RT[4] + g_RT[3]));
        Kc_save[4] = exp((g_RT[2] + g_RT[7]) - (g_RT[4] + g_RT[6]));
        Kc_save[5] = exp((g_RT[2] + g_RT[9]) - (g_RT[1] + g_RT[14]));
        Kc_save[6] = exp((g_RT[2] + g_RT[10]) - (g_RT[1] + g_RT[16]));
        Kc_save[7] = exp((g_RT[2] + g_RT[11]) - (g_RT[0] + g_RT[14]));
        Kc_save[8] = exp((g_RT[2] + g_RT[11]) - (g_RT[1] + g_RT[16]));
        Kc_save[9] = exp((g_RT[2] + g_RT[12]) - (g_RT[1] + g_RT[17]));
        Kc_save[10] = exp((g_RT[2] + g_RT[13]) - (g_RT[4] + g_RT[12]));
        Kc_save[11] = 1.0 / (refC) * exp((g_RT[2] + g_RT[14]) - (g_RT[15]));
        Kc_save[12] = exp((g_RT[2] + g_RT[16]) - (g_RT[4] + g_RT[14]));
        Kc_save[13] = exp((g_RT[2] + g_RT[16]) - (g_RT[1] + g_RT[15]));
        Kc_save[14] = exp((g_RT[2] + g_RT[17]) - (g_RT[4] + g_RT[16]));
        Kc_save[15] = exp((g_RT[2] + g_RT[18]) - (g_RT[4] + g_RT[17]));
        Kc_save[16] = exp((g_RT[2] + g_RT[19]) - (g_RT[4] + g_RT[17]));
        Kc_save[17] = exp((g_RT[2] + g_RT[20]) - (g_RT[4] + g_RT[18]));
        Kc_save[18] = exp((g_RT[2] + g_RT[20]) - (g_RT[4] + g_RT[19]));
        Kc_save[19] = exp((g_RT[2] + g_RT[21]) - (g_RT[9] + g_RT[14]));
        Kc_save[20] = exp((g_RT[2] + g_RT[22]) - (g_RT[1] + g_RT[27]));
        Kc_save[21] = exp((g_RT[2] + g_RT[22]) - (g_RT[4] + g_RT[21]));
        Kc_save[22] = exp((g_RT[2] + g_RT[22]) - (g_RT[14] + g_RT[10]));
        Kc_save[23] = exp((g_RT[2] + g_RT[23]) - (g_RT[1] + g_RT[28]));
        Kc_save[24] = exp((g_RT[2] + g_RT[24]) - (g_RT[12] + g_RT[16]));
        Kc_save[25] = exp((g_RT[2] + g_RT[25]) - (g_RT[12] + g_RT[17]));
        Kc_save[26] = exp((g_RT[2] + g_RT[26]) - (g_RT[4] + g_RT[25]));
        Kc_save[27] = refC * exp((g_RT[2] + g_RT[27]) - (g_RT[1] + 2 * g_RT[14]));
        Kc_save[28] = exp((g_RT[2] + g_RT[28]) - (g_RT[4] + g_RT[27]));
        Kc_save[29] = exp((g_RT[2] + g_RT[28]) - (g_RT[10] + g_RT[15]));
        Kc_save[30] = exp((g_RT[3] + g_RT[14]) - (g_RT[2] + g_RT[15]));
        Kc_save[31] = exp((g_RT[3] + g_RT[17]) - (g_RT[6] + g_RT[16]));
        Kc_save[32] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));
        Kc_save[33] = 1.0 / (refC) * exp((g_RT[1] + 2 * g_RT[3]) - (g_RT[6] + g_RT[3]));
        Kc_save[34] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));
        Kc_save[35] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[30]) - (g_RT[6] + g_RT[30]));
        Kc_save[36] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[31]) - (g_RT[6] + g_RT[31]));
        Kc_save[37] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));
        Kc_save[38] = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));
        Kc_save[39] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));
        Kc_save[40] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));
        Kc_save[41] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[15]) - (g_RT[0] + g_RT[15]));
        Kc_save[42] = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));
        Kc_save[43] = exp((g_RT[1] + g_RT[6]) - (g_RT[2] + g_RT[5]));
        Kc_save[44] = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));
        Kc_save[45] = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));
        Kc_save[46] = exp((g_RT[1] + g_RT[7]) - (g_RT[6] + g_RT[0]));
        Kc_save[47] = exp((g_RT[1] + g_RT[7]) - (g_RT[4] + g_RT[5]));
        Kc_save[48] = exp((g_RT[1] + g_RT[9]) - (g_RT[8] + g_RT[0]));
        Kc_save[49] = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[12]));
        Kc_save[50] = exp((g_RT[1] + g_RT[11]) - (g_RT[9] + g_RT[0]));
        Kc_save[51] = 1.0 / (refC) * exp((g_RT[1] + g_RT[12]) - (g_RT[13]));
        Kc_save[52] = exp((g_RT[1] + g_RT[13]) - (g_RT[12] + g_RT[0]));
        Kc_save[53] = 1.0 / (refC) * exp((g_RT[1] + g_RT[16]) - (g_RT[17]));
        Kc_save[54] = exp((g_RT[1] + g_RT[16]) - (g_RT[0] + g_RT[14]));
        Kc_save[55] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[18]));
        Kc_save[56] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[19]));
        Kc_save[57] = exp((g_RT[1] + g_RT[17]) - (g_RT[16] + g_RT[0]));
        Kc_save[58] = 1.0 / (refC) * exp((g_RT[1] + g_RT[18]) - (g_RT[20]));
        Kc_save[59] = exp((g_RT[1] + g_RT[18]) - (g_RT[0] + g_RT[17]));
        Kc_save[60] = exp((g_RT[1] + g_RT[18]) - (g_RT[4] + g_RT[12]));
        Kc_save[61] = exp((g_RT[1] + g_RT[18]) - (g_RT[11] + g_RT[5]));
        Kc_save[62] = 1.0 / (refC) * exp((g_RT[1] + g_RT[19]) - (g_RT[20]));
        Kc_save[63] = exp((g_RT[1] + g_RT[19]) - (g_RT[1] + g_RT[18]));
        Kc_save[64] = exp((g_RT[1] + g_RT[19]) - (g_RT[0] + g_RT[17]));
        Kc_save[65] = exp((g_RT[1] + g_RT[19]) - (g_RT[4] + g_RT[12]));
        Kc_save[66] = exp((g_RT[1] + g_RT[19]) - (g_RT[11] + g_RT[5]));
        Kc_save[67] = exp((g_RT[1] + g_RT[20]) - (g_RT[18] + g_RT[0]));
        Kc_save[68] = exp((g_RT[1] + g_RT[20]) - (g_RT[19] + g_RT[0]));
        Kc_save[69] = 1.0 / (refC) * exp((g_RT[1] + g_RT[21]) - (g_RT[22]));
        Kc_save[70] = 1.0 / (refC) * exp((g_RT[1] + g_RT[22]) - (g_RT[23]));
        Kc_save[71] = 1.0 / (refC) * exp((g_RT[1] + g_RT[23]) - (g_RT[24]));
        Kc_save[72] = exp((g_RT[1] + g_RT[23]) - (g_RT[0] + g_RT[22]));
        Kc_save[73] = 1.0 / (refC) * exp((g_RT[1] + g_RT[24]) - (g_RT[25]));
        Kc_save[74] = exp((g_RT[1] + g_RT[24]) - (g_RT[23] + g_RT[0]));
        Kc_save[75] = 1.0 / (refC) * exp((g_RT[1] + g_RT[25]) - (g_RT[26]));
        Kc_save[76] = exp((g_RT[1] + g_RT[25]) - (g_RT[0] + g_RT[24]));
        Kc_save[77] = exp((g_RT[1] + g_RT[26]) - (g_RT[25] + g_RT[0]));
        Kc_save[78] = exp((g_RT[1] + g_RT[27]) - (g_RT[11] + g_RT[14]));
        Kc_save[79] = exp((g_RT[1] + g_RT[28]) - (g_RT[27] + g_RT[0]));
        Kc_save[80] = exp((g_RT[1] + g_RT[28]) - (g_RT[12] + g_RT[14]));
        Kc_save[81] = exp((g_RT[1] + g_RT[29]) - (g_RT[1] + g_RT[28]));
        Kc_save[82] = 1.0 / (refC) * exp((g_RT[0] + g_RT[14]) - (g_RT[17]));
        Kc_save[83] = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));
        Kc_save[84] = 1.0 / (refC) * exp((2 * g_RT[4]) - (g_RT[7]));
        Kc_save[85] = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));
        Kc_save[86] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
        Kc_save[87] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
        Kc_save[88] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
        Kc_save[89] = exp((g_RT[4] + g_RT[8]) - (g_RT[1] + g_RT[14]));
        Kc_save[90] = exp((g_RT[4] + g_RT[9]) - (g_RT[1] + g_RT[16]));
        Kc_save[91] = exp((g_RT[4] + g_RT[10]) - (g_RT[1] + g_RT[17]));
        Kc_save[92] = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[5]));
        Kc_save[93] = exp((g_RT[4] + g_RT[11]) - (g_RT[1] + g_RT[17]));
        Kc_save[94] = 1.0 / (refC) * exp((g_RT[4] + g_RT[12]) - (g_RT[20]));
        Kc_save[95] = exp((g_RT[4] + g_RT[12]) - (g_RT[10] + g_RT[5]));
        Kc_save[96] = exp((g_RT[4] + g_RT[12]) - (g_RT[11] + g_RT[5]));
        Kc_save[97] = exp((g_RT[4] + g_RT[13]) - (g_RT[12] + g_RT[5]));
        Kc_save[98] = exp((g_RT[4] + g_RT[14]) - (g_RT[1] + g_RT[15]));
        Kc_save[99] = exp((g_RT[4] + g_RT[16]) - (g_RT[5] + g_RT[14]));
        Kc_save[100] = exp((g_RT[4] + g_RT[17]) - (g_RT[16] + g_RT[5]));
        Kc_save[101] = exp((g_RT[4] + g_RT[18]) - (g_RT[5] + g_RT[17]));
        Kc_save[102] = exp((g_RT[4] + g_RT[19]) - (g_RT[5] + g_RT[17]));
        Kc_save[103] = exp((g_RT[4] + g_RT[20]) - (g_RT[18] + g_RT[5]));
        Kc_save[104] = exp((g_RT[4] + g_RT[20]) - (g_RT[19] + g_RT[5]));
        Kc_save[105] = exp((g_RT[4] + g_RT[21]) - (g_RT[1] + g_RT[27]));
        Kc_save[106] = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[28]));
        Kc_save[107] = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[29]));
        Kc_save[108] = exp((g_RT[4] + g_RT[22]) - (g_RT[21] + g_RT[5]));
        Kc_save[109] = exp((g_RT[4] + g_RT[22]) - (g_RT[12] + g_RT[14]));
        Kc_save[110] = exp((g_RT[4] + g_RT[23]) - (g_RT[5] + g_RT[22]));
        Kc_save[111] = exp((g_RT[4] + g_RT[24]) - (g_RT[23] + g_RT[5]));
        Kc_save[112] = exp((g_RT[4] + g_RT[26]) - (g_RT[25] + g_RT[5]));
        Kc_save[113] = exp((g_RT[4] + g_RT[28]) - (g_RT[27] + g_RT[5]));
        Kc_save[114] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
        Kc_save[115] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
        Kc_save[116] = exp((g_RT[6] + g_RT[10]) - (g_RT[4] + g_RT[17]));
        Kc_save[117] = exp((g_RT[6] + g_RT[12]) - (g_RT[3] + g_RT[13]));
        Kc_save[118] = exp((g_RT[6] + g_RT[12]) - (g_RT[4] + g_RT[19]));
        Kc_save[119] = exp((g_RT[6] + g_RT[14]) - (g_RT[4] + g_RT[15]));
        Kc_save[120] = exp((g_RT[6] + g_RT[17]) - (g_RT[16] + g_RT[7]));
        Kc_save[121] = exp((g_RT[8] + g_RT[3]) - (g_RT[2] + g_RT[14]));
        Kc_save[122] = exp((g_RT[8] + g_RT[10]) - (g_RT[1] + g_RT[21]));
        Kc_save[123] = exp((g_RT[8] + g_RT[12]) - (g_RT[1] + g_RT[22]));
        Kc_save[124] = exp((g_RT[9] + g_RT[3]) - (g_RT[2] + g_RT[16]));
        Kc_save[125] = exp((g_RT[9] + g_RT[0]) - (g_RT[1] + g_RT[10]));
        Kc_save[126] = exp((g_RT[9] + g_RT[5]) - (g_RT[1] + g_RT[17]));
        Kc_save[127] = exp((g_RT[9] + g_RT[10]) - (g_RT[1] + g_RT[22]));
        Kc_save[128] = exp((g_RT[9] + g_RT[12]) - (g_RT[1] + g_RT[23]));
        Kc_save[129] = exp((g_RT[9] + g_RT[13]) - (g_RT[1] + g_RT[24]));
        Kc_save[130] = 1.0 / (refC) * exp((g_RT[9] + g_RT[14]) - (g_RT[27]));
        Kc_save[131] = exp((g_RT[9] + g_RT[15]) - (g_RT[16] + g_RT[14]));
        Kc_save[132] = exp((g_RT[9] + g_RT[17]) - (g_RT[1] + g_RT[28]));
        Kc_save[133] = exp((g_RT[9] + g_RT[27]) - (g_RT[14] + g_RT[22]));
        Kc_save[134] = exp((g_RT[10] + g_RT[3]) - (g_RT[4] + g_RT[16]));
        Kc_save[135] = exp((g_RT[10] + g_RT[0]) - (g_RT[1] + g_RT[12]));
        Kc_save[136] = exp((2 * g_RT[10]) - (g_RT[0] + g_RT[22]));
        Kc_save[137] = exp((g_RT[10] + g_RT[12]) - (g_RT[1] + g_RT[24]));
        Kc_save[138] = exp((g_RT[10] + g_RT[13]) - (2 * g_RT[12]));
        Kc_save[139] = 1.0 / (refC) * exp((g_RT[10] + g_RT[14]) - (g_RT[28]));
        Kc_save[140] = exp((g_RT[10] + g_RT[27]) - (g_RT[23] + g_RT[14]));
        Kc_save[141] = exp((g_RT[11] + g_RT[30]) - (g_RT[10] + g_RT[30]));
        Kc_save[142] = exp((g_RT[11] + g_RT[31]) - (g_RT[10] + g_RT[31]));
        Kc_save[143] = refC * exp((g_RT[11] + g_RT[3]) - (g_RT[1] + g_RT[4] + g_RT[14]));
        Kc_save[144] = exp((g_RT[11] + g_RT[3]) - (g_RT[14] + g_RT[5]));
        Kc_save[145] = exp((g_RT[11] + g_RT[0]) - (g_RT[12] + g_RT[1]));
        Kc_save[146] = 1.0 / (refC) * exp((g_RT[11] + g_RT[5]) - (g_RT[20]));
        Kc_save[147] = exp((g_RT[11] + g_RT[5]) - (g_RT[10] + g_RT[5]));
        Kc_save[148] = exp((g_RT[11] + g_RT[12]) - (g_RT[1] + g_RT[24]));
        Kc_save[149] = exp((g_RT[11] + g_RT[13]) - (2 * g_RT[12]));
        Kc_save[150] = exp((g_RT[11] + g_RT[14]) - (g_RT[10] + g_RT[14]));
        Kc_save[151] = exp((g_RT[11] + g_RT[15]) - (g_RT[10] + g_RT[15]));
        Kc_save[152] = exp((g_RT[11] + g_RT[15]) - (g_RT[14] + g_RT[17]));
        Kc_save[153] = exp((g_RT[11] + g_RT[26]) - (g_RT[12] + g_RT[25]));
        Kc_save[154] = exp((g_RT[12] + g_RT[3]) - (g_RT[2] + g_RT[19]));
        Kc_save[155] = exp((g_RT[12] + g_RT[3]) - (g_RT[4] + g_RT[17]));
        Kc_save[156] = exp((g_RT[12] + g_RT[7]) - (g_RT[6] + g_RT[13]));
        Kc_save[157] = 1.0 / (refC) * exp((2 * g_RT[12]) - (g_RT[26]));
        Kc_save[158] = exp((2 * g_RT[12]) - (g_RT[1] + g_RT[25]));
        Kc_save[159] = exp((g_RT[12] + g_RT[16]) - (g_RT[13] + g_RT[14]));
        Kc_save[160] = exp((g_RT[12] + g_RT[17]) - (g_RT[16] + g_RT[13]));
        Kc_save[161] = exp((g_RT[12] + g_RT[20]) - (g_RT[18] + g_RT[13]));
        Kc_save[162] = exp((g_RT[12] + g_RT[20]) - (g_RT[19] + g_RT[13]));
        Kc_save[163] = exp((g_RT[12] + g_RT[24]) - (g_RT[23] + g_RT[13]));
        Kc_save[164] = exp((g_RT[12] + g_RT[26]) - (g_RT[25] + g_RT[13]));
        Kc_save[165] = refC * exp((g_RT[16] + g_RT[5]) - (g_RT[1] + g_RT[14] + g_RT[5]));
        Kc_save[166] = refC * exp((g_RT[16]) - (g_RT[1] + g_RT[14]));
        Kc_save[167] = exp((g_RT[16] + g_RT[3]) - (g_RT[6] + g_RT[14]));
        Kc_save[168] = exp((g_RT[18] + g_RT[3]) - (g_RT[6] + g_RT[17]));
        Kc_save[169] = exp((g_RT[19] + g_RT[3]) - (g_RT[6] + g_RT[17]));
        Kc_save[170] = exp((g_RT[21] + g_RT[3]) - (g_RT[16] + g_RT[14]));
        Kc_save[171] = exp((g_RT[21] + g_RT[0]) - (g_RT[1] + g_RT[22]));
        Kc_save[172] = exp((g_RT[23] + g_RT[3]) - (g_RT[16] + g_RT[17]));
        Kc_save[173] = refC * exp((g_RT[24]) - (g_RT[0] + g_RT[22]));
        Kc_save[174] = exp((g_RT[25] + g_RT[3]) - (g_RT[6] + g_RT[24]));
        Kc_save[175] = refC * exp((g_RT[27] + g_RT[3]) - (g_RT[4] + 2 * g_RT[14]));
        Kc_save[176] = refC * exp((2 * g_RT[27]) - (2 * g_RT[14] + g_RT[22]));
    }

    /*reaction 1: 2 O + M <=> O2 + M */
    phi_f = sc[2]*sc[2];
    alpha = mixture + 1.4*sc[0] + 14.4*sc[5] + sc[13] + 0.75*sc[14] + 2.6*sc[15] + 2*sc[26] + -0.17*sc[31];
    k_f = alpha * k_f_save[0];
    q_f = phi_f * k_f;
    phi_r = sc[3];
    Kc = Kc_save[0];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 2 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 2: O + H + M <=> OH + M */
    phi_f = sc[2]*sc[1];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = alpha * k_f_save[1];
    q_f = phi_f * k_f;
    phi_r = sc[4];
    Kc = Kc_save[1];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 3: O + H2 <=> H + OH */
    phi_f = sc[2]*sc[0];
    k_f = k_f_save[2];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[4];
    Kc = Kc_save[2];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 4: O + HO2 <=> OH + O2 */
    phi_f = sc[2]*sc[6];
    k_f = k_f_save[3];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = Kc_save[3];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 5: O + H2O2 <=> OH + HO2 */
    phi_f = sc[2]*sc[7];
    k_f = k_f_save[4];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[6];
    Kc = Kc_save[4];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 6: O + CH <=> H + CO */
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

    /*reaction 7: O + CH2 <=> H + HCO */
    phi_f = sc[2]*sc[10];
    k_f = k_f_save[6];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_save[6];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 8: O + CH2(S) <=> H2 + CO */
    phi_f = sc[2]*sc[11];
    k_f = k_f_save[7];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[14];
    Kc = Kc_save[7];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 9: O + CH2(S) <=> H + HCO */
    phi_f = sc[2]*sc[11];
    k_f = k_f_save[8];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_save[8];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 10: O + CH3 <=> H + CH2O */
    phi_f = sc[2]*sc[12];
    k_f = k_f_save[9];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[9];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 11: O + CH4 <=> OH + CH3 */
    phi_f = sc[2]*sc[13];
    k_f = k_f_save[10];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = Kc_save[10];
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
    k_f = alpha * k_f_save[11];
    q_f = phi_f * k_f;
    phi_r = sc[15];
    Kc = Kc_save[11];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 13: O + HCO <=> OH + CO */
    phi_f = sc[2]*sc[16];
    k_f = k_f_save[12];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[14];
    Kc = Kc_save[12];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[16] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 14: O + HCO <=> H + CO2 */
    phi_f = sc[2]*sc[16];
    k_f = k_f_save[13];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[15];
    Kc = Kc_save[13];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[16] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 15: O + CH2O <=> OH + HCO */
    phi_f = sc[2]*sc[17];
    k_f = k_f_save[14];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[16];
    Kc = Kc_save[14];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 16: O + CH2OH <=> OH + CH2O */
    phi_f = sc[2]*sc[18];
    k_f = k_f_save[15];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[15];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 17: O + CH3O <=> OH + CH2O */
    phi_f = sc[2]*sc[19];
    k_f = k_f_save[16];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[16];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 18: O + CH3OH <=> OH + CH2OH */
    phi_f = sc[2]*sc[20];
    k_f = k_f_save[17];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[18];
    Kc = Kc_save[17];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[18] += 1 * qdot;

    /*reaction 19: O + CH3OH <=> OH + CH3O */
    phi_f = sc[2]*sc[20];
    k_f = k_f_save[18];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[19];
    Kc = Kc_save[18];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[19] += 1 * qdot;

    /*reaction 20: O + C2H <=> CH + CO */
    phi_f = sc[2]*sc[21];
    k_f = k_f_save[19];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[14];
    Kc = Kc_save[19];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[21] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 21: O + C2H2 <=> H + HCCO */
    phi_f = sc[2]*sc[22];
    k_f = k_f_save[20];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[27];
    Kc = Kc_save[20];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[27] += 1 * qdot;

    /*reaction 22: O + C2H2 <=> OH + C2H */
    phi_f = sc[2]*sc[22];
    k_f = k_f_save[21];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[21];
    Kc = Kc_save[21];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[21] += 1 * qdot;

    /*reaction 23: O + C2H2 <=> CO + CH2 */
    phi_f = sc[2]*sc[22];
    k_f = k_f_save[22];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[10];
    Kc = Kc_save[22];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 24: O + C2H3 <=> H + CH2CO */
    phi_f = sc[2]*sc[23];
    k_f = k_f_save[23];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[23];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[23] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[28] += 1 * qdot;

    /*reaction 25: O + C2H4 <=> CH3 + HCO */
    phi_f = sc[2]*sc[24];
    k_f = k_f_save[24];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[16];
    Kc = Kc_save[24];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[24] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 26: O + C2H5 <=> CH3 + CH2O */
    phi_f = sc[2]*sc[25];
    k_f = k_f_save[25];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[17];
    Kc = Kc_save[25];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[25] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 27: O + C2H6 <=> OH + C2H5 */
    phi_f = sc[2]*sc[26];
    k_f = k_f_save[26];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[25];
    Kc = Kc_save[26];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[26] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[25] += 1 * qdot;

    /*reaction 28: O + HCCO <=> H + 2 CO */
    phi_f = sc[2]*sc[27];
    k_f = k_f_save[27];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14]*sc[14];
    Kc = Kc_save[27];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[27] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[14] += 2 * qdot;

    /*reaction 29: O + CH2CO <=> OH + HCCO */
    phi_f = sc[2]*sc[28];
    k_f = k_f_save[28];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[27];
    Kc = Kc_save[28];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[28] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[27] += 1 * qdot;

    /*reaction 30: O + CH2CO <=> CH2 + CO2 */
    phi_f = sc[2]*sc[28];
    k_f = k_f_save[29];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[15];
    Kc = Kc_save[29];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[28] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 31: O2 + CO <=> O + CO2 */
    phi_f = sc[3]*sc[14];
    k_f = k_f_save[30];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[15];
    Kc = Kc_save[30];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 32: O2 + CH2O <=> HO2 + HCO */
    phi_f = sc[3]*sc[17];
    k_f = k_f_save[31];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[16];
    Kc = Kc_save[31];
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
    k_f = alpha * k_f_save[32];
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = Kc_save[32];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 34: H + 2 O2 <=> HO2 + O2 */
    phi_f = sc[1]*sc[3]*sc[3];
    k_f = k_f_save[33];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = Kc_save[33];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 2 * qdot;
    wdot[6] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 35: H + O2 + H2O <=> HO2 + H2O */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = k_f_save[34];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = Kc_save[34];
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
    k_f = k_f_save[35];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[30];
    Kc = Kc_save[35];
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
    k_f = k_f_save[36];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[31];
    Kc = Kc_save[36];
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
    k_f = k_f_save[37];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = Kc_save[37];
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
    k_f = alpha * k_f_save[38];
    q_f = phi_f * k_f;
    phi_r = sc[0];
    Kc = Kc_save[38];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    phi_f = sc[1]*sc[1]*sc[0];
    k_f = k_f_save[39];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[0];
    Kc = Kc_save[39];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[0] += 2 * qdot;

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    phi_f = sc[1]*sc[1]*sc[5];
    k_f = k_f_save[40];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[5];
    Kc = Kc_save[40];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    phi_f = sc[1]*sc[1]*sc[15];
    k_f = k_f_save[41];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[15];
    Kc = Kc_save[41];
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
    k_f = alpha * k_f_save[42];
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = Kc_save[42];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 44: H + HO2 <=> O + H2O */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[43];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = Kc_save[43];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 45: H + HO2 <=> O2 + H2 */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[44];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = Kc_save[44];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 46: H + HO2 <=> 2 OH */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[45];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = Kc_save[45];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[4] += 2 * qdot;

    /*reaction 47: H + H2O2 <=> HO2 + H2 */
    phi_f = sc[1]*sc[7];
    k_f = k_f_save[46];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = Kc_save[46];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 48: H + H2O2 <=> OH + H2O */
    phi_f = sc[1]*sc[7];
    k_f = k_f_save[47];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[5];
    Kc = Kc_save[47];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 49: H + CH <=> C + H2 */
    phi_f = sc[1]*sc[9];
    k_f = k_f_save[48];
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = Kc_save[48];
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
    k_f = k_f_save[49];
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
    phi_r = sc[12];
    Kc = Kc_save[49];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 51: H + CH2(S) <=> CH + H2 */
    phi_f = sc[1]*sc[11];
    k_f = k_f_save[50];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = Kc_save[50];
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
    k_f = k_f_save[51];
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
    phi_r = sc[13];
    Kc = Kc_save[51];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 53: H + CH4 <=> CH3 + H2 */
    phi_f = sc[1]*sc[13];
    k_f = k_f_save[52];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[0];
    Kc = Kc_save[52];
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
    k_f = k_f_save[53];
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
    phi_r = sc[17];
    Kc = Kc_save[53];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[16] -= 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 55: H + HCO <=> H2 + CO */
    phi_f = sc[1]*sc[16];
    k_f = k_f_save[54];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[14];
    Kc = Kc_save[54];
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
    k_f = k_f_save[55];
    redP = 1e-12 * alpha / k_f * 1.27e+32*exp(-4.82*tc[0]-3286.0048277048194905*invT);
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
    Kc = Kc_save[55];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[18] += 1 * qdot;

    /*reaction 57: H + CH2O (+M) <=> CH3O (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[56];
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
    phi_r = sc[19];
    Kc = Kc_save[56];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[19] += 1 * qdot;

    /*reaction 58: H + CH2O <=> HCO + H2 */
    phi_f = sc[1]*sc[17];
    k_f = k_f_save[57];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[0];
    Kc = Kc_save[57];
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
    k_f = k_f_save[58];
    redP = 1e-12 * alpha / k_f * 3e+31*exp(-4.8*tc[0]-1660.6149971555748834*invT);
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
    Kc = Kc_save[58];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[20] += 1 * qdot;

    /*reaction 60: H + CH2OH <=> H2 + CH2O */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[59];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[17];
    Kc = Kc_save[59];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 61: H + CH2OH <=> OH + CH3 */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[60];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = Kc_save[60];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 62: H + CH2OH <=> CH2(S) + H2O */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[61];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[61];
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
    k_f = k_f_save[62];
    redP = 1e-12 * alpha / k_f * 8.6e+28*exp(-4*tc[0]-1522.2304140592771091*invT);
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
    Kc = Kc_save[62];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[20] += 1 * qdot;

    /*reaction 64: H + CH3O <=> H + CH2OH */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[63];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[18];
    Kc = Kc_save[63];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[18] += 1 * qdot;

    /*reaction 65: H + CH3O <=> H2 + CH2O */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[64];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[17];
    Kc = Kc_save[64];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 66: H + CH3O <=> OH + CH3 */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[65];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = Kc_save[65];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 67: H + CH3O <=> CH2(S) + H2O */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[66];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[66];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 68: H + CH3OH <=> CH2OH + H2 */
    phi_f = sc[1]*sc[20];
    k_f = k_f_save[67];
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[0];
    Kc = Kc_save[67];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 69: H + CH3OH <=> CH3O + H2 */
    phi_f = sc[1]*sc[20];
    k_f = k_f_save[68];
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[0];
    Kc = Kc_save[68];
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
    k_f = k_f_save[69];
    redP = 1e-12 * alpha / k_f * 3.75e+33*exp(-4.8*tc[0]-956.11166502896742259*invT);
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
    Kc = Kc_save[69];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[21] -= 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 71: H + C2H2 (+M) <=> C2H3 (+M) */
    phi_f = sc[1]*sc[22];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[70];
    redP = 1e-12 * alpha / k_f * 3.8e+40*exp(-7.27*tc[0]-3633.2243271100760467*invT);
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
    Kc = Kc_save[70];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[23] += 1 * qdot;

    /*reaction 72: H + C2H3 (+M) <=> C2H4 (+M) */
    phi_f = sc[1]*sc[23];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[71];
    redP = 1e-12 * alpha / k_f * 1.4e+30*exp(-3.86*tc[0]-1670.6793304716693456*invT);
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
    Kc = Kc_save[71];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[23] -= 1 * qdot;
    wdot[24] += 1 * qdot;

    /*reaction 73: H + C2H3 <=> H2 + C2H2 */
    phi_f = sc[1]*sc[23];
    k_f = k_f_save[72];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = Kc_save[72];
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
    k_f = k_f_save[73];
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
    phi_r = sc[25];
    Kc = Kc_save[73];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[24] -= 1 * qdot;
    wdot[25] += 1 * qdot;

    /*reaction 75: H + C2H4 <=> C2H3 + H2 */
    phi_f = sc[1]*sc[24];
    k_f = k_f_save[74];
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[0];
    Kc = Kc_save[74];
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
    k_f = k_f_save[75];
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
    phi_r = sc[26];
    Kc = Kc_save[75];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[25] -= 1 * qdot;
    wdot[26] += 1 * qdot;

    /*reaction 77: H + C2H5 <=> H2 + C2H4 */
    phi_f = sc[1]*sc[25];
    k_f = k_f_save[76];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[24];
    Kc = Kc_save[76];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[25] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[24] += 1 * qdot;

    /*reaction 78: H + C2H6 <=> C2H5 + H2 */
    phi_f = sc[1]*sc[26];
    k_f = k_f_save[77];
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[0];
    Kc = Kc_save[77];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[26] -= 1 * qdot;
    wdot[25] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 79: H + HCCO <=> CH2(S) + CO */
    phi_f = sc[1]*sc[27];
    k_f = k_f_save[78];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[14];
    Kc = Kc_save[78];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[27] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 80: H + CH2CO <=> HCCO + H2 */
    phi_f = sc[1]*sc[28];
    k_f = k_f_save[79];
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[0];
    Kc = Kc_save[79];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[28] -= 1 * qdot;
    wdot[27] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 81: H + CH2CO <=> CH3 + CO */
    phi_f = sc[1]*sc[28];
    k_f = k_f_save[80];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[14];
    Kc = Kc_save[80];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[28] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 82: H + HCCOH <=> H + CH2CO */
    phi_f = sc[1]*sc[29];
    k_f = k_f_save[81];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[81];
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
    k_f = k_f_save[82];
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
    phi_r = sc[17];
    Kc = Kc_save[82];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 84: OH + H2 <=> H + H2O */
    phi_f = sc[4]*sc[0];
    k_f = k_f_save[83];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    Kc = Kc_save[83];
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
    k_f = k_f_save[84];
    redP = 1e-12 * alpha / k_f * 2.3e+18*exp(-0.9*tc[0]+855.46833186802348337*invT);
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
    Kc = Kc_save[84];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 2 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 86: 2 OH <=> O + H2O */
    phi_f = sc[4]*sc[4];
    k_f = k_f_save[85];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = Kc_save[85];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 2 * qdot;
    wdot[2] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 87: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = k_f_save[86];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = Kc_save[86];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    phi_f = sc[4]*sc[7];
    k_f = k_f_save[87];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = Kc_save[87];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 89: OH + H2O2 <=> HO2 + H2O */
    phi_f = sc[4]*sc[7];
    k_f = k_f_save[88];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = Kc_save[88];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 90: OH + C <=> H + CO */
    phi_f = sc[4]*sc[8];
    k_f = k_f_save[89];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[89];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 91: OH + CH <=> H + HCO */
    phi_f = sc[4]*sc[9];
    k_f = k_f_save[90];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_save[90];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[9] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 92: OH + CH2 <=> H + CH2O */
    phi_f = sc[4]*sc[10];
    k_f = k_f_save[91];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[91];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 93: OH + CH2 <=> CH + H2O */
    phi_f = sc[4]*sc[10];
    k_f = k_f_save[92];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[5];
    Kc = Kc_save[92];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 94: OH + CH2(S) <=> H + CH2O */
    phi_f = sc[4]*sc[11];
    k_f = k_f_save[93];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[93];
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
    k_f = k_f_save[94];
    redP = 1e-12 * alpha / k_f * 2.7e+38*exp(-6.3*tc[0]-1559.9716639946311716*invT);
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
    Kc = Kc_save[94];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[20] += 1 * qdot;

    /*reaction 96: OH + CH3 <=> CH2 + H2O */
    phi_f = sc[4]*sc[12];
    k_f = k_f_save[95];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[5];
    Kc = Kc_save[95];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 97: OH + CH3 <=> CH2(S) + H2O */
    phi_f = sc[4]*sc[12];
    k_f = k_f_save[96];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[96];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 98: OH + CH4 <=> CH3 + H2O */
    phi_f = sc[4]*sc[13];
    k_f = k_f_save[97];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[5];
    Kc = Kc_save[97];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 99: OH + CO <=> H + CO2 */
    phi_f = sc[4]*sc[14];
    k_f = k_f_save[98];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[15];
    Kc = Kc_save[98];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 100: OH + HCO <=> H2O + CO */
    phi_f = sc[4]*sc[16];
    k_f = k_f_save[99];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[14];
    Kc = Kc_save[99];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[16] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 101: OH + CH2O <=> HCO + H2O */
    phi_f = sc[4]*sc[17];
    k_f = k_f_save[100];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[5];
    Kc = Kc_save[100];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 102: OH + CH2OH <=> H2O + CH2O */
    phi_f = sc[4]*sc[18];
    k_f = k_f_save[101];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[17];
    Kc = Kc_save[101];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 103: OH + CH3O <=> H2O + CH2O */
    phi_f = sc[4]*sc[19];
    k_f = k_f_save[102];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[17];
    Kc = Kc_save[102];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 104: OH + CH3OH <=> CH2OH + H2O */
    phi_f = sc[4]*sc[20];
    k_f = k_f_save[103];
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[5];
    Kc = Kc_save[103];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 105: OH + CH3OH <=> CH3O + H2O */
    phi_f = sc[4]*sc[20];
    k_f = k_f_save[104];
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[5];
    Kc = Kc_save[104];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[19] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 106: OH + C2H <=> H + HCCO */
    phi_f = sc[4]*sc[21];
    k_f = k_f_save[105];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[27];
    Kc = Kc_save[105];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[21] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[27] += 1 * qdot;

    /*reaction 107: OH + C2H2 <=> H + CH2CO */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[106];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[106];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[28] += 1 * qdot;

    /*reaction 108: OH + C2H2 <=> H + HCCOH */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[107];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[29];
    Kc = Kc_save[107];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[29] += 1 * qdot;

    /*reaction 109: OH + C2H2 <=> C2H + H2O */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[108];
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[5];
    Kc = Kc_save[108];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[21] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 110: OH + C2H2 <=> CH3 + CO */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[109];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[14];
    Kc = Kc_save[109];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 111: OH + C2H3 <=> H2O + C2H2 */
    phi_f = sc[4]*sc[23];
    k_f = k_f_save[110];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[22];
    Kc = Kc_save[110];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[23] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 112: OH + C2H4 <=> C2H3 + H2O */
    phi_f = sc[4]*sc[24];
    k_f = k_f_save[111];
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[5];
    Kc = Kc_save[111];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[24] -= 1 * qdot;
    wdot[23] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 113: OH + C2H6 <=> C2H5 + H2O */
    phi_f = sc[4]*sc[26];
    k_f = k_f_save[112];
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[5];
    Kc = Kc_save[112];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[26] -= 1 * qdot;
    wdot[25] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 114: OH + CH2CO <=> HCCO + H2O */
    phi_f = sc[4]*sc[28];
    k_f = k_f_save[113];
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[5];
    Kc = Kc_save[113];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[28] -= 1 * qdot;
    wdot[27] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    phi_f = sc[6]*sc[6];
    k_f = k_f_save[114];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = Kc_save[114];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 2 * qdot;
    wdot[3] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 116: 2 HO2 <=> O2 + H2O2 */
    phi_f = sc[6]*sc[6];
    k_f = k_f_save[115];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = Kc_save[115];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 2 * qdot;
    wdot[3] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 117: HO2 + CH2 <=> OH + CH2O */
    phi_f = sc[6]*sc[10];
    k_f = k_f_save[116];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[116];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 118: HO2 + CH3 <=> O2 + CH4 */
    phi_f = sc[6]*sc[12];
    k_f = k_f_save[117];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[13];
    Kc = Kc_save[117];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 119: HO2 + CH3 <=> OH + CH3O */
    phi_f = sc[6]*sc[12];
    k_f = k_f_save[118];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[19];
    Kc = Kc_save[118];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[19] += 1 * qdot;

    /*reaction 120: HO2 + CO <=> OH + CO2 */
    phi_f = sc[6]*sc[14];
    k_f = k_f_save[119];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[15];
    Kc = Kc_save[119];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 121: HO2 + CH2O <=> HCO + H2O2 */
    phi_f = sc[6]*sc[17];
    k_f = k_f_save[120];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[7];
    Kc = Kc_save[120];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 122: C + O2 <=> O + CO */
    phi_f = sc[8]*sc[3];
    k_f = k_f_save[121];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[14];
    Kc = Kc_save[121];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 123: C + CH2 <=> H + C2H */
    phi_f = sc[8]*sc[10];
    k_f = k_f_save[122];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[21];
    Kc = Kc_save[122];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[21] += 1 * qdot;

    /*reaction 124: C + CH3 <=> H + C2H2 */
    phi_f = sc[8]*sc[12];
    k_f = k_f_save[123];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = Kc_save[123];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 125: CH + O2 <=> O + HCO */
    phi_f = sc[9]*sc[3];
    k_f = k_f_save[124];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[16];
    Kc = Kc_save[124];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 126: CH + H2 <=> H + CH2 */
    phi_f = sc[9]*sc[0];
    k_f = k_f_save[125];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[10];
    Kc = Kc_save[125];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 127: CH + H2O <=> H + CH2O */
    phi_f = sc[9]*sc[5];
    k_f = k_f_save[126];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[126];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 128: CH + CH2 <=> H + C2H2 */
    phi_f = sc[9]*sc[10];
    k_f = k_f_save[127];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = Kc_save[127];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 129: CH + CH3 <=> H + C2H3 */
    phi_f = sc[9]*sc[12];
    k_f = k_f_save[128];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[23];
    Kc = Kc_save[128];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[23] += 1 * qdot;

    /*reaction 130: CH + CH4 <=> H + C2H4 */
    phi_f = sc[9]*sc[13];
    k_f = k_f_save[129];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = Kc_save[129];
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
    k_f = k_f_save[130];
    redP = 1e-12 * alpha / k_f * 2.69e+28*exp(-3.74*tc[0]-974.227464997937318*invT);
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
    Kc = Kc_save[130];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[27] += 1 * qdot;

    /*reaction 132: CH + CO2 <=> HCO + CO */
    phi_f = sc[9]*sc[15];
    k_f = k_f_save[131];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[14];
    Kc = Kc_save[131];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 133: CH + CH2O <=> H + CH2CO */
    phi_f = sc[9]*sc[17];
    k_f = k_f_save[132];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[132];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[28] += 1 * qdot;

    /*reaction 134: CH + HCCO <=> CO + C2H2 */
    phi_f = sc[9]*sc[27];
    k_f = k_f_save[133];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[22];
    Kc = Kc_save[133];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[27] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 135: CH2 + O2 <=> OH + HCO */
    phi_f = sc[10]*sc[3];
    k_f = k_f_save[134];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[16];
    Kc = Kc_save[134];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 136: CH2 + H2 <=> H + CH3 */
    phi_f = sc[10]*sc[0];
    k_f = k_f_save[135];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = Kc_save[135];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 137: 2 CH2 <=> H2 + C2H2 */
    phi_f = sc[10]*sc[10];
    k_f = k_f_save[136];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = Kc_save[136];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 2 * qdot;
    wdot[0] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 138: CH2 + CH3 <=> H + C2H4 */
    phi_f = sc[10]*sc[12];
    k_f = k_f_save[137];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = Kc_save[137];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[24] += 1 * qdot;

    /*reaction 139: CH2 + CH4 <=> 2 CH3 */
    phi_f = sc[10]*sc[13];
    k_f = k_f_save[138];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[12];
    Kc = Kc_save[138];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[12] += 2 * qdot;

    /*reaction 140: CH2 + CO (+M) <=> CH2CO (+M) */
    phi_f = sc[10]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[139];
    redP = 1e-12 * alpha / k_f * 2.69e+33*exp(-5.11*tc[0]-3570.3222438844863973*invT);
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
    Kc = Kc_save[139];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[28] += 1 * qdot;

    /*reaction 141: CH2 + HCCO <=> C2H3 + CO */
    phi_f = sc[10]*sc[27];
    k_f = k_f_save[140];
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[14];
    Kc = Kc_save[140];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[27] -= 1 * qdot;
    wdot[23] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 142: CH2(S) + N2 <=> CH2 + N2 */
    phi_f = sc[11]*sc[30];
    k_f = k_f_save[141];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[30];
    Kc = Kc_save[141];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[30] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[30] += 1 * qdot;

    /*reaction 143: CH2(S) + AR <=> CH2 + AR */
    phi_f = sc[11]*sc[31];
    k_f = k_f_save[142];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[31];
    Kc = Kc_save[142];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[31] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[31] += 1 * qdot;

    /*reaction 144: CH2(S) + O2 <=> H + OH + CO */
    phi_f = sc[11]*sc[3];
    k_f = k_f_save[143];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[4]*sc[14];
    Kc = Kc_save[143];
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
    k_f = k_f_save[144];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[5];
    Kc = Kc_save[144];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 146: CH2(S) + H2 <=> CH3 + H */
    phi_f = sc[11]*sc[0];
    k_f = k_f_save[145];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[1];
    Kc = Kc_save[145];
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
    k_f = k_f_save[146];
    redP = 1e-12 * alpha / k_f * 2.7e+38*exp(-6.3*tc[0]-1559.9716639946311716*invT);
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
    Kc = Kc_save[146];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[20] += 1 * qdot;

    /*reaction 148: CH2(S) + H2O <=> CH2 + H2O */
    phi_f = sc[11]*sc[5];
    k_f = k_f_save[147];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[5];
    Kc = Kc_save[147];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 149: CH2(S) + CH3 <=> H + C2H4 */
    phi_f = sc[11]*sc[12];
    k_f = k_f_save[148];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = Kc_save[148];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[24] += 1 * qdot;

    /*reaction 150: CH2(S) + CH4 <=> 2 CH3 */
    phi_f = sc[11]*sc[13];
    k_f = k_f_save[149];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[12];
    Kc = Kc_save[149];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[12] += 2 * qdot;

    /*reaction 151: CH2(S) + CO <=> CH2 + CO */
    phi_f = sc[11]*sc[14];
    k_f = k_f_save[150];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[14];
    Kc = Kc_save[150];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 152: CH2(S) + CO2 <=> CH2 + CO2 */
    phi_f = sc[11]*sc[15];
    k_f = k_f_save[151];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[15];
    Kc = Kc_save[151];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 153: CH2(S) + CO2 <=> CO + CH2O */
    phi_f = sc[11]*sc[15];
    k_f = k_f_save[152];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[17];
    Kc = Kc_save[152];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 154: CH2(S) + C2H6 <=> CH3 + C2H5 */
    phi_f = sc[11]*sc[26];
    k_f = k_f_save[153];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[25];
    Kc = Kc_save[153];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[26] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[25] += 1 * qdot;

    /*reaction 155: CH3 + O2 <=> O + CH3O */
    phi_f = sc[12]*sc[3];
    k_f = k_f_save[154];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[19];
    Kc = Kc_save[154];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[19] += 1 * qdot;

    /*reaction 156: CH3 + O2 <=> OH + CH2O */
    phi_f = sc[12]*sc[3];
    k_f = k_f_save[155];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[155];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 157: CH3 + H2O2 <=> HO2 + CH4 */
    phi_f = sc[12]*sc[7];
    k_f = k_f_save[156];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = Kc_save[156];
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
    k_f = k_f_save[157];
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
    phi_r = sc[26];
    Kc = Kc_save[157];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 2 * qdot;
    wdot[26] += 1 * qdot;

    /*reaction 159: 2 CH3 <=> H + C2H5 */
    phi_f = sc[12]*sc[12];
    k_f = k_f_save[158];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[25];
    Kc = Kc_save[158];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 2 * qdot;
    wdot[1] += 1 * qdot;
    wdot[25] += 1 * qdot;

    /*reaction 160: CH3 + HCO <=> CH4 + CO */
    phi_f = sc[12]*sc[16];
    k_f = k_f_save[159];
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[14];
    Kc = Kc_save[159];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[16] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 161: CH3 + CH2O <=> HCO + CH4 */
    phi_f = sc[12]*sc[17];
    k_f = k_f_save[160];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[13];
    Kc = Kc_save[160];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 162: CH3 + CH3OH <=> CH2OH + CH4 */
    phi_f = sc[12]*sc[20];
    k_f = k_f_save[161];
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[13];
    Kc = Kc_save[161];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 163: CH3 + CH3OH <=> CH3O + CH4 */
    phi_f = sc[12]*sc[20];
    k_f = k_f_save[162];
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[13];
    Kc = Kc_save[162];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[19] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 164: CH3 + C2H4 <=> C2H3 + CH4 */
    phi_f = sc[12]*sc[24];
    k_f = k_f_save[163];
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[13];
    Kc = Kc_save[163];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[24] -= 1 * qdot;
    wdot[23] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 165: CH3 + C2H6 <=> C2H5 + CH4 */
    phi_f = sc[12]*sc[26];
    k_f = k_f_save[164];
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[13];
    Kc = Kc_save[164];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[26] -= 1 * qdot;
    wdot[25] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 166: HCO + H2O <=> H + CO + H2O */
    phi_f = sc[16]*sc[5];
    k_f = k_f_save[165];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14]*sc[5];
    Kc = Kc_save[165];
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
    k_f = alpha * k_f_save[166];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[166];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 168: HCO + O2 <=> HO2 + CO */
    phi_f = sc[16]*sc[3];
    k_f = k_f_save[167];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[14];
    Kc = Kc_save[167];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 169: CH2OH + O2 <=> HO2 + CH2O */
    phi_f = sc[18]*sc[3];
    k_f = k_f_save[168];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = Kc_save[168];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 170: CH3O + O2 <=> HO2 + CH2O */
    phi_f = sc[19]*sc[3];
    k_f = k_f_save[169];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = Kc_save[169];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 171: C2H + O2 <=> HCO + CO */
    phi_f = sc[21]*sc[3];
    k_f = k_f_save[170];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[14];
    Kc = Kc_save[170];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 172: C2H + H2 <=> H + C2H2 */
    phi_f = sc[21]*sc[0];
    k_f = k_f_save[171];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = Kc_save[171];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 173: C2H3 + O2 <=> HCO + CH2O */
    phi_f = sc[23]*sc[3];
    k_f = k_f_save[172];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[17];
    Kc = Kc_save[172];
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
    k_f = k_f_save[173];
    redP = 1e-6 * alpha / k_f * 7e+50*exp(-9.31*tc[0]-50251.216247259304509*invT);
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
    Kc = Kc_save[173];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[24] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 175: C2H5 + O2 <=> HO2 + C2H4 */
    phi_f = sc[25]*sc[3];
    k_f = k_f_save[174];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[24];
    Kc = Kc_save[174];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[25] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[24] += 1 * qdot;

    /*reaction 176: HCCO + O2 <=> OH + 2 CO */
    phi_f = sc[27]*sc[3];
    k_f = k_f_save[175];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[14]*sc[14];
    Kc = Kc_save[175];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[27] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[14] += 2 * qdot;

    /*reaction 177: 2 HCCO <=> 2 CO + C2H2 */
    phi_f = sc[27]*sc[27];
    k_f = k_f_save[176];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[14]*sc[22];
    Kc = Kc_save[176];
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

    double invT = 1.0 / tc[1];

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 32; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    if (T != T_save)
    {
        T_save = T;

        k_f_save[0] = 1e-12 * 1.2e+17*exp(-1*tc[0]);
        k_f_save[1] = 1e-12 * 5e+17*exp(-1*tc[0]);
        k_f_save[2] = 1e-06 * 50000*exp(2.67*tc[0]-3165.2328279116868543*invT);
        k_f_save[3] = 1e-06 * 2e+13;
        k_f_save[4] = 1e-06 * 9.63e+06*exp(2*tc[0]-2012.866663218878557*invT);
        k_f_save[5] = 1e-06 * 5.7e+13;
        k_f_save[6] = 1e-06 * 8e+13;
        k_f_save[7] = 1e-06 * 1.5e+13;
        k_f_save[8] = 1e-06 * 1.5e+13;
        k_f_save[9] = 1e-06 * 8.43e+13;
        k_f_save[10] = 1e-06 * 1.02e+09*exp(1.5*tc[0]-4327.6633259205891591*invT);
        k_f_save[11] = 1e-12 * 6.02e+14*exp(-1509.6499974141590883*invT);
        k_f_save[12] = 1e-06 * 3e+13;
        k_f_save[13] = 1e-06 * 3e+13;
        k_f_save[14] = 1e-06 * 3.9e+13*exp(-1781.3869969487077469*invT);
        k_f_save[15] = 1e-06 * 1e+13;
        k_f_save[16] = 1e-06 * 1e+13;
        k_f_save[17] = 1e-06 * 388000*exp(2.5*tc[0]-1559.9716639946311716*invT);
        k_f_save[18] = 1e-06 * 130000*exp(2.5*tc[0]-2516.0833290235987079*invT);
        k_f_save[19] = 1e-06 * 5e+13;
        k_f_save[20] = 1e-06 * 1.02e+07*exp(2*tc[0]-956.11166502896742259*invT);
        k_f_save[21] = 1e-06 * 4.6e+19*exp(-1.41*tc[0]-14568.122475046635373*invT);
        k_f_save[22] = 1e-06 * 1.02e+07*exp(2*tc[0]-956.11166502896742259*invT);
        k_f_save[23] = 1e-06 * 3e+13;
        k_f_save[24] = 1e-06 * 1.92e+07*exp(1.83*tc[0]-110.7076664770383303*invT);
        k_f_save[25] = 1e-06 * 1.32e+14;
        k_f_save[26] = 1e-06 * 8.98e+07*exp(1.92*tc[0]-2863.3028284288548093*invT);
        k_f_save[27] = 1e-06 * 1e+14;
        k_f_save[28] = 1e-06 * 1e+13*exp(-4025.733326437757114*invT);
        k_f_save[29] = 1e-06 * 1.75e+12*exp(-679.34249883637164658*invT);
        k_f_save[30] = 1e-06 * 2.5e+12*exp(-24053.756625465601246*invT);
        k_f_save[31] = 1e-06 * 1e+14*exp(-20128.666632188789663*invT);
        k_f_save[32] = 1e-12 * 2.8e+18*exp(-0.86*tc[0]);
        k_f_save[33] = 1e-12 * 3e+20*exp(-1.72*tc[0]);
        k_f_save[34] = 1e-12 * 9.38e+18*exp(-0.76*tc[0]);
        k_f_save[35] = 1e-12 * 3.75e+20*exp(-1.72*tc[0]);
        k_f_save[36] = 1e-12 * 7e+17*exp(-0.8*tc[0]);
        k_f_save[37] = 1e-06 * 8.3e+13*exp(-7252.8618042434254676*invT);
        k_f_save[38] = 1e-12 * 1e+18*exp(-1*tc[0]);
        k_f_save[39] = 1e-12 * 9e+16*exp(-0.6*tc[0]);
        k_f_save[40] = 1e-12 * 6e+19*exp(-1.25*tc[0]);
        k_f_save[41] = 1e-12 * 5.5e+20*exp(-2*tc[0]);
        k_f_save[42] = 1e-12 * 2.2e+22*exp(-2*tc[0]);
        k_f_save[43] = 1e-06 * 3.97e+12*exp(-337.65838275496690812*invT);
        k_f_save[44] = 1e-06 * 2.8e+13*exp(-537.43539907944057177*invT);
        k_f_save[45] = 1e-06 * 1.34e+14*exp(-319.54258278599701271*invT);
        k_f_save[46] = 1e-06 * 1.21e+07*exp(2*tc[0]-2616.7266621845424197*invT);
        k_f_save[47] = 1e-06 * 1e+13*exp(-1811.579996896990906*invT);
        k_f_save[48] = 1e-06 * 1.1e+14;
        k_f_save[49] = 1e-06 * 2.5e+16*exp(-0.8*tc[0]);
        k_f_save[50] = 1e-06 * 3e+13;
        k_f_save[51] = 1e-06 * 1.27e+16*exp(-0.63*tc[0]-192.73198300320765952*invT);
        k_f_save[52] = 1e-06 * 6.6e+08*exp(1.62*tc[0]-5454.8686573231616421*invT);
        k_f_save[53] = 1e-06 * 1.09e+12*exp(0.48*tc[0]+130.83633310922709825*invT);
        k_f_save[54] = 1e-06 * 7.34e+13;
        k_f_save[55] = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1811.579996896990906*invT);
        k_f_save[56] = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1308.3633310922712099*invT);
        k_f_save[57] = 1e-06 * 2.3e+10*exp(1.05*tc[0]-1648.0345805104568626*invT);
        k_f_save[58] = 1e-06 * 1.8e+13;
        k_f_save[59] = 1e-06 * 2e+13;
        k_f_save[60] = 1e-06 * 1.2e+13;
        k_f_save[61] = 1e-06 * 6e+12;
        k_f_save[62] = 1e-06 * 5e+13;
        k_f_save[63] = 1e-06 * 3.4e+06*exp(1.6*tc[0]);
        k_f_save[64] = 1e-06 * 2e+13;
        k_f_save[65] = 1e-06 * 3.2e+13;
        k_f_save[66] = 1e-06 * 1.6e+13;
        k_f_save[67] = 1e-06 * 1.7e+07*exp(2.1*tc[0]-2450.6651624689848177*invT);
        k_f_save[68] = 1e-06 * 4.2e+06*exp(2.1*tc[0]-2450.6651624689848177*invT);
        k_f_save[69] = 1e-06 * 1e+17*exp(-1*tc[0]);
        k_f_save[70] = 1e-06 * 5.6e+12*exp(-1207.7199979313272706*invT);
        k_f_save[71] = 1e-06 * 6.08e+12*exp(0.27*tc[0]-140.90066642532150354*invT);
        k_f_save[72] = 1e-06 * 3e+13;
        k_f_save[73] = 1e-06 * 1.08e+12*exp(0.454*tc[0]-915.85433176458980142*invT);
        k_f_save[74] = 1e-06 * 1.325e+06*exp(2.53*tc[0]-6159.3719894497689893*invT);
        k_f_save[75] = 1e-06 * 5.21e+17*exp(-0.99*tc[0]-795.08233197145705162*invT);
        k_f_save[76] = 1e-06 * 2e+12;
        k_f_save[77] = 1e-06 * 1.15e+08*exp(1.9*tc[0]-3789.2214935095394139*invT);
        k_f_save[78] = 1e-06 * 1e+14;
        k_f_save[79] = 1e-06 * 5e+13*exp(-4025.733326437757114*invT);
        k_f_save[80] = 1e-06 * 1.13e+13*exp(-1725.0267303785790318*invT);
        k_f_save[81] = 1e-06 * 1e+13;
        k_f_save[82] = 1e-06 * 4.3e+07*exp(1.5*tc[0]-40056.046598055690993*invT);
        k_f_save[83] = 1e-06 * 2.16e+08*exp(1.51*tc[0]-1726.0331637101885462*invT);
        k_f_save[84] = 1e-06 * 7.4e+13*exp(-0.37*tc[0]);
        k_f_save[85] = 1e-06 * 35700*exp(2.4*tc[0]+1061.7871648479585929*invT);
        k_f_save[86] = 1e-06 * 2.9e+13*exp(+251.60833290235981963*invT);
        k_f_save[87] = 1e-06 * 1.75e+12*exp(-161.02933305751031412*invT);
        k_f_save[88] = 1e-06 * 5.8e+14*exp(-4810.751325093120613*invT);
        k_f_save[89] = 1e-06 * 5e+13;
        k_f_save[90] = 1e-06 * 3e+13;
        k_f_save[91] = 1e-06 * 2e+13;
        k_f_save[92] = 1e-06 * 1.13e+07*exp(2*tc[0]-1509.6499974141590883*invT);
        k_f_save[93] = 1e-06 * 3e+13;
        k_f_save[94] = 1e-06 * 6.3e+13;
        k_f_save[95] = 1e-06 * 5.6e+07*exp(1.6*tc[0]-2727.4343286615808211*invT);
        k_f_save[96] = 1e-06 * 2.501e+13;
        k_f_save[97] = 1e-06 * 1e+08*exp(1.6*tc[0]-1570.0359973107254064*invT);
        k_f_save[98] = 1e-06 * 4.76e+07*exp(1.228*tc[0]-35.225166606330375885*invT);
        k_f_save[99] = 1e-06 * 5e+13;
        k_f_save[100] = 1e-06 * 3.43e+09*exp(1.18*tc[0]+224.93784961470970529*invT);
        k_f_save[101] = 1e-06 * 5e+12;
        k_f_save[102] = 1e-06 * 5e+12;
        k_f_save[103] = 1e-06 * 1.44e+06*exp(2*tc[0]+422.70199927596451062*invT);
        k_f_save[104] = 1e-06 * 6.3e+06*exp(2*tc[0]-754.82499870707954415*invT);
        k_f_save[105] = 1e-06 * 2e+13;
        k_f_save[106] = 1e-06 * 0.000218*exp(4.5*tc[0]+503.21666580471963925*invT);
        k_f_save[107] = 1e-06 * 504000*exp(2.3*tc[0]-6793.4249883637157836*invT);
        k_f_save[108] = 1e-06 * 3.37e+07*exp(2*tc[0]-7045.0333212660752906*invT);
        k_f_save[109] = 1e-06 * 0.000483*exp(4*tc[0]+1006.4333316094392785*invT);
        k_f_save[110] = 1e-06 * 5e+12;
        k_f_save[111] = 1e-06 * 3.6e+06*exp(2*tc[0]-1258.0416645117993539*invT);
        k_f_save[112] = 1e-06 * 3.54e+06*exp(2.12*tc[0]-437.79849925010609013*invT);
        k_f_save[113] = 1e-06 * 7.5e+12*exp(-1006.4333316094392785*invT);
        k_f_save[114] = 1e-06 * 1.3e+11*exp(+820.24316526169309327*invT);
        k_f_save[115] = 1e-06 * 4.2e+14*exp(-6038.5999896566363532*invT);
        k_f_save[116] = 1e-06 * 2e+13;
        k_f_save[117] = 1e-06 * 1e+12;
        k_f_save[118] = 1e-06 * 2e+13;
        k_f_save[119] = 1e-06 * 1.5e+14*exp(-11875.913312991384373*invT);
        k_f_save[120] = 1e-06 * 1e+12*exp(-4025.733326437757114*invT);
        k_f_save[121] = 1e-06 * 5.8e+13*exp(-289.85279950351855405*invT);
        k_f_save[122] = 1e-06 * 5e+13;
        k_f_save[123] = 1e-06 * 5e+13;
        k_f_save[124] = 1e-06 * 3.3e+13;
        k_f_save[125] = 1e-06 * 1.107e+08*exp(1.79*tc[0]-840.37183189388190385*invT);
        k_f_save[126] = 1e-06 * 5.71e+12*exp(+379.92858268256338761*invT);
        k_f_save[127] = 1e-06 * 4e+13;
        k_f_save[128] = 1e-06 * 3e+13;
        k_f_save[129] = 1e-06 * 6e+13;
        k_f_save[130] = 1e-06 * 5e+13;
        k_f_save[131] = 1e-06 * 3.4e+12*exp(-347.2194994052565562*invT);
        k_f_save[132] = 1e-06 * 9.46e+13*exp(+259.15658288943063781*invT);
        k_f_save[133] = 1e-06 * 5e+13;
        k_f_save[134] = 1e-06 * 1.32e+13*exp(-754.82499870707954415*invT);
        k_f_save[135] = 1e-06 * 500000*exp(2*tc[0]-3638.256493768123164*invT);
        k_f_save[136] = 1e-06 * 3.2e+13;
        k_f_save[137] = 1e-06 * 4e+13;
        k_f_save[138] = 1e-06 * 2.46e+06*exp(2*tc[0]-4161.6018262050320118*invT);
        k_f_save[139] = 1e-06 * 8.1e+11*exp(0.5*tc[0]-2269.5071627792858635*invT);
        k_f_save[140] = 1e-06 * 3e+13;
        k_f_save[141] = 1e-06 * 1.5e+13*exp(-301.92999948283181766*invT);
        k_f_save[142] = 1e-06 * 9e+12*exp(-301.92999948283181766*invT);
        k_f_save[143] = 1e-06 * 2.8e+13;
        k_f_save[144] = 1e-06 * 1.2e+13;
        k_f_save[145] = 1e-06 * 7e+13;
        k_f_save[146] = 1e-06 * 2e+13;
        k_f_save[147] = 1e-06 * 3e+13;
        k_f_save[148] = 1e-06 * 1.2e+13*exp(+286.83349950869023814*invT);
        k_f_save[149] = 1e-06 * 1.6e+13*exp(+286.83349950869023814*invT);
        k_f_save[150] = 1e-06 * 9e+12;
        k_f_save[151] = 1e-06 * 7e+12;
        k_f_save[152] = 1e-06 * 1.4e+13;
        k_f_save[153] = 1e-06 * 4e+13*exp(+276.76916619259583285*invT);
        k_f_save[154] = 1e-06 * 2.675e+13*exp(-14492.639975175927248*invT);
        k_f_save[155] = 1e-06 * 3.6e+10*exp(-4498.7569922941938785*invT);
        k_f_save[156] = 1e-06 * 24500*exp(2.47*tc[0]-2606.6623288684481849*invT);
        k_f_save[157] = 1e-06 * 2.12e+16*exp(-0.97*tc[0]-311.99433279892622295*invT);
        k_f_save[158] = 1e-06 * 4.99e+12*exp(0.1*tc[0]-5334.096657530029006*invT);
        k_f_save[159] = 1e-06 * 2.648e+13;
        k_f_save[160] = 1e-06 * 3320*exp(2.81*tc[0]-2948.8496616156576238*invT);
        k_f_save[161] = 1e-06 * 3e+07*exp(1.5*tc[0]-5001.973658098913802*invT);
        k_f_save[162] = 1e-06 * 1e+07*exp(1.5*tc[0]-5001.973658098913802*invT);
        k_f_save[163] = 1e-06 * 227000*exp(2*tc[0]-4629.5933254034207494*invT);
        k_f_save[164] = 1e-06 * 6.14e+06*exp(1.74*tc[0]-5258.6141576593208811*invT);
        k_f_save[165] = 1e-06 * 2.244e+18*exp(-1*tc[0]-8554.6833186802341515*invT);
        k_f_save[166] = 1e-06 * 1.87e+17*exp(-1*tc[0]-8554.6833186802341515*invT);
        k_f_save[167] = 1e-06 * 7.6e+12*exp(-201.28666632188787844*invT);
        k_f_save[168] = 1e-06 * 1.8e+13*exp(-452.89499922424772649*invT);
        k_f_save[169] = 1e-06 * 4.28e-13*exp(7.6*tc[0]+1776.3548302906606295*invT);
        k_f_save[170] = 1e-06 * 5e+13*exp(-754.82499870707954415*invT);
        k_f_save[171] = 1e-06 * 407000*exp(2.4*tc[0]-100.64333316094393922*invT);
        k_f_save[172] = 1e-06 * 3.98e+12*exp(+120.77199979313272138*invT);
        k_f_save[173] = 1 * 8e+12*exp(0.44*tc[0]-44670.543423484967207*invT);
        k_f_save[174] = 1e-06 * 8.4e+11*exp(-1949.9645799932889076*invT);
        k_f_save[175] = 1e-06 * 1.6e+12*exp(-429.7470325972306*invT);
        k_f_save[176] = 1e-06 * 1e+13;

        Kc_save[0] = 1.0 / (refC) * exp((2 * g_RT[2]) - (g_RT[3]));
        Kc_save[1] = 1.0 / (refC) * exp((g_RT[2] + g_RT[1]) - (g_RT[4]));
        Kc_save[2] = exp((g_RT[2] + g_RT[0]) - (g_RT[1] + g_RT[4]));
        Kc_save[3] = exp((g_RT[2] + g_RT[6]) - (g_RT[4] + g_RT[3]));
        Kc_save[4] = exp((g_RT[2] + g_RT[7]) - (g_RT[4] + g_RT[6]));
        Kc_save[5] = exp((g_RT[2] + g_RT[9]) - (g_RT[1] + g_RT[14]));
        Kc_save[6] = exp((g_RT[2] + g_RT[10]) - (g_RT[1] + g_RT[16]));
        Kc_save[7] = exp((g_RT[2] + g_RT[11]) - (g_RT[0] + g_RT[14]));
        Kc_save[8] = exp((g_RT[2] + g_RT[11]) - (g_RT[1] + g_RT[16]));
        Kc_save[9] = exp((g_RT[2] + g_RT[12]) - (g_RT[1] + g_RT[17]));
        Kc_save[10] = exp((g_RT[2] + g_RT[13]) - (g_RT[4] + g_RT[12]));
        Kc_save[11] = 1.0 / (refC) * exp((g_RT[2] + g_RT[14]) - (g_RT[15]));
        Kc_save[12] = exp((g_RT[2] + g_RT[16]) - (g_RT[4] + g_RT[14]));
        Kc_save[13] = exp((g_RT[2] + g_RT[16]) - (g_RT[1] + g_RT[15]));
        Kc_save[14] = exp((g_RT[2] + g_RT[17]) - (g_RT[4] + g_RT[16]));
        Kc_save[15] = exp((g_RT[2] + g_RT[18]) - (g_RT[4] + g_RT[17]));
        Kc_save[16] = exp((g_RT[2] + g_RT[19]) - (g_RT[4] + g_RT[17]));
        Kc_save[17] = exp((g_RT[2] + g_RT[20]) - (g_RT[4] + g_RT[18]));
        Kc_save[18] = exp((g_RT[2] + g_RT[20]) - (g_RT[4] + g_RT[19]));
        Kc_save[19] = exp((g_RT[2] + g_RT[21]) - (g_RT[9] + g_RT[14]));
        Kc_save[20] = exp((g_RT[2] + g_RT[22]) - (g_RT[1] + g_RT[27]));
        Kc_save[21] = exp((g_RT[2] + g_RT[22]) - (g_RT[4] + g_RT[21]));
        Kc_save[22] = exp((g_RT[2] + g_RT[22]) - (g_RT[14] + g_RT[10]));
        Kc_save[23] = exp((g_RT[2] + g_RT[23]) - (g_RT[1] + g_RT[28]));
        Kc_save[24] = exp((g_RT[2] + g_RT[24]) - (g_RT[12] + g_RT[16]));
        Kc_save[25] = exp((g_RT[2] + g_RT[25]) - (g_RT[12] + g_RT[17]));
        Kc_save[26] = exp((g_RT[2] + g_RT[26]) - (g_RT[4] + g_RT[25]));
        Kc_save[27] = refC * exp((g_RT[2] + g_RT[27]) - (g_RT[1] + 2 * g_RT[14]));
        Kc_save[28] = exp((g_RT[2] + g_RT[28]) - (g_RT[4] + g_RT[27]));
        Kc_save[29] = exp((g_RT[2] + g_RT[28]) - (g_RT[10] + g_RT[15]));
        Kc_save[30] = exp((g_RT[3] + g_RT[14]) - (g_RT[2] + g_RT[15]));
        Kc_save[31] = exp((g_RT[3] + g_RT[17]) - (g_RT[6] + g_RT[16]));
        Kc_save[32] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));
        Kc_save[33] = 1.0 / (refC) * exp((g_RT[1] + 2 * g_RT[3]) - (g_RT[6] + g_RT[3]));
        Kc_save[34] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));
        Kc_save[35] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[30]) - (g_RT[6] + g_RT[30]));
        Kc_save[36] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[31]) - (g_RT[6] + g_RT[31]));
        Kc_save[37] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));
        Kc_save[38] = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));
        Kc_save[39] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));
        Kc_save[40] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));
        Kc_save[41] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[15]) - (g_RT[0] + g_RT[15]));
        Kc_save[42] = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));
        Kc_save[43] = exp((g_RT[1] + g_RT[6]) - (g_RT[2] + g_RT[5]));
        Kc_save[44] = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));
        Kc_save[45] = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));
        Kc_save[46] = exp((g_RT[1] + g_RT[7]) - (g_RT[6] + g_RT[0]));
        Kc_save[47] = exp((g_RT[1] + g_RT[7]) - (g_RT[4] + g_RT[5]));
        Kc_save[48] = exp((g_RT[1] + g_RT[9]) - (g_RT[8] + g_RT[0]));
        Kc_save[49] = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[12]));
        Kc_save[50] = exp((g_RT[1] + g_RT[11]) - (g_RT[9] + g_RT[0]));
        Kc_save[51] = 1.0 / (refC) * exp((g_RT[1] + g_RT[12]) - (g_RT[13]));
        Kc_save[52] = exp((g_RT[1] + g_RT[13]) - (g_RT[12] + g_RT[0]));
        Kc_save[53] = 1.0 / (refC) * exp((g_RT[1] + g_RT[16]) - (g_RT[17]));
        Kc_save[54] = exp((g_RT[1] + g_RT[16]) - (g_RT[0] + g_RT[14]));
        Kc_save[55] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[18]));
        Kc_save[56] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[19]));
        Kc_save[57] = exp((g_RT[1] + g_RT[17]) - (g_RT[16] + g_RT[0]));
        Kc_save[58] = 1.0 / (refC) * exp((g_RT[1] + g_RT[18]) - (g_RT[20]));
        Kc_save[59] = exp((g_RT[1] + g_RT[18]) - (g_RT[0] + g_RT[17]));
        Kc_save[60] = exp((g_RT[1] + g_RT[18]) - (g_RT[4] + g_RT[12]));
        Kc_save[61] = exp((g_RT[1] + g_RT[18]) - (g_RT[11] + g_RT[5]));
        Kc_save[62] = 1.0 / (refC) * exp((g_RT[1] + g_RT[19]) - (g_RT[20]));
        Kc_save[63] = exp((g_RT[1] + g_RT[19]) - (g_RT[1] + g_RT[18]));
        Kc_save[64] = exp((g_RT[1] + g_RT[19]) - (g_RT[0] + g_RT[17]));
        Kc_save[65] = exp((g_RT[1] + g_RT[19]) - (g_RT[4] + g_RT[12]));
        Kc_save[66] = exp((g_RT[1] + g_RT[19]) - (g_RT[11] + g_RT[5]));
        Kc_save[67] = exp((g_RT[1] + g_RT[20]) - (g_RT[18] + g_RT[0]));
        Kc_save[68] = exp((g_RT[1] + g_RT[20]) - (g_RT[19] + g_RT[0]));
        Kc_save[69] = 1.0 / (refC) * exp((g_RT[1] + g_RT[21]) - (g_RT[22]));
        Kc_save[70] = 1.0 / (refC) * exp((g_RT[1] + g_RT[22]) - (g_RT[23]));
        Kc_save[71] = 1.0 / (refC) * exp((g_RT[1] + g_RT[23]) - (g_RT[24]));
        Kc_save[72] = exp((g_RT[1] + g_RT[23]) - (g_RT[0] + g_RT[22]));
        Kc_save[73] = 1.0 / (refC) * exp((g_RT[1] + g_RT[24]) - (g_RT[25]));
        Kc_save[74] = exp((g_RT[1] + g_RT[24]) - (g_RT[23] + g_RT[0]));
        Kc_save[75] = 1.0 / (refC) * exp((g_RT[1] + g_RT[25]) - (g_RT[26]));
        Kc_save[76] = exp((g_RT[1] + g_RT[25]) - (g_RT[0] + g_RT[24]));
        Kc_save[77] = exp((g_RT[1] + g_RT[26]) - (g_RT[25] + g_RT[0]));
        Kc_save[78] = exp((g_RT[1] + g_RT[27]) - (g_RT[11] + g_RT[14]));
        Kc_save[79] = exp((g_RT[1] + g_RT[28]) - (g_RT[27] + g_RT[0]));
        Kc_save[80] = exp((g_RT[1] + g_RT[28]) - (g_RT[12] + g_RT[14]));
        Kc_save[81] = exp((g_RT[1] + g_RT[29]) - (g_RT[1] + g_RT[28]));
        Kc_save[82] = 1.0 / (refC) * exp((g_RT[0] + g_RT[14]) - (g_RT[17]));
        Kc_save[83] = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));
        Kc_save[84] = 1.0 / (refC) * exp((2 * g_RT[4]) - (g_RT[7]));
        Kc_save[85] = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));
        Kc_save[86] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
        Kc_save[87] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
        Kc_save[88] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
        Kc_save[89] = exp((g_RT[4] + g_RT[8]) - (g_RT[1] + g_RT[14]));
        Kc_save[90] = exp((g_RT[4] + g_RT[9]) - (g_RT[1] + g_RT[16]));
        Kc_save[91] = exp((g_RT[4] + g_RT[10]) - (g_RT[1] + g_RT[17]));
        Kc_save[92] = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[5]));
        Kc_save[93] = exp((g_RT[4] + g_RT[11]) - (g_RT[1] + g_RT[17]));
        Kc_save[94] = 1.0 / (refC) * exp((g_RT[4] + g_RT[12]) - (g_RT[20]));
        Kc_save[95] = exp((g_RT[4] + g_RT[12]) - (g_RT[10] + g_RT[5]));
        Kc_save[96] = exp((g_RT[4] + g_RT[12]) - (g_RT[11] + g_RT[5]));
        Kc_save[97] = exp((g_RT[4] + g_RT[13]) - (g_RT[12] + g_RT[5]));
        Kc_save[98] = exp((g_RT[4] + g_RT[14]) - (g_RT[1] + g_RT[15]));
        Kc_save[99] = exp((g_RT[4] + g_RT[16]) - (g_RT[5] + g_RT[14]));
        Kc_save[100] = exp((g_RT[4] + g_RT[17]) - (g_RT[16] + g_RT[5]));
        Kc_save[101] = exp((g_RT[4] + g_RT[18]) - (g_RT[5] + g_RT[17]));
        Kc_save[102] = exp((g_RT[4] + g_RT[19]) - (g_RT[5] + g_RT[17]));
        Kc_save[103] = exp((g_RT[4] + g_RT[20]) - (g_RT[18] + g_RT[5]));
        Kc_save[104] = exp((g_RT[4] + g_RT[20]) - (g_RT[19] + g_RT[5]));
        Kc_save[105] = exp((g_RT[4] + g_RT[21]) - (g_RT[1] + g_RT[27]));
        Kc_save[106] = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[28]));
        Kc_save[107] = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[29]));
        Kc_save[108] = exp((g_RT[4] + g_RT[22]) - (g_RT[21] + g_RT[5]));
        Kc_save[109] = exp((g_RT[4] + g_RT[22]) - (g_RT[12] + g_RT[14]));
        Kc_save[110] = exp((g_RT[4] + g_RT[23]) - (g_RT[5] + g_RT[22]));
        Kc_save[111] = exp((g_RT[4] + g_RT[24]) - (g_RT[23] + g_RT[5]));
        Kc_save[112] = exp((g_RT[4] + g_RT[26]) - (g_RT[25] + g_RT[5]));
        Kc_save[113] = exp((g_RT[4] + g_RT[28]) - (g_RT[27] + g_RT[5]));
        Kc_save[114] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
        Kc_save[115] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
        Kc_save[116] = exp((g_RT[6] + g_RT[10]) - (g_RT[4] + g_RT[17]));
        Kc_save[117] = exp((g_RT[6] + g_RT[12]) - (g_RT[3] + g_RT[13]));
        Kc_save[118] = exp((g_RT[6] + g_RT[12]) - (g_RT[4] + g_RT[19]));
        Kc_save[119] = exp((g_RT[6] + g_RT[14]) - (g_RT[4] + g_RT[15]));
        Kc_save[120] = exp((g_RT[6] + g_RT[17]) - (g_RT[16] + g_RT[7]));
        Kc_save[121] = exp((g_RT[8] + g_RT[3]) - (g_RT[2] + g_RT[14]));
        Kc_save[122] = exp((g_RT[8] + g_RT[10]) - (g_RT[1] + g_RT[21]));
        Kc_save[123] = exp((g_RT[8] + g_RT[12]) - (g_RT[1] + g_RT[22]));
        Kc_save[124] = exp((g_RT[9] + g_RT[3]) - (g_RT[2] + g_RT[16]));
        Kc_save[125] = exp((g_RT[9] + g_RT[0]) - (g_RT[1] + g_RT[10]));
        Kc_save[126] = exp((g_RT[9] + g_RT[5]) - (g_RT[1] + g_RT[17]));
        Kc_save[127] = exp((g_RT[9] + g_RT[10]) - (g_RT[1] + g_RT[22]));
        Kc_save[128] = exp((g_RT[9] + g_RT[12]) - (g_RT[1] + g_RT[23]));
        Kc_save[129] = exp((g_RT[9] + g_RT[13]) - (g_RT[1] + g_RT[24]));
        Kc_save[130] = 1.0 / (refC) * exp((g_RT[9] + g_RT[14]) - (g_RT[27]));
        Kc_save[131] = exp((g_RT[9] + g_RT[15]) - (g_RT[16] + g_RT[14]));
        Kc_save[132] = exp((g_RT[9] + g_RT[17]) - (g_RT[1] + g_RT[28]));
        Kc_save[133] = exp((g_RT[9] + g_RT[27]) - (g_RT[14] + g_RT[22]));
        Kc_save[134] = exp((g_RT[10] + g_RT[3]) - (g_RT[4] + g_RT[16]));
        Kc_save[135] = exp((g_RT[10] + g_RT[0]) - (g_RT[1] + g_RT[12]));
        Kc_save[136] = exp((2 * g_RT[10]) - (g_RT[0] + g_RT[22]));
        Kc_save[137] = exp((g_RT[10] + g_RT[12]) - (g_RT[1] + g_RT[24]));
        Kc_save[138] = exp((g_RT[10] + g_RT[13]) - (2 * g_RT[12]));
        Kc_save[139] = 1.0 / (refC) * exp((g_RT[10] + g_RT[14]) - (g_RT[28]));
        Kc_save[140] = exp((g_RT[10] + g_RT[27]) - (g_RT[23] + g_RT[14]));
        Kc_save[141] = exp((g_RT[11] + g_RT[30]) - (g_RT[10] + g_RT[30]));
        Kc_save[142] = exp((g_RT[11] + g_RT[31]) - (g_RT[10] + g_RT[31]));
        Kc_save[143] = refC * exp((g_RT[11] + g_RT[3]) - (g_RT[1] + g_RT[4] + g_RT[14]));
        Kc_save[144] = exp((g_RT[11] + g_RT[3]) - (g_RT[14] + g_RT[5]));
        Kc_save[145] = exp((g_RT[11] + g_RT[0]) - (g_RT[12] + g_RT[1]));
        Kc_save[146] = 1.0 / (refC) * exp((g_RT[11] + g_RT[5]) - (g_RT[20]));
        Kc_save[147] = exp((g_RT[11] + g_RT[5]) - (g_RT[10] + g_RT[5]));
        Kc_save[148] = exp((g_RT[11] + g_RT[12]) - (g_RT[1] + g_RT[24]));
        Kc_save[149] = exp((g_RT[11] + g_RT[13]) - (2 * g_RT[12]));
        Kc_save[150] = exp((g_RT[11] + g_RT[14]) - (g_RT[10] + g_RT[14]));
        Kc_save[151] = exp((g_RT[11] + g_RT[15]) - (g_RT[10] + g_RT[15]));
        Kc_save[152] = exp((g_RT[11] + g_RT[15]) - (g_RT[14] + g_RT[17]));
        Kc_save[153] = exp((g_RT[11] + g_RT[26]) - (g_RT[12] + g_RT[25]));
        Kc_save[154] = exp((g_RT[12] + g_RT[3]) - (g_RT[2] + g_RT[19]));
        Kc_save[155] = exp((g_RT[12] + g_RT[3]) - (g_RT[4] + g_RT[17]));
        Kc_save[156] = exp((g_RT[12] + g_RT[7]) - (g_RT[6] + g_RT[13]));
        Kc_save[157] = 1.0 / (refC) * exp((2 * g_RT[12]) - (g_RT[26]));
        Kc_save[158] = exp((2 * g_RT[12]) - (g_RT[1] + g_RT[25]));
        Kc_save[159] = exp((g_RT[12] + g_RT[16]) - (g_RT[13] + g_RT[14]));
        Kc_save[160] = exp((g_RT[12] + g_RT[17]) - (g_RT[16] + g_RT[13]));
        Kc_save[161] = exp((g_RT[12] + g_RT[20]) - (g_RT[18] + g_RT[13]));
        Kc_save[162] = exp((g_RT[12] + g_RT[20]) - (g_RT[19] + g_RT[13]));
        Kc_save[163] = exp((g_RT[12] + g_RT[24]) - (g_RT[23] + g_RT[13]));
        Kc_save[164] = exp((g_RT[12] + g_RT[26]) - (g_RT[25] + g_RT[13]));
        Kc_save[165] = refC * exp((g_RT[16] + g_RT[5]) - (g_RT[1] + g_RT[14] + g_RT[5]));
        Kc_save[166] = refC * exp((g_RT[16]) - (g_RT[1] + g_RT[14]));
        Kc_save[167] = exp((g_RT[16] + g_RT[3]) - (g_RT[6] + g_RT[14]));
        Kc_save[168] = exp((g_RT[18] + g_RT[3]) - (g_RT[6] + g_RT[17]));
        Kc_save[169] = exp((g_RT[19] + g_RT[3]) - (g_RT[6] + g_RT[17]));
        Kc_save[170] = exp((g_RT[21] + g_RT[3]) - (g_RT[16] + g_RT[14]));
        Kc_save[171] = exp((g_RT[21] + g_RT[0]) - (g_RT[1] + g_RT[22]));
        Kc_save[172] = exp((g_RT[23] + g_RT[3]) - (g_RT[16] + g_RT[17]));
        Kc_save[173] = refC * exp((g_RT[24]) - (g_RT[0] + g_RT[22]));
        Kc_save[174] = exp((g_RT[25] + g_RT[3]) - (g_RT[6] + g_RT[24]));
        Kc_save[175] = refC * exp((g_RT[27] + g_RT[3]) - (g_RT[4] + 2 * g_RT[14]));
        Kc_save[176] = refC * exp((2 * g_RT[27]) - (2 * g_RT[14] + g_RT[22]));
    }

    /*reaction 1: 2 O + M <=> O2 + M */
    phi_f = sc[2]*sc[2];
    alpha = mixture + 1.4*sc[0] + 14.4*sc[5] + sc[13] + 0.75*sc[14] + 2.6*sc[15] + 2*sc[26] + -0.17*sc[31];
    k_f = alpha * k_f_save[0];
    q_f = phi_f * k_f;
    phi_r = sc[3];
    Kc = Kc_save[0];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[0] = q_f - q_r;

    /*reaction 2: O + H + M <=> OH + M */
    phi_f = sc[2]*sc[1];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = alpha * k_f_save[1];
    q_f = phi_f * k_f;
    phi_r = sc[4];
    Kc = Kc_save[1];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[1] = q_f - q_r;

    /*reaction 3: O + H2 <=> H + OH */
    phi_f = sc[2]*sc[0];
    k_f = k_f_save[2];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[4];
    Kc = Kc_save[2];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[2] = q_f - q_r;

    /*reaction 4: O + HO2 <=> OH + O2 */
    phi_f = sc[2]*sc[6];
    k_f = k_f_save[3];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = Kc_save[3];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[3] = q_f - q_r;

    /*reaction 5: O + H2O2 <=> OH + HO2 */
    phi_f = sc[2]*sc[7];
    k_f = k_f_save[4];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[6];
    Kc = Kc_save[4];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[4] = q_f - q_r;

    /*reaction 6: O + CH <=> H + CO */
    phi_f = sc[2]*sc[9];
    k_f = k_f_save[5];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[5];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[5] = q_f - q_r;

    /*reaction 7: O + CH2 <=> H + HCO */
    phi_f = sc[2]*sc[10];
    k_f = k_f_save[6];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_save[6];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[6] = q_f - q_r;

    /*reaction 8: O + CH2(S) <=> H2 + CO */
    phi_f = sc[2]*sc[11];
    k_f = k_f_save[7];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[14];
    Kc = Kc_save[7];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[7] = q_f - q_r;

    /*reaction 9: O + CH2(S) <=> H + HCO */
    phi_f = sc[2]*sc[11];
    k_f = k_f_save[8];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_save[8];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[8] = q_f - q_r;

    /*reaction 10: O + CH3 <=> H + CH2O */
    phi_f = sc[2]*sc[12];
    k_f = k_f_save[9];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[9];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[9] = q_f - q_r;

    /*reaction 11: O + CH4 <=> OH + CH3 */
    phi_f = sc[2]*sc[13];
    k_f = k_f_save[10];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = Kc_save[10];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[10] = q_f - q_r;

    /*reaction 12: O + CO + M <=> CO2 + M */
    phi_f = sc[2]*sc[14];
    alpha = mixture + sc[0] + 5*sc[3] + 5*sc[5] + sc[13] + 0.5*sc[14] + 2.5*sc[15] + 2*sc[26] + -0.5*sc[31];
    k_f = alpha * k_f_save[11];
    q_f = phi_f * k_f;
    phi_r = sc[15];
    Kc = Kc_save[11];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[11] = q_f - q_r;

    /*reaction 13: O + HCO <=> OH + CO */
    phi_f = sc[2]*sc[16];
    k_f = k_f_save[12];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[14];
    Kc = Kc_save[12];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[12] = q_f - q_r;

    /*reaction 14: O + HCO <=> H + CO2 */
    phi_f = sc[2]*sc[16];
    k_f = k_f_save[13];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[15];
    Kc = Kc_save[13];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[13] = q_f - q_r;

    /*reaction 15: O + CH2O <=> OH + HCO */
    phi_f = sc[2]*sc[17];
    k_f = k_f_save[14];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[16];
    Kc = Kc_save[14];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[14] = q_f - q_r;

    /*reaction 16: O + CH2OH <=> OH + CH2O */
    phi_f = sc[2]*sc[18];
    k_f = k_f_save[15];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[15];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[15] = q_f - q_r;

    /*reaction 17: O + CH3O <=> OH + CH2O */
    phi_f = sc[2]*sc[19];
    k_f = k_f_save[16];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[16];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[16] = q_f - q_r;

    /*reaction 18: O + CH3OH <=> OH + CH2OH */
    phi_f = sc[2]*sc[20];
    k_f = k_f_save[17];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[18];
    Kc = Kc_save[17];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[17] = q_f - q_r;

    /*reaction 19: O + CH3OH <=> OH + CH3O */
    phi_f = sc[2]*sc[20];
    k_f = k_f_save[18];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[19];
    Kc = Kc_save[18];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[18] = q_f - q_r;

    /*reaction 20: O + C2H <=> CH + CO */
    phi_f = sc[2]*sc[21];
    k_f = k_f_save[19];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[14];
    Kc = Kc_save[19];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[19] = q_f - q_r;

    /*reaction 21: O + C2H2 <=> H + HCCO */
    phi_f = sc[2]*sc[22];
    k_f = k_f_save[20];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[27];
    Kc = Kc_save[20];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[20] = q_f - q_r;

    /*reaction 22: O + C2H2 <=> OH + C2H */
    phi_f = sc[2]*sc[22];
    k_f = k_f_save[21];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[21];
    Kc = Kc_save[21];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[21] = q_f - q_r;

    /*reaction 23: O + C2H2 <=> CO + CH2 */
    phi_f = sc[2]*sc[22];
    k_f = k_f_save[22];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[10];
    Kc = Kc_save[22];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[22] = q_f - q_r;

    /*reaction 24: O + C2H3 <=> H + CH2CO */
    phi_f = sc[2]*sc[23];
    k_f = k_f_save[23];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[23];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[23] = q_f - q_r;

    /*reaction 25: O + C2H4 <=> CH3 + HCO */
    phi_f = sc[2]*sc[24];
    k_f = k_f_save[24];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[16];
    Kc = Kc_save[24];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[24] = q_f - q_r;

    /*reaction 26: O + C2H5 <=> CH3 + CH2O */
    phi_f = sc[2]*sc[25];
    k_f = k_f_save[25];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[17];
    Kc = Kc_save[25];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[25] = q_f - q_r;

    /*reaction 27: O + C2H6 <=> OH + C2H5 */
    phi_f = sc[2]*sc[26];
    k_f = k_f_save[26];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[25];
    Kc = Kc_save[26];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[26] = q_f - q_r;

    /*reaction 28: O + HCCO <=> H + 2 CO */
    phi_f = sc[2]*sc[27];
    k_f = k_f_save[27];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14]*sc[14];
    Kc = Kc_save[27];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[27] = q_f - q_r;

    /*reaction 29: O + CH2CO <=> OH + HCCO */
    phi_f = sc[2]*sc[28];
    k_f = k_f_save[28];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[27];
    Kc = Kc_save[28];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[28] = q_f - q_r;

    /*reaction 30: O + CH2CO <=> CH2 + CO2 */
    phi_f = sc[2]*sc[28];
    k_f = k_f_save[29];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[15];
    Kc = Kc_save[29];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[29] = q_f - q_r;

    /*reaction 31: O2 + CO <=> O + CO2 */
    phi_f = sc[3]*sc[14];
    k_f = k_f_save[30];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[15];
    Kc = Kc_save[30];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[30] = q_f - q_r;

    /*reaction 32: O2 + CH2O <=> HO2 + HCO */
    phi_f = sc[3]*sc[17];
    k_f = k_f_save[31];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[16];
    Kc = Kc_save[31];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[31] = q_f - q_r;

    /*reaction 33: H + O2 + M <=> HO2 + M */
    phi_f = sc[1]*sc[3];
    alpha = mixture + -1*sc[3] + -1*sc[5] + -0.25*sc[14] + 0.5*sc[15] + 0.5*sc[26] + -1*sc[30] + -1*sc[31];
    k_f = alpha * k_f_save[32];
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = Kc_save[32];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[32] = q_f - q_r;

    /*reaction 34: H + 2 O2 <=> HO2 + O2 */
    phi_f = sc[1]*sc[3]*sc[3];
    k_f = k_f_save[33];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = Kc_save[33];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[33] = q_f - q_r;

    /*reaction 35: H + O2 + H2O <=> HO2 + H2O */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = k_f_save[34];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = Kc_save[34];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[34] = q_f - q_r;

    /*reaction 36: H + O2 + N2 <=> HO2 + N2 */
    phi_f = sc[1]*sc[3]*sc[30];
    k_f = k_f_save[35];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[30];
    Kc = Kc_save[35];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[35] = q_f - q_r;

    /*reaction 37: H + O2 + AR <=> HO2 + AR */
    phi_f = sc[1]*sc[3]*sc[31];
    k_f = k_f_save[36];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[31];
    Kc = Kc_save[36];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[36] = q_f - q_r;

    /*reaction 38: H + O2 <=> O + OH */
    phi_f = sc[1]*sc[3];
    k_f = k_f_save[37];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = Kc_save[37];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[37] = q_f - q_r;

    /*reaction 39: 2 H + M <=> H2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[0] + -1*sc[5] + sc[13] + -1*sc[15] + 2*sc[26] + -0.37*sc[31];
    k_f = alpha * k_f_save[38];
    q_f = phi_f * k_f;
    phi_r = sc[0];
    Kc = Kc_save[38];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[38] = q_f - q_r;

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    phi_f = sc[1]*sc[1]*sc[0];
    k_f = k_f_save[39];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[0];
    Kc = Kc_save[39];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[39] = q_f - q_r;

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    phi_f = sc[1]*sc[1]*sc[5];
    k_f = k_f_save[40];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[5];
    Kc = Kc_save[40];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[40] = q_f - q_r;

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    phi_f = sc[1]*sc[1]*sc[15];
    k_f = k_f_save[41];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[15];
    Kc = Kc_save[41];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[41] = q_f - q_r;

    /*reaction 43: H + OH + M <=> H2O + M */
    phi_f = sc[1]*sc[4];
    alpha = mixture + -0.27*sc[0] + 2.65*sc[5] + sc[13] + 2*sc[26] + -0.62*sc[31];
    k_f = alpha * k_f_save[42];
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = Kc_save[42];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[42] = q_f - q_r;

    /*reaction 44: H + HO2 <=> O + H2O */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[43];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = Kc_save[43];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[43] = q_f - q_r;

    /*reaction 45: H + HO2 <=> O2 + H2 */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[44];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = Kc_save[44];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[44] = q_f - q_r;

    /*reaction 46: H + HO2 <=> 2 OH */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[45];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = Kc_save[45];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[45] = q_f - q_r;

    /*reaction 47: H + H2O2 <=> HO2 + H2 */
    phi_f = sc[1]*sc[7];
    k_f = k_f_save[46];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = Kc_save[46];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[46] = q_f - q_r;

    /*reaction 48: H + H2O2 <=> OH + H2O */
    phi_f = sc[1]*sc[7];
    k_f = k_f_save[47];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[5];
    Kc = Kc_save[47];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[47] = q_f - q_r;

    /*reaction 49: H + CH <=> C + H2 */
    phi_f = sc[1]*sc[9];
    k_f = k_f_save[48];
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = Kc_save[48];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[48] = q_f - q_r;

    /*reaction 50: H + CH2 (+M) <=> CH3 (+M) */
    phi_f = sc[1]*sc[10];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[49];
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
    phi_r = sc[12];
    Kc = Kc_save[49];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[49] = q_f - q_r;

    /*reaction 51: H + CH2(S) <=> CH + H2 */
    phi_f = sc[1]*sc[11];
    k_f = k_f_save[50];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = Kc_save[50];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[50] = q_f - q_r;

    /*reaction 52: H + CH3 (+M) <=> CH4 (+M) */
    phi_f = sc[1]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[51];
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
    phi_r = sc[13];
    Kc = Kc_save[51];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[51] = q_f - q_r;

    /*reaction 53: H + CH4 <=> CH3 + H2 */
    phi_f = sc[1]*sc[13];
    k_f = k_f_save[52];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[0];
    Kc = Kc_save[52];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[52] = q_f - q_r;

    /*reaction 54: H + HCO (+M) <=> CH2O (+M) */
    phi_f = sc[1]*sc[16];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[53];
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
    phi_r = sc[17];
    Kc = Kc_save[53];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[53] = q_f - q_r;

    /*reaction 55: H + HCO <=> H2 + CO */
    phi_f = sc[1]*sc[16];
    k_f = k_f_save[54];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[14];
    Kc = Kc_save[54];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[54] = q_f - q_r;

    /*reaction 56: H + CH2O (+M) <=> CH2OH (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[55];
    redP = 1e-12 * alpha / k_f * 1.27e+32*exp(-4.82*tc[0]-3286.0048277048194905*invT);
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
    Kc = Kc_save[55];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[55] = q_f - q_r;

    /*reaction 57: H + CH2O (+M) <=> CH3O (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[56];
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
    phi_r = sc[19];
    Kc = Kc_save[56];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[56] = q_f - q_r;

    /*reaction 58: H + CH2O <=> HCO + H2 */
    phi_f = sc[1]*sc[17];
    k_f = k_f_save[57];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[0];
    Kc = Kc_save[57];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[57] = q_f - q_r;

    /*reaction 59: H + CH2OH (+M) <=> CH3OH (+M) */
    phi_f = sc[1]*sc[18];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[58];
    redP = 1e-12 * alpha / k_f * 3e+31*exp(-4.8*tc[0]-1660.6149971555748834*invT);
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
    Kc = Kc_save[58];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[58] = q_f - q_r;

    /*reaction 60: H + CH2OH <=> H2 + CH2O */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[59];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[17];
    Kc = Kc_save[59];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[59] = q_f - q_r;

    /*reaction 61: H + CH2OH <=> OH + CH3 */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[60];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = Kc_save[60];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[60] = q_f - q_r;

    /*reaction 62: H + CH2OH <=> CH2(S) + H2O */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[61];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[61];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[61] = q_f - q_r;

    /*reaction 63: H + CH3O (+M) <=> CH3OH (+M) */
    phi_f = sc[1]*sc[19];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[62];
    redP = 1e-12 * alpha / k_f * 8.6e+28*exp(-4*tc[0]-1522.2304140592771091*invT);
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
    Kc = Kc_save[62];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[62] = q_f - q_r;

    /*reaction 64: H + CH3O <=> H + CH2OH */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[63];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[18];
    Kc = Kc_save[63];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[63] = q_f - q_r;

    /*reaction 65: H + CH3O <=> H2 + CH2O */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[64];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[17];
    Kc = Kc_save[64];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[64] = q_f - q_r;

    /*reaction 66: H + CH3O <=> OH + CH3 */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[65];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = Kc_save[65];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[65] = q_f - q_r;

    /*reaction 67: H + CH3O <=> CH2(S) + H2O */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[66];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[66];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[66] = q_f - q_r;

    /*reaction 68: H + CH3OH <=> CH2OH + H2 */
    phi_f = sc[1]*sc[20];
    k_f = k_f_save[67];
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[0];
    Kc = Kc_save[67];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[67] = q_f - q_r;

    /*reaction 69: H + CH3OH <=> CH3O + H2 */
    phi_f = sc[1]*sc[20];
    k_f = k_f_save[68];
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[0];
    Kc = Kc_save[68];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[68] = q_f - q_r;

    /*reaction 70: H + C2H (+M) <=> C2H2 (+M) */
    phi_f = sc[1]*sc[21];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[69];
    redP = 1e-12 * alpha / k_f * 3.75e+33*exp(-4.8*tc[0]-956.11166502896742259*invT);
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
    Kc = Kc_save[69];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[69] = q_f - q_r;

    /*reaction 71: H + C2H2 (+M) <=> C2H3 (+M) */
    phi_f = sc[1]*sc[22];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[70];
    redP = 1e-12 * alpha / k_f * 3.8e+40*exp(-7.27*tc[0]-3633.2243271100760467*invT);
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
    Kc = Kc_save[70];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[70] = q_f - q_r;

    /*reaction 72: H + C2H3 (+M) <=> C2H4 (+M) */
    phi_f = sc[1]*sc[23];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[71];
    redP = 1e-12 * alpha / k_f * 1.4e+30*exp(-3.86*tc[0]-1670.6793304716693456*invT);
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
    Kc = Kc_save[71];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[71] = q_f - q_r;

    /*reaction 73: H + C2H3 <=> H2 + C2H2 */
    phi_f = sc[1]*sc[23];
    k_f = k_f_save[72];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = Kc_save[72];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[72] = q_f - q_r;

    /*reaction 74: H + C2H4 (+M) <=> C2H5 (+M) */
    phi_f = sc[1]*sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[73];
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
    phi_r = sc[25];
    Kc = Kc_save[73];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[73] = q_f - q_r;

    /*reaction 75: H + C2H4 <=> C2H3 + H2 */
    phi_f = sc[1]*sc[24];
    k_f = k_f_save[74];
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[0];
    Kc = Kc_save[74];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[74] = q_f - q_r;

    /*reaction 76: H + C2H5 (+M) <=> C2H6 (+M) */
    phi_f = sc[1]*sc[25];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[75];
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
    phi_r = sc[26];
    Kc = Kc_save[75];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[75] = q_f - q_r;

    /*reaction 77: H + C2H5 <=> H2 + C2H4 */
    phi_f = sc[1]*sc[25];
    k_f = k_f_save[76];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[24];
    Kc = Kc_save[76];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[76] = q_f - q_r;

    /*reaction 78: H + C2H6 <=> C2H5 + H2 */
    phi_f = sc[1]*sc[26];
    k_f = k_f_save[77];
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[0];
    Kc = Kc_save[77];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[77] = q_f - q_r;

    /*reaction 79: H + HCCO <=> CH2(S) + CO */
    phi_f = sc[1]*sc[27];
    k_f = k_f_save[78];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[14];
    Kc = Kc_save[78];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[78] = q_f - q_r;

    /*reaction 80: H + CH2CO <=> HCCO + H2 */
    phi_f = sc[1]*sc[28];
    k_f = k_f_save[79];
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[0];
    Kc = Kc_save[79];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[79] = q_f - q_r;

    /*reaction 81: H + CH2CO <=> CH3 + CO */
    phi_f = sc[1]*sc[28];
    k_f = k_f_save[80];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[14];
    Kc = Kc_save[80];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[80] = q_f - q_r;

    /*reaction 82: H + HCCOH <=> H + CH2CO */
    phi_f = sc[1]*sc[29];
    k_f = k_f_save[81];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[81];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[81] = q_f - q_r;

    /*reaction 83: H2 + CO (+M) <=> CH2O (+M) */
    phi_f = sc[0]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[82];
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
    phi_r = sc[17];
    Kc = Kc_save[82];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[82] = q_f - q_r;

    /*reaction 84: OH + H2 <=> H + H2O */
    phi_f = sc[4]*sc[0];
    k_f = k_f_save[83];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    Kc = Kc_save[83];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[83] = q_f - q_r;

    /*reaction 85: 2 OH (+M) <=> H2O2 (+M) */
    phi_f = sc[4]*sc[4];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[84];
    redP = 1e-12 * alpha / k_f * 2.3e+18*exp(-0.9*tc[0]+855.46833186802348337*invT);
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
    Kc = Kc_save[84];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[84] = q_f - q_r;

    /*reaction 86: 2 OH <=> O + H2O */
    phi_f = sc[4]*sc[4];
    k_f = k_f_save[85];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = Kc_save[85];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[85] = q_f - q_r;

    /*reaction 87: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = k_f_save[86];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = Kc_save[86];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[86] = q_f - q_r;

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    phi_f = sc[4]*sc[7];
    k_f = k_f_save[87];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = Kc_save[87];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[87] = q_f - q_r;

    /*reaction 89: OH + H2O2 <=> HO2 + H2O */
    phi_f = sc[4]*sc[7];
    k_f = k_f_save[88];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = Kc_save[88];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[88] = q_f - q_r;

    /*reaction 90: OH + C <=> H + CO */
    phi_f = sc[4]*sc[8];
    k_f = k_f_save[89];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[89];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[89] = q_f - q_r;

    /*reaction 91: OH + CH <=> H + HCO */
    phi_f = sc[4]*sc[9];
    k_f = k_f_save[90];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_save[90];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[90] = q_f - q_r;

    /*reaction 92: OH + CH2 <=> H + CH2O */
    phi_f = sc[4]*sc[10];
    k_f = k_f_save[91];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[91];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[91] = q_f - q_r;

    /*reaction 93: OH + CH2 <=> CH + H2O */
    phi_f = sc[4]*sc[10];
    k_f = k_f_save[92];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[5];
    Kc = Kc_save[92];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[92] = q_f - q_r;

    /*reaction 94: OH + CH2(S) <=> H + CH2O */
    phi_f = sc[4]*sc[11];
    k_f = k_f_save[93];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[93];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[93] = q_f - q_r;

    /*reaction 95: OH + CH3 (+M) <=> CH3OH (+M) */
    phi_f = sc[4]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[94];
    redP = 1e-12 * alpha / k_f * 2.7e+38*exp(-6.3*tc[0]-1559.9716639946311716*invT);
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
    Kc = Kc_save[94];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[94] = q_f - q_r;

    /*reaction 96: OH + CH3 <=> CH2 + H2O */
    phi_f = sc[4]*sc[12];
    k_f = k_f_save[95];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[5];
    Kc = Kc_save[95];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[95] = q_f - q_r;

    /*reaction 97: OH + CH3 <=> CH2(S) + H2O */
    phi_f = sc[4]*sc[12];
    k_f = k_f_save[96];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[96];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[96] = q_f - q_r;

    /*reaction 98: OH + CH4 <=> CH3 + H2O */
    phi_f = sc[4]*sc[13];
    k_f = k_f_save[97];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[5];
    Kc = Kc_save[97];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[97] = q_f - q_r;

    /*reaction 99: OH + CO <=> H + CO2 */
    phi_f = sc[4]*sc[14];
    k_f = k_f_save[98];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[15];
    Kc = Kc_save[98];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[98] = q_f - q_r;

    /*reaction 100: OH + HCO <=> H2O + CO */
    phi_f = sc[4]*sc[16];
    k_f = k_f_save[99];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[14];
    Kc = Kc_save[99];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[99] = q_f - q_r;

    /*reaction 101: OH + CH2O <=> HCO + H2O */
    phi_f = sc[4]*sc[17];
    k_f = k_f_save[100];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[5];
    Kc = Kc_save[100];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[100] = q_f - q_r;

    /*reaction 102: OH + CH2OH <=> H2O + CH2O */
    phi_f = sc[4]*sc[18];
    k_f = k_f_save[101];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[17];
    Kc = Kc_save[101];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[101] = q_f - q_r;

    /*reaction 103: OH + CH3O <=> H2O + CH2O */
    phi_f = sc[4]*sc[19];
    k_f = k_f_save[102];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[17];
    Kc = Kc_save[102];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[102] = q_f - q_r;

    /*reaction 104: OH + CH3OH <=> CH2OH + H2O */
    phi_f = sc[4]*sc[20];
    k_f = k_f_save[103];
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[5];
    Kc = Kc_save[103];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[103] = q_f - q_r;

    /*reaction 105: OH + CH3OH <=> CH3O + H2O */
    phi_f = sc[4]*sc[20];
    k_f = k_f_save[104];
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[5];
    Kc = Kc_save[104];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[104] = q_f - q_r;

    /*reaction 106: OH + C2H <=> H + HCCO */
    phi_f = sc[4]*sc[21];
    k_f = k_f_save[105];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[27];
    Kc = Kc_save[105];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[105] = q_f - q_r;

    /*reaction 107: OH + C2H2 <=> H + CH2CO */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[106];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[106];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[106] = q_f - q_r;

    /*reaction 108: OH + C2H2 <=> H + HCCOH */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[107];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[29];
    Kc = Kc_save[107];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[107] = q_f - q_r;

    /*reaction 109: OH + C2H2 <=> C2H + H2O */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[108];
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[5];
    Kc = Kc_save[108];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[108] = q_f - q_r;

    /*reaction 110: OH + C2H2 <=> CH3 + CO */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[109];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[14];
    Kc = Kc_save[109];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[109] = q_f - q_r;

    /*reaction 111: OH + C2H3 <=> H2O + C2H2 */
    phi_f = sc[4]*sc[23];
    k_f = k_f_save[110];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[22];
    Kc = Kc_save[110];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[110] = q_f - q_r;

    /*reaction 112: OH + C2H4 <=> C2H3 + H2O */
    phi_f = sc[4]*sc[24];
    k_f = k_f_save[111];
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[5];
    Kc = Kc_save[111];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[111] = q_f - q_r;

    /*reaction 113: OH + C2H6 <=> C2H5 + H2O */
    phi_f = sc[4]*sc[26];
    k_f = k_f_save[112];
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[5];
    Kc = Kc_save[112];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[112] = q_f - q_r;

    /*reaction 114: OH + CH2CO <=> HCCO + H2O */
    phi_f = sc[4]*sc[28];
    k_f = k_f_save[113];
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[5];
    Kc = Kc_save[113];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[113] = q_f - q_r;

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    phi_f = sc[6]*sc[6];
    k_f = k_f_save[114];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = Kc_save[114];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[114] = q_f - q_r;

    /*reaction 116: 2 HO2 <=> O2 + H2O2 */
    phi_f = sc[6]*sc[6];
    k_f = k_f_save[115];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = Kc_save[115];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[115] = q_f - q_r;

    /*reaction 117: HO2 + CH2 <=> OH + CH2O */
    phi_f = sc[6]*sc[10];
    k_f = k_f_save[116];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[116];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[116] = q_f - q_r;

    /*reaction 118: HO2 + CH3 <=> O2 + CH4 */
    phi_f = sc[6]*sc[12];
    k_f = k_f_save[117];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[13];
    Kc = Kc_save[117];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[117] = q_f - q_r;

    /*reaction 119: HO2 + CH3 <=> OH + CH3O */
    phi_f = sc[6]*sc[12];
    k_f = k_f_save[118];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[19];
    Kc = Kc_save[118];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[118] = q_f - q_r;

    /*reaction 120: HO2 + CO <=> OH + CO2 */
    phi_f = sc[6]*sc[14];
    k_f = k_f_save[119];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[15];
    Kc = Kc_save[119];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[119] = q_f - q_r;

    /*reaction 121: HO2 + CH2O <=> HCO + H2O2 */
    phi_f = sc[6]*sc[17];
    k_f = k_f_save[120];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[7];
    Kc = Kc_save[120];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[120] = q_f - q_r;

    /*reaction 122: C + O2 <=> O + CO */
    phi_f = sc[8]*sc[3];
    k_f = k_f_save[121];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[14];
    Kc = Kc_save[121];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[121] = q_f - q_r;

    /*reaction 123: C + CH2 <=> H + C2H */
    phi_f = sc[8]*sc[10];
    k_f = k_f_save[122];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[21];
    Kc = Kc_save[122];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[122] = q_f - q_r;

    /*reaction 124: C + CH3 <=> H + C2H2 */
    phi_f = sc[8]*sc[12];
    k_f = k_f_save[123];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = Kc_save[123];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[123] = q_f - q_r;

    /*reaction 125: CH + O2 <=> O + HCO */
    phi_f = sc[9]*sc[3];
    k_f = k_f_save[124];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[16];
    Kc = Kc_save[124];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[124] = q_f - q_r;

    /*reaction 126: CH + H2 <=> H + CH2 */
    phi_f = sc[9]*sc[0];
    k_f = k_f_save[125];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[10];
    Kc = Kc_save[125];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[125] = q_f - q_r;

    /*reaction 127: CH + H2O <=> H + CH2O */
    phi_f = sc[9]*sc[5];
    k_f = k_f_save[126];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[126];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[126] = q_f - q_r;

    /*reaction 128: CH + CH2 <=> H + C2H2 */
    phi_f = sc[9]*sc[10];
    k_f = k_f_save[127];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = Kc_save[127];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[127] = q_f - q_r;

    /*reaction 129: CH + CH3 <=> H + C2H3 */
    phi_f = sc[9]*sc[12];
    k_f = k_f_save[128];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[23];
    Kc = Kc_save[128];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[128] = q_f - q_r;

    /*reaction 130: CH + CH4 <=> H + C2H4 */
    phi_f = sc[9]*sc[13];
    k_f = k_f_save[129];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = Kc_save[129];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[129] = q_f - q_r;

    /*reaction 131: CH + CO (+M) <=> HCCO (+M) */
    phi_f = sc[9]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[130];
    redP = 1e-12 * alpha / k_f * 2.69e+28*exp(-3.74*tc[0]-974.227464997937318*invT);
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
    Kc = Kc_save[130];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[130] = q_f - q_r;

    /*reaction 132: CH + CO2 <=> HCO + CO */
    phi_f = sc[9]*sc[15];
    k_f = k_f_save[131];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[14];
    Kc = Kc_save[131];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[131] = q_f - q_r;

    /*reaction 133: CH + CH2O <=> H + CH2CO */
    phi_f = sc[9]*sc[17];
    k_f = k_f_save[132];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[132];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[132] = q_f - q_r;

    /*reaction 134: CH + HCCO <=> CO + C2H2 */
    phi_f = sc[9]*sc[27];
    k_f = k_f_save[133];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[22];
    Kc = Kc_save[133];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[133] = q_f - q_r;

    /*reaction 135: CH2 + O2 <=> OH + HCO */
    phi_f = sc[10]*sc[3];
    k_f = k_f_save[134];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[16];
    Kc = Kc_save[134];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[134] = q_f - q_r;

    /*reaction 136: CH2 + H2 <=> H + CH3 */
    phi_f = sc[10]*sc[0];
    k_f = k_f_save[135];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = Kc_save[135];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[135] = q_f - q_r;

    /*reaction 137: 2 CH2 <=> H2 + C2H2 */
    phi_f = sc[10]*sc[10];
    k_f = k_f_save[136];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = Kc_save[136];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[136] = q_f - q_r;

    /*reaction 138: CH2 + CH3 <=> H + C2H4 */
    phi_f = sc[10]*sc[12];
    k_f = k_f_save[137];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = Kc_save[137];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[137] = q_f - q_r;

    /*reaction 139: CH2 + CH4 <=> 2 CH3 */
    phi_f = sc[10]*sc[13];
    k_f = k_f_save[138];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[12];
    Kc = Kc_save[138];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[138] = q_f - q_r;

    /*reaction 140: CH2 + CO (+M) <=> CH2CO (+M) */
    phi_f = sc[10]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[139];
    redP = 1e-12 * alpha / k_f * 2.69e+33*exp(-5.11*tc[0]-3570.3222438844863973*invT);
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
    Kc = Kc_save[139];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[139] = q_f - q_r;

    /*reaction 141: CH2 + HCCO <=> C2H3 + CO */
    phi_f = sc[10]*sc[27];
    k_f = k_f_save[140];
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[14];
    Kc = Kc_save[140];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[140] = q_f - q_r;

    /*reaction 142: CH2(S) + N2 <=> CH2 + N2 */
    phi_f = sc[11]*sc[30];
    k_f = k_f_save[141];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[30];
    Kc = Kc_save[141];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[141] = q_f - q_r;

    /*reaction 143: CH2(S) + AR <=> CH2 + AR */
    phi_f = sc[11]*sc[31];
    k_f = k_f_save[142];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[31];
    Kc = Kc_save[142];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[142] = q_f - q_r;

    /*reaction 144: CH2(S) + O2 <=> H + OH + CO */
    phi_f = sc[11]*sc[3];
    k_f = k_f_save[143];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[4]*sc[14];
    Kc = Kc_save[143];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[143] = q_f - q_r;

    /*reaction 145: CH2(S) + O2 <=> CO + H2O */
    phi_f = sc[11]*sc[3];
    k_f = k_f_save[144];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[5];
    Kc = Kc_save[144];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[144] = q_f - q_r;

    /*reaction 146: CH2(S) + H2 <=> CH3 + H */
    phi_f = sc[11]*sc[0];
    k_f = k_f_save[145];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[1];
    Kc = Kc_save[145];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[145] = q_f - q_r;

    /*reaction 147: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    phi_f = sc[11]*sc[5];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[146];
    redP = 1e-12 * alpha / k_f * 2.7e+38*exp(-6.3*tc[0]-1559.9716639946311716*invT);
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
    Kc = Kc_save[146];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[146] = q_f - q_r;

    /*reaction 148: CH2(S) + H2O <=> CH2 + H2O */
    phi_f = sc[11]*sc[5];
    k_f = k_f_save[147];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[5];
    Kc = Kc_save[147];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[147] = q_f - q_r;

    /*reaction 149: CH2(S) + CH3 <=> H + C2H4 */
    phi_f = sc[11]*sc[12];
    k_f = k_f_save[148];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = Kc_save[148];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[148] = q_f - q_r;

    /*reaction 150: CH2(S) + CH4 <=> 2 CH3 */
    phi_f = sc[11]*sc[13];
    k_f = k_f_save[149];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[12];
    Kc = Kc_save[149];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[149] = q_f - q_r;

    /*reaction 151: CH2(S) + CO <=> CH2 + CO */
    phi_f = sc[11]*sc[14];
    k_f = k_f_save[150];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[14];
    Kc = Kc_save[150];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[150] = q_f - q_r;

    /*reaction 152: CH2(S) + CO2 <=> CH2 + CO2 */
    phi_f = sc[11]*sc[15];
    k_f = k_f_save[151];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[15];
    Kc = Kc_save[151];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[151] = q_f - q_r;

    /*reaction 153: CH2(S) + CO2 <=> CO + CH2O */
    phi_f = sc[11]*sc[15];
    k_f = k_f_save[152];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[17];
    Kc = Kc_save[152];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[152] = q_f - q_r;

    /*reaction 154: CH2(S) + C2H6 <=> CH3 + C2H5 */
    phi_f = sc[11]*sc[26];
    k_f = k_f_save[153];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[25];
    Kc = Kc_save[153];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[153] = q_f - q_r;

    /*reaction 155: CH3 + O2 <=> O + CH3O */
    phi_f = sc[12]*sc[3];
    k_f = k_f_save[154];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[19];
    Kc = Kc_save[154];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[154] = q_f - q_r;

    /*reaction 156: CH3 + O2 <=> OH + CH2O */
    phi_f = sc[12]*sc[3];
    k_f = k_f_save[155];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[155];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[155] = q_f - q_r;

    /*reaction 157: CH3 + H2O2 <=> HO2 + CH4 */
    phi_f = sc[12]*sc[7];
    k_f = k_f_save[156];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = Kc_save[156];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[156] = q_f - q_r;

    /*reaction 158: 2 CH3 (+M) <=> C2H6 (+M) */
    phi_f = sc[12]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[157];
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
    phi_r = sc[26];
    Kc = Kc_save[157];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[157] = q_f - q_r;

    /*reaction 159: 2 CH3 <=> H + C2H5 */
    phi_f = sc[12]*sc[12];
    k_f = k_f_save[158];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[25];
    Kc = Kc_save[158];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[158] = q_f - q_r;

    /*reaction 160: CH3 + HCO <=> CH4 + CO */
    phi_f = sc[12]*sc[16];
    k_f = k_f_save[159];
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[14];
    Kc = Kc_save[159];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[159] = q_f - q_r;

    /*reaction 161: CH3 + CH2O <=> HCO + CH4 */
    phi_f = sc[12]*sc[17];
    k_f = k_f_save[160];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[13];
    Kc = Kc_save[160];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[160] = q_f - q_r;

    /*reaction 162: CH3 + CH3OH <=> CH2OH + CH4 */
    phi_f = sc[12]*sc[20];
    k_f = k_f_save[161];
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[13];
    Kc = Kc_save[161];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[161] = q_f - q_r;

    /*reaction 163: CH3 + CH3OH <=> CH3O + CH4 */
    phi_f = sc[12]*sc[20];
    k_f = k_f_save[162];
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[13];
    Kc = Kc_save[162];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[162] = q_f - q_r;

    /*reaction 164: CH3 + C2H4 <=> C2H3 + CH4 */
    phi_f = sc[12]*sc[24];
    k_f = k_f_save[163];
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[13];
    Kc = Kc_save[163];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[163] = q_f - q_r;

    /*reaction 165: CH3 + C2H6 <=> C2H5 + CH4 */
    phi_f = sc[12]*sc[26];
    k_f = k_f_save[164];
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[13];
    Kc = Kc_save[164];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[164] = q_f - q_r;

    /*reaction 166: HCO + H2O <=> H + CO + H2O */
    phi_f = sc[16]*sc[5];
    k_f = k_f_save[165];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14]*sc[5];
    Kc = Kc_save[165];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[165] = q_f - q_r;

    /*reaction 167: HCO + M <=> H + CO + M */
    phi_f = sc[16];
    alpha = mixture + sc[0] + -1*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = alpha * k_f_save[166];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[166];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[166] = q_f - q_r;

    /*reaction 168: HCO + O2 <=> HO2 + CO */
    phi_f = sc[16]*sc[3];
    k_f = k_f_save[167];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[14];
    Kc = Kc_save[167];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[167] = q_f - q_r;

    /*reaction 169: CH2OH + O2 <=> HO2 + CH2O */
    phi_f = sc[18]*sc[3];
    k_f = k_f_save[168];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = Kc_save[168];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[168] = q_f - q_r;

    /*reaction 170: CH3O + O2 <=> HO2 + CH2O */
    phi_f = sc[19]*sc[3];
    k_f = k_f_save[169];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = Kc_save[169];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[169] = q_f - q_r;

    /*reaction 171: C2H + O2 <=> HCO + CO */
    phi_f = sc[21]*sc[3];
    k_f = k_f_save[170];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[14];
    Kc = Kc_save[170];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[170] = q_f - q_r;

    /*reaction 172: C2H + H2 <=> H + C2H2 */
    phi_f = sc[21]*sc[0];
    k_f = k_f_save[171];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = Kc_save[171];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[171] = q_f - q_r;

    /*reaction 173: C2H3 + O2 <=> HCO + CH2O */
    phi_f = sc[23]*sc[3];
    k_f = k_f_save[172];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[17];
    Kc = Kc_save[172];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[172] = q_f - q_r;

    /*reaction 174: C2H4 (+M) <=> H2 + C2H2 (+M) */
    phi_f = sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[173];
    redP = 1e-6 * alpha / k_f * 7e+50*exp(-9.31*tc[0]-50251.216247259304509*invT);
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
    Kc = Kc_save[173];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[173] = q_f - q_r;

    /*reaction 175: C2H5 + O2 <=> HO2 + C2H4 */
    phi_f = sc[25]*sc[3];
    k_f = k_f_save[174];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[24];
    Kc = Kc_save[174];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[174] = q_f - q_r;

    /*reaction 176: HCCO + O2 <=> OH + 2 CO */
    phi_f = sc[27]*sc[3];
    k_f = k_f_save[175];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[14]*sc[14];
    Kc = Kc_save[175];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[175] = q_f - q_r;

    /*reaction 177: 2 HCCO <=> 2 CO + C2H2 */
    phi_f = sc[27]*sc[27];
    k_f = k_f_save[176];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[14]*sc[22];
    Kc = Kc_save[176];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[176] = q_f - q_r;

    return;
}


/*compute the progress rate for each reaction */
void progressRateFR(double * q_f, double * q_r, double * sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[32];                /*Gibbs free energy */
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
    for (id = 0; id < 32; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    if (T != T_save)
    {
        T_save = T;

        k_f_save[0] = 1e-12 * 1.2e+17*exp(-1*tc[0]);
        k_f_save[1] = 1e-12 * 5e+17*exp(-1*tc[0]);
        k_f_save[2] = 1e-06 * 50000*exp(2.67*tc[0]-3165.2328279116868543*invT);
        k_f_save[3] = 1e-06 * 2e+13;
        k_f_save[4] = 1e-06 * 9.63e+06*exp(2*tc[0]-2012.866663218878557*invT);
        k_f_save[5] = 1e-06 * 5.7e+13;
        k_f_save[6] = 1e-06 * 8e+13;
        k_f_save[7] = 1e-06 * 1.5e+13;
        k_f_save[8] = 1e-06 * 1.5e+13;
        k_f_save[9] = 1e-06 * 8.43e+13;
        k_f_save[10] = 1e-06 * 1.02e+09*exp(1.5*tc[0]-4327.6633259205891591*invT);
        k_f_save[11] = 1e-12 * 6.02e+14*exp(-1509.6499974141590883*invT);
        k_f_save[12] = 1e-06 * 3e+13;
        k_f_save[13] = 1e-06 * 3e+13;
        k_f_save[14] = 1e-06 * 3.9e+13*exp(-1781.3869969487077469*invT);
        k_f_save[15] = 1e-06 * 1e+13;
        k_f_save[16] = 1e-06 * 1e+13;
        k_f_save[17] = 1e-06 * 388000*exp(2.5*tc[0]-1559.9716639946311716*invT);
        k_f_save[18] = 1e-06 * 130000*exp(2.5*tc[0]-2516.0833290235987079*invT);
        k_f_save[19] = 1e-06 * 5e+13;
        k_f_save[20] = 1e-06 * 1.02e+07*exp(2*tc[0]-956.11166502896742259*invT);
        k_f_save[21] = 1e-06 * 4.6e+19*exp(-1.41*tc[0]-14568.122475046635373*invT);
        k_f_save[22] = 1e-06 * 1.02e+07*exp(2*tc[0]-956.11166502896742259*invT);
        k_f_save[23] = 1e-06 * 3e+13;
        k_f_save[24] = 1e-06 * 1.92e+07*exp(1.83*tc[0]-110.7076664770383303*invT);
        k_f_save[25] = 1e-06 * 1.32e+14;
        k_f_save[26] = 1e-06 * 8.98e+07*exp(1.92*tc[0]-2863.3028284288548093*invT);
        k_f_save[27] = 1e-06 * 1e+14;
        k_f_save[28] = 1e-06 * 1e+13*exp(-4025.733326437757114*invT);
        k_f_save[29] = 1e-06 * 1.75e+12*exp(-679.34249883637164658*invT);
        k_f_save[30] = 1e-06 * 2.5e+12*exp(-24053.756625465601246*invT);
        k_f_save[31] = 1e-06 * 1e+14*exp(-20128.666632188789663*invT);
        k_f_save[32] = 1e-12 * 2.8e+18*exp(-0.86*tc[0]);
        k_f_save[33] = 1e-12 * 3e+20*exp(-1.72*tc[0]);
        k_f_save[34] = 1e-12 * 9.38e+18*exp(-0.76*tc[0]);
        k_f_save[35] = 1e-12 * 3.75e+20*exp(-1.72*tc[0]);
        k_f_save[36] = 1e-12 * 7e+17*exp(-0.8*tc[0]);
        k_f_save[37] = 1e-06 * 8.3e+13*exp(-7252.8618042434254676*invT);
        k_f_save[38] = 1e-12 * 1e+18*exp(-1*tc[0]);
        k_f_save[39] = 1e-12 * 9e+16*exp(-0.6*tc[0]);
        k_f_save[40] = 1e-12 * 6e+19*exp(-1.25*tc[0]);
        k_f_save[41] = 1e-12 * 5.5e+20*exp(-2*tc[0]);
        k_f_save[42] = 1e-12 * 2.2e+22*exp(-2*tc[0]);
        k_f_save[43] = 1e-06 * 3.97e+12*exp(-337.65838275496690812*invT);
        k_f_save[44] = 1e-06 * 2.8e+13*exp(-537.43539907944057177*invT);
        k_f_save[45] = 1e-06 * 1.34e+14*exp(-319.54258278599701271*invT);
        k_f_save[46] = 1e-06 * 1.21e+07*exp(2*tc[0]-2616.7266621845424197*invT);
        k_f_save[47] = 1e-06 * 1e+13*exp(-1811.579996896990906*invT);
        k_f_save[48] = 1e-06 * 1.1e+14;
        k_f_save[49] = 1e-06 * 2.5e+16*exp(-0.8*tc[0]);
        k_f_save[50] = 1e-06 * 3e+13;
        k_f_save[51] = 1e-06 * 1.27e+16*exp(-0.63*tc[0]-192.73198300320765952*invT);
        k_f_save[52] = 1e-06 * 6.6e+08*exp(1.62*tc[0]-5454.8686573231616421*invT);
        k_f_save[53] = 1e-06 * 1.09e+12*exp(0.48*tc[0]+130.83633310922709825*invT);
        k_f_save[54] = 1e-06 * 7.34e+13;
        k_f_save[55] = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1811.579996896990906*invT);
        k_f_save[56] = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1308.3633310922712099*invT);
        k_f_save[57] = 1e-06 * 2.3e+10*exp(1.05*tc[0]-1648.0345805104568626*invT);
        k_f_save[58] = 1e-06 * 1.8e+13;
        k_f_save[59] = 1e-06 * 2e+13;
        k_f_save[60] = 1e-06 * 1.2e+13;
        k_f_save[61] = 1e-06 * 6e+12;
        k_f_save[62] = 1e-06 * 5e+13;
        k_f_save[63] = 1e-06 * 3.4e+06*exp(1.6*tc[0]);
        k_f_save[64] = 1e-06 * 2e+13;
        k_f_save[65] = 1e-06 * 3.2e+13;
        k_f_save[66] = 1e-06 * 1.6e+13;
        k_f_save[67] = 1e-06 * 1.7e+07*exp(2.1*tc[0]-2450.6651624689848177*invT);
        k_f_save[68] = 1e-06 * 4.2e+06*exp(2.1*tc[0]-2450.6651624689848177*invT);
        k_f_save[69] = 1e-06 * 1e+17*exp(-1*tc[0]);
        k_f_save[70] = 1e-06 * 5.6e+12*exp(-1207.7199979313272706*invT);
        k_f_save[71] = 1e-06 * 6.08e+12*exp(0.27*tc[0]-140.90066642532150354*invT);
        k_f_save[72] = 1e-06 * 3e+13;
        k_f_save[73] = 1e-06 * 1.08e+12*exp(0.454*tc[0]-915.85433176458980142*invT);
        k_f_save[74] = 1e-06 * 1.325e+06*exp(2.53*tc[0]-6159.3719894497689893*invT);
        k_f_save[75] = 1e-06 * 5.21e+17*exp(-0.99*tc[0]-795.08233197145705162*invT);
        k_f_save[76] = 1e-06 * 2e+12;
        k_f_save[77] = 1e-06 * 1.15e+08*exp(1.9*tc[0]-3789.2214935095394139*invT);
        k_f_save[78] = 1e-06 * 1e+14;
        k_f_save[79] = 1e-06 * 5e+13*exp(-4025.733326437757114*invT);
        k_f_save[80] = 1e-06 * 1.13e+13*exp(-1725.0267303785790318*invT);
        k_f_save[81] = 1e-06 * 1e+13;
        k_f_save[82] = 1e-06 * 4.3e+07*exp(1.5*tc[0]-40056.046598055690993*invT);
        k_f_save[83] = 1e-06 * 2.16e+08*exp(1.51*tc[0]-1726.0331637101885462*invT);
        k_f_save[84] = 1e-06 * 7.4e+13*exp(-0.37*tc[0]);
        k_f_save[85] = 1e-06 * 35700*exp(2.4*tc[0]+1061.7871648479585929*invT);
        k_f_save[86] = 1e-06 * 2.9e+13*exp(+251.60833290235981963*invT);
        k_f_save[87] = 1e-06 * 1.75e+12*exp(-161.02933305751031412*invT);
        k_f_save[88] = 1e-06 * 5.8e+14*exp(-4810.751325093120613*invT);
        k_f_save[89] = 1e-06 * 5e+13;
        k_f_save[90] = 1e-06 * 3e+13;
        k_f_save[91] = 1e-06 * 2e+13;
        k_f_save[92] = 1e-06 * 1.13e+07*exp(2*tc[0]-1509.6499974141590883*invT);
        k_f_save[93] = 1e-06 * 3e+13;
        k_f_save[94] = 1e-06 * 6.3e+13;
        k_f_save[95] = 1e-06 * 5.6e+07*exp(1.6*tc[0]-2727.4343286615808211*invT);
        k_f_save[96] = 1e-06 * 2.501e+13;
        k_f_save[97] = 1e-06 * 1e+08*exp(1.6*tc[0]-1570.0359973107254064*invT);
        k_f_save[98] = 1e-06 * 4.76e+07*exp(1.228*tc[0]-35.225166606330375885*invT);
        k_f_save[99] = 1e-06 * 5e+13;
        k_f_save[100] = 1e-06 * 3.43e+09*exp(1.18*tc[0]+224.93784961470970529*invT);
        k_f_save[101] = 1e-06 * 5e+12;
        k_f_save[102] = 1e-06 * 5e+12;
        k_f_save[103] = 1e-06 * 1.44e+06*exp(2*tc[0]+422.70199927596451062*invT);
        k_f_save[104] = 1e-06 * 6.3e+06*exp(2*tc[0]-754.82499870707954415*invT);
        k_f_save[105] = 1e-06 * 2e+13;
        k_f_save[106] = 1e-06 * 0.000218*exp(4.5*tc[0]+503.21666580471963925*invT);
        k_f_save[107] = 1e-06 * 504000*exp(2.3*tc[0]-6793.4249883637157836*invT);
        k_f_save[108] = 1e-06 * 3.37e+07*exp(2*tc[0]-7045.0333212660752906*invT);
        k_f_save[109] = 1e-06 * 0.000483*exp(4*tc[0]+1006.4333316094392785*invT);
        k_f_save[110] = 1e-06 * 5e+12;
        k_f_save[111] = 1e-06 * 3.6e+06*exp(2*tc[0]-1258.0416645117993539*invT);
        k_f_save[112] = 1e-06 * 3.54e+06*exp(2.12*tc[0]-437.79849925010609013*invT);
        k_f_save[113] = 1e-06 * 7.5e+12*exp(-1006.4333316094392785*invT);
        k_f_save[114] = 1e-06 * 1.3e+11*exp(+820.24316526169309327*invT);
        k_f_save[115] = 1e-06 * 4.2e+14*exp(-6038.5999896566363532*invT);
        k_f_save[116] = 1e-06 * 2e+13;
        k_f_save[117] = 1e-06 * 1e+12;
        k_f_save[118] = 1e-06 * 2e+13;
        k_f_save[119] = 1e-06 * 1.5e+14*exp(-11875.913312991384373*invT);
        k_f_save[120] = 1e-06 * 1e+12*exp(-4025.733326437757114*invT);
        k_f_save[121] = 1e-06 * 5.8e+13*exp(-289.85279950351855405*invT);
        k_f_save[122] = 1e-06 * 5e+13;
        k_f_save[123] = 1e-06 * 5e+13;
        k_f_save[124] = 1e-06 * 3.3e+13;
        k_f_save[125] = 1e-06 * 1.107e+08*exp(1.79*tc[0]-840.37183189388190385*invT);
        k_f_save[126] = 1e-06 * 5.71e+12*exp(+379.92858268256338761*invT);
        k_f_save[127] = 1e-06 * 4e+13;
        k_f_save[128] = 1e-06 * 3e+13;
        k_f_save[129] = 1e-06 * 6e+13;
        k_f_save[130] = 1e-06 * 5e+13;
        k_f_save[131] = 1e-06 * 3.4e+12*exp(-347.2194994052565562*invT);
        k_f_save[132] = 1e-06 * 9.46e+13*exp(+259.15658288943063781*invT);
        k_f_save[133] = 1e-06 * 5e+13;
        k_f_save[134] = 1e-06 * 1.32e+13*exp(-754.82499870707954415*invT);
        k_f_save[135] = 1e-06 * 500000*exp(2*tc[0]-3638.256493768123164*invT);
        k_f_save[136] = 1e-06 * 3.2e+13;
        k_f_save[137] = 1e-06 * 4e+13;
        k_f_save[138] = 1e-06 * 2.46e+06*exp(2*tc[0]-4161.6018262050320118*invT);
        k_f_save[139] = 1e-06 * 8.1e+11*exp(0.5*tc[0]-2269.5071627792858635*invT);
        k_f_save[140] = 1e-06 * 3e+13;
        k_f_save[141] = 1e-06 * 1.5e+13*exp(-301.92999948283181766*invT);
        k_f_save[142] = 1e-06 * 9e+12*exp(-301.92999948283181766*invT);
        k_f_save[143] = 1e-06 * 2.8e+13;
        k_f_save[144] = 1e-06 * 1.2e+13;
        k_f_save[145] = 1e-06 * 7e+13;
        k_f_save[146] = 1e-06 * 2e+13;
        k_f_save[147] = 1e-06 * 3e+13;
        k_f_save[148] = 1e-06 * 1.2e+13*exp(+286.83349950869023814*invT);
        k_f_save[149] = 1e-06 * 1.6e+13*exp(+286.83349950869023814*invT);
        k_f_save[150] = 1e-06 * 9e+12;
        k_f_save[151] = 1e-06 * 7e+12;
        k_f_save[152] = 1e-06 * 1.4e+13;
        k_f_save[153] = 1e-06 * 4e+13*exp(+276.76916619259583285*invT);
        k_f_save[154] = 1e-06 * 2.675e+13*exp(-14492.639975175927248*invT);
        k_f_save[155] = 1e-06 * 3.6e+10*exp(-4498.7569922941938785*invT);
        k_f_save[156] = 1e-06 * 24500*exp(2.47*tc[0]-2606.6623288684481849*invT);
        k_f_save[157] = 1e-06 * 2.12e+16*exp(-0.97*tc[0]-311.99433279892622295*invT);
        k_f_save[158] = 1e-06 * 4.99e+12*exp(0.1*tc[0]-5334.096657530029006*invT);
        k_f_save[159] = 1e-06 * 2.648e+13;
        k_f_save[160] = 1e-06 * 3320*exp(2.81*tc[0]-2948.8496616156576238*invT);
        k_f_save[161] = 1e-06 * 3e+07*exp(1.5*tc[0]-5001.973658098913802*invT);
        k_f_save[162] = 1e-06 * 1e+07*exp(1.5*tc[0]-5001.973658098913802*invT);
        k_f_save[163] = 1e-06 * 227000*exp(2*tc[0]-4629.5933254034207494*invT);
        k_f_save[164] = 1e-06 * 6.14e+06*exp(1.74*tc[0]-5258.6141576593208811*invT);
        k_f_save[165] = 1e-06 * 2.244e+18*exp(-1*tc[0]-8554.6833186802341515*invT);
        k_f_save[166] = 1e-06 * 1.87e+17*exp(-1*tc[0]-8554.6833186802341515*invT);
        k_f_save[167] = 1e-06 * 7.6e+12*exp(-201.28666632188787844*invT);
        k_f_save[168] = 1e-06 * 1.8e+13*exp(-452.89499922424772649*invT);
        k_f_save[169] = 1e-06 * 4.28e-13*exp(7.6*tc[0]+1776.3548302906606295*invT);
        k_f_save[170] = 1e-06 * 5e+13*exp(-754.82499870707954415*invT);
        k_f_save[171] = 1e-06 * 407000*exp(2.4*tc[0]-100.64333316094393922*invT);
        k_f_save[172] = 1e-06 * 3.98e+12*exp(+120.77199979313272138*invT);
        k_f_save[173] = 1 * 8e+12*exp(0.44*tc[0]-44670.543423484967207*invT);
        k_f_save[174] = 1e-06 * 8.4e+11*exp(-1949.9645799932889076*invT);
        k_f_save[175] = 1e-06 * 1.6e+12*exp(-429.7470325972306*invT);
        k_f_save[176] = 1e-06 * 1e+13;

        Kc_save[0] = 1.0 / (refC) * exp((2 * g_RT[2]) - (g_RT[3]));
        Kc_save[1] = 1.0 / (refC) * exp((g_RT[2] + g_RT[1]) - (g_RT[4]));
        Kc_save[2] = exp((g_RT[2] + g_RT[0]) - (g_RT[1] + g_RT[4]));
        Kc_save[3] = exp((g_RT[2] + g_RT[6]) - (g_RT[4] + g_RT[3]));
        Kc_save[4] = exp((g_RT[2] + g_RT[7]) - (g_RT[4] + g_RT[6]));
        Kc_save[5] = exp((g_RT[2] + g_RT[9]) - (g_RT[1] + g_RT[14]));
        Kc_save[6] = exp((g_RT[2] + g_RT[10]) - (g_RT[1] + g_RT[16]));
        Kc_save[7] = exp((g_RT[2] + g_RT[11]) - (g_RT[0] + g_RT[14]));
        Kc_save[8] = exp((g_RT[2] + g_RT[11]) - (g_RT[1] + g_RT[16]));
        Kc_save[9] = exp((g_RT[2] + g_RT[12]) - (g_RT[1] + g_RT[17]));
        Kc_save[10] = exp((g_RT[2] + g_RT[13]) - (g_RT[4] + g_RT[12]));
        Kc_save[11] = 1.0 / (refC) * exp((g_RT[2] + g_RT[14]) - (g_RT[15]));
        Kc_save[12] = exp((g_RT[2] + g_RT[16]) - (g_RT[4] + g_RT[14]));
        Kc_save[13] = exp((g_RT[2] + g_RT[16]) - (g_RT[1] + g_RT[15]));
        Kc_save[14] = exp((g_RT[2] + g_RT[17]) - (g_RT[4] + g_RT[16]));
        Kc_save[15] = exp((g_RT[2] + g_RT[18]) - (g_RT[4] + g_RT[17]));
        Kc_save[16] = exp((g_RT[2] + g_RT[19]) - (g_RT[4] + g_RT[17]));
        Kc_save[17] = exp((g_RT[2] + g_RT[20]) - (g_RT[4] + g_RT[18]));
        Kc_save[18] = exp((g_RT[2] + g_RT[20]) - (g_RT[4] + g_RT[19]));
        Kc_save[19] = exp((g_RT[2] + g_RT[21]) - (g_RT[9] + g_RT[14]));
        Kc_save[20] = exp((g_RT[2] + g_RT[22]) - (g_RT[1] + g_RT[27]));
        Kc_save[21] = exp((g_RT[2] + g_RT[22]) - (g_RT[4] + g_RT[21]));
        Kc_save[22] = exp((g_RT[2] + g_RT[22]) - (g_RT[14] + g_RT[10]));
        Kc_save[23] = exp((g_RT[2] + g_RT[23]) - (g_RT[1] + g_RT[28]));
        Kc_save[24] = exp((g_RT[2] + g_RT[24]) - (g_RT[12] + g_RT[16]));
        Kc_save[25] = exp((g_RT[2] + g_RT[25]) - (g_RT[12] + g_RT[17]));
        Kc_save[26] = exp((g_RT[2] + g_RT[26]) - (g_RT[4] + g_RT[25]));
        Kc_save[27] = refC * exp((g_RT[2] + g_RT[27]) - (g_RT[1] + 2 * g_RT[14]));
        Kc_save[28] = exp((g_RT[2] + g_RT[28]) - (g_RT[4] + g_RT[27]));
        Kc_save[29] = exp((g_RT[2] + g_RT[28]) - (g_RT[10] + g_RT[15]));
        Kc_save[30] = exp((g_RT[3] + g_RT[14]) - (g_RT[2] + g_RT[15]));
        Kc_save[31] = exp((g_RT[3] + g_RT[17]) - (g_RT[6] + g_RT[16]));
        Kc_save[32] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));
        Kc_save[33] = 1.0 / (refC) * exp((g_RT[1] + 2 * g_RT[3]) - (g_RT[6] + g_RT[3]));
        Kc_save[34] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));
        Kc_save[35] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[30]) - (g_RT[6] + g_RT[30]));
        Kc_save[36] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[31]) - (g_RT[6] + g_RT[31]));
        Kc_save[37] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));
        Kc_save[38] = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));
        Kc_save[39] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));
        Kc_save[40] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));
        Kc_save[41] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[15]) - (g_RT[0] + g_RT[15]));
        Kc_save[42] = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));
        Kc_save[43] = exp((g_RT[1] + g_RT[6]) - (g_RT[2] + g_RT[5]));
        Kc_save[44] = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));
        Kc_save[45] = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));
        Kc_save[46] = exp((g_RT[1] + g_RT[7]) - (g_RT[6] + g_RT[0]));
        Kc_save[47] = exp((g_RT[1] + g_RT[7]) - (g_RT[4] + g_RT[5]));
        Kc_save[48] = exp((g_RT[1] + g_RT[9]) - (g_RT[8] + g_RT[0]));
        Kc_save[49] = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[12]));
        Kc_save[50] = exp((g_RT[1] + g_RT[11]) - (g_RT[9] + g_RT[0]));
        Kc_save[51] = 1.0 / (refC) * exp((g_RT[1] + g_RT[12]) - (g_RT[13]));
        Kc_save[52] = exp((g_RT[1] + g_RT[13]) - (g_RT[12] + g_RT[0]));
        Kc_save[53] = 1.0 / (refC) * exp((g_RT[1] + g_RT[16]) - (g_RT[17]));
        Kc_save[54] = exp((g_RT[1] + g_RT[16]) - (g_RT[0] + g_RT[14]));
        Kc_save[55] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[18]));
        Kc_save[56] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[19]));
        Kc_save[57] = exp((g_RT[1] + g_RT[17]) - (g_RT[16] + g_RT[0]));
        Kc_save[58] = 1.0 / (refC) * exp((g_RT[1] + g_RT[18]) - (g_RT[20]));
        Kc_save[59] = exp((g_RT[1] + g_RT[18]) - (g_RT[0] + g_RT[17]));
        Kc_save[60] = exp((g_RT[1] + g_RT[18]) - (g_RT[4] + g_RT[12]));
        Kc_save[61] = exp((g_RT[1] + g_RT[18]) - (g_RT[11] + g_RT[5]));
        Kc_save[62] = 1.0 / (refC) * exp((g_RT[1] + g_RT[19]) - (g_RT[20]));
        Kc_save[63] = exp((g_RT[1] + g_RT[19]) - (g_RT[1] + g_RT[18]));
        Kc_save[64] = exp((g_RT[1] + g_RT[19]) - (g_RT[0] + g_RT[17]));
        Kc_save[65] = exp((g_RT[1] + g_RT[19]) - (g_RT[4] + g_RT[12]));
        Kc_save[66] = exp((g_RT[1] + g_RT[19]) - (g_RT[11] + g_RT[5]));
        Kc_save[67] = exp((g_RT[1] + g_RT[20]) - (g_RT[18] + g_RT[0]));
        Kc_save[68] = exp((g_RT[1] + g_RT[20]) - (g_RT[19] + g_RT[0]));
        Kc_save[69] = 1.0 / (refC) * exp((g_RT[1] + g_RT[21]) - (g_RT[22]));
        Kc_save[70] = 1.0 / (refC) * exp((g_RT[1] + g_RT[22]) - (g_RT[23]));
        Kc_save[71] = 1.0 / (refC) * exp((g_RT[1] + g_RT[23]) - (g_RT[24]));
        Kc_save[72] = exp((g_RT[1] + g_RT[23]) - (g_RT[0] + g_RT[22]));
        Kc_save[73] = 1.0 / (refC) * exp((g_RT[1] + g_RT[24]) - (g_RT[25]));
        Kc_save[74] = exp((g_RT[1] + g_RT[24]) - (g_RT[23] + g_RT[0]));
        Kc_save[75] = 1.0 / (refC) * exp((g_RT[1] + g_RT[25]) - (g_RT[26]));
        Kc_save[76] = exp((g_RT[1] + g_RT[25]) - (g_RT[0] + g_RT[24]));
        Kc_save[77] = exp((g_RT[1] + g_RT[26]) - (g_RT[25] + g_RT[0]));
        Kc_save[78] = exp((g_RT[1] + g_RT[27]) - (g_RT[11] + g_RT[14]));
        Kc_save[79] = exp((g_RT[1] + g_RT[28]) - (g_RT[27] + g_RT[0]));
        Kc_save[80] = exp((g_RT[1] + g_RT[28]) - (g_RT[12] + g_RT[14]));
        Kc_save[81] = exp((g_RT[1] + g_RT[29]) - (g_RT[1] + g_RT[28]));
        Kc_save[82] = 1.0 / (refC) * exp((g_RT[0] + g_RT[14]) - (g_RT[17]));
        Kc_save[83] = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));
        Kc_save[84] = 1.0 / (refC) * exp((2 * g_RT[4]) - (g_RT[7]));
        Kc_save[85] = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));
        Kc_save[86] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
        Kc_save[87] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
        Kc_save[88] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
        Kc_save[89] = exp((g_RT[4] + g_RT[8]) - (g_RT[1] + g_RT[14]));
        Kc_save[90] = exp((g_RT[4] + g_RT[9]) - (g_RT[1] + g_RT[16]));
        Kc_save[91] = exp((g_RT[4] + g_RT[10]) - (g_RT[1] + g_RT[17]));
        Kc_save[92] = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[5]));
        Kc_save[93] = exp((g_RT[4] + g_RT[11]) - (g_RT[1] + g_RT[17]));
        Kc_save[94] = 1.0 / (refC) * exp((g_RT[4] + g_RT[12]) - (g_RT[20]));
        Kc_save[95] = exp((g_RT[4] + g_RT[12]) - (g_RT[10] + g_RT[5]));
        Kc_save[96] = exp((g_RT[4] + g_RT[12]) - (g_RT[11] + g_RT[5]));
        Kc_save[97] = exp((g_RT[4] + g_RT[13]) - (g_RT[12] + g_RT[5]));
        Kc_save[98] = exp((g_RT[4] + g_RT[14]) - (g_RT[1] + g_RT[15]));
        Kc_save[99] = exp((g_RT[4] + g_RT[16]) - (g_RT[5] + g_RT[14]));
        Kc_save[100] = exp((g_RT[4] + g_RT[17]) - (g_RT[16] + g_RT[5]));
        Kc_save[101] = exp((g_RT[4] + g_RT[18]) - (g_RT[5] + g_RT[17]));
        Kc_save[102] = exp((g_RT[4] + g_RT[19]) - (g_RT[5] + g_RT[17]));
        Kc_save[103] = exp((g_RT[4] + g_RT[20]) - (g_RT[18] + g_RT[5]));
        Kc_save[104] = exp((g_RT[4] + g_RT[20]) - (g_RT[19] + g_RT[5]));
        Kc_save[105] = exp((g_RT[4] + g_RT[21]) - (g_RT[1] + g_RT[27]));
        Kc_save[106] = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[28]));
        Kc_save[107] = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[29]));
        Kc_save[108] = exp((g_RT[4] + g_RT[22]) - (g_RT[21] + g_RT[5]));
        Kc_save[109] = exp((g_RT[4] + g_RT[22]) - (g_RT[12] + g_RT[14]));
        Kc_save[110] = exp((g_RT[4] + g_RT[23]) - (g_RT[5] + g_RT[22]));
        Kc_save[111] = exp((g_RT[4] + g_RT[24]) - (g_RT[23] + g_RT[5]));
        Kc_save[112] = exp((g_RT[4] + g_RT[26]) - (g_RT[25] + g_RT[5]));
        Kc_save[113] = exp((g_RT[4] + g_RT[28]) - (g_RT[27] + g_RT[5]));
        Kc_save[114] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
        Kc_save[115] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
        Kc_save[116] = exp((g_RT[6] + g_RT[10]) - (g_RT[4] + g_RT[17]));
        Kc_save[117] = exp((g_RT[6] + g_RT[12]) - (g_RT[3] + g_RT[13]));
        Kc_save[118] = exp((g_RT[6] + g_RT[12]) - (g_RT[4] + g_RT[19]));
        Kc_save[119] = exp((g_RT[6] + g_RT[14]) - (g_RT[4] + g_RT[15]));
        Kc_save[120] = exp((g_RT[6] + g_RT[17]) - (g_RT[16] + g_RT[7]));
        Kc_save[121] = exp((g_RT[8] + g_RT[3]) - (g_RT[2] + g_RT[14]));
        Kc_save[122] = exp((g_RT[8] + g_RT[10]) - (g_RT[1] + g_RT[21]));
        Kc_save[123] = exp((g_RT[8] + g_RT[12]) - (g_RT[1] + g_RT[22]));
        Kc_save[124] = exp((g_RT[9] + g_RT[3]) - (g_RT[2] + g_RT[16]));
        Kc_save[125] = exp((g_RT[9] + g_RT[0]) - (g_RT[1] + g_RT[10]));
        Kc_save[126] = exp((g_RT[9] + g_RT[5]) - (g_RT[1] + g_RT[17]));
        Kc_save[127] = exp((g_RT[9] + g_RT[10]) - (g_RT[1] + g_RT[22]));
        Kc_save[128] = exp((g_RT[9] + g_RT[12]) - (g_RT[1] + g_RT[23]));
        Kc_save[129] = exp((g_RT[9] + g_RT[13]) - (g_RT[1] + g_RT[24]));
        Kc_save[130] = 1.0 / (refC) * exp((g_RT[9] + g_RT[14]) - (g_RT[27]));
        Kc_save[131] = exp((g_RT[9] + g_RT[15]) - (g_RT[16] + g_RT[14]));
        Kc_save[132] = exp((g_RT[9] + g_RT[17]) - (g_RT[1] + g_RT[28]));
        Kc_save[133] = exp((g_RT[9] + g_RT[27]) - (g_RT[14] + g_RT[22]));
        Kc_save[134] = exp((g_RT[10] + g_RT[3]) - (g_RT[4] + g_RT[16]));
        Kc_save[135] = exp((g_RT[10] + g_RT[0]) - (g_RT[1] + g_RT[12]));
        Kc_save[136] = exp((2 * g_RT[10]) - (g_RT[0] + g_RT[22]));
        Kc_save[137] = exp((g_RT[10] + g_RT[12]) - (g_RT[1] + g_RT[24]));
        Kc_save[138] = exp((g_RT[10] + g_RT[13]) - (2 * g_RT[12]));
        Kc_save[139] = 1.0 / (refC) * exp((g_RT[10] + g_RT[14]) - (g_RT[28]));
        Kc_save[140] = exp((g_RT[10] + g_RT[27]) - (g_RT[23] + g_RT[14]));
        Kc_save[141] = exp((g_RT[11] + g_RT[30]) - (g_RT[10] + g_RT[30]));
        Kc_save[142] = exp((g_RT[11] + g_RT[31]) - (g_RT[10] + g_RT[31]));
        Kc_save[143] = refC * exp((g_RT[11] + g_RT[3]) - (g_RT[1] + g_RT[4] + g_RT[14]));
        Kc_save[144] = exp((g_RT[11] + g_RT[3]) - (g_RT[14] + g_RT[5]));
        Kc_save[145] = exp((g_RT[11] + g_RT[0]) - (g_RT[12] + g_RT[1]));
        Kc_save[146] = 1.0 / (refC) * exp((g_RT[11] + g_RT[5]) - (g_RT[20]));
        Kc_save[147] = exp((g_RT[11] + g_RT[5]) - (g_RT[10] + g_RT[5]));
        Kc_save[148] = exp((g_RT[11] + g_RT[12]) - (g_RT[1] + g_RT[24]));
        Kc_save[149] = exp((g_RT[11] + g_RT[13]) - (2 * g_RT[12]));
        Kc_save[150] = exp((g_RT[11] + g_RT[14]) - (g_RT[10] + g_RT[14]));
        Kc_save[151] = exp((g_RT[11] + g_RT[15]) - (g_RT[10] + g_RT[15]));
        Kc_save[152] = exp((g_RT[11] + g_RT[15]) - (g_RT[14] + g_RT[17]));
        Kc_save[153] = exp((g_RT[11] + g_RT[26]) - (g_RT[12] + g_RT[25]));
        Kc_save[154] = exp((g_RT[12] + g_RT[3]) - (g_RT[2] + g_RT[19]));
        Kc_save[155] = exp((g_RT[12] + g_RT[3]) - (g_RT[4] + g_RT[17]));
        Kc_save[156] = exp((g_RT[12] + g_RT[7]) - (g_RT[6] + g_RT[13]));
        Kc_save[157] = 1.0 / (refC) * exp((2 * g_RT[12]) - (g_RT[26]));
        Kc_save[158] = exp((2 * g_RT[12]) - (g_RT[1] + g_RT[25]));
        Kc_save[159] = exp((g_RT[12] + g_RT[16]) - (g_RT[13] + g_RT[14]));
        Kc_save[160] = exp((g_RT[12] + g_RT[17]) - (g_RT[16] + g_RT[13]));
        Kc_save[161] = exp((g_RT[12] + g_RT[20]) - (g_RT[18] + g_RT[13]));
        Kc_save[162] = exp((g_RT[12] + g_RT[20]) - (g_RT[19] + g_RT[13]));
        Kc_save[163] = exp((g_RT[12] + g_RT[24]) - (g_RT[23] + g_RT[13]));
        Kc_save[164] = exp((g_RT[12] + g_RT[26]) - (g_RT[25] + g_RT[13]));
        Kc_save[165] = refC * exp((g_RT[16] + g_RT[5]) - (g_RT[1] + g_RT[14] + g_RT[5]));
        Kc_save[166] = refC * exp((g_RT[16]) - (g_RT[1] + g_RT[14]));
        Kc_save[167] = exp((g_RT[16] + g_RT[3]) - (g_RT[6] + g_RT[14]));
        Kc_save[168] = exp((g_RT[18] + g_RT[3]) - (g_RT[6] + g_RT[17]));
        Kc_save[169] = exp((g_RT[19] + g_RT[3]) - (g_RT[6] + g_RT[17]));
        Kc_save[170] = exp((g_RT[21] + g_RT[3]) - (g_RT[16] + g_RT[14]));
        Kc_save[171] = exp((g_RT[21] + g_RT[0]) - (g_RT[1] + g_RT[22]));
        Kc_save[172] = exp((g_RT[23] + g_RT[3]) - (g_RT[16] + g_RT[17]));
        Kc_save[173] = refC * exp((g_RT[24]) - (g_RT[0] + g_RT[22]));
        Kc_save[174] = exp((g_RT[25] + g_RT[3]) - (g_RT[6] + g_RT[24]));
        Kc_save[175] = refC * exp((g_RT[27] + g_RT[3]) - (g_RT[4] + 2 * g_RT[14]));
        Kc_save[176] = refC * exp((2 * g_RT[27]) - (2 * g_RT[14] + g_RT[22]));
    }

    /*reaction 1: 2 O + M <=> O2 + M */
    phi_f = sc[2]*sc[2];
    alpha = mixture + 1.4*sc[0] + 14.4*sc[5] + sc[13] + 0.75*sc[14] + 2.6*sc[15] + 2*sc[26] + -0.17*sc[31];
    k_f = alpha * k_f_save[0];
    q_f[0] = phi_f * k_f;
    phi_r = sc[3];
    Kc = Kc_save[0];
    k_r = k_f / Kc;
    q_r[0] = phi_r * k_r;

    /*reaction 2: O + H + M <=> OH + M */
    phi_f = sc[2]*sc[1];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = alpha * k_f_save[1];
    q_f[1] = phi_f * k_f;
    phi_r = sc[4];
    Kc = Kc_save[1];
    k_r = k_f / Kc;
    q_r[1] = phi_r * k_r;

    /*reaction 3: O + H2 <=> H + OH */
    phi_f = sc[2]*sc[0];
    k_f = k_f_save[2];
    q_f[2] = phi_f * k_f;
    phi_r = sc[1]*sc[4];
    Kc = Kc_save[2];
    k_r = k_f / Kc;
    q_r[2] = phi_r * k_r;

    /*reaction 4: O + HO2 <=> OH + O2 */
    phi_f = sc[2]*sc[6];
    k_f = k_f_save[3];
    q_f[3] = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = Kc_save[3];
    k_r = k_f / Kc;
    q_r[3] = phi_r * k_r;

    /*reaction 5: O + H2O2 <=> OH + HO2 */
    phi_f = sc[2]*sc[7];
    k_f = k_f_save[4];
    q_f[4] = phi_f * k_f;
    phi_r = sc[4]*sc[6];
    Kc = Kc_save[4];
    k_r = k_f / Kc;
    q_r[4] = phi_r * k_r;

    /*reaction 6: O + CH <=> H + CO */
    phi_f = sc[2]*sc[9];
    k_f = k_f_save[5];
    q_f[5] = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[5];
    k_r = k_f / Kc;
    q_r[5] = phi_r * k_r;

    /*reaction 7: O + CH2 <=> H + HCO */
    phi_f = sc[2]*sc[10];
    k_f = k_f_save[6];
    q_f[6] = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_save[6];
    k_r = k_f / Kc;
    q_r[6] = phi_r * k_r;

    /*reaction 8: O + CH2(S) <=> H2 + CO */
    phi_f = sc[2]*sc[11];
    k_f = k_f_save[7];
    q_f[7] = phi_f * k_f;
    phi_r = sc[0]*sc[14];
    Kc = Kc_save[7];
    k_r = k_f / Kc;
    q_r[7] = phi_r * k_r;

    /*reaction 9: O + CH2(S) <=> H + HCO */
    phi_f = sc[2]*sc[11];
    k_f = k_f_save[8];
    q_f[8] = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_save[8];
    k_r = k_f / Kc;
    q_r[8] = phi_r * k_r;

    /*reaction 10: O + CH3 <=> H + CH2O */
    phi_f = sc[2]*sc[12];
    k_f = k_f_save[9];
    q_f[9] = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[9];
    k_r = k_f / Kc;
    q_r[9] = phi_r * k_r;

    /*reaction 11: O + CH4 <=> OH + CH3 */
    phi_f = sc[2]*sc[13];
    k_f = k_f_save[10];
    q_f[10] = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = Kc_save[10];
    k_r = k_f / Kc;
    q_r[10] = phi_r * k_r;

    /*reaction 12: O + CO + M <=> CO2 + M */
    phi_f = sc[2]*sc[14];
    alpha = mixture + sc[0] + 5*sc[3] + 5*sc[5] + sc[13] + 0.5*sc[14] + 2.5*sc[15] + 2*sc[26] + -0.5*sc[31];
    k_f = alpha * k_f_save[11];
    q_f[11] = phi_f * k_f;
    phi_r = sc[15];
    Kc = Kc_save[11];
    k_r = k_f / Kc;
    q_r[11] = phi_r * k_r;

    /*reaction 13: O + HCO <=> OH + CO */
    phi_f = sc[2]*sc[16];
    k_f = k_f_save[12];
    q_f[12] = phi_f * k_f;
    phi_r = sc[4]*sc[14];
    Kc = Kc_save[12];
    k_r = k_f / Kc;
    q_r[12] = phi_r * k_r;

    /*reaction 14: O + HCO <=> H + CO2 */
    phi_f = sc[2]*sc[16];
    k_f = k_f_save[13];
    q_f[13] = phi_f * k_f;
    phi_r = sc[1]*sc[15];
    Kc = Kc_save[13];
    k_r = k_f / Kc;
    q_r[13] = phi_r * k_r;

    /*reaction 15: O + CH2O <=> OH + HCO */
    phi_f = sc[2]*sc[17];
    k_f = k_f_save[14];
    q_f[14] = phi_f * k_f;
    phi_r = sc[4]*sc[16];
    Kc = Kc_save[14];
    k_r = k_f / Kc;
    q_r[14] = phi_r * k_r;

    /*reaction 16: O + CH2OH <=> OH + CH2O */
    phi_f = sc[2]*sc[18];
    k_f = k_f_save[15];
    q_f[15] = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[15];
    k_r = k_f / Kc;
    q_r[15] = phi_r * k_r;

    /*reaction 17: O + CH3O <=> OH + CH2O */
    phi_f = sc[2]*sc[19];
    k_f = k_f_save[16];
    q_f[16] = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[16];
    k_r = k_f / Kc;
    q_r[16] = phi_r * k_r;

    /*reaction 18: O + CH3OH <=> OH + CH2OH */
    phi_f = sc[2]*sc[20];
    k_f = k_f_save[17];
    q_f[17] = phi_f * k_f;
    phi_r = sc[4]*sc[18];
    Kc = Kc_save[17];
    k_r = k_f / Kc;
    q_r[17] = phi_r * k_r;

    /*reaction 19: O + CH3OH <=> OH + CH3O */
    phi_f = sc[2]*sc[20];
    k_f = k_f_save[18];
    q_f[18] = phi_f * k_f;
    phi_r = sc[4]*sc[19];
    Kc = Kc_save[18];
    k_r = k_f / Kc;
    q_r[18] = phi_r * k_r;

    /*reaction 20: O + C2H <=> CH + CO */
    phi_f = sc[2]*sc[21];
    k_f = k_f_save[19];
    q_f[19] = phi_f * k_f;
    phi_r = sc[9]*sc[14];
    Kc = Kc_save[19];
    k_r = k_f / Kc;
    q_r[19] = phi_r * k_r;

    /*reaction 21: O + C2H2 <=> H + HCCO */
    phi_f = sc[2]*sc[22];
    k_f = k_f_save[20];
    q_f[20] = phi_f * k_f;
    phi_r = sc[1]*sc[27];
    Kc = Kc_save[20];
    k_r = k_f / Kc;
    q_r[20] = phi_r * k_r;

    /*reaction 22: O + C2H2 <=> OH + C2H */
    phi_f = sc[2]*sc[22];
    k_f = k_f_save[21];
    q_f[21] = phi_f * k_f;
    phi_r = sc[4]*sc[21];
    Kc = Kc_save[21];
    k_r = k_f / Kc;
    q_r[21] = phi_r * k_r;

    /*reaction 23: O + C2H2 <=> CO + CH2 */
    phi_f = sc[2]*sc[22];
    k_f = k_f_save[22];
    q_f[22] = phi_f * k_f;
    phi_r = sc[14]*sc[10];
    Kc = Kc_save[22];
    k_r = k_f / Kc;
    q_r[22] = phi_r * k_r;

    /*reaction 24: O + C2H3 <=> H + CH2CO */
    phi_f = sc[2]*sc[23];
    k_f = k_f_save[23];
    q_f[23] = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[23];
    k_r = k_f / Kc;
    q_r[23] = phi_r * k_r;

    /*reaction 25: O + C2H4 <=> CH3 + HCO */
    phi_f = sc[2]*sc[24];
    k_f = k_f_save[24];
    q_f[24] = phi_f * k_f;
    phi_r = sc[12]*sc[16];
    Kc = Kc_save[24];
    k_r = k_f / Kc;
    q_r[24] = phi_r * k_r;

    /*reaction 26: O + C2H5 <=> CH3 + CH2O */
    phi_f = sc[2]*sc[25];
    k_f = k_f_save[25];
    q_f[25] = phi_f * k_f;
    phi_r = sc[12]*sc[17];
    Kc = Kc_save[25];
    k_r = k_f / Kc;
    q_r[25] = phi_r * k_r;

    /*reaction 27: O + C2H6 <=> OH + C2H5 */
    phi_f = sc[2]*sc[26];
    k_f = k_f_save[26];
    q_f[26] = phi_f * k_f;
    phi_r = sc[4]*sc[25];
    Kc = Kc_save[26];
    k_r = k_f / Kc;
    q_r[26] = phi_r * k_r;

    /*reaction 28: O + HCCO <=> H + 2 CO */
    phi_f = sc[2]*sc[27];
    k_f = k_f_save[27];
    q_f[27] = phi_f * k_f;
    phi_r = sc[1]*sc[14]*sc[14];
    Kc = Kc_save[27];
    k_r = k_f / Kc;
    q_r[27] = phi_r * k_r;

    /*reaction 29: O + CH2CO <=> OH + HCCO */
    phi_f = sc[2]*sc[28];
    k_f = k_f_save[28];
    q_f[28] = phi_f * k_f;
    phi_r = sc[4]*sc[27];
    Kc = Kc_save[28];
    k_r = k_f / Kc;
    q_r[28] = phi_r * k_r;

    /*reaction 30: O + CH2CO <=> CH2 + CO2 */
    phi_f = sc[2]*sc[28];
    k_f = k_f_save[29];
    q_f[29] = phi_f * k_f;
    phi_r = sc[10]*sc[15];
    Kc = Kc_save[29];
    k_r = k_f / Kc;
    q_r[29] = phi_r * k_r;

    /*reaction 31: O2 + CO <=> O + CO2 */
    phi_f = sc[3]*sc[14];
    k_f = k_f_save[30];
    q_f[30] = phi_f * k_f;
    phi_r = sc[2]*sc[15];
    Kc = Kc_save[30];
    k_r = k_f / Kc;
    q_r[30] = phi_r * k_r;

    /*reaction 32: O2 + CH2O <=> HO2 + HCO */
    phi_f = sc[3]*sc[17];
    k_f = k_f_save[31];
    q_f[31] = phi_f * k_f;
    phi_r = sc[6]*sc[16];
    Kc = Kc_save[31];
    k_r = k_f / Kc;
    q_r[31] = phi_r * k_r;

    /*reaction 33: H + O2 + M <=> HO2 + M */
    phi_f = sc[1]*sc[3];
    alpha = mixture + -1*sc[3] + -1*sc[5] + -0.25*sc[14] + 0.5*sc[15] + 0.5*sc[26] + -1*sc[30] + -1*sc[31];
    k_f = alpha * k_f_save[32];
    q_f[32] = phi_f * k_f;
    phi_r = sc[6];
    Kc = Kc_save[32];
    k_r = k_f / Kc;
    q_r[32] = phi_r * k_r;

    /*reaction 34: H + 2 O2 <=> HO2 + O2 */
    phi_f = sc[1]*sc[3]*sc[3];
    k_f = k_f_save[33];
    q_f[33] = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = Kc_save[33];
    k_r = k_f / Kc;
    q_r[33] = phi_r * k_r;

    /*reaction 35: H + O2 + H2O <=> HO2 + H2O */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = k_f_save[34];
    q_f[34] = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = Kc_save[34];
    k_r = k_f / Kc;
    q_r[34] = phi_r * k_r;

    /*reaction 36: H + O2 + N2 <=> HO2 + N2 */
    phi_f = sc[1]*sc[3]*sc[30];
    k_f = k_f_save[35];
    q_f[35] = phi_f * k_f;
    phi_r = sc[6]*sc[30];
    Kc = Kc_save[35];
    k_r = k_f / Kc;
    q_r[35] = phi_r * k_r;

    /*reaction 37: H + O2 + AR <=> HO2 + AR */
    phi_f = sc[1]*sc[3]*sc[31];
    k_f = k_f_save[36];
    q_f[36] = phi_f * k_f;
    phi_r = sc[6]*sc[31];
    Kc = Kc_save[36];
    k_r = k_f / Kc;
    q_r[36] = phi_r * k_r;

    /*reaction 38: H + O2 <=> O + OH */
    phi_f = sc[1]*sc[3];
    k_f = k_f_save[37];
    q_f[37] = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = Kc_save[37];
    k_r = k_f / Kc;
    q_r[37] = phi_r * k_r;

    /*reaction 39: 2 H + M <=> H2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[0] + -1*sc[5] + sc[13] + -1*sc[15] + 2*sc[26] + -0.37*sc[31];
    k_f = alpha * k_f_save[38];
    q_f[38] = phi_f * k_f;
    phi_r = sc[0];
    Kc = Kc_save[38];
    k_r = k_f / Kc;
    q_r[38] = phi_r * k_r;

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    phi_f = sc[1]*sc[1]*sc[0];
    k_f = k_f_save[39];
    q_f[39] = phi_f * k_f;
    phi_r = sc[0]*sc[0];
    Kc = Kc_save[39];
    k_r = k_f / Kc;
    q_r[39] = phi_r * k_r;

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    phi_f = sc[1]*sc[1]*sc[5];
    k_f = k_f_save[40];
    q_f[40] = phi_f * k_f;
    phi_r = sc[0]*sc[5];
    Kc = Kc_save[40];
    k_r = k_f / Kc;
    q_r[40] = phi_r * k_r;

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    phi_f = sc[1]*sc[1]*sc[15];
    k_f = k_f_save[41];
    q_f[41] = phi_f * k_f;
    phi_r = sc[0]*sc[15];
    Kc = Kc_save[41];
    k_r = k_f / Kc;
    q_r[41] = phi_r * k_r;

    /*reaction 43: H + OH + M <=> H2O + M */
    phi_f = sc[1]*sc[4];
    alpha = mixture + -0.27*sc[0] + 2.65*sc[5] + sc[13] + 2*sc[26] + -0.62*sc[31];
    k_f = alpha * k_f_save[42];
    q_f[42] = phi_f * k_f;
    phi_r = sc[5];
    Kc = Kc_save[42];
    k_r = k_f / Kc;
    q_r[42] = phi_r * k_r;

    /*reaction 44: H + HO2 <=> O + H2O */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[43];
    q_f[43] = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = Kc_save[43];
    k_r = k_f / Kc;
    q_r[43] = phi_r * k_r;

    /*reaction 45: H + HO2 <=> O2 + H2 */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[44];
    q_f[44] = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = Kc_save[44];
    k_r = k_f / Kc;
    q_r[44] = phi_r * k_r;

    /*reaction 46: H + HO2 <=> 2 OH */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[45];
    q_f[45] = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = Kc_save[45];
    k_r = k_f / Kc;
    q_r[45] = phi_r * k_r;

    /*reaction 47: H + H2O2 <=> HO2 + H2 */
    phi_f = sc[1]*sc[7];
    k_f = k_f_save[46];
    q_f[46] = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = Kc_save[46];
    k_r = k_f / Kc;
    q_r[46] = phi_r * k_r;

    /*reaction 48: H + H2O2 <=> OH + H2O */
    phi_f = sc[1]*sc[7];
    k_f = k_f_save[47];
    q_f[47] = phi_f * k_f;
    phi_r = sc[4]*sc[5];
    Kc = Kc_save[47];
    k_r = k_f / Kc;
    q_r[47] = phi_r * k_r;

    /*reaction 49: H + CH <=> C + H2 */
    phi_f = sc[1]*sc[9];
    k_f = k_f_save[48];
    q_f[48] = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = Kc_save[48];
    k_r = k_f / Kc;
    q_r[48] = phi_r * k_r;

    /*reaction 50: H + CH2 (+M) <=> CH3 (+M) */
    phi_f = sc[1]*sc[10];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[49];
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
    q_f[49] = phi_f * k_f;
    phi_r = sc[12];
    Kc = Kc_save[49];
    k_r = k_f / Kc;
    q_r[49] = phi_r * k_r;

    /*reaction 51: H + CH2(S) <=> CH + H2 */
    phi_f = sc[1]*sc[11];
    k_f = k_f_save[50];
    q_f[50] = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = Kc_save[50];
    k_r = k_f / Kc;
    q_r[50] = phi_r * k_r;

    /*reaction 52: H + CH3 (+M) <=> CH4 (+M) */
    phi_f = sc[1]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[51];
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
    q_f[51] = phi_f * k_f;
    phi_r = sc[13];
    Kc = Kc_save[51];
    k_r = k_f / Kc;
    q_r[51] = phi_r * k_r;

    /*reaction 53: H + CH4 <=> CH3 + H2 */
    phi_f = sc[1]*sc[13];
    k_f = k_f_save[52];
    q_f[52] = phi_f * k_f;
    phi_r = sc[12]*sc[0];
    Kc = Kc_save[52];
    k_r = k_f / Kc;
    q_r[52] = phi_r * k_r;

    /*reaction 54: H + HCO (+M) <=> CH2O (+M) */
    phi_f = sc[1]*sc[16];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[53];
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
    q_f[53] = phi_f * k_f;
    phi_r = sc[17];
    Kc = Kc_save[53];
    k_r = k_f / Kc;
    q_r[53] = phi_r * k_r;

    /*reaction 55: H + HCO <=> H2 + CO */
    phi_f = sc[1]*sc[16];
    k_f = k_f_save[54];
    q_f[54] = phi_f * k_f;
    phi_r = sc[0]*sc[14];
    Kc = Kc_save[54];
    k_r = k_f / Kc;
    q_r[54] = phi_r * k_r;

    /*reaction 56: H + CH2O (+M) <=> CH2OH (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[55];
    redP = 1e-12 * alpha / k_f * 1.27e+32*exp(-4.82*tc[0]-3286.0048277048194905*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2813*exp(T/-103))+ (0.7187*exp(T/-1291))+ (exp(-4160/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[55] = phi_f * k_f;
    phi_r = sc[18];
    Kc = Kc_save[55];
    k_r = k_f / Kc;
    q_r[55] = phi_r * k_r;

    /*reaction 57: H + CH2O (+M) <=> CH3O (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[56];
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
    q_f[56] = phi_f * k_f;
    phi_r = sc[19];
    Kc = Kc_save[56];
    k_r = k_f / Kc;
    q_r[56] = phi_r * k_r;

    /*reaction 58: H + CH2O <=> HCO + H2 */
    phi_f = sc[1]*sc[17];
    k_f = k_f_save[57];
    q_f[57] = phi_f * k_f;
    phi_r = sc[16]*sc[0];
    Kc = Kc_save[57];
    k_r = k_f / Kc;
    q_r[57] = phi_r * k_r;

    /*reaction 59: H + CH2OH (+M) <=> CH3OH (+M) */
    phi_f = sc[1]*sc[18];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[58];
    redP = 1e-12 * alpha / k_f * 3e+31*exp(-4.8*tc[0]-1660.6149971555748834*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2321*exp(T/-338))+ (0.7679*exp(T/-1812))+ (exp(-5081/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[58] = phi_f * k_f;
    phi_r = sc[20];
    Kc = Kc_save[58];
    k_r = k_f / Kc;
    q_r[58] = phi_r * k_r;

    /*reaction 60: H + CH2OH <=> H2 + CH2O */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[59];
    q_f[59] = phi_f * k_f;
    phi_r = sc[0]*sc[17];
    Kc = Kc_save[59];
    k_r = k_f / Kc;
    q_r[59] = phi_r * k_r;

    /*reaction 61: H + CH2OH <=> OH + CH3 */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[60];
    q_f[60] = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = Kc_save[60];
    k_r = k_f / Kc;
    q_r[60] = phi_r * k_r;

    /*reaction 62: H + CH2OH <=> CH2(S) + H2O */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[61];
    q_f[61] = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[61];
    k_r = k_f / Kc;
    q_r[61] = phi_r * k_r;

    /*reaction 63: H + CH3O (+M) <=> CH3OH (+M) */
    phi_f = sc[1]*sc[19];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[62];
    redP = 1e-12 * alpha / k_f * 8.6e+28*exp(-4*tc[0]-1522.2304140592771091*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.1098*exp(T/-144))+ (0.8902*exp(T/-2838))+ (exp(-45569/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[62] = phi_f * k_f;
    phi_r = sc[20];
    Kc = Kc_save[62];
    k_r = k_f / Kc;
    q_r[62] = phi_r * k_r;

    /*reaction 64: H + CH3O <=> H + CH2OH */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[63];
    q_f[63] = phi_f * k_f;
    phi_r = sc[1]*sc[18];
    Kc = Kc_save[63];
    k_r = k_f / Kc;
    q_r[63] = phi_r * k_r;

    /*reaction 65: H + CH3O <=> H2 + CH2O */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[64];
    q_f[64] = phi_f * k_f;
    phi_r = sc[0]*sc[17];
    Kc = Kc_save[64];
    k_r = k_f / Kc;
    q_r[64] = phi_r * k_r;

    /*reaction 66: H + CH3O <=> OH + CH3 */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[65];
    q_f[65] = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = Kc_save[65];
    k_r = k_f / Kc;
    q_r[65] = phi_r * k_r;

    /*reaction 67: H + CH3O <=> CH2(S) + H2O */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[66];
    q_f[66] = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[66];
    k_r = k_f / Kc;
    q_r[66] = phi_r * k_r;

    /*reaction 68: H + CH3OH <=> CH2OH + H2 */
    phi_f = sc[1]*sc[20];
    k_f = k_f_save[67];
    q_f[67] = phi_f * k_f;
    phi_r = sc[18]*sc[0];
    Kc = Kc_save[67];
    k_r = k_f / Kc;
    q_r[67] = phi_r * k_r;

    /*reaction 69: H + CH3OH <=> CH3O + H2 */
    phi_f = sc[1]*sc[20];
    k_f = k_f_save[68];
    q_f[68] = phi_f * k_f;
    phi_r = sc[19]*sc[0];
    Kc = Kc_save[68];
    k_r = k_f / Kc;
    q_r[68] = phi_r * k_r;

    /*reaction 70: H + C2H (+M) <=> C2H2 (+M) */
    phi_f = sc[1]*sc[21];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[69];
    redP = 1e-12 * alpha / k_f * 3.75e+33*exp(-4.8*tc[0]-956.11166502896742259*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.3536*exp(T/-132))+ (0.6464*exp(T/-1315))+ (exp(-5566/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[69] = phi_f * k_f;
    phi_r = sc[22];
    Kc = Kc_save[69];
    k_r = k_f / Kc;
    q_r[69] = phi_r * k_r;

    /*reaction 71: H + C2H2 (+M) <=> C2H3 (+M) */
    phi_f = sc[1]*sc[22];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[70];
    redP = 1e-12 * alpha / k_f * 3.8e+40*exp(-7.27*tc[0]-3633.2243271100760467*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2493*exp(T/-98.5))+ (0.7507*exp(T/-1302))+ (exp(-4167/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[70] = phi_f * k_f;
    phi_r = sc[23];
    Kc = Kc_save[70];
    k_r = k_f / Kc;
    q_r[70] = phi_r * k_r;

    /*reaction 72: H + C2H3 (+M) <=> C2H4 (+M) */
    phi_f = sc[1]*sc[23];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[71];
    redP = 1e-12 * alpha / k_f * 1.4e+30*exp(-3.86*tc[0]-1670.6793304716693456*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.218*exp(T/-207.5))+ (0.782*exp(T/-2663))+ (exp(-6095/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[71] = phi_f * k_f;
    phi_r = sc[24];
    Kc = Kc_save[71];
    k_r = k_f / Kc;
    q_r[71] = phi_r * k_r;

    /*reaction 73: H + C2H3 <=> H2 + C2H2 */
    phi_f = sc[1]*sc[23];
    k_f = k_f_save[72];
    q_f[72] = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = Kc_save[72];
    k_r = k_f / Kc;
    q_r[72] = phi_r * k_r;

    /*reaction 74: H + C2H4 (+M) <=> C2H5 (+M) */
    phi_f = sc[1]*sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[73];
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
    q_f[73] = phi_f * k_f;
    phi_r = sc[25];
    Kc = Kc_save[73];
    k_r = k_f / Kc;
    q_r[73] = phi_r * k_r;

    /*reaction 75: H + C2H4 <=> C2H3 + H2 */
    phi_f = sc[1]*sc[24];
    k_f = k_f_save[74];
    q_f[74] = phi_f * k_f;
    phi_r = sc[23]*sc[0];
    Kc = Kc_save[74];
    k_r = k_f / Kc;
    q_r[74] = phi_r * k_r;

    /*reaction 76: H + C2H5 (+M) <=> C2H6 (+M) */
    phi_f = sc[1]*sc[25];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[75];
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
    q_f[75] = phi_f * k_f;
    phi_r = sc[26];
    Kc = Kc_save[75];
    k_r = k_f / Kc;
    q_r[75] = phi_r * k_r;

    /*reaction 77: H + C2H5 <=> H2 + C2H4 */
    phi_f = sc[1]*sc[25];
    k_f = k_f_save[76];
    q_f[76] = phi_f * k_f;
    phi_r = sc[0]*sc[24];
    Kc = Kc_save[76];
    k_r = k_f / Kc;
    q_r[76] = phi_r * k_r;

    /*reaction 78: H + C2H6 <=> C2H5 + H2 */
    phi_f = sc[1]*sc[26];
    k_f = k_f_save[77];
    q_f[77] = phi_f * k_f;
    phi_r = sc[25]*sc[0];
    Kc = Kc_save[77];
    k_r = k_f / Kc;
    q_r[77] = phi_r * k_r;

    /*reaction 79: H + HCCO <=> CH2(S) + CO */
    phi_f = sc[1]*sc[27];
    k_f = k_f_save[78];
    q_f[78] = phi_f * k_f;
    phi_r = sc[11]*sc[14];
    Kc = Kc_save[78];
    k_r = k_f / Kc;
    q_r[78] = phi_r * k_r;

    /*reaction 80: H + CH2CO <=> HCCO + H2 */
    phi_f = sc[1]*sc[28];
    k_f = k_f_save[79];
    q_f[79] = phi_f * k_f;
    phi_r = sc[27]*sc[0];
    Kc = Kc_save[79];
    k_r = k_f / Kc;
    q_r[79] = phi_r * k_r;

    /*reaction 81: H + CH2CO <=> CH3 + CO */
    phi_f = sc[1]*sc[28];
    k_f = k_f_save[80];
    q_f[80] = phi_f * k_f;
    phi_r = sc[12]*sc[14];
    Kc = Kc_save[80];
    k_r = k_f / Kc;
    q_r[80] = phi_r * k_r;

    /*reaction 82: H + HCCOH <=> H + CH2CO */
    phi_f = sc[1]*sc[29];
    k_f = k_f_save[81];
    q_f[81] = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[81];
    k_r = k_f / Kc;
    q_r[81] = phi_r * k_r;

    /*reaction 83: H2 + CO (+M) <=> CH2O (+M) */
    phi_f = sc[0]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[82];
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
    q_f[82] = phi_f * k_f;
    phi_r = sc[17];
    Kc = Kc_save[82];
    k_r = k_f / Kc;
    q_r[82] = phi_r * k_r;

    /*reaction 84: OH + H2 <=> H + H2O */
    phi_f = sc[4]*sc[0];
    k_f = k_f_save[83];
    q_f[83] = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    Kc = Kc_save[83];
    k_r = k_f / Kc;
    q_r[83] = phi_r * k_r;

    /*reaction 85: 2 OH (+M) <=> H2O2 (+M) */
    phi_f = sc[4]*sc[4];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[84];
    redP = 1e-12 * alpha / k_f * 2.3e+18*exp(-0.9*tc[0]+855.46833186802348337*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2654*exp(T/-94))+ (0.7346*exp(T/-1756))+ (exp(-5182/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[84] = phi_f * k_f;
    phi_r = sc[7];
    Kc = Kc_save[84];
    k_r = k_f / Kc;
    q_r[84] = phi_r * k_r;

    /*reaction 86: 2 OH <=> O + H2O */
    phi_f = sc[4]*sc[4];
    k_f = k_f_save[85];
    q_f[85] = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = Kc_save[85];
    k_r = k_f / Kc;
    q_r[85] = phi_r * k_r;

    /*reaction 87: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = k_f_save[86];
    q_f[86] = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = Kc_save[86];
    k_r = k_f / Kc;
    q_r[86] = phi_r * k_r;

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    phi_f = sc[4]*sc[7];
    k_f = k_f_save[87];
    q_f[87] = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = Kc_save[87];
    k_r = k_f / Kc;
    q_r[87] = phi_r * k_r;

    /*reaction 89: OH + H2O2 <=> HO2 + H2O */
    phi_f = sc[4]*sc[7];
    k_f = k_f_save[88];
    q_f[88] = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = Kc_save[88];
    k_r = k_f / Kc;
    q_r[88] = phi_r * k_r;

    /*reaction 90: OH + C <=> H + CO */
    phi_f = sc[4]*sc[8];
    k_f = k_f_save[89];
    q_f[89] = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[89];
    k_r = k_f / Kc;
    q_r[89] = phi_r * k_r;

    /*reaction 91: OH + CH <=> H + HCO */
    phi_f = sc[4]*sc[9];
    k_f = k_f_save[90];
    q_f[90] = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_save[90];
    k_r = k_f / Kc;
    q_r[90] = phi_r * k_r;

    /*reaction 92: OH + CH2 <=> H + CH2O */
    phi_f = sc[4]*sc[10];
    k_f = k_f_save[91];
    q_f[91] = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[91];
    k_r = k_f / Kc;
    q_r[91] = phi_r * k_r;

    /*reaction 93: OH + CH2 <=> CH + H2O */
    phi_f = sc[4]*sc[10];
    k_f = k_f_save[92];
    q_f[92] = phi_f * k_f;
    phi_r = sc[9]*sc[5];
    Kc = Kc_save[92];
    k_r = k_f / Kc;
    q_r[92] = phi_r * k_r;

    /*reaction 94: OH + CH2(S) <=> H + CH2O */
    phi_f = sc[4]*sc[11];
    k_f = k_f_save[93];
    q_f[93] = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[93];
    k_r = k_f / Kc;
    q_r[93] = phi_r * k_r;

    /*reaction 95: OH + CH3 (+M) <=> CH3OH (+M) */
    phi_f = sc[4]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[94];
    redP = 1e-12 * alpha / k_f * 2.7e+38*exp(-6.3*tc[0]-1559.9716639946311716*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.7895*exp(T/-83.5))+ (0.2105*exp(T/-5398))+ (exp(-8370/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[94] = phi_f * k_f;
    phi_r = sc[20];
    Kc = Kc_save[94];
    k_r = k_f / Kc;
    q_r[94] = phi_r * k_r;

    /*reaction 96: OH + CH3 <=> CH2 + H2O */
    phi_f = sc[4]*sc[12];
    k_f = k_f_save[95];
    q_f[95] = phi_f * k_f;
    phi_r = sc[10]*sc[5];
    Kc = Kc_save[95];
    k_r = k_f / Kc;
    q_r[95] = phi_r * k_r;

    /*reaction 97: OH + CH3 <=> CH2(S) + H2O */
    phi_f = sc[4]*sc[12];
    k_f = k_f_save[96];
    q_f[96] = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[96];
    k_r = k_f / Kc;
    q_r[96] = phi_r * k_r;

    /*reaction 98: OH + CH4 <=> CH3 + H2O */
    phi_f = sc[4]*sc[13];
    k_f = k_f_save[97];
    q_f[97] = phi_f * k_f;
    phi_r = sc[12]*sc[5];
    Kc = Kc_save[97];
    k_r = k_f / Kc;
    q_r[97] = phi_r * k_r;

    /*reaction 99: OH + CO <=> H + CO2 */
    phi_f = sc[4]*sc[14];
    k_f = k_f_save[98];
    q_f[98] = phi_f * k_f;
    phi_r = sc[1]*sc[15];
    Kc = Kc_save[98];
    k_r = k_f / Kc;
    q_r[98] = phi_r * k_r;

    /*reaction 100: OH + HCO <=> H2O + CO */
    phi_f = sc[4]*sc[16];
    k_f = k_f_save[99];
    q_f[99] = phi_f * k_f;
    phi_r = sc[5]*sc[14];
    Kc = Kc_save[99];
    k_r = k_f / Kc;
    q_r[99] = phi_r * k_r;

    /*reaction 101: OH + CH2O <=> HCO + H2O */
    phi_f = sc[4]*sc[17];
    k_f = k_f_save[100];
    q_f[100] = phi_f * k_f;
    phi_r = sc[16]*sc[5];
    Kc = Kc_save[100];
    k_r = k_f / Kc;
    q_r[100] = phi_r * k_r;

    /*reaction 102: OH + CH2OH <=> H2O + CH2O */
    phi_f = sc[4]*sc[18];
    k_f = k_f_save[101];
    q_f[101] = phi_f * k_f;
    phi_r = sc[5]*sc[17];
    Kc = Kc_save[101];
    k_r = k_f / Kc;
    q_r[101] = phi_r * k_r;

    /*reaction 103: OH + CH3O <=> H2O + CH2O */
    phi_f = sc[4]*sc[19];
    k_f = k_f_save[102];
    q_f[102] = phi_f * k_f;
    phi_r = sc[5]*sc[17];
    Kc = Kc_save[102];
    k_r = k_f / Kc;
    q_r[102] = phi_r * k_r;

    /*reaction 104: OH + CH3OH <=> CH2OH + H2O */
    phi_f = sc[4]*sc[20];
    k_f = k_f_save[103];
    q_f[103] = phi_f * k_f;
    phi_r = sc[18]*sc[5];
    Kc = Kc_save[103];
    k_r = k_f / Kc;
    q_r[103] = phi_r * k_r;

    /*reaction 105: OH + CH3OH <=> CH3O + H2O */
    phi_f = sc[4]*sc[20];
    k_f = k_f_save[104];
    q_f[104] = phi_f * k_f;
    phi_r = sc[19]*sc[5];
    Kc = Kc_save[104];
    k_r = k_f / Kc;
    q_r[104] = phi_r * k_r;

    /*reaction 106: OH + C2H <=> H + HCCO */
    phi_f = sc[4]*sc[21];
    k_f = k_f_save[105];
    q_f[105] = phi_f * k_f;
    phi_r = sc[1]*sc[27];
    Kc = Kc_save[105];
    k_r = k_f / Kc;
    q_r[105] = phi_r * k_r;

    /*reaction 107: OH + C2H2 <=> H + CH2CO */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[106];
    q_f[106] = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[106];
    k_r = k_f / Kc;
    q_r[106] = phi_r * k_r;

    /*reaction 108: OH + C2H2 <=> H + HCCOH */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[107];
    q_f[107] = phi_f * k_f;
    phi_r = sc[1]*sc[29];
    Kc = Kc_save[107];
    k_r = k_f / Kc;
    q_r[107] = phi_r * k_r;

    /*reaction 109: OH + C2H2 <=> C2H + H2O */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[108];
    q_f[108] = phi_f * k_f;
    phi_r = sc[21]*sc[5];
    Kc = Kc_save[108];
    k_r = k_f / Kc;
    q_r[108] = phi_r * k_r;

    /*reaction 110: OH + C2H2 <=> CH3 + CO */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[109];
    q_f[109] = phi_f * k_f;
    phi_r = sc[12]*sc[14];
    Kc = Kc_save[109];
    k_r = k_f / Kc;
    q_r[109] = phi_r * k_r;

    /*reaction 111: OH + C2H3 <=> H2O + C2H2 */
    phi_f = sc[4]*sc[23];
    k_f = k_f_save[110];
    q_f[110] = phi_f * k_f;
    phi_r = sc[5]*sc[22];
    Kc = Kc_save[110];
    k_r = k_f / Kc;
    q_r[110] = phi_r * k_r;

    /*reaction 112: OH + C2H4 <=> C2H3 + H2O */
    phi_f = sc[4]*sc[24];
    k_f = k_f_save[111];
    q_f[111] = phi_f * k_f;
    phi_r = sc[23]*sc[5];
    Kc = Kc_save[111];
    k_r = k_f / Kc;
    q_r[111] = phi_r * k_r;

    /*reaction 113: OH + C2H6 <=> C2H5 + H2O */
    phi_f = sc[4]*sc[26];
    k_f = k_f_save[112];
    q_f[112] = phi_f * k_f;
    phi_r = sc[25]*sc[5];
    Kc = Kc_save[112];
    k_r = k_f / Kc;
    q_r[112] = phi_r * k_r;

    /*reaction 114: OH + CH2CO <=> HCCO + H2O */
    phi_f = sc[4]*sc[28];
    k_f = k_f_save[113];
    q_f[113] = phi_f * k_f;
    phi_r = sc[27]*sc[5];
    Kc = Kc_save[113];
    k_r = k_f / Kc;
    q_r[113] = phi_r * k_r;

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    phi_f = sc[6]*sc[6];
    k_f = k_f_save[114];
    q_f[114] = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = Kc_save[114];
    k_r = k_f / Kc;
    q_r[114] = phi_r * k_r;

    /*reaction 116: 2 HO2 <=> O2 + H2O2 */
    phi_f = sc[6]*sc[6];
    k_f = k_f_save[115];
    q_f[115] = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = Kc_save[115];
    k_r = k_f / Kc;
    q_r[115] = phi_r * k_r;

    /*reaction 117: HO2 + CH2 <=> OH + CH2O */
    phi_f = sc[6]*sc[10];
    k_f = k_f_save[116];
    q_f[116] = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[116];
    k_r = k_f / Kc;
    q_r[116] = phi_r * k_r;

    /*reaction 118: HO2 + CH3 <=> O2 + CH4 */
    phi_f = sc[6]*sc[12];
    k_f = k_f_save[117];
    q_f[117] = phi_f * k_f;
    phi_r = sc[3]*sc[13];
    Kc = Kc_save[117];
    k_r = k_f / Kc;
    q_r[117] = phi_r * k_r;

    /*reaction 119: HO2 + CH3 <=> OH + CH3O */
    phi_f = sc[6]*sc[12];
    k_f = k_f_save[118];
    q_f[118] = phi_f * k_f;
    phi_r = sc[4]*sc[19];
    Kc = Kc_save[118];
    k_r = k_f / Kc;
    q_r[118] = phi_r * k_r;

    /*reaction 120: HO2 + CO <=> OH + CO2 */
    phi_f = sc[6]*sc[14];
    k_f = k_f_save[119];
    q_f[119] = phi_f * k_f;
    phi_r = sc[4]*sc[15];
    Kc = Kc_save[119];
    k_r = k_f / Kc;
    q_r[119] = phi_r * k_r;

    /*reaction 121: HO2 + CH2O <=> HCO + H2O2 */
    phi_f = sc[6]*sc[17];
    k_f = k_f_save[120];
    q_f[120] = phi_f * k_f;
    phi_r = sc[16]*sc[7];
    Kc = Kc_save[120];
    k_r = k_f / Kc;
    q_r[120] = phi_r * k_r;

    /*reaction 122: C + O2 <=> O + CO */
    phi_f = sc[8]*sc[3];
    k_f = k_f_save[121];
    q_f[121] = phi_f * k_f;
    phi_r = sc[2]*sc[14];
    Kc = Kc_save[121];
    k_r = k_f / Kc;
    q_r[121] = phi_r * k_r;

    /*reaction 123: C + CH2 <=> H + C2H */
    phi_f = sc[8]*sc[10];
    k_f = k_f_save[122];
    q_f[122] = phi_f * k_f;
    phi_r = sc[1]*sc[21];
    Kc = Kc_save[122];
    k_r = k_f / Kc;
    q_r[122] = phi_r * k_r;

    /*reaction 124: C + CH3 <=> H + C2H2 */
    phi_f = sc[8]*sc[12];
    k_f = k_f_save[123];
    q_f[123] = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = Kc_save[123];
    k_r = k_f / Kc;
    q_r[123] = phi_r * k_r;

    /*reaction 125: CH + O2 <=> O + HCO */
    phi_f = sc[9]*sc[3];
    k_f = k_f_save[124];
    q_f[124] = phi_f * k_f;
    phi_r = sc[2]*sc[16];
    Kc = Kc_save[124];
    k_r = k_f / Kc;
    q_r[124] = phi_r * k_r;

    /*reaction 126: CH + H2 <=> H + CH2 */
    phi_f = sc[9]*sc[0];
    k_f = k_f_save[125];
    q_f[125] = phi_f * k_f;
    phi_r = sc[1]*sc[10];
    Kc = Kc_save[125];
    k_r = k_f / Kc;
    q_r[125] = phi_r * k_r;

    /*reaction 127: CH + H2O <=> H + CH2O */
    phi_f = sc[9]*sc[5];
    k_f = k_f_save[126];
    q_f[126] = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[126];
    k_r = k_f / Kc;
    q_r[126] = phi_r * k_r;

    /*reaction 128: CH + CH2 <=> H + C2H2 */
    phi_f = sc[9]*sc[10];
    k_f = k_f_save[127];
    q_f[127] = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = Kc_save[127];
    k_r = k_f / Kc;
    q_r[127] = phi_r * k_r;

    /*reaction 129: CH + CH3 <=> H + C2H3 */
    phi_f = sc[9]*sc[12];
    k_f = k_f_save[128];
    q_f[128] = phi_f * k_f;
    phi_r = sc[1]*sc[23];
    Kc = Kc_save[128];
    k_r = k_f / Kc;
    q_r[128] = phi_r * k_r;

    /*reaction 130: CH + CH4 <=> H + C2H4 */
    phi_f = sc[9]*sc[13];
    k_f = k_f_save[129];
    q_f[129] = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = Kc_save[129];
    k_r = k_f / Kc;
    q_r[129] = phi_r * k_r;

    /*reaction 131: CH + CO (+M) <=> HCCO (+M) */
    phi_f = sc[9]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[130];
    redP = 1e-12 * alpha / k_f * 2.69e+28*exp(-3.74*tc[0]-974.227464997937318*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.4243*exp(T/-237))+ (0.5757*exp(T/-1652))+ (exp(-5069/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[130] = phi_f * k_f;
    phi_r = sc[27];
    Kc = Kc_save[130];
    k_r = k_f / Kc;
    q_r[130] = phi_r * k_r;

    /*reaction 132: CH + CO2 <=> HCO + CO */
    phi_f = sc[9]*sc[15];
    k_f = k_f_save[131];
    q_f[131] = phi_f * k_f;
    phi_r = sc[16]*sc[14];
    Kc = Kc_save[131];
    k_r = k_f / Kc;
    q_r[131] = phi_r * k_r;

    /*reaction 133: CH + CH2O <=> H + CH2CO */
    phi_f = sc[9]*sc[17];
    k_f = k_f_save[132];
    q_f[132] = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[132];
    k_r = k_f / Kc;
    q_r[132] = phi_r * k_r;

    /*reaction 134: CH + HCCO <=> CO + C2H2 */
    phi_f = sc[9]*sc[27];
    k_f = k_f_save[133];
    q_f[133] = phi_f * k_f;
    phi_r = sc[14]*sc[22];
    Kc = Kc_save[133];
    k_r = k_f / Kc;
    q_r[133] = phi_r * k_r;

    /*reaction 135: CH2 + O2 <=> OH + HCO */
    phi_f = sc[10]*sc[3];
    k_f = k_f_save[134];
    q_f[134] = phi_f * k_f;
    phi_r = sc[4]*sc[16];
    Kc = Kc_save[134];
    k_r = k_f / Kc;
    q_r[134] = phi_r * k_r;

    /*reaction 136: CH2 + H2 <=> H + CH3 */
    phi_f = sc[10]*sc[0];
    k_f = k_f_save[135];
    q_f[135] = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = Kc_save[135];
    k_r = k_f / Kc;
    q_r[135] = phi_r * k_r;

    /*reaction 137: 2 CH2 <=> H2 + C2H2 */
    phi_f = sc[10]*sc[10];
    k_f = k_f_save[136];
    q_f[136] = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = Kc_save[136];
    k_r = k_f / Kc;
    q_r[136] = phi_r * k_r;

    /*reaction 138: CH2 + CH3 <=> H + C2H4 */
    phi_f = sc[10]*sc[12];
    k_f = k_f_save[137];
    q_f[137] = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = Kc_save[137];
    k_r = k_f / Kc;
    q_r[137] = phi_r * k_r;

    /*reaction 139: CH2 + CH4 <=> 2 CH3 */
    phi_f = sc[10]*sc[13];
    k_f = k_f_save[138];
    q_f[138] = phi_f * k_f;
    phi_r = sc[12]*sc[12];
    Kc = Kc_save[138];
    k_r = k_f / Kc;
    q_r[138] = phi_r * k_r;

    /*reaction 140: CH2 + CO (+M) <=> CH2CO (+M) */
    phi_f = sc[10]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[139];
    redP = 1e-12 * alpha / k_f * 2.69e+33*exp(-5.11*tc[0]-3570.3222438844863973*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.4093*exp(T/-275))+ (0.5907*exp(T/-1226))+ (exp(-5185/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[139] = phi_f * k_f;
    phi_r = sc[28];
    Kc = Kc_save[139];
    k_r = k_f / Kc;
    q_r[139] = phi_r * k_r;

    /*reaction 141: CH2 + HCCO <=> C2H3 + CO */
    phi_f = sc[10]*sc[27];
    k_f = k_f_save[140];
    q_f[140] = phi_f * k_f;
    phi_r = sc[23]*sc[14];
    Kc = Kc_save[140];
    k_r = k_f / Kc;
    q_r[140] = phi_r * k_r;

    /*reaction 142: CH2(S) + N2 <=> CH2 + N2 */
    phi_f = sc[11]*sc[30];
    k_f = k_f_save[141];
    q_f[141] = phi_f * k_f;
    phi_r = sc[10]*sc[30];
    Kc = Kc_save[141];
    k_r = k_f / Kc;
    q_r[141] = phi_r * k_r;

    /*reaction 143: CH2(S) + AR <=> CH2 + AR */
    phi_f = sc[11]*sc[31];
    k_f = k_f_save[142];
    q_f[142] = phi_f * k_f;
    phi_r = sc[10]*sc[31];
    Kc = Kc_save[142];
    k_r = k_f / Kc;
    q_r[142] = phi_r * k_r;

    /*reaction 144: CH2(S) + O2 <=> H + OH + CO */
    phi_f = sc[11]*sc[3];
    k_f = k_f_save[143];
    q_f[143] = phi_f * k_f;
    phi_r = sc[1]*sc[4]*sc[14];
    Kc = Kc_save[143];
    k_r = k_f / Kc;
    q_r[143] = phi_r * k_r;

    /*reaction 145: CH2(S) + O2 <=> CO + H2O */
    phi_f = sc[11]*sc[3];
    k_f = k_f_save[144];
    q_f[144] = phi_f * k_f;
    phi_r = sc[14]*sc[5];
    Kc = Kc_save[144];
    k_r = k_f / Kc;
    q_r[144] = phi_r * k_r;

    /*reaction 146: CH2(S) + H2 <=> CH3 + H */
    phi_f = sc[11]*sc[0];
    k_f = k_f_save[145];
    q_f[145] = phi_f * k_f;
    phi_r = sc[12]*sc[1];
    Kc = Kc_save[145];
    k_r = k_f / Kc;
    q_r[145] = phi_r * k_r;

    /*reaction 147: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    phi_f = sc[11]*sc[5];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[146];
    redP = 1e-12 * alpha / k_f * 2.7e+38*exp(-6.3*tc[0]-1559.9716639946311716*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.8493*exp(T/-134))+ (0.1507*exp(T/-2383))+ (exp(-7265/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[146] = phi_f * k_f;
    phi_r = sc[20];
    Kc = Kc_save[146];
    k_r = k_f / Kc;
    q_r[146] = phi_r * k_r;

    /*reaction 148: CH2(S) + H2O <=> CH2 + H2O */
    phi_f = sc[11]*sc[5];
    k_f = k_f_save[147];
    q_f[147] = phi_f * k_f;
    phi_r = sc[10]*sc[5];
    Kc = Kc_save[147];
    k_r = k_f / Kc;
    q_r[147] = phi_r * k_r;

    /*reaction 149: CH2(S) + CH3 <=> H + C2H4 */
    phi_f = sc[11]*sc[12];
    k_f = k_f_save[148];
    q_f[148] = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = Kc_save[148];
    k_r = k_f / Kc;
    q_r[148] = phi_r * k_r;

    /*reaction 150: CH2(S) + CH4 <=> 2 CH3 */
    phi_f = sc[11]*sc[13];
    k_f = k_f_save[149];
    q_f[149] = phi_f * k_f;
    phi_r = sc[12]*sc[12];
    Kc = Kc_save[149];
    k_r = k_f / Kc;
    q_r[149] = phi_r * k_r;

    /*reaction 151: CH2(S) + CO <=> CH2 + CO */
    phi_f = sc[11]*sc[14];
    k_f = k_f_save[150];
    q_f[150] = phi_f * k_f;
    phi_r = sc[10]*sc[14];
    Kc = Kc_save[150];
    k_r = k_f / Kc;
    q_r[150] = phi_r * k_r;

    /*reaction 152: CH2(S) + CO2 <=> CH2 + CO2 */
    phi_f = sc[11]*sc[15];
    k_f = k_f_save[151];
    q_f[151] = phi_f * k_f;
    phi_r = sc[10]*sc[15];
    Kc = Kc_save[151];
    k_r = k_f / Kc;
    q_r[151] = phi_r * k_r;

    /*reaction 153: CH2(S) + CO2 <=> CO + CH2O */
    phi_f = sc[11]*sc[15];
    k_f = k_f_save[152];
    q_f[152] = phi_f * k_f;
    phi_r = sc[14]*sc[17];
    Kc = Kc_save[152];
    k_r = k_f / Kc;
    q_r[152] = phi_r * k_r;

    /*reaction 154: CH2(S) + C2H6 <=> CH3 + C2H5 */
    phi_f = sc[11]*sc[26];
    k_f = k_f_save[153];
    q_f[153] = phi_f * k_f;
    phi_r = sc[12]*sc[25];
    Kc = Kc_save[153];
    k_r = k_f / Kc;
    q_r[153] = phi_r * k_r;

    /*reaction 155: CH3 + O2 <=> O + CH3O */
    phi_f = sc[12]*sc[3];
    k_f = k_f_save[154];
    q_f[154] = phi_f * k_f;
    phi_r = sc[2]*sc[19];
    Kc = Kc_save[154];
    k_r = k_f / Kc;
    q_r[154] = phi_r * k_r;

    /*reaction 156: CH3 + O2 <=> OH + CH2O */
    phi_f = sc[12]*sc[3];
    k_f = k_f_save[155];
    q_f[155] = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[155];
    k_r = k_f / Kc;
    q_r[155] = phi_r * k_r;

    /*reaction 157: CH3 + H2O2 <=> HO2 + CH4 */
    phi_f = sc[12]*sc[7];
    k_f = k_f_save[156];
    q_f[156] = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = Kc_save[156];
    k_r = k_f / Kc;
    q_r[156] = phi_r * k_r;

    /*reaction 158: 2 CH3 (+M) <=> C2H6 (+M) */
    phi_f = sc[12]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[157];
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
    q_f[157] = phi_f * k_f;
    phi_r = sc[26];
    Kc = Kc_save[157];
    k_r = k_f / Kc;
    q_r[157] = phi_r * k_r;

    /*reaction 159: 2 CH3 <=> H + C2H5 */
    phi_f = sc[12]*sc[12];
    k_f = k_f_save[158];
    q_f[158] = phi_f * k_f;
    phi_r = sc[1]*sc[25];
    Kc = Kc_save[158];
    k_r = k_f / Kc;
    q_r[158] = phi_r * k_r;

    /*reaction 160: CH3 + HCO <=> CH4 + CO */
    phi_f = sc[12]*sc[16];
    k_f = k_f_save[159];
    q_f[159] = phi_f * k_f;
    phi_r = sc[13]*sc[14];
    Kc = Kc_save[159];
    k_r = k_f / Kc;
    q_r[159] = phi_r * k_r;

    /*reaction 161: CH3 + CH2O <=> HCO + CH4 */
    phi_f = sc[12]*sc[17];
    k_f = k_f_save[160];
    q_f[160] = phi_f * k_f;
    phi_r = sc[16]*sc[13];
    Kc = Kc_save[160];
    k_r = k_f / Kc;
    q_r[160] = phi_r * k_r;

    /*reaction 162: CH3 + CH3OH <=> CH2OH + CH4 */
    phi_f = sc[12]*sc[20];
    k_f = k_f_save[161];
    q_f[161] = phi_f * k_f;
    phi_r = sc[18]*sc[13];
    Kc = Kc_save[161];
    k_r = k_f / Kc;
    q_r[161] = phi_r * k_r;

    /*reaction 163: CH3 + CH3OH <=> CH3O + CH4 */
    phi_f = sc[12]*sc[20];
    k_f = k_f_save[162];
    q_f[162] = phi_f * k_f;
    phi_r = sc[19]*sc[13];
    Kc = Kc_save[162];
    k_r = k_f / Kc;
    q_r[162] = phi_r * k_r;

    /*reaction 164: CH3 + C2H4 <=> C2H3 + CH4 */
    phi_f = sc[12]*sc[24];
    k_f = k_f_save[163];
    q_f[163] = phi_f * k_f;
    phi_r = sc[23]*sc[13];
    Kc = Kc_save[163];
    k_r = k_f / Kc;
    q_r[163] = phi_r * k_r;

    /*reaction 165: CH3 + C2H6 <=> C2H5 + CH4 */
    phi_f = sc[12]*sc[26];
    k_f = k_f_save[164];
    q_f[164] = phi_f * k_f;
    phi_r = sc[25]*sc[13];
    Kc = Kc_save[164];
    k_r = k_f / Kc;
    q_r[164] = phi_r * k_r;

    /*reaction 166: HCO + H2O <=> H + CO + H2O */
    phi_f = sc[16]*sc[5];
    k_f = k_f_save[165];
    q_f[165] = phi_f * k_f;
    phi_r = sc[1]*sc[14]*sc[5];
    Kc = Kc_save[165];
    k_r = k_f / Kc;
    q_r[165] = phi_r * k_r;

    /*reaction 167: HCO + M <=> H + CO + M */
    phi_f = sc[16];
    alpha = mixture + sc[0] + -1*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = alpha * k_f_save[166];
    q_f[166] = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[166];
    k_r = k_f / Kc;
    q_r[166] = phi_r * k_r;

    /*reaction 168: HCO + O2 <=> HO2 + CO */
    phi_f = sc[16]*sc[3];
    k_f = k_f_save[167];
    q_f[167] = phi_f * k_f;
    phi_r = sc[6]*sc[14];
    Kc = Kc_save[167];
    k_r = k_f / Kc;
    q_r[167] = phi_r * k_r;

    /*reaction 169: CH2OH + O2 <=> HO2 + CH2O */
    phi_f = sc[18]*sc[3];
    k_f = k_f_save[168];
    q_f[168] = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = Kc_save[168];
    k_r = k_f / Kc;
    q_r[168] = phi_r * k_r;

    /*reaction 170: CH3O + O2 <=> HO2 + CH2O */
    phi_f = sc[19]*sc[3];
    k_f = k_f_save[169];
    q_f[169] = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = Kc_save[169];
    k_r = k_f / Kc;
    q_r[169] = phi_r * k_r;

    /*reaction 171: C2H + O2 <=> HCO + CO */
    phi_f = sc[21]*sc[3];
    k_f = k_f_save[170];
    q_f[170] = phi_f * k_f;
    phi_r = sc[16]*sc[14];
    Kc = Kc_save[170];
    k_r = k_f / Kc;
    q_r[170] = phi_r * k_r;

    /*reaction 172: C2H + H2 <=> H + C2H2 */
    phi_f = sc[21]*sc[0];
    k_f = k_f_save[171];
    q_f[171] = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = Kc_save[171];
    k_r = k_f / Kc;
    q_r[171] = phi_r * k_r;

    /*reaction 173: C2H3 + O2 <=> HCO + CH2O */
    phi_f = sc[23]*sc[3];
    k_f = k_f_save[172];
    q_f[172] = phi_f * k_f;
    phi_r = sc[16]*sc[17];
    Kc = Kc_save[172];
    k_r = k_f / Kc;
    q_r[172] = phi_r * k_r;

    /*reaction 174: C2H4 (+M) <=> H2 + C2H2 (+M) */
    phi_f = sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26] + -0.3*sc[31];
    k_f = k_f_save[173];
    redP = 1e-6 * alpha / k_f * 7e+50*exp(-9.31*tc[0]-50251.216247259304509*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2655*exp(T/-180))+ (0.7345*exp(T/-1035))+ (exp(-5417/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[173] = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = Kc_save[173];
    k_r = k_f / Kc;
    q_r[173] = phi_r * k_r;

    /*reaction 175: C2H5 + O2 <=> HO2 + C2H4 */
    phi_f = sc[25]*sc[3];
    k_f = k_f_save[174];
    q_f[174] = phi_f * k_f;
    phi_r = sc[6]*sc[24];
    Kc = Kc_save[174];
    k_r = k_f / Kc;
    q_r[174] = phi_r * k_r;

    /*reaction 176: HCCO + O2 <=> OH + 2 CO */
    phi_f = sc[27]*sc[3];
    k_f = k_f_save[175];
    q_f[175] = phi_f * k_f;
    phi_r = sc[4]*sc[14]*sc[14];
    Kc = Kc_save[175];
    k_r = k_f / Kc;
    q_r[175] = phi_r * k_r;

    /*reaction 177: 2 HCCO <=> 2 CO + C2H2 */
    phi_f = sc[27]*sc[27];
    k_f = k_f_save[176];
    q_f[176] = phi_f * k_f;
    phi_r = sc[14]*sc[14]*sc[22];
    Kc = Kc_save[176];
    k_r = k_f / Kc;
    q_r[176] = phi_r * k_r;

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
        /*species 7: H2O2 */
        species[7] =
            -1.770258210000000e+04 * invT
            +8.410619499999998e-01
            -4.276112690000000e+00 * tc[0]
            +2.714112085000000e-04 * tc[1]
            -2.788928350000000e-06 * tc[2]
            +1.798090108333333e-09 * tc[3]
            -4.312271815000000e-13 * tc[4];
        /*species 8: C */
        species[8] =
            +8.544388320000000e+04 * invT
            -1.977068930000000e+00
            -2.554239550000000e+00 * tc[0]
            +1.607688620000000e-04 * tc[1]
            -1.222987075000000e-07 * tc[2]
            +6.101957408333333e-11 * tc[3]
            -1.332607230000000e-14 * tc[4];
        /*species 9: CH */
        species[9] =
            +7.079729340000000e+04 * invT
            +1.405805570000000e+00
            -3.489816650000000e+00 * tc[0]
            -1.619177705000000e-04 * tc[1]
            +2.814984416666667e-07 * tc[2]
            -2.635144391666666e-10 * tc[3]
            +7.030453350000001e-14 * tc[4];
        /*species 10: CH2 */
        species[10] =
            +4.600404010000000e+04 * invT
            +2.200146820000000e+00
            -3.762678670000000e+00 * tc[0]
            -4.844360715000000e-04 * tc[1]
            -4.658164016666667e-07 * tc[2]
            +3.209092941666667e-10 * tc[3]
            -8.437085950000000e-14 * tc[4];
        /*species 11: CH2(S) */
        species[11] =
            +5.049681630000000e+04 * invT
            +4.967723077000000e+00
            -4.198604110000000e+00 * tc[0]
            +1.183307095000000e-03 * tc[1]
            -1.372160366666667e-06 * tc[2]
            +5.573466508333334e-10 * tc[3]
            -9.715736850000000e-14 * tc[4];
        /*species 12: CH3 */
        species[12] =
            +1.644499880000000e+04 * invT
            +2.069026070000000e+00
            -3.673590400000000e+00 * tc[0]
            -1.005475875000000e-03 * tc[1]
            -9.550364266666668e-07 * tc[2]
            +5.725978541666666e-10 * tc[3]
            -1.271928670000000e-13 * tc[4];
        /*species 13: CH4 */
        species[13] =
            -1.024664760000000e+04 * invT
            +9.791179889999999e+00
            -5.149876130000000e+00 * tc[0]
            +6.835489400000000e-03 * tc[1]
            -8.196676650000000e-06 * tc[2]
            +4.039525216666667e-09 * tc[3]
            -8.334697800000000e-13 * tc[4];
        /*species 14: CO */
        species[14] =
            -1.434408600000000e+04 * invT
            +7.112418999999992e-02
            -3.579533470000000e+00 * tc[0]
            +3.051768400000000e-04 * tc[1]
            -1.694690550000000e-07 * tc[2]
            -7.558382366666667e-11 * tc[3]
            +4.522122495000000e-14 * tc[4];
        /*species 15: CO2 */
        species[15] =
            -4.837196970000000e+04 * invT
            -7.544278700000000e+00
            -2.356773520000000e+00 * tc[0]
            -4.492298385000000e-03 * tc[1]
            +1.187260448333333e-06 * tc[2]
            -2.049325183333333e-10 * tc[3]
            +7.184977399999999e-15 * tc[4];
        /*species 16: HCO */
        species[16] =
            +3.839564960000000e+03 * invT
            +8.268134100000002e-01
            -4.221185840000000e+00 * tc[0]
            +1.621962660000000e-03 * tc[1]
            -2.296657433333333e-06 * tc[2]
            +1.109534108333333e-09 * tc[3]
            -2.168844325000000e-13 * tc[4];
        /*species 17: CH2O */
        species[17] =
            -1.430895670000000e+04 * invT
            +4.190910250000000e+00
            -4.793723150000000e+00 * tc[0]
            +4.954166845000000e-03 * tc[1]
            -6.220333466666666e-06 * tc[2]
            +3.160710508333333e-09 * tc[3]
            -6.588632600000000e-13 * tc[4];
        /*species 18: CH2OH */
        species[18] =
            -3.193913670000000e+03 * invT
            -1.609133250000000e+00
            -3.863889180000000e+00 * tc[0]
            -2.798361520000000e-03 * tc[1]
            -9.887863183333334e-07 * tc[2]
            +8.711001000000001e-10 * tc[3]
            -2.184836390000000e-13 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +9.786011000000000e+02 * invT
            -1.104597300000000e+01
            -2.106204000000000e+00 * tc[0]
            -3.608297500000000e-03 * tc[1]
            -8.897453333333333e-07 * tc[2]
            +6.148030000000000e-10 * tc[3]
            -1.037805000000000e-13 * tc[4];
        /*species 20: CH3OH */
        species[20] =
            -2.564276560000000e+04 * invT
            +7.219494050000001e+00
            -5.715395820000000e+00 * tc[0]
            +7.615456450000000e-03 * tc[1]
            -1.087401925000000e-05 * tc[2]
            +5.923390741666667e-09 * tc[3]
            -1.306763490000000e-12 * tc[4];
        /*species 21: C2H */
        species[21] =
            +6.683939320000001e+04 * invT
            -3.333307050000000e+00
            -2.889657330000000e+00 * tc[0]
            -6.704980550000000e-03 * tc[1]
            +4.746158350000000e-06 * tc[2]
            -2.456592041666667e-09 * tc[3]
            +5.466575550000000e-13 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            +2.642898070000000e+04 * invT
            -1.313102400600000e+01
            -8.086810940000000e-01 * tc[0]
            -1.168078145000000e-02 * tc[1]
            +5.919530250000000e-06 * tc[2]
            -2.334603641666667e-09 * tc[3]
            +4.250364870000000e-13 * tc[4];
        /*species 23: C2H3 */
        species[23] =
            +3.485984680000000e+04 * invT
            -5.298073800000000e+00
            -3.212466450000000e+00 * tc[0]
            -7.573958100000000e-04 * tc[1]
            -4.320156866666666e-06 * tc[2]
            +2.980482058333333e-09 * tc[3]
            -7.357543650000000e-13 * tc[4];
        /*species 24: C2H4 */
        species[24] =
            +5.089775930000000e+03 * invT
            -1.381294799999999e-01
            -3.959201480000000e+00 * tc[0]
            +3.785261235000000e-03 * tc[1]
            -9.516504866666667e-06 * tc[2]
            +5.763239608333333e-09 * tc[3]
            -1.349421865000000e-12 * tc[4];
        /*species 25: C2H5 */
        species[25] =
            +1.284162650000000e+04 * invT
            -4.007435600000004e-01
            -4.306465680000000e+00 * tc[0]
            +2.093294460000000e-03 * tc[1]
            -8.285713450000000e-06 * tc[2]
            +4.992721716666666e-09 * tc[3]
            -1.152545020000000e-12 * tc[4];
        /*species 26: C2H6 */
        species[26] =
            -1.152220550000000e+04 * invT
            +1.624601760000000e+00
            -4.291424920000000e+00 * tc[0]
            +2.750771350000000e-03 * tc[1]
            -9.990638133333334e-06 * tc[2]
            +5.903885708333334e-09 * tc[3]
            -1.343428855000000e-12 * tc[4];
        /*species 27: HCCO */
        species[27] =
            +2.005944900000000e+04 * invT
            -1.023869560000000e+01
            -2.251721400000000e+00 * tc[0]
            -8.827510500000000e-03 * tc[1]
            +3.954850166666666e-06 * tc[2]
            -1.439646583333334e-09 * tc[3]
            +2.533240550000000e-13 * tc[4];
        /*species 28: CH2CO */
        species[28] =
            -7.042918040000000e+03 * invT
            -1.007981170000000e+01
            -2.135836300000000e+00 * tc[0]
            -9.059436050000000e-03 * tc[1]
            +2.899124566666666e-06 * tc[2]
            -7.786646400000000e-10 * tc[3]
            +1.007288075000000e-13 * tc[4];
        /*species 29: HCCOH */
        species[29] =
            +8.031614300000000e+03 * invT
            -1.263194570000000e+01
            -1.242373300000000e+00 * tc[0]
            -1.553610050000000e-02 * tc[1]
            +8.477810666666667e-06 * tc[2]
            -3.594760916666667e-09 * tc[3]
            +7.007297000000000e-13 * tc[4];
        /*species 30: N2 */
        species[30] =
            -1.020899900000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041202000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427000000000e-13 * tc[4];
        /*species 31: AR */
        species[31] =
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
        /*species 7: H2O2 */
        species[7] =
            -1.786178770000000e+04 * invT
            +1.248846229999999e+00
            -4.165002850000000e+00 * tc[0]
            -2.454158470000000e-03 * tc[1]
            +3.168987083333333e-07 * tc[2]
            -3.093216550000000e-11 * tc[3]
            +1.439541525000000e-15 * tc[4];
        /*species 8: C */
        species[8] =
            +8.545129530000000e+04 * invT
            -2.308834850000000e+00
            -2.492668880000000e+00 * tc[0]
            -2.399446420000000e-05 * tc[1]
            +1.207225033333333e-08 * tc[2]
            -3.119091908333333e-12 * tc[3]
            +2.436389465000000e-16 * tc[4];
        /*species 9: CH */
        species[9] =
            +7.101243640000001e+04 * invT
            -2.606515260000000e+00
            -2.878464730000000e+00 * tc[0]
            -4.854568405000000e-04 * tc[1]
            -2.407427583333333e-08 * tc[2]
            +1.089065408333333e-11 * tc[3]
            -8.803969149999999e-16 * tc[4];
        /*species 10: CH2 */
        species[10] =
            +4.626360400000000e+04 * invT
            -3.297092110000000e+00
            -2.874101130000000e+00 * tc[0]
            -1.828196460000000e-03 * tc[1]
            +2.348243283333333e-07 * tc[2]
            -2.168162908333333e-11 * tc[3]
            +9.386378350000000e-16 * tc[4];
        /*species 11: CH2(S) */
        species[11] =
            +5.092599970000000e+04 * invT
            -6.334463270000000e+00
            -2.292038420000000e+00 * tc[0]
            -2.327943185000000e-03 * tc[1]
            +3.353199116666667e-07 * tc[2]
            -3.482550000000000e-11 * tc[3]
            +1.698581825000000e-15 * tc[4];
        /*species 12: CH3 */
        species[12] =
            +1.677558430000000e+04 * invT
            -6.194354070000000e+00
            -2.285717720000000e+00 * tc[0]
            -3.619950185000000e-03 * tc[1]
            +4.978572466666667e-07 * tc[2]
            -4.964038700000000e-11 * tc[3]
            +2.335771970000000e-15 * tc[4];
        /*species 13: CH4 */
        species[13] =
            -9.468344590000001e+03 * invT
            -1.836246650500000e+01
            -7.485149500000000e-02 * tc[0]
            -6.695473350000000e-03 * tc[1]
            +9.554763483333333e-07 * tc[2]
            -1.019104458333333e-10 * tc[3]
            +5.090761500000000e-15 * tc[4];
        /*species 14: CO */
        species[14] =
            -1.415187240000000e+04 * invT
            -5.103502110000000e+00
            -2.715185610000000e+00 * tc[0]
            -1.031263715000000e-03 * tc[1]
            +1.664709618333334e-07 * tc[2]
            -1.917108400000000e-11 * tc[3]
            +1.018238580000000e-15 * tc[4];
        /*species 15: CO2 */
        species[15] =
            -4.875916600000000e+04 * invT
            +1.585822230000000e+00
            -3.857460290000000e+00 * tc[0]
            -2.207185130000000e-03 * tc[1]
            +3.691356733333334e-07 * tc[2]
            -4.362418233333334e-11 * tc[3]
            +2.360420820000000e-15 * tc[4];
        /*species 16: HCO */
        species[16] =
            +4.011918150000000e+03 * invT
            -7.026170540000000e+00
            -2.772174380000000e+00 * tc[0]
            -2.478477630000000e-03 * tc[1]
            +4.140760216666667e-07 * tc[2]
            -4.909681483333334e-11 * tc[3]
            +2.667543555000000e-15 * tc[4];
        /*species 17: CH2O */
        species[17] =
            -1.399583230000000e+04 * invT
            -1.189563292000000e+01
            -1.760690080000000e+00 * tc[0]
            -4.600000410000000e-03 * tc[1]
            +7.370980216666666e-07 * tc[2]
            -8.386767666666666e-11 * tc[3]
            +4.419278200000001e-15 * tc[4];
        /*species 18: CH2OH */
        species[18] =
            -3.242506270000000e+03 * invT
            -2.117766459999999e+00
            -3.692665690000000e+00 * tc[0]
            -4.322883985000000e-03 * tc[1]
            +6.251685333333334e-07 * tc[2]
            -6.560288633333334e-11 * tc[3]
            +3.242771005000000e-15 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +1.278325200000000e+02 * invT
            +8.412240000000000e-01
            -3.770799000000000e+00 * tc[0]
            -3.935748500000000e-03 * tc[1]
            +4.427306666666667e-07 * tc[2]
            -3.287025833333333e-11 * tc[3]
            +1.056308000000000e-15 * tc[4];
        /*species 20: CH3OH */
        species[20] =
            -2.537487470000000e+04 * invT
            -1.271265439000000e+01
            -1.789707910000000e+00 * tc[0]
            -7.046914600000000e-03 * tc[1]
            +1.060834725000000e-06 * tc[2]
            -1.151425708333333e-10 * tc[3]
            +5.853011000000000e-15 * tc[4];
        /*species 21: C2H */
        species[21] =
            +6.712106500000000e+04 * invT
            -3.468088230000000e+00
            -3.167806520000000e+00 * tc[0]
            -2.376109510000000e-03 * tc[1]
            +3.063117950000000e-07 * tc[2]
            -2.534918766666666e-11 * tc[3]
            +8.861638500000000e-16 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            +2.593599920000000e+04 * invT
            +5.377850850000001e+00
            -4.147569640000000e+00 * tc[0]
            -2.980833320000000e-03 * tc[1]
            +3.954914200000000e-07 * tc[2]
            -3.895101425000000e-11 * tc[3]
            +1.806176065000000e-15 * tc[4];
        /*species 23: C2H3 */
        species[23] =
            +3.461287390000000e+04 * invT
            -4.770599780000000e+00
            -3.016724000000000e+00 * tc[0]
            -5.165114600000000e-03 * tc[1]
            +7.801372483333333e-07 * tc[2]
            -8.480274000000000e-11 * tc[3]
            +4.313035205000000e-15 * tc[4];
        /*species 24: C2H4 */
        species[24] =
            +4.939886140000000e+03 * invT
            -8.269258140000002e+00
            -2.036111160000000e+00 * tc[0]
            -7.322707550000000e-03 * tc[1]
            +1.118463191666667e-06 * tc[2]
            -1.226857691666667e-10 * tc[3]
            +6.285303050000000e-15 * tc[4];
        /*species 25: C2H5 */
        species[25] =
            +1.285752000000000e+04 * invT
            -1.150777788000000e+01
            -1.954656420000000e+00 * tc[0]
            -8.698636100000001e-03 * tc[1]
            +1.330344446666667e-06 * tc[2]
            -1.460147408333333e-10 * tc[3]
            +7.482078800000000e-15 * tc[4];
        /*species 26: C2H6 */
        species[26] =
            -1.142639320000000e+04 * invT
            -1.404372920000000e+01
            -1.071881500000000e+00 * tc[0]
            -1.084263385000000e-02 * tc[1]
            +1.670934450000000e-06 * tc[2]
            -1.845100008333333e-10 * tc[3]
            +9.500144500000000e-15 * tc[4];
        /*species 27: HCCO */
        species[27] =
            +1.932721500000000e+04 * invT
            +9.558465300000000e+00
            -5.628205800000000e+00 * tc[0]
            -2.042670050000000e-03 * tc[1]
            +2.655757833333333e-07 * tc[2]
            -2.385504333333333e-11 * tc[3]
            +9.703915999999999e-16 * tc[4];
        /*species 28: CH2CO */
        species[28] =
            -7.551053110000000e+03 * invT
            +3.879050115000000e+00
            -4.511297320000000e+00 * tc[0]
            -4.501798725000000e-03 * tc[1]
            +6.948993916666666e-07 * tc[2]
            -7.694549016666667e-11 * tc[3]
            +3.974191005000000e-15 * tc[4];
        /*species 29: HCCOH */
        species[29] =
            +7.264626000000000e+03 * invT
            +1.352560330000000e+01
            -5.923829100000000e+00 * tc[0]
            -3.396180000000000e-03 * tc[1]
            +4.276427333333333e-07 * tc[2]
            -3.748986750000000e-11 * tc[3]
            +1.497005050000000e-15 * tc[4];
        /*species 30: N2 */
        species[30] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439884000000000e-04 * tc[1]
            +9.474600000000001e-08 * tc[2]
            -8.414198333333333e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
        /*species 31: AR */
        species[31] =
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
        /*species 7: H2O2 */
        species[7] =
            -1.77025821e+04 * invT
            -1.58938050e-01
            -4.27611269e+00 * tc[0]
            +2.71411208e-04 * tc[1]
            -2.78892835e-06 * tc[2]
            +1.79809011e-09 * tc[3]
            -4.31227182e-13 * tc[4];
        /*species 8: C */
        species[8] =
            +8.54438832e+04 * invT
            -2.97706893e+00
            -2.55423955e+00 * tc[0]
            +1.60768862e-04 * tc[1]
            -1.22298707e-07 * tc[2]
            +6.10195741e-11 * tc[3]
            -1.33260723e-14 * tc[4];
        /*species 9: CH */
        species[9] =
            +7.07972934e+04 * invT
            +4.05805570e-01
            -3.48981665e+00 * tc[0]
            -1.61917771e-04 * tc[1]
            +2.81498442e-07 * tc[2]
            -2.63514439e-10 * tc[3]
            +7.03045335e-14 * tc[4];
        /*species 10: CH2 */
        species[10] =
            +4.60040401e+04 * invT
            +1.20014682e+00
            -3.76267867e+00 * tc[0]
            -4.84436072e-04 * tc[1]
            -4.65816402e-07 * tc[2]
            +3.20909294e-10 * tc[3]
            -8.43708595e-14 * tc[4];
        /*species 11: CH2(S) */
        species[11] =
            +5.04968163e+04 * invT
            +3.96772308e+00
            -4.19860411e+00 * tc[0]
            +1.18330710e-03 * tc[1]
            -1.37216037e-06 * tc[2]
            +5.57346651e-10 * tc[3]
            -9.71573685e-14 * tc[4];
        /*species 12: CH3 */
        species[12] =
            +1.64449988e+04 * invT
            +1.06902607e+00
            -3.67359040e+00 * tc[0]
            -1.00547588e-03 * tc[1]
            -9.55036427e-07 * tc[2]
            +5.72597854e-10 * tc[3]
            -1.27192867e-13 * tc[4];
        /*species 13: CH4 */
        species[13] =
            -1.02466476e+04 * invT
            +8.79117989e+00
            -5.14987613e+00 * tc[0]
            +6.83548940e-03 * tc[1]
            -8.19667665e-06 * tc[2]
            +4.03952522e-09 * tc[3]
            -8.33469780e-13 * tc[4];
        /*species 14: CO */
        species[14] =
            -1.43440860e+04 * invT
            -9.28875810e-01
            -3.57953347e+00 * tc[0]
            +3.05176840e-04 * tc[1]
            -1.69469055e-07 * tc[2]
            -7.55838237e-11 * tc[3]
            +4.52212249e-14 * tc[4];
        /*species 15: CO2 */
        species[15] =
            -4.83719697e+04 * invT
            -8.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
        /*species 16: HCO */
        species[16] =
            +3.83956496e+03 * invT
            -1.73186590e-01
            -4.22118584e+00 * tc[0]
            +1.62196266e-03 * tc[1]
            -2.29665743e-06 * tc[2]
            +1.10953411e-09 * tc[3]
            -2.16884432e-13 * tc[4];
        /*species 17: CH2O */
        species[17] =
            -1.43089567e+04 * invT
            +3.19091025e+00
            -4.79372315e+00 * tc[0]
            +4.95416684e-03 * tc[1]
            -6.22033347e-06 * tc[2]
            +3.16071051e-09 * tc[3]
            -6.58863260e-13 * tc[4];
        /*species 18: CH2OH */
        species[18] =
            -3.19391367e+03 * invT
            -2.60913325e+00
            -3.86388918e+00 * tc[0]
            -2.79836152e-03 * tc[1]
            -9.88786318e-07 * tc[2]
            +8.71100100e-10 * tc[3]
            -2.18483639e-13 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +9.78601100e+02 * invT
            -1.20459730e+01
            -2.10620400e+00 * tc[0]
            -3.60829750e-03 * tc[1]
            -8.89745333e-07 * tc[2]
            +6.14803000e-10 * tc[3]
            -1.03780500e-13 * tc[4];
        /*species 20: CH3OH */
        species[20] =
            -2.56427656e+04 * invT
            +6.21949405e+00
            -5.71539582e+00 * tc[0]
            +7.61545645e-03 * tc[1]
            -1.08740193e-05 * tc[2]
            +5.92339074e-09 * tc[3]
            -1.30676349e-12 * tc[4];
        /*species 21: C2H */
        species[21] =
            +6.68393932e+04 * invT
            -4.33330705e+00
            -2.88965733e+00 * tc[0]
            -6.70498055e-03 * tc[1]
            +4.74615835e-06 * tc[2]
            -2.45659204e-09 * tc[3]
            +5.46657555e-13 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            +2.64289807e+04 * invT
            -1.41310240e+01
            -8.08681094e-01 * tc[0]
            -1.16807815e-02 * tc[1]
            +5.91953025e-06 * tc[2]
            -2.33460364e-09 * tc[3]
            +4.25036487e-13 * tc[4];
        /*species 23: C2H3 */
        species[23] =
            +3.48598468e+04 * invT
            -6.29807380e+00
            -3.21246645e+00 * tc[0]
            -7.57395810e-04 * tc[1]
            -4.32015687e-06 * tc[2]
            +2.98048206e-09 * tc[3]
            -7.35754365e-13 * tc[4];
        /*species 24: C2H4 */
        species[24] =
            +5.08977593e+03 * invT
            -1.13812948e+00
            -3.95920148e+00 * tc[0]
            +3.78526124e-03 * tc[1]
            -9.51650487e-06 * tc[2]
            +5.76323961e-09 * tc[3]
            -1.34942187e-12 * tc[4];
        /*species 25: C2H5 */
        species[25] =
            +1.28416265e+04 * invT
            -1.40074356e+00
            -4.30646568e+00 * tc[0]
            +2.09329446e-03 * tc[1]
            -8.28571345e-06 * tc[2]
            +4.99272172e-09 * tc[3]
            -1.15254502e-12 * tc[4];
        /*species 26: C2H6 */
        species[26] =
            -1.15222055e+04 * invT
            +6.24601760e-01
            -4.29142492e+00 * tc[0]
            +2.75077135e-03 * tc[1]
            -9.99063813e-06 * tc[2]
            +5.90388571e-09 * tc[3]
            -1.34342886e-12 * tc[4];
        /*species 27: HCCO */
        species[27] =
            +2.00594490e+04 * invT
            -1.12386956e+01
            -2.25172140e+00 * tc[0]
            -8.82751050e-03 * tc[1]
            +3.95485017e-06 * tc[2]
            -1.43964658e-09 * tc[3]
            +2.53324055e-13 * tc[4];
        /*species 28: CH2CO */
        species[28] =
            -7.04291804e+03 * invT
            -1.10798117e+01
            -2.13583630e+00 * tc[0]
            -9.05943605e-03 * tc[1]
            +2.89912457e-06 * tc[2]
            -7.78664640e-10 * tc[3]
            +1.00728807e-13 * tc[4];
        /*species 29: HCCOH */
        species[29] =
            +8.03161430e+03 * invT
            -1.36319457e+01
            -1.24237330e+00 * tc[0]
            -1.55361005e-02 * tc[1]
            +8.47781067e-06 * tc[2]
            -3.59476092e-09 * tc[3]
            +7.00729700e-13 * tc[4];
        /*species 30: N2 */
        species[30] =
            -1.02089990e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
        /*species 31: AR */
        species[31] =
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
        /*species 7: H2O2 */
        species[7] =
            -1.78617877e+04 * invT
            +2.48846230e-01
            -4.16500285e+00 * tc[0]
            -2.45415847e-03 * tc[1]
            +3.16898708e-07 * tc[2]
            -3.09321655e-11 * tc[3]
            +1.43954153e-15 * tc[4];
        /*species 8: C */
        species[8] =
            +8.54512953e+04 * invT
            -3.30883485e+00
            -2.49266888e+00 * tc[0]
            -2.39944642e-05 * tc[1]
            +1.20722503e-08 * tc[2]
            -3.11909191e-12 * tc[3]
            +2.43638946e-16 * tc[4];
        /*species 9: CH */
        species[9] =
            +7.10124364e+04 * invT
            -3.60651526e+00
            -2.87846473e+00 * tc[0]
            -4.85456840e-04 * tc[1]
            -2.40742758e-08 * tc[2]
            +1.08906541e-11 * tc[3]
            -8.80396915e-16 * tc[4];
        /*species 10: CH2 */
        species[10] =
            +4.62636040e+04 * invT
            -4.29709211e+00
            -2.87410113e+00 * tc[0]
            -1.82819646e-03 * tc[1]
            +2.34824328e-07 * tc[2]
            -2.16816291e-11 * tc[3]
            +9.38637835e-16 * tc[4];
        /*species 11: CH2(S) */
        species[11] =
            +5.09259997e+04 * invT
            -7.33446327e+00
            -2.29203842e+00 * tc[0]
            -2.32794318e-03 * tc[1]
            +3.35319912e-07 * tc[2]
            -3.48255000e-11 * tc[3]
            +1.69858182e-15 * tc[4];
        /*species 12: CH3 */
        species[12] =
            +1.67755843e+04 * invT
            -7.19435407e+00
            -2.28571772e+00 * tc[0]
            -3.61995018e-03 * tc[1]
            +4.97857247e-07 * tc[2]
            -4.96403870e-11 * tc[3]
            +2.33577197e-15 * tc[4];
        /*species 13: CH4 */
        species[13] =
            -9.46834459e+03 * invT
            -1.93624665e+01
            -7.48514950e-02 * tc[0]
            -6.69547335e-03 * tc[1]
            +9.55476348e-07 * tc[2]
            -1.01910446e-10 * tc[3]
            +5.09076150e-15 * tc[4];
        /*species 14: CO */
        species[14] =
            -1.41518724e+04 * invT
            -6.10350211e+00
            -2.71518561e+00 * tc[0]
            -1.03126372e-03 * tc[1]
            +1.66470962e-07 * tc[2]
            -1.91710840e-11 * tc[3]
            +1.01823858e-15 * tc[4];
        /*species 15: CO2 */
        species[15] =
            -4.87591660e+04 * invT
            +5.85822230e-01
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
        /*species 16: HCO */
        species[16] =
            +4.01191815e+03 * invT
            -8.02617054e+00
            -2.77217438e+00 * tc[0]
            -2.47847763e-03 * tc[1]
            +4.14076022e-07 * tc[2]
            -4.90968148e-11 * tc[3]
            +2.66754356e-15 * tc[4];
        /*species 17: CH2O */
        species[17] =
            -1.39958323e+04 * invT
            -1.28956329e+01
            -1.76069008e+00 * tc[0]
            -4.60000041e-03 * tc[1]
            +7.37098022e-07 * tc[2]
            -8.38676767e-11 * tc[3]
            +4.41927820e-15 * tc[4];
        /*species 18: CH2OH */
        species[18] =
            -3.24250627e+03 * invT
            -3.11776646e+00
            -3.69266569e+00 * tc[0]
            -4.32288399e-03 * tc[1]
            +6.25168533e-07 * tc[2]
            -6.56028863e-11 * tc[3]
            +3.24277101e-15 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +1.27832520e+02 * invT
            -1.58776000e-01
            -3.77079900e+00 * tc[0]
            -3.93574850e-03 * tc[1]
            +4.42730667e-07 * tc[2]
            -3.28702583e-11 * tc[3]
            +1.05630800e-15 * tc[4];
        /*species 20: CH3OH */
        species[20] =
            -2.53748747e+04 * invT
            -1.37126544e+01
            -1.78970791e+00 * tc[0]
            -7.04691460e-03 * tc[1]
            +1.06083472e-06 * tc[2]
            -1.15142571e-10 * tc[3]
            +5.85301100e-15 * tc[4];
        /*species 21: C2H */
        species[21] =
            +6.71210650e+04 * invT
            -4.46808823e+00
            -3.16780652e+00 * tc[0]
            -2.37610951e-03 * tc[1]
            +3.06311795e-07 * tc[2]
            -2.53491877e-11 * tc[3]
            +8.86163850e-16 * tc[4];
        /*species 22: C2H2 */
        species[22] =
            +2.59359992e+04 * invT
            +4.37785085e+00
            -4.14756964e+00 * tc[0]
            -2.98083332e-03 * tc[1]
            +3.95491420e-07 * tc[2]
            -3.89510143e-11 * tc[3]
            +1.80617607e-15 * tc[4];
        /*species 23: C2H3 */
        species[23] =
            +3.46128739e+04 * invT
            -5.77059978e+00
            -3.01672400e+00 * tc[0]
            -5.16511460e-03 * tc[1]
            +7.80137248e-07 * tc[2]
            -8.48027400e-11 * tc[3]
            +4.31303520e-15 * tc[4];
        /*species 24: C2H4 */
        species[24] =
            +4.93988614e+03 * invT
            -9.26925814e+00
            -2.03611116e+00 * tc[0]
            -7.32270755e-03 * tc[1]
            +1.11846319e-06 * tc[2]
            -1.22685769e-10 * tc[3]
            +6.28530305e-15 * tc[4];
        /*species 25: C2H5 */
        species[25] =
            +1.28575200e+04 * invT
            -1.25077779e+01
            -1.95465642e+00 * tc[0]
            -8.69863610e-03 * tc[1]
            +1.33034445e-06 * tc[2]
            -1.46014741e-10 * tc[3]
            +7.48207880e-15 * tc[4];
        /*species 26: C2H6 */
        species[26] =
            -1.14263932e+04 * invT
            -1.50437292e+01
            -1.07188150e+00 * tc[0]
            -1.08426339e-02 * tc[1]
            +1.67093445e-06 * tc[2]
            -1.84510001e-10 * tc[3]
            +9.50014450e-15 * tc[4];
        /*species 27: HCCO */
        species[27] =
            +1.93272150e+04 * invT
            +8.55846530e+00
            -5.62820580e+00 * tc[0]
            -2.04267005e-03 * tc[1]
            +2.65575783e-07 * tc[2]
            -2.38550433e-11 * tc[3]
            +9.70391600e-16 * tc[4];
        /*species 28: CH2CO */
        species[28] =
            -7.55105311e+03 * invT
            +2.87905011e+00
            -4.51129732e+00 * tc[0]
            -4.50179872e-03 * tc[1]
            +6.94899392e-07 * tc[2]
            -7.69454902e-11 * tc[3]
            +3.97419100e-15 * tc[4];
        /*species 29: HCCOH */
        species[29] =
            +7.26462600e+03 * invT
            +1.25256033e+01
            -5.92382910e+00 * tc[0]
            -3.39618000e-03 * tc[1]
            +4.27642733e-07 * tc[2]
            -3.74898675e-11 * tc[3]
            +1.49700505e-15 * tc[4];
        /*species 30: N2 */
        species[30] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 31: AR */
        species[31] =
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
        /*species 7: H2O2 */
        species[7] =
            +3.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 * invT;
        /*species 8: C */
        species[8] =
            +1.55423955e+00
            -1.60768862e-04 * tc[1]
            +2.44597415e-07 * tc[2]
            -1.83058722e-10 * tc[3]
            +5.33042892e-14 * tc[4]
            +8.54438832e+04 * invT;
        /*species 9: CH */
        species[9] =
            +2.48981665e+00
            +1.61917771e-04 * tc[1]
            -5.62996883e-07 * tc[2]
            +7.90543317e-10 * tc[3]
            -2.81218134e-13 * tc[4]
            +7.07972934e+04 * invT;
        /*species 10: CH2 */
        species[10] =
            +2.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 * invT;
        /*species 11: CH2(S) */
        species[11] =
            +3.19860411e+00
            -1.18330710e-03 * tc[1]
            +2.74432073e-06 * tc[2]
            -1.67203995e-09 * tc[3]
            +3.88629474e-13 * tc[4]
            +5.04968163e+04 * invT;
        /*species 12: CH3 */
        species[12] =
            +2.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 * invT;
        /*species 13: CH4 */
        species[13] =
            +4.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 14: CO */
        species[14] =
            +2.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 15: CO2 */
        species[15] =
            +1.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 16: HCO */
        species[16] =
            +3.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 * invT;
        /*species 17: CH2O */
        species[17] =
            +3.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 * invT;
        /*species 18: CH2OH */
        species[18] =
            +2.86388918e+00
            +2.79836152e-03 * tc[1]
            +1.97757264e-06 * tc[2]
            -2.61330030e-09 * tc[3]
            +8.73934556e-13 * tc[4]
            -3.19391367e+03 * invT;
        /*species 19: CH3O */
        species[19] =
            +1.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122000e-13 * tc[4]
            +9.78601100e+02 * invT;
        /*species 20: CH3OH */
        species[20] =
            +4.71539582e+00
            -7.61545645e-03 * tc[1]
            +2.17480385e-05 * tc[2]
            -1.77701722e-08 * tc[3]
            +5.22705396e-12 * tc[4]
            -2.56427656e+04 * invT;
        /*species 21: C2H */
        species[21] =
            +1.88965733e+00
            +6.70498055e-03 * tc[1]
            -9.49231670e-06 * tc[2]
            +7.36977613e-09 * tc[3]
            -2.18663022e-12 * tc[4]
            +6.68393932e+04 * invT;
        /*species 22: C2H2 */
        species[22] =
            -1.91318906e-01
            +1.16807815e-02 * tc[1]
            -1.18390605e-05 * tc[2]
            +7.00381092e-09 * tc[3]
            -1.70014595e-12 * tc[4]
            +2.64289807e+04 * invT;
        /*species 23: C2H3 */
        species[23] =
            +2.21246645e+00
            +7.57395810e-04 * tc[1]
            +8.64031373e-06 * tc[2]
            -8.94144617e-09 * tc[3]
            +2.94301746e-12 * tc[4]
            +3.48598468e+04 * invT;
        /*species 24: C2H4 */
        species[24] =
            +2.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 * invT;
        /*species 25: C2H5 */
        species[25] =
            +3.30646568e+00
            -2.09329446e-03 * tc[1]
            +1.65714269e-05 * tc[2]
            -1.49781651e-08 * tc[3]
            +4.61018008e-12 * tc[4]
            +1.28416265e+04 * invT;
        /*species 26: C2H6 */
        species[26] =
            +3.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 * invT;
        /*species 27: HCCO */
        species[27] =
            +1.25172140e+00
            +8.82751050e-03 * tc[1]
            -7.90970033e-06 * tc[2]
            +4.31893975e-09 * tc[3]
            -1.01329622e-12 * tc[4]
            +2.00594490e+04 * invT;
        /*species 28: CH2CO */
        species[28] =
            +1.13583630e+00
            +9.05943605e-03 * tc[1]
            -5.79824913e-06 * tc[2]
            +2.33599392e-09 * tc[3]
            -4.02915230e-13 * tc[4]
            -7.04291804e+03 * invT;
        /*species 29: HCCOH */
        species[29] =
            +2.42373300e-01
            +1.55361005e-02 * tc[1]
            -1.69556213e-05 * tc[2]
            +1.07842828e-08 * tc[3]
            -2.80291880e-12 * tc[4]
            +8.03161430e+03 * invT;
        /*species 30: N2 */
        species[30] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
        /*species 31: AR */
        species[31] =
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
        /*species 7: H2O2 */
        species[7] =
            +3.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 * invT;
        /*species 8: C */
        species[8] =
            +1.49266888e+00
            +2.39944642e-05 * tc[1]
            -2.41445007e-08 * tc[2]
            +9.35727573e-12 * tc[3]
            -9.74555786e-16 * tc[4]
            +8.54512953e+04 * invT;
        /*species 9: CH */
        species[9] =
            +1.87846473e+00
            +4.85456840e-04 * tc[1]
            +4.81485517e-08 * tc[2]
            -3.26719623e-11 * tc[3]
            +3.52158766e-15 * tc[4]
            +7.10124364e+04 * invT;
        /*species 10: CH2 */
        species[10] =
            +1.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 * invT;
        /*species 11: CH2(S) */
        species[11] =
            +1.29203842e+00
            +2.32794318e-03 * tc[1]
            -6.70639823e-07 * tc[2]
            +1.04476500e-10 * tc[3]
            -6.79432730e-15 * tc[4]
            +5.09259997e+04 * invT;
        /*species 12: CH3 */
        species[12] =
            +1.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 * invT;
        /*species 13: CH4 */
        species[13] =
            -9.25148505e-01
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 14: CO */
        species[14] =
            +1.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 15: CO2 */
        species[15] =
            +2.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 16: HCO */
        species[16] =
            +1.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 * invT;
        /*species 17: CH2O */
        species[17] =
            +7.60690080e-01
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 * invT;
        /*species 18: CH2OH */
        species[18] =
            +2.69266569e+00
            +4.32288399e-03 * tc[1]
            -1.25033707e-06 * tc[2]
            +1.96808659e-10 * tc[3]
            -1.29710840e-14 * tc[4]
            -3.24250627e+03 * invT;
        /*species 19: CH3O */
        species[19] =
            +2.77079900e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832520e+02 * invT;
        /*species 20: CH3OH */
        species[20] =
            +7.89707910e-01
            +7.04691460e-03 * tc[1]
            -2.12166945e-06 * tc[2]
            +3.45427713e-10 * tc[3]
            -2.34120440e-14 * tc[4]
            -2.53748747e+04 * invT;
        /*species 21: C2H */
        species[21] =
            +2.16780652e+00
            +2.37610951e-03 * tc[1]
            -6.12623590e-07 * tc[2]
            +7.60475630e-11 * tc[3]
            -3.54465540e-15 * tc[4]
            +6.71210650e+04 * invT;
        /*species 22: C2H2 */
        species[22] =
            +3.14756964e+00
            +2.98083332e-03 * tc[1]
            -7.90982840e-07 * tc[2]
            +1.16853043e-10 * tc[3]
            -7.22470426e-15 * tc[4]
            +2.59359992e+04 * invT;
        /*species 23: C2H3 */
        species[23] =
            +2.01672400e+00
            +5.16511460e-03 * tc[1]
            -1.56027450e-06 * tc[2]
            +2.54408220e-10 * tc[3]
            -1.72521408e-14 * tc[4]
            +3.46128739e+04 * invT;
        /*species 24: C2H4 */
        species[24] =
            +1.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 * invT;
        /*species 25: C2H5 */
        species[25] =
            +9.54656420e-01
            +8.69863610e-03 * tc[1]
            -2.66068889e-06 * tc[2]
            +4.38044223e-10 * tc[3]
            -2.99283152e-14 * tc[4]
            +1.28575200e+04 * invT;
        /*species 26: C2H6 */
        species[26] =
            +7.18815000e-02
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 * invT;
        /*species 27: HCCO */
        species[27] =
            +4.62820580e+00
            +2.04267005e-03 * tc[1]
            -5.31151567e-07 * tc[2]
            +7.15651300e-11 * tc[3]
            -3.88156640e-15 * tc[4]
            +1.93272150e+04 * invT;
        /*species 28: CH2CO */
        species[28] =
            +3.51129732e+00
            +4.50179872e-03 * tc[1]
            -1.38979878e-06 * tc[2]
            +2.30836470e-10 * tc[3]
            -1.58967640e-14 * tc[4]
            -7.55105311e+03 * invT;
        /*species 29: HCCOH */
        species[29] =
            +4.92382910e+00
            +3.39618000e-03 * tc[1]
            -8.55285467e-07 * tc[2]
            +1.12469603e-10 * tc[3]
            -5.98802020e-15 * tc[4]
            +7.26462600e+03 * invT;
        /*species 30: N2 */
        species[30] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 31: AR */
        species[31] =
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
        /*species 7: H2O2 */
        species[7] =
            +4.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 * invT;
        /*species 8: C */
        species[8] =
            +2.55423955e+00
            -1.60768862e-04 * tc[1]
            +2.44597415e-07 * tc[2]
            -1.83058722e-10 * tc[3]
            +5.33042892e-14 * tc[4]
            +8.54438832e+04 * invT;
        /*species 9: CH */
        species[9] =
            +3.48981665e+00
            +1.61917771e-04 * tc[1]
            -5.62996883e-07 * tc[2]
            +7.90543317e-10 * tc[3]
            -2.81218134e-13 * tc[4]
            +7.07972934e+04 * invT;
        /*species 10: CH2 */
        species[10] =
            +3.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 * invT;
        /*species 11: CH2(S) */
        species[11] =
            +4.19860411e+00
            -1.18330710e-03 * tc[1]
            +2.74432073e-06 * tc[2]
            -1.67203995e-09 * tc[3]
            +3.88629474e-13 * tc[4]
            +5.04968163e+04 * invT;
        /*species 12: CH3 */
        species[12] =
            +3.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 * invT;
        /*species 13: CH4 */
        species[13] =
            +5.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 14: CO */
        species[14] =
            +3.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 15: CO2 */
        species[15] =
            +2.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 16: HCO */
        species[16] =
            +4.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 * invT;
        /*species 17: CH2O */
        species[17] =
            +4.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 * invT;
        /*species 18: CH2OH */
        species[18] =
            +3.86388918e+00
            +2.79836152e-03 * tc[1]
            +1.97757264e-06 * tc[2]
            -2.61330030e-09 * tc[3]
            +8.73934556e-13 * tc[4]
            -3.19391367e+03 * invT;
        /*species 19: CH3O */
        species[19] =
            +2.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122000e-13 * tc[4]
            +9.78601100e+02 * invT;
        /*species 20: CH3OH */
        species[20] =
            +5.71539582e+00
            -7.61545645e-03 * tc[1]
            +2.17480385e-05 * tc[2]
            -1.77701722e-08 * tc[3]
            +5.22705396e-12 * tc[4]
            -2.56427656e+04 * invT;
        /*species 21: C2H */
        species[21] =
            +2.88965733e+00
            +6.70498055e-03 * tc[1]
            -9.49231670e-06 * tc[2]
            +7.36977613e-09 * tc[3]
            -2.18663022e-12 * tc[4]
            +6.68393932e+04 * invT;
        /*species 22: C2H2 */
        species[22] =
            +8.08681094e-01
            +1.16807815e-02 * tc[1]
            -1.18390605e-05 * tc[2]
            +7.00381092e-09 * tc[3]
            -1.70014595e-12 * tc[4]
            +2.64289807e+04 * invT;
        /*species 23: C2H3 */
        species[23] =
            +3.21246645e+00
            +7.57395810e-04 * tc[1]
            +8.64031373e-06 * tc[2]
            -8.94144617e-09 * tc[3]
            +2.94301746e-12 * tc[4]
            +3.48598468e+04 * invT;
        /*species 24: C2H4 */
        species[24] =
            +3.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 * invT;
        /*species 25: C2H5 */
        species[25] =
            +4.30646568e+00
            -2.09329446e-03 * tc[1]
            +1.65714269e-05 * tc[2]
            -1.49781651e-08 * tc[3]
            +4.61018008e-12 * tc[4]
            +1.28416265e+04 * invT;
        /*species 26: C2H6 */
        species[26] =
            +4.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 * invT;
        /*species 27: HCCO */
        species[27] =
            +2.25172140e+00
            +8.82751050e-03 * tc[1]
            -7.90970033e-06 * tc[2]
            +4.31893975e-09 * tc[3]
            -1.01329622e-12 * tc[4]
            +2.00594490e+04 * invT;
        /*species 28: CH2CO */
        species[28] =
            +2.13583630e+00
            +9.05943605e-03 * tc[1]
            -5.79824913e-06 * tc[2]
            +2.33599392e-09 * tc[3]
            -4.02915230e-13 * tc[4]
            -7.04291804e+03 * invT;
        /*species 29: HCCOH */
        species[29] =
            +1.24237330e+00
            +1.55361005e-02 * tc[1]
            -1.69556213e-05 * tc[2]
            +1.07842828e-08 * tc[3]
            -2.80291880e-12 * tc[4]
            +8.03161430e+03 * invT;
        /*species 30: N2 */
        species[30] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
        /*species 31: AR */
        species[31] =
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
        /*species 7: H2O2 */
        species[7] =
            +4.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 * invT;
        /*species 8: C */
        species[8] =
            +2.49266888e+00
            +2.39944642e-05 * tc[1]
            -2.41445007e-08 * tc[2]
            +9.35727573e-12 * tc[3]
            -9.74555786e-16 * tc[4]
            +8.54512953e+04 * invT;
        /*species 9: CH */
        species[9] =
            +2.87846473e+00
            +4.85456840e-04 * tc[1]
            +4.81485517e-08 * tc[2]
            -3.26719623e-11 * tc[3]
            +3.52158766e-15 * tc[4]
            +7.10124364e+04 * invT;
        /*species 10: CH2 */
        species[10] =
            +2.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 * invT;
        /*species 11: CH2(S) */
        species[11] =
            +2.29203842e+00
            +2.32794318e-03 * tc[1]
            -6.70639823e-07 * tc[2]
            +1.04476500e-10 * tc[3]
            -6.79432730e-15 * tc[4]
            +5.09259997e+04 * invT;
        /*species 12: CH3 */
        species[12] =
            +2.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 * invT;
        /*species 13: CH4 */
        species[13] =
            +7.48514950e-02
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 14: CO */
        species[14] =
            +2.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 15: CO2 */
        species[15] =
            +3.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 16: HCO */
        species[16] =
            +2.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 * invT;
        /*species 17: CH2O */
        species[17] =
            +1.76069008e+00
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 * invT;
        /*species 18: CH2OH */
        species[18] =
            +3.69266569e+00
            +4.32288399e-03 * tc[1]
            -1.25033707e-06 * tc[2]
            +1.96808659e-10 * tc[3]
            -1.29710840e-14 * tc[4]
            -3.24250627e+03 * invT;
        /*species 19: CH3O */
        species[19] =
            +3.77079900e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832520e+02 * invT;
        /*species 20: CH3OH */
        species[20] =
            +1.78970791e+00
            +7.04691460e-03 * tc[1]
            -2.12166945e-06 * tc[2]
            +3.45427713e-10 * tc[3]
            -2.34120440e-14 * tc[4]
            -2.53748747e+04 * invT;
        /*species 21: C2H */
        species[21] =
            +3.16780652e+00
            +2.37610951e-03 * tc[1]
            -6.12623590e-07 * tc[2]
            +7.60475630e-11 * tc[3]
            -3.54465540e-15 * tc[4]
            +6.71210650e+04 * invT;
        /*species 22: C2H2 */
        species[22] =
            +4.14756964e+00
            +2.98083332e-03 * tc[1]
            -7.90982840e-07 * tc[2]
            +1.16853043e-10 * tc[3]
            -7.22470426e-15 * tc[4]
            +2.59359992e+04 * invT;
        /*species 23: C2H3 */
        species[23] =
            +3.01672400e+00
            +5.16511460e-03 * tc[1]
            -1.56027450e-06 * tc[2]
            +2.54408220e-10 * tc[3]
            -1.72521408e-14 * tc[4]
            +3.46128739e+04 * invT;
        /*species 24: C2H4 */
        species[24] =
            +2.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 * invT;
        /*species 25: C2H5 */
        species[25] =
            +1.95465642e+00
            +8.69863610e-03 * tc[1]
            -2.66068889e-06 * tc[2]
            +4.38044223e-10 * tc[3]
            -2.99283152e-14 * tc[4]
            +1.28575200e+04 * invT;
        /*species 26: C2H6 */
        species[26] =
            +1.07188150e+00
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 * invT;
        /*species 27: HCCO */
        species[27] =
            +5.62820580e+00
            +2.04267005e-03 * tc[1]
            -5.31151567e-07 * tc[2]
            +7.15651300e-11 * tc[3]
            -3.88156640e-15 * tc[4]
            +1.93272150e+04 * invT;
        /*species 28: CH2CO */
        species[28] =
            +4.51129732e+00
            +4.50179872e-03 * tc[1]
            -1.38979878e-06 * tc[2]
            +2.30836470e-10 * tc[3]
            -1.58967640e-14 * tc[4]
            -7.55105311e+03 * invT;
        /*species 29: HCCOH */
        species[29] =
            +5.92382910e+00
            +3.39618000e-03 * tc[1]
            -8.55285467e-07 * tc[2]
            +1.12469603e-10 * tc[3]
            -5.98802020e-15 * tc[4]
            +7.26462600e+03 * invT;
        /*species 30: N2 */
        species[30] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 31: AR */
        species[31] =
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
    wt[8] = 12.011150; /*C */
    wt[9] = 13.019120; /*CH */
    wt[10] = 14.027090; /*CH2 */
    wt[11] = 14.027090; /*CH2(S) */
    wt[12] = 15.035060; /*CH3 */
    wt[13] = 16.043030; /*CH4 */
    wt[14] = 28.010550; /*CO */
    wt[15] = 44.009950; /*CO2 */
    wt[16] = 29.018520; /*HCO */
    wt[17] = 30.026490; /*CH2O */
    wt[18] = 31.034460; /*CH2OH */
    wt[19] = 31.034460; /*CH3O */
    wt[20] = 32.042430; /*CH3OH */
    wt[21] = 25.030270; /*C2H */
    wt[22] = 26.038240; /*C2H2 */
    wt[23] = 27.046210; /*C2H3 */
    wt[24] = 28.054180; /*C2H4 */
    wt[25] = 29.062150; /*C2H5 */
    wt[26] = 30.070120; /*C2H6 */
    wt[27] = 41.029670; /*HCCO */
    wt[28] = 42.037640; /*CH2CO */
    wt[29] = 42.037640; /*HCCOH */
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
    y[8] = phi[8]*12.011150;   XW += y[8]; /*C */
    y[9] = phi[9]*13.019120;   XW += y[9]; /*CH */
    y[10] = phi[10]*14.027090;   XW += y[10]; /*CH2 */
    y[11] = phi[11]*14.027090;   XW += y[11]; /*CH2(S) */
    y[12] = phi[12]*15.035060;   XW += y[12]; /*CH3 */
    y[13] = phi[13]*16.043030;   XW += y[13]; /*CH4 */
    y[14] = phi[14]*28.010550;   XW += y[14]; /*CO */
    y[15] = phi[15]*44.009950;   XW += y[15]; /*CO2 */
    y[16] = phi[16]*29.018520;   XW += y[16]; /*HCO */
    y[17] = phi[17]*30.026490;   XW += y[17]; /*CH2O */
    y[18] = phi[18]*31.034460;   XW += y[18]; /*CH2OH */
    y[19] = phi[19]*31.034460;   XW += y[19]; /*CH3O */
    y[20] = phi[20]*32.042430;   XW += y[20]; /*CH3OH */
    y[21] = phi[21]*25.030270;   XW += y[21]; /*C2H */
    y[22] = phi[22]*26.038240;   XW += y[22]; /*C2H2 */
    y[23] = phi[23]*27.046210;   XW += y[23]; /*C2H3 */
    y[24] = phi[24]*28.054180;   XW += y[24]; /*C2H4 */
    y[25] = phi[25]*29.062150;   XW += y[25]; /*C2H5 */
    y[26] = phi[26]*30.070120;   XW += y[26]; /*C2H6 */
    y[27] = phi[27]*41.029670;   XW += y[27]; /*HCCO */
    y[28] = phi[28]*42.037640;   XW += y[28]; /*CH2CO */
    y[29] = phi[29]*42.037640;   XW += y[29]; /*HCCOH */
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
    phi[8] = y[8]/ 1.20111500e-02; /*C (wt in kg) */
    phi[9] = y[9]/ 1.30191200e-02; /*CH (wt in kg) */
    phi[10] = y[10]/ 1.40270900e-02; /*CH2 (wt in kg) */
    phi[11] = y[11]/ 1.40270900e-02; /*CH2(S) (wt in kg) */
    phi[12] = y[12]/ 1.50350600e-02; /*CH3 (wt in kg) */
    phi[13] = y[13]/ 1.60430300e-02; /*CH4 (wt in kg) */
    phi[14] = y[14]/ 2.80105500e-02; /*CO (wt in kg) */
    phi[15] = y[15]/ 4.40099500e-02; /*CO2 (wt in kg) */
    phi[16] = y[16]/ 2.90185200e-02; /*HCO (wt in kg) */
    phi[17] = y[17]/ 3.00264900e-02; /*CH2O (wt in kg) */
    phi[18] = y[18]/ 3.10344600e-02; /*CH2OH (wt in kg) */
    phi[19] = y[19]/ 3.10344600e-02; /*CH3O (wt in kg) */
    phi[20] = y[20]/ 3.20424300e-02; /*CH3OH (wt in kg) */
    phi[21] = y[21]/ 2.50302700e-02; /*C2H (wt in kg) */
    phi[22] = y[22]/ 2.60382400e-02; /*C2H2 (wt in kg) */
    phi[23] = y[23]/ 2.70462100e-02; /*C2H3 (wt in kg) */
    phi[24] = y[24]/ 2.80541800e-02; /*C2H4 (wt in kg) */
    phi[25] = y[25]/ 2.90621500e-02; /*C2H5 (wt in kg) */
    phi[26] = y[26]/ 3.00701200e-02; /*C2H6 (wt in kg) */
    phi[27] = y[27]/ 4.10296700e-02; /*HCCO (wt in kg) */
    phi[28] = y[28]/ 4.20376400e-02; /*CH2CO (wt in kg) */
    phi[29] = y[29]/ 4.20376400e-02; /*HCCOH (wt in kg) */
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
    y[8] = c[8] * 12.011150 / (*rho); 
    y[9] = c[9] * 13.019120 / (*rho); 
    y[10] = c[10] * 14.027090 / (*rho); 
    y[11] = c[11] * 14.027090 / (*rho); 
    y[12] = c[12] * 15.035060 / (*rho); 
    y[13] = c[13] * 16.043030 / (*rho); 
    y[14] = c[14] * 28.010550 / (*rho); 
    y[15] = c[15] * 44.009950 / (*rho); 
    y[16] = c[16] * 29.018520 / (*rho); 
    y[17] = c[17] * 30.026490 / (*rho); 
    y[18] = c[18] * 31.034460 / (*rho); 
    y[19] = c[19] * 31.034460 / (*rho); 
    y[20] = c[20] * 32.042430 / (*rho); 
    y[21] = c[21] * 25.030270 / (*rho); 
    y[22] = c[22] * 26.038240 / (*rho); 
    y[23] = c[23] * 27.046210 / (*rho); 
    y[24] = c[24] * 28.054180 / (*rho); 
    y[25] = c[25] * 29.062150 / (*rho); 
    y[26] = c[26] * 30.070120 / (*rho); 
    y[27] = c[27] * 41.029670 / (*rho); 
    y[28] = c[28] * 42.037640 / (*rho); 
    y[29] = c[29] * 42.037640 / (*rho); 
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
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENIMC EGTRANSETLENIMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENIMC egtransetlenimc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENIMC egtransetlenimc_
#endif
extern "C" { void egtransetLENIMC(int* LENIMC) {
  *LENIMC =          130;}}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENRMC EGTRANSETLENRMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENRMC egtransetlenrmc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENRMC egtransetlenrmc_
#endif
extern "C" { void egtransetLENRMC(int* LENRMC) {
  *LENRMC =        20576;}}
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
  *KK =           32;}}
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
  WT[           7] =   0.3401474022865295E+02;
  WT[           8] =   0.1201115036010742E+02;
  WT[           9] =   0.1301912033557892E+02;
  WT[          10] =   0.1402709031105042E+02;
  WT[          11] =   0.1402709031105042E+02;
  WT[          12] =   0.1503506028652191E+02;
  WT[          13] =   0.1604303026199341E+02;
  WT[          14] =   0.2801055049896240E+02;
  WT[          15] =   0.4400995063781738E+02;
  WT[          16] =   0.2901852047443390E+02;
  WT[          17] =   0.3002649044990540E+02;
  WT[          18] =   0.3103446042537689E+02;
  WT[          19] =   0.3103446042537689E+02;
  WT[          20] =   0.3204243040084839E+02;
  WT[          21] =   0.2503027069568634E+02;
  WT[          22] =   0.2603824067115784E+02;
  WT[          23] =   0.2704621064662933E+02;
  WT[          24] =   0.2805418062210083E+02;
  WT[          25] =   0.2906215059757233E+02;
  WT[          26] =   0.3007012057304382E+02;
  WT[          27] =   0.4102967083454132E+02;
  WT[          28] =   0.4203764081001282E+02;
  WT[          29] =   0.4203764081001282E+02;
  WT[          30] =   0.2801339912414551E+02;
  WT[          31] =   0.3994800186157227E+02;
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
  EPS[           7] =   0.1074000000000000E+03;
  EPS[           8] =   0.7140000000000001E+02;
  EPS[           9] =   0.8000000000000000E+02;
  EPS[          10] =   0.1440000000000000E+03;
  EPS[          11] =   0.1440000000000000E+03;
  EPS[          12] =   0.1440000000000000E+03;
  EPS[          13] =   0.1414000000000000E+03;
  EPS[          14] =   0.9809999999999999E+02;
  EPS[          15] =   0.2440000000000000E+03;
  EPS[          16] =   0.4980000000000000E+03;
  EPS[          17] =   0.4980000000000000E+03;
  EPS[          18] =   0.4170000000000000E+03;
  EPS[          19] =   0.4170000000000000E+03;
  EPS[          20] =   0.4818000000000000E+03;
  EPS[          21] =   0.2090000000000000E+03;
  EPS[          22] =   0.2090000000000000E+03;
  EPS[          23] =   0.2090000000000000E+03;
  EPS[          24] =   0.2808000000000000E+03;
  EPS[          25] =   0.2523000000000000E+03;
  EPS[          26] =   0.2523000000000000E+03;
  EPS[          27] =   0.1500000000000000E+03;
  EPS[          28] =   0.4360000000000000E+03;
  EPS[          29] =   0.4360000000000000E+03;
  EPS[          30] =   0.9753000000000000E+02;
  EPS[          31] =   0.1365000000000000E+03;
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
  SIG[           7] =   0.3458000000000000E+01;
  SIG[           8] =   0.3298000000000000E+01;
  SIG[           9] =   0.2750000000000000E+01;
  SIG[          10] =   0.3800000000000000E+01;
  SIG[          11] =   0.3800000000000000E+01;
  SIG[          12] =   0.3800000000000000E+01;
  SIG[          13] =   0.3746000000000000E+01;
  SIG[          14] =   0.3650000000000000E+01;
  SIG[          15] =   0.3763000000000000E+01;
  SIG[          16] =   0.3590000000000000E+01;
  SIG[          17] =   0.3590000000000000E+01;
  SIG[          18] =   0.3690000000000000E+01;
  SIG[          19] =   0.3690000000000000E+01;
  SIG[          20] =   0.3626000000000000E+01;
  SIG[          21] =   0.4100000000000000E+01;
  SIG[          22] =   0.4100000000000000E+01;
  SIG[          23] =   0.4100000000000000E+01;
  SIG[          24] =   0.3971000000000000E+01;
  SIG[          25] =   0.4302000000000000E+01;
  SIG[          26] =   0.4302000000000000E+01;
  SIG[          27] =   0.2500000000000000E+01;
  SIG[          28] =   0.3970000000000000E+01;
  SIG[          29] =   0.3970000000000000E+01;
  SIG[          30] =   0.3621000000000000E+01;
  SIG[          31] =   0.3330000000000000E+01;
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
  DIP[          15] =   0.0000000000000000E+00;
  DIP[          16] =   0.0000000000000000E+00;
  DIP[          17] =   0.0000000000000000E+00;
  DIP[          18] =   0.1700000000000000E+01;
  DIP[          19] =   0.1700000000000000E+01;
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
  DIP[          30] =   0.0000000000000000E+00;
  DIP[          31] =   0.0000000000000000E+00;
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
  POL[          10] =   0.0000000000000000E+00;
  POL[          11] =   0.0000000000000000E+00;
  POL[          12] =   0.0000000000000000E+00;
  POL[          13] =   0.2600000000000000E+01;
  POL[          14] =   0.1950000000000000E+01;
  POL[          15] =   0.2650000000000000E+01;
  POL[          16] =   0.0000000000000000E+00;
  POL[          17] =   0.0000000000000000E+00;
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
  POL[          29] =   0.0000000000000000E+00;
  POL[          30] =   0.1760000000000000E+01;
  POL[          31] =   0.0000000000000000E+00;
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
  ZROT[           7] =   0.3800000000000000E+01;
  ZROT[           8] =   0.0000000000000000E+00;
  ZROT[           9] =   0.0000000000000000E+00;
  ZROT[          10] =   0.0000000000000000E+00;
  ZROT[          11] =   0.0000000000000000E+00;
  ZROT[          12] =   0.0000000000000000E+00;
  ZROT[          13] =   0.1300000000000000E+02;
  ZROT[          14] =   0.1800000000000000E+01;
  ZROT[          15] =   0.2100000000000000E+01;
  ZROT[          16] =   0.0000000000000000E+00;
  ZROT[          17] =   0.2000000000000000E+01;
  ZROT[          18] =   0.2000000000000000E+01;
  ZROT[          19] =   0.2000000000000000E+01;
  ZROT[          20] =   0.1000000000000000E+01;
  ZROT[          21] =   0.2500000000000000E+01;
  ZROT[          22] =   0.2500000000000000E+01;
  ZROT[          23] =   0.1000000000000000E+01;
  ZROT[          24] =   0.1500000000000000E+01;
  ZROT[          25] =   0.1500000000000000E+01;
  ZROT[          26] =   0.1500000000000000E+01;
  ZROT[          27] =   0.1000000000000000E+01;
  ZROT[          28] =   0.2000000000000000E+01;
  ZROT[          29] =   0.2000000000000000E+01;
  ZROT[          30] =   0.4000000000000000E+01;
  ZROT[          31] =   0.0000000000000000E+00;
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
  NLIN[           7] =            2;
  NLIN[           8] =            0;
  NLIN[           9] =            1;
  NLIN[          10] =            1;
  NLIN[          11] =            1;
  NLIN[          12] =            1;
  NLIN[          13] =            2;
  NLIN[          14] =            1;
  NLIN[          15] =            1;
  NLIN[          16] =            2;
  NLIN[          17] =            2;
  NLIN[          18] =            2;
  NLIN[          19] =            2;
  NLIN[          20] =            2;
  NLIN[          21] =            1;
  NLIN[          22] =            1;
  NLIN[          23] =            2;
  NLIN[          24] =            2;
  NLIN[          25] =            2;
  NLIN[          26] =            2;
  NLIN[          27] =            2;
  NLIN[          28] =            2;
  NLIN[          29] =            2;
  NLIN[          30] =            1;
  NLIN[          31] =            0;
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
  COFLAM[          28] =   0.8827769530599365E+00;
  COFLAM[          29] =   0.1315528335457282E+01;
  COFLAM[          30] =   0.1916184484138214E-01;
  COFLAM[          31] =  -0.4416817199461963E-02;
  COFLAM[          32] =   0.2041557184183207E+01;
  COFLAM[          33] =   0.1706343688810633E+01;
  COFLAM[          34] =  -0.1443632622298301E+00;
  COFLAM[          35] =   0.6539115411757299E-02;
  COFLAM[          36] =   0.2080465726344590E+02;
  COFLAM[          37] =  -0.6239956140710875E+01;
  COFLAM[          38] =   0.9821413259074330E+00;
  COFLAM[          39] =  -0.4502534139719348E-01;
  COFLAM[          40] =   0.1291635112732534E+02;
  COFLAM[          41] =  -0.3737065816673101E+01;
  COFLAM[          42] =   0.7157937477351572E+00;
  COFLAM[          43] =  -0.3638374043162176E-01;
  COFLAM[          44] =   0.1893642052511076E+02;
  COFLAM[          45] =  -0.6509584574818159E+01;
  COFLAM[          46] =   0.1132853049549838E+01;
  COFLAM[          47] =  -0.5695791215151112E-01;
  COFLAM[          48] =   0.1399177056183472E+02;
  COFLAM[          49] =  -0.4641882989478738E+01;
  COFLAM[          50] =   0.9076435786312527E+00;
  COFLAM[          51] =  -0.4772397826828819E-01;
  COFLAM[          52] =   0.1330618431273573E+02;
  COFLAM[          53] =  -0.4960294456747341E+01;
  COFLAM[          54] =   0.1032808842591664E+01;
  COFLAM[          55] =  -0.5633567903356620E-01;
  COFLAM[          56] =   0.1187710072654253E+02;
  COFLAM[          57] =  -0.3154801252509626E+01;
  COFLAM[          58] =   0.6020483454882906E+00;
  COFLAM[          59] =  -0.3032714732778013E-01;
  COFLAM[          60] =  -0.1135070822782991E+02;
  COFLAM[          61] =   0.5875667873528579E+01;
  COFLAM[          62] =  -0.5677982250304914E+00;
  COFLAM[          63] =   0.2031670238952729E-01;
  COFLAM[          64] =   0.6296028277200143E+01;
  COFLAM[          65] =  -0.2225281214264024E+01;
  COFLAM[          66] =   0.6369149369376514E+00;
  COFLAM[          67] =  -0.3808451133860326E-01;
  COFLAM[          68] =   0.5384151635429418E+01;
  COFLAM[          69] =  -0.2389146833355944E+01;
  COFLAM[          70] =   0.7389861798813848E+00;
  COFLAM[          71] =  -0.4581386404459022E-01;
  COFLAM[          72] =   0.1061523737201134E+01;
  COFLAM[          73] =   0.4055448816410179E-01;
  COFLAM[          74] =   0.3317451351857446E+00;
  COFLAM[          75] =  -0.2409090777140087E-01;
  COFLAM[          76] =  -0.6138056391436082E+01;
  COFLAM[          77] =   0.2471262972505274E+01;
  COFLAM[          78] =   0.6476625884336293E-01;
  COFLAM[          79] =  -0.1455104151494988E-01;
  COFLAM[          80] =  -0.1838961899802688E+01;
  COFLAM[          81] =   0.5576367774184181E+00;
  COFLAM[          82] =   0.3510283465541862E+00;
  COFLAM[          83] =  -0.2853844971316656E-01;
  COFLAM[          84] =   0.1336229863097232E+02;
  COFLAM[          85] =  -0.4326103989562728E+01;
  COFLAM[          86] =   0.8367814315037552E+00;
  COFLAM[          87] =  -0.4380214396494905E-01;
  COFLAM[          88] =  -0.7691048523709211E+01;
  COFLAM[          89] =   0.4564166690336450E+01;
  COFLAM[          90] =  -0.4040787948018396E+00;
  COFLAM[          91] =   0.1405248077558339E-01;
  COFLAM[          92] =  -0.9095387619073980E+01;
  COFLAM[          93] =   0.4544249531780739E+01;
  COFLAM[          94] =  -0.3175822855243875E+00;
  COFLAM[          95] =   0.6570713165664151E-02;
  COFLAM[          96] =  -0.1460870939008838E+02;
  COFLAM[          97] =   0.6359722173433303E+01;
  COFLAM[          98] =  -0.5034536870651064E+00;
  COFLAM[          99] =   0.1259513381293373E-01;
  COFLAM[         100] =  -0.8941827721234478E+01;
  COFLAM[         101] =   0.4021583693700529E+01;
  COFLAM[         102] =  -0.1835686948063458E+00;
  COFLAM[         103] =  -0.1963330023789806E-02;
  COFLAM[         104] =  -0.1098227802311090E+02;
  COFLAM[         105] =   0.4703047406066884E+01;
  COFLAM[         106] =  -0.2517963879729279E+00;
  COFLAM[         107] =   0.1532853967214294E-03;
  COFLAM[         108] =  -0.5965774793396434E+01;
  COFLAM[         109] =   0.4394152739209455E+01;
  COFLAM[         110] =  -0.3976634471827024E+00;
  COFLAM[         111] =   0.1397721984656401E-01;
  COFLAM[         112] =  -0.8329425744536485E+01;
  COFLAM[         113] =   0.3971067920427290E+01;
  COFLAM[         114] =  -0.2213342783035633E+00;
  COFLAM[         115] =   0.1483723034163886E-02;
  COFLAM[         116] =  -0.5367731950169370E+01;
  COFLAM[         117] =   0.2983564111219091E+01;
  COFLAM[         118] =  -0.1147307150440449E+00;
  COFLAM[         119] =  -0.2123138953210028E-02;
  COFLAM[         120] =   0.1293004274651541E+02;
  COFLAM[         121] =  -0.3528374680486067E+01;
  COFLAM[         122] =   0.6455829015382131E+00;
  COFLAM[         123] =  -0.3194413600157287E-01;
  COFLAM[         124] =  -0.3166636127434923E+01;
  COFLAM[         125] =   0.3467381630049378E+01;
  COFLAM[         126] =  -0.3746257298439007E+00;
  COFLAM[         127] =   0.1658331946651994E-01;
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
  COFETA[          28] =  -0.1712754275668106E+02;
  COFETA[          29] =   0.2678088349030581E+01;
  COFETA[          30] =  -0.2721592407921878E+00;
  COFETA[          31] =   0.1214173232611457E-01;
  COFETA[          32] =  -0.1503046107299309E+02;
  COFETA[          33] =   0.1706343688810645E+01;
  COFETA[          34] =  -0.1443632622298313E+00;
  COFETA[          35] =   0.6539115411757332E-02;
  COFETA[          36] =  -0.1520334313770844E+02;
  COFETA[          37] =   0.1929024677614508E+01;
  COFETA[          38] =  -0.1738657445088065E+00;
  COFETA[          39] =   0.7841476915125538E-02;
  COFETA[          40] =  -0.2026303235679370E+02;
  COFETA[          41] =   0.3630408718312100E+01;
  COFETA[          42] =  -0.3952256976799819E+00;
  COFETA[          43] =   0.1745288966941290E-01;
  COFETA[          44] =  -0.2026303235679370E+02;
  COFETA[          45] =   0.3630408718312100E+01;
  COFETA[          46] =  -0.3952256976799819E+00;
  COFETA[          47] =   0.1745288966941290E-01;
  COFETA[          48] =  -0.2022833518474935E+02;
  COFETA[          49] =   0.3630408718312077E+01;
  COFETA[          50] =  -0.3952256976799788E+00;
  COFETA[          51] =   0.1745288966941275E-01;
  COFETA[          52] =  -0.2000457400541702E+02;
  COFETA[          53] =   0.3569542093003574E+01;
  COFETA[          54] =  -0.3874920392785122E+00;
  COFETA[          55] =   0.1712461411247914E-01;
  COFETA[          56] =  -0.1661561262594180E+02;
  COFETA[          57] =   0.2400975158113569E+01;
  COFETA[          58] =  -0.2357717790312281E+00;
  COFETA[          59] =   0.1054820948438182E-01;
  COFETA[          60] =  -0.2397057295682966E+02;
  COFETA[          61] =   0.5130426196036950E+01;
  COFETA[          62] =  -0.5724284704186094E+00;
  COFETA[          63] =   0.2440888721969576E-01;
  COFETA[          64] =  -0.1987003264957855E+02;
  COFETA[          65] =   0.2703811616927088E+01;
  COFETA[          66] =  -0.1672355189642018E+00;
  COFETA[          67] =   0.3212257118450056E-02;
  COFETA[          68] =  -0.1985295977766010E+02;
  COFETA[          69] =   0.2703811616927011E+01;
  COFETA[          70] =  -0.1672355189641908E+00;
  COFETA[          71] =   0.3212257118449540E-02;
  COFETA[          72] =  -0.1997384799268116E+02;
  COFETA[          73] =   0.2861750522819142E+01;
  COFETA[          74] =  -0.2024899143538781E+00;
  COFETA[          75] =   0.5362464436156140E-02;
  COFETA[          76] =  -0.1997384799268116E+02;
  COFETA[          77] =   0.2861750522819142E+01;
  COFETA[          78] =  -0.2024899143538781E+00;
  COFETA[          79] =   0.5362464436156140E-02;
  COFETA[          80] =  -0.2057883672241691E+02;
  COFETA[          81] =   0.3046121203848835E+01;
  COFETA[          82] =  -0.2180172142434544E+00;
  COFETA[          83] =   0.5665599148971810E-02;
  COFETA[          84] =  -0.2335627494271669E+02;
  COFETA[          85] =   0.4790351551626667E+01;
  COFETA[          86] =  -0.5364560275637629E+00;
  COFETA[          87] =   0.2318560947296586E-01;
  COFETA[          88] =  -0.2333653477658102E+02;
  COFETA[          89] =   0.4790351551626705E+01;
  COFETA[          90] =  -0.5364560275637680E+00;
  COFETA[          91] =   0.2318560947296609E-01;
  COFETA[          92] =  -0.2331754444859468E+02;
  COFETA[          93] =   0.4790351551626613E+01;
  COFETA[          94] =  -0.5364560275637551E+00;
  COFETA[          95] =   0.2318560947296548E-01;
  COFETA[          96] =  -0.2504838169296130E+02;
  COFETA[          97] =   0.5332796564390715E+01;
  COFETA[          98] =  -0.5890465746677508E+00;
  COFETA[          99] =   0.2473730756384840E-01;
  COFETA[         100] =  -0.2463014790579715E+02;
  COFETA[         101] =   0.5183142307718661E+01;
  COFETA[         102] =  -0.5771878965089483E+00;
  COFETA[         103] =   0.2453023678847775E-01;
  COFETA[         104] =  -0.2461310023284140E+02;
  COFETA[         105] =   0.5183142307718687E+01;
  COFETA[         106] =  -0.5771878965089522E+00;
  COFETA[         107] =   0.2453023678847794E-01;
  COFETA[         108] =  -0.1921922328200846E+02;
  COFETA[         109] =   0.3751322790269299E+01;
  COFETA[         110] =  -0.4102632112528644E+00;
  COFETA[         111] =   0.1807670334472438E-01;
  COFETA[         112] =  -0.2234339879666697E+02;
  COFETA[         113] =   0.3866553929075889E+01;
  COFETA[         114] =  -0.3419581218331971E+00;
  COFETA[         115] =   0.1173156132420011E-01;
  COFETA[         116] =  -0.2234339879666697E+02;
  COFETA[         117] =   0.3866553929075889E+01;
  COFETA[         118] =  -0.3419581218331971E+00;
  COFETA[         119] =   0.1173156132420011E-01;
  COFETA[         120] =  -0.1656563666406150E+02;
  COFETA[         121] =   0.2388167035581858E+01;
  COFETA[         122] =  -0.2341208182867099E+00;
  COFETA[         123] =   0.1047727172770704E-01;
  COFETA[         124] =  -0.1903691114465788E+02;
  COFETA[         125] =   0.3467381630049308E+01;
  COFETA[         126] =  -0.3746257298438905E+00;
  COFETA[         127] =   0.1658331946651944E-01;
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
  COFD[          28] =  -0.1271518733681371E+02;
  COFD[          29] =   0.2932402241140174E+01;
  COFD[          30] =  -0.1743585082612658E+00;
  COFD[          31] =   0.7919486191765343E-02;
  COFD[          32] =  -0.1081045105273500E+02;
  COFD[          33] =   0.2188985884544625E+01;
  COFD[          34] =  -0.7096785665772207E-01;
  COFD[          35] =   0.3131621700706094E-02;
  COFD[          36] =  -0.1092230647098016E+02;
  COFD[          37] =   0.2306200733264308E+01;
  COFD[          38] =  -0.8737245270736949E-01;
  COFD[          39] =   0.3897864029739334E-02;
  COFD[          40] =  -0.1355980676188247E+02;
  COFD[          41] =   0.3231182295054575E+01;
  COFD[          42] =  -0.2138406992439825E+00;
  COFD[          43] =   0.9659145053095021E-02;
  COFD[          44] =  -0.1355980676188247E+02;
  COFD[          45] =   0.3231182295054575E+01;
  COFD[          46] =  -0.2138406992439825E+00;
  COFD[          47] =   0.9659145053095021E-02;
  COFD[          48] =  -0.1358415095584097E+02;
  COFD[          49] =   0.3239845970896248E+01;
  COFD[          50] =  -0.2150500406377150E+00;
  COFD[          51] =   0.9715225271588263E-02;
  COFD[          52] =  -0.1353824314289343E+02;
  COFD[          53] =   0.3227352939332374E+01;
  COFD[          54] =  -0.2134245159235849E+00;
  COFD[          55] =   0.9644648995332943E-02;
  COFD[          56] =  -0.1268296956267843E+02;
  COFD[          57] =   0.2903428401599904E+01;
  COFD[          58] =  -0.1707667302851327E+00;
  COFD[          59] =   0.7771111802295727E-02;
  COFD[          60] =  -0.1564739122681956E+02;
  COFD[          61] =   0.4025825276486324E+01;
  COFD[          62] =  -0.3181195147607377E+00;
  COFD[          63] =   0.1422136468533644E-01;
  COFD[          64] =  -0.1746999275651299E+02;
  COFD[          65] =   0.4670556529920976E+01;
  COFD[          66] =  -0.3931788949833612E+00;
  COFD[          67] =   0.1709979614643272E-01;
  COFD[          68] =  -0.1747771970561251E+02;
  COFD[          69] =   0.4673159679696371E+01;
  COFD[          70] =  -0.3935109715079382E+00;
  COFD[          71] =   0.1711386921544412E-01;
  COFD[          72] =  -0.1797854335218850E+02;
  COFD[          73] =   0.4905020138101169E+01;
  COFD[          74] =  -0.4293084201751438E+00;
  COFD[          75] =   0.1891302999394448E-01;
  COFD[          76] =  -0.1797854335218850E+02;
  COFD[          77] =   0.4905020138101169E+01;
  COFD[          78] =  -0.4293084201751438E+00;
  COFD[          79] =   0.1891302999394448E-01;
  COFD[          80] =  -0.1783622363376704E+02;
  COFD[          81] =   0.4831292682846398E+01;
  COFD[          82] =  -0.4169745817203137E+00;
  COFD[          83] =   0.1825729345293071E-01;
  COFD[          84] =  -0.1527556886659549E+02;
  COFD[          85] =   0.3862563385245359E+01;
  COFD[          86] =  -0.2980840675934164E+00;
  COFD[          87] =   0.1340588281258234E-01;
  COFD[          88] =  -0.1528921665984086E+02;
  COFD[          89] =   0.3867752822044799E+01;
  COFD[          90] =  -0.2988043102116158E+00;
  COFD[          91] =   0.1343910541069478E-01;
  COFD[          92] =  -0.1530201153719761E+02;
  COFD[          93] =   0.3872623565106621E+01;
  COFD[          94] =  -0.2994802960421091E+00;
  COFD[          95] =   0.1347028558095244E-01;
  COFD[          96] =  -0.1646041046594534E+02;
  COFD[          97] =   0.4329097783536894E+01;
  COFD[          98] =  -0.3588538378541796E+00;
  COFD[          99] =   0.1604603061706265E-01;
  COFD[         100] =  -0.1615267929807825E+02;
  COFD[         101] =   0.4174121490551778E+01;
  COFD[         102] =  -0.3383475624805202E+00;
  COFD[         103] =   0.1513993577326728E-01;
  COFD[         104] =  -0.1616515483252968E+02;
  COFD[         105] =   0.4178910660883203E+01;
  COFD[         106] =  -0.3390062877759687E+00;
  COFD[         107] =   0.1517006334795817E-01;
  COFD[         108] =  -0.1309372900688130E+02;
  COFD[         109] =   0.3192057362682311E+01;
  COFD[         110] =  -0.2080508536154327E+00;
  COFD[         111] =   0.9377275278300783E-02;
  COFD[         112] =  -0.1769072264213409E+02;
  COFD[         113] =   0.4739685996132764E+01;
  COFD[         114] =  -0.4053730868412535E+00;
  COFD[         115] =   0.1776343130216211E-01;
  COFD[         116] =  -0.1769072264213409E+02;
  COFD[         117] =   0.4739685996132764E+01;
  COFD[         118] =  -0.4053730868412535E+00;
  COFD[         119] =   0.1776343130216211E-01;
  COFD[         120] =  -0.1266066516595179E+02;
  COFD[         121] =   0.2898076057146818E+01;
  COFD[         122] =  -0.1700453759318706E+00;
  COFD[         123] =   0.7738690296212197E-02;
  COFD[         124] =  -0.1342264266757377E+02;
  COFD[         125] =   0.3221777556990338E+01;
  COFD[         126] =  -0.2128673394230190E+00;
  COFD[         127] =   0.9627787551744238E-02;
  COFD[         128] =  -0.1201910564637032E+02;
  COFD[         129] =   0.3033139269467014E+01;
  COFD[         130] =  -0.1860150369684279E+00;
  COFD[         131] =   0.8361628520974778E-02;
  COFD[         132] =  -0.1519084593013382E+02;
  COFD[         133] =   0.4384306309333000E+01;
  COFD[         134] =  -0.3557550068895214E+00;
  COFD[         135] =   0.1548014912218549E-01;
  COFD[         136] =  -0.1531287529990816E+02;
  COFD[         137] =   0.4270888061293578E+01;
  COFD[         138] =  -0.3482328696611331E+00;
  COFD[         139] =   0.1545097630860032E-01;
  COFD[         140] =  -0.1761107854705398E+02;
  COFD[         141] =   0.5048254268555682E+01;
  COFD[         142] =  -0.4490719655111177E+00;
  COFD[         143] =   0.1981332179135272E-01;
  COFD[         144] =  -0.1533072430565730E+02;
  COFD[         145] =   0.4277567720643019E+01;
  COFD[         146] =  -0.3491407301790579E+00;
  COFD[         147] =   0.1549202933974097E-01;
  COFD[         148] =  -0.1651315089501800E+02;
  COFD[         149] =   0.4220701654162998E+01;
  COFD[         150] =  -0.2825947916575933E+00;
  COFD[         151] =   0.1018065117940972E-01;
  COFD[         152] =  -0.1762111333528786E+02;
  COFD[         153] =   0.5052180357519856E+01;
  COFD[         154] =  -0.4496001366665920E+00;
  COFD[         155] =   0.1983696581079775E-01;
  COFD[         156] =  -0.1763070827156596E+02;
  COFD[         157] =   0.5055937247475734E+01;
  COFD[         158] =  -0.4501055595390027E+00;
  COFD[         159] =   0.1985959199754691E-01;
  COFD[         160] =  -0.1560173804820947E+02;
  COFD[         161] =   0.4330818915646310E+01;
  COFD[         162] =  -0.3589811094371836E+00;
  COFD[         163] =   0.1604703157842258E-01;
  COFD[         164] =  -0.1524698785298783E+02;
  COFD[         165] =   0.4246332712416798E+01;
  COFD[         166] =  -0.3448955366037200E+00;
  COFD[         167] =   0.1530006748241336E-01;
  COFD[         168] =  -0.1798724510881690E+02;
  COFD[         169] =   0.5079566464178115E+01;
  COFD[         170] =  -0.4444140659462563E+00;
  COFD[         171] =   0.1923462920552580E-01;
  COFD[         172] =  -0.1798724510881690E+02;
  COFD[         173] =   0.5079566464178115E+01;
  COFD[         174] =  -0.4444140659462563E+00;
  COFD[         175] =   0.1923462920552580E-01;
  COFD[         176] =  -0.1801636950310948E+02;
  COFD[         177] =   0.5090088857460223E+01;
  COFD[         178] =  -0.4457552498826564E+00;
  COFD[         179] =   0.1929141580461785E-01;
  COFD[         180] =  -0.1790236652645725E+02;
  COFD[         181] =   0.5054241512059871E+01;
  COFD[         182] =  -0.4413163862361035E+00;
  COFD[         183] =   0.1910841720545356E-01;
  COFD[         184] =  -0.1719132602269712E+02;
  COFD[         185] =   0.4859367175853364E+01;
  COFD[         186] =  -0.4240580935572132E+00;
  COFD[         187] =   0.1870696403708276E-01;
  COFD[         188] =  -0.1725034548774479E+02;
  COFD[         189] =   0.4566921339972978E+01;
  COFD[         190] =  -0.3509192271095338E+00;
  COFD[         191] =   0.1402685086906688E-01;
  COFD[         192] =  -0.1540645026491837E+02;
  COFD[         193] =   0.3572427121876953E+01;
  COFD[         194] =  -0.1878763535992870E+00;
  COFD[         195] =   0.5591984032298197E-02;
  COFD[         196] =  -0.1539084008249366E+02;
  COFD[         197] =   0.3565006575777335E+01;
  COFD[         198] =  -0.1867580729432450E+00;
  COFD[         199] =   0.5536884066371989E-02;
  COFD[         200] =  -0.1561125740910447E+02;
  COFD[         201] =   0.3698144089263720E+01;
  COFD[         202] =  -0.2103080207621357E+00;
  COFD[         203] =   0.6792204226221277E-02;
  COFD[         204] =  -0.1561125740910447E+02;
  COFD[         205] =   0.3698144089263720E+01;
  COFD[         206] =  -0.2103080207621357E+00;
  COFD[         207] =   0.6792204226221277E-02;
  COFD[         208] =  -0.1545485195366714E+02;
  COFD[         209] =   0.3597796727187503E+01;
  COFD[         210] =  -0.1924739949851062E+00;
  COFD[         211] =   0.5844460086632715E-02;
  COFD[         212] =  -0.1774160477281574E+02;
  COFD[         213] =   0.4779667070939399E+01;
  COFD[         214] =  -0.3856864284506084E+00;
  COFD[         215] =   0.1583474705456999E-01;
  COFD[         216] =  -0.1774507254557799E+02;
  COFD[         217] =   0.4780159131489816E+01;
  COFD[         218] =  -0.3856587269477316E+00;
  COFD[         219] =   0.1582942762638978E-01;
  COFD[         220] =  -0.1774834356244655E+02;
  COFD[         221] =   0.4780626617651506E+01;
  COFD[         222] =  -0.3856327873318030E+00;
  COFD[         223] =   0.1582440706500046E-01;
  COFD[         224] =  -0.1694497223361135E+02;
  COFD[         225] =   0.4353290721870846E+01;
  COFD[         226] =  -0.3147441391485086E+00;
  COFD[         227] =   0.1210201943729390E-01;
  COFD[         228] =  -0.1738386031480564E+02;
  COFD[         229] =   0.4515894896789344E+01;
  COFD[         230] =  -0.3391856065220307E+00;
  COFD[         231] =   0.1329276720561771E-01;
  COFD[         232] =  -0.1738223879329547E+02;
  COFD[         233] =   0.4514299549461490E+01;
  COFD[         234] =  -0.3388651702152571E+00;
  COFD[         235] =   0.1327383213284630E-01;
  COFD[         236] =  -0.1695212975476949E+02;
  COFD[         237] =   0.4856357283707361E+01;
  COFD[         238] =  -0.4155359687875360E+00;
  COFD[         239] =   0.1799470588389149E-01;
  COFD[         240] =  -0.1510629945043326E+02;
  COFD[         241] =   0.3409059092179077E+01;
  COFD[         242] =  -0.1659094843453963E+00;
  COFD[         243] =   0.4585455445383278E-02;
  COFD[         244] =  -0.1510629945043326E+02;
  COFD[         245] =   0.3409059092179077E+01;
  COFD[         246] =  -0.1659094843453963E+00;
  COFD[         247] =   0.4585455445383278E-02;
  COFD[         248] =  -0.1711567017943538E+02;
  COFD[         249] =   0.4833264586465570E+01;
  COFD[         250] =  -0.4205908440542898E+00;
  COFD[         251] =   0.1855345683322859E-01;
  COFD[         252] =  -0.1752038279713604E+02;
  COFD[         253] =   0.4959639086677385E+01;
  COFD[         254] =  -0.4293974875408633E+00;
  COFD[         255] =   0.1860806087249651E-01;
  COFD[         256] =  -0.1093372460559999E+02;
  COFD[         257] =   0.2306133612968357E+01;
  COFD[         258] =  -0.8736427899565385E-01;
  COFD[         259] =   0.3897526119633928E-02;
  COFD[         260] =  -0.1531287529990816E+02;
  COFD[         261] =   0.4270888061293578E+01;
  COFD[         262] =  -0.3482328696611331E+00;
  COFD[         263] =   0.1545097630860032E-01;
  COFD[         264] =  -0.1356847565002614E+02;
  COFD[         265] =   0.3062175626208075E+01;
  COFD[         266] =  -0.1889691489381823E+00;
  COFD[         267] =   0.8454139853458284E-02;
  COFD[         268] =  -0.1511165027002312E+02;
  COFD[         269] =   0.3515441252015493E+01;
  COFD[         270] =  -0.2491198761221894E+00;
  COFD[         271] =   0.1111625360314578E-01;
  COFD[         272] =  -0.1358466752355666E+02;
  COFD[         273] =   0.3062671492876047E+01;
  COFD[         274] =  -0.1890385291581307E+00;
  COFD[         275] =   0.8457362732499073E-02;
  COFD[         276] =  -0.1906508600184398E+02;
  COFD[         277] =   0.4989055420187643E+01;
  COFD[         278] =  -0.4192832417141421E+00;
  COFD[         279] =   0.1761870656062123E-01;
  COFD[         280] =  -0.1513524421958843E+02;
  COFD[         281] =   0.3523297948981071E+01;
  COFD[         282] =  -0.2502118585975575E+00;
  COFD[         283] =   0.1116667641067510E-01;
  COFD[         284] =  -0.1515824683455671E+02;
  COFD[         285] =   0.3531011422568503E+01;
  COFD[         286] =  -0.2512839088019835E+00;
  COFD[         287] =   0.1121617798528330E-01;
  COFD[         288] =  -0.1340378974537447E+02;
  COFD[         289] =   0.2956605082389272E+01;
  COFD[         290] =  -0.1751292924768263E+00;
  COFD[         291] =   0.7848688796275457E-02;
  COFD[         292] =  -0.1352724646961571E+02;
  COFD[         293] =   0.3067755508215616E+01;
  COFD[         294] =  -0.1897498600053720E+00;
  COFD[         295] =   0.8490405512152714E-02;
  COFD[         296] =  -0.1558834457822415E+02;
  COFD[         297] =   0.3699567664433178E+01;
  COFD[         298] =  -0.2706504848826866E+00;
  COFD[         299] =   0.1194595865045622E-01;
  COFD[         300] =  -0.1558834457822415E+02;
  COFD[         301] =   0.3699567664433178E+01;
  COFD[         302] =  -0.2706504848826866E+00;
  COFD[         303] =   0.1194595865045622E-01;
  COFD[         304] =  -0.1561687122378518E+02;
  COFD[         305] =   0.3703886367486989E+01;
  COFD[         306] =  -0.2712410029515650E+00;
  COFD[         307] =   0.1197280722366742E-01;
  COFD[         308] =  -0.1556217965348612E+02;
  COFD[         309] =   0.3683628994574988E+01;
  COFD[         310] =  -0.2686556693855442E+00;
  COFD[         311] =   0.1186253607235680E-01;
  COFD[         312] =  -0.1478198899238034E+02;
  COFD[         313] =   0.3371092581023257E+01;
  COFD[         314] =  -0.2298567856863946E+00;
  COFD[         315] =   0.1025846107153025E-01;
  COFD[         316] =  -0.1793673861488134E+02;
  COFD[         317] =   0.4489875914439151E+01;
  COFD[         318] =  -0.3697556337020395E+00;
  COFD[         319] =   0.1609253879718407E-01;
  COFD[         320] =  -0.1915075646861156E+02;
  COFD[         321] =   0.4868838487330766E+01;
  COFD[         322] =  -0.4046210142562804E+00;
  COFD[         323] =   0.1701122981099066E-01;
  COFD[         324] =  -0.1914986770399175E+02;
  COFD[         325] =   0.4864899621234062E+01;
  COFD[         326] =  -0.4039220297382060E+00;
  COFD[         327] =   0.1697258735720344E-01;
  COFD[         328] =  -0.1878355293596396E+02;
  COFD[         329] =   0.4745970603113269E+01;
  COFD[         330] =  -0.3924372163577848E+00;
  COFD[         331] =   0.1663659729501466E-01;
  COFD[         332] =  -0.1878355293596396E+02;
  COFD[         333] =   0.4745970603113269E+01;
  COFD[         334] =  -0.3924372163577848E+00;
  COFD[         335] =   0.1663659729501466E-01;
  COFD[         336] =  -0.1907199325752651E+02;
  COFD[         337] =   0.4830462464213029E+01;
  COFD[         338] =  -0.3998168598650015E+00;
  COFD[         339] =   0.1680839848129470E-01;
  COFD[         340] =  -0.1728281584130743E+02;
  COFD[         341] =   0.4249575630985691E+01;
  COFD[         342] =  -0.3404513157409610E+00;
  COFD[         343] =   0.1490330472499406E-01;
  COFD[         344] =  -0.1731509452354737E+02;
  COFD[         345] =   0.4259553692469185E+01;
  COFD[         346] =  -0.3417709853226915E+00;
  COFD[         347] =   0.1496136243888912E-01;
  COFD[         348] =  -0.1734665788392033E+02;
  COFD[         349] =   0.4269427891180465E+01;
  COFD[         350] =  -0.3430766910038187E+00;
  COFD[         351] =   0.1501879627721415E-01;
  COFD[         352] =  -0.1800091724620567E+02;
  COFD[         353] =   0.4493072014688184E+01;
  COFD[         354] =  -0.3678451189139861E+00;
  COFD[         355] =   0.1591448159607913E-01;
  COFD[         356] =  -0.1802433609839966E+02;
  COFD[         357] =   0.4483646303442399E+01;
  COFD[         358] =  -0.3688151334706006E+00;
  COFD[         359] =   0.1604826293273014E-01;
  COFD[         360] =  -0.1805221831789452E+02;
  COFD[         361] =   0.4492279525932131E+01;
  COFD[         362] =  -0.3699242263001061E+00;
  COFD[         363] =   0.1609559465701492E-01;
  COFD[         364] =  -0.1575708126216198E+02;
  COFD[         365] =   0.3859037941929713E+01;
  COFD[         366] =  -0.2920671356487180E+00;
  COFD[         367] =   0.1290449436293418E-01;
  COFD[         368] =  -0.1900240647811995E+02;
  COFD[         369] =   0.4757008879827381E+01;
  COFD[         370] =  -0.3907894080963009E+00;
  COFD[         371] =   0.1643024015529323E-01;
  COFD[         372] =  -0.1900240647811995E+02;
  COFD[         373] =   0.4757008879827381E+01;
  COFD[         374] =  -0.3907894080963009E+00;
  COFD[         375] =   0.1643024015529323E-01;
  COFD[         376] =  -0.1474819919893963E+02;
  COFD[         377] =   0.3361311502126538E+01;
  COFD[         378] =  -0.2285465363406641E+00;
  COFD[         379] =   0.1019990809984005E-01;
  COFD[         380] =  -0.1581858256688188E+02;
  COFD[         381] =   0.3775899039344615E+01;
  COFD[         382] =  -0.2814687920927101E+00;
  COFD[         383] =   0.1245218573218981E-01;
  COFD[         384] =  -0.1270707915711888E+02;
  COFD[         385] =   0.2929619774221139E+01;
  COFD[         386] =  -0.1739670601708100E+00;
  COFD[         387] =   0.7901198812704345E-02;
  COFD[         388] =  -0.1761107854705398E+02;
  COFD[         389] =   0.5048254268555682E+01;
  COFD[         390] =  -0.4490719655111177E+00;
  COFD[         391] =   0.1981332179135272E-01;
  COFD[         392] =  -0.1511165027002312E+02;
  COFD[         393] =   0.3515441252015493E+01;
  COFD[         394] =  -0.2491198761221894E+00;
  COFD[         395] =   0.1111625360314578E-01;
  COFD[         396] =  -0.1605076758771229E+02;
  COFD[         397] =   0.3688226344868140E+01;
  COFD[         398] =  -0.2690947493495641E+00;
  COFD[         399] =   0.1187504521154671E-01;
  COFD[         400] =  -0.1509620843583920E+02;
  COFD[         401] =   0.3500304707531062E+01;
  COFD[         402] =  -0.2470160084537919E+00;
  COFD[         403] =   0.1101910392929726E-01;
  COFD[         404] =  -0.2029800570223752E+02;
  COFD[         405] =   0.5167991782317953E+01;
  COFD[         406] =  -0.4260628585538284E+00;
  COFD[         407] =   0.1725358239585169E-01;
  COFD[         408] =  -0.1605887758547131E+02;
  COFD[         409] =   0.3688399632843664E+01;
  COFD[         410] =  -0.2691183502617316E+00;
  COFD[         411] =   0.1187611401222301E-01;
  COFD[         412] =  -0.1606742328871227E+02;
  COFD[         413] =   0.3688898239072069E+01;
  COFD[         414] =  -0.2691862577212903E+00;
  COFD[         415] =   0.1187918929470605E-01;
  COFD[         416] =  -0.1480659595361166E+02;
  COFD[         417] =   0.3371224676621789E+01;
  COFD[         418] =  -0.2300778346811563E+00;
  COFD[         419] =   0.1027669190040931E-01;
  COFD[         420] =  -0.1516681611308053E+02;
  COFD[         421] =   0.3569087909107446E+01;
  COFD[         422] =  -0.2565755708070836E+00;
  COFD[         423] =   0.1146050720294671E-01;
  COFD[         424] =  -0.1707915522153321E+02;
  COFD[         425] =   0.4132970234506635E+01;
  COFD[         426] =  -0.3260547355537204E+00;
  COFD[         427] =   0.1431106723651207E-01;
  COFD[         428] =  -0.1707915522153321E+02;
  COFD[         429] =   0.4132970234506635E+01;
  COFD[         430] =  -0.3260547355537204E+00;
  COFD[         431] =   0.1431106723651207E-01;
  COFD[         432] =  -0.1707094794267795E+02;
  COFD[         433] =   0.4119810624341889E+01;
  COFD[         434] =  -0.3242844616145015E+00;
  COFD[         435] =   0.1423184749658152E-01;
  COFD[         436] =  -0.1696131964138645E+02;
  COFD[         437] =   0.4074167708212651E+01;
  COFD[         438] =  -0.3182576642303228E+00;
  COFD[         439] =   0.1396612551510442E-01;
  COFD[         440] =  -0.1582244348422915E+02;
  COFD[         441] =   0.3599211636009877E+01;
  COFD[         442] =  -0.2581898369286856E+00;
  COFD[         443] =   0.1143096575439008E-01;
  COFD[         444] =  -0.1868514242925209E+02;
  COFD[         445] =   0.4565512377700616E+01;
  COFD[         446] =  -0.3754538805580755E+00;
  COFD[         447] =   0.1617305773699574E-01;
  COFD[         448] =  -0.2031657960333580E+02;
  COFD[         449] =   0.5113614640951361E+01;
  COFD[         450] =  -0.4319293017089550E+00;
  COFD[         451] =   0.1802805482147992E-01;
  COFD[         452] =  -0.2032689236579443E+02;
  COFD[         453] =   0.5114221260960875E+01;
  COFD[         454] =  -0.4320136109530544E+00;
  COFD[         455] =   0.1803190494311217E-01;
  COFD[         456] =  -0.2016013281728147E+02;
  COFD[         457] =   0.5063566265537697E+01;
  COFD[         458] =  -0.4281515792120217E+00;
  COFD[         459] =   0.1797180686456706E-01;
  COFD[         460] =  -0.2016013281728147E+02;
  COFD[         461] =   0.5063566265537697E+01;
  COFD[         462] =  -0.4281515792120217E+00;
  COFD[         463] =   0.1797180686456706E-01;
  COFD[         464] =  -0.2028793707722566E+02;
  COFD[         465] =   0.5097377789913596E+01;
  COFD[         466] =  -0.4307531399728672E+00;
  COFD[         467] =   0.1801467570282063E-01;
  COFD[         468] =  -0.1816374974359966E+02;
  COFD[         469] =   0.4404539112241986E+01;
  COFD[         470] =  -0.3569514300098177E+00;
  COFD[         471] =   0.1547163489736868E-01;
  COFD[         472] =  -0.1817374730686181E+02;
  COFD[         473] =   0.4404209750065649E+01;
  COFD[         474] =  -0.3569181269134391E+00;
  COFD[         475] =   0.1547062305483025E-01;
  COFD[         476] =  -0.1818389332535593E+02;
  COFD[         477] =   0.4404158898744668E+01;
  COFD[         478] =  -0.3569180632138942E+00;
  COFD[         479] =   0.1547091340944365E-01;
  COFD[         480] =  -0.1901929468754062E+02;
  COFD[         481] =   0.4693839754776671E+01;
  COFD[         482] =  -0.3901000075815114E+00;
  COFD[         483] =   0.1672631492507185E-01;
  COFD[         484] =  -0.1881741117557303E+02;
  COFD[         485] =   0.4595519619722809E+01;
  COFD[         486] =  -0.3791053358880415E+00;
  COFD[         487] =   0.1632232422675412E-01;
  COFD[         488] =  -0.1882824574536174E+02;
  COFD[         489] =   0.4596254296231404E+01;
  COFD[         490] =  -0.3791928797888692E+00;
  COFD[         491] =   0.1632573982847240E-01;
  COFD[         492] =  -0.1687140720665078E+02;
  COFD[         493] =   0.4075868327717259E+01;
  COFD[         494] =  -0.3175542600319969E+00;
  COFD[         495] =   0.1389707797540003E-01;
  COFD[         496] =  -0.2009809718705050E+02;
  COFD[         497] =   0.4977640189149046E+01;
  COFD[         498] =  -0.4165323971276452E+00;
  COFD[         499] =   0.1744635825561870E-01;
  COFD[         500] =  -0.2009809718705050E+02;
  COFD[         501] =   0.4977640189149046E+01;
  COFD[         502] =  -0.4165323971276452E+00;
  COFD[         503] =   0.1744635825561870E-01;
  COFD[         504] =  -0.1578873902175925E+02;
  COFD[         505] =   0.3589258101555060E+01;
  COFD[         506] =  -0.2568718883472193E+00;
  COFD[         507] =   0.1137271971440306E-01;
  COFD[         508] =  -0.1677455230511400E+02;
  COFD[         509] =   0.3946214344521416E+01;
  COFD[         510] =  -0.3012613112287137E+00;
  COFD[         511] =   0.1321325422150532E-01;
  COFD[         512] =  -0.1093673058791694E+02;
  COFD[         513] =   0.2306118133758978E+01;
  COFD[         514] =  -0.8736238957367691E-01;
  COFD[         515] =   0.3897448280154413E-02;
  COFD[         516] =  -0.1533072430565730E+02;
  COFD[         517] =   0.4277567720643019E+01;
  COFD[         518] =  -0.3491407301790579E+00;
  COFD[         519] =   0.1549202933974097E-01;
  COFD[         520] =  -0.1358466752355666E+02;
  COFD[         521] =   0.3062671492876047E+01;
  COFD[         522] =  -0.1890385291581307E+00;
  COFD[         523] =   0.8457362732499073E-02;
  COFD[         524] =  -0.1509620843583920E+02;
  COFD[         525] =   0.3500304707531062E+01;
  COFD[         526] =  -0.2470160084537919E+00;
  COFD[         527] =   0.1101910392929726E-01;
  COFD[         528] =  -0.1359902342804023E+02;
  COFD[         529] =   0.3062175626208111E+01;
  COFD[         530] =  -0.1889691489381874E+00;
  COFD[         531] =   0.8454139853458534E-02;
  COFD[         532] =  -0.1908171031786331E+02;
  COFD[         533] =   0.4989633286602976E+01;
  COFD[         534] =  -0.4194071192530971E+00;
  COFD[         535] =   0.1762632910887997E-01;
  COFD[         536] =  -0.1511946285859442E+02;
  COFD[         537] =   0.3507927470352835E+01;
  COFD[         538] =  -0.2480755290441053E+00;
  COFD[         539] =   0.1106802953504252E-01;
  COFD[         540] =  -0.1514219804803715E+02;
  COFD[         541] =   0.3515441252015498E+01;
  COFD[         542] =  -0.2491198761221899E+00;
  COFD[         543] =   0.1111625360314580E-01;
  COFD[         544] =  -0.1342265746209028E+02;
  COFD[         545] =   0.2959193932129379E+01;
  COFD[         546] =  -0.1754915954841434E+00;
  COFD[         547] =   0.7865516993304794E-02;
  COFD[         548] =  -0.1354935397948868E+02;
  COFD[         549] =   0.3071473642095726E+01;
  COFD[         550] =  -0.1902700720448971E+00;
  COFD[         551] =   0.8514570131842562E-02;
  COFD[         552] =  -0.1559686009787643E+02;
  COFD[         553] =   0.3697254029246602E+01;
  COFD[         554] =  -0.2703332731709061E+00;
  COFD[         555] =   0.1193149792533076E-01;
  COFD[         556] =  -0.1559686009787643E+02;
  COFD[         557] =   0.3697254029246602E+01;
  COFD[         558] =  -0.2703332731709061E+00;
  COFD[         559] =   0.1193149792533076E-01;
  COFD[         560] =  -0.1562210625534589E+02;
  COFD[         561] =   0.3699988818188416E+01;
  COFD[         562] =  -0.2707081263160739E+00;
  COFD[         563] =   0.1194858184298051E-01;
  COFD[         564] =  -0.1556532980903788E+02;
  COFD[         565] =   0.3678661452138608E+01;
  COFD[         566] =  -0.2679777825820168E+00;
  COFD[         567] =   0.1183177506717497E-01;
  COFD[         568] =  -0.1476985219999908E+02;
  COFD[         569] =   0.3357751196108267E+01;
  COFD[         570] =  -0.2280046831004087E+00;
  COFD[         571] =   0.1017303667696902E-01;
  COFD[         572] =  -0.1792454027766008E+02;
  COFD[         573] =   0.4476326984714260E+01;
  COFD[         574] =  -0.3680191278135618E+00;
  COFD[         575] =   0.1601860672433683E-01;
  COFD[         576] =  -0.1918204714711579E+02;
  COFD[         577] =   0.4875615818554563E+01;
  COFD[         578] =  -0.4058262469897604E+00;
  COFD[         579] =   0.1707794755285167E-01;
  COFD[         580] =  -0.1918177070837826E+02;
  COFD[         581] =   0.4871873199573439E+01;
  COFD[         582] =  -0.4051602536432213E+00;
  COFD[         583] =   0.1704106539271711E-01;
  COFD[         584] =  -0.1880409131690931E+02;
  COFD[         585] =   0.4747896766189557E+01;
  COFD[         586] =  -0.3929513558466455E+00;
  COFD[         587] =   0.1667071927589917E-01;
  COFD[         588] =  -0.1880409131690931E+02;
  COFD[         589] =   0.4747896766189557E+01;
  COFD[         590] =  -0.3929513558466455E+00;
  COFD[         591] =   0.1667071927589917E-01;
  COFD[         592] =  -0.1910332751367901E+02;
  COFD[         593] =   0.4837040211978248E+01;
  COFD[         594] =  -0.4010050297336601E+00;
  COFD[         595] =   0.1687475885991967E-01;
  COFD[         596] =  -0.1726463880632943E+02;
  COFD[         597] =   0.4234775238868422E+01;
  COFD[         598] =  -0.3384933692140595E+00;
  COFD[         599] =   0.1481714586415791E-01;
  COFD[         600] =  -0.1729594517981932E+02;
  COFD[         601] =   0.4244246410490659E+01;
  COFD[         602] =  -0.3397463840897515E+00;
  COFD[         603] =   0.1487228746520592E-01;
  COFD[         604] =  -0.1732667376192460E+02;
  COFD[         605] =   0.4253674558730024E+01;
  COFD[         606] =  -0.3409934577271315E+00;
  COFD[         607] =   0.1492715703864712E-01;
  COFD[         608] =  -0.1799417924090561E+02;
  COFD[         609] =   0.4483237269191199E+01;
  COFD[         610] =  -0.3666457816700467E+00;
  COFD[         611] =   0.1586624138953086E-01;
  COFD[         612] =  -0.1800542957417654E+02;
  COFD[         613] =   0.4468504271448043E+01;
  COFD[         614] =  -0.3668692513471102E+00;
  COFD[         615] =   0.1596519376748332E-01;
  COFD[         616] =  -0.1803293840758895E+02;
  COFD[         617] =   0.4476898762718537E+01;
  COFD[         618] =  -0.3679481179513786E+00;
  COFD[         619] =   0.1601125468656841E-01;
  COFD[         620] =  -0.1574882548231690E+02;
  COFD[         621] =   0.3846525250144765E+01;
  COFD[         622] =  -0.2903598878615679E+00;
  COFD[         623] =   0.1282703770058140E-01;
  COFD[         624] =  -0.1902751852931795E+02;
  COFD[         625] =   0.4760378344803059E+01;
  COFD[         626] =  -0.3915558236239631E+00;
  COFD[         627] =   0.1647818289394700E-01;
  COFD[         628] =  -0.1902751852931795E+02;
  COFD[         629] =   0.4760378344803059E+01;
  COFD[         630] =  -0.3915558236239631E+00;
  COFD[         631] =   0.1647818289394700E-01;
  COFD[         632] =  -0.1473661449825940E+02;
  COFD[         633] =   0.3348204558833826E+01;
  COFD[         634] =  -0.2267271723233180E+00;
  COFD[         635] =   0.1011600240359858E-01;
  COFD[         636] =  -0.1580533640016843E+02;
  COFD[         637] =   0.3761415063014861E+01;
  COFD[         638] =  -0.2795015971789234E+00;
  COFD[         639] =   0.1236332420639666E-01;
  COFD[         640] =  -0.1818497883283899E+02;
  COFD[         641] =   0.5043624597310085E+01;
  COFD[         642] =  -0.4377066809210735E+00;
  COFD[         643] =   0.1887643866477002E-01;
  COFD[         644] =  -0.1651315089501800E+02;
  COFD[         645] =   0.4220701654162998E+01;
  COFD[         646] =  -0.2825947916575933E+00;
  COFD[         647] =   0.1018065117940972E-01;
  COFD[         648] =  -0.1906508600184398E+02;
  COFD[         649] =   0.4989055420187643E+01;
  COFD[         650] =  -0.4192832417141421E+00;
  COFD[         651] =   0.1761870656062123E-01;
  COFD[         652] =  -0.2029800570223752E+02;
  COFD[         653] =   0.5167991782317953E+01;
  COFD[         654] =  -0.4260628585538284E+00;
  COFD[         655] =   0.1725358239585169E-01;
  COFD[         656] =  -0.1908171031786331E+02;
  COFD[         657] =   0.4989633286602976E+01;
  COFD[         658] =  -0.4194071192530971E+00;
  COFD[         659] =   0.1762632910887997E-01;
  COFD[         660] =  -0.1180242969654581E+02;
  COFD[         661] =   0.8887910786873973E+00;
  COFD[         662] =   0.2461702200576275E+00;
  COFD[         663] =  -0.1607058638443385E-01;
  COFD[         664] =  -0.1969226895876598E+02;
  COFD[         665] =   0.4975023535625081E+01;
  COFD[         666] =  -0.4062730945223436E+00;
  COFD[         667] =   0.1660121367605663E-01;
  COFD[         668] =  -0.1967561834414118E+02;
  COFD[         669] =   0.4964896032502449E+01;
  COFD[         670] =  -0.4047371090295282E+00;
  COFD[         671] =   0.1652517237194395E-01;
  COFD[         672] =  -0.1896641926147093E+02;
  COFD[         673] =   0.4932217655116672E+01;
  COFD[         674] =  -0.4145721175243606E+00;
  COFD[         675] =   0.1752187325720189E-01;
  COFD[         676] =  -0.1898988022108976E+02;
  COFD[         677] =   0.4979356265300015E+01;
  COFD[         678] =  -0.4176235134136435E+00;
  COFD[         679] =   0.1752883066124088E-01;
  COFD[         680] =  -0.2040913503195636E+02;
  COFD[         681] =   0.5224092639655129E+01;
  COFD[         682] =  -0.4327353018434304E+00;
  COFD[         683] =   0.1753064056853202E-01;
  COFD[         684] =  -0.2040913503195636E+02;
  COFD[         685] =   0.5224092639655129E+01;
  COFD[         686] =  -0.4327353018434304E+00;
  COFD[         687] =   0.1753064056853202E-01;
  COFD[         688] =  -0.2042453573019913E+02;
  COFD[         689] =   0.5222435968901098E+01;
  COFD[         690] =  -0.4324987976330752E+00;
  COFD[         691] =   0.1751958312757915E-01;
  COFD[         692] =  -0.2052802769939090E+02;
  COFD[         693] =   0.5190809253464677E+01;
  COFD[         694] =  -0.4197156907449909E+00;
  COFD[         695] =   0.1661902490158712E-01;
  COFD[         696] =  -0.2031234004461199E+02;
  COFD[         697] =   0.5186312454417227E+01;
  COFD[         698] =  -0.4313742876779729E+00;
  COFD[         699] =   0.1759974844922337E-01;
  COFD[         700] =  -0.1924815690645235E+02;
  COFD[         701] =   0.4356731855020855E+01;
  COFD[         702] =  -0.2835834611821879E+00;
  COFD[         703] =   0.9647123135613308E-02;
  COFD[         704] =  -0.1757312087019186E+02;
  COFD[         705] =   0.3511204681563826E+01;
  COFD[         706] =  -0.1497966127823777E+00;
  COFD[         707] =   0.2966580974304371E-02;
  COFD[         708] =  -0.1757470287826831E+02;
  COFD[         709] =   0.3509436887842453E+01;
  COFD[         710] =  -0.1496102927025162E+00;
  COFD[         711] =   0.2961974572337444E-02;
  COFD[         712] =  -0.1619587514529141E+02;
  COFD[         713] =   0.2843235032975332E+01;
  COFD[         714] =  -0.5104314541884596E-01;
  COFD[         715] =  -0.1707937587486378E-02;
  COFD[         716] =  -0.1619587514529141E+02;
  COFD[         717] =   0.2843235032975332E+01;
  COFD[         718] =  -0.5104314541884596E-01;
  COFD[         719] =  -0.1707937587486378E-02;
  COFD[         720] =  -0.1783328418687383E+02;
  COFD[         721] =   0.3623064754744803E+01;
  COFD[         722] =  -0.1668883567851333E+00;
  COFD[         723] =   0.3813867862735125E-02;
  COFD[         724] =  -0.2027561345681607E+02;
  COFD[         725] =   0.4923011006260475E+01;
  COFD[         726] =  -0.3762758154685706E+00;
  COFD[         727] =   0.1440493997760856E-01;
  COFD[         728] =  -0.2024506268677845E+02;
  COFD[         729] =   0.4905829896937226E+01;
  COFD[         730] =  -0.3737705712203094E+00;
  COFD[         731] =   0.1428535217779750E-01;
  COFD[         732] =  -0.2021382909884265E+02;
  COFD[         733] =   0.4888560960839369E+01;
  COFD[         734] =  -0.3712523102934615E+00;
  COFD[         735] =   0.1416513273755347E-01;
  COFD[         736] =  -0.1988688758541057E+02;
  COFD[         737] =   0.4655348588622681E+01;
  COFD[         738] =  -0.3287200044949488E+00;
  COFD[         739] =   0.1186464366263145E-01;
  COFD[         740] =  -0.2010033725040245E+02;
  COFD[         741] =   0.4740608584978043E+01;
  COFD[         742] =  -0.3442903017689571E+00;
  COFD[         743] =   0.1271241178230842E-01;
  COFD[         744] =  -0.2007381993778366E+02;
  COFD[         745] =   0.4726146050811121E+01;
  COFD[         746] =  -0.3422109775805109E+00;
  COFD[         747] =   0.1261488904944177E-01;
  COFD[         748] =  -0.1993851215676389E+02;
  COFD[         749] =   0.5098663262971057E+01;
  COFD[         750] =  -0.4129849728390237E+00;
  COFD[         751] =   0.1653578043838391E-01;
  COFD[         752] =  -0.1836390465233753E+02;
  COFD[         753] =   0.3824966367676553E+01;
  COFD[         754] =  -0.1992927348787672E+00;
  COFD[         755] =   0.5478606239007274E-02;
  COFD[         756] =  -0.1836390465233753E+02;
  COFD[         757] =   0.3824966367676553E+01;
  COFD[         758] =  -0.1992927348787672E+00;
  COFD[         759] =   0.5478606239007274E-02;
  COFD[         760] =  -0.2024331300258610E+02;
  COFD[         761] =   0.5167058816761965E+01;
  COFD[         762] =  -0.4292618020154789E+00;
  COFD[         763] =   0.1752039422653532E-01;
  COFD[         764] =  -0.1991427613539219E+02;
  COFD[         765] =   0.4985074718206674E+01;
  COFD[         766] =  -0.3988046631693363E+00;
  COFD[         767] =   0.1592993698509921E-01;
  COFD[         768] =  -0.1271124360389947E+02;
  COFD[         769] =   0.2931047949505553E+01;
  COFD[         770] =  -0.1741679825526774E+00;
  COFD[         771] =   0.7910585398301936E-02;
  COFD[         772] =  -0.1762111333528786E+02;
  COFD[         773] =   0.5052180357519856E+01;
  COFD[         774] =  -0.4496001366665920E+00;
  COFD[         775] =   0.1983696581079775E-01;
  COFD[         776] =  -0.1513524421958843E+02;
  COFD[         777] =   0.3523297948981071E+01;
  COFD[         778] =  -0.2502118585975575E+00;
  COFD[         779] =   0.1116667641067510E-01;
  COFD[         780] =  -0.1605887758547131E+02;
  COFD[         781] =   0.3688399632843664E+01;
  COFD[         782] =  -0.2691183502617316E+00;
  COFD[         783] =   0.1187611401222301E-01;
  COFD[         784] =  -0.1511946285859442E+02;
  COFD[         785] =   0.3507927470352835E+01;
  COFD[         786] =  -0.2480755290441053E+00;
  COFD[         787] =   0.1106802953504252E-01;
  COFD[         788] =  -0.1969226895876598E+02;
  COFD[         789] =   0.4975023535625081E+01;
  COFD[         790] =  -0.4062730945223436E+00;
  COFD[         791] =   0.1660121367605663E-01;
  COFD[         792] =  -0.1606627473213039E+02;
  COFD[         793] =   0.3688226344868143E+01;
  COFD[         794] =  -0.2690947493495643E+00;
  COFD[         795] =   0.1187504521154673E-01;
  COFD[         796] =  -0.1607413035852147E+02;
  COFD[         797] =   0.3688389366224249E+01;
  COFD[         798] =  -0.2691169520025368E+00;
  COFD[         799] =   0.1187605069009256E-01;
  COFD[         800] =  -0.1482502680337964E+02;
  COFD[         801] =   0.3377285630938066E+01;
  COFD[         802] =  -0.2309189157634036E+00;
  COFD[         803] =   0.1031547111935213E-01;
  COFD[         804] =  -0.1519066098667452E+02;
  COFD[         805] =   0.3577339223376338E+01;
  COFD[         806] =  -0.2577222156910252E+00;
  COFD[         807] =   0.1151344824477554E-01;
  COFD[         808] =  -0.1709853282395253E+02;
  COFD[         809] =   0.4139002691011624E+01;
  COFD[         810] =  -0.3268662406958736E+00;
  COFD[         811] =   0.1434738242897968E-01;
  COFD[         812] =  -0.1709853282395253E+02;
  COFD[         813] =   0.4139002691011624E+01;
  COFD[         814] =  -0.3268662406958736E+00;
  COFD[         815] =   0.1434738242897968E-01;
  COFD[         816] =  -0.1709003104221693E+02;
  COFD[         817] =   0.4125627998285617E+01;
  COFD[         818] =  -0.3250670331906144E+00;
  COFD[         819] =   0.1426686739847163E-01;
  COFD[         820] =  -0.1698006191382480E+02;
  COFD[         821] =   0.4079755891690070E+01;
  COFD[         822] =  -0.3190093940831370E+00;
  COFD[         823] =   0.1399976445173436E-01;
  COFD[         824] =  -0.1583284545979817E+02;
  COFD[         825] =   0.3600569993944359E+01;
  COFD[         826] =  -0.2583773478248265E+00;
  COFD[         827] =   0.1143956672273590E-01;
  COFD[         828] =  -0.1868972720602865E+02;
  COFD[         829] =   0.4564003101584411E+01;
  COFD[         830] =  -0.3752914154261894E+00;
  COFD[         831] =   0.1616756912994680E-01;
  COFD[         832] =  -0.2032173485009814E+02;
  COFD[         833] =   0.5112656067238100E+01;
  COFD[         834] =  -0.4317894593025599E+00;
  COFD[         835] =   0.1802136540679988E-01;
  COFD[         836] =  -0.2033308000404864E+02;
  COFD[         837] =   0.5113689932626547E+01;
  COFD[         838] =  -0.4319400850127593E+00;
  COFD[         839] =   0.1802856188319298E-01;
  COFD[         840] =  -0.1994316533928780E+02;
  COFD[         841] =   0.4989293390342089E+01;
  COFD[         842] =  -0.4200959859596722E+00;
  COFD[         843] =   0.1768643162186107E-01;
  COFD[         844] =  -0.1994316533928780E+02;
  COFD[         845] =   0.4989293390342089E+01;
  COFD[         846] =  -0.4200959859596722E+00;
  COFD[         847] =   0.1768643162186107E-01;
  COFD[         848] =  -0.2029614681465456E+02;
  COFD[         849] =   0.5097689110913429E+01;
  COFD[         850] =  -0.4308107142909775E+00;
  COFD[         851] =   0.1801795640022364E-01;
  COFD[         852] =  -0.1817178289920526E+02;
  COFD[         853] =   0.4404990527806040E+01;
  COFD[         854] =  -0.3570007541989366E+00;
  COFD[         855] =   0.1547334532679651E-01;
  COFD[         856] =  -0.1818135452916891E+02;
  COFD[         857] =   0.4404445109750770E+01;
  COFD[         858] =  -0.3569414786935936E+00;
  COFD[         859] =   0.1547130685514495E-01;
  COFD[         860] =  -0.1819107613580373E+02;
  COFD[         861] =   0.4404180695515474E+01;
  COFD[         862] =  -0.3569157439691602E+00;
  COFD[         863] =   0.1547058259292524E-01;
  COFD[         864] =  -0.1902570873205879E+02;
  COFD[         865] =   0.4693487292163261E+01;
  COFD[         866] =  -0.3900471625097403E+00;
  COFD[         867] =   0.1672372332121293E-01;
  COFD[         868] =  -0.1882326100183751E+02;
  COFD[         869] =   0.4594957876626870E+01;
  COFD[         870] =  -0.3790360690538448E+00;
  COFD[         871] =   0.1631950746180167E-01;
  COFD[         872] =  -0.1883388798798641E+02;
  COFD[         873] =   0.4595580957493091E+01;
  COFD[         874] =  -0.3791127554800276E+00;
  COFD[         875] =   0.1632261913288707E-01;
  COFD[         876] =  -0.1687995772389570E+02;
  COFD[         877] =   0.4075816169135930E+01;
  COFD[         878] =  -0.3175458358276380E+00;
  COFD[         879] =   0.1389664535256806E-01;
  COFD[         880] =  -0.2011350568921456E+02;
  COFD[         881] =   0.4980977035855766E+01;
  COFD[         882] =  -0.4170734091977097E+00;
  COFD[         883] =   0.1747455638443835E-01;
  COFD[         884] =  -0.2011350568921456E+02;
  COFD[         885] =   0.4980977035855766E+01;
  COFD[         886] =  -0.4170734091977097E+00;
  COFD[         887] =   0.1747455638443835E-01;
  COFD[         888] =  -0.1579929195048582E+02;
  COFD[         889] =   0.3590679879226839E+01;
  COFD[         890] =  -0.2570681340587627E+00;
  COFD[         891] =   0.1138172060899465E-01;
  COFD[         892] =  -0.1677798202265061E+02;
  COFD[         893] =   0.3944105793182093E+01;
  COFD[         894] =  -0.3009768138615855E+00;
  COFD[         895] =   0.1320048849753941E-01;
  COFD[         896] =  -0.1271518733681371E+02;
  COFD[         897] =   0.2932402241140174E+01;
  COFD[         898] =  -0.1743585082612658E+00;
  COFD[         899] =   0.7919486191765343E-02;
  COFD[         900] =  -0.1763070827156596E+02;
  COFD[         901] =   0.5055937247475734E+01;
  COFD[         902] =  -0.4501055595390027E+00;
  COFD[         903] =   0.1985959199754691E-01;
  COFD[         904] =  -0.1515824683455671E+02;
  COFD[         905] =   0.3531011422568503E+01;
  COFD[         906] =  -0.2512839088019835E+00;
  COFD[         907] =   0.1121617798528330E-01;
  COFD[         908] =  -0.1606742328871227E+02;
  COFD[         909] =   0.3688898239072069E+01;
  COFD[         910] =  -0.2691862577212903E+00;
  COFD[         911] =   0.1187918929470605E-01;
  COFD[         912] =  -0.1514219804803715E+02;
  COFD[         913] =   0.3515441252015498E+01;
  COFD[         914] =  -0.2491198761221899E+00;
  COFD[         915] =   0.1111625360314580E-01;
  COFD[         916] =  -0.1967561834414118E+02;
  COFD[         917] =   0.4964896032502449E+01;
  COFD[         918] =  -0.4047371090295282E+00;
  COFD[         919] =   0.1652517237194395E-01;
  COFD[         920] =  -0.1607413035852147E+02;
  COFD[         921] =   0.3688389366224249E+01;
  COFD[         922] =  -0.2691169520025368E+00;
  COFD[         923] =   0.1187605069009256E-01;
  COFD[         924] =  -0.1608131536572632E+02;
  COFD[         925] =   0.3688226344868147E+01;
  COFD[         926] =  -0.2690947493495650E+00;
  COFD[         927] =   0.1187504521154676E-01;
  COFD[         928] =  -0.1484281133878113E+02;
  COFD[         929] =   0.3383163225244592E+01;
  COFD[         930] =  -0.2317345266445065E+00;
  COFD[         931] =   0.1035307519950314E-01;
  COFD[         932] =  -0.1521372216636978E+02;
  COFD[         933] =   0.3585353537093755E+01;
  COFD[         934] =  -0.2588358998442430E+00;
  COFD[         935] =   0.1156486668966244E-01;
  COFD[         936] =  -0.1711740515182854E+02;
  COFD[         937] =   0.4144924979490134E+01;
  COFD[         938] =  -0.3276629227930643E+00;
  COFD[         939] =   0.1438303438727467E-01;
  COFD[         940] =  -0.1711740515182854E+02;
  COFD[         941] =   0.4144924979490134E+01;
  COFD[         942] =  -0.3276629227930643E+00;
  COFD[         943] =   0.1438303438727467E-01;
  COFD[         944] =  -0.1710868501481421E+02;
  COFD[         945] =   0.4131370391002061E+01;
  COFD[         946] =  -0.3258395192522007E+00;
  COFD[         947] =   0.1430143624455670E-01;
  COFD[         948] =  -0.1699844377766270E+02;
  COFD[         949] =   0.4085300971938842E+01;
  COFD[         950] =  -0.3197553288341113E+00;
  COFD[         951] =   0.1403314440725374E-01;
  COFD[         952] =  -0.1584365035133258E+02;
  COFD[         953] =   0.3602236757402370E+01;
  COFD[         954] =  -0.2586073726747737E+00;
  COFD[         955] =   0.1145011549760728E-01;
  COFD[         956] =  -0.1869430855359004E+02;
  COFD[         957] =   0.4562621031385030E+01;
  COFD[         958] =  -0.3751419042160075E+00;
  COFD[         959] =   0.1616247820310002E-01;
  COFD[         960] =  -0.2032584895987123E+02;
  COFD[         961] =   0.5111360106328648E+01;
  COFD[         962] =  -0.4315969740328803E+00;
  COFD[         963] =   0.1801200837165252E-01;
  COFD[         964] =  -0.2033818516187319E+02;
  COFD[         965] =   0.5112803625509553E+01;
  COFD[         966] =  -0.4318111908846191E+00;
  COFD[         967] =   0.1802241389102129E-01;
  COFD[         968] =  -0.1994945707378163E+02;
  COFD[         969] =   0.4988909811915903E+01;
  COFD[         970] =  -0.4200448560138315E+00;
  COFD[         971] =   0.1768418875431581E-01;
  COFD[         972] =  -0.1994945707378163E+02;
  COFD[         973] =   0.4988909811915903E+01;
  COFD[         974] =  -0.4200448560138315E+00;
  COFD[         975] =   0.1768418875431581E-01;
  COFD[         976] =  -0.2030325152362526E+02;
  COFD[         977] =   0.5097636743004845E+01;
  COFD[         978] =  -0.4308113500403411E+00;
  COFD[         979] =   0.1801834516701545E-01;
  COFD[         980] =  -0.1817988211332206E+02;
  COFD[         981] =   0.4405585855409423E+01;
  COFD[         982] =  -0.3570675278892711E+00;
  COFD[         983] =   0.1547575275234058E-01;
  COFD[         984] =  -0.1818905402728113E+02;
  COFD[         985] =   0.4404836108862263E+01;
  COFD[         986] =  -0.3569836886410043E+00;
  COFD[         987] =   0.1547274326407077E-01;
  COFD[         988] =  -0.1819837498757435E+02;
  COFD[         989] =   0.4404368843395893E+01;
  COFD[         990] =  -0.3569335691748086E+00;
  COFD[         991] =   0.1547105514857662E-01;
  COFD[         992] =  -0.1903192515013230E+02;
  COFD[         993] =   0.4693162157411670E+01;
  COFD[         994] =  -0.3899939058159093E+00;
  COFD[         995] =   0.1672093460134929E-01;
  COFD[         996] =  -0.1882911310545702E+02;
  COFD[         997] =   0.4594509279764214E+01;
  COFD[         998] =  -0.3789782571863847E+00;
  COFD[         999] =   0.1631703818360878E-01;
  COFD[        1000] =  -0.1883954098397452E+02;
  COFD[        1001] =   0.4595024886620513E+01;
  COFD[        1002] =  -0.3790444797196689E+00;
  COFD[        1003] =   0.1631985649333879E-01;
  COFD[        1004] =  -0.1688882624224305E+02;
  COFD[        1005] =   0.4076043066784620E+01;
  COFD[        1006] =  -0.3175741894503459E+00;
  COFD[        1007] =   0.1389782472898128E-01;
  COFD[        1008] =  -0.2012780676364666E+02;
  COFD[        1009] =   0.4983954385219352E+01;
  COFD[        1010] =  -0.4175569796553835E+00;
  COFD[        1011] =   0.1749979298053702E-01;
  COFD[        1012] =  -0.2012780676364666E+02;
  COFD[        1013] =   0.4983954385219352E+01;
  COFD[        1014] =  -0.4175569796553835E+00;
  COFD[        1015] =   0.1749979298053702E-01;
  COFD[        1016] =  -0.1581023736073806E+02;
  COFD[        1017] =   0.3592405682636056E+01;
  COFD[        1018] =  -0.2573062871053550E+00;
  COFD[        1019] =   0.1139264138769127E-01;
  COFD[        1020] =  -0.1678197080138354E+02;
  COFD[        1021] =   0.3942379447222446E+01;
  COFD[        1022] =  -0.3007437798097645E+00;
  COFD[        1023] =   0.1319002761262226E-01;
  COFD[        1024] =  -0.1081045105273500E+02;
  COFD[        1025] =   0.2188985884544625E+01;
  COFD[        1026] =  -0.7096785665772207E-01;
  COFD[        1027] =   0.3131621700706094E-02;
  COFD[        1028] =  -0.1560173804820947E+02;
  COFD[        1029] =   0.4330818915646310E+01;
  COFD[        1030] =  -0.3589811094371836E+00;
  COFD[        1031] =   0.1604703157842258E-01;
  COFD[        1032] =  -0.1340378974537447E+02;
  COFD[        1033] =   0.2956605082389272E+01;
  COFD[        1034] =  -0.1751292924768263E+00;
  COFD[        1035] =   0.7848688796275457E-02;
  COFD[        1036] =  -0.1480659595361166E+02;
  COFD[        1037] =   0.3371224676621789E+01;
  COFD[        1038] =  -0.2300778346811563E+00;
  COFD[        1039] =   0.1027669190040931E-01;
  COFD[        1040] =  -0.1342265746209028E+02;
  COFD[        1041] =   0.2959193932129379E+01;
  COFD[        1042] =  -0.1754915954841434E+00;
  COFD[        1043] =   0.7865516993304794E-02;
  COFD[        1044] =  -0.1896641926147093E+02;
  COFD[        1045] =   0.4932217655116672E+01;
  COFD[        1046] =  -0.4145721175243606E+00;
  COFD[        1047] =   0.1752187325720189E-01;
  COFD[        1048] =  -0.1482502680337964E+02;
  COFD[        1049] =   0.3377285630938066E+01;
  COFD[        1050] =  -0.2309189157634036E+00;
  COFD[        1051] =   0.1031547111935213E-01;
  COFD[        1052] =  -0.1484281133878113E+02;
  COFD[        1053] =   0.3383163225244592E+01;
  COFD[        1054] =  -0.2317345266445065E+00;
  COFD[        1055] =   0.1035307519950314E-01;
  COFD[        1056] =  -0.1328345305948485E+02;
  COFD[        1057] =   0.2871644258385096E+01;
  COFD[        1058] =  -0.1641218343558189E+00;
  COFD[        1059] =   0.7373067026762586E-02;
  COFD[        1060] =  -0.1335130819487631E+02;
  COFD[        1061] =   0.2954148842887173E+01;
  COFD[        1062] =  -0.1747842853483490E+00;
  COFD[        1063] =   0.7832621744754471E-02;
  COFD[        1064] =  -0.1532762961286062E+02;
  COFD[        1065] =   0.3571504907965739E+01;
  COFD[        1066] =  -0.2545911293907537E+00;
  COFD[        1067] =   0.1127475678174522E-01;
  COFD[        1068] =  -0.1532762961286062E+02;
  COFD[        1069] =   0.3571504907965739E+01;
  COFD[        1070] =  -0.2545911293907537E+00;
  COFD[        1071] =   0.1127475678174522E-01;
  COFD[        1072] =  -0.1536331593809694E+02;
  COFD[        1073] =   0.3579933868091604E+01;
  COFD[        1074] =  -0.2557555059112962E+00;
  COFD[        1075] =   0.1132820669357713E-01;
  COFD[        1076] =  -0.1530099655236273E+02;
  COFD[        1077] =   0.3556645688039140E+01;
  COFD[        1078] =  -0.2527119705024404E+00;
  COFD[        1079] =   0.1119535838430676E-01;
  COFD[        1080] =  -0.1451842224029973E+02;
  COFD[        1081] =   0.3244170570593186E+01;
  COFD[        1082] =  -0.2131823036382577E+00;
  COFD[        1083] =   0.9527296769607143E-02;
  COFD[        1084] =  -0.1780941797497641E+02;
  COFD[        1085] =   0.4440466484221220E+01;
  COFD[        1086] =  -0.3655226891172044E+00;
  COFD[        1087] =   0.1600114172157665E-01;
  COFD[        1088] =  -0.1895961739996744E+02;
  COFD[        1089] =   0.4791450137258701E+01;
  COFD[        1090] =  -0.3969919759224265E+00;
  COFD[        1091] =   0.1677814534495883E-01;
  COFD[        1092] =  -0.1896276103155867E+02;
  COFD[        1093] =   0.4789835029851207E+01;
  COFD[        1094] =  -0.3966415331888314E+00;
  COFD[        1095] =   0.1675665490813235E-01;
  COFD[        1096] =  -0.1867823422906515E+02;
  COFD[        1097] =   0.4704907071510658E+01;
  COFD[        1098] =  -0.3897247750072582E+00;
  COFD[        1099] =   0.1662675912848450E-01;
  COFD[        1100] =  -0.1867823422906515E+02;
  COFD[        1101] =   0.4704907071510658E+01;
  COFD[        1102] =  -0.3897247750072582E+00;
  COFD[        1103] =   0.1662675912848450E-01;
  COFD[        1104] =  -0.1890046153942015E+02;
  COFD[        1105] =   0.4763373152612735E+01;
  COFD[        1106] =  -0.3936930155883911E+00;
  COFD[        1107] =   0.1664833359988295E-01;
  COFD[        1108] =  -0.1700832009372265E+02;
  COFD[        1109] =   0.4130721660167143E+01;
  COFD[        1110] =  -0.3260471778937963E+00;
  COFD[        1111] =   0.1432174203443516E-01;
  COFD[        1112] =  -0.1704125657308844E+02;
  COFD[        1113] =   0.4141672233114677E+01;
  COFD[        1114] =  -0.3275236148236659E+00;
  COFD[        1115] =   0.1438795853925400E-01;
  COFD[        1116] =  -0.1707305377045087E+02;
  COFD[        1117] =   0.4152317255109030E+01;
  COFD[        1118] =  -0.3289587183559789E+00;
  COFD[        1119] =   0.1445231595352617E-01;
  COFD[        1120] =  -0.1790606312284444E+02;
  COFD[        1121] =   0.4454152121123291E+01;
  COFD[        1122] =  -0.3650609129312846E+00;
  COFD[        1123] =   0.1588913670146698E-01;
  COFD[        1124] =  -0.1792078479501453E+02;
  COFD[        1125] =   0.4443307443203482E+01;
  COFD[        1126] =  -0.3657004955367218E+00;
  COFD[        1127] =   0.1600230799211683E-01;
  COFD[        1128] =  -0.1794869549063484E+02;
  COFD[        1129] =   0.4452595998578516E+01;
  COFD[        1130] =  -0.3669291177144341E+00;
  COFD[        1131] =   0.1605636779086480E-01;
  COFD[        1132] =  -0.1559790150312219E+02;
  COFD[        1133] =   0.3773944780684847E+01;
  COFD[        1134] =  -0.2821502865661450E+00;
  COFD[        1135] =   0.1252380559259624E-01;
  COFD[        1136] =  -0.1894980246607794E+02;
  COFD[        1137] =   0.4748233808709042E+01;
  COFD[        1138] =  -0.3930779606458064E+00;
  COFD[        1139] =   0.1667238331348419E-01;
  COFD[        1140] =  -0.1894980246607794E+02;
  COFD[        1141] =   0.4748233808709042E+01;
  COFD[        1142] =  -0.3930779606458064E+00;
  COFD[        1143] =   0.1667238331348419E-01;
  COFD[        1144] =  -0.1449021183214253E+02;
  COFD[        1145] =   0.3236447428837034E+01;
  COFD[        1146] =  -0.2121618734582983E+00;
  COFD[        1147] =   0.9482346564680543E-02;
  COFD[        1148] =  -0.1559892478263072E+02;
  COFD[        1149] =   0.3675877348865438E+01;
  COFD[        1150] =  -0.2694954491084075E+00;
  COFD[        1151] =   0.1197817600584641E-01;
  COFD[        1152] =  -0.1092230647098016E+02;
  COFD[        1153] =   0.2306200733264308E+01;
  COFD[        1154] =  -0.8737245270736949E-01;
  COFD[        1155] =   0.3897864029739334E-02;
  COFD[        1156] =  -0.1524698785298783E+02;
  COFD[        1157] =   0.4246332712416798E+01;
  COFD[        1158] =  -0.3448955366037200E+00;
  COFD[        1159] =   0.1530006748241336E-01;
  COFD[        1160] =  -0.1352724646961571E+02;
  COFD[        1161] =   0.3067755508215616E+01;
  COFD[        1162] =  -0.1897498600053720E+00;
  COFD[        1163] =   0.8490405512152714E-02;
  COFD[        1164] =  -0.1516681611308053E+02;
  COFD[        1165] =   0.3569087909107446E+01;
  COFD[        1166] =  -0.2565755708070836E+00;
  COFD[        1167] =   0.1146050720294671E-01;
  COFD[        1168] =  -0.1354935397948868E+02;
  COFD[        1169] =   0.3071473642095726E+01;
  COFD[        1170] =  -0.1902700720448971E+00;
  COFD[        1171] =   0.8514570131842562E-02;
  COFD[        1172] =  -0.1898988022108976E+02;
  COFD[        1173] =   0.4979356265300015E+01;
  COFD[        1174] =  -0.4176235134136435E+00;
  COFD[        1175] =   0.1752883066124088E-01;
  COFD[        1176] =  -0.1519066098667452E+02;
  COFD[        1177] =   0.3577339223376338E+01;
  COFD[        1178] =  -0.2577222156910252E+00;
  COFD[        1179] =   0.1151344824477554E-01;
  COFD[        1180] =  -0.1521372216636978E+02;
  COFD[        1181] =   0.3585353537093755E+01;
  COFD[        1182] =  -0.2588358998442430E+00;
  COFD[        1183] =   0.1156486668966244E-01;
  COFD[        1184] =  -0.1335130819487631E+02;
  COFD[        1185] =   0.2954148842887173E+01;
  COFD[        1186] =  -0.1747842853483490E+00;
  COFD[        1187] =   0.7832621744754471E-02;
  COFD[        1188] =  -0.1346540957089169E+02;
  COFD[        1189] =   0.3062175626208091E+01;
  COFD[        1190] =  -0.1889691489381845E+00;
  COFD[        1191] =   0.8454139853458399E-02;
  COFD[        1192] =  -0.1558104532401472E+02;
  COFD[        1193] =   0.3717715590428343E+01;
  COFD[        1194] =  -0.2731288491005951E+00;
  COFD[        1195] =   0.1205850196546502E-01;
  COFD[        1196] =  -0.1558104532401472E+02;
  COFD[        1197] =   0.3717715590428343E+01;
  COFD[        1198] =  -0.2731288491005951E+00;
  COFD[        1199] =   0.1205850196546502E-01;
  COFD[        1200] =  -0.1562066217821344E+02;
  COFD[        1201] =   0.3727389496240396E+01;
  COFD[        1202] =  -0.2744484704780327E+00;
  COFD[        1203] =   0.1211835981295997E-01;
  COFD[        1204] =  -0.1557193549697668E+02;
  COFD[        1205] =   0.3710292232761665E+01;
  COFD[        1206] =  -0.2722900119431717E+00;
  COFD[        1207] =   0.1202726621500630E-01;
  COFD[        1208] =  -0.1482895398372444E+02;
  COFD[        1209] =   0.3419996364037291E+01;
  COFD[        1210] =  -0.2366451614068416E+00;
  COFD[        1211] =   0.1057153967193523E-01;
  COFD[        1212] =  -0.1797275361165766E+02;
  COFD[        1213] =   0.4534581804401434E+01;
  COFD[        1214] =  -0.3754824832646733E+00;
  COFD[        1215] =   0.1633624099541172E-01;
  COFD[        1216] =  -0.1903795928403629E+02;
  COFD[        1217] =   0.4843815937780696E+01;
  COFD[        1218] =  -0.4001936970428022E+00;
  COFD[        1219] =   0.1676692845698968E-01;
  COFD[        1220] =  -0.1903547544589747E+02;
  COFD[        1221] =   0.4839452409427587E+01;
  COFD[        1222] =  -0.3994241694363118E+00;
  COFD[        1223] =   0.1672455395152373E-01;
  COFD[        1224] =  -0.1871000553285175E+02;
  COFD[        1225] =   0.4738678605504111E+01;
  COFD[        1226] =  -0.3905510969010333E+00;
  COFD[        1227] =   0.1651277300611445E-01;
  COFD[        1228] =  -0.1871000553285175E+02;
  COFD[        1229] =   0.4738678605504111E+01;
  COFD[        1230] =  -0.3905510969010333E+00;
  COFD[        1231] =   0.1651277300611445E-01;
  COFD[        1232] =  -0.1896065537612248E+02;
  COFD[        1233] =   0.4806906431965693E+01;
  COFD[        1234] =  -0.3955804840609516E+00;
  COFD[        1235] =   0.1657242891378949E-01;
  COFD[        1236] =  -0.1735339864186368E+02;
  COFD[        1237] =   0.4304534508536137E+01;
  COFD[        1238] =  -0.3477174523791359E+00;
  COFD[        1239] =   0.1522286508630585E-01;
  COFD[        1240] =  -0.1738771924301953E+02;
  COFD[        1241] =   0.4315709821275419E+01;
  COFD[        1242] =  -0.3491942943333771E+00;
  COFD[        1243] =   0.1528778889586191E-01;
  COFD[        1244] =  -0.1742090070168673E+02;
  COFD[        1245] =   0.4326592022152555E+01;
  COFD[        1246] =  -0.3506322276126634E+00;
  COFD[        1247] =   0.1535099522033430E-01;
  COFD[        1248] =  -0.1802582448642654E+02;
  COFD[        1249] =   0.4528207779823734E+01;
  COFD[        1250] =  -0.3721205839071864E+00;
  COFD[        1251] =   0.1608601965175457E-01;
  COFD[        1252] =  -0.1809136688643801E+02;
  COFD[        1253] =   0.4537201706428067E+01;
  COFD[        1254] =  -0.3756919101622764E+00;
  COFD[        1255] =   0.1634159264026368E-01;
  COFD[        1256] =  -0.1811958823661455E+02;
  COFD[        1257] =   0.4546265606449557E+01;
  COFD[        1258] =  -0.3768550606974422E+00;
  COFD[        1259] =   0.1639117743141738E-01;
  COFD[        1260] =  -0.1578342511116079E+02;
  COFD[        1261] =   0.3901749158409890E+01;
  COFD[        1262] =  -0.2978934623317541E+00;
  COFD[        1263] =   0.1316878577801081E-01;
  COFD[        1264] =  -0.1891456819097078E+02;
  COFD[        1265] =   0.4745553725181649E+01;
  COFD[        1266] =  -0.3882051404698895E+00;
  COFD[        1267] =   0.1626915360804150E-01;
  COFD[        1268] =  -0.1891456819097078E+02;
  COFD[        1269] =   0.4745553725181649E+01;
  COFD[        1270] =  -0.3882051404698895E+00;
  COFD[        1271] =   0.1626915360804150E-01;
  COFD[        1272] =  -0.1479325303072075E+02;
  COFD[        1273] =   0.3409403744671862E+01;
  COFD[        1274] =  -0.2352215543219259E+00;
  COFD[        1275] =   0.1050772674323951E-01;
  COFD[        1276] =  -0.1586061989522203E+02;
  COFD[        1277] =   0.3824799145446967E+01;
  COFD[        1278] =  -0.2881091362509902E+00;
  COFD[        1279] =   0.1275209876189350E-01;
  COFD[        1280] =  -0.1355980676188247E+02;
  COFD[        1281] =   0.3231182295054575E+01;
  COFD[        1282] =  -0.2138406992439825E+00;
  COFD[        1283] =   0.9659145053095021E-02;
  COFD[        1284] =  -0.1798724510881690E+02;
  COFD[        1285] =   0.5079566464178115E+01;
  COFD[        1286] =  -0.4444140659462563E+00;
  COFD[        1287] =   0.1923462920552580E-01;
  COFD[        1288] =  -0.1558834457822415E+02;
  COFD[        1289] =   0.3699567664433178E+01;
  COFD[        1290] =  -0.2706504848826866E+00;
  COFD[        1291] =   0.1194595865045622E-01;
  COFD[        1292] =  -0.1707915522153321E+02;
  COFD[        1293] =   0.4132970234506635E+01;
  COFD[        1294] =  -0.3260547355537204E+00;
  COFD[        1295] =   0.1431106723651207E-01;
  COFD[        1296] =  -0.1559686009787643E+02;
  COFD[        1297] =   0.3697254029246602E+01;
  COFD[        1298] =  -0.2703332731709061E+00;
  COFD[        1299] =   0.1193149792533076E-01;
  COFD[        1300] =  -0.2040913503195636E+02;
  COFD[        1301] =   0.5224092639655129E+01;
  COFD[        1302] =  -0.4327353018434304E+00;
  COFD[        1303] =   0.1753064056853202E-01;
  COFD[        1304] =  -0.1709853282395253E+02;
  COFD[        1305] =   0.4139002691011624E+01;
  COFD[        1306] =  -0.3268662406958736E+00;
  COFD[        1307] =   0.1434738242897968E-01;
  COFD[        1308] =  -0.1711740515182854E+02;
  COFD[        1309] =   0.4144924979490134E+01;
  COFD[        1310] =  -0.3276629227930643E+00;
  COFD[        1311] =   0.1438303438727467E-01;
  COFD[        1312] =  -0.1532762961286062E+02;
  COFD[        1313] =   0.3571504907965739E+01;
  COFD[        1314] =  -0.2545911293907537E+00;
  COFD[        1315] =   0.1127475678174522E-01;
  COFD[        1316] =  -0.1558104532401472E+02;
  COFD[        1317] =   0.3717715590428343E+01;
  COFD[        1318] =  -0.2731288491005951E+00;
  COFD[        1319] =   0.1205850196546502E-01;
  COFD[        1320] =  -0.1769521074851900E+02;
  COFD[        1321] =   0.4367699379901299E+01;
  COFD[        1322] =  -0.3537277912984029E+00;
  COFD[        1323] =   0.1539755051214757E-01;
  COFD[        1324] =  -0.1769521074851900E+02;
  COFD[        1325] =   0.4367699379901299E+01;
  COFD[        1326] =  -0.3537277912984029E+00;
  COFD[        1327] =   0.1539755051214757E-01;
  COFD[        1328] =  -0.1771391298979426E+02;
  COFD[        1329] =   0.4368348686795481E+01;
  COFD[        1330] =  -0.3538107613054961E+00;
  COFD[        1331] =   0.1540107028020538E-01;
  COFD[        1332] =  -0.1765053635794619E+02;
  COFD[        1333] =   0.4344873087030904E+01;
  COFD[        1334] =  -0.3508852749505070E+00;
  COFD[        1335] =   0.1527908501173896E-01;
  COFD[        1336] =  -0.1661775419297465E+02;
  COFD[        1337] =   0.3942050022877351E+01;
  COFD[        1338] =  -0.3009976504232300E+00;
  COFD[        1339] =   0.1321284000123418E-01;
  COFD[        1340] =  -0.1920160126704923E+02;
  COFD[        1341] =   0.4768586004236171E+01;
  COFD[        1342] =  -0.3938478869794211E+00;
  COFD[        1343] =   0.1663350528215831E-01;
  COFD[        1344] =  -0.2039806228291693E+02;
  COFD[        1345] =   0.5097933786196286E+01;
  COFD[        1346] =  -0.4183216752667500E+00;
  COFD[        1347] =   0.1696948938817924E-01;
  COFD[        1348] =  -0.2038049892591359E+02;
  COFD[        1349] =   0.5087331121491625E+01;
  COFD[        1350] =  -0.4167206268593121E+00;
  COFD[        1351] =   0.1689044439564830E-01;
  COFD[        1352] =  -0.2006425198525411E+02;
  COFD[        1353] =   0.4999965267895282E+01;
  COFD[        1354] =  -0.4103327875262441E+00;
  COFD[        1355] =   0.1680582455088142E-01;
  COFD[        1356] =  -0.2006425198525411E+02;
  COFD[        1357] =   0.4999965267895282E+01;
  COFD[        1358] =  -0.4103327875262441E+00;
  COFD[        1359] =   0.1680582455088142E-01;
  COFD[        1360] =  -0.2031516188866531E+02;
  COFD[        1361] =   0.5061917475702875E+01;
  COFD[        1362] =  -0.4141727881359089E+00;
  COFD[        1363] =   0.1680909956107790E-01;
  COFD[        1364] =  -0.1891736997617915E+02;
  COFD[        1365] =   0.4700704222975122E+01;
  COFD[        1366] =  -0.3902105374820116E+00;
  COFD[        1367] =   0.1669545640046241E-01;
  COFD[        1368] =  -0.1892750366101823E+02;
  COFD[        1369] =   0.4701319071264327E+01;
  COFD[        1370] =  -0.3901920115746619E+00;
  COFD[        1371] =   0.1669020608030539E-01;
  COFD[        1372] =  -0.1893716028173049E+02;
  COFD[        1373] =   0.4701920423433434E+01;
  COFD[        1374] =  -0.3901723473212924E+00;
  COFD[        1375] =   0.1668493495004967E-01;
  COFD[        1376] =  -0.1956675647234025E+02;
  COFD[        1377] =   0.4893594117189979E+01;
  COFD[        1378] =  -0.4083224878273328E+00;
  COFD[        1379] =   0.1719480931404233E-01;
  COFD[        1380] =  -0.1941109254354329E+02;
  COFD[        1381] =   0.4819578196751006E+01;
  COFD[        1382] =  -0.4008189828944452E+00;
  COFD[        1383] =   0.1695355424415370E-01;
  COFD[        1384] =  -0.1941401903486105E+02;
  COFD[        1385] =   0.4817662369355102E+01;
  COFD[        1386] =  -0.4004301095900044E+00;
  COFD[        1387] =   0.1693043453375954E-01;
  COFD[        1388] =  -0.1781814842212219E+02;
  COFD[        1389] =   0.4479878049884859E+01;
  COFD[        1390] =  -0.3677558817559772E+00;
  COFD[        1391] =   0.1598132627369455E-01;
  COFD[        1392] =  -0.1997313604118449E+02;
  COFD[        1393] =   0.4888134291131157E+01;
  COFD[        1394] =  -0.3919147663317576E+00;
  COFD[        1395] =   0.1584552869833823E-01;
  COFD[        1396] =  -0.1997313604118449E+02;
  COFD[        1397] =   0.4888134291131157E+01;
  COFD[        1398] =  -0.3919147663317576E+00;
  COFD[        1399] =   0.1584552869833823E-01;
  COFD[        1400] =  -0.1658740614019908E+02;
  COFD[        1401] =   0.3933646317245189E+01;
  COFD[        1402] =  -0.2999344650396615E+00;
  COFD[        1403] =   0.1316797448289397E-01;
  COFD[        1404] =  -0.1779539137381903E+02;
  COFD[        1405] =   0.4385582017594726E+01;
  COFD[        1406] =  -0.3562391667708982E+00;
  COFD[        1407] =   0.1551072662807571E-01;
  COFD[        1408] =  -0.1355980676188247E+02;
  COFD[        1409] =   0.3231182295054575E+01;
  COFD[        1410] =  -0.2138406992439825E+00;
  COFD[        1411] =   0.9659145053095021E-02;
  COFD[        1412] =  -0.1798724510881690E+02;
  COFD[        1413] =   0.5079566464178115E+01;
  COFD[        1414] =  -0.4444140659462563E+00;
  COFD[        1415] =   0.1923462920552580E-01;
  COFD[        1416] =  -0.1558834457822415E+02;
  COFD[        1417] =   0.3699567664433178E+01;
  COFD[        1418] =  -0.2706504848826866E+00;
  COFD[        1419] =   0.1194595865045622E-01;
  COFD[        1420] =  -0.1707915522153321E+02;
  COFD[        1421] =   0.4132970234506635E+01;
  COFD[        1422] =  -0.3260547355537204E+00;
  COFD[        1423] =   0.1431106723651207E-01;
  COFD[        1424] =  -0.1559686009787643E+02;
  COFD[        1425] =   0.3697254029246602E+01;
  COFD[        1426] =  -0.2703332731709061E+00;
  COFD[        1427] =   0.1193149792533076E-01;
  COFD[        1428] =  -0.2040913503195636E+02;
  COFD[        1429] =   0.5224092639655129E+01;
  COFD[        1430] =  -0.4327353018434304E+00;
  COFD[        1431] =   0.1753064056853202E-01;
  COFD[        1432] =  -0.1709853282395253E+02;
  COFD[        1433] =   0.4139002691011624E+01;
  COFD[        1434] =  -0.3268662406958736E+00;
  COFD[        1435] =   0.1434738242897968E-01;
  COFD[        1436] =  -0.1711740515182854E+02;
  COFD[        1437] =   0.4144924979490134E+01;
  COFD[        1438] =  -0.3276629227930643E+00;
  COFD[        1439] =   0.1438303438727467E-01;
  COFD[        1440] =  -0.1532762961286062E+02;
  COFD[        1441] =   0.3571504907965739E+01;
  COFD[        1442] =  -0.2545911293907537E+00;
  COFD[        1443] =   0.1127475678174522E-01;
  COFD[        1444] =  -0.1558104532401472E+02;
  COFD[        1445] =   0.3717715590428343E+01;
  COFD[        1446] =  -0.2731288491005951E+00;
  COFD[        1447] =   0.1205850196546502E-01;
  COFD[        1448] =  -0.1769521074851900E+02;
  COFD[        1449] =   0.4367699379901299E+01;
  COFD[        1450] =  -0.3537277912984029E+00;
  COFD[        1451] =   0.1539755051214757E-01;
  COFD[        1452] =  -0.1769521074851900E+02;
  COFD[        1453] =   0.4367699379901299E+01;
  COFD[        1454] =  -0.3537277912984029E+00;
  COFD[        1455] =   0.1539755051214757E-01;
  COFD[        1456] =  -0.1771391298979426E+02;
  COFD[        1457] =   0.4368348686795481E+01;
  COFD[        1458] =  -0.3538107613054961E+00;
  COFD[        1459] =   0.1540107028020538E-01;
  COFD[        1460] =  -0.1765053635794619E+02;
  COFD[        1461] =   0.4344873087030904E+01;
  COFD[        1462] =  -0.3508852749505070E+00;
  COFD[        1463] =   0.1527908501173896E-01;
  COFD[        1464] =  -0.1661775419297465E+02;
  COFD[        1465] =   0.3942050022877351E+01;
  COFD[        1466] =  -0.3009976504232300E+00;
  COFD[        1467] =   0.1321284000123418E-01;
  COFD[        1468] =  -0.1920160126704923E+02;
  COFD[        1469] =   0.4768586004236171E+01;
  COFD[        1470] =  -0.3938478869794211E+00;
  COFD[        1471] =   0.1663350528215831E-01;
  COFD[        1472] =  -0.2039806228291693E+02;
  COFD[        1473] =   0.5097933786196286E+01;
  COFD[        1474] =  -0.4183216752667500E+00;
  COFD[        1475] =   0.1696948938817924E-01;
  COFD[        1476] =  -0.2038049892591359E+02;
  COFD[        1477] =   0.5087331121491625E+01;
  COFD[        1478] =  -0.4167206268593121E+00;
  COFD[        1479] =   0.1689044439564830E-01;
  COFD[        1480] =  -0.2006425198525411E+02;
  COFD[        1481] =   0.4999965267895282E+01;
  COFD[        1482] =  -0.4103327875262441E+00;
  COFD[        1483] =   0.1680582455088142E-01;
  COFD[        1484] =  -0.2006425198525411E+02;
  COFD[        1485] =   0.4999965267895282E+01;
  COFD[        1486] =  -0.4103327875262441E+00;
  COFD[        1487] =   0.1680582455088142E-01;
  COFD[        1488] =  -0.2031516188866531E+02;
  COFD[        1489] =   0.5061917475702875E+01;
  COFD[        1490] =  -0.4141727881359089E+00;
  COFD[        1491] =   0.1680909956107790E-01;
  COFD[        1492] =  -0.1891736997617915E+02;
  COFD[        1493] =   0.4700704222975122E+01;
  COFD[        1494] =  -0.3902105374820116E+00;
  COFD[        1495] =   0.1669545640046241E-01;
  COFD[        1496] =  -0.1892750366101823E+02;
  COFD[        1497] =   0.4701319071264327E+01;
  COFD[        1498] =  -0.3901920115746619E+00;
  COFD[        1499] =   0.1669020608030539E-01;
  COFD[        1500] =  -0.1893716028173049E+02;
  COFD[        1501] =   0.4701920423433434E+01;
  COFD[        1502] =  -0.3901723473212924E+00;
  COFD[        1503] =   0.1668493495004967E-01;
  COFD[        1504] =  -0.1956675647234025E+02;
  COFD[        1505] =   0.4893594117189979E+01;
  COFD[        1506] =  -0.4083224878273328E+00;
  COFD[        1507] =   0.1719480931404233E-01;
  COFD[        1508] =  -0.1941109254354329E+02;
  COFD[        1509] =   0.4819578196751006E+01;
  COFD[        1510] =  -0.4008189828944452E+00;
  COFD[        1511] =   0.1695355424415370E-01;
  COFD[        1512] =  -0.1941401903486105E+02;
  COFD[        1513] =   0.4817662369355102E+01;
  COFD[        1514] =  -0.4004301095900044E+00;
  COFD[        1515] =   0.1693043453375954E-01;
  COFD[        1516] =  -0.1781814842212219E+02;
  COFD[        1517] =   0.4479878049884859E+01;
  COFD[        1518] =  -0.3677558817559772E+00;
  COFD[        1519] =   0.1598132627369455E-01;
  COFD[        1520] =  -0.1997313604118449E+02;
  COFD[        1521] =   0.4888134291131157E+01;
  COFD[        1522] =  -0.3919147663317576E+00;
  COFD[        1523] =   0.1584552869833823E-01;
  COFD[        1524] =  -0.1997313604118449E+02;
  COFD[        1525] =   0.4888134291131157E+01;
  COFD[        1526] =  -0.3919147663317576E+00;
  COFD[        1527] =   0.1584552869833823E-01;
  COFD[        1528] =  -0.1658740614019908E+02;
  COFD[        1529] =   0.3933646317245189E+01;
  COFD[        1530] =  -0.2999344650396615E+00;
  COFD[        1531] =   0.1316797448289397E-01;
  COFD[        1532] =  -0.1779539137381903E+02;
  COFD[        1533] =   0.4385582017594726E+01;
  COFD[        1534] =  -0.3562391667708982E+00;
  COFD[        1535] =   0.1551072662807571E-01;
  COFD[        1536] =  -0.1358415095584097E+02;
  COFD[        1537] =   0.3239845970896248E+01;
  COFD[        1538] =  -0.2150500406377150E+00;
  COFD[        1539] =   0.9715225271588263E-02;
  COFD[        1540] =  -0.1801636950310948E+02;
  COFD[        1541] =   0.5090088857460223E+01;
  COFD[        1542] =  -0.4457552498826564E+00;
  COFD[        1543] =   0.1929141580461785E-01;
  COFD[        1544] =  -0.1561687122378518E+02;
  COFD[        1545] =   0.3703886367486989E+01;
  COFD[        1546] =  -0.2712410029515650E+00;
  COFD[        1547] =   0.1197280722366742E-01;
  COFD[        1548] =  -0.1707094794267795E+02;
  COFD[        1549] =   0.4119810624341889E+01;
  COFD[        1550] =  -0.3242844616145015E+00;
  COFD[        1551] =   0.1423184749658152E-01;
  COFD[        1552] =  -0.1562210625534589E+02;
  COFD[        1553] =   0.3699988818188416E+01;
  COFD[        1554] =  -0.2707081263160739E+00;
  COFD[        1555] =   0.1194858184298051E-01;
  COFD[        1556] =  -0.2042453573019913E+02;
  COFD[        1557] =   0.5222435968901098E+01;
  COFD[        1558] =  -0.4324987976330752E+00;
  COFD[        1559] =   0.1751958312757915E-01;
  COFD[        1560] =  -0.1709003104221693E+02;
  COFD[        1561] =   0.4125627998285617E+01;
  COFD[        1562] =  -0.3250670331906144E+00;
  COFD[        1563] =   0.1426686739847163E-01;
  COFD[        1564] =  -0.1710868501481421E+02;
  COFD[        1565] =   0.4131370391002061E+01;
  COFD[        1566] =  -0.3258395192522007E+00;
  COFD[        1567] =   0.1430143624455670E-01;
  COFD[        1568] =  -0.1536331593809694E+02;
  COFD[        1569] =   0.3579933868091604E+01;
  COFD[        1570] =  -0.2557555059112962E+00;
  COFD[        1571] =   0.1132820669357713E-01;
  COFD[        1572] =  -0.1562066217821344E+02;
  COFD[        1573] =   0.3727389496240396E+01;
  COFD[        1574] =  -0.2744484704780327E+00;
  COFD[        1575] =   0.1211835981295997E-01;
  COFD[        1576] =  -0.1771391298979426E+02;
  COFD[        1577] =   0.4368348686795481E+01;
  COFD[        1578] =  -0.3538107613054961E+00;
  COFD[        1579] =   0.1540107028020538E-01;
  COFD[        1580] =  -0.1771391298979426E+02;
  COFD[        1581] =   0.4368348686795481E+01;
  COFD[        1582] =  -0.3538107613054961E+00;
  COFD[        1583] =   0.1540107028020538E-01;
  COFD[        1584] =  -0.1772990792056323E+02;
  COFD[        1585] =   0.4367699379901269E+01;
  COFD[        1586] =  -0.3537277912983985E+00;
  COFD[        1587] =   0.1539755051214736E-01;
  COFD[        1588] =  -0.1766459414605863E+02;
  COFD[        1589] =   0.4343242779501765E+01;
  COFD[        1590] =  -0.3506769601296933E+00;
  COFD[        1591] =   0.1527024825917116E-01;
  COFD[        1592] =  -0.1661214673917502E+02;
  COFD[        1593] =   0.3930381940303202E+01;
  COFD[        1594] =  -0.2994310262293675E+00;
  COFD[        1595] =   0.1314286546846328E-01;
  COFD[        1596] =  -0.1923136095279614E+02;
  COFD[        1597] =   0.4771652981887103E+01;
  COFD[        1598] =  -0.3944968044729790E+00;
  COFD[        1599] =   0.1667283377652949E-01;
  COFD[        1600] =  -0.2046646394524606E+02;
  COFD[        1601] =   0.5118752767847491E+01;
  COFD[        1602] =  -0.4214659652131379E+00;
  COFD[        1603] =   0.1712474984073264E-01;
  COFD[        1604] =  -0.2045023591963550E+02;
  COFD[        1605] =   0.5108641161759671E+01;
  COFD[        1606] =  -0.4199387136777197E+00;
  COFD[        1607] =   0.1704933199723114E-01;
  COFD[        1608] =  -0.2012688925964889E+02;
  COFD[        1609] =   0.5018247199767933E+01;
  COFD[        1610] =  -0.4131414445091446E+00;
  COFD[        1611] =   0.1694646037906594E-01;
  COFD[        1612] =  -0.2012688925964889E+02;
  COFD[        1613] =   0.5018247199767933E+01;
  COFD[        1614] =  -0.4131414445091446E+00;
  COFD[        1615] =   0.1694646037906594E-01;
  COFD[        1616] =  -0.2038511261124467E+02;
  COFD[        1617] =   0.5083152894379470E+01;
  COFD[        1618] =  -0.4173863035588838E+00;
  COFD[        1619] =   0.1696799786866059E-01;
  COFD[        1620] =  -0.1893404747219508E+02;
  COFD[        1621] =   0.4699655533259762E+01;
  COFD[        1622] =  -0.3902381208946921E+00;
  COFD[        1623] =   0.1670405772171351E-01;
  COFD[        1624] =  -0.1894432968047615E+02;
  COFD[        1625] =   0.4700246655831890E+01;
  COFD[        1626] =  -0.3902232238988751E+00;
  COFD[        1627] =   0.1669926672613195E-01;
  COFD[        1628] =  -0.1895414515078359E+02;
  COFD[        1629] =   0.4700828894373993E+01;
  COFD[        1630] =  -0.3902069149863984E+00;
  COFD[        1631] =   0.1669440360391662E-01;
  COFD[        1632] =  -0.1960184896363511E+02;
  COFD[        1633] =   0.4900346283733346E+01;
  COFD[        1634] =  -0.4094945328875742E+00;
  COFD[        1635] =   0.1725872606725634E-01;
  COFD[        1636] =  -0.1943925299002717E+02;
  COFD[        1637] =   0.4823356728877616E+01;
  COFD[        1638] =  -0.4015874038077185E+00;
  COFD[        1639] =   0.1699928308073373E-01;
  COFD[        1640] =  -0.1944256474839717E+02;
  COFD[        1641] =   0.4821522129341849E+01;
  COFD[        1642] =  -0.4012140535427826E+00;
  COFD[        1643] =   0.1697705724216803E-01;
  COFD[        1644] =  -0.1782506140103899E+02;
  COFD[        1645] =   0.4472550110816844E+01;
  COFD[        1646] =  -0.3668215949172294E+00;
  COFD[        1647] =   0.1594178999185863E-01;
  COFD[        1648] =  -0.2004560573876227E+02;
  COFD[        1649] =   0.4909854771274889E+01;
  COFD[        1650] =  -0.3952294049921926E+00;
  COFD[        1651] =   0.1601054023046091E-01;
  COFD[        1652] =  -0.2004560573876227E+02;
  COFD[        1653] =   0.4909854771274889E+01;
  COFD[        1654] =  -0.3952294049921926E+00;
  COFD[        1655] =   0.1601054023046091E-01;
  COFD[        1656] =  -0.1658229952875101E+02;
  COFD[        1657] =   0.3922184822536529E+01;
  COFD[        1658] =  -0.2983959485115160E+00;
  COFD[        1659] =   0.1309927190981370E-01;
  COFD[        1660] =  -0.1779710613097560E+02;
  COFD[        1661] =   0.4376269742117769E+01;
  COFD[        1662] =  -0.3550500828076509E+00;
  COFD[        1663] =   0.1546030915508793E-01;
  COFD[        1664] =  -0.1353824314289343E+02;
  COFD[        1665] =   0.3227352939332374E+01;
  COFD[        1666] =  -0.2134245159235849E+00;
  COFD[        1667] =   0.9644648995332943E-02;
  COFD[        1668] =  -0.1790236652645725E+02;
  COFD[        1669] =   0.5054241512059871E+01;
  COFD[        1670] =  -0.4413163862361035E+00;
  COFD[        1671] =   0.1910841720545356E-01;
  COFD[        1672] =  -0.1556217965348612E+02;
  COFD[        1673] =   0.3683628994574988E+01;
  COFD[        1674] =  -0.2686556693855442E+00;
  COFD[        1675] =   0.1186253607235680E-01;
  COFD[        1676] =  -0.1696131964138645E+02;
  COFD[        1677] =   0.4074167708212651E+01;
  COFD[        1678] =  -0.3182576642303228E+00;
  COFD[        1679] =   0.1396612551510442E-01;
  COFD[        1680] =  -0.1556532980903788E+02;
  COFD[        1681] =   0.3678661452138608E+01;
  COFD[        1682] =  -0.2679777825820168E+00;
  COFD[        1683] =   0.1183177506717497E-01;
  COFD[        1684] =  -0.2052802769939090E+02;
  COFD[        1685] =   0.5190809253464677E+01;
  COFD[        1686] =  -0.4197156907449909E+00;
  COFD[        1687] =   0.1661902490158712E-01;
  COFD[        1688] =  -0.1698006191382480E+02;
  COFD[        1689] =   0.4079755891690070E+01;
  COFD[        1690] =  -0.3190093940831370E+00;
  COFD[        1691] =   0.1399976445173436E-01;
  COFD[        1692] =  -0.1699844377766270E+02;
  COFD[        1693] =   0.4085300971938842E+01;
  COFD[        1694] =  -0.3197553288341113E+00;
  COFD[        1695] =   0.1403314440725374E-01;
  COFD[        1696] =  -0.1530099655236273E+02;
  COFD[        1697] =   0.3556645688039140E+01;
  COFD[        1698] =  -0.2527119705024404E+00;
  COFD[        1699] =   0.1119535838430676E-01;
  COFD[        1700] =  -0.1557193549697668E+02;
  COFD[        1701] =   0.3710292232761665E+01;
  COFD[        1702] =  -0.2722900119431717E+00;
  COFD[        1703] =   0.1202726621500630E-01;
  COFD[        1704] =  -0.1765053635794619E+02;
  COFD[        1705] =   0.4344873087030904E+01;
  COFD[        1706] =  -0.3508852749505070E+00;
  COFD[        1707] =   0.1527908501173896E-01;
  COFD[        1708] =  -0.1765053635794619E+02;
  COFD[        1709] =   0.4344873087030904E+01;
  COFD[        1710] =  -0.3508852749505070E+00;
  COFD[        1711] =   0.1527908501173896E-01;
  COFD[        1712] =  -0.1766459414605863E+02;
  COFD[        1713] =   0.4343242779501765E+01;
  COFD[        1714] =  -0.3506769601296933E+00;
  COFD[        1715] =   0.1527024825917116E-01;
  COFD[        1716] =  -0.1760280207090635E+02;
  COFD[        1717] =   0.4320037087778842E+01;
  COFD[        1718] =  -0.3477962090917468E+00;
  COFD[        1719] =   0.1515062790890059E-01;
  COFD[        1720] =  -0.1654504240089281E+02;
  COFD[        1721] =   0.3902990496709335E+01;
  COFD[        1722] =  -0.2959726841289695E+00;
  COFD[        1723] =   0.1299731318228341E-01;
  COFD[        1724] =  -0.1920918070133062E+02;
  COFD[        1725] =   0.4764713252836613E+01;
  COFD[        1726] =  -0.3942405534249802E+00;
  COFD[        1727] =   0.1668923743050566E-01;
  COFD[        1728] =  -0.2048733868393580E+02;
  COFD[        1729] =   0.5132889365248658E+01;
  COFD[        1730] =  -0.4243039026916984E+00;
  COFD[        1731] =   0.1728905581961421E-01;
  COFD[        1732] =  -0.2047278277859408E+02;
  COFD[        1733] =   0.5123414930187060E+01;
  COFD[        1734] =  -0.4228706903648477E+00;
  COFD[        1735] =   0.1721819830097161E-01;
  COFD[        1736] =  -0.2034321636694622E+02;
  COFD[        1737] =   0.5089562300279502E+01;
  COFD[        1738] =  -0.4211891334428361E+00;
  COFD[        1739] =   0.1724909778947470E-01;
  COFD[        1740] =  -0.2034321636694622E+02;
  COFD[        1741] =   0.5089562300279502E+01;
  COFD[        1742] =  -0.4211891334428361E+00;
  COFD[        1743] =   0.1724909778947470E-01;
  COFD[        1744] =  -0.2040021632891986E+02;
  COFD[        1745] =   0.5094681415071051E+01;
  COFD[        1746] =  -0.4198825806270313E+00;
  COFD[        1747] =   0.1711748447383598E-01;
  COFD[        1748] =  -0.1887960780725439E+02;
  COFD[        1749] =   0.4679550122415976E+01;
  COFD[        1750] =  -0.3880189808953476E+00;
  COFD[        1751] =   0.1662436514641822E-01;
  COFD[        1752] =  -0.1889022870587310E+02;
  COFD[        1753] =   0.4680203446785792E+01;
  COFD[        1754] =  -0.3880194721331493E+00;
  COFD[        1755] =   0.1662057070374977E-01;
  COFD[        1756] =  -0.1890039213728353E+02;
  COFD[        1757] =   0.4680852507758593E+01;
  COFD[        1758] =  -0.3880184127527163E+00;
  COFD[        1759] =   0.1661666470021389E-01;
  COFD[        1760] =  -0.1957015754946189E+02;
  COFD[        1761] =   0.4890650507507935E+01;
  COFD[        1762] =  -0.4088565156459115E+00;
  COFD[        1763] =   0.1725720932243744E-01;
  COFD[        1764] =  -0.1938496483487797E+02;
  COFD[        1765] =   0.4803651837482498E+01;
  COFD[        1766] =  -0.3995406103851766E+00;
  COFD[        1767] =   0.1693221082596063E-01;
  COFD[        1768] =  -0.1938900888397246E+02;
  COFD[        1769] =   0.4802051950424429E+01;
  COFD[        1770] =  -0.3992044152028926E+00;
  COFD[        1771] =   0.1691189193038190E-01;
  COFD[        1772] =  -0.1775786527319250E+02;
  COFD[        1773] =   0.4446009317231327E+01;
  COFD[        1774] =  -0.3635789359718016E+00;
  COFD[        1775] =   0.1580980006557884E-01;
  COFD[        1776] =  -0.2008186085747107E+02;
  COFD[        1777] =   0.4928935175864744E+01;
  COFD[        1778] =  -0.3987527366535738E+00;
  COFD[        1779] =   0.1620682050200429E-01;
  COFD[        1780] =  -0.2008186085747107E+02;
  COFD[        1781] =   0.4928935175864744E+01;
  COFD[        1782] =  -0.3987527366535738E+00;
  COFD[        1783] =   0.1620682050200429E-01;
  COFD[        1784] =  -0.1651845590523054E+02;
  COFD[        1785] =   0.3896132878391479E+01;
  COFD[        1786] =  -0.2951170694006329E+00;
  COFD[        1787] =   0.1296170338310411E-01;
  COFD[        1788] =  -0.1772757250170383E+02;
  COFD[        1789] =   0.4347725932428951E+01;
  COFD[        1790] =  -0.3515015798062593E+00;
  COFD[        1791] =   0.1531308129463464E-01;
  COFD[        1792] =  -0.1268296956267843E+02;
  COFD[        1793] =   0.2903428401599904E+01;
  COFD[        1794] =  -0.1707667302851327E+00;
  COFD[        1795] =   0.7771111802295727E-02;
  COFD[        1796] =  -0.1719132602269712E+02;
  COFD[        1797] =   0.4859367175853364E+01;
  COFD[        1798] =  -0.4240580935572132E+00;
  COFD[        1799] =   0.1870696403708276E-01;
  COFD[        1800] =  -0.1478198899238034E+02;
  COFD[        1801] =   0.3371092581023257E+01;
  COFD[        1802] =  -0.2298567856863946E+00;
  COFD[        1803] =   0.1025846107153025E-01;
  COFD[        1804] =  -0.1582244348422915E+02;
  COFD[        1805] =   0.3599211636009877E+01;
  COFD[        1806] =  -0.2581898369286856E+00;
  COFD[        1807] =   0.1143096575439008E-01;
  COFD[        1808] =  -0.1476985219999908E+02;
  COFD[        1809] =   0.3357751196108267E+01;
  COFD[        1810] =  -0.2280046831004087E+00;
  COFD[        1811] =   0.1017303667696902E-01;
  COFD[        1812] =  -0.2031234004461199E+02;
  COFD[        1813] =   0.5186312454417227E+01;
  COFD[        1814] =  -0.4313742876779729E+00;
  COFD[        1815] =   0.1759974844922337E-01;
  COFD[        1816] =  -0.1583284545979817E+02;
  COFD[        1817] =   0.3600569993944359E+01;
  COFD[        1818] =  -0.2583773478248265E+00;
  COFD[        1819] =   0.1143956672273590E-01;
  COFD[        1820] =  -0.1584365035133258E+02;
  COFD[        1821] =   0.3602236757402370E+01;
  COFD[        1822] =  -0.2586073726747737E+00;
  COFD[        1823] =   0.1145011549760728E-01;
  COFD[        1824] =  -0.1451842224029973E+02;
  COFD[        1825] =   0.3244170570593186E+01;
  COFD[        1826] =  -0.2131823036382577E+00;
  COFD[        1827] =   0.9527296769607143E-02;
  COFD[        1828] =  -0.1482895398372444E+02;
  COFD[        1829] =   0.3419996364037291E+01;
  COFD[        1830] =  -0.2366451614068416E+00;
  COFD[        1831] =   0.1057153967193523E-01;
  COFD[        1832] =  -0.1661775419297465E+02;
  COFD[        1833] =   0.3942050022877351E+01;
  COFD[        1834] =  -0.3009976504232300E+00;
  COFD[        1835] =   0.1321284000123418E-01;
  COFD[        1836] =  -0.1661775419297465E+02;
  COFD[        1837] =   0.3942050022877351E+01;
  COFD[        1838] =  -0.3009976504232300E+00;
  COFD[        1839] =   0.1321284000123418E-01;
  COFD[        1840] =  -0.1661214673917502E+02;
  COFD[        1841] =   0.3930381940303202E+01;
  COFD[        1842] =  -0.2994310262293675E+00;
  COFD[        1843] =   0.1314286546846328E-01;
  COFD[        1844] =  -0.1654504240089281E+02;
  COFD[        1845] =   0.3902990496709335E+01;
  COFD[        1846] =  -0.2959726841289695E+00;
  COFD[        1847] =   0.1299731318228341E-01;
  COFD[        1848] =  -0.1551400389810888E+02;
  COFD[        1849] =   0.3472777644794381E+01;
  COFD[        1850] =  -0.2416432013010369E+00;
  COFD[        1851] =   0.1070799930256162E-01;
  COFD[        1852] =  -0.1853079633836986E+02;
  COFD[        1853] =   0.4516962151189352E+01;
  COFD[        1854] =  -0.3708427592293742E+00;
  COFD[        1855] =   0.1604368344993687E-01;
  COFD[        1856] =  -0.2009021453031469E+02;
  COFD[        1857] =   0.5037634138971256E+01;
  COFD[        1858] =  -0.4241196608853510E+00;
  COFD[        1859] =   0.1777020759153531E-01;
  COFD[        1860] =  -0.2009778046496907E+02;
  COFD[        1861] =   0.5037181430781024E+01;
  COFD[        1862] =  -0.4240362270814515E+00;
  COFD[        1863] =   0.1776546650877937E-01;
  COFD[        1864] =  -0.2003301781547066E+02;
  COFD[        1865] =   0.5028625403684593E+01;
  COFD[        1866] =  -0.4259384280972174E+00;
  COFD[        1867] =   0.1797007932483027E-01;
  COFD[        1868] =  -0.2003301781547066E+02;
  COFD[        1869] =   0.5028625403684593E+01;
  COFD[        1870] =  -0.4259384280972174E+00;
  COFD[        1871] =   0.1797007932483027E-01;
  COFD[        1872] =  -0.2002461546597822E+02;
  COFD[        1873] =   0.5005680869193603E+01;
  COFD[        1874] =  -0.4206313238956306E+00;
  COFD[        1875] =   0.1764515597235788E-01;
  COFD[        1876] =  -0.1800841736946497E+02;
  COFD[        1877] =   0.4351921120886195E+01;
  COFD[        1878] =  -0.3517363610689986E+00;
  COFD[        1879] =   0.1531343604033216E-01;
  COFD[        1880] =  -0.1801976002924096E+02;
  COFD[        1881] =   0.4352352585107287E+01;
  COFD[        1882] =  -0.3517965558706893E+00;
  COFD[        1883] =   0.1531621848405548E-01;
  COFD[        1884] =  -0.1803148637733784E+02;
  COFD[        1885] =   0.4353159047512093E+01;
  COFD[        1886] =  -0.3519044208218940E+00;
  COFD[        1887] =   0.1532101211937375E-01;
  COFD[        1888] =  -0.1879103964706305E+02;
  COFD[        1889] =   0.4613167070471964E+01;
  COFD[        1890] =  -0.3811821259847005E+00;
  COFD[        1891] =   0.1640353922053438E-01;
  COFD[        1892] =  -0.1864571977877140E+02;
  COFD[        1893] =   0.4539069966095800E+01;
  COFD[        1894] =  -0.3735011850255855E+00;
  COFD[        1895] =   0.1615133908644331E-01;
  COFD[        1896] =  -0.1865835471859333E+02;
  COFD[        1897] =   0.4540734987509840E+01;
  COFD[        1898] =  -0.3737078963922989E+00;
  COFD[        1899] =   0.1615982140755382E-01;
  COFD[        1900] =  -0.1658530159401671E+02;
  COFD[        1901] =   0.3962484009621037E+01;
  COFD[        1902] =  -0.3034260162001460E+00;
  COFD[        1903] =   0.1330924504181773E-01;
  COFD[        1904] =  -0.1996446349112875E+02;
  COFD[        1905] =   0.4943053687420567E+01;
  COFD[        1906] =  -0.4143586003201034E+00;
  COFD[        1907] =   0.1744554338516522E-01;
  COFD[        1908] =  -0.1996446349112875E+02;
  COFD[        1909] =   0.4943053687420567E+01;
  COFD[        1910] =  -0.4143586003201034E+00;
  COFD[        1911] =   0.1744554338516522E-01;
  COFD[        1912] =  -0.1549017806948574E+02;
  COFD[        1913] =   0.3466859355125401E+01;
  COFD[        1914] =  -0.2408856054712291E+00;
  COFD[        1915] =   0.1067561900190753E-01;
  COFD[        1916] =  -0.1650781936547825E+02;
  COFD[        1917] =   0.3842122413055806E+01;
  COFD[        1918] =  -0.2882098997185555E+00;
  COFD[        1919] =   0.1266705263344362E-01;
  COFD[        1920] =  -0.1564739122681956E+02;
  COFD[        1921] =   0.4025825276486324E+01;
  COFD[        1922] =  -0.3181195147607377E+00;
  COFD[        1923] =   0.1422136468533644E-01;
  COFD[        1924] =  -0.1725034548774479E+02;
  COFD[        1925] =   0.4566921339972978E+01;
  COFD[        1926] =  -0.3509192271095338E+00;
  COFD[        1927] =   0.1402685086906688E-01;
  COFD[        1928] =  -0.1793673861488134E+02;
  COFD[        1929] =   0.4489875914439151E+01;
  COFD[        1930] =  -0.3697556337020395E+00;
  COFD[        1931] =   0.1609253879718407E-01;
  COFD[        1932] =  -0.1868514242925209E+02;
  COFD[        1933] =   0.4565512377700616E+01;
  COFD[        1934] =  -0.3754538805580755E+00;
  COFD[        1935] =   0.1617305773699574E-01;
  COFD[        1936] =  -0.1792454027766008E+02;
  COFD[        1937] =   0.4476326984714260E+01;
  COFD[        1938] =  -0.3680191278135618E+00;
  COFD[        1939] =   0.1601860672433683E-01;
  COFD[        1940] =  -0.1924815690645235E+02;
  COFD[        1941] =   0.4356731855020855E+01;
  COFD[        1942] =  -0.2835834611821879E+00;
  COFD[        1943] =   0.9647123135613308E-02;
  COFD[        1944] =  -0.1868972720602865E+02;
  COFD[        1945] =   0.4564003101584411E+01;
  COFD[        1946] =  -0.3752914154261894E+00;
  COFD[        1947] =   0.1616756912994680E-01;
  COFD[        1948] =  -0.1869430855359004E+02;
  COFD[        1949] =   0.4562621031385030E+01;
  COFD[        1950] =  -0.3751419042160075E+00;
  COFD[        1951] =   0.1616247820310002E-01;
  COFD[        1952] =  -0.1780941797497641E+02;
  COFD[        1953] =   0.4440466484221220E+01;
  COFD[        1954] =  -0.3655226891172044E+00;
  COFD[        1955] =   0.1600114172157665E-01;
  COFD[        1956] =  -0.1797275361165766E+02;
  COFD[        1957] =   0.4534581804401434E+01;
  COFD[        1958] =  -0.3754824832646733E+00;
  COFD[        1959] =   0.1633624099541172E-01;
  COFD[        1960] =  -0.1920160126704923E+02;
  COFD[        1961] =   0.4768586004236171E+01;
  COFD[        1962] =  -0.3938478869794211E+00;
  COFD[        1963] =   0.1663350528215831E-01;
  COFD[        1964] =  -0.1920160126704923E+02;
  COFD[        1965] =   0.4768586004236171E+01;
  COFD[        1966] =  -0.3938478869794211E+00;
  COFD[        1967] =   0.1663350528215831E-01;
  COFD[        1968] =  -0.1923136095279614E+02;
  COFD[        1969] =   0.4771652981887103E+01;
  COFD[        1970] =  -0.3944968044729790E+00;
  COFD[        1971] =   0.1667283377652949E-01;
  COFD[        1972] =  -0.1920918070133062E+02;
  COFD[        1973] =   0.4764713252836613E+01;
  COFD[        1974] =  -0.3942405534249802E+00;
  COFD[        1975] =   0.1668923743050566E-01;
  COFD[        1976] =  -0.1853079633836986E+02;
  COFD[        1977] =   0.4516962151189352E+01;
  COFD[        1978] =  -0.3708427592293742E+00;
  COFD[        1979] =   0.1604368344993687E-01;
  COFD[        1980] =  -0.2074682089064003E+02;
  COFD[        1981] =   0.5125913401864712E+01;
  COFD[        1982] =  -0.4300036496008145E+00;
  COFD[        1983] =   0.1780187794631543E-01;
  COFD[        1984] =  -0.2103480542370173E+02;
  COFD[        1985] =   0.5058663047911905E+01;
  COFD[        1986] =  -0.3955876986371539E+00;
  COFD[        1987] =   0.1531365269600822E-01;
  COFD[        1988] =  -0.2106203504941583E+02;
  COFD[        1989] =   0.5066201505062354E+01;
  COFD[        1990] =  -0.3966848046994042E+00;
  COFD[        1991] =   0.1536587788350550E-01;
  COFD[        1992] =  -0.2111193727907456E+02;
  COFD[        1993] =   0.5115004803331022E+01;
  COFD[        1994] =  -0.4073081193984433E+00;
  COFD[        1995] =   0.1598596435486438E-01;
  COFD[        1996] =  -0.2111193727907456E+02;
  COFD[        1997] =   0.5115004803331022E+01;
  COFD[        1998] =  -0.4073081193984433E+00;
  COFD[        1999] =   0.1598596435486438E-01;
  COFD[        2000] =  -0.2111512317179221E+02;
  COFD[        2001] =   0.5089049253783277E+01;
  COFD[        2002] =  -0.4010600389026057E+00;
  COFD[        2003] =   0.1560812400748400E-01;
  COFD[        2004] =  -0.2036669718093571E+02;
  COFD[        2005] =   0.5042701296018135E+01;
  COFD[        2006] =  -0.4228382373074129E+00;
  COFD[        2007] =   0.1763442522229875E-01;
  COFD[        2008] =  -0.2039137163285755E+02;
  COFD[        2009] =   0.5048551827610241E+01;
  COFD[        2010] =  -0.4237554975633671E+00;
  COFD[        2011] =   0.1768105924766743E-01;
  COFD[        2012] =  -0.2041424604605752E+02;
  COFD[        2013] =   0.5053848971414935E+01;
  COFD[        2014] =  -0.4245860964146182E+00;
  COFD[        2015] =   0.1772329162485299E-01;
  COFD[        2016] =  -0.2089328018920450E+02;
  COFD[        2017] =   0.5166930585415199E+01;
  COFD[        2018] =  -0.4306289589876920E+00;
  COFD[        2019] =   0.1764267219152353E-01;
  COFD[        2020] =  -0.2075862463488640E+02;
  COFD[        2021] =   0.5108417296497037E+01;
  COFD[        2022] =  -0.4262851069696628E+00;
  COFD[        2023] =   0.1758280414396813E-01;
  COFD[        2024] =  -0.2077840705218999E+02;
  COFD[        2025] =   0.5112879719806333E+01;
  COFD[        2026] =  -0.4269664397363305E+00;
  COFD[        2027] =   0.1761673548420319E-01;
  COFD[        2028] =  -0.1934676521684402E+02;
  COFD[        2029] =   0.4844117217038898E+01;
  COFD[        2030] =  -0.4050654922829996E+00;
  COFD[        2031] =   0.1718562636894347E-01;
  COFD[        2032] =  -0.2140522217493100E+02;
  COFD[        2033] =   0.5181866929752220E+01;
  COFD[        2034] =  -0.4179647553808900E+00;
  COFD[        2035] =   0.1653089926470964E-01;
  COFD[        2036] =  -0.2140522217493100E+02;
  COFD[        2037] =   0.5181866929752220E+01;
  COFD[        2038] =  -0.4179647553808900E+00;
  COFD[        2039] =   0.1653089926470964E-01;
  COFD[        2040] =  -0.1850244542172855E+02;
  COFD[        2041] =   0.4509546506902997E+01;
  COFD[        2042] =  -0.3699250247632380E+00;
  COFD[        2043] =   0.1600569104209367E-01;
  COFD[        2044] =  -0.1931016351705685E+02;
  COFD[        2045] =   0.4757139666603682E+01;
  COFD[        2046] =  -0.3961700036764624E+00;
  COFD[        2047] =   0.1690081368592412E-01;
  COFD[        2048] =  -0.1746999275651299E+02;
  COFD[        2049] =   0.4670556529920976E+01;
  COFD[        2050] =  -0.3931788949833612E+00;
  COFD[        2051] =   0.1709979614643272E-01;
  COFD[        2052] =  -0.1540645026491837E+02;
  COFD[        2053] =   0.3572427121876953E+01;
  COFD[        2054] =  -0.1878763535992870E+00;
  COFD[        2055] =   0.5591984032298197E-02;
  COFD[        2056] =  -0.1915075646861156E+02;
  COFD[        2057] =   0.4868838487330766E+01;
  COFD[        2058] =  -0.4046210142562804E+00;
  COFD[        2059] =   0.1701122981099066E-01;
  COFD[        2060] =  -0.2031657960333580E+02;
  COFD[        2061] =   0.5113614640951361E+01;
  COFD[        2062] =  -0.4319293017089550E+00;
  COFD[        2063] =   0.1802805482147992E-01;
  COFD[        2064] =  -0.1918204714711579E+02;
  COFD[        2065] =   0.4875615818554563E+01;
  COFD[        2066] =  -0.4058262469897604E+00;
  COFD[        2067] =   0.1707794755285167E-01;
  COFD[        2068] =  -0.1757312087019186E+02;
  COFD[        2069] =   0.3511204681563826E+01;
  COFD[        2070] =  -0.1497966127823777E+00;
  COFD[        2071] =   0.2966580974304371E-02;
  COFD[        2072] =  -0.2032173485009814E+02;
  COFD[        2073] =   0.5112656067238100E+01;
  COFD[        2074] =  -0.4317894593025599E+00;
  COFD[        2075] =   0.1802136540679988E-01;
  COFD[        2076] =  -0.2032584895987123E+02;
  COFD[        2077] =   0.5111360106328648E+01;
  COFD[        2078] =  -0.4315969740328803E+00;
  COFD[        2079] =   0.1801200837165252E-01;
  COFD[        2080] =  -0.1895961739996744E+02;
  COFD[        2081] =   0.4791450137258701E+01;
  COFD[        2082] =  -0.3969919759224265E+00;
  COFD[        2083] =   0.1677814534495883E-01;
  COFD[        2084] =  -0.1903795928403629E+02;
  COFD[        2085] =   0.4843815937780696E+01;
  COFD[        2086] =  -0.4001936970428022E+00;
  COFD[        2087] =   0.1676692845698968E-01;
  COFD[        2088] =  -0.2039806228291693E+02;
  COFD[        2089] =   0.5097933786196286E+01;
  COFD[        2090] =  -0.4183216752667500E+00;
  COFD[        2091] =   0.1696948938817924E-01;
  COFD[        2092] =  -0.2039806228291693E+02;
  COFD[        2093] =   0.5097933786196286E+01;
  COFD[        2094] =  -0.4183216752667500E+00;
  COFD[        2095] =   0.1696948938817924E-01;
  COFD[        2096] =  -0.2046646394524606E+02;
  COFD[        2097] =   0.5118752767847491E+01;
  COFD[        2098] =  -0.4214659652131379E+00;
  COFD[        2099] =   0.1712474984073264E-01;
  COFD[        2100] =  -0.2048733868393580E+02;
  COFD[        2101] =   0.5132889365248658E+01;
  COFD[        2102] =  -0.4243039026916984E+00;
  COFD[        2103] =   0.1728905581961421E-01;
  COFD[        2104] =  -0.2009021453031469E+02;
  COFD[        2105] =   0.5037634138971256E+01;
  COFD[        2106] =  -0.4241196608853510E+00;
  COFD[        2107] =   0.1777020759153531E-01;
  COFD[        2108] =  -0.2103480542370173E+02;
  COFD[        2109] =   0.5058663047911905E+01;
  COFD[        2110] =  -0.3955876986371539E+00;
  COFD[        2111] =   0.1531365269600822E-01;
  COFD[        2112] =  -0.1885553227873480E+02;
  COFD[        2113] =   0.3923352162852614E+01;
  COFD[        2114] =  -0.2115865709699218E+00;
  COFD[        2115] =   0.5939285843445688E-02;
  COFD[        2116] =  -0.1886381013480732E+02;
  COFD[        2117] =   0.3923279625734948E+01;
  COFD[        2118] =  -0.2115777534001486E+00;
  COFD[        2119] =   0.5938970665371634E-02;
  COFD[        2120] =  -0.1981646458797120E+02;
  COFD[        2121] =   0.4378886175261416E+01;
  COFD[        2122] =  -0.2813828834600260E+00;
  COFD[        2123] =   0.9370857616443682E-02;
  COFD[        2124] =  -0.1981646458797120E+02;
  COFD[        2125] =   0.4378886175261416E+01;
  COFD[        2126] =  -0.2813828834600260E+00;
  COFD[        2127] =   0.9370857616443682E-02;
  COFD[        2128] =  -0.1910884672535752E+02;
  COFD[        2129] =   0.4028503489934288E+01;
  COFD[        2130] =  -0.2273989668967882E+00;
  COFD[        2131] =   0.6707295013470504E-02;
  COFD[        2132] =  -0.2111949107097320E+02;
  COFD[        2133] =   0.5171726543065573E+01;
  COFD[        2134] =  -0.4171067949473953E+00;
  COFD[        2135] =   0.1650848848746980E-01;
  COFD[        2136] =  -0.2113317649345642E+02;
  COFD[        2137] =   0.5173175010557953E+01;
  COFD[        2138] =  -0.4173234769896741E+00;
  COFD[        2139] =   0.1651910472359747E-01;
  COFD[        2140] =  -0.2114401307743194E+02;
  COFD[        2141] =   0.5173607724136159E+01;
  COFD[        2142] =  -0.4173904988942925E+00;
  COFD[        2143] =   0.1652249796434591E-01;
  COFD[        2144] =  -0.2090933405715619E+02;
  COFD[        2145] =   0.4972999527197874E+01;
  COFD[        2146] =  -0.3786889263628220E+00;
  COFD[        2147] =   0.1435806171655258E-01;
  COFD[        2148] =  -0.2110606940488274E+02;
  COFD[        2149] =   0.5056386072237322E+01;
  COFD[        2150] =  -0.3941314241059736E+00;
  COFD[        2151] =   0.1520386988001248E-01;
  COFD[        2152] =  -0.2110836259086301E+02;
  COFD[        2153] =   0.5053674020415143E+01;
  COFD[        2154] =  -0.3937388243323189E+00;
  COFD[        2155] =   0.1518528106281032E-01;
  COFD[        2156] =  -0.2062895610140275E+02;
  COFD[        2157] =   0.5231225460685535E+01;
  COFD[        2158] =  -0.4370160277895190E+00;
  COFD[        2159] =   0.1784652320332632E-01;
  COFD[        2160] =  -0.1972933562163679E+02;
  COFD[        2161] =   0.4264735345000579E+01;
  COFD[        2162] =  -0.2639983274633736E+00;
  COFD[        2163] =   0.8526026392458589E-02;
  COFD[        2164] =  -0.1972933562163679E+02;
  COFD[        2165] =   0.4264735345000579E+01;
  COFD[        2166] =  -0.2639983274633736E+00;
  COFD[        2167] =   0.8526026392458589E-02;
  COFD[        2168] =  -0.2006577079410498E+02;
  COFD[        2169] =   0.5032536447034316E+01;
  COFD[        2170] =  -0.4235925902580250E+00;
  COFD[        2171] =   0.1775279409683131E-01;
  COFD[        2172] =  -0.2070496369832456E+02;
  COFD[        2173] =   0.5191072368701859E+01;
  COFD[        2174] =  -0.4346910627906122E+00;
  COFD[        2175] =   0.1785842811134005E-01;
  COFD[        2176] =  -0.1747771970561251E+02;
  COFD[        2177] =   0.4673159679696371E+01;
  COFD[        2178] =  -0.3935109715079382E+00;
  COFD[        2179] =   0.1711386921544412E-01;
  COFD[        2180] =  -0.1539084008249366E+02;
  COFD[        2181] =   0.3565006575777335E+01;
  COFD[        2182] =  -0.1867580729432450E+00;
  COFD[        2183] =   0.5536884066371989E-02;
  COFD[        2184] =  -0.1914986770399175E+02;
  COFD[        2185] =   0.4864899621234062E+01;
  COFD[        2186] =  -0.4039220297382060E+00;
  COFD[        2187] =   0.1697258735720344E-01;
  COFD[        2188] =  -0.2032689236579443E+02;
  COFD[        2189] =   0.5114221260960875E+01;
  COFD[        2190] =  -0.4320136109530544E+00;
  COFD[        2191] =   0.1803190494311217E-01;
  COFD[        2192] =  -0.1918177070837826E+02;
  COFD[        2193] =   0.4871873199573439E+01;
  COFD[        2194] =  -0.4051602536432213E+00;
  COFD[        2195] =   0.1704106539271711E-01;
  COFD[        2196] =  -0.1757470287826831E+02;
  COFD[        2197] =   0.3509436887842453E+01;
  COFD[        2198] =  -0.1496102927025162E+00;
  COFD[        2199] =   0.2961974572337444E-02;
  COFD[        2200] =  -0.2033308000404864E+02;
  COFD[        2201] =   0.5113689932626547E+01;
  COFD[        2202] =  -0.4319400850127593E+00;
  COFD[        2203] =   0.1802856188319298E-01;
  COFD[        2204] =  -0.2033818516187319E+02;
  COFD[        2205] =   0.5112803625509553E+01;
  COFD[        2206] =  -0.4318111908846191E+00;
  COFD[        2207] =   0.1802241389102129E-01;
  COFD[        2208] =  -0.1896276103155867E+02;
  COFD[        2209] =   0.4789835029851207E+01;
  COFD[        2210] =  -0.3966415331888314E+00;
  COFD[        2211] =   0.1675665490813235E-01;
  COFD[        2212] =  -0.1903547544589747E+02;
  COFD[        2213] =   0.4839452409427587E+01;
  COFD[        2214] =  -0.3994241694363118E+00;
  COFD[        2215] =   0.1672455395152373E-01;
  COFD[        2216] =  -0.2038049892591359E+02;
  COFD[        2217] =   0.5087331121491625E+01;
  COFD[        2218] =  -0.4167206268593121E+00;
  COFD[        2219] =   0.1689044439564830E-01;
  COFD[        2220] =  -0.2038049892591359E+02;
  COFD[        2221] =   0.5087331121491625E+01;
  COFD[        2222] =  -0.4167206268593121E+00;
  COFD[        2223] =   0.1689044439564830E-01;
  COFD[        2224] =  -0.2045023591963550E+02;
  COFD[        2225] =   0.5108641161759671E+01;
  COFD[        2226] =  -0.4199387136777197E+00;
  COFD[        2227] =   0.1704933199723114E-01;
  COFD[        2228] =  -0.2047278277859408E+02;
  COFD[        2229] =   0.5123414930187060E+01;
  COFD[        2230] =  -0.4228706903648477E+00;
  COFD[        2231] =   0.1721819830097161E-01;
  COFD[        2232] =  -0.2009778046496907E+02;
  COFD[        2233] =   0.5037181430781024E+01;
  COFD[        2234] =  -0.4240362270814515E+00;
  COFD[        2235] =   0.1776546650877937E-01;
  COFD[        2236] =  -0.2106203504941583E+02;
  COFD[        2237] =   0.5066201505062354E+01;
  COFD[        2238] =  -0.3966848046994042E+00;
  COFD[        2239] =   0.1536587788350550E-01;
  COFD[        2240] =  -0.1886381013480732E+02;
  COFD[        2241] =   0.3923279625734948E+01;
  COFD[        2242] =  -0.2115777534001486E+00;
  COFD[        2243] =   0.5938970665371634E-02;
  COFD[        2244] =  -0.1887260515065304E+02;
  COFD[        2245] =   0.3923352162852601E+01;
  COFD[        2246] =  -0.2115865709699200E+00;
  COFD[        2247] =   0.5939285843445600E-02;
  COFD[        2248] =  -0.1982650591715289E+02;
  COFD[        2249] =   0.4379424471801094E+01;
  COFD[        2250] =  -0.2814553408580349E+00;
  COFD[        2251] =   0.9373962429419786E-02;
  COFD[        2252] =  -0.1982650591715289E+02;
  COFD[        2253] =   0.4379424471801094E+01;
  COFD[        2254] =  -0.2814553408580349E+00;
  COFD[        2255] =   0.9373962429419786E-02;
  COFD[        2256] =  -0.1911879294336443E+02;
  COFD[        2257] =   0.4028925322234755E+01;
  COFD[        2258] =  -0.2274518266873066E+00;
  COFD[        2259] =   0.6709297609983591E-02;
  COFD[        2260] =  -0.2112276604025784E+02;
  COFD[        2261] =   0.5169665446870568E+01;
  COFD[        2262] =  -0.4168003974284392E+00;
  COFD[        2263] =   0.1649356899195412E-01;
  COFD[        2264] =  -0.2113855609134408E+02;
  COFD[        2265] =   0.5171981200778053E+01;
  COFD[        2266] =  -0.4171447495559704E+00;
  COFD[        2267] =   0.1651034132923998E-01;
  COFD[        2268] =  -0.2115144768087362E+02;
  COFD[        2269] =   0.5173261676775666E+01;
  COFD[        2270] =  -0.4173365706093256E+00;
  COFD[        2271] =   0.1651975240232191E-01;
  COFD[        2272] =  -0.2091775085800750E+02;
  COFD[        2273] =   0.4973029974129495E+01;
  COFD[        2274] =  -0.3786913370630364E+00;
  COFD[        2275] =   0.1435807083846824E-01;
  COFD[        2276] =  -0.2111907495038380E+02;
  COFD[        2277] =   0.5058381671032987E+01;
  COFD[        2278] =  -0.3944196800959486E+00;
  COFD[        2279] =   0.1521748580862574E-01;
  COFD[        2280] =  -0.2112313758034541E+02;
  COFD[        2281] =   0.5056389549388699E+01;
  COFD[        2282] =  -0.3941319269977342E+00;
  COFD[        2283] =   0.1520389366699413E-01;
  COFD[        2284] =  -0.2063989479595861E+02;
  COFD[        2285] =   0.5231645585324166E+01;
  COFD[        2286] =  -0.4370732404462373E+00;
  COFD[        2287] =   0.1784907337224841E-01;
  COFD[        2288] =  -0.1974668857323474E+02;
  COFD[        2289] =   0.4267776444888849E+01;
  COFD[        2290] =  -0.2644050809170580E+00;
  COFD[        2291] =   0.8543267577410659E-02;
  COFD[        2292] =  -0.1974668857323474E+02;
  COFD[        2293] =   0.4267776444888849E+01;
  COFD[        2294] =  -0.2644050809170580E+00;
  COFD[        2295] =   0.8543267577410659E-02;
  COFD[        2296] =  -0.2007323696432831E+02;
  COFD[        2297] =   0.5032033200707714E+01;
  COFD[        2298] =  -0.4235009227577490E+00;
  COFD[        2299] =   0.1774762291235686E-01;
  COFD[        2300] =  -0.2072332061257778E+02;
  COFD[        2301] =   0.5194960177957292E+01;
  COFD[        2302] =  -0.4352726149968372E+00;
  COFD[        2303] =   0.1788688372961993E-01;
  COFD[        2304] =  -0.1797854335218850E+02;
  COFD[        2305] =   0.4905020138101169E+01;
  COFD[        2306] =  -0.4293084201751438E+00;
  COFD[        2307] =   0.1891302999394448E-01;
  COFD[        2308] =  -0.1561125740910447E+02;
  COFD[        2309] =   0.3698144089263720E+01;
  COFD[        2310] =  -0.2103080207621357E+00;
  COFD[        2311] =   0.6792204226221277E-02;
  COFD[        2312] =  -0.1878355293596396E+02;
  COFD[        2313] =   0.4745970603113269E+01;
  COFD[        2314] =  -0.3924372163577848E+00;
  COFD[        2315] =   0.1663659729501466E-01;
  COFD[        2316] =  -0.2016013281728147E+02;
  COFD[        2317] =   0.5063566265537697E+01;
  COFD[        2318] =  -0.4281515792120217E+00;
  COFD[        2319] =   0.1797180686456706E-01;
  COFD[        2320] =  -0.1880409131690931E+02;
  COFD[        2321] =   0.4747896766189557E+01;
  COFD[        2322] =  -0.3929513558466455E+00;
  COFD[        2323] =   0.1667071927589917E-01;
  COFD[        2324] =  -0.1619587514529141E+02;
  COFD[        2325] =   0.2843235032975332E+01;
  COFD[        2326] =  -0.5104314541884596E-01;
  COFD[        2327] =  -0.1707937587486378E-02;
  COFD[        2328] =  -0.1994316533928780E+02;
  COFD[        2329] =   0.4989293390342089E+01;
  COFD[        2330] =  -0.4200959859596722E+00;
  COFD[        2331] =   0.1768643162186107E-01;
  COFD[        2332] =  -0.1994945707378163E+02;
  COFD[        2333] =   0.4988909811915903E+01;
  COFD[        2334] =  -0.4200448560138315E+00;
  COFD[        2335] =   0.1768418875431581E-01;
  COFD[        2336] =  -0.1867823422906515E+02;
  COFD[        2337] =   0.4704907071510658E+01;
  COFD[        2338] =  -0.3897247750072582E+00;
  COFD[        2339] =   0.1662675912848450E-01;
  COFD[        2340] =  -0.1871000553285175E+02;
  COFD[        2341] =   0.4738678605504111E+01;
  COFD[        2342] =  -0.3905510969010333E+00;
  COFD[        2343] =   0.1651277300611445E-01;
  COFD[        2344] =  -0.2006425198525411E+02;
  COFD[        2345] =   0.4999965267895282E+01;
  COFD[        2346] =  -0.4103327875262441E+00;
  COFD[        2347] =   0.1680582455088142E-01;
  COFD[        2348] =  -0.2006425198525411E+02;
  COFD[        2349] =   0.4999965267895282E+01;
  COFD[        2350] =  -0.4103327875262441E+00;
  COFD[        2351] =   0.1680582455088142E-01;
  COFD[        2352] =  -0.2012688925964889E+02;
  COFD[        2353] =   0.5018247199767933E+01;
  COFD[        2354] =  -0.4131414445091446E+00;
  COFD[        2355] =   0.1694646037906594E-01;
  COFD[        2356] =  -0.2034321636694622E+02;
  COFD[        2357] =   0.5089562300279502E+01;
  COFD[        2358] =  -0.4211891334428361E+00;
  COFD[        2359] =   0.1724909778947470E-01;
  COFD[        2360] =  -0.2003301781547066E+02;
  COFD[        2361] =   0.5028625403684593E+01;
  COFD[        2362] =  -0.4259384280972174E+00;
  COFD[        2363] =   0.1797007932483027E-01;
  COFD[        2364] =  -0.2111193727907456E+02;
  COFD[        2365] =   0.5115004803331022E+01;
  COFD[        2366] =  -0.4073081193984433E+00;
  COFD[        2367] =   0.1598596435486438E-01;
  COFD[        2368] =  -0.1981646458797120E+02;
  COFD[        2369] =   0.4378886175261416E+01;
  COFD[        2370] =  -0.2813828834600260E+00;
  COFD[        2371] =   0.9370857616443682E-02;
  COFD[        2372] =  -0.1982650591715289E+02;
  COFD[        2373] =   0.4379424471801094E+01;
  COFD[        2374] =  -0.2814553408580349E+00;
  COFD[        2375] =   0.9373962429419786E-02;
  COFD[        2376] =  -0.1906244143057827E+02;
  COFD[        2377] =   0.4049325847301449E+01;
  COFD[        2378] =  -0.2359825286728792E+00;
  COFD[        2379] =   0.7307807686876559E-02;
  COFD[        2380] =  -0.1906244143057827E+02;
  COFD[        2381] =   0.4049325847301449E+01;
  COFD[        2382] =  -0.2359825286728792E+00;
  COFD[        2383] =   0.7307807686876559E-02;
  COFD[        2384] =  -0.1998941407765696E+02;
  COFD[        2385] =   0.4450214098944700E+01;
  COFD[        2386] =  -0.2925300774906292E+00;
  COFD[        2387] =   0.9925402868644597E-02;
  COFD[        2388] =  -0.2111130556738677E+02;
  COFD[        2389] =   0.5217492138717765E+01;
  COFD[        2390] =  -0.4300187424364486E+00;
  COFD[        2391] =   0.1734092727984471E-01;
  COFD[        2392] =  -0.2112929323949077E+02;
  COFD[        2393] =   0.5220748890969238E+01;
  COFD[        2394] =  -0.4305058272023291E+00;
  COFD[        2395] =   0.1736476920050096E-01;
  COFD[        2396] =  -0.2114446337698859E+02;
  COFD[        2397] =   0.5222999206283971E+01;
  COFD[        2398] =  -0.4308431817261412E+00;
  COFD[        2399] =   0.1738131867299051E-01;
  COFD[        2400] =  -0.2114547150354106E+02;
  COFD[        2401] =   0.5122846710851252E+01;
  COFD[        2402] =  -0.4061105011633706E+00;
  COFD[        2403] =   0.1585327078701023E-01;
  COFD[        2404] =  -0.2123091991916349E+02;
  COFD[        2405] =   0.5158619088859661E+01;
  COFD[        2406] =  -0.4148430725527501E+00;
  COFD[        2407] =   0.1638779761015646E-01;
  COFD[        2408] =  -0.2123706489334075E+02;
  COFD[        2409] =   0.5157473227475789E+01;
  COFD[        2410] =  -0.4146762094312990E+00;
  COFD[        2411] =   0.1637983966492620E-01;
  COFD[        2412] =  -0.2025196923914840E+02;
  COFD[        2413] =   0.5111045666234173E+01;
  COFD[        2414] =  -0.4261812523465015E+00;
  COFD[        2415] =   0.1756072363253846E-01;
  COFD[        2416] =  -0.2039593029382244E+02;
  COFD[        2417] =   0.4595576419324821E+01;
  COFD[        2418] =  -0.3164280753520152E+00;
  COFD[        2419] =   0.1115468584869250E-01;
  COFD[        2420] =  -0.2039593029382244E+02;
  COFD[        2421] =   0.4595576419324821E+01;
  COFD[        2422] =  -0.3164280753520152E+00;
  COFD[        2423] =   0.1115468584869250E-01;
  COFD[        2424] =  -0.2000088490858467E+02;
  COFD[        2425] =   0.5021789747731799E+01;
  COFD[        2426] =  -0.4253480595637520E+00;
  COFD[        2427] =   0.1795667045275812E-01;
  COFD[        2428] =  -0.2043431836412989E+02;
  COFD[        2429] =   0.5116486808224121E+01;
  COFD[        2430] =  -0.4301998553589693E+00;
  COFD[        2431] =   0.1786838693497468E-01;
  COFD[        2432] =  -0.1797854335218850E+02;
  COFD[        2433] =   0.4905020138101169E+01;
  COFD[        2434] =  -0.4293084201751438E+00;
  COFD[        2435] =   0.1891302999394448E-01;
  COFD[        2436] =  -0.1561125740910447E+02;
  COFD[        2437] =   0.3698144089263720E+01;
  COFD[        2438] =  -0.2103080207621357E+00;
  COFD[        2439] =   0.6792204226221277E-02;
  COFD[        2440] =  -0.1878355293596396E+02;
  COFD[        2441] =   0.4745970603113269E+01;
  COFD[        2442] =  -0.3924372163577848E+00;
  COFD[        2443] =   0.1663659729501466E-01;
  COFD[        2444] =  -0.2016013281728147E+02;
  COFD[        2445] =   0.5063566265537697E+01;
  COFD[        2446] =  -0.4281515792120217E+00;
  COFD[        2447] =   0.1797180686456706E-01;
  COFD[        2448] =  -0.1880409131690931E+02;
  COFD[        2449] =   0.4747896766189557E+01;
  COFD[        2450] =  -0.3929513558466455E+00;
  COFD[        2451] =   0.1667071927589917E-01;
  COFD[        2452] =  -0.1619587514529141E+02;
  COFD[        2453] =   0.2843235032975332E+01;
  COFD[        2454] =  -0.5104314541884596E-01;
  COFD[        2455] =  -0.1707937587486378E-02;
  COFD[        2456] =  -0.1994316533928780E+02;
  COFD[        2457] =   0.4989293390342089E+01;
  COFD[        2458] =  -0.4200959859596722E+00;
  COFD[        2459] =   0.1768643162186107E-01;
  COFD[        2460] =  -0.1994945707378163E+02;
  COFD[        2461] =   0.4988909811915903E+01;
  COFD[        2462] =  -0.4200448560138315E+00;
  COFD[        2463] =   0.1768418875431581E-01;
  COFD[        2464] =  -0.1867823422906515E+02;
  COFD[        2465] =   0.4704907071510658E+01;
  COFD[        2466] =  -0.3897247750072582E+00;
  COFD[        2467] =   0.1662675912848450E-01;
  COFD[        2468] =  -0.1871000553285175E+02;
  COFD[        2469] =   0.4738678605504111E+01;
  COFD[        2470] =  -0.3905510969010333E+00;
  COFD[        2471] =   0.1651277300611445E-01;
  COFD[        2472] =  -0.2006425198525411E+02;
  COFD[        2473] =   0.4999965267895282E+01;
  COFD[        2474] =  -0.4103327875262441E+00;
  COFD[        2475] =   0.1680582455088142E-01;
  COFD[        2476] =  -0.2006425198525411E+02;
  COFD[        2477] =   0.4999965267895282E+01;
  COFD[        2478] =  -0.4103327875262441E+00;
  COFD[        2479] =   0.1680582455088142E-01;
  COFD[        2480] =  -0.2012688925964889E+02;
  COFD[        2481] =   0.5018247199767933E+01;
  COFD[        2482] =  -0.4131414445091446E+00;
  COFD[        2483] =   0.1694646037906594E-01;
  COFD[        2484] =  -0.2034321636694622E+02;
  COFD[        2485] =   0.5089562300279502E+01;
  COFD[        2486] =  -0.4211891334428361E+00;
  COFD[        2487] =   0.1724909778947470E-01;
  COFD[        2488] =  -0.2003301781547066E+02;
  COFD[        2489] =   0.5028625403684593E+01;
  COFD[        2490] =  -0.4259384280972174E+00;
  COFD[        2491] =   0.1797007932483027E-01;
  COFD[        2492] =  -0.2111193727907456E+02;
  COFD[        2493] =   0.5115004803331022E+01;
  COFD[        2494] =  -0.4073081193984433E+00;
  COFD[        2495] =   0.1598596435486438E-01;
  COFD[        2496] =  -0.1981646458797120E+02;
  COFD[        2497] =   0.4378886175261416E+01;
  COFD[        2498] =  -0.2813828834600260E+00;
  COFD[        2499] =   0.9370857616443682E-02;
  COFD[        2500] =  -0.1982650591715289E+02;
  COFD[        2501] =   0.4379424471801094E+01;
  COFD[        2502] =  -0.2814553408580349E+00;
  COFD[        2503] =   0.9373962429419786E-02;
  COFD[        2504] =  -0.1906244143057827E+02;
  COFD[        2505] =   0.4049325847301449E+01;
  COFD[        2506] =  -0.2359825286728792E+00;
  COFD[        2507] =   0.7307807686876559E-02;
  COFD[        2508] =  -0.1906244143057827E+02;
  COFD[        2509] =   0.4049325847301449E+01;
  COFD[        2510] =  -0.2359825286728792E+00;
  COFD[        2511] =   0.7307807686876559E-02;
  COFD[        2512] =  -0.1998941407765696E+02;
  COFD[        2513] =   0.4450214098944700E+01;
  COFD[        2514] =  -0.2925300774906292E+00;
  COFD[        2515] =   0.9925402868644597E-02;
  COFD[        2516] =  -0.2111130556738677E+02;
  COFD[        2517] =   0.5217492138717765E+01;
  COFD[        2518] =  -0.4300187424364486E+00;
  COFD[        2519] =   0.1734092727984471E-01;
  COFD[        2520] =  -0.2112929323949077E+02;
  COFD[        2521] =   0.5220748890969238E+01;
  COFD[        2522] =  -0.4305058272023291E+00;
  COFD[        2523] =   0.1736476920050096E-01;
  COFD[        2524] =  -0.2114446337698859E+02;
  COFD[        2525] =   0.5222999206283971E+01;
  COFD[        2526] =  -0.4308431817261412E+00;
  COFD[        2527] =   0.1738131867299051E-01;
  COFD[        2528] =  -0.2114547150354106E+02;
  COFD[        2529] =   0.5122846710851252E+01;
  COFD[        2530] =  -0.4061105011633706E+00;
  COFD[        2531] =   0.1585327078701023E-01;
  COFD[        2532] =  -0.2123091991916349E+02;
  COFD[        2533] =   0.5158619088859661E+01;
  COFD[        2534] =  -0.4148430725527501E+00;
  COFD[        2535] =   0.1638779761015646E-01;
  COFD[        2536] =  -0.2123706489334075E+02;
  COFD[        2537] =   0.5157473227475789E+01;
  COFD[        2538] =  -0.4146762094312990E+00;
  COFD[        2539] =   0.1637983966492620E-01;
  COFD[        2540] =  -0.2025196923914840E+02;
  COFD[        2541] =   0.5111045666234173E+01;
  COFD[        2542] =  -0.4261812523465015E+00;
  COFD[        2543] =   0.1756072363253846E-01;
  COFD[        2544] =  -0.2039593029382244E+02;
  COFD[        2545] =   0.4595576419324821E+01;
  COFD[        2546] =  -0.3164280753520152E+00;
  COFD[        2547] =   0.1115468584869250E-01;
  COFD[        2548] =  -0.2039593029382244E+02;
  COFD[        2549] =   0.4595576419324821E+01;
  COFD[        2550] =  -0.3164280753520152E+00;
  COFD[        2551] =   0.1115468584869250E-01;
  COFD[        2552] =  -0.2000088490858467E+02;
  COFD[        2553] =   0.5021789747731799E+01;
  COFD[        2554] =  -0.4253480595637520E+00;
  COFD[        2555] =   0.1795667045275812E-01;
  COFD[        2556] =  -0.2043431836412989E+02;
  COFD[        2557] =   0.5116486808224121E+01;
  COFD[        2558] =  -0.4301998553589693E+00;
  COFD[        2559] =   0.1786838693497468E-01;
  COFD[        2560] =  -0.1783622363376704E+02;
  COFD[        2561] =   0.4831292682846398E+01;
  COFD[        2562] =  -0.4169745817203137E+00;
  COFD[        2563] =   0.1825729345293071E-01;
  COFD[        2564] =  -0.1545485195366714E+02;
  COFD[        2565] =   0.3597796727187503E+01;
  COFD[        2566] =  -0.1924739949851062E+00;
  COFD[        2567] =   0.5844460086632715E-02;
  COFD[        2568] =  -0.1907199325752651E+02;
  COFD[        2569] =   0.4830462464213029E+01;
  COFD[        2570] =  -0.3998168598650015E+00;
  COFD[        2571] =   0.1680839848129470E-01;
  COFD[        2572] =  -0.2028793707722566E+02;
  COFD[        2573] =   0.5097377789913596E+01;
  COFD[        2574] =  -0.4307531399728672E+00;
  COFD[        2575] =   0.1801467570282063E-01;
  COFD[        2576] =  -0.1910332751367901E+02;
  COFD[        2577] =   0.4837040211978248E+01;
  COFD[        2578] =  -0.4010050297336601E+00;
  COFD[        2579] =   0.1687475885991967E-01;
  COFD[        2580] =  -0.1783328418687383E+02;
  COFD[        2581] =   0.3623064754744803E+01;
  COFD[        2582] =  -0.1668883567851333E+00;
  COFD[        2583] =   0.3813867862735125E-02;
  COFD[        2584] =  -0.2029614681465456E+02;
  COFD[        2585] =   0.5097689110913429E+01;
  COFD[        2586] =  -0.4308107142909775E+00;
  COFD[        2587] =   0.1801795640022364E-01;
  COFD[        2588] =  -0.2030325152362526E+02;
  COFD[        2589] =   0.5097636743004845E+01;
  COFD[        2590] =  -0.4308113500403411E+00;
  COFD[        2591] =   0.1801834516701545E-01;
  COFD[        2592] =  -0.1890046153942015E+02;
  COFD[        2593] =   0.4763373152612735E+01;
  COFD[        2594] =  -0.3936930155883911E+00;
  COFD[        2595] =   0.1664833359988295E-01;
  COFD[        2596] =  -0.1896065537612248E+02;
  COFD[        2597] =   0.4806906431965693E+01;
  COFD[        2598] =  -0.3955804840609516E+00;
  COFD[        2599] =   0.1657242891378949E-01;
  COFD[        2600] =  -0.2031516188866531E+02;
  COFD[        2601] =   0.5061917475702875E+01;
  COFD[        2602] =  -0.4141727881359089E+00;
  COFD[        2603] =   0.1680909956107790E-01;
  COFD[        2604] =  -0.2031516188866531E+02;
  COFD[        2605] =   0.5061917475702875E+01;
  COFD[        2606] =  -0.4141727881359089E+00;
  COFD[        2607] =   0.1680909956107790E-01;
  COFD[        2608] =  -0.2038511261124467E+02;
  COFD[        2609] =   0.5083152894379470E+01;
  COFD[        2610] =  -0.4173863035588838E+00;
  COFD[        2611] =   0.1696799786866059E-01;
  COFD[        2612] =  -0.2040021632891986E+02;
  COFD[        2613] =   0.5094681415071051E+01;
  COFD[        2614] =  -0.4198825806270313E+00;
  COFD[        2615] =   0.1711748447383598E-01;
  COFD[        2616] =  -0.2002461546597822E+02;
  COFD[        2617] =   0.5005680869193603E+01;
  COFD[        2618] =  -0.4206313238956306E+00;
  COFD[        2619] =   0.1764515597235788E-01;
  COFD[        2620] =  -0.2111512317179221E+02;
  COFD[        2621] =   0.5089049253783277E+01;
  COFD[        2622] =  -0.4010600389026057E+00;
  COFD[        2623] =   0.1560812400748400E-01;
  COFD[        2624] =  -0.1910884672535752E+02;
  COFD[        2625] =   0.4028503489934288E+01;
  COFD[        2626] =  -0.2273989668967882E+00;
  COFD[        2627] =   0.6707295013470504E-02;
  COFD[        2628] =  -0.1911879294336443E+02;
  COFD[        2629] =   0.4028925322234755E+01;
  COFD[        2630] =  -0.2274518266873066E+00;
  COFD[        2631] =   0.6709297609983591E-02;
  COFD[        2632] =  -0.1998941407765696E+02;
  COFD[        2633] =   0.4450214098944700E+01;
  COFD[        2634] =  -0.2925300774906292E+00;
  COFD[        2635] =   0.9925402868644597E-02;
  COFD[        2636] =  -0.1998941407765696E+02;
  COFD[        2637] =   0.4450214098944700E+01;
  COFD[        2638] =  -0.2925300774906292E+00;
  COFD[        2639] =   0.9925402868644597E-02;
  COFD[        2640] =  -0.1934717976162822E+02;
  COFD[        2641] =   0.4127136016250256E+01;
  COFD[        2642] =  -0.2423483032940467E+00;
  COFD[        2643] =   0.7437657210565075E-02;
  COFD[        2644] =  -0.2118617333461202E+02;
  COFD[        2645] =   0.5199233449520693E+01;
  COFD[        2646] =  -0.4222135111558194E+00;
  COFD[        2647] =   0.1679135338146060E-01;
  COFD[        2648] =  -0.2120583405822262E+02;
  COFD[        2649] =   0.5203145910746485E+01;
  COFD[        2650] =  -0.4227936349938661E+00;
  COFD[        2651] =   0.1681951784299868E-01;
  COFD[        2652] =  -0.2122260802133029E+02;
  COFD[        2653] =   0.5206031619048144E+01;
  COFD[        2654] =  -0.4232222706241738E+00;
  COFD[        2655] =   0.1684036379103880E-01;
  COFD[        2656] =  -0.2096319952980276E+02;
  COFD[        2657] =   0.4993526798744635E+01;
  COFD[        2658] =  -0.3826329210237235E+00;
  COFD[        2659] =   0.1457641079733719E-01;
  COFD[        2660] =  -0.2122143337952809E+02;
  COFD[        2661] =   0.5104877422506141E+01;
  COFD[        2662] =  -0.4022844500757584E+00;
  COFD[        2663] =   0.1563115659144747E-01;
  COFD[        2664] =  -0.2122916193717567E+02;
  COFD[        2665] =   0.5104386292705568E+01;
  COFD[        2666] =  -0.4022150641632897E+00;
  COFD[        2667] =   0.1562795996905004E-01;
  COFD[        2668] =  -0.2057956090814800E+02;
  COFD[        2669] =   0.5203149310742000E+01;
  COFD[        2670] =  -0.4339608258248750E+00;
  COFD[        2671] =   0.1773356769165475E-01;
  COFD[        2672] =  -0.1994912980553432E+02;
  COFD[        2673] =   0.4354262271960964E+01;
  COFD[        2674] =  -0.2776328852684010E+00;
  COFD[        2675] =   0.9191073490857130E-02;
  COFD[        2676] =  -0.1994912980553432E+02;
  COFD[        2677] =   0.4354262271960964E+01;
  COFD[        2678] =  -0.2776328852684010E+00;
  COFD[        2679] =   0.9191073490857130E-02;
  COFD[        2680] =  -0.1999915273024637E+02;
  COFD[        2681] =   0.5000095313030963E+01;
  COFD[        2682] =  -0.4200303081186300E+00;
  COFD[        2683] =   0.1762406687147381E-01;
  COFD[        2684] =  -0.2065593872227729E+02;
  COFD[        2685] =   0.5165539551318380E+01;
  COFD[        2686] =  -0.4322266808318962E+00;
  COFD[        2687] =   0.1778317395017445E-01;
  COFD[        2688] =  -0.1527556886659549E+02;
  COFD[        2689] =   0.3862563385245359E+01;
  COFD[        2690] =  -0.2980840675934164E+00;
  COFD[        2691] =   0.1340588281258234E-01;
  COFD[        2692] =  -0.1774160477281574E+02;
  COFD[        2693] =   0.4779667070939399E+01;
  COFD[        2694] =  -0.3856864284506084E+00;
  COFD[        2695] =   0.1583474705456999E-01;
  COFD[        2696] =  -0.1728281584130743E+02;
  COFD[        2697] =   0.4249575630985691E+01;
  COFD[        2698] =  -0.3404513157409610E+00;
  COFD[        2699] =   0.1490330472499406E-01;
  COFD[        2700] =  -0.1816374974359966E+02;
  COFD[        2701] =   0.4404539112241986E+01;
  COFD[        2702] =  -0.3569514300098177E+00;
  COFD[        2703] =   0.1547163489736868E-01;
  COFD[        2704] =  -0.1726463880632943E+02;
  COFD[        2705] =   0.4234775238868422E+01;
  COFD[        2706] =  -0.3384933692140595E+00;
  COFD[        2707] =   0.1481714586415791E-01;
  COFD[        2708] =  -0.2027561345681607E+02;
  COFD[        2709] =   0.4923011006260475E+01;
  COFD[        2710] =  -0.3762758154685706E+00;
  COFD[        2711] =   0.1440493997760856E-01;
  COFD[        2712] =  -0.1817178289920526E+02;
  COFD[        2713] =   0.4404990527806040E+01;
  COFD[        2714] =  -0.3570007541989366E+00;
  COFD[        2715] =   0.1547334532679651E-01;
  COFD[        2716] =  -0.1817988211332206E+02;
  COFD[        2717] =   0.4405585855409423E+01;
  COFD[        2718] =  -0.3570675278892711E+00;
  COFD[        2719] =   0.1547575275234058E-01;
  COFD[        2720] =  -0.1700832009372265E+02;
  COFD[        2721] =   0.4130721660167143E+01;
  COFD[        2722] =  -0.3260471778937963E+00;
  COFD[        2723] =   0.1432174203443516E-01;
  COFD[        2724] =  -0.1735339864186368E+02;
  COFD[        2725] =   0.4304534508536137E+01;
  COFD[        2726] =  -0.3477174523791359E+00;
  COFD[        2727] =   0.1522286508630585E-01;
  COFD[        2728] =  -0.1891736997617915E+02;
  COFD[        2729] =   0.4700704222975122E+01;
  COFD[        2730] =  -0.3902105374820116E+00;
  COFD[        2731] =   0.1669545640046241E-01;
  COFD[        2732] =  -0.1891736997617915E+02;
  COFD[        2733] =   0.4700704222975122E+01;
  COFD[        2734] =  -0.3902105374820116E+00;
  COFD[        2735] =   0.1669545640046241E-01;
  COFD[        2736] =  -0.1893404747219508E+02;
  COFD[        2737] =   0.4699655533259762E+01;
  COFD[        2738] =  -0.3902381208946921E+00;
  COFD[        2739] =   0.1670405772171351E-01;
  COFD[        2740] =  -0.1887960780725439E+02;
  COFD[        2741] =   0.4679550122415976E+01;
  COFD[        2742] =  -0.3880189808953476E+00;
  COFD[        2743] =   0.1662436514641822E-01;
  COFD[        2744] =  -0.1800841736946497E+02;
  COFD[        2745] =   0.4351921120886195E+01;
  COFD[        2746] =  -0.3517363610689986E+00;
  COFD[        2747] =   0.1531343604033216E-01;
  COFD[        2748] =  -0.2036669718093571E+02;
  COFD[        2749] =   0.5042701296018135E+01;
  COFD[        2750] =  -0.4228382373074129E+00;
  COFD[        2751] =   0.1763442522229875E-01;
  COFD[        2752] =  -0.2111949107097320E+02;
  COFD[        2753] =   0.5171726543065573E+01;
  COFD[        2754] =  -0.4171067949473953E+00;
  COFD[        2755] =   0.1650848848746980E-01;
  COFD[        2756] =  -0.2112276604025784E+02;
  COFD[        2757] =   0.5169665446870568E+01;
  COFD[        2758] =  -0.4168003974284392E+00;
  COFD[        2759] =   0.1649356899195412E-01;
  COFD[        2760] =  -0.2111130556738677E+02;
  COFD[        2761] =   0.5217492138717765E+01;
  COFD[        2762] =  -0.4300187424364486E+00;
  COFD[        2763] =   0.1734092727984471E-01;
  COFD[        2764] =  -0.2111130556738677E+02;
  COFD[        2765] =   0.5217492138717765E+01;
  COFD[        2766] =  -0.4300187424364486E+00;
  COFD[        2767] =   0.1734092727984471E-01;
  COFD[        2768] =  -0.2118617333461202E+02;
  COFD[        2769] =   0.5199233449520693E+01;
  COFD[        2770] =  -0.4222135111558194E+00;
  COFD[        2771] =   0.1679135338146060E-01;
  COFD[        2772] =  -0.2006339000470825E+02;
  COFD[        2773] =   0.4983892868408139E+01;
  COFD[        2774] =  -0.4201561978851934E+00;
  COFD[        2775] =   0.1772039680073912E-01;
  COFD[        2776] =  -0.2007277841647668E+02;
  COFD[        2777] =   0.4983697515500141E+01;
  COFD[        2778] =  -0.4201241733575779E+00;
  COFD[        2779] =   0.1771871500132603E-01;
  COFD[        2780] =  -0.2008090305086635E+02;
  COFD[        2781] =   0.4983141838534659E+01;
  COFD[        2782] =  -0.4200330807285494E+00;
  COFD[        2783] =   0.1771393120130671E-01;
  COFD[        2784] =  -0.2061603066740442E+02;
  COFD[        2785] =   0.5123762574292296E+01;
  COFD[        2786] =  -0.4301681676037199E+00;
  COFD[        2787] =   0.1782694722668276E-01;
  COFD[        2788] =  -0.2057806613901781E+02;
  COFD[        2789] =   0.5105514013792447E+01;
  COFD[        2790] =  -0.4311352749843362E+00;
  COFD[        2791] =   0.1800316675813045E-01;
  COFD[        2792] =  -0.2058028529774318E+02;
  COFD[        2793] =   0.5102835837876370E+01;
  COFD[        2794] =  -0.4307166347937537E+00;
  COFD[        2795] =   0.1798192512994230E-01;
  COFD[        2796] =  -0.1891361306283375E+02;
  COFD[        2797] =   0.4728417953878362E+01;
  COFD[        2798] =  -0.3935830815508194E+00;
  COFD[        2799] =   0.1683444422800744E-01;
  COFD[        2800] =  -0.2106259252673142E+02;
  COFD[        2801] =   0.5117958864455729E+01;
  COFD[        2802] =  -0.4136130650204043E+00;
  COFD[        2803] =   0.1648627810521946E-01;
  COFD[        2804] =  -0.2106259252673142E+02;
  COFD[        2805] =   0.5117958864455729E+01;
  COFD[        2806] =  -0.4136130650204043E+00;
  COFD[        2807] =   0.1648627810521946E-01;
  COFD[        2808] =  -0.1797950802803312E+02;
  COFD[        2809] =   0.4344100406483967E+01;
  COFD[        2810] =  -0.3507611460911521E+00;
  COFD[        2811] =   0.1527276328163101E-01;
  COFD[        2812] =  -0.1893874700238532E+02;
  COFD[        2813] =   0.4666666006155531E+01;
  COFD[        2814] =  -0.3875266304039654E+00;
  COFD[        2815] =   0.1665285657108702E-01;
  COFD[        2816] =  -0.1528921665984086E+02;
  COFD[        2817] =   0.3867752822044799E+01;
  COFD[        2818] =  -0.2988043102116158E+00;
  COFD[        2819] =   0.1343910541069478E-01;
  COFD[        2820] =  -0.1774507254557799E+02;
  COFD[        2821] =   0.4780159131489816E+01;
  COFD[        2822] =  -0.3856587269477316E+00;
  COFD[        2823] =   0.1582942762638978E-01;
  COFD[        2824] =  -0.1731509452354737E+02;
  COFD[        2825] =   0.4259553692469185E+01;
  COFD[        2826] =  -0.3417709853226915E+00;
  COFD[        2827] =   0.1496136243888912E-01;
  COFD[        2828] =  -0.1817374730686181E+02;
  COFD[        2829] =   0.4404209750065649E+01;
  COFD[        2830] =  -0.3569181269134391E+00;
  COFD[        2831] =   0.1547062305483025E-01;
  COFD[        2832] =  -0.1729594517981932E+02;
  COFD[        2833] =   0.4244246410490659E+01;
  COFD[        2834] =  -0.3397463840897515E+00;
  COFD[        2835] =   0.1487228746520592E-01;
  COFD[        2836] =  -0.2024506268677845E+02;
  COFD[        2837] =   0.4905829896937226E+01;
  COFD[        2838] =  -0.3737705712203094E+00;
  COFD[        2839] =   0.1428535217779750E-01;
  COFD[        2840] =  -0.1818135452916891E+02;
  COFD[        2841] =   0.4404445109750770E+01;
  COFD[        2842] =  -0.3569414786935936E+00;
  COFD[        2843] =   0.1547130685514495E-01;
  COFD[        2844] =  -0.1818905402728113E+02;
  COFD[        2845] =   0.4404836108862263E+01;
  COFD[        2846] =  -0.3569836886410043E+00;
  COFD[        2847] =   0.1547274326407077E-01;
  COFD[        2848] =  -0.1704125657308844E+02;
  COFD[        2849] =   0.4141672233114677E+01;
  COFD[        2850] =  -0.3275236148236659E+00;
  COFD[        2851] =   0.1438795853925400E-01;
  COFD[        2852] =  -0.1738771924301953E+02;
  COFD[        2853] =   0.4315709821275419E+01;
  COFD[        2854] =  -0.3491942943333771E+00;
  COFD[        2855] =   0.1528778889586191E-01;
  COFD[        2856] =  -0.1892750366101823E+02;
  COFD[        2857] =   0.4701319071264327E+01;
  COFD[        2858] =  -0.3901920115746619E+00;
  COFD[        2859] =   0.1669020608030539E-01;
  COFD[        2860] =  -0.1892750366101823E+02;
  COFD[        2861] =   0.4701319071264327E+01;
  COFD[        2862] =  -0.3901920115746619E+00;
  COFD[        2863] =   0.1669020608030539E-01;
  COFD[        2864] =  -0.1894432968047615E+02;
  COFD[        2865] =   0.4700246655831890E+01;
  COFD[        2866] =  -0.3902232238988751E+00;
  COFD[        2867] =   0.1669926672613195E-01;
  COFD[        2868] =  -0.1889022870587310E+02;
  COFD[        2869] =   0.4680203446785792E+01;
  COFD[        2870] =  -0.3880194721331493E+00;
  COFD[        2871] =   0.1662057070374977E-01;
  COFD[        2872] =  -0.1801976002924096E+02;
  COFD[        2873] =   0.4352352585107287E+01;
  COFD[        2874] =  -0.3517965558706893E+00;
  COFD[        2875] =   0.1531621848405548E-01;
  COFD[        2876] =  -0.2039137163285755E+02;
  COFD[        2877] =   0.5048551827610241E+01;
  COFD[        2878] =  -0.4237554975633671E+00;
  COFD[        2879] =   0.1768105924766743E-01;
  COFD[        2880] =  -0.2113317649345642E+02;
  COFD[        2881] =   0.5173175010557953E+01;
  COFD[        2882] =  -0.4173234769896741E+00;
  COFD[        2883] =   0.1651910472359747E-01;
  COFD[        2884] =  -0.2113855609134408E+02;
  COFD[        2885] =   0.5171981200778053E+01;
  COFD[        2886] =  -0.4171447495559704E+00;
  COFD[        2887] =   0.1651034132923998E-01;
  COFD[        2888] =  -0.2112929323949077E+02;
  COFD[        2889] =   0.5220748890969238E+01;
  COFD[        2890] =  -0.4305058272023291E+00;
  COFD[        2891] =   0.1736476920050096E-01;
  COFD[        2892] =  -0.2112929323949077E+02;
  COFD[        2893] =   0.5220748890969238E+01;
  COFD[        2894] =  -0.4305058272023291E+00;
  COFD[        2895] =   0.1736476920050096E-01;
  COFD[        2896] =  -0.2120583405822262E+02;
  COFD[        2897] =   0.5203145910746485E+01;
  COFD[        2898] =  -0.4227936349938661E+00;
  COFD[        2899] =   0.1681951784299868E-01;
  COFD[        2900] =  -0.2007277841647668E+02;
  COFD[        2901] =   0.4983697515500141E+01;
  COFD[        2902] =  -0.4201241733575779E+00;
  COFD[        2903] =   0.1771871500132603E-01;
  COFD[        2904] =  -0.2008313017084402E+02;
  COFD[        2905] =   0.4983892868408141E+01;
  COFD[        2906] =  -0.4201561978851936E+00;
  COFD[        2907] =   0.1772039680073913E-01;
  COFD[        2908] =  -0.2009217954894486E+02;
  COFD[        2909] =   0.4983712068565962E+01;
  COFD[        2910] =  -0.4201265590632418E+00;
  COFD[        2911] =   0.1771884028887708E-01;
  COFD[        2912] =  -0.2062921491889363E+02;
  COFD[        2913] =   0.5125110770555073E+01;
  COFD[        2914] =  -0.4303770100880492E+00;
  COFD[        2915] =   0.1783747461488234E-01;
  COFD[        2916] =  -0.2059391054107079E+02;
  COFD[        2917] =   0.5108068575679129E+01;
  COFD[        2918] =  -0.4315347958033334E+00;
  COFD[        2919] =   0.1802344633255497E-01;
  COFD[        2920] =  -0.2059735929037355E+02;
  COFD[        2921] =   0.5105897247055966E+01;
  COFD[        2922] =  -0.4311951951062819E+00;
  COFD[        2923] =   0.1800620766882382E-01;
  COFD[        2924] =  -0.1892576390190177E+02;
  COFD[        2925] =   0.4728500364939531E+01;
  COFD[        2926] =  -0.3936077379198002E+00;
  COFD[        2927] =   0.1683614402582818E-01;
  COFD[        2928] =  -0.2109960034604860E+02;
  COFD[        2929] =   0.5129133288825352E+01;
  COFD[        2930] =  -0.4152735461892595E+00;
  COFD[        2931] =   0.1656709485798127E-01;
  COFD[        2932] =  -0.2109960034604860E+02;
  COFD[        2933] =   0.5129133288825352E+01;
  COFD[        2934] =  -0.4152735461892595E+00;
  COFD[        2935] =   0.1656709485798127E-01;
  COFD[        2936] =  -0.1799101249332176E+02;
  COFD[        2937] =   0.4344595334851394E+01;
  COFD[        2938] =  -0.3508294711072030E+00;
  COFD[        2939] =   0.1527589170814564E-01;
  COFD[        2940] =  -0.1894826439081699E+02;
  COFD[        2941] =   0.4665904241842780E+01;
  COFD[        2942] =  -0.3874653193184607E+00;
  COFD[        2943] =   0.1665191589170439E-01;
  COFD[        2944] =  -0.1530201153719761E+02;
  COFD[        2945] =   0.3872623565106621E+01;
  COFD[        2946] =  -0.2994802960421091E+00;
  COFD[        2947] =   0.1347028558095244E-01;
  COFD[        2948] =  -0.1774834356244655E+02;
  COFD[        2949] =   0.4780626617651506E+01;
  COFD[        2950] =  -0.3856327873318030E+00;
  COFD[        2951] =   0.1582440706500046E-01;
  COFD[        2952] =  -0.1734665788392033E+02;
  COFD[        2953] =   0.4269427891180465E+01;
  COFD[        2954] =  -0.3430766910038187E+00;
  COFD[        2955] =   0.1501879627721415E-01;
  COFD[        2956] =  -0.1818389332535593E+02;
  COFD[        2957] =   0.4404158898744668E+01;
  COFD[        2958] =  -0.3569180632138942E+00;
  COFD[        2959] =   0.1547091340944365E-01;
  COFD[        2960] =  -0.1732667376192460E+02;
  COFD[        2961] =   0.4253674558730024E+01;
  COFD[        2962] =  -0.3409934577271315E+00;
  COFD[        2963] =   0.1492715703864712E-01;
  COFD[        2964] =  -0.2021382909884265E+02;
  COFD[        2965] =   0.4888560960839369E+01;
  COFD[        2966] =  -0.3712523102934615E+00;
  COFD[        2967] =   0.1416513273755347E-01;
  COFD[        2968] =  -0.1819107613580373E+02;
  COFD[        2969] =   0.4404180695515474E+01;
  COFD[        2970] =  -0.3569157439691602E+00;
  COFD[        2971] =   0.1547058259292524E-01;
  COFD[        2972] =  -0.1819837498757435E+02;
  COFD[        2973] =   0.4404368843395893E+01;
  COFD[        2974] =  -0.3569335691748086E+00;
  COFD[        2975] =   0.1547105514857662E-01;
  COFD[        2976] =  -0.1707305377045087E+02;
  COFD[        2977] =   0.4152317255109030E+01;
  COFD[        2978] =  -0.3289587183559789E+00;
  COFD[        2979] =   0.1445231595352617E-01;
  COFD[        2980] =  -0.1742090070168673E+02;
  COFD[        2981] =   0.4326592022152555E+01;
  COFD[        2982] =  -0.3506322276126634E+00;
  COFD[        2983] =   0.1535099522033430E-01;
  COFD[        2984] =  -0.1893716028173049E+02;
  COFD[        2985] =   0.4701920423433434E+01;
  COFD[        2986] =  -0.3901723473212924E+00;
  COFD[        2987] =   0.1668493495004967E-01;
  COFD[        2988] =  -0.1893716028173049E+02;
  COFD[        2989] =   0.4701920423433434E+01;
  COFD[        2990] =  -0.3901723473212924E+00;
  COFD[        2991] =   0.1668493495004967E-01;
  COFD[        2992] =  -0.1895414515078359E+02;
  COFD[        2993] =   0.4700828894373993E+01;
  COFD[        2994] =  -0.3902069149863984E+00;
  COFD[        2995] =   0.1669440360391662E-01;
  COFD[        2996] =  -0.1890039213728353E+02;
  COFD[        2997] =   0.4680852507758593E+01;
  COFD[        2998] =  -0.3880184127527163E+00;
  COFD[        2999] =   0.1661666470021389E-01;
  COFD[        3000] =  -0.1803148637733784E+02;
  COFD[        3001] =   0.4353159047512093E+01;
  COFD[        3002] =  -0.3519044208218940E+00;
  COFD[        3003] =   0.1532101211937375E-01;
  COFD[        3004] =  -0.2041424604605752E+02;
  COFD[        3005] =   0.5053848971414935E+01;
  COFD[        3006] =  -0.4245860964146182E+00;
  COFD[        3007] =   0.1772329162485299E-01;
  COFD[        3008] =  -0.2114401307743194E+02;
  COFD[        3009] =   0.5173607724136159E+01;
  COFD[        3010] =  -0.4173904988942925E+00;
  COFD[        3011] =   0.1652249796434591E-01;
  COFD[        3012] =  -0.2115144768087362E+02;
  COFD[        3013] =   0.5173261676775666E+01;
  COFD[        3014] =  -0.4173365706093256E+00;
  COFD[        3015] =   0.1651975240232191E-01;
  COFD[        3016] =  -0.2114446337698859E+02;
  COFD[        3017] =   0.5222999206283971E+01;
  COFD[        3018] =  -0.4308431817261412E+00;
  COFD[        3019] =   0.1738131867299051E-01;
  COFD[        3020] =  -0.2114446337698859E+02;
  COFD[        3021] =   0.5222999206283971E+01;
  COFD[        3022] =  -0.4308431817261412E+00;
  COFD[        3023] =   0.1738131867299051E-01;
  COFD[        3024] =  -0.2122260802133029E+02;
  COFD[        3025] =   0.5206031619048144E+01;
  COFD[        3026] =  -0.4232222706241738E+00;
  COFD[        3027] =   0.1684036379103880E-01;
  COFD[        3028] =  -0.2008090305086635E+02;
  COFD[        3029] =   0.4983141838534659E+01;
  COFD[        3030] =  -0.4200330807285494E+00;
  COFD[        3031] =   0.1771393120130671E-01;
  COFD[        3032] =  -0.2009217954894486E+02;
  COFD[        3033] =   0.4983712068565962E+01;
  COFD[        3034] =  -0.4201265590632418E+00;
  COFD[        3035] =   0.1771884028887708E-01;
  COFD[        3036] =  -0.2010212049883026E+02;
  COFD[        3037] =   0.4983892868408190E+01;
  COFD[        3038] =  -0.4201561978852008E+00;
  COFD[        3039] =   0.1772039680073947E-01;
  COFD[        3040] =  -0.2064026808809792E+02;
  COFD[        3041] =   0.5125729516119874E+01;
  COFD[        3042] =  -0.4304737187234308E+00;
  COFD[        3043] =   0.1784238527988932E-01;
  COFD[        3044] =  -0.2060778849004208E+02;
  COFD[        3045] =   0.5109962431730637E+01;
  COFD[        3046] =  -0.4318312376520599E+00;
  COFD[        3047] =   0.1803850360308551E-01;
  COFD[        3048] =  -0.2061243965594226E+02;
  COFD[        3049] =   0.5108287844008140E+01;
  COFD[        3050] =  -0.4315691026948430E+00;
  COFD[        3051] =   0.1802518831107918E-01;
  COFD[        3052] =  -0.1893741821656993E+02;
  COFD[        3053] =   0.4728594869696905E+01;
  COFD[        3054] =  -0.3936271603294507E+00;
  COFD[        3055] =   0.1683730982264022E-01;
  COFD[        3056] =  -0.2113363149383363E+02;
  COFD[        3057] =   0.5139253110457338E+01;
  COFD[        3058] =  -0.4167775115779367E+00;
  COFD[        3059] =   0.1664030318824075E-01;
  COFD[        3060] =  -0.2113363149383363E+02;
  COFD[        3061] =   0.5139253110457338E+01;
  COFD[        3062] =  -0.4167775115779367E+00;
  COFD[        3063] =   0.1664030318824075E-01;
  COFD[        3064] =  -0.1800289314001846E+02;
  COFD[        3065] =   0.4345462384297714E+01;
  COFD[        3066] =  -0.3509451050542303E+00;
  COFD[        3067] =   0.1528101631398434E-01;
  COFD[        3068] =  -0.1895758159603200E+02;
  COFD[        3069] =   0.4665279128753897E+01;
  COFD[        3070] =  -0.3874162488122672E+00;
  COFD[        3071] =   0.1665125384367453E-01;
  COFD[        3072] =  -0.1646041046594534E+02;
  COFD[        3073] =   0.4329097783536894E+01;
  COFD[        3074] =  -0.3588538378541796E+00;
  COFD[        3075] =   0.1604603061706265E-01;
  COFD[        3076] =  -0.1694497223361135E+02;
  COFD[        3077] =   0.4353290721870846E+01;
  COFD[        3078] =  -0.3147441391485086E+00;
  COFD[        3079] =   0.1210201943729390E-01;
  COFD[        3080] =  -0.1800091724620567E+02;
  COFD[        3081] =   0.4493072014688184E+01;
  COFD[        3082] =  -0.3678451189139861E+00;
  COFD[        3083] =   0.1591448159607913E-01;
  COFD[        3084] =  -0.1901929468754062E+02;
  COFD[        3085] =   0.4693839754776671E+01;
  COFD[        3086] =  -0.3901000075815114E+00;
  COFD[        3087] =   0.1672631492507185E-01;
  COFD[        3088] =  -0.1799417924090561E+02;
  COFD[        3089] =   0.4483237269191199E+01;
  COFD[        3090] =  -0.3666457816700467E+00;
  COFD[        3091] =   0.1586624138953086E-01;
  COFD[        3092] =  -0.1988688758541057E+02;
  COFD[        3093] =   0.4655348588622681E+01;
  COFD[        3094] =  -0.3287200044949488E+00;
  COFD[        3095] =   0.1186464366263145E-01;
  COFD[        3096] =  -0.1902570873205879E+02;
  COFD[        3097] =   0.4693487292163261E+01;
  COFD[        3098] =  -0.3900471625097403E+00;
  COFD[        3099] =   0.1672372332121293E-01;
  COFD[        3100] =  -0.1903192515013230E+02;
  COFD[        3101] =   0.4693162157411670E+01;
  COFD[        3102] =  -0.3899939058159093E+00;
  COFD[        3103] =   0.1672093460134929E-01;
  COFD[        3104] =  -0.1790606312284444E+02;
  COFD[        3105] =   0.4454152121123291E+01;
  COFD[        3106] =  -0.3650609129312846E+00;
  COFD[        3107] =   0.1588913670146698E-01;
  COFD[        3108] =  -0.1802582448642654E+02;
  COFD[        3109] =   0.4528207779823734E+01;
  COFD[        3110] =  -0.3721205839071864E+00;
  COFD[        3111] =   0.1608601965175457E-01;
  COFD[        3112] =  -0.1956675647234025E+02;
  COFD[        3113] =   0.4893594117189979E+01;
  COFD[        3114] =  -0.4083224878273328E+00;
  COFD[        3115] =   0.1719480931404233E-01;
  COFD[        3116] =  -0.1956675647234025E+02;
  COFD[        3117] =   0.4893594117189979E+01;
  COFD[        3118] =  -0.4083224878273328E+00;
  COFD[        3119] =   0.1719480931404233E-01;
  COFD[        3120] =  -0.1960184896363511E+02;
  COFD[        3121] =   0.4900346283733346E+01;
  COFD[        3122] =  -0.4094945328875742E+00;
  COFD[        3123] =   0.1725872606725634E-01;
  COFD[        3124] =  -0.1957015754946189E+02;
  COFD[        3125] =   0.4890650507507935E+01;
  COFD[        3126] =  -0.4088565156459115E+00;
  COFD[        3127] =   0.1725720932243744E-01;
  COFD[        3128] =  -0.1879103964706305E+02;
  COFD[        3129] =   0.4613167070471964E+01;
  COFD[        3130] =  -0.3811821259847005E+00;
  COFD[        3131] =   0.1640353922053438E-01;
  COFD[        3132] =  -0.2089328018920450E+02;
  COFD[        3133] =   0.5166930585415199E+01;
  COFD[        3134] =  -0.4306289589876920E+00;
  COFD[        3135] =   0.1764267219152353E-01;
  COFD[        3136] =  -0.2090933405715619E+02;
  COFD[        3137] =   0.4972999527197874E+01;
  COFD[        3138] =  -0.3786889263628220E+00;
  COFD[        3139] =   0.1435806171655258E-01;
  COFD[        3140] =  -0.2091775085800750E+02;
  COFD[        3141] =   0.4973029974129495E+01;
  COFD[        3142] =  -0.3786913370630364E+00;
  COFD[        3143] =   0.1435807083846824E-01;
  COFD[        3144] =  -0.2114547150354106E+02;
  COFD[        3145] =   0.5122846710851252E+01;
  COFD[        3146] =  -0.4061105011633706E+00;
  COFD[        3147] =   0.1585327078701023E-01;
  COFD[        3148] =  -0.2114547150354106E+02;
  COFD[        3149] =   0.5122846710851252E+01;
  COFD[        3150] =  -0.4061105011633706E+00;
  COFD[        3151] =   0.1585327078701023E-01;
  COFD[        3152] =  -0.2096319952980276E+02;
  COFD[        3153] =   0.4993526798744635E+01;
  COFD[        3154] =  -0.3826329210237235E+00;
  COFD[        3155] =   0.1457641079733719E-01;
  COFD[        3156] =  -0.2061603066740442E+02;
  COFD[        3157] =   0.5123762574292296E+01;
  COFD[        3158] =  -0.4301681676037199E+00;
  COFD[        3159] =   0.1782694722668276E-01;
  COFD[        3160] =  -0.2062921491889363E+02;
  COFD[        3161] =   0.5125110770555073E+01;
  COFD[        3162] =  -0.4303770100880492E+00;
  COFD[        3163] =   0.1783747461488234E-01;
  COFD[        3164] =  -0.2064026808809792E+02;
  COFD[        3165] =   0.5125729516119874E+01;
  COFD[        3166] =  -0.4304737187234308E+00;
  COFD[        3167] =   0.1784238527988932E-01;
  COFD[        3168] =  -0.2117733141272616E+02;
  COFD[        3169] =   0.5262507718327111E+01;
  COFD[        3170] =  -0.4398579014089885E+00;
  COFD[        3171] =   0.1792489830320274E-01;
  COFD[        3172] =  -0.2108189705348579E+02;
  COFD[        3173] =   0.5219096733932689E+01;
  COFD[        3174] =  -0.4371166395377749E+00;
  COFD[        3175] =   0.1791475976287751E-01;
  COFD[        3176] =  -0.2108550712837031E+02;
  COFD[        3177] =   0.5216945829110916E+01;
  COFD[        3178] =  -0.4367926562187473E+00;
  COFD[        3179] =   0.1789880186533132E-01;
  COFD[        3180] =  -0.1961481660694155E+02;
  COFD[        3181] =   0.4942832081837352E+01;
  COFD[        3182] =  -0.4152110398291045E+00;
  COFD[        3183] =   0.1752038373578458E-01;
  COFD[        3184] =  -0.2115552988307172E+02;
  COFD[        3185] =   0.5049408694314179E+01;
  COFD[        3186] =  -0.3939342769912307E+00;
  COFD[        3187] =   0.1522309753276939E-01;
  COFD[        3188] =  -0.2115552988307172E+02;
  COFD[        3189] =   0.5049408694314179E+01;
  COFD[        3190] =  -0.3939342769912307E+00;
  COFD[        3191] =   0.1522309753276939E-01;
  COFD[        3192] =  -0.1876345278233955E+02;
  COFD[        3193] =   0.4606235123260915E+01;
  COFD[        3194] =  -0.3803641770348614E+00;
  COFD[        3195] =   0.1637150008428244E-01;
  COFD[        3196] =  -0.1962450162993201E+02;
  COFD[        3197] =   0.4875751864229406E+01;
  COFD[        3198] =  -0.4087050411778894E+00;
  COFD[        3199] =   0.1732743672805724E-01;
  COFD[        3200] =  -0.1615267929807825E+02;
  COFD[        3201] =   0.4174121490551778E+01;
  COFD[        3202] =  -0.3383475624805202E+00;
  COFD[        3203] =   0.1513993577326728E-01;
  COFD[        3204] =  -0.1738386031480564E+02;
  COFD[        3205] =   0.4515894896789344E+01;
  COFD[        3206] =  -0.3391856065220307E+00;
  COFD[        3207] =   0.1329276720561771E-01;
  COFD[        3208] =  -0.1802433609839966E+02;
  COFD[        3209] =   0.4483646303442399E+01;
  COFD[        3210] =  -0.3688151334706006E+00;
  COFD[        3211] =   0.1604826293273014E-01;
  COFD[        3212] =  -0.1881741117557303E+02;
  COFD[        3213] =   0.4595519619722809E+01;
  COFD[        3214] =  -0.3791053358880415E+00;
  COFD[        3215] =   0.1632232422675412E-01;
  COFD[        3216] =  -0.1800542957417654E+02;
  COFD[        3217] =   0.4468504271448043E+01;
  COFD[        3218] =  -0.3668692513471102E+00;
  COFD[        3219] =   0.1596519376748332E-01;
  COFD[        3220] =  -0.2010033725040245E+02;
  COFD[        3221] =   0.4740608584978043E+01;
  COFD[        3222] =  -0.3442903017689571E+00;
  COFD[        3223] =   0.1271241178230842E-01;
  COFD[        3224] =  -0.1882326100183751E+02;
  COFD[        3225] =   0.4594957876626870E+01;
  COFD[        3226] =  -0.3790360690538448E+00;
  COFD[        3227] =   0.1631950746180167E-01;
  COFD[        3228] =  -0.1882911310545702E+02;
  COFD[        3229] =   0.4594509279764214E+01;
  COFD[        3230] =  -0.3789782571863847E+00;
  COFD[        3231] =   0.1631703818360878E-01;
  COFD[        3232] =  -0.1792078479501453E+02;
  COFD[        3233] =   0.4443307443203482E+01;
  COFD[        3234] =  -0.3657004955367218E+00;
  COFD[        3235] =   0.1600230799211683E-01;
  COFD[        3236] =  -0.1809136688643801E+02;
  COFD[        3237] =   0.4537201706428067E+01;
  COFD[        3238] =  -0.3756919101622764E+00;
  COFD[        3239] =   0.1634159264026368E-01;
  COFD[        3240] =  -0.1941109254354329E+02;
  COFD[        3241] =   0.4819578196751006E+01;
  COFD[        3242] =  -0.4008189828944452E+00;
  COFD[        3243] =   0.1695355424415370E-01;
  COFD[        3244] =  -0.1941109254354329E+02;
  COFD[        3245] =   0.4819578196751006E+01;
  COFD[        3246] =  -0.4008189828944452E+00;
  COFD[        3247] =   0.1695355424415370E-01;
  COFD[        3248] =  -0.1943925299002717E+02;
  COFD[        3249] =   0.4823356728877616E+01;
  COFD[        3250] =  -0.4015874038077185E+00;
  COFD[        3251] =   0.1699928308073373E-01;
  COFD[        3252] =  -0.1938496483487797E+02;
  COFD[        3253] =   0.4803651837482498E+01;
  COFD[        3254] =  -0.3995406103851766E+00;
  COFD[        3255] =   0.1693221082596063E-01;
  COFD[        3256] =  -0.1864571977877140E+02;
  COFD[        3257] =   0.4539069966095800E+01;
  COFD[        3258] =  -0.3735011850255855E+00;
  COFD[        3259] =   0.1615133908644331E-01;
  COFD[        3260] =  -0.2075862463488640E+02;
  COFD[        3261] =   0.5108417296497037E+01;
  COFD[        3262] =  -0.4262851069696628E+00;
  COFD[        3263] =   0.1758280414396813E-01;
  COFD[        3264] =  -0.2110606940488274E+02;
  COFD[        3265] =   0.5056386072237322E+01;
  COFD[        3266] =  -0.3941314241059736E+00;
  COFD[        3267] =   0.1520386988001248E-01;
  COFD[        3268] =  -0.2111907495038380E+02;
  COFD[        3269] =   0.5058381671032987E+01;
  COFD[        3270] =  -0.3944196800959486E+00;
  COFD[        3271] =   0.1521748580862574E-01;
  COFD[        3272] =  -0.2123091991916349E+02;
  COFD[        3273] =   0.5158619088859661E+01;
  COFD[        3274] =  -0.4148430725527501E+00;
  COFD[        3275] =   0.1638779761015646E-01;
  COFD[        3276] =  -0.2123091991916349E+02;
  COFD[        3277] =   0.5158619088859661E+01;
  COFD[        3278] =  -0.4148430725527501E+00;
  COFD[        3279] =   0.1638779761015646E-01;
  COFD[        3280] =  -0.2122143337952809E+02;
  COFD[        3281] =   0.5104877422506141E+01;
  COFD[        3282] =  -0.4022844500757584E+00;
  COFD[        3283] =   0.1563115659144747E-01;
  COFD[        3284] =  -0.2057806613901781E+02;
  COFD[        3285] =   0.5105514013792447E+01;
  COFD[        3286] =  -0.4311352749843362E+00;
  COFD[        3287] =   0.1800316675813045E-01;
  COFD[        3288] =  -0.2059391054107079E+02;
  COFD[        3289] =   0.5108068575679129E+01;
  COFD[        3290] =  -0.4315347958033334E+00;
  COFD[        3291] =   0.1802344633255497E-01;
  COFD[        3292] =  -0.2060778849004208E+02;
  COFD[        3293] =   0.5109962431730637E+01;
  COFD[        3294] =  -0.4318312376520599E+00;
  COFD[        3295] =   0.1803850360308551E-01;
  COFD[        3296] =  -0.2108189705348579E+02;
  COFD[        3297] =   0.5219096733932689E+01;
  COFD[        3298] =  -0.4371166395377749E+00;
  COFD[        3299] =   0.1791475976287751E-01;
  COFD[        3300] =  -0.2091185026658962E+02;
  COFD[        3301] =   0.5148560230790555E+01;
  COFD[        3302] =  -0.4310562392852861E+00;
  COFD[        3303] =   0.1777358978417338E-01;
  COFD[        3304] =  -0.2091965964592144E+02;
  COFD[        3305] =   0.5148262389934279E+01;
  COFD[        3306] =  -0.4310108894089904E+00;
  COFD[        3307] =   0.1777133583281057E-01;
  COFD[        3308] =  -0.1945124782431341E+02;
  COFD[        3309] =   0.4857428845983466E+01;
  COFD[        3310] =  -0.4062322300213044E+00;
  COFD[        3311] =   0.1721513119766667E-01;
  COFD[        3312] =  -0.2129367686674487E+02;
  COFD[        3313] =   0.5109799100487325E+01;
  COFD[        3314] =  -0.4062923797917856E+00;
  COFD[        3315] =   0.1593137641849987E-01;
  COFD[        3316] =  -0.2129367686674487E+02;
  COFD[        3317] =   0.5109799100487325E+01;
  COFD[        3318] =  -0.4062923797917856E+00;
  COFD[        3319] =   0.1593137641849987E-01;
  COFD[        3320] =  -0.1861866210245701E+02;
  COFD[        3321] =   0.4532060302179115E+01;
  COFD[        3322] =  -0.3726516669477872E+00;
  COFD[        3323] =   0.1611699967368039E-01;
  COFD[        3324] =  -0.1943393371703941E+02;
  COFD[        3325] =   0.4784664368892784E+01;
  COFD[        3326] =  -0.3991077435475853E+00;
  COFD[        3327] =   0.1700210403998681E-01;
  COFD[        3328] =  -0.1616515483252968E+02;
  COFD[        3329] =   0.4178910660883203E+01;
  COFD[        3330] =  -0.3390062877759687E+00;
  COFD[        3331] =   0.1517006334795817E-01;
  COFD[        3332] =  -0.1738223879329547E+02;
  COFD[        3333] =   0.4514299549461490E+01;
  COFD[        3334] =  -0.3388651702152571E+00;
  COFD[        3335] =   0.1327383213284630E-01;
  COFD[        3336] =  -0.1805221831789452E+02;
  COFD[        3337] =   0.4492279525932131E+01;
  COFD[        3338] =  -0.3699242263001061E+00;
  COFD[        3339] =   0.1609559465701492E-01;
  COFD[        3340] =  -0.1882824574536174E+02;
  COFD[        3341] =   0.4596254296231404E+01;
  COFD[        3342] =  -0.3791928797888692E+00;
  COFD[        3343] =   0.1632573982847240E-01;
  COFD[        3344] =  -0.1803293840758895E+02;
  COFD[        3345] =   0.4476898762718537E+01;
  COFD[        3346] =  -0.3679481179513786E+00;
  COFD[        3347] =   0.1601125468656841E-01;
  COFD[        3348] =  -0.2007381993778366E+02;
  COFD[        3349] =   0.4726146050811121E+01;
  COFD[        3350] =  -0.3422109775805109E+00;
  COFD[        3351] =   0.1261488904944177E-01;
  COFD[        3352] =  -0.1883388798798641E+02;
  COFD[        3353] =   0.4595580957493091E+01;
  COFD[        3354] =  -0.3791127554800276E+00;
  COFD[        3355] =   0.1632261913288707E-01;
  COFD[        3356] =  -0.1883954098397452E+02;
  COFD[        3357] =   0.4595024886620513E+01;
  COFD[        3358] =  -0.3790444797196689E+00;
  COFD[        3359] =   0.1631985649333879E-01;
  COFD[        3360] =  -0.1794869549063484E+02;
  COFD[        3361] =   0.4452595998578516E+01;
  COFD[        3362] =  -0.3669291177144341E+00;
  COFD[        3363] =   0.1605636779086480E-01;
  COFD[        3364] =  -0.1811958823661455E+02;
  COFD[        3365] =   0.4546265606449557E+01;
  COFD[        3366] =  -0.3768550606974422E+00;
  COFD[        3367] =   0.1639117743141738E-01;
  COFD[        3368] =  -0.1941401903486105E+02;
  COFD[        3369] =   0.4817662369355102E+01;
  COFD[        3370] =  -0.4004301095900044E+00;
  COFD[        3371] =   0.1693043453375954E-01;
  COFD[        3372] =  -0.1941401903486105E+02;
  COFD[        3373] =   0.4817662369355102E+01;
  COFD[        3374] =  -0.4004301095900044E+00;
  COFD[        3375] =   0.1693043453375954E-01;
  COFD[        3376] =  -0.1944256474839717E+02;
  COFD[        3377] =   0.4821522129341849E+01;
  COFD[        3378] =  -0.4012140535427826E+00;
  COFD[        3379] =   0.1697705724216803E-01;
  COFD[        3380] =  -0.1938900888397246E+02;
  COFD[        3381] =   0.4802051950424429E+01;
  COFD[        3382] =  -0.3992044152028926E+00;
  COFD[        3383] =   0.1691189193038190E-01;
  COFD[        3384] =  -0.1865835471859333E+02;
  COFD[        3385] =   0.4540734987509840E+01;
  COFD[        3386] =  -0.3737078963922989E+00;
  COFD[        3387] =   0.1615982140755382E-01;
  COFD[        3388] =  -0.2077840705218999E+02;
  COFD[        3389] =   0.5112879719806333E+01;
  COFD[        3390] =  -0.4269664397363305E+00;
  COFD[        3391] =   0.1761673548420319E-01;
  COFD[        3392] =  -0.2110836259086301E+02;
  COFD[        3393] =   0.5053674020415143E+01;
  COFD[        3394] =  -0.3937388243323189E+00;
  COFD[        3395] =   0.1518528106281032E-01;
  COFD[        3396] =  -0.2112313758034541E+02;
  COFD[        3397] =   0.5056389549388699E+01;
  COFD[        3398] =  -0.3941319269977342E+00;
  COFD[        3399] =   0.1520389366699413E-01;
  COFD[        3400] =  -0.2123706489334075E+02;
  COFD[        3401] =   0.5157473227475789E+01;
  COFD[        3402] =  -0.4146762094312990E+00;
  COFD[        3403] =   0.1637983966492620E-01;
  COFD[        3404] =  -0.2123706489334075E+02;
  COFD[        3405] =   0.5157473227475789E+01;
  COFD[        3406] =  -0.4146762094312990E+00;
  COFD[        3407] =   0.1637983966492620E-01;
  COFD[        3408] =  -0.2122916193717567E+02;
  COFD[        3409] =   0.5104386292705568E+01;
  COFD[        3410] =  -0.4022150641632897E+00;
  COFD[        3411] =   0.1562795996905004E-01;
  COFD[        3412] =  -0.2058028529774318E+02;
  COFD[        3413] =   0.5102835837876370E+01;
  COFD[        3414] =  -0.4307166347937537E+00;
  COFD[        3415] =   0.1798192512994230E-01;
  COFD[        3416] =  -0.2059735929037355E+02;
  COFD[        3417] =   0.5105897247055966E+01;
  COFD[        3418] =  -0.4311951951062819E+00;
  COFD[        3419] =   0.1800620766882382E-01;
  COFD[        3420] =  -0.2061243965594226E+02;
  COFD[        3421] =   0.5108287844008140E+01;
  COFD[        3422] =  -0.4315691026948430E+00;
  COFD[        3423] =   0.1802518831107918E-01;
  COFD[        3424] =  -0.2108550712837031E+02;
  COFD[        3425] =   0.5216945829110916E+01;
  COFD[        3426] =  -0.4367926562187473E+00;
  COFD[        3427] =   0.1789880186533132E-01;
  COFD[        3428] =  -0.2091965964592144E+02;
  COFD[        3429] =   0.5148262389934279E+01;
  COFD[        3430] =  -0.4310108894089904E+00;
  COFD[        3431] =   0.1777133583281057E-01;
  COFD[        3432] =  -0.2092889793954542E+02;
  COFD[        3433] =   0.5148560230790550E+01;
  COFD[        3434] =  -0.4310562392852854E+00;
  COFD[        3435] =   0.1777358978417334E-01;
  COFD[        3436] =  -0.1946053017195549E+02;
  COFD[        3437] =   0.4857004282574310E+01;
  COFD[        3438] =  -0.4061496791697473E+00;
  COFD[        3439] =   0.1721033304344743E-01;
  COFD[        3440] =  -0.2131778720446518E+02;
  COFD[        3441] =   0.5116107656502815E+01;
  COFD[        3442] =  -0.4072201243811285E+00;
  COFD[        3443] =   0.1597605640732582E-01;
  COFD[        3444] =  -0.2131778720446518E+02;
  COFD[        3445] =   0.5116107656502815E+01;
  COFD[        3446] =  -0.4072201243811285E+00;
  COFD[        3447] =   0.1597605640732582E-01;
  COFD[        3448] =  -0.1863141958605046E+02;
  COFD[        3449] =   0.4533772315478532E+01;
  COFD[        3450] =  -0.3728642283605221E+00;
  COFD[        3451] =   0.1612572322702682E-01;
  COFD[        3452] =  -0.1944451603679417E+02;
  COFD[        3453] =   0.4785108067913721E+01;
  COFD[        3454] =  -0.3991858464067674E+00;
  COFD[        3455] =   0.1700639619802289E-01;
  COFD[        3456] =  -0.1309372900688130E+02;
  COFD[        3457] =   0.3192057362682311E+01;
  COFD[        3458] =  -0.2080508536154327E+00;
  COFD[        3459] =   0.9377275278300783E-02;
  COFD[        3460] =  -0.1695212975476949E+02;
  COFD[        3461] =   0.4856357283707361E+01;
  COFD[        3462] =  -0.4155359687875360E+00;
  COFD[        3463] =   0.1799470588389149E-01;
  COFD[        3464] =  -0.1575708126216198E+02;
  COFD[        3465] =   0.3859037941929713E+01;
  COFD[        3466] =  -0.2920671356487180E+00;
  COFD[        3467] =   0.1290449436293418E-01;
  COFD[        3468] =  -0.1687140720665078E+02;
  COFD[        3469] =   0.4075868327717259E+01;
  COFD[        3470] =  -0.3175542600319969E+00;
  COFD[        3471] =   0.1389707797540003E-01;
  COFD[        3472] =  -0.1574882548231690E+02;
  COFD[        3473] =   0.3846525250144765E+01;
  COFD[        3474] =  -0.2903598878615679E+00;
  COFD[        3475] =   0.1282703770058140E-01;
  COFD[        3476] =  -0.1993851215676389E+02;
  COFD[        3477] =   0.5098663262971057E+01;
  COFD[        3478] =  -0.4129849728390237E+00;
  COFD[        3479] =   0.1653578043838391E-01;
  COFD[        3480] =  -0.1687995772389570E+02;
  COFD[        3481] =   0.4075816169135930E+01;
  COFD[        3482] =  -0.3175458358276380E+00;
  COFD[        3483] =   0.1389664535256806E-01;
  COFD[        3484] =  -0.1688882624224305E+02;
  COFD[        3485] =   0.4076043066784620E+01;
  COFD[        3486] =  -0.3175741894503459E+00;
  COFD[        3487] =   0.1389782472898128E-01;
  COFD[        3488] =  -0.1559790150312219E+02;
  COFD[        3489] =   0.3773944780684847E+01;
  COFD[        3490] =  -0.2821502865661450E+00;
  COFD[        3491] =   0.1252380559259624E-01;
  COFD[        3492] =  -0.1578342511116079E+02;
  COFD[        3493] =   0.3901749158409890E+01;
  COFD[        3494] =  -0.2978934623317541E+00;
  COFD[        3495] =   0.1316878577801081E-01;
  COFD[        3496] =  -0.1781814842212219E+02;
  COFD[        3497] =   0.4479878049884859E+01;
  COFD[        3498] =  -0.3677558817559772E+00;
  COFD[        3499] =   0.1598132627369455E-01;
  COFD[        3500] =  -0.1781814842212219E+02;
  COFD[        3501] =   0.4479878049884859E+01;
  COFD[        3502] =  -0.3677558817559772E+00;
  COFD[        3503] =   0.1598132627369455E-01;
  COFD[        3504] =  -0.1782506140103899E+02;
  COFD[        3505] =   0.4472550110816844E+01;
  COFD[        3506] =  -0.3668215949172294E+00;
  COFD[        3507] =   0.1594178999185863E-01;
  COFD[        3508] =  -0.1775786527319250E+02;
  COFD[        3509] =   0.4446009317231327E+01;
  COFD[        3510] =  -0.3635789359718016E+00;
  COFD[        3511] =   0.1580980006557884E-01;
  COFD[        3512] =  -0.1658530159401671E+02;
  COFD[        3513] =   0.3962484009621037E+01;
  COFD[        3514] =  -0.3034260162001460E+00;
  COFD[        3515] =   0.1330924504181773E-01;
  COFD[        3516] =  -0.1934676521684402E+02;
  COFD[        3517] =   0.4844117217038898E+01;
  COFD[        3518] =  -0.4050654922829996E+00;
  COFD[        3519] =   0.1718562636894347E-01;
  COFD[        3520] =  -0.2062895610140275E+02;
  COFD[        3521] =   0.5231225460685535E+01;
  COFD[        3522] =  -0.4370160277895190E+00;
  COFD[        3523] =   0.1784652320332632E-01;
  COFD[        3524] =  -0.2063989479595861E+02;
  COFD[        3525] =   0.5231645585324166E+01;
  COFD[        3526] =  -0.4370732404462373E+00;
  COFD[        3527] =   0.1784907337224841E-01;
  COFD[        3528] =  -0.2025196923914840E+02;
  COFD[        3529] =   0.5111045666234173E+01;
  COFD[        3530] =  -0.4261812523465015E+00;
  COFD[        3531] =   0.1756072363253846E-01;
  COFD[        3532] =  -0.2025196923914840E+02;
  COFD[        3533] =   0.5111045666234173E+01;
  COFD[        3534] =  -0.4261812523465015E+00;
  COFD[        3535] =   0.1756072363253846E-01;
  COFD[        3536] =  -0.2057956090814800E+02;
  COFD[        3537] =   0.5203149310742000E+01;
  COFD[        3538] =  -0.4339608258248750E+00;
  COFD[        3539] =   0.1773356769165475E-01;
  COFD[        3540] =  -0.1891361306283375E+02;
  COFD[        3541] =   0.4728417953878362E+01;
  COFD[        3542] =  -0.3935830815508194E+00;
  COFD[        3543] =   0.1683444422800744E-01;
  COFD[        3544] =  -0.1892576390190177E+02;
  COFD[        3545] =   0.4728500364939531E+01;
  COFD[        3546] =  -0.3936077379198002E+00;
  COFD[        3547] =   0.1683614402582818E-01;
  COFD[        3548] =  -0.1893741821656993E+02;
  COFD[        3549] =   0.4728594869696905E+01;
  COFD[        3550] =  -0.3936271603294507E+00;
  COFD[        3551] =   0.1683730982264022E-01;
  COFD[        3552] =  -0.1961481660694155E+02;
  COFD[        3553] =   0.4942832081837352E+01;
  COFD[        3554] =  -0.4152110398291045E+00;
  COFD[        3555] =   0.1752038373578458E-01;
  COFD[        3556] =  -0.1945124782431341E+02;
  COFD[        3557] =   0.4857428845983466E+01;
  COFD[        3558] =  -0.4062322300213044E+00;
  COFD[        3559] =   0.1721513119766667E-01;
  COFD[        3560] =  -0.1946053017195549E+02;
  COFD[        3561] =   0.4857004282574310E+01;
  COFD[        3562] =  -0.4061496791697473E+00;
  COFD[        3563] =   0.1721033304344743E-01;
  COFD[        3564] =  -0.1753045259567557E+02;
  COFD[        3565] =   0.4405877830669391E+01;
  COFD[        3566] =  -0.3571392522269409E+00;
  COFD[        3567] =   0.1548049973130148E-01;
  COFD[        3568] =  -0.2042383193691537E+02;
  COFD[        3569] =   0.5091192695446570E+01;
  COFD[        3570] =  -0.4210911032655343E+00;
  COFD[        3571] =   0.1723542670509821E-01;
  COFD[        3572] =  -0.2042383193691537E+02;
  COFD[        3573] =   0.5091192695446570E+01;
  COFD[        3574] =  -0.4210911032655343E+00;
  COFD[        3575] =   0.1723542670509821E-01;
  COFD[        3576] =  -0.1654369229918947E+02;
  COFD[        3577] =   0.3949972451339995E+01;
  COFD[        3578] =  -0.3017790380457315E+00;
  COFD[        3579] =   0.1323685590120937E-01;
  COFD[        3580] =  -0.1767426187212671E+02;
  COFD[        3581] =   0.4360850719063046E+01;
  COFD[        3582] =  -0.3528984604085298E+00;
  COFD[        3583] =   0.1536360662305269E-01;
  COFD[        3584] =  -0.1769072264213409E+02;
  COFD[        3585] =   0.4739685996132764E+01;
  COFD[        3586] =  -0.4053730868412535E+00;
  COFD[        3587] =   0.1776343130216211E-01;
  COFD[        3588] =  -0.1510629945043326E+02;
  COFD[        3589] =   0.3409059092179077E+01;
  COFD[        3590] =  -0.1659094843453963E+00;
  COFD[        3591] =   0.4585455445383278E-02;
  COFD[        3592] =  -0.1900240647811995E+02;
  COFD[        3593] =   0.4757008879827381E+01;
  COFD[        3594] =  -0.3907894080963009E+00;
  COFD[        3595] =   0.1643024015529323E-01;
  COFD[        3596] =  -0.2009809718705050E+02;
  COFD[        3597] =   0.4977640189149046E+01;
  COFD[        3598] =  -0.4165323971276452E+00;
  COFD[        3599] =   0.1744635825561870E-01;
  COFD[        3600] =  -0.1902751852931795E+02;
  COFD[        3601] =   0.4760378344803059E+01;
  COFD[        3602] =  -0.3915558236239631E+00;
  COFD[        3603] =   0.1647818289394700E-01;
  COFD[        3604] =  -0.1836390465233753E+02;
  COFD[        3605] =   0.3824966367676553E+01;
  COFD[        3606] =  -0.1992927348787672E+00;
  COFD[        3607] =   0.5478606239007274E-02;
  COFD[        3608] =  -0.2011350568921456E+02;
  COFD[        3609] =   0.4980977035855766E+01;
  COFD[        3610] =  -0.4170734091977097E+00;
  COFD[        3611] =   0.1747455638443835E-01;
  COFD[        3612] =  -0.2012780676364666E+02;
  COFD[        3613] =   0.4983954385219352E+01;
  COFD[        3614] =  -0.4175569796553835E+00;
  COFD[        3615] =   0.1749979298053702E-01;
  COFD[        3616] =  -0.1894980246607794E+02;
  COFD[        3617] =   0.4748233808709042E+01;
  COFD[        3618] =  -0.3930779606458064E+00;
  COFD[        3619] =   0.1667238331348419E-01;
  COFD[        3620] =  -0.1891456819097078E+02;
  COFD[        3621] =   0.4745553725181649E+01;
  COFD[        3622] =  -0.3882051404698895E+00;
  COFD[        3623] =   0.1626915360804150E-01;
  COFD[        3624] =  -0.1997313604118449E+02;
  COFD[        3625] =   0.4888134291131157E+01;
  COFD[        3626] =  -0.3919147663317576E+00;
  COFD[        3627] =   0.1584552869833823E-01;
  COFD[        3628] =  -0.1997313604118449E+02;
  COFD[        3629] =   0.4888134291131157E+01;
  COFD[        3630] =  -0.3919147663317576E+00;
  COFD[        3631] =   0.1584552869833823E-01;
  COFD[        3632] =  -0.2004560573876227E+02;
  COFD[        3633] =   0.4909854771274889E+01;
  COFD[        3634] =  -0.3952294049921926E+00;
  COFD[        3635] =   0.1601054023046091E-01;
  COFD[        3636] =  -0.2008186085747107E+02;
  COFD[        3637] =   0.4928935175864744E+01;
  COFD[        3638] =  -0.3987527366535738E+00;
  COFD[        3639] =   0.1620682050200429E-01;
  COFD[        3640] =  -0.1996446349112875E+02;
  COFD[        3641] =   0.4943053687420567E+01;
  COFD[        3642] =  -0.4143586003201034E+00;
  COFD[        3643] =   0.1744554338516522E-01;
  COFD[        3644] =  -0.2140522217493100E+02;
  COFD[        3645] =   0.5181866929752220E+01;
  COFD[        3646] =  -0.4179647553808900E+00;
  COFD[        3647] =   0.1653089926470964E-01;
  COFD[        3648] =  -0.1972933562163679E+02;
  COFD[        3649] =   0.4264735345000579E+01;
  COFD[        3650] =  -0.2639983274633736E+00;
  COFD[        3651] =   0.8526026392458589E-02;
  COFD[        3652] =  -0.1974668857323474E+02;
  COFD[        3653] =   0.4267776444888849E+01;
  COFD[        3654] =  -0.2644050809170580E+00;
  COFD[        3655] =   0.8543267577410659E-02;
  COFD[        3656] =  -0.2039593029382244E+02;
  COFD[        3657] =   0.4595576419324821E+01;
  COFD[        3658] =  -0.3164280753520152E+00;
  COFD[        3659] =   0.1115468584869250E-01;
  COFD[        3660] =  -0.2039593029382244E+02;
  COFD[        3661] =   0.4595576419324821E+01;
  COFD[        3662] =  -0.3164280753520152E+00;
  COFD[        3663] =   0.1115468584869250E-01;
  COFD[        3664] =  -0.1994912980553432E+02;
  COFD[        3665] =   0.4354262271960964E+01;
  COFD[        3666] =  -0.2776328852684010E+00;
  COFD[        3667] =   0.9191073490857130E-02;
  COFD[        3668] =  -0.2106259252673142E+02;
  COFD[        3669] =   0.5117958864455729E+01;
  COFD[        3670] =  -0.4136130650204043E+00;
  COFD[        3671] =   0.1648627810521946E-01;
  COFD[        3672] =  -0.2109960034604860E+02;
  COFD[        3673] =   0.5129133288825352E+01;
  COFD[        3674] =  -0.4152735461892595E+00;
  COFD[        3675] =   0.1656709485798127E-01;
  COFD[        3676] =  -0.2113363149383363E+02;
  COFD[        3677] =   0.5139253110457338E+01;
  COFD[        3678] =  -0.4167775115779367E+00;
  COFD[        3679] =   0.1664030318824075E-01;
  COFD[        3680] =  -0.2115552988307172E+02;
  COFD[        3681] =   0.5049408694314179E+01;
  COFD[        3682] =  -0.3939342769912307E+00;
  COFD[        3683] =   0.1522309753276939E-01;
  COFD[        3684] =  -0.2129367686674487E+02;
  COFD[        3685] =   0.5109799100487325E+01;
  COFD[        3686] =  -0.4062923797917856E+00;
  COFD[        3687] =   0.1593137641849987E-01;
  COFD[        3688] =  -0.2131778720446518E+02;
  COFD[        3689] =   0.5116107656502815E+01;
  COFD[        3690] =  -0.4072201243811285E+00;
  COFD[        3691] =   0.1597605640732582E-01;
  COFD[        3692] =  -0.2042383193691537E+02;
  COFD[        3693] =   0.5091192695446570E+01;
  COFD[        3694] =  -0.4210911032655343E+00;
  COFD[        3695] =   0.1723542670509821E-01;
  COFD[        3696] =  -0.2047887601043758E+02;
  COFD[        3697] =   0.4549054655103817E+01;
  COFD[        3698] =  -0.3083705960322820E+00;
  COFD[        3699] =   0.1072673844162708E-01;
  COFD[        3700] =  -0.2047887601043758E+02;
  COFD[        3701] =   0.4549054655103817E+01;
  COFD[        3702] =  -0.3083705960322820E+00;
  COFD[        3703] =   0.1072673844162708E-01;
  COFD[        3704] =  -0.1994082237937651E+02;
  COFD[        3705] =   0.4937902067410034E+01;
  COFD[        3706] =  -0.4137962628355027E+00;
  COFD[        3707] =   0.1742541544956688E-01;
  COFD[        3708] =  -0.2061978339529828E+02;
  COFD[        3709] =   0.5111890805302203E+01;
  COFD[        3710] =  -0.4278633892785161E+00;
  COFD[        3711] =   0.1769525721185838E-01;
  COFD[        3712] =  -0.1769072264213409E+02;
  COFD[        3713] =   0.4739685996132764E+01;
  COFD[        3714] =  -0.4053730868412535E+00;
  COFD[        3715] =   0.1776343130216211E-01;
  COFD[        3716] =  -0.1510629945043326E+02;
  COFD[        3717] =   0.3409059092179077E+01;
  COFD[        3718] =  -0.1659094843453963E+00;
  COFD[        3719] =   0.4585455445383278E-02;
  COFD[        3720] =  -0.1900240647811995E+02;
  COFD[        3721] =   0.4757008879827381E+01;
  COFD[        3722] =  -0.3907894080963009E+00;
  COFD[        3723] =   0.1643024015529323E-01;
  COFD[        3724] =  -0.2009809718705050E+02;
  COFD[        3725] =   0.4977640189149046E+01;
  COFD[        3726] =  -0.4165323971276452E+00;
  COFD[        3727] =   0.1744635825561870E-01;
  COFD[        3728] =  -0.1902751852931795E+02;
  COFD[        3729] =   0.4760378344803059E+01;
  COFD[        3730] =  -0.3915558236239631E+00;
  COFD[        3731] =   0.1647818289394700E-01;
  COFD[        3732] =  -0.1836390465233753E+02;
  COFD[        3733] =   0.3824966367676553E+01;
  COFD[        3734] =  -0.1992927348787672E+00;
  COFD[        3735] =   0.5478606239007274E-02;
  COFD[        3736] =  -0.2011350568921456E+02;
  COFD[        3737] =   0.4980977035855766E+01;
  COFD[        3738] =  -0.4170734091977097E+00;
  COFD[        3739] =   0.1747455638443835E-01;
  COFD[        3740] =  -0.2012780676364666E+02;
  COFD[        3741] =   0.4983954385219352E+01;
  COFD[        3742] =  -0.4175569796553835E+00;
  COFD[        3743] =   0.1749979298053702E-01;
  COFD[        3744] =  -0.1894980246607794E+02;
  COFD[        3745] =   0.4748233808709042E+01;
  COFD[        3746] =  -0.3930779606458064E+00;
  COFD[        3747] =   0.1667238331348419E-01;
  COFD[        3748] =  -0.1891456819097078E+02;
  COFD[        3749] =   0.4745553725181649E+01;
  COFD[        3750] =  -0.3882051404698895E+00;
  COFD[        3751] =   0.1626915360804150E-01;
  COFD[        3752] =  -0.1997313604118449E+02;
  COFD[        3753] =   0.4888134291131157E+01;
  COFD[        3754] =  -0.3919147663317576E+00;
  COFD[        3755] =   0.1584552869833823E-01;
  COFD[        3756] =  -0.1997313604118449E+02;
  COFD[        3757] =   0.4888134291131157E+01;
  COFD[        3758] =  -0.3919147663317576E+00;
  COFD[        3759] =   0.1584552869833823E-01;
  COFD[        3760] =  -0.2004560573876227E+02;
  COFD[        3761] =   0.4909854771274889E+01;
  COFD[        3762] =  -0.3952294049921926E+00;
  COFD[        3763] =   0.1601054023046091E-01;
  COFD[        3764] =  -0.2008186085747107E+02;
  COFD[        3765] =   0.4928935175864744E+01;
  COFD[        3766] =  -0.3987527366535738E+00;
  COFD[        3767] =   0.1620682050200429E-01;
  COFD[        3768] =  -0.1996446349112875E+02;
  COFD[        3769] =   0.4943053687420567E+01;
  COFD[        3770] =  -0.4143586003201034E+00;
  COFD[        3771] =   0.1744554338516522E-01;
  COFD[        3772] =  -0.2140522217493100E+02;
  COFD[        3773] =   0.5181866929752220E+01;
  COFD[        3774] =  -0.4179647553808900E+00;
  COFD[        3775] =   0.1653089926470964E-01;
  COFD[        3776] =  -0.1972933562163679E+02;
  COFD[        3777] =   0.4264735345000579E+01;
  COFD[        3778] =  -0.2639983274633736E+00;
  COFD[        3779] =   0.8526026392458589E-02;
  COFD[        3780] =  -0.1974668857323474E+02;
  COFD[        3781] =   0.4267776444888849E+01;
  COFD[        3782] =  -0.2644050809170580E+00;
  COFD[        3783] =   0.8543267577410659E-02;
  COFD[        3784] =  -0.2039593029382244E+02;
  COFD[        3785] =   0.4595576419324821E+01;
  COFD[        3786] =  -0.3164280753520152E+00;
  COFD[        3787] =   0.1115468584869250E-01;
  COFD[        3788] =  -0.2039593029382244E+02;
  COFD[        3789] =   0.4595576419324821E+01;
  COFD[        3790] =  -0.3164280753520152E+00;
  COFD[        3791] =   0.1115468584869250E-01;
  COFD[        3792] =  -0.1994912980553432E+02;
  COFD[        3793] =   0.4354262271960964E+01;
  COFD[        3794] =  -0.2776328852684010E+00;
  COFD[        3795] =   0.9191073490857130E-02;
  COFD[        3796] =  -0.2106259252673142E+02;
  COFD[        3797] =   0.5117958864455729E+01;
  COFD[        3798] =  -0.4136130650204043E+00;
  COFD[        3799] =   0.1648627810521946E-01;
  COFD[        3800] =  -0.2109960034604860E+02;
  COFD[        3801] =   0.5129133288825352E+01;
  COFD[        3802] =  -0.4152735461892595E+00;
  COFD[        3803] =   0.1656709485798127E-01;
  COFD[        3804] =  -0.2113363149383363E+02;
  COFD[        3805] =   0.5139253110457338E+01;
  COFD[        3806] =  -0.4167775115779367E+00;
  COFD[        3807] =   0.1664030318824075E-01;
  COFD[        3808] =  -0.2115552988307172E+02;
  COFD[        3809] =   0.5049408694314179E+01;
  COFD[        3810] =  -0.3939342769912307E+00;
  COFD[        3811] =   0.1522309753276939E-01;
  COFD[        3812] =  -0.2129367686674487E+02;
  COFD[        3813] =   0.5109799100487325E+01;
  COFD[        3814] =  -0.4062923797917856E+00;
  COFD[        3815] =   0.1593137641849987E-01;
  COFD[        3816] =  -0.2131778720446518E+02;
  COFD[        3817] =   0.5116107656502815E+01;
  COFD[        3818] =  -0.4072201243811285E+00;
  COFD[        3819] =   0.1597605640732582E-01;
  COFD[        3820] =  -0.2042383193691537E+02;
  COFD[        3821] =   0.5091192695446570E+01;
  COFD[        3822] =  -0.4210911032655343E+00;
  COFD[        3823] =   0.1723542670509821E-01;
  COFD[        3824] =  -0.2047887601043758E+02;
  COFD[        3825] =   0.4549054655103817E+01;
  COFD[        3826] =  -0.3083705960322820E+00;
  COFD[        3827] =   0.1072673844162708E-01;
  COFD[        3828] =  -0.2047887601043758E+02;
  COFD[        3829] =   0.4549054655103817E+01;
  COFD[        3830] =  -0.3083705960322820E+00;
  COFD[        3831] =   0.1072673844162708E-01;
  COFD[        3832] =  -0.1994082237937651E+02;
  COFD[        3833] =   0.4937902067410034E+01;
  COFD[        3834] =  -0.4137962628355027E+00;
  COFD[        3835] =   0.1742541544956688E-01;
  COFD[        3836] =  -0.2061978339529828E+02;
  COFD[        3837] =   0.5111890805302203E+01;
  COFD[        3838] =  -0.4278633892785161E+00;
  COFD[        3839] =   0.1769525721185838E-01;
  COFD[        3840] =  -0.1266066516595179E+02;
  COFD[        3841] =   0.2898076057146818E+01;
  COFD[        3842] =  -0.1700453759318706E+00;
  COFD[        3843] =   0.7738690296212197E-02;
  COFD[        3844] =  -0.1711567017943538E+02;
  COFD[        3845] =   0.4833264586465570E+01;
  COFD[        3846] =  -0.4205908440542898E+00;
  COFD[        3847] =   0.1855345683322859E-01;
  COFD[        3848] =  -0.1474819919893963E+02;
  COFD[        3849] =   0.3361311502126538E+01;
  COFD[        3850] =  -0.2285465363406641E+00;
  COFD[        3851] =   0.1019990809984005E-01;
  COFD[        3852] =  -0.1578873902175925E+02;
  COFD[        3853] =   0.3589258101555060E+01;
  COFD[        3854] =  -0.2568718883472193E+00;
  COFD[        3855] =   0.1137271971440306E-01;
  COFD[        3856] =  -0.1473661449825940E+02;
  COFD[        3857] =   0.3348204558833826E+01;
  COFD[        3858] =  -0.2267271723233180E+00;
  COFD[        3859] =   0.1011600240359858E-01;
  COFD[        3860] =  -0.2024331300258610E+02;
  COFD[        3861] =   0.5167058816761965E+01;
  COFD[        3862] =  -0.4292618020154789E+00;
  COFD[        3863] =   0.1752039422653532E-01;
  COFD[        3864] =  -0.1579929195048582E+02;
  COFD[        3865] =   0.3590679879226839E+01;
  COFD[        3866] =  -0.2570681340587627E+00;
  COFD[        3867] =   0.1138172060899465E-01;
  COFD[        3868] =  -0.1581023736073806E+02;
  COFD[        3869] =   0.3592405682636056E+01;
  COFD[        3870] =  -0.2573062871053550E+00;
  COFD[        3871] =   0.1139264138769127E-01;
  COFD[        3872] =  -0.1449021183214253E+02;
  COFD[        3873] =   0.3236447428837034E+01;
  COFD[        3874] =  -0.2121618734582983E+00;
  COFD[        3875] =   0.9482346564680543E-02;
  COFD[        3876] =  -0.1479325303072075E+02;
  COFD[        3877] =   0.3409403744671862E+01;
  COFD[        3878] =  -0.2352215543219259E+00;
  COFD[        3879] =   0.1050772674323951E-01;
  COFD[        3880] =  -0.1658740614019908E+02;
  COFD[        3881] =   0.3933646317245189E+01;
  COFD[        3882] =  -0.2999344650396615E+00;
  COFD[        3883] =   0.1316797448289397E-01;
  COFD[        3884] =  -0.1658740614019908E+02;
  COFD[        3885] =   0.3933646317245189E+01;
  COFD[        3886] =  -0.2999344650396615E+00;
  COFD[        3887] =   0.1316797448289397E-01;
  COFD[        3888] =  -0.1658229952875101E+02;
  COFD[        3889] =   0.3922184822536529E+01;
  COFD[        3890] =  -0.2983959485115160E+00;
  COFD[        3891] =   0.1309927190981370E-01;
  COFD[        3892] =  -0.1651845590523054E+02;
  COFD[        3893] =   0.3896132878391479E+01;
  COFD[        3894] =  -0.2951170694006329E+00;
  COFD[        3895] =   0.1296170338310411E-01;
  COFD[        3896] =  -0.1549017806948574E+02;
  COFD[        3897] =   0.3466859355125401E+01;
  COFD[        3898] =  -0.2408856054712291E+00;
  COFD[        3899] =   0.1067561900190753E-01;
  COFD[        3900] =  -0.1850244542172855E+02;
  COFD[        3901] =   0.4509546506902997E+01;
  COFD[        3902] =  -0.3699250247632380E+00;
  COFD[        3903] =   0.1600569104209367E-01;
  COFD[        3904] =  -0.2006577079410498E+02;
  COFD[        3905] =   0.5032536447034316E+01;
  COFD[        3906] =  -0.4235925902580250E+00;
  COFD[        3907] =   0.1775279409683131E-01;
  COFD[        3908] =  -0.2007323696432831E+02;
  COFD[        3909] =   0.5032033200707714E+01;
  COFD[        3910] =  -0.4235009227577490E+00;
  COFD[        3911] =   0.1774762291235686E-01;
  COFD[        3912] =  -0.2000088490858467E+02;
  COFD[        3913] =   0.5021789747731799E+01;
  COFD[        3914] =  -0.4253480595637520E+00;
  COFD[        3915] =   0.1795667045275812E-01;
  COFD[        3916] =  -0.2000088490858467E+02;
  COFD[        3917] =   0.5021789747731799E+01;
  COFD[        3918] =  -0.4253480595637520E+00;
  COFD[        3919] =   0.1795667045275812E-01;
  COFD[        3920] =  -0.1999915273024637E+02;
  COFD[        3921] =   0.5000095313030963E+01;
  COFD[        3922] =  -0.4200303081186300E+00;
  COFD[        3923] =   0.1762406687147381E-01;
  COFD[        3924] =  -0.1797950802803312E+02;
  COFD[        3925] =   0.4344100406483967E+01;
  COFD[        3926] =  -0.3507611460911521E+00;
  COFD[        3927] =   0.1527276328163101E-01;
  COFD[        3928] =  -0.1799101249332176E+02;
  COFD[        3929] =   0.4344595334851394E+01;
  COFD[        3930] =  -0.3508294711072030E+00;
  COFD[        3931] =   0.1527589170814564E-01;
  COFD[        3932] =  -0.1800289314001846E+02;
  COFD[        3933] =   0.4345462384297714E+01;
  COFD[        3934] =  -0.3509451050542303E+00;
  COFD[        3935] =   0.1528101631398434E-01;
  COFD[        3936] =  -0.1876345278233955E+02;
  COFD[        3937] =   0.4606235123260915E+01;
  COFD[        3938] =  -0.3803641770348614E+00;
  COFD[        3939] =   0.1637150008428244E-01;
  COFD[        3940] =  -0.1861866210245701E+02;
  COFD[        3941] =   0.4532060302179115E+01;
  COFD[        3942] =  -0.3726516669477872E+00;
  COFD[        3943] =   0.1611699967368039E-01;
  COFD[        3944] =  -0.1863141958605046E+02;
  COFD[        3945] =   0.4533772315478532E+01;
  COFD[        3946] =  -0.3728642283605221E+00;
  COFD[        3947] =   0.1612572322702682E-01;
  COFD[        3948] =  -0.1654369229918947E+02;
  COFD[        3949] =   0.3949972451339995E+01;
  COFD[        3950] =  -0.3017790380457315E+00;
  COFD[        3951] =   0.1323685590120937E-01;
  COFD[        3952] =  -0.1994082237937651E+02;
  COFD[        3953] =   0.4937902067410034E+01;
  COFD[        3954] =  -0.4137962628355027E+00;
  COFD[        3955] =   0.1742541544956688E-01;
  COFD[        3956] =  -0.1994082237937651E+02;
  COFD[        3957] =   0.4937902067410034E+01;
  COFD[        3958] =  -0.4137962628355027E+00;
  COFD[        3959] =   0.1742541544956688E-01;
  COFD[        3960] =  -0.1546612761154138E+02;
  COFD[        3961] =   0.3460858520880728E+01;
  COFD[        3962] =  -0.2401164641793465E+00;
  COFD[        3963] =   0.1064270570979806E-01;
  COFD[        3964] =  -0.1648483159780408E+02;
  COFD[        3965] =   0.3836871316958563E+01;
  COFD[        3966] =  -0.2875628541513582E+00;
  COFD[        3967] =   0.1264044562964702E-01;
  COFD[        3968] =  -0.1342264266757377E+02;
  COFD[        3969] =   0.3221777556990338E+01;
  COFD[        3970] =  -0.2128673394230190E+00;
  COFD[        3971] =   0.9627787551744238E-02;
  COFD[        3972] =  -0.1752038279713604E+02;
  COFD[        3973] =   0.4959639086677385E+01;
  COFD[        3974] =  -0.4293974875408633E+00;
  COFD[        3975] =   0.1860806087249651E-01;
  COFD[        3976] =  -0.1581858256688188E+02;
  COFD[        3977] =   0.3775899039344615E+01;
  COFD[        3978] =  -0.2814687920927101E+00;
  COFD[        3979] =   0.1245218573218981E-01;
  COFD[        3980] =  -0.1677455230511400E+02;
  COFD[        3981] =   0.3946214344521416E+01;
  COFD[        3982] =  -0.3012613112287137E+00;
  COFD[        3983] =   0.1321325422150532E-01;
  COFD[        3984] =  -0.1580533640016843E+02;
  COFD[        3985] =   0.3761415063014861E+01;
  COFD[        3986] =  -0.2795015971789234E+00;
  COFD[        3987] =   0.1236332420639666E-01;
  COFD[        3988] =  -0.1991427613539219E+02;
  COFD[        3989] =   0.4985074718206674E+01;
  COFD[        3990] =  -0.3988046631693363E+00;
  COFD[        3991] =   0.1592993698509921E-01;
  COFD[        3992] =  -0.1677798202265061E+02;
  COFD[        3993] =   0.3944105793182093E+01;
  COFD[        3994] =  -0.3009768138615855E+00;
  COFD[        3995] =   0.1320048849753941E-01;
  COFD[        3996] =  -0.1678197080138354E+02;
  COFD[        3997] =   0.3942379447222446E+01;
  COFD[        3998] =  -0.3007437798097645E+00;
  COFD[        3999] =   0.1319002761262226E-01;
  COFD[        4000] =  -0.1559892478263072E+02;
  COFD[        4001] =   0.3675877348865438E+01;
  COFD[        4002] =  -0.2694954491084075E+00;
  COFD[        4003] =   0.1197817600584641E-01;
  COFD[        4004] =  -0.1586061989522203E+02;
  COFD[        4005] =   0.3824799145446967E+01;
  COFD[        4006] =  -0.2881091362509902E+00;
  COFD[        4007] =   0.1275209876189350E-01;
  COFD[        4008] =  -0.1779539137381903E+02;
  COFD[        4009] =   0.4385582017594726E+01;
  COFD[        4010] =  -0.3562391667708982E+00;
  COFD[        4011] =   0.1551072662807571E-01;
  COFD[        4012] =  -0.1779539137381903E+02;
  COFD[        4013] =   0.4385582017594726E+01;
  COFD[        4014] =  -0.3562391667708982E+00;
  COFD[        4015] =   0.1551072662807571E-01;
  COFD[        4016] =  -0.1779710613097560E+02;
  COFD[        4017] =   0.4376269742117769E+01;
  COFD[        4018] =  -0.3550500828076509E+00;
  COFD[        4019] =   0.1546030915508793E-01;
  COFD[        4020] =  -0.1772757250170383E+02;
  COFD[        4021] =   0.4347725932428951E+01;
  COFD[        4022] =  -0.3515015798062593E+00;
  COFD[        4023] =   0.1531308129463464E-01;
  COFD[        4024] =  -0.1650781936547825E+02;
  COFD[        4025] =   0.3842122413055806E+01;
  COFD[        4026] =  -0.2882098997185555E+00;
  COFD[        4027] =   0.1266705263344362E-01;
  COFD[        4028] =  -0.1931016351705685E+02;
  COFD[        4029] =   0.4757139666603682E+01;
  COFD[        4030] =  -0.3961700036764624E+00;
  COFD[        4031] =   0.1690081368592412E-01;
  COFD[        4032] =  -0.2070496369832456E+02;
  COFD[        4033] =   0.5191072368701859E+01;
  COFD[        4034] =  -0.4346910627906122E+00;
  COFD[        4035] =   0.1785842811134005E-01;
  COFD[        4036] =  -0.2072332061257778E+02;
  COFD[        4037] =   0.5194960177957292E+01;
  COFD[        4038] =  -0.4352726149968372E+00;
  COFD[        4039] =   0.1788688372961993E-01;
  COFD[        4040] =  -0.2043431836412989E+02;
  COFD[        4041] =   0.5116486808224121E+01;
  COFD[        4042] =  -0.4301998553589693E+00;
  COFD[        4043] =   0.1786838693497468E-01;
  COFD[        4044] =  -0.2043431836412989E+02;
  COFD[        4045] =   0.5116486808224121E+01;
  COFD[        4046] =  -0.4301998553589693E+00;
  COFD[        4047] =   0.1786838693497468E-01;
  COFD[        4048] =  -0.2065593872227729E+02;
  COFD[        4049] =   0.5165539551318380E+01;
  COFD[        4050] =  -0.4322266808318962E+00;
  COFD[        4051] =   0.1778317395017445E-01;
  COFD[        4052] =  -0.1893874700238532E+02;
  COFD[        4053] =   0.4666666006155531E+01;
  COFD[        4054] =  -0.3875266304039654E+00;
  COFD[        4055] =   0.1665285657108702E-01;
  COFD[        4056] =  -0.1894826439081699E+02;
  COFD[        4057] =   0.4665904241842780E+01;
  COFD[        4058] =  -0.3874653193184607E+00;
  COFD[        4059] =   0.1665191589170439E-01;
  COFD[        4060] =  -0.1895758159603200E+02;
  COFD[        4061] =   0.4665279128753897E+01;
  COFD[        4062] =  -0.3874162488122672E+00;
  COFD[        4063] =   0.1665125384367453E-01;
  COFD[        4064] =  -0.1962450162993201E+02;
  COFD[        4065] =   0.4875751864229406E+01;
  COFD[        4066] =  -0.4087050411778894E+00;
  COFD[        4067] =   0.1732743672805724E-01;
  COFD[        4068] =  -0.1943393371703941E+02;
  COFD[        4069] =   0.4784664368892784E+01;
  COFD[        4070] =  -0.3991077435475853E+00;
  COFD[        4071] =   0.1700210403998681E-01;
  COFD[        4072] =  -0.1944451603679417E+02;
  COFD[        4073] =   0.4785108067913721E+01;
  COFD[        4074] =  -0.3991858464067674E+00;
  COFD[        4075] =   0.1700639619802289E-01;
  COFD[        4076] =  -0.1767426187212671E+02;
  COFD[        4077] =   0.4360850719063046E+01;
  COFD[        4078] =  -0.3528984604085298E+00;
  COFD[        4079] =   0.1536360662305269E-01;
  COFD[        4080] =  -0.2061978339529828E+02;
  COFD[        4081] =   0.5111890805302203E+01;
  COFD[        4082] =  -0.4278633892785161E+00;
  COFD[        4083] =   0.1769525721185838E-01;
  COFD[        4084] =  -0.2061978339529828E+02;
  COFD[        4085] =   0.5111890805302203E+01;
  COFD[        4086] =  -0.4278633892785161E+00;
  COFD[        4087] =   0.1769525721185838E-01;
  COFD[        4088] =  -0.1648483159780408E+02;
  COFD[        4089] =   0.3836871316958563E+01;
  COFD[        4090] =  -0.2875628541513582E+00;
  COFD[        4091] =   0.1264044562964702E-01;
  COFD[        4092] =  -0.1771336337770503E+02;
  COFD[        4093] =   0.4289097727284545E+01;
  COFD[        4094] =  -0.3450079724750131E+00;
  COFD[        4095] =   0.1508113614518283E-01;
};  }




#if 0




\\
\\
\\  This is the mechanism file
\\
\\
! GRI-Mech Version 1.2  released 11/16/94   CHEMKIN-II format
! See README file at anonymous FTP site GRI_MECH at CRVAX.SRI.COM or
!   WorldWideWeb home page at http://diesel.fsc.psu.edu/~gri_mech
!      for additional information, contacts, and disclaimer.
ELEMENTS                                                                        
O  H  C  N  AR                                                                  
END                                                                             
SPECIES                                                                         
H2      H       O       O2      OH      H2O     HO2     H2O2                    
C       CH      CH2     CH2(S)  CH3     CH4     CO      CO2                     
HCO     CH2O    CH2OH   CH3O    CH3OH   C2H     C2H2    C2H3                    
C2H4    C2H5    C2H6    HCCO    CH2CO   HCCOH   N2      AR                      
END                                                                             
! THERMO
! Insert GRI-Mech thermodynamics here or use in default file
! END
REACTIONS                                                                       
2O+M<=>O2+M                              1.200E+17   -1.000      0.00           
H2/ 2.40/ H2O/15.40/ CH4/ 2.00/ CO/ 1.75/ CO2/ 3.60/ C2H6/ 3.00/ AR/ 0.83/      
O+H+M<=>OH+M                             5.000E+17   -1.000      0.00           
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
O+H2<=>H+OH                              5.000E+04    2.670   6290.00           
O+HO2<=>OH+O2                            2.000E+13    0.000      0.00           
O+H2O2<=>OH+HO2                          9.630E+06    2.000   4000.00           
O+CH<=>H+CO                              5.700E+13    0.000      0.00           
O+CH2<=>H+HCO                            8.000E+13    0.000      0.00           
O+CH2(S)<=>H2+CO                         1.500E+13    0.000      0.00           
O+CH2(S)<=>H+HCO                         1.500E+13    0.000      0.00           
O+CH3<=>H+CH2O                           8.430E+13    0.000      0.00           
O+CH4<=>OH+CH3                           1.020E+09    1.500   8600.00           
O+CO+M<=>CO2+M                           6.020E+14    0.000   3000.00           
H2/2.00/ O2/6.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/3.50/ C2H6/3.00/ AR/0.50/    
O+HCO<=>OH+CO                            3.000E+13    0.000      0.00           
O+HCO<=>H+CO2                            3.000E+13    0.000      0.00           
O+CH2O<=>OH+HCO                          3.900E+13    0.000   3540.00           
O+CH2OH<=>OH+CH2O                        1.000E+13    0.000      0.00           
O+CH3O<=>OH+CH2O                         1.000E+13    0.000      0.00           
O+CH3OH<=>OH+CH2OH                       3.880E+05    2.500   3100.00           
O+CH3OH<=>OH+CH3O                        1.300E+05    2.500   5000.00           
O+C2H<=>CH+CO                            5.000E+13    0.000      0.00           
O+C2H2<=>H+HCCO                          1.020E+07    2.000   1900.00           
O+C2H2<=>OH+C2H                          4.600E+19   -1.410  28950.00           
O+C2H2<=>CO+CH2                          1.020E+07    2.000   1900.00           
O+C2H3<=>H+CH2CO                         3.000E+13    0.000      0.00           
O+C2H4<=>CH3+HCO                         1.920E+07    1.830    220.00           
O+C2H5<=>CH3+CH2O                        1.320E+14    0.000      0.00           
O+C2H6<=>OH+C2H5                         8.980E+07    1.920   5690.00           
O+HCCO<=>H+2CO                           1.000E+14    0.000      0.00           
O+CH2CO<=>OH+HCCO                        1.000E+13    0.000   8000.00           
O+CH2CO<=>CH2+CO2                        1.750E+12    0.000   1350.00           
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
H+HO2<=>O+H2O                            3.970E+12    0.000    671.00           
H+HO2<=>O2+H2                            2.800E+13    0.000   1068.00           
H+HO2<=>2OH                              1.340E+14    0.000    635.00           
H+H2O2<=>HO2+H2                          1.210E+07    2.000   5200.00           
H+H2O2<=>OH+H2O                          1.000E+13    0.000   3600.00           
H+CH<=>C+H2                              1.100E+14    0.000      0.00           
H+CH2(+M)<=>CH3(+M)                      2.500E+16   -0.800      0.00           
     LOW  /  3.200E+27   -3.140   1230.00/                                      
     TROE/  0.6800   78.00  1995.00  5590.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
H+CH2(S)<=>CH+H2                         3.000E+13    0.000      0.00           
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
H+CH2O(+M)<=>CH2OH(+M)                   5.400E+11    0.454   3600.00           
     LOW  /  1.270E+32   -4.820   6530.00/                                      
     TROE/  0.7187  103.00  1291.00  4160.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/                      
H+CH2O(+M)<=>CH3O(+M)                    5.400E+11    0.454   2600.00           
     LOW  /  2.200E+30   -4.800   5560.00/                                      
     TROE/  0.7580   94.00  1555.00  4200.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/                      
H+CH2O<=>HCO+H2                          2.300E+10    1.050   3275.00           
H+CH2OH(+M)<=>CH3OH(+M)                  1.800E+13    0.000      0.00           
     LOW  /  3.000E+31   -4.800   3300.00/                                      
     TROE/  0.7679  338.00  1812.00  5081.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/                      
H+CH2OH<=>H2+CH2O                        2.000E+13    0.000      0.00           
H+CH2OH<=>OH+CH3                         1.200E+13    0.000      0.00           
H+CH2OH<=>CH2(S)+H2O                     6.000E+12    0.000      0.00           
H+CH3O(+M)<=>CH3OH(+M)                   5.000E+13    0.000      0.00           
     LOW  /  8.600E+28   -4.000   3025.00/                                      
     TROE/  0.8902  144.00  2838.00 45569.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/                      
H+CH3O<=>H+CH2OH                         3.400E+06    1.600      0.00           
H+CH3O<=>H2+CH2O                         2.000E+13    0.000      0.00           
H+CH3O<=>OH+CH3                          3.200E+13    0.000      0.00           
H+CH3O<=>CH2(S)+H2O                      1.600E+13    0.000      0.00           
H+CH3OH<=>CH2OH+H2                       1.700E+07    2.100   4870.00           
H+CH3OH<=>CH3O+H2                        4.200E+06    2.100   4870.00           
H+C2H(+M)<=>C2H2(+M)                     1.000E+17   -1.000      0.00           
     LOW  /  3.750E+33   -4.800   1900.00/                                      
     TROE/  0.6464  132.00  1315.00  5566.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
H+C2H2(+M)<=>C2H3(+M)                    5.600E+12    0.000   2400.00           
     LOW  /  3.800E+40   -7.270   7220.00/                                      
     TROE/  0.7507   98.50  1302.00  4167.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
H+C2H3(+M)<=>C2H4(+M)                    6.080E+12    0.270    280.00           
     LOW  /  1.400E+30   -3.860   3320.00/                                      
     TROE/  0.7820  207.50  2663.00  6095.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
H+C2H3<=>H2+C2H2                         3.000E+13    0.000      0.00           
H+C2H4(+M)<=>C2H5(+M)                    1.080E+12    0.454   1820.00           
     LOW  /  1.200E+42   -7.620   6970.00/                                      
     TROE/  0.9753  210.00   984.00  4374.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
H+C2H4<=>C2H3+H2                         1.325E+06    2.530  12240.00           
H+C2H5(+M)<=>C2H6(+M)                    5.210E+17   -0.990   1580.00           
     LOW  /  1.990E+41   -7.080   6685.00/                                      
     TROE/  0.8422  125.00  2219.00  6882.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
H+C2H5<=>H2+C2H4                         2.000E+12    0.000      0.00           
H+C2H6<=>C2H5+H2                         1.150E+08    1.900   7530.00           
H+HCCO<=>CH2(S)+CO                       1.000E+14    0.000      0.00           
H+CH2CO<=>HCCO+H2                        5.000E+13    0.000   8000.00           
H+CH2CO<=>CH3+CO                         1.130E+13    0.000   3428.00           
H+HCCOH<=>H+CH2CO                        1.000E+13    0.000      0.00           
H2+CO(+M)<=>CH2O(+M)                     4.300E+07    1.500  79600.00           
     LOW  /  5.070E+27   -3.420  84350.00/                                      
     TROE/  0.9320  197.00  1540.00 10300.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
OH+H2<=>H+H2O                            2.160E+08    1.510   3430.00           
2OH(+M)<=>H2O2(+M)                       7.400E+13   -0.370      0.00           
     LOW  /  2.300E+18   -0.900  -1700.00/                                      
     TROE/  0.7346   94.00  1756.00  5182.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
2OH<=>O+H2O                              3.570E+04    2.400  -2110.00           
OH+HO2<=>O2+H2O                          2.900E+13    0.000   -500.00           
OH+H2O2<=>HO2+H2O                        1.750E+12    0.000    320.00           
 DUPLICATE                                                                      
OH+H2O2<=>HO2+H2O                        5.800E+14    0.000   9560.00           
 DUPLICATE                                                                      
OH+C<=>H+CO                              5.000E+13    0.000      0.00           
OH+CH<=>H+HCO                            3.000E+13    0.000      0.00           
OH+CH2<=>H+CH2O                          2.000E+13    0.000      0.00           
OH+CH2<=>CH+H2O                          1.130E+07    2.000   3000.00           
OH+CH2(S)<=>H+CH2O                       3.000E+13    0.000      0.00           
OH+CH3(+M)<=>CH3OH(+M)                   6.300E+13    0.000      0.00           
     LOW  /  2.700E+38   -6.300   3100.00/                                      
     TROE/  0.2105   83.50  5398.00  8370.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/                      
OH+CH3<=>CH2+H2O                         5.600E+07    1.600   5420.00           
OH+CH3<=>CH2(S)+H2O                      2.501E+13    0.000      0.00           
OH+CH4<=>CH3+H2O                         1.000E+08    1.600   3120.00           
OH+CO<=>H+CO2                            4.760E+07    1.228     70.00           
OH+HCO<=>H2O+CO                          5.000E+13    0.000      0.00           
OH+CH2O<=>HCO+H2O                        3.430E+09    1.180   -447.00           
OH+CH2OH<=>H2O+CH2O                      5.000E+12    0.000      0.00           
OH+CH3O<=>H2O+CH2O                       5.000E+12    0.000      0.00           
OH+CH3OH<=>CH2OH+H2O                     1.440E+06    2.000   -840.00           
OH+CH3OH<=>CH3O+H2O                      6.300E+06    2.000   1500.00           
OH+C2H<=>H+HCCO                          2.000E+13    0.000      0.00           
OH+C2H2<=>H+CH2CO                        2.180E-04    4.500  -1000.00           
OH+C2H2<=>H+HCCOH                        5.040E+05    2.300  13500.00           
OH+C2H2<=>C2H+H2O                        3.370E+07    2.000  14000.00           
OH+C2H2<=>CH3+CO                         4.830E-04    4.000  -2000.00           
OH+C2H3<=>H2O+C2H2                       5.000E+12    0.000      0.00           
OH+C2H4<=>C2H3+H2O                       3.600E+06    2.000   2500.00           
OH+C2H6<=>C2H5+H2O                       3.540E+06    2.120    870.00           
OH+CH2CO<=>HCCO+H2O                      7.500E+12    0.000   2000.00           
2HO2<=>O2+H2O2                           1.300E+11    0.000  -1630.00           
 DUPLICATE                                                                      
2HO2<=>O2+H2O2                           4.200E+14    0.000  12000.00           
 DUPLICATE                                                                      
HO2+CH2<=>OH+CH2O                        2.000E+13    0.000      0.00           
HO2+CH3<=>O2+CH4                         1.000E+12    0.000      0.00           
HO2+CH3<=>OH+CH3O                        2.000E+13    0.000      0.00           
HO2+CO<=>OH+CO2                          1.500E+14    0.000  23600.00           
HO2+CH2O<=>HCO+H2O2                      1.000E+12    0.000   8000.00           
C+O2<=>O+CO                              5.800E+13    0.000    576.00           
C+CH2<=>H+C2H                            5.000E+13    0.000      0.00           
C+CH3<=>H+C2H2                           5.000E+13    0.000      0.00           
CH+O2<=>O+HCO                            3.300E+13    0.000      0.00           
CH+H2<=>H+CH2                            1.107E+08    1.790   1670.00           
CH+H2O<=>H+CH2O                          5.710E+12    0.000   -755.00           
CH+CH2<=>H+C2H2                          4.000E+13    0.000      0.00           
CH+CH3<=>H+C2H3                          3.000E+13    0.000      0.00           
CH+CH4<=>H+C2H4                          6.000E+13    0.000      0.00           
CH+CO(+M)<=>HCCO(+M)                     5.000E+13    0.000      0.00           
     LOW  /  2.690E+28   -3.740   1936.00/                                      
     TROE/  0.5757  237.00  1652.00  5069.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
CH+CO2<=>HCO+CO                          3.400E+12    0.000    690.00           
CH+CH2O<=>H+CH2CO                        9.460E+13    0.000   -515.00           
CH+HCCO<=>CO+C2H2                        5.000E+13    0.000      0.00           
CH2+O2<=>OH+HCO                          1.320E+13    0.000   1500.00           
CH2+H2<=>H+CH3                           5.000E+05    2.000   7230.00           
2CH2<=>H2+C2H2                           3.200E+13    0.000      0.00           
CH2+CH3<=>H+C2H4                         4.000E+13    0.000      0.00           
CH2+CH4<=>2CH3                           2.460E+06    2.000   8270.00           
CH2+CO(+M)<=>CH2CO(+M)                   8.100E+11    0.500   4510.00           
     LOW  /  2.690E+33   -5.110   7095.00/                                      
     TROE/  0.5907  275.00  1226.00  5185.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
CH2+HCCO<=>C2H3+CO                       3.000E+13    0.000      0.00           
CH2(S)+N2<=>CH2+N2                       1.500E+13    0.000    600.00           
CH2(S)+AR<=>CH2+AR                       9.000E+12    0.000    600.00           
CH2(S)+O2<=>H+OH+CO                      2.800E+13    0.000      0.00           
CH2(S)+O2<=>CO+H2O                       1.200E+13    0.000      0.00           
CH2(S)+H2<=>CH3+H                        7.000E+13    0.000      0.00           
CH2(S)+H2O(+M)<=>CH3OH(+M)               2.000E+13    0.000      0.00           
     LOW  /  2.700E+38   -6.300   3100.00/                                      
     TROE/  0.1507  134.00  2383.00  7265.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/                      
CH2(S)+H2O<=>CH2+H2O                     3.000E+13    0.000      0.00           
CH2(S)+CH3<=>H+C2H4                      1.200E+13    0.000   -570.00           
CH2(S)+CH4<=>2CH3                        1.600E+13    0.000   -570.00           
CH2(S)+CO<=>CH2+CO                       9.000E+12    0.000      0.00           
CH2(S)+CO2<=>CH2+CO2                     7.000E+12    0.000      0.00           
CH2(S)+CO2<=>CO+CH2O                     1.400E+13    0.000      0.00           
CH2(S)+C2H6<=>CH3+C2H5                   4.000E+13    0.000   -550.00           
CH3+O2<=>O+CH3O                          2.675E+13    0.000  28800.00           
CH3+O2<=>OH+CH2O                         3.600E+10    0.000   8940.00           
CH3+H2O2<=>HO2+CH4                       2.450E+04    2.470   5180.00           
2CH3(+M)<=>C2H6(+M)                      2.120E+16   -0.970    620.00           
     LOW  /  1.770E+50   -9.670   6220.00/                                      
     TROE/  0.5325  151.00  1038.00  4970.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
2CH3<=>H+C2H5                            4.990E+12    0.100  10600.00           
CH3+HCO<=>CH4+CO                         2.648E+13    0.000      0.00           
CH3+CH2O<=>HCO+CH4                       3.320E+03    2.810   5860.00           
CH3+CH3OH<=>CH2OH+CH4                    3.000E+07    1.500   9940.00           
CH3+CH3OH<=>CH3O+CH4                     1.000E+07    1.500   9940.00           
CH3+C2H4<=>C2H3+CH4                      2.270E+05    2.000   9200.00           
CH3+C2H6<=>C2H5+CH4                      6.140E+06    1.740  10450.00           
HCO+H2O<=>H+CO+H2O                       2.244E+18   -1.000  17000.00           
HCO+M<=>H+CO+M                           1.870E+17   -1.000  17000.00           
H2/2.00/ H2O/0.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/                      
HCO+O2<=>HO2+CO                          7.600E+12    0.000    400.00           
CH2OH+O2<=>HO2+CH2O                      1.800E+13    0.000    900.00           
CH3O+O2<=>HO2+CH2O                       4.280E-13    7.600  -3530.00           
C2H+O2<=>HCO+CO                          5.000E+13    0.000   1500.00           
C2H+H2<=>H+C2H2                          4.070E+05    2.400    200.00           
C2H3+O2<=>HCO+CH2O                       3.980E+12    0.000   -240.00           
C2H4(+M)<=>H2+C2H2(+M)                   8.000E+12    0.440  88770.00           
     LOW  /  7.000E+50   -9.310  99860.00/                                      
     TROE/  0.7345  180.00  1035.00  5417.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
C2H5+O2<=>HO2+C2H4                       8.400E+11    0.000   3875.00           
HCCO+O2<=>OH+2CO                         1.600E+12    0.000    854.00           
2HCCO<=>2CO+C2H2                         1.000E+13    0.000      0.00           
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
