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
    *kk = 35;
    *ii = 217;
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

    /* H  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = ' ';

    /* C  */
    kname[ 2*lenkname + 0 ] = 'C';
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
    for (i=0; i<lenkname*35; i++) {
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

    /* C3H7  */
    kname[ 31*lenkname + 0 ] = 'C';
    kname[ 31*lenkname + 1 ] = '3';
    kname[ 31*lenkname + 2 ] = 'H';
    kname[ 31*lenkname + 3 ] = '7';
    kname[ 31*lenkname + 4 ] = ' ';

    /* C3H8  */
    kname[ 32*lenkname + 0 ] = 'C';
    kname[ 32*lenkname + 1 ] = '3';
    kname[ 32*lenkname + 2 ] = 'H';
    kname[ 32*lenkname + 3 ] = '8';
    kname[ 32*lenkname + 4 ] = ' ';

    /* CH2CHO  */
    kname[ 33*lenkname + 0 ] = 'C';
    kname[ 33*lenkname + 1 ] = 'H';
    kname[ 33*lenkname + 2 ] = '2';
    kname[ 33*lenkname + 3 ] = 'C';
    kname[ 33*lenkname + 4 ] = 'H';
    kname[ 33*lenkname + 5 ] = 'O';
    kname[ 33*lenkname + 6 ] = ' ';

    /* CH3CHO  */
    kname[ 34*lenkname + 0 ] = 'C';
    kname[ 34*lenkname + 1 ] = 'H';
    kname[ 34*lenkname + 2 ] = '3';
    kname[ 34*lenkname + 3 ] = 'C';
    kname[ 34*lenkname + 4 ] = 'H';
    kname[ 34*lenkname + 5 ] = 'O';
    kname[ 34*lenkname + 6 ] = ' ';

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
    XW += x[31]*43.089240; /*C3H7 */
    XW += x[32]*44.097210; /*C3H8 */
    XW += x[33]*43.045610; /*CH2CHO */
    XW += x[34]*44.053580; /*CH3CHO */
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
    YOW += y[31]/43.089240; /*C3H7 */
    YOW += y[32]/44.097210; /*C3H8 */
    YOW += y[33]/43.045610; /*CH2CHO */
    YOW += y[34]/44.053580; /*CH3CHO */
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
    W += c[31]*43.089240; /*C3H7 */
    W += c[32]*44.097210; /*C3H8 */
    W += c[33]*43.045610; /*CH2CHO */
    W += c[34]*44.053580; /*CH3CHO */

    for (id = 0; id < 35; ++id) {
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
    XW += x[31]*43.089240; /*C3H7 */
    XW += x[32]*44.097210; /*C3H8 */
    XW += x[33]*43.045610; /*CH2CHO */
    XW += x[34]*44.053580; /*CH3CHO */
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
    YOW += y[31]/43.089240; /*C3H7 */
    YOW += y[32]/44.097210; /*C3H8 */
    YOW += y[33]/43.045610; /*CH2CHO */
    YOW += y[34]/44.053580; /*CH3CHO */
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
    W += c[31]*43.089240; /*C3H7 */
    W += c[32]*44.097210; /*C3H8 */
    W += c[33]*43.045610; /*CH2CHO */
    W += c[34]*44.053580; /*CH3CHO */

    for (id = 0; id < 35; ++id) {
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
    YOW += y[31]/43.089240; /*C3H7 */
    YOW += y[32]/44.097210; /*C3H8 */
    YOW += y[33]/43.045610; /*CH2CHO */
    YOW += y[34]/44.053580; /*CH3CHO */
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
    XW += x[31]*43.089240; /*C3H7 */
    XW += x[32]*44.097210; /*C3H8 */
    XW += x[33]*43.045610; /*CH2CHO */
    XW += x[34]*44.053580; /*CH3CHO */
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
    W += c[31]*43.089240; /*C3H7 */
    W += c[32]*44.097210; /*C3H8 */
    W += c[33]*43.045610; /*CH2CHO */
    W += c[34]*44.053580; /*CH3CHO */

    for (id = 0; id < 35; ++id) {
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
    YOW += y[31]/43.089240; /*C3H7 */
    YOW += y[32]/44.097210; /*C3H8 */
    YOW += y[33]/43.045610; /*CH2CHO */
    YOW += y[34]/44.053580; /*CH3CHO */
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
    x[31] = y[31]/(43.089240*YOW); 
    x[32] = y[32]/(44.097210*YOW); 
    x[33] = y[33]/(43.045610*YOW); 
    x[34] = y[34]/(44.053580*YOW); 

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
    YOW += y[31]/43.089240; /*C3H7 */
    YOW += y[32]/44.097210; /*C3H8 */
    YOW += y[33]/43.045610; /*CH2CHO */
    YOW += y[34]/44.053580; /*CH3CHO */
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
    c[31] = PWORT * y[31]/43.089240; 
    c[32] = PWORT * y[32]/44.097210; 
    c[33] = PWORT * y[33]/43.045610; 
    c[34] = PWORT * y[34]/44.053580; 

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
    c[31] = (*rho) * y[31]/43.089240; 
    c[32] = (*rho) * y[32]/44.097210; 
    c[33] = (*rho) * y[33]/43.045610; 
    c[34] = (*rho) * y[34]/44.053580; 

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
    XW += x[31]*43.089240; /*C3H7 */
    XW += x[32]*44.097210; /*C3H8 */
    XW += x[33]*43.045610; /*CH2CHO */
    XW += x[34]*44.053580; /*CH3CHO */
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
    y[31] = x[31]*43.089240/XW; 
    y[32] = x[32]*44.097210/XW; 
    y[33] = x[33]*43.045610/XW; 
    y[34] = x[34]*44.053580/XW; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 35; ++id) {
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
    XW += x[31]*43.089240; /*C3H7 */
    XW += x[32]*44.097210; /*C3H8 */
    XW += x[33]*43.045610; /*CH2CHO */
    XW += x[34]*44.053580; /*CH3CHO */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 35; ++id) {
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
    for (id = 0; id < 35; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    for (id = 0; id < 35; ++id) {
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
    CW += c[31]*43.089240; /*C3H7 */
    CW += c[32]*44.097210; /*C3H8 */
    CW += c[33]*43.045610; /*CH2CHO */
    CW += c[34]*44.053580; /*CH3CHO */
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
    y[31] = c[31]*43.089240/CW; 
    y[32] = c[32]*44.097210/CW; 
    y[33] = c[33]*43.045610/CW; 
    y[34] = c[34]*44.053580/CW; 

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
    for (id = 0; id < 35; ++id) {
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
    for (id = 0; id < 35; ++id) {
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
    for (id = 0; id < 35; ++id) {
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
    for (id = 0; id < 35; ++id) {
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
    for (id = 0; id < 35; ++id) {
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
    for (id = 0; id < 35; ++id) {
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
    for (id = 0; id < 35; ++id) {
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
    cvms[31] *= 1.929602378691293e+06; /*C3H7 */
    cvms[32] *= 1.885495703696447e+06; /*C3H8 */
    cvms[33] *= 1.931558177477332e+06; /*CH2CHO */
    cvms[34] *= 1.887363070152301e+06; /*CH3CHO */
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
    cpms[31] *= 1.929602378691293e+06; /*C3H7 */
    cpms[32] *= 1.885495703696447e+06; /*C3H8 */
    cpms[33] *= 1.931558177477332e+06; /*CH2CHO */
    cpms[34] *= 1.887363070152301e+06; /*CH3CHO */
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
    ums[31] *= RT/43.089240; /*C3H7 */
    ums[32] *= RT/44.097210; /*C3H8 */
    ums[33] *= RT/43.045610; /*CH2CHO */
    ums[34] *= RT/44.053580; /*CH3CHO */
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
    hms[31] *= RT/43.089240; /*C3H7 */
    hms[32] *= RT/44.097210; /*C3H8 */
    hms[33] *= RT/43.045610; /*CH2CHO */
    hms[34] *= RT/44.053580; /*CH3CHO */
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
    gms[31] *= RT/43.089240; /*C3H7 */
    gms[32] *= RT/44.097210; /*C3H8 */
    gms[33] *= RT/43.045610; /*CH2CHO */
    gms[34] *= RT/44.053580; /*CH3CHO */
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
    ams[31] *= RT/43.089240; /*C3H7 */
    ams[32] *= RT/44.097210; /*C3H8 */
    ams[33] *= RT/43.045610; /*CH2CHO */
    ams[34] *= RT/44.053580; /*CH3CHO */
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
    sms[31] *= 1.929602378691293e+06; /*C3H7 */
    sms[32] *= 1.885495703696447e+06; /*C3H8 */
    sms[33] *= 1.931558177477332e+06; /*CH2CHO */
    sms[34] *= 1.887363070152301e+06; /*CH3CHO */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double *T, double *x, int * iwrk, double * rwrk, double * cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[35]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 35; ++id) {
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
    double cpor[35]; /* temporary storage */
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
    result += cpor[31]*y[31]/43.089240; /*C3H7 */
    result += cpor[32]*y[32]/44.097210; /*C3H8 */
    result += cpor[33]*y[33]/43.045610; /*CH2CHO */
    result += cpor[34]*y[34]/44.053580; /*CH3CHO */

    *cpbs = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double *T, double *x, int * iwrk, double * rwrk, double * cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[35]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 35; ++id) {
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
    double cvor[35]; /* temporary storage */
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
    result += cvor[31]*y[31]/43.089240; /*C3H7 */
    result += cvor[32]*y[32]/44.097210; /*C3H8 */
    result += cvor[33]*y[33]/43.045610; /*CH2CHO */
    result += cvor[34]*y[34]/44.053580; /*CH3CHO */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double *T, double *x, int * iwrk, double * rwrk, double * hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[35]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 35; ++id) {
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
    double hml[35]; /* temporary storage */
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
    result += y[31]*hml[31]/43.089240; /*C3H7 */
    result += y[32]*hml[32]/44.097210; /*C3H8 */
    result += y[33]*hml[33]/43.045610; /*CH2CHO */
    result += y[34]*hml[34]/44.053580; /*CH3CHO */

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double *T, double *x, int * iwrk, double * rwrk, double * ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[35]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 35; ++id) {
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
    double ums[35]; /* temporary energy array */
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
    result += y[31]*ums[31]/43.089240; /*C3H7 */
    result += y[32]*ums[32]/44.097210; /*C3H8 */
    result += y[33]*ums[33]/43.045610; /*CH2CHO */
    result += y[34]*ums[34]/44.053580; /*CH3CHO */

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
    double sor[35]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 35; ++id) {
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
    double sor[35]; /* temporary storage */
    double x[35]; /* need a ytx conversion */
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
    YOW += y[31]/43.089240; /*C3H7 */
    YOW += y[32]/44.097210; /*C3H8 */
    YOW += y[33]/43.045610; /*CH2CHO */
    YOW += y[34]/44.053580; /*CH3CHO */
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
    x[31] = y[31]/(43.089240*YOW); 
    x[32] = y[32]/(44.097210*YOW); 
    x[33] = y[33]/(43.045610*YOW); 
    x[34] = y[34]/(44.053580*YOW); 
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
    double gort[35]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 35; ++id) {
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
    double gort[35]; /* temporary storage */
    double x[35]; /* need a ytx conversion */
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
    YOW += y[31]/43.089240; /*C3H7 */
    YOW += y[32]/44.097210; /*C3H8 */
    YOW += y[33]/43.045610; /*CH2CHO */
    YOW += y[34]/44.053580; /*CH3CHO */
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
    x[31] = y[31]/(43.089240*YOW); 
    x[32] = y[32]/(44.097210*YOW); 
    x[33] = y[33]/(43.045610*YOW); 
    x[34] = y[34]/(44.053580*YOW); 
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
    double aort[35]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 35; ++id) {
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
    double aort[35]; /* temporary storage */
    double x[35]; /* need a ytx conversion */
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
    YOW += y[31]/43.089240; /*C3H7 */
    YOW += y[32]/44.097210; /*C3H8 */
    YOW += y[33]/43.045610; /*CH2CHO */
    YOW += y[34]/44.053580; /*CH3CHO */
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
    x[31] = y[31]/(43.089240*YOW); 
    x[32] = y[32]/(44.097210*YOW); 
    x[33] = y[33]/(43.045610*YOW); 
    x[34] = y[34]/(44.053580*YOW); 
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
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double * T, double * C, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 35; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 35; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[35]; /*temporary storage */
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
    YOW += y[31]/43.089240; /*C3H7 */
    YOW += y[32]/44.097210; /*C3H8 */
    YOW += y[33]/43.045610; /*CH2CHO */
    YOW += y[34]/44.053580; /*CH3CHO */
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
    c[31] = PWORT * y[31]/43.089240; 
    c[32] = PWORT * y[32]/44.097210; 
    c[33] = PWORT * y[33]/43.045610; 
    c[34] = PWORT * y[34]/44.053580; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 35; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[35]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 35; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 35; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[35]; /*temporary storage */
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
    c[31] = 1e6 * (*rho) * y[31]/43.089240; 
    c[32] = 1e6 * (*rho) * y[32]/44.097210; 
    c[33] = 1e6 * (*rho) * y[33]/43.045610; 
    c[34] = 1e6 * (*rho) * y[34]/44.053580; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 35; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[35]; /*temporary storage */
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
    XW += x[31]*43.089240; /*C3H7 */
    XW += x[32]*44.097210; /*C3H8 */
    XW += x[33]*43.045610; /*CH2CHO */
    XW += x[34]*44.053580; /*CH3CHO */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 35; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 35; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double * T, double * C, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 35; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 35; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 217; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double * P, double * T, double * x, int * iwrk, double * rwrk, double * q_f, double * q_r)
{
    int id; /*loop counter */
    double c[35]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 35; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 217; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[35]; /*temporary storage */
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
    YOW += y[31]/43.089240; /*C3H7 */
    YOW += y[32]/44.097210; /*C3H8 */
    YOW += y[33]/43.045610; /*CH2CHO */
    YOW += y[34]/44.053580; /*CH3CHO */
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
    c[31] = PWORT * y[31]/43.089240; 
    c[32] = PWORT * y[32]/44.097210; 
    c[33] = PWORT * y[33]/43.045610; 
    c[34] = PWORT * y[34]/44.053580; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 217; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[35]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 35; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 217; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[35]; /*temporary storage */
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
    c[31] = 1e6 * (*rho) * y[31]/43.089240; 
    c[32] = 1e6 * (*rho) * y[32]/44.097210; 
    c[33] = 1e6 * (*rho) * y[33]/43.045610; 
    c[34] = 1e6 * (*rho) * y[34]/44.053580; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 217; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[35]; /*temporary storage */
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
    XW += x[31]*43.089240; /*C3H7 */
    XW += x[32]*44.097210; /*C3H8 */
    XW += x[33]*43.045610; /*CH2CHO */
    XW += x[34]*44.053580; /*CH3CHO */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 35; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 217; ++id) {
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
    for (id = 0; id < 35 * kd; ++ id) {
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
    nuki[ 30 * kd + 35 ] += -1 ;
    nuki[ 6 * kd + 35 ] += +1 ;
    nuki[ 30 * kd + 35 ] += +1 ;

    /*reaction 37: H + O2 <=> O + OH */
    nuki[ 1 * kd + 36 ] += -1 ;
    nuki[ 3 * kd + 36 ] += -1 ;
    nuki[ 2 * kd + 36 ] += +1 ;
    nuki[ 4 * kd + 36 ] += +1 ;

    /*reaction 38: 2 H + M <=> H2 + M */
    nuki[ 1 * kd + 37 ] += -2 ;
    nuki[ 0 * kd + 37 ] += +1 ;

    /*reaction 39: 2 H + H2 <=> 2 H2 */
    nuki[ 1 * kd + 38 ] += -2 ;
    nuki[ 0 * kd + 38 ] += -1 ;
    nuki[ 0 * kd + 38 ] += +2 ;

    /*reaction 40: 2 H + H2O <=> H2 + H2O */
    nuki[ 1 * kd + 39 ] += -2 ;
    nuki[ 5 * kd + 39 ] += -1 ;
    nuki[ 0 * kd + 39 ] += +1 ;
    nuki[ 5 * kd + 39 ] += +1 ;

    /*reaction 41: 2 H + CO2 <=> H2 + CO2 */
    nuki[ 1 * kd + 40 ] += -2 ;
    nuki[ 15 * kd + 40 ] += -1 ;
    nuki[ 0 * kd + 40 ] += +1 ;
    nuki[ 15 * kd + 40 ] += +1 ;

    /*reaction 42: H + OH + M <=> H2O + M */
    nuki[ 1 * kd + 41 ] += -1 ;
    nuki[ 4 * kd + 41 ] += -1 ;
    nuki[ 5 * kd + 41 ] += +1 ;

    /*reaction 43: H + HO2 <=> O + H2O */
    nuki[ 1 * kd + 42 ] += -1 ;
    nuki[ 6 * kd + 42 ] += -1 ;
    nuki[ 2 * kd + 42 ] += +1 ;
    nuki[ 5 * kd + 42 ] += +1 ;

    /*reaction 44: H + HO2 <=> O2 + H2 */
    nuki[ 1 * kd + 43 ] += -1 ;
    nuki[ 6 * kd + 43 ] += -1 ;
    nuki[ 3 * kd + 43 ] += +1 ;
    nuki[ 0 * kd + 43 ] += +1 ;

    /*reaction 45: H + HO2 <=> 2 OH */
    nuki[ 1 * kd + 44 ] += -1 ;
    nuki[ 6 * kd + 44 ] += -1 ;
    nuki[ 4 * kd + 44 ] += +2 ;

    /*reaction 46: H + H2O2 <=> HO2 + H2 */
    nuki[ 1 * kd + 45 ] += -1 ;
    nuki[ 7 * kd + 45 ] += -1 ;
    nuki[ 6 * kd + 45 ] += +1 ;
    nuki[ 0 * kd + 45 ] += +1 ;

    /*reaction 47: H + H2O2 <=> OH + H2O */
    nuki[ 1 * kd + 46 ] += -1 ;
    nuki[ 7 * kd + 46 ] += -1 ;
    nuki[ 4 * kd + 46 ] += +1 ;
    nuki[ 5 * kd + 46 ] += +1 ;

    /*reaction 48: H + CH <=> C + H2 */
    nuki[ 1 * kd + 47 ] += -1 ;
    nuki[ 9 * kd + 47 ] += -1 ;
    nuki[ 8 * kd + 47 ] += +1 ;
    nuki[ 0 * kd + 47 ] += +1 ;

    /*reaction 49: H + CH2 (+M) <=> CH3 (+M) */
    nuki[ 1 * kd + 48 ] += -1 ;
    nuki[ 10 * kd + 48 ] += -1 ;
    nuki[ 12 * kd + 48 ] += +1 ;

    /*reaction 50: H + CH2(S) <=> CH + H2 */
    nuki[ 1 * kd + 49 ] += -1 ;
    nuki[ 11 * kd + 49 ] += -1 ;
    nuki[ 9 * kd + 49 ] += +1 ;
    nuki[ 0 * kd + 49 ] += +1 ;

    /*reaction 51: H + CH3 (+M) <=> CH4 (+M) */
    nuki[ 1 * kd + 50 ] += -1 ;
    nuki[ 12 * kd + 50 ] += -1 ;
    nuki[ 13 * kd + 50 ] += +1 ;

    /*reaction 52: H + CH4 <=> CH3 + H2 */
    nuki[ 1 * kd + 51 ] += -1 ;
    nuki[ 13 * kd + 51 ] += -1 ;
    nuki[ 12 * kd + 51 ] += +1 ;
    nuki[ 0 * kd + 51 ] += +1 ;

    /*reaction 53: H + HCO (+M) <=> CH2O (+M) */
    nuki[ 1 * kd + 52 ] += -1 ;
    nuki[ 16 * kd + 52 ] += -1 ;
    nuki[ 17 * kd + 52 ] += +1 ;

    /*reaction 54: H + HCO <=> H2 + CO */
    nuki[ 1 * kd + 53 ] += -1 ;
    nuki[ 16 * kd + 53 ] += -1 ;
    nuki[ 0 * kd + 53 ] += +1 ;
    nuki[ 14 * kd + 53 ] += +1 ;

    /*reaction 55: H + CH2O (+M) <=> CH2OH (+M) */
    nuki[ 1 * kd + 54 ] += -1 ;
    nuki[ 17 * kd + 54 ] += -1 ;
    nuki[ 18 * kd + 54 ] += +1 ;

    /*reaction 56: H + CH2O (+M) <=> CH3O (+M) */
    nuki[ 1 * kd + 55 ] += -1 ;
    nuki[ 17 * kd + 55 ] += -1 ;
    nuki[ 19 * kd + 55 ] += +1 ;

    /*reaction 57: H + CH2O <=> HCO + H2 */
    nuki[ 1 * kd + 56 ] += -1 ;
    nuki[ 17 * kd + 56 ] += -1 ;
    nuki[ 16 * kd + 56 ] += +1 ;
    nuki[ 0 * kd + 56 ] += +1 ;

    /*reaction 58: H + CH2OH (+M) <=> CH3OH (+M) */
    nuki[ 1 * kd + 57 ] += -1 ;
    nuki[ 18 * kd + 57 ] += -1 ;
    nuki[ 20 * kd + 57 ] += +1 ;

    /*reaction 59: H + CH2OH <=> H2 + CH2O */
    nuki[ 1 * kd + 58 ] += -1 ;
    nuki[ 18 * kd + 58 ] += -1 ;
    nuki[ 0 * kd + 58 ] += +1 ;
    nuki[ 17 * kd + 58 ] += +1 ;

    /*reaction 60: H + CH2OH <=> OH + CH3 */
    nuki[ 1 * kd + 59 ] += -1 ;
    nuki[ 18 * kd + 59 ] += -1 ;
    nuki[ 4 * kd + 59 ] += +1 ;
    nuki[ 12 * kd + 59 ] += +1 ;

    /*reaction 61: H + CH2OH <=> CH2(S) + H2O */
    nuki[ 1 * kd + 60 ] += -1 ;
    nuki[ 18 * kd + 60 ] += -1 ;
    nuki[ 11 * kd + 60 ] += +1 ;
    nuki[ 5 * kd + 60 ] += +1 ;

    /*reaction 62: H + CH3O (+M) <=> CH3OH (+M) */
    nuki[ 1 * kd + 61 ] += -1 ;
    nuki[ 19 * kd + 61 ] += -1 ;
    nuki[ 20 * kd + 61 ] += +1 ;

    /*reaction 63: H + CH3O <=> H + CH2OH */
    nuki[ 1 * kd + 62 ] += -1 ;
    nuki[ 19 * kd + 62 ] += -1 ;
    nuki[ 1 * kd + 62 ] += +1 ;
    nuki[ 18 * kd + 62 ] += +1 ;

    /*reaction 64: H + CH3O <=> H2 + CH2O */
    nuki[ 1 * kd + 63 ] += -1 ;
    nuki[ 19 * kd + 63 ] += -1 ;
    nuki[ 0 * kd + 63 ] += +1 ;
    nuki[ 17 * kd + 63 ] += +1 ;

    /*reaction 65: H + CH3O <=> OH + CH3 */
    nuki[ 1 * kd + 64 ] += -1 ;
    nuki[ 19 * kd + 64 ] += -1 ;
    nuki[ 4 * kd + 64 ] += +1 ;
    nuki[ 12 * kd + 64 ] += +1 ;

    /*reaction 66: H + CH3O <=> CH2(S) + H2O */
    nuki[ 1 * kd + 65 ] += -1 ;
    nuki[ 19 * kd + 65 ] += -1 ;
    nuki[ 11 * kd + 65 ] += +1 ;
    nuki[ 5 * kd + 65 ] += +1 ;

    /*reaction 67: H + CH3OH <=> CH2OH + H2 */
    nuki[ 1 * kd + 66 ] += -1 ;
    nuki[ 20 * kd + 66 ] += -1 ;
    nuki[ 18 * kd + 66 ] += +1 ;
    nuki[ 0 * kd + 66 ] += +1 ;

    /*reaction 68: H + CH3OH <=> CH3O + H2 */
    nuki[ 1 * kd + 67 ] += -1 ;
    nuki[ 20 * kd + 67 ] += -1 ;
    nuki[ 19 * kd + 67 ] += +1 ;
    nuki[ 0 * kd + 67 ] += +1 ;

    /*reaction 69: H + C2H (+M) <=> C2H2 (+M) */
    nuki[ 1 * kd + 68 ] += -1 ;
    nuki[ 21 * kd + 68 ] += -1 ;
    nuki[ 22 * kd + 68 ] += +1 ;

    /*reaction 70: H + C2H2 (+M) <=> C2H3 (+M) */
    nuki[ 1 * kd + 69 ] += -1 ;
    nuki[ 22 * kd + 69 ] += -1 ;
    nuki[ 23 * kd + 69 ] += +1 ;

    /*reaction 71: H + C2H3 (+M) <=> C2H4 (+M) */
    nuki[ 1 * kd + 70 ] += -1 ;
    nuki[ 23 * kd + 70 ] += -1 ;
    nuki[ 24 * kd + 70 ] += +1 ;

    /*reaction 72: H + C2H3 <=> H2 + C2H2 */
    nuki[ 1 * kd + 71 ] += -1 ;
    nuki[ 23 * kd + 71 ] += -1 ;
    nuki[ 0 * kd + 71 ] += +1 ;
    nuki[ 22 * kd + 71 ] += +1 ;

    /*reaction 73: H + C2H4 (+M) <=> C2H5 (+M) */
    nuki[ 1 * kd + 72 ] += -1 ;
    nuki[ 24 * kd + 72 ] += -1 ;
    nuki[ 25 * kd + 72 ] += +1 ;

    /*reaction 74: H + C2H4 <=> C2H3 + H2 */
    nuki[ 1 * kd + 73 ] += -1 ;
    nuki[ 24 * kd + 73 ] += -1 ;
    nuki[ 23 * kd + 73 ] += +1 ;
    nuki[ 0 * kd + 73 ] += +1 ;

    /*reaction 75: H + C2H5 (+M) <=> C2H6 (+M) */
    nuki[ 1 * kd + 74 ] += -1 ;
    nuki[ 25 * kd + 74 ] += -1 ;
    nuki[ 26 * kd + 74 ] += +1 ;

    /*reaction 76: H + C2H5 <=> H2 + C2H4 */
    nuki[ 1 * kd + 75 ] += -1 ;
    nuki[ 25 * kd + 75 ] += -1 ;
    nuki[ 0 * kd + 75 ] += +1 ;
    nuki[ 24 * kd + 75 ] += +1 ;

    /*reaction 77: H + C2H6 <=> C2H5 + H2 */
    nuki[ 1 * kd + 76 ] += -1 ;
    nuki[ 26 * kd + 76 ] += -1 ;
    nuki[ 25 * kd + 76 ] += +1 ;
    nuki[ 0 * kd + 76 ] += +1 ;

    /*reaction 78: H + HCCO <=> CH2(S) + CO */
    nuki[ 1 * kd + 77 ] += -1 ;
    nuki[ 27 * kd + 77 ] += -1 ;
    nuki[ 11 * kd + 77 ] += +1 ;
    nuki[ 14 * kd + 77 ] += +1 ;

    /*reaction 79: H + CH2CO <=> HCCO + H2 */
    nuki[ 1 * kd + 78 ] += -1 ;
    nuki[ 28 * kd + 78 ] += -1 ;
    nuki[ 27 * kd + 78 ] += +1 ;
    nuki[ 0 * kd + 78 ] += +1 ;

    /*reaction 80: H + CH2CO <=> CH3 + CO */
    nuki[ 1 * kd + 79 ] += -1 ;
    nuki[ 28 * kd + 79 ] += -1 ;
    nuki[ 12 * kd + 79 ] += +1 ;
    nuki[ 14 * kd + 79 ] += +1 ;

    /*reaction 81: H + HCCOH <=> H + CH2CO */
    nuki[ 1 * kd + 80 ] += -1 ;
    nuki[ 29 * kd + 80 ] += -1 ;
    nuki[ 1 * kd + 80 ] += +1 ;
    nuki[ 28 * kd + 80 ] += +1 ;

    /*reaction 82: H2 + CO (+M) <=> CH2O (+M) */
    nuki[ 0 * kd + 81 ] += -1 ;
    nuki[ 14 * kd + 81 ] += -1 ;
    nuki[ 17 * kd + 81 ] += +1 ;

    /*reaction 83: OH + H2 <=> H + H2O */
    nuki[ 4 * kd + 82 ] += -1 ;
    nuki[ 0 * kd + 82 ] += -1 ;
    nuki[ 1 * kd + 82 ] += +1 ;
    nuki[ 5 * kd + 82 ] += +1 ;

    /*reaction 84: 2 OH (+M) <=> H2O2 (+M) */
    nuki[ 4 * kd + 83 ] += -2 ;
    nuki[ 7 * kd + 83 ] += +1 ;

    /*reaction 85: 2 OH <=> O + H2O */
    nuki[ 4 * kd + 84 ] += -2 ;
    nuki[ 2 * kd + 84 ] += +1 ;
    nuki[ 5 * kd + 84 ] += +1 ;

    /*reaction 86: OH + HO2 <=> O2 + H2O */
    nuki[ 4 * kd + 85 ] += -1 ;
    nuki[ 6 * kd + 85 ] += -1 ;
    nuki[ 3 * kd + 85 ] += +1 ;
    nuki[ 5 * kd + 85 ] += +1 ;

    /*reaction 87: OH + H2O2 <=> HO2 + H2O */
    nuki[ 4 * kd + 86 ] += -1 ;
    nuki[ 7 * kd + 86 ] += -1 ;
    nuki[ 6 * kd + 86 ] += +1 ;
    nuki[ 5 * kd + 86 ] += +1 ;

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    nuki[ 4 * kd + 87 ] += -1 ;
    nuki[ 7 * kd + 87 ] += -1 ;
    nuki[ 6 * kd + 87 ] += +1 ;
    nuki[ 5 * kd + 87 ] += +1 ;

    /*reaction 89: OH + C <=> H + CO */
    nuki[ 4 * kd + 88 ] += -1 ;
    nuki[ 8 * kd + 88 ] += -1 ;
    nuki[ 1 * kd + 88 ] += +1 ;
    nuki[ 14 * kd + 88 ] += +1 ;

    /*reaction 90: OH + CH <=> H + HCO */
    nuki[ 4 * kd + 89 ] += -1 ;
    nuki[ 9 * kd + 89 ] += -1 ;
    nuki[ 1 * kd + 89 ] += +1 ;
    nuki[ 16 * kd + 89 ] += +1 ;

    /*reaction 91: OH + CH2 <=> H + CH2O */
    nuki[ 4 * kd + 90 ] += -1 ;
    nuki[ 10 * kd + 90 ] += -1 ;
    nuki[ 1 * kd + 90 ] += +1 ;
    nuki[ 17 * kd + 90 ] += +1 ;

    /*reaction 92: OH + CH2 <=> CH + H2O */
    nuki[ 4 * kd + 91 ] += -1 ;
    nuki[ 10 * kd + 91 ] += -1 ;
    nuki[ 9 * kd + 91 ] += +1 ;
    nuki[ 5 * kd + 91 ] += +1 ;

    /*reaction 93: OH + CH2(S) <=> H + CH2O */
    nuki[ 4 * kd + 92 ] += -1 ;
    nuki[ 11 * kd + 92 ] += -1 ;
    nuki[ 1 * kd + 92 ] += +1 ;
    nuki[ 17 * kd + 92 ] += +1 ;

    /*reaction 94: OH + CH3 (+M) <=> CH3OH (+M) */
    nuki[ 4 * kd + 93 ] += -1 ;
    nuki[ 12 * kd + 93 ] += -1 ;
    nuki[ 20 * kd + 93 ] += +1 ;

    /*reaction 95: OH + CH3 <=> CH2 + H2O */
    nuki[ 4 * kd + 94 ] += -1 ;
    nuki[ 12 * kd + 94 ] += -1 ;
    nuki[ 10 * kd + 94 ] += +1 ;
    nuki[ 5 * kd + 94 ] += +1 ;

    /*reaction 96: OH + CH3 <=> CH2(S) + H2O */
    nuki[ 4 * kd + 95 ] += -1 ;
    nuki[ 12 * kd + 95 ] += -1 ;
    nuki[ 11 * kd + 95 ] += +1 ;
    nuki[ 5 * kd + 95 ] += +1 ;

    /*reaction 97: OH + CH4 <=> CH3 + H2O */
    nuki[ 4 * kd + 96 ] += -1 ;
    nuki[ 13 * kd + 96 ] += -1 ;
    nuki[ 12 * kd + 96 ] += +1 ;
    nuki[ 5 * kd + 96 ] += +1 ;

    /*reaction 98: OH + CO <=> H + CO2 */
    nuki[ 4 * kd + 97 ] += -1 ;
    nuki[ 14 * kd + 97 ] += -1 ;
    nuki[ 1 * kd + 97 ] += +1 ;
    nuki[ 15 * kd + 97 ] += +1 ;

    /*reaction 99: OH + HCO <=> H2O + CO */
    nuki[ 4 * kd + 98 ] += -1 ;
    nuki[ 16 * kd + 98 ] += -1 ;
    nuki[ 5 * kd + 98 ] += +1 ;
    nuki[ 14 * kd + 98 ] += +1 ;

    /*reaction 100: OH + CH2O <=> HCO + H2O */
    nuki[ 4 * kd + 99 ] += -1 ;
    nuki[ 17 * kd + 99 ] += -1 ;
    nuki[ 16 * kd + 99 ] += +1 ;
    nuki[ 5 * kd + 99 ] += +1 ;

    /*reaction 101: OH + CH2OH <=> H2O + CH2O */
    nuki[ 4 * kd + 100 ] += -1 ;
    nuki[ 18 * kd + 100 ] += -1 ;
    nuki[ 5 * kd + 100 ] += +1 ;
    nuki[ 17 * kd + 100 ] += +1 ;

    /*reaction 102: OH + CH3O <=> H2O + CH2O */
    nuki[ 4 * kd + 101 ] += -1 ;
    nuki[ 19 * kd + 101 ] += -1 ;
    nuki[ 5 * kd + 101 ] += +1 ;
    nuki[ 17 * kd + 101 ] += +1 ;

    /*reaction 103: OH + CH3OH <=> CH2OH + H2O */
    nuki[ 4 * kd + 102 ] += -1 ;
    nuki[ 20 * kd + 102 ] += -1 ;
    nuki[ 18 * kd + 102 ] += +1 ;
    nuki[ 5 * kd + 102 ] += +1 ;

    /*reaction 104: OH + CH3OH <=> CH3O + H2O */
    nuki[ 4 * kd + 103 ] += -1 ;
    nuki[ 20 * kd + 103 ] += -1 ;
    nuki[ 19 * kd + 103 ] += +1 ;
    nuki[ 5 * kd + 103 ] += +1 ;

    /*reaction 105: OH + C2H <=> H + HCCO */
    nuki[ 4 * kd + 104 ] += -1 ;
    nuki[ 21 * kd + 104 ] += -1 ;
    nuki[ 1 * kd + 104 ] += +1 ;
    nuki[ 27 * kd + 104 ] += +1 ;

    /*reaction 106: OH + C2H2 <=> H + CH2CO */
    nuki[ 4 * kd + 105 ] += -1 ;
    nuki[ 22 * kd + 105 ] += -1 ;
    nuki[ 1 * kd + 105 ] += +1 ;
    nuki[ 28 * kd + 105 ] += +1 ;

    /*reaction 107: OH + C2H2 <=> H + HCCOH */
    nuki[ 4 * kd + 106 ] += -1 ;
    nuki[ 22 * kd + 106 ] += -1 ;
    nuki[ 1 * kd + 106 ] += +1 ;
    nuki[ 29 * kd + 106 ] += +1 ;

    /*reaction 108: OH + C2H2 <=> C2H + H2O */
    nuki[ 4 * kd + 107 ] += -1 ;
    nuki[ 22 * kd + 107 ] += -1 ;
    nuki[ 21 * kd + 107 ] += +1 ;
    nuki[ 5 * kd + 107 ] += +1 ;

    /*reaction 109: OH + C2H2 <=> CH3 + CO */
    nuki[ 4 * kd + 108 ] += -1 ;
    nuki[ 22 * kd + 108 ] += -1 ;
    nuki[ 12 * kd + 108 ] += +1 ;
    nuki[ 14 * kd + 108 ] += +1 ;

    /*reaction 110: OH + C2H3 <=> H2O + C2H2 */
    nuki[ 4 * kd + 109 ] += -1 ;
    nuki[ 23 * kd + 109 ] += -1 ;
    nuki[ 5 * kd + 109 ] += +1 ;
    nuki[ 22 * kd + 109 ] += +1 ;

    /*reaction 111: OH + C2H4 <=> C2H3 + H2O */
    nuki[ 4 * kd + 110 ] += -1 ;
    nuki[ 24 * kd + 110 ] += -1 ;
    nuki[ 23 * kd + 110 ] += +1 ;
    nuki[ 5 * kd + 110 ] += +1 ;

    /*reaction 112: OH + C2H6 <=> C2H5 + H2O */
    nuki[ 4 * kd + 111 ] += -1 ;
    nuki[ 26 * kd + 111 ] += -1 ;
    nuki[ 25 * kd + 111 ] += +1 ;
    nuki[ 5 * kd + 111 ] += +1 ;

    /*reaction 113: OH + CH2CO <=> HCCO + H2O */
    nuki[ 4 * kd + 112 ] += -1 ;
    nuki[ 28 * kd + 112 ] += -1 ;
    nuki[ 27 * kd + 112 ] += +1 ;
    nuki[ 5 * kd + 112 ] += +1 ;

    /*reaction 114: 2 HO2 <=> O2 + H2O2 */
    nuki[ 6 * kd + 113 ] += -2 ;
    nuki[ 3 * kd + 113 ] += +1 ;
    nuki[ 7 * kd + 113 ] += +1 ;

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    nuki[ 6 * kd + 114 ] += -2 ;
    nuki[ 3 * kd + 114 ] += +1 ;
    nuki[ 7 * kd + 114 ] += +1 ;

    /*reaction 116: HO2 + CH2 <=> OH + CH2O */
    nuki[ 6 * kd + 115 ] += -1 ;
    nuki[ 10 * kd + 115 ] += -1 ;
    nuki[ 4 * kd + 115 ] += +1 ;
    nuki[ 17 * kd + 115 ] += +1 ;

    /*reaction 117: HO2 + CH3 <=> O2 + CH4 */
    nuki[ 6 * kd + 116 ] += -1 ;
    nuki[ 12 * kd + 116 ] += -1 ;
    nuki[ 3 * kd + 116 ] += +1 ;
    nuki[ 13 * kd + 116 ] += +1 ;

    /*reaction 118: HO2 + CH3 <=> OH + CH3O */
    nuki[ 6 * kd + 117 ] += -1 ;
    nuki[ 12 * kd + 117 ] += -1 ;
    nuki[ 4 * kd + 117 ] += +1 ;
    nuki[ 19 * kd + 117 ] += +1 ;

    /*reaction 119: HO2 + CO <=> OH + CO2 */
    nuki[ 6 * kd + 118 ] += -1 ;
    nuki[ 14 * kd + 118 ] += -1 ;
    nuki[ 4 * kd + 118 ] += +1 ;
    nuki[ 15 * kd + 118 ] += +1 ;

    /*reaction 120: HO2 + CH2O <=> HCO + H2O2 */
    nuki[ 6 * kd + 119 ] += -1 ;
    nuki[ 17 * kd + 119 ] += -1 ;
    nuki[ 16 * kd + 119 ] += +1 ;
    nuki[ 7 * kd + 119 ] += +1 ;

    /*reaction 121: C + O2 <=> O + CO */
    nuki[ 8 * kd + 120 ] += -1 ;
    nuki[ 3 * kd + 120 ] += -1 ;
    nuki[ 2 * kd + 120 ] += +1 ;
    nuki[ 14 * kd + 120 ] += +1 ;

    /*reaction 122: C + CH2 <=> H + C2H */
    nuki[ 8 * kd + 121 ] += -1 ;
    nuki[ 10 * kd + 121 ] += -1 ;
    nuki[ 1 * kd + 121 ] += +1 ;
    nuki[ 21 * kd + 121 ] += +1 ;

    /*reaction 123: C + CH3 <=> H + C2H2 */
    nuki[ 8 * kd + 122 ] += -1 ;
    nuki[ 12 * kd + 122 ] += -1 ;
    nuki[ 1 * kd + 122 ] += +1 ;
    nuki[ 22 * kd + 122 ] += +1 ;

    /*reaction 124: CH + O2 <=> O + HCO */
    nuki[ 9 * kd + 123 ] += -1 ;
    nuki[ 3 * kd + 123 ] += -1 ;
    nuki[ 2 * kd + 123 ] += +1 ;
    nuki[ 16 * kd + 123 ] += +1 ;

    /*reaction 125: CH + H2 <=> H + CH2 */
    nuki[ 9 * kd + 124 ] += -1 ;
    nuki[ 0 * kd + 124 ] += -1 ;
    nuki[ 1 * kd + 124 ] += +1 ;
    nuki[ 10 * kd + 124 ] += +1 ;

    /*reaction 126: CH + H2O <=> H + CH2O */
    nuki[ 9 * kd + 125 ] += -1 ;
    nuki[ 5 * kd + 125 ] += -1 ;
    nuki[ 1 * kd + 125 ] += +1 ;
    nuki[ 17 * kd + 125 ] += +1 ;

    /*reaction 127: CH + CH2 <=> H + C2H2 */
    nuki[ 9 * kd + 126 ] += -1 ;
    nuki[ 10 * kd + 126 ] += -1 ;
    nuki[ 1 * kd + 126 ] += +1 ;
    nuki[ 22 * kd + 126 ] += +1 ;

    /*reaction 128: CH + CH3 <=> H + C2H3 */
    nuki[ 9 * kd + 127 ] += -1 ;
    nuki[ 12 * kd + 127 ] += -1 ;
    nuki[ 1 * kd + 127 ] += +1 ;
    nuki[ 23 * kd + 127 ] += +1 ;

    /*reaction 129: CH + CH4 <=> H + C2H4 */
    nuki[ 9 * kd + 128 ] += -1 ;
    nuki[ 13 * kd + 128 ] += -1 ;
    nuki[ 1 * kd + 128 ] += +1 ;
    nuki[ 24 * kd + 128 ] += +1 ;

    /*reaction 130: CH + CO (+M) <=> HCCO (+M) */
    nuki[ 9 * kd + 129 ] += -1 ;
    nuki[ 14 * kd + 129 ] += -1 ;
    nuki[ 27 * kd + 129 ] += +1 ;

    /*reaction 131: CH + CO2 <=> HCO + CO */
    nuki[ 9 * kd + 130 ] += -1 ;
    nuki[ 15 * kd + 130 ] += -1 ;
    nuki[ 16 * kd + 130 ] += +1 ;
    nuki[ 14 * kd + 130 ] += +1 ;

    /*reaction 132: CH + CH2O <=> H + CH2CO */
    nuki[ 9 * kd + 131 ] += -1 ;
    nuki[ 17 * kd + 131 ] += -1 ;
    nuki[ 1 * kd + 131 ] += +1 ;
    nuki[ 28 * kd + 131 ] += +1 ;

    /*reaction 133: CH + HCCO <=> CO + C2H2 */
    nuki[ 9 * kd + 132 ] += -1 ;
    nuki[ 27 * kd + 132 ] += -1 ;
    nuki[ 14 * kd + 132 ] += +1 ;
    nuki[ 22 * kd + 132 ] += +1 ;

    /*reaction 134: CH2 + O2 => OH + H + CO */
    nuki[ 10 * kd + 133 ] += -1 ;
    nuki[ 3 * kd + 133 ] += -1 ;
    nuki[ 4 * kd + 133 ] += +1 ;
    nuki[ 1 * kd + 133 ] += +1 ;
    nuki[ 14 * kd + 133 ] += +1 ;

    /*reaction 135: CH2 + H2 <=> H + CH3 */
    nuki[ 10 * kd + 134 ] += -1 ;
    nuki[ 0 * kd + 134 ] += -1 ;
    nuki[ 1 * kd + 134 ] += +1 ;
    nuki[ 12 * kd + 134 ] += +1 ;

    /*reaction 136: 2 CH2 <=> H2 + C2H2 */
    nuki[ 10 * kd + 135 ] += -2 ;
    nuki[ 0 * kd + 135 ] += +1 ;
    nuki[ 22 * kd + 135 ] += +1 ;

    /*reaction 137: CH2 + CH3 <=> H + C2H4 */
    nuki[ 10 * kd + 136 ] += -1 ;
    nuki[ 12 * kd + 136 ] += -1 ;
    nuki[ 1 * kd + 136 ] += +1 ;
    nuki[ 24 * kd + 136 ] += +1 ;

    /*reaction 138: CH2 + CH4 <=> 2 CH3 */
    nuki[ 10 * kd + 137 ] += -1 ;
    nuki[ 13 * kd + 137 ] += -1 ;
    nuki[ 12 * kd + 137 ] += +2 ;

    /*reaction 139: CH2 + CO (+M) <=> CH2CO (+M) */
    nuki[ 10 * kd + 138 ] += -1 ;
    nuki[ 14 * kd + 138 ] += -1 ;
    nuki[ 28 * kd + 138 ] += +1 ;

    /*reaction 140: CH2 + HCCO <=> C2H3 + CO */
    nuki[ 10 * kd + 139 ] += -1 ;
    nuki[ 27 * kd + 139 ] += -1 ;
    nuki[ 23 * kd + 139 ] += +1 ;
    nuki[ 14 * kd + 139 ] += +1 ;

    /*reaction 141: CH2(S) + N2 <=> CH2 + N2 */
    nuki[ 11 * kd + 140 ] += -1 ;
    nuki[ 30 * kd + 140 ] += -1 ;
    nuki[ 10 * kd + 140 ] += +1 ;
    nuki[ 30 * kd + 140 ] += +1 ;

    /*reaction 142: CH2(S) + O2 <=> H + OH + CO */
    nuki[ 11 * kd + 141 ] += -1 ;
    nuki[ 3 * kd + 141 ] += -1 ;
    nuki[ 1 * kd + 141 ] += +1 ;
    nuki[ 4 * kd + 141 ] += +1 ;
    nuki[ 14 * kd + 141 ] += +1 ;

    /*reaction 143: CH2(S) + O2 <=> CO + H2O */
    nuki[ 11 * kd + 142 ] += -1 ;
    nuki[ 3 * kd + 142 ] += -1 ;
    nuki[ 14 * kd + 142 ] += +1 ;
    nuki[ 5 * kd + 142 ] += +1 ;

    /*reaction 144: CH2(S) + H2 <=> CH3 + H */
    nuki[ 11 * kd + 143 ] += -1 ;
    nuki[ 0 * kd + 143 ] += -1 ;
    nuki[ 12 * kd + 143 ] += +1 ;
    nuki[ 1 * kd + 143 ] += +1 ;

    /*reaction 145: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    nuki[ 11 * kd + 144 ] += -1 ;
    nuki[ 5 * kd + 144 ] += -1 ;
    nuki[ 20 * kd + 144 ] += +1 ;

    /*reaction 146: CH2(S) + H2O <=> CH2 + H2O */
    nuki[ 11 * kd + 145 ] += -1 ;
    nuki[ 5 * kd + 145 ] += -1 ;
    nuki[ 10 * kd + 145 ] += +1 ;
    nuki[ 5 * kd + 145 ] += +1 ;

    /*reaction 147: CH2(S) + CH3 <=> H + C2H4 */
    nuki[ 11 * kd + 146 ] += -1 ;
    nuki[ 12 * kd + 146 ] += -1 ;
    nuki[ 1 * kd + 146 ] += +1 ;
    nuki[ 24 * kd + 146 ] += +1 ;

    /*reaction 148: CH2(S) + CH4 <=> 2 CH3 */
    nuki[ 11 * kd + 147 ] += -1 ;
    nuki[ 13 * kd + 147 ] += -1 ;
    nuki[ 12 * kd + 147 ] += +2 ;

    /*reaction 149: CH2(S) + CO <=> CH2 + CO */
    nuki[ 11 * kd + 148 ] += -1 ;
    nuki[ 14 * kd + 148 ] += -1 ;
    nuki[ 10 * kd + 148 ] += +1 ;
    nuki[ 14 * kd + 148 ] += +1 ;

    /*reaction 150: CH2(S) + CO2 <=> CH2 + CO2 */
    nuki[ 11 * kd + 149 ] += -1 ;
    nuki[ 15 * kd + 149 ] += -1 ;
    nuki[ 10 * kd + 149 ] += +1 ;
    nuki[ 15 * kd + 149 ] += +1 ;

    /*reaction 151: CH2(S) + CO2 <=> CO + CH2O */
    nuki[ 11 * kd + 150 ] += -1 ;
    nuki[ 15 * kd + 150 ] += -1 ;
    nuki[ 14 * kd + 150 ] += +1 ;
    nuki[ 17 * kd + 150 ] += +1 ;

    /*reaction 152: CH2(S) + C2H6 <=> CH3 + C2H5 */
    nuki[ 11 * kd + 151 ] += -1 ;
    nuki[ 26 * kd + 151 ] += -1 ;
    nuki[ 12 * kd + 151 ] += +1 ;
    nuki[ 25 * kd + 151 ] += +1 ;

    /*reaction 153: CH3 + O2 <=> O + CH3O */
    nuki[ 12 * kd + 152 ] += -1 ;
    nuki[ 3 * kd + 152 ] += -1 ;
    nuki[ 2 * kd + 152 ] += +1 ;
    nuki[ 19 * kd + 152 ] += +1 ;

    /*reaction 154: CH3 + O2 <=> OH + CH2O */
    nuki[ 12 * kd + 153 ] += -1 ;
    nuki[ 3 * kd + 153 ] += -1 ;
    nuki[ 4 * kd + 153 ] += +1 ;
    nuki[ 17 * kd + 153 ] += +1 ;

    /*reaction 155: CH3 + H2O2 <=> HO2 + CH4 */
    nuki[ 12 * kd + 154 ] += -1 ;
    nuki[ 7 * kd + 154 ] += -1 ;
    nuki[ 6 * kd + 154 ] += +1 ;
    nuki[ 13 * kd + 154 ] += +1 ;

    /*reaction 156: 2 CH3 (+M) <=> C2H6 (+M) */
    nuki[ 12 * kd + 155 ] += -2 ;
    nuki[ 26 * kd + 155 ] += +1 ;

    /*reaction 157: 2 CH3 <=> H + C2H5 */
    nuki[ 12 * kd + 156 ] += -2 ;
    nuki[ 1 * kd + 156 ] += +1 ;
    nuki[ 25 * kd + 156 ] += +1 ;

    /*reaction 158: CH3 + HCO <=> CH4 + CO */
    nuki[ 12 * kd + 157 ] += -1 ;
    nuki[ 16 * kd + 157 ] += -1 ;
    nuki[ 13 * kd + 157 ] += +1 ;
    nuki[ 14 * kd + 157 ] += +1 ;

    /*reaction 159: CH3 + CH2O <=> HCO + CH4 */
    nuki[ 12 * kd + 158 ] += -1 ;
    nuki[ 17 * kd + 158 ] += -1 ;
    nuki[ 16 * kd + 158 ] += +1 ;
    nuki[ 13 * kd + 158 ] += +1 ;

    /*reaction 160: CH3 + CH3OH <=> CH2OH + CH4 */
    nuki[ 12 * kd + 159 ] += -1 ;
    nuki[ 20 * kd + 159 ] += -1 ;
    nuki[ 18 * kd + 159 ] += +1 ;
    nuki[ 13 * kd + 159 ] += +1 ;

    /*reaction 161: CH3 + CH3OH <=> CH3O + CH4 */
    nuki[ 12 * kd + 160 ] += -1 ;
    nuki[ 20 * kd + 160 ] += -1 ;
    nuki[ 19 * kd + 160 ] += +1 ;
    nuki[ 13 * kd + 160 ] += +1 ;

    /*reaction 162: CH3 + C2H4 <=> C2H3 + CH4 */
    nuki[ 12 * kd + 161 ] += -1 ;
    nuki[ 24 * kd + 161 ] += -1 ;
    nuki[ 23 * kd + 161 ] += +1 ;
    nuki[ 13 * kd + 161 ] += +1 ;

    /*reaction 163: CH3 + C2H6 <=> C2H5 + CH4 */
    nuki[ 12 * kd + 162 ] += -1 ;
    nuki[ 26 * kd + 162 ] += -1 ;
    nuki[ 25 * kd + 162 ] += +1 ;
    nuki[ 13 * kd + 162 ] += +1 ;

    /*reaction 164: HCO + H2O <=> H + CO + H2O */
    nuki[ 16 * kd + 163 ] += -1 ;
    nuki[ 5 * kd + 163 ] += -1 ;
    nuki[ 1 * kd + 163 ] += +1 ;
    nuki[ 14 * kd + 163 ] += +1 ;
    nuki[ 5 * kd + 163 ] += +1 ;

    /*reaction 165: HCO + M <=> H + CO + M */
    nuki[ 16 * kd + 164 ] += -1 ;
    nuki[ 1 * kd + 164 ] += +1 ;
    nuki[ 14 * kd + 164 ] += +1 ;

    /*reaction 166: HCO + O2 <=> HO2 + CO */
    nuki[ 16 * kd + 165 ] += -1 ;
    nuki[ 3 * kd + 165 ] += -1 ;
    nuki[ 6 * kd + 165 ] += +1 ;
    nuki[ 14 * kd + 165 ] += +1 ;

    /*reaction 167: CH2OH + O2 <=> HO2 + CH2O */
    nuki[ 18 * kd + 166 ] += -1 ;
    nuki[ 3 * kd + 166 ] += -1 ;
    nuki[ 6 * kd + 166 ] += +1 ;
    nuki[ 17 * kd + 166 ] += +1 ;

    /*reaction 168: CH3O + O2 <=> HO2 + CH2O */
    nuki[ 19 * kd + 167 ] += -1 ;
    nuki[ 3 * kd + 167 ] += -1 ;
    nuki[ 6 * kd + 167 ] += +1 ;
    nuki[ 17 * kd + 167 ] += +1 ;

    /*reaction 169: C2H + O2 <=> HCO + CO */
    nuki[ 21 * kd + 168 ] += -1 ;
    nuki[ 3 * kd + 168 ] += -1 ;
    nuki[ 16 * kd + 168 ] += +1 ;
    nuki[ 14 * kd + 168 ] += +1 ;

    /*reaction 170: C2H + H2 <=> H + C2H2 */
    nuki[ 21 * kd + 169 ] += -1 ;
    nuki[ 0 * kd + 169 ] += -1 ;
    nuki[ 1 * kd + 169 ] += +1 ;
    nuki[ 22 * kd + 169 ] += +1 ;

    /*reaction 171: C2H3 + O2 <=> HCO + CH2O */
    nuki[ 23 * kd + 170 ] += -1 ;
    nuki[ 3 * kd + 170 ] += -1 ;
    nuki[ 16 * kd + 170 ] += +1 ;
    nuki[ 17 * kd + 170 ] += +1 ;

    /*reaction 172: C2H4 (+M) <=> H2 + C2H2 (+M) */
    nuki[ 24 * kd + 171 ] += -1 ;
    nuki[ 0 * kd + 171 ] += +1 ;
    nuki[ 22 * kd + 171 ] += +1 ;

    /*reaction 173: C2H5 + O2 <=> HO2 + C2H4 */
    nuki[ 25 * kd + 172 ] += -1 ;
    nuki[ 3 * kd + 172 ] += -1 ;
    nuki[ 6 * kd + 172 ] += +1 ;
    nuki[ 24 * kd + 172 ] += +1 ;

    /*reaction 174: HCCO + O2 <=> OH + 2 CO */
    nuki[ 27 * kd + 173 ] += -1 ;
    nuki[ 3 * kd + 173 ] += -1 ;
    nuki[ 4 * kd + 173 ] += +1 ;
    nuki[ 14 * kd + 173 ] += +2 ;

    /*reaction 175: 2 HCCO <=> 2 CO + C2H2 */
    nuki[ 27 * kd + 174 ] += -2 ;
    nuki[ 14 * kd + 174 ] += +2 ;
    nuki[ 22 * kd + 174 ] += +1 ;

    /*reaction 176: O + CH3 => H + H2 + CO */
    nuki[ 2 * kd + 175 ] += -1 ;
    nuki[ 12 * kd + 175 ] += -1 ;
    nuki[ 1 * kd + 175 ] += +1 ;
    nuki[ 0 * kd + 175 ] += +1 ;
    nuki[ 14 * kd + 175 ] += +1 ;

    /*reaction 177: O + C2H4 <=> H + CH2CHO */
    nuki[ 2 * kd + 176 ] += -1 ;
    nuki[ 24 * kd + 176 ] += -1 ;
    nuki[ 1 * kd + 176 ] += +1 ;
    nuki[ 33 * kd + 176 ] += +1 ;

    /*reaction 178: O + C2H5 <=> H + CH3CHO */
    nuki[ 2 * kd + 177 ] += -1 ;
    nuki[ 25 * kd + 177 ] += -1 ;
    nuki[ 1 * kd + 177 ] += +1 ;
    nuki[ 34 * kd + 177 ] += +1 ;

    /*reaction 179: OH + HO2 <=> O2 + H2O */
    nuki[ 4 * kd + 178 ] += -1 ;
    nuki[ 6 * kd + 178 ] += -1 ;
    nuki[ 3 * kd + 178 ] += +1 ;
    nuki[ 5 * kd + 178 ] += +1 ;

    /*reaction 180: OH + CH3 => H2 + CH2O */
    nuki[ 4 * kd + 179 ] += -1 ;
    nuki[ 12 * kd + 179 ] += -1 ;
    nuki[ 0 * kd + 179 ] += +1 ;
    nuki[ 17 * kd + 179 ] += +1 ;

    /*reaction 181: CH + H2 (+M) <=> CH3 (+M) */
    nuki[ 9 * kd + 180 ] += -1 ;
    nuki[ 0 * kd + 180 ] += -1 ;
    nuki[ 12 * kd + 180 ] += +1 ;

    /*reaction 182: CH2 + O2 => 2 H + CO2 */
    nuki[ 10 * kd + 181 ] += -1 ;
    nuki[ 3 * kd + 181 ] += -1 ;
    nuki[ 1 * kd + 181 ] += +2 ;
    nuki[ 15 * kd + 181 ] += +1 ;

    /*reaction 183: CH2 + O2 <=> O + CH2O */
    nuki[ 10 * kd + 182 ] += -1 ;
    nuki[ 3 * kd + 182 ] += -1 ;
    nuki[ 2 * kd + 182 ] += +1 ;
    nuki[ 17 * kd + 182 ] += +1 ;

    /*reaction 184: CH2 + CH2 => 2 H + C2H2 */
    nuki[ 10 * kd + 183 ] += -1 ;
    nuki[ 10 * kd + 183 ] += -1 ;
    nuki[ 1 * kd + 183 ] += +2 ;
    nuki[ 22 * kd + 183 ] += +1 ;

    /*reaction 185: CH2(S) + H2O => H2 + CH2O */
    nuki[ 11 * kd + 184 ] += -1 ;
    nuki[ 5 * kd + 184 ] += -1 ;
    nuki[ 0 * kd + 184 ] += +1 ;
    nuki[ 17 * kd + 184 ] += +1 ;

    /*reaction 186: C2H3 + O2 <=> O + CH2CHO */
    nuki[ 23 * kd + 185 ] += -1 ;
    nuki[ 3 * kd + 185 ] += -1 ;
    nuki[ 2 * kd + 185 ] += +1 ;
    nuki[ 33 * kd + 185 ] += +1 ;

    /*reaction 187: C2H3 + O2 <=> HO2 + C2H2 */
    nuki[ 23 * kd + 186 ] += -1 ;
    nuki[ 3 * kd + 186 ] += -1 ;
    nuki[ 6 * kd + 186 ] += +1 ;
    nuki[ 22 * kd + 186 ] += +1 ;

    /*reaction 188: O + CH3CHO <=> OH + CH2CHO */
    nuki[ 2 * kd + 187 ] += -1 ;
    nuki[ 34 * kd + 187 ] += -1 ;
    nuki[ 4 * kd + 187 ] += +1 ;
    nuki[ 33 * kd + 187 ] += +1 ;

    /*reaction 189: O + CH3CHO => OH + CH3 + CO */
    nuki[ 2 * kd + 188 ] += -1 ;
    nuki[ 34 * kd + 188 ] += -1 ;
    nuki[ 4 * kd + 188 ] += +1 ;
    nuki[ 12 * kd + 188 ] += +1 ;
    nuki[ 14 * kd + 188 ] += +1 ;

    /*reaction 190: O2 + CH3CHO => HO2 + CH3 + CO */
    nuki[ 3 * kd + 189 ] += -1 ;
    nuki[ 34 * kd + 189 ] += -1 ;
    nuki[ 6 * kd + 189 ] += +1 ;
    nuki[ 12 * kd + 189 ] += +1 ;
    nuki[ 14 * kd + 189 ] += +1 ;

    /*reaction 191: H + CH3CHO <=> CH2CHO + H2 */
    nuki[ 1 * kd + 190 ] += -1 ;
    nuki[ 34 * kd + 190 ] += -1 ;
    nuki[ 33 * kd + 190 ] += +1 ;
    nuki[ 0 * kd + 190 ] += +1 ;

    /*reaction 192: H + CH3CHO => CH3 + H2 + CO */
    nuki[ 1 * kd + 191 ] += -1 ;
    nuki[ 34 * kd + 191 ] += -1 ;
    nuki[ 12 * kd + 191 ] += +1 ;
    nuki[ 0 * kd + 191 ] += +1 ;
    nuki[ 14 * kd + 191 ] += +1 ;

    /*reaction 193: OH + CH3CHO => CH3 + H2O + CO */
    nuki[ 4 * kd + 192 ] += -1 ;
    nuki[ 34 * kd + 192 ] += -1 ;
    nuki[ 12 * kd + 192 ] += +1 ;
    nuki[ 5 * kd + 192 ] += +1 ;
    nuki[ 14 * kd + 192 ] += +1 ;

    /*reaction 194: HO2 + CH3CHO => CH3 + H2O2 + CO */
    nuki[ 6 * kd + 193 ] += -1 ;
    nuki[ 34 * kd + 193 ] += -1 ;
    nuki[ 12 * kd + 193 ] += +1 ;
    nuki[ 7 * kd + 193 ] += +1 ;
    nuki[ 14 * kd + 193 ] += +1 ;

    /*reaction 195: CH3 + CH3CHO => CH3 + CH4 + CO */
    nuki[ 12 * kd + 194 ] += -1 ;
    nuki[ 34 * kd + 194 ] += -1 ;
    nuki[ 12 * kd + 194 ] += +1 ;
    nuki[ 13 * kd + 194 ] += +1 ;
    nuki[ 14 * kd + 194 ] += +1 ;

    /*reaction 196: H + CH2CO (+M) <=> CH2CHO (+M) */
    nuki[ 1 * kd + 195 ] += -1 ;
    nuki[ 28 * kd + 195 ] += -1 ;
    nuki[ 33 * kd + 195 ] += +1 ;

    /*reaction 197: O + CH2CHO => H + CH2 + CO2 */
    nuki[ 2 * kd + 196 ] += -1 ;
    nuki[ 33 * kd + 196 ] += -1 ;
    nuki[ 1 * kd + 196 ] += +1 ;
    nuki[ 10 * kd + 196 ] += +1 ;
    nuki[ 15 * kd + 196 ] += +1 ;

    /*reaction 198: O2 + CH2CHO => OH + CO + CH2O */
    nuki[ 3 * kd + 197 ] += -1 ;
    nuki[ 33 * kd + 197 ] += -1 ;
    nuki[ 4 * kd + 197 ] += +1 ;
    nuki[ 14 * kd + 197 ] += +1 ;
    nuki[ 17 * kd + 197 ] += +1 ;

    /*reaction 199: O2 + CH2CHO => OH + 2 HCO */
    nuki[ 3 * kd + 198 ] += -1 ;
    nuki[ 33 * kd + 198 ] += -1 ;
    nuki[ 4 * kd + 198 ] += +1 ;
    nuki[ 16 * kd + 198 ] += +2 ;

    /*reaction 200: H + CH2CHO <=> CH3 + HCO */
    nuki[ 1 * kd + 199 ] += -1 ;
    nuki[ 33 * kd + 199 ] += -1 ;
    nuki[ 12 * kd + 199 ] += +1 ;
    nuki[ 16 * kd + 199 ] += +1 ;

    /*reaction 201: H + CH2CHO <=> CH2CO + H2 */
    nuki[ 1 * kd + 200 ] += -1 ;
    nuki[ 33 * kd + 200 ] += -1 ;
    nuki[ 28 * kd + 200 ] += +1 ;
    nuki[ 0 * kd + 200 ] += +1 ;

    /*reaction 202: OH + CH2CHO <=> H2O + CH2CO */
    nuki[ 4 * kd + 201 ] += -1 ;
    nuki[ 33 * kd + 201 ] += -1 ;
    nuki[ 5 * kd + 201 ] += +1 ;
    nuki[ 28 * kd + 201 ] += +1 ;

    /*reaction 203: OH + CH2CHO <=> HCO + CH2OH */
    nuki[ 4 * kd + 202 ] += -1 ;
    nuki[ 33 * kd + 202 ] += -1 ;
    nuki[ 16 * kd + 202 ] += +1 ;
    nuki[ 18 * kd + 202 ] += +1 ;

    /*reaction 204: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    nuki[ 12 * kd + 203 ] += -1 ;
    nuki[ 25 * kd + 203 ] += -1 ;
    nuki[ 32 * kd + 203 ] += +1 ;

    /*reaction 205: O + C3H8 <=> OH + C3H7 */
    nuki[ 2 * kd + 204 ] += -1 ;
    nuki[ 32 * kd + 204 ] += -1 ;
    nuki[ 4 * kd + 204 ] += +1 ;
    nuki[ 31 * kd + 204 ] += +1 ;

    /*reaction 206: H + C3H8 <=> C3H7 + H2 */
    nuki[ 1 * kd + 205 ] += -1 ;
    nuki[ 32 * kd + 205 ] += -1 ;
    nuki[ 31 * kd + 205 ] += +1 ;
    nuki[ 0 * kd + 205 ] += +1 ;

    /*reaction 207: OH + C3H8 <=> C3H7 + H2O */
    nuki[ 4 * kd + 206 ] += -1 ;
    nuki[ 32 * kd + 206 ] += -1 ;
    nuki[ 31 * kd + 206 ] += +1 ;
    nuki[ 5 * kd + 206 ] += +1 ;

    /*reaction 208: C3H7 + H2O2 <=> HO2 + C3H8 */
    nuki[ 31 * kd + 207 ] += -1 ;
    nuki[ 7 * kd + 207 ] += -1 ;
    nuki[ 6 * kd + 207 ] += +1 ;
    nuki[ 32 * kd + 207 ] += +1 ;

    /*reaction 209: CH3 + C3H8 <=> C3H7 + CH4 */
    nuki[ 12 * kd + 208 ] += -1 ;
    nuki[ 32 * kd + 208 ] += -1 ;
    nuki[ 31 * kd + 208 ] += +1 ;
    nuki[ 13 * kd + 208 ] += +1 ;

    /*reaction 210: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    nuki[ 12 * kd + 209 ] += -1 ;
    nuki[ 24 * kd + 209 ] += -1 ;
    nuki[ 31 * kd + 209 ] += +1 ;

    /*reaction 211: O + C3H7 <=> C2H5 + CH2O */
    nuki[ 2 * kd + 210 ] += -1 ;
    nuki[ 31 * kd + 210 ] += -1 ;
    nuki[ 25 * kd + 210 ] += +1 ;
    nuki[ 17 * kd + 210 ] += +1 ;

    /*reaction 212: H + C3H7 (+M) <=> C3H8 (+M) */
    nuki[ 1 * kd + 211 ] += -1 ;
    nuki[ 31 * kd + 211 ] += -1 ;
    nuki[ 32 * kd + 211 ] += +1 ;

    /*reaction 213: H + C3H7 <=> CH3 + C2H5 */
    nuki[ 1 * kd + 212 ] += -1 ;
    nuki[ 31 * kd + 212 ] += -1 ;
    nuki[ 12 * kd + 212 ] += +1 ;
    nuki[ 25 * kd + 212 ] += +1 ;

    /*reaction 214: OH + C3H7 <=> C2H5 + CH2OH */
    nuki[ 4 * kd + 213 ] += -1 ;
    nuki[ 31 * kd + 213 ] += -1 ;
    nuki[ 25 * kd + 213 ] += +1 ;
    nuki[ 18 * kd + 213 ] += +1 ;

    /*reaction 215: HO2 + C3H7 <=> O2 + C3H8 */
    nuki[ 6 * kd + 214 ] += -1 ;
    nuki[ 31 * kd + 214 ] += -1 ;
    nuki[ 3 * kd + 214 ] += +1 ;
    nuki[ 32 * kd + 214 ] += +1 ;

    /*reaction 216: HO2 + C3H7 => OH + C2H5 + CH2O */
    nuki[ 6 * kd + 215 ] += -1 ;
    nuki[ 31 * kd + 215 ] += -1 ;
    nuki[ 4 * kd + 215 ] += +1 ;
    nuki[ 25 * kd + 215 ] += +1 ;
    nuki[ 17 * kd + 215 ] += +1 ;

    /*reaction 217: CH3 + C3H7 <=> 2 C2H5 */
    nuki[ 12 * kd + 216 ] += -1 ;
    nuki[ 31 * kd + 216 ] += -1 ;
    nuki[ 25 * kd + 216 ] += +2 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < 4 * 35; ++ id) {
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

    /*C3H7 */
    ncf[ 31 * kd + 2 ] = 3; /*C */
    ncf[ 31 * kd + 1 ] = 7; /*H */

    /*C3H8 */
    ncf[ 32 * kd + 2 ] = 3; /*C */
    ncf[ 32 * kd + 1 ] = 8; /*H */

    /*CH2CHO */
    ncf[ 33 * kd + 0 ] = 1; /*O */
    ncf[ 33 * kd + 1 ] = 3; /*H */
    ncf[ 33 * kd + 2 ] = 2; /*C */

    /*CH3CHO */
    ncf[ 34 * kd + 2 ] = 2; /*C */
    ncf[ 34 * kd + 1 ] = 4; /*H */
    ncf[ 34 * kd + 0 ] = 1; /*O */

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

    /*reaction 37: H + O2 <=> O + OH */
    a[36] = 2.65e+16;
    b[36] = -0.6707;
    e[36] = 17041;

    /*reaction 38: 2 H + M <=> H2 + M */
    a[37] = 1e+18;
    b[37] = -1;
    e[37] = 0;

    /*reaction 39: 2 H + H2 <=> 2 H2 */
    a[38] = 9e+16;
    b[38] = -0.6;
    e[38] = 0;

    /*reaction 40: 2 H + H2O <=> H2 + H2O */
    a[39] = 6e+19;
    b[39] = -1.25;
    e[39] = 0;

    /*reaction 41: 2 H + CO2 <=> H2 + CO2 */
    a[40] = 5.5e+20;
    b[40] = -2;
    e[40] = 0;

    /*reaction 42: H + OH + M <=> H2O + M */
    a[41] = 2.2e+22;
    b[41] = -2;
    e[41] = 0;

    /*reaction 43: H + HO2 <=> O + H2O */
    a[42] = 3.97e+12;
    b[42] = 0;
    e[42] = 671;

    /*reaction 44: H + HO2 <=> O2 + H2 */
    a[43] = 4.48e+13;
    b[43] = 0;
    e[43] = 1068;

    /*reaction 45: H + HO2 <=> 2 OH */
    a[44] = 8.4e+13;
    b[44] = 0;
    e[44] = 635;

    /*reaction 46: H + H2O2 <=> HO2 + H2 */
    a[45] = 1.21e+07;
    b[45] = 2;
    e[45] = 5200;

    /*reaction 47: H + H2O2 <=> OH + H2O */
    a[46] = 1e+13;
    b[46] = 0;
    e[46] = 3600;

    /*reaction 48: H + CH <=> C + H2 */
    a[47] = 1.65e+14;
    b[47] = 0;
    e[47] = 0;

    /*reaction 49: H + CH2 (+M) <=> CH3 (+M) */
    a[48] = 6e+14;
    b[48] = 0;
    e[48] = 0;

    /*reaction 50: H + CH2(S) <=> CH + H2 */
    a[49] = 3e+13;
    b[49] = 0;
    e[49] = 0;

    /*reaction 51: H + CH3 (+M) <=> CH4 (+M) */
    a[50] = 1.39e+16;
    b[50] = -0.534;
    e[50] = 536;

    /*reaction 52: H + CH4 <=> CH3 + H2 */
    a[51] = 6.6e+08;
    b[51] = 1.62;
    e[51] = 10840;

    /*reaction 53: H + HCO (+M) <=> CH2O (+M) */
    a[52] = 1.09e+12;
    b[52] = 0.48;
    e[52] = -260;

    /*reaction 54: H + HCO <=> H2 + CO */
    a[53] = 7.34e+13;
    b[53] = 0;
    e[53] = 0;

    /*reaction 55: H + CH2O (+M) <=> CH2OH (+M) */
    a[54] = 5.4e+11;
    b[54] = 0.454;
    e[54] = 3600;

    /*reaction 56: H + CH2O (+M) <=> CH3O (+M) */
    a[55] = 5.4e+11;
    b[55] = 0.454;
    e[55] = 2600;

    /*reaction 57: H + CH2O <=> HCO + H2 */
    a[56] = 5.74e+07;
    b[56] = 1.9;
    e[56] = 2742;

    /*reaction 58: H + CH2OH (+M) <=> CH3OH (+M) */
    a[57] = 1.055e+12;
    b[57] = 0.5;
    e[57] = 86;

    /*reaction 59: H + CH2OH <=> H2 + CH2O */
    a[58] = 2e+13;
    b[58] = 0;
    e[58] = 0;

    /*reaction 60: H + CH2OH <=> OH + CH3 */
    a[59] = 1.65e+11;
    b[59] = 0.65;
    e[59] = -284;

    /*reaction 61: H + CH2OH <=> CH2(S) + H2O */
    a[60] = 3.28e+13;
    b[60] = -0.09;
    e[60] = 610;

    /*reaction 62: H + CH3O (+M) <=> CH3OH (+M) */
    a[61] = 2.43e+12;
    b[61] = 0.515;
    e[61] = 50;

    /*reaction 63: H + CH3O <=> H + CH2OH */
    a[62] = 4.15e+07;
    b[62] = 1.63;
    e[62] = 1924;

    /*reaction 64: H + CH3O <=> H2 + CH2O */
    a[63] = 2e+13;
    b[63] = 0;
    e[63] = 0;

    /*reaction 65: H + CH3O <=> OH + CH3 */
    a[64] = 1.5e+12;
    b[64] = 0.5;
    e[64] = -110;

    /*reaction 66: H + CH3O <=> CH2(S) + H2O */
    a[65] = 2.62e+14;
    b[65] = -0.23;
    e[65] = 1070;

    /*reaction 67: H + CH3OH <=> CH2OH + H2 */
    a[66] = 1.7e+07;
    b[66] = 2.1;
    e[66] = 4870;

    /*reaction 68: H + CH3OH <=> CH3O + H2 */
    a[67] = 4.2e+06;
    b[67] = 2.1;
    e[67] = 4870;

    /*reaction 69: H + C2H (+M) <=> C2H2 (+M) */
    a[68] = 1e+17;
    b[68] = -1;
    e[68] = 0;

    /*reaction 70: H + C2H2 (+M) <=> C2H3 (+M) */
    a[69] = 5.6e+12;
    b[69] = 0;
    e[69] = 2400;

    /*reaction 71: H + C2H3 (+M) <=> C2H4 (+M) */
    a[70] = 6.08e+12;
    b[70] = 0.27;
    e[70] = 280;

    /*reaction 72: H + C2H3 <=> H2 + C2H2 */
    a[71] = 3e+13;
    b[71] = 0;
    e[71] = 0;

    /*reaction 73: H + C2H4 (+M) <=> C2H5 (+M) */
    a[72] = 5.4e+11;
    b[72] = 0.454;
    e[72] = 1820;

    /*reaction 74: H + C2H4 <=> C2H3 + H2 */
    a[73] = 1.325e+06;
    b[73] = 2.53;
    e[73] = 12240;

    /*reaction 75: H + C2H5 (+M) <=> C2H6 (+M) */
    a[74] = 5.21e+17;
    b[74] = -0.99;
    e[74] = 1580;

    /*reaction 76: H + C2H5 <=> H2 + C2H4 */
    a[75] = 2e+12;
    b[75] = 0;
    e[75] = 0;

    /*reaction 77: H + C2H6 <=> C2H5 + H2 */
    a[76] = 1.15e+08;
    b[76] = 1.9;
    e[76] = 7530;

    /*reaction 78: H + HCCO <=> CH2(S) + CO */
    a[77] = 1e+14;
    b[77] = 0;
    e[77] = 0;

    /*reaction 79: H + CH2CO <=> HCCO + H2 */
    a[78] = 5e+13;
    b[78] = 0;
    e[78] = 8000;

    /*reaction 80: H + CH2CO <=> CH3 + CO */
    a[79] = 1.13e+13;
    b[79] = 0;
    e[79] = 3428;

    /*reaction 81: H + HCCOH <=> H + CH2CO */
    a[80] = 1e+13;
    b[80] = 0;
    e[80] = 0;

    /*reaction 82: H2 + CO (+M) <=> CH2O (+M) */
    a[81] = 4.3e+07;
    b[81] = 1.5;
    e[81] = 79600;

    /*reaction 83: OH + H2 <=> H + H2O */
    a[82] = 2.16e+08;
    b[82] = 1.51;
    e[82] = 3430;

    /*reaction 84: 2 OH (+M) <=> H2O2 (+M) */
    a[83] = 7.4e+13;
    b[83] = -0.37;
    e[83] = 0;

    /*reaction 85: 2 OH <=> O + H2O */
    a[84] = 35700;
    b[84] = 2.4;
    e[84] = -2110;

    /*reaction 86: OH + HO2 <=> O2 + H2O */
    a[85] = 1.45e+13;
    b[85] = 0;
    e[85] = -500;

    /*reaction 87: OH + H2O2 <=> HO2 + H2O */
    a[86] = 2e+12;
    b[86] = 0;
    e[86] = 427;

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    a[87] = 1.7e+18;
    b[87] = 0;
    e[87] = 29410;

    /*reaction 89: OH + C <=> H + CO */
    a[88] = 5e+13;
    b[88] = 0;
    e[88] = 0;

    /*reaction 90: OH + CH <=> H + HCO */
    a[89] = 3e+13;
    b[89] = 0;
    e[89] = 0;

    /*reaction 91: OH + CH2 <=> H + CH2O */
    a[90] = 2e+13;
    b[90] = 0;
    e[90] = 0;

    /*reaction 92: OH + CH2 <=> CH + H2O */
    a[91] = 1.13e+07;
    b[91] = 2;
    e[91] = 3000;

    /*reaction 93: OH + CH2(S) <=> H + CH2O */
    a[92] = 3e+13;
    b[92] = 0;
    e[92] = 0;

    /*reaction 94: OH + CH3 (+M) <=> CH3OH (+M) */
    a[93] = 2.79e+18;
    b[93] = -1.43;
    e[93] = 1330;

    /*reaction 95: OH + CH3 <=> CH2 + H2O */
    a[94] = 5.6e+07;
    b[94] = 1.6;
    e[94] = 5420;

    /*reaction 96: OH + CH3 <=> CH2(S) + H2O */
    a[95] = 6.44e+17;
    b[95] = -1.34;
    e[95] = 1417;

    /*reaction 97: OH + CH4 <=> CH3 + H2O */
    a[96] = 1e+08;
    b[96] = 1.6;
    e[96] = 3120;

    /*reaction 98: OH + CO <=> H + CO2 */
    a[97] = 4.76e+07;
    b[97] = 1.228;
    e[97] = 70;

    /*reaction 99: OH + HCO <=> H2O + CO */
    a[98] = 5e+13;
    b[98] = 0;
    e[98] = 0;

    /*reaction 100: OH + CH2O <=> HCO + H2O */
    a[99] = 3.43e+09;
    b[99] = 1.18;
    e[99] = -447;

    /*reaction 101: OH + CH2OH <=> H2O + CH2O */
    a[100] = 5e+12;
    b[100] = 0;
    e[100] = 0;

    /*reaction 102: OH + CH3O <=> H2O + CH2O */
    a[101] = 5e+12;
    b[101] = 0;
    e[101] = 0;

    /*reaction 103: OH + CH3OH <=> CH2OH + H2O */
    a[102] = 1.44e+06;
    b[102] = 2;
    e[102] = -840;

    /*reaction 104: OH + CH3OH <=> CH3O + H2O */
    a[103] = 6.3e+06;
    b[103] = 2;
    e[103] = 1500;

    /*reaction 105: OH + C2H <=> H + HCCO */
    a[104] = 2e+13;
    b[104] = 0;
    e[104] = 0;

    /*reaction 106: OH + C2H2 <=> H + CH2CO */
    a[105] = 0.000218;
    b[105] = 4.5;
    e[105] = -1000;

    /*reaction 107: OH + C2H2 <=> H + HCCOH */
    a[106] = 504000;
    b[106] = 2.3;
    e[106] = 13500;

    /*reaction 108: OH + C2H2 <=> C2H + H2O */
    a[107] = 3.37e+07;
    b[107] = 2;
    e[107] = 14000;

    /*reaction 109: OH + C2H2 <=> CH3 + CO */
    a[108] = 0.000483;
    b[108] = 4;
    e[108] = -2000;

    /*reaction 110: OH + C2H3 <=> H2O + C2H2 */
    a[109] = 5e+12;
    b[109] = 0;
    e[109] = 0;

    /*reaction 111: OH + C2H4 <=> C2H3 + H2O */
    a[110] = 3.6e+06;
    b[110] = 2;
    e[110] = 2500;

    /*reaction 112: OH + C2H6 <=> C2H5 + H2O */
    a[111] = 3.54e+06;
    b[111] = 2.12;
    e[111] = 870;

    /*reaction 113: OH + CH2CO <=> HCCO + H2O */
    a[112] = 7.5e+12;
    b[112] = 0;
    e[112] = 2000;

    /*reaction 114: 2 HO2 <=> O2 + H2O2 */
    a[113] = 1.3e+11;
    b[113] = 0;
    e[113] = -1630;

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    a[114] = 4.2e+14;
    b[114] = 0;
    e[114] = 12000;

    /*reaction 116: HO2 + CH2 <=> OH + CH2O */
    a[115] = 2e+13;
    b[115] = 0;
    e[115] = 0;

    /*reaction 117: HO2 + CH3 <=> O2 + CH4 */
    a[116] = 1e+12;
    b[116] = 0;
    e[116] = 0;

    /*reaction 118: HO2 + CH3 <=> OH + CH3O */
    a[117] = 3.78e+13;
    b[117] = 0;
    e[117] = 0;

    /*reaction 119: HO2 + CO <=> OH + CO2 */
    a[118] = 1.5e+14;
    b[118] = 0;
    e[118] = 23600;

    /*reaction 120: HO2 + CH2O <=> HCO + H2O2 */
    a[119] = 5.6e+06;
    b[119] = 2;
    e[119] = 12000;

    /*reaction 121: C + O2 <=> O + CO */
    a[120] = 5.8e+13;
    b[120] = 0;
    e[120] = 576;

    /*reaction 122: C + CH2 <=> H + C2H */
    a[121] = 5e+13;
    b[121] = 0;
    e[121] = 0;

    /*reaction 123: C + CH3 <=> H + C2H2 */
    a[122] = 5e+13;
    b[122] = 0;
    e[122] = 0;

    /*reaction 124: CH + O2 <=> O + HCO */
    a[123] = 6.71e+13;
    b[123] = 0;
    e[123] = 0;

    /*reaction 125: CH + H2 <=> H + CH2 */
    a[124] = 1.08e+14;
    b[124] = 0;
    e[124] = 3110;

    /*reaction 126: CH + H2O <=> H + CH2O */
    a[125] = 5.71e+12;
    b[125] = 0;
    e[125] = -755;

    /*reaction 127: CH + CH2 <=> H + C2H2 */
    a[126] = 4e+13;
    b[126] = 0;
    e[126] = 0;

    /*reaction 128: CH + CH3 <=> H + C2H3 */
    a[127] = 3e+13;
    b[127] = 0;
    e[127] = 0;

    /*reaction 129: CH + CH4 <=> H + C2H4 */
    a[128] = 6e+13;
    b[128] = 0;
    e[128] = 0;

    /*reaction 130: CH + CO (+M) <=> HCCO (+M) */
    a[129] = 5e+13;
    b[129] = 0;
    e[129] = 0;

    /*reaction 131: CH + CO2 <=> HCO + CO */
    a[130] = 1.9e+14;
    b[130] = 0;
    e[130] = 15792;

    /*reaction 132: CH + CH2O <=> H + CH2CO */
    a[131] = 9.46e+13;
    b[131] = 0;
    e[131] = -515;

    /*reaction 133: CH + HCCO <=> CO + C2H2 */
    a[132] = 5e+13;
    b[132] = 0;
    e[132] = 0;

    /*reaction 134: CH2 + O2 => OH + H + CO */
    a[133] = 5e+12;
    b[133] = 0;
    e[133] = 1500;

    /*reaction 135: CH2 + H2 <=> H + CH3 */
    a[134] = 500000;
    b[134] = 2;
    e[134] = 7230;

    /*reaction 136: 2 CH2 <=> H2 + C2H2 */
    a[135] = 1.6e+15;
    b[135] = 0;
    e[135] = 11944;

    /*reaction 137: CH2 + CH3 <=> H + C2H4 */
    a[136] = 4e+13;
    b[136] = 0;
    e[136] = 0;

    /*reaction 138: CH2 + CH4 <=> 2 CH3 */
    a[137] = 2.46e+06;
    b[137] = 2;
    e[137] = 8270;

    /*reaction 139: CH2 + CO (+M) <=> CH2CO (+M) */
    a[138] = 8.1e+11;
    b[138] = 0.5;
    e[138] = 4510;

    /*reaction 140: CH2 + HCCO <=> C2H3 + CO */
    a[139] = 3e+13;
    b[139] = 0;
    e[139] = 0;

    /*reaction 141: CH2(S) + N2 <=> CH2 + N2 */
    a[140] = 1.5e+13;
    b[140] = 0;
    e[140] = 600;

    /*reaction 142: CH2(S) + O2 <=> H + OH + CO */
    a[141] = 2.8e+13;
    b[141] = 0;
    e[141] = 0;

    /*reaction 143: CH2(S) + O2 <=> CO + H2O */
    a[142] = 1.2e+13;
    b[142] = 0;
    e[142] = 0;

    /*reaction 144: CH2(S) + H2 <=> CH3 + H */
    a[143] = 7e+13;
    b[143] = 0;
    e[143] = 0;

    /*reaction 145: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    a[144] = 4.82e+17;
    b[144] = -1.16;
    e[144] = 1145;

    /*reaction 146: CH2(S) + H2O <=> CH2 + H2O */
    a[145] = 3e+13;
    b[145] = 0;
    e[145] = 0;

    /*reaction 147: CH2(S) + CH3 <=> H + C2H4 */
    a[146] = 1.2e+13;
    b[146] = 0;
    e[146] = -570;

    /*reaction 148: CH2(S) + CH4 <=> 2 CH3 */
    a[147] = 1.6e+13;
    b[147] = 0;
    e[147] = -570;

    /*reaction 149: CH2(S) + CO <=> CH2 + CO */
    a[148] = 9e+12;
    b[148] = 0;
    e[148] = 0;

    /*reaction 150: CH2(S) + CO2 <=> CH2 + CO2 */
    a[149] = 7e+12;
    b[149] = 0;
    e[149] = 0;

    /*reaction 151: CH2(S) + CO2 <=> CO + CH2O */
    a[150] = 1.4e+13;
    b[150] = 0;
    e[150] = 0;

    /*reaction 152: CH2(S) + C2H6 <=> CH3 + C2H5 */
    a[151] = 4e+13;
    b[151] = 0;
    e[151] = -550;

    /*reaction 153: CH3 + O2 <=> O + CH3O */
    a[152] = 3.56e+13;
    b[152] = 0;
    e[152] = 30480;

    /*reaction 154: CH3 + O2 <=> OH + CH2O */
    a[153] = 2.31e+12;
    b[153] = 0;
    e[153] = 20315;

    /*reaction 155: CH3 + H2O2 <=> HO2 + CH4 */
    a[154] = 24500;
    b[154] = 2.47;
    e[154] = 5180;

    /*reaction 156: 2 CH3 (+M) <=> C2H6 (+M) */
    a[155] = 6.77e+16;
    b[155] = -1.18;
    e[155] = 654;

    /*reaction 157: 2 CH3 <=> H + C2H5 */
    a[156] = 6.84e+12;
    b[156] = 0.1;
    e[156] = 10600;

    /*reaction 158: CH3 + HCO <=> CH4 + CO */
    a[157] = 2.648e+13;
    b[157] = 0;
    e[157] = 0;

    /*reaction 159: CH3 + CH2O <=> HCO + CH4 */
    a[158] = 3320;
    b[158] = 2.81;
    e[158] = 5860;

    /*reaction 160: CH3 + CH3OH <=> CH2OH + CH4 */
    a[159] = 3e+07;
    b[159] = 1.5;
    e[159] = 9940;

    /*reaction 161: CH3 + CH3OH <=> CH3O + CH4 */
    a[160] = 1e+07;
    b[160] = 1.5;
    e[160] = 9940;

    /*reaction 162: CH3 + C2H4 <=> C2H3 + CH4 */
    a[161] = 227000;
    b[161] = 2;
    e[161] = 9200;

    /*reaction 163: CH3 + C2H6 <=> C2H5 + CH4 */
    a[162] = 6.14e+06;
    b[162] = 1.74;
    e[162] = 10450;

    /*reaction 164: HCO + H2O <=> H + CO + H2O */
    a[163] = 1.5e+18;
    b[163] = -1;
    e[163] = 17000;

    /*reaction 165: HCO + M <=> H + CO + M */
    a[164] = 1.87e+17;
    b[164] = -1;
    e[164] = 17000;

    /*reaction 166: HCO + O2 <=> HO2 + CO */
    a[165] = 1.345e+13;
    b[165] = 0;
    e[165] = 400;

    /*reaction 167: CH2OH + O2 <=> HO2 + CH2O */
    a[166] = 1.8e+13;
    b[166] = 0;
    e[166] = 900;

    /*reaction 168: CH3O + O2 <=> HO2 + CH2O */
    a[167] = 4.28e-13;
    b[167] = 7.6;
    e[167] = -3530;

    /*reaction 169: C2H + O2 <=> HCO + CO */
    a[168] = 1e+13;
    b[168] = 0;
    e[168] = -755;

    /*reaction 170: C2H + H2 <=> H + C2H2 */
    a[169] = 5.68e+10;
    b[169] = 0.9;
    e[169] = 1993;

    /*reaction 171: C2H3 + O2 <=> HCO + CH2O */
    a[170] = 4.58e+16;
    b[170] = -1.39;
    e[170] = 1015;

    /*reaction 172: C2H4 (+M) <=> H2 + C2H2 (+M) */
    a[171] = 8e+12;
    b[171] = 0.44;
    e[171] = 86770;

    /*reaction 173: C2H5 + O2 <=> HO2 + C2H4 */
    a[172] = 8.4e+11;
    b[172] = 0;
    e[172] = 3875;

    /*reaction 174: HCCO + O2 <=> OH + 2 CO */
    a[173] = 3.2e+12;
    b[173] = 0;
    e[173] = 854;

    /*reaction 175: 2 HCCO <=> 2 CO + C2H2 */
    a[174] = 1e+13;
    b[174] = 0;
    e[174] = 0;

    /*reaction 176: O + CH3 => H + H2 + CO */
    a[175] = 3.37e+13;
    b[175] = 0;
    e[175] = 0;

    /*reaction 177: O + C2H4 <=> H + CH2CHO */
    a[176] = 6.7e+06;
    b[176] = 1.83;
    e[176] = 220;

    /*reaction 178: O + C2H5 <=> H + CH3CHO */
    a[177] = 1.096e+14;
    b[177] = 0;
    e[177] = 0;

    /*reaction 179: OH + HO2 <=> O2 + H2O */
    a[178] = 5e+15;
    b[178] = 0;
    e[178] = 17330;

    /*reaction 180: OH + CH3 => H2 + CH2O */
    a[179] = 8e+09;
    b[179] = 0.5;
    e[179] = -1755;

    /*reaction 181: CH + H2 (+M) <=> CH3 (+M) */
    a[180] = 1.97e+12;
    b[180] = 0.43;
    e[180] = -370;

    /*reaction 182: CH2 + O2 => 2 H + CO2 */
    a[181] = 5.8e+12;
    b[181] = 0;
    e[181] = 1500;

    /*reaction 183: CH2 + O2 <=> O + CH2O */
    a[182] = 2.4e+12;
    b[182] = 0;
    e[182] = 1500;

    /*reaction 184: CH2 + CH2 => 2 H + C2H2 */
    a[183] = 2e+14;
    b[183] = 0;
    e[183] = 10989;

    /*reaction 185: CH2(S) + H2O => H2 + CH2O */
    a[184] = 6.82e+10;
    b[184] = 0.25;
    e[184] = -935;

    /*reaction 186: C2H3 + O2 <=> O + CH2CHO */
    a[185] = 3.03e+11;
    b[185] = 0.29;
    e[185] = 11;

    /*reaction 187: C2H3 + O2 <=> HO2 + C2H2 */
    a[186] = 1.337e+06;
    b[186] = 1.61;
    e[186] = -384;

    /*reaction 188: O + CH3CHO <=> OH + CH2CHO */
    a[187] = 5.84e+12;
    b[187] = 0;
    e[187] = 1808;

    /*reaction 189: O + CH3CHO => OH + CH3 + CO */
    a[188] = 5.84e+12;
    b[188] = 0;
    e[188] = 1808;

    /*reaction 190: O2 + CH3CHO => HO2 + CH3 + CO */
    a[189] = 3.01e+13;
    b[189] = 0;
    e[189] = 39150;

    /*reaction 191: H + CH3CHO <=> CH2CHO + H2 */
    a[190] = 2.05e+09;
    b[190] = 1.16;
    e[190] = 2405;

    /*reaction 192: H + CH3CHO => CH3 + H2 + CO */
    a[191] = 2.05e+09;
    b[191] = 1.16;
    e[191] = 2405;

    /*reaction 193: OH + CH3CHO => CH3 + H2O + CO */
    a[192] = 2.343e+10;
    b[192] = 0.73;
    e[192] = -1113;

    /*reaction 194: HO2 + CH3CHO => CH3 + H2O2 + CO */
    a[193] = 3.01e+12;
    b[193] = 0;
    e[193] = 11923;

    /*reaction 195: CH3 + CH3CHO => CH3 + CH4 + CO */
    a[194] = 2.72e+06;
    b[194] = 1.77;
    e[194] = 5920;

    /*reaction 196: H + CH2CO (+M) <=> CH2CHO (+M) */
    a[195] = 4.865e+11;
    b[195] = 0.422;
    e[195] = -1755;

    /*reaction 197: O + CH2CHO => H + CH2 + CO2 */
    a[196] = 1.5e+14;
    b[196] = 0;
    e[196] = 0;

    /*reaction 198: O2 + CH2CHO => OH + CO + CH2O */
    a[197] = 1.81e+10;
    b[197] = 0;
    e[197] = 0;

    /*reaction 199: O2 + CH2CHO => OH + 2 HCO */
    a[198] = 2.35e+10;
    b[198] = 0;
    e[198] = 0;

    /*reaction 200: H + CH2CHO <=> CH3 + HCO */
    a[199] = 2.2e+13;
    b[199] = 0;
    e[199] = 0;

    /*reaction 201: H + CH2CHO <=> CH2CO + H2 */
    a[200] = 1.1e+13;
    b[200] = 0;
    e[200] = 0;

    /*reaction 202: OH + CH2CHO <=> H2O + CH2CO */
    a[201] = 1.2e+13;
    b[201] = 0;
    e[201] = 0;

    /*reaction 203: OH + CH2CHO <=> HCO + CH2OH */
    a[202] = 3.01e+13;
    b[202] = 0;
    e[202] = 0;

    /*reaction 204: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    a[203] = 9.43e+12;
    b[203] = 0;
    e[203] = 0;

    /*reaction 205: O + C3H8 <=> OH + C3H7 */
    a[204] = 193000;
    b[204] = 2.68;
    e[204] = 3716;

    /*reaction 206: H + C3H8 <=> C3H7 + H2 */
    a[205] = 1.32e+06;
    b[205] = 2.54;
    e[205] = 6756;

    /*reaction 207: OH + C3H8 <=> C3H7 + H2O */
    a[206] = 3.16e+07;
    b[206] = 1.8;
    e[206] = 934;

    /*reaction 208: C3H7 + H2O2 <=> HO2 + C3H8 */
    a[207] = 378;
    b[207] = 2.72;
    e[207] = 1500;

    /*reaction 209: CH3 + C3H8 <=> C3H7 + CH4 */
    a[208] = 0.903;
    b[208] = 3.65;
    e[208] = 7154;

    /*reaction 210: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    a[209] = 2.55e+06;
    b[209] = 1.6;
    e[209] = 5700;

    /*reaction 211: O + C3H7 <=> C2H5 + CH2O */
    a[210] = 9.64e+13;
    b[210] = 0;
    e[210] = 0;

    /*reaction 212: H + C3H7 (+M) <=> C3H8 (+M) */
    a[211] = 3.613e+13;
    b[211] = 0;
    e[211] = 0;

    /*reaction 213: H + C3H7 <=> CH3 + C2H5 */
    a[212] = 4.06e+06;
    b[212] = 2.19;
    e[212] = 890;

    /*reaction 214: OH + C3H7 <=> C2H5 + CH2OH */
    a[213] = 2.41e+13;
    b[213] = 0;
    e[213] = 0;

    /*reaction 215: HO2 + C3H7 <=> O2 + C3H8 */
    a[214] = 2.55e+10;
    b[214] = 0.255;
    e[214] = -943;

    /*reaction 216: HO2 + C3H7 => OH + C2H5 + CH2O */
    a[215] = 2.41e+13;
    b[215] = 0;
    e[215] = 0;

    /*reaction 217: CH3 + C3H7 <=> 2 C2H5 */
    a[216] = 1.927e+13;
    b[216] = -0.32;
    e[216] = 0;

    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double * T, double * C, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[35]; /* temporary storage */

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

    /*reaction 37: H + O2 <=> O + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: 2 H + M <=> H2 + M */
    eqcon[37] *= 1e+06; 

    /*reaction 39: 2 H + H2 <=> 2 H2 */
    eqcon[38] *= 1e+06; 

    /*reaction 40: 2 H + H2O <=> H2 + H2O */
    eqcon[39] *= 1e+06; 

    /*reaction 41: 2 H + CO2 <=> H2 + CO2 */
    eqcon[40] *= 1e+06; 

    /*reaction 42: H + OH + M <=> H2O + M */
    eqcon[41] *= 1e+06; 

    /*reaction 43: H + HO2 <=> O + H2O */
    /*eqcon[42] *= 1;  */

    /*reaction 44: H + HO2 <=> O2 + H2 */
    /*eqcon[43] *= 1;  */

    /*reaction 45: H + HO2 <=> 2 OH */
    /*eqcon[44] *= 1;  */

    /*reaction 46: H + H2O2 <=> HO2 + H2 */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + H2O2 <=> OH + H2O */
    /*eqcon[46] *= 1;  */

    /*reaction 48: H + CH <=> C + H2 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[48] *= 1e+06; 

    /*reaction 50: H + CH2(S) <=> CH + H2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[50] *= 1e+06; 

    /*reaction 52: H + CH4 <=> CH3 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: H + HCO (+M) <=> CH2O (+M) */
    eqcon[52] *= 1e+06; 

    /*reaction 54: H + HCO <=> H2 + CO */
    /*eqcon[53] *= 1;  */

    /*reaction 55: H + CH2O (+M) <=> CH2OH (+M) */
    eqcon[54] *= 1e+06; 

    /*reaction 56: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[55] *= 1e+06; 

    /*reaction 57: H + CH2O <=> HCO + H2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: H + CH2OH (+M) <=> CH3OH (+M) */
    eqcon[57] *= 1e+06; 

    /*reaction 59: H + CH2OH <=> H2 + CH2O */
    /*eqcon[58] *= 1;  */

    /*reaction 60: H + CH2OH <=> OH + CH3 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: H + CH2OH <=> CH2(S) + H2O */
    /*eqcon[60] *= 1;  */

    /*reaction 62: H + CH3O (+M) <=> CH3OH (+M) */
    eqcon[61] *= 1e+06; 

    /*reaction 63: H + CH3O <=> H + CH2OH */
    /*eqcon[62] *= 1;  */

    /*reaction 64: H + CH3O <=> H2 + CH2O */
    /*eqcon[63] *= 1;  */

    /*reaction 65: H + CH3O <=> OH + CH3 */
    /*eqcon[64] *= 1;  */

    /*reaction 66: H + CH3O <=> CH2(S) + H2O */
    /*eqcon[65] *= 1;  */

    /*reaction 67: H + CH3OH <=> CH2OH + H2 */
    /*eqcon[66] *= 1;  */

    /*reaction 68: H + CH3OH <=> CH3O + H2 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: H + C2H (+M) <=> C2H2 (+M) */
    eqcon[68] *= 1e+06; 

    /*reaction 70: H + C2H2 (+M) <=> C2H3 (+M) */
    eqcon[69] *= 1e+06; 

    /*reaction 71: H + C2H3 (+M) <=> C2H4 (+M) */
    eqcon[70] *= 1e+06; 

    /*reaction 72: H + C2H3 <=> H2 + C2H2 */
    /*eqcon[71] *= 1;  */

    /*reaction 73: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[72] *= 1e+06; 

    /*reaction 74: H + C2H4 <=> C2H3 + H2 */
    /*eqcon[73] *= 1;  */

    /*reaction 75: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[74] *= 1e+06; 

    /*reaction 76: H + C2H5 <=> H2 + C2H4 */
    /*eqcon[75] *= 1;  */

    /*reaction 77: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: H + HCCO <=> CH2(S) + CO */
    /*eqcon[77] *= 1;  */

    /*reaction 79: H + CH2CO <=> HCCO + H2 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: H + CH2CO <=> CH3 + CO */
    /*eqcon[79] *= 1;  */

    /*reaction 81: H + HCCOH <=> H + CH2CO */
    /*eqcon[80] *= 1;  */

    /*reaction 82: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[81] *= 1e+06; 

    /*reaction 83: OH + H2 <=> H + H2O */
    /*eqcon[82] *= 1;  */

    /*reaction 84: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[83] *= 1e+06; 

    /*reaction 85: 2 OH <=> O + H2O */
    /*eqcon[84] *= 1;  */

    /*reaction 86: OH + HO2 <=> O2 + H2O */
    /*eqcon[85] *= 1;  */

    /*reaction 87: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[86] *= 1;  */

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[87] *= 1;  */

    /*reaction 89: OH + C <=> H + CO */
    /*eqcon[88] *= 1;  */

    /*reaction 90: OH + CH <=> H + HCO */
    /*eqcon[89] *= 1;  */

    /*reaction 91: OH + CH2 <=> H + CH2O */
    /*eqcon[90] *= 1;  */

    /*reaction 92: OH + CH2 <=> CH + H2O */
    /*eqcon[91] *= 1;  */

    /*reaction 93: OH + CH2(S) <=> H + CH2O */
    /*eqcon[92] *= 1;  */

    /*reaction 94: OH + CH3 (+M) <=> CH3OH (+M) */
    eqcon[93] *= 1e+06; 

    /*reaction 95: OH + CH3 <=> CH2 + H2O */
    /*eqcon[94] *= 1;  */

    /*reaction 96: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[95] *= 1;  */

    /*reaction 97: OH + CH4 <=> CH3 + H2O */
    /*eqcon[96] *= 1;  */

    /*reaction 98: OH + CO <=> H + CO2 */
    /*eqcon[97] *= 1;  */

    /*reaction 99: OH + HCO <=> H2O + CO */
    /*eqcon[98] *= 1;  */

    /*reaction 100: OH + CH2O <=> HCO + H2O */
    /*eqcon[99] *= 1;  */

    /*reaction 101: OH + CH2OH <=> H2O + CH2O */
    /*eqcon[100] *= 1;  */

    /*reaction 102: OH + CH3O <=> H2O + CH2O */
    /*eqcon[101] *= 1;  */

    /*reaction 103: OH + CH3OH <=> CH2OH + H2O */
    /*eqcon[102] *= 1;  */

    /*reaction 104: OH + CH3OH <=> CH3O + H2O */
    /*eqcon[103] *= 1;  */

    /*reaction 105: OH + C2H <=> H + HCCO */
    /*eqcon[104] *= 1;  */

    /*reaction 106: OH + C2H2 <=> H + CH2CO */
    /*eqcon[105] *= 1;  */

    /*reaction 107: OH + C2H2 <=> H + HCCOH */
    /*eqcon[106] *= 1;  */

    /*reaction 108: OH + C2H2 <=> C2H + H2O */
    /*eqcon[107] *= 1;  */

    /*reaction 109: OH + C2H2 <=> CH3 + CO */
    /*eqcon[108] *= 1;  */

    /*reaction 110: OH + C2H3 <=> H2O + C2H2 */
    /*eqcon[109] *= 1;  */

    /*reaction 111: OH + C2H4 <=> C2H3 + H2O */
    /*eqcon[110] *= 1;  */

    /*reaction 112: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[111] *= 1;  */

    /*reaction 113: OH + CH2CO <=> HCCO + H2O */
    /*eqcon[112] *= 1;  */

    /*reaction 114: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[113] *= 1;  */

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[114] *= 1;  */

    /*reaction 116: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[115] *= 1;  */

    /*reaction 117: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[116] *= 1;  */

    /*reaction 118: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[117] *= 1;  */

    /*reaction 119: HO2 + CO <=> OH + CO2 */
    /*eqcon[118] *= 1;  */

    /*reaction 120: HO2 + CH2O <=> HCO + H2O2 */
    /*eqcon[119] *= 1;  */

    /*reaction 121: C + O2 <=> O + CO */
    /*eqcon[120] *= 1;  */

    /*reaction 122: C + CH2 <=> H + C2H */
    /*eqcon[121] *= 1;  */

    /*reaction 123: C + CH3 <=> H + C2H2 */
    /*eqcon[122] *= 1;  */

    /*reaction 124: CH + O2 <=> O + HCO */
    /*eqcon[123] *= 1;  */

    /*reaction 125: CH + H2 <=> H + CH2 */
    /*eqcon[124] *= 1;  */

    /*reaction 126: CH + H2O <=> H + CH2O */
    /*eqcon[125] *= 1;  */

    /*reaction 127: CH + CH2 <=> H + C2H2 */
    /*eqcon[126] *= 1;  */

    /*reaction 128: CH + CH3 <=> H + C2H3 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: CH + CH4 <=> H + C2H4 */
    /*eqcon[128] *= 1;  */

    /*reaction 130: CH + CO (+M) <=> HCCO (+M) */
    eqcon[129] *= 1e+06; 

    /*reaction 131: CH + CO2 <=> HCO + CO */
    /*eqcon[130] *= 1;  */

    /*reaction 132: CH + CH2O <=> H + CH2CO */
    /*eqcon[131] *= 1;  */

    /*reaction 133: CH + HCCO <=> CO + C2H2 */
    /*eqcon[132] *= 1;  */

    /*reaction 134: CH2 + O2 => OH + H + CO */
    eqcon[133] *= 1e-06; 

    /*reaction 135: CH2 + H2 <=> H + CH3 */
    /*eqcon[134] *= 1;  */

    /*reaction 136: 2 CH2 <=> H2 + C2H2 */
    /*eqcon[135] *= 1;  */

    /*reaction 137: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[137] *= 1;  */

    /*reaction 139: CH2 + CO (+M) <=> CH2CO (+M) */
    eqcon[138] *= 1e+06; 

    /*reaction 140: CH2 + HCCO <=> C2H3 + CO */
    /*eqcon[139] *= 1;  */

    /*reaction 141: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[140] *= 1;  */

    /*reaction 142: CH2(S) + O2 <=> H + OH + CO */
    eqcon[141] *= 1e-06; 

    /*reaction 143: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[142] *= 1;  */

    /*reaction 144: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[143] *= 1;  */

    /*reaction 145: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    eqcon[144] *= 1e+06; 

    /*reaction 146: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[145] *= 1;  */

    /*reaction 147: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[146] *= 1;  */

    /*reaction 148: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[147] *= 1;  */

    /*reaction 149: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[148] *= 1;  */

    /*reaction 150: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[149] *= 1;  */

    /*reaction 151: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[150] *= 1;  */

    /*reaction 152: CH2(S) + C2H6 <=> CH3 + C2H5 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: CH3 + O2 <=> O + CH3O */
    /*eqcon[152] *= 1;  */

    /*reaction 154: CH3 + O2 <=> OH + CH2O */
    /*eqcon[153] *= 1;  */

    /*reaction 155: CH3 + H2O2 <=> HO2 + CH4 */
    /*eqcon[154] *= 1;  */

    /*reaction 156: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[155] *= 1e+06; 

    /*reaction 157: 2 CH3 <=> H + C2H5 */
    /*eqcon[156] *= 1;  */

    /*reaction 158: CH3 + HCO <=> CH4 + CO */
    /*eqcon[157] *= 1;  */

    /*reaction 159: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[158] *= 1;  */

    /*reaction 160: CH3 + CH3OH <=> CH2OH + CH4 */
    /*eqcon[159] *= 1;  */

    /*reaction 161: CH3 + CH3OH <=> CH3O + CH4 */
    /*eqcon[160] *= 1;  */

    /*reaction 162: CH3 + C2H4 <=> C2H3 + CH4 */
    /*eqcon[161] *= 1;  */

    /*reaction 163: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[162] *= 1;  */

    /*reaction 164: HCO + H2O <=> H + CO + H2O */
    eqcon[163] *= 1e-06; 

    /*reaction 165: HCO + M <=> H + CO + M */
    eqcon[164] *= 1e-06; 

    /*reaction 166: HCO + O2 <=> HO2 + CO */
    /*eqcon[165] *= 1;  */

    /*reaction 167: CH2OH + O2 <=> HO2 + CH2O */
    /*eqcon[166] *= 1;  */

    /*reaction 168: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[167] *= 1;  */

    /*reaction 169: C2H + O2 <=> HCO + CO */
    /*eqcon[168] *= 1;  */

    /*reaction 170: C2H + H2 <=> H + C2H2 */
    /*eqcon[169] *= 1;  */

    /*reaction 171: C2H3 + O2 <=> HCO + CH2O */
    /*eqcon[170] *= 1;  */

    /*reaction 172: C2H4 (+M) <=> H2 + C2H2 (+M) */
    eqcon[171] *= 1e-06; 

    /*reaction 173: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[172] *= 1;  */

    /*reaction 174: HCCO + O2 <=> OH + 2 CO */
    eqcon[173] *= 1e-06; 

    /*reaction 175: 2 HCCO <=> 2 CO + C2H2 */
    eqcon[174] *= 1e-06; 

    /*reaction 176: O + CH3 => H + H2 + CO */
    eqcon[175] *= 1e-06; 

    /*reaction 177: O + C2H4 <=> H + CH2CHO */
    /*eqcon[176] *= 1;  */

    /*reaction 178: O + C2H5 <=> H + CH3CHO */
    /*eqcon[177] *= 1;  */

    /*reaction 179: OH + HO2 <=> O2 + H2O */
    /*eqcon[178] *= 1;  */

    /*reaction 180: OH + CH3 => H2 + CH2O */
    /*eqcon[179] *= 1;  */

    /*reaction 181: CH + H2 (+M) <=> CH3 (+M) */
    eqcon[180] *= 1e+06; 

    /*reaction 182: CH2 + O2 => 2 H + CO2 */
    eqcon[181] *= 1e-06; 

    /*reaction 183: CH2 + O2 <=> O + CH2O */
    /*eqcon[182] *= 1;  */

    /*reaction 184: CH2 + CH2 => 2 H + C2H2 */
    eqcon[183] *= 1e-06; 

    /*reaction 185: CH2(S) + H2O => H2 + CH2O */
    /*eqcon[184] *= 1;  */

    /*reaction 186: C2H3 + O2 <=> O + CH2CHO */
    /*eqcon[185] *= 1;  */

    /*reaction 187: C2H3 + O2 <=> HO2 + C2H2 */
    /*eqcon[186] *= 1;  */

    /*reaction 188: O + CH3CHO <=> OH + CH2CHO */
    /*eqcon[187] *= 1;  */

    /*reaction 189: O + CH3CHO => OH + CH3 + CO */
    eqcon[188] *= 1e-06; 

    /*reaction 190: O2 + CH3CHO => HO2 + CH3 + CO */
    eqcon[189] *= 1e-06; 

    /*reaction 191: H + CH3CHO <=> CH2CHO + H2 */
    /*eqcon[190] *= 1;  */

    /*reaction 192: H + CH3CHO => CH3 + H2 + CO */
    eqcon[191] *= 1e-06; 

    /*reaction 193: OH + CH3CHO => CH3 + H2O + CO */
    eqcon[192] *= 1e-06; 

    /*reaction 194: HO2 + CH3CHO => CH3 + H2O2 + CO */
    eqcon[193] *= 1e-06; 

    /*reaction 195: CH3 + CH3CHO => CH3 + CH4 + CO */
    eqcon[194] *= 1e-06; 

    /*reaction 196: H + CH2CO (+M) <=> CH2CHO (+M) */
    eqcon[195] *= 1e+06; 

    /*reaction 197: O + CH2CHO => H + CH2 + CO2 */
    eqcon[196] *= 1e-06; 

    /*reaction 198: O2 + CH2CHO => OH + CO + CH2O */
    eqcon[197] *= 1e-06; 

    /*reaction 199: O2 + CH2CHO => OH + 2 HCO */
    eqcon[198] *= 1e-06; 

    /*reaction 200: H + CH2CHO <=> CH3 + HCO */
    /*eqcon[199] *= 1;  */

    /*reaction 201: H + CH2CHO <=> CH2CO + H2 */
    /*eqcon[200] *= 1;  */

    /*reaction 202: OH + CH2CHO <=> H2O + CH2CO */
    /*eqcon[201] *= 1;  */

    /*reaction 203: OH + CH2CHO <=> HCO + CH2OH */
    /*eqcon[202] *= 1;  */

    /*reaction 204: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    eqcon[203] *= 1e+06; 

    /*reaction 205: O + C3H8 <=> OH + C3H7 */
    /*eqcon[204] *= 1;  */

    /*reaction 206: H + C3H8 <=> C3H7 + H2 */
    /*eqcon[205] *= 1;  */

    /*reaction 207: OH + C3H8 <=> C3H7 + H2O */
    /*eqcon[206] *= 1;  */

    /*reaction 208: C3H7 + H2O2 <=> HO2 + C3H8 */
    /*eqcon[207] *= 1;  */

    /*reaction 209: CH3 + C3H8 <=> C3H7 + CH4 */
    /*eqcon[208] *= 1;  */

    /*reaction 210: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    eqcon[209] *= 1e+06; 

    /*reaction 211: O + C3H7 <=> C2H5 + CH2O */
    /*eqcon[210] *= 1;  */

    /*reaction 212: H + C3H7 (+M) <=> C3H8 (+M) */
    eqcon[211] *= 1e+06; 

    /*reaction 213: H + C3H7 <=> CH3 + C2H5 */
    /*eqcon[212] *= 1;  */

    /*reaction 214: OH + C3H7 <=> C2H5 + CH2OH */
    /*eqcon[213] *= 1;  */

    /*reaction 215: HO2 + C3H7 <=> O2 + C3H8 */
    /*eqcon[214] *= 1;  */

    /*reaction 216: HO2 + C3H7 => OH + C2H5 + CH2O */
    eqcon[215] *= 1e-06; 

    /*reaction 217: CH3 + C3H7 <=> 2 C2H5 */
    /*eqcon[216] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[35]; /* temporary storage */

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

    /*reaction 37: H + O2 <=> O + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: 2 H + M <=> H2 + M */
    eqcon[37] *= 1e+06; 

    /*reaction 39: 2 H + H2 <=> 2 H2 */
    eqcon[38] *= 1e+06; 

    /*reaction 40: 2 H + H2O <=> H2 + H2O */
    eqcon[39] *= 1e+06; 

    /*reaction 41: 2 H + CO2 <=> H2 + CO2 */
    eqcon[40] *= 1e+06; 

    /*reaction 42: H + OH + M <=> H2O + M */
    eqcon[41] *= 1e+06; 

    /*reaction 43: H + HO2 <=> O + H2O */
    /*eqcon[42] *= 1;  */

    /*reaction 44: H + HO2 <=> O2 + H2 */
    /*eqcon[43] *= 1;  */

    /*reaction 45: H + HO2 <=> 2 OH */
    /*eqcon[44] *= 1;  */

    /*reaction 46: H + H2O2 <=> HO2 + H2 */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + H2O2 <=> OH + H2O */
    /*eqcon[46] *= 1;  */

    /*reaction 48: H + CH <=> C + H2 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[48] *= 1e+06; 

    /*reaction 50: H + CH2(S) <=> CH + H2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[50] *= 1e+06; 

    /*reaction 52: H + CH4 <=> CH3 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: H + HCO (+M) <=> CH2O (+M) */
    eqcon[52] *= 1e+06; 

    /*reaction 54: H + HCO <=> H2 + CO */
    /*eqcon[53] *= 1;  */

    /*reaction 55: H + CH2O (+M) <=> CH2OH (+M) */
    eqcon[54] *= 1e+06; 

    /*reaction 56: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[55] *= 1e+06; 

    /*reaction 57: H + CH2O <=> HCO + H2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: H + CH2OH (+M) <=> CH3OH (+M) */
    eqcon[57] *= 1e+06; 

    /*reaction 59: H + CH2OH <=> H2 + CH2O */
    /*eqcon[58] *= 1;  */

    /*reaction 60: H + CH2OH <=> OH + CH3 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: H + CH2OH <=> CH2(S) + H2O */
    /*eqcon[60] *= 1;  */

    /*reaction 62: H + CH3O (+M) <=> CH3OH (+M) */
    eqcon[61] *= 1e+06; 

    /*reaction 63: H + CH3O <=> H + CH2OH */
    /*eqcon[62] *= 1;  */

    /*reaction 64: H + CH3O <=> H2 + CH2O */
    /*eqcon[63] *= 1;  */

    /*reaction 65: H + CH3O <=> OH + CH3 */
    /*eqcon[64] *= 1;  */

    /*reaction 66: H + CH3O <=> CH2(S) + H2O */
    /*eqcon[65] *= 1;  */

    /*reaction 67: H + CH3OH <=> CH2OH + H2 */
    /*eqcon[66] *= 1;  */

    /*reaction 68: H + CH3OH <=> CH3O + H2 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: H + C2H (+M) <=> C2H2 (+M) */
    eqcon[68] *= 1e+06; 

    /*reaction 70: H + C2H2 (+M) <=> C2H3 (+M) */
    eqcon[69] *= 1e+06; 

    /*reaction 71: H + C2H3 (+M) <=> C2H4 (+M) */
    eqcon[70] *= 1e+06; 

    /*reaction 72: H + C2H3 <=> H2 + C2H2 */
    /*eqcon[71] *= 1;  */

    /*reaction 73: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[72] *= 1e+06; 

    /*reaction 74: H + C2H4 <=> C2H3 + H2 */
    /*eqcon[73] *= 1;  */

    /*reaction 75: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[74] *= 1e+06; 

    /*reaction 76: H + C2H5 <=> H2 + C2H4 */
    /*eqcon[75] *= 1;  */

    /*reaction 77: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: H + HCCO <=> CH2(S) + CO */
    /*eqcon[77] *= 1;  */

    /*reaction 79: H + CH2CO <=> HCCO + H2 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: H + CH2CO <=> CH3 + CO */
    /*eqcon[79] *= 1;  */

    /*reaction 81: H + HCCOH <=> H + CH2CO */
    /*eqcon[80] *= 1;  */

    /*reaction 82: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[81] *= 1e+06; 

    /*reaction 83: OH + H2 <=> H + H2O */
    /*eqcon[82] *= 1;  */

    /*reaction 84: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[83] *= 1e+06; 

    /*reaction 85: 2 OH <=> O + H2O */
    /*eqcon[84] *= 1;  */

    /*reaction 86: OH + HO2 <=> O2 + H2O */
    /*eqcon[85] *= 1;  */

    /*reaction 87: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[86] *= 1;  */

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[87] *= 1;  */

    /*reaction 89: OH + C <=> H + CO */
    /*eqcon[88] *= 1;  */

    /*reaction 90: OH + CH <=> H + HCO */
    /*eqcon[89] *= 1;  */

    /*reaction 91: OH + CH2 <=> H + CH2O */
    /*eqcon[90] *= 1;  */

    /*reaction 92: OH + CH2 <=> CH + H2O */
    /*eqcon[91] *= 1;  */

    /*reaction 93: OH + CH2(S) <=> H + CH2O */
    /*eqcon[92] *= 1;  */

    /*reaction 94: OH + CH3 (+M) <=> CH3OH (+M) */
    eqcon[93] *= 1e+06; 

    /*reaction 95: OH + CH3 <=> CH2 + H2O */
    /*eqcon[94] *= 1;  */

    /*reaction 96: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[95] *= 1;  */

    /*reaction 97: OH + CH4 <=> CH3 + H2O */
    /*eqcon[96] *= 1;  */

    /*reaction 98: OH + CO <=> H + CO2 */
    /*eqcon[97] *= 1;  */

    /*reaction 99: OH + HCO <=> H2O + CO */
    /*eqcon[98] *= 1;  */

    /*reaction 100: OH + CH2O <=> HCO + H2O */
    /*eqcon[99] *= 1;  */

    /*reaction 101: OH + CH2OH <=> H2O + CH2O */
    /*eqcon[100] *= 1;  */

    /*reaction 102: OH + CH3O <=> H2O + CH2O */
    /*eqcon[101] *= 1;  */

    /*reaction 103: OH + CH3OH <=> CH2OH + H2O */
    /*eqcon[102] *= 1;  */

    /*reaction 104: OH + CH3OH <=> CH3O + H2O */
    /*eqcon[103] *= 1;  */

    /*reaction 105: OH + C2H <=> H + HCCO */
    /*eqcon[104] *= 1;  */

    /*reaction 106: OH + C2H2 <=> H + CH2CO */
    /*eqcon[105] *= 1;  */

    /*reaction 107: OH + C2H2 <=> H + HCCOH */
    /*eqcon[106] *= 1;  */

    /*reaction 108: OH + C2H2 <=> C2H + H2O */
    /*eqcon[107] *= 1;  */

    /*reaction 109: OH + C2H2 <=> CH3 + CO */
    /*eqcon[108] *= 1;  */

    /*reaction 110: OH + C2H3 <=> H2O + C2H2 */
    /*eqcon[109] *= 1;  */

    /*reaction 111: OH + C2H4 <=> C2H3 + H2O */
    /*eqcon[110] *= 1;  */

    /*reaction 112: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[111] *= 1;  */

    /*reaction 113: OH + CH2CO <=> HCCO + H2O */
    /*eqcon[112] *= 1;  */

    /*reaction 114: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[113] *= 1;  */

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[114] *= 1;  */

    /*reaction 116: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[115] *= 1;  */

    /*reaction 117: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[116] *= 1;  */

    /*reaction 118: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[117] *= 1;  */

    /*reaction 119: HO2 + CO <=> OH + CO2 */
    /*eqcon[118] *= 1;  */

    /*reaction 120: HO2 + CH2O <=> HCO + H2O2 */
    /*eqcon[119] *= 1;  */

    /*reaction 121: C + O2 <=> O + CO */
    /*eqcon[120] *= 1;  */

    /*reaction 122: C + CH2 <=> H + C2H */
    /*eqcon[121] *= 1;  */

    /*reaction 123: C + CH3 <=> H + C2H2 */
    /*eqcon[122] *= 1;  */

    /*reaction 124: CH + O2 <=> O + HCO */
    /*eqcon[123] *= 1;  */

    /*reaction 125: CH + H2 <=> H + CH2 */
    /*eqcon[124] *= 1;  */

    /*reaction 126: CH + H2O <=> H + CH2O */
    /*eqcon[125] *= 1;  */

    /*reaction 127: CH + CH2 <=> H + C2H2 */
    /*eqcon[126] *= 1;  */

    /*reaction 128: CH + CH3 <=> H + C2H3 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: CH + CH4 <=> H + C2H4 */
    /*eqcon[128] *= 1;  */

    /*reaction 130: CH + CO (+M) <=> HCCO (+M) */
    eqcon[129] *= 1e+06; 

    /*reaction 131: CH + CO2 <=> HCO + CO */
    /*eqcon[130] *= 1;  */

    /*reaction 132: CH + CH2O <=> H + CH2CO */
    /*eqcon[131] *= 1;  */

    /*reaction 133: CH + HCCO <=> CO + C2H2 */
    /*eqcon[132] *= 1;  */

    /*reaction 134: CH2 + O2 => OH + H + CO */
    eqcon[133] *= 1e-06; 

    /*reaction 135: CH2 + H2 <=> H + CH3 */
    /*eqcon[134] *= 1;  */

    /*reaction 136: 2 CH2 <=> H2 + C2H2 */
    /*eqcon[135] *= 1;  */

    /*reaction 137: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[137] *= 1;  */

    /*reaction 139: CH2 + CO (+M) <=> CH2CO (+M) */
    eqcon[138] *= 1e+06; 

    /*reaction 140: CH2 + HCCO <=> C2H3 + CO */
    /*eqcon[139] *= 1;  */

    /*reaction 141: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[140] *= 1;  */

    /*reaction 142: CH2(S) + O2 <=> H + OH + CO */
    eqcon[141] *= 1e-06; 

    /*reaction 143: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[142] *= 1;  */

    /*reaction 144: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[143] *= 1;  */

    /*reaction 145: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    eqcon[144] *= 1e+06; 

    /*reaction 146: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[145] *= 1;  */

    /*reaction 147: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[146] *= 1;  */

    /*reaction 148: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[147] *= 1;  */

    /*reaction 149: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[148] *= 1;  */

    /*reaction 150: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[149] *= 1;  */

    /*reaction 151: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[150] *= 1;  */

    /*reaction 152: CH2(S) + C2H6 <=> CH3 + C2H5 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: CH3 + O2 <=> O + CH3O */
    /*eqcon[152] *= 1;  */

    /*reaction 154: CH3 + O2 <=> OH + CH2O */
    /*eqcon[153] *= 1;  */

    /*reaction 155: CH3 + H2O2 <=> HO2 + CH4 */
    /*eqcon[154] *= 1;  */

    /*reaction 156: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[155] *= 1e+06; 

    /*reaction 157: 2 CH3 <=> H + C2H5 */
    /*eqcon[156] *= 1;  */

    /*reaction 158: CH3 + HCO <=> CH4 + CO */
    /*eqcon[157] *= 1;  */

    /*reaction 159: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[158] *= 1;  */

    /*reaction 160: CH3 + CH3OH <=> CH2OH + CH4 */
    /*eqcon[159] *= 1;  */

    /*reaction 161: CH3 + CH3OH <=> CH3O + CH4 */
    /*eqcon[160] *= 1;  */

    /*reaction 162: CH3 + C2H4 <=> C2H3 + CH4 */
    /*eqcon[161] *= 1;  */

    /*reaction 163: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[162] *= 1;  */

    /*reaction 164: HCO + H2O <=> H + CO + H2O */
    eqcon[163] *= 1e-06; 

    /*reaction 165: HCO + M <=> H + CO + M */
    eqcon[164] *= 1e-06; 

    /*reaction 166: HCO + O2 <=> HO2 + CO */
    /*eqcon[165] *= 1;  */

    /*reaction 167: CH2OH + O2 <=> HO2 + CH2O */
    /*eqcon[166] *= 1;  */

    /*reaction 168: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[167] *= 1;  */

    /*reaction 169: C2H + O2 <=> HCO + CO */
    /*eqcon[168] *= 1;  */

    /*reaction 170: C2H + H2 <=> H + C2H2 */
    /*eqcon[169] *= 1;  */

    /*reaction 171: C2H3 + O2 <=> HCO + CH2O */
    /*eqcon[170] *= 1;  */

    /*reaction 172: C2H4 (+M) <=> H2 + C2H2 (+M) */
    eqcon[171] *= 1e-06; 

    /*reaction 173: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[172] *= 1;  */

    /*reaction 174: HCCO + O2 <=> OH + 2 CO */
    eqcon[173] *= 1e-06; 

    /*reaction 175: 2 HCCO <=> 2 CO + C2H2 */
    eqcon[174] *= 1e-06; 

    /*reaction 176: O + CH3 => H + H2 + CO */
    eqcon[175] *= 1e-06; 

    /*reaction 177: O + C2H4 <=> H + CH2CHO */
    /*eqcon[176] *= 1;  */

    /*reaction 178: O + C2H5 <=> H + CH3CHO */
    /*eqcon[177] *= 1;  */

    /*reaction 179: OH + HO2 <=> O2 + H2O */
    /*eqcon[178] *= 1;  */

    /*reaction 180: OH + CH3 => H2 + CH2O */
    /*eqcon[179] *= 1;  */

    /*reaction 181: CH + H2 (+M) <=> CH3 (+M) */
    eqcon[180] *= 1e+06; 

    /*reaction 182: CH2 + O2 => 2 H + CO2 */
    eqcon[181] *= 1e-06; 

    /*reaction 183: CH2 + O2 <=> O + CH2O */
    /*eqcon[182] *= 1;  */

    /*reaction 184: CH2 + CH2 => 2 H + C2H2 */
    eqcon[183] *= 1e-06; 

    /*reaction 185: CH2(S) + H2O => H2 + CH2O */
    /*eqcon[184] *= 1;  */

    /*reaction 186: C2H3 + O2 <=> O + CH2CHO */
    /*eqcon[185] *= 1;  */

    /*reaction 187: C2H3 + O2 <=> HO2 + C2H2 */
    /*eqcon[186] *= 1;  */

    /*reaction 188: O + CH3CHO <=> OH + CH2CHO */
    /*eqcon[187] *= 1;  */

    /*reaction 189: O + CH3CHO => OH + CH3 + CO */
    eqcon[188] *= 1e-06; 

    /*reaction 190: O2 + CH3CHO => HO2 + CH3 + CO */
    eqcon[189] *= 1e-06; 

    /*reaction 191: H + CH3CHO <=> CH2CHO + H2 */
    /*eqcon[190] *= 1;  */

    /*reaction 192: H + CH3CHO => CH3 + H2 + CO */
    eqcon[191] *= 1e-06; 

    /*reaction 193: OH + CH3CHO => CH3 + H2O + CO */
    eqcon[192] *= 1e-06; 

    /*reaction 194: HO2 + CH3CHO => CH3 + H2O2 + CO */
    eqcon[193] *= 1e-06; 

    /*reaction 195: CH3 + CH3CHO => CH3 + CH4 + CO */
    eqcon[194] *= 1e-06; 

    /*reaction 196: H + CH2CO (+M) <=> CH2CHO (+M) */
    eqcon[195] *= 1e+06; 

    /*reaction 197: O + CH2CHO => H + CH2 + CO2 */
    eqcon[196] *= 1e-06; 

    /*reaction 198: O2 + CH2CHO => OH + CO + CH2O */
    eqcon[197] *= 1e-06; 

    /*reaction 199: O2 + CH2CHO => OH + 2 HCO */
    eqcon[198] *= 1e-06; 

    /*reaction 200: H + CH2CHO <=> CH3 + HCO */
    /*eqcon[199] *= 1;  */

    /*reaction 201: H + CH2CHO <=> CH2CO + H2 */
    /*eqcon[200] *= 1;  */

    /*reaction 202: OH + CH2CHO <=> H2O + CH2CO */
    /*eqcon[201] *= 1;  */

    /*reaction 203: OH + CH2CHO <=> HCO + CH2OH */
    /*eqcon[202] *= 1;  */

    /*reaction 204: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    eqcon[203] *= 1e+06; 

    /*reaction 205: O + C3H8 <=> OH + C3H7 */
    /*eqcon[204] *= 1;  */

    /*reaction 206: H + C3H8 <=> C3H7 + H2 */
    /*eqcon[205] *= 1;  */

    /*reaction 207: OH + C3H8 <=> C3H7 + H2O */
    /*eqcon[206] *= 1;  */

    /*reaction 208: C3H7 + H2O2 <=> HO2 + C3H8 */
    /*eqcon[207] *= 1;  */

    /*reaction 209: CH3 + C3H8 <=> C3H7 + CH4 */
    /*eqcon[208] *= 1;  */

    /*reaction 210: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    eqcon[209] *= 1e+06; 

    /*reaction 211: O + C3H7 <=> C2H5 + CH2O */
    /*eqcon[210] *= 1;  */

    /*reaction 212: H + C3H7 (+M) <=> C3H8 (+M) */
    eqcon[211] *= 1e+06; 

    /*reaction 213: H + C3H7 <=> CH3 + C2H5 */
    /*eqcon[212] *= 1;  */

    /*reaction 214: OH + C3H7 <=> C2H5 + CH2OH */
    /*eqcon[213] *= 1;  */

    /*reaction 215: HO2 + C3H7 <=> O2 + C3H8 */
    /*eqcon[214] *= 1;  */

    /*reaction 216: HO2 + C3H7 => OH + C2H5 + CH2O */
    eqcon[215] *= 1e-06; 

    /*reaction 217: CH3 + C3H7 <=> 2 C2H5 */
    /*eqcon[216] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[35]; /* temporary storage */

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

    /*reaction 37: H + O2 <=> O + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: 2 H + M <=> H2 + M */
    eqcon[37] *= 1e+06; 

    /*reaction 39: 2 H + H2 <=> 2 H2 */
    eqcon[38] *= 1e+06; 

    /*reaction 40: 2 H + H2O <=> H2 + H2O */
    eqcon[39] *= 1e+06; 

    /*reaction 41: 2 H + CO2 <=> H2 + CO2 */
    eqcon[40] *= 1e+06; 

    /*reaction 42: H + OH + M <=> H2O + M */
    eqcon[41] *= 1e+06; 

    /*reaction 43: H + HO2 <=> O + H2O */
    /*eqcon[42] *= 1;  */

    /*reaction 44: H + HO2 <=> O2 + H2 */
    /*eqcon[43] *= 1;  */

    /*reaction 45: H + HO2 <=> 2 OH */
    /*eqcon[44] *= 1;  */

    /*reaction 46: H + H2O2 <=> HO2 + H2 */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + H2O2 <=> OH + H2O */
    /*eqcon[46] *= 1;  */

    /*reaction 48: H + CH <=> C + H2 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[48] *= 1e+06; 

    /*reaction 50: H + CH2(S) <=> CH + H2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[50] *= 1e+06; 

    /*reaction 52: H + CH4 <=> CH3 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: H + HCO (+M) <=> CH2O (+M) */
    eqcon[52] *= 1e+06; 

    /*reaction 54: H + HCO <=> H2 + CO */
    /*eqcon[53] *= 1;  */

    /*reaction 55: H + CH2O (+M) <=> CH2OH (+M) */
    eqcon[54] *= 1e+06; 

    /*reaction 56: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[55] *= 1e+06; 

    /*reaction 57: H + CH2O <=> HCO + H2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: H + CH2OH (+M) <=> CH3OH (+M) */
    eqcon[57] *= 1e+06; 

    /*reaction 59: H + CH2OH <=> H2 + CH2O */
    /*eqcon[58] *= 1;  */

    /*reaction 60: H + CH2OH <=> OH + CH3 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: H + CH2OH <=> CH2(S) + H2O */
    /*eqcon[60] *= 1;  */

    /*reaction 62: H + CH3O (+M) <=> CH3OH (+M) */
    eqcon[61] *= 1e+06; 

    /*reaction 63: H + CH3O <=> H + CH2OH */
    /*eqcon[62] *= 1;  */

    /*reaction 64: H + CH3O <=> H2 + CH2O */
    /*eqcon[63] *= 1;  */

    /*reaction 65: H + CH3O <=> OH + CH3 */
    /*eqcon[64] *= 1;  */

    /*reaction 66: H + CH3O <=> CH2(S) + H2O */
    /*eqcon[65] *= 1;  */

    /*reaction 67: H + CH3OH <=> CH2OH + H2 */
    /*eqcon[66] *= 1;  */

    /*reaction 68: H + CH3OH <=> CH3O + H2 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: H + C2H (+M) <=> C2H2 (+M) */
    eqcon[68] *= 1e+06; 

    /*reaction 70: H + C2H2 (+M) <=> C2H3 (+M) */
    eqcon[69] *= 1e+06; 

    /*reaction 71: H + C2H3 (+M) <=> C2H4 (+M) */
    eqcon[70] *= 1e+06; 

    /*reaction 72: H + C2H3 <=> H2 + C2H2 */
    /*eqcon[71] *= 1;  */

    /*reaction 73: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[72] *= 1e+06; 

    /*reaction 74: H + C2H4 <=> C2H3 + H2 */
    /*eqcon[73] *= 1;  */

    /*reaction 75: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[74] *= 1e+06; 

    /*reaction 76: H + C2H5 <=> H2 + C2H4 */
    /*eqcon[75] *= 1;  */

    /*reaction 77: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: H + HCCO <=> CH2(S) + CO */
    /*eqcon[77] *= 1;  */

    /*reaction 79: H + CH2CO <=> HCCO + H2 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: H + CH2CO <=> CH3 + CO */
    /*eqcon[79] *= 1;  */

    /*reaction 81: H + HCCOH <=> H + CH2CO */
    /*eqcon[80] *= 1;  */

    /*reaction 82: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[81] *= 1e+06; 

    /*reaction 83: OH + H2 <=> H + H2O */
    /*eqcon[82] *= 1;  */

    /*reaction 84: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[83] *= 1e+06; 

    /*reaction 85: 2 OH <=> O + H2O */
    /*eqcon[84] *= 1;  */

    /*reaction 86: OH + HO2 <=> O2 + H2O */
    /*eqcon[85] *= 1;  */

    /*reaction 87: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[86] *= 1;  */

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[87] *= 1;  */

    /*reaction 89: OH + C <=> H + CO */
    /*eqcon[88] *= 1;  */

    /*reaction 90: OH + CH <=> H + HCO */
    /*eqcon[89] *= 1;  */

    /*reaction 91: OH + CH2 <=> H + CH2O */
    /*eqcon[90] *= 1;  */

    /*reaction 92: OH + CH2 <=> CH + H2O */
    /*eqcon[91] *= 1;  */

    /*reaction 93: OH + CH2(S) <=> H + CH2O */
    /*eqcon[92] *= 1;  */

    /*reaction 94: OH + CH3 (+M) <=> CH3OH (+M) */
    eqcon[93] *= 1e+06; 

    /*reaction 95: OH + CH3 <=> CH2 + H2O */
    /*eqcon[94] *= 1;  */

    /*reaction 96: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[95] *= 1;  */

    /*reaction 97: OH + CH4 <=> CH3 + H2O */
    /*eqcon[96] *= 1;  */

    /*reaction 98: OH + CO <=> H + CO2 */
    /*eqcon[97] *= 1;  */

    /*reaction 99: OH + HCO <=> H2O + CO */
    /*eqcon[98] *= 1;  */

    /*reaction 100: OH + CH2O <=> HCO + H2O */
    /*eqcon[99] *= 1;  */

    /*reaction 101: OH + CH2OH <=> H2O + CH2O */
    /*eqcon[100] *= 1;  */

    /*reaction 102: OH + CH3O <=> H2O + CH2O */
    /*eqcon[101] *= 1;  */

    /*reaction 103: OH + CH3OH <=> CH2OH + H2O */
    /*eqcon[102] *= 1;  */

    /*reaction 104: OH + CH3OH <=> CH3O + H2O */
    /*eqcon[103] *= 1;  */

    /*reaction 105: OH + C2H <=> H + HCCO */
    /*eqcon[104] *= 1;  */

    /*reaction 106: OH + C2H2 <=> H + CH2CO */
    /*eqcon[105] *= 1;  */

    /*reaction 107: OH + C2H2 <=> H + HCCOH */
    /*eqcon[106] *= 1;  */

    /*reaction 108: OH + C2H2 <=> C2H + H2O */
    /*eqcon[107] *= 1;  */

    /*reaction 109: OH + C2H2 <=> CH3 + CO */
    /*eqcon[108] *= 1;  */

    /*reaction 110: OH + C2H3 <=> H2O + C2H2 */
    /*eqcon[109] *= 1;  */

    /*reaction 111: OH + C2H4 <=> C2H3 + H2O */
    /*eqcon[110] *= 1;  */

    /*reaction 112: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[111] *= 1;  */

    /*reaction 113: OH + CH2CO <=> HCCO + H2O */
    /*eqcon[112] *= 1;  */

    /*reaction 114: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[113] *= 1;  */

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[114] *= 1;  */

    /*reaction 116: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[115] *= 1;  */

    /*reaction 117: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[116] *= 1;  */

    /*reaction 118: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[117] *= 1;  */

    /*reaction 119: HO2 + CO <=> OH + CO2 */
    /*eqcon[118] *= 1;  */

    /*reaction 120: HO2 + CH2O <=> HCO + H2O2 */
    /*eqcon[119] *= 1;  */

    /*reaction 121: C + O2 <=> O + CO */
    /*eqcon[120] *= 1;  */

    /*reaction 122: C + CH2 <=> H + C2H */
    /*eqcon[121] *= 1;  */

    /*reaction 123: C + CH3 <=> H + C2H2 */
    /*eqcon[122] *= 1;  */

    /*reaction 124: CH + O2 <=> O + HCO */
    /*eqcon[123] *= 1;  */

    /*reaction 125: CH + H2 <=> H + CH2 */
    /*eqcon[124] *= 1;  */

    /*reaction 126: CH + H2O <=> H + CH2O */
    /*eqcon[125] *= 1;  */

    /*reaction 127: CH + CH2 <=> H + C2H2 */
    /*eqcon[126] *= 1;  */

    /*reaction 128: CH + CH3 <=> H + C2H3 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: CH + CH4 <=> H + C2H4 */
    /*eqcon[128] *= 1;  */

    /*reaction 130: CH + CO (+M) <=> HCCO (+M) */
    eqcon[129] *= 1e+06; 

    /*reaction 131: CH + CO2 <=> HCO + CO */
    /*eqcon[130] *= 1;  */

    /*reaction 132: CH + CH2O <=> H + CH2CO */
    /*eqcon[131] *= 1;  */

    /*reaction 133: CH + HCCO <=> CO + C2H2 */
    /*eqcon[132] *= 1;  */

    /*reaction 134: CH2 + O2 => OH + H + CO */
    eqcon[133] *= 1e-06; 

    /*reaction 135: CH2 + H2 <=> H + CH3 */
    /*eqcon[134] *= 1;  */

    /*reaction 136: 2 CH2 <=> H2 + C2H2 */
    /*eqcon[135] *= 1;  */

    /*reaction 137: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[137] *= 1;  */

    /*reaction 139: CH2 + CO (+M) <=> CH2CO (+M) */
    eqcon[138] *= 1e+06; 

    /*reaction 140: CH2 + HCCO <=> C2H3 + CO */
    /*eqcon[139] *= 1;  */

    /*reaction 141: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[140] *= 1;  */

    /*reaction 142: CH2(S) + O2 <=> H + OH + CO */
    eqcon[141] *= 1e-06; 

    /*reaction 143: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[142] *= 1;  */

    /*reaction 144: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[143] *= 1;  */

    /*reaction 145: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    eqcon[144] *= 1e+06; 

    /*reaction 146: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[145] *= 1;  */

    /*reaction 147: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[146] *= 1;  */

    /*reaction 148: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[147] *= 1;  */

    /*reaction 149: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[148] *= 1;  */

    /*reaction 150: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[149] *= 1;  */

    /*reaction 151: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[150] *= 1;  */

    /*reaction 152: CH2(S) + C2H6 <=> CH3 + C2H5 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: CH3 + O2 <=> O + CH3O */
    /*eqcon[152] *= 1;  */

    /*reaction 154: CH3 + O2 <=> OH + CH2O */
    /*eqcon[153] *= 1;  */

    /*reaction 155: CH3 + H2O2 <=> HO2 + CH4 */
    /*eqcon[154] *= 1;  */

    /*reaction 156: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[155] *= 1e+06; 

    /*reaction 157: 2 CH3 <=> H + C2H5 */
    /*eqcon[156] *= 1;  */

    /*reaction 158: CH3 + HCO <=> CH4 + CO */
    /*eqcon[157] *= 1;  */

    /*reaction 159: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[158] *= 1;  */

    /*reaction 160: CH3 + CH3OH <=> CH2OH + CH4 */
    /*eqcon[159] *= 1;  */

    /*reaction 161: CH3 + CH3OH <=> CH3O + CH4 */
    /*eqcon[160] *= 1;  */

    /*reaction 162: CH3 + C2H4 <=> C2H3 + CH4 */
    /*eqcon[161] *= 1;  */

    /*reaction 163: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[162] *= 1;  */

    /*reaction 164: HCO + H2O <=> H + CO + H2O */
    eqcon[163] *= 1e-06; 

    /*reaction 165: HCO + M <=> H + CO + M */
    eqcon[164] *= 1e-06; 

    /*reaction 166: HCO + O2 <=> HO2 + CO */
    /*eqcon[165] *= 1;  */

    /*reaction 167: CH2OH + O2 <=> HO2 + CH2O */
    /*eqcon[166] *= 1;  */

    /*reaction 168: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[167] *= 1;  */

    /*reaction 169: C2H + O2 <=> HCO + CO */
    /*eqcon[168] *= 1;  */

    /*reaction 170: C2H + H2 <=> H + C2H2 */
    /*eqcon[169] *= 1;  */

    /*reaction 171: C2H3 + O2 <=> HCO + CH2O */
    /*eqcon[170] *= 1;  */

    /*reaction 172: C2H4 (+M) <=> H2 + C2H2 (+M) */
    eqcon[171] *= 1e-06; 

    /*reaction 173: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[172] *= 1;  */

    /*reaction 174: HCCO + O2 <=> OH + 2 CO */
    eqcon[173] *= 1e-06; 

    /*reaction 175: 2 HCCO <=> 2 CO + C2H2 */
    eqcon[174] *= 1e-06; 

    /*reaction 176: O + CH3 => H + H2 + CO */
    eqcon[175] *= 1e-06; 

    /*reaction 177: O + C2H4 <=> H + CH2CHO */
    /*eqcon[176] *= 1;  */

    /*reaction 178: O + C2H5 <=> H + CH3CHO */
    /*eqcon[177] *= 1;  */

    /*reaction 179: OH + HO2 <=> O2 + H2O */
    /*eqcon[178] *= 1;  */

    /*reaction 180: OH + CH3 => H2 + CH2O */
    /*eqcon[179] *= 1;  */

    /*reaction 181: CH + H2 (+M) <=> CH3 (+M) */
    eqcon[180] *= 1e+06; 

    /*reaction 182: CH2 + O2 => 2 H + CO2 */
    eqcon[181] *= 1e-06; 

    /*reaction 183: CH2 + O2 <=> O + CH2O */
    /*eqcon[182] *= 1;  */

    /*reaction 184: CH2 + CH2 => 2 H + C2H2 */
    eqcon[183] *= 1e-06; 

    /*reaction 185: CH2(S) + H2O => H2 + CH2O */
    /*eqcon[184] *= 1;  */

    /*reaction 186: C2H3 + O2 <=> O + CH2CHO */
    /*eqcon[185] *= 1;  */

    /*reaction 187: C2H3 + O2 <=> HO2 + C2H2 */
    /*eqcon[186] *= 1;  */

    /*reaction 188: O + CH3CHO <=> OH + CH2CHO */
    /*eqcon[187] *= 1;  */

    /*reaction 189: O + CH3CHO => OH + CH3 + CO */
    eqcon[188] *= 1e-06; 

    /*reaction 190: O2 + CH3CHO => HO2 + CH3 + CO */
    eqcon[189] *= 1e-06; 

    /*reaction 191: H + CH3CHO <=> CH2CHO + H2 */
    /*eqcon[190] *= 1;  */

    /*reaction 192: H + CH3CHO => CH3 + H2 + CO */
    eqcon[191] *= 1e-06; 

    /*reaction 193: OH + CH3CHO => CH3 + H2O + CO */
    eqcon[192] *= 1e-06; 

    /*reaction 194: HO2 + CH3CHO => CH3 + H2O2 + CO */
    eqcon[193] *= 1e-06; 

    /*reaction 195: CH3 + CH3CHO => CH3 + CH4 + CO */
    eqcon[194] *= 1e-06; 

    /*reaction 196: H + CH2CO (+M) <=> CH2CHO (+M) */
    eqcon[195] *= 1e+06; 

    /*reaction 197: O + CH2CHO => H + CH2 + CO2 */
    eqcon[196] *= 1e-06; 

    /*reaction 198: O2 + CH2CHO => OH + CO + CH2O */
    eqcon[197] *= 1e-06; 

    /*reaction 199: O2 + CH2CHO => OH + 2 HCO */
    eqcon[198] *= 1e-06; 

    /*reaction 200: H + CH2CHO <=> CH3 + HCO */
    /*eqcon[199] *= 1;  */

    /*reaction 201: H + CH2CHO <=> CH2CO + H2 */
    /*eqcon[200] *= 1;  */

    /*reaction 202: OH + CH2CHO <=> H2O + CH2CO */
    /*eqcon[201] *= 1;  */

    /*reaction 203: OH + CH2CHO <=> HCO + CH2OH */
    /*eqcon[202] *= 1;  */

    /*reaction 204: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    eqcon[203] *= 1e+06; 

    /*reaction 205: O + C3H8 <=> OH + C3H7 */
    /*eqcon[204] *= 1;  */

    /*reaction 206: H + C3H8 <=> C3H7 + H2 */
    /*eqcon[205] *= 1;  */

    /*reaction 207: OH + C3H8 <=> C3H7 + H2O */
    /*eqcon[206] *= 1;  */

    /*reaction 208: C3H7 + H2O2 <=> HO2 + C3H8 */
    /*eqcon[207] *= 1;  */

    /*reaction 209: CH3 + C3H8 <=> C3H7 + CH4 */
    /*eqcon[208] *= 1;  */

    /*reaction 210: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    eqcon[209] *= 1e+06; 

    /*reaction 211: O + C3H7 <=> C2H5 + CH2O */
    /*eqcon[210] *= 1;  */

    /*reaction 212: H + C3H7 (+M) <=> C3H8 (+M) */
    eqcon[211] *= 1e+06; 

    /*reaction 213: H + C3H7 <=> CH3 + C2H5 */
    /*eqcon[212] *= 1;  */

    /*reaction 214: OH + C3H7 <=> C2H5 + CH2OH */
    /*eqcon[213] *= 1;  */

    /*reaction 215: HO2 + C3H7 <=> O2 + C3H8 */
    /*eqcon[214] *= 1;  */

    /*reaction 216: HO2 + C3H7 => OH + C2H5 + CH2O */
    eqcon[215] *= 1e-06; 

    /*reaction 217: CH3 + C3H7 <=> 2 C2H5 */
    /*eqcon[216] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[35]; /* temporary storage */

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

    /*reaction 37: H + O2 <=> O + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: 2 H + M <=> H2 + M */
    eqcon[37] *= 1e+06; 

    /*reaction 39: 2 H + H2 <=> 2 H2 */
    eqcon[38] *= 1e+06; 

    /*reaction 40: 2 H + H2O <=> H2 + H2O */
    eqcon[39] *= 1e+06; 

    /*reaction 41: 2 H + CO2 <=> H2 + CO2 */
    eqcon[40] *= 1e+06; 

    /*reaction 42: H + OH + M <=> H2O + M */
    eqcon[41] *= 1e+06; 

    /*reaction 43: H + HO2 <=> O + H2O */
    /*eqcon[42] *= 1;  */

    /*reaction 44: H + HO2 <=> O2 + H2 */
    /*eqcon[43] *= 1;  */

    /*reaction 45: H + HO2 <=> 2 OH */
    /*eqcon[44] *= 1;  */

    /*reaction 46: H + H2O2 <=> HO2 + H2 */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + H2O2 <=> OH + H2O */
    /*eqcon[46] *= 1;  */

    /*reaction 48: H + CH <=> C + H2 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[48] *= 1e+06; 

    /*reaction 50: H + CH2(S) <=> CH + H2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[50] *= 1e+06; 

    /*reaction 52: H + CH4 <=> CH3 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: H + HCO (+M) <=> CH2O (+M) */
    eqcon[52] *= 1e+06; 

    /*reaction 54: H + HCO <=> H2 + CO */
    /*eqcon[53] *= 1;  */

    /*reaction 55: H + CH2O (+M) <=> CH2OH (+M) */
    eqcon[54] *= 1e+06; 

    /*reaction 56: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[55] *= 1e+06; 

    /*reaction 57: H + CH2O <=> HCO + H2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: H + CH2OH (+M) <=> CH3OH (+M) */
    eqcon[57] *= 1e+06; 

    /*reaction 59: H + CH2OH <=> H2 + CH2O */
    /*eqcon[58] *= 1;  */

    /*reaction 60: H + CH2OH <=> OH + CH3 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: H + CH2OH <=> CH2(S) + H2O */
    /*eqcon[60] *= 1;  */

    /*reaction 62: H + CH3O (+M) <=> CH3OH (+M) */
    eqcon[61] *= 1e+06; 

    /*reaction 63: H + CH3O <=> H + CH2OH */
    /*eqcon[62] *= 1;  */

    /*reaction 64: H + CH3O <=> H2 + CH2O */
    /*eqcon[63] *= 1;  */

    /*reaction 65: H + CH3O <=> OH + CH3 */
    /*eqcon[64] *= 1;  */

    /*reaction 66: H + CH3O <=> CH2(S) + H2O */
    /*eqcon[65] *= 1;  */

    /*reaction 67: H + CH3OH <=> CH2OH + H2 */
    /*eqcon[66] *= 1;  */

    /*reaction 68: H + CH3OH <=> CH3O + H2 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: H + C2H (+M) <=> C2H2 (+M) */
    eqcon[68] *= 1e+06; 

    /*reaction 70: H + C2H2 (+M) <=> C2H3 (+M) */
    eqcon[69] *= 1e+06; 

    /*reaction 71: H + C2H3 (+M) <=> C2H4 (+M) */
    eqcon[70] *= 1e+06; 

    /*reaction 72: H + C2H3 <=> H2 + C2H2 */
    /*eqcon[71] *= 1;  */

    /*reaction 73: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[72] *= 1e+06; 

    /*reaction 74: H + C2H4 <=> C2H3 + H2 */
    /*eqcon[73] *= 1;  */

    /*reaction 75: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[74] *= 1e+06; 

    /*reaction 76: H + C2H5 <=> H2 + C2H4 */
    /*eqcon[75] *= 1;  */

    /*reaction 77: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: H + HCCO <=> CH2(S) + CO */
    /*eqcon[77] *= 1;  */

    /*reaction 79: H + CH2CO <=> HCCO + H2 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: H + CH2CO <=> CH3 + CO */
    /*eqcon[79] *= 1;  */

    /*reaction 81: H + HCCOH <=> H + CH2CO */
    /*eqcon[80] *= 1;  */

    /*reaction 82: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[81] *= 1e+06; 

    /*reaction 83: OH + H2 <=> H + H2O */
    /*eqcon[82] *= 1;  */

    /*reaction 84: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[83] *= 1e+06; 

    /*reaction 85: 2 OH <=> O + H2O */
    /*eqcon[84] *= 1;  */

    /*reaction 86: OH + HO2 <=> O2 + H2O */
    /*eqcon[85] *= 1;  */

    /*reaction 87: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[86] *= 1;  */

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[87] *= 1;  */

    /*reaction 89: OH + C <=> H + CO */
    /*eqcon[88] *= 1;  */

    /*reaction 90: OH + CH <=> H + HCO */
    /*eqcon[89] *= 1;  */

    /*reaction 91: OH + CH2 <=> H + CH2O */
    /*eqcon[90] *= 1;  */

    /*reaction 92: OH + CH2 <=> CH + H2O */
    /*eqcon[91] *= 1;  */

    /*reaction 93: OH + CH2(S) <=> H + CH2O */
    /*eqcon[92] *= 1;  */

    /*reaction 94: OH + CH3 (+M) <=> CH3OH (+M) */
    eqcon[93] *= 1e+06; 

    /*reaction 95: OH + CH3 <=> CH2 + H2O */
    /*eqcon[94] *= 1;  */

    /*reaction 96: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[95] *= 1;  */

    /*reaction 97: OH + CH4 <=> CH3 + H2O */
    /*eqcon[96] *= 1;  */

    /*reaction 98: OH + CO <=> H + CO2 */
    /*eqcon[97] *= 1;  */

    /*reaction 99: OH + HCO <=> H2O + CO */
    /*eqcon[98] *= 1;  */

    /*reaction 100: OH + CH2O <=> HCO + H2O */
    /*eqcon[99] *= 1;  */

    /*reaction 101: OH + CH2OH <=> H2O + CH2O */
    /*eqcon[100] *= 1;  */

    /*reaction 102: OH + CH3O <=> H2O + CH2O */
    /*eqcon[101] *= 1;  */

    /*reaction 103: OH + CH3OH <=> CH2OH + H2O */
    /*eqcon[102] *= 1;  */

    /*reaction 104: OH + CH3OH <=> CH3O + H2O */
    /*eqcon[103] *= 1;  */

    /*reaction 105: OH + C2H <=> H + HCCO */
    /*eqcon[104] *= 1;  */

    /*reaction 106: OH + C2H2 <=> H + CH2CO */
    /*eqcon[105] *= 1;  */

    /*reaction 107: OH + C2H2 <=> H + HCCOH */
    /*eqcon[106] *= 1;  */

    /*reaction 108: OH + C2H2 <=> C2H + H2O */
    /*eqcon[107] *= 1;  */

    /*reaction 109: OH + C2H2 <=> CH3 + CO */
    /*eqcon[108] *= 1;  */

    /*reaction 110: OH + C2H3 <=> H2O + C2H2 */
    /*eqcon[109] *= 1;  */

    /*reaction 111: OH + C2H4 <=> C2H3 + H2O */
    /*eqcon[110] *= 1;  */

    /*reaction 112: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[111] *= 1;  */

    /*reaction 113: OH + CH2CO <=> HCCO + H2O */
    /*eqcon[112] *= 1;  */

    /*reaction 114: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[113] *= 1;  */

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[114] *= 1;  */

    /*reaction 116: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[115] *= 1;  */

    /*reaction 117: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[116] *= 1;  */

    /*reaction 118: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[117] *= 1;  */

    /*reaction 119: HO2 + CO <=> OH + CO2 */
    /*eqcon[118] *= 1;  */

    /*reaction 120: HO2 + CH2O <=> HCO + H2O2 */
    /*eqcon[119] *= 1;  */

    /*reaction 121: C + O2 <=> O + CO */
    /*eqcon[120] *= 1;  */

    /*reaction 122: C + CH2 <=> H + C2H */
    /*eqcon[121] *= 1;  */

    /*reaction 123: C + CH3 <=> H + C2H2 */
    /*eqcon[122] *= 1;  */

    /*reaction 124: CH + O2 <=> O + HCO */
    /*eqcon[123] *= 1;  */

    /*reaction 125: CH + H2 <=> H + CH2 */
    /*eqcon[124] *= 1;  */

    /*reaction 126: CH + H2O <=> H + CH2O */
    /*eqcon[125] *= 1;  */

    /*reaction 127: CH + CH2 <=> H + C2H2 */
    /*eqcon[126] *= 1;  */

    /*reaction 128: CH + CH3 <=> H + C2H3 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: CH + CH4 <=> H + C2H4 */
    /*eqcon[128] *= 1;  */

    /*reaction 130: CH + CO (+M) <=> HCCO (+M) */
    eqcon[129] *= 1e+06; 

    /*reaction 131: CH + CO2 <=> HCO + CO */
    /*eqcon[130] *= 1;  */

    /*reaction 132: CH + CH2O <=> H + CH2CO */
    /*eqcon[131] *= 1;  */

    /*reaction 133: CH + HCCO <=> CO + C2H2 */
    /*eqcon[132] *= 1;  */

    /*reaction 134: CH2 + O2 => OH + H + CO */
    eqcon[133] *= 1e-06; 

    /*reaction 135: CH2 + H2 <=> H + CH3 */
    /*eqcon[134] *= 1;  */

    /*reaction 136: 2 CH2 <=> H2 + C2H2 */
    /*eqcon[135] *= 1;  */

    /*reaction 137: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[137] *= 1;  */

    /*reaction 139: CH2 + CO (+M) <=> CH2CO (+M) */
    eqcon[138] *= 1e+06; 

    /*reaction 140: CH2 + HCCO <=> C2H3 + CO */
    /*eqcon[139] *= 1;  */

    /*reaction 141: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[140] *= 1;  */

    /*reaction 142: CH2(S) + O2 <=> H + OH + CO */
    eqcon[141] *= 1e-06; 

    /*reaction 143: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[142] *= 1;  */

    /*reaction 144: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[143] *= 1;  */

    /*reaction 145: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    eqcon[144] *= 1e+06; 

    /*reaction 146: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[145] *= 1;  */

    /*reaction 147: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[146] *= 1;  */

    /*reaction 148: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[147] *= 1;  */

    /*reaction 149: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[148] *= 1;  */

    /*reaction 150: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[149] *= 1;  */

    /*reaction 151: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[150] *= 1;  */

    /*reaction 152: CH2(S) + C2H6 <=> CH3 + C2H5 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: CH3 + O2 <=> O + CH3O */
    /*eqcon[152] *= 1;  */

    /*reaction 154: CH3 + O2 <=> OH + CH2O */
    /*eqcon[153] *= 1;  */

    /*reaction 155: CH3 + H2O2 <=> HO2 + CH4 */
    /*eqcon[154] *= 1;  */

    /*reaction 156: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[155] *= 1e+06; 

    /*reaction 157: 2 CH3 <=> H + C2H5 */
    /*eqcon[156] *= 1;  */

    /*reaction 158: CH3 + HCO <=> CH4 + CO */
    /*eqcon[157] *= 1;  */

    /*reaction 159: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[158] *= 1;  */

    /*reaction 160: CH3 + CH3OH <=> CH2OH + CH4 */
    /*eqcon[159] *= 1;  */

    /*reaction 161: CH3 + CH3OH <=> CH3O + CH4 */
    /*eqcon[160] *= 1;  */

    /*reaction 162: CH3 + C2H4 <=> C2H3 + CH4 */
    /*eqcon[161] *= 1;  */

    /*reaction 163: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[162] *= 1;  */

    /*reaction 164: HCO + H2O <=> H + CO + H2O */
    eqcon[163] *= 1e-06; 

    /*reaction 165: HCO + M <=> H + CO + M */
    eqcon[164] *= 1e-06; 

    /*reaction 166: HCO + O2 <=> HO2 + CO */
    /*eqcon[165] *= 1;  */

    /*reaction 167: CH2OH + O2 <=> HO2 + CH2O */
    /*eqcon[166] *= 1;  */

    /*reaction 168: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[167] *= 1;  */

    /*reaction 169: C2H + O2 <=> HCO + CO */
    /*eqcon[168] *= 1;  */

    /*reaction 170: C2H + H2 <=> H + C2H2 */
    /*eqcon[169] *= 1;  */

    /*reaction 171: C2H3 + O2 <=> HCO + CH2O */
    /*eqcon[170] *= 1;  */

    /*reaction 172: C2H4 (+M) <=> H2 + C2H2 (+M) */
    eqcon[171] *= 1e-06; 

    /*reaction 173: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[172] *= 1;  */

    /*reaction 174: HCCO + O2 <=> OH + 2 CO */
    eqcon[173] *= 1e-06; 

    /*reaction 175: 2 HCCO <=> 2 CO + C2H2 */
    eqcon[174] *= 1e-06; 

    /*reaction 176: O + CH3 => H + H2 + CO */
    eqcon[175] *= 1e-06; 

    /*reaction 177: O + C2H4 <=> H + CH2CHO */
    /*eqcon[176] *= 1;  */

    /*reaction 178: O + C2H5 <=> H + CH3CHO */
    /*eqcon[177] *= 1;  */

    /*reaction 179: OH + HO2 <=> O2 + H2O */
    /*eqcon[178] *= 1;  */

    /*reaction 180: OH + CH3 => H2 + CH2O */
    /*eqcon[179] *= 1;  */

    /*reaction 181: CH + H2 (+M) <=> CH3 (+M) */
    eqcon[180] *= 1e+06; 

    /*reaction 182: CH2 + O2 => 2 H + CO2 */
    eqcon[181] *= 1e-06; 

    /*reaction 183: CH2 + O2 <=> O + CH2O */
    /*eqcon[182] *= 1;  */

    /*reaction 184: CH2 + CH2 => 2 H + C2H2 */
    eqcon[183] *= 1e-06; 

    /*reaction 185: CH2(S) + H2O => H2 + CH2O */
    /*eqcon[184] *= 1;  */

    /*reaction 186: C2H3 + O2 <=> O + CH2CHO */
    /*eqcon[185] *= 1;  */

    /*reaction 187: C2H3 + O2 <=> HO2 + C2H2 */
    /*eqcon[186] *= 1;  */

    /*reaction 188: O + CH3CHO <=> OH + CH2CHO */
    /*eqcon[187] *= 1;  */

    /*reaction 189: O + CH3CHO => OH + CH3 + CO */
    eqcon[188] *= 1e-06; 

    /*reaction 190: O2 + CH3CHO => HO2 + CH3 + CO */
    eqcon[189] *= 1e-06; 

    /*reaction 191: H + CH3CHO <=> CH2CHO + H2 */
    /*eqcon[190] *= 1;  */

    /*reaction 192: H + CH3CHO => CH3 + H2 + CO */
    eqcon[191] *= 1e-06; 

    /*reaction 193: OH + CH3CHO => CH3 + H2O + CO */
    eqcon[192] *= 1e-06; 

    /*reaction 194: HO2 + CH3CHO => CH3 + H2O2 + CO */
    eqcon[193] *= 1e-06; 

    /*reaction 195: CH3 + CH3CHO => CH3 + CH4 + CO */
    eqcon[194] *= 1e-06; 

    /*reaction 196: H + CH2CO (+M) <=> CH2CHO (+M) */
    eqcon[195] *= 1e+06; 

    /*reaction 197: O + CH2CHO => H + CH2 + CO2 */
    eqcon[196] *= 1e-06; 

    /*reaction 198: O2 + CH2CHO => OH + CO + CH2O */
    eqcon[197] *= 1e-06; 

    /*reaction 199: O2 + CH2CHO => OH + 2 HCO */
    eqcon[198] *= 1e-06; 

    /*reaction 200: H + CH2CHO <=> CH3 + HCO */
    /*eqcon[199] *= 1;  */

    /*reaction 201: H + CH2CHO <=> CH2CO + H2 */
    /*eqcon[200] *= 1;  */

    /*reaction 202: OH + CH2CHO <=> H2O + CH2CO */
    /*eqcon[201] *= 1;  */

    /*reaction 203: OH + CH2CHO <=> HCO + CH2OH */
    /*eqcon[202] *= 1;  */

    /*reaction 204: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    eqcon[203] *= 1e+06; 

    /*reaction 205: O + C3H8 <=> OH + C3H7 */
    /*eqcon[204] *= 1;  */

    /*reaction 206: H + C3H8 <=> C3H7 + H2 */
    /*eqcon[205] *= 1;  */

    /*reaction 207: OH + C3H8 <=> C3H7 + H2O */
    /*eqcon[206] *= 1;  */

    /*reaction 208: C3H7 + H2O2 <=> HO2 + C3H8 */
    /*eqcon[207] *= 1;  */

    /*reaction 209: CH3 + C3H8 <=> C3H7 + CH4 */
    /*eqcon[208] *= 1;  */

    /*reaction 210: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    eqcon[209] *= 1e+06; 

    /*reaction 211: O + C3H7 <=> C2H5 + CH2O */
    /*eqcon[210] *= 1;  */

    /*reaction 212: H + C3H7 (+M) <=> C3H8 (+M) */
    eqcon[211] *= 1e+06; 

    /*reaction 213: H + C3H7 <=> CH3 + C2H5 */
    /*eqcon[212] *= 1;  */

    /*reaction 214: OH + C3H7 <=> C2H5 + CH2OH */
    /*eqcon[213] *= 1;  */

    /*reaction 215: HO2 + C3H7 <=> O2 + C3H8 */
    /*eqcon[214] *= 1;  */

    /*reaction 216: HO2 + C3H7 => OH + C2H5 + CH2O */
    eqcon[215] *= 1e-06; 

    /*reaction 217: CH3 + C3H7 <=> 2 C2H5 */
    /*eqcon[216] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[35]; /* temporary storage */

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

    /*reaction 37: H + O2 <=> O + OH */
    /*eqcon[36] *= 1;  */

    /*reaction 38: 2 H + M <=> H2 + M */
    eqcon[37] *= 1e+06; 

    /*reaction 39: 2 H + H2 <=> 2 H2 */
    eqcon[38] *= 1e+06; 

    /*reaction 40: 2 H + H2O <=> H2 + H2O */
    eqcon[39] *= 1e+06; 

    /*reaction 41: 2 H + CO2 <=> H2 + CO2 */
    eqcon[40] *= 1e+06; 

    /*reaction 42: H + OH + M <=> H2O + M */
    eqcon[41] *= 1e+06; 

    /*reaction 43: H + HO2 <=> O + H2O */
    /*eqcon[42] *= 1;  */

    /*reaction 44: H + HO2 <=> O2 + H2 */
    /*eqcon[43] *= 1;  */

    /*reaction 45: H + HO2 <=> 2 OH */
    /*eqcon[44] *= 1;  */

    /*reaction 46: H + H2O2 <=> HO2 + H2 */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + H2O2 <=> OH + H2O */
    /*eqcon[46] *= 1;  */

    /*reaction 48: H + CH <=> C + H2 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[48] *= 1e+06; 

    /*reaction 50: H + CH2(S) <=> CH + H2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[50] *= 1e+06; 

    /*reaction 52: H + CH4 <=> CH3 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: H + HCO (+M) <=> CH2O (+M) */
    eqcon[52] *= 1e+06; 

    /*reaction 54: H + HCO <=> H2 + CO */
    /*eqcon[53] *= 1;  */

    /*reaction 55: H + CH2O (+M) <=> CH2OH (+M) */
    eqcon[54] *= 1e+06; 

    /*reaction 56: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[55] *= 1e+06; 

    /*reaction 57: H + CH2O <=> HCO + H2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: H + CH2OH (+M) <=> CH3OH (+M) */
    eqcon[57] *= 1e+06; 

    /*reaction 59: H + CH2OH <=> H2 + CH2O */
    /*eqcon[58] *= 1;  */

    /*reaction 60: H + CH2OH <=> OH + CH3 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: H + CH2OH <=> CH2(S) + H2O */
    /*eqcon[60] *= 1;  */

    /*reaction 62: H + CH3O (+M) <=> CH3OH (+M) */
    eqcon[61] *= 1e+06; 

    /*reaction 63: H + CH3O <=> H + CH2OH */
    /*eqcon[62] *= 1;  */

    /*reaction 64: H + CH3O <=> H2 + CH2O */
    /*eqcon[63] *= 1;  */

    /*reaction 65: H + CH3O <=> OH + CH3 */
    /*eqcon[64] *= 1;  */

    /*reaction 66: H + CH3O <=> CH2(S) + H2O */
    /*eqcon[65] *= 1;  */

    /*reaction 67: H + CH3OH <=> CH2OH + H2 */
    /*eqcon[66] *= 1;  */

    /*reaction 68: H + CH3OH <=> CH3O + H2 */
    /*eqcon[67] *= 1;  */

    /*reaction 69: H + C2H (+M) <=> C2H2 (+M) */
    eqcon[68] *= 1e+06; 

    /*reaction 70: H + C2H2 (+M) <=> C2H3 (+M) */
    eqcon[69] *= 1e+06; 

    /*reaction 71: H + C2H3 (+M) <=> C2H4 (+M) */
    eqcon[70] *= 1e+06; 

    /*reaction 72: H + C2H3 <=> H2 + C2H2 */
    /*eqcon[71] *= 1;  */

    /*reaction 73: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[72] *= 1e+06; 

    /*reaction 74: H + C2H4 <=> C2H3 + H2 */
    /*eqcon[73] *= 1;  */

    /*reaction 75: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[74] *= 1e+06; 

    /*reaction 76: H + C2H5 <=> H2 + C2H4 */
    /*eqcon[75] *= 1;  */

    /*reaction 77: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: H + HCCO <=> CH2(S) + CO */
    /*eqcon[77] *= 1;  */

    /*reaction 79: H + CH2CO <=> HCCO + H2 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: H + CH2CO <=> CH3 + CO */
    /*eqcon[79] *= 1;  */

    /*reaction 81: H + HCCOH <=> H + CH2CO */
    /*eqcon[80] *= 1;  */

    /*reaction 82: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[81] *= 1e+06; 

    /*reaction 83: OH + H2 <=> H + H2O */
    /*eqcon[82] *= 1;  */

    /*reaction 84: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[83] *= 1e+06; 

    /*reaction 85: 2 OH <=> O + H2O */
    /*eqcon[84] *= 1;  */

    /*reaction 86: OH + HO2 <=> O2 + H2O */
    /*eqcon[85] *= 1;  */

    /*reaction 87: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[86] *= 1;  */

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[87] *= 1;  */

    /*reaction 89: OH + C <=> H + CO */
    /*eqcon[88] *= 1;  */

    /*reaction 90: OH + CH <=> H + HCO */
    /*eqcon[89] *= 1;  */

    /*reaction 91: OH + CH2 <=> H + CH2O */
    /*eqcon[90] *= 1;  */

    /*reaction 92: OH + CH2 <=> CH + H2O */
    /*eqcon[91] *= 1;  */

    /*reaction 93: OH + CH2(S) <=> H + CH2O */
    /*eqcon[92] *= 1;  */

    /*reaction 94: OH + CH3 (+M) <=> CH3OH (+M) */
    eqcon[93] *= 1e+06; 

    /*reaction 95: OH + CH3 <=> CH2 + H2O */
    /*eqcon[94] *= 1;  */

    /*reaction 96: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[95] *= 1;  */

    /*reaction 97: OH + CH4 <=> CH3 + H2O */
    /*eqcon[96] *= 1;  */

    /*reaction 98: OH + CO <=> H + CO2 */
    /*eqcon[97] *= 1;  */

    /*reaction 99: OH + HCO <=> H2O + CO */
    /*eqcon[98] *= 1;  */

    /*reaction 100: OH + CH2O <=> HCO + H2O */
    /*eqcon[99] *= 1;  */

    /*reaction 101: OH + CH2OH <=> H2O + CH2O */
    /*eqcon[100] *= 1;  */

    /*reaction 102: OH + CH3O <=> H2O + CH2O */
    /*eqcon[101] *= 1;  */

    /*reaction 103: OH + CH3OH <=> CH2OH + H2O */
    /*eqcon[102] *= 1;  */

    /*reaction 104: OH + CH3OH <=> CH3O + H2O */
    /*eqcon[103] *= 1;  */

    /*reaction 105: OH + C2H <=> H + HCCO */
    /*eqcon[104] *= 1;  */

    /*reaction 106: OH + C2H2 <=> H + CH2CO */
    /*eqcon[105] *= 1;  */

    /*reaction 107: OH + C2H2 <=> H + HCCOH */
    /*eqcon[106] *= 1;  */

    /*reaction 108: OH + C2H2 <=> C2H + H2O */
    /*eqcon[107] *= 1;  */

    /*reaction 109: OH + C2H2 <=> CH3 + CO */
    /*eqcon[108] *= 1;  */

    /*reaction 110: OH + C2H3 <=> H2O + C2H2 */
    /*eqcon[109] *= 1;  */

    /*reaction 111: OH + C2H4 <=> C2H3 + H2O */
    /*eqcon[110] *= 1;  */

    /*reaction 112: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[111] *= 1;  */

    /*reaction 113: OH + CH2CO <=> HCCO + H2O */
    /*eqcon[112] *= 1;  */

    /*reaction 114: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[113] *= 1;  */

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[114] *= 1;  */

    /*reaction 116: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[115] *= 1;  */

    /*reaction 117: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[116] *= 1;  */

    /*reaction 118: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[117] *= 1;  */

    /*reaction 119: HO2 + CO <=> OH + CO2 */
    /*eqcon[118] *= 1;  */

    /*reaction 120: HO2 + CH2O <=> HCO + H2O2 */
    /*eqcon[119] *= 1;  */

    /*reaction 121: C + O2 <=> O + CO */
    /*eqcon[120] *= 1;  */

    /*reaction 122: C + CH2 <=> H + C2H */
    /*eqcon[121] *= 1;  */

    /*reaction 123: C + CH3 <=> H + C2H2 */
    /*eqcon[122] *= 1;  */

    /*reaction 124: CH + O2 <=> O + HCO */
    /*eqcon[123] *= 1;  */

    /*reaction 125: CH + H2 <=> H + CH2 */
    /*eqcon[124] *= 1;  */

    /*reaction 126: CH + H2O <=> H + CH2O */
    /*eqcon[125] *= 1;  */

    /*reaction 127: CH + CH2 <=> H + C2H2 */
    /*eqcon[126] *= 1;  */

    /*reaction 128: CH + CH3 <=> H + C2H3 */
    /*eqcon[127] *= 1;  */

    /*reaction 129: CH + CH4 <=> H + C2H4 */
    /*eqcon[128] *= 1;  */

    /*reaction 130: CH + CO (+M) <=> HCCO (+M) */
    eqcon[129] *= 1e+06; 

    /*reaction 131: CH + CO2 <=> HCO + CO */
    /*eqcon[130] *= 1;  */

    /*reaction 132: CH + CH2O <=> H + CH2CO */
    /*eqcon[131] *= 1;  */

    /*reaction 133: CH + HCCO <=> CO + C2H2 */
    /*eqcon[132] *= 1;  */

    /*reaction 134: CH2 + O2 => OH + H + CO */
    eqcon[133] *= 1e-06; 

    /*reaction 135: CH2 + H2 <=> H + CH3 */
    /*eqcon[134] *= 1;  */

    /*reaction 136: 2 CH2 <=> H2 + C2H2 */
    /*eqcon[135] *= 1;  */

    /*reaction 137: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[136] *= 1;  */

    /*reaction 138: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[137] *= 1;  */

    /*reaction 139: CH2 + CO (+M) <=> CH2CO (+M) */
    eqcon[138] *= 1e+06; 

    /*reaction 140: CH2 + HCCO <=> C2H3 + CO */
    /*eqcon[139] *= 1;  */

    /*reaction 141: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[140] *= 1;  */

    /*reaction 142: CH2(S) + O2 <=> H + OH + CO */
    eqcon[141] *= 1e-06; 

    /*reaction 143: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[142] *= 1;  */

    /*reaction 144: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[143] *= 1;  */

    /*reaction 145: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    eqcon[144] *= 1e+06; 

    /*reaction 146: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[145] *= 1;  */

    /*reaction 147: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[146] *= 1;  */

    /*reaction 148: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[147] *= 1;  */

    /*reaction 149: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[148] *= 1;  */

    /*reaction 150: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[149] *= 1;  */

    /*reaction 151: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[150] *= 1;  */

    /*reaction 152: CH2(S) + C2H6 <=> CH3 + C2H5 */
    /*eqcon[151] *= 1;  */

    /*reaction 153: CH3 + O2 <=> O + CH3O */
    /*eqcon[152] *= 1;  */

    /*reaction 154: CH3 + O2 <=> OH + CH2O */
    /*eqcon[153] *= 1;  */

    /*reaction 155: CH3 + H2O2 <=> HO2 + CH4 */
    /*eqcon[154] *= 1;  */

    /*reaction 156: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[155] *= 1e+06; 

    /*reaction 157: 2 CH3 <=> H + C2H5 */
    /*eqcon[156] *= 1;  */

    /*reaction 158: CH3 + HCO <=> CH4 + CO */
    /*eqcon[157] *= 1;  */

    /*reaction 159: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[158] *= 1;  */

    /*reaction 160: CH3 + CH3OH <=> CH2OH + CH4 */
    /*eqcon[159] *= 1;  */

    /*reaction 161: CH3 + CH3OH <=> CH3O + CH4 */
    /*eqcon[160] *= 1;  */

    /*reaction 162: CH3 + C2H4 <=> C2H3 + CH4 */
    /*eqcon[161] *= 1;  */

    /*reaction 163: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[162] *= 1;  */

    /*reaction 164: HCO + H2O <=> H + CO + H2O */
    eqcon[163] *= 1e-06; 

    /*reaction 165: HCO + M <=> H + CO + M */
    eqcon[164] *= 1e-06; 

    /*reaction 166: HCO + O2 <=> HO2 + CO */
    /*eqcon[165] *= 1;  */

    /*reaction 167: CH2OH + O2 <=> HO2 + CH2O */
    /*eqcon[166] *= 1;  */

    /*reaction 168: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[167] *= 1;  */

    /*reaction 169: C2H + O2 <=> HCO + CO */
    /*eqcon[168] *= 1;  */

    /*reaction 170: C2H + H2 <=> H + C2H2 */
    /*eqcon[169] *= 1;  */

    /*reaction 171: C2H3 + O2 <=> HCO + CH2O */
    /*eqcon[170] *= 1;  */

    /*reaction 172: C2H4 (+M) <=> H2 + C2H2 (+M) */
    eqcon[171] *= 1e-06; 

    /*reaction 173: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[172] *= 1;  */

    /*reaction 174: HCCO + O2 <=> OH + 2 CO */
    eqcon[173] *= 1e-06; 

    /*reaction 175: 2 HCCO <=> 2 CO + C2H2 */
    eqcon[174] *= 1e-06; 

    /*reaction 176: O + CH3 => H + H2 + CO */
    eqcon[175] *= 1e-06; 

    /*reaction 177: O + C2H4 <=> H + CH2CHO */
    /*eqcon[176] *= 1;  */

    /*reaction 178: O + C2H5 <=> H + CH3CHO */
    /*eqcon[177] *= 1;  */

    /*reaction 179: OH + HO2 <=> O2 + H2O */
    /*eqcon[178] *= 1;  */

    /*reaction 180: OH + CH3 => H2 + CH2O */
    /*eqcon[179] *= 1;  */

    /*reaction 181: CH + H2 (+M) <=> CH3 (+M) */
    eqcon[180] *= 1e+06; 

    /*reaction 182: CH2 + O2 => 2 H + CO2 */
    eqcon[181] *= 1e-06; 

    /*reaction 183: CH2 + O2 <=> O + CH2O */
    /*eqcon[182] *= 1;  */

    /*reaction 184: CH2 + CH2 => 2 H + C2H2 */
    eqcon[183] *= 1e-06; 

    /*reaction 185: CH2(S) + H2O => H2 + CH2O */
    /*eqcon[184] *= 1;  */

    /*reaction 186: C2H3 + O2 <=> O + CH2CHO */
    /*eqcon[185] *= 1;  */

    /*reaction 187: C2H3 + O2 <=> HO2 + C2H2 */
    /*eqcon[186] *= 1;  */

    /*reaction 188: O + CH3CHO <=> OH + CH2CHO */
    /*eqcon[187] *= 1;  */

    /*reaction 189: O + CH3CHO => OH + CH3 + CO */
    eqcon[188] *= 1e-06; 

    /*reaction 190: O2 + CH3CHO => HO2 + CH3 + CO */
    eqcon[189] *= 1e-06; 

    /*reaction 191: H + CH3CHO <=> CH2CHO + H2 */
    /*eqcon[190] *= 1;  */

    /*reaction 192: H + CH3CHO => CH3 + H2 + CO */
    eqcon[191] *= 1e-06; 

    /*reaction 193: OH + CH3CHO => CH3 + H2O + CO */
    eqcon[192] *= 1e-06; 

    /*reaction 194: HO2 + CH3CHO => CH3 + H2O2 + CO */
    eqcon[193] *= 1e-06; 

    /*reaction 195: CH3 + CH3CHO => CH3 + CH4 + CO */
    eqcon[194] *= 1e-06; 

    /*reaction 196: H + CH2CO (+M) <=> CH2CHO (+M) */
    eqcon[195] *= 1e+06; 

    /*reaction 197: O + CH2CHO => H + CH2 + CO2 */
    eqcon[196] *= 1e-06; 

    /*reaction 198: O2 + CH2CHO => OH + CO + CH2O */
    eqcon[197] *= 1e-06; 

    /*reaction 199: O2 + CH2CHO => OH + 2 HCO */
    eqcon[198] *= 1e-06; 

    /*reaction 200: H + CH2CHO <=> CH3 + HCO */
    /*eqcon[199] *= 1;  */

    /*reaction 201: H + CH2CHO <=> CH2CO + H2 */
    /*eqcon[200] *= 1;  */

    /*reaction 202: OH + CH2CHO <=> H2O + CH2CO */
    /*eqcon[201] *= 1;  */

    /*reaction 203: OH + CH2CHO <=> HCO + CH2OH */
    /*eqcon[202] *= 1;  */

    /*reaction 204: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    eqcon[203] *= 1e+06; 

    /*reaction 205: O + C3H8 <=> OH + C3H7 */
    /*eqcon[204] *= 1;  */

    /*reaction 206: H + C3H8 <=> C3H7 + H2 */
    /*eqcon[205] *= 1;  */

    /*reaction 207: OH + C3H8 <=> C3H7 + H2O */
    /*eqcon[206] *= 1;  */

    /*reaction 208: C3H7 + H2O2 <=> HO2 + C3H8 */
    /*eqcon[207] *= 1;  */

    /*reaction 209: CH3 + C3H8 <=> C3H7 + CH4 */
    /*eqcon[208] *= 1;  */

    /*reaction 210: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    eqcon[209] *= 1e+06; 

    /*reaction 211: O + C3H7 <=> C2H5 + CH2O */
    /*eqcon[210] *= 1;  */

    /*reaction 212: H + C3H7 (+M) <=> C3H8 (+M) */
    eqcon[211] *= 1e+06; 

    /*reaction 213: H + C3H7 <=> CH3 + C2H5 */
    /*eqcon[212] *= 1;  */

    /*reaction 214: OH + C3H7 <=> C2H5 + CH2OH */
    /*eqcon[213] *= 1;  */

    /*reaction 215: HO2 + C3H7 <=> O2 + C3H8 */
    /*eqcon[214] *= 1;  */

    /*reaction 216: HO2 + C3H7 => OH + C2H5 + CH2O */
    eqcon[215] *= 1e-06; 

    /*reaction 217: CH3 + C3H7 <=> 2 C2H5 */
    /*eqcon[216] *= 1;  */
}

static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[217];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[217];
#ifdef _OPENMP
#pragma omp threadprivate(Kc_save)
#endif

/*compute the production rate for each species */
void productionRate(double * wdot, double * sc, double T)
{
    double qdot;

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[35];                /*Gibbs free energy */
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
    for (id = 0; id < 35; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*zero out wdot */
    for (id = 0; id < 35; ++id) {
        wdot[id] = 0.0;
    }

    if (T != T_save)
    {
        T_save = T;

        k_f_save[0] = 1e-12 * 1.2e+17*exp(-1*tc[0]);
        k_f_save[1] = 1e-12 * 5e+17*exp(-1*tc[0]);
        k_f_save[2] = 1e-06 * 38700*exp(2.7*tc[0]-3150.1363279375455022*invT);
        k_f_save[3] = 1e-06 * 2e+13;
        k_f_save[4] = 1e-06 * 9.63e+06*exp(2*tc[0]-2012.866663218878557*invT);
        k_f_save[5] = 1e-06 * 5.7e+13;
        k_f_save[6] = 1e-06 * 8e+13;
        k_f_save[7] = 1e-06 * 1.5e+13;
        k_f_save[8] = 1e-06 * 1.5e+13;
        k_f_save[9] = 1e-06 * 5.06e+13;
        k_f_save[10] = 1e-06 * 1.02e+09*exp(1.5*tc[0]-4327.6633259205891591*invT);
        k_f_save[11] = 1e-06 * 1.8e+10*exp(-1200.1717479442565946*invT);
        k_f_save[12] = 1e-06 * 3e+13;
        k_f_save[13] = 1e-06 * 3e+13;
        k_f_save[14] = 1e-06 * 3.9e+13*exp(-1781.3869969487077469*invT);
        k_f_save[15] = 1e-06 * 1e+13;
        k_f_save[16] = 1e-06 * 1e+13;
        k_f_save[17] = 1e-06 * 388000*exp(2.5*tc[0]-1559.9716639946311716*invT);
        k_f_save[18] = 1e-06 * 130000*exp(2.5*tc[0]-2516.0833290235987079*invT);
        k_f_save[19] = 1e-06 * 5e+13;
        k_f_save[20] = 1e-06 * 1.35e+07*exp(2*tc[0]-956.11166502896742259*invT);
        k_f_save[21] = 1e-06 * 4.6e+19*exp(-1.41*tc[0]-14568.122475046635373*invT);
        k_f_save[22] = 1e-06 * 6.94e+06*exp(2*tc[0]-956.11166502896742259*invT);
        k_f_save[23] = 1e-06 * 3e+13;
        k_f_save[24] = 1e-06 * 1.25e+07*exp(1.83*tc[0]-110.7076664770383303*invT);
        k_f_save[25] = 1e-06 * 2.24e+13;
        k_f_save[26] = 1e-06 * 8.98e+07*exp(1.92*tc[0]-2863.3028284288548093*invT);
        k_f_save[27] = 1e-06 * 1e+14;
        k_f_save[28] = 1e-06 * 1e+13*exp(-4025.733326437757114*invT);
        k_f_save[29] = 1e-06 * 1.75e+12*exp(-679.34249883637164658*invT);
        k_f_save[30] = 1e-06 * 2.5e+12*exp(-24053.756625465601246*invT);
        k_f_save[31] = 1e-06 * 1e+14*exp(-20128.666632188789663*invT);
        k_f_save[32] = 1e-12 * 2.8e+18*exp(-0.86*tc[0]);
        k_f_save[33] = 1e-12 * 2.08e+19*exp(-1.24*tc[0]);
        k_f_save[34] = 1e-12 * 1.126e+19*exp(-0.76*tc[0]);
        k_f_save[35] = 1e-12 * 2.6e+19*exp(-1.24*tc[0]);
        k_f_save[36] = 1e-06 * 2.65e+16*exp(-0.6707*tc[0]-8575.3152019782282878*invT);
        k_f_save[37] = 1e-12 * 1e+18*exp(-1*tc[0]);
        k_f_save[38] = 1e-12 * 9e+16*exp(-0.6*tc[0]);
        k_f_save[39] = 1e-12 * 6e+19*exp(-1.25*tc[0]);
        k_f_save[40] = 1e-12 * 5.5e+20*exp(-2*tc[0]);
        k_f_save[41] = 1e-12 * 2.2e+22*exp(-2*tc[0]);
        k_f_save[42] = 1e-06 * 3.97e+12*exp(-337.65838275496690812*invT);
        k_f_save[43] = 1e-06 * 4.48e+13*exp(-537.43539907944057177*invT);
        k_f_save[44] = 1e-06 * 8.4e+13*exp(-319.54258278599701271*invT);
        k_f_save[45] = 1e-06 * 1.21e+07*exp(2*tc[0]-2616.7266621845424197*invT);
        k_f_save[46] = 1e-06 * 1e+13*exp(-1811.579996896990906*invT);
        k_f_save[47] = 1e-06 * 1.65e+14;
        k_f_save[48] = 1e-06 * 6e+14;
        k_f_save[49] = 1e-06 * 3e+13;
        k_f_save[50] = 1e-06 * 1.39e+16*exp(-0.534*tc[0]-269.72413287132980031*invT);
        k_f_save[51] = 1e-06 * 6.6e+08*exp(1.62*tc[0]-5454.8686573231616421*invT);
        k_f_save[52] = 1e-06 * 1.09e+12*exp(0.48*tc[0]+130.83633310922709825*invT);
        k_f_save[53] = 1e-06 * 7.34e+13;
        k_f_save[54] = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1811.579996896990906*invT);
        k_f_save[55] = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1308.3633310922712099*invT);
        k_f_save[56] = 1e-06 * 5.74e+07*exp(1.9*tc[0]-1379.8200976365412771*invT);
        k_f_save[57] = 1e-06 * 1.055e+12*exp(0.5*tc[0]-43.276633259205894433*invT);
        k_f_save[58] = 1e-06 * 2e+13;
        k_f_save[59] = 1e-06 * 1.65e+11*exp(0.65*tc[0]+142.9135330885404187*invT);
        k_f_save[60] = 1e-06 * 3.28e+13*exp(-0.09*tc[0]-306.96216614087899188*invT);
        k_f_save[61] = 1e-06 * 2.43e+12*exp(0.515*tc[0]-25.160833290235984805*invT);
        k_f_save[62] = 1e-06 * 4.15e+07*exp(1.63*tc[0]-968.1888650082806862*invT);
        k_f_save[63] = 1e-06 * 2e+13;
        k_f_save[64] = 1e-06 * 1.5e+12*exp(0.5*tc[0]+55.35383323851916515*invT);
        k_f_save[65] = 1e-06 * 2.62e+14*exp(-0.23*tc[0]-538.44183241105008619*invT);
        k_f_save[66] = 1e-06 * 1.7e+07*exp(2.1*tc[0]-2450.6651624689848177*invT);
        k_f_save[67] = 1e-06 * 4.2e+06*exp(2.1*tc[0]-2450.6651624689848177*invT);
        k_f_save[68] = 1e-06 * 1e+17*exp(-1*tc[0]);
        k_f_save[69] = 1e-06 * 5.6e+12*exp(-1207.7199979313272706*invT);
        k_f_save[70] = 1e-06 * 6.08e+12*exp(0.27*tc[0]-140.90066642532150354*invT);
        k_f_save[71] = 1e-06 * 3e+13;
        k_f_save[72] = 1e-06 * 5.4e+11*exp(0.454*tc[0]-915.85433176458980142*invT);
        k_f_save[73] = 1e-06 * 1.325e+06*exp(2.53*tc[0]-6159.3719894497689893*invT);
        k_f_save[74] = 1e-06 * 5.21e+17*exp(-0.99*tc[0]-795.08233197145705162*invT);
        k_f_save[75] = 1e-06 * 2e+12;
        k_f_save[76] = 1e-06 * 1.15e+08*exp(1.9*tc[0]-3789.2214935095394139*invT);
        k_f_save[77] = 1e-06 * 1e+14;
        k_f_save[78] = 1e-06 * 5e+13*exp(-4025.733326437757114*invT);
        k_f_save[79] = 1e-06 * 1.13e+13*exp(-1725.0267303785790318*invT);
        k_f_save[80] = 1e-06 * 1e+13;
        k_f_save[81] = 1e-06 * 4.3e+07*exp(1.5*tc[0]-40056.046598055690993*invT);
        k_f_save[82] = 1e-06 * 2.16e+08*exp(1.51*tc[0]-1726.0331637101885462*invT);
        k_f_save[83] = 1e-06 * 7.4e+13*exp(-0.37*tc[0]);
        k_f_save[84] = 1e-06 * 35700*exp(2.4*tc[0]+1061.7871648479585929*invT);
        k_f_save[85] = 1e-06 * 1.45e+13*exp(+251.60833290235981963*invT);
        k_f_save[86] = 1e-06 * 2e+12*exp(-214.8735162986153*invT);
        k_f_save[87] = 1e-06 * 1.7e+18*exp(-14799.602141316805501*invT);
        k_f_save[88] = 1e-06 * 5e+13;
        k_f_save[89] = 1e-06 * 3e+13;
        k_f_save[90] = 1e-06 * 2e+13;
        k_f_save[91] = 1e-06 * 1.13e+07*exp(2*tc[0]-1509.6499974141590883*invT);
        k_f_save[92] = 1e-06 * 3e+13;
        k_f_save[93] = 1e-06 * 2.79e+18*exp(-1.43*tc[0]-669.27816552027718444*invT);
        k_f_save[94] = 1e-06 * 5.6e+07*exp(1.6*tc[0]-2727.4343286615808211*invT);
        k_f_save[95] = 1e-06 * 6.44e+17*exp(-1.34*tc[0]-713.05801544528776503*invT);
        k_f_save[96] = 1e-06 * 1e+08*exp(1.6*tc[0]-1570.0359973107254064*invT);
        k_f_save[97] = 1e-06 * 4.76e+07*exp(1.228*tc[0]-35.225166606330375885*invT);
        k_f_save[98] = 1e-06 * 5e+13;
        k_f_save[99] = 1e-06 * 3.43e+09*exp(1.18*tc[0]+224.93784961470970529*invT);
        k_f_save[100] = 1e-06 * 5e+12;
        k_f_save[101] = 1e-06 * 5e+12;
        k_f_save[102] = 1e-06 * 1.44e+06*exp(2*tc[0]+422.70199927596451062*invT);
        k_f_save[103] = 1e-06 * 6.3e+06*exp(2*tc[0]-754.82499870707954415*invT);
        k_f_save[104] = 1e-06 * 2e+13;
        k_f_save[105] = 1e-06 * 0.000218*exp(4.5*tc[0]+503.21666580471963925*invT);
        k_f_save[106] = 1e-06 * 504000*exp(2.3*tc[0]-6793.4249883637157836*invT);
        k_f_save[107] = 1e-06 * 3.37e+07*exp(2*tc[0]-7045.0333212660752906*invT);
        k_f_save[108] = 1e-06 * 0.000483*exp(4*tc[0]+1006.4333316094392785*invT);
        k_f_save[109] = 1e-06 * 5e+12;
        k_f_save[110] = 1e-06 * 3.6e+06*exp(2*tc[0]-1258.0416645117993539*invT);
        k_f_save[111] = 1e-06 * 3.54e+06*exp(2.12*tc[0]-437.79849925010609013*invT);
        k_f_save[112] = 1e-06 * 7.5e+12*exp(-1006.4333316094392785*invT);
        k_f_save[113] = 1e-06 * 1.3e+11*exp(+820.24316526169309327*invT);
        k_f_save[114] = 1e-06 * 4.2e+14*exp(-6038.5999896566363532*invT);
        k_f_save[115] = 1e-06 * 2e+13;
        k_f_save[116] = 1e-06 * 1e+12;
        k_f_save[117] = 1e-06 * 3.78e+13;
        k_f_save[118] = 1e-06 * 1.5e+14*exp(-11875.913312991384373*invT);
        k_f_save[119] = 1e-06 * 5.6e+06*exp(2*tc[0]-6038.5999896566363532*invT);
        k_f_save[120] = 1e-06 * 5.8e+13*exp(-289.85279950351855405*invT);
        k_f_save[121] = 1e-06 * 5e+13;
        k_f_save[122] = 1e-06 * 5e+13;
        k_f_save[123] = 1e-06 * 6.71e+13;
        k_f_save[124] = 1e-06 * 1.08e+14*exp(-1565.003830652678289*invT);
        k_f_save[125] = 1e-06 * 5.71e+12*exp(+379.92858268256338761*invT);
        k_f_save[126] = 1e-06 * 4e+13;
        k_f_save[127] = 1e-06 * 3e+13;
        k_f_save[128] = 1e-06 * 6e+13;
        k_f_save[129] = 1e-06 * 5e+13;
        k_f_save[130] = 1e-06 * 1.9e+14*exp(-7946.7975863881329133*invT);
        k_f_save[131] = 1e-06 * 9.46e+13*exp(+259.15658288943063781*invT);
        k_f_save[132] = 1e-06 * 5e+13;
        k_f_save[133] = 1e-06 * 5e+12*exp(-754.82499870707954415*invT);
        k_f_save[134] = 1e-06 * 500000*exp(2*tc[0]-3638.256493768123164*invT);
        k_f_save[135] = 1e-06 * 1.6e+15*exp(-6010.4198563715717682*invT);
        k_f_save[136] = 1e-06 * 4e+13;
        k_f_save[137] = 1e-06 * 2.46e+06*exp(2*tc[0]-4161.6018262050320118*invT);
        k_f_save[138] = 1e-06 * 8.1e+11*exp(0.5*tc[0]-2269.5071627792858635*invT);
        k_f_save[139] = 1e-06 * 3e+13;
        k_f_save[140] = 1e-06 * 1.5e+13*exp(-301.92999948283181766*invT);
        k_f_save[141] = 1e-06 * 2.8e+13;
        k_f_save[142] = 1e-06 * 1.2e+13;
        k_f_save[143] = 1e-06 * 7e+13;
        k_f_save[144] = 1e-06 * 4.82e+17*exp(-1.16*tc[0]-576.18308234640414867*invT);
        k_f_save[145] = 1e-06 * 3e+13;
        k_f_save[146] = 1e-06 * 1.2e+13*exp(+286.83349950869023814*invT);
        k_f_save[147] = 1e-06 * 1.6e+13*exp(+286.83349950869023814*invT);
        k_f_save[148] = 1e-06 * 9e+12;
        k_f_save[149] = 1e-06 * 7e+12;
        k_f_save[150] = 1e-06 * 1.4e+13;
        k_f_save[151] = 1e-06 * 4e+13*exp(+276.76916619259583285*invT);
        k_f_save[152] = 1e-06 * 3.56e+13*exp(-15338.0439737278557*invT);
        k_f_save[153] = 1e-06 * 2.31e+12*exp(-10222.846565822879711*invT);
        k_f_save[154] = 1e-06 * 24500*exp(2.47*tc[0]-2606.6623288684481849*invT);
        k_f_save[155] = 1e-06 * 6.77e+16*exp(-1.18*tc[0]-329.10369943628666078*invT);
        k_f_save[156] = 1e-06 * 6.84e+12*exp(0.1*tc[0]-5334.096657530029006*invT);
        k_f_save[157] = 1e-06 * 2.648e+13;
        k_f_save[158] = 1e-06 * 3320*exp(2.81*tc[0]-2948.8496616156576238*invT);
        k_f_save[159] = 1e-06 * 3e+07*exp(1.5*tc[0]-5001.973658098913802*invT);
        k_f_save[160] = 1e-06 * 1e+07*exp(1.5*tc[0]-5001.973658098913802*invT);
        k_f_save[161] = 1e-06 * 227000*exp(2*tc[0]-4629.5933254034207494*invT);
        k_f_save[162] = 1e-06 * 6.14e+06*exp(1.74*tc[0]-5258.6141576593208811*invT);
        k_f_save[163] = 1e-06 * 1.5e+18*exp(-1*tc[0]-8554.6833186802341515*invT);
        k_f_save[164] = 1e-06 * 1.87e+17*exp(-1*tc[0]-8554.6833186802341515*invT);
        k_f_save[165] = 1e-06 * 1.345e+13*exp(-201.28666632188787844*invT);
        k_f_save[166] = 1e-06 * 1.8e+13*exp(-452.89499922424772649*invT);
        k_f_save[167] = 1e-06 * 4.28e-13*exp(7.6*tc[0]+1776.3548302906606295*invT);
        k_f_save[168] = 1e-06 * 1e+13*exp(+379.92858268256338761*invT);
        k_f_save[169] = 1e-06 * 5.68e+10*exp(0.9*tc[0]-1002.9108149488064328*invT);
        k_f_save[170] = 1e-06 * 4.58e+16*exp(-1.39*tc[0]-510.76491579179048586*invT);
        k_f_save[171] = 1 * 8e+12*exp(0.44*tc[0]-43664.110091875525541*invT);
        k_f_save[172] = 1e-06 * 8.4e+11*exp(-1949.9645799932889076*invT);
        k_f_save[173] = 1e-06 * 3.2e+12*exp(-429.7470325972306*invT);
        k_f_save[174] = 1e-06 * 1e+13;
        k_f_save[175] = 1e-06 * 3.37e+13;
        k_f_save[176] = 1e-06 * 6.7e+06*exp(1.83*tc[0]-110.7076664770383303*invT);
        k_f_save[177] = 1e-06 * 1.096e+14;
        k_f_save[178] = 1e-06 * 5e+15*exp(-8720.7448183957912988*invT);
        k_f_save[179] = 1e-06 * 8e+09*exp(0.5*tc[0]+883.1452484872830837*invT);
        k_f_save[180] = 1e-06 * 1.97e+12*exp(0.43*tc[0]+186.1901663477462705*invT);
        k_f_save[181] = 1e-06 * 5.8e+12*exp(-754.82499870707954415*invT);
        k_f_save[182] = 1e-06 * 2.4e+12*exp(-754.82499870707954415*invT);
        k_f_save[183] = 1e-06 * 2e+14*exp(-5529.8479405280650099*invT);
        k_f_save[184] = 1e-06 * 6.82e+10*exp(0.25*tc[0]+470.50758252741292154*invT);
        k_f_save[185] = 1e-06 * 3.03e+11*exp(0.29*tc[0]-5.5353833238519163373*invT);
        k_f_save[186] = 1e-06 * 1.337e+06*exp(1.61*tc[0]+193.23519966901235989*invT);
        k_f_save[187] = 1e-06 * 5.84e+12*exp(-909.81573177493316962*invT);
        k_f_save[188] = 1e-06 * 5.84e+12*exp(-909.81573177493316962*invT);
        k_f_save[189] = 1e-06 * 3.01e+13*exp(-19700.932466254773317*invT);
        k_f_save[190] = 1e-06 * 2.05e+09*exp(1.16*tc[0]-1210.2360812603508293*invT);
        k_f_save[191] = 1e-06 * 2.05e+09*exp(1.16*tc[0]-1210.2360812603508293*invT);
        k_f_save[192] = 1e-06 * 2.343e+10*exp(0.73*tc[0]+560.08014904065305473*invT);
        k_f_save[193] = 1e-06 * 3.01e+12*exp(-5999.8523063896727763*invT);
        k_f_save[194] = 1e-06 * 2.72e+06*exp(1.77*tc[0]-2979.0426615639403281*invT);
        k_f_save[195] = 1e-06 * 4.865e+11*exp(0.422*tc[0]+883.1452484872830837*invT);
        k_f_save[196] = 1e-06 * 1.5e+14;
        k_f_save[197] = 1e-06 * 1.81e+10;
        k_f_save[198] = 1e-06 * 2.35e+10;
        k_f_save[199] = 1e-06 * 2.2e+13;
        k_f_save[200] = 1e-06 * 1.1e+13;
        k_f_save[201] = 1e-06 * 1.2e+13;
        k_f_save[202] = 1e-06 * 3.01e+13;
        k_f_save[203] = 1e-06 * 9.43e+12;
        k_f_save[204] = 1e-06 * 193000*exp(2.68*tc[0]-1869.9531301303384225*invT);
        k_f_save[205] = 1e-06 * 1.32e+06*exp(2.54*tc[0]-3399.7317941766864351*invT);
        k_f_save[206] = 1e-06 * 3.16e+07*exp(1.8*tc[0]-470.00436586160816432*invT);
        k_f_save[207] = 1e-06 * 378*exp(2.72*tc[0]-754.82499870707954415*invT);
        k_f_save[208] = 1e-06 * 0.903*exp(3.65*tc[0]-3600.0120271669643444*invT);
        k_f_save[209] = 1e-06 * 2.55e+06*exp(1.6*tc[0]-2868.3349950869019267*invT);
        k_f_save[210] = 1e-06 * 9.64e+13;
        k_f_save[211] = 1e-06 * 3.613e+13;
        k_f_save[212] = 1e-06 * 4.06e+06*exp(2.19*tc[0]-447.86283256620055226*invT);
        k_f_save[213] = 1e-06 * 2.41e+13;
        k_f_save[214] = 1e-06 * 2.55e+10*exp(0.255*tc[0]+474.53331585385063818*invT);
        k_f_save[215] = 1e-06 * 2.41e+13;
        k_f_save[216] = 1e-06 * 1.927e+13*exp(-0.32*tc[0]);

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
        Kc_save[36] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));
        Kc_save[37] = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));
        Kc_save[38] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));
        Kc_save[39] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));
        Kc_save[40] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[15]) - (g_RT[0] + g_RT[15]));
        Kc_save[41] = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));
        Kc_save[42] = exp((g_RT[1] + g_RT[6]) - (g_RT[2] + g_RT[5]));
        Kc_save[43] = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));
        Kc_save[44] = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));
        Kc_save[45] = exp((g_RT[1] + g_RT[7]) - (g_RT[6] + g_RT[0]));
        Kc_save[46] = exp((g_RT[1] + g_RT[7]) - (g_RT[4] + g_RT[5]));
        Kc_save[47] = exp((g_RT[1] + g_RT[9]) - (g_RT[8] + g_RT[0]));
        Kc_save[48] = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[12]));
        Kc_save[49] = exp((g_RT[1] + g_RT[11]) - (g_RT[9] + g_RT[0]));
        Kc_save[50] = 1.0 / (refC) * exp((g_RT[1] + g_RT[12]) - (g_RT[13]));
        Kc_save[51] = exp((g_RT[1] + g_RT[13]) - (g_RT[12] + g_RT[0]));
        Kc_save[52] = 1.0 / (refC) * exp((g_RT[1] + g_RT[16]) - (g_RT[17]));
        Kc_save[53] = exp((g_RT[1] + g_RT[16]) - (g_RT[0] + g_RT[14]));
        Kc_save[54] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[18]));
        Kc_save[55] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[19]));
        Kc_save[56] = exp((g_RT[1] + g_RT[17]) - (g_RT[16] + g_RT[0]));
        Kc_save[57] = 1.0 / (refC) * exp((g_RT[1] + g_RT[18]) - (g_RT[20]));
        Kc_save[58] = exp((g_RT[1] + g_RT[18]) - (g_RT[0] + g_RT[17]));
        Kc_save[59] = exp((g_RT[1] + g_RT[18]) - (g_RT[4] + g_RT[12]));
        Kc_save[60] = exp((g_RT[1] + g_RT[18]) - (g_RT[11] + g_RT[5]));
        Kc_save[61] = 1.0 / (refC) * exp((g_RT[1] + g_RT[19]) - (g_RT[20]));
        Kc_save[62] = exp((g_RT[1] + g_RT[19]) - (g_RT[1] + g_RT[18]));
        Kc_save[63] = exp((g_RT[1] + g_RT[19]) - (g_RT[0] + g_RT[17]));
        Kc_save[64] = exp((g_RT[1] + g_RT[19]) - (g_RT[4] + g_RT[12]));
        Kc_save[65] = exp((g_RT[1] + g_RT[19]) - (g_RT[11] + g_RT[5]));
        Kc_save[66] = exp((g_RT[1] + g_RT[20]) - (g_RT[18] + g_RT[0]));
        Kc_save[67] = exp((g_RT[1] + g_RT[20]) - (g_RT[19] + g_RT[0]));
        Kc_save[68] = 1.0 / (refC) * exp((g_RT[1] + g_RT[21]) - (g_RT[22]));
        Kc_save[69] = 1.0 / (refC) * exp((g_RT[1] + g_RT[22]) - (g_RT[23]));
        Kc_save[70] = 1.0 / (refC) * exp((g_RT[1] + g_RT[23]) - (g_RT[24]));
        Kc_save[71] = exp((g_RT[1] + g_RT[23]) - (g_RT[0] + g_RT[22]));
        Kc_save[72] = 1.0 / (refC) * exp((g_RT[1] + g_RT[24]) - (g_RT[25]));
        Kc_save[73] = exp((g_RT[1] + g_RT[24]) - (g_RT[23] + g_RT[0]));
        Kc_save[74] = 1.0 / (refC) * exp((g_RT[1] + g_RT[25]) - (g_RT[26]));
        Kc_save[75] = exp((g_RT[1] + g_RT[25]) - (g_RT[0] + g_RT[24]));
        Kc_save[76] = exp((g_RT[1] + g_RT[26]) - (g_RT[25] + g_RT[0]));
        Kc_save[77] = exp((g_RT[1] + g_RT[27]) - (g_RT[11] + g_RT[14]));
        Kc_save[78] = exp((g_RT[1] + g_RT[28]) - (g_RT[27] + g_RT[0]));
        Kc_save[79] = exp((g_RT[1] + g_RT[28]) - (g_RT[12] + g_RT[14]));
        Kc_save[80] = exp((g_RT[1] + g_RT[29]) - (g_RT[1] + g_RT[28]));
        Kc_save[81] = 1.0 / (refC) * exp((g_RT[0] + g_RT[14]) - (g_RT[17]));
        Kc_save[82] = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));
        Kc_save[83] = 1.0 / (refC) * exp((2 * g_RT[4]) - (g_RT[7]));
        Kc_save[84] = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));
        Kc_save[85] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
        Kc_save[86] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
        Kc_save[87] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
        Kc_save[88] = exp((g_RT[4] + g_RT[8]) - (g_RT[1] + g_RT[14]));
        Kc_save[89] = exp((g_RT[4] + g_RT[9]) - (g_RT[1] + g_RT[16]));
        Kc_save[90] = exp((g_RT[4] + g_RT[10]) - (g_RT[1] + g_RT[17]));
        Kc_save[91] = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[5]));
        Kc_save[92] = exp((g_RT[4] + g_RT[11]) - (g_RT[1] + g_RT[17]));
        Kc_save[93] = 1.0 / (refC) * exp((g_RT[4] + g_RT[12]) - (g_RT[20]));
        Kc_save[94] = exp((g_RT[4] + g_RT[12]) - (g_RT[10] + g_RT[5]));
        Kc_save[95] = exp((g_RT[4] + g_RT[12]) - (g_RT[11] + g_RT[5]));
        Kc_save[96] = exp((g_RT[4] + g_RT[13]) - (g_RT[12] + g_RT[5]));
        Kc_save[97] = exp((g_RT[4] + g_RT[14]) - (g_RT[1] + g_RT[15]));
        Kc_save[98] = exp((g_RT[4] + g_RT[16]) - (g_RT[5] + g_RT[14]));
        Kc_save[99] = exp((g_RT[4] + g_RT[17]) - (g_RT[16] + g_RT[5]));
        Kc_save[100] = exp((g_RT[4] + g_RT[18]) - (g_RT[5] + g_RT[17]));
        Kc_save[101] = exp((g_RT[4] + g_RT[19]) - (g_RT[5] + g_RT[17]));
        Kc_save[102] = exp((g_RT[4] + g_RT[20]) - (g_RT[18] + g_RT[5]));
        Kc_save[103] = exp((g_RT[4] + g_RT[20]) - (g_RT[19] + g_RT[5]));
        Kc_save[104] = exp((g_RT[4] + g_RT[21]) - (g_RT[1] + g_RT[27]));
        Kc_save[105] = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[28]));
        Kc_save[106] = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[29]));
        Kc_save[107] = exp((g_RT[4] + g_RT[22]) - (g_RT[21] + g_RT[5]));
        Kc_save[108] = exp((g_RT[4] + g_RT[22]) - (g_RT[12] + g_RT[14]));
        Kc_save[109] = exp((g_RT[4] + g_RT[23]) - (g_RT[5] + g_RT[22]));
        Kc_save[110] = exp((g_RT[4] + g_RT[24]) - (g_RT[23] + g_RT[5]));
        Kc_save[111] = exp((g_RT[4] + g_RT[26]) - (g_RT[25] + g_RT[5]));
        Kc_save[112] = exp((g_RT[4] + g_RT[28]) - (g_RT[27] + g_RT[5]));
        Kc_save[113] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
        Kc_save[114] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
        Kc_save[115] = exp((g_RT[6] + g_RT[10]) - (g_RT[4] + g_RT[17]));
        Kc_save[116] = exp((g_RT[6] + g_RT[12]) - (g_RT[3] + g_RT[13]));
        Kc_save[117] = exp((g_RT[6] + g_RT[12]) - (g_RT[4] + g_RT[19]));
        Kc_save[118] = exp((g_RT[6] + g_RT[14]) - (g_RT[4] + g_RT[15]));
        Kc_save[119] = exp((g_RT[6] + g_RT[17]) - (g_RT[16] + g_RT[7]));
        Kc_save[120] = exp((g_RT[8] + g_RT[3]) - (g_RT[2] + g_RT[14]));
        Kc_save[121] = exp((g_RT[8] + g_RT[10]) - (g_RT[1] + g_RT[21]));
        Kc_save[122] = exp((g_RT[8] + g_RT[12]) - (g_RT[1] + g_RT[22]));
        Kc_save[123] = exp((g_RT[9] + g_RT[3]) - (g_RT[2] + g_RT[16]));
        Kc_save[124] = exp((g_RT[9] + g_RT[0]) - (g_RT[1] + g_RT[10]));
        Kc_save[125] = exp((g_RT[9] + g_RT[5]) - (g_RT[1] + g_RT[17]));
        Kc_save[126] = exp((g_RT[9] + g_RT[10]) - (g_RT[1] + g_RT[22]));
        Kc_save[127] = exp((g_RT[9] + g_RT[12]) - (g_RT[1] + g_RT[23]));
        Kc_save[128] = exp((g_RT[9] + g_RT[13]) - (g_RT[1] + g_RT[24]));
        Kc_save[129] = 1.0 / (refC) * exp((g_RT[9] + g_RT[14]) - (g_RT[27]));
        Kc_save[130] = exp((g_RT[9] + g_RT[15]) - (g_RT[16] + g_RT[14]));
        Kc_save[131] = exp((g_RT[9] + g_RT[17]) - (g_RT[1] + g_RT[28]));
        Kc_save[132] = exp((g_RT[9] + g_RT[27]) - (g_RT[14] + g_RT[22]));
        Kc_save[133] = refC * exp((g_RT[10] + g_RT[3]) - (g_RT[4] + g_RT[1] + g_RT[14]));
        Kc_save[134] = exp((g_RT[10] + g_RT[0]) - (g_RT[1] + g_RT[12]));
        Kc_save[135] = exp((2 * g_RT[10]) - (g_RT[0] + g_RT[22]));
        Kc_save[136] = exp((g_RT[10] + g_RT[12]) - (g_RT[1] + g_RT[24]));
        Kc_save[137] = exp((g_RT[10] + g_RT[13]) - (2 * g_RT[12]));
        Kc_save[138] = 1.0 / (refC) * exp((g_RT[10] + g_RT[14]) - (g_RT[28]));
        Kc_save[139] = exp((g_RT[10] + g_RT[27]) - (g_RT[23] + g_RT[14]));
        Kc_save[140] = exp((g_RT[11] + g_RT[30]) - (g_RT[10] + g_RT[30]));
        Kc_save[141] = refC * exp((g_RT[11] + g_RT[3]) - (g_RT[1] + g_RT[4] + g_RT[14]));
        Kc_save[142] = exp((g_RT[11] + g_RT[3]) - (g_RT[14] + g_RT[5]));
        Kc_save[143] = exp((g_RT[11] + g_RT[0]) - (g_RT[12] + g_RT[1]));
        Kc_save[144] = 1.0 / (refC) * exp((g_RT[11] + g_RT[5]) - (g_RT[20]));
        Kc_save[145] = exp((g_RT[11] + g_RT[5]) - (g_RT[10] + g_RT[5]));
        Kc_save[146] = exp((g_RT[11] + g_RT[12]) - (g_RT[1] + g_RT[24]));
        Kc_save[147] = exp((g_RT[11] + g_RT[13]) - (2 * g_RT[12]));
        Kc_save[148] = exp((g_RT[11] + g_RT[14]) - (g_RT[10] + g_RT[14]));
        Kc_save[149] = exp((g_RT[11] + g_RT[15]) - (g_RT[10] + g_RT[15]));
        Kc_save[150] = exp((g_RT[11] + g_RT[15]) - (g_RT[14] + g_RT[17]));
        Kc_save[151] = exp((g_RT[11] + g_RT[26]) - (g_RT[12] + g_RT[25]));
        Kc_save[152] = exp((g_RT[12] + g_RT[3]) - (g_RT[2] + g_RT[19]));
        Kc_save[153] = exp((g_RT[12] + g_RT[3]) - (g_RT[4] + g_RT[17]));
        Kc_save[154] = exp((g_RT[12] + g_RT[7]) - (g_RT[6] + g_RT[13]));
        Kc_save[155] = 1.0 / (refC) * exp((2 * g_RT[12]) - (g_RT[26]));
        Kc_save[156] = exp((2 * g_RT[12]) - (g_RT[1] + g_RT[25]));
        Kc_save[157] = exp((g_RT[12] + g_RT[16]) - (g_RT[13] + g_RT[14]));
        Kc_save[158] = exp((g_RT[12] + g_RT[17]) - (g_RT[16] + g_RT[13]));
        Kc_save[159] = exp((g_RT[12] + g_RT[20]) - (g_RT[18] + g_RT[13]));
        Kc_save[160] = exp((g_RT[12] + g_RT[20]) - (g_RT[19] + g_RT[13]));
        Kc_save[161] = exp((g_RT[12] + g_RT[24]) - (g_RT[23] + g_RT[13]));
        Kc_save[162] = exp((g_RT[12] + g_RT[26]) - (g_RT[25] + g_RT[13]));
        Kc_save[163] = refC * exp((g_RT[16] + g_RT[5]) - (g_RT[1] + g_RT[14] + g_RT[5]));
        Kc_save[164] = refC * exp((g_RT[16]) - (g_RT[1] + g_RT[14]));
        Kc_save[165] = exp((g_RT[16] + g_RT[3]) - (g_RT[6] + g_RT[14]));
        Kc_save[166] = exp((g_RT[18] + g_RT[3]) - (g_RT[6] + g_RT[17]));
        Kc_save[167] = exp((g_RT[19] + g_RT[3]) - (g_RT[6] + g_RT[17]));
        Kc_save[168] = exp((g_RT[21] + g_RT[3]) - (g_RT[16] + g_RT[14]));
        Kc_save[169] = exp((g_RT[21] + g_RT[0]) - (g_RT[1] + g_RT[22]));
        Kc_save[170] = exp((g_RT[23] + g_RT[3]) - (g_RT[16] + g_RT[17]));
        Kc_save[171] = refC * exp((g_RT[24]) - (g_RT[0] + g_RT[22]));
        Kc_save[172] = exp((g_RT[25] + g_RT[3]) - (g_RT[6] + g_RT[24]));
        Kc_save[173] = refC * exp((g_RT[27] + g_RT[3]) - (g_RT[4] + 2 * g_RT[14]));
        Kc_save[174] = refC * exp((2 * g_RT[27]) - (2 * g_RT[14] + g_RT[22]));
        Kc_save[175] = refC * exp((g_RT[2] + g_RT[12]) - (g_RT[1] + g_RT[0] + g_RT[14]));
        Kc_save[176] = exp((g_RT[2] + g_RT[24]) - (g_RT[1] + g_RT[33]));
        Kc_save[177] = exp((g_RT[2] + g_RT[25]) - (g_RT[1] + g_RT[34]));
        Kc_save[178] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
        Kc_save[179] = exp((g_RT[4] + g_RT[12]) - (g_RT[0] + g_RT[17]));
        Kc_save[180] = 1.0 / (refC) * exp((g_RT[9] + g_RT[0]) - (g_RT[12]));
        Kc_save[181] = refC * exp((g_RT[10] + g_RT[3]) - (2 * g_RT[1] + g_RT[15]));
        Kc_save[182] = exp((g_RT[10] + g_RT[3]) - (g_RT[2] + g_RT[17]));
        Kc_save[183] = refC * exp((g_RT[10] + g_RT[10]) - (2 * g_RT[1] + g_RT[22]));
        Kc_save[184] = exp((g_RT[11] + g_RT[5]) - (g_RT[0] + g_RT[17]));
        Kc_save[185] = exp((g_RT[23] + g_RT[3]) - (g_RT[2] + g_RT[33]));
        Kc_save[186] = exp((g_RT[23] + g_RT[3]) - (g_RT[6] + g_RT[22]));
        Kc_save[187] = exp((g_RT[2] + g_RT[34]) - (g_RT[4] + g_RT[33]));
        Kc_save[188] = refC * exp((g_RT[2] + g_RT[34]) - (g_RT[4] + g_RT[12] + g_RT[14]));
        Kc_save[189] = refC * exp((g_RT[3] + g_RT[34]) - (g_RT[6] + g_RT[12] + g_RT[14]));
        Kc_save[190] = exp((g_RT[1] + g_RT[34]) - (g_RT[33] + g_RT[0]));
        Kc_save[191] = refC * exp((g_RT[1] + g_RT[34]) - (g_RT[12] + g_RT[0] + g_RT[14]));
        Kc_save[192] = refC * exp((g_RT[4] + g_RT[34]) - (g_RT[12] + g_RT[5] + g_RT[14]));
        Kc_save[193] = refC * exp((g_RT[6] + g_RT[34]) - (g_RT[12] + g_RT[7] + g_RT[14]));
        Kc_save[194] = refC * exp((g_RT[12] + g_RT[34]) - (g_RT[12] + g_RT[13] + g_RT[14]));
        Kc_save[195] = 1.0 / (refC) * exp((g_RT[1] + g_RT[28]) - (g_RT[33]));
        Kc_save[196] = refC * exp((g_RT[2] + g_RT[33]) - (g_RT[1] + g_RT[10] + g_RT[15]));
        Kc_save[197] = refC * exp((g_RT[3] + g_RT[33]) - (g_RT[4] + g_RT[14] + g_RT[17]));
        Kc_save[198] = refC * exp((g_RT[3] + g_RT[33]) - (g_RT[4] + 2 * g_RT[16]));
        Kc_save[199] = exp((g_RT[1] + g_RT[33]) - (g_RT[12] + g_RT[16]));
        Kc_save[200] = exp((g_RT[1] + g_RT[33]) - (g_RT[28] + g_RT[0]));
        Kc_save[201] = exp((g_RT[4] + g_RT[33]) - (g_RT[5] + g_RT[28]));
        Kc_save[202] = exp((g_RT[4] + g_RT[33]) - (g_RT[16] + g_RT[18]));
        Kc_save[203] = 1.0 / (refC) * exp((g_RT[12] + g_RT[25]) - (g_RT[32]));
        Kc_save[204] = exp((g_RT[2] + g_RT[32]) - (g_RT[4] + g_RT[31]));
        Kc_save[205] = exp((g_RT[1] + g_RT[32]) - (g_RT[31] + g_RT[0]));
        Kc_save[206] = exp((g_RT[4] + g_RT[32]) - (g_RT[31] + g_RT[5]));
        Kc_save[207] = exp((g_RT[31] + g_RT[7]) - (g_RT[6] + g_RT[32]));
        Kc_save[208] = exp((g_RT[12] + g_RT[32]) - (g_RT[31] + g_RT[13]));
        Kc_save[209] = 1.0 / (refC) * exp((g_RT[12] + g_RT[24]) - (g_RT[31]));
        Kc_save[210] = exp((g_RT[2] + g_RT[31]) - (g_RT[25] + g_RT[17]));
        Kc_save[211] = 1.0 / (refC) * exp((g_RT[1] + g_RT[31]) - (g_RT[32]));
        Kc_save[212] = exp((g_RT[1] + g_RT[31]) - (g_RT[12] + g_RT[25]));
        Kc_save[213] = exp((g_RT[4] + g_RT[31]) - (g_RT[25] + g_RT[18]));
        Kc_save[214] = exp((g_RT[6] + g_RT[31]) - (g_RT[3] + g_RT[32]));
        Kc_save[215] = refC * exp((g_RT[6] + g_RT[31]) - (g_RT[4] + g_RT[25] + g_RT[17]));
        Kc_save[216] = exp((g_RT[12] + g_RT[31]) - (2 * g_RT[25]));
    }

    /*reaction 1: 2 O + M <=> O2 + M */
    phi_f = sc[2]*sc[2];
    alpha = mixture + 1.4*sc[0] + 14.4*sc[5] + sc[13] + 0.75*sc[14] + 2.6*sc[15] + 2*sc[26];
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
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
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

    /*reaction 12: O + CO (+M) <=> CO2 (+M) */
    phi_f = sc[2]*sc[14];
    alpha = mixture + sc[0] + 5*sc[3] + 5*sc[5] + sc[13] + 0.5*sc[14] + 2.5*sc[15] + 2*sc[26];
    k_f = k_f_save[11];
    redP = 1e-12 * alpha / k_f * 6.02e+14*exp(-1509.6499974141590883*invT);
    F = redP / (1 + redP);
    k_f *= F;
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
    alpha = mixture + -1*sc[3] + -1*sc[5] + -0.25*sc[14] + 0.5*sc[15] + 0.5*sc[26] + -1*sc[30];
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

    /*reaction 37: H + O2 <=> O + OH */
    phi_f = sc[1]*sc[3];
    k_f = k_f_save[36];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = Kc_save[36];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 38: 2 H + M <=> H2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[0] + -1*sc[5] + sc[13] + -1*sc[15] + 2*sc[26];
    k_f = alpha * k_f_save[37];
    q_f = phi_f * k_f;
    phi_r = sc[0];
    Kc = Kc_save[37];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 39: 2 H + H2 <=> 2 H2 */
    phi_f = sc[1]*sc[1]*sc[0];
    k_f = k_f_save[38];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[0];
    Kc = Kc_save[38];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[0] += 2 * qdot;

    /*reaction 40: 2 H + H2O <=> H2 + H2O */
    phi_f = sc[1]*sc[1]*sc[5];
    k_f = k_f_save[39];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[5];
    Kc = Kc_save[39];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 41: 2 H + CO2 <=> H2 + CO2 */
    phi_f = sc[1]*sc[1]*sc[15];
    k_f = k_f_save[40];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[15];
    Kc = Kc_save[40];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 42: H + OH + M <=> H2O + M */
    phi_f = sc[1]*sc[4];
    alpha = mixture + -0.27*sc[0] + 2.65*sc[5] + sc[13] + 2*sc[26];
    k_f = alpha * k_f_save[41];
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = Kc_save[41];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 43: H + HO2 <=> O + H2O */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[42];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = Kc_save[42];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 44: H + HO2 <=> O2 + H2 */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[43];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = Kc_save[43];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 45: H + HO2 <=> 2 OH */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[44];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = Kc_save[44];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[4] += 2 * qdot;

    /*reaction 46: H + H2O2 <=> HO2 + H2 */
    phi_f = sc[1]*sc[7];
    k_f = k_f_save[45];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = Kc_save[45];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 47: H + H2O2 <=> OH + H2O */
    phi_f = sc[1]*sc[7];
    k_f = k_f_save[46];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[5];
    Kc = Kc_save[46];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 48: H + CH <=> C + H2 */
    phi_f = sc[1]*sc[9];
    k_f = k_f_save[47];
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = Kc_save[47];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[9] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 49: H + CH2 (+M) <=> CH3 (+M) */
    phi_f = sc[1]*sc[10];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[48];
    redP = 1e-12 * alpha / k_f * 1.04e+26*exp(-2.76*tc[0]-805.14666528755151376*invT);
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
    Kc = Kc_save[48];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 50: H + CH2(S) <=> CH + H2 */
    phi_f = sc[1]*sc[11];
    k_f = k_f_save[49];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = Kc_save[49];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 51: H + CH3 (+M) <=> CH4 (+M) */
    phi_f = sc[1]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + 2*sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[50];
    redP = 1e-12 * alpha / k_f * 2.62e+33*exp(-4.76*tc[0]-1227.8486645635159675*invT);
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
    Kc = Kc_save[50];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 52: H + CH4 <=> CH3 + H2 */
    phi_f = sc[1]*sc[13];
    k_f = k_f_save[51];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[0];
    Kc = Kc_save[51];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 53: H + HCO (+M) <=> CH2O (+M) */
    phi_f = sc[1]*sc[16];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[52];
    redP = 1e-12 * alpha / k_f * 2.47e+24*exp(-2.57*tc[0]-213.86708296700587084*invT);
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
    Kc = Kc_save[52];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[16] -= 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 54: H + HCO <=> H2 + CO */
    phi_f = sc[1]*sc[16];
    k_f = k_f_save[53];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[14];
    Kc = Kc_save[53];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[16] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 55: H + CH2O (+M) <=> CH2OH (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[54];
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
    Kc = Kc_save[54];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[18] += 1 * qdot;

    /*reaction 56: H + CH2O (+M) <=> CH3O (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[55];
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
    Kc = Kc_save[55];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[19] += 1 * qdot;

    /*reaction 57: H + CH2O <=> HCO + H2 */
    phi_f = sc[1]*sc[17];
    k_f = k_f_save[56];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[0];
    Kc = Kc_save[56];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 58: H + CH2OH (+M) <=> CH3OH (+M) */
    phi_f = sc[1]*sc[18];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[57];
    redP = 1e-12 * alpha / k_f * 4.36e+31*exp(-4.65*tc[0]-2556.3406622879761017*invT);
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
    Kc = Kc_save[57];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[20] += 1 * qdot;

    /*reaction 59: H + CH2OH <=> H2 + CH2O */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[58];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[17];
    Kc = Kc_save[58];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 60: H + CH2OH <=> OH + CH3 */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[59];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = Kc_save[59];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 61: H + CH2OH <=> CH2(S) + H2O */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[60];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[60];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 62: H + CH3O (+M) <=> CH3OH (+M) */
    phi_f = sc[1]*sc[19];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[61];
    redP = 1e-12 * alpha / k_f * 4.66e+41*exp(-7.44*tc[0]-7085.2906545304531392*invT);
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
    Kc = Kc_save[61];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[20] += 1 * qdot;

    /*reaction 63: H + CH3O <=> H + CH2OH */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[62];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[18];
    Kc = Kc_save[62];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[18] += 1 * qdot;

    /*reaction 64: H + CH3O <=> H2 + CH2O */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[63];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[17];
    Kc = Kc_save[63];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 65: H + CH3O <=> OH + CH3 */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[64];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = Kc_save[64];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 66: H + CH3O <=> CH2(S) + H2O */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[65];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[65];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 67: H + CH3OH <=> CH2OH + H2 */
    phi_f = sc[1]*sc[20];
    k_f = k_f_save[66];
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[0];
    Kc = Kc_save[66];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 68: H + CH3OH <=> CH3O + H2 */
    phi_f = sc[1]*sc[20];
    k_f = k_f_save[67];
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[0];
    Kc = Kc_save[67];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[19] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 69: H + C2H (+M) <=> C2H2 (+M) */
    phi_f = sc[1]*sc[21];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[68];
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
    Kc = Kc_save[68];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[21] -= 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 70: H + C2H2 (+M) <=> C2H3 (+M) */
    phi_f = sc[1]*sc[22];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[69];
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
    Kc = Kc_save[69];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[23] += 1 * qdot;

    /*reaction 71: H + C2H3 (+M) <=> C2H4 (+M) */
    phi_f = sc[1]*sc[23];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[70];
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
    Kc = Kc_save[70];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[23] -= 1 * qdot;
    wdot[24] += 1 * qdot;

    /*reaction 72: H + C2H3 <=> H2 + C2H2 */
    phi_f = sc[1]*sc[23];
    k_f = k_f_save[71];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = Kc_save[71];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[23] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 73: H + C2H4 (+M) <=> C2H5 (+M) */
    phi_f = sc[1]*sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[72];
    redP = 1e-12 * alpha / k_f * 6e+41*exp(-7.62*tc[0]-3507.4201606588962932*invT);
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
    Kc = Kc_save[72];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[24] -= 1 * qdot;
    wdot[25] += 1 * qdot;

    /*reaction 74: H + C2H4 <=> C2H3 + H2 */
    phi_f = sc[1]*sc[24];
    k_f = k_f_save[73];
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[0];
    Kc = Kc_save[73];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[24] -= 1 * qdot;
    wdot[23] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 75: H + C2H5 (+M) <=> C2H6 (+M) */
    phi_f = sc[1]*sc[25];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[74];
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
    Kc = Kc_save[74];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[25] -= 1 * qdot;
    wdot[26] += 1 * qdot;

    /*reaction 76: H + C2H5 <=> H2 + C2H4 */
    phi_f = sc[1]*sc[25];
    k_f = k_f_save[75];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[24];
    Kc = Kc_save[75];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[25] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[24] += 1 * qdot;

    /*reaction 77: H + C2H6 <=> C2H5 + H2 */
    phi_f = sc[1]*sc[26];
    k_f = k_f_save[76];
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[0];
    Kc = Kc_save[76];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[26] -= 1 * qdot;
    wdot[25] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 78: H + HCCO <=> CH2(S) + CO */
    phi_f = sc[1]*sc[27];
    k_f = k_f_save[77];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[14];
    Kc = Kc_save[77];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[27] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 79: H + CH2CO <=> HCCO + H2 */
    phi_f = sc[1]*sc[28];
    k_f = k_f_save[78];
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[0];
    Kc = Kc_save[78];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[28] -= 1 * qdot;
    wdot[27] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 80: H + CH2CO <=> CH3 + CO */
    phi_f = sc[1]*sc[28];
    k_f = k_f_save[79];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[14];
    Kc = Kc_save[79];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[28] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 81: H + HCCOH <=> H + CH2CO */
    phi_f = sc[1]*sc[29];
    k_f = k_f_save[80];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[80];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[29] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[28] += 1 * qdot;

    /*reaction 82: H2 + CO (+M) <=> CH2O (+M) */
    phi_f = sc[0]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[81];
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
    Kc = Kc_save[81];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 83: OH + H2 <=> H + H2O */
    phi_f = sc[4]*sc[0];
    k_f = k_f_save[82];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    Kc = Kc_save[82];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 84: 2 OH (+M) <=> H2O2 (+M) */
    phi_f = sc[4]*sc[4];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[83];
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
    Kc = Kc_save[83];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 2 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 85: 2 OH <=> O + H2O */
    phi_f = sc[4]*sc[4];
    k_f = k_f_save[84];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = Kc_save[84];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 2 * qdot;
    wdot[2] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 86: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = k_f_save[85];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = Kc_save[85];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 87: OH + H2O2 <=> HO2 + H2O */
    phi_f = sc[4]*sc[7];
    k_f = k_f_save[86];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = Kc_save[86];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[6] += 1 * qdot;
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

    /*reaction 89: OH + C <=> H + CO */
    phi_f = sc[4]*sc[8];
    k_f = k_f_save[88];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[88];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 90: OH + CH <=> H + HCO */
    phi_f = sc[4]*sc[9];
    k_f = k_f_save[89];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_save[89];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[9] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 91: OH + CH2 <=> H + CH2O */
    phi_f = sc[4]*sc[10];
    k_f = k_f_save[90];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[90];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 92: OH + CH2 <=> CH + H2O */
    phi_f = sc[4]*sc[10];
    k_f = k_f_save[91];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[5];
    Kc = Kc_save[91];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 93: OH + CH2(S) <=> H + CH2O */
    phi_f = sc[4]*sc[11];
    k_f = k_f_save[92];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[92];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[11] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 94: OH + CH3 (+M) <=> CH3OH (+M) */
    phi_f = sc[4]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[93];
    redP = 1e-12 * alpha / k_f * 4e+36*exp(-5.92*tc[0]-1580.1003306268198685*invT);
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
    Kc = Kc_save[93];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[20] += 1 * qdot;

    /*reaction 95: OH + CH3 <=> CH2 + H2O */
    phi_f = sc[4]*sc[12];
    k_f = k_f_save[94];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[5];
    Kc = Kc_save[94];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 96: OH + CH3 <=> CH2(S) + H2O */
    phi_f = sc[4]*sc[12];
    k_f = k_f_save[95];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[95];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[11] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 97: OH + CH4 <=> CH3 + H2O */
    phi_f = sc[4]*sc[13];
    k_f = k_f_save[96];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[5];
    Kc = Kc_save[96];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 98: OH + CO <=> H + CO2 */
    phi_f = sc[4]*sc[14];
    k_f = k_f_save[97];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[15];
    Kc = Kc_save[97];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 99: OH + HCO <=> H2O + CO */
    phi_f = sc[4]*sc[16];
    k_f = k_f_save[98];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[14];
    Kc = Kc_save[98];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[16] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 100: OH + CH2O <=> HCO + H2O */
    phi_f = sc[4]*sc[17];
    k_f = k_f_save[99];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[5];
    Kc = Kc_save[99];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 101: OH + CH2OH <=> H2O + CH2O */
    phi_f = sc[4]*sc[18];
    k_f = k_f_save[100];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[17];
    Kc = Kc_save[100];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[18] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 102: OH + CH3O <=> H2O + CH2O */
    phi_f = sc[4]*sc[19];
    k_f = k_f_save[101];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[17];
    Kc = Kc_save[101];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[19] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 103: OH + CH3OH <=> CH2OH + H2O */
    phi_f = sc[4]*sc[20];
    k_f = k_f_save[102];
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[5];
    Kc = Kc_save[102];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 104: OH + CH3OH <=> CH3O + H2O */
    phi_f = sc[4]*sc[20];
    k_f = k_f_save[103];
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[5];
    Kc = Kc_save[103];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[19] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 105: OH + C2H <=> H + HCCO */
    phi_f = sc[4]*sc[21];
    k_f = k_f_save[104];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[27];
    Kc = Kc_save[104];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[21] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[27] += 1 * qdot;

    /*reaction 106: OH + C2H2 <=> H + CH2CO */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[105];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[105];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[28] += 1 * qdot;

    /*reaction 107: OH + C2H2 <=> H + HCCOH */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[106];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[29];
    Kc = Kc_save[106];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[29] += 1 * qdot;

    /*reaction 108: OH + C2H2 <=> C2H + H2O */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[107];
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[5];
    Kc = Kc_save[107];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[21] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 109: OH + C2H2 <=> CH3 + CO */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[108];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[14];
    Kc = Kc_save[108];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[22] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 110: OH + C2H3 <=> H2O + C2H2 */
    phi_f = sc[4]*sc[23];
    k_f = k_f_save[109];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[22];
    Kc = Kc_save[109];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[23] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 111: OH + C2H4 <=> C2H3 + H2O */
    phi_f = sc[4]*sc[24];
    k_f = k_f_save[110];
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[5];
    Kc = Kc_save[110];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[24] -= 1 * qdot;
    wdot[23] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 112: OH + C2H6 <=> C2H5 + H2O */
    phi_f = sc[4]*sc[26];
    k_f = k_f_save[111];
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[5];
    Kc = Kc_save[111];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[26] -= 1 * qdot;
    wdot[25] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 113: OH + CH2CO <=> HCCO + H2O */
    phi_f = sc[4]*sc[28];
    k_f = k_f_save[112];
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[5];
    Kc = Kc_save[112];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[28] -= 1 * qdot;
    wdot[27] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 114: 2 HO2 <=> O2 + H2O2 */
    phi_f = sc[6]*sc[6];
    k_f = k_f_save[113];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = Kc_save[113];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 2 * qdot;
    wdot[3] += 1 * qdot;
    wdot[7] += 1 * qdot;

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

    /*reaction 116: HO2 + CH2 <=> OH + CH2O */
    phi_f = sc[6]*sc[10];
    k_f = k_f_save[115];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[115];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 117: HO2 + CH3 <=> O2 + CH4 */
    phi_f = sc[6]*sc[12];
    k_f = k_f_save[116];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[13];
    Kc = Kc_save[116];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 118: HO2 + CH3 <=> OH + CH3O */
    phi_f = sc[6]*sc[12];
    k_f = k_f_save[117];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[19];
    Kc = Kc_save[117];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[19] += 1 * qdot;

    /*reaction 119: HO2 + CO <=> OH + CO2 */
    phi_f = sc[6]*sc[14];
    k_f = k_f_save[118];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[15];
    Kc = Kc_save[118];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 120: HO2 + CH2O <=> HCO + H2O2 */
    phi_f = sc[6]*sc[17];
    k_f = k_f_save[119];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[7];
    Kc = Kc_save[119];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[7] += 1 * qdot;

    /*reaction 121: C + O2 <=> O + CO */
    phi_f = sc[8]*sc[3];
    k_f = k_f_save[120];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[14];
    Kc = Kc_save[120];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 122: C + CH2 <=> H + C2H */
    phi_f = sc[8]*sc[10];
    k_f = k_f_save[121];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[21];
    Kc = Kc_save[121];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[21] += 1 * qdot;

    /*reaction 123: C + CH3 <=> H + C2H2 */
    phi_f = sc[8]*sc[12];
    k_f = k_f_save[122];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = Kc_save[122];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 124: CH + O2 <=> O + HCO */
    phi_f = sc[9]*sc[3];
    k_f = k_f_save[123];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[16];
    Kc = Kc_save[123];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 125: CH + H2 <=> H + CH2 */
    phi_f = sc[9]*sc[0];
    k_f = k_f_save[124];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[10];
    Kc = Kc_save[124];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[10] += 1 * qdot;

    /*reaction 126: CH + H2O <=> H + CH2O */
    phi_f = sc[9]*sc[5];
    k_f = k_f_save[125];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[125];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 127: CH + CH2 <=> H + C2H2 */
    phi_f = sc[9]*sc[10];
    k_f = k_f_save[126];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = Kc_save[126];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 128: CH + CH3 <=> H + C2H3 */
    phi_f = sc[9]*sc[12];
    k_f = k_f_save[127];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[23];
    Kc = Kc_save[127];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[23] += 1 * qdot;

    /*reaction 129: CH + CH4 <=> H + C2H4 */
    phi_f = sc[9]*sc[13];
    k_f = k_f_save[128];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = Kc_save[128];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[24] += 1 * qdot;

    /*reaction 130: CH + CO (+M) <=> HCCO (+M) */
    phi_f = sc[9]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[129];
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
    Kc = Kc_save[129];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[27] += 1 * qdot;

    /*reaction 131: CH + CO2 <=> HCO + CO */
    phi_f = sc[9]*sc[15];
    k_f = k_f_save[130];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[14];
    Kc = Kc_save[130];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 132: CH + CH2O <=> H + CH2CO */
    phi_f = sc[9]*sc[17];
    k_f = k_f_save[131];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[131];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[28] += 1 * qdot;

    /*reaction 133: CH + HCCO <=> CO + C2H2 */
    phi_f = sc[9]*sc[27];
    k_f = k_f_save[132];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[22];
    Kc = Kc_save[132];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[27] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 134: CH2 + O2 => OH + H + CO */
    phi_f = sc[10]*sc[3];
    k_f = k_f_save[133];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 135: CH2 + H2 <=> H + CH3 */
    phi_f = sc[10]*sc[0];
    k_f = k_f_save[134];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = Kc_save[134];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 136: 2 CH2 <=> H2 + C2H2 */
    phi_f = sc[10]*sc[10];
    k_f = k_f_save[135];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = Kc_save[135];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 2 * qdot;
    wdot[0] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 137: CH2 + CH3 <=> H + C2H4 */
    phi_f = sc[10]*sc[12];
    k_f = k_f_save[136];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = Kc_save[136];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[24] += 1 * qdot;

    /*reaction 138: CH2 + CH4 <=> 2 CH3 */
    phi_f = sc[10]*sc[13];
    k_f = k_f_save[137];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[12];
    Kc = Kc_save[137];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[12] += 2 * qdot;

    /*reaction 139: CH2 + CO (+M) <=> CH2CO (+M) */
    phi_f = sc[10]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[138];
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
    Kc = Kc_save[138];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[28] += 1 * qdot;

    /*reaction 140: CH2 + HCCO <=> C2H3 + CO */
    phi_f = sc[10]*sc[27];
    k_f = k_f_save[139];
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[14];
    Kc = Kc_save[139];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[27] -= 1 * qdot;
    wdot[23] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 141: CH2(S) + N2 <=> CH2 + N2 */
    phi_f = sc[11]*sc[30];
    k_f = k_f_save[140];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[30];
    Kc = Kc_save[140];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[30] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[30] += 1 * qdot;

    /*reaction 142: CH2(S) + O2 <=> H + OH + CO */
    phi_f = sc[11]*sc[3];
    k_f = k_f_save[141];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[4]*sc[14];
    Kc = Kc_save[141];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 143: CH2(S) + O2 <=> CO + H2O */
    phi_f = sc[11]*sc[3];
    k_f = k_f_save[142];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[5];
    Kc = Kc_save[142];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 144: CH2(S) + H2 <=> CH3 + H */
    phi_f = sc[11]*sc[0];
    k_f = k_f_save[143];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[1];
    Kc = Kc_save[143];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 145: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    phi_f = sc[11]*sc[5];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[144];
    redP = 1e-12 * alpha / k_f * 1.88e+38*exp(-6.36*tc[0]-2536.2119956557871774*invT);
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
    Kc = Kc_save[144];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[20] += 1 * qdot;

    /*reaction 146: CH2(S) + H2O <=> CH2 + H2O */
    phi_f = sc[11]*sc[5];
    k_f = k_f_save[145];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[5];
    Kc = Kc_save[145];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 147: CH2(S) + CH3 <=> H + C2H4 */
    phi_f = sc[11]*sc[12];
    k_f = k_f_save[146];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = Kc_save[146];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[24] += 1 * qdot;

    /*reaction 148: CH2(S) + CH4 <=> 2 CH3 */
    phi_f = sc[11]*sc[13];
    k_f = k_f_save[147];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[12];
    Kc = Kc_save[147];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[12] += 2 * qdot;

    /*reaction 149: CH2(S) + CO <=> CH2 + CO */
    phi_f = sc[11]*sc[14];
    k_f = k_f_save[148];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[14];
    Kc = Kc_save[148];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[14] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 150: CH2(S) + CO2 <=> CH2 + CO2 */
    phi_f = sc[11]*sc[15];
    k_f = k_f_save[149];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[15];
    Kc = Kc_save[149];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 151: CH2(S) + CO2 <=> CO + CH2O */
    phi_f = sc[11]*sc[15];
    k_f = k_f_save[150];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[17];
    Kc = Kc_save[150];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[15] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 152: CH2(S) + C2H6 <=> CH3 + C2H5 */
    phi_f = sc[11]*sc[26];
    k_f = k_f_save[151];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[25];
    Kc = Kc_save[151];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[26] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[25] += 1 * qdot;

    /*reaction 153: CH3 + O2 <=> O + CH3O */
    phi_f = sc[12]*sc[3];
    k_f = k_f_save[152];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[19];
    Kc = Kc_save[152];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[19] += 1 * qdot;

    /*reaction 154: CH3 + O2 <=> OH + CH2O */
    phi_f = sc[12]*sc[3];
    k_f = k_f_save[153];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[153];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 155: CH3 + H2O2 <=> HO2 + CH4 */
    phi_f = sc[12]*sc[7];
    k_f = k_f_save[154];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = Kc_save[154];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 156: 2 CH3 (+M) <=> C2H6 (+M) */
    phi_f = sc[12]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[155];
    redP = 1e-12 * alpha / k_f * 3.4e+41*exp(-7.03*tc[0]-1389.8844309526357392*invT);
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
    Kc = Kc_save[155];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 2 * qdot;
    wdot[26] += 1 * qdot;

    /*reaction 157: 2 CH3 <=> H + C2H5 */
    phi_f = sc[12]*sc[12];
    k_f = k_f_save[156];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[25];
    Kc = Kc_save[156];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 2 * qdot;
    wdot[1] += 1 * qdot;
    wdot[25] += 1 * qdot;

    /*reaction 158: CH3 + HCO <=> CH4 + CO */
    phi_f = sc[12]*sc[16];
    k_f = k_f_save[157];
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[14];
    Kc = Kc_save[157];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[16] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 159: CH3 + CH2O <=> HCO + CH4 */
    phi_f = sc[12]*sc[17];
    k_f = k_f_save[158];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[13];
    Kc = Kc_save[158];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 160: CH3 + CH3OH <=> CH2OH + CH4 */
    phi_f = sc[12]*sc[20];
    k_f = k_f_save[159];
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[13];
    Kc = Kc_save[159];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[18] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 161: CH3 + CH3OH <=> CH3O + CH4 */
    phi_f = sc[12]*sc[20];
    k_f = k_f_save[160];
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[13];
    Kc = Kc_save[160];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[20] -= 1 * qdot;
    wdot[19] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 162: CH3 + C2H4 <=> C2H3 + CH4 */
    phi_f = sc[12]*sc[24];
    k_f = k_f_save[161];
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[13];
    Kc = Kc_save[161];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[24] -= 1 * qdot;
    wdot[23] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 163: CH3 + C2H6 <=> C2H5 + CH4 */
    phi_f = sc[12]*sc[26];
    k_f = k_f_save[162];
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[13];
    Kc = Kc_save[162];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[26] -= 1 * qdot;
    wdot[25] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 164: HCO + H2O <=> H + CO + H2O */
    phi_f = sc[16]*sc[5];
    k_f = k_f_save[163];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14]*sc[5];
    Kc = Kc_save[163];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 165: HCO + M <=> H + CO + M */
    phi_f = sc[16];
    alpha = mixture + sc[0] + -1*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = alpha * k_f_save[164];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[164];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 166: HCO + O2 <=> HO2 + CO */
    phi_f = sc[16]*sc[3];
    k_f = k_f_save[165];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[14];
    Kc = Kc_save[165];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 167: CH2OH + O2 <=> HO2 + CH2O */
    phi_f = sc[18]*sc[3];
    k_f = k_f_save[166];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = Kc_save[166];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[18] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 168: CH3O + O2 <=> HO2 + CH2O */
    phi_f = sc[19]*sc[3];
    k_f = k_f_save[167];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = Kc_save[167];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[19] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 169: C2H + O2 <=> HCO + CO */
    phi_f = sc[21]*sc[3];
    k_f = k_f_save[168];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[14];
    Kc = Kc_save[168];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 170: C2H + H2 <=> H + C2H2 */
    phi_f = sc[21]*sc[0];
    k_f = k_f_save[169];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = Kc_save[169];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[21] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 171: C2H3 + O2 <=> HCO + CH2O */
    phi_f = sc[23]*sc[3];
    k_f = k_f_save[170];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[17];
    Kc = Kc_save[170];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[23] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 172: C2H4 (+M) <=> H2 + C2H2 (+M) */
    phi_f = sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[171];
    redP = 1e-6 * alpha / k_f * 1.58e+51*exp(-9.3*tc[0]-49214.589915701588325*invT);
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
    Kc = Kc_save[171];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[24] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 173: C2H5 + O2 <=> HO2 + C2H4 */
    phi_f = sc[25]*sc[3];
    k_f = k_f_save[172];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[24];
    Kc = Kc_save[172];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[25] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[24] += 1 * qdot;

    /*reaction 174: HCCO + O2 <=> OH + 2 CO */
    phi_f = sc[27]*sc[3];
    k_f = k_f_save[173];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[14]*sc[14];
    Kc = Kc_save[173];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[27] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[14] += 2 * qdot;

    /*reaction 175: 2 HCCO <=> 2 CO + C2H2 */
    phi_f = sc[27]*sc[27];
    k_f = k_f_save[174];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[14]*sc[22];
    Kc = Kc_save[174];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[27] -= 2 * qdot;
    wdot[14] += 2 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 176: O + CH3 => H + H2 + CO */
    phi_f = sc[2]*sc[12];
    k_f = k_f_save[175];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 177: O + C2H4 <=> H + CH2CHO */
    phi_f = sc[2]*sc[24];
    k_f = k_f_save[176];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[33];
    Kc = Kc_save[176];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[24] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[33] += 1 * qdot;

    /*reaction 178: O + C2H5 <=> H + CH3CHO */
    phi_f = sc[2]*sc[25];
    k_f = k_f_save[177];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[34];
    Kc = Kc_save[177];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[25] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[34] += 1 * qdot;

    /*reaction 179: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = k_f_save[178];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = Kc_save[178];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 180: OH + CH3 => H2 + CH2O */
    phi_f = sc[4]*sc[12];
    k_f = k_f_save[179];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[12] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 181: CH + H2 (+M) <=> CH3 (+M) */
    phi_f = sc[9]*sc[0];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[180];
    redP = 1e-12 * alpha / k_f * 4.82e+25*exp(-2.8*tc[0]-296.89783282478464344*invT);
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
    Kc = Kc_save[180];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[12] += 1 * qdot;

    /*reaction 182: CH2 + O2 => 2 H + CO2 */
    phi_f = sc[10]*sc[3];
    k_f = k_f_save[181];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[1] += 2 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 183: CH2 + O2 <=> O + CH2O */
    phi_f = sc[10]*sc[3];
    k_f = k_f_save[182];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[17];
    Kc = Kc_save[182];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 184: CH2 + CH2 => 2 H + C2H2 */
    phi_f = sc[10]*sc[10];
    k_f = k_f_save[183];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[10] -= 1 * qdot;
    wdot[1] += 2 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 185: CH2(S) + H2O => H2 + CH2O */
    phi_f = sc[11]*sc[5];
    k_f = k_f_save[184];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 186: C2H3 + O2 <=> O + CH2CHO */
    phi_f = sc[23]*sc[3];
    k_f = k_f_save[185];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[33];
    Kc = Kc_save[185];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[23] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[33] += 1 * qdot;

    /*reaction 187: C2H3 + O2 <=> HO2 + C2H2 */
    phi_f = sc[23]*sc[3];
    k_f = k_f_save[186];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[22];
    Kc = Kc_save[186];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[23] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[22] += 1 * qdot;

    /*reaction 188: O + CH3CHO <=> OH + CH2CHO */
    phi_f = sc[2]*sc[34];
    k_f = k_f_save[187];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[33];
    Kc = Kc_save[187];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[34] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[33] += 1 * qdot;

    /*reaction 189: O + CH3CHO => OH + CH3 + CO */
    phi_f = sc[2]*sc[34];
    k_f = k_f_save[188];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[34] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 190: O2 + CH3CHO => HO2 + CH3 + CO */
    phi_f = sc[3]*sc[34];
    k_f = k_f_save[189];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[34] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 191: H + CH3CHO <=> CH2CHO + H2 */
    phi_f = sc[1]*sc[34];
    k_f = k_f_save[190];
    q_f = phi_f * k_f;
    phi_r = sc[33]*sc[0];
    Kc = Kc_save[190];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[34] -= 1 * qdot;
    wdot[33] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 192: H + CH3CHO => CH3 + H2 + CO */
    phi_f = sc[1]*sc[34];
    k_f = k_f_save[191];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[34] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 193: OH + CH3CHO => CH3 + H2O + CO */
    phi_f = sc[4]*sc[34];
    k_f = k_f_save[192];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[34] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 194: HO2 + CH3CHO => CH3 + H2O2 + CO */
    phi_f = sc[6]*sc[34];
    k_f = k_f_save[193];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[34] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 195: CH3 + CH3CHO => CH3 + CH4 + CO */
    phi_f = sc[12]*sc[34];
    k_f = k_f_save[194];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[34] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 196: H + CH2CO (+M) <=> CH2CHO (+M) */
    phi_f = sc[1]*sc[28];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[195];
    redP = 1e-12 * alpha / k_f * 1.012e+42*exp(-7.63*tc[0]-1939.3970300113896883*invT);
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
    phi_r = sc[33];
    Kc = Kc_save[195];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[28] -= 1 * qdot;
    wdot[33] += 1 * qdot;

    /*reaction 197: O + CH2CHO => H + CH2 + CO2 */
    phi_f = sc[2]*sc[33];
    k_f = k_f_save[196];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[33] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[15] += 1 * qdot;

    /*reaction 198: O2 + CH2CHO => OH + CO + CH2O */
    phi_f = sc[3]*sc[33];
    k_f = k_f_save[197];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[33] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 199: O2 + CH2CHO => OH + 2 HCO */
    phi_f = sc[3]*sc[33];
    k_f = k_f_save[198];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[33] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[16] += 2 * qdot;

    /*reaction 200: H + CH2CHO <=> CH3 + HCO */
    phi_f = sc[1]*sc[33];
    k_f = k_f_save[199];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[16];
    Kc = Kc_save[199];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[33] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 201: H + CH2CHO <=> CH2CO + H2 */
    phi_f = sc[1]*sc[33];
    k_f = k_f_save[200];
    q_f = phi_f * k_f;
    phi_r = sc[28]*sc[0];
    Kc = Kc_save[200];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[33] -= 1 * qdot;
    wdot[28] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 202: OH + CH2CHO <=> H2O + CH2CO */
    phi_f = sc[4]*sc[33];
    k_f = k_f_save[201];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[28];
    Kc = Kc_save[201];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[33] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[28] += 1 * qdot;

    /*reaction 203: OH + CH2CHO <=> HCO + CH2OH */
    phi_f = sc[4]*sc[33];
    k_f = k_f_save[202];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[18];
    Kc = Kc_save[202];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[33] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[18] += 1 * qdot;

    /*reaction 204: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    phi_f = sc[12]*sc[25];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[203];
    redP = 1e-12 * alpha / k_f * 2.71e+74*exp(-16.82*tc[0]-6574.525738738662767*invT);
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
    phi_r = sc[32];
    Kc = Kc_save[203];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[25] -= 1 * qdot;
    wdot[32] += 1 * qdot;

    /*reaction 205: O + C3H8 <=> OH + C3H7 */
    phi_f = sc[2]*sc[32];
    k_f = k_f_save[204];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[31];
    Kc = Kc_save[204];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[32] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[31] += 1 * qdot;

    /*reaction 206: H + C3H8 <=> C3H7 + H2 */
    phi_f = sc[1]*sc[32];
    k_f = k_f_save[205];
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[0];
    Kc = Kc_save[205];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[32] -= 1 * qdot;
    wdot[31] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 207: OH + C3H8 <=> C3H7 + H2O */
    phi_f = sc[4]*sc[32];
    k_f = k_f_save[206];
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[5];
    Kc = Kc_save[206];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[32] -= 1 * qdot;
    wdot[31] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 208: C3H7 + H2O2 <=> HO2 + C3H8 */
    phi_f = sc[31]*sc[7];
    k_f = k_f_save[207];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[32];
    Kc = Kc_save[207];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[31] -= 1 * qdot;
    wdot[7] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[32] += 1 * qdot;

    /*reaction 209: CH3 + C3H8 <=> C3H7 + CH4 */
    phi_f = sc[12]*sc[32];
    k_f = k_f_save[208];
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[13];
    Kc = Kc_save[208];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[32] -= 1 * qdot;
    wdot[31] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 210: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    phi_f = sc[12]*sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[209];
    redP = 1e-12 * alpha / k_f * 3e+63*exp(-14.6*tc[0]-9143.4468176717564347*invT);
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
    phi_r = sc[31];
    Kc = Kc_save[209];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[24] -= 1 * qdot;
    wdot[31] += 1 * qdot;

    /*reaction 211: O + C3H7 <=> C2H5 + CH2O */
    phi_f = sc[2]*sc[31];
    k_f = k_f_save[210];
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[17];
    Kc = Kc_save[210];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 1 * qdot;
    wdot[31] -= 1 * qdot;
    wdot[25] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 212: H + C3H7 (+M) <=> C3H8 (+M) */
    phi_f = sc[1]*sc[31];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[211];
    redP = 1e-12 * alpha / k_f * 4.42e+61*exp(-13.545*tc[0]-5715.0316735442011122*invT);
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
    phi_r = sc[32];
    Kc = Kc_save[211];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[31] -= 1 * qdot;
    wdot[32] += 1 * qdot;

    /*reaction 213: H + C3H7 <=> CH3 + C2H5 */
    phi_f = sc[1]*sc[31];
    k_f = k_f_save[212];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[25];
    Kc = Kc_save[212];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[31] -= 1 * qdot;
    wdot[12] += 1 * qdot;
    wdot[25] += 1 * qdot;

    /*reaction 214: OH + C3H7 <=> C2H5 + CH2OH */
    phi_f = sc[4]*sc[31];
    k_f = k_f_save[213];
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[18];
    Kc = Kc_save[213];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[31] -= 1 * qdot;
    wdot[25] += 1 * qdot;
    wdot[18] += 1 * qdot;

    /*reaction 215: HO2 + C3H7 <=> O2 + C3H8 */
    phi_f = sc[6]*sc[31];
    k_f = k_f_save[214];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[32];
    Kc = Kc_save[214];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[31] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[32] += 1 * qdot;

    /*reaction 216: HO2 + C3H7 => OH + C2H5 + CH2O */
    phi_f = sc[6]*sc[31];
    k_f = k_f_save[215];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[31] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[25] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 217: CH3 + C3H7 <=> 2 C2H5 */
    phi_f = sc[12]*sc[31];
    k_f = k_f_save[216];
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[25];
    Kc = Kc_save[216];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[31] -= 1 * qdot;
    wdot[25] += 2 * qdot;

    return;
}


/*compute the progress rate for each reaction */
void progressRate(double * qdot, double * sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[35];                /*Gibbs free energy */
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
    for (id = 0; id < 35; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    if (T != T_save)
    {
        T_save = T;

        k_f_save[0] = 1e-12 * 1.2e+17*exp(-1*tc[0]);
        k_f_save[1] = 1e-12 * 5e+17*exp(-1*tc[0]);
        k_f_save[2] = 1e-06 * 38700*exp(2.7*tc[0]-3150.1363279375455022*invT);
        k_f_save[3] = 1e-06 * 2e+13;
        k_f_save[4] = 1e-06 * 9.63e+06*exp(2*tc[0]-2012.866663218878557*invT);
        k_f_save[5] = 1e-06 * 5.7e+13;
        k_f_save[6] = 1e-06 * 8e+13;
        k_f_save[7] = 1e-06 * 1.5e+13;
        k_f_save[8] = 1e-06 * 1.5e+13;
        k_f_save[9] = 1e-06 * 5.06e+13;
        k_f_save[10] = 1e-06 * 1.02e+09*exp(1.5*tc[0]-4327.6633259205891591*invT);
        k_f_save[11] = 1e-06 * 1.8e+10*exp(-1200.1717479442565946*invT);
        k_f_save[12] = 1e-06 * 3e+13;
        k_f_save[13] = 1e-06 * 3e+13;
        k_f_save[14] = 1e-06 * 3.9e+13*exp(-1781.3869969487077469*invT);
        k_f_save[15] = 1e-06 * 1e+13;
        k_f_save[16] = 1e-06 * 1e+13;
        k_f_save[17] = 1e-06 * 388000*exp(2.5*tc[0]-1559.9716639946311716*invT);
        k_f_save[18] = 1e-06 * 130000*exp(2.5*tc[0]-2516.0833290235987079*invT);
        k_f_save[19] = 1e-06 * 5e+13;
        k_f_save[20] = 1e-06 * 1.35e+07*exp(2*tc[0]-956.11166502896742259*invT);
        k_f_save[21] = 1e-06 * 4.6e+19*exp(-1.41*tc[0]-14568.122475046635373*invT);
        k_f_save[22] = 1e-06 * 6.94e+06*exp(2*tc[0]-956.11166502896742259*invT);
        k_f_save[23] = 1e-06 * 3e+13;
        k_f_save[24] = 1e-06 * 1.25e+07*exp(1.83*tc[0]-110.7076664770383303*invT);
        k_f_save[25] = 1e-06 * 2.24e+13;
        k_f_save[26] = 1e-06 * 8.98e+07*exp(1.92*tc[0]-2863.3028284288548093*invT);
        k_f_save[27] = 1e-06 * 1e+14;
        k_f_save[28] = 1e-06 * 1e+13*exp(-4025.733326437757114*invT);
        k_f_save[29] = 1e-06 * 1.75e+12*exp(-679.34249883637164658*invT);
        k_f_save[30] = 1e-06 * 2.5e+12*exp(-24053.756625465601246*invT);
        k_f_save[31] = 1e-06 * 1e+14*exp(-20128.666632188789663*invT);
        k_f_save[32] = 1e-12 * 2.8e+18*exp(-0.86*tc[0]);
        k_f_save[33] = 1e-12 * 2.08e+19*exp(-1.24*tc[0]);
        k_f_save[34] = 1e-12 * 1.126e+19*exp(-0.76*tc[0]);
        k_f_save[35] = 1e-12 * 2.6e+19*exp(-1.24*tc[0]);
        k_f_save[36] = 1e-06 * 2.65e+16*exp(-0.6707*tc[0]-8575.3152019782282878*invT);
        k_f_save[37] = 1e-12 * 1e+18*exp(-1*tc[0]);
        k_f_save[38] = 1e-12 * 9e+16*exp(-0.6*tc[0]);
        k_f_save[39] = 1e-12 * 6e+19*exp(-1.25*tc[0]);
        k_f_save[40] = 1e-12 * 5.5e+20*exp(-2*tc[0]);
        k_f_save[41] = 1e-12 * 2.2e+22*exp(-2*tc[0]);
        k_f_save[42] = 1e-06 * 3.97e+12*exp(-337.65838275496690812*invT);
        k_f_save[43] = 1e-06 * 4.48e+13*exp(-537.43539907944057177*invT);
        k_f_save[44] = 1e-06 * 8.4e+13*exp(-319.54258278599701271*invT);
        k_f_save[45] = 1e-06 * 1.21e+07*exp(2*tc[0]-2616.7266621845424197*invT);
        k_f_save[46] = 1e-06 * 1e+13*exp(-1811.579996896990906*invT);
        k_f_save[47] = 1e-06 * 1.65e+14;
        k_f_save[48] = 1e-06 * 6e+14;
        k_f_save[49] = 1e-06 * 3e+13;
        k_f_save[50] = 1e-06 * 1.39e+16*exp(-0.534*tc[0]-269.72413287132980031*invT);
        k_f_save[51] = 1e-06 * 6.6e+08*exp(1.62*tc[0]-5454.8686573231616421*invT);
        k_f_save[52] = 1e-06 * 1.09e+12*exp(0.48*tc[0]+130.83633310922709825*invT);
        k_f_save[53] = 1e-06 * 7.34e+13;
        k_f_save[54] = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1811.579996896990906*invT);
        k_f_save[55] = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1308.3633310922712099*invT);
        k_f_save[56] = 1e-06 * 5.74e+07*exp(1.9*tc[0]-1379.8200976365412771*invT);
        k_f_save[57] = 1e-06 * 1.055e+12*exp(0.5*tc[0]-43.276633259205894433*invT);
        k_f_save[58] = 1e-06 * 2e+13;
        k_f_save[59] = 1e-06 * 1.65e+11*exp(0.65*tc[0]+142.9135330885404187*invT);
        k_f_save[60] = 1e-06 * 3.28e+13*exp(-0.09*tc[0]-306.96216614087899188*invT);
        k_f_save[61] = 1e-06 * 2.43e+12*exp(0.515*tc[0]-25.160833290235984805*invT);
        k_f_save[62] = 1e-06 * 4.15e+07*exp(1.63*tc[0]-968.1888650082806862*invT);
        k_f_save[63] = 1e-06 * 2e+13;
        k_f_save[64] = 1e-06 * 1.5e+12*exp(0.5*tc[0]+55.35383323851916515*invT);
        k_f_save[65] = 1e-06 * 2.62e+14*exp(-0.23*tc[0]-538.44183241105008619*invT);
        k_f_save[66] = 1e-06 * 1.7e+07*exp(2.1*tc[0]-2450.6651624689848177*invT);
        k_f_save[67] = 1e-06 * 4.2e+06*exp(2.1*tc[0]-2450.6651624689848177*invT);
        k_f_save[68] = 1e-06 * 1e+17*exp(-1*tc[0]);
        k_f_save[69] = 1e-06 * 5.6e+12*exp(-1207.7199979313272706*invT);
        k_f_save[70] = 1e-06 * 6.08e+12*exp(0.27*tc[0]-140.90066642532150354*invT);
        k_f_save[71] = 1e-06 * 3e+13;
        k_f_save[72] = 1e-06 * 5.4e+11*exp(0.454*tc[0]-915.85433176458980142*invT);
        k_f_save[73] = 1e-06 * 1.325e+06*exp(2.53*tc[0]-6159.3719894497689893*invT);
        k_f_save[74] = 1e-06 * 5.21e+17*exp(-0.99*tc[0]-795.08233197145705162*invT);
        k_f_save[75] = 1e-06 * 2e+12;
        k_f_save[76] = 1e-06 * 1.15e+08*exp(1.9*tc[0]-3789.2214935095394139*invT);
        k_f_save[77] = 1e-06 * 1e+14;
        k_f_save[78] = 1e-06 * 5e+13*exp(-4025.733326437757114*invT);
        k_f_save[79] = 1e-06 * 1.13e+13*exp(-1725.0267303785790318*invT);
        k_f_save[80] = 1e-06 * 1e+13;
        k_f_save[81] = 1e-06 * 4.3e+07*exp(1.5*tc[0]-40056.046598055690993*invT);
        k_f_save[82] = 1e-06 * 2.16e+08*exp(1.51*tc[0]-1726.0331637101885462*invT);
        k_f_save[83] = 1e-06 * 7.4e+13*exp(-0.37*tc[0]);
        k_f_save[84] = 1e-06 * 35700*exp(2.4*tc[0]+1061.7871648479585929*invT);
        k_f_save[85] = 1e-06 * 1.45e+13*exp(+251.60833290235981963*invT);
        k_f_save[86] = 1e-06 * 2e+12*exp(-214.8735162986153*invT);
        k_f_save[87] = 1e-06 * 1.7e+18*exp(-14799.602141316805501*invT);
        k_f_save[88] = 1e-06 * 5e+13;
        k_f_save[89] = 1e-06 * 3e+13;
        k_f_save[90] = 1e-06 * 2e+13;
        k_f_save[91] = 1e-06 * 1.13e+07*exp(2*tc[0]-1509.6499974141590883*invT);
        k_f_save[92] = 1e-06 * 3e+13;
        k_f_save[93] = 1e-06 * 2.79e+18*exp(-1.43*tc[0]-669.27816552027718444*invT);
        k_f_save[94] = 1e-06 * 5.6e+07*exp(1.6*tc[0]-2727.4343286615808211*invT);
        k_f_save[95] = 1e-06 * 6.44e+17*exp(-1.34*tc[0]-713.05801544528776503*invT);
        k_f_save[96] = 1e-06 * 1e+08*exp(1.6*tc[0]-1570.0359973107254064*invT);
        k_f_save[97] = 1e-06 * 4.76e+07*exp(1.228*tc[0]-35.225166606330375885*invT);
        k_f_save[98] = 1e-06 * 5e+13;
        k_f_save[99] = 1e-06 * 3.43e+09*exp(1.18*tc[0]+224.93784961470970529*invT);
        k_f_save[100] = 1e-06 * 5e+12;
        k_f_save[101] = 1e-06 * 5e+12;
        k_f_save[102] = 1e-06 * 1.44e+06*exp(2*tc[0]+422.70199927596451062*invT);
        k_f_save[103] = 1e-06 * 6.3e+06*exp(2*tc[0]-754.82499870707954415*invT);
        k_f_save[104] = 1e-06 * 2e+13;
        k_f_save[105] = 1e-06 * 0.000218*exp(4.5*tc[0]+503.21666580471963925*invT);
        k_f_save[106] = 1e-06 * 504000*exp(2.3*tc[0]-6793.4249883637157836*invT);
        k_f_save[107] = 1e-06 * 3.37e+07*exp(2*tc[0]-7045.0333212660752906*invT);
        k_f_save[108] = 1e-06 * 0.000483*exp(4*tc[0]+1006.4333316094392785*invT);
        k_f_save[109] = 1e-06 * 5e+12;
        k_f_save[110] = 1e-06 * 3.6e+06*exp(2*tc[0]-1258.0416645117993539*invT);
        k_f_save[111] = 1e-06 * 3.54e+06*exp(2.12*tc[0]-437.79849925010609013*invT);
        k_f_save[112] = 1e-06 * 7.5e+12*exp(-1006.4333316094392785*invT);
        k_f_save[113] = 1e-06 * 1.3e+11*exp(+820.24316526169309327*invT);
        k_f_save[114] = 1e-06 * 4.2e+14*exp(-6038.5999896566363532*invT);
        k_f_save[115] = 1e-06 * 2e+13;
        k_f_save[116] = 1e-06 * 1e+12;
        k_f_save[117] = 1e-06 * 3.78e+13;
        k_f_save[118] = 1e-06 * 1.5e+14*exp(-11875.913312991384373*invT);
        k_f_save[119] = 1e-06 * 5.6e+06*exp(2*tc[0]-6038.5999896566363532*invT);
        k_f_save[120] = 1e-06 * 5.8e+13*exp(-289.85279950351855405*invT);
        k_f_save[121] = 1e-06 * 5e+13;
        k_f_save[122] = 1e-06 * 5e+13;
        k_f_save[123] = 1e-06 * 6.71e+13;
        k_f_save[124] = 1e-06 * 1.08e+14*exp(-1565.003830652678289*invT);
        k_f_save[125] = 1e-06 * 5.71e+12*exp(+379.92858268256338761*invT);
        k_f_save[126] = 1e-06 * 4e+13;
        k_f_save[127] = 1e-06 * 3e+13;
        k_f_save[128] = 1e-06 * 6e+13;
        k_f_save[129] = 1e-06 * 5e+13;
        k_f_save[130] = 1e-06 * 1.9e+14*exp(-7946.7975863881329133*invT);
        k_f_save[131] = 1e-06 * 9.46e+13*exp(+259.15658288943063781*invT);
        k_f_save[132] = 1e-06 * 5e+13;
        k_f_save[133] = 1e-06 * 5e+12*exp(-754.82499870707954415*invT);
        k_f_save[134] = 1e-06 * 500000*exp(2*tc[0]-3638.256493768123164*invT);
        k_f_save[135] = 1e-06 * 1.6e+15*exp(-6010.4198563715717682*invT);
        k_f_save[136] = 1e-06 * 4e+13;
        k_f_save[137] = 1e-06 * 2.46e+06*exp(2*tc[0]-4161.6018262050320118*invT);
        k_f_save[138] = 1e-06 * 8.1e+11*exp(0.5*tc[0]-2269.5071627792858635*invT);
        k_f_save[139] = 1e-06 * 3e+13;
        k_f_save[140] = 1e-06 * 1.5e+13*exp(-301.92999948283181766*invT);
        k_f_save[141] = 1e-06 * 2.8e+13;
        k_f_save[142] = 1e-06 * 1.2e+13;
        k_f_save[143] = 1e-06 * 7e+13;
        k_f_save[144] = 1e-06 * 4.82e+17*exp(-1.16*tc[0]-576.18308234640414867*invT);
        k_f_save[145] = 1e-06 * 3e+13;
        k_f_save[146] = 1e-06 * 1.2e+13*exp(+286.83349950869023814*invT);
        k_f_save[147] = 1e-06 * 1.6e+13*exp(+286.83349950869023814*invT);
        k_f_save[148] = 1e-06 * 9e+12;
        k_f_save[149] = 1e-06 * 7e+12;
        k_f_save[150] = 1e-06 * 1.4e+13;
        k_f_save[151] = 1e-06 * 4e+13*exp(+276.76916619259583285*invT);
        k_f_save[152] = 1e-06 * 3.56e+13*exp(-15338.0439737278557*invT);
        k_f_save[153] = 1e-06 * 2.31e+12*exp(-10222.846565822879711*invT);
        k_f_save[154] = 1e-06 * 24500*exp(2.47*tc[0]-2606.6623288684481849*invT);
        k_f_save[155] = 1e-06 * 6.77e+16*exp(-1.18*tc[0]-329.10369943628666078*invT);
        k_f_save[156] = 1e-06 * 6.84e+12*exp(0.1*tc[0]-5334.096657530029006*invT);
        k_f_save[157] = 1e-06 * 2.648e+13;
        k_f_save[158] = 1e-06 * 3320*exp(2.81*tc[0]-2948.8496616156576238*invT);
        k_f_save[159] = 1e-06 * 3e+07*exp(1.5*tc[0]-5001.973658098913802*invT);
        k_f_save[160] = 1e-06 * 1e+07*exp(1.5*tc[0]-5001.973658098913802*invT);
        k_f_save[161] = 1e-06 * 227000*exp(2*tc[0]-4629.5933254034207494*invT);
        k_f_save[162] = 1e-06 * 6.14e+06*exp(1.74*tc[0]-5258.6141576593208811*invT);
        k_f_save[163] = 1e-06 * 1.5e+18*exp(-1*tc[0]-8554.6833186802341515*invT);
        k_f_save[164] = 1e-06 * 1.87e+17*exp(-1*tc[0]-8554.6833186802341515*invT);
        k_f_save[165] = 1e-06 * 1.345e+13*exp(-201.28666632188787844*invT);
        k_f_save[166] = 1e-06 * 1.8e+13*exp(-452.89499922424772649*invT);
        k_f_save[167] = 1e-06 * 4.28e-13*exp(7.6*tc[0]+1776.3548302906606295*invT);
        k_f_save[168] = 1e-06 * 1e+13*exp(+379.92858268256338761*invT);
        k_f_save[169] = 1e-06 * 5.68e+10*exp(0.9*tc[0]-1002.9108149488064328*invT);
        k_f_save[170] = 1e-06 * 4.58e+16*exp(-1.39*tc[0]-510.76491579179048586*invT);
        k_f_save[171] = 1 * 8e+12*exp(0.44*tc[0]-43664.110091875525541*invT);
        k_f_save[172] = 1e-06 * 8.4e+11*exp(-1949.9645799932889076*invT);
        k_f_save[173] = 1e-06 * 3.2e+12*exp(-429.7470325972306*invT);
        k_f_save[174] = 1e-06 * 1e+13;
        k_f_save[175] = 1e-06 * 3.37e+13;
        k_f_save[176] = 1e-06 * 6.7e+06*exp(1.83*tc[0]-110.7076664770383303*invT);
        k_f_save[177] = 1e-06 * 1.096e+14;
        k_f_save[178] = 1e-06 * 5e+15*exp(-8720.7448183957912988*invT);
        k_f_save[179] = 1e-06 * 8e+09*exp(0.5*tc[0]+883.1452484872830837*invT);
        k_f_save[180] = 1e-06 * 1.97e+12*exp(0.43*tc[0]+186.1901663477462705*invT);
        k_f_save[181] = 1e-06 * 5.8e+12*exp(-754.82499870707954415*invT);
        k_f_save[182] = 1e-06 * 2.4e+12*exp(-754.82499870707954415*invT);
        k_f_save[183] = 1e-06 * 2e+14*exp(-5529.8479405280650099*invT);
        k_f_save[184] = 1e-06 * 6.82e+10*exp(0.25*tc[0]+470.50758252741292154*invT);
        k_f_save[185] = 1e-06 * 3.03e+11*exp(0.29*tc[0]-5.5353833238519163373*invT);
        k_f_save[186] = 1e-06 * 1.337e+06*exp(1.61*tc[0]+193.23519966901235989*invT);
        k_f_save[187] = 1e-06 * 5.84e+12*exp(-909.81573177493316962*invT);
        k_f_save[188] = 1e-06 * 5.84e+12*exp(-909.81573177493316962*invT);
        k_f_save[189] = 1e-06 * 3.01e+13*exp(-19700.932466254773317*invT);
        k_f_save[190] = 1e-06 * 2.05e+09*exp(1.16*tc[0]-1210.2360812603508293*invT);
        k_f_save[191] = 1e-06 * 2.05e+09*exp(1.16*tc[0]-1210.2360812603508293*invT);
        k_f_save[192] = 1e-06 * 2.343e+10*exp(0.73*tc[0]+560.08014904065305473*invT);
        k_f_save[193] = 1e-06 * 3.01e+12*exp(-5999.8523063896727763*invT);
        k_f_save[194] = 1e-06 * 2.72e+06*exp(1.77*tc[0]-2979.0426615639403281*invT);
        k_f_save[195] = 1e-06 * 4.865e+11*exp(0.422*tc[0]+883.1452484872830837*invT);
        k_f_save[196] = 1e-06 * 1.5e+14;
        k_f_save[197] = 1e-06 * 1.81e+10;
        k_f_save[198] = 1e-06 * 2.35e+10;
        k_f_save[199] = 1e-06 * 2.2e+13;
        k_f_save[200] = 1e-06 * 1.1e+13;
        k_f_save[201] = 1e-06 * 1.2e+13;
        k_f_save[202] = 1e-06 * 3.01e+13;
        k_f_save[203] = 1e-06 * 9.43e+12;
        k_f_save[204] = 1e-06 * 193000*exp(2.68*tc[0]-1869.9531301303384225*invT);
        k_f_save[205] = 1e-06 * 1.32e+06*exp(2.54*tc[0]-3399.7317941766864351*invT);
        k_f_save[206] = 1e-06 * 3.16e+07*exp(1.8*tc[0]-470.00436586160816432*invT);
        k_f_save[207] = 1e-06 * 378*exp(2.72*tc[0]-754.82499870707954415*invT);
        k_f_save[208] = 1e-06 * 0.903*exp(3.65*tc[0]-3600.0120271669643444*invT);
        k_f_save[209] = 1e-06 * 2.55e+06*exp(1.6*tc[0]-2868.3349950869019267*invT);
        k_f_save[210] = 1e-06 * 9.64e+13;
        k_f_save[211] = 1e-06 * 3.613e+13;
        k_f_save[212] = 1e-06 * 4.06e+06*exp(2.19*tc[0]-447.86283256620055226*invT);
        k_f_save[213] = 1e-06 * 2.41e+13;
        k_f_save[214] = 1e-06 * 2.55e+10*exp(0.255*tc[0]+474.53331585385063818*invT);
        k_f_save[215] = 1e-06 * 2.41e+13;
        k_f_save[216] = 1e-06 * 1.927e+13*exp(-0.32*tc[0]);

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
        Kc_save[36] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));
        Kc_save[37] = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));
        Kc_save[38] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));
        Kc_save[39] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));
        Kc_save[40] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[15]) - (g_RT[0] + g_RT[15]));
        Kc_save[41] = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));
        Kc_save[42] = exp((g_RT[1] + g_RT[6]) - (g_RT[2] + g_RT[5]));
        Kc_save[43] = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));
        Kc_save[44] = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));
        Kc_save[45] = exp((g_RT[1] + g_RT[7]) - (g_RT[6] + g_RT[0]));
        Kc_save[46] = exp((g_RT[1] + g_RT[7]) - (g_RT[4] + g_RT[5]));
        Kc_save[47] = exp((g_RT[1] + g_RT[9]) - (g_RT[8] + g_RT[0]));
        Kc_save[48] = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[12]));
        Kc_save[49] = exp((g_RT[1] + g_RT[11]) - (g_RT[9] + g_RT[0]));
        Kc_save[50] = 1.0 / (refC) * exp((g_RT[1] + g_RT[12]) - (g_RT[13]));
        Kc_save[51] = exp((g_RT[1] + g_RT[13]) - (g_RT[12] + g_RT[0]));
        Kc_save[52] = 1.0 / (refC) * exp((g_RT[1] + g_RT[16]) - (g_RT[17]));
        Kc_save[53] = exp((g_RT[1] + g_RT[16]) - (g_RT[0] + g_RT[14]));
        Kc_save[54] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[18]));
        Kc_save[55] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[19]));
        Kc_save[56] = exp((g_RT[1] + g_RT[17]) - (g_RT[16] + g_RT[0]));
        Kc_save[57] = 1.0 / (refC) * exp((g_RT[1] + g_RT[18]) - (g_RT[20]));
        Kc_save[58] = exp((g_RT[1] + g_RT[18]) - (g_RT[0] + g_RT[17]));
        Kc_save[59] = exp((g_RT[1] + g_RT[18]) - (g_RT[4] + g_RT[12]));
        Kc_save[60] = exp((g_RT[1] + g_RT[18]) - (g_RT[11] + g_RT[5]));
        Kc_save[61] = 1.0 / (refC) * exp((g_RT[1] + g_RT[19]) - (g_RT[20]));
        Kc_save[62] = exp((g_RT[1] + g_RT[19]) - (g_RT[1] + g_RT[18]));
        Kc_save[63] = exp((g_RT[1] + g_RT[19]) - (g_RT[0] + g_RT[17]));
        Kc_save[64] = exp((g_RT[1] + g_RT[19]) - (g_RT[4] + g_RT[12]));
        Kc_save[65] = exp((g_RT[1] + g_RT[19]) - (g_RT[11] + g_RT[5]));
        Kc_save[66] = exp((g_RT[1] + g_RT[20]) - (g_RT[18] + g_RT[0]));
        Kc_save[67] = exp((g_RT[1] + g_RT[20]) - (g_RT[19] + g_RT[0]));
        Kc_save[68] = 1.0 / (refC) * exp((g_RT[1] + g_RT[21]) - (g_RT[22]));
        Kc_save[69] = 1.0 / (refC) * exp((g_RT[1] + g_RT[22]) - (g_RT[23]));
        Kc_save[70] = 1.0 / (refC) * exp((g_RT[1] + g_RT[23]) - (g_RT[24]));
        Kc_save[71] = exp((g_RT[1] + g_RT[23]) - (g_RT[0] + g_RT[22]));
        Kc_save[72] = 1.0 / (refC) * exp((g_RT[1] + g_RT[24]) - (g_RT[25]));
        Kc_save[73] = exp((g_RT[1] + g_RT[24]) - (g_RT[23] + g_RT[0]));
        Kc_save[74] = 1.0 / (refC) * exp((g_RT[1] + g_RT[25]) - (g_RT[26]));
        Kc_save[75] = exp((g_RT[1] + g_RT[25]) - (g_RT[0] + g_RT[24]));
        Kc_save[76] = exp((g_RT[1] + g_RT[26]) - (g_RT[25] + g_RT[0]));
        Kc_save[77] = exp((g_RT[1] + g_RT[27]) - (g_RT[11] + g_RT[14]));
        Kc_save[78] = exp((g_RT[1] + g_RT[28]) - (g_RT[27] + g_RT[0]));
        Kc_save[79] = exp((g_RT[1] + g_RT[28]) - (g_RT[12] + g_RT[14]));
        Kc_save[80] = exp((g_RT[1] + g_RT[29]) - (g_RT[1] + g_RT[28]));
        Kc_save[81] = 1.0 / (refC) * exp((g_RT[0] + g_RT[14]) - (g_RT[17]));
        Kc_save[82] = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));
        Kc_save[83] = 1.0 / (refC) * exp((2 * g_RT[4]) - (g_RT[7]));
        Kc_save[84] = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));
        Kc_save[85] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
        Kc_save[86] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
        Kc_save[87] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
        Kc_save[88] = exp((g_RT[4] + g_RT[8]) - (g_RT[1] + g_RT[14]));
        Kc_save[89] = exp((g_RT[4] + g_RT[9]) - (g_RT[1] + g_RT[16]));
        Kc_save[90] = exp((g_RT[4] + g_RT[10]) - (g_RT[1] + g_RT[17]));
        Kc_save[91] = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[5]));
        Kc_save[92] = exp((g_RT[4] + g_RT[11]) - (g_RT[1] + g_RT[17]));
        Kc_save[93] = 1.0 / (refC) * exp((g_RT[4] + g_RT[12]) - (g_RT[20]));
        Kc_save[94] = exp((g_RT[4] + g_RT[12]) - (g_RT[10] + g_RT[5]));
        Kc_save[95] = exp((g_RT[4] + g_RT[12]) - (g_RT[11] + g_RT[5]));
        Kc_save[96] = exp((g_RT[4] + g_RT[13]) - (g_RT[12] + g_RT[5]));
        Kc_save[97] = exp((g_RT[4] + g_RT[14]) - (g_RT[1] + g_RT[15]));
        Kc_save[98] = exp((g_RT[4] + g_RT[16]) - (g_RT[5] + g_RT[14]));
        Kc_save[99] = exp((g_RT[4] + g_RT[17]) - (g_RT[16] + g_RT[5]));
        Kc_save[100] = exp((g_RT[4] + g_RT[18]) - (g_RT[5] + g_RT[17]));
        Kc_save[101] = exp((g_RT[4] + g_RT[19]) - (g_RT[5] + g_RT[17]));
        Kc_save[102] = exp((g_RT[4] + g_RT[20]) - (g_RT[18] + g_RT[5]));
        Kc_save[103] = exp((g_RT[4] + g_RT[20]) - (g_RT[19] + g_RT[5]));
        Kc_save[104] = exp((g_RT[4] + g_RT[21]) - (g_RT[1] + g_RT[27]));
        Kc_save[105] = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[28]));
        Kc_save[106] = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[29]));
        Kc_save[107] = exp((g_RT[4] + g_RT[22]) - (g_RT[21] + g_RT[5]));
        Kc_save[108] = exp((g_RT[4] + g_RT[22]) - (g_RT[12] + g_RT[14]));
        Kc_save[109] = exp((g_RT[4] + g_RT[23]) - (g_RT[5] + g_RT[22]));
        Kc_save[110] = exp((g_RT[4] + g_RT[24]) - (g_RT[23] + g_RT[5]));
        Kc_save[111] = exp((g_RT[4] + g_RT[26]) - (g_RT[25] + g_RT[5]));
        Kc_save[112] = exp((g_RT[4] + g_RT[28]) - (g_RT[27] + g_RT[5]));
        Kc_save[113] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
        Kc_save[114] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
        Kc_save[115] = exp((g_RT[6] + g_RT[10]) - (g_RT[4] + g_RT[17]));
        Kc_save[116] = exp((g_RT[6] + g_RT[12]) - (g_RT[3] + g_RT[13]));
        Kc_save[117] = exp((g_RT[6] + g_RT[12]) - (g_RT[4] + g_RT[19]));
        Kc_save[118] = exp((g_RT[6] + g_RT[14]) - (g_RT[4] + g_RT[15]));
        Kc_save[119] = exp((g_RT[6] + g_RT[17]) - (g_RT[16] + g_RT[7]));
        Kc_save[120] = exp((g_RT[8] + g_RT[3]) - (g_RT[2] + g_RT[14]));
        Kc_save[121] = exp((g_RT[8] + g_RT[10]) - (g_RT[1] + g_RT[21]));
        Kc_save[122] = exp((g_RT[8] + g_RT[12]) - (g_RT[1] + g_RT[22]));
        Kc_save[123] = exp((g_RT[9] + g_RT[3]) - (g_RT[2] + g_RT[16]));
        Kc_save[124] = exp((g_RT[9] + g_RT[0]) - (g_RT[1] + g_RT[10]));
        Kc_save[125] = exp((g_RT[9] + g_RT[5]) - (g_RT[1] + g_RT[17]));
        Kc_save[126] = exp((g_RT[9] + g_RT[10]) - (g_RT[1] + g_RT[22]));
        Kc_save[127] = exp((g_RT[9] + g_RT[12]) - (g_RT[1] + g_RT[23]));
        Kc_save[128] = exp((g_RT[9] + g_RT[13]) - (g_RT[1] + g_RT[24]));
        Kc_save[129] = 1.0 / (refC) * exp((g_RT[9] + g_RT[14]) - (g_RT[27]));
        Kc_save[130] = exp((g_RT[9] + g_RT[15]) - (g_RT[16] + g_RT[14]));
        Kc_save[131] = exp((g_RT[9] + g_RT[17]) - (g_RT[1] + g_RT[28]));
        Kc_save[132] = exp((g_RT[9] + g_RT[27]) - (g_RT[14] + g_RT[22]));
        Kc_save[133] = refC * exp((g_RT[10] + g_RT[3]) - (g_RT[4] + g_RT[1] + g_RT[14]));
        Kc_save[134] = exp((g_RT[10] + g_RT[0]) - (g_RT[1] + g_RT[12]));
        Kc_save[135] = exp((2 * g_RT[10]) - (g_RT[0] + g_RT[22]));
        Kc_save[136] = exp((g_RT[10] + g_RT[12]) - (g_RT[1] + g_RT[24]));
        Kc_save[137] = exp((g_RT[10] + g_RT[13]) - (2 * g_RT[12]));
        Kc_save[138] = 1.0 / (refC) * exp((g_RT[10] + g_RT[14]) - (g_RT[28]));
        Kc_save[139] = exp((g_RT[10] + g_RT[27]) - (g_RT[23] + g_RT[14]));
        Kc_save[140] = exp((g_RT[11] + g_RT[30]) - (g_RT[10] + g_RT[30]));
        Kc_save[141] = refC * exp((g_RT[11] + g_RT[3]) - (g_RT[1] + g_RT[4] + g_RT[14]));
        Kc_save[142] = exp((g_RT[11] + g_RT[3]) - (g_RT[14] + g_RT[5]));
        Kc_save[143] = exp((g_RT[11] + g_RT[0]) - (g_RT[12] + g_RT[1]));
        Kc_save[144] = 1.0 / (refC) * exp((g_RT[11] + g_RT[5]) - (g_RT[20]));
        Kc_save[145] = exp((g_RT[11] + g_RT[5]) - (g_RT[10] + g_RT[5]));
        Kc_save[146] = exp((g_RT[11] + g_RT[12]) - (g_RT[1] + g_RT[24]));
        Kc_save[147] = exp((g_RT[11] + g_RT[13]) - (2 * g_RT[12]));
        Kc_save[148] = exp((g_RT[11] + g_RT[14]) - (g_RT[10] + g_RT[14]));
        Kc_save[149] = exp((g_RT[11] + g_RT[15]) - (g_RT[10] + g_RT[15]));
        Kc_save[150] = exp((g_RT[11] + g_RT[15]) - (g_RT[14] + g_RT[17]));
        Kc_save[151] = exp((g_RT[11] + g_RT[26]) - (g_RT[12] + g_RT[25]));
        Kc_save[152] = exp((g_RT[12] + g_RT[3]) - (g_RT[2] + g_RT[19]));
        Kc_save[153] = exp((g_RT[12] + g_RT[3]) - (g_RT[4] + g_RT[17]));
        Kc_save[154] = exp((g_RT[12] + g_RT[7]) - (g_RT[6] + g_RT[13]));
        Kc_save[155] = 1.0 / (refC) * exp((2 * g_RT[12]) - (g_RT[26]));
        Kc_save[156] = exp((2 * g_RT[12]) - (g_RT[1] + g_RT[25]));
        Kc_save[157] = exp((g_RT[12] + g_RT[16]) - (g_RT[13] + g_RT[14]));
        Kc_save[158] = exp((g_RT[12] + g_RT[17]) - (g_RT[16] + g_RT[13]));
        Kc_save[159] = exp((g_RT[12] + g_RT[20]) - (g_RT[18] + g_RT[13]));
        Kc_save[160] = exp((g_RT[12] + g_RT[20]) - (g_RT[19] + g_RT[13]));
        Kc_save[161] = exp((g_RT[12] + g_RT[24]) - (g_RT[23] + g_RT[13]));
        Kc_save[162] = exp((g_RT[12] + g_RT[26]) - (g_RT[25] + g_RT[13]));
        Kc_save[163] = refC * exp((g_RT[16] + g_RT[5]) - (g_RT[1] + g_RT[14] + g_RT[5]));
        Kc_save[164] = refC * exp((g_RT[16]) - (g_RT[1] + g_RT[14]));
        Kc_save[165] = exp((g_RT[16] + g_RT[3]) - (g_RT[6] + g_RT[14]));
        Kc_save[166] = exp((g_RT[18] + g_RT[3]) - (g_RT[6] + g_RT[17]));
        Kc_save[167] = exp((g_RT[19] + g_RT[3]) - (g_RT[6] + g_RT[17]));
        Kc_save[168] = exp((g_RT[21] + g_RT[3]) - (g_RT[16] + g_RT[14]));
        Kc_save[169] = exp((g_RT[21] + g_RT[0]) - (g_RT[1] + g_RT[22]));
        Kc_save[170] = exp((g_RT[23] + g_RT[3]) - (g_RT[16] + g_RT[17]));
        Kc_save[171] = refC * exp((g_RT[24]) - (g_RT[0] + g_RT[22]));
        Kc_save[172] = exp((g_RT[25] + g_RT[3]) - (g_RT[6] + g_RT[24]));
        Kc_save[173] = refC * exp((g_RT[27] + g_RT[3]) - (g_RT[4] + 2 * g_RT[14]));
        Kc_save[174] = refC * exp((2 * g_RT[27]) - (2 * g_RT[14] + g_RT[22]));
        Kc_save[175] = refC * exp((g_RT[2] + g_RT[12]) - (g_RT[1] + g_RT[0] + g_RT[14]));
        Kc_save[176] = exp((g_RT[2] + g_RT[24]) - (g_RT[1] + g_RT[33]));
        Kc_save[177] = exp((g_RT[2] + g_RT[25]) - (g_RT[1] + g_RT[34]));
        Kc_save[178] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
        Kc_save[179] = exp((g_RT[4] + g_RT[12]) - (g_RT[0] + g_RT[17]));
        Kc_save[180] = 1.0 / (refC) * exp((g_RT[9] + g_RT[0]) - (g_RT[12]));
        Kc_save[181] = refC * exp((g_RT[10] + g_RT[3]) - (2 * g_RT[1] + g_RT[15]));
        Kc_save[182] = exp((g_RT[10] + g_RT[3]) - (g_RT[2] + g_RT[17]));
        Kc_save[183] = refC * exp((g_RT[10] + g_RT[10]) - (2 * g_RT[1] + g_RT[22]));
        Kc_save[184] = exp((g_RT[11] + g_RT[5]) - (g_RT[0] + g_RT[17]));
        Kc_save[185] = exp((g_RT[23] + g_RT[3]) - (g_RT[2] + g_RT[33]));
        Kc_save[186] = exp((g_RT[23] + g_RT[3]) - (g_RT[6] + g_RT[22]));
        Kc_save[187] = exp((g_RT[2] + g_RT[34]) - (g_RT[4] + g_RT[33]));
        Kc_save[188] = refC * exp((g_RT[2] + g_RT[34]) - (g_RT[4] + g_RT[12] + g_RT[14]));
        Kc_save[189] = refC * exp((g_RT[3] + g_RT[34]) - (g_RT[6] + g_RT[12] + g_RT[14]));
        Kc_save[190] = exp((g_RT[1] + g_RT[34]) - (g_RT[33] + g_RT[0]));
        Kc_save[191] = refC * exp((g_RT[1] + g_RT[34]) - (g_RT[12] + g_RT[0] + g_RT[14]));
        Kc_save[192] = refC * exp((g_RT[4] + g_RT[34]) - (g_RT[12] + g_RT[5] + g_RT[14]));
        Kc_save[193] = refC * exp((g_RT[6] + g_RT[34]) - (g_RT[12] + g_RT[7] + g_RT[14]));
        Kc_save[194] = refC * exp((g_RT[12] + g_RT[34]) - (g_RT[12] + g_RT[13] + g_RT[14]));
        Kc_save[195] = 1.0 / (refC) * exp((g_RT[1] + g_RT[28]) - (g_RT[33]));
        Kc_save[196] = refC * exp((g_RT[2] + g_RT[33]) - (g_RT[1] + g_RT[10] + g_RT[15]));
        Kc_save[197] = refC * exp((g_RT[3] + g_RT[33]) - (g_RT[4] + g_RT[14] + g_RT[17]));
        Kc_save[198] = refC * exp((g_RT[3] + g_RT[33]) - (g_RT[4] + 2 * g_RT[16]));
        Kc_save[199] = exp((g_RT[1] + g_RT[33]) - (g_RT[12] + g_RT[16]));
        Kc_save[200] = exp((g_RT[1] + g_RT[33]) - (g_RT[28] + g_RT[0]));
        Kc_save[201] = exp((g_RT[4] + g_RT[33]) - (g_RT[5] + g_RT[28]));
        Kc_save[202] = exp((g_RT[4] + g_RT[33]) - (g_RT[16] + g_RT[18]));
        Kc_save[203] = 1.0 / (refC) * exp((g_RT[12] + g_RT[25]) - (g_RT[32]));
        Kc_save[204] = exp((g_RT[2] + g_RT[32]) - (g_RT[4] + g_RT[31]));
        Kc_save[205] = exp((g_RT[1] + g_RT[32]) - (g_RT[31] + g_RT[0]));
        Kc_save[206] = exp((g_RT[4] + g_RT[32]) - (g_RT[31] + g_RT[5]));
        Kc_save[207] = exp((g_RT[31] + g_RT[7]) - (g_RT[6] + g_RT[32]));
        Kc_save[208] = exp((g_RT[12] + g_RT[32]) - (g_RT[31] + g_RT[13]));
        Kc_save[209] = 1.0 / (refC) * exp((g_RT[12] + g_RT[24]) - (g_RT[31]));
        Kc_save[210] = exp((g_RT[2] + g_RT[31]) - (g_RT[25] + g_RT[17]));
        Kc_save[211] = 1.0 / (refC) * exp((g_RT[1] + g_RT[31]) - (g_RT[32]));
        Kc_save[212] = exp((g_RT[1] + g_RT[31]) - (g_RT[12] + g_RT[25]));
        Kc_save[213] = exp((g_RT[4] + g_RT[31]) - (g_RT[25] + g_RT[18]));
        Kc_save[214] = exp((g_RT[6] + g_RT[31]) - (g_RT[3] + g_RT[32]));
        Kc_save[215] = refC * exp((g_RT[6] + g_RT[31]) - (g_RT[4] + g_RT[25] + g_RT[17]));
        Kc_save[216] = exp((g_RT[12] + g_RT[31]) - (2 * g_RT[25]));
    }

    /*reaction 1: 2 O + M <=> O2 + M */
    phi_f = sc[2]*sc[2];
    alpha = mixture + 1.4*sc[0] + 14.4*sc[5] + sc[13] + 0.75*sc[14] + 2.6*sc[15] + 2*sc[26];
    k_f = alpha * k_f_save[0];
    q_f = phi_f * k_f;
    phi_r = sc[3];
    Kc = Kc_save[0];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[0] = q_f - q_r;

    /*reaction 2: O + H + M <=> OH + M */
    phi_f = sc[2]*sc[1];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
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

    /*reaction 12: O + CO (+M) <=> CO2 (+M) */
    phi_f = sc[2]*sc[14];
    alpha = mixture + sc[0] + 5*sc[3] + 5*sc[5] + sc[13] + 0.5*sc[14] + 2.5*sc[15] + 2*sc[26];
    k_f = k_f_save[11];
    redP = 1e-12 * alpha / k_f * 6.02e+14*exp(-1509.6499974141590883*invT);
    F = redP / (1 + redP);
    k_f *= F;
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
    alpha = mixture + -1*sc[3] + -1*sc[5] + -0.25*sc[14] + 0.5*sc[15] + 0.5*sc[26] + -1*sc[30];
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

    /*reaction 37: H + O2 <=> O + OH */
    phi_f = sc[1]*sc[3];
    k_f = k_f_save[36];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = Kc_save[36];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[36] = q_f - q_r;

    /*reaction 38: 2 H + M <=> H2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[0] + -1*sc[5] + sc[13] + -1*sc[15] + 2*sc[26];
    k_f = alpha * k_f_save[37];
    q_f = phi_f * k_f;
    phi_r = sc[0];
    Kc = Kc_save[37];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[37] = q_f - q_r;

    /*reaction 39: 2 H + H2 <=> 2 H2 */
    phi_f = sc[1]*sc[1]*sc[0];
    k_f = k_f_save[38];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[0];
    Kc = Kc_save[38];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[38] = q_f - q_r;

    /*reaction 40: 2 H + H2O <=> H2 + H2O */
    phi_f = sc[1]*sc[1]*sc[5];
    k_f = k_f_save[39];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[5];
    Kc = Kc_save[39];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[39] = q_f - q_r;

    /*reaction 41: 2 H + CO2 <=> H2 + CO2 */
    phi_f = sc[1]*sc[1]*sc[15];
    k_f = k_f_save[40];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[15];
    Kc = Kc_save[40];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[40] = q_f - q_r;

    /*reaction 42: H + OH + M <=> H2O + M */
    phi_f = sc[1]*sc[4];
    alpha = mixture + -0.27*sc[0] + 2.65*sc[5] + sc[13] + 2*sc[26];
    k_f = alpha * k_f_save[41];
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = Kc_save[41];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[41] = q_f - q_r;

    /*reaction 43: H + HO2 <=> O + H2O */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[42];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = Kc_save[42];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[42] = q_f - q_r;

    /*reaction 44: H + HO2 <=> O2 + H2 */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[43];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = Kc_save[43];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[43] = q_f - q_r;

    /*reaction 45: H + HO2 <=> 2 OH */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[44];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = Kc_save[44];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[44] = q_f - q_r;

    /*reaction 46: H + H2O2 <=> HO2 + H2 */
    phi_f = sc[1]*sc[7];
    k_f = k_f_save[45];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = Kc_save[45];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[45] = q_f - q_r;

    /*reaction 47: H + H2O2 <=> OH + H2O */
    phi_f = sc[1]*sc[7];
    k_f = k_f_save[46];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[5];
    Kc = Kc_save[46];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[46] = q_f - q_r;

    /*reaction 48: H + CH <=> C + H2 */
    phi_f = sc[1]*sc[9];
    k_f = k_f_save[47];
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = Kc_save[47];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[47] = q_f - q_r;

    /*reaction 49: H + CH2 (+M) <=> CH3 (+M) */
    phi_f = sc[1]*sc[10];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[48];
    redP = 1e-12 * alpha / k_f * 1.04e+26*exp(-2.76*tc[0]-805.14666528755151376*invT);
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
    Kc = Kc_save[48];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[48] = q_f - q_r;

    /*reaction 50: H + CH2(S) <=> CH + H2 */
    phi_f = sc[1]*sc[11];
    k_f = k_f_save[49];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = Kc_save[49];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[49] = q_f - q_r;

    /*reaction 51: H + CH3 (+M) <=> CH4 (+M) */
    phi_f = sc[1]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + 2*sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[50];
    redP = 1e-12 * alpha / k_f * 2.62e+33*exp(-4.76*tc[0]-1227.8486645635159675*invT);
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
    Kc = Kc_save[50];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[50] = q_f - q_r;

    /*reaction 52: H + CH4 <=> CH3 + H2 */
    phi_f = sc[1]*sc[13];
    k_f = k_f_save[51];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[0];
    Kc = Kc_save[51];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[51] = q_f - q_r;

    /*reaction 53: H + HCO (+M) <=> CH2O (+M) */
    phi_f = sc[1]*sc[16];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[52];
    redP = 1e-12 * alpha / k_f * 2.47e+24*exp(-2.57*tc[0]-213.86708296700587084*invT);
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
    Kc = Kc_save[52];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[52] = q_f - q_r;

    /*reaction 54: H + HCO <=> H2 + CO */
    phi_f = sc[1]*sc[16];
    k_f = k_f_save[53];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[14];
    Kc = Kc_save[53];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[53] = q_f - q_r;

    /*reaction 55: H + CH2O (+M) <=> CH2OH (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[54];
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
    Kc = Kc_save[54];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[54] = q_f - q_r;

    /*reaction 56: H + CH2O (+M) <=> CH3O (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[55];
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
    Kc = Kc_save[55];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[55] = q_f - q_r;

    /*reaction 57: H + CH2O <=> HCO + H2 */
    phi_f = sc[1]*sc[17];
    k_f = k_f_save[56];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[0];
    Kc = Kc_save[56];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[56] = q_f - q_r;

    /*reaction 58: H + CH2OH (+M) <=> CH3OH (+M) */
    phi_f = sc[1]*sc[18];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[57];
    redP = 1e-12 * alpha / k_f * 4.36e+31*exp(-4.65*tc[0]-2556.3406622879761017*invT);
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
    Kc = Kc_save[57];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[57] = q_f - q_r;

    /*reaction 59: H + CH2OH <=> H2 + CH2O */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[58];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[17];
    Kc = Kc_save[58];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[58] = q_f - q_r;

    /*reaction 60: H + CH2OH <=> OH + CH3 */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[59];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = Kc_save[59];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[59] = q_f - q_r;

    /*reaction 61: H + CH2OH <=> CH2(S) + H2O */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[60];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[60];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[60] = q_f - q_r;

    /*reaction 62: H + CH3O (+M) <=> CH3OH (+M) */
    phi_f = sc[1]*sc[19];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[61];
    redP = 1e-12 * alpha / k_f * 4.66e+41*exp(-7.44*tc[0]-7085.2906545304531392*invT);
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
    Kc = Kc_save[61];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[61] = q_f - q_r;

    /*reaction 63: H + CH3O <=> H + CH2OH */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[62];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[18];
    Kc = Kc_save[62];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[62] = q_f - q_r;

    /*reaction 64: H + CH3O <=> H2 + CH2O */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[63];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[17];
    Kc = Kc_save[63];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[63] = q_f - q_r;

    /*reaction 65: H + CH3O <=> OH + CH3 */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[64];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = Kc_save[64];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[64] = q_f - q_r;

    /*reaction 66: H + CH3O <=> CH2(S) + H2O */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[65];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[65];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[65] = q_f - q_r;

    /*reaction 67: H + CH3OH <=> CH2OH + H2 */
    phi_f = sc[1]*sc[20];
    k_f = k_f_save[66];
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[0];
    Kc = Kc_save[66];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[66] = q_f - q_r;

    /*reaction 68: H + CH3OH <=> CH3O + H2 */
    phi_f = sc[1]*sc[20];
    k_f = k_f_save[67];
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[0];
    Kc = Kc_save[67];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[67] = q_f - q_r;

    /*reaction 69: H + C2H (+M) <=> C2H2 (+M) */
    phi_f = sc[1]*sc[21];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[68];
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
    Kc = Kc_save[68];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[68] = q_f - q_r;

    /*reaction 70: H + C2H2 (+M) <=> C2H3 (+M) */
    phi_f = sc[1]*sc[22];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[69];
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
    Kc = Kc_save[69];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[69] = q_f - q_r;

    /*reaction 71: H + C2H3 (+M) <=> C2H4 (+M) */
    phi_f = sc[1]*sc[23];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[70];
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
    Kc = Kc_save[70];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[70] = q_f - q_r;

    /*reaction 72: H + C2H3 <=> H2 + C2H2 */
    phi_f = sc[1]*sc[23];
    k_f = k_f_save[71];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = Kc_save[71];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[71] = q_f - q_r;

    /*reaction 73: H + C2H4 (+M) <=> C2H5 (+M) */
    phi_f = sc[1]*sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[72];
    redP = 1e-12 * alpha / k_f * 6e+41*exp(-7.62*tc[0]-3507.4201606588962932*invT);
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
    Kc = Kc_save[72];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[72] = q_f - q_r;

    /*reaction 74: H + C2H4 <=> C2H3 + H2 */
    phi_f = sc[1]*sc[24];
    k_f = k_f_save[73];
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[0];
    Kc = Kc_save[73];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[73] = q_f - q_r;

    /*reaction 75: H + C2H5 (+M) <=> C2H6 (+M) */
    phi_f = sc[1]*sc[25];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[74];
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
    Kc = Kc_save[74];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[74] = q_f - q_r;

    /*reaction 76: H + C2H5 <=> H2 + C2H4 */
    phi_f = sc[1]*sc[25];
    k_f = k_f_save[75];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[24];
    Kc = Kc_save[75];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[75] = q_f - q_r;

    /*reaction 77: H + C2H6 <=> C2H5 + H2 */
    phi_f = sc[1]*sc[26];
    k_f = k_f_save[76];
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[0];
    Kc = Kc_save[76];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[76] = q_f - q_r;

    /*reaction 78: H + HCCO <=> CH2(S) + CO */
    phi_f = sc[1]*sc[27];
    k_f = k_f_save[77];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[14];
    Kc = Kc_save[77];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[77] = q_f - q_r;

    /*reaction 79: H + CH2CO <=> HCCO + H2 */
    phi_f = sc[1]*sc[28];
    k_f = k_f_save[78];
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[0];
    Kc = Kc_save[78];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[78] = q_f - q_r;

    /*reaction 80: H + CH2CO <=> CH3 + CO */
    phi_f = sc[1]*sc[28];
    k_f = k_f_save[79];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[14];
    Kc = Kc_save[79];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[79] = q_f - q_r;

    /*reaction 81: H + HCCOH <=> H + CH2CO */
    phi_f = sc[1]*sc[29];
    k_f = k_f_save[80];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[80];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[80] = q_f - q_r;

    /*reaction 82: H2 + CO (+M) <=> CH2O (+M) */
    phi_f = sc[0]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[81];
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
    Kc = Kc_save[81];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[81] = q_f - q_r;

    /*reaction 83: OH + H2 <=> H + H2O */
    phi_f = sc[4]*sc[0];
    k_f = k_f_save[82];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    Kc = Kc_save[82];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[82] = q_f - q_r;

    /*reaction 84: 2 OH (+M) <=> H2O2 (+M) */
    phi_f = sc[4]*sc[4];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[83];
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
    Kc = Kc_save[83];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[83] = q_f - q_r;

    /*reaction 85: 2 OH <=> O + H2O */
    phi_f = sc[4]*sc[4];
    k_f = k_f_save[84];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = Kc_save[84];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[84] = q_f - q_r;

    /*reaction 86: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = k_f_save[85];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = Kc_save[85];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[85] = q_f - q_r;

    /*reaction 87: OH + H2O2 <=> HO2 + H2O */
    phi_f = sc[4]*sc[7];
    k_f = k_f_save[86];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
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

    /*reaction 89: OH + C <=> H + CO */
    phi_f = sc[4]*sc[8];
    k_f = k_f_save[88];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[88];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[88] = q_f - q_r;

    /*reaction 90: OH + CH <=> H + HCO */
    phi_f = sc[4]*sc[9];
    k_f = k_f_save[89];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_save[89];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[89] = q_f - q_r;

    /*reaction 91: OH + CH2 <=> H + CH2O */
    phi_f = sc[4]*sc[10];
    k_f = k_f_save[90];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[90];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[90] = q_f - q_r;

    /*reaction 92: OH + CH2 <=> CH + H2O */
    phi_f = sc[4]*sc[10];
    k_f = k_f_save[91];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[5];
    Kc = Kc_save[91];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[91] = q_f - q_r;

    /*reaction 93: OH + CH2(S) <=> H + CH2O */
    phi_f = sc[4]*sc[11];
    k_f = k_f_save[92];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[92];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[92] = q_f - q_r;

    /*reaction 94: OH + CH3 (+M) <=> CH3OH (+M) */
    phi_f = sc[4]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[93];
    redP = 1e-12 * alpha / k_f * 4e+36*exp(-5.92*tc[0]-1580.1003306268198685*invT);
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
    Kc = Kc_save[93];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[93] = q_f - q_r;

    /*reaction 95: OH + CH3 <=> CH2 + H2O */
    phi_f = sc[4]*sc[12];
    k_f = k_f_save[94];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[5];
    Kc = Kc_save[94];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[94] = q_f - q_r;

    /*reaction 96: OH + CH3 <=> CH2(S) + H2O */
    phi_f = sc[4]*sc[12];
    k_f = k_f_save[95];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[95];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[95] = q_f - q_r;

    /*reaction 97: OH + CH4 <=> CH3 + H2O */
    phi_f = sc[4]*sc[13];
    k_f = k_f_save[96];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[5];
    Kc = Kc_save[96];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[96] = q_f - q_r;

    /*reaction 98: OH + CO <=> H + CO2 */
    phi_f = sc[4]*sc[14];
    k_f = k_f_save[97];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[15];
    Kc = Kc_save[97];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[97] = q_f - q_r;

    /*reaction 99: OH + HCO <=> H2O + CO */
    phi_f = sc[4]*sc[16];
    k_f = k_f_save[98];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[14];
    Kc = Kc_save[98];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[98] = q_f - q_r;

    /*reaction 100: OH + CH2O <=> HCO + H2O */
    phi_f = sc[4]*sc[17];
    k_f = k_f_save[99];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[5];
    Kc = Kc_save[99];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[99] = q_f - q_r;

    /*reaction 101: OH + CH2OH <=> H2O + CH2O */
    phi_f = sc[4]*sc[18];
    k_f = k_f_save[100];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[17];
    Kc = Kc_save[100];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[100] = q_f - q_r;

    /*reaction 102: OH + CH3O <=> H2O + CH2O */
    phi_f = sc[4]*sc[19];
    k_f = k_f_save[101];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[17];
    Kc = Kc_save[101];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[101] = q_f - q_r;

    /*reaction 103: OH + CH3OH <=> CH2OH + H2O */
    phi_f = sc[4]*sc[20];
    k_f = k_f_save[102];
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[5];
    Kc = Kc_save[102];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[102] = q_f - q_r;

    /*reaction 104: OH + CH3OH <=> CH3O + H2O */
    phi_f = sc[4]*sc[20];
    k_f = k_f_save[103];
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[5];
    Kc = Kc_save[103];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[103] = q_f - q_r;

    /*reaction 105: OH + C2H <=> H + HCCO */
    phi_f = sc[4]*sc[21];
    k_f = k_f_save[104];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[27];
    Kc = Kc_save[104];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[104] = q_f - q_r;

    /*reaction 106: OH + C2H2 <=> H + CH2CO */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[105];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[105];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[105] = q_f - q_r;

    /*reaction 107: OH + C2H2 <=> H + HCCOH */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[106];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[29];
    Kc = Kc_save[106];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[106] = q_f - q_r;

    /*reaction 108: OH + C2H2 <=> C2H + H2O */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[107];
    q_f = phi_f * k_f;
    phi_r = sc[21]*sc[5];
    Kc = Kc_save[107];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[107] = q_f - q_r;

    /*reaction 109: OH + C2H2 <=> CH3 + CO */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[108];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[14];
    Kc = Kc_save[108];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[108] = q_f - q_r;

    /*reaction 110: OH + C2H3 <=> H2O + C2H2 */
    phi_f = sc[4]*sc[23];
    k_f = k_f_save[109];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[22];
    Kc = Kc_save[109];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[109] = q_f - q_r;

    /*reaction 111: OH + C2H4 <=> C2H3 + H2O */
    phi_f = sc[4]*sc[24];
    k_f = k_f_save[110];
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[5];
    Kc = Kc_save[110];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[110] = q_f - q_r;

    /*reaction 112: OH + C2H6 <=> C2H5 + H2O */
    phi_f = sc[4]*sc[26];
    k_f = k_f_save[111];
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[5];
    Kc = Kc_save[111];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[111] = q_f - q_r;

    /*reaction 113: OH + CH2CO <=> HCCO + H2O */
    phi_f = sc[4]*sc[28];
    k_f = k_f_save[112];
    q_f = phi_f * k_f;
    phi_r = sc[27]*sc[5];
    Kc = Kc_save[112];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[112] = q_f - q_r;

    /*reaction 114: 2 HO2 <=> O2 + H2O2 */
    phi_f = sc[6]*sc[6];
    k_f = k_f_save[113];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[7];
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

    /*reaction 116: HO2 + CH2 <=> OH + CH2O */
    phi_f = sc[6]*sc[10];
    k_f = k_f_save[115];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[115];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[115] = q_f - q_r;

    /*reaction 117: HO2 + CH3 <=> O2 + CH4 */
    phi_f = sc[6]*sc[12];
    k_f = k_f_save[116];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[13];
    Kc = Kc_save[116];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[116] = q_f - q_r;

    /*reaction 118: HO2 + CH3 <=> OH + CH3O */
    phi_f = sc[6]*sc[12];
    k_f = k_f_save[117];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[19];
    Kc = Kc_save[117];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[117] = q_f - q_r;

    /*reaction 119: HO2 + CO <=> OH + CO2 */
    phi_f = sc[6]*sc[14];
    k_f = k_f_save[118];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[15];
    Kc = Kc_save[118];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[118] = q_f - q_r;

    /*reaction 120: HO2 + CH2O <=> HCO + H2O2 */
    phi_f = sc[6]*sc[17];
    k_f = k_f_save[119];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[7];
    Kc = Kc_save[119];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[119] = q_f - q_r;

    /*reaction 121: C + O2 <=> O + CO */
    phi_f = sc[8]*sc[3];
    k_f = k_f_save[120];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[14];
    Kc = Kc_save[120];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[120] = q_f - q_r;

    /*reaction 122: C + CH2 <=> H + C2H */
    phi_f = sc[8]*sc[10];
    k_f = k_f_save[121];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[21];
    Kc = Kc_save[121];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[121] = q_f - q_r;

    /*reaction 123: C + CH3 <=> H + C2H2 */
    phi_f = sc[8]*sc[12];
    k_f = k_f_save[122];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = Kc_save[122];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[122] = q_f - q_r;

    /*reaction 124: CH + O2 <=> O + HCO */
    phi_f = sc[9]*sc[3];
    k_f = k_f_save[123];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[16];
    Kc = Kc_save[123];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[123] = q_f - q_r;

    /*reaction 125: CH + H2 <=> H + CH2 */
    phi_f = sc[9]*sc[0];
    k_f = k_f_save[124];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[10];
    Kc = Kc_save[124];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[124] = q_f - q_r;

    /*reaction 126: CH + H2O <=> H + CH2O */
    phi_f = sc[9]*sc[5];
    k_f = k_f_save[125];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[125];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[125] = q_f - q_r;

    /*reaction 127: CH + CH2 <=> H + C2H2 */
    phi_f = sc[9]*sc[10];
    k_f = k_f_save[126];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = Kc_save[126];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[126] = q_f - q_r;

    /*reaction 128: CH + CH3 <=> H + C2H3 */
    phi_f = sc[9]*sc[12];
    k_f = k_f_save[127];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[23];
    Kc = Kc_save[127];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[127] = q_f - q_r;

    /*reaction 129: CH + CH4 <=> H + C2H4 */
    phi_f = sc[9]*sc[13];
    k_f = k_f_save[128];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = Kc_save[128];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[128] = q_f - q_r;

    /*reaction 130: CH + CO (+M) <=> HCCO (+M) */
    phi_f = sc[9]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[129];
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
    Kc = Kc_save[129];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[129] = q_f - q_r;

    /*reaction 131: CH + CO2 <=> HCO + CO */
    phi_f = sc[9]*sc[15];
    k_f = k_f_save[130];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[14];
    Kc = Kc_save[130];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[130] = q_f - q_r;

    /*reaction 132: CH + CH2O <=> H + CH2CO */
    phi_f = sc[9]*sc[17];
    k_f = k_f_save[131];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[131];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[131] = q_f - q_r;

    /*reaction 133: CH + HCCO <=> CO + C2H2 */
    phi_f = sc[9]*sc[27];
    k_f = k_f_save[132];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[22];
    Kc = Kc_save[132];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[132] = q_f - q_r;

    /*reaction 134: CH2 + O2 => OH + H + CO */
    phi_f = sc[10]*sc[3];
    k_f = k_f_save[133];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[133] = q_f - q_r;

    /*reaction 135: CH2 + H2 <=> H + CH3 */
    phi_f = sc[10]*sc[0];
    k_f = k_f_save[134];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = Kc_save[134];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[134] = q_f - q_r;

    /*reaction 136: 2 CH2 <=> H2 + C2H2 */
    phi_f = sc[10]*sc[10];
    k_f = k_f_save[135];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = Kc_save[135];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[135] = q_f - q_r;

    /*reaction 137: CH2 + CH3 <=> H + C2H4 */
    phi_f = sc[10]*sc[12];
    k_f = k_f_save[136];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = Kc_save[136];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[136] = q_f - q_r;

    /*reaction 138: CH2 + CH4 <=> 2 CH3 */
    phi_f = sc[10]*sc[13];
    k_f = k_f_save[137];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[12];
    Kc = Kc_save[137];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[137] = q_f - q_r;

    /*reaction 139: CH2 + CO (+M) <=> CH2CO (+M) */
    phi_f = sc[10]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[138];
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
    Kc = Kc_save[138];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[138] = q_f - q_r;

    /*reaction 140: CH2 + HCCO <=> C2H3 + CO */
    phi_f = sc[10]*sc[27];
    k_f = k_f_save[139];
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[14];
    Kc = Kc_save[139];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[139] = q_f - q_r;

    /*reaction 141: CH2(S) + N2 <=> CH2 + N2 */
    phi_f = sc[11]*sc[30];
    k_f = k_f_save[140];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[30];
    Kc = Kc_save[140];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[140] = q_f - q_r;

    /*reaction 142: CH2(S) + O2 <=> H + OH + CO */
    phi_f = sc[11]*sc[3];
    k_f = k_f_save[141];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[4]*sc[14];
    Kc = Kc_save[141];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[141] = q_f - q_r;

    /*reaction 143: CH2(S) + O2 <=> CO + H2O */
    phi_f = sc[11]*sc[3];
    k_f = k_f_save[142];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[5];
    Kc = Kc_save[142];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[142] = q_f - q_r;

    /*reaction 144: CH2(S) + H2 <=> CH3 + H */
    phi_f = sc[11]*sc[0];
    k_f = k_f_save[143];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[1];
    Kc = Kc_save[143];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[143] = q_f - q_r;

    /*reaction 145: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    phi_f = sc[11]*sc[5];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[144];
    redP = 1e-12 * alpha / k_f * 1.88e+38*exp(-6.36*tc[0]-2536.2119956557871774*invT);
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
    Kc = Kc_save[144];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[144] = q_f - q_r;

    /*reaction 146: CH2(S) + H2O <=> CH2 + H2O */
    phi_f = sc[11]*sc[5];
    k_f = k_f_save[145];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[5];
    Kc = Kc_save[145];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[145] = q_f - q_r;

    /*reaction 147: CH2(S) + CH3 <=> H + C2H4 */
    phi_f = sc[11]*sc[12];
    k_f = k_f_save[146];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = Kc_save[146];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[146] = q_f - q_r;

    /*reaction 148: CH2(S) + CH4 <=> 2 CH3 */
    phi_f = sc[11]*sc[13];
    k_f = k_f_save[147];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[12];
    Kc = Kc_save[147];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[147] = q_f - q_r;

    /*reaction 149: CH2(S) + CO <=> CH2 + CO */
    phi_f = sc[11]*sc[14];
    k_f = k_f_save[148];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[14];
    Kc = Kc_save[148];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[148] = q_f - q_r;

    /*reaction 150: CH2(S) + CO2 <=> CH2 + CO2 */
    phi_f = sc[11]*sc[15];
    k_f = k_f_save[149];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[15];
    Kc = Kc_save[149];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[149] = q_f - q_r;

    /*reaction 151: CH2(S) + CO2 <=> CO + CH2O */
    phi_f = sc[11]*sc[15];
    k_f = k_f_save[150];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[17];
    Kc = Kc_save[150];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[150] = q_f - q_r;

    /*reaction 152: CH2(S) + C2H6 <=> CH3 + C2H5 */
    phi_f = sc[11]*sc[26];
    k_f = k_f_save[151];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[25];
    Kc = Kc_save[151];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[151] = q_f - q_r;

    /*reaction 153: CH3 + O2 <=> O + CH3O */
    phi_f = sc[12]*sc[3];
    k_f = k_f_save[152];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[19];
    Kc = Kc_save[152];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[152] = q_f - q_r;

    /*reaction 154: CH3 + O2 <=> OH + CH2O */
    phi_f = sc[12]*sc[3];
    k_f = k_f_save[153];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[153];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[153] = q_f - q_r;

    /*reaction 155: CH3 + H2O2 <=> HO2 + CH4 */
    phi_f = sc[12]*sc[7];
    k_f = k_f_save[154];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = Kc_save[154];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[154] = q_f - q_r;

    /*reaction 156: 2 CH3 (+M) <=> C2H6 (+M) */
    phi_f = sc[12]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[155];
    redP = 1e-12 * alpha / k_f * 3.4e+41*exp(-7.03*tc[0]-1389.8844309526357392*invT);
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
    Kc = Kc_save[155];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[155] = q_f - q_r;

    /*reaction 157: 2 CH3 <=> H + C2H5 */
    phi_f = sc[12]*sc[12];
    k_f = k_f_save[156];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[25];
    Kc = Kc_save[156];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[156] = q_f - q_r;

    /*reaction 158: CH3 + HCO <=> CH4 + CO */
    phi_f = sc[12]*sc[16];
    k_f = k_f_save[157];
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[14];
    Kc = Kc_save[157];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[157] = q_f - q_r;

    /*reaction 159: CH3 + CH2O <=> HCO + CH4 */
    phi_f = sc[12]*sc[17];
    k_f = k_f_save[158];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[13];
    Kc = Kc_save[158];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[158] = q_f - q_r;

    /*reaction 160: CH3 + CH3OH <=> CH2OH + CH4 */
    phi_f = sc[12]*sc[20];
    k_f = k_f_save[159];
    q_f = phi_f * k_f;
    phi_r = sc[18]*sc[13];
    Kc = Kc_save[159];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[159] = q_f - q_r;

    /*reaction 161: CH3 + CH3OH <=> CH3O + CH4 */
    phi_f = sc[12]*sc[20];
    k_f = k_f_save[160];
    q_f = phi_f * k_f;
    phi_r = sc[19]*sc[13];
    Kc = Kc_save[160];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[160] = q_f - q_r;

    /*reaction 162: CH3 + C2H4 <=> C2H3 + CH4 */
    phi_f = sc[12]*sc[24];
    k_f = k_f_save[161];
    q_f = phi_f * k_f;
    phi_r = sc[23]*sc[13];
    Kc = Kc_save[161];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[161] = q_f - q_r;

    /*reaction 163: CH3 + C2H6 <=> C2H5 + CH4 */
    phi_f = sc[12]*sc[26];
    k_f = k_f_save[162];
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[13];
    Kc = Kc_save[162];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[162] = q_f - q_r;

    /*reaction 164: HCO + H2O <=> H + CO + H2O */
    phi_f = sc[16]*sc[5];
    k_f = k_f_save[163];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14]*sc[5];
    Kc = Kc_save[163];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[163] = q_f - q_r;

    /*reaction 165: HCO + M <=> H + CO + M */
    phi_f = sc[16];
    alpha = mixture + sc[0] + -1*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = alpha * k_f_save[164];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[164];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[164] = q_f - q_r;

    /*reaction 166: HCO + O2 <=> HO2 + CO */
    phi_f = sc[16]*sc[3];
    k_f = k_f_save[165];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[14];
    Kc = Kc_save[165];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[165] = q_f - q_r;

    /*reaction 167: CH2OH + O2 <=> HO2 + CH2O */
    phi_f = sc[18]*sc[3];
    k_f = k_f_save[166];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = Kc_save[166];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[166] = q_f - q_r;

    /*reaction 168: CH3O + O2 <=> HO2 + CH2O */
    phi_f = sc[19]*sc[3];
    k_f = k_f_save[167];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = Kc_save[167];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[167] = q_f - q_r;

    /*reaction 169: C2H + O2 <=> HCO + CO */
    phi_f = sc[21]*sc[3];
    k_f = k_f_save[168];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[14];
    Kc = Kc_save[168];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[168] = q_f - q_r;

    /*reaction 170: C2H + H2 <=> H + C2H2 */
    phi_f = sc[21]*sc[0];
    k_f = k_f_save[169];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = Kc_save[169];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[169] = q_f - q_r;

    /*reaction 171: C2H3 + O2 <=> HCO + CH2O */
    phi_f = sc[23]*sc[3];
    k_f = k_f_save[170];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[17];
    Kc = Kc_save[170];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[170] = q_f - q_r;

    /*reaction 172: C2H4 (+M) <=> H2 + C2H2 (+M) */
    phi_f = sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[171];
    redP = 1e-6 * alpha / k_f * 1.58e+51*exp(-9.3*tc[0]-49214.589915701588325*invT);
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
    Kc = Kc_save[171];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[171] = q_f - q_r;

    /*reaction 173: C2H5 + O2 <=> HO2 + C2H4 */
    phi_f = sc[25]*sc[3];
    k_f = k_f_save[172];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[24];
    Kc = Kc_save[172];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[172] = q_f - q_r;

    /*reaction 174: HCCO + O2 <=> OH + 2 CO */
    phi_f = sc[27]*sc[3];
    k_f = k_f_save[173];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[14]*sc[14];
    Kc = Kc_save[173];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[173] = q_f - q_r;

    /*reaction 175: 2 HCCO <=> 2 CO + C2H2 */
    phi_f = sc[27]*sc[27];
    k_f = k_f_save[174];
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[14]*sc[22];
    Kc = Kc_save[174];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[174] = q_f - q_r;

    /*reaction 176: O + CH3 => H + H2 + CO */
    phi_f = sc[2]*sc[12];
    k_f = k_f_save[175];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[175] = q_f - q_r;

    /*reaction 177: O + C2H4 <=> H + CH2CHO */
    phi_f = sc[2]*sc[24];
    k_f = k_f_save[176];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[33];
    Kc = Kc_save[176];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[176] = q_f - q_r;

    /*reaction 178: O + C2H5 <=> H + CH3CHO */
    phi_f = sc[2]*sc[25];
    k_f = k_f_save[177];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[34];
    Kc = Kc_save[177];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[177] = q_f - q_r;

    /*reaction 179: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = k_f_save[178];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = Kc_save[178];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[178] = q_f - q_r;

    /*reaction 180: OH + CH3 => H2 + CH2O */
    phi_f = sc[4]*sc[12];
    k_f = k_f_save[179];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[179] = q_f - q_r;

    /*reaction 181: CH + H2 (+M) <=> CH3 (+M) */
    phi_f = sc[9]*sc[0];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[180];
    redP = 1e-12 * alpha / k_f * 4.82e+25*exp(-2.8*tc[0]-296.89783282478464344*invT);
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
    Kc = Kc_save[180];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[180] = q_f - q_r;

    /*reaction 182: CH2 + O2 => 2 H + CO2 */
    phi_f = sc[10]*sc[3];
    k_f = k_f_save[181];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[181] = q_f - q_r;

    /*reaction 183: CH2 + O2 <=> O + CH2O */
    phi_f = sc[10]*sc[3];
    k_f = k_f_save[182];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[17];
    Kc = Kc_save[182];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[182] = q_f - q_r;

    /*reaction 184: CH2 + CH2 => 2 H + C2H2 */
    phi_f = sc[10]*sc[10];
    k_f = k_f_save[183];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[183] = q_f - q_r;

    /*reaction 185: CH2(S) + H2O => H2 + CH2O */
    phi_f = sc[11]*sc[5];
    k_f = k_f_save[184];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[184] = q_f - q_r;

    /*reaction 186: C2H3 + O2 <=> O + CH2CHO */
    phi_f = sc[23]*sc[3];
    k_f = k_f_save[185];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[33];
    Kc = Kc_save[185];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[185] = q_f - q_r;

    /*reaction 187: C2H3 + O2 <=> HO2 + C2H2 */
    phi_f = sc[23]*sc[3];
    k_f = k_f_save[186];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[22];
    Kc = Kc_save[186];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[186] = q_f - q_r;

    /*reaction 188: O + CH3CHO <=> OH + CH2CHO */
    phi_f = sc[2]*sc[34];
    k_f = k_f_save[187];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[33];
    Kc = Kc_save[187];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[187] = q_f - q_r;

    /*reaction 189: O + CH3CHO => OH + CH3 + CO */
    phi_f = sc[2]*sc[34];
    k_f = k_f_save[188];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[188] = q_f - q_r;

    /*reaction 190: O2 + CH3CHO => HO2 + CH3 + CO */
    phi_f = sc[3]*sc[34];
    k_f = k_f_save[189];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[189] = q_f - q_r;

    /*reaction 191: H + CH3CHO <=> CH2CHO + H2 */
    phi_f = sc[1]*sc[34];
    k_f = k_f_save[190];
    q_f = phi_f * k_f;
    phi_r = sc[33]*sc[0];
    Kc = Kc_save[190];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[190] = q_f - q_r;

    /*reaction 192: H + CH3CHO => CH3 + H2 + CO */
    phi_f = sc[1]*sc[34];
    k_f = k_f_save[191];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[191] = q_f - q_r;

    /*reaction 193: OH + CH3CHO => CH3 + H2O + CO */
    phi_f = sc[4]*sc[34];
    k_f = k_f_save[192];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[192] = q_f - q_r;

    /*reaction 194: HO2 + CH3CHO => CH3 + H2O2 + CO */
    phi_f = sc[6]*sc[34];
    k_f = k_f_save[193];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[193] = q_f - q_r;

    /*reaction 195: CH3 + CH3CHO => CH3 + CH4 + CO */
    phi_f = sc[12]*sc[34];
    k_f = k_f_save[194];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[194] = q_f - q_r;

    /*reaction 196: H + CH2CO (+M) <=> CH2CHO (+M) */
    phi_f = sc[1]*sc[28];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[195];
    redP = 1e-12 * alpha / k_f * 1.012e+42*exp(-7.63*tc[0]-1939.3970300113896883*invT);
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
    phi_r = sc[33];
    Kc = Kc_save[195];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[195] = q_f - q_r;

    /*reaction 197: O + CH2CHO => H + CH2 + CO2 */
    phi_f = sc[2]*sc[33];
    k_f = k_f_save[196];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[196] = q_f - q_r;

    /*reaction 198: O2 + CH2CHO => OH + CO + CH2O */
    phi_f = sc[3]*sc[33];
    k_f = k_f_save[197];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[197] = q_f - q_r;

    /*reaction 199: O2 + CH2CHO => OH + 2 HCO */
    phi_f = sc[3]*sc[33];
    k_f = k_f_save[198];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[198] = q_f - q_r;

    /*reaction 200: H + CH2CHO <=> CH3 + HCO */
    phi_f = sc[1]*sc[33];
    k_f = k_f_save[199];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[16];
    Kc = Kc_save[199];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[199] = q_f - q_r;

    /*reaction 201: H + CH2CHO <=> CH2CO + H2 */
    phi_f = sc[1]*sc[33];
    k_f = k_f_save[200];
    q_f = phi_f * k_f;
    phi_r = sc[28]*sc[0];
    Kc = Kc_save[200];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[200] = q_f - q_r;

    /*reaction 202: OH + CH2CHO <=> H2O + CH2CO */
    phi_f = sc[4]*sc[33];
    k_f = k_f_save[201];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[28];
    Kc = Kc_save[201];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[201] = q_f - q_r;

    /*reaction 203: OH + CH2CHO <=> HCO + CH2OH */
    phi_f = sc[4]*sc[33];
    k_f = k_f_save[202];
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[18];
    Kc = Kc_save[202];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[202] = q_f - q_r;

    /*reaction 204: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    phi_f = sc[12]*sc[25];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[203];
    redP = 1e-12 * alpha / k_f * 2.71e+74*exp(-16.82*tc[0]-6574.525738738662767*invT);
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
    phi_r = sc[32];
    Kc = Kc_save[203];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[203] = q_f - q_r;

    /*reaction 205: O + C3H8 <=> OH + C3H7 */
    phi_f = sc[2]*sc[32];
    k_f = k_f_save[204];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[31];
    Kc = Kc_save[204];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[204] = q_f - q_r;

    /*reaction 206: H + C3H8 <=> C3H7 + H2 */
    phi_f = sc[1]*sc[32];
    k_f = k_f_save[205];
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[0];
    Kc = Kc_save[205];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[205] = q_f - q_r;

    /*reaction 207: OH + C3H8 <=> C3H7 + H2O */
    phi_f = sc[4]*sc[32];
    k_f = k_f_save[206];
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[5];
    Kc = Kc_save[206];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[206] = q_f - q_r;

    /*reaction 208: C3H7 + H2O2 <=> HO2 + C3H8 */
    phi_f = sc[31]*sc[7];
    k_f = k_f_save[207];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[32];
    Kc = Kc_save[207];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[207] = q_f - q_r;

    /*reaction 209: CH3 + C3H8 <=> C3H7 + CH4 */
    phi_f = sc[12]*sc[32];
    k_f = k_f_save[208];
    q_f = phi_f * k_f;
    phi_r = sc[31]*sc[13];
    Kc = Kc_save[208];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[208] = q_f - q_r;

    /*reaction 210: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    phi_f = sc[12]*sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[209];
    redP = 1e-12 * alpha / k_f * 3e+63*exp(-14.6*tc[0]-9143.4468176717564347*invT);
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
    phi_r = sc[31];
    Kc = Kc_save[209];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[209] = q_f - q_r;

    /*reaction 211: O + C3H7 <=> C2H5 + CH2O */
    phi_f = sc[2]*sc[31];
    k_f = k_f_save[210];
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[17];
    Kc = Kc_save[210];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[210] = q_f - q_r;

    /*reaction 212: H + C3H7 (+M) <=> C3H8 (+M) */
    phi_f = sc[1]*sc[31];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[211];
    redP = 1e-12 * alpha / k_f * 4.42e+61*exp(-13.545*tc[0]-5715.0316735442011122*invT);
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
    phi_r = sc[32];
    Kc = Kc_save[211];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[211] = q_f - q_r;

    /*reaction 213: H + C3H7 <=> CH3 + C2H5 */
    phi_f = sc[1]*sc[31];
    k_f = k_f_save[212];
    q_f = phi_f * k_f;
    phi_r = sc[12]*sc[25];
    Kc = Kc_save[212];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[212] = q_f - q_r;

    /*reaction 214: OH + C3H7 <=> C2H5 + CH2OH */
    phi_f = sc[4]*sc[31];
    k_f = k_f_save[213];
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[18];
    Kc = Kc_save[213];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[213] = q_f - q_r;

    /*reaction 215: HO2 + C3H7 <=> O2 + C3H8 */
    phi_f = sc[6]*sc[31];
    k_f = k_f_save[214];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[32];
    Kc = Kc_save[214];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[214] = q_f - q_r;

    /*reaction 216: HO2 + C3H7 => OH + C2H5 + CH2O */
    phi_f = sc[6]*sc[31];
    k_f = k_f_save[215];
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[215] = q_f - q_r;

    /*reaction 217: CH3 + C3H7 <=> 2 C2H5 */
    phi_f = sc[12]*sc[31];
    k_f = k_f_save[216];
    q_f = phi_f * k_f;
    phi_r = sc[25]*sc[25];
    Kc = Kc_save[216];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[216] = q_f - q_r;

    return;
}


/*compute the progress rate for each reaction */
void progressRateFR(double * q_f, double * q_r, double * sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[35];                /*Gibbs free energy */
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
    for (id = 0; id < 35; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    if (T != T_save)
    {
        T_save = T;

        k_f_save[0] = 1e-12 * 1.2e+17*exp(-1*tc[0]);
        k_f_save[1] = 1e-12 * 5e+17*exp(-1*tc[0]);
        k_f_save[2] = 1e-06 * 38700*exp(2.7*tc[0]-3150.1363279375455022*invT);
        k_f_save[3] = 1e-06 * 2e+13;
        k_f_save[4] = 1e-06 * 9.63e+06*exp(2*tc[0]-2012.866663218878557*invT);
        k_f_save[5] = 1e-06 * 5.7e+13;
        k_f_save[6] = 1e-06 * 8e+13;
        k_f_save[7] = 1e-06 * 1.5e+13;
        k_f_save[8] = 1e-06 * 1.5e+13;
        k_f_save[9] = 1e-06 * 5.06e+13;
        k_f_save[10] = 1e-06 * 1.02e+09*exp(1.5*tc[0]-4327.6633259205891591*invT);
        k_f_save[11] = 1e-06 * 1.8e+10*exp(-1200.1717479442565946*invT);
        k_f_save[12] = 1e-06 * 3e+13;
        k_f_save[13] = 1e-06 * 3e+13;
        k_f_save[14] = 1e-06 * 3.9e+13*exp(-1781.3869969487077469*invT);
        k_f_save[15] = 1e-06 * 1e+13;
        k_f_save[16] = 1e-06 * 1e+13;
        k_f_save[17] = 1e-06 * 388000*exp(2.5*tc[0]-1559.9716639946311716*invT);
        k_f_save[18] = 1e-06 * 130000*exp(2.5*tc[0]-2516.0833290235987079*invT);
        k_f_save[19] = 1e-06 * 5e+13;
        k_f_save[20] = 1e-06 * 1.35e+07*exp(2*tc[0]-956.11166502896742259*invT);
        k_f_save[21] = 1e-06 * 4.6e+19*exp(-1.41*tc[0]-14568.122475046635373*invT);
        k_f_save[22] = 1e-06 * 6.94e+06*exp(2*tc[0]-956.11166502896742259*invT);
        k_f_save[23] = 1e-06 * 3e+13;
        k_f_save[24] = 1e-06 * 1.25e+07*exp(1.83*tc[0]-110.7076664770383303*invT);
        k_f_save[25] = 1e-06 * 2.24e+13;
        k_f_save[26] = 1e-06 * 8.98e+07*exp(1.92*tc[0]-2863.3028284288548093*invT);
        k_f_save[27] = 1e-06 * 1e+14;
        k_f_save[28] = 1e-06 * 1e+13*exp(-4025.733326437757114*invT);
        k_f_save[29] = 1e-06 * 1.75e+12*exp(-679.34249883637164658*invT);
        k_f_save[30] = 1e-06 * 2.5e+12*exp(-24053.756625465601246*invT);
        k_f_save[31] = 1e-06 * 1e+14*exp(-20128.666632188789663*invT);
        k_f_save[32] = 1e-12 * 2.8e+18*exp(-0.86*tc[0]);
        k_f_save[33] = 1e-12 * 2.08e+19*exp(-1.24*tc[0]);
        k_f_save[34] = 1e-12 * 1.126e+19*exp(-0.76*tc[0]);
        k_f_save[35] = 1e-12 * 2.6e+19*exp(-1.24*tc[0]);
        k_f_save[36] = 1e-06 * 2.65e+16*exp(-0.6707*tc[0]-8575.3152019782282878*invT);
        k_f_save[37] = 1e-12 * 1e+18*exp(-1*tc[0]);
        k_f_save[38] = 1e-12 * 9e+16*exp(-0.6*tc[0]);
        k_f_save[39] = 1e-12 * 6e+19*exp(-1.25*tc[0]);
        k_f_save[40] = 1e-12 * 5.5e+20*exp(-2*tc[0]);
        k_f_save[41] = 1e-12 * 2.2e+22*exp(-2*tc[0]);
        k_f_save[42] = 1e-06 * 3.97e+12*exp(-337.65838275496690812*invT);
        k_f_save[43] = 1e-06 * 4.48e+13*exp(-537.43539907944057177*invT);
        k_f_save[44] = 1e-06 * 8.4e+13*exp(-319.54258278599701271*invT);
        k_f_save[45] = 1e-06 * 1.21e+07*exp(2*tc[0]-2616.7266621845424197*invT);
        k_f_save[46] = 1e-06 * 1e+13*exp(-1811.579996896990906*invT);
        k_f_save[47] = 1e-06 * 1.65e+14;
        k_f_save[48] = 1e-06 * 6e+14;
        k_f_save[49] = 1e-06 * 3e+13;
        k_f_save[50] = 1e-06 * 1.39e+16*exp(-0.534*tc[0]-269.72413287132980031*invT);
        k_f_save[51] = 1e-06 * 6.6e+08*exp(1.62*tc[0]-5454.8686573231616421*invT);
        k_f_save[52] = 1e-06 * 1.09e+12*exp(0.48*tc[0]+130.83633310922709825*invT);
        k_f_save[53] = 1e-06 * 7.34e+13;
        k_f_save[54] = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1811.579996896990906*invT);
        k_f_save[55] = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1308.3633310922712099*invT);
        k_f_save[56] = 1e-06 * 5.74e+07*exp(1.9*tc[0]-1379.8200976365412771*invT);
        k_f_save[57] = 1e-06 * 1.055e+12*exp(0.5*tc[0]-43.276633259205894433*invT);
        k_f_save[58] = 1e-06 * 2e+13;
        k_f_save[59] = 1e-06 * 1.65e+11*exp(0.65*tc[0]+142.9135330885404187*invT);
        k_f_save[60] = 1e-06 * 3.28e+13*exp(-0.09*tc[0]-306.96216614087899188*invT);
        k_f_save[61] = 1e-06 * 2.43e+12*exp(0.515*tc[0]-25.160833290235984805*invT);
        k_f_save[62] = 1e-06 * 4.15e+07*exp(1.63*tc[0]-968.1888650082806862*invT);
        k_f_save[63] = 1e-06 * 2e+13;
        k_f_save[64] = 1e-06 * 1.5e+12*exp(0.5*tc[0]+55.35383323851916515*invT);
        k_f_save[65] = 1e-06 * 2.62e+14*exp(-0.23*tc[0]-538.44183241105008619*invT);
        k_f_save[66] = 1e-06 * 1.7e+07*exp(2.1*tc[0]-2450.6651624689848177*invT);
        k_f_save[67] = 1e-06 * 4.2e+06*exp(2.1*tc[0]-2450.6651624689848177*invT);
        k_f_save[68] = 1e-06 * 1e+17*exp(-1*tc[0]);
        k_f_save[69] = 1e-06 * 5.6e+12*exp(-1207.7199979313272706*invT);
        k_f_save[70] = 1e-06 * 6.08e+12*exp(0.27*tc[0]-140.90066642532150354*invT);
        k_f_save[71] = 1e-06 * 3e+13;
        k_f_save[72] = 1e-06 * 5.4e+11*exp(0.454*tc[0]-915.85433176458980142*invT);
        k_f_save[73] = 1e-06 * 1.325e+06*exp(2.53*tc[0]-6159.3719894497689893*invT);
        k_f_save[74] = 1e-06 * 5.21e+17*exp(-0.99*tc[0]-795.08233197145705162*invT);
        k_f_save[75] = 1e-06 * 2e+12;
        k_f_save[76] = 1e-06 * 1.15e+08*exp(1.9*tc[0]-3789.2214935095394139*invT);
        k_f_save[77] = 1e-06 * 1e+14;
        k_f_save[78] = 1e-06 * 5e+13*exp(-4025.733326437757114*invT);
        k_f_save[79] = 1e-06 * 1.13e+13*exp(-1725.0267303785790318*invT);
        k_f_save[80] = 1e-06 * 1e+13;
        k_f_save[81] = 1e-06 * 4.3e+07*exp(1.5*tc[0]-40056.046598055690993*invT);
        k_f_save[82] = 1e-06 * 2.16e+08*exp(1.51*tc[0]-1726.0331637101885462*invT);
        k_f_save[83] = 1e-06 * 7.4e+13*exp(-0.37*tc[0]);
        k_f_save[84] = 1e-06 * 35700*exp(2.4*tc[0]+1061.7871648479585929*invT);
        k_f_save[85] = 1e-06 * 1.45e+13*exp(+251.60833290235981963*invT);
        k_f_save[86] = 1e-06 * 2e+12*exp(-214.8735162986153*invT);
        k_f_save[87] = 1e-06 * 1.7e+18*exp(-14799.602141316805501*invT);
        k_f_save[88] = 1e-06 * 5e+13;
        k_f_save[89] = 1e-06 * 3e+13;
        k_f_save[90] = 1e-06 * 2e+13;
        k_f_save[91] = 1e-06 * 1.13e+07*exp(2*tc[0]-1509.6499974141590883*invT);
        k_f_save[92] = 1e-06 * 3e+13;
        k_f_save[93] = 1e-06 * 2.79e+18*exp(-1.43*tc[0]-669.27816552027718444*invT);
        k_f_save[94] = 1e-06 * 5.6e+07*exp(1.6*tc[0]-2727.4343286615808211*invT);
        k_f_save[95] = 1e-06 * 6.44e+17*exp(-1.34*tc[0]-713.05801544528776503*invT);
        k_f_save[96] = 1e-06 * 1e+08*exp(1.6*tc[0]-1570.0359973107254064*invT);
        k_f_save[97] = 1e-06 * 4.76e+07*exp(1.228*tc[0]-35.225166606330375885*invT);
        k_f_save[98] = 1e-06 * 5e+13;
        k_f_save[99] = 1e-06 * 3.43e+09*exp(1.18*tc[0]+224.93784961470970529*invT);
        k_f_save[100] = 1e-06 * 5e+12;
        k_f_save[101] = 1e-06 * 5e+12;
        k_f_save[102] = 1e-06 * 1.44e+06*exp(2*tc[0]+422.70199927596451062*invT);
        k_f_save[103] = 1e-06 * 6.3e+06*exp(2*tc[0]-754.82499870707954415*invT);
        k_f_save[104] = 1e-06 * 2e+13;
        k_f_save[105] = 1e-06 * 0.000218*exp(4.5*tc[0]+503.21666580471963925*invT);
        k_f_save[106] = 1e-06 * 504000*exp(2.3*tc[0]-6793.4249883637157836*invT);
        k_f_save[107] = 1e-06 * 3.37e+07*exp(2*tc[0]-7045.0333212660752906*invT);
        k_f_save[108] = 1e-06 * 0.000483*exp(4*tc[0]+1006.4333316094392785*invT);
        k_f_save[109] = 1e-06 * 5e+12;
        k_f_save[110] = 1e-06 * 3.6e+06*exp(2*tc[0]-1258.0416645117993539*invT);
        k_f_save[111] = 1e-06 * 3.54e+06*exp(2.12*tc[0]-437.79849925010609013*invT);
        k_f_save[112] = 1e-06 * 7.5e+12*exp(-1006.4333316094392785*invT);
        k_f_save[113] = 1e-06 * 1.3e+11*exp(+820.24316526169309327*invT);
        k_f_save[114] = 1e-06 * 4.2e+14*exp(-6038.5999896566363532*invT);
        k_f_save[115] = 1e-06 * 2e+13;
        k_f_save[116] = 1e-06 * 1e+12;
        k_f_save[117] = 1e-06 * 3.78e+13;
        k_f_save[118] = 1e-06 * 1.5e+14*exp(-11875.913312991384373*invT);
        k_f_save[119] = 1e-06 * 5.6e+06*exp(2*tc[0]-6038.5999896566363532*invT);
        k_f_save[120] = 1e-06 * 5.8e+13*exp(-289.85279950351855405*invT);
        k_f_save[121] = 1e-06 * 5e+13;
        k_f_save[122] = 1e-06 * 5e+13;
        k_f_save[123] = 1e-06 * 6.71e+13;
        k_f_save[124] = 1e-06 * 1.08e+14*exp(-1565.003830652678289*invT);
        k_f_save[125] = 1e-06 * 5.71e+12*exp(+379.92858268256338761*invT);
        k_f_save[126] = 1e-06 * 4e+13;
        k_f_save[127] = 1e-06 * 3e+13;
        k_f_save[128] = 1e-06 * 6e+13;
        k_f_save[129] = 1e-06 * 5e+13;
        k_f_save[130] = 1e-06 * 1.9e+14*exp(-7946.7975863881329133*invT);
        k_f_save[131] = 1e-06 * 9.46e+13*exp(+259.15658288943063781*invT);
        k_f_save[132] = 1e-06 * 5e+13;
        k_f_save[133] = 1e-06 * 5e+12*exp(-754.82499870707954415*invT);
        k_f_save[134] = 1e-06 * 500000*exp(2*tc[0]-3638.256493768123164*invT);
        k_f_save[135] = 1e-06 * 1.6e+15*exp(-6010.4198563715717682*invT);
        k_f_save[136] = 1e-06 * 4e+13;
        k_f_save[137] = 1e-06 * 2.46e+06*exp(2*tc[0]-4161.6018262050320118*invT);
        k_f_save[138] = 1e-06 * 8.1e+11*exp(0.5*tc[0]-2269.5071627792858635*invT);
        k_f_save[139] = 1e-06 * 3e+13;
        k_f_save[140] = 1e-06 * 1.5e+13*exp(-301.92999948283181766*invT);
        k_f_save[141] = 1e-06 * 2.8e+13;
        k_f_save[142] = 1e-06 * 1.2e+13;
        k_f_save[143] = 1e-06 * 7e+13;
        k_f_save[144] = 1e-06 * 4.82e+17*exp(-1.16*tc[0]-576.18308234640414867*invT);
        k_f_save[145] = 1e-06 * 3e+13;
        k_f_save[146] = 1e-06 * 1.2e+13*exp(+286.83349950869023814*invT);
        k_f_save[147] = 1e-06 * 1.6e+13*exp(+286.83349950869023814*invT);
        k_f_save[148] = 1e-06 * 9e+12;
        k_f_save[149] = 1e-06 * 7e+12;
        k_f_save[150] = 1e-06 * 1.4e+13;
        k_f_save[151] = 1e-06 * 4e+13*exp(+276.76916619259583285*invT);
        k_f_save[152] = 1e-06 * 3.56e+13*exp(-15338.0439737278557*invT);
        k_f_save[153] = 1e-06 * 2.31e+12*exp(-10222.846565822879711*invT);
        k_f_save[154] = 1e-06 * 24500*exp(2.47*tc[0]-2606.6623288684481849*invT);
        k_f_save[155] = 1e-06 * 6.77e+16*exp(-1.18*tc[0]-329.10369943628666078*invT);
        k_f_save[156] = 1e-06 * 6.84e+12*exp(0.1*tc[0]-5334.096657530029006*invT);
        k_f_save[157] = 1e-06 * 2.648e+13;
        k_f_save[158] = 1e-06 * 3320*exp(2.81*tc[0]-2948.8496616156576238*invT);
        k_f_save[159] = 1e-06 * 3e+07*exp(1.5*tc[0]-5001.973658098913802*invT);
        k_f_save[160] = 1e-06 * 1e+07*exp(1.5*tc[0]-5001.973658098913802*invT);
        k_f_save[161] = 1e-06 * 227000*exp(2*tc[0]-4629.5933254034207494*invT);
        k_f_save[162] = 1e-06 * 6.14e+06*exp(1.74*tc[0]-5258.6141576593208811*invT);
        k_f_save[163] = 1e-06 * 1.5e+18*exp(-1*tc[0]-8554.6833186802341515*invT);
        k_f_save[164] = 1e-06 * 1.87e+17*exp(-1*tc[0]-8554.6833186802341515*invT);
        k_f_save[165] = 1e-06 * 1.345e+13*exp(-201.28666632188787844*invT);
        k_f_save[166] = 1e-06 * 1.8e+13*exp(-452.89499922424772649*invT);
        k_f_save[167] = 1e-06 * 4.28e-13*exp(7.6*tc[0]+1776.3548302906606295*invT);
        k_f_save[168] = 1e-06 * 1e+13*exp(+379.92858268256338761*invT);
        k_f_save[169] = 1e-06 * 5.68e+10*exp(0.9*tc[0]-1002.9108149488064328*invT);
        k_f_save[170] = 1e-06 * 4.58e+16*exp(-1.39*tc[0]-510.76491579179048586*invT);
        k_f_save[171] = 1 * 8e+12*exp(0.44*tc[0]-43664.110091875525541*invT);
        k_f_save[172] = 1e-06 * 8.4e+11*exp(-1949.9645799932889076*invT);
        k_f_save[173] = 1e-06 * 3.2e+12*exp(-429.7470325972306*invT);
        k_f_save[174] = 1e-06 * 1e+13;
        k_f_save[175] = 1e-06 * 3.37e+13;
        k_f_save[176] = 1e-06 * 6.7e+06*exp(1.83*tc[0]-110.7076664770383303*invT);
        k_f_save[177] = 1e-06 * 1.096e+14;
        k_f_save[178] = 1e-06 * 5e+15*exp(-8720.7448183957912988*invT);
        k_f_save[179] = 1e-06 * 8e+09*exp(0.5*tc[0]+883.1452484872830837*invT);
        k_f_save[180] = 1e-06 * 1.97e+12*exp(0.43*tc[0]+186.1901663477462705*invT);
        k_f_save[181] = 1e-06 * 5.8e+12*exp(-754.82499870707954415*invT);
        k_f_save[182] = 1e-06 * 2.4e+12*exp(-754.82499870707954415*invT);
        k_f_save[183] = 1e-06 * 2e+14*exp(-5529.8479405280650099*invT);
        k_f_save[184] = 1e-06 * 6.82e+10*exp(0.25*tc[0]+470.50758252741292154*invT);
        k_f_save[185] = 1e-06 * 3.03e+11*exp(0.29*tc[0]-5.5353833238519163373*invT);
        k_f_save[186] = 1e-06 * 1.337e+06*exp(1.61*tc[0]+193.23519966901235989*invT);
        k_f_save[187] = 1e-06 * 5.84e+12*exp(-909.81573177493316962*invT);
        k_f_save[188] = 1e-06 * 5.84e+12*exp(-909.81573177493316962*invT);
        k_f_save[189] = 1e-06 * 3.01e+13*exp(-19700.932466254773317*invT);
        k_f_save[190] = 1e-06 * 2.05e+09*exp(1.16*tc[0]-1210.2360812603508293*invT);
        k_f_save[191] = 1e-06 * 2.05e+09*exp(1.16*tc[0]-1210.2360812603508293*invT);
        k_f_save[192] = 1e-06 * 2.343e+10*exp(0.73*tc[0]+560.08014904065305473*invT);
        k_f_save[193] = 1e-06 * 3.01e+12*exp(-5999.8523063896727763*invT);
        k_f_save[194] = 1e-06 * 2.72e+06*exp(1.77*tc[0]-2979.0426615639403281*invT);
        k_f_save[195] = 1e-06 * 4.865e+11*exp(0.422*tc[0]+883.1452484872830837*invT);
        k_f_save[196] = 1e-06 * 1.5e+14;
        k_f_save[197] = 1e-06 * 1.81e+10;
        k_f_save[198] = 1e-06 * 2.35e+10;
        k_f_save[199] = 1e-06 * 2.2e+13;
        k_f_save[200] = 1e-06 * 1.1e+13;
        k_f_save[201] = 1e-06 * 1.2e+13;
        k_f_save[202] = 1e-06 * 3.01e+13;
        k_f_save[203] = 1e-06 * 9.43e+12;
        k_f_save[204] = 1e-06 * 193000*exp(2.68*tc[0]-1869.9531301303384225*invT);
        k_f_save[205] = 1e-06 * 1.32e+06*exp(2.54*tc[0]-3399.7317941766864351*invT);
        k_f_save[206] = 1e-06 * 3.16e+07*exp(1.8*tc[0]-470.00436586160816432*invT);
        k_f_save[207] = 1e-06 * 378*exp(2.72*tc[0]-754.82499870707954415*invT);
        k_f_save[208] = 1e-06 * 0.903*exp(3.65*tc[0]-3600.0120271669643444*invT);
        k_f_save[209] = 1e-06 * 2.55e+06*exp(1.6*tc[0]-2868.3349950869019267*invT);
        k_f_save[210] = 1e-06 * 9.64e+13;
        k_f_save[211] = 1e-06 * 3.613e+13;
        k_f_save[212] = 1e-06 * 4.06e+06*exp(2.19*tc[0]-447.86283256620055226*invT);
        k_f_save[213] = 1e-06 * 2.41e+13;
        k_f_save[214] = 1e-06 * 2.55e+10*exp(0.255*tc[0]+474.53331585385063818*invT);
        k_f_save[215] = 1e-06 * 2.41e+13;
        k_f_save[216] = 1e-06 * 1.927e+13*exp(-0.32*tc[0]);

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
        Kc_save[36] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));
        Kc_save[37] = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));
        Kc_save[38] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));
        Kc_save[39] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));
        Kc_save[40] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[15]) - (g_RT[0] + g_RT[15]));
        Kc_save[41] = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));
        Kc_save[42] = exp((g_RT[1] + g_RT[6]) - (g_RT[2] + g_RT[5]));
        Kc_save[43] = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));
        Kc_save[44] = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));
        Kc_save[45] = exp((g_RT[1] + g_RT[7]) - (g_RT[6] + g_RT[0]));
        Kc_save[46] = exp((g_RT[1] + g_RT[7]) - (g_RT[4] + g_RT[5]));
        Kc_save[47] = exp((g_RT[1] + g_RT[9]) - (g_RT[8] + g_RT[0]));
        Kc_save[48] = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[12]));
        Kc_save[49] = exp((g_RT[1] + g_RT[11]) - (g_RT[9] + g_RT[0]));
        Kc_save[50] = 1.0 / (refC) * exp((g_RT[1] + g_RT[12]) - (g_RT[13]));
        Kc_save[51] = exp((g_RT[1] + g_RT[13]) - (g_RT[12] + g_RT[0]));
        Kc_save[52] = 1.0 / (refC) * exp((g_RT[1] + g_RT[16]) - (g_RT[17]));
        Kc_save[53] = exp((g_RT[1] + g_RT[16]) - (g_RT[0] + g_RT[14]));
        Kc_save[54] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[18]));
        Kc_save[55] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[19]));
        Kc_save[56] = exp((g_RT[1] + g_RT[17]) - (g_RT[16] + g_RT[0]));
        Kc_save[57] = 1.0 / (refC) * exp((g_RT[1] + g_RT[18]) - (g_RT[20]));
        Kc_save[58] = exp((g_RT[1] + g_RT[18]) - (g_RT[0] + g_RT[17]));
        Kc_save[59] = exp((g_RT[1] + g_RT[18]) - (g_RT[4] + g_RT[12]));
        Kc_save[60] = exp((g_RT[1] + g_RT[18]) - (g_RT[11] + g_RT[5]));
        Kc_save[61] = 1.0 / (refC) * exp((g_RT[1] + g_RT[19]) - (g_RT[20]));
        Kc_save[62] = exp((g_RT[1] + g_RT[19]) - (g_RT[1] + g_RT[18]));
        Kc_save[63] = exp((g_RT[1] + g_RT[19]) - (g_RT[0] + g_RT[17]));
        Kc_save[64] = exp((g_RT[1] + g_RT[19]) - (g_RT[4] + g_RT[12]));
        Kc_save[65] = exp((g_RT[1] + g_RT[19]) - (g_RT[11] + g_RT[5]));
        Kc_save[66] = exp((g_RT[1] + g_RT[20]) - (g_RT[18] + g_RT[0]));
        Kc_save[67] = exp((g_RT[1] + g_RT[20]) - (g_RT[19] + g_RT[0]));
        Kc_save[68] = 1.0 / (refC) * exp((g_RT[1] + g_RT[21]) - (g_RT[22]));
        Kc_save[69] = 1.0 / (refC) * exp((g_RT[1] + g_RT[22]) - (g_RT[23]));
        Kc_save[70] = 1.0 / (refC) * exp((g_RT[1] + g_RT[23]) - (g_RT[24]));
        Kc_save[71] = exp((g_RT[1] + g_RT[23]) - (g_RT[0] + g_RT[22]));
        Kc_save[72] = 1.0 / (refC) * exp((g_RT[1] + g_RT[24]) - (g_RT[25]));
        Kc_save[73] = exp((g_RT[1] + g_RT[24]) - (g_RT[23] + g_RT[0]));
        Kc_save[74] = 1.0 / (refC) * exp((g_RT[1] + g_RT[25]) - (g_RT[26]));
        Kc_save[75] = exp((g_RT[1] + g_RT[25]) - (g_RT[0] + g_RT[24]));
        Kc_save[76] = exp((g_RT[1] + g_RT[26]) - (g_RT[25] + g_RT[0]));
        Kc_save[77] = exp((g_RT[1] + g_RT[27]) - (g_RT[11] + g_RT[14]));
        Kc_save[78] = exp((g_RT[1] + g_RT[28]) - (g_RT[27] + g_RT[0]));
        Kc_save[79] = exp((g_RT[1] + g_RT[28]) - (g_RT[12] + g_RT[14]));
        Kc_save[80] = exp((g_RT[1] + g_RT[29]) - (g_RT[1] + g_RT[28]));
        Kc_save[81] = 1.0 / (refC) * exp((g_RT[0] + g_RT[14]) - (g_RT[17]));
        Kc_save[82] = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));
        Kc_save[83] = 1.0 / (refC) * exp((2 * g_RT[4]) - (g_RT[7]));
        Kc_save[84] = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));
        Kc_save[85] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
        Kc_save[86] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
        Kc_save[87] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
        Kc_save[88] = exp((g_RT[4] + g_RT[8]) - (g_RT[1] + g_RT[14]));
        Kc_save[89] = exp((g_RT[4] + g_RT[9]) - (g_RT[1] + g_RT[16]));
        Kc_save[90] = exp((g_RT[4] + g_RT[10]) - (g_RT[1] + g_RT[17]));
        Kc_save[91] = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[5]));
        Kc_save[92] = exp((g_RT[4] + g_RT[11]) - (g_RT[1] + g_RT[17]));
        Kc_save[93] = 1.0 / (refC) * exp((g_RT[4] + g_RT[12]) - (g_RT[20]));
        Kc_save[94] = exp((g_RT[4] + g_RT[12]) - (g_RT[10] + g_RT[5]));
        Kc_save[95] = exp((g_RT[4] + g_RT[12]) - (g_RT[11] + g_RT[5]));
        Kc_save[96] = exp((g_RT[4] + g_RT[13]) - (g_RT[12] + g_RT[5]));
        Kc_save[97] = exp((g_RT[4] + g_RT[14]) - (g_RT[1] + g_RT[15]));
        Kc_save[98] = exp((g_RT[4] + g_RT[16]) - (g_RT[5] + g_RT[14]));
        Kc_save[99] = exp((g_RT[4] + g_RT[17]) - (g_RT[16] + g_RT[5]));
        Kc_save[100] = exp((g_RT[4] + g_RT[18]) - (g_RT[5] + g_RT[17]));
        Kc_save[101] = exp((g_RT[4] + g_RT[19]) - (g_RT[5] + g_RT[17]));
        Kc_save[102] = exp((g_RT[4] + g_RT[20]) - (g_RT[18] + g_RT[5]));
        Kc_save[103] = exp((g_RT[4] + g_RT[20]) - (g_RT[19] + g_RT[5]));
        Kc_save[104] = exp((g_RT[4] + g_RT[21]) - (g_RT[1] + g_RT[27]));
        Kc_save[105] = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[28]));
        Kc_save[106] = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[29]));
        Kc_save[107] = exp((g_RT[4] + g_RT[22]) - (g_RT[21] + g_RT[5]));
        Kc_save[108] = exp((g_RT[4] + g_RT[22]) - (g_RT[12] + g_RT[14]));
        Kc_save[109] = exp((g_RT[4] + g_RT[23]) - (g_RT[5] + g_RT[22]));
        Kc_save[110] = exp((g_RT[4] + g_RT[24]) - (g_RT[23] + g_RT[5]));
        Kc_save[111] = exp((g_RT[4] + g_RT[26]) - (g_RT[25] + g_RT[5]));
        Kc_save[112] = exp((g_RT[4] + g_RT[28]) - (g_RT[27] + g_RT[5]));
        Kc_save[113] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
        Kc_save[114] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
        Kc_save[115] = exp((g_RT[6] + g_RT[10]) - (g_RT[4] + g_RT[17]));
        Kc_save[116] = exp((g_RT[6] + g_RT[12]) - (g_RT[3] + g_RT[13]));
        Kc_save[117] = exp((g_RT[6] + g_RT[12]) - (g_RT[4] + g_RT[19]));
        Kc_save[118] = exp((g_RT[6] + g_RT[14]) - (g_RT[4] + g_RT[15]));
        Kc_save[119] = exp((g_RT[6] + g_RT[17]) - (g_RT[16] + g_RT[7]));
        Kc_save[120] = exp((g_RT[8] + g_RT[3]) - (g_RT[2] + g_RT[14]));
        Kc_save[121] = exp((g_RT[8] + g_RT[10]) - (g_RT[1] + g_RT[21]));
        Kc_save[122] = exp((g_RT[8] + g_RT[12]) - (g_RT[1] + g_RT[22]));
        Kc_save[123] = exp((g_RT[9] + g_RT[3]) - (g_RT[2] + g_RT[16]));
        Kc_save[124] = exp((g_RT[9] + g_RT[0]) - (g_RT[1] + g_RT[10]));
        Kc_save[125] = exp((g_RT[9] + g_RT[5]) - (g_RT[1] + g_RT[17]));
        Kc_save[126] = exp((g_RT[9] + g_RT[10]) - (g_RT[1] + g_RT[22]));
        Kc_save[127] = exp((g_RT[9] + g_RT[12]) - (g_RT[1] + g_RT[23]));
        Kc_save[128] = exp((g_RT[9] + g_RT[13]) - (g_RT[1] + g_RT[24]));
        Kc_save[129] = 1.0 / (refC) * exp((g_RT[9] + g_RT[14]) - (g_RT[27]));
        Kc_save[130] = exp((g_RT[9] + g_RT[15]) - (g_RT[16] + g_RT[14]));
        Kc_save[131] = exp((g_RT[9] + g_RT[17]) - (g_RT[1] + g_RT[28]));
        Kc_save[132] = exp((g_RT[9] + g_RT[27]) - (g_RT[14] + g_RT[22]));
        Kc_save[133] = refC * exp((g_RT[10] + g_RT[3]) - (g_RT[4] + g_RT[1] + g_RT[14]));
        Kc_save[134] = exp((g_RT[10] + g_RT[0]) - (g_RT[1] + g_RT[12]));
        Kc_save[135] = exp((2 * g_RT[10]) - (g_RT[0] + g_RT[22]));
        Kc_save[136] = exp((g_RT[10] + g_RT[12]) - (g_RT[1] + g_RT[24]));
        Kc_save[137] = exp((g_RT[10] + g_RT[13]) - (2 * g_RT[12]));
        Kc_save[138] = 1.0 / (refC) * exp((g_RT[10] + g_RT[14]) - (g_RT[28]));
        Kc_save[139] = exp((g_RT[10] + g_RT[27]) - (g_RT[23] + g_RT[14]));
        Kc_save[140] = exp((g_RT[11] + g_RT[30]) - (g_RT[10] + g_RT[30]));
        Kc_save[141] = refC * exp((g_RT[11] + g_RT[3]) - (g_RT[1] + g_RT[4] + g_RT[14]));
        Kc_save[142] = exp((g_RT[11] + g_RT[3]) - (g_RT[14] + g_RT[5]));
        Kc_save[143] = exp((g_RT[11] + g_RT[0]) - (g_RT[12] + g_RT[1]));
        Kc_save[144] = 1.0 / (refC) * exp((g_RT[11] + g_RT[5]) - (g_RT[20]));
        Kc_save[145] = exp((g_RT[11] + g_RT[5]) - (g_RT[10] + g_RT[5]));
        Kc_save[146] = exp((g_RT[11] + g_RT[12]) - (g_RT[1] + g_RT[24]));
        Kc_save[147] = exp((g_RT[11] + g_RT[13]) - (2 * g_RT[12]));
        Kc_save[148] = exp((g_RT[11] + g_RT[14]) - (g_RT[10] + g_RT[14]));
        Kc_save[149] = exp((g_RT[11] + g_RT[15]) - (g_RT[10] + g_RT[15]));
        Kc_save[150] = exp((g_RT[11] + g_RT[15]) - (g_RT[14] + g_RT[17]));
        Kc_save[151] = exp((g_RT[11] + g_RT[26]) - (g_RT[12] + g_RT[25]));
        Kc_save[152] = exp((g_RT[12] + g_RT[3]) - (g_RT[2] + g_RT[19]));
        Kc_save[153] = exp((g_RT[12] + g_RT[3]) - (g_RT[4] + g_RT[17]));
        Kc_save[154] = exp((g_RT[12] + g_RT[7]) - (g_RT[6] + g_RT[13]));
        Kc_save[155] = 1.0 / (refC) * exp((2 * g_RT[12]) - (g_RT[26]));
        Kc_save[156] = exp((2 * g_RT[12]) - (g_RT[1] + g_RT[25]));
        Kc_save[157] = exp((g_RT[12] + g_RT[16]) - (g_RT[13] + g_RT[14]));
        Kc_save[158] = exp((g_RT[12] + g_RT[17]) - (g_RT[16] + g_RT[13]));
        Kc_save[159] = exp((g_RT[12] + g_RT[20]) - (g_RT[18] + g_RT[13]));
        Kc_save[160] = exp((g_RT[12] + g_RT[20]) - (g_RT[19] + g_RT[13]));
        Kc_save[161] = exp((g_RT[12] + g_RT[24]) - (g_RT[23] + g_RT[13]));
        Kc_save[162] = exp((g_RT[12] + g_RT[26]) - (g_RT[25] + g_RT[13]));
        Kc_save[163] = refC * exp((g_RT[16] + g_RT[5]) - (g_RT[1] + g_RT[14] + g_RT[5]));
        Kc_save[164] = refC * exp((g_RT[16]) - (g_RT[1] + g_RT[14]));
        Kc_save[165] = exp((g_RT[16] + g_RT[3]) - (g_RT[6] + g_RT[14]));
        Kc_save[166] = exp((g_RT[18] + g_RT[3]) - (g_RT[6] + g_RT[17]));
        Kc_save[167] = exp((g_RT[19] + g_RT[3]) - (g_RT[6] + g_RT[17]));
        Kc_save[168] = exp((g_RT[21] + g_RT[3]) - (g_RT[16] + g_RT[14]));
        Kc_save[169] = exp((g_RT[21] + g_RT[0]) - (g_RT[1] + g_RT[22]));
        Kc_save[170] = exp((g_RT[23] + g_RT[3]) - (g_RT[16] + g_RT[17]));
        Kc_save[171] = refC * exp((g_RT[24]) - (g_RT[0] + g_RT[22]));
        Kc_save[172] = exp((g_RT[25] + g_RT[3]) - (g_RT[6] + g_RT[24]));
        Kc_save[173] = refC * exp((g_RT[27] + g_RT[3]) - (g_RT[4] + 2 * g_RT[14]));
        Kc_save[174] = refC * exp((2 * g_RT[27]) - (2 * g_RT[14] + g_RT[22]));
        Kc_save[175] = refC * exp((g_RT[2] + g_RT[12]) - (g_RT[1] + g_RT[0] + g_RT[14]));
        Kc_save[176] = exp((g_RT[2] + g_RT[24]) - (g_RT[1] + g_RT[33]));
        Kc_save[177] = exp((g_RT[2] + g_RT[25]) - (g_RT[1] + g_RT[34]));
        Kc_save[178] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
        Kc_save[179] = exp((g_RT[4] + g_RT[12]) - (g_RT[0] + g_RT[17]));
        Kc_save[180] = 1.0 / (refC) * exp((g_RT[9] + g_RT[0]) - (g_RT[12]));
        Kc_save[181] = refC * exp((g_RT[10] + g_RT[3]) - (2 * g_RT[1] + g_RT[15]));
        Kc_save[182] = exp((g_RT[10] + g_RT[3]) - (g_RT[2] + g_RT[17]));
        Kc_save[183] = refC * exp((g_RT[10] + g_RT[10]) - (2 * g_RT[1] + g_RT[22]));
        Kc_save[184] = exp((g_RT[11] + g_RT[5]) - (g_RT[0] + g_RT[17]));
        Kc_save[185] = exp((g_RT[23] + g_RT[3]) - (g_RT[2] + g_RT[33]));
        Kc_save[186] = exp((g_RT[23] + g_RT[3]) - (g_RT[6] + g_RT[22]));
        Kc_save[187] = exp((g_RT[2] + g_RT[34]) - (g_RT[4] + g_RT[33]));
        Kc_save[188] = refC * exp((g_RT[2] + g_RT[34]) - (g_RT[4] + g_RT[12] + g_RT[14]));
        Kc_save[189] = refC * exp((g_RT[3] + g_RT[34]) - (g_RT[6] + g_RT[12] + g_RT[14]));
        Kc_save[190] = exp((g_RT[1] + g_RT[34]) - (g_RT[33] + g_RT[0]));
        Kc_save[191] = refC * exp((g_RT[1] + g_RT[34]) - (g_RT[12] + g_RT[0] + g_RT[14]));
        Kc_save[192] = refC * exp((g_RT[4] + g_RT[34]) - (g_RT[12] + g_RT[5] + g_RT[14]));
        Kc_save[193] = refC * exp((g_RT[6] + g_RT[34]) - (g_RT[12] + g_RT[7] + g_RT[14]));
        Kc_save[194] = refC * exp((g_RT[12] + g_RT[34]) - (g_RT[12] + g_RT[13] + g_RT[14]));
        Kc_save[195] = 1.0 / (refC) * exp((g_RT[1] + g_RT[28]) - (g_RT[33]));
        Kc_save[196] = refC * exp((g_RT[2] + g_RT[33]) - (g_RT[1] + g_RT[10] + g_RT[15]));
        Kc_save[197] = refC * exp((g_RT[3] + g_RT[33]) - (g_RT[4] + g_RT[14] + g_RT[17]));
        Kc_save[198] = refC * exp((g_RT[3] + g_RT[33]) - (g_RT[4] + 2 * g_RT[16]));
        Kc_save[199] = exp((g_RT[1] + g_RT[33]) - (g_RT[12] + g_RT[16]));
        Kc_save[200] = exp((g_RT[1] + g_RT[33]) - (g_RT[28] + g_RT[0]));
        Kc_save[201] = exp((g_RT[4] + g_RT[33]) - (g_RT[5] + g_RT[28]));
        Kc_save[202] = exp((g_RT[4] + g_RT[33]) - (g_RT[16] + g_RT[18]));
        Kc_save[203] = 1.0 / (refC) * exp((g_RT[12] + g_RT[25]) - (g_RT[32]));
        Kc_save[204] = exp((g_RT[2] + g_RT[32]) - (g_RT[4] + g_RT[31]));
        Kc_save[205] = exp((g_RT[1] + g_RT[32]) - (g_RT[31] + g_RT[0]));
        Kc_save[206] = exp((g_RT[4] + g_RT[32]) - (g_RT[31] + g_RT[5]));
        Kc_save[207] = exp((g_RT[31] + g_RT[7]) - (g_RT[6] + g_RT[32]));
        Kc_save[208] = exp((g_RT[12] + g_RT[32]) - (g_RT[31] + g_RT[13]));
        Kc_save[209] = 1.0 / (refC) * exp((g_RT[12] + g_RT[24]) - (g_RT[31]));
        Kc_save[210] = exp((g_RT[2] + g_RT[31]) - (g_RT[25] + g_RT[17]));
        Kc_save[211] = 1.0 / (refC) * exp((g_RT[1] + g_RT[31]) - (g_RT[32]));
        Kc_save[212] = exp((g_RT[1] + g_RT[31]) - (g_RT[12] + g_RT[25]));
        Kc_save[213] = exp((g_RT[4] + g_RT[31]) - (g_RT[25] + g_RT[18]));
        Kc_save[214] = exp((g_RT[6] + g_RT[31]) - (g_RT[3] + g_RT[32]));
        Kc_save[215] = refC * exp((g_RT[6] + g_RT[31]) - (g_RT[4] + g_RT[25] + g_RT[17]));
        Kc_save[216] = exp((g_RT[12] + g_RT[31]) - (2 * g_RT[25]));
    }

    /*reaction 1: 2 O + M <=> O2 + M */
    phi_f = sc[2]*sc[2];
    alpha = mixture + 1.4*sc[0] + 14.4*sc[5] + sc[13] + 0.75*sc[14] + 2.6*sc[15] + 2*sc[26];
    k_f = alpha * k_f_save[0];
    q_f[0] = phi_f * k_f;
    phi_r = sc[3];
    Kc = Kc_save[0];
    k_r = k_f / Kc;
    q_r[0] = phi_r * k_r;

    /*reaction 2: O + H + M <=> OH + M */
    phi_f = sc[2]*sc[1];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
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

    /*reaction 12: O + CO (+M) <=> CO2 (+M) */
    phi_f = sc[2]*sc[14];
    alpha = mixture + sc[0] + 5*sc[3] + 5*sc[5] + sc[13] + 0.5*sc[14] + 2.5*sc[15] + 2*sc[26];
    k_f = k_f_save[11];
    redP = 1e-12 * alpha / k_f * 6.02e+14*exp(-1509.6499974141590883*invT);
    F = redP / (1 + redP);
    k_f *= F;
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
    alpha = mixture + -1*sc[3] + -1*sc[5] + -0.25*sc[14] + 0.5*sc[15] + 0.5*sc[26] + -1*sc[30];
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

    /*reaction 37: H + O2 <=> O + OH */
    phi_f = sc[1]*sc[3];
    k_f = k_f_save[36];
    q_f[36] = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = Kc_save[36];
    k_r = k_f / Kc;
    q_r[36] = phi_r * k_r;

    /*reaction 38: 2 H + M <=> H2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[0] + -1*sc[5] + sc[13] + -1*sc[15] + 2*sc[26];
    k_f = alpha * k_f_save[37];
    q_f[37] = phi_f * k_f;
    phi_r = sc[0];
    Kc = Kc_save[37];
    k_r = k_f / Kc;
    q_r[37] = phi_r * k_r;

    /*reaction 39: 2 H + H2 <=> 2 H2 */
    phi_f = sc[1]*sc[1]*sc[0];
    k_f = k_f_save[38];
    q_f[38] = phi_f * k_f;
    phi_r = sc[0]*sc[0];
    Kc = Kc_save[38];
    k_r = k_f / Kc;
    q_r[38] = phi_r * k_r;

    /*reaction 40: 2 H + H2O <=> H2 + H2O */
    phi_f = sc[1]*sc[1]*sc[5];
    k_f = k_f_save[39];
    q_f[39] = phi_f * k_f;
    phi_r = sc[0]*sc[5];
    Kc = Kc_save[39];
    k_r = k_f / Kc;
    q_r[39] = phi_r * k_r;

    /*reaction 41: 2 H + CO2 <=> H2 + CO2 */
    phi_f = sc[1]*sc[1]*sc[15];
    k_f = k_f_save[40];
    q_f[40] = phi_f * k_f;
    phi_r = sc[0]*sc[15];
    Kc = Kc_save[40];
    k_r = k_f / Kc;
    q_r[40] = phi_r * k_r;

    /*reaction 42: H + OH + M <=> H2O + M */
    phi_f = sc[1]*sc[4];
    alpha = mixture + -0.27*sc[0] + 2.65*sc[5] + sc[13] + 2*sc[26];
    k_f = alpha * k_f_save[41];
    q_f[41] = phi_f * k_f;
    phi_r = sc[5];
    Kc = Kc_save[41];
    k_r = k_f / Kc;
    q_r[41] = phi_r * k_r;

    /*reaction 43: H + HO2 <=> O + H2O */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[42];
    q_f[42] = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = Kc_save[42];
    k_r = k_f / Kc;
    q_r[42] = phi_r * k_r;

    /*reaction 44: H + HO2 <=> O2 + H2 */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[43];
    q_f[43] = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = Kc_save[43];
    k_r = k_f / Kc;
    q_r[43] = phi_r * k_r;

    /*reaction 45: H + HO2 <=> 2 OH */
    phi_f = sc[1]*sc[6];
    k_f = k_f_save[44];
    q_f[44] = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = Kc_save[44];
    k_r = k_f / Kc;
    q_r[44] = phi_r * k_r;

    /*reaction 46: H + H2O2 <=> HO2 + H2 */
    phi_f = sc[1]*sc[7];
    k_f = k_f_save[45];
    q_f[45] = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = Kc_save[45];
    k_r = k_f / Kc;
    q_r[45] = phi_r * k_r;

    /*reaction 47: H + H2O2 <=> OH + H2O */
    phi_f = sc[1]*sc[7];
    k_f = k_f_save[46];
    q_f[46] = phi_f * k_f;
    phi_r = sc[4]*sc[5];
    Kc = Kc_save[46];
    k_r = k_f / Kc;
    q_r[46] = phi_r * k_r;

    /*reaction 48: H + CH <=> C + H2 */
    phi_f = sc[1]*sc[9];
    k_f = k_f_save[47];
    q_f[47] = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = Kc_save[47];
    k_r = k_f / Kc;
    q_r[47] = phi_r * k_r;

    /*reaction 49: H + CH2 (+M) <=> CH3 (+M) */
    phi_f = sc[1]*sc[10];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[48];
    redP = 1e-12 * alpha / k_f * 1.04e+26*exp(-2.76*tc[0]-805.14666528755151376*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.438*exp(T/-91))+ (0.562*exp(T/-5836))+ (exp(-8552/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[48] = phi_f * k_f;
    phi_r = sc[12];
    Kc = Kc_save[48];
    k_r = k_f / Kc;
    q_r[48] = phi_r * k_r;

    /*reaction 50: H + CH2(S) <=> CH + H2 */
    phi_f = sc[1]*sc[11];
    k_f = k_f_save[49];
    q_f[49] = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = Kc_save[49];
    k_r = k_f / Kc;
    q_r[49] = phi_r * k_r;

    /*reaction 51: H + CH3 (+M) <=> CH4 (+M) */
    phi_f = sc[1]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + 2*sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[50];
    redP = 1e-12 * alpha / k_f * 2.62e+33*exp(-4.76*tc[0]-1227.8486645635159675*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.217*exp(T/-74))+ (0.783*exp(T/-2941))+ (exp(-6964/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[50] = phi_f * k_f;
    phi_r = sc[13];
    Kc = Kc_save[50];
    k_r = k_f / Kc;
    q_r[50] = phi_r * k_r;

    /*reaction 52: H + CH4 <=> CH3 + H2 */
    phi_f = sc[1]*sc[13];
    k_f = k_f_save[51];
    q_f[51] = phi_f * k_f;
    phi_r = sc[12]*sc[0];
    Kc = Kc_save[51];
    k_r = k_f / Kc;
    q_r[51] = phi_r * k_r;

    /*reaction 53: H + HCO (+M) <=> CH2O (+M) */
    phi_f = sc[1]*sc[16];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[52];
    redP = 1e-12 * alpha / k_f * 2.47e+24*exp(-2.57*tc[0]-213.86708296700587084*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2176*exp(T/-271))+ (0.7824*exp(T/-2755))+ (exp(-6570/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[52] = phi_f * k_f;
    phi_r = sc[17];
    Kc = Kc_save[52];
    k_r = k_f / Kc;
    q_r[52] = phi_r * k_r;

    /*reaction 54: H + HCO <=> H2 + CO */
    phi_f = sc[1]*sc[16];
    k_f = k_f_save[53];
    q_f[53] = phi_f * k_f;
    phi_r = sc[0]*sc[14];
    Kc = Kc_save[53];
    k_r = k_f / Kc;
    q_r[53] = phi_r * k_r;

    /*reaction 55: H + CH2O (+M) <=> CH2OH (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[54];
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
    q_f[54] = phi_f * k_f;
    phi_r = sc[18];
    Kc = Kc_save[54];
    k_r = k_f / Kc;
    q_r[54] = phi_r * k_r;

    /*reaction 56: H + CH2O (+M) <=> CH3O (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[55];
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
    q_f[55] = phi_f * k_f;
    phi_r = sc[19];
    Kc = Kc_save[55];
    k_r = k_f / Kc;
    q_r[55] = phi_r * k_r;

    /*reaction 57: H + CH2O <=> HCO + H2 */
    phi_f = sc[1]*sc[17];
    k_f = k_f_save[56];
    q_f[56] = phi_f * k_f;
    phi_r = sc[16]*sc[0];
    Kc = Kc_save[56];
    k_r = k_f / Kc;
    q_r[56] = phi_r * k_r;

    /*reaction 58: H + CH2OH (+M) <=> CH3OH (+M) */
    phi_f = sc[1]*sc[18];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[57];
    redP = 1e-12 * alpha / k_f * 4.36e+31*exp(-4.65*tc[0]-2556.3406622879761017*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.4*exp(T/-100))+ (0.6*exp(T/-90000))+ (exp(-10000/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[57] = phi_f * k_f;
    phi_r = sc[20];
    Kc = Kc_save[57];
    k_r = k_f / Kc;
    q_r[57] = phi_r * k_r;

    /*reaction 59: H + CH2OH <=> H2 + CH2O */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[58];
    q_f[58] = phi_f * k_f;
    phi_r = sc[0]*sc[17];
    Kc = Kc_save[58];
    k_r = k_f / Kc;
    q_r[58] = phi_r * k_r;

    /*reaction 60: H + CH2OH <=> OH + CH3 */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[59];
    q_f[59] = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = Kc_save[59];
    k_r = k_f / Kc;
    q_r[59] = phi_r * k_r;

    /*reaction 61: H + CH2OH <=> CH2(S) + H2O */
    phi_f = sc[1]*sc[18];
    k_f = k_f_save[60];
    q_f[60] = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[60];
    k_r = k_f / Kc;
    q_r[60] = phi_r * k_r;

    /*reaction 62: H + CH3O (+M) <=> CH3OH (+M) */
    phi_f = sc[1]*sc[19];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[61];
    redP = 1e-12 * alpha / k_f * 4.66e+41*exp(-7.44*tc[0]-7085.2906545304531392*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.3*exp(T/-100))+ (0.7*exp(T/-90000))+ (exp(-10000/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[61] = phi_f * k_f;
    phi_r = sc[20];
    Kc = Kc_save[61];
    k_r = k_f / Kc;
    q_r[61] = phi_r * k_r;

    /*reaction 63: H + CH3O <=> H + CH2OH */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[62];
    q_f[62] = phi_f * k_f;
    phi_r = sc[1]*sc[18];
    Kc = Kc_save[62];
    k_r = k_f / Kc;
    q_r[62] = phi_r * k_r;

    /*reaction 64: H + CH3O <=> H2 + CH2O */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[63];
    q_f[63] = phi_f * k_f;
    phi_r = sc[0]*sc[17];
    Kc = Kc_save[63];
    k_r = k_f / Kc;
    q_r[63] = phi_r * k_r;

    /*reaction 65: H + CH3O <=> OH + CH3 */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[64];
    q_f[64] = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = Kc_save[64];
    k_r = k_f / Kc;
    q_r[64] = phi_r * k_r;

    /*reaction 66: H + CH3O <=> CH2(S) + H2O */
    phi_f = sc[1]*sc[19];
    k_f = k_f_save[65];
    q_f[65] = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[65];
    k_r = k_f / Kc;
    q_r[65] = phi_r * k_r;

    /*reaction 67: H + CH3OH <=> CH2OH + H2 */
    phi_f = sc[1]*sc[20];
    k_f = k_f_save[66];
    q_f[66] = phi_f * k_f;
    phi_r = sc[18]*sc[0];
    Kc = Kc_save[66];
    k_r = k_f / Kc;
    q_r[66] = phi_r * k_r;

    /*reaction 68: H + CH3OH <=> CH3O + H2 */
    phi_f = sc[1]*sc[20];
    k_f = k_f_save[67];
    q_f[67] = phi_f * k_f;
    phi_r = sc[19]*sc[0];
    Kc = Kc_save[67];
    k_r = k_f / Kc;
    q_r[67] = phi_r * k_r;

    /*reaction 69: H + C2H (+M) <=> C2H2 (+M) */
    phi_f = sc[1]*sc[21];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[68];
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
    q_f[68] = phi_f * k_f;
    phi_r = sc[22];
    Kc = Kc_save[68];
    k_r = k_f / Kc;
    q_r[68] = phi_r * k_r;

    /*reaction 70: H + C2H2 (+M) <=> C2H3 (+M) */
    phi_f = sc[1]*sc[22];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[69];
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
    q_f[69] = phi_f * k_f;
    phi_r = sc[23];
    Kc = Kc_save[69];
    k_r = k_f / Kc;
    q_r[69] = phi_r * k_r;

    /*reaction 71: H + C2H3 (+M) <=> C2H4 (+M) */
    phi_f = sc[1]*sc[23];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[70];
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
    q_f[70] = phi_f * k_f;
    phi_r = sc[24];
    Kc = Kc_save[70];
    k_r = k_f / Kc;
    q_r[70] = phi_r * k_r;

    /*reaction 72: H + C2H3 <=> H2 + C2H2 */
    phi_f = sc[1]*sc[23];
    k_f = k_f_save[71];
    q_f[71] = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = Kc_save[71];
    k_r = k_f / Kc;
    q_r[71] = phi_r * k_r;

    /*reaction 73: H + C2H4 (+M) <=> C2H5 (+M) */
    phi_f = sc[1]*sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[72];
    redP = 1e-12 * alpha / k_f * 6e+41*exp(-7.62*tc[0]-3507.4201606588962932*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.0247*exp(T/-210))+ (0.9753*exp(T/-984))+ (exp(-4374/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[72] = phi_f * k_f;
    phi_r = sc[25];
    Kc = Kc_save[72];
    k_r = k_f / Kc;
    q_r[72] = phi_r * k_r;

    /*reaction 74: H + C2H4 <=> C2H3 + H2 */
    phi_f = sc[1]*sc[24];
    k_f = k_f_save[73];
    q_f[73] = phi_f * k_f;
    phi_r = sc[23]*sc[0];
    Kc = Kc_save[73];
    k_r = k_f / Kc;
    q_r[73] = phi_r * k_r;

    /*reaction 75: H + C2H5 (+M) <=> C2H6 (+M) */
    phi_f = sc[1]*sc[25];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[74];
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
    q_f[74] = phi_f * k_f;
    phi_r = sc[26];
    Kc = Kc_save[74];
    k_r = k_f / Kc;
    q_r[74] = phi_r * k_r;

    /*reaction 76: H + C2H5 <=> H2 + C2H4 */
    phi_f = sc[1]*sc[25];
    k_f = k_f_save[75];
    q_f[75] = phi_f * k_f;
    phi_r = sc[0]*sc[24];
    Kc = Kc_save[75];
    k_r = k_f / Kc;
    q_r[75] = phi_r * k_r;

    /*reaction 77: H + C2H6 <=> C2H5 + H2 */
    phi_f = sc[1]*sc[26];
    k_f = k_f_save[76];
    q_f[76] = phi_f * k_f;
    phi_r = sc[25]*sc[0];
    Kc = Kc_save[76];
    k_r = k_f / Kc;
    q_r[76] = phi_r * k_r;

    /*reaction 78: H + HCCO <=> CH2(S) + CO */
    phi_f = sc[1]*sc[27];
    k_f = k_f_save[77];
    q_f[77] = phi_f * k_f;
    phi_r = sc[11]*sc[14];
    Kc = Kc_save[77];
    k_r = k_f / Kc;
    q_r[77] = phi_r * k_r;

    /*reaction 79: H + CH2CO <=> HCCO + H2 */
    phi_f = sc[1]*sc[28];
    k_f = k_f_save[78];
    q_f[78] = phi_f * k_f;
    phi_r = sc[27]*sc[0];
    Kc = Kc_save[78];
    k_r = k_f / Kc;
    q_r[78] = phi_r * k_r;

    /*reaction 80: H + CH2CO <=> CH3 + CO */
    phi_f = sc[1]*sc[28];
    k_f = k_f_save[79];
    q_f[79] = phi_f * k_f;
    phi_r = sc[12]*sc[14];
    Kc = Kc_save[79];
    k_r = k_f / Kc;
    q_r[79] = phi_r * k_r;

    /*reaction 81: H + HCCOH <=> H + CH2CO */
    phi_f = sc[1]*sc[29];
    k_f = k_f_save[80];
    q_f[80] = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[80];
    k_r = k_f / Kc;
    q_r[80] = phi_r * k_r;

    /*reaction 82: H2 + CO (+M) <=> CH2O (+M) */
    phi_f = sc[0]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[81];
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
    q_f[81] = phi_f * k_f;
    phi_r = sc[17];
    Kc = Kc_save[81];
    k_r = k_f / Kc;
    q_r[81] = phi_r * k_r;

    /*reaction 83: OH + H2 <=> H + H2O */
    phi_f = sc[4]*sc[0];
    k_f = k_f_save[82];
    q_f[82] = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    Kc = Kc_save[82];
    k_r = k_f / Kc;
    q_r[82] = phi_r * k_r;

    /*reaction 84: 2 OH (+M) <=> H2O2 (+M) */
    phi_f = sc[4]*sc[4];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[83];
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
    q_f[83] = phi_f * k_f;
    phi_r = sc[7];
    Kc = Kc_save[83];
    k_r = k_f / Kc;
    q_r[83] = phi_r * k_r;

    /*reaction 85: 2 OH <=> O + H2O */
    phi_f = sc[4]*sc[4];
    k_f = k_f_save[84];
    q_f[84] = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = Kc_save[84];
    k_r = k_f / Kc;
    q_r[84] = phi_r * k_r;

    /*reaction 86: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = k_f_save[85];
    q_f[85] = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = Kc_save[85];
    k_r = k_f / Kc;
    q_r[85] = phi_r * k_r;

    /*reaction 87: OH + H2O2 <=> HO2 + H2O */
    phi_f = sc[4]*sc[7];
    k_f = k_f_save[86];
    q_f[86] = phi_f * k_f;
    phi_r = sc[6]*sc[5];
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

    /*reaction 89: OH + C <=> H + CO */
    phi_f = sc[4]*sc[8];
    k_f = k_f_save[88];
    q_f[88] = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[88];
    k_r = k_f / Kc;
    q_r[88] = phi_r * k_r;

    /*reaction 90: OH + CH <=> H + HCO */
    phi_f = sc[4]*sc[9];
    k_f = k_f_save[89];
    q_f[89] = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_save[89];
    k_r = k_f / Kc;
    q_r[89] = phi_r * k_r;

    /*reaction 91: OH + CH2 <=> H + CH2O */
    phi_f = sc[4]*sc[10];
    k_f = k_f_save[90];
    q_f[90] = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[90];
    k_r = k_f / Kc;
    q_r[90] = phi_r * k_r;

    /*reaction 92: OH + CH2 <=> CH + H2O */
    phi_f = sc[4]*sc[10];
    k_f = k_f_save[91];
    q_f[91] = phi_f * k_f;
    phi_r = sc[9]*sc[5];
    Kc = Kc_save[91];
    k_r = k_f / Kc;
    q_r[91] = phi_r * k_r;

    /*reaction 93: OH + CH2(S) <=> H + CH2O */
    phi_f = sc[4]*sc[11];
    k_f = k_f_save[92];
    q_f[92] = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[92];
    k_r = k_f / Kc;
    q_r[92] = phi_r * k_r;

    /*reaction 94: OH + CH3 (+M) <=> CH3OH (+M) */
    phi_f = sc[4]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[93];
    redP = 1e-12 * alpha / k_f * 4e+36*exp(-5.92*tc[0]-1580.1003306268198685*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.588*exp(T/-195))+ (0.412*exp(T/-5900))+ (exp(-6394/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[93] = phi_f * k_f;
    phi_r = sc[20];
    Kc = Kc_save[93];
    k_r = k_f / Kc;
    q_r[93] = phi_r * k_r;

    /*reaction 95: OH + CH3 <=> CH2 + H2O */
    phi_f = sc[4]*sc[12];
    k_f = k_f_save[94];
    q_f[94] = phi_f * k_f;
    phi_r = sc[10]*sc[5];
    Kc = Kc_save[94];
    k_r = k_f / Kc;
    q_r[94] = phi_r * k_r;

    /*reaction 96: OH + CH3 <=> CH2(S) + H2O */
    phi_f = sc[4]*sc[12];
    k_f = k_f_save[95];
    q_f[95] = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_save[95];
    k_r = k_f / Kc;
    q_r[95] = phi_r * k_r;

    /*reaction 97: OH + CH4 <=> CH3 + H2O */
    phi_f = sc[4]*sc[13];
    k_f = k_f_save[96];
    q_f[96] = phi_f * k_f;
    phi_r = sc[12]*sc[5];
    Kc = Kc_save[96];
    k_r = k_f / Kc;
    q_r[96] = phi_r * k_r;

    /*reaction 98: OH + CO <=> H + CO2 */
    phi_f = sc[4]*sc[14];
    k_f = k_f_save[97];
    q_f[97] = phi_f * k_f;
    phi_r = sc[1]*sc[15];
    Kc = Kc_save[97];
    k_r = k_f / Kc;
    q_r[97] = phi_r * k_r;

    /*reaction 99: OH + HCO <=> H2O + CO */
    phi_f = sc[4]*sc[16];
    k_f = k_f_save[98];
    q_f[98] = phi_f * k_f;
    phi_r = sc[5]*sc[14];
    Kc = Kc_save[98];
    k_r = k_f / Kc;
    q_r[98] = phi_r * k_r;

    /*reaction 100: OH + CH2O <=> HCO + H2O */
    phi_f = sc[4]*sc[17];
    k_f = k_f_save[99];
    q_f[99] = phi_f * k_f;
    phi_r = sc[16]*sc[5];
    Kc = Kc_save[99];
    k_r = k_f / Kc;
    q_r[99] = phi_r * k_r;

    /*reaction 101: OH + CH2OH <=> H2O + CH2O */
    phi_f = sc[4]*sc[18];
    k_f = k_f_save[100];
    q_f[100] = phi_f * k_f;
    phi_r = sc[5]*sc[17];
    Kc = Kc_save[100];
    k_r = k_f / Kc;
    q_r[100] = phi_r * k_r;

    /*reaction 102: OH + CH3O <=> H2O + CH2O */
    phi_f = sc[4]*sc[19];
    k_f = k_f_save[101];
    q_f[101] = phi_f * k_f;
    phi_r = sc[5]*sc[17];
    Kc = Kc_save[101];
    k_r = k_f / Kc;
    q_r[101] = phi_r * k_r;

    /*reaction 103: OH + CH3OH <=> CH2OH + H2O */
    phi_f = sc[4]*sc[20];
    k_f = k_f_save[102];
    q_f[102] = phi_f * k_f;
    phi_r = sc[18]*sc[5];
    Kc = Kc_save[102];
    k_r = k_f / Kc;
    q_r[102] = phi_r * k_r;

    /*reaction 104: OH + CH3OH <=> CH3O + H2O */
    phi_f = sc[4]*sc[20];
    k_f = k_f_save[103];
    q_f[103] = phi_f * k_f;
    phi_r = sc[19]*sc[5];
    Kc = Kc_save[103];
    k_r = k_f / Kc;
    q_r[103] = phi_r * k_r;

    /*reaction 105: OH + C2H <=> H + HCCO */
    phi_f = sc[4]*sc[21];
    k_f = k_f_save[104];
    q_f[104] = phi_f * k_f;
    phi_r = sc[1]*sc[27];
    Kc = Kc_save[104];
    k_r = k_f / Kc;
    q_r[104] = phi_r * k_r;

    /*reaction 106: OH + C2H2 <=> H + CH2CO */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[105];
    q_f[105] = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[105];
    k_r = k_f / Kc;
    q_r[105] = phi_r * k_r;

    /*reaction 107: OH + C2H2 <=> H + HCCOH */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[106];
    q_f[106] = phi_f * k_f;
    phi_r = sc[1]*sc[29];
    Kc = Kc_save[106];
    k_r = k_f / Kc;
    q_r[106] = phi_r * k_r;

    /*reaction 108: OH + C2H2 <=> C2H + H2O */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[107];
    q_f[107] = phi_f * k_f;
    phi_r = sc[21]*sc[5];
    Kc = Kc_save[107];
    k_r = k_f / Kc;
    q_r[107] = phi_r * k_r;

    /*reaction 109: OH + C2H2 <=> CH3 + CO */
    phi_f = sc[4]*sc[22];
    k_f = k_f_save[108];
    q_f[108] = phi_f * k_f;
    phi_r = sc[12]*sc[14];
    Kc = Kc_save[108];
    k_r = k_f / Kc;
    q_r[108] = phi_r * k_r;

    /*reaction 110: OH + C2H3 <=> H2O + C2H2 */
    phi_f = sc[4]*sc[23];
    k_f = k_f_save[109];
    q_f[109] = phi_f * k_f;
    phi_r = sc[5]*sc[22];
    Kc = Kc_save[109];
    k_r = k_f / Kc;
    q_r[109] = phi_r * k_r;

    /*reaction 111: OH + C2H4 <=> C2H3 + H2O */
    phi_f = sc[4]*sc[24];
    k_f = k_f_save[110];
    q_f[110] = phi_f * k_f;
    phi_r = sc[23]*sc[5];
    Kc = Kc_save[110];
    k_r = k_f / Kc;
    q_r[110] = phi_r * k_r;

    /*reaction 112: OH + C2H6 <=> C2H5 + H2O */
    phi_f = sc[4]*sc[26];
    k_f = k_f_save[111];
    q_f[111] = phi_f * k_f;
    phi_r = sc[25]*sc[5];
    Kc = Kc_save[111];
    k_r = k_f / Kc;
    q_r[111] = phi_r * k_r;

    /*reaction 113: OH + CH2CO <=> HCCO + H2O */
    phi_f = sc[4]*sc[28];
    k_f = k_f_save[112];
    q_f[112] = phi_f * k_f;
    phi_r = sc[27]*sc[5];
    Kc = Kc_save[112];
    k_r = k_f / Kc;
    q_r[112] = phi_r * k_r;

    /*reaction 114: 2 HO2 <=> O2 + H2O2 */
    phi_f = sc[6]*sc[6];
    k_f = k_f_save[113];
    q_f[113] = phi_f * k_f;
    phi_r = sc[3]*sc[7];
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

    /*reaction 116: HO2 + CH2 <=> OH + CH2O */
    phi_f = sc[6]*sc[10];
    k_f = k_f_save[115];
    q_f[115] = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[115];
    k_r = k_f / Kc;
    q_r[115] = phi_r * k_r;

    /*reaction 117: HO2 + CH3 <=> O2 + CH4 */
    phi_f = sc[6]*sc[12];
    k_f = k_f_save[116];
    q_f[116] = phi_f * k_f;
    phi_r = sc[3]*sc[13];
    Kc = Kc_save[116];
    k_r = k_f / Kc;
    q_r[116] = phi_r * k_r;

    /*reaction 118: HO2 + CH3 <=> OH + CH3O */
    phi_f = sc[6]*sc[12];
    k_f = k_f_save[117];
    q_f[117] = phi_f * k_f;
    phi_r = sc[4]*sc[19];
    Kc = Kc_save[117];
    k_r = k_f / Kc;
    q_r[117] = phi_r * k_r;

    /*reaction 119: HO2 + CO <=> OH + CO2 */
    phi_f = sc[6]*sc[14];
    k_f = k_f_save[118];
    q_f[118] = phi_f * k_f;
    phi_r = sc[4]*sc[15];
    Kc = Kc_save[118];
    k_r = k_f / Kc;
    q_r[118] = phi_r * k_r;

    /*reaction 120: HO2 + CH2O <=> HCO + H2O2 */
    phi_f = sc[6]*sc[17];
    k_f = k_f_save[119];
    q_f[119] = phi_f * k_f;
    phi_r = sc[16]*sc[7];
    Kc = Kc_save[119];
    k_r = k_f / Kc;
    q_r[119] = phi_r * k_r;

    /*reaction 121: C + O2 <=> O + CO */
    phi_f = sc[8]*sc[3];
    k_f = k_f_save[120];
    q_f[120] = phi_f * k_f;
    phi_r = sc[2]*sc[14];
    Kc = Kc_save[120];
    k_r = k_f / Kc;
    q_r[120] = phi_r * k_r;

    /*reaction 122: C + CH2 <=> H + C2H */
    phi_f = sc[8]*sc[10];
    k_f = k_f_save[121];
    q_f[121] = phi_f * k_f;
    phi_r = sc[1]*sc[21];
    Kc = Kc_save[121];
    k_r = k_f / Kc;
    q_r[121] = phi_r * k_r;

    /*reaction 123: C + CH3 <=> H + C2H2 */
    phi_f = sc[8]*sc[12];
    k_f = k_f_save[122];
    q_f[122] = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = Kc_save[122];
    k_r = k_f / Kc;
    q_r[122] = phi_r * k_r;

    /*reaction 124: CH + O2 <=> O + HCO */
    phi_f = sc[9]*sc[3];
    k_f = k_f_save[123];
    q_f[123] = phi_f * k_f;
    phi_r = sc[2]*sc[16];
    Kc = Kc_save[123];
    k_r = k_f / Kc;
    q_r[123] = phi_r * k_r;

    /*reaction 125: CH + H2 <=> H + CH2 */
    phi_f = sc[9]*sc[0];
    k_f = k_f_save[124];
    q_f[124] = phi_f * k_f;
    phi_r = sc[1]*sc[10];
    Kc = Kc_save[124];
    k_r = k_f / Kc;
    q_r[124] = phi_r * k_r;

    /*reaction 126: CH + H2O <=> H + CH2O */
    phi_f = sc[9]*sc[5];
    k_f = k_f_save[125];
    q_f[125] = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_save[125];
    k_r = k_f / Kc;
    q_r[125] = phi_r * k_r;

    /*reaction 127: CH + CH2 <=> H + C2H2 */
    phi_f = sc[9]*sc[10];
    k_f = k_f_save[126];
    q_f[126] = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = Kc_save[126];
    k_r = k_f / Kc;
    q_r[126] = phi_r * k_r;

    /*reaction 128: CH + CH3 <=> H + C2H3 */
    phi_f = sc[9]*sc[12];
    k_f = k_f_save[127];
    q_f[127] = phi_f * k_f;
    phi_r = sc[1]*sc[23];
    Kc = Kc_save[127];
    k_r = k_f / Kc;
    q_r[127] = phi_r * k_r;

    /*reaction 129: CH + CH4 <=> H + C2H4 */
    phi_f = sc[9]*sc[13];
    k_f = k_f_save[128];
    q_f[128] = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = Kc_save[128];
    k_r = k_f / Kc;
    q_r[128] = phi_r * k_r;

    /*reaction 130: CH + CO (+M) <=> HCCO (+M) */
    phi_f = sc[9]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[129];
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
    q_f[129] = phi_f * k_f;
    phi_r = sc[27];
    Kc = Kc_save[129];
    k_r = k_f / Kc;
    q_r[129] = phi_r * k_r;

    /*reaction 131: CH + CO2 <=> HCO + CO */
    phi_f = sc[9]*sc[15];
    k_f = k_f_save[130];
    q_f[130] = phi_f * k_f;
    phi_r = sc[16]*sc[14];
    Kc = Kc_save[130];
    k_r = k_f / Kc;
    q_r[130] = phi_r * k_r;

    /*reaction 132: CH + CH2O <=> H + CH2CO */
    phi_f = sc[9]*sc[17];
    k_f = k_f_save[131];
    q_f[131] = phi_f * k_f;
    phi_r = sc[1]*sc[28];
    Kc = Kc_save[131];
    k_r = k_f / Kc;
    q_r[131] = phi_r * k_r;

    /*reaction 133: CH + HCCO <=> CO + C2H2 */
    phi_f = sc[9]*sc[27];
    k_f = k_f_save[132];
    q_f[132] = phi_f * k_f;
    phi_r = sc[14]*sc[22];
    Kc = Kc_save[132];
    k_r = k_f / Kc;
    q_r[132] = phi_r * k_r;

    /*reaction 134: CH2 + O2 => OH + H + CO */
    phi_f = sc[10]*sc[3];
    k_f = k_f_save[133];
    q_f[133] = phi_f * k_f;
    q_r[133] = 0.0;

    /*reaction 135: CH2 + H2 <=> H + CH3 */
    phi_f = sc[10]*sc[0];
    k_f = k_f_save[134];
    q_f[134] = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = Kc_save[134];
    k_r = k_f / Kc;
    q_r[134] = phi_r * k_r;

    /*reaction 136: 2 CH2 <=> H2 + C2H2 */
    phi_f = sc[10]*sc[10];
    k_f = k_f_save[135];
    q_f[135] = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = Kc_save[135];
    k_r = k_f / Kc;
    q_r[135] = phi_r * k_r;

    /*reaction 137: CH2 + CH3 <=> H + C2H4 */
    phi_f = sc[10]*sc[12];
    k_f = k_f_save[136];
    q_f[136] = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = Kc_save[136];
    k_r = k_f / Kc;
    q_r[136] = phi_r * k_r;

    /*reaction 138: CH2 + CH4 <=> 2 CH3 */
    phi_f = sc[10]*sc[13];
    k_f = k_f_save[137];
    q_f[137] = phi_f * k_f;
    phi_r = sc[12]*sc[12];
    Kc = Kc_save[137];
    k_r = k_f / Kc;
    q_r[137] = phi_r * k_r;

    /*reaction 139: CH2 + CO (+M) <=> CH2CO (+M) */
    phi_f = sc[10]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[138];
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
    q_f[138] = phi_f * k_f;
    phi_r = sc[28];
    Kc = Kc_save[138];
    k_r = k_f / Kc;
    q_r[138] = phi_r * k_r;

    /*reaction 140: CH2 + HCCO <=> C2H3 + CO */
    phi_f = sc[10]*sc[27];
    k_f = k_f_save[139];
    q_f[139] = phi_f * k_f;
    phi_r = sc[23]*sc[14];
    Kc = Kc_save[139];
    k_r = k_f / Kc;
    q_r[139] = phi_r * k_r;

    /*reaction 141: CH2(S) + N2 <=> CH2 + N2 */
    phi_f = sc[11]*sc[30];
    k_f = k_f_save[140];
    q_f[140] = phi_f * k_f;
    phi_r = sc[10]*sc[30];
    Kc = Kc_save[140];
    k_r = k_f / Kc;
    q_r[140] = phi_r * k_r;

    /*reaction 142: CH2(S) + O2 <=> H + OH + CO */
    phi_f = sc[11]*sc[3];
    k_f = k_f_save[141];
    q_f[141] = phi_f * k_f;
    phi_r = sc[1]*sc[4]*sc[14];
    Kc = Kc_save[141];
    k_r = k_f / Kc;
    q_r[141] = phi_r * k_r;

    /*reaction 143: CH2(S) + O2 <=> CO + H2O */
    phi_f = sc[11]*sc[3];
    k_f = k_f_save[142];
    q_f[142] = phi_f * k_f;
    phi_r = sc[14]*sc[5];
    Kc = Kc_save[142];
    k_r = k_f / Kc;
    q_r[142] = phi_r * k_r;

    /*reaction 144: CH2(S) + H2 <=> CH3 + H */
    phi_f = sc[11]*sc[0];
    k_f = k_f_save[143];
    q_f[143] = phi_f * k_f;
    phi_r = sc[12]*sc[1];
    Kc = Kc_save[143];
    k_r = k_f / Kc;
    q_r[143] = phi_r * k_r;

    /*reaction 145: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    phi_f = sc[11]*sc[5];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[144];
    redP = 1e-12 * alpha / k_f * 1.88e+38*exp(-6.36*tc[0]-2536.2119956557871774*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.3973*exp(T/-208))+ (0.6027*exp(T/-3922))+ (exp(-10180/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[144] = phi_f * k_f;
    phi_r = sc[20];
    Kc = Kc_save[144];
    k_r = k_f / Kc;
    q_r[144] = phi_r * k_r;

    /*reaction 146: CH2(S) + H2O <=> CH2 + H2O */
    phi_f = sc[11]*sc[5];
    k_f = k_f_save[145];
    q_f[145] = phi_f * k_f;
    phi_r = sc[10]*sc[5];
    Kc = Kc_save[145];
    k_r = k_f / Kc;
    q_r[145] = phi_r * k_r;

    /*reaction 147: CH2(S) + CH3 <=> H + C2H4 */
    phi_f = sc[11]*sc[12];
    k_f = k_f_save[146];
    q_f[146] = phi_f * k_f;
    phi_r = sc[1]*sc[24];
    Kc = Kc_save[146];
    k_r = k_f / Kc;
    q_r[146] = phi_r * k_r;

    /*reaction 148: CH2(S) + CH4 <=> 2 CH3 */
    phi_f = sc[11]*sc[13];
    k_f = k_f_save[147];
    q_f[147] = phi_f * k_f;
    phi_r = sc[12]*sc[12];
    Kc = Kc_save[147];
    k_r = k_f / Kc;
    q_r[147] = phi_r * k_r;

    /*reaction 149: CH2(S) + CO <=> CH2 + CO */
    phi_f = sc[11]*sc[14];
    k_f = k_f_save[148];
    q_f[148] = phi_f * k_f;
    phi_r = sc[10]*sc[14];
    Kc = Kc_save[148];
    k_r = k_f / Kc;
    q_r[148] = phi_r * k_r;

    /*reaction 150: CH2(S) + CO2 <=> CH2 + CO2 */
    phi_f = sc[11]*sc[15];
    k_f = k_f_save[149];
    q_f[149] = phi_f * k_f;
    phi_r = sc[10]*sc[15];
    Kc = Kc_save[149];
    k_r = k_f / Kc;
    q_r[149] = phi_r * k_r;

    /*reaction 151: CH2(S) + CO2 <=> CO + CH2O */
    phi_f = sc[11]*sc[15];
    k_f = k_f_save[150];
    q_f[150] = phi_f * k_f;
    phi_r = sc[14]*sc[17];
    Kc = Kc_save[150];
    k_r = k_f / Kc;
    q_r[150] = phi_r * k_r;

    /*reaction 152: CH2(S) + C2H6 <=> CH3 + C2H5 */
    phi_f = sc[11]*sc[26];
    k_f = k_f_save[151];
    q_f[151] = phi_f * k_f;
    phi_r = sc[12]*sc[25];
    Kc = Kc_save[151];
    k_r = k_f / Kc;
    q_r[151] = phi_r * k_r;

    /*reaction 153: CH3 + O2 <=> O + CH3O */
    phi_f = sc[12]*sc[3];
    k_f = k_f_save[152];
    q_f[152] = phi_f * k_f;
    phi_r = sc[2]*sc[19];
    Kc = Kc_save[152];
    k_r = k_f / Kc;
    q_r[152] = phi_r * k_r;

    /*reaction 154: CH3 + O2 <=> OH + CH2O */
    phi_f = sc[12]*sc[3];
    k_f = k_f_save[153];
    q_f[153] = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_save[153];
    k_r = k_f / Kc;
    q_r[153] = phi_r * k_r;

    /*reaction 155: CH3 + H2O2 <=> HO2 + CH4 */
    phi_f = sc[12]*sc[7];
    k_f = k_f_save[154];
    q_f[154] = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = Kc_save[154];
    k_r = k_f / Kc;
    q_r[154] = phi_r * k_r;

    /*reaction 156: 2 CH3 (+M) <=> C2H6 (+M) */
    phi_f = sc[12]*sc[12];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[155];
    redP = 1e-12 * alpha / k_f * 3.4e+41*exp(-7.03*tc[0]-1389.8844309526357392*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.381*exp(T/-73.2))+ (0.619*exp(T/-1180))+ (exp(-9999/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[155] = phi_f * k_f;
    phi_r = sc[26];
    Kc = Kc_save[155];
    k_r = k_f / Kc;
    q_r[155] = phi_r * k_r;

    /*reaction 157: 2 CH3 <=> H + C2H5 */
    phi_f = sc[12]*sc[12];
    k_f = k_f_save[156];
    q_f[156] = phi_f * k_f;
    phi_r = sc[1]*sc[25];
    Kc = Kc_save[156];
    k_r = k_f / Kc;
    q_r[156] = phi_r * k_r;

    /*reaction 158: CH3 + HCO <=> CH4 + CO */
    phi_f = sc[12]*sc[16];
    k_f = k_f_save[157];
    q_f[157] = phi_f * k_f;
    phi_r = sc[13]*sc[14];
    Kc = Kc_save[157];
    k_r = k_f / Kc;
    q_r[157] = phi_r * k_r;

    /*reaction 159: CH3 + CH2O <=> HCO + CH4 */
    phi_f = sc[12]*sc[17];
    k_f = k_f_save[158];
    q_f[158] = phi_f * k_f;
    phi_r = sc[16]*sc[13];
    Kc = Kc_save[158];
    k_r = k_f / Kc;
    q_r[158] = phi_r * k_r;

    /*reaction 160: CH3 + CH3OH <=> CH2OH + CH4 */
    phi_f = sc[12]*sc[20];
    k_f = k_f_save[159];
    q_f[159] = phi_f * k_f;
    phi_r = sc[18]*sc[13];
    Kc = Kc_save[159];
    k_r = k_f / Kc;
    q_r[159] = phi_r * k_r;

    /*reaction 161: CH3 + CH3OH <=> CH3O + CH4 */
    phi_f = sc[12]*sc[20];
    k_f = k_f_save[160];
    q_f[160] = phi_f * k_f;
    phi_r = sc[19]*sc[13];
    Kc = Kc_save[160];
    k_r = k_f / Kc;
    q_r[160] = phi_r * k_r;

    /*reaction 162: CH3 + C2H4 <=> C2H3 + CH4 */
    phi_f = sc[12]*sc[24];
    k_f = k_f_save[161];
    q_f[161] = phi_f * k_f;
    phi_r = sc[23]*sc[13];
    Kc = Kc_save[161];
    k_r = k_f / Kc;
    q_r[161] = phi_r * k_r;

    /*reaction 163: CH3 + C2H6 <=> C2H5 + CH4 */
    phi_f = sc[12]*sc[26];
    k_f = k_f_save[162];
    q_f[162] = phi_f * k_f;
    phi_r = sc[25]*sc[13];
    Kc = Kc_save[162];
    k_r = k_f / Kc;
    q_r[162] = phi_r * k_r;

    /*reaction 164: HCO + H2O <=> H + CO + H2O */
    phi_f = sc[16]*sc[5];
    k_f = k_f_save[163];
    q_f[163] = phi_f * k_f;
    phi_r = sc[1]*sc[14]*sc[5];
    Kc = Kc_save[163];
    k_r = k_f / Kc;
    q_r[163] = phi_r * k_r;

    /*reaction 165: HCO + M <=> H + CO + M */
    phi_f = sc[16];
    alpha = mixture + sc[0] + -1*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = alpha * k_f_save[164];
    q_f[164] = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_save[164];
    k_r = k_f / Kc;
    q_r[164] = phi_r * k_r;

    /*reaction 166: HCO + O2 <=> HO2 + CO */
    phi_f = sc[16]*sc[3];
    k_f = k_f_save[165];
    q_f[165] = phi_f * k_f;
    phi_r = sc[6]*sc[14];
    Kc = Kc_save[165];
    k_r = k_f / Kc;
    q_r[165] = phi_r * k_r;

    /*reaction 167: CH2OH + O2 <=> HO2 + CH2O */
    phi_f = sc[18]*sc[3];
    k_f = k_f_save[166];
    q_f[166] = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = Kc_save[166];
    k_r = k_f / Kc;
    q_r[166] = phi_r * k_r;

    /*reaction 168: CH3O + O2 <=> HO2 + CH2O */
    phi_f = sc[19]*sc[3];
    k_f = k_f_save[167];
    q_f[167] = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = Kc_save[167];
    k_r = k_f / Kc;
    q_r[167] = phi_r * k_r;

    /*reaction 169: C2H + O2 <=> HCO + CO */
    phi_f = sc[21]*sc[3];
    k_f = k_f_save[168];
    q_f[168] = phi_f * k_f;
    phi_r = sc[16]*sc[14];
    Kc = Kc_save[168];
    k_r = k_f / Kc;
    q_r[168] = phi_r * k_r;

    /*reaction 170: C2H + H2 <=> H + C2H2 */
    phi_f = sc[21]*sc[0];
    k_f = k_f_save[169];
    q_f[169] = phi_f * k_f;
    phi_r = sc[1]*sc[22];
    Kc = Kc_save[169];
    k_r = k_f / Kc;
    q_r[169] = phi_r * k_r;

    /*reaction 171: C2H3 + O2 <=> HCO + CH2O */
    phi_f = sc[23]*sc[3];
    k_f = k_f_save[170];
    q_f[170] = phi_f * k_f;
    phi_r = sc[16]*sc[17];
    Kc = Kc_save[170];
    k_r = k_f / Kc;
    q_r[170] = phi_r * k_r;

    /*reaction 172: C2H4 (+M) <=> H2 + C2H2 (+M) */
    phi_f = sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[171];
    redP = 1e-6 * alpha / k_f * 1.58e+51*exp(-9.3*tc[0]-49214.589915701588325*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.2655*exp(T/-180))+ (0.7345*exp(T/-1035))+ (exp(-5417/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[171] = phi_f * k_f;
    phi_r = sc[0]*sc[22];
    Kc = Kc_save[171];
    k_r = k_f / Kc;
    q_r[171] = phi_r * k_r;

    /*reaction 173: C2H5 + O2 <=> HO2 + C2H4 */
    phi_f = sc[25]*sc[3];
    k_f = k_f_save[172];
    q_f[172] = phi_f * k_f;
    phi_r = sc[6]*sc[24];
    Kc = Kc_save[172];
    k_r = k_f / Kc;
    q_r[172] = phi_r * k_r;

    /*reaction 174: HCCO + O2 <=> OH + 2 CO */
    phi_f = sc[27]*sc[3];
    k_f = k_f_save[173];
    q_f[173] = phi_f * k_f;
    phi_r = sc[4]*sc[14]*sc[14];
    Kc = Kc_save[173];
    k_r = k_f / Kc;
    q_r[173] = phi_r * k_r;

    /*reaction 175: 2 HCCO <=> 2 CO + C2H2 */
    phi_f = sc[27]*sc[27];
    k_f = k_f_save[174];
    q_f[174] = phi_f * k_f;
    phi_r = sc[14]*sc[14]*sc[22];
    Kc = Kc_save[174];
    k_r = k_f / Kc;
    q_r[174] = phi_r * k_r;

    /*reaction 176: O + CH3 => H + H2 + CO */
    phi_f = sc[2]*sc[12];
    k_f = k_f_save[175];
    q_f[175] = phi_f * k_f;
    q_r[175] = 0.0;

    /*reaction 177: O + C2H4 <=> H + CH2CHO */
    phi_f = sc[2]*sc[24];
    k_f = k_f_save[176];
    q_f[176] = phi_f * k_f;
    phi_r = sc[1]*sc[33];
    Kc = Kc_save[176];
    k_r = k_f / Kc;
    q_r[176] = phi_r * k_r;

    /*reaction 178: O + C2H5 <=> H + CH3CHO */
    phi_f = sc[2]*sc[25];
    k_f = k_f_save[177];
    q_f[177] = phi_f * k_f;
    phi_r = sc[1]*sc[34];
    Kc = Kc_save[177];
    k_r = k_f / Kc;
    q_r[177] = phi_r * k_r;

    /*reaction 179: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = k_f_save[178];
    q_f[178] = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = Kc_save[178];
    k_r = k_f / Kc;
    q_r[178] = phi_r * k_r;

    /*reaction 180: OH + CH3 => H2 + CH2O */
    phi_f = sc[4]*sc[12];
    k_f = k_f_save[179];
    q_f[179] = phi_f * k_f;
    q_r[179] = 0.0;

    /*reaction 181: CH + H2 (+M) <=> CH3 (+M) */
    phi_f = sc[9]*sc[0];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[180];
    redP = 1e-12 * alpha / k_f * 4.82e+25*exp(-2.8*tc[0]-296.89783282478464344*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.422*exp(T/-122))+ (0.578*exp(T/-2535))+ (exp(-9365/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[180] = phi_f * k_f;
    phi_r = sc[12];
    Kc = Kc_save[180];
    k_r = k_f / Kc;
    q_r[180] = phi_r * k_r;

    /*reaction 182: CH2 + O2 => 2 H + CO2 */
    phi_f = sc[10]*sc[3];
    k_f = k_f_save[181];
    q_f[181] = phi_f * k_f;
    q_r[181] = 0.0;

    /*reaction 183: CH2 + O2 <=> O + CH2O */
    phi_f = sc[10]*sc[3];
    k_f = k_f_save[182];
    q_f[182] = phi_f * k_f;
    phi_r = sc[2]*sc[17];
    Kc = Kc_save[182];
    k_r = k_f / Kc;
    q_r[182] = phi_r * k_r;

    /*reaction 184: CH2 + CH2 => 2 H + C2H2 */
    phi_f = sc[10]*sc[10];
    k_f = k_f_save[183];
    q_f[183] = phi_f * k_f;
    q_r[183] = 0.0;

    /*reaction 185: CH2(S) + H2O => H2 + CH2O */
    phi_f = sc[11]*sc[5];
    k_f = k_f_save[184];
    q_f[184] = phi_f * k_f;
    q_r[184] = 0.0;

    /*reaction 186: C2H3 + O2 <=> O + CH2CHO */
    phi_f = sc[23]*sc[3];
    k_f = k_f_save[185];
    q_f[185] = phi_f * k_f;
    phi_r = sc[2]*sc[33];
    Kc = Kc_save[185];
    k_r = k_f / Kc;
    q_r[185] = phi_r * k_r;

    /*reaction 187: C2H3 + O2 <=> HO2 + C2H2 */
    phi_f = sc[23]*sc[3];
    k_f = k_f_save[186];
    q_f[186] = phi_f * k_f;
    phi_r = sc[6]*sc[22];
    Kc = Kc_save[186];
    k_r = k_f / Kc;
    q_r[186] = phi_r * k_r;

    /*reaction 188: O + CH3CHO <=> OH + CH2CHO */
    phi_f = sc[2]*sc[34];
    k_f = k_f_save[187];
    q_f[187] = phi_f * k_f;
    phi_r = sc[4]*sc[33];
    Kc = Kc_save[187];
    k_r = k_f / Kc;
    q_r[187] = phi_r * k_r;

    /*reaction 189: O + CH3CHO => OH + CH3 + CO */
    phi_f = sc[2]*sc[34];
    k_f = k_f_save[188];
    q_f[188] = phi_f * k_f;
    q_r[188] = 0.0;

    /*reaction 190: O2 + CH3CHO => HO2 + CH3 + CO */
    phi_f = sc[3]*sc[34];
    k_f = k_f_save[189];
    q_f[189] = phi_f * k_f;
    q_r[189] = 0.0;

    /*reaction 191: H + CH3CHO <=> CH2CHO + H2 */
    phi_f = sc[1]*sc[34];
    k_f = k_f_save[190];
    q_f[190] = phi_f * k_f;
    phi_r = sc[33]*sc[0];
    Kc = Kc_save[190];
    k_r = k_f / Kc;
    q_r[190] = phi_r * k_r;

    /*reaction 192: H + CH3CHO => CH3 + H2 + CO */
    phi_f = sc[1]*sc[34];
    k_f = k_f_save[191];
    q_f[191] = phi_f * k_f;
    q_r[191] = 0.0;

    /*reaction 193: OH + CH3CHO => CH3 + H2O + CO */
    phi_f = sc[4]*sc[34];
    k_f = k_f_save[192];
    q_f[192] = phi_f * k_f;
    q_r[192] = 0.0;

    /*reaction 194: HO2 + CH3CHO => CH3 + H2O2 + CO */
    phi_f = sc[6]*sc[34];
    k_f = k_f_save[193];
    q_f[193] = phi_f * k_f;
    q_r[193] = 0.0;

    /*reaction 195: CH3 + CH3CHO => CH3 + CH4 + CO */
    phi_f = sc[12]*sc[34];
    k_f = k_f_save[194];
    q_f[194] = phi_f * k_f;
    q_r[194] = 0.0;

    /*reaction 196: H + CH2CO (+M) <=> CH2CHO (+M) */
    phi_f = sc[1]*sc[28];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[195];
    redP = 1e-12 * alpha / k_f * 1.012e+42*exp(-7.63*tc[0]-1939.3970300113896883*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.535*exp(T/-201))+ (0.465*exp(T/-1773))+ (exp(-5333/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[195] = phi_f * k_f;
    phi_r = sc[33];
    Kc = Kc_save[195];
    k_r = k_f / Kc;
    q_r[195] = phi_r * k_r;

    /*reaction 197: O + CH2CHO => H + CH2 + CO2 */
    phi_f = sc[2]*sc[33];
    k_f = k_f_save[196];
    q_f[196] = phi_f * k_f;
    q_r[196] = 0.0;

    /*reaction 198: O2 + CH2CHO => OH + CO + CH2O */
    phi_f = sc[3]*sc[33];
    k_f = k_f_save[197];
    q_f[197] = phi_f * k_f;
    q_r[197] = 0.0;

    /*reaction 199: O2 + CH2CHO => OH + 2 HCO */
    phi_f = sc[3]*sc[33];
    k_f = k_f_save[198];
    q_f[198] = phi_f * k_f;
    q_r[198] = 0.0;

    /*reaction 200: H + CH2CHO <=> CH3 + HCO */
    phi_f = sc[1]*sc[33];
    k_f = k_f_save[199];
    q_f[199] = phi_f * k_f;
    phi_r = sc[12]*sc[16];
    Kc = Kc_save[199];
    k_r = k_f / Kc;
    q_r[199] = phi_r * k_r;

    /*reaction 201: H + CH2CHO <=> CH2CO + H2 */
    phi_f = sc[1]*sc[33];
    k_f = k_f_save[200];
    q_f[200] = phi_f * k_f;
    phi_r = sc[28]*sc[0];
    Kc = Kc_save[200];
    k_r = k_f / Kc;
    q_r[200] = phi_r * k_r;

    /*reaction 202: OH + CH2CHO <=> H2O + CH2CO */
    phi_f = sc[4]*sc[33];
    k_f = k_f_save[201];
    q_f[201] = phi_f * k_f;
    phi_r = sc[5]*sc[28];
    Kc = Kc_save[201];
    k_r = k_f / Kc;
    q_r[201] = phi_r * k_r;

    /*reaction 203: OH + CH2CHO <=> HCO + CH2OH */
    phi_f = sc[4]*sc[33];
    k_f = k_f_save[202];
    q_f[202] = phi_f * k_f;
    phi_r = sc[16]*sc[18];
    Kc = Kc_save[202];
    k_r = k_f / Kc;
    q_r[202] = phi_r * k_r;

    /*reaction 204: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    phi_f = sc[12]*sc[25];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[203];
    redP = 1e-12 * alpha / k_f * 2.71e+74*exp(-16.82*tc[0]-6574.525738738662767*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.8473*exp(T/-291))+ (0.1527*exp(T/-2742))+ (exp(-7748/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[203] = phi_f * k_f;
    phi_r = sc[32];
    Kc = Kc_save[203];
    k_r = k_f / Kc;
    q_r[203] = phi_r * k_r;

    /*reaction 205: O + C3H8 <=> OH + C3H7 */
    phi_f = sc[2]*sc[32];
    k_f = k_f_save[204];
    q_f[204] = phi_f * k_f;
    phi_r = sc[4]*sc[31];
    Kc = Kc_save[204];
    k_r = k_f / Kc;
    q_r[204] = phi_r * k_r;

    /*reaction 206: H + C3H8 <=> C3H7 + H2 */
    phi_f = sc[1]*sc[32];
    k_f = k_f_save[205];
    q_f[205] = phi_f * k_f;
    phi_r = sc[31]*sc[0];
    Kc = Kc_save[205];
    k_r = k_f / Kc;
    q_r[205] = phi_r * k_r;

    /*reaction 207: OH + C3H8 <=> C3H7 + H2O */
    phi_f = sc[4]*sc[32];
    k_f = k_f_save[206];
    q_f[206] = phi_f * k_f;
    phi_r = sc[31]*sc[5];
    Kc = Kc_save[206];
    k_r = k_f / Kc;
    q_r[206] = phi_r * k_r;

    /*reaction 208: C3H7 + H2O2 <=> HO2 + C3H8 */
    phi_f = sc[31]*sc[7];
    k_f = k_f_save[207];
    q_f[207] = phi_f * k_f;
    phi_r = sc[6]*sc[32];
    Kc = Kc_save[207];
    k_r = k_f / Kc;
    q_r[207] = phi_r * k_r;

    /*reaction 209: CH3 + C3H8 <=> C3H7 + CH4 */
    phi_f = sc[12]*sc[32];
    k_f = k_f_save[208];
    q_f[208] = phi_f * k_f;
    phi_r = sc[31]*sc[13];
    Kc = Kc_save[208];
    k_r = k_f / Kc;
    q_r[208] = phi_r * k_r;

    /*reaction 210: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    phi_f = sc[12]*sc[24];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[209];
    redP = 1e-12 * alpha / k_f * 3e+63*exp(-14.6*tc[0]-9143.4468176717564347*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.8106*exp(T/-277))+ (0.1894*exp(T/-8748))+ (exp(-7891/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[209] = phi_f * k_f;
    phi_r = sc[31];
    Kc = Kc_save[209];
    k_r = k_f / Kc;
    q_r[209] = phi_r * k_r;

    /*reaction 211: O + C3H7 <=> C2H5 + CH2O */
    phi_f = sc[2]*sc[31];
    k_f = k_f_save[210];
    q_f[210] = phi_f * k_f;
    phi_r = sc[25]*sc[17];
    Kc = Kc_save[210];
    k_r = k_f / Kc;
    q_r[210] = phi_r * k_r;

    /*reaction 212: H + C3H7 (+M) <=> C3H8 (+M) */
    phi_f = sc[1]*sc[31];
    alpha = mixture + sc[0] + 5*sc[5] + sc[13] + 0.5*sc[14] + sc[15] + 2*sc[26];
    k_f = k_f_save[211];
    redP = 1e-12 * alpha / k_f * 4.42e+61*exp(-13.545*tc[0]-5715.0316735442011122*invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10((0.685*exp(T/-369))+ (0.315*exp(T/-3285))+ (exp(-6667/T)));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[211] = phi_f * k_f;
    phi_r = sc[32];
    Kc = Kc_save[211];
    k_r = k_f / Kc;
    q_r[211] = phi_r * k_r;

    /*reaction 213: H + C3H7 <=> CH3 + C2H5 */
    phi_f = sc[1]*sc[31];
    k_f = k_f_save[212];
    q_f[212] = phi_f * k_f;
    phi_r = sc[12]*sc[25];
    Kc = Kc_save[212];
    k_r = k_f / Kc;
    q_r[212] = phi_r * k_r;

    /*reaction 214: OH + C3H7 <=> C2H5 + CH2OH */
    phi_f = sc[4]*sc[31];
    k_f = k_f_save[213];
    q_f[213] = phi_f * k_f;
    phi_r = sc[25]*sc[18];
    Kc = Kc_save[213];
    k_r = k_f / Kc;
    q_r[213] = phi_r * k_r;

    /*reaction 215: HO2 + C3H7 <=> O2 + C3H8 */
    phi_f = sc[6]*sc[31];
    k_f = k_f_save[214];
    q_f[214] = phi_f * k_f;
    phi_r = sc[3]*sc[32];
    Kc = Kc_save[214];
    k_r = k_f / Kc;
    q_r[214] = phi_r * k_r;

    /*reaction 216: HO2 + C3H7 => OH + C2H5 + CH2O */
    phi_f = sc[6]*sc[31];
    k_f = k_f_save[215];
    q_f[215] = phi_f * k_f;
    q_r[215] = 0.0;

    /*reaction 217: CH3 + C3H7 <=> 2 C2H5 */
    phi_f = sc[12]*sc[31];
    k_f = k_f_save[216];
    q_f[216] = phi_f * k_f;
    phi_r = sc[25]*sc[25];
    Kc = Kc_save[216];
    k_r = k_f / Kc;
    q_r[216] = phi_r * k_r;

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
    kc[35] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[30]) - (g_RT[6] + g_RT[30]));

    /*reaction 37: H + O2 <=> O + OH */
    kc[36] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));

    /*reaction 38: 2 H + M <=> H2 + M */
    kc[37] = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));

    /*reaction 39: 2 H + H2 <=> 2 H2 */
    kc[38] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));

    /*reaction 40: 2 H + H2O <=> H2 + H2O */
    kc[39] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));

    /*reaction 41: 2 H + CO2 <=> H2 + CO2 */
    kc[40] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[15]) - (g_RT[0] + g_RT[15]));

    /*reaction 42: H + OH + M <=> H2O + M */
    kc[41] = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));

    /*reaction 43: H + HO2 <=> O + H2O */
    kc[42] = exp((g_RT[1] + g_RT[6]) - (g_RT[2] + g_RT[5]));

    /*reaction 44: H + HO2 <=> O2 + H2 */
    kc[43] = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));

    /*reaction 45: H + HO2 <=> 2 OH */
    kc[44] = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));

    /*reaction 46: H + H2O2 <=> HO2 + H2 */
    kc[45] = exp((g_RT[1] + g_RT[7]) - (g_RT[6] + g_RT[0]));

    /*reaction 47: H + H2O2 <=> OH + H2O */
    kc[46] = exp((g_RT[1] + g_RT[7]) - (g_RT[4] + g_RT[5]));

    /*reaction 48: H + CH <=> C + H2 */
    kc[47] = exp((g_RT[1] + g_RT[9]) - (g_RT[8] + g_RT[0]));

    /*reaction 49: H + CH2 (+M) <=> CH3 (+M) */
    kc[48] = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[12]));

    /*reaction 50: H + CH2(S) <=> CH + H2 */
    kc[49] = exp((g_RT[1] + g_RT[11]) - (g_RT[9] + g_RT[0]));

    /*reaction 51: H + CH3 (+M) <=> CH4 (+M) */
    kc[50] = 1.0 / (refC) * exp((g_RT[1] + g_RT[12]) - (g_RT[13]));

    /*reaction 52: H + CH4 <=> CH3 + H2 */
    kc[51] = exp((g_RT[1] + g_RT[13]) - (g_RT[12] + g_RT[0]));

    /*reaction 53: H + HCO (+M) <=> CH2O (+M) */
    kc[52] = 1.0 / (refC) * exp((g_RT[1] + g_RT[16]) - (g_RT[17]));

    /*reaction 54: H + HCO <=> H2 + CO */
    kc[53] = exp((g_RT[1] + g_RT[16]) - (g_RT[0] + g_RT[14]));

    /*reaction 55: H + CH2O (+M) <=> CH2OH (+M) */
    kc[54] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[18]));

    /*reaction 56: H + CH2O (+M) <=> CH3O (+M) */
    kc[55] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[19]));

    /*reaction 57: H + CH2O <=> HCO + H2 */
    kc[56] = exp((g_RT[1] + g_RT[17]) - (g_RT[16] + g_RT[0]));

    /*reaction 58: H + CH2OH (+M) <=> CH3OH (+M) */
    kc[57] = 1.0 / (refC) * exp((g_RT[1] + g_RT[18]) - (g_RT[20]));

    /*reaction 59: H + CH2OH <=> H2 + CH2O */
    kc[58] = exp((g_RT[1] + g_RT[18]) - (g_RT[0] + g_RT[17]));

    /*reaction 60: H + CH2OH <=> OH + CH3 */
    kc[59] = exp((g_RT[1] + g_RT[18]) - (g_RT[4] + g_RT[12]));

    /*reaction 61: H + CH2OH <=> CH2(S) + H2O */
    kc[60] = exp((g_RT[1] + g_RT[18]) - (g_RT[11] + g_RT[5]));

    /*reaction 62: H + CH3O (+M) <=> CH3OH (+M) */
    kc[61] = 1.0 / (refC) * exp((g_RT[1] + g_RT[19]) - (g_RT[20]));

    /*reaction 63: H + CH3O <=> H + CH2OH */
    kc[62] = exp((g_RT[1] + g_RT[19]) - (g_RT[1] + g_RT[18]));

    /*reaction 64: H + CH3O <=> H2 + CH2O */
    kc[63] = exp((g_RT[1] + g_RT[19]) - (g_RT[0] + g_RT[17]));

    /*reaction 65: H + CH3O <=> OH + CH3 */
    kc[64] = exp((g_RT[1] + g_RT[19]) - (g_RT[4] + g_RT[12]));

    /*reaction 66: H + CH3O <=> CH2(S) + H2O */
    kc[65] = exp((g_RT[1] + g_RT[19]) - (g_RT[11] + g_RT[5]));

    /*reaction 67: H + CH3OH <=> CH2OH + H2 */
    kc[66] = exp((g_RT[1] + g_RT[20]) - (g_RT[18] + g_RT[0]));

    /*reaction 68: H + CH3OH <=> CH3O + H2 */
    kc[67] = exp((g_RT[1] + g_RT[20]) - (g_RT[19] + g_RT[0]));

    /*reaction 69: H + C2H (+M) <=> C2H2 (+M) */
    kc[68] = 1.0 / (refC) * exp((g_RT[1] + g_RT[21]) - (g_RT[22]));

    /*reaction 70: H + C2H2 (+M) <=> C2H3 (+M) */
    kc[69] = 1.0 / (refC) * exp((g_RT[1] + g_RT[22]) - (g_RT[23]));

    /*reaction 71: H + C2H3 (+M) <=> C2H4 (+M) */
    kc[70] = 1.0 / (refC) * exp((g_RT[1] + g_RT[23]) - (g_RT[24]));

    /*reaction 72: H + C2H3 <=> H2 + C2H2 */
    kc[71] = exp((g_RT[1] + g_RT[23]) - (g_RT[0] + g_RT[22]));

    /*reaction 73: H + C2H4 (+M) <=> C2H5 (+M) */
    kc[72] = 1.0 / (refC) * exp((g_RT[1] + g_RT[24]) - (g_RT[25]));

    /*reaction 74: H + C2H4 <=> C2H3 + H2 */
    kc[73] = exp((g_RT[1] + g_RT[24]) - (g_RT[23] + g_RT[0]));

    /*reaction 75: H + C2H5 (+M) <=> C2H6 (+M) */
    kc[74] = 1.0 / (refC) * exp((g_RT[1] + g_RT[25]) - (g_RT[26]));

    /*reaction 76: H + C2H5 <=> H2 + C2H4 */
    kc[75] = exp((g_RT[1] + g_RT[25]) - (g_RT[0] + g_RT[24]));

    /*reaction 77: H + C2H6 <=> C2H5 + H2 */
    kc[76] = exp((g_RT[1] + g_RT[26]) - (g_RT[25] + g_RT[0]));

    /*reaction 78: H + HCCO <=> CH2(S) + CO */
    kc[77] = exp((g_RT[1] + g_RT[27]) - (g_RT[11] + g_RT[14]));

    /*reaction 79: H + CH2CO <=> HCCO + H2 */
    kc[78] = exp((g_RT[1] + g_RT[28]) - (g_RT[27] + g_RT[0]));

    /*reaction 80: H + CH2CO <=> CH3 + CO */
    kc[79] = exp((g_RT[1] + g_RT[28]) - (g_RT[12] + g_RT[14]));

    /*reaction 81: H + HCCOH <=> H + CH2CO */
    kc[80] = exp((g_RT[1] + g_RT[29]) - (g_RT[1] + g_RT[28]));

    /*reaction 82: H2 + CO (+M) <=> CH2O (+M) */
    kc[81] = 1.0 / (refC) * exp((g_RT[0] + g_RT[14]) - (g_RT[17]));

    /*reaction 83: OH + H2 <=> H + H2O */
    kc[82] = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));

    /*reaction 84: 2 OH (+M) <=> H2O2 (+M) */
    kc[83] = 1.0 / (refC) * exp((2 * g_RT[4]) - (g_RT[7]));

    /*reaction 85: 2 OH <=> O + H2O */
    kc[84] = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));

    /*reaction 86: OH + HO2 <=> O2 + H2O */
    kc[85] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));

    /*reaction 87: OH + H2O2 <=> HO2 + H2O */
    kc[86] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));

    /*reaction 88: OH + H2O2 <=> HO2 + H2O */
    kc[87] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));

    /*reaction 89: OH + C <=> H + CO */
    kc[88] = exp((g_RT[4] + g_RT[8]) - (g_RT[1] + g_RT[14]));

    /*reaction 90: OH + CH <=> H + HCO */
    kc[89] = exp((g_RT[4] + g_RT[9]) - (g_RT[1] + g_RT[16]));

    /*reaction 91: OH + CH2 <=> H + CH2O */
    kc[90] = exp((g_RT[4] + g_RT[10]) - (g_RT[1] + g_RT[17]));

    /*reaction 92: OH + CH2 <=> CH + H2O */
    kc[91] = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[5]));

    /*reaction 93: OH + CH2(S) <=> H + CH2O */
    kc[92] = exp((g_RT[4] + g_RT[11]) - (g_RT[1] + g_RT[17]));

    /*reaction 94: OH + CH3 (+M) <=> CH3OH (+M) */
    kc[93] = 1.0 / (refC) * exp((g_RT[4] + g_RT[12]) - (g_RT[20]));

    /*reaction 95: OH + CH3 <=> CH2 + H2O */
    kc[94] = exp((g_RT[4] + g_RT[12]) - (g_RT[10] + g_RT[5]));

    /*reaction 96: OH + CH3 <=> CH2(S) + H2O */
    kc[95] = exp((g_RT[4] + g_RT[12]) - (g_RT[11] + g_RT[5]));

    /*reaction 97: OH + CH4 <=> CH3 + H2O */
    kc[96] = exp((g_RT[4] + g_RT[13]) - (g_RT[12] + g_RT[5]));

    /*reaction 98: OH + CO <=> H + CO2 */
    kc[97] = exp((g_RT[4] + g_RT[14]) - (g_RT[1] + g_RT[15]));

    /*reaction 99: OH + HCO <=> H2O + CO */
    kc[98] = exp((g_RT[4] + g_RT[16]) - (g_RT[5] + g_RT[14]));

    /*reaction 100: OH + CH2O <=> HCO + H2O */
    kc[99] = exp((g_RT[4] + g_RT[17]) - (g_RT[16] + g_RT[5]));

    /*reaction 101: OH + CH2OH <=> H2O + CH2O */
    kc[100] = exp((g_RT[4] + g_RT[18]) - (g_RT[5] + g_RT[17]));

    /*reaction 102: OH + CH3O <=> H2O + CH2O */
    kc[101] = exp((g_RT[4] + g_RT[19]) - (g_RT[5] + g_RT[17]));

    /*reaction 103: OH + CH3OH <=> CH2OH + H2O */
    kc[102] = exp((g_RT[4] + g_RT[20]) - (g_RT[18] + g_RT[5]));

    /*reaction 104: OH + CH3OH <=> CH3O + H2O */
    kc[103] = exp((g_RT[4] + g_RT[20]) - (g_RT[19] + g_RT[5]));

    /*reaction 105: OH + C2H <=> H + HCCO */
    kc[104] = exp((g_RT[4] + g_RT[21]) - (g_RT[1] + g_RT[27]));

    /*reaction 106: OH + C2H2 <=> H + CH2CO */
    kc[105] = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[28]));

    /*reaction 107: OH + C2H2 <=> H + HCCOH */
    kc[106] = exp((g_RT[4] + g_RT[22]) - (g_RT[1] + g_RT[29]));

    /*reaction 108: OH + C2H2 <=> C2H + H2O */
    kc[107] = exp((g_RT[4] + g_RT[22]) - (g_RT[21] + g_RT[5]));

    /*reaction 109: OH + C2H2 <=> CH3 + CO */
    kc[108] = exp((g_RT[4] + g_RT[22]) - (g_RT[12] + g_RT[14]));

    /*reaction 110: OH + C2H3 <=> H2O + C2H2 */
    kc[109] = exp((g_RT[4] + g_RT[23]) - (g_RT[5] + g_RT[22]));

    /*reaction 111: OH + C2H4 <=> C2H3 + H2O */
    kc[110] = exp((g_RT[4] + g_RT[24]) - (g_RT[23] + g_RT[5]));

    /*reaction 112: OH + C2H6 <=> C2H5 + H2O */
    kc[111] = exp((g_RT[4] + g_RT[26]) - (g_RT[25] + g_RT[5]));

    /*reaction 113: OH + CH2CO <=> HCCO + H2O */
    kc[112] = exp((g_RT[4] + g_RT[28]) - (g_RT[27] + g_RT[5]));

    /*reaction 114: 2 HO2 <=> O2 + H2O2 */
    kc[113] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));

    /*reaction 115: 2 HO2 <=> O2 + H2O2 */
    kc[114] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));

    /*reaction 116: HO2 + CH2 <=> OH + CH2O */
    kc[115] = exp((g_RT[6] + g_RT[10]) - (g_RT[4] + g_RT[17]));

    /*reaction 117: HO2 + CH3 <=> O2 + CH4 */
    kc[116] = exp((g_RT[6] + g_RT[12]) - (g_RT[3] + g_RT[13]));

    /*reaction 118: HO2 + CH3 <=> OH + CH3O */
    kc[117] = exp((g_RT[6] + g_RT[12]) - (g_RT[4] + g_RT[19]));

    /*reaction 119: HO2 + CO <=> OH + CO2 */
    kc[118] = exp((g_RT[6] + g_RT[14]) - (g_RT[4] + g_RT[15]));

    /*reaction 120: HO2 + CH2O <=> HCO + H2O2 */
    kc[119] = exp((g_RT[6] + g_RT[17]) - (g_RT[16] + g_RT[7]));

    /*reaction 121: C + O2 <=> O + CO */
    kc[120] = exp((g_RT[8] + g_RT[3]) - (g_RT[2] + g_RT[14]));

    /*reaction 122: C + CH2 <=> H + C2H */
    kc[121] = exp((g_RT[8] + g_RT[10]) - (g_RT[1] + g_RT[21]));

    /*reaction 123: C + CH3 <=> H + C2H2 */
    kc[122] = exp((g_RT[8] + g_RT[12]) - (g_RT[1] + g_RT[22]));

    /*reaction 124: CH + O2 <=> O + HCO */
    kc[123] = exp((g_RT[9] + g_RT[3]) - (g_RT[2] + g_RT[16]));

    /*reaction 125: CH + H2 <=> H + CH2 */
    kc[124] = exp((g_RT[9] + g_RT[0]) - (g_RT[1] + g_RT[10]));

    /*reaction 126: CH + H2O <=> H + CH2O */
    kc[125] = exp((g_RT[9] + g_RT[5]) - (g_RT[1] + g_RT[17]));

    /*reaction 127: CH + CH2 <=> H + C2H2 */
    kc[126] = exp((g_RT[9] + g_RT[10]) - (g_RT[1] + g_RT[22]));

    /*reaction 128: CH + CH3 <=> H + C2H3 */
    kc[127] = exp((g_RT[9] + g_RT[12]) - (g_RT[1] + g_RT[23]));

    /*reaction 129: CH + CH4 <=> H + C2H4 */
    kc[128] = exp((g_RT[9] + g_RT[13]) - (g_RT[1] + g_RT[24]));

    /*reaction 130: CH + CO (+M) <=> HCCO (+M) */
    kc[129] = 1.0 / (refC) * exp((g_RT[9] + g_RT[14]) - (g_RT[27]));

    /*reaction 131: CH + CO2 <=> HCO + CO */
    kc[130] = exp((g_RT[9] + g_RT[15]) - (g_RT[16] + g_RT[14]));

    /*reaction 132: CH + CH2O <=> H + CH2CO */
    kc[131] = exp((g_RT[9] + g_RT[17]) - (g_RT[1] + g_RT[28]));

    /*reaction 133: CH + HCCO <=> CO + C2H2 */
    kc[132] = exp((g_RT[9] + g_RT[27]) - (g_RT[14] + g_RT[22]));

    /*reaction 134: CH2 + O2 => OH + H + CO */
    kc[133] = refC * exp((g_RT[10] + g_RT[3]) - (g_RT[4] + g_RT[1] + g_RT[14]));

    /*reaction 135: CH2 + H2 <=> H + CH3 */
    kc[134] = exp((g_RT[10] + g_RT[0]) - (g_RT[1] + g_RT[12]));

    /*reaction 136: 2 CH2 <=> H2 + C2H2 */
    kc[135] = exp((2 * g_RT[10]) - (g_RT[0] + g_RT[22]));

    /*reaction 137: CH2 + CH3 <=> H + C2H4 */
    kc[136] = exp((g_RT[10] + g_RT[12]) - (g_RT[1] + g_RT[24]));

    /*reaction 138: CH2 + CH4 <=> 2 CH3 */
    kc[137] = exp((g_RT[10] + g_RT[13]) - (2 * g_RT[12]));

    /*reaction 139: CH2 + CO (+M) <=> CH2CO (+M) */
    kc[138] = 1.0 / (refC) * exp((g_RT[10] + g_RT[14]) - (g_RT[28]));

    /*reaction 140: CH2 + HCCO <=> C2H3 + CO */
    kc[139] = exp((g_RT[10] + g_RT[27]) - (g_RT[23] + g_RT[14]));

    /*reaction 141: CH2(S) + N2 <=> CH2 + N2 */
    kc[140] = exp((g_RT[11] + g_RT[30]) - (g_RT[10] + g_RT[30]));

    /*reaction 142: CH2(S) + O2 <=> H + OH + CO */
    kc[141] = refC * exp((g_RT[11] + g_RT[3]) - (g_RT[1] + g_RT[4] + g_RT[14]));

    /*reaction 143: CH2(S) + O2 <=> CO + H2O */
    kc[142] = exp((g_RT[11] + g_RT[3]) - (g_RT[14] + g_RT[5]));

    /*reaction 144: CH2(S) + H2 <=> CH3 + H */
    kc[143] = exp((g_RT[11] + g_RT[0]) - (g_RT[12] + g_RT[1]));

    /*reaction 145: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    kc[144] = 1.0 / (refC) * exp((g_RT[11] + g_RT[5]) - (g_RT[20]));

    /*reaction 146: CH2(S) + H2O <=> CH2 + H2O */
    kc[145] = exp((g_RT[11] + g_RT[5]) - (g_RT[10] + g_RT[5]));

    /*reaction 147: CH2(S) + CH3 <=> H + C2H4 */
    kc[146] = exp((g_RT[11] + g_RT[12]) - (g_RT[1] + g_RT[24]));

    /*reaction 148: CH2(S) + CH4 <=> 2 CH3 */
    kc[147] = exp((g_RT[11] + g_RT[13]) - (2 * g_RT[12]));

    /*reaction 149: CH2(S) + CO <=> CH2 + CO */
    kc[148] = exp((g_RT[11] + g_RT[14]) - (g_RT[10] + g_RT[14]));

    /*reaction 150: CH2(S) + CO2 <=> CH2 + CO2 */
    kc[149] = exp((g_RT[11] + g_RT[15]) - (g_RT[10] + g_RT[15]));

    /*reaction 151: CH2(S) + CO2 <=> CO + CH2O */
    kc[150] = exp((g_RT[11] + g_RT[15]) - (g_RT[14] + g_RT[17]));

    /*reaction 152: CH2(S) + C2H6 <=> CH3 + C2H5 */
    kc[151] = exp((g_RT[11] + g_RT[26]) - (g_RT[12] + g_RT[25]));

    /*reaction 153: CH3 + O2 <=> O + CH3O */
    kc[152] = exp((g_RT[12] + g_RT[3]) - (g_RT[2] + g_RT[19]));

    /*reaction 154: CH3 + O2 <=> OH + CH2O */
    kc[153] = exp((g_RT[12] + g_RT[3]) - (g_RT[4] + g_RT[17]));

    /*reaction 155: CH3 + H2O2 <=> HO2 + CH4 */
    kc[154] = exp((g_RT[12] + g_RT[7]) - (g_RT[6] + g_RT[13]));

    /*reaction 156: 2 CH3 (+M) <=> C2H6 (+M) */
    kc[155] = 1.0 / (refC) * exp((2 * g_RT[12]) - (g_RT[26]));

    /*reaction 157: 2 CH3 <=> H + C2H5 */
    kc[156] = exp((2 * g_RT[12]) - (g_RT[1] + g_RT[25]));

    /*reaction 158: CH3 + HCO <=> CH4 + CO */
    kc[157] = exp((g_RT[12] + g_RT[16]) - (g_RT[13] + g_RT[14]));

    /*reaction 159: CH3 + CH2O <=> HCO + CH4 */
    kc[158] = exp((g_RT[12] + g_RT[17]) - (g_RT[16] + g_RT[13]));

    /*reaction 160: CH3 + CH3OH <=> CH2OH + CH4 */
    kc[159] = exp((g_RT[12] + g_RT[20]) - (g_RT[18] + g_RT[13]));

    /*reaction 161: CH3 + CH3OH <=> CH3O + CH4 */
    kc[160] = exp((g_RT[12] + g_RT[20]) - (g_RT[19] + g_RT[13]));

    /*reaction 162: CH3 + C2H4 <=> C2H3 + CH4 */
    kc[161] = exp((g_RT[12] + g_RT[24]) - (g_RT[23] + g_RT[13]));

    /*reaction 163: CH3 + C2H6 <=> C2H5 + CH4 */
    kc[162] = exp((g_RT[12] + g_RT[26]) - (g_RT[25] + g_RT[13]));

    /*reaction 164: HCO + H2O <=> H + CO + H2O */
    kc[163] = refC * exp((g_RT[16] + g_RT[5]) - (g_RT[1] + g_RT[14] + g_RT[5]));

    /*reaction 165: HCO + M <=> H + CO + M */
    kc[164] = refC * exp((g_RT[16]) - (g_RT[1] + g_RT[14]));

    /*reaction 166: HCO + O2 <=> HO2 + CO */
    kc[165] = exp((g_RT[16] + g_RT[3]) - (g_RT[6] + g_RT[14]));

    /*reaction 167: CH2OH + O2 <=> HO2 + CH2O */
    kc[166] = exp((g_RT[18] + g_RT[3]) - (g_RT[6] + g_RT[17]));

    /*reaction 168: CH3O + O2 <=> HO2 + CH2O */
    kc[167] = exp((g_RT[19] + g_RT[3]) - (g_RT[6] + g_RT[17]));

    /*reaction 169: C2H + O2 <=> HCO + CO */
    kc[168] = exp((g_RT[21] + g_RT[3]) - (g_RT[16] + g_RT[14]));

    /*reaction 170: C2H + H2 <=> H + C2H2 */
    kc[169] = exp((g_RT[21] + g_RT[0]) - (g_RT[1] + g_RT[22]));

    /*reaction 171: C2H3 + O2 <=> HCO + CH2O */
    kc[170] = exp((g_RT[23] + g_RT[3]) - (g_RT[16] + g_RT[17]));

    /*reaction 172: C2H4 (+M) <=> H2 + C2H2 (+M) */
    kc[171] = refC * exp((g_RT[24]) - (g_RT[0] + g_RT[22]));

    /*reaction 173: C2H5 + O2 <=> HO2 + C2H4 */
    kc[172] = exp((g_RT[25] + g_RT[3]) - (g_RT[6] + g_RT[24]));

    /*reaction 174: HCCO + O2 <=> OH + 2 CO */
    kc[173] = refC * exp((g_RT[27] + g_RT[3]) - (g_RT[4] + 2 * g_RT[14]));

    /*reaction 175: 2 HCCO <=> 2 CO + C2H2 */
    kc[174] = refC * exp((2 * g_RT[27]) - (2 * g_RT[14] + g_RT[22]));

    /*reaction 176: O + CH3 => H + H2 + CO */
    kc[175] = refC * exp((g_RT[2] + g_RT[12]) - (g_RT[1] + g_RT[0] + g_RT[14]));

    /*reaction 177: O + C2H4 <=> H + CH2CHO */
    kc[176] = exp((g_RT[2] + g_RT[24]) - (g_RT[1] + g_RT[33]));

    /*reaction 178: O + C2H5 <=> H + CH3CHO */
    kc[177] = exp((g_RT[2] + g_RT[25]) - (g_RT[1] + g_RT[34]));

    /*reaction 179: OH + HO2 <=> O2 + H2O */
    kc[178] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));

    /*reaction 180: OH + CH3 => H2 + CH2O */
    kc[179] = exp((g_RT[4] + g_RT[12]) - (g_RT[0] + g_RT[17]));

    /*reaction 181: CH + H2 (+M) <=> CH3 (+M) */
    kc[180] = 1.0 / (refC) * exp((g_RT[9] + g_RT[0]) - (g_RT[12]));

    /*reaction 182: CH2 + O2 => 2 H + CO2 */
    kc[181] = refC * exp((g_RT[10] + g_RT[3]) - (2 * g_RT[1] + g_RT[15]));

    /*reaction 183: CH2 + O2 <=> O + CH2O */
    kc[182] = exp((g_RT[10] + g_RT[3]) - (g_RT[2] + g_RT[17]));

    /*reaction 184: CH2 + CH2 => 2 H + C2H2 */
    kc[183] = refC * exp((g_RT[10] + g_RT[10]) - (2 * g_RT[1] + g_RT[22]));

    /*reaction 185: CH2(S) + H2O => H2 + CH2O */
    kc[184] = exp((g_RT[11] + g_RT[5]) - (g_RT[0] + g_RT[17]));

    /*reaction 186: C2H3 + O2 <=> O + CH2CHO */
    kc[185] = exp((g_RT[23] + g_RT[3]) - (g_RT[2] + g_RT[33]));

    /*reaction 187: C2H3 + O2 <=> HO2 + C2H2 */
    kc[186] = exp((g_RT[23] + g_RT[3]) - (g_RT[6] + g_RT[22]));

    /*reaction 188: O + CH3CHO <=> OH + CH2CHO */
    kc[187] = exp((g_RT[2] + g_RT[34]) - (g_RT[4] + g_RT[33]));

    /*reaction 189: O + CH3CHO => OH + CH3 + CO */
    kc[188] = refC * exp((g_RT[2] + g_RT[34]) - (g_RT[4] + g_RT[12] + g_RT[14]));

    /*reaction 190: O2 + CH3CHO => HO2 + CH3 + CO */
    kc[189] = refC * exp((g_RT[3] + g_RT[34]) - (g_RT[6] + g_RT[12] + g_RT[14]));

    /*reaction 191: H + CH3CHO <=> CH2CHO + H2 */
    kc[190] = exp((g_RT[1] + g_RT[34]) - (g_RT[33] + g_RT[0]));

    /*reaction 192: H + CH3CHO => CH3 + H2 + CO */
    kc[191] = refC * exp((g_RT[1] + g_RT[34]) - (g_RT[12] + g_RT[0] + g_RT[14]));

    /*reaction 193: OH + CH3CHO => CH3 + H2O + CO */
    kc[192] = refC * exp((g_RT[4] + g_RT[34]) - (g_RT[12] + g_RT[5] + g_RT[14]));

    /*reaction 194: HO2 + CH3CHO => CH3 + H2O2 + CO */
    kc[193] = refC * exp((g_RT[6] + g_RT[34]) - (g_RT[12] + g_RT[7] + g_RT[14]));

    /*reaction 195: CH3 + CH3CHO => CH3 + CH4 + CO */
    kc[194] = refC * exp((g_RT[12] + g_RT[34]) - (g_RT[12] + g_RT[13] + g_RT[14]));

    /*reaction 196: H + CH2CO (+M) <=> CH2CHO (+M) */
    kc[195] = 1.0 / (refC) * exp((g_RT[1] + g_RT[28]) - (g_RT[33]));

    /*reaction 197: O + CH2CHO => H + CH2 + CO2 */
    kc[196] = refC * exp((g_RT[2] + g_RT[33]) - (g_RT[1] + g_RT[10] + g_RT[15]));

    /*reaction 198: O2 + CH2CHO => OH + CO + CH2O */
    kc[197] = refC * exp((g_RT[3] + g_RT[33]) - (g_RT[4] + g_RT[14] + g_RT[17]));

    /*reaction 199: O2 + CH2CHO => OH + 2 HCO */
    kc[198] = refC * exp((g_RT[3] + g_RT[33]) - (g_RT[4] + 2 * g_RT[16]));

    /*reaction 200: H + CH2CHO <=> CH3 + HCO */
    kc[199] = exp((g_RT[1] + g_RT[33]) - (g_RT[12] + g_RT[16]));

    /*reaction 201: H + CH2CHO <=> CH2CO + H2 */
    kc[200] = exp((g_RT[1] + g_RT[33]) - (g_RT[28] + g_RT[0]));

    /*reaction 202: OH + CH2CHO <=> H2O + CH2CO */
    kc[201] = exp((g_RT[4] + g_RT[33]) - (g_RT[5] + g_RT[28]));

    /*reaction 203: OH + CH2CHO <=> HCO + CH2OH */
    kc[202] = exp((g_RT[4] + g_RT[33]) - (g_RT[16] + g_RT[18]));

    /*reaction 204: CH3 + C2H5 (+M) <=> C3H8 (+M) */
    kc[203] = 1.0 / (refC) * exp((g_RT[12] + g_RT[25]) - (g_RT[32]));

    /*reaction 205: O + C3H8 <=> OH + C3H7 */
    kc[204] = exp((g_RT[2] + g_RT[32]) - (g_RT[4] + g_RT[31]));

    /*reaction 206: H + C3H8 <=> C3H7 + H2 */
    kc[205] = exp((g_RT[1] + g_RT[32]) - (g_RT[31] + g_RT[0]));

    /*reaction 207: OH + C3H8 <=> C3H7 + H2O */
    kc[206] = exp((g_RT[4] + g_RT[32]) - (g_RT[31] + g_RT[5]));

    /*reaction 208: C3H7 + H2O2 <=> HO2 + C3H8 */
    kc[207] = exp((g_RT[31] + g_RT[7]) - (g_RT[6] + g_RT[32]));

    /*reaction 209: CH3 + C3H8 <=> C3H7 + CH4 */
    kc[208] = exp((g_RT[12] + g_RT[32]) - (g_RT[31] + g_RT[13]));

    /*reaction 210: CH3 + C2H4 (+M) <=> C3H7 (+M) */
    kc[209] = 1.0 / (refC) * exp((g_RT[12] + g_RT[24]) - (g_RT[31]));

    /*reaction 211: O + C3H7 <=> C2H5 + CH2O */
    kc[210] = exp((g_RT[2] + g_RT[31]) - (g_RT[25] + g_RT[17]));

    /*reaction 212: H + C3H7 (+M) <=> C3H8 (+M) */
    kc[211] = 1.0 / (refC) * exp((g_RT[1] + g_RT[31]) - (g_RT[32]));

    /*reaction 213: H + C3H7 <=> CH3 + C2H5 */
    kc[212] = exp((g_RT[1] + g_RT[31]) - (g_RT[12] + g_RT[25]));

    /*reaction 214: OH + C3H7 <=> C2H5 + CH2OH */
    kc[213] = exp((g_RT[4] + g_RT[31]) - (g_RT[25] + g_RT[18]));

    /*reaction 215: HO2 + C3H7 <=> O2 + C3H8 */
    kc[214] = exp((g_RT[6] + g_RT[31]) - (g_RT[3] + g_RT[32]));

    /*reaction 216: HO2 + C3H7 => OH + C2H5 + CH2O */
    kc[215] = refC * exp((g_RT[6] + g_RT[31]) - (g_RT[4] + g_RT[25] + g_RT[17]));

    /*reaction 217: CH3 + C3H7 <=> 2 C2H5 */
    kc[216] = exp((g_RT[12] + g_RT[31]) - (2 * g_RT[25]));

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
        /*species 31: C3H7 */
        species[31] =
            +1.063186300000000e+04 * invT
            -2.007100720000000e+01
            -1.051551800000000e+00 * tc[0]
            -1.299599000000000e-02 * tc[1]
            -3.966756666666666e-07 * tc[2]
            +1.634130750000000e-09 * tc[3]
            -4.686623500000000e-13 * tc[4];
        /*species 32: C3H8 */
        species[32] =
            -1.395852000000000e+04 * invT
            -1.826813719000000e+01
            -9.335538100000000e-01 * tc[0]
            -1.321228950000000e-02 * tc[1]
            -1.017662116666667e-06 * tc[2]
            +1.831458250000000e-09 * tc[3]
            -4.757462650000000e-13 * tc[4];
        /*species 33: CH2CHO */
        species[33] =
            +1.521476600000000e+03 * invT
            -6.149227999999999e+00
            -3.409062000000000e+00 * tc[0]
            -5.369287000000000e-03 * tc[1]
            -3.152486666666667e-07 * tc[2]
            +5.965485833333333e-10 * tc[3]
            -1.433692500000000e-13 * tc[4];
        /*species 34: CH3CHO */
        species[34] =
            -2.157287800000000e+04 * invT
            +6.264436000000000e-01
            -4.729459500000000e+00 * tc[0]
            +1.596642900000000e-03 * tc[1]
            -7.922486833333334e-06 * tc[2]
            +4.788217583333333e-09 * tc[3]
            -1.096555600000000e-12 * tc[4];
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
        /*species 31: C3H7 */
        species[31] =
            +8.298433600000000e+03 * invT
            +2.318287870000000e+01
            -7.702698700000000e+00 * tc[0]
            -8.022101500000000e-03 * tc[1]
            +8.805536666666666e-07 * tc[2]
            -6.358215833333333e-11 * tc[3]
            +1.969614200000000e-15 * tc[4];
        /*species 32: C3H8 */
        species[32] =
            -1.646751600000000e+04 * invT
            +2.542648580000000e+01
            -7.534136800000000e+00 * tc[0]
            -9.436119499999999e-03 * tc[1]
            +1.045308183333333e-06 * tc[2]
            -7.622970750000001e-11 * tc[3]
            +2.391903450000000e-15 * tc[4];
        /*species 33: CH2CHO */
        species[33] =
            +4.903218000000000e+02 * invT
            +1.102092100000000e+01
            -5.975670000000000e+00 * tc[0]
            -4.065295500000000e-03 * tc[1]
            +4.572706666666667e-07 * tc[2]
            -3.391920000000000e-11 * tc[3]
            +1.088008500000000e-15 * tc[4];
        /*species 34: CH3CHO */
        species[34] =
            -2.259312200000000e+04 * invT
            +8.884902499999999e+00
            -5.404110800000000e+00 * tc[0]
            -5.861529500000000e-03 * tc[1]
            +7.043856166666666e-07 * tc[2]
            -5.697704250000000e-11 * tc[3]
            +2.049243150000000e-15 * tc[4];
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
        /*species 31: C3H7 */
        species[31] =
            +1.06318630e+04 * invT
            -2.10710072e+01
            -1.05155180e+00 * tc[0]
            -1.29959900e-02 * tc[1]
            -3.96675667e-07 * tc[2]
            +1.63413075e-09 * tc[3]
            -4.68662350e-13 * tc[4];
        /*species 32: C3H8 */
        species[32] =
            -1.39585200e+04 * invT
            -1.92681372e+01
            -9.33553810e-01 * tc[0]
            -1.32122895e-02 * tc[1]
            -1.01766212e-06 * tc[2]
            +1.83145825e-09 * tc[3]
            -4.75746265e-13 * tc[4];
        /*species 33: CH2CHO */
        species[33] =
            +1.52147660e+03 * invT
            -7.14922800e+00
            -3.40906200e+00 * tc[0]
            -5.36928700e-03 * tc[1]
            -3.15248667e-07 * tc[2]
            +5.96548583e-10 * tc[3]
            -1.43369250e-13 * tc[4];
        /*species 34: CH3CHO */
        species[34] =
            -2.15728780e+04 * invT
            -3.73556400e-01
            -4.72945950e+00 * tc[0]
            +1.59664290e-03 * tc[1]
            -7.92248683e-06 * tc[2]
            +4.78821758e-09 * tc[3]
            -1.09655560e-12 * tc[4];
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
        /*species 31: C3H7 */
        species[31] =
            +8.29843360e+03 * invT
            +2.21828787e+01
            -7.70269870e+00 * tc[0]
            -8.02210150e-03 * tc[1]
            +8.80553667e-07 * tc[2]
            -6.35821583e-11 * tc[3]
            +1.96961420e-15 * tc[4];
        /*species 32: C3H8 */
        species[32] =
            -1.64675160e+04 * invT
            +2.44264858e+01
            -7.53413680e+00 * tc[0]
            -9.43611950e-03 * tc[1]
            +1.04530818e-06 * tc[2]
            -7.62297075e-11 * tc[3]
            +2.39190345e-15 * tc[4];
        /*species 33: CH2CHO */
        species[33] =
            +4.90321800e+02 * invT
            +1.00209210e+01
            -5.97567000e+00 * tc[0]
            -4.06529550e-03 * tc[1]
            +4.57270667e-07 * tc[2]
            -3.39192000e-11 * tc[3]
            +1.08800850e-15 * tc[4];
        /*species 34: CH3CHO */
        species[34] =
            -2.25931220e+04 * invT
            +7.88490250e+00
            -5.40411080e+00 * tc[0]
            -5.86152950e-03 * tc[1]
            +7.04385617e-07 * tc[2]
            -5.69770425e-11 * tc[3]
            +2.04924315e-15 * tc[4];
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
        /*species 31: C3H7 */
        species[31] =
            +5.15518000e-02
            +2.59919800e-02 * tc[1]
            +2.38005400e-06 * tc[2]
            -1.96095690e-08 * tc[3]
            +9.37324700e-12 * tc[4];
        /*species 32: C3H8 */
        species[32] =
            -6.64461900e-02
            +2.64245790e-02 * tc[1]
            +6.10597270e-06 * tc[2]
            -2.19774990e-08 * tc[3]
            +9.51492530e-12 * tc[4];
        /*species 33: CH2CHO */
        species[33] =
            +2.40906200e+00
            +1.07385740e-02 * tc[1]
            +1.89149200e-06 * tc[2]
            -7.15858300e-09 * tc[3]
            +2.86738500e-12 * tc[4];
        /*species 34: CH3CHO */
        species[34] =
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
        /*species 30: N2 */
        species[30] =
            +1.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 31: C3H7 */
        species[31] =
            +6.70269870e+00
            +1.60442030e-02 * tc[1]
            -5.28332200e-06 * tc[2]
            +7.62985900e-10 * tc[3]
            -3.93922840e-14 * tc[4];
        /*species 32: C3H8 */
        species[32] =
            +6.53413680e+00
            +1.88722390e-02 * tc[1]
            -6.27184910e-06 * tc[2]
            +9.14756490e-10 * tc[3]
            -4.78380690e-14 * tc[4];
        /*species 33: CH2CHO */
        species[33] =
            +4.97567000e+00
            +8.13059100e-03 * tc[1]
            -2.74362400e-06 * tc[2]
            +4.07030400e-10 * tc[3]
            -2.17601700e-14 * tc[4];
        /*species 34: CH3CHO */
        species[34] =
            +4.40411080e+00
            +1.17230590e-02 * tc[1]
            -4.22631370e-06 * tc[2]
            +6.83724510e-10 * tc[3]
            -4.09848630e-14 * tc[4];
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
        /*species 31: C3H7 */
        species[31] =
            +1.05155180e+00
            +2.59919800e-02 * tc[1]
            +2.38005400e-06 * tc[2]
            -1.96095690e-08 * tc[3]
            +9.37324700e-12 * tc[4];
        /*species 32: C3H8 */
        species[32] =
            +9.33553810e-01
            +2.64245790e-02 * tc[1]
            +6.10597270e-06 * tc[2]
            -2.19774990e-08 * tc[3]
            +9.51492530e-12 * tc[4];
        /*species 33: CH2CHO */
        species[33] =
            +3.40906200e+00
            +1.07385740e-02 * tc[1]
            +1.89149200e-06 * tc[2]
            -7.15858300e-09 * tc[3]
            +2.86738500e-12 * tc[4];
        /*species 34: CH3CHO */
        species[34] =
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
        /*species 30: N2 */
        species[30] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 31: C3H7 */
        species[31] =
            +7.70269870e+00
            +1.60442030e-02 * tc[1]
            -5.28332200e-06 * tc[2]
            +7.62985900e-10 * tc[3]
            -3.93922840e-14 * tc[4];
        /*species 32: C3H8 */
        species[32] =
            +7.53413680e+00
            +1.88722390e-02 * tc[1]
            -6.27184910e-06 * tc[2]
            +9.14756490e-10 * tc[3]
            -4.78380690e-14 * tc[4];
        /*species 33: CH2CHO */
        species[33] =
            +5.97567000e+00
            +8.13059100e-03 * tc[1]
            -2.74362400e-06 * tc[2]
            +4.07030400e-10 * tc[3]
            -2.17601700e-14 * tc[4];
        /*species 34: CH3CHO */
        species[34] =
            +5.40411080e+00
            +1.17230590e-02 * tc[1]
            -4.22631370e-06 * tc[2]
            +6.83724510e-10 * tc[3]
            -4.09848630e-14 * tc[4];
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
        /*species 31: C3H7 */
        species[31] =
            +5.15518000e-02
            +1.29959900e-02 * tc[1]
            +7.93351333e-07 * tc[2]
            -4.90239225e-09 * tc[3]
            +1.87464940e-12 * tc[4]
            +1.06318630e+04 * invT;
        /*species 32: C3H8 */
        species[32] =
            -6.64461900e-02
            +1.32122895e-02 * tc[1]
            +2.03532423e-06 * tc[2]
            -5.49437475e-09 * tc[3]
            +1.90298506e-12 * tc[4]
            -1.39585200e+04 * invT;
        /*species 33: CH2CHO */
        species[33] =
            +2.40906200e+00
            +5.36928700e-03 * tc[1]
            +6.30497333e-07 * tc[2]
            -1.78964575e-09 * tc[3]
            +5.73477000e-13 * tc[4]
            +1.52147660e+03 * invT;
        /*species 34: CH3CHO */
        species[34] =
            +3.72945950e+00
            -1.59664290e-03 * tc[1]
            +1.58449737e-05 * tc[2]
            -1.43646527e-08 * tc[3]
            +4.38622240e-12 * tc[4]
            -2.15728780e+04 * invT;
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
        /*species 31: C3H7 */
        species[31] =
            +6.70269870e+00
            +8.02210150e-03 * tc[1]
            -1.76110733e-06 * tc[2]
            +1.90746475e-10 * tc[3]
            -7.87845680e-15 * tc[4]
            +8.29843360e+03 * invT;
        /*species 32: C3H8 */
        species[32] =
            +6.53413680e+00
            +9.43611950e-03 * tc[1]
            -2.09061637e-06 * tc[2]
            +2.28689123e-10 * tc[3]
            -9.56761380e-15 * tc[4]
            -1.64675160e+04 * invT;
        /*species 33: CH2CHO */
        species[33] =
            +4.97567000e+00
            +4.06529550e-03 * tc[1]
            -9.14541333e-07 * tc[2]
            +1.01757600e-10 * tc[3]
            -4.35203400e-15 * tc[4]
            +4.90321800e+02 * invT;
        /*species 34: CH3CHO */
        species[34] =
            +4.40411080e+00
            +5.86152950e-03 * tc[1]
            -1.40877123e-06 * tc[2]
            +1.70931128e-10 * tc[3]
            -8.19697260e-15 * tc[4]
            -2.25931220e+04 * invT;
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
        /*species 31: C3H7 */
        species[31] =
            +1.05155180e+00
            +1.29959900e-02 * tc[1]
            +7.93351333e-07 * tc[2]
            -4.90239225e-09 * tc[3]
            +1.87464940e-12 * tc[4]
            +1.06318630e+04 * invT;
        /*species 32: C3H8 */
        species[32] =
            +9.33553810e-01
            +1.32122895e-02 * tc[1]
            +2.03532423e-06 * tc[2]
            -5.49437475e-09 * tc[3]
            +1.90298506e-12 * tc[4]
            -1.39585200e+04 * invT;
        /*species 33: CH2CHO */
        species[33] =
            +3.40906200e+00
            +5.36928700e-03 * tc[1]
            +6.30497333e-07 * tc[2]
            -1.78964575e-09 * tc[3]
            +5.73477000e-13 * tc[4]
            +1.52147660e+03 * invT;
        /*species 34: CH3CHO */
        species[34] =
            +4.72945950e+00
            -1.59664290e-03 * tc[1]
            +1.58449737e-05 * tc[2]
            -1.43646527e-08 * tc[3]
            +4.38622240e-12 * tc[4]
            -2.15728780e+04 * invT;
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
        /*species 31: C3H7 */
        species[31] =
            +7.70269870e+00
            +8.02210150e-03 * tc[1]
            -1.76110733e-06 * tc[2]
            +1.90746475e-10 * tc[3]
            -7.87845680e-15 * tc[4]
            +8.29843360e+03 * invT;
        /*species 32: C3H8 */
        species[32] =
            +7.53413680e+00
            +9.43611950e-03 * tc[1]
            -2.09061637e-06 * tc[2]
            +2.28689123e-10 * tc[3]
            -9.56761380e-15 * tc[4]
            -1.64675160e+04 * invT;
        /*species 33: CH2CHO */
        species[33] =
            +5.97567000e+00
            +4.06529550e-03 * tc[1]
            -9.14541333e-07 * tc[2]
            +1.01757600e-10 * tc[3]
            -4.35203400e-15 * tc[4]
            +4.90321800e+02 * invT;
        /*species 34: CH3CHO */
        species[34] =
            +5.40411080e+00
            +5.86152950e-03 * tc[1]
            -1.40877123e-06 * tc[2]
            +1.70931128e-10 * tc[3]
            -8.19697260e-15 * tc[4]
            -2.25931220e+04 * invT;
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
        /*species 31: C3H7 */
        species[31] =
            +1.05155180e+00 * tc[0]
            +2.59919800e-02 * tc[1]
            +1.19002700e-06 * tc[2]
            -6.53652300e-09 * tc[3]
            +2.34331175e-12 * tc[4]
            +2.11225590e+01 ;
        /*species 32: C3H8 */
        species[32] =
            +9.33553810e-01 * tc[0]
            +2.64245790e-02 * tc[1]
            +3.05298635e-06 * tc[2]
            -7.32583300e-09 * tc[3]
            +2.37873132e-12 * tc[4]
            +1.92016910e+01 ;
        /*species 33: CH2CHO */
        species[33] =
            +3.40906200e+00 * tc[0]
            +1.07385740e-02 * tc[1]
            +9.45746000e-07 * tc[2]
            -2.38619433e-09 * tc[3]
            +7.16846250e-13 * tc[4]
            +9.55829000e+00 ;
        /*species 34: CH3CHO */
        species[34] =
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
        /*species 30: N2 */
        species[30] =
            +2.92664000e+00 * tc[0]
            +1.48797680e-03 * tc[1]
            -2.84238000e-07 * tc[2]
            +3.36567933e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
        /*species 31: C3H7 */
        species[31] =
            +7.70269870e+00 * tc[0]
            +1.60442030e-02 * tc[1]
            -2.64166100e-06 * tc[2]
            +2.54328633e-10 * tc[3]
            -9.84807100e-15 * tc[4]
            -1.54801800e+01 ;
        /*species 32: C3H8 */
        species[32] =
            +7.53413680e+00 * tc[0]
            +1.88722390e-02 * tc[1]
            -3.13592455e-06 * tc[2]
            +3.04918830e-10 * tc[3]
            -1.19595173e-14 * tc[4]
            -1.78923490e+01 ;
        /*species 33: CH2CHO */
        species[33] =
            +5.97567000e+00 * tc[0]
            +8.13059100e-03 * tc[1]
            -1.37181200e-06 * tc[2]
            +1.35676800e-10 * tc[3]
            -5.44004250e-15 * tc[4]
            -5.04525100e+00 ;
        /*species 34: CH3CHO */
        species[34] =
            +5.40411080e+00 * tc[0]
            +1.17230590e-02 * tc[1]
            -2.11315685e-06 * tc[2]
            +2.27908170e-10 * tc[3]
            -1.02462158e-14 * tc[4]
            -3.48079170e+00 ;
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
    wt[31] = 43.089240; /*C3H7 */
    wt[32] = 44.097210; /*C3H8 */
    wt[33] = 43.045610; /*CH2CHO */
    wt[34] = 44.053580; /*CH3CHO */

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
    y[31] = phi[31]*43.089240;   XW += y[31]; /*C3H7 */
    y[32] = phi[32]*44.097210;   XW += y[32]; /*C3H8 */
    y[33] = phi[33]*43.045610;   XW += y[33]; /*CH2CHO */
    y[34] = phi[34]*44.053580;   XW += y[34]; /*CH3CHO */
    for (id = 0; id < 35; ++id) {
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
    phi[31] = y[31]/ 4.30892400e-02; /*C3H7 (wt in kg) */
    phi[32] = y[32]/ 4.40972100e-02; /*C3H8 (wt in kg) */
    phi[33] = y[33]/ 4.30456100e-02; /*CH2CHO (wt in kg) */
    phi[34] = y[34]/ 4.40535800e-02; /*CH3CHO (wt in kg) */

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
    y[31] = c[31] * 43.089240 / (*rho); 
    y[32] = c[32] * 44.097210 / (*rho); 
    y[33] = c[33] * 43.045610 / (*rho); 
    y[34] = c[34] * 44.053580 / (*rho); 

    return;
}


/*ddebdf compatible right hand side of CV burner */
/*rwrk[0] and rwrk[1] should contain rho and ene respectively */
/*working variable phi contains specific mole numbers */
void fecvrhs_(double * time, double * phi, double * phidot, double * rwrk, int * iwrk)
{
    double rho,ene; /*CV Parameters */
    double y[35], wdot[35]; /*temporary storage */
    int i; /*Loop counter */
    double temperature,pressure; /*temporary var */
    rho = rwrk[0];
    ene = rwrk[1];
    fephity_(phi, iwrk, rwrk, y);
    feeytt_(&ene, y, iwrk, rwrk, &temperature);
    CKPY(&rho, &temperature,  y, iwrk, rwrk, &pressure);
    CKWYP(&pressure, &temperature,  y, iwrk, rwrk, wdot);
    for (i=0; i<35; ++i) phidot[i] = wdot[i] / (rho/1000.0); 

    return;
}


/*returns the dimensionality of the cv burner (number of species) */
int fecvdim_()
{
    return 35;
}


/*ddebdf compatible right hand side of ZND solver */
/*rwrk[0] : scaling factor for pressure */
/*rwrk[1] : preshock density (g/cc)  */
/*rwrk[2] : detonation velocity (cm/s)  */
/*solution vector: [P; rho; y0 ... ylast]  */
void fezndrhs_(double * time, double * z, double * zdot, double * rwrk, int * iwrk)
{
    double psc,rho1,udet; /*ZND Parameters */
    double wt[35], hms[35], wdot[35]; /*temporary storage */
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
    for (i=0; i<35; ++i) {
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
    return 38;
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
    if (strcmp(s1, "C3H7")==0) return 31; 
    if (strcmp(s1, "C3H8")==0) return 32; 
    if (strcmp(s1, "CH2CHO")==0) return 33; 
    if (strcmp(s1, "CH3CHO")==0) return 34; 
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
    if (sn==31) return "C3H7"; 
    if (sn==32) return "C3H8"; 
    if (sn==33) return "CH2CHO"; 
    if (sn==34) return "CH3CHO"; 
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
! GRI-Mech Version 3.0 3/12/99  CHEMKIN-II format
! See README30 file at anonymous FTP site unix.sri.com, directory gri;
! WorldWideWeb home page http://www.me.berkeley.edu/gri_mech/ or
! through http://www.gri.org , under 'Basic  Research', 
! for additional information, contacts, and disclaimer
ELEMENTS
O  H  C  N
END
SPECIES
H2      H       O       O2      OH      H2O     HO2     H2O2    
C       CH      CH2     CH2(S)  CH3     CH4     CO      CO2     
HCO     CH2O    CH2OH   CH3O    CH3OH   C2H     C2H2    C2H3    
C2H4    C2H5    C2H6    HCCO    CH2CO   HCCOH   N2      
C3H7    C3H8    CH2CHO  CH3CHO
END
!THERMO
! Insert GRI-Mech thermodynamics here or use in default file
!END
REACTIONS
2O+M<=>O2+M                              1.200E+17   -1.000        .00
H2/ 2.40/ H2O/15.40/ CH4/ 2.00/ CO/ 1.75/ CO2/ 3.60/ C2H6/ 3.00/ 
O+H+M<=>OH+M                             5.000E+17   -1.000        .00
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
O+H2<=>H+OH                              3.870E+04    2.700    6260.00
O+HO2<=>OH+O2                            2.000E+13     .000        .00
O+H2O2<=>OH+HO2                          9.630E+06    2.000    4000.00
O+CH<=>H+CO                              5.700E+13     .000        .00
O+CH2<=>H+HCO                            8.000E+13     .000        .00
O+CH2(S)<=>H2+CO                         1.500E+13     .000        .00
O+CH2(S)<=>H+HCO                         1.500E+13     .000        .00
O+CH3<=>H+CH2O                           5.060E+13     .000        .00
O+CH4<=>OH+CH3                           1.020E+09    1.500    8600.00
O+CO(+M)<=>CO2(+M)                       1.800E+10     .000    2385.00
   LOW/ 6.020E+14     .000    3000.00/
H2/2.00/ O2/6.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/3.50/ C2H6/3.00/ 
O+HCO<=>OH+CO                            3.000E+13     .000        .00
O+HCO<=>H+CO2                            3.000E+13     .000        .00
O+CH2O<=>OH+HCO                          3.900E+13     .000    3540.00
O+CH2OH<=>OH+CH2O                        1.000E+13     .000        .00
O+CH3O<=>OH+CH2O                         1.000E+13     .000        .00
O+CH3OH<=>OH+CH2OH                       3.880E+05    2.500    3100.00
O+CH3OH<=>OH+CH3O                        1.300E+05    2.500    5000.00
O+C2H<=>CH+CO                            5.000E+13     .000        .00
O+C2H2<=>H+HCCO                          1.350E+07    2.000    1900.00
O+C2H2<=>OH+C2H                          4.600E+19   -1.410   28950.00
O+C2H2<=>CO+CH2                          6.940E+06    2.000    1900.00
O+C2H3<=>H+CH2CO                         3.000E+13     .000        .00
O+C2H4<=>CH3+HCO                         1.250E+07    1.830     220.00
O+C2H5<=>CH3+CH2O                        2.240E+13     .000        .00
O+C2H6<=>OH+C2H5                         8.980E+07    1.920    5690.00
O+HCCO<=>H+2CO                           1.000E+14     .000        .00
O+CH2CO<=>OH+HCCO                        1.000E+13     .000    8000.00
O+CH2CO<=>CH2+CO2                        1.750E+12     .000    1350.00
O2+CO<=>O+CO2                            2.500E+12     .000   47800.00
O2+CH2O<=>HO2+HCO                        1.000E+14     .000   40000.00
H+O2+M<=>HO2+M                           2.800E+18    -.860        .00
O2/ .00/ H2O/ .00/ CO/ .75/ CO2/1.50/ C2H6/1.50/  N2/  .00/
H+2O2<=>HO2+O2                           2.080E+19   -1.240        .00
H+O2+H2O<=>HO2+H2O                       11.26E+18    -.760        .00
H+O2+N2<=>HO2+N2                         2.600E+19   -1.240        .00
H+O2<=>O+OH                              2.650E+16    -.6707  17041.00
2H+M<=>H2+M                              1.000E+18   -1.000        .00
H2/ .00/ H2O/ .00/ CH4/2.00/ CO2/ .00/ C2H6/3.00/ 
2H+H2<=>2H2                              9.000E+16    -.600        .00
2H+H2O<=>H2+H2O                          6.000E+19   -1.250        .00
2H+CO2<=>H2+CO2                          5.500E+20   -2.000        .00
H+OH+M<=>H2O+M                           2.200E+22   -2.000        .00
H2/ .73/ H2O/3.65/ CH4/2.00/ C2H6/3.00/ 
H+HO2<=>O+H2O                            3.970E+12     .000     671.00
H+HO2<=>O2+H2                            4.480E+13     .000    1068.00
H+HO2<=>2OH                              0.840E+14     .000     635.00
H+H2O2<=>HO2+H2                          1.210E+07    2.000    5200.00
H+H2O2<=>OH+H2O                          1.000E+13     .000    3600.00
H+CH<=>C+H2                              1.650E+14     .000        .00
H+CH2(+M)<=>CH3(+M)                      6.000E+14     .000        .00
     LOW  /  1.040E+26   -2.760   1600.00/
     TROE/   .5620  91.00  5836.00  8552.00/
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
H+CH2(S)<=>CH+H2                         3.000E+13     .000        .00
H+CH3(+M)<=>CH4(+M)                      13.90E+15    -.534     536.00
     LOW  /  2.620E+33   -4.760   2440.00/
     TROE/   .7830   74.00  2941.00  6964.00 /
H2/2.00/ H2O/6.00/ CH4/3.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
H+CH4<=>CH3+H2                           6.600E+08    1.620   10840.00
H+HCO(+M)<=>CH2O(+M)                     1.090E+12     .480    -260.00
     LOW  /  2.470E+24   -2.570    425.00/
     TROE/   .7824  271.00  2755.00  6570.00 /
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
H+HCO<=>H2+CO                            7.340E+13     .000        .00
H+CH2O(+M)<=>CH2OH(+M)                   5.400E+11     .454    3600.00
     LOW  /  1.270E+32   -4.820   6530.00/
     TROE/   .7187  103.00  1291.00  4160.00 /
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
H+CH2O(+M)<=>CH3O(+M)                    5.400E+11     .454    2600.00
     LOW  /  2.200E+30   -4.800   5560.00/
     TROE/   .7580   94.00  1555.00  4200.00 /
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
H+CH2O<=>HCO+H2                          5.740E+07    1.900    2742.00
H+CH2OH(+M)<=>CH3OH(+M)                  1.055E+12     .500      86.00
     LOW  /  4.360E+31   -4.650   5080.00/
     TROE/   .600  100.00  90000.0  10000.0 /
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
H+CH2OH<=>H2+CH2O                        2.000E+13     .000        .00
H+CH2OH<=>OH+CH3                         1.650E+11     .650    -284.00
H+CH2OH<=>CH2(S)+H2O                     3.280E+13    -.090     610.00
H+CH3O(+M)<=>CH3OH(+M)                   2.430E+12     .515      50.00
     LOW  /  4.660E+41   -7.440   14080.0/
     TROE/   .700  100.00  90000.0 10000.00 /
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
H+CH3O<=>H+CH2OH                         4.150E+07    1.630    1924.00
H+CH3O<=>H2+CH2O                         2.000E+13     .000        .00
H+CH3O<=>OH+CH3                          1.500E+12     .500    -110.00
H+CH3O<=>CH2(S)+H2O                      2.620E+14    -.230    1070.00
H+CH3OH<=>CH2OH+H2                       1.700E+07    2.100    4870.00
H+CH3OH<=>CH3O+H2                        4.200E+06    2.100    4870.00
H+C2H(+M)<=>C2H2(+M)                     1.000E+17   -1.000        .00
     LOW  /  3.750E+33   -4.800   1900.00/
     TROE/   .6464  132.00  1315.00  5566.00 /
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
H+C2H2(+M)<=>C2H3(+M)                    5.600E+12     .000    2400.00
     LOW  /  3.800E+40   -7.270   7220.00/
     TROE/   .7507   98.50  1302.00  4167.00 /
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
H+C2H3(+M)<=>C2H4(+M)                    6.080E+12     .270     280.00
     LOW  /  1.400E+30   -3.860   3320.00/
     TROE/   .7820  207.50  2663.00  6095.00 /
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
H+C2H3<=>H2+C2H2                         3.000E+13     .000        .00
H+C2H4(+M)<=>C2H5(+M)                    0.540E+12     .454    1820.00
     LOW  /  0.600E+42   -7.620   6970.00/
     TROE/   .9753  210.00   984.00  4374.00 /
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
H+C2H4<=>C2H3+H2                         1.325E+06    2.530   12240.00
H+C2H5(+M)<=>C2H6(+M)                    5.210E+17    -.990    1580.00
     LOW  /  1.990E+41   -7.080   6685.00/
     TROE/   .8422  125.00  2219.00  6882.00 /
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
H+C2H5<=>H2+C2H4                         2.000E+12     .000        .00
H+C2H6<=>C2H5+H2                         1.150E+08    1.900    7530.00
H+HCCO<=>CH2(S)+CO                       1.000E+14     .000        .00
H+CH2CO<=>HCCO+H2                        5.000E+13     .000    8000.00
H+CH2CO<=>CH3+CO                         1.130E+13     .000    3428.00
H+HCCOH<=>H+CH2CO                        1.000E+13     .000        .00
H2+CO(+M)<=>CH2O(+M)                     4.300E+07    1.500   79600.00
     LOW  /  5.070E+27   -3.420  84350.00/
     TROE/   .9320  197.00  1540.00 10300.00 /
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
OH+H2<=>H+H2O                            2.160E+08    1.510    3430.00
2OH(+M)<=>H2O2(+M)                       7.400E+13    -.370        .00
     LOW  /  2.300E+18    -.900  -1700.00/
     TROE/   .7346   94.00  1756.00  5182.00 /
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
2OH<=>O+H2O                              3.570E+04    2.400   -2110.00
OH+HO2<=>O2+H2O                          1.450E+13     .000    -500.00
 DUPLICATE
OH+H2O2<=>HO2+H2O                        2.000E+12     .000     427.00
 DUPLICATE
OH+H2O2<=>HO2+H2O                        1.700E+18     .000   29410.00
 DUPLICATE
OH+C<=>H+CO                              5.000E+13     .000        .00
OH+CH<=>H+HCO                            3.000E+13     .000        .00
OH+CH2<=>H+CH2O                          2.000E+13     .000        .00
OH+CH2<=>CH+H2O                          1.130E+07    2.000    3000.00
OH+CH2(S)<=>H+CH2O                       3.000E+13     .000        .00
OH+CH3(+M)<=>CH3OH(+M)                   2.790E+18   -1.430    1330.00
     LOW  /  4.000E+36   -5.920   3140.00/
     TROE/   .4120  195.0  5900.00  6394.00/ 
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
OH+CH3<=>CH2+H2O                         5.600E+07    1.600    5420.00
OH+CH3<=>CH2(S)+H2O                      6.440E+17   -1.340    1417.00
OH+CH4<=>CH3+H2O                         1.000E+08    1.600    3120.00
OH+CO<=>H+CO2                            4.760E+07    1.228      70.00
OH+HCO<=>H2O+CO                          5.000E+13     .000        .00
OH+CH2O<=>HCO+H2O                        3.430E+09    1.180    -447.00
OH+CH2OH<=>H2O+CH2O                      5.000E+12     .000        .00
OH+CH3O<=>H2O+CH2O                       5.000E+12     .000        .00
OH+CH3OH<=>CH2OH+H2O                     1.440E+06    2.000    -840.00
OH+CH3OH<=>CH3O+H2O                      6.300E+06    2.000    1500.00
OH+C2H<=>H+HCCO                          2.000E+13     .000        .00
OH+C2H2<=>H+CH2CO                        2.180E-04    4.500   -1000.00
OH+C2H2<=>H+HCCOH                        5.040E+05    2.300   13500.00
OH+C2H2<=>C2H+H2O                        3.370E+07    2.000   14000.00
OH+C2H2<=>CH3+CO                         4.830E-04    4.000   -2000.00
OH+C2H3<=>H2O+C2H2                       5.000E+12     .000        .00
OH+C2H4<=>C2H3+H2O                       3.600E+06    2.000    2500.00
OH+C2H6<=>C2H5+H2O                       3.540E+06    2.120     870.00
OH+CH2CO<=>HCCO+H2O                      7.500E+12     .000    2000.00
2HO2<=>O2+H2O2                           1.300E+11     .000   -1630.00
 DUPLICATE
2HO2<=>O2+H2O2                           4.200E+14     .000   12000.00
 DUPLICATE
HO2+CH2<=>OH+CH2O                        2.000E+13     .000        .00
HO2+CH3<=>O2+CH4                         1.000E+12     .000        .00
HO2+CH3<=>OH+CH3O                        3.780E+13     .000        .00
HO2+CO<=>OH+CO2                          1.500E+14     .000   23600.00
HO2+CH2O<=>HCO+H2O2                      5.600E+06    2.000   12000.00
C+O2<=>O+CO                              5.800E+13     .000     576.00
C+CH2<=>H+C2H                            5.000E+13     .000        .00
C+CH3<=>H+C2H2                           5.000E+13     .000        .00
CH+O2<=>O+HCO                            6.710E+13     .000        .00
CH+H2<=>H+CH2                            1.080E+14     .000    3110.00
CH+H2O<=>H+CH2O                          5.710E+12     .000    -755.00
CH+CH2<=>H+C2H2                          4.000E+13     .000        .00
CH+CH3<=>H+C2H3                          3.000E+13     .000        .00
CH+CH4<=>H+C2H4                          6.000E+13     .000        .00
CH+CO(+M)<=>HCCO(+M)                     5.000E+13     .000        .00
     LOW  /  2.690E+28   -3.740   1936.00/
     TROE/   .5757  237.00  1652.00  5069.00 /
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
CH+CO2<=>HCO+CO                          1.900E+14     .000   15792.00
CH+CH2O<=>H+CH2CO                        9.460E+13     .000    -515.00
CH+HCCO<=>CO+C2H2                        5.000E+13     .000        .00
CH2+O2=>OH+H+CO                          5.000E+12     .000    1500.00
CH2+H2<=>H+CH3                           5.000E+05    2.000    7230.00
2CH2<=>H2+C2H2                           1.600E+15     .000   11944.00
CH2+CH3<=>H+C2H4                         4.000E+13     .000        .00
CH2+CH4<=>2CH3                           2.460E+06    2.000    8270.00
CH2+CO(+M)<=>CH2CO(+M)                   8.100E+11     .500    4510.00
     LOW  /  2.690E+33   -5.110   7095.00/
     TROE/   .5907  275.00  1226.00  5185.00 /
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
CH2+HCCO<=>C2H3+CO                       3.000E+13     .000        .00
CH2(S)+N2<=>CH2+N2                       1.500E+13     .000     600.00
CH2(S)+O2<=>H+OH+CO                      2.800E+13     .000        .00
CH2(S)+O2<=>CO+H2O                       1.200E+13     .000        .00
CH2(S)+H2<=>CH3+H                        7.000E+13     .000        .00
CH2(S)+H2O(+M)<=>CH3OH(+M)               4.820E+17   -1.160    1145.00
     LOW  /  1.880E+38   -6.360   5040.00/
     TROE/   .6027  208.00  3922.00  10180.0 /
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
CH2(S)+H2O<=>CH2+H2O                     3.000E+13     .000        .00
CH2(S)+CH3<=>H+C2H4                      1.200E+13     .000    -570.00
CH2(S)+CH4<=>2CH3                        1.600E+13     .000    -570.00
CH2(S)+CO<=>CH2+CO                       9.000E+12     .000        .00
CH2(S)+CO2<=>CH2+CO2                     7.000E+12     .000        .00
CH2(S)+CO2<=>CO+CH2O                     1.400E+13     .000        .00
CH2(S)+C2H6<=>CH3+C2H5                   4.000E+13     .000    -550.00
CH3+O2<=>O+CH3O                          3.560E+13     .000   30480.00
CH3+O2<=>OH+CH2O                         2.310E+12     .000   20315.00
CH3+H2O2<=>HO2+CH4                       2.450E+04    2.470    5180.00
2CH3(+M)<=>C2H6(+M)                      6.770E+16   -1.180     654.00
     LOW  /  3.400E+41   -7.030   2762.00/
     TROE/   .6190  73.20  1180.00  9999.00 /
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
2CH3<=>H+C2H5                            6.840E+12     .100   10600.00
CH3+HCO<=>CH4+CO                         2.648E+13     .000        .00
CH3+CH2O<=>HCO+CH4                       3.320E+03    2.810    5860.00
CH3+CH3OH<=>CH2OH+CH4                    3.000E+07    1.500    9940.00
CH3+CH3OH<=>CH3O+CH4                     1.000E+07    1.500    9940.00
CH3+C2H4<=>C2H3+CH4                      2.270E+05    2.000    9200.00
CH3+C2H6<=>C2H5+CH4                      6.140E+06    1.740   10450.00
HCO+H2O<=>H+CO+H2O                       1.500E+18   -1.000   17000.00
HCO+M<=>H+CO+M                           1.870E+17   -1.000   17000.00
H2/2.00/ H2O/ .00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
HCO+O2<=>HO2+CO                          13.45E+12     .000     400.00
CH2OH+O2<=>HO2+CH2O                      1.800E+13     .000     900.00
CH3O+O2<=>HO2+CH2O                       4.280E-13    7.600   -3530.00
C2H+O2<=>HCO+CO                          1.000E+13     .000    -755.00
C2H+H2<=>H+C2H2                          5.680E+10    0.900    1993.00
C2H3+O2<=>HCO+CH2O                       4.580E+16   -1.390    1015.00
C2H4(+M)<=>H2+C2H2(+M)                   8.000E+12     .440   86770.00
     LOW  /  1.580E+51   -9.300  97800.00/
     TROE/   .7345  180.00  1035.00  5417.00 /
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
C2H5+O2<=>HO2+C2H4                       8.400E+11     .000    3875.00
HCCO+O2<=>OH+2CO                         3.200E+12     .000     854.00
2HCCO<=>2CO+C2H2                         1.000E+13     .000        .00
O+CH3=>H+H2+CO                           3.370E+13     .000        .00
O+C2H4<=>H+CH2CHO                        6.700E+06    1.830     220.00
O+C2H5<=>H+CH3CHO                        1.096E+14     .000        .00
OH+HO2<=>O2+H2O                          0.500E+16     .000   17330.00
  DUPLICATE
OH+CH3=>H2+CH2O                          8.000E+09     .500   -1755.00
CH+H2(+M)<=>CH3(+M)                      1.970E+12     .430    -370.00
   LOW/ 4.820E+25  -2.80  590.0 /
   TROE/ .578  122.0  2535.0  9365.0 /
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
CH2+O2=>2H+CO2                           5.800E+12     .000    1500.00
CH2+O2<=>O+CH2O                          2.400E+12     .000    1500.00
CH2+CH2=>2H+C2H2                         2.000E+14     .000   10989.00
CH2(S)+H2O=>H2+CH2O                      6.820E+10     .250    -935.00
C2H3+O2<=>O+CH2CHO                       3.030E+11     .290      11.00
C2H3+O2<=>HO2+C2H2                       1.337E+06    1.610    -384.00
O+CH3CHO<=>OH+CH2CHO                     5.840E+12     .000    1808.00
O+CH3CHO=>OH+CH3+CO                      5.840E+12     .000    1808.00
O2+CH3CHO=>HO2+CH3+CO                    3.010E+13     .000   39150.00
H+CH3CHO<=>CH2CHO+H2                     2.050E+09    1.160    2405.00
H+CH3CHO=>CH3+H2+CO                      2.050E+09    1.160    2405.00
OH+CH3CHO=>CH3+H2O+CO                    2.343E+10    0.730   -1113.00
HO2+CH3CHO=>CH3+H2O2+CO                  3.010E+12     .000   11923.00
CH3+CH3CHO=>CH3+CH4+CO                   2.720E+06    1.770    5920.00
H+CH2CO(+M)<=>CH2CHO(+M)                 4.865E+11    0.422   -1755.00
    LOW/ 1.012E+42  -7.63  3854.0/
    TROE/ 0.465  201.0  1773.0  5333.0 /
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
O+CH2CHO=>H+CH2+CO2                      1.500E+14     .000       .00
O2+CH2CHO=>OH+CO+CH2O                    1.810E+10     .000       .00
O2+CH2CHO=>OH+2HCO                       2.350E+10     .000       .00
H+CH2CHO<=>CH3+HCO                       2.200E+13     .000       .00
H+CH2CHO<=>CH2CO+H2                      1.100E+13     .000       .00
OH+CH2CHO<=>H2O+CH2CO                    1.200E+13     .000       .00
OH+CH2CHO<=>HCO+CH2OH                    3.010E+13     .000       .00
CH3+C2H5(+M)<=>C3H8(+M)                  .9430E+13     .000       .00
     LOW/ 2.710E+74  -16.82  13065.0 /
     TROE/ .1527  291.0  2742.0  7748.0 / 
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
O+C3H8<=>OH+C3H7                         1.930E+05    2.680   3716.00
H+C3H8<=>C3H7+H2                         1.320E+06    2.540   6756.00
OH+C3H8<=>C3H7+H2O                       3.160E+07    1.800    934.00
C3H7+H2O2<=>HO2+C3H8                     3.780E+02    2.720   1500.00
CH3+C3H8<=>C3H7+CH4                      0.903E+00    3.650   7154.00
CH3+C2H4(+M)<=>C3H7(+M)                  2.550E+06    1.600   5700.00
      LOW/ 3.00E+63  -14.6  18170./
      TROE/ .1894  277.0  8748.0  7891.0 / 
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
O+C3H7<=>C2H5+CH2O                       9.640E+13     .000       .00
H+C3H7(+M)<=>C3H8(+M)                    3.613E+13     .000       .00
      LOW/ 4.420E+61  -13.545  11357.0/
      TROE/ .315  369.0  3285.0  6667.0 / 
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ 
H+C3H7<=>CH3+C2H5                        4.060E+06    2.190    890.00
OH+C3H7<=>C2H5+CH2OH                     2.410E+13     .000       .00
HO2+C3H7<=>O2+C3H8                       2.550E+10    0.255   -943.00
HO2+C3H7=>OH+C2H5+CH2O                   2.410E+13     .000       .00
CH3+C3H7<=>2C2H5                         1.927E+13   -0.320       .00
END

\\
\\
\\  This is the therm file
\\
\\
THERMO 
   300.000  1000.000  5000.000
! GRI-Mech Version 3.0 Thermodynamics released 3/12/99
! NASA Polynomial format for CHEMKIN-II
! see README file for disclaimer
O                 L 1/90O   1   00   00   00G   200.000  3500.000  1000.000    1
 2.56942078E+00-8.59741137E-05 4.19484589E-08-1.00177799E-11 1.22833691E-15    2
 2.92175791E+04 4.78433864E+00 3.16826710E+00-3.27931884E-03 6.64306396E-06    3
-6.12806624E-09 2.11265971E-12 2.91222592E+04 2.05193346E+00                   4
O2                TPIS89O   2   00   00   00G   200.000  3500.000  1000.000    1
 3.28253784E+00 1.48308754E-03-7.57966669E-07 2.09470555E-10-2.16717794E-14    2
-1.08845772E+03 5.45323129E+00 3.78245636E+00-2.99673416E-03 9.84730201E-06    3
-9.68129509E-09 3.24372837E-12-1.06394356E+03 3.65767573E+00                   4
H                 L 7/88H   1   00   00   00G   200.000  3500.000  1000.000    1
 2.50000001E+00-2.30842973E-11 1.61561948E-14-4.73515235E-18 4.98197357E-22    2
 2.54736599E+04-4.46682914E-01 2.50000000E+00 7.05332819E-13-1.99591964E-15    3
 2.30081632E-18-9.27732332E-22 2.54736599E+04-4.46682853E-01                   4
H2                TPIS78H   2   00   00   00G   200.000  3500.000  1000.000    1
 3.33727920E+00-4.94024731E-05 4.99456778E-07-1.79566394E-10 2.00255376E-14    2
-9.50158922E+02-3.20502331E+00 2.34433112E+00 7.98052075E-03-1.94781510E-05    3
 2.01572094E-08-7.37611761E-12-9.17935173E+02 6.83010238E-01                   4
OH                RUS 78O   1H   1   00   00G   200.000  3500.000  1000.000    1
 3.09288767E+00 5.48429716E-04 1.26505228E-07-8.79461556E-11 1.17412376E-14    2
 3.85865700E+03 4.47669610E+00 3.99201543E+00-2.40131752E-03 4.61793841E-06    3
-3.88113333E-09 1.36411470E-12 3.61508056E+03-1.03925458E-01                   4
H2O               L 8/89H   2O   1   00   00G   200.000  3500.000  1000.000    1
 3.03399249E+00 2.17691804E-03-1.64072518E-07-9.70419870E-11 1.68200992E-14    2
-3.00042971E+04 4.96677010E+00 4.19864056E+00-2.03643410E-03 6.52040211E-06    3
-5.48797062E-09 1.77197817E-12-3.02937267E+04-8.49032208E-01                   4
HO2               L 5/89H   1O   2   00   00G   200.000  3500.000  1000.000    1
 4.01721090E+00 2.23982013E-03-6.33658150E-07 1.14246370E-10-1.07908535E-14    2
 1.11856713E+02 3.78510215E+00 4.30179801E+00-4.74912051E-03 2.11582891E-05    3
-2.42763894E-08 9.29225124E-12 2.94808040E+02 3.71666245E+00                   4
H2O2              L 7/88H   2O   2   00   00G   200.000  3500.000  1000.000    1
 4.16500285E+00 4.90831694E-03-1.90139225E-06 3.71185986E-10-2.87908305E-14    2
-1.78617877E+04 2.91615662E+00 4.27611269E+00-5.42822417E-04 1.67335701E-05    3
-2.15770813E-08 8.62454363E-12-1.77025821E+04 3.43505074E+00                   4
C                 L11/88C   1   00   00   00G   200.000  3500.000  1000.000    1
 2.49266888E+00 4.79889284E-05-7.24335020E-08 3.74291029E-11-4.87277893E-15    2
 8.54512953E+04 4.80150373E+00 2.55423955E+00-3.21537724E-04 7.33792245E-07    3
-7.32234889E-10 2.66521446E-13 8.54438832E+04 4.53130848E+00                   4
CH                TPIS79C   1H   1   00   00G   200.000  3500.000  1000.000    1
 2.87846473E+00 9.70913681E-04 1.44445655E-07-1.30687849E-10 1.76079383E-14    2
 7.10124364E+04 5.48497999E+00 3.48981665E+00 3.23835541E-04-1.68899065E-06    3
 3.16217327E-09-1.40609067E-12 7.07972934E+04 2.08401108E+00                   4
CH2               L S/93C   1H   2   00   00G   200.000  3500.000  1000.000    1
 2.87410113E+00 3.65639292E-03-1.40894597E-06 2.60179549E-10-1.87727567E-14    2
 4.62636040E+04 6.17119324E+00 3.76267867E+00 9.68872143E-04 2.79489841E-06    3
-3.85091153E-09 1.68741719E-12 4.60040401E+04 1.56253185E+00                   4
CH2(S)            L S/93C   1H   2   00   00G   200.000  3500.000  1000.000    1
 2.29203842E+00 4.65588637E-03-2.01191947E-06 4.17906000E-10-3.39716365E-14    2
 5.09259997E+04 8.62650169E+00 4.19860411E+00-2.36661419E-03 8.23296220E-06    3
-6.68815981E-09 1.94314737E-12 5.04968163E+04-7.69118967E-01                   4
CH3               L11/89C   1H   3   00   00G   200.000  3500.000  1000.000    1
 2.28571772E+00 7.23990037E-03-2.98714348E-06 5.95684644E-10-4.67154394E-14    2
 1.67755843E+04 8.48007179E+00 3.67359040E+00 2.01095175E-03 5.73021856E-06    3
-6.87117425E-09 2.54385734E-12 1.64449988E+04 1.60456433E+00                   4
CH4               L 8/88C   1H   4   00   00G   200.000  3500.000  1000.000    1
 7.48514950E-02 1.33909467E-02-5.73285809E-06 1.22292535E-09-1.01815230E-13    2
-9.46834459E+03 1.84373180E+01 5.14987613E+00-1.36709788E-02 4.91800599E-05    3
-4.84743026E-08 1.66693956E-11-1.02466476E+04-4.64130376E+00                   4
CO                TPIS79C   1O   1   00   00G   200.000  3500.000  1000.000    1
 2.71518561E+00 2.06252743E-03-9.98825771E-07 2.30053008E-10-2.03647716E-14    2
-1.41518724E+04 7.81868772E+00 3.57953347E+00-6.10353680E-04 1.01681433E-06    3
 9.07005884E-10-9.04424499E-13-1.43440860E+04 3.50840928E+00                   4
CO2               L 7/88C   1O   2   00   00G   200.000  3500.000  1000.000    1
 3.85746029E+00 4.41437026E-03-2.21481404E-06 5.23490188E-10-4.72084164E-14    2
-4.87591660E+04 2.27163806E+00 2.35677352E+00 8.98459677E-03-7.12356269E-06    3
 2.45919022E-09-1.43699548E-13-4.83719697E+04 9.90105222E+00                   4
HCO               L12/89H   1C   1O   1   00G   200.000  3500.000  1000.000    1
 2.77217438E+00 4.95695526E-03-2.48445613E-06 5.89161778E-10-5.33508711E-14    2
 4.01191815E+03 9.79834492E+00 4.22118584E+00-3.24392532E-03 1.37799446E-05    3
-1.33144093E-08 4.33768865E-12 3.83956496E+03 3.39437243E+00                   4
CH2O              L 8/88H   2C   1O   1   00G   200.000  3500.000  1000.000    1
 1.76069008E+00 9.20000082E-03-4.42258813E-06 1.00641212E-09-8.83855640E-14    2
-1.39958323E+04 1.36563230E+01 4.79372315E+00-9.90833369E-03 3.73220008E-05    3
-3.79285261E-08 1.31772652E-11-1.43089567E+04 6.02812900E-01                   4
CH2OH             GUNL93C   1H   3O   1   00G   200.000  3500.000  1000.000    1
 3.69266569E+00 8.64576797E-03-3.75101120E-06 7.87234636E-10-6.48554201E-14    2
-3.24250627E+03 5.81043215E+00 3.86388918E+00 5.59672304E-03 5.93271791E-06    3
-1.04532012E-08 4.36967278E-12-3.19391367E+03 5.47302243E+00                   4
CH3O              121686C   1H   3O   1     G  0300.00   3000.00   1000.000    1
 0.03770799E+02 0.07871497E-01-0.02656384E-04 0.03944431E-08-0.02112616E-12    2
 0.12783252E+03 0.02929575E+02 0.02106204E+02 0.07216595E-01 0.05338472E-04    3
-0.07377636E-07 0.02075610E-10 0.09786011E+04 0.13152177E+02                   4
CH3OH             L 8/88C   1H   4O   1   00G   200.000  3500.000  1000.000    1
 1.78970791E+00 1.40938292E-02-6.36500835E-06 1.38171085E-09-1.17060220E-13    2
-2.53748747E+04 1.45023623E+01 5.71539582E+00-1.52309129E-02 6.52441155E-05    3
-7.10806889E-08 2.61352698E-11-2.56427656E+04-1.50409823E+00                   4
C2H               L 1/91C   2H   1   00   00G   200.000  3500.000  1000.000    1
 3.16780652E+00 4.75221902E-03-1.83787077E-06 3.04190252E-10-1.77232770E-14    2
 6.71210650E+04 6.63589475E+00 2.88965733E+00 1.34099611E-02-2.84769501E-05    3
 2.94791045E-08-1.09331511E-11 6.68393932E+04 6.22296438E+00                   4
C2H2              L 1/91C   2H   2   00   00G   200.000  3500.000  1000.000    1
 4.14756964E+00 5.96166664E-03-2.37294852E-06 4.67412171E-10-3.61235213E-14    2
 2.59359992E+04-1.23028121E+00 8.08681094E-01 2.33615629E-02-3.55171815E-05    3
 2.80152437E-08-8.50072974E-12 2.64289807E+04 1.39397051E+01                   4
C2H3              L 2/92C   2H   3   00   00G   200.000  3500.000  1000.000    1
 3.01672400E+00 1.03302292E-02-4.68082349E-06 1.01763288E-09-8.62607041E-14    2
 3.46128739E+04 7.78732378E+00 3.21246645E+00 1.51479162E-03 2.59209412E-05    3
-3.57657847E-08 1.47150873E-11 3.48598468E+04 8.51054025E+00                   4
C2H4              L 1/91C   2H   4   00   00G   200.000  3500.000  1000.000    1
 2.03611116E+00 1.46454151E-02-6.71077915E-06 1.47222923E-09-1.25706061E-13    2
 4.93988614E+03 1.03053693E+01 3.95920148E+00-7.57052247E-03 5.70990292E-05    3
-6.91588753E-08 2.69884373E-11 5.08977593E+03 4.09733096E+00                   4
C2H5              L12/92C   2H   5   00   00G   200.000  3500.000  1000.000    1
 1.95465642E+00 1.73972722E-02-7.98206668E-06 1.75217689E-09-1.49641576E-13    2
 1.28575200E+04 1.34624343E+01 4.30646568E+00-4.18658892E-03 4.97142807E-05    3
-5.99126606E-08 2.30509004E-11 1.28416265E+04 4.70720924E+00                   4
C2H6              L 8/88C   2H   6   00   00G   200.000  3500.000  1000.000    1
 1.07188150E+00 2.16852677E-02-1.00256067E-05 2.21412001E-09-1.90002890E-13    2
-1.14263932E+04 1.51156107E+01 4.29142492E+00-5.50154270E-03 5.99438288E-05    3
-7.08466285E-08 2.68685771E-11-1.15222055E+04 2.66682316E+00                   4
CH2CO             L 5/90C   2H   2O   1   00G   200.000  3500.000  1000.000    1
 4.51129732E+00 9.00359745E-03-4.16939635E-06 9.23345882E-10-7.94838201E-14    2
-7.55105311E+03 6.32247205E-01 2.13583630E+00 1.81188721E-02-1.73947474E-05    3
 9.34397568E-09-2.01457615E-12-7.04291804E+03 1.22156480E+01                   4
HCCO              SRIC91H   1C   2O   1     G  0300.00   4000.00   1000.000    1
 0.56282058E+01 0.40853401E-02-0.15934547E-05 0.28626052E-09-0.19407832E-13    2
 0.19327215E+05-0.39302595E+01 0.22517214E+01 0.17655021E-01-0.23729101E-04    3
 0.17275759E-07-0.50664811E-11 0.20059449E+05 0.12490417E+02                   4
HCCOH              SRI91C   2O   1H   20   0G   300.000  5000.000  1000.000    1
 0.59238291E+01 0.67923600E-02-0.25658564E-05 0.44987841E-09-0.29940101E-13    2
 0.72646260E+04-0.76017742E+01 0.12423733E+01 0.31072201E-01-0.50866864E-04    3
 0.43137131E-07-0.14014594E-10 0.80316143E+04 0.13874319E+02                   4
H2CN               41687H   2C   1N   1     G  0300.00   4000.000  1000.000    1
 0.52097030E+01 0.29692911E-02-0.28555891E-06-0.16355500E-09 0.30432589E-13    2
 0.27677109E+05-0.44444780E+01 0.28516610E+01 0.56952331E-02 0.10711400E-05    3
-0.16226120E-08-0.23511081E-12 0.28637820E+05 0.89927511E+01                   4
HCN               GRI/98H   1C   1N   1    0G   200.000  6000.000  1000.000    1
 0.38022392E+01 0.31464228E-02-0.10632185E-05 0.16619757E-09-0.97997570E-14    2
 0.14407292E+05 0.15754601E+01 0.22589886E+01 0.10051170E-01-0.13351763E-04    3
 0.10092349E-07-0.30089028E-11 0.14712633E+05 0.89164419E+01                   4
HNO               And93 H   1N   1O   1    0G   200.000  6000.000  1000.000    1
 0.29792509E+01 0.34944059E-02-0.78549778E-06 0.57479594E-10-0.19335916E-15    2
 0.11750582E+05 0.86063728E+01 0.45334916E+01-0.56696171E-02 0.18473207E-04    3
-0.17137094E-07 0.55454573E-11 0.11548297E+05 0.17498417E+01                   4
N                 L 6/88N   1   00   00   00G  200.000  6000.000   1000.000    1
 0.24159429E+01 0.17489065E-03-0.11902369E-06 0.30226245E-10-0.20360982E-14    2
 0.56133773E+05 0.46496096E+01 0.25000000E+01 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 0.56104637E+05 0.41939087E+01                   4
NNH               T07/93N   2H   1   00   00G   200.000  6000.000  1000.000    1
 0.37667544E+01 0.28915082E-02-0.10416620E-05 0.16842594E-09-0.10091896E-13    2
 0.28650697E+05 0.44705067E+01 0.43446927E+01-0.48497072E-02 0.20059459E-04    3
-0.21726464E-07 0.79469539E-11 0.28791973E+05 0.29779410E+01                   4
N2O               L 7/88N   2O   1    0    0G   200.000  6000.000  1000.000    1
 0.48230729E+01 0.26270251E-02-0.95850874E-06 0.16000712E-09-0.97752303E-14    2
 0.80734048E+04-0.22017207E+01 0.22571502E+01 0.11304728E-01-0.13671319E-04    3
 0.96819806E-08-0.29307182E-11 0.87417744E+04 0.10757992E+02                   4
NH                And94 N   1H   1    0    0G   200.000  6000.000  1000.000    1
 0.27836928E+01 0.13298430E-02-0.42478047E-06 0.78348501E-10-0.55044470E-14    2
 0.42120848E+05 0.57407799E+01 0.34929085E+01 0.31179198E-03-0.14890484E-05    3
 0.24816442E-08-0.10356967E-11 0.41880629E+05 0.18483278E+01                   4
NH2               And89 N   1H   2    0    0G   200.000  6000.000  1000.000    1
 0.28347421E+01 0.32073082E-02-0.93390804E-06 0.13702953E-09-0.79206144E-14    2
 0.22171957E+05 0.65204163E+01 0.42040029E+01-0.21061385E-02 0.71068348E-05    3
-0.56115197E-08 0.16440717E-11 0.21885910E+05-0.14184248E+00                   4
NH3               J 6/77N   1H   3    0    0G   200.000  6000.000  1000.000    1
 0.26344521E+01 0.56662560E-02-0.17278676E-05 0.23867161E-09-0.12578786E-13    2
-0.65446958E+04 0.65662928E+01 0.42860274E+01-0.46605230E-02 0.21718513E-04    3
-0.22808887E-07 0.82638046E-11-0.67417285E+04-0.62537277E+00                   4
NO                RUS 78N   1O   1    0    0G   200.000  6000.000  1000.000    1
 0.32606056E+01 0.11911043E-02-0.42917048E-06 0.69457669E-10-0.40336099E-14    2
 0.99209746E+04 0.63693027E+01 0.42184763E+01-0.46389760E-02 0.11041022E-04    3
-0.93361354E-08 0.28035770E-11 0.98446230E+04 0.22808464E+01                   4
NO2               L 7/88N   1O   2    0    0G   200.000  6000.000  1000.000    1
 0.48847542E+01 0.21723956E-02-0.82806906E-06 0.15747510E-09-0.10510895E-13    2
 0.23164983E+04-0.11741695E+00 0.39440312E+01-0.15854290E-02 0.16657812E-04    3
-0.20475426E-07 0.78350564E-11 0.28966179E+04 0.63119917E+01                   4
HCNO              BDEA94H   1N   1C   1O   1G   300.000  5000.000  1382.000    1
 6.59860456E+00 3.02778626E-03-1.07704346E-06 1.71666528E-10-1.01439391E-14    2
 1.79661339E+04-1.03306599E+01 2.64727989E+00 1.27505342E-02-1.04794236E-05    3
 4.41432836E-09-7.57521466E-13 1.92990252E+04 1.07332972E+01                   4
HOCN              BDEA94H   1N   1C   1O   1G   300.000  5000.000  1368.000    1
 5.89784885E+00 3.16789393E-03-1.11801064E-06 1.77243144E-10-1.04339177E-14    2
-3.70653331E+03-6.18167825E+00 3.78604952E+00 6.88667922E-03-3.21487864E-06    3
 5.17195767E-10 1.19360788E-14-2.82698400E+03 5.63292162E+00                   4
HNCO              BDEA94H   1N   1C   1O   1G   300.000  5000.000  1478.000    1
 6.22395134E+00 3.17864004E-03-1.09378755E-06 1.70735163E-10-9.95021955E-15    2
-1.66599344E+04-8.38224741E+00 3.63096317E+00 7.30282357E-03-2.28050003E-06    3
-6.61271298E-10 3.62235752E-13-1.55873636E+04 6.19457727E+00                   4
NCO               EA 93 N   1C   1O   1    0G   200.000  6000.000  1000.000    1
 0.51521845E+01 0.23051761E-02-0.88033153E-06 0.14789098E-09-0.90977996E-14    2
 0.14004123E+05-0.25442660E+01 0.28269308E+01 0.88051688E-02-0.83866134E-05    3
 0.48016964E-08-0.13313595E-11 0.14682477E+05 0.95504646E+01                   4
CN                HBH92 C   1N   1    0    0G   200.000  6000.000  1000.000    1
 0.37459805E+01 0.43450775E-04 0.29705984E-06-0.68651806E-10 0.44134173E-14    2
 0.51536188E+05 0.27867601E+01 0.36129351E+01-0.95551327E-03 0.21442977E-05    3
-0.31516323E-09-0.46430356E-12 0.51708340E+05 0.39804995E+01                   4
HCNN              SRI/94C   1N   2H   10   0G   300.000  5000.000  1000.000    1
 0.58946362E+01 0.39895959E-02-0.15982380E-05 0.29249395E-09-0.20094686E-13    2
 0.53452941E+05-0.51030502E+01 0.25243194E+01 0.15960619E-01-0.18816354E-04    3
 0.12125540E-07-0.32357378E-11 0.54261984E+05 0.11675870E+02                   4
N2                121286N   2               G   300.000  5000.000  1000.000    1
 0.02926640E+02 0.14879768E-02-0.05684760E-05 0.10097038E-09-0.06753351E-13    2
-0.09227977E+04 0.05980528E+02 0.03298677E+02 0.14082404E-02-0.03963222E-04    3
 0.05641515E-07-0.02444854E-10-0.10208999E+04 0.03950372E+02                   4
AR                120186AR  1               G   300.000  5000.000  1000.000    1
 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-0.07453750E+04 0.04366000E+02 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-0.07453750E+04 0.04366000E+02                   4
C3H8              L 4/85C   3H   8    0    0G   300.000  5000.000  1000.00     1
 0.75341368E+01 0.18872239E-01-0.62718491E-05 0.91475649E-09-0.47838069E-13    2
-0.16467516E+05-0.17892349E+02 0.93355381E+00 0.26424579E-01 0.61059727E-05    3
-0.21977499E-07 0.95149253E-11-0.13958520E+05 0.19201691E+02                   4
C3H7              L 9/84C   3H   7    0    0G   300.000  5000.000  1000.00     1
 0.77026987E+01 0.16044203E-01-0.52833220E-05 0.76298590E-09-0.39392284E-13    2
 0.82984336E+04-0.15480180E+02 0.10515518E+01 0.25991980E-01 0.23800540E-05    3
-0.19609569E-07 0.93732470E-11 0.10631863E+05 0.21122559E+02                   4
CH3CHO            L 8/88C   2H   4O   1    0G   200.000  6000.000  1000.00     1
 0.54041108E+01 0.11723059E-01-0.42263137E-05 0.68372451E-09-0.40984863E-13    2
-0.22593122E+05-0.34807917E+01 0.47294595E+01-0.31932858E-02 0.47534921E-04    3
-0.57458611E-07 0.21931112E-10-0.21572878E+05 0.41030159E+01                   4
CH2CHO            SAND86O   1H   3C   2     G  0300.00   5000.00  1000.00      1
 0.05975670E+02 0.08130591E-01-0.02743624E-04 0.04070304E-08-0.02176017E-12    2
 0.04903218E+04-0.05045251E+02 0.03409062E+02 0.10738574E-01 0.01891492E-04    3
-0.07158583E-07 0.02867385E-10 0.15214766E+04 0.09558290E+02                   4
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
CH2CHO             2   436.000     3.970     0.000     0.000     2.000
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
