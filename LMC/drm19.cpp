/*  -*- C -*-  */
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
void CKINIT(int * leniwk, int * lenrwk, int * lencwk, int * linc, int * lout, int * ickwrk, double * rckwrk, char * cckwrk );
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double * rval, int * kerr, int lenline);
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
void CKINIT(int * leniwk, int * lenrwk, int * lencwk, int * linc, int * lout, int * ickwrk, double * rckwrk, char * cckwrk )
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
    kname[ 0*lenkname     ] = 'O';
    kname[ 0*lenkname + 1 ] = ' ';

    /* H  */
    kname[ 1*lenkname     ] = 'H';
    kname[ 1*lenkname + 1 ] = ' ';

    /* C  */
    kname[ 2*lenkname     ] = 'C';
    kname[ 2*lenkname + 1 ] = ' ';

    /* N  */
    kname[ 3*lenkname     ] = 'N';
    kname[ 3*lenkname + 1 ] = ' ';

    /* AR  */
    kname[ 4*lenkname     ] = 'A';
    kname[ 4*lenkname + 1 ] = 'R';
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
    kname[ 0*lenkname     ] = 'H';
    kname[ 0*lenkname + 1 ] = '2';
    kname[ 0*lenkname + 2 ] = ' ';

    /* H  */
    kname[ 1*lenkname     ] = 'H';
    kname[ 1*lenkname + 1 ] = ' ';

    /* O  */
    kname[ 2*lenkname     ] = 'O';
    kname[ 2*lenkname + 1 ] = ' ';

    /* O2  */
    kname[ 3*lenkname     ] = 'O';
    kname[ 3*lenkname + 1 ] = '2';
    kname[ 3*lenkname + 2 ] = ' ';

    /* OH  */
    kname[ 4*lenkname     ] = 'O';
    kname[ 4*lenkname + 1 ] = 'H';
    kname[ 4*lenkname + 2 ] = ' ';

    /* H2O  */
    kname[ 5*lenkname     ] = 'H';
    kname[ 5*lenkname + 1 ] = '2';
    kname[ 5*lenkname + 2 ] = 'O';
    kname[ 5*lenkname + 3 ] = ' ';

    /* HO2  */
    kname[ 6*lenkname     ] = 'H';
    kname[ 6*lenkname + 1 ] = 'O';
    kname[ 6*lenkname + 2 ] = '2';
    kname[ 6*lenkname + 3 ] = ' ';

    /* CH2  */
    kname[ 7*lenkname     ] = 'C';
    kname[ 7*lenkname + 1 ] = 'H';
    kname[ 7*lenkname + 2 ] = '2';
    kname[ 7*lenkname + 3 ] = ' ';

    /* CH2(S)  */
    kname[ 8*lenkname     ] = 'C';
    kname[ 8*lenkname + 1 ] = 'H';
    kname[ 8*lenkname + 2 ] = '2';
    kname[ 8*lenkname + 3 ] = '(';
    kname[ 8*lenkname + 4 ] = 'S';
    kname[ 8*lenkname + 5 ] = ')';
    kname[ 8*lenkname + 6 ] = ' ';

    /* CH3  */
    kname[ 9*lenkname     ] = 'C';
    kname[ 9*lenkname + 1 ] = 'H';
    kname[ 9*lenkname + 2 ] = '3';
    kname[ 9*lenkname + 3 ] = ' ';

    /* CH4  */
    kname[ 10*lenkname     ] = 'C';
    kname[ 10*lenkname + 1 ] = 'H';
    kname[ 10*lenkname + 2 ] = '4';
    kname[ 10*lenkname + 3 ] = ' ';

    /* CO  */
    kname[ 11*lenkname     ] = 'C';
    kname[ 11*lenkname + 1 ] = 'O';
    kname[ 11*lenkname + 2 ] = ' ';

    /* CO2  */
    kname[ 12*lenkname     ] = 'C';
    kname[ 12*lenkname + 1 ] = 'O';
    kname[ 12*lenkname + 2 ] = '2';
    kname[ 12*lenkname + 3 ] = ' ';

    /* HCO  */
    kname[ 13*lenkname     ] = 'H';
    kname[ 13*lenkname + 1 ] = 'C';
    kname[ 13*lenkname + 2 ] = 'O';
    kname[ 13*lenkname + 3 ] = ' ';

    /* CH2O  */
    kname[ 14*lenkname     ] = 'C';
    kname[ 14*lenkname + 1 ] = 'H';
    kname[ 14*lenkname + 2 ] = '2';
    kname[ 14*lenkname + 3 ] = 'O';
    kname[ 14*lenkname + 4 ] = ' ';

    /* CH3O  */
    kname[ 15*lenkname     ] = 'C';
    kname[ 15*lenkname + 1 ] = 'H';
    kname[ 15*lenkname + 2 ] = '3';
    kname[ 15*lenkname + 3 ] = 'O';
    kname[ 15*lenkname + 4 ] = ' ';

    /* C2H4  */
    kname[ 16*lenkname     ] = 'C';
    kname[ 16*lenkname + 1 ] = '2';
    kname[ 16*lenkname + 2 ] = 'H';
    kname[ 16*lenkname + 3 ] = '4';
    kname[ 16*lenkname + 4 ] = ' ';

    /* C2H5  */
    kname[ 17*lenkname     ] = 'C';
    kname[ 17*lenkname + 1 ] = '2';
    kname[ 17*lenkname + 2 ] = 'H';
    kname[ 17*lenkname + 3 ] = '5';
    kname[ 17*lenkname + 4 ] = ' ';

    /* C2H6  */
    kname[ 18*lenkname     ] = 'C';
    kname[ 18*lenkname + 1 ] = '2';
    kname[ 18*lenkname + 2 ] = 'H';
    kname[ 18*lenkname + 3 ] = '6';
    kname[ 18*lenkname + 4 ] = ' ';

    /* N2  */
    kname[ 19*lenkname     ] = 'N';
    kname[ 19*lenkname + 1 ] = '2';
    kname[ 19*lenkname + 2 ] = ' ';

    /* AR  */
    kname[ 20*lenkname     ] = 'A';
    kname[ 20*lenkname + 1 ] = 'R';
    kname[ 20*lenkname + 2 ] = ' ';

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
    XW += x[7]*14.026940; /*CH2 */
    XW += x[8]*14.026940; /*CH2(S) */
    XW += x[9]*15.034910; /*CH3 */
    XW += x[10]*16.042880; /*CH4 */
    XW += x[11]*28.010400; /*CO */
    XW += x[12]*44.009800; /*CO2 */
    XW += x[13]*29.018370; /*HCO */
    XW += x[14]*30.026340; /*CH2O */
    XW += x[15]*31.034310; /*CH3O */
    XW += x[16]*28.053880; /*C2H4 */
    XW += x[17]*29.061850; /*C2H5 */
    XW += x[18]*30.069820; /*C2H6 */
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
    YOW += y[7]/14.026940; /*CH2 */
    YOW += y[8]/14.026940; /*CH2(S) */
    YOW += y[9]/15.034910; /*CH3 */
    YOW += y[10]/16.042880; /*CH4 */
    YOW += y[11]/28.010400; /*CO */
    YOW += y[12]/44.009800; /*CO2 */
    YOW += y[13]/29.018370; /*HCO */
    YOW += y[14]/30.026340; /*CH2O */
    YOW += y[15]/31.034310; /*CH3O */
    YOW += y[16]/28.053880; /*C2H4 */
    YOW += y[17]/29.061850; /*C2H5 */
    YOW += y[18]/30.069820; /*C2H6 */
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
    W += c[7]*14.026940; /*CH2 */
    W += c[8]*14.026940; /*CH2(S) */
    W += c[9]*15.034910; /*CH3 */
    W += c[10]*16.042880; /*CH4 */
    W += c[11]*28.010400; /*CO */
    W += c[12]*44.009800; /*CO2 */
    W += c[13]*29.018370; /*HCO */
    W += c[14]*30.026340; /*CH2O */
    W += c[15]*31.034310; /*CH3O */
    W += c[16]*28.053880; /*C2H4 */
    W += c[17]*29.061850; /*C2H5 */
    W += c[18]*30.069820; /*C2H6 */
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
    XW += x[7]*14.026940; /*CH2 */
    XW += x[8]*14.026940; /*CH2(S) */
    XW += x[9]*15.034910; /*CH3 */
    XW += x[10]*16.042880; /*CH4 */
    XW += x[11]*28.010400; /*CO */
    XW += x[12]*44.009800; /*CO2 */
    XW += x[13]*29.018370; /*HCO */
    XW += x[14]*30.026340; /*CH2O */
    XW += x[15]*31.034310; /*CH3O */
    XW += x[16]*28.053880; /*C2H4 */
    XW += x[17]*29.061850; /*C2H5 */
    XW += x[18]*30.069820; /*C2H6 */
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
    YOW += y[7]/14.026940; /*CH2 */
    YOW += y[8]/14.026940; /*CH2(S) */
    YOW += y[9]/15.034910; /*CH3 */
    YOW += y[10]/16.042880; /*CH4 */
    YOW += y[11]/28.010400; /*CO */
    YOW += y[12]/44.009800; /*CO2 */
    YOW += y[13]/29.018370; /*HCO */
    YOW += y[14]/30.026340; /*CH2O */
    YOW += y[15]/31.034310; /*CH3O */
    YOW += y[16]/28.053880; /*C2H4 */
    YOW += y[17]/29.061850; /*C2H5 */
    YOW += y[18]/30.069820; /*C2H6 */
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
    W += c[7]*14.026940; /*CH2 */
    W += c[8]*14.026940; /*CH2(S) */
    W += c[9]*15.034910; /*CH3 */
    W += c[10]*16.042880; /*CH4 */
    W += c[11]*28.010400; /*CO */
    W += c[12]*44.009800; /*CO2 */
    W += c[13]*29.018370; /*HCO */
    W += c[14]*30.026340; /*CH2O */
    W += c[15]*31.034310; /*CH3O */
    W += c[16]*28.053880; /*C2H4 */
    W += c[17]*29.061850; /*C2H5 */
    W += c[18]*30.069820; /*C2H6 */
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
    YOW += y[7]/14.026940; /*CH2 */
    YOW += y[8]/14.026940; /*CH2(S) */
    YOW += y[9]/15.034910; /*CH3 */
    YOW += y[10]/16.042880; /*CH4 */
    YOW += y[11]/28.010400; /*CO */
    YOW += y[12]/44.009800; /*CO2 */
    YOW += y[13]/29.018370; /*HCO */
    YOW += y[14]/30.026340; /*CH2O */
    YOW += y[15]/31.034310; /*CH3O */
    YOW += y[16]/28.053880; /*C2H4 */
    YOW += y[17]/29.061850; /*C2H5 */
    YOW += y[18]/30.069820; /*C2H6 */
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
    XW += x[7]*14.026940; /*CH2 */
    XW += x[8]*14.026940; /*CH2(S) */
    XW += x[9]*15.034910; /*CH3 */
    XW += x[10]*16.042880; /*CH4 */
    XW += x[11]*28.010400; /*CO */
    XW += x[12]*44.009800; /*CO2 */
    XW += x[13]*29.018370; /*HCO */
    XW += x[14]*30.026340; /*CH2O */
    XW += x[15]*31.034310; /*CH3O */
    XW += x[16]*28.053880; /*C2H4 */
    XW += x[17]*29.061850; /*C2H5 */
    XW += x[18]*30.069820; /*C2H6 */
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
    W += c[7]*14.026940; /*CH2 */
    W += c[8]*14.026940; /*CH2(S) */
    W += c[9]*15.034910; /*CH3 */
    W += c[10]*16.042880; /*CH4 */
    W += c[11]*28.010400; /*CO */
    W += c[12]*44.009800; /*CO2 */
    W += c[13]*29.018370; /*HCO */
    W += c[14]*30.026340; /*CH2O */
    W += c[15]*31.034310; /*CH3O */
    W += c[16]*28.053880; /*C2H4 */
    W += c[17]*29.061850; /*C2H5 */
    W += c[18]*30.069820; /*C2H6 */
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
    YOW += y[7]/14.026940; /*CH2 */
    YOW += y[8]/14.026940; /*CH2(S) */
    YOW += y[9]/15.034910; /*CH3 */
    YOW += y[10]/16.042880; /*CH4 */
    YOW += y[11]/28.010400; /*CO */
    YOW += y[12]/44.009800; /*CO2 */
    YOW += y[13]/29.018370; /*HCO */
    YOW += y[14]/30.026340; /*CH2O */
    YOW += y[15]/31.034310; /*CH3O */
    YOW += y[16]/28.053880; /*C2H4 */
    YOW += y[17]/29.061850; /*C2H5 */
    YOW += y[18]/30.069820; /*C2H6 */
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
    x[7] = y[7]/(14.026940*YOW); 
    x[8] = y[8]/(14.026940*YOW); 
    x[9] = y[9]/(15.034910*YOW); 
    x[10] = y[10]/(16.042880*YOW); 
    x[11] = y[11]/(28.010400*YOW); 
    x[12] = y[12]/(44.009800*YOW); 
    x[13] = y[13]/(29.018370*YOW); 
    x[14] = y[14]/(30.026340*YOW); 
    x[15] = y[15]/(31.034310*YOW); 
    x[16] = y[16]/(28.053880*YOW); 
    x[17] = y[17]/(29.061850*YOW); 
    x[18] = y[18]/(30.069820*YOW); 
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
    YOW += y[7]/14.026940; /*CH2 */
    YOW += y[8]/14.026940; /*CH2(S) */
    YOW += y[9]/15.034910; /*CH3 */
    YOW += y[10]/16.042880; /*CH4 */
    YOW += y[11]/28.010400; /*CO */
    YOW += y[12]/44.009800; /*CO2 */
    YOW += y[13]/29.018370; /*HCO */
    YOW += y[14]/30.026340; /*CH2O */
    YOW += y[15]/31.034310; /*CH3O */
    YOW += y[16]/28.053880; /*C2H4 */
    YOW += y[17]/29.061850; /*C2H5 */
    YOW += y[18]/30.069820; /*C2H6 */
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
    c[7] = PWORT * y[7]/14.026940; 
    c[8] = PWORT * y[8]/14.026940; 
    c[9] = PWORT * y[9]/15.034910; 
    c[10] = PWORT * y[10]/16.042880; 
    c[11] = PWORT * y[11]/28.010400; 
    c[12] = PWORT * y[12]/44.009800; 
    c[13] = PWORT * y[13]/29.018370; 
    c[14] = PWORT * y[14]/30.026340; 
    c[15] = PWORT * y[15]/31.034310; 
    c[16] = PWORT * y[16]/28.053880; 
    c[17] = PWORT * y[17]/29.061850; 
    c[18] = PWORT * y[18]/30.069820; 
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
    c[7] = (*rho) * y[7]/14.026940; 
    c[8] = (*rho) * y[8]/14.026940; 
    c[9] = (*rho) * y[9]/15.034910; 
    c[10] = (*rho) * y[10]/16.042880; 
    c[11] = (*rho) * y[11]/28.010400; 
    c[12] = (*rho) * y[12]/44.009800; 
    c[13] = (*rho) * y[13]/29.018370; 
    c[14] = (*rho) * y[14]/30.026340; 
    c[15] = (*rho) * y[15]/31.034310; 
    c[16] = (*rho) * y[16]/28.053880; 
    c[17] = (*rho) * y[17]/29.061850; 
    c[18] = (*rho) * y[18]/30.069820; 
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
    XW += x[7]*14.026940; /*CH2 */
    XW += x[8]*14.026940; /*CH2(S) */
    XW += x[9]*15.034910; /*CH3 */
    XW += x[10]*16.042880; /*CH4 */
    XW += x[11]*28.010400; /*CO */
    XW += x[12]*44.009800; /*CO2 */
    XW += x[13]*29.018370; /*HCO */
    XW += x[14]*30.026340; /*CH2O */
    XW += x[15]*31.034310; /*CH3O */
    XW += x[16]*28.053880; /*C2H4 */
    XW += x[17]*29.061850; /*C2H5 */
    XW += x[18]*30.069820; /*C2H6 */
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
    y[7] = x[7]*14.026940/XW; 
    y[8] = x[8]*14.026940/XW; 
    y[9] = x[9]*15.034910/XW; 
    y[10] = x[10]*16.042880/XW; 
    y[11] = x[11]*28.010400/XW; 
    y[12] = x[12]*44.009800/XW; 
    y[13] = x[13]*29.018370/XW; 
    y[14] = x[14]*30.026340/XW; 
    y[15] = x[15]*31.034310/XW; 
    y[16] = x[16]*28.053880/XW; 
    y[17] = x[17]*29.061850/XW; 
    y[18] = x[18]*30.069820/XW; 
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
    XW += x[7]*14.026940; /*CH2 */
    XW += x[8]*14.026940; /*CH2(S) */
    XW += x[9]*15.034910; /*CH3 */
    XW += x[10]*16.042880; /*CH4 */
    XW += x[11]*28.010400; /*CO */
    XW += x[12]*44.009800; /*CO2 */
    XW += x[13]*29.018370; /*HCO */
    XW += x[14]*30.026340; /*CH2O */
    XW += x[15]*31.034310; /*CH3O */
    XW += x[16]*28.053880; /*C2H4 */
    XW += x[17]*29.061850; /*C2H5 */
    XW += x[18]*30.069820; /*C2H6 */
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
    CW += c[7]*14.026940; /*CH2 */
    CW += c[8]*14.026940; /*CH2(S) */
    CW += c[9]*15.034910; /*CH3 */
    CW += c[10]*16.042880; /*CH4 */
    CW += c[11]*28.010400; /*CO */
    CW += c[12]*44.009800; /*CO2 */
    CW += c[13]*29.018370; /*HCO */
    CW += c[14]*30.026340; /*CH2O */
    CW += c[15]*31.034310; /*CH3O */
    CW += c[16]*28.053880; /*C2H4 */
    CW += c[17]*29.061850; /*C2H5 */
    CW += c[18]*30.069820; /*C2H6 */
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
    y[7] = c[7]*14.026940/CW; 
    y[8] = c[8]*14.026940/CW; 
    y[9] = c[9]*15.034910/CW; 
    y[10] = c[10]*16.042880/CW; 
    y[11] = c[11]*28.010400/CW; 
    y[12] = c[12]*44.009800/CW; 
    y[13] = c[13]*29.018370/CW; 
    y[14] = c[14]*30.026340/CW; 
    y[15] = c[15]*31.034310/CW; 
    y[16] = c[16]*28.053880/CW; 
    y[17] = c[17]*29.061850/CW; 
    y[18] = c[18]*30.069820/CW; 
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
    cvms[0] *= 41241306.784924; /*H2 */
    cvms[1] *= 82482613.569848; /*H */
    cvms[2] *= 5196444.866683; /*O */
    cvms[3] *= 2598222.433341; /*O2 */
    cvms[4] *= 4888468.940230; /*OH */
    cvms[5] *= 4614955.920899; /*H2O */
    cvms[6] *= 2518877.187922; /*HO2 */
    cvms[7] *= 5927165.867966; /*CH2 */
    cvms[8] *= 5927165.867966; /*CH2(S) */
    cvms[9] *= 5529796.985815; /*CH3 */
    cvms[10] *= 5182361.271792; /*CH4 */
    cvms[11] *= 2968183.246223; /*CO */
    cvms[12] *= 1889124.694954; /*CO2 */
    cvms[13] *= 2865081.670680; /*HCO */
    cvms[14] *= 2768902.237169; /*CH2O */
    cvms[15] *= 2678970.468491; /*CH3O */
    cvms[16] *= 2963582.933983; /*C2H4 */
    cvms[17] *= 2860795.166171; /*C2H5 */
    cvms[18] *= 2764898.492908; /*C2H6 */
    cvms[19] *= 2967865.378712; /*N2 */
    cvms[20] *= 2081205.567237; /*AR */
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
    cpms[7] *= 5927165.867966; /*CH2 */
    cpms[8] *= 5927165.867966; /*CH2(S) */
    cpms[9] *= 5529796.985815; /*CH3 */
    cpms[10] *= 5182361.271792; /*CH4 */
    cpms[11] *= 2968183.246223; /*CO */
    cpms[12] *= 1889124.694954; /*CO2 */
    cpms[13] *= 2865081.670680; /*HCO */
    cpms[14] *= 2768902.237169; /*CH2O */
    cpms[15] *= 2678970.468491; /*CH3O */
    cpms[16] *= 2963582.933983; /*C2H4 */
    cpms[17] *= 2860795.166171; /*C2H5 */
    cpms[18] *= 2764898.492908; /*C2H6 */
    cpms[19] *= 2967865.378712; /*N2 */
    cpms[20] *= 2081205.567237; /*AR */
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
    ums[7] *= RT/14.026940; /*CH2 */
    ums[8] *= RT/14.026940; /*CH2(S) */
    ums[9] *= RT/15.034910; /*CH3 */
    ums[10] *= RT/16.042880; /*CH4 */
    ums[11] *= RT/28.010400; /*CO */
    ums[12] *= RT/44.009800; /*CO2 */
    ums[13] *= RT/29.018370; /*HCO */
    ums[14] *= RT/30.026340; /*CH2O */
    ums[15] *= RT/31.034310; /*CH3O */
    ums[16] *= RT/28.053880; /*C2H4 */
    ums[17] *= RT/29.061850; /*C2H5 */
    ums[18] *= RT/30.069820; /*C2H6 */
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
    hms[7] *= RT/14.026940; /*CH2 */
    hms[8] *= RT/14.026940; /*CH2(S) */
    hms[9] *= RT/15.034910; /*CH3 */
    hms[10] *= RT/16.042880; /*CH4 */
    hms[11] *= RT/28.010400; /*CO */
    hms[12] *= RT/44.009800; /*CO2 */
    hms[13] *= RT/29.018370; /*HCO */
    hms[14] *= RT/30.026340; /*CH2O */
    hms[15] *= RT/31.034310; /*CH3O */
    hms[16] *= RT/28.053880; /*C2H4 */
    hms[17] *= RT/29.061850; /*C2H5 */
    hms[18] *= RT/30.069820; /*C2H6 */
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
    gms[7] *= RT/14.026940; /*CH2 */
    gms[8] *= RT/14.026940; /*CH2(S) */
    gms[9] *= RT/15.034910; /*CH3 */
    gms[10] *= RT/16.042880; /*CH4 */
    gms[11] *= RT/28.010400; /*CO */
    gms[12] *= RT/44.009800; /*CO2 */
    gms[13] *= RT/29.018370; /*HCO */
    gms[14] *= RT/30.026340; /*CH2O */
    gms[15] *= RT/31.034310; /*CH3O */
    gms[16] *= RT/28.053880; /*C2H4 */
    gms[17] *= RT/29.061850; /*C2H5 */
    gms[18] *= RT/30.069820; /*C2H6 */
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
    ams[7] *= RT/14.026940; /*CH2 */
    ams[8] *= RT/14.026940; /*CH2(S) */
    ams[9] *= RT/15.034910; /*CH3 */
    ams[10] *= RT/16.042880; /*CH4 */
    ams[11] *= RT/28.010400; /*CO */
    ams[12] *= RT/44.009800; /*CO2 */
    ams[13] *= RT/29.018370; /*HCO */
    ams[14] *= RT/30.026340; /*CH2O */
    ams[15] *= RT/31.034310; /*CH3O */
    ams[16] *= RT/28.053880; /*C2H4 */
    ams[17] *= RT/29.061850; /*C2H5 */
    ams[18] *= RT/30.069820; /*C2H6 */
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
    sms[0] *= 41241306.784924; /*H2 */
    sms[1] *= 82482613.569848; /*H */
    sms[2] *= 5196444.866683; /*O */
    sms[3] *= 2598222.433341; /*O2 */
    sms[4] *= 4888468.940230; /*OH */
    sms[5] *= 4614955.920899; /*H2O */
    sms[6] *= 2518877.187922; /*HO2 */
    sms[7] *= 5927165.867966; /*CH2 */
    sms[8] *= 5927165.867966; /*CH2(S) */
    sms[9] *= 5529796.985815; /*CH3 */
    sms[10] *= 5182361.271792; /*CH4 */
    sms[11] *= 2968183.246223; /*CO */
    sms[12] *= 1889124.694954; /*CO2 */
    sms[13] *= 2865081.670680; /*HCO */
    sms[14] *= 2768902.237169; /*CH2O */
    sms[15] *= 2678970.468491; /*CH3O */
    sms[16] *= 2963582.933983; /*C2H4 */
    sms[17] *= 2860795.166171; /*C2H5 */
    sms[18] *= 2764898.492908; /*C2H6 */
    sms[19] *= 2967865.378712; /*N2 */
    sms[20] *= 2081205.567237; /*AR */
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
    result += cpor[0]*y[0]/2.01594; /*H2 */
    result += cpor[1]*y[1]/1.00797; /*H */
    result += cpor[2]*y[2]/15.9994; /*O */
    result += cpor[3]*y[3]/31.9988; /*O2 */
    result += cpor[4]*y[4]/17.0074; /*OH */
    result += cpor[5]*y[5]/18.0153; /*H2O */
    result += cpor[6]*y[6]/33.0068; /*HO2 */
    result += cpor[7]*y[7]/14.0269; /*CH2 */
    result += cpor[8]*y[8]/14.0269; /*CH2(S) */
    result += cpor[9]*y[9]/15.0349; /*CH3 */
    result += cpor[10]*y[10]/16.0429; /*CH4 */
    result += cpor[11]*y[11]/28.0104; /*CO */
    result += cpor[12]*y[12]/44.0098; /*CO2 */
    result += cpor[13]*y[13]/29.0184; /*HCO */
    result += cpor[14]*y[14]/30.0263; /*CH2O */
    result += cpor[15]*y[15]/31.0343; /*CH3O */
    result += cpor[16]*y[16]/28.0539; /*C2H4 */
    result += cpor[17]*y[17]/29.0618; /*C2H5 */
    result += cpor[18]*y[18]/30.0698; /*C2H6 */
    result += cpor[19]*y[19]/28.0134; /*N2 */
    result += cpor[20]*y[20]/39.948; /*AR */

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
    result += cvor[0]*y[0]/2.01594; /*H2 */
    result += cvor[1]*y[1]/1.00797; /*H */
    result += cvor[2]*y[2]/15.9994; /*O */
    result += cvor[3]*y[3]/31.9988; /*O2 */
    result += cvor[4]*y[4]/17.0074; /*OH */
    result += cvor[5]*y[5]/18.0153; /*H2O */
    result += cvor[6]*y[6]/33.0068; /*HO2 */
    result += cvor[7]*y[7]/14.0269; /*CH2 */
    result += cvor[8]*y[8]/14.0269; /*CH2(S) */
    result += cvor[9]*y[9]/15.0349; /*CH3 */
    result += cvor[10]*y[10]/16.0429; /*CH4 */
    result += cvor[11]*y[11]/28.0104; /*CO */
    result += cvor[12]*y[12]/44.0098; /*CO2 */
    result += cvor[13]*y[13]/29.0184; /*HCO */
    result += cvor[14]*y[14]/30.0263; /*CH2O */
    result += cvor[15]*y[15]/31.0343; /*CH3O */
    result += cvor[16]*y[16]/28.0539; /*C2H4 */
    result += cvor[17]*y[17]/29.0618; /*C2H5 */
    result += cvor[18]*y[18]/30.0698; /*C2H6 */
    result += cvor[19]*y[19]/28.0134; /*N2 */
    result += cvor[20]*y[20]/39.948; /*AR */

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
    result += y[7]*hml[7]/14.026940; /*CH2 */
    result += y[8]*hml[8]/14.026940; /*CH2(S) */
    result += y[9]*hml[9]/15.034910; /*CH3 */
    result += y[10]*hml[10]/16.042880; /*CH4 */
    result += y[11]*hml[11]/28.010400; /*CO */
    result += y[12]*hml[12]/44.009800; /*CO2 */
    result += y[13]*hml[13]/29.018370; /*HCO */
    result += y[14]*hml[14]/30.026340; /*CH2O */
    result += y[15]*hml[15]/31.034310; /*CH3O */
    result += y[16]*hml[16]/28.053880; /*C2H4 */
    result += y[17]*hml[17]/29.061850; /*C2H5 */
    result += y[18]*hml[18]/30.069820; /*C2H6 */
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
    result += y[7]*ums[7]/14.026940; /*CH2 */
    result += y[8]*ums[8]/14.026940; /*CH2(S) */
    result += y[9]*ums[9]/15.034910; /*CH3 */
    result += y[10]*ums[10]/16.042880; /*CH4 */
    result += y[11]*ums[11]/28.010400; /*CO */
    result += y[12]*ums[12]/44.009800; /*CO2 */
    result += y[13]*ums[13]/29.018370; /*HCO */
    result += y[14]*ums[14]/30.026340; /*CH2O */
    result += y[15]*ums[15]/31.034310; /*CH3O */
    result += y[16]*ums[16]/28.053880; /*C2H4 */
    result += y[17]*ums[17]/29.061850; /*C2H5 */
    result += y[18]*ums[18]/30.069820; /*C2H6 */
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
    YOW += y[7]/14.026940; /*CH2 */
    YOW += y[8]/14.026940; /*CH2(S) */
    YOW += y[9]/15.034910; /*CH3 */
    YOW += y[10]/16.042880; /*CH4 */
    YOW += y[11]/28.010400; /*CO */
    YOW += y[12]/44.009800; /*CO2 */
    YOW += y[13]/29.018370; /*HCO */
    YOW += y[14]/30.026340; /*CH2O */
    YOW += y[15]/31.034310; /*CH3O */
    YOW += y[16]/28.053880; /*C2H4 */
    YOW += y[17]/29.061850; /*C2H5 */
    YOW += y[18]/30.069820; /*C2H6 */
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
    x[7] = y[7]/(14.026940*YOW); 
    x[8] = y[8]/(14.026940*YOW); 
    x[9] = y[9]/(15.034910*YOW); 
    x[10] = y[10]/(16.042880*YOW); 
    x[11] = y[11]/(28.010400*YOW); 
    x[12] = y[12]/(44.009800*YOW); 
    x[13] = y[13]/(29.018370*YOW); 
    x[14] = y[14]/(30.026340*YOW); 
    x[15] = y[15]/(31.034310*YOW); 
    x[16] = y[16]/(28.053880*YOW); 
    x[17] = y[17]/(29.061850*YOW); 
    x[18] = y[18]/(30.069820*YOW); 
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
    YOW += y[7]/14.026940; /*CH2 */
    YOW += y[8]/14.026940; /*CH2(S) */
    YOW += y[9]/15.034910; /*CH3 */
    YOW += y[10]/16.042880; /*CH4 */
    YOW += y[11]/28.010400; /*CO */
    YOW += y[12]/44.009800; /*CO2 */
    YOW += y[13]/29.018370; /*HCO */
    YOW += y[14]/30.026340; /*CH2O */
    YOW += y[15]/31.034310; /*CH3O */
    YOW += y[16]/28.053880; /*C2H4 */
    YOW += y[17]/29.061850; /*C2H5 */
    YOW += y[18]/30.069820; /*C2H6 */
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
    x[7] = y[7]/(14.026940*YOW); 
    x[8] = y[8]/(14.026940*YOW); 
    x[9] = y[9]/(15.034910*YOW); 
    x[10] = y[10]/(16.042880*YOW); 
    x[11] = y[11]/(28.010400*YOW); 
    x[12] = y[12]/(44.009800*YOW); 
    x[13] = y[13]/(29.018370*YOW); 
    x[14] = y[14]/(30.026340*YOW); 
    x[15] = y[15]/(31.034310*YOW); 
    x[16] = y[16]/(28.053880*YOW); 
    x[17] = y[17]/(29.061850*YOW); 
    x[18] = y[18]/(30.069820*YOW); 
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
    YOW += y[7]/14.026940; /*CH2 */
    YOW += y[8]/14.026940; /*CH2(S) */
    YOW += y[9]/15.034910; /*CH3 */
    YOW += y[10]/16.042880; /*CH4 */
    YOW += y[11]/28.010400; /*CO */
    YOW += y[12]/44.009800; /*CO2 */
    YOW += y[13]/29.018370; /*HCO */
    YOW += y[14]/30.026340; /*CH2O */
    YOW += y[15]/31.034310; /*CH3O */
    YOW += y[16]/28.053880; /*C2H4 */
    YOW += y[17]/29.061850; /*C2H5 */
    YOW += y[18]/30.069820; /*C2H6 */
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
    x[7] = y[7]/(14.026940*YOW); 
    x[8] = y[8]/(14.026940*YOW); 
    x[9] = y[9]/(15.034910*YOW); 
    x[10] = y[10]/(16.042880*YOW); 
    x[11] = y[11]/(28.010400*YOW); 
    x[12] = y[12]/(44.009800*YOW); 
    x[13] = y[13]/(29.018370*YOW); 
    x[14] = y[14]/(30.026340*YOW); 
    x[15] = y[15]/(31.034310*YOW); 
    x[16] = y[16]/(28.053880*YOW); 
    x[17] = y[17]/(29.061850*YOW); 
    x[18] = y[18]/(30.069820*YOW); 
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
    YOW += y[7]/14.026940; /*CH2 */
    YOW += y[8]/14.026940; /*CH2(S) */
    YOW += y[9]/15.034910; /*CH3 */
    YOW += y[10]/16.042880; /*CH4 */
    YOW += y[11]/28.010400; /*CO */
    YOW += y[12]/44.009800; /*CO2 */
    YOW += y[13]/29.018370; /*HCO */
    YOW += y[14]/30.026340; /*CH2O */
    YOW += y[15]/31.034310; /*CH3O */
    YOW += y[16]/28.053880; /*C2H4 */
    YOW += y[17]/29.061850; /*C2H5 */
    YOW += y[18]/30.069820; /*C2H6 */
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
    c[7] = PWORT * y[7]/14.026940; 
    c[8] = PWORT * y[8]/14.026940; 
    c[9] = PWORT * y[9]/15.034910; 
    c[10] = PWORT * y[10]/16.042880; 
    c[11] = PWORT * y[11]/28.010400; 
    c[12] = PWORT * y[12]/44.009800; 
    c[13] = PWORT * y[13]/29.018370; 
    c[14] = PWORT * y[14]/30.026340; 
    c[15] = PWORT * y[15]/31.034310; 
    c[16] = PWORT * y[16]/28.053880; 
    c[17] = PWORT * y[17]/29.061850; 
    c[18] = PWORT * y[18]/30.069820; 
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
    c[7] = 1e6 * (*rho) * y[7]/14.026940; 
    c[8] = 1e6 * (*rho) * y[8]/14.026940; 
    c[9] = 1e6 * (*rho) * y[9]/15.034910; 
    c[10] = 1e6 * (*rho) * y[10]/16.042880; 
    c[11] = 1e6 * (*rho) * y[11]/28.010400; 
    c[12] = 1e6 * (*rho) * y[12]/44.009800; 
    c[13] = 1e6 * (*rho) * y[13]/29.018370; 
    c[14] = 1e6 * (*rho) * y[14]/30.026340; 
    c[15] = 1e6 * (*rho) * y[15]/31.034310; 
    c[16] = 1e6 * (*rho) * y[16]/28.053880; 
    c[17] = 1e6 * (*rho) * y[17]/29.061850; 
    c[18] = 1e6 * (*rho) * y[18]/30.069820; 
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
    XW += x[7]*14.026940; /*CH2 */
    XW += x[8]*14.026940; /*CH2(S) */
    XW += x[9]*15.034910; /*CH3 */
    XW += x[10]*16.042880; /*CH4 */
    XW += x[11]*28.010400; /*CO */
    XW += x[12]*44.009800; /*CO2 */
    XW += x[13]*29.018370; /*HCO */
    XW += x[14]*30.026340; /*CH2O */
    XW += x[15]*31.034310; /*CH3O */
    XW += x[16]*28.053880; /*C2H4 */
    XW += x[17]*29.061850; /*C2H5 */
    XW += x[18]*30.069820; /*C2H6 */
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


/*Returns the rate of progress for each reaction */
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
    YOW += y[7]/14.026940; /*CH2 */
    YOW += y[8]/14.026940; /*CH2(S) */
    YOW += y[9]/15.034910; /*CH3 */
    YOW += y[10]/16.042880; /*CH4 */
    YOW += y[11]/28.010400; /*CO */
    YOW += y[12]/44.009800; /*CO2 */
    YOW += y[13]/29.018370; /*HCO */
    YOW += y[14]/30.026340; /*CH2O */
    YOW += y[15]/31.034310; /*CH3O */
    YOW += y[16]/28.053880; /*C2H4 */
    YOW += y[17]/29.061850; /*C2H5 */
    YOW += y[18]/30.069820; /*C2H6 */
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
    c[7] = PWORT * y[7]/14.026940; 
    c[8] = PWORT * y[8]/14.026940; 
    c[9] = PWORT * y[9]/15.034910; 
    c[10] = PWORT * y[10]/16.042880; 
    c[11] = PWORT * y[11]/28.010400; 
    c[12] = PWORT * y[12]/44.009800; 
    c[13] = PWORT * y[13]/29.018370; 
    c[14] = PWORT * y[14]/30.026340; 
    c[15] = PWORT * y[15]/31.034310; 
    c[16] = PWORT * y[16]/28.053880; 
    c[17] = PWORT * y[17]/29.061850; 
    c[18] = PWORT * y[18]/30.069820; 
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
    c[7] = 1e6 * (*rho) * y[7]/14.026940; 
    c[8] = 1e6 * (*rho) * y[8]/14.026940; 
    c[9] = 1e6 * (*rho) * y[9]/15.034910; 
    c[10] = 1e6 * (*rho) * y[10]/16.042880; 
    c[11] = 1e6 * (*rho) * y[11]/28.010400; 
    c[12] = 1e6 * (*rho) * y[12]/44.009800; 
    c[13] = 1e6 * (*rho) * y[13]/29.018370; 
    c[14] = 1e6 * (*rho) * y[14]/30.026340; 
    c[15] = 1e6 * (*rho) * y[15]/31.034310; 
    c[16] = 1e6 * (*rho) * y[16]/28.053880; 
    c[17] = 1e6 * (*rho) * y[17]/29.061850; 
    c[18] = 1e6 * (*rho) * y[18]/30.069820; 
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
    XW += x[7]*14.026940; /*CH2 */
    XW += x[8]*14.026940; /*CH2(S) */
    XW += x[9]*15.034910; /*CH3 */
    XW += x[10]*16.042880; /*CH4 */
    XW += x[11]*28.010400; /*CO */
    XW += x[12]*44.009800; /*CO2 */
    XW += x[13]*29.018370; /*HCO */
    XW += x[14]*30.026340; /*CH2O */
    XW += x[15]*31.034310; /*CH3O */
    XW += x[16]*28.053880; /*C2H4 */
    XW += x[17]*29.061850; /*C2H5 */
    XW += x[18]*30.069820; /*C2H6 */
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
    for (id = 0; id < 21 * 84; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: O + H + M <=> OH + M */
    nuki[ 2 * kd     ] += -1 ;
    nuki[     kd     ] += -1 ;
    nuki[ 4 * kd     ] += +1 ;

    /*reaction 2: O + H2 <=> H + OH */
    nuki[ 2 * kd + 1 ] += -1 ;
    nuki[          1 ] += -1 ;
    nuki[     kd + 1 ] += +1 ;
    nuki[ 4 * kd + 1 ] += +1 ;

    /*reaction 3: O + HO2 <=> OH + O2 */
    nuki[ 2 * kd + 2 ] += -1 ;
    nuki[ 6 * kd + 2 ] += -1 ;
    nuki[ 4 * kd + 2 ] += +1 ;
    nuki[ 3 * kd + 2 ] += +1 ;

    /*reaction 4: O + CH2 <=> H + HCO */
    nuki[ 2 * kd + 3 ] += -1 ;
    nuki[ 7 * kd + 3 ] += -1 ;
    nuki[     kd + 3 ] += +1 ;
    nuki[ 13 * kd + 3 ] += +1 ;

    /*reaction 5: O + CH2(S) <=> H + HCO */
    nuki[ 2 * kd + 4 ] += -1 ;
    nuki[ 8 * kd + 4 ] += -1 ;
    nuki[     kd + 4 ] += +1 ;
    nuki[ 13 * kd + 4 ] += +1 ;

    /*reaction 6: O + CH3 <=> H + CH2O */
    nuki[ 2 * kd + 5 ] += -1 ;
    nuki[ 9 * kd + 5 ] += -1 ;
    nuki[     kd + 5 ] += +1 ;
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
    nuki[     kd + 9 ] += +1 ;
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
    nuki[     kd + 16 ] += -1 ;
    nuki[ 3 * kd + 16 ] += -1 ;
    nuki[ 6 * kd + 16 ] += +1 ;

    /*reaction 18: H + 2 O2 <=> HO2 + O2 */
    nuki[     kd + 17 ] += -1 ;
    nuki[ 3 * kd + 17 ] += -2 ;
    nuki[ 6 * kd + 17 ] += +1 ;
    nuki[ 3 * kd + 17 ] += +1 ;

    /*reaction 19: H + O2 + H2O <=> HO2 + H2O */
    nuki[     kd + 18 ] += -1 ;
    nuki[ 3 * kd + 18 ] += -1 ;
    nuki[ 5 * kd + 18 ] += -1 ;
    nuki[ 6 * kd + 18 ] += +1 ;
    nuki[ 5 * kd + 18 ] += +1 ;

    /*reaction 20: H + O2 + N2 <=> HO2 + N2 */
    nuki[     kd + 19 ] += -1 ;
    nuki[ 3 * kd + 19 ] += -1 ;
    nuki[ 19 * kd + 19 ] += -1 ;
    nuki[ 6 * kd + 19 ] += +1 ;
    nuki[ 19 * kd + 19 ] += +1 ;

    /*reaction 21: H + O2 + AR <=> HO2 + AR */
    nuki[     kd + 20 ] += -1 ;
    nuki[ 3 * kd + 20 ] += -1 ;
    nuki[ 20 * kd + 20 ] += -1 ;
    nuki[ 6 * kd + 20 ] += +1 ;
    nuki[ 20 * kd + 20 ] += +1 ;

    /*reaction 22: H + O2 <=> O + OH */
    nuki[     kd + 21 ] += -1 ;
    nuki[ 3 * kd + 21 ] += -1 ;
    nuki[ 2 * kd + 21 ] += +1 ;
    nuki[ 4 * kd + 21 ] += +1 ;

    /*reaction 23: 2 H + M <=> H2 + M */
    nuki[     kd + 22 ] += -2 ;
    nuki[          22 ] += +1 ;

    /*reaction 24: 2 H + H2 <=> 2 H2 */
    nuki[     kd + 23 ] += -2 ;
    nuki[          23 ] += -1 ;
    nuki[          23 ] += +2 ;

    /*reaction 25: 2 H + H2O <=> H2 + H2O */
    nuki[     kd + 24 ] += -2 ;
    nuki[ 5 * kd + 24 ] += -1 ;
    nuki[          24 ] += +1 ;
    nuki[ 5 * kd + 24 ] += +1 ;

    /*reaction 26: 2 H + CO2 <=> H2 + CO2 */
    nuki[     kd + 25 ] += -2 ;
    nuki[ 12 * kd + 25 ] += -1 ;
    nuki[          25 ] += +1 ;
    nuki[ 12 * kd + 25 ] += +1 ;

    /*reaction 27: H + OH + M <=> H2O + M */
    nuki[     kd + 26 ] += -1 ;
    nuki[ 4 * kd + 26 ] += -1 ;
    nuki[ 5 * kd + 26 ] += +1 ;

    /*reaction 28: H + HO2 <=> O2 + H2 */
    nuki[     kd + 27 ] += -1 ;
    nuki[ 6 * kd + 27 ] += -1 ;
    nuki[ 3 * kd + 27 ] += +1 ;
    nuki[          27 ] += +1 ;

    /*reaction 29: H + HO2 <=> 2 OH */
    nuki[     kd + 28 ] += -1 ;
    nuki[ 6 * kd + 28 ] += -1 ;
    nuki[ 4 * kd + 28 ] += +2 ;

    /*reaction 30: H + CH2 (+M) <=> CH3 (+M) */
    nuki[     kd + 29 ] += -1 ;
    nuki[ 7 * kd + 29 ] += -1 ;
    nuki[ 9 * kd + 29 ] += +1 ;

    /*reaction 31: H + CH3 (+M) <=> CH4 (+M) */
    nuki[     kd + 30 ] += -1 ;
    nuki[ 9 * kd + 30 ] += -1 ;
    nuki[ 10 * kd + 30 ] += +1 ;

    /*reaction 32: H + CH4 <=> CH3 + H2 */
    nuki[     kd + 31 ] += -1 ;
    nuki[ 10 * kd + 31 ] += -1 ;
    nuki[ 9 * kd + 31 ] += +1 ;
    nuki[          31 ] += +1 ;

    /*reaction 33: H + HCO (+M) <=> CH2O (+M) */
    nuki[     kd + 32 ] += -1 ;
    nuki[ 13 * kd + 32 ] += -1 ;
    nuki[ 14 * kd + 32 ] += +1 ;

    /*reaction 34: H + HCO <=> H2 + CO */
    nuki[     kd + 33 ] += -1 ;
    nuki[ 13 * kd + 33 ] += -1 ;
    nuki[          33 ] += +1 ;
    nuki[ 11 * kd + 33 ] += +1 ;

    /*reaction 35: H + CH2O (+M) <=> CH3O (+M) */
    nuki[     kd + 34 ] += -1 ;
    nuki[ 14 * kd + 34 ] += -1 ;
    nuki[ 15 * kd + 34 ] += +1 ;

    /*reaction 36: H + CH2O <=> HCO + H2 */
    nuki[     kd + 35 ] += -1 ;
    nuki[ 14 * kd + 35 ] += -1 ;
    nuki[ 13 * kd + 35 ] += +1 ;
    nuki[          35 ] += +1 ;

    /*reaction 37: H + CH3O <=> OH + CH3 */
    nuki[     kd + 36 ] += -1 ;
    nuki[ 15 * kd + 36 ] += -1 ;
    nuki[ 4 * kd + 36 ] += +1 ;
    nuki[ 9 * kd + 36 ] += +1 ;

    /*reaction 38: H + C2H4 (+M) <=> C2H5 (+M) */
    nuki[     kd + 37 ] += -1 ;
    nuki[ 16 * kd + 37 ] += -1 ;
    nuki[ 17 * kd + 37 ] += +1 ;

    /*reaction 39: H + C2H5 (+M) <=> C2H6 (+M) */
    nuki[     kd + 38 ] += -1 ;
    nuki[ 17 * kd + 38 ] += -1 ;
    nuki[ 18 * kd + 38 ] += +1 ;

    /*reaction 40: H + C2H6 <=> C2H5 + H2 */
    nuki[     kd + 39 ] += -1 ;
    nuki[ 18 * kd + 39 ] += -1 ;
    nuki[ 17 * kd + 39 ] += +1 ;
    nuki[          39 ] += +1 ;

    /*reaction 41: H2 + CO (+M) <=> CH2O (+M) */
    nuki[          40 ] += -1 ;
    nuki[ 11 * kd + 40 ] += -1 ;
    nuki[ 14 * kd + 40 ] += +1 ;

    /*reaction 42: OH + H2 <=> H + H2O */
    nuki[ 4 * kd + 41 ] += -1 ;
    nuki[          41 ] += -1 ;
    nuki[     kd + 41 ] += +1 ;
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
    nuki[     kd + 44 ] += +1 ;
    nuki[ 14 * kd + 44 ] += +1 ;

    /*reaction 46: OH + CH2(S) <=> H + CH2O */
    nuki[ 4 * kd + 45 ] += -1 ;
    nuki[ 8 * kd + 45 ] += -1 ;
    nuki[     kd + 45 ] += +1 ;
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
    nuki[     kd + 49 ] += +1 ;
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
    nuki[          58 ] += -1 ;
    nuki[     kd + 58 ] += +1 ;
    nuki[ 9 * kd + 58 ] += +1 ;

    /*reaction 60: CH2 + CH3 <=> H + C2H4 */
    nuki[ 7 * kd + 59 ] += -1 ;
    nuki[ 9 * kd + 59 ] += -1 ;
    nuki[     kd + 59 ] += +1 ;
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
    nuki[     kd + 63 ] += +1 ;
    nuki[ 4 * kd + 63 ] += +1 ;
    nuki[ 11 * kd + 63 ] += +1 ;

    /*reaction 65: CH2(S) + O2 <=> CO + H2O */
    nuki[ 8 * kd + 64 ] += -1 ;
    nuki[ 3 * kd + 64 ] += -1 ;
    nuki[ 11 * kd + 64 ] += +1 ;
    nuki[ 5 * kd + 64 ] += +1 ;

    /*reaction 66: CH2(S) + H2 <=> CH3 + H */
    nuki[ 8 * kd + 65 ] += -1 ;
    nuki[          65 ] += -1 ;
    nuki[ 9 * kd + 65 ] += +1 ;
    nuki[     kd + 65 ] += +1 ;

    /*reaction 67: CH2(S) + H2O <=> CH2 + H2O */
    nuki[ 8 * kd + 66 ] += -1 ;
    nuki[ 5 * kd + 66 ] += -1 ;
    nuki[ 7 * kd + 66 ] += +1 ;
    nuki[ 5 * kd + 66 ] += +1 ;

    /*reaction 68: CH2(S) + CH3 <=> H + C2H4 */
    nuki[ 8 * kd + 67 ] += -1 ;
    nuki[ 9 * kd + 67 ] += -1 ;
    nuki[     kd + 67 ] += +1 ;
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
    nuki[     kd + 75 ] += +1 ;
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
    nuki[     kd + 79 ] += +1 ;
    nuki[ 11 * kd + 79 ] += +1 ;
    nuki[ 5 * kd + 79 ] += +1 ;

    /*reaction 81: HCO + M <=> H + CO + M */
    nuki[ 13 * kd + 80 ] += -1 ;
    nuki[     kd + 80 ] += +1 ;
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
    for (id = 0; id < kd * 21; ++ id) {
         ncf[id] = 0; 
    }

    /*H2 */
    ncf[          1 ] = 2; /*H */

    /*H */
    ncf[     kd + 1 ] = 1; /*H */

    /*O */
    ncf[ 2 * kd     ] = 1; /*O */

    /*O2 */
    ncf[ 3 * kd     ] = 2; /*O */

    /*OH */
    ncf[ 4 * kd     ] = 1; /*O */
    ncf[ 4 * kd + 1 ] = 1; /*H */

    /*H2O */
    ncf[ 5 * kd + 1 ] = 2; /*H */
    ncf[ 5 * kd     ] = 1; /*O */

    /*HO2 */
    ncf[ 6 * kd + 1 ] = 1; /*H */
    ncf[ 6 * kd     ] = 2; /*O */

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
    ncf[ 11 * kd     ] = 1; /*O */

    /*CO2 */
    ncf[ 12 * kd + 2 ] = 1; /*C */
    ncf[ 12 * kd     ] = 2; /*O */

    /*HCO */
    ncf[ 13 * kd + 1 ] = 1; /*H */
    ncf[ 13 * kd + 2 ] = 1; /*C */
    ncf[ 13 * kd     ] = 1; /*O */

    /*CH2O */
    ncf[ 14 * kd + 1 ] = 2; /*H */
    ncf[ 14 * kd + 2 ] = 1; /*C */
    ncf[ 14 * kd     ] = 1; /*O */

    /*CH3O */
    ncf[ 15 * kd + 2 ] = 1; /*C */
    ncf[ 15 * kd + 1 ] = 3; /*H */
    ncf[ 15 * kd     ] = 1; /*O */

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


/*compute the production rate for each species */
void productionRate(double * wdot, double * sc, double T)
{
    double qdot;

    static double T_old = -1, k_f_old[84], Kc_old[84];

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

    if (T != T_old)
    {
        T_old = T;

        k_f_old[0] = 1e-12 * 5e+17*exp(-1*tc[0]);
        k_f_old[1] = 1e-06 * 50000*exp(2.67*tc[0]-3165.58/tc[1]);
        k_f_old[2] = 1e-06 * 2e+13;
        k_f_old[3] = 1e-06 * 8e+13;
        k_f_old[4] = 1e-06 * 1.5e+13;
        k_f_old[5] = 1e-06 * 8.43e+13;
        k_f_old[6] = 1e-06 * 1.02e+09*exp(1.5*tc[0]-4328.13/tc[1]);
        k_f_old[7] = 1e-12 * 6.02e+14*exp(-1509.81/tc[1]);
        k_f_old[8] = 1e-06 * 3e+13;
        k_f_old[9] = 1e-06 * 3e+13;
        k_f_old[10] = 1e-06 * 3.9e+13*exp(-1781.58/tc[1]);
        k_f_old[11] = 1e-06 * 1.92e+07*exp(1.83*tc[0]-110.72/tc[1]);
        k_f_old[12] = 1e-06 * 1.32e+14;
        k_f_old[13] = 1e-06 * 8.98e+07*exp(1.92*tc[0]-2863.61/tc[1]);
        k_f_old[14] = 1e-06 * 2.5e+12*exp(-24056.4/tc[1]);
        k_f_old[15] = 1e-06 * 1e+14*exp(-20130.9/tc[1]);
        k_f_old[16] = 1e-12 * 2.8e+18*exp(-0.86*tc[0]);
        k_f_old[17] = 1e-12 * 3e+20*exp(-1.72*tc[0]);
        k_f_old[18] = 1e-12 * 9.38e+18*exp(-0.76*tc[0]);
        k_f_old[19] = 1e-12 * 3.75e+20*exp(-1.72*tc[0]);
        k_f_old[20] = 1e-12 * 7e+17*exp(-0.8*tc[0]);
        k_f_old[21] = 1e-06 * 8.3e+13*exp(-7253.65/tc[1]);
        k_f_old[22] = 1e-12 * 1e+18*exp(-1*tc[0]);
        k_f_old[23] = 1e-12 * 9e+16*exp(-0.6*tc[0]);
        k_f_old[24] = 1e-12 * 6e+19*exp(-1.25*tc[0]);
        k_f_old[25] = 1e-12 * 5.5e+20*exp(-2*tc[0]);
        k_f_old[26] = 1e-12 * 2.2e+22*exp(-2*tc[0]);
        k_f_old[27] = 1e-06 * 2.8e+13*exp(-537.494/tc[1]);
        k_f_old[28] = 1e-06 * 1.34e+14*exp(-319.577/tc[1]);
        k_f_old[29] = 1e-06 * 2.5e+16*exp(-0.8*tc[0]);
        k_f_old[30] = 1e-06 * 1.27e+16*exp(-0.63*tc[0]-192.753/tc[1]);
        k_f_old[31] = 1e-06 * 6.6e+08*exp(1.62*tc[0]-5455.46/tc[1]);
        k_f_old[32] = 1e-06 * 1.09e+12*exp(0.48*tc[0]+130.851/tc[1]);
        k_f_old[33] = 1e-06 * 7.34e+13;
        k_f_old[34] = 1e-06 * 5.4e+11*exp(0.454*tc[0]-1308.51/tc[1]);
        k_f_old[35] = 1e-06 * 2.3e+10*exp(1.05*tc[0]-1648.21/tc[1]);
        k_f_old[36] = 1e-06 * 3.2e+13;
        k_f_old[37] = 1e-06 * 1.08e+12*exp(0.454*tc[0]-915.954/tc[1]);
        k_f_old[38] = 1e-06 * 5.21e+17*exp(-0.99*tc[0]-795.169/tc[1]);
        k_f_old[39] = 1e-06 * 1.15e+08*exp(1.9*tc[0]-3789.63/tc[1]);
        k_f_old[40] = 1e-06 * 4.3e+07*exp(1.5*tc[0]-40060.4/tc[1]);
        k_f_old[41] = 1e-06 * 2.16e+08*exp(1.51*tc[0]-1726.22/tc[1]);
        k_f_old[42] = 1e-06 * 35700*exp(2.4*tc[0]+1061.9/tc[1]);
        k_f_old[43] = 1e-06 * 2.9e+13*exp(+251.636/tc[1]);
        k_f_old[44] = 1e-06 * 2e+13;
        k_f_old[45] = 1e-06 * 3e+13;
        k_f_old[46] = 1e-06 * 5.6e+07*exp(1.6*tc[0]-2727.73/tc[1]);
        k_f_old[47] = 1e-06 * 2.501e+13;
        k_f_old[48] = 1e-06 * 1e+08*exp(1.6*tc[0]-1570.21/tc[1]);
        k_f_old[49] = 1e-06 * 4.76e+07*exp(1.228*tc[0]-35.229/tc[1]);
        k_f_old[50] = 1e-06 * 5e+13;
        k_f_old[51] = 1e-06 * 3.43e+09*exp(1.18*tc[0]+224.962/tc[1]);
        k_f_old[52] = 1e-06 * 3.54e+06*exp(2.12*tc[0]-437.846/tc[1]);
        k_f_old[53] = 1e-06 * 2e+13;
        k_f_old[54] = 1e-06 * 1e+12;
        k_f_old[55] = 1e-06 * 2e+13;
        k_f_old[56] = 1e-06 * 1.5e+14*exp(-11877.2/tc[1]);
        k_f_old[57] = 1e-06 * 1.32e+13*exp(-754.907/tc[1]);
        k_f_old[58] = 1e-06 * 500000*exp(2*tc[0]-3638.65/tc[1]);
        k_f_old[59] = 1e-06 * 4e+13;
        k_f_old[60] = 1e-06 * 2.46e+06*exp(2*tc[0]-4162.05/tc[1]);
        k_f_old[61] = 1e-06 * 1.5e+13*exp(-301.963/tc[1]);
        k_f_old[62] = 1e-06 * 9e+12*exp(-301.963/tc[1]);
        k_f_old[63] = 1e-06 * 2.8e+13;
        k_f_old[64] = 1e-06 * 1.2e+13;
        k_f_old[65] = 1e-06 * 7e+13;
        k_f_old[66] = 1e-06 * 3e+13;
        k_f_old[67] = 1e-06 * 1.2e+13*exp(+286.865/tc[1]);
        k_f_old[68] = 1e-06 * 1.6e+13*exp(+286.865/tc[1]);
        k_f_old[69] = 1e-06 * 9e+12;
        k_f_old[70] = 1e-06 * 7e+12;
        k_f_old[71] = 1e-06 * 1.4e+13;
        k_f_old[72] = 1e-06 * 2.675e+13*exp(-14494.2/tc[1]);
        k_f_old[73] = 1e-06 * 3.6e+10*exp(-4499.25/tc[1]);
        k_f_old[74] = 1e-06 * 2.12e+16*exp(-0.97*tc[0]-312.028/tc[1]);
        k_f_old[75] = 1e-06 * 4.99e+12*exp(0.1*tc[0]-5334.68/tc[1]);
        k_f_old[76] = 1e-06 * 2.648e+13;
        k_f_old[77] = 1e-06 * 3320*exp(2.81*tc[0]-2949.17/tc[1]);
        k_f_old[78] = 1e-06 * 6.14e+06*exp(1.74*tc[0]-5259.18/tc[1]);
        k_f_old[79] = 1e-06 * 2.244e+18*exp(-1*tc[0]-8555.61/tc[1]);
        k_f_old[80] = 1e-06 * 1.87e+17*exp(-1*tc[0]-8555.61/tc[1]);
        k_f_old[81] = 1e-06 * 7.6e+12*exp(-201.309/tc[1]);
        k_f_old[82] = 1e-06 * 4.28e-13*exp(7.6*tc[0]+1776.55/tc[1]);
        k_f_old[83] = 1e-06 * 8.4e+11*exp(-1950.18/tc[1]);

        Kc_old[0] = 1.0 / (refC) * exp((g_RT[2] + g_RT[1]) - (g_RT[4]));
        Kc_old[1] = exp((g_RT[2] + g_RT[0]) - (g_RT[1] + g_RT[4]));
        Kc_old[2] = exp((g_RT[2] + g_RT[6]) - (g_RT[4] + g_RT[3]));
        Kc_old[3] = exp((g_RT[2] + g_RT[7]) - (g_RT[1] + g_RT[13]));
        Kc_old[4] = exp((g_RT[2] + g_RT[8]) - (g_RT[1] + g_RT[13]));
        Kc_old[5] = exp((g_RT[2] + g_RT[9]) - (g_RT[1] + g_RT[14]));
        Kc_old[6] = exp((g_RT[2] + g_RT[10]) - (g_RT[4] + g_RT[9]));
        Kc_old[7] = 1.0 / (refC) * exp((g_RT[2] + g_RT[11]) - (g_RT[12]));
        Kc_old[8] = exp((g_RT[2] + g_RT[13]) - (g_RT[4] + g_RT[11]));
        Kc_old[9] = exp((g_RT[2] + g_RT[13]) - (g_RT[1] + g_RT[12]));
        Kc_old[10] = exp((g_RT[2] + g_RT[14]) - (g_RT[4] + g_RT[13]));
        Kc_old[11] = exp((g_RT[2] + g_RT[16]) - (g_RT[9] + g_RT[13]));
        Kc_old[12] = exp((g_RT[2] + g_RT[17]) - (g_RT[9] + g_RT[14]));
        Kc_old[13] = exp((g_RT[2] + g_RT[18]) - (g_RT[4] + g_RT[17]));
        Kc_old[14] = exp((g_RT[3] + g_RT[11]) - (g_RT[2] + g_RT[12]));
        Kc_old[15] = exp((g_RT[3] + g_RT[14]) - (g_RT[6] + g_RT[13]));
        Kc_old[16] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));
        Kc_old[17] = 1.0 / (refC) * exp((g_RT[1] + 2 * g_RT[3]) - (g_RT[6] + g_RT[3]));
        Kc_old[18] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));
        Kc_old[19] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[19]) - (g_RT[6] + g_RT[19]));
        Kc_old[20] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[20]) - (g_RT[6] + g_RT[20]));
        Kc_old[21] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));
        Kc_old[22] = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));
        Kc_old[23] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));
        Kc_old[24] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));
        Kc_old[25] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[12]) - (g_RT[0] + g_RT[12]));
        Kc_old[26] = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));
        Kc_old[27] = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));
        Kc_old[28] = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));
        Kc_old[29] = 1.0 / (refC) * exp((g_RT[1] + g_RT[7]) - (g_RT[9]));
        Kc_old[30] = 1.0 / (refC) * exp((g_RT[1] + g_RT[9]) - (g_RT[10]));
        Kc_old[31] = exp((g_RT[1] + g_RT[10]) - (g_RT[9] + g_RT[0]));
        Kc_old[32] = 1.0 / (refC) * exp((g_RT[1] + g_RT[13]) - (g_RT[14]));
        Kc_old[33] = exp((g_RT[1] + g_RT[13]) - (g_RT[0] + g_RT[11]));
        Kc_old[34] = 1.0 / (refC) * exp((g_RT[1] + g_RT[14]) - (g_RT[15]));
        Kc_old[35] = exp((g_RT[1] + g_RT[14]) - (g_RT[13] + g_RT[0]));
        Kc_old[36] = exp((g_RT[1] + g_RT[15]) - (g_RT[4] + g_RT[9]));
        Kc_old[37] = 1.0 / (refC) * exp((g_RT[1] + g_RT[16]) - (g_RT[17]));
        Kc_old[38] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[18]));
        Kc_old[39] = exp((g_RT[1] + g_RT[18]) - (g_RT[17] + g_RT[0]));
        Kc_old[40] = 1.0 / (refC) * exp((g_RT[0] + g_RT[11]) - (g_RT[14]));
        Kc_old[41] = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));
        Kc_old[42] = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));
        Kc_old[43] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
        Kc_old[44] = exp((g_RT[4] + g_RT[7]) - (g_RT[1] + g_RT[14]));
        Kc_old[45] = exp((g_RT[4] + g_RT[8]) - (g_RT[1] + g_RT[14]));
        Kc_old[46] = exp((g_RT[4] + g_RT[9]) - (g_RT[7] + g_RT[5]));
        Kc_old[47] = exp((g_RT[4] + g_RT[9]) - (g_RT[8] + g_RT[5]));
        Kc_old[48] = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[5]));
        Kc_old[49] = exp((g_RT[4] + g_RT[11]) - (g_RT[1] + g_RT[12]));
        Kc_old[50] = exp((g_RT[4] + g_RT[13]) - (g_RT[5] + g_RT[11]));
        Kc_old[51] = exp((g_RT[4] + g_RT[14]) - (g_RT[13] + g_RT[5]));
        Kc_old[52] = exp((g_RT[4] + g_RT[18]) - (g_RT[17] + g_RT[5]));
        Kc_old[53] = exp((g_RT[6] + g_RT[7]) - (g_RT[4] + g_RT[14]));
        Kc_old[54] = exp((g_RT[6] + g_RT[9]) - (g_RT[3] + g_RT[10]));
        Kc_old[55] = exp((g_RT[6] + g_RT[9]) - (g_RT[4] + g_RT[15]));
        Kc_old[56] = exp((g_RT[6] + g_RT[11]) - (g_RT[4] + g_RT[12]));
        Kc_old[57] = exp((g_RT[7] + g_RT[3]) - (g_RT[4] + g_RT[13]));
        Kc_old[58] = exp((g_RT[7] + g_RT[0]) - (g_RT[1] + g_RT[9]));
        Kc_old[59] = exp((g_RT[7] + g_RT[9]) - (g_RT[1] + g_RT[16]));
        Kc_old[60] = exp((g_RT[7] + g_RT[10]) - (2 * g_RT[9]));
        Kc_old[61] = exp((g_RT[8] + g_RT[19]) - (g_RT[7] + g_RT[19]));
        Kc_old[62] = exp((g_RT[8] + g_RT[20]) - (g_RT[7] + g_RT[20]));
        Kc_old[63] = refC * exp((g_RT[8] + g_RT[3]) - (g_RT[1] + g_RT[4] + g_RT[11]));
        Kc_old[64] = exp((g_RT[8] + g_RT[3]) - (g_RT[11] + g_RT[5]));
        Kc_old[65] = exp((g_RT[8] + g_RT[0]) - (g_RT[9] + g_RT[1]));
        Kc_old[66] = exp((g_RT[8] + g_RT[5]) - (g_RT[7] + g_RT[5]));
        Kc_old[67] = exp((g_RT[8] + g_RT[9]) - (g_RT[1] + g_RT[16]));
        Kc_old[68] = exp((g_RT[8] + g_RT[10]) - (2 * g_RT[9]));
        Kc_old[69] = exp((g_RT[8] + g_RT[11]) - (g_RT[7] + g_RT[11]));
        Kc_old[70] = exp((g_RT[8] + g_RT[12]) - (g_RT[7] + g_RT[12]));
        Kc_old[71] = exp((g_RT[8] + g_RT[12]) - (g_RT[11] + g_RT[14]));
        Kc_old[72] = exp((g_RT[9] + g_RT[3]) - (g_RT[2] + g_RT[15]));
        Kc_old[73] = exp((g_RT[9] + g_RT[3]) - (g_RT[4] + g_RT[14]));
        Kc_old[74] = 1.0 / (refC) * exp((2 * g_RT[9]) - (g_RT[18]));
        Kc_old[75] = exp((2 * g_RT[9]) - (g_RT[1] + g_RT[17]));
        Kc_old[76] = exp((g_RT[9] + g_RT[13]) - (g_RT[10] + g_RT[11]));
        Kc_old[77] = exp((g_RT[9] + g_RT[14]) - (g_RT[13] + g_RT[10]));
        Kc_old[78] = exp((g_RT[9] + g_RT[18]) - (g_RT[17] + g_RT[10]));
        Kc_old[79] = refC * exp((g_RT[13] + g_RT[5]) - (g_RT[1] + g_RT[11] + g_RT[5]));
        Kc_old[80] = refC * exp((g_RT[13]) - (g_RT[1] + g_RT[11]));
        Kc_old[81] = exp((g_RT[13] + g_RT[3]) - (g_RT[6] + g_RT[11]));
        Kc_old[82] = exp((g_RT[15] + g_RT[3]) - (g_RT[6] + g_RT[14]));
        Kc_old[83] = exp((g_RT[17] + g_RT[3]) - (g_RT[6] + g_RT[16]));
    }

    /*reaction 1: O + H + M <=> OH + M */
    phi_f = sc[2]*sc[1];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = alpha * k_f_old[0];
    q_f = phi_f * k_f;
    phi_r = sc[4];
    Kc = Kc_old[0];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -=     qdot;
    wdot[1] -=     qdot;
    wdot[4] +=     qdot;

    /*reaction 2: O + H2 <=> H + OH */
    phi_f = sc[2]*sc[0];
    k_f = k_f_old[1];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[4];
    Kc = Kc_old[1];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -=     qdot;
    wdot[0] -=     qdot;
    wdot[1] +=     qdot;
    wdot[4] +=     qdot;

    /*reaction 3: O + HO2 <=> OH + O2 */
    phi_f = sc[2]*sc[6];
    k_f = k_f_old[2];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = Kc_old[2];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -=     qdot;
    wdot[6] -=     qdot;
    wdot[4] +=     qdot;
    wdot[3] +=     qdot;

    /*reaction 4: O + CH2 <=> H + HCO */
    phi_f = sc[2]*sc[7];
    k_f = k_f_old[3];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[13];
    Kc = Kc_old[3];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -=     qdot;
    wdot[7] -=     qdot;
    wdot[1] +=     qdot;
    wdot[13] +=     qdot;

    /*reaction 5: O + CH2(S) <=> H + HCO */
    phi_f = sc[2]*sc[8];
    k_f = k_f_old[4];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[13];
    Kc = Kc_old[4];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -=     qdot;
    wdot[8] -=     qdot;
    wdot[1] +=     qdot;
    wdot[13] +=     qdot;

    /*reaction 6: O + CH3 <=> H + CH2O */
    phi_f = sc[2]*sc[9];
    k_f = k_f_old[5];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_old[5];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -=     qdot;
    wdot[9] -=     qdot;
    wdot[1] +=     qdot;
    wdot[14] +=     qdot;

    /*reaction 7: O + CH4 <=> OH + CH3 */
    phi_f = sc[2]*sc[10];
    k_f = k_f_old[6];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = Kc_old[6];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -=     qdot;
    wdot[10] -=     qdot;
    wdot[4] +=     qdot;
    wdot[9] +=     qdot;

    /*reaction 8: O + CO + M <=> CO2 + M */
    phi_f = sc[2]*sc[11];
    alpha = mixture + sc[0] + 5*sc[3] + 5*sc[5] + sc[10] + 0.5*sc[11] + 2.5*sc[12] + 2*sc[18] + -0.5*sc[20];
    k_f = alpha * k_f_old[7];
    q_f = phi_f * k_f;
    phi_r = sc[12];
    Kc = Kc_old[7];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -=     qdot;
    wdot[11] -=     qdot;
    wdot[12] +=     qdot;

    /*reaction 9: O + HCO <=> OH + CO */
    phi_f = sc[2]*sc[13];
    k_f = k_f_old[8];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[11];
    Kc = Kc_old[8];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -=     qdot;
    wdot[13] -=     qdot;
    wdot[4] +=     qdot;
    wdot[11] +=     qdot;

    /*reaction 10: O + HCO <=> H + CO2 */
    phi_f = sc[2]*sc[13];
    k_f = k_f_old[9];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = Kc_old[9];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -=     qdot;
    wdot[13] -=     qdot;
    wdot[1] +=     qdot;
    wdot[12] +=     qdot;

    /*reaction 11: O + CH2O <=> OH + HCO */
    phi_f = sc[2]*sc[14];
    k_f = k_f_old[10];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[13];
    Kc = Kc_old[10];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -=     qdot;
    wdot[14] -=     qdot;
    wdot[4] +=     qdot;
    wdot[13] +=     qdot;

    /*reaction 12: O + C2H4 <=> CH3 + HCO */
    phi_f = sc[2]*sc[16];
    k_f = k_f_old[11];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[13];
    Kc = Kc_old[11];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -=     qdot;
    wdot[16] -=     qdot;
    wdot[9] +=     qdot;
    wdot[13] +=     qdot;

    /*reaction 13: O + C2H5 <=> CH3 + CH2O */
    phi_f = sc[2]*sc[17];
    k_f = k_f_old[12];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[14];
    Kc = Kc_old[12];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -=     qdot;
    wdot[17] -=     qdot;
    wdot[9] +=     qdot;
    wdot[14] +=     qdot;

    /*reaction 14: O + C2H6 <=> OH + C2H5 */
    phi_f = sc[2]*sc[18];
    k_f = k_f_old[13];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = Kc_old[13];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -=     qdot;
    wdot[18] -=     qdot;
    wdot[4] +=     qdot;
    wdot[17] +=     qdot;

    /*reaction 15: O2 + CO <=> O + CO2 */
    phi_f = sc[3]*sc[11];
    k_f = k_f_old[14];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[12];
    Kc = Kc_old[14];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -=     qdot;
    wdot[11] -=     qdot;
    wdot[2] +=     qdot;
    wdot[12] +=     qdot;

    /*reaction 16: O2 + CH2O <=> HO2 + HCO */
    phi_f = sc[3]*sc[14];
    k_f = k_f_old[15];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = Kc_old[15];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -=     qdot;
    wdot[14] -=     qdot;
    wdot[6] +=     qdot;
    wdot[13] +=     qdot;

    /*reaction 17: H + O2 + M <=> HO2 + M */
    phi_f = sc[1]*sc[3];
    alpha = mixture + -1*sc[3] + -1*sc[5] + -0.25*sc[11] + 0.5*sc[12] + 0.5*sc[18] + -1*sc[19] + -1*sc[20];
    k_f = alpha * k_f_old[16];
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = Kc_old[16];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[3] -=     qdot;
    wdot[6] +=     qdot;

    /*reaction 18: H + 2 O2 <=> HO2 + O2 */
    phi_f = sc[1]*sc[3]*sc[3];
    k_f = k_f_old[17];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = Kc_old[17];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[3] -= 2 * qdot;
    wdot[6] +=     qdot;
    wdot[3] +=     qdot;

    /*reaction 19: H + O2 + H2O <=> HO2 + H2O */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = k_f_old[18];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = Kc_old[18];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[3] -=     qdot;
    wdot[5] -=     qdot;
    wdot[6] +=     qdot;
    wdot[5] +=     qdot;

    /*reaction 20: H + O2 + N2 <=> HO2 + N2 */
    phi_f = sc[1]*sc[3]*sc[19];
    k_f = k_f_old[19];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[19];
    Kc = Kc_old[19];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[3] -=     qdot;
    wdot[19] -=     qdot;
    wdot[6] +=     qdot;
    wdot[19] +=     qdot;

    /*reaction 21: H + O2 + AR <=> HO2 + AR */
    phi_f = sc[1]*sc[3]*sc[20];
    k_f = k_f_old[20];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[20];
    Kc = Kc_old[20];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[3] -=     qdot;
    wdot[20] -=     qdot;
    wdot[6] +=     qdot;
    wdot[20] +=     qdot;

    /*reaction 22: H + O2 <=> O + OH */
    phi_f = sc[1]*sc[3];
    k_f = k_f_old[21];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = Kc_old[21];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[3] -=     qdot;
    wdot[2] +=     qdot;
    wdot[4] +=     qdot;

    /*reaction 23: 2 H + M <=> H2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[0] + -1*sc[5] + sc[10] + -1*sc[12] + 2*sc[18] + -0.37*sc[20];
    k_f = alpha * k_f_old[22];
    q_f = phi_f * k_f;
    phi_r = sc[0];
    Kc = Kc_old[22];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[0] +=     qdot;

    /*reaction 24: 2 H + H2 <=> 2 H2 */
    phi_f = sc[1]*sc[1]*sc[0];
    k_f = k_f_old[23];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[0];
    Kc = Kc_old[23];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[0] -=     qdot;
    wdot[0] += 2 * qdot;

    /*reaction 25: 2 H + H2O <=> H2 + H2O */
    phi_f = sc[1]*sc[1]*sc[5];
    k_f = k_f_old[24];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[5];
    Kc = Kc_old[24];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[5] -=     qdot;
    wdot[0] +=     qdot;
    wdot[5] +=     qdot;

    /*reaction 26: 2 H + CO2 <=> H2 + CO2 */
    phi_f = sc[1]*sc[1]*sc[12];
    k_f = k_f_old[25];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[12];
    Kc = Kc_old[25];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[12] -=     qdot;
    wdot[0] +=     qdot;
    wdot[12] +=     qdot;

    /*reaction 27: H + OH + M <=> H2O + M */
    phi_f = sc[1]*sc[4];
    alpha = mixture + -0.27*sc[0] + 2.65*sc[5] + sc[10] + 2*sc[18] + -0.62*sc[20];
    k_f = alpha * k_f_old[26];
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = Kc_old[26];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[4] -=     qdot;
    wdot[5] +=     qdot;

    /*reaction 28: H + HO2 <=> O2 + H2 */
    phi_f = sc[1]*sc[6];
    k_f = k_f_old[27];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = Kc_old[27];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[6] -=     qdot;
    wdot[3] +=     qdot;
    wdot[0] +=     qdot;

    /*reaction 29: H + HO2 <=> 2 OH */
    phi_f = sc[1]*sc[6];
    k_f = k_f_old[28];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = Kc_old[28];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[6] -=     qdot;
    wdot[4] += 2 * qdot;

    /*reaction 30: H + CH2 (+M) <=> CH3 (+M) */
    phi_f = sc[1]*sc[7];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_old[29];
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
    phi_r = sc[9];
    Kc = Kc_old[29];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[7] -=     qdot;
    wdot[9] +=     qdot;

    /*reaction 31: H + CH3 (+M) <=> CH4 (+M) */
    phi_f = sc[1]*sc[9];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_old[30];
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
    phi_r = sc[10];
    Kc = Kc_old[30];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[9] -=     qdot;
    wdot[10] +=     qdot;

    /*reaction 32: H + CH4 <=> CH3 + H2 */
    phi_f = sc[1]*sc[10];
    k_f = k_f_old[31];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = Kc_old[31];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[10] -=     qdot;
    wdot[9] +=     qdot;
    wdot[0] +=     qdot;

    /*reaction 33: H + HCO (+M) <=> CH2O (+M) */
    phi_f = sc[1]*sc[13];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_old[32];
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
    phi_r = sc[14];
    Kc = Kc_old[32];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[13] -=     qdot;
    wdot[14] +=     qdot;

    /*reaction 34: H + HCO <=> H2 + CO */
    phi_f = sc[1]*sc[13];
    k_f = k_f_old[33];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[11];
    Kc = Kc_old[33];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[13] -=     qdot;
    wdot[0] +=     qdot;
    wdot[11] +=     qdot;

    /*reaction 35: H + CH2O (+M) <=> CH3O (+M) */
    phi_f = sc[1]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18];
    k_f = k_f_old[34];
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
    phi_r = sc[15];
    Kc = Kc_old[34];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[14] -=     qdot;
    wdot[15] +=     qdot;

    /*reaction 36: H + CH2O <=> HCO + H2 */
    phi_f = sc[1]*sc[14];
    k_f = k_f_old[35];
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[0];
    Kc = Kc_old[35];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[14] -=     qdot;
    wdot[13] +=     qdot;
    wdot[0] +=     qdot;

    /*reaction 37: H + CH3O <=> OH + CH3 */
    phi_f = sc[1]*sc[15];
    k_f = k_f_old[36];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = Kc_old[36];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[15] -=     qdot;
    wdot[4] +=     qdot;
    wdot[9] +=     qdot;

    /*reaction 38: H + C2H4 (+M) <=> C2H5 (+M) */
    phi_f = sc[1]*sc[16];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_old[37];
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
    phi_r = sc[17];
    Kc = Kc_old[37];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[16] -=     qdot;
    wdot[17] +=     qdot;

    /*reaction 39: H + C2H5 (+M) <=> C2H6 (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_old[38];
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
    phi_r = sc[18];
    Kc = Kc_old[38];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[17] -=     qdot;
    wdot[18] +=     qdot;

    /*reaction 40: H + C2H6 <=> C2H5 + H2 */
    phi_f = sc[1]*sc[18];
    k_f = k_f_old[39];
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[0];
    Kc = Kc_old[39];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -=     qdot;
    wdot[18] -=     qdot;
    wdot[17] +=     qdot;
    wdot[0] +=     qdot;

    /*reaction 41: H2 + CO (+M) <=> CH2O (+M) */
    phi_f = sc[0]*sc[11];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_old[40];
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
    phi_r = sc[14];
    Kc = Kc_old[40];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -=     qdot;
    wdot[11] -=     qdot;
    wdot[14] +=     qdot;

    /*reaction 42: OH + H2 <=> H + H2O */
    phi_f = sc[4]*sc[0];
    k_f = k_f_old[41];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    Kc = Kc_old[41];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -=     qdot;
    wdot[0] -=     qdot;
    wdot[1] +=     qdot;
    wdot[5] +=     qdot;

    /*reaction 43: 2 OH <=> O + H2O */
    phi_f = sc[4]*sc[4];
    k_f = k_f_old[42];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = Kc_old[42];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 2 * qdot;
    wdot[2] +=     qdot;
    wdot[5] +=     qdot;

    /*reaction 44: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = k_f_old[43];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = Kc_old[43];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -=     qdot;
    wdot[6] -=     qdot;
    wdot[3] +=     qdot;
    wdot[5] +=     qdot;

    /*reaction 45: OH + CH2 <=> H + CH2O */
    phi_f = sc[4]*sc[7];
    k_f = k_f_old[44];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_old[44];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -=     qdot;
    wdot[7] -=     qdot;
    wdot[1] +=     qdot;
    wdot[14] +=     qdot;

    /*reaction 46: OH + CH2(S) <=> H + CH2O */
    phi_f = sc[4]*sc[8];
    k_f = k_f_old[45];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = Kc_old[45];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -=     qdot;
    wdot[8] -=     qdot;
    wdot[1] +=     qdot;
    wdot[14] +=     qdot;

    /*reaction 47: OH + CH3 <=> CH2 + H2O */
    phi_f = sc[4]*sc[9];
    k_f = k_f_old[46];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[5];
    Kc = Kc_old[46];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -=     qdot;
    wdot[9] -=     qdot;
    wdot[7] +=     qdot;
    wdot[5] +=     qdot;

    /*reaction 48: OH + CH3 <=> CH2(S) + H2O */
    phi_f = sc[4]*sc[9];
    k_f = k_f_old[47];
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[5];
    Kc = Kc_old[47];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -=     qdot;
    wdot[9] -=     qdot;
    wdot[8] +=     qdot;
    wdot[5] +=     qdot;

    /*reaction 49: OH + CH4 <=> CH3 + H2O */
    phi_f = sc[4]*sc[10];
    k_f = k_f_old[48];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[5];
    Kc = Kc_old[48];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -=     qdot;
    wdot[10] -=     qdot;
    wdot[9] +=     qdot;
    wdot[5] +=     qdot;

    /*reaction 50: OH + CO <=> H + CO2 */
    phi_f = sc[4]*sc[11];
    k_f = k_f_old[49];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = Kc_old[49];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -=     qdot;
    wdot[11] -=     qdot;
    wdot[1] +=     qdot;
    wdot[12] +=     qdot;

    /*reaction 51: OH + HCO <=> H2O + CO */
    phi_f = sc[4]*sc[13];
    k_f = k_f_old[50];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[11];
    Kc = Kc_old[50];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -=     qdot;
    wdot[13] -=     qdot;
    wdot[5] +=     qdot;
    wdot[11] +=     qdot;

    /*reaction 52: OH + CH2O <=> HCO + H2O */
    phi_f = sc[4]*sc[14];
    k_f = k_f_old[51];
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[5];
    Kc = Kc_old[51];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -=     qdot;
    wdot[14] -=     qdot;
    wdot[13] +=     qdot;
    wdot[5] +=     qdot;

    /*reaction 53: OH + C2H6 <=> C2H5 + H2O */
    phi_f = sc[4]*sc[18];
    k_f = k_f_old[52];
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[5];
    Kc = Kc_old[52];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -=     qdot;
    wdot[18] -=     qdot;
    wdot[17] +=     qdot;
    wdot[5] +=     qdot;

    /*reaction 54: HO2 + CH2 <=> OH + CH2O */
    phi_f = sc[6]*sc[7];
    k_f = k_f_old[53];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[14];
    Kc = Kc_old[53];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -=     qdot;
    wdot[7] -=     qdot;
    wdot[4] +=     qdot;
    wdot[14] +=     qdot;

    /*reaction 55: HO2 + CH3 <=> O2 + CH4 */
    phi_f = sc[6]*sc[9];
    k_f = k_f_old[54];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[10];
    Kc = Kc_old[54];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -=     qdot;
    wdot[9] -=     qdot;
    wdot[3] +=     qdot;
    wdot[10] +=     qdot;

    /*reaction 56: HO2 + CH3 <=> OH + CH3O */
    phi_f = sc[6]*sc[9];
    k_f = k_f_old[55];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[15];
    Kc = Kc_old[55];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -=     qdot;
    wdot[9] -=     qdot;
    wdot[4] +=     qdot;
    wdot[15] +=     qdot;

    /*reaction 57: HO2 + CO <=> OH + CO2 */
    phi_f = sc[6]*sc[11];
    k_f = k_f_old[56];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = Kc_old[56];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -=     qdot;
    wdot[11] -=     qdot;
    wdot[4] +=     qdot;
    wdot[12] +=     qdot;

    /*reaction 58: CH2 + O2 <=> OH + HCO */
    phi_f = sc[7]*sc[3];
    k_f = k_f_old[57];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[13];
    Kc = Kc_old[57];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -=     qdot;
    wdot[3] -=     qdot;
    wdot[4] +=     qdot;
    wdot[13] +=     qdot;

    /*reaction 59: CH2 + H2 <=> H + CH3 */
    phi_f = sc[7]*sc[0];
    k_f = k_f_old[58];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[9];
    Kc = Kc_old[58];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -=     qdot;
    wdot[0] -=     qdot;
    wdot[1] +=     qdot;
    wdot[9] +=     qdot;

    /*reaction 60: CH2 + CH3 <=> H + C2H4 */
    phi_f = sc[7]*sc[9];
    k_f = k_f_old[59];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_old[59];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -=     qdot;
    wdot[9] -=     qdot;
    wdot[1] +=     qdot;
    wdot[16] +=     qdot;

    /*reaction 61: CH2 + CH4 <=> 2 CH3 */
    phi_f = sc[7]*sc[10];
    k_f = k_f_old[60];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[9];
    Kc = Kc_old[60];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -=     qdot;
    wdot[10] -=     qdot;
    wdot[9] += 2 * qdot;

    /*reaction 62: CH2(S) + N2 <=> CH2 + N2 */
    phi_f = sc[8]*sc[19];
    k_f = k_f_old[61];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[19];
    Kc = Kc_old[61];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -=     qdot;
    wdot[19] -=     qdot;
    wdot[7] +=     qdot;
    wdot[19] +=     qdot;

    /*reaction 63: CH2(S) + AR <=> CH2 + AR */
    phi_f = sc[8]*sc[20];
    k_f = k_f_old[62];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[20];
    Kc = Kc_old[62];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -=     qdot;
    wdot[20] -=     qdot;
    wdot[7] +=     qdot;
    wdot[20] +=     qdot;

    /*reaction 64: CH2(S) + O2 <=> H + OH + CO */
    phi_f = sc[8]*sc[3];
    k_f = k_f_old[63];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[4]*sc[11];
    Kc = Kc_old[63];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -=     qdot;
    wdot[3] -=     qdot;
    wdot[1] +=     qdot;
    wdot[4] +=     qdot;
    wdot[11] +=     qdot;

    /*reaction 65: CH2(S) + O2 <=> CO + H2O */
    phi_f = sc[8]*sc[3];
    k_f = k_f_old[64];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = Kc_old[64];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -=     qdot;
    wdot[3] -=     qdot;
    wdot[11] +=     qdot;
    wdot[5] +=     qdot;

    /*reaction 66: CH2(S) + H2 <=> CH3 + H */
    phi_f = sc[8]*sc[0];
    k_f = k_f_old[65];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[1];
    Kc = Kc_old[65];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -=     qdot;
    wdot[0] -=     qdot;
    wdot[9] +=     qdot;
    wdot[1] +=     qdot;

    /*reaction 67: CH2(S) + H2O <=> CH2 + H2O */
    phi_f = sc[8]*sc[5];
    k_f = k_f_old[66];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[5];
    Kc = Kc_old[66];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -=     qdot;
    wdot[5] -=     qdot;
    wdot[7] +=     qdot;
    wdot[5] +=     qdot;

    /*reaction 68: CH2(S) + CH3 <=> H + C2H4 */
    phi_f = sc[8]*sc[9];
    k_f = k_f_old[67];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = Kc_old[67];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -=     qdot;
    wdot[9] -=     qdot;
    wdot[1] +=     qdot;
    wdot[16] +=     qdot;

    /*reaction 69: CH2(S) + CH4 <=> 2 CH3 */
    phi_f = sc[8]*sc[10];
    k_f = k_f_old[68];
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[9];
    Kc = Kc_old[68];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -=     qdot;
    wdot[10] -=     qdot;
    wdot[9] += 2 * qdot;

    /*reaction 70: CH2(S) + CO <=> CH2 + CO */
    phi_f = sc[8]*sc[11];
    k_f = k_f_old[69];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[11];
    Kc = Kc_old[69];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -=     qdot;
    wdot[11] -=     qdot;
    wdot[7] +=     qdot;
    wdot[11] +=     qdot;

    /*reaction 71: CH2(S) + CO2 <=> CH2 + CO2 */
    phi_f = sc[8]*sc[12];
    k_f = k_f_old[70];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[12];
    Kc = Kc_old[70];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -=     qdot;
    wdot[12] -=     qdot;
    wdot[7] +=     qdot;
    wdot[12] +=     qdot;

    /*reaction 72: CH2(S) + CO2 <=> CO + CH2O */
    phi_f = sc[8]*sc[12];
    k_f = k_f_old[71];
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[14];
    Kc = Kc_old[71];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -=     qdot;
    wdot[12] -=     qdot;
    wdot[11] +=     qdot;
    wdot[14] +=     qdot;

    /*reaction 73: CH3 + O2 <=> O + CH3O */
    phi_f = sc[9]*sc[3];
    k_f = k_f_old[72];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[15];
    Kc = Kc_old[72];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -=     qdot;
    wdot[3] -=     qdot;
    wdot[2] +=     qdot;
    wdot[15] +=     qdot;

    /*reaction 74: CH3 + O2 <=> OH + CH2O */
    phi_f = sc[9]*sc[3];
    k_f = k_f_old[73];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[14];
    Kc = Kc_old[73];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -=     qdot;
    wdot[3] -=     qdot;
    wdot[4] +=     qdot;
    wdot[14] +=     qdot;

    /*reaction 75: 2 CH3 (+M) <=> C2H6 (+M) */
    phi_f = sc[9]*sc[9];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = k_f_old[74];
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
    phi_r = sc[18];
    Kc = Kc_old[74];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 2 * qdot;
    wdot[18] +=     qdot;

    /*reaction 76: 2 CH3 <=> H + C2H5 */
    phi_f = sc[9]*sc[9];
    k_f = k_f_old[75];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = Kc_old[75];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 2 * qdot;
    wdot[1] +=     qdot;
    wdot[17] +=     qdot;

    /*reaction 77: CH3 + HCO <=> CH4 + CO */
    phi_f = sc[9]*sc[13];
    k_f = k_f_old[76];
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[11];
    Kc = Kc_old[76];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -=     qdot;
    wdot[13] -=     qdot;
    wdot[10] +=     qdot;
    wdot[11] +=     qdot;

    /*reaction 78: CH3 + CH2O <=> HCO + CH4 */
    phi_f = sc[9]*sc[14];
    k_f = k_f_old[77];
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[10];
    Kc = Kc_old[77];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -=     qdot;
    wdot[14] -=     qdot;
    wdot[13] +=     qdot;
    wdot[10] +=     qdot;

    /*reaction 79: CH3 + C2H6 <=> C2H5 + CH4 */
    phi_f = sc[9]*sc[18];
    k_f = k_f_old[78];
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[10];
    Kc = Kc_old[78];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -=     qdot;
    wdot[18] -=     qdot;
    wdot[17] +=     qdot;
    wdot[10] +=     qdot;

    /*reaction 80: HCO + H2O <=> H + CO + H2O */
    phi_f = sc[13]*sc[5];
    k_f = k_f_old[79];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[11]*sc[5];
    Kc = Kc_old[79];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -=     qdot;
    wdot[5] -=     qdot;
    wdot[1] +=     qdot;
    wdot[11] +=     qdot;
    wdot[5] +=     qdot;

    /*reaction 81: HCO + M <=> H + CO + M */
    phi_f = sc[13];
    alpha = mixture + sc[0] + -1*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18];
    k_f = k_f_old[80];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[11];
    Kc = Kc_old[80];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -=     qdot;
    wdot[1] +=     qdot;
    wdot[11] +=     qdot;

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    phi_f = sc[13]*sc[3];
    k_f = k_f_old[81];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[11];
    Kc = Kc_old[81];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -=     qdot;
    wdot[3] -=     qdot;
    wdot[6] +=     qdot;
    wdot[11] +=     qdot;

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    phi_f = sc[15]*sc[3];
    k_f = k_f_old[82];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[14];
    Kc = Kc_old[82];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -=     qdot;
    wdot[3] -=     qdot;
    wdot[6] +=     qdot;
    wdot[14] +=     qdot;

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    phi_f = sc[17]*sc[3];
    k_f = k_f_old[83];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[16];
    Kc = Kc_old[83];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[17] -=     qdot;
    wdot[3] -=     qdot;
    wdot[6] +=     qdot;
    wdot[16] +=     qdot;

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

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 21; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*reaction 1: O + H + M <=> OH + M */
    phi_f = sc[2]*sc[1];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = 1e-12 * alpha * 5e+17*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[4];
    Kc = 1.0 / (refC) * exp((g_RT[2] + g_RT[1]) - (g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[0] = q_f - q_r;

    /*reaction 2: O + H2 <=> H + OH */
    phi_f = sc[2]*sc[0];
    k_f = 1e-06 * 50000*exp(2.67*tc[0]-3165.58/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[4];
    Kc = exp((g_RT[2] + g_RT[0]) - (g_RT[1] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[1] = q_f - q_r;

    /*reaction 3: O + HO2 <=> OH + O2 */
    phi_f = sc[2]*sc[6];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = exp((g_RT[2] + g_RT[6]) - (g_RT[4] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[2] = q_f - q_r;

    /*reaction 4: O + CH2 <=> H + HCO */
    phi_f = sc[2]*sc[7];
    k_f = 1e-06 * 8e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[13];
    Kc = exp((g_RT[2] + g_RT[7]) - (g_RT[1] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[3] = q_f - q_r;

    /*reaction 5: O + CH2(S) <=> H + HCO */
    phi_f = sc[2]*sc[8];
    k_f = 1e-06 * 1.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[13];
    Kc = exp((g_RT[2] + g_RT[8]) - (g_RT[1] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[4] = q_f - q_r;

    /*reaction 6: O + CH3 <=> H + CH2O */
    phi_f = sc[2]*sc[9];
    k_f = 1e-06 * 8.43e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = exp((g_RT[2] + g_RT[9]) - (g_RT[1] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[5] = q_f - q_r;

    /*reaction 7: O + CH4 <=> OH + CH3 */
    phi_f = sc[2]*sc[10];
    k_f = 1e-06 * 1.02e+09*exp(1.5*tc[0]-4328.13/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = exp((g_RT[2] + g_RT[10]) - (g_RT[4] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[6] = q_f - q_r;

    /*reaction 8: O + CO + M <=> CO2 + M */
    phi_f = sc[2]*sc[11];
    alpha = mixture + sc[0] + 5*sc[3] + 5*sc[5] + sc[10] + 0.5*sc[11] + 2.5*sc[12] + 2*sc[18] + -0.5*sc[20];
    k_f = 1e-12 * alpha * 6.02e+14*exp(-1509.81/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[12];
    Kc = 1.0 / (refC) * exp((g_RT[2] + g_RT[11]) - (g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[7] = q_f - q_r;

    /*reaction 9: O + HCO <=> OH + CO */
    phi_f = sc[2]*sc[13];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[11];
    Kc = exp((g_RT[2] + g_RT[13]) - (g_RT[4] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[8] = q_f - q_r;

    /*reaction 10: O + HCO <=> H + CO2 */
    phi_f = sc[2]*sc[13];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = exp((g_RT[2] + g_RT[13]) - (g_RT[1] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[9] = q_f - q_r;

    /*reaction 11: O + CH2O <=> OH + HCO */
    phi_f = sc[2]*sc[14];
    k_f = 1e-06 * 3.9e+13*exp(-1781.58/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[13];
    Kc = exp((g_RT[2] + g_RT[14]) - (g_RT[4] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[10] = q_f - q_r;

    /*reaction 12: O + C2H4 <=> CH3 + HCO */
    phi_f = sc[2]*sc[16];
    k_f = 1e-06 * 1.92e+07*exp(1.83*tc[0]-110.72/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[13];
    Kc = exp((g_RT[2] + g_RT[16]) - (g_RT[9] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[11] = q_f - q_r;

    /*reaction 13: O + C2H5 <=> CH3 + CH2O */
    phi_f = sc[2]*sc[17];
    k_f = 1e-06 * 1.32e+14;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[14];
    Kc = exp((g_RT[2] + g_RT[17]) - (g_RT[9] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[12] = q_f - q_r;

    /*reaction 14: O + C2H6 <=> OH + C2H5 */
    phi_f = sc[2]*sc[18];
    k_f = 1e-06 * 8.98e+07*exp(1.92*tc[0]-2863.61/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = exp((g_RT[2] + g_RT[18]) - (g_RT[4] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[13] = q_f - q_r;

    /*reaction 15: O2 + CO <=> O + CO2 */
    phi_f = sc[3]*sc[11];
    k_f = 1e-06 * 2.5e+12*exp(-24056.4/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[12];
    Kc = exp((g_RT[3] + g_RT[11]) - (g_RT[2] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[14] = q_f - q_r;

    /*reaction 16: O2 + CH2O <=> HO2 + HCO */
    phi_f = sc[3]*sc[14];
    k_f = 1e-06 * 1e+14*exp(-20130.9/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = exp((g_RT[3] + g_RT[14]) - (g_RT[6] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[15] = q_f - q_r;

    /*reaction 17: H + O2 + M <=> HO2 + M */
    phi_f = sc[1]*sc[3];
    alpha = mixture + -1*sc[3] + -1*sc[5] + -0.25*sc[11] + 0.5*sc[12] + 0.5*sc[18] + -1*sc[19] + -1*sc[20];
    k_f = 1e-12 * alpha * 2.8e+18*exp(-0.86*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[16] = q_f - q_r;

    /*reaction 18: H + 2 O2 <=> HO2 + O2 */
    phi_f = sc[1]*sc[3]*sc[3];
    k_f = 1e-12 * 3e+20*exp(-1.72*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = 1.0 / (refC) * exp((g_RT[1] + 2 * g_RT[3]) - (g_RT[6] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[17] = q_f - q_r;

    /*reaction 19: H + O2 + H2O <=> HO2 + H2O */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = 1e-12 * 9.38e+18*exp(-0.76*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[18] = q_f - q_r;

    /*reaction 20: H + O2 + N2 <=> HO2 + N2 */
    phi_f = sc[1]*sc[3]*sc[19];
    k_f = 1e-12 * 3.75e+20*exp(-1.72*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[19];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[19]) - (g_RT[6] + g_RT[19]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[19] = q_f - q_r;

    /*reaction 21: H + O2 + AR <=> HO2 + AR */
    phi_f = sc[1]*sc[3]*sc[20];
    k_f = 1e-12 * 7e+17*exp(-0.8*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[20];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[20]) - (g_RT[6] + g_RT[20]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[20] = q_f - q_r;

    /*reaction 22: H + O2 <=> O + OH */
    phi_f = sc[1]*sc[3];
    k_f = 1e-06 * 8.3e+13*exp(-7253.65/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[21] = q_f - q_r;

    /*reaction 23: 2 H + M <=> H2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[0] + -1*sc[5] + sc[10] + -1*sc[12] + 2*sc[18] + -0.37*sc[20];
    k_f = 1e-12 * alpha * 1e+18*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[22] = q_f - q_r;

    /*reaction 24: 2 H + H2 <=> 2 H2 */
    phi_f = sc[1]*sc[1]*sc[0];
    k_f = 1e-12 * 9e+16*exp(-0.6*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[0];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[23] = q_f - q_r;

    /*reaction 25: 2 H + H2O <=> H2 + H2O */
    phi_f = sc[1]*sc[1]*sc[5];
    k_f = 1e-12 * 6e+19*exp(-1.25*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[5];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[24] = q_f - q_r;

    /*reaction 26: 2 H + CO2 <=> H2 + CO2 */
    phi_f = sc[1]*sc[1]*sc[12];
    k_f = 1e-12 * 5.5e+20*exp(-2*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[12];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[12]) - (g_RT[0] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[25] = q_f - q_r;

    /*reaction 27: H + OH + M <=> H2O + M */
    phi_f = sc[1]*sc[4];
    alpha = mixture + -0.27*sc[0] + 2.65*sc[5] + sc[10] + 2*sc[18] + -0.62*sc[20];
    k_f = 1e-12 * alpha * 2.2e+22*exp(-2*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[26] = q_f - q_r;

    /*reaction 28: H + HO2 <=> O2 + H2 */
    phi_f = sc[1]*sc[6];
    k_f = 1e-06 * 2.8e+13*exp(-537.494/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[27] = q_f - q_r;

    /*reaction 29: H + HO2 <=> 2 OH */
    phi_f = sc[1]*sc[6];
    k_f = 1e-06 * 1.34e+14*exp(-319.577/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[28] = q_f - q_r;

    /*reaction 30: H + CH2 (+M) <=> CH3 (+M) */
    phi_f = sc[1]*sc[7];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
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
    phi_r = sc[9];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[7]) - (g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[29] = q_f - q_r;

    /*reaction 31: H + CH3 (+M) <=> CH4 (+M) */
    phi_f = sc[1]*sc[9];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
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
    phi_r = sc[10];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[9]) - (g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[30] = q_f - q_r;

    /*reaction 32: H + CH4 <=> CH3 + H2 */
    phi_f = sc[1]*sc[10];
    k_f = 1e-06 * 6.6e+08*exp(1.62*tc[0]-5455.46/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = exp((g_RT[1] + g_RT[10]) - (g_RT[9] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[31] = q_f - q_r;

    /*reaction 33: H + HCO (+M) <=> CH2O (+M) */
    phi_f = sc[1]*sc[13];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
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
    phi_r = sc[14];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[13]) - (g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[32] = q_f - q_r;

    /*reaction 34: H + HCO <=> H2 + CO */
    phi_f = sc[1]*sc[13];
    k_f = 1e-06 * 7.34e+13;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[11];
    Kc = exp((g_RT[1] + g_RT[13]) - (g_RT[0] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[33] = q_f - q_r;

    /*reaction 35: H + CH2O (+M) <=> CH3O (+M) */
    phi_f = sc[1]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18];
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
    phi_r = sc[15];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[14]) - (g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[34] = q_f - q_r;

    /*reaction 36: H + CH2O <=> HCO + H2 */
    phi_f = sc[1]*sc[14];
    k_f = 1e-06 * 2.3e+10*exp(1.05*tc[0]-1648.21/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[0];
    Kc = exp((g_RT[1] + g_RT[14]) - (g_RT[13] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[35] = q_f - q_r;

    /*reaction 37: H + CH3O <=> OH + CH3 */
    phi_f = sc[1]*sc[15];
    k_f = 1e-06 * 3.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = exp((g_RT[1] + g_RT[15]) - (g_RT[4] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[36] = q_f - q_r;

    /*reaction 38: H + C2H4 (+M) <=> C2H5 (+M) */
    phi_f = sc[1]*sc[16];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
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
    phi_r = sc[17];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[16]) - (g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[37] = q_f - q_r;

    /*reaction 39: H + C2H5 (+M) <=> C2H6 (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
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
    phi_r = sc[18];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[18]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[38] = q_f - q_r;

    /*reaction 40: H + C2H6 <=> C2H5 + H2 */
    phi_f = sc[1]*sc[18];
    k_f = 1e-06 * 1.15e+08*exp(1.9*tc[0]-3789.63/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[0];
    Kc = exp((g_RT[1] + g_RT[18]) - (g_RT[17] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[39] = q_f - q_r;

    /*reaction 41: H2 + CO (+M) <=> CH2O (+M) */
    phi_f = sc[0]*sc[11];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
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
    phi_r = sc[14];
    Kc = 1.0 / (refC) * exp((g_RT[0] + g_RT[11]) - (g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[40] = q_f - q_r;

    /*reaction 42: OH + H2 <=> H + H2O */
    phi_f = sc[4]*sc[0];
    k_f = 1e-06 * 2.16e+08*exp(1.51*tc[0]-1726.22/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    Kc = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[41] = q_f - q_r;

    /*reaction 43: 2 OH <=> O + H2O */
    phi_f = sc[4]*sc[4];
    k_f = 1e-06 * 35700*exp(2.4*tc[0]+1061.9/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[42] = q_f - q_r;

    /*reaction 44: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = 1e-06 * 2.9e+13*exp(+251.636/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[43] = q_f - q_r;

    /*reaction 45: OH + CH2 <=> H + CH2O */
    phi_f = sc[4]*sc[7];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = exp((g_RT[4] + g_RT[7]) - (g_RT[1] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[44] = q_f - q_r;

    /*reaction 46: OH + CH2(S) <=> H + CH2O */
    phi_f = sc[4]*sc[8];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = exp((g_RT[4] + g_RT[8]) - (g_RT[1] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[45] = q_f - q_r;

    /*reaction 47: OH + CH3 <=> CH2 + H2O */
    phi_f = sc[4]*sc[9];
    k_f = 1e-06 * 5.6e+07*exp(1.6*tc[0]-2727.73/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[5];
    Kc = exp((g_RT[4] + g_RT[9]) - (g_RT[7] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[46] = q_f - q_r;

    /*reaction 48: OH + CH3 <=> CH2(S) + H2O */
    phi_f = sc[4]*sc[9];
    k_f = 1e-06 * 2.501e+13;
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[5];
    Kc = exp((g_RT[4] + g_RT[9]) - (g_RT[8] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[47] = q_f - q_r;

    /*reaction 49: OH + CH4 <=> CH3 + H2O */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 1e+08*exp(1.6*tc[0]-1570.21/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[5];
    Kc = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[48] = q_f - q_r;

    /*reaction 50: OH + CO <=> H + CO2 */
    phi_f = sc[4]*sc[11];
    k_f = 1e-06 * 4.76e+07*exp(1.228*tc[0]-35.229/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = exp((g_RT[4] + g_RT[11]) - (g_RT[1] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[49] = q_f - q_r;

    /*reaction 51: OH + HCO <=> H2O + CO */
    phi_f = sc[4]*sc[13];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[11];
    Kc = exp((g_RT[4] + g_RT[13]) - (g_RT[5] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[50] = q_f - q_r;

    /*reaction 52: OH + CH2O <=> HCO + H2O */
    phi_f = sc[4]*sc[14];
    k_f = 1e-06 * 3.43e+09*exp(1.18*tc[0]+224.962/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[5];
    Kc = exp((g_RT[4] + g_RT[14]) - (g_RT[13] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[51] = q_f - q_r;

    /*reaction 53: OH + C2H6 <=> C2H5 + H2O */
    phi_f = sc[4]*sc[18];
    k_f = 1e-06 * 3.54e+06*exp(2.12*tc[0]-437.846/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[5];
    Kc = exp((g_RT[4] + g_RT[18]) - (g_RT[17] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[52] = q_f - q_r;

    /*reaction 54: HO2 + CH2 <=> OH + CH2O */
    phi_f = sc[6]*sc[7];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[14];
    Kc = exp((g_RT[6] + g_RT[7]) - (g_RT[4] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[53] = q_f - q_r;

    /*reaction 55: HO2 + CH3 <=> O2 + CH4 */
    phi_f = sc[6]*sc[9];
    k_f = 1e-06 * 1e+12;
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[10];
    Kc = exp((g_RT[6] + g_RT[9]) - (g_RT[3] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[54] = q_f - q_r;

    /*reaction 56: HO2 + CH3 <=> OH + CH3O */
    phi_f = sc[6]*sc[9];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[15];
    Kc = exp((g_RT[6] + g_RT[9]) - (g_RT[4] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[55] = q_f - q_r;

    /*reaction 57: HO2 + CO <=> OH + CO2 */
    phi_f = sc[6]*sc[11];
    k_f = 1e-06 * 1.5e+14*exp(-11877.2/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = exp((g_RT[6] + g_RT[11]) - (g_RT[4] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[56] = q_f - q_r;

    /*reaction 58: CH2 + O2 <=> OH + HCO */
    phi_f = sc[7]*sc[3];
    k_f = 1e-06 * 1.32e+13*exp(-754.907/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[13];
    Kc = exp((g_RT[7] + g_RT[3]) - (g_RT[4] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[57] = q_f - q_r;

    /*reaction 59: CH2 + H2 <=> H + CH3 */
    phi_f = sc[7]*sc[0];
    k_f = 1e-06 * 500000*exp(2*tc[0]-3638.65/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[9];
    Kc = exp((g_RT[7] + g_RT[0]) - (g_RT[1] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[58] = q_f - q_r;

    /*reaction 60: CH2 + CH3 <=> H + C2H4 */
    phi_f = sc[7]*sc[9];
    k_f = 1e-06 * 4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = exp((g_RT[7] + g_RT[9]) - (g_RT[1] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[59] = q_f - q_r;

    /*reaction 61: CH2 + CH4 <=> 2 CH3 */
    phi_f = sc[7]*sc[10];
    k_f = 1e-06 * 2.46e+06*exp(2*tc[0]-4162.05/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[9];
    Kc = exp((g_RT[7] + g_RT[10]) - (2 * g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[60] = q_f - q_r;

    /*reaction 62: CH2(S) + N2 <=> CH2 + N2 */
    phi_f = sc[8]*sc[19];
    k_f = 1e-06 * 1.5e+13*exp(-301.963/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[19];
    Kc = exp((g_RT[8] + g_RT[19]) - (g_RT[7] + g_RT[19]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[61] = q_f - q_r;

    /*reaction 63: CH2(S) + AR <=> CH2 + AR */
    phi_f = sc[8]*sc[20];
    k_f = 1e-06 * 9e+12*exp(-301.963/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[20];
    Kc = exp((g_RT[8] + g_RT[20]) - (g_RT[7] + g_RT[20]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[62] = q_f - q_r;

    /*reaction 64: CH2(S) + O2 <=> H + OH + CO */
    phi_f = sc[8]*sc[3];
    k_f = 1e-06 * 2.8e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[4]*sc[11];
    Kc = refC * exp((g_RT[8] + g_RT[3]) - (g_RT[1] + g_RT[4] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[63] = q_f - q_r;

    /*reaction 65: CH2(S) + O2 <=> CO + H2O */
    phi_f = sc[8]*sc[3];
    k_f = 1e-06 * 1.2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = exp((g_RT[8] + g_RT[3]) - (g_RT[11] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[64] = q_f - q_r;

    /*reaction 66: CH2(S) + H2 <=> CH3 + H */
    phi_f = sc[8]*sc[0];
    k_f = 1e-06 * 7e+13;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[1];
    Kc = exp((g_RT[8] + g_RT[0]) - (g_RT[9] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[65] = q_f - q_r;

    /*reaction 67: CH2(S) + H2O <=> CH2 + H2O */
    phi_f = sc[8]*sc[5];
    k_f = 1e-06 * 3e+13;
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[5];
    Kc = exp((g_RT[8] + g_RT[5]) - (g_RT[7] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[66] = q_f - q_r;

    /*reaction 68: CH2(S) + CH3 <=> H + C2H4 */
    phi_f = sc[8]*sc[9];
    k_f = 1e-06 * 1.2e+13*exp(+286.865/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = exp((g_RT[8] + g_RT[9]) - (g_RT[1] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[67] = q_f - q_r;

    /*reaction 69: CH2(S) + CH4 <=> 2 CH3 */
    phi_f = sc[8]*sc[10];
    k_f = 1e-06 * 1.6e+13*exp(+286.865/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[9];
    Kc = exp((g_RT[8] + g_RT[10]) - (2 * g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[68] = q_f - q_r;

    /*reaction 70: CH2(S) + CO <=> CH2 + CO */
    phi_f = sc[8]*sc[11];
    k_f = 1e-06 * 9e+12;
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[11];
    Kc = exp((g_RT[8] + g_RT[11]) - (g_RT[7] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[69] = q_f - q_r;

    /*reaction 71: CH2(S) + CO2 <=> CH2 + CO2 */
    phi_f = sc[8]*sc[12];
    k_f = 1e-06 * 7e+12;
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[12];
    Kc = exp((g_RT[8] + g_RT[12]) - (g_RT[7] + g_RT[12]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[70] = q_f - q_r;

    /*reaction 72: CH2(S) + CO2 <=> CO + CH2O */
    phi_f = sc[8]*sc[12];
    k_f = 1e-06 * 1.4e+13;
    q_f = phi_f * k_f;
    phi_r = sc[11]*sc[14];
    Kc = exp((g_RT[8] + g_RT[12]) - (g_RT[11] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[71] = q_f - q_r;

    /*reaction 73: CH3 + O2 <=> O + CH3O */
    phi_f = sc[9]*sc[3];
    k_f = 1e-06 * 2.675e+13*exp(-14494.2/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[15];
    Kc = exp((g_RT[9] + g_RT[3]) - (g_RT[2] + g_RT[15]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[72] = q_f - q_r;

    /*reaction 74: CH3 + O2 <=> OH + CH2O */
    phi_f = sc[9]*sc[3];
    k_f = 1e-06 * 3.6e+10*exp(-4499.25/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[14];
    Kc = exp((g_RT[9] + g_RT[3]) - (g_RT[4] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[73] = q_f - q_r;

    /*reaction 75: 2 CH3 (+M) <=> C2H6 (+M) */
    phi_f = sc[9]*sc[9];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
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
    phi_r = sc[18];
    Kc = 1.0 / (refC) * exp((2 * g_RT[9]) - (g_RT[18]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[74] = q_f - q_r;

    /*reaction 76: 2 CH3 <=> H + C2H5 */
    phi_f = sc[9]*sc[9];
    k_f = 1e-06 * 4.99e+12*exp(0.1*tc[0]-5334.68/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = exp((2 * g_RT[9]) - (g_RT[1] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[75] = q_f - q_r;

    /*reaction 77: CH3 + HCO <=> CH4 + CO */
    phi_f = sc[9]*sc[13];
    k_f = 1e-06 * 2.648e+13;
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[11];
    Kc = exp((g_RT[9] + g_RT[13]) - (g_RT[10] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[76] = q_f - q_r;

    /*reaction 78: CH3 + CH2O <=> HCO + CH4 */
    phi_f = sc[9]*sc[14];
    k_f = 1e-06 * 3320*exp(2.81*tc[0]-2949.17/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[10];
    Kc = exp((g_RT[9] + g_RT[14]) - (g_RT[13] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[77] = q_f - q_r;

    /*reaction 79: CH3 + C2H6 <=> C2H5 + CH4 */
    phi_f = sc[9]*sc[18];
    k_f = 1e-06 * 6.14e+06*exp(1.74*tc[0]-5259.18/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[10];
    Kc = exp((g_RT[9] + g_RT[18]) - (g_RT[17] + g_RT[10]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[78] = q_f - q_r;

    /*reaction 80: HCO + H2O <=> H + CO + H2O */
    phi_f = sc[13]*sc[5];
    k_f = 1e-06 * 2.244e+18*exp(-1*tc[0]-8555.61/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[11]*sc[5];
    Kc = refC * exp((g_RT[13] + g_RT[5]) - (g_RT[1] + g_RT[11] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[79] = q_f - q_r;

    /*reaction 81: HCO + M <=> H + CO + M */
    phi_f = sc[13];
    alpha = mixture + sc[0] + -1*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18];
    k_f = 1e-06 * alpha * 1.87e+17*exp(-1*tc[0]-8555.61/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[11];
    Kc = refC * exp((g_RT[13]) - (g_RT[1] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[80] = q_f - q_r;

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    phi_f = sc[13]*sc[3];
    k_f = 1e-06 * 7.6e+12*exp(-201.309/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[11];
    Kc = exp((g_RT[13] + g_RT[3]) - (g_RT[6] + g_RT[11]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[81] = q_f - q_r;

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    phi_f = sc[15]*sc[3];
    k_f = 1e-06 * 4.28e-13*exp(7.6*tc[0]+1776.55/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[14];
    Kc = exp((g_RT[15] + g_RT[3]) - (g_RT[6] + g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[82] = q_f - q_r;

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    phi_f = sc[17]*sc[3];
    k_f = 1e-06 * 8.4e+11*exp(-1950.18/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[16];
    Kc = exp((g_RT[17] + g_RT[3]) - (g_RT[6] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[83] = q_f - q_r;

    return;
}


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

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 21; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*reaction 1: O + H + M <=> OH + M */
    phi_f = sc[2]*sc[1];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
    k_f = 1e-12 * alpha * 5e+17*exp(-1*tc[0]);
    q_f[0] = phi_f * k_f;
    phi_r = sc[4];
    Kc = 1.0 / (refC) * exp((g_RT[2] + g_RT[1]) - (g_RT[4]));
    k_r = k_f / Kc;
    q_r[0] = phi_r * k_r;

    /*reaction 2: O + H2 <=> H + OH */
    phi_f = sc[2]*sc[0];
    k_f = 1e-06 * 50000*exp(2.67*tc[0]-3165.58/tc[1]);
    q_f[1] = phi_f * k_f;
    phi_r = sc[1]*sc[4];
    Kc = exp((g_RT[2] + g_RT[0]) - (g_RT[1] + g_RT[4]));
    k_r = k_f / Kc;
    q_r[1] = phi_r * k_r;

    /*reaction 3: O + HO2 <=> OH + O2 */
    phi_f = sc[2]*sc[6];
    k_f = 1e-06 * 2e+13;
    q_f[2] = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = exp((g_RT[2] + g_RT[6]) - (g_RT[4] + g_RT[3]));
    k_r = k_f / Kc;
    q_r[2] = phi_r * k_r;

    /*reaction 4: O + CH2 <=> H + HCO */
    phi_f = sc[2]*sc[7];
    k_f = 1e-06 * 8e+13;
    q_f[3] = phi_f * k_f;
    phi_r = sc[1]*sc[13];
    Kc = exp((g_RT[2] + g_RT[7]) - (g_RT[1] + g_RT[13]));
    k_r = k_f / Kc;
    q_r[3] = phi_r * k_r;

    /*reaction 5: O + CH2(S) <=> H + HCO */
    phi_f = sc[2]*sc[8];
    k_f = 1e-06 * 1.5e+13;
    q_f[4] = phi_f * k_f;
    phi_r = sc[1]*sc[13];
    Kc = exp((g_RT[2] + g_RT[8]) - (g_RT[1] + g_RT[13]));
    k_r = k_f / Kc;
    q_r[4] = phi_r * k_r;

    /*reaction 6: O + CH3 <=> H + CH2O */
    phi_f = sc[2]*sc[9];
    k_f = 1e-06 * 8.43e+13;
    q_f[5] = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = exp((g_RT[2] + g_RT[9]) - (g_RT[1] + g_RT[14]));
    k_r = k_f / Kc;
    q_r[5] = phi_r * k_r;

    /*reaction 7: O + CH4 <=> OH + CH3 */
    phi_f = sc[2]*sc[10];
    k_f = 1e-06 * 1.02e+09*exp(1.5*tc[0]-4328.13/tc[1]);
    q_f[6] = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = exp((g_RT[2] + g_RT[10]) - (g_RT[4] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[6] = phi_r * k_r;

    /*reaction 8: O + CO + M <=> CO2 + M */
    phi_f = sc[2]*sc[11];
    alpha = mixture + sc[0] + 5*sc[3] + 5*sc[5] + sc[10] + 0.5*sc[11] + 2.5*sc[12] + 2*sc[18] + -0.5*sc[20];
    k_f = 1e-12 * alpha * 6.02e+14*exp(-1509.81/tc[1]);
    q_f[7] = phi_f * k_f;
    phi_r = sc[12];
    Kc = 1.0 / (refC) * exp((g_RT[2] + g_RT[11]) - (g_RT[12]));
    k_r = k_f / Kc;
    q_r[7] = phi_r * k_r;

    /*reaction 9: O + HCO <=> OH + CO */
    phi_f = sc[2]*sc[13];
    k_f = 1e-06 * 3e+13;
    q_f[8] = phi_f * k_f;
    phi_r = sc[4]*sc[11];
    Kc = exp((g_RT[2] + g_RT[13]) - (g_RT[4] + g_RT[11]));
    k_r = k_f / Kc;
    q_r[8] = phi_r * k_r;

    /*reaction 10: O + HCO <=> H + CO2 */
    phi_f = sc[2]*sc[13];
    k_f = 1e-06 * 3e+13;
    q_f[9] = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = exp((g_RT[2] + g_RT[13]) - (g_RT[1] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[9] = phi_r * k_r;

    /*reaction 11: O + CH2O <=> OH + HCO */
    phi_f = sc[2]*sc[14];
    k_f = 1e-06 * 3.9e+13*exp(-1781.58/tc[1]);
    q_f[10] = phi_f * k_f;
    phi_r = sc[4]*sc[13];
    Kc = exp((g_RT[2] + g_RT[14]) - (g_RT[4] + g_RT[13]));
    k_r = k_f / Kc;
    q_r[10] = phi_r * k_r;

    /*reaction 12: O + C2H4 <=> CH3 + HCO */
    phi_f = sc[2]*sc[16];
    k_f = 1e-06 * 1.92e+07*exp(1.83*tc[0]-110.72/tc[1]);
    q_f[11] = phi_f * k_f;
    phi_r = sc[9]*sc[13];
    Kc = exp((g_RT[2] + g_RT[16]) - (g_RT[9] + g_RT[13]));
    k_r = k_f / Kc;
    q_r[11] = phi_r * k_r;

    /*reaction 13: O + C2H5 <=> CH3 + CH2O */
    phi_f = sc[2]*sc[17];
    k_f = 1e-06 * 1.32e+14;
    q_f[12] = phi_f * k_f;
    phi_r = sc[9]*sc[14];
    Kc = exp((g_RT[2] + g_RT[17]) - (g_RT[9] + g_RT[14]));
    k_r = k_f / Kc;
    q_r[12] = phi_r * k_r;

    /*reaction 14: O + C2H6 <=> OH + C2H5 */
    phi_f = sc[2]*sc[18];
    k_f = 1e-06 * 8.98e+07*exp(1.92*tc[0]-2863.61/tc[1]);
    q_f[13] = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = exp((g_RT[2] + g_RT[18]) - (g_RT[4] + g_RT[17]));
    k_r = k_f / Kc;
    q_r[13] = phi_r * k_r;

    /*reaction 15: O2 + CO <=> O + CO2 */
    phi_f = sc[3]*sc[11];
    k_f = 1e-06 * 2.5e+12*exp(-24056.4/tc[1]);
    q_f[14] = phi_f * k_f;
    phi_r = sc[2]*sc[12];
    Kc = exp((g_RT[3] + g_RT[11]) - (g_RT[2] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[14] = phi_r * k_r;

    /*reaction 16: O2 + CH2O <=> HO2 + HCO */
    phi_f = sc[3]*sc[14];
    k_f = 1e-06 * 1e+14*exp(-20130.9/tc[1]);
    q_f[15] = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = exp((g_RT[3] + g_RT[14]) - (g_RT[6] + g_RT[13]));
    k_r = k_f / Kc;
    q_r[15] = phi_r * k_r;

    /*reaction 17: H + O2 + M <=> HO2 + M */
    phi_f = sc[1]*sc[3];
    alpha = mixture + -1*sc[3] + -1*sc[5] + -0.25*sc[11] + 0.5*sc[12] + 0.5*sc[18] + -1*sc[19] + -1*sc[20];
    k_f = 1e-12 * alpha * 2.8e+18*exp(-0.86*tc[0]);
    q_f[16] = phi_f * k_f;
    phi_r = sc[6];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));
    k_r = k_f / Kc;
    q_r[16] = phi_r * k_r;

    /*reaction 18: H + 2 O2 <=> HO2 + O2 */
    phi_f = sc[1]*sc[3]*sc[3];
    k_f = 1e-12 * 3e+20*exp(-1.72*tc[0]);
    q_f[17] = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = 1.0 / (refC) * exp((g_RT[1] + 2 * g_RT[3]) - (g_RT[6] + g_RT[3]));
    k_r = k_f / Kc;
    q_r[17] = phi_r * k_r;

    /*reaction 19: H + O2 + H2O <=> HO2 + H2O */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = 1e-12 * 9.38e+18*exp(-0.76*tc[0]);
    q_f[18] = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r[18] = phi_r * k_r;

    /*reaction 20: H + O2 + N2 <=> HO2 + N2 */
    phi_f = sc[1]*sc[3]*sc[19];
    k_f = 1e-12 * 3.75e+20*exp(-1.72*tc[0]);
    q_f[19] = phi_f * k_f;
    phi_r = sc[6]*sc[19];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[19]) - (g_RT[6] + g_RT[19]));
    k_r = k_f / Kc;
    q_r[19] = phi_r * k_r;

    /*reaction 21: H + O2 + AR <=> HO2 + AR */
    phi_f = sc[1]*sc[3]*sc[20];
    k_f = 1e-12 * 7e+17*exp(-0.8*tc[0]);
    q_f[20] = phi_f * k_f;
    phi_r = sc[6]*sc[20];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[20]) - (g_RT[6] + g_RT[20]));
    k_r = k_f / Kc;
    q_r[20] = phi_r * k_r;

    /*reaction 22: H + O2 <=> O + OH */
    phi_f = sc[1]*sc[3];
    k_f = 1e-06 * 8.3e+13*exp(-7253.65/tc[1]);
    q_f[21] = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));
    k_r = k_f / Kc;
    q_r[21] = phi_r * k_r;

    /*reaction 23: 2 H + M <=> H2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[0] + -1*sc[5] + sc[10] + -1*sc[12] + 2*sc[18] + -0.37*sc[20];
    k_f = 1e-12 * alpha * 1e+18*exp(-1*tc[0]);
    q_f[22] = phi_f * k_f;
    phi_r = sc[0];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));
    k_r = k_f / Kc;
    q_r[22] = phi_r * k_r;

    /*reaction 24: 2 H + H2 <=> 2 H2 */
    phi_f = sc[1]*sc[1]*sc[0];
    k_f = 1e-12 * 9e+16*exp(-0.6*tc[0]);
    q_f[23] = phi_f * k_f;
    phi_r = sc[0]*sc[0];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));
    k_r = k_f / Kc;
    q_r[23] = phi_r * k_r;

    /*reaction 25: 2 H + H2O <=> H2 + H2O */
    phi_f = sc[1]*sc[1]*sc[5];
    k_f = 1e-12 * 6e+19*exp(-1.25*tc[0]);
    q_f[24] = phi_f * k_f;
    phi_r = sc[0]*sc[5];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));
    k_r = k_f / Kc;
    q_r[24] = phi_r * k_r;

    /*reaction 26: 2 H + CO2 <=> H2 + CO2 */
    phi_f = sc[1]*sc[1]*sc[12];
    k_f = 1e-12 * 5.5e+20*exp(-2*tc[0]);
    q_f[25] = phi_f * k_f;
    phi_r = sc[0]*sc[12];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[12]) - (g_RT[0] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[25] = phi_r * k_r;

    /*reaction 27: H + OH + M <=> H2O + M */
    phi_f = sc[1]*sc[4];
    alpha = mixture + -0.27*sc[0] + 2.65*sc[5] + sc[10] + 2*sc[18] + -0.62*sc[20];
    k_f = 1e-12 * alpha * 2.2e+22*exp(-2*tc[0]);
    q_f[26] = phi_f * k_f;
    phi_r = sc[5];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));
    k_r = k_f / Kc;
    q_r[26] = phi_r * k_r;

    /*reaction 28: H + HO2 <=> O2 + H2 */
    phi_f = sc[1]*sc[6];
    k_f = 1e-06 * 2.8e+13*exp(-537.494/tc[1]);
    q_f[27] = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[27] = phi_r * k_r;

    /*reaction 29: H + HO2 <=> 2 OH */
    phi_f = sc[1]*sc[6];
    k_f = 1e-06 * 1.34e+14*exp(-319.577/tc[1]);
    q_f[28] = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));
    k_r = k_f / Kc;
    q_r[28] = phi_r * k_r;

    /*reaction 30: H + CH2 (+M) <=> CH3 (+M) */
    phi_f = sc[1]*sc[7];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
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
    q_f[29] = phi_f * k_f;
    phi_r = sc[9];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[7]) - (g_RT[9]));
    k_r = k_f / Kc;
    q_r[29] = phi_r * k_r;

    /*reaction 31: H + CH3 (+M) <=> CH4 (+M) */
    phi_f = sc[1]*sc[9];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
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
    q_f[30] = phi_f * k_f;
    phi_r = sc[10];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[9]) - (g_RT[10]));
    k_r = k_f / Kc;
    q_r[20] = phi_r * k_r;

    /*reaction 32: H + CH4 <=> CH3 + H2 */
    phi_f = sc[1]*sc[10];
    k_f = 1e-06 * 6.6e+08*exp(1.62*tc[0]-5455.46/tc[1]);
    q_f[31] = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = exp((g_RT[1] + g_RT[10]) - (g_RT[9] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[31] = phi_r * k_r;

    /*reaction 33: H + HCO (+M) <=> CH2O (+M) */
    phi_f = sc[1]*sc[13];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
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
    q_f[32] = phi_f * k_f;
    phi_r = sc[14];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[13]) - (g_RT[14]));
    k_r = k_f / Kc;
    q_r[32] = phi_r * k_r;

    /*reaction 34: H + HCO <=> H2 + CO */
    phi_f = sc[1]*sc[13];
    k_f = 1e-06 * 7.34e+13;
    q_f[33] = phi_f * k_f;
    phi_r = sc[0]*sc[11];
    Kc = exp((g_RT[1] + g_RT[13]) - (g_RT[0] + g_RT[11]));
    k_r = k_f / Kc;
    q_r[33] = phi_r * k_r;

    /*reaction 35: H + CH2O (+M) <=> CH3O (+M) */
    phi_f = sc[1]*sc[14];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18];
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
    q_f[34] = phi_f * k_f;
    phi_r = sc[15];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[14]) - (g_RT[15]));
    k_r = k_f / Kc;
    q_r[34] = phi_r * k_r;

    /*reaction 36: H + CH2O <=> HCO + H2 */
    phi_f = sc[1]*sc[14];
    k_f = 1e-06 * 2.3e+10*exp(1.05*tc[0]-1648.21/tc[1]);
    q_f[35] = phi_f * k_f;
    phi_r = sc[13]*sc[0];
    Kc = exp((g_RT[1] + g_RT[14]) - (g_RT[13] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[35] = phi_r * k_r;

    /*reaction 37: H + CH3O <=> OH + CH3 */
    phi_f = sc[1]*sc[15];
    k_f = 1e-06 * 3.2e+13;
    q_f[36] = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = exp((g_RT[1] + g_RT[15]) - (g_RT[4] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[36] = phi_r * k_r;

    /*reaction 38: H + C2H4 (+M) <=> C2H5 (+M) */
    phi_f = sc[1]*sc[16];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
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
    q_f[37] = phi_f * k_f;
    phi_r = sc[17];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[16]) - (g_RT[17]));
    k_r = k_f / Kc;
    q_r[37] = phi_r * k_r;

    /*reaction 39: H + C2H5 (+M) <=> C2H6 (+M) */
    phi_f = sc[1]*sc[17];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
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
    q_f[38] = phi_f * k_f;
    phi_r = sc[18];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[18]));
    k_r = k_f / Kc;
    q_r[38] = phi_r * k_r;

    /*reaction 40: H + C2H6 <=> C2H5 + H2 */
    phi_f = sc[1]*sc[18];
    k_f = 1e-06 * 1.15e+08*exp(1.9*tc[0]-3789.63/tc[1]);
    q_f[39] = phi_f * k_f;
    phi_r = sc[17]*sc[0];
    Kc = exp((g_RT[1] + g_RT[18]) - (g_RT[17] + g_RT[0]));
    k_r = k_f / Kc;
    q_r[39] = phi_r * k_r;

    /*reaction 41: H2 + CO (+M) <=> CH2O (+M) */
    phi_f = sc[0]*sc[11];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
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
    q_f[40] = phi_f * k_f;
    phi_r = sc[14];
    Kc = 1.0 / (refC) * exp((g_RT[0] + g_RT[11]) - (g_RT[14]));
    k_r = k_f / Kc;
    q_r[40] = phi_r * k_r;

    /*reaction 42: OH + H2 <=> H + H2O */
    phi_f = sc[4]*sc[0];
    k_f = 1e-06 * 2.16e+08*exp(1.51*tc[0]-1726.22/tc[1]);
    q_f[41] = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    Kc = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));
    k_r = k_f / Kc;
    q_r[41] = phi_r * k_r;

    /*reaction 43: 2 OH <=> O + H2O */
    phi_f = sc[4]*sc[4];
    k_f = 1e-06 * 35700*exp(2.4*tc[0]+1061.9/tc[1]);
    q_f[42] = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));
    k_r = k_f / Kc;
    q_r[42] = phi_r * k_r;

    /*reaction 44: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = 1e-06 * 2.9e+13*exp(+251.636/tc[1]);
    q_f[43] = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
    k_r = k_f / Kc;
    q_r[43] = phi_r * k_r;

    /*reaction 45: OH + CH2 <=> H + CH2O */
    phi_f = sc[4]*sc[7];
    k_f = 1e-06 * 2e+13;
    q_f[44] = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = exp((g_RT[4] + g_RT[7]) - (g_RT[1] + g_RT[14]));
    k_r = k_f / Kc;
    q_r[44] = phi_r * k_r;

    /*reaction 46: OH + CH2(S) <=> H + CH2O */
    phi_f = sc[4]*sc[8];
    k_f = 1e-06 * 3e+13;
    q_f[45] = phi_f * k_f;
    phi_r = sc[1]*sc[14];
    Kc = exp((g_RT[4] + g_RT[8]) - (g_RT[1] + g_RT[14]));
    k_r = k_f / Kc;
    q_r[45] = phi_r * k_r;

    /*reaction 47: OH + CH3 <=> CH2 + H2O */
    phi_f = sc[4]*sc[9];
    k_f = 1e-06 * 5.6e+07*exp(1.6*tc[0]-2727.73/tc[1]);
    q_f[46] = phi_f * k_f;
    phi_r = sc[7]*sc[5];
    Kc = exp((g_RT[4] + g_RT[9]) - (g_RT[7] + g_RT[5]));
    k_r = k_f / Kc;
    q_r[46] = phi_r * k_r;

    /*reaction 48: OH + CH3 <=> CH2(S) + H2O */
    phi_f = sc[4]*sc[9];
    k_f = 1e-06 * 2.501e+13;
    q_f[47] = phi_f * k_f;
    phi_r = sc[8]*sc[5];
    Kc = exp((g_RT[4] + g_RT[9]) - (g_RT[8] + g_RT[5]));
    k_r = k_f / Kc;
    q_r[47] = phi_r * k_r;

    /*reaction 49: OH + CH4 <=> CH3 + H2O */
    phi_f = sc[4]*sc[10];
    k_f = 1e-06 * 1e+08*exp(1.6*tc[0]-1570.21/tc[1]);
    q_f[48] = phi_f * k_f;
    phi_r = sc[9]*sc[5];
    Kc = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[5]));
    k_r = k_f / Kc;
    q_r[48] = phi_r * k_r;

    /*reaction 50: OH + CO <=> H + CO2 */
    phi_f = sc[4]*sc[11];
    k_f = 1e-06 * 4.76e+07*exp(1.228*tc[0]-35.229/tc[1]);
    q_f[49] = phi_f * k_f;
    phi_r = sc[1]*sc[12];
    Kc = exp((g_RT[4] + g_RT[11]) - (g_RT[1] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[49] = phi_r * k_r;

    /*reaction 51: OH + HCO <=> H2O + CO */
    phi_f = sc[4]*sc[13];
    k_f = 1e-06 * 5e+13;
    q_f[50] = phi_f * k_f;
    phi_r = sc[5]*sc[11];
    Kc = exp((g_RT[4] + g_RT[13]) - (g_RT[5] + g_RT[11]));
    k_r = k_f / Kc;
    q_r[50] = phi_r * k_r;

    /*reaction 52: OH + CH2O <=> HCO + H2O */
    phi_f = sc[4]*sc[14];
    k_f = 1e-06 * 3.43e+09*exp(1.18*tc[0]+224.962/tc[1]);
    q_f[51] = phi_f * k_f;
    phi_r = sc[13]*sc[5];
    Kc = exp((g_RT[4] + g_RT[14]) - (g_RT[13] + g_RT[5]));
    k_r = k_f / Kc;
    q_r[51] = phi_r * k_r;

    /*reaction 53: OH + C2H6 <=> C2H5 + H2O */
    phi_f = sc[4]*sc[18];
    k_f = 1e-06 * 3.54e+06*exp(2.12*tc[0]-437.846/tc[1]);
    q_f[52] = phi_f * k_f;
    phi_r = sc[17]*sc[5];
    Kc = exp((g_RT[4] + g_RT[18]) - (g_RT[17] + g_RT[5]));
    k_r = k_f / Kc;
    q_r[52] = phi_r * k_r;

    /*reaction 54: HO2 + CH2 <=> OH + CH2O */
    phi_f = sc[6]*sc[7];
    k_f = 1e-06 * 2e+13;
    q_f[53] = phi_f * k_f;
    phi_r = sc[4]*sc[14];
    Kc = exp((g_RT[6] + g_RT[7]) - (g_RT[4] + g_RT[14]));
    k_r = k_f / Kc;
    q_r[53] = phi_r * k_r;

    /*reaction 55: HO2 + CH3 <=> O2 + CH4 */
    phi_f = sc[6]*sc[9];
    k_f = 1e-06 * 1e+12;
    q_f[54] = phi_f * k_f;
    phi_r = sc[3]*sc[10];
    Kc = exp((g_RT[6] + g_RT[9]) - (g_RT[3] + g_RT[10]));
    k_r = k_f / Kc;
    q_r[54] = phi_r * k_r;

    /*reaction 56: HO2 + CH3 <=> OH + CH3O */
    phi_f = sc[6]*sc[9];
    k_f = 1e-06 * 2e+13;
    q_f[55] = phi_f * k_f;
    phi_r = sc[4]*sc[15];
    Kc = exp((g_RT[6] + g_RT[9]) - (g_RT[4] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[55] = phi_r * k_r;

    /*reaction 57: HO2 + CO <=> OH + CO2 */
    phi_f = sc[6]*sc[11];
    k_f = 1e-06 * 1.5e+14*exp(-11877.2/tc[1]);
    q_f[56] = phi_f * k_f;
    phi_r = sc[4]*sc[12];
    Kc = exp((g_RT[6] + g_RT[11]) - (g_RT[4] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[56] = phi_r * k_r;

    /*reaction 58: CH2 + O2 <=> OH + HCO */
    phi_f = sc[7]*sc[3];
    k_f = 1e-06 * 1.32e+13*exp(-754.907/tc[1]);
    q_f[57] = phi_f * k_f;
    phi_r = sc[4]*sc[13];
    Kc = exp((g_RT[7] + g_RT[3]) - (g_RT[4] + g_RT[13]));
    k_r = k_f / Kc;
    q_r[57] = phi_r * k_r;

    /*reaction 59: CH2 + H2 <=> H + CH3 */
    phi_f = sc[7]*sc[0];
    k_f = 1e-06 * 500000*exp(2*tc[0]-3638.65/tc[1]);
    q_f[58] = phi_f * k_f;
    phi_r = sc[1]*sc[9];
    Kc = exp((g_RT[7] + g_RT[0]) - (g_RT[1] + g_RT[9]));
    k_r = k_f / Kc;
    q_r[58] = phi_r * k_r;

    /*reaction 60: CH2 + CH3 <=> H + C2H4 */
    phi_f = sc[7]*sc[9];
    k_f = 1e-06 * 4e+13;
    q_f[59] = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = exp((g_RT[7] + g_RT[9]) - (g_RT[1] + g_RT[16]));
    k_r = k_f / Kc;
    q_r[59] = phi_r * k_r;

    /*reaction 61: CH2 + CH4 <=> 2 CH3 */
    phi_f = sc[7]*sc[10];
    k_f = 1e-06 * 2.46e+06*exp(2*tc[0]-4162.05/tc[1]);
    q_f[60] = phi_f * k_f;
    phi_r = sc[9]*sc[9];
    Kc = exp((g_RT[7] + g_RT[10]) - (2 * g_RT[9]));
    k_r = k_f / Kc;
    q_r[60] = phi_r * k_r;

    /*reaction 62: CH2(S) + N2 <=> CH2 + N2 */
    phi_f = sc[8]*sc[19];
    k_f = 1e-06 * 1.5e+13*exp(-301.963/tc[1]);
    q_f[61] = phi_f * k_f;
    phi_r = sc[7]*sc[19];
    Kc = exp((g_RT[8] + g_RT[19]) - (g_RT[7] + g_RT[19]));
    k_r = k_f / Kc;
    q_r[61] = phi_r * k_r;

    /*reaction 63: CH2(S) + AR <=> CH2 + AR */
    phi_f = sc[8]*sc[20];
    k_f = 1e-06 * 9e+12*exp(-301.963/tc[1]);
    q_f[62] = phi_f * k_f;
    phi_r = sc[7]*sc[20];
    Kc = exp((g_RT[8] + g_RT[20]) - (g_RT[7] + g_RT[20]));
    k_r = k_f / Kc;
    q_r[62] = phi_r * k_r;

    /*reaction 64: CH2(S) + O2 <=> H + OH + CO */
    phi_f = sc[8]*sc[3];
    k_f = 1e-06 * 2.8e+13;
    q_f[63] = phi_f * k_f;
    phi_r = sc[1]*sc[4]*sc[11];
    Kc = refC * exp((g_RT[8] + g_RT[3]) - (g_RT[1] + g_RT[4] + g_RT[11]));
    k_r = k_f / Kc;
    q_r[63] = phi_r * k_r;

    /*reaction 65: CH2(S) + O2 <=> CO + H2O */
    phi_f = sc[8]*sc[3];
    k_f = 1e-06 * 1.2e+13;
    q_f[64] = phi_f * k_f;
    phi_r = sc[11]*sc[5];
    Kc = exp((g_RT[8] + g_RT[3]) - (g_RT[11] + g_RT[5]));
    k_r = k_f / Kc;
    q_r[64] = phi_r * k_r;

    /*reaction 66: CH2(S) + H2 <=> CH3 + H */
    phi_f = sc[8]*sc[0];
    k_f = 1e-06 * 7e+13;
    q_f[65] = phi_f * k_f;
    phi_r = sc[9]*sc[1];
    Kc = exp((g_RT[8] + g_RT[0]) - (g_RT[9] + g_RT[1]));
    k_r = k_f / Kc;
    q_r[65] = phi_r * k_r;

    /*reaction 67: CH2(S) + H2O <=> CH2 + H2O */
    phi_f = sc[8]*sc[5];
    k_f = 1e-06 * 3e+13;
    q_f[66] = phi_f * k_f;
    phi_r = sc[7]*sc[5];
    Kc = exp((g_RT[8] + g_RT[5]) - (g_RT[7] + g_RT[5]));
    k_r = k_f / Kc;
    q_r[66] = phi_r * k_r;

    /*reaction 68: CH2(S) + CH3 <=> H + C2H4 */
    phi_f = sc[8]*sc[9];
    k_f = 1e-06 * 1.2e+13*exp(+286.865/tc[1]);
    q_f[67] = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = exp((g_RT[8] + g_RT[9]) - (g_RT[1] + g_RT[16]));
    k_r = k_f / Kc;
    q_r[67] = phi_r * k_r;

    /*reaction 69: CH2(S) + CH4 <=> 2 CH3 */
    phi_f = sc[8]*sc[10];
    k_f = 1e-06 * 1.6e+13*exp(+286.865/tc[1]);
    q_f[68] = phi_f * k_f;
    phi_r = sc[9]*sc[9];
    Kc = exp((g_RT[8] + g_RT[10]) - (2 * g_RT[9]));
    k_r = k_f / Kc;
    q_r[68] = phi_r * k_r;

    /*reaction 70: CH2(S) + CO <=> CH2 + CO */
    phi_f = sc[8]*sc[11];
    k_f = 1e-06 * 9e+12;
    q_f[69] = phi_f * k_f;
    phi_r = sc[7]*sc[11];
    Kc = exp((g_RT[8] + g_RT[11]) - (g_RT[7] + g_RT[11]));
    k_r = k_f / Kc;
    q_r[69] = phi_r * k_r;

    /*reaction 71: CH2(S) + CO2 <=> CH2 + CO2 */
    phi_f = sc[8]*sc[12];
    k_f = 1e-06 * 7e+12;
    q_f[70] = phi_f * k_f;
    phi_r = sc[7]*sc[12];
    Kc = exp((g_RT[8] + g_RT[12]) - (g_RT[7] + g_RT[12]));
    k_r = k_f / Kc;
    q_r[70] = phi_r * k_r;

    /*reaction 72: CH2(S) + CO2 <=> CO + CH2O */
    phi_f = sc[8]*sc[12];
    k_f = 1e-06 * 1.4e+13;
    q_f[71] = phi_f * k_f;
    phi_r = sc[11]*sc[14];
    Kc = exp((g_RT[8] + g_RT[12]) - (g_RT[11] + g_RT[14]));
    k_r = k_f / Kc;
    q_r[71] = phi_r * k_r;

    /*reaction 73: CH3 + O2 <=> O + CH3O */
    phi_f = sc[9]*sc[3];
    k_f = 1e-06 * 2.675e+13*exp(-14494.2/tc[1]);
    q_f[72] = phi_f * k_f;
    phi_r = sc[2]*sc[15];
    Kc = exp((g_RT[9] + g_RT[3]) - (g_RT[2] + g_RT[15]));
    k_r = k_f / Kc;
    q_r[72] = phi_r * k_r;

    /*reaction 74: CH3 + O2 <=> OH + CH2O */
    phi_f = sc[9]*sc[3];
    k_f = 1e-06 * 3.6e+10*exp(-4499.25/tc[1]);
    q_f[73] = phi_f * k_f;
    phi_r = sc[4]*sc[14];
    Kc = exp((g_RT[9] + g_RT[3]) - (g_RT[4] + g_RT[14]));
    k_r = k_f / Kc;
    q_r[73] = phi_r * k_r;

    /*reaction 75: 2 CH3 (+M) <=> C2H6 (+M) */
    phi_f = sc[9]*sc[9];
    alpha = mixture + sc[0] + 5*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18] + -0.3*sc[20];
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
    q_f[74] = phi_f * k_f;
    phi_r = sc[18];
    Kc = 1.0 / (refC) * exp((2 * g_RT[9]) - (g_RT[18]));
    k_r = k_f / Kc;
    q_r[74] = phi_r * k_r;

    /*reaction 76: 2 CH3 <=> H + C2H5 */
    phi_f = sc[9]*sc[9];
    k_f = 1e-06 * 4.99e+12*exp(0.1*tc[0]-5334.68/tc[1]);
    q_f[75] = phi_f * k_f;
    phi_r = sc[1]*sc[17];
    Kc = exp((2 * g_RT[9]) - (g_RT[1] + g_RT[17]));
    k_r = k_f / Kc;
    q_r[75] = phi_r * k_r;

    /*reaction 77: CH3 + HCO <=> CH4 + CO */
    phi_f = sc[9]*sc[13];
    k_f = 1e-06 * 2.648e+13;
    q_f[76] = phi_f * k_f;
    phi_r = sc[10]*sc[11];
    Kc = exp((g_RT[9] + g_RT[13]) - (g_RT[10] + g_RT[11]));
    k_r = k_f / Kc;
    q_r[76] = phi_r * k_r;

    /*reaction 78: CH3 + CH2O <=> HCO + CH4 */
    phi_f = sc[9]*sc[14];
    k_f = 1e-06 * 3320*exp(2.81*tc[0]-2949.17/tc[1]);
    q_f[77] = phi_f * k_f;
    phi_r = sc[13]*sc[10];
    Kc = exp((g_RT[9] + g_RT[14]) - (g_RT[13] + g_RT[10]));
    k_r = k_f / Kc;
    q_r[77] = phi_r * k_r;

    /*reaction 79: CH3 + C2H6 <=> C2H5 + CH4 */
    phi_f = sc[9]*sc[18];
    k_f = 1e-06 * 6.14e+06*exp(1.74*tc[0]-5259.18/tc[1]);
    q_f[78] = phi_f * k_f;
    phi_r = sc[17]*sc[10];
    Kc = exp((g_RT[9] + g_RT[18]) - (g_RT[17] + g_RT[10]));
    k_r = k_f / Kc;
    q_r[78] = phi_r * k_r;

    /*reaction 80: HCO + H2O <=> H + CO + H2O */
    phi_f = sc[13]*sc[5];
    k_f = 1e-06 * 2.244e+18*exp(-1*tc[0]-8555.61/tc[1]);
    q_f[79] = phi_f * k_f;
    phi_r = sc[1]*sc[11]*sc[5];
    Kc = refC * exp((g_RT[13] + g_RT[5]) - (g_RT[1] + g_RT[11] + g_RT[5]));
    k_r = k_f / Kc;
    q_r[79] = phi_r * k_r;

    /*reaction 81: HCO + M <=> H + CO + M */
    phi_f = sc[13];
    alpha = mixture + sc[0] + -1*sc[5] + sc[10] + 0.5*sc[11] + sc[12] + 2*sc[18];
    k_f = 1e-06 * alpha * 1.87e+17*exp(-1*tc[0]-8555.61/tc[1]);
    q_f[80] = phi_f * k_f;
    phi_r = sc[1]*sc[11];
    Kc = refC * exp((g_RT[13]) - (g_RT[1] + g_RT[11]));
    k_r = k_f / Kc;
    q_r[80] = phi_r * k_r;

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    phi_f = sc[13]*sc[3];
    k_f = 1e-06 * 7.6e+12*exp(-201.309/tc[1]);
    q_f[81] = phi_f * k_f;
    phi_r = sc[6]*sc[11];
    Kc = exp((g_RT[13] + g_RT[3]) - (g_RT[6] + g_RT[11]));
    k_r = k_f / Kc;
    q_r[81] = phi_r * k_r;

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    phi_f = sc[15]*sc[3];
    k_f = 1e-06 * 4.28e-13*exp(7.6*tc[0]+1776.55/tc[1]);
    q_f[82] = phi_f * k_f;
    phi_r = sc[6]*sc[14];
    Kc = exp((g_RT[15] + g_RT[3]) - (g_RT[6] + g_RT[14]));
    k_r = k_f / Kc;
    q_r[82] = phi_r * k_r;

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    phi_f = sc[17]*sc[3];
    k_f = 1e-06 * 8.4e+11*exp(-1950.18/tc[1]);
    q_f[83] = phi_f * k_f;
    phi_r = sc[6]*sc[16];
    Kc = exp((g_RT[17] + g_RT[3]) - (g_RT[6] + g_RT[16]));
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
    double T = tc[1];

    static double T_old = -1, species_old[21];

    if (T == T_old)
    {
        for (int i = 0; i < 21; i++)
            species[i] = species_old[i];
        return;
    }

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
        /*species 7: CH2 */
        species[7] =
            +4.60040401e+04 / tc[1]
            +2.20014682e+00
            -3.76267867e+00 * tc[0]
            -4.84436072e-04 * tc[1]
            -4.65816402e-07 * tc[2]
            +3.20909294e-10 * tc[3]
            -8.43708595e-14 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +5.04968163e+04 / tc[1]
            +4.96772308e+00
            -4.19860411e+00 * tc[0]
            +1.18330710e-03 * tc[1]
            -1.37216037e-06 * tc[2]
            +5.57346651e-10 * tc[3]
            -9.71573685e-14 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +1.64449988e+04 / tc[1]
            +2.06902607e+00
            -3.67359040e+00 * tc[0]
            -1.00547588e-03 * tc[1]
            -9.55036427e-07 * tc[2]
            +5.72597854e-10 * tc[3]
            -1.27192867e-13 * tc[4];
        /*species 10: CH4 */
        species[10] =
            -1.02466476e+04 / tc[1]
            +9.79117989e+00
            -5.14987613e+00 * tc[0]
            +6.83548940e-03 * tc[1]
            -8.19667665e-06 * tc[2]
            +4.03952522e-09 * tc[3]
            -8.33469780e-13 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.43440860e+04 / tc[1]
            +7.11241900e-02
            -3.57953347e+00 * tc[0]
            +3.05176840e-04 * tc[1]
            -1.69469055e-07 * tc[2]
            -7.55838237e-11 * tc[3]
            +4.52212249e-14 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.83719697e+04 / tc[1]
            -7.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
        /*species 13: HCO */
        species[13] =
            +3.83956496e+03 / tc[1]
            +8.26813410e-01
            -4.22118584e+00 * tc[0]
            +1.62196266e-03 * tc[1]
            -2.29665743e-06 * tc[2]
            +1.10953411e-09 * tc[3]
            -2.16884432e-13 * tc[4];
        /*species 14: CH2O */
        species[14] =
            -1.43089567e+04 / tc[1]
            +4.19091025e+00
            -4.79372315e+00 * tc[0]
            +4.95416684e-03 * tc[1]
            -6.22033347e-06 * tc[2]
            +3.16071051e-09 * tc[3]
            -6.58863260e-13 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +9.78601100e+02 / tc[1]
            -1.10459730e+01
            -2.10620400e+00 * tc[0]
            -3.60829750e-03 * tc[1]
            -8.89745333e-07 * tc[2]
            +6.14803000e-10 * tc[3]
            -1.03780500e-13 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +5.08977593e+03 / tc[1]
            -1.38129480e-01
            -3.95920148e+00 * tc[0]
            +3.78526124e-03 * tc[1]
            -9.51650487e-06 * tc[2]
            +5.76323961e-09 * tc[3]
            -1.34942187e-12 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +1.28416265e+04 / tc[1]
            -4.00743560e-01
            -4.30646568e+00 * tc[0]
            +2.09329446e-03 * tc[1]
            -8.28571345e-06 * tc[2]
            +4.99272172e-09 * tc[3]
            -1.15254502e-12 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            -1.15222055e+04 / tc[1]
            +1.62460176e+00
            -4.29142492e+00 * tc[0]
            +2.75077135e-03 * tc[1]
            -9.99063813e-06 * tc[2]
            +5.90388571e-09 * tc[3]
            -1.34342886e-12 * tc[4];
        /*species 19: N2 */
        species[19] =
            -1.02089990e+03 / tc[1]
            -6.51695000e-01
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
        /*species 20: AR */
        species[20] =
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
        /*species 7: CH2 */
        species[7] =
            +4.62636040e+04 / tc[1]
            -3.29709211e+00
            -2.87410113e+00 * tc[0]
            -1.82819646e-03 * tc[1]
            +2.34824328e-07 * tc[2]
            -2.16816291e-11 * tc[3]
            +9.38637835e-16 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +5.09259997e+04 / tc[1]
            -6.33446327e+00
            -2.29203842e+00 * tc[0]
            -2.32794318e-03 * tc[1]
            +3.35319912e-07 * tc[2]
            -3.48255000e-11 * tc[3]
            +1.69858182e-15 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +1.67755843e+04 / tc[1]
            -6.19435407e+00
            -2.28571772e+00 * tc[0]
            -3.61995018e-03 * tc[1]
            +4.97857247e-07 * tc[2]
            -4.96403870e-11 * tc[3]
            +2.33577197e-15 * tc[4];
        /*species 10: CH4 */
        species[10] =
            -9.46834459e+03 / tc[1]
            -1.83624665e+01
            -7.48514950e-02 * tc[0]
            -6.69547335e-03 * tc[1]
            +9.55476348e-07 * tc[2]
            -1.01910446e-10 * tc[3]
            +5.09076150e-15 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.41518724e+04 / tc[1]
            -5.10350211e+00
            -2.71518561e+00 * tc[0]
            -1.03126372e-03 * tc[1]
            +1.66470962e-07 * tc[2]
            -1.91710840e-11 * tc[3]
            +1.01823858e-15 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.87591660e+04 / tc[1]
            +1.58582223e+00
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
        /*species 13: HCO */
        species[13] =
            +4.01191815e+03 / tc[1]
            -7.02617054e+00
            -2.77217438e+00 * tc[0]
            -2.47847763e-03 * tc[1]
            +4.14076022e-07 * tc[2]
            -4.90968148e-11 * tc[3]
            +2.66754356e-15 * tc[4];
        /*species 14: CH2O */
        species[14] =
            -1.39958323e+04 / tc[1]
            -1.18956329e+01
            -1.76069008e+00 * tc[0]
            -4.60000041e-03 * tc[1]
            +7.37098022e-07 * tc[2]
            -8.38676767e-11 * tc[3]
            +4.41927820e-15 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +1.27832520e+02 / tc[1]
            +8.41224000e-01
            -3.77079900e+00 * tc[0]
            -3.93574850e-03 * tc[1]
            +4.42730667e-07 * tc[2]
            -3.28702583e-11 * tc[3]
            +1.05630800e-15 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +4.93988614e+03 / tc[1]
            -8.26925814e+00
            -2.03611116e+00 * tc[0]
            -7.32270755e-03 * tc[1]
            +1.11846319e-06 * tc[2]
            -1.22685769e-10 * tc[3]
            +6.28530305e-15 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +1.28575200e+04 / tc[1]
            -1.15077779e+01
            -1.95465642e+00 * tc[0]
            -8.69863610e-03 * tc[1]
            +1.33034445e-06 * tc[2]
            -1.46014741e-10 * tc[3]
            +7.48207880e-15 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            -1.14263932e+04 / tc[1]
            -1.40437292e+01
            -1.07188150e+00 * tc[0]
            -1.08426339e-02 * tc[1]
            +1.67093445e-06 * tc[2]
            -1.84510001e-10 * tc[3]
            +9.50014450e-15 * tc[4];
        /*species 19: N2 */
        species[19] =
            -9.22797700e+02 / tc[1]
            -3.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 20: AR */
        species[20] =
            -7.45375000e+02 / tc[1]
            -1.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
    }
    T_old = T;
    for (int i = 0; i < 21; i++)
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
        /*species 7: CH2 */
        species[7] =
            +4.60040401e+04 / tc[1]
            +1.20014682e+00
            -3.76267867e+00 * tc[0]
            -4.84436072e-04 * tc[1]
            -4.65816402e-07 * tc[2]
            +3.20909294e-10 * tc[3]
            -8.43708595e-14 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +5.04968163e+04 / tc[1]
            +3.96772308e+00
            -4.19860411e+00 * tc[0]
            +1.18330710e-03 * tc[1]
            -1.37216037e-06 * tc[2]
            +5.57346651e-10 * tc[3]
            -9.71573685e-14 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +1.64449988e+04 / tc[1]
            +1.06902607e+00
            -3.67359040e+00 * tc[0]
            -1.00547588e-03 * tc[1]
            -9.55036427e-07 * tc[2]
            +5.72597854e-10 * tc[3]
            -1.27192867e-13 * tc[4];
        /*species 10: CH4 */
        species[10] =
            -1.02466476e+04 / tc[1]
            +8.79117989e+00
            -5.14987613e+00 * tc[0]
            +6.83548940e-03 * tc[1]
            -8.19667665e-06 * tc[2]
            +4.03952522e-09 * tc[3]
            -8.33469780e-13 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.43440860e+04 / tc[1]
            -9.28875810e-01
            -3.57953347e+00 * tc[0]
            +3.05176840e-04 * tc[1]
            -1.69469055e-07 * tc[2]
            -7.55838237e-11 * tc[3]
            +4.52212249e-14 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.83719697e+04 / tc[1]
            -8.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
        /*species 13: HCO */
        species[13] =
            +3.83956496e+03 / tc[1]
            -1.73186590e-01
            -4.22118584e+00 * tc[0]
            +1.62196266e-03 * tc[1]
            -2.29665743e-06 * tc[2]
            +1.10953411e-09 * tc[3]
            -2.16884432e-13 * tc[4];
        /*species 14: CH2O */
        species[14] =
            -1.43089567e+04 / tc[1]
            +3.19091025e+00
            -4.79372315e+00 * tc[0]
            +4.95416684e-03 * tc[1]
            -6.22033347e-06 * tc[2]
            +3.16071051e-09 * tc[3]
            -6.58863260e-13 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +9.78601100e+02 / tc[1]
            -1.20459730e+01
            -2.10620400e+00 * tc[0]
            -3.60829750e-03 * tc[1]
            -8.89745333e-07 * tc[2]
            +6.14803000e-10 * tc[3]
            -1.03780500e-13 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +5.08977593e+03 / tc[1]
            -1.13812948e+00
            -3.95920148e+00 * tc[0]
            +3.78526124e-03 * tc[1]
            -9.51650487e-06 * tc[2]
            +5.76323961e-09 * tc[3]
            -1.34942187e-12 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +1.28416265e+04 / tc[1]
            -1.40074356e+00
            -4.30646568e+00 * tc[0]
            +2.09329446e-03 * tc[1]
            -8.28571345e-06 * tc[2]
            +4.99272172e-09 * tc[3]
            -1.15254502e-12 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            -1.15222055e+04 / tc[1]
            +6.24601760e-01
            -4.29142492e+00 * tc[0]
            +2.75077135e-03 * tc[1]
            -9.99063813e-06 * tc[2]
            +5.90388571e-09 * tc[3]
            -1.34342886e-12 * tc[4];
        /*species 19: N2 */
        species[19] =
            -1.02089990e+03 / tc[1]
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
        /*species 20: AR */
        species[20] =
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
        /*species 7: CH2 */
        species[7] =
            +4.62636040e+04 / tc[1]
            -4.29709211e+00
            -2.87410113e+00 * tc[0]
            -1.82819646e-03 * tc[1]
            +2.34824328e-07 * tc[2]
            -2.16816291e-11 * tc[3]
            +9.38637835e-16 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +5.09259997e+04 / tc[1]
            -7.33446327e+00
            -2.29203842e+00 * tc[0]
            -2.32794318e-03 * tc[1]
            +3.35319912e-07 * tc[2]
            -3.48255000e-11 * tc[3]
            +1.69858182e-15 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +1.67755843e+04 / tc[1]
            -7.19435407e+00
            -2.28571772e+00 * tc[0]
            -3.61995018e-03 * tc[1]
            +4.97857247e-07 * tc[2]
            -4.96403870e-11 * tc[3]
            +2.33577197e-15 * tc[4];
        /*species 10: CH4 */
        species[10] =
            -9.46834459e+03 / tc[1]
            -1.93624665e+01
            -7.48514950e-02 * tc[0]
            -6.69547335e-03 * tc[1]
            +9.55476348e-07 * tc[2]
            -1.01910446e-10 * tc[3]
            +5.09076150e-15 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.41518724e+04 / tc[1]
            -6.10350211e+00
            -2.71518561e+00 * tc[0]
            -1.03126372e-03 * tc[1]
            +1.66470962e-07 * tc[2]
            -1.91710840e-11 * tc[3]
            +1.01823858e-15 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.87591660e+04 / tc[1]
            +5.85822230e-01
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
        /*species 13: HCO */
        species[13] =
            +4.01191815e+03 / tc[1]
            -8.02617054e+00
            -2.77217438e+00 * tc[0]
            -2.47847763e-03 * tc[1]
            +4.14076022e-07 * tc[2]
            -4.90968148e-11 * tc[3]
            +2.66754356e-15 * tc[4];
        /*species 14: CH2O */
        species[14] =
            -1.39958323e+04 / tc[1]
            -1.28956329e+01
            -1.76069008e+00 * tc[0]
            -4.60000041e-03 * tc[1]
            +7.37098022e-07 * tc[2]
            -8.38676767e-11 * tc[3]
            +4.41927820e-15 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +1.27832520e+02 / tc[1]
            -1.58776000e-01
            -3.77079900e+00 * tc[0]
            -3.93574850e-03 * tc[1]
            +4.42730667e-07 * tc[2]
            -3.28702583e-11 * tc[3]
            +1.05630800e-15 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +4.93988614e+03 / tc[1]
            -9.26925814e+00
            -2.03611116e+00 * tc[0]
            -7.32270755e-03 * tc[1]
            +1.11846319e-06 * tc[2]
            -1.22685769e-10 * tc[3]
            +6.28530305e-15 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +1.28575200e+04 / tc[1]
            -1.25077779e+01
            -1.95465642e+00 * tc[0]
            -8.69863610e-03 * tc[1]
            +1.33034445e-06 * tc[2]
            -1.46014741e-10 * tc[3]
            +7.48207880e-15 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            -1.14263932e+04 / tc[1]
            -1.50437292e+01
            -1.07188150e+00 * tc[0]
            -1.08426339e-02 * tc[1]
            +1.67093445e-06 * tc[2]
            -1.84510001e-10 * tc[3]
            +9.50014450e-15 * tc[4];
        /*species 19: N2 */
        species[19] =
            -9.22797700e+02 / tc[1]
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 20: AR */
        species[20] =
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
        /*species 7: CH2 */
        species[7] =
            +2.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 / tc[1];
        /*species 8: CH2(S) */
        species[8] =
            +3.19860411e+00
            -1.18330710e-03 * tc[1]
            +2.74432073e-06 * tc[2]
            -1.67203995e-09 * tc[3]
            +3.88629474e-13 * tc[4]
            +5.04968163e+04 / tc[1];
        /*species 9: CH3 */
        species[9] =
            +2.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 / tc[1];
        /*species 10: CH4 */
        species[10] =
            +4.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 / tc[1];
        /*species 11: CO */
        species[11] =
            +2.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 / tc[1];
        /*species 12: CO2 */
        species[12] =
            +1.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 / tc[1];
        /*species 13: HCO */
        species[13] =
            +3.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 / tc[1];
        /*species 14: CH2O */
        species[14] =
            +3.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 / tc[1];
        /*species 15: CH3O */
        species[15] =
            +1.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122000e-13 * tc[4]
            +9.78601100e+02 / tc[1];
        /*species 16: C2H4 */
        species[16] =
            +2.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 / tc[1];
        /*species 17: C2H5 */
        species[17] =
            +3.30646568e+00
            -2.09329446e-03 * tc[1]
            +1.65714269e-05 * tc[2]
            -1.49781651e-08 * tc[3]
            +4.61018008e-12 * tc[4]
            +1.28416265e+04 / tc[1];
        /*species 18: C2H6 */
        species[18] =
            +3.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 / tc[1];
        /*species 19: N2 */
        species[19] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 / tc[1];
        /*species 20: AR */
        species[20] =
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
        /*species 7: CH2 */
        species[7] =
            +1.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 / tc[1];
        /*species 8: CH2(S) */
        species[8] =
            +1.29203842e+00
            +2.32794318e-03 * tc[1]
            -6.70639823e-07 * tc[2]
            +1.04476500e-10 * tc[3]
            -6.79432730e-15 * tc[4]
            +5.09259997e+04 / tc[1];
        /*species 9: CH3 */
        species[9] =
            +1.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 / tc[1];
        /*species 10: CH4 */
        species[10] =
            -9.25148505e-01
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 / tc[1];
        /*species 11: CO */
        species[11] =
            +1.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 / tc[1];
        /*species 12: CO2 */
        species[12] =
            +2.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 / tc[1];
        /*species 13: HCO */
        species[13] =
            +1.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 / tc[1];
        /*species 14: CH2O */
        species[14] =
            +7.60690080e-01
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 / tc[1];
        /*species 15: CH3O */
        species[15] =
            +2.77079900e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832520e+02 / tc[1];
        /*species 16: C2H4 */
        species[16] =
            +1.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 / tc[1];
        /*species 17: C2H5 */
        species[17] =
            +9.54656420e-01
            +8.69863610e-03 * tc[1]
            -2.66068889e-06 * tc[2]
            +4.38044223e-10 * tc[3]
            -2.99283152e-14 * tc[4]
            +1.28575200e+04 / tc[1];
        /*species 18: C2H6 */
        species[18] =
            +7.18815000e-02
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 / tc[1];
        /*species 19: N2 */
        species[19] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 / tc[1];
        /*species 20: AR */
        species[20] =
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
        /*species 7: CH2 */
        species[7] =
            +3.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 / tc[1];
        /*species 8: CH2(S) */
        species[8] =
            +4.19860411e+00
            -1.18330710e-03 * tc[1]
            +2.74432073e-06 * tc[2]
            -1.67203995e-09 * tc[3]
            +3.88629474e-13 * tc[4]
            +5.04968163e+04 / tc[1];
        /*species 9: CH3 */
        species[9] =
            +3.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 / tc[1];
        /*species 10: CH4 */
        species[10] =
            +5.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 / tc[1];
        /*species 11: CO */
        species[11] =
            +3.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 / tc[1];
        /*species 12: CO2 */
        species[12] =
            +2.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 / tc[1];
        /*species 13: HCO */
        species[13] =
            +4.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 / tc[1];
        /*species 14: CH2O */
        species[14] =
            +4.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 / tc[1];
        /*species 15: CH3O */
        species[15] =
            +2.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122000e-13 * tc[4]
            +9.78601100e+02 / tc[1];
        /*species 16: C2H4 */
        species[16] =
            +3.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 / tc[1];
        /*species 17: C2H5 */
        species[17] =
            +4.30646568e+00
            -2.09329446e-03 * tc[1]
            +1.65714269e-05 * tc[2]
            -1.49781651e-08 * tc[3]
            +4.61018008e-12 * tc[4]
            +1.28416265e+04 / tc[1];
        /*species 18: C2H6 */
        species[18] =
            +4.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 / tc[1];
        /*species 19: N2 */
        species[19] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 / tc[1];
        /*species 20: AR */
        species[20] =
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
        /*species 7: CH2 */
        species[7] =
            +2.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 / tc[1];
        /*species 8: CH2(S) */
        species[8] =
            +2.29203842e+00
            +2.32794318e-03 * tc[1]
            -6.70639823e-07 * tc[2]
            +1.04476500e-10 * tc[3]
            -6.79432730e-15 * tc[4]
            +5.09259997e+04 / tc[1];
        /*species 9: CH3 */
        species[9] =
            +2.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 / tc[1];
        /*species 10: CH4 */
        species[10] =
            +7.48514950e-02
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 / tc[1];
        /*species 11: CO */
        species[11] =
            +2.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 / tc[1];
        /*species 12: CO2 */
        species[12] =
            +3.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 / tc[1];
        /*species 13: HCO */
        species[13] =
            +2.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 / tc[1];
        /*species 14: CH2O */
        species[14] =
            +1.76069008e+00
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 / tc[1];
        /*species 15: CH3O */
        species[15] =
            +3.77079900e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832520e+02 / tc[1];
        /*species 16: C2H4 */
        species[16] =
            +2.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 / tc[1];
        /*species 17: C2H5 */
        species[17] =
            +1.95465642e+00
            +8.69863610e-03 * tc[1]
            -2.66068889e-06 * tc[2]
            +4.38044223e-10 * tc[3]
            -2.99283152e-14 * tc[4]
            +1.28575200e+04 / tc[1];
        /*species 18: C2H6 */
        species[18] =
            +1.07188150e+00
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 / tc[1];
        /*species 19: N2 */
        species[19] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 / tc[1];
        /*species 20: AR */
        species[20] =
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
    wt[7] = 14.026940; /*CH2 */
    wt[8] = 14.026940; /*CH2(S) */
    wt[9] = 15.034910; /*CH3 */
    wt[10] = 16.042880; /*CH4 */
    wt[11] = 28.010400; /*CO */
    wt[12] = 44.009800; /*CO2 */
    wt[13] = 29.018370; /*HCO */
    wt[14] = 30.026340; /*CH2O */
    wt[15] = 31.034310; /*CH3O */
    wt[16] = 28.053880; /*C2H4 */
    wt[17] = 29.061850; /*C2H5 */
    wt[18] = 30.069820; /*C2H6 */
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
    y[7] = phi[7]*14.026940;   XW += y[7]; /*CH2 */
    y[8] = phi[8]*14.026940;   XW += y[8]; /*CH2(S) */
    y[9] = phi[9]*15.034910;   XW += y[9]; /*CH3 */
    y[10] = phi[10]*16.042880;   XW += y[10]; /*CH4 */
    y[11] = phi[11]*28.010400;   XW += y[11]; /*CO */
    y[12] = phi[12]*44.009800;   XW += y[12]; /*CO2 */
    y[13] = phi[13]*29.018370;   XW += y[13]; /*HCO */
    y[14] = phi[14]*30.026340;   XW += y[14]; /*CH2O */
    y[15] = phi[15]*31.034310;   XW += y[15]; /*CH3O */
    y[16] = phi[16]*28.053880;   XW += y[16]; /*C2H4 */
    y[17] = phi[17]*29.061850;   XW += y[17]; /*C2H5 */
    y[18] = phi[18]*30.069820;   XW += y[18]; /*C2H6 */
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
    phi[7] = y[7]/ 1.40269400e-02; /*CH2 (wt in kg) */
    phi[8] = y[8]/ 1.40269400e-02; /*CH2(S) (wt in kg) */
    phi[9] = y[9]/ 1.50349100e-02; /*CH3 (wt in kg) */
    phi[10] = y[10]/ 1.60428800e-02; /*CH4 (wt in kg) */
    phi[11] = y[11]/ 2.80104000e-02; /*CO (wt in kg) */
    phi[12] = y[12]/ 4.40098000e-02; /*CO2 (wt in kg) */
    phi[13] = y[13]/ 2.90183700e-02; /*HCO (wt in kg) */
    phi[14] = y[14]/ 3.00263400e-02; /*CH2O (wt in kg) */
    phi[15] = y[15]/ 3.10343100e-02; /*CH3O (wt in kg) */
    phi[16] = y[16]/ 2.80538800e-02; /*C2H4 (wt in kg) */
    phi[17] = y[17]/ 2.90618500e-02; /*C2H5 (wt in kg) */
    phi[18] = y[18]/ 3.00698200e-02; /*C2H6 (wt in kg) */
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
    y[7] = c[7] * 14.026940 / (*rho); 
    y[8] = c[8] * 14.026940 / (*rho); 
    y[9] = c[9] * 15.034910 / (*rho); 
    y[10] = c[10] * 16.042880 / (*rho); 
    y[11] = c[11] * 28.010400 / (*rho); 
    y[12] = c[12] * 44.009800 / (*rho); 
    y[13] = c[13] * 29.018370 / (*rho); 
    y[14] = c[14] * 30.026340 / (*rho); 
    y[15] = c[15] * 31.034310 / (*rho); 
    y[16] = c[16] * 28.053880 / (*rho); 
    y[17] = c[17] * 29.061850 / (*rho); 
    y[18] = c[18] * 30.069820 / (*rho); 
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
