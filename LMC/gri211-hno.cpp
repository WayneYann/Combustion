
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
void CKSYME(char * cckwrk, int * lout, char * kname, int * kerr, int lencck, int lenkname);
void CKSYMS(char * cckwrk, int * lout, char * kname, int * kerr, int lencck, int lenkname);
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
    *mm = 4;
    *kk = 18;
    *ii = 68;
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
void CKSYME(char * cckwrk, int * lout, char * kname, int * kerr, int lencck, int lenkname )
{
    int i; /*Loop Counter */
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
void CKSYMS(char * cckwrk, int * lout, char * kname, int * kerr, int lencck, int lenkname )
{
    int i; /*Loop Counter */
    /*clear kname */
    for (i=0; i<lenkname*18; i++) {
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

    /* N  */
    kname[ 8*lenkname + 0 ] = 'N';
    kname[ 8*lenkname + 1 ] = ' ';

    /* NH  */
    kname[ 9*lenkname + 0 ] = 'N';
    kname[ 9*lenkname + 1 ] = 'H';
    kname[ 9*lenkname + 2 ] = ' ';

    /* NH2  */
    kname[ 10*lenkname + 0 ] = 'N';
    kname[ 10*lenkname + 1 ] = 'H';
    kname[ 10*lenkname + 2 ] = '2';
    kname[ 10*lenkname + 3 ] = ' ';

    /* NH3  */
    kname[ 11*lenkname + 0 ] = 'N';
    kname[ 11*lenkname + 1 ] = 'H';
    kname[ 11*lenkname + 2 ] = '3';
    kname[ 11*lenkname + 3 ] = ' ';

    /* NNH  */
    kname[ 12*lenkname + 0 ] = 'N';
    kname[ 12*lenkname + 1 ] = 'N';
    kname[ 12*lenkname + 2 ] = 'H';
    kname[ 12*lenkname + 3 ] = ' ';

    /* NO  */
    kname[ 13*lenkname + 0 ] = 'N';
    kname[ 13*lenkname + 1 ] = 'O';
    kname[ 13*lenkname + 2 ] = ' ';

    /* NO2  */
    kname[ 14*lenkname + 0 ] = 'N';
    kname[ 14*lenkname + 1 ] = 'O';
    kname[ 14*lenkname + 2 ] = '2';
    kname[ 14*lenkname + 3 ] = ' ';

    /* N2O  */
    kname[ 15*lenkname + 0 ] = 'N';
    kname[ 15*lenkname + 1 ] = '2';
    kname[ 15*lenkname + 2 ] = 'O';
    kname[ 15*lenkname + 3 ] = ' ';

    /* HNO  */
    kname[ 16*lenkname + 0 ] = 'H';
    kname[ 16*lenkname + 1 ] = 'N';
    kname[ 16*lenkname + 2 ] = 'O';
    kname[ 16*lenkname + 3 ] = ' ';

    /* N2  */
    kname[ 17*lenkname + 0 ] = 'N';
    kname[ 17*lenkname + 1 ] = '2';
    kname[ 17*lenkname + 2 ] = ' ';

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
    XW += x[8]*14.006700; /*N */
    XW += x[9]*15.014670; /*NH */
    XW += x[10]*16.022640; /*NH2 */
    XW += x[11]*17.030610; /*NH3 */
    XW += x[12]*29.021370; /*NNH */
    XW += x[13]*30.006100; /*NO */
    XW += x[14]*46.005500; /*NO2 */
    XW += x[15]*44.012800; /*N2O */
    XW += x[16]*31.014070; /*HNO */
    XW += x[17]*28.013400; /*N2 */
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
    YOW += y[8]/14.006700; /*N */
    YOW += y[9]/15.014670; /*NH */
    YOW += y[10]/16.022640; /*NH2 */
    YOW += y[11]/17.030610; /*NH3 */
    YOW += y[12]/29.021370; /*NNH */
    YOW += y[13]/30.006100; /*NO */
    YOW += y[14]/46.005500; /*NO2 */
    YOW += y[15]/44.012800; /*N2O */
    YOW += y[16]/31.014070; /*HNO */
    YOW += y[17]/28.013400; /*N2 */
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
    W += c[8]*14.006700; /*N */
    W += c[9]*15.014670; /*NH */
    W += c[10]*16.022640; /*NH2 */
    W += c[11]*17.030610; /*NH3 */
    W += c[12]*29.021370; /*NNH */
    W += c[13]*30.006100; /*NO */
    W += c[14]*46.005500; /*NO2 */
    W += c[15]*44.012800; /*N2O */
    W += c[16]*31.014070; /*HNO */
    W += c[17]*28.013400; /*N2 */

    for (id = 0; id < 18; ++id) {
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
    XW += x[8]*14.006700; /*N */
    XW += x[9]*15.014670; /*NH */
    XW += x[10]*16.022640; /*NH2 */
    XW += x[11]*17.030610; /*NH3 */
    XW += x[12]*29.021370; /*NNH */
    XW += x[13]*30.006100; /*NO */
    XW += x[14]*46.005500; /*NO2 */
    XW += x[15]*44.012800; /*N2O */
    XW += x[16]*31.014070; /*HNO */
    XW += x[17]*28.013400; /*N2 */
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
    YOW += y[8]/14.006700; /*N */
    YOW += y[9]/15.014670; /*NH */
    YOW += y[10]/16.022640; /*NH2 */
    YOW += y[11]/17.030610; /*NH3 */
    YOW += y[12]/29.021370; /*NNH */
    YOW += y[13]/30.006100; /*NO */
    YOW += y[14]/46.005500; /*NO2 */
    YOW += y[15]/44.012800; /*N2O */
    YOW += y[16]/31.014070; /*HNO */
    YOW += y[17]/28.013400; /*N2 */
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
    W += c[8]*14.006700; /*N */
    W += c[9]*15.014670; /*NH */
    W += c[10]*16.022640; /*NH2 */
    W += c[11]*17.030610; /*NH3 */
    W += c[12]*29.021370; /*NNH */
    W += c[13]*30.006100; /*NO */
    W += c[14]*46.005500; /*NO2 */
    W += c[15]*44.012800; /*N2O */
    W += c[16]*31.014070; /*HNO */
    W += c[17]*28.013400; /*N2 */

    for (id = 0; id < 18; ++id) {
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
    YOW += y[8]/14.006700; /*N */
    YOW += y[9]/15.014670; /*NH */
    YOW += y[10]/16.022640; /*NH2 */
    YOW += y[11]/17.030610; /*NH3 */
    YOW += y[12]/29.021370; /*NNH */
    YOW += y[13]/30.006100; /*NO */
    YOW += y[14]/46.005500; /*NO2 */
    YOW += y[15]/44.012800; /*N2O */
    YOW += y[16]/31.014070; /*HNO */
    YOW += y[17]/28.013400; /*N2 */
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
    XW += x[8]*14.006700; /*N */
    XW += x[9]*15.014670; /*NH */
    XW += x[10]*16.022640; /*NH2 */
    XW += x[11]*17.030610; /*NH3 */
    XW += x[12]*29.021370; /*NNH */
    XW += x[13]*30.006100; /*NO */
    XW += x[14]*46.005500; /*NO2 */
    XW += x[15]*44.012800; /*N2O */
    XW += x[16]*31.014070; /*HNO */
    XW += x[17]*28.013400; /*N2 */
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
    W += c[8]*14.006700; /*N */
    W += c[9]*15.014670; /*NH */
    W += c[10]*16.022640; /*NH2 */
    W += c[11]*17.030610; /*NH3 */
    W += c[12]*29.021370; /*NNH */
    W += c[13]*30.006100; /*NO */
    W += c[14]*46.005500; /*NO2 */
    W += c[15]*44.012800; /*N2O */
    W += c[16]*31.014070; /*HNO */
    W += c[17]*28.013400; /*N2 */

    for (id = 0; id < 18; ++id) {
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
    YOW += y[8]/14.006700; /*N */
    YOW += y[9]/15.014670; /*NH */
    YOW += y[10]/16.022640; /*NH2 */
    YOW += y[11]/17.030610; /*NH3 */
    YOW += y[12]/29.021370; /*NNH */
    YOW += y[13]/30.006100; /*NO */
    YOW += y[14]/46.005500; /*NO2 */
    YOW += y[15]/44.012800; /*N2O */
    YOW += y[16]/31.014070; /*HNO */
    YOW += y[17]/28.013400; /*N2 */
    /*Now compute conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(14.006700*YOW); 
    x[9] = y[9]/(15.014670*YOW); 
    x[10] = y[10]/(16.022640*YOW); 
    x[11] = y[11]/(17.030610*YOW); 
    x[12] = y[12]/(29.021370*YOW); 
    x[13] = y[13]/(30.006100*YOW); 
    x[14] = y[14]/(46.005500*YOW); 
    x[15] = y[15]/(44.012800*YOW); 
    x[16] = y[16]/(31.014070*YOW); 
    x[17] = y[17]/(28.013400*YOW); 

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
    YOW += y[8]/14.006700; /*N */
    YOW += y[9]/15.014670; /*NH */
    YOW += y[10]/16.022640; /*NH2 */
    YOW += y[11]/17.030610; /*NH3 */
    YOW += y[12]/29.021370; /*NNH */
    YOW += y[13]/30.006100; /*NO */
    YOW += y[14]/46.005500; /*NO2 */
    YOW += y[15]/44.012800; /*N2O */
    YOW += y[16]/31.014070; /*HNO */
    YOW += y[17]/28.013400; /*N2 */
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
    c[8] = PWORT * y[8]/14.006700; 
    c[9] = PWORT * y[9]/15.014670; 
    c[10] = PWORT * y[10]/16.022640; 
    c[11] = PWORT * y[11]/17.030610; 
    c[12] = PWORT * y[12]/29.021370; 
    c[13] = PWORT * y[13]/30.006100; 
    c[14] = PWORT * y[14]/46.005500; 
    c[15] = PWORT * y[15]/44.012800; 
    c[16] = PWORT * y[16]/31.014070; 
    c[17] = PWORT * y[17]/28.013400; 

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
    c[8] = (*rho) * y[8]/14.006700; 
    c[9] = (*rho) * y[9]/15.014670; 
    c[10] = (*rho) * y[10]/16.022640; 
    c[11] = (*rho) * y[11]/17.030610; 
    c[12] = (*rho) * y[12]/29.021370; 
    c[13] = (*rho) * y[13]/30.006100; 
    c[14] = (*rho) * y[14]/46.005500; 
    c[15] = (*rho) * y[15]/44.012800; 
    c[16] = (*rho) * y[16]/31.014070; 
    c[17] = (*rho) * y[17]/28.013400; 

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
    XW += x[8]*14.006700; /*N */
    XW += x[9]*15.014670; /*NH */
    XW += x[10]*16.022640; /*NH2 */
    XW += x[11]*17.030610; /*NH3 */
    XW += x[12]*29.021370; /*NNH */
    XW += x[13]*30.006100; /*NO */
    XW += x[14]*46.005500; /*NO2 */
    XW += x[15]*44.012800; /*N2O */
    XW += x[16]*31.014070; /*HNO */
    XW += x[17]*28.013400; /*N2 */
    /*Now compute conversion */
    y[0] = x[0]*2.015940/XW; 
    y[1] = x[1]*1.007970/XW; 
    y[2] = x[2]*15.999400/XW; 
    y[3] = x[3]*31.998800/XW; 
    y[4] = x[4]*17.007370/XW; 
    y[5] = x[5]*18.015340/XW; 
    y[6] = x[6]*33.006770/XW; 
    y[7] = x[7]*34.014740/XW; 
    y[8] = x[8]*14.006700/XW; 
    y[9] = x[9]*15.014670/XW; 
    y[10] = x[10]*16.022640/XW; 
    y[11] = x[11]*17.030610/XW; 
    y[12] = x[12]*29.021370/XW; 
    y[13] = x[13]*30.006100/XW; 
    y[14] = x[14]*46.005500/XW; 
    y[15] = x[15]*44.012800/XW; 
    y[16] = x[16]*31.014070/XW; 
    y[17] = x[17]*28.013400/XW; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.314e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 18; ++id) {
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
    XW += x[8]*14.006700; /*N */
    XW += x[9]*15.014670; /*NH */
    XW += x[10]*16.022640; /*NH2 */
    XW += x[11]*17.030610; /*NH3 */
    XW += x[12]*29.021370; /*NNH */
    XW += x[13]*30.006100; /*NO */
    XW += x[14]*46.005500; /*NO2 */
    XW += x[15]*44.012800; /*N2O */
    XW += x[16]*31.014070; /*HNO */
    XW += x[17]*28.013400; /*N2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 18; ++id) {
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
    for (id = 0; id < 18; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    for (id = 0; id < 18; ++id) {
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
    CW += c[8]*14.006700; /*N */
    CW += c[9]*15.014670; /*NH */
    CW += c[10]*16.022640; /*NH2 */
    CW += c[11]*17.030610; /*NH3 */
    CW += c[12]*29.021370; /*NNH */
    CW += c[13]*30.006100; /*NO */
    CW += c[14]*46.005500; /*NO2 */
    CW += c[15]*44.012800; /*N2O */
    CW += c[16]*31.014070; /*HNO */
    CW += c[17]*28.013400; /*N2 */
    /*Now compute conversion */
    y[0] = c[0]*2.015940/CW; 
    y[1] = c[1]*1.007970/CW; 
    y[2] = c[2]*15.999400/CW; 
    y[3] = c[3]*31.998800/CW; 
    y[4] = c[4]*17.007370/CW; 
    y[5] = c[5]*18.015340/CW; 
    y[6] = c[6]*33.006770/CW; 
    y[7] = c[7]*34.014740/CW; 
    y[8] = c[8]*14.006700/CW; 
    y[9] = c[9]*15.014670/CW; 
    y[10] = c[10]*16.022640/CW; 
    y[11] = c[11]*17.030610/CW; 
    y[12] = c[12]*29.021370/CW; 
    y[13] = c[13]*30.006100/CW; 
    y[14] = c[14]*46.005500/CW; 
    y[15] = c[15]*44.012800/CW; 
    y[16] = c[16]*31.014070/CW; 
    y[17] = c[17]*28.013400/CW; 

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
    for (id = 0; id < 18; ++id) {
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
    for (id = 0; id < 18; ++id) {
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
    for (id = 0; id < 18; ++id) {
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
    for (id = 0; id < 18; ++id) {
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
    for (id = 0; id < 18; ++id) {
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
    for (id = 0; id < 18; ++id) {
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
    for (id = 0; id < 18; ++id) {
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
    cvms[8] *= 5935730.757423; /*N */
    cvms[9] *= 5537251.234959; /*NH */
    cvms[10] *= 5188907.695611; /*NH2 */
    cvms[11] *= 4881798.127020; /*NH3 */
    cvms[12] *= 2864785.501167; /*NNH */
    cvms[13] *= 2770769.943445; /*NO */
    cvms[14] *= 1807175.229049; /*NO2 */
    cvms[15] *= 1888995.928457; /*N2O */
    cvms[16] *= 2680718.783442; /*HNO */
    cvms[17] *= 2967865.378712; /*N2 */
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
    cpms[8] *= 5935730.757423; /*N */
    cpms[9] *= 5537251.234959; /*NH */
    cpms[10] *= 5188907.695611; /*NH2 */
    cpms[11] *= 4881798.127020; /*NH3 */
    cpms[12] *= 2864785.501167; /*NNH */
    cpms[13] *= 2770769.943445; /*NO */
    cpms[14] *= 1807175.229049; /*NO2 */
    cpms[15] *= 1888995.928457; /*N2O */
    cpms[16] *= 2680718.783442; /*HNO */
    cpms[17] *= 2967865.378712; /*N2 */
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
    ums[8] *= RT/14.006700; /*N */
    ums[9] *= RT/15.014670; /*NH */
    ums[10] *= RT/16.022640; /*NH2 */
    ums[11] *= RT/17.030610; /*NH3 */
    ums[12] *= RT/29.021370; /*NNH */
    ums[13] *= RT/30.006100; /*NO */
    ums[14] *= RT/46.005500; /*NO2 */
    ums[15] *= RT/44.012800; /*N2O */
    ums[16] *= RT/31.014070; /*HNO */
    ums[17] *= RT/28.013400; /*N2 */
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
    hms[8] *= RT/14.006700; /*N */
    hms[9] *= RT/15.014670; /*NH */
    hms[10] *= RT/16.022640; /*NH2 */
    hms[11] *= RT/17.030610; /*NH3 */
    hms[12] *= RT/29.021370; /*NNH */
    hms[13] *= RT/30.006100; /*NO */
    hms[14] *= RT/46.005500; /*NO2 */
    hms[15] *= RT/44.012800; /*N2O */
    hms[16] *= RT/31.014070; /*HNO */
    hms[17] *= RT/28.013400; /*N2 */
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
    gms[8] *= RT/14.006700; /*N */
    gms[9] *= RT/15.014670; /*NH */
    gms[10] *= RT/16.022640; /*NH2 */
    gms[11] *= RT/17.030610; /*NH3 */
    gms[12] *= RT/29.021370; /*NNH */
    gms[13] *= RT/30.006100; /*NO */
    gms[14] *= RT/46.005500; /*NO2 */
    gms[15] *= RT/44.012800; /*N2O */
    gms[16] *= RT/31.014070; /*HNO */
    gms[17] *= RT/28.013400; /*N2 */
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
    ams[8] *= RT/14.006700; /*N */
    ams[9] *= RT/15.014670; /*NH */
    ams[10] *= RT/16.022640; /*NH2 */
    ams[11] *= RT/17.030610; /*NH3 */
    ams[12] *= RT/29.021370; /*NNH */
    ams[13] *= RT/30.006100; /*NO */
    ams[14] *= RT/46.005500; /*NO2 */
    ams[15] *= RT/44.012800; /*N2O */
    ams[16] *= RT/31.014070; /*HNO */
    ams[17] *= RT/28.013400; /*N2 */
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
    sms[8] *= 5935730.757423; /*N */
    sms[9] *= 5537251.234959; /*NH */
    sms[10] *= 5188907.695611; /*NH2 */
    sms[11] *= 4881798.127020; /*NH3 */
    sms[12] *= 2864785.501167; /*NNH */
    sms[13] *= 2770769.943445; /*NO */
    sms[14] *= 1807175.229049; /*NO2 */
    sms[15] *= 1888995.928457; /*N2O */
    sms[16] *= 2680718.783442; /*HNO */
    sms[17] *= 2967865.378712; /*N2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double *T, double *x, int * iwrk, double * rwrk, double * cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[18]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 18; ++id) {
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
    double cpor[18]; /* temporary storage */
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
    result += cpor[8]*y[8]/14.0067; /*N */
    result += cpor[9]*y[9]/15.0147; /*NH */
    result += cpor[10]*y[10]/16.0226; /*NH2 */
    result += cpor[11]*y[11]/17.0306; /*NH3 */
    result += cpor[12]*y[12]/29.0214; /*NNH */
    result += cpor[13]*y[13]/30.0061; /*NO */
    result += cpor[14]*y[14]/46.0055; /*NO2 */
    result += cpor[15]*y[15]/44.0128; /*N2O */
    result += cpor[16]*y[16]/31.0141; /*HNO */
    result += cpor[17]*y[17]/28.0134; /*N2 */

    *cpbs = result * 8.314e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double *T, double *x, int * iwrk, double * rwrk, double * cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[18]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 18; ++id) {
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
    double cvor[18]; /* temporary storage */
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
    result += cvor[8]*y[8]/14.0067; /*N */
    result += cvor[9]*y[9]/15.0147; /*NH */
    result += cvor[10]*y[10]/16.0226; /*NH2 */
    result += cvor[11]*y[11]/17.0306; /*NH3 */
    result += cvor[12]*y[12]/29.0214; /*NNH */
    result += cvor[13]*y[13]/30.0061; /*NO */
    result += cvor[14]*y[14]/46.0055; /*NO2 */
    result += cvor[15]*y[15]/44.0128; /*N2O */
    result += cvor[16]*y[16]/31.0141; /*HNO */
    result += cvor[17]*y[17]/28.0134; /*N2 */

    *cvbs = result * 8.314e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double *T, double *x, int * iwrk, double * rwrk, double * hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[18]; /* temporary storage */
    double RT = 8.314e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 18; ++id) {
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
    double hml[18]; /* temporary storage */
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
    result += y[8]*hml[8]/14.006700; /*N */
    result += y[9]*hml[9]/15.014670; /*NH */
    result += y[10]*hml[10]/16.022640; /*NH2 */
    result += y[11]*hml[11]/17.030610; /*NH3 */
    result += y[12]*hml[12]/29.021370; /*NNH */
    result += y[13]*hml[13]/30.006100; /*NO */
    result += y[14]*hml[14]/46.005500; /*NO2 */
    result += y[15]*hml[15]/44.012800; /*N2O */
    result += y[16]*hml[16]/31.014070; /*HNO */
    result += y[17]*hml[17]/28.013400; /*N2 */

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double *T, double *x, int * iwrk, double * rwrk, double * ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[18]; /* temporary energy array */
    double RT = 8.314e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 18; ++id) {
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
    double ums[18]; /* temporary energy array */
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
    result += y[8]*ums[8]/14.006700; /*N */
    result += y[9]*ums[9]/15.014670; /*NH */
    result += y[10]*ums[10]/16.022640; /*NH2 */
    result += y[11]*ums[11]/17.030610; /*NH3 */
    result += y[12]*ums[12]/29.021370; /*NNH */
    result += y[13]*ums[13]/30.006100; /*NO */
    result += y[14]*ums[14]/46.005500; /*NO2 */
    result += y[15]*ums[15]/44.012800; /*N2O */
    result += y[16]*ums[16]/31.014070; /*HNO */
    result += y[17]*ums[17]/28.013400; /*N2 */

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
    double sor[18]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 18; ++id) {
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
    double sor[18]; /* temporary storage */
    double x[18]; /* need a ytx conversion */
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
    YOW += y[8]/14.006700; /*N */
    YOW += y[9]/15.014670; /*NH */
    YOW += y[10]/16.022640; /*NH2 */
    YOW += y[11]/17.030610; /*NH3 */
    YOW += y[12]/29.021370; /*NNH */
    YOW += y[13]/30.006100; /*NO */
    YOW += y[14]/46.005500; /*NO2 */
    YOW += y[15]/44.012800; /*N2O */
    YOW += y[16]/31.014070; /*HNO */
    YOW += y[17]/28.013400; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(14.006700*YOW); 
    x[9] = y[9]/(15.014670*YOW); 
    x[10] = y[10]/(16.022640*YOW); 
    x[11] = y[11]/(17.030610*YOW); 
    x[12] = y[12]/(29.021370*YOW); 
    x[13] = y[13]/(30.006100*YOW); 
    x[14] = y[14]/(46.005500*YOW); 
    x[15] = y[15]/(44.012800*YOW); 
    x[16] = y[16]/(31.014070*YOW); 
    x[17] = y[17]/(28.013400*YOW); 
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
    double gort[18]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 18; ++id) {
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
    double gort[18]; /* temporary storage */
    double x[18]; /* need a ytx conversion */
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
    YOW += y[8]/14.006700; /*N */
    YOW += y[9]/15.014670; /*NH */
    YOW += y[10]/16.022640; /*NH2 */
    YOW += y[11]/17.030610; /*NH3 */
    YOW += y[12]/29.021370; /*NNH */
    YOW += y[13]/30.006100; /*NO */
    YOW += y[14]/46.005500; /*NO2 */
    YOW += y[15]/44.012800; /*N2O */
    YOW += y[16]/31.014070; /*HNO */
    YOW += y[17]/28.013400; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(14.006700*YOW); 
    x[9] = y[9]/(15.014670*YOW); 
    x[10] = y[10]/(16.022640*YOW); 
    x[11] = y[11]/(17.030610*YOW); 
    x[12] = y[12]/(29.021370*YOW); 
    x[13] = y[13]/(30.006100*YOW); 
    x[14] = y[14]/(46.005500*YOW); 
    x[15] = y[15]/(44.012800*YOW); 
    x[16] = y[16]/(31.014070*YOW); 
    x[17] = y[17]/(28.013400*YOW); 
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
    double aort[18]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 18; ++id) {
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
    double aort[18]; /* temporary storage */
    double x[18]; /* need a ytx conversion */
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
    YOW += y[8]/14.006700; /*N */
    YOW += y[9]/15.014670; /*NH */
    YOW += y[10]/16.022640; /*NH2 */
    YOW += y[11]/17.030610; /*NH3 */
    YOW += y[12]/29.021370; /*NNH */
    YOW += y[13]/30.006100; /*NO */
    YOW += y[14]/46.005500; /*NO2 */
    YOW += y[15]/44.012800; /*N2O */
    YOW += y[16]/31.014070; /*HNO */
    YOW += y[17]/28.013400; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(14.006700*YOW); 
    x[9] = y[9]/(15.014670*YOW); 
    x[10] = y[10]/(16.022640*YOW); 
    x[11] = y[11]/(17.030610*YOW); 
    x[12] = y[12]/(29.021370*YOW); 
    x[13] = y[13]/(30.006100*YOW); 
    x[14] = y[14]/(46.005500*YOW); 
    x[15] = y[15]/(44.012800*YOW); 
    x[16] = y[16]/(31.014070*YOW); 
    x[17] = y[17]/(28.013400*YOW); 
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
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double * T, double * C, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 18; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 18; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[18]; /*temporary storage */
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
    YOW += y[8]/14.006700; /*N */
    YOW += y[9]/15.014670; /*NH */
    YOW += y[10]/16.022640; /*NH2 */
    YOW += y[11]/17.030610; /*NH3 */
    YOW += y[12]/29.021370; /*NNH */
    YOW += y[13]/30.006100; /*NO */
    YOW += y[14]/46.005500; /*NO2 */
    YOW += y[15]/44.012800; /*N2O */
    YOW += y[16]/31.014070; /*HNO */
    YOW += y[17]/28.013400; /*N2 */
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
    c[8] = PWORT * y[8]/14.006700; 
    c[9] = PWORT * y[9]/15.014670; 
    c[10] = PWORT * y[10]/16.022640; 
    c[11] = PWORT * y[11]/17.030610; 
    c[12] = PWORT * y[12]/29.021370; 
    c[13] = PWORT * y[13]/30.006100; 
    c[14] = PWORT * y[14]/46.005500; 
    c[15] = PWORT * y[15]/44.012800; 
    c[16] = PWORT * y[16]/31.014070; 
    c[17] = PWORT * y[17]/28.013400; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 18; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[18]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.314e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 18; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 18; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[18]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]/2.015940; 
    c[1] = 1e6 * (*rho) * y[1]/1.007970; 
    c[2] = 1e6 * (*rho) * y[2]/15.999400; 
    c[3] = 1e6 * (*rho) * y[3]/31.998800; 
    c[4] = 1e6 * (*rho) * y[4]/17.007370; 
    c[5] = 1e6 * (*rho) * y[5]/18.015340; 
    c[6] = 1e6 * (*rho) * y[6]/33.006770; 
    c[7] = 1e6 * (*rho) * y[7]/34.014740; 
    c[8] = 1e6 * (*rho) * y[8]/14.006700; 
    c[9] = 1e6 * (*rho) * y[9]/15.014670; 
    c[10] = 1e6 * (*rho) * y[10]/16.022640; 
    c[11] = 1e6 * (*rho) * y[11]/17.030610; 
    c[12] = 1e6 * (*rho) * y[12]/29.021370; 
    c[13] = 1e6 * (*rho) * y[13]/30.006100; 
    c[14] = 1e6 * (*rho) * y[14]/46.005500; 
    c[15] = 1e6 * (*rho) * y[15]/44.012800; 
    c[16] = 1e6 * (*rho) * y[16]/31.014070; 
    c[17] = 1e6 * (*rho) * y[17]/28.013400; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 18; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[18]; /*temporary storage */
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
    XW += x[8]*14.006700; /*N */
    XW += x[9]*15.014670; /*NH */
    XW += x[10]*16.022640; /*NH2 */
    XW += x[11]*17.030610; /*NH3 */
    XW += x[12]*29.021370; /*NNH */
    XW += x[13]*30.006100; /*NO */
    XW += x[14]*46.005500; /*NO2 */
    XW += x[15]*44.012800; /*N2O */
    XW += x[16]*31.014070; /*HNO */
    XW += x[17]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 18; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 18; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double * T, double * C, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 18; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 18; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 68; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[18]; /*temporary storage */
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
    YOW += y[8]/14.006700; /*N */
    YOW += y[9]/15.014670; /*NH */
    YOW += y[10]/16.022640; /*NH2 */
    YOW += y[11]/17.030610; /*NH3 */
    YOW += y[12]/29.021370; /*NNH */
    YOW += y[13]/30.006100; /*NO */
    YOW += y[14]/46.005500; /*NO2 */
    YOW += y[15]/44.012800; /*N2O */
    YOW += y[16]/31.014070; /*HNO */
    YOW += y[17]/28.013400; /*N2 */
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
    c[8] = PWORT * y[8]/14.006700; 
    c[9] = PWORT * y[9]/15.014670; 
    c[10] = PWORT * y[10]/16.022640; 
    c[11] = PWORT * y[11]/17.030610; 
    c[12] = PWORT * y[12]/29.021370; 
    c[13] = PWORT * y[13]/30.006100; 
    c[14] = PWORT * y[14]/46.005500; 
    c[15] = PWORT * y[15]/44.012800; 
    c[16] = PWORT * y[16]/31.014070; 
    c[17] = PWORT * y[17]/28.013400; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 68; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[18]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.314e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 18; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 68; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[18]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]/2.015940; 
    c[1] = 1e6 * (*rho) * y[1]/1.007970; 
    c[2] = 1e6 * (*rho) * y[2]/15.999400; 
    c[3] = 1e6 * (*rho) * y[3]/31.998800; 
    c[4] = 1e6 * (*rho) * y[4]/17.007370; 
    c[5] = 1e6 * (*rho) * y[5]/18.015340; 
    c[6] = 1e6 * (*rho) * y[6]/33.006770; 
    c[7] = 1e6 * (*rho) * y[7]/34.014740; 
    c[8] = 1e6 * (*rho) * y[8]/14.006700; 
    c[9] = 1e6 * (*rho) * y[9]/15.014670; 
    c[10] = 1e6 * (*rho) * y[10]/16.022640; 
    c[11] = 1e6 * (*rho) * y[11]/17.030610; 
    c[12] = 1e6 * (*rho) * y[12]/29.021370; 
    c[13] = 1e6 * (*rho) * y[13]/30.006100; 
    c[14] = 1e6 * (*rho) * y[14]/46.005500; 
    c[15] = 1e6 * (*rho) * y[15]/44.012800; 
    c[16] = 1e6 * (*rho) * y[16]/31.014070; 
    c[17] = 1e6 * (*rho) * y[17]/28.013400; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 68; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[18]; /*temporary storage */
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
    XW += x[8]*14.006700; /*N */
    XW += x[9]*15.014670; /*NH */
    XW += x[10]*16.022640; /*NH2 */
    XW += x[11]*17.030610; /*NH3 */
    XW += x[12]*29.021370; /*NNH */
    XW += x[13]*30.006100; /*NO */
    XW += x[14]*46.005500; /*NO2 */
    XW += x[15]*44.012800; /*N2O */
    XW += x[16]*31.014070; /*HNO */
    XW += x[17]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 18; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 68; ++id) {
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
    for (id = 0; id < 18 * 68; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: 2 O + M <=> O2 + M */
    nuki[ 2 * kd + 0 ] = -2 ;
    nuki[ 3 * kd + 0 ] = +1 ;

    /*reaction 2: O + H + M <=> OH + M */
    nuki[ 2 * kd + 1 ] = -1 ;
    nuki[ 1 * kd + 1 ] = -1 ;
    nuki[ 4 * kd + 1 ] = +1 ;

    /*reaction 3: O + H2 <=> H + OH */
    nuki[ 2 * kd + 2 ] = -1 ;
    nuki[ 0 * kd + 2 ] = -1 ;
    nuki[ 1 * kd + 2 ] = +1 ;
    nuki[ 4 * kd + 2 ] = +1 ;

    /*reaction 4: O + HO2 <=> OH + O2 */
    nuki[ 2 * kd + 3 ] = -1 ;
    nuki[ 6 * kd + 3 ] = -1 ;
    nuki[ 4 * kd + 3 ] = +1 ;
    nuki[ 3 * kd + 3 ] = +1 ;

    /*reaction 5: O + H2O2 <=> OH + HO2 */
    nuki[ 2 * kd + 4 ] = -1 ;
    nuki[ 7 * kd + 4 ] = -1 ;
    nuki[ 4 * kd + 4 ] = +1 ;
    nuki[ 6 * kd + 4 ] = +1 ;

    /*reaction 6: H + O2 + M <=> HO2 + M */
    nuki[ 1 * kd + 5 ] = -1 ;
    nuki[ 3 * kd + 5 ] = -1 ;
    nuki[ 6 * kd + 5 ] = +1 ;

    /*reaction 7: H + 2 O2 <=> HO2 + O2 */
    nuki[ 1 * kd + 6 ] = -1 ;
    nuki[ 3 * kd + 6 ] = -2 ;
    nuki[ 6 * kd + 6 ] = +1 ;
    nuki[ 3 * kd + 6 ] = +1 ;

    /*reaction 8: H + O2 + H2O <=> HO2 + H2O */
    nuki[ 1 * kd + 7 ] = -1 ;
    nuki[ 3 * kd + 7 ] = -1 ;
    nuki[ 5 * kd + 7 ] = -1 ;
    nuki[ 6 * kd + 7 ] = +1 ;
    nuki[ 5 * kd + 7 ] = +1 ;

    /*reaction 9: H + O2 + N2 <=> HO2 + N2 */
    nuki[ 1 * kd + 8 ] = -1 ;
    nuki[ 3 * kd + 8 ] = -1 ;
    nuki[ 17 * kd + 8 ] = -1 ;
    nuki[ 6 * kd + 8 ] = +1 ;
    nuki[ 17 * kd + 8 ] = +1 ;

    /*reaction 10: H + O2 <=> O + OH */
    nuki[ 1 * kd + 9 ] = -1 ;
    nuki[ 3 * kd + 9 ] = -1 ;
    nuki[ 2 * kd + 9 ] = +1 ;
    nuki[ 4 * kd + 9 ] = +1 ;

    /*reaction 11: 2 H + M <=> H2 + M */
    nuki[ 1 * kd + 10 ] = -2 ;
    nuki[ 0 * kd + 10 ] = +1 ;

    /*reaction 12: 2 H + H2 <=> 2 H2 */
    nuki[ 1 * kd + 11 ] = -2 ;
    nuki[ 0 * kd + 11 ] = -1 ;
    nuki[ 0 * kd + 11 ] = +2 ;

    /*reaction 13: 2 H + H2O <=> H2 + H2O */
    nuki[ 1 * kd + 12 ] = -2 ;
    nuki[ 5 * kd + 12 ] = -1 ;
    nuki[ 0 * kd + 12 ] = +1 ;
    nuki[ 5 * kd + 12 ] = +1 ;

    /*reaction 14: H + OH + M <=> H2O + M */
    nuki[ 1 * kd + 13 ] = -1 ;
    nuki[ 4 * kd + 13 ] = -1 ;
    nuki[ 5 * kd + 13 ] = +1 ;

    /*reaction 15: H + HO2 <=> O + H2O */
    nuki[ 1 * kd + 14 ] = -1 ;
    nuki[ 6 * kd + 14 ] = -1 ;
    nuki[ 2 * kd + 14 ] = +1 ;
    nuki[ 5 * kd + 14 ] = +1 ;

    /*reaction 16: H + HO2 <=> O2 + H2 */
    nuki[ 1 * kd + 15 ] = -1 ;
    nuki[ 6 * kd + 15 ] = -1 ;
    nuki[ 3 * kd + 15 ] = +1 ;
    nuki[ 0 * kd + 15 ] = +1 ;

    /*reaction 17: H + HO2 <=> 2 OH */
    nuki[ 1 * kd + 16 ] = -1 ;
    nuki[ 6 * kd + 16 ] = -1 ;
    nuki[ 4 * kd + 16 ] = +2 ;

    /*reaction 18: H + H2O2 <=> HO2 + H2 */
    nuki[ 1 * kd + 17 ] = -1 ;
    nuki[ 7 * kd + 17 ] = -1 ;
    nuki[ 6 * kd + 17 ] = +1 ;
    nuki[ 0 * kd + 17 ] = +1 ;

    /*reaction 19: H + H2O2 <=> OH + H2O */
    nuki[ 1 * kd + 18 ] = -1 ;
    nuki[ 7 * kd + 18 ] = -1 ;
    nuki[ 4 * kd + 18 ] = +1 ;
    nuki[ 5 * kd + 18 ] = +1 ;

    /*reaction 20: OH + H2 <=> H + H2O */
    nuki[ 4 * kd + 19 ] = -1 ;
    nuki[ 0 * kd + 19 ] = -1 ;
    nuki[ 1 * kd + 19 ] = +1 ;
    nuki[ 5 * kd + 19 ] = +1 ;

    /*reaction 21: 2 OH (+M) <=> H2O2 (+M) */
    nuki[ 4 * kd + 20 ] = -2 ;
    nuki[ 7 * kd + 20 ] = +1 ;

    /*reaction 22: 2 OH <=> O + H2O */
    nuki[ 4 * kd + 21 ] = -2 ;
    nuki[ 2 * kd + 21 ] = +1 ;
    nuki[ 5 * kd + 21 ] = +1 ;

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    nuki[ 4 * kd + 22 ] = -1 ;
    nuki[ 6 * kd + 22 ] = -1 ;
    nuki[ 3 * kd + 22 ] = +1 ;
    nuki[ 5 * kd + 22 ] = +1 ;

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    nuki[ 4 * kd + 23 ] = -1 ;
    nuki[ 7 * kd + 23 ] = -1 ;
    nuki[ 6 * kd + 23 ] = +1 ;
    nuki[ 5 * kd + 23 ] = +1 ;

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    nuki[ 4 * kd + 24 ] = -1 ;
    nuki[ 7 * kd + 24 ] = -1 ;
    nuki[ 6 * kd + 24 ] = +1 ;
    nuki[ 5 * kd + 24 ] = +1 ;

    /*reaction 26: 2 HO2 <=> O2 + H2O2 */
    nuki[ 6 * kd + 25 ] = -2 ;
    nuki[ 3 * kd + 25 ] = +1 ;
    nuki[ 7 * kd + 25 ] = +1 ;

    /*reaction 27: 2 HO2 <=> O2 + H2O2 */
    nuki[ 6 * kd + 26 ] = -2 ;
    nuki[ 3 * kd + 26 ] = +1 ;
    nuki[ 7 * kd + 26 ] = +1 ;

    /*reaction 28: N + NO <=> N2 + O */
    nuki[ 8 * kd + 27 ] = -1 ;
    nuki[ 13 * kd + 27 ] = -1 ;
    nuki[ 17 * kd + 27 ] = +1 ;
    nuki[ 2 * kd + 27 ] = +1 ;

    /*reaction 29: N + O2 <=> NO + O */
    nuki[ 8 * kd + 28 ] = -1 ;
    nuki[ 3 * kd + 28 ] = -1 ;
    nuki[ 13 * kd + 28 ] = +1 ;
    nuki[ 2 * kd + 28 ] = +1 ;

    /*reaction 30: N + OH <=> NO + H */
    nuki[ 8 * kd + 29 ] = -1 ;
    nuki[ 4 * kd + 29 ] = -1 ;
    nuki[ 13 * kd + 29 ] = +1 ;
    nuki[ 1 * kd + 29 ] = +1 ;

    /*reaction 31: N2O + O <=> N2 + O2 */
    nuki[ 15 * kd + 30 ] = -1 ;
    nuki[ 2 * kd + 30 ] = -1 ;
    nuki[ 17 * kd + 30 ] = +1 ;
    nuki[ 3 * kd + 30 ] = +1 ;

    /*reaction 32: N2O + O <=> 2 NO */
    nuki[ 15 * kd + 31 ] = -1 ;
    nuki[ 2 * kd + 31 ] = -1 ;
    nuki[ 13 * kd + 31 ] = +2 ;

    /*reaction 33: N2O + H <=> N2 + OH */
    nuki[ 15 * kd + 32 ] = -1 ;
    nuki[ 1 * kd + 32 ] = -1 ;
    nuki[ 17 * kd + 32 ] = +1 ;
    nuki[ 4 * kd + 32 ] = +1 ;

    /*reaction 34: N2O + OH <=> N2 + HO2 */
    nuki[ 15 * kd + 33 ] = -1 ;
    nuki[ 4 * kd + 33 ] = -1 ;
    nuki[ 17 * kd + 33 ] = +1 ;
    nuki[ 6 * kd + 33 ] = +1 ;

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    nuki[ 15 * kd + 34 ] = -1 ;
    nuki[ 17 * kd + 34 ] = +1 ;
    nuki[ 2 * kd + 34 ] = +1 ;

    /*reaction 36: HO2 + NO <=> NO2 + OH */
    nuki[ 6 * kd + 35 ] = -1 ;
    nuki[ 13 * kd + 35 ] = -1 ;
    nuki[ 14 * kd + 35 ] = +1 ;
    nuki[ 4 * kd + 35 ] = +1 ;

    /*reaction 37: NO + O + M <=> NO2 + M */
    nuki[ 13 * kd + 36 ] = -1 ;
    nuki[ 2 * kd + 36 ] = -1 ;
    nuki[ 14 * kd + 36 ] = +1 ;

    /*reaction 38: NO2 + O <=> NO + O2 */
    nuki[ 14 * kd + 37 ] = -1 ;
    nuki[ 2 * kd + 37 ] = -1 ;
    nuki[ 13 * kd + 37 ] = +1 ;
    nuki[ 3 * kd + 37 ] = +1 ;

    /*reaction 39: NO2 + H <=> NO + OH */
    nuki[ 14 * kd + 38 ] = -1 ;
    nuki[ 1 * kd + 38 ] = -1 ;
    nuki[ 13 * kd + 38 ] = +1 ;
    nuki[ 4 * kd + 38 ] = +1 ;

    /*reaction 40: NH + O <=> NO + H */
    nuki[ 9 * kd + 39 ] = -1 ;
    nuki[ 2 * kd + 39 ] = -1 ;
    nuki[ 13 * kd + 39 ] = +1 ;
    nuki[ 1 * kd + 39 ] = +1 ;

    /*reaction 41: NH + H <=> N + H2 */
    nuki[ 9 * kd + 40 ] = -1 ;
    nuki[ 1 * kd + 40 ] = -1 ;
    nuki[ 8 * kd + 40 ] = +1 ;
    nuki[ 0 * kd + 40 ] = +1 ;

    /*reaction 42: NH + OH <=> HNO + H */
    nuki[ 9 * kd + 41 ] = -1 ;
    nuki[ 4 * kd + 41 ] = -1 ;
    nuki[ 16 * kd + 41 ] = +1 ;
    nuki[ 1 * kd + 41 ] = +1 ;

    /*reaction 43: NH + OH <=> N + H2O */
    nuki[ 9 * kd + 42 ] = -1 ;
    nuki[ 4 * kd + 42 ] = -1 ;
    nuki[ 8 * kd + 42 ] = +1 ;
    nuki[ 5 * kd + 42 ] = +1 ;

    /*reaction 44: NH + O2 <=> HNO + O */
    nuki[ 9 * kd + 43 ] = -1 ;
    nuki[ 3 * kd + 43 ] = -1 ;
    nuki[ 16 * kd + 43 ] = +1 ;
    nuki[ 2 * kd + 43 ] = +1 ;

    /*reaction 45: NH + O2 <=> NO + OH */
    nuki[ 9 * kd + 44 ] = -1 ;
    nuki[ 3 * kd + 44 ] = -1 ;
    nuki[ 13 * kd + 44 ] = +1 ;
    nuki[ 4 * kd + 44 ] = +1 ;

    /*reaction 46: NH + N <=> N2 + H */
    nuki[ 9 * kd + 45 ] = -1 ;
    nuki[ 8 * kd + 45 ] = -1 ;
    nuki[ 17 * kd + 45 ] = +1 ;
    nuki[ 1 * kd + 45 ] = +1 ;

    /*reaction 47: NH + H2O <=> HNO + H2 */
    nuki[ 9 * kd + 46 ] = -1 ;
    nuki[ 5 * kd + 46 ] = -1 ;
    nuki[ 16 * kd + 46 ] = +1 ;
    nuki[ 0 * kd + 46 ] = +1 ;

    /*reaction 48: NH + NO <=> N2 + OH */
    nuki[ 9 * kd + 47 ] = -1 ;
    nuki[ 13 * kd + 47 ] = -1 ;
    nuki[ 17 * kd + 47 ] = +1 ;
    nuki[ 4 * kd + 47 ] = +1 ;

    /*reaction 49: NH + NO <=> N2O + H */
    nuki[ 9 * kd + 48 ] = -1 ;
    nuki[ 13 * kd + 48 ] = -1 ;
    nuki[ 15 * kd + 48 ] = +1 ;
    nuki[ 1 * kd + 48 ] = +1 ;

    /*reaction 50: NH2 + O <=> OH + NH */
    nuki[ 10 * kd + 49 ] = -1 ;
    nuki[ 2 * kd + 49 ] = -1 ;
    nuki[ 4 * kd + 49 ] = +1 ;
    nuki[ 9 * kd + 49 ] = +1 ;

    /*reaction 51: NH2 + O <=> H + HNO */
    nuki[ 10 * kd + 50 ] = -1 ;
    nuki[ 2 * kd + 50 ] = -1 ;
    nuki[ 1 * kd + 50 ] = +1 ;
    nuki[ 16 * kd + 50 ] = +1 ;

    /*reaction 52: NH2 + H <=> NH + H2 */
    nuki[ 10 * kd + 51 ] = -1 ;
    nuki[ 1 * kd + 51 ] = -1 ;
    nuki[ 9 * kd + 51 ] = +1 ;
    nuki[ 0 * kd + 51 ] = +1 ;

    /*reaction 53: NH2 + OH <=> NH + H2O */
    nuki[ 10 * kd + 52 ] = -1 ;
    nuki[ 4 * kd + 52 ] = -1 ;
    nuki[ 9 * kd + 52 ] = +1 ;
    nuki[ 5 * kd + 52 ] = +1 ;

    /*reaction 54: NNH <=> N2 + H */
    nuki[ 12 * kd + 53 ] = -1 ;
    nuki[ 17 * kd + 53 ] = +1 ;
    nuki[ 1 * kd + 53 ] = +1 ;

    /*reaction 55: NNH + M <=> N2 + H + M */
    nuki[ 12 * kd + 54 ] = -1 ;
    nuki[ 17 * kd + 54 ] = +1 ;
    nuki[ 1 * kd + 54 ] = +1 ;

    /*reaction 56: NNH + O2 <=> HO2 + N2 */
    nuki[ 12 * kd + 55 ] = -1 ;
    nuki[ 3 * kd + 55 ] = -1 ;
    nuki[ 6 * kd + 55 ] = +1 ;
    nuki[ 17 * kd + 55 ] = +1 ;

    /*reaction 57: NNH + O <=> OH + N2 */
    nuki[ 12 * kd + 56 ] = -1 ;
    nuki[ 2 * kd + 56 ] = -1 ;
    nuki[ 4 * kd + 56 ] = +1 ;
    nuki[ 17 * kd + 56 ] = +1 ;

    /*reaction 58: NNH + O <=> NH + NO */
    nuki[ 12 * kd + 57 ] = -1 ;
    nuki[ 2 * kd + 57 ] = -1 ;
    nuki[ 9 * kd + 57 ] = +1 ;
    nuki[ 13 * kd + 57 ] = +1 ;

    /*reaction 59: NNH + H <=> H2 + N2 */
    nuki[ 12 * kd + 58 ] = -1 ;
    nuki[ 1 * kd + 58 ] = -1 ;
    nuki[ 0 * kd + 58 ] = +1 ;
    nuki[ 17 * kd + 58 ] = +1 ;

    /*reaction 60: NNH + OH <=> H2O + N2 */
    nuki[ 12 * kd + 59 ] = -1 ;
    nuki[ 4 * kd + 59 ] = -1 ;
    nuki[ 5 * kd + 59 ] = +1 ;
    nuki[ 17 * kd + 59 ] = +1 ;

    /*reaction 61: H + NO + M <=> HNO + M */
    nuki[ 1 * kd + 60 ] = -1 ;
    nuki[ 13 * kd + 60 ] = -1 ;
    nuki[ 16 * kd + 60 ] = +1 ;

    /*reaction 62: HNO + O <=> NO + OH */
    nuki[ 16 * kd + 61 ] = -1 ;
    nuki[ 2 * kd + 61 ] = -1 ;
    nuki[ 13 * kd + 61 ] = +1 ;
    nuki[ 4 * kd + 61 ] = +1 ;

    /*reaction 63: HNO + H <=> H2 + NO */
    nuki[ 16 * kd + 62 ] = -1 ;
    nuki[ 1 * kd + 62 ] = -1 ;
    nuki[ 0 * kd + 62 ] = +1 ;
    nuki[ 13 * kd + 62 ] = +1 ;

    /*reaction 64: HNO + OH <=> NO + H2O */
    nuki[ 16 * kd + 63 ] = -1 ;
    nuki[ 4 * kd + 63 ] = -1 ;
    nuki[ 13 * kd + 63 ] = +1 ;
    nuki[ 5 * kd + 63 ] = +1 ;

    /*reaction 65: HNO + O2 <=> HO2 + NO */
    nuki[ 16 * kd + 64 ] = -1 ;
    nuki[ 3 * kd + 64 ] = -1 ;
    nuki[ 6 * kd + 64 ] = +1 ;
    nuki[ 13 * kd + 64 ] = +1 ;

    /*reaction 66: NH3 + H <=> NH2 + H2 */
    nuki[ 11 * kd + 65 ] = -1 ;
    nuki[ 1 * kd + 65 ] = -1 ;
    nuki[ 10 * kd + 65 ] = +1 ;
    nuki[ 0 * kd + 65 ] = +1 ;

    /*reaction 67: NH3 + OH <=> NH2 + H2O */
    nuki[ 11 * kd + 66 ] = -1 ;
    nuki[ 4 * kd + 66 ] = -1 ;
    nuki[ 10 * kd + 66 ] = +1 ;
    nuki[ 5 * kd + 66 ] = +1 ;

    /*reaction 68: NH3 + O <=> NH2 + OH */
    nuki[ 11 * kd + 67 ] = -1 ;
    nuki[ 2 * kd + 67 ] = -1 ;
    nuki[ 10 * kd + 67 ] = +1 ;
    nuki[ 4 * kd + 67 ] = +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < 4 * 18; ++ id) {
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

    /*N */
    ncf[ 8 * kd + 3 ] = 1; /*N */

    /*NH */
    ncf[ 9 * kd + 3 ] = 1; /*N */
    ncf[ 9 * kd + 1 ] = 1; /*H */

    /*NH2 */
    ncf[ 10 * kd + 3 ] = 1; /*N */
    ncf[ 10 * kd + 1 ] = 2; /*H */

    /*NH3 */
    ncf[ 11 * kd + 3 ] = 1; /*N */
    ncf[ 11 * kd + 1 ] = 3; /*H */

    /*NNH */
    ncf[ 12 * kd + 3 ] = 2; /*N */
    ncf[ 12 * kd + 1 ] = 1; /*H */

    /*NO */
    ncf[ 13 * kd + 3 ] = 1; /*N */
    ncf[ 13 * kd + 0 ] = 1; /*O */

    /*NO2 */
    ncf[ 14 * kd + 3 ] = 1; /*N */
    ncf[ 14 * kd + 0 ] = 2; /*O */

    /*N2O */
    ncf[ 15 * kd + 3 ] = 2; /*N */
    ncf[ 15 * kd + 0 ] = 1; /*O */

    /*HNO */
    ncf[ 16 * kd + 1 ] = 1; /*H */
    ncf[ 16 * kd + 3 ] = 1; /*N */
    ncf[ 16 * kd + 0 ] = 1; /*O */

    /*N2 */
    ncf[ 17 * kd + 3 ] = 2; /*N */

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

    /*reaction 6: H + O2 + M <=> HO2 + M */
    a[5] = 2.8e+18;
    b[5] = -0.86;
    e[5] = 0;

    /*reaction 7: H + 2 O2 <=> HO2 + O2 */
    a[6] = 3e+20;
    b[6] = -1.72;
    e[6] = 0;

    /*reaction 8: H + O2 + H2O <=> HO2 + H2O */
    a[7] = 9.38e+18;
    b[7] = -0.76;
    e[7] = 0;

    /*reaction 9: H + O2 + N2 <=> HO2 + N2 */
    a[8] = 3.75e+20;
    b[8] = -1.72;
    e[8] = 0;

    /*reaction 10: H + O2 <=> O + OH */
    a[9] = 8.3e+13;
    b[9] = 0;
    e[9] = 14413;

    /*reaction 11: 2 H + M <=> H2 + M */
    a[10] = 1e+18;
    b[10] = -1;
    e[10] = 0;

    /*reaction 12: 2 H + H2 <=> 2 H2 */
    a[11] = 9e+16;
    b[11] = -0.6;
    e[11] = 0;

    /*reaction 13: 2 H + H2O <=> H2 + H2O */
    a[12] = 6e+19;
    b[12] = -1.25;
    e[12] = 0;

    /*reaction 14: H + OH + M <=> H2O + M */
    a[13] = 2.2e+22;
    b[13] = -2;
    e[13] = 0;

    /*reaction 15: H + HO2 <=> O + H2O */
    a[14] = 3.97e+12;
    b[14] = 0;
    e[14] = 671;

    /*reaction 16: H + HO2 <=> O2 + H2 */
    a[15] = 2.8e+13;
    b[15] = 0;
    e[15] = 1068;

    /*reaction 17: H + HO2 <=> 2 OH */
    a[16] = 1.34e+14;
    b[16] = 0;
    e[16] = 635;

    /*reaction 18: H + H2O2 <=> HO2 + H2 */
    a[17] = 1.21e+07;
    b[17] = 2;
    e[17] = 5200;

    /*reaction 19: H + H2O2 <=> OH + H2O */
    a[18] = 1e+13;
    b[18] = 0;
    e[18] = 3600;

    /*reaction 20: OH + H2 <=> H + H2O */
    a[19] = 2.16e+08;
    b[19] = 1.51;
    e[19] = 3430;

    /*reaction 21: 2 OH (+M) <=> H2O2 (+M) */
    a[20] = 7.4e+13;
    b[20] = -0.37;
    e[20] = 0;

    /*reaction 22: 2 OH <=> O + H2O */
    a[21] = 35700;
    b[21] = 2.4;
    e[21] = -2110;

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    a[22] = 2.9e+13;
    b[22] = 0;
    e[22] = -500;

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    a[23] = 1.75e+12;
    b[23] = 0;
    e[23] = 320;

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    a[24] = 5.8e+14;
    b[24] = 0;
    e[24] = 9560;

    /*reaction 26: 2 HO2 <=> O2 + H2O2 */
    a[25] = 1.3e+11;
    b[25] = 0;
    e[25] = -1630;

    /*reaction 27: 2 HO2 <=> O2 + H2O2 */
    a[26] = 4.2e+14;
    b[26] = 0;
    e[26] = 12000;

    /*reaction 28: N + NO <=> N2 + O */
    a[27] = 3.5e+13;
    b[27] = 0;
    e[27] = 330;

    /*reaction 29: N + O2 <=> NO + O */
    a[28] = 2.65e+12;
    b[28] = 0;
    e[28] = 6400;

    /*reaction 30: N + OH <=> NO + H */
    a[29] = 7.333e+13;
    b[29] = 0;
    e[29] = 1120;

    /*reaction 31: N2O + O <=> N2 + O2 */
    a[30] = 1.4e+12;
    b[30] = 0;
    e[30] = 10810;

    /*reaction 32: N2O + O <=> 2 NO */
    a[31] = 2.9e+13;
    b[31] = 0;
    e[31] = 23150;

    /*reaction 33: N2O + H <=> N2 + OH */
    a[32] = 4.4e+14;
    b[32] = 0;
    e[32] = 18880;

    /*reaction 34: N2O + OH <=> N2 + HO2 */
    a[33] = 2e+12;
    b[33] = 0;
    e[33] = 21060;

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    a[34] = 1.3e+11;
    b[34] = 0;
    e[34] = 59620;

    /*reaction 36: HO2 + NO <=> NO2 + OH */
    a[35] = 2.11e+12;
    b[35] = 0;
    e[35] = -480;

    /*reaction 37: NO + O + M <=> NO2 + M */
    a[36] = 1.06e+20;
    b[36] = -1.41;
    e[36] = 0;

    /*reaction 38: NO2 + O <=> NO + O2 */
    a[37] = 3.9e+12;
    b[37] = 0;
    e[37] = -240;

    /*reaction 39: NO2 + H <=> NO + OH */
    a[38] = 1.32e+14;
    b[38] = 0;
    e[38] = 360;

    /*reaction 40: NH + O <=> NO + H */
    a[39] = 5e+13;
    b[39] = 0;
    e[39] = 0;

    /*reaction 41: NH + H <=> N + H2 */
    a[40] = 3.2e+13;
    b[40] = 0;
    e[40] = 330;

    /*reaction 42: NH + OH <=> HNO + H */
    a[41] = 2e+13;
    b[41] = 0;
    e[41] = 0;

    /*reaction 43: NH + OH <=> N + H2O */
    a[42] = 2e+09;
    b[42] = 1.2;
    e[42] = 0;

    /*reaction 44: NH + O2 <=> HNO + O */
    a[43] = 461000;
    b[43] = 2;
    e[43] = 6500;

    /*reaction 45: NH + O2 <=> NO + OH */
    a[44] = 1.28e+06;
    b[44] = 1.5;
    e[44] = 100;

    /*reaction 46: NH + N <=> N2 + H */
    a[45] = 1.5e+13;
    b[45] = 0;
    e[45] = 0;

    /*reaction 47: NH + H2O <=> HNO + H2 */
    a[46] = 2e+13;
    b[46] = 0;
    e[46] = 13850;

    /*reaction 48: NH + NO <=> N2 + OH */
    a[47] = 2.16e+13;
    b[47] = -0.23;
    e[47] = 0;

    /*reaction 49: NH + NO <=> N2O + H */
    a[48] = 4.16e+14;
    b[48] = -0.45;
    e[48] = 0;

    /*reaction 50: NH2 + O <=> OH + NH */
    a[49] = 7e+12;
    b[49] = 0;
    e[49] = 0;

    /*reaction 51: NH2 + O <=> H + HNO */
    a[50] = 4.6e+13;
    b[50] = 0;
    e[50] = 0;

    /*reaction 52: NH2 + H <=> NH + H2 */
    a[51] = 4e+13;
    b[51] = 0;
    e[51] = 3650;

    /*reaction 53: NH2 + OH <=> NH + H2O */
    a[52] = 9e+07;
    b[52] = 1.5;
    e[52] = -460;

    /*reaction 54: NNH <=> N2 + H */
    a[53] = 3.3e+08;
    b[53] = 0;
    e[53] = 0;

    /*reaction 55: NNH + M <=> N2 + H + M */
    a[54] = 1.3e+14;
    b[54] = -0.11;
    e[54] = 4980;

    /*reaction 56: NNH + O2 <=> HO2 + N2 */
    a[55] = 5e+12;
    b[55] = 0;
    e[55] = 0;

    /*reaction 57: NNH + O <=> OH + N2 */
    a[56] = 2.5e+13;
    b[56] = 0;
    e[56] = 0;

    /*reaction 58: NNH + O <=> NH + NO */
    a[57] = 7e+13;
    b[57] = 0;
    e[57] = 0;

    /*reaction 59: NNH + H <=> H2 + N2 */
    a[58] = 5e+13;
    b[58] = 0;
    e[58] = 0;

    /*reaction 60: NNH + OH <=> H2O + N2 */
    a[59] = 2e+13;
    b[59] = 0;
    e[59] = 0;

    /*reaction 61: H + NO + M <=> HNO + M */
    a[60] = 8.95e+19;
    b[60] = -1.32;
    e[60] = 740;

    /*reaction 62: HNO + O <=> NO + OH */
    a[61] = 2.5e+13;
    b[61] = 0;
    e[61] = 0;

    /*reaction 63: HNO + H <=> H2 + NO */
    a[62] = 4.5e+11;
    b[62] = 0.72;
    e[62] = 660;

    /*reaction 64: HNO + OH <=> NO + H2O */
    a[63] = 1.3e+07;
    b[63] = 1.9;
    e[63] = -950;

    /*reaction 65: HNO + O2 <=> HO2 + NO */
    a[64] = 1e+13;
    b[64] = 0;
    e[64] = 13000;

    /*reaction 66: NH3 + H <=> NH2 + H2 */
    a[65] = 540000;
    b[65] = 2.4;
    e[65] = 9915;

    /*reaction 67: NH3 + OH <=> NH2 + H2O */
    a[66] = 5e+07;
    b[66] = 1.6;
    e[66] = 955;

    /*reaction 68: NH3 + O <=> NH2 + OH */
    a[67] = 9.4e+06;
    b[67] = 1.94;
    e[67] = 6460;

    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double * T, double * C, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[18]; /* temporary storage */

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

    /*reaction 6: H + O2 + M <=> HO2 + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + 2 O2 <=> HO2 + O2 */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 + H2O <=> HO2 + H2O */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H + O2 + N2 <=> HO2 + N2 */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 <=> O + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: 2 H + M <=> H2 + M */
    eqcon[10] *= 1e+06; 

    /*reaction 12: 2 H + H2 <=> 2 H2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: 2 H + H2O <=> H2 + H2O */
    eqcon[12] *= 1e+06; 

    /*reaction 14: H + OH + M <=> H2O + M */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + HO2 <=> O + H2O */
    /*eqcon[14] *= 1;  */

    /*reaction 16: H + HO2 <=> O2 + H2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + HO2 <=> 2 OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H + H2O2 <=> HO2 + H2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + H2O2 <=> OH + H2O */
    /*eqcon[18] *= 1;  */

    /*reaction 20: OH + H2 <=> H + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[20] *= 1e+06; 

    /*reaction 22: 2 OH <=> O + H2O */
    /*eqcon[21] *= 1;  */

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[23] *= 1;  */

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: N + NO <=> N2 + O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: N + O2 <=> NO + O */
    /*eqcon[28] *= 1;  */

    /*reaction 30: N + OH <=> NO + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: N2O + O <=> N2 + O2 */
    /*eqcon[30] *= 1;  */

    /*reaction 32: N2O + O <=> 2 NO */
    /*eqcon[31] *= 1;  */

    /*reaction 33: N2O + H <=> N2 + OH */
    /*eqcon[32] *= 1;  */

    /*reaction 34: N2O + OH <=> N2 + HO2 */
    /*eqcon[33] *= 1;  */

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    eqcon[34] *= 1e-06; 

    /*reaction 36: HO2 + NO <=> NO2 + OH */
    /*eqcon[35] *= 1;  */

    /*reaction 37: NO + O + M <=> NO2 + M */
    eqcon[36] *= 1e+06; 

    /*reaction 38: NO2 + O <=> NO + O2 */
    /*eqcon[37] *= 1;  */

    /*reaction 39: NO2 + H <=> NO + OH */
    /*eqcon[38] *= 1;  */

    /*reaction 40: NH + O <=> NO + H */
    /*eqcon[39] *= 1;  */

    /*reaction 41: NH + H <=> N + H2 */
    /*eqcon[40] *= 1;  */

    /*reaction 42: NH + OH <=> HNO + H */
    /*eqcon[41] *= 1;  */

    /*reaction 43: NH + OH <=> N + H2O */
    /*eqcon[42] *= 1;  */

    /*reaction 44: NH + O2 <=> HNO + O */
    /*eqcon[43] *= 1;  */

    /*reaction 45: NH + O2 <=> NO + OH */
    /*eqcon[44] *= 1;  */

    /*reaction 46: NH + N <=> N2 + H */
    /*eqcon[45] *= 1;  */

    /*reaction 47: NH + H2O <=> HNO + H2 */
    /*eqcon[46] *= 1;  */

    /*reaction 48: NH + NO <=> N2 + OH */
    /*eqcon[47] *= 1;  */

    /*reaction 49: NH + NO <=> N2O + H */
    /*eqcon[48] *= 1;  */

    /*reaction 50: NH2 + O <=> OH + NH */
    /*eqcon[49] *= 1;  */

    /*reaction 51: NH2 + O <=> H + HNO */
    /*eqcon[50] *= 1;  */

    /*reaction 52: NH2 + H <=> NH + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: NH2 + OH <=> NH + H2O */
    /*eqcon[52] *= 1;  */

    /*reaction 54: NNH <=> N2 + H */
    eqcon[53] *= 1e-06; 

    /*reaction 55: NNH + M <=> N2 + H + M */
    eqcon[54] *= 1e-06; 

    /*reaction 56: NNH + O2 <=> HO2 + N2 */
    /*eqcon[55] *= 1;  */

    /*reaction 57: NNH + O <=> OH + N2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: NNH + O <=> NH + NO */
    /*eqcon[57] *= 1;  */

    /*reaction 59: NNH + H <=> H2 + N2 */
    /*eqcon[58] *= 1;  */

    /*reaction 60: NNH + OH <=> H2O + N2 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: H + NO + M <=> HNO + M */
    eqcon[60] *= 1e+06; 

    /*reaction 62: HNO + O <=> NO + OH */
    /*eqcon[61] *= 1;  */

    /*reaction 63: HNO + H <=> H2 + NO */
    /*eqcon[62] *= 1;  */

    /*reaction 64: HNO + OH <=> NO + H2O */
    /*eqcon[63] *= 1;  */

    /*reaction 65: HNO + O2 <=> HO2 + NO */
    /*eqcon[64] *= 1;  */

    /*reaction 66: NH3 + H <=> NH2 + H2 */
    /*eqcon[65] *= 1;  */

    /*reaction 67: NH3 + OH <=> NH2 + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: NH3 + O <=> NH2 + OH */
    /*eqcon[67] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[18]; /* temporary storage */

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

    /*reaction 6: H + O2 + M <=> HO2 + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + 2 O2 <=> HO2 + O2 */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 + H2O <=> HO2 + H2O */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H + O2 + N2 <=> HO2 + N2 */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 <=> O + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: 2 H + M <=> H2 + M */
    eqcon[10] *= 1e+06; 

    /*reaction 12: 2 H + H2 <=> 2 H2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: 2 H + H2O <=> H2 + H2O */
    eqcon[12] *= 1e+06; 

    /*reaction 14: H + OH + M <=> H2O + M */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + HO2 <=> O + H2O */
    /*eqcon[14] *= 1;  */

    /*reaction 16: H + HO2 <=> O2 + H2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + HO2 <=> 2 OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H + H2O2 <=> HO2 + H2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + H2O2 <=> OH + H2O */
    /*eqcon[18] *= 1;  */

    /*reaction 20: OH + H2 <=> H + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[20] *= 1e+06; 

    /*reaction 22: 2 OH <=> O + H2O */
    /*eqcon[21] *= 1;  */

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[23] *= 1;  */

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: N + NO <=> N2 + O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: N + O2 <=> NO + O */
    /*eqcon[28] *= 1;  */

    /*reaction 30: N + OH <=> NO + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: N2O + O <=> N2 + O2 */
    /*eqcon[30] *= 1;  */

    /*reaction 32: N2O + O <=> 2 NO */
    /*eqcon[31] *= 1;  */

    /*reaction 33: N2O + H <=> N2 + OH */
    /*eqcon[32] *= 1;  */

    /*reaction 34: N2O + OH <=> N2 + HO2 */
    /*eqcon[33] *= 1;  */

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    eqcon[34] *= 1e-06; 

    /*reaction 36: HO2 + NO <=> NO2 + OH */
    /*eqcon[35] *= 1;  */

    /*reaction 37: NO + O + M <=> NO2 + M */
    eqcon[36] *= 1e+06; 

    /*reaction 38: NO2 + O <=> NO + O2 */
    /*eqcon[37] *= 1;  */

    /*reaction 39: NO2 + H <=> NO + OH */
    /*eqcon[38] *= 1;  */

    /*reaction 40: NH + O <=> NO + H */
    /*eqcon[39] *= 1;  */

    /*reaction 41: NH + H <=> N + H2 */
    /*eqcon[40] *= 1;  */

    /*reaction 42: NH + OH <=> HNO + H */
    /*eqcon[41] *= 1;  */

    /*reaction 43: NH + OH <=> N + H2O */
    /*eqcon[42] *= 1;  */

    /*reaction 44: NH + O2 <=> HNO + O */
    /*eqcon[43] *= 1;  */

    /*reaction 45: NH + O2 <=> NO + OH */
    /*eqcon[44] *= 1;  */

    /*reaction 46: NH + N <=> N2 + H */
    /*eqcon[45] *= 1;  */

    /*reaction 47: NH + H2O <=> HNO + H2 */
    /*eqcon[46] *= 1;  */

    /*reaction 48: NH + NO <=> N2 + OH */
    /*eqcon[47] *= 1;  */

    /*reaction 49: NH + NO <=> N2O + H */
    /*eqcon[48] *= 1;  */

    /*reaction 50: NH2 + O <=> OH + NH */
    /*eqcon[49] *= 1;  */

    /*reaction 51: NH2 + O <=> H + HNO */
    /*eqcon[50] *= 1;  */

    /*reaction 52: NH2 + H <=> NH + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: NH2 + OH <=> NH + H2O */
    /*eqcon[52] *= 1;  */

    /*reaction 54: NNH <=> N2 + H */
    eqcon[53] *= 1e-06; 

    /*reaction 55: NNH + M <=> N2 + H + M */
    eqcon[54] *= 1e-06; 

    /*reaction 56: NNH + O2 <=> HO2 + N2 */
    /*eqcon[55] *= 1;  */

    /*reaction 57: NNH + O <=> OH + N2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: NNH + O <=> NH + NO */
    /*eqcon[57] *= 1;  */

    /*reaction 59: NNH + H <=> H2 + N2 */
    /*eqcon[58] *= 1;  */

    /*reaction 60: NNH + OH <=> H2O + N2 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: H + NO + M <=> HNO + M */
    eqcon[60] *= 1e+06; 

    /*reaction 62: HNO + O <=> NO + OH */
    /*eqcon[61] *= 1;  */

    /*reaction 63: HNO + H <=> H2 + NO */
    /*eqcon[62] *= 1;  */

    /*reaction 64: HNO + OH <=> NO + H2O */
    /*eqcon[63] *= 1;  */

    /*reaction 65: HNO + O2 <=> HO2 + NO */
    /*eqcon[64] *= 1;  */

    /*reaction 66: NH3 + H <=> NH2 + H2 */
    /*eqcon[65] *= 1;  */

    /*reaction 67: NH3 + OH <=> NH2 + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: NH3 + O <=> NH2 + OH */
    /*eqcon[67] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[18]; /* temporary storage */

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

    /*reaction 6: H + O2 + M <=> HO2 + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + 2 O2 <=> HO2 + O2 */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 + H2O <=> HO2 + H2O */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H + O2 + N2 <=> HO2 + N2 */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 <=> O + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: 2 H + M <=> H2 + M */
    eqcon[10] *= 1e+06; 

    /*reaction 12: 2 H + H2 <=> 2 H2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: 2 H + H2O <=> H2 + H2O */
    eqcon[12] *= 1e+06; 

    /*reaction 14: H + OH + M <=> H2O + M */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + HO2 <=> O + H2O */
    /*eqcon[14] *= 1;  */

    /*reaction 16: H + HO2 <=> O2 + H2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + HO2 <=> 2 OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H + H2O2 <=> HO2 + H2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + H2O2 <=> OH + H2O */
    /*eqcon[18] *= 1;  */

    /*reaction 20: OH + H2 <=> H + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[20] *= 1e+06; 

    /*reaction 22: 2 OH <=> O + H2O */
    /*eqcon[21] *= 1;  */

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[23] *= 1;  */

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: N + NO <=> N2 + O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: N + O2 <=> NO + O */
    /*eqcon[28] *= 1;  */

    /*reaction 30: N + OH <=> NO + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: N2O + O <=> N2 + O2 */
    /*eqcon[30] *= 1;  */

    /*reaction 32: N2O + O <=> 2 NO */
    /*eqcon[31] *= 1;  */

    /*reaction 33: N2O + H <=> N2 + OH */
    /*eqcon[32] *= 1;  */

    /*reaction 34: N2O + OH <=> N2 + HO2 */
    /*eqcon[33] *= 1;  */

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    eqcon[34] *= 1e-06; 

    /*reaction 36: HO2 + NO <=> NO2 + OH */
    /*eqcon[35] *= 1;  */

    /*reaction 37: NO + O + M <=> NO2 + M */
    eqcon[36] *= 1e+06; 

    /*reaction 38: NO2 + O <=> NO + O2 */
    /*eqcon[37] *= 1;  */

    /*reaction 39: NO2 + H <=> NO + OH */
    /*eqcon[38] *= 1;  */

    /*reaction 40: NH + O <=> NO + H */
    /*eqcon[39] *= 1;  */

    /*reaction 41: NH + H <=> N + H2 */
    /*eqcon[40] *= 1;  */

    /*reaction 42: NH + OH <=> HNO + H */
    /*eqcon[41] *= 1;  */

    /*reaction 43: NH + OH <=> N + H2O */
    /*eqcon[42] *= 1;  */

    /*reaction 44: NH + O2 <=> HNO + O */
    /*eqcon[43] *= 1;  */

    /*reaction 45: NH + O2 <=> NO + OH */
    /*eqcon[44] *= 1;  */

    /*reaction 46: NH + N <=> N2 + H */
    /*eqcon[45] *= 1;  */

    /*reaction 47: NH + H2O <=> HNO + H2 */
    /*eqcon[46] *= 1;  */

    /*reaction 48: NH + NO <=> N2 + OH */
    /*eqcon[47] *= 1;  */

    /*reaction 49: NH + NO <=> N2O + H */
    /*eqcon[48] *= 1;  */

    /*reaction 50: NH2 + O <=> OH + NH */
    /*eqcon[49] *= 1;  */

    /*reaction 51: NH2 + O <=> H + HNO */
    /*eqcon[50] *= 1;  */

    /*reaction 52: NH2 + H <=> NH + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: NH2 + OH <=> NH + H2O */
    /*eqcon[52] *= 1;  */

    /*reaction 54: NNH <=> N2 + H */
    eqcon[53] *= 1e-06; 

    /*reaction 55: NNH + M <=> N2 + H + M */
    eqcon[54] *= 1e-06; 

    /*reaction 56: NNH + O2 <=> HO2 + N2 */
    /*eqcon[55] *= 1;  */

    /*reaction 57: NNH + O <=> OH + N2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: NNH + O <=> NH + NO */
    /*eqcon[57] *= 1;  */

    /*reaction 59: NNH + H <=> H2 + N2 */
    /*eqcon[58] *= 1;  */

    /*reaction 60: NNH + OH <=> H2O + N2 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: H + NO + M <=> HNO + M */
    eqcon[60] *= 1e+06; 

    /*reaction 62: HNO + O <=> NO + OH */
    /*eqcon[61] *= 1;  */

    /*reaction 63: HNO + H <=> H2 + NO */
    /*eqcon[62] *= 1;  */

    /*reaction 64: HNO + OH <=> NO + H2O */
    /*eqcon[63] *= 1;  */

    /*reaction 65: HNO + O2 <=> HO2 + NO */
    /*eqcon[64] *= 1;  */

    /*reaction 66: NH3 + H <=> NH2 + H2 */
    /*eqcon[65] *= 1;  */

    /*reaction 67: NH3 + OH <=> NH2 + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: NH3 + O <=> NH2 + OH */
    /*eqcon[67] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[18]; /* temporary storage */

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

    /*reaction 6: H + O2 + M <=> HO2 + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + 2 O2 <=> HO2 + O2 */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 + H2O <=> HO2 + H2O */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H + O2 + N2 <=> HO2 + N2 */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 <=> O + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: 2 H + M <=> H2 + M */
    eqcon[10] *= 1e+06; 

    /*reaction 12: 2 H + H2 <=> 2 H2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: 2 H + H2O <=> H2 + H2O */
    eqcon[12] *= 1e+06; 

    /*reaction 14: H + OH + M <=> H2O + M */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + HO2 <=> O + H2O */
    /*eqcon[14] *= 1;  */

    /*reaction 16: H + HO2 <=> O2 + H2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + HO2 <=> 2 OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H + H2O2 <=> HO2 + H2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + H2O2 <=> OH + H2O */
    /*eqcon[18] *= 1;  */

    /*reaction 20: OH + H2 <=> H + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[20] *= 1e+06; 

    /*reaction 22: 2 OH <=> O + H2O */
    /*eqcon[21] *= 1;  */

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[23] *= 1;  */

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: N + NO <=> N2 + O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: N + O2 <=> NO + O */
    /*eqcon[28] *= 1;  */

    /*reaction 30: N + OH <=> NO + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: N2O + O <=> N2 + O2 */
    /*eqcon[30] *= 1;  */

    /*reaction 32: N2O + O <=> 2 NO */
    /*eqcon[31] *= 1;  */

    /*reaction 33: N2O + H <=> N2 + OH */
    /*eqcon[32] *= 1;  */

    /*reaction 34: N2O + OH <=> N2 + HO2 */
    /*eqcon[33] *= 1;  */

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    eqcon[34] *= 1e-06; 

    /*reaction 36: HO2 + NO <=> NO2 + OH */
    /*eqcon[35] *= 1;  */

    /*reaction 37: NO + O + M <=> NO2 + M */
    eqcon[36] *= 1e+06; 

    /*reaction 38: NO2 + O <=> NO + O2 */
    /*eqcon[37] *= 1;  */

    /*reaction 39: NO2 + H <=> NO + OH */
    /*eqcon[38] *= 1;  */

    /*reaction 40: NH + O <=> NO + H */
    /*eqcon[39] *= 1;  */

    /*reaction 41: NH + H <=> N + H2 */
    /*eqcon[40] *= 1;  */

    /*reaction 42: NH + OH <=> HNO + H */
    /*eqcon[41] *= 1;  */

    /*reaction 43: NH + OH <=> N + H2O */
    /*eqcon[42] *= 1;  */

    /*reaction 44: NH + O2 <=> HNO + O */
    /*eqcon[43] *= 1;  */

    /*reaction 45: NH + O2 <=> NO + OH */
    /*eqcon[44] *= 1;  */

    /*reaction 46: NH + N <=> N2 + H */
    /*eqcon[45] *= 1;  */

    /*reaction 47: NH + H2O <=> HNO + H2 */
    /*eqcon[46] *= 1;  */

    /*reaction 48: NH + NO <=> N2 + OH */
    /*eqcon[47] *= 1;  */

    /*reaction 49: NH + NO <=> N2O + H */
    /*eqcon[48] *= 1;  */

    /*reaction 50: NH2 + O <=> OH + NH */
    /*eqcon[49] *= 1;  */

    /*reaction 51: NH2 + O <=> H + HNO */
    /*eqcon[50] *= 1;  */

    /*reaction 52: NH2 + H <=> NH + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: NH2 + OH <=> NH + H2O */
    /*eqcon[52] *= 1;  */

    /*reaction 54: NNH <=> N2 + H */
    eqcon[53] *= 1e-06; 

    /*reaction 55: NNH + M <=> N2 + H + M */
    eqcon[54] *= 1e-06; 

    /*reaction 56: NNH + O2 <=> HO2 + N2 */
    /*eqcon[55] *= 1;  */

    /*reaction 57: NNH + O <=> OH + N2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: NNH + O <=> NH + NO */
    /*eqcon[57] *= 1;  */

    /*reaction 59: NNH + H <=> H2 + N2 */
    /*eqcon[58] *= 1;  */

    /*reaction 60: NNH + OH <=> H2O + N2 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: H + NO + M <=> HNO + M */
    eqcon[60] *= 1e+06; 

    /*reaction 62: HNO + O <=> NO + OH */
    /*eqcon[61] *= 1;  */

    /*reaction 63: HNO + H <=> H2 + NO */
    /*eqcon[62] *= 1;  */

    /*reaction 64: HNO + OH <=> NO + H2O */
    /*eqcon[63] *= 1;  */

    /*reaction 65: HNO + O2 <=> HO2 + NO */
    /*eqcon[64] *= 1;  */

    /*reaction 66: NH3 + H <=> NH2 + H2 */
    /*eqcon[65] *= 1;  */

    /*reaction 67: NH3 + OH <=> NH2 + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: NH3 + O <=> NH2 + OH */
    /*eqcon[67] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[18]; /* temporary storage */

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

    /*reaction 6: H + O2 + M <=> HO2 + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + 2 O2 <=> HO2 + O2 */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + O2 + H2O <=> HO2 + H2O */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H + O2 + N2 <=> HO2 + N2 */
    eqcon[8] *= 1e+06; 

    /*reaction 10: H + O2 <=> O + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: 2 H + M <=> H2 + M */
    eqcon[10] *= 1e+06; 

    /*reaction 12: 2 H + H2 <=> 2 H2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: 2 H + H2O <=> H2 + H2O */
    eqcon[12] *= 1e+06; 

    /*reaction 14: H + OH + M <=> H2O + M */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + HO2 <=> O + H2O */
    /*eqcon[14] *= 1;  */

    /*reaction 16: H + HO2 <=> O2 + H2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + HO2 <=> 2 OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H + H2O2 <=> HO2 + H2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + H2O2 <=> OH + H2O */
    /*eqcon[18] *= 1;  */

    /*reaction 20: OH + H2 <=> H + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[20] *= 1e+06; 

    /*reaction 22: 2 OH <=> O + H2O */
    /*eqcon[21] *= 1;  */

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[23] *= 1;  */

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: N + NO <=> N2 + O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: N + O2 <=> NO + O */
    /*eqcon[28] *= 1;  */

    /*reaction 30: N + OH <=> NO + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: N2O + O <=> N2 + O2 */
    /*eqcon[30] *= 1;  */

    /*reaction 32: N2O + O <=> 2 NO */
    /*eqcon[31] *= 1;  */

    /*reaction 33: N2O + H <=> N2 + OH */
    /*eqcon[32] *= 1;  */

    /*reaction 34: N2O + OH <=> N2 + HO2 */
    /*eqcon[33] *= 1;  */

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    eqcon[34] *= 1e-06; 

    /*reaction 36: HO2 + NO <=> NO2 + OH */
    /*eqcon[35] *= 1;  */

    /*reaction 37: NO + O + M <=> NO2 + M */
    eqcon[36] *= 1e+06; 

    /*reaction 38: NO2 + O <=> NO + O2 */
    /*eqcon[37] *= 1;  */

    /*reaction 39: NO2 + H <=> NO + OH */
    /*eqcon[38] *= 1;  */

    /*reaction 40: NH + O <=> NO + H */
    /*eqcon[39] *= 1;  */

    /*reaction 41: NH + H <=> N + H2 */
    /*eqcon[40] *= 1;  */

    /*reaction 42: NH + OH <=> HNO + H */
    /*eqcon[41] *= 1;  */

    /*reaction 43: NH + OH <=> N + H2O */
    /*eqcon[42] *= 1;  */

    /*reaction 44: NH + O2 <=> HNO + O */
    /*eqcon[43] *= 1;  */

    /*reaction 45: NH + O2 <=> NO + OH */
    /*eqcon[44] *= 1;  */

    /*reaction 46: NH + N <=> N2 + H */
    /*eqcon[45] *= 1;  */

    /*reaction 47: NH + H2O <=> HNO + H2 */
    /*eqcon[46] *= 1;  */

    /*reaction 48: NH + NO <=> N2 + OH */
    /*eqcon[47] *= 1;  */

    /*reaction 49: NH + NO <=> N2O + H */
    /*eqcon[48] *= 1;  */

    /*reaction 50: NH2 + O <=> OH + NH */
    /*eqcon[49] *= 1;  */

    /*reaction 51: NH2 + O <=> H + HNO */
    /*eqcon[50] *= 1;  */

    /*reaction 52: NH2 + H <=> NH + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: NH2 + OH <=> NH + H2O */
    /*eqcon[52] *= 1;  */

    /*reaction 54: NNH <=> N2 + H */
    eqcon[53] *= 1e-06; 

    /*reaction 55: NNH + M <=> N2 + H + M */
    eqcon[54] *= 1e-06; 

    /*reaction 56: NNH + O2 <=> HO2 + N2 */
    /*eqcon[55] *= 1;  */

    /*reaction 57: NNH + O <=> OH + N2 */
    /*eqcon[56] *= 1;  */

    /*reaction 58: NNH + O <=> NH + NO */
    /*eqcon[57] *= 1;  */

    /*reaction 59: NNH + H <=> H2 + N2 */
    /*eqcon[58] *= 1;  */

    /*reaction 60: NNH + OH <=> H2O + N2 */
    /*eqcon[59] *= 1;  */

    /*reaction 61: H + NO + M <=> HNO + M */
    eqcon[60] *= 1e+06; 

    /*reaction 62: HNO + O <=> NO + OH */
    /*eqcon[61] *= 1;  */

    /*reaction 63: HNO + H <=> H2 + NO */
    /*eqcon[62] *= 1;  */

    /*reaction 64: HNO + OH <=> NO + H2O */
    /*eqcon[63] *= 1;  */

    /*reaction 65: HNO + O2 <=> HO2 + NO */
    /*eqcon[64] *= 1;  */

    /*reaction 66: NH3 + H <=> NH2 + H2 */
    /*eqcon[65] *= 1;  */

    /*reaction 67: NH3 + OH <=> NH2 + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: NH3 + O <=> NH2 + OH */
    /*eqcon[67] *= 1;  */
}


/*compute the production rate for each species */
void productionRate(double * wdot, double * sc, double T)
{
    double qdot;

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[18];                /*Gibbs free energy */
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
    for (id = 0; id < 18; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*zero out wdot */
    for (id = 0; id < 18; ++id) {
        wdot[id] = 0.0;
    }

    /*reaction 1: 2 O + M <=> O2 + M */
    phi_f = sc[2]*sc[2];
    alpha = mixture + 1.4*sc[0] + 14.4*sc[5];
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
    alpha = mixture + sc[0] + 5*sc[5];
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

    /*reaction 6: H + O2 + M <=> HO2 + M */
    phi_f = sc[1]*sc[3];
    alpha = mixture + -1*sc[3] + -1*sc[5] + -1*sc[17];
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

    /*reaction 7: H + 2 O2 <=> HO2 + O2 */
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

    /*reaction 8: H + O2 + H2O <=> HO2 + H2O */
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

    /*reaction 9: H + O2 + N2 <=> HO2 + N2 */
    phi_f = sc[1]*sc[3]*sc[17];
    k_f = 1e-12 * 3.75e+20*exp(-1.72*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[17]) - (g_RT[6] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[17] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 10: H + O2 <=> O + OH */
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

    /*reaction 11: 2 H + M <=> H2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[0] + -1*sc[5];
    k_f = 1e-12 * alpha * 1e+18*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 12: 2 H + H2 <=> 2 H2 */
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

    /*reaction 13: 2 H + H2O <=> H2 + H2O */
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

    /*reaction 14: H + OH + M <=> H2O + M */
    phi_f = sc[1]*sc[4];
    alpha = mixture + -0.27*sc[0] + 2.65*sc[5];
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

    /*reaction 15: H + HO2 <=> O + H2O */
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

    /*reaction 16: H + HO2 <=> O2 + H2 */
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

    /*reaction 17: H + HO2 <=> 2 OH */
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

    /*reaction 18: H + H2O2 <=> HO2 + H2 */
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

    /*reaction 19: H + H2O2 <=> OH + H2O */
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

    /*reaction 20: OH + H2 <=> H + H2O */
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

    /*reaction 21: 2 OH (+M) <=> H2O2 (+M) */
    phi_f = sc[4]*sc[4];
    alpha = mixture + sc[0] + 5*sc[5];
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

    /*reaction 22: 2 OH <=> O + H2O */
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

    /*reaction 23: OH + HO2 <=> O2 + H2O */
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

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
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

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
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

    /*reaction 26: 2 HO2 <=> O2 + H2O2 */
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

    /*reaction 27: 2 HO2 <=> O2 + H2O2 */
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

    /*reaction 28: N + NO <=> N2 + O */
    phi_f = sc[8]*sc[13];
    k_f = 1e-06 * 3.5e+13*exp(-166.08/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[2];
    Kc = exp((g_RT[8] + g_RT[13]) - (g_RT[17] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 29: N + O2 <=> NO + O */
    phi_f = sc[8]*sc[3];
    k_f = 1e-06 * 2.65e+12*exp(-3220.94/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[2];
    Kc = exp((g_RT[8] + g_RT[3]) - (g_RT[13] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 30: N + OH <=> NO + H */
    phi_f = sc[8]*sc[4];
    k_f = 1e-06 * 7.333e+13*exp(-563.664/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[1];
    Kc = exp((g_RT[8] + g_RT[4]) - (g_RT[13] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[8] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 31: N2O + O <=> N2 + O2 */
    phi_f = sc[15]*sc[2];
    k_f = 1e-06 * 1.4e+12*exp(-5440.36/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[3];
    Kc = exp((g_RT[15] + g_RT[2]) - (g_RT[17] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 32: N2O + O <=> 2 NO */
    phi_f = sc[15]*sc[2];
    k_f = 1e-06 * 2.9e+13*exp(-11650.7/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[13];
    Kc = exp((g_RT[15] + g_RT[2]) - (2 * g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[13] += 2 * qdot;

    /*reaction 33: N2O + H <=> N2 + OH */
    phi_f = sc[15]*sc[1];
    k_f = 1e-06 * 4.4e+14*exp(-9501.76/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[4];
    Kc = exp((g_RT[15] + g_RT[1]) - (g_RT[17] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 34: N2O + OH <=> N2 + HO2 */
    phi_f = sc[15]*sc[4];
    k_f = 1e-06 * 2e+12*exp(-10598.9/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[6];
    Kc = exp((g_RT[15] + g_RT[4]) - (g_RT[17] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    phi_f = sc[15];
    alpha = mixture + sc[0] + 5*sc[5];
    k_f = 1 * 1.3e+11*exp(-30005/tc[1]);
    redP = 1e-12 * alpha / k_f * 6.2e+14*exp(-28233.5/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[2];
    Kc = refC * exp((g_RT[15]) - (g_RT[17] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[15] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 36: HO2 + NO <=> NO2 + OH */
    phi_f = sc[6]*sc[13];
    k_f = 1e-06 * 2.11e+12*exp(+241.57/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[4];
    Kc = exp((g_RT[6] + g_RT[13]) - (g_RT[14] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[14] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 37: NO + O + M <=> NO2 + M */
    phi_f = sc[13]*sc[2];
    alpha = mixture + sc[0] + 5*sc[5];
    k_f = 1e-12 * alpha * 1.06e+20*exp(-1.41*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[14];
    Kc = 1.0 / (refC) * exp((g_RT[13] + g_RT[2]) - (g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[13] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[14] += 1 * qdot;

    /*reaction 38: NO2 + O <=> NO + O2 */
    phi_f = sc[14]*sc[2];
    k_f = 1e-06 * 3.9e+12*exp(+120.785/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[3];
    Kc = exp((g_RT[14] + g_RT[2]) - (g_RT[13] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[14] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 39: NO2 + H <=> NO + OH */
    phi_f = sc[14]*sc[1];
    k_f = 1e-06 * 1.32e+14*exp(-181.178/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[4];
    Kc = exp((g_RT[14] + g_RT[1]) - (g_RT[13] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[14] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 40: NH + O <=> NO + H */
    phi_f = sc[9]*sc[2];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[1];
    Kc = exp((g_RT[9] + g_RT[2]) - (g_RT[13] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 41: NH + H <=> N + H2 */
    phi_f = sc[9]*sc[1];
    k_f = 1e-06 * 3.2e+13*exp(-166.08/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = exp((g_RT[9] + g_RT[1]) - (g_RT[8] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 42: NH + OH <=> HNO + H */
    phi_f = sc[9]*sc[4];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[1];
    Kc = exp((g_RT[9] + g_RT[4]) - (g_RT[16] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 43: NH + OH <=> N + H2O */
    phi_f = sc[9]*sc[4];
    k_f = 1e-06 * 2e+09*exp(1.2*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[5];
    Kc = exp((g_RT[9] + g_RT[4]) - (g_RT[8] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[8] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 44: NH + O2 <=> HNO + O */
    phi_f = sc[9]*sc[3];
    k_f = 1e-06 * 461000*exp(2*tc[0]-3271.26/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[2];
    Kc = exp((g_RT[9] + g_RT[3]) - (g_RT[16] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 45: NH + O2 <=> NO + OH */
    phi_f = sc[9]*sc[3];
    k_f = 1e-06 * 1.28e+06*exp(1.5*tc[0]-50.3271/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[4];
    Kc = exp((g_RT[9] + g_RT[3]) - (g_RT[13] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 46: NH + N <=> N2 + H */
    phi_f = sc[9]*sc[8];
    k_f = 1e-06 * 1.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[1];
    Kc = exp((g_RT[9] + g_RT[8]) - (g_RT[17] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[8] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 47: NH + H2O <=> HNO + H2 */
    phi_f = sc[9]*sc[5];
    k_f = 1e-06 * 2e+13*exp(-6970.31/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[0];
    Kc = exp((g_RT[9] + g_RT[5]) - (g_RT[16] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[16] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 48: NH + NO <=> N2 + OH */
    phi_f = sc[9]*sc[13];
    k_f = 1e-06 * 2.16e+13*exp(-0.23*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[4];
    Kc = exp((g_RT[9] + g_RT[13]) - (g_RT[17] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 49: NH + NO <=> N2O + H */
    phi_f = sc[9]*sc[13];
    k_f = 1e-06 * 4.16e+14*exp(-0.45*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[1];
    Kc = exp((g_RT[9] + g_RT[13]) - (g_RT[15] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[9] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[15] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 50: NH2 + O <=> OH + NH */
    phi_f = sc[10]*sc[2];
    k_f = 1e-06 * 7e+12;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = exp((g_RT[10] + g_RT[2]) - (g_RT[4] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[9] += 1 * qdot;

    /*reaction 51: NH2 + O <=> H + HNO */
    phi_f = sc[10]*sc[2];
    k_f = 1e-06 * 4.6e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = exp((g_RT[10] + g_RT[2]) - (g_RT[1] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 52: NH2 + H <=> NH + H2 */
    phi_f = sc[10]*sc[1];
    k_f = 1e-06 * 4e+13*exp(-1836.94/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = exp((g_RT[10] + g_RT[1]) - (g_RT[9] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 53: NH2 + OH <=> NH + H2O */
    phi_f = sc[10]*sc[4];
    k_f = 1e-06 * 9e+07*exp(1.5*tc[0]+231.505/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[5];
    Kc = exp((g_RT[10] + g_RT[4]) - (g_RT[9] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[10] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 54: NNH <=> N2 + H */
    phi_f = sc[12];
    k_f = 1 * 3.3e+08;
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[1];
    Kc = refC * exp((g_RT[12]) - (g_RT[17] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 55: NNH + M <=> N2 + H + M */
    phi_f = sc[12];
    alpha = mixture + sc[0] + 5*sc[5];
    k_f = 1e-06 * alpha * 1.3e+14*exp(-0.11*tc[0]-2506.29/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[1];
    Kc = refC * exp((g_RT[12]) - (g_RT[17] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[17] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 56: NNH + O2 <=> HO2 + N2 */
    phi_f = sc[12]*sc[3];
    k_f = 1e-06 * 5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = exp((g_RT[12] + g_RT[3]) - (g_RT[6] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 57: NNH + O <=> OH + N2 */
    phi_f = sc[12]*sc[2];
    k_f = 1e-06 * 2.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = exp((g_RT[12] + g_RT[2]) - (g_RT[4] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 58: NNH + O <=> NH + NO */
    phi_f = sc[12]*sc[2];
    k_f = 1e-06 * 7e+13;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[13];
    Kc = exp((g_RT[12] + g_RT[2]) - (g_RT[9] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[9] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 59: NNH + H <=> H2 + N2 */
    phi_f = sc[12]*sc[1];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[17];
    Kc = exp((g_RT[12] + g_RT[1]) - (g_RT[0] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 60: NNH + OH <=> H2O + N2 */
    phi_f = sc[12]*sc[4];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[17];
    Kc = exp((g_RT[12] + g_RT[4]) - (g_RT[5] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[12] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[17] += 1 * qdot;

    /*reaction 61: H + NO + M <=> HNO + M */
    phi_f = sc[1]*sc[13];
    alpha = mixture + sc[0] + 5*sc[5];
    k_f = 1e-12 * alpha * 8.95e+19*exp(-1.32*tc[0]-372.421/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[13]) - (g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 1 * qdot;
    wdot[13] -= 1 * qdot;
    wdot[16] += 1 * qdot;

    /*reaction 62: HNO + O <=> NO + OH */
    phi_f = sc[16]*sc[2];
    k_f = 1e-06 * 2.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[4];
    Kc = exp((g_RT[16] + g_RT[2]) - (g_RT[13] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[4] += 1 * qdot;

    /*reaction 63: HNO + H <=> H2 + NO */
    phi_f = sc[16]*sc[1];
    k_f = 1e-06 * 4.5e+11*exp(0.72*tc[0]-332.159/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[13];
    Kc = exp((g_RT[16] + g_RT[1]) - (g_RT[0] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 64: HNO + OH <=> NO + H2O */
    phi_f = sc[16]*sc[4];
    k_f = 1e-06 * 1.3e+07*exp(1.9*tc[0]+478.108/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[5];
    Kc = exp((g_RT[16] + g_RT[4]) - (g_RT[13] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[13] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 65: HNO + O2 <=> HO2 + NO */
    phi_f = sc[16]*sc[3];
    k_f = 1e-06 * 1e+13*exp(-6542.53/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = exp((g_RT[16] + g_RT[3]) - (g_RT[6] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[16] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[13] += 1 * qdot;

    /*reaction 66: NH3 + H <=> NH2 + H2 */
    phi_f = sc[11]*sc[1];
    k_f = 1e-06 * 540000*exp(2.4*tc[0]-4989.93/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[0];
    Kc = exp((g_RT[11] + g_RT[1]) - (g_RT[10] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 67: NH3 + OH <=> NH2 + H2O */
    phi_f = sc[11]*sc[4];
    k_f = 1e-06 * 5e+07*exp(1.6*tc[0]-480.624/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[5];
    Kc = exp((g_RT[11] + g_RT[4]) - (g_RT[10] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 68: NH3 + O <=> NH2 + OH */
    phi_f = sc[11]*sc[2];
    k_f = 1e-06 * 9.4e+06*exp(1.94*tc[0]-3251.13/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[4];
    Kc = exp((g_RT[11] + g_RT[2]) - (g_RT[10] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[11] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[10] += 1 * qdot;
    wdot[4] += 1 * qdot;

    return;
}


/*compute the progress rate for each reaction */
void progressRate(double * qdot, double * sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[18];                /*Gibbs free energy */
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
    for (id = 0; id < 18; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*reaction 1: 2 O + M <=> O2 + M */
    phi_f = sc[2]*sc[2];
    alpha = mixture + 1.4*sc[0] + 14.4*sc[5];
    k_f = 1e-12 * alpha * 1.2e+17*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[3];
    Kc = 1.0 / (refC) * exp((2 * g_RT[2]) - (g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[0] = q_f - q_r;

    /*reaction 2: O + H + M <=> OH + M */
    phi_f = sc[2]*sc[1];
    alpha = mixture + sc[0] + 5*sc[5];
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

    /*reaction 6: H + O2 + M <=> HO2 + M */
    phi_f = sc[1]*sc[3];
    alpha = mixture + -1*sc[3] + -1*sc[5] + -1*sc[17];
    k_f = 1e-12 * alpha * 2.8e+18*exp(-0.86*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[5] = q_f - q_r;

    /*reaction 7: H + 2 O2 <=> HO2 + O2 */
    phi_f = sc[1]*sc[3]*sc[3];
    k_f = 1e-12 * 3e+20*exp(-1.72*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = 1.0 / (refC) * exp((g_RT[1] + 2 * g_RT[3]) - (g_RT[6] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[6] = q_f - q_r;

    /*reaction 8: H + O2 + H2O <=> HO2 + H2O */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = 1e-12 * 9.38e+18*exp(-0.76*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[7] = q_f - q_r;

    /*reaction 9: H + O2 + N2 <=> HO2 + N2 */
    phi_f = sc[1]*sc[3]*sc[17];
    k_f = 1e-12 * 3.75e+20*exp(-1.72*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[17]) - (g_RT[6] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[8] = q_f - q_r;

    /*reaction 10: H + O2 <=> O + OH */
    phi_f = sc[1]*sc[3];
    k_f = 1e-06 * 8.3e+13*exp(-7253.65/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[9] = q_f - q_r;

    /*reaction 11: 2 H + M <=> H2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[0] + -1*sc[5];
    k_f = 1e-12 * alpha * 1e+18*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[10] = q_f - q_r;

    /*reaction 12: 2 H + H2 <=> 2 H2 */
    phi_f = sc[1]*sc[1]*sc[0];
    k_f = 1e-12 * 9e+16*exp(-0.6*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[0];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[11] = q_f - q_r;

    /*reaction 13: 2 H + H2O <=> H2 + H2O */
    phi_f = sc[1]*sc[1]*sc[5];
    k_f = 1e-12 * 6e+19*exp(-1.25*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[5];
    Kc = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[12] = q_f - q_r;

    /*reaction 14: H + OH + M <=> H2O + M */
    phi_f = sc[1]*sc[4];
    alpha = mixture + -0.27*sc[0] + 2.65*sc[5];
    k_f = 1e-12 * alpha * 2.2e+22*exp(-2*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[13] = q_f - q_r;

    /*reaction 15: H + HO2 <=> O + H2O */
    phi_f = sc[1]*sc[6];
    k_f = 1e-06 * 3.97e+12*exp(-337.695/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = exp((g_RT[1] + g_RT[6]) - (g_RT[2] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[14] = q_f - q_r;

    /*reaction 16: H + HO2 <=> O2 + H2 */
    phi_f = sc[1]*sc[6];
    k_f = 1e-06 * 2.8e+13*exp(-537.494/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[15] = q_f - q_r;

    /*reaction 17: H + HO2 <=> 2 OH */
    phi_f = sc[1]*sc[6];
    k_f = 1e-06 * 1.34e+14*exp(-319.577/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[16] = q_f - q_r;

    /*reaction 18: H + H2O2 <=> HO2 + H2 */
    phi_f = sc[1]*sc[7];
    k_f = 1e-06 * 1.21e+07*exp(2*tc[0]-2617.01/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = exp((g_RT[1] + g_RT[7]) - (g_RT[6] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[17] = q_f - q_r;

    /*reaction 19: H + H2O2 <=> OH + H2O */
    phi_f = sc[1]*sc[7];
    k_f = 1e-06 * 1e+13*exp(-1811.78/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[5];
    Kc = exp((g_RT[1] + g_RT[7]) - (g_RT[4] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[18] = q_f - q_r;

    /*reaction 20: OH + H2 <=> H + H2O */
    phi_f = sc[4]*sc[0];
    k_f = 1e-06 * 2.16e+08*exp(1.51*tc[0]-1726.22/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    Kc = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[19] = q_f - q_r;

    /*reaction 21: 2 OH (+M) <=> H2O2 (+M) */
    phi_f = sc[4]*sc[4];
    alpha = mixture + sc[0] + 5*sc[5];
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
    qdot[20] = q_f - q_r;

    /*reaction 22: 2 OH <=> O + H2O */
    phi_f = sc[4]*sc[4];
    k_f = 1e-06 * 35700*exp(2.4*tc[0]+1061.9/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[21] = q_f - q_r;

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = 1e-06 * 2.9e+13*exp(+251.636/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[22] = q_f - q_r;

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    phi_f = sc[4]*sc[7];
    k_f = 1e-06 * 1.75e+12*exp(-161.047/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[23] = q_f - q_r;

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    phi_f = sc[4]*sc[7];
    k_f = 1e-06 * 5.8e+14*exp(-4811.27/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[24] = q_f - q_r;

    /*reaction 26: 2 HO2 <=> O2 + H2O2 */
    phi_f = sc[6]*sc[6];
    k_f = 1e-06 * 1.3e+11*exp(+820.332/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[25] = q_f - q_r;

    /*reaction 27: 2 HO2 <=> O2 + H2O2 */
    phi_f = sc[6]*sc[6];
    k_f = 1e-06 * 4.2e+14*exp(-6039.26/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[26] = q_f - q_r;

    /*reaction 28: N + NO <=> N2 + O */
    phi_f = sc[8]*sc[13];
    k_f = 1e-06 * 3.5e+13*exp(-166.08/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[2];
    Kc = exp((g_RT[8] + g_RT[13]) - (g_RT[17] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[27] = q_f - q_r;

    /*reaction 29: N + O2 <=> NO + O */
    phi_f = sc[8]*sc[3];
    k_f = 1e-06 * 2.65e+12*exp(-3220.94/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[2];
    Kc = exp((g_RT[8] + g_RT[3]) - (g_RT[13] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[28] = q_f - q_r;

    /*reaction 30: N + OH <=> NO + H */
    phi_f = sc[8]*sc[4];
    k_f = 1e-06 * 7.333e+13*exp(-563.664/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[1];
    Kc = exp((g_RT[8] + g_RT[4]) - (g_RT[13] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[29] = q_f - q_r;

    /*reaction 31: N2O + O <=> N2 + O2 */
    phi_f = sc[15]*sc[2];
    k_f = 1e-06 * 1.4e+12*exp(-5440.36/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[3];
    Kc = exp((g_RT[15] + g_RT[2]) - (g_RT[17] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[30] = q_f - q_r;

    /*reaction 32: N2O + O <=> 2 NO */
    phi_f = sc[15]*sc[2];
    k_f = 1e-06 * 2.9e+13*exp(-11650.7/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[13];
    Kc = exp((g_RT[15] + g_RT[2]) - (2 * g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[31] = q_f - q_r;

    /*reaction 33: N2O + H <=> N2 + OH */
    phi_f = sc[15]*sc[1];
    k_f = 1e-06 * 4.4e+14*exp(-9501.76/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[4];
    Kc = exp((g_RT[15] + g_RT[1]) - (g_RT[17] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[32] = q_f - q_r;

    /*reaction 34: N2O + OH <=> N2 + HO2 */
    phi_f = sc[15]*sc[4];
    k_f = 1e-06 * 2e+12*exp(-10598.9/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[6];
    Kc = exp((g_RT[15] + g_RT[4]) - (g_RT[17] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[33] = q_f - q_r;

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    phi_f = sc[15];
    alpha = mixture + sc[0] + 5*sc[5];
    k_f = 1 * 1.3e+11*exp(-30005/tc[1]);
    redP = 1e-12 * alpha / k_f * 6.2e+14*exp(-28233.5/tc[1]);
    F = redP / (1 + redP);
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[2];
    Kc = refC * exp((g_RT[15]) - (g_RT[17] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[34] = q_f - q_r;

    /*reaction 36: HO2 + NO <=> NO2 + OH */
    phi_f = sc[6]*sc[13];
    k_f = 1e-06 * 2.11e+12*exp(+241.57/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[14]*sc[4];
    Kc = exp((g_RT[6] + g_RT[13]) - (g_RT[14] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[35] = q_f - q_r;

    /*reaction 37: NO + O + M <=> NO2 + M */
    phi_f = sc[13]*sc[2];
    alpha = mixture + sc[0] + 5*sc[5];
    k_f = 1e-12 * alpha * 1.06e+20*exp(-1.41*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[14];
    Kc = 1.0 / (refC) * exp((g_RT[13] + g_RT[2]) - (g_RT[14]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[36] = q_f - q_r;

    /*reaction 38: NO2 + O <=> NO + O2 */
    phi_f = sc[14]*sc[2];
    k_f = 1e-06 * 3.9e+12*exp(+120.785/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[3];
    Kc = exp((g_RT[14] + g_RT[2]) - (g_RT[13] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[37] = q_f - q_r;

    /*reaction 39: NO2 + H <=> NO + OH */
    phi_f = sc[14]*sc[1];
    k_f = 1e-06 * 1.32e+14*exp(-181.178/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[4];
    Kc = exp((g_RT[14] + g_RT[1]) - (g_RT[13] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[38] = q_f - q_r;

    /*reaction 40: NH + O <=> NO + H */
    phi_f = sc[9]*sc[2];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[1];
    Kc = exp((g_RT[9] + g_RT[2]) - (g_RT[13] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[39] = q_f - q_r;

    /*reaction 41: NH + H <=> N + H2 */
    phi_f = sc[9]*sc[1];
    k_f = 1e-06 * 3.2e+13*exp(-166.08/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[0];
    Kc = exp((g_RT[9] + g_RT[1]) - (g_RT[8] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[40] = q_f - q_r;

    /*reaction 42: NH + OH <=> HNO + H */
    phi_f = sc[9]*sc[4];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[1];
    Kc = exp((g_RT[9] + g_RT[4]) - (g_RT[16] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[41] = q_f - q_r;

    /*reaction 43: NH + OH <=> N + H2O */
    phi_f = sc[9]*sc[4];
    k_f = 1e-06 * 2e+09*exp(1.2*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[8]*sc[5];
    Kc = exp((g_RT[9] + g_RT[4]) - (g_RT[8] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[42] = q_f - q_r;

    /*reaction 44: NH + O2 <=> HNO + O */
    phi_f = sc[9]*sc[3];
    k_f = 1e-06 * 461000*exp(2*tc[0]-3271.26/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[2];
    Kc = exp((g_RT[9] + g_RT[3]) - (g_RT[16] + g_RT[2]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[43] = q_f - q_r;

    /*reaction 45: NH + O2 <=> NO + OH */
    phi_f = sc[9]*sc[3];
    k_f = 1e-06 * 1.28e+06*exp(1.5*tc[0]-50.3271/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[4];
    Kc = exp((g_RT[9] + g_RT[3]) - (g_RT[13] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[44] = q_f - q_r;

    /*reaction 46: NH + N <=> N2 + H */
    phi_f = sc[9]*sc[8];
    k_f = 1e-06 * 1.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[1];
    Kc = exp((g_RT[9] + g_RT[8]) - (g_RT[17] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[45] = q_f - q_r;

    /*reaction 47: NH + H2O <=> HNO + H2 */
    phi_f = sc[9]*sc[5];
    k_f = 1e-06 * 2e+13*exp(-6970.31/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16]*sc[0];
    Kc = exp((g_RT[9] + g_RT[5]) - (g_RT[16] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[46] = q_f - q_r;

    /*reaction 48: NH + NO <=> N2 + OH */
    phi_f = sc[9]*sc[13];
    k_f = 1e-06 * 2.16e+13*exp(-0.23*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[4];
    Kc = exp((g_RT[9] + g_RT[13]) - (g_RT[17] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[47] = q_f - q_r;

    /*reaction 49: NH + NO <=> N2O + H */
    phi_f = sc[9]*sc[13];
    k_f = 1e-06 * 4.16e+14*exp(-0.45*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[15]*sc[1];
    Kc = exp((g_RT[9] + g_RT[13]) - (g_RT[15] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[48] = q_f - q_r;

    /*reaction 50: NH2 + O <=> OH + NH */
    phi_f = sc[10]*sc[2];
    k_f = 1e-06 * 7e+12;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[9];
    Kc = exp((g_RT[10] + g_RT[2]) - (g_RT[4] + g_RT[9]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[49] = q_f - q_r;

    /*reaction 51: NH2 + O <=> H + HNO */
    phi_f = sc[10]*sc[2];
    k_f = 1e-06 * 4.6e+13;
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[16];
    Kc = exp((g_RT[10] + g_RT[2]) - (g_RT[1] + g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[50] = q_f - q_r;

    /*reaction 52: NH2 + H <=> NH + H2 */
    phi_f = sc[10]*sc[1];
    k_f = 1e-06 * 4e+13*exp(-1836.94/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[0];
    Kc = exp((g_RT[10] + g_RT[1]) - (g_RT[9] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[51] = q_f - q_r;

    /*reaction 53: NH2 + OH <=> NH + H2O */
    phi_f = sc[10]*sc[4];
    k_f = 1e-06 * 9e+07*exp(1.5*tc[0]+231.505/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[5];
    Kc = exp((g_RT[10] + g_RT[4]) - (g_RT[9] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[52] = q_f - q_r;

    /*reaction 54: NNH <=> N2 + H */
    phi_f = sc[12];
    k_f = 1 * 3.3e+08;
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[1];
    Kc = refC * exp((g_RT[12]) - (g_RT[17] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[53] = q_f - q_r;

    /*reaction 55: NNH + M <=> N2 + H + M */
    phi_f = sc[12];
    alpha = mixture + sc[0] + 5*sc[5];
    k_f = 1e-06 * alpha * 1.3e+14*exp(-0.11*tc[0]-2506.29/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[17]*sc[1];
    Kc = refC * exp((g_RT[12]) - (g_RT[17] + g_RT[1]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[54] = q_f - q_r;

    /*reaction 56: NNH + O2 <=> HO2 + N2 */
    phi_f = sc[12]*sc[3];
    k_f = 1e-06 * 5e+12;
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[17];
    Kc = exp((g_RT[12] + g_RT[3]) - (g_RT[6] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[55] = q_f - q_r;

    /*reaction 57: NNH + O <=> OH + N2 */
    phi_f = sc[12]*sc[2];
    k_f = 1e-06 * 2.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[17];
    Kc = exp((g_RT[12] + g_RT[2]) - (g_RT[4] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[56] = q_f - q_r;

    /*reaction 58: NNH + O <=> NH + NO */
    phi_f = sc[12]*sc[2];
    k_f = 1e-06 * 7e+13;
    q_f = phi_f * k_f;
    phi_r = sc[9]*sc[13];
    Kc = exp((g_RT[12] + g_RT[2]) - (g_RT[9] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[57] = q_f - q_r;

    /*reaction 59: NNH + H <=> H2 + N2 */
    phi_f = sc[12]*sc[1];
    k_f = 1e-06 * 5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[17];
    Kc = exp((g_RT[12] + g_RT[1]) - (g_RT[0] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[58] = q_f - q_r;

    /*reaction 60: NNH + OH <=> H2O + N2 */
    phi_f = sc[12]*sc[4];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[17];
    Kc = exp((g_RT[12] + g_RT[4]) - (g_RT[5] + g_RT[17]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[59] = q_f - q_r;

    /*reaction 61: H + NO + M <=> HNO + M */
    phi_f = sc[1]*sc[13];
    alpha = mixture + sc[0] + 5*sc[5];
    k_f = 1e-12 * alpha * 8.95e+19*exp(-1.32*tc[0]-372.421/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[16];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[13]) - (g_RT[16]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[60] = q_f - q_r;

    /*reaction 62: HNO + O <=> NO + OH */
    phi_f = sc[16]*sc[2];
    k_f = 1e-06 * 2.5e+13;
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[4];
    Kc = exp((g_RT[16] + g_RT[2]) - (g_RT[13] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[61] = q_f - q_r;

    /*reaction 63: HNO + H <=> H2 + NO */
    phi_f = sc[16]*sc[1];
    k_f = 1e-06 * 4.5e+11*exp(0.72*tc[0]-332.159/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[13];
    Kc = exp((g_RT[16] + g_RT[1]) - (g_RT[0] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[62] = q_f - q_r;

    /*reaction 64: HNO + OH <=> NO + H2O */
    phi_f = sc[16]*sc[4];
    k_f = 1e-06 * 1.3e+07*exp(1.9*tc[0]+478.108/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[13]*sc[5];
    Kc = exp((g_RT[16] + g_RT[4]) - (g_RT[13] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[63] = q_f - q_r;

    /*reaction 65: HNO + O2 <=> HO2 + NO */
    phi_f = sc[16]*sc[3];
    k_f = 1e-06 * 1e+13*exp(-6542.53/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[13];
    Kc = exp((g_RT[16] + g_RT[3]) - (g_RT[6] + g_RT[13]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[64] = q_f - q_r;

    /*reaction 66: NH3 + H <=> NH2 + H2 */
    phi_f = sc[11]*sc[1];
    k_f = 1e-06 * 540000*exp(2.4*tc[0]-4989.93/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[0];
    Kc = exp((g_RT[11] + g_RT[1]) - (g_RT[10] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[65] = q_f - q_r;

    /*reaction 67: NH3 + OH <=> NH2 + H2O */
    phi_f = sc[11]*sc[4];
    k_f = 1e-06 * 5e+07*exp(1.6*tc[0]-480.624/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[5];
    Kc = exp((g_RT[11] + g_RT[4]) - (g_RT[10] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[66] = q_f - q_r;

    /*reaction 68: NH3 + O <=> NH2 + OH */
    phi_f = sc[11]*sc[2];
    k_f = 1e-06 * 9.4e+06*exp(1.94*tc[0]-3251.13/tc[1]);
    q_f = phi_f * k_f;
    phi_r = sc[10]*sc[4];
    Kc = exp((g_RT[11] + g_RT[2]) - (g_RT[10] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[67] = q_f - q_r;

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

    /*reaction 6: H + O2 + M <=> HO2 + M */
    kc[5] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));

    /*reaction 7: H + 2 O2 <=> HO2 + O2 */
    kc[6] = 1.0 / (refC) * exp((g_RT[1] + 2 * g_RT[3]) - (g_RT[6] + g_RT[3]));

    /*reaction 8: H + O2 + H2O <=> HO2 + H2O */
    kc[7] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));

    /*reaction 9: H + O2 + N2 <=> HO2 + N2 */
    kc[8] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[17]) - (g_RT[6] + g_RT[17]));

    /*reaction 10: H + O2 <=> O + OH */
    kc[9] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));

    /*reaction 11: 2 H + M <=> H2 + M */
    kc[10] = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));

    /*reaction 12: 2 H + H2 <=> 2 H2 */
    kc[11] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));

    /*reaction 13: 2 H + H2O <=> H2 + H2O */
    kc[12] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));

    /*reaction 14: H + OH + M <=> H2O + M */
    kc[13] = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));

    /*reaction 15: H + HO2 <=> O + H2O */
    kc[14] = exp((g_RT[1] + g_RT[6]) - (g_RT[2] + g_RT[5]));

    /*reaction 16: H + HO2 <=> O2 + H2 */
    kc[15] = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));

    /*reaction 17: H + HO2 <=> 2 OH */
    kc[16] = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));

    /*reaction 18: H + H2O2 <=> HO2 + H2 */
    kc[17] = exp((g_RT[1] + g_RT[7]) - (g_RT[6] + g_RT[0]));

    /*reaction 19: H + H2O2 <=> OH + H2O */
    kc[18] = exp((g_RT[1] + g_RT[7]) - (g_RT[4] + g_RT[5]));

    /*reaction 20: OH + H2 <=> H + H2O */
    kc[19] = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));

    /*reaction 21: 2 OH (+M) <=> H2O2 (+M) */
    kc[20] = 1.0 / (refC) * exp((2 * g_RT[4]) - (g_RT[7]));

    /*reaction 22: 2 OH <=> O + H2O */
    kc[21] = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    kc[22] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    kc[23] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    kc[24] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));

    /*reaction 26: 2 HO2 <=> O2 + H2O2 */
    kc[25] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));

    /*reaction 27: 2 HO2 <=> O2 + H2O2 */
    kc[26] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));

    /*reaction 28: N + NO <=> N2 + O */
    kc[27] = exp((g_RT[8] + g_RT[13]) - (g_RT[17] + g_RT[2]));

    /*reaction 29: N + O2 <=> NO + O */
    kc[28] = exp((g_RT[8] + g_RT[3]) - (g_RT[13] + g_RT[2]));

    /*reaction 30: N + OH <=> NO + H */
    kc[29] = exp((g_RT[8] + g_RT[4]) - (g_RT[13] + g_RT[1]));

    /*reaction 31: N2O + O <=> N2 + O2 */
    kc[30] = exp((g_RT[15] + g_RT[2]) - (g_RT[17] + g_RT[3]));

    /*reaction 32: N2O + O <=> 2 NO */
    kc[31] = exp((g_RT[15] + g_RT[2]) - (2 * g_RT[13]));

    /*reaction 33: N2O + H <=> N2 + OH */
    kc[32] = exp((g_RT[15] + g_RT[1]) - (g_RT[17] + g_RT[4]));

    /*reaction 34: N2O + OH <=> N2 + HO2 */
    kc[33] = exp((g_RT[15] + g_RT[4]) - (g_RT[17] + g_RT[6]));

    /*reaction 35: N2O (+M) <=> N2 + O (+M) */
    kc[34] = refC * exp((g_RT[15]) - (g_RT[17] + g_RT[2]));

    /*reaction 36: HO2 + NO <=> NO2 + OH */
    kc[35] = exp((g_RT[6] + g_RT[13]) - (g_RT[14] + g_RT[4]));

    /*reaction 37: NO + O + M <=> NO2 + M */
    kc[36] = 1.0 / (refC) * exp((g_RT[13] + g_RT[2]) - (g_RT[14]));

    /*reaction 38: NO2 + O <=> NO + O2 */
    kc[37] = exp((g_RT[14] + g_RT[2]) - (g_RT[13] + g_RT[3]));

    /*reaction 39: NO2 + H <=> NO + OH */
    kc[38] = exp((g_RT[14] + g_RT[1]) - (g_RT[13] + g_RT[4]));

    /*reaction 40: NH + O <=> NO + H */
    kc[39] = exp((g_RT[9] + g_RT[2]) - (g_RT[13] + g_RT[1]));

    /*reaction 41: NH + H <=> N + H2 */
    kc[40] = exp((g_RT[9] + g_RT[1]) - (g_RT[8] + g_RT[0]));

    /*reaction 42: NH + OH <=> HNO + H */
    kc[41] = exp((g_RT[9] + g_RT[4]) - (g_RT[16] + g_RT[1]));

    /*reaction 43: NH + OH <=> N + H2O */
    kc[42] = exp((g_RT[9] + g_RT[4]) - (g_RT[8] + g_RT[5]));

    /*reaction 44: NH + O2 <=> HNO + O */
    kc[43] = exp((g_RT[9] + g_RT[3]) - (g_RT[16] + g_RT[2]));

    /*reaction 45: NH + O2 <=> NO + OH */
    kc[44] = exp((g_RT[9] + g_RT[3]) - (g_RT[13] + g_RT[4]));

    /*reaction 46: NH + N <=> N2 + H */
    kc[45] = exp((g_RT[9] + g_RT[8]) - (g_RT[17] + g_RT[1]));

    /*reaction 47: NH + H2O <=> HNO + H2 */
    kc[46] = exp((g_RT[9] + g_RT[5]) - (g_RT[16] + g_RT[0]));

    /*reaction 48: NH + NO <=> N2 + OH */
    kc[47] = exp((g_RT[9] + g_RT[13]) - (g_RT[17] + g_RT[4]));

    /*reaction 49: NH + NO <=> N2O + H */
    kc[48] = exp((g_RT[9] + g_RT[13]) - (g_RT[15] + g_RT[1]));

    /*reaction 50: NH2 + O <=> OH + NH */
    kc[49] = exp((g_RT[10] + g_RT[2]) - (g_RT[4] + g_RT[9]));

    /*reaction 51: NH2 + O <=> H + HNO */
    kc[50] = exp((g_RT[10] + g_RT[2]) - (g_RT[1] + g_RT[16]));

    /*reaction 52: NH2 + H <=> NH + H2 */
    kc[51] = exp((g_RT[10] + g_RT[1]) - (g_RT[9] + g_RT[0]));

    /*reaction 53: NH2 + OH <=> NH + H2O */
    kc[52] = exp((g_RT[10] + g_RT[4]) - (g_RT[9] + g_RT[5]));

    /*reaction 54: NNH <=> N2 + H */
    kc[53] = refC * exp((g_RT[12]) - (g_RT[17] + g_RT[1]));

    /*reaction 55: NNH + M <=> N2 + H + M */
    kc[54] = refC * exp((g_RT[12]) - (g_RT[17] + g_RT[1]));

    /*reaction 56: NNH + O2 <=> HO2 + N2 */
    kc[55] = exp((g_RT[12] + g_RT[3]) - (g_RT[6] + g_RT[17]));

    /*reaction 57: NNH + O <=> OH + N2 */
    kc[56] = exp((g_RT[12] + g_RT[2]) - (g_RT[4] + g_RT[17]));

    /*reaction 58: NNH + O <=> NH + NO */
    kc[57] = exp((g_RT[12] + g_RT[2]) - (g_RT[9] + g_RT[13]));

    /*reaction 59: NNH + H <=> H2 + N2 */
    kc[58] = exp((g_RT[12] + g_RT[1]) - (g_RT[0] + g_RT[17]));

    /*reaction 60: NNH + OH <=> H2O + N2 */
    kc[59] = exp((g_RT[12] + g_RT[4]) - (g_RT[5] + g_RT[17]));

    /*reaction 61: H + NO + M <=> HNO + M */
    kc[60] = 1.0 / (refC) * exp((g_RT[1] + g_RT[13]) - (g_RT[16]));

    /*reaction 62: HNO + O <=> NO + OH */
    kc[61] = exp((g_RT[16] + g_RT[2]) - (g_RT[13] + g_RT[4]));

    /*reaction 63: HNO + H <=> H2 + NO */
    kc[62] = exp((g_RT[16] + g_RT[1]) - (g_RT[0] + g_RT[13]));

    /*reaction 64: HNO + OH <=> NO + H2O */
    kc[63] = exp((g_RT[16] + g_RT[4]) - (g_RT[13] + g_RT[5]));

    /*reaction 65: HNO + O2 <=> HO2 + NO */
    kc[64] = exp((g_RT[16] + g_RT[3]) - (g_RT[6] + g_RT[13]));

    /*reaction 66: NH3 + H <=> NH2 + H2 */
    kc[65] = exp((g_RT[11] + g_RT[1]) - (g_RT[10] + g_RT[0]));

    /*reaction 67: NH3 + OH <=> NH2 + H2O */
    kc[66] = exp((g_RT[11] + g_RT[4]) - (g_RT[10] + g_RT[5]));

    /*reaction 68: NH3 + O <=> NH2 + OH */
    kc[67] = exp((g_RT[11] + g_RT[2]) - (g_RT[10] + g_RT[4]));

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
        /*species 8: N */
        species[8] =
            +5.61046370e+04 / tc[1]
            -1.69390870e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 9: NH */
        species[9] =
            +4.18806290e+04 / tc[1]
            +1.64458070e+00
            -3.49290850e+00 * tc[0]
            -1.55895990e-04 * tc[1]
            +2.48174733e-07 * tc[2]
            -2.06803683e-10 * tc[3]
            +5.17848350e-14 * tc[4];
        /*species 10: NH2 */
        species[10] =
            +2.18859100e+04 / tc[1]
            +4.34584538e+00
            -4.20400290e+00 * tc[0]
            +1.05306925e-03 * tc[1]
            -1.18447247e-06 * tc[2]
            +4.67626642e-10 * tc[3]
            -8.22035850e-14 * tc[4];
        /*species 11: NH3 */
        species[11] =
            -6.74172850e+03 / tc[1]
            +4.91140017e+00
            -4.28602740e+00 * tc[0]
            +2.33026150e-03 * tc[1]
            -3.61975217e-06 * tc[2]
            +1.90074058e-09 * tc[3]
            -4.13190230e-13 * tc[4];
        /*species 12: NNH */
        species[12] =
            +2.87919730e+04 / tc[1]
            +1.36675170e+00
            -4.34469270e+00 * tc[0]
            +2.42485360e-03 * tc[1]
            -3.34324317e-06 * tc[2]
            +1.81053867e-09 * tc[3]
            -3.97347695e-13 * tc[4];
        /*species 13: NO */
        species[13] =
            +9.84462300e+03 / tc[1]
            +1.93762990e+00
            -4.21847630e+00 * tc[0]
            +2.31948800e-03 * tc[1]
            -1.84017033e-06 * tc[2]
            +7.78011283e-10 * tc[3]
            -1.40178850e-13 * tc[4];
        /*species 14: NO2 */
        species[14] =
            +2.89661790e+03 / tc[1]
            -2.36796050e+00
            -3.94403120e+00 * tc[0]
            +7.92714500e-04 * tc[1]
            -2.77630200e-06 * tc[2]
            +1.70628550e-09 * tc[3]
            -3.91752820e-13 * tc[4];
        /*species 15: N2O */
        species[15] =
            +8.74177440e+03 / tc[1]
            -8.50084180e+00
            -2.25715020e+00 * tc[0]
            -5.65236400e-03 * tc[1]
            +2.27855317e-06 * tc[2]
            -8.06831717e-10 * tc[3]
            +1.46535910e-13 * tc[4];
        /*species 16: HNO */
        species[16] =
            +1.15482970e+04 / tc[1]
            +2.78364990e+00
            -4.53349160e+00 * tc[0]
            +2.83480855e-03 * tc[1]
            -3.07886783e-06 * tc[2]
            +1.42809117e-09 * tc[3]
            -2.77272865e-13 * tc[4];
        /*species 17: N2 */
        species[17] =
            -1.02089990e+03 / tc[1]
            -6.51695000e-01
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
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
        /*species 8: N */
        species[8] =
            +5.61337730e+04 / tc[1]
            -2.23366670e+00
            -2.41594290e+00 * tc[0]
            -8.74453250e-05 * tc[1]
            +1.98372817e-08 * tc[2]
            -2.51885375e-12 * tc[3]
            +1.01804910e-16 * tc[4];
        /*species 9: NH */
        species[9] =
            +4.21208480e+04 / tc[1]
            -2.95708710e+00
            -2.78369280e+00 * tc[0]
            -6.64921500e-04 * tc[1]
            +7.07967450e-08 * tc[2]
            -6.52904175e-12 * tc[3]
            +2.75222350e-16 * tc[4];
        /*species 10: NH2 */
        species[10] =
            +2.21719570e+04 / tc[1]
            -3.68567420e+00
            -2.83474210e+00 * tc[0]
            -1.60365410e-03 * tc[1]
            +1.55651340e-07 * tc[2]
            -1.14191275e-11 * tc[3]
            +3.96030720e-16 * tc[4];
        /*species 11: NH3 */
        species[11] =
            -6.54469580e+03 / tc[1]
            -3.93184070e+00
            -2.63445210e+00 * tc[0]
            -2.83312800e-03 * tc[1]
            +2.87977933e-07 * tc[2]
            -1.98893008e-11 * tc[3]
            +6.28939300e-16 * tc[4];
        /*species 12: NNH */
        species[12] =
            +2.86506970e+04 / tc[1]
            -7.03752300e-01
            -3.76675440e+00 * tc[0]
            -1.44575410e-03 * tc[1]
            +1.73610333e-07 * tc[2]
            -1.40354950e-11 * tc[3]
            +5.04594800e-16 * tc[4];
        /*species 13: NO */
        species[13] =
            +9.92097460e+03 / tc[1]
            -3.10869710e+00
            -3.26060560e+00 * tc[0]
            -5.95552150e-04 * tc[1]
            +7.15284133e-08 * tc[2]
            -5.78813908e-12 * tc[3]
            +2.01680495e-16 * tc[4];
        /*species 14: NO2 */
        species[14] =
            +2.31649830e+03 / tc[1]
            +5.00217115e+00
            -4.88475420e+00 * tc[0]
            -1.08619780e-03 * tc[1]
            +1.38011510e-07 * tc[2]
            -1.31229250e-11 * tc[3]
            +5.25544750e-16 * tc[4];
        /*species 15: N2O */
        species[15] =
            +8.07340480e+03 / tc[1]
            +7.02479360e+00
            -4.82307290e+00 * tc[0]
            -1.31351255e-03 * tc[1]
            +1.59751457e-07 * tc[2]
            -1.33339267e-11 * tc[3]
            +4.88761515e-16 * tc[4];
        /*species 16: HNO */
        species[16] =
            +1.17505820e+04 / tc[1]
            -5.62712190e+00
            -2.97925090e+00 * tc[0]
            -1.74720295e-03 * tc[1]
            +1.30916297e-07 * tc[2]
            -4.78996617e-12 * tc[3]
            +9.66795800e-18 * tc[4];
        /*species 17: N2 */
        species[17] =
            -9.22797700e+02 / tc[1]
            -3.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
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
        /*species 8: N */
        species[8] =
            +5.61046370e+04 / tc[1]
            -2.69390870e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 9: NH */
        species[9] =
            +4.18806290e+04 / tc[1]
            +6.44580700e-01
            -3.49290850e+00 * tc[0]
            -1.55895990e-04 * tc[1]
            +2.48174733e-07 * tc[2]
            -2.06803683e-10 * tc[3]
            +5.17848350e-14 * tc[4];
        /*species 10: NH2 */
        species[10] =
            +2.18859100e+04 / tc[1]
            +3.34584538e+00
            -4.20400290e+00 * tc[0]
            +1.05306925e-03 * tc[1]
            -1.18447247e-06 * tc[2]
            +4.67626642e-10 * tc[3]
            -8.22035850e-14 * tc[4];
        /*species 11: NH3 */
        species[11] =
            -6.74172850e+03 / tc[1]
            +3.91140017e+00
            -4.28602740e+00 * tc[0]
            +2.33026150e-03 * tc[1]
            -3.61975217e-06 * tc[2]
            +1.90074058e-09 * tc[3]
            -4.13190230e-13 * tc[4];
        /*species 12: NNH */
        species[12] =
            +2.87919730e+04 / tc[1]
            +3.66751700e-01
            -4.34469270e+00 * tc[0]
            +2.42485360e-03 * tc[1]
            -3.34324317e-06 * tc[2]
            +1.81053867e-09 * tc[3]
            -3.97347695e-13 * tc[4];
        /*species 13: NO */
        species[13] =
            +9.84462300e+03 / tc[1]
            +9.37629900e-01
            -4.21847630e+00 * tc[0]
            +2.31948800e-03 * tc[1]
            -1.84017033e-06 * tc[2]
            +7.78011283e-10 * tc[3]
            -1.40178850e-13 * tc[4];
        /*species 14: NO2 */
        species[14] =
            +2.89661790e+03 / tc[1]
            -3.36796050e+00
            -3.94403120e+00 * tc[0]
            +7.92714500e-04 * tc[1]
            -2.77630200e-06 * tc[2]
            +1.70628550e-09 * tc[3]
            -3.91752820e-13 * tc[4];
        /*species 15: N2O */
        species[15] =
            +8.74177440e+03 / tc[1]
            -9.50084180e+00
            -2.25715020e+00 * tc[0]
            -5.65236400e-03 * tc[1]
            +2.27855317e-06 * tc[2]
            -8.06831717e-10 * tc[3]
            +1.46535910e-13 * tc[4];
        /*species 16: HNO */
        species[16] =
            +1.15482970e+04 / tc[1]
            +1.78364990e+00
            -4.53349160e+00 * tc[0]
            +2.83480855e-03 * tc[1]
            -3.07886783e-06 * tc[2]
            +1.42809117e-09 * tc[3]
            -2.77272865e-13 * tc[4];
        /*species 17: N2 */
        species[17] =
            -1.02089990e+03 / tc[1]
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
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
        /*species 8: N */
        species[8] =
            +5.61337730e+04 / tc[1]
            -3.23366670e+00
            -2.41594290e+00 * tc[0]
            -8.74453250e-05 * tc[1]
            +1.98372817e-08 * tc[2]
            -2.51885375e-12 * tc[3]
            +1.01804910e-16 * tc[4];
        /*species 9: NH */
        species[9] =
            +4.21208480e+04 / tc[1]
            -3.95708710e+00
            -2.78369280e+00 * tc[0]
            -6.64921500e-04 * tc[1]
            +7.07967450e-08 * tc[2]
            -6.52904175e-12 * tc[3]
            +2.75222350e-16 * tc[4];
        /*species 10: NH2 */
        species[10] =
            +2.21719570e+04 / tc[1]
            -4.68567420e+00
            -2.83474210e+00 * tc[0]
            -1.60365410e-03 * tc[1]
            +1.55651340e-07 * tc[2]
            -1.14191275e-11 * tc[3]
            +3.96030720e-16 * tc[4];
        /*species 11: NH3 */
        species[11] =
            -6.54469580e+03 / tc[1]
            -4.93184070e+00
            -2.63445210e+00 * tc[0]
            -2.83312800e-03 * tc[1]
            +2.87977933e-07 * tc[2]
            -1.98893008e-11 * tc[3]
            +6.28939300e-16 * tc[4];
        /*species 12: NNH */
        species[12] =
            +2.86506970e+04 / tc[1]
            -1.70375230e+00
            -3.76675440e+00 * tc[0]
            -1.44575410e-03 * tc[1]
            +1.73610333e-07 * tc[2]
            -1.40354950e-11 * tc[3]
            +5.04594800e-16 * tc[4];
        /*species 13: NO */
        species[13] =
            +9.92097460e+03 / tc[1]
            -4.10869710e+00
            -3.26060560e+00 * tc[0]
            -5.95552150e-04 * tc[1]
            +7.15284133e-08 * tc[2]
            -5.78813908e-12 * tc[3]
            +2.01680495e-16 * tc[4];
        /*species 14: NO2 */
        species[14] =
            +2.31649830e+03 / tc[1]
            +4.00217115e+00
            -4.88475420e+00 * tc[0]
            -1.08619780e-03 * tc[1]
            +1.38011510e-07 * tc[2]
            -1.31229250e-11 * tc[3]
            +5.25544750e-16 * tc[4];
        /*species 15: N2O */
        species[15] =
            +8.07340480e+03 / tc[1]
            +6.02479360e+00
            -4.82307290e+00 * tc[0]
            -1.31351255e-03 * tc[1]
            +1.59751457e-07 * tc[2]
            -1.33339267e-11 * tc[3]
            +4.88761515e-16 * tc[4];
        /*species 16: HNO */
        species[16] =
            +1.17505820e+04 / tc[1]
            -6.62712190e+00
            -2.97925090e+00 * tc[0]
            -1.74720295e-03 * tc[1]
            +1.30916297e-07 * tc[2]
            -4.78996617e-12 * tc[3]
            +9.66795800e-18 * tc[4];
        /*species 17: N2 */
        species[17] =
            -9.22797700e+02 / tc[1]
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
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
        /*species 8: N */
        species[8] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 9: NH */
        species[9] =
            +2.49290850e+00
            +3.11791980e-04 * tc[1]
            -1.48904840e-06 * tc[2]
            +2.48164420e-09 * tc[3]
            -1.03569670e-12 * tc[4];
        /*species 10: NH2 */
        species[10] =
            +3.20400290e+00
            -2.10613850e-03 * tc[1]
            +7.10683480e-06 * tc[2]
            -5.61151970e-09 * tc[3]
            +1.64407170e-12 * tc[4];
        /*species 11: NH3 */
        species[11] =
            +3.28602740e+00
            -4.66052300e-03 * tc[1]
            +2.17185130e-05 * tc[2]
            -2.28088870e-08 * tc[3]
            +8.26380460e-12 * tc[4];
        /*species 12: NNH */
        species[12] =
            +3.34469270e+00
            -4.84970720e-03 * tc[1]
            +2.00594590e-05 * tc[2]
            -2.17264640e-08 * tc[3]
            +7.94695390e-12 * tc[4];
        /*species 13: NO */
        species[13] =
            +3.21847630e+00
            -4.63897600e-03 * tc[1]
            +1.10410220e-05 * tc[2]
            -9.33613540e-09 * tc[3]
            +2.80357700e-12 * tc[4];
        /*species 14: NO2 */
        species[14] =
            +2.94403120e+00
            -1.58542900e-03 * tc[1]
            +1.66578120e-05 * tc[2]
            -2.04754260e-08 * tc[3]
            +7.83505640e-12 * tc[4];
        /*species 15: N2O */
        species[15] =
            +1.25715020e+00
            +1.13047280e-02 * tc[1]
            -1.36713190e-05 * tc[2]
            +9.68198060e-09 * tc[3]
            -2.93071820e-12 * tc[4];
        /*species 16: HNO */
        species[16] =
            +3.53349160e+00
            -5.66961710e-03 * tc[1]
            +1.84732070e-05 * tc[2]
            -1.71370940e-08 * tc[3]
            +5.54545730e-12 * tc[4];
        /*species 17: N2 */
        species[17] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
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
        /*species 8: N */
        species[8] =
            +1.41594290e+00
            +1.74890650e-04 * tc[1]
            -1.19023690e-07 * tc[2]
            +3.02262450e-11 * tc[3]
            -2.03609820e-15 * tc[4];
        /*species 9: NH */
        species[9] =
            +1.78369280e+00
            +1.32984300e-03 * tc[1]
            -4.24780470e-07 * tc[2]
            +7.83485010e-11 * tc[3]
            -5.50444700e-15 * tc[4];
        /*species 10: NH2 */
        species[10] =
            +1.83474210e+00
            +3.20730820e-03 * tc[1]
            -9.33908040e-07 * tc[2]
            +1.37029530e-10 * tc[3]
            -7.92061440e-15 * tc[4];
        /*species 11: NH3 */
        species[11] =
            +1.63445210e+00
            +5.66625600e-03 * tc[1]
            -1.72786760e-06 * tc[2]
            +2.38671610e-10 * tc[3]
            -1.25787860e-14 * tc[4];
        /*species 12: NNH */
        species[12] =
            +2.76675440e+00
            +2.89150820e-03 * tc[1]
            -1.04166200e-06 * tc[2]
            +1.68425940e-10 * tc[3]
            -1.00918960e-14 * tc[4];
        /*species 13: NO */
        species[13] =
            +2.26060560e+00
            +1.19110430e-03 * tc[1]
            -4.29170480e-07 * tc[2]
            +6.94576690e-11 * tc[3]
            -4.03360990e-15 * tc[4];
        /*species 14: NO2 */
        species[14] =
            +3.88475420e+00
            +2.17239560e-03 * tc[1]
            -8.28069060e-07 * tc[2]
            +1.57475100e-10 * tc[3]
            -1.05108950e-14 * tc[4];
        /*species 15: N2O */
        species[15] =
            +3.82307290e+00
            +2.62702510e-03 * tc[1]
            -9.58508740e-07 * tc[2]
            +1.60007120e-10 * tc[3]
            -9.77523030e-15 * tc[4];
        /*species 16: HNO */
        species[16] =
            +1.97925090e+00
            +3.49440590e-03 * tc[1]
            -7.85497780e-07 * tc[2]
            +5.74795940e-11 * tc[3]
            -1.93359160e-16 * tc[4];
        /*species 17: N2 */
        species[17] =
            +1.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
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
        /*species 8: N */
        species[8] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 9: NH */
        species[9] =
            +3.49290850e+00
            +3.11791980e-04 * tc[1]
            -1.48904840e-06 * tc[2]
            +2.48164420e-09 * tc[3]
            -1.03569670e-12 * tc[4];
        /*species 10: NH2 */
        species[10] =
            +4.20400290e+00
            -2.10613850e-03 * tc[1]
            +7.10683480e-06 * tc[2]
            -5.61151970e-09 * tc[3]
            +1.64407170e-12 * tc[4];
        /*species 11: NH3 */
        species[11] =
            +4.28602740e+00
            -4.66052300e-03 * tc[1]
            +2.17185130e-05 * tc[2]
            -2.28088870e-08 * tc[3]
            +8.26380460e-12 * tc[4];
        /*species 12: NNH */
        species[12] =
            +4.34469270e+00
            -4.84970720e-03 * tc[1]
            +2.00594590e-05 * tc[2]
            -2.17264640e-08 * tc[3]
            +7.94695390e-12 * tc[4];
        /*species 13: NO */
        species[13] =
            +4.21847630e+00
            -4.63897600e-03 * tc[1]
            +1.10410220e-05 * tc[2]
            -9.33613540e-09 * tc[3]
            +2.80357700e-12 * tc[4];
        /*species 14: NO2 */
        species[14] =
            +3.94403120e+00
            -1.58542900e-03 * tc[1]
            +1.66578120e-05 * tc[2]
            -2.04754260e-08 * tc[3]
            +7.83505640e-12 * tc[4];
        /*species 15: N2O */
        species[15] =
            +2.25715020e+00
            +1.13047280e-02 * tc[1]
            -1.36713190e-05 * tc[2]
            +9.68198060e-09 * tc[3]
            -2.93071820e-12 * tc[4];
        /*species 16: HNO */
        species[16] =
            +4.53349160e+00
            -5.66961710e-03 * tc[1]
            +1.84732070e-05 * tc[2]
            -1.71370940e-08 * tc[3]
            +5.54545730e-12 * tc[4];
        /*species 17: N2 */
        species[17] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
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
        /*species 8: N */
        species[8] =
            +2.41594290e+00
            +1.74890650e-04 * tc[1]
            -1.19023690e-07 * tc[2]
            +3.02262450e-11 * tc[3]
            -2.03609820e-15 * tc[4];
        /*species 9: NH */
        species[9] =
            +2.78369280e+00
            +1.32984300e-03 * tc[1]
            -4.24780470e-07 * tc[2]
            +7.83485010e-11 * tc[3]
            -5.50444700e-15 * tc[4];
        /*species 10: NH2 */
        species[10] =
            +2.83474210e+00
            +3.20730820e-03 * tc[1]
            -9.33908040e-07 * tc[2]
            +1.37029530e-10 * tc[3]
            -7.92061440e-15 * tc[4];
        /*species 11: NH3 */
        species[11] =
            +2.63445210e+00
            +5.66625600e-03 * tc[1]
            -1.72786760e-06 * tc[2]
            +2.38671610e-10 * tc[3]
            -1.25787860e-14 * tc[4];
        /*species 12: NNH */
        species[12] =
            +3.76675440e+00
            +2.89150820e-03 * tc[1]
            -1.04166200e-06 * tc[2]
            +1.68425940e-10 * tc[3]
            -1.00918960e-14 * tc[4];
        /*species 13: NO */
        species[13] =
            +3.26060560e+00
            +1.19110430e-03 * tc[1]
            -4.29170480e-07 * tc[2]
            +6.94576690e-11 * tc[3]
            -4.03360990e-15 * tc[4];
        /*species 14: NO2 */
        species[14] =
            +4.88475420e+00
            +2.17239560e-03 * tc[1]
            -8.28069060e-07 * tc[2]
            +1.57475100e-10 * tc[3]
            -1.05108950e-14 * tc[4];
        /*species 15: N2O */
        species[15] =
            +4.82307290e+00
            +2.62702510e-03 * tc[1]
            -9.58508740e-07 * tc[2]
            +1.60007120e-10 * tc[3]
            -9.77523030e-15 * tc[4];
        /*species 16: HNO */
        species[16] =
            +2.97925090e+00
            +3.49440590e-03 * tc[1]
            -7.85497780e-07 * tc[2]
            +5.74795940e-11 * tc[3]
            -1.93359160e-16 * tc[4];
        /*species 17: N2 */
        species[17] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
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
        /*species 8: N */
        species[8] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +5.61046370e+04 / tc[1];
        /*species 9: NH */
        species[9] =
            +2.49290850e+00
            +1.55895990e-04 * tc[1]
            -4.96349467e-07 * tc[2]
            +6.20411050e-10 * tc[3]
            -2.07139340e-13 * tc[4]
            +4.18806290e+04 / tc[1];
        /*species 10: NH2 */
        species[10] =
            +3.20400290e+00
            -1.05306925e-03 * tc[1]
            +2.36894493e-06 * tc[2]
            -1.40287992e-09 * tc[3]
            +3.28814340e-13 * tc[4]
            +2.18859100e+04 / tc[1];
        /*species 11: NH3 */
        species[11] =
            +3.28602740e+00
            -2.33026150e-03 * tc[1]
            +7.23950433e-06 * tc[2]
            -5.70222175e-09 * tc[3]
            +1.65276092e-12 * tc[4]
            -6.74172850e+03 / tc[1];
        /*species 12: NNH */
        species[12] =
            +3.34469270e+00
            -2.42485360e-03 * tc[1]
            +6.68648633e-06 * tc[2]
            -5.43161600e-09 * tc[3]
            +1.58939078e-12 * tc[4]
            +2.87919730e+04 / tc[1];
        /*species 13: NO */
        species[13] =
            +3.21847630e+00
            -2.31948800e-03 * tc[1]
            +3.68034067e-06 * tc[2]
            -2.33403385e-09 * tc[3]
            +5.60715400e-13 * tc[4]
            +9.84462300e+03 / tc[1];
        /*species 14: NO2 */
        species[14] =
            +2.94403120e+00
            -7.92714500e-04 * tc[1]
            +5.55260400e-06 * tc[2]
            -5.11885650e-09 * tc[3]
            +1.56701128e-12 * tc[4]
            +2.89661790e+03 / tc[1];
        /*species 15: N2O */
        species[15] =
            +1.25715020e+00
            +5.65236400e-03 * tc[1]
            -4.55710633e-06 * tc[2]
            +2.42049515e-09 * tc[3]
            -5.86143640e-13 * tc[4]
            +8.74177440e+03 / tc[1];
        /*species 16: HNO */
        species[16] =
            +3.53349160e+00
            -2.83480855e-03 * tc[1]
            +6.15773567e-06 * tc[2]
            -4.28427350e-09 * tc[3]
            +1.10909146e-12 * tc[4]
            +1.15482970e+04 / tc[1];
        /*species 17: N2 */
        species[17] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 / tc[1];
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
        /*species 8: N */
        species[8] =
            +1.41594290e+00
            +8.74453250e-05 * tc[1]
            -3.96745633e-08 * tc[2]
            +7.55656125e-12 * tc[3]
            -4.07219640e-16 * tc[4]
            +5.61337730e+04 / tc[1];
        /*species 9: NH */
        species[9] =
            +1.78369280e+00
            +6.64921500e-04 * tc[1]
            -1.41593490e-07 * tc[2]
            +1.95871253e-11 * tc[3]
            -1.10088940e-15 * tc[4]
            +4.21208480e+04 / tc[1];
        /*species 10: NH2 */
        species[10] =
            +1.83474210e+00
            +1.60365410e-03 * tc[1]
            -3.11302680e-07 * tc[2]
            +3.42573825e-11 * tc[3]
            -1.58412288e-15 * tc[4]
            +2.21719570e+04 / tc[1];
        /*species 11: NH3 */
        species[11] =
            +1.63445210e+00
            +2.83312800e-03 * tc[1]
            -5.75955867e-07 * tc[2]
            +5.96679025e-11 * tc[3]
            -2.51575720e-15 * tc[4]
            -6.54469580e+03 / tc[1];
        /*species 12: NNH */
        species[12] =
            +2.76675440e+00
            +1.44575410e-03 * tc[1]
            -3.47220667e-07 * tc[2]
            +4.21064850e-11 * tc[3]
            -2.01837920e-15 * tc[4]
            +2.86506970e+04 / tc[1];
        /*species 13: NO */
        species[13] =
            +2.26060560e+00
            +5.95552150e-04 * tc[1]
            -1.43056827e-07 * tc[2]
            +1.73644173e-11 * tc[3]
            -8.06721980e-16 * tc[4]
            +9.92097460e+03 / tc[1];
        /*species 14: NO2 */
        species[14] =
            +3.88475420e+00
            +1.08619780e-03 * tc[1]
            -2.76023020e-07 * tc[2]
            +3.93687750e-11 * tc[3]
            -2.10217900e-15 * tc[4]
            +2.31649830e+03 / tc[1];
        /*species 15: N2O */
        species[15] =
            +3.82307290e+00
            +1.31351255e-03 * tc[1]
            -3.19502913e-07 * tc[2]
            +4.00017800e-11 * tc[3]
            -1.95504606e-15 * tc[4]
            +8.07340480e+03 / tc[1];
        /*species 16: HNO */
        species[16] =
            +1.97925090e+00
            +1.74720295e-03 * tc[1]
            -2.61832593e-07 * tc[2]
            +1.43698985e-11 * tc[3]
            -3.86718320e-17 * tc[4]
            +1.17505820e+04 / tc[1];
        /*species 17: N2 */
        species[17] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 / tc[1];
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
        /*species 8: N */
        species[8] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +5.61046370e+04 / tc[1];
        /*species 9: NH */
        species[9] =
            +3.49290850e+00
            +1.55895990e-04 * tc[1]
            -4.96349467e-07 * tc[2]
            +6.20411050e-10 * tc[3]
            -2.07139340e-13 * tc[4]
            +4.18806290e+04 / tc[1];
        /*species 10: NH2 */
        species[10] =
            +4.20400290e+00
            -1.05306925e-03 * tc[1]
            +2.36894493e-06 * tc[2]
            -1.40287992e-09 * tc[3]
            +3.28814340e-13 * tc[4]
            +2.18859100e+04 / tc[1];
        /*species 11: NH3 */
        species[11] =
            +4.28602740e+00
            -2.33026150e-03 * tc[1]
            +7.23950433e-06 * tc[2]
            -5.70222175e-09 * tc[3]
            +1.65276092e-12 * tc[4]
            -6.74172850e+03 / tc[1];
        /*species 12: NNH */
        species[12] =
            +4.34469270e+00
            -2.42485360e-03 * tc[1]
            +6.68648633e-06 * tc[2]
            -5.43161600e-09 * tc[3]
            +1.58939078e-12 * tc[4]
            +2.87919730e+04 / tc[1];
        /*species 13: NO */
        species[13] =
            +4.21847630e+00
            -2.31948800e-03 * tc[1]
            +3.68034067e-06 * tc[2]
            -2.33403385e-09 * tc[3]
            +5.60715400e-13 * tc[4]
            +9.84462300e+03 / tc[1];
        /*species 14: NO2 */
        species[14] =
            +3.94403120e+00
            -7.92714500e-04 * tc[1]
            +5.55260400e-06 * tc[2]
            -5.11885650e-09 * tc[3]
            +1.56701128e-12 * tc[4]
            +2.89661790e+03 / tc[1];
        /*species 15: N2O */
        species[15] =
            +2.25715020e+00
            +5.65236400e-03 * tc[1]
            -4.55710633e-06 * tc[2]
            +2.42049515e-09 * tc[3]
            -5.86143640e-13 * tc[4]
            +8.74177440e+03 / tc[1];
        /*species 16: HNO */
        species[16] =
            +4.53349160e+00
            -2.83480855e-03 * tc[1]
            +6.15773567e-06 * tc[2]
            -4.28427350e-09 * tc[3]
            +1.10909146e-12 * tc[4]
            +1.15482970e+04 / tc[1];
        /*species 17: N2 */
        species[17] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 / tc[1];
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
        /*species 8: N */
        species[8] =
            +2.41594290e+00
            +8.74453250e-05 * tc[1]
            -3.96745633e-08 * tc[2]
            +7.55656125e-12 * tc[3]
            -4.07219640e-16 * tc[4]
            +5.61337730e+04 / tc[1];
        /*species 9: NH */
        species[9] =
            +2.78369280e+00
            +6.64921500e-04 * tc[1]
            -1.41593490e-07 * tc[2]
            +1.95871253e-11 * tc[3]
            -1.10088940e-15 * tc[4]
            +4.21208480e+04 / tc[1];
        /*species 10: NH2 */
        species[10] =
            +2.83474210e+00
            +1.60365410e-03 * tc[1]
            -3.11302680e-07 * tc[2]
            +3.42573825e-11 * tc[3]
            -1.58412288e-15 * tc[4]
            +2.21719570e+04 / tc[1];
        /*species 11: NH3 */
        species[11] =
            +2.63445210e+00
            +2.83312800e-03 * tc[1]
            -5.75955867e-07 * tc[2]
            +5.96679025e-11 * tc[3]
            -2.51575720e-15 * tc[4]
            -6.54469580e+03 / tc[1];
        /*species 12: NNH */
        species[12] =
            +3.76675440e+00
            +1.44575410e-03 * tc[1]
            -3.47220667e-07 * tc[2]
            +4.21064850e-11 * tc[3]
            -2.01837920e-15 * tc[4]
            +2.86506970e+04 / tc[1];
        /*species 13: NO */
        species[13] =
            +3.26060560e+00
            +5.95552150e-04 * tc[1]
            -1.43056827e-07 * tc[2]
            +1.73644173e-11 * tc[3]
            -8.06721980e-16 * tc[4]
            +9.92097460e+03 / tc[1];
        /*species 14: NO2 */
        species[14] =
            +4.88475420e+00
            +1.08619780e-03 * tc[1]
            -2.76023020e-07 * tc[2]
            +3.93687750e-11 * tc[3]
            -2.10217900e-15 * tc[4]
            +2.31649830e+03 / tc[1];
        /*species 15: N2O */
        species[15] =
            +4.82307290e+00
            +1.31351255e-03 * tc[1]
            -3.19502913e-07 * tc[2]
            +4.00017800e-11 * tc[3]
            -1.95504606e-15 * tc[4]
            +8.07340480e+03 / tc[1];
        /*species 16: HNO */
        species[16] =
            +2.97925090e+00
            +1.74720295e-03 * tc[1]
            -2.61832593e-07 * tc[2]
            +1.43698985e-11 * tc[3]
            -3.86718320e-17 * tc[4]
            +1.17505820e+04 / tc[1];
        /*species 17: N2 */
        species[17] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 / tc[1];
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
        /*species 8: N */
        species[8] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.19390870e+00 ;
        /*species 9: NH */
        species[9] =
            +3.49290850e+00 * tc[0]
            +3.11791980e-04 * tc[1]
            -7.44524200e-07 * tc[2]
            +8.27214733e-10 * tc[3]
            -2.58924175e-13 * tc[4]
            +1.84832780e+00 ;
        /*species 10: NH2 */
        species[10] =
            +4.20400290e+00 * tc[0]
            -2.10613850e-03 * tc[1]
            +3.55341740e-06 * tc[2]
            -1.87050657e-09 * tc[3]
            +4.11017925e-13 * tc[4]
            -1.41842480e-01 ;
        /*species 11: NH3 */
        species[11] =
            +4.28602740e+00 * tc[0]
            -4.66052300e-03 * tc[1]
            +1.08592565e-05 * tc[2]
            -7.60296233e-09 * tc[3]
            +2.06595115e-12 * tc[4]
            -6.25372770e-01 ;
        /*species 12: NNH */
        species[12] =
            +4.34469270e+00 * tc[0]
            -4.84970720e-03 * tc[1]
            +1.00297295e-05 * tc[2]
            -7.24215467e-09 * tc[3]
            +1.98673848e-12 * tc[4]
            +2.97794100e+00 ;
        /*species 13: NO */
        species[13] =
            +4.21847630e+00 * tc[0]
            -4.63897600e-03 * tc[1]
            +5.52051100e-06 * tc[2]
            -3.11204513e-09 * tc[3]
            +7.00894250e-13 * tc[4]
            +2.28084640e+00 ;
        /*species 14: NO2 */
        species[14] =
            +3.94403120e+00 * tc[0]
            -1.58542900e-03 * tc[1]
            +8.32890600e-06 * tc[2]
            -6.82514200e-09 * tc[3]
            +1.95876410e-12 * tc[4]
            +6.31199170e+00 ;
        /*species 15: N2O */
        species[15] =
            +2.25715020e+00 * tc[0]
            +1.13047280e-02 * tc[1]
            -6.83565950e-06 * tc[2]
            +3.22732687e-09 * tc[3]
            -7.32679550e-13 * tc[4]
            +1.07579920e+01 ;
        /*species 16: HNO */
        species[16] =
            +4.53349160e+00 * tc[0]
            -5.66961710e-03 * tc[1]
            +9.23660350e-06 * tc[2]
            -5.71236467e-09 * tc[3]
            +1.38636433e-12 * tc[4]
            +1.74984170e+00 ;
        /*species 17: N2 */
        species[17] =
            +3.29867700e+00 * tc[0]
            +1.40824040e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213500e-13 * tc[4]
            +3.95037200e+00 ;
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
        /*species 8: N */
        species[8] =
            +2.41594290e+00 * tc[0]
            +1.74890650e-04 * tc[1]
            -5.95118450e-08 * tc[2]
            +1.00754150e-11 * tc[3]
            -5.09024550e-16 * tc[4]
            +4.64960960e+00 ;
        /*species 9: NH */
        species[9] =
            +2.78369280e+00 * tc[0]
            +1.32984300e-03 * tc[1]
            -2.12390235e-07 * tc[2]
            +2.61161670e-11 * tc[3]
            -1.37611175e-15 * tc[4]
            +5.74077990e+00 ;
        /*species 10: NH2 */
        species[10] =
            +2.83474210e+00 * tc[0]
            +3.20730820e-03 * tc[1]
            -4.66954020e-07 * tc[2]
            +4.56765100e-11 * tc[3]
            -1.98015360e-15 * tc[4]
            +6.52041630e+00 ;
        /*species 11: NH3 */
        species[11] =
            +2.63445210e+00 * tc[0]
            +5.66625600e-03 * tc[1]
            -8.63933800e-07 * tc[2]
            +7.95572033e-11 * tc[3]
            -3.14469650e-15 * tc[4]
            +6.56629280e+00 ;
        /*species 12: NNH */
        species[12] =
            +3.76675440e+00 * tc[0]
            +2.89150820e-03 * tc[1]
            -5.20831000e-07 * tc[2]
            +5.61419800e-11 * tc[3]
            -2.52297400e-15 * tc[4]
            +4.47050670e+00 ;
        /*species 13: NO */
        species[13] =
            +3.26060560e+00 * tc[0]
            +1.19110430e-03 * tc[1]
            -2.14585240e-07 * tc[2]
            +2.31525563e-11 * tc[3]
            -1.00840247e-15 * tc[4]
            +6.36930270e+00 ;
        /*species 14: NO2 */
        species[14] =
            +4.88475420e+00 * tc[0]
            +2.17239560e-03 * tc[1]
            -4.14034530e-07 * tc[2]
            +5.24917000e-11 * tc[3]
            -2.62772375e-15 * tc[4]
            -1.17416950e-01 ;
        /*species 15: N2O */
        species[15] =
            +4.82307290e+00 * tc[0]
            +2.62702510e-03 * tc[1]
            -4.79254370e-07 * tc[2]
            +5.33357067e-11 * tc[3]
            -2.44380757e-15 * tc[4]
            -2.20172070e+00 ;
        /*species 16: HNO */
        species[16] =
            +2.97925090e+00 * tc[0]
            +3.49440590e-03 * tc[1]
            -3.92748890e-07 * tc[2]
            +1.91598647e-11 * tc[3]
            -4.83397900e-17 * tc[4]
            +8.60637280e+00 ;
        /*species 17: N2 */
        species[17] =
            +2.92664000e+00 * tc[0]
            +1.48797680e-03 * tc[1]
            -2.84238000e-07 * tc[2]
            +3.36567933e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
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
    wt[8] = 14.006700; /*N */
    wt[9] = 15.014670; /*NH */
    wt[10] = 16.022640; /*NH2 */
    wt[11] = 17.030610; /*NH3 */
    wt[12] = 29.021370; /*NNH */
    wt[13] = 30.006100; /*NO */
    wt[14] = 46.005500; /*NO2 */
    wt[15] = 44.012800; /*N2O */
    wt[16] = 31.014070; /*HNO */
    wt[17] = 28.013400; /*N2 */

    return;
}


/*get temperature given internal energy in mass units and mass fracs */
int feeytt_(double * e, double * y, int * iwrk, double * rwrk, double * t)
{
    const int maxiter = 50;
    const double tol  = 0.001;
    double ein  = *e;
    double tmin = 300; // max lower bound for thermo def
    double tmax = 3500; // min upper bound for thermo def
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
    y[8] = phi[8]*14.006700;   XW += y[8]; /*N */
    y[9] = phi[9]*15.014670;   XW += y[9]; /*NH */
    y[10] = phi[10]*16.022640;   XW += y[10]; /*NH2 */
    y[11] = phi[11]*17.030610;   XW += y[11]; /*NH3 */
    y[12] = phi[12]*29.021370;   XW += y[12]; /*NNH */
    y[13] = phi[13]*30.006100;   XW += y[13]; /*NO */
    y[14] = phi[14]*46.005500;   XW += y[14]; /*NO2 */
    y[15] = phi[15]*44.012800;   XW += y[15]; /*N2O */
    y[16] = phi[16]*31.014070;   XW += y[16]; /*HNO */
    y[17] = phi[17]*28.013400;   XW += y[17]; /*N2 */
    for (id = 0; id < 18; ++id) {
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
    phi[8] = y[8]/ 1.40067000e-02; /*N (wt in kg) */
    phi[9] = y[9]/ 1.50146700e-02; /*NH (wt in kg) */
    phi[10] = y[10]/ 1.60226400e-02; /*NH2 (wt in kg) */
    phi[11] = y[11]/ 1.70306100e-02; /*NH3 (wt in kg) */
    phi[12] = y[12]/ 2.90213700e-02; /*NNH (wt in kg) */
    phi[13] = y[13]/ 3.00061000e-02; /*NO (wt in kg) */
    phi[14] = y[14]/ 4.60055000e-02; /*NO2 (wt in kg) */
    phi[15] = y[15]/ 4.40128000e-02; /*N2O (wt in kg) */
    phi[16] = y[16]/ 3.10140700e-02; /*HNO (wt in kg) */
    phi[17] = y[17]/ 2.80134000e-02; /*N2 (wt in kg) */

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
    y[8] = c[8] * 14.006700 / (*rho); 
    y[9] = c[9] * 15.014670 / (*rho); 
    y[10] = c[10] * 16.022640 / (*rho); 
    y[11] = c[11] * 17.030610 / (*rho); 
    y[12] = c[12] * 29.021370 / (*rho); 
    y[13] = c[13] * 30.006100 / (*rho); 
    y[14] = c[14] * 46.005500 / (*rho); 
    y[15] = c[15] * 44.012800 / (*rho); 
    y[16] = c[16] * 31.014070 / (*rho); 
    y[17] = c[17] * 28.013400 / (*rho); 

    return;
}


/*ddebdf compatible right hand side of CV burner */
/*rwrk[0] and rwrk[1] should contain rho and ene respectively */
/*working variable phi contains specific mole numbers */
void fecvrhs_(double * time, double * phi, double * phidot, double * rwrk, int * iwrk)
{
    double rho,ene; /*CV Parameters */
    double y[18], wdot[18]; /*temporary storage */
    int i; /*Loop counter */
    double temperature,pressure; /*temporary var */
    rho = rwrk[0];
    ene = rwrk[1];
    fephity_(phi, iwrk, rwrk, y);
    feeytt_(&ene, y, iwrk, rwrk, &temperature);
    CKPY(&rho, &temperature,  y, iwrk, rwrk, &pressure);
    CKWYP(&pressure, &temperature,  y, iwrk, rwrk, wdot);
    for (i=0; i<18; ++i) phidot[i] = wdot[i] / (rho/1000.0); 

    return;
}


/*returns the dimensionality of the cv burner (number of species) */
int fecvdim_()
{
    return 18;
}


/*ddebdf compatible right hand side of ZND solver */
/*rwrk[0] : scaling factor for pressure */
/*rwrk[1] : preshock density (g/cc)  */
/*rwrk[2] : detonation velocity (cm/s)  */
/*solution vector: [P; rho; y0 ... ylast]  */
void fezndrhs_(double * time, double * z, double * zdot, double * rwrk, int * iwrk)
{
    double psc,rho1,udet; /*ZND Parameters */
    double wt[18], hms[18], wdot[18]; /*temporary storage */
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
    for (i=0; i<18; ++i) {
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
    return 21;
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
    if (strcmp(s1, "N")==0) return 8; 
    if (strcmp(s1, "NH")==0) return 9; 
    if (strcmp(s1, "NH2")==0) return 10; 
    if (strcmp(s1, "NH3")==0) return 11; 
    if (strcmp(s1, "NNH")==0) return 12; 
    if (strcmp(s1, "NO")==0) return 13; 
    if (strcmp(s1, "NO2")==0) return 14; 
    if (strcmp(s1, "N2O")==0) return 15; 
    if (strcmp(s1, "HNO")==0) return 16; 
    if (strcmp(s1, "N2")==0) return 17; 
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
    if (sn==8) return "N"; 
    if (sn==9) return "NH"; 
    if (sn==10) return "NH2"; 
    if (sn==11) return "NH3"; 
    if (sn==12) return "NNH"; 
    if (sn==13) return "NO"; 
    if (sn==14) return "NO2"; 
    if (sn==15) return "N2O"; 
    if (sn==16) return "HNO"; 
    if (sn==17) return "N2"; 
    /*species name not found */
    return "NOTFOUND";
}
