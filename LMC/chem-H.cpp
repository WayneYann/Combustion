
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


/*A few mechanism parameters */
void CKINDX(int * iwrk, double * rwrk, int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 3;
    *kk = 9;
    *ii = 27;
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
    for (i=0; i<lenkname*9; i++) {
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

    /* N2  */
    kname[ 8*lenkname + 0 ] = 'N';
    kname[ 8*lenkname + 1 ] = '2';
    kname[ 8*lenkname + 2 ] = ' ';

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
    XW += x[8]*28.013400; /*N2 */
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
    YOW += y[8]/28.013400; /*N2 */
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
    W += c[8]*28.013400; /*N2 */

    for (id = 0; id < 9; ++id) {
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
    XW += x[8]*28.013400; /*N2 */
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
    YOW += y[8]/28.013400; /*N2 */
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
    W += c[8]*28.013400; /*N2 */

    for (id = 0; id < 9; ++id) {
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
    YOW += y[8]/28.013400; /*N2 */
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
    XW += x[8]*28.013400; /*N2 */
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
    W += c[8]*28.013400; /*N2 */

    for (id = 0; id < 9; ++id) {
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
    YOW += y[8]/28.013400; /*N2 */
    /*Now compute conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.013400*YOW); 

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
    YOW += y[8]/28.013400; /*N2 */
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
    c[8] = PWORT * y[8]/28.013400; 

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
    c[8] = (*rho) * y[8]/28.013400; 

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
    XW += x[8]*28.013400; /*N2 */
    /*Now compute conversion */
    y[0] = x[0]*2.015940/XW; 
    y[1] = x[1]*1.007970/XW; 
    y[2] = x[2]*15.999400/XW; 
    y[3] = x[3]*31.998800/XW; 
    y[4] = x[4]*17.007370/XW; 
    y[5] = x[5]*18.015340/XW; 
    y[6] = x[6]*33.006770/XW; 
    y[7] = x[7]*34.014740/XW; 
    y[8] = x[8]*28.013400/XW; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.314e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 9; ++id) {
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
    XW += x[8]*28.013400; /*N2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    for (id = 0; id < 9; ++id) {
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
    CW += c[8]*28.013400; /*N2 */
    /*Now compute conversion */
    y[0] = c[0]*2.015940/CW; 
    y[1] = c[1]*1.007970/CW; 
    y[2] = c[2]*15.999400/CW; 
    y[3] = c[3]*31.998800/CW; 
    y[4] = c[4]*17.007370/CW; 
    y[5] = c[5]*18.015340/CW; 
    y[6] = c[6]*33.006770/CW; 
    y[7] = c[7]*34.014740/CW; 
    y[8] = c[8]*28.013400/CW; 

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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    cvms[8] *= 2967865.378712; /*N2 */
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
    cpms[8] *= 2967865.378712; /*N2 */
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
    ums[8] *= RT/28.013400; /*N2 */
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
    hms[8] *= RT/28.013400; /*N2 */
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
    gms[8] *= RT/28.013400; /*N2 */
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
    ams[8] *= RT/28.013400; /*N2 */
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
    sms[8] *= 2967865.378712; /*N2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double *T, double *x, int * iwrk, double * rwrk, double * cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[9]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
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
    double cpor[9]; /* temporary storage */
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
    result += cpor[8]*y[8]/28.0134; /*N2 */

    *cpbs = result * 8.314e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double *T, double *x, int * iwrk, double * rwrk, double * cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[9]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
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
    double cvor[9]; /* temporary storage */
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
    result += cvor[8]*y[8]/28.0134; /*N2 */

    *cvbs = result * 8.314e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double *T, double *x, int * iwrk, double * rwrk, double * hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[9]; /* temporary storage */
    double RT = 8.314e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
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
    double hml[9]; /* temporary storage */
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
    result += y[8]*hml[8]/28.013400; /*N2 */

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double *T, double *x, int * iwrk, double * rwrk, double * ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[9]; /* temporary energy array */
    double RT = 8.314e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
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
    double ums[9]; /* temporary energy array */
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
    result += y[8]*ums[8]/28.013400; /*N2 */

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
    double sor[9]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 9; ++id) {
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
    double sor[9]; /* temporary storage */
    double x[9]; /* need a ytx conversion */
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
    YOW += y[8]/28.013400; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.013400*YOW); 
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
    double gort[9]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 9; ++id) {
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
    double gort[9]; /* temporary storage */
    double x[9]; /* need a ytx conversion */
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
    YOW += y[8]/28.013400; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.013400*YOW); 
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
    double aort[9]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 9; ++id) {
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
    double aort[9]; /* temporary storage */
    double x[9]; /* need a ytx conversion */
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
    YOW += y[8]/28.013400; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.013400*YOW); 
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
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double * T, double * C, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 9; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
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
    YOW += y[8]/28.013400; /*N2 */
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
    c[8] = PWORT * y[8]/28.013400; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.314e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]/2.015940; 
    c[1] = 1e6 * (*rho) * y[1]/1.007970; 
    c[2] = 1e6 * (*rho) * y[2]/15.999400; 
    c[3] = 1e6 * (*rho) * y[3]/31.998800; 
    c[4] = 1e6 * (*rho) * y[4]/17.007370; 
    c[5] = 1e6 * (*rho) * y[5]/18.015340; 
    c[6] = 1e6 * (*rho) * y[6]/33.006770; 
    c[7] = 1e6 * (*rho) * y[7]/34.014740; 
    c[8] = 1e6 * (*rho) * y[8]/28.013400; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
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
    XW += x[8]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double * T, double * C, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 9; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 27; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
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
    YOW += y[8]/28.013400; /*N2 */
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
    c[8] = PWORT * y[8]/28.013400; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 27; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.314e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 27; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]/2.015940; 
    c[1] = 1e6 * (*rho) * y[1]/1.007970; 
    c[2] = 1e6 * (*rho) * y[2]/15.999400; 
    c[3] = 1e6 * (*rho) * y[3]/31.998800; 
    c[4] = 1e6 * (*rho) * y[4]/17.007370; 
    c[5] = 1e6 * (*rho) * y[5]/18.015340; 
    c[6] = 1e6 * (*rho) * y[6]/33.006770; 
    c[7] = 1e6 * (*rho) * y[7]/34.014740; 
    c[8] = 1e6 * (*rho) * y[8]/28.013400; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 27; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
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
    XW += x[8]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 27; ++id) {
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
    for (id = 0; id < 9 * 27; ++ id) {
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
    nuki[ 8 * kd + 8 ] = -1 ;
    nuki[ 6 * kd + 8 ] = +1 ;
    nuki[ 8 * kd + 8 ] = +1 ;

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
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < 3 * 9; ++ id) {
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

    /*N2 */
    ncf[ 8 * kd + 2 ] = 2; /*N */

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

    return;
}



/*Returns the equil constants for each reaction */
void CKEQC(double * T, double * C, int * iwrk, double * rwrk, double * eqcon)
{
   double tT = *T; /*temporary temperature */
   double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
   double gort[9]; /* temporary storage */
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
}

/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * P, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
   double tT = *T; /*temporary temperature */
   double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
   double gort[9]; /* temporary storage */
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
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * P, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

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
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

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
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

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
}


/*compute the production rate for each species */
void productionRate(double * wdot, double * sc, double T)
{
    double qdot;

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[9];                /*Gibbs free energy */
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
    double refC  = 101325 / 8.314 / T;
    double irefC = 1.0 / refC;
    double Tinv  = 1.0 / T;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 9; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*zero out wdot */
    for (id = 0; id < 9; ++id) {
        wdot[id] = 0.0;
    }

    /*reaction 1: 2 O + M <=> O2 + M */
    phi_f = sc[2]*sc[2];
    alpha = mixture + 1.4*sc[0] + 14.4*sc[5];
    k_f = 1e-12 * alpha * 1.2e+17*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[3];
    Kc = irefC * exp((2 * g_RT[2]) - (g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= 2 * qdot;
    wdot[3] += qdot;

    /*reaction 2: O + H + M <=> OH + M */
    phi_f = sc[2]*sc[1];
    alpha = mixture + sc[0] + 5*sc[5];
    k_f = 1e-12 * alpha * 5e+17*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[4];
    Kc = irefC * exp((g_RT[2] + g_RT[1]) - (g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= qdot;
    wdot[1] -= qdot;
    wdot[4] += qdot;

    /*reaction 3: O + H2 <=> H + OH */
    phi_f = sc[2]*sc[0];
    k_f = 1e-06 * 50000*exp(2.67*tc[0]-3165.58*Tinv);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[4];
    Kc = exp((g_RT[2] + g_RT[0]) - (g_RT[1] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= qdot;
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[4] += qdot;

    /*reaction 4: O + HO2 <=> OH + O2 */
    phi_f = sc[2]*sc[6];
    k_f = 1e-06 * 2e+13;
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[3];
    Kc = exp((g_RT[2] + g_RT[6]) - (g_RT[4] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= qdot;
    wdot[6] -= qdot;
    wdot[4] += qdot;
    wdot[3] += qdot;

    /*reaction 5: O + H2O2 <=> OH + HO2 */
    phi_f = sc[2]*sc[7];
    k_f = 1e-06 * 9.63e+06*exp(2*tc[0]-2013.09*Tinv);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[6];
    Kc = exp((g_RT[2] + g_RT[7]) - (g_RT[4] + g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[2] -= qdot;
    wdot[7] -= qdot;
    wdot[4] += qdot;
    wdot[6] += qdot;

    /*reaction 6: H + O2 + M <=> HO2 + M */
    phi_f = sc[1]*sc[3];
    alpha = mixture + -1*sc[3] + -1*sc[5] + -1*sc[8];
    k_f = 1e-12 * alpha * 2.8e+18*exp(-0.86*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = irefC * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;

    /*reaction 7: H + 2 O2 <=> HO2 + O2 */
    phi_f = sc[1]*sc[3]*sc[3];
    k_f = 1e-12 * 3e+20*exp(-1.72*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[3];
    Kc = irefC * exp((g_RT[1] + 2 * g_RT[3]) - (g_RT[6] + g_RT[3]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= qdot;
    wdot[3] -= 2 * qdot;
    wdot[6] += qdot;
    wdot[3] += qdot;

    /*reaction 8: H + O2 + H2O <=> HO2 + H2O */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = 1e-12 * 9.38e+18*exp(-0.76*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = irefC * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[5] += qdot;

    /*reaction 9: H + O2 + N2 <=> HO2 + N2 */
    phi_f = sc[1]*sc[3]*sc[8];
    k_f = 1e-12 * 3.75e+20*exp(-1.72*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[8];
    Kc = irefC * exp((g_RT[1] + g_RT[3] + g_RT[8]) - (g_RT[6] + g_RT[8]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[8] -= qdot;
    wdot[6] += qdot;
    wdot[8] += qdot;

    /*reaction 10: H + O2 <=> O + OH */
    phi_f = sc[1]*sc[3];
    k_f = 1e-06 * 8.3e+13*exp(-7253.65*Tinv);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[4];
    Kc = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[2] += qdot;
    wdot[4] += qdot;

    /*reaction 11: 2 H + M <=> H2 + M */
    phi_f = sc[1]*sc[1];
    alpha = mixture + -1*sc[0] + -1*sc[5];
    k_f = 1e-12 * alpha * 1e+18*exp(-1*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0];
    Kc = irefC * exp((2 * g_RT[1]) - (g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[0] += qdot;

    /*reaction 12: 2 H + H2 <=> 2 H2 */
    phi_f = sc[1]*sc[1]*sc[0];
    k_f = 1e-12 * 9e+16*exp(-0.6*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[0];
    Kc = irefC * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[0] -= qdot;
    wdot[0] += 2 * qdot;

    /*reaction 13: 2 H + H2O <=> H2 + H2O */
    phi_f = sc[1]*sc[1]*sc[5];
    k_f = 1e-12 * 6e+19*exp(-1.25*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[5];
    Kc = irefC * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= 2 * qdot;
    wdot[5] -= qdot;
    wdot[0] += qdot;
    wdot[5] += qdot;

    /*reaction 14: H + OH + M <=> H2O + M */
    phi_f = sc[1]*sc[4];
    alpha = mixture + -0.27*sc[0] + 2.65*sc[5];
    k_f = 1e-12 * alpha * 2.2e+22*exp(-2*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = irefC * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    /*reaction 15: H + HO2 <=> O + H2O */
    phi_f = sc[1]*sc[6];
    k_f = 1e-06 * 3.97e+12*exp(-337.695*Tinv);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = exp((g_RT[1] + g_RT[6]) - (g_RT[2] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= qdot;
    wdot[6] -= qdot;
    wdot[2] += qdot;
    wdot[5] += qdot;

    /*reaction 16: H + HO2 <=> O2 + H2 */
    phi_f = sc[1]*sc[6];
    k_f = 1e-06 * 2.8e+13*exp(-537.494*Tinv);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[0];
    Kc = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= qdot;
    wdot[6] -= qdot;
    wdot[3] += qdot;
    wdot[0] += qdot;

    /*reaction 17: H + HO2 <=> 2 OH */
    phi_f = sc[1]*sc[6];
    k_f = 1e-06 * 1.34e+14*exp(-319.577*Tinv);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[4];
    Kc = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= qdot;
    wdot[6] -= qdot;
    wdot[4] += 2 * qdot;

    /*reaction 18: H + H2O2 <=> HO2 + H2 */
    phi_f = sc[1]*sc[7];
    k_f = 1e-06 * 1.21e+07*exp(2*tc[0]-2617.01*Tinv);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = exp((g_RT[1] + g_RT[7]) - (g_RT[6] + g_RT[0]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= qdot;
    wdot[7] -= qdot;
    wdot[6] += qdot;
    wdot[0] += qdot;

    /*reaction 19: H + H2O2 <=> OH + H2O */
    phi_f = sc[1]*sc[7];
    k_f = 1e-06 * 1e+13*exp(-1811.78*Tinv);
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[5];
    Kc = exp((g_RT[1] + g_RT[7]) - (g_RT[4] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[1] -= qdot;
    wdot[7] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;

    /*reaction 20: OH + H2 <=> H + H2O */
    phi_f = sc[4]*sc[0];
    k_f = 1e-06 * 2.16e+08*exp(1.51*tc[0]-1726.22*Tinv);
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    Kc = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= qdot;
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[5] += qdot;

    /*reaction 21: 2 OH (+M) <=> H2O2 (+M) */
    phi_f = sc[4]*sc[4];
    alpha = mixture + sc[0] + 5*sc[5];
    k_f = 1e-06 * 7.4e+13*exp(-0.37*tc[0]);
    redP = 1e-12 * alpha / k_f * 2.3e+18*exp(-0.9*tc[0]+855.561*Tinv);
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
    Kc = irefC * exp((2 * g_RT[4]) - (g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 2 * qdot;
    wdot[7] += qdot;

    /*reaction 22: 2 OH <=> O + H2O */
    phi_f = sc[4]*sc[4];
    k_f = 1e-06 * 35700*exp(2.4*tc[0]+1061.9*Tinv);
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 2 * qdot;
    wdot[2] += qdot;
    wdot[5] += qdot;

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    phi_f = sc[4]*sc[6];
    k_f = 1e-06 * 2.9e+13*exp(+251.636*Tinv);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= qdot;
    wdot[6] -= qdot;
    wdot[3] += qdot;
    wdot[5] += qdot;

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    phi_f = sc[4]*sc[7];
    k_f = 1e-06 * 1.75e+12*exp(-161.047*Tinv);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= qdot;
    wdot[7] -= qdot;
    wdot[6] += qdot;
    wdot[5] += qdot;

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    phi_f = sc[4]*sc[7];
    k_f = 1e-06 * 5.8e+14*exp(-4811.27*Tinv);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[5];
    Kc = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= qdot;
    wdot[7] -= qdot;
    wdot[6] += qdot;
    wdot[5] += qdot;

    /*reaction 26: 2 HO2 <=> O2 + H2O2 */
    phi_f = sc[6]*sc[6];
    k_f = 1e-06 * 1.3e+11*exp(+820.332*Tinv);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 2 * qdot;
    wdot[3] += qdot;
    wdot[7] += qdot;

    /*reaction 27: 2 HO2 <=> O2 + H2O2 */
    phi_f = sc[6]*sc[6];
    k_f = 1e-06 * 4.2e+14*exp(-6039.26*Tinv);
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[7];
    Kc = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 2 * qdot;
    wdot[3] += qdot;
    wdot[7] += qdot;

    return;
}


/*compute the progress rate for each reaction */
void progressRate(double * qdot, double * sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[9];                /*Gibbs free energy */
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
    for (id = 0; id < 9; ++id) {
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
    alpha = mixture + -1*sc[3] + -1*sc[5] + -1*sc[8];
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
    phi_f = sc[1]*sc[3]*sc[8];
    k_f = 1e-12 * 3.75e+20*exp(-1.72*tc[0]);
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[8];
    Kc = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[8]) - (g_RT[6] + g_RT[8]));
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
    kc[8] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[8]) - (g_RT[6] + g_RT[8]));

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

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void gibbs(double * species, double * tc)
{
    static double T_save = -1;
    static double species_save[9];

    if (T_save == tc[1])
    {
        for (int i = 0; i < 9; i++)
            species[i] = species_save[i];
        return;
    }

    /*temperature */
    double T    = tc[1];
    double Tinv = 1.0/T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            -1.01252100e+03 * Tinv
            +6.59221800e+00
            -3.29812400e+00 * tc[0]
            -4.12472100e-04 * tc[1]
            +1.35716917e-07 * tc[2]
            +7.89619500e-12 * tc[3]
            -2.06743600e-14 * tc[4];
        /*species 1: H */
        species[1] =
            +2.54716300e+04 * Tinv
            +2.96011760e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 2: O */
        species[2] =
            +2.91476400e+04 * Tinv
            -1.75660000e-02
            -2.94642900e+00 * tc[0]
            +8.19083000e-04 * tc[1]
            -4.03505333e-07 * tc[2]
            +1.33570250e-10 * tc[3]
            -1.94534800e-14 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.00524900e+03 * Tinv
            -2.82180200e+00
            -3.21293600e+00 * tc[0]
            -5.63743000e-04 * tc[1]
            +9.59358333e-08 * tc[2]
            -1.09489750e-10 * tc[3]
            +4.38427700e-14 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.60678200e+03 * Tinv
            +2.27840600e+00
            -3.63726600e+00 * tc[0]
            -9.25455000e-05 * tc[1]
            +2.79360833e-07 * tc[2]
            -1.98933583e-10 * tc[3]
            +4.21572100e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.02081100e+04 * Tinv
            +7.96609000e-01
            -3.38684200e+00 * tc[0]
            -1.73749100e-03 * tc[1]
            +1.05911600e-06 * tc[2]
            -5.80715083e-10 * tc[3]
            +1.25329400e-13 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.76227400e+02 * Tinv
            -6.24276100e+00
            -2.97996300e+00 * tc[0]
            -2.49834850e-03 * tc[1]
            +6.31832833e-07 * tc[2]
            -1.96182667e-10 * tc[3]
            +4.04451200e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.76631500e+04 * Tinv
            -3.39660900e+00
            -3.38875400e+00 * tc[0]
            -3.28461300e-03 * tc[1]
            +2.47502167e-08 * tc[2]
            +3.85483833e-10 * tc[3]
            -1.23575750e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.02090000e+03 * Tinv
            -6.51695000e-01
            -3.29867700e+00 * tc[0]
            -7.04120000e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242750e-13 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -8.35034000e+02 * Tinv
            +4.34653300e+00
            -2.99142300e+00 * tc[0]
            -3.50032200e-04 * tc[1]
            +9.38971500e-09 * tc[2]
            +7.69298167e-13 * tc[3]
            -7.91376000e-17 * tc[4];
        /*species 1: H */
        species[1] =
            +2.54716300e+04 * Tinv
            +2.96011760e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 2: O */
        species[2] =
            +2.92308000e+04 * Tinv
            -2.37824800e+00
            -2.54206000e+00 * tc[0]
            +1.37753100e-05 * tc[1]
            +5.17133833e-10 * tc[2]
            -3.79255583e-13 * tc[3]
            +2.18402600e-17 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.23393000e+03 * Tinv
            +5.08412000e-01
            -3.69757800e+00 * tc[0]
            -3.06759850e-04 * tc[1]
            +2.09807000e-08 * tc[2]
            -1.47940083e-12 * tc[3]
            +5.68217500e-17 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.88688800e+03 * Tinv
            -2.71298200e+00
            -2.88273000e+00 * tc[0]
            -5.06987000e-04 * tc[1]
            +3.79479500e-08 * tc[2]
            -1.81223667e-12 * tc[3]
            +2.56315250e-17 * tc[4];
        /*species 5: H2O */
        species[5] =
            -2.98992100e+04 * Tinv
            -4.19067100e+00
            -2.67214600e+00 * tc[0]
            -1.52814650e-03 * tc[1]
            +1.45504333e-07 * tc[2]
            -1.00083000e-11 * tc[3]
            +3.19580900e-16 * tc[4];
        /*species 6: HO2 */
        species[6] =
            -1.57972700e+02 * Tinv
            +5.96162000e-01
            -4.07219100e+00 * tc[0]
            -1.06564800e-03 * tc[1]
            +8.84690833e-08 * tc[2]
            -5.09355750e-12 * tc[3]
            +1.42058250e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.80069600e+04 * Tinv
            +4.07203000e+00
            -4.57316700e+00 * tc[0]
            -2.16806800e-03 * tc[1]
            +2.45781500e-07 * tc[2]
            -1.95742000e-11 * tc[3]
            +7.15827000e-16 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.22797700e+02 * Tinv
            -3.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988500e-04 * tc[1]
            +9.47460167e-08 * tc[2]
            -8.41420000e-12 * tc[3]
            +3.37667550e-16 * tc[4];
    }

    T_save = T;

    for (int i = 0; i < 9; i++)
      species_save[i] = species[i];

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
            -1.01252100e+03 / tc[1]
            +5.59221800e+00
            -3.29812400e+00 * tc[0]
            -4.12472100e-04 * tc[1]
            +1.35716917e-07 * tc[2]
            +7.89619500e-12 * tc[3]
            -2.06743600e-14 * tc[4];
        /*species 1: H */
        species[1] =
            +2.54716300e+04 / tc[1]
            +1.96011760e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 2: O */
        species[2] =
            +2.91476400e+04 / tc[1]
            -1.01756600e+00
            -2.94642900e+00 * tc[0]
            +8.19083000e-04 * tc[1]
            -4.03505333e-07 * tc[2]
            +1.33570250e-10 * tc[3]
            -1.94534800e-14 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.00524900e+03 / tc[1]
            -3.82180200e+00
            -3.21293600e+00 * tc[0]
            -5.63743000e-04 * tc[1]
            +9.59358333e-08 * tc[2]
            -1.09489750e-10 * tc[3]
            +4.38427700e-14 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.60678200e+03 / tc[1]
            +1.27840600e+00
            -3.63726600e+00 * tc[0]
            -9.25455000e-05 * tc[1]
            +2.79360833e-07 * tc[2]
            -1.98933583e-10 * tc[3]
            +4.21572100e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.02081100e+04 / tc[1]
            -2.03391000e-01
            -3.38684200e+00 * tc[0]
            -1.73749100e-03 * tc[1]
            +1.05911600e-06 * tc[2]
            -5.80715083e-10 * tc[3]
            +1.25329400e-13 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.76227400e+02 / tc[1]
            -7.24276100e+00
            -2.97996300e+00 * tc[0]
            -2.49834850e-03 * tc[1]
            +6.31832833e-07 * tc[2]
            -1.96182667e-10 * tc[3]
            +4.04451200e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.76631500e+04 / tc[1]
            -4.39660900e+00
            -3.38875400e+00 * tc[0]
            -3.28461300e-03 * tc[1]
            +2.47502167e-08 * tc[2]
            +3.85483833e-10 * tc[3]
            -1.23575750e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.02090000e+03 / tc[1]
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120000e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242750e-13 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -8.35034000e+02 / tc[1]
            +3.34653300e+00
            -2.99142300e+00 * tc[0]
            -3.50032200e-04 * tc[1]
            +9.38971500e-09 * tc[2]
            +7.69298167e-13 * tc[3]
            -7.91376000e-17 * tc[4];
        /*species 1: H */
        species[1] =
            +2.54716300e+04 / tc[1]
            +1.96011760e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 2: O */
        species[2] =
            +2.92308000e+04 / tc[1]
            -3.37824800e+00
            -2.54206000e+00 * tc[0]
            +1.37753100e-05 * tc[1]
            +5.17133833e-10 * tc[2]
            -3.79255583e-13 * tc[3]
            +2.18402600e-17 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.23393000e+03 / tc[1]
            -4.91588000e-01
            -3.69757800e+00 * tc[0]
            -3.06759850e-04 * tc[1]
            +2.09807000e-08 * tc[2]
            -1.47940083e-12 * tc[3]
            +5.68217500e-17 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.88688800e+03 / tc[1]
            -3.71298200e+00
            -2.88273000e+00 * tc[0]
            -5.06987000e-04 * tc[1]
            +3.79479500e-08 * tc[2]
            -1.81223667e-12 * tc[3]
            +2.56315250e-17 * tc[4];
        /*species 5: H2O */
        species[5] =
            -2.98992100e+04 / tc[1]
            -5.19067100e+00
            -2.67214600e+00 * tc[0]
            -1.52814650e-03 * tc[1]
            +1.45504333e-07 * tc[2]
            -1.00083000e-11 * tc[3]
            +3.19580900e-16 * tc[4];
        /*species 6: HO2 */
        species[6] =
            -1.57972700e+02 / tc[1]
            -4.03838000e-01
            -4.07219100e+00 * tc[0]
            -1.06564800e-03 * tc[1]
            +8.84690833e-08 * tc[2]
            -5.09355750e-12 * tc[3]
            +1.42058250e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.80069600e+04 / tc[1]
            +3.07203000e+00
            -4.57316700e+00 * tc[0]
            -2.16806800e-03 * tc[1]
            +2.45781500e-07 * tc[2]
            -1.95742000e-11 * tc[3]
            +7.15827000e-16 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.22797700e+02 / tc[1]
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988500e-04 * tc[1]
            +9.47460167e-08 * tc[2]
            -8.41420000e-12 * tc[3]
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
            +2.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 1: H */
        species[1] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 2: O */
        species[2] =
            +1.94642900e+00
            -1.63816600e-03 * tc[1]
            +2.42103200e-06 * tc[2]
            -1.60284300e-09 * tc[3]
            +3.89069600e-13 * tc[4];
        /*species 3: O2 */
        species[3] =
            +2.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
        /*species 4: OH */
        species[4] =
            +2.63726600e+00
            +1.85091000e-04 * tc[1]
            -1.67616500e-06 * tc[2]
            +2.38720300e-09 * tc[3]
            -8.43144200e-13 * tc[4];
        /*species 5: H2O */
        species[5] =
            +2.38684200e+00
            +3.47498200e-03 * tc[1]
            -6.35469600e-06 * tc[2]
            +6.96858100e-09 * tc[3]
            -2.50658800e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.97996300e+00
            +4.99669700e-03 * tc[1]
            -3.79099700e-06 * tc[2]
            +2.35419200e-09 * tc[3]
            -8.08902400e-13 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +2.38875400e+00
            +6.56922600e-03 * tc[1]
            -1.48501300e-07 * tc[2]
            -4.62580600e-09 * tc[3]
            +2.47151500e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +1.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 1: H */
        species[1] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 2: O */
        species[2] =
            +1.54206000e+00
            -2.75506200e-05 * tc[1]
            -3.10280300e-09 * tc[2]
            +4.55106700e-12 * tc[3]
            -4.36805200e-16 * tc[4];
        /*species 3: O2 */
        species[3] =
            +2.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
        /*species 4: OH */
        species[4] =
            +1.88273000e+00
            +1.01397400e-03 * tc[1]
            -2.27687700e-07 * tc[2]
            +2.17468400e-11 * tc[3]
            -5.12630500e-16 * tc[4];
        /*species 5: H2O */
        species[5] =
            +1.67214600e+00
            +3.05629300e-03 * tc[1]
            -8.73026000e-07 * tc[2]
            +1.20099600e-10 * tc[3]
            -6.39161800e-15 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.07219100e+00
            +2.13129600e-03 * tc[1]
            -5.30814500e-07 * tc[2]
            +6.11226900e-11 * tc[3]
            -2.84116500e-15 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.57316700e+00
            +4.33613600e-03 * tc[1]
            -1.47468900e-06 * tc[2]
            +2.34890400e-10 * tc[3]
            -1.43165400e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
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
            +3.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 1: H */
        species[1] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 2: O */
        species[2] =
            +2.94642900e+00
            -1.63816600e-03 * tc[1]
            +2.42103200e-06 * tc[2]
            -1.60284300e-09 * tc[3]
            +3.89069600e-13 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.63726600e+00
            +1.85091000e-04 * tc[1]
            -1.67616500e-06 * tc[2]
            +2.38720300e-09 * tc[3]
            -8.43144200e-13 * tc[4];
        /*species 5: H2O */
        species[5] =
            +3.38684200e+00
            +3.47498200e-03 * tc[1]
            -6.35469600e-06 * tc[2]
            +6.96858100e-09 * tc[3]
            -2.50658800e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +2.97996300e+00
            +4.99669700e-03 * tc[1]
            -3.79099700e-06 * tc[2]
            +2.35419200e-09 * tc[3]
            -8.08902400e-13 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.38875400e+00
            +6.56922600e-03 * tc[1]
            -1.48501300e-07 * tc[2]
            -4.62580600e-09 * tc[3]
            +2.47151500e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 1: H */
        species[1] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 2: O */
        species[2] =
            +2.54206000e+00
            -2.75506200e-05 * tc[1]
            -3.10280300e-09 * tc[2]
            +4.55106700e-12 * tc[3]
            -4.36805200e-16 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
        /*species 4: OH */
        species[4] =
            +2.88273000e+00
            +1.01397400e-03 * tc[1]
            -2.27687700e-07 * tc[2]
            +2.17468400e-11 * tc[3]
            -5.12630500e-16 * tc[4];
        /*species 5: H2O */
        species[5] =
            +2.67214600e+00
            +3.05629300e-03 * tc[1]
            -8.73026000e-07 * tc[2]
            +1.20099600e-10 * tc[3]
            -6.39161800e-15 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.07219100e+00
            +2.13129600e-03 * tc[1]
            -5.30814500e-07 * tc[2]
            +6.11226900e-11 * tc[3]
            -2.84116500e-15 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.57316700e+00
            +4.33613600e-03 * tc[1]
            -1.47468900e-06 * tc[2]
            +2.34890400e-10 * tc[3]
            -1.43165400e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
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
            +2.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 / tc[1];
        /*species 1: H */
        species[1] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 / tc[1];
        /*species 2: O */
        species[2] =
            +1.94642900e+00
            -8.19083000e-04 * tc[1]
            +8.07010667e-07 * tc[2]
            -4.00710750e-10 * tc[3]
            +7.78139200e-14 * tc[4]
            +2.91476400e+04 / tc[1];
        /*species 3: O2 */
        species[3] =
            +2.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 / tc[1];
        /*species 4: OH */
        species[4] =
            +2.63726600e+00
            +9.25455000e-05 * tc[1]
            -5.58721667e-07 * tc[2]
            +5.96800750e-10 * tc[3]
            -1.68628840e-13 * tc[4]
            +3.60678200e+03 / tc[1];
        /*species 5: H2O */
        species[5] =
            +2.38684200e+00
            +1.73749100e-03 * tc[1]
            -2.11823200e-06 * tc[2]
            +1.74214525e-09 * tc[3]
            -5.01317600e-13 * tc[4]
            -3.02081100e+04 / tc[1];
        /*species 6: HO2 */
        species[6] =
            +1.97996300e+00
            +2.49834850e-03 * tc[1]
            -1.26366567e-06 * tc[2]
            +5.88548000e-10 * tc[3]
            -1.61780480e-13 * tc[4]
            +1.76227400e+02 / tc[1];
        /*species 7: H2O2 */
        species[7] =
            +2.38875400e+00
            +3.28461300e-03 * tc[1]
            -4.95004333e-08 * tc[2]
            -1.15645150e-09 * tc[3]
            +4.94303000e-13 * tc[4]
            -1.76631500e+04 / tc[1];
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 / tc[1];
    } else {
        /*species 0: H2 */
        species[0] =
            +1.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 / tc[1];
        /*species 1: H */
        species[1] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 / tc[1];
        /*species 2: O */
        species[2] =
            +1.54206000e+00
            -1.37753100e-05 * tc[1]
            -1.03426767e-09 * tc[2]
            +1.13776675e-12 * tc[3]
            -8.73610400e-17 * tc[4]
            +2.92308000e+04 / tc[1];
        /*species 3: O2 */
        species[3] =
            +2.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 / tc[1];
        /*species 4: OH */
        species[4] =
            +1.88273000e+00
            +5.06987000e-04 * tc[1]
            -7.58959000e-08 * tc[2]
            +5.43671000e-12 * tc[3]
            -1.02526100e-16 * tc[4]
            +3.88688800e+03 / tc[1];
        /*species 5: H2O */
        species[5] =
            +1.67214600e+00
            +1.52814650e-03 * tc[1]
            -2.91008667e-07 * tc[2]
            +3.00249000e-11 * tc[3]
            -1.27832360e-15 * tc[4]
            -2.98992100e+04 / tc[1];
        /*species 6: HO2 */
        species[6] =
            +3.07219100e+00
            +1.06564800e-03 * tc[1]
            -1.76938167e-07 * tc[2]
            +1.52806725e-11 * tc[3]
            -5.68233000e-16 * tc[4]
            -1.57972700e+02 / tc[1];
        /*species 7: H2O2 */
        species[7] =
            +3.57316700e+00
            +2.16806800e-03 * tc[1]
            -4.91563000e-07 * tc[2]
            +5.87226000e-11 * tc[3]
            -2.86330800e-15 * tc[4]
            -1.80069600e+04 / tc[1];
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
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
            +3.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 / tc[1];
        /*species 1: H */
        species[1] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 / tc[1];
        /*species 2: O */
        species[2] =
            +2.94642900e+00
            -8.19083000e-04 * tc[1]
            +8.07010667e-07 * tc[2]
            -4.00710750e-10 * tc[3]
            +7.78139200e-14 * tc[4]
            +2.91476400e+04 / tc[1];
        /*species 3: O2 */
        species[3] =
            +3.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 / tc[1];
        /*species 4: OH */
        species[4] =
            +3.63726600e+00
            +9.25455000e-05 * tc[1]
            -5.58721667e-07 * tc[2]
            +5.96800750e-10 * tc[3]
            -1.68628840e-13 * tc[4]
            +3.60678200e+03 / tc[1];
        /*species 5: H2O */
        species[5] =
            +3.38684200e+00
            +1.73749100e-03 * tc[1]
            -2.11823200e-06 * tc[2]
            +1.74214525e-09 * tc[3]
            -5.01317600e-13 * tc[4]
            -3.02081100e+04 / tc[1];
        /*species 6: HO2 */
        species[6] =
            +2.97996300e+00
            +2.49834850e-03 * tc[1]
            -1.26366567e-06 * tc[2]
            +5.88548000e-10 * tc[3]
            -1.61780480e-13 * tc[4]
            +1.76227400e+02 / tc[1];
        /*species 7: H2O2 */
        species[7] =
            +3.38875400e+00
            +3.28461300e-03 * tc[1]
            -4.95004333e-08 * tc[2]
            -1.15645150e-09 * tc[3]
            +4.94303000e-13 * tc[4]
            -1.76631500e+04 / tc[1];
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 / tc[1];
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 / tc[1];
        /*species 1: H */
        species[1] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 / tc[1];
        /*species 2: O */
        species[2] =
            +2.54206000e+00
            -1.37753100e-05 * tc[1]
            -1.03426767e-09 * tc[2]
            +1.13776675e-12 * tc[3]
            -8.73610400e-17 * tc[4]
            +2.92308000e+04 / tc[1];
        /*species 3: O2 */
        species[3] =
            +3.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 / tc[1];
        /*species 4: OH */
        species[4] =
            +2.88273000e+00
            +5.06987000e-04 * tc[1]
            -7.58959000e-08 * tc[2]
            +5.43671000e-12 * tc[3]
            -1.02526100e-16 * tc[4]
            +3.88688800e+03 / tc[1];
        /*species 5: H2O */
        species[5] =
            +2.67214600e+00
            +1.52814650e-03 * tc[1]
            -2.91008667e-07 * tc[2]
            +3.00249000e-11 * tc[3]
            -1.27832360e-15 * tc[4]
            -2.98992100e+04 / tc[1];
        /*species 6: HO2 */
        species[6] =
            +4.07219100e+00
            +1.06564800e-03 * tc[1]
            -1.76938167e-07 * tc[2]
            +1.52806725e-11 * tc[3]
            -5.68233000e-16 * tc[4]
            -1.57972700e+02 / tc[1];
        /*species 7: H2O2 */
        species[7] =
            +4.57316700e+00
            +2.16806800e-03 * tc[1]
            -4.91563000e-07 * tc[2]
            +5.87226000e-11 * tc[3]
            -2.86330800e-15 * tc[4]
            -1.80069600e+04 / tc[1];
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
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
            +3.29812400e+00 * tc[0]
            +8.24944200e-04 * tc[1]
            -4.07150750e-07 * tc[2]
            -3.15847800e-11 * tc[3]
            +1.03371800e-13 * tc[4]
            -3.29409400e+00 ;
        /*species 1: H */
        species[1] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117600e-01 ;
        /*species 2: O */
        species[2] =
            +2.94642900e+00 * tc[0]
            -1.63816600e-03 * tc[1]
            +1.21051600e-06 * tc[2]
            -5.34281000e-10 * tc[3]
            +9.72674000e-14 * tc[4]
            +2.96399500e+00 ;
        /*species 3: O2 */
        species[3] =
            +3.21293600e+00 * tc[0]
            +1.12748600e-03 * tc[1]
            -2.87807500e-07 * tc[2]
            +4.37959000e-10 * tc[3]
            -2.19213850e-13 * tc[4]
            +6.03473800e+00 ;
        /*species 4: OH */
        species[4] =
            +3.63726600e+00 * tc[0]
            +1.85091000e-04 * tc[1]
            -8.38082500e-07 * tc[2]
            +7.95734333e-10 * tc[3]
            -2.10786050e-13 * tc[4]
            +1.35886000e+00 ;
        /*species 5: H2O */
        species[5] =
            +3.38684200e+00 * tc[0]
            +3.47498200e-03 * tc[1]
            -3.17734800e-06 * tc[2]
            +2.32286033e-09 * tc[3]
            -6.26647000e-13 * tc[4]
            +2.59023300e+00 ;
        /*species 6: HO2 */
        species[6] =
            +2.97996300e+00 * tc[0]
            +4.99669700e-03 * tc[1]
            -1.89549850e-06 * tc[2]
            +7.84730667e-10 * tc[3]
            -2.02225600e-13 * tc[4]
            +9.22272400e+00 ;
        /*species 7: H2O2 */
        species[7] =
            +3.38875400e+00 * tc[0]
            +6.56922600e-03 * tc[1]
            -7.42506500e-08 * tc[2]
            -1.54193533e-09 * tc[3]
            +6.17878750e-13 * tc[4]
            +6.78536300e+00 ;
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00 * tc[0]
            +1.40824000e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213750e-13 * tc[4]
            +3.95037200e+00 ;
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142300e+00 * tc[0]
            +7.00064400e-04 * tc[1]
            -2.81691450e-08 * tc[2]
            -3.07719267e-12 * tc[3]
            +3.95688000e-16 * tc[4]
            -1.35511000e+00 ;
        /*species 1: H */
        species[1] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117600e-01 ;
        /*species 2: O */
        species[2] =
            +2.54206000e+00 * tc[0]
            -2.75506200e-05 * tc[1]
            -1.55140150e-09 * tc[2]
            +1.51702233e-12 * tc[3]
            -1.09201300e-16 * tc[4]
            +4.92030800e+00 ;
        /*species 3: O2 */
        species[3] =
            +3.69757800e+00 * tc[0]
            +6.13519700e-04 * tc[1]
            -6.29421000e-08 * tc[2]
            +5.91760333e-12 * tc[3]
            -2.84108750e-16 * tc[4]
            +3.18916600e+00 ;
        /*species 4: OH */
        species[4] =
            +2.88273000e+00 * tc[0]
            +1.01397400e-03 * tc[1]
            -1.13843850e-07 * tc[2]
            +7.24894667e-12 * tc[3]
            -1.28157625e-16 * tc[4]
            +5.59571200e+00 ;
        /*species 5: H2O */
        species[5] =
            +2.67214600e+00 * tc[0]
            +3.05629300e-03 * tc[1]
            -4.36513000e-07 * tc[2]
            +4.00332000e-11 * tc[3]
            -1.59790450e-15 * tc[4]
            +6.86281700e+00 ;
        /*species 6: HO2 */
        species[6] =
            +4.07219100e+00 * tc[0]
            +2.13129600e-03 * tc[1]
            -2.65407250e-07 * tc[2]
            +2.03742300e-11 * tc[3]
            -7.10291250e-16 * tc[4]
            +3.47602900e+00 ;
        /*species 7: H2O2 */
        species[7] =
            +4.57316700e+00 * tc[0]
            +4.33613600e-03 * tc[1]
            -7.37344500e-07 * tc[2]
            +7.82968000e-11 * tc[3]
            -3.57913500e-15 * tc[4]
            +5.01137000e-01 ;
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00 * tc[0]
            +1.48797700e-03 * tc[1]
            -2.84238050e-07 * tc[2]
            +3.36568000e-11 * tc[3]
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
    wt[8] = 28.013400; /*N2 */

    return;
}


/*get temperature given internal energy in mass units and mass fracs */
int feeytt_(double * e, double * y, int * iwrk, double * rwrk, double * t)
{
    const int maxiter = 50;
    const double tol  = 0.001;
    double ein  = *e;
    double tmin = 300; // max lower bound for thermo def
    double tmax = 5000; // min upper bound for thermo def
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
    y[8] = phi[8]*28.013400;   XW += y[8]; /*N2 */
    for (id = 0; id < 9; ++id) {
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
    phi[8] = y[8]/ 2.80134000e-02; /*N2 (wt in kg) */

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
    y[8] = c[8] * 28.013400 / (*rho); 

    return;
}


/*ddebdf compatible right hand side of CV burner */
/*rwrk[0] and rwrk[1] should contain rho and ene respectively */
/*working variable phi contains specific mole numbers */
void fecvrhs_(double * time, double * phi, double * phidot, double * rwrk, int * iwrk)
{
    double rho,ene; /*CV Parameters */
    double y[9], wdot[9]; /*temporary storage */
    int i; /*Loop counter */
    double temperature,pressure; /*temporary var */
    rho = rwrk[0];
    ene = rwrk[1];
    fephity_(phi, iwrk, rwrk, y);
    feeytt_(&ene, y, iwrk, rwrk, &temperature);
    CKPY(&rho, &temperature,  y, iwrk, rwrk, &pressure);
    CKWYP(&pressure, &temperature,  y, iwrk, rwrk, wdot);
    for (i=0; i<9; ++i) phidot[i] = wdot[i] / (rho/1000.0); 

    return;
}


/*returns the dimensionality of the cv burner (number of species) */
int fecvdim_()
{
    return 9;
}


/*ddebdf compatible right hand side of ZND solver */
/*rwrk[0] : scaling factor for pressure */
/*rwrk[1] : preshock density (g/cc)  */
/*rwrk[2] : detonation velocity (cm/s)  */
/*solution vector: [P; rho; y0 ... ylast]  */
void fezndrhs_(double * time, double * z, double * zdot, double * rwrk, int * iwrk)
{
    double psc,rho1,udet; /*ZND Parameters */
    double wt[9], hms[9], wdot[9]; /*temporary storage */
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
    for (i=0; i<9; ++i) {
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
    return 12;
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
    if (strcmp(s1, "N2")==0) return 8; 
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
    if (sn==8) return "N2"; 
    /*species name not found */
    return "NOTFOUND";
}

/* End of file  */
