From - Thu Feb 12 08:37:22 2004
X-Mozilla-Status: 0001
X-Mozilla-Status2: 00000000
Return-path: <MSDay@lbl.gov>
Received: from mta1.lbl.gov (mta1.lbl.gov [128.3.41.24])
 by imapc.lbl.gov (iPlanet Messaging Server 5.2 HotFix 1.21 (built Sep  8
 2003)) with ESMTP id <0HSY00EI848XMP@imapc.lbl.gov> for mjlijewski@lbl.gov;
 Wed, 11 Feb 2004 16:34:09 -0800 (PST)
Received: from mta1.lbl.gov (localhost [127.0.0.1])
	by mta1.lbl.gov (8.12.10/8.12.10) with ESMTP id i1C0Y6qc009660	for
 <mjlijewski@imapc.lbl.gov>; Wed, 11 Feb 2004 16:34:07 -0800 (PST)
Received: from lbl.gov (hedorah.lbl.gov [128.3.5.24])
	by mta1.lbl.gov (8.12.10/8.12.10) with ESMTP id i1C0Y6OW009656	for
 <MJLijewski@lbl.gov>; Wed, 11 Feb 2004 16:34:06 -0800 (PST)
Date: Wed, 11 Feb 2004 16:35:04 -0800
From: Marc Day <MSDay@lbl.gov>
Subject: Fuego file
To: Mike Lijewski <MJLijewski@lbl.gov>
Message-id: <402ACA38.3090007@lbl.gov>
Organization: Lawrence Berkeley National Laboratory
MIME-version: 1.0
Content-type: multipart/mixed; boundary=------------030905060905030504010201
X-Accept-Language: en,pdf
User-Agent: Mozilla/5.0 (Windows; U; Windows NT 5.0; en-US; rv:1.4)
 Gecko/20030624
Original-recipient: rfc822;mjlijewski@imapc.lbl.gov

This is a multi-part message in MIME format.
--------------030905060905030504010201
Content-Type: text/plain; charset=us-ascii; format=flowed
Content-Transfer-Encoding: 7bit

Mike,
Attached you'll find the fuego output for the chem-CH4-2step.inp
mech (but with AR removed).  Can you add this to the LMC capability?
-M


--------------030905060905030504010201
Content-Type: text/plain;
 name="chem-CH4-2step.c"
Content-Transfer-Encoding: 7bit
Content-Disposition: inline;
 filename="chem-CH4-2step.c"

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
void fgindx_(int * iwrk, double *rwrk, int * mm, int * kk, int * ii, int * nfit );
void fgxnum_(char * line, int * nexp, int * lout, int * nval, double * rval, int * kerr, int lenline);
void fgsnum_(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * rval, int * kerr, int lenline, int lenkray);
void fgsyms_(char * cckwrk, int * lout, char * kname, int * kerr, int lencck, int lenkname);
void fgrp_(int * ickwrk, double * rckwrk, double * ru, double * ruc, double * pa);
void fgpx_(double * rho, double * T, double * x, int * iwrk, double *rwrk, double * P);
void fgpy_(double * rho, double * T, double * y, int * iwrk, double *rwrk, double * P);
void fgpc_(double * rho, double * T, double * c, int * iwrk, double *rwrk, double * P);
void fgrhox_(double * P, double * T, double * x, int * iwrk, double *rwrk, double * rho);
void fgrhoy_(double * P, double * T, double * y, int * iwrk, double *rwrk, double * rho);
void fgrhoc_(double * P, double * T, double * c, int * iwrk, double *rwrk, double * rho);
void fgwt_(int * iwrk, double *rwrk, double * wt);
void fgmmwy_(double * y, int * iwrk, double * rwrk, double * wtm);
void fgmmwx_(double * x, int * iwrk, double * rwrk, double * wtm);
void fgmmwc_(double * c, int * iwrk, double * rwrk, double * wtm);
void fgytx_(double * y, int * iwrk, double * rwrk, double * x);
void fgytcp_(double * P, double * T, double * y, int * iwrk, double * rwrk, double * c);
void fgytcr_(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * c);
void fgxty_(double * x, int * iwrk, double * rwrk, double * y);
void fgxtcp_(double * P, double * T, double * x, int * iwrk, double * rwrk, double * c);
void fgxtcr_(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * c);
void fgctx_(double * c, int * iwrk, double * rwrk, double * x);
void fgcty_(double * c, int * iwrk, double * rwrk, double * y);
void fgcpor_(double * T, int * iwrk, double * rwrk, double * cpor);
void fghort_(double * T, int * iwrk, double * rwrk, double * hort);
void fgsor_(double * T, int * iwrk, double * rwrk, double * sor);
void fgcvml_(double * T, int * iwrk, double * rwrk, double * cvml);
void fgcpml_(double * T, int * iwrk, double * rwrk, double * cvml);
void fguml_(double * T, int * iwrk, double * rwrk, double * uml);
void fghml_(double * T, int * iwrk, double * rwrk, double * uml);
void fggml_(double * T, int * iwrk, double * rwrk, double * gml);
void fgaml_(double * T, int * iwrk, double * rwrk, double * aml);
void fgsml_(double * T, int * iwrk, double * rwrk, double * sml);
void fgcvms_(double * T, int * iwrk, double * rwrk, double * cvms);
void fgcpms_(double * T, int * iwrk, double * rwrk, double * cvms);
void fgums_(double * T, int * iwrk, double * rwrk, double * ums);
void fghms_(double * T, int * iwrk, double * rwrk, double * ums);
void fggms_(double * T, int * iwrk, double * rwrk, double * gms);
void fgams_(double * T, int * iwrk, double * rwrk, double * ams);
void fgsms_(double * T, int * iwrk, double * rwrk, double * sms);
void fgcpbl_(double * T, double * x, int * iwrk, double * rwrk, double * cpbl);
void fgcpbs_(double * T, double * y, int * iwrk, double * rwrk, double * cpbs);
void fgcvbl_(double * T, double * x, int * iwrk, double * rwrk, double * cpbl);
void fgcvbs_(double * T, double * y, int * iwrk, double * rwrk, double * cpbs);
void fghbml_(double * T, double * x, int * iwrk, double * rwrk, double * hbml);
void fghbms_(double * T, double * y, int * iwrk, double * rwrk, double * hbms);
void fgubml_(double * T, double * x, int * iwrk, double * rwrk, double * ubml);
void fgubms_(double * T, double * y, int * iwrk, double * rwrk, double * ubms);
void fgsbml_(double * P, double * T, double * x, int * iwrk, double * rwrk, double * sbml);
void fgsbms_(double * P, double * T, double * y, int * iwrk, double * rwrk, double * sbms);
void fggbml_(double * P, double * T, double * x, int * iwrk, double * rwrk, double * gbml);
void fggbms_(double * P, double * T, double * y, int * iwrk, double * rwrk, double * gbms);
void fgabml_(double * P, double * T, double * x, int * iwrk, double * rwrk, double * abml);
void fgabms_(double * P, double * T, double * y, int * iwrk, double * rwrk, double * abms);
void fgwc_(double * T, double * C, int * iwrk, double *rwrk, double * wdot);
void fgwyp_(double * P, double * T, double * y, int * iwrk, double *rwrk, double * wdot);
void fgwxp_(double * P, double * T, double * x, int * iwrk, double *rwrk, double * wdot);
void fgwyr_(double * rho, double * T, double * y, int * iwrk, double *rwrk, double * wdot);
void fgwxr_(double * rho, double * T, double * x, int * iwrk, double *rwrk, double * wdot);
void fgqc_(double * T, double * C, int * iwrk, double *rwrk, double * qdot);
void fgqyp_(double * P, double * T, double * y, int * iwrk, double *rwrk, double * qdot);
void fgqxp_(double * P, double * T, double * x, int * iwrk, double *rwrk, double * qdot);
void fgqyr_(double * rho, double * T, double * y, int * iwrk, double *rwrk, double * qdot);
void fgqxr_(double * rho, double * T, double * x, int * iwrk, double *rwrk, double * qdot);
void fgnu_(int * kdim, int * iwrk, double *rwrk, int * nuki);
void fgncf_(int * mdim, int * iwrk, double *rwrk, int * ncf);
void fgabe_(int * iwrk, double *rwrk, double * a, double * b, double * e );
void fgeqc_(double * T, double * C , int * iwrk, double *rwrk, double * eqcon );
void fgeqyp_(double * P, double * T, double * y, int * iwrk, double *rwrk, double * eqcon);
void fgeqxp_(double * P, double * T, double * x, int * iwrk, double *rwrk, double * eqcon);
void fgeqyr_(double * rho, double * T, double * y, int * iwrk, double *rwrk, double * eqcon);
void fgeqxr_(double * rho, double * T, double * x, int * iwrk, double *rwrk, double * eqcon);
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
void fgindx_(int * iwrk, double * rwrk, int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 4;
    *kk = 6;
    *ii = 3;
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
void fgxnum_(char * line, int * nexp, int * lout, int * nval, double * rval, int * kerr, int lenline )
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
void fgsnum_(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * rval, int * kerr, int lenline, int lenkray)
{
    /*Not done yet ... */
}


/* Returns the char strings of species names */
void fgsyms_(char * cckwrk, int * lout, char * kname, int * kerr, int lencck, int lenkname )
{
    int i; /*Loop Counter */
    /*clear kname */
    for (i=0; i<lenkname*6; i++) {
        kname[i] = ' ';
    }

    /* O2  */
    kname[ 0*lenkname + 0 ] = 'O';
    kname[ 0*lenkname + 1 ] = '2';
    kname[ 0*lenkname + 2 ] = ' ';

    /* H2O  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = '2';
    kname[ 1*lenkname + 2 ] = 'O';
    kname[ 1*lenkname + 3 ] = ' ';

    /* CH4  */
    kname[ 2*lenkname + 0 ] = 'C';
    kname[ 2*lenkname + 1 ] = 'H';
    kname[ 2*lenkname + 2 ] = '4';
    kname[ 2*lenkname + 3 ] = ' ';

    /* CO  */
    kname[ 3*lenkname + 0 ] = 'C';
    kname[ 3*lenkname + 1 ] = 'O';
    kname[ 3*lenkname + 2 ] = ' ';

    /* CO2  */
    kname[ 4*lenkname + 0 ] = 'C';
    kname[ 4*lenkname + 1 ] = 'O';
    kname[ 4*lenkname + 2 ] = '2';
    kname[ 4*lenkname + 3 ] = ' ';

    /* N2  */
    kname[ 5*lenkname + 0 ] = 'N';
    kname[ 5*lenkname + 1 ] = '2';
    kname[ 5*lenkname + 2 ] = ' ';

}


/* Returns R, Rc, Patm */
void fgrp_(int * ickwrk, double * rckwrk, double * ru, double * ruc, double * pa)
{
     *ru  = 8.314e+07; 
     *ruc = 1.987; 
     *pa  = 1.01325e+06; 
}


/*Compute P = rhoRT/W(x) */
void fgpx_(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * P)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*31.998800; /*O2 */
    XW += x[1]*18.015340; /*H2O */
    XW += x[2]*16.042880; /*CH4 */
    XW += x[3]*28.010400; /*CO */
    XW += x[4]*44.009800; /*CO2 */
    XW += x[5]*28.013400; /*N2 */
    *P = *rho * 8.314e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void fgpy_(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]/31.998800; /*O2 */
    YOW += y[1]/18.015340; /*H2O */
    YOW += y[2]/16.042880; /*CH4 */
    YOW += y[3]/28.010400; /*CO */
    YOW += y[4]/44.009800; /*CO2 */
    YOW += y[5]/28.013400; /*N2 */
    *P = *rho * 8.314e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(c) */
void fgpc_(double * rho, double * T, double * c, int * iwrk, double * rwrk, double * P)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*31.998800; /*O2 */
    W += c[1]*18.015340; /*H2O */
    W += c[2]*16.042880; /*CH4 */
    W += c[3]*28.010400; /*CO */
    W += c[4]*44.009800; /*CO2 */
    W += c[5]*28.013400; /*N2 */

    for (id = 0; id < 6; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.314e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void fgrhox_(double * P, double * T, double * x, int * iwrk, double * rwrk, double * rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*31.998800; /*O2 */
    XW += x[1]*18.015340; /*H2O */
    XW += x[2]*16.042880; /*CH4 */
    XW += x[3]*28.010400; /*CO */
    XW += x[4]*44.009800; /*CO2 */
    XW += x[5]*28.013400; /*N2 */
    *rho = *P * XW / (8.314e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void fgrhoy_(double * P, double * T, double * y, int * iwrk, double * rwrk, double * rho)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]/31.998800; /*O2 */
    YOW += y[1]/18.015340; /*H2O */
    YOW += y[2]/16.042880; /*CH4 */
    YOW += y[3]/28.010400; /*CO */
    YOW += y[4]/44.009800; /*CO2 */
    YOW += y[5]/28.013400; /*N2 */
    *rho = *P / (8.314e+07 * (*T) * YOW); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(c)/(R*T) */
void fgrhoc_(double * P, double * T, double * c, int * iwrk, double * rwrk, double * rho)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*31.998800; /*O2 */
    W += c[1]*18.015340; /*H2O */
    W += c[2]*16.042880; /*CH4 */
    W += c[3]*28.010400; /*CO */
    W += c[4]*44.009800; /*CO2 */
    W += c[5]*28.013400; /*N2 */

    for (id = 0; id < 6; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.314e+07); /*rho = PW/(R*T) */

    return;
}


/*get molecular weight for all species */
void fgwt_(int * iwrk, double * rwrk, double * wt)
{
    molecularWeight(wt);
}


/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
void fgmmwy_(double *y, int * iwrk, double * rwrk, double * wtm)
{
    double YOW = 0;/* see Eq 3 in CK Manual */
    YOW += y[0]/31.998800; /*O2 */
    YOW += y[1]/18.015340; /*H2O */
    YOW += y[2]/16.042880; /*CH4 */
    YOW += y[3]/28.010400; /*CO */
    YOW += y[4]/44.009800; /*CO2 */
    YOW += y[5]/28.013400; /*N2 */
    *wtm = 1.0 / YOW;

    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
void fgmmwx_(double *x, int * iwrk, double * rwrk, double * wtm)
{
    double XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*31.998800; /*O2 */
    XW += x[1]*18.015340; /*H2O */
    XW += x[2]*16.042880; /*CH4 */
    XW += x[3]*28.010400; /*CO */
    XW += x[4]*44.009800; /*CO2 */
    XW += x[5]*28.013400; /*N2 */
    *wtm = XW;

    return;
}


/*given c[species]: molar concentration */
/*returns mean molecular weight (gm/mole) */
void fgmmwc_(double *c, int * iwrk, double * rwrk, double * wtm)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*31.998800; /*O2 */
    W += c[1]*18.015340; /*H2O */
    W += c[2]*16.042880; /*CH4 */
    W += c[3]*28.010400; /*CO */
    W += c[4]*44.009800; /*CO2 */
    W += c[5]*28.013400; /*N2 */

    for (id = 0; id < 6; ++id) {
        sumC += c[id];
    }
    /* CK provides no guard against divison by zero */
    *wtm = W/sumC;

    return;
}


/*convert y[species] (mass fracs) to x[species] (mole fracs) */
void fgytx_(double * y, int * iwrk, double * rwrk, double * x)
{
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/31.998800; /*O2 */
    YOW += y[1]/18.015340; /*H2O */
    YOW += y[2]/16.042880; /*CH4 */
    YOW += y[3]/28.010400; /*CO */
    YOW += y[4]/44.009800; /*CO2 */
    YOW += y[5]/28.013400; /*N2 */
    /*Now compute conversion */
    x[0] = y[0]/(31.998800*YOW); 
    x[1] = y[1]/(18.015340*YOW); 
    x[2] = y[2]/(16.042880*YOW); 
    x[3] = y[3]/(28.010400*YOW); 
    x[4] = y[4]/(44.009800*YOW); 
    x[5] = y[5]/(28.013400*YOW); 

    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void fgytcp_(double * P, double * T, double * y, int * iwrk, double * rwrk, double * c)
{
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/31.998800; /*O2 */
    YOW += y[1]/18.015340; /*H2O */
    YOW += y[2]/16.042880; /*CH4 */
    YOW += y[3]/28.010400; /*CO */
    YOW += y[4]/44.009800; /*CO2 */
    YOW += y[5]/28.013400; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.314e+07 * (*T)); 
    /*Now compute conversion */
    c[0] = PWORT * y[0]/31.998800; 
    c[1] = PWORT * y[1]/18.015340; 
    c[2] = PWORT * y[2]/16.042880; 
    c[3] = PWORT * y[3]/28.010400; 
    c[4] = PWORT * y[4]/44.009800; 
    c[5] = PWORT * y[5]/28.013400; 

    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void fgytcr_(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * c)
{
    /*See Eq 8 (Temperature not used) */
    c[0] = (*rho) * y[0]/31.998800; 
    c[1] = (*rho) * y[1]/18.015340; 
    c[2] = (*rho) * y[2]/16.042880; 
    c[3] = (*rho) * y[3]/28.010400; 
    c[4] = (*rho) * y[4]/44.009800; 
    c[5] = (*rho) * y[5]/28.013400; 

    return;
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
void fgxty_(double * x, int * iwrk, double * rwrk, double * y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*31.998800; /*O2 */
    XW += x[1]*18.015340; /*H2O */
    XW += x[2]*16.042880; /*CH4 */
    XW += x[3]*28.010400; /*CO */
    XW += x[4]*44.009800; /*CO2 */
    XW += x[5]*28.013400; /*N2 */
    /*Now compute conversion */
    y[0] = x[0]*31.998800/XW; 
    y[1] = x[1]*18.015340/XW; 
    y[2] = x[2]*16.042880/XW; 
    y[3] = x[3]*28.010400/XW; 
    y[4] = x[4]*44.009800/XW; 
    y[5] = x[5]*28.013400/XW; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void fgxtcp_(double * P, double * T, double * x, int * iwrk, double * rwrk, double * c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.314e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 6; ++id) {
        c[id] = x[id]*PORT;
    }

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void fgxtcr_(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * c)
{
    int id; /*loop counter */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*31.998800; /*O2 */
    XW += x[1]*18.015340; /*H2O */
    XW += x[2]*16.042880; /*CH4 */
    XW += x[3]*28.010400; /*CO */
    XW += x[4]*44.009800; /*CO2 */
    XW += x[5]*28.013400; /*N2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 6; ++id) {
        c[id] = x[id]*ROW;
    }

    return;
}


/*convert c[species] (molar conc) to x[species] (mole fracs) */
void fgctx_(double * c, int * iwrk, double * rwrk, double * x)
{
    int id; /*loop counter */
    double sumC = 0; 

    /*compute sum of c  */
    for (id = 0; id < 6; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    for (id = 0; id < 6; ++id) {
        x[id] = c[id]/sumC;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
void fgcty_(double * c, int * iwrk, double * rwrk, double * y)
{
    double CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*31.998800; /*O2 */
    CW += c[1]*18.015340; /*H2O */
    CW += c[2]*16.042880; /*CH4 */
    CW += c[3]*28.010400; /*CO */
    CW += c[4]*44.009800; /*CO2 */
    CW += c[5]*28.013400; /*N2 */
    /*Now compute conversion */
    y[0] = c[0]*31.998800/CW; 
    y[1] = c[1]*18.015340/CW; 
    y[2] = c[2]*16.042880/CW; 
    y[3] = c[3]*28.010400/CW; 
    y[4] = c[4]*44.009800/CW; 
    y[5] = c[5]*28.013400/CW; 

    return;
}


/*get Cp/R as a function of T  */
/*for all species (Eq 19) */
void fgcpor_(double *T, int * iwrk, double * rwrk, double * cpor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpor, tc);
}


/*get H/RT as a function of T  */
/*for all species (Eq 20) */
void fghort_(double *T, int * iwrk, double * rwrk, double * hort)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEnthalpy(hort, tc);
}


/*get S/R as a function of T  */
/*for all species (Eq 21) */
void fgsor_(double *T, int * iwrk, double * rwrk, double * sor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sor, tc);
}


/*get specific heat at constant volume as a function  */
/*of T for all species (molar units) */
void fgcvml_(double *T, int * iwrk, double * rwrk, double * cvml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 6; ++id) {
        cvml[id] *= 8.314e+07;
    }
}


/*get specific heat at constant pressure as a  */
/*function of T for all species (molar units) */
void fgcpml_(double *T, int * iwrk, double * rwrk, double * cpml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 6; ++id) {
        cpml[id] *= 8.314e+07;
    }
}


/*get internal energy as a function  */
/*of T for all species (molar units) */
void fguml_(double *T, int * iwrk, double * rwrk, double * uml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 6; ++id) {
        uml[id] *= RT;
    }
}


/*get enthalpy as a function  */
/*of T for all species (molar units) */
void fghml_(double *T, int * iwrk, double * rwrk, double * hml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 6; ++id) {
        hml[id] *= RT;
    }
}


/*get standard-state Gibbs energy as a function  */
/*of T for all species (molar units) */
void fggml_(double *T, int * iwrk, double * rwrk, double * gml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    gibbs(gml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 6; ++id) {
        gml[id] *= RT;
    }
}


/*get standard-state Helmholtz free energy as a  */
/*function of T for all species (molar units) */
void fgaml_(double *T, int * iwrk, double * rwrk, double * aml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    helmholtz(aml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 6; ++id) {
        aml[id] *= RT;
    }
}


/*Returns the standard-state entropies in molar units */
void fgsml_(double *T, int * iwrk, double * rwrk, double * sml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 6; ++id) {
        sml[id] *= 8.314e+07;
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
void fgcvms_(double *T, int * iwrk, double * rwrk, double * cvms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 2598222.433341; /*O2 */
    cvms[1] *= 4614955.920899; /*H2O */
    cvms[2] *= 5182361.271792; /*CH4 */
    cvms[3] *= 2968183.246223; /*CO */
    cvms[4] *= 1889124.694954; /*CO2 */
    cvms[5] *= 2967865.378712; /*N2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
void fgcpms_(double *T, int * iwrk, double * rwrk, double * cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 2598222.433341; /*O2 */
    cpms[1] *= 4614955.920899; /*H2O */
    cpms[2] *= 5182361.271792; /*CH4 */
    cpms[3] *= 2968183.246223; /*CO */
    cpms[4] *= 1889124.694954; /*CO2 */
    cpms[5] *= 2967865.378712; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
void fgums_(double *T, int * iwrk, double * rwrk, double * ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    ums[0] *= RT/31.998800; /*O2 */
    ums[1] *= RT/18.015340; /*H2O */
    ums[2] *= RT/16.042880; /*CH4 */
    ums[3] *= RT/28.010400; /*CO */
    ums[4] *= RT/44.009800; /*CO2 */
    ums[5] *= RT/28.013400; /*N2 */
}


/*Returns enthalpy in mass units (Eq 27.) */
void fghms_(double *T, int * iwrk, double * rwrk, double * hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    speciesEnthalpy(hms, tc);
    hms[0] *= RT/31.998800; /*O2 */
    hms[1] *= RT/18.015340; /*H2O */
    hms[2] *= RT/16.042880; /*CH4 */
    hms[3] *= RT/28.010400; /*CO */
    hms[4] *= RT/44.009800; /*CO2 */
    hms[5] *= RT/28.013400; /*N2 */
}


/*Returns gibbs in mass units (Eq 31.) */
void fggms_(double *T, int * iwrk, double * rwrk, double * gms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    gibbs(gms, tc);
    gms[0] *= RT/31.998800; /*O2 */
    gms[1] *= RT/18.015340; /*H2O */
    gms[2] *= RT/16.042880; /*CH4 */
    gms[3] *= RT/28.010400; /*CO */
    gms[4] *= RT/44.009800; /*CO2 */
    gms[5] *= RT/28.013400; /*N2 */
}


/*Returns helmholtz in mass units (Eq 32.) */
void fgams_(double *T, int * iwrk, double * rwrk, double * ams)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    helmholtz(ams, tc);
    ams[0] *= RT/31.998800; /*O2 */
    ams[1] *= RT/18.015340; /*H2O */
    ams[2] *= RT/16.042880; /*CH4 */
    ams[3] *= RT/28.010400; /*CO */
    ams[4] *= RT/44.009800; /*CO2 */
    ams[5] *= RT/28.013400; /*N2 */
}


/*Returns the entropies in mass units (Eq 28.) */
void fgsms_(double *T, int * iwrk, double * rwrk, double * sms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 2598222.433341; /*O2 */
    sms[1] *= 4614955.920899; /*H2O */
    sms[2] *= 5182361.271792; /*CH4 */
    sms[3] *= 2968183.246223; /*CO */
    sms[4] *= 1889124.694954; /*CO2 */
    sms[5] *= 2967865.378712; /*N2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void fgcpbl_(double *T, double *x, int * iwrk, double * rwrk, double * cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[6]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 6; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.314e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
void fgcpbs_(double *T, double *y, int * iwrk, double * rwrk, double * cpbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[6]; /* temporary storage */
    cp_R(cpor, tc);
    /*multiply by y/molecularweight */
    result += cpor[0]*y[0]/31.9988; /*O2 */
    result += cpor[1]*y[1]/18.0153; /*H2O */
    result += cpor[2]*y[2]/16.0429; /*CH4 */
    result += cpor[3]*y[3]/28.0104; /*CO */
    result += cpor[4]*y[4]/44.0098; /*CO2 */
    result += cpor[5]*y[5]/28.0134; /*N2 */

    *cpbs = result * 8.314e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void fgcvbl_(double *T, double *x, int * iwrk, double * rwrk, double * cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[6]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 6; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.314e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
void fgcvbs_(double *T, double *y, int * iwrk, double * rwrk, double * cvbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[6]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]/31.9988; /*O2 */
    result += cvor[1]*y[1]/18.0153; /*H2O */
    result += cvor[2]*y[2]/16.0429; /*CH4 */
    result += cvor[3]*y[3]/28.0104; /*CO */
    result += cvor[4]*y[4]/44.0098; /*CO2 */
    result += cvor[5]*y[5]/28.0134; /*N2 */

    *cvbs = result * 8.314e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void fghbml_(double *T, double *x, int * iwrk, double * rwrk, double * hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[6]; /* temporary storage */
    double RT = 8.314e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 6; ++id) {
        result += x[id]*hml[id];
    }

    *hbml = result * RT;
}


/*Returns mean enthalpy of mixture in mass units */
void fghbms_(double *T, double *y, int * iwrk, double * rwrk, double * hbms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[6]; /* temporary storage */
    double RT = 8.314e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*hml[0]/31.998800; /*O2 */
    result += y[1]*hml[1]/18.015340; /*H2O */
    result += y[2]*hml[2]/16.042880; /*CH4 */
    result += y[3]*hml[3]/28.010400; /*CO */
    result += y[4]*hml[4]/44.009800; /*CO2 */
    result += y[5]*hml[5]/28.013400; /*N2 */

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void fgubml_(double *T, double *x, int * iwrk, double * rwrk, double * ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[6]; /* temporary energy array */
    double RT = 8.314e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 6; ++id) {
        result += x[id]*uml[id];
    }

    *ubml = result * RT;
}


/*get mean internal energy in mass units */
void fgubms_(double *T, double *y, int * iwrk, double * rwrk, double * ubms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double ums[6]; /* temporary energy array */
    double RT = 8.314e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]/31.998800; /*O2 */
    result += y[1]*ums[1]/18.015340; /*H2O */
    result += y[2]*ums[2]/16.042880; /*CH4 */
    result += y[3]*ums[3]/28.010400; /*CO */
    result += y[4]*ums[4]/44.009800; /*CO2 */
    result += y[5]*ums[5]/28.013400; /*N2 */

    *ubms = result * RT;
}


/*get mixture entropy in molar units */
void fgsbml_(double *P, double *T, double *x, int * iwrk, double * rwrk, double * sbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[6]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 6; ++id) {
        result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
    }

    *sbml = result * 8.314e+07;
}


/*get mixture entropy in mass units */
void fgsbms_(double *P, double *T, double *y, int * iwrk, double * rwrk, double * sbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[6]; /* temporary storage */
    double x[6]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/31.998800; /*O2 */
    YOW += y[1]/18.015340; /*H2O */
    YOW += y[2]/16.042880; /*CH4 */
    YOW += y[3]/28.010400; /*CO */
    YOW += y[4]/44.009800; /*CO2 */
    YOW += y[5]/28.013400; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(31.998800*YOW); 
    x[1] = y[1]/(18.015340*YOW); 
    x[2] = y[2]/(16.042880*YOW); 
    x[3] = y[3]/(28.010400*YOW); 
    x[4] = y[4]/(44.009800*YOW); 
    x[5] = y[5]/(28.013400*YOW); 
    speciesEntropy(sor, tc);
    /*Perform computation in Eq 42 and 43 */
    result += x[0]*(sor[0]-log((x[0]+1e-100))-logPratio);
    result += x[1]*(sor[1]-log((x[1]+1e-100))-logPratio);
    result += x[2]*(sor[2]-log((x[2]+1e-100))-logPratio);
    result += x[3]*(sor[3]-log((x[3]+1e-100))-logPratio);
    result += x[4]*(sor[4]-log((x[4]+1e-100))-logPratio);
    result += x[5]*(sor[5]-log((x[5]+1e-100))-logPratio);
    /*Scale by R/W */
    *sbms = result * 8.314e+07 * YOW;
}


/*Returns mean gibbs free energy in molar units */
void fggbml_(double *P, double *T, double *x, int * iwrk, double * rwrk, double * gbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    double gort[6]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 6; ++id) {
        result += x[id]*(gort[id]+log((x[id]+1e-100))+logPratio);
    }

    *gbml = result * RT;
}


/*Returns mixture gibbs free energy in mass units */
void fggbms_(double *P, double *T, double *y, int * iwrk, double * rwrk, double * gbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    double gort[6]; /* temporary storage */
    double x[6]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/31.998800; /*O2 */
    YOW += y[1]/18.015340; /*H2O */
    YOW += y[2]/16.042880; /*CH4 */
    YOW += y[3]/28.010400; /*CO */
    YOW += y[4]/44.009800; /*CO2 */
    YOW += y[5]/28.013400; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(31.998800*YOW); 
    x[1] = y[1]/(18.015340*YOW); 
    x[2] = y[2]/(16.042880*YOW); 
    x[3] = y[3]/(28.010400*YOW); 
    x[4] = y[4]/(44.009800*YOW); 
    x[5] = y[5]/(28.013400*YOW); 
    gibbs(gort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(gort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(gort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(gort[2]+log((x[2]+1e-100))+logPratio);
    result += x[3]*(gort[3]+log((x[3]+1e-100))+logPratio);
    result += x[4]*(gort[4]+log((x[4]+1e-100))+logPratio);
    result += x[5]*(gort[5]+log((x[5]+1e-100))+logPratio);
    /*Scale by RT/W */
    *gbms = result * RT * YOW;
}


/*Returns mean helmholtz free energy in molar units */
void fgabml_(double *P, double *T, double *x, int * iwrk, double * rwrk, double * abml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    double aort[6]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 6; ++id) {
        result += x[id]*(aort[id]+log((x[id]+1e-100))+logPratio);
    }

    *abml = result * RT;
}


/*Returns mixture helmholtz free energy in mass units */
void fgabms_(double *P, double *T, double *y, int * iwrk, double * rwrk, double * abms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.314e+07*tT; /*R*T */
    double aort[6]; /* temporary storage */
    double x[6]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/31.998800; /*O2 */
    YOW += y[1]/18.015340; /*H2O */
    YOW += y[2]/16.042880; /*CH4 */
    YOW += y[3]/28.010400; /*CO */
    YOW += y[4]/44.009800; /*CO2 */
    YOW += y[5]/28.013400; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(31.998800*YOW); 
    x[1] = y[1]/(18.015340*YOW); 
    x[2] = y[2]/(16.042880*YOW); 
    x[3] = y[3]/(28.010400*YOW); 
    x[4] = y[4]/(44.009800*YOW); 
    x[5] = y[5]/(28.013400*YOW); 
    helmholtz(aort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(aort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(aort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(aort[2]+log((x[2]+1e-100))+logPratio);
    result += x[3]*(aort[3]+log((x[3]+1e-100))+logPratio);
    result += x[4]*(aort[4]+log((x[4]+1e-100))+logPratio);
    result += x[5]*(aort[5]+log((x[5]+1e-100))+logPratio);
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void fgwc_(double * T, double * C, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 6; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 6; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void fgwyp_(double * P, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[6]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/31.998800; /*O2 */
    YOW += y[1]/18.015340; /*H2O */
    YOW += y[2]/16.042880; /*CH4 */
    YOW += y[3]/28.010400; /*CO */
    YOW += y[4]/44.009800; /*CO2 */
    YOW += y[5]/28.013400; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.314e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]/31.998800; 
    c[1] = PWORT * y[1]/18.015340; 
    c[2] = PWORT * y[2]/16.042880; 
    c[3] = PWORT * y[3]/28.010400; 
    c[4] = PWORT * y[4]/44.009800; 
    c[5] = PWORT * y[5]/28.013400; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 6; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void fgwxp_(double * P, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[6]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.314e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 6; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 6; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void fgwyr_(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[6]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]/31.998800; 
    c[1] = 1e6 * (*rho) * y[1]/18.015340; 
    c[2] = 1e6 * (*rho) * y[2]/16.042880; 
    c[3] = 1e6 * (*rho) * y[3]/28.010400; 
    c[4] = 1e6 * (*rho) * y[4]/44.009800; 
    c[5] = 1e6 * (*rho) * y[5]/28.013400; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 6; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void fgwxr_(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * wdot)
{
    int id; /*loop counter */
    double c[6]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*31.998800; /*O2 */
    XW += x[1]*18.015340; /*H2O */
    XW += x[2]*16.042880; /*CH4 */
    XW += x[3]*28.010400; /*CO */
    XW += x[4]*44.009800; /*CO2 */
    XW += x[5]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 6; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 6; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void fgqc_(double * T, double * C, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 6; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 6; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 3; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void fgqyp_(double * P, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[6]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]/31.998800; /*O2 */
    YOW += y[1]/18.015340; /*H2O */
    YOW += y[2]/16.042880; /*CH4 */
    YOW += y[3]/28.010400; /*CO */
    YOW += y[4]/44.009800; /*CO2 */
    YOW += y[5]/28.013400; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.314e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]/31.998800; 
    c[1] = PWORT * y[1]/18.015340; 
    c[2] = PWORT * y[2]/16.042880; 
    c[3] = PWORT * y[3]/28.010400; 
    c[4] = PWORT * y[4]/44.009800; 
    c[5] = PWORT * y[5]/28.013400; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void fgqxp_(double * P, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[6]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.314e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 6; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void fgqyr_(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[6]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]/31.998800; 
    c[1] = 1e6 * (*rho) * y[1]/18.015340; 
    c[2] = 1e6 * (*rho) * y[2]/16.042880; 
    c[3] = 1e6 * (*rho) * y[3]/28.010400; 
    c[4] = 1e6 * (*rho) * y[4]/44.009800; 
    c[5] = 1e6 * (*rho) * y[5]/28.013400; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void fgqxr_(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * qdot)
{
    int id; /*loop counter */
    double c[6]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*31.998800; /*O2 */
    XW += x[1]*18.015340; /*H2O */
    XW += x[2]*16.042880; /*CH4 */
    XW += x[3]*28.010400; /*CO */
    XW += x[4]*44.009800; /*CO2 */
    XW += x[5]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 6; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 3; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the stoichiometric coefficients */
/*of the reaction mechanism. (Eq 50) */
void fgnu_(int * kdim, int * iwrk, double * rwrk, int * nuki)
{
    int id; /*loop counter */
    int kd = (*kdim); 
    /*Zero nuki */
    for (id = 0; id < 6 * 3; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: 2 CH4 + 3 O2 => 2 CO + 4 H2O */
    nuki[ 2 * kd + 0 ] = -2 ;
    nuki[ 0 * kd + 0 ] = -3 ;
    nuki[ 3 * kd + 0 ] = +2 ;
    nuki[ 1 * kd + 0 ] = +4 ;

    /*reaction 2: 2 CO + O2 => 2 CO2 */
    nuki[ 3 * kd + 1 ] = -2 ;
    nuki[ 0 * kd + 1 ] = -1 ;
    nuki[ 4 * kd + 1 ] = +2 ;

    /*reaction 3: 2 CO2 => 2 CO + O2 */
    nuki[ 4 * kd + 2 ] = -2 ;
    nuki[ 3 * kd + 2 ] = +2 ;
    nuki[ 0 * kd + 2 ] = +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void fgncf_(int * mdim, int * iwrk, double * rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < 4 * 6; ++ id) {
         ncf[id] = 0; 
    }

    /*O2 */
    ncf[ 0 * kd + 0 ] = 2; /*O */

    /*H2O */
    ncf[ 1 * kd + 1 ] = 2; /*H */
    ncf[ 1 * kd + 0 ] = 1; /*O */

    /*CH4 */
    ncf[ 2 * kd + 2 ] = 1; /*C */
    ncf[ 2 * kd + 1 ] = 4; /*H */

    /*CO */
    ncf[ 3 * kd + 2 ] = 1; /*C */
    ncf[ 3 * kd + 0 ] = 1; /*O */

    /*CO2 */
    ncf[ 4 * kd + 2 ] = 1; /*C */
    ncf[ 4 * kd + 0 ] = 2; /*O */

    /*N2 */
    ncf[ 5 * kd + 3 ] = 2; /*N */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void fgabe_(int * iwrk, double * rwrk, double * a, double * b, double * e)
{

    /*reaction 1: 2 CH4 + 3 O2 => 2 CO + 4 H2O */
    a[0] = 1.545e+14;
    b[0] = 0.5;
    e[0] = 39895;

    /*reaction 2: 2 CO + O2 => 2 CO2 */
    a[1] = 1.991e+14;
    b[1] = 0;
    e[1] = 40000;

    /*reaction 3: 2 CO2 => 2 CO + O2 */
    a[2] = 2.5e+08;
    b[2] = 0;
    e[2] = 40000;

    return;
}


/*Returns the equil constants for each reaction */
void fgeqc_(double * T, double * C, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[6]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 CH4 + 3 O2 => 2 CO + 4 H2O */
    eqcon[0] *= 1e-06; 

    /*reaction 2: 2 CO + O2 => 2 CO2 */
    eqcon[1] *= 1e+06; 

    /*reaction 3: 2 CO2 => 2 CO + O2 */
    eqcon[2] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void fgeqyp_(double * P, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[6]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 CH4 + 3 O2 => 2 CO + 4 H2O */
    eqcon[0] *= 1e-06; 

    /*reaction 2: 2 CO + O2 => 2 CO2 */
    eqcon[1] *= 1e+06; 

    /*reaction 3: 2 CO2 => 2 CO + O2 */
    eqcon[2] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void fgeqxp_(double * P, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[6]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 CH4 + 3 O2 => 2 CO + 4 H2O */
    eqcon[0] *= 1e-06; 

    /*reaction 2: 2 CO + O2 => 2 CO2 */
    eqcon[1] *= 1e+06; 

    /*reaction 3: 2 CO2 => 2 CO + O2 */
    eqcon[2] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void fgeqyr_(double * rho, double * T, double * y, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[6]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 CH4 + 3 O2 => 2 CO + 4 H2O */
    eqcon[0] *= 1e-06; 

    /*reaction 2: 2 CO + O2 => 2 CO2 */
    eqcon[1] *= 1e+06; 

    /*reaction 3: 2 CO2 => 2 CO + O2 */
    eqcon[2] *= 1e-06; 
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void fgeqxr_(double * rho, double * T, double * x, int * iwrk, double * rwrk, double * eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[6]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 CH4 + 3 O2 => 2 CO + 4 H2O */
    eqcon[0] *= 1e-06; 

    /*reaction 2: 2 CO + O2 => 2 CO2 */
    eqcon[1] *= 1e+06; 

    /*reaction 3: 2 CO2 => 2 CO + O2 */
    eqcon[2] *= 1e-06; 
}


/*compute the production rate for each species */
void productionRate(double * wdot, double * sc, double T)
{
    double qdot;

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[6];                /*Gibbs free energy */
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
    for (id = 0; id < 6; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*zero out wdot */
    for (id = 0; id < 6; ++id) {
        wdot[id] = 0.0;
    }

    /*reaction 1: 2 CH4 + 3 O2 => 2 CO + 4 H2O */
    phi_f = sc[2]*sc[2]*sc[0]*sc[0]*sc[0];
    k_f = 1e-24 * 1.545e+14*exp(0.5*tc[0]-20078/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[2] -= 2 * qdot;
    wdot[0] -= 3 * qdot;
    wdot[3] += 2 * qdot;
    wdot[1] += 4 * qdot;

    /*reaction 2: 2 CO + O2 => 2 CO2 */
    phi_f = sc[3]*sc[3]*sc[0];
    k_f = 1e-12 * 1.991e+14*exp(-20130.9/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[3] -= 2 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[4] += 2 * qdot;

    /*reaction 3: 2 CO2 => 2 CO + O2 */
    phi_f = sc[4]*sc[4];
    k_f = 1e-06 * 2.5e+08*exp(-20130.9/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot = q_f - q_r;
    wdot[4] -= 2 * qdot;
    wdot[3] += 2 * qdot;
    wdot[0] += 1 * qdot;

    return;
}


/*compute the progress rate for each reaction */
void progressRate(double * qdot, double * sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[6];                /*Gibbs free energy */
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
    for (id = 0; id < 6; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    /*reaction 1: 2 CH4 + 3 O2 => 2 CO + 4 H2O */
    phi_f = sc[2]*sc[2]*sc[0]*sc[0]*sc[0];
    k_f = 1e-24 * 1.545e+14*exp(0.5*tc[0]-20078/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[0] = q_f - q_r;

    /*reaction 2: 2 CO + O2 => 2 CO2 */
    phi_f = sc[3]*sc[3]*sc[0];
    k_f = 1e-12 * 1.991e+14*exp(-20130.9/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[1] = q_f - q_r;

    /*reaction 3: 2 CO2 => 2 CO + O2 */
    phi_f = sc[4]*sc[4];
    k_f = 1e-06 * 2.5e+08*exp(-20130.9/tc[1]);
    q_f = phi_f * k_f;
    q_r = 0.0;
    qdot[2] = q_f - q_r;

    return;
}


/*compute the equilibrium constants for each reaction */
void equilibriumConstants(double *kc, double * g_RT, double T)
{
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.314 / T;

    /*reaction 1: 2 CH4 + 3 O2 => 2 CO + 4 H2O */
    kc[0] = refC * exp((2 * g_RT[2] + 3 * g_RT[0]) - (2 * g_RT[3] + 4 * g_RT[1]));

    /*reaction 2: 2 CO + O2 => 2 CO2 */
    kc[1] = 1.0 / (refC) * exp((2 * g_RT[3] + g_RT[0]) - (2 * g_RT[4]));

    /*reaction 3: 2 CO2 => 2 CO + O2 */
    kc[2] = refC * exp((2 * g_RT[4]) - (2 * g_RT[3] + g_RT[0]));

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
        /*species 0: O2 */
        species[0] =
            -1.06394356e+03 / tc[1]
            +1.24780630e-01
            -3.78245636e+00 * tc[0]
            +1.49836708e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774591e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 1: H2O */
        species[1] =
            -3.02937267e+04 / tc[1]
            +5.04767277e+00
            -4.19864056e+00 * tc[0]
            +1.01821705e-03 * tc[1]
            -1.08673369e-06 * tc[2]
            +4.57330885e-10 * tc[3]
            -8.85989085e-14 * tc[4];
        /*species 2: CH4 */
        species[2] =
            -1.02466476e+04 / tc[1]
            +9.79117989e+00
            -5.14987613e+00 * tc[0]
            +6.83548940e-03 * tc[1]
            -8.19667665e-06 * tc[2]
            +4.03952522e-09 * tc[3]
            -8.33469780e-13 * tc[4];
        /*species 3: CO */
        species[3] =
            -1.43440860e+04 / tc[1]
            +7.11241900e-02
            -3.57953347e+00 * tc[0]
            +3.05176840e-04 * tc[1]
            -1.69469055e-07 * tc[2]
            -7.55838237e-11 * tc[3]
            +4.52212249e-14 * tc[4];
        /*species 4: CO2 */
        species[4] =
            -4.83719697e+04 / tc[1]
            -7.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
        /*species 5: N2 */
        species[5] =
            -1.02089990e+03 / tc[1]
            -6.51695000e-01
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
    } else {
        /*species 0: O2 */
        species[0] =
            -1.08845772e+03 / tc[1]
            -2.17069345e+00
            -3.28253784e+00 * tc[0]
            -7.41543770e-04 * tc[1]
            +1.26327778e-07 * tc[2]
            -1.74558796e-11 * tc[3]
            +1.08358897e-15 * tc[4];
        /*species 1: H2O */
        species[1] =
            -3.00042971e+04 / tc[1]
            -1.93277761e+00
            -3.03399249e+00 * tc[0]
            -1.08845902e-03 * tc[1]
            +2.73454197e-08 * tc[2]
            +8.08683225e-12 * tc[3]
            -8.41004960e-16 * tc[4];
        /*species 2: CH4 */
        species[2] =
            -9.46834459e+03 / tc[1]
            -1.83624665e+01
            -7.48514950e-02 * tc[0]
            -6.69547335e-03 * tc[1]
            +9.55476348e-07 * tc[2]
            -1.01910446e-10 * tc[3]
            +5.09076150e-15 * tc[4];
        /*species 3: CO */
        species[3] =
            -1.41518724e+04 / tc[1]
            -5.10350211e+00
            -2.71518561e+00 * tc[0]
            -1.03126372e-03 * tc[1]
            +1.66470962e-07 * tc[2]
            -1.91710840e-11 * tc[3]
            +1.01823858e-15 * tc[4];
        /*species 4: CO2 */
        species[4] =
            -4.87591660e+04 / tc[1]
            +1.58582223e+00
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
        /*species 5: N2 */
        species[5] =
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
        /*species 0: O2 */
        species[0] =
            -1.06394356e+03 / tc[1]
            -8.75219370e-01
            -3.78245636e+00 * tc[0]
            +1.49836708e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774591e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 1: H2O */
        species[1] =
            -3.02937267e+04 / tc[1]
            +4.04767277e+00
            -4.19864056e+00 * tc[0]
            +1.01821705e-03 * tc[1]
            -1.08673369e-06 * tc[2]
            +4.57330885e-10 * tc[3]
            -8.85989085e-14 * tc[4];
        /*species 2: CH4 */
        species[2] =
            -1.02466476e+04 / tc[1]
            +8.79117989e+00
            -5.14987613e+00 * tc[0]
            +6.83548940e-03 * tc[1]
            -8.19667665e-06 * tc[2]
            +4.03952522e-09 * tc[3]
            -8.33469780e-13 * tc[4];
        /*species 3: CO */
        species[3] =
            -1.43440860e+04 / tc[1]
            -9.28875810e-01
            -3.57953347e+00 * tc[0]
            +3.05176840e-04 * tc[1]
            -1.69469055e-07 * tc[2]
            -7.55838237e-11 * tc[3]
            +4.52212249e-14 * tc[4];
        /*species 4: CO2 */
        species[4] =
            -4.83719697e+04 / tc[1]
            -8.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
        /*species 5: N2 */
        species[5] =
            -1.02089990e+03 / tc[1]
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
    } else {
        /*species 0: O2 */
        species[0] =
            -1.08845772e+03 / tc[1]
            -3.17069345e+00
            -3.28253784e+00 * tc[0]
            -7.41543770e-04 * tc[1]
            +1.26327778e-07 * tc[2]
            -1.74558796e-11 * tc[3]
            +1.08358897e-15 * tc[4];
        /*species 1: H2O */
        species[1] =
            -3.00042971e+04 / tc[1]
            -2.93277761e+00
            -3.03399249e+00 * tc[0]
            -1.08845902e-03 * tc[1]
            +2.73454197e-08 * tc[2]
            +8.08683225e-12 * tc[3]
            -8.41004960e-16 * tc[4];
        /*species 2: CH4 */
        species[2] =
            -9.46834459e+03 / tc[1]
            -1.93624665e+01
            -7.48514950e-02 * tc[0]
            -6.69547335e-03 * tc[1]
            +9.55476348e-07 * tc[2]
            -1.01910446e-10 * tc[3]
            +5.09076150e-15 * tc[4];
        /*species 3: CO */
        species[3] =
            -1.41518724e+04 / tc[1]
            -6.10350211e+00
            -2.71518561e+00 * tc[0]
            -1.03126372e-03 * tc[1]
            +1.66470962e-07 * tc[2]
            -1.91710840e-11 * tc[3]
            +1.01823858e-15 * tc[4];
        /*species 4: CO2 */
        species[4] =
            -4.87591660e+04 / tc[1]
            +5.85822230e-01
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
        /*species 5: N2 */
        species[5] =
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
        /*species 0: O2 */
        species[0] =
            +2.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 1: H2O */
        species[1] =
            +3.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 2: CH4 */
        species[2] =
            +4.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 3: CO */
        species[3] =
            +2.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 4: CO2 */
        species[4] =
            +1.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 5: N2 */
        species[5] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
    } else {
        /*species 0: O2 */
        species[0] =
            +2.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 1: H2O */
        species[1] =
            +2.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 2: CH4 */
        species[2] =
            -9.25148505e-01
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 3: CO */
        species[3] =
            +1.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 4: CO2 */
        species[4] =
            +2.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 5: N2 */
        species[5] =
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
        /*species 0: O2 */
        species[0] =
            +3.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 1: H2O */
        species[1] =
            +4.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 2: CH4 */
        species[2] =
            +5.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 3: CO */
        species[3] =
            +3.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 4: CO2 */
        species[4] =
            +2.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 5: N2 */
        species[5] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
    } else {
        /*species 0: O2 */
        species[0] =
            +3.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 1: H2O */
        species[1] =
            +3.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 2: CH4 */
        species[2] =
            +7.48514950e-02
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 3: CO */
        species[3] =
            +2.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 4: CO2 */
        species[4] =
            +3.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 5: N2 */
        species[5] =
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
        /*species 0: O2 */
        species[0] =
            +2.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 / tc[1];
        /*species 1: H2O */
        species[1] =
            +3.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 / tc[1];
        /*species 2: CH4 */
        species[2] =
            +4.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 / tc[1];
        /*species 3: CO */
        species[3] =
            +2.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 / tc[1];
        /*species 4: CO2 */
        species[4] =
            +1.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 / tc[1];
        /*species 5: N2 */
        species[5] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 / tc[1];
    } else {
        /*species 0: O2 */
        species[0] =
            +2.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 / tc[1];
        /*species 1: H2O */
        species[1] =
            +2.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 / tc[1];
        /*species 2: CH4 */
        species[2] =
            -9.25148505e-01
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 / tc[1];
        /*species 3: CO */
        species[3] =
            +1.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 / tc[1];
        /*species 4: CO2 */
        species[4] =
            +2.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 / tc[1];
        /*species 5: N2 */
        species[5] =
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
        /*species 0: O2 */
        species[0] =
            +3.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 / tc[1];
        /*species 1: H2O */
        species[1] =
            +4.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 / tc[1];
        /*species 2: CH4 */
        species[2] =
            +5.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 / tc[1];
        /*species 3: CO */
        species[3] =
            +3.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 / tc[1];
        /*species 4: CO2 */
        species[4] =
            +2.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 / tc[1];
        /*species 5: N2 */
        species[5] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 / tc[1];
    } else {
        /*species 0: O2 */
        species[0] =
            +3.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 / tc[1];
        /*species 1: H2O */
        species[1] =
            +3.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 / tc[1];
        /*species 2: CH4 */
        species[2] =
            +7.48514950e-02
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 / tc[1];
        /*species 3: CO */
        species[3] =
            +2.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 / tc[1];
        /*species 4: CO2 */
        species[4] =
            +3.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 / tc[1];
        /*species 5: N2 */
        species[5] =
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
        /*species 0: O2 */
        species[0] =
            +3.78245636e+00 * tc[0]
            -2.99673416e-03 * tc[1]
            +4.92365101e-06 * tc[2]
            -3.22709836e-09 * tc[3]
            +8.10932092e-13 * tc[4]
            +3.65767573e+00 ;
        /*species 1: H2O */
        species[1] =
            +4.19864056e+00 * tc[0]
            -2.03643410e-03 * tc[1]
            +3.26020105e-06 * tc[2]
            -1.82932354e-09 * tc[3]
            +4.42994543e-13 * tc[4]
            -8.49032208e-01 ;
        /*species 2: CH4 */
        species[2] =
            +5.14987613e+00 * tc[0]
            -1.36709788e-02 * tc[1]
            +2.45900299e-05 * tc[2]
            -1.61581009e-08 * tc[3]
            +4.16734890e-12 * tc[4]
            -4.64130376e+00 ;
        /*species 3: CO */
        species[3] =
            +3.57953347e+00 * tc[0]
            -6.10353680e-04 * tc[1]
            +5.08407165e-07 * tc[2]
            +3.02335295e-10 * tc[3]
            -2.26106125e-13 * tc[4]
            +3.50840928e+00 ;
        /*species 4: CO2 */
        species[4] =
            +2.35677352e+00 * tc[0]
            +8.98459677e-03 * tc[1]
            -3.56178134e-06 * tc[2]
            +8.19730073e-10 * tc[3]
            -3.59248870e-14 * tc[4]
            +9.90105222e+00 ;
        /*species 5: N2 */
        species[5] =
            +3.29867700e+00 * tc[0]
            +1.40824040e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213500e-13 * tc[4]
            +3.95037200e+00 ;
    } else {
        /*species 0: O2 */
        species[0] =
            +3.28253784e+00 * tc[0]
            +1.48308754e-03 * tc[1]
            -3.78983334e-07 * tc[2]
            +6.98235183e-11 * tc[3]
            -5.41794485e-15 * tc[4]
            +5.45323129e+00 ;
        /*species 1: H2O */
        species[1] =
            +3.03399249e+00 * tc[0]
            +2.17691804e-03 * tc[1]
            -8.20362590e-08 * tc[2]
            -3.23473290e-11 * tc[3]
            +4.20502480e-15 * tc[4]
            +4.96677010e+00 ;
        /*species 2: CH4 */
        species[2] =
            +7.48514950e-02 * tc[0]
            +1.33909467e-02 * tc[1]
            -2.86642905e-06 * tc[2]
            +4.07641783e-10 * tc[3]
            -2.54538075e-14 * tc[4]
            +1.84373180e+01 ;
        /*species 3: CO */
        species[3] =
            +2.71518561e+00 * tc[0]
            +2.06252743e-03 * tc[1]
            -4.99412886e-07 * tc[2]
            +7.66843360e-11 * tc[3]
            -5.09119290e-15 * tc[4]
            +7.81868772e+00 ;
        /*species 4: CO2 */
        species[4] =
            +3.85746029e+00 * tc[0]
            +4.41437026e-03 * tc[1]
            -1.10740702e-06 * tc[2]
            +1.74496729e-10 * tc[3]
            -1.18021041e-14 * tc[4]
            +2.27163806e+00 ;
        /*species 5: N2 */
        species[5] =
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
    wt[0] = 31.998800; /*O2 */
    wt[1] = 18.015340; /*H2O */
    wt[2] = 16.042880; /*CH4 */
    wt[3] = 28.010400; /*CO */
    wt[4] = 44.009800; /*CO2 */
    wt[5] = 28.013400; /*N2 */

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
    fgubms_(&tmin, y, iwrk, rwrk, &emin);
    fgubms_(&tmax, y, iwrk, rwrk, &emax);
    if (ein < emin) {
        /*Linear Extrapolation below tmin */
        fgcvbs_(&tmin, y, iwrk, rwrk, &cv);
        *t = tmin - (emin-ein)/cv;
        return 1;
    }
    if (ein > emax) {
        /*Linear Extrapolation above tmax */
        fgcvbs_(&tmax, y, iwrk, rwrk, &cv);
        *t = tmax - (emax-ein)/cv;
        return 1;
    }
    t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);
    for (i = 0; i < maxiter; ++i) {
        fgubms_(&t1,y,iwrk,rwrk,&e1);
        fgcvbs_(&t1,y,iwrk,rwrk,&cv);
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
    y[0] = phi[0]*31.998800;   XW += y[0]; /*O2 */
    y[1] = phi[1]*18.015340;   XW += y[1]; /*H2O */
    y[2] = phi[2]*16.042880;   XW += y[2]; /*CH4 */
    y[3] = phi[3]*28.010400;   XW += y[3]; /*CO */
    y[4] = phi[4]*44.009800;   XW += y[4]; /*CO2 */
    y[5] = phi[5]*28.013400;   XW += y[5]; /*N2 */
    for (id = 0; id < 6; ++id) {
        y[id] = y[id]/XW;
    }

    return;
}


/*convert y[species] (mass fracs) to phi[species] (specific mole num) */
void feytphi_(double * y, int * iwrk, double * rwrk, double * phi)
{
    phi[0] = y[0]/ 3.19988000e-02; /*O2 (wt in kg) */
    phi[1] = y[1]/ 1.80153400e-02; /*H2O (wt in kg) */
    phi[2] = y[2]/ 1.60428800e-02; /*CH4 (wt in kg) */
    phi[3] = y[3]/ 2.80104000e-02; /*CO (wt in kg) */
    phi[4] = y[4]/ 4.40098000e-02; /*CO2 (wt in kg) */
    phi[5] = y[5]/ 2.80134000e-02; /*N2 (wt in kg) */

    return;
}


/*reverse of ytcr, useful for rate computations */
void fectyr_(double * c, double * rho, int * iwrk, double * rwrk, double * y)
{
    y[0] = c[0] * 31.998800 / (*rho); 
    y[1] = c[1] * 18.015340 / (*rho); 
    y[2] = c[2] * 16.042880 / (*rho); 
    y[3] = c[3] * 28.010400 / (*rho); 
    y[4] = c[4] * 44.009800 / (*rho); 
    y[5] = c[5] * 28.013400 / (*rho); 

    return;
}


/*ddebdf compatible right hand side of CV burner */
/*rwrk[0] and rwrk[1] should contain rho and ene respectively */
/*working variable phi contains specific mole numbers */
void fecvrhs_(double * time, double * phi, double * phidot, double * rwrk, int * iwrk)
{
    double rho,ene; /*CV Parameters */
    double y[6], wdot[6]; /*temporary storage */
    int i; /*Loop counter */
    double temperature,pressure; /*temporary var */
    rho = rwrk[0];
    ene = rwrk[1];
    fephity_(phi, iwrk, rwrk, y);
    feeytt_(&ene, y, iwrk, rwrk, &temperature);
    fgpy_(&rho, &temperature,  y, iwrk, rwrk, &pressure);
    fgwyp_(&pressure, &temperature,  y, iwrk, rwrk, wdot);
    for (i=0; i<6; ++i) phidot[i] = wdot[i] / (rho/1000.0); 

    return;
}


/*returns the dimensionality of the cv burner (number of species) */
int fecvdim_()
{
    return 6;
}


/*ddebdf compatible right hand side of ZND solver */
/*rwrk[0] : scaling factor for pressure */
/*rwrk[1] : preshock density (g/cc)  */
/*rwrk[2] : detonation velocity (cm/s)  */
/*solution vector: [P; rho; y0 ... ylast]  */
void fezndrhs_(double * time, double * z, double * zdot, double * rwrk, int * iwrk)
{
    double psc,rho1,udet; /*ZND Parameters */
    double wt[6], hms[6], wdot[6]; /*temporary storage */
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

    fgmmwy_(y, 0, 0, &wtm);

    T = p * wtm / rho / ru;

    uvel = (rho1 * udet)/ rho;

    fgcpbs_(&T, y, 0, 0, &cp);
    fgcvbs_(&T, y, 0, 0, &cv);
    gam = cp/cv;

    son = sqrt(fabs(gam*ru*T/wtm));
    xm = uvel/son;

    fghms_(&T, 0, 0, hms);
    fgwt_(0, 0, wt);
    fgwyp_(&p, &T, y, 0, 0, wdot);

    sum = 0.0;
    for (i=0; i<6; ++i) {
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
    return 9;
}


/*returns the name of the source mechanism file  */
char* femechfile_()
{
    return "";
}


/*returns the species number */
int fesymnum_(const char* s1)
{
    if (strcmp(s1, "O2")==0) return 0; 
    if (strcmp(s1, "H2O")==0) return 1; 
    if (strcmp(s1, "CH4")==0) return 2; 
    if (strcmp(s1, "CO")==0) return 3; 
    if (strcmp(s1, "CO2")==0) return 4; 
    if (strcmp(s1, "N2")==0) return 5; 
    /*species name not found */
    return -1;
}


/*returns the species name */
char* fesymname_(int sn)
{
    if (sn==0) return "O2"; 
    if (sn==1) return "H2O"; 
    if (sn==2) return "CH4"; 
    if (sn==3) return "CO"; 
    if (sn==4) return "CO2"; 
    if (sn==5) return "N2"; 
    /*species name not found */
    return "NOTFOUND";
}

/* End of file  */

--------------030905060905030504010201--

