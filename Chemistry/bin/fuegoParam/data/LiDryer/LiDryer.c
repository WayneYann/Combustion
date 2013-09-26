
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
#define DWDOT DWDOT
#define VCKHMS VCKHMS
#define VCKPY VCKPY
#define VCKWYR VCKWYR
#define VCKYTX VCKYTX
#define GET_T_GIVEN_EY GET_T_GIVEN_EY
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
#define DWDOT dwdot
#define VCKHMS vckhms
#define VCKPY vckpy
#define VCKWYR vckwyr
#define VCKYTX vckytx
#define GET_T_GIVEN_EY get_t_given_ey
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
#define DWDOT dwdot_
#define VCKHMS vckhms_
#define VCKPY vckpy_
#define VCKWYR vckwyr_
#define VCKYTX vckytx_
#define GET_T_GIVEN_EY get_t_given_ey_
#endif

/*function declarations */
void molecularWeight(double * restrict  wt);
void gibbs(double * restrict  species, double * restrict  tc);
void helmholtz(double * restrict  species, double * restrict  tc);
void speciesInternalEnergy(double * restrict  species, double * restrict  tc);
void speciesEnthalpy(double * restrict  species, double * restrict  tc);
void speciesEntropy(double * restrict  species, double * restrict  tc);
void cp_R(double * restrict  species, double * restrict  tc);
void cv_R(double * restrict  species, double * restrict  tc);
void equilibriumConstants(double * restrict  kc, double * restrict  g_RT, double T);
void productionRate(double * restrict  wdot, double * restrict  sc, double T);
void progressRate(double * restrict  qdot, double * restrict  speciesConc, double T);
void progressRateFR(double * restrict  q_f, double * restrict  q_r, double * restrict  speciesConc, double T);
void CKINIT();
void CKINDX(int * iwrk, double * restrict rwrk, int * mm, int * kk, int * ii, int * nfit );
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double * restrict  rval, int * kerr, int lenline);
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * restrict  rval, int * kerr, int lenline, int lenkray);
void CKSYME(int * kname, int * lenkname);
void CKSYMS(int * kname, int * lenkname);
void CKRP(int * ickwrk, double * restrict  rckwrk, double * restrict  ru, double * restrict  ruc, double * restrict  pa);
void CKPX(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict rwrk, double * restrict  P);
void CKPY(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  P);
void CKPC(double * restrict  rho, double * restrict  T, double * restrict  c, int * iwrk, double * restrict rwrk, double * restrict  P);
void CKRHOX(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict rwrk, double * restrict  rho);
void CKRHOY(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  rho);
void CKRHOC(double * restrict  P, double * restrict  T, double * restrict  c, int * iwrk, double * restrict rwrk, double * restrict  rho);
void CKWT(int * iwrk, double * restrict rwrk, double * restrict  wt);
void CKMMWY(double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  wtm);
void CKMMWX(double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  wtm);
void CKMMWC(double * restrict  c, int * iwrk, double * restrict  rwrk, double * restrict  wtm);
void CKYTX(double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  x);
void CKYTCP(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  c);
void CKYTCR(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  c);
void CKXTY(double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  y);
void CKXTCP(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  c);
void CKXTCR(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  c);
void CKCTX(double * restrict  c, int * iwrk, double * restrict  rwrk, double * restrict  x);
void CKCTY(double * restrict  c, int * iwrk, double * restrict  rwrk, double * restrict  y);
void CKCPOR(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  cpor);
void CKHORT(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  hort);
void CKSOR(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  sor);
void CKCVML(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  cvml);
void CKCPML(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  cvml);
void CKUML(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  uml);
void CKHML(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  uml);
void CKGML(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  gml);
void CKAML(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  aml);
void CKSML(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  sml);
void CKCVMS(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  cvms);
void CKCPMS(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  cvms);
void CKUMS(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  ums);
void CKHMS(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  ums);
void CKGMS(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  gms);
void CKAMS(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  ams);
void CKSMS(double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  sms);
void CKCPBL(double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  cpbl);
void CKCPBS(double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  cpbs);
void CKCVBL(double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  cpbl);
void CKCVBS(double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  cpbs);
void CKHBML(double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  hbml);
void CKHBMS(double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  hbms);
void CKUBML(double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  ubml);
void CKUBMS(double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  ubms);
void CKSBML(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  sbml);
void CKSBMS(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  sbms);
void CKGBML(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  gbml);
void CKGBMS(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  gbms);
void CKABML(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  abml);
void CKABMS(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  abms);
void CKWC(double * restrict  T, double * restrict  C, int * iwrk, double * restrict rwrk, double * restrict  wdot);
void CKWYP(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  wdot);
void CKWXP(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict rwrk, double * restrict  wdot);
void CKWYR(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  wdot);
void CKWXR(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict rwrk, double * restrict  wdot);
void CKQC(double * restrict  T, double * restrict  C, int * iwrk, double * restrict rwrk, double * restrict  qdot);
void CKKFKR(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict rwrk, double * restrict  q_f, double * restrict  q_r);
void CKQYP(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  qdot);
void CKQXP(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict rwrk, double * restrict  qdot);
void CKQYR(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  qdot);
void CKQXR(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict rwrk, double * restrict  qdot);
void CKNU(int * kdim, int * iwrk, double * restrict rwrk, int * nuki);
void CKNCF(int * mdim, int * iwrk, double * restrict rwrk, int * ncf);
void CKABE(int * iwrk, double * restrict rwrk, double * restrict  a, double * restrict  b, double * restrict  e );
void CKEQC(double * restrict  T, double * restrict  C , int * iwrk, double * restrict rwrk, double * restrict  eqcon );
void CKEQYP(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  eqcon);
void CKEQXP(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict rwrk, double * restrict  eqcon);
void CKEQYR(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  eqcon);
void CKEQXR(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict rwrk, double * restrict  eqcon);
void DWDOT(double * restrict  J, double * restrict  sc, double * T, int * consP);
void aJacobian(double * restrict J, double * restrict sc, double T, int consP);
void dcvpRdT(double * restrict  species, double * restrict  tc);
void GET_T_GIVEN_EY(double * restrict  e, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  t, int *ierr);
/*vector version */
void vproductionRate(int npt, double * restrict wdot, double * restrict c, double * restrict T);
void VCKHMS(int * restrict np, double * restrict  T, int * iwrk, double * restrict  rwrk, double * restrict  ums);
void VCKPY(int * restrict np, double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict rwrk, double * restrict  P);
void VCKWYR(int * restrict np, double * restrict rho, double * restrict T,
            double * restrict y, int * restrict iwrk, double * restrict rwrk,
            double * restrict wdot);
void VCKYTX(int * restrict np, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  x);

/* Inverse molecular weights */
static const double imw[9] = {
    1.0 / 2.015940,  /*H2 */
    1.0 / 31.998800,  /*O2 */
    1.0 / 18.015340,  /*H2O */
    1.0 / 1.007970,  /*H */
    1.0 / 15.999400,  /*O */
    1.0 / 17.007370,  /*OH */
    1.0 / 33.006770,  /*HO2 */
    1.0 / 34.014740,  /*H2O2 */
    1.0 / 28.013400};  /*N2 */

struct ReactionData {
    double activation_units, prefactor_units, phase_units;
    double fwd_A,fwd_beta,fwd_Ea;
    double low_A,low_beta,low_Ea;
    double rev_A,rev_beta,rev_Ea;
    int is_PD, troe_len, sri_len;
    double troe_a,troe_Ts, troe_Tss, troe_Tsss;
    double sri_a, sri_b, sri_c, sri_d, sri_e;
};

static struct ReactionData R[21], R_DEF[21];

struct ReactionData* GetReactionData(int id)
{
    if (id<0 || id>=21) {
        printf("GetReactionData: Bad reaction id = \%d",id);
        abort();
    };
    return &(R[id]);
}

struct ReactionData* GetDefaultReactionData(int id)
{
    if (id<0 || id>=21) {
        printf("GetDefaultReactionData: Bad reaction id = \%d",id);
        abort();
    };
    return &(R_DEF[id]);
}


void SetReactionData(int id, const struct ReactionData * rhs)
{
    if (id<0 || id>=175) {
        printf("SetReactionData: Bad reaction id = %d",id);
        abort();
    }
    struct ReactionData* r = &(R[id]);

    r->fwd_A = rhs->fwd_A;
    r->fwd_beta = rhs->fwd_beta;
    r->fwd_Ea = rhs->fwd_Ea;

    r->low_A = rhs->low_A;
    r->low_beta = rhs->low_beta;
    r->low_Ea = rhs->low_Ea;

    r->rev_A = rhs->rev_A;
    r->rev_beta = rhs->rev_beta;
    r->rev_Ea = rhs->rev_Ea;

    r->is_PD = rhs->is_PD;
    r->troe_len = rhs->troe_len;
    r->sri_len = rhs->sri_len;

    r->troe_a = rhs->troe_a;
    r->troe_Ts = rhs->troe_Ts;
    r->troe_Tss = rhs->troe_Tss;
    r->troe_Tsss = rhs->troe_Tsss;

    r->sri_a = rhs->sri_a;
    r->sri_b = rhs->sri_b;
    r->sri_c = rhs->sri_c;
    r->sri_d = rhs->sri_d;
    r->sri_e = rhs->sri_e;
}

/* Initializes static database */
void CKINIT()
{
    R[0].fwd_A     = 3.547e+15;
    R[0].fwd_beta  = -0.406;
    R[0].fwd_Ea    = 16599;
    R[0].prefactor_units = 1e-06;
    R[0].activation_units = 0.503217;
    R[0].phase_units = 1e-12;

    R_DEF[0] = R[0];

    R[1].fwd_A     = 50800;
    R[1].fwd_beta  = 2.67;
    R[1].fwd_Ea    = 6290;
    R[1].prefactor_units = 1e-06;
    R[1].activation_units = 0.503217;
    R[1].phase_units = 1e-12;

    R_DEF[1] = R[1];

    R[2].fwd_A     = 2.16e+08;
    R[2].fwd_beta  = 1.51;
    R[2].fwd_Ea    = 3430;
    R[2].prefactor_units = 1e-06;
    R[2].activation_units = 0.503217;
    R[2].phase_units = 1e-12;

    R_DEF[2] = R[2];

    R[3].fwd_A     = 2.97e+06;
    R[3].fwd_beta  = 2.02;
    R[3].fwd_Ea    = 13400;
    R[3].prefactor_units = 1e-06;
    R[3].activation_units = 0.503217;
    R[3].phase_units = 1e-12;

    R_DEF[3] = R[3];

    R[4].fwd_A     = 4.577e+19;
    R[4].fwd_beta  = -1.4;
    R[4].fwd_Ea    = 104380;
    R[4].prefactor_units = 1e-06;
    R[4].activation_units = 0.503217;
    R[4].phase_units = 1e-6;

    R_DEF[4] = R[4];

    R[5].fwd_A     = 6.165e+15;
    R[5].fwd_beta  = -0.5;
    R[5].fwd_Ea    = 0;
    R[5].prefactor_units = 1e-12;
    R[5].activation_units = 0.503217;
    R[5].phase_units = 1e-12;

    R_DEF[5] = R[5];

    R[6].fwd_A     = 4.714e+18;
    R[6].fwd_beta  = -1;
    R[6].fwd_Ea    = 0;
    R[6].prefactor_units = 1e-12;
    R[6].activation_units = 0.503217;
    R[6].phase_units = 1e-12;

    R_DEF[6] = R[6];

    R[7].fwd_A     = 3.8e+22;
    R[7].fwd_beta  = -2;
    R[7].fwd_Ea    = 0;
    R[7].prefactor_units = 1e-12;
    R[7].activation_units = 0.503217;
    R[7].phase_units = 1e-12;

    R_DEF[7] = R[7];

    R[8].fwd_A     = 1.475e+12;
    R[8].fwd_beta  = 0.6;
    R[8].fwd_Ea    = 0;
    R[8].low_A     = 6.366e+20;
    R[8].low_beta  = -1.72;
    R[8].low_Ea    = 524.8;
    R[8].troe_a    = 0.8;
    R[8].troe_Tsss = 1e-30;
    R[8].troe_Ts   = 1e+30;
    R[8].troe_len  = 3;
    R[8].is_PD = 1;
    R[8].prefactor_units = 1e-06;
    R[8].activation_units = 0.503217;
    R[8].phase_units = 1e-12;

    R_DEF[8] = R[8];

    R[9].fwd_A     = 1.66e+13;
    R[9].fwd_beta  = 0;
    R[9].fwd_Ea    = 823;
    R[9].prefactor_units = 1e-06;
    R[9].activation_units = 0.503217;
    R[9].phase_units = 1e-12;

    R_DEF[9] = R[9];

    R[10].fwd_A     = 7.079e+13;
    R[10].fwd_beta  = 0;
    R[10].fwd_Ea    = 295;
    R[10].prefactor_units = 1e-06;
    R[10].activation_units = 0.503217;
    R[10].phase_units = 1e-12;

    R_DEF[10] = R[10];

    R[11].fwd_A     = 3.25e+13;
    R[11].fwd_beta  = 0;
    R[11].fwd_Ea    = 0;
    R[11].prefactor_units = 1e-06;
    R[11].activation_units = 0.503217;
    R[11].phase_units = 1e-12;

    R_DEF[11] = R[11];

    R[12].fwd_A     = 2.89e+13;
    R[12].fwd_beta  = 0;
    R[12].fwd_Ea    = -497;
    R[12].prefactor_units = 1e-06;
    R[12].activation_units = 0.503217;
    R[12].phase_units = 1e-12;

    R_DEF[12] = R[12];

    R[13].fwd_A     = 4.2e+14;
    R[13].fwd_beta  = 0;
    R[13].fwd_Ea    = 11982;
    R[13].prefactor_units = 1e-06;
    R[13].activation_units = 0.503217;
    R[13].phase_units = 1e-12;

    R_DEF[13] = R[13];

    R[14].fwd_A     = 1.3e+11;
    R[14].fwd_beta  = 0;
    R[14].fwd_Ea    = -1629.3;
    R[14].prefactor_units = 1e-06;
    R[14].activation_units = 0.503217;
    R[14].phase_units = 1e-12;

    R_DEF[14] = R[14];

    R[15].fwd_A     = 2.951e+14;
    R[15].fwd_beta  = 0;
    R[15].fwd_Ea    = 48430;
    R[15].low_A     = 1.202e+17;
    R[15].low_beta  = 0;
    R[15].low_Ea    = 45500;
    R[15].troe_a    = 0.5;
    R[15].troe_Tsss = 1e-30;
    R[15].troe_Ts   = 1e+30;
    R[15].troe_len  = 3;
    R[15].is_PD = 1;
    R[15].prefactor_units = 1;
    R[15].activation_units = 0.503217;
    R[15].phase_units = 1e-6;

    R_DEF[15] = R[15];

    R[16].fwd_A     = 2.41e+13;
    R[16].fwd_beta  = 0;
    R[16].fwd_Ea    = 3970;
    R[16].prefactor_units = 1e-06;
    R[16].activation_units = 0.503217;
    R[16].phase_units = 1e-12;

    R_DEF[16] = R[16];

    R[17].fwd_A     = 4.82e+13;
    R[17].fwd_beta  = 0;
    R[17].fwd_Ea    = 7950;
    R[17].prefactor_units = 1e-06;
    R[17].activation_units = 0.503217;
    R[17].phase_units = 1e-12;

    R_DEF[17] = R[17];

    R[18].fwd_A     = 9.55e+06;
    R[18].fwd_beta  = 2;
    R[18].fwd_Ea    = 3970;
    R[18].prefactor_units = 1e-06;
    R[18].activation_units = 0.503217;
    R[18].phase_units = 1e-12;

    R_DEF[18] = R[18];

    R[19].fwd_A     = 1e+12;
    R[19].fwd_beta  = 0;
    R[19].fwd_Ea    = 0;
    R[19].prefactor_units = 1e-06;
    R[19].activation_units = 0.503217;
    R[19].phase_units = 1e-12;

    R_DEF[19] = R[19];

    R[20].fwd_A     = 5.8e+14;
    R[20].fwd_beta  = 0;
    R[20].fwd_Ea    = 9557;
    R[20].prefactor_units = 1e-06;
    R[20].activation_units = 0.503217;
    R[20].phase_units = 1e-12;

    R_DEF[20] = R[20];

}


/*A few mechanism parameters */
void CKINDX(int * iwrk, double * restrict  rwrk, int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 3;
    *kk = 9;
    *ii = 21;
    *nfit = -1; /*Why do you need this anyway ?  */
}


/* ckxnum... for parsing strings  */
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double * restrict  rval, int * kerr, int lenline )
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
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * restrict  rval, int * kerr, int lenline, int lenkray)
{
    /*Not done yet ... */
}


/* Returns the char strings of element names */
void CKSYME(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*3; i++) {
        kname[i] = ' ';
    }

    /* H  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = ' ';

    /* O  */
    kname[ 1*lenkname + 0 ] = 'O';
    kname[ 1*lenkname + 1 ] = ' ';

    /* N  */
    kname[ 2*lenkname + 0 ] = 'N';
    kname[ 2*lenkname + 1 ] = ' ';

}


/* Returns the char strings of species names */
void CKSYMS(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*9; i++) {
        kname[i] = ' ';
    }

    /* H2  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = '2';
    kname[ 0*lenkname + 2 ] = ' ';

    /* O2  */
    kname[ 1*lenkname + 0 ] = 'O';
    kname[ 1*lenkname + 1 ] = '2';
    kname[ 1*lenkname + 2 ] = ' ';

    /* H2O  */
    kname[ 2*lenkname + 0 ] = 'H';
    kname[ 2*lenkname + 1 ] = '2';
    kname[ 2*lenkname + 2 ] = 'O';
    kname[ 2*lenkname + 3 ] = ' ';

    /* H  */
    kname[ 3*lenkname + 0 ] = 'H';
    kname[ 3*lenkname + 1 ] = ' ';

    /* O  */
    kname[ 4*lenkname + 0 ] = 'O';
    kname[ 4*lenkname + 1 ] = ' ';

    /* OH  */
    kname[ 5*lenkname + 0 ] = 'O';
    kname[ 5*lenkname + 1 ] = 'H';
    kname[ 5*lenkname + 2 ] = ' ';

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
void CKRP(int * ickwrk, double * restrict  rckwrk, double * restrict  ru, double * restrict  ruc, double * restrict  pa)
{
     *ru  = 8.31451e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}


/*Compute P = rhoRT/W(x) */
void CKPX(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  P)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void CKPY(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void VCKPY(int * restrict np, double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  P)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<9; n++) {
        for (int i=0; i<(*np); i++) {
            YOW[i] += y[n*(*np)+i] * imw[n];
        }
    }

    for (int i=0; i<(*np); i++) {
        P[i] = rho[i] * 8.31451e+07 * T[i] * YOW[i]; /*P = rho*R*T/W */
    }

    return;
}


/*Compute P = rhoRT/W(c) */
void CKPC(double * restrict  rho, double * restrict  T, double * restrict  c, int * iwrk, double * restrict  rwrk, double * restrict  P)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*18.015340; /*H2O */
    W += c[3]*1.007970; /*H */
    W += c[4]*15.999400; /*O */
    W += c[5]*17.007370; /*OH */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */

    for (id = 0; id < 9; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31451e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  rho)
{
    double YOW = 0;
    double tmp[9];

    for (int i = 0; i < 9; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
    {
        YOW += tmp[i];
    }

    *rho = *P / (83145100 * (*T) * YOW);/*rho = P*W/(R*T) */
    return;
}


/*Compute rho = P*W(c)/(R*T) */
void CKRHOC(double * restrict  P, double * restrict  T, double * restrict  c, int * iwrk, double * restrict  rwrk, double * restrict  rho)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*18.015340; /*H2O */
    W += c[3]*1.007970; /*H */
    W += c[4]*15.999400; /*O */
    W += c[5]*17.007370; /*OH */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */

    for (id = 0; id < 9; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.31451e+07); /*rho = PW/(R*T) */

    return;
}


/*get molecular weight for all species */
void CKWT(int * iwrk, double * restrict  rwrk, double * restrict  wt)
{
    molecularWeight(wt);
}


/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWY(double * restrict y, int * iwrk, double * restrict  rwrk, double * restrict  wtm)
{
    double YOW = 0;
    double tmp[9];

    for (int i = 0; i < 9; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
    {
        YOW += tmp[i];
    }

    *wtm = 1.0 / YOW;
    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWX(double * restrict x, int * iwrk, double * restrict  rwrk, double * restrict  wtm)
{
    double XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    *wtm = XW;

    return;
}


/*given c[species]: molar concentration */
/*returns mean molecular weight (gm/mole) */
void CKMMWC(double * restrict c, int * iwrk, double * restrict  rwrk, double * restrict  wtm)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*18.015340; /*H2O */
    W += c[3]*1.007970; /*H */
    W += c[4]*15.999400; /*O */
    W += c[5]*17.007370; /*OH */
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
void CKYTX(double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  x)
{
    double YOW = 0;
    double tmp[9];

    for (int i = 0; i < 9; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
    {
        YOW += tmp[i];
    }

    double YOWINV = 1.0/YOW;

    for (int i = 0; i < 9; i++)
    {
        x[i] = y[i]*imw[i]*YOWINV;
    }
    return;
}


/*convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs) */
void VCKYTX(int * restrict np, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  x)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<9; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<9; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCP(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  c)
{
    double YOW = 0;
    double PWORT;

    /*Compute inverse of mean molecular wt first */
    for (int i = 0; i < 9; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 9; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCR(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  c)
{
    for (int i = 0; i < 9; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
void CKXTY(double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    /*Now compute conversion */
    y[0] = x[0]*2.015940/XW; 
    y[1] = x[1]*31.998800/XW; 
    y[2] = x[2]*18.015340/XW; 
    y[3] = x[3]*1.007970/XW; 
    y[4] = x[4]*15.999400/XW; 
    y[5] = x[5]*17.007370/XW; 
    y[6] = x[6]*33.006770/XW; 
    y[7] = x[7]*34.014740/XW; 
    y[8] = x[8]*28.013400/XW; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*PORT;
    }

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCR(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  c)
{
    int id; /*loop counter */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
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
void CKCTX(double * restrict  c, int * iwrk, double * restrict  rwrk, double * restrict  x)
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
void CKCTY(double * restrict  c, int * iwrk, double * restrict  rwrk, double * restrict  y)
{
    double CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*2.015940; /*H2 */
    CW += c[1]*31.998800; /*O2 */
    CW += c[2]*18.015340; /*H2O */
    CW += c[3]*1.007970; /*H */
    CW += c[4]*15.999400; /*O */
    CW += c[5]*17.007370; /*OH */
    CW += c[6]*33.006770; /*HO2 */
    CW += c[7]*34.014740; /*H2O2 */
    CW += c[8]*28.013400; /*N2 */
    /*Now compute conversion */
    y[0] = c[0]*2.015940/CW; 
    y[1] = c[1]*31.998800/CW; 
    y[2] = c[2]*18.015340/CW; 
    y[3] = c[3]*1.007970/CW; 
    y[4] = c[4]*15.999400/CW; 
    y[5] = c[5]*17.007370/CW; 
    y[6] = c[6]*33.006770/CW; 
    y[7] = c[7]*34.014740/CW; 
    y[8] = c[8]*28.013400/CW; 

    return;
}


/*get Cp/R as a function of T  */
/*for all species (Eq 19) */
void CKCPOR(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  cpor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpor, tc);
}


/*get H/RT as a function of T  */
/*for all species (Eq 20) */
void CKHORT(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  hort)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEnthalpy(hort, tc);
}


/*get S/R as a function of T  */
/*for all species (Eq 21) */
void CKSOR(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  sor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sor, tc);
}


/*get specific heat at constant volume as a function  */
/*of T for all species (molar units) */
void CKCVML(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  cvml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        cvml[id] *= 8.31451e+07;
    }
}


/*get specific heat at constant pressure as a  */
/*function of T for all species (molar units) */
void CKCPML(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  cpml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        cpml[id] *= 8.31451e+07;
    }
}


/*get internal energy as a function  */
/*of T for all species (molar units) */
void CKUML(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  uml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        uml[id] *= RT;
    }
}


/*get enthalpy as a function  */
/*of T for all species (molar units) */
void CKHML(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  hml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        hml[id] *= RT;
    }
}


/*get standard-state Gibbs energy as a function  */
/*of T for all species (molar units) */
void CKGML(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  gml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        gml[id] *= RT;
    }
}


/*get standard-state Helmholtz free energy as a  */
/*function of T for all species (molar units) */
void CKAML(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  aml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(aml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        aml[id] *= RT;
    }
}


/*Returns the standard-state entropies in molar units */
void CKSML(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  sml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        sml[id] *= 8.31451e+07;
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
void CKCVMS(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  cvms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 4.124383662212169e+07; /*H2 */
    cvms[1] *= 2.598381814318037e+06; /*O2 */
    cvms[2] *= 4.615239012974499e+06; /*H2O */
    cvms[3] *= 8.248767324424338e+07; /*H */
    cvms[4] *= 5.196763628636074e+06; /*O */
    cvms[5] *= 4.888768810227566e+06; /*OH */
    cvms[6] *= 2.519031701678171e+06; /*HO2 */
    cvms[7] *= 2.444384405113783e+06; /*H2O2 */
    cvms[8] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
void CKCPMS(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 4.124383662212169e+07; /*H2 */
    cpms[1] *= 2.598381814318037e+06; /*O2 */
    cpms[2] *= 4.615239012974499e+06; /*H2O */
    cpms[3] *= 8.248767324424338e+07; /*H */
    cpms[4] *= 5.196763628636074e+06; /*O */
    cpms[5] *= 4.888768810227566e+06; /*OH */
    cpms[6] *= 2.519031701678171e+06; /*HO2 */
    cpms[7] *= 2.444384405113783e+06; /*H2O2 */
    cpms[8] *= 2.968047434442088e+06; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    for (int i = 0; i < 9; i++)
    {
        ums[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void CKHMS(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hms, tc);
    for (int i = 0; i < 9; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int * restrict np, double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  hms)
{
    double tc[5], h[9];

    for (int i=0; i<(*np); i++) {
        tc[0] = 0.0;
        tc[1] = T[i];
        tc[2] = T[i]*T[i];
        tc[3] = T[i]*T[i]*T[i];
        tc[4] = T[i]*T[i]*T[i]*T[i];

        speciesEnthalpy(h, tc);

        hms[0*(*np)+i] = h[0];
        hms[1*(*np)+i] = h[1];
        hms[2*(*np)+i] = h[2];
        hms[3*(*np)+i] = h[3];
        hms[4*(*np)+i] = h[4];
        hms[5*(*np)+i] = h[5];
        hms[6*(*np)+i] = h[6];
        hms[7*(*np)+i] = h[7];
        hms[8*(*np)+i] = h[8];
    }

    for (int n=0; n<9; n++) {
        for (int i=0; i<(*np); i++) {
            hms[n*(*np)+i] *= 8.31451e+07 * T[i] * imw[n];
        }
    }
}


/*Returns gibbs in mass units (Eq 31.) */
void CKGMS(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  gms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gms, tc);
    for (int i = 0; i < 9; i++)
    {
        gms[i] *= RT*imw[i];
    }
}


/*Returns helmholtz in mass units (Eq 32.) */
void CKAMS(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  ams)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(ams, tc);
    for (int i = 0; i < 9; i++)
    {
        ams[i] *= RT*imw[i];
    }
}


/*Returns the entropies in mass units (Eq 28.) */
void CKSMS(double * restrict T, int * iwrk, double * restrict  rwrk, double * restrict  sms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 4.124383662212169e+07; /*H2 */
    sms[1] *= 2.598381814318037e+06; /*O2 */
    sms[2] *= 4.615239012974499e+06; /*H2O */
    sms[3] *= 8.248767324424338e+07; /*H */
    sms[4] *= 5.196763628636074e+06; /*O */
    sms[5] *= 4.888768810227566e+06; /*OH */
    sms[6] *= 2.519031701678171e+06; /*HO2 */
    sms[7] *= 2.444384405113783e+06; /*H2O2 */
    sms[8] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double * restrict T, double * restrict x, int * iwrk, double * restrict  rwrk, double * restrict  cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[9]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
void CKCPBS(double * restrict T, double * restrict y, int * iwrk, double * restrict  rwrk, double * restrict  cpbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[9], tresult[9]; /* temporary storage */
    cp_R(cpor, tc);
    for (int i = 0; i < 9; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 9; i++)
    {
        result += tresult[i];
    }

    *cpbs = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double * restrict T, double * restrict x, int * iwrk, double * restrict  rwrk, double * restrict  cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[9]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
void CKCVBS(double * restrict T, double * restrict y, int * iwrk, double * restrict  rwrk, double * restrict  cvbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[9]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H2 */
    result += cvor[1]*y[1]*imw[1]; /*O2 */
    result += cvor[2]*y[2]*imw[2]; /*H2O */
    result += cvor[3]*y[3]*imw[3]; /*H */
    result += cvor[4]*y[4]*imw[4]; /*O */
    result += cvor[5]*y[5]*imw[5]; /*OH */
    result += cvor[6]*y[6]*imw[6]; /*HO2 */
    result += cvor[7]*y[7]*imw[7]; /*H2O2 */
    result += cvor[8]*y[8]*imw[8]; /*N2 */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double * restrict T, double * restrict x, int * iwrk, double * restrict  rwrk, double * restrict  hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[9]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
        result += x[id]*hml[id];
    }

    *hbml = result * RT;
}


/*Returns mean enthalpy of mixture in mass units */
void CKHBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict  rwrk, double * restrict  hbms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[9], tmp[9]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    int id;
    for (id = 0; id < 9; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 9; ++id) {
        result += tmp[id];
    }

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double * restrict T, double * restrict x, int * iwrk, double * restrict  rwrk, double * restrict  ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[9]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
        result += x[id]*uml[id];
    }

    *ubml = result * RT;
}


/*get mean internal energy in mass units */
void CKUBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict  rwrk, double * restrict  ubms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double ums[9]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H2 */
    result += y[1]*ums[1]*imw[1]; /*O2 */
    result += y[2]*ums[2]*imw[2]; /*H2O */
    result += y[3]*ums[3]*imw[3]; /*H */
    result += y[4]*ums[4]*imw[4]; /*O */
    result += y[5]*ums[5]*imw[5]; /*OH */
    result += y[6]*ums[6]*imw[6]; /*HO2 */
    result += y[7]*ums[7]*imw[7]; /*H2O2 */
    result += y[8]*ums[8]*imw[8]; /*N2 */

    *ubms = result * RT;
}


/*get mixture entropy in molar units */
void CKSBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict  rwrk, double * restrict  sbml)
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

    *sbml = result * 8.31451e+07;
}


/*get mixture entropy in mass units */
void CKSBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict  rwrk, double * restrict  sbms)
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
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(18.015340*YOW); 
    x[3] = y[3]/(1.007970*YOW); 
    x[4] = y[4]/(15.999400*YOW); 
    x[5] = y[5]/(17.007370*YOW); 
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
    *sbms = result * 8.31451e+07 * YOW;
}


/*Returns mean gibbs free energy in molar units */
void CKGBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict  rwrk, double * restrict  gbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
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
void CKGBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict  rwrk, double * restrict  gbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double gort[9]; /* temporary storage */
    double x[9]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(18.015340*YOW); 
    x[3] = y[3]/(1.007970*YOW); 
    x[4] = y[4]/(15.999400*YOW); 
    x[5] = y[5]/(17.007370*YOW); 
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
void CKABML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict  rwrk, double * restrict  abml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
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
void CKABMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict  rwrk, double * restrict  abms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double aort[9]; /* temporary storage */
    double x[9]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(18.015340*YOW); 
    x[3] = y[3]/(1.007970*YOW); 
    x[4] = y[4]/(15.999400*YOW); 
    x[5] = y[5]/(17.007370*YOW); 
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
void CKWC(double * restrict  T, double * restrict  C, int * iwrk, double * restrict  rwrk, double * restrict  wdot)
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
void CKWYP(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  wdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]*imw[0]; 
    c[1] = PWORT * y[1]*imw[1]; 
    c[2] = PWORT * y[2]*imw[2]; 
    c[3] = PWORT * y[3]*imw[3]; 
    c[4] = PWORT * y[4]*imw[4]; 
    c[5] = PWORT * y[5]*imw[5]; 
    c[6] = PWORT * y[6]*imw[6]; 
    c[7] = PWORT * y[7]*imw[7]; 
    c[8] = PWORT * y[8]*imw[8]; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  wdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

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
void CKWYR(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  wdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 
    c[3] = 1e6 * (*rho) * y[3]*imw[3]; 
    c[4] = 1e6 * (*rho) * y[4]*imw[4]; 
    c[5] = 1e6 * (*rho) * y[5]*imw[5]; 
    c[6] = 1e6 * (*rho) * y[6]*imw[6]; 
    c[7] = 1e6 * (*rho) * y[7]*imw[7]; 
    c[8] = 1e6 * (*rho) * y[8]*imw[8]; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int * restrict np, double * restrict rho, double * restrict T,
	    double * restrict y, int * restrict iwrk, double * restrict rwrk,
	    double * restrict wdot)
{
    double c[9*(*np)]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<9; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<9*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  wdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
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
void CKQC(double * restrict  T, double * restrict  C, int * iwrk, double * restrict  rwrk, double * restrict  qdot)
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

    for (id = 0; id < 21; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  q_f, double * restrict  q_r)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]*imw[0]; 
    c[1] = PWORT * y[1]*imw[1]; 
    c[2] = PWORT * y[2]*imw[2]; 
    c[3] = PWORT * y[3]*imw[3]; 
    c[4] = PWORT * y[4]*imw[4]; 
    c[5] = PWORT * y[5]*imw[5]; 
    c[6] = PWORT * y[6]*imw[6]; 
    c[7] = PWORT * y[7]*imw[7]; 
    c[8] = PWORT * y[8]*imw[8]; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 
    c[3] = 1e6 * (*rho) * y[3]*imw[3]; 
    c[4] = 1e6 * (*rho) * y[4]*imw[4]; 
    c[5] = 1e6 * (*rho) * y[5]*imw[5]; 
    c[6] = 1e6 * (*rho) * y[6]*imw[6]; 
    c[7] = 1e6 * (*rho) * y[7]*imw[7]; 
    c[8] = 1e6 * (*rho) * y[8]*imw[8]; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
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
    for (id = 0; id < 21; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the stoichiometric coefficients */
/*of the reaction mechanism. (Eq 50) */
void CKNU(int * kdim, int * iwrk, double * restrict  rwrk, int * nuki)
{
    int id; /*loop counter */
    int kd = (*kdim); 
    /*Zero nuki */
    for (id = 0; id < 9 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: H + O2 <=> O + OH */
    nuki[ 3 * kd + 0 ] += -1 ;
    nuki[ 1 * kd + 0 ] += -1 ;
    nuki[ 4 * kd + 0 ] += +1 ;
    nuki[ 5 * kd + 0 ] += +1 ;

    /*reaction 2: O + H2 <=> H + OH */
    nuki[ 4 * kd + 1 ] += -1 ;
    nuki[ 0 * kd + 1 ] += -1 ;
    nuki[ 3 * kd + 1 ] += +1 ;
    nuki[ 5 * kd + 1 ] += +1 ;

    /*reaction 3: H2 + OH <=> H2O + H */
    nuki[ 0 * kd + 2 ] += -1 ;
    nuki[ 5 * kd + 2 ] += -1 ;
    nuki[ 2 * kd + 2 ] += +1 ;
    nuki[ 3 * kd + 2 ] += +1 ;

    /*reaction 4: O + H2O <=> OH + OH */
    nuki[ 4 * kd + 3 ] += -1 ;
    nuki[ 2 * kd + 3 ] += -1 ;
    nuki[ 5 * kd + 3 ] += +1 ;
    nuki[ 5 * kd + 3 ] += +1 ;

    /*reaction 5: H2 + M <=> H + H + M */
    nuki[ 0 * kd + 4 ] += -1 ;
    nuki[ 3 * kd + 4 ] += +1 ;
    nuki[ 3 * kd + 4 ] += +1 ;

    /*reaction 6: O + O + M <=> O2 + M */
    nuki[ 4 * kd + 5 ] += -1 ;
    nuki[ 4 * kd + 5 ] += -1 ;
    nuki[ 1 * kd + 5 ] += +1 ;

    /*reaction 7: O + H + M <=> OH + M */
    nuki[ 4 * kd + 6 ] += -1 ;
    nuki[ 3 * kd + 6 ] += -1 ;
    nuki[ 5 * kd + 6 ] += +1 ;

    /*reaction 8: H + OH + M <=> H2O + M */
    nuki[ 3 * kd + 7 ] += -1 ;
    nuki[ 5 * kd + 7 ] += -1 ;
    nuki[ 2 * kd + 7 ] += +1 ;

    /*reaction 9: H + O2 (+M) <=> HO2 (+M) */
    nuki[ 3 * kd + 8 ] += -1 ;
    nuki[ 1 * kd + 8 ] += -1 ;
    nuki[ 6 * kd + 8 ] += +1 ;

    /*reaction 10: HO2 + H <=> H2 + O2 */
    nuki[ 6 * kd + 9 ] += -1 ;
    nuki[ 3 * kd + 9 ] += -1 ;
    nuki[ 0 * kd + 9 ] += +1 ;
    nuki[ 1 * kd + 9 ] += +1 ;

    /*reaction 11: HO2 + H <=> OH + OH */
    nuki[ 6 * kd + 10 ] += -1 ;
    nuki[ 3 * kd + 10 ] += -1 ;
    nuki[ 5 * kd + 10 ] += +1 ;
    nuki[ 5 * kd + 10 ] += +1 ;

    /*reaction 12: HO2 + O <=> O2 + OH */
    nuki[ 6 * kd + 11 ] += -1 ;
    nuki[ 4 * kd + 11 ] += -1 ;
    nuki[ 1 * kd + 11 ] += +1 ;
    nuki[ 5 * kd + 11 ] += +1 ;

    /*reaction 13: HO2 + OH <=> H2O + O2 */
    nuki[ 6 * kd + 12 ] += -1 ;
    nuki[ 5 * kd + 12 ] += -1 ;
    nuki[ 2 * kd + 12 ] += +1 ;
    nuki[ 1 * kd + 12 ] += +1 ;

    /*reaction 14: HO2 + HO2 <=> H2O2 + O2 */
    nuki[ 6 * kd + 13 ] += -1 ;
    nuki[ 6 * kd + 13 ] += -1 ;
    nuki[ 7 * kd + 13 ] += +1 ;
    nuki[ 1 * kd + 13 ] += +1 ;

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    nuki[ 6 * kd + 14 ] += -1 ;
    nuki[ 6 * kd + 14 ] += -1 ;
    nuki[ 7 * kd + 14 ] += +1 ;
    nuki[ 1 * kd + 14 ] += +1 ;

    /*reaction 16: H2O2 (+M) <=> OH + OH (+M) */
    nuki[ 7 * kd + 15 ] += -1 ;
    nuki[ 5 * kd + 15 ] += +1 ;
    nuki[ 5 * kd + 15 ] += +1 ;

    /*reaction 17: H2O2 + H <=> H2O + OH */
    nuki[ 7 * kd + 16 ] += -1 ;
    nuki[ 3 * kd + 16 ] += -1 ;
    nuki[ 2 * kd + 16 ] += +1 ;
    nuki[ 5 * kd + 16 ] += +1 ;

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    nuki[ 7 * kd + 17 ] += -1 ;
    nuki[ 3 * kd + 17 ] += -1 ;
    nuki[ 6 * kd + 17 ] += +1 ;
    nuki[ 0 * kd + 17 ] += +1 ;

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    nuki[ 7 * kd + 18 ] += -1 ;
    nuki[ 4 * kd + 18 ] += -1 ;
    nuki[ 5 * kd + 18 ] += +1 ;
    nuki[ 6 * kd + 18 ] += +1 ;

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    nuki[ 7 * kd + 19 ] += -1 ;
    nuki[ 5 * kd + 19 ] += -1 ;
    nuki[ 6 * kd + 19 ] += +1 ;
    nuki[ 2 * kd + 19 ] += +1 ;

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    nuki[ 7 * kd + 20 ] += -1 ;
    nuki[ 5 * kd + 20 ] += -1 ;
    nuki[ 6 * kd + 20 ] += +1 ;
    nuki[ 2 * kd + 20 ] += +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * restrict  rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < 3 * 9; ++ id) {
         ncf[id] = 0; 
    }

    /*H2 */
    ncf[ 0 * kd + 0 ] = 2; /*H */

    /*O2 */
    ncf[ 1 * kd + 1 ] = 2; /*O */

    /*H2O */
    ncf[ 2 * kd + 0 ] = 2; /*H */
    ncf[ 2 * kd + 1 ] = 1; /*O */

    /*H */
    ncf[ 3 * kd + 0 ] = 1; /*H */

    /*O */
    ncf[ 4 * kd + 1 ] = 1; /*O */

    /*OH */
    ncf[ 5 * kd + 1 ] = 1; /*O */
    ncf[ 5 * kd + 0 ] = 1; /*H */

    /*HO2 */
    ncf[ 6 * kd + 0 ] = 1; /*H */
    ncf[ 6 * kd + 1 ] = 2; /*O */

    /*H2O2 */
    ncf[ 7 * kd + 0 ] = 2; /*H */
    ncf[ 7 * kd + 1 ] = 2; /*O */

    /*N2 */
    ncf[ 8 * kd + 2 ] = 2; /*N */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE(int * iwrk, double * restrict  rwrk, double * restrict  a, double * restrict  b, double * restrict  e)
{

    /*reaction 1: H + O2 <=> O + OH */
    a[0] = 3.547e+15;
    b[0] = -0.406;
    e[0] = 16599;

    /*reaction 2: O + H2 <=> H + OH */
    a[1] = 50800;
    b[1] = 2.67;
    e[1] = 6290;

    /*reaction 3: H2 + OH <=> H2O + H */
    a[2] = 2.16e+08;
    b[2] = 1.51;
    e[2] = 3430;

    /*reaction 4: O + H2O <=> OH + OH */
    a[3] = 2.97e+06;
    b[3] = 2.02;
    e[3] = 13400;

    /*reaction 5: H2 + M <=> H + H + M */
    a[4] = 4.577e+19;
    b[4] = -1.4;
    e[4] = 104380;

    /*reaction 6: O + O + M <=> O2 + M */
    a[5] = 6.165e+15;
    b[5] = -0.5;
    e[5] = 0;

    /*reaction 7: O + H + M <=> OH + M */
    a[6] = 4.714e+18;
    b[6] = -1;
    e[6] = 0;

    /*reaction 8: H + OH + M <=> H2O + M */
    a[7] = 3.8e+22;
    b[7] = -2;
    e[7] = 0;

    /*reaction 9: H + O2 (+M) <=> HO2 (+M) */
    a[8] = 1.475e+12;
    b[8] = 0.6;
    e[8] = 0;

    /*reaction 10: HO2 + H <=> H2 + O2 */
    a[9] = 1.66e+13;
    b[9] = 0;
    e[9] = 823;

    /*reaction 11: HO2 + H <=> OH + OH */
    a[10] = 7.079e+13;
    b[10] = 0;
    e[10] = 295;

    /*reaction 12: HO2 + O <=> O2 + OH */
    a[11] = 3.25e+13;
    b[11] = 0;
    e[11] = 0;

    /*reaction 13: HO2 + OH <=> H2O + O2 */
    a[12] = 2.89e+13;
    b[12] = 0;
    e[12] = -497;

    /*reaction 14: HO2 + HO2 <=> H2O2 + O2 */
    a[13] = 4.2e+14;
    b[13] = 0;
    e[13] = 11982;

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    a[14] = 1.3e+11;
    b[14] = 0;
    e[14] = -1629.3;

    /*reaction 16: H2O2 (+M) <=> OH + OH (+M) */
    a[15] = 2.951e+14;
    b[15] = 0;
    e[15] = 48430;

    /*reaction 17: H2O2 + H <=> H2O + OH */
    a[16] = 2.41e+13;
    b[16] = 0;
    e[16] = 3970;

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    a[17] = 4.82e+13;
    b[17] = 0;
    e[17] = 7950;

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    a[18] = 9.55e+06;
    b[18] = 2;
    e[18] = 3970;

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    a[19] = 1e+12;
    b[19] = 0;
    e[19] = 0;

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    a[20] = 5.8e+14;
    b[20] = 0;
    e[20] = 9557;

    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double * restrict  T, double * restrict  C, int * iwrk, double * restrict  rwrk, double * restrict  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> O + OH */
    /*eqcon[0] *= 1;  */

    /*reaction 2: O + H2 <=> H + OH */
    /*eqcon[1] *= 1;  */

    /*reaction 3: H2 + OH <=> H2O + H */
    /*eqcon[2] *= 1;  */

    /*reaction 4: O + H2O <=> OH + OH */
    /*eqcon[3] *= 1;  */

    /*reaction 5: H2 + M <=> H + H + M */
    eqcon[4] *= 1e-06; 

    /*reaction 6: O + O + M <=> O2 + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + H + M <=> OH + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + OH + M <=> H2O + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H + O2 (+M) <=> HO2 (+M) */
    eqcon[8] *= 1e+06; 

    /*reaction 10: HO2 + H <=> H2 + O2 */
    /*eqcon[9] *= 1;  */

    /*reaction 11: HO2 + H <=> OH + OH */
    /*eqcon[10] *= 1;  */

    /*reaction 12: HO2 + O <=> O2 + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + OH <=> H2O + O2 */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[15] *= 1e-06; 

    /*reaction 17: H2O2 + H <=> H2O + OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[20] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * restrict  P, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> O + OH */
    /*eqcon[0] *= 1;  */

    /*reaction 2: O + H2 <=> H + OH */
    /*eqcon[1] *= 1;  */

    /*reaction 3: H2 + OH <=> H2O + H */
    /*eqcon[2] *= 1;  */

    /*reaction 4: O + H2O <=> OH + OH */
    /*eqcon[3] *= 1;  */

    /*reaction 5: H2 + M <=> H + H + M */
    eqcon[4] *= 1e-06; 

    /*reaction 6: O + O + M <=> O2 + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + H + M <=> OH + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + OH + M <=> H2O + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H + O2 (+M) <=> HO2 (+M) */
    eqcon[8] *= 1e+06; 

    /*reaction 10: HO2 + H <=> H2 + O2 */
    /*eqcon[9] *= 1;  */

    /*reaction 11: HO2 + H <=> OH + OH */
    /*eqcon[10] *= 1;  */

    /*reaction 12: HO2 + O <=> O2 + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + OH <=> H2O + O2 */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[15] *= 1e-06; 

    /*reaction 17: H2O2 + H <=> H2O + OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[20] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * restrict  P, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> O + OH */
    /*eqcon[0] *= 1;  */

    /*reaction 2: O + H2 <=> H + OH */
    /*eqcon[1] *= 1;  */

    /*reaction 3: H2 + OH <=> H2O + H */
    /*eqcon[2] *= 1;  */

    /*reaction 4: O + H2O <=> OH + OH */
    /*eqcon[3] *= 1;  */

    /*reaction 5: H2 + M <=> H + H + M */
    eqcon[4] *= 1e-06; 

    /*reaction 6: O + O + M <=> O2 + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + H + M <=> OH + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + OH + M <=> H2O + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H + O2 (+M) <=> HO2 (+M) */
    eqcon[8] *= 1e+06; 

    /*reaction 10: HO2 + H <=> H2 + O2 */
    /*eqcon[9] *= 1;  */

    /*reaction 11: HO2 + H <=> OH + OH */
    /*eqcon[10] *= 1;  */

    /*reaction 12: HO2 + O <=> O2 + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + OH <=> H2O + O2 */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[15] *= 1e-06; 

    /*reaction 17: H2O2 + H <=> H2O + OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[20] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * restrict  rho, double * restrict  T, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> O + OH */
    /*eqcon[0] *= 1;  */

    /*reaction 2: O + H2 <=> H + OH */
    /*eqcon[1] *= 1;  */

    /*reaction 3: H2 + OH <=> H2O + H */
    /*eqcon[2] *= 1;  */

    /*reaction 4: O + H2O <=> OH + OH */
    /*eqcon[3] *= 1;  */

    /*reaction 5: H2 + M <=> H + H + M */
    eqcon[4] *= 1e-06; 

    /*reaction 6: O + O + M <=> O2 + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + H + M <=> OH + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + OH + M <=> H2O + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H + O2 (+M) <=> HO2 (+M) */
    eqcon[8] *= 1e+06; 

    /*reaction 10: HO2 + H <=> H2 + O2 */
    /*eqcon[9] *= 1;  */

    /*reaction 11: HO2 + H <=> OH + OH */
    /*eqcon[10] *= 1;  */

    /*reaction 12: HO2 + O <=> O2 + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + OH <=> H2O + O2 */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[15] *= 1e-06; 

    /*reaction 17: H2O2 + H <=> H2O + OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[20] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * restrict  rho, double * restrict  T, double * restrict  x, int * iwrk, double * restrict  rwrk, double * restrict  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 <=> O + OH */
    /*eqcon[0] *= 1;  */

    /*reaction 2: O + H2 <=> H + OH */
    /*eqcon[1] *= 1;  */

    /*reaction 3: H2 + OH <=> H2O + H */
    /*eqcon[2] *= 1;  */

    /*reaction 4: O + H2O <=> OH + OH */
    /*eqcon[3] *= 1;  */

    /*reaction 5: H2 + M <=> H + H + M */
    eqcon[4] *= 1e-06; 

    /*reaction 6: O + O + M <=> O2 + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + H + M <=> OH + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + OH + M <=> H2O + M */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H + O2 (+M) <=> HO2 (+M) */
    eqcon[8] *= 1e+06; 

    /*reaction 10: HO2 + H <=> H2 + O2 */
    /*eqcon[9] *= 1;  */

    /*reaction 11: HO2 + H <=> OH + OH */
    /*eqcon[10] *= 1;  */

    /*reaction 12: HO2 + O <=> O2 + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + OH <=> H2O + O2 */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[15] *= 1e-06; 

    /*reaction 17: H2O2 + H <=> H2O + OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[20] *= 1;  */
}

static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[21];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[21];
#ifdef _OPENMP
#pragma omp threadprivate(Kc_save)
#endif


/*compute the production rate for each species */
void productionRate(double * restrict  wdot, double * restrict  sc, double T)
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

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */

    double invT = 1.0 / tc[1];

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;
    double refCinv = 1 / refC;

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

    struct ReactionData* Ri;
    if (T != T_save)
    {
        T_save = T;

        for (int i=0; i<21; ++i) {
            Ri = &(R[i]);
            k_f_save[i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[0] - Ri->activation_units * Ri->fwd_Ea * invT);
        };
        Kc_save[0] = exp((g_RT[3] + g_RT[1]) - (g_RT[4] + g_RT[5]));
        Kc_save[1] = exp((g_RT[4] + g_RT[0]) - (g_RT[3] + g_RT[5]));
        Kc_save[2] = exp((g_RT[0] + g_RT[5]) - (g_RT[2] + g_RT[3]));
        Kc_save[3] = exp((g_RT[4] + g_RT[2]) - (g_RT[5] + g_RT[5]));
        Kc_save[4] = refC * exp((g_RT[0]) - (g_RT[3] + g_RT[3]));
        Kc_save[5] = 1.0 / (refC) * exp((g_RT[4] + g_RT[4]) - (g_RT[1]));
        Kc_save[6] = 1.0 / (refC) * exp((g_RT[4] + g_RT[3]) - (g_RT[5]));
        Kc_save[7] = 1.0 / (refC) * exp((g_RT[3] + g_RT[5]) - (g_RT[2]));
        Kc_save[8] = 1.0 / (refC) * exp((g_RT[3] + g_RT[1]) - (g_RT[6]));
        Kc_save[9] = exp((g_RT[6] + g_RT[3]) - (g_RT[0] + g_RT[1]));
        Kc_save[10] = exp((g_RT[6] + g_RT[3]) - (g_RT[5] + g_RT[5]));
        Kc_save[11] = exp((g_RT[6] + g_RT[4]) - (g_RT[1] + g_RT[5]));
        Kc_save[12] = exp((g_RT[6] + g_RT[5]) - (g_RT[2] + g_RT[1]));
        Kc_save[13] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[1]));
        Kc_save[14] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[1]));
        Kc_save[15] = refC * exp((g_RT[7]) - (g_RT[5] + g_RT[5]));
        Kc_save[16] = exp((g_RT[7] + g_RT[3]) - (g_RT[2] + g_RT[5]));
        Kc_save[17] = exp((g_RT[7] + g_RT[3]) - (g_RT[6] + g_RT[0]));
        Kc_save[18] = exp((g_RT[7] + g_RT[4]) - (g_RT[5] + g_RT[6]));
        Kc_save[19] = exp((g_RT[7] + g_RT[5]) - (g_RT[6] + g_RT[2]));
        Kc_save[20] = exp((g_RT[7] + g_RT[5]) - (g_RT[6] + g_RT[2]));
    }

    /*reaction 1: H + O2 <=> O + OH */
    phi_f = sc[3]*sc[1];
    k_f = k_f_save[0];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[5];
    Kc = Kc_save[0];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[4] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 2: O + H2 <=> H + OH */
    phi_f = sc[4]*sc[0];
    k_f = k_f_save[1];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = Kc_save[1];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[0] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 3: H2 + OH <=> H2O + H */
    phi_f = sc[0]*sc[5];
    k_f = k_f_save[2];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[3];
    Kc = Kc_save[2];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 4: O + H2O <=> OH + OH */
    phi_f = sc[4]*sc[2];
    k_f = k_f_save[3];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[5];
    Kc = Kc_save[3];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[2] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 5: H2 + M <=> H + H + M */
    phi_f = sc[0];
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    k_f = alpha * k_f_save[4];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[3];
    Kc = Kc_save[4];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[0] -= 1 * qdot;
    wdot[3] += 1 * qdot;
    wdot[3] += 1 * qdot;

    /*reaction 6: O + O + M <=> O2 + M */
    phi_f = sc[4]*sc[4];
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    k_f = alpha * k_f_save[5];
    q_f = phi_f * k_f;
    phi_r = sc[1];
    Kc = Kc_save[5];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 7: O + H + M <=> OH + M */
    phi_f = sc[4]*sc[3];
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    k_f = alpha * k_f_save[6];
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = Kc_save[6];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[4] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 8: H + OH + M <=> H2O + M */
    phi_f = sc[3]*sc[5];
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    k_f = alpha * k_f_save[7];
    q_f = phi_f * k_f;
    phi_r = sc[2];
    Kc = Kc_save[7];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 9: H + O2 (+M) <=> HO2 (+M) */
    phi_f = sc[3]*sc[1];
    alpha = mixture + sc[0] + 10*sc[2] + -0.22*sc[1];
    k_f = k_f_save[8];
    Ri = &(R[8]);
    redP = alpha / k_f * Ri->phase_units * Ri->low_A * exp(Ri->low_beta * tc[0] - Ri->activation_units * Ri->low_Ea *invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (fabs(Ri->troe_Tsss) > 1.e-100 ? (1-Ri->troe_a)*exp(-T/Ri->troe_Tsss) : 0) 
        + (fabs(Ri->troe_Ts) > 1.e-100 ? Ri->troe_a * exp(-T/Ri->troe_Ts) : 0) 
        + (Ri->troe_len == 4 ? exp(-Ri->troe_Tss * invT) : 0) );
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = Kc_save[8];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[3] -= 1 * qdot;
    wdot[1] -= 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 10: HO2 + H <=> H2 + O2 */
    phi_f = sc[6]*sc[3];
    k_f = k_f_save[9];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[1];
    Kc = Kc_save[9];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[0] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 11: HO2 + H <=> OH + OH */
    phi_f = sc[6]*sc[3];
    k_f = k_f_save[10];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[5];
    Kc = Kc_save[10];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 12: HO2 + O <=> O2 + OH */
    phi_f = sc[6]*sc[4];
    k_f = k_f_save[11];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    Kc = Kc_save[11];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[1] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 13: HO2 + OH <=> H2O + O2 */
    phi_f = sc[6]*sc[5];
    k_f = k_f_save[12];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[1];
    Kc = Kc_save[12];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 14: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[6]*sc[6];
    k_f = k_f_save[13];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[1];
    Kc = Kc_save[13];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[6]*sc[6];
    k_f = k_f_save[14];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[1];
    Kc = Kc_save[14];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[6] -= 1 * qdot;
    wdot[6] -= 1 * qdot;
    wdot[7] += 1 * qdot;
    wdot[1] += 1 * qdot;

    /*reaction 16: H2O2 (+M) <=> OH + OH (+M) */
    phi_f = sc[7];
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    k_f = k_f_save[15];
    Ri = &(R[15]);
    redP = alpha / k_f * Ri->phase_units * Ri->low_A * exp(Ri->low_beta * tc[0] - Ri->activation_units * Ri->low_Ea *invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (fabs(Ri->troe_Tsss) > 1.e-100 ? (1-Ri->troe_a)*exp(-T/Ri->troe_Tsss) : 0) 
        + (fabs(Ri->troe_Ts) > 1.e-100 ? Ri->troe_a * exp(-T/Ri->troe_Ts) : 0) 
        + (Ri->troe_len == 4 ? exp(-Ri->troe_Tss * invT) : 0) );
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[5];
    Kc = Kc_save[15];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 17: H2O2 + H <=> H2O + OH */
    phi_f = sc[7]*sc[3];
    k_f = k_f_save[16];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = Kc_save[16];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[2] += 1 * qdot;
    wdot[5] += 1 * qdot;

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    phi_f = sc[7]*sc[3];
    k_f = k_f_save[17];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = Kc_save[17];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[3] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[0] += 1 * qdot;

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    phi_f = sc[7]*sc[4];
    k_f = k_f_save[18];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[6];
    Kc = Kc_save[18];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[4] -= 1 * qdot;
    wdot[5] += 1 * qdot;
    wdot[6] += 1 * qdot;

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    phi_f = sc[7]*sc[5];
    k_f = k_f_save[19];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[2];
    Kc = Kc_save[19];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[2] += 1 * qdot;

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    phi_f = sc[7]*sc[5];
    k_f = k_f_save[20];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[2];
    Kc = Kc_save[20];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot = q_f - q_r;
    wdot[7] -= 1 * qdot;
    wdot[5] -= 1 * qdot;
    wdot[6] += 1 * qdot;
    wdot[2] += 1 * qdot;

    return;
}


/*compute the production rate for each species */
void vproductionRate(int npt, double * restrict wdot, double * restrict sc, double * restrict T)
{
    double k_f_s[21][npt], Kc_s[21][npt], mixture[npt], g_RT[9*npt];
    double tc[5*npt], invT[npt];
    struct ReactionData* Ri;

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        tc[0*npt+i] = log(T[i]);
        tc[1*npt+i] = T[i];
        tc[2*npt+i] = T[i]*T[i];
        tc[3*npt+i] = T[i]*T[i]*T[i];
        tc[4*npt+i] = T[i]*T[i]*T[i]*T[i];
        invT[i] = 1.0 / T[i];
    }

#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        Ri = &(R[i]);
        k_f_s[0][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[1][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[2][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[3][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[4][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[5][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[6][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[7][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[8][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[9][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[10][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[11][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[12][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[13][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[14][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[15][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[16][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[17][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[18][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[19][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        k_f_s[20][i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
    }

    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        double tg[5], g[9];
        tg[0] = tc[0*npt+i];
        tg[1] = tc[1*npt+i];
        tg[2] = tc[2*npt+i];
        tg[3] = tc[3*npt+i];
        tg[4] = tc[4*npt+i];

        gibbs(g, tg);

        g_RT[0*npt+i] = g[0];
        g_RT[1*npt+i] = g[1];
        g_RT[2*npt+i] = g[2];
        g_RT[3*npt+i] = g[3];
        g_RT[4*npt+i] = g[4];
        g_RT[5*npt+i] = g[5];
        g_RT[6*npt+i] = g[6];
        g_RT[7*npt+i] = g[7];
        g_RT[8*npt+i] = g[8];
    }

#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        double refC = 101325. / 8.31451 / T[i];

        Kc_s[0][i] = exp((g_RT[3*npt+i] + g_RT[1*npt+i]) - (g_RT[4*npt+i] + g_RT[5*npt+i]));
        Kc_s[1][i] = exp((g_RT[4*npt+i] + g_RT[0*npt+i]) - (g_RT[3*npt+i] + g_RT[5*npt+i]));
        Kc_s[2][i] = exp((g_RT[0*npt+i] + g_RT[5*npt+i]) - (g_RT[2*npt+i] + g_RT[3*npt+i]));
        Kc_s[3][i] = exp((g_RT[4*npt+i] + g_RT[2*npt+i]) - (g_RT[5*npt+i] + g_RT[5*npt+i]));
        Kc_s[4][i] = refC * exp((g_RT[0*npt+i]) - (g_RT[3*npt+i] + g_RT[3*npt+i]));
        Kc_s[5][i] = 1.0 / (refC) * exp((g_RT[4*npt+i] + g_RT[4*npt+i]) - (g_RT[1*npt+i]));
        Kc_s[6][i] = 1.0 / (refC) * exp((g_RT[4*npt+i] + g_RT[3*npt+i]) - (g_RT[5*npt+i]));
        Kc_s[7][i] = 1.0 / (refC) * exp((g_RT[3*npt+i] + g_RT[5*npt+i]) - (g_RT[2*npt+i]));
        Kc_s[8][i] = 1.0 / (refC) * exp((g_RT[3*npt+i] + g_RT[1*npt+i]) - (g_RT[6*npt+i]));
        Kc_s[9][i] = exp((g_RT[6*npt+i] + g_RT[3*npt+i]) - (g_RT[0*npt+i] + g_RT[1*npt+i]));
        Kc_s[10][i] = exp((g_RT[6*npt+i] + g_RT[3*npt+i]) - (g_RT[5*npt+i] + g_RT[5*npt+i]));
        Kc_s[11][i] = exp((g_RT[6*npt+i] + g_RT[4*npt+i]) - (g_RT[1*npt+i] + g_RT[5*npt+i]));
        Kc_s[12][i] = exp((g_RT[6*npt+i] + g_RT[5*npt+i]) - (g_RT[2*npt+i] + g_RT[1*npt+i]));
        Kc_s[13][i] = exp((g_RT[6*npt+i] + g_RT[6*npt+i]) - (g_RT[7*npt+i] + g_RT[1*npt+i]));
        Kc_s[14][i] = exp((g_RT[6*npt+i] + g_RT[6*npt+i]) - (g_RT[7*npt+i] + g_RT[1*npt+i]));
        Kc_s[15][i] = refC * exp((g_RT[7*npt+i]) - (g_RT[5*npt+i] + g_RT[5*npt+i]));
        Kc_s[16][i] = exp((g_RT[7*npt+i] + g_RT[3*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i]));
        Kc_s[17][i] = exp((g_RT[7*npt+i] + g_RT[3*npt+i]) - (g_RT[6*npt+i] + g_RT[0*npt+i]));
        Kc_s[18][i] = exp((g_RT[7*npt+i] + g_RT[4*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[19][i] = exp((g_RT[7*npt+i] + g_RT[5*npt+i]) - (g_RT[6*npt+i] + g_RT[2*npt+i]));
        Kc_s[20][i] = exp((g_RT[7*npt+i] + g_RT[5*npt+i]) - (g_RT[6*npt+i] + g_RT[2*npt+i]));
    }

    for (int i=0; i<npt; i++) {
        mixture[i] = 0.0;
    }

    for (int n=0; n<9; n++) {
        for (int i=0; i<npt; i++) {
            mixture[i] += sc[n*npt+i];
            wdot[n*npt+i] = 0.0;
        }
    }

#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        double alpha, redP, F, logPred, logFcent;
        double troe_c, troe_n, troe, F_troe;
        double X, F_src;

        /*reaction 1: H + O2 <=> O + OH */
        phi_f = sc[3*npt+i]*sc[1*npt+i];
        k_f = k_f_s[0][i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[5*npt+i];
        Kc = Kc_s[0][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= 1 * qdot;
        wdot[1*npt+i] -= 1 * qdot;
        wdot[4*npt+i] += 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;

        /*reaction 2: O + H2 <=> H + OH */
        phi_f = sc[4*npt+i]*sc[0*npt+i];
        k_f = k_f_s[1][i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[5*npt+i];
        Kc = Kc_s[1][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= 1 * qdot;
        wdot[0*npt+i] -= 1 * qdot;
        wdot[3*npt+i] += 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;

        /*reaction 3: H2 + OH <=> H2O + H */
        phi_f = sc[0*npt+i]*sc[5*npt+i];
        k_f = k_f_s[2][i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[3*npt+i];
        Kc = Kc_s[2][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= 1 * qdot;
        wdot[5*npt+i] -= 1 * qdot;
        wdot[2*npt+i] += 1 * qdot;
        wdot[3*npt+i] += 1 * qdot;

        /*reaction 4: O + H2O <=> OH + OH */
        phi_f = sc[4*npt+i]*sc[2*npt+i];
        k_f = k_f_s[3][i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[5*npt+i];
        Kc = Kc_s[3][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= 1 * qdot;
        wdot[2*npt+i] -= 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;

        /*reaction 5: H2 + M <=> H + H + M */
        phi_f = sc[0*npt+i];
        alpha = mixture[i] + 1.5*sc[0*npt+i] + 11*sc[2*npt+i];
        k_f = alpha * k_f_s[4][i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[3*npt+i];
        Kc = Kc_s[4][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= 1 * qdot;
        wdot[3*npt+i] += 1 * qdot;
        wdot[3*npt+i] += 1 * qdot;

        /*reaction 6: O + O + M <=> O2 + M */
        phi_f = sc[4*npt+i]*sc[4*npt+i];
        alpha = mixture[i] + 1.5*sc[0*npt+i] + 11*sc[2*npt+i];
        k_f = alpha * k_f_s[5][i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i];
        Kc = Kc_s[5][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= 1 * qdot;
        wdot[4*npt+i] -= 1 * qdot;
        wdot[1*npt+i] += 1 * qdot;

        /*reaction 7: O + H + M <=> OH + M */
        phi_f = sc[4*npt+i]*sc[3*npt+i];
        alpha = mixture[i] + 1.5*sc[0*npt+i] + 11*sc[2*npt+i];
        k_f = alpha * k_f_s[6][i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i];
        Kc = Kc_s[6][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= 1 * qdot;
        wdot[3*npt+i] -= 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;

        /*reaction 8: H + OH + M <=> H2O + M */
        phi_f = sc[3*npt+i]*sc[5*npt+i];
        alpha = mixture[i] + 1.5*sc[0*npt+i] + 11*sc[2*npt+i];
        k_f = alpha * k_f_s[7][i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i];
        Kc = Kc_s[7][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= 1 * qdot;
        wdot[5*npt+i] -= 1 * qdot;
        wdot[2*npt+i] += 1 * qdot;

        /*reaction 9: H + O2 (+M) <=> HO2 (+M) */
        phi_f = sc[3*npt+i]*sc[1*npt+i];
        alpha = mixture[i] + sc[0*npt+i] + 10*sc[2*npt+i] + -0.22*sc[1*npt+i];
        Ri = &(R[8]);;
        redP = alpha / k_f * Ri->phase_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(Ri->troe_Tsss) > 1.e-100 ? (1-Ri->troe_a)*exp(-T[i]/Ri->troe_Tsss) : 0) 
            + (fabs(Ri->troe_Ts) > 1.e-100 ? Ri->troe_a * exp(-T[i]/Ri->troe_Ts) : 0) 
            + (Ri->troe_len == 4 ? exp(-Ri->troe_Tss * invT[i]) : 0) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10, logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i];
        Kc = Kc_s[8][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= 1 * qdot;
        wdot[1*npt+i] -= 1 * qdot;
        wdot[6*npt+i] += 1 * qdot;

        /*reaction 10: HO2 + H <=> H2 + O2 */
        phi_f = sc[6*npt+i]*sc[3*npt+i];
        k_f = k_f_s[9][i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[1*npt+i];
        Kc = Kc_s[9][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= 1 * qdot;
        wdot[3*npt+i] -= 1 * qdot;
        wdot[0*npt+i] += 1 * qdot;
        wdot[1*npt+i] += 1 * qdot;

        /*reaction 11: HO2 + H <=> OH + OH */
        phi_f = sc[6*npt+i]*sc[3*npt+i];
        k_f = k_f_s[10][i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[5*npt+i];
        Kc = Kc_s[10][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= 1 * qdot;
        wdot[3*npt+i] -= 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;

        /*reaction 12: HO2 + O <=> O2 + OH */
        phi_f = sc[6*npt+i]*sc[4*npt+i];
        k_f = k_f_s[11][i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[5*npt+i];
        Kc = Kc_s[11][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= 1 * qdot;
        wdot[4*npt+i] -= 1 * qdot;
        wdot[1*npt+i] += 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;

        /*reaction 13: HO2 + OH <=> H2O + O2 */
        phi_f = sc[6*npt+i]*sc[5*npt+i];
        k_f = k_f_s[12][i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[1*npt+i];
        Kc = Kc_s[12][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= 1 * qdot;
        wdot[5*npt+i] -= 1 * qdot;
        wdot[2*npt+i] += 1 * qdot;
        wdot[1*npt+i] += 1 * qdot;

        /*reaction 14: HO2 + HO2 <=> H2O2 + O2 */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[13][i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[1*npt+i];
        Kc = Kc_s[13][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= 1 * qdot;
        wdot[6*npt+i] -= 1 * qdot;
        wdot[7*npt+i] += 1 * qdot;
        wdot[1*npt+i] += 1 * qdot;

        /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[14][i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[1*npt+i];
        Kc = Kc_s[14][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= 1 * qdot;
        wdot[6*npt+i] -= 1 * qdot;
        wdot[7*npt+i] += 1 * qdot;
        wdot[1*npt+i] += 1 * qdot;

        /*reaction 16: H2O2 (+M) <=> OH + OH (+M) */
        phi_f = sc[7*npt+i];
        alpha = mixture[i] + 1.5*sc[0*npt+i] + 11*sc[2*npt+i];
        Ri = &(R[15]);;
        redP = alpha / k_f * Ri->phase_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[i] - Ri->activation_units * Ri->fwd_Ea * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(Ri->troe_Tsss) > 1.e-100 ? (1-Ri->troe_a)*exp(-T[i]/Ri->troe_Tsss) : 0) 
            + (fabs(Ri->troe_Ts) > 1.e-100 ? Ri->troe_a * exp(-T[i]/Ri->troe_Ts) : 0) 
            + (Ri->troe_len == 4 ? exp(-Ri->troe_Tss * invT[i]) : 0) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10, logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[5*npt+i];
        Kc = Kc_s[15][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;

        /*reaction 17: H2O2 + H <=> H2O + OH */
        phi_f = sc[7*npt+i]*sc[3*npt+i];
        k_f = k_f_s[16][i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[5*npt+i];
        Kc = Kc_s[16][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= 1 * qdot;
        wdot[3*npt+i] -= 1 * qdot;
        wdot[2*npt+i] += 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;

        /*reaction 18: H2O2 + H <=> HO2 + H2 */
        phi_f = sc[7*npt+i]*sc[3*npt+i];
        k_f = k_f_s[17][i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[0*npt+i];
        Kc = Kc_s[17][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= 1 * qdot;
        wdot[3*npt+i] -= 1 * qdot;
        wdot[6*npt+i] += 1 * qdot;
        wdot[0*npt+i] += 1 * qdot;

        /*reaction 19: H2O2 + O <=> OH + HO2 */
        phi_f = sc[7*npt+i]*sc[4*npt+i];
        k_f = k_f_s[18][i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[18][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= 1 * qdot;
        wdot[4*npt+i] -= 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;
        wdot[6*npt+i] += 1 * qdot;

        /*reaction 20: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[7*npt+i]*sc[5*npt+i];
        k_f = k_f_s[19][i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[2*npt+i];
        Kc = Kc_s[19][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= 1 * qdot;
        wdot[5*npt+i] -= 1 * qdot;
        wdot[6*npt+i] += 1 * qdot;
        wdot[2*npt+i] += 1 * qdot;

        /*reaction 21: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[7*npt+i]*sc[5*npt+i];
        k_f = k_f_s[20][i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[2*npt+i];
        Kc = Kc_s[20][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= 1 * qdot;
        wdot[5*npt+i] -= 1 * qdot;
        wdot[6*npt+i] += 1 * qdot;
        wdot[2*npt+i] += 1 * qdot;
    }
}

/*compute the reaction Jacobian */
void DWDOT(double * restrict  J, double * restrict  sc, double * Tp, int * consP)
{
    double c[9];

    for (int k=0; k<9; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    for (int k=0; k<9; k++) {
        J[90+k] *= 1.e-6;
    }

    /* dTdot/d[X] */
    for (int k=0; k<9; k++) {
        J[k*10+9] *= 1.e6;
    }

    return;
}

/*compute the reaction Jacobian */
void aJacobian(double * restrict J, double * restrict sc, double T, int consP)
{
    for (int i=0; i<100; i++) {
        J[i] = 0.0;
    }

    double wdot[9];
    for (int k=0; k<9; k++) {
        wdot[k] = 0.0;
    }

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];
    double invT2 = invT * invT;

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;
    double refCinv = 1 / refC;

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int k = 0; k < 9; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[9];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[9];
    speciesEnthalpy(h_RT, tc);

    double q, q_nocor, Corr, alpha;
    double dlnk0dT, dqdT;
    double dqdci, dcdc_fac, dqdc[9];
    double Pr, fPr, F, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    struct ReactionData* Ri;

    double Kc[21], k_f[21], k_r[21], k_0[21], dlnkfdT[21], dlnKcdT[21], dkrdT[21], phi_f[21], phi_r[21];

    phi_f[0] = sc[3]*sc[1];
    phi_r[0] = sc[4]*sc[5];
    Kc[0] = exp((g_RT[3] + g_RT[1]) - (g_RT[4] + g_RT[5]));
    dlnKcdT[0] = invT * (-(h_RT[3] + h_RT[1]) + (h_RT[4] + h_RT[5]));
    phi_f[1] = sc[0]*sc[4];
    phi_r[1] = sc[3]*sc[5];
    Kc[1] = exp((g_RT[0] + g_RT[4]) - (g_RT[3] + g_RT[5]));
    dlnKcdT[1] = invT * (-(h_RT[0] + h_RT[4]) + (h_RT[3] + h_RT[5]));
    phi_f[2] = sc[0]*sc[5];
    phi_r[2] = sc[3]*sc[2];
    Kc[2] = exp((g_RT[0] + g_RT[5]) - (g_RT[3] + g_RT[2]));
    dlnKcdT[2] = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[3] + h_RT[2]));
    phi_f[3] = sc[2]*sc[4];
    phi_r[3] = sc[5]*sc[5];
    Kc[3] = exp((g_RT[2] + g_RT[4]) - (g_RT[5] + g_RT[5]));
    dlnKcdT[3] = invT * (-(h_RT[2] + h_RT[4]) + (2*h_RT[5]));
    phi_f[4] = sc[0];
    phi_r[4] = sc[3]*sc[3];
    Kc[4] = refC * exp((g_RT[0]) - (g_RT[3] + g_RT[3]));
    dlnKcdT[4] = invT * (-(h_RT[0]) + (2*h_RT[3]) - 1);
    phi_f[5] = sc[4]*sc[4];
    phi_r[5] = sc[1];
    Kc[5] = refCinv * exp((g_RT[4] + g_RT[4]) - (g_RT[1]));
    dlnKcdT[5] = invT * (-(2*h_RT[4]) + (h_RT[1]) + 1);
    phi_f[6] = sc[3]*sc[4];
    phi_r[6] = sc[5];
    Kc[6] = refCinv * exp((g_RT[3] + g_RT[4]) - (g_RT[5]));
    dlnKcdT[6] = invT * (-(h_RT[3] + h_RT[4]) + (h_RT[5]) + 1);
    phi_f[7] = sc[3]*sc[5];
    phi_r[7] = sc[2];
    Kc[7] = refCinv * exp((g_RT[3] + g_RT[5]) - (g_RT[2]));
    dlnKcdT[7] = invT * (-(h_RT[3] + h_RT[5]) + (h_RT[2]) + 1);
    phi_f[8] = sc[3]*sc[1];
    phi_r[8] = sc[6];
    Kc[8] = refCinv * exp((g_RT[3] + g_RT[1]) - (g_RT[6]));
    dlnKcdT[8] = invT * (-(h_RT[3] + h_RT[1]) + (h_RT[6]) + 1);
    phi_f[9] = sc[3]*sc[6];
    phi_r[9] = sc[0]*sc[1];
    Kc[9] = exp((g_RT[3] + g_RT[6]) - (g_RT[0] + g_RT[1]));
    dlnKcdT[9] = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[0] + h_RT[1]));
    phi_f[10] = sc[3]*sc[6];
    phi_r[10] = sc[5]*sc[5];
    Kc[10] = exp((g_RT[3] + g_RT[6]) - (g_RT[5] + g_RT[5]));
    dlnKcdT[10] = invT * (-(h_RT[3] + h_RT[6]) + (2*h_RT[5]));
    phi_f[11] = sc[6]*sc[4];
    phi_r[11] = sc[1]*sc[5];
    Kc[11] = exp((g_RT[6] + g_RT[4]) - (g_RT[1] + g_RT[5]));
    dlnKcdT[11] = invT * (-(h_RT[6] + h_RT[4]) + (h_RT[1] + h_RT[5]));
    phi_f[12] = sc[6]*sc[5];
    phi_r[12] = sc[2]*sc[1];
    Kc[12] = exp((g_RT[6] + g_RT[5]) - (g_RT[2] + g_RT[1]));
    dlnKcdT[12] = invT * (-(h_RT[6] + h_RT[5]) + (h_RT[2] + h_RT[1]));
    phi_f[13] = sc[6]*sc[6];
    phi_r[13] = sc[7]*sc[1];
    Kc[13] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[1]));
    dlnKcdT[13] = invT * (-(2*h_RT[6]) + (h_RT[7] + h_RT[1]));
    phi_f[14] = sc[6]*sc[6];
    phi_r[14] = sc[7]*sc[1];
    Kc[14] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[1]));
    dlnKcdT[14] = invT * (-(2*h_RT[6]) + (h_RT[7] + h_RT[1]));
    phi_f[15] = sc[7];
    phi_r[15] = sc[5]*sc[5];
    Kc[15] = refC * exp((g_RT[7]) - (g_RT[5] + g_RT[5]));
    dlnKcdT[15] = invT * (-(h_RT[7]) + (2*h_RT[5]) - 1);
    phi_f[16] = sc[3]*sc[7];
    phi_r[16] = sc[2]*sc[5];
    Kc[16] = exp((g_RT[3] + g_RT[7]) - (g_RT[2] + g_RT[5]));
    dlnKcdT[16] = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[2] + h_RT[5]));
    phi_f[17] = sc[3]*sc[7];
    phi_r[17] = sc[0]*sc[6];
    Kc[17] = exp((g_RT[3] + g_RT[7]) - (g_RT[0] + g_RT[6]));
    dlnKcdT[17] = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[0] + h_RT[6]));
    phi_f[18] = sc[7]*sc[4];
    phi_r[18] = sc[6]*sc[5];
    Kc[18] = exp((g_RT[7] + g_RT[4]) - (g_RT[6] + g_RT[5]));
    dlnKcdT[18] = invT * (-(h_RT[7] + h_RT[4]) + (h_RT[6] + h_RT[5]));
    phi_f[19] = sc[7]*sc[5];
    phi_r[19] = sc[2]*sc[6];
    Kc[19] = exp((g_RT[7] + g_RT[5]) - (g_RT[2] + g_RT[6]));
    dlnKcdT[19] = invT * (-(h_RT[7] + h_RT[5]) + (h_RT[2] + h_RT[6]));
    phi_f[20] = sc[7]*sc[5];
    phi_r[20] = sc[2]*sc[6];
    Kc[20] = exp((g_RT[7] + g_RT[5]) - (g_RT[2] + g_RT[6]));
    dlnKcdT[20] = invT * (-(h_RT[7] + h_RT[5]) + (h_RT[2] + h_RT[6]));

    for (int i=0; i<21; ++i) {
        Ri = &(R[i]);
        k_f[i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[0] - Ri->activation_units * Ri->fwd_Ea * invT);
        dlnkfdT[i] = Ri->fwd_beta * invT + Ri->activation_units * Ri->fwd_Ea * invT2;
        k_r[i] = k_f[i] / Kc[i];
        dkrdT[i] = (dlnkfdT[i] - dlnKcdT[i])*k_r[i];
        k_0[i] = (Ri->is_PD ? Ri->low_A * exp(Ri->low_beta * tc[0] - Ri->activation_units * Ri->low_Ea * invT) : 0);
    }

    /*reaction 1: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    Ri = &(R[0]);
    /* rate of progress */
    q = k_f[0]*phi_f[0] - k_r[0]*phi_r[0];
    dqdT = dlnkfdT[0]*k_f[0]*phi_f[0] - dkrdT[0]*phi_r[0];
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[3] -= q; /* H */
    wdot[4] += q; /* O */
    wdot[5] += q; /* OH */
    /* d()/d[O2] */
    dqdci =  + k_f[0] *sc[3];
    J[11] += -1 * dqdci;          /* dwdot[O2]/d[O2] */
    J[13] += -1 * dqdci;          /* dwdot[H]/d[O2] */
    J[14] += 1 * dqdci;           /* dwdot[O]/d[O2] */
    J[15] += 1 * dqdci;           /* dwdot[OH]/d[O2] */
    /* d()/d[H] */
    dqdci =  + k_f[0] *sc[1];
    J[31] += -1 * dqdci;          /* dwdot[O2]/d[H] */
    J[33] += -1 * dqdci;          /* dwdot[H]/d[H] */
    J[34] += 1 * dqdci;           /* dwdot[O]/d[H] */
    J[35] += 1 * dqdci;           /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r[0]*sc[5];
    J[41] += -1 * dqdci;          /* dwdot[O2]/d[O] */
    J[43] += -1 * dqdci;          /* dwdot[H]/d[O] */
    J[44] += 1 * dqdci;           /* dwdot[O]/d[O] */
    J[45] += 1 * dqdci;           /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r[0]*sc[4];
    J[51] += -1 * dqdci;          /* dwdot[O2]/d[OH] */
    J[53] += -1 * dqdci;          /* dwdot[H]/d[OH] */
    J[54] += 1 * dqdci;           /* dwdot[O]/d[OH] */
    J[55] += 1 * dqdci;           /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[91] += -1 * dqdT;           /* dwdot[O2]/dT */
    J[93] += -1 * dqdT;           /* dwdot[H]/dT */
    J[94] += 1 * dqdT;            /* dwdot[O]/dT */
    J[95] += 1 * dqdT;            /* dwdot[OH]/dT */

    /*reaction 2: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    Ri = &(R[1]);
    /* rate of progress */
    q = k_f[1]*phi_f[1] - k_r[1]*phi_r[1];
    dqdT = dlnkfdT[1]*k_f[1]*phi_f[1] - dkrdT[1]*phi_r[1];
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[3] += q; /* H */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f[1] *sc[4];
    J[0] += -1 * dqdci;           /* dwdot[H2]/d[H2] */
    J[3] += 1 * dqdci;            /* dwdot[H]/d[H2] */
    J[4] += -1 * dqdci;           /* dwdot[O]/d[H2] */
    J[5] += 1 * dqdci;            /* dwdot[OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r[1]*sc[5];
    J[30] += -1 * dqdci;          /* dwdot[H2]/d[H] */
    J[33] += 1 * dqdci;           /* dwdot[H]/d[H] */
    J[34] += -1 * dqdci;          /* dwdot[O]/d[H] */
    J[35] += 1 * dqdci;           /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f[1] *sc[0];
    J[40] += -1 * dqdci;          /* dwdot[H2]/d[O] */
    J[43] += 1 * dqdci;           /* dwdot[H]/d[O] */
    J[44] += -1 * dqdci;          /* dwdot[O]/d[O] */
    J[45] += 1 * dqdci;           /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r[1]*sc[3];
    J[50] += -1 * dqdci;          /* dwdot[H2]/d[OH] */
    J[53] += 1 * dqdci;           /* dwdot[H]/d[OH] */
    J[54] += -1 * dqdci;          /* dwdot[O]/d[OH] */
    J[55] += 1 * dqdci;           /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[90] += -1 * dqdT;           /* dwdot[H2]/dT */
    J[93] += 1 * dqdT;            /* dwdot[H]/dT */
    J[94] += -1 * dqdT;           /* dwdot[O]/dT */
    J[95] += 1 * dqdT;            /* dwdot[OH]/dT */

    /*reaction 3: H2 + OH <=> H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    Ri = &(R[2]);
    /* rate of progress */
    q = k_f[2]*phi_f[2] - k_r[2]*phi_r[2];
    dqdT = dlnkfdT[2]*k_f[2]*phi_f[2] - dkrdT[2]*phi_r[2];
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[2] += q; /* H2O */
    wdot[3] += q; /* H */
    wdot[5] -= q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f[2] *sc[5];
    J[0] += -1 * dqdci;           /* dwdot[H2]/d[H2] */
    J[2] += 1 * dqdci;            /* dwdot[H2O]/d[H2] */
    J[3] += 1 * dqdci;            /* dwdot[H]/d[H2] */
    J[5] += -1 * dqdci;           /* dwdot[OH]/d[H2] */
    /* d()/d[H2O] */
    dqdci =  - k_r[2]*sc[3];
    J[20] += -1 * dqdci;          /* dwdot[H2]/d[H2O] */
    J[22] += 1 * dqdci;           /* dwdot[H2O]/d[H2O] */
    J[23] += 1 * dqdci;           /* dwdot[H]/d[H2O] */
    J[25] += -1 * dqdci;          /* dwdot[OH]/d[H2O] */
    /* d()/d[H] */
    dqdci =  - k_r[2]*sc[2];
    J[30] += -1 * dqdci;          /* dwdot[H2]/d[H] */
    J[32] += 1 * dqdci;           /* dwdot[H2O]/d[H] */
    J[33] += 1 * dqdci;           /* dwdot[H]/d[H] */
    J[35] += -1 * dqdci;          /* dwdot[OH]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f[2] *sc[0];
    J[50] += -1 * dqdci;          /* dwdot[H2]/d[OH] */
    J[52] += 1 * dqdci;           /* dwdot[H2O]/d[OH] */
    J[53] += 1 * dqdci;           /* dwdot[H]/d[OH] */
    J[55] += -1 * dqdci;          /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[90] += -1 * dqdT;           /* dwdot[H2]/dT */
    J[92] += 1 * dqdT;            /* dwdot[H2O]/dT */
    J[93] += 1 * dqdT;            /* dwdot[H]/dT */
    J[95] += -1 * dqdT;           /* dwdot[OH]/dT */

    /*reaction 4: O + H2O <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    Ri = &(R[3]);
    /* rate of progress */
    q = k_f[3]*phi_f[3] - k_r[3]*phi_r[3];
    dqdT = dlnkfdT[3]*k_f[3]*phi_f[3] - dkrdT[3]*phi_r[3];
    /* update wdot */
    wdot[2] -= q; /* H2O */
    wdot[4] -= q; /* O */
    wdot[5] += 2 * q; /* OH */
    /* d()/d[H2O] */
    dqdci =  + k_f[3] *sc[4];
    J[22] += -1 * dqdci;          /* dwdot[H2O]/d[H2O] */
    J[24] += -1 * dqdci;          /* dwdot[O]/d[H2O] */
    J[25] += 2 * dqdci;           /* dwdot[OH]/d[H2O] */
    /* d()/d[O] */
    dqdci =  + k_f[3] *sc[2];
    J[42] += -1 * dqdci;          /* dwdot[H2O]/d[O] */
    J[44] += -1 * dqdci;          /* dwdot[O]/d[O] */
    J[45] += 2 * dqdci;           /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r[3]*2*sc[5];
    J[52] += -1 * dqdci;          /* dwdot[H2O]/d[OH] */
    J[54] += -1 * dqdci;          /* dwdot[O]/d[OH] */
    J[55] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[92] += -1 * dqdT;           /* dwdot[H2O]/dT */
    J[94] += -1 * dqdT;           /* dwdot[O]/dT */
    J[95] += 2 * dqdT;            /* dwdot[OH]/dT */

    /*reaction 5: H2 + M <=> H + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    Ri = &(R[4]);
    /* 3-body correction factor */
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    /* rate of progress */
    q_nocor = k_f[4]*phi_f[4] - k_r[4]*phi_r[4];
    q = alpha * q_nocor;
    dqdT *= alpha * (dlnkfdT[4]*k_f[4]*phi_f[4] - dkrdT[4]*phi_r[4]);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[3] += 2 * q; /* H */
    /* for convenience */
    k_f[4] *= alpha;
    k_r[4] *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci =  1.5*q_nocor + k_f[4] ;
        J[0] += -1 * dqdci;           /* dwdot[H2]/d[H2] */
        J[3] += 2 * dqdci;            /* dwdot[H]/d[H2] */
        /* d()/d[H2O] */
        dqdci =  11*q_nocor;
        J[20] += -1 * dqdci;          /* dwdot[H2]/d[H2O] */
        J[23] += 2 * dqdci;           /* dwdot[H]/d[H2O] */
        /* d()/d[H] */
        dqdci =  - k_r[4]*2*sc[3];
        J[30] += -1 * dqdci;          /* dwdot[H2]/d[H] */
        J[33] += 2 * dqdci;           /* dwdot[H]/d[H] */
    } else {
        dqdc[0] = 2.5*q_nocor + k_f[4] ;
        dqdc[1] = q_nocor;
        dqdc[2] = 12*q_nocor;
        dqdc[3] = q_nocor - k_r[4]*2*sc[3];
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; ++k) {
            J[10*k+0] += -1 * dqdc[k];
            J[10*k+3] += 2 * dqdc[k];
        }
    }
    J[90] += -1 * dqdT; /* dwdot[H2]/dT */
    J[93] += 2 * dqdT; /* dwdot[H]/dT */

    /*reaction 6: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    Ri = &(R[5]);
    /* 3-body correction factor */
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    /* rate of progress */
    q_nocor = k_f[5]*phi_f[5] - k_r[5]*phi_r[5];
    q = alpha * q_nocor;
    dqdT *= alpha * (dlnkfdT[5]*k_f[5]*phi_f[5] - dkrdT[5]*phi_r[5]);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[4] -= 2 * q; /* O */
    /* for convenience */
    k_f[5] *= alpha;
    k_r[5] *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci =  1.5*q_nocor;
        J[1] += 1 * dqdci;            /* dwdot[O2]/d[H2] */
        J[4] += -2 * dqdci;           /* dwdot[O]/d[H2] */
        /* d()/d[O2] */
        dqdci =  - k_r[5];
        J[11] += 1 * dqdci;           /* dwdot[O2]/d[O2] */
        J[14] += -2 * dqdci;          /* dwdot[O]/d[O2] */
        /* d()/d[H2O] */
        dqdci =  11*q_nocor;
        J[21] += 1 * dqdci;           /* dwdot[O2]/d[H2O] */
        J[24] += -2 * dqdci;          /* dwdot[O]/d[H2O] */
        /* d()/d[O] */
        dqdci =  + k_f[5] *2*sc[4];
        J[41] += 1 * dqdci;           /* dwdot[O2]/d[O] */
        J[44] += -2 * dqdci;          /* dwdot[O]/d[O] */
    } else {
        dqdc[0] = 2.5*q_nocor;
        dqdc[1] = q_nocor - k_r[5];
        dqdc[2] = 12*q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor + k_f[5] *2*sc[4];
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; ++k) {
            J[10*k+1] += 1 * dqdc[k];
            J[10*k+4] += -2 * dqdc[k];
        }
    }
    J[91] += 1 * dqdT; /* dwdot[O2]/dT */
    J[94] += -2 * dqdT; /* dwdot[O]/dT */

    /*reaction 7: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    Ri = &(R[6]);
    /* 3-body correction factor */
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    /* rate of progress */
    q_nocor = k_f[6]*phi_f[6] - k_r[6]*phi_r[6];
    q = alpha * q_nocor;
    dqdT *= alpha * (dlnkfdT[6]*k_f[6]*phi_f[6] - dkrdT[6]*phi_r[6]);
    /* update wdot */
    wdot[3] -= q; /* H */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    /* for convenience */
    k_f[6] *= alpha;
    k_r[6] *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci =  1.5*q_nocor;
        J[3] += -1 * dqdci;           /* dwdot[H]/d[H2] */
        J[4] += -1 * dqdci;           /* dwdot[O]/d[H2] */
        J[5] += 1 * dqdci;            /* dwdot[OH]/d[H2] */
        /* d()/d[H2O] */
        dqdci =  11*q_nocor;
        J[23] += -1 * dqdci;          /* dwdot[H]/d[H2O] */
        J[24] += -1 * dqdci;          /* dwdot[O]/d[H2O] */
        J[25] += 1 * dqdci;           /* dwdot[OH]/d[H2O] */
        /* d()/d[H] */
        dqdci =  + k_f[6] *sc[4];
        J[33] += -1 * dqdci;          /* dwdot[H]/d[H] */
        J[34] += -1 * dqdci;          /* dwdot[O]/d[H] */
        J[35] += 1 * dqdci;           /* dwdot[OH]/d[H] */
        /* d()/d[O] */
        dqdci =  + k_f[6] *sc[3];
        J[43] += -1 * dqdci;          /* dwdot[H]/d[O] */
        J[44] += -1 * dqdci;          /* dwdot[O]/d[O] */
        J[45] += 1 * dqdci;           /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r[6];
        J[53] += -1 * dqdci;          /* dwdot[H]/d[OH] */
        J[54] += -1 * dqdci;          /* dwdot[O]/d[OH] */
        J[55] += 1 * dqdci;           /* dwdot[OH]/d[OH] */
    } else {
        dqdc[0] = 2.5*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = 12*q_nocor;
        dqdc[3] = q_nocor + k_f[6] *sc[4];
        dqdc[4] = q_nocor + k_f[6] *sc[3];
        dqdc[5] = q_nocor - k_r[6];
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; ++k) {
            J[10*k+3] += -1 * dqdc[k];
            J[10*k+4] += -1 * dqdc[k];
            J[10*k+5] += 1 * dqdc[k];
        }
    }
    J[93] += -1 * dqdT; /* dwdot[H]/dT */
    J[94] += -1 * dqdT; /* dwdot[O]/dT */
    J[95] += 1 * dqdT; /* dwdot[OH]/dT */

    /*reaction 8: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    Ri = &(R[7]);
    /* 3-body correction factor */
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    /* rate of progress */
    q_nocor = k_f[7]*phi_f[7] - k_r[7]*phi_r[7];
    q = alpha * q_nocor;
    dqdT *= alpha * (dlnkfdT[7]*k_f[7]*phi_f[7] - dkrdT[7]*phi_r[7]);
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[3] -= q; /* H */
    wdot[5] -= q; /* OH */
    /* for convenience */
    k_f[7] *= alpha;
    k_r[7] *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci =  1.5*q_nocor;
        J[2] += 1 * dqdci;            /* dwdot[H2O]/d[H2] */
        J[3] += -1 * dqdci;           /* dwdot[H]/d[H2] */
        J[5] += -1 * dqdci;           /* dwdot[OH]/d[H2] */
        /* d()/d[H2O] */
        dqdci =  11*q_nocor - k_r[7];
        J[22] += 1 * dqdci;           /* dwdot[H2O]/d[H2O] */
        J[23] += -1 * dqdci;          /* dwdot[H]/d[H2O] */
        J[25] += -1 * dqdci;          /* dwdot[OH]/d[H2O] */
        /* d()/d[H] */
        dqdci =  + k_f[7] *sc[5];
        J[32] += 1 * dqdci;           /* dwdot[H2O]/d[H] */
        J[33] += -1 * dqdci;          /* dwdot[H]/d[H] */
        J[35] += -1 * dqdci;          /* dwdot[OH]/d[H] */
        /* d()/d[OH] */
        dqdci =  + k_f[7] *sc[3];
        J[52] += 1 * dqdci;           /* dwdot[H2O]/d[OH] */
        J[53] += -1 * dqdci;          /* dwdot[H]/d[OH] */
        J[55] += -1 * dqdci;          /* dwdot[OH]/d[OH] */
    } else {
        dqdc[0] = 2.5*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = 12*q_nocor - k_r[7];
        dqdc[3] = q_nocor + k_f[7] *sc[5];
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor + k_f[7] *sc[3];
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; ++k) {
            J[10*k+2] += 1 * dqdc[k];
            J[10*k+3] += -1 * dqdc[k];
            J[10*k+5] += -1 * dqdc[k];
        }
    }
    J[92] += 1 * dqdT; /* dwdot[H2O]/dT */
    J[93] += -1 * dqdT; /* dwdot[H]/dT */
    J[95] += -1 * dqdT; /* dwdot[OH]/dT */

    /*reaction 9: H + O2 (+M) <=> HO2 (+M) */
    /*a pressure-fall-off reaction */
    Ri = &(R[8]);
    /* 3-body correction factor */
    alpha = mixture + sc[0] + 10*sc[2] + -0.22*sc[1];
    /* pressure-fall-off */
    Pr = Ri->phase_units * alpha / k_f[8] * k_0[8];
    fPr = Pr / (1.0+Pr);
    dlnk0dT = Ri->low_beta * invT + Ri->activation_units * Ri->low_Ea * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT[8]);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(Ri->troe_Tsss) > 1.e-100 ? (1-Ri->troe_a)*exp(-T/Ri->troe_Tsss) : 0);
    Fcent2 = (fabs(Ri->troe_Ts) > 1.e-100 ? Ri->troe_a * exp(-T/Ri->troe_Ts) : 0);
    Fcent3 = (Ri->troe_len == 4 ? exp(-Ri->troe_Tss * invT) : 0);
    Fcent = Fcent1 + Fcent2 + (Ri->troe_len == 4 ? Fcent3 : 0);
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( (fabs(Ri->troe_Tsss) > 1.e-100 ? -Fcent1/Ri->troe_Tsss : 0) + (fabs(Ri->troe_Ts) > 1.e-100 ? -Fcent2/Ri->troe_Ts : 0) + (Ri->troe_len == 4 ? Fcent3*Ri->troe_Tss*invT2 : 0) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* rate of progress */
    q_nocor = k_f[8]*phi_f[8] - k_r[8]*phi_r[8];
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr * (dlnkfdT[8]*k_f[8]*phi_f[8] - dkrdT[8]*phi_r[8]) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[3] -= q; /* H */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f[8] *= Corr;
    k_r[8] *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = dcdc_fac;
        J[1] += -1 * dqdci;           /* dwdot[O2]/d[H2] */
        J[3] += -1 * dqdci;           /* dwdot[H]/d[H2] */
        J[6] += 1 * dqdci;            /* dwdot[HO2]/d[H2] */
        /* d()/d[O2] */
        dqdci =  -0.22*dcdc_fac + k_f[8] *sc[3];
        J[11] += -1 * dqdci;          /* dwdot[O2]/d[O2] */
        J[13] += -1 * dqdci;          /* dwdot[H]/d[O2] */
        J[16] += 1 * dqdci;           /* dwdot[HO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci =  10*dcdc_fac;
        J[21] += -1 * dqdci;          /* dwdot[O2]/d[H2O] */
        J[23] += -1 * dqdci;          /* dwdot[H]/d[H2O] */
        J[26] += 1 * dqdci;           /* dwdot[HO2]/d[H2O] */
        /* d()/d[H] */
        dqdci =  + k_f[8] *sc[1];
        J[31] += -1 * dqdci;          /* dwdot[O2]/d[H] */
        J[33] += -1 * dqdci;          /* dwdot[H]/d[H] */
        J[36] += 1 * dqdci;           /* dwdot[HO2]/d[H] */
        /* d()/d[HO2] */
        dqdci =  - k_r[8];
        J[61] += -1 * dqdci;          /* dwdot[O2]/d[HO2] */
        J[63] += -1 * dqdci;          /* dwdot[H]/d[HO2] */
        J[66] += 1 * dqdci;           /* dwdot[HO2]/d[HO2] */
    } else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = 0.78*dcdc_fac + k_f[8] *sc[3];
        dqdc[2] = 11*dcdc_fac;
        dqdc[3] = dcdc_fac + k_f[8] *sc[1];
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac - k_r[8];
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        for (int k=0; k<9; ++k) {
            J[10*k+1] += -1 * dqdc[k];
            J[10*k+3] += -1 * dqdc[k];
            J[10*k+6] += 1 * dqdc[k];
        }
    }
    J[91] += -1 * dqdT; /* dwdot[O2]/dT */
    J[93] += -1 * dqdT; /* dwdot[H]/dT */
    J[96] += 1 * dqdT; /* dwdot[HO2]/dT */

    /*reaction 10: HO2 + H <=> H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    Ri = &(R[9]);
    /* rate of progress */
    q = k_f[9]*phi_f[9] - k_r[9]*phi_r[9];
    dqdT = dlnkfdT[9]*k_f[9]*phi_f[9] - dkrdT[9]*phi_r[9];
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] += q; /* O2 */
    wdot[3] -= q; /* H */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H2] */
    dqdci =  - k_r[9]*sc[1];
    J[0] += 1 * dqdci;            /* dwdot[H2]/d[H2] */
    J[1] += 1 * dqdci;            /* dwdot[O2]/d[H2] */
    J[3] += -1 * dqdci;           /* dwdot[H]/d[H2] */
    J[6] += -1 * dqdci;           /* dwdot[HO2]/d[H2] */
    /* d()/d[O2] */
    dqdci =  - k_r[9]*sc[0];
    J[10] += 1 * dqdci;           /* dwdot[H2]/d[O2] */
    J[11] += 1 * dqdci;           /* dwdot[O2]/d[O2] */
    J[13] += -1 * dqdci;          /* dwdot[H]/d[O2] */
    J[16] += -1 * dqdci;          /* dwdot[HO2]/d[O2] */
    /* d()/d[H] */
    dqdci =  + k_f[9] *sc[6];
    J[30] += 1 * dqdci;           /* dwdot[H2]/d[H] */
    J[31] += 1 * dqdci;           /* dwdot[O2]/d[H] */
    J[33] += -1 * dqdci;          /* dwdot[H]/d[H] */
    J[36] += -1 * dqdci;          /* dwdot[HO2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  + k_f[9] *sc[3];
    J[60] += 1 * dqdci;           /* dwdot[H2]/d[HO2] */
    J[61] += 1 * dqdci;           /* dwdot[O2]/d[HO2] */
    J[63] += -1 * dqdci;          /* dwdot[H]/d[HO2] */
    J[66] += -1 * dqdci;          /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[90] += 1 * dqdT;            /* dwdot[H2]/dT */
    J[91] += 1 * dqdT;            /* dwdot[O2]/dT */
    J[93] += -1 * dqdT;           /* dwdot[H]/dT */
    J[96] += -1 * dqdT;           /* dwdot[HO2]/dT */

    /*reaction 11: HO2 + H <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    Ri = &(R[10]);
    /* rate of progress */
    q = k_f[10]*phi_f[10] - k_r[10]*phi_r[10];
    dqdT = dlnkfdT[10]*k_f[10]*phi_f[10] - dkrdT[10]*phi_r[10];
    /* update wdot */
    wdot[3] -= q; /* H */
    wdot[5] += 2 * q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f[10] *sc[6];
    J[33] += -1 * dqdci;          /* dwdot[H]/d[H] */
    J[35] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[36] += -1 * dqdci;          /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r[10]*2*sc[5];
    J[53] += -1 * dqdci;          /* dwdot[H]/d[OH] */
    J[55] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[56] += -1 * dqdci;          /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f[10] *sc[3];
    J[63] += -1 * dqdci;          /* dwdot[H]/d[HO2] */
    J[65] += 2 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[66] += -1 * dqdci;          /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[93] += -1 * dqdT;           /* dwdot[H]/dT */
    J[95] += 2 * dqdT;            /* dwdot[OH]/dT */
    J[96] += -1 * dqdT;           /* dwdot[HO2]/dT */

    /*reaction 12: HO2 + O <=> O2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    Ri = &(R[11]);
    /* rate of progress */
    q = k_f[11]*phi_f[11] - k_r[11]*phi_r[11];
    dqdT = dlnkfdT[11]*k_f[11]*phi_f[11] - dkrdT[11]*phi_r[11];
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r[11]*sc[5];
    J[11] += 1 * dqdci;           /* dwdot[O2]/d[O2] */
    J[14] += -1 * dqdci;          /* dwdot[O]/d[O2] */
    J[15] += 1 * dqdci;           /* dwdot[OH]/d[O2] */
    J[16] += -1 * dqdci;          /* dwdot[HO2]/d[O2] */
    /* d()/d[O] */
    dqdci =  + k_f[11] *sc[6];
    J[41] += 1 * dqdci;           /* dwdot[O2]/d[O] */
    J[44] += -1 * dqdci;          /* dwdot[O]/d[O] */
    J[45] += 1 * dqdci;           /* dwdot[OH]/d[O] */
    J[46] += -1 * dqdci;          /* dwdot[HO2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r[11]*sc[1];
    J[51] += 1 * dqdci;           /* dwdot[O2]/d[OH] */
    J[54] += -1 * dqdci;          /* dwdot[O]/d[OH] */
    J[55] += 1 * dqdci;           /* dwdot[OH]/d[OH] */
    J[56] += -1 * dqdci;          /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f[11] *sc[4];
    J[61] += 1 * dqdci;           /* dwdot[O2]/d[HO2] */
    J[64] += -1 * dqdci;          /* dwdot[O]/d[HO2] */
    J[65] += 1 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[66] += -1 * dqdci;          /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] += 1 * dqdT;            /* dwdot[O2]/dT */
    J[94] += -1 * dqdT;           /* dwdot[O]/dT */
    J[95] += 1 * dqdT;            /* dwdot[OH]/dT */
    J[96] += -1 * dqdT;           /* dwdot[HO2]/dT */

    /*reaction 13: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    Ri = &(R[12]);
    /* rate of progress */
    q = k_f[12]*phi_f[12] - k_r[12]*phi_r[12];
    dqdT = dlnkfdT[12]*k_f[12]*phi_f[12] - dkrdT[12]*phi_r[12];
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[2] += q; /* H2O */
    wdot[5] -= q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r[12]*sc[2];
    J[11] += 1 * dqdci;           /* dwdot[O2]/d[O2] */
    J[12] += 1 * dqdci;           /* dwdot[H2O]/d[O2] */
    J[15] += -1 * dqdci;          /* dwdot[OH]/d[O2] */
    J[16] += -1 * dqdci;          /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r[12]*sc[1];
    J[21] += 1 * dqdci;           /* dwdot[O2]/d[H2O] */
    J[22] += 1 * dqdci;           /* dwdot[H2O]/d[H2O] */
    J[25] += -1 * dqdci;          /* dwdot[OH]/d[H2O] */
    J[26] += -1 * dqdci;          /* dwdot[HO2]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f[12] *sc[6];
    J[51] += 1 * dqdci;           /* dwdot[O2]/d[OH] */
    J[52] += 1 * dqdci;           /* dwdot[H2O]/d[OH] */
    J[55] += -1 * dqdci;          /* dwdot[OH]/d[OH] */
    J[56] += -1 * dqdci;          /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f[12] *sc[5];
    J[61] += 1 * dqdci;           /* dwdot[O2]/d[HO2] */
    J[62] += 1 * dqdci;           /* dwdot[H2O]/d[HO2] */
    J[65] += -1 * dqdci;          /* dwdot[OH]/d[HO2] */
    J[66] += -1 * dqdci;          /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] += 1 * dqdT;            /* dwdot[O2]/dT */
    J[92] += 1 * dqdT;            /* dwdot[H2O]/dT */
    J[95] += -1 * dqdT;           /* dwdot[OH]/dT */
    J[96] += -1 * dqdT;           /* dwdot[HO2]/dT */

    /*reaction 14: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    Ri = &(R[13]);
    /* rate of progress */
    q = k_f[13]*phi_f[13] - k_r[13]*phi_r[13];
    dqdT = dlnkfdT[13]*k_f[13]*phi_f[13] - dkrdT[13]*phi_r[13];
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r[13]*sc[7];
    J[11] += 1 * dqdci;           /* dwdot[O2]/d[O2] */
    J[16] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[17] += 1 * dqdci;           /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f[13] *2*sc[6];
    J[61] += 1 * dqdci;           /* dwdot[O2]/d[HO2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[67] += 1 * dqdci;           /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r[13]*sc[1];
    J[71] += 1 * dqdci;           /* dwdot[O2]/d[H2O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[77] += 1 * dqdci;           /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[91] += 1 * dqdT;            /* dwdot[O2]/dT */
    J[96] += -2 * dqdT;           /* dwdot[HO2]/dT */
    J[97] += 1 * dqdT;            /* dwdot[H2O2]/dT */

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    Ri = &(R[14]);
    /* rate of progress */
    q = k_f[14]*phi_f[14] - k_r[14]*phi_r[14];
    dqdT = dlnkfdT[14]*k_f[14]*phi_f[14] - dkrdT[14]*phi_r[14];
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r[14]*sc[7];
    J[11] += 1 * dqdci;           /* dwdot[O2]/d[O2] */
    J[16] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[17] += 1 * dqdci;           /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f[14] *2*sc[6];
    J[61] += 1 * dqdci;           /* dwdot[O2]/d[HO2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[67] += 1 * dqdci;           /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r[14]*sc[1];
    J[71] += 1 * dqdci;           /* dwdot[O2]/d[H2O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[77] += 1 * dqdci;           /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[91] += 1 * dqdT;            /* dwdot[O2]/dT */
    J[96] += -2 * dqdT;           /* dwdot[HO2]/dT */
    J[97] += 1 * dqdT;            /* dwdot[H2O2]/dT */

    /*reaction 16: H2O2 (+M) <=> OH + OH (+M) */
    /*a pressure-fall-off reaction */
    Ri = &(R[15]);
    /* 3-body correction factor */
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    /* pressure-fall-off */
    Pr = Ri->phase_units * alpha / k_f[15] * k_0[15];
    fPr = Pr / (1.0+Pr);
    dlnk0dT = Ri->low_beta * invT + Ri->activation_units * Ri->low_Ea * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT[15]);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(Ri->troe_Tsss) > 1.e-100 ? (1-Ri->troe_a)*exp(-T/Ri->troe_Tsss) : 0);
    Fcent2 = (fabs(Ri->troe_Ts) > 1.e-100 ? Ri->troe_a * exp(-T/Ri->troe_Ts) : 0);
    Fcent3 = (Ri->troe_len == 4 ? exp(-Ri->troe_Tss * invT) : 0);
    Fcent = Fcent1 + Fcent2 + (Ri->troe_len == 4 ? Fcent3 : 0);
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( (fabs(Ri->troe_Tsss) > 1.e-100 ? -Fcent1/Ri->troe_Tsss : 0) + (fabs(Ri->troe_Ts) > 1.e-100 ? -Fcent2/Ri->troe_Ts : 0) + (Ri->troe_len == 4 ? Fcent3*Ri->troe_Tss*invT2 : 0) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* rate of progress */
    q_nocor = k_f[15]*phi_f[15] - k_r[15]*phi_r[15];
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr * (dlnkfdT[15]*k_f[15]*phi_f[15] - dkrdT[15]*phi_r[15]) + dlnCorrdT*q;
    /* update wdot */
    wdot[5] += 2 * q; /* OH */
    wdot[7] -= q; /* H2O2 */
    /* for convenience */
    k_f[15] *= Corr;
    k_r[15] *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci =  1.5*dcdc_fac;
        J[5] += 2 * dqdci;            /* dwdot[OH]/d[H2] */
        J[7] += -1 * dqdci;           /* dwdot[H2O2]/d[H2] */
        /* d()/d[H2O] */
        dqdci =  11*dcdc_fac;
        J[25] += 2 * dqdci;           /* dwdot[OH]/d[H2O] */
        J[27] += -1 * dqdci;          /* dwdot[H2O2]/d[H2O] */
        /* d()/d[OH] */
        dqdci =  - k_r[15]*2*sc[5];
        J[55] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
        J[57] += -1 * dqdci;          /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O2] */
        dqdci =  + k_f[15] ;
        J[75] += 2 * dqdci;           /* dwdot[OH]/d[H2O2] */
        J[77] += -1 * dqdci;          /* dwdot[H2O2]/d[H2O2] */
    } else {
        dqdc[0] = 2.5*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = 12*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac - k_r[15]*2*sc[5];
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac + k_f[15] ;
        dqdc[8] = dcdc_fac;
        for (int k=0; k<9; ++k) {
            J[10*k+5] += 2 * dqdc[k];
            J[10*k+7] += -1 * dqdc[k];
        }
    }
    J[95] += 2 * dqdT; /* dwdot[OH]/dT */
    J[97] += -1 * dqdT; /* dwdot[H2O2]/dT */

    /*reaction 17: H2O2 + H <=> H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    Ri = &(R[16]);
    /* rate of progress */
    q = k_f[16]*phi_f[16] - k_r[16]*phi_r[16];
    dqdT = dlnkfdT[16]*k_f[16]*phi_f[16] - dkrdT[16]*phi_r[16];
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[3] -= q; /* H */
    wdot[5] += q; /* OH */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2O] */
    dqdci =  - k_r[16]*sc[5];
    J[22] += 1 * dqdci;           /* dwdot[H2O]/d[H2O] */
    J[23] += -1 * dqdci;          /* dwdot[H]/d[H2O] */
    J[25] += 1 * dqdci;           /* dwdot[OH]/d[H2O] */
    J[27] += -1 * dqdci;          /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H] */
    dqdci =  + k_f[16] *sc[7];
    J[32] += 1 * dqdci;           /* dwdot[H2O]/d[H] */
    J[33] += -1 * dqdci;          /* dwdot[H]/d[H] */
    J[35] += 1 * dqdci;           /* dwdot[OH]/d[H] */
    J[37] += -1 * dqdci;          /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r[16]*sc[2];
    J[52] += 1 * dqdci;           /* dwdot[H2O]/d[OH] */
    J[53] += -1 * dqdci;          /* dwdot[H]/d[OH] */
    J[55] += 1 * dqdci;           /* dwdot[OH]/d[OH] */
    J[57] += -1 * dqdci;          /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O2] */
    dqdci =  + k_f[16] *sc[3];
    J[72] += 1 * dqdci;           /* dwdot[H2O]/d[H2O2] */
    J[73] += -1 * dqdci;          /* dwdot[H]/d[H2O2] */
    J[75] += 1 * dqdci;           /* dwdot[OH]/d[H2O2] */
    J[77] += -1 * dqdci;          /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[92] += 1 * dqdT;            /* dwdot[H2O]/dT */
    J[93] += -1 * dqdT;           /* dwdot[H]/dT */
    J[95] += 1 * dqdT;            /* dwdot[OH]/dT */
    J[97] += -1 * dqdT;           /* dwdot[H2O2]/dT */

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    Ri = &(R[17]);
    /* rate of progress */
    q = k_f[17]*phi_f[17] - k_r[17]*phi_r[17];
    dqdT = dlnkfdT[17]*k_f[17]*phi_f[17] - dkrdT[17]*phi_r[17];
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[3] -= q; /* H */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2] */
    dqdci =  - k_r[17]*sc[6];
    J[0] += 1 * dqdci;            /* dwdot[H2]/d[H2] */
    J[3] += -1 * dqdci;           /* dwdot[H]/d[H2] */
    J[6] += 1 * dqdci;            /* dwdot[HO2]/d[H2] */
    J[7] += -1 * dqdci;           /* dwdot[H2O2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f[17] *sc[7];
    J[30] += 1 * dqdci;           /* dwdot[H2]/d[H] */
    J[33] += -1 * dqdci;          /* dwdot[H]/d[H] */
    J[36] += 1 * dqdci;           /* dwdot[HO2]/d[H] */
    J[37] += -1 * dqdci;          /* dwdot[H2O2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  - k_r[17]*sc[0];
    J[60] += 1 * dqdci;           /* dwdot[H2]/d[HO2] */
    J[63] += -1 * dqdci;          /* dwdot[H]/d[HO2] */
    J[66] += 1 * dqdci;           /* dwdot[HO2]/d[HO2] */
    J[67] += -1 * dqdci;          /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f[17] *sc[3];
    J[70] += 1 * dqdci;           /* dwdot[H2]/d[H2O2] */
    J[73] += -1 * dqdci;          /* dwdot[H]/d[H2O2] */
    J[76] += 1 * dqdci;           /* dwdot[HO2]/d[H2O2] */
    J[77] += -1 * dqdci;          /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[90] += 1 * dqdT;            /* dwdot[H2]/dT */
    J[93] += -1 * dqdT;           /* dwdot[H]/dT */
    J[96] += 1 * dqdT;            /* dwdot[HO2]/dT */
    J[97] += -1 * dqdT;           /* dwdot[H2O2]/dT */

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    Ri = &(R[18]);
    /* rate of progress */
    q = k_f[18]*phi_f[18] - k_r[18]*phi_r[18];
    dqdT = dlnkfdT[18]*k_f[18]*phi_f[18] - dkrdT[18]*phi_r[18];
    /* update wdot */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[O] */
    dqdci =  + k_f[18] *sc[7];
    J[44] += -1 * dqdci;          /* dwdot[O]/d[O] */
    J[45] += 1 * dqdci;           /* dwdot[OH]/d[O] */
    J[46] += 1 * dqdci;           /* dwdot[HO2]/d[O] */
    J[47] += -1 * dqdci;          /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r[18]*sc[6];
    J[54] += -1 * dqdci;          /* dwdot[O]/d[OH] */
    J[55] += 1 * dqdci;           /* dwdot[OH]/d[OH] */
    J[56] += 1 * dqdci;           /* dwdot[HO2]/d[OH] */
    J[57] += -1 * dqdci;          /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r[18]*sc[5];
    J[64] += -1 * dqdci;          /* dwdot[O]/d[HO2] */
    J[65] += 1 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[66] += 1 * dqdci;           /* dwdot[HO2]/d[HO2] */
    J[67] += -1 * dqdci;          /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f[18] *sc[4];
    J[74] += -1 * dqdci;          /* dwdot[O]/d[H2O2] */
    J[75] += 1 * dqdci;           /* dwdot[OH]/d[H2O2] */
    J[76] += 1 * dqdci;           /* dwdot[HO2]/d[H2O2] */
    J[77] += -1 * dqdci;          /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[94] += -1 * dqdT;           /* dwdot[O]/dT */
    J[95] += 1 * dqdT;            /* dwdot[OH]/dT */
    J[96] += 1 * dqdT;            /* dwdot[HO2]/dT */
    J[97] += -1 * dqdT;           /* dwdot[H2O2]/dT */

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    Ri = &(R[19]);
    /* rate of progress */
    q = k_f[19]*phi_f[19] - k_r[19]*phi_r[19];
    dqdT = dlnkfdT[19]*k_f[19]*phi_f[19] - dkrdT[19]*phi_r[19];
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2O] */
    dqdci =  - k_r[19]*sc[6];
    J[22] += 1 * dqdci;           /* dwdot[H2O]/d[H2O] */
    J[25] += -1 * dqdci;          /* dwdot[OH]/d[H2O] */
    J[26] += 1 * dqdci;           /* dwdot[HO2]/d[H2O] */
    J[27] += -1 * dqdci;          /* dwdot[H2O2]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f[19] *sc[7];
    J[52] += 1 * dqdci;           /* dwdot[H2O]/d[OH] */
    J[55] += -1 * dqdci;          /* dwdot[OH]/d[OH] */
    J[56] += 1 * dqdci;           /* dwdot[HO2]/d[OH] */
    J[57] += -1 * dqdci;          /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r[19]*sc[2];
    J[62] += 1 * dqdci;           /* dwdot[H2O]/d[HO2] */
    J[65] += -1 * dqdci;          /* dwdot[OH]/d[HO2] */
    J[66] += 1 * dqdci;           /* dwdot[HO2]/d[HO2] */
    J[67] += -1 * dqdci;          /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f[19] *sc[5];
    J[72] += 1 * dqdci;           /* dwdot[H2O]/d[H2O2] */
    J[75] += -1 * dqdci;          /* dwdot[OH]/d[H2O2] */
    J[76] += 1 * dqdci;           /* dwdot[HO2]/d[H2O2] */
    J[77] += -1 * dqdci;          /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[92] += 1 * dqdT;            /* dwdot[H2O]/dT */
    J[95] += -1 * dqdT;           /* dwdot[OH]/dT */
    J[96] += 1 * dqdT;            /* dwdot[HO2]/dT */
    J[97] += -1 * dqdT;           /* dwdot[H2O2]/dT */

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    Ri = &(R[20]);
    /* rate of progress */
    q = k_f[20]*phi_f[20] - k_r[20]*phi_r[20];
    dqdT = dlnkfdT[20]*k_f[20]*phi_f[20] - dkrdT[20]*phi_r[20];
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2O] */
    dqdci =  - k_r[20]*sc[6];
    J[22] += 1 * dqdci;           /* dwdot[H2O]/d[H2O] */
    J[25] += -1 * dqdci;          /* dwdot[OH]/d[H2O] */
    J[26] += 1 * dqdci;           /* dwdot[HO2]/d[H2O] */
    J[27] += -1 * dqdci;          /* dwdot[H2O2]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f[20] *sc[7];
    J[52] += 1 * dqdci;           /* dwdot[H2O]/d[OH] */
    J[55] += -1 * dqdci;          /* dwdot[OH]/d[OH] */
    J[56] += 1 * dqdci;           /* dwdot[HO2]/d[OH] */
    J[57] += -1 * dqdci;          /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r[20]*sc[2];
    J[62] += 1 * dqdci;           /* dwdot[H2O]/d[HO2] */
    J[65] += -1 * dqdci;          /* dwdot[OH]/d[HO2] */
    J[66] += 1 * dqdci;           /* dwdot[HO2]/d[HO2] */
    J[67] += -1 * dqdci;          /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f[20] *sc[5];
    J[72] += 1 * dqdci;           /* dwdot[H2O]/d[H2O2] */
    J[75] += -1 * dqdci;          /* dwdot[OH]/d[H2O2] */
    J[76] += 1 * dqdci;           /* dwdot[HO2]/d[H2O2] */
    J[77] += -1 * dqdci;          /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[92] += 1 * dqdT;            /* dwdot[H2O]/dT */
    J[95] += -1 * dqdT;           /* dwdot[OH]/dT */
    J[96] += 1 * dqdT;            /* dwdot[HO2]/dT */
    J[97] += -1 * dqdT;           /* dwdot[H2O2]/dT */

    double c_R[9], dcRdT[9], e_RT[9];
    double * eh_RT;
    if (consP) {
        cp_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        eh_RT = &h_RT[0];
    }
    else {
        cv_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        speciesInternalEnergy(e_RT, tc);
        eh_RT = &e_RT[0];
    }

    double cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 9; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[90+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 9; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 9; ++m) {
            dehmixdc += eh_RT[m]*J[k*10+m];
        }
        J[k*10+9] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[99] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}


/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void dcvpRdT(double * restrict  species, double * restrict  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +8.24944174e-04
            -1.62860306e-06 * tc[1]
            -2.84263030e-10 * tc[2]
            +1.65394890e-12 * tc[3];
        /*species 1: O2 */
        species[1] =
            +1.12748635e-03
            -1.15123009e-06 * tc[1]
            +3.94163169e-09 * tc[2]
            -3.50742157e-12 * tc[3];
        /*species 2: H2O */
        species[2] =
            +3.47498246e-03
            -1.27093927e-05 * tc[1]
            +2.09057438e-08 * tc[2]
            -1.00263539e-11 * tc[3];
        /*species 3: H */
        species[3] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 4: O */
        species[4] =
            -1.63816649e-03
            +4.84206340e-06 * tc[1]
            -4.80852957e-09 * tc[2]
            +1.55627854e-12 * tc[3];
        /*species 5: OH */
        species[5] =
            -3.22544939e-03
            +1.30552938e-05 * tc[1]
            -1.73956093e-08 * tc[2]
            +8.24949516e-12 * tc[3];
        /*species 6: HO2 */
        species[6] =
            -4.74912051e-03
            +4.23165782e-05 * tc[1]
            -7.28291682e-08 * tc[2]
            +3.71690050e-11 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +6.56922581e-03
            -2.97002516e-07 * tc[1]
            -1.38774166e-08 * tc[2]
            +9.88605900e-12 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.40824000e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77942000e-12 * tc[3];
    } else {
        /*species 0: H2 */
        species[0] =
            +7.00064411e-04
            -1.12676574e-07 * tc[1]
            -2.76947345e-11 * tc[2]
            +6.33100716e-15 * tc[3];
        /*species 1: O2 */
        species[1] =
            +6.13519689e-04
            -2.51768398e-07 * tc[1]
            +5.32584444e-11 * tc[2]
            -4.54574124e-15 * tc[3];
        /*species 2: H2O */
        species[2] =
            +3.05629289e-03
            -1.74605202e-06 * tc[1]
            +3.60298917e-10 * tc[2]
            -2.55664715e-14 * tc[3];
        /*species 3: H */
        species[3] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 4: O */
        species[4] =
            -2.75506191e-05
            -6.20560670e-09 * tc[1]
            +1.36532023e-11 * tc[2]
            -1.74722060e-15 * tc[3];
        /*species 5: OH */
        species[5] =
            +1.05650448e-03
            -5.18165516e-07 * tc[1]
            +9.15656022e-11 * tc[2]
            -5.32783504e-15 * tc[3];
        /*species 6: HO2 */
        species[6] =
            +2.23982013e-03
            -1.26731630e-06 * tc[1]
            +3.42739110e-10 * tc[2]
            -4.31634140e-14 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +4.33613639e-03
            -2.94937764e-06 * tc[1]
            +7.04671071e-10 * tc[2]
            -5.72661424e-14 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.48797700e-03
            -1.13695220e-06 * tc[1]
            +3.02911200e-10 * tc[2]
            -2.70134040e-14 * tc[3];
    }
    return;
}


/*compute the progress rate for each reaction */
void progressRate(double * restrict  qdot, double * restrict  sc, double T)
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

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */

    double invT = 1.0 / tc[1];

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;
    double refCinv = 1 / refC;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 9; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    struct ReactionData* Ri;
    if (T != T_save)
    {
        T_save = T;

        for (int i=0; i<21; ++i) {
            Ri = &(R[i]);
            k_f_save[i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[0] - Ri->activation_units * Ri->fwd_Ea * invT);
        };

        Kc_save[0] = exp((g_RT[3] + g_RT[1]) - (g_RT[4] + g_RT[5]));
        Kc_save[1] = exp((g_RT[4] + g_RT[0]) - (g_RT[3] + g_RT[5]));
        Kc_save[2] = exp((g_RT[0] + g_RT[5]) - (g_RT[2] + g_RT[3]));
        Kc_save[3] = exp((g_RT[4] + g_RT[2]) - (g_RT[5] + g_RT[5]));
        Kc_save[4] = refC * exp((g_RT[0]) - (g_RT[3] + g_RT[3]));
        Kc_save[5] = 1.0 / (refC) * exp((g_RT[4] + g_RT[4]) - (g_RT[1]));
        Kc_save[6] = 1.0 / (refC) * exp((g_RT[4] + g_RT[3]) - (g_RT[5]));
        Kc_save[7] = 1.0 / (refC) * exp((g_RT[3] + g_RT[5]) - (g_RT[2]));
        Kc_save[8] = 1.0 / (refC) * exp((g_RT[3] + g_RT[1]) - (g_RT[6]));
        Kc_save[9] = exp((g_RT[6] + g_RT[3]) - (g_RT[0] + g_RT[1]));
        Kc_save[10] = exp((g_RT[6] + g_RT[3]) - (g_RT[5] + g_RT[5]));
        Kc_save[11] = exp((g_RT[6] + g_RT[4]) - (g_RT[1] + g_RT[5]));
        Kc_save[12] = exp((g_RT[6] + g_RT[5]) - (g_RT[2] + g_RT[1]));
        Kc_save[13] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[1]));
        Kc_save[14] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[1]));
        Kc_save[15] = refC * exp((g_RT[7]) - (g_RT[5] + g_RT[5]));
        Kc_save[16] = exp((g_RT[7] + g_RT[3]) - (g_RT[2] + g_RT[5]));
        Kc_save[17] = exp((g_RT[7] + g_RT[3]) - (g_RT[6] + g_RT[0]));
        Kc_save[18] = exp((g_RT[7] + g_RT[4]) - (g_RT[5] + g_RT[6]));
        Kc_save[19] = exp((g_RT[7] + g_RT[5]) - (g_RT[6] + g_RT[2]));
        Kc_save[20] = exp((g_RT[7] + g_RT[5]) - (g_RT[6] + g_RT[2]));
    }

    /*reaction 1: H + O2 <=> O + OH */
    phi_f = sc[3]*sc[1];
    k_f = k_f_save[0];
    q_f = phi_f * k_f;
    phi_r = sc[4]*sc[5];
    Kc = Kc_save[0];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[0] = q_f - q_r;

    /*reaction 2: O + H2 <=> H + OH */
    phi_f = sc[4]*sc[0];
    k_f = k_f_save[1];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = Kc_save[1];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[1] = q_f - q_r;

    /*reaction 3: H2 + OH <=> H2O + H */
    phi_f = sc[0]*sc[5];
    k_f = k_f_save[2];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[3];
    Kc = Kc_save[2];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[2] = q_f - q_r;

    /*reaction 4: O + H2O <=> OH + OH */
    phi_f = sc[4]*sc[2];
    k_f = k_f_save[3];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[5];
    Kc = Kc_save[3];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[3] = q_f - q_r;

    /*reaction 5: H2 + M <=> H + H + M */
    phi_f = sc[0];
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    k_f = alpha * k_f_save[4];
    q_f = phi_f * k_f;
    phi_r = sc[3]*sc[3];
    Kc = Kc_save[4];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[4] = q_f - q_r;

    /*reaction 6: O + O + M <=> O2 + M */
    phi_f = sc[4]*sc[4];
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    k_f = alpha * k_f_save[5];
    q_f = phi_f * k_f;
    phi_r = sc[1];
    Kc = Kc_save[5];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[5] = q_f - q_r;

    /*reaction 7: O + H + M <=> OH + M */
    phi_f = sc[4]*sc[3];
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    k_f = alpha * k_f_save[6];
    q_f = phi_f * k_f;
    phi_r = sc[5];
    Kc = Kc_save[6];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[6] = q_f - q_r;

    /*reaction 8: H + OH + M <=> H2O + M */
    phi_f = sc[3]*sc[5];
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    k_f = alpha * k_f_save[7];
    q_f = phi_f * k_f;
    phi_r = sc[2];
    Kc = Kc_save[7];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[7] = q_f - q_r;

    /*reaction 9: H + O2 (+M) <=> HO2 (+M) */
    phi_f = sc[3]*sc[1];
    alpha = mixture + sc[0] + 10*sc[2] + -0.22*sc[1];
    k_f = k_f_save[8];
    Ri = &(R[8]);
    redP = alpha / k_f * Ri->phase_units * Ri->low_A * exp(Ri->low_beta * tc[0] - Ri->activation_units * Ri->low_Ea *invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (fabs(Ri->troe_Tsss) > 1.e-100 ? (1-Ri->troe_a)*exp(-T/Ri->troe_Tsss) : 0) 
        + (fabs(Ri->troe_Ts) > 1.e-100 ? Ri->troe_a * exp(-T/Ri->troe_Ts) : 0) 
        + (Ri->troe_len == 4 ? exp(-Ri->troe_Tss * invT) : 0) );
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[6];
    Kc = Kc_save[8];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[8] = q_f - q_r;

    /*reaction 10: HO2 + H <=> H2 + O2 */
    phi_f = sc[6]*sc[3];
    k_f = k_f_save[9];
    q_f = phi_f * k_f;
    phi_r = sc[0]*sc[1];
    Kc = Kc_save[9];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[9] = q_f - q_r;

    /*reaction 11: HO2 + H <=> OH + OH */
    phi_f = sc[6]*sc[3];
    k_f = k_f_save[10];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[5];
    Kc = Kc_save[10];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[10] = q_f - q_r;

    /*reaction 12: HO2 + O <=> O2 + OH */
    phi_f = sc[6]*sc[4];
    k_f = k_f_save[11];
    q_f = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    Kc = Kc_save[11];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[11] = q_f - q_r;

    /*reaction 13: HO2 + OH <=> H2O + O2 */
    phi_f = sc[6]*sc[5];
    k_f = k_f_save[12];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[1];
    Kc = Kc_save[12];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[12] = q_f - q_r;

    /*reaction 14: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[6]*sc[6];
    k_f = k_f_save[13];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[1];
    Kc = Kc_save[13];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[13] = q_f - q_r;

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[6]*sc[6];
    k_f = k_f_save[14];
    q_f = phi_f * k_f;
    phi_r = sc[7]*sc[1];
    Kc = Kc_save[14];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[14] = q_f - q_r;

    /*reaction 16: H2O2 (+M) <=> OH + OH (+M) */
    phi_f = sc[7];
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    k_f = k_f_save[15];
    Ri = &(R[15]);
    redP = alpha / k_f * Ri->phase_units * Ri->low_A * exp(Ri->low_beta * tc[0] - Ri->activation_units * Ri->low_Ea *invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (fabs(Ri->troe_Tsss) > 1.e-100 ? (1-Ri->troe_a)*exp(-T/Ri->troe_Tsss) : 0) 
        + (fabs(Ri->troe_Ts) > 1.e-100 ? Ri->troe_a * exp(-T/Ri->troe_Ts) : 0) 
        + (Ri->troe_len == 4 ? exp(-Ri->troe_Tss * invT) : 0) );
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[5];
    Kc = Kc_save[15];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[15] = q_f - q_r;

    /*reaction 17: H2O2 + H <=> H2O + OH */
    phi_f = sc[7]*sc[3];
    k_f = k_f_save[16];
    q_f = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = Kc_save[16];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[16] = q_f - q_r;

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    phi_f = sc[7]*sc[3];
    k_f = k_f_save[17];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = Kc_save[17];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[17] = q_f - q_r;

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    phi_f = sc[7]*sc[4];
    k_f = k_f_save[18];
    q_f = phi_f * k_f;
    phi_r = sc[5]*sc[6];
    Kc = Kc_save[18];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[18] = q_f - q_r;

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    phi_f = sc[7]*sc[5];
    k_f = k_f_save[19];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[2];
    Kc = Kc_save[19];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[19] = q_f - q_r;

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    phi_f = sc[7]*sc[5];
    k_f = k_f_save[20];
    q_f = phi_f * k_f;
    phi_r = sc[6]*sc[2];
    Kc = Kc_save[20];
    k_r = k_f / Kc;
    q_r = phi_r * k_r;
    qdot[20] = q_f - q_r;

    return;
}


/*compute the progress rate for each reaction */
void progressRateFR(double * restrict  q_f, double * restrict  q_r, double * restrict  sc, double T)
{

    int id; /*loop counter */
    double mixture;                 /*mixture concentration */
    double g_RT[9];                /*Gibbs free energy */
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

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */

    double invT = 1.0 / tc[1];

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    /*compute the mixture concentration */
    mixture = 0.0;
    for (id = 0; id < 9; ++id) {
        mixture += sc[id];
    }

    /*compute the Gibbs free energy */
    gibbs(g_RT, tc);

    struct ReactionData* Ri;
    if (T != T_save)
    {
        T_save = T;

        for (int i=0; i<21; ++i) {
            Ri = &(R[i]);
            k_f_save[i] = Ri->prefactor_units * Ri->fwd_A * exp(Ri->fwd_beta * tc[0] - Ri->activation_units * Ri->fwd_Ea * invT);
        };
        Kc_save[0] = exp((g_RT[3] + g_RT[1]) - (g_RT[4] + g_RT[5]));
        Kc_save[1] = exp((g_RT[4] + g_RT[0]) - (g_RT[3] + g_RT[5]));
        Kc_save[2] = exp((g_RT[0] + g_RT[5]) - (g_RT[2] + g_RT[3]));
        Kc_save[3] = exp((g_RT[4] + g_RT[2]) - (g_RT[5] + g_RT[5]));
        Kc_save[4] = refC * exp((g_RT[0]) - (g_RT[3] + g_RT[3]));
        Kc_save[5] = 1.0 / (refC) * exp((g_RT[4] + g_RT[4]) - (g_RT[1]));
        Kc_save[6] = 1.0 / (refC) * exp((g_RT[4] + g_RT[3]) - (g_RT[5]));
        Kc_save[7] = 1.0 / (refC) * exp((g_RT[3] + g_RT[5]) - (g_RT[2]));
        Kc_save[8] = 1.0 / (refC) * exp((g_RT[3] + g_RT[1]) - (g_RT[6]));
        Kc_save[9] = exp((g_RT[6] + g_RT[3]) - (g_RT[0] + g_RT[1]));
        Kc_save[10] = exp((g_RT[6] + g_RT[3]) - (g_RT[5] + g_RT[5]));
        Kc_save[11] = exp((g_RT[6] + g_RT[4]) - (g_RT[1] + g_RT[5]));
        Kc_save[12] = exp((g_RT[6] + g_RT[5]) - (g_RT[2] + g_RT[1]));
        Kc_save[13] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[1]));
        Kc_save[14] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[1]));
        Kc_save[15] = refC * exp((g_RT[7]) - (g_RT[5] + g_RT[5]));
        Kc_save[16] = exp((g_RT[7] + g_RT[3]) - (g_RT[2] + g_RT[5]));
        Kc_save[17] = exp((g_RT[7] + g_RT[3]) - (g_RT[6] + g_RT[0]));
        Kc_save[18] = exp((g_RT[7] + g_RT[4]) - (g_RT[5] + g_RT[6]));
        Kc_save[19] = exp((g_RT[7] + g_RT[5]) - (g_RT[6] + g_RT[2]));
        Kc_save[20] = exp((g_RT[7] + g_RT[5]) - (g_RT[6] + g_RT[2]));
    }

    /*reaction 1: H + O2 <=> O + OH */
    phi_f = sc[3]*sc[1];
    k_f = k_f_save[0];
    q_f[0] = phi_f * k_f;
    phi_r = sc[4]*sc[5];
    Kc = Kc_save[0];
    k_r = k_f / Kc;
    q_r[0] = phi_r * k_r;

    /*reaction 2: O + H2 <=> H + OH */
    phi_f = sc[4]*sc[0];
    k_f = k_f_save[1];
    q_f[1] = phi_f * k_f;
    phi_r = sc[3]*sc[5];
    Kc = Kc_save[1];
    k_r = k_f / Kc;
    q_r[1] = phi_r * k_r;

    /*reaction 3: H2 + OH <=> H2O + H */
    phi_f = sc[0]*sc[5];
    k_f = k_f_save[2];
    q_f[2] = phi_f * k_f;
    phi_r = sc[2]*sc[3];
    Kc = Kc_save[2];
    k_r = k_f / Kc;
    q_r[2] = phi_r * k_r;

    /*reaction 4: O + H2O <=> OH + OH */
    phi_f = sc[4]*sc[2];
    k_f = k_f_save[3];
    q_f[3] = phi_f * k_f;
    phi_r = sc[5]*sc[5];
    Kc = Kc_save[3];
    k_r = k_f / Kc;
    q_r[3] = phi_r * k_r;

    /*reaction 5: H2 + M <=> H + H + M */
    phi_f = sc[0];
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    k_f = alpha * k_f_save[4];
    q_f[4] = phi_f * k_f;
    phi_r = sc[3]*sc[3];
    Kc = Kc_save[4];
    k_r = k_f / Kc;
    q_r[4] = phi_r * k_r;

    /*reaction 6: O + O + M <=> O2 + M */
    phi_f = sc[4]*sc[4];
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    k_f = alpha * k_f_save[5];
    q_f[5] = phi_f * k_f;
    phi_r = sc[1];
    Kc = Kc_save[5];
    k_r = k_f / Kc;
    q_r[5] = phi_r * k_r;

    /*reaction 7: O + H + M <=> OH + M */
    phi_f = sc[4]*sc[3];
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    k_f = alpha * k_f_save[6];
    q_f[6] = phi_f * k_f;
    phi_r = sc[5];
    Kc = Kc_save[6];
    k_r = k_f / Kc;
    q_r[6] = phi_r * k_r;

    /*reaction 8: H + OH + M <=> H2O + M */
    phi_f = sc[3]*sc[5];
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    k_f = alpha * k_f_save[7];
    q_f[7] = phi_f * k_f;
    phi_r = sc[2];
    Kc = Kc_save[7];
    k_r = k_f / Kc;
    q_r[7] = phi_r * k_r;

    /*reaction 9: H + O2 (+M) <=> HO2 (+M) */
    phi_f = sc[3]*sc[1];
    alpha = mixture + sc[0] + 10*sc[2] + -0.22*sc[1];
    k_f = k_f_save[8];
    Ri = &(R[8]);
    redP = alpha / k_f * Ri->phase_units * Ri->low_A * exp(Ri->low_beta * tc[0] - Ri->activation_units * Ri->low_Ea *invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (fabs(Ri->troe_Tsss) > 1.e-100 ? (1-Ri->troe_a)*exp(-T/Ri->troe_Tsss) : 0) 
        + (fabs(Ri->troe_Ts) > 1.e-100 ? Ri->troe_a * exp(-T/Ri->troe_Ts) : 0) 
        + (Ri->troe_len == 4 ? exp(-Ri->troe_Tss * invT) : 0) );
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[8] = phi_f * k_f;
    phi_r = sc[6];
    Kc = Kc_save[8];
    k_r = k_f / Kc;
    q_r[8] = phi_r * k_r;

    /*reaction 10: HO2 + H <=> H2 + O2 */
    phi_f = sc[6]*sc[3];
    k_f = k_f_save[9];
    q_f[9] = phi_f * k_f;
    phi_r = sc[0]*sc[1];
    Kc = Kc_save[9];
    k_r = k_f / Kc;
    q_r[9] = phi_r * k_r;

    /*reaction 11: HO2 + H <=> OH + OH */
    phi_f = sc[6]*sc[3];
    k_f = k_f_save[10];
    q_f[10] = phi_f * k_f;
    phi_r = sc[5]*sc[5];
    Kc = Kc_save[10];
    k_r = k_f / Kc;
    q_r[10] = phi_r * k_r;

    /*reaction 12: HO2 + O <=> O2 + OH */
    phi_f = sc[6]*sc[4];
    k_f = k_f_save[11];
    q_f[11] = phi_f * k_f;
    phi_r = sc[1]*sc[5];
    Kc = Kc_save[11];
    k_r = k_f / Kc;
    q_r[11] = phi_r * k_r;

    /*reaction 13: HO2 + OH <=> H2O + O2 */
    phi_f = sc[6]*sc[5];
    k_f = k_f_save[12];
    q_f[12] = phi_f * k_f;
    phi_r = sc[2]*sc[1];
    Kc = Kc_save[12];
    k_r = k_f / Kc;
    q_r[12] = phi_r * k_r;

    /*reaction 14: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[6]*sc[6];
    k_f = k_f_save[13];
    q_f[13] = phi_f * k_f;
    phi_r = sc[7]*sc[1];
    Kc = Kc_save[13];
    k_r = k_f / Kc;
    q_r[13] = phi_r * k_r;

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    phi_f = sc[6]*sc[6];
    k_f = k_f_save[14];
    q_f[14] = phi_f * k_f;
    phi_r = sc[7]*sc[1];
    Kc = Kc_save[14];
    k_r = k_f / Kc;
    q_r[14] = phi_r * k_r;

    /*reaction 16: H2O2 (+M) <=> OH + OH (+M) */
    phi_f = sc[7];
    alpha = mixture + 1.5*sc[0] + 11*sc[2];
    k_f = k_f_save[15];
    Ri = &(R[15]);
    redP = alpha / k_f * Ri->phase_units * Ri->low_A * exp(Ri->low_beta * tc[0] - Ri->activation_units * Ri->low_Ea *invT);
    F = redP / (1 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (fabs(Ri->troe_Tsss) > 1.e-100 ? (1-Ri->troe_a)*exp(-T/Ri->troe_Tsss) : 0) 
        + (fabs(Ri->troe_Ts) > 1.e-100 ? Ri->troe_a * exp(-T/Ri->troe_Ts) : 0) 
        + (Ri->troe_len == 4 ? exp(-Ri->troe_Tss * invT) : 0) );
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    F *= F_troe;
    k_f *= F;
    q_f[15] = phi_f * k_f;
    phi_r = sc[5]*sc[5];
    Kc = Kc_save[15];
    k_r = k_f / Kc;
    q_r[15] = phi_r * k_r;

    /*reaction 17: H2O2 + H <=> H2O + OH */
    phi_f = sc[7]*sc[3];
    k_f = k_f_save[16];
    q_f[16] = phi_f * k_f;
    phi_r = sc[2]*sc[5];
    Kc = Kc_save[16];
    k_r = k_f / Kc;
    q_r[16] = phi_r * k_r;

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    phi_f = sc[7]*sc[3];
    k_f = k_f_save[17];
    q_f[17] = phi_f * k_f;
    phi_r = sc[6]*sc[0];
    Kc = Kc_save[17];
    k_r = k_f / Kc;
    q_r[17] = phi_r * k_r;

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    phi_f = sc[7]*sc[4];
    k_f = k_f_save[18];
    q_f[18] = phi_f * k_f;
    phi_r = sc[5]*sc[6];
    Kc = Kc_save[18];
    k_r = k_f / Kc;
    q_r[18] = phi_r * k_r;

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    phi_f = sc[7]*sc[5];
    k_f = k_f_save[19];
    q_f[19] = phi_f * k_f;
    phi_r = sc[6]*sc[2];
    Kc = Kc_save[19];
    k_r = k_f / Kc;
    q_r[19] = phi_r * k_r;

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    phi_f = sc[7]*sc[5];
    k_f = k_f_save[20];
    q_f[20] = phi_f * k_f;
    phi_r = sc[6]*sc[2];
    Kc = Kc_save[20];
    k_r = k_f / Kc;
    q_r[20] = phi_r * k_r;

    return;
}


/*compute the equilibrium constants for each reaction */
void equilibriumConstants(double * restrict kc, double * restrict  g_RT, double T)
{
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    /*reaction 1: H + O2 <=> O + OH */
    kc[0] = exp((g_RT[3] + g_RT[1]) - (g_RT[4] + g_RT[5]));

    /*reaction 2: O + H2 <=> H + OH */
    kc[1] = exp((g_RT[4] + g_RT[0]) - (g_RT[3] + g_RT[5]));

    /*reaction 3: H2 + OH <=> H2O + H */
    kc[2] = exp((g_RT[0] + g_RT[5]) - (g_RT[2] + g_RT[3]));

    /*reaction 4: O + H2O <=> OH + OH */
    kc[3] = exp((g_RT[4] + g_RT[2]) - (g_RT[5] + g_RT[5]));

    /*reaction 5: H2 + M <=> H + H + M */
    kc[4] = refC * exp((g_RT[0]) - (g_RT[3] + g_RT[3]));

    /*reaction 6: O + O + M <=> O2 + M */
    kc[5] = 1.0 / (refC) * exp((g_RT[4] + g_RT[4]) - (g_RT[1]));

    /*reaction 7: O + H + M <=> OH + M */
    kc[6] = 1.0 / (refC) * exp((g_RT[4] + g_RT[3]) - (g_RT[5]));

    /*reaction 8: H + OH + M <=> H2O + M */
    kc[7] = 1.0 / (refC) * exp((g_RT[3] + g_RT[5]) - (g_RT[2]));

    /*reaction 9: H + O2 (+M) <=> HO2 (+M) */
    kc[8] = 1.0 / (refC) * exp((g_RT[3] + g_RT[1]) - (g_RT[6]));

    /*reaction 10: HO2 + H <=> H2 + O2 */
    kc[9] = exp((g_RT[6] + g_RT[3]) - (g_RT[0] + g_RT[1]));

    /*reaction 11: HO2 + H <=> OH + OH */
    kc[10] = exp((g_RT[6] + g_RT[3]) - (g_RT[5] + g_RT[5]));

    /*reaction 12: HO2 + O <=> O2 + OH */
    kc[11] = exp((g_RT[6] + g_RT[4]) - (g_RT[1] + g_RT[5]));

    /*reaction 13: HO2 + OH <=> H2O + O2 */
    kc[12] = exp((g_RT[6] + g_RT[5]) - (g_RT[2] + g_RT[1]));

    /*reaction 14: HO2 + HO2 <=> H2O2 + O2 */
    kc[13] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[1]));

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    kc[14] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[1]));

    /*reaction 16: H2O2 (+M) <=> OH + OH (+M) */
    kc[15] = refC * exp((g_RT[7]) - (g_RT[5] + g_RT[5]));

    /*reaction 17: H2O2 + H <=> H2O + OH */
    kc[16] = exp((g_RT[7] + g_RT[3]) - (g_RT[2] + g_RT[5]));

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    kc[17] = exp((g_RT[7] + g_RT[3]) - (g_RT[6] + g_RT[0]));

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    kc[18] = exp((g_RT[7] + g_RT[4]) - (g_RT[5] + g_RT[6]));

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    kc[19] = exp((g_RT[7] + g_RT[5]) - (g_RT[6] + g_RT[2]));

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    kc[20] = exp((g_RT[7] + g_RT[5]) - (g_RT[6] + g_RT[2]));

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void gibbs(double * restrict  species, double * restrict  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            -1.012520870000000e+03 * invT
            +6.592218400000000e+00
            -3.298124310000000e+00 * tc[0]
            -4.124720870000000e-04 * tc[1]
            +1.357169215000000e-07 * tc[2]
            +7.896195275000000e-12 * tc[3]
            -2.067436120000000e-14 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.005249020000000e+03 * invT
            -2.821801190000000e+00
            -3.212936400000000e+00 * tc[0]
            -5.637431750000000e-04 * tc[1]
            +9.593584116666666e-08 * tc[2]
            -1.094897691666667e-10 * tc[3]
            +4.384276960000000e-14 * tc[4];
        /*species 2: H2O */
        species[2] =
            -3.020811330000000e+04 * invT
            +7.966096399999998e-01
            -3.386842490000000e+00 * tc[0]
            -1.737491230000000e-03 * tc[1]
            +1.059116055000000e-06 * tc[2]
            -5.807151058333333e-10 * tc[3]
            +1.253294235000000e-13 * tc[4];
        /*species 3: H */
        species[3] =
            +2.547162700000000e+04 * invT
            +2.960117608000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.914764450000000e+04 * invT
            -1.756619999999964e-02
            -2.946428780000000e+00 * tc[0]
            +8.190832450000000e-04 * tc[1]
            -4.035052833333333e-07 * tc[2]
            +1.335702658333333e-10 * tc[3]
            -1.945348180000000e-14 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.346309130000000e+03 * invT
            +4.815738570000000e+00
            -4.125305610000000e+00 * tc[0]
            +1.612724695000000e-03 * tc[1]
            -1.087941151666667e-06 * tc[2]
            +4.832113691666666e-10 * tc[3]
            -1.031186895000000e-13 * tc[4];
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
            -1.766314650000000e+04 * invT
            -3.396609550000000e+00
            -3.388753650000000e+00 * tc[0]
            -3.284612905000000e-03 * tc[1]
            +2.475020966666667e-08 * tc[2]
            +3.854837933333333e-10 * tc[3]
            -1.235757375000000e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.020900000000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041200000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427500000000e-13 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -8.350339970000000e+02 * invT
            +4.346533540000000e+00
            -2.991423370000000e+00 * tc[0]
            -3.500322055000000e-04 * tc[1]
            +9.389714483333333e-09 * tc[2]
            +7.692981816666667e-13 * tc[3]
            -7.913758950000000e-17 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.233930180000000e+03 * invT
            +5.084126000000002e-01
            -3.697578190000000e+00 * tc[0]
            -3.067598445000000e-04 * tc[1]
            +2.098069983333333e-08 * tc[2]
            -1.479401233333333e-12 * tc[3]
            +5.682176550000000e-17 * tc[4];
        /*species 2: H2O */
        species[2] =
            -2.989920900000000e+04 * invT
            -4.190671200000001e+00
            -2.672145610000000e+00 * tc[0]
            -1.528146445000000e-03 * tc[1]
            +1.455043351666667e-07 * tc[2]
            -1.000830325000000e-11 * tc[3]
            +3.195808935000000e-16 * tc[4];
        /*species 3: H */
        species[3] =
            +2.547162700000000e+04 * invT
            +2.960117638000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.923080270000000e+04 * invT
            -2.378248450000000e+00
            -2.542059660000000e+00 * tc[0]
            +1.377530955000000e-05 * tc[1]
            +5.171338916666667e-10 * tc[2]
            -3.792556183333333e-13 * tc[3]
            +2.184025750000000e-17 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.683628750000000e+03 * invT
            -2.836911870000000e+00
            -2.864728860000000e+00 * tc[0]
            -5.282522400000000e-04 * tc[1]
            +4.318045966666667e-08 * tc[2]
            -2.543488950000000e-12 * tc[3]
            +6.659793800000000e-17 * tc[4];
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
            -1.800696090000000e+04 * invT
            +4.072029891000000e+00
            -4.573166850000000e+00 * tc[0]
            -2.168068195000000e-03 * tc[1]
            +2.457814700000000e-07 * tc[2]
            -1.957419641666667e-11 * tc[3]
            +7.158267800000000e-16 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439885000000000e-04 * tc[1]
            +9.474601666666666e-08 * tc[2]
            -8.414199999999999e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
    }
    return;
}


/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void helmholtz(double * restrict  species, double * restrict  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            -1.01252087e+03 * invT
            +5.59221840e+00
            -3.29812431e+00 * tc[0]
            -4.12472087e-04 * tc[1]
            +1.35716922e-07 * tc[2]
            +7.89619527e-12 * tc[3]
            -2.06743612e-14 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.00524902e+03 * invT
            -3.82180119e+00
            -3.21293640e+00 * tc[0]
            -5.63743175e-04 * tc[1]
            +9.59358412e-08 * tc[2]
            -1.09489769e-10 * tc[3]
            +4.38427696e-14 * tc[4];
        /*species 2: H2O */
        species[2] =
            -3.02081133e+04 * invT
            -2.03390360e-01
            -3.38684249e+00 * tc[0]
            -1.73749123e-03 * tc[1]
            +1.05911606e-06 * tc[2]
            -5.80715106e-10 * tc[3]
            +1.25329424e-13 * tc[4];
        /*species 3: H */
        species[3] =
            +2.54716270e+04 * invT
            +1.96011761e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.91476445e+04 * invT
            -1.01756620e+00
            -2.94642878e+00 * tc[0]
            +8.19083245e-04 * tc[1]
            -4.03505283e-07 * tc[2]
            +1.33570266e-10 * tc[3]
            -1.94534818e-14 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.34630913e+03 * invT
            +3.81573857e+00
            -4.12530561e+00 * tc[0]
            +1.61272470e-03 * tc[1]
            -1.08794115e-06 * tc[2]
            +4.83211369e-10 * tc[3]
            -1.03118689e-13 * tc[4];
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
            -1.76631465e+04 * invT
            -4.39660955e+00
            -3.38875365e+00 * tc[0]
            -3.28461290e-03 * tc[1]
            +2.47502097e-08 * tc[2]
            +3.85483793e-10 * tc[3]
            -1.23575738e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.02090000e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120000e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242750e-13 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -8.35033997e+02 * invT
            +3.34653354e+00
            -2.99142337e+00 * tc[0]
            -3.50032206e-04 * tc[1]
            +9.38971448e-09 * tc[2]
            +7.69298182e-13 * tc[3]
            -7.91375895e-17 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.23393018e+03 * invT
            -4.91587400e-01
            -3.69757819e+00 * tc[0]
            -3.06759845e-04 * tc[1]
            +2.09806998e-08 * tc[2]
            -1.47940123e-12 * tc[3]
            +5.68217655e-17 * tc[4];
        /*species 2: H2O */
        species[2] =
            -2.98992090e+04 * invT
            -5.19067120e+00
            -2.67214561e+00 * tc[0]
            -1.52814644e-03 * tc[1]
            +1.45504335e-07 * tc[2]
            -1.00083033e-11 * tc[3]
            +3.19580894e-16 * tc[4];
        /*species 3: H */
        species[3] =
            +2.54716270e+04 * invT
            +1.96011764e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.92308027e+04 * invT
            -3.37824845e+00
            -2.54205966e+00 * tc[0]
            +1.37753096e-05 * tc[1]
            +5.17133892e-10 * tc[2]
            -3.79255618e-13 * tc[3]
            +2.18402575e-17 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.68362875e+03 * invT
            -3.83691187e+00
            -2.86472886e+00 * tc[0]
            -5.28252240e-04 * tc[1]
            +4.31804597e-08 * tc[2]
            -2.54348895e-12 * tc[3]
            +6.65979380e-17 * tc[4];
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
            -1.80069609e+04 * invT
            +3.07202989e+00
            -4.57316685e+00 * tc[0]
            -2.16806820e-03 * tc[1]
            +2.45781470e-07 * tc[2]
            -1.95741964e-11 * tc[3]
            +7.15826780e-16 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.22797700e+02 * invT
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
void cv_R(double * restrict  species, double * restrict  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +2.29812431e+00
            +8.24944174e-04 * tc[1]
            -8.14301529e-07 * tc[2]
            -9.47543433e-11 * tc[3]
            +4.13487224e-13 * tc[4];
        /*species 1: O2 */
        species[1] =
            +2.21293640e+00
            +1.12748635e-03 * tc[1]
            -5.75615047e-07 * tc[2]
            +1.31387723e-09 * tc[3]
            -8.76855392e-13 * tc[4];
        /*species 2: H2O */
        species[2] =
            +2.38684249e+00
            +3.47498246e-03 * tc[1]
            -6.35469633e-06 * tc[2]
            +6.96858127e-09 * tc[3]
            -2.50658847e-12 * tc[4];
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +1.94642878e+00
            -1.63816649e-03 * tc[1]
            +2.42103170e-06 * tc[2]
            -1.60284319e-09 * tc[3]
            +3.89069636e-13 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +2.38875365e+00
            +6.56922581e-03 * tc[1]
            -1.48501258e-07 * tc[2]
            -4.62580552e-09 * tc[3]
            +2.47151475e-12 * tc[4];
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
            +1.99142337e+00
            +7.00064411e-04 * tc[1]
            -5.63382869e-08 * tc[2]
            -9.23157818e-12 * tc[3]
            +1.58275179e-15 * tc[4];
        /*species 1: O2 */
        species[1] =
            +2.69757819e+00
            +6.13519689e-04 * tc[1]
            -1.25884199e-07 * tc[2]
            +1.77528148e-11 * tc[3]
            -1.13643531e-15 * tc[4];
        /*species 2: H2O */
        species[2] =
            +1.67214561e+00
            +3.05629289e-03 * tc[1]
            -8.73026011e-07 * tc[2]
            +1.20099639e-10 * tc[3]
            -6.39161787e-15 * tc[4];
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +1.54205966e+00
            -2.75506191e-05 * tc[1]
            -3.10280335e-09 * tc[2]
            +4.55106742e-12 * tc[3]
            -4.36805150e-16 * tc[4];
        /*species 5: OH */
        species[5] =
            +1.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.57316685e+00
            +4.33613639e-03 * tc[1]
            -1.47468882e-06 * tc[2]
            +2.34890357e-10 * tc[3]
            -1.43165356e-14 * tc[4];
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
void cp_R(double * restrict  species, double * restrict  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +3.29812431e+00
            +8.24944174e-04 * tc[1]
            -8.14301529e-07 * tc[2]
            -9.47543433e-11 * tc[3]
            +4.13487224e-13 * tc[4];
        /*species 1: O2 */
        species[1] =
            +3.21293640e+00
            +1.12748635e-03 * tc[1]
            -5.75615047e-07 * tc[2]
            +1.31387723e-09 * tc[3]
            -8.76855392e-13 * tc[4];
        /*species 2: H2O */
        species[2] =
            +3.38684249e+00
            +3.47498246e-03 * tc[1]
            -6.35469633e-06 * tc[2]
            +6.96858127e-09 * tc[3]
            -2.50658847e-12 * tc[4];
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.94642878e+00
            -1.63816649e-03 * tc[1]
            +2.42103170e-06 * tc[2]
            -1.60284319e-09 * tc[3]
            +3.89069636e-13 * tc[4];
        /*species 5: OH */
        species[5] =
            +4.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.38875365e+00
            +6.56922581e-03 * tc[1]
            -1.48501258e-07 * tc[2]
            -4.62580552e-09 * tc[3]
            +2.47151475e-12 * tc[4];
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
            +2.99142337e+00
            +7.00064411e-04 * tc[1]
            -5.63382869e-08 * tc[2]
            -9.23157818e-12 * tc[3]
            +1.58275179e-15 * tc[4];
        /*species 1: O2 */
        species[1] =
            +3.69757819e+00
            +6.13519689e-04 * tc[1]
            -1.25884199e-07 * tc[2]
            +1.77528148e-11 * tc[3]
            -1.13643531e-15 * tc[4];
        /*species 2: H2O */
        species[2] =
            +2.67214561e+00
            +3.05629289e-03 * tc[1]
            -8.73026011e-07 * tc[2]
            +1.20099639e-10 * tc[3]
            -6.39161787e-15 * tc[4];
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.54205966e+00
            -2.75506191e-05 * tc[1]
            -3.10280335e-09 * tc[2]
            +4.55106742e-12 * tc[3]
            -4.36805150e-16 * tc[4];
        /*species 5: OH */
        species[5] =
            +2.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.57316685e+00
            +4.33613639e-03 * tc[1]
            -1.47468882e-06 * tc[2]
            +2.34890357e-10 * tc[3]
            -1.43165356e-14 * tc[4];
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
void speciesInternalEnergy(double * restrict  species, double * restrict  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +2.29812431e+00
            +4.12472087e-04 * tc[1]
            -2.71433843e-07 * tc[2]
            -2.36885858e-11 * tc[3]
            +8.26974448e-14 * tc[4]
            -1.01252087e+03 * invT;
        /*species 1: O2 */
        species[1] =
            +2.21293640e+00
            +5.63743175e-04 * tc[1]
            -1.91871682e-07 * tc[2]
            +3.28469308e-10 * tc[3]
            -1.75371078e-13 * tc[4]
            -1.00524902e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +2.38684249e+00
            +1.73749123e-03 * tc[1]
            -2.11823211e-06 * tc[2]
            +1.74214532e-09 * tc[3]
            -5.01317694e-13 * tc[4]
            -3.02081133e+04 * invT;
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +1.94642878e+00
            -8.19083245e-04 * tc[1]
            +8.07010567e-07 * tc[2]
            -4.00710797e-10 * tc[3]
            +7.78139272e-14 * tc[4]
            +2.91476445e+04 * invT;
        /*species 5: OH */
        species[5] =
            +3.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.34630913e+03 * invT;
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
            +2.38875365e+00
            +3.28461290e-03 * tc[1]
            -4.95004193e-08 * tc[2]
            -1.15645138e-09 * tc[3]
            +4.94302950e-13 * tc[4]
            -1.76631465e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +1.99142337e+00
            +3.50032206e-04 * tc[1]
            -1.87794290e-08 * tc[2]
            -2.30789455e-12 * tc[3]
            +3.16550358e-16 * tc[4]
            -8.35033997e+02 * invT;
        /*species 1: O2 */
        species[1] =
            +2.69757819e+00
            +3.06759845e-04 * tc[1]
            -4.19613997e-08 * tc[2]
            +4.43820370e-12 * tc[3]
            -2.27287062e-16 * tc[4]
            -1.23393018e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +1.67214561e+00
            +1.52814644e-03 * tc[1]
            -2.91008670e-07 * tc[2]
            +3.00249098e-11 * tc[3]
            -1.27832357e-15 * tc[4]
            -2.98992090e+04 * invT;
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +1.54205966e+00
            -1.37753096e-05 * tc[1]
            -1.03426778e-09 * tc[2]
            +1.13776685e-12 * tc[3]
            -8.73610300e-17 * tc[4]
            +2.92308027e+04 * invT;
        /*species 5: OH */
        species[5] =
            +1.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.68362875e+03 * invT;
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
            +3.57316685e+00
            +2.16806820e-03 * tc[1]
            -4.91562940e-07 * tc[2]
            +5.87225893e-11 * tc[3]
            -2.86330712e-15 * tc[4]
            -1.80069609e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEnthalpy(double * restrict  species, double * restrict  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +3.29812431e+00
            +4.12472087e-04 * tc[1]
            -2.71433843e-07 * tc[2]
            -2.36885858e-11 * tc[3]
            +8.26974448e-14 * tc[4]
            -1.01252087e+03 * invT;
        /*species 1: O2 */
        species[1] =
            +3.21293640e+00
            +5.63743175e-04 * tc[1]
            -1.91871682e-07 * tc[2]
            +3.28469308e-10 * tc[3]
            -1.75371078e-13 * tc[4]
            -1.00524902e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +3.38684249e+00
            +1.73749123e-03 * tc[1]
            -2.11823211e-06 * tc[2]
            +1.74214532e-09 * tc[3]
            -5.01317694e-13 * tc[4]
            -3.02081133e+04 * invT;
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +2.94642878e+00
            -8.19083245e-04 * tc[1]
            +8.07010567e-07 * tc[2]
            -4.00710797e-10 * tc[3]
            +7.78139272e-14 * tc[4]
            +2.91476445e+04 * invT;
        /*species 5: OH */
        species[5] =
            +4.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.34630913e+03 * invT;
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
            +3.38875365e+00
            +3.28461290e-03 * tc[1]
            -4.95004193e-08 * tc[2]
            -1.15645138e-09 * tc[3]
            +4.94302950e-13 * tc[4]
            -1.76631465e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142337e+00
            +3.50032206e-04 * tc[1]
            -1.87794290e-08 * tc[2]
            -2.30789455e-12 * tc[3]
            +3.16550358e-16 * tc[4]
            -8.35033997e+02 * invT;
        /*species 1: O2 */
        species[1] =
            +3.69757819e+00
            +3.06759845e-04 * tc[1]
            -4.19613997e-08 * tc[2]
            +4.43820370e-12 * tc[3]
            -2.27287062e-16 * tc[4]
            -1.23393018e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +2.67214561e+00
            +1.52814644e-03 * tc[1]
            -2.91008670e-07 * tc[2]
            +3.00249098e-11 * tc[3]
            -1.27832357e-15 * tc[4]
            -2.98992090e+04 * invT;
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +2.54205966e+00
            -1.37753096e-05 * tc[1]
            -1.03426778e-09 * tc[2]
            +1.13776685e-12 * tc[3]
            -8.73610300e-17 * tc[4]
            +2.92308027e+04 * invT;
        /*species 5: OH */
        species[5] =
            +2.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.68362875e+03 * invT;
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
            +4.57316685e+00
            +2.16806820e-03 * tc[1]
            -4.91562940e-07 * tc[2]
            +5.87225893e-11 * tc[3]
            -2.86330712e-15 * tc[4]
            -1.80069609e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEntropy(double * restrict  species, double * restrict  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +3.29812431e+00 * tc[0]
            +8.24944174e-04 * tc[1]
            -4.07150765e-07 * tc[2]
            -3.15847811e-11 * tc[3]
            +1.03371806e-13 * tc[4]
            -3.29409409e+00 ;
        /*species 1: O2 */
        species[1] =
            +3.21293640e+00 * tc[0]
            +1.12748635e-03 * tc[1]
            -2.87807523e-07 * tc[2]
            +4.37959077e-10 * tc[3]
            -2.19213848e-13 * tc[4]
            +6.03473759e+00 ;
        /*species 2: H2O */
        species[2] =
            +3.38684249e+00 * tc[0]
            +3.47498246e-03 * tc[1]
            -3.17734817e-06 * tc[2]
            +2.32286042e-09 * tc[3]
            -6.26647117e-13 * tc[4]
            +2.59023285e+00 ;
        /*species 3: H */
        species[3] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117608e-01 ;
        /*species 4: O */
        species[4] =
            +2.94642878e+00 * tc[0]
            -1.63816649e-03 * tc[1]
            +1.21051585e-06 * tc[2]
            -5.34281063e-10 * tc[3]
            +9.72674090e-14 * tc[4]
            +2.96399498e+00 ;
        /*species 5: OH */
        species[5] =
            +4.12530561e+00 * tc[0]
            -3.22544939e-03 * tc[1]
            +3.26382346e-06 * tc[2]
            -1.93284548e-09 * tc[3]
            +5.15593447e-13 * tc[4]
            -6.90432960e-01 ;
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
            +3.38875365e+00 * tc[0]
            +6.56922581e-03 * tc[1]
            -7.42506290e-08 * tc[2]
            -1.54193517e-09 * tc[3]
            +6.17878688e-13 * tc[4]
            +6.78536320e+00 ;
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
            +2.99142337e+00 * tc[0]
            +7.00064411e-04 * tc[1]
            -2.81691434e-08 * tc[2]
            -3.07719273e-12 * tc[3]
            +3.95687948e-16 * tc[4]
            -1.35511017e+00 ;
        /*species 1: O2 */
        species[1] =
            +3.69757819e+00 * tc[0]
            +6.13519689e-04 * tc[1]
            -6.29420995e-08 * tc[2]
            +5.91760493e-12 * tc[3]
            -2.84108828e-16 * tc[4]
            +3.18916559e+00 ;
        /*species 2: H2O */
        species[2] =
            +2.67214561e+00 * tc[0]
            +3.05629289e-03 * tc[1]
            -4.36513005e-07 * tc[2]
            +4.00332130e-11 * tc[3]
            -1.59790447e-15 * tc[4]
            +6.86281681e+00 ;
        /*species 3: H */
        species[3] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117638e-01 ;
        /*species 4: O */
        species[4] =
            +2.54205966e+00 * tc[0]
            -2.75506191e-05 * tc[1]
            -1.55140167e-09 * tc[2]
            +1.51702247e-12 * tc[3]
            -1.09201287e-16 * tc[4]
            +4.92030811e+00 ;
        /*species 5: OH */
        species[5] =
            +2.86472886e+00 * tc[0]
            +1.05650448e-03 * tc[1]
            -1.29541379e-07 * tc[2]
            +1.01739558e-11 * tc[3]
            -3.32989690e-16 * tc[4]
            +5.70164073e+00 ;
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
            +4.57316685e+00 * tc[0]
            +4.33613639e-03 * tc[1]
            -7.37344410e-07 * tc[2]
            +7.82967857e-11 * tc[3]
            -3.57913390e-15 * tc[4]
            +5.01136959e-01 ;
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
void molecularWeight(double * restrict  wt)
{
    wt[0] = 2.015940; /*H2 */
    wt[1] = 31.998800; /*O2 */
    wt[2] = 18.015340; /*H2O */
    wt[3] = 1.007970; /*H */
    wt[4] = 15.999400; /*O */
    wt[5] = 17.007370; /*OH */
    wt[6] = 33.006770; /*HO2 */
    wt[7] = 34.014740; /*H2O2 */
    wt[8] = 28.013400; /*N2 */

    return;
}
/* get temperature given internal energy in mass units and mass fracs */
void GET_T_GIVEN_EY(double * restrict  e, double * restrict  y, int * iwrk, double * restrict  rwrk, double * restrict  t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const double tol  = 1.e-12;
#else
    const int maxiter = 200;
    const double tol  = 1.e-6;
#endif
    double ein  = *e;
    double tmin = 250;/*max lower bound for thermo def */
    double tmax = 3500;/*min upper bound for thermo def */
    double e1,emin,emax,cv,t1,dt;
    int i;/* loop counter */
    CKUBMS(&tmin, y, iwrk, rwrk, &emin);
    CKUBMS(&tmax, y, iwrk, rwrk, &emax);
    if (ein < emin) {
        /*Linear Extrapolation below tmin */
        CKCVBS(&tmin, y, iwrk, rwrk, &cv);
        *t = tmin - (emin-ein)/cv;
        *ierr = 1;
        return;
    }
    if (ein > emax) {
        /*Linear Extrapolation above tmax */
        CKCVBS(&tmax, y, iwrk, rwrk, &cv);
        *t = tmax - (emax-ein)/cv;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKUBMS(&t1,y,iwrk,rwrk,&e1);
        CKCVBS(&t1,y,iwrk,rwrk,&cv);
        dt = (ein - e1) / cv;
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
    }
    *t = t1;
    *ierr = 0;
    return;
}

/* End of file  */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENIMC EGTRANSETLENIMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENIMC egtransetlenimc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENIMC egtransetlenimc_
#endif
void egtransetLENIMC(int* LENIMC) {
  *LENIMC =           38;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENRMC EGTRANSETLENRMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENRMC egtransetlenrmc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENRMC egtransetlenrmc_
#endif
void egtransetLENRMC(int* LENRMC) {
  *LENRMC =         1854;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNO EGTRANSETNO
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNO egtransetno
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNO egtransetno_
#endif
void egtransetNO(int* NO) {
  *NO =            4;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetKK EGTRANSETKK
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetKK egtransetkk
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetKK egtransetkk_
#endif
void egtransetKK(int* KK) {
  *KK =            9;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLITE EGTRANSETNLITE
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLITE egtransetnlite
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLITE egtransetnlite_
#endif
void egtransetNLITE(int* NLITE) {
  *NLITE =            2;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetPATM EGTRANSETPATM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetPATM egtransetpatm
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetPATM egtransetpatm_
#endif
void egtransetPATM(double* PATM) {
  *PATM =   0.1013250000000000E+07;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetWT EGTRANSETWT
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetWT egtransetwt
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetWT egtransetwt_
#endif
void egtransetWT(double* WT) {
  WT[           0] =   0.2015939950942993E+01;
  WT[           1] =   0.3199880027770996E+02;
  WT[           2] =   0.1801534008979797E+02;
  WT[           3] =   0.1007969975471497E+01;
  WT[           4] =   0.1599940013885498E+02;
  WT[           5] =   0.1700737011432648E+02;
  WT[           6] =   0.3300677025318146E+02;
  WT[           7] =   0.3401474022865295E+02;
  WT[           8] =   0.2801339912414551E+02;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetEPS EGTRANSETEPS
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetEPS egtranseteps
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetEPS egtranseteps_
#endif
void egtransetEPS(double* EPS) {
  EPS[           0] =   0.3800000000000000E+02;
  EPS[           1] =   0.1074000000000000E+03;
  EPS[           2] =   0.5724000000000000E+03;
  EPS[           3] =   0.1450000000000000E+03;
  EPS[           4] =   0.8000000000000000E+02;
  EPS[           5] =   0.8000000000000000E+02;
  EPS[           6] =   0.1074000000000000E+03;
  EPS[           7] =   0.1074000000000000E+03;
  EPS[           8] =   0.9753000000000000E+02;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetSIG EGTRANSETSIG
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetSIG egtransetsig
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetSIG egtransetsig_
#endif
void egtransetSIG(double* SIG) {
  SIG[           0] =   0.2920000000000000E+01;
  SIG[           1] =   0.3458000000000000E+01;
  SIG[           2] =   0.2605000000000000E+01;
  SIG[           3] =   0.2050000000000000E+01;
  SIG[           4] =   0.2750000000000000E+01;
  SIG[           5] =   0.2750000000000000E+01;
  SIG[           6] =   0.3458000000000000E+01;
  SIG[           7] =   0.3458000000000000E+01;
  SIG[           8] =   0.3621000000000000E+01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetDIP EGTRANSETDIP
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetDIP egtransetdip
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetDIP egtransetdip_
#endif
void egtransetDIP(double* DIP) {
  DIP[           0] =   0.0000000000000000E+00;
  DIP[           1] =   0.0000000000000000E+00;
  DIP[           2] =   0.1844000000000000E+01;
  DIP[           3] =   0.0000000000000000E+00;
  DIP[           4] =   0.0000000000000000E+00;
  DIP[           5] =   0.0000000000000000E+00;
  DIP[           6] =   0.0000000000000000E+00;
  DIP[           7] =   0.0000000000000000E+00;
  DIP[           8] =   0.0000000000000000E+00;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetPOL EGTRANSETPOL
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetPOL egtransetpol
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetPOL egtransetpol_
#endif
void egtransetPOL(double* POL) {
  POL[           0] =   0.7900000000000000E+00;
  POL[           1] =   0.1600000000000000E+01;
  POL[           2] =   0.0000000000000000E+00;
  POL[           3] =   0.0000000000000000E+00;
  POL[           4] =   0.0000000000000000E+00;
  POL[           5] =   0.0000000000000000E+00;
  POL[           6] =   0.0000000000000000E+00;
  POL[           7] =   0.0000000000000000E+00;
  POL[           8] =   0.1760000000000000E+01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetZROT EGTRANSETZROT
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetZROT egtransetzrot
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetZROT egtransetzrot_
#endif
void egtransetZROT(double* ZROT) {
  ZROT[           0] =   0.2800000000000000E+03;
  ZROT[           1] =   0.3800000000000000E+01;
  ZROT[           2] =   0.4000000000000000E+01;
  ZROT[           3] =   0.0000000000000000E+00;
  ZROT[           4] =   0.0000000000000000E+00;
  ZROT[           5] =   0.0000000000000000E+00;
  ZROT[           6] =   0.1000000000000000E+01;
  ZROT[           7] =   0.3800000000000000E+01;
  ZROT[           8] =   0.4000000000000000E+01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLIN EGTRANSETNLIN
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLIN egtransetnlin
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLIN egtransetnlin_
#endif
void egtransetNLIN(int* NLIN) {
  NLIN[           0] =            1;
  NLIN[           1] =            1;
  NLIN[           2] =            2;
  NLIN[           3] =            0;
  NLIN[           4] =            0;
  NLIN[           5] =            1;
  NLIN[           6] =            2;
  NLIN[           7] =            2;
  NLIN[           8] =            1;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFLAM EGTRANSETCOFLAM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFLAM egtransetcoflam
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFLAM egtransetcoflam_
#endif
void egtransetCOFLAM(double* COFLAM) {
  COFLAM[           0] =   0.1109466217088744E+02;
  COFLAM[           1] =  -0.1314982879973400E+01;
  COFLAM[           2] =   0.2434842345478272E+00;
  COFLAM[           3] =  -0.8971821506339107E-02;
  COFLAM[           4] =  -0.2513896927262683E+01;
  COFLAM[           5] =   0.3151661994436593E+01;
  COFLAM[           6] =  -0.3099623183029560E+00;
  COFLAM[           7] =   0.1344793787858907E-01;
  COFLAM[           8] =   0.2212500250653218E+02;
  COFLAM[           9] =  -0.8451477268785496E+01;
  COFLAM[          10] =   0.1459333704606090E+01;
  COFLAM[          11] =  -0.7286069710088973E-01;
  COFLAM[          12] =  -0.3253788579319486E+00;
  COFLAM[          13] =   0.3416397862210875E+01;
  COFLAM[          14] =  -0.3631104489846262E+00;
  COFLAM[          15] =   0.1585986202702918E-01;
  COFLAM[          16] =   0.1969666526185763E+01;
  COFLAM[          17] =   0.1801396425547300E+01;
  COFLAM[          18] =  -0.1549102264837655E+00;
  COFLAM[          19] =   0.6908759721669425E-02;
  COFLAM[          20] =   0.1605449323860262E+02;
  COFLAM[          21] =  -0.4103244048729492E+01;
  COFLAM[          22] =   0.6631513498156665E+00;
  COFLAM[          23] =  -0.2977137989389932E-01;
  COFLAM[          24] =   0.5546401577805573E+00;
  COFLAM[          25] =   0.1591057931808813E+01;
  COFLAM[          26] =  -0.5282455808284543E-01;
  COFLAM[          27] =   0.4072391521895438E-03;
  COFLAM[          28] =   0.1486267947622502E+01;
  COFLAM[          29] =   0.1062271348472054E+01;
  COFLAM[          30] =   0.5716895740566978E-01;
  COFLAM[          31] =  -0.6382577562940780E-02;
  COFLAM[          32] =   0.1154289596065436E+02;
  COFLAM[          33] =  -0.2911524559760463E+01;
  COFLAM[          34] =   0.5546581827235442E+00;
  COFLAM[          35] =  -0.2750103005663747E-01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFETA EGTRANSETCOFETA
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFETA egtransetcofeta
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFETA egtransetcofeta_
#endif
void egtransetCOFETA(double* COFETA) {
  COFETA[           0] =  -0.1376710086380533E+02;
  COFETA[           1] =   0.9708665581008417E+00;
  COFETA[           2] =  -0.4534923959308444E-01;
  COFETA[           3] =   0.2096011921433090E-02;
  COFETA[           4] =  -0.1681080110872507E+02;
  COFETA[           5] =   0.2522528725237932E+01;
  COFETA[           6] =  -0.2490712798240462E+00;
  COFETA[           7] =   0.1100615806447165E-01;
  COFETA[           8] =  -0.1187800764707307E+02;
  COFETA[           9] =  -0.7882519505482808E+00;
  COFETA[          10] =   0.3341408170058896E+00;
  COFETA[          11] =  -0.1986366361647418E-01;
  COFETA[          12] =  -0.1987529414716897E+02;
  COFETA[          13] =   0.3416397862210910E+01;
  COFETA[          14] =  -0.3631104489846305E+00;
  COFETA[          15] =   0.1585986202702936E-01;
  COFETA[          16] =  -0.1481563591580279E+02;
  COFETA[          17] =   0.1801396425547327E+01;
  COFETA[          18] =  -0.1549102264837694E+00;
  COFETA[          19] =   0.6908759721669611E-02;
  COFETA[          20] =  -0.1478508813778873E+02;
  COFETA[          21] =   0.1801396425547305E+01;
  COFETA[          22] =  -0.1549102264837663E+00;
  COFETA[          23] =   0.6908759721669463E-02;
  COFETA[          24] =  -0.1679529396430702E+02;
  COFETA[          25] =   0.2522528725237947E+01;
  COFETA[          26] =  -0.2490712798240481E+00;
  COFETA[          27] =   0.1100615806447173E-01;
  COFETA[          28] =  -0.1678025333071109E+02;
  COFETA[          29] =   0.2522528725237946E+01;
  COFETA[          30] =  -0.2490712798240483E+00;
  COFETA[          31] =   0.1100615806447175E-01;
  COFETA[          32] =  -0.1626172843728175E+02;
  COFETA[          33] =   0.2251740453876138E+01;
  COFETA[          34] =  -0.2138340893699383E+00;
  COFETA[          35] =   0.9477823154448534E-02;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFD EGTRANSETCOFD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFD egtransetcofd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFD egtransetcofd_
#endif
void egtransetCOFD(double* COFD) {
  COFD[           0] =  -0.1023073719095340E+02;
  COFD[           1] =   0.2153598714195162E+01;
  COFD[           2] =  -0.6969019063693850E-01;
  COFD[           3] =   0.3233961733522977E-02;
  COFD[           4] =  -0.1228690337278460E+02;
  COFD[           5] =   0.2739817607547901E+01;
  COFD[           6] =  -0.1455907617836413E+00;
  COFD[           7] =   0.6496697467578658E-02;
  COFD[           8] =  -0.1787213075631945E+02;
  COFD[           9] =   0.4905112324365878E+01;
  COFD[          10] =  -0.4173872143518051E+00;
  COFD[          11] =   0.1788845489573668E-01;
  COFD[          12] =  -0.1168685516750758E+02;
  COFD[          13] =   0.2883726644585007E+01;
  COFD[          14] =  -0.1637777890869061E+00;
  COFD[          15] =   0.7265872474830623E-02;
  COFD[          16] =  -0.1060109579664357E+02;
  COFD[          17] =   0.2157129122185176E+01;
  COFD[          18] =  -0.6524729611183182E-01;
  COFD[          19] =   0.2809597810849289E-02;
  COFD[          20] =  -0.1060412419807472E+02;
  COFD[          21] =   0.2157123706881953E+01;
  COFD[          22] =  -0.6524690383610432E-01;
  COFD[          23] =   0.2809593789707914E-02;
  COFD[          24] =  -0.1229066279603360E+02;
  COFD[          25] =   0.2741057082889793E+01;
  COFD[          26] =  -0.1457627022193994E+00;
  COFD[          27] =   0.6504615244946633E-02;
  COFD[          28] =  -0.1229422250798061E+02;
  COFD[          29] =   0.2742232459895579E+01;
  COFD[          30] =  -0.1459257489141948E+00;
  COFD[          31] =   0.6512123400236565E-02;
  COFD[          32] =  -0.1217511489930491E+02;
  COFD[          33] =   0.2679286720576265E+01;
  COFD[          34] =  -0.1373997446113497E+00;
  COFD[          35] =   0.6125379451516252E-02;
  COFD[          36] =  -0.1228690337278460E+02;
  COFD[          37] =   0.2739817607547901E+01;
  COFD[          38] =  -0.1455907617836413E+00;
  COFD[          39] =   0.6496697467578658E-02;
  COFD[          40] =  -0.1579169675646239E+02;
  COFD[          41] =   0.3572143437285479E+01;
  COFD[          42] =  -0.2518469828462104E+00;
  COFD[          43] =   0.1102533331592793E-01;
  COFD[          44] =  -0.2036133619470493E+02;
  COFD[          45] =   0.5195864695910879E+01;
  COFD[          46] =  -0.4301216528920454E+00;
  COFD[          47] =   0.1744936825492251E-01;
  COFD[          48] =  -0.1719486110604737E+02;
  COFD[          49] =   0.4862600572193149E+01;
  COFD[          50] =  -0.4215909171603991E+00;
  COFD[          51] =   0.1846373197001472E-01;
  COFD[          52] =  -0.1474363832488901E+02;
  COFD[          53] =   0.3350090901251298E+01;
  COFD[          54] =  -0.2245382243795846E+00;
  COFD[          55] =   0.9906528312288899E-02;
  COFD[          56] =  -0.1473449002110577E+02;
  COFD[          57] =   0.3337793989426398E+01;
  COFD[          58] =  -0.2228575541864579E+00;
  COFD[          59] =   0.9830229599517820E-02;
  COFD[          60] =  -0.1579979030842146E+02;
  COFD[          61] =   0.3572309323030401E+01;
  COFD[          62] =  -0.2518694694768392E+00;
  COFD[          63] =   0.1102634618303224E-01;
  COFD[          64] =  -0.1580828869487550E+02;
  COFD[          65] =   0.3572786632028933E+01;
  COFD[          66] =  -0.2519341709914386E+00;
  COFD[          67] =   0.1102926053643803E-01;
  COFD[          68] =  -0.1544409203507008E+02;
  COFD[          69] =   0.3434913447661099E+01;
  COFD[          70] =  -0.2339977102148624E+00;
  COFD[          71] =   0.1025033359000777E-01;
  COFD[          72] =  -0.1787213075631945E+02;
  COFD[          73] =   0.4905112324365878E+01;
  COFD[          74] =  -0.4173872143518051E+00;
  COFD[          75] =   0.1788845489573668E-01;
  COFD[          76] =  -0.2036133619470493E+02;
  COFD[          77] =   0.5195864695910879E+01;
  COFD[          78] =  -0.4301216528920454E+00;
  COFD[          79] =   0.1744936825492251E-01;
  COFD[          80] =  -0.1301206458003408E+02;
  COFD[          81] =   0.1429168452568504E+01;
  COFD[          82] =   0.1661557715508861E+00;
  COFD[          83] =  -0.1214321823404827E-01;
  COFD[          84] =  -0.1695183404426207E+02;
  COFD[          85] =   0.4416203740445888E+01;
  COFD[          86] =  -0.3114890192543383E+00;
  COFD[          87] =   0.1159663710619628E-01;
  COFD[          88] =  -0.1894546586601941E+02;
  COFD[          89] =   0.4935720676293722E+01;
  COFD[          90] =  -0.4114057447467022E+00;
  COFD[          91] =   0.1723315986098057E-01;
  COFD[          92] =  -0.1896004007077217E+02;
  COFD[          93] =   0.4935377774711531E+01;
  COFD[          94] =  -0.4113926173034602E+00;
  COFD[          95] =   0.1723402811988101E-01;
  COFD[          96] =  -0.1973436691750415E+02;
  COFD[          97] =   0.4993125184848788E+01;
  COFD[          98] =  -0.4088531920998837E+00;
  COFD[          99] =   0.1672325900844039E-01;
  COFD[         100] =  -0.1972379377029397E+02;
  COFD[         101] =   0.4985700637038792E+01;
  COFD[         102] =  -0.4077156220815392E+00;
  COFD[         103] =   0.1666668649763390E-01;
  COFD[         104] =  -0.2025975101746627E+02;
  COFD[         105] =   0.5176212790334279E+01;
  COFD[         106] =  -0.4308686680573706E+00;
  COFD[         107] =   0.1761066268522406E-01;
  COFD[         108] =  -0.1168685516750758E+02;
  COFD[         109] =   0.2883726644585007E+01;
  COFD[         110] =  -0.1637777890869061E+00;
  COFD[         111] =   0.7265872474830623E-02;
  COFD[         112] =  -0.1719486110604737E+02;
  COFD[         113] =   0.4862600572193149E+01;
  COFD[         114] =  -0.4215909171603991E+00;
  COFD[         115] =   0.1846373197001472E-01;
  COFD[         116] =  -0.1695183404426207E+02;
  COFD[         117] =   0.4416203740445888E+01;
  COFD[         118] =  -0.3114890192543383E+00;
  COFD[         119] =   0.1159663710619628E-01;
  COFD[         120] =  -0.1476537827842222E+02;
  COFD[         121] =   0.4195448036781197E+01;
  COFD[         122] =  -0.3279615676730872E+00;
  COFD[         123] =   0.1412388797372668E-01;
  COFD[         124] =  -0.1500301598770682E+02;
  COFD[         125] =   0.4131920691035950E+01;
  COFD[         126] =  -0.3275337549968449E+00;
  COFD[         127] =   0.1442760802547180E-01;
  COFD[         128] =  -0.1502025959280452E+02;
  COFD[         129] =   0.4138327950805945E+01;
  COFD[         130] =  -0.3284006010219734E+00;
  COFD[         131] =   0.1446660099762999E-01;
  COFD[         132] =  -0.1720445512239349E+02;
  COFD[         133] =   0.4866330968714804E+01;
  COFD[         134] =  -0.4220901924876364E+00;
  COFD[         135] =   0.1848595807162496E-01;
  COFD[         136] =  -0.1721362788263434E+02;
  COFD[         137] =   0.4869900420159440E+01;
  COFD[         138] =  -0.4225679381202458E+00;
  COFD[         139] =   0.1850722613298993E-01;
  COFD[         140] =  -0.1676897875889624E+02;
  COFD[         141] =   0.4677755283143664E+01;
  COFD[         142] =  -0.3974479445519867E+00;
  COFD[         143] =   0.1741111717763365E-01;
  COFD[         144] =  -0.1060109579664357E+02;
  COFD[         145] =   0.2157129122185176E+01;
  COFD[         146] =  -0.6524729611183182E-01;
  COFD[         147] =   0.2809597810849289E-02;
  COFD[         148] =  -0.1474363832488901E+02;
  COFD[         149] =   0.3350090901251298E+01;
  COFD[         150] =  -0.2245382243795846E+00;
  COFD[         151] =   0.9906528312288899E-02;
  COFD[         152] =  -0.1894546586601941E+02;
  COFD[         153] =   0.4935720676293722E+01;
  COFD[         154] =  -0.4114057447467022E+00;
  COFD[         155] =   0.1723315986098057E-01;
  COFD[         156] =  -0.1500301598770682E+02;
  COFD[         157] =   0.4131920691035950E+01;
  COFD[         158] =  -0.3275337549968449E+00;
  COFD[         159] =   0.1442760802547180E-01;
  COFD[         160] =  -0.1329503778398344E+02;
  COFD[         161] =   0.2938989816956702E+01;
  COFD[         162] =  -0.1706012465018490E+00;
  COFD[         163] =   0.7547358484820122E-02;
  COFD[         164] =  -0.1331106228772027E+02;
  COFD[         165] =   0.2939408772392007E+01;
  COFD[         166] =  -0.1706589698617190E+00;
  COFD[         167] =   0.7549998035548285E-02;
  COFD[         168] =  -0.1476396639560408E+02;
  COFD[         169] =   0.3356474034011327E+01;
  COFD[         170] =  -0.2254105976116209E+00;
  COFD[         171] =   0.9946130954817655E-02;
  COFD[         172] =  -0.1478376315874751E+02;
  COFD[         173] =   0.3362741024796286E+01;
  COFD[         174] =  -0.2262670743508222E+00;
  COFD[         175] =   0.9985011202484037E-02;
  COFD[         176] =  -0.1442662508885140E+02;
  COFD[         177] =   0.3216422191096490E+01;
  COFD[         178] =  -0.2069561315576079E+00;
  COFD[         179] =   0.9135288598631700E-02;
  COFD[         180] =  -0.1060412419807472E+02;
  COFD[         181] =   0.2157123706881953E+01;
  COFD[         182] =  -0.6524690383610432E-01;
  COFD[         183] =   0.2809593789707914E-02;
  COFD[         184] =  -0.1473449002110577E+02;
  COFD[         185] =   0.3337793989426398E+01;
  COFD[         186] =  -0.2228575541864579E+00;
  COFD[         187] =   0.9830229599517820E-02;
  COFD[         188] =  -0.1896004007077217E+02;
  COFD[         189] =   0.4935377774711531E+01;
  COFD[         190] =  -0.4113926173034602E+00;
  COFD[         191] =   0.1723402811988101E-01;
  COFD[         192] =  -0.1502025959280452E+02;
  COFD[         193] =   0.4138327950805945E+01;
  COFD[         194] =  -0.3284006010219734E+00;
  COFD[         195] =   0.1446660099762999E-01;
  COFD[         196] =  -0.1331106228772027E+02;
  COFD[         197] =   0.2939408772392007E+01;
  COFD[         198] =  -0.1706589698617190E+00;
  COFD[         199] =   0.7549998035548285E-02;
  COFD[         200] =  -0.1332558556199739E+02;
  COFD[         201] =   0.2938989816956671E+01;
  COFD[         202] =  -0.1706012465018444E+00;
  COFD[         203] =   0.7547358484819904E-02;
  COFD[         204] =  -0.1475457475365462E+02;
  COFD[         205] =   0.3343986593397501E+01;
  COFD[         206] =  -0.2237039347662408E+00;
  COFD[         207] =   0.9868653787965939E-02;
  COFD[         208] =  -0.1477418610290300E+02;
  COFD[         209] =   0.3350090901251285E+01;
  COFD[         210] =  -0.2245382243795827E+00;
  COFD[         211] =   0.9906528312288795E-02;
  COFD[         212] =  -0.1442063092189782E+02;
  COFD[         213] =   0.3205844724257626E+01;
  COFD[         214] =  -0.2055144159347264E+00;
  COFD[         215] =   0.9070008933822527E-02;
  COFD[         216] =  -0.1229066279603360E+02;
  COFD[         217] =   0.2741057082889793E+01;
  COFD[         218] =  -0.1457627022193994E+00;
  COFD[         219] =   0.6504615244946633E-02;
  COFD[         220] =  -0.1579979030842146E+02;
  COFD[         221] =   0.3572309323030401E+01;
  COFD[         222] =  -0.2518694694768392E+00;
  COFD[         223] =   0.1102634618303224E-01;
  COFD[         224] =  -0.1973436691750415E+02;
  COFD[         225] =   0.4993125184848788E+01;
  COFD[         226] =  -0.4088531920998837E+00;
  COFD[         227] =   0.1672325900844039E-01;
  COFD[         228] =  -0.1720445512239349E+02;
  COFD[         229] =   0.4866330968714804E+01;
  COFD[         230] =  -0.4220901924876364E+00;
  COFD[         231] =   0.1848595807162496E-01;
  COFD[         232] =  -0.1476396639560408E+02;
  COFD[         233] =   0.3356474034011327E+01;
  COFD[         234] =  -0.2254105976116209E+00;
  COFD[         235] =   0.9946130954817655E-02;
  COFD[         236] =  -0.1475457475365462E+02;
  COFD[         237] =   0.3343986593397501E+01;
  COFD[         238] =  -0.2237039347662408E+00;
  COFD[         239] =   0.9868653787965939E-02;
  COFD[         240] =  -0.1580720390088048E+02;
  COFD[         241] =   0.3572143437285478E+01;
  COFD[         242] =  -0.2518469828462103E+00;
  COFD[         243] =   0.1102533331592793E-01;
  COFD[         244] =  -0.1581504405580379E+02;
  COFD[         245] =   0.3572299494956542E+01;
  COFD[         246] =  -0.2518681372333372E+00;
  COFD[         247] =   0.1102628617469197E-01;
  COFD[         248] =  -0.1545394137084459E+02;
  COFD[         249] =   0.3436021618273004E+01;
  COFD[         250] =  -0.2341477712817018E+00;
  COFD[         251] =   0.1025708506984392E-01;
  COFD[         252] =  -0.1229422250798061E+02;
  COFD[         253] =   0.2742232459895579E+01;
  COFD[         254] =  -0.1459257489141948E+00;
  COFD[         255] =   0.6512123400236565E-02;
  COFD[         256] =  -0.1580828869487550E+02;
  COFD[         257] =   0.3572786632028933E+01;
  COFD[         258] =  -0.2519341709914386E+00;
  COFD[         259] =   0.1102926053643803E-01;
  COFD[         260] =  -0.1972379377029397E+02;
  COFD[         261] =   0.4985700637038792E+01;
  COFD[         262] =  -0.4077156220815392E+00;
  COFD[         263] =   0.1666668649763390E-01;
  COFD[         264] =  -0.1721362788263434E+02;
  COFD[         265] =   0.4869900420159440E+01;
  COFD[         266] =  -0.4225679381202458E+00;
  COFD[         267] =   0.1850722613298993E-01;
  COFD[         268] =  -0.1478376315874751E+02;
  COFD[         269] =   0.3362741024796286E+01;
  COFD[         270] =  -0.2262670743508222E+00;
  COFD[         271] =   0.9985011202484037E-02;
  COFD[         272] =  -0.1477418610290300E+02;
  COFD[         273] =   0.3350090901251285E+01;
  COFD[         274] =  -0.2245382243795827E+00;
  COFD[         275] =   0.9906528312288795E-02;
  COFD[         276] =  -0.1581504405580379E+02;
  COFD[         277] =   0.3572299494956542E+01;
  COFD[         278] =  -0.2518681372333372E+00;
  COFD[         279] =   0.1102628617469197E-01;
  COFD[         280] =  -0.1582224453447645E+02;
  COFD[         281] =   0.3572143437285498E+01;
  COFD[         282] =  -0.2518469828462133E+00;
  COFD[         283] =   0.1102533331592807E-01;
  COFD[         284] =  -0.1546403419937241E+02;
  COFD[         285] =   0.3437367417375052E+01;
  COFD[         286] =  -0.2343299630280416E+00;
  COFD[         287] =   0.1026528034463288E-01;
  COFD[         288] =  -0.1217511489930491E+02;
  COFD[         289] =   0.2679286720576265E+01;
  COFD[         290] =  -0.1373997446113497E+00;
  COFD[         291] =   0.6125379451516252E-02;
  COFD[         292] =  -0.1544409203507008E+02;
  COFD[         293] =   0.3434913447661099E+01;
  COFD[         294] =  -0.2339977102148624E+00;
  COFD[         295] =   0.1025033359000777E-01;
  COFD[         296] =  -0.2025975101746627E+02;
  COFD[         297] =   0.5176212790334279E+01;
  COFD[         298] =  -0.4308686680573706E+00;
  COFD[         299] =   0.1761066268522406E-01;
  COFD[         300] =  -0.1676897875889624E+02;
  COFD[         301] =   0.4677755283143664E+01;
  COFD[         302] =  -0.3974479445519867E+00;
  COFD[         303] =   0.1741111717763365E-01;
  COFD[         304] =  -0.1442662508885140E+02;
  COFD[         305] =   0.3216422191096490E+01;
  COFD[         306] =  -0.2069561315576079E+00;
  COFD[         307] =   0.9135288598631700E-02;
  COFD[         308] =  -0.1442063092189782E+02;
  COFD[         309] =   0.3205844724257626E+01;
  COFD[         310] =  -0.2055144159347264E+00;
  COFD[         311] =   0.9070008933822527E-02;
  COFD[         312] =  -0.1545394137084459E+02;
  COFD[         313] =   0.3436021618273004E+01;
  COFD[         314] =  -0.2341477712817018E+00;
  COFD[         315] =   0.1025708506984392E-01;
  COFD[         316] =  -0.1546403419937241E+02;
  COFD[         317] =   0.3437367417375052E+01;
  COFD[         318] =  -0.2343299630280416E+00;
  COFD[         319] =   0.1026528034463288E-01;
  COFD[         320] =  -0.1521415383275665E+02;
  COFD[         321] =   0.3348053449783496E+01;
  COFD[         322] =  -0.2233657260417595E+00;
  COFD[         323] =   0.9817787279109837E-02;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetKTDIF EGTRANSETKTDIF
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetKTDIF egtransetktdif
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetKTDIF egtransetktdif_
#endif
void egtransetKTDIF(int* KTDIF) {
  KTDIF[           0] =            1;
  KTDIF[           1] =            4;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFTD EGTRANSETCOFTD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFTD egtransetcoftd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFTD egtransetcoftd_
#endif
void egtransetCOFTD(double* COFTD) {
  COFTD[           0] =   0.0000000000000000E+00;
  COFTD[           1] =   0.0000000000000000E+00;
  COFTD[           2] =   0.0000000000000000E+00;
  COFTD[           3] =   0.0000000000000000E+00;
  COFTD[           4] =   0.4427392097103151E+00;
  COFTD[           5] =   0.7117708739452678E-04;
  COFTD[           6] =  -0.3847681435586731E-07;
  COFTD[           7] =   0.6863235673780089E-11;
  COFTD[           8] =   0.6020321200392946E-01;
  COFTD[           9] =   0.5615615962686471E-03;
  COFTD[          10] =  -0.2553727790421687E-06;
  COFTD[          11] =   0.3633898304582181E-10;
  COFTD[          12] =  -0.1525347878079855E+00;
  COFTD[          13] =  -0.5464040665780899E-04;
  COFTD[          14] =   0.2934125089337968E-07;
  COFTD[          15] =  -0.4870919731196246E-11;
  COFTD[          16] =   0.4155834523029710E+00;
  COFTD[          17] =   0.1097383911871691E-04;
  COFTD[          18] =  -0.3960224834087406E-08;
  COFTD[          19] =   0.1144145307552153E-11;
  COFTD[          20] =   0.4219325599835782E+00;
  COFTD[          21] =   0.1114149277732201E-04;
  COFTD[          22] =  -0.4020727469049531E-08;
  COFTD[          23] =   0.1161625074178193E-11;
  COFTD[          24] =   0.4444526949057779E+00;
  COFTD[          25] =   0.7145255629999482E-04;
  COFTD[          26] =  -0.3862572696699677E-07;
  COFTD[          27] =   0.6889797705021213E-11;
  COFTD[          28] =   0.4460703094919803E+00;
  COFTD[          29] =   0.7171261254133862E-04;
  COFTD[          30] =  -0.3876630782084383E-07;
  COFTD[          31] =   0.6914873573367535E-11;
  COFTD[          32] =   0.4452620885099525E+00;
  COFTD[          33] =   0.4946972054970773E-04;
  COFTD[          34] =  -0.2630235132337235E-07;
  COFTD[          35] =   0.4903063330400376E-11;
  COFTD[          36] =   0.1525347878079855E+00;
  COFTD[          37] =   0.5464040665780899E-04;
  COFTD[          38] =  -0.2934125089337968E-07;
  COFTD[          39] =   0.4870919731196246E-11;
  COFTD[          40] =   0.2204829568680314E+00;
  COFTD[          41] =   0.4801643237618346E-03;
  COFTD[          42] =  -0.2329279629943256E-06;
  COFTD[          43] =   0.3464704623498651E-10;
  COFTD[          44] =  -0.1418836567013430E+00;
  COFTD[          45] =   0.7665588601033144E-03;
  COFTD[          46] =  -0.3065500230258140E-06;
  COFTD[          47] =   0.4029595271283201E-10;
  COFTD[          48] =   0.0000000000000000E+00;
  COFTD[          49] =   0.0000000000000000E+00;
  COFTD[          50] =   0.0000000000000000E+00;
  COFTD[          51] =   0.0000000000000000E+00;
  COFTD[          52] =   0.2700102625474920E+00;
  COFTD[          53] =   0.3615551214404962E-03;
  COFTD[          54] =  -0.1807447676772846E-06;
  COFTD[          55] =   0.2753212716087208E-10;
  COFTD[          56] =   0.2720417770528547E+00;
  COFTD[          57] =   0.3642754049836650E-03;
  COFTD[          58] =  -0.1821046627192010E-06;
  COFTD[          59] =   0.2773927453061686E-10;
  COFTD[          60] =   0.2209079675460805E+00;
  COFTD[          61] =   0.4810899053474402E-03;
  COFTD[          62] =  -0.2333769631025201E-06;
  COFTD[          63] =   0.3471383309607501E-10;
  COFTD[          64] =   0.2213085142120510E+00;
  COFTD[          65] =   0.4819622095914184E-03;
  COFTD[          66] =  -0.2338001183446032E-06;
  COFTD[          67] =   0.3477677564298334E-10;
  COFTD[          68] =   0.2407445350007590E+00;
  COFTD[          69] =   0.4453434831756282E-03;
  COFTD[          70] =  -0.2181738909148746E-06;
  COFTD[          71] =   0.3269585309435789E-10;
};




#if 0




\\
\\
\\  This is the mechanism file
\\
\\
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
! H2/O2 oxidation reaction mechanism --
! (c) Li, Zhao, Kazakov, and Dryer, Princeton University, 2003.
!
!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! HOW TO USE THIS MECHANISM:
!
! Due to 
! (1) limitations of CHEMKIN-II format (specifically, an inability to implement
!     temperature-dependent collision efficiencies in falloff reactions)
! and
! (2) lack of fundamental understanding of the mixing rules for the falloff 
!     reactions with the bath gases that have different broadening factors,
!
! the present implementation represents a compromise (approximate) formulation.
!
! As a consequence, PRIOR TO ITS USE IN THE CALCULATIONS, THIS FILE HAS TO BE
! MODIFIED. DEPENDING ON WHAT BATH GAS (DILUTANT) IS MOST ABUNDANT IN YOUR SYSTEM
! (THE PRESENT CHOICES ARE N2, AR, OR HE),  YOU  SHOULD UNCOMMENT THE CORRESPONDING
! BLOCK FOR THE REACTION H+O2(+M)=HO2(+M), AND COMMENT THE BLOCK FOR OTHER DILUTANT(S).
! AS GIVEN, THE MAIN DILUTANT IS SET TO BE N2.
!
! 
! HOW TO REFERENCE THIS MECHANISM:
!
! Li, J., Zhao, Z., Kazakov, A., and Dryer, F.L. "An Updated Comprehensive Kinetic Model
! of Hydrogen Combustion", Int. J. Chem. Kinet. 2004 (in press).
!
!
! HOW TO CONTACT THE AUTHORS:
!
!    Prof. Frederick L. Dryer 
!    D-329-D Engineering Quadrangle 
!    Mechanical and Aerospace Engineering 
!    Princeton University 
!    Princeton, NJ 08544-5263 
!    Phone: 609-258-5206 
!    Lab:    609-258-0316 
!    FAX:    609-258-1939
!    Email: fldryer@Princeton.EDU
! 
!**********************************************************************************************
! Development notes:
!
!The following H2/O2 mechanism is based on Mueller et al's (Int.J.Chem.Kinet.1999,31:113)
!Changes:
!
!1.update the standard heat of formation of OH at 0K to 8.85kcal/mol (Ruscic et al, 
!  J. Phys. Chem. A, 2002, 106:2727)
!
!2.update the rate constant of H+O2=O+OH as proposed by Hessler (J. Phys. Chem. A, 1998,
!  102:4517)
!
!3.update the low-pressure-limit rate constant of H+O2(+M)=HO2(+M) with bath gases: H2, 
!  O2, N2, AR, HE, H2O as proposed by Michael et al (J. Phys. Chem. A, 2002,106:5297).
!  The third-body efficiency of H2, O2, and H2O are taken as the average value over 
!  the temperature range of 300-3000K. 
!  The Fc in Troe's form with N2 and AR/HE as bath gas are different, so the fall-off 
!  kinetics is expressed in two sets, for N2 and AR/HE, respectively.
! 
!4.for all other recombination reactions, assume the third-body efficiency of HE is
!  the same as AR.
!
!5.modify the A factor of the rate constant of H+OH+M=H2O+M to 3.8E+22.
!
!END OF NOTES
!**********************************************************************************************
ELEMENTS
H O N
END
!C AR HE

SPECIES
H2 O2 H2O H O OH HO2 H2O2 N2
END
! AR HE CO CO2

THERMO ALL
300.0 1000.0 5000.0
H                 120186H   1               G  0300.00   5000.00  1000.00      1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
 2.54716270E+04-4.60117638E-01 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 2.54716270E+04-4.60117608E-01                   4
O                 120186O   1               G  0300.00   5000.00  1000.00      1
 2.54205966E+00-2.75506191E-05-3.10280335E-09 4.55106742E-12-4.36805150E-16    2
 2.92308027E+04 4.92030811E+00 2.94642878E+00-1.63816649E-03 2.42103170E-06    3
-1.60284319E-09 3.89069636E-13 2.91476445E+04 2.96399498E+00                   4
OH                S 9/01O   1H   1    0    0G   200.000  6000.000 1000.        1
 2.86472886E+00 1.05650448E-03-2.59082758E-07 3.05218674E-11-1.33195876E-15    2
 3.68362875E+03 5.70164073E+00 4.12530561E+00-3.22544939E-03 6.52764691E-06    3
-5.79853643E-09 2.06237379E-12 3.34630913E+03-6.90432960E-01 4.51532273E+03    4
H2                121286H   2               G  0300.00   5000.00  1000.00      1
 2.99142337E+00 7.00064411E-04-5.63382869E-08-9.23157818E-12 1.58275179E-15    2
-8.35033997E+02-1.35511017E+00 3.29812431E+00 8.24944174E-04-8.14301529E-07    3
-9.47543433E-11 4.13487224E-13-1.01252087E+03-3.29409409E+00                   4
O2                121386O   2               G  0300.00   5000.00  1000.00      1
 3.69757819E+00 6.13519689E-04-1.25884199E-07 1.77528148E-11-1.13643531E-15    2
-1.23393018E+03 3.18916559E+00 3.21293640E+00 1.12748635E-03-5.75615047E-07    3
 1.31387723E-09-8.76855392E-13-1.00524902E+03 6.03473759E+00                   4
H2O                20387H   2O   1          G  0300.00   5000.00  1000.00      1
 2.67214561E+00 3.05629289E-03-8.73026011E-07 1.20099639E-10-6.39161787E-15    2
-2.98992090E+04 6.86281681E+00 3.38684249E+00 3.47498246E-03-6.35469633E-06    3
 6.96858127E-09-2.50658847E-12-3.02081133E+04 2.59023285E+00                   4
HO2               L 5/89H   1O   2   00   00G   200.000  3500.000  1000.000    1
 4.01721090E+00 2.23982013E-03-6.33658150E-07 1.14246370E-10-1.07908535E-14    2
 1.11856713E+02 3.78510215E+00 4.30179801E+00-4.74912051E-03 2.11582891E-05    3
-2.42763894E-08 9.29225124E-12 2.94808040E+02 3.71666245E+00 1.00021620E+04    4
H2O2              120186H   2O   2          G  0300.00   5000.00  1000.00      1
 4.57316685E+00 4.33613639E-03-1.47468882E-06 2.34890357E-10-1.43165356E-14    2
-1.80069609E+04 5.01136959E-01 3.38875365E+00 6.56922581E-03-1.48501258E-07    3
-4.62580552E-09 2.47151475E-12-1.76631465E+04 6.78536320E+00                   4
N2                121286N   2               G  0300.00   5000.00  1000.00      1
 0.02926640E+02 0.01487977E-01-0.05684761E-05 0.01009704E-08-0.06753351E-13    2
-0.09227977E+04 0.05980528E+02 0.03298677E+02 0.01408240E-01-0.03963222E-04    3
 0.05641515E-07-0.02444855E-10-0.01020900E+05 0.03950372E+02                   4
END

!AR                120186AR  1               G  0300.00   5000.00  1000.00      1
! 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
!-0.07453750E+04 0.04366001E+02 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
! 0.00000000E+00 0.00000000E+00-0.07453750E+04 0.04366001E+02                   4
!HE                120186HE  1               G  0300.00   5000.00  1000.00      1
! 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
!-0.07453750E+04 0.09153489E+01 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
! 0.00000000E+00 0.00000000E+00-0.07453750E+04 0.09153488E+01                   4
!CO                121286C   1O   1          G  0300.00   5000.00  1000.00      1
! 0.03025078E+02 0.01442689E-01-0.05630828E-05 0.01018581E-08-0.06910952E-13    2
!-0.01426835E+06 0.06108218E+02 0.03262452E+02 0.01511941E-01-0.03881755E-04    3
! 0.05581944E-07-0.02474951E-10-0.01431054E+06 0.04848897E+02                   4
!CO2               121286C   1O   2          G  0300.00   5000.00  1000.00      1
! 0.04453623E+02 0.03140169E-01-0.01278411E-04 0.02393997E-08-0.01669033E-12    2
!-0.04896696E+06-0.09553959E+01 0.02275725E+02 0.09922072E-01-0.01040911E-03    3
! 0.06866687E-07-0.02117280E-10-0.04837314E+06 0.01018849E+03                   4


REACTIONS

!H2-O2 Chain Reactions

! Hessler, J. Phys. Chem. A, 102:4517 (1998)
H+O2=O+OH                 3.547e+15 -0.406  1.6599E+4

! Sutherland et al., 21st Symposium, p. 929 (1986)
O+H2=H+OH                 0.508E+05  2.67  0.629E+04   

! Michael and Sutherland, J. Phys. Chem. 92:3853 (1988)
H2+OH=H2O+H               0.216E+09  1.51  0.343E+04

! Sutherland et al., 23rd Symposium, p. 51 (1990)
O+H2O=OH+OH               2.97e+06   2.02  1.34e+4


!H2-O2 Dissociation Reactions

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
H2+M=H+H+M                4.577E+19 -1.40  1.0438E+05
   H2/2.5/ H2O/12/
!   CO/1.9/ CO2/3.8/
!   AR/0.0/ HE/0.0/

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
!H2+AR=H+H+AR              5.84e18   -1.1   1.0438E+05 			

!H2+HE=H+H+HE              5.84e18   -1.1   1.0438E+05

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
O+O+M=O2+M                6.165E+15 -0.50  0.000E+00 
   H2/2.5/ H2O/12/
!   AR/0.0/  HE/0.0/
!   CO/1.9/ CO2/3.8/

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
!O+O+AR=O2+AR              1.886E+13 0.00  -1.788E+03

!O+O+HE=O2+HE              1.886E+13 0.00  -1.788E+03

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
O+H+M=OH+M                4.714E+18 -1.00  0.000E+00
   H2/2.5/ H2O/12/
!   AR/0.75/ HE/0.75/
!   CO/1.9/ CO2/3.8/

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
!H+OH+M=H2O+M              2.212E+22 -2.00  0.000E+00
H+OH+M=H2O+M               3.800E+22 -2.00  0.000E+00  
   H2/2.5/ H2O/12/   
!   AR/0.38/ HE/0.38/  
!   CO/1.9/ CO2/3.8/


!Formation and Consumption of HO2

! Cobos et al., J. Phys. Chem. 89:342 (1985) for kinf
! Michael, et al., J. Phys. Chem. A, 106:5297 (2002) for k0

!******************************************************************************
! MAIN BATH GAS IS N2 (comment this reaction otherwise)
!
H+O2(+M)=HO2(+M)      1.475E+12  0.60  0.00E+00
    LOW/6.366E+20  -1.72  5.248E+02/
    TROE/0.8  1E-30  1E+30/
    H2/2.0/ H2O/11./ O2/0.78/
!CO/1.9/ CO2/3.8/
    
!******************************************************************************
! MAIN BATH GAS IS AR OR HE (comment this reaction otherwise)
!
!H+O2(+M)=HO2(+M)      1.475E+12  0.60  0.00E+00
!    LOW/9.042E+19  -1.50  4.922E+02/
!    TROE/0.5 1E-30  1E+30/
!    H2/3.0/ H2O/16/ O2/1.1/ CO/2.7/ CO2/5.4/ HE/1.2/

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986) [modified]
HO2+H=H2+O2               1.66E+13   0.00   0.823E+03

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986) [modified]
HO2+H=OH+OH               7.079E+13   0.00   2.95E+02

! Baulch et al., J. Phys. Chem. Ref Data, 21:411 (1992)
HO2+O=O2+OH               0.325E+14  0.00   0.00E+00   

! Keyser, J. Phys. Chem. 92:1193 (1988)
HO2+OH=H2O+O2             2.890E+13  0.00 -4.970E+02


!Formation and Consumption of H2O2

! Hippler et al., J. Chem. Phys. 93:1755 (1990)
HO2+HO2=H2O2+O2            4.200e+14  0.00  1.1982e+04
  DUPLICATE
HO2+HO2=H2O2+O2            1.300e+11  0.00 -1.6293e+3
  DUPLICATE

! Brouwer et al., J. Chem. Phys. 86:6171 (1987) for kinf
! Warnatz, J. in Combustion chemistry (1984) for k0
H2O2(+M)=OH+OH(+M)         2.951e+14   0.00  4.843E+04 
  LOW/1.202E+17  0.00  4.55E+04/
  TROE/0.5 1E-30 1E+30/
  H2/2.5/ H2O/12/          
!  CO/1.9/ CO2/3.8/
!  AR/0.64/ HE/0.64/

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
H2O2+H=H2O+OH             0.241E+14  0.00  0.397E+04

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
H2O2+H=HO2+H2             0.482E+14  0.00  0.795E+04  

! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
H2O2+O=OH+HO2 		  9.550E+06  2.00  3.970E+03 

! Hippler and Troe, J. Chem. Phys. Lett. 192:333 (1992)
H2O2+OH=HO2+H2O           1.000E+12  0.00  0.000	
    DUPLICATE 
H2O2+OH=HO2+H2O           5.800E+14  0.00  9.557E+03
    DUPLICATE

END

\\
\\
\\  This is the therm file
\\
\\
THERMO ALL
300.0 1000.0 5000.0
H                 120186H   1               G  0300.00   5000.00  1000.00      1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
 2.54716270E+04-4.60117638E-01 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 2.54716270E+04-4.60117608E-01                   4
O                 120186O   1               G  0300.00   5000.00  1000.00      1
 2.54205966E+00-2.75506191E-05-3.10280335E-09 4.55106742E-12-4.36805150E-16    2
 2.92308027E+04 4.92030811E+00 2.94642878E+00-1.63816649E-03 2.42103170E-06    3
-1.60284319E-09 3.89069636E-13 2.91476445E+04 2.96399498E+00                   4
OH                S 9/01O   1H   1    0    0G   200.000  6000.000 1000.        1
 2.86472886E+00 1.05650448E-03-2.59082758E-07 3.05218674E-11-1.33195876E-15    2
 3.68362875E+03 5.70164073E+00 4.12530561E+00-3.22544939E-03 6.52764691E-06    3
-5.79853643E-09 2.06237379E-12 3.34630913E+03-6.90432960E-01 4.51532273E+03    4
H2                121286H   2               G  0300.00   5000.00  1000.00      1
 2.99142337E+00 7.00064411E-04-5.63382869E-08-9.23157818E-12 1.58275179E-15    2
-8.35033997E+02-1.35511017E+00 3.29812431E+00 8.24944174E-04-8.14301529E-07    3
-9.47543433E-11 4.13487224E-13-1.01252087E+03-3.29409409E+00                   4
O2                121386O   2               G  0300.00   5000.00  1000.00      1
 3.69757819E+00 6.13519689E-04-1.25884199E-07 1.77528148E-11-1.13643531E-15    2
-1.23393018E+03 3.18916559E+00 3.21293640E+00 1.12748635E-03-5.75615047E-07    3
 1.31387723E-09-8.76855392E-13-1.00524902E+03 6.03473759E+00                   4
H2O                20387H   2O   1          G  0300.00   5000.00  1000.00      1
 2.67214561E+00 3.05629289E-03-8.73026011E-07 1.20099639E-10-6.39161787E-15    2
-2.98992090E+04 6.86281681E+00 3.38684249E+00 3.47498246E-03-6.35469633E-06    3
 6.96858127E-09-2.50658847E-12-3.02081133E+04 2.59023285E+00                   4
HO2               L 5/89H   1O   2   00   00G   200.000  3500.000  1000.000    1
 4.01721090E+00 2.23982013E-03-6.33658150E-07 1.14246370E-10-1.07908535E-14    2
 1.11856713E+02 3.78510215E+00 4.30179801E+00-4.74912051E-03 2.11582891E-05    3
-2.42763894E-08 9.29225124E-12 2.94808040E+02 3.71666245E+00 1.00021620E+04    4
H2O2              120186H   2O   2          G  0300.00   5000.00  1000.00      1
 4.57316685E+00 4.33613639E-03-1.47468882E-06 2.34890357E-10-1.43165356E-14    2
-1.80069609E+04 5.01136959E-01 3.38875365E+00 6.56922581E-03-1.48501258E-07    3
-4.62580552E-09 2.47151475E-12-1.76631465E+04 6.78536320E+00                   4
N2                121286N   2               G  0300.00   5000.00  1000.00      1
 0.02926640E+02 0.01487977E-01-0.05684761E-05 0.01009704E-08-0.06753351E-13    2
-0.09227977E+04 0.05980528E+02 0.03298677E+02 0.01408240E-01-0.03963222E-04    3
 0.05641515E-07-0.02444855E-10-0.01020900E+05 0.03950372E+02                   4
END

\\
\\
\\  This is the tran file
\\
\\
                                                                                
AR                 0   136.500     3.330     0.000     0.000     0.000          
C                  0    71.400     3.298     0.000     0.000     0.000 ! *      
CH                 1    80.000     2.750     0.000     0.000     0.000          
CH2                1   144.000     3.800     0.000     0.000     0.000          
CH2*               1   144.000     3.800     0.000     0.000     0.000          
CH3                1   144.000     3.800     0.000     0.000     0.000          
CH4                2   141.400     3.746     0.000     2.600    13.000          
CO                 1    98.100     3.650     0.000     1.950     1.800          
CO2                1   244.000     3.763     0.000     2.650     2.100          
HCO                2   498.000     3.590     0.000     0.000     0.000          
CH2O               2   498.000     3.590     0.000     0.000     2.000          
CH2OH              2   417.000     3.690     1.700     0.000     2.000          
CH3O               2   417.000     3.690     1.700     0.000     2.000          
CH3OH              2   481.800     3.626     0.000     0.000     1.000 ! SVE    
                                                                                
C2                 1    97.530     3.621     0.000     1.760     4.000          
C2O                1   232.400     3.828     0.000     0.000     1.000 ! *      
C2H                1   209.000     4.100     0.000     0.000     2.500          
C2H2               1   209.000     4.100     0.000     0.000     2.500          
H2CC               2   209.000     4.100     0.000     0.000     2.500
C2H3               2   209.000     4.100     0.000     0.000     1.000 ! *      
C2H4               2   280.800     3.971     0.000     0.000     1.500          
C2H5               2   252.300     4.302     0.000     0.000     1.500          
C2H6               2   252.300     4.302     0.000     0.000     1.500          
HCCO               2   150.000     2.500     0.000     0.000     1.000 ! *      
HCCOH              2   436.000     3.970     0.000     0.000     2.000          
CH2CO              2   436.000     3.970     0.000     0.000     2.000          
CH2CHO             2   436.000     3.970     0.000     0.000     2.000          
C2H2OH             2   224.700     4.162     0.000     0.000     1.000 ! *      
                                                                                
C3H2               2   209.000     4.100     0.000     0.000     1.000 ! *      
C3H3               2   252.000     4.760     0.000     0.000     1.000 ! JAM    
aC3H4              1   252.000     4.760     0.000     0.000     1.000          
pC3H4              1   252.000     4.760     0.000     0.000     1.000          
cC3H4              1   252.000     4.760     0.000     0.000     1.000
                                                                                
C4H                1   357.000     5.180     0.000     0.000     1.000          
C4H2               1   357.000     5.180     0.000     0.000     1.000          
H2C4O              2   357.000     5.180     0.000     0.000     1.000 ! JAM    
C4H2OH             2   224.700     4.162     0.000     0.000     1.000 ! *      
i-C4H3             2   357.000     5.180     0.000     0.000     1.000 ! JAM    
n-C4H3             2   357.000     5.180     0.000     0.000     1.000 ! JAM    
C4H4               2   357.000     5.180     0.000     0.000     1.000 ! JAM    
i-C4H5             2   357.000     5.180     0.000     0.000     1.000 ! JAM    
n-C4H5             2   357.000     5.180     0.000     0.000     1.000 ! JAM    
C4H6               2   357.000     5.180     0.000     0.000     1.000         
C4H612             2   357.000     5.180     0.000     0.000     1.000 
                                                                                
C5H2               1   357.000     5.180     0.000     0.000     1.000          
C5H3               1   357.000     5.180     0.000     0.000     1.000          
C5H5               1   357.000     5.180     0.000     0.000     1.000          
C5H6               1   357.000     5.180     0.000     0.000     1.000          
                                                                                
C6H                1   357.000     5.180     0.000     0.000     1.000          
C6H2               1   357.000     5.180     0.000     0.000     1.000          
C6H3               2   357.000     5.180     0.000     0.000     1.000  !       
l-C6H4             2   412.300     5.349     0.000     0.000     1.000  !(JAM)  
n-C6H5             2   412.300     5.349     0.000     0.000     1.000  !(JAM)  
i-C6H5             2   412.300     5.349     0.000     0.000     1.000  !(JAM)  
l-C6H6             2   412.300     5.349     0.000     0.000     1.000  !(SVE)  
n-C6H7             2   412.300     5.349     0.000     0.000     1.000  !(JAM)  
i-C6H7             2   412.300     5.349     0.000     0.000     1.000  !(JAM)  
C6H8               2   412.300     5.349     0.000     0.000     1.000  !(JAM)  
                                                                                
HE                 0    10.200     2.576     0.000     0.000     0.000 ! *      
H                  0   145.000     2.050     0.000     0.000     0.000          
H2                 1    38.000     2.920     0.000     0.790   280.000          
H2O                2   572.400     2.605     1.844     0.000     4.000          
H2O2               2   107.400     3.458     0.000     0.000     3.800          
HO2                2   107.400     3.458     0.000     0.000     1.000 ! *      
N2                 1    97.530     3.621     0.000     1.760     4.000          
O                  0    80.000     2.750     0.000     0.000     0.000          
O2                 1   107.400     3.458     0.000     1.600     3.800          
OH                 1    80.000     2.750     0.000     0.000     0.000          
                                                                                
                                                                                
The Lennard-Jones parameters of polycyclic aromatic hydrocarbons were estimated 
based on the critical temperature and pressure. See H. Wang and M. Frenklach,   
"Transport Properties of Polycyclic Aromatic Hydrocarbons for Flame Modeling."  
Combustion and Flame, 96:163-170 (1994)                                         
                                                                                
c-C6H4             2   464.8       5.29      0.00     10.32      0.000  !  benze
A1                 2   464.8       5.29      0.00     10.32      0.000  !  benze
A1-                2   464.8       5.29      0.00     10.32      0.000  !  benze
c-C6H7             2   464.8       5.29      0.00     10.32      0.000  !  benze
C5H4O              2   464.8       5.29      0.00     10.32      0.000  !  benze
C5H5O              2   464.8       5.29      0.00     10.32      0.000  !  benze
C5H4OH             2   464.8       5.29      0.00     10.32      0.000  !  benze
C6H5O              2   464.8       5.29      0.00     10.32      0.000  !  benze
C6H5OH             2   464.8       5.29      0.00     10.32      0.000  !  benze

aC3H5              2   266.800     4.982     0.000     0.000     1.000
CH3CCH2            2   266.800     4.982     0.000     0.000     1.000
CH3CHCH            2   266.800     4.982     0.000     0.000     1.000
C3H6               2   266.800     4.982     0.000     0.000     1.000
C3H7               2   266.800     4.982     0.000     0.000     1.000
C4H6               2   357.000     5.180     0.000     0.000     1.000
iC3H7              2   266.800     4.982     0.000     0.000     1.000
nC3H7              2   266.800     4.982     0.000     0.000     1.000
C3H8               2   266.800     4.982     0.000     0.000     1.000
C4H                1   357.000     5.180     0.000     0.000     1.000
C4H2               1   357.000     5.180     0.000     0.000     1.000
C4H2OH             2   224.700     4.162     0.000     0.000     1.000 ! *
iC4H5              2   357.000     5.176     0.000     0.000     1.000
C4H6               2   357.000     5.176     0.000     0.000     1.000
C4H7               2   357.000     5.176     0.000     0.000     1.000
iC4H7              2   357.000     5.176     0.000     0.000     1.000
C4H81              2   357.000     5.176     0.000     0.000     1.000
C4H82              2   357.000     5.176     0.000     0.000     1.000
iC4H8              2   357.000     5.176     0.000     0.000     1.000
tC4H9              2   357.000     5.176     0.000     0.000     1.000
iC4H9              2   357.000     5.176     0.000     0.000     1.000
pC4H9              2   357.000     5.176     0.000     0.000     1.000
sC4H9              2   357.000     5.176     0.000     0.000     1.000
C4H10              2   357.000     5.176     0.000     0.000     1.000
iC4H10             2   357.000     5.176     0.000     0.000     1.000
CH3COCH3           2   357.000     5.176     0.000     0.000     1.000
C2H3CHO            2   357.000     5.176     0.000     0.000     1.000
iC4H7O             2   450.000     5.500     0.000     0.000     1.000 ! JAM
CH3CHO             2   436.000     3.970     0.000     0.000     2.000
CH3CO              2   436.000     3.970     0.000     0.000     2.000
                                                                                
 1-15: Species name                                                             
 16-80: Molecular parameters                                                    
        molecule index: 0 = atom, 1= linear molec.                              
                        2 = nonlinear molec.                                    
        L-J potential well depth, e/kb (K)                                      
        L-J collision diameter, s,                                             
        Dipole moment, f, Debye                                                 
        Polarizability, `, 2                                                   
        Rotational relaxation number, Zrot at 298K                              
        Comments                                                                
                                                                                

#endif
