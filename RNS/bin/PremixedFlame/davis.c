
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
#define CKAWT CKAWT
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
#define GET_REACTION_MAP GET_REACTION_MAP
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
#define CKAWT ckawt
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
#define GET_REACTION_MAP get_reaction_map
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
#define CKAWT ckawt_
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
#define GET_REACTION_MAP get_reaction_map_
#endif

/*function declarations */
void atomicWeight(double * restrict awt);
void molecularWeight(double * restrict wt);
void gibbs(double * restrict species, double * restrict tc);
void helmholtz(double * restrict species, double * restrict tc);
void speciesInternalEnergy(double * restrict species, double * restrict tc);
void speciesEnthalpy(double * restrict species, double * restrict tc);
void speciesEntropy(double * restrict species, double * restrict tc);
void cp_R(double * restrict species, double * restrict tc);
void cv_R(double * restrict species, double * restrict tc);
void equilibriumConstants(double * restrict kc, double * restrict g_RT, double T);
void productionRate(double * restrict wdot, double * restrict sc, double T);
void comp_k_f(double * restrict tc, double invT, double * restrict k_f);
void comp_Kc(double * restrict tc, double invT, double * restrict Kc);
void comp_qfqr(double * restrict q_f, double * restrict q_r, double * restrict sc, double * restrict tc, double invT);
void progressRate(double * restrict qdot, double * restrict speciesConc, double T);
void progressRateFR(double * restrict q_f, double * restrict q_r, double * restrict speciesConc, double T);
void CKINIT();
void CKINDX(int * iwrk, double * restrict rwrk, int * mm, int * kk, int * ii, int * nfit );
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double * restrict rval, int * kerr, int lenline);
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * restrict rval, int * kerr, int lenline, int lenkray);
void CKSYME(int * kname, int * lenkname);
void CKSYMS(int * kname, int * lenkname);
void CKRP(int * ickwrk, double * restrict rckwrk, double * restrict ru, double * restrict ruc, double * restrict pa);
void CKPX(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict P);
void CKPY(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P);
void CKPC(double * restrict rho, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict P);
void CKRHOX(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict rho);
void CKRHOY(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict rho);
void CKRHOC(double * restrict P, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict rho);
void CKWT(int * iwrk, double * restrict rwrk, double * restrict wt);
void CKAWT(int * iwrk, double * restrict rwrk, double * restrict awt);
void CKMMWY(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wtm);
void CKMMWX(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wtm);
void CKMMWC(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict wtm);
void CKYTX(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x);
void CKYTCP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c);
void CKYTCR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c);
void CKXTY(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict y);
void CKXTCP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c);
void CKXTCR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c);
void CKCTX(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict x);
void CKCTY(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict y);
void CKCPOR(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpor);
void CKHORT(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hort);
void CKSOR(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sor);
void CKCVML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvml);
void CKCPML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvml);
void CKUML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict uml);
void CKHML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict uml);
void CKGML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gml);
void CKAML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict aml);
void CKSML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sml);
void CKCVMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvms);
void CKCPMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvms);
void CKUMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums);
void CKHMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums);
void CKGMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gms);
void CKAMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ams);
void CKSMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sms);
void CKCPBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl);
void CKCPBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cpbs);
void CKCVBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl);
void CKCVBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cpbs);
void CKHBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict hbml);
void CKHBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict hbms);
void CKUBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict ubml);
void CKUBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict ubms);
void CKSBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict sbml);
void CKSBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict sbms);
void CKGBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict gbml);
void CKGBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict gbms);
void CKABML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict abml);
void CKABMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict abms);
void CKWC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict wdot);
void CKWYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot);
void CKWXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot);
void CKWYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot);
void CKWXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot);
void CKQC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKKFKR(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict q_f, double * restrict q_r);
void CKQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKNU(int * kdim, int * iwrk, double * restrict rwrk, int * nuki);
void CKNCF(int * mdim, int * iwrk, double * restrict rwrk, int * ncf);
void CKABE(int * iwrk, double * restrict rwrk, double * restrict a, double * restrict b, double * restrict e );
void CKEQC(double * restrict T, double * restrict C , int * iwrk, double * restrict rwrk, double * restrict eqcon );
void CKEQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon);
void CKEQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon);
void CKEQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon);
void CKEQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon);
void DWDOT(double * restrict J, double * restrict sc, double * restrict T, int * consP);
void aJacobian(double * restrict J, double * restrict sc, double T, int consP);
void dcvpRdT(double * restrict species, double * restrict tc);
void GET_T_GIVEN_EY(double * restrict e, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int *ierr);
void GET_REACTION_MAP(int * restrict rmap);
/*vector version */
void vproductionRate(int npt, double * restrict wdot, double * restrict c, double * restrict T);
void VCKHMS(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums);
void VCKPY(int * restrict np, double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P);
void VCKWYR(int * restrict np, double * restrict rho, double * restrict T,
            double * restrict y, int * restrict iwrk, double * restrict rwrk,
            double * restrict wdot);
void VCKYTX(int * restrict np, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x);
void vcomp_k_f(int npt, double * restrict k_f_s, double * restrict tc, double * restrict invT);
void vcomp_gibbs(int npt, double * restrict g_RT, double * restrict tc);
void vcomp_Kc(int npt, double * restrict Kc_s, double * restrict g_RT, double * restrict invT);
void vcomp_wdot(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,
                double * restrict k_f_s, double * restrict Kc_s,
                double * restrict tc, double * restrict invT, double * restrict T);

/* Inverse molecular weights */
static const double imw[14] = {
    1.0 / 2.015940,  /*H2 */
    1.0 / 1.007970,  /*H */
    1.0 / 39.948000,  /*AR */
    1.0 / 28.013400,  /*N2 */
    1.0 / 4.002600,  /*HE */
    1.0 / 15.999400,  /*O */
    1.0 / 17.007370,  /*OH */
    1.0 / 29.018520,  /*HCO */
    1.0 / 33.006770,  /*HO2 */
    1.0 / 18.015340,  /*H2O */
    1.0 / 28.010550,  /*CO */
    1.0 / 31.998800,  /*O2 */
    1.0 / 34.014740,  /*H2O2 */
    1.0 / 44.009950};  /*CO2 */

struct ReactionData {
    double fwd_A,fwd_beta,fwd_Ea;
    double low_A,low_beta,low_Ea;
    double rev_A,rev_beta,rev_Ea;
    double troe_a,troe_Ts, troe_Tss, troe_Tsss;
    double sri_a, sri_b, sri_c, sri_d, sri_e;
    double activation_units, prefactor_units, phase_units;
    int is_PD, troe_len, sri_len;
};

static struct ReactionData R[38], R_DEF[38];

static double fwd_A[38], fwd_beta[38], fwd_Ea[38];
static double low_A[38], low_beta[38], low_Ea[38];
static double rev_A[38], rev_beta[38], rev_Ea[38];
static double troe_a[38],troe_Ts[38], troe_Tss[38], troe_Tsss[38];
static double sri_a[38], sri_b[38], sri_c[38], sri_d[38], sri_e[38];
static double activation_units[38], prefactor_units[38], phase_units[38];
static int is_PD[38], troe_len[38], sri_len[38];
static int rxn_map[38] = {8,9,10,11,3,12,13,14,4,5,6,0,15,1,16,17,18,19,20,21,22,23,24,25,26,27,2,28,29,30,31,32,33,34,35,7,36,37};

void GET_REACTION_MAP(int *rmap)
{
    for (int i=0; i<38; ++i) {
        rmap[i] = rxn_map[i];
    }
}

struct ReactionData* GetReactionData(int id)
{
    if (id<0 || id>=38) {
        printf("GetReactionData: Bad reaction id = %d",id);
        abort();
    };
    return &(R[rxn_map[id]]);
}

struct ReactionData* GetDefaultReactionData(int id)
{
    if (id<0 || id>=38) {
        printf("GetDefaultReactionData: Bad reaction id = %d",id);
        abort();
    };
    return &(R_DEF[rxn_map[id]]);
}

void CopyReactionDataToTranspose(int i, const struct ReactionData * rhs)
{
    fwd_A[i]    = rhs->fwd_A;
    fwd_beta[i] = rhs->fwd_beta;
    fwd_Ea[i]   = rhs->fwd_Ea;
    low_A[i]    = rhs->low_A;
    low_beta[i] = rhs->low_beta;
    low_Ea[i]   = rhs->low_Ea;
    rev_A[i]    = rhs->rev_A;
    rev_beta[i] = rhs->rev_beta;
    rev_Ea[i]   = rhs->rev_Ea;
    troe_a[i]    = rhs->troe_a;
    troe_Ts[i]   = rhs->troe_Ts;
    troe_Tss[i]  = rhs->troe_Tss;
    troe_Tsss[i] = rhs->troe_Tsss;
    sri_a[i] = rhs->sri_a;
    sri_b[i] = rhs->sri_b;
    sri_c[i] = rhs->sri_c;
    sri_d[i] = rhs->sri_d;
    sri_e[i] = rhs->sri_e;
    activation_units[i] = rhs->activation_units;
    prefactor_units[i]  = rhs->prefactor_units;
    phase_units[i]      = rhs->phase_units;
    is_PD[i]    = rhs->is_PD;
    troe_len[i] = rhs->troe_len;
    sri_len[i]  = rhs->sri_len;
}

void SetReactionData(int id, const struct ReactionData * rhs)
{
    if (id<0 || id>=38) {
        printf("SetReactionData: Bad reaction id = %d",id);
        abort();
    }
    R[rxn_map[id]] = *rhs;
    CopyReactionDataToTranspose(rxn_map[id],rhs);
}

/* Initializes static database */
void CKINIT()
{
    R[0].fwd_A     = 5116000000000;
    R[0].fwd_beta  = 0.44;
    R[0].fwd_Ea    = 0;
    R[0].low_A     = 6.328e+19;
    R[0].low_beta  = -1.3999999999999999;
    R[0].low_Ea    = 0;
    R[0].troe_a    = 0.5;
    R[0].troe_Tsss = 1.0000000000000001e-30;
    R[0].troe_Ts   = 1e+30;
    R[0].troe_len  = 3;
    R[0].prefactor_units  = 1.0000000000000002e-06;
    R[0].activation_units = 0.50321666580471969;
    R[0].phase_units      = 1e-12;
    R[0].is_PD = 1;

    R[1].fwd_A     = 111000000000000;
    R[1].fwd_beta  = -0.37;
    R[1].fwd_Ea    = 0;
    R[1].low_A     = 2.01e+17;
    R[1].low_beta  = -0.58399999999999996;
    R[1].low_Ea    = -2293;
    R[1].troe_a    = 0.73460000000000003;
    R[1].troe_Tsss = 94;
    R[1].troe_Ts   = 1756;
    R[1].troe_Tss  = 5182;
    R[1].troe_len  = 4;
    R[1].prefactor_units  = 1.0000000000000002e-06;
    R[1].activation_units = 0.50321666580471969;
    R[1].phase_units      = 1e-12;
    R[1].is_PD = 1;

    R[2].fwd_A     = 13620000000;
    R[2].fwd_beta  = 0;
    R[2].fwd_Ea    = 2384;
    R[2].low_A     = 1.1729999999999999e+24;
    R[2].low_beta  = -2.79;
    R[2].low_Ea    = 4191;
    R[2].prefactor_units  = 1.0000000000000002e-06;
    R[2].activation_units = 0.50321666580471969;
    R[2].phase_units      = 1e-12;
    R[2].is_PD = 1;

    R[3].fwd_A     = 1.78e+18;
    R[3].fwd_beta  = -1;
    R[3].fwd_Ea    = 0;
    R[3].prefactor_units  = 1.0000000000000002e-12;
    R[3].activation_units = 0.50321666580471969;
    R[3].phase_units      = 1e-12;
    R[3].is_PD = 0;

    R[4].fwd_A     = 4.4e+22;
    R[4].fwd_beta  = -2;
    R[4].fwd_Ea    = 0;
    R[4].prefactor_units  = 1.0000000000000002e-12;
    R[4].activation_units = 0.50321666580471969;
    R[4].phase_units      = 1e-12;
    R[4].is_PD = 0;

    R[5].fwd_A     = 9.428e+18;
    R[5].fwd_beta  = -1;
    R[5].fwd_Ea    = 0;
    R[5].prefactor_units  = 1.0000000000000002e-12;
    R[5].activation_units = 0.50321666580471969;
    R[5].phase_units      = 1e-12;
    R[5].is_PD = 0;

    R[6].fwd_A     = 1.2e+17;
    R[6].fwd_beta  = -1;
    R[6].fwd_Ea    = 0;
    R[6].prefactor_units  = 1.0000000000000002e-12;
    R[6].activation_units = 0.50321666580471969;
    R[6].phase_units      = 1e-12;
    R[6].is_PD = 0;

    R[7].fwd_A     = 1.87e+17;
    R[7].fwd_beta  = -1;
    R[7].fwd_Ea    = 17000;
    R[7].prefactor_units  = 1.0000000000000002e-06;
    R[7].activation_units = 0.50321666580471969;
    R[7].phase_units      = 1e-6;
    R[7].is_PD = 0;

    R[8].fwd_A     = 26440000000000000;
    R[8].fwd_beta  = -0.67069999999999996;
    R[8].fwd_Ea    = 17041;
    R[8].prefactor_units  = 1.0000000000000002e-06;
    R[8].activation_units = 0.50321666580471969;
    R[8].phase_units      = 1e-12;
    R[8].is_PD = 0;

    R[9].fwd_A     = 45890;
    R[9].fwd_beta  = 2.7000000000000002;
    R[9].fwd_Ea    = 6260;
    R[9].prefactor_units  = 1.0000000000000002e-06;
    R[9].activation_units = 0.50321666580471969;
    R[9].phase_units      = 1e-12;
    R[9].is_PD = 0;

    R[10].fwd_A     = 173400000;
    R[10].fwd_beta  = 1.51;
    R[10].fwd_Ea    = 3430;
    R[10].prefactor_units  = 1.0000000000000002e-06;
    R[10].activation_units = 0.50321666580471969;
    R[10].phase_units      = 1e-12;
    R[10].is_PD = 0;

    R[11].fwd_A     = 39730;
    R[11].fwd_beta  = 2.3999999999999999;
    R[11].fwd_Ea    = -2110;
    R[11].prefactor_units  = 1.0000000000000002e-06;
    R[11].activation_units = 0.50321666580471969;
    R[11].phase_units      = 1e-12;
    R[11].is_PD = 0;

    R[12].fwd_A     = 90000000000000000;
    R[12].fwd_beta  = -0.59999999999999998;
    R[12].fwd_Ea    = 0;
    R[12].prefactor_units  = 1.0000000000000002e-12;
    R[12].activation_units = 0.50321666580471969;
    R[12].phase_units      = 1e-18;
    R[12].is_PD = 0;

    R[13].fwd_A     = 5.624e+19;
    R[13].fwd_beta  = -1.25;
    R[13].fwd_Ea    = 0;
    R[13].prefactor_units  = 1.0000000000000002e-12;
    R[13].activation_units = 0.50321666580471969;
    R[13].phase_units      = 1e-18;
    R[13].is_PD = 0;

    R[14].fwd_A     = 5.5e+20;
    R[14].fwd_beta  = -2;
    R[14].fwd_Ea    = 0;
    R[14].prefactor_units  = 1.0000000000000002e-12;
    R[14].activation_units = 0.50321666580471969;
    R[14].phase_units      = 1e-18;
    R[14].is_PD = 0;

    R[15].fwd_A     = 591600;
    R[15].fwd_beta  = 2.4329999999999998;
    R[15].fwd_Ea    = 53502;
    R[15].prefactor_units  = 1.0000000000000002e-06;
    R[15].activation_units = 0.50321666580471969;
    R[15].phase_units      = 1e-12;
    R[15].is_PD = 0;

    R[16].fwd_A     = 3970000000000;
    R[16].fwd_beta  = 0;
    R[16].fwd_Ea    = 671;
    R[16].prefactor_units  = 1.0000000000000002e-06;
    R[16].activation_units = 0.50321666580471969;
    R[16].phase_units      = 1e-12;
    R[16].is_PD = 0;

    R[17].fwd_A     = 74850000000000;
    R[17].fwd_beta  = 0;
    R[17].fwd_Ea    = 295;
    R[17].prefactor_units  = 1.0000000000000002e-06;
    R[17].activation_units = 0.50321666580471969;
    R[17].phase_units      = 1e-12;
    R[17].is_PD = 0;

    R[18].fwd_A     = 40000000000000;
    R[18].fwd_beta  = 0;
    R[18].fwd_Ea    = 0;
    R[18].prefactor_units  = 1.0000000000000002e-06;
    R[18].activation_units = 0.50321666580471969;
    R[18].phase_units      = 1e-12;
    R[18].is_PD = 0;

    R[19].fwd_A     = 23750000000000;
    R[19].fwd_beta  = 0;
    R[19].fwd_Ea    = -500;
    R[19].prefactor_units  = 1.0000000000000002e-06;
    R[19].activation_units = 0.50321666580471969;
    R[19].phase_units      = 1e-12;
    R[19].is_PD = 0;

    R[20].fwd_A     = 10000000000000000;
    R[20].fwd_beta  = 0;
    R[20].fwd_Ea    = 17330;
    R[20].prefactor_units  = 1.0000000000000002e-06;
    R[20].activation_units = 0.50321666580471969;
    R[20].phase_units      = 1e-12;
    R[20].is_PD = 0;

    R[21].fwd_A     = 130000000000;
    R[21].fwd_beta  = 0;
    R[21].fwd_Ea    = -1630;
    R[21].prefactor_units  = 1.0000000000000002e-06;
    R[21].activation_units = 0.50321666580471969;
    R[21].phase_units      = 1e-12;
    R[21].is_PD = 0;

    R[22].fwd_A     = 365800000000000;
    R[22].fwd_beta  = 0;
    R[22].fwd_Ea    = 12000;
    R[22].prefactor_units  = 1.0000000000000002e-06;
    R[22].activation_units = 0.50321666580471969;
    R[22].phase_units      = 1e-12;
    R[22].is_PD = 0;

    R[23].fwd_A     = 6050000;
    R[23].fwd_beta  = 2;
    R[23].fwd_Ea    = 5200;
    R[23].prefactor_units  = 1.0000000000000002e-06;
    R[23].activation_units = 0.50321666580471969;
    R[23].phase_units      = 1e-12;
    R[23].is_PD = 0;

    R[24].fwd_A     = 24100000000000;
    R[24].fwd_beta  = 0;
    R[24].fwd_Ea    = 3970;
    R[24].prefactor_units  = 1.0000000000000002e-06;
    R[24].activation_units = 0.50321666580471969;
    R[24].phase_units      = 1e-12;
    R[24].is_PD = 0;

    R[25].fwd_A     = 9630000;
    R[25].fwd_beta  = 2;
    R[25].fwd_Ea    = 3970;
    R[25].prefactor_units  = 1.0000000000000002e-06;
    R[25].activation_units = 0.50321666580471969;
    R[25].phase_units      = 1e-12;
    R[25].is_PD = 0;

    R[26].fwd_A     = 2000000000000;
    R[26].fwd_beta  = 0;
    R[26].fwd_Ea    = 427;
    R[26].prefactor_units  = 1.0000000000000002e-06;
    R[26].activation_units = 0.50321666580471969;
    R[26].phase_units      = 1e-12;
    R[26].is_PD = 0;

    R[27].fwd_A     = 2.6700000000000001e+41;
    R[27].fwd_beta  = -7;
    R[27].fwd_Ea    = 37600;
    R[27].prefactor_units  = 1.0000000000000002e-06;
    R[27].activation_units = 0.50321666580471969;
    R[27].phase_units      = 1e-12;
    R[27].is_PD = 0;

    R[28].fwd_A     = 800000000000;
    R[28].fwd_beta  = 0.14000000000000001;
    R[28].fwd_Ea    = 7352;
    R[28].prefactor_units  = 1.0000000000000002e-06;
    R[28].activation_units = 0.50321666580471969;
    R[28].phase_units      = 1e-12;
    R[28].is_PD = 0;

    R[29].fwd_A     = 87840000000;
    R[29].fwd_beta  = 0.029999999999999999;
    R[29].fwd_Ea    = -16;
    R[29].prefactor_units  = 1.0000000000000002e-06;
    R[29].activation_units = 0.50321666580471969;
    R[29].phase_units      = 1e-12;
    R[29].is_PD = 0;

    R[30].fwd_A     = 1119000000000;
    R[30].fwd_beta  = 0;
    R[30].fwd_Ea    = 47700;
    R[30].prefactor_units  = 1.0000000000000002e-06;
    R[30].activation_units = 0.50321666580471969;
    R[30].phase_units      = 1e-12;
    R[30].is_PD = 0;

    R[31].fwd_A     = 30100000000000;
    R[31].fwd_beta  = 0;
    R[31].fwd_Ea    = 23000;
    R[31].prefactor_units  = 1.0000000000000002e-06;
    R[31].activation_units = 0.50321666580471969;
    R[31].phase_units      = 1e-12;
    R[31].is_PD = 0;

    R[32].fwd_A     = 120000000000000;
    R[32].fwd_beta  = 0;
    R[32].fwd_Ea    = 0;
    R[32].prefactor_units  = 1.0000000000000002e-06;
    R[32].activation_units = 0.50321666580471969;
    R[32].phase_units      = 1e-12;
    R[32].is_PD = 0;

    R[33].fwd_A     = 30000000000000;
    R[33].fwd_beta  = 0;
    R[33].fwd_Ea    = 0;
    R[33].prefactor_units  = 1.0000000000000002e-06;
    R[33].activation_units = 0.50321666580471969;
    R[33].phase_units      = 1e-12;
    R[33].is_PD = 0;

    R[34].fwd_A     = 30000000000000;
    R[34].fwd_beta  = 0;
    R[34].fwd_Ea    = 0;
    R[34].prefactor_units  = 1.0000000000000002e-06;
    R[34].activation_units = 0.50321666580471969;
    R[34].phase_units      = 1e-12;
    R[34].is_PD = 0;

    R[35].fwd_A     = 30200000000000;
    R[35].fwd_beta  = 0;
    R[35].fwd_Ea    = 0;
    R[35].prefactor_units  = 1.0000000000000002e-06;
    R[35].activation_units = 0.50321666580471969;
    R[35].phase_units      = 1e-12;
    R[35].is_PD = 0;

    R[36].fwd_A     = 2.244e+18;
    R[36].fwd_beta  = -1;
    R[36].fwd_Ea    = 17000;
    R[36].prefactor_units  = 1.0000000000000002e-06;
    R[36].activation_units = 0.50321666580471969;
    R[36].phase_units      = 1e-12;
    R[36].is_PD = 0;

    R[37].fwd_A     = 12040000000;
    R[37].fwd_beta  = 0.80700000000000005;
    R[37].fwd_Ea    = -727;
    R[37].prefactor_units  = 1.0000000000000002e-06;
    R[37].activation_units = 0.50321666580471969;
    R[37].phase_units      = 1e-12;
    R[37].is_PD = 0;

    for (int i=0; i<38; i++)
    {
        R_DEF[i] = R[i];
        CopyReactionDataToTranspose(i,&(R[i]));
    }

}


/*A few mechanism parameters */
void CKINDX(int * iwrk, double * restrict rwrk, int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 6;
    *kk = 14;
    *ii = 38;
    *nfit = -1; /*Why do you need this anyway ?  */
}


/* ckxnum... for parsing strings  */
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double * restrict rval, int * kerr, int lenline )
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
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * restrict rval, int * kerr, int lenline, int lenkray)
{
    /*Not done yet ... */
}


/* Returns the char strings of element names */
void CKSYME(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*6; i++) {
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

    /* HE  */
    kname[ 5*lenkname + 0 ] = 'H';
    kname[ 5*lenkname + 1 ] = 'E';
    kname[ 5*lenkname + 2 ] = ' ';

}


/* Returns the char strings of species names */
void CKSYMS(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*14; i++) {
        kname[i] = ' ';
    }

    /* H2  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = '2';
    kname[ 0*lenkname + 2 ] = ' ';

    /* H  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = ' ';

    /* AR  */
    kname[ 2*lenkname + 0 ] = 'A';
    kname[ 2*lenkname + 1 ] = 'R';
    kname[ 2*lenkname + 2 ] = ' ';

    /* N2  */
    kname[ 3*lenkname + 0 ] = 'N';
    kname[ 3*lenkname + 1 ] = '2';
    kname[ 3*lenkname + 2 ] = ' ';

    /* HE  */
    kname[ 4*lenkname + 0 ] = 'H';
    kname[ 4*lenkname + 1 ] = 'E';
    kname[ 4*lenkname + 2 ] = ' ';

    /* O  */
    kname[ 5*lenkname + 0 ] = 'O';
    kname[ 5*lenkname + 1 ] = ' ';

    /* OH  */
    kname[ 6*lenkname + 0 ] = 'O';
    kname[ 6*lenkname + 1 ] = 'H';
    kname[ 6*lenkname + 2 ] = ' ';

    /* HCO  */
    kname[ 7*lenkname + 0 ] = 'H';
    kname[ 7*lenkname + 1 ] = 'C';
    kname[ 7*lenkname + 2 ] = 'O';
    kname[ 7*lenkname + 3 ] = ' ';

    /* HO2  */
    kname[ 8*lenkname + 0 ] = 'H';
    kname[ 8*lenkname + 1 ] = 'O';
    kname[ 8*lenkname + 2 ] = '2';
    kname[ 8*lenkname + 3 ] = ' ';

    /* H2O  */
    kname[ 9*lenkname + 0 ] = 'H';
    kname[ 9*lenkname + 1 ] = '2';
    kname[ 9*lenkname + 2 ] = 'O';
    kname[ 9*lenkname + 3 ] = ' ';

    /* CO  */
    kname[ 10*lenkname + 0 ] = 'C';
    kname[ 10*lenkname + 1 ] = 'O';
    kname[ 10*lenkname + 2 ] = ' ';

    /* O2  */
    kname[ 11*lenkname + 0 ] = 'O';
    kname[ 11*lenkname + 1 ] = '2';
    kname[ 11*lenkname + 2 ] = ' ';

    /* H2O2  */
    kname[ 12*lenkname + 0 ] = 'H';
    kname[ 12*lenkname + 1 ] = '2';
    kname[ 12*lenkname + 2 ] = 'O';
    kname[ 12*lenkname + 3 ] = '2';
    kname[ 12*lenkname + 4 ] = ' ';

    /* CO2  */
    kname[ 13*lenkname + 0 ] = 'C';
    kname[ 13*lenkname + 1 ] = 'O';
    kname[ 13*lenkname + 2 ] = '2';
    kname[ 13*lenkname + 3 ] = ' ';

}


/* Returns R, Rc, Patm */
void CKRP(int * ickwrk, double * restrict rckwrk, double * restrict ru, double * restrict ruc, double * restrict pa)
{
     *ru  = 8.31451e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}


/*Compute P = rhoRT/W(x) */
void CKPX(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict P)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*39.948000; /*AR */
    XW += x[3]*28.013400; /*N2 */
    XW += x[4]*4.002600; /*HE */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*17.007370; /*OH */
    XW += x[7]*29.018520; /*HCO */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*18.015340; /*H2O */
    XW += x[10]*28.010550; /*CO */
    XW += x[11]*31.998800; /*O2 */
    XW += x[12]*34.014740; /*H2O2 */
    XW += x[13]*44.009950; /*CO2 */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void CKPY(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*AR */
    YOW += y[3]*imw[3]; /*N2 */
    YOW += y[4]*imw[4]; /*HE */
    YOW += y[5]*imw[5]; /*O */
    YOW += y[6]*imw[6]; /*OH */
    YOW += y[7]*imw[7]; /*HCO */
    YOW += y[8]*imw[8]; /*HO2 */
    YOW += y[9]*imw[9]; /*H2O */
    YOW += y[10]*imw[10]; /*CO */
    YOW += y[11]*imw[11]; /*O2 */
    YOW += y[12]*imw[12]; /*H2O2 */
    YOW += y[13]*imw[13]; /*CO2 */
    *P = *rho * 8.31451e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void VCKPY(int * restrict np, double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<14; n++) {
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
void CKPC(double * restrict rho, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict P)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*1.007970; /*H */
    W += c[2]*39.948000; /*AR */
    W += c[3]*28.013400; /*N2 */
    W += c[4]*4.002600; /*HE */
    W += c[5]*15.999400; /*O */
    W += c[6]*17.007370; /*OH */
    W += c[7]*29.018520; /*HCO */
    W += c[8]*33.006770; /*HO2 */
    W += c[9]*18.015340; /*H2O */
    W += c[10]*28.010550; /*CO */
    W += c[11]*31.998800; /*O2 */
    W += c[12]*34.014740; /*H2O2 */
    W += c[13]*44.009950; /*CO2 */

    for (id = 0; id < 14; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31451e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*39.948000; /*AR */
    XW += x[3]*28.013400; /*N2 */
    XW += x[4]*4.002600; /*HE */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*17.007370; /*OH */
    XW += x[7]*29.018520; /*HCO */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*18.015340; /*H2O */
    XW += x[10]*28.010550; /*CO */
    XW += x[11]*31.998800; /*O2 */
    XW += x[12]*34.014740; /*H2O2 */
    XW += x[13]*44.009950; /*CO2 */
    *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    double YOW = 0;
    double tmp[14];

    for (int i = 0; i < 14; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 14; i++)
    {
        YOW += tmp[i];
    }

    *rho = *P / (8.31451e+07 * (*T) * YOW);/*rho = P*W/(R*T) */
    return;
}


/*Compute rho = P*W(c)/(R*T) */
void CKRHOC(double * restrict P, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*1.007970; /*H */
    W += c[2]*39.948000; /*AR */
    W += c[3]*28.013400; /*N2 */
    W += c[4]*4.002600; /*HE */
    W += c[5]*15.999400; /*O */
    W += c[6]*17.007370; /*OH */
    W += c[7]*29.018520; /*HCO */
    W += c[8]*33.006770; /*HO2 */
    W += c[9]*18.015340; /*H2O */
    W += c[10]*28.010550; /*CO */
    W += c[11]*31.998800; /*O2 */
    W += c[12]*34.014740; /*H2O2 */
    W += c[13]*44.009950; /*CO2 */

    for (id = 0; id < 14; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.31451e+07); /*rho = PW/(R*T) */

    return;
}


/*get molecular weight for all species */
void CKWT(int * iwrk, double * restrict rwrk, double * restrict wt)
{
    molecularWeight(wt);
}


/*get atomic weight for all elements */
void CKAWT(int * iwrk, double * restrict rwrk, double * restrict awt)
{
    atomicWeight(awt);
}


/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWY(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wtm)
{
    double YOW = 0;
    double tmp[14];

    for (int i = 0; i < 14; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 14; i++)
    {
        YOW += tmp[i];
    }

    *wtm = 1.0 / YOW;
    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWX(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wtm)
{
    double XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*39.948000; /*AR */
    XW += x[3]*28.013400; /*N2 */
    XW += x[4]*4.002600; /*HE */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*17.007370; /*OH */
    XW += x[7]*29.018520; /*HCO */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*18.015340; /*H2O */
    XW += x[10]*28.010550; /*CO */
    XW += x[11]*31.998800; /*O2 */
    XW += x[12]*34.014740; /*H2O2 */
    XW += x[13]*44.009950; /*CO2 */
    *wtm = XW;

    return;
}


/*given c[species]: molar concentration */
/*returns mean molecular weight (gm/mole) */
void CKMMWC(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict wtm)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*1.007970; /*H */
    W += c[2]*39.948000; /*AR */
    W += c[3]*28.013400; /*N2 */
    W += c[4]*4.002600; /*HE */
    W += c[5]*15.999400; /*O */
    W += c[6]*17.007370; /*OH */
    W += c[7]*29.018520; /*HCO */
    W += c[8]*33.006770; /*HO2 */
    W += c[9]*18.015340; /*H2O */
    W += c[10]*28.010550; /*CO */
    W += c[11]*31.998800; /*O2 */
    W += c[12]*34.014740; /*H2O2 */
    W += c[13]*44.009950; /*CO2 */

    for (id = 0; id < 14; ++id) {
        sumC += c[id];
    }
    /* CK provides no guard against divison by zero */
    *wtm = W/sumC;

    return;
}


/*convert y[species] (mass fracs) to x[species] (mole fracs) */
void CKYTX(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x)
{
    double YOW = 0;
    double tmp[14];

    for (int i = 0; i < 14; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 14; i++)
    {
        YOW += tmp[i];
    }

    double YOWINV = 1.0/YOW;

    for (int i = 0; i < 14; i++)
    {
        x[i] = y[i]*imw[i]*YOWINV;
    }
    return;
}


/*convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs) */
void VCKYTX(int * restrict np, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<14; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<14; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c)
{
    double YOW = 0;
    double PWORT;

    /*Compute inverse of mean molecular wt first */
    for (int i = 0; i < 14; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 14; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 14; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c)
{
    for (int i = 0; i < 14; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
void CKXTY(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*39.948000; /*AR */
    XW += x[3]*28.013400; /*N2 */
    XW += x[4]*4.002600; /*HE */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*17.007370; /*OH */
    XW += x[7]*29.018520; /*HCO */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*18.015340; /*H2O */
    XW += x[10]*28.010550; /*CO */
    XW += x[11]*31.998800; /*O2 */
    XW += x[12]*34.014740; /*H2O2 */
    XW += x[13]*44.009950; /*CO2 */
    /*Now compute conversion */
    double XWinv = 1.0/XW;
    y[0] = x[0]*2.015940*XWinv; 
    y[1] = x[1]*1.007970*XWinv; 
    y[2] = x[2]*39.948000*XWinv; 
    y[3] = x[3]*28.013400*XWinv; 
    y[4] = x[4]*4.002600*XWinv; 
    y[5] = x[5]*15.999400*XWinv; 
    y[6] = x[6]*17.007370*XWinv; 
    y[7] = x[7]*29.018520*XWinv; 
    y[8] = x[8]*33.006770*XWinv; 
    y[9] = x[9]*18.015340*XWinv; 
    y[10] = x[10]*28.010550*XWinv; 
    y[11] = x[11]*31.998800*XWinv; 
    y[12] = x[12]*34.014740*XWinv; 
    y[13] = x[13]*44.009950*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 14; ++id) {
        c[id] = x[id]*PORT;
    }

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c)
{
    int id; /*loop counter */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*39.948000; /*AR */
    XW += x[3]*28.013400; /*N2 */
    XW += x[4]*4.002600; /*HE */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*17.007370; /*OH */
    XW += x[7]*29.018520; /*HCO */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*18.015340; /*H2O */
    XW += x[10]*28.010550; /*CO */
    XW += x[11]*31.998800; /*O2 */
    XW += x[12]*34.014740; /*H2O2 */
    XW += x[13]*44.009950; /*CO2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 14; ++id) {
        c[id] = x[id]*ROW;
    }

    return;
}


/*convert c[species] (molar conc) to x[species] (mole fracs) */
void CKCTX(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict x)
{
    int id; /*loop counter */
    double sumC = 0; 

    /*compute sum of c  */
    for (id = 0; id < 14; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    double sumCinv = 1.0/sumC;
    for (id = 0; id < 14; ++id) {
        x[id] = c[id]*sumCinv;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
void CKCTY(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict y)
{
    double CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*2.015940; /*H2 */
    CW += c[1]*1.007970; /*H */
    CW += c[2]*39.948000; /*AR */
    CW += c[3]*28.013400; /*N2 */
    CW += c[4]*4.002600; /*HE */
    CW += c[5]*15.999400; /*O */
    CW += c[6]*17.007370; /*OH */
    CW += c[7]*29.018520; /*HCO */
    CW += c[8]*33.006770; /*HO2 */
    CW += c[9]*18.015340; /*H2O */
    CW += c[10]*28.010550; /*CO */
    CW += c[11]*31.998800; /*O2 */
    CW += c[12]*34.014740; /*H2O2 */
    CW += c[13]*44.009950; /*CO2 */
    /*Now compute conversion */
    double CWinv = 1.0/CW;
    y[0] = c[0]*2.015940*CWinv; 
    y[1] = c[1]*1.007970*CWinv; 
    y[2] = c[2]*39.948000*CWinv; 
    y[3] = c[3]*28.013400*CWinv; 
    y[4] = c[4]*4.002600*CWinv; 
    y[5] = c[5]*15.999400*CWinv; 
    y[6] = c[6]*17.007370*CWinv; 
    y[7] = c[7]*29.018520*CWinv; 
    y[8] = c[8]*33.006770*CWinv; 
    y[9] = c[9]*18.015340*CWinv; 
    y[10] = c[10]*28.010550*CWinv; 
    y[11] = c[11]*31.998800*CWinv; 
    y[12] = c[12]*34.014740*CWinv; 
    y[13] = c[13]*44.009950*CWinv; 

    return;
}


/*get Cp/R as a function of T  */
/*for all species (Eq 19) */
void CKCPOR(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpor, tc);
}


/*get H/RT as a function of T  */
/*for all species (Eq 20) */
void CKHORT(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hort)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEnthalpy(hort, tc);
}


/*get S/R as a function of T  */
/*for all species (Eq 21) */
void CKSOR(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sor, tc);
}


/*get specific heat at constant volume as a function  */
/*of T for all species (molar units) */
void CKCVML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
        cvml[id] *= 8.31451e+07;
    }
}


/*get specific heat at constant pressure as a  */
/*function of T for all species (molar units) */
void CKCPML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
        cpml[id] *= 8.31451e+07;
    }
}


/*get internal energy as a function  */
/*of T for all species (molar units) */
void CKUML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict uml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
        uml[id] *= RT;
    }
}


/*get enthalpy as a function  */
/*of T for all species (molar units) */
void CKHML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
        hml[id] *= RT;
    }
}


/*get standard-state Gibbs energy as a function  */
/*of T for all species (molar units) */
void CKGML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
        gml[id] *= RT;
    }
}


/*get standard-state Helmholtz free energy as a  */
/*function of T for all species (molar units) */
void CKAML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict aml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(aml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
        aml[id] *= RT;
    }
}


/*Returns the standard-state entropies in molar units */
void CKSML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
        sml[id] *= 8.31451e+07;
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
void CKCVMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 4.124383662212169e+07; /*H2 */
    cvms[1] *= 8.248767324424338e+07; /*H */
    cvms[2] *= 2.081333233203164e+06; /*AR */
    cvms[3] *= 2.968047434442088e+06; /*N2 */
    cvms[4] *= 2.077277269774646e+07; /*HE */
    cvms[5] *= 5.196763628636074e+06; /*O */
    cvms[6] *= 4.888768810227566e+06; /*OH */
    cvms[7] *= 2.865242610581105e+06; /*HCO */
    cvms[8] *= 2.519031701678171e+06; /*HO2 */
    cvms[9] *= 4.615239012974499e+06; /*H2O */
    cvms[10] *= 2.968349425484326e+06; /*CO */
    cvms[11] *= 2.598381814318037e+06; /*O2 */
    cvms[12] *= 2.444384405113783e+06; /*H2O2 */
    cvms[13] *= 1.889234139098090e+06; /*CO2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
void CKCPMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 4.124383662212169e+07; /*H2 */
    cpms[1] *= 8.248767324424338e+07; /*H */
    cpms[2] *= 2.081333233203164e+06; /*AR */
    cpms[3] *= 2.968047434442088e+06; /*N2 */
    cpms[4] *= 2.077277269774646e+07; /*HE */
    cpms[5] *= 5.196763628636074e+06; /*O */
    cpms[6] *= 4.888768810227566e+06; /*OH */
    cpms[7] *= 2.865242610581105e+06; /*HCO */
    cpms[8] *= 2.519031701678171e+06; /*HO2 */
    cpms[9] *= 4.615239012974499e+06; /*H2O */
    cpms[10] *= 2.968349425484326e+06; /*CO */
    cpms[11] *= 2.598381814318037e+06; /*O2 */
    cpms[12] *= 2.444384405113783e+06; /*H2O2 */
    cpms[13] *= 1.889234139098090e+06; /*CO2 */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    for (int i = 0; i < 14; i++)
    {
        ums[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void CKHMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hms, tc);
    for (int i = 0; i < 14; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hms)
{
    double tc[5], h[14];

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
        hms[9*(*np)+i] = h[9];
        hms[10*(*np)+i] = h[10];
        hms[11*(*np)+i] = h[11];
        hms[12*(*np)+i] = h[12];
        hms[13*(*np)+i] = h[13];
    }

    for (int n=0; n<14; n++) {
        for (int i=0; i<(*np); i++) {
            hms[n*(*np)+i] *= 8.31451e+07 * T[i] * imw[n];
        }
    }
}


/*Returns gibbs in mass units (Eq 31.) */
void CKGMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gms, tc);
    for (int i = 0; i < 14; i++)
    {
        gms[i] *= RT*imw[i];
    }
}


/*Returns helmholtz in mass units (Eq 32.) */
void CKAMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ams)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(ams, tc);
    for (int i = 0; i < 14; i++)
    {
        ams[i] *= RT*imw[i];
    }
}


/*Returns the entropies in mass units (Eq 28.) */
void CKSMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 4.124383662212169e+07; /*H2 */
    sms[1] *= 8.248767324424338e+07; /*H */
    sms[2] *= 2.081333233203164e+06; /*AR */
    sms[3] *= 2.968047434442088e+06; /*N2 */
    sms[4] *= 2.077277269774646e+07; /*HE */
    sms[5] *= 5.196763628636074e+06; /*O */
    sms[6] *= 4.888768810227566e+06; /*OH */
    sms[7] *= 2.865242610581105e+06; /*HCO */
    sms[8] *= 2.519031701678171e+06; /*HO2 */
    sms[9] *= 4.615239012974499e+06; /*H2O */
    sms[10] *= 2.968349425484326e+06; /*CO */
    sms[11] *= 2.598381814318037e+06; /*O2 */
    sms[12] *= 2.444384405113783e+06; /*H2O2 */
    sms[13] *= 1.889234139098090e+06; /*CO2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[14]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 14; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
void CKCPBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cpbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[14], tresult[14]; /* temporary storage */
    cp_R(cpor, tc);
    for (int i = 0; i < 14; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 14; i++)
    {
        result += tresult[i];
    }

    *cpbs = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[14]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 14; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
void CKCVBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cvbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[14]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H2 */
    result += cvor[1]*y[1]*imw[1]; /*H */
    result += cvor[2]*y[2]*imw[2]; /*AR */
    result += cvor[3]*y[3]*imw[3]; /*N2 */
    result += cvor[4]*y[4]*imw[4]; /*HE */
    result += cvor[5]*y[5]*imw[5]; /*O */
    result += cvor[6]*y[6]*imw[6]; /*OH */
    result += cvor[7]*y[7]*imw[7]; /*HCO */
    result += cvor[8]*y[8]*imw[8]; /*HO2 */
    result += cvor[9]*y[9]*imw[9]; /*H2O */
    result += cvor[10]*y[10]*imw[10]; /*CO */
    result += cvor[11]*y[11]*imw[11]; /*O2 */
    result += cvor[12]*y[12]*imw[12]; /*H2O2 */
    result += cvor[13]*y[13]*imw[13]; /*CO2 */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[14]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 14; ++id) {
        result += x[id]*hml[id];
    }

    *hbml = result * RT;
}


/*Returns mean enthalpy of mixture in mass units */
void CKHBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict hbms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[14], tmp[14]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    int id;
    for (id = 0; id < 14; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 14; ++id) {
        result += tmp[id];
    }

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[14]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 14; ++id) {
        result += x[id]*uml[id];
    }

    *ubml = result * RT;
}


/*get mean internal energy in mass units */
void CKUBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict ubms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double ums[14]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H2 */
    result += y[1]*ums[1]*imw[1]; /*H */
    result += y[2]*ums[2]*imw[2]; /*AR */
    result += y[3]*ums[3]*imw[3]; /*N2 */
    result += y[4]*ums[4]*imw[4]; /*HE */
    result += y[5]*ums[5]*imw[5]; /*O */
    result += y[6]*ums[6]*imw[6]; /*OH */
    result += y[7]*ums[7]*imw[7]; /*HCO */
    result += y[8]*ums[8]*imw[8]; /*HO2 */
    result += y[9]*ums[9]*imw[9]; /*H2O */
    result += y[10]*ums[10]*imw[10]; /*CO */
    result += y[11]*ums[11]*imw[11]; /*O2 */
    result += y[12]*ums[12]*imw[12]; /*H2O2 */
    result += y[13]*ums[13]*imw[13]; /*CO2 */

    *ubms = result * RT;
}


/*get mixture entropy in molar units */
void CKSBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict sbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[14]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 14; ++id) {
        result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
    }

    *sbml = result * 8.31451e+07;
}


/*get mixture entropy in mass units */
void CKSBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict sbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[14]; /* temporary storage */
    double x[14]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*AR */
    YOW += y[3]*imw[3]; /*N2 */
    YOW += y[4]*imw[4]; /*HE */
    YOW += y[5]*imw[5]; /*O */
    YOW += y[6]*imw[6]; /*OH */
    YOW += y[7]*imw[7]; /*HCO */
    YOW += y[8]*imw[8]; /*HO2 */
    YOW += y[9]*imw[9]; /*H2O */
    YOW += y[10]*imw[10]; /*CO */
    YOW += y[11]*imw[11]; /*O2 */
    YOW += y[12]*imw[12]; /*H2O2 */
    YOW += y[13]*imw[13]; /*CO2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(39.948000*YOW); 
    x[3] = y[3]/(28.013400*YOW); 
    x[4] = y[4]/(4.002600*YOW); 
    x[5] = y[5]/(15.999400*YOW); 
    x[6] = y[6]/(17.007370*YOW); 
    x[7] = y[7]/(29.018520*YOW); 
    x[8] = y[8]/(33.006770*YOW); 
    x[9] = y[9]/(18.015340*YOW); 
    x[10] = y[10]/(28.010550*YOW); 
    x[11] = y[11]/(31.998800*YOW); 
    x[12] = y[12]/(34.014740*YOW); 
    x[13] = y[13]/(44.009950*YOW); 
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
    /*Scale by R/W */
    *sbms = result * 8.31451e+07 * YOW;
}


/*Returns mean gibbs free energy in molar units */
void CKGBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict gbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double gort[14]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 14; ++id) {
        result += x[id]*(gort[id]+log((x[id]+1e-100))+logPratio);
    }

    *gbml = result * RT;
}


/*Returns mixture gibbs free energy in mass units */
void CKGBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict gbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double gort[14]; /* temporary storage */
    double x[14]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*AR */
    YOW += y[3]*imw[3]; /*N2 */
    YOW += y[4]*imw[4]; /*HE */
    YOW += y[5]*imw[5]; /*O */
    YOW += y[6]*imw[6]; /*OH */
    YOW += y[7]*imw[7]; /*HCO */
    YOW += y[8]*imw[8]; /*HO2 */
    YOW += y[9]*imw[9]; /*H2O */
    YOW += y[10]*imw[10]; /*CO */
    YOW += y[11]*imw[11]; /*O2 */
    YOW += y[12]*imw[12]; /*H2O2 */
    YOW += y[13]*imw[13]; /*CO2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(39.948000*YOW); 
    x[3] = y[3]/(28.013400*YOW); 
    x[4] = y[4]/(4.002600*YOW); 
    x[5] = y[5]/(15.999400*YOW); 
    x[6] = y[6]/(17.007370*YOW); 
    x[7] = y[7]/(29.018520*YOW); 
    x[8] = y[8]/(33.006770*YOW); 
    x[9] = y[9]/(18.015340*YOW); 
    x[10] = y[10]/(28.010550*YOW); 
    x[11] = y[11]/(31.998800*YOW); 
    x[12] = y[12]/(34.014740*YOW); 
    x[13] = y[13]/(44.009950*YOW); 
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
    /*Scale by RT/W */
    *gbms = result * RT * YOW;
}


/*Returns mean helmholtz free energy in molar units */
void CKABML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict abml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double aort[14]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 14; ++id) {
        result += x[id]*(aort[id]+log((x[id]+1e-100))+logPratio);
    }

    *abml = result * RT;
}


/*Returns mixture helmholtz free energy in mass units */
void CKABMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict abms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double aort[14]; /* temporary storage */
    double x[14]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*AR */
    YOW += y[3]*imw[3]; /*N2 */
    YOW += y[4]*imw[4]; /*HE */
    YOW += y[5]*imw[5]; /*O */
    YOW += y[6]*imw[6]; /*OH */
    YOW += y[7]*imw[7]; /*HCO */
    YOW += y[8]*imw[8]; /*HO2 */
    YOW += y[9]*imw[9]; /*H2O */
    YOW += y[10]*imw[10]; /*CO */
    YOW += y[11]*imw[11]; /*O2 */
    YOW += y[12]*imw[12]; /*H2O2 */
    YOW += y[13]*imw[13]; /*CO2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(39.948000*YOW); 
    x[3] = y[3]/(28.013400*YOW); 
    x[4] = y[4]/(4.002600*YOW); 
    x[5] = y[5]/(15.999400*YOW); 
    x[6] = y[6]/(17.007370*YOW); 
    x[7] = y[7]/(29.018520*YOW); 
    x[8] = y[8]/(33.006770*YOW); 
    x[9] = y[9]/(18.015340*YOW); 
    x[10] = y[10]/(28.010550*YOW); 
    x[11] = y[11]/(31.998800*YOW); 
    x[12] = y[12]/(34.014740*YOW); 
    x[13] = y[13]/(44.009950*YOW); 
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
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 14; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[14]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*AR */
    YOW += y[3]*imw[3]; /*N2 */
    YOW += y[4]*imw[4]; /*HE */
    YOW += y[5]*imw[5]; /*O */
    YOW += y[6]*imw[6]; /*OH */
    YOW += y[7]*imw[7]; /*HCO */
    YOW += y[8]*imw[8]; /*HO2 */
    YOW += y[9]*imw[9]; /*H2O */
    YOW += y[10]*imw[10]; /*CO */
    YOW += y[11]*imw[11]; /*O2 */
    YOW += y[12]*imw[12]; /*H2O2 */
    YOW += y[13]*imw[13]; /*CO2 */
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
    c[9] = PWORT * y[9]*imw[9]; 
    c[10] = PWORT * y[10]*imw[10]; 
    c[11] = PWORT * y[11]*imw[11]; 
    c[12] = PWORT * y[12]*imw[12]; 
    c[13] = PWORT * y[13]*imw[13]; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[14]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 14; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[14]; /*temporary storage */
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
    c[9] = 1e6 * (*rho) * y[9]*imw[9]; 
    c[10] = 1e6 * (*rho) * y[10]*imw[10]; 
    c[11] = 1e6 * (*rho) * y[11]*imw[11]; 
    c[12] = 1e6 * (*rho) * y[12]*imw[12]; 
    c[13] = 1e6 * (*rho) * y[13]*imw[13]; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int * restrict np, double * restrict rho, double * restrict T,
	    double * restrict y, int * restrict iwrk, double * restrict rwrk,
	    double * restrict wdot)
{
    double c[14*(*np)]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<14; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<14*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[14]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*39.948000; /*AR */
    XW += x[3]*28.013400; /*N2 */
    XW += x[4]*4.002600; /*HE */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*17.007370; /*OH */
    XW += x[7]*29.018520; /*HCO */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*18.015340; /*H2O */
    XW += x[10]*28.010550; /*CO */
    XW += x[11]*31.998800; /*O2 */
    XW += x[12]*34.014740; /*H2O2 */
    XW += x[13]*44.009950; /*CO2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 14; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 14; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 38; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict q_f, double * restrict q_r)
{
    int id; /*loop counter */
    double c[14]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 14; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 38; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[14]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*AR */
    YOW += y[3]*imw[3]; /*N2 */
    YOW += y[4]*imw[4]; /*HE */
    YOW += y[5]*imw[5]; /*O */
    YOW += y[6]*imw[6]; /*OH */
    YOW += y[7]*imw[7]; /*HCO */
    YOW += y[8]*imw[8]; /*HO2 */
    YOW += y[9]*imw[9]; /*H2O */
    YOW += y[10]*imw[10]; /*CO */
    YOW += y[11]*imw[11]; /*O2 */
    YOW += y[12]*imw[12]; /*H2O2 */
    YOW += y[13]*imw[13]; /*CO2 */
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
    c[9] = PWORT * y[9]*imw[9]; 
    c[10] = PWORT * y[10]*imw[10]; 
    c[11] = PWORT * y[11]*imw[11]; 
    c[12] = PWORT * y[12]*imw[12]; 
    c[13] = PWORT * y[13]*imw[13]; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 38; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[14]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 14; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 38; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[14]; /*temporary storage */
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
    c[9] = 1e6 * (*rho) * y[9]*imw[9]; 
    c[10] = 1e6 * (*rho) * y[10]*imw[10]; 
    c[11] = 1e6 * (*rho) * y[11]*imw[11]; 
    c[12] = 1e6 * (*rho) * y[12]*imw[12]; 
    c[13] = 1e6 * (*rho) * y[13]*imw[13]; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 38; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[14]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*39.948000; /*AR */
    XW += x[3]*28.013400; /*N2 */
    XW += x[4]*4.002600; /*HE */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*17.007370; /*OH */
    XW += x[7]*29.018520; /*HCO */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*18.015340; /*H2O */
    XW += x[10]*28.010550; /*CO */
    XW += x[11]*31.998800; /*O2 */
    XW += x[12]*34.014740; /*H2O2 */
    XW += x[13]*44.009950; /*CO2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 14; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 38; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the stoichiometric coefficients */
/*of the reaction mechanism. (Eq 50) */
void CKNU(int * kdim, int * iwrk, double * restrict rwrk, int * nuki)
{
    int id; /*loop counter */
    int kd = (*kdim); 
    /*Zero nuki */
    for (id = 0; id < 14 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    nuki[ 1 * kd + 0 ] += -1 ;
    nuki[ 11 * kd + 0 ] += -1 ;
    nuki[ 8 * kd + 0 ] += +1 ;

    /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
    nuki[ 6 * kd + 1 ] += -1 ;
    nuki[ 6 * kd + 1 ] += -1 ;
    nuki[ 12 * kd + 1 ] += +1 ;

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    nuki[ 10 * kd + 2 ] += -1 ;
    nuki[ 5 * kd + 2 ] += -1 ;
    nuki[ 13 * kd + 2 ] += +1 ;

    /*reaction 4: H + H + M <=> H2 + M */
    nuki[ 1 * kd + 3 ] += -1 ;
    nuki[ 1 * kd + 3 ] += -1 ;
    nuki[ 0 * kd + 3 ] += +1 ;

    /*reaction 5: H + OH + M <=> H2O + M */
    nuki[ 1 * kd + 4 ] += -1 ;
    nuki[ 6 * kd + 4 ] += -1 ;
    nuki[ 9 * kd + 4 ] += +1 ;

    /*reaction 6: O + H + M <=> OH + M */
    nuki[ 5 * kd + 5 ] += -1 ;
    nuki[ 1 * kd + 5 ] += -1 ;
    nuki[ 6 * kd + 5 ] += +1 ;

    /*reaction 7: O + O + M <=> O2 + M */
    nuki[ 5 * kd + 6 ] += -1 ;
    nuki[ 5 * kd + 6 ] += -1 ;
    nuki[ 11 * kd + 6 ] += +1 ;

    /*reaction 8: HCO + M <=> CO + H + M */
    nuki[ 7 * kd + 7 ] += -1 ;
    nuki[ 10 * kd + 7 ] += +1 ;
    nuki[ 1 * kd + 7 ] += +1 ;

    /*reaction 9: H + O2 <=> O + OH */
    nuki[ 1 * kd + 8 ] += -1 ;
    nuki[ 11 * kd + 8 ] += -1 ;
    nuki[ 5 * kd + 8 ] += +1 ;
    nuki[ 6 * kd + 8 ] += +1 ;

    /*reaction 10: O + H2 <=> H + OH */
    nuki[ 5 * kd + 9 ] += -1 ;
    nuki[ 0 * kd + 9 ] += -1 ;
    nuki[ 1 * kd + 9 ] += +1 ;
    nuki[ 6 * kd + 9 ] += +1 ;

    /*reaction 11: OH + H2 <=> H + H2O */
    nuki[ 6 * kd + 10 ] += -1 ;
    nuki[ 0 * kd + 10 ] += -1 ;
    nuki[ 1 * kd + 10 ] += +1 ;
    nuki[ 9 * kd + 10 ] += +1 ;

    /*reaction 12: OH + OH <=> O + H2O */
    nuki[ 6 * kd + 11 ] += -1 ;
    nuki[ 6 * kd + 11 ] += -1 ;
    nuki[ 5 * kd + 11 ] += +1 ;
    nuki[ 9 * kd + 11 ] += +1 ;

    /*reaction 13: H + H + H2 <=> H2 + H2 */
    nuki[ 1 * kd + 12 ] += -1 ;
    nuki[ 1 * kd + 12 ] += -1 ;
    nuki[ 0 * kd + 12 ] += -1 ;
    nuki[ 0 * kd + 12 ] += +1 ;
    nuki[ 0 * kd + 12 ] += +1 ;

    /*reaction 14: H + H + H2O <=> H2 + H2O */
    nuki[ 1 * kd + 13 ] += -1 ;
    nuki[ 1 * kd + 13 ] += -1 ;
    nuki[ 9 * kd + 13 ] += -1 ;
    nuki[ 0 * kd + 13 ] += +1 ;
    nuki[ 9 * kd + 13 ] += +1 ;

    /*reaction 15: H + H + CO2 <=> H2 + CO2 */
    nuki[ 1 * kd + 14 ] += -1 ;
    nuki[ 1 * kd + 14 ] += -1 ;
    nuki[ 13 * kd + 14 ] += -1 ;
    nuki[ 0 * kd + 14 ] += +1 ;
    nuki[ 13 * kd + 14 ] += +1 ;

    /*reaction 16: H2 + O2 <=> HO2 + H */
    nuki[ 0 * kd + 15 ] += -1 ;
    nuki[ 11 * kd + 15 ] += -1 ;
    nuki[ 8 * kd + 15 ] += +1 ;
    nuki[ 1 * kd + 15 ] += +1 ;

    /*reaction 17: HO2 + H <=> O + H2O */
    nuki[ 8 * kd + 16 ] += -1 ;
    nuki[ 1 * kd + 16 ] += -1 ;
    nuki[ 5 * kd + 16 ] += +1 ;
    nuki[ 9 * kd + 16 ] += +1 ;

    /*reaction 18: HO2 + H <=> OH + OH */
    nuki[ 8 * kd + 17 ] += -1 ;
    nuki[ 1 * kd + 17 ] += -1 ;
    nuki[ 6 * kd + 17 ] += +1 ;
    nuki[ 6 * kd + 17 ] += +1 ;

    /*reaction 19: HO2 + O <=> OH + O2 */
    nuki[ 8 * kd + 18 ] += -1 ;
    nuki[ 5 * kd + 18 ] += -1 ;
    nuki[ 6 * kd + 18 ] += +1 ;
    nuki[ 11 * kd + 18 ] += +1 ;

    /*reaction 20: HO2 + OH <=> O2 + H2O */
    nuki[ 8 * kd + 19 ] += -1 ;
    nuki[ 6 * kd + 19 ] += -1 ;
    nuki[ 11 * kd + 19 ] += +1 ;
    nuki[ 9 * kd + 19 ] += +1 ;

    /*reaction 21: HO2 + OH <=> O2 + H2O */
    nuki[ 8 * kd + 20 ] += -1 ;
    nuki[ 6 * kd + 20 ] += -1 ;
    nuki[ 11 * kd + 20 ] += +1 ;
    nuki[ 9 * kd + 20 ] += +1 ;

    /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
    nuki[ 8 * kd + 21 ] += -1 ;
    nuki[ 8 * kd + 21 ] += -1 ;
    nuki[ 11 * kd + 21 ] += +1 ;
    nuki[ 12 * kd + 21 ] += +1 ;

    /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
    nuki[ 8 * kd + 22 ] += -1 ;
    nuki[ 8 * kd + 22 ] += -1 ;
    nuki[ 11 * kd + 22 ] += +1 ;
    nuki[ 12 * kd + 22 ] += +1 ;

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    nuki[ 12 * kd + 23 ] += -1 ;
    nuki[ 1 * kd + 23 ] += -1 ;
    nuki[ 8 * kd + 23 ] += +1 ;
    nuki[ 0 * kd + 23 ] += +1 ;

    /*reaction 25: H2O2 + H <=> OH + H2O */
    nuki[ 12 * kd + 24 ] += -1 ;
    nuki[ 1 * kd + 24 ] += -1 ;
    nuki[ 6 * kd + 24 ] += +1 ;
    nuki[ 9 * kd + 24 ] += +1 ;

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    nuki[ 12 * kd + 25 ] += -1 ;
    nuki[ 5 * kd + 25 ] += -1 ;
    nuki[ 6 * kd + 25 ] += +1 ;
    nuki[ 8 * kd + 25 ] += +1 ;

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    nuki[ 12 * kd + 26 ] += -1 ;
    nuki[ 6 * kd + 26 ] += -1 ;
    nuki[ 8 * kd + 26 ] += +1 ;
    nuki[ 9 * kd + 26 ] += +1 ;

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    nuki[ 12 * kd + 27 ] += -1 ;
    nuki[ 6 * kd + 27 ] += -1 ;
    nuki[ 8 * kd + 27 ] += +1 ;
    nuki[ 9 * kd + 27 ] += +1 ;

    /*reaction 29: CO + OH <=> CO2 + H */
    nuki[ 10 * kd + 28 ] += -1 ;
    nuki[ 6 * kd + 28 ] += -1 ;
    nuki[ 13 * kd + 28 ] += +1 ;
    nuki[ 1 * kd + 28 ] += +1 ;

    /*reaction 30: CO + OH <=> CO2 + H */
    nuki[ 10 * kd + 29 ] += -1 ;
    nuki[ 6 * kd + 29 ] += -1 ;
    nuki[ 13 * kd + 29 ] += +1 ;
    nuki[ 1 * kd + 29 ] += +1 ;

    /*reaction 31: CO + O2 <=> CO2 + O */
    nuki[ 10 * kd + 30 ] += -1 ;
    nuki[ 11 * kd + 30 ] += -1 ;
    nuki[ 13 * kd + 30 ] += +1 ;
    nuki[ 5 * kd + 30 ] += +1 ;

    /*reaction 32: CO + HO2 <=> CO2 + OH */
    nuki[ 10 * kd + 31 ] += -1 ;
    nuki[ 8 * kd + 31 ] += -1 ;
    nuki[ 13 * kd + 31 ] += +1 ;
    nuki[ 6 * kd + 31 ] += +1 ;

    /*reaction 33: HCO + H <=> CO + H2 */
    nuki[ 7 * kd + 32 ] += -1 ;
    nuki[ 1 * kd + 32 ] += -1 ;
    nuki[ 10 * kd + 32 ] += +1 ;
    nuki[ 0 * kd + 32 ] += +1 ;

    /*reaction 34: HCO + O <=> CO + OH */
    nuki[ 7 * kd + 33 ] += -1 ;
    nuki[ 5 * kd + 33 ] += -1 ;
    nuki[ 10 * kd + 33 ] += +1 ;
    nuki[ 6 * kd + 33 ] += +1 ;

    /*reaction 35: HCO + O <=> CO2 + H */
    nuki[ 7 * kd + 34 ] += -1 ;
    nuki[ 5 * kd + 34 ] += -1 ;
    nuki[ 13 * kd + 34 ] += +1 ;
    nuki[ 1 * kd + 34 ] += +1 ;

    /*reaction 36: HCO + OH <=> CO + H2O */
    nuki[ 7 * kd + 35 ] += -1 ;
    nuki[ 6 * kd + 35 ] += -1 ;
    nuki[ 10 * kd + 35 ] += +1 ;
    nuki[ 9 * kd + 35 ] += +1 ;

    /*reaction 37: HCO + H2O <=> CO + H + H2O */
    nuki[ 7 * kd + 36 ] += -1 ;
    nuki[ 9 * kd + 36 ] += -1 ;
    nuki[ 10 * kd + 36 ] += +1 ;
    nuki[ 1 * kd + 36 ] += +1 ;
    nuki[ 9 * kd + 36 ] += +1 ;

    /*reaction 38: HCO + O2 <=> CO + HO2 */
    nuki[ 7 * kd + 37 ] += -1 ;
    nuki[ 11 * kd + 37 ] += -1 ;
    nuki[ 10 * kd + 37 ] += +1 ;
    nuki[ 8 * kd + 37 ] += +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * restrict rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < kd * 14; ++ id) {
         ncf[id] = 0; 
    }

    /*H2 */
    ncf[ 0 * kd + 1 ] = 2; /*H */

    /*H */
    ncf[ 1 * kd + 1 ] = 1; /*H */

    /*AR */
    ncf[ 2 * kd + 4 ] = 1; /*AR */

    /*N2 */
    ncf[ 3 * kd + 3 ] = 2; /*N */

    /*HE */
    ncf[ 4 * kd + 5 ] = 1; /*HE */

    /*O */
    ncf[ 5 * kd + 0 ] = 1; /*O */

    /*OH */
    ncf[ 6 * kd + 0 ] = 1; /*O */
    ncf[ 6 * kd + 1 ] = 1; /*H */

    /*HCO */
    ncf[ 7 * kd + 1 ] = 1; /*H */
    ncf[ 7 * kd + 2 ] = 1; /*C */
    ncf[ 7 * kd + 0 ] = 1; /*O */

    /*HO2 */
    ncf[ 8 * kd + 1 ] = 1; /*H */
    ncf[ 8 * kd + 0 ] = 2; /*O */

    /*H2O */
    ncf[ 9 * kd + 1 ] = 2; /*H */
    ncf[ 9 * kd + 0 ] = 1; /*O */

    /*CO */
    ncf[ 10 * kd + 2 ] = 1; /*C */
    ncf[ 10 * kd + 0 ] = 1; /*O */

    /*O2 */
    ncf[ 11 * kd + 0 ] = 2; /*O */

    /*H2O2 */
    ncf[ 12 * kd + 1 ] = 2; /*H */
    ncf[ 12 * kd + 0 ] = 2; /*O */

    /*CO2 */
    ncf[ 13 * kd + 2 ] = 1; /*C */
    ncf[ 13 * kd + 0 ] = 2; /*O */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE(int * iwrk, double * restrict rwrk, double * restrict a, double * restrict b, double * restrict e)
{
    for (int i=0; i<38; ++i) {
        a[i] = fwd_A[i];
        b[i] = fwd_beta[i];
        e[i] = fwd_Ea[i];
    }

    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[14]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + M <=> H2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + OH + M <=> H2O + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + O + M <=> O2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: HCO + M <=> CO + H + M */
    eqcon[7] *= 1e-06; 

    /*reaction 9: H + O2 <=> O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: OH + H2 <=> H + H2O */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: H + H + H2 <=> H2 + H2 */
    eqcon[12] *= 1e+06; 

    /*reaction 14: H + H + H2O <=> H2 + H2O */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + H + CO2 <=> H2 + CO2 */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H2 + O2 <=> HO2 + H */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + H <=> O + H2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + O <=> OH + O2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> O2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> O2 + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + H <=> OH + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: CO + OH <=> CO2 + H */
    /*eqcon[28] *= 1;  */

    /*reaction 30: CO + OH <=> CO2 + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: CO + O2 <=> CO2 + O */
    /*eqcon[30] *= 1;  */

    /*reaction 32: CO + HO2 <=> CO2 + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: HCO + H <=> CO + H2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: HCO + O <=> CO + OH */
    /*eqcon[33] *= 1;  */

    /*reaction 35: HCO + O <=> CO2 + H */
    /*eqcon[34] *= 1;  */

    /*reaction 36: HCO + OH <=> CO + H2O */
    /*eqcon[35] *= 1;  */

    /*reaction 37: HCO + H2O <=> CO + H + H2O */
    eqcon[36] *= 1e-06; 

    /*reaction 38: HCO + O2 <=> CO + HO2 */
    /*eqcon[37] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[14]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + M <=> H2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + OH + M <=> H2O + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + O + M <=> O2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: HCO + M <=> CO + H + M */
    eqcon[7] *= 1e-06; 

    /*reaction 9: H + O2 <=> O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: OH + H2 <=> H + H2O */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: H + H + H2 <=> H2 + H2 */
    eqcon[12] *= 1e+06; 

    /*reaction 14: H + H + H2O <=> H2 + H2O */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + H + CO2 <=> H2 + CO2 */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H2 + O2 <=> HO2 + H */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + H <=> O + H2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + O <=> OH + O2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> O2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> O2 + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + H <=> OH + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: CO + OH <=> CO2 + H */
    /*eqcon[28] *= 1;  */

    /*reaction 30: CO + OH <=> CO2 + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: CO + O2 <=> CO2 + O */
    /*eqcon[30] *= 1;  */

    /*reaction 32: CO + HO2 <=> CO2 + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: HCO + H <=> CO + H2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: HCO + O <=> CO + OH */
    /*eqcon[33] *= 1;  */

    /*reaction 35: HCO + O <=> CO2 + H */
    /*eqcon[34] *= 1;  */

    /*reaction 36: HCO + OH <=> CO + H2O */
    /*eqcon[35] *= 1;  */

    /*reaction 37: HCO + H2O <=> CO + H + H2O */
    eqcon[36] *= 1e-06; 

    /*reaction 38: HCO + O2 <=> CO + HO2 */
    /*eqcon[37] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[14]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + M <=> H2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + OH + M <=> H2O + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + O + M <=> O2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: HCO + M <=> CO + H + M */
    eqcon[7] *= 1e-06; 

    /*reaction 9: H + O2 <=> O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: OH + H2 <=> H + H2O */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: H + H + H2 <=> H2 + H2 */
    eqcon[12] *= 1e+06; 

    /*reaction 14: H + H + H2O <=> H2 + H2O */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + H + CO2 <=> H2 + CO2 */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H2 + O2 <=> HO2 + H */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + H <=> O + H2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + O <=> OH + O2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> O2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> O2 + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + H <=> OH + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: CO + OH <=> CO2 + H */
    /*eqcon[28] *= 1;  */

    /*reaction 30: CO + OH <=> CO2 + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: CO + O2 <=> CO2 + O */
    /*eqcon[30] *= 1;  */

    /*reaction 32: CO + HO2 <=> CO2 + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: HCO + H <=> CO + H2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: HCO + O <=> CO + OH */
    /*eqcon[33] *= 1;  */

    /*reaction 35: HCO + O <=> CO2 + H */
    /*eqcon[34] *= 1;  */

    /*reaction 36: HCO + OH <=> CO + H2O */
    /*eqcon[35] *= 1;  */

    /*reaction 37: HCO + H2O <=> CO + H + H2O */
    eqcon[36] *= 1e-06; 

    /*reaction 38: HCO + O2 <=> CO + HO2 */
    /*eqcon[37] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[14]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + M <=> H2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + OH + M <=> H2O + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + O + M <=> O2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: HCO + M <=> CO + H + M */
    eqcon[7] *= 1e-06; 

    /*reaction 9: H + O2 <=> O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: OH + H2 <=> H + H2O */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: H + H + H2 <=> H2 + H2 */
    eqcon[12] *= 1e+06; 

    /*reaction 14: H + H + H2O <=> H2 + H2O */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + H + CO2 <=> H2 + CO2 */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H2 + O2 <=> HO2 + H */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + H <=> O + H2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + O <=> OH + O2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> O2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> O2 + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + H <=> OH + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: CO + OH <=> CO2 + H */
    /*eqcon[28] *= 1;  */

    /*reaction 30: CO + OH <=> CO2 + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: CO + O2 <=> CO2 + O */
    /*eqcon[30] *= 1;  */

    /*reaction 32: CO + HO2 <=> CO2 + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: HCO + H <=> CO + H2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: HCO + O <=> CO + OH */
    /*eqcon[33] *= 1;  */

    /*reaction 35: HCO + O <=> CO2 + H */
    /*eqcon[34] *= 1;  */

    /*reaction 36: HCO + OH <=> CO + H2O */
    /*eqcon[35] *= 1;  */

    /*reaction 37: HCO + H2O <=> CO + H + H2O */
    eqcon[36] *= 1e-06; 

    /*reaction 38: HCO + O2 <=> CO + HO2 */
    /*eqcon[37] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[14]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + M <=> H2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + OH + M <=> H2O + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + O + M <=> O2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: HCO + M <=> CO + H + M */
    eqcon[7] *= 1e-06; 

    /*reaction 9: H + O2 <=> O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: OH + H2 <=> H + H2O */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: H + H + H2 <=> H2 + H2 */
    eqcon[12] *= 1e+06; 

    /*reaction 14: H + H + H2O <=> H2 + H2O */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + H + CO2 <=> H2 + CO2 */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H2 + O2 <=> HO2 + H */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + H <=> O + H2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + O <=> OH + O2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> O2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> O2 + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + H <=> OH + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: CO + OH <=> CO2 + H */
    /*eqcon[28] *= 1;  */

    /*reaction 30: CO + OH <=> CO2 + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: CO + O2 <=> CO2 + O */
    /*eqcon[30] *= 1;  */

    /*reaction 32: CO + HO2 <=> CO2 + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: HCO + H <=> CO + H2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: HCO + O <=> CO + OH */
    /*eqcon[33] *= 1;  */

    /*reaction 35: HCO + O <=> CO2 + H */
    /*eqcon[34] *= 1;  */

    /*reaction 36: HCO + OH <=> CO + H2O */
    /*eqcon[35] *= 1;  */

    /*reaction 37: HCO + H2O <=> CO + H + H2O */
    eqcon[36] *= 1e-06; 

    /*reaction 38: HCO + O2 <=> CO + HO2 */
    /*eqcon[37] *= 1;  */
}

static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[38];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[38];
#ifdef _OPENMP
#pragma omp threadprivate(Kc_save)
#endif


/*compute the production rate for each species */
void productionRate(double * restrict wdot, double * restrict sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    double qdot, q_f[38], q_r[38];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 14; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[1] -= qdot;
    wdot[8] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[1]-q_r[1];
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[2]-q_r[2];
    wdot[5] -= qdot;
    wdot[10] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[1] -= qdot;

    qdot = q_f[4]-q_r[4];
    wdot[1] -= qdot;
    wdot[6] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[1] -= qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[6]-q_r[6];
    wdot[5] -= qdot;
    wdot[5] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[1] += qdot;
    wdot[7] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[1] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[9]-q_r[9];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[6] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[12]-q_r[12];
    wdot[0] -= qdot;
    wdot[0] += qdot;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[1] -= qdot;

    qdot = q_f[13]-q_r[13];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[1] -= qdot;
    wdot[9] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[14]-q_r[14];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[1] -= qdot;
    wdot[13] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[8] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[16]-q_r[16];
    wdot[1] -= qdot;
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[17]-q_r[17];
    wdot[1] -= qdot;
    wdot[6] += qdot;
    wdot[6] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[18]-q_r[18];
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[8] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[19]-q_r[19];
    wdot[6] -= qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;
    wdot[11] += qdot;

    qdot = q_f[20]-q_r[20];
    wdot[6] -= qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;
    wdot[11] += qdot;

    qdot = q_f[21]-q_r[21];
    wdot[8] -= qdot;
    wdot[8] -= qdot;
    wdot[11] += qdot;
    wdot[12] += qdot;

    qdot = q_f[22]-q_r[22];
    wdot[8] -= qdot;
    wdot[8] -= qdot;
    wdot[11] += qdot;
    wdot[12] += qdot;

    qdot = q_f[23]-q_r[23];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[8] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[24]-q_r[24];
    wdot[1] -= qdot;
    wdot[6] += qdot;
    wdot[9] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[25]-q_r[25];
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[8] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[26]-q_r[26];
    wdot[6] -= qdot;
    wdot[8] += qdot;
    wdot[9] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[27]-q_r[27];
    wdot[6] -= qdot;
    wdot[8] += qdot;
    wdot[9] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[28]-q_r[28];
    wdot[1] += qdot;
    wdot[6] -= qdot;
    wdot[10] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[29]-q_r[29];
    wdot[1] += qdot;
    wdot[6] -= qdot;
    wdot[10] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[30]-q_r[30];
    wdot[5] += qdot;
    wdot[10] -= qdot;
    wdot[11] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[31]-q_r[31];
    wdot[6] += qdot;
    wdot[8] -= qdot;
    wdot[10] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[32]-q_r[32];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[7] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[33]-q_r[33];
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[34]-q_r[34];
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[7] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[35]-q_r[35];
    wdot[6] -= qdot;
    wdot[7] -= qdot;
    wdot[9] += qdot;
    wdot[10] += qdot;

    qdot = q_f[36]-q_r[36];
    wdot[1] += qdot;
    wdot[7] -= qdot;
    wdot[9] -= qdot;
    wdot[9] += qdot;
    wdot[10] += qdot;

    qdot = q_f[37]-q_r[37];
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[10] += qdot;
    wdot[11] -= qdot;

    return;
}

void comp_k_f(double * restrict tc, double invT, double * restrict k_f)
{
#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<38; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(double * restrict tc, double invT, double * restrict Kc)
{
    /*compute the Gibbs free energy */
    double g_RT[14];
    gibbs(g_RT, tc);

    Kc[0] = g_RT[1] - g_RT[8] + g_RT[11];
    Kc[1] = g_RT[6] + g_RT[6] - g_RT[12];
    Kc[2] = g_RT[5] + g_RT[10] - g_RT[13];
    Kc[3] = -g_RT[0] + g_RT[1] + g_RT[1];
    Kc[4] = g_RT[1] + g_RT[6] - g_RT[9];
    Kc[5] = g_RT[1] + g_RT[5] - g_RT[6];
    Kc[6] = g_RT[5] + g_RT[5] - g_RT[11];
    Kc[7] = -g_RT[1] + g_RT[7] - g_RT[10];
    Kc[8] = g_RT[1] - g_RT[5] - g_RT[6] + g_RT[11];
    Kc[9] = g_RT[0] - g_RT[1] + g_RT[5] - g_RT[6];
    Kc[10] = g_RT[0] - g_RT[1] + g_RT[6] - g_RT[9];
    Kc[11] = -g_RT[5] + g_RT[6] + g_RT[6] - g_RT[9];
    Kc[12] = g_RT[0] - g_RT[0] - g_RT[0] + g_RT[1] + g_RT[1];
    Kc[13] = -g_RT[0] + g_RT[1] + g_RT[1] + g_RT[9] - g_RT[9];
    Kc[14] = -g_RT[0] + g_RT[1] + g_RT[1] + g_RT[13] - g_RT[13];
    Kc[15] = g_RT[0] - g_RT[1] - g_RT[8] + g_RT[11];
    Kc[16] = g_RT[1] - g_RT[5] + g_RT[8] - g_RT[9];
    Kc[17] = g_RT[1] - g_RT[6] - g_RT[6] + g_RT[8];
    Kc[18] = g_RT[5] - g_RT[6] + g_RT[8] - g_RT[11];
    Kc[19] = g_RT[6] + g_RT[8] - g_RT[9] - g_RT[11];
    Kc[20] = g_RT[6] + g_RT[8] - g_RT[9] - g_RT[11];
    Kc[21] = g_RT[8] + g_RT[8] - g_RT[11] - g_RT[12];
    Kc[22] = g_RT[8] + g_RT[8] - g_RT[11] - g_RT[12];
    Kc[23] = -g_RT[0] + g_RT[1] - g_RT[8] + g_RT[12];
    Kc[24] = g_RT[1] - g_RT[6] - g_RT[9] + g_RT[12];
    Kc[25] = g_RT[5] - g_RT[6] - g_RT[8] + g_RT[12];
    Kc[26] = g_RT[6] - g_RT[8] - g_RT[9] + g_RT[12];
    Kc[27] = g_RT[6] - g_RT[8] - g_RT[9] + g_RT[12];
    Kc[28] = -g_RT[1] + g_RT[6] + g_RT[10] - g_RT[13];
    Kc[29] = -g_RT[1] + g_RT[6] + g_RT[10] - g_RT[13];
    Kc[30] = -g_RT[5] + g_RT[10] + g_RT[11] - g_RT[13];
    Kc[31] = -g_RT[6] + g_RT[8] + g_RT[10] - g_RT[13];
    Kc[32] = -g_RT[0] + g_RT[1] + g_RT[7] - g_RT[10];
    Kc[33] = g_RT[5] - g_RT[6] + g_RT[7] - g_RT[10];
    Kc[34] = -g_RT[1] + g_RT[5] + g_RT[7] - g_RT[13];
    Kc[35] = g_RT[6] + g_RT[7] - g_RT[9] - g_RT[10];
    Kc[36] = -g_RT[1] + g_RT[7] + g_RT[9] - g_RT[9] - g_RT[10];
    Kc[37] = g_RT[7] - g_RT[8] - g_RT[10] + g_RT[11];

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<38; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 * invT;
    double refCinv = 1 / refC;

    Kc[0] *= refCinv;
    Kc[1] *= refCinv;
    Kc[2] *= refCinv;
    Kc[3] *= refCinv;
    Kc[4] *= refCinv;
    Kc[5] *= refCinv;
    Kc[6] *= refCinv;
    Kc[7] *= refC;
    Kc[12] *= refCinv;
    Kc[13] *= refCinv;
    Kc[14] *= refCinv;
    Kc[36] *= refC;

    return;
}

void comp_qfqr(double * restrict qf, double * restrict qr, double * restrict sc, double * restrict tc, double invT)
{

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    qf[0] = sc[1]*sc[11];
    qr[0] = sc[8];

    /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
    qf[1] = sc[6]*sc[6];
    qr[1] = sc[12];

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    qf[2] = sc[5]*sc[10];
    qr[2] = sc[13];

    /*reaction 4: H + H + M <=> H2 + M */
    qf[3] = sc[1]*sc[1];
    qr[3] = sc[0];

    /*reaction 5: H + OH + M <=> H2O + M */
    qf[4] = sc[1]*sc[6];
    qr[4] = sc[9];

    /*reaction 6: O + H + M <=> OH + M */
    qf[5] = sc[1]*sc[5];
    qr[5] = sc[6];

    /*reaction 7: O + O + M <=> O2 + M */
    qf[6] = sc[5]*sc[5];
    qr[6] = sc[11];

    /*reaction 8: HCO + M <=> CO + H + M */
    qf[7] = sc[7];
    qr[7] = sc[1]*sc[10];

    /*reaction 9: H + O2 <=> O + OH */
    qf[8] = sc[1]*sc[11];
    qr[8] = sc[5]*sc[6];

    /*reaction 10: O + H2 <=> H + OH */
    qf[9] = sc[0]*sc[5];
    qr[9] = sc[1]*sc[6];

    /*reaction 11: OH + H2 <=> H + H2O */
    qf[10] = sc[0]*sc[6];
    qr[10] = sc[1]*sc[9];

    /*reaction 12: OH + OH <=> O + H2O */
    qf[11] = sc[6]*sc[6];
    qr[11] = sc[5]*sc[9];

    /*reaction 13: H + H + H2 <=> H2 + H2 */
    qf[12] = sc[0]*sc[1]*sc[1];
    qr[12] = sc[0]*sc[0];

    /*reaction 14: H + H + H2O <=> H2 + H2O */
    qf[13] = sc[1]*sc[1]*sc[9];
    qr[13] = sc[0]*sc[9];

    /*reaction 15: H + H + CO2 <=> H2 + CO2 */
    qf[14] = sc[1]*sc[1]*sc[13];
    qr[14] = sc[0]*sc[13];

    /*reaction 16: H2 + O2 <=> HO2 + H */
    qf[15] = sc[0]*sc[11];
    qr[15] = sc[1]*sc[8];

    /*reaction 17: HO2 + H <=> O + H2O */
    qf[16] = sc[1]*sc[8];
    qr[16] = sc[5]*sc[9];

    /*reaction 18: HO2 + H <=> OH + OH */
    qf[17] = sc[1]*sc[8];
    qr[17] = sc[6]*sc[6];

    /*reaction 19: HO2 + O <=> OH + O2 */
    qf[18] = sc[5]*sc[8];
    qr[18] = sc[6]*sc[11];

    /*reaction 20: HO2 + OH <=> O2 + H2O */
    qf[19] = sc[6]*sc[8];
    qr[19] = sc[9]*sc[11];

    /*reaction 21: HO2 + OH <=> O2 + H2O */
    qf[20] = sc[6]*sc[8];
    qr[20] = sc[9]*sc[11];

    /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
    qf[21] = sc[8]*sc[8];
    qr[21] = sc[11]*sc[12];

    /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
    qf[22] = sc[8]*sc[8];
    qr[22] = sc[11]*sc[12];

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    qf[23] = sc[1]*sc[12];
    qr[23] = sc[0]*sc[8];

    /*reaction 25: H2O2 + H <=> OH + H2O */
    qf[24] = sc[1]*sc[12];
    qr[24] = sc[6]*sc[9];

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    qf[25] = sc[5]*sc[12];
    qr[25] = sc[6]*sc[8];

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    qf[26] = sc[6]*sc[12];
    qr[26] = sc[8]*sc[9];

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    qf[27] = sc[6]*sc[12];
    qr[27] = sc[8]*sc[9];

    /*reaction 29: CO + OH <=> CO2 + H */
    qf[28] = sc[6]*sc[10];
    qr[28] = sc[1]*sc[13];

    /*reaction 30: CO + OH <=> CO2 + H */
    qf[29] = sc[6]*sc[10];
    qr[29] = sc[1]*sc[13];

    /*reaction 31: CO + O2 <=> CO2 + O */
    qf[30] = sc[10]*sc[11];
    qr[30] = sc[5]*sc[13];

    /*reaction 32: CO + HO2 <=> CO2 + OH */
    qf[31] = sc[8]*sc[10];
    qr[31] = sc[6]*sc[13];

    /*reaction 33: HCO + H <=> CO + H2 */
    qf[32] = sc[1]*sc[7];
    qr[32] = sc[0]*sc[10];

    /*reaction 34: HCO + O <=> CO + OH */
    qf[33] = sc[5]*sc[7];
    qr[33] = sc[6]*sc[10];

    /*reaction 35: HCO + O <=> CO2 + H */
    qf[34] = sc[5]*sc[7];
    qr[34] = sc[1]*sc[13];

    /*reaction 36: HCO + OH <=> CO + H2O */
    qf[35] = sc[6]*sc[7];
    qr[35] = sc[9]*sc[10];

    /*reaction 37: HCO + H2O <=> CO + H + H2O */
    qf[36] = sc[7]*sc[9];
    qr[36] = sc[1]*sc[9]*sc[10];

    /*reaction 38: HCO + O2 <=> CO + HO2 */
    qf[37] = sc[7]*sc[11];
    qr[37] = sc[8]*sc[10];

    double T = tc[1];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 14; ++i) {
        mixture += sc[i];
    }

    double Corr[38];
    for (int i = 0; i < 38; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        double alpha[2];
        alpha[0] = mixture - 0.15000000000000002*sc[11] + 10.890000000000001*sc[9] + 0.09000000000000008*sc[10] + 1.1800000000000002*sc[13] - 0.59999999999999998*sc[2] - 0.54000000000000004*sc[4] - 0.25*sc[0];
        alpha[1] = mixture + sc[0] + 5*sc[9] + 0.75*sc[10] + 2.6000000000000001*sc[13] - 0.30000000000000004*sc[2] - 0.30000000000000004*sc[4];
        for (int i=0; i<2; i++)
        {
            double redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;
            redP = alpha[i-0] / k_f_save[i] * phase_units[i] * low_A[i] * exp(low_beta[i] * tc[0] - activation_units[i] * low_Ea[i] *invT);
            F = redP / (1.0 + redP);
            logPred = log10(redP);
            logFcent = log10(
                (fabs(troe_Tsss[i]) > 1.e-100 ? (1.-troe_a[i])*exp(-T/troe_Tsss[i]) : 0.) 
                + (fabs(troe_Ts[i]) > 1.e-100 ? troe_a[i] * exp(-T/troe_Ts[i]) : 0.) 
                + (troe_len[i] == 4 ? exp(-troe_Tss[i] * invT) : 0.) );
            troe_c = -.4 - .67 * logFcent;
            troe_n = .75 - 1.27 * logFcent;
            troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
            F_troe = pow(10., logFcent / (1.0 + troe*troe));
            Corr[i] = F * F_troe;
        }
    }

    /* Lindemann */
    {
        double alpha;
        alpha = mixture + sc[0] + 11*sc[9] + 0.75*sc[10] + 2.6000000000000001*sc[13] - 0.30000000000000004*sc[2] - 0.30000000000000004*sc[4];
        double redP = alpha / k_f_save[2] * phase_units[2] * low_A[2] * exp(low_beta[2] * tc[0] - activation_units[2] * low_Ea[2] * invT);
        Corr[2] = redP / (1. + redP);
    }

    /* simple three-body correction */
    {
        double alpha;
        alpha = mixture - sc[0] - sc[9] - sc[13] - 0.37*sc[2] - 0.37*sc[4];
        Corr[3] = alpha;
        alpha = mixture + sc[0] + 5.2999999999999998*sc[9] + 0.75*sc[10] + 2.6000000000000001*sc[13] - 0.62*sc[2] - 0.62*sc[4];
        Corr[4] = alpha;
        alpha = mixture + sc[0] + 11*sc[9] + 0.75*sc[10] + 2.6000000000000001*sc[13] - 0.30000000000000004*sc[2] - 0.30000000000000004*sc[4];
        Corr[5] = alpha;
        alpha = mixture + 1.3999999999999999*sc[0] + 14.4*sc[9] + 0.75*sc[10] + 2.6000000000000001*sc[13] - 0.17000000000000004*sc[2] - 0.17000000000000004*sc[4];
        Corr[6] = alpha;
        alpha = mixture + sc[0] - sc[9] + 0.75*sc[10] + 2.6000000000000001*sc[13];
        Corr[7] = alpha;
    }

    for (int i=0; i<38; i++)
    {
        qf[i] *= Corr[i] * k_f_save[i];
        qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];
    }

    return;
}


/*compute the production rate for each species */
void vproductionRate(int npt, double * restrict wdot, double * restrict sc, double * restrict T)
{
    double k_f_s[38*npt], Kc_s[38*npt], mixture[npt], g_RT[14*npt];
    double tc[5*npt], invT[npt];

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

    for (int i=0; i<npt; i++) {
        mixture[i] = 0.0;
    }

    for (int n=0; n<14; n++) {
        for (int i=0; i<npt; i++) {
            mixture[i] += sc[n*npt+i];
            wdot[n*npt+i] = 0.0;
        }
    }

    vcomp_k_f(npt, k_f_s, tc, invT);

    vcomp_gibbs(npt, g_RT, tc);

    vcomp_Kc(npt, Kc_s, g_RT, invT);

    vcomp_wdot(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
}

void vcomp_k_f(int npt, double * restrict k_f_s, double * restrict tc, double * restrict invT)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        k_f_s[0*npt+i] = prefactor_units[0] * fwd_A[0] * exp(fwd_beta[0] * tc[i] - activation_units[0] * fwd_Ea[0] * invT[i]);
        k_f_s[1*npt+i] = prefactor_units[1] * fwd_A[1] * exp(fwd_beta[1] * tc[i] - activation_units[1] * fwd_Ea[1] * invT[i]);
        k_f_s[2*npt+i] = prefactor_units[2] * fwd_A[2] * exp(fwd_beta[2] * tc[i] - activation_units[2] * fwd_Ea[2] * invT[i]);
        k_f_s[3*npt+i] = prefactor_units[3] * fwd_A[3] * exp(fwd_beta[3] * tc[i] - activation_units[3] * fwd_Ea[3] * invT[i]);
        k_f_s[4*npt+i] = prefactor_units[4] * fwd_A[4] * exp(fwd_beta[4] * tc[i] - activation_units[4] * fwd_Ea[4] * invT[i]);
        k_f_s[5*npt+i] = prefactor_units[5] * fwd_A[5] * exp(fwd_beta[5] * tc[i] - activation_units[5] * fwd_Ea[5] * invT[i]);
        k_f_s[6*npt+i] = prefactor_units[6] * fwd_A[6] * exp(fwd_beta[6] * tc[i] - activation_units[6] * fwd_Ea[6] * invT[i]);
        k_f_s[7*npt+i] = prefactor_units[7] * fwd_A[7] * exp(fwd_beta[7] * tc[i] - activation_units[7] * fwd_Ea[7] * invT[i]);
        k_f_s[8*npt+i] = prefactor_units[8] * fwd_A[8] * exp(fwd_beta[8] * tc[i] - activation_units[8] * fwd_Ea[8] * invT[i]);
        k_f_s[9*npt+i] = prefactor_units[9] * fwd_A[9] * exp(fwd_beta[9] * tc[i] - activation_units[9] * fwd_Ea[9] * invT[i]);
        k_f_s[10*npt+i] = prefactor_units[10] * fwd_A[10] * exp(fwd_beta[10] * tc[i] - activation_units[10] * fwd_Ea[10] * invT[i]);
        k_f_s[11*npt+i] = prefactor_units[11] * fwd_A[11] * exp(fwd_beta[11] * tc[i] - activation_units[11] * fwd_Ea[11] * invT[i]);
        k_f_s[12*npt+i] = prefactor_units[12] * fwd_A[12] * exp(fwd_beta[12] * tc[i] - activation_units[12] * fwd_Ea[12] * invT[i]);
        k_f_s[13*npt+i] = prefactor_units[13] * fwd_A[13] * exp(fwd_beta[13] * tc[i] - activation_units[13] * fwd_Ea[13] * invT[i]);
        k_f_s[14*npt+i] = prefactor_units[14] * fwd_A[14] * exp(fwd_beta[14] * tc[i] - activation_units[14] * fwd_Ea[14] * invT[i]);
        k_f_s[15*npt+i] = prefactor_units[15] * fwd_A[15] * exp(fwd_beta[15] * tc[i] - activation_units[15] * fwd_Ea[15] * invT[i]);
        k_f_s[16*npt+i] = prefactor_units[16] * fwd_A[16] * exp(fwd_beta[16] * tc[i] - activation_units[16] * fwd_Ea[16] * invT[i]);
        k_f_s[17*npt+i] = prefactor_units[17] * fwd_A[17] * exp(fwd_beta[17] * tc[i] - activation_units[17] * fwd_Ea[17] * invT[i]);
        k_f_s[18*npt+i] = prefactor_units[18] * fwd_A[18] * exp(fwd_beta[18] * tc[i] - activation_units[18] * fwd_Ea[18] * invT[i]);
        k_f_s[19*npt+i] = prefactor_units[19] * fwd_A[19] * exp(fwd_beta[19] * tc[i] - activation_units[19] * fwd_Ea[19] * invT[i]);
        k_f_s[20*npt+i] = prefactor_units[20] * fwd_A[20] * exp(fwd_beta[20] * tc[i] - activation_units[20] * fwd_Ea[20] * invT[i]);
        k_f_s[21*npt+i] = prefactor_units[21] * fwd_A[21] * exp(fwd_beta[21] * tc[i] - activation_units[21] * fwd_Ea[21] * invT[i]);
        k_f_s[22*npt+i] = prefactor_units[22] * fwd_A[22] * exp(fwd_beta[22] * tc[i] - activation_units[22] * fwd_Ea[22] * invT[i]);
        k_f_s[23*npt+i] = prefactor_units[23] * fwd_A[23] * exp(fwd_beta[23] * tc[i] - activation_units[23] * fwd_Ea[23] * invT[i]);
        k_f_s[24*npt+i] = prefactor_units[24] * fwd_A[24] * exp(fwd_beta[24] * tc[i] - activation_units[24] * fwd_Ea[24] * invT[i]);
        k_f_s[25*npt+i] = prefactor_units[25] * fwd_A[25] * exp(fwd_beta[25] * tc[i] - activation_units[25] * fwd_Ea[25] * invT[i]);
        k_f_s[26*npt+i] = prefactor_units[26] * fwd_A[26] * exp(fwd_beta[26] * tc[i] - activation_units[26] * fwd_Ea[26] * invT[i]);
        k_f_s[27*npt+i] = prefactor_units[27] * fwd_A[27] * exp(fwd_beta[27] * tc[i] - activation_units[27] * fwd_Ea[27] * invT[i]);
        k_f_s[28*npt+i] = prefactor_units[28] * fwd_A[28] * exp(fwd_beta[28] * tc[i] - activation_units[28] * fwd_Ea[28] * invT[i]);
        k_f_s[29*npt+i] = prefactor_units[29] * fwd_A[29] * exp(fwd_beta[29] * tc[i] - activation_units[29] * fwd_Ea[29] * invT[i]);
        k_f_s[30*npt+i] = prefactor_units[30] * fwd_A[30] * exp(fwd_beta[30] * tc[i] - activation_units[30] * fwd_Ea[30] * invT[i]);
        k_f_s[31*npt+i] = prefactor_units[31] * fwd_A[31] * exp(fwd_beta[31] * tc[i] - activation_units[31] * fwd_Ea[31] * invT[i]);
        k_f_s[32*npt+i] = prefactor_units[32] * fwd_A[32] * exp(fwd_beta[32] * tc[i] - activation_units[32] * fwd_Ea[32] * invT[i]);
        k_f_s[33*npt+i] = prefactor_units[33] * fwd_A[33] * exp(fwd_beta[33] * tc[i] - activation_units[33] * fwd_Ea[33] * invT[i]);
        k_f_s[34*npt+i] = prefactor_units[34] * fwd_A[34] * exp(fwd_beta[34] * tc[i] - activation_units[34] * fwd_Ea[34] * invT[i]);
        k_f_s[35*npt+i] = prefactor_units[35] * fwd_A[35] * exp(fwd_beta[35] * tc[i] - activation_units[35] * fwd_Ea[35] * invT[i]);
        k_f_s[36*npt+i] = prefactor_units[36] * fwd_A[36] * exp(fwd_beta[36] * tc[i] - activation_units[36] * fwd_Ea[36] * invT[i]);
        k_f_s[37*npt+i] = prefactor_units[37] * fwd_A[37] * exp(fwd_beta[37] * tc[i] - activation_units[37] * fwd_Ea[37] * invT[i]);
    }
}

void vcomp_gibbs(int npt, double * restrict g_RT, double * restrict tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        double tg[5], g[14];
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
        g_RT[9*npt+i] = g[9];
        g_RT[10*npt+i] = g[10];
        g_RT[11*npt+i] = g[11];
        g_RT[12*npt+i] = g[12];
        g_RT[13*npt+i] = g[13];
    }
}

void vcomp_Kc(int npt, double * restrict Kc_s, double * restrict g_RT, double * restrict invT)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        double refC = (101325. / 8.31451) * invT[i];
        double refCinv = 1.0 / refC;

        Kc_s[0*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[11*npt+i]) - (g_RT[8*npt+i]));
        Kc_s[1*npt+i] = refCinv * exp((g_RT[6*npt+i] + g_RT[6*npt+i]) - (g_RT[12*npt+i]));
        Kc_s[2*npt+i] = refCinv * exp((g_RT[5*npt+i] + g_RT[10*npt+i]) - (g_RT[13*npt+i]));
        Kc_s[3*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[1*npt+i]) - (g_RT[0*npt+i]));
        Kc_s[4*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[6*npt+i]) - (g_RT[9*npt+i]));
        Kc_s[5*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[5*npt+i]) - (g_RT[6*npt+i]));
        Kc_s[6*npt+i] = refCinv * exp((g_RT[5*npt+i] + g_RT[5*npt+i]) - (g_RT[11*npt+i]));
        Kc_s[7*npt+i] = refC * exp((g_RT[7*npt+i]) - (g_RT[1*npt+i] + g_RT[10*npt+i]));
        Kc_s[8*npt+i] = exp((g_RT[1*npt+i] + g_RT[11*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[9*npt+i] = exp((g_RT[0*npt+i] + g_RT[5*npt+i]) - (g_RT[1*npt+i] + g_RT[6*npt+i]));
        Kc_s[10*npt+i] = exp((g_RT[0*npt+i] + g_RT[6*npt+i]) - (g_RT[1*npt+i] + g_RT[9*npt+i]));
        Kc_s[11*npt+i] = exp((g_RT[6*npt+i] + g_RT[6*npt+i]) - (g_RT[5*npt+i] + g_RT[9*npt+i]));
        Kc_s[12*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[1*npt+i] + g_RT[1*npt+i]) - (g_RT[0*npt+i] + g_RT[0*npt+i]));
        Kc_s[13*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[1*npt+i] + g_RT[9*npt+i]) - (g_RT[0*npt+i] + g_RT[9*npt+i]));
        Kc_s[14*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[1*npt+i] + g_RT[13*npt+i]) - (g_RT[0*npt+i] + g_RT[13*npt+i]));
        Kc_s[15*npt+i] = exp((g_RT[0*npt+i] + g_RT[11*npt+i]) - (g_RT[1*npt+i] + g_RT[8*npt+i]));
        Kc_s[16*npt+i] = exp((g_RT[1*npt+i] + g_RT[8*npt+i]) - (g_RT[5*npt+i] + g_RT[9*npt+i]));
        Kc_s[17*npt+i] = exp((g_RT[1*npt+i] + g_RT[8*npt+i]) - (g_RT[6*npt+i] + g_RT[6*npt+i]));
        Kc_s[18*npt+i] = exp((g_RT[5*npt+i] + g_RT[8*npt+i]) - (g_RT[6*npt+i] + g_RT[11*npt+i]));
        Kc_s[19*npt+i] = exp((g_RT[6*npt+i] + g_RT[8*npt+i]) - (g_RT[9*npt+i] + g_RT[11*npt+i]));
        Kc_s[20*npt+i] = exp((g_RT[6*npt+i] + g_RT[8*npt+i]) - (g_RT[9*npt+i] + g_RT[11*npt+i]));
        Kc_s[21*npt+i] = exp((g_RT[8*npt+i] + g_RT[8*npt+i]) - (g_RT[11*npt+i] + g_RT[12*npt+i]));
        Kc_s[22*npt+i] = exp((g_RT[8*npt+i] + g_RT[8*npt+i]) - (g_RT[11*npt+i] + g_RT[12*npt+i]));
        Kc_s[23*npt+i] = exp((g_RT[1*npt+i] + g_RT[12*npt+i]) - (g_RT[0*npt+i] + g_RT[8*npt+i]));
        Kc_s[24*npt+i] = exp((g_RT[1*npt+i] + g_RT[12*npt+i]) - (g_RT[6*npt+i] + g_RT[9*npt+i]));
        Kc_s[25*npt+i] = exp((g_RT[5*npt+i] + g_RT[12*npt+i]) - (g_RT[6*npt+i] + g_RT[8*npt+i]));
        Kc_s[26*npt+i] = exp((g_RT[6*npt+i] + g_RT[12*npt+i]) - (g_RT[8*npt+i] + g_RT[9*npt+i]));
        Kc_s[27*npt+i] = exp((g_RT[6*npt+i] + g_RT[12*npt+i]) - (g_RT[8*npt+i] + g_RT[9*npt+i]));
        Kc_s[28*npt+i] = exp((g_RT[6*npt+i] + g_RT[10*npt+i]) - (g_RT[1*npt+i] + g_RT[13*npt+i]));
        Kc_s[29*npt+i] = exp((g_RT[6*npt+i] + g_RT[10*npt+i]) - (g_RT[1*npt+i] + g_RT[13*npt+i]));
        Kc_s[30*npt+i] = exp((g_RT[10*npt+i] + g_RT[11*npt+i]) - (g_RT[5*npt+i] + g_RT[13*npt+i]));
        Kc_s[31*npt+i] = exp((g_RT[8*npt+i] + g_RT[10*npt+i]) - (g_RT[6*npt+i] + g_RT[13*npt+i]));
        Kc_s[32*npt+i] = exp((g_RT[1*npt+i] + g_RT[7*npt+i]) - (g_RT[0*npt+i] + g_RT[10*npt+i]));
        Kc_s[33*npt+i] = exp((g_RT[5*npt+i] + g_RT[7*npt+i]) - (g_RT[6*npt+i] + g_RT[10*npt+i]));
        Kc_s[34*npt+i] = exp((g_RT[5*npt+i] + g_RT[7*npt+i]) - (g_RT[1*npt+i] + g_RT[13*npt+i]));
        Kc_s[35*npt+i] = exp((g_RT[6*npt+i] + g_RT[7*npt+i]) - (g_RT[9*npt+i] + g_RT[10*npt+i]));
        Kc_s[36*npt+i] = refC * exp((g_RT[7*npt+i] + g_RT[9*npt+i]) - (g_RT[1*npt+i] + g_RT[9*npt+i] + g_RT[10*npt+i]));
        Kc_s[37*npt+i] = exp((g_RT[7*npt+i] + g_RT[11*npt+i]) - (g_RT[8*npt+i] + g_RT[10*npt+i]));
    }
}

void vcomp_wdot(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,
		double * restrict k_f_s, double * restrict Kc_s,
		double * restrict tc, double * restrict invT, double * restrict T)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        double alpha;
        double redP, F;
        double logPred;
        double logFcent, troe_c, troe_n, troe, F_troe;

        /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
        phi_f = sc[1*npt+i]*sc[11*npt+i];
        alpha = mixture[i] + -0.15000000000000002*sc[11*npt+i] + 10.890000000000001*sc[9*npt+i] + 0.09000000000000008*sc[10*npt+i] + 1.1800000000000002*sc[13*npt+i] + -0.59999999999999998*sc[2*npt+i] + -0.54000000000000004*sc[4*npt+i] + -0.25*sc[0*npt+i];
        k_f = k_f_s[0*npt+i];
        redP = alpha / k_f * phase_units[0] * low_A[0] * exp(low_beta[0] * tc[i] - activation_units[0] * low_Ea[0] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[0]) > 1.e-100 ? (1.-troe_a[0])*exp(-T[i]/troe_Tsss[0]) : 0.) 
            + (fabs(troe_Ts[0]) > 1.e-100 ? troe_a[0] * exp(-T[i]/troe_Ts[0]) : 0.) 
            + (troe_len[0] == 4 ? exp(-troe_Tss[0] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i];
        Kc = Kc_s[0*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        alpha = mixture[i] + sc[0*npt+i] + 5*sc[9*npt+i] + 0.75*sc[10*npt+i] + 2.6000000000000001*sc[13*npt+i] + -0.30000000000000004*sc[2*npt+i] + -0.30000000000000004*sc[4*npt+i];
        k_f = k_f_s[1*npt+i];
        redP = alpha / k_f * phase_units[1] * low_A[1] * exp(low_beta[1] * tc[i] - activation_units[1] * low_Ea[1] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[1]) > 1.e-100 ? (1.-troe_a[1])*exp(-T[i]/troe_Tsss[1]) : 0.) 
            + (fabs(troe_Ts[1]) > 1.e-100 ? troe_a[1] * exp(-T[i]/troe_Ts[1]) : 0.) 
            + (troe_len[1] == 4 ? exp(-troe_Tss[1] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[12*npt+i];
        Kc = Kc_s[1*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 3: CO + O (+M) <=> CO2 (+M) */
        phi_f = sc[5*npt+i]*sc[10*npt+i];
        alpha = mixture[i] + sc[0*npt+i] + 11*sc[9*npt+i] + 0.75*sc[10*npt+i] + 2.6000000000000001*sc[13*npt+i] + -0.30000000000000004*sc[2*npt+i] + -0.30000000000000004*sc[4*npt+i];
        k_f = k_f_s[2*npt+i];
        redP = alpha / k_f * phase_units[2] * low_A[2] * exp(low_beta[2] * tc[i] - activation_units[2] * low_Ea[2] * invT[i]);
        F = redP / (1 + redP);
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i];
        Kc = Kc_s[2*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 4: H + H + M <=> H2 + M */
        phi_f = sc[1*npt+i]*sc[1*npt+i];
        alpha = mixture[i] + -sc[0*npt+i] + -sc[9*npt+i] + -sc[13*npt+i] + -0.37*sc[2*npt+i] + -0.37*sc[4*npt+i];
        k_f = alpha * k_f_s[3*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[1*npt+i] -= qdot;

        /*reaction 5: H + OH + M <=> H2O + M */
        phi_f = sc[1*npt+i]*sc[6*npt+i];
        alpha = mixture[i] + sc[0*npt+i] + 5.2999999999999998*sc[9*npt+i] + 0.75*sc[10*npt+i] + 2.6000000000000001*sc[13*npt+i] + -0.62*sc[2*npt+i] + -0.62*sc[4*npt+i];
        k_f = alpha * k_f_s[4*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i];
        Kc = Kc_s[4*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 6: O + H + M <=> OH + M */
        phi_f = sc[1*npt+i]*sc[5*npt+i];
        alpha = mixture[i] + sc[0*npt+i] + 11*sc[9*npt+i] + 0.75*sc[10*npt+i] + 2.6000000000000001*sc[13*npt+i] + -0.30000000000000004*sc[2*npt+i] + -0.30000000000000004*sc[4*npt+i];
        k_f = alpha * k_f_s[5*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i];
        Kc = Kc_s[5*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 7: O + O + M <=> O2 + M */
        phi_f = sc[5*npt+i]*sc[5*npt+i];
        alpha = mixture[i] + 1.3999999999999999*sc[0*npt+i] + 14.4*sc[9*npt+i] + 0.75*sc[10*npt+i] + 2.6000000000000001*sc[13*npt+i] + -0.17000000000000004*sc[2*npt+i] + -0.17000000000000004*sc[4*npt+i];
        k_f = alpha * k_f_s[6*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i];
        Kc = Kc_s[6*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 8: HCO + M <=> CO + H + M */
        phi_f = sc[7*npt+i];
        alpha = mixture[i] + sc[0*npt+i] + -sc[9*npt+i] + 0.75*sc[10*npt+i] + 2.6000000000000001*sc[13*npt+i];
        k_f = alpha * k_f_s[7*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[10*npt+i];
        Kc = Kc_s[7*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 9: H + O2 <=> O + OH */
        phi_f = sc[1*npt+i]*sc[11*npt+i];
        k_f = k_f_s[8*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[8*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 10: O + H2 <=> H + OH */
        phi_f = sc[0*npt+i]*sc[5*npt+i];
        k_f = k_f_s[9*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[6*npt+i];
        Kc = Kc_s[9*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 11: OH + H2 <=> H + H2O */
        phi_f = sc[0*npt+i]*sc[6*npt+i];
        k_f = k_f_s[10*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[9*npt+i];
        Kc = Kc_s[10*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 12: OH + OH <=> O + H2O */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[11*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[9*npt+i];
        Kc = Kc_s[11*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 13: H + H + H2 <=> H2 + H2 */
        phi_f = sc[0*npt+i]*sc[1*npt+i]*sc[1*npt+i];
        k_f = k_f_s[12*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[0*npt+i];
        Kc = Kc_s[12*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[0*npt+i] += qdot;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[1*npt+i] -= qdot;

        /*reaction 14: H + H + H2O <=> H2 + H2O */
        phi_f = sc[1*npt+i]*sc[1*npt+i]*sc[9*npt+i];
        k_f = k_f_s[13*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[9*npt+i];
        Kc = Kc_s[13*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[1*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 15: H + H + CO2 <=> H2 + CO2 */
        phi_f = sc[1*npt+i]*sc[1*npt+i]*sc[13*npt+i];
        k_f = k_f_s[14*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[13*npt+i];
        Kc = Kc_s[14*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[1*npt+i] -= qdot;
        wdot[13*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 16: H2 + O2 <=> HO2 + H */
        phi_f = sc[0*npt+i]*sc[11*npt+i];
        k_f = k_f_s[15*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[8*npt+i];
        Kc = Kc_s[15*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 17: HO2 + H <=> O + H2O */
        phi_f = sc[1*npt+i]*sc[8*npt+i];
        k_f = k_f_s[16*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[9*npt+i];
        Kc = Kc_s[16*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 18: HO2 + H <=> OH + OH */
        phi_f = sc[1*npt+i]*sc[8*npt+i];
        k_f = k_f_s[17*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[6*npt+i];
        Kc = Kc_s[17*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;

        /*reaction 19: HO2 + O <=> OH + O2 */
        phi_f = sc[5*npt+i]*sc[8*npt+i];
        k_f = k_f_s[18*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[11*npt+i];
        Kc = Kc_s[18*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 20: HO2 + OH <=> O2 + H2O */
        phi_f = sc[6*npt+i]*sc[8*npt+i];
        k_f = k_f_s[19*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[11*npt+i];
        Kc = Kc_s[19*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 21: HO2 + OH <=> O2 + H2O */
        phi_f = sc[6*npt+i]*sc[8*npt+i];
        k_f = k_f_s[20*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[11*npt+i];
        Kc = Kc_s[20*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
        phi_f = sc[8*npt+i]*sc[8*npt+i];
        k_f = k_f_s[21*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[12*npt+i];
        Kc = Kc_s[21*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[8*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
        phi_f = sc[8*npt+i]*sc[8*npt+i];
        k_f = k_f_s[22*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[12*npt+i];
        Kc = Kc_s[22*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[8*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 24: H2O2 + H <=> HO2 + H2 */
        phi_f = sc[1*npt+i]*sc[12*npt+i];
        k_f = k_f_s[23*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[8*npt+i];
        Kc = Kc_s[23*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 25: H2O2 + H <=> OH + H2O */
        phi_f = sc[1*npt+i]*sc[12*npt+i];
        k_f = k_f_s[24*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[9*npt+i];
        Kc = Kc_s[24*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 26: H2O2 + O <=> OH + HO2 */
        phi_f = sc[5*npt+i]*sc[12*npt+i];
        k_f = k_f_s[25*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[8*npt+i];
        Kc = Kc_s[25*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 27: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[6*npt+i]*sc[12*npt+i];
        k_f = k_f_s[26*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[9*npt+i];
        Kc = Kc_s[26*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 28: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[6*npt+i]*sc[12*npt+i];
        k_f = k_f_s[27*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[9*npt+i];
        Kc = Kc_s[27*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 29: CO + OH <=> CO2 + H */
        phi_f = sc[6*npt+i]*sc[10*npt+i];
        k_f = k_f_s[28*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[13*npt+i];
        Kc = Kc_s[28*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 30: CO + OH <=> CO2 + H */
        phi_f = sc[6*npt+i]*sc[10*npt+i];
        k_f = k_f_s[29*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[13*npt+i];
        Kc = Kc_s[29*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 31: CO + O2 <=> CO2 + O */
        phi_f = sc[10*npt+i]*sc[11*npt+i];
        k_f = k_f_s[30*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[13*npt+i];
        Kc = Kc_s[30*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 32: CO + HO2 <=> CO2 + OH */
        phi_f = sc[8*npt+i]*sc[10*npt+i];
        k_f = k_f_s[31*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[13*npt+i];
        Kc = Kc_s[31*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 33: HCO + H <=> CO + H2 */
        phi_f = sc[1*npt+i]*sc[7*npt+i];
        k_f = k_f_s[32*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[10*npt+i];
        Kc = Kc_s[32*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 34: HCO + O <=> CO + OH */
        phi_f = sc[5*npt+i]*sc[7*npt+i];
        k_f = k_f_s[33*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[10*npt+i];
        Kc = Kc_s[33*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 35: HCO + O <=> CO2 + H */
        phi_f = sc[5*npt+i]*sc[7*npt+i];
        k_f = k_f_s[34*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[13*npt+i];
        Kc = Kc_s[34*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 36: HCO + OH <=> CO + H2O */
        phi_f = sc[6*npt+i]*sc[7*npt+i];
        k_f = k_f_s[35*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[10*npt+i];
        Kc = Kc_s[35*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 37: HCO + H2O <=> CO + H + H2O */
        phi_f = sc[7*npt+i]*sc[9*npt+i];
        k_f = k_f_s[36*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[9*npt+i]*sc[10*npt+i];
        Kc = Kc_s[36*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 38: HCO + O2 <=> CO + HO2 */
        phi_f = sc[7*npt+i]*sc[11*npt+i];
        k_f = k_f_s[37*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[10*npt+i];
        Kc = Kc_s[37*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
    }
}

/*compute the reaction Jacobian */
void DWDOT(double * restrict J, double * restrict sc, double * restrict Tp, int * consP)
{
    double c[14];

    for (int k=0; k<14; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    for (int k=0; k<14; k++) {
        J[210+k] *= 1.e-6;
    }

    /* dTdot/d[X] */
    for (int k=0; k<14; k++) {
        J[k*15+14] *= 1.e6;
    }

    return;
}

/*compute the reaction Jacobian */
void aJacobian(double * restrict J, double * restrict sc, double T, int consP)
{
    for (int i=0; i<225; i++) {
        J[i] = 0.0;
    }

    double wdot[14];
    for (int k=0; k<14; k++) {
        wdot[k] = 0.0;
    }

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];
    double invT2 = invT * invT;

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;
    double refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int k = 0; k < 14; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[14];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[14];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[14];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture - 0.15000000000000002*sc[11] + 10.890000000000001*sc[9] + 0.09000000000000008*sc[10] + 1.1800000000000002*sc[13] - 0.59999999999999998*sc[2] - 0.54000000000000004*sc[4] - 0.25*sc[0];
    /* forward */
    phi_f = sc[1]*sc[11];
    k_f = prefactor_units[0] * fwd_A[0]
                * exp(fwd_beta[0] * tc[0] - activation_units[0] * fwd_Ea[0] * invT);
    dlnkfdT = fwd_beta[0] * invT + activation_units[0] * fwd_Ea[0] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[0] * exp(low_beta[0] * tc[0] - activation_units[0] * low_Ea[0] * invT);
    Pr = phase_units[0] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[0] * invT + activation_units[0] * low_Ea[0] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[0]) > 1.e-100 ? (1.-troe_a[0])*exp(-T/troe_Tsss[0]) : 0.);
    Fcent2 = (fabs(troe_Ts[0]) > 1.e-100 ? troe_a[0] * exp(-T/troe_Ts[0]) : 0.);
    Fcent3 = (troe_len[0] == 4 ? exp(-troe_Tss[0] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[0]) > 1.e-100 ? -Fcent1/troe_Tsss[0] : 0.)
      + (fabs(troe_Ts[0]) > 1.e-100 ? -Fcent2/troe_Ts[0] : 0.)
      + (troe_len[0] == 4 ? Fcent3*troe_Tss[0]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[8];
    Kc = refCinv * exp(g_RT[1] - g_RT[8] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[11]) + (h_RT[8]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[8] += q; /* HO2 */
    wdot[11] -= q; /* O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci =  -0.25*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[8] += dqdci;                /* dwdot[HO2]/d[H2] */
        J[11] -= dqdci;               /* dwdot[O2]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[11];
        J[16] -= dqdci;               /* dwdot[H]/d[H] */
        J[23] += dqdci;               /* dwdot[HO2]/d[H] */
        J[26] -= dqdci;               /* dwdot[O2]/d[H] */
        /* d()/d[AR] */
        dqdci =  -0.59999999999999998*dcdc_fac;
        J[31] -= dqdci;               /* dwdot[H]/d[AR] */
        J[38] += dqdci;               /* dwdot[HO2]/d[AR] */
        J[41] -= dqdci;               /* dwdot[O2]/d[AR] */
        /* d()/d[HE] */
        dqdci =  -0.54000000000000004*dcdc_fac;
        J[61] -= dqdci;               /* dwdot[H]/d[HE] */
        J[68] += dqdci;               /* dwdot[HO2]/d[HE] */
        J[71] -= dqdci;               /* dwdot[O2]/d[HE] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[121] -= dqdci;              /* dwdot[H]/d[HO2] */
        J[128] += dqdci;              /* dwdot[HO2]/d[HO2] */
        J[131] -= dqdci;              /* dwdot[O2]/d[HO2] */
        /* d()/d[H2O] */
        dqdci =  10.890000000000001*dcdc_fac;
        J[136] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[143] += dqdci;              /* dwdot[HO2]/d[H2O] */
        J[146] -= dqdci;              /* dwdot[O2]/d[H2O] */
        /* d()/d[CO] */
        dqdci =  0.09000000000000008*dcdc_fac;
        J[151] -= dqdci;              /* dwdot[H]/d[CO] */
        J[158] += dqdci;              /* dwdot[HO2]/d[CO] */
        J[161] -= dqdci;              /* dwdot[O2]/d[CO] */
        /* d()/d[O2] */
        dqdci =  -0.15000000000000002*dcdc_fac + k_f*sc[1];
        J[166] -= dqdci;              /* dwdot[H]/d[O2] */
        J[173] += dqdci;              /* dwdot[HO2]/d[O2] */
        J[176] -= dqdci;              /* dwdot[O2]/d[O2] */
        /* d()/d[CO2] */
        dqdci =  1.1800000000000002*dcdc_fac;
        J[196] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[203] += dqdci;              /* dwdot[HO2]/d[CO2] */
        J[206] -= dqdci;              /* dwdot[O2]/d[CO2] */
    }
    else {
        dqdc[0] =  0.75*dcdc_fac;
        dqdc[1] =  dcdc_fac + k_f*sc[11];
        dqdc[2] =  0.40000000000000002*dcdc_fac;
        dqdc[3] =  dcdc_fac;
        dqdc[4] =  0.46000000000000002*dcdc_fac;
        dqdc[5] =  dcdc_fac;
        dqdc[6] =  dcdc_fac;
        dqdc[7] =  dcdc_fac;
        dqdc[8] =  dcdc_fac - k_r;
        dqdc[9] =  11.890000000000001*dcdc_fac;
        dqdc[10] =  1.0900000000000001*dcdc_fac;
        dqdc[11] =  0.84999999999999998*dcdc_fac + k_f*sc[1];
        dqdc[12] =  dcdc_fac;
        dqdc[13] =  2.1800000000000002*dcdc_fac;
        for (int k=0; k<14; k++) {
            J[15*k+1] -= dqdc[k];
            J[15*k+8] += dqdc[k];
            J[15*k+11] -= dqdc[k];
        }
    }
    J[211] -= dqdT; /* dwdot[H]/dT */
    J[218] += dqdT; /* dwdot[HO2]/dT */
    J[221] -= dqdT; /* dwdot[O2]/dT */

    /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + sc[0] + 5*sc[9] + 0.75*sc[10] + 2.6000000000000001*sc[13] - 0.30000000000000004*sc[2] - 0.30000000000000004*sc[4];
    /* forward */
    phi_f = sc[6]*sc[6];
    k_f = prefactor_units[1] * fwd_A[1]
                * exp(fwd_beta[1] * tc[0] - activation_units[1] * fwd_Ea[1] * invT);
    dlnkfdT = fwd_beta[1] * invT + activation_units[1] * fwd_Ea[1] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[1] * exp(low_beta[1] * tc[0] - activation_units[1] * low_Ea[1] * invT);
    Pr = phase_units[1] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[1] * invT + activation_units[1] * low_Ea[1] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[1]) > 1.e-100 ? (1.-troe_a[1])*exp(-T/troe_Tsss[1]) : 0.);
    Fcent2 = (fabs(troe_Ts[1]) > 1.e-100 ? troe_a[1] * exp(-T/troe_Ts[1]) : 0.);
    Fcent3 = (troe_len[1] == 4 ? exp(-troe_Tss[1] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[1]) > 1.e-100 ? -Fcent1/troe_Tsss[1] : 0.)
      + (fabs(troe_Ts[1]) > 1.e-100 ? -Fcent2/troe_Ts[1] : 0.)
      + (troe_len[1] == 4 ? Fcent3*troe_Tss[1]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[12];
    Kc = refCinv * exp(g_RT[6] + g_RT[6] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[6]) + (h_RT[12]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[6] -= 2 * q; /* OH */
    wdot[12] += q; /* H2O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci =  dcdc_fac;
        J[6] += -2 * dqdci;           /* dwdot[OH]/d[H2] */
        J[12] += dqdci;               /* dwdot[H2O2]/d[H2] */
        /* d()/d[AR] */
        dqdci =  -0.30000000000000004*dcdc_fac;
        J[36] += -2 * dqdci;          /* dwdot[OH]/d[AR] */
        J[42] += dqdci;               /* dwdot[H2O2]/d[AR] */
        /* d()/d[HE] */
        dqdci =  -0.30000000000000004*dcdc_fac;
        J[66] += -2 * dqdci;          /* dwdot[OH]/d[HE] */
        J[72] += dqdci;               /* dwdot[H2O2]/d[HE] */
        /* d()/d[OH] */
        dqdci =  + k_f*2*sc[6];
        J[96] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
        J[102] += dqdci;              /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O] */
        dqdci =  5*dcdc_fac;
        J[141] += -2 * dqdci;         /* dwdot[OH]/d[H2O] */
        J[147] += dqdci;              /* dwdot[H2O2]/d[H2O] */
        /* d()/d[CO] */
        dqdci =  0.75*dcdc_fac;
        J[156] += -2 * dqdci;         /* dwdot[OH]/d[CO] */
        J[162] += dqdci;              /* dwdot[H2O2]/d[CO] */
        /* d()/d[H2O2] */
        dqdci =  - k_r;
        J[186] += -2 * dqdci;         /* dwdot[OH]/d[H2O2] */
        J[192] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
        /* d()/d[CO2] */
        dqdci =  2.6000000000000001*dcdc_fac;
        J[201] += -2 * dqdci;         /* dwdot[OH]/d[CO2] */
        J[207] += dqdci;              /* dwdot[H2O2]/d[CO2] */
    }
    else {
        dqdc[0] =  2*dcdc_fac;
        dqdc[1] =  dcdc_fac;
        dqdc[2] =  0.69999999999999996*dcdc_fac;
        dqdc[3] =  dcdc_fac;
        dqdc[4] =  0.69999999999999996*dcdc_fac;
        dqdc[5] =  dcdc_fac;
        dqdc[6] =  dcdc_fac + k_f*2*sc[6];
        dqdc[7] =  dcdc_fac;
        dqdc[8] =  dcdc_fac;
        dqdc[9] =  6*dcdc_fac;
        dqdc[10] =  1.75*dcdc_fac;
        dqdc[11] =  dcdc_fac;
        dqdc[12] =  dcdc_fac - k_r;
        dqdc[13] =  3.6000000000000001*dcdc_fac;
        for (int k=0; k<14; k++) {
            J[15*k+6] += -2 * dqdc[k];
            J[15*k+12] += dqdc[k];
        }
    }
    J[216] += -2 * dqdT; /* dwdot[OH]/dT */
    J[222] += dqdT; /* dwdot[H2O2]/dT */

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + sc[0] + 11*sc[9] + 0.75*sc[10] + 2.6000000000000001*sc[13] - 0.30000000000000004*sc[2] - 0.30000000000000004*sc[4];
    /* forward */
    phi_f = sc[5]*sc[10];
    k_f = prefactor_units[2] * fwd_A[2]
                * exp(fwd_beta[2] * tc[0] - activation_units[2] * fwd_Ea[2] * invT);
    dlnkfdT = fwd_beta[2] * invT + activation_units[2] * fwd_Ea[2] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[2] * exp(low_beta[2] * tc[0] - activation_units[2] * low_Ea[2] * invT);
    Pr = phase_units[2] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[2] * invT + activation_units[2] * low_Ea[2] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Lindemann form */
    F = 1.0;
    dlogFdlogPr = 0.0;
    dlogFdT = 0.0;
    /* reverse */
    phi_r = sc[13];
    Kc = refCinv * exp(g_RT[5] + g_RT[10] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[10]) + (h_RT[13]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[10] -= q; /* CO */
    wdot[13] += q; /* CO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci =  dcdc_fac;
        J[5] -= dqdci;                /* dwdot[O]/d[H2] */
        J[10] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[13] += dqdci;               /* dwdot[CO2]/d[H2] */
        /* d()/d[AR] */
        dqdci =  -0.30000000000000004*dcdc_fac;
        J[35] -= dqdci;               /* dwdot[O]/d[AR] */
        J[40] -= dqdci;               /* dwdot[CO]/d[AR] */
        J[43] += dqdci;               /* dwdot[CO2]/d[AR] */
        /* d()/d[HE] */
        dqdci =  -0.30000000000000004*dcdc_fac;
        J[65] -= dqdci;               /* dwdot[O]/d[HE] */
        J[70] -= dqdci;               /* dwdot[CO]/d[HE] */
        J[73] += dqdci;               /* dwdot[CO2]/d[HE] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[10];
        J[80] -= dqdci;               /* dwdot[O]/d[O] */
        J[85] -= dqdci;               /* dwdot[CO]/d[O] */
        J[88] += dqdci;               /* dwdot[CO2]/d[O] */
        /* d()/d[H2O] */
        dqdci =  11*dcdc_fac;
        J[140] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[145] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[148] += dqdci;              /* dwdot[CO2]/d[H2O] */
        /* d()/d[CO] */
        dqdci =  0.75*dcdc_fac + k_f*sc[5];
        J[155] -= dqdci;              /* dwdot[O]/d[CO] */
        J[160] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[163] += dqdci;              /* dwdot[CO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci =  2.6000000000000001*dcdc_fac - k_r;
        J[200] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[205] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[208] += dqdci;              /* dwdot[CO2]/d[CO2] */
    }
    else {
        dqdc[0] =  2*dcdc_fac;
        dqdc[1] =  dcdc_fac;
        dqdc[2] =  0.69999999999999996*dcdc_fac;
        dqdc[3] =  dcdc_fac;
        dqdc[4] =  0.69999999999999996*dcdc_fac;
        dqdc[5] =  dcdc_fac + k_f*sc[10];
        dqdc[6] =  dcdc_fac;
        dqdc[7] =  dcdc_fac;
        dqdc[8] =  dcdc_fac;
        dqdc[9] =  12*dcdc_fac;
        dqdc[10] =  1.75*dcdc_fac + k_f*sc[5];
        dqdc[11] =  dcdc_fac;
        dqdc[12] =  dcdc_fac;
        dqdc[13] =  3.6000000000000001*dcdc_fac - k_r;
        for (int k=0; k<14; k++) {
            J[15*k+5] -= dqdc[k];
            J[15*k+10] -= dqdc[k];
            J[15*k+13] += dqdc[k];
        }
    }
    J[215] -= dqdT; /* dwdot[O]/dT */
    J[220] -= dqdT; /* dwdot[CO]/dT */
    J[223] += dqdT; /* dwdot[CO2]/dT */

    /*reaction 4: H + H + M <=> H2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture - sc[0] - sc[9] - sc[13] - 0.37*sc[2] - 0.37*sc[4];
    /* forward */
    phi_f = sc[1]*sc[1];
    k_f = prefactor_units[3] * fwd_A[3]
                * exp(fwd_beta[3] * tc[0] - activation_units[3] * fwd_Ea[3] * invT);
    dlnkfdT = fwd_beta[3] * invT + activation_units[3] * fwd_Ea[3] * invT2;
    /* reverse */
    phi_r = sc[0];
    Kc = refCinv * exp(-g_RT[0] + g_RT[1] + g_RT[1]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[1]) + (h_RT[0]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci =  -q_nocor - k_r;
        J[0] += dqdci;                /* dwdot[H2]/d[H2] */
        J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*2*sc[1];
        J[15] += dqdci;               /* dwdot[H2]/d[H] */
        J[16] += -2 * dqdci;          /* dwdot[H]/d[H] */
        /* d()/d[AR] */
        dqdci =  -0.37*q_nocor;
        J[30] += dqdci;               /* dwdot[H2]/d[AR] */
        J[31] += -2 * dqdci;          /* dwdot[H]/d[AR] */
        /* d()/d[HE] */
        dqdci =  -0.37*q_nocor;
        J[60] += dqdci;               /* dwdot[H2]/d[HE] */
        J[61] += -2 * dqdci;          /* dwdot[H]/d[HE] */
        /* d()/d[H2O] */
        dqdci =  -q_nocor;
        J[135] += dqdci;              /* dwdot[H2]/d[H2O] */
        J[136] += -2 * dqdci;         /* dwdot[H]/d[H2O] */
        /* d()/d[CO2] */
        dqdci =  -q_nocor;
        J[195] += dqdci;              /* dwdot[H2]/d[CO2] */
        J[196] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
    }
    else {
        dqdc[0] =  - k_r;
        dqdc[1] =  q_nocor + k_f*2*sc[1];
        dqdc[2] =  0.63*q_nocor;
        dqdc[3] =  q_nocor;
        dqdc[4] =  0.63*q_nocor;
        dqdc[5] =  q_nocor;
        dqdc[6] =  q_nocor;
        dqdc[7] =  q_nocor;
        dqdc[8] =  q_nocor;
        dqdc[10] =  q_nocor;
        dqdc[11] =  q_nocor;
        dqdc[12] =  q_nocor;
        for (int k=0; k<14; k++) {
            J[15*k+0] += dqdc[k];
            J[15*k+1] += -2 * dqdc[k];
        }
    }
    J[210] += dqdT; /* dwdot[H2]/dT */
    J[211] += -2 * dqdT; /* dwdot[H]/dT */

    /*reaction 5: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + sc[0] + 5.2999999999999998*sc[9] + 0.75*sc[10] + 2.6000000000000001*sc[13] - 0.62*sc[2] - 0.62*sc[4];
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = prefactor_units[4] * fwd_A[4]
                * exp(fwd_beta[4] * tc[0] - activation_units[4] * fwd_Ea[4] * invT);
    dlnkfdT = fwd_beta[4] * invT + activation_units[4] * fwd_Ea[4] * invT2;
    /* reverse */
    phi_r = sc[9];
    Kc = refCinv * exp(g_RT[1] + g_RT[6] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (h_RT[9]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[6] -= q; /* OH */
    wdot[9] += q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci =  q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[6] -= dqdci;                /* dwdot[OH]/d[H2] */
        J[9] += dqdci;                /* dwdot[H2O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[6];
        J[16] -= dqdci;               /* dwdot[H]/d[H] */
        J[21] -= dqdci;               /* dwdot[OH]/d[H] */
        J[24] += dqdci;               /* dwdot[H2O]/d[H] */
        /* d()/d[AR] */
        dqdci =  -0.62*q_nocor;
        J[31] -= dqdci;               /* dwdot[H]/d[AR] */
        J[36] -= dqdci;               /* dwdot[OH]/d[AR] */
        J[39] += dqdci;               /* dwdot[H2O]/d[AR] */
        /* d()/d[HE] */
        dqdci =  -0.62*q_nocor;
        J[61] -= dqdci;               /* dwdot[H]/d[HE] */
        J[66] -= dqdci;               /* dwdot[OH]/d[HE] */
        J[69] += dqdci;               /* dwdot[H2O]/d[HE] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[1];
        J[91] -= dqdci;               /* dwdot[H]/d[OH] */
        J[96] -= dqdci;               /* dwdot[OH]/d[OH] */
        J[99] += dqdci;               /* dwdot[H2O]/d[OH] */
        /* d()/d[H2O] */
        dqdci =  5.2999999999999998*q_nocor - k_r;
        J[136] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[141] -= dqdci;              /* dwdot[OH]/d[H2O] */
        J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
        /* d()/d[CO] */
        dqdci =  0.75*q_nocor;
        J[151] -= dqdci;              /* dwdot[H]/d[CO] */
        J[156] -= dqdci;              /* dwdot[OH]/d[CO] */
        J[159] += dqdci;              /* dwdot[H2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci =  2.6000000000000001*q_nocor;
        J[196] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[201] -= dqdci;              /* dwdot[OH]/d[CO2] */
        J[204] += dqdci;              /* dwdot[H2O]/d[CO2] */
    }
    else {
        dqdc[0] =  2*q_nocor;
        dqdc[1] =  q_nocor + k_f*sc[6];
        dqdc[2] =  0.38*q_nocor;
        dqdc[3] =  q_nocor;
        dqdc[4] =  0.38*q_nocor;
        dqdc[5] =  q_nocor;
        dqdc[6] =  q_nocor + k_f*sc[1];
        dqdc[7] =  q_nocor;
        dqdc[8] =  q_nocor;
        dqdc[9] =  6.2999999999999998*q_nocor - k_r;
        dqdc[10] =  1.75*q_nocor;
        dqdc[11] =  q_nocor;
        dqdc[12] =  q_nocor;
        dqdc[13] =  3.6000000000000001*q_nocor;
        for (int k=0; k<14; k++) {
            J[15*k+1] -= dqdc[k];
            J[15*k+6] -= dqdc[k];
            J[15*k+9] += dqdc[k];
        }
    }
    J[211] -= dqdT; /* dwdot[H]/dT */
    J[216] -= dqdT; /* dwdot[OH]/dT */
    J[219] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 6: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + sc[0] + 11*sc[9] + 0.75*sc[10] + 2.6000000000000001*sc[13] - 0.30000000000000004*sc[2] - 0.30000000000000004*sc[4];
    /* forward */
    phi_f = sc[1]*sc[5];
    k_f = prefactor_units[5] * fwd_A[5]
                * exp(fwd_beta[5] * tc[0] - activation_units[5] * fwd_Ea[5] * invT);
    dlnkfdT = fwd_beta[5] * invT + activation_units[5] * fwd_Ea[5] * invT2;
    /* reverse */
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[5]) + (h_RT[6]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[5] -= q; /* O */
    wdot[6] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci =  q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[5] -= dqdci;                /* dwdot[O]/d[H2] */
        J[6] += dqdci;                /* dwdot[OH]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[5];
        J[16] -= dqdci;               /* dwdot[H]/d[H] */
        J[20] -= dqdci;               /* dwdot[O]/d[H] */
        J[21] += dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[AR] */
        dqdci =  -0.30000000000000004*q_nocor;
        J[31] -= dqdci;               /* dwdot[H]/d[AR] */
        J[35] -= dqdci;               /* dwdot[O]/d[AR] */
        J[36] += dqdci;               /* dwdot[OH]/d[AR] */
        /* d()/d[HE] */
        dqdci =  -0.30000000000000004*q_nocor;
        J[61] -= dqdci;               /* dwdot[H]/d[HE] */
        J[65] -= dqdci;               /* dwdot[O]/d[HE] */
        J[66] += dqdci;               /* dwdot[OH]/d[HE] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[1];
        J[76] -= dqdci;               /* dwdot[H]/d[O] */
        J[80] -= dqdci;               /* dwdot[O]/d[O] */
        J[81] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[91] -= dqdci;               /* dwdot[H]/d[OH] */
        J[95] -= dqdci;               /* dwdot[O]/d[OH] */
        J[96] += dqdci;               /* dwdot[OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci =  11*q_nocor;
        J[136] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[140] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[141] += dqdci;              /* dwdot[OH]/d[H2O] */
        /* d()/d[CO] */
        dqdci =  0.75*q_nocor;
        J[151] -= dqdci;              /* dwdot[H]/d[CO] */
        J[155] -= dqdci;              /* dwdot[O]/d[CO] */
        J[156] += dqdci;              /* dwdot[OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci =  2.6000000000000001*q_nocor;
        J[196] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[200] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[201] += dqdci;              /* dwdot[OH]/d[CO2] */
    }
    else {
        dqdc[0] =  2*q_nocor;
        dqdc[1] =  q_nocor + k_f*sc[5];
        dqdc[2] =  0.69999999999999996*q_nocor;
        dqdc[3] =  q_nocor;
        dqdc[4] =  0.69999999999999996*q_nocor;
        dqdc[5] =  q_nocor + k_f*sc[1];
        dqdc[6] =  q_nocor - k_r;
        dqdc[7] =  q_nocor;
        dqdc[8] =  q_nocor;
        dqdc[9] =  12*q_nocor;
        dqdc[10] =  1.75*q_nocor;
        dqdc[11] =  q_nocor;
        dqdc[12] =  q_nocor;
        dqdc[13] =  3.6000000000000001*q_nocor;
        for (int k=0; k<14; k++) {
            J[15*k+1] -= dqdc[k];
            J[15*k+5] -= dqdc[k];
            J[15*k+6] += dqdc[k];
        }
    }
    J[211] -= dqdT; /* dwdot[H]/dT */
    J[215] -= dqdT; /* dwdot[O]/dT */
    J[216] += dqdT; /* dwdot[OH]/dT */

    /*reaction 7: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + 1.3999999999999999*sc[0] + 14.4*sc[9] + 0.75*sc[10] + 2.6000000000000001*sc[13] - 0.17000000000000004*sc[2] - 0.17000000000000004*sc[4];
    /* forward */
    phi_f = sc[5]*sc[5];
    k_f = prefactor_units[6] * fwd_A[6]
                * exp(fwd_beta[6] * tc[0] - activation_units[6] * fwd_Ea[6] * invT);
    dlnkfdT = fwd_beta[6] * invT + activation_units[6] * fwd_Ea[6] * invT2;
    /* reverse */
    phi_r = sc[11];
    Kc = refCinv * exp(g_RT[5] + g_RT[5] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[5]) + (h_RT[11]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= 2 * q; /* O */
    wdot[11] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci =  1.3999999999999999*q_nocor;
        J[5] += -2 * dqdci;           /* dwdot[O]/d[H2] */
        J[11] += dqdci;               /* dwdot[O2]/d[H2] */
        /* d()/d[AR] */
        dqdci =  -0.17000000000000004*q_nocor;
        J[35] += -2 * dqdci;          /* dwdot[O]/d[AR] */
        J[41] += dqdci;               /* dwdot[O2]/d[AR] */
        /* d()/d[HE] */
        dqdci =  -0.17000000000000004*q_nocor;
        J[65] += -2 * dqdci;          /* dwdot[O]/d[HE] */
        J[71] += dqdci;               /* dwdot[O2]/d[HE] */
        /* d()/d[O] */
        dqdci =  + k_f*2*sc[5];
        J[80] += -2 * dqdci;          /* dwdot[O]/d[O] */
        J[86] += dqdci;               /* dwdot[O2]/d[O] */
        /* d()/d[H2O] */
        dqdci =  14.4*q_nocor;
        J[140] += -2 * dqdci;         /* dwdot[O]/d[H2O] */
        J[146] += dqdci;              /* dwdot[O2]/d[H2O] */
        /* d()/d[CO] */
        dqdci =  0.75*q_nocor;
        J[155] += -2 * dqdci;         /* dwdot[O]/d[CO] */
        J[161] += dqdci;              /* dwdot[O2]/d[CO] */
        /* d()/d[O2] */
        dqdci =  - k_r;
        J[170] += -2 * dqdci;         /* dwdot[O]/d[O2] */
        J[176] += dqdci;              /* dwdot[O2]/d[O2] */
        /* d()/d[CO2] */
        dqdci =  2.6000000000000001*q_nocor;
        J[200] += -2 * dqdci;         /* dwdot[O]/d[CO2] */
        J[206] += dqdci;              /* dwdot[O2]/d[CO2] */
    }
    else {
        dqdc[0] =  2.3999999999999999*q_nocor;
        dqdc[1] =  q_nocor;
        dqdc[2] =  0.82999999999999996*q_nocor;
        dqdc[3] =  q_nocor;
        dqdc[4] =  0.82999999999999996*q_nocor;
        dqdc[5] =  q_nocor + k_f*2*sc[5];
        dqdc[6] =  q_nocor;
        dqdc[7] =  q_nocor;
        dqdc[8] =  q_nocor;
        dqdc[9] =  15.4*q_nocor;
        dqdc[10] =  1.75*q_nocor;
        dqdc[11] =  q_nocor - k_r;
        dqdc[12] =  q_nocor;
        dqdc[13] =  3.6000000000000001*q_nocor;
        for (int k=0; k<14; k++) {
            J[15*k+5] += -2 * dqdc[k];
            J[15*k+11] += dqdc[k];
        }
    }
    J[215] += -2 * dqdT; /* dwdot[O]/dT */
    J[221] += dqdT; /* dwdot[O2]/dT */

    /*reaction 8: HCO + M <=> CO + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + sc[0] - sc[9] + 0.75*sc[10] + 2.6000000000000001*sc[13];
    /* forward */
    phi_f = sc[7];
    k_f = prefactor_units[7] * fwd_A[7]
                * exp(fwd_beta[7] * tc[0] - activation_units[7] * fwd_Ea[7] * invT);
    dlnkfdT = fwd_beta[7] * invT + activation_units[7] * fwd_Ea[7] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[10];
    Kc = refC * exp(-g_RT[1] + g_RT[7] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7]) + (h_RT[1] + h_RT[10]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[7] -= q; /* HCO */
    wdot[10] += q; /* CO */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci =  q_nocor;
        J[1] += dqdci;                /* dwdot[H]/d[H2] */
        J[7] -= dqdci;                /* dwdot[HCO]/d[H2] */
        J[10] += dqdci;               /* dwdot[CO]/d[H2] */
        /* d()/d[H] */
        dqdci =  - k_r*sc[10];
        J[16] += dqdci;               /* dwdot[H]/d[H] */
        J[22] -= dqdci;               /* dwdot[HCO]/d[H] */
        J[25] += dqdci;               /* dwdot[CO]/d[H] */
        /* d()/d[HCO] */
        dqdci =  + k_f;
        J[106] += dqdci;              /* dwdot[H]/d[HCO] */
        J[112] -= dqdci;              /* dwdot[HCO]/d[HCO] */
        J[115] += dqdci;              /* dwdot[CO]/d[HCO] */
        /* d()/d[H2O] */
        dqdci =  -q_nocor;
        J[136] += dqdci;              /* dwdot[H]/d[H2O] */
        J[142] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        J[145] += dqdci;              /* dwdot[CO]/d[H2O] */
        /* d()/d[CO] */
        dqdci =  0.75*q_nocor - k_r*sc[1];
        J[151] += dqdci;              /* dwdot[H]/d[CO] */
        J[157] -= dqdci;              /* dwdot[HCO]/d[CO] */
        J[160] += dqdci;              /* dwdot[CO]/d[CO] */
        /* d()/d[CO2] */
        dqdci =  2.6000000000000001*q_nocor;
        J[196] += dqdci;              /* dwdot[H]/d[CO2] */
        J[202] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        J[205] += dqdci;              /* dwdot[CO]/d[CO2] */
    }
    else {
        dqdc[0] =  2*q_nocor;
        dqdc[1] =  q_nocor - k_r*sc[10];
        dqdc[2] =  q_nocor;
        dqdc[3] =  q_nocor;
        dqdc[4] =  q_nocor;
        dqdc[5] =  q_nocor;
        dqdc[6] =  q_nocor;
        dqdc[7] =  q_nocor + k_f;
        dqdc[8] =  q_nocor;
        dqdc[10] =  1.75*q_nocor - k_r*sc[1];
        dqdc[11] =  q_nocor;
        dqdc[12] =  q_nocor;
        dqdc[13] =  3.6000000000000001*q_nocor;
        for (int k=0; k<14; k++) {
            J[15*k+1] += dqdc[k];
            J[15*k+7] -= dqdc[k];
            J[15*k+10] += dqdc[k];
        }
    }
    J[211] += dqdT; /* dwdot[H]/dT */
    J[217] -= dqdT; /* dwdot[HCO]/dT */
    J[220] += dqdT; /* dwdot[CO]/dT */

    /*reaction 9: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[11];
    k_f = prefactor_units[8] * fwd_A[8]
                * exp(fwd_beta[8] * tc[0] - activation_units[8] * fwd_Ea[8] * invT);
    dlnkfdT = fwd_beta[8] * invT + activation_units[8] * fwd_Ea[8] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = exp(g_RT[1] - g_RT[5] - g_RT[6] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[11]) + (h_RT[5] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[5] += q; /* O */
    wdot[6] += q; /* OH */
    wdot[11] -= q; /* O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[11];
    J[16] -= dqdci;               /* dwdot[H]/d[H] */
    J[20] += dqdci;               /* dwdot[O]/d[H] */
    J[21] += dqdci;               /* dwdot[OH]/d[H] */
    J[26] -= dqdci;               /* dwdot[O2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[6];
    J[76] -= dqdci;               /* dwdot[H]/d[O] */
    J[80] += dqdci;               /* dwdot[O]/d[O] */
    J[81] += dqdci;               /* dwdot[OH]/d[O] */
    J[86] -= dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[5];
    J[91] -= dqdci;               /* dwdot[H]/d[OH] */
    J[95] += dqdci;               /* dwdot[O]/d[OH] */
    J[96] += dqdci;               /* dwdot[OH]/d[OH] */
    J[101] -= dqdci;              /* dwdot[O2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1];
    J[166] -= dqdci;              /* dwdot[H]/d[O2] */
    J[170] += dqdci;              /* dwdot[O]/d[O2] */
    J[171] += dqdci;              /* dwdot[OH]/d[O2] */
    J[176] -= dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[211] -= dqdT;               /* dwdot[H]/dT */
    J[215] += dqdT;               /* dwdot[O]/dT */
    J[216] += dqdT;               /* dwdot[OH]/dT */
    J[221] -= dqdT;               /* dwdot[O2]/dT */

    /*reaction 10: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = prefactor_units[9] * fwd_A[9]
                * exp(fwd_beta[9] * tc[0] - activation_units[9] * fwd_Ea[9] * invT);
    dlnkfdT = fwd_beta[9] * invT + activation_units[9] * fwd_Ea[9] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[6];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[1] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[5] -= q; /* O */
    wdot[6] += q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[5];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[5] -= dqdci;                /* dwdot[O]/d[H2] */
    J[6] += dqdci;                /* dwdot[OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[6];
    J[15] -= dqdci;               /* dwdot[H2]/d[H] */
    J[16] += dqdci;               /* dwdot[H]/d[H] */
    J[20] -= dqdci;               /* dwdot[O]/d[H] */
    J[21] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[0];
    J[75] -= dqdci;               /* dwdot[H2]/d[O] */
    J[76] += dqdci;               /* dwdot[H]/d[O] */
    J[80] -= dqdci;               /* dwdot[O]/d[O] */
    J[81] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[90] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[91] += dqdci;               /* dwdot[H]/d[OH] */
    J[95] -= dqdci;               /* dwdot[O]/d[OH] */
    J[96] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[210] -= dqdT;               /* dwdot[H2]/dT */
    J[211] += dqdT;               /* dwdot[H]/dT */
    J[215] -= dqdT;               /* dwdot[O]/dT */
    J[216] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 11: OH + H2 <=> H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[6];
    k_f = prefactor_units[10] * fwd_A[10]
                * exp(fwd_beta[10] * tc[0] - activation_units[10] * fwd_Ea[10] * invT);
    dlnkfdT = fwd_beta[10] * invT + activation_units[10] * fwd_Ea[10] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[9];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[6] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (h_RT[1] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[6] -= q; /* OH */
    wdot[9] += q; /* H2O */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[6];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[6] -= dqdci;                /* dwdot[OH]/d[H2] */
    J[9] += dqdci;                /* dwdot[H2O]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[15] -= dqdci;               /* dwdot[H2]/d[H] */
    J[16] += dqdci;               /* dwdot[H]/d[H] */
    J[21] -= dqdci;               /* dwdot[OH]/d[H] */
    J[24] += dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[0];
    J[90] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[91] += dqdci;               /* dwdot[H]/d[OH] */
    J[96] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[99] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[135] -= dqdci;              /* dwdot[H2]/d[H2O] */
    J[136] += dqdci;              /* dwdot[H]/d[H2O] */
    J[141] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[210] -= dqdT;               /* dwdot[H2]/dT */
    J[211] += dqdT;               /* dwdot[H]/dT */
    J[216] -= dqdT;               /* dwdot[OH]/dT */
    J[219] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 12: OH + OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[6];
    k_f = prefactor_units[11] * fwd_A[11]
                * exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
    dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[6]) + (h_RT[5] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* O */
    wdot[6] -= 2 * q; /* OH */
    wdot[9] += q; /* H2O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[9];
    J[80] += dqdci;               /* dwdot[O]/d[O] */
    J[81] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[84] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2*sc[6];
    J[95] += dqdci;               /* dwdot[O]/d[OH] */
    J[96] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[99] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[140] += dqdci;              /* dwdot[O]/d[H2O] */
    J[141] += -2 * dqdci;         /* dwdot[OH]/d[H2O] */
    J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[215] += dqdT;               /* dwdot[O]/dT */
    J[216] += -2 * dqdT;          /* dwdot[OH]/dT */
    J[219] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 13: H + H + H2 <=> H2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[1]*sc[1];
    k_f = prefactor_units[12] * fwd_A[12]
                * exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
    dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[0];
    Kc = refCinv * exp(g_RT[0] - g_RT[0] - g_RT[0] + g_RT[1] + g_RT[1]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + 2*h_RT[1]) + (2*h_RT[0]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[1]*sc[1] - k_r*2*sc[0];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[0]*2*sc[1];
    J[15] += dqdci;               /* dwdot[H2]/d[H] */
    J[16] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/dT */
    J[210] += dqdT;               /* dwdot[H2]/dT */
    J[211] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 14: H + H + H2O <=> H2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[1]*sc[9];
    k_f = prefactor_units[13] * fwd_A[13]
                * exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
    dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[9];
    Kc = refCinv * exp(-g_RT[0] + g_RT[1] + g_RT[1] + g_RT[9] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[1] + h_RT[9]) + (h_RT[0] + h_RT[9]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[9];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*2*sc[1]*sc[9];
    J[15] += dqdci;               /* dwdot[H2]/d[H] */
    J[16] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[1]*sc[1] - k_r*sc[0];
    J[135] += dqdci;              /* dwdot[H2]/d[H2O] */
    J[136] += -2 * dqdci;         /* dwdot[H]/d[H2O] */
    /* d()/dT */
    J[210] += dqdT;               /* dwdot[H2]/dT */
    J[211] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 15: H + H + CO2 <=> H2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[1]*sc[13];
    k_f = prefactor_units[14] * fwd_A[14]
                * exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
    dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[13];
    Kc = refCinv * exp(-g_RT[0] + g_RT[1] + g_RT[1] + g_RT[13] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[1] + h_RT[13]) + (h_RT[0] + h_RT[13]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[13];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*2*sc[1]*sc[13];
    J[15] += dqdci;               /* dwdot[H2]/d[H] */
    J[16] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[1]*sc[1] - k_r*sc[0];
    J[195] += dqdci;              /* dwdot[H2]/d[CO2] */
    J[196] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
    /* d()/dT */
    J[210] += dqdT;               /* dwdot[H2]/dT */
    J[211] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 16: H2 + O2 <=> HO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[11];
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[8];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[8] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[11]) + (h_RT[1] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[8] += q; /* HO2 */
    wdot[11] -= q; /* O2 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[11];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[8] += dqdci;                /* dwdot[HO2]/d[H2] */
    J[11] -= dqdci;               /* dwdot[O2]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[8];
    J[15] -= dqdci;               /* dwdot[H2]/d[H] */
    J[16] += dqdci;               /* dwdot[H]/d[H] */
    J[23] += dqdci;               /* dwdot[HO2]/d[H] */
    J[26] -= dqdci;               /* dwdot[O2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[1];
    J[120] -= dqdci;              /* dwdot[H2]/d[HO2] */
    J[121] += dqdci;              /* dwdot[H]/d[HO2] */
    J[128] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[131] -= dqdci;              /* dwdot[O2]/d[HO2] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[0];
    J[165] -= dqdci;              /* dwdot[H2]/d[O2] */
    J[166] += dqdci;              /* dwdot[H]/d[O2] */
    J[173] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[176] -= dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[210] -= dqdT;               /* dwdot[H2]/dT */
    J[211] += dqdT;               /* dwdot[H]/dT */
    J[218] += dqdT;               /* dwdot[HO2]/dT */
    J[221] -= dqdT;               /* dwdot[O2]/dT */

    /*reaction 17: HO2 + H <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[8];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = exp(g_RT[1] - g_RT[5] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[8]) + (h_RT[5] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[5] += q; /* O */
    wdot[8] -= q; /* HO2 */
    wdot[9] += q; /* H2O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[8];
    J[16] -= dqdci;               /* dwdot[H]/d[H] */
    J[20] += dqdci;               /* dwdot[O]/d[H] */
    J[23] -= dqdci;               /* dwdot[HO2]/d[H] */
    J[24] += dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[9];
    J[76] -= dqdci;               /* dwdot[H]/d[O] */
    J[80] += dqdci;               /* dwdot[O]/d[O] */
    J[83] -= dqdci;               /* dwdot[HO2]/d[O] */
    J[84] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[121] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[125] += dqdci;              /* dwdot[O]/d[HO2] */
    J[128] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[129] += dqdci;              /* dwdot[H2O]/d[HO2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[136] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[140] += dqdci;              /* dwdot[O]/d[H2O] */
    J[143] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[211] -= dqdT;               /* dwdot[H]/dT */
    J[215] += dqdT;               /* dwdot[O]/dT */
    J[218] -= dqdT;               /* dwdot[HO2]/dT */
    J[219] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[8];
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[6];
    Kc = exp(g_RT[1] - g_RT[6] - g_RT[6] + g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[8]) + (2*h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[6] += 2 * q; /* OH */
    wdot[8] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[8];
    J[16] -= dqdci;               /* dwdot[H]/d[H] */
    J[21] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[23] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2*sc[6];
    J[91] -= dqdci;               /* dwdot[H]/d[OH] */
    J[96] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[98] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[121] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[126] += 2 * dqdci;          /* dwdot[OH]/d[HO2] */
    J[128] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[211] -= dqdT;               /* dwdot[H]/dT */
    J[216] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[218] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 19: HO2 + O <=> OH + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[8];
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[11];
    Kc = exp(g_RT[5] - g_RT[6] + g_RT[8] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[8]) + (h_RT[6] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[6] += q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[11] += q; /* O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[8];
    J[80] -= dqdci;               /* dwdot[O]/d[O] */
    J[81] += dqdci;               /* dwdot[OH]/d[O] */
    J[83] -= dqdci;               /* dwdot[HO2]/d[O] */
    J[86] += dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[11];
    J[95] -= dqdci;               /* dwdot[O]/d[OH] */
    J[96] += dqdci;               /* dwdot[OH]/d[OH] */
    J[98] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[101] += dqdci;              /* dwdot[O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[5];
    J[125] -= dqdci;              /* dwdot[O]/d[HO2] */
    J[126] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[128] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[131] += dqdci;              /* dwdot[O2]/d[HO2] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[6];
    J[170] -= dqdci;              /* dwdot[O]/d[O2] */
    J[171] += dqdci;              /* dwdot[OH]/d[O2] */
    J[173] -= dqdci;              /* dwdot[HO2]/d[O2] */
    J[176] += dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[215] -= dqdT;               /* dwdot[O]/dT */
    J[216] += dqdT;               /* dwdot[OH]/dT */
    J[218] -= dqdT;               /* dwdot[HO2]/dT */
    J[221] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 20: HO2 + OH <=> O2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[8];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[11];
    Kc = exp(g_RT[6] + g_RT[8] - g_RT[9] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[8]) + (h_RT[9] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] -= q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[9] += q; /* H2O */
    wdot[11] += q; /* O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[96] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[98] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[99] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[101] += dqdci;              /* dwdot[O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[6];
    J[126] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[128] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[129] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[131] += dqdci;              /* dwdot[O2]/d[HO2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[11];
    J[141] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[143] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[146] += dqdci;              /* dwdot[O2]/d[H2O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[9];
    J[171] -= dqdci;              /* dwdot[OH]/d[O2] */
    J[173] -= dqdci;              /* dwdot[HO2]/d[O2] */
    J[174] += dqdci;              /* dwdot[H2O]/d[O2] */
    J[176] += dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[216] -= dqdT;               /* dwdot[OH]/dT */
    J[218] -= dqdT;               /* dwdot[HO2]/dT */
    J[219] += dqdT;               /* dwdot[H2O]/dT */
    J[221] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 21: HO2 + OH <=> O2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[8];
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[11];
    Kc = exp(g_RT[6] + g_RT[8] - g_RT[9] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[8]) + (h_RT[9] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] -= q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[9] += q; /* H2O */
    wdot[11] += q; /* O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[96] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[98] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[99] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[101] += dqdci;              /* dwdot[O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[6];
    J[126] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[128] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[129] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[131] += dqdci;              /* dwdot[O2]/d[HO2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[11];
    J[141] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[143] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[146] += dqdci;              /* dwdot[O2]/d[H2O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[9];
    J[171] -= dqdci;              /* dwdot[OH]/d[O2] */
    J[173] -= dqdci;              /* dwdot[HO2]/d[O2] */
    J[174] += dqdci;              /* dwdot[H2O]/d[O2] */
    J[176] += dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[216] -= dqdT;               /* dwdot[OH]/dT */
    J[218] -= dqdT;               /* dwdot[HO2]/dT */
    J[219] += dqdT;               /* dwdot[H2O]/dT */
    J[221] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[8];
    k_f = prefactor_units[21] * fwd_A[21]
                * exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
    dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[12];
    Kc = exp(g_RT[8] + g_RT[8] - g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[8]) + (h_RT[11] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] -= 2 * q; /* HO2 */
    wdot[11] += q; /* O2 */
    wdot[12] += q; /* H2O2 */
    /* d()/d[HO2] */
    dqdci =  + k_f*2*sc[8];
    J[128] += -2 * dqdci;         /* dwdot[HO2]/d[HO2] */
    J[131] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[132] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[12];
    J[173] += -2 * dqdci;         /* dwdot[HO2]/d[O2] */
    J[176] += dqdci;              /* dwdot[O2]/d[O2] */
    J[177] += dqdci;              /* dwdot[H2O2]/d[O2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[11];
    J[188] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[191] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[192] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[218] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[221] += dqdT;               /* dwdot[O2]/dT */
    J[222] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[8];
    k_f = prefactor_units[22] * fwd_A[22]
                * exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
    dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[12];
    Kc = exp(g_RT[8] + g_RT[8] - g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[8]) + (h_RT[11] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] -= 2 * q; /* HO2 */
    wdot[11] += q; /* O2 */
    wdot[12] += q; /* H2O2 */
    /* d()/d[HO2] */
    dqdci =  + k_f*2*sc[8];
    J[128] += -2 * dqdci;         /* dwdot[HO2]/d[HO2] */
    J[131] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[132] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[12];
    J[173] += -2 * dqdci;         /* dwdot[HO2]/d[O2] */
    J[176] += dqdci;              /* dwdot[O2]/d[O2] */
    J[177] += dqdci;              /* dwdot[H2O2]/d[O2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[11];
    J[188] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[191] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[192] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[218] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[221] += dqdT;               /* dwdot[O2]/dT */
    J[222] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[12];
    k_f = prefactor_units[23] * fwd_A[23]
                * exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
    dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[8];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[8] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[12]) + (h_RT[0] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[8] += q; /* HO2 */
    wdot[12] -= q; /* H2O2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[8];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[8] += dqdci;                /* dwdot[HO2]/d[H2] */
    J[12] -= dqdci;               /* dwdot[H2O2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[12];
    J[15] += dqdci;               /* dwdot[H2]/d[H] */
    J[16] -= dqdci;               /* dwdot[H]/d[H] */
    J[23] += dqdci;               /* dwdot[HO2]/d[H] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[0];
    J[120] += dqdci;              /* dwdot[H2]/d[HO2] */
    J[121] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[128] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[132] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[180] += dqdci;              /* dwdot[H2]/d[H2O2] */
    J[181] -= dqdci;              /* dwdot[H]/d[H2O2] */
    J[188] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[192] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[210] += dqdT;               /* dwdot[H2]/dT */
    J[211] -= dqdT;               /* dwdot[H]/dT */
    J[218] += dqdT;               /* dwdot[HO2]/dT */
    J[222] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 25: H2O2 + H <=> OH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[12];
    k_f = prefactor_units[24] * fwd_A[24]
                * exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
    dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[9];
    Kc = exp(g_RT[1] - g_RT[6] - g_RT[9] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[12]) + (h_RT[6] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[6] += q; /* OH */
    wdot[9] += q; /* H2O */
    wdot[12] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[12];
    J[16] -= dqdci;               /* dwdot[H]/d[H] */
    J[21] += dqdci;               /* dwdot[OH]/d[H] */
    J[24] += dqdci;               /* dwdot[H2O]/d[H] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[9];
    J[91] -= dqdci;               /* dwdot[H]/d[OH] */
    J[96] += dqdci;               /* dwdot[OH]/d[OH] */
    J[99] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[102] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[136] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[141] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[147] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[181] -= dqdci;              /* dwdot[H]/d[H2O2] */
    J[186] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[189] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[192] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[211] -= dqdT;               /* dwdot[H]/dT */
    J[216] += dqdT;               /* dwdot[OH]/dT */
    J[219] += dqdT;               /* dwdot[H2O]/dT */
    J[222] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[12];
    k_f = prefactor_units[25] * fwd_A[25]
                * exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
    dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[8];
    Kc = exp(g_RT[5] - g_RT[6] - g_RT[8] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[12]) + (h_RT[6] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[6] += q; /* OH */
    wdot[8] += q; /* HO2 */
    wdot[12] -= q; /* H2O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[12];
    J[80] -= dqdci;               /* dwdot[O]/d[O] */
    J[81] += dqdci;               /* dwdot[OH]/d[O] */
    J[83] += dqdci;               /* dwdot[HO2]/d[O] */
    J[87] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[8];
    J[95] -= dqdci;               /* dwdot[O]/d[OH] */
    J[96] += dqdci;               /* dwdot[OH]/d[OH] */
    J[98] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[102] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[6];
    J[125] -= dqdci;              /* dwdot[O]/d[HO2] */
    J[126] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[128] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[132] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[5];
    J[185] -= dqdci;              /* dwdot[O]/d[H2O2] */
    J[186] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[188] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[192] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[215] -= dqdT;               /* dwdot[O]/dT */
    J[216] += dqdT;               /* dwdot[OH]/dT */
    J[218] += dqdT;               /* dwdot[HO2]/dT */
    J[222] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[12];
    k_f = prefactor_units[26] * fwd_A[26]
                * exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
    dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[9];
    Kc = exp(g_RT[6] - g_RT[8] - g_RT[9] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[8] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] -= q; /* OH */
    wdot[8] += q; /* HO2 */
    wdot[9] += q; /* H2O */
    wdot[12] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[12];
    J[96] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[98] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[99] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[102] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[9];
    J[126] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[128] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[129] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[132] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[8];
    J[141] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[143] += dqdci;              /* dwdot[HO2]/d[H2O] */
    J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[147] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[6];
    J[186] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[188] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[189] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[192] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[216] -= dqdT;               /* dwdot[OH]/dT */
    J[218] += dqdT;               /* dwdot[HO2]/dT */
    J[219] += dqdT;               /* dwdot[H2O]/dT */
    J[222] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[12];
    k_f = prefactor_units[27] * fwd_A[27]
                * exp(fwd_beta[27] * tc[0] - activation_units[27] * fwd_Ea[27] * invT);
    dlnkfdT = fwd_beta[27] * invT + activation_units[27] * fwd_Ea[27] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[9];
    Kc = exp(g_RT[6] - g_RT[8] - g_RT[9] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[8] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] -= q; /* OH */
    wdot[8] += q; /* HO2 */
    wdot[9] += q; /* H2O */
    wdot[12] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[12];
    J[96] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[98] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[99] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[102] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[9];
    J[126] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[128] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[129] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[132] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[8];
    J[141] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[143] += dqdci;              /* dwdot[HO2]/d[H2O] */
    J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[147] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[6];
    J[186] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[188] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[189] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[192] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[216] -= dqdT;               /* dwdot[OH]/dT */
    J[218] += dqdT;               /* dwdot[HO2]/dT */
    J[219] += dqdT;               /* dwdot[H2O]/dT */
    J[222] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 29: CO + OH <=> CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[10];
    k_f = prefactor_units[28] * fwd_A[28]
                * exp(fwd_beta[28] * tc[0] - activation_units[28] * fwd_Ea[28] * invT);
    dlnkfdT = fwd_beta[28] * invT + activation_units[28] * fwd_Ea[28] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[13];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[10] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[10]) + (h_RT[1] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[6] -= q; /* OH */
    wdot[10] -= q; /* CO */
    wdot[13] += q; /* CO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[16] += dqdci;               /* dwdot[H]/d[H] */
    J[21] -= dqdci;               /* dwdot[OH]/d[H] */
    J[25] -= dqdci;               /* dwdot[CO]/d[H] */
    J[28] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[91] += dqdci;               /* dwdot[H]/d[OH] */
    J[96] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[100] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[103] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[6];
    J[151] += dqdci;              /* dwdot[H]/d[CO] */
    J[156] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[160] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[163] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[196] += dqdci;              /* dwdot[H]/d[CO2] */
    J[201] -= dqdci;              /* dwdot[OH]/d[CO2] */
    J[205] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[208] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[211] += dqdT;               /* dwdot[H]/dT */
    J[216] -= dqdT;               /* dwdot[OH]/dT */
    J[220] -= dqdT;               /* dwdot[CO]/dT */
    J[223] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 30: CO + OH <=> CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[10];
    k_f = prefactor_units[29] * fwd_A[29]
                * exp(fwd_beta[29] * tc[0] - activation_units[29] * fwd_Ea[29] * invT);
    dlnkfdT = fwd_beta[29] * invT + activation_units[29] * fwd_Ea[29] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[13];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[10] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[10]) + (h_RT[1] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[6] -= q; /* OH */
    wdot[10] -= q; /* CO */
    wdot[13] += q; /* CO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[16] += dqdci;               /* dwdot[H]/d[H] */
    J[21] -= dqdci;               /* dwdot[OH]/d[H] */
    J[25] -= dqdci;               /* dwdot[CO]/d[H] */
    J[28] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[91] += dqdci;               /* dwdot[H]/d[OH] */
    J[96] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[100] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[103] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[6];
    J[151] += dqdci;              /* dwdot[H]/d[CO] */
    J[156] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[160] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[163] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[196] += dqdci;              /* dwdot[H]/d[CO2] */
    J[201] -= dqdci;              /* dwdot[OH]/d[CO2] */
    J[205] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[208] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[211] += dqdT;               /* dwdot[H]/dT */
    J[216] -= dqdT;               /* dwdot[OH]/dT */
    J[220] -= dqdT;               /* dwdot[CO]/dT */
    J[223] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 31: CO + O2 <=> CO2 + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[11];
    k_f = prefactor_units[30] * fwd_A[30]
                * exp(fwd_beta[30] * tc[0] - activation_units[30] * fwd_Ea[30] * invT);
    dlnkfdT = fwd_beta[30] * invT + activation_units[30] * fwd_Ea[30] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[13];
    Kc = exp(-g_RT[5] + g_RT[10] + g_RT[11] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[11]) + (h_RT[5] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* O */
    wdot[10] -= q; /* CO */
    wdot[11] -= q; /* O2 */
    wdot[13] += q; /* CO2 */
    /* d()/d[O] */
    dqdci =  - k_r*sc[13];
    J[80] += dqdci;               /* dwdot[O]/d[O] */
    J[85] -= dqdci;               /* dwdot[CO]/d[O] */
    J[86] -= dqdci;               /* dwdot[O2]/d[O] */
    J[88] += dqdci;               /* dwdot[CO2]/d[O] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[11];
    J[155] += dqdci;              /* dwdot[O]/d[CO] */
    J[160] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[161] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[163] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[10];
    J[170] += dqdci;              /* dwdot[O]/d[O2] */
    J[175] -= dqdci;              /* dwdot[CO]/d[O2] */
    J[176] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[178] += dqdci;              /* dwdot[CO2]/d[O2] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[5];
    J[200] += dqdci;              /* dwdot[O]/d[CO2] */
    J[205] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[206] -= dqdci;              /* dwdot[O2]/d[CO2] */
    J[208] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[215] += dqdT;               /* dwdot[O]/dT */
    J[220] -= dqdT;               /* dwdot[CO]/dT */
    J[221] -= dqdT;               /* dwdot[O2]/dT */
    J[223] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 32: CO + HO2 <=> CO2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[10];
    k_f = prefactor_units[31] * fwd_A[31]
                * exp(fwd_beta[31] * tc[0] - activation_units[31] * fwd_Ea[31] * invT);
    dlnkfdT = fwd_beta[31] * invT + activation_units[31] * fwd_Ea[31] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[13];
    Kc = exp(-g_RT[6] + g_RT[8] + g_RT[10] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[10]) + (h_RT[6] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] += q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[10] -= q; /* CO */
    wdot[13] += q; /* CO2 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[13];
    J[96] += dqdci;               /* dwdot[OH]/d[OH] */
    J[98] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[100] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[103] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[10];
    J[126] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[128] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[130] -= dqdci;              /* dwdot[CO]/d[HO2] */
    J[133] += dqdci;              /* dwdot[CO2]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[8];
    J[156] += dqdci;              /* dwdot[OH]/d[CO] */
    J[158] -= dqdci;              /* dwdot[HO2]/d[CO] */
    J[160] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[163] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[6];
    J[201] += dqdci;              /* dwdot[OH]/d[CO2] */
    J[203] -= dqdci;              /* dwdot[HO2]/d[CO2] */
    J[205] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[208] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[216] += dqdT;               /* dwdot[OH]/dT */
    J[218] -= dqdT;               /* dwdot[HO2]/dT */
    J[220] -= dqdT;               /* dwdot[CO]/dT */
    J[223] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 33: HCO + H <=> CO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = prefactor_units[32] * fwd_A[32]
                * exp(fwd_beta[32] * tc[0] - activation_units[32] * fwd_Ea[32] * invT);
    dlnkfdT = fwd_beta[32] * invT + activation_units[32] * fwd_Ea[32] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[10];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[7] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[7]) + (h_RT[0] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[7] -= q; /* HCO */
    wdot[10] += q; /* CO */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[10];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[7] -= dqdci;                /* dwdot[HCO]/d[H2] */
    J[10] += dqdci;               /* dwdot[CO]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[15] += dqdci;               /* dwdot[H2]/d[H] */
    J[16] -= dqdci;               /* dwdot[H]/d[H] */
    J[22] -= dqdci;               /* dwdot[HCO]/d[H] */
    J[25] += dqdci;               /* dwdot[CO]/d[H] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[1];
    J[105] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[106] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[112] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[115] += dqdci;              /* dwdot[CO]/d[HCO] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[0];
    J[150] += dqdci;              /* dwdot[H2]/d[CO] */
    J[151] -= dqdci;              /* dwdot[H]/d[CO] */
    J[157] -= dqdci;              /* dwdot[HCO]/d[CO] */
    J[160] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[210] += dqdT;               /* dwdot[H2]/dT */
    J[211] -= dqdT;               /* dwdot[H]/dT */
    J[217] -= dqdT;               /* dwdot[HCO]/dT */
    J[220] += dqdT;               /* dwdot[CO]/dT */

    /*reaction 34: HCO + O <=> CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = prefactor_units[33] * fwd_A[33]
                * exp(fwd_beta[33] * tc[0] - activation_units[33] * fwd_Ea[33] * invT);
    dlnkfdT = fwd_beta[33] * invT + activation_units[33] * fwd_Ea[33] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[10];
    Kc = exp(g_RT[5] - g_RT[6] + g_RT[7] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[6] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[6] += q; /* OH */
    wdot[7] -= q; /* HCO */
    wdot[10] += q; /* CO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[80] -= dqdci;               /* dwdot[O]/d[O] */
    J[81] += dqdci;               /* dwdot[OH]/d[O] */
    J[82] -= dqdci;               /* dwdot[HCO]/d[O] */
    J[85] += dqdci;               /* dwdot[CO]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[10];
    J[95] -= dqdci;               /* dwdot[O]/d[OH] */
    J[96] += dqdci;               /* dwdot[OH]/d[OH] */
    J[97] -= dqdci;               /* dwdot[HCO]/d[OH] */
    J[100] += dqdci;              /* dwdot[CO]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[5];
    J[110] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[111] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[112] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[115] += dqdci;              /* dwdot[CO]/d[HCO] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[6];
    J[155] -= dqdci;              /* dwdot[O]/d[CO] */
    J[156] += dqdci;              /* dwdot[OH]/d[CO] */
    J[157] -= dqdci;              /* dwdot[HCO]/d[CO] */
    J[160] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[215] -= dqdT;               /* dwdot[O]/dT */
    J[216] += dqdT;               /* dwdot[OH]/dT */
    J[217] -= dqdT;               /* dwdot[HCO]/dT */
    J[220] += dqdT;               /* dwdot[CO]/dT */

    /*reaction 35: HCO + O <=> CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = prefactor_units[34] * fwd_A[34]
                * exp(fwd_beta[34] * tc[0] - activation_units[34] * fwd_Ea[34] * invT);
    dlnkfdT = fwd_beta[34] * invT + activation_units[34] * fwd_Ea[34] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[13];
    Kc = exp(-g_RT[1] + g_RT[5] + g_RT[7] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[1] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[5] -= q; /* O */
    wdot[7] -= q; /* HCO */
    wdot[13] += q; /* CO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[16] += dqdci;               /* dwdot[H]/d[H] */
    J[20] -= dqdci;               /* dwdot[O]/d[H] */
    J[22] -= dqdci;               /* dwdot[HCO]/d[H] */
    J[28] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[76] += dqdci;               /* dwdot[H]/d[O] */
    J[80] -= dqdci;               /* dwdot[O]/d[O] */
    J[82] -= dqdci;               /* dwdot[HCO]/d[O] */
    J[88] += dqdci;               /* dwdot[CO2]/d[O] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[5];
    J[106] += dqdci;              /* dwdot[H]/d[HCO] */
    J[110] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[112] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[118] += dqdci;              /* dwdot[CO2]/d[HCO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[196] += dqdci;              /* dwdot[H]/d[CO2] */
    J[200] -= dqdci;              /* dwdot[O]/d[CO2] */
    J[202] -= dqdci;              /* dwdot[HCO]/d[CO2] */
    J[208] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[211] += dqdT;               /* dwdot[H]/dT */
    J[215] -= dqdT;               /* dwdot[O]/dT */
    J[217] -= dqdT;               /* dwdot[HCO]/dT */
    J[223] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 36: HCO + OH <=> CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[7];
    k_f = prefactor_units[35] * fwd_A[35]
                * exp(fwd_beta[35] * tc[0] - activation_units[35] * fwd_Ea[35] * invT);
    dlnkfdT = fwd_beta[35] * invT + activation_units[35] * fwd_Ea[35] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[10];
    Kc = exp(g_RT[6] + g_RT[7] - g_RT[9] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[7]) + (h_RT[9] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] -= q; /* OH */
    wdot[7] -= q; /* HCO */
    wdot[9] += q; /* H2O */
    wdot[10] += q; /* CO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[96] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[97] -= dqdci;               /* dwdot[HCO]/d[OH] */
    J[99] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[100] += dqdci;              /* dwdot[CO]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[6];
    J[111] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[112] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[114] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[115] += dqdci;              /* dwdot[CO]/d[HCO] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[10];
    J[141] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[142] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[145] += dqdci;              /* dwdot[CO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[9];
    J[156] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[157] -= dqdci;              /* dwdot[HCO]/d[CO] */
    J[159] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[160] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[216] -= dqdT;               /* dwdot[OH]/dT */
    J[217] -= dqdT;               /* dwdot[HCO]/dT */
    J[219] += dqdT;               /* dwdot[H2O]/dT */
    J[220] += dqdT;               /* dwdot[CO]/dT */

    /*reaction 37: HCO + H2O <=> CO + H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[9];
    k_f = prefactor_units[36] * fwd_A[36]
                * exp(fwd_beta[36] * tc[0] - activation_units[36] * fwd_Ea[36] * invT);
    dlnkfdT = fwd_beta[36] * invT + activation_units[36] * fwd_Ea[36] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[9]*sc[10];
    Kc = refC * exp(-g_RT[1] + g_RT[7] + g_RT[9] - g_RT[9] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[9]) + (h_RT[1] + h_RT[9] + h_RT[10]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[7] -= q; /* HCO */
    wdot[10] += q; /* CO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9]*sc[10];
    J[16] += dqdci;               /* dwdot[H]/d[H] */
    J[22] -= dqdci;               /* dwdot[HCO]/d[H] */
    J[25] += dqdci;               /* dwdot[CO]/d[H] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[9];
    J[106] += dqdci;              /* dwdot[H]/d[HCO] */
    J[112] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[115] += dqdci;              /* dwdot[CO]/d[HCO] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[7] - k_r*sc[1]*sc[10];
    J[136] += dqdci;              /* dwdot[H]/d[H2O] */
    J[142] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    J[145] += dqdci;              /* dwdot[CO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*sc[9];
    J[151] += dqdci;              /* dwdot[H]/d[CO] */
    J[157] -= dqdci;              /* dwdot[HCO]/d[CO] */
    J[160] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[211] += dqdT;               /* dwdot[H]/dT */
    J[217] -= dqdT;               /* dwdot[HCO]/dT */
    J[220] += dqdT;               /* dwdot[CO]/dT */

    /*reaction 38: HCO + O2 <=> CO + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[11];
    k_f = prefactor_units[37] * fwd_A[37]
                * exp(fwd_beta[37] * tc[0] - activation_units[37] * fwd_Ea[37] * invT);
    dlnkfdT = fwd_beta[37] * invT + activation_units[37] * fwd_Ea[37] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[10];
    Kc = exp(g_RT[7] - g_RT[8] - g_RT[10] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[11]) + (h_RT[8] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] -= q; /* HCO */
    wdot[8] += q; /* HO2 */
    wdot[10] += q; /* CO */
    wdot[11] -= q; /* O2 */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[11];
    J[112] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[113] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[115] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[116] -= dqdci;              /* dwdot[O2]/d[HCO] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[10];
    J[127] -= dqdci;              /* dwdot[HCO]/d[HO2] */
    J[128] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[130] += dqdci;              /* dwdot[CO]/d[HO2] */
    J[131] -= dqdci;              /* dwdot[O2]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[8];
    J[157] -= dqdci;              /* dwdot[HCO]/d[CO] */
    J[158] += dqdci;              /* dwdot[HO2]/d[CO] */
    J[160] += dqdci;              /* dwdot[CO]/d[CO] */
    J[161] -= dqdci;              /* dwdot[O2]/d[CO] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[7];
    J[172] -= dqdci;              /* dwdot[HCO]/d[O2] */
    J[173] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[175] += dqdci;              /* dwdot[CO]/d[O2] */
    J[176] -= dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[217] -= dqdT;               /* dwdot[HCO]/dT */
    J[218] += dqdT;               /* dwdot[HO2]/dT */
    J[220] += dqdT;               /* dwdot[CO]/dT */
    J[221] -= dqdT;               /* dwdot[O2]/dT */

    double c_R[14], dcRdT[14], e_RT[14];
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
    for (int k = 0; k < 14; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[210+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 14; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 14; ++m) {
            dehmixdc += eh_RT[m]*J[k*15+m];
        }
        J[k*15+14] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[224] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}


/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void dcvpRdT(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +7.98052075e-03
            -3.89563020e-05 * tc[1]
            +6.04716282e-08 * tc[2]
            -2.95044704e-11 * tc[3];
        /*species 1: H */
        species[1] =
            +7.05332819e-13
            -3.99183928e-15 * tc[1]
            +6.90244896e-18 * tc[2]
            -3.71092933e-21 * tc[3];
        /*species 2: AR */
        species[2] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 3: N2 */
        species[3] =
            +1.40824040e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77941600e-12 * tc[3];
        /*species 4: HE */
        species[4] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 5: O */
        species[5] =
            -3.27931884e-03
            +1.32861279e-05 * tc[1]
            -1.83841987e-08 * tc[2]
            +8.45063884e-12 * tc[3];
        /*species 6: OH */
        species[6] =
            -3.22544939e-03
            +1.30552938e-05 * tc[1]
            -1.73956093e-08 * tc[2]
            +8.24949516e-12 * tc[3];
        /*species 7: HCO */
        species[7] =
            -3.24392532e-03
            +2.75598892e-05 * tc[1]
            -3.99432279e-08 * tc[2]
            +1.73507546e-11 * tc[3];
        /*species 8: HO2 */
        species[8] =
            -4.74912051e-03
            +4.23165782e-05 * tc[1]
            -7.28291682e-08 * tc[2]
            +3.71690050e-11 * tc[3];
        /*species 9: H2O */
        species[9] =
            -2.03643410e-03
            +1.30408042e-05 * tc[1]
            -1.64639119e-08 * tc[2]
            +7.08791268e-12 * tc[3];
        /*species 10: CO */
        species[10] =
            -6.10353680e-04
            +2.03362866e-06 * tc[1]
            +2.72101765e-09 * tc[2]
            -3.61769800e-12 * tc[3];
        /*species 11: O2 */
        species[11] =
            -2.99673416e-03
            +1.96946040e-05 * tc[1]
            -2.90438853e-08 * tc[2]
            +1.29749135e-11 * tc[3];
        /*species 12: H2O2 */
        species[12] =
            -5.42822417e-04
            +3.34671402e-05 * tc[1]
            -6.47312439e-08 * tc[2]
            +3.44981745e-11 * tc[3];
        /*species 13: CO2 */
        species[13] =
            +8.98459677e-03
            -1.42471254e-05 * tc[1]
            +7.37757066e-09 * tc[2]
            -5.74798192e-13 * tc[3];
    } else {
        /*species 0: H2 */
        species[0] =
            -4.94024731e-05
            +9.98913556e-07 * tc[1]
            -5.38699182e-10 * tc[2]
            +8.01021504e-14 * tc[3];
        /*species 1: H */
        species[1] =
            -2.30842973e-11
            +3.23123896e-14 * tc[1]
            -1.42054571e-17 * tc[2]
            +1.99278943e-21 * tc[3];
        /*species 2: AR */
        species[2] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 3: N2 */
        species[3] =
            +1.48797680e-03
            -1.13695200e-06 * tc[1]
            +3.02911140e-10 * tc[2]
            -2.70134040e-14 * tc[3];
        /*species 4: HE */
        species[4] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 5: O */
        species[5] =
            -8.59741137e-05
            +8.38969178e-08 * tc[1]
            -3.00533397e-11 * tc[2]
            +4.91334764e-15 * tc[3];
        /*species 6: OH */
        species[6] =
            +1.05650448e-03
            -5.18165516e-07 * tc[1]
            +9.15656022e-11 * tc[2]
            -5.32783504e-15 * tc[3];
        /*species 7: HCO */
        species[7] =
            +4.95695526e-03
            -4.96891226e-06 * tc[1]
            +1.76748533e-09 * tc[2]
            -2.13403484e-13 * tc[3];
        /*species 8: HO2 */
        species[8] =
            +2.23982013e-03
            -1.26731630e-06 * tc[1]
            +3.42739110e-10 * tc[2]
            -4.31634140e-14 * tc[3];
        /*species 9: H2O */
        species[9] =
            +2.17691804e-03
            -3.28145036e-07 * tc[1]
            -2.91125961e-10 * tc[2]
            +6.72803968e-14 * tc[3];
        /*species 10: CO */
        species[10] =
            +2.06252743e-03
            -1.99765154e-06 * tc[1]
            +6.90159024e-10 * tc[2]
            -8.14590864e-14 * tc[3];
        /*species 11: O2 */
        species[11] =
            +1.48308754e-03
            -1.51593334e-06 * tc[1]
            +6.28411665e-10 * tc[2]
            -8.66871176e-14 * tc[3];
        /*species 12: H2O2 */
        species[12] =
            +4.90831694e-03
            -3.80278450e-06 * tc[1]
            +1.11355796e-09 * tc[2]
            -1.15163322e-13 * tc[3];
        /*species 13: CO2 */
        species[13] =
            +4.41437026e-03
            -4.42962808e-06 * tc[1]
            +1.57047056e-09 * tc[2]
            -1.88833666e-13 * tc[3];
    }
    return;
}


/*compute the progress rate for each reaction */
void progressRate(double * restrict qdot, double * restrict sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    double q_f[38], q_r[38];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 38; ++i) {
        qdot[i] = q_f[i] - q_r[i];
    }

    return;
}


/*compute the progress rate for each reaction */
void progressRateFR(double * restrict q_f, double * restrict q_r, double * restrict sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    comp_qfqr(q_f, q_r, sc, tc, invT);

    return;
}


/*compute the equilibrium constants for each reaction */
void equilibriumConstants(double * restrict kc, double * restrict g_RT, double T)
{
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    kc[0] = 1.0 / (refC) * exp((g_RT[1] + g_RT[11]) - (g_RT[8]));

    /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
    kc[1] = 1.0 / (refC) * exp((g_RT[6] + g_RT[6]) - (g_RT[12]));

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    kc[2] = 1.0 / (refC) * exp((g_RT[10] + g_RT[5]) - (g_RT[13]));

    /*reaction 4: H + H + M <=> H2 + M */
    kc[3] = 1.0 / (refC) * exp((g_RT[1] + g_RT[1]) - (g_RT[0]));

    /*reaction 5: H + OH + M <=> H2O + M */
    kc[4] = 1.0 / (refC) * exp((g_RT[1] + g_RT[6]) - (g_RT[9]));

    /*reaction 6: O + H + M <=> OH + M */
    kc[5] = 1.0 / (refC) * exp((g_RT[5] + g_RT[1]) - (g_RT[6]));

    /*reaction 7: O + O + M <=> O2 + M */
    kc[6] = 1.0 / (refC) * exp((g_RT[5] + g_RT[5]) - (g_RT[11]));

    /*reaction 8: HCO + M <=> CO + H + M */
    kc[7] = refC * exp((g_RT[7]) - (g_RT[10] + g_RT[1]));

    /*reaction 9: H + O2 <=> O + OH */
    kc[8] = exp((g_RT[1] + g_RT[11]) - (g_RT[5] + g_RT[6]));

    /*reaction 10: O + H2 <=> H + OH */
    kc[9] = exp((g_RT[5] + g_RT[0]) - (g_RT[1] + g_RT[6]));

    /*reaction 11: OH + H2 <=> H + H2O */
    kc[10] = exp((g_RT[6] + g_RT[0]) - (g_RT[1] + g_RT[9]));

    /*reaction 12: OH + OH <=> O + H2O */
    kc[11] = exp((g_RT[6] + g_RT[6]) - (g_RT[5] + g_RT[9]));

    /*reaction 13: H + H + H2 <=> H2 + H2 */
    kc[12] = 1.0 / (refC) * exp((g_RT[1] + g_RT[1] + g_RT[0]) - (g_RT[0] + g_RT[0]));

    /*reaction 14: H + H + H2O <=> H2 + H2O */
    kc[13] = 1.0 / (refC) * exp((g_RT[1] + g_RT[1] + g_RT[9]) - (g_RT[0] + g_RT[9]));

    /*reaction 15: H + H + CO2 <=> H2 + CO2 */
    kc[14] = 1.0 / (refC) * exp((g_RT[1] + g_RT[1] + g_RT[13]) - (g_RT[0] + g_RT[13]));

    /*reaction 16: H2 + O2 <=> HO2 + H */
    kc[15] = exp((g_RT[0] + g_RT[11]) - (g_RT[8] + g_RT[1]));

    /*reaction 17: HO2 + H <=> O + H2O */
    kc[16] = exp((g_RT[8] + g_RT[1]) - (g_RT[5] + g_RT[9]));

    /*reaction 18: HO2 + H <=> OH + OH */
    kc[17] = exp((g_RT[8] + g_RT[1]) - (g_RT[6] + g_RT[6]));

    /*reaction 19: HO2 + O <=> OH + O2 */
    kc[18] = exp((g_RT[8] + g_RT[5]) - (g_RT[6] + g_RT[11]));

    /*reaction 20: HO2 + OH <=> O2 + H2O */
    kc[19] = exp((g_RT[8] + g_RT[6]) - (g_RT[11] + g_RT[9]));

    /*reaction 21: HO2 + OH <=> O2 + H2O */
    kc[20] = exp((g_RT[8] + g_RT[6]) - (g_RT[11] + g_RT[9]));

    /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
    kc[21] = exp((g_RT[8] + g_RT[8]) - (g_RT[11] + g_RT[12]));

    /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
    kc[22] = exp((g_RT[8] + g_RT[8]) - (g_RT[11] + g_RT[12]));

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    kc[23] = exp((g_RT[12] + g_RT[1]) - (g_RT[8] + g_RT[0]));

    /*reaction 25: H2O2 + H <=> OH + H2O */
    kc[24] = exp((g_RT[12] + g_RT[1]) - (g_RT[6] + g_RT[9]));

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    kc[25] = exp((g_RT[12] + g_RT[5]) - (g_RT[6] + g_RT[8]));

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    kc[26] = exp((g_RT[12] + g_RT[6]) - (g_RT[8] + g_RT[9]));

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    kc[27] = exp((g_RT[12] + g_RT[6]) - (g_RT[8] + g_RT[9]));

    /*reaction 29: CO + OH <=> CO2 + H */
    kc[28] = exp((g_RT[10] + g_RT[6]) - (g_RT[13] + g_RT[1]));

    /*reaction 30: CO + OH <=> CO2 + H */
    kc[29] = exp((g_RT[10] + g_RT[6]) - (g_RT[13] + g_RT[1]));

    /*reaction 31: CO + O2 <=> CO2 + O */
    kc[30] = exp((g_RT[10] + g_RT[11]) - (g_RT[13] + g_RT[5]));

    /*reaction 32: CO + HO2 <=> CO2 + OH */
    kc[31] = exp((g_RT[10] + g_RT[8]) - (g_RT[13] + g_RT[6]));

    /*reaction 33: HCO + H <=> CO + H2 */
    kc[32] = exp((g_RT[7] + g_RT[1]) - (g_RT[10] + g_RT[0]));

    /*reaction 34: HCO + O <=> CO + OH */
    kc[33] = exp((g_RT[7] + g_RT[5]) - (g_RT[10] + g_RT[6]));

    /*reaction 35: HCO + O <=> CO2 + H */
    kc[34] = exp((g_RT[7] + g_RT[5]) - (g_RT[13] + g_RT[1]));

    /*reaction 36: HCO + OH <=> CO + H2O */
    kc[35] = exp((g_RT[7] + g_RT[6]) - (g_RT[10] + g_RT[9]));

    /*reaction 37: HCO + H2O <=> CO + H + H2O */
    kc[36] = refC * exp((g_RT[7] + g_RT[9]) - (g_RT[10] + g_RT[1] + g_RT[9]));

    /*reaction 38: HCO + O2 <=> CO + HO2 */
    kc[37] = exp((g_RT[7] + g_RT[11]) - (g_RT[10] + g_RT[8]));

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void gibbs(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

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
        /*species 2: AR */
        species[2] =
            -7.453750000000000e+02 * invT
            -1.866000000000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 3: N2 */
        species[3] =
            -1.020899900000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041202000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427000000000e-13 * tc[4];
        /*species 4: HE */
        species[4] =
            -7.453750000000000e+02 * invT
            +1.571276026000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 5: O */
        species[5] =
            +2.912225920000000e+04 * invT
            +1.116333640000000e+00
            -3.168267100000000e+00 * tc[0]
            +1.639659420000000e-03 * tc[1]
            -1.107177326666667e-06 * tc[2]
            +5.106721866666666e-10 * tc[3]
            -1.056329855000000e-13 * tc[4];
        /*species 6: OH */
        species[6] =
            +3.381538120000000e+03 * invT
            +4.815738570000000e+00
            -4.125305610000000e+00 * tc[0]
            +1.612724695000000e-03 * tc[1]
            -1.087941151666667e-06 * tc[2]
            +4.832113691666666e-10 * tc[3]
            -1.031186895000000e-13 * tc[4];
        /*species 7: HCO */
        species[7] =
            +3.839564960000000e+03 * invT
            +8.268134100000002e-01
            -4.221185840000000e+00 * tc[0]
            +1.621962660000000e-03 * tc[1]
            -2.296657433333333e-06 * tc[2]
            +1.109534108333333e-09 * tc[3]
            -2.168844325000000e-13 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +2.948080400000000e+02 * invT
            +5.851355599999999e-01
            -4.301798010000000e+00 * tc[0]
            +2.374560255000000e-03 * tc[1]
            -3.526381516666666e-06 * tc[2]
            +2.023032450000000e-09 * tc[3]
            -4.646125620000001e-13 * tc[4];
        /*species 9: H2O */
        species[9] =
            -3.029372670000000e+04 * invT
            +5.047672768000000e+00
            -4.198640560000000e+00 * tc[0]
            +1.018217050000000e-03 * tc[1]
            -1.086733685000000e-06 * tc[2]
            +4.573308850000000e-10 * tc[3]
            -8.859890850000000e-14 * tc[4];
        /*species 10: CO */
        species[10] =
            -1.434408600000000e+04 * invT
            +7.112418999999992e-02
            -3.579533470000000e+00 * tc[0]
            +3.051768400000000e-04 * tc[1]
            -1.694690550000000e-07 * tc[2]
            -7.558382366666667e-11 * tc[3]
            +4.522122495000000e-14 * tc[4];
        /*species 11: O2 */
        species[11] =
            -1.063943560000000e+03 * invT
            +1.247806300000001e-01
            -3.782456360000000e+00 * tc[0]
            +1.498367080000000e-03 * tc[1]
            -1.641217001666667e-06 * tc[2]
            +8.067745908333334e-10 * tc[3]
            -1.621864185000000e-13 * tc[4];
        /*species 12: H2O2 */
        species[12] =
            -1.770258210000000e+04 * invT
            +8.410619499999998e-01
            -4.276112690000000e+00 * tc[0]
            +2.714112085000000e-04 * tc[1]
            -2.788928350000000e-06 * tc[2]
            +1.798090108333333e-09 * tc[3]
            -4.312271815000000e-13 * tc[4];
        /*species 13: CO2 */
        species[13] =
            -4.837196970000000e+04 * invT
            -7.544278700000000e+00
            -2.356773520000000e+00 * tc[0]
            -4.492298385000000e-03 * tc[1]
            +1.187260448333333e-06 * tc[2]
            -2.049325183333333e-10 * tc[3]
            +7.184977399999999e-15 * tc[4];
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
        /*species 2: AR */
        species[2] =
            -7.453750000000000e+02 * invT
            -1.866000000000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 3: N2 */
        species[3] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439884000000000e-04 * tc[1]
            +9.474600000000001e-08 * tc[2]
            -8.414198333333333e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
        /*species 4: HE */
        species[4] =
            -7.453750000000000e+02 * invT
            +1.571276026000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 5: O */
        species[5] =
            +2.921757910000000e+04 * invT
            -2.214917859999999e+00
            -2.569420780000000e+00 * tc[0]
            +4.298705685000000e-05 * tc[1]
            -6.991409816666667e-09 * tc[2]
            +8.348149916666666e-13 * tc[3]
            -6.141684549999999e-17 * tc[4];
        /*species 6: OH */
        species[6] =
            +3.718857740000000e+03 * invT
            -2.836911870000000e+00
            -2.864728860000000e+00 * tc[0]
            -5.282522400000000e-04 * tc[1]
            +4.318045966666667e-08 * tc[2]
            -2.543488950000000e-12 * tc[3]
            +6.659793800000000e-17 * tc[4];
        /*species 7: HCO */
        species[7] =
            +4.011918150000000e+03 * invT
            -7.026170540000000e+00
            -2.772174380000000e+00 * tc[0]
            -2.478477630000000e-03 * tc[1]
            +4.140760216666667e-07 * tc[2]
            -4.909681483333334e-11 * tc[3]
            +2.667543555000000e-15 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +1.118567130000000e+02 * invT
            +2.321087500000001e-01
            -4.017210900000000e+00 * tc[0]
            -1.119910065000000e-03 * tc[1]
            +1.056096916666667e-07 * tc[2]
            -9.520530833333334e-12 * tc[3]
            +5.395426750000000e-16 * tc[4];
        /*species 9: H2O */
        species[9] =
            -3.000429710000000e+04 * invT
            -1.932777610000000e+00
            -3.033992490000000e+00 * tc[0]
            -1.088459020000000e-03 * tc[1]
            +2.734541966666666e-08 * tc[2]
            +8.086832250000000e-12 * tc[3]
            -8.410049600000000e-16 * tc[4];
        /*species 10: CO */
        species[10] =
            -1.415187240000000e+04 * invT
            -5.103502110000000e+00
            -2.715185610000000e+00 * tc[0]
            -1.031263715000000e-03 * tc[1]
            +1.664709618333334e-07 * tc[2]
            -1.917108400000000e-11 * tc[3]
            +1.018238580000000e-15 * tc[4];
        /*species 11: O2 */
        species[11] =
            -1.088457720000000e+03 * invT
            -2.170693450000000e+00
            -3.282537840000000e+00 * tc[0]
            -7.415437700000000e-04 * tc[1]
            +1.263277781666667e-07 * tc[2]
            -1.745587958333333e-11 * tc[3]
            +1.083588970000000e-15 * tc[4];
        /*species 12: H2O2 */
        species[12] =
            -1.786178770000000e+04 * invT
            +1.248846229999999e+00
            -4.165002850000000e+00 * tc[0]
            -2.454158470000000e-03 * tc[1]
            +3.168987083333333e-07 * tc[2]
            -3.093216550000000e-11 * tc[3]
            +1.439541525000000e-15 * tc[4];
        /*species 13: CO2 */
        species[13] =
            -4.875916600000000e+04 * invT
            +1.585822230000000e+00
            -3.857460290000000e+00 * tc[0]
            -2.207185130000000e-03 * tc[1]
            +3.691356733333334e-07 * tc[2]
            -4.362418233333334e-11 * tc[3]
            +2.360420820000000e-15 * tc[4];
    }
    return;
}


/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void helmholtz(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

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
        /*species 2: AR */
        species[2] =
            -7.45375000e+02 * invT
            -2.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 3: N2 */
        species[3] =
            -1.02089990e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
        /*species 4: HE */
        species[4] =
            -7.45375000e+02 * invT
            +5.71276026e-01
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 5: O */
        species[5] =
            +2.91222592e+04 * invT
            +1.16333640e-01
            -3.16826710e+00 * tc[0]
            +1.63965942e-03 * tc[1]
            -1.10717733e-06 * tc[2]
            +5.10672187e-10 * tc[3]
            -1.05632985e-13 * tc[4];
        /*species 6: OH */
        species[6] =
            +3.38153812e+03 * invT
            +3.81573857e+00
            -4.12530561e+00 * tc[0]
            +1.61272470e-03 * tc[1]
            -1.08794115e-06 * tc[2]
            +4.83211369e-10 * tc[3]
            -1.03118689e-13 * tc[4];
        /*species 7: HCO */
        species[7] =
            +3.83956496e+03 * invT
            -1.73186590e-01
            -4.22118584e+00 * tc[0]
            +1.62196266e-03 * tc[1]
            -2.29665743e-06 * tc[2]
            +1.10953411e-09 * tc[3]
            -2.16884432e-13 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +2.94808040e+02 * invT
            -4.14864440e-01
            -4.30179801e+00 * tc[0]
            +2.37456025e-03 * tc[1]
            -3.52638152e-06 * tc[2]
            +2.02303245e-09 * tc[3]
            -4.64612562e-13 * tc[4];
        /*species 9: H2O */
        species[9] =
            -3.02937267e+04 * invT
            +4.04767277e+00
            -4.19864056e+00 * tc[0]
            +1.01821705e-03 * tc[1]
            -1.08673369e-06 * tc[2]
            +4.57330885e-10 * tc[3]
            -8.85989085e-14 * tc[4];
        /*species 10: CO */
        species[10] =
            -1.43440860e+04 * invT
            -9.28875810e-01
            -3.57953347e+00 * tc[0]
            +3.05176840e-04 * tc[1]
            -1.69469055e-07 * tc[2]
            -7.55838237e-11 * tc[3]
            +4.52212249e-14 * tc[4];
        /*species 11: O2 */
        species[11] =
            -1.06394356e+03 * invT
            -8.75219370e-01
            -3.78245636e+00 * tc[0]
            +1.49836708e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774591e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 12: H2O2 */
        species[12] =
            -1.77025821e+04 * invT
            -1.58938050e-01
            -4.27611269e+00 * tc[0]
            +2.71411208e-04 * tc[1]
            -2.78892835e-06 * tc[2]
            +1.79809011e-09 * tc[3]
            -4.31227182e-13 * tc[4];
        /*species 13: CO2 */
        species[13] =
            -4.83719697e+04 * invT
            -8.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
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
        /*species 2: AR */
        species[2] =
            -7.45375000e+02 * invT
            -2.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 3: N2 */
        species[3] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 4: HE */
        species[4] =
            -7.45375000e+02 * invT
            +5.71276026e-01
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 5: O */
        species[5] =
            +2.92175791e+04 * invT
            -3.21491786e+00
            -2.56942078e+00 * tc[0]
            +4.29870569e-05 * tc[1]
            -6.99140982e-09 * tc[2]
            +8.34814992e-13 * tc[3]
            -6.14168455e-17 * tc[4];
        /*species 6: OH */
        species[6] =
            +3.71885774e+03 * invT
            -3.83691187e+00
            -2.86472886e+00 * tc[0]
            -5.28252240e-04 * tc[1]
            +4.31804597e-08 * tc[2]
            -2.54348895e-12 * tc[3]
            +6.65979380e-17 * tc[4];
        /*species 7: HCO */
        species[7] =
            +4.01191815e+03 * invT
            -8.02617054e+00
            -2.77217438e+00 * tc[0]
            -2.47847763e-03 * tc[1]
            +4.14076022e-07 * tc[2]
            -4.90968148e-11 * tc[3]
            +2.66754356e-15 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +1.11856713e+02 * invT
            -7.67891250e-01
            -4.01721090e+00 * tc[0]
            -1.11991006e-03 * tc[1]
            +1.05609692e-07 * tc[2]
            -9.52053083e-12 * tc[3]
            +5.39542675e-16 * tc[4];
        /*species 9: H2O */
        species[9] =
            -3.00042971e+04 * invT
            -2.93277761e+00
            -3.03399249e+00 * tc[0]
            -1.08845902e-03 * tc[1]
            +2.73454197e-08 * tc[2]
            +8.08683225e-12 * tc[3]
            -8.41004960e-16 * tc[4];
        /*species 10: CO */
        species[10] =
            -1.41518724e+04 * invT
            -6.10350211e+00
            -2.71518561e+00 * tc[0]
            -1.03126372e-03 * tc[1]
            +1.66470962e-07 * tc[2]
            -1.91710840e-11 * tc[3]
            +1.01823858e-15 * tc[4];
        /*species 11: O2 */
        species[11] =
            -1.08845772e+03 * invT
            -3.17069345e+00
            -3.28253784e+00 * tc[0]
            -7.41543770e-04 * tc[1]
            +1.26327778e-07 * tc[2]
            -1.74558796e-11 * tc[3]
            +1.08358897e-15 * tc[4];
        /*species 12: H2O2 */
        species[12] =
            -1.78617877e+04 * invT
            +2.48846230e-01
            -4.16500285e+00 * tc[0]
            -2.45415847e-03 * tc[1]
            +3.16898708e-07 * tc[2]
            -3.09321655e-11 * tc[3]
            +1.43954153e-15 * tc[4];
        /*species 13: CO2 */
        species[13] =
            -4.87591660e+04 * invT
            +5.85822230e-01
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
    }
    return;
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cv_R(double * restrict species, double * restrict tc)
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
        /*species 2: AR */
        species[2] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 3: N2 */
        species[3] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 4: HE */
        species[4] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 5: O */
        species[5] =
            +2.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 6: OH */
        species[6] =
            +3.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 7: HCO */
        species[7] =
            +3.22118584e+00
            -3.24392532e-03 * tc[1]
            +1.37799446e-05 * tc[2]
            -1.33144093e-08 * tc[3]
            +4.33768865e-12 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 9: H2O */
        species[9] =
            +3.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 10: CO */
        species[10] =
            +2.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 11: O2 */
        species[11] =
            +2.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 12: H2O2 */
        species[12] =
            +3.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 13: CO2 */
        species[13] =
            +1.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
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
        /*species 2: AR */
        species[2] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 3: N2 */
        species[3] =
            +1.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 4: HE */
        species[4] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 5: O */
        species[5] =
            +1.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 6: OH */
        species[6] =
            +1.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 7: HCO */
        species[7] =
            +1.77217438e+00
            +4.95695526e-03 * tc[1]
            -2.48445613e-06 * tc[2]
            +5.89161778e-10 * tc[3]
            -5.33508711e-14 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 9: H2O */
        species[9] =
            +2.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 10: CO */
        species[10] =
            +1.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 11: O2 */
        species[11] =
            +2.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 12: H2O2 */
        species[12] =
            +3.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 13: CO2 */
        species[13] =
            +2.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cp_R(double * restrict species, double * restrict tc)
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
        /*species 2: AR */
        species[2] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 3: N2 */
        species[3] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 4: HE */
        species[4] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 5: O */
        species[5] =
            +3.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 6: OH */
        species[6] =
            +4.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 7: HCO */
        species[7] =
            +4.22118584e+00
            -3.24392532e-03 * tc[1]
            +1.37799446e-05 * tc[2]
            -1.33144093e-08 * tc[3]
            +4.33768865e-12 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 9: H2O */
        species[9] =
            +4.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 10: CO */
        species[10] =
            +3.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 11: O2 */
        species[11] =
            +3.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 12: H2O2 */
        species[12] =
            +4.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 13: CO2 */
        species[13] =
            +2.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
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
        /*species 2: AR */
        species[2] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 3: N2 */
        species[3] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 4: HE */
        species[4] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 5: O */
        species[5] =
            +2.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 6: OH */
        species[6] =
            +2.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 7: HCO */
        species[7] =
            +2.77217438e+00
            +4.95695526e-03 * tc[1]
            -2.48445613e-06 * tc[2]
            +5.89161778e-10 * tc[3]
            -5.33508711e-14 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 9: H2O */
        species[9] =
            +3.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 10: CO */
        species[10] =
            +2.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 11: O2 */
        species[11] =
            +3.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 12: H2O2 */
        species[12] =
            +4.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 13: CO2 */
        species[13] =
            +3.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesInternalEnergy(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

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
        /*species 2: AR */
        species[2] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 3: N2 */
        species[3] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
        /*species 4: HE */
        species[4] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 5: O */
        species[5] =
            +2.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 6: OH */
        species[6] =
            +3.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.38153812e+03 * invT;
        /*species 7: HCO */
        species[7] =
            +3.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 * invT;
        /*species 8: HO2 */
        species[8] =
            +3.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 9: H2O */
        species[9] =
            +3.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 10: CO */
        species[10] =
            +2.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 11: O2 */
        species[11] =
            +2.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 12: H2O2 */
        species[12] =
            +3.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 * invT;
        /*species 13: CO2 */
        species[13] =
            +1.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
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
        /*species 2: AR */
        species[2] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 3: N2 */
        species[3] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 4: HE */
        species[4] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 5: O */
        species[5] =
            +1.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 6: OH */
        species[6] =
            +1.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.71885774e+03 * invT;
        /*species 7: HCO */
        species[7] =
            +1.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 * invT;
        /*species 8: HO2 */
        species[8] =
            +3.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 9: H2O */
        species[9] =
            +2.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 10: CO */
        species[10] =
            +1.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 11: O2 */
        species[11] =
            +2.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 12: H2O2 */
        species[12] =
            +3.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 * invT;
        /*species 13: CO2 */
        species[13] =
            +2.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEnthalpy(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

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
        /*species 2: AR */
        species[2] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 3: N2 */
        species[3] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
        /*species 4: HE */
        species[4] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 5: O */
        species[5] =
            +3.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 6: OH */
        species[6] =
            +4.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.38153812e+03 * invT;
        /*species 7: HCO */
        species[7] =
            +4.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 * invT;
        /*species 8: HO2 */
        species[8] =
            +4.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 9: H2O */
        species[9] =
            +4.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 10: CO */
        species[10] =
            +3.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 11: O2 */
        species[11] =
            +3.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 12: H2O2 */
        species[12] =
            +4.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 * invT;
        /*species 13: CO2 */
        species[13] =
            +2.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
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
        /*species 2: AR */
        species[2] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 3: N2 */
        species[3] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 4: HE */
        species[4] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 5: O */
        species[5] =
            +2.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 6: OH */
        species[6] =
            +2.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.71885774e+03 * invT;
        /*species 7: HCO */
        species[7] =
            +2.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 * invT;
        /*species 8: HO2 */
        species[8] =
            +4.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 9: H2O */
        species[9] =
            +3.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 10: CO */
        species[10] =
            +2.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 11: O2 */
        species[11] =
            +3.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 12: H2O2 */
        species[12] =
            +4.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 * invT;
        /*species 13: CO2 */
        species[13] =
            +3.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEntropy(double * restrict species, double * restrict tc)
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
        /*species 2: AR */
        species[2] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.36600000e+00 ;
        /*species 3: N2 */
        species[3] =
            +3.29867700e+00 * tc[0]
            +1.40824040e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213500e-13 * tc[4]
            +3.95037200e+00 ;
        /*species 4: HE */
        species[4] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +9.28723974e-01 ;
        /*species 5: O */
        species[5] =
            +3.16826710e+00 * tc[0]
            -3.27931884e-03 * tc[1]
            +3.32153198e-06 * tc[2]
            -2.04268875e-09 * tc[3]
            +5.28164927e-13 * tc[4]
            +2.05193346e+00 ;
        /*species 6: OH */
        species[6] =
            +4.12530561e+00 * tc[0]
            -3.22544939e-03 * tc[1]
            +3.26382346e-06 * tc[2]
            -1.93284548e-09 * tc[3]
            +5.15593447e-13 * tc[4]
            -6.90432960e-01 ;
        /*species 7: HCO */
        species[7] =
            +4.22118584e+00 * tc[0]
            -3.24392532e-03 * tc[1]
            +6.88997230e-06 * tc[2]
            -4.43813643e-09 * tc[3]
            +1.08442216e-12 * tc[4]
            +3.39437243e+00 ;
        /*species 8: HO2 */
        species[8] =
            +4.30179801e+00 * tc[0]
            -4.74912051e-03 * tc[1]
            +1.05791445e-05 * tc[2]
            -8.09212980e-09 * tc[3]
            +2.32306281e-12 * tc[4]
            +3.71666245e+00 ;
        /*species 9: H2O */
        species[9] =
            +4.19864056e+00 * tc[0]
            -2.03643410e-03 * tc[1]
            +3.26020105e-06 * tc[2]
            -1.82932354e-09 * tc[3]
            +4.42994543e-13 * tc[4]
            -8.49032208e-01 ;
        /*species 10: CO */
        species[10] =
            +3.57953347e+00 * tc[0]
            -6.10353680e-04 * tc[1]
            +5.08407165e-07 * tc[2]
            +3.02335295e-10 * tc[3]
            -2.26106125e-13 * tc[4]
            +3.50840928e+00 ;
        /*species 11: O2 */
        species[11] =
            +3.78245636e+00 * tc[0]
            -2.99673416e-03 * tc[1]
            +4.92365101e-06 * tc[2]
            -3.22709836e-09 * tc[3]
            +8.10932092e-13 * tc[4]
            +3.65767573e+00 ;
        /*species 12: H2O2 */
        species[12] =
            +4.27611269e+00 * tc[0]
            -5.42822417e-04 * tc[1]
            +8.36678505e-06 * tc[2]
            -7.19236043e-09 * tc[3]
            +2.15613591e-12 * tc[4]
            +3.43505074e+00 ;
        /*species 13: CO2 */
        species[13] =
            +2.35677352e+00 * tc[0]
            +8.98459677e-03 * tc[1]
            -3.56178134e-06 * tc[2]
            +8.19730073e-10 * tc[3]
            -3.59248870e-14 * tc[4]
            +9.90105222e+00 ;
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
        /*species 2: AR */
        species[2] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.36600000e+00 ;
        /*species 3: N2 */
        species[3] =
            +2.92664000e+00 * tc[0]
            +1.48797680e-03 * tc[1]
            -2.84238000e-07 * tc[2]
            +3.36567933e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
        /*species 4: HE */
        species[4] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +9.28723974e-01 ;
        /*species 5: O */
        species[5] =
            +2.56942078e+00 * tc[0]
            -8.59741137e-05 * tc[1]
            +2.09742295e-08 * tc[2]
            -3.33925997e-12 * tc[3]
            +3.07084227e-16 * tc[4]
            +4.78433864e+00 ;
        /*species 6: OH */
        species[6] =
            +2.86472886e+00 * tc[0]
            +1.05650448e-03 * tc[1]
            -1.29541379e-07 * tc[2]
            +1.01739558e-11 * tc[3]
            -3.32989690e-16 * tc[4]
            +5.70164073e+00 ;
        /*species 7: HCO */
        species[7] =
            +2.77217438e+00 * tc[0]
            +4.95695526e-03 * tc[1]
            -1.24222806e-06 * tc[2]
            +1.96387259e-10 * tc[3]
            -1.33377178e-14 * tc[4]
            +9.79834492e+00 ;
        /*species 8: HO2 */
        species[8] =
            +4.01721090e+00 * tc[0]
            +2.23982013e-03 * tc[1]
            -3.16829075e-07 * tc[2]
            +3.80821233e-11 * tc[3]
            -2.69771337e-15 * tc[4]
            +3.78510215e+00 ;
        /*species 9: H2O */
        species[9] =
            +3.03399249e+00 * tc[0]
            +2.17691804e-03 * tc[1]
            -8.20362590e-08 * tc[2]
            -3.23473290e-11 * tc[3]
            +4.20502480e-15 * tc[4]
            +4.96677010e+00 ;
        /*species 10: CO */
        species[10] =
            +2.71518561e+00 * tc[0]
            +2.06252743e-03 * tc[1]
            -4.99412886e-07 * tc[2]
            +7.66843360e-11 * tc[3]
            -5.09119290e-15 * tc[4]
            +7.81868772e+00 ;
        /*species 11: O2 */
        species[11] =
            +3.28253784e+00 * tc[0]
            +1.48308754e-03 * tc[1]
            -3.78983334e-07 * tc[2]
            +6.98235183e-11 * tc[3]
            -5.41794485e-15 * tc[4]
            +5.45323129e+00 ;
        /*species 12: H2O2 */
        species[12] =
            +4.16500285e+00 * tc[0]
            +4.90831694e-03 * tc[1]
            -9.50696125e-07 * tc[2]
            +1.23728662e-10 * tc[3]
            -7.19770763e-15 * tc[4]
            +2.91615662e+00 ;
        /*species 13: CO2 */
        species[13] =
            +3.85746029e+00 * tc[0]
            +4.41437026e-03 * tc[1]
            -1.10740702e-06 * tc[2]
            +1.74496729e-10 * tc[3]
            -1.18021041e-14 * tc[4]
            +2.27163806e+00 ;
    }
    return;
}


/*save molecular weights into array */
void molecularWeight(double * restrict wt)
{
    wt[0] = 2.015940; /*H2 */
    wt[1] = 1.007970; /*H */
    wt[2] = 39.948000; /*AR */
    wt[3] = 28.013400; /*N2 */
    wt[4] = 4.002600; /*HE */
    wt[5] = 15.999400; /*O */
    wt[6] = 17.007370; /*OH */
    wt[7] = 29.018520; /*HCO */
    wt[8] = 33.006770; /*HO2 */
    wt[9] = 18.015340; /*H2O */
    wt[10] = 28.010550; /*CO */
    wt[11] = 31.998800; /*O2 */
    wt[12] = 34.014740; /*H2O2 */
    wt[13] = 44.009950; /*CO2 */

    return;
}


/*save atomic weights into array */
void atomicWeight(double * restrict awt)
{
    awt[0] = 15.999400; /*O */
    awt[1] = 1.007970; /*H */
    awt[2] = 12.011150; /*C */
    awt[3] = 14.006700; /*N */
    awt[4] = 39.948000; /*AR */
    awt[5] = 4.002600; /*HE */

    return;
}
/* get temperature given internal energy in mass units and mass fracs */
void GET_T_GIVEN_EY(double * restrict e, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int * ierr)
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
  *LENIMC =           59;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENRMC EGTRANSETLENRMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENRMC egtransetlenrmc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENRMC egtransetlenrmc_
#endif
void egtransetLENRMC(int* LENRMC) {
  *LENRMC =         4270;}
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
  *KK =           14;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLITE EGTRANSETNLITE
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLITE egtransetnlite
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLITE egtransetnlite_
#endif
void egtransetNLITE(int* NLITE) {
  *NLITE =            3;}
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
  WT[           1] =   0.1007969975471497E+01;
  WT[           2] =   0.3994800186157227E+02;
  WT[           3] =   0.2801339912414551E+02;
  WT[           4] =   0.4002600193023682E+01;
  WT[           5] =   0.1599940013885498E+02;
  WT[           6] =   0.1700737011432648E+02;
  WT[           7] =   0.2901852047443390E+02;
  WT[           8] =   0.3300677025318146E+02;
  WT[           9] =   0.1801534008979797E+02;
  WT[          10] =   0.2801055049896240E+02;
  WT[          11] =   0.3199880027770996E+02;
  WT[          12] =   0.3401474022865295E+02;
  WT[          13] =   0.4400995063781738E+02;
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
  EPS[           1] =   0.1450000000000000E+03;
  EPS[           2] =   0.1365000000000000E+03;
  EPS[           3] =   0.9753000000000000E+02;
  EPS[           4] =   0.1020000000000000E+02;
  EPS[           5] =   0.8000000000000000E+02;
  EPS[           6] =   0.8000000000000000E+02;
  EPS[           7] =   0.4980000000000000E+03;
  EPS[           8] =   0.1074000000000000E+03;
  EPS[           9] =   0.5724000000000000E+03;
  EPS[          10] =   0.9809999999999999E+02;
  EPS[          11] =   0.1074000000000000E+03;
  EPS[          12] =   0.1074000000000000E+03;
  EPS[          13] =   0.2440000000000000E+03;
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
  SIG[           1] =   0.2050000000000000E+01;
  SIG[           2] =   0.3330000000000000E+01;
  SIG[           3] =   0.3621000000000000E+01;
  SIG[           4] =   0.2576000000000000E+01;
  SIG[           5] =   0.2750000000000000E+01;
  SIG[           6] =   0.2750000000000000E+01;
  SIG[           7] =   0.3590000000000000E+01;
  SIG[           8] =   0.3458000000000000E+01;
  SIG[           9] =   0.2605000000000000E+01;
  SIG[          10] =   0.3650000000000000E+01;
  SIG[          11] =   0.3458000000000000E+01;
  SIG[          12] =   0.3458000000000000E+01;
  SIG[          13] =   0.3763000000000000E+01;
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
  DIP[           2] =   0.0000000000000000E+00;
  DIP[           3] =   0.0000000000000000E+00;
  DIP[           4] =   0.0000000000000000E+00;
  DIP[           5] =   0.0000000000000000E+00;
  DIP[           6] =   0.0000000000000000E+00;
  DIP[           7] =   0.0000000000000000E+00;
  DIP[           8] =   0.0000000000000000E+00;
  DIP[           9] =   0.1844000000000000E+01;
  DIP[          10] =   0.0000000000000000E+00;
  DIP[          11] =   0.0000000000000000E+00;
  DIP[          12] =   0.0000000000000000E+00;
  DIP[          13] =   0.0000000000000000E+00;
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
  POL[           1] =   0.0000000000000000E+00;
  POL[           2] =   0.0000000000000000E+00;
  POL[           3] =   0.1760000000000000E+01;
  POL[           4] =   0.0000000000000000E+00;
  POL[           5] =   0.0000000000000000E+00;
  POL[           6] =   0.0000000000000000E+00;
  POL[           7] =   0.0000000000000000E+00;
  POL[           8] =   0.0000000000000000E+00;
  POL[           9] =   0.0000000000000000E+00;
  POL[          10] =   0.1950000000000000E+01;
  POL[          11] =   0.1600000000000000E+01;
  POL[          12] =   0.0000000000000000E+00;
  POL[          13] =   0.2650000000000000E+01;
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
  ZROT[           1] =   0.0000000000000000E+00;
  ZROT[           2] =   0.0000000000000000E+00;
  ZROT[           3] =   0.4000000000000000E+01;
  ZROT[           4] =   0.0000000000000000E+00;
  ZROT[           5] =   0.0000000000000000E+00;
  ZROT[           6] =   0.0000000000000000E+00;
  ZROT[           7] =   0.0000000000000000E+00;
  ZROT[           8] =   0.1000000000000000E+01;
  ZROT[           9] =   0.4000000000000000E+01;
  ZROT[          10] =   0.1800000000000000E+01;
  ZROT[          11] =   0.3800000000000000E+01;
  ZROT[          12] =   0.3800000000000000E+01;
  ZROT[          13] =   0.2100000000000000E+01;
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
  NLIN[           1] =            0;
  NLIN[           2] =            0;
  NLIN[           3] =            1;
  NLIN[           4] =            0;
  NLIN[           5] =            0;
  NLIN[           6] =            1;
  NLIN[           7] =            2;
  NLIN[           8] =            2;
  NLIN[           9] =            2;
  NLIN[          10] =            1;
  NLIN[          11] =            1;
  NLIN[          12] =            2;
  NLIN[          13] =            1;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFLAM EGTRANSETCOFLAM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFLAM egtransetcoflam
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFLAM egtransetcoflam_
#endif
void egtransetCOFLAM(double* COFLAM) {
  COFLAM[           0] =   0.1158096158594203E+02;
  COFLAM[           1] =  -0.1520398050330763E+01;
  COFLAM[           2] =   0.2722793392212857E+00;
  COFLAM[           3] =  -0.1031257389909708E-01;
  COFLAM[           4] =  -0.3253788579319486E+00;
  COFLAM[           5] =   0.3416397862210875E+01;
  COFLAM[           6] =  -0.3631104489846262E+00;
  COFLAM[           7] =   0.1585986202702918E-01;
  COFLAM[           8] =  -0.2731910414923424E+01;
  COFLAM[           9] =   0.3271507387309670E+01;
  COFLAM[          10] =  -0.3454038971012363E+00;
  COFLAM[          11] =   0.1513951680124476E-01;
  COFLAM[          12] =   0.1154286034252058E+02;
  COFLAM[          13] =  -0.2911509179774523E+01;
  COFLAM[          14] =   0.5546559973428081E+00;
  COFLAM[          15] =  -0.2750092774174382E-01;
  COFLAM[          16] =   0.7017839888449308E+01;
  COFLAM[          17] =   0.2176250912976264E+00;
  COFLAM[          18] =   0.5619507686207265E-01;
  COFLAM[          19] =  -0.2368322714922462E-02;
  COFLAM[          20] =   0.1969666526185763E+01;
  COFLAM[          21] =   0.1801396425547300E+01;
  COFLAM[          22] =  -0.1549102264837655E+00;
  COFLAM[          23] =   0.6908759721669425E-02;
  COFLAM[          24] =   0.1605449323860262E+02;
  COFLAM[          25] =  -0.4103244048729492E+01;
  COFLAM[          26] =   0.6631513498156665E+00;
  COFLAM[          27] =  -0.2977137989389932E-01;
  COFLAM[          28] =   0.1709383690534594E+01;
  COFLAM[          29] =  -0.1709779525397148E+00;
  COFLAM[          30] =   0.3320736523343362E+00;
  COFLAM[          31] =  -0.2309480476620785E-01;
  COFLAM[          32] =   0.5546401577805573E+00;
  COFLAM[          33] =   0.1591057931808813E+01;
  COFLAM[          34] =  -0.5282455808284543E-01;
  COFLAM[          35] =   0.4072391521895438E-03;
  COFLAM[          36] =   0.2277162937237345E+02;
  COFLAM[          37] =  -0.8704976392893332E+01;
  COFLAM[          38] =   0.1490812066077175E+01;
  COFLAM[          39] =  -0.7406350595204310E-01;
  COFLAM[          40] =   0.9917329347278134E+01;
  COFLAM[          41] =  -0.2280409546557303E+01;
  COFLAM[          42] =   0.4727591838750459E+00;
  COFLAM[          43] =  -0.2399081985632417E-01;
  COFLAM[          44] =  -0.3013713897782378E+01;
  COFLAM[          45] =   0.3375529164503225E+01;
  COFLAM[          46] =  -0.3433040143356826E+00;
  COFLAM[          47] =   0.1509986872028326E-01;
  COFLAM[          48] =   0.6252374647013528E+00;
  COFLAM[          49] =   0.1431934364158212E+01;
  COFLAM[          50] =   0.1750712355480039E-02;
  COFLAM[          51] =  -0.3554756889670176E-02;
  COFLAM[          52] =  -0.1239079949009015E+02;
  COFLAM[          53] =   0.6341842149327256E+01;
  COFLAM[          54] =  -0.6370214552008658E+00;
  COFLAM[          55] =   0.2372277065304921E-01;
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
  COFETA[           4] =  -0.1987529414716897E+02;
  COFETA[           5] =   0.3416397862210910E+01;
  COFETA[           6] =  -0.3631104489846305E+00;
  COFETA[           7] =   0.1585986202702936E-01;
  COFETA[           8] =  -0.1860218543214680E+02;
  COFETA[           9] =   0.3271507387309783E+01;
  COFETA[          10] =  -0.3454038971012522E+00;
  COFETA[          11] =   0.1513951680124550E-01;
  COFETA[          12] =  -0.1626172843728175E+02;
  COFETA[          13] =   0.2251740453876138E+01;
  COFETA[          14] =  -0.2138340893699383E+00;
  COFETA[          15] =   0.9477823154448534E-02;
  COFETA[          16] =  -0.1115306958556883E+02;
  COFETA[          17] =   0.2176250912977139E+00;
  COFETA[          18] =   0.5619507686206007E-01;
  COFETA[          19] =  -0.2368322714921862E-02;
  COFETA[          20] =  -0.1481563591580279E+02;
  COFETA[          21] =   0.1801396425547327E+01;
  COFETA[          22] =  -0.1549102264837694E+00;
  COFETA[          23] =   0.6908759721669611E-02;
  COFETA[          24] =  -0.1478508813778873E+02;
  COFETA[          25] =   0.1801396425547305E+01;
  COFETA[          26] =  -0.1549102264837663E+00;
  COFETA[          27] =   0.6908759721669463E-02;
  COFETA[          28] =  -0.2114891508052214E+02;
  COFETA[          29] =   0.3277849124444586E+01;
  COFETA[          30] =  -0.2525897277867938E+00;
  COFETA[          31] =   0.7417162247591944E-02;
  COFETA[          32] =  -0.1679529396430702E+02;
  COFETA[          33] =   0.2522528725237947E+01;
  COFETA[          34] =  -0.2490712798240481E+00;
  COFETA[          35] =   0.1100615806447173E-01;
  COFETA[          36] =  -0.1187800764707307E+02;
  COFETA[          37] =  -0.7882519505482808E+00;
  COFETA[          38] =   0.3341408170058896E+00;
  COFETA[          39] =  -0.1986366361647418E-01;
  COFETA[          40] =  -0.1631251619137886E+02;
  COFETA[          41] =   0.2264906486512875E+01;
  COFETA[          42] =  -0.2155365004850580E+00;
  COFETA[          43] =   0.9551179268412450E-02;
  COFETA[          44] =  -0.1681080110872507E+02;
  COFETA[          45] =   0.2522528725237932E+01;
  COFETA[          46] =  -0.2490712798240462E+00;
  COFETA[          47] =   0.1100615806447165E-01;
  COFETA[          48] =  -0.1678025333071109E+02;
  COFETA[          49] =   0.2522528725237946E+01;
  COFETA[          50] =  -0.2490712798240483E+00;
  COFETA[          51] =   0.1100615806447175E-01;
  COFETA[          52] =  -0.2364265995035323E+02;
  COFETA[          53] =   0.4983971927098391E+01;
  COFETA[          54] =  -0.5507551343971877E+00;
  COFETA[          55] =   0.2334595082438979E-01;
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
  COFD[           4] =  -0.1168685516750758E+02;
  COFD[           5] =   0.2883726644585007E+01;
  COFD[           6] =  -0.1637777890869061E+00;
  COFD[           7] =   0.7265872474830623E-02;
  COFD[           8] =  -0.1272604230538407E+02;
  COFD[           9] =   0.2921765441423414E+01;
  COFD[          10] =  -0.1699072228926505E+00;
  COFD[          11] =   0.7582709264066988E-02;
  COFD[          12] =  -0.1217511489930491E+02;
  COFD[          13] =   0.2679286720576265E+01;
  COFD[          14] =  -0.1373997446113497E+00;
  COFD[          15] =   0.6125379451516252E-02;
  COFD[          16] =  -0.9853589402094190E+01;
  COFD[          17] =   0.2052962092552261E+01;
  COFD[          18] =  -0.5895182333307221E-01;
  COFD[          19] =   0.2891951845562504E-02;
  COFD[          20] =  -0.1060109579664357E+02;
  COFD[          21] =   0.2157129122185176E+01;
  COFD[          22] =  -0.6524729611183182E-01;
  COFD[          23] =   0.2809597810849289E-02;
  COFD[          24] =  -0.1060412419807472E+02;
  COFD[          25] =   0.2157123706881953E+01;
  COFD[          26] =  -0.6524690383610432E-01;
  COFD[          27] =   0.2809593789707914E-02;
  COFD[          28] =  -0.1752884957692982E+02;
  COFD[          29] =   0.4691328454817521E+01;
  COFD[          30] =  -0.3953771017966438E+00;
  COFD[          31] =   0.1716337377732627E-01;
  COFD[          32] =  -0.1229066279603360E+02;
  COFD[          33] =   0.2741057082889793E+01;
  COFD[          34] =  -0.1457627022193994E+00;
  COFD[          35] =   0.6504615244946633E-02;
  COFD[          36] =  -0.1787213075631945E+02;
  COFD[          37] =   0.4905112324365878E+01;
  COFD[          38] =  -0.4173872143518051E+00;
  COFD[          39] =   0.1788845489573668E-01;
  COFD[          40] =  -0.1220110902734887E+02;
  COFD[          41] =   0.2686277536441436E+01;
  COFD[          42] =  -0.1383620679476532E+00;
  COFD[          43] =   0.6169543276698330E-02;
  COFD[          44] =  -0.1228690337278460E+02;
  COFD[          45] =   0.2739817607547901E+01;
  COFD[          46] =  -0.1455907617836413E+00;
  COFD[          47] =   0.6496697467578658E-02;
  COFD[          48] =  -0.1229422250798061E+02;
  COFD[          49] =   0.2742232459895579E+01;
  COFD[          50] =  -0.1459257489141948E+00;
  COFD[          51] =   0.6512123400236565E-02;
  COFD[          52] =  -0.1530635401605573E+02;
  COFD[          53] =   0.3873896336632167E+01;
  COFD[          54] =  -0.2956113676790514E+00;
  COFD[          55] =   0.1311329899072374E-01;
  COFD[          56] =  -0.1168685516750758E+02;
  COFD[          57] =   0.2883726644585007E+01;
  COFD[          58] =  -0.1637777890869061E+00;
  COFD[          59] =   0.7265872474830623E-02;
  COFD[          60] =  -0.1476537827842222E+02;
  COFD[          61] =   0.4195448036781197E+01;
  COFD[          62] =  -0.3279615676730872E+00;
  COFD[          63] =   0.1412388797372668E-01;
  COFD[          64] =  -0.1715574553730659E+02;
  COFD[          65] =   0.4798805397829151E+01;
  COFD[          66] =  -0.4058992577741924E+00;
  COFD[          67] =   0.1747030892545845E-01;
  COFD[          68] =  -0.1676897875889624E+02;
  COFD[          69] =   0.4677755283143664E+01;
  COFD[          70] =  -0.3974479445519867E+00;
  COFD[          71] =   0.1741111717763365E-01;
  COFD[          72] =  -0.9735381521569142E+01;
  COFD[          73] =   0.2191202791663548E+01;
  COFD[          74] =  -0.7485252920875622E-01;
  COFD[          75] =   0.3470152548075659E-02;
  COFD[          76] =  -0.1500301598770682E+02;
  COFD[          77] =   0.4131920691035950E+01;
  COFD[          78] =  -0.3275337549968449E+00;
  COFD[          79] =   0.1442760802547180E-01;
  COFD[          80] =  -0.1502025959280452E+02;
  COFD[          81] =   0.4138327950805945E+01;
  COFD[          82] =  -0.3284006010219734E+00;
  COFD[          83] =   0.1446660099762999E-01;
  COFD[          84] =  -0.1633601357384808E+02;
  COFD[          85] =   0.3987313912522149E+01;
  COFD[          86] =  -0.2492282800992027E+00;
  COFD[          87] =   0.8598651953545469E-02;
  COFD[          88] =  -0.1720445512239349E+02;
  COFD[          89] =   0.4866330968714804E+01;
  COFD[          90] =  -0.4220901924876364E+00;
  COFD[          91] =   0.1848595807162496E-01;
  COFD[          92] =  -0.1695183404426207E+02;
  COFD[          93] =   0.4416203740445888E+01;
  COFD[          94] =  -0.3114890192543383E+00;
  COFD[          95] =   0.1159663710619628E-01;
  COFD[          96] =  -0.1683623819222381E+02;
  COFD[          97] =   0.4700302341393257E+01;
  COFD[          98] =  -0.4004143455401054E+00;
  COFD[          99] =   0.1754115750370506E-01;
  COFD[         100] =  -0.1719486110604737E+02;
  COFD[         101] =   0.4862600572193149E+01;
  COFD[         102] =  -0.4215909171603991E+00;
  COFD[         103] =   0.1846373197001472E-01;
  COFD[         104] =  -0.1721362788263434E+02;
  COFD[         105] =   0.4869900420159440E+01;
  COFD[         106] =  -0.4225679381202458E+00;
  COFD[         107] =   0.1850722613298993E-01;
  COFD[         108] =  -0.1788417846762182E+02;
  COFD[         109] =   0.4849004827267381E+01;
  COFD[         110] =  -0.3925012395929342E+00;
  COFD[         111] =   0.1605802238455311E-01;
  COFD[         112] =  -0.1272604230538407E+02;
  COFD[         113] =   0.2921765441423414E+01;
  COFD[         114] =  -0.1699072228926505E+00;
  COFD[         115] =   0.7582709264066988E-02;
  COFD[         116] =  -0.1715574553730659E+02;
  COFD[         117] =   0.4798805397829151E+01;
  COFD[         118] =  -0.4058992577741924E+00;
  COFD[         119] =   0.1747030892545845E-01;
  COFD[         120] =  -0.1736905261068727E+02;
  COFD[         121] =   0.4134132180231587E+01;
  COFD[         122] =  -0.3219126846923269E+00;
  COFD[         123] =   0.1394110637546756E-01;
  COFD[         124] =  -0.1621324217553288E+02;
  COFD[         125] =   0.3714519551995751E+01;
  COFD[         126] =  -0.2693102767279618E+00;
  COFD[         127] =   0.1173858202719274E-01;
  COFD[         128] =  -0.1075685106039131E+02;
  COFD[         129] =   0.2105974066569993E+01;
  COFD[         130] =  -0.6310245347805579E-01;
  COFD[         131] =   0.2931474019944176E-02;
  COFD[         132] =  -0.1555067129792547E+02;
  COFD[         133] =   0.3655895323186071E+01;
  COFD[         134] =  -0.2636354398546018E+00;
  COFD[         135] =   0.1157316255067020E-01;
  COFD[         136] =  -0.1553803545411258E+02;
  COFD[         137] =   0.3641692220020330E+01;
  COFD[         138] =  -0.2617120813083521E+00;
  COFD[         139] =   0.1148659824340386E-01;
  COFD[         140] =  -0.2062980903668024E+02;
  COFD[         141] =   0.5158372106011019E+01;
  COFD[         142] =  -0.4299697617186866E+00;
  COFD[         143] =   0.1763220992933028E-01;
  COFD[         144] =  -0.1647873164232925E+02;
  COFD[         145] =   0.3810732050841141E+01;
  COFD[         146] =  -0.2812626132670425E+00;
  COFD[         147] =   0.1223411282905946E-01;
  COFD[         148] =  -0.2001704718429840E+02;
  COFD[         149] =   0.5030685873326798E+01;
  COFD[         150] =  -0.4055157556084015E+00;
  COFD[         151] =   0.1625731432523394E-01;
  COFD[         152] =  -0.1623589255135331E+02;
  COFD[         153] =   0.3719609230867837E+01;
  COFD[         154] =  -0.2699320567627367E+00;
  COFD[         155] =   0.1176389135749383E-01;
  COFD[         156] =  -0.1647508018746298E+02;
  COFD[         157] =   0.3812741716466225E+01;
  COFD[         158] =  -0.2815324383266347E+00;
  COFD[         159] =   0.1224615497227192E-01;
  COFD[         160] =  -0.1648290476852054E+02;
  COFD[         161] =   0.3809087920338710E+01;
  COFD[         162] =  -0.2810417778226483E+00;
  COFD[         163] =   0.1222425349132379E-01;
  COFD[         164] =  -0.1912503774058016E+02;
  COFD[         165] =   0.4674260376378655E+01;
  COFD[         166] =  -0.3838726146843580E+00;
  COFD[         167] =   0.1629606485550395E-01;
  COFD[         168] =  -0.1217511489930491E+02;
  COFD[         169] =   0.2679286720576265E+01;
  COFD[         170] =  -0.1373997446113497E+00;
  COFD[         171] =   0.6125379451516252E-02;
  COFD[         172] =  -0.1676897875889624E+02;
  COFD[         173] =   0.4677755283143664E+01;
  COFD[         174] =  -0.3974479445519867E+00;
  COFD[         175] =   0.1741111717763365E-01;
  COFD[         176] =  -0.1621324217553288E+02;
  COFD[         177] =   0.3714519551995751E+01;
  COFD[         178] =  -0.2693102767279618E+00;
  COFD[         179] =   0.1173858202719274E-01;
  COFD[         180] =  -0.1521415383275665E+02;
  COFD[         181] =   0.3348053449783496E+01;
  COFD[         182] =  -0.2233657260417595E+00;
  COFD[         183] =   0.9817787279109837E-02;
  COFD[         184] =  -0.1014817848789097E+02;
  COFD[         185] =   0.1823570562174538E+01;
  COFD[         186] =  -0.2296236173006269E-01;
  COFD[         187] =   0.1036389912916652E-02;
  COFD[         188] =  -0.1442662508885140E+02;
  COFD[         189] =   0.3216422191096490E+01;
  COFD[         190] =  -0.2069561315576079E+00;
  COFD[         191] =   0.9135288598631700E-02;
  COFD[         192] =  -0.1442063092189782E+02;
  COFD[         193] =   0.3205844724257626E+01;
  COFD[         194] =  -0.2055144159347264E+00;
  COFD[         195] =   0.9070008933822527E-02;
  COFD[         196] =  -0.1997245750726936E+02;
  COFD[         197] =   0.4990620273390539E+01;
  COFD[         198] =  -0.4173575849700564E+00;
  COFD[         199] =   0.1744556746387847E-01;
  COFD[         200] =  -0.1545394137084459E+02;
  COFD[         201] =   0.3436021618273004E+01;
  COFD[         202] =  -0.2341477712817018E+00;
  COFD[         203] =   0.1025708506984392E-01;
  COFD[         204] =  -0.2025975101746627E+02;
  COFD[         205] =   0.5176212790334279E+01;
  COFD[         206] =  -0.4308686680573706E+00;
  COFD[         207] =   0.1761066268522406E-01;
  COFD[         208] =  -0.1523800917587685E+02;
  COFD[         209] =   0.3353970687690228E+01;
  COFD[         210] =  -0.2241227404004539E+00;
  COFD[         211] =   0.9850105182191860E-02;
  COFD[         212] =  -0.1544409203507008E+02;
  COFD[         213] =   0.3434913447661099E+01;
  COFD[         214] =  -0.2339977102148624E+00;
  COFD[         215] =   0.1025033359000777E-01;
  COFD[         216] =  -0.1546403419937241E+02;
  COFD[         217] =   0.3437367417375052E+01;
  COFD[         218] =  -0.2343299630280416E+00;
  COFD[         219] =   0.1026528034463288E-01;
  COFD[         220] =  -0.1811180157901078E+02;
  COFD[         221] =   0.4336077093491970E+01;
  COFD[         222] =  -0.3444079746717316E+00;
  COFD[         223] =   0.1476188131006474E-01;
  COFD[         224] =  -0.9853589402094190E+01;
  COFD[         225] =   0.2052962092552261E+01;
  COFD[         226] =  -0.5895182333307221E-01;
  COFD[         227] =   0.2891951845562504E-02;
  COFD[         228] =  -0.9735381521569142E+01;
  COFD[         229] =   0.2191202791663548E+01;
  COFD[         230] =  -0.7485252920875622E-01;
  COFD[         231] =   0.3470152548075659E-02;
  COFD[         232] =  -0.1075685106039131E+02;
  COFD[         233] =   0.2105974066569993E+01;
  COFD[         234] =  -0.6310245347805579E-01;
  COFD[         235] =   0.2931474019944176E-02;
  COFD[         236] =  -0.1014817848789097E+02;
  COFD[         237] =   0.1823570562174538E+01;
  COFD[         238] =  -0.2296236173006269E-01;
  COFD[         239] =   0.1036389912916652E-02;
  COFD[         240] =  -0.7722235809097614E+01;
  COFD[         241] =   0.1138839639005951E+01;
  COFD[         242] =   0.7226014852152088E-01;
  COFD[         243] =  -0.3322612720077926E-02;
  COFD[         244] =  -0.9681129070187980E+01;
  COFD[         245] =   0.1775605691412164E+01;
  COFD[         246] =  -0.1628295103092079E-01;
  COFD[         247] =   0.7262086859513526E-03;
  COFD[         248] =  -0.9686347735719096E+01;
  COFD[         249] =   0.1775508180457379E+01;
  COFD[         250] =  -0.1627025530906716E-01;
  COFD[         251] =   0.7256647582828849E-03;
  COFD[         252] =  -0.1303387396358854E+02;
  COFD[         253] =   0.2940290699291584E+01;
  COFD[         254] =  -0.1728777131903053E+00;
  COFD[         255] =   0.7735387047566533E-02;
  COFD[         256] =  -0.1028179959593585E+02;
  COFD[         257] =   0.1896257251270889E+01;
  COFD[         258] =  -0.3329823450147470E-01;
  COFD[         259] =   0.1524623974194515E-02;
  COFD[         260] =  -0.1286980431258231E+02;
  COFD[         261] =   0.3009638621722941E+01;
  COFD[         262] =  -0.1808900918568948E+00;
  COFD[         263] =   0.8039983754952261E-02;
  COFD[         264] =  -0.1015561153728631E+02;
  COFD[         265] =   0.1822574904230457E+01;
  COFD[         266] =  -0.2281022330021956E-01;
  COFD[         267] =   0.1028887956167319E-02;
  COFD[         268] =  -0.1028043354360960E+02;
  COFD[         269] =   0.1896303471491212E+01;
  COFD[         270] =  -0.3330438429816188E-01;
  COFD[         271] =   0.1524890767209287E-02;
  COFD[         272] =  -0.1028308834746900E+02;
  COFD[         273] =   0.1896214006095980E+01;
  COFD[         274] =  -0.3329247904451228E-01;
  COFD[         275] =   0.1524374377170872E-02;
  COFD[         276] =  -0.1084949698568986E+02;
  COFD[         277] =   0.2037615396660156E+01;
  COFD[         278] =  -0.4955359664787870E-01;
  COFD[         279] =   0.2122039387077741E-02;
  COFD[         280] =  -0.1060109579664357E+02;
  COFD[         281] =   0.2157129122185176E+01;
  COFD[         282] =  -0.6524729611183182E-01;
  COFD[         283] =   0.2809597810849289E-02;
  COFD[         284] =  -0.1500301598770682E+02;
  COFD[         285] =   0.4131920691035950E+01;
  COFD[         286] =  -0.3275337549968449E+00;
  COFD[         287] =   0.1442760802547180E-01;
  COFD[         288] =  -0.1555067129792547E+02;
  COFD[         289] =   0.3655895323186071E+01;
  COFD[         290] =  -0.2636354398546018E+00;
  COFD[         291] =   0.1157316255067020E-01;
  COFD[         292] =  -0.1442662508885140E+02;
  COFD[         293] =   0.3216422191096490E+01;
  COFD[         294] =  -0.2069561315576079E+00;
  COFD[         295] =   0.9135288598631700E-02;
  COFD[         296] =  -0.9681129070187980E+01;
  COFD[         297] =   0.1775605691412164E+01;
  COFD[         298] =  -0.1628295103092079E-01;
  COFD[         299] =   0.7262086859513526E-03;
  COFD[         300] =  -0.1329503778398344E+02;
  COFD[         301] =   0.2938989816956702E+01;
  COFD[         302] =  -0.1706012465018490E+00;
  COFD[         303] =   0.7547358484820122E-02;
  COFD[         304] =  -0.1331106228772027E+02;
  COFD[         305] =   0.2939408772392007E+01;
  COFD[         306] =  -0.1706589698617190E+00;
  COFD[         307] =   0.7549998035548285E-02;
  COFD[         308] =  -0.1908503263157613E+02;
  COFD[         309] =   0.4839712718287958E+01;
  COFD[         310] =  -0.4003349576619368E+00;
  COFD[         311] =   0.1680186521574959E-01;
  COFD[         312] =  -0.1476396639560408E+02;
  COFD[         313] =   0.3356474034011327E+01;
  COFD[         314] =  -0.2254105976116209E+00;
  COFD[         315] =   0.9946130954817655E-02;
  COFD[         316] =  -0.1894546586601941E+02;
  COFD[         317] =   0.4935720676293722E+01;
  COFD[         318] =  -0.4114057447467022E+00;
  COFD[         319] =   0.1723315986098057E-01;
  COFD[         320] =  -0.1445980203722702E+02;
  COFD[         321] =   0.3225883257074325E+01;
  COFD[         322] =  -0.2082125398549377E+00;
  COFD[         323] =   0.9190899948945185E-02;
  COFD[         324] =  -0.1474363832488901E+02;
  COFD[         325] =   0.3350090901251298E+01;
  COFD[         326] =  -0.2245382243795846E+00;
  COFD[         327] =   0.9906528312288899E-02;
  COFD[         328] =  -0.1478376315874751E+02;
  COFD[         329] =   0.3362741024796286E+01;
  COFD[         330] =  -0.2262670743508222E+00;
  COFD[         331] =   0.9985011202484037E-02;
  COFD[         332] =  -0.1758398460398111E+02;
  COFD[         333] =   0.4332475641237926E+01;
  COFD[         334] =  -0.3464887767008653E+00;
  COFD[         335] =   0.1495283324430825E-01;
  COFD[         336] =  -0.1060412419807472E+02;
  COFD[         337] =   0.2157123706881953E+01;
  COFD[         338] =  -0.6524690383610432E-01;
  COFD[         339] =   0.2809593789707914E-02;
  COFD[         340] =  -0.1502025959280452E+02;
  COFD[         341] =   0.4138327950805945E+01;
  COFD[         342] =  -0.3284006010219734E+00;
  COFD[         343] =   0.1446660099762999E-01;
  COFD[         344] =  -0.1553803545411258E+02;
  COFD[         345] =   0.3641692220020330E+01;
  COFD[         346] =  -0.2617120813083521E+00;
  COFD[         347] =   0.1148659824340386E-01;
  COFD[         348] =  -0.1442063092189782E+02;
  COFD[         349] =   0.3205844724257626E+01;
  COFD[         350] =  -0.2055144159347264E+00;
  COFD[         351] =   0.9070008933822527E-02;
  COFD[         352] =  -0.9686347735719096E+01;
  COFD[         353] =   0.1775508180457379E+01;
  COFD[         354] =  -0.1627025530906716E-01;
  COFD[         355] =   0.7256647582828849E-03;
  COFD[         356] =  -0.1331106228772027E+02;
  COFD[         357] =   0.2939408772392007E+01;
  COFD[         358] =  -0.1706589698617190E+00;
  COFD[         359] =   0.7549998035548285E-02;
  COFD[         360] =  -0.1332558556199739E+02;
  COFD[         361] =   0.2938989816956671E+01;
  COFD[         362] =  -0.1706012465018444E+00;
  COFD[         363] =   0.7547358484819904E-02;
  COFD[         364] =  -0.1910414995755510E+02;
  COFD[         365] =   0.4841043235938669E+01;
  COFD[         366] =  -0.4007327593260240E+00;
  COFD[         367] =   0.1682891976067088E-01;
  COFD[         368] =  -0.1475457475365462E+02;
  COFD[         369] =   0.3343986593397501E+01;
  COFD[         370] =  -0.2237039347662408E+00;
  COFD[         371] =   0.9868653787965939E-02;
  COFD[         372] =  -0.1896004007077217E+02;
  COFD[         373] =   0.4935377774711531E+01;
  COFD[         374] =  -0.4113926173034602E+00;
  COFD[         375] =   0.1723402811988101E-01;
  COFD[         376] =  -0.1445334512020511E+02;
  COFD[         377] =   0.3215111284442823E+01;
  COFD[         378] =  -0.2067440027231415E+00;
  COFD[         379] =   0.9124392199487649E-02;
  COFD[         380] =  -0.1473449002110577E+02;
  COFD[         381] =   0.3337793989426398E+01;
  COFD[         382] =  -0.2228575541864579E+00;
  COFD[         383] =   0.9830229599517820E-02;
  COFD[         384] =  -0.1477418610290300E+02;
  COFD[         385] =   0.3350090901251285E+01;
  COFD[         386] =  -0.2245382243795827E+00;
  COFD[         387] =   0.9906528312288795E-02;
  COFD[         388] =  -0.1757223862273510E+02;
  COFD[         389] =   0.4319080008767561E+01;
  COFD[         390] =  -0.3447676561461133E+00;
  COFD[         391] =   0.1487930120557579E-01;
  COFD[         392] =  -0.1752884957692982E+02;
  COFD[         393] =   0.4691328454817521E+01;
  COFD[         394] =  -0.3953771017966438E+00;
  COFD[         395] =   0.1716337377732627E-01;
  COFD[         396] =  -0.1633601357384808E+02;
  COFD[         397] =   0.3987313912522149E+01;
  COFD[         398] =  -0.2492282800992027E+00;
  COFD[         399] =   0.8598651953545469E-02;
  COFD[         400] =  -0.2062980903668024E+02;
  COFD[         401] =   0.5158372106011019E+01;
  COFD[         402] =  -0.4299697617186866E+00;
  COFD[         403] =   0.1763220992933028E-01;
  COFD[         404] =  -0.1997245750726936E+02;
  COFD[         405] =   0.4990620273390539E+01;
  COFD[         406] =  -0.4173575849700564E+00;
  COFD[         407] =   0.1744556746387847E-01;
  COFD[         408] =  -0.1303387396358854E+02;
  COFD[         409] =   0.2940290699291584E+01;
  COFD[         410] =  -0.1728777131903053E+00;
  COFD[         411] =   0.7735387047566533E-02;
  COFD[         412] =  -0.1908503263157613E+02;
  COFD[         413] =   0.4839712718287958E+01;
  COFD[         414] =  -0.4003349576619368E+00;
  COFD[         415] =   0.1680186521574959E-01;
  COFD[         416] =  -0.1910414995755510E+02;
  COFD[         417] =   0.4841043235938669E+01;
  COFD[         418] =  -0.4007327593260240E+00;
  COFD[         419] =   0.1682891976067088E-01;
  COFD[         420] =  -0.1950109402915842E+02;
  COFD[         421] =   0.4213160796246491E+01;
  COFD[         422] =  -0.2546858655973295E+00;
  COFD[         423] =   0.8062896988504849E-02;
  COFD[         424] =  -0.2019519349819304E+02;
  COFD[         425] =   0.5057098278246105E+01;
  COFD[         426] =  -0.4237036359144079E+00;
  COFD[         427] =   0.1763107794907143E-01;
  COFD[         428] =  -0.1833537066207706E+02;
  COFD[         429] =   0.3854730356423516E+01;
  COFD[         430] =  -0.2010584934954492E+00;
  COFD[         431] =   0.5500040907376769E-02;
  COFD[         432] =  -0.1999976389851343E+02;
  COFD[         433] =   0.4996976782821175E+01;
  COFD[         434] =  -0.4180681383611913E+00;
  COFD[         435] =   0.1747184649902579E-01;
  COFD[         436] =  -0.2018972661108984E+02;
  COFD[         437] =   0.5057917524708468E+01;
  COFD[         438] =  -0.4238228321005105E+00;
  COFD[         439] =   0.1763675334011968E-01;
  COFD[         440] =  -0.2019989177195713E+02;
  COFD[         441] =   0.5056064258174870E+01;
  COFD[         442] =  -0.4235500744230125E+00;
  COFD[         443] =   0.1762363757106131E-01;
  COFD[         444] =  -0.2121396758009720E+02;
  COFD[         445] =   0.5139008041762467E+01;
  COFD[         446] =  -0.4075088009634611E+00;
  COFD[         447] =   0.1589915057303567E-01;
  COFD[         448] =  -0.1229066279603360E+02;
  COFD[         449] =   0.2741057082889793E+01;
  COFD[         450] =  -0.1457627022193994E+00;
  COFD[         451] =   0.6504615244946633E-02;
  COFD[         452] =  -0.1720445512239349E+02;
  COFD[         453] =   0.4866330968714804E+01;
  COFD[         454] =  -0.4220901924876364E+00;
  COFD[         455] =   0.1848595807162496E-01;
  COFD[         456] =  -0.1647873164232925E+02;
  COFD[         457] =   0.3810732050841141E+01;
  COFD[         458] =  -0.2812626132670425E+00;
  COFD[         459] =   0.1223411282905946E-01;
  COFD[         460] =  -0.1545394137084459E+02;
  COFD[         461] =   0.3436021618273004E+01;
  COFD[         462] =  -0.2341477712817018E+00;
  COFD[         463] =   0.1025708506984392E-01;
  COFD[         464] =  -0.1028179959593585E+02;
  COFD[         465] =   0.1896257251270889E+01;
  COFD[         466] =  -0.3329823450147470E-01;
  COFD[         467] =   0.1524623974194515E-02;
  COFD[         468] =  -0.1476396639560408E+02;
  COFD[         469] =   0.3356474034011327E+01;
  COFD[         470] =  -0.2254105976116209E+00;
  COFD[         471] =   0.9946130954817655E-02;
  COFD[         472] =  -0.1475457475365462E+02;
  COFD[         473] =   0.3343986593397501E+01;
  COFD[         474] =  -0.2237039347662408E+00;
  COFD[         475] =   0.9868653787965939E-02;
  COFD[         476] =  -0.2019519349819304E+02;
  COFD[         477] =   0.5057098278246105E+01;
  COFD[         478] =  -0.4237036359144079E+00;
  COFD[         479] =   0.1763107794907143E-01;
  COFD[         480] =  -0.1580720390088048E+02;
  COFD[         481] =   0.3572143437285478E+01;
  COFD[         482] =  -0.2518469828462103E+00;
  COFD[         483] =   0.1102533331592793E-01;
  COFD[         484] =  -0.1973436691750415E+02;
  COFD[         485] =   0.4993125184848788E+01;
  COFD[         486] =  -0.4088531920998837E+00;
  COFD[         487] =   0.1672325900844039E-01;
  COFD[         488] =  -0.1548438513232288E+02;
  COFD[         489] =   0.3444569041797159E+01;
  COFD[         490] =  -0.2352646456759844E+00;
  COFD[         491] =   0.1030578796275302E-01;
  COFD[         492] =  -0.1579979030842146E+02;
  COFD[         493] =   0.3572309323030401E+01;
  COFD[         494] =  -0.2518694694768392E+00;
  COFD[         495] =   0.1102634618303224E-01;
  COFD[         496] =  -0.1581504405580379E+02;
  COFD[         497] =   0.3572299494956542E+01;
  COFD[         498] =  -0.2518681372333372E+00;
  COFD[         499] =   0.1102628617469197E-01;
  COFD[         500] =  -0.1847775599474948E+02;
  COFD[         501] =   0.4468042284926762E+01;
  COFD[         502] =  -0.3608995423836463E+00;
  COFD[         503] =   0.1545260601935584E-01;
  COFD[         504] =  -0.1787213075631945E+02;
  COFD[         505] =   0.4905112324365878E+01;
  COFD[         506] =  -0.4173872143518051E+00;
  COFD[         507] =   0.1788845489573668E-01;
  COFD[         508] =  -0.1695183404426207E+02;
  COFD[         509] =   0.4416203740445888E+01;
  COFD[         510] =  -0.3114890192543383E+00;
  COFD[         511] =   0.1159663710619628E-01;
  COFD[         512] =  -0.2001704718429840E+02;
  COFD[         513] =   0.5030685873326798E+01;
  COFD[         514] =  -0.4055157556084015E+00;
  COFD[         515] =   0.1625731432523394E-01;
  COFD[         516] =  -0.2025975101746627E+02;
  COFD[         517] =   0.5176212790334279E+01;
  COFD[         518] =  -0.4308686680573706E+00;
  COFD[         519] =   0.1761066268522406E-01;
  COFD[         520] =  -0.1286980431258231E+02;
  COFD[         521] =   0.3009638621722941E+01;
  COFD[         522] =  -0.1808900918568948E+00;
  COFD[         523] =   0.8039983754952261E-02;
  COFD[         524] =  -0.1894546586601941E+02;
  COFD[         525] =   0.4935720676293722E+01;
  COFD[         526] =  -0.4114057447467022E+00;
  COFD[         527] =   0.1723315986098057E-01;
  COFD[         528] =  -0.1896004007077217E+02;
  COFD[         529] =   0.4935377774711531E+01;
  COFD[         530] =  -0.4113926173034602E+00;
  COFD[         531] =   0.1723402811988101E-01;
  COFD[         532] =  -0.1833537066207706E+02;
  COFD[         533] =   0.3854730356423516E+01;
  COFD[         534] =  -0.2010584934954492E+00;
  COFD[         535] =   0.5500040907376769E-02;
  COFD[         536] =  -0.1973436691750415E+02;
  COFD[         537] =   0.4993125184848788E+01;
  COFD[         538] =  -0.4088531920998837E+00;
  COFD[         539] =   0.1672325900844039E-01;
  COFD[         540] =  -0.1301206458003408E+02;
  COFD[         541] =   0.1429168452568504E+01;
  COFD[         542] =   0.1661557715508861E+00;
  COFD[         543] =  -0.1214321823404827E-01;
  COFD[         544] =  -0.2032965701164541E+02;
  COFD[         545] =   0.5194406293826699E+01;
  COFD[         546] =  -0.4326101189773075E+00;
  COFD[         547] =   0.1766176568164014E-01;
  COFD[         548] =  -0.2036133619470493E+02;
  COFD[         549] =   0.5195864695910879E+01;
  COFD[         550] =  -0.4301216528920454E+00;
  COFD[         551] =   0.1744936825492251E-01;
  COFD[         552] =  -0.1972379377029397E+02;
  COFD[         553] =   0.4985700637038792E+01;
  COFD[         554] =  -0.4077156220815392E+00;
  COFD[         555] =   0.1666668649763390E-01;
  COFD[         556] =  -0.1964903327620521E+02;
  COFD[         557] =   0.4536783758398668E+01;
  COFD[         558] =  -0.3103487526523180E+00;
  COFD[         559] =   0.1096426285225869E-01;
  COFD[         560] =  -0.1220110902734887E+02;
  COFD[         561] =   0.2686277536441436E+01;
  COFD[         562] =  -0.1383620679476532E+00;
  COFD[         563] =   0.6169543276698330E-02;
  COFD[         564] =  -0.1683623819222381E+02;
  COFD[         565] =   0.4700302341393257E+01;
  COFD[         566] =  -0.4004143455401054E+00;
  COFD[         567] =   0.1754115750370506E-01;
  COFD[         568] =  -0.1623589255135331E+02;
  COFD[         569] =   0.3719609230867837E+01;
  COFD[         570] =  -0.2699320567627367E+00;
  COFD[         571] =   0.1176389135749383E-01;
  COFD[         572] =  -0.1523800917587685E+02;
  COFD[         573] =   0.3353970687690228E+01;
  COFD[         574] =  -0.2241227404004539E+00;
  COFD[         575] =   0.9850105182191860E-02;
  COFD[         576] =  -0.1015561153728631E+02;
  COFD[         577] =   0.1822574904230457E+01;
  COFD[         578] =  -0.2281022330021956E-01;
  COFD[         579] =   0.1028887956167319E-02;
  COFD[         580] =  -0.1445980203722702E+02;
  COFD[         581] =   0.3225883257074325E+01;
  COFD[         582] =  -0.2082125398549377E+00;
  COFD[         583] =   0.9190899948945185E-02;
  COFD[         584] =  -0.1445334512020511E+02;
  COFD[         585] =   0.3215111284442823E+01;
  COFD[         586] =  -0.2067440027231415E+00;
  COFD[         587] =   0.9124392199487649E-02;
  COFD[         588] =  -0.1999976389851343E+02;
  COFD[         589] =   0.4996976782821175E+01;
  COFD[         590] =  -0.4180681383611913E+00;
  COFD[         591] =   0.1747184649902579E-01;
  COFD[         592] =  -0.1548438513232288E+02;
  COFD[         593] =   0.3444569041797159E+01;
  COFD[         594] =  -0.2352646456759844E+00;
  COFD[         595] =   0.1030578796275302E-01;
  COFD[         596] =  -0.2032965701164541E+02;
  COFD[         597] =   0.5194406293826699E+01;
  COFD[         598] =  -0.4326101189773075E+00;
  COFD[         599] =   0.1766176568164014E-01;
  COFD[         600] =  -0.1526178871149328E+02;
  COFD[         601] =   0.3359866584086094E+01;
  COFD[         602] =  -0.2248764571855301E+00;
  COFD[         603] =   0.9882257175669460E-02;
  COFD[         604] =  -0.1547465173805414E+02;
  COFD[         605] =   0.3443508682096921E+01;
  COFD[         606] =  -0.2351210220554292E+00;
  COFD[         607] =   0.1029932464766231E-01;
  COFD[         608] =  -0.1549437107750523E+02;
  COFD[         609] =   0.3445870797498957E+01;
  COFD[         610] =  -0.2354409162022220E+00;
  COFD[         611] =   0.1031371856736254E-01;
  COFD[         612] =  -0.1814306967104452E+02;
  COFD[         613] =   0.4344729601783576E+01;
  COFD[         614] =  -0.3455001737088529E+00;
  COFD[         615] =   0.1480806222798380E-01;
  COFD[         616] =  -0.1228690337278460E+02;
  COFD[         617] =   0.2739817607547901E+01;
  COFD[         618] =  -0.1455907617836413E+00;
  COFD[         619] =   0.6496697467578658E-02;
  COFD[         620] =  -0.1719486110604737E+02;
  COFD[         621] =   0.4862600572193149E+01;
  COFD[         622] =  -0.4215909171603991E+00;
  COFD[         623] =   0.1846373197001472E-01;
  COFD[         624] =  -0.1647508018746298E+02;
  COFD[         625] =   0.3812741716466225E+01;
  COFD[         626] =  -0.2815324383266347E+00;
  COFD[         627] =   0.1224615497227192E-01;
  COFD[         628] =  -0.1544409203507008E+02;
  COFD[         629] =   0.3434913447661099E+01;
  COFD[         630] =  -0.2339977102148624E+00;
  COFD[         631] =   0.1025033359000777E-01;
  COFD[         632] =  -0.1028043354360960E+02;
  COFD[         633] =   0.1896303471491212E+01;
  COFD[         634] =  -0.3330438429816188E-01;
  COFD[         635] =   0.1524890767209287E-02;
  COFD[         636] =  -0.1474363832488901E+02;
  COFD[         637] =   0.3350090901251298E+01;
  COFD[         638] =  -0.2245382243795846E+00;
  COFD[         639] =   0.9906528312288899E-02;
  COFD[         640] =  -0.1473449002110577E+02;
  COFD[         641] =   0.3337793989426398E+01;
  COFD[         642] =  -0.2228575541864579E+00;
  COFD[         643] =   0.9830229599517820E-02;
  COFD[         644] =  -0.2018972661108984E+02;
  COFD[         645] =   0.5057917524708468E+01;
  COFD[         646] =  -0.4238228321005105E+00;
  COFD[         647] =   0.1763675334011968E-01;
  COFD[         648] =  -0.1579979030842146E+02;
  COFD[         649] =   0.3572309323030401E+01;
  COFD[         650] =  -0.2518694694768392E+00;
  COFD[         651] =   0.1102634618303224E-01;
  COFD[         652] =  -0.2036133619470493E+02;
  COFD[         653] =   0.5195864695910879E+01;
  COFD[         654] =  -0.4301216528920454E+00;
  COFD[         655] =   0.1744936825492251E-01;
  COFD[         656] =  -0.1547465173805414E+02;
  COFD[         657] =   0.3443508682096921E+01;
  COFD[         658] =  -0.2351210220554292E+00;
  COFD[         659] =   0.1029932464766231E-01;
  COFD[         660] =  -0.1579169675646239E+02;
  COFD[         661] =   0.3572143437285479E+01;
  COFD[         662] =  -0.2518469828462104E+00;
  COFD[         663] =   0.1102533331592793E-01;
  COFD[         664] =  -0.1580828869487550E+02;
  COFD[         665] =   0.3572786632028933E+01;
  COFD[         666] =  -0.2519341709914386E+00;
  COFD[         667] =   0.1102926053643803E-01;
  COFD[         668] =  -0.1847553110880651E+02;
  COFD[         669] =   0.4470592201543791E+01;
  COFD[         670] =  -0.3612138681134173E+00;
  COFD[         671] =   0.1546543427860560E-01;
  COFD[         672] =  -0.1229422250798061E+02;
  COFD[         673] =   0.2742232459895579E+01;
  COFD[         674] =  -0.1459257489141948E+00;
  COFD[         675] =   0.6512123400236565E-02;
  COFD[         676] =  -0.1721362788263434E+02;
  COFD[         677] =   0.4869900420159440E+01;
  COFD[         678] =  -0.4225679381202458E+00;
  COFD[         679] =   0.1850722613298993E-01;
  COFD[         680] =  -0.1648290476852054E+02;
  COFD[         681] =   0.3809087920338710E+01;
  COFD[         682] =  -0.2810417778226483E+00;
  COFD[         683] =   0.1222425349132379E-01;
  COFD[         684] =  -0.1546403419937241E+02;
  COFD[         685] =   0.3437367417375052E+01;
  COFD[         686] =  -0.2343299630280416E+00;
  COFD[         687] =   0.1026528034463288E-01;
  COFD[         688] =  -0.1028308834746900E+02;
  COFD[         689] =   0.1896214006095980E+01;
  COFD[         690] =  -0.3329247904451228E-01;
  COFD[         691] =   0.1524374377170872E-02;
  COFD[         692] =  -0.1478376315874751E+02;
  COFD[         693] =   0.3362741024796286E+01;
  COFD[         694] =  -0.2262670743508222E+00;
  COFD[         695] =   0.9985011202484037E-02;
  COFD[         696] =  -0.1477418610290300E+02;
  COFD[         697] =   0.3350090901251285E+01;
  COFD[         698] =  -0.2245382243795827E+00;
  COFD[         699] =   0.9906528312288795E-02;
  COFD[         700] =  -0.2019989177195713E+02;
  COFD[         701] =   0.5056064258174870E+01;
  COFD[         702] =  -0.4235500744230125E+00;
  COFD[         703] =   0.1762363757106131E-01;
  COFD[         704] =  -0.1581504405580379E+02;
  COFD[         705] =   0.3572299494956542E+01;
  COFD[         706] =  -0.2518681372333372E+00;
  COFD[         707] =   0.1102628617469197E-01;
  COFD[         708] =  -0.1972379377029397E+02;
  COFD[         709] =   0.4985700637038792E+01;
  COFD[         710] =  -0.4077156220815392E+00;
  COFD[         711] =   0.1666668649763390E-01;
  COFD[         712] =  -0.1549437107750523E+02;
  COFD[         713] =   0.3445870797498957E+01;
  COFD[         714] =  -0.2354409162022220E+00;
  COFD[         715] =   0.1031371856736254E-01;
  COFD[         716] =  -0.1580828869487550E+02;
  COFD[         717] =   0.3572786632028933E+01;
  COFD[         718] =  -0.2519341709914386E+00;
  COFD[         719] =   0.1102926053643803E-01;
  COFD[         720] =  -0.1582224453447645E+02;
  COFD[         721] =   0.3572143437285498E+01;
  COFD[         722] =  -0.2518469828462133E+00;
  COFD[         723] =   0.1102533331592807E-01;
  COFD[         724] =  -0.1848021362972200E+02;
  COFD[         725] =   0.4465723773985839E+01;
  COFD[         726] =  -0.3606133849438846E+00;
  COFD[         727] =   0.1544091120986748E-01;
  COFD[         728] =  -0.1530635401605573E+02;
  COFD[         729] =   0.3873896336632167E+01;
  COFD[         730] =  -0.2956113676790514E+00;
  COFD[         731] =   0.1311329899072374E-01;
  COFD[         732] =  -0.1788417846762182E+02;
  COFD[         733] =   0.4849004827267381E+01;
  COFD[         734] =  -0.3925012395929342E+00;
  COFD[         735] =   0.1605802238455311E-01;
  COFD[         736] =  -0.1912503774058016E+02;
  COFD[         737] =   0.4674260376378655E+01;
  COFD[         738] =  -0.3838726146843580E+00;
  COFD[         739] =   0.1629606485550395E-01;
  COFD[         740] =  -0.1811180157901078E+02;
  COFD[         741] =   0.4336077093491970E+01;
  COFD[         742] =  -0.3444079746717316E+00;
  COFD[         743] =   0.1476188131006474E-01;
  COFD[         744] =  -0.1084949698568986E+02;
  COFD[         745] =   0.2037615396660156E+01;
  COFD[         746] =  -0.4955359664787870E-01;
  COFD[         747] =   0.2122039387077741E-02;
  COFD[         748] =  -0.1758398460398111E+02;
  COFD[         749] =   0.4332475641237926E+01;
  COFD[         750] =  -0.3464887767008653E+00;
  COFD[         751] =   0.1495283324430825E-01;
  COFD[         752] =  -0.1757223862273510E+02;
  COFD[         753] =   0.4319080008767561E+01;
  COFD[         754] =  -0.3447676561461133E+00;
  COFD[         755] =   0.1487930120557579E-01;
  COFD[         756] =  -0.2121396758009720E+02;
  COFD[         757] =   0.5139008041762467E+01;
  COFD[         758] =  -0.4075088009634611E+00;
  COFD[         759] =   0.1589915057303567E-01;
  COFD[         760] =  -0.1847775599474948E+02;
  COFD[         761] =   0.4468042284926762E+01;
  COFD[         762] =  -0.3608995423836463E+00;
  COFD[         763] =   0.1545260601935584E-01;
  COFD[         764] =  -0.1964903327620521E+02;
  COFD[         765] =   0.4536783758398668E+01;
  COFD[         766] =  -0.3103487526523180E+00;
  COFD[         767] =   0.1096426285225869E-01;
  COFD[         768] =  -0.1814306967104452E+02;
  COFD[         769] =   0.4344729601783576E+01;
  COFD[         770] =  -0.3455001737088529E+00;
  COFD[         771] =   0.1480806222798380E-01;
  COFD[         772] =  -0.1847553110880651E+02;
  COFD[         773] =   0.4470592201543791E+01;
  COFD[         774] =  -0.3612138681134173E+00;
  COFD[         775] =   0.1546543427860560E-01;
  COFD[         776] =  -0.1848021362972200E+02;
  COFD[         777] =   0.4465723773985839E+01;
  COFD[         778] =  -0.3606133849438846E+00;
  COFD[         779] =   0.1544091120986748E-01;
  COFD[         780] =  -0.2069000320825284E+02;
  COFD[         781] =   0.5101158792846945E+01;
  COFD[         782] =  -0.4264320198306029E+00;
  COFD[         783] =   0.1763100216755543E-01;
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
  KTDIF[           1] =            2;
  KTDIF[           2] =            5;
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
  COFTD[           4] =  -0.1525347878079855E+00;
  COFTD[           5] =  -0.5464040665780899E-04;
  COFTD[           6] =   0.2934125089337968E-07;
  COFTD[           7] =  -0.4870919731196246E-11;
  COFTD[           8] =   0.4225303552740193E+00;
  COFTD[           9] =   0.1320842793755543E-03;
  COFTD[          10] =  -0.7122224283164592E-07;
  COFTD[          11] =   0.1195161063668316E-10;
  COFTD[          12] =   0.4452620885099525E+00;
  COFTD[          13] =   0.4946972054970773E-04;
  COFTD[          14] =  -0.2630235132337235E-07;
  COFTD[          15] =   0.4903063330400376E-11;
  COFTD[          16] =   0.1616137301763599E+00;
  COFTD[          17] =   0.4741552334083435E-04;
  COFTD[          18] =  -0.1671152373050399E-07;
  COFTD[          19] =  -0.1889821705016664E-11;
  COFTD[          20] =   0.4155834523029710E+00;
  COFTD[          21] =   0.1097383911871691E-04;
  COFTD[          22] =  -0.3960224834087406E-08;
  COFTD[          23] =   0.1144145307552153E-11;
  COFTD[          24] =   0.4219325599835782E+00;
  COFTD[          25] =   0.1114149277732201E-04;
  COFTD[          26] =  -0.4020727469049531E-08;
  COFTD[          27] =   0.1161625074178193E-11;
  COFTD[          28] =   0.1613923904161360E+00;
  COFTD[          29] =   0.5010841669717118E-03;
  COFTD[          30] =  -0.2382739816291534E-06;
  COFTD[          31] =   0.3493444507237758E-10;
  COFTD[          32] =   0.4444526949057779E+00;
  COFTD[          33] =   0.7145255629999482E-04;
  COFTD[          34] =  -0.3862572696699677E-07;
  COFTD[          35] =   0.6889797705021213E-11;
  COFTD[          36] =   0.6020321200392946E-01;
  COFTD[          37] =   0.5615615962686471E-03;
  COFTD[          38] =  -0.2553727790421687E-06;
  COFTD[          39] =   0.3633898304582181E-10;
  COFTD[          40] =   0.4446537425501394E+00;
  COFTD[          41] =   0.5066317405316880E-04;
  COFTD[          42] =  -0.2698209736758623E-07;
  COFTD[          43] =   0.5012898785410712E-11;
  COFTD[          44] =   0.4427392097103151E+00;
  COFTD[          45] =   0.7117708739452678E-04;
  COFTD[          46] =  -0.3847681435586731E-07;
  COFTD[          47] =   0.6863235673780089E-11;
  COFTD[          48] =   0.4460703094919803E+00;
  COFTD[          49] =   0.7171261254133862E-04;
  COFTD[          50] =  -0.3876630782084383E-07;
  COFTD[          51] =   0.6914873573367535E-11;
  COFTD[          52] =   0.3257425690224397E+00;
  COFTD[          53] =   0.3036334352955890E-03;
  COFTD[          54] =  -0.1552903446031957E-06;
  COFTD[          55] =   0.2414664576246333E-10;
  COFTD[          56] =   0.1525347878079855E+00;
  COFTD[          57] =   0.5464040665780899E-04;
  COFTD[          58] =  -0.2934125089337968E-07;
  COFTD[          59] =   0.4870919731196246E-11;
  COFTD[          60] =   0.0000000000000000E+00;
  COFTD[          61] =   0.0000000000000000E+00;
  COFTD[          62] =   0.0000000000000000E+00;
  COFTD[          63] =   0.0000000000000000E+00;
  COFTD[          64] =   0.1654293333730542E+00;
  COFTD[          65] =   0.5612389631132157E-03;
  COFTD[          66] =  -0.2656505647884942E-06;
  COFTD[          67] =   0.3882296208775095E-10;
  COFTD[          68] =   0.2407445350007590E+00;
  COFTD[          69] =   0.4453434831756282E-03;
  COFTD[          70] =  -0.2181738909148746E-06;
  COFTD[          71] =   0.3269585309435789E-10;
  COFTD[          72] =   0.3407625358054880E+00;
  COFTD[          73] =  -0.4040578770176091E-04;
  COFTD[          74] =   0.3278795395318293E-07;
  COFTD[          75] =  -0.6270938141265901E-11;
  COFTD[          76] =   0.2700102625474920E+00;
  COFTD[          77] =   0.3615551214404962E-03;
  COFTD[          78] =  -0.1807447676772846E-06;
  COFTD[          79] =   0.2753212716087208E-10;
  COFTD[          80] =   0.2720417770528547E+00;
  COFTD[          81] =   0.3642754049836650E-03;
  COFTD[          82] =  -0.1821046627192010E-06;
  COFTD[          83] =   0.2773927453061686E-10;
  COFTD[          84] =  -0.1246478985813586E+00;
  COFTD[          85] =   0.7965256655680803E-03;
  COFTD[          86] =  -0.3249988027840275E-06;
  COFTD[          87] =   0.4325167626647755E-10;
  COFTD[          88] =   0.2209079675460805E+00;
  COFTD[          89] =   0.4810899053474402E-03;
  COFTD[          90] =  -0.2333769631025201E-06;
  COFTD[          91] =   0.3471383309607501E-10;
  COFTD[          92] =  -0.1418836567013430E+00;
  COFTD[          93] =   0.7665588601033144E-03;
  COFTD[          94] =  -0.3065500230258140E-06;
  COFTD[          95] =   0.4029595271283201E-10;
  COFTD[          96] =   0.2394100540729324E+00;
  COFTD[          97] =   0.4471972131282707E-03;
  COFTD[          98] =  -0.2189517205098697E-06;
  COFTD[          99] =   0.3279735363108610E-10;
  COFTD[         100] =   0.2204829568680314E+00;
  COFTD[         101] =   0.4801643237618346E-03;
  COFTD[         102] =  -0.2329279629943256E-06;
  COFTD[         103] =   0.3464704623498651E-10;
  COFTD[         104] =   0.2213085142120510E+00;
  COFTD[         105] =   0.4819622095914184E-03;
  COFTD[         106] =  -0.2338001183446032E-06;
  COFTD[         107] =   0.3477677564298334E-10;
  COFTD[         108] =   0.2443704234105951E-01;
  COFTD[         109] =   0.7182425471481816E-03;
  COFTD[         110] =  -0.3197185256879899E-06;
  COFTD[         111] =   0.4488287145198055E-10;
  COFTD[         112] =  -0.1616137301763599E+00;
  COFTD[         113] =  -0.4741552334083435E-04;
  COFTD[         114] =   0.1671152373050399E-07;
  COFTD[         115] =   0.1889821705016664E-11;
  COFTD[         116] =  -0.3407625358054880E+00;
  COFTD[         117] =   0.4040578770176091E-04;
  COFTD[         118] =  -0.3278795395318293E-07;
  COFTD[         119] =   0.6270938141265901E-11;
  COFTD[         120] =   0.4663152819246087E+00;
  COFTD[         121] =  -0.5601505800988478E-04;
  COFTD[         122] =   0.4659876633602610E-07;
  COFTD[         123] =  -0.9136462914858121E-11;
  COFTD[         124] =   0.4220100375895677E+00;
  COFTD[         125] =  -0.4140425208641108E-04;
  COFTD[         126] =   0.4387516156683957E-07;
  COFTD[         127] =  -0.1028602437920772E-10;
  COFTD[         128] =   0.0000000000000000E+00;
  COFTD[         129] =   0.0000000000000000E+00;
  COFTD[         130] =   0.0000000000000000E+00;
  COFTD[         131] =   0.0000000000000000E+00;
  COFTD[         132] =   0.3315880233940311E+00;
  COFTD[         133] =  -0.1963882557865052E-04;
  COFTD[         134] =   0.3023888393001825E-07;
  COFTD[         135] =  -0.8449980090305466E-11;
  COFTD[         136] =   0.3422032136897760E+00;
  COFTD[         137] =  -0.2026752702741960E-04;
  COFTD[         138] =   0.3120692705763941E-07;
  COFTD[         139] =  -0.8720490905912508E-11;
  COFTD[         140] =   0.3566291429282130E+00;
  COFTD[         141] =   0.1061168823365192E-03;
  COFTD[         142] =  -0.5729173426194156E-07;
  COFTD[         143] =   0.9652647160627432E-11;
  COFTD[         144] =   0.4436492507591140E+00;
  COFTD[         145] =  -0.4874846422125629E-04;
  COFTD[         146] =   0.4697186596249824E-07;
  COFTD[         147] =  -0.1035687579663234E-10;
  COFTD[         148] =   0.2849835858178697E+00;
  COFTD[         149] =   0.1154600116315077E-03;
  COFTD[         150] =  -0.6171979338326076E-07;
  COFTD[         151] =   0.1015042212845429E-10;
  COFTD[         152] =   0.4221715234173643E+00;
  COFTD[         153] =  -0.4177501100679457E-04;
  COFTD[         154] =   0.4397262274959381E-07;
  COFTD[         155] =  -0.1026729305000021E-10;
  COFTD[         156] =   0.4402209436618440E+00;
  COFTD[         157] =  -0.4837175963856001E-04;
  COFTD[         158] =   0.4660889007292850E-07;
  COFTD[         159] =  -0.1027684286354751E-10;
  COFTD[         160] =   0.4468957655172876E+00;
  COFTD[         161] =  -0.4910519334513472E-04;
  COFTD[         162] =   0.4731559438265286E-07;
  COFTD[         163] =  -0.1043266483507838E-10;
  COFTD[         164] =   0.4596633923456365E+00;
  COFTD[         165] =  -0.1747709148267824E-04;
  COFTD[         166] =   0.1428880756642898E-07;
  COFTD[         167] =  -0.2036106279600706E-11;
};




#if 0




\\
\\
\\  This is the mechanism file
\\
\\
!***************************************************************************
!****************************************************************************
!          Revised H2/CO high temperature Combustion Mechanism
!
!              Scott Davis, Ameya Joshi, and Hai Wang
!Department of Mechanical Engineering, University of Delaware, Neark, DE 19716
!
!                                 January 2003
!
!****************************************************************************
!
! Reference sources can be found at the end of the file.
!
!****************************************************************************
!
! Please contact Scott Davis at daviss@exponent.com or
! Hai Wang at hwang@me.udel.edu for questions and comments
!     
!============================================================================  
ELEMENTS
O  H  C  N  AR HE
END
SPECIES
!
!First Two Species Must be H2 and H for Modified Burner Stabilized Flame Calc.
!
H2  H
!
AR  N2  HE
!
! radicals in order of increasing molecular weight:
!
O  OH  HCO HO2
!
! stable molecules in order of increasing molecular weight:
!
H2O  CO  O2  H2O2  CO2
!
REACTIONS
!
! Reactions of H2/O2
!
 H+O2 = O+OH             2.644E+16  -0.6707 17041.00 !GRI3.0 * 1.00
 O+H2 = H+OH             4.589E+04   2.700   6260.00 !GRI3.0 * 1.19
 OH+H2 = H+H2O           1.734E+08   1.510   3430.00 !GRI3.0 * 0.80
 OH+OH = O+H2O           3.973E+04   2.400  -2110.00 !GRI3.0 * 1.11
 H+H+M = H2+M            1.780E+18  -1.000      0.00 !GRI3.0 * 1.78
   H2/0.0/ H2O/0.0/ CO2/0.0/ AR/0.63/ HE/0.63/
 H+H+H2 = H2+H2          9.000E+16  -0.600      0.00 !GRI3.0
 H+H+H2O = H2+H2O        5.624E+19  -1.250      0.00 !GRI3.0 * 0.94
 H+H+CO2 = H2+CO2        5.500E+20  -2.000      0.00 !GRI3.0
 H+OH+M = H2O+M          4.400E+22  -2.000      0.00 !GRI3.0 * 2.00
   H2/2.0/ H2O/6.30/ CO/1.75/ CO2/3.6/  AR/0.38/ HE/0.38/
 O+H+M = OH+M            9.428E+18  -1.000      0.00 !86TSA/HAM * 2.00
   H2/2.0/ H2O/12.0/ CO/1.75/ CO2/3.6/  AR/0.7/ HE/0.7/
 O+O+M = O2+M            1.200E+17  -1.000      0.00 !GRI3.0
   H2/2.4/ H2O/15.4/  CO/1.75/ CO2/3.6/  AR/0.83/ HE/0.83/
 H+O2(+M) = HO2(+M)      5.116E+12   0.440      0.00 !00TROE - Based on M=N2 * 1.10
   LOW/6.328E+19  -1.400  0.00/
   TROE/0.5  1E-30  1E+30/
   O2/0.85/  H2O/11.89/ CO/1.09/ CO2/2.18/ AR/0.40/ HE/0.46/ H2/0.75/
!  O2/0.75/  H2O/12.0/ CO/1.2/ CO2/2.4/ AR/0.53/ HE/0.53/
! H+O2(+M) = HO2(+M)    4.651E+12   0.440      0.00 !00TROE - Based on M=AR
!   LOW/7.490E+18  -1.200  0.00/
!   TROE/0.5  1E-30  1E+30/
! H+O2(+M) = HO2(+M)  4.651E+12   0.440      0.00 !00TROE - Based on M=H2O
!   LOW/5.753E+20  -1.400  0.00/                      !10xN2
!   TROE/0.0 345.0 10 345.0 /                         !FSC
 H2+O2 = HO2+H           5.916E+05   2.433  53502.00 !00MIC/SUT * 0.80
 OH+OH(+M) = H2O2(+M)    1.110E+14   -0.370     0.00 !88ZEL/EWI * 1.50
   LOW  /  2.010E+17   -0.584  -2293.00/     !Fit 88ZEL/EWI and 92BAU/COB
   TROE/  0.7346   94.00  1756.00  5182.00 / !H2O=6xN2 88ZEL/EWI
   H2/2.0/ H2O/6.00/ CO/1.75/ CO2/3.6/ AR/0.7/ HE/0.7/
!
! Reactions of HO2
!
 HO2+H = O+H2O                 3.970E+12    0.000    671.00 !GRI3.0
 HO2+H = OH+OH                 7.485E+13    0.000    295.00 !99MUE/KIM * 1.06
 HO2+O = OH+O2                 4.000E+13    0.000      0.00 !GRI3.0 * 2.00
 HO2+OH = O2+H2O               2.375E+13    0.000   -500.00 !88KEY * 0.82
  DUPLICATE
 HO2+OH = O2+H2O               1.000E+16    0.000  17330.00 !95HIP/NEU
  DUPLICATE
 HO2+HO2 = O2+H2O2             1.300E+11    0.000  -1630.00 !90HIP/TRO
  DUPLICATE
 HO2+HO2 = O2+H2O2             3.658E+14    0.000  12000.00 !90HIP/TRO * 0.87
  DUPLICATE
!
! Reactions of H2O2
!
 H2O2+H = HO2+H2               6.050E+06    2.000   5200.00 !GRI3.0 * 0.50
 H2O2+H = OH+H2O               2.410E+13    0.000   3970.00 !86TSA/HAM 
 H2O2+O = OH+HO2               9.630E+06    2.000   3970.00 !86TSA/HAM
 H2O2+OH = HO2+H2O             2.000E+12    0.000    427.00 !95HIP/NEU
  DUPLICATE
 H2O2+OH = HO2+H2O             2.670E+41   -7.000  37600.00 !Refit95HIP/NEU
  DUPLICATE                                                 !2.2E14 MAX K
!
! Reactions of CO/CO2
!
 CO+O(+M)=CO2(+M)              1.362E+10    0.000   2384.00 !99MUE/KIM * 0.76
   LOW/1.173E+24 -2.79  4191./
   H2/2.0/ H2O/12/ CO/1.75/ CO2/3.6/ AR/0.7/ HE/0.7/
! CO+OH = CO2+H                 4.760E+07    1.228     70.00 !GRI3.0
  CO+OH=CO2+H                  8.000E+11    0.140   7352.00 !This Work * 0.83
DUP
  CO+OH=CO2+H                  8.784E+10   0.030    -16.00 !          * 1.20
DUP
 CO+O2 = CO2+O                 1.119E+12    0.000  47700.00 !86TSA/HAM * 0.44
 CO+HO2 = CO2+OH               3.010E+13    0.000  23000.00 !99MUE/KIM
!
! Reactions of HCO
!
 HCO+H = CO+H2                 1.200E+14    0.000      0.00 !02FRI/DAV * 1.00
 HCO+O = CO+OH                 3.000E+13    0.000      0.00 !GRI3.0 
 HCO+O = CO2+H                 3.000E+13    0.000      0.00 !GRI3.0
 HCO+OH = CO+H2O               3.020E+13    0.000      0.00 !86TSA/HAM
 HCO+M = CO+H+M                1.870E+17   -1.000  17000.00 !02FRI/DAV * 2.00
    H2/2.0/ H2O/0.0/ CO/1.75/ CO2/3.6/ 
 HCO+H2O = CO+H+H2O            2.244E+18   -1.000  17000.00 !12xM * 2.00
 HCO+O2 = CO+HO2               1.204E+10    0.807   -727.00 !96HSU/MEB
END

\\
\\
\\  This is the therm file
\\
\\
THERMO
   298.000  1000.000  5000.000
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
OH                S 9/01O   1H   1    0    0G   200.000  6000.000   1000.00    1
 2.86472886E+00 1.05650448E-03-2.59082758E-07 3.05218674E-11-1.33195876E-15    2
 3.71885774E+03 5.70164073E+00 4.12530561E+00-3.22544939E-03 6.52764691E-06    3
-5.79853643E-09 2.06237379E-12 3.38153812E+03-6.90432960E-01 4.51532273E+03    4
<<<<<<< HEAD
=======
OH-OLD            RUS 78O   1H   1   00   00G   200.000  3500.000  1000.000    1
 3.09288767E+00 5.48429716E-04 1.26505228E-07-8.79461556E-11 1.17412376E-14    2
 3.85865700E+03 4.47669610E+00 3.99201543E+00-2.40131752E-03 4.61793841E-06    3
-3.88113333E-09 1.36411470E-12 3.61508056E+03-1.03925458E-01 8.81310600E+03    4
>>>>>>> 96047fcebbc048e76be07737f734966d562c3001
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
CH2*              L S/93C   1H   2   00   00G   200.000  3500.000  1000.000    1
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
C2O               121286C   2O   1          G  0300.00   5000.00  1000.00      1
 0.04849809E+02 0.02947585E-01-0.01090729E-04 0.01792562E-08-0.01115758E-12    2
 0.03282055E+06-0.06453226E+01 0.03368851E+02 0.08241803E-01-0.08765145E-04    3
 0.05569262E-07-0.01540009E-10 0.03317081E+06 0.06713314E+02                   4
C2H               L 1/91C   2H   1   00   00G   200.000  3500.000  1000.000    1
 3.16780652E+00 4.75221902E-03-1.83787077E-06 3.04190252E-10-1.77232770E-14    2
 6.71210650E+04 6.63589475E+00 2.88965733E+00 1.34099611E-02-2.84769501E-05    3
 2.94791045E-08-1.09331511E-11 6.68393932E+04 6.22296438E+00 1.04544720E+04    4
C2H2              L 1/91C   2H   2   00   00G   200.000  3500.000  1000.000    1
 4.14756964E+00 5.96166664E-03-2.37294852E-06 4.67412171E-10-3.61235213E-14    2
 2.59359992E+04-1.23028121E+00 8.08681094E-01 2.33615629E-02-3.55171815E-05    3
 2.80152437E-08-8.50072974E-12 2.64289807E+04 1.39397051E+01 1.00058390E+04    4
H2CC              L12/89H   2C   2    0    0G   200.000  6000.000  1000.000    1
 0.42780340E+01 0.47562804E-02-0.16301009E-05 0.25462806E-09-0.14886379E-13    2
 0.48316688E+05 0.64023701E+00 0.32815483E+01 0.69764791E-02-0.23855244E-05    3
-0.12104432E-08 0.98189545E-12 0.48621794E+05 0.59203910E+01 0.49887266E+05    4
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
CH2CHO            T04/83O   1H   3C   2    0G   300.000  5000.000              1
 0.59756699E+01 0.81305914E-02-0.27436245E-05 0.40703041E-09-0.21760171E-13    2
 0.49032178E+03-0.50320879E+01 0.34090624E+01 0.10738574E-01 0.18914925E-05    3
-0.71585831E-08 0.28673851E-11 0.15214766E+04 0.95714535E+01 0.30474436E+04    4
CH3CO             T 9/92C   2H   3O   1    0G   200.000  6000.0    1000.0      1
 0.59447731E+01 0.78667205E-02-0.28865882E-05 0.47270875E-09-0.28599861E-13    2
-0.37873075E+04-0.50136751E+01 0.41634257E+01-0.23261610E-03 0.34267820E-04    3
-0.44105227E-07 0.17275612E-10-0.26574529E+04 0.73468280E+01-0.12027167E+04    4
CH3CHO            L 8/88C   2H   4O   1    0G   200.000  6000.0    1000.0      1
 0.54041108E+01 0.11723059E-01-0.42263137E-05 0.68372451E-09-0.40984863E-13    2
-0.22593122E+05-0.34807917E+01 0.47294595E+01-0.31932858E-02 0.47534921E-04    3
-0.57458611E-07 0.21931112E-10-0.21572878E+05 0.41030159E+01-0.19987949E+05    4
HCCO              SRIC91H   1C   2O   1     G  0300.00   4000.00  1000.00      1
 0.56282058E+01 0.40853401E-02-0.15934547E-05 0.28626052E-09-0.19407832E-13    2
 0.19327215E+05-0.39302595E+01 0.22517214E+01 0.17655021E-01-0.23729101E-04    3
 0.17275759E-07-0.50664811E-11 0.20059449E+05 0.12490417E+02                   4
HCCOH              SRI91C   2O   1H   20   0G   300.000  5000.000   1000.0     1
 0.59238291E+01 0.67923600E-02-0.25658564E-05 0.44987841E-09-0.29940101E-13    2
 0.72646260E+04-0.76017742E+01 0.12423733E+01 0.31072201E-01-0.50866864E-04    3
 0.43137131E-07-0.14014594E-10 0.80316143E+04 0.13874319E+02                   4
CH2OCH            T 9/92C   2H   3O   1    0G   298.150  3000.000  1000.0      1
 0.48131470E+00 0.20711914E-01-0.12693155E-04 0.34579642E-08-0.35399703E-12    2
 0.15648642E+05 0.34629876E+02 0.10854772E+01 0.12845259E-01 0.24138660E-05    3
-0.44642672E-08-0.29381916E-12 0.15910655E+05 0.33395312E+02 0.16817588E+05    4
CH2OCH2           T 6/92C   2H   4O   1    0G   298.150  3000.0    1000.0      1
 0.54887641E+01 0.12046190E-01-0.43336931E-05 0.70028311E-09-0.41949088E-13    2
-0.91804251E+04-0.70799605E+01 0.37590532E+01-0.94412180E-02 0.80309721E-04    3
-0.10080788E-06 0.40039921E-10-0.75608143E+04 0.78497475E+01-0.63304657E+04    4
C3H2              UT8/99C   3H   2          G   200.000  3500.00 1000.00       1
 0.71572642E+01 0.55660834E-02-0.22944369E-05 0.44922385E-09-0.34017182E-13    2
 0.66001539E+05-0.11347860E+02 0.29415061E+01 0.27682947E-01-0.46382396E-04    3
 0.39372626E-07-0.12770546E-10 0.66671312E+05 0.79103537E+01                   4
cC3H2             121686C   3H   2          G  0300.00   5000.00  1000.00      1
 0.06530853E+02 0.05870316E-01-0.01720777E-04 0.02127498E-08-0.08291910E-13    2
 0.05115214E+06-0.01122728E+03 0.02691077E+02 0.01480366E+00-0.03250551E-04    3
-0.08644363E-07 0.05284878E-10 0.05219072E+06 0.08757391E+02                   4
C3H3              T 5/97C   3H   3    0    0G   200.000  6000.000              1
 7.14221880E+00 7.61902005E-03-2.67459950E-06 4.24914801E-10-2.51475415E-14    2
 3.89087427E+04-1.25848436E+01 1.35110927E+00 3.27411223E-02-4.73827135E-05    3
 3.76309808E-08-1.18540923E-11 4.01057783E+04 1.52058924E+01 4.16139977E+04    4
aC3H4             L 8/89C   3H   4    0    0G   200.000  6000.000              1
 0.63168722E+01 0.11133728E-01-0.39629378E-05 0.63564238E-09-0.37875540E-13    2
 0.20117495E+05-0.10995766E+02 0.26130445E+01 0.12122575E-01 0.18539880E-04    3
-0.34525149E-07 0.15335079E-10 0.21541567E+05 0.10226139E+02 0.22962267E+05    4
pC3H4             T 2/90H   4C   3    0    0G   200.000  6000.000              1
 0.60252400E+01 0.11336542E-01-0.40223391E-05 0.64376063E-09-0.38299635E-13    2
 0.19620942E+05-0.86043785E+01 0.26803869E+01 0.15799651E-01 0.25070596E-05    3
-0.13657623E-07 0.66154285E-11 0.20802374E+05 0.98769351E+01 0.22302059E+05    4
cC3H4             T12/81C   3H   4    0    0G   300.000  5000.000              1
 0.66999931E+01 0.10357372E-01-0.34551167E-05 0.50652949E-09-0.26682276E-13    2
 0.30199051E+05-0.13378770E+02-0.24621047E-01 0.23197215E-01-0.18474357E-05    3
-0.15927593E-07 0.86846155E-11 0.32334137E+05 0.22729762E+02 0.3332728 E+05    4
aC3H5             PD5/98C   3H   5    0    0G   300.000  3000.000              1
 0.65007877E+01 0.14324731E-01-0.56781632E-05 0.11080801E-08-0.90363887E-13    2
 0.17482449E+05-0.11243050E+02 0.13631835E+01 0.19813821E-01 0.12497060E-04    3
-0.33355555E-07 0.15846571E-10 0.19245629E+05 0.17173214E+02                   4
CH3CCH2           PD5/98C   3H   5    0    0G   300.000  3000.000              1
 0.54255528E+01 0.15511072E-01-0.56678350E-05 0.79224388E-09-0.16878034E-13    2
 0.27843027E+05-0.33527184E+01 0.17329209E+01 0.22394620E-01-0.51490611E-05    3
-0.67596466E-08 0.38253211E-11 0.29040498E+05 0.16568878E+02                   4
CH3CHCH           PD5/98C   3H   5    0    0G   300.000  3000.000              1
 0.53725281E+01 0.15780509E-01-0.59922850E-05 0.93089664E-09-0.36550966E-13    2
 0.29614760E+05-0.34186478E+01 0.91372931E+00 0.26432343E-01-0.11758950E-04    3
-0.23035678E-08 0.27715488E-11 0.30916867E+05 0.19989269E+02                   4
C3H6              120186C   3H   6          G  0300.00   5000.00  1000.00      1
 0.06732257E+02 0.01490834E+00-0.04949899E-04 0.07212022E-08-0.03766204E-12    2
-0.09235703E+04-0.01331335E+03 0.01493307E+02 0.02092518E+00 0.04486794E-04    3
-0.01668912E-06 0.07158146E-10 0.01074826E+05 0.01614534E+03                   4
nC3H7             P11/94C   3H   7    0    0G   300.000  3000.000              1
 0.77097479E+01 0.16031485E-01-0.52720238E-05 0.75888352E-09-0.38862719E-13    2
 0.79762236E+04-0.15515297E+02 0.10491173E+01 0.26008973E-01 0.23542516E-05    3
-0.19595132E-07 0.93720207E-11 0.10312346E+05 0.21136034E+02                   4
iC3H7             P11/94C   3H   7    0    0G   300.000  3000.000              1
 0.65192741E+01 0.17220104E-01-0.57364217E-05 0.84130732E-09-0.44565913E-13    2
 0.73227193E+04-0.90830215E+01 0.14449199E+01 0.20999112E-01 0.77036222E-05    3
-0.18476253E-07 0.71282962E-11 0.94223724E+04 0.20116317E+02                   4
C3H8              P11/94C   3H   8    0    0G   300.000  3000.000              1
 0.75244152E+01 0.18898282E-01-0.62921041E-05 0.92161457E-09-0.48684478E-13    2
-0.16564394E+05-0.17838375E+02 0.92851093E+00 0.26460566E-01 0.60332446E-05    3
-0.21914953E-07 0.94961544E-11-0.14057907E+05 0.19225538E+02                   4
CH2CHCO           T05/99C   3H   3O   1    0G   200.000  6000.0    1000.0      1
 6.95842227E+00 1.07193211E-02-3.85218494E-06 6.22009064E-10-3.72401640E-14    2
 5.64826498E+03-1.14745786E+01 3.21169467E+00 1.18422105E-02 1.67462582E-05    3
-3.06947176E-08 1.33048816E-11 7.12815750E+03 1.00881663E+01 8.70564832E+03    4
C2H3CHO           T 6/92C   3H   4O   1    0G   298.150  3000.0    1000.0      1
 0.48353180E+01 0.19772601E-01-0.10426628E-04 0.26525803E-08-0.26278207E-12    2
-0.11557837E+05 0.18853144E+01 0.11529584E+01 0.28040214E-01-0.15072153E-04    3
 0.15905842E-08 0.84930371E-12-0.10417694E+05 0.21453279E+02-0.89572567E+04    4
CH3COCH3          T 5/92C   3H   6O   1    0G   200.000  6000.000  1000.0      1
 0.72975991E+01 0.17566207E-01-0.63170456E-05 0.10203086E-08-0.61094016E-13    2
-0.29817680E+05-0.12756981E+02 0.55557943E+01-0.28365428E-02 0.70568945E-04    3
-0.87810488E-07 0.34028266E-10-0.28113337E+05 0.23226600E+01-0.26116945E+05    4
CH3CHOCH2         T 6/92C   3H   6O   1    0G   298.150  3000.0    1000.0      1
 0.86900558E+01 0.16020987E-01-0.53971753E-05 0.79941542E-09-0.42656366E-13    2
-0.15420691E+05-0.22485016E+02 0.48733836E+00 0.28519690E-01 0.30096162E-05    3
-0.22652642E-07 0.10706728E-10-0.12556434E+05 0.22605270E+02-0.11156446E+05    4
CH3CH2CHO         T 9/92C   3H   6O   1    0G   273.150  5000.000  1000.0      1
 0.33137982E+01 0.26619606E-01-0.10475596E-04 0.18815334E-08-0.12761310E-12    2
-0.25459603E+05 0.96608447E+01 0.76044596E+01-0.86403564E-02 0.73930097E-04    3
-0.79687398E-07 0.28004927E-10-0.25489789E+05-0.67643691E+01-0.23097645E+05    4
C4H               P 1/93C   4H   1    0    0G   300.000  3000.000              1
 0.77697593E+01 0.49829976E-02-0.17628546E-05 0.28144284E-09-0.16689869E-13    2
 0.94345900E+05-0.14165274E+02 0.13186295E+01 0.38582956E-01-0.71385623E-04    3
 0.65356359E-07-0.22617666E-10 0.95456106E+05 0.15567583E+02                   4
C4H2              P 1/93C   4H   2    0    0G   300.000  3000.000              1
 0.86637708E+01 0.67247189E-02-0.23593397E-05 0.37506380E-09-0.22230940E-13    2
 0.53252275E+05-0.21093503E+02-0.39201030E+00 0.51937565E-01-0.91737340E-04    3
 0.80471986E-07-0.26898218E-10 0.54845266E+05 0.20957794E+02                   4
C4H4              H6W/94C   4H   4    0    0G   300.000  3000.000              1
 0.66507092E+01 0.16129434E-01-0.71938875E-05 0.14981787E-08-0.11864110E-12    2
 0.31195992E+05-0.97952118E+01-0.19152479E+01 0.52750878E-01-0.71655944E-04    3
 0.55072423E-07-0.17286228E-10 0.32978504E+05 0.31419983E+02                   4
n-C4H3            H6W/94C   4H   3    0    0G   300.000  3000.000              1
 0.54328279E+01 0.16860981E-01-0.94313109E-05 0.25703895E-08-0.27456309E-12    2
 0.61600680E+05-0.15673981E+01-0.31684113E+00 0.46912100E-01-0.68093810E-04    3
 0.53179921E-07-0.16523005E-10 0.62476199E+05 0.24622559E+02                   4
i-C4H3            AB1/93C   4H   3    0    0G   300.000  3000.000              1
 0.90978165E+01 0.92207119E-02-0.33878441E-05 0.49160498E-09-0.14529780E-13    2
 0.56600574E+05-0.19802597E+02 0.20830412E+01 0.40834274E-01-0.62159685E-04    3
 0.51679358E-07-0.17029184E-10 0.58005129E+05 0.13617462E+02                   4
n-C4H5            H6W/94C   4H   5    0    0G   300.000  3000.000              1
 0.98501978E+01 0.10779008E-01-0.13672125E-05-0.77200535E-09 0.18366314E-12    2
 0.38840301E+05-0.26001846E+02 0.16305321E+00 0.39830137E-01-0.34000128E-04    3
 0.15147233E-07-0.24665825E-11 0.41429766E+05 0.23536163E+02                   4
i-C4H5            H6W/94C   4H   5    0    0G   300.000  3000.000              1
 0.10229092E+02 0.94850138E-02-0.90406445E-07-0.12596100E-08 0.24781468E-12    2
 0.34642812E+05-0.28564529E+02-0.19932900E-01 0.38005672E-01-0.27559450E-04    3
 0.77835551E-08 0.40209383E-12 0.37496223E+05 0.24394241E+02                   4
C4H5-2            H6W/94C   4H   5    0    0G   300.000  3000.000              1
 1.45381710E+01-8.56770560E-03 2.35595240E-05-1.36763790E-08 2.44369270E-12    2
 3.32590950E+04-4.53694970E+01 2.96962800E+00 2.44422450E-02-9.12514240E-06    3
-4.24668710E-18 1.63047280E-21 3.55033160E+04 1.20360510E+01 3.73930550E+04    4
C4H6              H6W/94C   4H   6    0    0G   300.000  3000.000              1
 0.88673134E+01 0.14918670E-01-0.31548716E-05-0.41841330E-09 0.15761258E-12    2
 0.91338516E+04-0.23328171E+02 0.11284465E+00 0.34369022E-01-0.11107392E-04    3
-0.92106660E-08 0.62065179E-11 0.11802270E+05 0.23089996E+02                   4
C4H612            A 8/83C   4H   6    0    0G   300.     3000.     1000.0      1
  0.1781557E+02 -0.4257502E-02  0.1051185E-04 -0.4473844E-08  0.5848138E-12    2
  0.1267342E+05 -0.6982662E+02  0.1023467E+01  0.3495919E-01 -0.2200905E-04    3
  0.6942272E-08 -0.7879187E-12  0.1811799E+05  0.1975066E+02  0.1950807E+05    4
C4H6-2            A 8/83C   4H   6    0    0G   300.     3000.     1000.0      1
  9.0338133E+00  8.2124510E-03  7.1753952E-06 -5.8834334E-09  1.0343915E-12    2
  1.4335068E+04 -2.0985762E+01  2.1373338E+00  2.6486229E-02 -9.0568711E-06    3
 -5.5386397E-19  2.1281884E-22  1.5710902E+04  1.3529426E+01  1.7488676E+04    4
n-C4H7            H6W/94C   4H   7    0    0G   300.000  3000.000              1
 0.11963392E+02 0.11425305E-01 0.78948909E-06-0.19858872E-08 0.36873645E-12    2
 0.16962977E+05-0.37542908E+02 0.28698254E+00 0.36964495E-01-0.86277441E-05    3
-0.15051821E-07 0.89891263E-11 0.20551301E+05 0.24484467E+02                   4
C4H7              AM1/94C   4H   7    0    0G   300.000  3000.000  1000.000    1
 0.11963392E+02 0.11425305E-01 0.78948909E-06-0.19858872E-08 0.36873645E-12    2
 0.16962977E+05-0.37542908E+02 0.28698254E+00 0.36964495E-01-0.86277441E-05    3
-0.15051821E-07 0.89891263E-11 0.20551301E+05 0.24484467E+02                   4
iC4H7             P11/94C   4H   7    0    0G   300.000  3000.000              1
 0.74491956E+01 0.22630504E-01-0.88095014E-05 0.14336478E-08-0.73247269E-13    2
 0.11196182E+05-0.11947779E+02 0.34512660E+01 0.24686039E-01 0.52359514E-05    3
-0.16130826E-07 0.53881687E-11 0.12783361E+05 0.11080150E+02                   4
C4H81             T 6/83C   4H   8    0    0G   300.000  5000.000              1
 0.20535841E+01 0.34350507E-01-0.15883197E-04 0.33089662E-08-0.25361045E-12    2
-0.21397231E+04 0.15543201E+02 0.11811380E+01 0.30853380E-01 0.50865247E-05    3
-0.24654888E-07 0.11110193E-10-0.17904004E+04 0.21062469E+02                   4
C4H82             T 6/83C   4H   8    0    0G   300.000  5000.00               1
 0.82797676E+00 0.35864539E-01-0.16634498E-04 0.34732759E-08-0.26657398E-12    2
-0.30521033E+04 0.21342545E+02 0.12594252E+01 0.27808424E-01 0.87013932E-05    3
-0.24402205E-07 0.98977710E-11-0.29647742E+04 0.20501129E+02                   4
iC4H8             T 6/83H   8C   4    0    0G   300.000  5000.0                1
 0.44609470E+01 0.29611487E-01-0.13077129E-04 0.26571934E-08-0.20134713E-12    2
-0.50066758E+04 0.10671549E+01 0.26471405E+01 0.25902957E-01 0.81985354E-05    3
-0.22193259E-07 0.88958580E-11-0.40373069E+04 0.12676388E+02                   4
pC4H9             P11/94C   4H   9    0    0G   300.000  3000.000              1
 0.90131759E+01 0.23992501E-01-0.89488794E-05 0.15311024E-08-0.98216680E-13    2
 0.47050163E+04-0.20568160E+02 0.22128506E+01 0.32808150E-01 0.35445091E-05    3
-0.24586397E-07 0.11510570E-10 0.71406715E+04 0.17148115E+02                   4
sC4H9             P11/94C   4H   9    0    0G   300.000  3000.000              1
 0.94263839E+01 0.21918998E-01-0.72868375E-05 0.10630334E-08-0.55649464E-13    2
 0.31965874E+04-0.22406051E+02 0.69428423E+00 0.33113346E-01 0.62942577E-05    3
-0.27025274E-07 0.11989315E-10 0.64175654E+04 0.26279789E+02                   4
sC4H9b            T07/95C   4H   9    0    0G   200.000  6000.000  1000.0      1
 0.88057265E+01 0.23630381E-01-0.84564737E-05 0.13612584E-08-0.81313232E-13    2
 0.37941169E+04-0.19996770E+02 0.46457042E+01 0.79313214E-02 0.70027013E-04    3
-0.95973349E-07 0.38628890E-10 0.62341181E+04 0.79642808E+01 0.84190169E+04    4
tC4H9             P11/94C   4H   9    0    0G   300.000  3000.000              1
 0.76607261E+01 0.23879414E-01-0.80890353E-05 0.12057521E-08-0.65009814E-13    2
 0.16207623E+04-0.14800281E+02 0.96167553E+00 0.25735856E-01 0.15609033E-04    3
-0.26656519E-07 0.89418010E-11 0.46564412E+04 0.24805366E+02                   4
iC4H9             P11/94C   4H   9    0    0G   300.000  3000.000              1
 0.36943491E+01 0.36043526E-01-0.18406555E-04 0.46356137E-08-0.45862968E-12    2
 0.53753792E+04 0.68486993E+01-0.52798452E+00 0.50164629E-01-0.37679302E-04    3
 0.17505938E-07-0.39549774E-11 0.64430984E+04 0.28114795E+02                   4
C4H10             P11/94C   4H  10    0    0G   300.000  3000.000              1
 0.10526774E+02 0.23590738E-01-0.78522480E-05 0.11448408E-08-0.59827703E-13    2
-0.20479223E+05-0.32198579E+02 0.15685419E+01 0.34652278E-01 0.68168129E-05    3
-0.27995097E-07 0.12307742E-10-0.17129977E+05 0.17908045E+02                   4
iC4H10            P11/94C   4H  10    0    0G   300.000  3000.000              1
 0.10846169E+02 0.23338389E-01-0.77833962E-05 0.11393807E-08-0.59918289E-13    2
-0.21669854E+05-0.35870573E+02 0.54109489E+00 0.37860301E-01 0.55459804E-05    3
-0.30500110E-07 0.14033357E-10-0.17977644E+05 0.21150935E+02                   4
H2C4O             120189H   2C   4O   1     G  0300.00   4000.00  1000.00      1
 0.01026888E+03 0.04896164E-01-0.04885081E-05-0.02708566E-08 0.05107013E-12    2
 0.02346903E+06-0.02815985E+03 0.04810971E+02 0.01313999E+00 0.09865073E-05    3
-0.06120720E-07 0.01640003E-10 0.02545803E+06 0.02113424E+02                   4
C4H4O             T03/97C   4H   4O   1    0G   200.000  6000.0    1000.0      1
 9.38935003E+00 1.40291241E-02-5.07755110E-06 8.24137332E-10-4.95319963E-14    2
-8.68241814E+03-2.79162920E+01 8.47469463E-01 1.31773796E-02 5.99735901E-05    3
-9.71562904E-08 4.22733796E-11-5.36785445E+03 2.14945172E+01-4.17166616E+03    4
C4H6O25           T 3/97C   4H   6O   1    0G   200.000  5000.000  1000.0      1
 8.60658242E+00 2.08310051E-02-8.42229481E-06 1.56717640E-09-1.09391202E-13    2
-1.76177415E+04-2.32464750E+01 2.67053463E+00 4.92586420E-03 8.86967406E-05    3
-1.26219194E-07 5.23991321E-11-1.46572472E+04 1.45722395E+01-1.30831522E+04    4
C4H6O23           T 3/97C   4H   6O   1    0G   200.000  5000.000  1000.0      1
 8.60658242E+00 2.08310051E-02-8.42229481E-06 1.56717640E-09-1.09391202E-13    2
-1.32392815E+04-2.32464750E+01 2.67053463E+00 4.92586420E-03 8.86967406E-05    3
-1.26219194E-07 5.23991321E-11-1.02787872E+04 1.45722395E+01-1.30831522E+04    4
C2H3CHOCH2        A 8/83C   4H   6O   1    0G   300.     3000.     1000.0      1
-4.72093360E+00 3.91413780E-02-6.52872650E-06-7.68209500E-09 2.51473310E-12    2
 1.75352252E+03 5.17190420E+01 7.97985440E-01 3.44034320E-02-1.24598510E-05    3
-5.18062790E-18 1.99359540E-21-6.48927540E+02 2.18896980E+01 1.00654250E+03    4
CH3CHCHCHO        T 5/92C   4H   6O   1    0G   298.150  3000.0    1000.0      1
 1.98794540E+01-2.09130550E-02 4.45360508E-05-2.60374870E-08 4.86836120E-12    2
-1.95278768E+04-6.87200320E+01-1.55577660E+00 4.09640630E-02-1.69868810E-05    3
-6.00928140E-18 2.31368530E-21-1.41394920E+04 3.74707580E+01-1.29340710E+04    4
CH3CHCHCO         T03/97C   4H   5O   1    0G   200.000  6000.000  1000.0      1
 8.90967920E+00 1.34364140E-02-7.62977390E-07-1.69114810E-09 2.95540440E-13    2
 1.48898740E+03-1.79662460E+01-1.08199860E+00 3.64929760E-02-1.52255950E-05    3
-5.62607170E-18 2.16113750E-21 3.56713230E+03 3.27142550E+01 4.73074990E+03    4
CH2CHCHCHO        T03/97C   4H   5O   1    0G   200.000  6000.000  1000.0      1
 8.90967920E+00 1.34364140E-02-7.62977390E-07-1.69114810E-09 2.95540440E-13    2
 1.48898740E+03-1.79662460E+01-1.08199860E+00 3.64929760E-02-1.52255950E-05    3
-5.62607170E-18 2.16113750E-21 3.56713230E+03 3.27142550E+01 4.73074990E+03    4
C5H2               20587C   5H   2          G  0300.00   5000.00  1000.00      1
 0.01132917E+03 0.07424057E-01-0.02628189E-04 0.04082541E-08-0.02301333E-12    2
 0.07878706E+06-0.03617117E+03 0.03062322E+02 0.02709998E+00-0.01009170E-03    3
-0.01272745E-06 0.09167219E-10 0.08114969E+06 0.07071078E+02                   4
C5H3               20387C   5H   3          G  0300.00   5000.00  1000.00      1
 0.01078762E+03 0.09539619E-01-0.03206745E-04 0.04733323E-08-0.02512135E-12    2
 0.06392904E+06-0.03005444E+03 0.04328720E+02 0.02352480E+00-0.05856723E-04    3
-0.01215449E-06 0.07726478E-10 0.06588531E+06 0.04173259E+02                   4
C5H5              L 7/89C   5H   5    0    0G   200.000  6000.000              1
 0.10844066e+02 0.15392837e-01-0.55630421e-05 0.90189371e-09-0.54156531e-13    2
 0.26900566e+05-0.35254948e+02-0.95903718e+00 0.31396859e-01 0.26723794e-04    3
-0.68941872e-07 0.33301856e-10 0.30729120e+05 0.29072816e+02 0.31954258e+05    4
C5H6            C-P10/85C   5H   6    0    0G   298.150  5000.000              1
 0.10624320E+02 0.17735448E-01-0.62330446E-05 0.97308317E-09-0.55500130E-13    2
 0.10772188E+05-0.35773422E+02-0.28978958E+01 0.43484777E-01-0.33511005E-05    3
-0.31103756E-07 0.16912444E-10 0.15084742E+05 0.36894760E+02 0.16068486E+05    4
C5H4OH            L 8/89C   5H   5O   1    0G   200.000  6000.000              1
 0.13367912e+02 0.15205785e-01-0.54592258e-05 0.88134866e-09-0.52774454e-13    2
 0.38411506e+04-0.45920839e+02-0.12822236e+01 0.49041160e-01-0.13688997e-04    3
-0.29133858e-07 0.19006964e-10 0.80087098e+04 0.30798358e+02 0.96365992e+04    4
C5H4O             P 1/93C   5H   4O   1    0G   300.000  3000.000              1
 0.47927242E+01 0.29221680E-01-0.15997486E-04 0.42069049E-08-0.42815179E-12    2
 0.22849286E+04-0.30131893E+01-0.23915355E+01 0.47363680E-01-0.30728171E-04    3
 0.78031552E-08-0.25145729E-12 0.43740152E+04 0.34594337E+02                   4
C5H5O             L 7/89C   5O   1H   5    0G   200.000  6000.000              1
 0.12606422e+02 0.16747260e-01-0.61098574e-05 0.99676557e-09-0.60113201e-13    2
 0.39313455e+04-0.42604277e+02 0.23042835e+00 0.32322691e-01 0.28900443e-04    3
-0.70679977e-07 0.33406891e-10 0.80753082e+04 0.25330974e+02                   4
cC5H8             T03/97C   5H   8O   0    0G   200.000  6000.000  1000.0      1
 0.77244792E+01 0.28322316E-01-0.11545236E-04 0.21540815E-08-0.15054178E-12    2
-0.78261573E+03-0.19769698E+02 0.26898140E+01 0.20954550E-02 0.11303687E-03    3
-0.15408070E-06 0.62763658E-10 0.23139663E+04 0.15294056E+02 0.39328836E+04    4
lC5H9             T03/97C   5H   9O   0    0G   200.000  6000.000  1000.0      1
 0.20313000E+02 0.10869880E-01-0.19063805E-05 0.00000000E+00 0.00000000E+00    2
 0.94061603E+04-0.82533815E+02 0.11430827E+01 0.44350789E-01-0.17825470E-04    3
 0.00000000E+00 0.00000000E+00 0.16967656E+05 0.24181940E+02 0.19122233E+05    4
cC5H9             T03/97C   5H   9O   0    0G   200.000  6000.000  1000.0      1
 0.11406802E+02 0.22563988E-01-0.70235595E-05 0.11321968E-08-0.73438204E-13    2
 0.75268769E+04-0.39636280E+02 0.29427128E+00 0.13823374E-01 0.90847653E-04    3
-0.13008694E-06 0.53051811E-10 0.12565712E+05 0.27389773E+02 0.13838458E+05    4
C6H2              P 1/93C   6H   2    0    0G   300.000  3000.000              1
 0.13226281E+02 0.73904302E-02-0.22715381E-05 0.25875217E-09-0.55356741E-14    2
 0.80565258E+05-0.41201176E+02-0.15932624E+01 0.80530145E-01-0.14800649E-03    3
 0.13300031E-06-0.45332313E-10 0.83273227E+05 0.27980873E+02                   4
C6H               P 1/93C   6H   1    0    0G   300.000  3000.000              1
 0.12370055E+02 0.52177699E-02-0.16885009E-05 0.25807149E-09-0.15472851E-13    2
 0.12158739E+06-0.34952797E+02-0.25630299E+00 0.63793827E-01-0.11440118E-03    3
 0.10136744E-06-0.34361855E-10 0.12408855E+06 0.24930750E+02                   4
l-C6H4            H6W/94C   6H   4    0    0G   300.000  3000.000              1
 0.12715182E+02 0.13839662E-01-0.43765440E-05 0.31541636E-09 0.46619026E-13    2
 0.57031148E+05-0.39464600E+02 0.29590225E+00 0.58053318E-01-0.67766756E-04    3
 0.43376762E-07-0.11418864E-10 0.60001371E+05 0.22318970E+02                   4
c-C6H4            H6W/94C   6H   4    0    0G   300.000  3000.000              1
 0.13849209E+02 0.78807920E-02 0.18243836E-05-0.21169166E-08 0.37459977E-12    2
 0.47446340E+05-0.50404953E+02-0.30991268E+01 0.54030564E-01-0.40839004E-04    3
 0.10738837E-07 0.98078490E-12 0.52205711E+05 0.37415207E+02                   4
n-C6H5            H6W/94C   6H   5    0    0G   300.000  3000.000              1
 0.16070068E+02 0.81899539E-02 0.17325165E-05-0.20624185E-08 0.36292345E-12    2
 0.64616867E+05-0.56163742E+02-0.61135769E+00 0.65082610E-01-0.78262397E-04    3
 0.53030828E-07-0.14946683E-10 0.68805375E+05 0.27635468E+02                   4
l-C6H6            H6W/94C   6H   6    0    0G   300.000  3000.000              1
 0.17584442E+02 0.64486600E-02 0.48933980E-05-0.34696221E-08 0.56150749E-12    2
 0.34111988E+05-0.66017838E+02-0.10170622E+01 0.61794821E-01-0.59461061E-04    3
 0.31873491E-07-0.71717693E-11 0.39202707E+05 0.29460373E+02                   4
c-C6H7            H6W/94C   6H   7    0    0G   300.000  3000.000              1
 0.19996841E+02 0.11189543E-02 0.11649756E-04-0.62779471E-08 0.94939508E-12    2
 0.16730059E+05-0.83746933E+02-0.30328493E+01 0.50804518E-01-0.69150292E-05    3
-0.29715974E-07 0.16296353E-10 0.23895383E+05 0.38909180E+02                   4
n-C6H7            H6W/94C   6H   7    0    0G   300.000  3000.000              1
 0.22577469E+02-0.30737517E-02 0.14225234E-04-0.69880848E-08 0.10232874E-11    2
 0.41228980E+05-0.91568619E+02 0.13248032E+00 0.57103366E-01-0.43712644E-04    3
 0.15538603E-07-0.12976356E-11 0.47730512E+05 0.25339081E+02                   4
C6H8              H6W/94C   6H   8    0    0G   300.000  3000.000              1
 0.28481979E+02-0.15702948E-01 0.26771697E-04-0.11780109E-07 0.16573427E-11    2
 0.93346445E+04-0.12500226E+03 0.15850439E+01 0.40215142E-01 0.78439543E-05    3
-0.38761325E-07 0.18545207E-10 0.17949613E+05 0.19112625E+02                   4
C6H3              H6W/94C   6H   3    0    0G   300.000  3000.000              1
 0.58188343E+01 0.27933408E-01-0.17825427E-04 0.53702536E-08-0.61707627E-12    2
 0.85188250E+05-0.92147827E+00 0.11790619E+01 0.55547360E-01-0.73076168E-04    3
 0.52076736E-07-0.15046964E-10 0.85647312E+05 0.19179199E+02                   4
i-C6H5            H6W/94C   6H   5    0    0G   300.000  3000.000              1
 0.22501663E+02-0.81009977E-02 0.15955695E-04-0.72310371E-08 0.10310424E-11    2
 0.58473410E+05-0.91224777E+02-0.78585434E+00 0.60221825E-01-0.62890264E-04    3
 0.36310730E-07-0.87000259E-11 0.64942270E+05 0.28658905E+02                   4
i-C6H7            H6W/94C   6H   7    0    0G   300.000  3000.000              1
 0.20481506E+02 0.79439697E-03 0.11450761E-04-0.60991177E-08 0.91756724E-12    2
 0.37728426E+05-0.81812073E+02-0.17099094E+01 0.62486034E-01-0.54290707E-04    3
 0.26959682E-07-0.58999090E-11 0.44086621E+05 0.33344772E+02                   4
A1                H6W/94C   6H   6    0    0G   300.000  3000.000              1
 0.17246994E+02 0.38420164E-02 0.82776232E-05-0.48961120E-08 0.76064545E-12    2
 0.26646055E+04-0.71945175E+02-0.48998680E+01 0.59806932E-01-0.36710087E-04    3
 0.32740399E-08 0.37600886E-11 0.91824570E+04 0.44095642E+02                   4
A1-               H6W/94C   6H   5    0    0G   300.000  3000.000              1
 0.14493439E+02 0.75712688E-02 0.37894542E-05-0.30769500E-08 0.51347820E-12    2
 0.33189977E+05-0.54288940E+02-0.49076147E+01 0.59790771E-01-0.45639827E-04    3
 0.14964993E-07-0.91767826E-12 0.38733410E+05 0.46567780E+02                   4
C6H5O             L12/84C   6H   5O   1    0G   300.000  5000.000              1
 0.13833984E+02 0.17618403E-01-0.60696257E-05 0.91988173E-09-0.50449181E-13    2
-0.69212549E+03-0.50392990E+02-0.18219433E+01 0.48122510E-01-0.46792302E-05    3
-0.34018594E-07 0.18649637E-10 0.42429180E+04 0.33526199E+02 0.57367379E+04    4
C6H5OH            L 4/84C   6H   6O   1    0G   300.000  5000.000              1
 0.14912073E+02 0.18378135E-01-0.61983128E-05 0.91983221E-09-0.49209565E-13    2
-0.18375199E+05-0.55924103E+02-0.16956539E+01 0.52271299E-01-0.72024050E-05    3
-0.35859603E-07 0.20449073E-10-0.13284121E+05 0.32542160E+02-0.11594207E+05    4
lC6H9             T 2/92C   6H   9O    0   0G   200.000  3000.000 1000.0       1
 0.23165919E+02 0.10813608E-01-0.17638168E-05 0.00000000E+00 0.00000000E+00    2
 0.11162402E+05-0.98600332E+02 0.31671271E+00 0.52069818E-01-0.21965057E-04    3
 0.00000000E+00 0.00000000E+00 0.19926824E+05 0.27879902E+02 0.22141533E+05    4
cC6H9             T 2/92C   6H   9O    0   0G   200.000  3000.000 1000.0       1
 0.26295828E+02 0.86828857E-02-0.15770376E-05 0.00000000E+00 0.00000000E+00    2
 0.20863563E+04-0.12573825E+03-0.35714300E+01 0.61696043E-01-0.26928803E-04    3
 0.00000000E+00 0.00000000E+00 0.13657039E+05 0.39986250E+02 0.15096500E+05    4
cC6H8             T03/97C   6H   8O   0    0G   200.000  6000.000  1000.0      1
 0.11779870E+02 0.25519980E-01-0.92666947E-05 0.15068122E-08-0.90658701E-13    2
 0.65486686E+04-0.41618805E+02 0.17265319E+01 0.14887612E-01 0.94809230E-04    3
-0.14083394E-06 0.58859873E-10 0.11021297E+05 0.19130886E+02 0.12784878E+05    4
A1CH3             L 6/87C   7H   8    0    0G   300.000  3000.000              1
 0.12940034E+02 0.26691287E-01-0.96838505E-05 0.15738629E-08-0.94663601E-13    2
-0.69764908E+03-0.46728785E+02 0.16152663E+01 0.21099438E-01 0.85366018E-04    3
-0.13261066E-06 0.55956604E-10 0.40756300E+04 0.20282210E+02 0.60135835E+04    4
A1CH2             L 6/87C   7H   7    0    0G   300.000  3000.000              1
 0.14043980E+02 0.23493873E-01-0.85375367E-05 0.13890841E-08-0.83614420E-13    2
 0.18564203E+05-0.51665589E+02 0.48111540E+00 0.38512832E-01 0.32861492E-04    3
-0.76972721E-07 0.35423068E-10 0.23307027E+05 0.23548820E+02 0.25317186E+05    4
A1C2H3            HW /94C   8H   8    0    0G   300.000  3000.000              1
 0.11303213E+02 0.33709887E-01-0.13208885E-04 0.21140962E-08-0.87311377E-13    2
 0.11725388E+05-0.34737919E+02-0.38678493E+01 0.67947865E-01-0.25230333E-04    3
-0.18017145E-07 0.12998470E-10 0.16200269E+05 0.45271770E+02                   4
cC6H10            T03/97C   6H  10O   0    0G   200.000  6000.000  1000.0      1
 0.11773904E+02 0.30947360E-01-0.11234330E-04 0.18262494E-08-0.10985119E-12    2
-0.72028376E+04-0.42658688E+02 0.23662378E+01 0.10681712E-01 0.11822112E-03    3
-0.16567854E-06 0.67612802E-10-0.24824973E+04 0.16769357E+02-0.55324968E+03    4
A1C2H             HW /94C   8H   6    0    0G   300.000  3000.000  1000.0      1
 0.24090759E+02 0.78232400E-03 0.11453964E-04-0.61620504E-08 0.93346685E-12    2
 0.27429445E+05-0.10499631E+03-0.52645016E+01 0.84511042E-01-0.76597848E-04    3
 0.33216978E-07-0.47673063E-11 0.35566242E+05 0.46378815E+02                   4
A2                HW /94C  10H   8    0    0G   300.000  3000.000  1000.0      1
 0.36468643E+02-0.15419513E-01 0.30160038E-04-0.13700120E-07 0.19582730E-11    2
 0.35091445E+04-0.17329489E+03-0.94505043E+01 0.11137849E+00-0.10345667E-03    3
 0.52800392E-07-0.11804439E-10 0.16695594E+05 0.65187668E+02                   4
N2                121286N   2               G  0300.00   5000.00  1000.00      1
 0.02926640E+02 0.14879768E-02-0.05684760E-05 0.10097038E-09-0.06753351E-13    2
-0.09227977E+04 0.05980528E+02 0.03298677E+02 0.14082404E-02-0.03963222E-04    3
 0.05641515E-07-0.02444854E-10-0.10208999E+04 0.03950372E+02                   4
AR                120186AR  1               G  0300.00   5000.00  1000.00      1
 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-0.07453750E+04 0.04366000E+02 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-0.07453750E+04 0.04366000E+02                   4
HE                L10/90HE  1    0    0    0G   200.000  6000.000 1000.        1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-7.45375000E+02 9.28723974E-01 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-7.45375000E+02 9.28723974E-01 0.00000000E+00    4
NE                L10/92NE  1    0    0    0G   200.000  6000.000   1000.00    1
 0.25000000E+01 0.             0.             0.             0.                2
-0.74537500E+03 0.33553227E+01 0.25000000E+01 0.             0.                3
 0.             0.            -0.74537498E+03 0.33553227E+01 0.00000000E+00    4
CH2CHCOCH3        T 3/97C   4H   6O   1    0G   200.000  3000.0    1000.0      1
 1.98794540E+01-2.09130550E-02 4.45360580E-05-2.60374870E-08 4.86836120E-12    2
-1.90786168E+04-6.97265750E+01-1.55577660E+00 4.09640630E-02-1.69868810E-05    3
-6.00928140E-18 2.31368530E-21-1.49447258E+04 3.64642160E+01-1.66079520E+04    4
CH2CHCH2CHO       T 5/92C   4H   6O   1    0G   298.150  3000.0    1000.0      1
 1.98794540E+01-2.09130550E-02 4.45360508E-05-2.60374870E-08 4.86836120E-12    2
-1.58539966E+04-6.71095639E+01-1.55577660E+00 4.09640630E-02-1.69868810E-05    3
-6.00928140E-18 2.31368530E-21-1.04656118E+04 3.90812260E+01-1.29340710E+04    4
ENDOFDATA

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
aC3H5              2   266.800     4.982     0.000     0.000     1.000
CH3CCH2            2   266.800     4.982     0.000     0.000     1.000
CH3CHCH            2   266.800     4.982     0.000     0.000     1.000
C3H6               2   266.800     4.982     0.000     0.000     1.000
iC3H7              2   266.800     4.982     0.000     0.000     1.000
nC3H7              2   266.800     4.982     0.000     0.000     1.000
C3H8               2   266.800     4.982     0.000     0.000     1.000
CH2OCH2            1   252.000     4.760     0.000     0.000     1.000
CH2OCH             1   252.000     4.760     0.000     0.000     1.000
CH3CH2CHO          1   252.000     4.760     0.000     0.000     1.000
                                                                                
C4H                1   357.000     5.180     0.000     0.000     1.000          
C4H2               1   357.000     5.180     0.000     0.000     1.000          
H2C4O              2   357.000     5.180     0.000     0.000     1.000 ! JAM    
C4H2OH             2   224.700     4.162     0.000     0.000     1.000 ! *      
i-C4H3             2   357.000     5.180     0.000     0.000     1.000 ! JAM    
n-C4H3             2   357.000     5.180     0.000     0.000     1.000 ! JAM    
C4H4               2   357.000     5.180     0.000     0.000     1.000 ! JAM    
i-C4H5             2   357.000     5.180     0.000     0.000     1.000 ! JAM    
n-C4H5             2   357.000     5.180     0.000     0.000     1.000 ! JAM    
C4H5-2             2   357.000     5.180     0.000     0.000     1.000 !
C4H6               2   357.000     5.180     0.000     0.000     1.000         
C4H6-2             2   357.000     5.180     0.000     0.000     1.000
C4H612             2   357.000     5.180     0.000     0.000     1.000 
C4H7               2   357.000     5.180     0.000     0.000     1.000
iC4H7              2   357.000     5.180     0.000     0.000     1.000
C4H81              2   357.000     5.180     0.000     0.000     1.000
C4H82              2   357.000     5.180     0.000     0.000     1.000
iC4H8              2   357.000     5.180     0.000     0.000     1.000
tC4H9              2   357.000     5.180     0.000     0.000     1.000
iC4H9              2   357.000     5.180     0.000     0.000     1.000
pC4H9              2   357.000     5.180     0.000     0.000     1.000
sC4H9              2   357.000     5.180     0.000     0.000     1.000
C4H10              2   357.000     5.180     0.000     0.000     1.000
iC4H10             2   357.000     5.180     0.000     0.000     1.000
CH3CHOCH2          2   357.000     5.180     0.000     0.000     1.000
                                                                                
C5H2               1   357.000     5.180     0.000     0.000     1.000          
C5H3               1   357.000     5.180     0.000     0.000     1.000          
C5H5               1   357.000     5.180     0.000     0.000     1.000          
C5H6               1   357.000     5.180     0.000     0.000     1.000          
C4H6O25            1   357.000     5.180     0.000     0.000     1.000
C4H6O23            1   357.000     5.180     0.000     0.000     1.000
C4H4O              1   357.000     5.180     0.000     0.000     1.000
CH2CHCO            1   357.000     5.180     0.000     0.000     1.000
CH3CHOCH2          1   357.000     5.180     0.000     0.000     1.000
CH2CHCHCHO         1   357.000     5.180     0.000     0.000     1.000
CH3CHCHCO          1   357.000     5.180     0.000     0.000     1.000
C2H3CHOCH2         1   357.000     5.180     0.000     0.000     1.000
CH3CHCHCHO         1   357.000     5.180     0.000     0.000     1.000
                                                                                
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
                                                                                
                                                                                
!The Lennard-Jones parameters of polycyclic aromatic hydrocarbons were estimated 
!based on the critical temperature and pressure. See H. Wang and M. Frenklach,   
!"Transport Properties of Polycyclic Aromatic Hydrocarbons for Flame Modeling."  
!Combustion and Flame, 96:163-170 (1994)                                         
!                                                                                
!c-C6H4             2   464.8       5.29      0.00     10.32      0.000  !  benze
!A1                 2   464.8       5.29      0.00     10.32      0.000  !  benze
!A1-                2   464.8       5.29      0.00     10.32      0.000  !  benze
!A1CH3              2   495.3       5.68      0.43     12.30      1.000  !
!A1C2H3             2   546.2       6.00      0.13     15.00      1.000  !
!A1CH2              2   495.3       5.68      0.43     12.30      1.000  !
!A1C2H              2   535.6       5.72      0.77     12.00      1.000  !
!A2                 2   630.4       6.18      0.00     16.50      1.000  !
!c-C6H7             2   464.8       5.29      0.00     10.32      0.000  !  benze
!C5H4O              2   464.8       5.29      0.00     10.32      0.000  !  benze
!C5H5O              2   464.8       5.29      0.00     10.32      0.000  !  benze
!C5H4OH             2   464.8       5.29      0.00     10.32      0.000  !  benze
!C6H5O              2   464.8       5.29      0.00     10.32      0.000  !  benze
!C6H5OH             2   464.8       5.29      0.00     10.32      0.000  !  benze
!
!CH3COCH3           2   357.000     5.176     0.000     0.000     1.000
!C2H3CHO            2   357.000     5.176     0.000     0.000     1.000
!iC4H7O             2   450.000     5.500     0.000     0.000     1.000 ! JAM
!CH3CHO             2   436.000     3.970     0.000     0.000     2.000
!CH3CO              2   436.000     3.970     0.000     0.000     2.000
! 1-15: Species name                                                             
! 16-80: Molecular parameters                                                    
!        molecule index: 0 = atom, 1= linear molec.                              
!                        2 = nonlinear molec.                                    
!        L-J potential well depth, e/kb (K)                                      
!        L-J collision diameter, s,                                              
!        Dipole moment, f, Debye                                                 
!        Polarizability, `,  2                                                   
!        Rotational relaxation number, Zrot at 298K                              
!        Comments                                                                
!END
!H HE   -9.66994265100   2.10026266000   -0.07705964500   0.00546112600 ! Middha et al, Proc. Comb. Inst., Vol. 29
!        0.93003284006   0.08015000695   -0.00947327267   0.00063458775
!        0.87637862374   0.10238278295   -0.01480299828   0.00098803605
!        1.06001553391  -0.05992751365    0.01026504945  -0.00073450868
!H H2  -11.74984983000   3.15068443400   -0.25747189600   0.01589155500 ! Middha et al, Proc. Comb. Inst., Vol. 29
!        0.68564849197   0.15339038119   -0.01367350360   0.00032210353
!        0.67794933764   0.13747070299   -0.01105868963   0.00016576525
!        0.65119387990   0.08215979427   -0.00438944167  -0.00027358319
!H2 HE -12.75127347000   3.42444798700   -0.28472577300   0.01593170100 ! Middha et al, Proc. Comb. Inst., Vol. 29
!        0.59534394712   0.20780905298   -0.02484767627   0.00100127755
!        0.67714740207   0.13579177260   -0.01142533631   0.00018885880
!        0.65191818266   0.08365135783   -0.00445010155  -0.00025986744
!H AR   -9.05107284400   1.61614185700   -0.00287779500   0.00130541500 ! AIChE 2002
!        0.68819287418   0.15342306998   -0.01769949486   0.00088795976
!        0.69683113086   0.17157988788   -0.02435683182   0.00136547873
!        0.67026727852   0.11433417790   -0.01504070012   0.00051875848
!N H2  -11.06296595000   2.35003553100   -0.10371499000   0.00580309100 ! Stallcop et al, Phys. Rev. A, 64, Art. 042722
!        1.31576376016  -0.13458345098    0.02492533153  -0.00119098283
!        6.99070003974  -2.59649971924    0.37720200602  -0.01788555212
!        1.28549476326  -0.20526851249    0.03796449261  -0.00230740895
!N N2  -14.50976666000   3.27038987700   -0.22411274000   0.01070436600 ! Stallcop et al, Phys. Rev. A, 64, Art. 042722
!        1.28549476326  -0.20526851249    0.03796449261  -0.00230740895
!        1.31470940230  -0.12581573177    0.02093874922  -0.00091239288
!        1.18530004024  -0.00109904512   -0.00679031552   0.00082436382
!H O2  -11.04103178000   2.40427694900   -0.10279690200   0.00532644300 ! Stallcop et al, Phys. Rev. A, 64, Art. 042722
!        1.29254489763  -0.18498887975    0.04271031056  -0.00260823132
!        1.99269990199  -0.27220069165    0.02059506743   0.00016444098
!        1.27566027457  -0.22786150854    0.04449949355  -0.00274619496
!O O2  -14.60250025000   3.29049804400   -0.22351565500   0.01068641100 ! Stallcop et al, Phys. Rev. A, 64, Art. 042722
!        1.29827044958  -0.17304622424    0.03635310801  -0.00206726443
!        9.96160341198  -3.95288550918    0.58275322682  -0.02828930035
!        1.26594984746  -0.24767207966    0.04946972829  -0.00300929402
!H N2  -13.27028844000   3.51865269300   -0.29664901800   0.01643138100 ! Stallcop et al, J. Chem. Phys., 97, 3431 (1992)
!        1.33864596568  -0.08545398502    0.00922905086   0.00004406488 
!       -2.20299987672   1.59160087079   -0.25339378410   0.01364477086
!        1.27615648812  -0.22593374251    0.04641828877  -0.00310151225
!H2 N2 -10.99943193000   2.20257995900   -0.08115516500   0.00440608700 ! Stallcop et al, Phys. Rev. A, 62, Art. 062709
!        1.31648435027  -0.13202102692    0.02416204045  -0.00120258929
!        1.94229781522  -0.41910587564    0.06711915359  -0.00336850295
!        1.29072759648  -0.19301215992    0.03398993578  -0.00196835629
!N2 N2 -16.51750614000   4.05271572500   -0.34593622800   0.01671006600 ! Stallcop et al, Phys. Rev. A, 62, Art. 062709
!        1.35301607939  -0.05405097875   -0.00113356164   0.00059799157
!        3.61489786826  -1.08630677454    0.15481851859  -0.00723438967
!        1.25620046061  -0.26881554497    0.05577645068  -0.00340833729
!H2 H2  -9.96095484000   2.05602189500   -0.06497689600   0.00413678100 ! Stallcop et al, J. Thermophys. Heat Tra., 12, 514 (1998)
!        1.32208755845  -0.12074925804    0.02204710917  -0.00105955055
!        3.63140155962  -1.13979323643    0.16942078618  -0.00812678015
!        1.29811062932  -0.17814060954    0.02956707397  -0.00170200851
!ENDDIFF

#endif
