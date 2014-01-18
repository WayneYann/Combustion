
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
static const double imw[9] = {
    1.0 / 2.015940,  /*H2 */
    1.0 / 1.007970,  /*H */
    1.0 / 15.999400,  /*O */
    1.0 / 31.998800,  /*O2 */
    1.0 / 17.007370,  /*OH */
    1.0 / 18.015340,  /*H2O */
    1.0 / 33.006770,  /*HO2 */
    1.0 / 34.014740,  /*H2O2 */
    1.0 / 28.013400};  /*N2 */

struct ReactionData {
    double fwd_A,fwd_beta,fwd_Ea;
    double low_A,low_beta,low_Ea;
    double rev_A,rev_beta,rev_Ea;
    double troe_a,troe_Ts, troe_Tss, troe_Tsss;
    double sri_a, sri_b, sri_c, sri_d, sri_e;
    double activation_units, prefactor_units, phase_units;
    int is_PD, troe_len, sri_len;
};

static struct ReactionData R[27], R_DEF[27];

static double fwd_A[27], fwd_beta[27], fwd_Ea[27];
static double low_A[27], low_beta[27], low_Ea[27];
static double rev_A[27], rev_beta[27], rev_Ea[27];
static double troe_a[27],troe_Ts[27], troe_Tss[27], troe_Tsss[27];
static double sri_a[27], sri_b[27], sri_c[27], sri_d[27], sri_e[27];
static double activation_units[27], prefactor_units[27], phase_units[27];
static int is_PD[27], troe_len[27], sri_len[27];
static int rxn_map[27] = {1,2,6,7,8,3,9,10,11,12,4,13,14,5,15,16,17,18,19,20,0,21,22,23,24,25,26};

void GET_REACTION_MAP(int *rmap)
{
    for (int i=0; i<27; ++i) {
        rmap[i] = rxn_map[i];
    }
}

struct ReactionData* GetReactionData(int id)
{
    if (id<0 || id>=27) {
        printf("GetReactionData: Bad reaction id = %d",id);
        abort();
    };
    return &(R[rxn_map[id]]);
}

struct ReactionData* GetDefaultReactionData(int id)
{
    if (id<0 || id>=27) {
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
    if (id<0 || id>=27) {
        printf("SetReactionData: Bad reaction id = %d",id);
        abort();
    }
    R[rxn_map[id]] = *rhs;
    CopyReactionDataToTranspose(rxn_map[id],rhs);
}

/* Initializes static database */
void CKINIT()
{
    R[0].fwd_A     = 74000000000000;
    R[0].fwd_beta  = -0.37;
    R[0].fwd_Ea    = 0;
    R[0].low_A     = 2.3e+18;
    R[0].low_beta  = -0.90000000000000002;
    R[0].low_Ea    = -1700;
    R[0].troe_a    = 0.73460000000000003;
    R[0].troe_Tsss = 94;
    R[0].troe_Ts   = 1756;
    R[0].troe_Tss  = 5182;
    R[0].troe_len  = 4;
    R[0].prefactor_units  = 1.0000000000000002e-06;
    R[0].activation_units = 0.50321666580471969;
    R[0].phase_units      = 1e-12;
    R[0].is_PD = 1;

    R[1].fwd_A     = 1.2e+17;
    R[1].fwd_beta  = -1;
    R[1].fwd_Ea    = 0;
    R[1].prefactor_units  = 1.0000000000000002e-12;
    R[1].activation_units = 0.50321666580471969;
    R[1].phase_units      = 1e-12;
    R[1].is_PD = 0;

    R[2].fwd_A     = 5e+17;
    R[2].fwd_beta  = -1;
    R[2].fwd_Ea    = 0;
    R[2].prefactor_units  = 1.0000000000000002e-12;
    R[2].activation_units = 0.50321666580471969;
    R[2].phase_units      = 1e-12;
    R[2].is_PD = 0;

    R[3].fwd_A     = 2.8e+18;
    R[3].fwd_beta  = -0.85999999999999999;
    R[3].fwd_Ea    = 0;
    R[3].prefactor_units  = 1.0000000000000002e-12;
    R[3].activation_units = 0.50321666580471969;
    R[3].phase_units      = 1e-12;
    R[3].is_PD = 0;

    R[4].fwd_A     = 1e+18;
    R[4].fwd_beta  = -1;
    R[4].fwd_Ea    = 0;
    R[4].prefactor_units  = 1.0000000000000002e-12;
    R[4].activation_units = 0.50321666580471969;
    R[4].phase_units      = 1e-12;
    R[4].is_PD = 0;

    R[5].fwd_A     = 2.2e+22;
    R[5].fwd_beta  = -2;
    R[5].fwd_Ea    = 0;
    R[5].prefactor_units  = 1.0000000000000002e-12;
    R[5].activation_units = 0.50321666580471969;
    R[5].phase_units      = 1e-12;
    R[5].is_PD = 0;

    R[6].fwd_A     = 50000;
    R[6].fwd_beta  = 2.6699999999999999;
    R[6].fwd_Ea    = 6290;
    R[6].prefactor_units  = 1.0000000000000002e-06;
    R[6].activation_units = 0.50321666580471969;
    R[6].phase_units      = 1e-12;
    R[6].is_PD = 0;

    R[7].fwd_A     = 20000000000000;
    R[7].fwd_beta  = 0;
    R[7].fwd_Ea    = 0;
    R[7].prefactor_units  = 1.0000000000000002e-06;
    R[7].activation_units = 0.50321666580471969;
    R[7].phase_units      = 1e-12;
    R[7].is_PD = 0;

    R[8].fwd_A     = 9630000;
    R[8].fwd_beta  = 2;
    R[8].fwd_Ea    = 4000;
    R[8].prefactor_units  = 1.0000000000000002e-06;
    R[8].activation_units = 0.50321666580471969;
    R[8].phase_units      = 1e-12;
    R[8].is_PD = 0;

    R[9].fwd_A     = 3e+20;
    R[9].fwd_beta  = -1.72;
    R[9].fwd_Ea    = 0;
    R[9].prefactor_units  = 1.0000000000000002e-12;
    R[9].activation_units = 0.50321666580471969;
    R[9].phase_units      = 1e-18;
    R[9].is_PD = 0;

    R[10].fwd_A     = 9.38e+18;
    R[10].fwd_beta  = -0.76000000000000001;
    R[10].fwd_Ea    = 0;
    R[10].prefactor_units  = 1.0000000000000002e-12;
    R[10].activation_units = 0.50321666580471969;
    R[10].phase_units      = 1e-18;
    R[10].is_PD = 0;

    R[11].fwd_A     = 3.75e+20;
    R[11].fwd_beta  = -1.72;
    R[11].fwd_Ea    = 0;
    R[11].prefactor_units  = 1.0000000000000002e-12;
    R[11].activation_units = 0.50321666580471969;
    R[11].phase_units      = 1e-18;
    R[11].is_PD = 0;

    R[12].fwd_A     = 83000000000000;
    R[12].fwd_beta  = 0;
    R[12].fwd_Ea    = 14413;
    R[12].prefactor_units  = 1.0000000000000002e-06;
    R[12].activation_units = 0.50321666580471969;
    R[12].phase_units      = 1e-12;
    R[12].is_PD = 0;

    R[13].fwd_A     = 90000000000000000;
    R[13].fwd_beta  = -0.59999999999999998;
    R[13].fwd_Ea    = 0;
    R[13].prefactor_units  = 1.0000000000000002e-12;
    R[13].activation_units = 0.50321666580471969;
    R[13].phase_units      = 1e-18;
    R[13].is_PD = 0;

    R[14].fwd_A     = 6e+19;
    R[14].fwd_beta  = -1.25;
    R[14].fwd_Ea    = 0;
    R[14].prefactor_units  = 1.0000000000000002e-12;
    R[14].activation_units = 0.50321666580471969;
    R[14].phase_units      = 1e-18;
    R[14].is_PD = 0;

    R[15].fwd_A     = 3970000000000;
    R[15].fwd_beta  = 0;
    R[15].fwd_Ea    = 671;
    R[15].prefactor_units  = 1.0000000000000002e-06;
    R[15].activation_units = 0.50321666580471969;
    R[15].phase_units      = 1e-12;
    R[15].is_PD = 0;

    R[16].fwd_A     = 28000000000000;
    R[16].fwd_beta  = 0;
    R[16].fwd_Ea    = 1068;
    R[16].prefactor_units  = 1.0000000000000002e-06;
    R[16].activation_units = 0.50321666580471969;
    R[16].phase_units      = 1e-12;
    R[16].is_PD = 0;

    R[17].fwd_A     = 134000000000000;
    R[17].fwd_beta  = 0;
    R[17].fwd_Ea    = 635;
    R[17].prefactor_units  = 1.0000000000000002e-06;
    R[17].activation_units = 0.50321666580471969;
    R[17].phase_units      = 1e-12;
    R[17].is_PD = 0;

    R[18].fwd_A     = 12100000;
    R[18].fwd_beta  = 2;
    R[18].fwd_Ea    = 5200;
    R[18].prefactor_units  = 1.0000000000000002e-06;
    R[18].activation_units = 0.50321666580471969;
    R[18].phase_units      = 1e-12;
    R[18].is_PD = 0;

    R[19].fwd_A     = 10000000000000;
    R[19].fwd_beta  = 0;
    R[19].fwd_Ea    = 3600;
    R[19].prefactor_units  = 1.0000000000000002e-06;
    R[19].activation_units = 0.50321666580471969;
    R[19].phase_units      = 1e-12;
    R[19].is_PD = 0;

    R[20].fwd_A     = 216000000;
    R[20].fwd_beta  = 1.51;
    R[20].fwd_Ea    = 3430;
    R[20].prefactor_units  = 1.0000000000000002e-06;
    R[20].activation_units = 0.50321666580471969;
    R[20].phase_units      = 1e-12;
    R[20].is_PD = 0;

    R[21].fwd_A     = 35700;
    R[21].fwd_beta  = 2.3999999999999999;
    R[21].fwd_Ea    = -2110;
    R[21].prefactor_units  = 1.0000000000000002e-06;
    R[21].activation_units = 0.50321666580471969;
    R[21].phase_units      = 1e-12;
    R[21].is_PD = 0;

    R[22].fwd_A     = 29000000000000;
    R[22].fwd_beta  = 0;
    R[22].fwd_Ea    = -500;
    R[22].prefactor_units  = 1.0000000000000002e-06;
    R[22].activation_units = 0.50321666580471969;
    R[22].phase_units      = 1e-12;
    R[22].is_PD = 0;

    R[23].fwd_A     = 1750000000000;
    R[23].fwd_beta  = 0;
    R[23].fwd_Ea    = 320;
    R[23].prefactor_units  = 1.0000000000000002e-06;
    R[23].activation_units = 0.50321666580471969;
    R[23].phase_units      = 1e-12;
    R[23].is_PD = 0;

    R[24].fwd_A     = 580000000000000;
    R[24].fwd_beta  = 0;
    R[24].fwd_Ea    = 9560;
    R[24].prefactor_units  = 1.0000000000000002e-06;
    R[24].activation_units = 0.50321666580471969;
    R[24].phase_units      = 1e-12;
    R[24].is_PD = 0;

    R[25].fwd_A     = 130000000000;
    R[25].fwd_beta  = 0;
    R[25].fwd_Ea    = -1630;
    R[25].prefactor_units  = 1.0000000000000002e-06;
    R[25].activation_units = 0.50321666580471969;
    R[25].phase_units      = 1e-12;
    R[25].is_PD = 0;

    R[26].fwd_A     = 420000000000000;
    R[26].fwd_beta  = 0;
    R[26].fwd_Ea    = 12000;
    R[26].prefactor_units  = 1.0000000000000002e-06;
    R[26].activation_units = 0.50321666580471969;
    R[26].phase_units      = 1e-12;
    R[26].is_PD = 0;

    for (int i=0; i<27; i++)
    {
        R_DEF[i] = R[i];
        CopyReactionDataToTranspose(i,&(R[i]));
    }

}


/*A few mechanism parameters */
void CKINDX(int * iwrk, double * restrict rwrk, int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 3;
    *kk = 9;
    *ii = 27;
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
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void CKPY(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*O2 */
    YOW += y[4]*imw[4]; /*OH */
    YOW += y[5]*imw[5]; /*H2O */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
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
void CKPC(double * restrict rho, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict P)
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
    *P = *rho * 8.31451e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict rho)
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
    *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict rho)
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
void CKMMWX(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wtm)
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
void CKMMWC(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict wtm)
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
void CKYTX(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x)
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
void VCKYTX(int * restrict np, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x)
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
void CKYTCP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c)
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
void CKYTCR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c)
{
    for (int i = 0; i < 9; i++)
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
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    /*Now compute conversion */
    double XWinv = 1.0/XW;
    y[0] = x[0]*2.015940*XWinv; 
    y[1] = x[1]*1.007970*XWinv; 
    y[2] = x[2]*15.999400*XWinv; 
    y[3] = x[3]*31.998800*XWinv; 
    y[4] = x[4]*17.007370*XWinv; 
    y[5] = x[5]*18.015340*XWinv; 
    y[6] = x[6]*33.006770*XWinv; 
    y[7] = x[7]*34.014740*XWinv; 
    y[8] = x[8]*28.013400*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c)
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
void CKXTCR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c)
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
void CKCTX(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict x)
{
    int id; /*loop counter */
    double sumC = 0; 

    /*compute sum of c  */
    for (id = 0; id < 9; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    double sumCinv = 1.0/sumC;
    for (id = 0; id < 9; ++id) {
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
    CW += c[2]*15.999400; /*O */
    CW += c[3]*31.998800; /*O2 */
    CW += c[4]*17.007370; /*OH */
    CW += c[5]*18.015340; /*H2O */
    CW += c[6]*33.006770; /*HO2 */
    CW += c[7]*34.014740; /*H2O2 */
    CW += c[8]*28.013400; /*N2 */
    /*Now compute conversion */
    double CWinv = 1.0/CW;
    y[0] = c[0]*2.015940*CWinv; 
    y[1] = c[1]*1.007970*CWinv; 
    y[2] = c[2]*15.999400*CWinv; 
    y[3] = c[3]*31.998800*CWinv; 
    y[4] = c[4]*17.007370*CWinv; 
    y[5] = c[5]*18.015340*CWinv; 
    y[6] = c[6]*33.006770*CWinv; 
    y[7] = c[7]*34.014740*CWinv; 
    y[8] = c[8]*28.013400*CWinv; 

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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    cvms[2] *= 5.196763628636074e+06; /*O */
    cvms[3] *= 2.598381814318037e+06; /*O2 */
    cvms[4] *= 4.888768810227566e+06; /*OH */
    cvms[5] *= 4.615239012974499e+06; /*H2O */
    cvms[6] *= 2.519031701678171e+06; /*HO2 */
    cvms[7] *= 2.444384405113783e+06; /*H2O2 */
    cvms[8] *= 2.968047434442088e+06; /*N2 */
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
    cpms[2] *= 5.196763628636074e+06; /*O */
    cpms[3] *= 2.598381814318037e+06; /*O2 */
    cpms[4] *= 4.888768810227566e+06; /*OH */
    cpms[5] *= 4.615239012974499e+06; /*H2O */
    cpms[6] *= 2.519031701678171e+06; /*HO2 */
    cpms[7] *= 2.444384405113783e+06; /*H2O2 */
    cpms[8] *= 2.968047434442088e+06; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums)
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
void CKHMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hms)
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
void VCKHMS(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hms)
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
void CKGMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gms)
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
void CKAMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ams)
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
void CKSMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sms)
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
    sms[8] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl)
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
void CKCPBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cpbs)
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
void CKCVBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cvbl)
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
void CKCVBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cvbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[9]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H2 */
    result += cvor[1]*y[1]*imw[1]; /*H */
    result += cvor[2]*y[2]*imw[2]; /*O */
    result += cvor[3]*y[3]*imw[3]; /*O2 */
    result += cvor[4]*y[4]*imw[4]; /*OH */
    result += cvor[5]*y[5]*imw[5]; /*H2O */
    result += cvor[6]*y[6]*imw[6]; /*HO2 */
    result += cvor[7]*y[7]*imw[7]; /*H2O2 */
    result += cvor[8]*y[8]*imw[8]; /*N2 */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict hbml)
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
void CKHBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict hbms)
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
void CKUBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict ubml)
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
void CKUBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict ubms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double ums[9]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H2 */
    result += y[1]*ums[1]*imw[1]; /*H */
    result += y[2]*ums[2]*imw[2]; /*O */
    result += y[3]*ums[3]*imw[3]; /*O2 */
    result += y[4]*ums[4]*imw[4]; /*OH */
    result += y[5]*ums[5]*imw[5]; /*H2O */
    result += y[6]*ums[6]*imw[6]; /*HO2 */
    result += y[7]*ums[7]*imw[7]; /*H2O2 */
    result += y[8]*ums[8]*imw[8]; /*N2 */

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
    double sor[9]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 9; ++id) {
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
    double sor[9]; /* temporary storage */
    double x[9]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*O2 */
    YOW += y[4]*imw[4]; /*OH */
    YOW += y[5]*imw[5]; /*H2O */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
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
void CKGBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict gbms)
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
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*O2 */
    YOW += y[4]*imw[4]; /*OH */
    YOW += y[5]*imw[5]; /*H2O */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
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
void CKABML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict abml)
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
void CKABMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict abms)
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
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*O2 */
    YOW += y[4]*imw[4]; /*OH */
    YOW += y[5]*imw[5]; /*H2O */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
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
void CKWC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict wdot)
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
void CKWYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*O2 */
    YOW += y[4]*imw[4]; /*OH */
    YOW += y[5]*imw[5]; /*H2O */
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
void CKWXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
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
void CKWYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
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
void CKWXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
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
void CKQC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict qdot)
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
/*Given P, T, and mole fractions */
void CKKFKR(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict q_f, double * restrict q_r)
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
    for (id = 0; id < 27; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*O2 */
    YOW += y[4]*imw[4]; /*OH */
    YOW += y[5]*imw[5]; /*H2O */
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
    for (id = 0; id < 27; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
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
    for (id = 0; id < 27; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
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
    for (id = 0; id < 27; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
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
void CKNU(int * kdim, int * iwrk, double * restrict rwrk, int * nuki)
{
    int id; /*loop counter */
    int kd = (*kdim); 
    /*Zero nuki */
    for (id = 0; id < 9 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: 2 OH (+M) <=> H2O2 (+M) */
    nuki[ 4 * kd + 0 ] += -2 ;
    nuki[ 7 * kd + 0 ] += +1 ;

    /*reaction 2: 2 O + M <=> O2 + M */
    nuki[ 2 * kd + 1 ] += -2 ;
    nuki[ 3 * kd + 1 ] += +1 ;

    /*reaction 3: O + H + M <=> OH + M */
    nuki[ 2 * kd + 2 ] += -1 ;
    nuki[ 1 * kd + 2 ] += -1 ;
    nuki[ 4 * kd + 2 ] += +1 ;

    /*reaction 4: H + O2 + M <=> HO2 + M */
    nuki[ 1 * kd + 3 ] += -1 ;
    nuki[ 3 * kd + 3 ] += -1 ;
    nuki[ 6 * kd + 3 ] += +1 ;

    /*reaction 5: 2 H + M <=> H2 + M */
    nuki[ 1 * kd + 4 ] += -2 ;
    nuki[ 0 * kd + 4 ] += +1 ;

    /*reaction 6: H + OH + M <=> H2O + M */
    nuki[ 1 * kd + 5 ] += -1 ;
    nuki[ 4 * kd + 5 ] += -1 ;
    nuki[ 5 * kd + 5 ] += +1 ;

    /*reaction 7: O + H2 <=> H + OH */
    nuki[ 2 * kd + 6 ] += -1 ;
    nuki[ 0 * kd + 6 ] += -1 ;
    nuki[ 1 * kd + 6 ] += +1 ;
    nuki[ 4 * kd + 6 ] += +1 ;

    /*reaction 8: O + HO2 <=> OH + O2 */
    nuki[ 2 * kd + 7 ] += -1 ;
    nuki[ 6 * kd + 7 ] += -1 ;
    nuki[ 4 * kd + 7 ] += +1 ;
    nuki[ 3 * kd + 7 ] += +1 ;

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    nuki[ 2 * kd + 8 ] += -1 ;
    nuki[ 7 * kd + 8 ] += -1 ;
    nuki[ 4 * kd + 8 ] += +1 ;
    nuki[ 6 * kd + 8 ] += +1 ;

    /*reaction 10: H + 2 O2 <=> HO2 + O2 */
    nuki[ 1 * kd + 9 ] += -1 ;
    nuki[ 3 * kd + 9 ] += -2 ;
    nuki[ 6 * kd + 9 ] += +1 ;
    nuki[ 3 * kd + 9 ] += +1 ;

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    nuki[ 1 * kd + 10 ] += -1 ;
    nuki[ 3 * kd + 10 ] += -1 ;
    nuki[ 5 * kd + 10 ] += -1 ;
    nuki[ 6 * kd + 10 ] += +1 ;
    nuki[ 5 * kd + 10 ] += +1 ;

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    nuki[ 1 * kd + 11 ] += -1 ;
    nuki[ 3 * kd + 11 ] += -1 ;
    nuki[ 8 * kd + 11 ] += -1 ;
    nuki[ 6 * kd + 11 ] += +1 ;
    nuki[ 8 * kd + 11 ] += +1 ;

    /*reaction 13: H + O2 <=> O + OH */
    nuki[ 1 * kd + 12 ] += -1 ;
    nuki[ 3 * kd + 12 ] += -1 ;
    nuki[ 2 * kd + 12 ] += +1 ;
    nuki[ 4 * kd + 12 ] += +1 ;

    /*reaction 14: 2 H + H2 <=> 2 H2 */
    nuki[ 1 * kd + 13 ] += -2 ;
    nuki[ 0 * kd + 13 ] += -1 ;
    nuki[ 0 * kd + 13 ] += +2 ;

    /*reaction 15: 2 H + H2O <=> H2 + H2O */
    nuki[ 1 * kd + 14 ] += -2 ;
    nuki[ 5 * kd + 14 ] += -1 ;
    nuki[ 0 * kd + 14 ] += +1 ;
    nuki[ 5 * kd + 14 ] += +1 ;

    /*reaction 16: H + HO2 <=> O + H2O */
    nuki[ 1 * kd + 15 ] += -1 ;
    nuki[ 6 * kd + 15 ] += -1 ;
    nuki[ 2 * kd + 15 ] += +1 ;
    nuki[ 5 * kd + 15 ] += +1 ;

    /*reaction 17: H + HO2 <=> O2 + H2 */
    nuki[ 1 * kd + 16 ] += -1 ;
    nuki[ 6 * kd + 16 ] += -1 ;
    nuki[ 3 * kd + 16 ] += +1 ;
    nuki[ 0 * kd + 16 ] += +1 ;

    /*reaction 18: H + HO2 <=> 2 OH */
    nuki[ 1 * kd + 17 ] += -1 ;
    nuki[ 6 * kd + 17 ] += -1 ;
    nuki[ 4 * kd + 17 ] += +2 ;

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    nuki[ 1 * kd + 18 ] += -1 ;
    nuki[ 7 * kd + 18 ] += -1 ;
    nuki[ 6 * kd + 18 ] += +1 ;
    nuki[ 0 * kd + 18 ] += +1 ;

    /*reaction 20: H + H2O2 <=> OH + H2O */
    nuki[ 1 * kd + 19 ] += -1 ;
    nuki[ 7 * kd + 19 ] += -1 ;
    nuki[ 4 * kd + 19 ] += +1 ;
    nuki[ 5 * kd + 19 ] += +1 ;

    /*reaction 21: OH + H2 <=> H + H2O */
    nuki[ 4 * kd + 20 ] += -1 ;
    nuki[ 0 * kd + 20 ] += -1 ;
    nuki[ 1 * kd + 20 ] += +1 ;
    nuki[ 5 * kd + 20 ] += +1 ;

    /*reaction 22: 2 OH <=> O + H2O */
    nuki[ 4 * kd + 21 ] += -2 ;
    nuki[ 2 * kd + 21 ] += +1 ;
    nuki[ 5 * kd + 21 ] += +1 ;

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    nuki[ 4 * kd + 22 ] += -1 ;
    nuki[ 6 * kd + 22 ] += -1 ;
    nuki[ 3 * kd + 22 ] += +1 ;
    nuki[ 5 * kd + 22 ] += +1 ;

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    nuki[ 4 * kd + 23 ] += -1 ;
    nuki[ 7 * kd + 23 ] += -1 ;
    nuki[ 6 * kd + 23 ] += +1 ;
    nuki[ 5 * kd + 23 ] += +1 ;

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    nuki[ 4 * kd + 24 ] += -1 ;
    nuki[ 7 * kd + 24 ] += -1 ;
    nuki[ 6 * kd + 24 ] += +1 ;
    nuki[ 5 * kd + 24 ] += +1 ;

    /*reaction 26: 2 HO2 <=> O2 + H2O2 */
    nuki[ 6 * kd + 25 ] += -2 ;
    nuki[ 3 * kd + 25 ] += +1 ;
    nuki[ 7 * kd + 25 ] += +1 ;

    /*reaction 27: 2 HO2 <=> O2 + H2O2 */
    nuki[ 6 * kd + 26 ] += -2 ;
    nuki[ 3 * kd + 26 ] += +1 ;
    nuki[ 7 * kd + 26 ] += +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * restrict rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < kd * 9; ++ id) {
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
void CKABE(int * iwrk, double * restrict rwrk, double * restrict a, double * restrict b, double * restrict e)
{
    for (int i=0; i<27; ++i) {
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
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: 2 O + M <=> O2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: O + H + M <=> OH + M */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + O2 + M <=> HO2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: 2 H + M <=> H2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + OH + M <=> H2O + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + H2 <=> H + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + HO2 <=> OH + O2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: H + 2 O2 <=> HO2 + O2 */
    eqcon[9] *= 1e+06; 

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    eqcon[10] *= 1e+06; 

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: H + O2 <=> O + OH */
    /*eqcon[12] *= 1;  */

    /*reaction 14: 2 H + H2 <=> 2 H2 */
    eqcon[13] *= 1e+06; 

    /*reaction 15: 2 H + H2O <=> H2 + H2O */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H + HO2 <=> O + H2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + HO2 <=> O2 + H2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H + HO2 <=> 2 OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H + H2O2 <=> OH + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: OH + H2 <=> H + H2O */
    /*eqcon[20] *= 1;  */

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
void CKEQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: 2 O + M <=> O2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: O + H + M <=> OH + M */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + O2 + M <=> HO2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: 2 H + M <=> H2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + OH + M <=> H2O + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + H2 <=> H + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + HO2 <=> OH + O2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: H + 2 O2 <=> HO2 + O2 */
    eqcon[9] *= 1e+06; 

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    eqcon[10] *= 1e+06; 

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: H + O2 <=> O + OH */
    /*eqcon[12] *= 1;  */

    /*reaction 14: 2 H + H2 <=> 2 H2 */
    eqcon[13] *= 1e+06; 

    /*reaction 15: 2 H + H2O <=> H2 + H2O */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H + HO2 <=> O + H2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + HO2 <=> O2 + H2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H + HO2 <=> 2 OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H + H2O2 <=> OH + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: OH + H2 <=> H + H2O */
    /*eqcon[20] *= 1;  */

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
void CKEQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: 2 O + M <=> O2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: O + H + M <=> OH + M */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + O2 + M <=> HO2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: 2 H + M <=> H2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + OH + M <=> H2O + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + H2 <=> H + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + HO2 <=> OH + O2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: H + 2 O2 <=> HO2 + O2 */
    eqcon[9] *= 1e+06; 

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    eqcon[10] *= 1e+06; 

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: H + O2 <=> O + OH */
    /*eqcon[12] *= 1;  */

    /*reaction 14: 2 H + H2 <=> 2 H2 */
    eqcon[13] *= 1e+06; 

    /*reaction 15: 2 H + H2O <=> H2 + H2O */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H + HO2 <=> O + H2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + HO2 <=> O2 + H2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H + HO2 <=> 2 OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H + H2O2 <=> OH + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: OH + H2 <=> H + H2O */
    /*eqcon[20] *= 1;  */

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
void CKEQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: 2 O + M <=> O2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: O + H + M <=> OH + M */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + O2 + M <=> HO2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: 2 H + M <=> H2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + OH + M <=> H2O + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + H2 <=> H + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + HO2 <=> OH + O2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: H + 2 O2 <=> HO2 + O2 */
    eqcon[9] *= 1e+06; 

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    eqcon[10] *= 1e+06; 

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: H + O2 <=> O + OH */
    /*eqcon[12] *= 1;  */

    /*reaction 14: 2 H + H2 <=> 2 H2 */
    eqcon[13] *= 1e+06; 

    /*reaction 15: 2 H + H2O <=> H2 + H2O */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H + HO2 <=> O + H2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + HO2 <=> O2 + H2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H + HO2 <=> 2 OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H + H2O2 <=> OH + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: OH + H2 <=> H + H2O */
    /*eqcon[20] *= 1;  */

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
void CKEQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: 2 O + M <=> O2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: O + H + M <=> OH + M */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + O2 + M <=> HO2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: 2 H + M <=> H2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + OH + M <=> H2O + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + H2 <=> H + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + HO2 <=> OH + O2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: H + 2 O2 <=> HO2 + O2 */
    eqcon[9] *= 1e+06; 

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    eqcon[10] *= 1e+06; 

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: H + O2 <=> O + OH */
    /*eqcon[12] *= 1;  */

    /*reaction 14: 2 H + H2 <=> 2 H2 */
    eqcon[13] *= 1e+06; 

    /*reaction 15: 2 H + H2O <=> H2 + H2O */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H + HO2 <=> O + H2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + HO2 <=> O2 + H2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H + HO2 <=> 2 OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H + H2O2 <=> OH + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: OH + H2 <=> H + H2O */
    /*eqcon[20] *= 1;  */

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

static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[27];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[27];
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

    double qdot, q_f[27], q_r[27];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 9; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[4] -= 2 * qdot;
    wdot[7] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[2] -= 2 * qdot;
    wdot[3] += qdot;

    qdot = q_f[2]-q_r[2];
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[0] += qdot;
    wdot[1] -= 2 * qdot;

    qdot = q_f[5]-q_r[5];
    wdot[1] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[6]-q_r[6];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[8]-q_r[8];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[9]-q_r[9];
    wdot[1] -= qdot;
    wdot[3] -= 2 * qdot;
    wdot[3] += qdot;
    wdot[6] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[8] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[12]-q_r[12];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[13]-q_r[13];
    wdot[0] -= qdot;
    wdot[0] += 2 * qdot;
    wdot[1] -= 2 * qdot;

    qdot = q_f[14]-q_r[14];
    wdot[0] += qdot;
    wdot[1] -= 2 * qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[16]-q_r[16];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[17]-q_r[17];
    wdot[1] -= qdot;
    wdot[4] += 2 * qdot;
    wdot[6] -= qdot;

    qdot = q_f[18]-q_r[18];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[19]-q_r[19];
    wdot[1] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[21]-q_r[21];
    wdot[2] += qdot;
    wdot[4] -= 2 * qdot;
    wdot[5] += qdot;

    qdot = q_f[22]-q_r[22];
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[23]-q_r[23];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[24]-q_r[24];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[25]-q_r[25];
    wdot[3] += qdot;
    wdot[6] -= 2 * qdot;
    wdot[7] += qdot;

    qdot = q_f[26]-q_r[26];
    wdot[3] += qdot;
    wdot[6] -= 2 * qdot;
    wdot[7] += qdot;

    return;
}

void comp_k_f(double * restrict tc, double invT, double * restrict k_f)
{
#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<27; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(double * restrict tc, double invT, double * restrict Kc)
{
    /*compute the Gibbs free energy */
    double g_RT[9];
    gibbs(g_RT, tc);

    Kc[0] = 2*g_RT[4] - g_RT[7];
    Kc[1] = 2*g_RT[2] - g_RT[3];
    Kc[2] = g_RT[1] + g_RT[2] - g_RT[4];
    Kc[3] = g_RT[1] + g_RT[3] - g_RT[6];
    Kc[4] = -g_RT[0] + 2*g_RT[1];
    Kc[5] = g_RT[1] + g_RT[4] - g_RT[5];
    Kc[6] = g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4];
    Kc[7] = g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6];
    Kc[8] = g_RT[2] - g_RT[4] - g_RT[6] + g_RT[7];
    Kc[9] = g_RT[1] + 2*g_RT[3] - g_RT[3] - g_RT[6];
    Kc[10] = g_RT[1] + g_RT[3] + g_RT[5] - g_RT[5] - g_RT[6];
    Kc[11] = g_RT[1] + g_RT[3] - g_RT[6] + g_RT[8] - g_RT[8];
    Kc[12] = g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4];
    Kc[13] = g_RT[0] - 2*g_RT[0] + 2*g_RT[1];
    Kc[14] = -g_RT[0] + 2*g_RT[1] + g_RT[5] - g_RT[5];
    Kc[15] = g_RT[1] - g_RT[2] - g_RT[5] + g_RT[6];
    Kc[16] = -g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6];
    Kc[17] = g_RT[1] - 2*g_RT[4] + g_RT[6];
    Kc[18] = -g_RT[0] + g_RT[1] - g_RT[6] + g_RT[7];
    Kc[19] = g_RT[1] - g_RT[4] - g_RT[5] + g_RT[7];
    Kc[20] = g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5];
    Kc[21] = -g_RT[2] + 2*g_RT[4] - g_RT[5];
    Kc[22] = -g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6];
    Kc[23] = g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7];
    Kc[24] = g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7];
    Kc[25] = -g_RT[3] + 2*g_RT[6] - g_RT[7];
    Kc[26] = -g_RT[3] + 2*g_RT[6] - g_RT[7];

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<27; ++i) {
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
    Kc[9] *= refCinv;
    Kc[10] *= refCinv;
    Kc[11] *= refCinv;
    Kc[13] *= refCinv;
    Kc[14] *= refCinv;

    return;
}

void comp_qfqr(double * restrict qf, double * restrict qr, double * restrict sc, double * restrict tc, double invT)
{

    /*reaction 1: 2 OH (+M) <=> H2O2 (+M) */
    qf[0] = sc[4]*sc[4];
    qr[0] = sc[7];

    /*reaction 2: 2 O + M <=> O2 + M */
    qf[1] = sc[2]*sc[2];
    qr[1] = sc[3];

    /*reaction 3: O + H + M <=> OH + M */
    qf[2] = sc[1]*sc[2];
    qr[2] = sc[4];

    /*reaction 4: H + O2 + M <=> HO2 + M */
    qf[3] = sc[1]*sc[3];
    qr[3] = sc[6];

    /*reaction 5: 2 H + M <=> H2 + M */
    qf[4] = sc[1]*sc[1];
    qr[4] = sc[0];

    /*reaction 6: H + OH + M <=> H2O + M */
    qf[5] = sc[1]*sc[4];
    qr[5] = sc[5];

    /*reaction 7: O + H2 <=> H + OH */
    qf[6] = sc[0]*sc[2];
    qr[6] = sc[1]*sc[4];

    /*reaction 8: O + HO2 <=> OH + O2 */
    qf[7] = sc[2]*sc[6];
    qr[7] = sc[3]*sc[4];

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    qf[8] = sc[2]*sc[7];
    qr[8] = sc[4]*sc[6];

    /*reaction 10: H + 2 O2 <=> HO2 + O2 */
    qf[9] = sc[1]*sc[3]*sc[3];
    qr[9] = sc[3]*sc[6];

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    qf[10] = sc[1]*sc[3]*sc[5];
    qr[10] = sc[5]*sc[6];

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    qf[11] = sc[1]*sc[3]*sc[8];
    qr[11] = sc[6]*sc[8];

    /*reaction 13: H + O2 <=> O + OH */
    qf[12] = sc[1]*sc[3];
    qr[12] = sc[2]*sc[4];

    /*reaction 14: 2 H + H2 <=> 2 H2 */
    qf[13] = sc[0]*sc[1]*sc[1];
    qr[13] = sc[0]*sc[0];

    /*reaction 15: 2 H + H2O <=> H2 + H2O */
    qf[14] = sc[1]*sc[1]*sc[5];
    qr[14] = sc[0]*sc[5];

    /*reaction 16: H + HO2 <=> O + H2O */
    qf[15] = sc[1]*sc[6];
    qr[15] = sc[2]*sc[5];

    /*reaction 17: H + HO2 <=> O2 + H2 */
    qf[16] = sc[1]*sc[6];
    qr[16] = sc[0]*sc[3];

    /*reaction 18: H + HO2 <=> 2 OH */
    qf[17] = sc[1]*sc[6];
    qr[17] = sc[4]*sc[4];

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    qf[18] = sc[1]*sc[7];
    qr[18] = sc[0]*sc[6];

    /*reaction 20: H + H2O2 <=> OH + H2O */
    qf[19] = sc[1]*sc[7];
    qr[19] = sc[4]*sc[5];

    /*reaction 21: OH + H2 <=> H + H2O */
    qf[20] = sc[0]*sc[4];
    qr[20] = sc[1]*sc[5];

    /*reaction 22: 2 OH <=> O + H2O */
    qf[21] = sc[4]*sc[4];
    qr[21] = sc[2]*sc[5];

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    qf[22] = sc[4]*sc[6];
    qr[22] = sc[3]*sc[5];

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    qf[23] = sc[4]*sc[7];
    qr[23] = sc[5]*sc[6];

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    qf[24] = sc[4]*sc[7];
    qr[24] = sc[5]*sc[6];

    /*reaction 26: 2 HO2 <=> O2 + H2O2 */
    qf[25] = sc[6]*sc[6];
    qr[25] = sc[3]*sc[7];

    /*reaction 27: 2 HO2 <=> O2 + H2O2 */
    qf[26] = sc[6]*sc[6];
    qr[26] = sc[3]*sc[7];

    double T = tc[1];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 9; ++i) {
        mixture += sc[i];
    }

    double Corr[27];
    for (int i = 0; i < 27; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        double alpha[1];
        alpha[0] = mixture + sc[0] + 5*sc[5];
        for (int i=0; i<1; i++)
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

    /* simple three-body correction */
    {
        double alpha;
        alpha = mixture + 1.3999999999999999*sc[0] + 14.4*sc[5];
        Corr[1] = alpha;
        alpha = mixture + sc[0] + 5*sc[5];
        Corr[2] = alpha;
        alpha = mixture - sc[3] - sc[5] - sc[8];
        Corr[3] = alpha;
        alpha = mixture - sc[0] - sc[5];
        Corr[4] = alpha;
        alpha = mixture - 0.27000000000000002*sc[0] + 2.6499999999999999*sc[5];
        Corr[5] = alpha;
    }

    for (int i=0; i<27; i++)
    {
        qf[i] *= Corr[i] * k_f_save[i];
        qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];
    }

    return;
}


/*compute the production rate for each species */
void vproductionRate(int npt, double * restrict wdot, double * restrict sc, double * restrict T)
{
    double k_f_s[27*npt], Kc_s[27*npt], mixture[npt], g_RT[9*npt];
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

    for (int n=0; n<9; n++) {
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
    }
}

void vcomp_gibbs(int npt, double * restrict g_RT, double * restrict tc)
{
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

        Kc_s[0*npt+i] = refCinv * exp((2 * g_RT[4*npt+i]) - (g_RT[7*npt+i]));
        Kc_s[1*npt+i] = refCinv * exp((2 * g_RT[2*npt+i]) - (g_RT[3*npt+i]));
        Kc_s[2*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[2*npt+i]) - (g_RT[4*npt+i]));
        Kc_s[3*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[6*npt+i]));
        Kc_s[4*npt+i] = refCinv * exp((2 * g_RT[1*npt+i]) - (g_RT[0*npt+i]));
        Kc_s[5*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[4*npt+i]) - (g_RT[5*npt+i]));
        Kc_s[6*npt+i] = exp((g_RT[0*npt+i] + g_RT[2*npt+i]) - (g_RT[1*npt+i] + g_RT[4*npt+i]));
        Kc_s[7*npt+i] = exp((g_RT[2*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[4*npt+i]));
        Kc_s[8*npt+i] = exp((g_RT[2*npt+i] + g_RT[7*npt+i]) - (g_RT[4*npt+i] + g_RT[6*npt+i]));
        Kc_s[9*npt+i] = refCinv * exp((g_RT[1*npt+i] + 2 * g_RT[3*npt+i]) - (g_RT[3*npt+i] + g_RT[6*npt+i]));
        Kc_s[10*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i] + g_RT[5*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[11*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i] + g_RT[8*npt+i]) - (g_RT[6*npt+i] + g_RT[8*npt+i]));
        Kc_s[12*npt+i] = exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[2*npt+i] + g_RT[4*npt+i]));
        Kc_s[13*npt+i] = refCinv * exp((g_RT[0*npt+i] + 2 * g_RT[1*npt+i]) - (2 * g_RT[0*npt+i]));
        Kc_s[14*npt+i] = refCinv * exp((2 * g_RT[1*npt+i] + g_RT[5*npt+i]) - (g_RT[0*npt+i] + g_RT[5*npt+i]));
        Kc_s[15*npt+i] = exp((g_RT[1*npt+i] + g_RT[6*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i]));
        Kc_s[16*npt+i] = exp((g_RT[1*npt+i] + g_RT[6*npt+i]) - (g_RT[0*npt+i] + g_RT[3*npt+i]));
        Kc_s[17*npt+i] = exp((g_RT[1*npt+i] + g_RT[6*npt+i]) - (2 * g_RT[4*npt+i]));
        Kc_s[18*npt+i] = exp((g_RT[1*npt+i] + g_RT[7*npt+i]) - (g_RT[0*npt+i] + g_RT[6*npt+i]));
        Kc_s[19*npt+i] = exp((g_RT[1*npt+i] + g_RT[7*npt+i]) - (g_RT[4*npt+i] + g_RT[5*npt+i]));
        Kc_s[20*npt+i] = exp((g_RT[0*npt+i] + g_RT[4*npt+i]) - (g_RT[1*npt+i] + g_RT[5*npt+i]));
        Kc_s[21*npt+i] = exp((2 * g_RT[4*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i]));
        Kc_s[22*npt+i] = exp((g_RT[4*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[5*npt+i]));
        Kc_s[23*npt+i] = exp((g_RT[4*npt+i] + g_RT[7*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[24*npt+i] = exp((g_RT[4*npt+i] + g_RT[7*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[25*npt+i] = exp((2 * g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[7*npt+i]));
        Kc_s[26*npt+i] = exp((2 * g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[7*npt+i]));
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

        /*reaction 1: 2 OH (+M) <=> H2O2 (+M) */
        phi_f = sc[4*npt+i]*sc[4*npt+i];
        alpha = mixture[i] + sc[0*npt+i] + 5*sc[5*npt+i];
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
        phi_r = sc[7*npt+i];
        Kc = Kc_s[0*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= 2 * qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 2: 2 O + M <=> O2 + M */
        phi_f = sc[2*npt+i]*sc[2*npt+i];
        alpha = mixture[i] + 1.3999999999999999*sc[0*npt+i] + 14.4*sc[5*npt+i];
        k_f = alpha * k_f_s[1*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i];
        Kc = Kc_s[1*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= 2 * qdot;
        wdot[3*npt+i] += qdot;

        /*reaction 3: O + H + M <=> OH + M */
        phi_f = sc[1*npt+i]*sc[2*npt+i];
        alpha = mixture[i] + sc[0*npt+i] + 5*sc[5*npt+i];
        k_f = alpha * k_f_s[2*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i];
        Kc = Kc_s[2*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 4: H + O2 + M <=> HO2 + M */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        alpha = mixture[i] + -sc[3*npt+i] + -sc[5*npt+i] + -sc[8*npt+i];
        k_f = alpha * k_f_s[3*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 5: 2 H + M <=> H2 + M */
        phi_f = sc[1*npt+i]*sc[1*npt+i];
        alpha = mixture[i] + -sc[0*npt+i] + -sc[5*npt+i];
        k_f = alpha * k_f_s[4*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i];
        Kc = Kc_s[4*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= 2 * qdot;

        /*reaction 6: H + OH + M <=> H2O + M */
        phi_f = sc[1*npt+i]*sc[4*npt+i];
        alpha = mixture[i] + -0.27000000000000002*sc[0*npt+i] + 2.6499999999999999*sc[5*npt+i];
        k_f = alpha * k_f_s[5*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i];
        Kc = Kc_s[5*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 7: O + H2 <=> H + OH */
        phi_f = sc[0*npt+i]*sc[2*npt+i];
        k_f = k_f_s[6*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[4*npt+i];
        Kc = Kc_s[6*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 8: O + HO2 <=> OH + O2 */
        phi_f = sc[2*npt+i]*sc[6*npt+i];
        k_f = k_f_s[7*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[4*npt+i];
        Kc = Kc_s[7*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 9: O + H2O2 <=> OH + HO2 */
        phi_f = sc[2*npt+i]*sc[7*npt+i];
        k_f = k_f_s[8*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[6*npt+i];
        Kc = Kc_s[8*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 10: H + 2 O2 <=> HO2 + O2 */
        phi_f = sc[1*npt+i]*sc[3*npt+i]*sc[3*npt+i];
        k_f = k_f_s[9*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[6*npt+i];
        Kc = Kc_s[9*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= 2 * qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
        phi_f = sc[1*npt+i]*sc[3*npt+i]*sc[5*npt+i];
        k_f = k_f_s[10*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[10*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
        phi_f = sc[1*npt+i]*sc[3*npt+i]*sc[8*npt+i];
        k_f = k_f_s[11*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[8*npt+i];
        Kc = Kc_s[11*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 13: H + O2 <=> O + OH */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        k_f = k_f_s[12*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[4*npt+i];
        Kc = Kc_s[12*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 14: 2 H + H2 <=> 2 H2 */
        phi_f = sc[0*npt+i]*sc[1*npt+i]*sc[1*npt+i];
        k_f = k_f_s[13*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[0*npt+i];
        Kc = Kc_s[13*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[0*npt+i] += 2 * qdot;
        wdot[1*npt+i] -= 2 * qdot;

        /*reaction 15: 2 H + H2O <=> H2 + H2O */
        phi_f = sc[1*npt+i]*sc[1*npt+i]*sc[5*npt+i];
        k_f = k_f_s[14*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[5*npt+i];
        Kc = Kc_s[14*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= 2 * qdot;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 16: H + HO2 <=> O + H2O */
        phi_f = sc[1*npt+i]*sc[6*npt+i];
        k_f = k_f_s[15*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[5*npt+i];
        Kc = Kc_s[15*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 17: H + HO2 <=> O2 + H2 */
        phi_f = sc[1*npt+i]*sc[6*npt+i];
        k_f = k_f_s[16*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[3*npt+i];
        Kc = Kc_s[16*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 18: H + HO2 <=> 2 OH */
        phi_f = sc[1*npt+i]*sc[6*npt+i];
        k_f = k_f_s[17*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[4*npt+i];
        Kc = Kc_s[17*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += 2 * qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 19: H + H2O2 <=> HO2 + H2 */
        phi_f = sc[1*npt+i]*sc[7*npt+i];
        k_f = k_f_s[18*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[6*npt+i];
        Kc = Kc_s[18*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 20: H + H2O2 <=> OH + H2O */
        phi_f = sc[1*npt+i]*sc[7*npt+i];
        k_f = k_f_s[19*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[5*npt+i];
        Kc = Kc_s[19*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 21: OH + H2 <=> H + H2O */
        phi_f = sc[0*npt+i]*sc[4*npt+i];
        k_f = k_f_s[20*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[5*npt+i];
        Kc = Kc_s[20*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 22: 2 OH <=> O + H2O */
        phi_f = sc[4*npt+i]*sc[4*npt+i];
        k_f = k_f_s[21*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[5*npt+i];
        Kc = Kc_s[21*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= 2 * qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 23: OH + HO2 <=> O2 + H2O */
        phi_f = sc[4*npt+i]*sc[6*npt+i];
        k_f = k_f_s[22*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[5*npt+i];
        Kc = Kc_s[22*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 24: OH + H2O2 <=> HO2 + H2O */
        phi_f = sc[4*npt+i]*sc[7*npt+i];
        k_f = k_f_s[23*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[23*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 25: OH + H2O2 <=> HO2 + H2O */
        phi_f = sc[4*npt+i]*sc[7*npt+i];
        k_f = k_f_s[24*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[24*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 26: 2 HO2 <=> O2 + H2O2 */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[25*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[7*npt+i];
        Kc = Kc_s[25*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= 2 * qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 27: 2 HO2 <=> O2 + H2O2 */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[26*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[7*npt+i];
        Kc = Kc_s[26*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= 2 * qdot;
        wdot[7*npt+i] += qdot;
    }
}

/*compute the reaction Jacobian */
void DWDOT(double * restrict J, double * restrict sc, double * restrict Tp, int * consP)
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
    double refCinv = 1.0 / refC;

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

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[9];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    /*reaction 1: 2 OH (+M) <=> H2O2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + sc[0] + 5*sc[5];
    /* forward */
    phi_f = sc[4]*sc[4];
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
    phi_r = sc[7];
    Kc = refCinv * exp(2*g_RT[4] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[4]) + (h_RT[7]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[4] -= 2 * q; /* OH */
    wdot[7] += q; /* H2O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci =  dcdc_fac;
        J[4] += -2 * dqdci;           /* dwdot[OH]/d[H2] */
        J[7] += dqdci;                /* dwdot[H2O2]/d[H2] */
        /* d()/d[OH] */
        dqdci =  + k_f*2*sc[4];
        J[44] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
        J[47] += dqdci;               /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O] */
        dqdci =  5*dcdc_fac;
        J[54] += -2 * dqdci;          /* dwdot[OH]/d[H2O] */
        J[57] += dqdci;               /* dwdot[H2O2]/d[H2O] */
        /* d()/d[H2O2] */
        dqdci =  - k_r;
        J[74] += -2 * dqdci;          /* dwdot[OH]/d[H2O2] */
        J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    }
    else {
        dqdc[0] =  2*dcdc_fac;
        dqdc[1] =  dcdc_fac;
        dqdc[2] =  dcdc_fac;
        dqdc[3] =  dcdc_fac;
        dqdc[4] =  dcdc_fac + k_f*2*sc[4];
        dqdc[5] =  6*dcdc_fac;
        dqdc[6] =  dcdc_fac;
        dqdc[7] =  dcdc_fac - k_r;
        dqdc[8] =  dcdc_fac;
        for (int k=0; k<9; k++) {
            J[10*k+4] += -2 * dqdc[k];
            J[10*k+7] += dqdc[k];
        }
    }
    J[94] += -2 * dqdT; /* dwdot[OH]/dT */
    J[97] += dqdT; /* dwdot[H2O2]/dT */

    /*reaction 2: 2 O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + 1.3999999999999999*sc[0] + 14.4*sc[5];
    /* forward */
    phi_f = sc[2]*sc[2];
    k_f = prefactor_units[1] * fwd_A[1]
                * exp(fwd_beta[1] * tc[0] - activation_units[1] * fwd_Ea[1] * invT);
    dlnkfdT = fwd_beta[1] * invT + activation_units[1] * fwd_Ea[1] * invT2;
    /* reverse */
    phi_r = sc[3];
    Kc = refCinv * exp(2*g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[2]) + (h_RT[3]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[3] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci =  1.3999999999999999*q_nocor;
        J[2] += -2 * dqdci;           /* dwdot[O]/d[H2] */
        J[3] += dqdci;                /* dwdot[O2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*2*sc[2];
        J[22] += -2 * dqdci;          /* dwdot[O]/d[O] */
        J[23] += dqdci;               /* dwdot[O2]/d[O] */
        /* d()/d[O2] */
        dqdci =  - k_r;
        J[32] += -2 * dqdci;          /* dwdot[O]/d[O2] */
        J[33] += dqdci;               /* dwdot[O2]/d[O2] */
        /* d()/d[H2O] */
        dqdci =  14.4*q_nocor;
        J[52] += -2 * dqdci;          /* dwdot[O]/d[H2O] */
        J[53] += dqdci;               /* dwdot[O2]/d[H2O] */
    }
    else {
        dqdc[0] =  2.3999999999999999*q_nocor;
        dqdc[1] =  q_nocor;
        dqdc[2] =  q_nocor + k_f*2*sc[2];
        dqdc[3] =  q_nocor - k_r;
        dqdc[4] =  q_nocor;
        dqdc[5] =  15.4*q_nocor;
        dqdc[6] =  q_nocor;
        dqdc[7] =  q_nocor;
        dqdc[8] =  q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+2] += -2 * dqdc[k];
            J[10*k+3] += dqdc[k];
        }
    }
    J[92] += -2 * dqdT; /* dwdot[O]/dT */
    J[93] += dqdT; /* dwdot[O2]/dT */

    /*reaction 3: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + sc[0] + 5*sc[5];
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = prefactor_units[2] * fwd_A[2]
                * exp(fwd_beta[2] * tc[0] - activation_units[2] * fwd_Ea[2] * invT);
    dlnkfdT = fwd_beta[2] * invT + activation_units[2] * fwd_Ea[2] * invT2;
    /* reverse */
    phi_r = sc[4];
    Kc = refCinv * exp(g_RT[1] + g_RT[2] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[4]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci =  q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[2] -= dqdci;                /* dwdot[O]/d[H2] */
        J[4] += dqdci;                /* dwdot[OH]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[2];
        J[11] -= dqdci;               /* dwdot[H]/d[H] */
        J[12] -= dqdci;               /* dwdot[O]/d[H] */
        J[14] += dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[1];
        J[21] -= dqdci;               /* dwdot[H]/d[O] */
        J[22] -= dqdci;               /* dwdot[O]/d[O] */
        J[24] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[41] -= dqdci;               /* dwdot[H]/d[OH] */
        J[42] -= dqdci;               /* dwdot[O]/d[OH] */
        J[44] += dqdci;               /* dwdot[OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci =  5*q_nocor;
        J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[52] -= dqdci;               /* dwdot[O]/d[H2O] */
        J[54] += dqdci;               /* dwdot[OH]/d[H2O] */
    }
    else {
        dqdc[0] =  2*q_nocor;
        dqdc[1] =  q_nocor + k_f*sc[2];
        dqdc[2] =  q_nocor + k_f*sc[1];
        dqdc[3] =  q_nocor;
        dqdc[4] =  q_nocor - k_r;
        dqdc[5] =  6*q_nocor;
        dqdc[6] =  q_nocor;
        dqdc[7] =  q_nocor;
        dqdc[8] =  q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+1] -= dqdc[k];
            J[10*k+2] -= dqdc[k];
            J[10*k+4] += dqdc[k];
        }
    }
    J[91] -= dqdT; /* dwdot[H]/dT */
    J[92] -= dqdT; /* dwdot[O]/dT */
    J[94] += dqdT; /* dwdot[OH]/dT */

    /*reaction 4: H + O2 + M <=> HO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture - sc[3] - sc[5] - sc[8];
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[3] * fwd_A[3]
                * exp(fwd_beta[3] * tc[0] - activation_units[3] * fwd_Ea[3] * invT);
    dlnkfdT = fwd_beta[3] * invT + activation_units[3] * fwd_Ea[3] * invT2;
    /* reverse */
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[6]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[3];
        J[11] -= dqdci;               /* dwdot[H]/d[H] */
        J[13] -= dqdci;               /* dwdot[O2]/d[H] */
        J[16] += dqdci;               /* dwdot[HO2]/d[H] */
        /* d()/d[O2] */
        dqdci =  -q_nocor + k_f*sc[1];
        J[31] -= dqdci;               /* dwdot[H]/d[O2] */
        J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[36] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci =  -q_nocor;
        J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[53] -= dqdci;               /* dwdot[O2]/d[H2O] */
        J[56] += dqdci;               /* dwdot[HO2]/d[H2O] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
        J[63] -= dqdci;               /* dwdot[O2]/d[HO2] */
        J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
        /* d()/d[N2] */
        dqdci =  -q_nocor;
        J[81] -= dqdci;               /* dwdot[H]/d[N2] */
        J[83] -= dqdci;               /* dwdot[O2]/d[N2] */
        J[86] += dqdci;               /* dwdot[HO2]/d[N2] */
    }
    else {
        dqdc[0] =  q_nocor;
        dqdc[1] =  q_nocor + k_f*sc[3];
        dqdc[2] =  q_nocor;
        dqdc[3] =  + k_f*sc[1];
        dqdc[4] =  q_nocor;
        dqdc[6] =  q_nocor - k_r;
        dqdc[7] =  q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+1] -= dqdc[k];
            J[10*k+3] -= dqdc[k];
            J[10*k+6] += dqdc[k];
        }
    }
    J[91] -= dqdT; /* dwdot[H]/dT */
    J[93] -= dqdT; /* dwdot[O2]/dT */
    J[96] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 5: 2 H + M <=> H2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture - sc[0] - sc[5];
    /* forward */
    phi_f = sc[1]*sc[1];
    k_f = prefactor_units[4] * fwd_A[4]
                * exp(fwd_beta[4] * tc[0] - activation_units[4] * fwd_Ea[4] * invT);
    dlnkfdT = fwd_beta[4] * invT + activation_units[4] * fwd_Ea[4] * invT2;
    /* reverse */
    phi_r = sc[0];
    Kc = refCinv * exp(-g_RT[0] + 2*g_RT[1]);
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
        J[10] += dqdci;               /* dwdot[H2]/d[H] */
        J[11] += -2 * dqdci;          /* dwdot[H]/d[H] */
        /* d()/d[H2O] */
        dqdci =  -q_nocor;
        J[50] += dqdci;               /* dwdot[H2]/d[H2O] */
        J[51] += -2 * dqdci;          /* dwdot[H]/d[H2O] */
    }
    else {
        dqdc[0] =  - k_r;
        dqdc[1] =  q_nocor + k_f*2*sc[1];
        dqdc[2] =  q_nocor;
        dqdc[3] =  q_nocor;
        dqdc[4] =  q_nocor;
        dqdc[6] =  q_nocor;
        dqdc[7] =  q_nocor;
        dqdc[8] =  q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+0] += dqdc[k];
            J[10*k+1] += -2 * dqdc[k];
        }
    }
    J[90] += dqdT; /* dwdot[H2]/dT */
    J[91] += -2 * dqdT; /* dwdot[H]/dT */

    /*reaction 6: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture - 0.27000000000000002*sc[0] + 2.6499999999999999*sc[5];
    /* forward */
    phi_f = sc[1]*sc[4];
    k_f = prefactor_units[5] * fwd_A[5]
                * exp(fwd_beta[5] * tc[0] - activation_units[5] * fwd_Ea[5] * invT);
    dlnkfdT = fwd_beta[5] * invT + activation_units[5] * fwd_Ea[5] * invT2;
    /* reverse */
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[1] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[4]) + (h_RT[5]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci =  -0.27000000000000002*q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[4] -= dqdci;                /* dwdot[OH]/d[H2] */
        J[5] += dqdci;                /* dwdot[H2O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[4];
        J[11] -= dqdci;               /* dwdot[H]/d[H] */
        J[14] -= dqdci;               /* dwdot[OH]/d[H] */
        J[15] += dqdci;               /* dwdot[H2O]/d[H] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[1];
        J[41] -= dqdci;               /* dwdot[H]/d[OH] */
        J[44] -= dqdci;               /* dwdot[OH]/d[OH] */
        J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
        /* d()/d[H2O] */
        dqdci =  2.6499999999999999*q_nocor - k_r;
        J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
        J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    }
    else {
        dqdc[0] =  0.72999999999999998*q_nocor;
        dqdc[1] =  q_nocor + k_f*sc[4];
        dqdc[2] =  q_nocor;
        dqdc[3] =  q_nocor;
        dqdc[4] =  q_nocor + k_f*sc[1];
        dqdc[5] =  3.6499999999999999*q_nocor - k_r;
        dqdc[6] =  q_nocor;
        dqdc[7] =  q_nocor;
        dqdc[8] =  q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+1] -= dqdc[k];
            J[10*k+4] -= dqdc[k];
            J[10*k+5] += dqdc[k];
        }
    }
    J[91] -= dqdT; /* dwdot[H]/dT */
    J[94] -= dqdT; /* dwdot[OH]/dT */
    J[95] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 7: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = prefactor_units[6] * fwd_A[6]
                * exp(fwd_beta[6] * tc[0] - activation_units[6] * fwd_Ea[6] * invT);
    dlnkfdT = fwd_beta[6] * invT + activation_units[6] * fwd_Ea[6] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[4];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[1] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[2] -= dqdci;                /* dwdot[O]/d[H2] */
    J[4] += dqdci;                /* dwdot[OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4];
    J[10] -= dqdci;               /* dwdot[H2]/d[H] */
    J[11] += dqdci;               /* dwdot[H]/d[H] */
    J[12] -= dqdci;               /* dwdot[O]/d[H] */
    J[14] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[0];
    J[20] -= dqdci;               /* dwdot[H2]/d[O] */
    J[21] += dqdci;               /* dwdot[H]/d[O] */
    J[22] -= dqdci;               /* dwdot[O]/d[O] */
    J[24] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[40] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[41] += dqdci;               /* dwdot[H]/d[OH] */
    J[42] -= dqdci;               /* dwdot[O]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[90] -= dqdT;                /* dwdot[H2]/dT */
    J[91] += dqdT;                /* dwdot[H]/dT */
    J[92] -= dqdT;                /* dwdot[O]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */

    /*reaction 8: O + HO2 <=> OH + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = prefactor_units[7] * fwd_A[7]
                * exp(fwd_beta[7] * tc[0] - activation_units[7] * fwd_Ea[7] * invT);
    dlnkfdT = fwd_beta[7] * invT + activation_units[7] * fwd_Ea[7] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[4];
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[3] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[22] -= dqdci;               /* dwdot[O]/d[O] */
    J[23] += dqdci;               /* dwdot[O2]/d[O] */
    J[24] += dqdci;               /* dwdot[OH]/d[O] */
    J[26] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[4];
    J[32] -= dqdci;               /* dwdot[O]/d[O2] */
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[34] += dqdci;               /* dwdot[OH]/d[O2] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[3];
    J[42] -= dqdci;               /* dwdot[O]/d[OH] */
    J[43] += dqdci;               /* dwdot[O2]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    J[46] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[62] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[64] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[92] -= dqdT;                /* dwdot[O]/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = prefactor_units[8] * fwd_A[8]
                * exp(fwd_beta[8] * tc[0] - activation_units[8] * fwd_Ea[8] * invT);
    dlnkfdT = fwd_beta[8] * invT + activation_units[8] * fwd_Ea[8] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[6];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[4] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[22] -= dqdci;               /* dwdot[O]/d[O] */
    J[24] += dqdci;               /* dwdot[OH]/d[O] */
    J[26] += dqdci;               /* dwdot[HO2]/d[O] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[6];
    J[42] -= dqdci;               /* dwdot[O]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[62] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[64] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[72] -= dqdci;               /* dwdot[O]/d[H2O2] */
    J[74] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[92] -= dqdT;                /* dwdot[O]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 10: H + 2 O2 <=> HO2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[3];
    k_f = prefactor_units[9] * fwd_A[9]
                * exp(fwd_beta[9] * tc[0] - activation_units[9] * fwd_Ea[9] * invT);
    dlnkfdT = fwd_beta[9] * invT + activation_units[9] * fwd_Ea[9] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[6];
    Kc = refCinv * exp(g_RT[1] + 2*g_RT[3] - g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + 2*h_RT[3]) + (h_RT[3] + h_RT[6]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[3];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[13] -= dqdci;               /* dwdot[O2]/d[H] */
    J[16] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*2*sc[3] - k_r*sc[6];
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[3];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[93] -= dqdT;                /* dwdot[O2]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = prefactor_units[10] * fwd_A[10]
                * exp(fwd_beta[10] * tc[0] - activation_units[10] * fwd_Ea[10] * invT);
    dlnkfdT = fwd_beta[10] * invT + activation_units[10] * fwd_Ea[10] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] + g_RT[5] - g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[5]) + (h_RT[5] + h_RT[6]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[5];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[13] -= dqdci;               /* dwdot[O2]/d[H] */
    J[16] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[5];
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[53] -= dqdci;               /* dwdot[O2]/d[H2O] */
    J[56] += dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[93] -= dqdT;                /* dwdot[O2]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[8];
    k_f = prefactor_units[11] * fwd_A[11]
                * exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
    dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[8];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[8] - g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[8]) + (h_RT[6] + h_RT[8]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[8];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[13] -= dqdci;               /* dwdot[O2]/d[H] */
    J[16] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[8];
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[8];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[81] -= dqdci;               /* dwdot[H]/d[N2] */
    J[83] -= dqdci;               /* dwdot[O2]/d[N2] */
    J[86] += dqdci;               /* dwdot[HO2]/d[N2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[93] -= dqdT;                /* dwdot[O2]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */

    /*reaction 13: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[12] * fwd_A[12]
                * exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
    dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[4];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[2] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[12] += dqdci;               /* dwdot[O]/d[H] */
    J[13] -= dqdci;               /* dwdot[O2]/d[H] */
    J[14] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[4];
    J[21] -= dqdci;               /* dwdot[H]/d[O] */
    J[22] += dqdci;               /* dwdot[O]/d[O] */
    J[23] -= dqdci;               /* dwdot[O2]/d[O] */
    J[24] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1];
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[32] += dqdci;               /* dwdot[O]/d[O2] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[34] += dqdci;               /* dwdot[OH]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[41] -= dqdci;               /* dwdot[H]/d[OH] */
    J[42] += dqdci;               /* dwdot[O]/d[OH] */
    J[43] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[92] += dqdT;                /* dwdot[O]/dT */
    J[93] -= dqdT;                /* dwdot[O2]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */

    /*reaction 14: 2 H + H2 <=> 2 H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[1]*sc[1];
    k_f = prefactor_units[13] * fwd_A[13]
                * exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
    dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[0];
    Kc = refCinv * exp(g_RT[0] - 2*g_RT[0] + 2*g_RT[1]);
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
    J[10] += dqdci;               /* dwdot[H2]/d[H] */
    J[11] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] += -2 * dqdT;           /* dwdot[H]/dT */

    /*reaction 15: 2 H + H2O <=> H2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[1]*sc[5];
    k_f = prefactor_units[14] * fwd_A[14]
                * exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
    dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[5];
    Kc = refCinv * exp(-g_RT[0] + 2*g_RT[1] + g_RT[5] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[1] + h_RT[5]) + (h_RT[0] + h_RT[5]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[5];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*2*sc[1]*sc[5];
    J[10] += dqdci;               /* dwdot[H2]/d[H] */
    J[11] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[1]*sc[1] - k_r*sc[0];
    J[50] += dqdci;               /* dwdot[H2]/d[H2O] */
    J[51] += -2 * dqdci;          /* dwdot[H]/d[H2O] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] += -2 * dqdT;           /* dwdot[H]/dT */

    /*reaction 16: H + HO2 <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = exp(g_RT[1] - g_RT[2] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* O */
    wdot[5] += q; /* H2O */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[12] += dqdci;               /* dwdot[O]/d[H] */
    J[15] += dqdci;               /* dwdot[H2O]/d[H] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[5];
    J[21] -= dqdci;               /* dwdot[H]/d[O] */
    J[22] += dqdci;               /* dwdot[O]/d[O] */
    J[25] += dqdci;               /* dwdot[H2O]/d[O] */
    J[26] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[52] += dqdci;               /* dwdot[O]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[62] += dqdci;               /* dwdot[O]/d[HO2] */
    J[65] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[92] += dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 17: H + HO2 <=> O2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[3] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[3] += dqdci;                /* dwdot[O2]/d[H2] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[10] += dqdci;               /* dwdot[H2]/d[H] */
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[13] += dqdci;               /* dwdot[O2]/d[H] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[0];
    J[30] += dqdci;               /* dwdot[H2]/d[O2] */
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[60] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 18: H + HO2 <=> 2 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[4];
    Kc = exp(g_RT[1] - 2*g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (2*h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += 2 * q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[14] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2*sc[4];
    J[41] -= dqdci;               /* dwdot[H]/d[OH] */
    J[44] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[46] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[64] += 2 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[94] += 2 * dqdT;            /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[6];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[7]) + (h_RT[0] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[6];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[6] += dqdci;                /* dwdot[HO2]/d[H2] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[10] += dqdci;               /* dwdot[H2]/d[H] */
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[16] += dqdci;               /* dwdot[HO2]/d[H] */
    J[17] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[0];
    J[60] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[70] += dqdci;               /* dwdot[H2]/d[H2O2] */
    J[71] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 20: H + H2O2 <=> OH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[1] - g_RT[4] - g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[7]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[14] += dqdci;               /* dwdot[OH]/d[H] */
    J[15] += dqdci;               /* dwdot[H2O]/d[H] */
    J[17] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[5];
    J[41] -= dqdci;               /* dwdot[H]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[4];
    J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[54] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[71] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[74] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 21: OH + H2 <=> H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[4];
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[4]) + (h_RT[1] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[4];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[4] -= dqdci;                /* dwdot[OH]/d[H2] */
    J[5] += dqdci;                /* dwdot[H2O]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5];
    J[10] -= dqdci;               /* dwdot[H2]/d[H] */
    J[11] += dqdci;               /* dwdot[H]/d[H] */
    J[14] -= dqdci;               /* dwdot[OH]/d[H] */
    J[15] += dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[0];
    J[40] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[41] += dqdci;               /* dwdot[H]/d[OH] */
    J[44] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[50] -= dqdci;               /* dwdot[H2]/d[H2O] */
    J[51] += dqdci;               /* dwdot[H]/d[H2O] */
    J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[90] -= dqdT;                /* dwdot[H2]/dT */
    J[91] += dqdT;                /* dwdot[H]/dT */
    J[94] -= dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */

    /*reaction 22: 2 OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[4];
    k_f = prefactor_units[21] * fwd_A[21]
                * exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
    dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = exp(-g_RT[2] + 2*g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[4]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[4] -= 2 * q; /* OH */
    wdot[5] += q; /* H2O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[5];
    J[22] += dqdci;               /* dwdot[O]/d[O] */
    J[24] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[25] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2*sc[4];
    J[42] += dqdci;               /* dwdot[O]/d[OH] */
    J[44] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[52] += dqdci;               /* dwdot[O]/d[H2O] */
    J[54] += -2 * dqdci;          /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[92] += dqdT;                /* dwdot[O]/dT */
    J[94] += -2 * dqdT;           /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[6];
    k_f = prefactor_units[22] * fwd_A[22]
                * exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
    dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(-g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[6]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[5];
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[34] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[35] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[43] += dqdci;               /* dwdot[O2]/d[OH] */
    J[44] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[46] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[53] += dqdci;               /* dwdot[O2]/d[H2O] */
    J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[64] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[65] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[94] -= dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = prefactor_units[23] * fwd_A[23]
                * exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
    dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[5] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[44] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[46] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[56] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[64] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[65] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[4];
    J[74] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[94] -= dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = prefactor_units[24] * fwd_A[24]
                * exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
    dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[5] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[44] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[46] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[56] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[64] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[65] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[4];
    J[74] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[94] -= dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 26: 2 HO2 <=> O2 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[6];
    k_f = prefactor_units[25] * fwd_A[25]
                * exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
    dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[7];
    Kc = exp(-g_RT[3] + 2*g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[6]) + (h_RT[3] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[37] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2*sc[6];
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[67] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[3];
    J[73] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[96] += -2 * dqdT;           /* dwdot[HO2]/dT */
    J[97] += dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 27: 2 HO2 <=> O2 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[6];
    k_f = prefactor_units[26] * fwd_A[26]
                * exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
    dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[7];
    Kc = exp(-g_RT[3] + 2*g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[6]) + (h_RT[3] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[37] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2*sc[6];
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[67] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[3];
    J[73] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[96] += -2 * dqdT;           /* dwdot[HO2]/dT */
    J[97] += dqdT;                /* dwdot[H2O2]/dT */

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
        /*species 2: O */
        species[2] =
            -3.27931884e-03
            +1.32861279e-05 * tc[1]
            -1.83841987e-08 * tc[2]
            +8.45063884e-12 * tc[3];
        /*species 3: O2 */
        species[3] =
            -2.99673416e-03
            +1.96946040e-05 * tc[1]
            -2.90438853e-08 * tc[2]
            +1.29749135e-11 * tc[3];
        /*species 4: OH */
        species[4] =
            -2.40131752e-03
            +9.23587682e-06 * tc[1]
            -1.16434000e-08 * tc[2]
            +5.45645880e-12 * tc[3];
        /*species 5: H2O */
        species[5] =
            -2.03643410e-03
            +1.30408042e-05 * tc[1]
            -1.64639119e-08 * tc[2]
            +7.08791268e-12 * tc[3];
        /*species 6: HO2 */
        species[6] =
            -4.74912051e-03
            +4.23165782e-05 * tc[1]
            -7.28291682e-08 * tc[2]
            +3.71690050e-11 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            -5.42822417e-04
            +3.34671402e-05 * tc[1]
            -6.47312439e-08 * tc[2]
            +3.44981745e-11 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.40824040e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77941600e-12 * tc[3];
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
        /*species 2: O */
        species[2] =
            -8.59741137e-05
            +8.38969178e-08 * tc[1]
            -3.00533397e-11 * tc[2]
            +4.91334764e-15 * tc[3];
        /*species 3: O2 */
        species[3] =
            +1.48308754e-03
            -1.51593334e-06 * tc[1]
            +6.28411665e-10 * tc[2]
            -8.66871176e-14 * tc[3];
        /*species 4: OH */
        species[4] =
            +5.48429716e-04
            +2.53010456e-07 * tc[1]
            -2.63838467e-10 * tc[2]
            +4.69649504e-14 * tc[3];
        /*species 5: H2O */
        species[5] =
            +2.17691804e-03
            -3.28145036e-07 * tc[1]
            -2.91125961e-10 * tc[2]
            +6.72803968e-14 * tc[3];
        /*species 6: HO2 */
        species[6] =
            +2.23982013e-03
            -1.26731630e-06 * tc[1]
            +3.42739110e-10 * tc[2]
            -4.31634140e-14 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +4.90831694e-03
            -3.80278450e-06 * tc[1]
            +1.11355796e-09 * tc[2]
            -1.15163322e-13 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.48797680e-03
            -1.13695200e-06 * tc[1]
            +3.02911140e-10 * tc[2]
            -2.70134040e-14 * tc[3];
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

    double q_f[27], q_r[27];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 27; ++i) {
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

    /*reaction 1: 2 OH (+M) <=> H2O2 (+M) */
    kc[0] = 1.0 / (refC) * exp((2 * g_RT[4]) - (g_RT[7]));

    /*reaction 2: 2 O + M <=> O2 + M */
    kc[1] = 1.0 / (refC) * exp((2 * g_RT[2]) - (g_RT[3]));

    /*reaction 3: O + H + M <=> OH + M */
    kc[2] = 1.0 / (refC) * exp((g_RT[2] + g_RT[1]) - (g_RT[4]));

    /*reaction 4: H + O2 + M <=> HO2 + M */
    kc[3] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));

    /*reaction 5: 2 H + M <=> H2 + M */
    kc[4] = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));

    /*reaction 6: H + OH + M <=> H2O + M */
    kc[5] = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));

    /*reaction 7: O + H2 <=> H + OH */
    kc[6] = exp((g_RT[2] + g_RT[0]) - (g_RT[1] + g_RT[4]));

    /*reaction 8: O + HO2 <=> OH + O2 */
    kc[7] = exp((g_RT[2] + g_RT[6]) - (g_RT[4] + g_RT[3]));

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    kc[8] = exp((g_RT[2] + g_RT[7]) - (g_RT[4] + g_RT[6]));

    /*reaction 10: H + 2 O2 <=> HO2 + O2 */
    kc[9] = 1.0 / (refC) * exp((g_RT[1] + 2 * g_RT[3]) - (g_RT[6] + g_RT[3]));

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    kc[10] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    kc[11] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[8]) - (g_RT[6] + g_RT[8]));

    /*reaction 13: H + O2 <=> O + OH */
    kc[12] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));

    /*reaction 14: 2 H + H2 <=> 2 H2 */
    kc[13] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));

    /*reaction 15: 2 H + H2O <=> H2 + H2O */
    kc[14] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));

    /*reaction 16: H + HO2 <=> O + H2O */
    kc[15] = exp((g_RT[1] + g_RT[6]) - (g_RT[2] + g_RT[5]));

    /*reaction 17: H + HO2 <=> O2 + H2 */
    kc[16] = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));

    /*reaction 18: H + HO2 <=> 2 OH */
    kc[17] = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    kc[18] = exp((g_RT[1] + g_RT[7]) - (g_RT[6] + g_RT[0]));

    /*reaction 20: H + H2O2 <=> OH + H2O */
    kc[19] = exp((g_RT[1] + g_RT[7]) - (g_RT[4] + g_RT[5]));

    /*reaction 21: OH + H2 <=> H + H2O */
    kc[20] = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));

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
        /*species 8: N2 */
        species[8] =
            -1.020899900000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041202000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427000000000e-13 * tc[4];
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
        /*species 8: N2 */
        species[8] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439884000000000e-04 * tc[1]
            +9.474600000000001e-08 * tc[2]
            -8.414198333333333e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
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
        /*species 8: N2 */
        species[8] =
            -1.02089990e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
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
        /*species 8: N2 */
        species[8] =
            -9.22797700e+02 * invT
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
        /*species 8: N2 */
        species[8] =
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
        /*species 8: N2 */
        species[8] =
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
        /*species 8: N2 */
        species[8] =
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
        /*species 8: N2 */
        species[8] =
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
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
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
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
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
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
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
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
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
        /*species 8: N2 */
        species[8] =
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
        /*species 8: N2 */
        species[8] =
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
void molecularWeight(double * restrict wt)
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


/*save atomic weights into array */
void atomicWeight(double * restrict awt)
{
    awt[0] = 15.999400; /*O */
    awt[1] = 1.007970; /*H */
    awt[2] = 14.006700; /*N */

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
  WT[           1] =   0.1007969975471497E+01;
  WT[           2] =   0.1599940013885498E+02;
  WT[           3] =   0.3199880027770996E+02;
  WT[           4] =   0.1700737011432648E+02;
  WT[           5] =   0.1801534008979797E+02;
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
  EPS[           1] =   0.1450000000000000E+03;
  EPS[           2] =   0.8000000000000000E+02;
  EPS[           3] =   0.1074000000000000E+03;
  EPS[           4] =   0.8000000000000000E+02;
  EPS[           5] =   0.5724000000000000E+03;
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
  SIG[           1] =   0.2050000000000000E+01;
  SIG[           2] =   0.2750000000000000E+01;
  SIG[           3] =   0.3458000000000000E+01;
  SIG[           4] =   0.2750000000000000E+01;
  SIG[           5] =   0.2605000000000000E+01;
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
  DIP[           2] =   0.0000000000000000E+00;
  DIP[           3] =   0.0000000000000000E+00;
  DIP[           4] =   0.0000000000000000E+00;
  DIP[           5] =   0.1844000000000000E+01;
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
  POL[           1] =   0.0000000000000000E+00;
  POL[           2] =   0.0000000000000000E+00;
  POL[           3] =   0.1600000000000000E+01;
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
  ZROT[           1] =   0.0000000000000000E+00;
  ZROT[           2] =   0.0000000000000000E+00;
  ZROT[           3] =   0.3800000000000000E+01;
  ZROT[           4] =   0.0000000000000000E+00;
  ZROT[           5] =   0.4000000000000000E+01;
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
  NLIN[           1] =            0;
  NLIN[           2] =            0;
  NLIN[           3] =            1;
  NLIN[           4] =            1;
  NLIN[           5] =            2;
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
  COFLAM[           0] =   0.1158096158594203E+02;
  COFLAM[           1] =  -0.1520398050330763E+01;
  COFLAM[           2] =   0.2722793392212857E+00;
  COFLAM[           3] =  -0.1031257389909708E-01;
  COFLAM[           4] =  -0.3253788579319486E+00;
  COFLAM[           5] =   0.3416397862210875E+01;
  COFLAM[           6] =  -0.3631104489846262E+00;
  COFLAM[           7] =   0.1585986202702918E-01;
  COFLAM[           8] =   0.1969666526185763E+01;
  COFLAM[           9] =   0.1801396425547300E+01;
  COFLAM[          10] =  -0.1549102264837655E+00;
  COFLAM[          11] =   0.6908759721669425E-02;
  COFLAM[          12] =  -0.3013713897782378E+01;
  COFLAM[          13] =   0.3375529164503225E+01;
  COFLAM[          14] =  -0.3433040143356826E+00;
  COFLAM[          15] =   0.1509986872028326E-01;
  COFLAM[          16] =   0.1534175890276220E+02;
  COFLAM[          17] =  -0.3776583656836550E+01;
  COFLAM[          18] =   0.6131005943656342E+00;
  COFLAM[          19] =  -0.2721052763006557E-01;
  COFLAM[          20] =   0.2277162937237345E+02;
  COFLAM[          21] =  -0.8704976392893332E+01;
  COFLAM[          22] =   0.1490812066077175E+01;
  COFLAM[          23] =  -0.7406350595204310E-01;
  COFLAM[          24] =   0.5546401577805573E+00;
  COFLAM[          25] =   0.1591057931808813E+01;
  COFLAM[          26] =  -0.5282455808284543E-01;
  COFLAM[          27] =   0.4072391521895438E-03;
  COFLAM[          28] =   0.6252374647013528E+00;
  COFLAM[          29] =   0.1431934364158212E+01;
  COFLAM[          30] =   0.1750712355480039E-02;
  COFLAM[          31] =  -0.3554756889670176E-02;
  COFLAM[          32] =   0.1154286034252058E+02;
  COFLAM[          33] =  -0.2911509179774523E+01;
  COFLAM[          34] =   0.5546559973428081E+00;
  COFLAM[          35] =  -0.2750092774174382E-01;
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
  COFETA[           8] =  -0.1481563591580279E+02;
  COFETA[           9] =   0.1801396425547327E+01;
  COFETA[          10] =  -0.1549102264837694E+00;
  COFETA[          11] =   0.6908759721669611E-02;
  COFETA[          12] =  -0.1681080110872507E+02;
  COFETA[          13] =   0.2522528725237932E+01;
  COFETA[          14] =  -0.2490712798240462E+00;
  COFETA[          15] =   0.1100615806447165E-01;
  COFETA[          16] =  -0.1478508813778873E+02;
  COFETA[          17] =   0.1801396425547305E+01;
  COFETA[          18] =  -0.1549102264837663E+00;
  COFETA[          19] =   0.6908759721669463E-02;
  COFETA[          20] =  -0.1187800764707307E+02;
  COFETA[          21] =  -0.7882519505482808E+00;
  COFETA[          22] =   0.3341408170058896E+00;
  COFETA[          23] =  -0.1986366361647418E-01;
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
  COFD[           4] =  -0.1168685516750758E+02;
  COFD[           5] =   0.2883726644585007E+01;
  COFD[           6] =  -0.1637777890869061E+00;
  COFD[           7] =   0.7265872474830623E-02;
  COFD[           8] =  -0.1060109579664357E+02;
  COFD[           9] =   0.2157129122185176E+01;
  COFD[          10] =  -0.6524729611183182E-01;
  COFD[          11] =   0.2809597810849289E-02;
  COFD[          12] =  -0.1228690337278460E+02;
  COFD[          13] =   0.2739817607547901E+01;
  COFD[          14] =  -0.1455907617836413E+00;
  COFD[          15] =   0.6496697467578658E-02;
  COFD[          16] =  -0.1060412419807472E+02;
  COFD[          17] =   0.2157123706881953E+01;
  COFD[          18] =  -0.6524690383610432E-01;
  COFD[          19] =   0.2809593789707914E-02;
  COFD[          20] =  -0.1787213075631945E+02;
  COFD[          21] =   0.4905112324365878E+01;
  COFD[          22] =  -0.4173872143518051E+00;
  COFD[          23] =   0.1788845489573668E-01;
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
  COFD[          36] =  -0.1168685516750758E+02;
  COFD[          37] =   0.2883726644585007E+01;
  COFD[          38] =  -0.1637777890869061E+00;
  COFD[          39] =   0.7265872474830623E-02;
  COFD[          40] =  -0.1476537827842222E+02;
  COFD[          41] =   0.4195448036781197E+01;
  COFD[          42] =  -0.3279615676730872E+00;
  COFD[          43] =   0.1412388797372668E-01;
  COFD[          44] =  -0.1500301598770682E+02;
  COFD[          45] =   0.4131920691035950E+01;
  COFD[          46] =  -0.3275337549968449E+00;
  COFD[          47] =   0.1442760802547180E-01;
  COFD[          48] =  -0.1719486110604737E+02;
  COFD[          49] =   0.4862600572193149E+01;
  COFD[          50] =  -0.4215909171603991E+00;
  COFD[          51] =   0.1846373197001472E-01;
  COFD[          52] =  -0.1502025959280452E+02;
  COFD[          53] =   0.4138327950805945E+01;
  COFD[          54] =  -0.3284006010219734E+00;
  COFD[          55] =   0.1446660099762999E-01;
  COFD[          56] =  -0.1695183404426207E+02;
  COFD[          57] =   0.4416203740445888E+01;
  COFD[          58] =  -0.3114890192543383E+00;
  COFD[          59] =   0.1159663710619628E-01;
  COFD[          60] =  -0.1720445512239349E+02;
  COFD[          61] =   0.4866330968714804E+01;
  COFD[          62] =  -0.4220901924876364E+00;
  COFD[          63] =   0.1848595807162496E-01;
  COFD[          64] =  -0.1721362788263434E+02;
  COFD[          65] =   0.4869900420159440E+01;
  COFD[          66] =  -0.4225679381202458E+00;
  COFD[          67] =   0.1850722613298993E-01;
  COFD[          68] =  -0.1676897875889624E+02;
  COFD[          69] =   0.4677755283143664E+01;
  COFD[          70] =  -0.3974479445519867E+00;
  COFD[          71] =   0.1741111717763365E-01;
  COFD[          72] =  -0.1060109579664357E+02;
  COFD[          73] =   0.2157129122185176E+01;
  COFD[          74] =  -0.6524729611183182E-01;
  COFD[          75] =   0.2809597810849289E-02;
  COFD[          76] =  -0.1500301598770682E+02;
  COFD[          77] =   0.4131920691035950E+01;
  COFD[          78] =  -0.3275337549968449E+00;
  COFD[          79] =   0.1442760802547180E-01;
  COFD[          80] =  -0.1329503778398344E+02;
  COFD[          81] =   0.2938989816956702E+01;
  COFD[          82] =  -0.1706012465018490E+00;
  COFD[          83] =   0.7547358484820122E-02;
  COFD[          84] =  -0.1474363832488901E+02;
  COFD[          85] =   0.3350090901251298E+01;
  COFD[          86] =  -0.2245382243795846E+00;
  COFD[          87] =   0.9906528312288899E-02;
  COFD[          88] =  -0.1331106228772027E+02;
  COFD[          89] =   0.2939408772392007E+01;
  COFD[          90] =  -0.1706589698617190E+00;
  COFD[          91] =   0.7549998035548285E-02;
  COFD[          92] =  -0.1894546586601941E+02;
  COFD[          93] =   0.4935720676293722E+01;
  COFD[          94] =  -0.4114057447467022E+00;
  COFD[          95] =   0.1723315986098057E-01;
  COFD[          96] =  -0.1476396639560408E+02;
  COFD[          97] =   0.3356474034011327E+01;
  COFD[          98] =  -0.2254105976116209E+00;
  COFD[          99] =   0.9946130954817655E-02;
  COFD[         100] =  -0.1478376315874751E+02;
  COFD[         101] =   0.3362741024796286E+01;
  COFD[         102] =  -0.2262670743508222E+00;
  COFD[         103] =   0.9985011202484037E-02;
  COFD[         104] =  -0.1442662508885140E+02;
  COFD[         105] =   0.3216422191096490E+01;
  COFD[         106] =  -0.2069561315576079E+00;
  COFD[         107] =   0.9135288598631700E-02;
  COFD[         108] =  -0.1228690337278460E+02;
  COFD[         109] =   0.2739817607547901E+01;
  COFD[         110] =  -0.1455907617836413E+00;
  COFD[         111] =   0.6496697467578658E-02;
  COFD[         112] =  -0.1719486110604737E+02;
  COFD[         113] =   0.4862600572193149E+01;
  COFD[         114] =  -0.4215909171603991E+00;
  COFD[         115] =   0.1846373197001472E-01;
  COFD[         116] =  -0.1474363832488901E+02;
  COFD[         117] =   0.3350090901251298E+01;
  COFD[         118] =  -0.2245382243795846E+00;
  COFD[         119] =   0.9906528312288899E-02;
  COFD[         120] =  -0.1579169675646239E+02;
  COFD[         121] =   0.3572143437285479E+01;
  COFD[         122] =  -0.2518469828462104E+00;
  COFD[         123] =   0.1102533331592793E-01;
  COFD[         124] =  -0.1473449002110577E+02;
  COFD[         125] =   0.3337793989426398E+01;
  COFD[         126] =  -0.2228575541864579E+00;
  COFD[         127] =   0.9830229599517820E-02;
  COFD[         128] =  -0.2036133619470493E+02;
  COFD[         129] =   0.5195864695910879E+01;
  COFD[         130] =  -0.4301216528920454E+00;
  COFD[         131] =   0.1744936825492251E-01;
  COFD[         132] =  -0.1579979030842146E+02;
  COFD[         133] =   0.3572309323030401E+01;
  COFD[         134] =  -0.2518694694768392E+00;
  COFD[         135] =   0.1102634618303224E-01;
  COFD[         136] =  -0.1580828869487550E+02;
  COFD[         137] =   0.3572786632028933E+01;
  COFD[         138] =  -0.2519341709914386E+00;
  COFD[         139] =   0.1102926053643803E-01;
  COFD[         140] =  -0.1544409203507008E+02;
  COFD[         141] =   0.3434913447661099E+01;
  COFD[         142] =  -0.2339977102148624E+00;
  COFD[         143] =   0.1025033359000777E-01;
  COFD[         144] =  -0.1060412419807472E+02;
  COFD[         145] =   0.2157123706881953E+01;
  COFD[         146] =  -0.6524690383610432E-01;
  COFD[         147] =   0.2809593789707914E-02;
  COFD[         148] =  -0.1502025959280452E+02;
  COFD[         149] =   0.4138327950805945E+01;
  COFD[         150] =  -0.3284006010219734E+00;
  COFD[         151] =   0.1446660099762999E-01;
  COFD[         152] =  -0.1331106228772027E+02;
  COFD[         153] =   0.2939408772392007E+01;
  COFD[         154] =  -0.1706589698617190E+00;
  COFD[         155] =   0.7549998035548285E-02;
  COFD[         156] =  -0.1473449002110577E+02;
  COFD[         157] =   0.3337793989426398E+01;
  COFD[         158] =  -0.2228575541864579E+00;
  COFD[         159] =   0.9830229599517820E-02;
  COFD[         160] =  -0.1332558556199739E+02;
  COFD[         161] =   0.2938989816956671E+01;
  COFD[         162] =  -0.1706012465018444E+00;
  COFD[         163] =   0.7547358484819904E-02;
  COFD[         164] =  -0.1896004007077217E+02;
  COFD[         165] =   0.4935377774711531E+01;
  COFD[         166] =  -0.4113926173034602E+00;
  COFD[         167] =   0.1723402811988101E-01;
  COFD[         168] =  -0.1475457475365462E+02;
  COFD[         169] =   0.3343986593397501E+01;
  COFD[         170] =  -0.2237039347662408E+00;
  COFD[         171] =   0.9868653787965939E-02;
  COFD[         172] =  -0.1477418610290300E+02;
  COFD[         173] =   0.3350090901251285E+01;
  COFD[         174] =  -0.2245382243795827E+00;
  COFD[         175] =   0.9906528312288795E-02;
  COFD[         176] =  -0.1442063092189782E+02;
  COFD[         177] =   0.3205844724257626E+01;
  COFD[         178] =  -0.2055144159347264E+00;
  COFD[         179] =   0.9070008933822527E-02;
  COFD[         180] =  -0.1787213075631945E+02;
  COFD[         181] =   0.4905112324365878E+01;
  COFD[         182] =  -0.4173872143518051E+00;
  COFD[         183] =   0.1788845489573668E-01;
  COFD[         184] =  -0.1695183404426207E+02;
  COFD[         185] =   0.4416203740445888E+01;
  COFD[         186] =  -0.3114890192543383E+00;
  COFD[         187] =   0.1159663710619628E-01;
  COFD[         188] =  -0.1894546586601941E+02;
  COFD[         189] =   0.4935720676293722E+01;
  COFD[         190] =  -0.4114057447467022E+00;
  COFD[         191] =   0.1723315986098057E-01;
  COFD[         192] =  -0.2036133619470493E+02;
  COFD[         193] =   0.5195864695910879E+01;
  COFD[         194] =  -0.4301216528920454E+00;
  COFD[         195] =   0.1744936825492251E-01;
  COFD[         196] =  -0.1896004007077217E+02;
  COFD[         197] =   0.4935377774711531E+01;
  COFD[         198] =  -0.4113926173034602E+00;
  COFD[         199] =   0.1723402811988101E-01;
  COFD[         200] =  -0.1301206458003408E+02;
  COFD[         201] =   0.1429168452568504E+01;
  COFD[         202] =   0.1661557715508861E+00;
  COFD[         203] =  -0.1214321823404827E-01;
  COFD[         204] =  -0.1973436691750415E+02;
  COFD[         205] =   0.4993125184848788E+01;
  COFD[         206] =  -0.4088531920998837E+00;
  COFD[         207] =   0.1672325900844039E-01;
  COFD[         208] =  -0.1972379377029397E+02;
  COFD[         209] =   0.4985700637038792E+01;
  COFD[         210] =  -0.4077156220815392E+00;
  COFD[         211] =   0.1666668649763390E-01;
  COFD[         212] =  -0.2025975101746627E+02;
  COFD[         213] =   0.5176212790334279E+01;
  COFD[         214] =  -0.4308686680573706E+00;
  COFD[         215] =   0.1761066268522406E-01;
  COFD[         216] =  -0.1229066279603360E+02;
  COFD[         217] =   0.2741057082889793E+01;
  COFD[         218] =  -0.1457627022193994E+00;
  COFD[         219] =   0.6504615244946633E-02;
  COFD[         220] =  -0.1720445512239349E+02;
  COFD[         221] =   0.4866330968714804E+01;
  COFD[         222] =  -0.4220901924876364E+00;
  COFD[         223] =   0.1848595807162496E-01;
  COFD[         224] =  -0.1476396639560408E+02;
  COFD[         225] =   0.3356474034011327E+01;
  COFD[         226] =  -0.2254105976116209E+00;
  COFD[         227] =   0.9946130954817655E-02;
  COFD[         228] =  -0.1579979030842146E+02;
  COFD[         229] =   0.3572309323030401E+01;
  COFD[         230] =  -0.2518694694768392E+00;
  COFD[         231] =   0.1102634618303224E-01;
  COFD[         232] =  -0.1475457475365462E+02;
  COFD[         233] =   0.3343986593397501E+01;
  COFD[         234] =  -0.2237039347662408E+00;
  COFD[         235] =   0.9868653787965939E-02;
  COFD[         236] =  -0.1973436691750415E+02;
  COFD[         237] =   0.4993125184848788E+01;
  COFD[         238] =  -0.4088531920998837E+00;
  COFD[         239] =   0.1672325900844039E-01;
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
  COFD[         256] =  -0.1721362788263434E+02;
  COFD[         257] =   0.4869900420159440E+01;
  COFD[         258] =  -0.4225679381202458E+00;
  COFD[         259] =   0.1850722613298993E-01;
  COFD[         260] =  -0.1478376315874751E+02;
  COFD[         261] =   0.3362741024796286E+01;
  COFD[         262] =  -0.2262670743508222E+00;
  COFD[         263] =   0.9985011202484037E-02;
  COFD[         264] =  -0.1580828869487550E+02;
  COFD[         265] =   0.3572786632028933E+01;
  COFD[         266] =  -0.2519341709914386E+00;
  COFD[         267] =   0.1102926053643803E-01;
  COFD[         268] =  -0.1477418610290300E+02;
  COFD[         269] =   0.3350090901251285E+01;
  COFD[         270] =  -0.2245382243795827E+00;
  COFD[         271] =   0.9906528312288795E-02;
  COFD[         272] =  -0.1972379377029397E+02;
  COFD[         273] =   0.4985700637038792E+01;
  COFD[         274] =  -0.4077156220815392E+00;
  COFD[         275] =   0.1666668649763390E-01;
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
  COFD[         292] =  -0.1676897875889624E+02;
  COFD[         293] =   0.4677755283143664E+01;
  COFD[         294] =  -0.3974479445519867E+00;
  COFD[         295] =   0.1741111717763365E-01;
  COFD[         296] =  -0.1442662508885140E+02;
  COFD[         297] =   0.3216422191096490E+01;
  COFD[         298] =  -0.2069561315576079E+00;
  COFD[         299] =   0.9135288598631700E-02;
  COFD[         300] =  -0.1544409203507008E+02;
  COFD[         301] =   0.3434913447661099E+01;
  COFD[         302] =  -0.2339977102148624E+00;
  COFD[         303] =   0.1025033359000777E-01;
  COFD[         304] =  -0.1442063092189782E+02;
  COFD[         305] =   0.3205844724257626E+01;
  COFD[         306] =  -0.2055144159347264E+00;
  COFD[         307] =   0.9070008933822527E-02;
  COFD[         308] =  -0.2025975101746627E+02;
  COFD[         309] =   0.5176212790334279E+01;
  COFD[         310] =  -0.4308686680573706E+00;
  COFD[         311] =   0.1761066268522406E-01;
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
  KTDIF[           1] =            2;
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
  COFTD[           8] =   0.4155834523029710E+00;
  COFTD[           9] =   0.1097383911871691E-04;
  COFTD[          10] =  -0.3960224834087406E-08;
  COFTD[          11] =   0.1144145307552153E-11;
  COFTD[          12] =   0.4427392097103151E+00;
  COFTD[          13] =   0.7117708739452678E-04;
  COFTD[          14] =  -0.3847681435586731E-07;
  COFTD[          15] =   0.6863235673780089E-11;
  COFTD[          16] =   0.4219325599835782E+00;
  COFTD[          17] =   0.1114149277732201E-04;
  COFTD[          18] =  -0.4020727469049531E-08;
  COFTD[          19] =   0.1161625074178193E-11;
  COFTD[          20] =   0.6020321200392946E-01;
  COFTD[          21] =   0.5615615962686471E-03;
  COFTD[          22] =  -0.2553727790421687E-06;
  COFTD[          23] =   0.3633898304582181E-10;
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
  COFTD[          40] =   0.0000000000000000E+00;
  COFTD[          41] =   0.0000000000000000E+00;
  COFTD[          42] =   0.0000000000000000E+00;
  COFTD[          43] =   0.0000000000000000E+00;
  COFTD[          44] =   0.2700102625474920E+00;
  COFTD[          45] =   0.3615551214404962E-03;
  COFTD[          46] =  -0.1807447676772846E-06;
  COFTD[          47] =   0.2753212716087208E-10;
  COFTD[          48] =   0.2204829568680314E+00;
  COFTD[          49] =   0.4801643237618346E-03;
  COFTD[          50] =  -0.2329279629943256E-06;
  COFTD[          51] =   0.3464704623498651E-10;
  COFTD[          52] =   0.2720417770528547E+00;
  COFTD[          53] =   0.3642754049836650E-03;
  COFTD[          54] =  -0.1821046627192010E-06;
  COFTD[          55] =   0.2773927453061686E-10;
  COFTD[          56] =  -0.1418836567013430E+00;
  COFTD[          57] =   0.7665588601033144E-03;
  COFTD[          58] =  -0.3065500230258140E-06;
  COFTD[          59] =   0.4029595271283201E-10;
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
! GRI-Mech Version 2.11 11/3/95  1st release 9/6/95  CHEMKIN-II format
! RETAINING ONLY THE RELEVANT HYDROGEN-OXYGEN CHEMISTRY
ELEMENTS
O  H  N
END
SPECIES
H2      H       O       O2      OH      H2O     HO2     H2O2    
N2
END
!THERMO
! Insert GRI-Mech thermodynamics here or use in default file
!END
REACTIONS
2O+M<=>O2+M                              1.200E+17   -1.000        .00
H2/ 2.40/ H2O/15.40/ 
O+H+M<=>OH+M                             5.000E+17   -1.000        .00
H2/2.00/ H2O/6.00/
O+H2<=>H+OH                              5.000E+04    2.670    6290.00
O+HO2<=>OH+O2                            2.000E+13     .000        .00
O+H2O2<=>OH+HO2                          9.630E+06    2.000    4000.00
H+O2+M<=>HO2+M                           2.800E+18    -.860        .00
O2/ .00/ H2O/ .00/ N2/ .00/
H+2O2<=>HO2+O2                           3.000E+20   -1.720        .00
H+O2+H2O<=>HO2+H2O                       9.380E+18    -.760        .00
H+O2+N2<=>HO2+N2                         3.750E+20   -1.720        .00
H+O2<=>O+OH                              8.300E+13     .000   14413.00
2H+M<=>H2+M                              1.000E+18   -1.000        .00
H2/ .00/ H2O/ .00/
2H+H2<=>2H2                              9.000E+16    -.600        .00
2H+H2O<=>H2+H2O                          6.000E+19   -1.250        .00
H+OH+M<=>H2O+M                           2.200E+22   -2.000        .00
H2/ .73/ H2O/3.65/
H+HO2<=>O+H2O                            3.970E+12     .000     671.00
H+HO2<=>O2+H2                            2.800E+13     .000    1068.00
H+HO2<=>2OH                              1.340E+14     .000     635.00
H+H2O2<=>HO2+H2                          1.210E+07    2.000    5200.00
H+H2O2<=>OH+H2O                          1.000E+13     .000    3600.00
OH+H2<=>H+H2O                            2.160E+08    1.510    3430.00
2OH(+M)<=>H2O2(+M)                       7.400E+13    -.370        .00
     LOW  /  2.300E+18    -.900  -1700.00/
     TROE/   .7346   94.00  1756.00  5182.00 /
H2/2.00/ H2O/6.00/
2OH<=>O+H2O                              3.570E+04    2.400   -2110.00
OH+HO2<=>O2+H2O                          2.900E+13     .000    -500.00
OH+H2O2<=>HO2+H2O                        1.750E+12     .000     320.00
 DUPLICATE
OH+H2O2<=>HO2+H2O                        5.800E+14     .000    9560.00
 DUPLICATE
2HO2<=>O2+H2O2                           1.300E+11     .000   -1630.00
 DUPLICATE
2HO2<=>O2+H2O2                           4.200E+14     .000   12000.00
 DUPLICATE
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
