
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#if defined(BL_FORT_USE_UPPERCASE)
#define CKINDX CKINDX
#define CKINIT CKINIT
#define CKFINALIZE CKFINALIZE
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
#define CKFINALIZE ckfinalize
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
#define CKFINALIZE ckfinalize_
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
void CKFINALIZE();
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
void vcomp_wdot_1_50(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,
                double * restrict k_f_s, double * restrict Kc_s,
                double * restrict tc, double * restrict invT, double * restrict T);
void vcomp_wdot_51_100(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,
                double * restrict k_f_s, double * restrict Kc_s,
                double * restrict tc, double * restrict invT, double * restrict T);
void vcomp_wdot_101_104(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,
                double * restrict k_f_s, double * restrict Kc_s,
                double * restrict tc, double * restrict invT, double * restrict T);

/* Inverse molecular weights */
static const double imw[24] = {
    1.0 / 2.015940,  /*H2 */
    1.0 / 1.007970,  /*H */
    1.0 / 15.999400,  /*O */
    1.0 / 31.998800,  /*O2 */
    1.0 / 17.007370,  /*OH */
    1.0 / 18.015340,  /*H2O */
    1.0 / 33.006770,  /*HO2 */
    1.0 / 34.014740,  /*H2O2 */
    1.0 / 14.027090,  /*CH2 */
    1.0 / 14.027090,  /*CH2(S) */
    1.0 / 15.035060,  /*CH3 */
    1.0 / 16.043030,  /*CH4 */
    1.0 / 28.010550,  /*CO */
    1.0 / 44.009950,  /*CO2 */
    1.0 / 29.018520,  /*HCO */
    1.0 / 30.026490,  /*CH2O */
    1.0 / 31.034460,  /*CH3O */
    1.0 / 26.038240,  /*C2H2 */
    1.0 / 27.046210,  /*C2H3 */
    1.0 / 28.054180,  /*C2H4 */
    1.0 / 29.062150,  /*C2H5 */
    1.0 / 30.070120,  /*C2H6 */
    1.0 / 28.013400,  /*N2 */
    1.0 / 39.948000};  /*AR */



static double fwd_A[104], fwd_beta[104], fwd_Ea[104];
static double low_A[104], low_beta[104], low_Ea[104];
static double rev_A[104], rev_beta[104], rev_Ea[104];
static double troe_a[104],troe_Ts[104], troe_Tss[104], troe_Tsss[104];
static double sri_a[104], sri_b[104], sri_c[104], sri_d[104], sri_e[104];
static double activation_units[104], prefactor_units[104], phase_units[104];
static int is_PD[104], troe_len[104], sri_len[104], nTB[104], *TBid[104];
static double *TB[104];

static double fwd_A_DEF[104], fwd_beta_DEF[104], fwd_Ea_DEF[104];
static double low_A_DEF[104], low_beta_DEF[104], low_Ea_DEF[104];
static double rev_A_DEF[104], rev_beta_DEF[104], rev_Ea_DEF[104];
static double troe_a_DEF[104],troe_Ts_DEF[104], troe_Tss_DEF[104], troe_Tsss_DEF[104];
static double sri_a_DEF[104], sri_b_DEF[104], sri_c_DEF[104], sri_d_DEF[104], sri_e_DEF[104];
static double activation_units_DEF[104], prefactor_units_DEF[104], phase_units_DEF[104];
static int is_PD_DEF[104], troe_len_DEF[104], sri_len_DEF[104], nTB_DEF[104], *TBid_DEF[104];
static double *TB_DEF[104];
static int rxn_map[104] = {12,18,19,20,21,22,23,13,24,25,26,27,28,29,30,31,32,33,14,34,35,36,37,38,15,39,40,41,16,42,43,44,0,1,45,2,46,3,47,48,4,5,49,6,50,7,51,8,52,9,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,10,94,95,96,97,98,99,17,100,101,102,11,103};

void GET_REACTION_MAP(int *rmap)
{
    for (int i=0; i<104; ++i) {
        rmap[i] = rxn_map[i];
    }
}


#include <ReactionData.H>
double* GetParamPtr(int                reaction_id,
                    REACTION_PARAMETER param_id,
                    int                species_id,
                    int                get_default)
{
  double* ret = 0;
  if (reaction_id<0 || reaction_id>=104) {
    printf("Bad reaction id = %d",reaction_id);
    abort();
  };
  int mrid = rxn_map[reaction_id];

  if (param_id == THIRD_BODY) {
    if (species_id<0 || species_id>=24) {
      printf("GetParamPtr: Bad species id = %d",species_id);
      abort();
    }
    if (get_default) {
      for (int i=0; i<nTB_DEF[mrid]; ++i) {
        if (species_id == TBid_DEF[mrid][i]) {
          ret = &(TB_DEF[mrid][i]);
        }
      }
    }
    else {
      for (int i=0; i<nTB[mrid]; ++i) {
        if (species_id == TBid[mrid][i]) {
          ret = &(TB[mrid][i]);
        }
      }
    }
    if (ret == 0) {
      printf("GetParamPtr: No TB for reaction id = %d",reaction_id);
      abort();
    }
  }
  else {
    if (     param_id == FWD_A)     {ret = (get_default ? &(fwd_A_DEF[mrid]) : &(fwd_A[mrid]));}
      else if (param_id == FWD_BETA)  {ret = (get_default ? &(fwd_beta_DEF[mrid]) : &(fwd_beta[mrid]));}
      else if (param_id == FWD_EA)    {ret = (get_default ? &(fwd_Ea_DEF[mrid]) : &(fwd_Ea[mrid]));}
      else if (param_id == LOW_A)     {ret = (get_default ? &(low_A_DEF[mrid]) : &(low_A[mrid]));}
      else if (param_id == LOW_BETA)  {ret = (get_default ? &(low_beta_DEF[mrid]) : &(low_beta[mrid]));}
      else if (param_id == LOW_EA)    {ret = (get_default ? &(low_Ea_DEF[mrid]) : &(low_Ea[mrid]));}
      else if (param_id == REV_A)     {ret = (get_default ? &(rev_A_DEF[mrid]) : &(rev_A[mrid]));}
      else if (param_id == REV_BETA)  {ret = (get_default ? &(rev_beta_DEF[mrid]) : &(rev_beta[mrid]));}
      else if (param_id == REV_EA)    {ret = (get_default ? &(rev_Ea_DEF[mrid]) : &(rev_Ea[mrid]));}
      else if (param_id == TROE_A)    {ret = (get_default ? &(troe_a_DEF[mrid]) : &(troe_a[mrid]));}
      else if (param_id == TROE_TS)   {ret = (get_default ? &(troe_Ts_DEF[mrid]) : &(troe_Ts[mrid]));}
      else if (param_id == TROE_TSS)  {ret = (get_default ? &(troe_Tss_DEF[mrid]) : &(troe_Tss[mrid]));}
      else if (param_id == TROE_TSSS) {ret = (get_default ? &(troe_Tsss_DEF[mrid]) : &(troe_Tsss[mrid]));}
      else if (param_id == SRI_A)     {ret = (get_default ? &(sri_a_DEF[mrid]) : &(sri_a[mrid]));}
      else if (param_id == SRI_B)     {ret = (get_default ? &(sri_b_DEF[mrid]) : &(sri_b[mrid]));}
      else if (param_id == SRI_C)     {ret = (get_default ? &(sri_c_DEF[mrid]) : &(sri_c[mrid]));}
      else if (param_id == SRI_D)     {ret = (get_default ? &(sri_d_DEF[mrid]) : &(sri_d[mrid]));}
      else if (param_id == SRI_E)     {ret = (get_default ? &(sri_e_DEF[mrid]) : &(sri_e[mrid]));}
    else {
      printf("GetParamPtr: Unknown parameter id");
      abort();
    }
  }
  return ret;
}

void ResetAllParametersToDefault()
{
    for (int i=0; i<104; i++) {
        if (nTB[i] != 0) {
            nTB[i] = 0;
            free(TB[i]);
            free(TBid[i]);
        }

        fwd_A[i]    = fwd_A_DEF[i];
        fwd_beta[i] = fwd_beta_DEF[i];
        fwd_Ea[i]   = fwd_Ea_DEF[i];

        low_A[i]    = low_A_DEF[i];
        low_beta[i] = low_beta_DEF[i];
        low_Ea[i]   = low_Ea_DEF[i];

        rev_A[i]    = rev_A_DEF[i];
        rev_beta[i] = rev_beta_DEF[i];
        rev_Ea[i]   = rev_Ea_DEF[i];

        troe_a[i]    = troe_a_DEF[i];
        troe_Ts[i]   = troe_Ts_DEF[i];
        troe_Tss[i]  = troe_Tss_DEF[i];
        troe_Tsss[i] = troe_Tsss_DEF[i];

        sri_a[i] = sri_a_DEF[i];
        sri_b[i] = sri_b_DEF[i];
        sri_c[i] = sri_c_DEF[i];
        sri_d[i] = sri_d_DEF[i];
        sri_e[i] = sri_e_DEF[i];

        is_PD[i]    = is_PD_DEF[i];
        troe_len[i] = troe_len_DEF[i];
        sri_len[i]  = sri_len_DEF[i];

        activation_units[i] = activation_units_DEF[i];
        prefactor_units[i]  = prefactor_units_DEF[i];
        phase_units[i]      = phase_units_DEF[i];

        nTB[i]  = nTB_DEF[i];
        if (nTB[i] != 0) {
           TB[i] = (double *) malloc(sizeof(double) * nTB[i]);
           TBid[i] = (int *) malloc(sizeof(int) * nTB[i]);
           for (int j=0; j<nTB[i]; j++) {
             TB[i][j] = TB_DEF[i][j];
             TBid[i][j] = TBid_DEF[i][j];
           }
        }
    }
}

void SetAllDefaults()
{
    for (int i=0; i<104; i++) {
        if (nTB_DEF[i] != 0) {
            nTB_DEF[i] = 0;
            free(TB_DEF[i]);
            free(TBid_DEF[i]);
        }

        fwd_A_DEF[i]    = fwd_A[i];
        fwd_beta_DEF[i] = fwd_beta[i];
        fwd_Ea_DEF[i]   = fwd_Ea[i];

        low_A_DEF[i]    = low_A[i];
        low_beta_DEF[i] = low_beta[i];
        low_Ea_DEF[i]   = low_Ea[i];

        rev_A_DEF[i]    = rev_A[i];
        rev_beta_DEF[i] = rev_beta[i];
        rev_Ea_DEF[i]   = rev_Ea[i];

        troe_a_DEF[i]    = troe_a[i];
        troe_Ts_DEF[i]   = troe_Ts[i];
        troe_Tss_DEF[i]  = troe_Tss[i];
        troe_Tsss_DEF[i] = troe_Tsss[i];

        sri_a_DEF[i] = sri_a[i];
        sri_b_DEF[i] = sri_b[i];
        sri_c_DEF[i] = sri_c[i];
        sri_d_DEF[i] = sri_d[i];
        sri_e_DEF[i] = sri_e[i];

        is_PD_DEF[i]    = is_PD[i];
        troe_len_DEF[i] = troe_len[i];
        sri_len_DEF[i]  = sri_len[i];

        activation_units_DEF[i] = activation_units[i];
        prefactor_units_DEF[i]  = prefactor_units[i];
        phase_units_DEF[i]      = phase_units[i];

        nTB_DEF[i]  = nTB[i];
        if (nTB_DEF[i] != 0) {
           TB_DEF[i] = (double *) malloc(sizeof(double) * nTB_DEF[i]);
           TBid_DEF[i] = (int *) malloc(sizeof(int) * nTB_DEF[i]);
           for (int j=0; j<nTB_DEF[i]; j++) {
             TB_DEF[i][j] = TB[i][j];
             TBid_DEF[i][j] = TBid[i][j];
           }
        }
    }
}

/* Finalizes parameter database */
void CKFINALIZE()
{
  for (int i=0; i<104; ++i) {
    free(TB[i]); TB[i] = 0; 
    free(TBid[i]); TBid[i] = 0;
    nTB[i] = 0;

    free(TB_DEF[i]); TB_DEF[i] = 0; 
    free(TBid_DEF[i]); TBid_DEF[i] = 0;
    nTB_DEF[i] = 0;
  }
}

/* Initializes parameter database */
void CKINIT()
{
    // (0):  O + H + M <=> OH + M
    fwd_A[12]     = 5e+17;
    fwd_beta[12]  = -1;
    fwd_Ea[12]    = 0;
    prefactor_units[12]  = 1.0000000000000002e-12;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = 1e-12;
    is_PD[12] = 0;
    nTB[12] = 7;
    TB[12] = (double *) malloc(7 * sizeof(double));
    TBid[12] = (int *) malloc(7 * sizeof(int));
    TBid[12][0] = 0; TB[12][0] = 2; // H2
    TBid[12][1] = 5; TB[12][1] = 6; // H2O
    TBid[12][2] = 11; TB[12][2] = 2; // CH4
    TBid[12][3] = 12; TB[12][3] = 1.5; // CO
    TBid[12][4] = 13; TB[12][4] = 2; // CO2
    TBid[12][5] = 21; TB[12][5] = 3; // C2H6
    TBid[12][6] = 23; TB[12][6] = 0.69999999999999996; // AR

    // (1):  O + H2 <=> H + OH
    fwd_A[18]     = 50000;
    fwd_beta[18]  = 2.6699999999999999;
    fwd_Ea[18]    = 6290;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = 1e-12;
    is_PD[18] = 0;
    nTB[18] = 0;

    // (2):  O + HO2 <=> OH + O2
    fwd_A[19]     = 20000000000000;
    fwd_beta[19]  = 0;
    fwd_Ea[19]    = 0;
    prefactor_units[19]  = 1.0000000000000002e-06;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = 1e-12;
    is_PD[19] = 0;
    nTB[19] = 0;

    // (3):  O + CH2 <=> H + HCO
    fwd_A[20]     = 80000000000000;
    fwd_beta[20]  = 0;
    fwd_Ea[20]    = 0;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = 1e-12;
    is_PD[20] = 0;
    nTB[20] = 0;

    // (4):  O + CH2(S) <=> H + HCO
    fwd_A[21]     = 15000000000000;
    fwd_beta[21]  = 0;
    fwd_Ea[21]    = 0;
    prefactor_units[21]  = 1.0000000000000002e-06;
    activation_units[21] = 0.50321666580471969;
    phase_units[21]      = 1e-12;
    is_PD[21] = 0;
    nTB[21] = 0;

    // (5):  O + CH3 <=> H + CH2O
    fwd_A[22]     = 84300000000000;
    fwd_beta[22]  = 0;
    fwd_Ea[22]    = 0;
    prefactor_units[22]  = 1.0000000000000002e-06;
    activation_units[22] = 0.50321666580471969;
    phase_units[22]      = 1e-12;
    is_PD[22] = 0;
    nTB[22] = 0;

    // (6):  O + CH4 <=> OH + CH3
    fwd_A[23]     = 1020000000;
    fwd_beta[23]  = 1.5;
    fwd_Ea[23]    = 8600;
    prefactor_units[23]  = 1.0000000000000002e-06;
    activation_units[23] = 0.50321666580471969;
    phase_units[23]      = 1e-12;
    is_PD[23] = 0;
    nTB[23] = 0;

    // (7):  O + CO + M <=> CO2 + M
    fwd_A[13]     = 602000000000000;
    fwd_beta[13]  = 0;
    fwd_Ea[13]    = 3000;
    prefactor_units[13]  = 1.0000000000000002e-12;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = 1e-12;
    is_PD[13] = 0;
    nTB[13] = 8;
    TB[13] = (double *) malloc(8 * sizeof(double));
    TBid[13] = (int *) malloc(8 * sizeof(int));
    TBid[13][0] = 0; TB[13][0] = 2; // H2
    TBid[13][1] = 3; TB[13][1] = 6; // O2
    TBid[13][2] = 5; TB[13][2] = 6; // H2O
    TBid[13][3] = 11; TB[13][3] = 2; // CH4
    TBid[13][4] = 12; TB[13][4] = 1.5; // CO
    TBid[13][5] = 13; TB[13][5] = 3.5; // CO2
    TBid[13][6] = 21; TB[13][6] = 3; // C2H6
    TBid[13][7] = 23; TB[13][7] = 0.5; // AR

    // (8):  O + HCO <=> OH + CO
    fwd_A[24]     = 30000000000000;
    fwd_beta[24]  = 0;
    fwd_Ea[24]    = 0;
    prefactor_units[24]  = 1.0000000000000002e-06;
    activation_units[24] = 0.50321666580471969;
    phase_units[24]      = 1e-12;
    is_PD[24] = 0;
    nTB[24] = 0;

    // (9):  O + HCO <=> H + CO2
    fwd_A[25]     = 30000000000000;
    fwd_beta[25]  = 0;
    fwd_Ea[25]    = 0;
    prefactor_units[25]  = 1.0000000000000002e-06;
    activation_units[25] = 0.50321666580471969;
    phase_units[25]      = 1e-12;
    is_PD[25] = 0;
    nTB[25] = 0;

    // (10):  O + CH2O <=> OH + HCO
    fwd_A[26]     = 39000000000000;
    fwd_beta[26]  = 0;
    fwd_Ea[26]    = 3540;
    prefactor_units[26]  = 1.0000000000000002e-06;
    activation_units[26] = 0.50321666580471969;
    phase_units[26]      = 1e-12;
    is_PD[26] = 0;
    nTB[26] = 0;

    // (11):  O + C2H2 <=> CH2(S) + CO
    fwd_A[27]     = 10200000;
    fwd_beta[27]  = 2;
    fwd_Ea[27]    = 1900;
    prefactor_units[27]  = 1.0000000000000002e-06;
    activation_units[27] = 0.50321666580471969;
    phase_units[27]      = 1e-12;
    is_PD[27] = 0;
    nTB[27] = 0;

    // (12):  O + C2H2 <=> CO + CH2
    fwd_A[28]     = 10200000;
    fwd_beta[28]  = 2;
    fwd_Ea[28]    = 1900;
    prefactor_units[28]  = 1.0000000000000002e-06;
    activation_units[28] = 0.50321666580471969;
    phase_units[28]      = 1e-12;
    is_PD[28] = 0;
    nTB[28] = 0;

    // (13):  O + C2H4 <=> CH3 + HCO
    fwd_A[29]     = 19200000;
    fwd_beta[29]  = 1.8300000000000001;
    fwd_Ea[29]    = 220;
    prefactor_units[29]  = 1.0000000000000002e-06;
    activation_units[29] = 0.50321666580471969;
    phase_units[29]      = 1e-12;
    is_PD[29] = 0;
    nTB[29] = 0;

    // (14):  O + C2H5 <=> CH3 + CH2O
    fwd_A[30]     = 132000000000000;
    fwd_beta[30]  = 0;
    fwd_Ea[30]    = 0;
    prefactor_units[30]  = 1.0000000000000002e-06;
    activation_units[30] = 0.50321666580471969;
    phase_units[30]      = 1e-12;
    is_PD[30] = 0;
    nTB[30] = 0;

    // (15):  O + C2H6 <=> OH + C2H5
    fwd_A[31]     = 89800000;
    fwd_beta[31]  = 1.9199999999999999;
    fwd_Ea[31]    = 5690;
    prefactor_units[31]  = 1.0000000000000002e-06;
    activation_units[31] = 0.50321666580471969;
    phase_units[31]      = 1e-12;
    is_PD[31] = 0;
    nTB[31] = 0;

    // (16):  O2 + CO <=> O + CO2
    fwd_A[32]     = 2500000000000;
    fwd_beta[32]  = 0;
    fwd_Ea[32]    = 47800;
    prefactor_units[32]  = 1.0000000000000002e-06;
    activation_units[32] = 0.50321666580471969;
    phase_units[32]      = 1e-12;
    is_PD[32] = 0;
    nTB[32] = 0;

    // (17):  O2 + CH2O <=> HO2 + HCO
    fwd_A[33]     = 100000000000000;
    fwd_beta[33]  = 0;
    fwd_Ea[33]    = 40000;
    prefactor_units[33]  = 1.0000000000000002e-06;
    activation_units[33] = 0.50321666580471969;
    phase_units[33]      = 1e-12;
    is_PD[33] = 0;
    nTB[33] = 0;

    // (18):  H + O2 + M <=> HO2 + M
    fwd_A[14]     = 2.8e+18;
    fwd_beta[14]  = -0.85999999999999999;
    fwd_Ea[14]    = 0;
    prefactor_units[14]  = 1.0000000000000002e-12;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = 1e-12;
    is_PD[14] = 0;
    nTB[14] = 7;
    TB[14] = (double *) malloc(7 * sizeof(double));
    TBid[14] = (int *) malloc(7 * sizeof(int));
    TBid[14][0] = 3; TB[14][0] = 0; // O2
    TBid[14][1] = 5; TB[14][1] = 0; // H2O
    TBid[14][2] = 12; TB[14][2] = 0.75; // CO
    TBid[14][3] = 13; TB[14][3] = 1.5; // CO2
    TBid[14][4] = 21; TB[14][4] = 1.5; // C2H6
    TBid[14][5] = 22; TB[14][5] = 0; // N2
    TBid[14][6] = 23; TB[14][6] = 0; // AR

    // (19):  H + 2 O2 <=> HO2 + O2
    fwd_A[34]     = 3e+20;
    fwd_beta[34]  = -1.72;
    fwd_Ea[34]    = 0;
    prefactor_units[34]  = 1.0000000000000002e-12;
    activation_units[34] = 0.50321666580471969;
    phase_units[34]      = 1e-18;
    is_PD[34] = 0;
    nTB[34] = 0;

    // (20):  H + O2 + H2O <=> HO2 + H2O
    fwd_A[35]     = 9.38e+18;
    fwd_beta[35]  = -0.76000000000000001;
    fwd_Ea[35]    = 0;
    prefactor_units[35]  = 1.0000000000000002e-12;
    activation_units[35] = 0.50321666580471969;
    phase_units[35]      = 1e-18;
    is_PD[35] = 0;
    nTB[35] = 0;

    // (21):  H + O2 + N2 <=> HO2 + N2
    fwd_A[36]     = 3.75e+20;
    fwd_beta[36]  = -1.72;
    fwd_Ea[36]    = 0;
    prefactor_units[36]  = 1.0000000000000002e-12;
    activation_units[36] = 0.50321666580471969;
    phase_units[36]      = 1e-18;
    is_PD[36] = 0;
    nTB[36] = 0;

    // (22):  H + O2 + AR <=> HO2 + AR
    fwd_A[37]     = 7e+17;
    fwd_beta[37]  = -0.80000000000000004;
    fwd_Ea[37]    = 0;
    prefactor_units[37]  = 1.0000000000000002e-12;
    activation_units[37] = 0.50321666580471969;
    phase_units[37]      = 1e-18;
    is_PD[37] = 0;
    nTB[37] = 0;

    // (23):  H + O2 <=> O + OH
    fwd_A[38]     = 83000000000000;
    fwd_beta[38]  = 0;
    fwd_Ea[38]    = 14413;
    prefactor_units[38]  = 1.0000000000000002e-06;
    activation_units[38] = 0.50321666580471969;
    phase_units[38]      = 1e-12;
    is_PD[38] = 0;
    nTB[38] = 0;

    // (24):  2 H + M <=> H2 + M
    fwd_A[15]     = 1e+18;
    fwd_beta[15]  = -1;
    fwd_Ea[15]    = 0;
    prefactor_units[15]  = 1.0000000000000002e-12;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = 1e-12;
    is_PD[15] = 0;
    nTB[15] = 6;
    TB[15] = (double *) malloc(6 * sizeof(double));
    TBid[15] = (int *) malloc(6 * sizeof(int));
    TBid[15][0] = 0; TB[15][0] = 0; // H2
    TBid[15][1] = 5; TB[15][1] = 0; // H2O
    TBid[15][2] = 11; TB[15][2] = 2; // CH4
    TBid[15][3] = 13; TB[15][3] = 0; // CO2
    TBid[15][4] = 21; TB[15][4] = 3; // C2H6
    TBid[15][5] = 23; TB[15][5] = 0.63; // AR

    // (25):  2 H + H2 <=> 2 H2
    fwd_A[39]     = 90000000000000000;
    fwd_beta[39]  = -0.59999999999999998;
    fwd_Ea[39]    = 0;
    prefactor_units[39]  = 1.0000000000000002e-12;
    activation_units[39] = 0.50321666580471969;
    phase_units[39]      = 1e-18;
    is_PD[39] = 0;
    nTB[39] = 0;

    // (26):  2 H + H2O <=> H2 + H2O
    fwd_A[40]     = 6e+19;
    fwd_beta[40]  = -1.25;
    fwd_Ea[40]    = 0;
    prefactor_units[40]  = 1.0000000000000002e-12;
    activation_units[40] = 0.50321666580471969;
    phase_units[40]      = 1e-18;
    is_PD[40] = 0;
    nTB[40] = 0;

    // (27):  2 H + CO2 <=> H2 + CO2
    fwd_A[41]     = 5.5e+20;
    fwd_beta[41]  = -2;
    fwd_Ea[41]    = 0;
    prefactor_units[41]  = 1.0000000000000002e-12;
    activation_units[41] = 0.50321666580471969;
    phase_units[41]      = 1e-18;
    is_PD[41] = 0;
    nTB[41] = 0;

    // (28):  H + OH + M <=> H2O + M
    fwd_A[16]     = 2.2e+22;
    fwd_beta[16]  = -2;
    fwd_Ea[16]    = 0;
    prefactor_units[16]  = 1.0000000000000002e-12;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = 1e-12;
    is_PD[16] = 0;
    nTB[16] = 5;
    TB[16] = (double *) malloc(5 * sizeof(double));
    TBid[16] = (int *) malloc(5 * sizeof(int));
    TBid[16][0] = 0; TB[16][0] = 0.72999999999999998; // H2
    TBid[16][1] = 5; TB[16][1] = 3.6499999999999999; // H2O
    TBid[16][2] = 11; TB[16][2] = 2; // CH4
    TBid[16][3] = 21; TB[16][3] = 3; // C2H6
    TBid[16][4] = 23; TB[16][4] = 0.38; // AR

    // (29):  H + HO2 <=> O2 + H2
    fwd_A[42]     = 28000000000000;
    fwd_beta[42]  = 0;
    fwd_Ea[42]    = 1068;
    prefactor_units[42]  = 1.0000000000000002e-06;
    activation_units[42] = 0.50321666580471969;
    phase_units[42]      = 1e-12;
    is_PD[42] = 0;
    nTB[42] = 0;

    // (30):  H + HO2 <=> 2 OH
    fwd_A[43]     = 134000000000000;
    fwd_beta[43]  = 0;
    fwd_Ea[43]    = 635;
    prefactor_units[43]  = 1.0000000000000002e-06;
    activation_units[43] = 0.50321666580471969;
    phase_units[43]      = 1e-12;
    is_PD[43] = 0;
    nTB[43] = 0;

    // (31):  H + H2O2 <=> HO2 + H2
    fwd_A[44]     = 12100000;
    fwd_beta[44]  = 2;
    fwd_Ea[44]    = 5200;
    prefactor_units[44]  = 1.0000000000000002e-06;
    activation_units[44] = 0.50321666580471969;
    phase_units[44]      = 1e-12;
    is_PD[44] = 0;
    nTB[44] = 0;

    // (32):  H + CH2 (+M) <=> CH3 (+M)
    fwd_A[0]     = 25000000000000000;
    fwd_beta[0]  = -0.80000000000000004;
    fwd_Ea[0]    = 0;
    low_A[0]     = 3.2000000000000002e+27;
    low_beta[0]  = -3.1400000000000001;
    low_Ea[0]    = 1230;
    troe_a[0]    = 0.68000000000000005;
    troe_Tsss[0] = 78;
    troe_Ts[0]   = 1995;
    troe_Tss[0]  = 5590;
    troe_len[0]  = 4;
    prefactor_units[0]  = 1.0000000000000002e-06;
    activation_units[0] = 0.50321666580471969;
    phase_units[0]      = 1e-12;
    is_PD[0] = 1;
    nTB[0] = 7;
    TB[0] = (double *) malloc(7 * sizeof(double));
    TBid[0] = (int *) malloc(7 * sizeof(int));
    TBid[0][0] = 0; TB[0][0] = 2; // H2
    TBid[0][1] = 5; TB[0][1] = 6; // H2O
    TBid[0][2] = 11; TB[0][2] = 2; // CH4
    TBid[0][3] = 12; TB[0][3] = 1.5; // CO
    TBid[0][4] = 13; TB[0][4] = 2; // CO2
    TBid[0][5] = 21; TB[0][5] = 3; // C2H6
    TBid[0][6] = 23; TB[0][6] = 0.69999999999999996; // AR

    // (33):  H + CH3 (+M) <=> CH4 (+M)
    fwd_A[1]     = 12700000000000000;
    fwd_beta[1]  = -0.63;
    fwd_Ea[1]    = 383;
    low_A[1]     = 2.4769999999999999e+33;
    low_beta[1]  = -4.7599999999999998;
    low_Ea[1]    = 2440;
    troe_a[1]    = 0.78300000000000003;
    troe_Tsss[1] = 74;
    troe_Ts[1]   = 2941;
    troe_Tss[1]  = 6964;
    troe_len[1]  = 4;
    prefactor_units[1]  = 1.0000000000000002e-06;
    activation_units[1] = 0.50321666580471969;
    phase_units[1]      = 1e-12;
    is_PD[1] = 1;
    nTB[1] = 7;
    TB[1] = (double *) malloc(7 * sizeof(double));
    TBid[1] = (int *) malloc(7 * sizeof(int));
    TBid[1][0] = 0; TB[1][0] = 2; // H2
    TBid[1][1] = 5; TB[1][1] = 6; // H2O
    TBid[1][2] = 11; TB[1][2] = 2; // CH4
    TBid[1][3] = 12; TB[1][3] = 1.5; // CO
    TBid[1][4] = 13; TB[1][4] = 2; // CO2
    TBid[1][5] = 21; TB[1][5] = 3; // C2H6
    TBid[1][6] = 23; TB[1][6] = 0.69999999999999996; // AR

    // (34):  H + CH4 <=> CH3 + H2
    fwd_A[45]     = 660000000;
    fwd_beta[45]  = 1.6200000000000001;
    fwd_Ea[45]    = 10840;
    prefactor_units[45]  = 1.0000000000000002e-06;
    activation_units[45] = 0.50321666580471969;
    phase_units[45]      = 1e-12;
    is_PD[45] = 0;
    nTB[45] = 0;

    // (35):  H + HCO (+M) <=> CH2O (+M)
    fwd_A[2]     = 1090000000000;
    fwd_beta[2]  = 0.47999999999999998;
    fwd_Ea[2]    = -260;
    low_A[2]     = 1.35e+24;
    low_beta[2]  = -2.5699999999999998;
    low_Ea[2]    = 1425;
    troe_a[2]    = 0.78239999999999998;
    troe_Tsss[2] = 271;
    troe_Ts[2]   = 2755;
    troe_Tss[2]  = 6570;
    troe_len[2]  = 4;
    prefactor_units[2]  = 1.0000000000000002e-06;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = 1e-12;
    is_PD[2] = 1;
    nTB[2] = 7;
    TB[2] = (double *) malloc(7 * sizeof(double));
    TBid[2] = (int *) malloc(7 * sizeof(int));
    TBid[2][0] = 0; TB[2][0] = 2; // H2
    TBid[2][1] = 5; TB[2][1] = 6; // H2O
    TBid[2][2] = 11; TB[2][2] = 2; // CH4
    TBid[2][3] = 12; TB[2][3] = 1.5; // CO
    TBid[2][4] = 13; TB[2][4] = 2; // CO2
    TBid[2][5] = 21; TB[2][5] = 3; // C2H6
    TBid[2][6] = 23; TB[2][6] = 0.69999999999999996; // AR

    // (36):  H + HCO <=> H2 + CO
    fwd_A[46]     = 73400000000000;
    fwd_beta[46]  = 0;
    fwd_Ea[46]    = 0;
    prefactor_units[46]  = 1.0000000000000002e-06;
    activation_units[46] = 0.50321666580471969;
    phase_units[46]      = 1e-12;
    is_PD[46] = 0;
    nTB[46] = 0;

    // (37):  H + CH2O (+M) <=> CH3O (+M)
    fwd_A[3]     = 540000000000;
    fwd_beta[3]  = 0.45400000000000001;
    fwd_Ea[3]    = 2600;
    low_A[3]     = 2.2e+30;
    low_beta[3]  = -4.7999999999999998;
    low_Ea[3]    = 5560;
    troe_a[3]    = 0.75800000000000001;
    troe_Tsss[3] = 94;
    troe_Ts[3]   = 1555;
    troe_Tss[3]  = 4200;
    troe_len[3]  = 4;
    prefactor_units[3]  = 1.0000000000000002e-06;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = 1e-12;
    is_PD[3] = 1;
    nTB[3] = 6;
    TB[3] = (double *) malloc(6 * sizeof(double));
    TBid[3] = (int *) malloc(6 * sizeof(int));
    TBid[3][0] = 0; TB[3][0] = 2; // H2
    TBid[3][1] = 5; TB[3][1] = 6; // H2O
    TBid[3][2] = 11; TB[3][2] = 2; // CH4
    TBid[3][3] = 12; TB[3][3] = 1.5; // CO
    TBid[3][4] = 13; TB[3][4] = 2; // CO2
    TBid[3][5] = 21; TB[3][5] = 3; // C2H6

    // (38):  H + CH2O <=> HCO + H2
    fwd_A[47]     = 23000000000;
    fwd_beta[47]  = 1.05;
    fwd_Ea[47]    = 3275;
    prefactor_units[47]  = 1.0000000000000002e-06;
    activation_units[47] = 0.50321666580471969;
    phase_units[47]      = 1e-12;
    is_PD[47] = 0;
    nTB[47] = 0;

    // (39):  H + CH3O <=> OH + CH3
    fwd_A[48]     = 32000000000000;
    fwd_beta[48]  = 0;
    fwd_Ea[48]    = 0;
    prefactor_units[48]  = 1.0000000000000002e-06;
    activation_units[48] = 0.50321666580471969;
    phase_units[48]      = 1e-12;
    is_PD[48] = 0;
    nTB[48] = 0;

    // (40):  H + C2H2 (+M) <=> C2H3 (+M)
    fwd_A[4]     = 5600000000000;
    fwd_beta[4]  = 0;
    fwd_Ea[4]    = 2400;
    low_A[4]     = 3.8e+40;
    low_beta[4]  = -7.2699999999999996;
    low_Ea[4]    = 7220;
    troe_a[4]    = 0.75070000000000003;
    troe_Tsss[4] = 98.5;
    troe_Ts[4]   = 1302;
    troe_Tss[4]  = 4167;
    troe_len[4]  = 4;
    prefactor_units[4]  = 1.0000000000000002e-06;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = 1e-12;
    is_PD[4] = 1;
    nTB[4] = 7;
    TB[4] = (double *) malloc(7 * sizeof(double));
    TBid[4] = (int *) malloc(7 * sizeof(int));
    TBid[4][0] = 0; TB[4][0] = 2; // H2
    TBid[4][1] = 5; TB[4][1] = 6; // H2O
    TBid[4][2] = 11; TB[4][2] = 2; // CH4
    TBid[4][3] = 12; TB[4][3] = 1.5; // CO
    TBid[4][4] = 13; TB[4][4] = 2; // CO2
    TBid[4][5] = 21; TB[4][5] = 3; // C2H6
    TBid[4][6] = 23; TB[4][6] = 0.69999999999999996; // AR

    // (41):  H + C2H3 (+M) <=> C2H4 (+M)
    fwd_A[5]     = 6080000000000;
    fwd_beta[5]  = 0.27000000000000002;
    fwd_Ea[5]    = 280;
    low_A[5]     = 1.3999999999999999e+30;
    low_beta[5]  = -3.8599999999999999;
    low_Ea[5]    = 3320;
    troe_a[5]    = 0.78200000000000003;
    troe_Tsss[5] = 207.5;
    troe_Ts[5]   = 2663;
    troe_Tss[5]  = 6095;
    troe_len[5]  = 4;
    prefactor_units[5]  = 1.0000000000000002e-06;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = 1e-12;
    is_PD[5] = 1;
    nTB[5] = 7;
    TB[5] = (double *) malloc(7 * sizeof(double));
    TBid[5] = (int *) malloc(7 * sizeof(int));
    TBid[5][0] = 0; TB[5][0] = 2; // H2
    TBid[5][1] = 5; TB[5][1] = 6; // H2O
    TBid[5][2] = 11; TB[5][2] = 2; // CH4
    TBid[5][3] = 12; TB[5][3] = 1.5; // CO
    TBid[5][4] = 13; TB[5][4] = 2; // CO2
    TBid[5][5] = 21; TB[5][5] = 3; // C2H6
    TBid[5][6] = 23; TB[5][6] = 0.69999999999999996; // AR

    // (42):  H + C2H3 <=> H2 + C2H2
    fwd_A[49]     = 30000000000000;
    fwd_beta[49]  = 0;
    fwd_Ea[49]    = 0;
    prefactor_units[49]  = 1.0000000000000002e-06;
    activation_units[49] = 0.50321666580471969;
    phase_units[49]      = 1e-12;
    is_PD[49] = 0;
    nTB[49] = 0;

    // (43):  H + C2H4 (+M) <=> C2H5 (+M)
    fwd_A[6]     = 1080000000000;
    fwd_beta[6]  = 0.45400000000000001;
    fwd_Ea[6]    = 1820;
    low_A[6]     = 1.1999999999999999e+42;
    low_beta[6]  = -7.6200000000000001;
    low_Ea[6]    = 6970;
    troe_a[6]    = 0.97529999999999994;
    troe_Tsss[6] = 210;
    troe_Ts[6]   = 984;
    troe_Tss[6]  = 4374;
    troe_len[6]  = 4;
    prefactor_units[6]  = 1.0000000000000002e-06;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = 1e-12;
    is_PD[6] = 1;
    nTB[6] = 7;
    TB[6] = (double *) malloc(7 * sizeof(double));
    TBid[6] = (int *) malloc(7 * sizeof(int));
    TBid[6][0] = 0; TB[6][0] = 2; // H2
    TBid[6][1] = 5; TB[6][1] = 6; // H2O
    TBid[6][2] = 11; TB[6][2] = 2; // CH4
    TBid[6][3] = 12; TB[6][3] = 1.5; // CO
    TBid[6][4] = 13; TB[6][4] = 2; // CO2
    TBid[6][5] = 21; TB[6][5] = 3; // C2H6
    TBid[6][6] = 23; TB[6][6] = 0.69999999999999996; // AR

    // (44):  H + C2H4 <=> C2H3 + H2
    fwd_A[50]     = 1325000;
    fwd_beta[50]  = 2.5299999999999998;
    fwd_Ea[50]    = 12240;
    prefactor_units[50]  = 1.0000000000000002e-06;
    activation_units[50] = 0.50321666580471969;
    phase_units[50]      = 1e-12;
    is_PD[50] = 0;
    nTB[50] = 0;

    // (45):  H + C2H5 (+M) <=> C2H6 (+M)
    fwd_A[7]     = 5.21e+17;
    fwd_beta[7]  = -0.98999999999999999;
    fwd_Ea[7]    = 1580;
    low_A[7]     = 1.9900000000000001e+41;
    low_beta[7]  = -7.0800000000000001;
    low_Ea[7]    = 6685;
    troe_a[7]    = 0.84219999999999995;
    troe_Tsss[7] = 125;
    troe_Ts[7]   = 2219;
    troe_Tss[7]  = 6882;
    troe_len[7]  = 4;
    prefactor_units[7]  = 1.0000000000000002e-06;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = 1e-12;
    is_PD[7] = 1;
    nTB[7] = 7;
    TB[7] = (double *) malloc(7 * sizeof(double));
    TBid[7] = (int *) malloc(7 * sizeof(int));
    TBid[7][0] = 0; TB[7][0] = 2; // H2
    TBid[7][1] = 5; TB[7][1] = 6; // H2O
    TBid[7][2] = 11; TB[7][2] = 2; // CH4
    TBid[7][3] = 12; TB[7][3] = 1.5; // CO
    TBid[7][4] = 13; TB[7][4] = 2; // CO2
    TBid[7][5] = 21; TB[7][5] = 3; // C2H6
    TBid[7][6] = 23; TB[7][6] = 0.69999999999999996; // AR

    // (46):  H + C2H6 <=> C2H5 + H2
    fwd_A[51]     = 115000000;
    fwd_beta[51]  = 1.8999999999999999;
    fwd_Ea[51]    = 7530;
    prefactor_units[51]  = 1.0000000000000002e-06;
    activation_units[51] = 0.50321666580471969;
    phase_units[51]      = 1e-12;
    is_PD[51] = 0;
    nTB[51] = 0;

    // (47):  H2 + CO (+M) <=> CH2O (+M)
    fwd_A[8]     = 43000000;
    fwd_beta[8]  = 1.5;
    fwd_Ea[8]    = 79600;
    low_A[8]     = 5.0699999999999998e+27;
    low_beta[8]  = -3.4199999999999999;
    low_Ea[8]    = 84350;
    troe_a[8]    = 0.93200000000000005;
    troe_Tsss[8] = 197;
    troe_Ts[8]   = 1540;
    troe_Tss[8]  = 10300;
    troe_len[8]  = 4;
    prefactor_units[8]  = 1.0000000000000002e-06;
    activation_units[8] = 0.50321666580471969;
    phase_units[8]      = 1e-12;
    is_PD[8] = 1;
    nTB[8] = 7;
    TB[8] = (double *) malloc(7 * sizeof(double));
    TBid[8] = (int *) malloc(7 * sizeof(int));
    TBid[8][0] = 0; TB[8][0] = 2; // H2
    TBid[8][1] = 5; TB[8][1] = 6; // H2O
    TBid[8][2] = 11; TB[8][2] = 2; // CH4
    TBid[8][3] = 12; TB[8][3] = 1.5; // CO
    TBid[8][4] = 13; TB[8][4] = 2; // CO2
    TBid[8][5] = 21; TB[8][5] = 3; // C2H6
    TBid[8][6] = 23; TB[8][6] = 0.69999999999999996; // AR

    // (48):  OH + H2 <=> H + H2O
    fwd_A[52]     = 216000000;
    fwd_beta[52]  = 1.51;
    fwd_Ea[52]    = 3430;
    prefactor_units[52]  = 1.0000000000000002e-06;
    activation_units[52] = 0.50321666580471969;
    phase_units[52]      = 1e-12;
    is_PD[52] = 0;
    nTB[52] = 0;

    // (49):  2 OH (+M) <=> H2O2 (+M)
    fwd_A[9]     = 74000000000000;
    fwd_beta[9]  = -0.37;
    fwd_Ea[9]    = 0;
    low_A[9]     = 2.3e+18;
    low_beta[9]  = -0.90000000000000002;
    low_Ea[9]    = -1700;
    troe_a[9]    = 0.73460000000000003;
    troe_Tsss[9] = 94;
    troe_Ts[9]   = 1756;
    troe_Tss[9]  = 5182;
    troe_len[9]  = 4;
    prefactor_units[9]  = 1.0000000000000002e-06;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = 1e-12;
    is_PD[9] = 1;
    nTB[9] = 7;
    TB[9] = (double *) malloc(7 * sizeof(double));
    TBid[9] = (int *) malloc(7 * sizeof(int));
    TBid[9][0] = 0; TB[9][0] = 2; // H2
    TBid[9][1] = 5; TB[9][1] = 6; // H2O
    TBid[9][2] = 11; TB[9][2] = 2; // CH4
    TBid[9][3] = 12; TB[9][3] = 1.5; // CO
    TBid[9][4] = 13; TB[9][4] = 2; // CO2
    TBid[9][5] = 21; TB[9][5] = 3; // C2H6
    TBid[9][6] = 23; TB[9][6] = 0.69999999999999996; // AR

    // (50):  2 OH <=> O + H2O
    fwd_A[53]     = 35700;
    fwd_beta[53]  = 2.3999999999999999;
    fwd_Ea[53]    = -2110;
    prefactor_units[53]  = 1.0000000000000002e-06;
    activation_units[53] = 0.50321666580471969;
    phase_units[53]      = 1e-12;
    is_PD[53] = 0;
    nTB[53] = 0;

    // (51):  OH + HO2 <=> O2 + H2O
    fwd_A[54]     = 29000000000000;
    fwd_beta[54]  = 0;
    fwd_Ea[54]    = -500;
    prefactor_units[54]  = 1.0000000000000002e-06;
    activation_units[54] = 0.50321666580471969;
    phase_units[54]      = 1e-12;
    is_PD[54] = 0;
    nTB[54] = 0;

    // (52):  OH + H2O2 <=> HO2 + H2O
    fwd_A[55]     = 580000000000000;
    fwd_beta[55]  = 0;
    fwd_Ea[55]    = 9560;
    prefactor_units[55]  = 1.0000000000000002e-06;
    activation_units[55] = 0.50321666580471969;
    phase_units[55]      = 1e-12;
    is_PD[55] = 0;
    nTB[55] = 0;

    // (53):  OH + CH2 <=> H + CH2O
    fwd_A[56]     = 20000000000000;
    fwd_beta[56]  = 0;
    fwd_Ea[56]    = 0;
    prefactor_units[56]  = 1.0000000000000002e-06;
    activation_units[56] = 0.50321666580471969;
    phase_units[56]      = 1e-12;
    is_PD[56] = 0;
    nTB[56] = 0;

    // (54):  OH + CH2(S) <=> H + CH2O
    fwd_A[57]     = 30000000000000;
    fwd_beta[57]  = 0;
    fwd_Ea[57]    = 0;
    prefactor_units[57]  = 1.0000000000000002e-06;
    activation_units[57] = 0.50321666580471969;
    phase_units[57]      = 1e-12;
    is_PD[57] = 0;
    nTB[57] = 0;

    // (55):  OH + CH3 <=> CH2 + H2O
    fwd_A[58]     = 56000000;
    fwd_beta[58]  = 1.6000000000000001;
    fwd_Ea[58]    = 5420;
    prefactor_units[58]  = 1.0000000000000002e-06;
    activation_units[58] = 0.50321666580471969;
    phase_units[58]      = 1e-12;
    is_PD[58] = 0;
    nTB[58] = 0;

    // (56):  OH + CH3 <=> CH2(S) + H2O
    fwd_A[59]     = 25010000000000;
    fwd_beta[59]  = 0;
    fwd_Ea[59]    = 0;
    prefactor_units[59]  = 1.0000000000000002e-06;
    activation_units[59] = 0.50321666580471969;
    phase_units[59]      = 1e-12;
    is_PD[59] = 0;
    nTB[59] = 0;

    // (57):  OH + CH4 <=> CH3 + H2O
    fwd_A[60]     = 100000000;
    fwd_beta[60]  = 1.6000000000000001;
    fwd_Ea[60]    = 3120;
    prefactor_units[60]  = 1.0000000000000002e-06;
    activation_units[60] = 0.50321666580471969;
    phase_units[60]      = 1e-12;
    is_PD[60] = 0;
    nTB[60] = 0;

    // (58):  OH + CO <=> H + CO2
    fwd_A[61]     = 47600000;
    fwd_beta[61]  = 1.228;
    fwd_Ea[61]    = 70;
    prefactor_units[61]  = 1.0000000000000002e-06;
    activation_units[61] = 0.50321666580471969;
    phase_units[61]      = 1e-12;
    is_PD[61] = 0;
    nTB[61] = 0;

    // (59):  OH + HCO <=> H2O + CO
    fwd_A[62]     = 50000000000000;
    fwd_beta[62]  = 0;
    fwd_Ea[62]    = 0;
    prefactor_units[62]  = 1.0000000000000002e-06;
    activation_units[62] = 0.50321666580471969;
    phase_units[62]      = 1e-12;
    is_PD[62] = 0;
    nTB[62] = 0;

    // (60):  OH + CH2O <=> HCO + H2O
    fwd_A[63]     = 3430000000;
    fwd_beta[63]  = 1.1799999999999999;
    fwd_Ea[63]    = -447;
    prefactor_units[63]  = 1.0000000000000002e-06;
    activation_units[63] = 0.50321666580471969;
    phase_units[63]      = 1e-12;
    is_PD[63] = 0;
    nTB[63] = 0;

    // (61):  OH + C2H2 <=> CH3 + CO
    fwd_A[64]     = 0.00048299999999999998;
    fwd_beta[64]  = 4;
    fwd_Ea[64]    = -2000;
    prefactor_units[64]  = 1.0000000000000002e-06;
    activation_units[64] = 0.50321666580471969;
    phase_units[64]      = 1e-12;
    is_PD[64] = 0;
    nTB[64] = 0;

    // (62):  OH + C2H3 <=> H2O + C2H2
    fwd_A[65]     = 5000000000000;
    fwd_beta[65]  = 0;
    fwd_Ea[65]    = 0;
    prefactor_units[65]  = 1.0000000000000002e-06;
    activation_units[65] = 0.50321666580471969;
    phase_units[65]      = 1e-12;
    is_PD[65] = 0;
    nTB[65] = 0;

    // (63):  OH + C2H4 <=> C2H3 + H2O
    fwd_A[66]     = 3600000;
    fwd_beta[66]  = 2;
    fwd_Ea[66]    = 2500;
    prefactor_units[66]  = 1.0000000000000002e-06;
    activation_units[66] = 0.50321666580471969;
    phase_units[66]      = 1e-12;
    is_PD[66] = 0;
    nTB[66] = 0;

    // (64):  OH + C2H6 <=> C2H5 + H2O
    fwd_A[67]     = 3540000;
    fwd_beta[67]  = 2.1200000000000001;
    fwd_Ea[67]    = 870;
    prefactor_units[67]  = 1.0000000000000002e-06;
    activation_units[67] = 0.50321666580471969;
    phase_units[67]      = 1e-12;
    is_PD[67] = 0;
    nTB[67] = 0;

    // (65):  2 HO2 <=> O2 + H2O2
    fwd_A[68]     = 130000000000;
    fwd_beta[68]  = 0;
    fwd_Ea[68]    = -1630;
    prefactor_units[68]  = 1.0000000000000002e-06;
    activation_units[68] = 0.50321666580471969;
    phase_units[68]      = 1e-12;
    is_PD[68] = 0;
    nTB[68] = 0;

    // (66):  2 HO2 <=> O2 + H2O2
    fwd_A[69]     = 420000000000000;
    fwd_beta[69]  = 0;
    fwd_Ea[69]    = 12000;
    prefactor_units[69]  = 1.0000000000000002e-06;
    activation_units[69] = 0.50321666580471969;
    phase_units[69]      = 1e-12;
    is_PD[69] = 0;
    nTB[69] = 0;

    // (67):  HO2 + CH2 <=> OH + CH2O
    fwd_A[70]     = 20000000000000;
    fwd_beta[70]  = 0;
    fwd_Ea[70]    = 0;
    prefactor_units[70]  = 1.0000000000000002e-06;
    activation_units[70] = 0.50321666580471969;
    phase_units[70]      = 1e-12;
    is_PD[70] = 0;
    nTB[70] = 0;

    // (68):  HO2 + CH3 <=> O2 + CH4
    fwd_A[71]     = 1000000000000;
    fwd_beta[71]  = 0;
    fwd_Ea[71]    = 0;
    prefactor_units[71]  = 1.0000000000000002e-06;
    activation_units[71] = 0.50321666580471969;
    phase_units[71]      = 1e-12;
    is_PD[71] = 0;
    nTB[71] = 0;

    // (69):  HO2 + CH3 <=> OH + CH3O
    fwd_A[72]     = 20000000000000;
    fwd_beta[72]  = 0;
    fwd_Ea[72]    = 0;
    prefactor_units[72]  = 1.0000000000000002e-06;
    activation_units[72] = 0.50321666580471969;
    phase_units[72]      = 1e-12;
    is_PD[72] = 0;
    nTB[72] = 0;

    // (70):  HO2 + CO <=> OH + CO2
    fwd_A[73]     = 150000000000000;
    fwd_beta[73]  = 0;
    fwd_Ea[73]    = 23600;
    prefactor_units[73]  = 1.0000000000000002e-06;
    activation_units[73] = 0.50321666580471969;
    phase_units[73]      = 1e-12;
    is_PD[73] = 0;
    nTB[73] = 0;

    // (71):  HO2 + CH2O <=> HCO + H2O2
    fwd_A[74]     = 1000000000000;
    fwd_beta[74]  = 0;
    fwd_Ea[74]    = 8000;
    prefactor_units[74]  = 1.0000000000000002e-06;
    activation_units[74] = 0.50321666580471969;
    phase_units[74]      = 1e-12;
    is_PD[74] = 0;
    nTB[74] = 0;

    // (72):  CH2 + O2 <=> OH + HCO
    fwd_A[75]     = 13200000000000;
    fwd_beta[75]  = 0;
    fwd_Ea[75]    = 1500;
    prefactor_units[75]  = 1.0000000000000002e-06;
    activation_units[75] = 0.50321666580471969;
    phase_units[75]      = 1e-12;
    is_PD[75] = 0;
    nTB[75] = 0;

    // (73):  CH2 + H2 <=> H + CH3
    fwd_A[76]     = 500000;
    fwd_beta[76]  = 2;
    fwd_Ea[76]    = 7230;
    prefactor_units[76]  = 1.0000000000000002e-06;
    activation_units[76] = 0.50321666580471969;
    phase_units[76]      = 1e-12;
    is_PD[76] = 0;
    nTB[76] = 0;

    // (74):  2 CH2 <=> H2 + C2H2
    fwd_A[77]     = 32000000000000;
    fwd_beta[77]  = 0;
    fwd_Ea[77]    = 0;
    prefactor_units[77]  = 1.0000000000000002e-06;
    activation_units[77] = 0.50321666580471969;
    phase_units[77]      = 1e-12;
    is_PD[77] = 0;
    nTB[77] = 0;

    // (75):  CH2 + CH3 <=> H + C2H4
    fwd_A[78]     = 40000000000000;
    fwd_beta[78]  = 0;
    fwd_Ea[78]    = 0;
    prefactor_units[78]  = 1.0000000000000002e-06;
    activation_units[78] = 0.50321666580471969;
    phase_units[78]      = 1e-12;
    is_PD[78] = 0;
    nTB[78] = 0;

    // (76):  CH2 + CH4 <=> 2 CH3
    fwd_A[79]     = 2460000;
    fwd_beta[79]  = 2;
    fwd_Ea[79]    = 8270;
    prefactor_units[79]  = 1.0000000000000002e-06;
    activation_units[79] = 0.50321666580471969;
    phase_units[79]      = 1e-12;
    is_PD[79] = 0;
    nTB[79] = 0;

    // (77):  CH2(S) + N2 <=> CH2 + N2
    fwd_A[80]     = 15000000000000;
    fwd_beta[80]  = 0;
    fwd_Ea[80]    = 600;
    prefactor_units[80]  = 1.0000000000000002e-06;
    activation_units[80] = 0.50321666580471969;
    phase_units[80]      = 1e-12;
    is_PD[80] = 0;
    nTB[80] = 0;

    // (78):  CH2(S) + AR <=> CH2 + AR
    fwd_A[81]     = 9000000000000;
    fwd_beta[81]  = 0;
    fwd_Ea[81]    = 600;
    prefactor_units[81]  = 1.0000000000000002e-06;
    activation_units[81] = 0.50321666580471969;
    phase_units[81]      = 1e-12;
    is_PD[81] = 0;
    nTB[81] = 0;

    // (79):  CH2(S) + O2 <=> H + OH + CO
    fwd_A[82]     = 28000000000000;
    fwd_beta[82]  = 0;
    fwd_Ea[82]    = 0;
    prefactor_units[82]  = 1.0000000000000002e-06;
    activation_units[82] = 0.50321666580471969;
    phase_units[82]      = 1e-12;
    is_PD[82] = 0;
    nTB[82] = 0;

    // (80):  CH2(S) + O2 <=> CO + H2O
    fwd_A[83]     = 12000000000000;
    fwd_beta[83]  = 0;
    fwd_Ea[83]    = 0;
    prefactor_units[83]  = 1.0000000000000002e-06;
    activation_units[83] = 0.50321666580471969;
    phase_units[83]      = 1e-12;
    is_PD[83] = 0;
    nTB[83] = 0;

    // (81):  CH2(S) + H2 <=> CH3 + H
    fwd_A[84]     = 70000000000000;
    fwd_beta[84]  = 0;
    fwd_Ea[84]    = 0;
    prefactor_units[84]  = 1.0000000000000002e-06;
    activation_units[84] = 0.50321666580471969;
    phase_units[84]      = 1e-12;
    is_PD[84] = 0;
    nTB[84] = 0;

    // (82):  CH2(S) + H2O <=> CH2 + H2O
    fwd_A[85]     = 30000000000000;
    fwd_beta[85]  = 0;
    fwd_Ea[85]    = 0;
    prefactor_units[85]  = 1.0000000000000002e-06;
    activation_units[85] = 0.50321666580471969;
    phase_units[85]      = 1e-12;
    is_PD[85] = 0;
    nTB[85] = 0;

    // (83):  CH2(S) + CH3 <=> H + C2H4
    fwd_A[86]     = 12000000000000;
    fwd_beta[86]  = 0;
    fwd_Ea[86]    = -570;
    prefactor_units[86]  = 1.0000000000000002e-06;
    activation_units[86] = 0.50321666580471969;
    phase_units[86]      = 1e-12;
    is_PD[86] = 0;
    nTB[86] = 0;

    // (84):  CH2(S) + CH4 <=> 2 CH3
    fwd_A[87]     = 16000000000000;
    fwd_beta[87]  = 0;
    fwd_Ea[87]    = -570;
    prefactor_units[87]  = 1.0000000000000002e-06;
    activation_units[87] = 0.50321666580471969;
    phase_units[87]      = 1e-12;
    is_PD[87] = 0;
    nTB[87] = 0;

    // (85):  CH2(S) + CO <=> CH2 + CO
    fwd_A[88]     = 9000000000000;
    fwd_beta[88]  = 0;
    fwd_Ea[88]    = 0;
    prefactor_units[88]  = 1.0000000000000002e-06;
    activation_units[88] = 0.50321666580471969;
    phase_units[88]      = 1e-12;
    is_PD[88] = 0;
    nTB[88] = 0;

    // (86):  CH2(S) + CO2 <=> CH2 + CO2
    fwd_A[89]     = 7000000000000;
    fwd_beta[89]  = 0;
    fwd_Ea[89]    = 0;
    prefactor_units[89]  = 1.0000000000000002e-06;
    activation_units[89] = 0.50321666580471969;
    phase_units[89]      = 1e-12;
    is_PD[89] = 0;
    nTB[89] = 0;

    // (87):  CH2(S) + CO2 <=> CO + CH2O
    fwd_A[90]     = 14000000000000;
    fwd_beta[90]  = 0;
    fwd_Ea[90]    = 0;
    prefactor_units[90]  = 1.0000000000000002e-06;
    activation_units[90] = 0.50321666580471969;
    phase_units[90]      = 1e-12;
    is_PD[90] = 0;
    nTB[90] = 0;

    // (88):  CH3 + O2 <=> O + CH3O
    fwd_A[91]     = 26750000000000;
    fwd_beta[91]  = 0;
    fwd_Ea[91]    = 28800;
    prefactor_units[91]  = 1.0000000000000002e-06;
    activation_units[91] = 0.50321666580471969;
    phase_units[91]      = 1e-12;
    is_PD[91] = 0;
    nTB[91] = 0;

    // (89):  CH3 + O2 <=> OH + CH2O
    fwd_A[92]     = 36000000000;
    fwd_beta[92]  = 0;
    fwd_Ea[92]    = 8940;
    prefactor_units[92]  = 1.0000000000000002e-06;
    activation_units[92] = 0.50321666580471969;
    phase_units[92]      = 1e-12;
    is_PD[92] = 0;
    nTB[92] = 0;

    // (90):  CH3 + H2O2 <=> HO2 + CH4
    fwd_A[93]     = 24500;
    fwd_beta[93]  = 2.4700000000000002;
    fwd_Ea[93]    = 5180;
    prefactor_units[93]  = 1.0000000000000002e-06;
    activation_units[93] = 0.50321666580471969;
    phase_units[93]      = 1e-12;
    is_PD[93] = 0;
    nTB[93] = 0;

    // (91):  2 CH3 (+M) <=> C2H6 (+M)
    fwd_A[10]     = 21200000000000000;
    fwd_beta[10]  = -0.96999999999999997;
    fwd_Ea[10]    = 620;
    low_A[10]     = 1.7700000000000001e+50;
    low_beta[10]  = -9.6699999999999999;
    low_Ea[10]    = 6220;
    troe_a[10]    = 0.53249999999999997;
    troe_Tsss[10] = 151;
    troe_Ts[10]   = 1038;
    troe_Tss[10]  = 4970;
    troe_len[10]  = 4;
    prefactor_units[10]  = 1.0000000000000002e-06;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = 1e-12;
    is_PD[10] = 1;
    nTB[10] = 7;
    TB[10] = (double *) malloc(7 * sizeof(double));
    TBid[10] = (int *) malloc(7 * sizeof(int));
    TBid[10][0] = 0; TB[10][0] = 2; // H2
    TBid[10][1] = 5; TB[10][1] = 6; // H2O
    TBid[10][2] = 11; TB[10][2] = 2; // CH4
    TBid[10][3] = 12; TB[10][3] = 1.5; // CO
    TBid[10][4] = 13; TB[10][4] = 2; // CO2
    TBid[10][5] = 21; TB[10][5] = 3; // C2H6
    TBid[10][6] = 23; TB[10][6] = 0.69999999999999996; // AR

    // (92):  2 CH3 <=> H + C2H5
    fwd_A[94]     = 4990000000000;
    fwd_beta[94]  = 0.10000000000000001;
    fwd_Ea[94]    = 10600;
    prefactor_units[94]  = 1.0000000000000002e-06;
    activation_units[94] = 0.50321666580471969;
    phase_units[94]      = 1e-12;
    is_PD[94] = 0;
    nTB[94] = 0;

    // (93):  CH3 + HCO <=> CH4 + CO
    fwd_A[95]     = 26480000000000;
    fwd_beta[95]  = 0;
    fwd_Ea[95]    = 0;
    prefactor_units[95]  = 1.0000000000000002e-06;
    activation_units[95] = 0.50321666580471969;
    phase_units[95]      = 1e-12;
    is_PD[95] = 0;
    nTB[95] = 0;

    // (94):  CH3 + CH2O <=> HCO + CH4
    fwd_A[96]     = 3320;
    fwd_beta[96]  = 2.8100000000000001;
    fwd_Ea[96]    = 5860;
    prefactor_units[96]  = 1.0000000000000002e-06;
    activation_units[96] = 0.50321666580471969;
    phase_units[96]      = 1e-12;
    is_PD[96] = 0;
    nTB[96] = 0;

    // (95):  CH3 + C2H4 <=> C2H3 + CH4
    fwd_A[97]     = 227000;
    fwd_beta[97]  = 2;
    fwd_Ea[97]    = 9200;
    prefactor_units[97]  = 1.0000000000000002e-06;
    activation_units[97] = 0.50321666580471969;
    phase_units[97]      = 1e-12;
    is_PD[97] = 0;
    nTB[97] = 0;

    // (96):  CH3 + C2H6 <=> C2H5 + CH4
    fwd_A[98]     = 6140000;
    fwd_beta[98]  = 1.74;
    fwd_Ea[98]    = 10450;
    prefactor_units[98]  = 1.0000000000000002e-06;
    activation_units[98] = 0.50321666580471969;
    phase_units[98]      = 1e-12;
    is_PD[98] = 0;
    nTB[98] = 0;

    // (97):  HCO + H2O <=> H + CO + H2O
    fwd_A[99]     = 2.244e+18;
    fwd_beta[99]  = -1;
    fwd_Ea[99]    = 17000;
    prefactor_units[99]  = 1.0000000000000002e-06;
    activation_units[99] = 0.50321666580471969;
    phase_units[99]      = 1e-12;
    is_PD[99] = 0;
    nTB[99] = 0;

    // (98):  HCO + M <=> H + CO + M
    fwd_A[17]     = 1.87e+17;
    fwd_beta[17]  = -1;
    fwd_Ea[17]    = 17000;
    prefactor_units[17]  = 1.0000000000000002e-06;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = 1e-6;
    is_PD[17] = 0;
    nTB[17] = 6;
    TB[17] = (double *) malloc(6 * sizeof(double));
    TBid[17] = (int *) malloc(6 * sizeof(int));
    TBid[17][0] = 0; TB[17][0] = 2; // H2
    TBid[17][1] = 5; TB[17][1] = 0; // H2O
    TBid[17][2] = 11; TB[17][2] = 2; // CH4
    TBid[17][3] = 12; TB[17][3] = 1.5; // CO
    TBid[17][4] = 13; TB[17][4] = 2; // CO2
    TBid[17][5] = 21; TB[17][5] = 3; // C2H6

    // (99):  HCO + O2 <=> HO2 + CO
    fwd_A[100]     = 7600000000000;
    fwd_beta[100]  = 0;
    fwd_Ea[100]    = 400;
    prefactor_units[100]  = 1.0000000000000002e-06;
    activation_units[100] = 0.50321666580471969;
    phase_units[100]      = 1e-12;
    is_PD[100] = 0;
    nTB[100] = 0;

    // (100):  CH3O + O2 <=> HO2 + CH2O
    fwd_A[101]     = 4.2799999999999999e-13;
    fwd_beta[101]  = 7.5999999999999996;
    fwd_Ea[101]    = -3530;
    prefactor_units[101]  = 1.0000000000000002e-06;
    activation_units[101] = 0.50321666580471969;
    phase_units[101]      = 1e-12;
    is_PD[101] = 0;
    nTB[101] = 0;

    // (101):  C2H3 + O2 <=> HCO + CH2O
    fwd_A[102]     = 3980000000000;
    fwd_beta[102]  = 0;
    fwd_Ea[102]    = -240;
    prefactor_units[102]  = 1.0000000000000002e-06;
    activation_units[102] = 0.50321666580471969;
    phase_units[102]      = 1e-12;
    is_PD[102] = 0;
    nTB[102] = 0;

    // (102):  C2H4 (+M) <=> H2 + C2H2 (+M)
    fwd_A[11]     = 8000000000000;
    fwd_beta[11]  = 0.44;
    fwd_Ea[11]    = 88770;
    low_A[11]     = 7.0000000000000001e+50;
    low_beta[11]  = -9.3100000000000005;
    low_Ea[11]    = 99860;
    troe_a[11]    = 0.73450000000000004;
    troe_Tsss[11] = 180;
    troe_Ts[11]   = 1035;
    troe_Tss[11]  = 5417;
    troe_len[11]  = 4;
    prefactor_units[11]  = 1;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = 1e-6;
    is_PD[11] = 1;
    nTB[11] = 7;
    TB[11] = (double *) malloc(7 * sizeof(double));
    TBid[11] = (int *) malloc(7 * sizeof(int));
    TBid[11][0] = 0; TB[11][0] = 2; // H2
    TBid[11][1] = 5; TB[11][1] = 6; // H2O
    TBid[11][2] = 11; TB[11][2] = 2; // CH4
    TBid[11][3] = 12; TB[11][3] = 1.5; // CO
    TBid[11][4] = 13; TB[11][4] = 2; // CO2
    TBid[11][5] = 21; TB[11][5] = 3; // C2H6
    TBid[11][6] = 23; TB[11][6] = 0.69999999999999996; // AR

    // (103):  C2H5 + O2 <=> HO2 + C2H4
    fwd_A[103]     = 840000000000;
    fwd_beta[103]  = 0;
    fwd_Ea[103]    = 3875;
    prefactor_units[103]  = 1.0000000000000002e-06;
    activation_units[103] = 0.50321666580471969;
    phase_units[103]      = 1e-12;
    is_PD[103] = 0;
    nTB[103] = 0;

    SetAllDefaults();
}



/*A few mechanism parameters */
void CKINDX(int * iwrk, double * restrict rwrk, int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 5;
    *kk = 24;
    *ii = 104;
    *nfit = -1; /*Why do you need this anyway ?  */
}


/* strtok_r: re-entrant (threadsafe) version of strtok, helper function for tokenizing strings  */
char *strtok_r(char *s, const char *delim, char **save_ptr)
{
    char *token;

    if (s == NULL)
        s = *save_ptr;

    /* Scan leading delimiters.  */
    s += strspn (s, delim);
    if (*s == '\0')
    {
        *save_ptr = s;
        return NULL;
    }

    /* Find the end of the token.  */
    token = s;
    s = strpbrk (token, delim);
    if (s == NULL)
        /* This token finishes the string.  */
        *save_ptr = __rawmemchr (token, '\0');
    else
    {
        /* Terminate the token and make *SAVE_PTR point past it.  */
        *s = '\0';
        *save_ptr = s + 1;
    }
    return token;
}



/* ckxnum... for parsing strings  */
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double * restrict rval, int * kerr, int lenline )
{
    int n,i; /*Loop Counters */
    char cstr[1000];
    char *saveptr;
    char *p; /*String Tokens */
    /* Strip Comments  */
    for (i=0; i<lenline; ++i) {
        if (line[i]=='!') {
            break;
        }
        cstr[i] = line[i];
    }
    cstr[i] = '\0';

    p = strtok_r(cstr," ", &saveptr);
    if (!p) {
        *nval = 0;
        *kerr = 1;
        return;
    }
    for (n=0; n<*nexp; ++n) {
        rval[n] = atof(p);
        p = strtok_r(NULL, " ", &saveptr);
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
    for (i=0; i<lenkname*24; i++) {
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

    /* CH2  */
    kname[ 8*lenkname + 0 ] = 'C';
    kname[ 8*lenkname + 1 ] = 'H';
    kname[ 8*lenkname + 2 ] = '2';
    kname[ 8*lenkname + 3 ] = ' ';

    /* CH2(S)  */
    kname[ 9*lenkname + 0 ] = 'C';
    kname[ 9*lenkname + 1 ] = 'H';
    kname[ 9*lenkname + 2 ] = '2';
    kname[ 9*lenkname + 3 ] = '(';
    kname[ 9*lenkname + 4 ] = 'S';
    kname[ 9*lenkname + 5 ] = ')';
    kname[ 9*lenkname + 6 ] = ' ';

    /* CH3  */
    kname[ 10*lenkname + 0 ] = 'C';
    kname[ 10*lenkname + 1 ] = 'H';
    kname[ 10*lenkname + 2 ] = '3';
    kname[ 10*lenkname + 3 ] = ' ';

    /* CH4  */
    kname[ 11*lenkname + 0 ] = 'C';
    kname[ 11*lenkname + 1 ] = 'H';
    kname[ 11*lenkname + 2 ] = '4';
    kname[ 11*lenkname + 3 ] = ' ';

    /* CO  */
    kname[ 12*lenkname + 0 ] = 'C';
    kname[ 12*lenkname + 1 ] = 'O';
    kname[ 12*lenkname + 2 ] = ' ';

    /* CO2  */
    kname[ 13*lenkname + 0 ] = 'C';
    kname[ 13*lenkname + 1 ] = 'O';
    kname[ 13*lenkname + 2 ] = '2';
    kname[ 13*lenkname + 3 ] = ' ';

    /* HCO  */
    kname[ 14*lenkname + 0 ] = 'H';
    kname[ 14*lenkname + 1 ] = 'C';
    kname[ 14*lenkname + 2 ] = 'O';
    kname[ 14*lenkname + 3 ] = ' ';

    /* CH2O  */
    kname[ 15*lenkname + 0 ] = 'C';
    kname[ 15*lenkname + 1 ] = 'H';
    kname[ 15*lenkname + 2 ] = '2';
    kname[ 15*lenkname + 3 ] = 'O';
    kname[ 15*lenkname + 4 ] = ' ';

    /* CH3O  */
    kname[ 16*lenkname + 0 ] = 'C';
    kname[ 16*lenkname + 1 ] = 'H';
    kname[ 16*lenkname + 2 ] = '3';
    kname[ 16*lenkname + 3 ] = 'O';
    kname[ 16*lenkname + 4 ] = ' ';

    /* C2H2  */
    kname[ 17*lenkname + 0 ] = 'C';
    kname[ 17*lenkname + 1 ] = '2';
    kname[ 17*lenkname + 2 ] = 'H';
    kname[ 17*lenkname + 3 ] = '2';
    kname[ 17*lenkname + 4 ] = ' ';

    /* C2H3  */
    kname[ 18*lenkname + 0 ] = 'C';
    kname[ 18*lenkname + 1 ] = '2';
    kname[ 18*lenkname + 2 ] = 'H';
    kname[ 18*lenkname + 3 ] = '3';
    kname[ 18*lenkname + 4 ] = ' ';

    /* C2H4  */
    kname[ 19*lenkname + 0 ] = 'C';
    kname[ 19*lenkname + 1 ] = '2';
    kname[ 19*lenkname + 2 ] = 'H';
    kname[ 19*lenkname + 3 ] = '4';
    kname[ 19*lenkname + 4 ] = ' ';

    /* C2H5  */
    kname[ 20*lenkname + 0 ] = 'C';
    kname[ 20*lenkname + 1 ] = '2';
    kname[ 20*lenkname + 2 ] = 'H';
    kname[ 20*lenkname + 3 ] = '5';
    kname[ 20*lenkname + 4 ] = ' ';

    /* C2H6  */
    kname[ 21*lenkname + 0 ] = 'C';
    kname[ 21*lenkname + 1 ] = '2';
    kname[ 21*lenkname + 2 ] = 'H';
    kname[ 21*lenkname + 3 ] = '6';
    kname[ 21*lenkname + 4 ] = ' ';

    /* N2  */
    kname[ 22*lenkname + 0 ] = 'N';
    kname[ 22*lenkname + 1 ] = '2';
    kname[ 22*lenkname + 2 ] = ' ';

    /* AR  */
    kname[ 23*lenkname + 0 ] = 'A';
    kname[ 23*lenkname + 1 ] = 'R';
    kname[ 23*lenkname + 2 ] = ' ';

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
    XW += x[8]*14.027090; /*CH2 */
    XW += x[9]*14.027090; /*CH2(S) */
    XW += x[10]*15.035060; /*CH3 */
    XW += x[11]*16.043030; /*CH4 */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*44.009950; /*CO2 */
    XW += x[14]*29.018520; /*HCO */
    XW += x[15]*30.026490; /*CH2O */
    XW += x[16]*31.034460; /*CH3O */
    XW += x[17]*26.038240; /*C2H2 */
    XW += x[18]*27.046210; /*C2H3 */
    XW += x[19]*28.054180; /*C2H4 */
    XW += x[20]*29.062150; /*C2H5 */
    XW += x[21]*30.070120; /*C2H6 */
    XW += x[22]*28.013400; /*N2 */
    XW += x[23]*39.948000; /*AR */
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
    YOW += y[8]*imw[8]; /*CH2 */
    YOW += y[9]*imw[9]; /*CH2(S) */
    YOW += y[10]*imw[10]; /*CH3 */
    YOW += y[11]*imw[11]; /*CH4 */
    YOW += y[12]*imw[12]; /*CO */
    YOW += y[13]*imw[13]; /*CO2 */
    YOW += y[14]*imw[14]; /*HCO */
    YOW += y[15]*imw[15]; /*CH2O */
    YOW += y[16]*imw[16]; /*CH3O */
    YOW += y[17]*imw[17]; /*C2H2 */
    YOW += y[18]*imw[18]; /*C2H3 */
    YOW += y[19]*imw[19]; /*C2H4 */
    YOW += y[20]*imw[20]; /*C2H5 */
    YOW += y[21]*imw[21]; /*C2H6 */
    YOW += y[22]*imw[22]; /*N2 */
    YOW += y[23]*imw[23]; /*AR */
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

    for (int n=0; n<24; n++) {
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
    W += c[8]*14.027090; /*CH2 */
    W += c[9]*14.027090; /*CH2(S) */
    W += c[10]*15.035060; /*CH3 */
    W += c[11]*16.043030; /*CH4 */
    W += c[12]*28.010550; /*CO */
    W += c[13]*44.009950; /*CO2 */
    W += c[14]*29.018520; /*HCO */
    W += c[15]*30.026490; /*CH2O */
    W += c[16]*31.034460; /*CH3O */
    W += c[17]*26.038240; /*C2H2 */
    W += c[18]*27.046210; /*C2H3 */
    W += c[19]*28.054180; /*C2H4 */
    W += c[20]*29.062150; /*C2H5 */
    W += c[21]*30.070120; /*C2H6 */
    W += c[22]*28.013400; /*N2 */
    W += c[23]*39.948000; /*AR */

    for (id = 0; id < 24; ++id) {
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
    XW += x[8]*14.027090; /*CH2 */
    XW += x[9]*14.027090; /*CH2(S) */
    XW += x[10]*15.035060; /*CH3 */
    XW += x[11]*16.043030; /*CH4 */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*44.009950; /*CO2 */
    XW += x[14]*29.018520; /*HCO */
    XW += x[15]*30.026490; /*CH2O */
    XW += x[16]*31.034460; /*CH3O */
    XW += x[17]*26.038240; /*C2H2 */
    XW += x[18]*27.046210; /*C2H3 */
    XW += x[19]*28.054180; /*C2H4 */
    XW += x[20]*29.062150; /*C2H5 */
    XW += x[21]*30.070120; /*C2H6 */
    XW += x[22]*28.013400; /*N2 */
    XW += x[23]*39.948000; /*AR */
    *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    double YOW = 0;
    double tmp[24];

    for (int i = 0; i < 24; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 24; i++)
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
    W += c[8]*14.027090; /*CH2 */
    W += c[9]*14.027090; /*CH2(S) */
    W += c[10]*15.035060; /*CH3 */
    W += c[11]*16.043030; /*CH4 */
    W += c[12]*28.010550; /*CO */
    W += c[13]*44.009950; /*CO2 */
    W += c[14]*29.018520; /*HCO */
    W += c[15]*30.026490; /*CH2O */
    W += c[16]*31.034460; /*CH3O */
    W += c[17]*26.038240; /*C2H2 */
    W += c[18]*27.046210; /*C2H3 */
    W += c[19]*28.054180; /*C2H4 */
    W += c[20]*29.062150; /*C2H5 */
    W += c[21]*30.070120; /*C2H6 */
    W += c[22]*28.013400; /*N2 */
    W += c[23]*39.948000; /*AR */

    for (id = 0; id < 24; ++id) {
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
    double tmp[24];

    for (int i = 0; i < 24; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 24; i++)
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
    XW += x[8]*14.027090; /*CH2 */
    XW += x[9]*14.027090; /*CH2(S) */
    XW += x[10]*15.035060; /*CH3 */
    XW += x[11]*16.043030; /*CH4 */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*44.009950; /*CO2 */
    XW += x[14]*29.018520; /*HCO */
    XW += x[15]*30.026490; /*CH2O */
    XW += x[16]*31.034460; /*CH3O */
    XW += x[17]*26.038240; /*C2H2 */
    XW += x[18]*27.046210; /*C2H3 */
    XW += x[19]*28.054180; /*C2H4 */
    XW += x[20]*29.062150; /*C2H5 */
    XW += x[21]*30.070120; /*C2H6 */
    XW += x[22]*28.013400; /*N2 */
    XW += x[23]*39.948000; /*AR */
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
    W += c[8]*14.027090; /*CH2 */
    W += c[9]*14.027090; /*CH2(S) */
    W += c[10]*15.035060; /*CH3 */
    W += c[11]*16.043030; /*CH4 */
    W += c[12]*28.010550; /*CO */
    W += c[13]*44.009950; /*CO2 */
    W += c[14]*29.018520; /*HCO */
    W += c[15]*30.026490; /*CH2O */
    W += c[16]*31.034460; /*CH3O */
    W += c[17]*26.038240; /*C2H2 */
    W += c[18]*27.046210; /*C2H3 */
    W += c[19]*28.054180; /*C2H4 */
    W += c[20]*29.062150; /*C2H5 */
    W += c[21]*30.070120; /*C2H6 */
    W += c[22]*28.013400; /*N2 */
    W += c[23]*39.948000; /*AR */

    for (id = 0; id < 24; ++id) {
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
    double tmp[24];

    for (int i = 0; i < 24; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 24; i++)
    {
        YOW += tmp[i];
    }

    double YOWINV = 1.0/YOW;

    for (int i = 0; i < 24; i++)
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

    for (int n=0; n<24; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<24; n++) {
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
    for (int i = 0; i < 24; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 24; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 24; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c)
{
    for (int i = 0; i < 24; i++)
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
    XW += x[8]*14.027090; /*CH2 */
    XW += x[9]*14.027090; /*CH2(S) */
    XW += x[10]*15.035060; /*CH3 */
    XW += x[11]*16.043030; /*CH4 */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*44.009950; /*CO2 */
    XW += x[14]*29.018520; /*HCO */
    XW += x[15]*30.026490; /*CH2O */
    XW += x[16]*31.034460; /*CH3O */
    XW += x[17]*26.038240; /*C2H2 */
    XW += x[18]*27.046210; /*C2H3 */
    XW += x[19]*28.054180; /*C2H4 */
    XW += x[20]*29.062150; /*C2H5 */
    XW += x[21]*30.070120; /*C2H6 */
    XW += x[22]*28.013400; /*N2 */
    XW += x[23]*39.948000; /*AR */
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
    y[8] = x[8]*14.027090*XWinv; 
    y[9] = x[9]*14.027090*XWinv; 
    y[10] = x[10]*15.035060*XWinv; 
    y[11] = x[11]*16.043030*XWinv; 
    y[12] = x[12]*28.010550*XWinv; 
    y[13] = x[13]*44.009950*XWinv; 
    y[14] = x[14]*29.018520*XWinv; 
    y[15] = x[15]*30.026490*XWinv; 
    y[16] = x[16]*31.034460*XWinv; 
    y[17] = x[17]*26.038240*XWinv; 
    y[18] = x[18]*27.046210*XWinv; 
    y[19] = x[19]*28.054180*XWinv; 
    y[20] = x[20]*29.062150*XWinv; 
    y[21] = x[21]*30.070120*XWinv; 
    y[22] = x[22]*28.013400*XWinv; 
    y[23] = x[23]*39.948000*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 24; ++id) {
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
    XW += x[8]*14.027090; /*CH2 */
    XW += x[9]*14.027090; /*CH2(S) */
    XW += x[10]*15.035060; /*CH3 */
    XW += x[11]*16.043030; /*CH4 */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*44.009950; /*CO2 */
    XW += x[14]*29.018520; /*HCO */
    XW += x[15]*30.026490; /*CH2O */
    XW += x[16]*31.034460; /*CH3O */
    XW += x[17]*26.038240; /*C2H2 */
    XW += x[18]*27.046210; /*C2H3 */
    XW += x[19]*28.054180; /*C2H4 */
    XW += x[20]*29.062150; /*C2H5 */
    XW += x[21]*30.070120; /*C2H6 */
    XW += x[22]*28.013400; /*N2 */
    XW += x[23]*39.948000; /*AR */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 24; ++id) {
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
    for (id = 0; id < 24; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    double sumCinv = 1.0/sumC;
    for (id = 0; id < 24; ++id) {
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
    CW += c[8]*14.027090; /*CH2 */
    CW += c[9]*14.027090; /*CH2(S) */
    CW += c[10]*15.035060; /*CH3 */
    CW += c[11]*16.043030; /*CH4 */
    CW += c[12]*28.010550; /*CO */
    CW += c[13]*44.009950; /*CO2 */
    CW += c[14]*29.018520; /*HCO */
    CW += c[15]*30.026490; /*CH2O */
    CW += c[16]*31.034460; /*CH3O */
    CW += c[17]*26.038240; /*C2H2 */
    CW += c[18]*27.046210; /*C2H3 */
    CW += c[19]*28.054180; /*C2H4 */
    CW += c[20]*29.062150; /*C2H5 */
    CW += c[21]*30.070120; /*C2H6 */
    CW += c[22]*28.013400; /*N2 */
    CW += c[23]*39.948000; /*AR */
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
    y[8] = c[8]*14.027090*CWinv; 
    y[9] = c[9]*14.027090*CWinv; 
    y[10] = c[10]*15.035060*CWinv; 
    y[11] = c[11]*16.043030*CWinv; 
    y[12] = c[12]*28.010550*CWinv; 
    y[13] = c[13]*44.009950*CWinv; 
    y[14] = c[14]*29.018520*CWinv; 
    y[15] = c[15]*30.026490*CWinv; 
    y[16] = c[16]*31.034460*CWinv; 
    y[17] = c[17]*26.038240*CWinv; 
    y[18] = c[18]*27.046210*CWinv; 
    y[19] = c[19]*28.054180*CWinv; 
    y[20] = c[20]*29.062150*CWinv; 
    y[21] = c[21]*30.070120*CWinv; 
    y[22] = c[22]*28.013400*CWinv; 
    y[23] = c[23]*39.948000*CWinv; 

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
    for (id = 0; id < 24; ++id) {
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
    for (id = 0; id < 24; ++id) {
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
    for (id = 0; id < 24; ++id) {
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
    for (id = 0; id < 24; ++id) {
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
    for (id = 0; id < 24; ++id) {
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
    for (id = 0; id < 24; ++id) {
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
    for (id = 0; id < 24; ++id) {
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
    cvms[8] *= 5.927466067445207e+06; /*CH2 */
    cvms[9] *= 5.927466067445207e+06; /*CH2(S) */
    cvms[10] *= 5.530081023953346e+06; /*CH3 */
    cvms[11] *= 5.182630712527496e+06; /*CH4 */
    cvms[12] *= 2.968349425484326e+06; /*CO */
    cvms[13] *= 1.889234139098090e+06; /*CO2 */
    cvms[14] *= 2.865242610581105e+06; /*HCO */
    cvms[15] *= 2.769058254894261e+06; /*CH2O */
    cvms[16] *= 2.679121853578248e+06; /*CH3O */
    cvms[17] *= 3.193192012977835e+06; /*C2H2 */
    cvms[18] *= 3.074186734481467e+06; /*C2H3 */
    cvms[19] *= 2.963733033722604e+06; /*C2H4 */
    cvms[20] *= 2.860941121011349e+06; /*C2H5 */
    cvms[21] *= 2.765040511976673e+06; /*C2H6 */
    cvms[22] *= 2.968047434442088e+06; /*N2 */
    cvms[23] *= 2.081333233203164e+06; /*AR */
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
    cpms[8] *= 5.927466067445207e+06; /*CH2 */
    cpms[9] *= 5.927466067445207e+06; /*CH2(S) */
    cpms[10] *= 5.530081023953346e+06; /*CH3 */
    cpms[11] *= 5.182630712527496e+06; /*CH4 */
    cpms[12] *= 2.968349425484326e+06; /*CO */
    cpms[13] *= 1.889234139098090e+06; /*CO2 */
    cpms[14] *= 2.865242610581105e+06; /*HCO */
    cpms[15] *= 2.769058254894261e+06; /*CH2O */
    cpms[16] *= 2.679121853578248e+06; /*CH3O */
    cpms[17] *= 3.193192012977835e+06; /*C2H2 */
    cpms[18] *= 3.074186734481467e+06; /*C2H3 */
    cpms[19] *= 2.963733033722604e+06; /*C2H4 */
    cpms[20] *= 2.860941121011349e+06; /*C2H5 */
    cpms[21] *= 2.765040511976673e+06; /*C2H6 */
    cpms[22] *= 2.968047434442088e+06; /*N2 */
    cpms[23] *= 2.081333233203164e+06; /*AR */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    for (int i = 0; i < 24; i++)
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
    for (int i = 0; i < 24; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hms)
{
    double tc[5], h[24];

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
        hms[14*(*np)+i] = h[14];
        hms[15*(*np)+i] = h[15];
        hms[16*(*np)+i] = h[16];
        hms[17*(*np)+i] = h[17];
        hms[18*(*np)+i] = h[18];
        hms[19*(*np)+i] = h[19];
        hms[20*(*np)+i] = h[20];
        hms[21*(*np)+i] = h[21];
        hms[22*(*np)+i] = h[22];
        hms[23*(*np)+i] = h[23];
    }

    for (int n=0; n<24; n++) {
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
    for (int i = 0; i < 24; i++)
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
    for (int i = 0; i < 24; i++)
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
    sms[8] *= 5.927466067445207e+06; /*CH2 */
    sms[9] *= 5.927466067445207e+06; /*CH2(S) */
    sms[10] *= 5.530081023953346e+06; /*CH3 */
    sms[11] *= 5.182630712527496e+06; /*CH4 */
    sms[12] *= 2.968349425484326e+06; /*CO */
    sms[13] *= 1.889234139098090e+06; /*CO2 */
    sms[14] *= 2.865242610581105e+06; /*HCO */
    sms[15] *= 2.769058254894261e+06; /*CH2O */
    sms[16] *= 2.679121853578248e+06; /*CH3O */
    sms[17] *= 3.193192012977835e+06; /*C2H2 */
    sms[18] *= 3.074186734481467e+06; /*C2H3 */
    sms[19] *= 2.963733033722604e+06; /*C2H4 */
    sms[20] *= 2.860941121011349e+06; /*C2H5 */
    sms[21] *= 2.765040511976673e+06; /*C2H6 */
    sms[22] *= 2.968047434442088e+06; /*N2 */
    sms[23] *= 2.081333233203164e+06; /*AR */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[24]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 24; ++id) {
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
    double cpor[24], tresult[24]; /* temporary storage */
    cp_R(cpor, tc);
    for (int i = 0; i < 24; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 24; i++)
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
    double cvor[24]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 24; ++id) {
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
    double cvor[24]; /* temporary storage */
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
    result += cvor[8]*y[8]*imw[8]; /*CH2 */
    result += cvor[9]*y[9]*imw[9]; /*CH2(S) */
    result += cvor[10]*y[10]*imw[10]; /*CH3 */
    result += cvor[11]*y[11]*imw[11]; /*CH4 */
    result += cvor[12]*y[12]*imw[12]; /*CO */
    result += cvor[13]*y[13]*imw[13]; /*CO2 */
    result += cvor[14]*y[14]*imw[14]; /*HCO */
    result += cvor[15]*y[15]*imw[15]; /*CH2O */
    result += cvor[16]*y[16]*imw[16]; /*CH3O */
    result += cvor[17]*y[17]*imw[17]; /*C2H2 */
    result += cvor[18]*y[18]*imw[18]; /*C2H3 */
    result += cvor[19]*y[19]*imw[19]; /*C2H4 */
    result += cvor[20]*y[20]*imw[20]; /*C2H5 */
    result += cvor[21]*y[21]*imw[21]; /*C2H6 */
    result += cvor[22]*y[22]*imw[22]; /*N2 */
    result += cvor[23]*y[23]*imw[23]; /*AR */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[24]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 24; ++id) {
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
    double hml[24], tmp[24]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    int id;
    for (id = 0; id < 24; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 24; ++id) {
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
    double uml[24]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 24; ++id) {
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
    double ums[24]; /* temporary energy array */
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
    result += y[8]*ums[8]*imw[8]; /*CH2 */
    result += y[9]*ums[9]*imw[9]; /*CH2(S) */
    result += y[10]*ums[10]*imw[10]; /*CH3 */
    result += y[11]*ums[11]*imw[11]; /*CH4 */
    result += y[12]*ums[12]*imw[12]; /*CO */
    result += y[13]*ums[13]*imw[13]; /*CO2 */
    result += y[14]*ums[14]*imw[14]; /*HCO */
    result += y[15]*ums[15]*imw[15]; /*CH2O */
    result += y[16]*ums[16]*imw[16]; /*CH3O */
    result += y[17]*ums[17]*imw[17]; /*C2H2 */
    result += y[18]*ums[18]*imw[18]; /*C2H3 */
    result += y[19]*ums[19]*imw[19]; /*C2H4 */
    result += y[20]*ums[20]*imw[20]; /*C2H5 */
    result += y[21]*ums[21]*imw[21]; /*C2H6 */
    result += y[22]*ums[22]*imw[22]; /*N2 */
    result += y[23]*ums[23]*imw[23]; /*AR */

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
    double sor[24]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 24; ++id) {
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
    double sor[24]; /* temporary storage */
    double x[24]; /* need a ytx conversion */
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
    YOW += y[8]*imw[8]; /*CH2 */
    YOW += y[9]*imw[9]; /*CH2(S) */
    YOW += y[10]*imw[10]; /*CH3 */
    YOW += y[11]*imw[11]; /*CH4 */
    YOW += y[12]*imw[12]; /*CO */
    YOW += y[13]*imw[13]; /*CO2 */
    YOW += y[14]*imw[14]; /*HCO */
    YOW += y[15]*imw[15]; /*CH2O */
    YOW += y[16]*imw[16]; /*CH3O */
    YOW += y[17]*imw[17]; /*C2H2 */
    YOW += y[18]*imw[18]; /*C2H3 */
    YOW += y[19]*imw[19]; /*C2H4 */
    YOW += y[20]*imw[20]; /*C2H5 */
    YOW += y[21]*imw[21]; /*C2H6 */
    YOW += y[22]*imw[22]; /*N2 */
    YOW += y[23]*imw[23]; /*AR */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(14.027090*YOW); 
    x[9] = y[9]/(14.027090*YOW); 
    x[10] = y[10]/(15.035060*YOW); 
    x[11] = y[11]/(16.043030*YOW); 
    x[12] = y[12]/(28.010550*YOW); 
    x[13] = y[13]/(44.009950*YOW); 
    x[14] = y[14]/(29.018520*YOW); 
    x[15] = y[15]/(30.026490*YOW); 
    x[16] = y[16]/(31.034460*YOW); 
    x[17] = y[17]/(26.038240*YOW); 
    x[18] = y[18]/(27.046210*YOW); 
    x[19] = y[19]/(28.054180*YOW); 
    x[20] = y[20]/(29.062150*YOW); 
    x[21] = y[21]/(30.070120*YOW); 
    x[22] = y[22]/(28.013400*YOW); 
    x[23] = y[23]/(39.948000*YOW); 
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
    double gort[24]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 24; ++id) {
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
    double gort[24]; /* temporary storage */
    double x[24]; /* need a ytx conversion */
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
    YOW += y[8]*imw[8]; /*CH2 */
    YOW += y[9]*imw[9]; /*CH2(S) */
    YOW += y[10]*imw[10]; /*CH3 */
    YOW += y[11]*imw[11]; /*CH4 */
    YOW += y[12]*imw[12]; /*CO */
    YOW += y[13]*imw[13]; /*CO2 */
    YOW += y[14]*imw[14]; /*HCO */
    YOW += y[15]*imw[15]; /*CH2O */
    YOW += y[16]*imw[16]; /*CH3O */
    YOW += y[17]*imw[17]; /*C2H2 */
    YOW += y[18]*imw[18]; /*C2H3 */
    YOW += y[19]*imw[19]; /*C2H4 */
    YOW += y[20]*imw[20]; /*C2H5 */
    YOW += y[21]*imw[21]; /*C2H6 */
    YOW += y[22]*imw[22]; /*N2 */
    YOW += y[23]*imw[23]; /*AR */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(14.027090*YOW); 
    x[9] = y[9]/(14.027090*YOW); 
    x[10] = y[10]/(15.035060*YOW); 
    x[11] = y[11]/(16.043030*YOW); 
    x[12] = y[12]/(28.010550*YOW); 
    x[13] = y[13]/(44.009950*YOW); 
    x[14] = y[14]/(29.018520*YOW); 
    x[15] = y[15]/(30.026490*YOW); 
    x[16] = y[16]/(31.034460*YOW); 
    x[17] = y[17]/(26.038240*YOW); 
    x[18] = y[18]/(27.046210*YOW); 
    x[19] = y[19]/(28.054180*YOW); 
    x[20] = y[20]/(29.062150*YOW); 
    x[21] = y[21]/(30.070120*YOW); 
    x[22] = y[22]/(28.013400*YOW); 
    x[23] = y[23]/(39.948000*YOW); 
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
    double aort[24]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 24; ++id) {
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
    double aort[24]; /* temporary storage */
    double x[24]; /* need a ytx conversion */
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
    YOW += y[8]*imw[8]; /*CH2 */
    YOW += y[9]*imw[9]; /*CH2(S) */
    YOW += y[10]*imw[10]; /*CH3 */
    YOW += y[11]*imw[11]; /*CH4 */
    YOW += y[12]*imw[12]; /*CO */
    YOW += y[13]*imw[13]; /*CO2 */
    YOW += y[14]*imw[14]; /*HCO */
    YOW += y[15]*imw[15]; /*CH2O */
    YOW += y[16]*imw[16]; /*CH3O */
    YOW += y[17]*imw[17]; /*C2H2 */
    YOW += y[18]*imw[18]; /*C2H3 */
    YOW += y[19]*imw[19]; /*C2H4 */
    YOW += y[20]*imw[20]; /*C2H5 */
    YOW += y[21]*imw[21]; /*C2H6 */
    YOW += y[22]*imw[22]; /*N2 */
    YOW += y[23]*imw[23]; /*AR */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(14.027090*YOW); 
    x[9] = y[9]/(14.027090*YOW); 
    x[10] = y[10]/(15.035060*YOW); 
    x[11] = y[11]/(16.043030*YOW); 
    x[12] = y[12]/(28.010550*YOW); 
    x[13] = y[13]/(44.009950*YOW); 
    x[14] = y[14]/(29.018520*YOW); 
    x[15] = y[15]/(30.026490*YOW); 
    x[16] = y[16]/(31.034460*YOW); 
    x[17] = y[17]/(26.038240*YOW); 
    x[18] = y[18]/(27.046210*YOW); 
    x[19] = y[19]/(28.054180*YOW); 
    x[20] = y[20]/(29.062150*YOW); 
    x[21] = y[21]/(30.070120*YOW); 
    x[22] = y[22]/(28.013400*YOW); 
    x[23] = y[23]/(39.948000*YOW); 
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
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 24; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 24; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[24]; /*temporary storage */
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
    YOW += y[8]*imw[8]; /*CH2 */
    YOW += y[9]*imw[9]; /*CH2(S) */
    YOW += y[10]*imw[10]; /*CH3 */
    YOW += y[11]*imw[11]; /*CH4 */
    YOW += y[12]*imw[12]; /*CO */
    YOW += y[13]*imw[13]; /*CO2 */
    YOW += y[14]*imw[14]; /*HCO */
    YOW += y[15]*imw[15]; /*CH2O */
    YOW += y[16]*imw[16]; /*CH3O */
    YOW += y[17]*imw[17]; /*C2H2 */
    YOW += y[18]*imw[18]; /*C2H3 */
    YOW += y[19]*imw[19]; /*C2H4 */
    YOW += y[20]*imw[20]; /*C2H5 */
    YOW += y[21]*imw[21]; /*C2H6 */
    YOW += y[22]*imw[22]; /*N2 */
    YOW += y[23]*imw[23]; /*AR */
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
    c[14] = PWORT * y[14]*imw[14]; 
    c[15] = PWORT * y[15]*imw[15]; 
    c[16] = PWORT * y[16]*imw[16]; 
    c[17] = PWORT * y[17]*imw[17]; 
    c[18] = PWORT * y[18]*imw[18]; 
    c[19] = PWORT * y[19]*imw[19]; 
    c[20] = PWORT * y[20]*imw[20]; 
    c[21] = PWORT * y[21]*imw[21]; 
    c[22] = PWORT * y[22]*imw[22]; 
    c[23] = PWORT * y[23]*imw[23]; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 24; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[24]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 24; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 24; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[24]; /*temporary storage */
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
    c[14] = 1e6 * (*rho) * y[14]*imw[14]; 
    c[15] = 1e6 * (*rho) * y[15]*imw[15]; 
    c[16] = 1e6 * (*rho) * y[16]*imw[16]; 
    c[17] = 1e6 * (*rho) * y[17]*imw[17]; 
    c[18] = 1e6 * (*rho) * y[18]*imw[18]; 
    c[19] = 1e6 * (*rho) * y[19]*imw[19]; 
    c[20] = 1e6 * (*rho) * y[20]*imw[20]; 
    c[21] = 1e6 * (*rho) * y[21]*imw[21]; 
    c[22] = 1e6 * (*rho) * y[22]*imw[22]; 
    c[23] = 1e6 * (*rho) * y[23]*imw[23]; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 24; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int * restrict np, double * restrict rho, double * restrict T,
	    double * restrict y, int * restrict iwrk, double * restrict rwrk,
	    double * restrict wdot)
{
    double c[24*(*np)]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<24; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<24*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[24]; /*temporary storage */
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
    XW += x[8]*14.027090; /*CH2 */
    XW += x[9]*14.027090; /*CH2(S) */
    XW += x[10]*15.035060; /*CH3 */
    XW += x[11]*16.043030; /*CH4 */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*44.009950; /*CO2 */
    XW += x[14]*29.018520; /*HCO */
    XW += x[15]*30.026490; /*CH2O */
    XW += x[16]*31.034460; /*CH3O */
    XW += x[17]*26.038240; /*C2H2 */
    XW += x[18]*27.046210; /*C2H3 */
    XW += x[19]*28.054180; /*C2H4 */
    XW += x[20]*29.062150; /*C2H5 */
    XW += x[21]*30.070120; /*C2H6 */
    XW += x[22]*28.013400; /*N2 */
    XW += x[23]*39.948000; /*AR */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 24; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 24; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 24; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 24; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 104; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict q_f, double * restrict q_r)
{
    int id; /*loop counter */
    double c[24]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 24; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 104; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[24]; /*temporary storage */
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
    YOW += y[8]*imw[8]; /*CH2 */
    YOW += y[9]*imw[9]; /*CH2(S) */
    YOW += y[10]*imw[10]; /*CH3 */
    YOW += y[11]*imw[11]; /*CH4 */
    YOW += y[12]*imw[12]; /*CO */
    YOW += y[13]*imw[13]; /*CO2 */
    YOW += y[14]*imw[14]; /*HCO */
    YOW += y[15]*imw[15]; /*CH2O */
    YOW += y[16]*imw[16]; /*CH3O */
    YOW += y[17]*imw[17]; /*C2H2 */
    YOW += y[18]*imw[18]; /*C2H3 */
    YOW += y[19]*imw[19]; /*C2H4 */
    YOW += y[20]*imw[20]; /*C2H5 */
    YOW += y[21]*imw[21]; /*C2H6 */
    YOW += y[22]*imw[22]; /*N2 */
    YOW += y[23]*imw[23]; /*AR */
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
    c[14] = PWORT * y[14]*imw[14]; 
    c[15] = PWORT * y[15]*imw[15]; 
    c[16] = PWORT * y[16]*imw[16]; 
    c[17] = PWORT * y[17]*imw[17]; 
    c[18] = PWORT * y[18]*imw[18]; 
    c[19] = PWORT * y[19]*imw[19]; 
    c[20] = PWORT * y[20]*imw[20]; 
    c[21] = PWORT * y[21]*imw[21]; 
    c[22] = PWORT * y[22]*imw[22]; 
    c[23] = PWORT * y[23]*imw[23]; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 104; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[24]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 24; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 104; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[24]; /*temporary storage */
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
    c[14] = 1e6 * (*rho) * y[14]*imw[14]; 
    c[15] = 1e6 * (*rho) * y[15]*imw[15]; 
    c[16] = 1e6 * (*rho) * y[16]*imw[16]; 
    c[17] = 1e6 * (*rho) * y[17]*imw[17]; 
    c[18] = 1e6 * (*rho) * y[18]*imw[18]; 
    c[19] = 1e6 * (*rho) * y[19]*imw[19]; 
    c[20] = 1e6 * (*rho) * y[20]*imw[20]; 
    c[21] = 1e6 * (*rho) * y[21]*imw[21]; 
    c[22] = 1e6 * (*rho) * y[22]*imw[22]; 
    c[23] = 1e6 * (*rho) * y[23]*imw[23]; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 104; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[24]; /*temporary storage */
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
    XW += x[8]*14.027090; /*CH2 */
    XW += x[9]*14.027090; /*CH2(S) */
    XW += x[10]*15.035060; /*CH3 */
    XW += x[11]*16.043030; /*CH4 */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*44.009950; /*CO2 */
    XW += x[14]*29.018520; /*HCO */
    XW += x[15]*30.026490; /*CH2O */
    XW += x[16]*31.034460; /*CH3O */
    XW += x[17]*26.038240; /*C2H2 */
    XW += x[18]*27.046210; /*C2H3 */
    XW += x[19]*28.054180; /*C2H4 */
    XW += x[20]*29.062150; /*C2H5 */
    XW += x[21]*30.070120; /*C2H6 */
    XW += x[22]*28.013400; /*N2 */
    XW += x[23]*39.948000; /*AR */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 24; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 104; ++id) {
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
    for (id = 0; id < 24 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
    nuki[ 1 * kd + 0 ] += -1 ;
    nuki[ 8 * kd + 0 ] += -1 ;
    nuki[ 10 * kd + 0 ] += +1 ;

    /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
    nuki[ 1 * kd + 1 ] += -1 ;
    nuki[ 10 * kd + 1 ] += -1 ;
    nuki[ 11 * kd + 1 ] += +1 ;

    /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
    nuki[ 1 * kd + 2 ] += -1 ;
    nuki[ 14 * kd + 2 ] += -1 ;
    nuki[ 15 * kd + 2 ] += +1 ;

    /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
    nuki[ 1 * kd + 3 ] += -1 ;
    nuki[ 15 * kd + 3 ] += -1 ;
    nuki[ 16 * kd + 3 ] += +1 ;

    /*reaction 5: H + C2H2 (+M) <=> C2H3 (+M) */
    nuki[ 1 * kd + 4 ] += -1 ;
    nuki[ 17 * kd + 4 ] += -1 ;
    nuki[ 18 * kd + 4 ] += +1 ;

    /*reaction 6: H + C2H3 (+M) <=> C2H4 (+M) */
    nuki[ 1 * kd + 5 ] += -1 ;
    nuki[ 18 * kd + 5 ] += -1 ;
    nuki[ 19 * kd + 5 ] += +1 ;

    /*reaction 7: H + C2H4 (+M) <=> C2H5 (+M) */
    nuki[ 1 * kd + 6 ] += -1 ;
    nuki[ 19 * kd + 6 ] += -1 ;
    nuki[ 20 * kd + 6 ] += +1 ;

    /*reaction 8: H + C2H5 (+M) <=> C2H6 (+M) */
    nuki[ 1 * kd + 7 ] += -1 ;
    nuki[ 20 * kd + 7 ] += -1 ;
    nuki[ 21 * kd + 7 ] += +1 ;

    /*reaction 9: H2 + CO (+M) <=> CH2O (+M) */
    nuki[ 0 * kd + 8 ] += -1 ;
    nuki[ 12 * kd + 8 ] += -1 ;
    nuki[ 15 * kd + 8 ] += +1 ;

    /*reaction 10: 2 OH (+M) <=> H2O2 (+M) */
    nuki[ 4 * kd + 9 ] += -2 ;
    nuki[ 7 * kd + 9 ] += +1 ;

    /*reaction 11: 2 CH3 (+M) <=> C2H6 (+M) */
    nuki[ 10 * kd + 10 ] += -2 ;
    nuki[ 21 * kd + 10 ] += +1 ;

    /*reaction 12: C2H4 (+M) <=> H2 + C2H2 (+M) */
    nuki[ 19 * kd + 11 ] += -1 ;
    nuki[ 0 * kd + 11 ] += +1 ;
    nuki[ 17 * kd + 11 ] += +1 ;

    /*reaction 13: O + H + M <=> OH + M */
    nuki[ 2 * kd + 12 ] += -1 ;
    nuki[ 1 * kd + 12 ] += -1 ;
    nuki[ 4 * kd + 12 ] += +1 ;

    /*reaction 14: O + CO + M <=> CO2 + M */
    nuki[ 2 * kd + 13 ] += -1 ;
    nuki[ 12 * kd + 13 ] += -1 ;
    nuki[ 13 * kd + 13 ] += +1 ;

    /*reaction 15: H + O2 + M <=> HO2 + M */
    nuki[ 1 * kd + 14 ] += -1 ;
    nuki[ 3 * kd + 14 ] += -1 ;
    nuki[ 6 * kd + 14 ] += +1 ;

    /*reaction 16: 2 H + M <=> H2 + M */
    nuki[ 1 * kd + 15 ] += -2 ;
    nuki[ 0 * kd + 15 ] += +1 ;

    /*reaction 17: H + OH + M <=> H2O + M */
    nuki[ 1 * kd + 16 ] += -1 ;
    nuki[ 4 * kd + 16 ] += -1 ;
    nuki[ 5 * kd + 16 ] += +1 ;

    /*reaction 18: HCO + M <=> H + CO + M */
    nuki[ 14 * kd + 17 ] += -1 ;
    nuki[ 1 * kd + 17 ] += +1 ;
    nuki[ 12 * kd + 17 ] += +1 ;

    /*reaction 19: O + H2 <=> H + OH */
    nuki[ 2 * kd + 18 ] += -1 ;
    nuki[ 0 * kd + 18 ] += -1 ;
    nuki[ 1 * kd + 18 ] += +1 ;
    nuki[ 4 * kd + 18 ] += +1 ;

    /*reaction 20: O + HO2 <=> OH + O2 */
    nuki[ 2 * kd + 19 ] += -1 ;
    nuki[ 6 * kd + 19 ] += -1 ;
    nuki[ 4 * kd + 19 ] += +1 ;
    nuki[ 3 * kd + 19 ] += +1 ;

    /*reaction 21: O + CH2 <=> H + HCO */
    nuki[ 2 * kd + 20 ] += -1 ;
    nuki[ 8 * kd + 20 ] += -1 ;
    nuki[ 1 * kd + 20 ] += +1 ;
    nuki[ 14 * kd + 20 ] += +1 ;

    /*reaction 22: O + CH2(S) <=> H + HCO */
    nuki[ 2 * kd + 21 ] += -1 ;
    nuki[ 9 * kd + 21 ] += -1 ;
    nuki[ 1 * kd + 21 ] += +1 ;
    nuki[ 14 * kd + 21 ] += +1 ;

    /*reaction 23: O + CH3 <=> H + CH2O */
    nuki[ 2 * kd + 22 ] += -1 ;
    nuki[ 10 * kd + 22 ] += -1 ;
    nuki[ 1 * kd + 22 ] += +1 ;
    nuki[ 15 * kd + 22 ] += +1 ;

    /*reaction 24: O + CH4 <=> OH + CH3 */
    nuki[ 2 * kd + 23 ] += -1 ;
    nuki[ 11 * kd + 23 ] += -1 ;
    nuki[ 4 * kd + 23 ] += +1 ;
    nuki[ 10 * kd + 23 ] += +1 ;

    /*reaction 25: O + HCO <=> OH + CO */
    nuki[ 2 * kd + 24 ] += -1 ;
    nuki[ 14 * kd + 24 ] += -1 ;
    nuki[ 4 * kd + 24 ] += +1 ;
    nuki[ 12 * kd + 24 ] += +1 ;

    /*reaction 26: O + HCO <=> H + CO2 */
    nuki[ 2 * kd + 25 ] += -1 ;
    nuki[ 14 * kd + 25 ] += -1 ;
    nuki[ 1 * kd + 25 ] += +1 ;
    nuki[ 13 * kd + 25 ] += +1 ;

    /*reaction 27: O + CH2O <=> OH + HCO */
    nuki[ 2 * kd + 26 ] += -1 ;
    nuki[ 15 * kd + 26 ] += -1 ;
    nuki[ 4 * kd + 26 ] += +1 ;
    nuki[ 14 * kd + 26 ] += +1 ;

    /*reaction 28: O + C2H2 <=> CH2(S) + CO */
    nuki[ 2 * kd + 27 ] += -1 ;
    nuki[ 17 * kd + 27 ] += -1 ;
    nuki[ 9 * kd + 27 ] += +1 ;
    nuki[ 12 * kd + 27 ] += +1 ;

    /*reaction 29: O + C2H2 <=> CO + CH2 */
    nuki[ 2 * kd + 28 ] += -1 ;
    nuki[ 17 * kd + 28 ] += -1 ;
    nuki[ 12 * kd + 28 ] += +1 ;
    nuki[ 8 * kd + 28 ] += +1 ;

    /*reaction 30: O + C2H4 <=> CH3 + HCO */
    nuki[ 2 * kd + 29 ] += -1 ;
    nuki[ 19 * kd + 29 ] += -1 ;
    nuki[ 10 * kd + 29 ] += +1 ;
    nuki[ 14 * kd + 29 ] += +1 ;

    /*reaction 31: O + C2H5 <=> CH3 + CH2O */
    nuki[ 2 * kd + 30 ] += -1 ;
    nuki[ 20 * kd + 30 ] += -1 ;
    nuki[ 10 * kd + 30 ] += +1 ;
    nuki[ 15 * kd + 30 ] += +1 ;

    /*reaction 32: O + C2H6 <=> OH + C2H5 */
    nuki[ 2 * kd + 31 ] += -1 ;
    nuki[ 21 * kd + 31 ] += -1 ;
    nuki[ 4 * kd + 31 ] += +1 ;
    nuki[ 20 * kd + 31 ] += +1 ;

    /*reaction 33: O2 + CO <=> O + CO2 */
    nuki[ 3 * kd + 32 ] += -1 ;
    nuki[ 12 * kd + 32 ] += -1 ;
    nuki[ 2 * kd + 32 ] += +1 ;
    nuki[ 13 * kd + 32 ] += +1 ;

    /*reaction 34: O2 + CH2O <=> HO2 + HCO */
    nuki[ 3 * kd + 33 ] += -1 ;
    nuki[ 15 * kd + 33 ] += -1 ;
    nuki[ 6 * kd + 33 ] += +1 ;
    nuki[ 14 * kd + 33 ] += +1 ;

    /*reaction 35: H + 2 O2 <=> HO2 + O2 */
    nuki[ 1 * kd + 34 ] += -1 ;
    nuki[ 3 * kd + 34 ] += -2 ;
    nuki[ 6 * kd + 34 ] += +1 ;
    nuki[ 3 * kd + 34 ] += +1 ;

    /*reaction 36: H + O2 + H2O <=> HO2 + H2O */
    nuki[ 1 * kd + 35 ] += -1 ;
    nuki[ 3 * kd + 35 ] += -1 ;
    nuki[ 5 * kd + 35 ] += -1 ;
    nuki[ 6 * kd + 35 ] += +1 ;
    nuki[ 5 * kd + 35 ] += +1 ;

    /*reaction 37: H + O2 + N2 <=> HO2 + N2 */
    nuki[ 1 * kd + 36 ] += -1 ;
    nuki[ 3 * kd + 36 ] += -1 ;
    nuki[ 22 * kd + 36 ] += -1 ;
    nuki[ 6 * kd + 36 ] += +1 ;
    nuki[ 22 * kd + 36 ] += +1 ;

    /*reaction 38: H + O2 + AR <=> HO2 + AR */
    nuki[ 1 * kd + 37 ] += -1 ;
    nuki[ 3 * kd + 37 ] += -1 ;
    nuki[ 23 * kd + 37 ] += -1 ;
    nuki[ 6 * kd + 37 ] += +1 ;
    nuki[ 23 * kd + 37 ] += +1 ;

    /*reaction 39: H + O2 <=> O + OH */
    nuki[ 1 * kd + 38 ] += -1 ;
    nuki[ 3 * kd + 38 ] += -1 ;
    nuki[ 2 * kd + 38 ] += +1 ;
    nuki[ 4 * kd + 38 ] += +1 ;

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
    nuki[ 13 * kd + 41 ] += -1 ;
    nuki[ 0 * kd + 41 ] += +1 ;
    nuki[ 13 * kd + 41 ] += +1 ;

    /*reaction 43: H + HO2 <=> O2 + H2 */
    nuki[ 1 * kd + 42 ] += -1 ;
    nuki[ 6 * kd + 42 ] += -1 ;
    nuki[ 3 * kd + 42 ] += +1 ;
    nuki[ 0 * kd + 42 ] += +1 ;

    /*reaction 44: H + HO2 <=> 2 OH */
    nuki[ 1 * kd + 43 ] += -1 ;
    nuki[ 6 * kd + 43 ] += -1 ;
    nuki[ 4 * kd + 43 ] += +2 ;

    /*reaction 45: H + H2O2 <=> HO2 + H2 */
    nuki[ 1 * kd + 44 ] += -1 ;
    nuki[ 7 * kd + 44 ] += -1 ;
    nuki[ 6 * kd + 44 ] += +1 ;
    nuki[ 0 * kd + 44 ] += +1 ;

    /*reaction 46: H + CH4 <=> CH3 + H2 */
    nuki[ 1 * kd + 45 ] += -1 ;
    nuki[ 11 * kd + 45 ] += -1 ;
    nuki[ 10 * kd + 45 ] += +1 ;
    nuki[ 0 * kd + 45 ] += +1 ;

    /*reaction 47: H + HCO <=> H2 + CO */
    nuki[ 1 * kd + 46 ] += -1 ;
    nuki[ 14 * kd + 46 ] += -1 ;
    nuki[ 0 * kd + 46 ] += +1 ;
    nuki[ 12 * kd + 46 ] += +1 ;

    /*reaction 48: H + CH2O <=> HCO + H2 */
    nuki[ 1 * kd + 47 ] += -1 ;
    nuki[ 15 * kd + 47 ] += -1 ;
    nuki[ 14 * kd + 47 ] += +1 ;
    nuki[ 0 * kd + 47 ] += +1 ;

    /*reaction 49: H + CH3O <=> OH + CH3 */
    nuki[ 1 * kd + 48 ] += -1 ;
    nuki[ 16 * kd + 48 ] += -1 ;
    nuki[ 4 * kd + 48 ] += +1 ;
    nuki[ 10 * kd + 48 ] += +1 ;

    /*reaction 50: H + C2H3 <=> H2 + C2H2 */
    nuki[ 1 * kd + 49 ] += -1 ;
    nuki[ 18 * kd + 49 ] += -1 ;
    nuki[ 0 * kd + 49 ] += +1 ;
    nuki[ 17 * kd + 49 ] += +1 ;

    /*reaction 51: H + C2H4 <=> C2H3 + H2 */
    nuki[ 1 * kd + 50 ] += -1 ;
    nuki[ 19 * kd + 50 ] += -1 ;
    nuki[ 18 * kd + 50 ] += +1 ;
    nuki[ 0 * kd + 50 ] += +1 ;

    /*reaction 52: H + C2H6 <=> C2H5 + H2 */
    nuki[ 1 * kd + 51 ] += -1 ;
    nuki[ 21 * kd + 51 ] += -1 ;
    nuki[ 20 * kd + 51 ] += +1 ;
    nuki[ 0 * kd + 51 ] += +1 ;

    /*reaction 53: OH + H2 <=> H + H2O */
    nuki[ 4 * kd + 52 ] += -1 ;
    nuki[ 0 * kd + 52 ] += -1 ;
    nuki[ 1 * kd + 52 ] += +1 ;
    nuki[ 5 * kd + 52 ] += +1 ;

    /*reaction 54: 2 OH <=> O + H2O */
    nuki[ 4 * kd + 53 ] += -2 ;
    nuki[ 2 * kd + 53 ] += +1 ;
    nuki[ 5 * kd + 53 ] += +1 ;

    /*reaction 55: OH + HO2 <=> O2 + H2O */
    nuki[ 4 * kd + 54 ] += -1 ;
    nuki[ 6 * kd + 54 ] += -1 ;
    nuki[ 3 * kd + 54 ] += +1 ;
    nuki[ 5 * kd + 54 ] += +1 ;

    /*reaction 56: OH + H2O2 <=> HO2 + H2O */
    nuki[ 4 * kd + 55 ] += -1 ;
    nuki[ 7 * kd + 55 ] += -1 ;
    nuki[ 6 * kd + 55 ] += +1 ;
    nuki[ 5 * kd + 55 ] += +1 ;

    /*reaction 57: OH + CH2 <=> H + CH2O */
    nuki[ 4 * kd + 56 ] += -1 ;
    nuki[ 8 * kd + 56 ] += -1 ;
    nuki[ 1 * kd + 56 ] += +1 ;
    nuki[ 15 * kd + 56 ] += +1 ;

    /*reaction 58: OH + CH2(S) <=> H + CH2O */
    nuki[ 4 * kd + 57 ] += -1 ;
    nuki[ 9 * kd + 57 ] += -1 ;
    nuki[ 1 * kd + 57 ] += +1 ;
    nuki[ 15 * kd + 57 ] += +1 ;

    /*reaction 59: OH + CH3 <=> CH2 + H2O */
    nuki[ 4 * kd + 58 ] += -1 ;
    nuki[ 10 * kd + 58 ] += -1 ;
    nuki[ 8 * kd + 58 ] += +1 ;
    nuki[ 5 * kd + 58 ] += +1 ;

    /*reaction 60: OH + CH3 <=> CH2(S) + H2O */
    nuki[ 4 * kd + 59 ] += -1 ;
    nuki[ 10 * kd + 59 ] += -1 ;
    nuki[ 9 * kd + 59 ] += +1 ;
    nuki[ 5 * kd + 59 ] += +1 ;

    /*reaction 61: OH + CH4 <=> CH3 + H2O */
    nuki[ 4 * kd + 60 ] += -1 ;
    nuki[ 11 * kd + 60 ] += -1 ;
    nuki[ 10 * kd + 60 ] += +1 ;
    nuki[ 5 * kd + 60 ] += +1 ;

    /*reaction 62: OH + CO <=> H + CO2 */
    nuki[ 4 * kd + 61 ] += -1 ;
    nuki[ 12 * kd + 61 ] += -1 ;
    nuki[ 1 * kd + 61 ] += +1 ;
    nuki[ 13 * kd + 61 ] += +1 ;

    /*reaction 63: OH + HCO <=> H2O + CO */
    nuki[ 4 * kd + 62 ] += -1 ;
    nuki[ 14 * kd + 62 ] += -1 ;
    nuki[ 5 * kd + 62 ] += +1 ;
    nuki[ 12 * kd + 62 ] += +1 ;

    /*reaction 64: OH + CH2O <=> HCO + H2O */
    nuki[ 4 * kd + 63 ] += -1 ;
    nuki[ 15 * kd + 63 ] += -1 ;
    nuki[ 14 * kd + 63 ] += +1 ;
    nuki[ 5 * kd + 63 ] += +1 ;

    /*reaction 65: OH + C2H2 <=> CH3 + CO */
    nuki[ 4 * kd + 64 ] += -1 ;
    nuki[ 17 * kd + 64 ] += -1 ;
    nuki[ 10 * kd + 64 ] += +1 ;
    nuki[ 12 * kd + 64 ] += +1 ;

    /*reaction 66: OH + C2H3 <=> H2O + C2H2 */
    nuki[ 4 * kd + 65 ] += -1 ;
    nuki[ 18 * kd + 65 ] += -1 ;
    nuki[ 5 * kd + 65 ] += +1 ;
    nuki[ 17 * kd + 65 ] += +1 ;

    /*reaction 67: OH + C2H4 <=> C2H3 + H2O */
    nuki[ 4 * kd + 66 ] += -1 ;
    nuki[ 19 * kd + 66 ] += -1 ;
    nuki[ 18 * kd + 66 ] += +1 ;
    nuki[ 5 * kd + 66 ] += +1 ;

    /*reaction 68: OH + C2H6 <=> C2H5 + H2O */
    nuki[ 4 * kd + 67 ] += -1 ;
    nuki[ 21 * kd + 67 ] += -1 ;
    nuki[ 20 * kd + 67 ] += +1 ;
    nuki[ 5 * kd + 67 ] += +1 ;

    /*reaction 69: 2 HO2 <=> O2 + H2O2 */
    nuki[ 6 * kd + 68 ] += -2 ;
    nuki[ 3 * kd + 68 ] += +1 ;
    nuki[ 7 * kd + 68 ] += +1 ;

    /*reaction 70: 2 HO2 <=> O2 + H2O2 */
    nuki[ 6 * kd + 69 ] += -2 ;
    nuki[ 3 * kd + 69 ] += +1 ;
    nuki[ 7 * kd + 69 ] += +1 ;

    /*reaction 71: HO2 + CH2 <=> OH + CH2O */
    nuki[ 6 * kd + 70 ] += -1 ;
    nuki[ 8 * kd + 70 ] += -1 ;
    nuki[ 4 * kd + 70 ] += +1 ;
    nuki[ 15 * kd + 70 ] += +1 ;

    /*reaction 72: HO2 + CH3 <=> O2 + CH4 */
    nuki[ 6 * kd + 71 ] += -1 ;
    nuki[ 10 * kd + 71 ] += -1 ;
    nuki[ 3 * kd + 71 ] += +1 ;
    nuki[ 11 * kd + 71 ] += +1 ;

    /*reaction 73: HO2 + CH3 <=> OH + CH3O */
    nuki[ 6 * kd + 72 ] += -1 ;
    nuki[ 10 * kd + 72 ] += -1 ;
    nuki[ 4 * kd + 72 ] += +1 ;
    nuki[ 16 * kd + 72 ] += +1 ;

    /*reaction 74: HO2 + CO <=> OH + CO2 */
    nuki[ 6 * kd + 73 ] += -1 ;
    nuki[ 12 * kd + 73 ] += -1 ;
    nuki[ 4 * kd + 73 ] += +1 ;
    nuki[ 13 * kd + 73 ] += +1 ;

    /*reaction 75: HO2 + CH2O <=> HCO + H2O2 */
    nuki[ 6 * kd + 74 ] += -1 ;
    nuki[ 15 * kd + 74 ] += -1 ;
    nuki[ 14 * kd + 74 ] += +1 ;
    nuki[ 7 * kd + 74 ] += +1 ;

    /*reaction 76: CH2 + O2 <=> OH + HCO */
    nuki[ 8 * kd + 75 ] += -1 ;
    nuki[ 3 * kd + 75 ] += -1 ;
    nuki[ 4 * kd + 75 ] += +1 ;
    nuki[ 14 * kd + 75 ] += +1 ;

    /*reaction 77: CH2 + H2 <=> H + CH3 */
    nuki[ 8 * kd + 76 ] += -1 ;
    nuki[ 0 * kd + 76 ] += -1 ;
    nuki[ 1 * kd + 76 ] += +1 ;
    nuki[ 10 * kd + 76 ] += +1 ;

    /*reaction 78: 2 CH2 <=> H2 + C2H2 */
    nuki[ 8 * kd + 77 ] += -2 ;
    nuki[ 0 * kd + 77 ] += +1 ;
    nuki[ 17 * kd + 77 ] += +1 ;

    /*reaction 79: CH2 + CH3 <=> H + C2H4 */
    nuki[ 8 * kd + 78 ] += -1 ;
    nuki[ 10 * kd + 78 ] += -1 ;
    nuki[ 1 * kd + 78 ] += +1 ;
    nuki[ 19 * kd + 78 ] += +1 ;

    /*reaction 80: CH2 + CH4 <=> 2 CH3 */
    nuki[ 8 * kd + 79 ] += -1 ;
    nuki[ 11 * kd + 79 ] += -1 ;
    nuki[ 10 * kd + 79 ] += +2 ;

    /*reaction 81: CH2(S) + N2 <=> CH2 + N2 */
    nuki[ 9 * kd + 80 ] += -1 ;
    nuki[ 22 * kd + 80 ] += -1 ;
    nuki[ 8 * kd + 80 ] += +1 ;
    nuki[ 22 * kd + 80 ] += +1 ;

    /*reaction 82: CH2(S) + AR <=> CH2 + AR */
    nuki[ 9 * kd + 81 ] += -1 ;
    nuki[ 23 * kd + 81 ] += -1 ;
    nuki[ 8 * kd + 81 ] += +1 ;
    nuki[ 23 * kd + 81 ] += +1 ;

    /*reaction 83: CH2(S) + O2 <=> H + OH + CO */
    nuki[ 9 * kd + 82 ] += -1 ;
    nuki[ 3 * kd + 82 ] += -1 ;
    nuki[ 1 * kd + 82 ] += +1 ;
    nuki[ 4 * kd + 82 ] += +1 ;
    nuki[ 12 * kd + 82 ] += +1 ;

    /*reaction 84: CH2(S) + O2 <=> CO + H2O */
    nuki[ 9 * kd + 83 ] += -1 ;
    nuki[ 3 * kd + 83 ] += -1 ;
    nuki[ 12 * kd + 83 ] += +1 ;
    nuki[ 5 * kd + 83 ] += +1 ;

    /*reaction 85: CH2(S) + H2 <=> CH3 + H */
    nuki[ 9 * kd + 84 ] += -1 ;
    nuki[ 0 * kd + 84 ] += -1 ;
    nuki[ 10 * kd + 84 ] += +1 ;
    nuki[ 1 * kd + 84 ] += +1 ;

    /*reaction 86: CH2(S) + H2O <=> CH2 + H2O */
    nuki[ 9 * kd + 85 ] += -1 ;
    nuki[ 5 * kd + 85 ] += -1 ;
    nuki[ 8 * kd + 85 ] += +1 ;
    nuki[ 5 * kd + 85 ] += +1 ;

    /*reaction 87: CH2(S) + CH3 <=> H + C2H4 */
    nuki[ 9 * kd + 86 ] += -1 ;
    nuki[ 10 * kd + 86 ] += -1 ;
    nuki[ 1 * kd + 86 ] += +1 ;
    nuki[ 19 * kd + 86 ] += +1 ;

    /*reaction 88: CH2(S) + CH4 <=> 2 CH3 */
    nuki[ 9 * kd + 87 ] += -1 ;
    nuki[ 11 * kd + 87 ] += -1 ;
    nuki[ 10 * kd + 87 ] += +2 ;

    /*reaction 89: CH2(S) + CO <=> CH2 + CO */
    nuki[ 9 * kd + 88 ] += -1 ;
    nuki[ 12 * kd + 88 ] += -1 ;
    nuki[ 8 * kd + 88 ] += +1 ;
    nuki[ 12 * kd + 88 ] += +1 ;

    /*reaction 90: CH2(S) + CO2 <=> CH2 + CO2 */
    nuki[ 9 * kd + 89 ] += -1 ;
    nuki[ 13 * kd + 89 ] += -1 ;
    nuki[ 8 * kd + 89 ] += +1 ;
    nuki[ 13 * kd + 89 ] += +1 ;

    /*reaction 91: CH2(S) + CO2 <=> CO + CH2O */
    nuki[ 9 * kd + 90 ] += -1 ;
    nuki[ 13 * kd + 90 ] += -1 ;
    nuki[ 12 * kd + 90 ] += +1 ;
    nuki[ 15 * kd + 90 ] += +1 ;

    /*reaction 92: CH3 + O2 <=> O + CH3O */
    nuki[ 10 * kd + 91 ] += -1 ;
    nuki[ 3 * kd + 91 ] += -1 ;
    nuki[ 2 * kd + 91 ] += +1 ;
    nuki[ 16 * kd + 91 ] += +1 ;

    /*reaction 93: CH3 + O2 <=> OH + CH2O */
    nuki[ 10 * kd + 92 ] += -1 ;
    nuki[ 3 * kd + 92 ] += -1 ;
    nuki[ 4 * kd + 92 ] += +1 ;
    nuki[ 15 * kd + 92 ] += +1 ;

    /*reaction 94: CH3 + H2O2 <=> HO2 + CH4 */
    nuki[ 10 * kd + 93 ] += -1 ;
    nuki[ 7 * kd + 93 ] += -1 ;
    nuki[ 6 * kd + 93 ] += +1 ;
    nuki[ 11 * kd + 93 ] += +1 ;

    /*reaction 95: 2 CH3 <=> H + C2H5 */
    nuki[ 10 * kd + 94 ] += -2 ;
    nuki[ 1 * kd + 94 ] += +1 ;
    nuki[ 20 * kd + 94 ] += +1 ;

    /*reaction 96: CH3 + HCO <=> CH4 + CO */
    nuki[ 10 * kd + 95 ] += -1 ;
    nuki[ 14 * kd + 95 ] += -1 ;
    nuki[ 11 * kd + 95 ] += +1 ;
    nuki[ 12 * kd + 95 ] += +1 ;

    /*reaction 97: CH3 + CH2O <=> HCO + CH4 */
    nuki[ 10 * kd + 96 ] += -1 ;
    nuki[ 15 * kd + 96 ] += -1 ;
    nuki[ 14 * kd + 96 ] += +1 ;
    nuki[ 11 * kd + 96 ] += +1 ;

    /*reaction 98: CH3 + C2H4 <=> C2H3 + CH4 */
    nuki[ 10 * kd + 97 ] += -1 ;
    nuki[ 19 * kd + 97 ] += -1 ;
    nuki[ 18 * kd + 97 ] += +1 ;
    nuki[ 11 * kd + 97 ] += +1 ;

    /*reaction 99: CH3 + C2H6 <=> C2H5 + CH4 */
    nuki[ 10 * kd + 98 ] += -1 ;
    nuki[ 21 * kd + 98 ] += -1 ;
    nuki[ 20 * kd + 98 ] += +1 ;
    nuki[ 11 * kd + 98 ] += +1 ;

    /*reaction 100: HCO + H2O <=> H + CO + H2O */
    nuki[ 14 * kd + 99 ] += -1 ;
    nuki[ 5 * kd + 99 ] += -1 ;
    nuki[ 1 * kd + 99 ] += +1 ;
    nuki[ 12 * kd + 99 ] += +1 ;
    nuki[ 5 * kd + 99 ] += +1 ;

    /*reaction 101: HCO + O2 <=> HO2 + CO */
    nuki[ 14 * kd + 100 ] += -1 ;
    nuki[ 3 * kd + 100 ] += -1 ;
    nuki[ 6 * kd + 100 ] += +1 ;
    nuki[ 12 * kd + 100 ] += +1 ;

    /*reaction 102: CH3O + O2 <=> HO2 + CH2O */
    nuki[ 16 * kd + 101 ] += -1 ;
    nuki[ 3 * kd + 101 ] += -1 ;
    nuki[ 6 * kd + 101 ] += +1 ;
    nuki[ 15 * kd + 101 ] += +1 ;

    /*reaction 103: C2H3 + O2 <=> HCO + CH2O */
    nuki[ 18 * kd + 102 ] += -1 ;
    nuki[ 3 * kd + 102 ] += -1 ;
    nuki[ 14 * kd + 102 ] += +1 ;
    nuki[ 15 * kd + 102 ] += +1 ;

    /*reaction 104: C2H5 + O2 <=> HO2 + C2H4 */
    nuki[ 20 * kd + 103 ] += -1 ;
    nuki[ 3 * kd + 103 ] += -1 ;
    nuki[ 6 * kd + 103 ] += +1 ;
    nuki[ 19 * kd + 103 ] += +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * restrict rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < kd * 24; ++ id) {
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

    /*CH2 */
    ncf[ 8 * kd + 2 ] = 1; /*C */
    ncf[ 8 * kd + 1 ] = 2; /*H */

    /*CH2(S) */
    ncf[ 9 * kd + 2 ] = 1; /*C */
    ncf[ 9 * kd + 1 ] = 2; /*H */

    /*CH3 */
    ncf[ 10 * kd + 2 ] = 1; /*C */
    ncf[ 10 * kd + 1 ] = 3; /*H */

    /*CH4 */
    ncf[ 11 * kd + 2 ] = 1; /*C */
    ncf[ 11 * kd + 1 ] = 4; /*H */

    /*CO */
    ncf[ 12 * kd + 2 ] = 1; /*C */
    ncf[ 12 * kd + 0 ] = 1; /*O */

    /*CO2 */
    ncf[ 13 * kd + 2 ] = 1; /*C */
    ncf[ 13 * kd + 0 ] = 2; /*O */

    /*HCO */
    ncf[ 14 * kd + 1 ] = 1; /*H */
    ncf[ 14 * kd + 2 ] = 1; /*C */
    ncf[ 14 * kd + 0 ] = 1; /*O */

    /*CH2O */
    ncf[ 15 * kd + 1 ] = 2; /*H */
    ncf[ 15 * kd + 2 ] = 1; /*C */
    ncf[ 15 * kd + 0 ] = 1; /*O */

    /*CH3O */
    ncf[ 16 * kd + 2 ] = 1; /*C */
    ncf[ 16 * kd + 1 ] = 3; /*H */
    ncf[ 16 * kd + 0 ] = 1; /*O */

    /*C2H2 */
    ncf[ 17 * kd + 2 ] = 2; /*C */
    ncf[ 17 * kd + 1 ] = 2; /*H */

    /*C2H3 */
    ncf[ 18 * kd + 2 ] = 2; /*C */
    ncf[ 18 * kd + 1 ] = 3; /*H */

    /*C2H4 */
    ncf[ 19 * kd + 2 ] = 2; /*C */
    ncf[ 19 * kd + 1 ] = 4; /*H */

    /*C2H5 */
    ncf[ 20 * kd + 2 ] = 2; /*C */
    ncf[ 20 * kd + 1 ] = 5; /*H */

    /*C2H6 */
    ncf[ 21 * kd + 2 ] = 2; /*C */
    ncf[ 21 * kd + 1 ] = 6; /*H */

    /*N2 */
    ncf[ 22 * kd + 3 ] = 2; /*N */

    /*AR */
    ncf[ 23 * kd + 4 ] = 1; /*AR */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE(int * iwrk, double * restrict rwrk, double * restrict a, double * restrict b, double * restrict e)
{
    for (int i=0; i<104; ++i) {
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
    double gort[24]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + C2H2 (+M) <=> C2H3 (+M) */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + C2H3 (+M) <=> C2H4 (+M) */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[8] *= 1e+06; 

    /*reaction 10: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[9] *= 1e+06; 

    /*reaction 11: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[10] *= 1e+06; 

    /*reaction 12: C2H4 (+M) <=> H2 + C2H2 (+M) */
    eqcon[11] *= 1e-06; 

    /*reaction 13: O + H + M <=> OH + M */
    eqcon[12] *= 1e+06; 

    /*reaction 14: O + CO + M <=> CO2 + M */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + O2 + M <=> HO2 + M */
    eqcon[14] *= 1e+06; 

    /*reaction 16: 2 H + M <=> H2 + M */
    eqcon[15] *= 1e+06; 

    /*reaction 17: H + OH + M <=> H2O + M */
    eqcon[16] *= 1e+06; 

    /*reaction 18: HCO + M <=> H + CO + M */
    eqcon[17] *= 1e-06; 

    /*reaction 19: O + H2 <=> H + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: O + HO2 <=> OH + O2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: O + CH2 <=> H + HCO */
    /*eqcon[20] *= 1;  */

    /*reaction 22: O + CH2(S) <=> H + HCO */
    /*eqcon[21] *= 1;  */

    /*reaction 23: O + CH3 <=> H + CH2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: O + CH4 <=> OH + CH3 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: O + HCO <=> OH + CO */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O + HCO <=> H + CO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: O + CH2O <=> OH + HCO */
    /*eqcon[26] *= 1;  */

    /*reaction 28: O + C2H2 <=> CH2(S) + CO */
    /*eqcon[27] *= 1;  */

    /*reaction 29: O + C2H2 <=> CO + CH2 */
    /*eqcon[28] *= 1;  */

    /*reaction 30: O + C2H4 <=> CH3 + HCO */
    /*eqcon[29] *= 1;  */

    /*reaction 31: O + C2H5 <=> CH3 + CH2O */
    /*eqcon[30] *= 1;  */

    /*reaction 32: O + C2H6 <=> OH + C2H5 */
    /*eqcon[31] *= 1;  */

    /*reaction 33: O2 + CO <=> O + CO2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: O2 + CH2O <=> HO2 + HCO */
    /*eqcon[33] *= 1;  */

    /*reaction 35: H + 2 O2 <=> HO2 + O2 */
    eqcon[34] *= 1e+06; 

    /*reaction 36: H + O2 + H2O <=> HO2 + H2O */
    eqcon[35] *= 1e+06; 

    /*reaction 37: H + O2 + N2 <=> HO2 + N2 */
    eqcon[36] *= 1e+06; 

    /*reaction 38: H + O2 + AR <=> HO2 + AR */
    eqcon[37] *= 1e+06; 

    /*reaction 39: H + O2 <=> O + OH */
    /*eqcon[38] *= 1;  */

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    eqcon[39] *= 1e+06; 

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    eqcon[40] *= 1e+06; 

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    eqcon[41] *= 1e+06; 

    /*reaction 43: H + HO2 <=> O2 + H2 */
    /*eqcon[42] *= 1;  */

    /*reaction 44: H + HO2 <=> 2 OH */
    /*eqcon[43] *= 1;  */

    /*reaction 45: H + H2O2 <=> HO2 + H2 */
    /*eqcon[44] *= 1;  */

    /*reaction 46: H + CH4 <=> CH3 + H2 */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + HCO <=> H2 + CO */
    /*eqcon[46] *= 1;  */

    /*reaction 48: H + CH2O <=> HCO + H2 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: H + CH3O <=> OH + CH3 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: H + C2H3 <=> H2 + C2H2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: H + C2H4 <=> C2H3 + H2 */
    /*eqcon[50] *= 1;  */

    /*reaction 52: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: OH + H2 <=> H + H2O */
    /*eqcon[52] *= 1;  */

    /*reaction 54: 2 OH <=> O + H2O */
    /*eqcon[53] *= 1;  */

    /*reaction 55: OH + HO2 <=> O2 + H2O */
    /*eqcon[54] *= 1;  */

    /*reaction 56: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: OH + CH2 <=> H + CH2O */
    /*eqcon[56] *= 1;  */

    /*reaction 58: OH + CH2(S) <=> H + CH2O */
    /*eqcon[57] *= 1;  */

    /*reaction 59: OH + CH3 <=> CH2 + H2O */
    /*eqcon[58] *= 1;  */

    /*reaction 60: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[59] *= 1;  */

    /*reaction 61: OH + CH4 <=> CH3 + H2O */
    /*eqcon[60] *= 1;  */

    /*reaction 62: OH + CO <=> H + CO2 */
    /*eqcon[61] *= 1;  */

    /*reaction 63: OH + HCO <=> H2O + CO */
    /*eqcon[62] *= 1;  */

    /*reaction 64: OH + CH2O <=> HCO + H2O */
    /*eqcon[63] *= 1;  */

    /*reaction 65: OH + C2H2 <=> CH3 + CO */
    /*eqcon[64] *= 1;  */

    /*reaction 66: OH + C2H3 <=> H2O + C2H2 */
    /*eqcon[65] *= 1;  */

    /*reaction 67: OH + C2H4 <=> C2H3 + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[67] *= 1;  */

    /*reaction 69: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[68] *= 1;  */

    /*reaction 70: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[69] *= 1;  */

    /*reaction 71: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[70] *= 1;  */

    /*reaction 72: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[71] *= 1;  */

    /*reaction 73: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[72] *= 1;  */

    /*reaction 74: HO2 + CO <=> OH + CO2 */
    /*eqcon[73] *= 1;  */

    /*reaction 75: HO2 + CH2O <=> HCO + H2O2 */
    /*eqcon[74] *= 1;  */

    /*reaction 76: CH2 + O2 <=> OH + HCO */
    /*eqcon[75] *= 1;  */

    /*reaction 77: CH2 + H2 <=> H + CH3 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: 2 CH2 <=> H2 + C2H2 */
    /*eqcon[77] *= 1;  */

    /*reaction 79: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[79] *= 1;  */

    /*reaction 81: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[80] *= 1;  */

    /*reaction 82: CH2(S) + AR <=> CH2 + AR */
    /*eqcon[81] *= 1;  */

    /*reaction 83: CH2(S) + O2 <=> H + OH + CO */
    eqcon[82] *= 1e-06; 

    /*reaction 84: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[83] *= 1;  */

    /*reaction 85: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[84] *= 1;  */

    /*reaction 86: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[85] *= 1;  */

    /*reaction 87: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[86] *= 1;  */

    /*reaction 88: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[87] *= 1;  */

    /*reaction 89: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[88] *= 1;  */

    /*reaction 90: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[89] *= 1;  */

    /*reaction 91: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[90] *= 1;  */

    /*reaction 92: CH3 + O2 <=> O + CH3O */
    /*eqcon[91] *= 1;  */

    /*reaction 93: CH3 + O2 <=> OH + CH2O */
    /*eqcon[92] *= 1;  */

    /*reaction 94: CH3 + H2O2 <=> HO2 + CH4 */
    /*eqcon[93] *= 1;  */

    /*reaction 95: 2 CH3 <=> H + C2H5 */
    /*eqcon[94] *= 1;  */

    /*reaction 96: CH3 + HCO <=> CH4 + CO */
    /*eqcon[95] *= 1;  */

    /*reaction 97: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[96] *= 1;  */

    /*reaction 98: CH3 + C2H4 <=> C2H3 + CH4 */
    /*eqcon[97] *= 1;  */

    /*reaction 99: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[98] *= 1;  */

    /*reaction 100: HCO + H2O <=> H + CO + H2O */
    eqcon[99] *= 1e-06; 

    /*reaction 101: HCO + O2 <=> HO2 + CO */
    /*eqcon[100] *= 1;  */

    /*reaction 102: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[101] *= 1;  */

    /*reaction 103: C2H3 + O2 <=> HCO + CH2O */
    /*eqcon[102] *= 1;  */

    /*reaction 104: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[103] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[24]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + C2H2 (+M) <=> C2H3 (+M) */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + C2H3 (+M) <=> C2H4 (+M) */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[8] *= 1e+06; 

    /*reaction 10: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[9] *= 1e+06; 

    /*reaction 11: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[10] *= 1e+06; 

    /*reaction 12: C2H4 (+M) <=> H2 + C2H2 (+M) */
    eqcon[11] *= 1e-06; 

    /*reaction 13: O + H + M <=> OH + M */
    eqcon[12] *= 1e+06; 

    /*reaction 14: O + CO + M <=> CO2 + M */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + O2 + M <=> HO2 + M */
    eqcon[14] *= 1e+06; 

    /*reaction 16: 2 H + M <=> H2 + M */
    eqcon[15] *= 1e+06; 

    /*reaction 17: H + OH + M <=> H2O + M */
    eqcon[16] *= 1e+06; 

    /*reaction 18: HCO + M <=> H + CO + M */
    eqcon[17] *= 1e-06; 

    /*reaction 19: O + H2 <=> H + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: O + HO2 <=> OH + O2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: O + CH2 <=> H + HCO */
    /*eqcon[20] *= 1;  */

    /*reaction 22: O + CH2(S) <=> H + HCO */
    /*eqcon[21] *= 1;  */

    /*reaction 23: O + CH3 <=> H + CH2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: O + CH4 <=> OH + CH3 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: O + HCO <=> OH + CO */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O + HCO <=> H + CO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: O + CH2O <=> OH + HCO */
    /*eqcon[26] *= 1;  */

    /*reaction 28: O + C2H2 <=> CH2(S) + CO */
    /*eqcon[27] *= 1;  */

    /*reaction 29: O + C2H2 <=> CO + CH2 */
    /*eqcon[28] *= 1;  */

    /*reaction 30: O + C2H4 <=> CH3 + HCO */
    /*eqcon[29] *= 1;  */

    /*reaction 31: O + C2H5 <=> CH3 + CH2O */
    /*eqcon[30] *= 1;  */

    /*reaction 32: O + C2H6 <=> OH + C2H5 */
    /*eqcon[31] *= 1;  */

    /*reaction 33: O2 + CO <=> O + CO2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: O2 + CH2O <=> HO2 + HCO */
    /*eqcon[33] *= 1;  */

    /*reaction 35: H + 2 O2 <=> HO2 + O2 */
    eqcon[34] *= 1e+06; 

    /*reaction 36: H + O2 + H2O <=> HO2 + H2O */
    eqcon[35] *= 1e+06; 

    /*reaction 37: H + O2 + N2 <=> HO2 + N2 */
    eqcon[36] *= 1e+06; 

    /*reaction 38: H + O2 + AR <=> HO2 + AR */
    eqcon[37] *= 1e+06; 

    /*reaction 39: H + O2 <=> O + OH */
    /*eqcon[38] *= 1;  */

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    eqcon[39] *= 1e+06; 

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    eqcon[40] *= 1e+06; 

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    eqcon[41] *= 1e+06; 

    /*reaction 43: H + HO2 <=> O2 + H2 */
    /*eqcon[42] *= 1;  */

    /*reaction 44: H + HO2 <=> 2 OH */
    /*eqcon[43] *= 1;  */

    /*reaction 45: H + H2O2 <=> HO2 + H2 */
    /*eqcon[44] *= 1;  */

    /*reaction 46: H + CH4 <=> CH3 + H2 */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + HCO <=> H2 + CO */
    /*eqcon[46] *= 1;  */

    /*reaction 48: H + CH2O <=> HCO + H2 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: H + CH3O <=> OH + CH3 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: H + C2H3 <=> H2 + C2H2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: H + C2H4 <=> C2H3 + H2 */
    /*eqcon[50] *= 1;  */

    /*reaction 52: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: OH + H2 <=> H + H2O */
    /*eqcon[52] *= 1;  */

    /*reaction 54: 2 OH <=> O + H2O */
    /*eqcon[53] *= 1;  */

    /*reaction 55: OH + HO2 <=> O2 + H2O */
    /*eqcon[54] *= 1;  */

    /*reaction 56: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: OH + CH2 <=> H + CH2O */
    /*eqcon[56] *= 1;  */

    /*reaction 58: OH + CH2(S) <=> H + CH2O */
    /*eqcon[57] *= 1;  */

    /*reaction 59: OH + CH3 <=> CH2 + H2O */
    /*eqcon[58] *= 1;  */

    /*reaction 60: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[59] *= 1;  */

    /*reaction 61: OH + CH4 <=> CH3 + H2O */
    /*eqcon[60] *= 1;  */

    /*reaction 62: OH + CO <=> H + CO2 */
    /*eqcon[61] *= 1;  */

    /*reaction 63: OH + HCO <=> H2O + CO */
    /*eqcon[62] *= 1;  */

    /*reaction 64: OH + CH2O <=> HCO + H2O */
    /*eqcon[63] *= 1;  */

    /*reaction 65: OH + C2H2 <=> CH3 + CO */
    /*eqcon[64] *= 1;  */

    /*reaction 66: OH + C2H3 <=> H2O + C2H2 */
    /*eqcon[65] *= 1;  */

    /*reaction 67: OH + C2H4 <=> C2H3 + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[67] *= 1;  */

    /*reaction 69: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[68] *= 1;  */

    /*reaction 70: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[69] *= 1;  */

    /*reaction 71: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[70] *= 1;  */

    /*reaction 72: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[71] *= 1;  */

    /*reaction 73: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[72] *= 1;  */

    /*reaction 74: HO2 + CO <=> OH + CO2 */
    /*eqcon[73] *= 1;  */

    /*reaction 75: HO2 + CH2O <=> HCO + H2O2 */
    /*eqcon[74] *= 1;  */

    /*reaction 76: CH2 + O2 <=> OH + HCO */
    /*eqcon[75] *= 1;  */

    /*reaction 77: CH2 + H2 <=> H + CH3 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: 2 CH2 <=> H2 + C2H2 */
    /*eqcon[77] *= 1;  */

    /*reaction 79: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[79] *= 1;  */

    /*reaction 81: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[80] *= 1;  */

    /*reaction 82: CH2(S) + AR <=> CH2 + AR */
    /*eqcon[81] *= 1;  */

    /*reaction 83: CH2(S) + O2 <=> H + OH + CO */
    eqcon[82] *= 1e-06; 

    /*reaction 84: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[83] *= 1;  */

    /*reaction 85: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[84] *= 1;  */

    /*reaction 86: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[85] *= 1;  */

    /*reaction 87: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[86] *= 1;  */

    /*reaction 88: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[87] *= 1;  */

    /*reaction 89: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[88] *= 1;  */

    /*reaction 90: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[89] *= 1;  */

    /*reaction 91: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[90] *= 1;  */

    /*reaction 92: CH3 + O2 <=> O + CH3O */
    /*eqcon[91] *= 1;  */

    /*reaction 93: CH3 + O2 <=> OH + CH2O */
    /*eqcon[92] *= 1;  */

    /*reaction 94: CH3 + H2O2 <=> HO2 + CH4 */
    /*eqcon[93] *= 1;  */

    /*reaction 95: 2 CH3 <=> H + C2H5 */
    /*eqcon[94] *= 1;  */

    /*reaction 96: CH3 + HCO <=> CH4 + CO */
    /*eqcon[95] *= 1;  */

    /*reaction 97: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[96] *= 1;  */

    /*reaction 98: CH3 + C2H4 <=> C2H3 + CH4 */
    /*eqcon[97] *= 1;  */

    /*reaction 99: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[98] *= 1;  */

    /*reaction 100: HCO + H2O <=> H + CO + H2O */
    eqcon[99] *= 1e-06; 

    /*reaction 101: HCO + O2 <=> HO2 + CO */
    /*eqcon[100] *= 1;  */

    /*reaction 102: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[101] *= 1;  */

    /*reaction 103: C2H3 + O2 <=> HCO + CH2O */
    /*eqcon[102] *= 1;  */

    /*reaction 104: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[103] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[24]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + C2H2 (+M) <=> C2H3 (+M) */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + C2H3 (+M) <=> C2H4 (+M) */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[8] *= 1e+06; 

    /*reaction 10: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[9] *= 1e+06; 

    /*reaction 11: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[10] *= 1e+06; 

    /*reaction 12: C2H4 (+M) <=> H2 + C2H2 (+M) */
    eqcon[11] *= 1e-06; 

    /*reaction 13: O + H + M <=> OH + M */
    eqcon[12] *= 1e+06; 

    /*reaction 14: O + CO + M <=> CO2 + M */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + O2 + M <=> HO2 + M */
    eqcon[14] *= 1e+06; 

    /*reaction 16: 2 H + M <=> H2 + M */
    eqcon[15] *= 1e+06; 

    /*reaction 17: H + OH + M <=> H2O + M */
    eqcon[16] *= 1e+06; 

    /*reaction 18: HCO + M <=> H + CO + M */
    eqcon[17] *= 1e-06; 

    /*reaction 19: O + H2 <=> H + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: O + HO2 <=> OH + O2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: O + CH2 <=> H + HCO */
    /*eqcon[20] *= 1;  */

    /*reaction 22: O + CH2(S) <=> H + HCO */
    /*eqcon[21] *= 1;  */

    /*reaction 23: O + CH3 <=> H + CH2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: O + CH4 <=> OH + CH3 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: O + HCO <=> OH + CO */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O + HCO <=> H + CO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: O + CH2O <=> OH + HCO */
    /*eqcon[26] *= 1;  */

    /*reaction 28: O + C2H2 <=> CH2(S) + CO */
    /*eqcon[27] *= 1;  */

    /*reaction 29: O + C2H2 <=> CO + CH2 */
    /*eqcon[28] *= 1;  */

    /*reaction 30: O + C2H4 <=> CH3 + HCO */
    /*eqcon[29] *= 1;  */

    /*reaction 31: O + C2H5 <=> CH3 + CH2O */
    /*eqcon[30] *= 1;  */

    /*reaction 32: O + C2H6 <=> OH + C2H5 */
    /*eqcon[31] *= 1;  */

    /*reaction 33: O2 + CO <=> O + CO2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: O2 + CH2O <=> HO2 + HCO */
    /*eqcon[33] *= 1;  */

    /*reaction 35: H + 2 O2 <=> HO2 + O2 */
    eqcon[34] *= 1e+06; 

    /*reaction 36: H + O2 + H2O <=> HO2 + H2O */
    eqcon[35] *= 1e+06; 

    /*reaction 37: H + O2 + N2 <=> HO2 + N2 */
    eqcon[36] *= 1e+06; 

    /*reaction 38: H + O2 + AR <=> HO2 + AR */
    eqcon[37] *= 1e+06; 

    /*reaction 39: H + O2 <=> O + OH */
    /*eqcon[38] *= 1;  */

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    eqcon[39] *= 1e+06; 

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    eqcon[40] *= 1e+06; 

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    eqcon[41] *= 1e+06; 

    /*reaction 43: H + HO2 <=> O2 + H2 */
    /*eqcon[42] *= 1;  */

    /*reaction 44: H + HO2 <=> 2 OH */
    /*eqcon[43] *= 1;  */

    /*reaction 45: H + H2O2 <=> HO2 + H2 */
    /*eqcon[44] *= 1;  */

    /*reaction 46: H + CH4 <=> CH3 + H2 */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + HCO <=> H2 + CO */
    /*eqcon[46] *= 1;  */

    /*reaction 48: H + CH2O <=> HCO + H2 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: H + CH3O <=> OH + CH3 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: H + C2H3 <=> H2 + C2H2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: H + C2H4 <=> C2H3 + H2 */
    /*eqcon[50] *= 1;  */

    /*reaction 52: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: OH + H2 <=> H + H2O */
    /*eqcon[52] *= 1;  */

    /*reaction 54: 2 OH <=> O + H2O */
    /*eqcon[53] *= 1;  */

    /*reaction 55: OH + HO2 <=> O2 + H2O */
    /*eqcon[54] *= 1;  */

    /*reaction 56: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: OH + CH2 <=> H + CH2O */
    /*eqcon[56] *= 1;  */

    /*reaction 58: OH + CH2(S) <=> H + CH2O */
    /*eqcon[57] *= 1;  */

    /*reaction 59: OH + CH3 <=> CH2 + H2O */
    /*eqcon[58] *= 1;  */

    /*reaction 60: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[59] *= 1;  */

    /*reaction 61: OH + CH4 <=> CH3 + H2O */
    /*eqcon[60] *= 1;  */

    /*reaction 62: OH + CO <=> H + CO2 */
    /*eqcon[61] *= 1;  */

    /*reaction 63: OH + HCO <=> H2O + CO */
    /*eqcon[62] *= 1;  */

    /*reaction 64: OH + CH2O <=> HCO + H2O */
    /*eqcon[63] *= 1;  */

    /*reaction 65: OH + C2H2 <=> CH3 + CO */
    /*eqcon[64] *= 1;  */

    /*reaction 66: OH + C2H3 <=> H2O + C2H2 */
    /*eqcon[65] *= 1;  */

    /*reaction 67: OH + C2H4 <=> C2H3 + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[67] *= 1;  */

    /*reaction 69: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[68] *= 1;  */

    /*reaction 70: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[69] *= 1;  */

    /*reaction 71: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[70] *= 1;  */

    /*reaction 72: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[71] *= 1;  */

    /*reaction 73: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[72] *= 1;  */

    /*reaction 74: HO2 + CO <=> OH + CO2 */
    /*eqcon[73] *= 1;  */

    /*reaction 75: HO2 + CH2O <=> HCO + H2O2 */
    /*eqcon[74] *= 1;  */

    /*reaction 76: CH2 + O2 <=> OH + HCO */
    /*eqcon[75] *= 1;  */

    /*reaction 77: CH2 + H2 <=> H + CH3 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: 2 CH2 <=> H2 + C2H2 */
    /*eqcon[77] *= 1;  */

    /*reaction 79: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[79] *= 1;  */

    /*reaction 81: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[80] *= 1;  */

    /*reaction 82: CH2(S) + AR <=> CH2 + AR */
    /*eqcon[81] *= 1;  */

    /*reaction 83: CH2(S) + O2 <=> H + OH + CO */
    eqcon[82] *= 1e-06; 

    /*reaction 84: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[83] *= 1;  */

    /*reaction 85: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[84] *= 1;  */

    /*reaction 86: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[85] *= 1;  */

    /*reaction 87: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[86] *= 1;  */

    /*reaction 88: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[87] *= 1;  */

    /*reaction 89: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[88] *= 1;  */

    /*reaction 90: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[89] *= 1;  */

    /*reaction 91: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[90] *= 1;  */

    /*reaction 92: CH3 + O2 <=> O + CH3O */
    /*eqcon[91] *= 1;  */

    /*reaction 93: CH3 + O2 <=> OH + CH2O */
    /*eqcon[92] *= 1;  */

    /*reaction 94: CH3 + H2O2 <=> HO2 + CH4 */
    /*eqcon[93] *= 1;  */

    /*reaction 95: 2 CH3 <=> H + C2H5 */
    /*eqcon[94] *= 1;  */

    /*reaction 96: CH3 + HCO <=> CH4 + CO */
    /*eqcon[95] *= 1;  */

    /*reaction 97: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[96] *= 1;  */

    /*reaction 98: CH3 + C2H4 <=> C2H3 + CH4 */
    /*eqcon[97] *= 1;  */

    /*reaction 99: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[98] *= 1;  */

    /*reaction 100: HCO + H2O <=> H + CO + H2O */
    eqcon[99] *= 1e-06; 

    /*reaction 101: HCO + O2 <=> HO2 + CO */
    /*eqcon[100] *= 1;  */

    /*reaction 102: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[101] *= 1;  */

    /*reaction 103: C2H3 + O2 <=> HCO + CH2O */
    /*eqcon[102] *= 1;  */

    /*reaction 104: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[103] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[24]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + C2H2 (+M) <=> C2H3 (+M) */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + C2H3 (+M) <=> C2H4 (+M) */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[8] *= 1e+06; 

    /*reaction 10: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[9] *= 1e+06; 

    /*reaction 11: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[10] *= 1e+06; 

    /*reaction 12: C2H4 (+M) <=> H2 + C2H2 (+M) */
    eqcon[11] *= 1e-06; 

    /*reaction 13: O + H + M <=> OH + M */
    eqcon[12] *= 1e+06; 

    /*reaction 14: O + CO + M <=> CO2 + M */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + O2 + M <=> HO2 + M */
    eqcon[14] *= 1e+06; 

    /*reaction 16: 2 H + M <=> H2 + M */
    eqcon[15] *= 1e+06; 

    /*reaction 17: H + OH + M <=> H2O + M */
    eqcon[16] *= 1e+06; 

    /*reaction 18: HCO + M <=> H + CO + M */
    eqcon[17] *= 1e-06; 

    /*reaction 19: O + H2 <=> H + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: O + HO2 <=> OH + O2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: O + CH2 <=> H + HCO */
    /*eqcon[20] *= 1;  */

    /*reaction 22: O + CH2(S) <=> H + HCO */
    /*eqcon[21] *= 1;  */

    /*reaction 23: O + CH3 <=> H + CH2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: O + CH4 <=> OH + CH3 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: O + HCO <=> OH + CO */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O + HCO <=> H + CO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: O + CH2O <=> OH + HCO */
    /*eqcon[26] *= 1;  */

    /*reaction 28: O + C2H2 <=> CH2(S) + CO */
    /*eqcon[27] *= 1;  */

    /*reaction 29: O + C2H2 <=> CO + CH2 */
    /*eqcon[28] *= 1;  */

    /*reaction 30: O + C2H4 <=> CH3 + HCO */
    /*eqcon[29] *= 1;  */

    /*reaction 31: O + C2H5 <=> CH3 + CH2O */
    /*eqcon[30] *= 1;  */

    /*reaction 32: O + C2H6 <=> OH + C2H5 */
    /*eqcon[31] *= 1;  */

    /*reaction 33: O2 + CO <=> O + CO2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: O2 + CH2O <=> HO2 + HCO */
    /*eqcon[33] *= 1;  */

    /*reaction 35: H + 2 O2 <=> HO2 + O2 */
    eqcon[34] *= 1e+06; 

    /*reaction 36: H + O2 + H2O <=> HO2 + H2O */
    eqcon[35] *= 1e+06; 

    /*reaction 37: H + O2 + N2 <=> HO2 + N2 */
    eqcon[36] *= 1e+06; 

    /*reaction 38: H + O2 + AR <=> HO2 + AR */
    eqcon[37] *= 1e+06; 

    /*reaction 39: H + O2 <=> O + OH */
    /*eqcon[38] *= 1;  */

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    eqcon[39] *= 1e+06; 

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    eqcon[40] *= 1e+06; 

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    eqcon[41] *= 1e+06; 

    /*reaction 43: H + HO2 <=> O2 + H2 */
    /*eqcon[42] *= 1;  */

    /*reaction 44: H + HO2 <=> 2 OH */
    /*eqcon[43] *= 1;  */

    /*reaction 45: H + H2O2 <=> HO2 + H2 */
    /*eqcon[44] *= 1;  */

    /*reaction 46: H + CH4 <=> CH3 + H2 */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + HCO <=> H2 + CO */
    /*eqcon[46] *= 1;  */

    /*reaction 48: H + CH2O <=> HCO + H2 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: H + CH3O <=> OH + CH3 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: H + C2H3 <=> H2 + C2H2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: H + C2H4 <=> C2H3 + H2 */
    /*eqcon[50] *= 1;  */

    /*reaction 52: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: OH + H2 <=> H + H2O */
    /*eqcon[52] *= 1;  */

    /*reaction 54: 2 OH <=> O + H2O */
    /*eqcon[53] *= 1;  */

    /*reaction 55: OH + HO2 <=> O2 + H2O */
    /*eqcon[54] *= 1;  */

    /*reaction 56: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: OH + CH2 <=> H + CH2O */
    /*eqcon[56] *= 1;  */

    /*reaction 58: OH + CH2(S) <=> H + CH2O */
    /*eqcon[57] *= 1;  */

    /*reaction 59: OH + CH3 <=> CH2 + H2O */
    /*eqcon[58] *= 1;  */

    /*reaction 60: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[59] *= 1;  */

    /*reaction 61: OH + CH4 <=> CH3 + H2O */
    /*eqcon[60] *= 1;  */

    /*reaction 62: OH + CO <=> H + CO2 */
    /*eqcon[61] *= 1;  */

    /*reaction 63: OH + HCO <=> H2O + CO */
    /*eqcon[62] *= 1;  */

    /*reaction 64: OH + CH2O <=> HCO + H2O */
    /*eqcon[63] *= 1;  */

    /*reaction 65: OH + C2H2 <=> CH3 + CO */
    /*eqcon[64] *= 1;  */

    /*reaction 66: OH + C2H3 <=> H2O + C2H2 */
    /*eqcon[65] *= 1;  */

    /*reaction 67: OH + C2H4 <=> C2H3 + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[67] *= 1;  */

    /*reaction 69: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[68] *= 1;  */

    /*reaction 70: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[69] *= 1;  */

    /*reaction 71: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[70] *= 1;  */

    /*reaction 72: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[71] *= 1;  */

    /*reaction 73: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[72] *= 1;  */

    /*reaction 74: HO2 + CO <=> OH + CO2 */
    /*eqcon[73] *= 1;  */

    /*reaction 75: HO2 + CH2O <=> HCO + H2O2 */
    /*eqcon[74] *= 1;  */

    /*reaction 76: CH2 + O2 <=> OH + HCO */
    /*eqcon[75] *= 1;  */

    /*reaction 77: CH2 + H2 <=> H + CH3 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: 2 CH2 <=> H2 + C2H2 */
    /*eqcon[77] *= 1;  */

    /*reaction 79: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[79] *= 1;  */

    /*reaction 81: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[80] *= 1;  */

    /*reaction 82: CH2(S) + AR <=> CH2 + AR */
    /*eqcon[81] *= 1;  */

    /*reaction 83: CH2(S) + O2 <=> H + OH + CO */
    eqcon[82] *= 1e-06; 

    /*reaction 84: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[83] *= 1;  */

    /*reaction 85: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[84] *= 1;  */

    /*reaction 86: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[85] *= 1;  */

    /*reaction 87: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[86] *= 1;  */

    /*reaction 88: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[87] *= 1;  */

    /*reaction 89: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[88] *= 1;  */

    /*reaction 90: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[89] *= 1;  */

    /*reaction 91: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[90] *= 1;  */

    /*reaction 92: CH3 + O2 <=> O + CH3O */
    /*eqcon[91] *= 1;  */

    /*reaction 93: CH3 + O2 <=> OH + CH2O */
    /*eqcon[92] *= 1;  */

    /*reaction 94: CH3 + H2O2 <=> HO2 + CH4 */
    /*eqcon[93] *= 1;  */

    /*reaction 95: 2 CH3 <=> H + C2H5 */
    /*eqcon[94] *= 1;  */

    /*reaction 96: CH3 + HCO <=> CH4 + CO */
    /*eqcon[95] *= 1;  */

    /*reaction 97: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[96] *= 1;  */

    /*reaction 98: CH3 + C2H4 <=> C2H3 + CH4 */
    /*eqcon[97] *= 1;  */

    /*reaction 99: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[98] *= 1;  */

    /*reaction 100: HCO + H2O <=> H + CO + H2O */
    eqcon[99] *= 1e-06; 

    /*reaction 101: HCO + O2 <=> HO2 + CO */
    /*eqcon[100] *= 1;  */

    /*reaction 102: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[101] *= 1;  */

    /*reaction 103: C2H3 + O2 <=> HCO + CH2O */
    /*eqcon[102] *= 1;  */

    /*reaction 104: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[103] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[24]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + C2H2 (+M) <=> C2H3 (+M) */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + C2H3 (+M) <=> C2H4 (+M) */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + C2H4 (+M) <=> C2H5 (+M) */
    eqcon[6] *= 1e+06; 

    /*reaction 8: H + C2H5 (+M) <=> C2H6 (+M) */
    eqcon[7] *= 1e+06; 

    /*reaction 9: H2 + CO (+M) <=> CH2O (+M) */
    eqcon[8] *= 1e+06; 

    /*reaction 10: 2 OH (+M) <=> H2O2 (+M) */
    eqcon[9] *= 1e+06; 

    /*reaction 11: 2 CH3 (+M) <=> C2H6 (+M) */
    eqcon[10] *= 1e+06; 

    /*reaction 12: C2H4 (+M) <=> H2 + C2H2 (+M) */
    eqcon[11] *= 1e-06; 

    /*reaction 13: O + H + M <=> OH + M */
    eqcon[12] *= 1e+06; 

    /*reaction 14: O + CO + M <=> CO2 + M */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + O2 + M <=> HO2 + M */
    eqcon[14] *= 1e+06; 

    /*reaction 16: 2 H + M <=> H2 + M */
    eqcon[15] *= 1e+06; 

    /*reaction 17: H + OH + M <=> H2O + M */
    eqcon[16] *= 1e+06; 

    /*reaction 18: HCO + M <=> H + CO + M */
    eqcon[17] *= 1e-06; 

    /*reaction 19: O + H2 <=> H + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: O + HO2 <=> OH + O2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: O + CH2 <=> H + HCO */
    /*eqcon[20] *= 1;  */

    /*reaction 22: O + CH2(S) <=> H + HCO */
    /*eqcon[21] *= 1;  */

    /*reaction 23: O + CH3 <=> H + CH2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: O + CH4 <=> OH + CH3 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: O + HCO <=> OH + CO */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O + HCO <=> H + CO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: O + CH2O <=> OH + HCO */
    /*eqcon[26] *= 1;  */

    /*reaction 28: O + C2H2 <=> CH2(S) + CO */
    /*eqcon[27] *= 1;  */

    /*reaction 29: O + C2H2 <=> CO + CH2 */
    /*eqcon[28] *= 1;  */

    /*reaction 30: O + C2H4 <=> CH3 + HCO */
    /*eqcon[29] *= 1;  */

    /*reaction 31: O + C2H5 <=> CH3 + CH2O */
    /*eqcon[30] *= 1;  */

    /*reaction 32: O + C2H6 <=> OH + C2H5 */
    /*eqcon[31] *= 1;  */

    /*reaction 33: O2 + CO <=> O + CO2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: O2 + CH2O <=> HO2 + HCO */
    /*eqcon[33] *= 1;  */

    /*reaction 35: H + 2 O2 <=> HO2 + O2 */
    eqcon[34] *= 1e+06; 

    /*reaction 36: H + O2 + H2O <=> HO2 + H2O */
    eqcon[35] *= 1e+06; 

    /*reaction 37: H + O2 + N2 <=> HO2 + N2 */
    eqcon[36] *= 1e+06; 

    /*reaction 38: H + O2 + AR <=> HO2 + AR */
    eqcon[37] *= 1e+06; 

    /*reaction 39: H + O2 <=> O + OH */
    /*eqcon[38] *= 1;  */

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    eqcon[39] *= 1e+06; 

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    eqcon[40] *= 1e+06; 

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    eqcon[41] *= 1e+06; 

    /*reaction 43: H + HO2 <=> O2 + H2 */
    /*eqcon[42] *= 1;  */

    /*reaction 44: H + HO2 <=> 2 OH */
    /*eqcon[43] *= 1;  */

    /*reaction 45: H + H2O2 <=> HO2 + H2 */
    /*eqcon[44] *= 1;  */

    /*reaction 46: H + CH4 <=> CH3 + H2 */
    /*eqcon[45] *= 1;  */

    /*reaction 47: H + HCO <=> H2 + CO */
    /*eqcon[46] *= 1;  */

    /*reaction 48: H + CH2O <=> HCO + H2 */
    /*eqcon[47] *= 1;  */

    /*reaction 49: H + CH3O <=> OH + CH3 */
    /*eqcon[48] *= 1;  */

    /*reaction 50: H + C2H3 <=> H2 + C2H2 */
    /*eqcon[49] *= 1;  */

    /*reaction 51: H + C2H4 <=> C2H3 + H2 */
    /*eqcon[50] *= 1;  */

    /*reaction 52: H + C2H6 <=> C2H5 + H2 */
    /*eqcon[51] *= 1;  */

    /*reaction 53: OH + H2 <=> H + H2O */
    /*eqcon[52] *= 1;  */

    /*reaction 54: 2 OH <=> O + H2O */
    /*eqcon[53] *= 1;  */

    /*reaction 55: OH + HO2 <=> O2 + H2O */
    /*eqcon[54] *= 1;  */

    /*reaction 56: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[55] *= 1;  */

    /*reaction 57: OH + CH2 <=> H + CH2O */
    /*eqcon[56] *= 1;  */

    /*reaction 58: OH + CH2(S) <=> H + CH2O */
    /*eqcon[57] *= 1;  */

    /*reaction 59: OH + CH3 <=> CH2 + H2O */
    /*eqcon[58] *= 1;  */

    /*reaction 60: OH + CH3 <=> CH2(S) + H2O */
    /*eqcon[59] *= 1;  */

    /*reaction 61: OH + CH4 <=> CH3 + H2O */
    /*eqcon[60] *= 1;  */

    /*reaction 62: OH + CO <=> H + CO2 */
    /*eqcon[61] *= 1;  */

    /*reaction 63: OH + HCO <=> H2O + CO */
    /*eqcon[62] *= 1;  */

    /*reaction 64: OH + CH2O <=> HCO + H2O */
    /*eqcon[63] *= 1;  */

    /*reaction 65: OH + C2H2 <=> CH3 + CO */
    /*eqcon[64] *= 1;  */

    /*reaction 66: OH + C2H3 <=> H2O + C2H2 */
    /*eqcon[65] *= 1;  */

    /*reaction 67: OH + C2H4 <=> C2H3 + H2O */
    /*eqcon[66] *= 1;  */

    /*reaction 68: OH + C2H6 <=> C2H5 + H2O */
    /*eqcon[67] *= 1;  */

    /*reaction 69: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[68] *= 1;  */

    /*reaction 70: 2 HO2 <=> O2 + H2O2 */
    /*eqcon[69] *= 1;  */

    /*reaction 71: HO2 + CH2 <=> OH + CH2O */
    /*eqcon[70] *= 1;  */

    /*reaction 72: HO2 + CH3 <=> O2 + CH4 */
    /*eqcon[71] *= 1;  */

    /*reaction 73: HO2 + CH3 <=> OH + CH3O */
    /*eqcon[72] *= 1;  */

    /*reaction 74: HO2 + CO <=> OH + CO2 */
    /*eqcon[73] *= 1;  */

    /*reaction 75: HO2 + CH2O <=> HCO + H2O2 */
    /*eqcon[74] *= 1;  */

    /*reaction 76: CH2 + O2 <=> OH + HCO */
    /*eqcon[75] *= 1;  */

    /*reaction 77: CH2 + H2 <=> H + CH3 */
    /*eqcon[76] *= 1;  */

    /*reaction 78: 2 CH2 <=> H2 + C2H2 */
    /*eqcon[77] *= 1;  */

    /*reaction 79: CH2 + CH3 <=> H + C2H4 */
    /*eqcon[78] *= 1;  */

    /*reaction 80: CH2 + CH4 <=> 2 CH3 */
    /*eqcon[79] *= 1;  */

    /*reaction 81: CH2(S) + N2 <=> CH2 + N2 */
    /*eqcon[80] *= 1;  */

    /*reaction 82: CH2(S) + AR <=> CH2 + AR */
    /*eqcon[81] *= 1;  */

    /*reaction 83: CH2(S) + O2 <=> H + OH + CO */
    eqcon[82] *= 1e-06; 

    /*reaction 84: CH2(S) + O2 <=> CO + H2O */
    /*eqcon[83] *= 1;  */

    /*reaction 85: CH2(S) + H2 <=> CH3 + H */
    /*eqcon[84] *= 1;  */

    /*reaction 86: CH2(S) + H2O <=> CH2 + H2O */
    /*eqcon[85] *= 1;  */

    /*reaction 87: CH2(S) + CH3 <=> H + C2H4 */
    /*eqcon[86] *= 1;  */

    /*reaction 88: CH2(S) + CH4 <=> 2 CH3 */
    /*eqcon[87] *= 1;  */

    /*reaction 89: CH2(S) + CO <=> CH2 + CO */
    /*eqcon[88] *= 1;  */

    /*reaction 90: CH2(S) + CO2 <=> CH2 + CO2 */
    /*eqcon[89] *= 1;  */

    /*reaction 91: CH2(S) + CO2 <=> CO + CH2O */
    /*eqcon[90] *= 1;  */

    /*reaction 92: CH3 + O2 <=> O + CH3O */
    /*eqcon[91] *= 1;  */

    /*reaction 93: CH3 + O2 <=> OH + CH2O */
    /*eqcon[92] *= 1;  */

    /*reaction 94: CH3 + H2O2 <=> HO2 + CH4 */
    /*eqcon[93] *= 1;  */

    /*reaction 95: 2 CH3 <=> H + C2H5 */
    /*eqcon[94] *= 1;  */

    /*reaction 96: CH3 + HCO <=> CH4 + CO */
    /*eqcon[95] *= 1;  */

    /*reaction 97: CH3 + CH2O <=> HCO + CH4 */
    /*eqcon[96] *= 1;  */

    /*reaction 98: CH3 + C2H4 <=> C2H3 + CH4 */
    /*eqcon[97] *= 1;  */

    /*reaction 99: CH3 + C2H6 <=> C2H5 + CH4 */
    /*eqcon[98] *= 1;  */

    /*reaction 100: HCO + H2O <=> H + CO + H2O */
    eqcon[99] *= 1e-06; 

    /*reaction 101: HCO + O2 <=> HO2 + CO */
    /*eqcon[100] *= 1;  */

    /*reaction 102: CH3O + O2 <=> HO2 + CH2O */
    /*eqcon[101] *= 1;  */

    /*reaction 103: C2H3 + O2 <=> HCO + CH2O */
    /*eqcon[102] *= 1;  */

    /*reaction 104: C2H5 + O2 <=> HO2 + C2H4 */
    /*eqcon[103] *= 1;  */
}

static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[104];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[104];
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

    double qdot, q_f[104], q_r[104];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 24; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[1] -= qdot;
    wdot[8] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[1] -= qdot;
    wdot[10] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[2]-q_r[2];
    wdot[1] -= qdot;
    wdot[14] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[1] -= qdot;
    wdot[15] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[1] -= qdot;
    wdot[17] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[1] -= qdot;
    wdot[18] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[6]-q_r[6];
    wdot[1] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[1] -= qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[0] -= qdot;
    wdot[12] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[9]-q_r[9];
    wdot[4] -= 2 * qdot;
    wdot[7] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[10] -= 2 * qdot;
    wdot[21] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[0] += qdot;
    wdot[17] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[12]-q_r[12];
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[13]-q_r[13];
    wdot[2] -= qdot;
    wdot[12] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[14]-q_r[14];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[0] += qdot;
    wdot[1] -= 2 * qdot;

    qdot = q_f[16]-q_r[16];
    wdot[1] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[17]-q_r[17];
    wdot[1] += qdot;
    wdot[12] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[18]-q_r[18];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[19]-q_r[19];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[8] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[21]-q_r[21];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[9] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[22]-q_r[22];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[10] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[23]-q_r[23];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[10] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[24]-q_r[24];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[12] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[25]-q_r[25];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[26]-q_r[26];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[14] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[27]-q_r[27];
    wdot[2] -= qdot;
    wdot[9] += qdot;
    wdot[12] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[28]-q_r[28];
    wdot[2] -= qdot;
    wdot[8] += qdot;
    wdot[12] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[29]-q_r[29];
    wdot[2] -= qdot;
    wdot[10] += qdot;
    wdot[14] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[30]-q_r[30];
    wdot[2] -= qdot;
    wdot[10] += qdot;
    wdot[15] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[31]-q_r[31];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[20] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[32]-q_r[32];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[12] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[33]-q_r[33];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[14] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[34]-q_r[34];
    wdot[1] -= qdot;
    wdot[3] -= 2 * qdot;
    wdot[3] += qdot;
    wdot[6] += qdot;

    qdot = q_f[35]-q_r[35];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;

    qdot = q_f[36]-q_r[36];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[22] -= qdot;
    wdot[22] += qdot;

    qdot = q_f[37]-q_r[37];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[23] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[38]-q_r[38];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[39]-q_r[39];
    wdot[0] -= qdot;
    wdot[0] += 2 * qdot;
    wdot[1] -= 2 * qdot;

    qdot = q_f[40]-q_r[40];
    wdot[0] += qdot;
    wdot[1] -= 2 * qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[41]-q_r[41];
    wdot[0] += qdot;
    wdot[1] -= 2 * qdot;
    wdot[13] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[42]-q_r[42];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[43]-q_r[43];
    wdot[1] -= qdot;
    wdot[4] += 2 * qdot;
    wdot[6] -= qdot;

    qdot = q_f[44]-q_r[44];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[45]-q_r[45];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[10] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[46]-q_r[46];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[12] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[47]-q_r[47];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[14] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[48]-q_r[48];
    wdot[1] -= qdot;
    wdot[4] += qdot;
    wdot[10] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[49]-q_r[49];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[50]-q_r[50];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[18] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[51]-q_r[51];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[20] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[52]-q_r[52];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[53]-q_r[53];
    wdot[2] += qdot;
    wdot[4] -= 2 * qdot;
    wdot[5] += qdot;

    qdot = q_f[54]-q_r[54];
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[55]-q_r[55];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[56]-q_r[56];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[8] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[57]-q_r[57];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[9] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[58]-q_r[58];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[8] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[59]-q_r[59];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[60]-q_r[60];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[10] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[61]-q_r[61];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[12] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[62]-q_r[62];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[12] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[63]-q_r[63];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[14] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[64]-q_r[64];
    wdot[4] -= qdot;
    wdot[10] += qdot;
    wdot[12] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[65]-q_r[65];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[66]-q_r[66];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[18] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[67]-q_r[67];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[20] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[68]-q_r[68];
    wdot[3] += qdot;
    wdot[6] -= 2 * qdot;
    wdot[7] += qdot;

    qdot = q_f[69]-q_r[69];
    wdot[3] += qdot;
    wdot[6] -= 2 * qdot;
    wdot[7] += qdot;

    qdot = q_f[70]-q_r[70];
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[8] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[71]-q_r[71];
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[10] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[72]-q_r[72];
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[10] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[73]-q_r[73];
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[12] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[74]-q_r[74];
    wdot[6] -= qdot;
    wdot[7] += qdot;
    wdot[14] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[75]-q_r[75];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[8] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[76]-q_r[76];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[8] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[77]-q_r[77];
    wdot[0] += qdot;
    wdot[8] -= 2 * qdot;
    wdot[17] += qdot;

    qdot = q_f[78]-q_r[78];
    wdot[1] += qdot;
    wdot[8] -= qdot;
    wdot[10] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[79]-q_r[79];
    wdot[8] -= qdot;
    wdot[10] += 2 * qdot;
    wdot[11] -= qdot;

    qdot = q_f[80]-q_r[80];
    wdot[8] += qdot;
    wdot[9] -= qdot;
    wdot[22] -= qdot;
    wdot[22] += qdot;

    qdot = q_f[81]-q_r[81];
    wdot[8] += qdot;
    wdot[9] -= qdot;
    wdot[23] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[82]-q_r[82];
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[9] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[83]-q_r[83];
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[9] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[84]-q_r[84];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[85]-q_r[85];
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[8] += qdot;
    wdot[9] -= qdot;

    qdot = q_f[86]-q_r[86];
    wdot[1] += qdot;
    wdot[9] -= qdot;
    wdot[10] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[87]-q_r[87];
    wdot[9] -= qdot;
    wdot[10] += 2 * qdot;
    wdot[11] -= qdot;

    qdot = q_f[88]-q_r[88];
    wdot[8] += qdot;
    wdot[9] -= qdot;
    wdot[12] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[89]-q_r[89];
    wdot[8] += qdot;
    wdot[9] -= qdot;
    wdot[13] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[90]-q_r[90];
    wdot[9] -= qdot;
    wdot[12] += qdot;
    wdot[13] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[91]-q_r[91];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[10] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[92]-q_r[92];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[10] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[93]-q_r[93];
    wdot[6] += qdot;
    wdot[7] -= qdot;
    wdot[10] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[94]-q_r[94];
    wdot[1] += qdot;
    wdot[10] -= 2 * qdot;
    wdot[20] += qdot;

    qdot = q_f[95]-q_r[95];
    wdot[10] -= qdot;
    wdot[11] += qdot;
    wdot[12] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[96]-q_r[96];
    wdot[10] -= qdot;
    wdot[11] += qdot;
    wdot[14] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[97]-q_r[97];
    wdot[10] -= qdot;
    wdot[11] += qdot;
    wdot[18] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[98]-q_r[98];
    wdot[10] -= qdot;
    wdot[11] += qdot;
    wdot[20] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[99]-q_r[99];
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[12] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[100]-q_r[100];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[12] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[101]-q_r[101];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[15] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[102]-q_r[102];
    wdot[3] -= qdot;
    wdot[14] += qdot;
    wdot[15] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[103]-q_r[103];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[19] += qdot;
    wdot[20] -= qdot;

    return;
}

void comp_k_f(double * restrict tc, double invT, double * restrict k_f)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<104; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(double * restrict tc, double invT, double * restrict Kc)
{
    /*compute the Gibbs free energy */
    double g_RT[24];
    gibbs(g_RT, tc);

    Kc[0] = g_RT[1] + g_RT[8] - g_RT[10];
    Kc[1] = g_RT[1] + g_RT[10] - g_RT[11];
    Kc[2] = g_RT[1] + g_RT[14] - g_RT[15];
    Kc[3] = g_RT[1] + g_RT[15] - g_RT[16];
    Kc[4] = g_RT[1] + g_RT[17] - g_RT[18];
    Kc[5] = g_RT[1] + g_RT[18] - g_RT[19];
    Kc[6] = g_RT[1] + g_RT[19] - g_RT[20];
    Kc[7] = g_RT[1] + g_RT[20] - g_RT[21];
    Kc[8] = g_RT[0] + g_RT[12] - g_RT[15];
    Kc[9] = 2*g_RT[4] - g_RT[7];
    Kc[10] = 2*g_RT[10] - g_RT[21];
    Kc[11] = -g_RT[0] - g_RT[17] + g_RT[19];
    Kc[12] = g_RT[1] + g_RT[2] - g_RT[4];
    Kc[13] = g_RT[2] + g_RT[12] - g_RT[13];
    Kc[14] = g_RT[1] + g_RT[3] - g_RT[6];
    Kc[15] = -g_RT[0] + 2*g_RT[1];
    Kc[16] = g_RT[1] + g_RT[4] - g_RT[5];
    Kc[17] = -g_RT[1] - g_RT[12] + g_RT[14];
    Kc[18] = g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4];
    Kc[19] = g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6];
    Kc[20] = -g_RT[1] + g_RT[2] + g_RT[8] - g_RT[14];
    Kc[21] = -g_RT[1] + g_RT[2] + g_RT[9] - g_RT[14];
    Kc[22] = -g_RT[1] + g_RT[2] + g_RT[10] - g_RT[15];
    Kc[23] = g_RT[2] - g_RT[4] - g_RT[10] + g_RT[11];
    Kc[24] = g_RT[2] - g_RT[4] - g_RT[12] + g_RT[14];
    Kc[25] = -g_RT[1] + g_RT[2] - g_RT[13] + g_RT[14];
    Kc[26] = g_RT[2] - g_RT[4] - g_RT[14] + g_RT[15];
    Kc[27] = g_RT[2] - g_RT[9] - g_RT[12] + g_RT[17];
    Kc[28] = g_RT[2] - g_RT[8] - g_RT[12] + g_RT[17];
    Kc[29] = g_RT[2] - g_RT[10] - g_RT[14] + g_RT[19];
    Kc[30] = g_RT[2] - g_RT[10] - g_RT[15] + g_RT[20];
    Kc[31] = g_RT[2] - g_RT[4] - g_RT[20] + g_RT[21];
    Kc[32] = -g_RT[2] + g_RT[3] + g_RT[12] - g_RT[13];
    Kc[33] = g_RT[3] - g_RT[6] - g_RT[14] + g_RT[15];
    Kc[34] = g_RT[1] + 2*g_RT[3] - g_RT[3] - g_RT[6];
    Kc[35] = g_RT[1] + g_RT[3] + g_RT[5] - g_RT[5] - g_RT[6];
    Kc[36] = g_RT[1] + g_RT[3] - g_RT[6] + g_RT[22] - g_RT[22];
    Kc[37] = g_RT[1] + g_RT[3] - g_RT[6] + g_RT[23] - g_RT[23];
    Kc[38] = g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4];
    Kc[39] = g_RT[0] - 2*g_RT[0] + 2*g_RT[1];
    Kc[40] = -g_RT[0] + 2*g_RT[1] + g_RT[5] - g_RT[5];
    Kc[41] = -g_RT[0] + 2*g_RT[1] + g_RT[13] - g_RT[13];
    Kc[42] = -g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6];
    Kc[43] = g_RT[1] - 2*g_RT[4] + g_RT[6];
    Kc[44] = -g_RT[0] + g_RT[1] - g_RT[6] + g_RT[7];
    Kc[45] = -g_RT[0] + g_RT[1] - g_RT[10] + g_RT[11];
    Kc[46] = -g_RT[0] + g_RT[1] - g_RT[12] + g_RT[14];
    Kc[47] = -g_RT[0] + g_RT[1] - g_RT[14] + g_RT[15];
    Kc[48] = g_RT[1] - g_RT[4] - g_RT[10] + g_RT[16];
    Kc[49] = -g_RT[0] + g_RT[1] - g_RT[17] + g_RT[18];
    Kc[50] = -g_RT[0] + g_RT[1] - g_RT[18] + g_RT[19];
    Kc[51] = -g_RT[0] + g_RT[1] - g_RT[20] + g_RT[21];
    Kc[52] = g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5];
    Kc[53] = -g_RT[2] + 2*g_RT[4] - g_RT[5];
    Kc[54] = -g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6];
    Kc[55] = g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7];
    Kc[56] = -g_RT[1] + g_RT[4] + g_RT[8] - g_RT[15];
    Kc[57] = -g_RT[1] + g_RT[4] + g_RT[9] - g_RT[15];
    Kc[58] = g_RT[4] - g_RT[5] - g_RT[8] + g_RT[10];
    Kc[59] = g_RT[4] - g_RT[5] - g_RT[9] + g_RT[10];
    Kc[60] = g_RT[4] - g_RT[5] - g_RT[10] + g_RT[11];
    Kc[61] = -g_RT[1] + g_RT[4] + g_RT[12] - g_RT[13];
    Kc[62] = g_RT[4] - g_RT[5] - g_RT[12] + g_RT[14];
    Kc[63] = g_RT[4] - g_RT[5] - g_RT[14] + g_RT[15];
    Kc[64] = g_RT[4] - g_RT[10] - g_RT[12] + g_RT[17];
    Kc[65] = g_RT[4] - g_RT[5] - g_RT[17] + g_RT[18];
    Kc[66] = g_RT[4] - g_RT[5] - g_RT[18] + g_RT[19];
    Kc[67] = g_RT[4] - g_RT[5] - g_RT[20] + g_RT[21];
    Kc[68] = -g_RT[3] + 2*g_RT[6] - g_RT[7];
    Kc[69] = -g_RT[3] + 2*g_RT[6] - g_RT[7];
    Kc[70] = -g_RT[4] + g_RT[6] + g_RT[8] - g_RT[15];
    Kc[71] = -g_RT[3] + g_RT[6] + g_RT[10] - g_RT[11];
    Kc[72] = -g_RT[4] + g_RT[6] + g_RT[10] - g_RT[16];
    Kc[73] = -g_RT[4] + g_RT[6] + g_RT[12] - g_RT[13];
    Kc[74] = g_RT[6] - g_RT[7] - g_RT[14] + g_RT[15];
    Kc[75] = g_RT[3] - g_RT[4] + g_RT[8] - g_RT[14];
    Kc[76] = g_RT[0] - g_RT[1] + g_RT[8] - g_RT[10];
    Kc[77] = -g_RT[0] + 2*g_RT[8] - g_RT[17];
    Kc[78] = -g_RT[1] + g_RT[8] + g_RT[10] - g_RT[19];
    Kc[79] = g_RT[8] - 2*g_RT[10] + g_RT[11];
    Kc[80] = -g_RT[8] + g_RT[9] + g_RT[22] - g_RT[22];
    Kc[81] = -g_RT[8] + g_RT[9] + g_RT[23] - g_RT[23];
    Kc[82] = -g_RT[1] + g_RT[3] - g_RT[4] + g_RT[9] - g_RT[12];
    Kc[83] = g_RT[3] - g_RT[5] + g_RT[9] - g_RT[12];
    Kc[84] = g_RT[0] - g_RT[1] + g_RT[9] - g_RT[10];
    Kc[85] = g_RT[5] - g_RT[5] - g_RT[8] + g_RT[9];
    Kc[86] = -g_RT[1] + g_RT[9] + g_RT[10] - g_RT[19];
    Kc[87] = g_RT[9] - 2*g_RT[10] + g_RT[11];
    Kc[88] = -g_RT[8] + g_RT[9] + g_RT[12] - g_RT[12];
    Kc[89] = -g_RT[8] + g_RT[9] + g_RT[13] - g_RT[13];
    Kc[90] = g_RT[9] - g_RT[12] + g_RT[13] - g_RT[15];
    Kc[91] = -g_RT[2] + g_RT[3] + g_RT[10] - g_RT[16];
    Kc[92] = g_RT[3] - g_RT[4] + g_RT[10] - g_RT[15];
    Kc[93] = -g_RT[6] + g_RT[7] + g_RT[10] - g_RT[11];
    Kc[94] = -g_RT[1] + 2*g_RT[10] - g_RT[20];
    Kc[95] = g_RT[10] - g_RT[11] - g_RT[12] + g_RT[14];
    Kc[96] = g_RT[10] - g_RT[11] - g_RT[14] + g_RT[15];
    Kc[97] = g_RT[10] - g_RT[11] - g_RT[18] + g_RT[19];
    Kc[98] = g_RT[10] - g_RT[11] - g_RT[20] + g_RT[21];
    Kc[99] = -g_RT[1] + g_RT[5] - g_RT[5] - g_RT[12] + g_RT[14];
    Kc[100] = g_RT[3] - g_RT[6] - g_RT[12] + g_RT[14];
    Kc[101] = g_RT[3] - g_RT[6] - g_RT[15] + g_RT[16];
    Kc[102] = g_RT[3] - g_RT[14] - g_RT[15] + g_RT[18];
    Kc[103] = g_RT[3] - g_RT[6] - g_RT[19] + g_RT[20];

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<104; ++i) {
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
    Kc[7] *= refCinv;
    Kc[8] *= refCinv;
    Kc[9] *= refCinv;
    Kc[10] *= refCinv;
    Kc[11] *= refC;
    Kc[12] *= refCinv;
    Kc[13] *= refCinv;
    Kc[14] *= refCinv;
    Kc[15] *= refCinv;
    Kc[16] *= refCinv;
    Kc[17] *= refC;
    Kc[34] *= refCinv;
    Kc[35] *= refCinv;
    Kc[36] *= refCinv;
    Kc[37] *= refCinv;
    Kc[39] *= refCinv;
    Kc[40] *= refCinv;
    Kc[41] *= refCinv;
    Kc[82] *= refC;
    Kc[99] *= refC;

    return;
}

void comp_qfqr(double * restrict qf, double * restrict qr, double * restrict sc, double * restrict tc, double invT)
{

    /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
    qf[0] = sc[1]*sc[8];
    qr[0] = sc[10];

    /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
    qf[1] = sc[1]*sc[10];
    qr[1] = sc[11];

    /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
    qf[2] = sc[1]*sc[14];
    qr[2] = sc[15];

    /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
    qf[3] = sc[1]*sc[15];
    qr[3] = sc[16];

    /*reaction 5: H + C2H2 (+M) <=> C2H3 (+M) */
    qf[4] = sc[1]*sc[17];
    qr[4] = sc[18];

    /*reaction 6: H + C2H3 (+M) <=> C2H4 (+M) */
    qf[5] = sc[1]*sc[18];
    qr[5] = sc[19];

    /*reaction 7: H + C2H4 (+M) <=> C2H5 (+M) */
    qf[6] = sc[1]*sc[19];
    qr[6] = sc[20];

    /*reaction 8: H + C2H5 (+M) <=> C2H6 (+M) */
    qf[7] = sc[1]*sc[20];
    qr[7] = sc[21];

    /*reaction 9: H2 + CO (+M) <=> CH2O (+M) */
    qf[8] = sc[0]*sc[12];
    qr[8] = sc[15];

    /*reaction 10: 2 OH (+M) <=> H2O2 (+M) */
    qf[9] = sc[4]*sc[4];
    qr[9] = sc[7];

    /*reaction 11: 2 CH3 (+M) <=> C2H6 (+M) */
    qf[10] = sc[10]*sc[10];
    qr[10] = sc[21];

    /*reaction 12: C2H4 (+M) <=> H2 + C2H2 (+M) */
    qf[11] = sc[19];
    qr[11] = sc[0]*sc[17];

    /*reaction 13: O + H + M <=> OH + M */
    qf[12] = sc[1]*sc[2];
    qr[12] = sc[4];

    /*reaction 14: O + CO + M <=> CO2 + M */
    qf[13] = sc[2]*sc[12];
    qr[13] = sc[13];

    /*reaction 15: H + O2 + M <=> HO2 + M */
    qf[14] = sc[1]*sc[3];
    qr[14] = sc[6];

    /*reaction 16: 2 H + M <=> H2 + M */
    qf[15] = sc[1]*sc[1];
    qr[15] = sc[0];

    /*reaction 17: H + OH + M <=> H2O + M */
    qf[16] = sc[1]*sc[4];
    qr[16] = sc[5];

    /*reaction 18: HCO + M <=> H + CO + M */
    qf[17] = sc[14];
    qr[17] = sc[1]*sc[12];

    /*reaction 19: O + H2 <=> H + OH */
    qf[18] = sc[0]*sc[2];
    qr[18] = sc[1]*sc[4];

    /*reaction 20: O + HO2 <=> OH + O2 */
    qf[19] = sc[2]*sc[6];
    qr[19] = sc[3]*sc[4];

    /*reaction 21: O + CH2 <=> H + HCO */
    qf[20] = sc[2]*sc[8];
    qr[20] = sc[1]*sc[14];

    /*reaction 22: O + CH2(S) <=> H + HCO */
    qf[21] = sc[2]*sc[9];
    qr[21] = sc[1]*sc[14];

    /*reaction 23: O + CH3 <=> H + CH2O */
    qf[22] = sc[2]*sc[10];
    qr[22] = sc[1]*sc[15];

    /*reaction 24: O + CH4 <=> OH + CH3 */
    qf[23] = sc[2]*sc[11];
    qr[23] = sc[4]*sc[10];

    /*reaction 25: O + HCO <=> OH + CO */
    qf[24] = sc[2]*sc[14];
    qr[24] = sc[4]*sc[12];

    /*reaction 26: O + HCO <=> H + CO2 */
    qf[25] = sc[2]*sc[14];
    qr[25] = sc[1]*sc[13];

    /*reaction 27: O + CH2O <=> OH + HCO */
    qf[26] = sc[2]*sc[15];
    qr[26] = sc[4]*sc[14];

    /*reaction 28: O + C2H2 <=> CH2(S) + CO */
    qf[27] = sc[2]*sc[17];
    qr[27] = sc[9]*sc[12];

    /*reaction 29: O + C2H2 <=> CO + CH2 */
    qf[28] = sc[2]*sc[17];
    qr[28] = sc[8]*sc[12];

    /*reaction 30: O + C2H4 <=> CH3 + HCO */
    qf[29] = sc[2]*sc[19];
    qr[29] = sc[10]*sc[14];

    /*reaction 31: O + C2H5 <=> CH3 + CH2O */
    qf[30] = sc[2]*sc[20];
    qr[30] = sc[10]*sc[15];

    /*reaction 32: O + C2H6 <=> OH + C2H5 */
    qf[31] = sc[2]*sc[21];
    qr[31] = sc[4]*sc[20];

    /*reaction 33: O2 + CO <=> O + CO2 */
    qf[32] = sc[3]*sc[12];
    qr[32] = sc[2]*sc[13];

    /*reaction 34: O2 + CH2O <=> HO2 + HCO */
    qf[33] = sc[3]*sc[15];
    qr[33] = sc[6]*sc[14];

    /*reaction 35: H + 2 O2 <=> HO2 + O2 */
    qf[34] = sc[1]*sc[3]*sc[3];
    qr[34] = sc[3]*sc[6];

    /*reaction 36: H + O2 + H2O <=> HO2 + H2O */
    qf[35] = sc[1]*sc[3]*sc[5];
    qr[35] = sc[5]*sc[6];

    /*reaction 37: H + O2 + N2 <=> HO2 + N2 */
    qf[36] = sc[1]*sc[3]*sc[22];
    qr[36] = sc[6]*sc[22];

    /*reaction 38: H + O2 + AR <=> HO2 + AR */
    qf[37] = sc[1]*sc[3]*sc[23];
    qr[37] = sc[6]*sc[23];

    /*reaction 39: H + O2 <=> O + OH */
    qf[38] = sc[1]*sc[3];
    qr[38] = sc[2]*sc[4];

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    qf[39] = sc[0]*sc[1]*sc[1];
    qr[39] = sc[0]*sc[0];

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    qf[40] = sc[1]*sc[1]*sc[5];
    qr[40] = sc[0]*sc[5];

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    qf[41] = sc[1]*sc[1]*sc[13];
    qr[41] = sc[0]*sc[13];

    /*reaction 43: H + HO2 <=> O2 + H2 */
    qf[42] = sc[1]*sc[6];
    qr[42] = sc[0]*sc[3];

    /*reaction 44: H + HO2 <=> 2 OH */
    qf[43] = sc[1]*sc[6];
    qr[43] = sc[4]*sc[4];

    /*reaction 45: H + H2O2 <=> HO2 + H2 */
    qf[44] = sc[1]*sc[7];
    qr[44] = sc[0]*sc[6];

    /*reaction 46: H + CH4 <=> CH3 + H2 */
    qf[45] = sc[1]*sc[11];
    qr[45] = sc[0]*sc[10];

    /*reaction 47: H + HCO <=> H2 + CO */
    qf[46] = sc[1]*sc[14];
    qr[46] = sc[0]*sc[12];

    /*reaction 48: H + CH2O <=> HCO + H2 */
    qf[47] = sc[1]*sc[15];
    qr[47] = sc[0]*sc[14];

    /*reaction 49: H + CH3O <=> OH + CH3 */
    qf[48] = sc[1]*sc[16];
    qr[48] = sc[4]*sc[10];

    /*reaction 50: H + C2H3 <=> H2 + C2H2 */
    qf[49] = sc[1]*sc[18];
    qr[49] = sc[0]*sc[17];

    /*reaction 51: H + C2H4 <=> C2H3 + H2 */
    qf[50] = sc[1]*sc[19];
    qr[50] = sc[0]*sc[18];

    /*reaction 52: H + C2H6 <=> C2H5 + H2 */
    qf[51] = sc[1]*sc[21];
    qr[51] = sc[0]*sc[20];

    /*reaction 53: OH + H2 <=> H + H2O */
    qf[52] = sc[0]*sc[4];
    qr[52] = sc[1]*sc[5];

    /*reaction 54: 2 OH <=> O + H2O */
    qf[53] = sc[4]*sc[4];
    qr[53] = sc[2]*sc[5];

    /*reaction 55: OH + HO2 <=> O2 + H2O */
    qf[54] = sc[4]*sc[6];
    qr[54] = sc[3]*sc[5];

    /*reaction 56: OH + H2O2 <=> HO2 + H2O */
    qf[55] = sc[4]*sc[7];
    qr[55] = sc[5]*sc[6];

    /*reaction 57: OH + CH2 <=> H + CH2O */
    qf[56] = sc[4]*sc[8];
    qr[56] = sc[1]*sc[15];

    /*reaction 58: OH + CH2(S) <=> H + CH2O */
    qf[57] = sc[4]*sc[9];
    qr[57] = sc[1]*sc[15];

    /*reaction 59: OH + CH3 <=> CH2 + H2O */
    qf[58] = sc[4]*sc[10];
    qr[58] = sc[5]*sc[8];

    /*reaction 60: OH + CH3 <=> CH2(S) + H2O */
    qf[59] = sc[4]*sc[10];
    qr[59] = sc[5]*sc[9];

    /*reaction 61: OH + CH4 <=> CH3 + H2O */
    qf[60] = sc[4]*sc[11];
    qr[60] = sc[5]*sc[10];

    /*reaction 62: OH + CO <=> H + CO2 */
    qf[61] = sc[4]*sc[12];
    qr[61] = sc[1]*sc[13];

    /*reaction 63: OH + HCO <=> H2O + CO */
    qf[62] = sc[4]*sc[14];
    qr[62] = sc[5]*sc[12];

    /*reaction 64: OH + CH2O <=> HCO + H2O */
    qf[63] = sc[4]*sc[15];
    qr[63] = sc[5]*sc[14];

    /*reaction 65: OH + C2H2 <=> CH3 + CO */
    qf[64] = sc[4]*sc[17];
    qr[64] = sc[10]*sc[12];

    /*reaction 66: OH + C2H3 <=> H2O + C2H2 */
    qf[65] = sc[4]*sc[18];
    qr[65] = sc[5]*sc[17];

    /*reaction 67: OH + C2H4 <=> C2H3 + H2O */
    qf[66] = sc[4]*sc[19];
    qr[66] = sc[5]*sc[18];

    /*reaction 68: OH + C2H6 <=> C2H5 + H2O */
    qf[67] = sc[4]*sc[21];
    qr[67] = sc[5]*sc[20];

    /*reaction 69: 2 HO2 <=> O2 + H2O2 */
    qf[68] = sc[6]*sc[6];
    qr[68] = sc[3]*sc[7];

    /*reaction 70: 2 HO2 <=> O2 + H2O2 */
    qf[69] = sc[6]*sc[6];
    qr[69] = sc[3]*sc[7];

    /*reaction 71: HO2 + CH2 <=> OH + CH2O */
    qf[70] = sc[6]*sc[8];
    qr[70] = sc[4]*sc[15];

    /*reaction 72: HO2 + CH3 <=> O2 + CH4 */
    qf[71] = sc[6]*sc[10];
    qr[71] = sc[3]*sc[11];

    /*reaction 73: HO2 + CH3 <=> OH + CH3O */
    qf[72] = sc[6]*sc[10];
    qr[72] = sc[4]*sc[16];

    /*reaction 74: HO2 + CO <=> OH + CO2 */
    qf[73] = sc[6]*sc[12];
    qr[73] = sc[4]*sc[13];

    /*reaction 75: HO2 + CH2O <=> HCO + H2O2 */
    qf[74] = sc[6]*sc[15];
    qr[74] = sc[7]*sc[14];

    /*reaction 76: CH2 + O2 <=> OH + HCO */
    qf[75] = sc[3]*sc[8];
    qr[75] = sc[4]*sc[14];

    /*reaction 77: CH2 + H2 <=> H + CH3 */
    qf[76] = sc[0]*sc[8];
    qr[76] = sc[1]*sc[10];

    /*reaction 78: 2 CH2 <=> H2 + C2H2 */
    qf[77] = sc[8]*sc[8];
    qr[77] = sc[0]*sc[17];

    /*reaction 79: CH2 + CH3 <=> H + C2H4 */
    qf[78] = sc[8]*sc[10];
    qr[78] = sc[1]*sc[19];

    /*reaction 80: CH2 + CH4 <=> 2 CH3 */
    qf[79] = sc[8]*sc[11];
    qr[79] = sc[10]*sc[10];

    /*reaction 81: CH2(S) + N2 <=> CH2 + N2 */
    qf[80] = sc[9]*sc[22];
    qr[80] = sc[8]*sc[22];

    /*reaction 82: CH2(S) + AR <=> CH2 + AR */
    qf[81] = sc[9]*sc[23];
    qr[81] = sc[8]*sc[23];

    /*reaction 83: CH2(S) + O2 <=> H + OH + CO */
    qf[82] = sc[3]*sc[9];
    qr[82] = sc[1]*sc[4]*sc[12];

    /*reaction 84: CH2(S) + O2 <=> CO + H2O */
    qf[83] = sc[3]*sc[9];
    qr[83] = sc[5]*sc[12];

    /*reaction 85: CH2(S) + H2 <=> CH3 + H */
    qf[84] = sc[0]*sc[9];
    qr[84] = sc[1]*sc[10];

    /*reaction 86: CH2(S) + H2O <=> CH2 + H2O */
    qf[85] = sc[5]*sc[9];
    qr[85] = sc[5]*sc[8];

    /*reaction 87: CH2(S) + CH3 <=> H + C2H4 */
    qf[86] = sc[9]*sc[10];
    qr[86] = sc[1]*sc[19];

    /*reaction 88: CH2(S) + CH4 <=> 2 CH3 */
    qf[87] = sc[9]*sc[11];
    qr[87] = sc[10]*sc[10];

    /*reaction 89: CH2(S) + CO <=> CH2 + CO */
    qf[88] = sc[9]*sc[12];
    qr[88] = sc[8]*sc[12];

    /*reaction 90: CH2(S) + CO2 <=> CH2 + CO2 */
    qf[89] = sc[9]*sc[13];
    qr[89] = sc[8]*sc[13];

    /*reaction 91: CH2(S) + CO2 <=> CO + CH2O */
    qf[90] = sc[9]*sc[13];
    qr[90] = sc[12]*sc[15];

    /*reaction 92: CH3 + O2 <=> O + CH3O */
    qf[91] = sc[3]*sc[10];
    qr[91] = sc[2]*sc[16];

    /*reaction 93: CH3 + O2 <=> OH + CH2O */
    qf[92] = sc[3]*sc[10];
    qr[92] = sc[4]*sc[15];

    /*reaction 94: CH3 + H2O2 <=> HO2 + CH4 */
    qf[93] = sc[7]*sc[10];
    qr[93] = sc[6]*sc[11];

    /*reaction 95: 2 CH3 <=> H + C2H5 */
    qf[94] = sc[10]*sc[10];
    qr[94] = sc[1]*sc[20];

    /*reaction 96: CH3 + HCO <=> CH4 + CO */
    qf[95] = sc[10]*sc[14];
    qr[95] = sc[11]*sc[12];

    /*reaction 97: CH3 + CH2O <=> HCO + CH4 */
    qf[96] = sc[10]*sc[15];
    qr[96] = sc[11]*sc[14];

    /*reaction 98: CH3 + C2H4 <=> C2H3 + CH4 */
    qf[97] = sc[10]*sc[19];
    qr[97] = sc[11]*sc[18];

    /*reaction 99: CH3 + C2H6 <=> C2H5 + CH4 */
    qf[98] = sc[10]*sc[21];
    qr[98] = sc[11]*sc[20];

    /*reaction 100: HCO + H2O <=> H + CO + H2O */
    qf[99] = sc[5]*sc[14];
    qr[99] = sc[1]*sc[5]*sc[12];

    /*reaction 101: HCO + O2 <=> HO2 + CO */
    qf[100] = sc[3]*sc[14];
    qr[100] = sc[6]*sc[12];

    /*reaction 102: CH3O + O2 <=> HO2 + CH2O */
    qf[101] = sc[3]*sc[16];
    qr[101] = sc[6]*sc[15];

    /*reaction 103: C2H3 + O2 <=> HCO + CH2O */
    qf[102] = sc[3]*sc[18];
    qr[102] = sc[14]*sc[15];

    /*reaction 104: C2H5 + O2 <=> HO2 + C2H4 */
    qf[103] = sc[3]*sc[20];
    qr[103] = sc[6]*sc[19];

    double T = tc[1];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 24; ++i) {
        mixture += sc[i];
    }

    double Corr[104];
    for (int i = 0; i < 104; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        double alpha[12];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[5] + (TB[0][2] - 1)*sc[11] + (TB[0][3] - 1)*sc[12] + (TB[0][4] - 1)*sc[13] + (TB[0][5] - 1)*sc[21] + (TB[0][6] - 1)*sc[23];
        alpha[1] = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[5] + (TB[1][2] - 1)*sc[11] + (TB[1][3] - 1)*sc[12] + (TB[1][4] - 1)*sc[13] + (TB[1][5] - 1)*sc[21] + (TB[1][6] - 1)*sc[23];
        alpha[2] = mixture + (TB[2][0] - 1)*sc[0] + (TB[2][1] - 1)*sc[5] + (TB[2][2] - 1)*sc[11] + (TB[2][3] - 1)*sc[12] + (TB[2][4] - 1)*sc[13] + (TB[2][5] - 1)*sc[21] + (TB[2][6] - 1)*sc[23];
        alpha[3] = mixture + (TB[3][0] - 1)*sc[0] + (TB[3][1] - 1)*sc[5] + (TB[3][2] - 1)*sc[11] + (TB[3][3] - 1)*sc[12] + (TB[3][4] - 1)*sc[13] + (TB[3][5] - 1)*sc[21];
        alpha[4] = mixture + (TB[4][0] - 1)*sc[0] + (TB[4][1] - 1)*sc[5] + (TB[4][2] - 1)*sc[11] + (TB[4][3] - 1)*sc[12] + (TB[4][4] - 1)*sc[13] + (TB[4][5] - 1)*sc[21] + (TB[4][6] - 1)*sc[23];
        alpha[5] = mixture + (TB[5][0] - 1)*sc[0] + (TB[5][1] - 1)*sc[5] + (TB[5][2] - 1)*sc[11] + (TB[5][3] - 1)*sc[12] + (TB[5][4] - 1)*sc[13] + (TB[5][5] - 1)*sc[21] + (TB[5][6] - 1)*sc[23];
        alpha[6] = mixture + (TB[6][0] - 1)*sc[0] + (TB[6][1] - 1)*sc[5] + (TB[6][2] - 1)*sc[11] + (TB[6][3] - 1)*sc[12] + (TB[6][4] - 1)*sc[13] + (TB[6][5] - 1)*sc[21] + (TB[6][6] - 1)*sc[23];
        alpha[7] = mixture + (TB[7][0] - 1)*sc[0] + (TB[7][1] - 1)*sc[5] + (TB[7][2] - 1)*sc[11] + (TB[7][3] - 1)*sc[12] + (TB[7][4] - 1)*sc[13] + (TB[7][5] - 1)*sc[21] + (TB[7][6] - 1)*sc[23];
        alpha[8] = mixture + (TB[8][0] - 1)*sc[0] + (TB[8][1] - 1)*sc[5] + (TB[8][2] - 1)*sc[11] + (TB[8][3] - 1)*sc[12] + (TB[8][4] - 1)*sc[13] + (TB[8][5] - 1)*sc[21] + (TB[8][6] - 1)*sc[23];
        alpha[9] = mixture + (TB[9][0] - 1)*sc[0] + (TB[9][1] - 1)*sc[5] + (TB[9][2] - 1)*sc[11] + (TB[9][3] - 1)*sc[12] + (TB[9][4] - 1)*sc[13] + (TB[9][5] - 1)*sc[21] + (TB[9][6] - 1)*sc[23];
        alpha[10] = mixture + (TB[10][0] - 1)*sc[0] + (TB[10][1] - 1)*sc[5] + (TB[10][2] - 1)*sc[11] + (TB[10][3] - 1)*sc[12] + (TB[10][4] - 1)*sc[13] + (TB[10][5] - 1)*sc[21] + (TB[10][6] - 1)*sc[23];
        alpha[11] = mixture + (TB[11][0] - 1)*sc[0] + (TB[11][1] - 1)*sc[5] + (TB[11][2] - 1)*sc[11] + (TB[11][3] - 1)*sc[12] + (TB[11][4] - 1)*sc[13] + (TB[11][5] - 1)*sc[21] + (TB[11][6] - 1)*sc[23];
#ifdef __INTEL_COMPILER
         #pragma simd
#endif
        for (int i=0; i<12; i++)
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
        alpha = mixture + (TB[12][0] - 1)*sc[0] + (TB[12][1] - 1)*sc[5] + (TB[12][2] - 1)*sc[11] + (TB[12][3] - 1)*sc[12] + (TB[12][4] - 1)*sc[13] + (TB[12][5] - 1)*sc[21] + (TB[12][6] - 1)*sc[23];
        Corr[12] = alpha;
        alpha = mixture + (TB[13][0] - 1)*sc[0] + (TB[13][1] - 1)*sc[3] + (TB[13][2] - 1)*sc[5] + (TB[13][3] - 1)*sc[11] + (TB[13][4] - 1)*sc[12] + (TB[13][5] - 1)*sc[13] + (TB[13][6] - 1)*sc[21] + (TB[13][7] - 1)*sc[23];
        Corr[13] = alpha;
        alpha = mixture + (TB[14][0] - 1)*sc[3] + (TB[14][1] - 1)*sc[5] + (TB[14][2] - 1)*sc[12] + (TB[14][3] - 1)*sc[13] + (TB[14][4] - 1)*sc[21] + (TB[14][5] - 1)*sc[22] + (TB[14][6] - 1)*sc[23];
        Corr[14] = alpha;
        alpha = mixture + (TB[15][0] - 1)*sc[0] + (TB[15][1] - 1)*sc[5] + (TB[15][2] - 1)*sc[11] + (TB[15][3] - 1)*sc[13] + (TB[15][4] - 1)*sc[21] + (TB[15][5] - 1)*sc[23];
        Corr[15] = alpha;
        alpha = mixture + (TB[16][0] - 1)*sc[0] + (TB[16][1] - 1)*sc[5] + (TB[16][2] - 1)*sc[11] + (TB[16][3] - 1)*sc[21] + (TB[16][4] - 1)*sc[23];
        Corr[16] = alpha;
        alpha = mixture + (TB[17][0] - 1)*sc[0] + (TB[17][1] - 1)*sc[5] + (TB[17][2] - 1)*sc[11] + (TB[17][3] - 1)*sc[12] + (TB[17][4] - 1)*sc[13] + (TB[17][5] - 1)*sc[21];
        Corr[17] = alpha;
    }

    for (int i=0; i<104; i++)
    {
        qf[i] *= Corr[i] * k_f_save[i];
        qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];
    }

    return;
}


/*compute the production rate for each species */
void vproductionRate(int npt, double * restrict wdot, double * restrict sc, double * restrict T)
{
    double k_f_s[104*npt], Kc_s[104*npt], mixture[npt], g_RT[24*npt];
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

    for (int n=0; n<24; n++) {
        for (int i=0; i<npt; i++) {
            mixture[i] += sc[n*npt+i];
            wdot[n*npt+i] = 0.0;
        }
    }

    vcomp_k_f(npt, k_f_s, tc, invT);

    vcomp_gibbs(npt, g_RT, tc);

    vcomp_Kc(npt, Kc_s, g_RT, invT);

    vcomp_wdot_1_50(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
    vcomp_wdot_51_100(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
    vcomp_wdot_101_104(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
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
        k_f_s[38*npt+i] = prefactor_units[38] * fwd_A[38] * exp(fwd_beta[38] * tc[i] - activation_units[38] * fwd_Ea[38] * invT[i]);
        k_f_s[39*npt+i] = prefactor_units[39] * fwd_A[39] * exp(fwd_beta[39] * tc[i] - activation_units[39] * fwd_Ea[39] * invT[i]);
        k_f_s[40*npt+i] = prefactor_units[40] * fwd_A[40] * exp(fwd_beta[40] * tc[i] - activation_units[40] * fwd_Ea[40] * invT[i]);
        k_f_s[41*npt+i] = prefactor_units[41] * fwd_A[41] * exp(fwd_beta[41] * tc[i] - activation_units[41] * fwd_Ea[41] * invT[i]);
        k_f_s[42*npt+i] = prefactor_units[42] * fwd_A[42] * exp(fwd_beta[42] * tc[i] - activation_units[42] * fwd_Ea[42] * invT[i]);
        k_f_s[43*npt+i] = prefactor_units[43] * fwd_A[43] * exp(fwd_beta[43] * tc[i] - activation_units[43] * fwd_Ea[43] * invT[i]);
        k_f_s[44*npt+i] = prefactor_units[44] * fwd_A[44] * exp(fwd_beta[44] * tc[i] - activation_units[44] * fwd_Ea[44] * invT[i]);
        k_f_s[45*npt+i] = prefactor_units[45] * fwd_A[45] * exp(fwd_beta[45] * tc[i] - activation_units[45] * fwd_Ea[45] * invT[i]);
        k_f_s[46*npt+i] = prefactor_units[46] * fwd_A[46] * exp(fwd_beta[46] * tc[i] - activation_units[46] * fwd_Ea[46] * invT[i]);
        k_f_s[47*npt+i] = prefactor_units[47] * fwd_A[47] * exp(fwd_beta[47] * tc[i] - activation_units[47] * fwd_Ea[47] * invT[i]);
        k_f_s[48*npt+i] = prefactor_units[48] * fwd_A[48] * exp(fwd_beta[48] * tc[i] - activation_units[48] * fwd_Ea[48] * invT[i]);
        k_f_s[49*npt+i] = prefactor_units[49] * fwd_A[49] * exp(fwd_beta[49] * tc[i] - activation_units[49] * fwd_Ea[49] * invT[i]);
        k_f_s[50*npt+i] = prefactor_units[50] * fwd_A[50] * exp(fwd_beta[50] * tc[i] - activation_units[50] * fwd_Ea[50] * invT[i]);
        k_f_s[51*npt+i] = prefactor_units[51] * fwd_A[51] * exp(fwd_beta[51] * tc[i] - activation_units[51] * fwd_Ea[51] * invT[i]);
        k_f_s[52*npt+i] = prefactor_units[52] * fwd_A[52] * exp(fwd_beta[52] * tc[i] - activation_units[52] * fwd_Ea[52] * invT[i]);
        k_f_s[53*npt+i] = prefactor_units[53] * fwd_A[53] * exp(fwd_beta[53] * tc[i] - activation_units[53] * fwd_Ea[53] * invT[i]);
        k_f_s[54*npt+i] = prefactor_units[54] * fwd_A[54] * exp(fwd_beta[54] * tc[i] - activation_units[54] * fwd_Ea[54] * invT[i]);
        k_f_s[55*npt+i] = prefactor_units[55] * fwd_A[55] * exp(fwd_beta[55] * tc[i] - activation_units[55] * fwd_Ea[55] * invT[i]);
        k_f_s[56*npt+i] = prefactor_units[56] * fwd_A[56] * exp(fwd_beta[56] * tc[i] - activation_units[56] * fwd_Ea[56] * invT[i]);
        k_f_s[57*npt+i] = prefactor_units[57] * fwd_A[57] * exp(fwd_beta[57] * tc[i] - activation_units[57] * fwd_Ea[57] * invT[i]);
        k_f_s[58*npt+i] = prefactor_units[58] * fwd_A[58] * exp(fwd_beta[58] * tc[i] - activation_units[58] * fwd_Ea[58] * invT[i]);
        k_f_s[59*npt+i] = prefactor_units[59] * fwd_A[59] * exp(fwd_beta[59] * tc[i] - activation_units[59] * fwd_Ea[59] * invT[i]);
        k_f_s[60*npt+i] = prefactor_units[60] * fwd_A[60] * exp(fwd_beta[60] * tc[i] - activation_units[60] * fwd_Ea[60] * invT[i]);
        k_f_s[61*npt+i] = prefactor_units[61] * fwd_A[61] * exp(fwd_beta[61] * tc[i] - activation_units[61] * fwd_Ea[61] * invT[i]);
        k_f_s[62*npt+i] = prefactor_units[62] * fwd_A[62] * exp(fwd_beta[62] * tc[i] - activation_units[62] * fwd_Ea[62] * invT[i]);
        k_f_s[63*npt+i] = prefactor_units[63] * fwd_A[63] * exp(fwd_beta[63] * tc[i] - activation_units[63] * fwd_Ea[63] * invT[i]);
        k_f_s[64*npt+i] = prefactor_units[64] * fwd_A[64] * exp(fwd_beta[64] * tc[i] - activation_units[64] * fwd_Ea[64] * invT[i]);
        k_f_s[65*npt+i] = prefactor_units[65] * fwd_A[65] * exp(fwd_beta[65] * tc[i] - activation_units[65] * fwd_Ea[65] * invT[i]);
        k_f_s[66*npt+i] = prefactor_units[66] * fwd_A[66] * exp(fwd_beta[66] * tc[i] - activation_units[66] * fwd_Ea[66] * invT[i]);
        k_f_s[67*npt+i] = prefactor_units[67] * fwd_A[67] * exp(fwd_beta[67] * tc[i] - activation_units[67] * fwd_Ea[67] * invT[i]);
        k_f_s[68*npt+i] = prefactor_units[68] * fwd_A[68] * exp(fwd_beta[68] * tc[i] - activation_units[68] * fwd_Ea[68] * invT[i]);
        k_f_s[69*npt+i] = prefactor_units[69] * fwd_A[69] * exp(fwd_beta[69] * tc[i] - activation_units[69] * fwd_Ea[69] * invT[i]);
        k_f_s[70*npt+i] = prefactor_units[70] * fwd_A[70] * exp(fwd_beta[70] * tc[i] - activation_units[70] * fwd_Ea[70] * invT[i]);
        k_f_s[71*npt+i] = prefactor_units[71] * fwd_A[71] * exp(fwd_beta[71] * tc[i] - activation_units[71] * fwd_Ea[71] * invT[i]);
        k_f_s[72*npt+i] = prefactor_units[72] * fwd_A[72] * exp(fwd_beta[72] * tc[i] - activation_units[72] * fwd_Ea[72] * invT[i]);
        k_f_s[73*npt+i] = prefactor_units[73] * fwd_A[73] * exp(fwd_beta[73] * tc[i] - activation_units[73] * fwd_Ea[73] * invT[i]);
        k_f_s[74*npt+i] = prefactor_units[74] * fwd_A[74] * exp(fwd_beta[74] * tc[i] - activation_units[74] * fwd_Ea[74] * invT[i]);
        k_f_s[75*npt+i] = prefactor_units[75] * fwd_A[75] * exp(fwd_beta[75] * tc[i] - activation_units[75] * fwd_Ea[75] * invT[i]);
        k_f_s[76*npt+i] = prefactor_units[76] * fwd_A[76] * exp(fwd_beta[76] * tc[i] - activation_units[76] * fwd_Ea[76] * invT[i]);
        k_f_s[77*npt+i] = prefactor_units[77] * fwd_A[77] * exp(fwd_beta[77] * tc[i] - activation_units[77] * fwd_Ea[77] * invT[i]);
        k_f_s[78*npt+i] = prefactor_units[78] * fwd_A[78] * exp(fwd_beta[78] * tc[i] - activation_units[78] * fwd_Ea[78] * invT[i]);
        k_f_s[79*npt+i] = prefactor_units[79] * fwd_A[79] * exp(fwd_beta[79] * tc[i] - activation_units[79] * fwd_Ea[79] * invT[i]);
        k_f_s[80*npt+i] = prefactor_units[80] * fwd_A[80] * exp(fwd_beta[80] * tc[i] - activation_units[80] * fwd_Ea[80] * invT[i]);
        k_f_s[81*npt+i] = prefactor_units[81] * fwd_A[81] * exp(fwd_beta[81] * tc[i] - activation_units[81] * fwd_Ea[81] * invT[i]);
        k_f_s[82*npt+i] = prefactor_units[82] * fwd_A[82] * exp(fwd_beta[82] * tc[i] - activation_units[82] * fwd_Ea[82] * invT[i]);
        k_f_s[83*npt+i] = prefactor_units[83] * fwd_A[83] * exp(fwd_beta[83] * tc[i] - activation_units[83] * fwd_Ea[83] * invT[i]);
        k_f_s[84*npt+i] = prefactor_units[84] * fwd_A[84] * exp(fwd_beta[84] * tc[i] - activation_units[84] * fwd_Ea[84] * invT[i]);
        k_f_s[85*npt+i] = prefactor_units[85] * fwd_A[85] * exp(fwd_beta[85] * tc[i] - activation_units[85] * fwd_Ea[85] * invT[i]);
        k_f_s[86*npt+i] = prefactor_units[86] * fwd_A[86] * exp(fwd_beta[86] * tc[i] - activation_units[86] * fwd_Ea[86] * invT[i]);
        k_f_s[87*npt+i] = prefactor_units[87] * fwd_A[87] * exp(fwd_beta[87] * tc[i] - activation_units[87] * fwd_Ea[87] * invT[i]);
        k_f_s[88*npt+i] = prefactor_units[88] * fwd_A[88] * exp(fwd_beta[88] * tc[i] - activation_units[88] * fwd_Ea[88] * invT[i]);
        k_f_s[89*npt+i] = prefactor_units[89] * fwd_A[89] * exp(fwd_beta[89] * tc[i] - activation_units[89] * fwd_Ea[89] * invT[i]);
        k_f_s[90*npt+i] = prefactor_units[90] * fwd_A[90] * exp(fwd_beta[90] * tc[i] - activation_units[90] * fwd_Ea[90] * invT[i]);
        k_f_s[91*npt+i] = prefactor_units[91] * fwd_A[91] * exp(fwd_beta[91] * tc[i] - activation_units[91] * fwd_Ea[91] * invT[i]);
        k_f_s[92*npt+i] = prefactor_units[92] * fwd_A[92] * exp(fwd_beta[92] * tc[i] - activation_units[92] * fwd_Ea[92] * invT[i]);
        k_f_s[93*npt+i] = prefactor_units[93] * fwd_A[93] * exp(fwd_beta[93] * tc[i] - activation_units[93] * fwd_Ea[93] * invT[i]);
        k_f_s[94*npt+i] = prefactor_units[94] * fwd_A[94] * exp(fwd_beta[94] * tc[i] - activation_units[94] * fwd_Ea[94] * invT[i]);
        k_f_s[95*npt+i] = prefactor_units[95] * fwd_A[95] * exp(fwd_beta[95] * tc[i] - activation_units[95] * fwd_Ea[95] * invT[i]);
        k_f_s[96*npt+i] = prefactor_units[96] * fwd_A[96] * exp(fwd_beta[96] * tc[i] - activation_units[96] * fwd_Ea[96] * invT[i]);
        k_f_s[97*npt+i] = prefactor_units[97] * fwd_A[97] * exp(fwd_beta[97] * tc[i] - activation_units[97] * fwd_Ea[97] * invT[i]);
        k_f_s[98*npt+i] = prefactor_units[98] * fwd_A[98] * exp(fwd_beta[98] * tc[i] - activation_units[98] * fwd_Ea[98] * invT[i]);
        k_f_s[99*npt+i] = prefactor_units[99] * fwd_A[99] * exp(fwd_beta[99] * tc[i] - activation_units[99] * fwd_Ea[99] * invT[i]);
        k_f_s[100*npt+i] = prefactor_units[100] * fwd_A[100] * exp(fwd_beta[100] * tc[i] - activation_units[100] * fwd_Ea[100] * invT[i]);
        k_f_s[101*npt+i] = prefactor_units[101] * fwd_A[101] * exp(fwd_beta[101] * tc[i] - activation_units[101] * fwd_Ea[101] * invT[i]);
        k_f_s[102*npt+i] = prefactor_units[102] * fwd_A[102] * exp(fwd_beta[102] * tc[i] - activation_units[102] * fwd_Ea[102] * invT[i]);
        k_f_s[103*npt+i] = prefactor_units[103] * fwd_A[103] * exp(fwd_beta[103] * tc[i] - activation_units[103] * fwd_Ea[103] * invT[i]);
    }
}

void vcomp_gibbs(int npt, double * restrict g_RT, double * restrict tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        double tg[5], g[24];
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
        g_RT[14*npt+i] = g[14];
        g_RT[15*npt+i] = g[15];
        g_RT[16*npt+i] = g[16];
        g_RT[17*npt+i] = g[17];
        g_RT[18*npt+i] = g[18];
        g_RT[19*npt+i] = g[19];
        g_RT[20*npt+i] = g[20];
        g_RT[21*npt+i] = g[21];
        g_RT[22*npt+i] = g[22];
        g_RT[23*npt+i] = g[23];
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

        Kc_s[0*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[8*npt+i]) - (g_RT[10*npt+i]));
        Kc_s[1*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[10*npt+i]) - (g_RT[11*npt+i]));
        Kc_s[2*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[14*npt+i]) - (g_RT[15*npt+i]));
        Kc_s[3*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[15*npt+i]) - (g_RT[16*npt+i]));
        Kc_s[4*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[17*npt+i]) - (g_RT[18*npt+i]));
        Kc_s[5*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[18*npt+i]) - (g_RT[19*npt+i]));
        Kc_s[6*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[19*npt+i]) - (g_RT[20*npt+i]));
        Kc_s[7*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[20*npt+i]) - (g_RT[21*npt+i]));
        Kc_s[8*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[12*npt+i]) - (g_RT[15*npt+i]));
        Kc_s[9*npt+i] = refCinv * exp((2 * g_RT[4*npt+i]) - (g_RT[7*npt+i]));
        Kc_s[10*npt+i] = refCinv * exp((2 * g_RT[10*npt+i]) - (g_RT[21*npt+i]));
        Kc_s[11*npt+i] = refC * exp((g_RT[19*npt+i]) - (g_RT[0*npt+i] + g_RT[17*npt+i]));
        Kc_s[12*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[2*npt+i]) - (g_RT[4*npt+i]));
        Kc_s[13*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[12*npt+i]) - (g_RT[13*npt+i]));
        Kc_s[14*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[6*npt+i]));
        Kc_s[15*npt+i] = refCinv * exp((2 * g_RT[1*npt+i]) - (g_RT[0*npt+i]));
        Kc_s[16*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[4*npt+i]) - (g_RT[5*npt+i]));
        Kc_s[17*npt+i] = refC * exp((g_RT[14*npt+i]) - (g_RT[1*npt+i] + g_RT[12*npt+i]));
        Kc_s[18*npt+i] = exp((g_RT[0*npt+i] + g_RT[2*npt+i]) - (g_RT[1*npt+i] + g_RT[4*npt+i]));
        Kc_s[19*npt+i] = exp((g_RT[2*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[4*npt+i]));
        Kc_s[20*npt+i] = exp((g_RT[2*npt+i] + g_RT[8*npt+i]) - (g_RT[1*npt+i] + g_RT[14*npt+i]));
        Kc_s[21*npt+i] = exp((g_RT[2*npt+i] + g_RT[9*npt+i]) - (g_RT[1*npt+i] + g_RT[14*npt+i]));
        Kc_s[22*npt+i] = exp((g_RT[2*npt+i] + g_RT[10*npt+i]) - (g_RT[1*npt+i] + g_RT[15*npt+i]));
        Kc_s[23*npt+i] = exp((g_RT[2*npt+i] + g_RT[11*npt+i]) - (g_RT[4*npt+i] + g_RT[10*npt+i]));
        Kc_s[24*npt+i] = exp((g_RT[2*npt+i] + g_RT[14*npt+i]) - (g_RT[4*npt+i] + g_RT[12*npt+i]));
        Kc_s[25*npt+i] = exp((g_RT[2*npt+i] + g_RT[14*npt+i]) - (g_RT[1*npt+i] + g_RT[13*npt+i]));
        Kc_s[26*npt+i] = exp((g_RT[2*npt+i] + g_RT[15*npt+i]) - (g_RT[4*npt+i] + g_RT[14*npt+i]));
        Kc_s[27*npt+i] = exp((g_RT[2*npt+i] + g_RT[17*npt+i]) - (g_RT[9*npt+i] + g_RT[12*npt+i]));
        Kc_s[28*npt+i] = exp((g_RT[2*npt+i] + g_RT[17*npt+i]) - (g_RT[8*npt+i] + g_RT[12*npt+i]));
        Kc_s[29*npt+i] = exp((g_RT[2*npt+i] + g_RT[19*npt+i]) - (g_RT[10*npt+i] + g_RT[14*npt+i]));
        Kc_s[30*npt+i] = exp((g_RT[2*npt+i] + g_RT[20*npt+i]) - (g_RT[10*npt+i] + g_RT[15*npt+i]));
        Kc_s[31*npt+i] = exp((g_RT[2*npt+i] + g_RT[21*npt+i]) - (g_RT[4*npt+i] + g_RT[20*npt+i]));
        Kc_s[32*npt+i] = exp((g_RT[3*npt+i] + g_RT[12*npt+i]) - (g_RT[2*npt+i] + g_RT[13*npt+i]));
        Kc_s[33*npt+i] = exp((g_RT[3*npt+i] + g_RT[15*npt+i]) - (g_RT[6*npt+i] + g_RT[14*npt+i]));
        Kc_s[34*npt+i] = refCinv * exp((g_RT[1*npt+i] + 2 * g_RT[3*npt+i]) - (g_RT[3*npt+i] + g_RT[6*npt+i]));
        Kc_s[35*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i] + g_RT[5*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[36*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i] + g_RT[22*npt+i]) - (g_RT[6*npt+i] + g_RT[22*npt+i]));
        Kc_s[37*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i] + g_RT[23*npt+i]) - (g_RT[6*npt+i] + g_RT[23*npt+i]));
        Kc_s[38*npt+i] = exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[2*npt+i] + g_RT[4*npt+i]));
        Kc_s[39*npt+i] = refCinv * exp((g_RT[0*npt+i] + 2 * g_RT[1*npt+i]) - (2 * g_RT[0*npt+i]));
        Kc_s[40*npt+i] = refCinv * exp((2 * g_RT[1*npt+i] + g_RT[5*npt+i]) - (g_RT[0*npt+i] + g_RT[5*npt+i]));
        Kc_s[41*npt+i] = refCinv * exp((2 * g_RT[1*npt+i] + g_RT[13*npt+i]) - (g_RT[0*npt+i] + g_RT[13*npt+i]));
        Kc_s[42*npt+i] = exp((g_RT[1*npt+i] + g_RT[6*npt+i]) - (g_RT[0*npt+i] + g_RT[3*npt+i]));
        Kc_s[43*npt+i] = exp((g_RT[1*npt+i] + g_RT[6*npt+i]) - (2 * g_RT[4*npt+i]));
        Kc_s[44*npt+i] = exp((g_RT[1*npt+i] + g_RT[7*npt+i]) - (g_RT[0*npt+i] + g_RT[6*npt+i]));
        Kc_s[45*npt+i] = exp((g_RT[1*npt+i] + g_RT[11*npt+i]) - (g_RT[0*npt+i] + g_RT[10*npt+i]));
        Kc_s[46*npt+i] = exp((g_RT[1*npt+i] + g_RT[14*npt+i]) - (g_RT[0*npt+i] + g_RT[12*npt+i]));
        Kc_s[47*npt+i] = exp((g_RT[1*npt+i] + g_RT[15*npt+i]) - (g_RT[0*npt+i] + g_RT[14*npt+i]));
        Kc_s[48*npt+i] = exp((g_RT[1*npt+i] + g_RT[16*npt+i]) - (g_RT[4*npt+i] + g_RT[10*npt+i]));
        Kc_s[49*npt+i] = exp((g_RT[1*npt+i] + g_RT[18*npt+i]) - (g_RT[0*npt+i] + g_RT[17*npt+i]));
        Kc_s[50*npt+i] = exp((g_RT[1*npt+i] + g_RT[19*npt+i]) - (g_RT[0*npt+i] + g_RT[18*npt+i]));
        Kc_s[51*npt+i] = exp((g_RT[1*npt+i] + g_RT[21*npt+i]) - (g_RT[0*npt+i] + g_RT[20*npt+i]));
        Kc_s[52*npt+i] = exp((g_RT[0*npt+i] + g_RT[4*npt+i]) - (g_RT[1*npt+i] + g_RT[5*npt+i]));
        Kc_s[53*npt+i] = exp((2 * g_RT[4*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i]));
        Kc_s[54*npt+i] = exp((g_RT[4*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[5*npt+i]));
        Kc_s[55*npt+i] = exp((g_RT[4*npt+i] + g_RT[7*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[56*npt+i] = exp((g_RT[4*npt+i] + g_RT[8*npt+i]) - (g_RT[1*npt+i] + g_RT[15*npt+i]));
        Kc_s[57*npt+i] = exp((g_RT[4*npt+i] + g_RT[9*npt+i]) - (g_RT[1*npt+i] + g_RT[15*npt+i]));
        Kc_s[58*npt+i] = exp((g_RT[4*npt+i] + g_RT[10*npt+i]) - (g_RT[5*npt+i] + g_RT[8*npt+i]));
        Kc_s[59*npt+i] = exp((g_RT[4*npt+i] + g_RT[10*npt+i]) - (g_RT[5*npt+i] + g_RT[9*npt+i]));
        Kc_s[60*npt+i] = exp((g_RT[4*npt+i] + g_RT[11*npt+i]) - (g_RT[5*npt+i] + g_RT[10*npt+i]));
        Kc_s[61*npt+i] = exp((g_RT[4*npt+i] + g_RT[12*npt+i]) - (g_RT[1*npt+i] + g_RT[13*npt+i]));
        Kc_s[62*npt+i] = exp((g_RT[4*npt+i] + g_RT[14*npt+i]) - (g_RT[5*npt+i] + g_RT[12*npt+i]));
        Kc_s[63*npt+i] = exp((g_RT[4*npt+i] + g_RT[15*npt+i]) - (g_RT[5*npt+i] + g_RT[14*npt+i]));
        Kc_s[64*npt+i] = exp((g_RT[4*npt+i] + g_RT[17*npt+i]) - (g_RT[10*npt+i] + g_RT[12*npt+i]));
        Kc_s[65*npt+i] = exp((g_RT[4*npt+i] + g_RT[18*npt+i]) - (g_RT[5*npt+i] + g_RT[17*npt+i]));
        Kc_s[66*npt+i] = exp((g_RT[4*npt+i] + g_RT[19*npt+i]) - (g_RT[5*npt+i] + g_RT[18*npt+i]));
        Kc_s[67*npt+i] = exp((g_RT[4*npt+i] + g_RT[21*npt+i]) - (g_RT[5*npt+i] + g_RT[20*npt+i]));
        Kc_s[68*npt+i] = exp((2 * g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[7*npt+i]));
        Kc_s[69*npt+i] = exp((2 * g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[7*npt+i]));
        Kc_s[70*npt+i] = exp((g_RT[6*npt+i] + g_RT[8*npt+i]) - (g_RT[4*npt+i] + g_RT[15*npt+i]));
        Kc_s[71*npt+i] = exp((g_RT[6*npt+i] + g_RT[10*npt+i]) - (g_RT[3*npt+i] + g_RT[11*npt+i]));
        Kc_s[72*npt+i] = exp((g_RT[6*npt+i] + g_RT[10*npt+i]) - (g_RT[4*npt+i] + g_RT[16*npt+i]));
        Kc_s[73*npt+i] = exp((g_RT[6*npt+i] + g_RT[12*npt+i]) - (g_RT[4*npt+i] + g_RT[13*npt+i]));
        Kc_s[74*npt+i] = exp((g_RT[6*npt+i] + g_RT[15*npt+i]) - (g_RT[7*npt+i] + g_RT[14*npt+i]));
        Kc_s[75*npt+i] = exp((g_RT[3*npt+i] + g_RT[8*npt+i]) - (g_RT[4*npt+i] + g_RT[14*npt+i]));
        Kc_s[76*npt+i] = exp((g_RT[0*npt+i] + g_RT[8*npt+i]) - (g_RT[1*npt+i] + g_RT[10*npt+i]));
        Kc_s[77*npt+i] = exp((2 * g_RT[8*npt+i]) - (g_RT[0*npt+i] + g_RT[17*npt+i]));
        Kc_s[78*npt+i] = exp((g_RT[8*npt+i] + g_RT[10*npt+i]) - (g_RT[1*npt+i] + g_RT[19*npt+i]));
        Kc_s[79*npt+i] = exp((g_RT[8*npt+i] + g_RT[11*npt+i]) - (2 * g_RT[10*npt+i]));
        Kc_s[80*npt+i] = exp((g_RT[9*npt+i] + g_RT[22*npt+i]) - (g_RT[8*npt+i] + g_RT[22*npt+i]));
        Kc_s[81*npt+i] = exp((g_RT[9*npt+i] + g_RT[23*npt+i]) - (g_RT[8*npt+i] + g_RT[23*npt+i]));
        Kc_s[82*npt+i] = refC * exp((g_RT[3*npt+i] + g_RT[9*npt+i]) - (g_RT[1*npt+i] + g_RT[4*npt+i] + g_RT[12*npt+i]));
        Kc_s[83*npt+i] = exp((g_RT[3*npt+i] + g_RT[9*npt+i]) - (g_RT[5*npt+i] + g_RT[12*npt+i]));
        Kc_s[84*npt+i] = exp((g_RT[0*npt+i] + g_RT[9*npt+i]) - (g_RT[1*npt+i] + g_RT[10*npt+i]));
        Kc_s[85*npt+i] = exp((g_RT[5*npt+i] + g_RT[9*npt+i]) - (g_RT[5*npt+i] + g_RT[8*npt+i]));
        Kc_s[86*npt+i] = exp((g_RT[9*npt+i] + g_RT[10*npt+i]) - (g_RT[1*npt+i] + g_RT[19*npt+i]));
        Kc_s[87*npt+i] = exp((g_RT[9*npt+i] + g_RT[11*npt+i]) - (2 * g_RT[10*npt+i]));
        Kc_s[88*npt+i] = exp((g_RT[9*npt+i] + g_RT[12*npt+i]) - (g_RT[8*npt+i] + g_RT[12*npt+i]));
        Kc_s[89*npt+i] = exp((g_RT[9*npt+i] + g_RT[13*npt+i]) - (g_RT[8*npt+i] + g_RT[13*npt+i]));
        Kc_s[90*npt+i] = exp((g_RT[9*npt+i] + g_RT[13*npt+i]) - (g_RT[12*npt+i] + g_RT[15*npt+i]));
        Kc_s[91*npt+i] = exp((g_RT[3*npt+i] + g_RT[10*npt+i]) - (g_RT[2*npt+i] + g_RT[16*npt+i]));
        Kc_s[92*npt+i] = exp((g_RT[3*npt+i] + g_RT[10*npt+i]) - (g_RT[4*npt+i] + g_RT[15*npt+i]));
        Kc_s[93*npt+i] = exp((g_RT[7*npt+i] + g_RT[10*npt+i]) - (g_RT[6*npt+i] + g_RT[11*npt+i]));
        Kc_s[94*npt+i] = exp((2 * g_RT[10*npt+i]) - (g_RT[1*npt+i] + g_RT[20*npt+i]));
        Kc_s[95*npt+i] = exp((g_RT[10*npt+i] + g_RT[14*npt+i]) - (g_RT[11*npt+i] + g_RT[12*npt+i]));
        Kc_s[96*npt+i] = exp((g_RT[10*npt+i] + g_RT[15*npt+i]) - (g_RT[11*npt+i] + g_RT[14*npt+i]));
        Kc_s[97*npt+i] = exp((g_RT[10*npt+i] + g_RT[19*npt+i]) - (g_RT[11*npt+i] + g_RT[18*npt+i]));
        Kc_s[98*npt+i] = exp((g_RT[10*npt+i] + g_RT[21*npt+i]) - (g_RT[11*npt+i] + g_RT[20*npt+i]));
        Kc_s[99*npt+i] = refC * exp((g_RT[5*npt+i] + g_RT[14*npt+i]) - (g_RT[1*npt+i] + g_RT[5*npt+i] + g_RT[12*npt+i]));
        Kc_s[100*npt+i] = exp((g_RT[3*npt+i] + g_RT[14*npt+i]) - (g_RT[6*npt+i] + g_RT[12*npt+i]));
        Kc_s[101*npt+i] = exp((g_RT[3*npt+i] + g_RT[16*npt+i]) - (g_RT[6*npt+i] + g_RT[15*npt+i]));
        Kc_s[102*npt+i] = exp((g_RT[3*npt+i] + g_RT[18*npt+i]) - (g_RT[14*npt+i] + g_RT[15*npt+i]));
        Kc_s[103*npt+i] = exp((g_RT[3*npt+i] + g_RT[20*npt+i]) - (g_RT[6*npt+i] + g_RT[19*npt+i]));
    }
}

void vcomp_wdot_1_50(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,
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

        /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
        phi_f = sc[1*npt+i]*sc[8*npt+i];
        alpha = mixture[i] + (TB[0][0] - 1)*sc[0*npt+i] + (TB[0][1] - 1)*sc[5*npt+i] + (TB[0][2] - 1)*sc[11*npt+i] + (TB[0][3] - 1)*sc[12*npt+i] + (TB[0][4] - 1)*sc[13*npt+i] + (TB[0][5] - 1)*sc[21*npt+i] + (TB[0][6] - 1)*sc[23*npt+i];
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
        phi_r = sc[10*npt+i];
        Kc = Kc_s[0*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
        phi_f = sc[1*npt+i]*sc[10*npt+i];
        alpha = mixture[i] + (TB[1][0] - 1)*sc[0*npt+i] + (TB[1][1] - 1)*sc[5*npt+i] + (TB[1][2] - 1)*sc[11*npt+i] + (TB[1][3] - 1)*sc[12*npt+i] + (TB[1][4] - 1)*sc[13*npt+i] + (TB[1][5] - 1)*sc[21*npt+i] + (TB[1][6] - 1)*sc[23*npt+i];
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
        phi_r = sc[11*npt+i];
        Kc = Kc_s[1*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
        phi_f = sc[1*npt+i]*sc[14*npt+i];
        alpha = mixture[i] + (TB[2][0] - 1)*sc[0*npt+i] + (TB[2][1] - 1)*sc[5*npt+i] + (TB[2][2] - 1)*sc[11*npt+i] + (TB[2][3] - 1)*sc[12*npt+i] + (TB[2][4] - 1)*sc[13*npt+i] + (TB[2][5] - 1)*sc[21*npt+i] + (TB[2][6] - 1)*sc[23*npt+i];
        k_f = k_f_s[2*npt+i];
        redP = alpha / k_f * phase_units[2] * low_A[2] * exp(low_beta[2] * tc[i] - activation_units[2] * low_Ea[2] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[2]) > 1.e-100 ? (1.-troe_a[2])*exp(-T[i]/troe_Tsss[2]) : 0.) 
            + (fabs(troe_Ts[2]) > 1.e-100 ? troe_a[2] * exp(-T[i]/troe_Ts[2]) : 0.) 
            + (troe_len[2] == 4 ? exp(-troe_Tss[2] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[15*npt+i];
        Kc = Kc_s[2*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
        phi_f = sc[1*npt+i]*sc[15*npt+i];
        alpha = mixture[i] + (TB[3][0] - 1)*sc[0*npt+i] + (TB[3][1] - 1)*sc[5*npt+i] + (TB[3][2] - 1)*sc[11*npt+i] + (TB[3][3] - 1)*sc[12*npt+i] + (TB[3][4] - 1)*sc[13*npt+i] + (TB[3][5] - 1)*sc[21*npt+i];
        k_f = k_f_s[3*npt+i];
        redP = alpha / k_f * phase_units[3] * low_A[3] * exp(low_beta[3] * tc[i] - activation_units[3] * low_Ea[3] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[3]) > 1.e-100 ? (1.-troe_a[3])*exp(-T[i]/troe_Tsss[3]) : 0.) 
            + (fabs(troe_Ts[3]) > 1.e-100 ? troe_a[3] * exp(-T[i]/troe_Ts[3]) : 0.) 
            + (troe_len[3] == 4 ? exp(-troe_Tss[3] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[16*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[15*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 5: H + C2H2 (+M) <=> C2H3 (+M) */
        phi_f = sc[1*npt+i]*sc[17*npt+i];
        alpha = mixture[i] + (TB[4][0] - 1)*sc[0*npt+i] + (TB[4][1] - 1)*sc[5*npt+i] + (TB[4][2] - 1)*sc[11*npt+i] + (TB[4][3] - 1)*sc[12*npt+i] + (TB[4][4] - 1)*sc[13*npt+i] + (TB[4][5] - 1)*sc[21*npt+i] + (TB[4][6] - 1)*sc[23*npt+i];
        k_f = k_f_s[4*npt+i];
        redP = alpha / k_f * phase_units[4] * low_A[4] * exp(low_beta[4] * tc[i] - activation_units[4] * low_Ea[4] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[4]) > 1.e-100 ? (1.-troe_a[4])*exp(-T[i]/troe_Tsss[4]) : 0.) 
            + (fabs(troe_Ts[4]) > 1.e-100 ? troe_a[4] * exp(-T[i]/troe_Ts[4]) : 0.) 
            + (troe_len[4] == 4 ? exp(-troe_Tss[4] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[18*npt+i];
        Kc = Kc_s[4*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[17*npt+i] -= qdot;
        wdot[18*npt+i] += qdot;

        /*reaction 6: H + C2H3 (+M) <=> C2H4 (+M) */
        phi_f = sc[1*npt+i]*sc[18*npt+i];
        alpha = mixture[i] + (TB[5][0] - 1)*sc[0*npt+i] + (TB[5][1] - 1)*sc[5*npt+i] + (TB[5][2] - 1)*sc[11*npt+i] + (TB[5][3] - 1)*sc[12*npt+i] + (TB[5][4] - 1)*sc[13*npt+i] + (TB[5][5] - 1)*sc[21*npt+i] + (TB[5][6] - 1)*sc[23*npt+i];
        k_f = k_f_s[5*npt+i];
        redP = alpha / k_f * phase_units[5] * low_A[5] * exp(low_beta[5] * tc[i] - activation_units[5] * low_Ea[5] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[5]) > 1.e-100 ? (1.-troe_a[5])*exp(-T[i]/troe_Tsss[5]) : 0.) 
            + (fabs(troe_Ts[5]) > 1.e-100 ? troe_a[5] * exp(-T[i]/troe_Ts[5]) : 0.) 
            + (troe_len[5] == 4 ? exp(-troe_Tss[5] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[19*npt+i];
        Kc = Kc_s[5*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[18*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 7: H + C2H4 (+M) <=> C2H5 (+M) */
        phi_f = sc[1*npt+i]*sc[19*npt+i];
        alpha = mixture[i] + (TB[6][0] - 1)*sc[0*npt+i] + (TB[6][1] - 1)*sc[5*npt+i] + (TB[6][2] - 1)*sc[11*npt+i] + (TB[6][3] - 1)*sc[12*npt+i] + (TB[6][4] - 1)*sc[13*npt+i] + (TB[6][5] - 1)*sc[21*npt+i] + (TB[6][6] - 1)*sc[23*npt+i];
        k_f = k_f_s[6*npt+i];
        redP = alpha / k_f * phase_units[6] * low_A[6] * exp(low_beta[6] * tc[i] - activation_units[6] * low_Ea[6] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[6]) > 1.e-100 ? (1.-troe_a[6])*exp(-T[i]/troe_Tsss[6]) : 0.) 
            + (fabs(troe_Ts[6]) > 1.e-100 ? troe_a[6] * exp(-T[i]/troe_Ts[6]) : 0.) 
            + (troe_len[6] == 4 ? exp(-troe_Tss[6] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[20*npt+i];
        Kc = Kc_s[6*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 8: H + C2H5 (+M) <=> C2H6 (+M) */
        phi_f = sc[1*npt+i]*sc[20*npt+i];
        alpha = mixture[i] + (TB[7][0] - 1)*sc[0*npt+i] + (TB[7][1] - 1)*sc[5*npt+i] + (TB[7][2] - 1)*sc[11*npt+i] + (TB[7][3] - 1)*sc[12*npt+i] + (TB[7][4] - 1)*sc[13*npt+i] + (TB[7][5] - 1)*sc[21*npt+i] + (TB[7][6] - 1)*sc[23*npt+i];
        k_f = k_f_s[7*npt+i];
        redP = alpha / k_f * phase_units[7] * low_A[7] * exp(low_beta[7] * tc[i] - activation_units[7] * low_Ea[7] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[7]) > 1.e-100 ? (1.-troe_a[7])*exp(-T[i]/troe_Tsss[7]) : 0.) 
            + (fabs(troe_Ts[7]) > 1.e-100 ? troe_a[7] * exp(-T[i]/troe_Ts[7]) : 0.) 
            + (troe_len[7] == 4 ? exp(-troe_Tss[7] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[21*npt+i];
        Kc = Kc_s[7*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[20*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;

        /*reaction 9: H2 + CO (+M) <=> CH2O (+M) */
        phi_f = sc[0*npt+i]*sc[12*npt+i];
        alpha = mixture[i] + (TB[8][0] - 1)*sc[0*npt+i] + (TB[8][1] - 1)*sc[5*npt+i] + (TB[8][2] - 1)*sc[11*npt+i] + (TB[8][3] - 1)*sc[12*npt+i] + (TB[8][4] - 1)*sc[13*npt+i] + (TB[8][5] - 1)*sc[21*npt+i] + (TB[8][6] - 1)*sc[23*npt+i];
        k_f = k_f_s[8*npt+i];
        redP = alpha / k_f * phase_units[8] * low_A[8] * exp(low_beta[8] * tc[i] - activation_units[8] * low_Ea[8] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[8]) > 1.e-100 ? (1.-troe_a[8])*exp(-T[i]/troe_Tsss[8]) : 0.) 
            + (fabs(troe_Ts[8]) > 1.e-100 ? troe_a[8] * exp(-T[i]/troe_Ts[8]) : 0.) 
            + (troe_len[8] == 4 ? exp(-troe_Tss[8] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[15*npt+i];
        Kc = Kc_s[8*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 10: 2 OH (+M) <=> H2O2 (+M) */
        phi_f = sc[4*npt+i]*sc[4*npt+i];
        alpha = mixture[i] + (TB[9][0] - 1)*sc[0*npt+i] + (TB[9][1] - 1)*sc[5*npt+i] + (TB[9][2] - 1)*sc[11*npt+i] + (TB[9][3] - 1)*sc[12*npt+i] + (TB[9][4] - 1)*sc[13*npt+i] + (TB[9][5] - 1)*sc[21*npt+i] + (TB[9][6] - 1)*sc[23*npt+i];
        k_f = k_f_s[9*npt+i];
        redP = alpha / k_f * phase_units[9] * low_A[9] * exp(low_beta[9] * tc[i] - activation_units[9] * low_Ea[9] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[9]) > 1.e-100 ? (1.-troe_a[9])*exp(-T[i]/troe_Tsss[9]) : 0.) 
            + (fabs(troe_Ts[9]) > 1.e-100 ? troe_a[9] * exp(-T[i]/troe_Ts[9]) : 0.) 
            + (troe_len[9] == 4 ? exp(-troe_Tss[9] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i];
        Kc = Kc_s[9*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= 2 * qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 11: 2 CH3 (+M) <=> C2H6 (+M) */
        phi_f = sc[10*npt+i]*sc[10*npt+i];
        alpha = mixture[i] + (TB[10][0] - 1)*sc[0*npt+i] + (TB[10][1] - 1)*sc[5*npt+i] + (TB[10][2] - 1)*sc[11*npt+i] + (TB[10][3] - 1)*sc[12*npt+i] + (TB[10][4] - 1)*sc[13*npt+i] + (TB[10][5] - 1)*sc[21*npt+i] + (TB[10][6] - 1)*sc[23*npt+i];
        k_f = k_f_s[10*npt+i];
        redP = alpha / k_f * phase_units[10] * low_A[10] * exp(low_beta[10] * tc[i] - activation_units[10] * low_Ea[10] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[10]) > 1.e-100 ? (1.-troe_a[10])*exp(-T[i]/troe_Tsss[10]) : 0.) 
            + (fabs(troe_Ts[10]) > 1.e-100 ? troe_a[10] * exp(-T[i]/troe_Ts[10]) : 0.) 
            + (troe_len[10] == 4 ? exp(-troe_Tss[10] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[21*npt+i];
        Kc = Kc_s[10*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= 2 * qdot;
        wdot[21*npt+i] += qdot;

        /*reaction 12: C2H4 (+M) <=> H2 + C2H2 (+M) */
        phi_f = sc[19*npt+i];
        alpha = mixture[i] + (TB[11][0] - 1)*sc[0*npt+i] + (TB[11][1] - 1)*sc[5*npt+i] + (TB[11][2] - 1)*sc[11*npt+i] + (TB[11][3] - 1)*sc[12*npt+i] + (TB[11][4] - 1)*sc[13*npt+i] + (TB[11][5] - 1)*sc[21*npt+i] + (TB[11][6] - 1)*sc[23*npt+i];
        k_f = k_f_s[11*npt+i];
        redP = alpha / k_f * phase_units[11] * low_A[11] * exp(low_beta[11] * tc[i] - activation_units[11] * low_Ea[11] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[11]) > 1.e-100 ? (1.-troe_a[11])*exp(-T[i]/troe_Tsss[11]) : 0.) 
            + (fabs(troe_Ts[11]) > 1.e-100 ? troe_a[11] * exp(-T[i]/troe_Ts[11]) : 0.) 
            + (troe_len[11] == 4 ? exp(-troe_Tss[11] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[17*npt+i];
        Kc = Kc_s[11*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 13: O + H + M <=> OH + M */
        phi_f = sc[1*npt+i]*sc[2*npt+i];
        alpha = mixture[i] + (TB[12][0] - 1)*sc[0*npt+i] + (TB[12][1] - 1)*sc[5*npt+i] + (TB[12][2] - 1)*sc[11*npt+i] + (TB[12][3] - 1)*sc[12*npt+i] + (TB[12][4] - 1)*sc[13*npt+i] + (TB[12][5] - 1)*sc[21*npt+i] + (TB[12][6] - 1)*sc[23*npt+i];
        k_f = alpha * k_f_s[12*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i];
        Kc = Kc_s[12*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 14: O + CO + M <=> CO2 + M */
        phi_f = sc[2*npt+i]*sc[12*npt+i];
        alpha = mixture[i] + (TB[13][0] - 1)*sc[0*npt+i] + (TB[13][1] - 1)*sc[3*npt+i] + (TB[13][2] - 1)*sc[5*npt+i] + (TB[13][3] - 1)*sc[11*npt+i] + (TB[13][4] - 1)*sc[12*npt+i] + (TB[13][5] - 1)*sc[13*npt+i] + (TB[13][6] - 1)*sc[21*npt+i] + (TB[13][7] - 1)*sc[23*npt+i];
        k_f = alpha * k_f_s[13*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i];
        Kc = Kc_s[13*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 15: H + O2 + M <=> HO2 + M */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        alpha = mixture[i] + (TB[14][0] - 1)*sc[3*npt+i] + (TB[14][1] - 1)*sc[5*npt+i] + (TB[14][2] - 1)*sc[12*npt+i] + (TB[14][3] - 1)*sc[13*npt+i] + (TB[14][4] - 1)*sc[21*npt+i] + (TB[14][5] - 1)*sc[22*npt+i] + (TB[14][6] - 1)*sc[23*npt+i];
        k_f = alpha * k_f_s[14*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i];
        Kc = Kc_s[14*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 16: 2 H + M <=> H2 + M */
        phi_f = sc[1*npt+i]*sc[1*npt+i];
        alpha = mixture[i] + (TB[15][0] - 1)*sc[0*npt+i] + (TB[15][1] - 1)*sc[5*npt+i] + (TB[15][2] - 1)*sc[11*npt+i] + (TB[15][3] - 1)*sc[13*npt+i] + (TB[15][4] - 1)*sc[21*npt+i] + (TB[15][5] - 1)*sc[23*npt+i];
        k_f = alpha * k_f_s[15*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i];
        Kc = Kc_s[15*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= 2 * qdot;

        /*reaction 17: H + OH + M <=> H2O + M */
        phi_f = sc[1*npt+i]*sc[4*npt+i];
        alpha = mixture[i] + (TB[16][0] - 1)*sc[0*npt+i] + (TB[16][1] - 1)*sc[5*npt+i] + (TB[16][2] - 1)*sc[11*npt+i] + (TB[16][3] - 1)*sc[21*npt+i] + (TB[16][4] - 1)*sc[23*npt+i];
        k_f = alpha * k_f_s[16*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i];
        Kc = Kc_s[16*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 18: HCO + M <=> H + CO + M */
        phi_f = sc[14*npt+i];
        alpha = mixture[i] + (TB[17][0] - 1)*sc[0*npt+i] + (TB[17][1] - 1)*sc[5*npt+i] + (TB[17][2] - 1)*sc[11*npt+i] + (TB[17][3] - 1)*sc[12*npt+i] + (TB[17][4] - 1)*sc[13*npt+i] + (TB[17][5] - 1)*sc[21*npt+i];
        k_f = alpha * k_f_s[17*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[12*npt+i];
        Kc = Kc_s[17*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 19: O + H2 <=> H + OH */
        phi_f = sc[0*npt+i]*sc[2*npt+i];
        k_f = k_f_s[18*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[4*npt+i];
        Kc = Kc_s[18*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 20: O + HO2 <=> OH + O2 */
        phi_f = sc[2*npt+i]*sc[6*npt+i];
        k_f = k_f_s[19*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[4*npt+i];
        Kc = Kc_s[19*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 21: O + CH2 <=> H + HCO */
        phi_f = sc[2*npt+i]*sc[8*npt+i];
        k_f = k_f_s[20*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[14*npt+i];
        Kc = Kc_s[20*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 22: O + CH2(S) <=> H + HCO */
        phi_f = sc[2*npt+i]*sc[9*npt+i];
        k_f = k_f_s[21*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[14*npt+i];
        Kc = Kc_s[21*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 23: O + CH3 <=> H + CH2O */
        phi_f = sc[2*npt+i]*sc[10*npt+i];
        k_f = k_f_s[22*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[15*npt+i];
        Kc = Kc_s[22*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 24: O + CH4 <=> OH + CH3 */
        phi_f = sc[2*npt+i]*sc[11*npt+i];
        k_f = k_f_s[23*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[10*npt+i];
        Kc = Kc_s[23*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 25: O + HCO <=> OH + CO */
        phi_f = sc[2*npt+i]*sc[14*npt+i];
        k_f = k_f_s[24*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[12*npt+i];
        Kc = Kc_s[24*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 26: O + HCO <=> H + CO2 */
        phi_f = sc[2*npt+i]*sc[14*npt+i];
        k_f = k_f_s[25*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[13*npt+i];
        Kc = Kc_s[25*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 27: O + CH2O <=> OH + HCO */
        phi_f = sc[2*npt+i]*sc[15*npt+i];
        k_f = k_f_s[26*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[14*npt+i];
        Kc = Kc_s[26*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 28: O + C2H2 <=> CH2(S) + CO */
        phi_f = sc[2*npt+i]*sc[17*npt+i];
        k_f = k_f_s[27*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[12*npt+i];
        Kc = Kc_s[27*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 29: O + C2H2 <=> CO + CH2 */
        phi_f = sc[2*npt+i]*sc[17*npt+i];
        k_f = k_f_s[28*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[12*npt+i];
        Kc = Kc_s[28*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 30: O + C2H4 <=> CH3 + HCO */
        phi_f = sc[2*npt+i]*sc[19*npt+i];
        k_f = k_f_s[29*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[14*npt+i];
        Kc = Kc_s[29*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 31: O + C2H5 <=> CH3 + CH2O */
        phi_f = sc[2*npt+i]*sc[20*npt+i];
        k_f = k_f_s[30*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[15*npt+i];
        Kc = Kc_s[30*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 32: O + C2H6 <=> OH + C2H5 */
        phi_f = sc[2*npt+i]*sc[21*npt+i];
        k_f = k_f_s[31*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[20*npt+i];
        Kc = Kc_s[31*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[20*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 33: O2 + CO <=> O + CO2 */
        phi_f = sc[3*npt+i]*sc[12*npt+i];
        k_f = k_f_s[32*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[13*npt+i];
        Kc = Kc_s[32*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 34: O2 + CH2O <=> HO2 + HCO */
        phi_f = sc[3*npt+i]*sc[15*npt+i];
        k_f = k_f_s[33*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[14*npt+i];
        Kc = Kc_s[33*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 35: H + 2 O2 <=> HO2 + O2 */
        phi_f = sc[1*npt+i]*sc[3*npt+i]*sc[3*npt+i];
        k_f = k_f_s[34*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[6*npt+i];
        Kc = Kc_s[34*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= 2 * qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 36: H + O2 + H2O <=> HO2 + H2O */
        phi_f = sc[1*npt+i]*sc[3*npt+i]*sc[5*npt+i];
        k_f = k_f_s[35*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[35*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 37: H + O2 + N2 <=> HO2 + N2 */
        phi_f = sc[1*npt+i]*sc[3*npt+i]*sc[22*npt+i];
        k_f = k_f_s[36*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[22*npt+i];
        Kc = Kc_s[36*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;
        wdot[22*npt+i] += qdot;

        /*reaction 38: H + O2 + AR <=> HO2 + AR */
        phi_f = sc[1*npt+i]*sc[3*npt+i]*sc[23*npt+i];
        k_f = k_f_s[37*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[23*npt+i];
        Kc = Kc_s[37*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[23*npt+i] -= qdot;
        wdot[23*npt+i] += qdot;

        /*reaction 39: H + O2 <=> O + OH */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        k_f = k_f_s[38*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[4*npt+i];
        Kc = Kc_s[38*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 40: 2 H + H2 <=> 2 H2 */
        phi_f = sc[0*npt+i]*sc[1*npt+i]*sc[1*npt+i];
        k_f = k_f_s[39*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[0*npt+i];
        Kc = Kc_s[39*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[0*npt+i] += 2 * qdot;
        wdot[1*npt+i] -= 2 * qdot;

        /*reaction 41: 2 H + H2O <=> H2 + H2O */
        phi_f = sc[1*npt+i]*sc[1*npt+i]*sc[5*npt+i];
        k_f = k_f_s[40*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[5*npt+i];
        Kc = Kc_s[40*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= 2 * qdot;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
        phi_f = sc[1*npt+i]*sc[1*npt+i]*sc[13*npt+i];
        k_f = k_f_s[41*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[13*npt+i];
        Kc = Kc_s[41*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= 2 * qdot;
        wdot[13*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 43: H + HO2 <=> O2 + H2 */
        phi_f = sc[1*npt+i]*sc[6*npt+i];
        k_f = k_f_s[42*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[3*npt+i];
        Kc = Kc_s[42*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 44: H + HO2 <=> 2 OH */
        phi_f = sc[1*npt+i]*sc[6*npt+i];
        k_f = k_f_s[43*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[4*npt+i];
        Kc = Kc_s[43*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += 2 * qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 45: H + H2O2 <=> HO2 + H2 */
        phi_f = sc[1*npt+i]*sc[7*npt+i];
        k_f = k_f_s[44*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[6*npt+i];
        Kc = Kc_s[44*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 46: H + CH4 <=> CH3 + H2 */
        phi_f = sc[1*npt+i]*sc[11*npt+i];
        k_f = k_f_s[45*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[10*npt+i];
        Kc = Kc_s[45*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 47: H + HCO <=> H2 + CO */
        phi_f = sc[1*npt+i]*sc[14*npt+i];
        k_f = k_f_s[46*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[12*npt+i];
        Kc = Kc_s[46*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 48: H + CH2O <=> HCO + H2 */
        phi_f = sc[1*npt+i]*sc[15*npt+i];
        k_f = k_f_s[47*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[14*npt+i];
        Kc = Kc_s[47*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 49: H + CH3O <=> OH + CH3 */
        phi_f = sc[1*npt+i]*sc[16*npt+i];
        k_f = k_f_s[48*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[10*npt+i];
        Kc = Kc_s[48*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 50: H + C2H3 <=> H2 + C2H2 */
        phi_f = sc[1*npt+i]*sc[18*npt+i];
        k_f = k_f_s[49*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[17*npt+i];
        Kc = Kc_s[49*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;
    }
}

void vcomp_wdot_51_100(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,
		double * restrict k_f_s, double * restrict Kc_s,
		double * restrict tc, double * restrict invT, double * restrict T)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;

        /*reaction 51: H + C2H4 <=> C2H3 + H2 */
        phi_f = sc[1*npt+i]*sc[19*npt+i];
        k_f = k_f_s[50*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[18*npt+i];
        Kc = Kc_s[50*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[18*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 52: H + C2H6 <=> C2H5 + H2 */
        phi_f = sc[1*npt+i]*sc[21*npt+i];
        k_f = k_f_s[51*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[20*npt+i];
        Kc = Kc_s[51*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 53: OH + H2 <=> H + H2O */
        phi_f = sc[0*npt+i]*sc[4*npt+i];
        k_f = k_f_s[52*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[5*npt+i];
        Kc = Kc_s[52*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 54: 2 OH <=> O + H2O */
        phi_f = sc[4*npt+i]*sc[4*npt+i];
        k_f = k_f_s[53*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[5*npt+i];
        Kc = Kc_s[53*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= 2 * qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 55: OH + HO2 <=> O2 + H2O */
        phi_f = sc[4*npt+i]*sc[6*npt+i];
        k_f = k_f_s[54*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[5*npt+i];
        Kc = Kc_s[54*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 56: OH + H2O2 <=> HO2 + H2O */
        phi_f = sc[4*npt+i]*sc[7*npt+i];
        k_f = k_f_s[55*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[55*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 57: OH + CH2 <=> H + CH2O */
        phi_f = sc[4*npt+i]*sc[8*npt+i];
        k_f = k_f_s[56*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[15*npt+i];
        Kc = Kc_s[56*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 58: OH + CH2(S) <=> H + CH2O */
        phi_f = sc[4*npt+i]*sc[9*npt+i];
        k_f = k_f_s[57*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[15*npt+i];
        Kc = Kc_s[57*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 59: OH + CH3 <=> CH2 + H2O */
        phi_f = sc[4*npt+i]*sc[10*npt+i];
        k_f = k_f_s[58*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[8*npt+i];
        Kc = Kc_s[58*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 60: OH + CH3 <=> CH2(S) + H2O */
        phi_f = sc[4*npt+i]*sc[10*npt+i];
        k_f = k_f_s[59*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[9*npt+i];
        Kc = Kc_s[59*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 61: OH + CH4 <=> CH3 + H2O */
        phi_f = sc[4*npt+i]*sc[11*npt+i];
        k_f = k_f_s[60*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[10*npt+i];
        Kc = Kc_s[60*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 62: OH + CO <=> H + CO2 */
        phi_f = sc[4*npt+i]*sc[12*npt+i];
        k_f = k_f_s[61*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[13*npt+i];
        Kc = Kc_s[61*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 63: OH + HCO <=> H2O + CO */
        phi_f = sc[4*npt+i]*sc[14*npt+i];
        k_f = k_f_s[62*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[12*npt+i];
        Kc = Kc_s[62*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 64: OH + CH2O <=> HCO + H2O */
        phi_f = sc[4*npt+i]*sc[15*npt+i];
        k_f = k_f_s[63*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[14*npt+i];
        Kc = Kc_s[63*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 65: OH + C2H2 <=> CH3 + CO */
        phi_f = sc[4*npt+i]*sc[17*npt+i];
        k_f = k_f_s[64*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[12*npt+i];
        Kc = Kc_s[64*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 66: OH + C2H3 <=> H2O + C2H2 */
        phi_f = sc[4*npt+i]*sc[18*npt+i];
        k_f = k_f_s[65*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[17*npt+i];
        Kc = Kc_s[65*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 67: OH + C2H4 <=> C2H3 + H2O */
        phi_f = sc[4*npt+i]*sc[19*npt+i];
        k_f = k_f_s[66*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[18*npt+i];
        Kc = Kc_s[66*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 68: OH + C2H6 <=> C2H5 + H2O */
        phi_f = sc[4*npt+i]*sc[21*npt+i];
        k_f = k_f_s[67*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[20*npt+i];
        Kc = Kc_s[67*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[20*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 69: 2 HO2 <=> O2 + H2O2 */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[68*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[7*npt+i];
        Kc = Kc_s[68*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= 2 * qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 70: 2 HO2 <=> O2 + H2O2 */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[69*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[7*npt+i];
        Kc = Kc_s[69*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= 2 * qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 71: HO2 + CH2 <=> OH + CH2O */
        phi_f = sc[6*npt+i]*sc[8*npt+i];
        k_f = k_f_s[70*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[15*npt+i];
        Kc = Kc_s[70*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 72: HO2 + CH3 <=> O2 + CH4 */
        phi_f = sc[6*npt+i]*sc[10*npt+i];
        k_f = k_f_s[71*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[11*npt+i];
        Kc = Kc_s[71*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 73: HO2 + CH3 <=> OH + CH3O */
        phi_f = sc[6*npt+i]*sc[10*npt+i];
        k_f = k_f_s[72*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[16*npt+i];
        Kc = Kc_s[72*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 74: HO2 + CO <=> OH + CO2 */
        phi_f = sc[6*npt+i]*sc[12*npt+i];
        k_f = k_f_s[73*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[13*npt+i];
        Kc = Kc_s[73*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 75: HO2 + CH2O <=> HCO + H2O2 */
        phi_f = sc[6*npt+i]*sc[15*npt+i];
        k_f = k_f_s[74*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[14*npt+i];
        Kc = Kc_s[74*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 76: CH2 + O2 <=> OH + HCO */
        phi_f = sc[3*npt+i]*sc[8*npt+i];
        k_f = k_f_s[75*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[14*npt+i];
        Kc = Kc_s[75*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 77: CH2 + H2 <=> H + CH3 */
        phi_f = sc[0*npt+i]*sc[8*npt+i];
        k_f = k_f_s[76*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[10*npt+i];
        Kc = Kc_s[76*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 78: 2 CH2 <=> H2 + C2H2 */
        phi_f = sc[8*npt+i]*sc[8*npt+i];
        k_f = k_f_s[77*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[17*npt+i];
        Kc = Kc_s[77*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[8*npt+i] -= 2 * qdot;
        wdot[17*npt+i] += qdot;

        /*reaction 79: CH2 + CH3 <=> H + C2H4 */
        phi_f = sc[8*npt+i]*sc[10*npt+i];
        k_f = k_f_s[78*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[19*npt+i];
        Kc = Kc_s[78*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 80: CH2 + CH4 <=> 2 CH3 */
        phi_f = sc[8*npt+i]*sc[11*npt+i];
        k_f = k_f_s[79*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[10*npt+i];
        Kc = Kc_s[79*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[8*npt+i] -= qdot;
        wdot[10*npt+i] += 2 * qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 81: CH2(S) + N2 <=> CH2 + N2 */
        phi_f = sc[9*npt+i]*sc[22*npt+i];
        k_f = k_f_s[80*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[22*npt+i];
        Kc = Kc_s[80*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[22*npt+i] -= qdot;
        wdot[22*npt+i] += qdot;

        /*reaction 82: CH2(S) + AR <=> CH2 + AR */
        phi_f = sc[9*npt+i]*sc[23*npt+i];
        k_f = k_f_s[81*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[23*npt+i];
        Kc = Kc_s[81*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[23*npt+i] -= qdot;
        wdot[23*npt+i] += qdot;

        /*reaction 83: CH2(S) + O2 <=> H + OH + CO */
        phi_f = sc[3*npt+i]*sc[9*npt+i];
        k_f = k_f_s[82*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[4*npt+i]*sc[12*npt+i];
        Kc = Kc_s[82*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 84: CH2(S) + O2 <=> CO + H2O */
        phi_f = sc[3*npt+i]*sc[9*npt+i];
        k_f = k_f_s[83*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[12*npt+i];
        Kc = Kc_s[83*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 85: CH2(S) + H2 <=> CH3 + H */
        phi_f = sc[0*npt+i]*sc[9*npt+i];
        k_f = k_f_s[84*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[10*npt+i];
        Kc = Kc_s[84*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 86: CH2(S) + H2O <=> CH2 + H2O */
        phi_f = sc[5*npt+i]*sc[9*npt+i];
        k_f = k_f_s[85*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[8*npt+i];
        Kc = Kc_s[85*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;

        /*reaction 87: CH2(S) + CH3 <=> H + C2H4 */
        phi_f = sc[9*npt+i]*sc[10*npt+i];
        k_f = k_f_s[86*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[19*npt+i];
        Kc = Kc_s[86*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 88: CH2(S) + CH4 <=> 2 CH3 */
        phi_f = sc[9*npt+i]*sc[11*npt+i];
        k_f = k_f_s[87*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[10*npt+i];
        Kc = Kc_s[87*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] -= qdot;
        wdot[10*npt+i] += 2 * qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 89: CH2(S) + CO <=> CH2 + CO */
        phi_f = sc[9*npt+i]*sc[12*npt+i];
        k_f = k_f_s[88*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[12*npt+i];
        Kc = Kc_s[88*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 90: CH2(S) + CO2 <=> CH2 + CO2 */
        phi_f = sc[9*npt+i]*sc[13*npt+i];
        k_f = k_f_s[89*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[13*npt+i];
        Kc = Kc_s[89*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[13*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 91: CH2(S) + CO2 <=> CO + CH2O */
        phi_f = sc[9*npt+i]*sc[13*npt+i];
        k_f = k_f_s[90*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[12*npt+i]*sc[15*npt+i];
        Kc = Kc_s[90*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 92: CH3 + O2 <=> O + CH3O */
        phi_f = sc[3*npt+i]*sc[10*npt+i];
        k_f = k_f_s[91*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[16*npt+i];
        Kc = Kc_s[91*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 93: CH3 + O2 <=> OH + CH2O */
        phi_f = sc[3*npt+i]*sc[10*npt+i];
        k_f = k_f_s[92*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[15*npt+i];
        Kc = Kc_s[92*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 94: CH3 + H2O2 <=> HO2 + CH4 */
        phi_f = sc[7*npt+i]*sc[10*npt+i];
        k_f = k_f_s[93*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[11*npt+i];
        Kc = Kc_s[93*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 95: 2 CH3 <=> H + C2H5 */
        phi_f = sc[10*npt+i]*sc[10*npt+i];
        k_f = k_f_s[94*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[20*npt+i];
        Kc = Kc_s[94*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[10*npt+i] -= 2 * qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 96: CH3 + HCO <=> CH4 + CO */
        phi_f = sc[10*npt+i]*sc[14*npt+i];
        k_f = k_f_s[95*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[12*npt+i];
        Kc = Kc_s[95*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 97: CH3 + CH2O <=> HCO + CH4 */
        phi_f = sc[10*npt+i]*sc[15*npt+i];
        k_f = k_f_s[96*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[14*npt+i];
        Kc = Kc_s[96*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 98: CH3 + C2H4 <=> C2H3 + CH4 */
        phi_f = sc[10*npt+i]*sc[19*npt+i];
        k_f = k_f_s[97*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[18*npt+i];
        Kc = Kc_s[97*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 99: CH3 + C2H6 <=> C2H5 + CH4 */
        phi_f = sc[10*npt+i]*sc[21*npt+i];
        k_f = k_f_s[98*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[20*npt+i];
        Kc = Kc_s[98*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[20*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 100: HCO + H2O <=> H + CO + H2O */
        phi_f = sc[5*npt+i]*sc[14*npt+i];
        k_f = k_f_s[99*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[5*npt+i]*sc[12*npt+i];
        Kc = Kc_s[99*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;
    }
}

void vcomp_wdot_101_104(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,
		double * restrict k_f_s, double * restrict Kc_s,
		double * restrict tc, double * restrict invT, double * restrict T)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;

        /*reaction 101: HCO + O2 <=> HO2 + CO */
        phi_f = sc[3*npt+i]*sc[14*npt+i];
        k_f = k_f_s[100*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[12*npt+i];
        Kc = Kc_s[100*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 102: CH3O + O2 <=> HO2 + CH2O */
        phi_f = sc[3*npt+i]*sc[16*npt+i];
        k_f = k_f_s[101*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[15*npt+i];
        Kc = Kc_s[101*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 103: C2H3 + O2 <=> HCO + CH2O */
        phi_f = sc[3*npt+i]*sc[18*npt+i];
        k_f = k_f_s[102*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[14*npt+i]*sc[15*npt+i];
        Kc = Kc_s[102*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 104: C2H5 + O2 <=> HO2 + C2H4 */
        phi_f = sc[3*npt+i]*sc[20*npt+i];
        k_f = k_f_s[103*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[19*npt+i];
        Kc = Kc_s[103*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[19*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;
    }
}


/*compute the reaction Jacobian */
void DWDOT(double * restrict J, double * restrict sc, double * restrict Tp, int * consP)
{
    double c[24];

    for (int k=0; k<24; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    for (int k=0; k<24; k++) {
        J[600+k] *= 1.e-6;
    }

    /* dTdot/d[X] */
    for (int k=0; k<24; k++) {
        J[k*25+24] *= 1.e6;
    }

    return;
}

/*compute the reaction Jacobian */
void aJacobian(double * restrict J, double * restrict sc, double T, int consP)
{
    for (int i=0; i<625; i++) {
        J[i] = 0.0;
    }

    double wdot[24];
    for (int k=0; k<24; k++) {
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
    for (int k = 0; k < 24; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[24];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[24];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[24];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[5] + (TB[0][2] - 1)*sc[11] + (TB[0][3] - 1)*sc[12] + (TB[0][4] - 1)*sc[13] + (TB[0][5] - 1)*sc[21] + (TB[0][6] - 1)*sc[23];
    /* forward */
    phi_f = sc[1]*sc[8];
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
    phi_r = sc[10];
    Kc = refCinv * exp(g_RT[1] + g_RT[8] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[8]) + (h_RT[10]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[8] -= q; /* CH2 */
    wdot[10] += q; /* CH3 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[0][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[8] -= dqdci;                /* dwdot[CH2]/d[H2] */
        J[10] += dqdci;               /* dwdot[CH3]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[8];
        J[26] -= dqdci;               /* dwdot[H]/d[H] */
        J[33] -= dqdci;               /* dwdot[CH2]/d[H] */
        J[35] += dqdci;               /* dwdot[CH3]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[0][1] - 1)*dcdc_fac;
        J[126] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[133] -= dqdci;              /* dwdot[CH2]/d[H2O] */
        J[135] += dqdci;              /* dwdot[CH3]/d[H2O] */
        /* d()/d[CH2] */
        dqdci =  + k_f*sc[1];
        J[201] -= dqdci;              /* dwdot[H]/d[CH2] */
        J[208] -= dqdci;              /* dwdot[CH2]/d[CH2] */
        J[210] += dqdci;              /* dwdot[CH3]/d[CH2] */
        /* d()/d[CH3] */
        dqdci =  - k_r;
        J[251] -= dqdci;              /* dwdot[H]/d[CH3] */
        J[258] -= dqdci;              /* dwdot[CH2]/d[CH3] */
        J[260] += dqdci;              /* dwdot[CH3]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (TB[0][2] - 1)*dcdc_fac;
        J[276] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[283] -= dqdci;              /* dwdot[CH2]/d[CH4] */
        J[285] += dqdci;              /* dwdot[CH3]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[0][3] - 1)*dcdc_fac;
        J[301] -= dqdci;              /* dwdot[H]/d[CO] */
        J[308] -= dqdci;              /* dwdot[CH2]/d[CO] */
        J[310] += dqdci;              /* dwdot[CH3]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[0][4] - 1)*dcdc_fac;
        J[326] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[333] -= dqdci;              /* dwdot[CH2]/d[CO2] */
        J[335] += dqdci;              /* dwdot[CH3]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[0][5] - 1)*dcdc_fac;
        J[526] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[533] -= dqdci;              /* dwdot[CH2]/d[C2H6] */
        J[535] += dqdci;              /* dwdot[CH3]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[0][6] - 1)*dcdc_fac;
        J[576] -= dqdci;              /* dwdot[H]/d[AR] */
        J[583] -= dqdci;              /* dwdot[CH2]/d[AR] */
        J[585] += dqdci;              /* dwdot[CH3]/d[AR] */
    }
    else {
        dqdc[0] = TB[0][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[8];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[0][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac + k_f*sc[1];
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac - k_r;
        dqdc[11] = TB[0][2]*dcdc_fac;
        dqdc[12] = TB[0][3]*dcdc_fac;
        dqdc[13] = TB[0][4]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = TB[0][5]*dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = TB[0][6]*dcdc_fac;
        for (int k=0; k<24; k++) {
            J[25*k+1] -= dqdc[k];
            J[25*k+8] -= dqdc[k];
            J[25*k+10] += dqdc[k];
        }
    }
    J[601] -= dqdT; /* dwdot[H]/dT */
    J[608] -= dqdT; /* dwdot[CH2]/dT */
    J[610] += dqdT; /* dwdot[CH3]/dT */

    /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[5] + (TB[1][2] - 1)*sc[11] + (TB[1][3] - 1)*sc[12] + (TB[1][4] - 1)*sc[13] + (TB[1][5] - 1)*sc[21] + (TB[1][6] - 1)*sc[23];
    /* forward */
    phi_f = sc[1]*sc[10];
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
    phi_r = sc[11];
    Kc = refCinv * exp(g_RT[1] + g_RT[10] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (h_RT[11]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[10] -= q; /* CH3 */
    wdot[11] += q; /* CH4 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[1][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[10] -= dqdci;               /* dwdot[CH3]/d[H2] */
        J[11] += dqdci;               /* dwdot[CH4]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[10];
        J[26] -= dqdci;               /* dwdot[H]/d[H] */
        J[35] -= dqdci;               /* dwdot[CH3]/d[H] */
        J[36] += dqdci;               /* dwdot[CH4]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[1][1] - 1)*dcdc_fac;
        J[126] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[135] -= dqdci;              /* dwdot[CH3]/d[H2O] */
        J[136] += dqdci;              /* dwdot[CH4]/d[H2O] */
        /* d()/d[CH3] */
        dqdci =  + k_f*sc[1];
        J[251] -= dqdci;              /* dwdot[H]/d[CH3] */
        J[260] -= dqdci;              /* dwdot[CH3]/d[CH3] */
        J[261] += dqdci;              /* dwdot[CH4]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (TB[1][2] - 1)*dcdc_fac - k_r;
        J[276] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[285] -= dqdci;              /* dwdot[CH3]/d[CH4] */
        J[286] += dqdci;              /* dwdot[CH4]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[1][3] - 1)*dcdc_fac;
        J[301] -= dqdci;              /* dwdot[H]/d[CO] */
        J[310] -= dqdci;              /* dwdot[CH3]/d[CO] */
        J[311] += dqdci;              /* dwdot[CH4]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[1][4] - 1)*dcdc_fac;
        J[326] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[335] -= dqdci;              /* dwdot[CH3]/d[CO2] */
        J[336] += dqdci;              /* dwdot[CH4]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[1][5] - 1)*dcdc_fac;
        J[526] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[535] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
        J[536] += dqdci;              /* dwdot[CH4]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[1][6] - 1)*dcdc_fac;
        J[576] -= dqdci;              /* dwdot[H]/d[AR] */
        J[585] -= dqdci;              /* dwdot[CH3]/d[AR] */
        J[586] += dqdci;              /* dwdot[CH4]/d[AR] */
    }
    else {
        dqdc[0] = TB[1][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[10];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[1][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac + k_f*sc[1];
        dqdc[11] = TB[1][2]*dcdc_fac - k_r;
        dqdc[12] = TB[1][3]*dcdc_fac;
        dqdc[13] = TB[1][4]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = TB[1][5]*dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = TB[1][6]*dcdc_fac;
        for (int k=0; k<24; k++) {
            J[25*k+1] -= dqdc[k];
            J[25*k+10] -= dqdc[k];
            J[25*k+11] += dqdc[k];
        }
    }
    J[601] -= dqdT; /* dwdot[H]/dT */
    J[610] -= dqdT; /* dwdot[CH3]/dT */
    J[611] += dqdT; /* dwdot[CH4]/dT */

    /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[2][0] - 1)*sc[0] + (TB[2][1] - 1)*sc[5] + (TB[2][2] - 1)*sc[11] + (TB[2][3] - 1)*sc[12] + (TB[2][4] - 1)*sc[13] + (TB[2][5] - 1)*sc[21] + (TB[2][6] - 1)*sc[23];
    /* forward */
    phi_f = sc[1]*sc[14];
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
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[2]) > 1.e-100 ? (1.-troe_a[2])*exp(-T/troe_Tsss[2]) : 0.);
    Fcent2 = (fabs(troe_Ts[2]) > 1.e-100 ? troe_a[2] * exp(-T/troe_Ts[2]) : 0.);
    Fcent3 = (troe_len[2] == 4 ? exp(-troe_Tss[2] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[2]) > 1.e-100 ? -Fcent1/troe_Tsss[2] : 0.)
      + (fabs(troe_Ts[2]) > 1.e-100 ? -Fcent2/troe_Ts[2] : 0.)
      + (troe_len[2] == 4 ? Fcent3*troe_Tss[2]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[15];
    Kc = refCinv * exp(g_RT[1] + g_RT[14] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[14]) + (h_RT[15]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[14] -= q; /* HCO */
    wdot[15] += q; /* CH2O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[2][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[14] -= dqdci;               /* dwdot[HCO]/d[H2] */
        J[15] += dqdci;               /* dwdot[CH2O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[14];
        J[26] -= dqdci;               /* dwdot[H]/d[H] */
        J[39] -= dqdci;               /* dwdot[HCO]/d[H] */
        J[40] += dqdci;               /* dwdot[CH2O]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[2][1] - 1)*dcdc_fac;
        J[126] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[139] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        J[140] += dqdci;              /* dwdot[CH2O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[2][2] - 1)*dcdc_fac;
        J[276] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[289] -= dqdci;              /* dwdot[HCO]/d[CH4] */
        J[290] += dqdci;              /* dwdot[CH2O]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[2][3] - 1)*dcdc_fac;
        J[301] -= dqdci;              /* dwdot[H]/d[CO] */
        J[314] -= dqdci;              /* dwdot[HCO]/d[CO] */
        J[315] += dqdci;              /* dwdot[CH2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[2][4] - 1)*dcdc_fac;
        J[326] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[339] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        J[340] += dqdci;              /* dwdot[CH2O]/d[CO2] */
        /* d()/d[HCO] */
        dqdci =  + k_f*sc[1];
        J[351] -= dqdci;              /* dwdot[H]/d[HCO] */
        J[364] -= dqdci;              /* dwdot[HCO]/d[HCO] */
        J[365] += dqdci;              /* dwdot[CH2O]/d[HCO] */
        /* d()/d[CH2O] */
        dqdci =  - k_r;
        J[376] -= dqdci;              /* dwdot[H]/d[CH2O] */
        J[389] -= dqdci;              /* dwdot[HCO]/d[CH2O] */
        J[390] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
        /* d()/d[C2H6] */
        dqdci = (TB[2][5] - 1)*dcdc_fac;
        J[526] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[539] -= dqdci;              /* dwdot[HCO]/d[C2H6] */
        J[540] += dqdci;              /* dwdot[CH2O]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[2][6] - 1)*dcdc_fac;
        J[576] -= dqdci;              /* dwdot[H]/d[AR] */
        J[589] -= dqdci;              /* dwdot[HCO]/d[AR] */
        J[590] += dqdci;              /* dwdot[CH2O]/d[AR] */
    }
    else {
        dqdc[0] = TB[2][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[14];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[2][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[2][2]*dcdc_fac;
        dqdc[12] = TB[2][3]*dcdc_fac;
        dqdc[13] = TB[2][4]*dcdc_fac;
        dqdc[14] = dcdc_fac + k_f*sc[1];
        dqdc[15] = dcdc_fac - k_r;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = TB[2][5]*dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = TB[2][6]*dcdc_fac;
        for (int k=0; k<24; k++) {
            J[25*k+1] -= dqdc[k];
            J[25*k+14] -= dqdc[k];
            J[25*k+15] += dqdc[k];
        }
    }
    J[601] -= dqdT; /* dwdot[H]/dT */
    J[614] -= dqdT; /* dwdot[HCO]/dT */
    J[615] += dqdT; /* dwdot[CH2O]/dT */

    /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[3][0] - 1)*sc[0] + (TB[3][1] - 1)*sc[5] + (TB[3][2] - 1)*sc[11] + (TB[3][3] - 1)*sc[12] + (TB[3][4] - 1)*sc[13] + (TB[3][5] - 1)*sc[21];
    /* forward */
    phi_f = sc[1]*sc[15];
    k_f = prefactor_units[3] * fwd_A[3]
                * exp(fwd_beta[3] * tc[0] - activation_units[3] * fwd_Ea[3] * invT);
    dlnkfdT = fwd_beta[3] * invT + activation_units[3] * fwd_Ea[3] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[3] * exp(low_beta[3] * tc[0] - activation_units[3] * low_Ea[3] * invT);
    Pr = phase_units[3] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[3] * invT + activation_units[3] * low_Ea[3] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[3]) > 1.e-100 ? (1.-troe_a[3])*exp(-T/troe_Tsss[3]) : 0.);
    Fcent2 = (fabs(troe_Ts[3]) > 1.e-100 ? troe_a[3] * exp(-T/troe_Ts[3]) : 0.);
    Fcent3 = (troe_len[3] == 4 ? exp(-troe_Tss[3] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[3]) > 1.e-100 ? -Fcent1/troe_Tsss[3] : 0.)
      + (fabs(troe_Ts[3]) > 1.e-100 ? -Fcent2/troe_Ts[3] : 0.)
      + (troe_len[3] == 4 ? Fcent3*troe_Tss[3]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[16];
    Kc = refCinv * exp(g_RT[1] + g_RT[15] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[15]) + (h_RT[16]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[15] -= q; /* CH2O */
    wdot[16] += q; /* CH3O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[3][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[15] -= dqdci;               /* dwdot[CH2O]/d[H2] */
        J[16] += dqdci;               /* dwdot[CH3O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[15];
        J[26] -= dqdci;               /* dwdot[H]/d[H] */
        J[40] -= dqdci;               /* dwdot[CH2O]/d[H] */
        J[41] += dqdci;               /* dwdot[CH3O]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[3][1] - 1)*dcdc_fac;
        J[126] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[140] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
        J[141] += dqdci;              /* dwdot[CH3O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[3][2] - 1)*dcdc_fac;
        J[276] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[290] -= dqdci;              /* dwdot[CH2O]/d[CH4] */
        J[291] += dqdci;              /* dwdot[CH3O]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[3][3] - 1)*dcdc_fac;
        J[301] -= dqdci;              /* dwdot[H]/d[CO] */
        J[315] -= dqdci;              /* dwdot[CH2O]/d[CO] */
        J[316] += dqdci;              /* dwdot[CH3O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[3][4] - 1)*dcdc_fac;
        J[326] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[340] -= dqdci;              /* dwdot[CH2O]/d[CO2] */
        J[341] += dqdci;              /* dwdot[CH3O]/d[CO2] */
        /* d()/d[CH2O] */
        dqdci =  + k_f*sc[1];
        J[376] -= dqdci;              /* dwdot[H]/d[CH2O] */
        J[390] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
        J[391] += dqdci;              /* dwdot[CH3O]/d[CH2O] */
        /* d()/d[CH3O] */
        dqdci =  - k_r;
        J[401] -= dqdci;              /* dwdot[H]/d[CH3O] */
        J[415] -= dqdci;              /* dwdot[CH2O]/d[CH3O] */
        J[416] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
        /* d()/d[C2H6] */
        dqdci = (TB[3][5] - 1)*dcdc_fac;
        J[526] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[540] -= dqdci;              /* dwdot[CH2O]/d[C2H6] */
        J[541] += dqdci;              /* dwdot[CH3O]/d[C2H6] */
    }
    else {
        dqdc[0] = TB[3][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[15];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[3][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[3][2]*dcdc_fac;
        dqdc[12] = TB[3][3]*dcdc_fac;
        dqdc[13] = TB[3][4]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac + k_f*sc[1];
        dqdc[16] = dcdc_fac - k_r;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = TB[3][5]*dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        for (int k=0; k<24; k++) {
            J[25*k+1] -= dqdc[k];
            J[25*k+15] -= dqdc[k];
            J[25*k+16] += dqdc[k];
        }
    }
    J[601] -= dqdT; /* dwdot[H]/dT */
    J[615] -= dqdT; /* dwdot[CH2O]/dT */
    J[616] += dqdT; /* dwdot[CH3O]/dT */

    /*reaction 5: H + C2H2 (+M) <=> C2H3 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[4][0] - 1)*sc[0] + (TB[4][1] - 1)*sc[5] + (TB[4][2] - 1)*sc[11] + (TB[4][3] - 1)*sc[12] + (TB[4][4] - 1)*sc[13] + (TB[4][5] - 1)*sc[21] + (TB[4][6] - 1)*sc[23];
    /* forward */
    phi_f = sc[1]*sc[17];
    k_f = prefactor_units[4] * fwd_A[4]
                * exp(fwd_beta[4] * tc[0] - activation_units[4] * fwd_Ea[4] * invT);
    dlnkfdT = fwd_beta[4] * invT + activation_units[4] * fwd_Ea[4] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[4] * exp(low_beta[4] * tc[0] - activation_units[4] * low_Ea[4] * invT);
    Pr = phase_units[4] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[4] * invT + activation_units[4] * low_Ea[4] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[4]) > 1.e-100 ? (1.-troe_a[4])*exp(-T/troe_Tsss[4]) : 0.);
    Fcent2 = (fabs(troe_Ts[4]) > 1.e-100 ? troe_a[4] * exp(-T/troe_Ts[4]) : 0.);
    Fcent3 = (troe_len[4] == 4 ? exp(-troe_Tss[4] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[4]) > 1.e-100 ? -Fcent1/troe_Tsss[4] : 0.)
      + (fabs(troe_Ts[4]) > 1.e-100 ? -Fcent2/troe_Ts[4] : 0.)
      + (troe_len[4] == 4 ? Fcent3*troe_Tss[4]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[18];
    Kc = refCinv * exp(g_RT[1] + g_RT[17] - g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[17]) + (h_RT[18]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[17] -= q; /* C2H2 */
    wdot[18] += q; /* C2H3 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[4][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[17] -= dqdci;               /* dwdot[C2H2]/d[H2] */
        J[18] += dqdci;               /* dwdot[C2H3]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[17];
        J[26] -= dqdci;               /* dwdot[H]/d[H] */
        J[42] -= dqdci;               /* dwdot[C2H2]/d[H] */
        J[43] += dqdci;               /* dwdot[C2H3]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[4][1] - 1)*dcdc_fac;
        J[126] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[142] -= dqdci;              /* dwdot[C2H2]/d[H2O] */
        J[143] += dqdci;              /* dwdot[C2H3]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[4][2] - 1)*dcdc_fac;
        J[276] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[292] -= dqdci;              /* dwdot[C2H2]/d[CH4] */
        J[293] += dqdci;              /* dwdot[C2H3]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[4][3] - 1)*dcdc_fac;
        J[301] -= dqdci;              /* dwdot[H]/d[CO] */
        J[317] -= dqdci;              /* dwdot[C2H2]/d[CO] */
        J[318] += dqdci;              /* dwdot[C2H3]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[4][4] - 1)*dcdc_fac;
        J[326] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[342] -= dqdci;              /* dwdot[C2H2]/d[CO2] */
        J[343] += dqdci;              /* dwdot[C2H3]/d[CO2] */
        /* d()/d[C2H2] */
        dqdci =  + k_f*sc[1];
        J[426] -= dqdci;              /* dwdot[H]/d[C2H2] */
        J[442] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
        J[443] += dqdci;              /* dwdot[C2H3]/d[C2H2] */
        /* d()/d[C2H3] */
        dqdci =  - k_r;
        J[451] -= dqdci;              /* dwdot[H]/d[C2H3] */
        J[467] -= dqdci;              /* dwdot[C2H2]/d[C2H3] */
        J[468] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
        /* d()/d[C2H6] */
        dqdci = (TB[4][5] - 1)*dcdc_fac;
        J[526] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[542] -= dqdci;              /* dwdot[C2H2]/d[C2H6] */
        J[543] += dqdci;              /* dwdot[C2H3]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[4][6] - 1)*dcdc_fac;
        J[576] -= dqdci;              /* dwdot[H]/d[AR] */
        J[592] -= dqdci;              /* dwdot[C2H2]/d[AR] */
        J[593] += dqdci;              /* dwdot[C2H3]/d[AR] */
    }
    else {
        dqdc[0] = TB[4][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[17];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[4][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[4][2]*dcdc_fac;
        dqdc[12] = TB[4][3]*dcdc_fac;
        dqdc[13] = TB[4][4]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac + k_f*sc[1];
        dqdc[18] = dcdc_fac - k_r;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = TB[4][5]*dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = TB[4][6]*dcdc_fac;
        for (int k=0; k<24; k++) {
            J[25*k+1] -= dqdc[k];
            J[25*k+17] -= dqdc[k];
            J[25*k+18] += dqdc[k];
        }
    }
    J[601] -= dqdT; /* dwdot[H]/dT */
    J[617] -= dqdT; /* dwdot[C2H2]/dT */
    J[618] += dqdT; /* dwdot[C2H3]/dT */

    /*reaction 6: H + C2H3 (+M) <=> C2H4 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[5][0] - 1)*sc[0] + (TB[5][1] - 1)*sc[5] + (TB[5][2] - 1)*sc[11] + (TB[5][3] - 1)*sc[12] + (TB[5][4] - 1)*sc[13] + (TB[5][5] - 1)*sc[21] + (TB[5][6] - 1)*sc[23];
    /* forward */
    phi_f = sc[1]*sc[18];
    k_f = prefactor_units[5] * fwd_A[5]
                * exp(fwd_beta[5] * tc[0] - activation_units[5] * fwd_Ea[5] * invT);
    dlnkfdT = fwd_beta[5] * invT + activation_units[5] * fwd_Ea[5] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[5] * exp(low_beta[5] * tc[0] - activation_units[5] * low_Ea[5] * invT);
    Pr = phase_units[5] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[5] * invT + activation_units[5] * low_Ea[5] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[5]) > 1.e-100 ? (1.-troe_a[5])*exp(-T/troe_Tsss[5]) : 0.);
    Fcent2 = (fabs(troe_Ts[5]) > 1.e-100 ? troe_a[5] * exp(-T/troe_Ts[5]) : 0.);
    Fcent3 = (troe_len[5] == 4 ? exp(-troe_Tss[5] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[5]) > 1.e-100 ? -Fcent1/troe_Tsss[5] : 0.)
      + (fabs(troe_Ts[5]) > 1.e-100 ? -Fcent2/troe_Ts[5] : 0.)
      + (troe_len[5] == 4 ? Fcent3*troe_Tss[5]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[19];
    Kc = refCinv * exp(g_RT[1] + g_RT[18] - g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[18]) + (h_RT[19]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[18] -= q; /* C2H3 */
    wdot[19] += q; /* C2H4 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[5][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[18] -= dqdci;               /* dwdot[C2H3]/d[H2] */
        J[19] += dqdci;               /* dwdot[C2H4]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[18];
        J[26] -= dqdci;               /* dwdot[H]/d[H] */
        J[43] -= dqdci;               /* dwdot[C2H3]/d[H] */
        J[44] += dqdci;               /* dwdot[C2H4]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[5][1] - 1)*dcdc_fac;
        J[126] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[143] -= dqdci;              /* dwdot[C2H3]/d[H2O] */
        J[144] += dqdci;              /* dwdot[C2H4]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[5][2] - 1)*dcdc_fac;
        J[276] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[293] -= dqdci;              /* dwdot[C2H3]/d[CH4] */
        J[294] += dqdci;              /* dwdot[C2H4]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[5][3] - 1)*dcdc_fac;
        J[301] -= dqdci;              /* dwdot[H]/d[CO] */
        J[318] -= dqdci;              /* dwdot[C2H3]/d[CO] */
        J[319] += dqdci;              /* dwdot[C2H4]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[5][4] - 1)*dcdc_fac;
        J[326] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[343] -= dqdci;              /* dwdot[C2H3]/d[CO2] */
        J[344] += dqdci;              /* dwdot[C2H4]/d[CO2] */
        /* d()/d[C2H3] */
        dqdci =  + k_f*sc[1];
        J[451] -= dqdci;              /* dwdot[H]/d[C2H3] */
        J[468] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
        J[469] += dqdci;              /* dwdot[C2H4]/d[C2H3] */
        /* d()/d[C2H4] */
        dqdci =  - k_r;
        J[476] -= dqdci;              /* dwdot[H]/d[C2H4] */
        J[493] -= dqdci;              /* dwdot[C2H3]/d[C2H4] */
        J[494] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
        /* d()/d[C2H6] */
        dqdci = (TB[5][5] - 1)*dcdc_fac;
        J[526] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[543] -= dqdci;              /* dwdot[C2H3]/d[C2H6] */
        J[544] += dqdci;              /* dwdot[C2H4]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[5][6] - 1)*dcdc_fac;
        J[576] -= dqdci;              /* dwdot[H]/d[AR] */
        J[593] -= dqdci;              /* dwdot[C2H3]/d[AR] */
        J[594] += dqdci;              /* dwdot[C2H4]/d[AR] */
    }
    else {
        dqdc[0] = TB[5][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[18];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[5][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[5][2]*dcdc_fac;
        dqdc[12] = TB[5][3]*dcdc_fac;
        dqdc[13] = TB[5][4]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac + k_f*sc[1];
        dqdc[19] = dcdc_fac - k_r;
        dqdc[20] = dcdc_fac;
        dqdc[21] = TB[5][5]*dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = TB[5][6]*dcdc_fac;
        for (int k=0; k<24; k++) {
            J[25*k+1] -= dqdc[k];
            J[25*k+18] -= dqdc[k];
            J[25*k+19] += dqdc[k];
        }
    }
    J[601] -= dqdT; /* dwdot[H]/dT */
    J[618] -= dqdT; /* dwdot[C2H3]/dT */
    J[619] += dqdT; /* dwdot[C2H4]/dT */

    /*reaction 7: H + C2H4 (+M) <=> C2H5 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[6][0] - 1)*sc[0] + (TB[6][1] - 1)*sc[5] + (TB[6][2] - 1)*sc[11] + (TB[6][3] - 1)*sc[12] + (TB[6][4] - 1)*sc[13] + (TB[6][5] - 1)*sc[21] + (TB[6][6] - 1)*sc[23];
    /* forward */
    phi_f = sc[1]*sc[19];
    k_f = prefactor_units[6] * fwd_A[6]
                * exp(fwd_beta[6] * tc[0] - activation_units[6] * fwd_Ea[6] * invT);
    dlnkfdT = fwd_beta[6] * invT + activation_units[6] * fwd_Ea[6] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[6] * exp(low_beta[6] * tc[0] - activation_units[6] * low_Ea[6] * invT);
    Pr = phase_units[6] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[6] * invT + activation_units[6] * low_Ea[6] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[6]) > 1.e-100 ? (1.-troe_a[6])*exp(-T/troe_Tsss[6]) : 0.);
    Fcent2 = (fabs(troe_Ts[6]) > 1.e-100 ? troe_a[6] * exp(-T/troe_Ts[6]) : 0.);
    Fcent3 = (troe_len[6] == 4 ? exp(-troe_Tss[6] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[6]) > 1.e-100 ? -Fcent1/troe_Tsss[6] : 0.)
      + (fabs(troe_Ts[6]) > 1.e-100 ? -Fcent2/troe_Ts[6] : 0.)
      + (troe_len[6] == 4 ? Fcent3*troe_Tss[6]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[20];
    Kc = refCinv * exp(g_RT[1] + g_RT[19] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[19]) + (h_RT[20]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[19] -= q; /* C2H4 */
    wdot[20] += q; /* C2H5 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[6][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[19] -= dqdci;               /* dwdot[C2H4]/d[H2] */
        J[20] += dqdci;               /* dwdot[C2H5]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[19];
        J[26] -= dqdci;               /* dwdot[H]/d[H] */
        J[44] -= dqdci;               /* dwdot[C2H4]/d[H] */
        J[45] += dqdci;               /* dwdot[C2H5]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[6][1] - 1)*dcdc_fac;
        J[126] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[144] -= dqdci;              /* dwdot[C2H4]/d[H2O] */
        J[145] += dqdci;              /* dwdot[C2H5]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[6][2] - 1)*dcdc_fac;
        J[276] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[294] -= dqdci;              /* dwdot[C2H4]/d[CH4] */
        J[295] += dqdci;              /* dwdot[C2H5]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[6][3] - 1)*dcdc_fac;
        J[301] -= dqdci;              /* dwdot[H]/d[CO] */
        J[319] -= dqdci;              /* dwdot[C2H4]/d[CO] */
        J[320] += dqdci;              /* dwdot[C2H5]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[6][4] - 1)*dcdc_fac;
        J[326] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[344] -= dqdci;              /* dwdot[C2H4]/d[CO2] */
        J[345] += dqdci;              /* dwdot[C2H5]/d[CO2] */
        /* d()/d[C2H4] */
        dqdci =  + k_f*sc[1];
        J[476] -= dqdci;              /* dwdot[H]/d[C2H4] */
        J[494] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
        J[495] += dqdci;              /* dwdot[C2H5]/d[C2H4] */
        /* d()/d[C2H5] */
        dqdci =  - k_r;
        J[501] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[519] -= dqdci;              /* dwdot[C2H4]/d[C2H5] */
        J[520] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
        /* d()/d[C2H6] */
        dqdci = (TB[6][5] - 1)*dcdc_fac;
        J[526] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[544] -= dqdci;              /* dwdot[C2H4]/d[C2H6] */
        J[545] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[6][6] - 1)*dcdc_fac;
        J[576] -= dqdci;              /* dwdot[H]/d[AR] */
        J[594] -= dqdci;              /* dwdot[C2H4]/d[AR] */
        J[595] += dqdci;              /* dwdot[C2H5]/d[AR] */
    }
    else {
        dqdc[0] = TB[6][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[19];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[6][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[6][2]*dcdc_fac;
        dqdc[12] = TB[6][3]*dcdc_fac;
        dqdc[13] = TB[6][4]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac + k_f*sc[1];
        dqdc[20] = dcdc_fac - k_r;
        dqdc[21] = TB[6][5]*dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = TB[6][6]*dcdc_fac;
        for (int k=0; k<24; k++) {
            J[25*k+1] -= dqdc[k];
            J[25*k+19] -= dqdc[k];
            J[25*k+20] += dqdc[k];
        }
    }
    J[601] -= dqdT; /* dwdot[H]/dT */
    J[619] -= dqdT; /* dwdot[C2H4]/dT */
    J[620] += dqdT; /* dwdot[C2H5]/dT */

    /*reaction 8: H + C2H5 (+M) <=> C2H6 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[7][0] - 1)*sc[0] + (TB[7][1] - 1)*sc[5] + (TB[7][2] - 1)*sc[11] + (TB[7][3] - 1)*sc[12] + (TB[7][4] - 1)*sc[13] + (TB[7][5] - 1)*sc[21] + (TB[7][6] - 1)*sc[23];
    /* forward */
    phi_f = sc[1]*sc[20];
    k_f = prefactor_units[7] * fwd_A[7]
                * exp(fwd_beta[7] * tc[0] - activation_units[7] * fwd_Ea[7] * invT);
    dlnkfdT = fwd_beta[7] * invT + activation_units[7] * fwd_Ea[7] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[7] * exp(low_beta[7] * tc[0] - activation_units[7] * low_Ea[7] * invT);
    Pr = phase_units[7] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[7] * invT + activation_units[7] * low_Ea[7] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[7]) > 1.e-100 ? (1.-troe_a[7])*exp(-T/troe_Tsss[7]) : 0.);
    Fcent2 = (fabs(troe_Ts[7]) > 1.e-100 ? troe_a[7] * exp(-T/troe_Ts[7]) : 0.);
    Fcent3 = (troe_len[7] == 4 ? exp(-troe_Tss[7] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[7]) > 1.e-100 ? -Fcent1/troe_Tsss[7] : 0.)
      + (fabs(troe_Ts[7]) > 1.e-100 ? -Fcent2/troe_Ts[7] : 0.)
      + (troe_len[7] == 4 ? Fcent3*troe_Tss[7]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[21];
    Kc = refCinv * exp(g_RT[1] + g_RT[20] - g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[20]) + (h_RT[21]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[20] -= q; /* C2H5 */
    wdot[21] += q; /* C2H6 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[7][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[20] -= dqdci;               /* dwdot[C2H5]/d[H2] */
        J[21] += dqdci;               /* dwdot[C2H6]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[20];
        J[26] -= dqdci;               /* dwdot[H]/d[H] */
        J[45] -= dqdci;               /* dwdot[C2H5]/d[H] */
        J[46] += dqdci;               /* dwdot[C2H6]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[7][1] - 1)*dcdc_fac;
        J[126] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[145] -= dqdci;              /* dwdot[C2H5]/d[H2O] */
        J[146] += dqdci;              /* dwdot[C2H6]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[7][2] - 1)*dcdc_fac;
        J[276] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[295] -= dqdci;              /* dwdot[C2H5]/d[CH4] */
        J[296] += dqdci;              /* dwdot[C2H6]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[7][3] - 1)*dcdc_fac;
        J[301] -= dqdci;              /* dwdot[H]/d[CO] */
        J[320] -= dqdci;              /* dwdot[C2H5]/d[CO] */
        J[321] += dqdci;              /* dwdot[C2H6]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[7][4] - 1)*dcdc_fac;
        J[326] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[345] -= dqdci;              /* dwdot[C2H5]/d[CO2] */
        J[346] += dqdci;              /* dwdot[C2H6]/d[CO2] */
        /* d()/d[C2H5] */
        dqdci =  + k_f*sc[1];
        J[501] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[520] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
        J[521] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
        /* d()/d[C2H6] */
        dqdci = (TB[7][5] - 1)*dcdc_fac - k_r;
        J[526] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[545] -= dqdci;              /* dwdot[C2H5]/d[C2H6] */
        J[546] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[7][6] - 1)*dcdc_fac;
        J[576] -= dqdci;              /* dwdot[H]/d[AR] */
        J[595] -= dqdci;              /* dwdot[C2H5]/d[AR] */
        J[596] += dqdci;              /* dwdot[C2H6]/d[AR] */
    }
    else {
        dqdc[0] = TB[7][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[20];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[7][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[7][2]*dcdc_fac;
        dqdc[12] = TB[7][3]*dcdc_fac;
        dqdc[13] = TB[7][4]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac + k_f*sc[1];
        dqdc[21] = TB[7][5]*dcdc_fac - k_r;
        dqdc[22] = dcdc_fac;
        dqdc[23] = TB[7][6]*dcdc_fac;
        for (int k=0; k<24; k++) {
            J[25*k+1] -= dqdc[k];
            J[25*k+20] -= dqdc[k];
            J[25*k+21] += dqdc[k];
        }
    }
    J[601] -= dqdT; /* dwdot[H]/dT */
    J[620] -= dqdT; /* dwdot[C2H5]/dT */
    J[621] += dqdT; /* dwdot[C2H6]/dT */

    /*reaction 9: H2 + CO (+M) <=> CH2O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[8][0] - 1)*sc[0] + (TB[8][1] - 1)*sc[5] + (TB[8][2] - 1)*sc[11] + (TB[8][3] - 1)*sc[12] + (TB[8][4] - 1)*sc[13] + (TB[8][5] - 1)*sc[21] + (TB[8][6] - 1)*sc[23];
    /* forward */
    phi_f = sc[0]*sc[12];
    k_f = prefactor_units[8] * fwd_A[8]
                * exp(fwd_beta[8] * tc[0] - activation_units[8] * fwd_Ea[8] * invT);
    dlnkfdT = fwd_beta[8] * invT + activation_units[8] * fwd_Ea[8] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[8] * exp(low_beta[8] * tc[0] - activation_units[8] * low_Ea[8] * invT);
    Pr = phase_units[8] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[8] * invT + activation_units[8] * low_Ea[8] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[8]) > 1.e-100 ? (1.-troe_a[8])*exp(-T/troe_Tsss[8]) : 0.);
    Fcent2 = (fabs(troe_Ts[8]) > 1.e-100 ? troe_a[8] * exp(-T/troe_Ts[8]) : 0.);
    Fcent3 = (troe_len[8] == 4 ? exp(-troe_Tss[8] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[8]) > 1.e-100 ? -Fcent1/troe_Tsss[8] : 0.)
      + (fabs(troe_Ts[8]) > 1.e-100 ? -Fcent2/troe_Ts[8] : 0.)
      + (troe_len[8] == 4 ? Fcent3*troe_Tss[8]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[15];
    Kc = refCinv * exp(g_RT[0] + g_RT[12] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[12]) + (h_RT[15]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[12] -= q; /* CO */
    wdot[15] += q; /* CH2O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[8][0] - 1)*dcdc_fac + k_f*sc[12];
        J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
        J[12] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[15] += dqdci;               /* dwdot[CH2O]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[8][1] - 1)*dcdc_fac;
        J[125] -= dqdci;              /* dwdot[H2]/d[H2O] */
        J[137] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[140] += dqdci;              /* dwdot[CH2O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[8][2] - 1)*dcdc_fac;
        J[275] -= dqdci;              /* dwdot[H2]/d[CH4] */
        J[287] -= dqdci;              /* dwdot[CO]/d[CH4] */
        J[290] += dqdci;              /* dwdot[CH2O]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[8][3] - 1)*dcdc_fac + k_f*sc[0];
        J[300] -= dqdci;              /* dwdot[H2]/d[CO] */
        J[312] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[315] += dqdci;              /* dwdot[CH2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[8][4] - 1)*dcdc_fac;
        J[325] -= dqdci;              /* dwdot[H2]/d[CO2] */
        J[337] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[340] += dqdci;              /* dwdot[CH2O]/d[CO2] */
        /* d()/d[CH2O] */
        dqdci =  - k_r;
        J[375] -= dqdci;              /* dwdot[H2]/d[CH2O] */
        J[387] -= dqdci;              /* dwdot[CO]/d[CH2O] */
        J[390] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
        /* d()/d[C2H6] */
        dqdci = (TB[8][5] - 1)*dcdc_fac;
        J[525] -= dqdci;              /* dwdot[H2]/d[C2H6] */
        J[537] -= dqdci;              /* dwdot[CO]/d[C2H6] */
        J[540] += dqdci;              /* dwdot[CH2O]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[8][6] - 1)*dcdc_fac;
        J[575] -= dqdci;              /* dwdot[H2]/d[AR] */
        J[587] -= dqdci;              /* dwdot[CO]/d[AR] */
        J[590] += dqdci;              /* dwdot[CH2O]/d[AR] */
    }
    else {
        dqdc[0] = TB[8][0]*dcdc_fac + k_f*sc[12];
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[8][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[8][2]*dcdc_fac;
        dqdc[12] = TB[8][3]*dcdc_fac + k_f*sc[0];
        dqdc[13] = TB[8][4]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac - k_r;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = TB[8][5]*dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = TB[8][6]*dcdc_fac;
        for (int k=0; k<24; k++) {
            J[25*k+0] -= dqdc[k];
            J[25*k+12] -= dqdc[k];
            J[25*k+15] += dqdc[k];
        }
    }
    J[600] -= dqdT; /* dwdot[H2]/dT */
    J[612] -= dqdT; /* dwdot[CO]/dT */
    J[615] += dqdT; /* dwdot[CH2O]/dT */

    /*reaction 10: 2 OH (+M) <=> H2O2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[9][0] - 1)*sc[0] + (TB[9][1] - 1)*sc[5] + (TB[9][2] - 1)*sc[11] + (TB[9][3] - 1)*sc[12] + (TB[9][4] - 1)*sc[13] + (TB[9][5] - 1)*sc[21] + (TB[9][6] - 1)*sc[23];
    /* forward */
    phi_f = sc[4]*sc[4];
    k_f = prefactor_units[9] * fwd_A[9]
                * exp(fwd_beta[9] * tc[0] - activation_units[9] * fwd_Ea[9] * invT);
    dlnkfdT = fwd_beta[9] * invT + activation_units[9] * fwd_Ea[9] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[9] * exp(low_beta[9] * tc[0] - activation_units[9] * low_Ea[9] * invT);
    Pr = phase_units[9] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[9] * invT + activation_units[9] * low_Ea[9] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[9]) > 1.e-100 ? (1.-troe_a[9])*exp(-T/troe_Tsss[9]) : 0.);
    Fcent2 = (fabs(troe_Ts[9]) > 1.e-100 ? troe_a[9] * exp(-T/troe_Ts[9]) : 0.);
    Fcent3 = (troe_len[9] == 4 ? exp(-troe_Tss[9] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[9]) > 1.e-100 ? -Fcent1/troe_Tsss[9] : 0.)
      + (fabs(troe_Ts[9]) > 1.e-100 ? -Fcent2/troe_Ts[9] : 0.)
      + (troe_len[9] == 4 ? Fcent3*troe_Tss[9]*invT2 : 0.) );
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
        dqdci = (TB[9][0] - 1)*dcdc_fac;
        J[4] += -2 * dqdci;           /* dwdot[OH]/d[H2] */
        J[7] += dqdci;                /* dwdot[H2O2]/d[H2] */
        /* d()/d[OH] */
        dqdci =  + k_f*2*sc[4];
        J[104] += -2 * dqdci;         /* dwdot[OH]/d[OH] */
        J[107] += dqdci;              /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[9][1] - 1)*dcdc_fac;
        J[129] += -2 * dqdci;         /* dwdot[OH]/d[H2O] */
        J[132] += dqdci;              /* dwdot[H2O2]/d[H2O] */
        /* d()/d[H2O2] */
        dqdci =  - k_r;
        J[179] += -2 * dqdci;         /* dwdot[OH]/d[H2O2] */
        J[182] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
        /* d()/d[CH4] */
        dqdci = (TB[9][2] - 1)*dcdc_fac;
        J[279] += -2 * dqdci;         /* dwdot[OH]/d[CH4] */
        J[282] += dqdci;              /* dwdot[H2O2]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[9][3] - 1)*dcdc_fac;
        J[304] += -2 * dqdci;         /* dwdot[OH]/d[CO] */
        J[307] += dqdci;              /* dwdot[H2O2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[9][4] - 1)*dcdc_fac;
        J[329] += -2 * dqdci;         /* dwdot[OH]/d[CO2] */
        J[332] += dqdci;              /* dwdot[H2O2]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[9][5] - 1)*dcdc_fac;
        J[529] += -2 * dqdci;         /* dwdot[OH]/d[C2H6] */
        J[532] += dqdci;              /* dwdot[H2O2]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[9][6] - 1)*dcdc_fac;
        J[579] += -2 * dqdci;         /* dwdot[OH]/d[AR] */
        J[582] += dqdci;              /* dwdot[H2O2]/d[AR] */
    }
    else {
        dqdc[0] = TB[9][0]*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac + k_f*2*sc[4];
        dqdc[5] = TB[9][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac - k_r;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[9][2]*dcdc_fac;
        dqdc[12] = TB[9][3]*dcdc_fac;
        dqdc[13] = TB[9][4]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = TB[9][5]*dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = TB[9][6]*dcdc_fac;
        for (int k=0; k<24; k++) {
            J[25*k+4] += -2 * dqdc[k];
            J[25*k+7] += dqdc[k];
        }
    }
    J[604] += -2 * dqdT; /* dwdot[OH]/dT */
    J[607] += dqdT; /* dwdot[H2O2]/dT */

    /*reaction 11: 2 CH3 (+M) <=> C2H6 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[10][0] - 1)*sc[0] + (TB[10][1] - 1)*sc[5] + (TB[10][2] - 1)*sc[11] + (TB[10][3] - 1)*sc[12] + (TB[10][4] - 1)*sc[13] + (TB[10][5] - 1)*sc[21] + (TB[10][6] - 1)*sc[23];
    /* forward */
    phi_f = sc[10]*sc[10];
    k_f = prefactor_units[10] * fwd_A[10]
                * exp(fwd_beta[10] * tc[0] - activation_units[10] * fwd_Ea[10] * invT);
    dlnkfdT = fwd_beta[10] * invT + activation_units[10] * fwd_Ea[10] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[10] * exp(low_beta[10] * tc[0] - activation_units[10] * low_Ea[10] * invT);
    Pr = phase_units[10] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[10] * invT + activation_units[10] * low_Ea[10] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[10]) > 1.e-100 ? (1.-troe_a[10])*exp(-T/troe_Tsss[10]) : 0.);
    Fcent2 = (fabs(troe_Ts[10]) > 1.e-100 ? troe_a[10] * exp(-T/troe_Ts[10]) : 0.);
    Fcent3 = (troe_len[10] == 4 ? exp(-troe_Tss[10] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[10]) > 1.e-100 ? -Fcent1/troe_Tsss[10] : 0.)
      + (fabs(troe_Ts[10]) > 1.e-100 ? -Fcent2/troe_Ts[10] : 0.)
      + (troe_len[10] == 4 ? Fcent3*troe_Tss[10]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[21];
    Kc = refCinv * exp(2*g_RT[10] - g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[10]) + (h_RT[21]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[10] -= 2 * q; /* CH3 */
    wdot[21] += q; /* C2H6 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[10][0] - 1)*dcdc_fac;
        J[10] += -2 * dqdci;          /* dwdot[CH3]/d[H2] */
        J[21] += dqdci;               /* dwdot[C2H6]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[10][1] - 1)*dcdc_fac;
        J[135] += -2 * dqdci;         /* dwdot[CH3]/d[H2O] */
        J[146] += dqdci;              /* dwdot[C2H6]/d[H2O] */
        /* d()/d[CH3] */
        dqdci =  + k_f*2*sc[10];
        J[260] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
        J[271] += dqdci;              /* dwdot[C2H6]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (TB[10][2] - 1)*dcdc_fac;
        J[285] += -2 * dqdci;         /* dwdot[CH3]/d[CH4] */
        J[296] += dqdci;              /* dwdot[C2H6]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[10][3] - 1)*dcdc_fac;
        J[310] += -2 * dqdci;         /* dwdot[CH3]/d[CO] */
        J[321] += dqdci;              /* dwdot[C2H6]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[10][4] - 1)*dcdc_fac;
        J[335] += -2 * dqdci;         /* dwdot[CH3]/d[CO2] */
        J[346] += dqdci;              /* dwdot[C2H6]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[10][5] - 1)*dcdc_fac - k_r;
        J[535] += -2 * dqdci;         /* dwdot[CH3]/d[C2H6] */
        J[546] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[10][6] - 1)*dcdc_fac;
        J[585] += -2 * dqdci;         /* dwdot[CH3]/d[AR] */
        J[596] += dqdci;              /* dwdot[C2H6]/d[AR] */
    }
    else {
        dqdc[0] = TB[10][0]*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[10][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac + k_f*2*sc[10];
        dqdc[11] = TB[10][2]*dcdc_fac;
        dqdc[12] = TB[10][3]*dcdc_fac;
        dqdc[13] = TB[10][4]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = TB[10][5]*dcdc_fac - k_r;
        dqdc[22] = dcdc_fac;
        dqdc[23] = TB[10][6]*dcdc_fac;
        for (int k=0; k<24; k++) {
            J[25*k+10] += -2 * dqdc[k];
            J[25*k+21] += dqdc[k];
        }
    }
    J[610] += -2 * dqdT; /* dwdot[CH3]/dT */
    J[621] += dqdT; /* dwdot[C2H6]/dT */

    /*reaction 12: C2H4 (+M) <=> H2 + C2H2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[11][0] - 1)*sc[0] + (TB[11][1] - 1)*sc[5] + (TB[11][2] - 1)*sc[11] + (TB[11][3] - 1)*sc[12] + (TB[11][4] - 1)*sc[13] + (TB[11][5] - 1)*sc[21] + (TB[11][6] - 1)*sc[23];
    /* forward */
    phi_f = sc[19];
    k_f = prefactor_units[11] * fwd_A[11]
                * exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
    dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[11] * exp(low_beta[11] * tc[0] - activation_units[11] * low_Ea[11] * invT);
    Pr = phase_units[11] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[11] * invT + activation_units[11] * low_Ea[11] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[11]) > 1.e-100 ? (1.-troe_a[11])*exp(-T/troe_Tsss[11]) : 0.);
    Fcent2 = (fabs(troe_Ts[11]) > 1.e-100 ? troe_a[11] * exp(-T/troe_Ts[11]) : 0.);
    Fcent3 = (troe_len[11] == 4 ? exp(-troe_Tss[11] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[11]) > 1.e-100 ? -Fcent1/troe_Tsss[11] : 0.)
      + (fabs(troe_Ts[11]) > 1.e-100 ? -Fcent2/troe_Ts[11] : 0.)
      + (troe_len[11] == 4 ? Fcent3*troe_Tss[11]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[0]*sc[17];
    Kc = refC * exp(-g_RT[0] - g_RT[17] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[19]) + (h_RT[0] + h_RT[17]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[17] += q; /* C2H2 */
    wdot[19] -= q; /* C2H4 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[11][0] - 1)*dcdc_fac - k_r*sc[17];
        J[0] += dqdci;                /* dwdot[H2]/d[H2] */
        J[17] += dqdci;               /* dwdot[C2H2]/d[H2] */
        J[19] -= dqdci;               /* dwdot[C2H4]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[11][1] - 1)*dcdc_fac;
        J[125] += dqdci;              /* dwdot[H2]/d[H2O] */
        J[142] += dqdci;              /* dwdot[C2H2]/d[H2O] */
        J[144] -= dqdci;              /* dwdot[C2H4]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[11][2] - 1)*dcdc_fac;
        J[275] += dqdci;              /* dwdot[H2]/d[CH4] */
        J[292] += dqdci;              /* dwdot[C2H2]/d[CH4] */
        J[294] -= dqdci;              /* dwdot[C2H4]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[11][3] - 1)*dcdc_fac;
        J[300] += dqdci;              /* dwdot[H2]/d[CO] */
        J[317] += dqdci;              /* dwdot[C2H2]/d[CO] */
        J[319] -= dqdci;              /* dwdot[C2H4]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[11][4] - 1)*dcdc_fac;
        J[325] += dqdci;              /* dwdot[H2]/d[CO2] */
        J[342] += dqdci;              /* dwdot[C2H2]/d[CO2] */
        J[344] -= dqdci;              /* dwdot[C2H4]/d[CO2] */
        /* d()/d[C2H2] */
        dqdci =  - k_r*sc[0];
        J[425] += dqdci;              /* dwdot[H2]/d[C2H2] */
        J[442] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
        J[444] -= dqdci;              /* dwdot[C2H4]/d[C2H2] */
        /* d()/d[C2H4] */
        dqdci =  + k_f;
        J[475] += dqdci;              /* dwdot[H2]/d[C2H4] */
        J[492] += dqdci;              /* dwdot[C2H2]/d[C2H4] */
        J[494] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
        /* d()/d[C2H6] */
        dqdci = (TB[11][5] - 1)*dcdc_fac;
        J[525] += dqdci;              /* dwdot[H2]/d[C2H6] */
        J[542] += dqdci;              /* dwdot[C2H2]/d[C2H6] */
        J[544] -= dqdci;              /* dwdot[C2H4]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[11][6] - 1)*dcdc_fac;
        J[575] += dqdci;              /* dwdot[H2]/d[AR] */
        J[592] += dqdci;              /* dwdot[C2H2]/d[AR] */
        J[594] -= dqdci;              /* dwdot[C2H4]/d[AR] */
    }
    else {
        dqdc[0] = TB[11][0]*dcdc_fac - k_r*sc[17];
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[11][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[11][2]*dcdc_fac;
        dqdc[12] = TB[11][3]*dcdc_fac;
        dqdc[13] = TB[11][4]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac - k_r*sc[0];
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac + k_f;
        dqdc[20] = dcdc_fac;
        dqdc[21] = TB[11][5]*dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = TB[11][6]*dcdc_fac;
        for (int k=0; k<24; k++) {
            J[25*k+0] += dqdc[k];
            J[25*k+17] += dqdc[k];
            J[25*k+19] -= dqdc[k];
        }
    }
    J[600] += dqdT; /* dwdot[H2]/dT */
    J[617] += dqdT; /* dwdot[C2H2]/dT */
    J[619] -= dqdT; /* dwdot[C2H4]/dT */

    /*reaction 13: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[12][0] - 1)*sc[0] + (TB[12][1] - 1)*sc[5] + (TB[12][2] - 1)*sc[11] + (TB[12][3] - 1)*sc[12] + (TB[12][4] - 1)*sc[13] + (TB[12][5] - 1)*sc[21] + (TB[12][6] - 1)*sc[23];
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = prefactor_units[12] * fwd_A[12]
                * exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
    dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
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
        dqdci = (TB[12][0] - 1)*q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[2] -= dqdci;                /* dwdot[O]/d[H2] */
        J[4] += dqdci;                /* dwdot[OH]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[2];
        J[26] -= dqdci;               /* dwdot[H]/d[H] */
        J[27] -= dqdci;               /* dwdot[O]/d[H] */
        J[29] += dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[1];
        J[51] -= dqdci;               /* dwdot[H]/d[O] */
        J[52] -= dqdci;               /* dwdot[O]/d[O] */
        J[54] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[101] -= dqdci;              /* dwdot[H]/d[OH] */
        J[102] -= dqdci;              /* dwdot[O]/d[OH] */
        J[104] += dqdci;              /* dwdot[OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[12][1] - 1)*q_nocor;
        J[126] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[127] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[129] += dqdci;              /* dwdot[OH]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[12][2] - 1)*q_nocor;
        J[276] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[277] -= dqdci;              /* dwdot[O]/d[CH4] */
        J[279] += dqdci;              /* dwdot[OH]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[12][3] - 1)*q_nocor;
        J[301] -= dqdci;              /* dwdot[H]/d[CO] */
        J[302] -= dqdci;              /* dwdot[O]/d[CO] */
        J[304] += dqdci;              /* dwdot[OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[12][4] - 1)*q_nocor;
        J[326] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[327] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[329] += dqdci;              /* dwdot[OH]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[12][5] - 1)*q_nocor;
        J[526] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[527] -= dqdci;              /* dwdot[O]/d[C2H6] */
        J[529] += dqdci;              /* dwdot[OH]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[12][6] - 1)*q_nocor;
        J[576] -= dqdci;              /* dwdot[H]/d[AR] */
        J[577] -= dqdci;              /* dwdot[O]/d[AR] */
        J[579] += dqdci;              /* dwdot[OH]/d[AR] */
    }
    else {
        dqdc[0] = TB[12][0];
        dqdc[1] = dcdc_fac + k_f*sc[2];
        dqdc[2] = dcdc_fac + k_f*sc[1];
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac - k_r;
        dqdc[5] = TB[12][1];
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[12][2];
        dqdc[12] = TB[12][3];
        dqdc[13] = TB[12][4];
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = TB[12][5];
        dqdc[22] = dcdc_fac;
        dqdc[23] = TB[12][6];
        for (int k=0; k<24; k++) {
            J[25*k+1] -= dqdc[k];
            J[25*k+2] -= dqdc[k];
            J[25*k+4] += dqdc[k];
        }
    }
    J[601] -= dqdT; /* dwdot[H]/dT */
    J[602] -= dqdT; /* dwdot[O]/dT */
    J[604] += dqdT; /* dwdot[OH]/dT */

    /*reaction 14: O + CO + M <=> CO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[13][0] - 1)*sc[0] + (TB[13][1] - 1)*sc[3] + (TB[13][2] - 1)*sc[5] + (TB[13][3] - 1)*sc[11] + (TB[13][4] - 1)*sc[12] + (TB[13][5] - 1)*sc[13] + (TB[13][6] - 1)*sc[21] + (TB[13][7] - 1)*sc[23];
    /* forward */
    phi_f = sc[2]*sc[12];
    k_f = prefactor_units[13] * fwd_A[13]
                * exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
    dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
    /* reverse */
    phi_r = sc[13];
    Kc = refCinv * exp(g_RT[2] + g_RT[12] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[12]) + (h_RT[13]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[12] -= q; /* CO */
    wdot[13] += q; /* CO2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[13][0] - 1)*q_nocor;
        J[2] -= dqdci;                /* dwdot[O]/d[H2] */
        J[12] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[13] += dqdci;               /* dwdot[CO2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[12];
        J[52] -= dqdci;               /* dwdot[O]/d[O] */
        J[62] -= dqdci;               /* dwdot[CO]/d[O] */
        J[63] += dqdci;               /* dwdot[CO2]/d[O] */
        /* d()/d[O2] */
        dqdci = (TB[13][1] - 1)*q_nocor;
        J[77] -= dqdci;               /* dwdot[O]/d[O2] */
        J[87] -= dqdci;               /* dwdot[CO]/d[O2] */
        J[88] += dqdci;               /* dwdot[CO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[13][2] - 1)*q_nocor;
        J[127] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[137] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[138] += dqdci;              /* dwdot[CO2]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[13][3] - 1)*q_nocor;
        J[277] -= dqdci;              /* dwdot[O]/d[CH4] */
        J[287] -= dqdci;              /* dwdot[CO]/d[CH4] */
        J[288] += dqdci;              /* dwdot[CO2]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[13][4] - 1)*q_nocor + k_f*sc[2];
        J[302] -= dqdci;              /* dwdot[O]/d[CO] */
        J[312] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[313] += dqdci;              /* dwdot[CO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[13][5] - 1)*q_nocor - k_r;
        J[327] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[337] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[338] += dqdci;              /* dwdot[CO2]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[13][6] - 1)*q_nocor;
        J[527] -= dqdci;              /* dwdot[O]/d[C2H6] */
        J[537] -= dqdci;              /* dwdot[CO]/d[C2H6] */
        J[538] += dqdci;              /* dwdot[CO2]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[13][7] - 1)*q_nocor;
        J[577] -= dqdci;              /* dwdot[O]/d[AR] */
        J[587] -= dqdci;              /* dwdot[CO]/d[AR] */
        J[588] += dqdci;              /* dwdot[CO2]/d[AR] */
    }
    else {
        dqdc[0] = TB[13][0];
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac + k_f*sc[12];
        dqdc[3] = TB[13][1];
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[13][2];
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[13][3];
        dqdc[12] = TB[13][4] + k_f*sc[2];
        dqdc[13] = TB[13][5] - k_r;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = TB[13][6];
        dqdc[22] = dcdc_fac;
        dqdc[23] = TB[13][7];
        for (int k=0; k<24; k++) {
            J[25*k+2] -= dqdc[k];
            J[25*k+12] -= dqdc[k];
            J[25*k+13] += dqdc[k];
        }
    }
    J[602] -= dqdT; /* dwdot[O]/dT */
    J[612] -= dqdT; /* dwdot[CO]/dT */
    J[613] += dqdT; /* dwdot[CO2]/dT */

    /*reaction 15: H + O2 + M <=> HO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[14][0] - 1)*sc[3] + (TB[14][1] - 1)*sc[5] + (TB[14][2] - 1)*sc[12] + (TB[14][3] - 1)*sc[13] + (TB[14][4] - 1)*sc[21] + (TB[14][5] - 1)*sc[22] + (TB[14][6] - 1)*sc[23];
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[14] * fwd_A[14]
                * exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
    dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
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
        J[26] -= dqdci;               /* dwdot[H]/d[H] */
        J[28] -= dqdci;               /* dwdot[O2]/d[H] */
        J[31] += dqdci;               /* dwdot[HO2]/d[H] */
        /* d()/d[O2] */
        dqdci = (TB[14][0] - 1)*q_nocor + k_f*sc[1];
        J[76] -= dqdci;               /* dwdot[H]/d[O2] */
        J[78] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[81] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[14][1] - 1)*q_nocor;
        J[126] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[128] -= dqdci;              /* dwdot[O2]/d[H2O] */
        J[131] += dqdci;              /* dwdot[HO2]/d[H2O] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[151] -= dqdci;              /* dwdot[H]/d[HO2] */
        J[153] -= dqdci;              /* dwdot[O2]/d[HO2] */
        J[156] += dqdci;              /* dwdot[HO2]/d[HO2] */
        /* d()/d[CO] */
        dqdci = (TB[14][2] - 1)*q_nocor;
        J[301] -= dqdci;              /* dwdot[H]/d[CO] */
        J[303] -= dqdci;              /* dwdot[O2]/d[CO] */
        J[306] += dqdci;              /* dwdot[HO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[14][3] - 1)*q_nocor;
        J[326] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[328] -= dqdci;              /* dwdot[O2]/d[CO2] */
        J[331] += dqdci;              /* dwdot[HO2]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[14][4] - 1)*q_nocor;
        J[526] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[528] -= dqdci;              /* dwdot[O2]/d[C2H6] */
        J[531] += dqdci;              /* dwdot[HO2]/d[C2H6] */
        /* d()/d[N2] */
        dqdci = (TB[14][5] - 1)*q_nocor;
        J[551] -= dqdci;              /* dwdot[H]/d[N2] */
        J[553] -= dqdci;              /* dwdot[O2]/d[N2] */
        J[556] += dqdci;              /* dwdot[HO2]/d[N2] */
        /* d()/d[AR] */
        dqdci = (TB[14][6] - 1)*q_nocor;
        J[576] -= dqdci;              /* dwdot[H]/d[AR] */
        J[578] -= dqdci;              /* dwdot[O2]/d[AR] */
        J[581] += dqdci;              /* dwdot[HO2]/d[AR] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[3];
        dqdc[2] = dcdc_fac;
        dqdc[3] = TB[14][0] + k_f*sc[1];
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[14][1];
        dqdc[6] = dcdc_fac - k_r;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[14][2];
        dqdc[13] = TB[14][3];
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = TB[14][4];
        dqdc[22] = TB[14][5];
        dqdc[23] = TB[14][6];
        for (int k=0; k<24; k++) {
            J[25*k+1] -= dqdc[k];
            J[25*k+3] -= dqdc[k];
            J[25*k+6] += dqdc[k];
        }
    }
    J[601] -= dqdT; /* dwdot[H]/dT */
    J[603] -= dqdT; /* dwdot[O2]/dT */
    J[606] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 16: 2 H + M <=> H2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[15][0] - 1)*sc[0] + (TB[15][1] - 1)*sc[5] + (TB[15][2] - 1)*sc[11] + (TB[15][3] - 1)*sc[13] + (TB[15][4] - 1)*sc[21] + (TB[15][5] - 1)*sc[23];
    /* forward */
    phi_f = sc[1]*sc[1];
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
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
        dqdci = (TB[15][0] - 1)*q_nocor - k_r;
        J[0] += dqdci;                /* dwdot[H2]/d[H2] */
        J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*2*sc[1];
        J[25] += dqdci;               /* dwdot[H2]/d[H] */
        J[26] += -2 * dqdci;          /* dwdot[H]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[15][1] - 1)*q_nocor;
        J[125] += dqdci;              /* dwdot[H2]/d[H2O] */
        J[126] += -2 * dqdci;         /* dwdot[H]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[15][2] - 1)*q_nocor;
        J[275] += dqdci;              /* dwdot[H2]/d[CH4] */
        J[276] += -2 * dqdci;         /* dwdot[H]/d[CH4] */
        /* d()/d[CO2] */
        dqdci = (TB[15][3] - 1)*q_nocor;
        J[325] += dqdci;              /* dwdot[H2]/d[CO2] */
        J[326] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[15][4] - 1)*q_nocor;
        J[525] += dqdci;              /* dwdot[H2]/d[C2H6] */
        J[526] += -2 * dqdci;         /* dwdot[H]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[15][5] - 1)*q_nocor;
        J[575] += dqdci;              /* dwdot[H2]/d[AR] */
        J[576] += -2 * dqdci;         /* dwdot[H]/d[AR] */
    }
    else {
        dqdc[0] = TB[15][0] - k_r;
        dqdc[1] = dcdc_fac + k_f*2*sc[1];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[15][1];
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[15][2];
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[15][3];
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = TB[15][4];
        dqdc[22] = dcdc_fac;
        dqdc[23] = TB[15][5];
        for (int k=0; k<24; k++) {
            J[25*k+0] += dqdc[k];
            J[25*k+1] += -2 * dqdc[k];
        }
    }
    J[600] += dqdT; /* dwdot[H2]/dT */
    J[601] += -2 * dqdT; /* dwdot[H]/dT */

    /*reaction 17: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[16][0] - 1)*sc[0] + (TB[16][1] - 1)*sc[5] + (TB[16][2] - 1)*sc[11] + (TB[16][3] - 1)*sc[21] + (TB[16][4] - 1)*sc[23];
    /* forward */
    phi_f = sc[1]*sc[4];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
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
        dqdci = (TB[16][0] - 1)*q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[4] -= dqdci;                /* dwdot[OH]/d[H2] */
        J[5] += dqdci;                /* dwdot[H2O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[4];
        J[26] -= dqdci;               /* dwdot[H]/d[H] */
        J[29] -= dqdci;               /* dwdot[OH]/d[H] */
        J[30] += dqdci;               /* dwdot[H2O]/d[H] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[1];
        J[101] -= dqdci;              /* dwdot[H]/d[OH] */
        J[104] -= dqdci;              /* dwdot[OH]/d[OH] */
        J[105] += dqdci;              /* dwdot[H2O]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[16][1] - 1)*q_nocor - k_r;
        J[126] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[129] -= dqdci;              /* dwdot[OH]/d[H2O] */
        J[130] += dqdci;              /* dwdot[H2O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[16][2] - 1)*q_nocor;
        J[276] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[279] -= dqdci;              /* dwdot[OH]/d[CH4] */
        J[280] += dqdci;              /* dwdot[H2O]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (TB[16][3] - 1)*q_nocor;
        J[526] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[529] -= dqdci;              /* dwdot[OH]/d[C2H6] */
        J[530] += dqdci;              /* dwdot[H2O]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[16][4] - 1)*q_nocor;
        J[576] -= dqdci;              /* dwdot[H]/d[AR] */
        J[579] -= dqdci;              /* dwdot[OH]/d[AR] */
        J[580] += dqdci;              /* dwdot[H2O]/d[AR] */
    }
    else {
        dqdc[0] = TB[16][0];
        dqdc[1] = dcdc_fac + k_f*sc[4];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac + k_f*sc[1];
        dqdc[5] = TB[16][1] - k_r;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[16][2];
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = TB[16][3];
        dqdc[22] = dcdc_fac;
        dqdc[23] = TB[16][4];
        for (int k=0; k<24; k++) {
            J[25*k+1] -= dqdc[k];
            J[25*k+4] -= dqdc[k];
            J[25*k+5] += dqdc[k];
        }
    }
    J[601] -= dqdT; /* dwdot[H]/dT */
    J[604] -= dqdT; /* dwdot[OH]/dT */
    J[605] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 18: HCO + M <=> H + CO + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[17][0] - 1)*sc[0] + (TB[17][1] - 1)*sc[5] + (TB[17][2] - 1)*sc[11] + (TB[17][3] - 1)*sc[12] + (TB[17][4] - 1)*sc[13] + (TB[17][5] - 1)*sc[21];
    /* forward */
    phi_f = sc[14];
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[12];
    Kc = refC * exp(-g_RT[1] - g_RT[12] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[14]) + (h_RT[1] + h_RT[12]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[12] += q; /* CO */
    wdot[14] -= q; /* HCO */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[17][0] - 1)*q_nocor;
        J[1] += dqdci;                /* dwdot[H]/d[H2] */
        J[12] += dqdci;               /* dwdot[CO]/d[H2] */
        J[14] -= dqdci;               /* dwdot[HCO]/d[H2] */
        /* d()/d[H] */
        dqdci =  - k_r*sc[12];
        J[26] += dqdci;               /* dwdot[H]/d[H] */
        J[37] += dqdci;               /* dwdot[CO]/d[H] */
        J[39] -= dqdci;               /* dwdot[HCO]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[17][1] - 1)*q_nocor;
        J[126] += dqdci;              /* dwdot[H]/d[H2O] */
        J[137] += dqdci;              /* dwdot[CO]/d[H2O] */
        J[139] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[17][2] - 1)*q_nocor;
        J[276] += dqdci;              /* dwdot[H]/d[CH4] */
        J[287] += dqdci;              /* dwdot[CO]/d[CH4] */
        J[289] -= dqdci;              /* dwdot[HCO]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[17][3] - 1)*q_nocor - k_r*sc[1];
        J[301] += dqdci;              /* dwdot[H]/d[CO] */
        J[312] += dqdci;              /* dwdot[CO]/d[CO] */
        J[314] -= dqdci;              /* dwdot[HCO]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[17][4] - 1)*q_nocor;
        J[326] += dqdci;              /* dwdot[H]/d[CO2] */
        J[337] += dqdci;              /* dwdot[CO]/d[CO2] */
        J[339] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        /* d()/d[HCO] */
        dqdci =  + k_f;
        J[351] += dqdci;              /* dwdot[H]/d[HCO] */
        J[362] += dqdci;              /* dwdot[CO]/d[HCO] */
        J[364] -= dqdci;              /* dwdot[HCO]/d[HCO] */
        /* d()/d[C2H6] */
        dqdci = (TB[17][5] - 1)*q_nocor;
        J[526] += dqdci;              /* dwdot[H]/d[C2H6] */
        J[537] += dqdci;              /* dwdot[CO]/d[C2H6] */
        J[539] -= dqdci;              /* dwdot[HCO]/d[C2H6] */
    }
    else {
        dqdc[0] = TB[17][0];
        dqdc[1] = dcdc_fac - k_r*sc[12];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[17][1];
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[17][2];
        dqdc[12] = TB[17][3] - k_r*sc[1];
        dqdc[13] = TB[17][4];
        dqdc[14] = dcdc_fac + k_f;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = TB[17][5];
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        for (int k=0; k<24; k++) {
            J[25*k+1] += dqdc[k];
            J[25*k+12] += dqdc[k];
            J[25*k+14] -= dqdc[k];
        }
    }
    J[601] += dqdT; /* dwdot[H]/dT */
    J[612] += dqdT; /* dwdot[CO]/dT */
    J[614] -= dqdT; /* dwdot[HCO]/dT */

    /*reaction 19: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
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
    J[25] -= dqdci;               /* dwdot[H2]/d[H] */
    J[26] += dqdci;               /* dwdot[H]/d[H] */
    J[27] -= dqdci;               /* dwdot[O]/d[H] */
    J[29] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[0];
    J[50] -= dqdci;               /* dwdot[H2]/d[O] */
    J[51] += dqdci;               /* dwdot[H]/d[O] */
    J[52] -= dqdci;               /* dwdot[O]/d[O] */
    J[54] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[100] -= dqdci;              /* dwdot[H2]/d[OH] */
    J[101] += dqdci;              /* dwdot[H]/d[OH] */
    J[102] -= dqdci;              /* dwdot[O]/d[OH] */
    J[104] += dqdci;              /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[600] -= dqdT;               /* dwdot[H2]/dT */
    J[601] += dqdT;               /* dwdot[H]/dT */
    J[602] -= dqdT;               /* dwdot[O]/dT */
    J[604] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 20: O + HO2 <=> OH + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
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
    J[52] -= dqdci;               /* dwdot[O]/d[O] */
    J[53] += dqdci;               /* dwdot[O2]/d[O] */
    J[54] += dqdci;               /* dwdot[OH]/d[O] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[4];
    J[77] -= dqdci;               /* dwdot[O]/d[O2] */
    J[78] += dqdci;               /* dwdot[O2]/d[O2] */
    J[79] += dqdci;               /* dwdot[OH]/d[O2] */
    J[81] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[3];
    J[102] -= dqdci;              /* dwdot[O]/d[OH] */
    J[103] += dqdci;              /* dwdot[O2]/d[OH] */
    J[104] += dqdci;              /* dwdot[OH]/d[OH] */
    J[106] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[152] -= dqdci;              /* dwdot[O]/d[HO2] */
    J[153] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[154] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[156] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[602] -= dqdT;               /* dwdot[O]/dT */
    J[603] += dqdT;               /* dwdot[O2]/dT */
    J[604] += dqdT;               /* dwdot[OH]/dT */
    J[606] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 21: O + CH2 <=> H + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[8];
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[14];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[8] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[8]) + (h_RT[1] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[8] -= q; /* CH2 */
    wdot[14] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[26] += dqdci;               /* dwdot[H]/d[H] */
    J[27] -= dqdci;               /* dwdot[O]/d[H] */
    J[33] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[39] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[8];
    J[51] += dqdci;               /* dwdot[H]/d[O] */
    J[52] -= dqdci;               /* dwdot[O]/d[O] */
    J[58] -= dqdci;               /* dwdot[CH2]/d[O] */
    J[64] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[2];
    J[201] += dqdci;              /* dwdot[H]/d[CH2] */
    J[202] -= dqdci;              /* dwdot[O]/d[CH2] */
    J[208] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[214] += dqdci;              /* dwdot[HCO]/d[CH2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[1];
    J[351] += dqdci;              /* dwdot[H]/d[HCO] */
    J[352] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[358] -= dqdci;              /* dwdot[CH2]/d[HCO] */
    J[364] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[601] += dqdT;               /* dwdot[H]/dT */
    J[602] -= dqdT;               /* dwdot[O]/dT */
    J[608] -= dqdT;               /* dwdot[CH2]/dT */
    J[614] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 22: O + CH2(S) <=> H + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[9];
    k_f = prefactor_units[21] * fwd_A[21]
                * exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
    dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[14];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[9] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[9]) + (h_RT[1] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[9] -= q; /* CH2(S) */
    wdot[14] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[26] += dqdci;               /* dwdot[H]/d[H] */
    J[27] -= dqdci;               /* dwdot[O]/d[H] */
    J[34] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[39] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[9];
    J[51] += dqdci;               /* dwdot[H]/d[O] */
    J[52] -= dqdci;               /* dwdot[O]/d[O] */
    J[59] -= dqdci;               /* dwdot[CH2(S)]/d[O] */
    J[64] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[2];
    J[226] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[227] -= dqdci;              /* dwdot[O]/d[CH2(S)] */
    J[234] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[239] += dqdci;              /* dwdot[HCO]/d[CH2(S)] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[1];
    J[351] += dqdci;              /* dwdot[H]/d[HCO] */
    J[352] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[359] -= dqdci;              /* dwdot[CH2(S)]/d[HCO] */
    J[364] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[601] += dqdT;               /* dwdot[H]/dT */
    J[602] -= dqdT;               /* dwdot[O]/dT */
    J[609] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[614] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 23: O + CH3 <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[10];
    k_f = prefactor_units[22] * fwd_A[22]
                * exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
    dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[15];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[10] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[10]) + (h_RT[1] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[10] -= q; /* CH3 */
    wdot[15] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[15];
    J[26] += dqdci;               /* dwdot[H]/d[H] */
    J[27] -= dqdci;               /* dwdot[O]/d[H] */
    J[35] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[40] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[10];
    J[51] += dqdci;               /* dwdot[H]/d[O] */
    J[52] -= dqdci;               /* dwdot[O]/d[O] */
    J[60] -= dqdci;               /* dwdot[CH3]/d[O] */
    J[65] += dqdci;               /* dwdot[CH2O]/d[O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[2];
    J[251] += dqdci;              /* dwdot[H]/d[CH3] */
    J[252] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[260] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[265] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[376] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[377] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[385] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[390] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[601] += dqdT;               /* dwdot[H]/dT */
    J[602] -= dqdT;               /* dwdot[O]/dT */
    J[610] -= dqdT;               /* dwdot[CH3]/dT */
    J[615] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 24: O + CH4 <=> OH + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[11];
    k_f = prefactor_units[23] * fwd_A[23]
                * exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
    dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[10];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[10] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[11]) + (h_RT[4] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[10] += q; /* CH3 */
    wdot[11] -= q; /* CH4 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[11];
    J[52] -= dqdci;               /* dwdot[O]/d[O] */
    J[54] += dqdci;               /* dwdot[OH]/d[O] */
    J[60] += dqdci;               /* dwdot[CH3]/d[O] */
    J[61] -= dqdci;               /* dwdot[CH4]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[10];
    J[102] -= dqdci;              /* dwdot[O]/d[OH] */
    J[104] += dqdci;              /* dwdot[OH]/d[OH] */
    J[110] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[111] -= dqdci;              /* dwdot[CH4]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[4];
    J[252] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[254] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[260] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[261] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[2];
    J[277] -= dqdci;              /* dwdot[O]/d[CH4] */
    J[279] += dqdci;              /* dwdot[OH]/d[CH4] */
    J[285] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[286] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[602] -= dqdT;               /* dwdot[O]/dT */
    J[604] += dqdT;               /* dwdot[OH]/dT */
    J[610] += dqdT;               /* dwdot[CH3]/dT */
    J[611] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 25: O + HCO <=> OH + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[14];
    k_f = prefactor_units[24] * fwd_A[24]
                * exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
    dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[12];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[12] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[14]) + (h_RT[4] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[12] += q; /* CO */
    wdot[14] -= q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[14];
    J[52] -= dqdci;               /* dwdot[O]/d[O] */
    J[54] += dqdci;               /* dwdot[OH]/d[O] */
    J[62] += dqdci;               /* dwdot[CO]/d[O] */
    J[64] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[12];
    J[102] -= dqdci;              /* dwdot[O]/d[OH] */
    J[104] += dqdci;              /* dwdot[OH]/d[OH] */
    J[112] += dqdci;              /* dwdot[CO]/d[OH] */
    J[114] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[4];
    J[302] -= dqdci;              /* dwdot[O]/d[CO] */
    J[304] += dqdci;              /* dwdot[OH]/d[CO] */
    J[312] += dqdci;              /* dwdot[CO]/d[CO] */
    J[314] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[2];
    J[352] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[354] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[362] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[364] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[602] -= dqdT;               /* dwdot[O]/dT */
    J[604] += dqdT;               /* dwdot[OH]/dT */
    J[612] += dqdT;               /* dwdot[CO]/dT */
    J[614] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 26: O + HCO <=> H + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[14];
    k_f = prefactor_units[25] * fwd_A[25]
                * exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
    dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[13];
    Kc = exp(-g_RT[1] + g_RT[2] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[14]) + (h_RT[1] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[13] += q; /* CO2 */
    wdot[14] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[26] += dqdci;               /* dwdot[H]/d[H] */
    J[27] -= dqdci;               /* dwdot[O]/d[H] */
    J[38] += dqdci;               /* dwdot[CO2]/d[H] */
    J[39] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[14];
    J[51] += dqdci;               /* dwdot[H]/d[O] */
    J[52] -= dqdci;               /* dwdot[O]/d[O] */
    J[63] += dqdci;               /* dwdot[CO2]/d[O] */
    J[64] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[326] += dqdci;              /* dwdot[H]/d[CO2] */
    J[327] -= dqdci;              /* dwdot[O]/d[CO2] */
    J[338] += dqdci;              /* dwdot[CO2]/d[CO2] */
    J[339] -= dqdci;              /* dwdot[HCO]/d[CO2] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[2];
    J[351] += dqdci;              /* dwdot[H]/d[HCO] */
    J[352] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[363] += dqdci;              /* dwdot[CO2]/d[HCO] */
    J[364] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[601] += dqdT;               /* dwdot[H]/dT */
    J[602] -= dqdT;               /* dwdot[O]/dT */
    J[613] += dqdT;               /* dwdot[CO2]/dT */
    J[614] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 27: O + CH2O <=> OH + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[15];
    k_f = prefactor_units[26] * fwd_A[26]
                * exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
    dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[14];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[14] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[15]) + (h_RT[4] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[14] += q; /* HCO */
    wdot[15] -= q; /* CH2O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[15];
    J[52] -= dqdci;               /* dwdot[O]/d[O] */
    J[54] += dqdci;               /* dwdot[OH]/d[O] */
    J[64] += dqdci;               /* dwdot[HCO]/d[O] */
    J[65] -= dqdci;               /* dwdot[CH2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[14];
    J[102] -= dqdci;              /* dwdot[O]/d[OH] */
    J[104] += dqdci;              /* dwdot[OH]/d[OH] */
    J[114] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[115] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[4];
    J[352] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[354] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[364] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[365] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[2];
    J[377] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[379] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[389] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[390] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[602] -= dqdT;               /* dwdot[O]/dT */
    J[604] += dqdT;               /* dwdot[OH]/dT */
    J[614] += dqdT;               /* dwdot[HCO]/dT */
    J[615] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 28: O + C2H2 <=> CH2(S) + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[17];
    k_f = prefactor_units[27] * fwd_A[27]
                * exp(fwd_beta[27] * tc[0] - activation_units[27] * fwd_Ea[27] * invT);
    dlnkfdT = fwd_beta[27] * invT + activation_units[27] * fwd_Ea[27] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[12];
    Kc = exp(g_RT[2] - g_RT[9] - g_RT[12] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[17]) + (h_RT[9] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[9] += q; /* CH2(S) */
    wdot[12] += q; /* CO */
    wdot[17] -= q; /* C2H2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[17];
    J[52] -= dqdci;               /* dwdot[O]/d[O] */
    J[59] += dqdci;               /* dwdot[CH2(S)]/d[O] */
    J[62] += dqdci;               /* dwdot[CO]/d[O] */
    J[67] -= dqdci;               /* dwdot[C2H2]/d[O] */
    /* d()/d[CH2(S)] */
    dqdci =  - k_r*sc[12];
    J[227] -= dqdci;              /* dwdot[O]/d[CH2(S)] */
    J[234] += dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[237] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    J[242] -= dqdci;              /* dwdot[C2H2]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[9];
    J[302] -= dqdci;              /* dwdot[O]/d[CO] */
    J[309] += dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[312] += dqdci;              /* dwdot[CO]/d[CO] */
    J[317] -= dqdci;              /* dwdot[C2H2]/d[CO] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[2];
    J[427] -= dqdci;              /* dwdot[O]/d[C2H2] */
    J[434] += dqdci;              /* dwdot[CH2(S)]/d[C2H2] */
    J[437] += dqdci;              /* dwdot[CO]/d[C2H2] */
    J[442] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    /* d()/dT */
    J[602] -= dqdT;               /* dwdot[O]/dT */
    J[609] += dqdT;               /* dwdot[CH2(S)]/dT */
    J[612] += dqdT;               /* dwdot[CO]/dT */
    J[617] -= dqdT;               /* dwdot[C2H2]/dT */

    /*reaction 29: O + C2H2 <=> CO + CH2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[17];
    k_f = prefactor_units[28] * fwd_A[28]
                * exp(fwd_beta[28] * tc[0] - activation_units[28] * fwd_Ea[28] * invT);
    dlnkfdT = fwd_beta[28] * invT + activation_units[28] * fwd_Ea[28] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[12];
    Kc = exp(g_RT[2] - g_RT[8] - g_RT[12] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[17]) + (h_RT[8] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[8] += q; /* CH2 */
    wdot[12] += q; /* CO */
    wdot[17] -= q; /* C2H2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[17];
    J[52] -= dqdci;               /* dwdot[O]/d[O] */
    J[58] += dqdci;               /* dwdot[CH2]/d[O] */
    J[62] += dqdci;               /* dwdot[CO]/d[O] */
    J[67] -= dqdci;               /* dwdot[C2H2]/d[O] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[12];
    J[202] -= dqdci;              /* dwdot[O]/d[CH2] */
    J[208] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[212] += dqdci;              /* dwdot[CO]/d[CH2] */
    J[217] -= dqdci;              /* dwdot[C2H2]/d[CH2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[8];
    J[302] -= dqdci;              /* dwdot[O]/d[CO] */
    J[308] += dqdci;              /* dwdot[CH2]/d[CO] */
    J[312] += dqdci;              /* dwdot[CO]/d[CO] */
    J[317] -= dqdci;              /* dwdot[C2H2]/d[CO] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[2];
    J[427] -= dqdci;              /* dwdot[O]/d[C2H2] */
    J[433] += dqdci;              /* dwdot[CH2]/d[C2H2] */
    J[437] += dqdci;              /* dwdot[CO]/d[C2H2] */
    J[442] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    /* d()/dT */
    J[602] -= dqdT;               /* dwdot[O]/dT */
    J[608] += dqdT;               /* dwdot[CH2]/dT */
    J[612] += dqdT;               /* dwdot[CO]/dT */
    J[617] -= dqdT;               /* dwdot[C2H2]/dT */

    /*reaction 30: O + C2H4 <=> CH3 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[19];
    k_f = prefactor_units[29] * fwd_A[29]
                * exp(fwd_beta[29] * tc[0] - activation_units[29] * fwd_Ea[29] * invT);
    dlnkfdT = fwd_beta[29] * invT + activation_units[29] * fwd_Ea[29] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[14];
    Kc = exp(g_RT[2] - g_RT[10] - g_RT[14] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[19]) + (h_RT[10] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[10] += q; /* CH3 */
    wdot[14] += q; /* HCO */
    wdot[19] -= q; /* C2H4 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[19];
    J[52] -= dqdci;               /* dwdot[O]/d[O] */
    J[60] += dqdci;               /* dwdot[CH3]/d[O] */
    J[64] += dqdci;               /* dwdot[HCO]/d[O] */
    J[69] -= dqdci;               /* dwdot[C2H4]/d[O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[14];
    J[252] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[260] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[264] += dqdci;              /* dwdot[HCO]/d[CH3] */
    J[269] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[10];
    J[352] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[360] += dqdci;              /* dwdot[CH3]/d[HCO] */
    J[364] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[369] -= dqdci;              /* dwdot[C2H4]/d[HCO] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[2];
    J[477] -= dqdci;              /* dwdot[O]/d[C2H4] */
    J[485] += dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[489] += dqdci;              /* dwdot[HCO]/d[C2H4] */
    J[494] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[602] -= dqdT;               /* dwdot[O]/dT */
    J[610] += dqdT;               /* dwdot[CH3]/dT */
    J[614] += dqdT;               /* dwdot[HCO]/dT */
    J[619] -= dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 31: O + C2H5 <=> CH3 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[20];
    k_f = prefactor_units[30] * fwd_A[30]
                * exp(fwd_beta[30] * tc[0] - activation_units[30] * fwd_Ea[30] * invT);
    dlnkfdT = fwd_beta[30] * invT + activation_units[30] * fwd_Ea[30] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[15];
    Kc = exp(g_RT[2] - g_RT[10] - g_RT[15] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[20]) + (h_RT[10] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[10] += q; /* CH3 */
    wdot[15] += q; /* CH2O */
    wdot[20] -= q; /* C2H5 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[20];
    J[52] -= dqdci;               /* dwdot[O]/d[O] */
    J[60] += dqdci;               /* dwdot[CH3]/d[O] */
    J[65] += dqdci;               /* dwdot[CH2O]/d[O] */
    J[70] -= dqdci;               /* dwdot[C2H5]/d[O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[15];
    J[252] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[260] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[265] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[270] -= dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[10];
    J[377] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[385] += dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[390] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[395] -= dqdci;              /* dwdot[C2H5]/d[CH2O] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[2];
    J[502] -= dqdci;              /* dwdot[O]/d[C2H5] */
    J[510] += dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[515] += dqdci;              /* dwdot[CH2O]/d[C2H5] */
    J[520] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[602] -= dqdT;               /* dwdot[O]/dT */
    J[610] += dqdT;               /* dwdot[CH3]/dT */
    J[615] += dqdT;               /* dwdot[CH2O]/dT */
    J[620] -= dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 32: O + C2H6 <=> OH + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[21];
    k_f = prefactor_units[31] * fwd_A[31]
                * exp(fwd_beta[31] * tc[0] - activation_units[31] * fwd_Ea[31] * invT);
    dlnkfdT = fwd_beta[31] * invT + activation_units[31] * fwd_Ea[31] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[20];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[20] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[21]) + (h_RT[4] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[20] += q; /* C2H5 */
    wdot[21] -= q; /* C2H6 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[21];
    J[52] -= dqdci;               /* dwdot[O]/d[O] */
    J[54] += dqdci;               /* dwdot[OH]/d[O] */
    J[70] += dqdci;               /* dwdot[C2H5]/d[O] */
    J[71] -= dqdci;               /* dwdot[C2H6]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[20];
    J[102] -= dqdci;              /* dwdot[O]/d[OH] */
    J[104] += dqdci;              /* dwdot[OH]/d[OH] */
    J[120] += dqdci;              /* dwdot[C2H5]/d[OH] */
    J[121] -= dqdci;              /* dwdot[C2H6]/d[OH] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[4];
    J[502] -= dqdci;              /* dwdot[O]/d[C2H5] */
    J[504] += dqdci;              /* dwdot[OH]/d[C2H5] */
    J[520] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[521] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[2];
    J[527] -= dqdci;              /* dwdot[O]/d[C2H6] */
    J[529] += dqdci;              /* dwdot[OH]/d[C2H6] */
    J[545] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[546] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[602] -= dqdT;               /* dwdot[O]/dT */
    J[604] += dqdT;               /* dwdot[OH]/dT */
    J[620] += dqdT;               /* dwdot[C2H5]/dT */
    J[621] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 33: O2 + CO <=> O + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[12];
    k_f = prefactor_units[32] * fwd_A[32]
                * exp(fwd_beta[32] * tc[0] - activation_units[32] * fwd_Ea[32] * invT);
    dlnkfdT = fwd_beta[32] * invT + activation_units[32] * fwd_Ea[32] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[13];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[12] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[12]) + (h_RT[2] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[12] -= q; /* CO */
    wdot[13] += q; /* CO2 */
    /* d()/d[O] */
    dqdci =  - k_r*sc[13];
    J[52] += dqdci;               /* dwdot[O]/d[O] */
    J[53] -= dqdci;               /* dwdot[O2]/d[O] */
    J[62] -= dqdci;               /* dwdot[CO]/d[O] */
    J[63] += dqdci;               /* dwdot[CO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[12];
    J[77] += dqdci;               /* dwdot[O]/d[O2] */
    J[78] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[87] -= dqdci;               /* dwdot[CO]/d[O2] */
    J[88] += dqdci;               /* dwdot[CO2]/d[O2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[3];
    J[302] += dqdci;              /* dwdot[O]/d[CO] */
    J[303] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[312] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[313] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[2];
    J[327] += dqdci;              /* dwdot[O]/d[CO2] */
    J[328] -= dqdci;              /* dwdot[O2]/d[CO2] */
    J[337] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[338] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[602] += dqdT;               /* dwdot[O]/dT */
    J[603] -= dqdT;               /* dwdot[O2]/dT */
    J[612] -= dqdT;               /* dwdot[CO]/dT */
    J[613] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 34: O2 + CH2O <=> HO2 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[15];
    k_f = prefactor_units[33] * fwd_A[33]
                * exp(fwd_beta[33] * tc[0] - activation_units[33] * fwd_Ea[33] * invT);
    dlnkfdT = fwd_beta[33] * invT + activation_units[33] * fwd_Ea[33] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[14];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[14] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[15]) + (h_RT[6] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[14] += q; /* HCO */
    wdot[15] -= q; /* CH2O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[15];
    J[78] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[81] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[89] += dqdci;               /* dwdot[HCO]/d[O2] */
    J[90] -= dqdci;               /* dwdot[CH2O]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[14];
    J[153] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[156] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[164] += dqdci;              /* dwdot[HCO]/d[HO2] */
    J[165] -= dqdci;              /* dwdot[CH2O]/d[HO2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[6];
    J[353] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[356] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[364] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[365] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[3];
    J[378] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[381] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[389] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[390] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[603] -= dqdT;               /* dwdot[O2]/dT */
    J[606] += dqdT;               /* dwdot[HO2]/dT */
    J[614] += dqdT;               /* dwdot[HCO]/dT */
    J[615] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 35: H + 2 O2 <=> HO2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[3];
    k_f = prefactor_units[34] * fwd_A[34]
                * exp(fwd_beta[34] * tc[0] - activation_units[34] * fwd_Ea[34] * invT);
    dlnkfdT = fwd_beta[34] * invT + activation_units[34] * fwd_Ea[34] * invT2;
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
    J[26] -= dqdci;               /* dwdot[H]/d[H] */
    J[28] -= dqdci;               /* dwdot[O2]/d[H] */
    J[31] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*2*sc[3] - k_r*sc[6];
    J[76] -= dqdci;               /* dwdot[H]/d[O2] */
    J[78] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[81] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[3];
    J[151] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[153] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[156] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[601] -= dqdT;               /* dwdot[H]/dT */
    J[603] -= dqdT;               /* dwdot[O2]/dT */
    J[606] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 36: H + O2 + H2O <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = prefactor_units[35] * fwd_A[35]
                * exp(fwd_beta[35] * tc[0] - activation_units[35] * fwd_Ea[35] * invT);
    dlnkfdT = fwd_beta[35] * invT + activation_units[35] * fwd_Ea[35] * invT2;
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
    J[26] -= dqdci;               /* dwdot[H]/d[H] */
    J[28] -= dqdci;               /* dwdot[O2]/d[H] */
    J[31] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[5];
    J[76] -= dqdci;               /* dwdot[H]/d[O2] */
    J[78] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[81] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[126] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[128] -= dqdci;              /* dwdot[O2]/d[H2O] */
    J[131] += dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[151] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[153] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[156] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[601] -= dqdT;               /* dwdot[H]/dT */
    J[603] -= dqdT;               /* dwdot[O2]/dT */
    J[606] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 37: H + O2 + N2 <=> HO2 + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[22];
    k_f = prefactor_units[36] * fwd_A[36]
                * exp(fwd_beta[36] * tc[0] - activation_units[36] * fwd_Ea[36] * invT);
    dlnkfdT = fwd_beta[36] * invT + activation_units[36] * fwd_Ea[36] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[22];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[22] - g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[22]) + (h_RT[6] + h_RT[22]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[22];
    J[26] -= dqdci;               /* dwdot[H]/d[H] */
    J[28] -= dqdci;               /* dwdot[O2]/d[H] */
    J[31] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[22];
    J[76] -= dqdci;               /* dwdot[H]/d[O2] */
    J[78] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[81] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[22];
    J[151] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[153] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[156] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[551] -= dqdci;              /* dwdot[H]/d[N2] */
    J[553] -= dqdci;              /* dwdot[O2]/d[N2] */
    J[556] += dqdci;              /* dwdot[HO2]/d[N2] */
    /* d()/dT */
    J[601] -= dqdT;               /* dwdot[H]/dT */
    J[603] -= dqdT;               /* dwdot[O2]/dT */
    J[606] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 38: H + O2 + AR <=> HO2 + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[23];
    k_f = prefactor_units[37] * fwd_A[37]
                * exp(fwd_beta[37] * tc[0] - activation_units[37] * fwd_Ea[37] * invT);
    dlnkfdT = fwd_beta[37] * invT + activation_units[37] * fwd_Ea[37] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[23];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[23] - g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[23]) + (h_RT[6] + h_RT[23]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[23];
    J[26] -= dqdci;               /* dwdot[H]/d[H] */
    J[28] -= dqdci;               /* dwdot[O2]/d[H] */
    J[31] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[23];
    J[76] -= dqdci;               /* dwdot[H]/d[O2] */
    J[78] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[81] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[23];
    J[151] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[153] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[156] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/d[AR] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[576] -= dqdci;              /* dwdot[H]/d[AR] */
    J[578] -= dqdci;              /* dwdot[O2]/d[AR] */
    J[581] += dqdci;              /* dwdot[HO2]/d[AR] */
    /* d()/dT */
    J[601] -= dqdT;               /* dwdot[H]/dT */
    J[603] -= dqdT;               /* dwdot[O2]/dT */
    J[606] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 39: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[38] * fwd_A[38]
                * exp(fwd_beta[38] * tc[0] - activation_units[38] * fwd_Ea[38] * invT);
    dlnkfdT = fwd_beta[38] * invT + activation_units[38] * fwd_Ea[38] * invT2;
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
    J[26] -= dqdci;               /* dwdot[H]/d[H] */
    J[27] += dqdci;               /* dwdot[O]/d[H] */
    J[28] -= dqdci;               /* dwdot[O2]/d[H] */
    J[29] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[4];
    J[51] -= dqdci;               /* dwdot[H]/d[O] */
    J[52] += dqdci;               /* dwdot[O]/d[O] */
    J[53] -= dqdci;               /* dwdot[O2]/d[O] */
    J[54] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1];
    J[76] -= dqdci;               /* dwdot[H]/d[O2] */
    J[77] += dqdci;               /* dwdot[O]/d[O2] */
    J[78] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[79] += dqdci;               /* dwdot[OH]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[101] -= dqdci;              /* dwdot[H]/d[OH] */
    J[102] += dqdci;              /* dwdot[O]/d[OH] */
    J[103] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[104] += dqdci;              /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[601] -= dqdT;               /* dwdot[H]/dT */
    J[602] += dqdT;               /* dwdot[O]/dT */
    J[603] -= dqdT;               /* dwdot[O2]/dT */
    J[604] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[1]*sc[1];
    k_f = prefactor_units[39] * fwd_A[39]
                * exp(fwd_beta[39] * tc[0] - activation_units[39] * fwd_Ea[39] * invT);
    dlnkfdT = fwd_beta[39] * invT + activation_units[39] * fwd_Ea[39] * invT2;
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
    J[25] += dqdci;               /* dwdot[H2]/d[H] */
    J[26] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/dT */
    J[600] += dqdT;               /* dwdot[H2]/dT */
    J[601] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[1]*sc[5];
    k_f = prefactor_units[40] * fwd_A[40]
                * exp(fwd_beta[40] * tc[0] - activation_units[40] * fwd_Ea[40] * invT);
    dlnkfdT = fwd_beta[40] * invT + activation_units[40] * fwd_Ea[40] * invT2;
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
    J[25] += dqdci;               /* dwdot[H2]/d[H] */
    J[26] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[1]*sc[1] - k_r*sc[0];
    J[125] += dqdci;              /* dwdot[H2]/d[H2O] */
    J[126] += -2 * dqdci;         /* dwdot[H]/d[H2O] */
    /* d()/dT */
    J[600] += dqdT;               /* dwdot[H2]/dT */
    J[601] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[1]*sc[13];
    k_f = prefactor_units[41] * fwd_A[41]
                * exp(fwd_beta[41] * tc[0] - activation_units[41] * fwd_Ea[41] * invT);
    dlnkfdT = fwd_beta[41] * invT + activation_units[41] * fwd_Ea[41] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[13];
    Kc = refCinv * exp(-g_RT[0] + 2*g_RT[1] + g_RT[13] - g_RT[13]);
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
    J[25] += dqdci;               /* dwdot[H2]/d[H] */
    J[26] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[1]*sc[1] - k_r*sc[0];
    J[325] += dqdci;              /* dwdot[H2]/d[CO2] */
    J[326] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
    /* d()/dT */
    J[600] += dqdT;               /* dwdot[H2]/dT */
    J[601] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 43: H + HO2 <=> O2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = prefactor_units[42] * fwd_A[42]
                * exp(fwd_beta[42] * tc[0] - activation_units[42] * fwd_Ea[42] * invT);
    dlnkfdT = fwd_beta[42] * invT + activation_units[42] * fwd_Ea[42] * invT2;
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
    J[25] += dqdci;               /* dwdot[H2]/d[H] */
    J[26] -= dqdci;               /* dwdot[H]/d[H] */
    J[28] += dqdci;               /* dwdot[O2]/d[H] */
    J[31] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[0];
    J[75] += dqdci;               /* dwdot[H2]/d[O2] */
    J[76] -= dqdci;               /* dwdot[H]/d[O2] */
    J[78] += dqdci;               /* dwdot[O2]/d[O2] */
    J[81] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[150] += dqdci;              /* dwdot[H2]/d[HO2] */
    J[151] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[153] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[156] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[600] += dqdT;               /* dwdot[H2]/dT */
    J[601] -= dqdT;               /* dwdot[H]/dT */
    J[603] += dqdT;               /* dwdot[O2]/dT */
    J[606] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 44: H + HO2 <=> 2 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = prefactor_units[43] * fwd_A[43]
                * exp(fwd_beta[43] * tc[0] - activation_units[43] * fwd_Ea[43] * invT);
    dlnkfdT = fwd_beta[43] * invT + activation_units[43] * fwd_Ea[43] * invT2;
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
    J[26] -= dqdci;               /* dwdot[H]/d[H] */
    J[29] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[31] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2*sc[4];
    J[101] -= dqdci;              /* dwdot[H]/d[OH] */
    J[104] += 2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[106] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[151] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[154] += 2 * dqdci;          /* dwdot[OH]/d[HO2] */
    J[156] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[601] -= dqdT;               /* dwdot[H]/dT */
    J[604] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[606] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 45: H + H2O2 <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = prefactor_units[44] * fwd_A[44]
                * exp(fwd_beta[44] * tc[0] - activation_units[44] * fwd_Ea[44] * invT);
    dlnkfdT = fwd_beta[44] * invT + activation_units[44] * fwd_Ea[44] * invT2;
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
    J[25] += dqdci;               /* dwdot[H2]/d[H] */
    J[26] -= dqdci;               /* dwdot[H]/d[H] */
    J[31] += dqdci;               /* dwdot[HO2]/d[H] */
    J[32] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[0];
    J[150] += dqdci;              /* dwdot[H2]/d[HO2] */
    J[151] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[156] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[157] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[175] += dqdci;              /* dwdot[H2]/d[H2O2] */
    J[176] -= dqdci;              /* dwdot[H]/d[H2O2] */
    J[181] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[182] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[600] += dqdT;               /* dwdot[H2]/dT */
    J[601] -= dqdT;               /* dwdot[H]/dT */
    J[606] += dqdT;               /* dwdot[HO2]/dT */
    J[607] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 46: H + CH4 <=> CH3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[11];
    k_f = prefactor_units[45] * fwd_A[45]
                * exp(fwd_beta[45] * tc[0] - activation_units[45] * fwd_Ea[45] * invT);
    dlnkfdT = fwd_beta[45] * invT + activation_units[45] * fwd_Ea[45] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[10];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[10] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[11]) + (h_RT[0] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[10] += q; /* CH3 */
    wdot[11] -= q; /* CH4 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[10];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[10] += dqdci;               /* dwdot[CH3]/d[H2] */
    J[11] -= dqdci;               /* dwdot[CH4]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[11];
    J[25] += dqdci;               /* dwdot[H2]/d[H] */
    J[26] -= dqdci;               /* dwdot[H]/d[H] */
    J[35] += dqdci;               /* dwdot[CH3]/d[H] */
    J[36] -= dqdci;               /* dwdot[CH4]/d[H] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[0];
    J[250] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[251] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[260] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[261] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[1];
    J[275] += dqdci;              /* dwdot[H2]/d[CH4] */
    J[276] -= dqdci;              /* dwdot[H]/d[CH4] */
    J[285] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[286] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[600] += dqdT;               /* dwdot[H2]/dT */
    J[601] -= dqdT;               /* dwdot[H]/dT */
    J[610] += dqdT;               /* dwdot[CH3]/dT */
    J[611] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 47: H + HCO <=> H2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[14];
    k_f = prefactor_units[46] * fwd_A[46]
                * exp(fwd_beta[46] * tc[0] - activation_units[46] * fwd_Ea[46] * invT);
    dlnkfdT = fwd_beta[46] * invT + activation_units[46] * fwd_Ea[46] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[12];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[12] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[14]) + (h_RT[0] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[12] += q; /* CO */
    wdot[14] -= q; /* HCO */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[12];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[12] += dqdci;               /* dwdot[CO]/d[H2] */
    J[14] -= dqdci;               /* dwdot[HCO]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[14];
    J[25] += dqdci;               /* dwdot[H2]/d[H] */
    J[26] -= dqdci;               /* dwdot[H]/d[H] */
    J[37] += dqdci;               /* dwdot[CO]/d[H] */
    J[39] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[0];
    J[300] += dqdci;              /* dwdot[H2]/d[CO] */
    J[301] -= dqdci;              /* dwdot[H]/d[CO] */
    J[312] += dqdci;              /* dwdot[CO]/d[CO] */
    J[314] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[1];
    J[350] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[351] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[362] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[364] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[600] += dqdT;               /* dwdot[H2]/dT */
    J[601] -= dqdT;               /* dwdot[H]/dT */
    J[612] += dqdT;               /* dwdot[CO]/dT */
    J[614] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 48: H + CH2O <=> HCO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[15];
    k_f = prefactor_units[47] * fwd_A[47]
                * exp(fwd_beta[47] * tc[0] - activation_units[47] * fwd_Ea[47] * invT);
    dlnkfdT = fwd_beta[47] * invT + activation_units[47] * fwd_Ea[47] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[14];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[14] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[15]) + (h_RT[0] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[14] += q; /* HCO */
    wdot[15] -= q; /* CH2O */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[14];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[14] += dqdci;               /* dwdot[HCO]/d[H2] */
    J[15] -= dqdci;               /* dwdot[CH2O]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[15];
    J[25] += dqdci;               /* dwdot[H2]/d[H] */
    J[26] -= dqdci;               /* dwdot[H]/d[H] */
    J[39] += dqdci;               /* dwdot[HCO]/d[H] */
    J[40] -= dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[0];
    J[350] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[351] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[364] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[365] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[1];
    J[375] += dqdci;              /* dwdot[H2]/d[CH2O] */
    J[376] -= dqdci;              /* dwdot[H]/d[CH2O] */
    J[389] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[390] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[600] += dqdT;               /* dwdot[H2]/dT */
    J[601] -= dqdT;               /* dwdot[H]/dT */
    J[614] += dqdT;               /* dwdot[HCO]/dT */
    J[615] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 49: H + CH3O <=> OH + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[16];
    k_f = prefactor_units[48] * fwd_A[48]
                * exp(fwd_beta[48] * tc[0] - activation_units[48] * fwd_Ea[48] * invT);
    dlnkfdT = fwd_beta[48] * invT + activation_units[48] * fwd_Ea[48] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[10];
    Kc = exp(g_RT[1] - g_RT[4] - g_RT[10] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[16]) + (h_RT[4] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += q; /* OH */
    wdot[10] += q; /* CH3 */
    wdot[16] -= q; /* CH3O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[16];
    J[26] -= dqdci;               /* dwdot[H]/d[H] */
    J[29] += dqdci;               /* dwdot[OH]/d[H] */
    J[35] += dqdci;               /* dwdot[CH3]/d[H] */
    J[41] -= dqdci;               /* dwdot[CH3O]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[10];
    J[101] -= dqdci;              /* dwdot[H]/d[OH] */
    J[104] += dqdci;              /* dwdot[OH]/d[OH] */
    J[110] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[116] -= dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[4];
    J[251] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[254] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[260] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[266] -= dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[1];
    J[401] -= dqdci;              /* dwdot[H]/d[CH3O] */
    J[404] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[410] += dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[416] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[601] -= dqdT;               /* dwdot[H]/dT */
    J[604] += dqdT;               /* dwdot[OH]/dT */
    J[610] += dqdT;               /* dwdot[CH3]/dT */
    J[616] -= dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 50: H + C2H3 <=> H2 + C2H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[18];
    k_f = prefactor_units[49] * fwd_A[49]
                * exp(fwd_beta[49] * tc[0] - activation_units[49] * fwd_Ea[49] * invT);
    dlnkfdT = fwd_beta[49] * invT + activation_units[49] * fwd_Ea[49] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[17];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[18]) + (h_RT[0] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[17] += q; /* C2H2 */
    wdot[18] -= q; /* C2H3 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[17];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[17] += dqdci;               /* dwdot[C2H2]/d[H2] */
    J[18] -= dqdci;               /* dwdot[C2H3]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[18];
    J[25] += dqdci;               /* dwdot[H2]/d[H] */
    J[26] -= dqdci;               /* dwdot[H]/d[H] */
    J[42] += dqdci;               /* dwdot[C2H2]/d[H] */
    J[43] -= dqdci;               /* dwdot[C2H3]/d[H] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[0];
    J[425] += dqdci;              /* dwdot[H2]/d[C2H2] */
    J[426] -= dqdci;              /* dwdot[H]/d[C2H2] */
    J[442] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[443] -= dqdci;              /* dwdot[C2H3]/d[C2H2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[1];
    J[450] += dqdci;              /* dwdot[H2]/d[C2H3] */
    J[451] -= dqdci;              /* dwdot[H]/d[C2H3] */
    J[467] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
    J[468] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[600] += dqdT;               /* dwdot[H2]/dT */
    J[601] -= dqdT;               /* dwdot[H]/dT */
    J[617] += dqdT;               /* dwdot[C2H2]/dT */
    J[618] -= dqdT;               /* dwdot[C2H3]/dT */

    /*reaction 51: H + C2H4 <=> C2H3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[19];
    k_f = prefactor_units[50] * fwd_A[50]
                * exp(fwd_beta[50] * tc[0] - activation_units[50] * fwd_Ea[50] * invT);
    dlnkfdT = fwd_beta[50] * invT + activation_units[50] * fwd_Ea[50] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[18];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[18] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[19]) + (h_RT[0] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[18] += q; /* C2H3 */
    wdot[19] -= q; /* C2H4 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[18];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[18] += dqdci;               /* dwdot[C2H3]/d[H2] */
    J[19] -= dqdci;               /* dwdot[C2H4]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[19];
    J[25] += dqdci;               /* dwdot[H2]/d[H] */
    J[26] -= dqdci;               /* dwdot[H]/d[H] */
    J[43] += dqdci;               /* dwdot[C2H3]/d[H] */
    J[44] -= dqdci;               /* dwdot[C2H4]/d[H] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[0];
    J[450] += dqdci;              /* dwdot[H2]/d[C2H3] */
    J[451] -= dqdci;              /* dwdot[H]/d[C2H3] */
    J[468] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[469] -= dqdci;              /* dwdot[C2H4]/d[C2H3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[1];
    J[475] += dqdci;              /* dwdot[H2]/d[C2H4] */
    J[476] -= dqdci;              /* dwdot[H]/d[C2H4] */
    J[493] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    J[494] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[600] += dqdT;               /* dwdot[H2]/dT */
    J[601] -= dqdT;               /* dwdot[H]/dT */
    J[618] += dqdT;               /* dwdot[C2H3]/dT */
    J[619] -= dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 52: H + C2H6 <=> C2H5 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[21];
    k_f = prefactor_units[51] * fwd_A[51]
                * exp(fwd_beta[51] * tc[0] - activation_units[51] * fwd_Ea[51] * invT);
    dlnkfdT = fwd_beta[51] * invT + activation_units[51] * fwd_Ea[51] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[20];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[20] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[21]) + (h_RT[0] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[20] += q; /* C2H5 */
    wdot[21] -= q; /* C2H6 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[20];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[20] += dqdci;               /* dwdot[C2H5]/d[H2] */
    J[21] -= dqdci;               /* dwdot[C2H6]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[21];
    J[25] += dqdci;               /* dwdot[H2]/d[H] */
    J[26] -= dqdci;               /* dwdot[H]/d[H] */
    J[45] += dqdci;               /* dwdot[C2H5]/d[H] */
    J[46] -= dqdci;               /* dwdot[C2H6]/d[H] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[0];
    J[500] += dqdci;              /* dwdot[H2]/d[C2H5] */
    J[501] -= dqdci;              /* dwdot[H]/d[C2H5] */
    J[520] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[521] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[1];
    J[525] += dqdci;              /* dwdot[H2]/d[C2H6] */
    J[526] -= dqdci;              /* dwdot[H]/d[C2H6] */
    J[545] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[546] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[600] += dqdT;               /* dwdot[H2]/dT */
    J[601] -= dqdT;               /* dwdot[H]/dT */
    J[620] += dqdT;               /* dwdot[C2H5]/dT */
    J[621] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 53: OH + H2 <=> H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[4];
    k_f = prefactor_units[52] * fwd_A[52]
                * exp(fwd_beta[52] * tc[0] - activation_units[52] * fwd_Ea[52] * invT);
    dlnkfdT = fwd_beta[52] * invT + activation_units[52] * fwd_Ea[52] * invT2;
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
    J[25] -= dqdci;               /* dwdot[H2]/d[H] */
    J[26] += dqdci;               /* dwdot[H]/d[H] */
    J[29] -= dqdci;               /* dwdot[OH]/d[H] */
    J[30] += dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[0];
    J[100] -= dqdci;              /* dwdot[H2]/d[OH] */
    J[101] += dqdci;              /* dwdot[H]/d[OH] */
    J[104] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[105] += dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[125] -= dqdci;              /* dwdot[H2]/d[H2O] */
    J[126] += dqdci;              /* dwdot[H]/d[H2O] */
    J[129] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[130] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[600] -= dqdT;               /* dwdot[H2]/dT */
    J[601] += dqdT;               /* dwdot[H]/dT */
    J[604] -= dqdT;               /* dwdot[OH]/dT */
    J[605] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 54: 2 OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[4];
    k_f = prefactor_units[53] * fwd_A[53]
                * exp(fwd_beta[53] * tc[0] - activation_units[53] * fwd_Ea[53] * invT);
    dlnkfdT = fwd_beta[53] * invT + activation_units[53] * fwd_Ea[53] * invT2;
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
    J[52] += dqdci;               /* dwdot[O]/d[O] */
    J[54] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2*sc[4];
    J[102] += dqdci;              /* dwdot[O]/d[OH] */
    J[104] += -2 * dqdci;         /* dwdot[OH]/d[OH] */
    J[105] += dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[127] += dqdci;              /* dwdot[O]/d[H2O] */
    J[129] += -2 * dqdci;         /* dwdot[OH]/d[H2O] */
    J[130] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[602] += dqdT;               /* dwdot[O]/dT */
    J[604] += -2 * dqdT;          /* dwdot[OH]/dT */
    J[605] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 55: OH + HO2 <=> O2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[6];
    k_f = prefactor_units[54] * fwd_A[54]
                * exp(fwd_beta[54] * tc[0] - activation_units[54] * fwd_Ea[54] * invT);
    dlnkfdT = fwd_beta[54] * invT + activation_units[54] * fwd_Ea[54] * invT2;
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
    J[78] += dqdci;               /* dwdot[O2]/d[O2] */
    J[79] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[80] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[81] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[103] += dqdci;              /* dwdot[O2]/d[OH] */
    J[104] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[105] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[106] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[128] += dqdci;              /* dwdot[O2]/d[H2O] */
    J[129] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[130] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[131] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[153] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[154] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[155] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[156] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[603] += dqdT;               /* dwdot[O2]/dT */
    J[604] -= dqdT;               /* dwdot[OH]/dT */
    J[605] += dqdT;               /* dwdot[H2O]/dT */
    J[606] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 56: OH + H2O2 <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = prefactor_units[55] * fwd_A[55]
                * exp(fwd_beta[55] * tc[0] - activation_units[55] * fwd_Ea[55] * invT);
    dlnkfdT = fwd_beta[55] * invT + activation_units[55] * fwd_Ea[55] * invT2;
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
    J[104] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[105] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[106] += dqdci;              /* dwdot[HO2]/d[OH] */
    J[107] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[129] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[130] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[131] += dqdci;              /* dwdot[HO2]/d[H2O] */
    J[132] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[154] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[155] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[156] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[157] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[4];
    J[179] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[180] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[181] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[182] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[604] -= dqdT;               /* dwdot[OH]/dT */
    J[605] += dqdT;               /* dwdot[H2O]/dT */
    J[606] += dqdT;               /* dwdot[HO2]/dT */
    J[607] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 57: OH + CH2 <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[8];
    k_f = prefactor_units[56] * fwd_A[56]
                * exp(fwd_beta[56] * tc[0] - activation_units[56] * fwd_Ea[56] * invT);
    dlnkfdT = fwd_beta[56] * invT + activation_units[56] * fwd_Ea[56] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[15];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[8] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[8]) + (h_RT[1] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[8] -= q; /* CH2 */
    wdot[15] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[15];
    J[26] += dqdci;               /* dwdot[H]/d[H] */
    J[29] -= dqdci;               /* dwdot[OH]/d[H] */
    J[33] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[40] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[101] += dqdci;              /* dwdot[H]/d[OH] */
    J[104] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[108] -= dqdci;              /* dwdot[CH2]/d[OH] */
    J[115] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[4];
    J[201] += dqdci;              /* dwdot[H]/d[CH2] */
    J[204] -= dqdci;              /* dwdot[OH]/d[CH2] */
    J[208] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[215] += dqdci;              /* dwdot[CH2O]/d[CH2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[376] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[379] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[383] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    J[390] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[601] += dqdT;               /* dwdot[H]/dT */
    J[604] -= dqdT;               /* dwdot[OH]/dT */
    J[608] -= dqdT;               /* dwdot[CH2]/dT */
    J[615] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 58: OH + CH2(S) <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[9];
    k_f = prefactor_units[57] * fwd_A[57]
                * exp(fwd_beta[57] * tc[0] - activation_units[57] * fwd_Ea[57] * invT);
    dlnkfdT = fwd_beta[57] * invT + activation_units[57] * fwd_Ea[57] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[15];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[9] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[9]) + (h_RT[1] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[9] -= q; /* CH2(S) */
    wdot[15] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[15];
    J[26] += dqdci;               /* dwdot[H]/d[H] */
    J[29] -= dqdci;               /* dwdot[OH]/d[H] */
    J[34] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[40] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[9];
    J[101] += dqdci;              /* dwdot[H]/d[OH] */
    J[104] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[109] -= dqdci;              /* dwdot[CH2(S)]/d[OH] */
    J[115] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[4];
    J[226] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[229] -= dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[234] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[240] += dqdci;              /* dwdot[CH2O]/d[CH2(S)] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[376] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[379] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[384] -= dqdci;              /* dwdot[CH2(S)]/d[CH2O] */
    J[390] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[601] += dqdT;               /* dwdot[H]/dT */
    J[604] -= dqdT;               /* dwdot[OH]/dT */
    J[609] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[615] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 59: OH + CH3 <=> CH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[10];
    k_f = prefactor_units[58] * fwd_A[58]
                * exp(fwd_beta[58] * tc[0] - activation_units[58] * fwd_Ea[58] * invT);
    dlnkfdT = fwd_beta[58] * invT + activation_units[58] * fwd_Ea[58] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[8];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[8] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[10]) + (h_RT[5] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[8] += q; /* CH2 */
    wdot[10] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[104] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[105] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[108] += dqdci;              /* dwdot[CH2]/d[OH] */
    J[110] -= dqdci;              /* dwdot[CH3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[8];
    J[129] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[130] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[133] += dqdci;              /* dwdot[CH2]/d[H2O] */
    J[135] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[5];
    J[204] -= dqdci;              /* dwdot[OH]/d[CH2] */
    J[205] += dqdci;              /* dwdot[H2O]/d[CH2] */
    J[208] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[210] -= dqdci;              /* dwdot[CH3]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[254] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[255] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[258] += dqdci;              /* dwdot[CH2]/d[CH3] */
    J[260] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[604] -= dqdT;               /* dwdot[OH]/dT */
    J[605] += dqdT;               /* dwdot[H2O]/dT */
    J[608] += dqdT;               /* dwdot[CH2]/dT */
    J[610] -= dqdT;               /* dwdot[CH3]/dT */

    /*reaction 60: OH + CH3 <=> CH2(S) + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[10];
    k_f = prefactor_units[59] * fwd_A[59]
                * exp(fwd_beta[59] * tc[0] - activation_units[59] * fwd_Ea[59] * invT);
    dlnkfdT = fwd_beta[59] * invT + activation_units[59] * fwd_Ea[59] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[10]) + (h_RT[5] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[9] += q; /* CH2(S) */
    wdot[10] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[104] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[105] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[109] += dqdci;              /* dwdot[CH2(S)]/d[OH] */
    J[110] -= dqdci;              /* dwdot[CH3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[9];
    J[129] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[130] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[134] += dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[135] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    /* d()/d[CH2(S)] */
    dqdci =  - k_r*sc[5];
    J[229] -= dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[230] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[234] += dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[235] -= dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[254] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[255] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[259] += dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[260] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[604] -= dqdT;               /* dwdot[OH]/dT */
    J[605] += dqdT;               /* dwdot[H2O]/dT */
    J[609] += dqdT;               /* dwdot[CH2(S)]/dT */
    J[610] -= dqdT;               /* dwdot[CH3]/dT */

    /*reaction 61: OH + CH4 <=> CH3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[11];
    k_f = prefactor_units[60] * fwd_A[60]
                * exp(fwd_beta[60] * tc[0] - activation_units[60] * fwd_Ea[60] * invT);
    dlnkfdT = fwd_beta[60] * invT + activation_units[60] * fwd_Ea[60] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[10];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[10] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[11]) + (h_RT[5] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[10] += q; /* CH3 */
    wdot[11] -= q; /* CH4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[104] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[105] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[110] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[111] -= dqdci;              /* dwdot[CH4]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[10];
    J[129] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[130] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[135] += dqdci;              /* dwdot[CH3]/d[H2O] */
    J[136] -= dqdci;              /* dwdot[CH4]/d[H2O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[5];
    J[254] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[255] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[260] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[261] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[4];
    J[279] -= dqdci;              /* dwdot[OH]/d[CH4] */
    J[280] += dqdci;              /* dwdot[H2O]/d[CH4] */
    J[285] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[286] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[604] -= dqdT;               /* dwdot[OH]/dT */
    J[605] += dqdT;               /* dwdot[H2O]/dT */
    J[610] += dqdT;               /* dwdot[CH3]/dT */
    J[611] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 62: OH + CO <=> H + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[12];
    k_f = prefactor_units[61] * fwd_A[61]
                * exp(fwd_beta[61] * tc[0] - activation_units[61] * fwd_Ea[61] * invT);
    dlnkfdT = fwd_beta[61] * invT + activation_units[61] * fwd_Ea[61] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[13];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[12] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[12]) + (h_RT[1] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[12] -= q; /* CO */
    wdot[13] += q; /* CO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[26] += dqdci;               /* dwdot[H]/d[H] */
    J[29] -= dqdci;               /* dwdot[OH]/d[H] */
    J[37] -= dqdci;               /* dwdot[CO]/d[H] */
    J[38] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[12];
    J[101] += dqdci;              /* dwdot[H]/d[OH] */
    J[104] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[112] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[113] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[4];
    J[301] += dqdci;              /* dwdot[H]/d[CO] */
    J[304] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[312] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[313] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[326] += dqdci;              /* dwdot[H]/d[CO2] */
    J[329] -= dqdci;              /* dwdot[OH]/d[CO2] */
    J[337] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[338] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[601] += dqdT;               /* dwdot[H]/dT */
    J[604] -= dqdT;               /* dwdot[OH]/dT */
    J[612] -= dqdT;               /* dwdot[CO]/dT */
    J[613] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 63: OH + HCO <=> H2O + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[14];
    k_f = prefactor_units[62] * fwd_A[62]
                * exp(fwd_beta[62] * tc[0] - activation_units[62] * fwd_Ea[62] * invT);
    dlnkfdT = fwd_beta[62] * invT + activation_units[62] * fwd_Ea[62] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[12];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[12] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[14]) + (h_RT[5] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[12] += q; /* CO */
    wdot[14] -= q; /* HCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[14];
    J[104] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[105] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[112] += dqdci;              /* dwdot[CO]/d[OH] */
    J[114] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[12];
    J[129] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[130] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[137] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[139] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[5];
    J[304] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[305] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[312] += dqdci;              /* dwdot[CO]/d[CO] */
    J[314] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[4];
    J[354] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[355] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[362] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[364] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[604] -= dqdT;               /* dwdot[OH]/dT */
    J[605] += dqdT;               /* dwdot[H2O]/dT */
    J[612] += dqdT;               /* dwdot[CO]/dT */
    J[614] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 64: OH + CH2O <=> HCO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[15];
    k_f = prefactor_units[63] * fwd_A[63]
                * exp(fwd_beta[63] * tc[0] - activation_units[63] * fwd_Ea[63] * invT);
    dlnkfdT = fwd_beta[63] * invT + activation_units[63] * fwd_Ea[63] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[14];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[14] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[15]) + (h_RT[5] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[14] += q; /* HCO */
    wdot[15] -= q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[15];
    J[104] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[105] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[114] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[115] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[14];
    J[129] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[130] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[139] += dqdci;              /* dwdot[HCO]/d[H2O] */
    J[140] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[5];
    J[354] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[355] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[364] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[365] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[4];
    J[379] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[380] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[389] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[390] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[604] -= dqdT;               /* dwdot[OH]/dT */
    J[605] += dqdT;               /* dwdot[H2O]/dT */
    J[614] += dqdT;               /* dwdot[HCO]/dT */
    J[615] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 65: OH + C2H2 <=> CH3 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[17];
    k_f = prefactor_units[64] * fwd_A[64]
                * exp(fwd_beta[64] * tc[0] - activation_units[64] * fwd_Ea[64] * invT);
    dlnkfdT = fwd_beta[64] * invT + activation_units[64] * fwd_Ea[64] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[12];
    Kc = exp(g_RT[4] - g_RT[10] - g_RT[12] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[17]) + (h_RT[10] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[10] += q; /* CH3 */
    wdot[12] += q; /* CO */
    wdot[17] -= q; /* C2H2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[17];
    J[104] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[110] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[112] += dqdci;              /* dwdot[CO]/d[OH] */
    J[117] -= dqdci;              /* dwdot[C2H2]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[12];
    J[254] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[260] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[262] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[267] -= dqdci;              /* dwdot[C2H2]/d[CH3] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[10];
    J[304] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[310] += dqdci;              /* dwdot[CH3]/d[CO] */
    J[312] += dqdci;              /* dwdot[CO]/d[CO] */
    J[317] -= dqdci;              /* dwdot[C2H2]/d[CO] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[4];
    J[429] -= dqdci;              /* dwdot[OH]/d[C2H2] */
    J[435] += dqdci;              /* dwdot[CH3]/d[C2H2] */
    J[437] += dqdci;              /* dwdot[CO]/d[C2H2] */
    J[442] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    /* d()/dT */
    J[604] -= dqdT;               /* dwdot[OH]/dT */
    J[610] += dqdT;               /* dwdot[CH3]/dT */
    J[612] += dqdT;               /* dwdot[CO]/dT */
    J[617] -= dqdT;               /* dwdot[C2H2]/dT */

    /*reaction 66: OH + C2H3 <=> H2O + C2H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[18];
    k_f = prefactor_units[65] * fwd_A[65]
                * exp(fwd_beta[65] * tc[0] - activation_units[65] * fwd_Ea[65] * invT);
    dlnkfdT = fwd_beta[65] * invT + activation_units[65] * fwd_Ea[65] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[17];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[18]) + (h_RT[5] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[17] += q; /* C2H2 */
    wdot[18] -= q; /* C2H3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[18];
    J[104] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[105] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[117] += dqdci;              /* dwdot[C2H2]/d[OH] */
    J[118] -= dqdci;              /* dwdot[C2H3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[17];
    J[129] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[130] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[142] += dqdci;              /* dwdot[C2H2]/d[H2O] */
    J[143] -= dqdci;              /* dwdot[C2H3]/d[H2O] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[5];
    J[429] -= dqdci;              /* dwdot[OH]/d[C2H2] */
    J[430] += dqdci;              /* dwdot[H2O]/d[C2H2] */
    J[442] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[443] -= dqdci;              /* dwdot[C2H3]/d[C2H2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[4];
    J[454] -= dqdci;              /* dwdot[OH]/d[C2H3] */
    J[455] += dqdci;              /* dwdot[H2O]/d[C2H3] */
    J[467] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
    J[468] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[604] -= dqdT;               /* dwdot[OH]/dT */
    J[605] += dqdT;               /* dwdot[H2O]/dT */
    J[617] += dqdT;               /* dwdot[C2H2]/dT */
    J[618] -= dqdT;               /* dwdot[C2H3]/dT */

    /*reaction 67: OH + C2H4 <=> C2H3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[19];
    k_f = prefactor_units[66] * fwd_A[66]
                * exp(fwd_beta[66] * tc[0] - activation_units[66] * fwd_Ea[66] * invT);
    dlnkfdT = fwd_beta[66] * invT + activation_units[66] * fwd_Ea[66] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[18];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[18] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[19]) + (h_RT[5] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[18] += q; /* C2H3 */
    wdot[19] -= q; /* C2H4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[19];
    J[104] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[105] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[118] += dqdci;              /* dwdot[C2H3]/d[OH] */
    J[119] -= dqdci;              /* dwdot[C2H4]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[18];
    J[129] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[130] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[143] += dqdci;              /* dwdot[C2H3]/d[H2O] */
    J[144] -= dqdci;              /* dwdot[C2H4]/d[H2O] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[5];
    J[454] -= dqdci;              /* dwdot[OH]/d[C2H3] */
    J[455] += dqdci;              /* dwdot[H2O]/d[C2H3] */
    J[468] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[469] -= dqdci;              /* dwdot[C2H4]/d[C2H3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[4];
    J[479] -= dqdci;              /* dwdot[OH]/d[C2H4] */
    J[480] += dqdci;              /* dwdot[H2O]/d[C2H4] */
    J[493] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    J[494] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[604] -= dqdT;               /* dwdot[OH]/dT */
    J[605] += dqdT;               /* dwdot[H2O]/dT */
    J[618] += dqdT;               /* dwdot[C2H3]/dT */
    J[619] -= dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 68: OH + C2H6 <=> C2H5 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[21];
    k_f = prefactor_units[67] * fwd_A[67]
                * exp(fwd_beta[67] * tc[0] - activation_units[67] * fwd_Ea[67] * invT);
    dlnkfdT = fwd_beta[67] * invT + activation_units[67] * fwd_Ea[67] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[20];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[20] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[21]) + (h_RT[5] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[20] += q; /* C2H5 */
    wdot[21] -= q; /* C2H6 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[21];
    J[104] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[105] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[120] += dqdci;              /* dwdot[C2H5]/d[OH] */
    J[121] -= dqdci;              /* dwdot[C2H6]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[20];
    J[129] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[130] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[145] += dqdci;              /* dwdot[C2H5]/d[H2O] */
    J[146] -= dqdci;              /* dwdot[C2H6]/d[H2O] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[5];
    J[504] -= dqdci;              /* dwdot[OH]/d[C2H5] */
    J[505] += dqdci;              /* dwdot[H2O]/d[C2H5] */
    J[520] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[521] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[4];
    J[529] -= dqdci;              /* dwdot[OH]/d[C2H6] */
    J[530] += dqdci;              /* dwdot[H2O]/d[C2H6] */
    J[545] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[546] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[604] -= dqdT;               /* dwdot[OH]/dT */
    J[605] += dqdT;               /* dwdot[H2O]/dT */
    J[620] += dqdT;               /* dwdot[C2H5]/dT */
    J[621] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 69: 2 HO2 <=> O2 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[6];
    k_f = prefactor_units[68] * fwd_A[68]
                * exp(fwd_beta[68] * tc[0] - activation_units[68] * fwd_Ea[68] * invT);
    dlnkfdT = fwd_beta[68] * invT + activation_units[68] * fwd_Ea[68] * invT2;
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
    J[78] += dqdci;               /* dwdot[O2]/d[O2] */
    J[81] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[82] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2*sc[6];
    J[153] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[156] += -2 * dqdci;         /* dwdot[HO2]/d[HO2] */
    J[157] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[3];
    J[178] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[181] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[182] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[603] += dqdT;               /* dwdot[O2]/dT */
    J[606] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[607] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 70: 2 HO2 <=> O2 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[6];
    k_f = prefactor_units[69] * fwd_A[69]
                * exp(fwd_beta[69] * tc[0] - activation_units[69] * fwd_Ea[69] * invT);
    dlnkfdT = fwd_beta[69] * invT + activation_units[69] * fwd_Ea[69] * invT2;
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
    J[78] += dqdci;               /* dwdot[O2]/d[O2] */
    J[81] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[82] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2*sc[6];
    J[153] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[156] += -2 * dqdci;         /* dwdot[HO2]/d[HO2] */
    J[157] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[3];
    J[178] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[181] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[182] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[603] += dqdT;               /* dwdot[O2]/dT */
    J[606] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[607] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 71: HO2 + CH2 <=> OH + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[8];
    k_f = prefactor_units[70] * fwd_A[70]
                * exp(fwd_beta[70] * tc[0] - activation_units[70] * fwd_Ea[70] * invT);
    dlnkfdT = fwd_beta[70] * invT + activation_units[70] * fwd_Ea[70] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[15];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[8] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[8]) + (h_RT[4] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[8] -= q; /* CH2 */
    wdot[15] += q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[15];
    J[104] += dqdci;              /* dwdot[OH]/d[OH] */
    J[106] -= dqdci;              /* dwdot[HO2]/d[OH] */
    J[108] -= dqdci;              /* dwdot[CH2]/d[OH] */
    J[115] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[8];
    J[154] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[156] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[158] -= dqdci;              /* dwdot[CH2]/d[HO2] */
    J[165] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[6];
    J[204] += dqdci;              /* dwdot[OH]/d[CH2] */
    J[206] -= dqdci;              /* dwdot[HO2]/d[CH2] */
    J[208] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[215] += dqdci;              /* dwdot[CH2O]/d[CH2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[4];
    J[379] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[381] -= dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[383] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    J[390] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[604] += dqdT;               /* dwdot[OH]/dT */
    J[606] -= dqdT;               /* dwdot[HO2]/dT */
    J[608] -= dqdT;               /* dwdot[CH2]/dT */
    J[615] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 72: HO2 + CH3 <=> O2 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[10];
    k_f = prefactor_units[71] * fwd_A[71]
                * exp(fwd_beta[71] * tc[0] - activation_units[71] * fwd_Ea[71] * invT);
    dlnkfdT = fwd_beta[71] * invT + activation_units[71] * fwd_Ea[71] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[11];
    Kc = exp(-g_RT[3] + g_RT[6] + g_RT[10] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[10]) + (h_RT[3] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    wdot[10] -= q; /* CH3 */
    wdot[11] += q; /* CH4 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[11];
    J[78] += dqdci;               /* dwdot[O2]/d[O2] */
    J[81] -= dqdci;               /* dwdot[HO2]/d[O2] */
    J[85] -= dqdci;               /* dwdot[CH3]/d[O2] */
    J[86] += dqdci;               /* dwdot[CH4]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[10];
    J[153] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[156] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[160] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[161] += dqdci;              /* dwdot[CH4]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[6];
    J[253] += dqdci;              /* dwdot[O2]/d[CH3] */
    J[256] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[260] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[261] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[3];
    J[278] += dqdci;              /* dwdot[O2]/d[CH4] */
    J[281] -= dqdci;              /* dwdot[HO2]/d[CH4] */
    J[285] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[286] += dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[603] += dqdT;               /* dwdot[O2]/dT */
    J[606] -= dqdT;               /* dwdot[HO2]/dT */
    J[610] -= dqdT;               /* dwdot[CH3]/dT */
    J[611] += dqdT;               /* dwdot[CH4]/dT */

    /*reaction 73: HO2 + CH3 <=> OH + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[10];
    k_f = prefactor_units[72] * fwd_A[72]
                * exp(fwd_beta[72] * tc[0] - activation_units[72] * fwd_Ea[72] * invT);
    dlnkfdT = fwd_beta[72] * invT + activation_units[72] * fwd_Ea[72] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[16];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[10] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[10]) + (h_RT[4] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[10] -= q; /* CH3 */
    wdot[16] += q; /* CH3O */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[16];
    J[104] += dqdci;              /* dwdot[OH]/d[OH] */
    J[106] -= dqdci;              /* dwdot[HO2]/d[OH] */
    J[110] -= dqdci;              /* dwdot[CH3]/d[OH] */
    J[116] += dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[10];
    J[154] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[156] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[160] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[166] += dqdci;              /* dwdot[CH3O]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[6];
    J[254] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[256] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[260] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[266] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[4];
    J[404] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[406] -= dqdci;              /* dwdot[HO2]/d[CH3O] */
    J[410] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[416] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[604] += dqdT;               /* dwdot[OH]/dT */
    J[606] -= dqdT;               /* dwdot[HO2]/dT */
    J[610] -= dqdT;               /* dwdot[CH3]/dT */
    J[616] += dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 74: HO2 + CO <=> OH + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[12];
    k_f = prefactor_units[73] * fwd_A[73]
                * exp(fwd_beta[73] * tc[0] - activation_units[73] * fwd_Ea[73] * invT);
    dlnkfdT = fwd_beta[73] * invT + activation_units[73] * fwd_Ea[73] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[13];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[12] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[4] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[12] -= q; /* CO */
    wdot[13] += q; /* CO2 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[13];
    J[104] += dqdci;              /* dwdot[OH]/d[OH] */
    J[106] -= dqdci;              /* dwdot[HO2]/d[OH] */
    J[112] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[113] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[12];
    J[154] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[156] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[162] -= dqdci;              /* dwdot[CO]/d[HO2] */
    J[163] += dqdci;              /* dwdot[CO2]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[6];
    J[304] += dqdci;              /* dwdot[OH]/d[CO] */
    J[306] -= dqdci;              /* dwdot[HO2]/d[CO] */
    J[312] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[313] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[4];
    J[329] += dqdci;              /* dwdot[OH]/d[CO2] */
    J[331] -= dqdci;              /* dwdot[HO2]/d[CO2] */
    J[337] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[338] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[604] += dqdT;               /* dwdot[OH]/dT */
    J[606] -= dqdT;               /* dwdot[HO2]/dT */
    J[612] -= dqdT;               /* dwdot[CO]/dT */
    J[613] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 75: HO2 + CH2O <=> HCO + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[15];
    k_f = prefactor_units[74] * fwd_A[74]
                * exp(fwd_beta[74] * tc[0] - activation_units[74] * fwd_Ea[74] * invT);
    dlnkfdT = fwd_beta[74] * invT + activation_units[74] * fwd_Ea[74] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[14];
    Kc = exp(g_RT[6] - g_RT[7] - g_RT[14] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[15]) + (h_RT[7] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] -= q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    wdot[14] += q; /* HCO */
    wdot[15] -= q; /* CH2O */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[15];
    J[156] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[157] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[164] += dqdci;              /* dwdot[HCO]/d[HO2] */
    J[165] -= dqdci;              /* dwdot[CH2O]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[14];
    J[181] -= dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[182] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[189] += dqdci;              /* dwdot[HCO]/d[H2O2] */
    J[190] -= dqdci;              /* dwdot[CH2O]/d[H2O2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[7];
    J[356] -= dqdci;              /* dwdot[HO2]/d[HCO] */
    J[357] += dqdci;              /* dwdot[H2O2]/d[HCO] */
    J[364] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[365] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[6];
    J[381] -= dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[382] += dqdci;              /* dwdot[H2O2]/d[CH2O] */
    J[389] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[390] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[606] -= dqdT;               /* dwdot[HO2]/dT */
    J[607] += dqdT;               /* dwdot[H2O2]/dT */
    J[614] += dqdT;               /* dwdot[HCO]/dT */
    J[615] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 76: CH2 + O2 <=> OH + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[8];
    k_f = prefactor_units[75] * fwd_A[75]
                * exp(fwd_beta[75] * tc[0] - activation_units[75] * fwd_Ea[75] * invT);
    dlnkfdT = fwd_beta[75] * invT + activation_units[75] * fwd_Ea[75] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[14];
    Kc = exp(g_RT[3] - g_RT[4] + g_RT[8] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[4] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[8] -= q; /* CH2 */
    wdot[14] += q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[8];
    J[78] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[79] += dqdci;               /* dwdot[OH]/d[O2] */
    J[83] -= dqdci;               /* dwdot[CH2]/d[O2] */
    J[89] += dqdci;               /* dwdot[HCO]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[14];
    J[103] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[104] += dqdci;              /* dwdot[OH]/d[OH] */
    J[108] -= dqdci;              /* dwdot[CH2]/d[OH] */
    J[114] += dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[3];
    J[203] -= dqdci;              /* dwdot[O2]/d[CH2] */
    J[204] += dqdci;              /* dwdot[OH]/d[CH2] */
    J[208] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[214] += dqdci;              /* dwdot[HCO]/d[CH2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[4];
    J[353] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[354] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[358] -= dqdci;              /* dwdot[CH2]/d[HCO] */
    J[364] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[603] -= dqdT;               /* dwdot[O2]/dT */
    J[604] += dqdT;               /* dwdot[OH]/dT */
    J[608] -= dqdT;               /* dwdot[CH2]/dT */
    J[614] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 77: CH2 + H2 <=> H + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[8];
    k_f = prefactor_units[76] * fwd_A[76]
                * exp(fwd_beta[76] * tc[0] - activation_units[76] * fwd_Ea[76] * invT);
    dlnkfdT = fwd_beta[76] * invT + activation_units[76] * fwd_Ea[76] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[10];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[8] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[8]) + (h_RT[1] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[8] -= q; /* CH2 */
    wdot[10] += q; /* CH3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[8];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[8] -= dqdci;                /* dwdot[CH2]/d[H2] */
    J[10] += dqdci;               /* dwdot[CH3]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[10];
    J[25] -= dqdci;               /* dwdot[H2]/d[H] */
    J[26] += dqdci;               /* dwdot[H]/d[H] */
    J[33] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[35] += dqdci;               /* dwdot[CH3]/d[H] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[0];
    J[200] -= dqdci;              /* dwdot[H2]/d[CH2] */
    J[201] += dqdci;              /* dwdot[H]/d[CH2] */
    J[208] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[210] += dqdci;              /* dwdot[CH3]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[1];
    J[250] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[251] += dqdci;              /* dwdot[H]/d[CH3] */
    J[258] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[260] += dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[600] -= dqdT;               /* dwdot[H2]/dT */
    J[601] += dqdT;               /* dwdot[H]/dT */
    J[608] -= dqdT;               /* dwdot[CH2]/dT */
    J[610] += dqdT;               /* dwdot[CH3]/dT */

    /*reaction 78: 2 CH2 <=> H2 + C2H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[8];
    k_f = prefactor_units[77] * fwd_A[77]
                * exp(fwd_beta[77] * tc[0] - activation_units[77] * fwd_Ea[77] * invT);
    dlnkfdT = fwd_beta[77] * invT + activation_units[77] * fwd_Ea[77] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[17];
    Kc = exp(-g_RT[0] + 2*g_RT[8] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[8]) + (h_RT[0] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[8] -= 2 * q; /* CH2 */
    wdot[17] += q; /* C2H2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[17];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[8] += -2 * dqdci;           /* dwdot[CH2]/d[H2] */
    J[17] += dqdci;               /* dwdot[C2H2]/d[H2] */
    /* d()/d[CH2] */
    dqdci =  + k_f*2*sc[8];
    J[200] += dqdci;              /* dwdot[H2]/d[CH2] */
    J[208] += -2 * dqdci;         /* dwdot[CH2]/d[CH2] */
    J[217] += dqdci;              /* dwdot[C2H2]/d[CH2] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[0];
    J[425] += dqdci;              /* dwdot[H2]/d[C2H2] */
    J[433] += -2 * dqdci;         /* dwdot[CH2]/d[C2H2] */
    J[442] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    /* d()/dT */
    J[600] += dqdT;               /* dwdot[H2]/dT */
    J[608] += -2 * dqdT;          /* dwdot[CH2]/dT */
    J[617] += dqdT;               /* dwdot[C2H2]/dT */

    /*reaction 79: CH2 + CH3 <=> H + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[10];
    k_f = prefactor_units[78] * fwd_A[78]
                * exp(fwd_beta[78] * tc[0] - activation_units[78] * fwd_Ea[78] * invT);
    dlnkfdT = fwd_beta[78] * invT + activation_units[78] * fwd_Ea[78] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[19];
    Kc = exp(-g_RT[1] + g_RT[8] + g_RT[10] - g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[10]) + (h_RT[1] + h_RT[19]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[8] -= q; /* CH2 */
    wdot[10] -= q; /* CH3 */
    wdot[19] += q; /* C2H4 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[19];
    J[26] += dqdci;               /* dwdot[H]/d[H] */
    J[33] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[35] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[44] += dqdci;               /* dwdot[C2H4]/d[H] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[10];
    J[201] += dqdci;              /* dwdot[H]/d[CH2] */
    J[208] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[210] -= dqdci;              /* dwdot[CH3]/d[CH2] */
    J[219] += dqdci;              /* dwdot[C2H4]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[8];
    J[251] += dqdci;              /* dwdot[H]/d[CH3] */
    J[258] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[260] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[269] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[1];
    J[476] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[483] -= dqdci;              /* dwdot[CH2]/d[C2H4] */
    J[485] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[494] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[601] += dqdT;               /* dwdot[H]/dT */
    J[608] -= dqdT;               /* dwdot[CH2]/dT */
    J[610] -= dqdT;               /* dwdot[CH3]/dT */
    J[619] += dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 80: CH2 + CH4 <=> 2 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[11];
    k_f = prefactor_units[79] * fwd_A[79]
                * exp(fwd_beta[79] * tc[0] - activation_units[79] * fwd_Ea[79] * invT);
    dlnkfdT = fwd_beta[79] * invT + activation_units[79] * fwd_Ea[79] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[10];
    Kc = exp(g_RT[8] - 2*g_RT[10] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[11]) + (2*h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] -= q; /* CH2 */
    wdot[10] += 2 * q; /* CH3 */
    wdot[11] -= q; /* CH4 */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[11];
    J[208] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[210] += 2 * dqdci;          /* dwdot[CH3]/d[CH2] */
    J[211] -= dqdci;              /* dwdot[CH4]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*2*sc[10];
    J[258] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[260] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
    J[261] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[8];
    J[283] -= dqdci;              /* dwdot[CH2]/d[CH4] */
    J[285] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[286] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[608] -= dqdT;               /* dwdot[CH2]/dT */
    J[610] += 2 * dqdT;           /* dwdot[CH3]/dT */
    J[611] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 81: CH2(S) + N2 <=> CH2 + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[22];
    k_f = prefactor_units[80] * fwd_A[80]
                * exp(fwd_beta[80] * tc[0] - activation_units[80] * fwd_Ea[80] * invT);
    dlnkfdT = fwd_beta[80] * invT + activation_units[80] * fwd_Ea[80] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[22];
    Kc = exp(-g_RT[8] + g_RT[9] + g_RT[22] - g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[22]) + (h_RT[8] + h_RT[22]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] += q; /* CH2 */
    wdot[9] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[22];
    J[208] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[209] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[22];
    J[233] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[234] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[9] - k_r*sc[8];
    J[558] += dqdci;              /* dwdot[CH2]/d[N2] */
    J[559] -= dqdci;              /* dwdot[CH2(S)]/d[N2] */
    /* d()/dT */
    J[608] += dqdT;               /* dwdot[CH2]/dT */
    J[609] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 82: CH2(S) + AR <=> CH2 + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[23];
    k_f = prefactor_units[81] * fwd_A[81]
                * exp(fwd_beta[81] * tc[0] - activation_units[81] * fwd_Ea[81] * invT);
    dlnkfdT = fwd_beta[81] * invT + activation_units[81] * fwd_Ea[81] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[23];
    Kc = exp(-g_RT[8] + g_RT[9] + g_RT[23] - g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[23]) + (h_RT[8] + h_RT[23]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] += q; /* CH2 */
    wdot[9] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[23];
    J[208] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[209] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[23];
    J[233] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[234] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[AR] */
    dqdci =  + k_f*sc[9] - k_r*sc[8];
    J[583] += dqdci;              /* dwdot[CH2]/d[AR] */
    J[584] -= dqdci;              /* dwdot[CH2(S)]/d[AR] */
    /* d()/dT */
    J[608] += dqdT;               /* dwdot[CH2]/dT */
    J[609] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 83: CH2(S) + O2 <=> H + OH + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[9];
    k_f = prefactor_units[82] * fwd_A[82]
                * exp(fwd_beta[82] * tc[0] - activation_units[82] * fwd_Ea[82] * invT);
    dlnkfdT = fwd_beta[82] * invT + activation_units[82] * fwd_Ea[82] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[4]*sc[12];
    Kc = refC * exp(-g_RT[1] + g_RT[3] - g_RT[4] + g_RT[9] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[9]) + (h_RT[1] + h_RT[4] + h_RT[12]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[9] -= q; /* CH2(S) */
    wdot[12] += q; /* CO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4]*sc[12];
    J[26] += dqdci;               /* dwdot[H]/d[H] */
    J[28] -= dqdci;               /* dwdot[O2]/d[H] */
    J[29] += dqdci;               /* dwdot[OH]/d[H] */
    J[34] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[37] += dqdci;               /* dwdot[CO]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[9];
    J[76] += dqdci;               /* dwdot[H]/d[O2] */
    J[78] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[79] += dqdci;               /* dwdot[OH]/d[O2] */
    J[84] -= dqdci;               /* dwdot[CH2(S)]/d[O2] */
    J[87] += dqdci;               /* dwdot[CO]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1]*sc[12];
    J[101] += dqdci;              /* dwdot[H]/d[OH] */
    J[103] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[104] += dqdci;              /* dwdot[OH]/d[OH] */
    J[109] -= dqdci;              /* dwdot[CH2(S)]/d[OH] */
    J[112] += dqdci;              /* dwdot[CO]/d[OH] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[3];
    J[226] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[228] -= dqdci;              /* dwdot[O2]/d[CH2(S)] */
    J[229] += dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[234] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[237] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*sc[4];
    J[301] += dqdci;              /* dwdot[H]/d[CO] */
    J[303] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[304] += dqdci;              /* dwdot[OH]/d[CO] */
    J[309] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[312] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[601] += dqdT;               /* dwdot[H]/dT */
    J[603] -= dqdT;               /* dwdot[O2]/dT */
    J[604] += dqdT;               /* dwdot[OH]/dT */
    J[609] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[612] += dqdT;               /* dwdot[CO]/dT */

    /*reaction 84: CH2(S) + O2 <=> CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[9];
    k_f = prefactor_units[83] * fwd_A[83]
                * exp(fwd_beta[83] * tc[0] - activation_units[83] * fwd_Ea[83] * invT);
    dlnkfdT = fwd_beta[83] * invT + activation_units[83] * fwd_Ea[83] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[12];
    Kc = exp(g_RT[3] - g_RT[5] + g_RT[9] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[9]) + (h_RT[5] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[5] += q; /* H2O */
    wdot[9] -= q; /* CH2(S) */
    wdot[12] += q; /* CO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[9];
    J[78] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[80] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[84] -= dqdci;               /* dwdot[CH2(S)]/d[O2] */
    J[87] += dqdci;               /* dwdot[CO]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[12];
    J[128] -= dqdci;              /* dwdot[O2]/d[H2O] */
    J[130] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[134] -= dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[137] += dqdci;              /* dwdot[CO]/d[H2O] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[3];
    J[228] -= dqdci;              /* dwdot[O2]/d[CH2(S)] */
    J[230] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[234] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[237] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[5];
    J[303] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[305] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[309] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[312] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[603] -= dqdT;               /* dwdot[O2]/dT */
    J[605] += dqdT;               /* dwdot[H2O]/dT */
    J[609] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[612] += dqdT;               /* dwdot[CO]/dT */

    /*reaction 85: CH2(S) + H2 <=> CH3 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[9];
    k_f = prefactor_units[84] * fwd_A[84]
                * exp(fwd_beta[84] * tc[0] - activation_units[84] * fwd_Ea[84] * invT);
    dlnkfdT = fwd_beta[84] * invT + activation_units[84] * fwd_Ea[84] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[10];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[9] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[9]) + (h_RT[1] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[9] -= q; /* CH2(S) */
    wdot[10] += q; /* CH3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[9];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[9] -= dqdci;                /* dwdot[CH2(S)]/d[H2] */
    J[10] += dqdci;               /* dwdot[CH3]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[10];
    J[25] -= dqdci;               /* dwdot[H2]/d[H] */
    J[26] += dqdci;               /* dwdot[H]/d[H] */
    J[34] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[35] += dqdci;               /* dwdot[CH3]/d[H] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[0];
    J[225] -= dqdci;              /* dwdot[H2]/d[CH2(S)] */
    J[226] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[234] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[235] += dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[1];
    J[250] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[251] += dqdci;              /* dwdot[H]/d[CH3] */
    J[259] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[260] += dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[600] -= dqdT;               /* dwdot[H2]/dT */
    J[601] += dqdT;               /* dwdot[H]/dT */
    J[609] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[610] += dqdT;               /* dwdot[CH3]/dT */

    /*reaction 86: CH2(S) + H2O <=> CH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[9];
    k_f = prefactor_units[85] * fwd_A[85]
                * exp(fwd_beta[85] * tc[0] - activation_units[85] * fwd_Ea[85] * invT);
    dlnkfdT = fwd_beta[85] * invT + activation_units[85] * fwd_Ea[85] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[8];
    Kc = exp(g_RT[5] - g_RT[5] - g_RT[8] + g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[9]) + (h_RT[5] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] += q; /* CH2 */
    wdot[9] -= q; /* CH2(S) */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[9] - k_r*sc[8];
    J[133] += dqdci;              /* dwdot[CH2]/d[H2O] */
    J[134] -= dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[5];
    J[208] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[209] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[5];
    J[233] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[234] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/dT */
    J[608] += dqdT;               /* dwdot[CH2]/dT */
    J[609] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 87: CH2(S) + CH3 <=> H + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[10];
    k_f = prefactor_units[86] * fwd_A[86]
                * exp(fwd_beta[86] * tc[0] - activation_units[86] * fwd_Ea[86] * invT);
    dlnkfdT = fwd_beta[86] * invT + activation_units[86] * fwd_Ea[86] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[19];
    Kc = exp(-g_RT[1] + g_RT[9] + g_RT[10] - g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[10]) + (h_RT[1] + h_RT[19]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[9] -= q; /* CH2(S) */
    wdot[10] -= q; /* CH3 */
    wdot[19] += q; /* C2H4 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[19];
    J[26] += dqdci;               /* dwdot[H]/d[H] */
    J[34] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[35] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[44] += dqdci;               /* dwdot[C2H4]/d[H] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[10];
    J[226] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[234] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[235] -= dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    J[244] += dqdci;              /* dwdot[C2H4]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[9];
    J[251] += dqdci;              /* dwdot[H]/d[CH3] */
    J[259] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[260] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[269] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[1];
    J[476] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[484] -= dqdci;              /* dwdot[CH2(S)]/d[C2H4] */
    J[485] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[494] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[601] += dqdT;               /* dwdot[H]/dT */
    J[609] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[610] -= dqdT;               /* dwdot[CH3]/dT */
    J[619] += dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 88: CH2(S) + CH4 <=> 2 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[11];
    k_f = prefactor_units[87] * fwd_A[87]
                * exp(fwd_beta[87] * tc[0] - activation_units[87] * fwd_Ea[87] * invT);
    dlnkfdT = fwd_beta[87] * invT + activation_units[87] * fwd_Ea[87] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[10];
    Kc = exp(g_RT[9] - 2*g_RT[10] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[11]) + (2*h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] -= q; /* CH2(S) */
    wdot[10] += 2 * q; /* CH3 */
    wdot[11] -= q; /* CH4 */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[11];
    J[234] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[235] += 2 * dqdci;          /* dwdot[CH3]/d[CH2(S)] */
    J[236] -= dqdci;              /* dwdot[CH4]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*2*sc[10];
    J[259] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[260] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
    J[261] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[9];
    J[284] -= dqdci;              /* dwdot[CH2(S)]/d[CH4] */
    J[285] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[286] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[609] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[610] += 2 * dqdT;           /* dwdot[CH3]/dT */
    J[611] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 89: CH2(S) + CO <=> CH2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[12];
    k_f = prefactor_units[88] * fwd_A[88]
                * exp(fwd_beta[88] * tc[0] - activation_units[88] * fwd_Ea[88] * invT);
    dlnkfdT = fwd_beta[88] * invT + activation_units[88] * fwd_Ea[88] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[12];
    Kc = exp(-g_RT[8] + g_RT[9] + g_RT[12] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[12]) + (h_RT[8] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] += q; /* CH2 */
    wdot[9] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[12];
    J[208] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[209] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[12];
    J[233] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[234] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[9] - k_r*sc[8];
    J[308] += dqdci;              /* dwdot[CH2]/d[CO] */
    J[309] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    /* d()/dT */
    J[608] += dqdT;               /* dwdot[CH2]/dT */
    J[609] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 90: CH2(S) + CO2 <=> CH2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[13];
    k_f = prefactor_units[89] * fwd_A[89]
                * exp(fwd_beta[89] * tc[0] - activation_units[89] * fwd_Ea[89] * invT);
    dlnkfdT = fwd_beta[89] * invT + activation_units[89] * fwd_Ea[89] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[13];
    Kc = exp(-g_RT[8] + g_RT[9] + g_RT[13] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[13]) + (h_RT[8] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] += q; /* CH2 */
    wdot[9] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[13];
    J[208] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[209] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[13];
    J[233] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[234] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[9] - k_r*sc[8];
    J[333] += dqdci;              /* dwdot[CH2]/d[CO2] */
    J[334] -= dqdci;              /* dwdot[CH2(S)]/d[CO2] */
    /* d()/dT */
    J[608] += dqdT;               /* dwdot[CH2]/dT */
    J[609] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 91: CH2(S) + CO2 <=> CO + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[13];
    k_f = prefactor_units[90] * fwd_A[90]
                * exp(fwd_beta[90] * tc[0] - activation_units[90] * fwd_Ea[90] * invT);
    dlnkfdT = fwd_beta[90] * invT + activation_units[90] * fwd_Ea[90] * invT2;
    /* reverse */
    phi_r = sc[12]*sc[15];
    Kc = exp(g_RT[9] - g_RT[12] + g_RT[13] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[13]) + (h_RT[12] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] -= q; /* CH2(S) */
    wdot[12] += q; /* CO */
    wdot[13] -= q; /* CO2 */
    wdot[15] += q; /* CH2O */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[13];
    J[234] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[237] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    J[238] -= dqdci;              /* dwdot[CO2]/d[CH2(S)] */
    J[240] += dqdci;              /* dwdot[CH2O]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[15];
    J[309] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[312] += dqdci;              /* dwdot[CO]/d[CO] */
    J[313] -= dqdci;              /* dwdot[CO2]/d[CO] */
    J[315] += dqdci;              /* dwdot[CH2O]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[9];
    J[334] -= dqdci;              /* dwdot[CH2(S)]/d[CO2] */
    J[337] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[338] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    J[340] += dqdci;              /* dwdot[CH2O]/d[CO2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[12];
    J[384] -= dqdci;              /* dwdot[CH2(S)]/d[CH2O] */
    J[387] += dqdci;              /* dwdot[CO]/d[CH2O] */
    J[388] -= dqdci;              /* dwdot[CO2]/d[CH2O] */
    J[390] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[609] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[612] += dqdT;               /* dwdot[CO]/dT */
    J[613] -= dqdT;               /* dwdot[CO2]/dT */
    J[615] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 92: CH3 + O2 <=> O + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[10];
    k_f = prefactor_units[91] * fwd_A[91]
                * exp(fwd_beta[91] * tc[0] - activation_units[91] * fwd_Ea[91] * invT);
    dlnkfdT = fwd_beta[91] * invT + activation_units[91] * fwd_Ea[91] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[16];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[10] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[10]) + (h_RT[2] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[10] -= q; /* CH3 */
    wdot[16] += q; /* CH3O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[16];
    J[52] += dqdci;               /* dwdot[O]/d[O] */
    J[53] -= dqdci;               /* dwdot[O2]/d[O] */
    J[60] -= dqdci;               /* dwdot[CH3]/d[O] */
    J[66] += dqdci;               /* dwdot[CH3O]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[10];
    J[77] += dqdci;               /* dwdot[O]/d[O2] */
    J[78] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[85] -= dqdci;               /* dwdot[CH3]/d[O2] */
    J[91] += dqdci;               /* dwdot[CH3O]/d[O2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[252] += dqdci;              /* dwdot[O]/d[CH3] */
    J[253] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[260] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[266] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[2];
    J[402] += dqdci;              /* dwdot[O]/d[CH3O] */
    J[403] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    J[410] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[416] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[602] += dqdT;               /* dwdot[O]/dT */
    J[603] -= dqdT;               /* dwdot[O2]/dT */
    J[610] -= dqdT;               /* dwdot[CH3]/dT */
    J[616] += dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 93: CH3 + O2 <=> OH + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[10];
    k_f = prefactor_units[92] * fwd_A[92]
                * exp(fwd_beta[92] * tc[0] - activation_units[92] * fwd_Ea[92] * invT);
    dlnkfdT = fwd_beta[92] * invT + activation_units[92] * fwd_Ea[92] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[15];
    Kc = exp(g_RT[3] - g_RT[4] + g_RT[10] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[10]) + (h_RT[4] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[10] -= q; /* CH3 */
    wdot[15] += q; /* CH2O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[10];
    J[78] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[79] += dqdci;               /* dwdot[OH]/d[O2] */
    J[85] -= dqdci;               /* dwdot[CH3]/d[O2] */
    J[90] += dqdci;               /* dwdot[CH2O]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[15];
    J[103] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[104] += dqdci;              /* dwdot[OH]/d[OH] */
    J[110] -= dqdci;              /* dwdot[CH3]/d[OH] */
    J[115] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[253] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[254] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[260] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[265] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[4];
    J[378] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[379] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[385] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[390] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[603] -= dqdT;               /* dwdot[O2]/dT */
    J[604] += dqdT;               /* dwdot[OH]/dT */
    J[610] -= dqdT;               /* dwdot[CH3]/dT */
    J[615] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 94: CH3 + H2O2 <=> HO2 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[10];
    k_f = prefactor_units[93] * fwd_A[93]
                * exp(fwd_beta[93] * tc[0] - activation_units[93] * fwd_Ea[93] * invT);
    dlnkfdT = fwd_beta[93] * invT + activation_units[93] * fwd_Ea[93] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[11];
    Kc = exp(-g_RT[6] + g_RT[7] + g_RT[10] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[10]) + (h_RT[6] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    wdot[10] -= q; /* CH3 */
    wdot[11] += q; /* CH4 */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[11];
    J[156] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[157] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[160] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[161] += dqdci;              /* dwdot[CH4]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[10];
    J[181] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[182] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[185] -= dqdci;              /* dwdot[CH3]/d[H2O2] */
    J[186] += dqdci;              /* dwdot[CH4]/d[H2O2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[7];
    J[256] += dqdci;              /* dwdot[HO2]/d[CH3] */
    J[257] -= dqdci;              /* dwdot[H2O2]/d[CH3] */
    J[260] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[261] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[6];
    J[281] += dqdci;              /* dwdot[HO2]/d[CH4] */
    J[282] -= dqdci;              /* dwdot[H2O2]/d[CH4] */
    J[285] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[286] += dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[606] += dqdT;               /* dwdot[HO2]/dT */
    J[607] -= dqdT;               /* dwdot[H2O2]/dT */
    J[610] -= dqdT;               /* dwdot[CH3]/dT */
    J[611] += dqdT;               /* dwdot[CH4]/dT */

    /*reaction 95: 2 CH3 <=> H + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[10];
    k_f = prefactor_units[94] * fwd_A[94]
                * exp(fwd_beta[94] * tc[0] - activation_units[94] * fwd_Ea[94] * invT);
    dlnkfdT = fwd_beta[94] * invT + activation_units[94] * fwd_Ea[94] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[20];
    Kc = exp(-g_RT[1] + 2*g_RT[10] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[10]) + (h_RT[1] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[10] -= 2 * q; /* CH3 */
    wdot[20] += q; /* C2H5 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[20];
    J[26] += dqdci;               /* dwdot[H]/d[H] */
    J[35] += -2 * dqdci;          /* dwdot[CH3]/d[H] */
    J[45] += dqdci;               /* dwdot[C2H5]/d[H] */
    /* d()/d[CH3] */
    dqdci =  + k_f*2*sc[10];
    J[251] += dqdci;              /* dwdot[H]/d[CH3] */
    J[260] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
    J[270] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[1];
    J[501] += dqdci;              /* dwdot[H]/d[C2H5] */
    J[510] += -2 * dqdci;         /* dwdot[CH3]/d[C2H5] */
    J[520] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[601] += dqdT;               /* dwdot[H]/dT */
    J[610] += -2 * dqdT;          /* dwdot[CH3]/dT */
    J[620] += dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 96: CH3 + HCO <=> CH4 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[14];
    k_f = prefactor_units[95] * fwd_A[95]
                * exp(fwd_beta[95] * tc[0] - activation_units[95] * fwd_Ea[95] * invT);
    dlnkfdT = fwd_beta[95] * invT + activation_units[95] * fwd_Ea[95] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[12];
    Kc = exp(g_RT[10] - g_RT[11] - g_RT[12] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[14]) + (h_RT[11] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] -= q; /* CH3 */
    wdot[11] += q; /* CH4 */
    wdot[12] += q; /* CO */
    wdot[14] -= q; /* HCO */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[14];
    J[260] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[261] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[262] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[264] -= dqdci;              /* dwdot[HCO]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[12];
    J[285] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[286] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[287] += dqdci;              /* dwdot[CO]/d[CH4] */
    J[289] -= dqdci;              /* dwdot[HCO]/d[CH4] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[11];
    J[310] -= dqdci;              /* dwdot[CH3]/d[CO] */
    J[311] += dqdci;              /* dwdot[CH4]/d[CO] */
    J[312] += dqdci;              /* dwdot[CO]/d[CO] */
    J[314] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[10];
    J[360] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[361] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[362] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[364] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[610] -= dqdT;               /* dwdot[CH3]/dT */
    J[611] += dqdT;               /* dwdot[CH4]/dT */
    J[612] += dqdT;               /* dwdot[CO]/dT */
    J[614] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 97: CH3 + CH2O <=> HCO + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[15];
    k_f = prefactor_units[96] * fwd_A[96]
                * exp(fwd_beta[96] * tc[0] - activation_units[96] * fwd_Ea[96] * invT);
    dlnkfdT = fwd_beta[96] * invT + activation_units[96] * fwd_Ea[96] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[14];
    Kc = exp(g_RT[10] - g_RT[11] - g_RT[14] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[15]) + (h_RT[11] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] -= q; /* CH3 */
    wdot[11] += q; /* CH4 */
    wdot[14] += q; /* HCO */
    wdot[15] -= q; /* CH2O */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[15];
    J[260] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[261] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[264] += dqdci;              /* dwdot[HCO]/d[CH3] */
    J[265] -= dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[14];
    J[285] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[286] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[289] += dqdci;              /* dwdot[HCO]/d[CH4] */
    J[290] -= dqdci;              /* dwdot[CH2O]/d[CH4] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[11];
    J[360] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[361] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[364] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[365] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[10];
    J[385] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[386] += dqdci;              /* dwdot[CH4]/d[CH2O] */
    J[389] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[390] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[610] -= dqdT;               /* dwdot[CH3]/dT */
    J[611] += dqdT;               /* dwdot[CH4]/dT */
    J[614] += dqdT;               /* dwdot[HCO]/dT */
    J[615] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 98: CH3 + C2H4 <=> C2H3 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[19];
    k_f = prefactor_units[97] * fwd_A[97]
                * exp(fwd_beta[97] * tc[0] - activation_units[97] * fwd_Ea[97] * invT);
    dlnkfdT = fwd_beta[97] * invT + activation_units[97] * fwd_Ea[97] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[18];
    Kc = exp(g_RT[10] - g_RT[11] - g_RT[18] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[19]) + (h_RT[11] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] -= q; /* CH3 */
    wdot[11] += q; /* CH4 */
    wdot[18] += q; /* C2H3 */
    wdot[19] -= q; /* C2H4 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[19];
    J[260] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[261] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[268] += dqdci;              /* dwdot[C2H3]/d[CH3] */
    J[269] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[18];
    J[285] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[286] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[293] += dqdci;              /* dwdot[C2H3]/d[CH4] */
    J[294] -= dqdci;              /* dwdot[C2H4]/d[CH4] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[11];
    J[460] -= dqdci;              /* dwdot[CH3]/d[C2H3] */
    J[461] += dqdci;              /* dwdot[CH4]/d[C2H3] */
    J[468] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[469] -= dqdci;              /* dwdot[C2H4]/d[C2H3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[10];
    J[485] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[486] += dqdci;              /* dwdot[CH4]/d[C2H4] */
    J[493] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    J[494] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[610] -= dqdT;               /* dwdot[CH3]/dT */
    J[611] += dqdT;               /* dwdot[CH4]/dT */
    J[618] += dqdT;               /* dwdot[C2H3]/dT */
    J[619] -= dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 99: CH3 + C2H6 <=> C2H5 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[21];
    k_f = prefactor_units[98] * fwd_A[98]
                * exp(fwd_beta[98] * tc[0] - activation_units[98] * fwd_Ea[98] * invT);
    dlnkfdT = fwd_beta[98] * invT + activation_units[98] * fwd_Ea[98] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[20];
    Kc = exp(g_RT[10] - g_RT[11] - g_RT[20] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[21]) + (h_RT[11] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] -= q; /* CH3 */
    wdot[11] += q; /* CH4 */
    wdot[20] += q; /* C2H5 */
    wdot[21] -= q; /* C2H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[21];
    J[260] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[261] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[270] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    J[271] -= dqdci;              /* dwdot[C2H6]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[20];
    J[285] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[286] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[295] += dqdci;              /* dwdot[C2H5]/d[CH4] */
    J[296] -= dqdci;              /* dwdot[C2H6]/d[CH4] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[11];
    J[510] -= dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[511] += dqdci;              /* dwdot[CH4]/d[C2H5] */
    J[520] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[521] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[10];
    J[535] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
    J[536] += dqdci;              /* dwdot[CH4]/d[C2H6] */
    J[545] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[546] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[610] -= dqdT;               /* dwdot[CH3]/dT */
    J[611] += dqdT;               /* dwdot[CH4]/dT */
    J[620] += dqdT;               /* dwdot[C2H5]/dT */
    J[621] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 100: HCO + H2O <=> H + CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[14];
    k_f = prefactor_units[99] * fwd_A[99]
                * exp(fwd_beta[99] * tc[0] - activation_units[99] * fwd_Ea[99] * invT);
    dlnkfdT = fwd_beta[99] * invT + activation_units[99] * fwd_Ea[99] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5]*sc[12];
    Kc = refC * exp(-g_RT[1] + g_RT[5] - g_RT[5] - g_RT[12] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[14]) + (h_RT[1] + h_RT[5] + h_RT[12]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[12] += q; /* CO */
    wdot[14] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5]*sc[12];
    J[26] += dqdci;               /* dwdot[H]/d[H] */
    J[37] += dqdci;               /* dwdot[CO]/d[H] */
    J[39] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[14] - k_r*sc[1]*sc[12];
    J[126] += dqdci;              /* dwdot[H]/d[H2O] */
    J[137] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[139] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*sc[5];
    J[301] += dqdci;              /* dwdot[H]/d[CO] */
    J[312] += dqdci;              /* dwdot[CO]/d[CO] */
    J[314] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[5];
    J[351] += dqdci;              /* dwdot[H]/d[HCO] */
    J[362] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[364] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[601] += dqdT;               /* dwdot[H]/dT */
    J[612] += dqdT;               /* dwdot[CO]/dT */
    J[614] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 101: HCO + O2 <=> HO2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[14];
    k_f = prefactor_units[100] * fwd_A[100]
                * exp(fwd_beta[100] * tc[0] - activation_units[100] * fwd_Ea[100] * invT);
    dlnkfdT = fwd_beta[100] * invT + activation_units[100] * fwd_Ea[100] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[12];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[12] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[14]) + (h_RT[6] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[12] += q; /* CO */
    wdot[14] -= q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[14];
    J[78] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[81] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[87] += dqdci;               /* dwdot[CO]/d[O2] */
    J[89] -= dqdci;               /* dwdot[HCO]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[12];
    J[153] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[156] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[162] += dqdci;              /* dwdot[CO]/d[HO2] */
    J[164] -= dqdci;              /* dwdot[HCO]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[6];
    J[303] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[306] += dqdci;              /* dwdot[HO2]/d[CO] */
    J[312] += dqdci;              /* dwdot[CO]/d[CO] */
    J[314] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[3];
    J[353] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[356] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[362] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[364] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[603] -= dqdT;               /* dwdot[O2]/dT */
    J[606] += dqdT;               /* dwdot[HO2]/dT */
    J[612] += dqdT;               /* dwdot[CO]/dT */
    J[614] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 102: CH3O + O2 <=> HO2 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[16];
    k_f = prefactor_units[101] * fwd_A[101]
                * exp(fwd_beta[101] * tc[0] - activation_units[101] * fwd_Ea[101] * invT);
    dlnkfdT = fwd_beta[101] * invT + activation_units[101] * fwd_Ea[101] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[15];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[15] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[16]) + (h_RT[6] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[15] += q; /* CH2O */
    wdot[16] -= q; /* CH3O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[16];
    J[78] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[81] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[90] += dqdci;               /* dwdot[CH2O]/d[O2] */
    J[91] -= dqdci;               /* dwdot[CH3O]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[15];
    J[153] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[156] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[165] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[166] -= dqdci;              /* dwdot[CH3O]/d[HO2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[6];
    J[378] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[381] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[390] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[391] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[3];
    J[403] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    J[406] += dqdci;              /* dwdot[HO2]/d[CH3O] */
    J[415] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[416] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[603] -= dqdT;               /* dwdot[O2]/dT */
    J[606] += dqdT;               /* dwdot[HO2]/dT */
    J[615] += dqdT;               /* dwdot[CH2O]/dT */
    J[616] -= dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 103: C2H3 + O2 <=> HCO + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[18];
    k_f = prefactor_units[102] * fwd_A[102]
                * exp(fwd_beta[102] * tc[0] - activation_units[102] * fwd_Ea[102] * invT);
    dlnkfdT = fwd_beta[102] * invT + activation_units[102] * fwd_Ea[102] * invT2;
    /* reverse */
    phi_r = sc[14]*sc[15];
    Kc = exp(g_RT[3] - g_RT[14] - g_RT[15] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[18]) + (h_RT[14] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[14] += q; /* HCO */
    wdot[15] += q; /* CH2O */
    wdot[18] -= q; /* C2H3 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[18];
    J[78] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[89] += dqdci;               /* dwdot[HCO]/d[O2] */
    J[90] += dqdci;               /* dwdot[CH2O]/d[O2] */
    J[93] -= dqdci;               /* dwdot[C2H3]/d[O2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[15];
    J[353] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[364] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[365] += dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[368] -= dqdci;              /* dwdot[C2H3]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[14];
    J[378] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[389] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[390] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[393] -= dqdci;              /* dwdot[C2H3]/d[CH2O] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[3];
    J[453] -= dqdci;              /* dwdot[O2]/d[C2H3] */
    J[464] += dqdci;              /* dwdot[HCO]/d[C2H3] */
    J[465] += dqdci;              /* dwdot[CH2O]/d[C2H3] */
    J[468] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[603] -= dqdT;               /* dwdot[O2]/dT */
    J[614] += dqdT;               /* dwdot[HCO]/dT */
    J[615] += dqdT;               /* dwdot[CH2O]/dT */
    J[618] -= dqdT;               /* dwdot[C2H3]/dT */

    /*reaction 104: C2H5 + O2 <=> HO2 + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[20];
    k_f = prefactor_units[103] * fwd_A[103]
                * exp(fwd_beta[103] * tc[0] - activation_units[103] * fwd_Ea[103] * invT);
    dlnkfdT = fwd_beta[103] * invT + activation_units[103] * fwd_Ea[103] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[19];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[19] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[20]) + (h_RT[6] + h_RT[19]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[19] += q; /* C2H4 */
    wdot[20] -= q; /* C2H5 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[20];
    J[78] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[81] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[94] += dqdci;               /* dwdot[C2H4]/d[O2] */
    J[95] -= dqdci;               /* dwdot[C2H5]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[19];
    J[153] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[156] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[169] += dqdci;              /* dwdot[C2H4]/d[HO2] */
    J[170] -= dqdci;              /* dwdot[C2H5]/d[HO2] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[6];
    J[478] -= dqdci;              /* dwdot[O2]/d[C2H4] */
    J[481] += dqdci;              /* dwdot[HO2]/d[C2H4] */
    J[494] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[495] -= dqdci;              /* dwdot[C2H5]/d[C2H4] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[3];
    J[503] -= dqdci;              /* dwdot[O2]/d[C2H5] */
    J[506] += dqdci;              /* dwdot[HO2]/d[C2H5] */
    J[519] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    J[520] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[603] -= dqdT;               /* dwdot[O2]/dT */
    J[606] += dqdT;               /* dwdot[HO2]/dT */
    J[619] += dqdT;               /* dwdot[C2H4]/dT */
    J[620] -= dqdT;               /* dwdot[C2H5]/dT */

    double c_R[24], dcRdT[24], e_RT[24];
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
    for (int k = 0; k < 24; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[600+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 24; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 24; ++m) {
            dehmixdc += eh_RT[m]*J[k*25+m];
        }
        J[k*25+24] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[624] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
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
        /*species 8: CH2 */
        species[8] =
            +9.68872143e-04
            +5.58979682e-06 * tc[1]
            -1.15527346e-08 * tc[2]
            +6.74966876e-12 * tc[3];
        /*species 9: CH2(S) */
        species[9] =
            -2.36661419e-03
            +1.64659244e-05 * tc[1]
            -2.00644794e-08 * tc[2]
            +7.77258948e-12 * tc[3];
        /*species 10: CH3 */
        species[10] =
            +2.01095175e-03
            +1.14604371e-05 * tc[1]
            -2.06135228e-08 * tc[2]
            +1.01754294e-11 * tc[3];
        /*species 11: CH4 */
        species[11] =
            -1.36709788e-02
            +9.83601198e-05 * tc[1]
            -1.45422908e-07 * tc[2]
            +6.66775824e-11 * tc[3];
        /*species 12: CO */
        species[12] =
            -6.10353680e-04
            +2.03362866e-06 * tc[1]
            +2.72101765e-09 * tc[2]
            -3.61769800e-12 * tc[3];
        /*species 13: CO2 */
        species[13] =
            +8.98459677e-03
            -1.42471254e-05 * tc[1]
            +7.37757066e-09 * tc[2]
            -5.74798192e-13 * tc[3];
        /*species 14: HCO */
        species[14] =
            -3.24392532e-03
            +2.75598892e-05 * tc[1]
            -3.99432279e-08 * tc[2]
            +1.73507546e-11 * tc[3];
        /*species 15: CH2O */
        species[15] =
            -9.90833369e-03
            +7.46440016e-05 * tc[1]
            -1.13785578e-07 * tc[2]
            +5.27090608e-11 * tc[3];
        /*species 16: CH3O */
        species[16] =
            +7.21659500e-03
            +1.06769440e-05 * tc[1]
            -2.21329080e-08 * tc[2]
            +8.30244000e-12 * tc[3];
        /*species 17: C2H2 */
        species[17] =
            +2.33615629e-02
            -7.10343630e-05 * tc[1]
            +8.40457311e-08 * tc[2]
            -3.40029190e-11 * tc[3];
        /*species 18: C2H3 */
        species[18] =
            +1.51479162e-03
            +5.18418824e-05 * tc[1]
            -1.07297354e-07 * tc[2]
            +5.88603492e-11 * tc[3];
        /*species 19: C2H4 */
        species[19] =
            -7.57052247e-03
            +1.14198058e-04 * tc[1]
            -2.07476626e-07 * tc[2]
            +1.07953749e-10 * tc[3];
        /*species 20: C2H5 */
        species[20] =
            -4.18658892e-03
            +9.94285614e-05 * tc[1]
            -1.79737982e-07 * tc[2]
            +9.22036016e-11 * tc[3];
        /*species 21: C2H6 */
        species[21] =
            -5.50154270e-03
            +1.19887658e-04 * tc[1]
            -2.12539886e-07 * tc[2]
            +1.07474308e-10 * tc[3];
        /*species 22: N2 */
        species[22] =
            +1.40824040e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77941600e-12 * tc[3];
        /*species 23: AR */
        species[23] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
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
        /*species 8: CH2 */
        species[8] =
            +3.65639292e-03
            -2.81789194e-06 * tc[1]
            +7.80538647e-10 * tc[2]
            -7.50910268e-14 * tc[3];
        /*species 9: CH2(S) */
        species[9] =
            +4.65588637e-03
            -4.02383894e-06 * tc[1]
            +1.25371800e-09 * tc[2]
            -1.35886546e-13 * tc[3];
        /*species 10: CH3 */
        species[10] =
            +7.23990037e-03
            -5.97428696e-06 * tc[1]
            +1.78705393e-09 * tc[2]
            -1.86861758e-13 * tc[3];
        /*species 11: CH4 */
        species[11] =
            +1.33909467e-02
            -1.14657162e-05 * tc[1]
            +3.66877605e-09 * tc[2]
            -4.07260920e-13 * tc[3];
        /*species 12: CO */
        species[12] =
            +2.06252743e-03
            -1.99765154e-06 * tc[1]
            +6.90159024e-10 * tc[2]
            -8.14590864e-14 * tc[3];
        /*species 13: CO2 */
        species[13] =
            +4.41437026e-03
            -4.42962808e-06 * tc[1]
            +1.57047056e-09 * tc[2]
            -1.88833666e-13 * tc[3];
        /*species 14: HCO */
        species[14] =
            +4.95695526e-03
            -4.96891226e-06 * tc[1]
            +1.76748533e-09 * tc[2]
            -2.13403484e-13 * tc[3];
        /*species 15: CH2O */
        species[15] =
            +9.20000082e-03
            -8.84517626e-06 * tc[1]
            +3.01923636e-09 * tc[2]
            -3.53542256e-13 * tc[3];
        /*species 16: CH3O */
        species[16] =
            +7.87149700e-03
            -5.31276800e-06 * tc[1]
            +1.18332930e-09 * tc[2]
            -8.45046400e-14 * tc[3];
        /*species 17: C2H2 */
        species[17] =
            +5.96166664e-03
            -4.74589704e-06 * tc[1]
            +1.40223651e-09 * tc[2]
            -1.44494085e-13 * tc[3];
        /*species 18: C2H3 */
        species[18] =
            +1.03302292e-02
            -9.36164698e-06 * tc[1]
            +3.05289864e-09 * tc[2]
            -3.45042816e-13 * tc[3];
        /*species 19: C2H4 */
        species[19] =
            +1.46454151e-02
            -1.34215583e-05 * tc[1]
            +4.41668769e-09 * tc[2]
            -5.02824244e-13 * tc[3];
        /*species 20: C2H5 */
        species[20] =
            +1.73972722e-02
            -1.59641334e-05 * tc[1]
            +5.25653067e-09 * tc[2]
            -5.98566304e-13 * tc[3];
        /*species 21: C2H6 */
        species[21] =
            +2.16852677e-02
            -2.00512134e-05 * tc[1]
            +6.64236003e-09 * tc[2]
            -7.60011560e-13 * tc[3];
        /*species 22: N2 */
        species[22] =
            +1.48797680e-03
            -1.13695200e-06 * tc[1]
            +3.02911140e-10 * tc[2]
            -2.70134040e-14 * tc[3];
        /*species 23: AR */
        species[23] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
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

    double q_f[104], q_r[104];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 104; ++i) {
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

    /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
    kc[0] = 1.0 / (refC) * exp((g_RT[1] + g_RT[8]) - (g_RT[10]));

    /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
    kc[1] = 1.0 / (refC) * exp((g_RT[1] + g_RT[10]) - (g_RT[11]));

    /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
    kc[2] = 1.0 / (refC) * exp((g_RT[1] + g_RT[14]) - (g_RT[15]));

    /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
    kc[3] = 1.0 / (refC) * exp((g_RT[1] + g_RT[15]) - (g_RT[16]));

    /*reaction 5: H + C2H2 (+M) <=> C2H3 (+M) */
    kc[4] = 1.0 / (refC) * exp((g_RT[1] + g_RT[17]) - (g_RT[18]));

    /*reaction 6: H + C2H3 (+M) <=> C2H4 (+M) */
    kc[5] = 1.0 / (refC) * exp((g_RT[1] + g_RT[18]) - (g_RT[19]));

    /*reaction 7: H + C2H4 (+M) <=> C2H5 (+M) */
    kc[6] = 1.0 / (refC) * exp((g_RT[1] + g_RT[19]) - (g_RT[20]));

    /*reaction 8: H + C2H5 (+M) <=> C2H6 (+M) */
    kc[7] = 1.0 / (refC) * exp((g_RT[1] + g_RT[20]) - (g_RT[21]));

    /*reaction 9: H2 + CO (+M) <=> CH2O (+M) */
    kc[8] = 1.0 / (refC) * exp((g_RT[0] + g_RT[12]) - (g_RT[15]));

    /*reaction 10: 2 OH (+M) <=> H2O2 (+M) */
    kc[9] = 1.0 / (refC) * exp((2 * g_RT[4]) - (g_RT[7]));

    /*reaction 11: 2 CH3 (+M) <=> C2H6 (+M) */
    kc[10] = 1.0 / (refC) * exp((2 * g_RT[10]) - (g_RT[21]));

    /*reaction 12: C2H4 (+M) <=> H2 + C2H2 (+M) */
    kc[11] = refC * exp((g_RT[19]) - (g_RT[0] + g_RT[17]));

    /*reaction 13: O + H + M <=> OH + M */
    kc[12] = 1.0 / (refC) * exp((g_RT[2] + g_RT[1]) - (g_RT[4]));

    /*reaction 14: O + CO + M <=> CO2 + M */
    kc[13] = 1.0 / (refC) * exp((g_RT[2] + g_RT[12]) - (g_RT[13]));

    /*reaction 15: H + O2 + M <=> HO2 + M */
    kc[14] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));

    /*reaction 16: 2 H + M <=> H2 + M */
    kc[15] = 1.0 / (refC) * exp((2 * g_RT[1]) - (g_RT[0]));

    /*reaction 17: H + OH + M <=> H2O + M */
    kc[16] = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));

    /*reaction 18: HCO + M <=> H + CO + M */
    kc[17] = refC * exp((g_RT[14]) - (g_RT[1] + g_RT[12]));

    /*reaction 19: O + H2 <=> H + OH */
    kc[18] = exp((g_RT[2] + g_RT[0]) - (g_RT[1] + g_RT[4]));

    /*reaction 20: O + HO2 <=> OH + O2 */
    kc[19] = exp((g_RT[2] + g_RT[6]) - (g_RT[4] + g_RT[3]));

    /*reaction 21: O + CH2 <=> H + HCO */
    kc[20] = exp((g_RT[2] + g_RT[8]) - (g_RT[1] + g_RT[14]));

    /*reaction 22: O + CH2(S) <=> H + HCO */
    kc[21] = exp((g_RT[2] + g_RT[9]) - (g_RT[1] + g_RT[14]));

    /*reaction 23: O + CH3 <=> H + CH2O */
    kc[22] = exp((g_RT[2] + g_RT[10]) - (g_RT[1] + g_RT[15]));

    /*reaction 24: O + CH4 <=> OH + CH3 */
    kc[23] = exp((g_RT[2] + g_RT[11]) - (g_RT[4] + g_RT[10]));

    /*reaction 25: O + HCO <=> OH + CO */
    kc[24] = exp((g_RT[2] + g_RT[14]) - (g_RT[4] + g_RT[12]));

    /*reaction 26: O + HCO <=> H + CO2 */
    kc[25] = exp((g_RT[2] + g_RT[14]) - (g_RT[1] + g_RT[13]));

    /*reaction 27: O + CH2O <=> OH + HCO */
    kc[26] = exp((g_RT[2] + g_RT[15]) - (g_RT[4] + g_RT[14]));

    /*reaction 28: O + C2H2 <=> CH2(S) + CO */
    kc[27] = exp((g_RT[2] + g_RT[17]) - (g_RT[9] + g_RT[12]));

    /*reaction 29: O + C2H2 <=> CO + CH2 */
    kc[28] = exp((g_RT[2] + g_RT[17]) - (g_RT[12] + g_RT[8]));

    /*reaction 30: O + C2H4 <=> CH3 + HCO */
    kc[29] = exp((g_RT[2] + g_RT[19]) - (g_RT[10] + g_RT[14]));

    /*reaction 31: O + C2H5 <=> CH3 + CH2O */
    kc[30] = exp((g_RT[2] + g_RT[20]) - (g_RT[10] + g_RT[15]));

    /*reaction 32: O + C2H6 <=> OH + C2H5 */
    kc[31] = exp((g_RT[2] + g_RT[21]) - (g_RT[4] + g_RT[20]));

    /*reaction 33: O2 + CO <=> O + CO2 */
    kc[32] = exp((g_RT[3] + g_RT[12]) - (g_RT[2] + g_RT[13]));

    /*reaction 34: O2 + CH2O <=> HO2 + HCO */
    kc[33] = exp((g_RT[3] + g_RT[15]) - (g_RT[6] + g_RT[14]));

    /*reaction 35: H + 2 O2 <=> HO2 + O2 */
    kc[34] = 1.0 / (refC) * exp((g_RT[1] + 2 * g_RT[3]) - (g_RT[6] + g_RT[3]));

    /*reaction 36: H + O2 + H2O <=> HO2 + H2O */
    kc[35] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));

    /*reaction 37: H + O2 + N2 <=> HO2 + N2 */
    kc[36] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[22]) - (g_RT[6] + g_RT[22]));

    /*reaction 38: H + O2 + AR <=> HO2 + AR */
    kc[37] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[23]) - (g_RT[6] + g_RT[23]));

    /*reaction 39: H + O2 <=> O + OH */
    kc[38] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));

    /*reaction 40: 2 H + H2 <=> 2 H2 */
    kc[39] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[0]) - (2 * g_RT[0]));

    /*reaction 41: 2 H + H2O <=> H2 + H2O */
    kc[40] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));

    /*reaction 42: 2 H + CO2 <=> H2 + CO2 */
    kc[41] = 1.0 / (refC) * exp((2 * g_RT[1] + g_RT[13]) - (g_RT[0] + g_RT[13]));

    /*reaction 43: H + HO2 <=> O2 + H2 */
    kc[42] = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));

    /*reaction 44: H + HO2 <=> 2 OH */
    kc[43] = exp((g_RT[1] + g_RT[6]) - (2 * g_RT[4]));

    /*reaction 45: H + H2O2 <=> HO2 + H2 */
    kc[44] = exp((g_RT[1] + g_RT[7]) - (g_RT[6] + g_RT[0]));

    /*reaction 46: H + CH4 <=> CH3 + H2 */
    kc[45] = exp((g_RT[1] + g_RT[11]) - (g_RT[10] + g_RT[0]));

    /*reaction 47: H + HCO <=> H2 + CO */
    kc[46] = exp((g_RT[1] + g_RT[14]) - (g_RT[0] + g_RT[12]));

    /*reaction 48: H + CH2O <=> HCO + H2 */
    kc[47] = exp((g_RT[1] + g_RT[15]) - (g_RT[14] + g_RT[0]));

    /*reaction 49: H + CH3O <=> OH + CH3 */
    kc[48] = exp((g_RT[1] + g_RT[16]) - (g_RT[4] + g_RT[10]));

    /*reaction 50: H + C2H3 <=> H2 + C2H2 */
    kc[49] = exp((g_RT[1] + g_RT[18]) - (g_RT[0] + g_RT[17]));

    /*reaction 51: H + C2H4 <=> C2H3 + H2 */
    kc[50] = exp((g_RT[1] + g_RT[19]) - (g_RT[18] + g_RT[0]));

    /*reaction 52: H + C2H6 <=> C2H5 + H2 */
    kc[51] = exp((g_RT[1] + g_RT[21]) - (g_RT[20] + g_RT[0]));

    /*reaction 53: OH + H2 <=> H + H2O */
    kc[52] = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));

    /*reaction 54: 2 OH <=> O + H2O */
    kc[53] = exp((2 * g_RT[4]) - (g_RT[2] + g_RT[5]));

    /*reaction 55: OH + HO2 <=> O2 + H2O */
    kc[54] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));

    /*reaction 56: OH + H2O2 <=> HO2 + H2O */
    kc[55] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));

    /*reaction 57: OH + CH2 <=> H + CH2O */
    kc[56] = exp((g_RT[4] + g_RT[8]) - (g_RT[1] + g_RT[15]));

    /*reaction 58: OH + CH2(S) <=> H + CH2O */
    kc[57] = exp((g_RT[4] + g_RT[9]) - (g_RT[1] + g_RT[15]));

    /*reaction 59: OH + CH3 <=> CH2 + H2O */
    kc[58] = exp((g_RT[4] + g_RT[10]) - (g_RT[8] + g_RT[5]));

    /*reaction 60: OH + CH3 <=> CH2(S) + H2O */
    kc[59] = exp((g_RT[4] + g_RT[10]) - (g_RT[9] + g_RT[5]));

    /*reaction 61: OH + CH4 <=> CH3 + H2O */
    kc[60] = exp((g_RT[4] + g_RT[11]) - (g_RT[10] + g_RT[5]));

    /*reaction 62: OH + CO <=> H + CO2 */
    kc[61] = exp((g_RT[4] + g_RT[12]) - (g_RT[1] + g_RT[13]));

    /*reaction 63: OH + HCO <=> H2O + CO */
    kc[62] = exp((g_RT[4] + g_RT[14]) - (g_RT[5] + g_RT[12]));

    /*reaction 64: OH + CH2O <=> HCO + H2O */
    kc[63] = exp((g_RT[4] + g_RT[15]) - (g_RT[14] + g_RT[5]));

    /*reaction 65: OH + C2H2 <=> CH3 + CO */
    kc[64] = exp((g_RT[4] + g_RT[17]) - (g_RT[10] + g_RT[12]));

    /*reaction 66: OH + C2H3 <=> H2O + C2H2 */
    kc[65] = exp((g_RT[4] + g_RT[18]) - (g_RT[5] + g_RT[17]));

    /*reaction 67: OH + C2H4 <=> C2H3 + H2O */
    kc[66] = exp((g_RT[4] + g_RT[19]) - (g_RT[18] + g_RT[5]));

    /*reaction 68: OH + C2H6 <=> C2H5 + H2O */
    kc[67] = exp((g_RT[4] + g_RT[21]) - (g_RT[20] + g_RT[5]));

    /*reaction 69: 2 HO2 <=> O2 + H2O2 */
    kc[68] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));

    /*reaction 70: 2 HO2 <=> O2 + H2O2 */
    kc[69] = exp((2 * g_RT[6]) - (g_RT[3] + g_RT[7]));

    /*reaction 71: HO2 + CH2 <=> OH + CH2O */
    kc[70] = exp((g_RT[6] + g_RT[8]) - (g_RT[4] + g_RT[15]));

    /*reaction 72: HO2 + CH3 <=> O2 + CH4 */
    kc[71] = exp((g_RT[6] + g_RT[10]) - (g_RT[3] + g_RT[11]));

    /*reaction 73: HO2 + CH3 <=> OH + CH3O */
    kc[72] = exp((g_RT[6] + g_RT[10]) - (g_RT[4] + g_RT[16]));

    /*reaction 74: HO2 + CO <=> OH + CO2 */
    kc[73] = exp((g_RT[6] + g_RT[12]) - (g_RT[4] + g_RT[13]));

    /*reaction 75: HO2 + CH2O <=> HCO + H2O2 */
    kc[74] = exp((g_RT[6] + g_RT[15]) - (g_RT[14] + g_RT[7]));

    /*reaction 76: CH2 + O2 <=> OH + HCO */
    kc[75] = exp((g_RT[8] + g_RT[3]) - (g_RT[4] + g_RT[14]));

    /*reaction 77: CH2 + H2 <=> H + CH3 */
    kc[76] = exp((g_RT[8] + g_RT[0]) - (g_RT[1] + g_RT[10]));

    /*reaction 78: 2 CH2 <=> H2 + C2H2 */
    kc[77] = exp((2 * g_RT[8]) - (g_RT[0] + g_RT[17]));

    /*reaction 79: CH2 + CH3 <=> H + C2H4 */
    kc[78] = exp((g_RT[8] + g_RT[10]) - (g_RT[1] + g_RT[19]));

    /*reaction 80: CH2 + CH4 <=> 2 CH3 */
    kc[79] = exp((g_RT[8] + g_RT[11]) - (2 * g_RT[10]));

    /*reaction 81: CH2(S) + N2 <=> CH2 + N2 */
    kc[80] = exp((g_RT[9] + g_RT[22]) - (g_RT[8] + g_RT[22]));

    /*reaction 82: CH2(S) + AR <=> CH2 + AR */
    kc[81] = exp((g_RT[9] + g_RT[23]) - (g_RT[8] + g_RT[23]));

    /*reaction 83: CH2(S) + O2 <=> H + OH + CO */
    kc[82] = refC * exp((g_RT[9] + g_RT[3]) - (g_RT[1] + g_RT[4] + g_RT[12]));

    /*reaction 84: CH2(S) + O2 <=> CO + H2O */
    kc[83] = exp((g_RT[9] + g_RT[3]) - (g_RT[12] + g_RT[5]));

    /*reaction 85: CH2(S) + H2 <=> CH3 + H */
    kc[84] = exp((g_RT[9] + g_RT[0]) - (g_RT[10] + g_RT[1]));

    /*reaction 86: CH2(S) + H2O <=> CH2 + H2O */
    kc[85] = exp((g_RT[9] + g_RT[5]) - (g_RT[8] + g_RT[5]));

    /*reaction 87: CH2(S) + CH3 <=> H + C2H4 */
    kc[86] = exp((g_RT[9] + g_RT[10]) - (g_RT[1] + g_RT[19]));

    /*reaction 88: CH2(S) + CH4 <=> 2 CH3 */
    kc[87] = exp((g_RT[9] + g_RT[11]) - (2 * g_RT[10]));

    /*reaction 89: CH2(S) + CO <=> CH2 + CO */
    kc[88] = exp((g_RT[9] + g_RT[12]) - (g_RT[8] + g_RT[12]));

    /*reaction 90: CH2(S) + CO2 <=> CH2 + CO2 */
    kc[89] = exp((g_RT[9] + g_RT[13]) - (g_RT[8] + g_RT[13]));

    /*reaction 91: CH2(S) + CO2 <=> CO + CH2O */
    kc[90] = exp((g_RT[9] + g_RT[13]) - (g_RT[12] + g_RT[15]));

    /*reaction 92: CH3 + O2 <=> O + CH3O */
    kc[91] = exp((g_RT[10] + g_RT[3]) - (g_RT[2] + g_RT[16]));

    /*reaction 93: CH3 + O2 <=> OH + CH2O */
    kc[92] = exp((g_RT[10] + g_RT[3]) - (g_RT[4] + g_RT[15]));

    /*reaction 94: CH3 + H2O2 <=> HO2 + CH4 */
    kc[93] = exp((g_RT[10] + g_RT[7]) - (g_RT[6] + g_RT[11]));

    /*reaction 95: 2 CH3 <=> H + C2H5 */
    kc[94] = exp((2 * g_RT[10]) - (g_RT[1] + g_RT[20]));

    /*reaction 96: CH3 + HCO <=> CH4 + CO */
    kc[95] = exp((g_RT[10] + g_RT[14]) - (g_RT[11] + g_RT[12]));

    /*reaction 97: CH3 + CH2O <=> HCO + CH4 */
    kc[96] = exp((g_RT[10] + g_RT[15]) - (g_RT[14] + g_RT[11]));

    /*reaction 98: CH3 + C2H4 <=> C2H3 + CH4 */
    kc[97] = exp((g_RT[10] + g_RT[19]) - (g_RT[18] + g_RT[11]));

    /*reaction 99: CH3 + C2H6 <=> C2H5 + CH4 */
    kc[98] = exp((g_RT[10] + g_RT[21]) - (g_RT[20] + g_RT[11]));

    /*reaction 100: HCO + H2O <=> H + CO + H2O */
    kc[99] = refC * exp((g_RT[14] + g_RT[5]) - (g_RT[1] + g_RT[12] + g_RT[5]));

    /*reaction 101: HCO + O2 <=> HO2 + CO */
    kc[100] = exp((g_RT[14] + g_RT[3]) - (g_RT[6] + g_RT[12]));

    /*reaction 102: CH3O + O2 <=> HO2 + CH2O */
    kc[101] = exp((g_RT[16] + g_RT[3]) - (g_RT[6] + g_RT[15]));

    /*reaction 103: C2H3 + O2 <=> HCO + CH2O */
    kc[102] = exp((g_RT[18] + g_RT[3]) - (g_RT[14] + g_RT[15]));

    /*reaction 104: C2H5 + O2 <=> HO2 + C2H4 */
    kc[103] = exp((g_RT[20] + g_RT[3]) - (g_RT[6] + g_RT[19]));

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
        /*species 8: CH2 */
        species[8] =
            +4.600404010000000e+04 * invT
            +2.200146820000000e+00
            -3.762678670000000e+00 * tc[0]
            -4.844360715000000e-04 * tc[1]
            -4.658164016666667e-07 * tc[2]
            +3.209092941666667e-10 * tc[3]
            -8.437085950000000e-14 * tc[4];
        /*species 9: CH2(S) */
        species[9] =
            +5.049681630000000e+04 * invT
            +4.967723077000000e+00
            -4.198604110000000e+00 * tc[0]
            +1.183307095000000e-03 * tc[1]
            -1.372160366666667e-06 * tc[2]
            +5.573466508333334e-10 * tc[3]
            -9.715736850000000e-14 * tc[4];
        /*species 10: CH3 */
        species[10] =
            +1.644499880000000e+04 * invT
            +2.069026070000000e+00
            -3.673590400000000e+00 * tc[0]
            -1.005475875000000e-03 * tc[1]
            -9.550364266666668e-07 * tc[2]
            +5.725978541666666e-10 * tc[3]
            -1.271928670000000e-13 * tc[4];
        /*species 11: CH4 */
        species[11] =
            -1.024664760000000e+04 * invT
            +9.791179889999999e+00
            -5.149876130000000e+00 * tc[0]
            +6.835489400000000e-03 * tc[1]
            -8.196676650000000e-06 * tc[2]
            +4.039525216666667e-09 * tc[3]
            -8.334697800000000e-13 * tc[4];
        /*species 12: CO */
        species[12] =
            -1.434408600000000e+04 * invT
            +7.112418999999992e-02
            -3.579533470000000e+00 * tc[0]
            +3.051768400000000e-04 * tc[1]
            -1.694690550000000e-07 * tc[2]
            -7.558382366666667e-11 * tc[3]
            +4.522122495000000e-14 * tc[4];
        /*species 13: CO2 */
        species[13] =
            -4.837196970000000e+04 * invT
            -7.544278700000000e+00
            -2.356773520000000e+00 * tc[0]
            -4.492298385000000e-03 * tc[1]
            +1.187260448333333e-06 * tc[2]
            -2.049325183333333e-10 * tc[3]
            +7.184977399999999e-15 * tc[4];
        /*species 14: HCO */
        species[14] =
            +3.839564960000000e+03 * invT
            +8.268134100000002e-01
            -4.221185840000000e+00 * tc[0]
            +1.621962660000000e-03 * tc[1]
            -2.296657433333333e-06 * tc[2]
            +1.109534108333333e-09 * tc[3]
            -2.168844325000000e-13 * tc[4];
        /*species 15: CH2O */
        species[15] =
            -1.430895670000000e+04 * invT
            +4.190910250000000e+00
            -4.793723150000000e+00 * tc[0]
            +4.954166845000000e-03 * tc[1]
            -6.220333466666666e-06 * tc[2]
            +3.160710508333333e-09 * tc[3]
            -6.588632600000000e-13 * tc[4];
        /*species 16: CH3O */
        species[16] =
            +9.786011000000000e+02 * invT
            -1.104597300000000e+01
            -2.106204000000000e+00 * tc[0]
            -3.608297500000000e-03 * tc[1]
            -8.897453333333333e-07 * tc[2]
            +6.148030000000000e-10 * tc[3]
            -1.037805000000000e-13 * tc[4];
        /*species 17: C2H2 */
        species[17] =
            +2.642898070000000e+04 * invT
            -1.313102400600000e+01
            -8.086810940000000e-01 * tc[0]
            -1.168078145000000e-02 * tc[1]
            +5.919530250000000e-06 * tc[2]
            -2.334603641666667e-09 * tc[3]
            +4.250364870000000e-13 * tc[4];
        /*species 18: C2H3 */
        species[18] =
            +3.485984680000000e+04 * invT
            -5.298073800000000e+00
            -3.212466450000000e+00 * tc[0]
            -7.573958100000000e-04 * tc[1]
            -4.320156866666666e-06 * tc[2]
            +2.980482058333333e-09 * tc[3]
            -7.357543650000000e-13 * tc[4];
        /*species 19: C2H4 */
        species[19] =
            +5.089775930000000e+03 * invT
            -1.381294799999999e-01
            -3.959201480000000e+00 * tc[0]
            +3.785261235000000e-03 * tc[1]
            -9.516504866666667e-06 * tc[2]
            +5.763239608333333e-09 * tc[3]
            -1.349421865000000e-12 * tc[4];
        /*species 20: C2H5 */
        species[20] =
            +1.284162650000000e+04 * invT
            -4.007435600000004e-01
            -4.306465680000000e+00 * tc[0]
            +2.093294460000000e-03 * tc[1]
            -8.285713450000000e-06 * tc[2]
            +4.992721716666666e-09 * tc[3]
            -1.152545020000000e-12 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            -1.152220550000000e+04 * invT
            +1.624601760000000e+00
            -4.291424920000000e+00 * tc[0]
            +2.750771350000000e-03 * tc[1]
            -9.990638133333334e-06 * tc[2]
            +5.903885708333334e-09 * tc[3]
            -1.343428855000000e-12 * tc[4];
        /*species 22: N2 */
        species[22] =
            -1.020899900000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041202000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427000000000e-13 * tc[4];
        /*species 23: AR */
        species[23] =
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
        /*species 8: CH2 */
        species[8] =
            +4.626360400000000e+04 * invT
            -3.297092110000000e+00
            -2.874101130000000e+00 * tc[0]
            -1.828196460000000e-03 * tc[1]
            +2.348243283333333e-07 * tc[2]
            -2.168162908333333e-11 * tc[3]
            +9.386378350000000e-16 * tc[4];
        /*species 9: CH2(S) */
        species[9] =
            +5.092599970000000e+04 * invT
            -6.334463270000000e+00
            -2.292038420000000e+00 * tc[0]
            -2.327943185000000e-03 * tc[1]
            +3.353199116666667e-07 * tc[2]
            -3.482550000000000e-11 * tc[3]
            +1.698581825000000e-15 * tc[4];
        /*species 10: CH3 */
        species[10] =
            +1.677558430000000e+04 * invT
            -6.194354070000000e+00
            -2.285717720000000e+00 * tc[0]
            -3.619950185000000e-03 * tc[1]
            +4.978572466666667e-07 * tc[2]
            -4.964038700000000e-11 * tc[3]
            +2.335771970000000e-15 * tc[4];
        /*species 11: CH4 */
        species[11] =
            -9.468344590000001e+03 * invT
            -1.836246650500000e+01
            -7.485149500000000e-02 * tc[0]
            -6.695473350000000e-03 * tc[1]
            +9.554763483333333e-07 * tc[2]
            -1.019104458333333e-10 * tc[3]
            +5.090761500000000e-15 * tc[4];
        /*species 12: CO */
        species[12] =
            -1.415187240000000e+04 * invT
            -5.103502110000000e+00
            -2.715185610000000e+00 * tc[0]
            -1.031263715000000e-03 * tc[1]
            +1.664709618333334e-07 * tc[2]
            -1.917108400000000e-11 * tc[3]
            +1.018238580000000e-15 * tc[4];
        /*species 13: CO2 */
        species[13] =
            -4.875916600000000e+04 * invT
            +1.585822230000000e+00
            -3.857460290000000e+00 * tc[0]
            -2.207185130000000e-03 * tc[1]
            +3.691356733333334e-07 * tc[2]
            -4.362418233333334e-11 * tc[3]
            +2.360420820000000e-15 * tc[4];
        /*species 14: HCO */
        species[14] =
            +4.011918150000000e+03 * invT
            -7.026170540000000e+00
            -2.772174380000000e+00 * tc[0]
            -2.478477630000000e-03 * tc[1]
            +4.140760216666667e-07 * tc[2]
            -4.909681483333334e-11 * tc[3]
            +2.667543555000000e-15 * tc[4];
        /*species 15: CH2O */
        species[15] =
            -1.399583230000000e+04 * invT
            -1.189563292000000e+01
            -1.760690080000000e+00 * tc[0]
            -4.600000410000000e-03 * tc[1]
            +7.370980216666666e-07 * tc[2]
            -8.386767666666666e-11 * tc[3]
            +4.419278200000001e-15 * tc[4];
        /*species 16: CH3O */
        species[16] =
            +1.278325200000000e+02 * invT
            +8.412240000000000e-01
            -3.770799000000000e+00 * tc[0]
            -3.935748500000000e-03 * tc[1]
            +4.427306666666667e-07 * tc[2]
            -3.287025833333333e-11 * tc[3]
            +1.056308000000000e-15 * tc[4];
        /*species 17: C2H2 */
        species[17] =
            +2.593599920000000e+04 * invT
            +5.377850850000001e+00
            -4.147569640000000e+00 * tc[0]
            -2.980833320000000e-03 * tc[1]
            +3.954914200000000e-07 * tc[2]
            -3.895101425000000e-11 * tc[3]
            +1.806176065000000e-15 * tc[4];
        /*species 18: C2H3 */
        species[18] =
            +3.461287390000000e+04 * invT
            -4.770599780000000e+00
            -3.016724000000000e+00 * tc[0]
            -5.165114600000000e-03 * tc[1]
            +7.801372483333333e-07 * tc[2]
            -8.480274000000000e-11 * tc[3]
            +4.313035205000000e-15 * tc[4];
        /*species 19: C2H4 */
        species[19] =
            +4.939886140000000e+03 * invT
            -8.269258140000002e+00
            -2.036111160000000e+00 * tc[0]
            -7.322707550000000e-03 * tc[1]
            +1.118463191666667e-06 * tc[2]
            -1.226857691666667e-10 * tc[3]
            +6.285303050000000e-15 * tc[4];
        /*species 20: C2H5 */
        species[20] =
            +1.285752000000000e+04 * invT
            -1.150777788000000e+01
            -1.954656420000000e+00 * tc[0]
            -8.698636100000001e-03 * tc[1]
            +1.330344446666667e-06 * tc[2]
            -1.460147408333333e-10 * tc[3]
            +7.482078800000000e-15 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            -1.142639320000000e+04 * invT
            -1.404372920000000e+01
            -1.071881500000000e+00 * tc[0]
            -1.084263385000000e-02 * tc[1]
            +1.670934450000000e-06 * tc[2]
            -1.845100008333333e-10 * tc[3]
            +9.500144500000000e-15 * tc[4];
        /*species 22: N2 */
        species[22] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439884000000000e-04 * tc[1]
            +9.474600000000001e-08 * tc[2]
            -8.414198333333333e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
        /*species 23: AR */
        species[23] =
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
        /*species 8: CH2 */
        species[8] =
            +4.60040401e+04 * invT
            +1.20014682e+00
            -3.76267867e+00 * tc[0]
            -4.84436072e-04 * tc[1]
            -4.65816402e-07 * tc[2]
            +3.20909294e-10 * tc[3]
            -8.43708595e-14 * tc[4];
        /*species 9: CH2(S) */
        species[9] =
            +5.04968163e+04 * invT
            +3.96772308e+00
            -4.19860411e+00 * tc[0]
            +1.18330710e-03 * tc[1]
            -1.37216037e-06 * tc[2]
            +5.57346651e-10 * tc[3]
            -9.71573685e-14 * tc[4];
        /*species 10: CH3 */
        species[10] =
            +1.64449988e+04 * invT
            +1.06902607e+00
            -3.67359040e+00 * tc[0]
            -1.00547588e-03 * tc[1]
            -9.55036427e-07 * tc[2]
            +5.72597854e-10 * tc[3]
            -1.27192867e-13 * tc[4];
        /*species 11: CH4 */
        species[11] =
            -1.02466476e+04 * invT
            +8.79117989e+00
            -5.14987613e+00 * tc[0]
            +6.83548940e-03 * tc[1]
            -8.19667665e-06 * tc[2]
            +4.03952522e-09 * tc[3]
            -8.33469780e-13 * tc[4];
        /*species 12: CO */
        species[12] =
            -1.43440860e+04 * invT
            -9.28875810e-01
            -3.57953347e+00 * tc[0]
            +3.05176840e-04 * tc[1]
            -1.69469055e-07 * tc[2]
            -7.55838237e-11 * tc[3]
            +4.52212249e-14 * tc[4];
        /*species 13: CO2 */
        species[13] =
            -4.83719697e+04 * invT
            -8.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
        /*species 14: HCO */
        species[14] =
            +3.83956496e+03 * invT
            -1.73186590e-01
            -4.22118584e+00 * tc[0]
            +1.62196266e-03 * tc[1]
            -2.29665743e-06 * tc[2]
            +1.10953411e-09 * tc[3]
            -2.16884432e-13 * tc[4];
        /*species 15: CH2O */
        species[15] =
            -1.43089567e+04 * invT
            +3.19091025e+00
            -4.79372315e+00 * tc[0]
            +4.95416684e-03 * tc[1]
            -6.22033347e-06 * tc[2]
            +3.16071051e-09 * tc[3]
            -6.58863260e-13 * tc[4];
        /*species 16: CH3O */
        species[16] =
            +9.78601100e+02 * invT
            -1.20459730e+01
            -2.10620400e+00 * tc[0]
            -3.60829750e-03 * tc[1]
            -8.89745333e-07 * tc[2]
            +6.14803000e-10 * tc[3]
            -1.03780500e-13 * tc[4];
        /*species 17: C2H2 */
        species[17] =
            +2.64289807e+04 * invT
            -1.41310240e+01
            -8.08681094e-01 * tc[0]
            -1.16807815e-02 * tc[1]
            +5.91953025e-06 * tc[2]
            -2.33460364e-09 * tc[3]
            +4.25036487e-13 * tc[4];
        /*species 18: C2H3 */
        species[18] =
            +3.48598468e+04 * invT
            -6.29807380e+00
            -3.21246645e+00 * tc[0]
            -7.57395810e-04 * tc[1]
            -4.32015687e-06 * tc[2]
            +2.98048206e-09 * tc[3]
            -7.35754365e-13 * tc[4];
        /*species 19: C2H4 */
        species[19] =
            +5.08977593e+03 * invT
            -1.13812948e+00
            -3.95920148e+00 * tc[0]
            +3.78526124e-03 * tc[1]
            -9.51650487e-06 * tc[2]
            +5.76323961e-09 * tc[3]
            -1.34942187e-12 * tc[4];
        /*species 20: C2H5 */
        species[20] =
            +1.28416265e+04 * invT
            -1.40074356e+00
            -4.30646568e+00 * tc[0]
            +2.09329446e-03 * tc[1]
            -8.28571345e-06 * tc[2]
            +4.99272172e-09 * tc[3]
            -1.15254502e-12 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            -1.15222055e+04 * invT
            +6.24601760e-01
            -4.29142492e+00 * tc[0]
            +2.75077135e-03 * tc[1]
            -9.99063813e-06 * tc[2]
            +5.90388571e-09 * tc[3]
            -1.34342886e-12 * tc[4];
        /*species 22: N2 */
        species[22] =
            -1.02089990e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
        /*species 23: AR */
        species[23] =
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
        /*species 8: CH2 */
        species[8] =
            +4.62636040e+04 * invT
            -4.29709211e+00
            -2.87410113e+00 * tc[0]
            -1.82819646e-03 * tc[1]
            +2.34824328e-07 * tc[2]
            -2.16816291e-11 * tc[3]
            +9.38637835e-16 * tc[4];
        /*species 9: CH2(S) */
        species[9] =
            +5.09259997e+04 * invT
            -7.33446327e+00
            -2.29203842e+00 * tc[0]
            -2.32794318e-03 * tc[1]
            +3.35319912e-07 * tc[2]
            -3.48255000e-11 * tc[3]
            +1.69858182e-15 * tc[4];
        /*species 10: CH3 */
        species[10] =
            +1.67755843e+04 * invT
            -7.19435407e+00
            -2.28571772e+00 * tc[0]
            -3.61995018e-03 * tc[1]
            +4.97857247e-07 * tc[2]
            -4.96403870e-11 * tc[3]
            +2.33577197e-15 * tc[4];
        /*species 11: CH4 */
        species[11] =
            -9.46834459e+03 * invT
            -1.93624665e+01
            -7.48514950e-02 * tc[0]
            -6.69547335e-03 * tc[1]
            +9.55476348e-07 * tc[2]
            -1.01910446e-10 * tc[3]
            +5.09076150e-15 * tc[4];
        /*species 12: CO */
        species[12] =
            -1.41518724e+04 * invT
            -6.10350211e+00
            -2.71518561e+00 * tc[0]
            -1.03126372e-03 * tc[1]
            +1.66470962e-07 * tc[2]
            -1.91710840e-11 * tc[3]
            +1.01823858e-15 * tc[4];
        /*species 13: CO2 */
        species[13] =
            -4.87591660e+04 * invT
            +5.85822230e-01
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
        /*species 14: HCO */
        species[14] =
            +4.01191815e+03 * invT
            -8.02617054e+00
            -2.77217438e+00 * tc[0]
            -2.47847763e-03 * tc[1]
            +4.14076022e-07 * tc[2]
            -4.90968148e-11 * tc[3]
            +2.66754356e-15 * tc[4];
        /*species 15: CH2O */
        species[15] =
            -1.39958323e+04 * invT
            -1.28956329e+01
            -1.76069008e+00 * tc[0]
            -4.60000041e-03 * tc[1]
            +7.37098022e-07 * tc[2]
            -8.38676767e-11 * tc[3]
            +4.41927820e-15 * tc[4];
        /*species 16: CH3O */
        species[16] =
            +1.27832520e+02 * invT
            -1.58776000e-01
            -3.77079900e+00 * tc[0]
            -3.93574850e-03 * tc[1]
            +4.42730667e-07 * tc[2]
            -3.28702583e-11 * tc[3]
            +1.05630800e-15 * tc[4];
        /*species 17: C2H2 */
        species[17] =
            +2.59359992e+04 * invT
            +4.37785085e+00
            -4.14756964e+00 * tc[0]
            -2.98083332e-03 * tc[1]
            +3.95491420e-07 * tc[2]
            -3.89510143e-11 * tc[3]
            +1.80617607e-15 * tc[4];
        /*species 18: C2H3 */
        species[18] =
            +3.46128739e+04 * invT
            -5.77059978e+00
            -3.01672400e+00 * tc[0]
            -5.16511460e-03 * tc[1]
            +7.80137248e-07 * tc[2]
            -8.48027400e-11 * tc[3]
            +4.31303520e-15 * tc[4];
        /*species 19: C2H4 */
        species[19] =
            +4.93988614e+03 * invT
            -9.26925814e+00
            -2.03611116e+00 * tc[0]
            -7.32270755e-03 * tc[1]
            +1.11846319e-06 * tc[2]
            -1.22685769e-10 * tc[3]
            +6.28530305e-15 * tc[4];
        /*species 20: C2H5 */
        species[20] =
            +1.28575200e+04 * invT
            -1.25077779e+01
            -1.95465642e+00 * tc[0]
            -8.69863610e-03 * tc[1]
            +1.33034445e-06 * tc[2]
            -1.46014741e-10 * tc[3]
            +7.48207880e-15 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            -1.14263932e+04 * invT
            -1.50437292e+01
            -1.07188150e+00 * tc[0]
            -1.08426339e-02 * tc[1]
            +1.67093445e-06 * tc[2]
            -1.84510001e-10 * tc[3]
            +9.50014450e-15 * tc[4];
        /*species 22: N2 */
        species[22] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 23: AR */
        species[23] =
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
        /*species 8: CH2 */
        species[8] =
            +2.76267867e+00
            +9.68872143e-04 * tc[1]
            +2.79489841e-06 * tc[2]
            -3.85091153e-09 * tc[3]
            +1.68741719e-12 * tc[4];
        /*species 9: CH2(S) */
        species[9] =
            +3.19860411e+00
            -2.36661419e-03 * tc[1]
            +8.23296220e-06 * tc[2]
            -6.68815981e-09 * tc[3]
            +1.94314737e-12 * tc[4];
        /*species 10: CH3 */
        species[10] =
            +2.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 11: CH4 */
        species[11] =
            +4.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 12: CO */
        species[12] =
            +2.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 13: CO2 */
        species[13] =
            +1.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 14: HCO */
        species[14] =
            +3.22118584e+00
            -3.24392532e-03 * tc[1]
            +1.37799446e-05 * tc[2]
            -1.33144093e-08 * tc[3]
            +4.33768865e-12 * tc[4];
        /*species 15: CH2O */
        species[15] =
            +3.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 16: CH3O */
        species[16] =
            +1.10620400e+00
            +7.21659500e-03 * tc[1]
            +5.33847200e-06 * tc[2]
            -7.37763600e-09 * tc[3]
            +2.07561000e-12 * tc[4];
        /*species 17: C2H2 */
        species[17] =
            -1.91318906e-01
            +2.33615629e-02 * tc[1]
            -3.55171815e-05 * tc[2]
            +2.80152437e-08 * tc[3]
            -8.50072974e-12 * tc[4];
        /*species 18: C2H3 */
        species[18] =
            +2.21246645e+00
            +1.51479162e-03 * tc[1]
            +2.59209412e-05 * tc[2]
            -3.57657847e-08 * tc[3]
            +1.47150873e-11 * tc[4];
        /*species 19: C2H4 */
        species[19] =
            +2.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 20: C2H5 */
        species[20] =
            +3.30646568e+00
            -4.18658892e-03 * tc[1]
            +4.97142807e-05 * tc[2]
            -5.99126606e-08 * tc[3]
            +2.30509004e-11 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            +3.29142492e+00
            -5.50154270e-03 * tc[1]
            +5.99438288e-05 * tc[2]
            -7.08466285e-08 * tc[3]
            +2.68685771e-11 * tc[4];
        /*species 22: N2 */
        species[22] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 23: AR */
        species[23] =
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
        /*species 8: CH2 */
        species[8] =
            +1.87410113e+00
            +3.65639292e-03 * tc[1]
            -1.40894597e-06 * tc[2]
            +2.60179549e-10 * tc[3]
            -1.87727567e-14 * tc[4];
        /*species 9: CH2(S) */
        species[9] =
            +1.29203842e+00
            +4.65588637e-03 * tc[1]
            -2.01191947e-06 * tc[2]
            +4.17906000e-10 * tc[3]
            -3.39716365e-14 * tc[4];
        /*species 10: CH3 */
        species[10] =
            +1.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 11: CH4 */
        species[11] =
            -9.25148505e-01
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 12: CO */
        species[12] =
            +1.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 13: CO2 */
        species[13] =
            +2.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 14: HCO */
        species[14] =
            +1.77217438e+00
            +4.95695526e-03 * tc[1]
            -2.48445613e-06 * tc[2]
            +5.89161778e-10 * tc[3]
            -5.33508711e-14 * tc[4];
        /*species 15: CH2O */
        species[15] =
            +7.60690080e-01
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 16: CH3O */
        species[16] =
            +2.77079900e+00
            +7.87149700e-03 * tc[1]
            -2.65638400e-06 * tc[2]
            +3.94443100e-10 * tc[3]
            -2.11261600e-14 * tc[4];
        /*species 17: C2H2 */
        species[17] =
            +3.14756964e+00
            +5.96166664e-03 * tc[1]
            -2.37294852e-06 * tc[2]
            +4.67412171e-10 * tc[3]
            -3.61235213e-14 * tc[4];
        /*species 18: C2H3 */
        species[18] =
            +2.01672400e+00
            +1.03302292e-02 * tc[1]
            -4.68082349e-06 * tc[2]
            +1.01763288e-09 * tc[3]
            -8.62607041e-14 * tc[4];
        /*species 19: C2H4 */
        species[19] =
            +1.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 20: C2H5 */
        species[20] =
            +9.54656420e-01
            +1.73972722e-02 * tc[1]
            -7.98206668e-06 * tc[2]
            +1.75217689e-09 * tc[3]
            -1.49641576e-13 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            +7.18815000e-02
            +2.16852677e-02 * tc[1]
            -1.00256067e-05 * tc[2]
            +2.21412001e-09 * tc[3]
            -1.90002890e-13 * tc[4];
        /*species 22: N2 */
        species[22] =
            +1.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 23: AR */
        species[23] =
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
        /*species 8: CH2 */
        species[8] =
            +3.76267867e+00
            +9.68872143e-04 * tc[1]
            +2.79489841e-06 * tc[2]
            -3.85091153e-09 * tc[3]
            +1.68741719e-12 * tc[4];
        /*species 9: CH2(S) */
        species[9] =
            +4.19860411e+00
            -2.36661419e-03 * tc[1]
            +8.23296220e-06 * tc[2]
            -6.68815981e-09 * tc[3]
            +1.94314737e-12 * tc[4];
        /*species 10: CH3 */
        species[10] =
            +3.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 11: CH4 */
        species[11] =
            +5.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 12: CO */
        species[12] =
            +3.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 13: CO2 */
        species[13] =
            +2.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 14: HCO */
        species[14] =
            +4.22118584e+00
            -3.24392532e-03 * tc[1]
            +1.37799446e-05 * tc[2]
            -1.33144093e-08 * tc[3]
            +4.33768865e-12 * tc[4];
        /*species 15: CH2O */
        species[15] =
            +4.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 16: CH3O */
        species[16] =
            +2.10620400e+00
            +7.21659500e-03 * tc[1]
            +5.33847200e-06 * tc[2]
            -7.37763600e-09 * tc[3]
            +2.07561000e-12 * tc[4];
        /*species 17: C2H2 */
        species[17] =
            +8.08681094e-01
            +2.33615629e-02 * tc[1]
            -3.55171815e-05 * tc[2]
            +2.80152437e-08 * tc[3]
            -8.50072974e-12 * tc[4];
        /*species 18: C2H3 */
        species[18] =
            +3.21246645e+00
            +1.51479162e-03 * tc[1]
            +2.59209412e-05 * tc[2]
            -3.57657847e-08 * tc[3]
            +1.47150873e-11 * tc[4];
        /*species 19: C2H4 */
        species[19] =
            +3.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 20: C2H5 */
        species[20] =
            +4.30646568e+00
            -4.18658892e-03 * tc[1]
            +4.97142807e-05 * tc[2]
            -5.99126606e-08 * tc[3]
            +2.30509004e-11 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            +4.29142492e+00
            -5.50154270e-03 * tc[1]
            +5.99438288e-05 * tc[2]
            -7.08466285e-08 * tc[3]
            +2.68685771e-11 * tc[4];
        /*species 22: N2 */
        species[22] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 23: AR */
        species[23] =
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
        /*species 8: CH2 */
        species[8] =
            +2.87410113e+00
            +3.65639292e-03 * tc[1]
            -1.40894597e-06 * tc[2]
            +2.60179549e-10 * tc[3]
            -1.87727567e-14 * tc[4];
        /*species 9: CH2(S) */
        species[9] =
            +2.29203842e+00
            +4.65588637e-03 * tc[1]
            -2.01191947e-06 * tc[2]
            +4.17906000e-10 * tc[3]
            -3.39716365e-14 * tc[4];
        /*species 10: CH3 */
        species[10] =
            +2.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 11: CH4 */
        species[11] =
            +7.48514950e-02
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 12: CO */
        species[12] =
            +2.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 13: CO2 */
        species[13] =
            +3.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 14: HCO */
        species[14] =
            +2.77217438e+00
            +4.95695526e-03 * tc[1]
            -2.48445613e-06 * tc[2]
            +5.89161778e-10 * tc[3]
            -5.33508711e-14 * tc[4];
        /*species 15: CH2O */
        species[15] =
            +1.76069008e+00
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 16: CH3O */
        species[16] =
            +3.77079900e+00
            +7.87149700e-03 * tc[1]
            -2.65638400e-06 * tc[2]
            +3.94443100e-10 * tc[3]
            -2.11261600e-14 * tc[4];
        /*species 17: C2H2 */
        species[17] =
            +4.14756964e+00
            +5.96166664e-03 * tc[1]
            -2.37294852e-06 * tc[2]
            +4.67412171e-10 * tc[3]
            -3.61235213e-14 * tc[4];
        /*species 18: C2H3 */
        species[18] =
            +3.01672400e+00
            +1.03302292e-02 * tc[1]
            -4.68082349e-06 * tc[2]
            +1.01763288e-09 * tc[3]
            -8.62607041e-14 * tc[4];
        /*species 19: C2H4 */
        species[19] =
            +2.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 20: C2H5 */
        species[20] =
            +1.95465642e+00
            +1.73972722e-02 * tc[1]
            -7.98206668e-06 * tc[2]
            +1.75217689e-09 * tc[3]
            -1.49641576e-13 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            +1.07188150e+00
            +2.16852677e-02 * tc[1]
            -1.00256067e-05 * tc[2]
            +2.21412001e-09 * tc[3]
            -1.90002890e-13 * tc[4];
        /*species 22: N2 */
        species[22] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 23: AR */
        species[23] =
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
        /*species 8: CH2 */
        species[8] =
            +2.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 * invT;
        /*species 9: CH2(S) */
        species[9] =
            +3.19860411e+00
            -1.18330710e-03 * tc[1]
            +2.74432073e-06 * tc[2]
            -1.67203995e-09 * tc[3]
            +3.88629474e-13 * tc[4]
            +5.04968163e+04 * invT;
        /*species 10: CH3 */
        species[10] =
            +2.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 * invT;
        /*species 11: CH4 */
        species[11] =
            +4.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 12: CO */
        species[12] =
            +2.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 13: CO2 */
        species[13] =
            +1.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 14: HCO */
        species[14] =
            +3.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 * invT;
        /*species 15: CH2O */
        species[15] =
            +3.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 * invT;
        /*species 16: CH3O */
        species[16] =
            +1.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122000e-13 * tc[4]
            +9.78601100e+02 * invT;
        /*species 17: C2H2 */
        species[17] =
            -1.91318906e-01
            +1.16807815e-02 * tc[1]
            -1.18390605e-05 * tc[2]
            +7.00381092e-09 * tc[3]
            -1.70014595e-12 * tc[4]
            +2.64289807e+04 * invT;
        /*species 18: C2H3 */
        species[18] =
            +2.21246645e+00
            +7.57395810e-04 * tc[1]
            +8.64031373e-06 * tc[2]
            -8.94144617e-09 * tc[3]
            +2.94301746e-12 * tc[4]
            +3.48598468e+04 * invT;
        /*species 19: C2H4 */
        species[19] =
            +2.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 * invT;
        /*species 20: C2H5 */
        species[20] =
            +3.30646568e+00
            -2.09329446e-03 * tc[1]
            +1.65714269e-05 * tc[2]
            -1.49781651e-08 * tc[3]
            +4.61018008e-12 * tc[4]
            +1.28416265e+04 * invT;
        /*species 21: C2H6 */
        species[21] =
            +3.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 * invT;
        /*species 22: N2 */
        species[22] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
        /*species 23: AR */
        species[23] =
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
        /*species 8: CH2 */
        species[8] =
            +1.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 * invT;
        /*species 9: CH2(S) */
        species[9] =
            +1.29203842e+00
            +2.32794318e-03 * tc[1]
            -6.70639823e-07 * tc[2]
            +1.04476500e-10 * tc[3]
            -6.79432730e-15 * tc[4]
            +5.09259997e+04 * invT;
        /*species 10: CH3 */
        species[10] =
            +1.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 * invT;
        /*species 11: CH4 */
        species[11] =
            -9.25148505e-01
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 12: CO */
        species[12] =
            +1.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 13: CO2 */
        species[13] =
            +2.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 14: HCO */
        species[14] =
            +1.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 * invT;
        /*species 15: CH2O */
        species[15] =
            +7.60690080e-01
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 * invT;
        /*species 16: CH3O */
        species[16] =
            +2.77079900e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832520e+02 * invT;
        /*species 17: C2H2 */
        species[17] =
            +3.14756964e+00
            +2.98083332e-03 * tc[1]
            -7.90982840e-07 * tc[2]
            +1.16853043e-10 * tc[3]
            -7.22470426e-15 * tc[4]
            +2.59359992e+04 * invT;
        /*species 18: C2H3 */
        species[18] =
            +2.01672400e+00
            +5.16511460e-03 * tc[1]
            -1.56027450e-06 * tc[2]
            +2.54408220e-10 * tc[3]
            -1.72521408e-14 * tc[4]
            +3.46128739e+04 * invT;
        /*species 19: C2H4 */
        species[19] =
            +1.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 * invT;
        /*species 20: C2H5 */
        species[20] =
            +9.54656420e-01
            +8.69863610e-03 * tc[1]
            -2.66068889e-06 * tc[2]
            +4.38044223e-10 * tc[3]
            -2.99283152e-14 * tc[4]
            +1.28575200e+04 * invT;
        /*species 21: C2H6 */
        species[21] =
            +7.18815000e-02
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 * invT;
        /*species 22: N2 */
        species[22] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 23: AR */
        species[23] =
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
        /*species 8: CH2 */
        species[8] =
            +3.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 * invT;
        /*species 9: CH2(S) */
        species[9] =
            +4.19860411e+00
            -1.18330710e-03 * tc[1]
            +2.74432073e-06 * tc[2]
            -1.67203995e-09 * tc[3]
            +3.88629474e-13 * tc[4]
            +5.04968163e+04 * invT;
        /*species 10: CH3 */
        species[10] =
            +3.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 * invT;
        /*species 11: CH4 */
        species[11] =
            +5.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 12: CO */
        species[12] =
            +3.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 13: CO2 */
        species[13] =
            +2.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 14: HCO */
        species[14] =
            +4.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 * invT;
        /*species 15: CH2O */
        species[15] =
            +4.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 * invT;
        /*species 16: CH3O */
        species[16] =
            +2.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122000e-13 * tc[4]
            +9.78601100e+02 * invT;
        /*species 17: C2H2 */
        species[17] =
            +8.08681094e-01
            +1.16807815e-02 * tc[1]
            -1.18390605e-05 * tc[2]
            +7.00381092e-09 * tc[3]
            -1.70014595e-12 * tc[4]
            +2.64289807e+04 * invT;
        /*species 18: C2H3 */
        species[18] =
            +3.21246645e+00
            +7.57395810e-04 * tc[1]
            +8.64031373e-06 * tc[2]
            -8.94144617e-09 * tc[3]
            +2.94301746e-12 * tc[4]
            +3.48598468e+04 * invT;
        /*species 19: C2H4 */
        species[19] =
            +3.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 * invT;
        /*species 20: C2H5 */
        species[20] =
            +4.30646568e+00
            -2.09329446e-03 * tc[1]
            +1.65714269e-05 * tc[2]
            -1.49781651e-08 * tc[3]
            +4.61018008e-12 * tc[4]
            +1.28416265e+04 * invT;
        /*species 21: C2H6 */
        species[21] =
            +4.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 * invT;
        /*species 22: N2 */
        species[22] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
        /*species 23: AR */
        species[23] =
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
        /*species 8: CH2 */
        species[8] =
            +2.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 * invT;
        /*species 9: CH2(S) */
        species[9] =
            +2.29203842e+00
            +2.32794318e-03 * tc[1]
            -6.70639823e-07 * tc[2]
            +1.04476500e-10 * tc[3]
            -6.79432730e-15 * tc[4]
            +5.09259997e+04 * invT;
        /*species 10: CH3 */
        species[10] =
            +2.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 * invT;
        /*species 11: CH4 */
        species[11] =
            +7.48514950e-02
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 12: CO */
        species[12] =
            +2.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 13: CO2 */
        species[13] =
            +3.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 14: HCO */
        species[14] =
            +2.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 * invT;
        /*species 15: CH2O */
        species[15] =
            +1.76069008e+00
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 * invT;
        /*species 16: CH3O */
        species[16] =
            +3.77079900e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832520e+02 * invT;
        /*species 17: C2H2 */
        species[17] =
            +4.14756964e+00
            +2.98083332e-03 * tc[1]
            -7.90982840e-07 * tc[2]
            +1.16853043e-10 * tc[3]
            -7.22470426e-15 * tc[4]
            +2.59359992e+04 * invT;
        /*species 18: C2H3 */
        species[18] =
            +3.01672400e+00
            +5.16511460e-03 * tc[1]
            -1.56027450e-06 * tc[2]
            +2.54408220e-10 * tc[3]
            -1.72521408e-14 * tc[4]
            +3.46128739e+04 * invT;
        /*species 19: C2H4 */
        species[19] =
            +2.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 * invT;
        /*species 20: C2H5 */
        species[20] =
            +1.95465642e+00
            +8.69863610e-03 * tc[1]
            -2.66068889e-06 * tc[2]
            +4.38044223e-10 * tc[3]
            -2.99283152e-14 * tc[4]
            +1.28575200e+04 * invT;
        /*species 21: C2H6 */
        species[21] =
            +1.07188150e+00
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 * invT;
        /*species 22: N2 */
        species[22] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 23: AR */
        species[23] =
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
        /*species 8: CH2 */
        species[8] =
            +3.76267867e+00 * tc[0]
            +9.68872143e-04 * tc[1]
            +1.39744921e-06 * tc[2]
            -1.28363718e-09 * tc[3]
            +4.21854298e-13 * tc[4]
            +1.56253185e+00 ;
        /*species 9: CH2(S) */
        species[9] =
            +4.19860411e+00 * tc[0]
            -2.36661419e-03 * tc[1]
            +4.11648110e-06 * tc[2]
            -2.22938660e-09 * tc[3]
            +4.85786843e-13 * tc[4]
            -7.69118967e-01 ;
        /*species 10: CH3 */
        species[10] =
            +3.67359040e+00 * tc[0]
            +2.01095175e-03 * tc[1]
            +2.86510928e-06 * tc[2]
            -2.29039142e-09 * tc[3]
            +6.35964335e-13 * tc[4]
            +1.60456433e+00 ;
        /*species 11: CH4 */
        species[11] =
            +5.14987613e+00 * tc[0]
            -1.36709788e-02 * tc[1]
            +2.45900299e-05 * tc[2]
            -1.61581009e-08 * tc[3]
            +4.16734890e-12 * tc[4]
            -4.64130376e+00 ;
        /*species 12: CO */
        species[12] =
            +3.57953347e+00 * tc[0]
            -6.10353680e-04 * tc[1]
            +5.08407165e-07 * tc[2]
            +3.02335295e-10 * tc[3]
            -2.26106125e-13 * tc[4]
            +3.50840928e+00 ;
        /*species 13: CO2 */
        species[13] =
            +2.35677352e+00 * tc[0]
            +8.98459677e-03 * tc[1]
            -3.56178134e-06 * tc[2]
            +8.19730073e-10 * tc[3]
            -3.59248870e-14 * tc[4]
            +9.90105222e+00 ;
        /*species 14: HCO */
        species[14] =
            +4.22118584e+00 * tc[0]
            -3.24392532e-03 * tc[1]
            +6.88997230e-06 * tc[2]
            -4.43813643e-09 * tc[3]
            +1.08442216e-12 * tc[4]
            +3.39437243e+00 ;
        /*species 15: CH2O */
        species[15] =
            +4.79372315e+00 * tc[0]
            -9.90833369e-03 * tc[1]
            +1.86610004e-05 * tc[2]
            -1.26428420e-08 * tc[3]
            +3.29431630e-12 * tc[4]
            +6.02812900e-01 ;
        /*species 16: CH3O */
        species[16] =
            +2.10620400e+00 * tc[0]
            +7.21659500e-03 * tc[1]
            +2.66923600e-06 * tc[2]
            -2.45921200e-09 * tc[3]
            +5.18902500e-13 * tc[4]
            +1.31521770e+01 ;
        /*species 17: C2H2 */
        species[17] =
            +8.08681094e-01 * tc[0]
            +2.33615629e-02 * tc[1]
            -1.77585907e-05 * tc[2]
            +9.33841457e-09 * tc[3]
            -2.12518243e-12 * tc[4]
            +1.39397051e+01 ;
        /*species 18: C2H3 */
        species[18] =
            +3.21246645e+00 * tc[0]
            +1.51479162e-03 * tc[1]
            +1.29604706e-05 * tc[2]
            -1.19219282e-08 * tc[3]
            +3.67877182e-12 * tc[4]
            +8.51054025e+00 ;
        /*species 19: C2H4 */
        species[19] =
            +3.95920148e+00 * tc[0]
            -7.57052247e-03 * tc[1]
            +2.85495146e-05 * tc[2]
            -2.30529584e-08 * tc[3]
            +6.74710933e-12 * tc[4]
            +4.09733096e+00 ;
        /*species 20: C2H5 */
        species[20] =
            +4.30646568e+00 * tc[0]
            -4.18658892e-03 * tc[1]
            +2.48571403e-05 * tc[2]
            -1.99708869e-08 * tc[3]
            +5.76272510e-12 * tc[4]
            +4.70720924e+00 ;
        /*species 21: C2H6 */
        species[21] =
            +4.29142492e+00 * tc[0]
            -5.50154270e-03 * tc[1]
            +2.99719144e-05 * tc[2]
            -2.36155428e-08 * tc[3]
            +6.71714427e-12 * tc[4]
            +2.66682316e+00 ;
        /*species 22: N2 */
        species[22] =
            +3.29867700e+00 * tc[0]
            +1.40824040e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213500e-13 * tc[4]
            +3.95037200e+00 ;
        /*species 23: AR */
        species[23] =
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
        /*species 8: CH2 */
        species[8] =
            +2.87410113e+00 * tc[0]
            +3.65639292e-03 * tc[1]
            -7.04472985e-07 * tc[2]
            +8.67265163e-11 * tc[3]
            -4.69318918e-15 * tc[4]
            +6.17119324e+00 ;
        /*species 9: CH2(S) */
        species[9] =
            +2.29203842e+00 * tc[0]
            +4.65588637e-03 * tc[1]
            -1.00595973e-06 * tc[2]
            +1.39302000e-10 * tc[3]
            -8.49290912e-15 * tc[4]
            +8.62650169e+00 ;
        /*species 10: CH3 */
        species[10] =
            +2.28571772e+00 * tc[0]
            +7.23990037e-03 * tc[1]
            -1.49357174e-06 * tc[2]
            +1.98561548e-10 * tc[3]
            -1.16788599e-14 * tc[4]
            +8.48007179e+00 ;
        /*species 11: CH4 */
        species[11] =
            +7.48514950e-02 * tc[0]
            +1.33909467e-02 * tc[1]
            -2.86642905e-06 * tc[2]
            +4.07641783e-10 * tc[3]
            -2.54538075e-14 * tc[4]
            +1.84373180e+01 ;
        /*species 12: CO */
        species[12] =
            +2.71518561e+00 * tc[0]
            +2.06252743e-03 * tc[1]
            -4.99412886e-07 * tc[2]
            +7.66843360e-11 * tc[3]
            -5.09119290e-15 * tc[4]
            +7.81868772e+00 ;
        /*species 13: CO2 */
        species[13] =
            +3.85746029e+00 * tc[0]
            +4.41437026e-03 * tc[1]
            -1.10740702e-06 * tc[2]
            +1.74496729e-10 * tc[3]
            -1.18021041e-14 * tc[4]
            +2.27163806e+00 ;
        /*species 14: HCO */
        species[14] =
            +2.77217438e+00 * tc[0]
            +4.95695526e-03 * tc[1]
            -1.24222806e-06 * tc[2]
            +1.96387259e-10 * tc[3]
            -1.33377178e-14 * tc[4]
            +9.79834492e+00 ;
        /*species 15: CH2O */
        species[15] =
            +1.76069008e+00 * tc[0]
            +9.20000082e-03 * tc[1]
            -2.21129406e-06 * tc[2]
            +3.35470707e-10 * tc[3]
            -2.20963910e-14 * tc[4]
            +1.36563230e+01 ;
        /*species 16: CH3O */
        species[16] =
            +3.77079900e+00 * tc[0]
            +7.87149700e-03 * tc[1]
            -1.32819200e-06 * tc[2]
            +1.31481033e-10 * tc[3]
            -5.28154000e-15 * tc[4]
            +2.92957500e+00 ;
        /*species 17: C2H2 */
        species[17] =
            +4.14756964e+00 * tc[0]
            +5.96166664e-03 * tc[1]
            -1.18647426e-06 * tc[2]
            +1.55804057e-10 * tc[3]
            -9.03088033e-15 * tc[4]
            -1.23028121e+00 ;
        /*species 18: C2H3 */
        species[18] =
            +3.01672400e+00 * tc[0]
            +1.03302292e-02 * tc[1]
            -2.34041174e-06 * tc[2]
            +3.39210960e-10 * tc[3]
            -2.15651760e-14 * tc[4]
            +7.78732378e+00 ;
        /*species 19: C2H4 */
        species[19] =
            +2.03611116e+00 * tc[0]
            +1.46454151e-02 * tc[1]
            -3.35538958e-06 * tc[2]
            +4.90743077e-10 * tc[3]
            -3.14265152e-14 * tc[4]
            +1.03053693e+01 ;
        /*species 20: C2H5 */
        species[20] =
            +1.95465642e+00 * tc[0]
            +1.73972722e-02 * tc[1]
            -3.99103334e-06 * tc[2]
            +5.84058963e-10 * tc[3]
            -3.74103940e-14 * tc[4]
            +1.34624343e+01 ;
        /*species 21: C2H6 */
        species[21] =
            +1.07188150e+00 * tc[0]
            +2.16852677e-02 * tc[1]
            -5.01280335e-06 * tc[2]
            +7.38040003e-10 * tc[3]
            -4.75007225e-14 * tc[4]
            +1.51156107e+01 ;
        /*species 22: N2 */
        species[22] =
            +2.92664000e+00 * tc[0]
            +1.48797680e-03 * tc[1]
            -2.84238000e-07 * tc[2]
            +3.36567933e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
        /*species 23: AR */
        species[23] =
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
    wt[8] = 14.027090; /*CH2 */
    wt[9] = 14.027090; /*CH2(S) */
    wt[10] = 15.035060; /*CH3 */
    wt[11] = 16.043030; /*CH4 */
    wt[12] = 28.010550; /*CO */
    wt[13] = 44.009950; /*CO2 */
    wt[14] = 29.018520; /*HCO */
    wt[15] = 30.026490; /*CH2O */
    wt[16] = 31.034460; /*CH3O */
    wt[17] = 26.038240; /*C2H2 */
    wt[18] = 27.046210; /*C2H3 */
    wt[19] = 28.054180; /*C2H4 */
    wt[20] = 29.062150; /*C2H5 */
    wt[21] = 30.070120; /*C2H6 */
    wt[22] = 28.013400; /*N2 */
    wt[23] = 39.948000; /*AR */

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
    double tmax = 4000;/*min upper bound for thermo def */
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
  *LENIMC =           98;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENRMC EGTRANSETLENRMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENRMC egtransetlenrmc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENRMC egtransetlenrmc_
#endif
void egtransetLENRMC(int* LENRMC) {
  *LENRMC =        11784;}
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
  *KK =           24;}
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
  WT[           8] =   0.1402709031105042E+02;
  WT[           9] =   0.1402709031105042E+02;
  WT[          10] =   0.1503506028652191E+02;
  WT[          11] =   0.1604303026199341E+02;
  WT[          12] =   0.2801055049896240E+02;
  WT[          13] =   0.4400995063781738E+02;
  WT[          14] =   0.2901852047443390E+02;
  WT[          15] =   0.3002649044990540E+02;
  WT[          16] =   0.3103446042537689E+02;
  WT[          17] =   0.2603824067115784E+02;
  WT[          18] =   0.2704621064662933E+02;
  WT[          19] =   0.2805418062210083E+02;
  WT[          20] =   0.2906215059757233E+02;
  WT[          21] =   0.3007012057304382E+02;
  WT[          22] =   0.2801339912414551E+02;
  WT[          23] =   0.3994800186157227E+02;
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
  EPS[           8] =   0.1440000000000000E+03;
  EPS[           9] =   0.1440000000000000E+03;
  EPS[          10] =   0.1440000000000000E+03;
  EPS[          11] =   0.1414000000000000E+03;
  EPS[          12] =   0.9809999999999999E+02;
  EPS[          13] =   0.2440000000000000E+03;
  EPS[          14] =   0.4980000000000000E+03;
  EPS[          15] =   0.4980000000000000E+03;
  EPS[          16] =   0.4170000000000000E+03;
  EPS[          17] =   0.2090000000000000E+03;
  EPS[          18] =   0.2090000000000000E+03;
  EPS[          19] =   0.2808000000000000E+03;
  EPS[          20] =   0.2523000000000000E+03;
  EPS[          21] =   0.2523000000000000E+03;
  EPS[          22] =   0.9753000000000000E+02;
  EPS[          23] =   0.1365000000000000E+03;
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
  SIG[           8] =   0.3800000000000000E+01;
  SIG[           9] =   0.3800000000000000E+01;
  SIG[          10] =   0.3800000000000000E+01;
  SIG[          11] =   0.3746000000000000E+01;
  SIG[          12] =   0.3650000000000000E+01;
  SIG[          13] =   0.3763000000000000E+01;
  SIG[          14] =   0.3590000000000000E+01;
  SIG[          15] =   0.3590000000000000E+01;
  SIG[          16] =   0.3690000000000000E+01;
  SIG[          17] =   0.4100000000000000E+01;
  SIG[          18] =   0.4100000000000000E+01;
  SIG[          19] =   0.3971000000000000E+01;
  SIG[          20] =   0.4302000000000000E+01;
  SIG[          21] =   0.4302000000000000E+01;
  SIG[          22] =   0.3621000000000000E+01;
  SIG[          23] =   0.3330000000000000E+01;
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
  DIP[           9] =   0.0000000000000000E+00;
  DIP[          10] =   0.0000000000000000E+00;
  DIP[          11] =   0.0000000000000000E+00;
  DIP[          12] =   0.0000000000000000E+00;
  DIP[          13] =   0.0000000000000000E+00;
  DIP[          14] =   0.0000000000000000E+00;
  DIP[          15] =   0.0000000000000000E+00;
  DIP[          16] =   0.1700000000000000E+01;
  DIP[          17] =   0.0000000000000000E+00;
  DIP[          18] =   0.0000000000000000E+00;
  DIP[          19] =   0.0000000000000000E+00;
  DIP[          20] =   0.0000000000000000E+00;
  DIP[          21] =   0.0000000000000000E+00;
  DIP[          22] =   0.0000000000000000E+00;
  DIP[          23] =   0.0000000000000000E+00;
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
  POL[           8] =   0.0000000000000000E+00;
  POL[           9] =   0.0000000000000000E+00;
  POL[          10] =   0.0000000000000000E+00;
  POL[          11] =   0.2600000000000000E+01;
  POL[          12] =   0.1950000000000000E+01;
  POL[          13] =   0.2650000000000000E+01;
  POL[          14] =   0.0000000000000000E+00;
  POL[          15] =   0.0000000000000000E+00;
  POL[          16] =   0.0000000000000000E+00;
  POL[          17] =   0.0000000000000000E+00;
  POL[          18] =   0.0000000000000000E+00;
  POL[          19] =   0.0000000000000000E+00;
  POL[          20] =   0.0000000000000000E+00;
  POL[          21] =   0.0000000000000000E+00;
  POL[          22] =   0.1760000000000000E+01;
  POL[          23] =   0.0000000000000000E+00;
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
  ZROT[           8] =   0.0000000000000000E+00;
  ZROT[           9] =   0.0000000000000000E+00;
  ZROT[          10] =   0.0000000000000000E+00;
  ZROT[          11] =   0.1300000000000000E+02;
  ZROT[          12] =   0.1800000000000000E+01;
  ZROT[          13] =   0.2100000000000000E+01;
  ZROT[          14] =   0.0000000000000000E+00;
  ZROT[          15] =   0.2000000000000000E+01;
  ZROT[          16] =   0.2000000000000000E+01;
  ZROT[          17] =   0.2500000000000000E+01;
  ZROT[          18] =   0.1000000000000000E+01;
  ZROT[          19] =   0.1500000000000000E+01;
  ZROT[          20] =   0.1500000000000000E+01;
  ZROT[          21] =   0.1500000000000000E+01;
  ZROT[          22] =   0.4000000000000000E+01;
  ZROT[          23] =   0.0000000000000000E+00;
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
  NLIN[           9] =            1;
  NLIN[          10] =            1;
  NLIN[          11] =            2;
  NLIN[          12] =            1;
  NLIN[          13] =            1;
  NLIN[          14] =            2;
  NLIN[          15] =            2;
  NLIN[          16] =            2;
  NLIN[          17] =            1;
  NLIN[          18] =            2;
  NLIN[          19] =            2;
  NLIN[          20] =            2;
  NLIN[          21] =            2;
  NLIN[          22] =            1;
  NLIN[          23] =            0;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFLAM EGTRANSETCOFLAM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFLAM egtransetcoflam
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFLAM egtransetcoflam_
#endif
void egtransetCOFLAM(double* COFLAM) {
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
  COFLAM[          32] =   0.1291635112732534E+02;
  COFLAM[          33] =  -0.3737065816673101E+01;
  COFLAM[          34] =   0.7157937477351572E+00;
  COFLAM[          35] =  -0.3638374043162176E-01;
  COFLAM[          36] =   0.1893642052511076E+02;
  COFLAM[          37] =  -0.6509584574818159E+01;
  COFLAM[          38] =   0.1132853049549838E+01;
  COFLAM[          39] =  -0.5695791215151112E-01;
  COFLAM[          40] =   0.1399177056183472E+02;
  COFLAM[          41] =  -0.4641882989478738E+01;
  COFLAM[          42] =   0.9076435786312527E+00;
  COFLAM[          43] =  -0.4772397826828819E-01;
  COFLAM[          44] =   0.1330618431273573E+02;
  COFLAM[          45] =  -0.4960294456747341E+01;
  COFLAM[          46] =   0.1032808842591664E+01;
  COFLAM[          47] =  -0.5633567903356620E-01;
  COFLAM[          48] =   0.1187710072654253E+02;
  COFLAM[          49] =  -0.3154801252509626E+01;
  COFLAM[          50] =   0.6020483454882906E+00;
  COFLAM[          51] =  -0.3032714732778013E-01;
  COFLAM[          52] =  -0.1135070822782991E+02;
  COFLAM[          53] =   0.5875667873528579E+01;
  COFLAM[          54] =  -0.5677982250304914E+00;
  COFLAM[          55] =   0.2031670238952729E-01;
  COFLAM[          56] =   0.6296028277200143E+01;
  COFLAM[          57] =  -0.2225281214264024E+01;
  COFLAM[          58] =   0.6369149369376514E+00;
  COFLAM[          59] =  -0.3808451133860326E-01;
  COFLAM[          60] =   0.5384151635429418E+01;
  COFLAM[          61] =  -0.2389146833355944E+01;
  COFLAM[          62] =   0.7389861798813848E+00;
  COFLAM[          63] =  -0.4581386404459022E-01;
  COFLAM[          64] =  -0.6138056391436082E+01;
  COFLAM[          65] =   0.2471262972505274E+01;
  COFLAM[          66] =   0.6476625884336293E-01;
  COFLAM[          67] =  -0.1455104151494988E-01;
  COFLAM[          68] =  -0.7691048523709211E+01;
  COFLAM[          69] =   0.4564166690336450E+01;
  COFLAM[          70] =  -0.4040787948018396E+00;
  COFLAM[          71] =   0.1405248077558339E-01;
  COFLAM[          72] =  -0.9095387619073980E+01;
  COFLAM[          73] =   0.4544249531780739E+01;
  COFLAM[          74] =  -0.3175822855243875E+00;
  COFLAM[          75] =   0.6570713165664151E-02;
  COFLAM[          76] =  -0.1460870939008838E+02;
  COFLAM[          77] =   0.6359722173433303E+01;
  COFLAM[          78] =  -0.5034536870651064E+00;
  COFLAM[          79] =   0.1259513381293373E-01;
  COFLAM[          80] =  -0.8941827721234478E+01;
  COFLAM[          81] =   0.4021583693700529E+01;
  COFLAM[          82] =  -0.1835686948063458E+00;
  COFLAM[          83] =  -0.1963330023789806E-02;
  COFLAM[          84] =  -0.1098227802311090E+02;
  COFLAM[          85] =   0.4703047406066884E+01;
  COFLAM[          86] =  -0.2517963879729279E+00;
  COFLAM[          87] =   0.1532853967214294E-03;
  COFLAM[          88] =   0.1293004274651541E+02;
  COFLAM[          89] =  -0.3528374680486067E+01;
  COFLAM[          90] =   0.6455829015382131E+00;
  COFLAM[          91] =  -0.3194413600157287E-01;
  COFLAM[          92] =  -0.3166636127434923E+01;
  COFLAM[          93] =   0.3467381630049378E+01;
  COFLAM[          94] =  -0.3746257298439007E+00;
  COFLAM[          95] =   0.1658331946651994E-01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFETA EGTRANSETCOFETA
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFETA egtransetcofeta
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFETA egtransetcofeta_
#endif
void egtransetCOFETA(double* COFETA) {
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
  COFETA[          32] =  -0.2026303235679370E+02;
  COFETA[          33] =   0.3630408718312100E+01;
  COFETA[          34] =  -0.3952256976799819E+00;
  COFETA[          35] =   0.1745288966941290E-01;
  COFETA[          36] =  -0.2026303235679370E+02;
  COFETA[          37] =   0.3630408718312100E+01;
  COFETA[          38] =  -0.3952256976799819E+00;
  COFETA[          39] =   0.1745288966941290E-01;
  COFETA[          40] =  -0.2022833518474935E+02;
  COFETA[          41] =   0.3630408718312077E+01;
  COFETA[          42] =  -0.3952256976799788E+00;
  COFETA[          43] =   0.1745288966941275E-01;
  COFETA[          44] =  -0.2000457400541702E+02;
  COFETA[          45] =   0.3569542093003574E+01;
  COFETA[          46] =  -0.3874920392785122E+00;
  COFETA[          47] =   0.1712461411247914E-01;
  COFETA[          48] =  -0.1661561262594180E+02;
  COFETA[          49] =   0.2400975158113569E+01;
  COFETA[          50] =  -0.2357717790312281E+00;
  COFETA[          51] =   0.1054820948438182E-01;
  COFETA[          52] =  -0.2397057295682966E+02;
  COFETA[          53] =   0.5130426196036950E+01;
  COFETA[          54] =  -0.5724284704186094E+00;
  COFETA[          55] =   0.2440888721969576E-01;
  COFETA[          56] =  -0.1987003264957855E+02;
  COFETA[          57] =   0.2703811616927088E+01;
  COFETA[          58] =  -0.1672355189642018E+00;
  COFETA[          59] =   0.3212257118450056E-02;
  COFETA[          60] =  -0.1985295977766010E+02;
  COFETA[          61] =   0.2703811616927011E+01;
  COFETA[          62] =  -0.1672355189641908E+00;
  COFETA[          63] =   0.3212257118449540E-02;
  COFETA[          64] =  -0.1997384799268116E+02;
  COFETA[          65] =   0.2861750522819142E+01;
  COFETA[          66] =  -0.2024899143538781E+00;
  COFETA[          67] =   0.5362464436156140E-02;
  COFETA[          68] =  -0.2333653477658102E+02;
  COFETA[          69] =   0.4790351551626705E+01;
  COFETA[          70] =  -0.5364560275637680E+00;
  COFETA[          71] =   0.2318560947296609E-01;
  COFETA[          72] =  -0.2331754444859468E+02;
  COFETA[          73] =   0.4790351551626613E+01;
  COFETA[          74] =  -0.5364560275637551E+00;
  COFETA[          75] =   0.2318560947296548E-01;
  COFETA[          76] =  -0.2504838169296130E+02;
  COFETA[          77] =   0.5332796564390715E+01;
  COFETA[          78] =  -0.5890465746677508E+00;
  COFETA[          79] =   0.2473730756384840E-01;
  COFETA[          80] =  -0.2463014790579715E+02;
  COFETA[          81] =   0.5183142307718661E+01;
  COFETA[          82] =  -0.5771878965089483E+00;
  COFETA[          83] =   0.2453023678847775E-01;
  COFETA[          84] =  -0.2461310023284140E+02;
  COFETA[          85] =   0.5183142307718687E+01;
  COFETA[          86] =  -0.5771878965089522E+00;
  COFETA[          87] =   0.2453023678847794E-01;
  COFETA[          88] =  -0.1656563666406150E+02;
  COFETA[          89] =   0.2388167035581858E+01;
  COFETA[          90] =  -0.2341208182867099E+00;
  COFETA[          91] =   0.1047727172770704E-01;
  COFETA[          92] =  -0.1903691114465788E+02;
  COFETA[          93] =   0.3467381630049308E+01;
  COFETA[          94] =  -0.3746257298438905E+00;
  COFETA[          95] =   0.1658331946651944E-01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFD EGTRANSETCOFD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFD egtransetcofd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFD egtransetcofd_
#endif
void egtransetCOFD(double* COFD) {
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
  COFD[          32] =  -0.1355980676188247E+02;
  COFD[          33] =   0.3231182295054575E+01;
  COFD[          34] =  -0.2138406992439825E+00;
  COFD[          35] =   0.9659145053095021E-02;
  COFD[          36] =  -0.1355980676188247E+02;
  COFD[          37] =   0.3231182295054575E+01;
  COFD[          38] =  -0.2138406992439825E+00;
  COFD[          39] =   0.9659145053095021E-02;
  COFD[          40] =  -0.1358415095584097E+02;
  COFD[          41] =   0.3239845970896248E+01;
  COFD[          42] =  -0.2150500406377150E+00;
  COFD[          43] =   0.9715225271588263E-02;
  COFD[          44] =  -0.1353824314289343E+02;
  COFD[          45] =   0.3227352939332374E+01;
  COFD[          46] =  -0.2134245159235849E+00;
  COFD[          47] =   0.9644648995332943E-02;
  COFD[          48] =  -0.1268296956267843E+02;
  COFD[          49] =   0.2903428401599904E+01;
  COFD[          50] =  -0.1707667302851327E+00;
  COFD[          51] =   0.7771111802295727E-02;
  COFD[          52] =  -0.1564739122681956E+02;
  COFD[          53] =   0.4025825276486324E+01;
  COFD[          54] =  -0.3181195147607377E+00;
  COFD[          55] =   0.1422136468533644E-01;
  COFD[          56] =  -0.1746999275651299E+02;
  COFD[          57] =   0.4670556529920976E+01;
  COFD[          58] =  -0.3931788949833612E+00;
  COFD[          59] =   0.1709979614643272E-01;
  COFD[          60] =  -0.1747771970561251E+02;
  COFD[          61] =   0.4673159679696371E+01;
  COFD[          62] =  -0.3935109715079382E+00;
  COFD[          63] =   0.1711386921544412E-01;
  COFD[          64] =  -0.1797854335218850E+02;
  COFD[          65] =   0.4905020138101169E+01;
  COFD[          66] =  -0.4293084201751438E+00;
  COFD[          67] =   0.1891302999394448E-01;
  COFD[          68] =  -0.1528921665984086E+02;
  COFD[          69] =   0.3867752822044799E+01;
  COFD[          70] =  -0.2988043102116158E+00;
  COFD[          71] =   0.1343910541069478E-01;
  COFD[          72] =  -0.1530201153719761E+02;
  COFD[          73] =   0.3872623565106621E+01;
  COFD[          74] =  -0.2994802960421091E+00;
  COFD[          75] =   0.1347028558095244E-01;
  COFD[          76] =  -0.1646041046594534E+02;
  COFD[          77] =   0.4329097783536894E+01;
  COFD[          78] =  -0.3588538378541796E+00;
  COFD[          79] =   0.1604603061706265E-01;
  COFD[          80] =  -0.1615267929807825E+02;
  COFD[          81] =   0.4174121490551778E+01;
  COFD[          82] =  -0.3383475624805202E+00;
  COFD[          83] =   0.1513993577326728E-01;
  COFD[          84] =  -0.1616515483252968E+02;
  COFD[          85] =   0.4178910660883203E+01;
  COFD[          86] =  -0.3390062877759687E+00;
  COFD[          87] =   0.1517006334795817E-01;
  COFD[          88] =  -0.1266066516595179E+02;
  COFD[          89] =   0.2898076057146818E+01;
  COFD[          90] =  -0.1700453759318706E+00;
  COFD[          91] =   0.7738690296212197E-02;
  COFD[          92] =  -0.1342264266757377E+02;
  COFD[          93] =   0.3221777556990338E+01;
  COFD[          94] =  -0.2128673394230190E+00;
  COFD[          95] =   0.9627787551744238E-02;
  COFD[          96] =  -0.1201910564637032E+02;
  COFD[          97] =   0.3033139269467014E+01;
  COFD[          98] =  -0.1860150369684279E+00;
  COFD[          99] =   0.8361628520974778E-02;
  COFD[         100] =  -0.1519084593013382E+02;
  COFD[         101] =   0.4384306309333000E+01;
  COFD[         102] =  -0.3557550068895214E+00;
  COFD[         103] =   0.1548014912218549E-01;
  COFD[         104] =  -0.1531287529990816E+02;
  COFD[         105] =   0.4270888061293578E+01;
  COFD[         106] =  -0.3482328696611331E+00;
  COFD[         107] =   0.1545097630860032E-01;
  COFD[         108] =  -0.1761107854705398E+02;
  COFD[         109] =   0.5048254268555682E+01;
  COFD[         110] =  -0.4490719655111177E+00;
  COFD[         111] =   0.1981332179135272E-01;
  COFD[         112] =  -0.1533072430565730E+02;
  COFD[         113] =   0.4277567720643019E+01;
  COFD[         114] =  -0.3491407301790579E+00;
  COFD[         115] =   0.1549202933974097E-01;
  COFD[         116] =  -0.1651315089501800E+02;
  COFD[         117] =   0.4220701654162998E+01;
  COFD[         118] =  -0.2825947916575933E+00;
  COFD[         119] =   0.1018065117940972E-01;
  COFD[         120] =  -0.1762111333528786E+02;
  COFD[         121] =   0.5052180357519856E+01;
  COFD[         122] =  -0.4496001366665920E+00;
  COFD[         123] =   0.1983696581079775E-01;
  COFD[         124] =  -0.1763070827156596E+02;
  COFD[         125] =   0.5055937247475734E+01;
  COFD[         126] =  -0.4501055595390027E+00;
  COFD[         127] =   0.1985959199754691E-01;
  COFD[         128] =  -0.1798724510881690E+02;
  COFD[         129] =   0.5079566464178115E+01;
  COFD[         130] =  -0.4444140659462563E+00;
  COFD[         131] =   0.1923462920552580E-01;
  COFD[         132] =  -0.1798724510881690E+02;
  COFD[         133] =   0.5079566464178115E+01;
  COFD[         134] =  -0.4444140659462563E+00;
  COFD[         135] =   0.1923462920552580E-01;
  COFD[         136] =  -0.1801636950310948E+02;
  COFD[         137] =   0.5090088857460223E+01;
  COFD[         138] =  -0.4457552498826564E+00;
  COFD[         139] =   0.1929141580461785E-01;
  COFD[         140] =  -0.1790236652645725E+02;
  COFD[         141] =   0.5054241512059871E+01;
  COFD[         142] =  -0.4413163862361035E+00;
  COFD[         143] =   0.1910841720545356E-01;
  COFD[         144] =  -0.1719132602269712E+02;
  COFD[         145] =   0.4859367175853364E+01;
  COFD[         146] =  -0.4240580935572132E+00;
  COFD[         147] =   0.1870696403708276E-01;
  COFD[         148] =  -0.1725034548774479E+02;
  COFD[         149] =   0.4566921339972978E+01;
  COFD[         150] =  -0.3509192271095338E+00;
  COFD[         151] =   0.1402685086906688E-01;
  COFD[         152] =  -0.1540645026491837E+02;
  COFD[         153] =   0.3572427121876953E+01;
  COFD[         154] =  -0.1878763535992870E+00;
  COFD[         155] =   0.5591984032298197E-02;
  COFD[         156] =  -0.1539084008249366E+02;
  COFD[         157] =   0.3565006575777335E+01;
  COFD[         158] =  -0.1867580729432450E+00;
  COFD[         159] =   0.5536884066371989E-02;
  COFD[         160] =  -0.1561125740910447E+02;
  COFD[         161] =   0.3698144089263720E+01;
  COFD[         162] =  -0.2103080207621357E+00;
  COFD[         163] =   0.6792204226221277E-02;
  COFD[         164] =  -0.1774507254557799E+02;
  COFD[         165] =   0.4780159131489816E+01;
  COFD[         166] =  -0.3856587269477316E+00;
  COFD[         167] =   0.1582942762638978E-01;
  COFD[         168] =  -0.1774834356244655E+02;
  COFD[         169] =   0.4780626617651506E+01;
  COFD[         170] =  -0.3856327873318030E+00;
  COFD[         171] =   0.1582440706500046E-01;
  COFD[         172] =  -0.1694497223361135E+02;
  COFD[         173] =   0.4353290721870846E+01;
  COFD[         174] =  -0.3147441391485086E+00;
  COFD[         175] =   0.1210201943729390E-01;
  COFD[         176] =  -0.1738386031480564E+02;
  COFD[         177] =   0.4515894896789344E+01;
  COFD[         178] =  -0.3391856065220307E+00;
  COFD[         179] =   0.1329276720561771E-01;
  COFD[         180] =  -0.1738223879329547E+02;
  COFD[         181] =   0.4514299549461490E+01;
  COFD[         182] =  -0.3388651702152571E+00;
  COFD[         183] =   0.1327383213284630E-01;
  COFD[         184] =  -0.1711567017943538E+02;
  COFD[         185] =   0.4833264586465570E+01;
  COFD[         186] =  -0.4205908440542898E+00;
  COFD[         187] =   0.1855345683322859E-01;
  COFD[         188] =  -0.1752038279713604E+02;
  COFD[         189] =   0.4959639086677385E+01;
  COFD[         190] =  -0.4293974875408633E+00;
  COFD[         191] =   0.1860806087249651E-01;
  COFD[         192] =  -0.1093372460559999E+02;
  COFD[         193] =   0.2306133612968357E+01;
  COFD[         194] =  -0.8736427899565385E-01;
  COFD[         195] =   0.3897526119633928E-02;
  COFD[         196] =  -0.1531287529990816E+02;
  COFD[         197] =   0.4270888061293578E+01;
  COFD[         198] =  -0.3482328696611331E+00;
  COFD[         199] =   0.1545097630860032E-01;
  COFD[         200] =  -0.1356847565002614E+02;
  COFD[         201] =   0.3062175626208075E+01;
  COFD[         202] =  -0.1889691489381823E+00;
  COFD[         203] =   0.8454139853458284E-02;
  COFD[         204] =  -0.1511165027002312E+02;
  COFD[         205] =   0.3515441252015493E+01;
  COFD[         206] =  -0.2491198761221894E+00;
  COFD[         207] =   0.1111625360314578E-01;
  COFD[         208] =  -0.1358466752355666E+02;
  COFD[         209] =   0.3062671492876047E+01;
  COFD[         210] =  -0.1890385291581307E+00;
  COFD[         211] =   0.8457362732499073E-02;
  COFD[         212] =  -0.1906508600184398E+02;
  COFD[         213] =   0.4989055420187643E+01;
  COFD[         214] =  -0.4192832417141421E+00;
  COFD[         215] =   0.1761870656062123E-01;
  COFD[         216] =  -0.1513524421958843E+02;
  COFD[         217] =   0.3523297948981071E+01;
  COFD[         218] =  -0.2502118585975575E+00;
  COFD[         219] =   0.1116667641067510E-01;
  COFD[         220] =  -0.1515824683455671E+02;
  COFD[         221] =   0.3531011422568503E+01;
  COFD[         222] =  -0.2512839088019835E+00;
  COFD[         223] =   0.1121617798528330E-01;
  COFD[         224] =  -0.1558834457822415E+02;
  COFD[         225] =   0.3699567664433178E+01;
  COFD[         226] =  -0.2706504848826866E+00;
  COFD[         227] =   0.1194595865045622E-01;
  COFD[         228] =  -0.1558834457822415E+02;
  COFD[         229] =   0.3699567664433178E+01;
  COFD[         230] =  -0.2706504848826866E+00;
  COFD[         231] =   0.1194595865045622E-01;
  COFD[         232] =  -0.1561687122378518E+02;
  COFD[         233] =   0.3703886367486989E+01;
  COFD[         234] =  -0.2712410029515650E+00;
  COFD[         235] =   0.1197280722366742E-01;
  COFD[         236] =  -0.1556217965348612E+02;
  COFD[         237] =   0.3683628994574988E+01;
  COFD[         238] =  -0.2686556693855442E+00;
  COFD[         239] =   0.1186253607235680E-01;
  COFD[         240] =  -0.1478198899238034E+02;
  COFD[         241] =   0.3371092581023257E+01;
  COFD[         242] =  -0.2298567856863946E+00;
  COFD[         243] =   0.1025846107153025E-01;
  COFD[         244] =  -0.1793673861488134E+02;
  COFD[         245] =   0.4489875914439151E+01;
  COFD[         246] =  -0.3697556337020395E+00;
  COFD[         247] =   0.1609253879718407E-01;
  COFD[         248] =  -0.1915075646861156E+02;
  COFD[         249] =   0.4868838487330766E+01;
  COFD[         250] =  -0.4046210142562804E+00;
  COFD[         251] =   0.1701122981099066E-01;
  COFD[         252] =  -0.1914986770399175E+02;
  COFD[         253] =   0.4864899621234062E+01;
  COFD[         254] =  -0.4039220297382060E+00;
  COFD[         255] =   0.1697258735720344E-01;
  COFD[         256] =  -0.1878355293596396E+02;
  COFD[         257] =   0.4745970603113269E+01;
  COFD[         258] =  -0.3924372163577848E+00;
  COFD[         259] =   0.1663659729501466E-01;
  COFD[         260] =  -0.1731509452354737E+02;
  COFD[         261] =   0.4259553692469185E+01;
  COFD[         262] =  -0.3417709853226915E+00;
  COFD[         263] =   0.1496136243888912E-01;
  COFD[         264] =  -0.1734665788392033E+02;
  COFD[         265] =   0.4269427891180465E+01;
  COFD[         266] =  -0.3430766910038187E+00;
  COFD[         267] =   0.1501879627721415E-01;
  COFD[         268] =  -0.1800091724620567E+02;
  COFD[         269] =   0.4493072014688184E+01;
  COFD[         270] =  -0.3678451189139861E+00;
  COFD[         271] =   0.1591448159607913E-01;
  COFD[         272] =  -0.1802433609839966E+02;
  COFD[         273] =   0.4483646303442399E+01;
  COFD[         274] =  -0.3688151334706006E+00;
  COFD[         275] =   0.1604826293273014E-01;
  COFD[         276] =  -0.1805221831789452E+02;
  COFD[         277] =   0.4492279525932131E+01;
  COFD[         278] =  -0.3699242263001061E+00;
  COFD[         279] =   0.1609559465701492E-01;
  COFD[         280] =  -0.1474819919893963E+02;
  COFD[         281] =   0.3361311502126538E+01;
  COFD[         282] =  -0.2285465363406641E+00;
  COFD[         283] =   0.1019990809984005E-01;
  COFD[         284] =  -0.1581858256688188E+02;
  COFD[         285] =   0.3775899039344615E+01;
  COFD[         286] =  -0.2814687920927101E+00;
  COFD[         287] =   0.1245218573218981E-01;
  COFD[         288] =  -0.1270707915711888E+02;
  COFD[         289] =   0.2929619774221139E+01;
  COFD[         290] =  -0.1739670601708100E+00;
  COFD[         291] =   0.7901198812704345E-02;
  COFD[         292] =  -0.1761107854705398E+02;
  COFD[         293] =   0.5048254268555682E+01;
  COFD[         294] =  -0.4490719655111177E+00;
  COFD[         295] =   0.1981332179135272E-01;
  COFD[         296] =  -0.1511165027002312E+02;
  COFD[         297] =   0.3515441252015493E+01;
  COFD[         298] =  -0.2491198761221894E+00;
  COFD[         299] =   0.1111625360314578E-01;
  COFD[         300] =  -0.1605076758771229E+02;
  COFD[         301] =   0.3688226344868140E+01;
  COFD[         302] =  -0.2690947493495641E+00;
  COFD[         303] =   0.1187504521154671E-01;
  COFD[         304] =  -0.1509620843583920E+02;
  COFD[         305] =   0.3500304707531062E+01;
  COFD[         306] =  -0.2470160084537919E+00;
  COFD[         307] =   0.1101910392929726E-01;
  COFD[         308] =  -0.2029800570223752E+02;
  COFD[         309] =   0.5167991782317953E+01;
  COFD[         310] =  -0.4260628585538284E+00;
  COFD[         311] =   0.1725358239585169E-01;
  COFD[         312] =  -0.1605887758547131E+02;
  COFD[         313] =   0.3688399632843664E+01;
  COFD[         314] =  -0.2691183502617316E+00;
  COFD[         315] =   0.1187611401222301E-01;
  COFD[         316] =  -0.1606742328871227E+02;
  COFD[         317] =   0.3688898239072069E+01;
  COFD[         318] =  -0.2691862577212903E+00;
  COFD[         319] =   0.1187918929470605E-01;
  COFD[         320] =  -0.1707915522153321E+02;
  COFD[         321] =   0.4132970234506635E+01;
  COFD[         322] =  -0.3260547355537204E+00;
  COFD[         323] =   0.1431106723651207E-01;
  COFD[         324] =  -0.1707915522153321E+02;
  COFD[         325] =   0.4132970234506635E+01;
  COFD[         326] =  -0.3260547355537204E+00;
  COFD[         327] =   0.1431106723651207E-01;
  COFD[         328] =  -0.1707094794267795E+02;
  COFD[         329] =   0.4119810624341889E+01;
  COFD[         330] =  -0.3242844616145015E+00;
  COFD[         331] =   0.1423184749658152E-01;
  COFD[         332] =  -0.1696131964138645E+02;
  COFD[         333] =   0.4074167708212651E+01;
  COFD[         334] =  -0.3182576642303228E+00;
  COFD[         335] =   0.1396612551510442E-01;
  COFD[         336] =  -0.1582244348422915E+02;
  COFD[         337] =   0.3599211636009877E+01;
  COFD[         338] =  -0.2581898369286856E+00;
  COFD[         339] =   0.1143096575439008E-01;
  COFD[         340] =  -0.1868514242925209E+02;
  COFD[         341] =   0.4565512377700616E+01;
  COFD[         342] =  -0.3754538805580755E+00;
  COFD[         343] =   0.1617305773699574E-01;
  COFD[         344] =  -0.2031657960333580E+02;
  COFD[         345] =   0.5113614640951361E+01;
  COFD[         346] =  -0.4319293017089550E+00;
  COFD[         347] =   0.1802805482147992E-01;
  COFD[         348] =  -0.2032689236579443E+02;
  COFD[         349] =   0.5114221260960875E+01;
  COFD[         350] =  -0.4320136109530544E+00;
  COFD[         351] =   0.1803190494311217E-01;
  COFD[         352] =  -0.2016013281728147E+02;
  COFD[         353] =   0.5063566265537697E+01;
  COFD[         354] =  -0.4281515792120217E+00;
  COFD[         355] =   0.1797180686456706E-01;
  COFD[         356] =  -0.1817374730686181E+02;
  COFD[         357] =   0.4404209750065649E+01;
  COFD[         358] =  -0.3569181269134391E+00;
  COFD[         359] =   0.1547062305483025E-01;
  COFD[         360] =  -0.1818389332535593E+02;
  COFD[         361] =   0.4404158898744668E+01;
  COFD[         362] =  -0.3569180632138942E+00;
  COFD[         363] =   0.1547091340944365E-01;
  COFD[         364] =  -0.1901929468754062E+02;
  COFD[         365] =   0.4693839754776671E+01;
  COFD[         366] =  -0.3901000075815114E+00;
  COFD[         367] =   0.1672631492507185E-01;
  COFD[         368] =  -0.1881741117557303E+02;
  COFD[         369] =   0.4595519619722809E+01;
  COFD[         370] =  -0.3791053358880415E+00;
  COFD[         371] =   0.1632232422675412E-01;
  COFD[         372] =  -0.1882824574536174E+02;
  COFD[         373] =   0.4596254296231404E+01;
  COFD[         374] =  -0.3791928797888692E+00;
  COFD[         375] =   0.1632573982847240E-01;
  COFD[         376] =  -0.1578873902175925E+02;
  COFD[         377] =   0.3589258101555060E+01;
  COFD[         378] =  -0.2568718883472193E+00;
  COFD[         379] =   0.1137271971440306E-01;
  COFD[         380] =  -0.1677455230511400E+02;
  COFD[         381] =   0.3946214344521416E+01;
  COFD[         382] =  -0.3012613112287137E+00;
  COFD[         383] =   0.1321325422150532E-01;
  COFD[         384] =  -0.1093673058791694E+02;
  COFD[         385] =   0.2306118133758978E+01;
  COFD[         386] =  -0.8736238957367691E-01;
  COFD[         387] =   0.3897448280154413E-02;
  COFD[         388] =  -0.1533072430565730E+02;
  COFD[         389] =   0.4277567720643019E+01;
  COFD[         390] =  -0.3491407301790579E+00;
  COFD[         391] =   0.1549202933974097E-01;
  COFD[         392] =  -0.1358466752355666E+02;
  COFD[         393] =   0.3062671492876047E+01;
  COFD[         394] =  -0.1890385291581307E+00;
  COFD[         395] =   0.8457362732499073E-02;
  COFD[         396] =  -0.1509620843583920E+02;
  COFD[         397] =   0.3500304707531062E+01;
  COFD[         398] =  -0.2470160084537919E+00;
  COFD[         399] =   0.1101910392929726E-01;
  COFD[         400] =  -0.1359902342804023E+02;
  COFD[         401] =   0.3062175626208111E+01;
  COFD[         402] =  -0.1889691489381874E+00;
  COFD[         403] =   0.8454139853458534E-02;
  COFD[         404] =  -0.1908171031786331E+02;
  COFD[         405] =   0.4989633286602976E+01;
  COFD[         406] =  -0.4194071192530971E+00;
  COFD[         407] =   0.1762632910887997E-01;
  COFD[         408] =  -0.1511946285859442E+02;
  COFD[         409] =   0.3507927470352835E+01;
  COFD[         410] =  -0.2480755290441053E+00;
  COFD[         411] =   0.1106802953504252E-01;
  COFD[         412] =  -0.1514219804803715E+02;
  COFD[         413] =   0.3515441252015498E+01;
  COFD[         414] =  -0.2491198761221899E+00;
  COFD[         415] =   0.1111625360314580E-01;
  COFD[         416] =  -0.1559686009787643E+02;
  COFD[         417] =   0.3697254029246602E+01;
  COFD[         418] =  -0.2703332731709061E+00;
  COFD[         419] =   0.1193149792533076E-01;
  COFD[         420] =  -0.1559686009787643E+02;
  COFD[         421] =   0.3697254029246602E+01;
  COFD[         422] =  -0.2703332731709061E+00;
  COFD[         423] =   0.1193149792533076E-01;
  COFD[         424] =  -0.1562210625534589E+02;
  COFD[         425] =   0.3699988818188416E+01;
  COFD[         426] =  -0.2707081263160739E+00;
  COFD[         427] =   0.1194858184298051E-01;
  COFD[         428] =  -0.1556532980903788E+02;
  COFD[         429] =   0.3678661452138608E+01;
  COFD[         430] =  -0.2679777825820168E+00;
  COFD[         431] =   0.1183177506717497E-01;
  COFD[         432] =  -0.1476985219999908E+02;
  COFD[         433] =   0.3357751196108267E+01;
  COFD[         434] =  -0.2280046831004087E+00;
  COFD[         435] =   0.1017303667696902E-01;
  COFD[         436] =  -0.1792454027766008E+02;
  COFD[         437] =   0.4476326984714260E+01;
  COFD[         438] =  -0.3680191278135618E+00;
  COFD[         439] =   0.1601860672433683E-01;
  COFD[         440] =  -0.1918204714711579E+02;
  COFD[         441] =   0.4875615818554563E+01;
  COFD[         442] =  -0.4058262469897604E+00;
  COFD[         443] =   0.1707794755285167E-01;
  COFD[         444] =  -0.1918177070837826E+02;
  COFD[         445] =   0.4871873199573439E+01;
  COFD[         446] =  -0.4051602536432213E+00;
  COFD[         447] =   0.1704106539271711E-01;
  COFD[         448] =  -0.1880409131690931E+02;
  COFD[         449] =   0.4747896766189557E+01;
  COFD[         450] =  -0.3929513558466455E+00;
  COFD[         451] =   0.1667071927589917E-01;
  COFD[         452] =  -0.1729594517981932E+02;
  COFD[         453] =   0.4244246410490659E+01;
  COFD[         454] =  -0.3397463840897515E+00;
  COFD[         455] =   0.1487228746520592E-01;
  COFD[         456] =  -0.1732667376192460E+02;
  COFD[         457] =   0.4253674558730024E+01;
  COFD[         458] =  -0.3409934577271315E+00;
  COFD[         459] =   0.1492715703864712E-01;
  COFD[         460] =  -0.1799417924090561E+02;
  COFD[         461] =   0.4483237269191199E+01;
  COFD[         462] =  -0.3666457816700467E+00;
  COFD[         463] =   0.1586624138953086E-01;
  COFD[         464] =  -0.1800542957417654E+02;
  COFD[         465] =   0.4468504271448043E+01;
  COFD[         466] =  -0.3668692513471102E+00;
  COFD[         467] =   0.1596519376748332E-01;
  COFD[         468] =  -0.1803293840758895E+02;
  COFD[         469] =   0.4476898762718537E+01;
  COFD[         470] =  -0.3679481179513786E+00;
  COFD[         471] =   0.1601125468656841E-01;
  COFD[         472] =  -0.1473661449825940E+02;
  COFD[         473] =   0.3348204558833826E+01;
  COFD[         474] =  -0.2267271723233180E+00;
  COFD[         475] =   0.1011600240359858E-01;
  COFD[         476] =  -0.1580533640016843E+02;
  COFD[         477] =   0.3761415063014861E+01;
  COFD[         478] =  -0.2795015971789234E+00;
  COFD[         479] =   0.1236332420639666E-01;
  COFD[         480] =  -0.1818497883283899E+02;
  COFD[         481] =   0.5043624597310085E+01;
  COFD[         482] =  -0.4377066809210735E+00;
  COFD[         483] =   0.1887643866477002E-01;
  COFD[         484] =  -0.1651315089501800E+02;
  COFD[         485] =   0.4220701654162998E+01;
  COFD[         486] =  -0.2825947916575933E+00;
  COFD[         487] =   0.1018065117940972E-01;
  COFD[         488] =  -0.1906508600184398E+02;
  COFD[         489] =   0.4989055420187643E+01;
  COFD[         490] =  -0.4192832417141421E+00;
  COFD[         491] =   0.1761870656062123E-01;
  COFD[         492] =  -0.2029800570223752E+02;
  COFD[         493] =   0.5167991782317953E+01;
  COFD[         494] =  -0.4260628585538284E+00;
  COFD[         495] =   0.1725358239585169E-01;
  COFD[         496] =  -0.1908171031786331E+02;
  COFD[         497] =   0.4989633286602976E+01;
  COFD[         498] =  -0.4194071192530971E+00;
  COFD[         499] =   0.1762632910887997E-01;
  COFD[         500] =  -0.1180242969654581E+02;
  COFD[         501] =   0.8887910786873973E+00;
  COFD[         502] =   0.2461702200576275E+00;
  COFD[         503] =  -0.1607058638443385E-01;
  COFD[         504] =  -0.1969226895876598E+02;
  COFD[         505] =   0.4975023535625081E+01;
  COFD[         506] =  -0.4062730945223436E+00;
  COFD[         507] =   0.1660121367605663E-01;
  COFD[         508] =  -0.1967561834414118E+02;
  COFD[         509] =   0.4964896032502449E+01;
  COFD[         510] =  -0.4047371090295282E+00;
  COFD[         511] =   0.1652517237194395E-01;
  COFD[         512] =  -0.2040913503195636E+02;
  COFD[         513] =   0.5224092639655129E+01;
  COFD[         514] =  -0.4327353018434304E+00;
  COFD[         515] =   0.1753064056853202E-01;
  COFD[         516] =  -0.2040913503195636E+02;
  COFD[         517] =   0.5224092639655129E+01;
  COFD[         518] =  -0.4327353018434304E+00;
  COFD[         519] =   0.1753064056853202E-01;
  COFD[         520] =  -0.2042453573019913E+02;
  COFD[         521] =   0.5222435968901098E+01;
  COFD[         522] =  -0.4324987976330752E+00;
  COFD[         523] =   0.1751958312757915E-01;
  COFD[         524] =  -0.2052802769939090E+02;
  COFD[         525] =   0.5190809253464677E+01;
  COFD[         526] =  -0.4197156907449909E+00;
  COFD[         527] =   0.1661902490158712E-01;
  COFD[         528] =  -0.2031234004461199E+02;
  COFD[         529] =   0.5186312454417227E+01;
  COFD[         530] =  -0.4313742876779729E+00;
  COFD[         531] =   0.1759974844922337E-01;
  COFD[         532] =  -0.1924815690645235E+02;
  COFD[         533] =   0.4356731855020855E+01;
  COFD[         534] =  -0.2835834611821879E+00;
  COFD[         535] =   0.9647123135613308E-02;
  COFD[         536] =  -0.1757312087019186E+02;
  COFD[         537] =   0.3511204681563826E+01;
  COFD[         538] =  -0.1497966127823777E+00;
  COFD[         539] =   0.2966580974304371E-02;
  COFD[         540] =  -0.1757470287826831E+02;
  COFD[         541] =   0.3509436887842453E+01;
  COFD[         542] =  -0.1496102927025162E+00;
  COFD[         543] =   0.2961974572337444E-02;
  COFD[         544] =  -0.1619587514529141E+02;
  COFD[         545] =   0.2843235032975332E+01;
  COFD[         546] =  -0.5104314541884596E-01;
  COFD[         547] =  -0.1707937587486378E-02;
  COFD[         548] =  -0.2024506268677845E+02;
  COFD[         549] =   0.4905829896937226E+01;
  COFD[         550] =  -0.3737705712203094E+00;
  COFD[         551] =   0.1428535217779750E-01;
  COFD[         552] =  -0.2021382909884265E+02;
  COFD[         553] =   0.4888560960839369E+01;
  COFD[         554] =  -0.3712523102934615E+00;
  COFD[         555] =   0.1416513273755347E-01;
  COFD[         556] =  -0.1988688758541057E+02;
  COFD[         557] =   0.4655348588622681E+01;
  COFD[         558] =  -0.3287200044949488E+00;
  COFD[         559] =   0.1186464366263145E-01;
  COFD[         560] =  -0.2010033725040245E+02;
  COFD[         561] =   0.4740608584978043E+01;
  COFD[         562] =  -0.3442903017689571E+00;
  COFD[         563] =   0.1271241178230842E-01;
  COFD[         564] =  -0.2007381993778366E+02;
  COFD[         565] =   0.4726146050811121E+01;
  COFD[         566] =  -0.3422109775805109E+00;
  COFD[         567] =   0.1261488904944177E-01;
  COFD[         568] =  -0.2024331300258610E+02;
  COFD[         569] =   0.5167058816761965E+01;
  COFD[         570] =  -0.4292618020154789E+00;
  COFD[         571] =   0.1752039422653532E-01;
  COFD[         572] =  -0.1991427613539219E+02;
  COFD[         573] =   0.4985074718206674E+01;
  COFD[         574] =  -0.3988046631693363E+00;
  COFD[         575] =   0.1592993698509921E-01;
  COFD[         576] =  -0.1271124360389947E+02;
  COFD[         577] =   0.2931047949505553E+01;
  COFD[         578] =  -0.1741679825526774E+00;
  COFD[         579] =   0.7910585398301936E-02;
  COFD[         580] =  -0.1762111333528786E+02;
  COFD[         581] =   0.5052180357519856E+01;
  COFD[         582] =  -0.4496001366665920E+00;
  COFD[         583] =   0.1983696581079775E-01;
  COFD[         584] =  -0.1513524421958843E+02;
  COFD[         585] =   0.3523297948981071E+01;
  COFD[         586] =  -0.2502118585975575E+00;
  COFD[         587] =   0.1116667641067510E-01;
  COFD[         588] =  -0.1605887758547131E+02;
  COFD[         589] =   0.3688399632843664E+01;
  COFD[         590] =  -0.2691183502617316E+00;
  COFD[         591] =   0.1187611401222301E-01;
  COFD[         592] =  -0.1511946285859442E+02;
  COFD[         593] =   0.3507927470352835E+01;
  COFD[         594] =  -0.2480755290441053E+00;
  COFD[         595] =   0.1106802953504252E-01;
  COFD[         596] =  -0.1969226895876598E+02;
  COFD[         597] =   0.4975023535625081E+01;
  COFD[         598] =  -0.4062730945223436E+00;
  COFD[         599] =   0.1660121367605663E-01;
  COFD[         600] =  -0.1606627473213039E+02;
  COFD[         601] =   0.3688226344868143E+01;
  COFD[         602] =  -0.2690947493495643E+00;
  COFD[         603] =   0.1187504521154673E-01;
  COFD[         604] =  -0.1607413035852147E+02;
  COFD[         605] =   0.3688389366224249E+01;
  COFD[         606] =  -0.2691169520025368E+00;
  COFD[         607] =   0.1187605069009256E-01;
  COFD[         608] =  -0.1709853282395253E+02;
  COFD[         609] =   0.4139002691011624E+01;
  COFD[         610] =  -0.3268662406958736E+00;
  COFD[         611] =   0.1434738242897968E-01;
  COFD[         612] =  -0.1709853282395253E+02;
  COFD[         613] =   0.4139002691011624E+01;
  COFD[         614] =  -0.3268662406958736E+00;
  COFD[         615] =   0.1434738242897968E-01;
  COFD[         616] =  -0.1709003104221693E+02;
  COFD[         617] =   0.4125627998285617E+01;
  COFD[         618] =  -0.3250670331906144E+00;
  COFD[         619] =   0.1426686739847163E-01;
  COFD[         620] =  -0.1698006191382480E+02;
  COFD[         621] =   0.4079755891690070E+01;
  COFD[         622] =  -0.3190093940831370E+00;
  COFD[         623] =   0.1399976445173436E-01;
  COFD[         624] =  -0.1583284545979817E+02;
  COFD[         625] =   0.3600569993944359E+01;
  COFD[         626] =  -0.2583773478248265E+00;
  COFD[         627] =   0.1143956672273590E-01;
  COFD[         628] =  -0.1868972720602865E+02;
  COFD[         629] =   0.4564003101584411E+01;
  COFD[         630] =  -0.3752914154261894E+00;
  COFD[         631] =   0.1616756912994680E-01;
  COFD[         632] =  -0.2032173485009814E+02;
  COFD[         633] =   0.5112656067238100E+01;
  COFD[         634] =  -0.4317894593025599E+00;
  COFD[         635] =   0.1802136540679988E-01;
  COFD[         636] =  -0.2033308000404864E+02;
  COFD[         637] =   0.5113689932626547E+01;
  COFD[         638] =  -0.4319400850127593E+00;
  COFD[         639] =   0.1802856188319298E-01;
  COFD[         640] =  -0.1994316533928780E+02;
  COFD[         641] =   0.4989293390342089E+01;
  COFD[         642] =  -0.4200959859596722E+00;
  COFD[         643] =   0.1768643162186107E-01;
  COFD[         644] =  -0.1818135452916891E+02;
  COFD[         645] =   0.4404445109750770E+01;
  COFD[         646] =  -0.3569414786935936E+00;
  COFD[         647] =   0.1547130685514495E-01;
  COFD[         648] =  -0.1819107613580373E+02;
  COFD[         649] =   0.4404180695515474E+01;
  COFD[         650] =  -0.3569157439691602E+00;
  COFD[         651] =   0.1547058259292524E-01;
  COFD[         652] =  -0.1902570873205879E+02;
  COFD[         653] =   0.4693487292163261E+01;
  COFD[         654] =  -0.3900471625097403E+00;
  COFD[         655] =   0.1672372332121293E-01;
  COFD[         656] =  -0.1882326100183751E+02;
  COFD[         657] =   0.4594957876626870E+01;
  COFD[         658] =  -0.3790360690538448E+00;
  COFD[         659] =   0.1631950746180167E-01;
  COFD[         660] =  -0.1883388798798641E+02;
  COFD[         661] =   0.4595580957493091E+01;
  COFD[         662] =  -0.3791127554800276E+00;
  COFD[         663] =   0.1632261913288707E-01;
  COFD[         664] =  -0.1579929195048582E+02;
  COFD[         665] =   0.3590679879226839E+01;
  COFD[         666] =  -0.2570681340587627E+00;
  COFD[         667] =   0.1138172060899465E-01;
  COFD[         668] =  -0.1677798202265061E+02;
  COFD[         669] =   0.3944105793182093E+01;
  COFD[         670] =  -0.3009768138615855E+00;
  COFD[         671] =   0.1320048849753941E-01;
  COFD[         672] =  -0.1271518733681371E+02;
  COFD[         673] =   0.2932402241140174E+01;
  COFD[         674] =  -0.1743585082612658E+00;
  COFD[         675] =   0.7919486191765343E-02;
  COFD[         676] =  -0.1763070827156596E+02;
  COFD[         677] =   0.5055937247475734E+01;
  COFD[         678] =  -0.4501055595390027E+00;
  COFD[         679] =   0.1985959199754691E-01;
  COFD[         680] =  -0.1515824683455671E+02;
  COFD[         681] =   0.3531011422568503E+01;
  COFD[         682] =  -0.2512839088019835E+00;
  COFD[         683] =   0.1121617798528330E-01;
  COFD[         684] =  -0.1606742328871227E+02;
  COFD[         685] =   0.3688898239072069E+01;
  COFD[         686] =  -0.2691862577212903E+00;
  COFD[         687] =   0.1187918929470605E-01;
  COFD[         688] =  -0.1514219804803715E+02;
  COFD[         689] =   0.3515441252015498E+01;
  COFD[         690] =  -0.2491198761221899E+00;
  COFD[         691] =   0.1111625360314580E-01;
  COFD[         692] =  -0.1967561834414118E+02;
  COFD[         693] =   0.4964896032502449E+01;
  COFD[         694] =  -0.4047371090295282E+00;
  COFD[         695] =   0.1652517237194395E-01;
  COFD[         696] =  -0.1607413035852147E+02;
  COFD[         697] =   0.3688389366224249E+01;
  COFD[         698] =  -0.2691169520025368E+00;
  COFD[         699] =   0.1187605069009256E-01;
  COFD[         700] =  -0.1608131536572632E+02;
  COFD[         701] =   0.3688226344868147E+01;
  COFD[         702] =  -0.2690947493495650E+00;
  COFD[         703] =   0.1187504521154676E-01;
  COFD[         704] =  -0.1711740515182854E+02;
  COFD[         705] =   0.4144924979490134E+01;
  COFD[         706] =  -0.3276629227930643E+00;
  COFD[         707] =   0.1438303438727467E-01;
  COFD[         708] =  -0.1711740515182854E+02;
  COFD[         709] =   0.4144924979490134E+01;
  COFD[         710] =  -0.3276629227930643E+00;
  COFD[         711] =   0.1438303438727467E-01;
  COFD[         712] =  -0.1710868501481421E+02;
  COFD[         713] =   0.4131370391002061E+01;
  COFD[         714] =  -0.3258395192522007E+00;
  COFD[         715] =   0.1430143624455670E-01;
  COFD[         716] =  -0.1699844377766270E+02;
  COFD[         717] =   0.4085300971938842E+01;
  COFD[         718] =  -0.3197553288341113E+00;
  COFD[         719] =   0.1403314440725374E-01;
  COFD[         720] =  -0.1584365035133258E+02;
  COFD[         721] =   0.3602236757402370E+01;
  COFD[         722] =  -0.2586073726747737E+00;
  COFD[         723] =   0.1145011549760728E-01;
  COFD[         724] =  -0.1869430855359004E+02;
  COFD[         725] =   0.4562621031385030E+01;
  COFD[         726] =  -0.3751419042160075E+00;
  COFD[         727] =   0.1616247820310002E-01;
  COFD[         728] =  -0.2032584895987123E+02;
  COFD[         729] =   0.5111360106328648E+01;
  COFD[         730] =  -0.4315969740328803E+00;
  COFD[         731] =   0.1801200837165252E-01;
  COFD[         732] =  -0.2033818516187319E+02;
  COFD[         733] =   0.5112803625509553E+01;
  COFD[         734] =  -0.4318111908846191E+00;
  COFD[         735] =   0.1802241389102129E-01;
  COFD[         736] =  -0.1994945707378163E+02;
  COFD[         737] =   0.4988909811915903E+01;
  COFD[         738] =  -0.4200448560138315E+00;
  COFD[         739] =   0.1768418875431581E-01;
  COFD[         740] =  -0.1818905402728113E+02;
  COFD[         741] =   0.4404836108862263E+01;
  COFD[         742] =  -0.3569836886410043E+00;
  COFD[         743] =   0.1547274326407077E-01;
  COFD[         744] =  -0.1819837498757435E+02;
  COFD[         745] =   0.4404368843395893E+01;
  COFD[         746] =  -0.3569335691748086E+00;
  COFD[         747] =   0.1547105514857662E-01;
  COFD[         748] =  -0.1903192515013230E+02;
  COFD[         749] =   0.4693162157411670E+01;
  COFD[         750] =  -0.3899939058159093E+00;
  COFD[         751] =   0.1672093460134929E-01;
  COFD[         752] =  -0.1882911310545702E+02;
  COFD[         753] =   0.4594509279764214E+01;
  COFD[         754] =  -0.3789782571863847E+00;
  COFD[         755] =   0.1631703818360878E-01;
  COFD[         756] =  -0.1883954098397452E+02;
  COFD[         757] =   0.4595024886620513E+01;
  COFD[         758] =  -0.3790444797196689E+00;
  COFD[         759] =   0.1631985649333879E-01;
  COFD[         760] =  -0.1581023736073806E+02;
  COFD[         761] =   0.3592405682636056E+01;
  COFD[         762] =  -0.2573062871053550E+00;
  COFD[         763] =   0.1139264138769127E-01;
  COFD[         764] =  -0.1678197080138354E+02;
  COFD[         765] =   0.3942379447222446E+01;
  COFD[         766] =  -0.3007437798097645E+00;
  COFD[         767] =   0.1319002761262226E-01;
  COFD[         768] =  -0.1355980676188247E+02;
  COFD[         769] =   0.3231182295054575E+01;
  COFD[         770] =  -0.2138406992439825E+00;
  COFD[         771] =   0.9659145053095021E-02;
  COFD[         772] =  -0.1798724510881690E+02;
  COFD[         773] =   0.5079566464178115E+01;
  COFD[         774] =  -0.4444140659462563E+00;
  COFD[         775] =   0.1923462920552580E-01;
  COFD[         776] =  -0.1558834457822415E+02;
  COFD[         777] =   0.3699567664433178E+01;
  COFD[         778] =  -0.2706504848826866E+00;
  COFD[         779] =   0.1194595865045622E-01;
  COFD[         780] =  -0.1707915522153321E+02;
  COFD[         781] =   0.4132970234506635E+01;
  COFD[         782] =  -0.3260547355537204E+00;
  COFD[         783] =   0.1431106723651207E-01;
  COFD[         784] =  -0.1559686009787643E+02;
  COFD[         785] =   0.3697254029246602E+01;
  COFD[         786] =  -0.2703332731709061E+00;
  COFD[         787] =   0.1193149792533076E-01;
  COFD[         788] =  -0.2040913503195636E+02;
  COFD[         789] =   0.5224092639655129E+01;
  COFD[         790] =  -0.4327353018434304E+00;
  COFD[         791] =   0.1753064056853202E-01;
  COFD[         792] =  -0.1709853282395253E+02;
  COFD[         793] =   0.4139002691011624E+01;
  COFD[         794] =  -0.3268662406958736E+00;
  COFD[         795] =   0.1434738242897968E-01;
  COFD[         796] =  -0.1711740515182854E+02;
  COFD[         797] =   0.4144924979490134E+01;
  COFD[         798] =  -0.3276629227930643E+00;
  COFD[         799] =   0.1438303438727467E-01;
  COFD[         800] =  -0.1769521074851900E+02;
  COFD[         801] =   0.4367699379901299E+01;
  COFD[         802] =  -0.3537277912984029E+00;
  COFD[         803] =   0.1539755051214757E-01;
  COFD[         804] =  -0.1769521074851900E+02;
  COFD[         805] =   0.4367699379901299E+01;
  COFD[         806] =  -0.3537277912984029E+00;
  COFD[         807] =   0.1539755051214757E-01;
  COFD[         808] =  -0.1771391298979426E+02;
  COFD[         809] =   0.4368348686795481E+01;
  COFD[         810] =  -0.3538107613054961E+00;
  COFD[         811] =   0.1540107028020538E-01;
  COFD[         812] =  -0.1765053635794619E+02;
  COFD[         813] =   0.4344873087030904E+01;
  COFD[         814] =  -0.3508852749505070E+00;
  COFD[         815] =   0.1527908501173896E-01;
  COFD[         816] =  -0.1661775419297465E+02;
  COFD[         817] =   0.3942050022877351E+01;
  COFD[         818] =  -0.3009976504232300E+00;
  COFD[         819] =   0.1321284000123418E-01;
  COFD[         820] =  -0.1920160126704923E+02;
  COFD[         821] =   0.4768586004236171E+01;
  COFD[         822] =  -0.3938478869794211E+00;
  COFD[         823] =   0.1663350528215831E-01;
  COFD[         824] =  -0.2039806228291693E+02;
  COFD[         825] =   0.5097933786196286E+01;
  COFD[         826] =  -0.4183216752667500E+00;
  COFD[         827] =   0.1696948938817924E-01;
  COFD[         828] =  -0.2038049892591359E+02;
  COFD[         829] =   0.5087331121491625E+01;
  COFD[         830] =  -0.4167206268593121E+00;
  COFD[         831] =   0.1689044439564830E-01;
  COFD[         832] =  -0.2006425198525411E+02;
  COFD[         833] =   0.4999965267895282E+01;
  COFD[         834] =  -0.4103327875262441E+00;
  COFD[         835] =   0.1680582455088142E-01;
  COFD[         836] =  -0.1892750366101823E+02;
  COFD[         837] =   0.4701319071264327E+01;
  COFD[         838] =  -0.3901920115746619E+00;
  COFD[         839] =   0.1669020608030539E-01;
  COFD[         840] =  -0.1893716028173049E+02;
  COFD[         841] =   0.4701920423433434E+01;
  COFD[         842] =  -0.3901723473212924E+00;
  COFD[         843] =   0.1668493495004967E-01;
  COFD[         844] =  -0.1956675647234025E+02;
  COFD[         845] =   0.4893594117189979E+01;
  COFD[         846] =  -0.4083224878273328E+00;
  COFD[         847] =   0.1719480931404233E-01;
  COFD[         848] =  -0.1941109254354329E+02;
  COFD[         849] =   0.4819578196751006E+01;
  COFD[         850] =  -0.4008189828944452E+00;
  COFD[         851] =   0.1695355424415370E-01;
  COFD[         852] =  -0.1941401903486105E+02;
  COFD[         853] =   0.4817662369355102E+01;
  COFD[         854] =  -0.4004301095900044E+00;
  COFD[         855] =   0.1693043453375954E-01;
  COFD[         856] =  -0.1658740614019908E+02;
  COFD[         857] =   0.3933646317245189E+01;
  COFD[         858] =  -0.2999344650396615E+00;
  COFD[         859] =   0.1316797448289397E-01;
  COFD[         860] =  -0.1779539137381903E+02;
  COFD[         861] =   0.4385582017594726E+01;
  COFD[         862] =  -0.3562391667708982E+00;
  COFD[         863] =   0.1551072662807571E-01;
  COFD[         864] =  -0.1355980676188247E+02;
  COFD[         865] =   0.3231182295054575E+01;
  COFD[         866] =  -0.2138406992439825E+00;
  COFD[         867] =   0.9659145053095021E-02;
  COFD[         868] =  -0.1798724510881690E+02;
  COFD[         869] =   0.5079566464178115E+01;
  COFD[         870] =  -0.4444140659462563E+00;
  COFD[         871] =   0.1923462920552580E-01;
  COFD[         872] =  -0.1558834457822415E+02;
  COFD[         873] =   0.3699567664433178E+01;
  COFD[         874] =  -0.2706504848826866E+00;
  COFD[         875] =   0.1194595865045622E-01;
  COFD[         876] =  -0.1707915522153321E+02;
  COFD[         877] =   0.4132970234506635E+01;
  COFD[         878] =  -0.3260547355537204E+00;
  COFD[         879] =   0.1431106723651207E-01;
  COFD[         880] =  -0.1559686009787643E+02;
  COFD[         881] =   0.3697254029246602E+01;
  COFD[         882] =  -0.2703332731709061E+00;
  COFD[         883] =   0.1193149792533076E-01;
  COFD[         884] =  -0.2040913503195636E+02;
  COFD[         885] =   0.5224092639655129E+01;
  COFD[         886] =  -0.4327353018434304E+00;
  COFD[         887] =   0.1753064056853202E-01;
  COFD[         888] =  -0.1709853282395253E+02;
  COFD[         889] =   0.4139002691011624E+01;
  COFD[         890] =  -0.3268662406958736E+00;
  COFD[         891] =   0.1434738242897968E-01;
  COFD[         892] =  -0.1711740515182854E+02;
  COFD[         893] =   0.4144924979490134E+01;
  COFD[         894] =  -0.3276629227930643E+00;
  COFD[         895] =   0.1438303438727467E-01;
  COFD[         896] =  -0.1769521074851900E+02;
  COFD[         897] =   0.4367699379901299E+01;
  COFD[         898] =  -0.3537277912984029E+00;
  COFD[         899] =   0.1539755051214757E-01;
  COFD[         900] =  -0.1769521074851900E+02;
  COFD[         901] =   0.4367699379901299E+01;
  COFD[         902] =  -0.3537277912984029E+00;
  COFD[         903] =   0.1539755051214757E-01;
  COFD[         904] =  -0.1771391298979426E+02;
  COFD[         905] =   0.4368348686795481E+01;
  COFD[         906] =  -0.3538107613054961E+00;
  COFD[         907] =   0.1540107028020538E-01;
  COFD[         908] =  -0.1765053635794619E+02;
  COFD[         909] =   0.4344873087030904E+01;
  COFD[         910] =  -0.3508852749505070E+00;
  COFD[         911] =   0.1527908501173896E-01;
  COFD[         912] =  -0.1661775419297465E+02;
  COFD[         913] =   0.3942050022877351E+01;
  COFD[         914] =  -0.3009976504232300E+00;
  COFD[         915] =   0.1321284000123418E-01;
  COFD[         916] =  -0.1920160126704923E+02;
  COFD[         917] =   0.4768586004236171E+01;
  COFD[         918] =  -0.3938478869794211E+00;
  COFD[         919] =   0.1663350528215831E-01;
  COFD[         920] =  -0.2039806228291693E+02;
  COFD[         921] =   0.5097933786196286E+01;
  COFD[         922] =  -0.4183216752667500E+00;
  COFD[         923] =   0.1696948938817924E-01;
  COFD[         924] =  -0.2038049892591359E+02;
  COFD[         925] =   0.5087331121491625E+01;
  COFD[         926] =  -0.4167206268593121E+00;
  COFD[         927] =   0.1689044439564830E-01;
  COFD[         928] =  -0.2006425198525411E+02;
  COFD[         929] =   0.4999965267895282E+01;
  COFD[         930] =  -0.4103327875262441E+00;
  COFD[         931] =   0.1680582455088142E-01;
  COFD[         932] =  -0.1892750366101823E+02;
  COFD[         933] =   0.4701319071264327E+01;
  COFD[         934] =  -0.3901920115746619E+00;
  COFD[         935] =   0.1669020608030539E-01;
  COFD[         936] =  -0.1893716028173049E+02;
  COFD[         937] =   0.4701920423433434E+01;
  COFD[         938] =  -0.3901723473212924E+00;
  COFD[         939] =   0.1668493495004967E-01;
  COFD[         940] =  -0.1956675647234025E+02;
  COFD[         941] =   0.4893594117189979E+01;
  COFD[         942] =  -0.4083224878273328E+00;
  COFD[         943] =   0.1719480931404233E-01;
  COFD[         944] =  -0.1941109254354329E+02;
  COFD[         945] =   0.4819578196751006E+01;
  COFD[         946] =  -0.4008189828944452E+00;
  COFD[         947] =   0.1695355424415370E-01;
  COFD[         948] =  -0.1941401903486105E+02;
  COFD[         949] =   0.4817662369355102E+01;
  COFD[         950] =  -0.4004301095900044E+00;
  COFD[         951] =   0.1693043453375954E-01;
  COFD[         952] =  -0.1658740614019908E+02;
  COFD[         953] =   0.3933646317245189E+01;
  COFD[         954] =  -0.2999344650396615E+00;
  COFD[         955] =   0.1316797448289397E-01;
  COFD[         956] =  -0.1779539137381903E+02;
  COFD[         957] =   0.4385582017594726E+01;
  COFD[         958] =  -0.3562391667708982E+00;
  COFD[         959] =   0.1551072662807571E-01;
  COFD[         960] =  -0.1358415095584097E+02;
  COFD[         961] =   0.3239845970896248E+01;
  COFD[         962] =  -0.2150500406377150E+00;
  COFD[         963] =   0.9715225271588263E-02;
  COFD[         964] =  -0.1801636950310948E+02;
  COFD[         965] =   0.5090088857460223E+01;
  COFD[         966] =  -0.4457552498826564E+00;
  COFD[         967] =   0.1929141580461785E-01;
  COFD[         968] =  -0.1561687122378518E+02;
  COFD[         969] =   0.3703886367486989E+01;
  COFD[         970] =  -0.2712410029515650E+00;
  COFD[         971] =   0.1197280722366742E-01;
  COFD[         972] =  -0.1707094794267795E+02;
  COFD[         973] =   0.4119810624341889E+01;
  COFD[         974] =  -0.3242844616145015E+00;
  COFD[         975] =   0.1423184749658152E-01;
  COFD[         976] =  -0.1562210625534589E+02;
  COFD[         977] =   0.3699988818188416E+01;
  COFD[         978] =  -0.2707081263160739E+00;
  COFD[         979] =   0.1194858184298051E-01;
  COFD[         980] =  -0.2042453573019913E+02;
  COFD[         981] =   0.5222435968901098E+01;
  COFD[         982] =  -0.4324987976330752E+00;
  COFD[         983] =   0.1751958312757915E-01;
  COFD[         984] =  -0.1709003104221693E+02;
  COFD[         985] =   0.4125627998285617E+01;
  COFD[         986] =  -0.3250670331906144E+00;
  COFD[         987] =   0.1426686739847163E-01;
  COFD[         988] =  -0.1710868501481421E+02;
  COFD[         989] =   0.4131370391002061E+01;
  COFD[         990] =  -0.3258395192522007E+00;
  COFD[         991] =   0.1430143624455670E-01;
  COFD[         992] =  -0.1771391298979426E+02;
  COFD[         993] =   0.4368348686795481E+01;
  COFD[         994] =  -0.3538107613054961E+00;
  COFD[         995] =   0.1540107028020538E-01;
  COFD[         996] =  -0.1771391298979426E+02;
  COFD[         997] =   0.4368348686795481E+01;
  COFD[         998] =  -0.3538107613054961E+00;
  COFD[         999] =   0.1540107028020538E-01;
  COFD[        1000] =  -0.1772990792056323E+02;
  COFD[        1001] =   0.4367699379901269E+01;
  COFD[        1002] =  -0.3537277912983985E+00;
  COFD[        1003] =   0.1539755051214736E-01;
  COFD[        1004] =  -0.1766459414605863E+02;
  COFD[        1005] =   0.4343242779501765E+01;
  COFD[        1006] =  -0.3506769601296933E+00;
  COFD[        1007] =   0.1527024825917116E-01;
  COFD[        1008] =  -0.1661214673917502E+02;
  COFD[        1009] =   0.3930381940303202E+01;
  COFD[        1010] =  -0.2994310262293675E+00;
  COFD[        1011] =   0.1314286546846328E-01;
  COFD[        1012] =  -0.1923136095279614E+02;
  COFD[        1013] =   0.4771652981887103E+01;
  COFD[        1014] =  -0.3944968044729790E+00;
  COFD[        1015] =   0.1667283377652949E-01;
  COFD[        1016] =  -0.2046646394524606E+02;
  COFD[        1017] =   0.5118752767847491E+01;
  COFD[        1018] =  -0.4214659652131379E+00;
  COFD[        1019] =   0.1712474984073264E-01;
  COFD[        1020] =  -0.2045023591963550E+02;
  COFD[        1021] =   0.5108641161759671E+01;
  COFD[        1022] =  -0.4199387136777197E+00;
  COFD[        1023] =   0.1704933199723114E-01;
  COFD[        1024] =  -0.2012688925964889E+02;
  COFD[        1025] =   0.5018247199767933E+01;
  COFD[        1026] =  -0.4131414445091446E+00;
  COFD[        1027] =   0.1694646037906594E-01;
  COFD[        1028] =  -0.1894432968047615E+02;
  COFD[        1029] =   0.4700246655831890E+01;
  COFD[        1030] =  -0.3902232238988751E+00;
  COFD[        1031] =   0.1669926672613195E-01;
  COFD[        1032] =  -0.1895414515078359E+02;
  COFD[        1033] =   0.4700828894373993E+01;
  COFD[        1034] =  -0.3902069149863984E+00;
  COFD[        1035] =   0.1669440360391662E-01;
  COFD[        1036] =  -0.1960184896363511E+02;
  COFD[        1037] =   0.4900346283733346E+01;
  COFD[        1038] =  -0.4094945328875742E+00;
  COFD[        1039] =   0.1725872606725634E-01;
  COFD[        1040] =  -0.1943925299002717E+02;
  COFD[        1041] =   0.4823356728877616E+01;
  COFD[        1042] =  -0.4015874038077185E+00;
  COFD[        1043] =   0.1699928308073373E-01;
  COFD[        1044] =  -0.1944256474839717E+02;
  COFD[        1045] =   0.4821522129341849E+01;
  COFD[        1046] =  -0.4012140535427826E+00;
  COFD[        1047] =   0.1697705724216803E-01;
  COFD[        1048] =  -0.1658229952875101E+02;
  COFD[        1049] =   0.3922184822536529E+01;
  COFD[        1050] =  -0.2983959485115160E+00;
  COFD[        1051] =   0.1309927190981370E-01;
  COFD[        1052] =  -0.1779710613097560E+02;
  COFD[        1053] =   0.4376269742117769E+01;
  COFD[        1054] =  -0.3550500828076509E+00;
  COFD[        1055] =   0.1546030915508793E-01;
  COFD[        1056] =  -0.1353824314289343E+02;
  COFD[        1057] =   0.3227352939332374E+01;
  COFD[        1058] =  -0.2134245159235849E+00;
  COFD[        1059] =   0.9644648995332943E-02;
  COFD[        1060] =  -0.1790236652645725E+02;
  COFD[        1061] =   0.5054241512059871E+01;
  COFD[        1062] =  -0.4413163862361035E+00;
  COFD[        1063] =   0.1910841720545356E-01;
  COFD[        1064] =  -0.1556217965348612E+02;
  COFD[        1065] =   0.3683628994574988E+01;
  COFD[        1066] =  -0.2686556693855442E+00;
  COFD[        1067] =   0.1186253607235680E-01;
  COFD[        1068] =  -0.1696131964138645E+02;
  COFD[        1069] =   0.4074167708212651E+01;
  COFD[        1070] =  -0.3182576642303228E+00;
  COFD[        1071] =   0.1396612551510442E-01;
  COFD[        1072] =  -0.1556532980903788E+02;
  COFD[        1073] =   0.3678661452138608E+01;
  COFD[        1074] =  -0.2679777825820168E+00;
  COFD[        1075] =   0.1183177506717497E-01;
  COFD[        1076] =  -0.2052802769939090E+02;
  COFD[        1077] =   0.5190809253464677E+01;
  COFD[        1078] =  -0.4197156907449909E+00;
  COFD[        1079] =   0.1661902490158712E-01;
  COFD[        1080] =  -0.1698006191382480E+02;
  COFD[        1081] =   0.4079755891690070E+01;
  COFD[        1082] =  -0.3190093940831370E+00;
  COFD[        1083] =   0.1399976445173436E-01;
  COFD[        1084] =  -0.1699844377766270E+02;
  COFD[        1085] =   0.4085300971938842E+01;
  COFD[        1086] =  -0.3197553288341113E+00;
  COFD[        1087] =   0.1403314440725374E-01;
  COFD[        1088] =  -0.1765053635794619E+02;
  COFD[        1089] =   0.4344873087030904E+01;
  COFD[        1090] =  -0.3508852749505070E+00;
  COFD[        1091] =   0.1527908501173896E-01;
  COFD[        1092] =  -0.1765053635794619E+02;
  COFD[        1093] =   0.4344873087030904E+01;
  COFD[        1094] =  -0.3508852749505070E+00;
  COFD[        1095] =   0.1527908501173896E-01;
  COFD[        1096] =  -0.1766459414605863E+02;
  COFD[        1097] =   0.4343242779501765E+01;
  COFD[        1098] =  -0.3506769601296933E+00;
  COFD[        1099] =   0.1527024825917116E-01;
  COFD[        1100] =  -0.1760280207090635E+02;
  COFD[        1101] =   0.4320037087778842E+01;
  COFD[        1102] =  -0.3477962090917468E+00;
  COFD[        1103] =   0.1515062790890059E-01;
  COFD[        1104] =  -0.1654504240089281E+02;
  COFD[        1105] =   0.3902990496709335E+01;
  COFD[        1106] =  -0.2959726841289695E+00;
  COFD[        1107] =   0.1299731318228341E-01;
  COFD[        1108] =  -0.1920918070133062E+02;
  COFD[        1109] =   0.4764713252836613E+01;
  COFD[        1110] =  -0.3942405534249802E+00;
  COFD[        1111] =   0.1668923743050566E-01;
  COFD[        1112] =  -0.2048733868393580E+02;
  COFD[        1113] =   0.5132889365248658E+01;
  COFD[        1114] =  -0.4243039026916984E+00;
  COFD[        1115] =   0.1728905581961421E-01;
  COFD[        1116] =  -0.2047278277859408E+02;
  COFD[        1117] =   0.5123414930187060E+01;
  COFD[        1118] =  -0.4228706903648477E+00;
  COFD[        1119] =   0.1721819830097161E-01;
  COFD[        1120] =  -0.2034321636694622E+02;
  COFD[        1121] =   0.5089562300279502E+01;
  COFD[        1122] =  -0.4211891334428361E+00;
  COFD[        1123] =   0.1724909778947470E-01;
  COFD[        1124] =  -0.1889022870587310E+02;
  COFD[        1125] =   0.4680203446785792E+01;
  COFD[        1126] =  -0.3880194721331493E+00;
  COFD[        1127] =   0.1662057070374977E-01;
  COFD[        1128] =  -0.1890039213728353E+02;
  COFD[        1129] =   0.4680852507758593E+01;
  COFD[        1130] =  -0.3880184127527163E+00;
  COFD[        1131] =   0.1661666470021389E-01;
  COFD[        1132] =  -0.1957015754946189E+02;
  COFD[        1133] =   0.4890650507507935E+01;
  COFD[        1134] =  -0.4088565156459115E+00;
  COFD[        1135] =   0.1725720932243744E-01;
  COFD[        1136] =  -0.1938496483487797E+02;
  COFD[        1137] =   0.4803651837482498E+01;
  COFD[        1138] =  -0.3995406103851766E+00;
  COFD[        1139] =   0.1693221082596063E-01;
  COFD[        1140] =  -0.1938900888397246E+02;
  COFD[        1141] =   0.4802051950424429E+01;
  COFD[        1142] =  -0.3992044152028926E+00;
  COFD[        1143] =   0.1691189193038190E-01;
  COFD[        1144] =  -0.1651845590523054E+02;
  COFD[        1145] =   0.3896132878391479E+01;
  COFD[        1146] =  -0.2951170694006329E+00;
  COFD[        1147] =   0.1296170338310411E-01;
  COFD[        1148] =  -0.1772757250170383E+02;
  COFD[        1149] =   0.4347725932428951E+01;
  COFD[        1150] =  -0.3515015798062593E+00;
  COFD[        1151] =   0.1531308129463464E-01;
  COFD[        1152] =  -0.1268296956267843E+02;
  COFD[        1153] =   0.2903428401599904E+01;
  COFD[        1154] =  -0.1707667302851327E+00;
  COFD[        1155] =   0.7771111802295727E-02;
  COFD[        1156] =  -0.1719132602269712E+02;
  COFD[        1157] =   0.4859367175853364E+01;
  COFD[        1158] =  -0.4240580935572132E+00;
  COFD[        1159] =   0.1870696403708276E-01;
  COFD[        1160] =  -0.1478198899238034E+02;
  COFD[        1161] =   0.3371092581023257E+01;
  COFD[        1162] =  -0.2298567856863946E+00;
  COFD[        1163] =   0.1025846107153025E-01;
  COFD[        1164] =  -0.1582244348422915E+02;
  COFD[        1165] =   0.3599211636009877E+01;
  COFD[        1166] =  -0.2581898369286856E+00;
  COFD[        1167] =   0.1143096575439008E-01;
  COFD[        1168] =  -0.1476985219999908E+02;
  COFD[        1169] =   0.3357751196108267E+01;
  COFD[        1170] =  -0.2280046831004087E+00;
  COFD[        1171] =   0.1017303667696902E-01;
  COFD[        1172] =  -0.2031234004461199E+02;
  COFD[        1173] =   0.5186312454417227E+01;
  COFD[        1174] =  -0.4313742876779729E+00;
  COFD[        1175] =   0.1759974844922337E-01;
  COFD[        1176] =  -0.1583284545979817E+02;
  COFD[        1177] =   0.3600569993944359E+01;
  COFD[        1178] =  -0.2583773478248265E+00;
  COFD[        1179] =   0.1143956672273590E-01;
  COFD[        1180] =  -0.1584365035133258E+02;
  COFD[        1181] =   0.3602236757402370E+01;
  COFD[        1182] =  -0.2586073726747737E+00;
  COFD[        1183] =   0.1145011549760728E-01;
  COFD[        1184] =  -0.1661775419297465E+02;
  COFD[        1185] =   0.3942050022877351E+01;
  COFD[        1186] =  -0.3009976504232300E+00;
  COFD[        1187] =   0.1321284000123418E-01;
  COFD[        1188] =  -0.1661775419297465E+02;
  COFD[        1189] =   0.3942050022877351E+01;
  COFD[        1190] =  -0.3009976504232300E+00;
  COFD[        1191] =   0.1321284000123418E-01;
  COFD[        1192] =  -0.1661214673917502E+02;
  COFD[        1193] =   0.3930381940303202E+01;
  COFD[        1194] =  -0.2994310262293675E+00;
  COFD[        1195] =   0.1314286546846328E-01;
  COFD[        1196] =  -0.1654504240089281E+02;
  COFD[        1197] =   0.3902990496709335E+01;
  COFD[        1198] =  -0.2959726841289695E+00;
  COFD[        1199] =   0.1299731318228341E-01;
  COFD[        1200] =  -0.1551400389810888E+02;
  COFD[        1201] =   0.3472777644794381E+01;
  COFD[        1202] =  -0.2416432013010369E+00;
  COFD[        1203] =   0.1070799930256162E-01;
  COFD[        1204] =  -0.1853079633836986E+02;
  COFD[        1205] =   0.4516962151189352E+01;
  COFD[        1206] =  -0.3708427592293742E+00;
  COFD[        1207] =   0.1604368344993687E-01;
  COFD[        1208] =  -0.2009021453031469E+02;
  COFD[        1209] =   0.5037634138971256E+01;
  COFD[        1210] =  -0.4241196608853510E+00;
  COFD[        1211] =   0.1777020759153531E-01;
  COFD[        1212] =  -0.2009778046496907E+02;
  COFD[        1213] =   0.5037181430781024E+01;
  COFD[        1214] =  -0.4240362270814515E+00;
  COFD[        1215] =   0.1776546650877937E-01;
  COFD[        1216] =  -0.2003301781547066E+02;
  COFD[        1217] =   0.5028625403684593E+01;
  COFD[        1218] =  -0.4259384280972174E+00;
  COFD[        1219] =   0.1797007932483027E-01;
  COFD[        1220] =  -0.1801976002924096E+02;
  COFD[        1221] =   0.4352352585107287E+01;
  COFD[        1222] =  -0.3517965558706893E+00;
  COFD[        1223] =   0.1531621848405548E-01;
  COFD[        1224] =  -0.1803148637733784E+02;
  COFD[        1225] =   0.4353159047512093E+01;
  COFD[        1226] =  -0.3519044208218940E+00;
  COFD[        1227] =   0.1532101211937375E-01;
  COFD[        1228] =  -0.1879103964706305E+02;
  COFD[        1229] =   0.4613167070471964E+01;
  COFD[        1230] =  -0.3811821259847005E+00;
  COFD[        1231] =   0.1640353922053438E-01;
  COFD[        1232] =  -0.1864571977877140E+02;
  COFD[        1233] =   0.4539069966095800E+01;
  COFD[        1234] =  -0.3735011850255855E+00;
  COFD[        1235] =   0.1615133908644331E-01;
  COFD[        1236] =  -0.1865835471859333E+02;
  COFD[        1237] =   0.4540734987509840E+01;
  COFD[        1238] =  -0.3737078963922989E+00;
  COFD[        1239] =   0.1615982140755382E-01;
  COFD[        1240] =  -0.1549017806948574E+02;
  COFD[        1241] =   0.3466859355125401E+01;
  COFD[        1242] =  -0.2408856054712291E+00;
  COFD[        1243] =   0.1067561900190753E-01;
  COFD[        1244] =  -0.1650781936547825E+02;
  COFD[        1245] =   0.3842122413055806E+01;
  COFD[        1246] =  -0.2882098997185555E+00;
  COFD[        1247] =   0.1266705263344362E-01;
  COFD[        1248] =  -0.1564739122681956E+02;
  COFD[        1249] =   0.4025825276486324E+01;
  COFD[        1250] =  -0.3181195147607377E+00;
  COFD[        1251] =   0.1422136468533644E-01;
  COFD[        1252] =  -0.1725034548774479E+02;
  COFD[        1253] =   0.4566921339972978E+01;
  COFD[        1254] =  -0.3509192271095338E+00;
  COFD[        1255] =   0.1402685086906688E-01;
  COFD[        1256] =  -0.1793673861488134E+02;
  COFD[        1257] =   0.4489875914439151E+01;
  COFD[        1258] =  -0.3697556337020395E+00;
  COFD[        1259] =   0.1609253879718407E-01;
  COFD[        1260] =  -0.1868514242925209E+02;
  COFD[        1261] =   0.4565512377700616E+01;
  COFD[        1262] =  -0.3754538805580755E+00;
  COFD[        1263] =   0.1617305773699574E-01;
  COFD[        1264] =  -0.1792454027766008E+02;
  COFD[        1265] =   0.4476326984714260E+01;
  COFD[        1266] =  -0.3680191278135618E+00;
  COFD[        1267] =   0.1601860672433683E-01;
  COFD[        1268] =  -0.1924815690645235E+02;
  COFD[        1269] =   0.4356731855020855E+01;
  COFD[        1270] =  -0.2835834611821879E+00;
  COFD[        1271] =   0.9647123135613308E-02;
  COFD[        1272] =  -0.1868972720602865E+02;
  COFD[        1273] =   0.4564003101584411E+01;
  COFD[        1274] =  -0.3752914154261894E+00;
  COFD[        1275] =   0.1616756912994680E-01;
  COFD[        1276] =  -0.1869430855359004E+02;
  COFD[        1277] =   0.4562621031385030E+01;
  COFD[        1278] =  -0.3751419042160075E+00;
  COFD[        1279] =   0.1616247820310002E-01;
  COFD[        1280] =  -0.1920160126704923E+02;
  COFD[        1281] =   0.4768586004236171E+01;
  COFD[        1282] =  -0.3938478869794211E+00;
  COFD[        1283] =   0.1663350528215831E-01;
  COFD[        1284] =  -0.1920160126704923E+02;
  COFD[        1285] =   0.4768586004236171E+01;
  COFD[        1286] =  -0.3938478869794211E+00;
  COFD[        1287] =   0.1663350528215831E-01;
  COFD[        1288] =  -0.1923136095279614E+02;
  COFD[        1289] =   0.4771652981887103E+01;
  COFD[        1290] =  -0.3944968044729790E+00;
  COFD[        1291] =   0.1667283377652949E-01;
  COFD[        1292] =  -0.1920918070133062E+02;
  COFD[        1293] =   0.4764713252836613E+01;
  COFD[        1294] =  -0.3942405534249802E+00;
  COFD[        1295] =   0.1668923743050566E-01;
  COFD[        1296] =  -0.1853079633836986E+02;
  COFD[        1297] =   0.4516962151189352E+01;
  COFD[        1298] =  -0.3708427592293742E+00;
  COFD[        1299] =   0.1604368344993687E-01;
  COFD[        1300] =  -0.2074682089064003E+02;
  COFD[        1301] =   0.5125913401864712E+01;
  COFD[        1302] =  -0.4300036496008145E+00;
  COFD[        1303] =   0.1780187794631543E-01;
  COFD[        1304] =  -0.2103480542370173E+02;
  COFD[        1305] =   0.5058663047911905E+01;
  COFD[        1306] =  -0.3955876986371539E+00;
  COFD[        1307] =   0.1531365269600822E-01;
  COFD[        1308] =  -0.2106203504941583E+02;
  COFD[        1309] =   0.5066201505062354E+01;
  COFD[        1310] =  -0.3966848046994042E+00;
  COFD[        1311] =   0.1536587788350550E-01;
  COFD[        1312] =  -0.2111193727907456E+02;
  COFD[        1313] =   0.5115004803331022E+01;
  COFD[        1314] =  -0.4073081193984433E+00;
  COFD[        1315] =   0.1598596435486438E-01;
  COFD[        1316] =  -0.2039137163285755E+02;
  COFD[        1317] =   0.5048551827610241E+01;
  COFD[        1318] =  -0.4237554975633671E+00;
  COFD[        1319] =   0.1768105924766743E-01;
  COFD[        1320] =  -0.2041424604605752E+02;
  COFD[        1321] =   0.5053848971414935E+01;
  COFD[        1322] =  -0.4245860964146182E+00;
  COFD[        1323] =   0.1772329162485299E-01;
  COFD[        1324] =  -0.2089328018920450E+02;
  COFD[        1325] =   0.5166930585415199E+01;
  COFD[        1326] =  -0.4306289589876920E+00;
  COFD[        1327] =   0.1764267219152353E-01;
  COFD[        1328] =  -0.2075862463488640E+02;
  COFD[        1329] =   0.5108417296497037E+01;
  COFD[        1330] =  -0.4262851069696628E+00;
  COFD[        1331] =   0.1758280414396813E-01;
  COFD[        1332] =  -0.2077840705218999E+02;
  COFD[        1333] =   0.5112879719806333E+01;
  COFD[        1334] =  -0.4269664397363305E+00;
  COFD[        1335] =   0.1761673548420319E-01;
  COFD[        1336] =  -0.1850244542172855E+02;
  COFD[        1337] =   0.4509546506902997E+01;
  COFD[        1338] =  -0.3699250247632380E+00;
  COFD[        1339] =   0.1600569104209367E-01;
  COFD[        1340] =  -0.1931016351705685E+02;
  COFD[        1341] =   0.4757139666603682E+01;
  COFD[        1342] =  -0.3961700036764624E+00;
  COFD[        1343] =   0.1690081368592412E-01;
  COFD[        1344] =  -0.1746999275651299E+02;
  COFD[        1345] =   0.4670556529920976E+01;
  COFD[        1346] =  -0.3931788949833612E+00;
  COFD[        1347] =   0.1709979614643272E-01;
  COFD[        1348] =  -0.1540645026491837E+02;
  COFD[        1349] =   0.3572427121876953E+01;
  COFD[        1350] =  -0.1878763535992870E+00;
  COFD[        1351] =   0.5591984032298197E-02;
  COFD[        1352] =  -0.1915075646861156E+02;
  COFD[        1353] =   0.4868838487330766E+01;
  COFD[        1354] =  -0.4046210142562804E+00;
  COFD[        1355] =   0.1701122981099066E-01;
  COFD[        1356] =  -0.2031657960333580E+02;
  COFD[        1357] =   0.5113614640951361E+01;
  COFD[        1358] =  -0.4319293017089550E+00;
  COFD[        1359] =   0.1802805482147992E-01;
  COFD[        1360] =  -0.1918204714711579E+02;
  COFD[        1361] =   0.4875615818554563E+01;
  COFD[        1362] =  -0.4058262469897604E+00;
  COFD[        1363] =   0.1707794755285167E-01;
  COFD[        1364] =  -0.1757312087019186E+02;
  COFD[        1365] =   0.3511204681563826E+01;
  COFD[        1366] =  -0.1497966127823777E+00;
  COFD[        1367] =   0.2966580974304371E-02;
  COFD[        1368] =  -0.2032173485009814E+02;
  COFD[        1369] =   0.5112656067238100E+01;
  COFD[        1370] =  -0.4317894593025599E+00;
  COFD[        1371] =   0.1802136540679988E-01;
  COFD[        1372] =  -0.2032584895987123E+02;
  COFD[        1373] =   0.5111360106328648E+01;
  COFD[        1374] =  -0.4315969740328803E+00;
  COFD[        1375] =   0.1801200837165252E-01;
  COFD[        1376] =  -0.2039806228291693E+02;
  COFD[        1377] =   0.5097933786196286E+01;
  COFD[        1378] =  -0.4183216752667500E+00;
  COFD[        1379] =   0.1696948938817924E-01;
  COFD[        1380] =  -0.2039806228291693E+02;
  COFD[        1381] =   0.5097933786196286E+01;
  COFD[        1382] =  -0.4183216752667500E+00;
  COFD[        1383] =   0.1696948938817924E-01;
  COFD[        1384] =  -0.2046646394524606E+02;
  COFD[        1385] =   0.5118752767847491E+01;
  COFD[        1386] =  -0.4214659652131379E+00;
  COFD[        1387] =   0.1712474984073264E-01;
  COFD[        1388] =  -0.2048733868393580E+02;
  COFD[        1389] =   0.5132889365248658E+01;
  COFD[        1390] =  -0.4243039026916984E+00;
  COFD[        1391] =   0.1728905581961421E-01;
  COFD[        1392] =  -0.2009021453031469E+02;
  COFD[        1393] =   0.5037634138971256E+01;
  COFD[        1394] =  -0.4241196608853510E+00;
  COFD[        1395] =   0.1777020759153531E-01;
  COFD[        1396] =  -0.2103480542370173E+02;
  COFD[        1397] =   0.5058663047911905E+01;
  COFD[        1398] =  -0.3955876986371539E+00;
  COFD[        1399] =   0.1531365269600822E-01;
  COFD[        1400] =  -0.1885553227873480E+02;
  COFD[        1401] =   0.3923352162852614E+01;
  COFD[        1402] =  -0.2115865709699218E+00;
  COFD[        1403] =   0.5939285843445688E-02;
  COFD[        1404] =  -0.1886381013480732E+02;
  COFD[        1405] =   0.3923279625734948E+01;
  COFD[        1406] =  -0.2115777534001486E+00;
  COFD[        1407] =   0.5938970665371634E-02;
  COFD[        1408] =  -0.1981646458797120E+02;
  COFD[        1409] =   0.4378886175261416E+01;
  COFD[        1410] =  -0.2813828834600260E+00;
  COFD[        1411] =   0.9370857616443682E-02;
  COFD[        1412] =  -0.2113317649345642E+02;
  COFD[        1413] =   0.5173175010557953E+01;
  COFD[        1414] =  -0.4173234769896741E+00;
  COFD[        1415] =   0.1651910472359747E-01;
  COFD[        1416] =  -0.2114401307743194E+02;
  COFD[        1417] =   0.5173607724136159E+01;
  COFD[        1418] =  -0.4173904988942925E+00;
  COFD[        1419] =   0.1652249796434591E-01;
  COFD[        1420] =  -0.2090933405715619E+02;
  COFD[        1421] =   0.4972999527197874E+01;
  COFD[        1422] =  -0.3786889263628220E+00;
  COFD[        1423] =   0.1435806171655258E-01;
  COFD[        1424] =  -0.2110606940488274E+02;
  COFD[        1425] =   0.5056386072237322E+01;
  COFD[        1426] =  -0.3941314241059736E+00;
  COFD[        1427] =   0.1520386988001248E-01;
  COFD[        1428] =  -0.2110836259086301E+02;
  COFD[        1429] =   0.5053674020415143E+01;
  COFD[        1430] =  -0.3937388243323189E+00;
  COFD[        1431] =   0.1518528106281032E-01;
  COFD[        1432] =  -0.2006577079410498E+02;
  COFD[        1433] =   0.5032536447034316E+01;
  COFD[        1434] =  -0.4235925902580250E+00;
  COFD[        1435] =   0.1775279409683131E-01;
  COFD[        1436] =  -0.2070496369832456E+02;
  COFD[        1437] =   0.5191072368701859E+01;
  COFD[        1438] =  -0.4346910627906122E+00;
  COFD[        1439] =   0.1785842811134005E-01;
  COFD[        1440] =  -0.1747771970561251E+02;
  COFD[        1441] =   0.4673159679696371E+01;
  COFD[        1442] =  -0.3935109715079382E+00;
  COFD[        1443] =   0.1711386921544412E-01;
  COFD[        1444] =  -0.1539084008249366E+02;
  COFD[        1445] =   0.3565006575777335E+01;
  COFD[        1446] =  -0.1867580729432450E+00;
  COFD[        1447] =   0.5536884066371989E-02;
  COFD[        1448] =  -0.1914986770399175E+02;
  COFD[        1449] =   0.4864899621234062E+01;
  COFD[        1450] =  -0.4039220297382060E+00;
  COFD[        1451] =   0.1697258735720344E-01;
  COFD[        1452] =  -0.2032689236579443E+02;
  COFD[        1453] =   0.5114221260960875E+01;
  COFD[        1454] =  -0.4320136109530544E+00;
  COFD[        1455] =   0.1803190494311217E-01;
  COFD[        1456] =  -0.1918177070837826E+02;
  COFD[        1457] =   0.4871873199573439E+01;
  COFD[        1458] =  -0.4051602536432213E+00;
  COFD[        1459] =   0.1704106539271711E-01;
  COFD[        1460] =  -0.1757470287826831E+02;
  COFD[        1461] =   0.3509436887842453E+01;
  COFD[        1462] =  -0.1496102927025162E+00;
  COFD[        1463] =   0.2961974572337444E-02;
  COFD[        1464] =  -0.2033308000404864E+02;
  COFD[        1465] =   0.5113689932626547E+01;
  COFD[        1466] =  -0.4319400850127593E+00;
  COFD[        1467] =   0.1802856188319298E-01;
  COFD[        1468] =  -0.2033818516187319E+02;
  COFD[        1469] =   0.5112803625509553E+01;
  COFD[        1470] =  -0.4318111908846191E+00;
  COFD[        1471] =   0.1802241389102129E-01;
  COFD[        1472] =  -0.2038049892591359E+02;
  COFD[        1473] =   0.5087331121491625E+01;
  COFD[        1474] =  -0.4167206268593121E+00;
  COFD[        1475] =   0.1689044439564830E-01;
  COFD[        1476] =  -0.2038049892591359E+02;
  COFD[        1477] =   0.5087331121491625E+01;
  COFD[        1478] =  -0.4167206268593121E+00;
  COFD[        1479] =   0.1689044439564830E-01;
  COFD[        1480] =  -0.2045023591963550E+02;
  COFD[        1481] =   0.5108641161759671E+01;
  COFD[        1482] =  -0.4199387136777197E+00;
  COFD[        1483] =   0.1704933199723114E-01;
  COFD[        1484] =  -0.2047278277859408E+02;
  COFD[        1485] =   0.5123414930187060E+01;
  COFD[        1486] =  -0.4228706903648477E+00;
  COFD[        1487] =   0.1721819830097161E-01;
  COFD[        1488] =  -0.2009778046496907E+02;
  COFD[        1489] =   0.5037181430781024E+01;
  COFD[        1490] =  -0.4240362270814515E+00;
  COFD[        1491] =   0.1776546650877937E-01;
  COFD[        1492] =  -0.2106203504941583E+02;
  COFD[        1493] =   0.5066201505062354E+01;
  COFD[        1494] =  -0.3966848046994042E+00;
  COFD[        1495] =   0.1536587788350550E-01;
  COFD[        1496] =  -0.1886381013480732E+02;
  COFD[        1497] =   0.3923279625734948E+01;
  COFD[        1498] =  -0.2115777534001486E+00;
  COFD[        1499] =   0.5938970665371634E-02;
  COFD[        1500] =  -0.1887260515065304E+02;
  COFD[        1501] =   0.3923352162852601E+01;
  COFD[        1502] =  -0.2115865709699200E+00;
  COFD[        1503] =   0.5939285843445600E-02;
  COFD[        1504] =  -0.1982650591715289E+02;
  COFD[        1505] =   0.4379424471801094E+01;
  COFD[        1506] =  -0.2814553408580349E+00;
  COFD[        1507] =   0.9373962429419786E-02;
  COFD[        1508] =  -0.2113855609134408E+02;
  COFD[        1509] =   0.5171981200778053E+01;
  COFD[        1510] =  -0.4171447495559704E+00;
  COFD[        1511] =   0.1651034132923998E-01;
  COFD[        1512] =  -0.2115144768087362E+02;
  COFD[        1513] =   0.5173261676775666E+01;
  COFD[        1514] =  -0.4173365706093256E+00;
  COFD[        1515] =   0.1651975240232191E-01;
  COFD[        1516] =  -0.2091775085800750E+02;
  COFD[        1517] =   0.4973029974129495E+01;
  COFD[        1518] =  -0.3786913370630364E+00;
  COFD[        1519] =   0.1435807083846824E-01;
  COFD[        1520] =  -0.2111907495038380E+02;
  COFD[        1521] =   0.5058381671032987E+01;
  COFD[        1522] =  -0.3944196800959486E+00;
  COFD[        1523] =   0.1521748580862574E-01;
  COFD[        1524] =  -0.2112313758034541E+02;
  COFD[        1525] =   0.5056389549388699E+01;
  COFD[        1526] =  -0.3941319269977342E+00;
  COFD[        1527] =   0.1520389366699413E-01;
  COFD[        1528] =  -0.2007323696432831E+02;
  COFD[        1529] =   0.5032033200707714E+01;
  COFD[        1530] =  -0.4235009227577490E+00;
  COFD[        1531] =   0.1774762291235686E-01;
  COFD[        1532] =  -0.2072332061257778E+02;
  COFD[        1533] =   0.5194960177957292E+01;
  COFD[        1534] =  -0.4352726149968372E+00;
  COFD[        1535] =   0.1788688372961993E-01;
  COFD[        1536] =  -0.1797854335218850E+02;
  COFD[        1537] =   0.4905020138101169E+01;
  COFD[        1538] =  -0.4293084201751438E+00;
  COFD[        1539] =   0.1891302999394448E-01;
  COFD[        1540] =  -0.1561125740910447E+02;
  COFD[        1541] =   0.3698144089263720E+01;
  COFD[        1542] =  -0.2103080207621357E+00;
  COFD[        1543] =   0.6792204226221277E-02;
  COFD[        1544] =  -0.1878355293596396E+02;
  COFD[        1545] =   0.4745970603113269E+01;
  COFD[        1546] =  -0.3924372163577848E+00;
  COFD[        1547] =   0.1663659729501466E-01;
  COFD[        1548] =  -0.2016013281728147E+02;
  COFD[        1549] =   0.5063566265537697E+01;
  COFD[        1550] =  -0.4281515792120217E+00;
  COFD[        1551] =   0.1797180686456706E-01;
  COFD[        1552] =  -0.1880409131690931E+02;
  COFD[        1553] =   0.4747896766189557E+01;
  COFD[        1554] =  -0.3929513558466455E+00;
  COFD[        1555] =   0.1667071927589917E-01;
  COFD[        1556] =  -0.1619587514529141E+02;
  COFD[        1557] =   0.2843235032975332E+01;
  COFD[        1558] =  -0.5104314541884596E-01;
  COFD[        1559] =  -0.1707937587486378E-02;
  COFD[        1560] =  -0.1994316533928780E+02;
  COFD[        1561] =   0.4989293390342089E+01;
  COFD[        1562] =  -0.4200959859596722E+00;
  COFD[        1563] =   0.1768643162186107E-01;
  COFD[        1564] =  -0.1994945707378163E+02;
  COFD[        1565] =   0.4988909811915903E+01;
  COFD[        1566] =  -0.4200448560138315E+00;
  COFD[        1567] =   0.1768418875431581E-01;
  COFD[        1568] =  -0.2006425198525411E+02;
  COFD[        1569] =   0.4999965267895282E+01;
  COFD[        1570] =  -0.4103327875262441E+00;
  COFD[        1571] =   0.1680582455088142E-01;
  COFD[        1572] =  -0.2006425198525411E+02;
  COFD[        1573] =   0.4999965267895282E+01;
  COFD[        1574] =  -0.4103327875262441E+00;
  COFD[        1575] =   0.1680582455088142E-01;
  COFD[        1576] =  -0.2012688925964889E+02;
  COFD[        1577] =   0.5018247199767933E+01;
  COFD[        1578] =  -0.4131414445091446E+00;
  COFD[        1579] =   0.1694646037906594E-01;
  COFD[        1580] =  -0.2034321636694622E+02;
  COFD[        1581] =   0.5089562300279502E+01;
  COFD[        1582] =  -0.4211891334428361E+00;
  COFD[        1583] =   0.1724909778947470E-01;
  COFD[        1584] =  -0.2003301781547066E+02;
  COFD[        1585] =   0.5028625403684593E+01;
  COFD[        1586] =  -0.4259384280972174E+00;
  COFD[        1587] =   0.1797007932483027E-01;
  COFD[        1588] =  -0.2111193727907456E+02;
  COFD[        1589] =   0.5115004803331022E+01;
  COFD[        1590] =  -0.4073081193984433E+00;
  COFD[        1591] =   0.1598596435486438E-01;
  COFD[        1592] =  -0.1981646458797120E+02;
  COFD[        1593] =   0.4378886175261416E+01;
  COFD[        1594] =  -0.2813828834600260E+00;
  COFD[        1595] =   0.9370857616443682E-02;
  COFD[        1596] =  -0.1982650591715289E+02;
  COFD[        1597] =   0.4379424471801094E+01;
  COFD[        1598] =  -0.2814553408580349E+00;
  COFD[        1599] =   0.9373962429419786E-02;
  COFD[        1600] =  -0.1906244143057827E+02;
  COFD[        1601] =   0.4049325847301449E+01;
  COFD[        1602] =  -0.2359825286728792E+00;
  COFD[        1603] =   0.7307807686876559E-02;
  COFD[        1604] =  -0.2112929323949077E+02;
  COFD[        1605] =   0.5220748890969238E+01;
  COFD[        1606] =  -0.4305058272023291E+00;
  COFD[        1607] =   0.1736476920050096E-01;
  COFD[        1608] =  -0.2114446337698859E+02;
  COFD[        1609] =   0.5222999206283971E+01;
  COFD[        1610] =  -0.4308431817261412E+00;
  COFD[        1611] =   0.1738131867299051E-01;
  COFD[        1612] =  -0.2114547150354106E+02;
  COFD[        1613] =   0.5122846710851252E+01;
  COFD[        1614] =  -0.4061105011633706E+00;
  COFD[        1615] =   0.1585327078701023E-01;
  COFD[        1616] =  -0.2123091991916349E+02;
  COFD[        1617] =   0.5158619088859661E+01;
  COFD[        1618] =  -0.4148430725527501E+00;
  COFD[        1619] =   0.1638779761015646E-01;
  COFD[        1620] =  -0.2123706489334075E+02;
  COFD[        1621] =   0.5157473227475789E+01;
  COFD[        1622] =  -0.4146762094312990E+00;
  COFD[        1623] =   0.1637983966492620E-01;
  COFD[        1624] =  -0.2000088490858467E+02;
  COFD[        1625] =   0.5021789747731799E+01;
  COFD[        1626] =  -0.4253480595637520E+00;
  COFD[        1627] =   0.1795667045275812E-01;
  COFD[        1628] =  -0.2043431836412989E+02;
  COFD[        1629] =   0.5116486808224121E+01;
  COFD[        1630] =  -0.4301998553589693E+00;
  COFD[        1631] =   0.1786838693497468E-01;
  COFD[        1632] =  -0.1528921665984086E+02;
  COFD[        1633] =   0.3867752822044799E+01;
  COFD[        1634] =  -0.2988043102116158E+00;
  COFD[        1635] =   0.1343910541069478E-01;
  COFD[        1636] =  -0.1774507254557799E+02;
  COFD[        1637] =   0.4780159131489816E+01;
  COFD[        1638] =  -0.3856587269477316E+00;
  COFD[        1639] =   0.1582942762638978E-01;
  COFD[        1640] =  -0.1731509452354737E+02;
  COFD[        1641] =   0.4259553692469185E+01;
  COFD[        1642] =  -0.3417709853226915E+00;
  COFD[        1643] =   0.1496136243888912E-01;
  COFD[        1644] =  -0.1817374730686181E+02;
  COFD[        1645] =   0.4404209750065649E+01;
  COFD[        1646] =  -0.3569181269134391E+00;
  COFD[        1647] =   0.1547062305483025E-01;
  COFD[        1648] =  -0.1729594517981932E+02;
  COFD[        1649] =   0.4244246410490659E+01;
  COFD[        1650] =  -0.3397463840897515E+00;
  COFD[        1651] =   0.1487228746520592E-01;
  COFD[        1652] =  -0.2024506268677845E+02;
  COFD[        1653] =   0.4905829896937226E+01;
  COFD[        1654] =  -0.3737705712203094E+00;
  COFD[        1655] =   0.1428535217779750E-01;
  COFD[        1656] =  -0.1818135452916891E+02;
  COFD[        1657] =   0.4404445109750770E+01;
  COFD[        1658] =  -0.3569414786935936E+00;
  COFD[        1659] =   0.1547130685514495E-01;
  COFD[        1660] =  -0.1818905402728113E+02;
  COFD[        1661] =   0.4404836108862263E+01;
  COFD[        1662] =  -0.3569836886410043E+00;
  COFD[        1663] =   0.1547274326407077E-01;
  COFD[        1664] =  -0.1892750366101823E+02;
  COFD[        1665] =   0.4701319071264327E+01;
  COFD[        1666] =  -0.3901920115746619E+00;
  COFD[        1667] =   0.1669020608030539E-01;
  COFD[        1668] =  -0.1892750366101823E+02;
  COFD[        1669] =   0.4701319071264327E+01;
  COFD[        1670] =  -0.3901920115746619E+00;
  COFD[        1671] =   0.1669020608030539E-01;
  COFD[        1672] =  -0.1894432968047615E+02;
  COFD[        1673] =   0.4700246655831890E+01;
  COFD[        1674] =  -0.3902232238988751E+00;
  COFD[        1675] =   0.1669926672613195E-01;
  COFD[        1676] =  -0.1889022870587310E+02;
  COFD[        1677] =   0.4680203446785792E+01;
  COFD[        1678] =  -0.3880194721331493E+00;
  COFD[        1679] =   0.1662057070374977E-01;
  COFD[        1680] =  -0.1801976002924096E+02;
  COFD[        1681] =   0.4352352585107287E+01;
  COFD[        1682] =  -0.3517965558706893E+00;
  COFD[        1683] =   0.1531621848405548E-01;
  COFD[        1684] =  -0.2039137163285755E+02;
  COFD[        1685] =   0.5048551827610241E+01;
  COFD[        1686] =  -0.4237554975633671E+00;
  COFD[        1687] =   0.1768105924766743E-01;
  COFD[        1688] =  -0.2113317649345642E+02;
  COFD[        1689] =   0.5173175010557953E+01;
  COFD[        1690] =  -0.4173234769896741E+00;
  COFD[        1691] =   0.1651910472359747E-01;
  COFD[        1692] =  -0.2113855609134408E+02;
  COFD[        1693] =   0.5171981200778053E+01;
  COFD[        1694] =  -0.4171447495559704E+00;
  COFD[        1695] =   0.1651034132923998E-01;
  COFD[        1696] =  -0.2112929323949077E+02;
  COFD[        1697] =   0.5220748890969238E+01;
  COFD[        1698] =  -0.4305058272023291E+00;
  COFD[        1699] =   0.1736476920050096E-01;
  COFD[        1700] =  -0.2008313017084402E+02;
  COFD[        1701] =   0.4983892868408141E+01;
  COFD[        1702] =  -0.4201561978851936E+00;
  COFD[        1703] =   0.1772039680073913E-01;
  COFD[        1704] =  -0.2009217954894486E+02;
  COFD[        1705] =   0.4983712068565962E+01;
  COFD[        1706] =  -0.4201265590632418E+00;
  COFD[        1707] =   0.1771884028887708E-01;
  COFD[        1708] =  -0.2062921491889363E+02;
  COFD[        1709] =   0.5125110770555073E+01;
  COFD[        1710] =  -0.4303770100880492E+00;
  COFD[        1711] =   0.1783747461488234E-01;
  COFD[        1712] =  -0.2059391054107079E+02;
  COFD[        1713] =   0.5108068575679129E+01;
  COFD[        1714] =  -0.4315347958033334E+00;
  COFD[        1715] =   0.1802344633255497E-01;
  COFD[        1716] =  -0.2059735929037355E+02;
  COFD[        1717] =   0.5105897247055966E+01;
  COFD[        1718] =  -0.4311951951062819E+00;
  COFD[        1719] =   0.1800620766882382E-01;
  COFD[        1720] =  -0.1799101249332176E+02;
  COFD[        1721] =   0.4344595334851394E+01;
  COFD[        1722] =  -0.3508294711072030E+00;
  COFD[        1723] =   0.1527589170814564E-01;
  COFD[        1724] =  -0.1894826439081699E+02;
  COFD[        1725] =   0.4665904241842780E+01;
  COFD[        1726] =  -0.3874653193184607E+00;
  COFD[        1727] =   0.1665191589170439E-01;
  COFD[        1728] =  -0.1530201153719761E+02;
  COFD[        1729] =   0.3872623565106621E+01;
  COFD[        1730] =  -0.2994802960421091E+00;
  COFD[        1731] =   0.1347028558095244E-01;
  COFD[        1732] =  -0.1774834356244655E+02;
  COFD[        1733] =   0.4780626617651506E+01;
  COFD[        1734] =  -0.3856327873318030E+00;
  COFD[        1735] =   0.1582440706500046E-01;
  COFD[        1736] =  -0.1734665788392033E+02;
  COFD[        1737] =   0.4269427891180465E+01;
  COFD[        1738] =  -0.3430766910038187E+00;
  COFD[        1739] =   0.1501879627721415E-01;
  COFD[        1740] =  -0.1818389332535593E+02;
  COFD[        1741] =   0.4404158898744668E+01;
  COFD[        1742] =  -0.3569180632138942E+00;
  COFD[        1743] =   0.1547091340944365E-01;
  COFD[        1744] =  -0.1732667376192460E+02;
  COFD[        1745] =   0.4253674558730024E+01;
  COFD[        1746] =  -0.3409934577271315E+00;
  COFD[        1747] =   0.1492715703864712E-01;
  COFD[        1748] =  -0.2021382909884265E+02;
  COFD[        1749] =   0.4888560960839369E+01;
  COFD[        1750] =  -0.3712523102934615E+00;
  COFD[        1751] =   0.1416513273755347E-01;
  COFD[        1752] =  -0.1819107613580373E+02;
  COFD[        1753] =   0.4404180695515474E+01;
  COFD[        1754] =  -0.3569157439691602E+00;
  COFD[        1755] =   0.1547058259292524E-01;
  COFD[        1756] =  -0.1819837498757435E+02;
  COFD[        1757] =   0.4404368843395893E+01;
  COFD[        1758] =  -0.3569335691748086E+00;
  COFD[        1759] =   0.1547105514857662E-01;
  COFD[        1760] =  -0.1893716028173049E+02;
  COFD[        1761] =   0.4701920423433434E+01;
  COFD[        1762] =  -0.3901723473212924E+00;
  COFD[        1763] =   0.1668493495004967E-01;
  COFD[        1764] =  -0.1893716028173049E+02;
  COFD[        1765] =   0.4701920423433434E+01;
  COFD[        1766] =  -0.3901723473212924E+00;
  COFD[        1767] =   0.1668493495004967E-01;
  COFD[        1768] =  -0.1895414515078359E+02;
  COFD[        1769] =   0.4700828894373993E+01;
  COFD[        1770] =  -0.3902069149863984E+00;
  COFD[        1771] =   0.1669440360391662E-01;
  COFD[        1772] =  -0.1890039213728353E+02;
  COFD[        1773] =   0.4680852507758593E+01;
  COFD[        1774] =  -0.3880184127527163E+00;
  COFD[        1775] =   0.1661666470021389E-01;
  COFD[        1776] =  -0.1803148637733784E+02;
  COFD[        1777] =   0.4353159047512093E+01;
  COFD[        1778] =  -0.3519044208218940E+00;
  COFD[        1779] =   0.1532101211937375E-01;
  COFD[        1780] =  -0.2041424604605752E+02;
  COFD[        1781] =   0.5053848971414935E+01;
  COFD[        1782] =  -0.4245860964146182E+00;
  COFD[        1783] =   0.1772329162485299E-01;
  COFD[        1784] =  -0.2114401307743194E+02;
  COFD[        1785] =   0.5173607724136159E+01;
  COFD[        1786] =  -0.4173904988942925E+00;
  COFD[        1787] =   0.1652249796434591E-01;
  COFD[        1788] =  -0.2115144768087362E+02;
  COFD[        1789] =   0.5173261676775666E+01;
  COFD[        1790] =  -0.4173365706093256E+00;
  COFD[        1791] =   0.1651975240232191E-01;
  COFD[        1792] =  -0.2114446337698859E+02;
  COFD[        1793] =   0.5222999206283971E+01;
  COFD[        1794] =  -0.4308431817261412E+00;
  COFD[        1795] =   0.1738131867299051E-01;
  COFD[        1796] =  -0.2009217954894486E+02;
  COFD[        1797] =   0.4983712068565962E+01;
  COFD[        1798] =  -0.4201265590632418E+00;
  COFD[        1799] =   0.1771884028887708E-01;
  COFD[        1800] =  -0.2010212049883026E+02;
  COFD[        1801] =   0.4983892868408190E+01;
  COFD[        1802] =  -0.4201561978852008E+00;
  COFD[        1803] =   0.1772039680073947E-01;
  COFD[        1804] =  -0.2064026808809792E+02;
  COFD[        1805] =   0.5125729516119874E+01;
  COFD[        1806] =  -0.4304737187234308E+00;
  COFD[        1807] =   0.1784238527988932E-01;
  COFD[        1808] =  -0.2060778849004208E+02;
  COFD[        1809] =   0.5109962431730637E+01;
  COFD[        1810] =  -0.4318312376520599E+00;
  COFD[        1811] =   0.1803850360308551E-01;
  COFD[        1812] =  -0.2061243965594226E+02;
  COFD[        1813] =   0.5108287844008140E+01;
  COFD[        1814] =  -0.4315691026948430E+00;
  COFD[        1815] =   0.1802518831107918E-01;
  COFD[        1816] =  -0.1800289314001846E+02;
  COFD[        1817] =   0.4345462384297714E+01;
  COFD[        1818] =  -0.3509451050542303E+00;
  COFD[        1819] =   0.1528101631398434E-01;
  COFD[        1820] =  -0.1895758159603200E+02;
  COFD[        1821] =   0.4665279128753897E+01;
  COFD[        1822] =  -0.3874162488122672E+00;
  COFD[        1823] =   0.1665125384367453E-01;
  COFD[        1824] =  -0.1646041046594534E+02;
  COFD[        1825] =   0.4329097783536894E+01;
  COFD[        1826] =  -0.3588538378541796E+00;
  COFD[        1827] =   0.1604603061706265E-01;
  COFD[        1828] =  -0.1694497223361135E+02;
  COFD[        1829] =   0.4353290721870846E+01;
  COFD[        1830] =  -0.3147441391485086E+00;
  COFD[        1831] =   0.1210201943729390E-01;
  COFD[        1832] =  -0.1800091724620567E+02;
  COFD[        1833] =   0.4493072014688184E+01;
  COFD[        1834] =  -0.3678451189139861E+00;
  COFD[        1835] =   0.1591448159607913E-01;
  COFD[        1836] =  -0.1901929468754062E+02;
  COFD[        1837] =   0.4693839754776671E+01;
  COFD[        1838] =  -0.3901000075815114E+00;
  COFD[        1839] =   0.1672631492507185E-01;
  COFD[        1840] =  -0.1799417924090561E+02;
  COFD[        1841] =   0.4483237269191199E+01;
  COFD[        1842] =  -0.3666457816700467E+00;
  COFD[        1843] =   0.1586624138953086E-01;
  COFD[        1844] =  -0.1988688758541057E+02;
  COFD[        1845] =   0.4655348588622681E+01;
  COFD[        1846] =  -0.3287200044949488E+00;
  COFD[        1847] =   0.1186464366263145E-01;
  COFD[        1848] =  -0.1902570873205879E+02;
  COFD[        1849] =   0.4693487292163261E+01;
  COFD[        1850] =  -0.3900471625097403E+00;
  COFD[        1851] =   0.1672372332121293E-01;
  COFD[        1852] =  -0.1903192515013230E+02;
  COFD[        1853] =   0.4693162157411670E+01;
  COFD[        1854] =  -0.3899939058159093E+00;
  COFD[        1855] =   0.1672093460134929E-01;
  COFD[        1856] =  -0.1956675647234025E+02;
  COFD[        1857] =   0.4893594117189979E+01;
  COFD[        1858] =  -0.4083224878273328E+00;
  COFD[        1859] =   0.1719480931404233E-01;
  COFD[        1860] =  -0.1956675647234025E+02;
  COFD[        1861] =   0.4893594117189979E+01;
  COFD[        1862] =  -0.4083224878273328E+00;
  COFD[        1863] =   0.1719480931404233E-01;
  COFD[        1864] =  -0.1960184896363511E+02;
  COFD[        1865] =   0.4900346283733346E+01;
  COFD[        1866] =  -0.4094945328875742E+00;
  COFD[        1867] =   0.1725872606725634E-01;
  COFD[        1868] =  -0.1957015754946189E+02;
  COFD[        1869] =   0.4890650507507935E+01;
  COFD[        1870] =  -0.4088565156459115E+00;
  COFD[        1871] =   0.1725720932243744E-01;
  COFD[        1872] =  -0.1879103964706305E+02;
  COFD[        1873] =   0.4613167070471964E+01;
  COFD[        1874] =  -0.3811821259847005E+00;
  COFD[        1875] =   0.1640353922053438E-01;
  COFD[        1876] =  -0.2089328018920450E+02;
  COFD[        1877] =   0.5166930585415199E+01;
  COFD[        1878] =  -0.4306289589876920E+00;
  COFD[        1879] =   0.1764267219152353E-01;
  COFD[        1880] =  -0.2090933405715619E+02;
  COFD[        1881] =   0.4972999527197874E+01;
  COFD[        1882] =  -0.3786889263628220E+00;
  COFD[        1883] =   0.1435806171655258E-01;
  COFD[        1884] =  -0.2091775085800750E+02;
  COFD[        1885] =   0.4973029974129495E+01;
  COFD[        1886] =  -0.3786913370630364E+00;
  COFD[        1887] =   0.1435807083846824E-01;
  COFD[        1888] =  -0.2114547150354106E+02;
  COFD[        1889] =   0.5122846710851252E+01;
  COFD[        1890] =  -0.4061105011633706E+00;
  COFD[        1891] =   0.1585327078701023E-01;
  COFD[        1892] =  -0.2062921491889363E+02;
  COFD[        1893] =   0.5125110770555073E+01;
  COFD[        1894] =  -0.4303770100880492E+00;
  COFD[        1895] =   0.1783747461488234E-01;
  COFD[        1896] =  -0.2064026808809792E+02;
  COFD[        1897] =   0.5125729516119874E+01;
  COFD[        1898] =  -0.4304737187234308E+00;
  COFD[        1899] =   0.1784238527988932E-01;
  COFD[        1900] =  -0.2117733141272616E+02;
  COFD[        1901] =   0.5262507718327111E+01;
  COFD[        1902] =  -0.4398579014089885E+00;
  COFD[        1903] =   0.1792489830320274E-01;
  COFD[        1904] =  -0.2108189705348579E+02;
  COFD[        1905] =   0.5219096733932689E+01;
  COFD[        1906] =  -0.4371166395377749E+00;
  COFD[        1907] =   0.1791475976287751E-01;
  COFD[        1908] =  -0.2108550712837031E+02;
  COFD[        1909] =   0.5216945829110916E+01;
  COFD[        1910] =  -0.4367926562187473E+00;
  COFD[        1911] =   0.1789880186533132E-01;
  COFD[        1912] =  -0.1876345278233955E+02;
  COFD[        1913] =   0.4606235123260915E+01;
  COFD[        1914] =  -0.3803641770348614E+00;
  COFD[        1915] =   0.1637150008428244E-01;
  COFD[        1916] =  -0.1962450162993201E+02;
  COFD[        1917] =   0.4875751864229406E+01;
  COFD[        1918] =  -0.4087050411778894E+00;
  COFD[        1919] =   0.1732743672805724E-01;
  COFD[        1920] =  -0.1615267929807825E+02;
  COFD[        1921] =   0.4174121490551778E+01;
  COFD[        1922] =  -0.3383475624805202E+00;
  COFD[        1923] =   0.1513993577326728E-01;
  COFD[        1924] =  -0.1738386031480564E+02;
  COFD[        1925] =   0.4515894896789344E+01;
  COFD[        1926] =  -0.3391856065220307E+00;
  COFD[        1927] =   0.1329276720561771E-01;
  COFD[        1928] =  -0.1802433609839966E+02;
  COFD[        1929] =   0.4483646303442399E+01;
  COFD[        1930] =  -0.3688151334706006E+00;
  COFD[        1931] =   0.1604826293273014E-01;
  COFD[        1932] =  -0.1881741117557303E+02;
  COFD[        1933] =   0.4595519619722809E+01;
  COFD[        1934] =  -0.3791053358880415E+00;
  COFD[        1935] =   0.1632232422675412E-01;
  COFD[        1936] =  -0.1800542957417654E+02;
  COFD[        1937] =   0.4468504271448043E+01;
  COFD[        1938] =  -0.3668692513471102E+00;
  COFD[        1939] =   0.1596519376748332E-01;
  COFD[        1940] =  -0.2010033725040245E+02;
  COFD[        1941] =   0.4740608584978043E+01;
  COFD[        1942] =  -0.3442903017689571E+00;
  COFD[        1943] =   0.1271241178230842E-01;
  COFD[        1944] =  -0.1882326100183751E+02;
  COFD[        1945] =   0.4594957876626870E+01;
  COFD[        1946] =  -0.3790360690538448E+00;
  COFD[        1947] =   0.1631950746180167E-01;
  COFD[        1948] =  -0.1882911310545702E+02;
  COFD[        1949] =   0.4594509279764214E+01;
  COFD[        1950] =  -0.3789782571863847E+00;
  COFD[        1951] =   0.1631703818360878E-01;
  COFD[        1952] =  -0.1941109254354329E+02;
  COFD[        1953] =   0.4819578196751006E+01;
  COFD[        1954] =  -0.4008189828944452E+00;
  COFD[        1955] =   0.1695355424415370E-01;
  COFD[        1956] =  -0.1941109254354329E+02;
  COFD[        1957] =   0.4819578196751006E+01;
  COFD[        1958] =  -0.4008189828944452E+00;
  COFD[        1959] =   0.1695355424415370E-01;
  COFD[        1960] =  -0.1943925299002717E+02;
  COFD[        1961] =   0.4823356728877616E+01;
  COFD[        1962] =  -0.4015874038077185E+00;
  COFD[        1963] =   0.1699928308073373E-01;
  COFD[        1964] =  -0.1938496483487797E+02;
  COFD[        1965] =   0.4803651837482498E+01;
  COFD[        1966] =  -0.3995406103851766E+00;
  COFD[        1967] =   0.1693221082596063E-01;
  COFD[        1968] =  -0.1864571977877140E+02;
  COFD[        1969] =   0.4539069966095800E+01;
  COFD[        1970] =  -0.3735011850255855E+00;
  COFD[        1971] =   0.1615133908644331E-01;
  COFD[        1972] =  -0.2075862463488640E+02;
  COFD[        1973] =   0.5108417296497037E+01;
  COFD[        1974] =  -0.4262851069696628E+00;
  COFD[        1975] =   0.1758280414396813E-01;
  COFD[        1976] =  -0.2110606940488274E+02;
  COFD[        1977] =   0.5056386072237322E+01;
  COFD[        1978] =  -0.3941314241059736E+00;
  COFD[        1979] =   0.1520386988001248E-01;
  COFD[        1980] =  -0.2111907495038380E+02;
  COFD[        1981] =   0.5058381671032987E+01;
  COFD[        1982] =  -0.3944196800959486E+00;
  COFD[        1983] =   0.1521748580862574E-01;
  COFD[        1984] =  -0.2123091991916349E+02;
  COFD[        1985] =   0.5158619088859661E+01;
  COFD[        1986] =  -0.4148430725527501E+00;
  COFD[        1987] =   0.1638779761015646E-01;
  COFD[        1988] =  -0.2059391054107079E+02;
  COFD[        1989] =   0.5108068575679129E+01;
  COFD[        1990] =  -0.4315347958033334E+00;
  COFD[        1991] =   0.1802344633255497E-01;
  COFD[        1992] =  -0.2060778849004208E+02;
  COFD[        1993] =   0.5109962431730637E+01;
  COFD[        1994] =  -0.4318312376520599E+00;
  COFD[        1995] =   0.1803850360308551E-01;
  COFD[        1996] =  -0.2108189705348579E+02;
  COFD[        1997] =   0.5219096733932689E+01;
  COFD[        1998] =  -0.4371166395377749E+00;
  COFD[        1999] =   0.1791475976287751E-01;
  COFD[        2000] =  -0.2091185026658962E+02;
  COFD[        2001] =   0.5148560230790555E+01;
  COFD[        2002] =  -0.4310562392852861E+00;
  COFD[        2003] =   0.1777358978417338E-01;
  COFD[        2004] =  -0.2091965964592144E+02;
  COFD[        2005] =   0.5148262389934279E+01;
  COFD[        2006] =  -0.4310108894089904E+00;
  COFD[        2007] =   0.1777133583281057E-01;
  COFD[        2008] =  -0.1861866210245701E+02;
  COFD[        2009] =   0.4532060302179115E+01;
  COFD[        2010] =  -0.3726516669477872E+00;
  COFD[        2011] =   0.1611699967368039E-01;
  COFD[        2012] =  -0.1943393371703941E+02;
  COFD[        2013] =   0.4784664368892784E+01;
  COFD[        2014] =  -0.3991077435475853E+00;
  COFD[        2015] =   0.1700210403998681E-01;
  COFD[        2016] =  -0.1616515483252968E+02;
  COFD[        2017] =   0.4178910660883203E+01;
  COFD[        2018] =  -0.3390062877759687E+00;
  COFD[        2019] =   0.1517006334795817E-01;
  COFD[        2020] =  -0.1738223879329547E+02;
  COFD[        2021] =   0.4514299549461490E+01;
  COFD[        2022] =  -0.3388651702152571E+00;
  COFD[        2023] =   0.1327383213284630E-01;
  COFD[        2024] =  -0.1805221831789452E+02;
  COFD[        2025] =   0.4492279525932131E+01;
  COFD[        2026] =  -0.3699242263001061E+00;
  COFD[        2027] =   0.1609559465701492E-01;
  COFD[        2028] =  -0.1882824574536174E+02;
  COFD[        2029] =   0.4596254296231404E+01;
  COFD[        2030] =  -0.3791928797888692E+00;
  COFD[        2031] =   0.1632573982847240E-01;
  COFD[        2032] =  -0.1803293840758895E+02;
  COFD[        2033] =   0.4476898762718537E+01;
  COFD[        2034] =  -0.3679481179513786E+00;
  COFD[        2035] =   0.1601125468656841E-01;
  COFD[        2036] =  -0.2007381993778366E+02;
  COFD[        2037] =   0.4726146050811121E+01;
  COFD[        2038] =  -0.3422109775805109E+00;
  COFD[        2039] =   0.1261488904944177E-01;
  COFD[        2040] =  -0.1883388798798641E+02;
  COFD[        2041] =   0.4595580957493091E+01;
  COFD[        2042] =  -0.3791127554800276E+00;
  COFD[        2043] =   0.1632261913288707E-01;
  COFD[        2044] =  -0.1883954098397452E+02;
  COFD[        2045] =   0.4595024886620513E+01;
  COFD[        2046] =  -0.3790444797196689E+00;
  COFD[        2047] =   0.1631985649333879E-01;
  COFD[        2048] =  -0.1941401903486105E+02;
  COFD[        2049] =   0.4817662369355102E+01;
  COFD[        2050] =  -0.4004301095900044E+00;
  COFD[        2051] =   0.1693043453375954E-01;
  COFD[        2052] =  -0.1941401903486105E+02;
  COFD[        2053] =   0.4817662369355102E+01;
  COFD[        2054] =  -0.4004301095900044E+00;
  COFD[        2055] =   0.1693043453375954E-01;
  COFD[        2056] =  -0.1944256474839717E+02;
  COFD[        2057] =   0.4821522129341849E+01;
  COFD[        2058] =  -0.4012140535427826E+00;
  COFD[        2059] =   0.1697705724216803E-01;
  COFD[        2060] =  -0.1938900888397246E+02;
  COFD[        2061] =   0.4802051950424429E+01;
  COFD[        2062] =  -0.3992044152028926E+00;
  COFD[        2063] =   0.1691189193038190E-01;
  COFD[        2064] =  -0.1865835471859333E+02;
  COFD[        2065] =   0.4540734987509840E+01;
  COFD[        2066] =  -0.3737078963922989E+00;
  COFD[        2067] =   0.1615982140755382E-01;
  COFD[        2068] =  -0.2077840705218999E+02;
  COFD[        2069] =   0.5112879719806333E+01;
  COFD[        2070] =  -0.4269664397363305E+00;
  COFD[        2071] =   0.1761673548420319E-01;
  COFD[        2072] =  -0.2110836259086301E+02;
  COFD[        2073] =   0.5053674020415143E+01;
  COFD[        2074] =  -0.3937388243323189E+00;
  COFD[        2075] =   0.1518528106281032E-01;
  COFD[        2076] =  -0.2112313758034541E+02;
  COFD[        2077] =   0.5056389549388699E+01;
  COFD[        2078] =  -0.3941319269977342E+00;
  COFD[        2079] =   0.1520389366699413E-01;
  COFD[        2080] =  -0.2123706489334075E+02;
  COFD[        2081] =   0.5157473227475789E+01;
  COFD[        2082] =  -0.4146762094312990E+00;
  COFD[        2083] =   0.1637983966492620E-01;
  COFD[        2084] =  -0.2059735929037355E+02;
  COFD[        2085] =   0.5105897247055966E+01;
  COFD[        2086] =  -0.4311951951062819E+00;
  COFD[        2087] =   0.1800620766882382E-01;
  COFD[        2088] =  -0.2061243965594226E+02;
  COFD[        2089] =   0.5108287844008140E+01;
  COFD[        2090] =  -0.4315691026948430E+00;
  COFD[        2091] =   0.1802518831107918E-01;
  COFD[        2092] =  -0.2108550712837031E+02;
  COFD[        2093] =   0.5216945829110916E+01;
  COFD[        2094] =  -0.4367926562187473E+00;
  COFD[        2095] =   0.1789880186533132E-01;
  COFD[        2096] =  -0.2091965964592144E+02;
  COFD[        2097] =   0.5148262389934279E+01;
  COFD[        2098] =  -0.4310108894089904E+00;
  COFD[        2099] =   0.1777133583281057E-01;
  COFD[        2100] =  -0.2092889793954542E+02;
  COFD[        2101] =   0.5148560230790550E+01;
  COFD[        2102] =  -0.4310562392852854E+00;
  COFD[        2103] =   0.1777358978417334E-01;
  COFD[        2104] =  -0.1863141958605046E+02;
  COFD[        2105] =   0.4533772315478532E+01;
  COFD[        2106] =  -0.3728642283605221E+00;
  COFD[        2107] =   0.1612572322702682E-01;
  COFD[        2108] =  -0.1944451603679417E+02;
  COFD[        2109] =   0.4785108067913721E+01;
  COFD[        2110] =  -0.3991858464067674E+00;
  COFD[        2111] =   0.1700639619802289E-01;
  COFD[        2112] =  -0.1266066516595179E+02;
  COFD[        2113] =   0.2898076057146818E+01;
  COFD[        2114] =  -0.1700453759318706E+00;
  COFD[        2115] =   0.7738690296212197E-02;
  COFD[        2116] =  -0.1711567017943538E+02;
  COFD[        2117] =   0.4833264586465570E+01;
  COFD[        2118] =  -0.4205908440542898E+00;
  COFD[        2119] =   0.1855345683322859E-01;
  COFD[        2120] =  -0.1474819919893963E+02;
  COFD[        2121] =   0.3361311502126538E+01;
  COFD[        2122] =  -0.2285465363406641E+00;
  COFD[        2123] =   0.1019990809984005E-01;
  COFD[        2124] =  -0.1578873902175925E+02;
  COFD[        2125] =   0.3589258101555060E+01;
  COFD[        2126] =  -0.2568718883472193E+00;
  COFD[        2127] =   0.1137271971440306E-01;
  COFD[        2128] =  -0.1473661449825940E+02;
  COFD[        2129] =   0.3348204558833826E+01;
  COFD[        2130] =  -0.2267271723233180E+00;
  COFD[        2131] =   0.1011600240359858E-01;
  COFD[        2132] =  -0.2024331300258610E+02;
  COFD[        2133] =   0.5167058816761965E+01;
  COFD[        2134] =  -0.4292618020154789E+00;
  COFD[        2135] =   0.1752039422653532E-01;
  COFD[        2136] =  -0.1579929195048582E+02;
  COFD[        2137] =   0.3590679879226839E+01;
  COFD[        2138] =  -0.2570681340587627E+00;
  COFD[        2139] =   0.1138172060899465E-01;
  COFD[        2140] =  -0.1581023736073806E+02;
  COFD[        2141] =   0.3592405682636056E+01;
  COFD[        2142] =  -0.2573062871053550E+00;
  COFD[        2143] =   0.1139264138769127E-01;
  COFD[        2144] =  -0.1658740614019908E+02;
  COFD[        2145] =   0.3933646317245189E+01;
  COFD[        2146] =  -0.2999344650396615E+00;
  COFD[        2147] =   0.1316797448289397E-01;
  COFD[        2148] =  -0.1658740614019908E+02;
  COFD[        2149] =   0.3933646317245189E+01;
  COFD[        2150] =  -0.2999344650396615E+00;
  COFD[        2151] =   0.1316797448289397E-01;
  COFD[        2152] =  -0.1658229952875101E+02;
  COFD[        2153] =   0.3922184822536529E+01;
  COFD[        2154] =  -0.2983959485115160E+00;
  COFD[        2155] =   0.1309927190981370E-01;
  COFD[        2156] =  -0.1651845590523054E+02;
  COFD[        2157] =   0.3896132878391479E+01;
  COFD[        2158] =  -0.2951170694006329E+00;
  COFD[        2159] =   0.1296170338310411E-01;
  COFD[        2160] =  -0.1549017806948574E+02;
  COFD[        2161] =   0.3466859355125401E+01;
  COFD[        2162] =  -0.2408856054712291E+00;
  COFD[        2163] =   0.1067561900190753E-01;
  COFD[        2164] =  -0.1850244542172855E+02;
  COFD[        2165] =   0.4509546506902997E+01;
  COFD[        2166] =  -0.3699250247632380E+00;
  COFD[        2167] =   0.1600569104209367E-01;
  COFD[        2168] =  -0.2006577079410498E+02;
  COFD[        2169] =   0.5032536447034316E+01;
  COFD[        2170] =  -0.4235925902580250E+00;
  COFD[        2171] =   0.1775279409683131E-01;
  COFD[        2172] =  -0.2007323696432831E+02;
  COFD[        2173] =   0.5032033200707714E+01;
  COFD[        2174] =  -0.4235009227577490E+00;
  COFD[        2175] =   0.1774762291235686E-01;
  COFD[        2176] =  -0.2000088490858467E+02;
  COFD[        2177] =   0.5021789747731799E+01;
  COFD[        2178] =  -0.4253480595637520E+00;
  COFD[        2179] =   0.1795667045275812E-01;
  COFD[        2180] =  -0.1799101249332176E+02;
  COFD[        2181] =   0.4344595334851394E+01;
  COFD[        2182] =  -0.3508294711072030E+00;
  COFD[        2183] =   0.1527589170814564E-01;
  COFD[        2184] =  -0.1800289314001846E+02;
  COFD[        2185] =   0.4345462384297714E+01;
  COFD[        2186] =  -0.3509451050542303E+00;
  COFD[        2187] =   0.1528101631398434E-01;
  COFD[        2188] =  -0.1876345278233955E+02;
  COFD[        2189] =   0.4606235123260915E+01;
  COFD[        2190] =  -0.3803641770348614E+00;
  COFD[        2191] =   0.1637150008428244E-01;
  COFD[        2192] =  -0.1861866210245701E+02;
  COFD[        2193] =   0.4532060302179115E+01;
  COFD[        2194] =  -0.3726516669477872E+00;
  COFD[        2195] =   0.1611699967368039E-01;
  COFD[        2196] =  -0.1863141958605046E+02;
  COFD[        2197] =   0.4533772315478532E+01;
  COFD[        2198] =  -0.3728642283605221E+00;
  COFD[        2199] =   0.1612572322702682E-01;
  COFD[        2200] =  -0.1546612761154138E+02;
  COFD[        2201] =   0.3460858520880728E+01;
  COFD[        2202] =  -0.2401164641793465E+00;
  COFD[        2203] =   0.1064270570979806E-01;
  COFD[        2204] =  -0.1648483159780408E+02;
  COFD[        2205] =   0.3836871316958563E+01;
  COFD[        2206] =  -0.2875628541513582E+00;
  COFD[        2207] =   0.1264044562964702E-01;
  COFD[        2208] =  -0.1342264266757377E+02;
  COFD[        2209] =   0.3221777556990338E+01;
  COFD[        2210] =  -0.2128673394230190E+00;
  COFD[        2211] =   0.9627787551744238E-02;
  COFD[        2212] =  -0.1752038279713604E+02;
  COFD[        2213] =   0.4959639086677385E+01;
  COFD[        2214] =  -0.4293974875408633E+00;
  COFD[        2215] =   0.1860806087249651E-01;
  COFD[        2216] =  -0.1581858256688188E+02;
  COFD[        2217] =   0.3775899039344615E+01;
  COFD[        2218] =  -0.2814687920927101E+00;
  COFD[        2219] =   0.1245218573218981E-01;
  COFD[        2220] =  -0.1677455230511400E+02;
  COFD[        2221] =   0.3946214344521416E+01;
  COFD[        2222] =  -0.3012613112287137E+00;
  COFD[        2223] =   0.1321325422150532E-01;
  COFD[        2224] =  -0.1580533640016843E+02;
  COFD[        2225] =   0.3761415063014861E+01;
  COFD[        2226] =  -0.2795015971789234E+00;
  COFD[        2227] =   0.1236332420639666E-01;
  COFD[        2228] =  -0.1991427613539219E+02;
  COFD[        2229] =   0.4985074718206674E+01;
  COFD[        2230] =  -0.3988046631693363E+00;
  COFD[        2231] =   0.1592993698509921E-01;
  COFD[        2232] =  -0.1677798202265061E+02;
  COFD[        2233] =   0.3944105793182093E+01;
  COFD[        2234] =  -0.3009768138615855E+00;
  COFD[        2235] =   0.1320048849753941E-01;
  COFD[        2236] =  -0.1678197080138354E+02;
  COFD[        2237] =   0.3942379447222446E+01;
  COFD[        2238] =  -0.3007437798097645E+00;
  COFD[        2239] =   0.1319002761262226E-01;
  COFD[        2240] =  -0.1779539137381903E+02;
  COFD[        2241] =   0.4385582017594726E+01;
  COFD[        2242] =  -0.3562391667708982E+00;
  COFD[        2243] =   0.1551072662807571E-01;
  COFD[        2244] =  -0.1779539137381903E+02;
  COFD[        2245] =   0.4385582017594726E+01;
  COFD[        2246] =  -0.3562391667708982E+00;
  COFD[        2247] =   0.1551072662807571E-01;
  COFD[        2248] =  -0.1779710613097560E+02;
  COFD[        2249] =   0.4376269742117769E+01;
  COFD[        2250] =  -0.3550500828076509E+00;
  COFD[        2251] =   0.1546030915508793E-01;
  COFD[        2252] =  -0.1772757250170383E+02;
  COFD[        2253] =   0.4347725932428951E+01;
  COFD[        2254] =  -0.3515015798062593E+00;
  COFD[        2255] =   0.1531308129463464E-01;
  COFD[        2256] =  -0.1650781936547825E+02;
  COFD[        2257] =   0.3842122413055806E+01;
  COFD[        2258] =  -0.2882098997185555E+00;
  COFD[        2259] =   0.1266705263344362E-01;
  COFD[        2260] =  -0.1931016351705685E+02;
  COFD[        2261] =   0.4757139666603682E+01;
  COFD[        2262] =  -0.3961700036764624E+00;
  COFD[        2263] =   0.1690081368592412E-01;
  COFD[        2264] =  -0.2070496369832456E+02;
  COFD[        2265] =   0.5191072368701859E+01;
  COFD[        2266] =  -0.4346910627906122E+00;
  COFD[        2267] =   0.1785842811134005E-01;
  COFD[        2268] =  -0.2072332061257778E+02;
  COFD[        2269] =   0.5194960177957292E+01;
  COFD[        2270] =  -0.4352726149968372E+00;
  COFD[        2271] =   0.1788688372961993E-01;
  COFD[        2272] =  -0.2043431836412989E+02;
  COFD[        2273] =   0.5116486808224121E+01;
  COFD[        2274] =  -0.4301998553589693E+00;
  COFD[        2275] =   0.1786838693497468E-01;
  COFD[        2276] =  -0.1894826439081699E+02;
  COFD[        2277] =   0.4665904241842780E+01;
  COFD[        2278] =  -0.3874653193184607E+00;
  COFD[        2279] =   0.1665191589170439E-01;
  COFD[        2280] =  -0.1895758159603200E+02;
  COFD[        2281] =   0.4665279128753897E+01;
  COFD[        2282] =  -0.3874162488122672E+00;
  COFD[        2283] =   0.1665125384367453E-01;
  COFD[        2284] =  -0.1962450162993201E+02;
  COFD[        2285] =   0.4875751864229406E+01;
  COFD[        2286] =  -0.4087050411778894E+00;
  COFD[        2287] =   0.1732743672805724E-01;
  COFD[        2288] =  -0.1943393371703941E+02;
  COFD[        2289] =   0.4784664368892784E+01;
  COFD[        2290] =  -0.3991077435475853E+00;
  COFD[        2291] =   0.1700210403998681E-01;
  COFD[        2292] =  -0.1944451603679417E+02;
  COFD[        2293] =   0.4785108067913721E+01;
  COFD[        2294] =  -0.3991858464067674E+00;
  COFD[        2295] =   0.1700639619802289E-01;
  COFD[        2296] =  -0.1648483159780408E+02;
  COFD[        2297] =   0.3836871316958563E+01;
  COFD[        2298] =  -0.2875628541513582E+00;
  COFD[        2299] =   0.1264044562964702E-01;
  COFD[        2300] =  -0.1771336337770503E+02;
  COFD[        2301] =   0.4289097727284545E+01;
  COFD[        2302] =  -0.3450079724750131E+00;
  COFD[        2303] =   0.1508113614518283E-01;
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
  COFTD[           4] =  -0.1441522353365304E+00;
  COFTD[           5] =  -0.7999936505724083E-04;
  COFTD[           6] =   0.4897074975189896E-07;
  COFTD[           7] =  -0.9142773649236618E-11;
  COFTD[           8] =   0.4066826061623353E+00;
  COFTD[           9] =   0.3847052789609525E-04;
  COFTD[          10] =  -0.2548469497314269E-07;
  COFTD[          11] =   0.5863025071898190E-11;
  COFTD[          12] =   0.4265800666136179E+00;
  COFTD[          13] =   0.1204072844552778E-03;
  COFTD[          14] =  -0.7672988784371693E-07;
  COFTD[          15] =   0.1520903572120744E-10;
  COFTD[          16] =   0.4128957304916273E+00;
  COFTD[          17] =   0.3905826430087406E-04;
  COFTD[          18] =  -0.2587403933152670E-07;
  COFTD[          19] =   0.5952597881665513E-11;
  COFTD[          20] =   0.2274730213366172E-01;
  COFTD[          21] =   0.6730786443559050E-03;
  COFTD[          22] =  -0.3409357669584640E-06;
  COFTD[          23] =   0.5484991439904711E-10;
  COFTD[          24] =   0.4282310128428898E+00;
  COFTD[          25] =   0.1208732836141846E-03;
  COFTD[          26] =  -0.7702684714613741E-07;
  COFTD[          27] =   0.1526789759061949E-10;
  COFTD[          28] =   0.4297895875586662E+00;
  COFTD[          29] =   0.1213132098175750E-03;
  COFTD[          30] =  -0.7730719138277033E-07;
  COFTD[          31] =   0.1532346610021863E-10;
  COFTD[          32] =   0.3247471361061310E+00;
  COFTD[          33] =   0.1777985650733727E-03;
  COFTD[          34] =  -0.1089347458861849E-06;
  COFTD[          35] =   0.2035959056112114E-10;
  COFTD[          36] =   0.3247471361061310E+00;
  COFTD[          37] =   0.1777985650733727E-03;
  COFTD[          38] =  -0.1089347458861849E-06;
  COFTD[          39] =   0.2035959056112114E-10;
  COFTD[          40] =   0.3311912921591189E+00;
  COFTD[          41] =   0.1813267307504234E-03;
  COFTD[          42] =  -0.1110964046786237E-06;
  COFTD[          43] =   0.2076359836954698E-10;
  COFTD[          44] =   0.3395573510152542E+00;
  COFTD[          45] =   0.1793350524644130E-03;
  COFTD[          46] =  -0.1101357188562818E-06;
  COFTD[          47] =   0.2064272631295983E-10;
  COFTD[          48] =   0.4306056706094065E+00;
  COFTD[          49] =   0.9359619868800346E-04;
  COFTD[          50] =  -0.6039837349070145E-07;
  COFTD[          51] =   0.1231151900773763E-10;
  COFTD[          52] =   0.2931916398206432E+00;
  COFTD[          53] =   0.4014300382755085E-03;
  COFTD[          54] =  -0.2307057834480938E-06;
  COFTD[          55] =   0.4051766204108762E-10;
  COFTD[          56] =   0.1221198807690243E+00;
  COFTD[          57] =   0.6183736619549319E-03;
  COFTD[          58] =  -0.3284226181678354E-06;
  COFTD[          59] =   0.5446035619319622E-10;
  COFTD[          60] =   0.1226934830045553E+00;
  COFTD[          61] =   0.6212781891511368E-03;
  COFTD[          62] =  -0.3299652330704521E-06;
  COFTD[          63] =   0.5471615878539925E-10;
  COFTD[          64] =   0.1403144190321918E+00;
  COFTD[          65] =   0.6012660048574510E-03;
  COFTD[          66] =  -0.3219150906946780E-06;
  COFTD[          67] =   0.5366790096442973E-10;
  COFTD[          68] =   0.3056133384265672E+00;
  COFTD[          69] =   0.3245059139327069E-03;
  COFTD[          70] =  -0.1898895909246885E-06;
  COFTD[          71] =   0.3386634966135091E-10;
  COFTD[          72] =   0.3073923760794651E+00;
  COFTD[          73] =   0.3263949291257170E-03;
  COFTD[          74] =  -0.1909949770111978E-06;
  COFTD[          75] =   0.3406349259864528E-10;
  COFTD[          76] =   0.2490175873186386E+00;
  COFTD[          77] =   0.4290366083999712E-03;
  COFTD[          78] =  -0.2426686389852732E-06;
  COFTD[          79] =   0.4208014064702138E-10;
  COFTD[          80] =   0.2727597106400999E+00;
  COFTD[          81] =   0.3944027517533329E-03;
  COFTD[          82] =  -0.2258005402743759E-06;
  COFTD[          83] =   0.3953256687340045E-10;
  COFTD[          84] =   0.2740370684853061E+00;
  COFTD[          85] =   0.3962497747170277E-03;
  COFTD[          86] =  -0.2268579841721363E-06;
  COFTD[          87] =   0.3971770137995339E-10;
  COFTD[          88] =   0.4313313898218810E+00;
  COFTD[          89] =   0.9205368774727828E-04;
  COFTD[          90] =  -0.5945097223307069E-07;
  COFTD[          91] =   0.1214380118036443E-10;
  COFTD[          92] =   0.4010129334861676E+00;
  COFTD[          93] =   0.1972528443628177E-03;
  COFTD[          94] =  -0.1216981610244689E-06;
  COFTD[          95] =   0.2294088734064920E-10;
  COFTD[          96] =   0.1441522353365304E+00;
  COFTD[          97] =   0.7999936505724083E-04;
  COFTD[          98] =  -0.4897074975189896E-07;
  COFTD[          99] =   0.9142773649236618E-11;
  COFTD[         100] =   0.0000000000000000E+00;
  COFTD[         101] =   0.0000000000000000E+00;
  COFTD[         102] =   0.0000000000000000E+00;
  COFTD[         103] =   0.0000000000000000E+00;
  COFTD[         104] =   0.2352832283423578E+00;
  COFTD[         105] =   0.4656706350120388E-03;
  COFTD[         106] =  -0.2609398458577575E-06;
  COFTD[         107] =   0.4492718580048779E-10;
  COFTD[         108] =   0.1798404103253150E+00;
  COFTD[         109] =   0.6017229463172420E-03;
  COFTD[         110] =  -0.3264339193657546E-06;
  COFTD[         111] =   0.5491123423829806E-10;
  COFTD[         112] =   0.2370534621354604E+00;
  COFTD[         113] =   0.4691742672103970E-03;
  COFTD[         114] =  -0.2629031159826994E-06;
  COFTD[         115] =   0.4526521083989819E-10;
  COFTD[         116] =  -0.1743526120296758E+00;
  COFTD[         117] =   0.8622469282596517E-03;
  COFTD[         118] =  -0.3795455129739271E-06;
  COFTD[         119] =   0.5602621271895701E-10;
  COFTD[         120] =   0.1801870770056730E+00;
  COFTD[         121] =   0.6028828485656729E-03;
  COFTD[         122] =  -0.3270631648338786E-06;
  COFTD[         123] =   0.5501708305866783E-10;
  COFTD[         124] =   0.1805137892277234E+00;
  COFTD[         125] =   0.6039759857559981E-03;
  COFTD[         126] =  -0.3276561903444169E-06;
  COFTD[         127] =   0.5511683912195231E-10;
  COFTD[         128] =   0.9907533094923988E-01;
  COFTD[         129] =   0.6442014322917300E-03;
  COFTD[         130] =  -0.3384859791884526E-06;
  COFTD[         131] =   0.5573567872490878E-10;
  COFTD[         132] =   0.9907533094923988E-01;
  COFTD[         133] =   0.6442014322917300E-03;
  COFTD[         134] =  -0.3384859791884526E-06;
  COFTD[         135] =   0.5573567872490878E-10;
  COFTD[         136] =   0.1000392102974500E+00;
  COFTD[         137] =   0.6504687084211580E-03;
  COFTD[         138] =  -0.3417790254177462E-06;
  COFTD[         139] =   0.5627791733432544E-10;
  COFTD[         140] =   0.1051242233178346E+00;
  COFTD[         141] =   0.6506660047178345E-03;
  COFTD[         142] =  -0.3425645641009918E-06;
  COFTD[         143] =   0.5648041610405423E-10;
  COFTD[         144] =   0.2001200095273127E+00;
  COFTD[         145] =   0.5647937464189248E-03;
  COFTD[         146] =  -0.3094455087464756E-06;
  COFTD[         147] =   0.5241393748626840E-10;
  COFTD[         148] =  -0.2003084371216687E-01;
  COFTD[         149] =   0.8504401712698103E-03;
  COFTD[         150] =  -0.4210644955951721E-06;
  COFTD[         151] =   0.6679597524516173E-10;
  COFTD[         152] =  -0.1609811732564022E+00;
  COFTD[         153] =   0.9038076297512394E-03;
  COFTD[         154] =  -0.4069279669787991E-06;
  COFTD[         155] =   0.6092022913992956E-10;
  COFTD[         156] =  -0.1613574734420234E+00;
  COFTD[         157] =   0.9059203176635061E-03;
  COFTD[         158] =  -0.4078791780205072E-06;
  COFTD[         159] =   0.6106263270843215E-10;
  COFTD[         160] =  -0.1312444259932660E+00;
  COFTD[         161] =   0.9039014427379127E-03;
  COFTD[         162] =  -0.4178315332949202E-06;
  COFTD[         163] =   0.6357257058426966E-10;
  COFTD[         164] =   0.1632460981654722E-01;
  COFTD[         165] =   0.7901334428207009E-03;
  COFTD[         166] =  -0.3982924858072905E-06;
  COFTD[         167] =   0.6388514753031584E-10;
  COFTD[         168] =   0.1637184926975502E-01;
  COFTD[         169] =   0.7924198969668966E-03;
  COFTD[         170] =  -0.3994450474585499E-06;
  COFTD[         171] =   0.6407001562035587E-10;
  COFTD[         172] =  -0.5087437787577148E-01;
  COFTD[         173] =   0.8543426436721789E-03;
  COFTD[         174] =  -0.4159264807346257E-06;
  COFTD[         175] =   0.6530633032370990E-10;
  COFTD[         176] =  -0.2716895748891503E-01;
  COFTD[         177] =   0.8372331893523684E-03;
  COFTD[         178] =  -0.4128876639137287E-06;
  COFTD[         179] =   0.6534052395309522E-10;
  COFTD[         180] =  -0.2723227826735267E-01;
  COFTD[         181] =   0.8391844698645221E-03;
  COFTD[         182] =  -0.4138499521537908E-06;
  COFTD[         183] =   0.6549280851690021E-10;
  COFTD[         184] =   0.2015217542366260E+00;
  COFTD[         185] =   0.5627441294507066E-03;
  COFTD[         186] =  -0.3085192624064295E-06;
  COFTD[         187] =   0.5228060239663243E-10;
  COFTD[         188] =   0.1221940303626176E+00;
  COFTD[         189] =   0.6903211781262630E-03;
  COFTD[         190] =  -0.3648449022916740E-06;
  COFTD[         191] =   0.6030549194456049E-10;
};




#if 0




\\
\\
\\  This is the mechanism file
\\
\\
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!     
! Reduced version of GRI-MECH 1.2. 22 species ( + N2, AR); 104 reactions. !     
!                                 PennState,  Dec, 1994.                  !     
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!     
ELEMENTS                                                                        
O  H  C  N  AR                                                                  
END                                                                             
SPECIES                                                                         
H2      H       O       O2      OH      H2O     HO2     H2O2                    
CH2     CH2(S)  CH3     CH4     CO      CO2     HCO     CH2O                    
CH3O    C2H2    C2H3    C2H4    C2H5    C2H6                                    
N2      AR                                                                      
END                                                                             
REACTIONS                                                                       
O+H+M<=>OH+M                             5.000E+17   -1.000      0.00           
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
O+H2<=>H+OH                              5.000E+04    2.670   6290.00           
O+HO2<=>OH+O2                            2.000E+13    0.000      0.00           
O+CH2<=>H+HCO                            8.000E+13    0.000      0.00           
O+CH2(S)<=>H+HCO                         1.500E+13    0.000      0.00           
O+CH3<=>H+CH2O                           8.430E+13    0.000      0.00           
O+CH4<=>OH+CH3                           1.020E+09    1.500   8600.00           
O+CO+M<=>CO2+M                           6.020E+14    0.000   3000.00           
H2/2.00/ O2/6.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/3.50/ C2H6/3.00/ AR/0.50/    
O+HCO<=>OH+CO                            3.000E+13    0.000      0.00           
O+HCO<=>H+CO2                            3.000E+13    0.000      0.00           
O+CH2O<=>OH+HCO                          3.900E+13    0.000   3540.00           
O+C2H2<=>CH2(S)+CO                       1.020E+07    2.000   1900.00           
O+C2H2<=>CO+CH2                          1.020E+07    2.000   1900.00           
O+C2H4<=>CH3+HCO                         1.920E+07    1.830    220.00           
O+C2H5<=>CH3+CH2O                        1.320E+14    0.000      0.00           
O+C2H6<=>OH+C2H5                         8.980E+07    1.920   5690.00           
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
H+HO2<=>O2+H2                            2.800E+13    0.000   1068.00           
H+HO2<=>2OH                              1.340E+14    0.000    635.00           
H+H2O2<=>HO2+H2                          1.210E+07    2.000   5200.00           
H+CH2(+M)<=>CH3(+M)                      2.500E+16   -0.800      0.00           
     LOW  /  3.200E+27   -3.140   1230.00/                                      
     TROE/  0.6800   78.00  1995.00  5590.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
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
H+CH2O(+M)<=>CH3O(+M)                    5.400E+11    0.454   2600.00           
     LOW  /  2.200E+30   -4.800   5560.00/                                      
     TROE/  0.7580   94.00  1555.00  4200.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/                      
H+CH2O<=>HCO+H2                          2.300E+10    1.050   3275.00           
H+CH3O<=>OH+CH3                          3.200E+13    0.000      0.00           
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
H+C2H6<=>C2H5+H2                         1.150E+08    1.900   7530.00           
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
OH+H2O2<=>HO2+H2O                        5.800E+14    0.000   9560.00           
OH+CH2<=>H+CH2O                          2.000E+13    0.000      0.00           
OH+CH2(S)<=>H+CH2O                       3.000E+13    0.000      0.00           
OH+CH3<=>CH2+H2O                         5.600E+07    1.600   5420.00           
OH+CH3<=>CH2(S)+H2O                      2.501E+13    0.000      0.00           
OH+CH4<=>CH3+H2O                         1.000E+08    1.600   3120.00           
OH+CO<=>H+CO2                            4.760E+07    1.228     70.00           
OH+HCO<=>H2O+CO                          5.000E+13    0.000      0.00           
OH+CH2O<=>HCO+H2O                        3.430E+09    1.180   -447.00           
OH+C2H2<=>CH3+CO                         4.830E-04    4.000  -2000.00           
OH+C2H3<=>H2O+C2H2                       5.000E+12    0.000      0.00           
OH+C2H4<=>C2H3+H2O                       3.600E+06    2.000   2500.00           
OH+C2H6<=>C2H5+H2O                       3.540E+06    2.120    870.00           
2HO2<=>O2+H2O2                           1.300E+11    0.000  -1630.00           
 DUPLICATE                                                                      
2HO2<=>O2+H2O2                           4.200E+14    0.000  12000.00           
 DUPLICATE                                                                      
HO2+CH2<=>OH+CH2O                        2.000E+13    0.000      0.00           
HO2+CH3<=>O2+CH4                         1.000E+12    0.000      0.00           
HO2+CH3<=>OH+CH3O                        2.000E+13    0.000      0.00           
HO2+CO<=>OH+CO2                          1.500E+14    0.000  23600.00           
HO2+CH2O<=>HCO+H2O2                      1.000E+12    0.000   8000.00           
CH2+O2<=>OH+HCO                          1.320E+13    0.000   1500.00           
CH2+H2<=>H+CH3                           5.000E+05    2.000   7230.00           
2CH2<=>H2+C2H2                           3.200E+13    0.000      0.00           
CH2+CH3<=>H+C2H4                         4.000E+13    0.000      0.00           
CH2+CH4<=>2CH3                           2.460E+06    2.000   8270.00           
CH2(S)+N2<=>CH2+N2                       1.500E+13    0.000    600.00           
CH2(S)+AR<=>CH2+AR                       9.000E+12    0.000    600.00           
CH2(S)+O2<=>H+OH+CO                      2.800E+13    0.000      0.00           
CH2(S)+O2<=>CO+H2O                       1.200E+13    0.000      0.00           
CH2(S)+H2<=>CH3+H                        7.000E+13    0.000      0.00           
CH2(S)+H2O<=>CH2+H2O                     3.000E+13    0.000      0.00           
CH2(S)+CH3<=>H+C2H4                      1.200E+13    0.000   -570.00           
CH2(S)+CH4<=>2CH3                        1.600E+13    0.000   -570.00           
CH2(S)+CO<=>CH2+CO                       9.000E+12    0.000      0.00           
CH2(S)+CO2<=>CH2+CO2                     7.000E+12    0.000      0.00           
CH2(S)+CO2<=>CO+CH2O                     1.400E+13    0.000      0.00           
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
CH3+C2H4<=>C2H3+CH4                      2.270E+05    2.000   9200.00           
CH3+C2H6<=>C2H5+CH4                      6.140E+06    1.740  10450.00           
HCO+H2O<=>H+CO+H2O                       2.244E+18   -1.000  17000.00           
HCO+M<=>H+CO+M                           1.870E+17   -1.000  17000.00           
H2/2.00/ H2O/0.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/                      
HCO+O2<=>HO2+CO                          7.600E+12    0.000    400.00           
CH3O+O2<=>HO2+CH2O                       4.280E-13    7.600  -3530.00           
C2H3+O2<=>HCO+CH2O                       3.980E+12    0.000   -240.00           
C2H4(+M)<=>H2+C2H2(+M)                   8.000E+12    0.440  88770.00           
     LOW  /  7.000E+50   -9.310  99860.00/                                      
     TROE/  0.7345  180.00  1035.00  5417.00 /                                  
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ C2H6/3.00/ AR/0.70/             
C2H5+O2<=>HO2+C2H4                       8.400E+11    0.000   3875.00           
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
